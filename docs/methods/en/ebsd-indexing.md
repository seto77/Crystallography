This is a draft written by Claude from the code (2026-11-05). Not yet reviewed by the author.

# EBSD indexing, orientation search, and geometry calibration

Target commit: `9ffcf47` (main). Files: `EBSD/EbsdBandDetector.cs`, `EBSD/EbsdIndexer.cs`, `EBSD/EbsdRadonIndexer.cs`, `EBSD/EbsdDictionaryIndexer.cs`, `EBSD/EbsdZoneAxisIndexer.cs`, `EBSD/EbsdOrientationSearch.cs`, `EBSD/EbsdGeometryCalibrator.cs`, `EBSD/EbsdPatternProjector.cs` (`EbsdPatternScorer`).

Orientations are rotation matrices $R$ (crystal → sample frame, $\mathbf v^{\rm s}=R\mathbf v^{\rm c}$), and the geometry is `EbsdDetectorGeometry` ([ebsd-geometry.md](ebsd-geometry.md)). The reflection input for every path is `Crystal.VectorOfG_KikuchiLine` ($\mathbf g$ in the crystal Cartesian frame, $d$ [nm], `RelativeIntensity` = kinematical intensity).

---

## 0. Overview

| Path | Input | Objective | Implementation |
|---|---|---|---|
| Band detection → pair-angle | Measured image | Angular residuals of the assignment, width likelihood | `EbsdBandDetector.Detect` → `EbsdIndexer.Index` |
| Radon template matching | Measured image | SNR-normalized sum of Radon evidence | `EbsdBandDetector.ComputeRadonMap` → `EbsdRadonIndexer.Index` |
| Dictionary matching | Measured image + master pattern | ZNCC after robust preprocessing | `EbsdDictionaryIndexer.Index` |
| Zone axis | Band intersections picked by hand | Angular residual + preference for low indices | `EbsdZoneAxisIndexer.Index` |
| Integrated (GUI Find) | The above + master pattern | Composite rank (Radon z + ZNCC) | `EbsdOrientationSearch.Run` |
| Geometry calibration | Measured + master pattern + orientation | ZNCC | `EbsdGeometryCalibrator.Run` |

---

## 1. Common: reflection catalog

`EbsdIndexer.BuildCatalog` / `EbsdRadonIndexer.BuildCatalog`:

- The sign of $\hat{\mathbf g} = \mathbf g/|\mathbf g|$ is made positive in lexicographic order ($z>0$, then $y>0$, then $x>0$), and those with $\hat{\mathbf g}\cdot\hat{\mathbf g}' > 0.99995$ (about 0.6°) are merged into one node ($\pm\mathbf g$ and harmonic reflections).
- `EbsdIndexer`: a node holds the list of reflections, the list of $d$, and the maximum intensity; the top `maxCatalogNodes` (default 120) by intensity.
- `EbsdRadonIndexer`: a node holds a representative (strongest) reflection and a weight $w$. The top `maxNodes` (`DefaultMaxNodes`) by intensity. The weight is $w=\min(\sqrt I,\ 3\,\mathrm{median}\sqrt I)$ if `weightExponent` = 0.5 (default), otherwise $w=\mathrm{clip}((I/\mathrm{median}I)^{p},\,0.5,\,2)$.

---

## 2. Band detection (`EbsdBandDetector.Detect`)

Procedure in the doc: downscale to a long side of 256 px → Gaussian background division → normalization → Radon (sum/√N normalization) → butterfly convolution with a width bank → MAD threshold + in-image line-distance NMS → local θ/ρ re-search → edges and width from the transverse profile.

| Stage | Content | Value |
|---|---|---|
| Downscale | Long side `WorkLongSide` | 256 px (position constants are converted with `WS = 256/256`) |
| Background | Division by a wide Gaussian (normalized convolution excluding invalid pixels from the weights) | $\sigma = 0.20\times$short side (changed from 0.1 in 260724Cl) |
| Smoothing, normalization | Gaussian → mean 0, variance 1 → ±5σ clip | $\sigma = 0.8$ |
| Radon | $\theta\in[0,180°)$, $\rho$ is the signed distance from the image center [px]. Line integral = bilinear, sum in 1 px steps / √N | $\Delta\theta = 0.5°$, $\Delta\rho = 1$ px |
| butterfly | $K_w(t) = G_\sigma(t) - \tfrac12G_\sigma(t-w/2) - \tfrac12G_\sigma(t+w/2)$, $\sigma=\max(0.75, 0.12w)$, mean 0, L2 norm 1. At each $(\theta,\rho)$ the response of the width with the largest $\lvert s\rvert$ is kept **with its sign** (bright bands positive, dark bands negative) | $w\in\{3,5,7,10,14,20,28,38\}$ px |
| θ smoothing | 3-tap filter adding only neighbors of the same sign as the center ($\rho$ flipped between $\theta$ and $\theta+180°$) | 1 pass |
| Threshold | $\lvert\text{response}\rvert > \mathrm{median} + 4.5\,\mathrm{MAD}$ and a local maximum along $\rho$ | |
| NMS | Same sign: in-image RMS distance 6 px. Opposite sign and nearly parallel ($\Delta\theta<10°$): $\max(6,\ 0.75\max w)$ (sidelobe suppression). Opposite sign, near-parallel, score ratio ≥ 0.7 is tagged as a "bright/dark lobe pair" | |
| Local re-search | $\rho$ over $\pm2$ px in 0.2 px steps, $\theta$ also in the neighborhood | |
| Duplicate removal | Same sign 5 px, opposite sign near-parallel $\max(5, 0.75\max w)$ | |
| Line verification gate | Discards bands with `CenterQuality` < (0.3 if Score ≥ 17, otherwise 0.5), or (not a lobe pair, Score < 17 and `EdgeLinearity` ≥ 0.55), or `WidthQuality` ≤ 0.05 | |

Output `EbsdBand` (native pixel coordinates): normalized center line $A\,\mathrm{col}+B\,\mathrm{row}+C=0$ ($A^2+B^2=1$), anchor points on the center line, points on both edges (each anchor ± half width), width, Score (in MAD units), CenterQuality (median of narrow-window anchor projections), WidthQuality, and others. `maxBands` default 12.

`ComputeRadonMap` shares the same steps from preprocessing through θ smoothing, does no discrete detection, and returns the map of $|\text{response}|$, `EbsdRadonMap.Abs`, together with robust null statistics ($\mu_0$ = median, $\sigma_0 = 1.4826\,\mathrm{MAD}$).

---

## 3. Pair-angle indexing (`EbsdIndexer.Index`)

### 3.1 Observations

For each band, the plane normal (lab) $\mathbf g^{\rm lab}$ = `LineToLabNormal`, and in the sample frame $\hat{\mathbf m} = R_x(\tau)\mathbf g^{\rm lab}$ (sign undetermined). $\sin\theta_B = |\hat{\mathbf g}\cdot\mathbf P_e|/|\mathbf P_e|$ is computed for each edge point; if there are 2 or more points and WidthQuality > 0.05, the median $\sin\theta_B$ and $\sigma = \max(1.4826\,\mathrm{MAD},\ 0.2\sin\theta_B)$. Weight $w_o = 0.3 + 0.7\,\text{CenterQuality}$.

### 3.2 Procedure

1. Theoretical pair table: the angles $\arccos|\hat{\mathbf g}_p\cdot\hat{\mathbf g}_q|\in(8°, 85°)$ between pairs of catalog directions, sorted by angle.
2. Observations: the top 8 by weight. For each theoretical pair within ± `pairTolDeg` (1.5°) of the angle of an observed pair $(i,j)$, the Wahba problem is solved for 2 correspondence swaps × 4 signs of the observed normals.
3. **Wahba / Kabsch** (`SolveWahba`): $H=\sum w\,\mathbf m\mathbf g^{\mathsf T}$, $H = U\Sigma V^{\mathsf T}$,

$$ R = U\,\mathrm{diag}\bigl(1,1,\det(UV^{\mathsf T})\bigr)\,V^{\mathsf T}\qquad(\text{undetermined if }\sigma_2<10^{-10}). $$

4. Scoring (`ScoreAndAssign`): pairs whose angle $\psi = \arccos|\hat{\mathbf m}\cdot R\hat{\mathbf g}|$ between the prediction $R\hat{\mathbf g}_n$ and the observation $\hat{\mathbf m}_o$ is below `assignTolDeg` (2.5°) are greedily assigned **one-to-one** in ascending order of residual.

$$ \text{Score} = \sum_{\rm assigned} w_o\,\Bigl(1-\frac{\psi}{\psi_{\rm tol}}\Bigr)\,f_w,\qquad
f_w = 0.4 + 0.6\exp\!\Bigl(-\tfrac12\min(r,4)^2\Bigr),\quad r = \min_m\frac{|\sin\theta_B - \lambda/(2d_m)|}{\sigma}. $$

   ($f_w=1$ if there is no $\sin\theta_B$.) $\lambda$ in nm. The assignment is the $(hkl)$ of the harmonic reflection with the closest $d_m$. Solutions with 3 or more assignments become candidates.
5. The top 60 are refined up to twice by "weighted Kabsch with the observed normal signs aligned to the assignment" (updated only when the score increases).
6. Duplicates with misorientation below 3° (without symmetry) are removed, keeping the top `maxCandidates` (10). `AngularRmsDeg` is also returned.

---

## 4. Radon template matching (`EbsdRadonIndexer.Index`)

### 4.1 Objective

For each predicted band center line $k$ intersecting the detector under orientation $R$ (and geometry), from the sample-frame $R\hat{\mathbf g}_k$ → lab → the same coefficients $(a,b,c)$ [mm] as `LabNormalToLine`: $\theta_k=\operatorname{atan2}(b,a)$, $\rho_k = -c/(s\sqrt{a^2+b^2})\times\text{Scale}$ (work px).

$$ \text{Score}(R) = \frac{\sum_k w_k\,\psi\bigl(E(\theta_k,\rho_k)-\mu_0\bigr)}{\sqrt{\sum_k w_k^2}\;(\sigma_0\ \text{if cap}=0)},\qquad
\psi(e) = \begin{cases}\max(e,-\sigma_0) & \text{cap}=0\\ \text{cap}\cdot\tanh\bigl(z/\text{cap}\bigr)\ (z>0),\ \max(z,-1)\ (z\le0),\ z=e/\sigma_0 & \text{cap}>0\end{cases} $$

(`ScoreExactCore`). $E$ is the bilinear value of the evidence map ($\rho$ outside the map gives 0). Among closely spaced predicted lines ($\Delta\theta<2°$ and $\Delta\rho<5$ px, taking into account $\theta$ periodicity and $\rho$ flip), only the first one in intensity order scores (prevention of double scoring). Invalid if $(\sum w)^2/\sum w^2 < 4$.

### 4.2 Search

- Seeds ①: local maxima of the evidence map (above $\mu_0+3\sigma_0$, NMS) taken as pseudo-bands for pair-angle + Kabsch (§3). Comment: because there were cases where grid + NM alone could not reach the ridge (θ width ~1°).
- Seeds ②: the SO(3) coarse grid $R = R_{\rm sphere}(i)R_z(\phi)$ (step `coarseStepDeg` = 3°) is scored with the top `coarseNodes` = 20 nodes by intensity and a **dilated map** (window max over $\theta\pm$`dilTheta` rows × $\rho\pm$`dilRho` columns, nearest-neighbor lookup). The dilation width is the amount a predicted line moves under a half-step rotation: `dilRho` $=\max(3,\lceil 0.7\sqrt{W_w^2+H_w^2}\,\delta/2\rceil)$, `dilTheta` $=\max(2,\lceil(\delta/2)/\Delta\theta\rceil)$. Top `maxGridSeeds` (52).
- Refinement: the exact score is maximized with Nelder–Mead. Variables are $[\omega_x,\omega_y,\omega_z]$ (sample-frame left perturbation, deg), plus $(\Delta Z,\Delta X,\Delta Y)$ [mm] if geometry search is enabled.
- Reduction of equivalents: if $Q=R_1^{\mathsf T}R_2$ preserves the set of catalog directions (including ±), the two are indistinguishable in the pattern and are reduced to one representative. Candidates with different geometries are kept separately (260922Cl).
- Geometry search (260922Cl, disabled by default): `zSearchRangeMm` (searches DetZ on a grid; the default step is the dilation width converted back to mm), `xSearchRangeMm`, `xyRefineRangeMm` (moves X and Y only during refinement), `coarseGeometryToleranceMm` (added to the dilation width).

---

## 5. Dictionary matching (`EbsdDictionaryIndexer.Index`)

### 5.1 Preprocessing and similarity (`EbsdPatternScorer`)

- `Downsample`: box downscaling to long side `targetLongSide` (lower bound 8 px).
- `RobustPreprocess` (applied to **both** the measured and the simulated pattern): wide background Gaussian with $\sigma = 0.1\times$short side, $v = \ln\bigl(\max(x, 0.01f)/\max(\mathrm{bg}, f)\bigr)$ ($f = \max(10^{-10}, 0.05\bar x)$) → DoG ($\sigma_1=1.5$, $\sigma_2=6$) → standardization → $\tanh(z/3)$ → re-standardization. Doc: a dynamical single slice has a heavy-tailed distribution whose variance is dominated by bright zone-axis spots, and with raw ZNCC the correct solution loses to false orientations.
- `RobustPreprocessFast`: approximates the Gaussian by 3 box passes ($r = \mathrm{round}((\sqrt{4\sigma^2+1}-1)/2)$, combined variance $\sigma^2$), with SIMD (`Vector.Log`), reuse of work buffers, and sequential execution.
- ZNCC: the reference $\mathbf a$ has mean 0 and variance 1 (population variance), $\mathrm{ZNCC} = \sum a_i(b_i-\bar b)/(n\,\sigma_b)$, 0 if $\sigma_b<10^{-12}$. `ZnccSimd` differs only in the summation order (calibration only).

### 5.2 Three-stage structure

1. **Coarse stage**: Fibonacci sphere with $n_s = \max(64,\lfloor4\pi/\delta^2\rfloor)$ points × $n_\phi=\max(16,\lfloor360/\delta\rfloor)$ in-plane rotations ($\delta$ = `coarseStepDeg`, 3° in the GUI path). The square grid uses the in-plane rotation decomposition projection ([ebsd-master-pattern.md](ebsd-master-pattern.md) §4.2); the hexagonal grid uses the normal `Project`. The GUI always uses `thoroughCoarse = true` (96 px, `RobustPreprocessFast`). The doc says `false` (48 px, lightweight preprocessing: mean-normalized log-ratio → box DoG 5×5 − 21×21 → tanh → normalization) is for comparison in the validation harness. NMS discarding misorientations below $1.5\delta$ keeps the top `coarseKeep` (64) basins.
2. **Middle stage**: the top ones are re-scored with full robust preprocessing (96 px), keeping the top `refineKeep` (12). Note in the doc: with `thoroughCoarse = true`, the reference, resolution and projector are the same as in the coarse stage, so this is "re-scoring with the same pipeline".
3. **Fine stage**: Nelder–Mead (sample-frame left perturbation; initial step $\delta/2$ with 120 evaluations → 0.5° with 80 evaluations) → misorientation NMS → top `maxCandidates`. The returned Score is NaN (adding Radon z and the composite rank is done by the integration layer).

Reduction by symmetry (`properSymmetries`): within an equivalence class $\{RS\}$, only the representative with maximal $\operatorname{tr}$ (= smallest rotation angle) is evaluated. For $R=r_0R_z$, $\operatorname{tr}(R_z M)$ is compared in closed form for $M=r_0$ and $M=Sr_0$. Correction in the doc (260725Cl): the tie band $10^{-12}$ rescues only strictly symmetry-invariant points, and coverage is maintained by the mechanism that "another grid point lies near the orbit partner". For this reason `EbsdOrientationSearch` limits FZ exclusion to a single proper rotation (C2).

---

## 6. Indexing from zone axes (`EbsdZoneAxisIndexer.Index`)

- Catalog: irreducible $[uvw]$ are formed from integer cross products of reflections and limited to `maxIndex` (4) or less. Weight = sum of the relative intensities of the reflections belonging to that zone. The top `maxNodes` (300), in both ± directions. Doc: allowing high indices makes it possible to explain any point (for Botallackite, $[10\,\bar2\,5]$ was assigned).
- Observations: picked points $(\mathrm{col},\mathrm{row})$ → **exit direction** $\hat{\mathbf o} = -$`PixelToSampleDirection`.
- Seeds: exhaustive search over point pairs whose angle matches a catalog pair within `toleranceDeg` (2°), followed by Kabsch. If indices are entered by hand, the assignment of that point is restricted to $\pm[uvw]$; if 2 or more points are entered by hand, the pairs are restricted to those.
- Scoring (`Evaluate`): each point is assigned to the closest (largest dot product) rotated catalog direction, and

$$ \text{Score} = \sum_{\rm assigned}\Bigl[\Bigl(1-\frac{\psi}{\psi_{\rm tol}}\Bigr) + 0.5\,\frac{W_c}{W_{\max}}\Bigr] + 2\,n_{\rm assigned}. $$

  Solutions in which different points are assigned to the same zone axis, or in which a hand-entered point is not assigned, are rejected. Returns `RmsDeg` and `RmsPx` (pixel RMS between the predicted position `SampleDirectionToPixel` and the picked position).
- Optional (`refineGeometry`): with the assignment fixed, 3 orientation + 3 detector center ($\Delta X_d,\Delta Y_d,\Delta Z_d$) parameters are optimized simultaneously with Nelder–Mead, minimizing the pixel residual. With fewer than 4 points, the geometry is not changed.

---

## 7. Integrated search (`EbsdOrientationSearch.Run`)

1. Engine: dictionary matching if `useDictionary`, otherwise Radon matching. Dictionary candidates have the Radon z attached afterwards.
2. Composite rank (when `context` is present): each candidate is projected, the robust ZNCC is computed, **standardized within the candidate set**, and

$$ \text{combo} = z_{\rm Radon} + 0.5\cdot\mathrm{clip}(z_{\rm ZNCC},\,\pm2). $$

   Radon scoring uses the evidence saturation cap = `SaturateCap` (8) (doc: it degrades on its own, so it is used only together with the composite).
3. Refinement: only the composite top 1 is refined by ZNCC (Nelder–Mead, initial step 0.25°, 120 evaluations). If the drop in Radon z exceeds `GuardMaxZDrop` (0.2), the result is not adopted (false-convergence guard).
4. Finishing: the final orientation is re-polished with **the same objective and the same steps** as the last stage of geometry calibration (`EbsdGeometryCalibrator.CreateFinalScorer`; Nelder–Mead with initial step 0.7° and 150 evaluations → `OrientationPolishStepDeg` = 0.1° with 100 evaluations). Doc: with different objectives, the orientation went back and forth by about 1° when repeating Find and Calibrate. The same guard applies to finishing.
5. Configuration that also searches the geometry (`EbsdSearchOptions.GeometryRobust`, 260922Cl): `ZSearchRangeMm = 12`, `XSearchRangeMm = 3` (step 3), `CalibrateTopCandidates = 5` (the top 5 are lightly calibrated with 5 starting points and 1 stage, then re-sorted by post-calibration ZNCC), `CalibrateIncludePlain = true` (the top results of the plain search with fixed geometry are also sent to calibration), `DictionaryZHypotheses = 1`, `DictionaryMergeRadon = true`, `SaturateCap = 4`, `RadonMaxGridSeeds = 150`. When candidates have been calibrated, the fine adjustments of §7.3 and 7.4 are not performed (doc: they shifted an orientation that calibration had matched to 0.03° by 0.66°).

`EbsdMatchingContext`: geometry, master pattern, the plane used for matching (MC-weighted composition or a single slice), the downscaled reference, raster dimensions, the current orientation, the full-resolution display values `DisplayReference`, and the FWHM of the high-pass filter applied to the simulation side, `SimFlattenFwhmPx` (in measured-image px).

Main values of `EbsdSearchOptions.Default`: Radon coarse step 3°, coarse nodes 20, `RadonMaxGridSeeds` 52, dictionary coarse step 3°, `DictionaryCoarseKeep` 64, `DictionaryRefineKeep` 12, `DictionaryThoroughCoarse = true`, `ZnccCoef` 0.5, `ZnccClip` 2, `GuardMaxZDrop` 0.2, geometry search 0.

---

## 8. Geometry calibration (`EbsdGeometryCalibrator.Run`)

### 8.1 Variables and constraints

- Variables: physical in-plane PC coordinates $(u_F, v_F)$ [mm], $\ln DD$, 3 orientation parameters (sample-frame left perturbation). **$\delta$ (DetTilt) is fixed** (doc: for a single pattern, DetTilt and the $x$ rotation of the orientation form a gauge freedom).
- Conversion back to geometry: `FromPatternCenter(u_F, v_F, e^{\ln DD}, δ)` ([ebsd-geometry.md](ebsd-geometry.md) §3.4, question G-I1).
- Soft bounds: if any of $|\Delta u_F| > 0.25\,W_{\rm phys}$, $|\Delta v_F| > 0.25\,H_{\rm phys}$, $|\Delta\ln DD| > 0.35$ holds relative to the starting point, the penalty $10 + |\Delta u|/W + |\Delta v|/H + |\Delta\ln DD|$ is returned.

### 8.2 Objective

$-\mathrm{ZNCC}$. The reference is `DisplayReference` (the displayed values, reflecting the user's "flatten background" setting), box-downscaled to each stage's long side and normalized. The simulation side is projected with `EbsdPatternProjector`, and if `SimFlattenFwhmPx` > 0, a high-pass filter with the same FWHM (converted by the stage's downscaling factor) is applied (`SubtractBoxBackground`, which subtracts a Gaussian background approximated by 3 box passes).

### 8.3 Procedure

- 3 resolution stages: long side 160 → 480 → full resolution (`StageLongSides` = [160, 480, 0]). Stage convergence criteria `(Tol, XTol)` = (1e-6, 1e-2), (1e-7, 3e-3), (1e-9, 3e-4). Evaluation limits (orientation / geometry / 6 variables) = (150, 120, 600), (200, 180, 900), (300, 250, 1500).
- Calibration from one starting point: alternating method (① orientation with geometry fixed [initial step 0.7°], ② geometry with orientation fixed [step = 1 % of detector width and height, $\ln DD$ 0.02]) for at most `MaxRounds` = 20 rounds, stopping when the improvement in one round is below `ZnccTolerance` = 1e-4 → orientation finishing → simultaneous Nelder–Mead over 6 variables (in the last stage only, `NelderMeadRestart` restarts up to 3 times with the step reduced to 1/4).
- Multi-start: `StartCount` = 40 points at stage 0. Offsets are a Halton sequence on $[-1,1]^3$ ([0] is the current geometry) × (PC: 8 % of detector width and height, $\ln DD$: 0.08). No random numbers are used, so the same input gives the same result. The top $\max(6,$ number of CPUs$)$ are carried over to stage 1, and 1 point to stage 2.
- Diagnostics: index of the best starting point, ZNCC difference between best and worst, number of starting points within 1e-3 of the best, and the spread of PC and DD within that cluster (`FlatU/V/Dd`).

---

## 9. Implementation questions

### I-I1 The doc and the value of `SaturateCap` do not match

The doc of `EbsdOrientationSearch.SaturateCap` says "determined by a coefficient sweep (toward the middle of the plateau 0.4-1.0)", but the value is 8 (4 in `GeometryRobust`). It is unclear whether the plateau range refers to a different normalization (e.g. ψ not in units of $\sigma_0$) or the doc is outdated.

### I-I2 Duplicate removal in pair-angle indexing does not consider symmetry

The last stage of `EbsdIndexer.Index` regards candidates with `MisorientationDeg(c.Rotation, current.Rotation)` (without `properSymmetries`) below 3° as duplicates. Symmetry-equivalent orientations remain as separate candidates and can consume the top `maxCandidates` slots (the doc states explicitly that "symmetry equivalents may be included", so this can also be read as by design).

### I-I3 The starting point of geometry calibration depends on the sign convention of `FromPatternCenter`

As in §8.1, calibration builds the geometry from $(u_F, v_F, \ln DD)$ with `FromPatternCenter`. Under the condition of G-I1 in [ebsd-geometry.md](ebsd-geometry.md) (configurations with $\mathbf n\cdot\mathbf C<0$), starting point [0] (which should be "the current geometry itself") does not coincide with the current geometry.

---

## 10. Items to confirm

- **I-1** The catalog intensity `RelativeIntensity` is the kinematical intensity (definition of `VectorOfG_KikuchiLine`). Both Radon and pair-angle set weights independently of dynamical intensities. Validity for crystals with heavy elements and low symmetry.
- **I-2** The comments say that the band detection constants (background $\sigma$ 0.2×short side, threshold 4.5 MAD, NMS 6 px, line verification 0.5/0.3/17/0.55/0.05) were determined by quantitative evaluation on a particular set of measured images. Should the conditions of those images (detector, resolution, binning) be documented?
- **I-3** $\lambda/(2d)$ in the width likelihood is the kinematical Bragg angle. The difference between the band edge position (Kikuchi line position) and the apparent width of the dynamical band (also related to the discussion of $A(E)$ in [ebsd-pattern-composition.md](ebsd-pattern-composition.md) §4.2).
- **I-4** Dictionary matching projects a single MC-weighted composite plane (or a single slice). The effect on matching accuracy of the variation of weights with position on the detector (the approximation of [ebsd-pattern-composition.md](ebsd-pattern-composition.md) §5).
- **I-5** The soft bounds, the multi-start widths, and the per-stage convergence criteria of geometry calibration are empirical values from the author's instrument (comments). Validity for other detectors.
- **I-6** The basis for the coefficients in the zone-axis scoring formula ($+0.5W/W_{\max}$, $+2n$).
- **I-7** In the Radon equivalence test ($Q$ preserves the set of catalog directions), the catalog identifies ±, so the reduction may be by a group larger than the crystal's point group (Laue group + operations that are kinematically indistinguishable). Is this as intended for polar crystals?
