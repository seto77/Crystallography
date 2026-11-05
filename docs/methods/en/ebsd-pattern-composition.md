This is a draft written by Claude from the code (2026-11-05). Not yet reviewed by the author.

# EBSD pattern composition: master pattern × MC weights

Target commit: `9ffcf47` (main). Files: `EBSD/EbsdPatternComposer.cs`, `EBSD/EbsdMonteCarloDistribution.cs` (`ComposeGlobalWeightedPattern`). Prerequisites: [ebsd-geometry.md](ebsd-geometry.md) (pixel → direction), [ebsd-master-pattern.md](ebsd-master-pattern.md) (meaning of the planes $M_{e,d}$), [ebsd-monte-carlo.md](ebsd-monte-carlo.md) (MC distribution).

`EbsdPatternComposer` is the display composer moved out of `FormEBSD`, and does not depend on the GUI. Its differences from `EbsdPatternProjector` (for matching, single slice) are as in the doc: this one handles a raster of the whole field of view, MC-weighted composition over all slices, and the current orientation of the crystal.

---

## 1. Notation

- $M_{e,d}(\hat{\mathbf d})$: the value of the plane at energy $E_e$ and depth (path length) $t_d$, interpolated along the pixel's line of sight (viewing direction) $\hat{\mathbf d}$ (crystal frame). Interpreted as the cumulative response $M(t_d)=\int_0^{t_d}R\,dz$ (doc of `DepthIntervals`).
- $\Delta M_{e,d} = \max\bigl(0,\ M_{e,d}-M_{e,d-1}\bigr)$ ($M_{e,-1}=0$). $\Delta t_d$ from `ComputeDepthIntervals`.
- $W_{e,d}(i)$: the (energy × depth) weight of pixel $i$ (§3).
- $\bar M_{e,d}$: the **full-sphere** directional mean of the plane (mean of the means over the $+Z$ and $-Z$ hemispheres; §4.3).
- $A_e$: the loss-dependent contrast factor (§4.1). $f_A(i)$: the fraction of electrons whose source is in the amorphous layer (§4.2).

---

## 2. Single-slice models

After building the lookup table (`BuildLookupTable`, independent of energy and depth), the specified slice is written.

| Model | Method | Pixel value |
|---|---|---|
| 0 current | `ApplySingleSliceModel0` | $M_{e,d}$ (= model 1 with factor 1) |
| 1 globally normalized | `ApplySingleSliceModel1` | $c_k\,M_{e,d}$, $c_k = 1/\sum_{\rm cells}(M^{+}_k + M^{-}_k)$ (sum over all cells of both hemispheres; 0 if ≤ 0) |
| 2 differential | `ApplySingleSliceModel2` | $\max(0,\,M_{e,d}-M_{e,d-1})$ (**not divided by $\Delta t$**) |

The square-grid Rosca–Lambert grid is equal-area, so the cell sum in $c_k$ is proportional to the full-sphere integral. If either plane has fewer than `GridSize²` entries, all pixels are 0.

---

## 3. Per-pixel weights $W_{e,d}(i)$ (`EvaluatePathLengthWeights`, depth mapping A2)

1. The pixel's **exit direction** (sample frame) is built with `ExitRay` ([ebsd-geometry.md](ebsd-geometry.md) §3.6). Pixels with $o_z\le0$ (below the sample surface) are 0.
2. $\mu = o_z/|\mathbf o|$, bin coordinates $(b_x,b_y)$ = `DirectionToBinCoords(o)` (same map as the MC distribution).
3. 4 × 4 taps of a **cubic B-spline** (`BinSplineTaps`, $t = b-\lfloor b\rfloor$):

$$ q(t) = \tfrac16\bigl[(1-t)^3,\ 3t^3-6t^2+4,\ -3t^3+3t^2+3t+1,\ t^3\bigr],\quad \text{taps } \lfloor b\rfloor-1,\dots,\lfloor b\rfloor+2 \text{ (only the indices are clamped at the edges)}. $$

   $\sum q = 1$, $\sum q\cdot(\text{tap position}) = b$ (reproduces first order exactly), $C^2$. It is an approximation, not an interpolation (the value at a bin center is $(v_{i-1}+4v_i+v_{i+1})/6$). $b$ is only limited to $[-3, n_b+2]$, not clamped.
4. Weight of tap $t$: $\omega_t = c_t\,s_b$, $c_t = q_x q_y$. For model 2, $s_b = F_b$ (electron fraction); for model 0/1, $s_b = [F_b>0]$. Taps with $\omega_t\le0$ are skipped.

$$ N(E_e) = \sum_t \omega_t G_{b_t}(E_e),\qquad \lambda(E_e) = \frac{\sum_t\omega_t G_{b_t}(E_e)\,\lambda_{b_t}(E_e)}{N(E_e)}\ \ (10^6\ \mathrm{nm} \text{ if } N=0). $$

5. $W_{e,d}$ from `FillPathLengthWeights(N, λ, μ, depths, Δt, sliceMass)` ([ebsd-monte-carlo.md](ebsd-monte-carlo.md) §5.5). Model 2 uses `sliceMass = true` (interval mass, total $\sum_e N_e$); model 0/1 uses `false` (density × interval width).
6. Amorphous fraction $f_A = \sum_t\omega_t f_{A,b_t}/\sum_t\omega_t$.

Note in the doc: the old version interpolated "weight arrays built with the bin-center $\mu$", but $\chi$ varies by about 11° even within one 16×16 exit-hemisphere bin, so **the parameters ($G$, $\lambda$, $F$) are interpolated first and then converted with the pixel's $\mu$**. Even at $\mu=1$ this does not agree with the old version (nonlinear).

---

## 4. MC-weighted models (`ApplyWeightedModel0/1/2`)

Precondition check: `dist.MatchesGridOf(mp)` (the energy and depth grid values agree to a relative $10^{-9}$). If they do not, an exception is thrown (260727Cl; because a mismatched grid would silently look up the weights of a different slice).

### 4.1 Basic form with A(E) disabled and no amorphous layer

$$
\begin{aligned}
\text{model 0:}\quad V_i &= \sum_{e,d} W_{e,d}(i)\,M_{e,d}(\hat{\mathbf d}_i),\\
\text{model 1:}\quad V_i &= \sum_{e,d} W_{e,d}(i)\,c_k\,M_{e,d}(\hat{\mathbf d}_i),\\
\text{model 2:}\quad V_i &= \Omega_i\sum_{e,d} W_{e,d}(i)\,\frac{\Delta M_{e,d}(\hat{\mathbf d}_i)}{\Delta t_d},\qquad \Omega_i = \cos^3\alpha_i .
\end{aligned}
$$

- In model 2, $\Delta M/\Delta t$ is the interval-mean response $\bar R$ (the MC-side weight is an interval mass, so it is not divided; 260726Cl).
- $\Omega_i$ is the $\cos^3\alpha$ part of the pixel's solid angle $d\Omega/dA = D/|\mathbf r|^3 = \cos^3\alpha/D^2$ ($\cos\alpha = D/|\mathbf r|$, $D$ = `ExitRay.PlaneDistance`). The exit-hemisphere bins have equal solid angle, so the "per solid angle" weight is multiplied by the pixel's solid angle. The constant $1/D^2$ is dropped (the display auto-stretches, and ZNCC is invariant to constant factors). **Not applied to model 0/1**.

### 4.2 Loss-dependent contrast factor $A(E)$ (260920Cl onward)

$$ A_e = \exp\!\Bigl(-\frac{\max(0, E_0-E_e)}{E_c}\Bigr),\qquad E_0 = \texttt{BeamEnergyKeV}\ (\max_e E_e \text{ if NaN}),\quad E_c = \texttt{CoherenceLossDecayKeV}. $$

Disabled if $E_c$ is NaN or ≤ 0 (all slices at full contrast). If all values are ≥ $1-10^{-12}$, `null` (same as disabled). The formula exists in only one place, `CoherenceFactors` (on the premise that display composition and the ZNCC objective use the same weights).

Composition when enabled (`CohAccum`, 260921Cl):

$$
V^{\rm coh}_i = \underbrace{\sum W\bar M}_{D_c}
\;+\; \bar A\,\Bigl(\sum W\Bigr)\,\frac{\sum W A\,I - \sum W A\bar M}{\sum W A},
$$

where $I$ is $M$ / $c_kM$ / $\Delta M/\Delta t$ depending on the model, $W$ is $W_{e,d}(i)$ ($W c_k$ for model 1), and the sums run over $(e,d)$. That is, $V = \sum W\bar M + \bar A\,(\sum W)\,(\langle I\rangle_A - \langle\bar M\rangle_A)$, where $\langle\cdot\rangle_A$ is the normalized mean weighted by $WA$. The first term is the brightness with $A$ disabled (the unmodulated component), and the second term scales the contrast by $\bar A$ uniformly over the whole pattern.

$\bar A$ (`MeanCoherentFraction`) is a single scalar independent of the pixel, obtained by sampling the **detector plane** at 48 × 48 points and multiplying each point's weight by the solid angle $\cos^3\alpha$:

$$ \bar A = \frac{\sum_{\rm pts}\sum_{e,d}\cos^3\alpha\;W_{e,d}\,A_e}{\sum_{\rm pts}\sum_{e,d}\cos^3\alpha\;W_{e,d}} $$

($W c_k$ for model 1). Sums are taken per row and then over rows in order, so the result does not depend on the parallel ordering. The range is `EbsdRasterView.HalfWidth/HalfHeight` (the physical dimensions of the detector), independent of the display pan and zoom.

For distributions with `energyHistogram = true`, the per-pixel effective $A$ ($G_A/G$ interpolated with the same 4 × 4 taps, `ComputeElectronWiseCoherence`; 260927Cl) is used instead of the nodal $A_e$.

History in the doc (summary): the old version added what was lost through $A$ back as a single scalar pedestal, but because $\sum W A$ carries an 8×8 imprint per pixel, broad dark lines at the bin spacing appeared after background flattening. The formula above removes the brightness imprint. The in-plane variation of the contrast (the difference in the fraction of low-loss electrons with take-off angle) is deliberately discarded.

What the doc says about the default value of $E_c$: an empirical value calibrated on a single measured Si 20 kV pattern (0.7 → 0.8 keV), not a constant independent of material and voltage. It also says the result for a second material (Botallackite) cannot establish a material dependence within the error (numerical details in the doc of `CoherenceLossDecayKeV`).

### 4.3 Amorphous layer (260919Cl)

Electrons with sources within the layer contribute as a uniform component without Kikuchi modulation:

$$ V_i = \Omega_i\Bigl[(1-f_A)\,V^{\rm coh}_i + f_A\sum_{e,d}W_{e,d}\,\bar M^{\rm ref}_e\Bigr] \quad(\text{for model 0/1, }\Omega_i=1\text{; for model 1, } W \text{ is multiplied by } c_k). $$

- $\bar M^{\rm ref}_e$ (`sphereMeans`): the directional mean per plane (excluding cells with value 0; for model 2, the mean of $\Delta M/\Delta t$), averaged over depth with the all-bin depth weights (`GlobalDepthWeights` / `GlobalDepthSliceWeights`) into a "reference value per energy", then averaged over $+Z$ and $-Z$ (`CollapseToEnergyReference`). A per-hemisphere value would produce a step at the equator, so the full-sphere mean is used.
- In contrast, $\bar M_{e,d}$ of §4.2 (`planeMeans`) is the full-sphere mean per plane, not collapsed onto energy (because the $A$ pedestal must cancel the coherent term one-to-one in $(e,d)$).
- Plane means are cached when `(MasterPattern, Distribution, differential)` are the same.

### 4.4 Background (optional output of model 2)

If passed to `ApplyWeightedModel2(..., background)`, it also writes $B_i = \Omega_i[(1-f_A)\sum W\bar M + f_A\sum W\bar M^{\rm ref}]$ (the Kikuchi modulation removed).

### 4.5 Diagnostic output

`ComputeEffectiveCoherenceModel2`: returns the coefficient of $A$ that the model 2 composition applies at pixel $i$ and node $e$, together with $\bar A$ (so that the external reference implementation `tools/EbsdProfileFit` can reproduce the same discretization; the tool is not in this repository).

---

## 5. Reduction to a single plane for matching (`EbsdMonteCarloDistribution.ComposeGlobalWeightedPattern`)

For ZNCC matching and geometry calibration, weights are not built per pixel; the weights are averaged into one set, the planes are reduced to a single plane, and it is projected with `EbsdPatternProjector`:

$$ \tilde w_{e,d} = \frac{A_e\,\bar W_{e,d}}{\sum A\bar W},\qquad
\tilde M^{\pm}(\text{cell}) = \sum_{e,d}\frac{\tilde w_{e,d}}{\Delta t_d}\,\max\bigl(0, M^{\pm}_{e,d}-M^{\pm}_{e,d-1}\bigr)\quad(\text{skipped if } \tilde w<10^{-15}). $$

- $\bar W$ is the argument `sliceWeights` (recommended: `EbsdPatternComposer.ComputeDetectorAverageSliceWeights` = pixel average over the detector image, subsampled with `stride` 4, $\cos^3\alpha$ applied to each pixel, pixels below the horizon counted in the denominator as 0). If null, the sum of `BinAbsoluteSliceWeights` over all bins (average over the whole exit hemisphere).
- The $(1-A)$ pedestal is not added (ZNCC is invariant to constant offsets and constant factors).
- Note in the doc: this is an approximation that "averages the weights first, then reduces to one plane and projects"; the final calibration of $E_c$ should be done with per-pixel composition.

---

## 6. Implementation questions

### P-I1 The weights of $\bar A$ and of the pixel values are not consistent for model 0/1

`MeanCoherentFraction` multiplies the weight of each sample point by $\cos^3\alpha$ even for model 0/1, whereas `ApplyWeightedModel0/1` does not multiply the pixel values by $\Omega_i$. The doc defines "W is each model's pixel weight × the pixel's solid angle", so this may be read as intentional, but whereas for model 2 it agrees with "the pixel average of the display composition weights", for model 0/1 it does not.

### P-I2 Single-slice model 2 does not divide by $\Delta t$

`ApplySingleSliceModel2` writes $\max(0, M_d - M_{d-1})$ as is. The weighted model 2 and `ComposeGlobalWeightedPattern` apply $/\Delta t_d$ (260726Cl). Single-slice display does not compare across slices, so this is not an apparent problem, but when switching between slices on an unequally spaced depth grid and comparing them, the brightness changes in proportion to the interval width.

---

## 7. Items to confirm

- **P-1** The physical standing of model 0 (summing the cumulative $M$ with density × interval-width weights) and model 1. The doc gives only the names (current / globally normalized master), without explaining which approximation they correspond to.
- **P-2** $\Omega_i=\cos^3\alpha$ in model 2 corresponds to "the number of electrons received by the detector per pixel", but does not include the angular dependence of phosphor emission or the camera sensitivity. Is background flattening assumed when comparing with measurements?
- **P-3** The default $E_c$ of $A(E)$ (commented as 0.8 keV in the UI field), and whether recalibration is needed when `energyHistogram = true` (the doc of `EbsdMonteCarloDistribution` says "with true, E_c must be recalibrated").
- **P-4** `MeanCoherentFraction` uses fixed sampling at 48 × 48 points. Is this sufficient when the detector aspect ratio is large?
- **P-5** `BinSplineTaps` is an approximation (smoothing the bin values by about one bin), so sharp angular dependence of the exit distribution (e.g. anisotropy of the BSE yield due to channeling) cannot be represented. Is it correct to understand that this is not a problem because the MC distribution does not include crystallinity (amorphous MC)?
