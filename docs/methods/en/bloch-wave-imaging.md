Draft written by Claude from the code (2026-11-05). Not yet reviewed by the author.

# Image simulation with the Bloch-wave method — HRTEM, STEM (ADF/TDS), potential images

Scope: `GetHRTEMImage` (L3709–3834), `RunSTEM` / `StemDoWork` (L2478–3634) and `GetPotentialImage` (L3640–3690) in `Diffraction/BetheMethod.cs`.
The eigenproblem, beam selection and `getU` conventions are in [bloch-wave.md](bloch-wave.md) (nm, $k = 1/\lambda$, plane waves $e^{+2\pi i\mathbf k\cdot\mathbf r}$, angles in radians).
For STEM-EDX (ionisation channels: matrix elements and normalisation) see [alchemi-ionization.md](alchemi-ionization.md) and `Diffraction/Ionization/README.md`; only the parts EDX shares with TDS are mentioned here.

---

## 1. HRTEM `GetHRTEMImage`

### 1.1 Inputs

| Argument | Unit / meaning (from the doc comments) |
|---|---|
| `BlochNum` | number of Bloch waves (beams) N |
| `AccVol` | kV |
| `rot` | crystal orientation |
| `aperture (R, X, Y)` | objective-aperture semi-angle and shift. The doc says "mrad", but `ExtractInsideBeams` uses `2 sin(R/2)/λ`, so **the value reaching this function is in radians** (→ item 1) |
| `size`, `res` | pixels, nm/pixel |
| `cs` | spherical aberration $C_s$ [nm] |
| `beta` | illumination semi-angle $\beta$ [rad] (spatial coherence) |
| `delta` | focal spread $\delta = C_c\,\Delta V/V$ [nm] (temporal coherence) |
| `thicknesses`, `defocusses` | arrays in nm |
| `quasiMode` | true: quasi-coherent (default) / false: TCC |
| `native` | true and native available → C++ grid solver |

260610Cl: the wavelength is computed from the argument `AccVol`, not the field `AccVoltage` (previously λ = ∞ → all-NaN images on a fresh instance).

### 1.2 Procedure

1. For each thickness, `GetDifractedBeamAmpriltudes(BlochNum, AccVol, rot, t)` ([bloch-wave.md](bloch-wave.md) §8; the eigendecomposition is reused after the first call) gives exit amplitudes $\psi_\mathbf g$.
2. `ExtractInsideBeams`: $|\mathbf g_{xy} - \mathbf c| < 2\sin(R/2)/\lambda$, $\mathbf c = (2\sin(X/2)/\lambda,\ 2\sin(Y/2)/\lambda)$; all beams if $R = \infty$.
3. The incident tilt $\mathbf k_{xy}$ is fixed at 0 (`var k = new PointD(0, 0)`), so $\mathbf q = \mathbf g_{xy}$.

### 1.3 Aberration and envelope functions

$$
\chi(q) = \pi\lambda q^2\Big(\Delta f + \tfrac12 C_s\lambda^2 q^2\Big),\qquad \text{transfer} = e^{-i\chi(q)}
$$

($\Delta f$ = `defocus`. With this sign, Scherzer defocus for $C_s > 0$ is on the $\Delta f < 0$ side. → item 2)

**Quasi-coherent** (envelopes applied independently per beam):

$$
L(q) = e^{-i\chi(q)}\;
\exp\!\big[-\pi^2\beta^2 q^2(\Delta f + \lambda^2 C_s q^2)^2\big]\;
\exp\!\big[-\tfrac12\pi^2\lambda^2\delta^2 q^4\big]
$$

$$
I(\mathbf r) = \Big|\sum_\mathbf g \psi_\mathbf g\,L(q_\mathbf g)\,e^{2\pi i\mathbf g\cdot\mathbf r}\Big|^2
$$

The second factor is spatial and the third temporal coherence (the comments at L3734–3735 have them the other way round → [bloch-wave.md](bloch-wave.md) Q3).

**TCC** (transmission cross-coefficient), over pairs $g \le h$:

$$
T(\mathbf q_1,\mathbf q_2) = e^{-i\chi(q_1)}e^{+i\chi(q_2)}
\exp\!\big[-\pi^2\beta^2\,|\Delta f(\mathbf q_1-\mathbf q_2) + \lambda^2C_s(q_1^2\mathbf q_1 - q_2^2\mathbf q_2)|^2\big]
\exp\!\big[-\tfrac12\pi^2\lambda^2\delta^2(q_1^2-q_2^2)^2\big]
$$

Coefficient $|\psi_g|^2$ for $g = h$ and $2\psi_g\psi_h^*$ for $g<h$; $I(\mathbf r) = \big|\mathrm{Re}\sum_{g\le h} (\ldots)\,T\,e^{2\pi i(\mathbf g-\mathbf h)\cdot\mathbf r}\big|$. Since $T(\mathbf q_2,\mathbf q_1) = T(\mathbf q_1,\mathbf q_2)^*$, taking the real part over half the pairs equals the full double sum.

### 1.4 Pixel coordinates

$\mathbf r = \big((x - c_X)\,\text{res} + s_X,\ -(y - c_Y)\,\text{res} + s_Y\big)$ with $c$ = pixel count/2 and $\mathbf s = R\,(\mathbf a+\mathbf b+\mathbf c)/2$ (half the body diagonal rotated by `Crystal.RotationMatrix`). The native version receives the same grid as start + step and advances along rows by phase-rotation recurrence, re-anchoring with `Exp` every 256 pixels (260610Cl, difference ~1e-13). Images are not normalised (raw $|\cdot|^2$ for quasi-coherent, raw $|\mathrm{Re}|$ for TCC).

---

## 2. STEM `RunSTEM` → `StemDoWork`

### 2.1 Inputs

| Argument | Unit / meaning |
|---|---|
| `maxNumOfBloch` | N |
| `voltage` | kV |
| `cs` | $C_s$ [nm] |
| `delta` | focal spread $\delta$ [nm] |
| `sliceThickness` | target depth step for the inelastic (TDS/EDX) integral [nm] |
| `imageSize`, `resolution` | pixels, nm/pixel |
| `sourceSize` | source size: **half width at half maximum** of the Gaussian blur applied to the image [nm] (`GaussianBlurFast(…, sourceSize/resolution)`, parameter name `hwhm`) |
| `baseRotation`, `thicknesses`, `defocusses` | orientation, nm, nm |
| `beamDirections` | incident directions, treated as a **square grid** (side $\sqrt{\text{Length}}$); only the inscribed circle is computed |
| `convergenceAngle` | convergence semi-angle $\alpha$ [rad] |
| `detAngleInner`, `detAngleOuter` | ADF detector inner/outer angles [rad] |
| `solver`, `thread` | [bloch-wave.md](bloch-wave.md) §6 |
| `ionizations` | STEM-EDX channels (empty → no EDX) |

`RunSTEM` returns immediately if busy (260711Cl). Arrays are copied before being handed to the worker (260801Cl). The result is published once on success as `ResultStem` (`StemSimulationResult`); `ResultSTEM` is a compatibility view.

### 2.2 Bloch-wave part

- `Find_gVectors` at the central direction, then **re-sorted by $|\mathbf g|^2$** (`OrderBy(e => e.Vec.Length2)`). Incident vector $\psi_0 = (-1, 0, \ldots)$. The potential matrix (local form, full-solid-angle TDS absorption) is built once.
- For each incident direction $\mathbf k$ (inside the circle): `reset_gVectors` → `getEigenMatrix` → EVD. Transmission coefficients

  $$\mathrm{tc}_\mathbf g(\mathbf k, t) = e^{2\pi i (P_\mathbf g/2)\,t}\sum_j C^{(j)}_\mathbf g e^{2\pi i\gamma^{(j)}t}\alpha_j$$

  are kept for the user thicknesses $t$ and all slice depths $z$.
- The 260711Cl "all-slice materialisation" path (native Eigen and favourable memory estimate; env `RECIPRO_STEM_ALLSLICE` forces it) builds tc for all depths with one GEMM right after the EVD and does not keep eigenvectors. Otherwise eigenvectors pre-multiplied by α are kept and tc is regenerated per slice by `GenerateTC1B`. Same formula either way.

### 2.3 Lens, aperture and detector functions (L2730–2756)

$$
L(\mathbf k, \mathbf k+\mathbf q; \Delta f) = \exp\!\Big\{-i\big[\chi(k) - \chi(|\mathbf k+\mathbf q|)\big]\Big\}\,
\exp\!\Big\{-\big[\pi\lambda\delta\,(|\mathbf k+\mathbf q|^2 - k^2)\big]^2\Big\}
$$

$\chi$ as in §1.3 (comment: see Mr. Otsuka's note "瀬戸先生・質問3への回答(瀬戸追記).docx"). The temporal exponent lacks the $\tfrac12$ used in HRTEM (→ [bloch-wave.md](bloch-wave.md) Q2). Spatial coherence is not in the lens function; it enters as the final real-space Gaussian blur (`sourceSize`).

- Condenser aperture: $A(\mathbf k) = [\,|\mathbf k_{xy}|^2 \le (\sin\alpha/\lambda)^2\,]$
- Detector: $D(\mathbf k_f) = [\,(\sin\theta_\mathrm{in}/\lambda)^2 \le |\mathbf k_{f,xy}|^2 \le (\sin\theta_\mathrm{out}/\lambda)^2\,]$

$\mathbf k_{xy}$ is the XY component of `getVecK0` (the surface-parallel component equals the vacuum one).

### 2.4 Elastic (coherent ADF)

- Set `qList`: beam differences $\mathbf g_i - \mathbf g_j$ for which some in-circle direction $\mathbf k$ satisfies $A(\mathbf k)$ and $A(\mathbf k+\mathbf q)$; sorted by $|\mathbf q|^2$.
- Amplitudes at $\mathbf k+\mathbf q$ are **bilinearly interpolated** from the 4 neighbouring grid directions (`list[kIndex]` stores indices `n0..n3` and weights `r`). The map to grid coordinates is $d = q_x\cdot(r_\mathrm{pix}-0.5)/\max_k k_x + (r_\mathrm{pix}-0.5)$ (y sign-flipped), assuming `beamDirections` is a square grid uniform in $\mathbf k_{xy}$. Entries are dropped unless all 4 neighbours are inside the circle.

$$
I_\mathrm{el}(\mathbf q; t, \Delta f) = \sum_{\mathbf k} L(\mathbf k,\mathbf k+\mathbf q;\Delta f)
\sum_{\substack{\mathbf g:\ D(\mathbf k+\mathbf g)\\ \mathbf g-\mathbf q\in\text{Beams}}}
\mathrm{tc}_\mathbf g(\mathbf k,t)\,\Big[\sum_{m=0}^{3} r_m\,\mathrm{tc}_{\mathbf g-\mathbf q}(\mathbf k_m,t)\Big]^*
$$

(interference of the two components exiting in the same direction $\mathbf k+\mathbf g = (\mathbf k+\mathbf q)+(\mathbf g-\mathbf q)$.)

### 2.5 TDS (incoherent ADF)

The flux removed by the absorptive potential and scattered into the detector range is integrated over depth:

$$
I_\mathrm{TDS}(\mathbf q; t, \Delta f) = \sum_{z_s \le t} \frac{2\pi}{k_\mathrm{vac}}\,\Delta z\;
\sum_\mathbf k L(\mathbf k,\mathbf k+\mathbf q;\Delta f)\;
\boldsymbol\psi(\mathbf k+\mathbf q, z_s)^\dagger\,\mathsf U'_\mathbf q\,\boldsymbol\psi(\mathbf k, z_s)
$$

$$
(\mathsf U'_\mathbf q)_{ij} = \overline{\texttt{getU}(\mathbf q + \mathbf g_i - \mathbf g_j;\ \theta_\mathrm{in},\theta_\mathrm{out}).\texttt{Imag}}
$$

- Annular `getU` (local form, [bloch-wave.md](bloch-wave.md) §2.4): $f'$ integrated only over detector angles [inner, outer]. The non-local form is commented out.
- Complex conjugate (L2895, comment: "taking the conjugate somehow works nicely"): `getU` uses the phase $e^{+2\pi i\mathbf G\cdot\mathbf r}$, so `getU(G)` is the textbook $U'_{-\mathbf G}$; because $f'$ is real, its conjugate is $U'_{\mathbf G}$. The Fourier component coupling $\psi^*_{\mathbf k+\mathbf q,i}\,\psi_{\mathbf k,j}$ is $U'_\mathbf G$ (textbook sign) with $\mathbf G = \mathbf q + \mathbf g_i - \mathbf g_j$, so the conjugate is the correct sign fix (Claude's analysis; it appears to be necessary rather than empirical).
- The factor $2\pi/k_\mathrm{vac}$ matches the absorption rate $d|\psi|^2/dz = -(2\pi/K)U'|\psi|^2$ ([bloch-wave.md](bloch-wave.md) §2.6).
- **Depth steps**: interval $[T_{t-1}, T_t]$ ($T_{-1} = 0$) split into $n = \max(1, \lfloor \Delta T/\text{sliceThickness}\rfloor)$ equal slices, evaluated at the **right end** of each, $z_s = T_{t-1} + e\,\Delta z$ ($e = 1..n$). Right Riemann sum, error $O(\Delta z)$. Accumulated over thickness ($I(T_t) = I(T_{t-1}) + \ldots$).
- Code paths (native `STEM_InelasticQ`, per-q U streaming `RECIPRO_STEM_USTREAM`, managed fallback) differ only in summation order; the formula is the same (260610Cl, 260711Cl comments).

### 2.6 Real-space images

$$
I_X(\mathbf R) = \frac{1}{r_\mathrm{pix}^2}\,\mathrm{Re}\sum_\mathbf q I_X(\mathbf q)\,e^{2\pi i\mathbf q\cdot\mathbf R},\qquad
\mathbf R = \big(\text{res}(x-c_X),\ -\text{res}(y-c_Y)\big) + \mathbf s_{xy}
$$

- $r_\mathrm{pix}$ = radius of the direction grid in pixels. $r_\mathrm{pix}^2$ is the "existing convention for relative intensity display", not an absolute normalisation (260802Cl: EDX was changed to divide by the number of q = 0 entries = number of in-aperture directions; ADF was left unchanged for regression compatibility).
- **Real part** (260802Cl; previously `.Magnitude`). The comment's reasoning: $I(-\mathbf q) = I(\mathbf q)^*$ is broken by $O(h^2)$ through the asymmetric interpolation, the resulting imaginary part is numerical error, and orthogonal projection removes it. The Hermitian residual (`MeasureHermitian`), the largest discarded imaginary part and the most negative pixel are recorded as diagnostics only; ADF never throws (EDX throws above a tolerance and symmetrises ±q).
- Then, if `sourceSize > 0`, a Gaussian blur (HWHM = sourceSize/resolution pixels), negative values clamped to 0, and `Both = Ela + TDS`.

### 2.7 Shared with STEM-EDX

EDX uses the same transmission coefficients, q entries, lens function and accumulation as TDS, replacing $\mathsf U'_\mathbf q$ with `FillIonizationUq` ($U_\mathrm{ion}$, negative phase built directly, cache separate from `uDictionary`) (`AccumulateInelasticQ`, 260802Cl). EDX requires the native helpers (`NotSupportedException` on the managed path).

---

## 3. Potential image `GetPotentialImage` (L3640–3690)

Synthesises `Ureal` ($\hat U_\mathbf g$) and `Uimag` ($U'_\mathbf g$) of the given beams in real space:

$$
S_\mathrm{cry}(\mathbf r) = \sum_\mathbf g \hat U_\mathbf g\,e^{2\pi i\mathbf g_{xy}\cdot\mathbf r},\qquad
S_\mathrm{th}(\mathbf r) = \sum_\mathbf g U'_\mathbf g\,e^{2\pi i\mathbf g_{xy}\cdot\mathbf r}
$$

Four output images: with `phase = true` $|S_\mathrm{cry}|, \arg S_\mathrm{cry}, |S_\mathrm{th}|, \arg S_\mathrm{th}$; otherwise real and imaginary parts. Coordinates

$$
\mathbf r = \big(-(x-c_X)\,\text{res} + s_X,\ -(H-1-y-c_Y)\,\text{res} + s_Y\big)
$$

are sign-reversed in both X and Y relative to HRTEM/STEM. Combined with `getU` returning the textbook $U_{-\mathbf g}$ ([bloch-wave.md](bloch-wave.md) §2.1), $S(\mathbf r) = U(-\mathbf r)$, so the reversed coordinates yield $U$ in the correct orientation (Claude's analysis). However Y uses $(H-1-y)$, a 1-pixel offset from HRTEM/STEM ($y$) (→ item 4). Unit nm⁻² ($U$ itself).

---

## Items to confirm

1. The `GetHRTEMImage` doc says the aperture semi-angle is in "mrad", but `ExtractInsideBeams` uses `2 sin(R/2)/λ` (radians). Does the caller (ReciPro) convert mrad→rad?
2. Defocus sign: $\chi = \pi\lambda q^2(\Delta f + \tfrac12 C_s\lambda^2q^2)$, so negative $\Delta f$ is underfocus. Does the UI sign match?
3. Missing $\tfrac12$ in the STEM temporal envelope ([bloch-wave.md](bloch-wave.md) Q2). Fix the definition of `delta` (rms or 1/e half width).
4. `GetPotentialImage` uses $(H-1-y)$ for Y, a 1-pixel offset relative to HRTEM/STEM ($y$). Does this matter for overlays?
5. The STEM TDS depth integral is a right Riemann sum (L2578). The midpoint rule would give $O(\Delta z^2)$; is the default `sliceThickness` small enough?
6. The bilinear interpolation assumes a square grid uniform in $\mathbf k_{xy}$ (`maxK = max k_x`, `coeff2`). Does the caller (ReciPro) generate `beamDirections` that way (an equal-angle grid on the sphere is non-uniform by the sine)?
7. The ADF normalisation $1/r_\mathrm{pix}^2$ differs from EDX (number of q = 0 entries) by a ratio of about $\pi/1.05^2$ (260802Cl comment). Are ADF and EDX ever compared on the same scale?
8. HRTEM fixes the incident tilt $\mathbf k_{xy} = 0$ (L3716). Is tilted illumination intentionally unsupported?
9. The STEM direction grid extends beyond the convergence angle (L2791 comment: "not the convergence angle, note (5% larger)"). The grid mapping uses `maxK` (largest computed $k_x$) while the aperture test uses $\sin\alpha/\lambda$. Is this consistent with how the caller builds the grid?

(9 items to confirm. Implementation questions are listed as Q1–Q3 in [bloch-wave.md](bloch-wave.md).)
