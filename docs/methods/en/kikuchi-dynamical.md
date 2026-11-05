This is a draft written by Claude from the code (2026-11-05). Not yet reviewed by the author.

# Dynamical Kikuchi diffraction — EBSD master patterns and TEM Kikuchi band profiles

This document covers the two Bloch-wave (Bethe) implementations of Kikuchi diffraction.

| Path | Purpose | Entry point |
|---|---|---|
| A. EBSD master pattern | Intensity of backscattered electrons per exit direction (full sphere, crystal-fixed) | `BetheMethod.RunEBSDNew` → `ebsdNew_DoWork` (`Diffraction/BetheMethod.cs`) |
| B. 1D Kikuchi band profile | Profile across a transmission (forward) Kikuchi band (two-beam / systematic row) | `KikuchiProfileCalculator.ComputeProfile(s)` (`Diffraction/Kikuchi/`) |

- For the general form of the potential matrix `getU`, the eigenproblem `getEigenMatrix` and the beam search `Find_gVectors`, see [bloch-wave.md](bloch-wave.md). This document describes what is specific to Kikuchi/EBSD.
- Scattering factors $f^{e}$, TDS absorptive form factors $f'$ (`FactorImaginary*`) and the Debye–Waller factor: see [scattering-factors.md](scattering-factors.md).
- Expansion of the unit cell and structure-factor conventions: [structure-factor.md](structure-factor.md).
- The EBSD Monte Carlo (depth/energy distributions), projection onto the detector, the master-pattern grid (Rosca–Lambert) and its symmetry operations live in `EBSD/` (separate document). Here only the quantities exchanged with it are defined.

---

## 1. Symbols and conventions (both paths)

| Symbol | Meaning | Units / convention | Implementation |
|---|---|---|---|
| $\mathbf g,\mathbf h$ | Reciprocal-lattice vectors | nm⁻¹, **no 2π** ($\lvert\mathbf g\rvert = 1/d$) | `Beam.Vec`; `Crystal.MatrixInverseTransposed` multiplied by `BaseRotation` |
| $s^2$ | Scattering-factor argument | $s^2 = \lvert\mathbf q\rvert^2/4$ [nm⁻²] ($s=\sin\theta/\lambda$) | `s2` in `getU` and in the kernels |
| $k_\mathrm{vac}$ | Vacuum wave number | nm⁻¹ (relativistically corrected, no 2π) | `UniversalConstants.Convert.EnergyToElectronWaveNumber(kV)` |
| $\gamma_\mathrm{rel}$ | Relativistic factor | $1 + eE/(m_0c^2)$ | in `getU`, and in the σ computation of `ebsdNew_DoWork` |
| $U_\mathbf g = U^{r}_\mathbf g + iU'_\mathbf g$ | Fourier coefficient of the optical potential | nm⁻². $U_\mathbf g = \gamma_\mathrm{rel}F_\mathbf g/(\pi V)$, $V$ [nm³] | `(Real, Imag)` of `getU(...)`; `Imag` is the absorption (TDS) |
| $U_0$ | Mean-inner-potential term | $U^{r}_{000}$ | `u0 = getU(AccVoltage).Real.Real` |
| Phase convention | Atomic-position factor | $e^{+2\pi i\,\mathbf g\cdot\mathbf r_n}$ (`getU`, `CreatePhaseFactors`) | sign changed to + on 20240524 (comment in `getU`) |
| $\hat{\mathbf n}$ | Surface normal (**pointing into the sample**) | unit vector, default `Surface = (0,0,−1)` | `BetheMethod.Surface` |
| $\mathbf K_0$ | Incident wave vector inside the crystal | $\lvert\mathbf K_0\rvert^2 = k_\mathrm{vac}^2 + U_0$, tangential component conserved | `getVecK0` |
| $Q_\mathbf g,\ P_\mathbf g$ | Coefficients of the Bethe equation | $Q_\mathbf g = \lvert\mathbf K_0\rvert^2 - \lvert\mathbf K_0+\mathbf g\rvert^2$, $P_\mathbf g = 2\hat{\mathbf n}\cdot(\mathbf K_0+\mathbf g)$ | `getQ`, `getP`, `getQP` |
| $\gamma_j,\ C_{\mathbf g j}$ | Eigenvalues / eigenvectors | $\gamma$ [nm⁻¹]. Wave field $\psi_\mathbf g(z) = \sum_j C_{\mathbf g j}\alpha_j e^{2\pi i\gamma_j z}$; decays when $\operatorname{Im}\gamma_j>0$ | `eigenValues`, `eigenVectors` (column-major: `[j*bLen+g]` = $C_{\mathbf g j}$) |
| $\alpha_j$ | Excitation amplitudes | $\alpha = C^{-1}\boldsymbol\psi(0)$, $\boldsymbol\psi(0)=\mathbf e_0$ (1 only for (000)) | `PartialPivLuSolve` / `LU().Solve` |
| $z,\ t$ | Depth | nm, distance along $\hat{\mathbf n}$ | `Thicknesses` |

Since 260829Cl the eigenproblem matrix (see [bloch-wave.md](bloch-wave.md)) divides row $\mathbf g$ by $P_\mathbf g$:

$$A_{\mathbf g\mathbf h} = \frac{1}{P_\mathbf g}\Bigl[Q_\mathbf g\,\delta_{\mathbf g\mathbf h} + \tilde U_{\mathbf g\mathbf h}\Bigr],\qquad
\tilde U_{\mathbf g\mathbf h} = \begin{cases} iU'_{0} & (\mathbf g=\mathbf h)\\ U_{\mathbf h-\mathbf g}^{\texttt{getU}} & (\mathbf g\neq\mathbf h)\end{cases}$$

(`getEigenMatrix`; row = $\mathbf g$, column = $\mathbf h$, element `getU(b[col] − b[row])`. The real part $U_0$ of the diagonal is folded into the refraction of $\mathbf K_0$ and is not in the matrix.) Because `getU` sums with $e^{+2\pi i\mathbf q\cdot\mathbf r}$, $U^{\texttt{getU}}_{\mathbf h-\mathbf g}$ equals $V_{\mathbf g-\mathbf h}$ of the expansion $V(\mathbf r)=\sum_\mathbf q V_\mathbf q e^{2\pi i\mathbf q\cdot\mathbf r}$. The source and TDS matrices below use the same layout (row $\mathbf g$, column $\mathbf h$ holds $V_{\mathbf g-\mathbf h}$).

Common thickness-integration kernel (`EBSDSolverManaged`, `ComputeTDSMatrixBackground`, `KikuchiProfileCalculator.Fjj`):

$$\lambda_{jj'} = 2\pi i\,(\gamma_j - \gamma_{j'}^{*}),\qquad
F_{jj'}(t) = \int_0^t e^{\lambda_{jj'}z}\,dz = \frac{e^{\lambda_{jj'}t}-1}{\lambda_{jj'}}\quad(F=t\ \text{if}\ \lvert\lambda_{jj'}\rvert^2<10^{-30}).$$

$\lambda_{j'j}=\lambda_{jj'}^{*}$ and $F_{j'j}=F_{jj'}^{*}$.

---

## 2. EBSD master pattern (`BetheMethod.RunEBSDNew`)

### 2.1 Inputs and defaults

`RunEBSDNew(int maxNumOfBloch, double[] voltages, Matrix3D rotation, double[] thickness, Vector3DBase[] beamDirections, Solver solver = Solver.Auto, int thread = 1, bool useNonLocalAbsorption = false, bool includeTDSBackground = false, bool includeAbsorbedFluxBackground = false)`

| Argument / property | Meaning | Default / notes |
|---|---|---|
| `maxNumOfBloch` | Number of beams (Bloch waves) $N$ | `Find_gVectors` picks $N+1$ and drops the trailing tie group ([bloch-wave.md](bloch-wave.md)) |
| `voltages` | Accelerating voltages (= electron energies) [kV] | the whole calculation is redone per voltage |
| `rotation` | Crystal rotation `BaseRotation` | The EBSD side (`EBSD/EBSD.cs`) passes `Matrix3D.IdentityMatrix` → directions are given in the crystal-fixed Cartesian frame |
| `thickness` | Depths $t_k$ [nm] | returns the **cumulative** intensity over $0\ldots t_k$ for each $t_k$ (§2.8) |
| `beamDirections` | Exit directions $\hat{\mathbf d}$ (unit vectors) | Must be **two square arrays (+Z and −Z hemispheres) concatenated**, length $2G^2$, otherwise `InvalidOperationException`. Invalid hexagonal cells are `null` |
| `solver` | Eigen solver | `Auto`: `Eigen_Eigen` if native is enabled, else `Eigen_MKL`; both override the thread count with `ProcessorCount` (260712Cl). `MtxExp_*` falls to the MathNet EVD path in EBSD |
| `useNonLocalAbsorption` | Use the non-local absorptive potential $U'(\mathbf g,\mathbf h)$ | false. §2.7 |
| `includeTDSBackground` | **Replace** the local point source by the non-local source $(2\pi/k)\psi^\dagger U'_\mathrm{back}\psi$ | false. Despite the name it is not added (260919Cl). §2.5 |
| `includeAbsorbedFluxBackground` | Re-inject the flux lost to absorption as a diffuse background | false. §2.6 |
| `SourceSmearingSigmaNm` (static) | Gaussian smearing width σ [nm] of the local source | `NaN` (off). A finite value also forces `IncludeTDSBackground` to true (trial feature, 260919Cl). §2.5(c) |
| `EbsdRepresentativeDirectionBlockSize` | Block width over which neighbouring directions share a beam set | property default 3. The EBSD side uses 1 for grid width $G\le512$, 3 for $\le1024$, 5 above (`EBSD/EBSD.cs`). §2.3 |
| `UseLocalSurfacePerBeamDirection` | Local-surface flag of the old path | default false, the EBSD side also passes false. The new path always uses the local surface $\hat{\mathbf n}=-\hat{\mathbf d}$ (below) |
| `ElasticIonModel` / `AbsorptionElasticModel` (static) | Ion model of the elastic factor / elastic factor inside the absorptive integrand | `Neutral` / `TemariHybrid`. [scattering-factors.md](scattering-factors.md) |

### 2.2 Geometry — exit directions via reciprocity

For each exit direction $\hat{\mathbf d}$ the code solves the "reversed wave" entering the crystal from the detector side (reciprocity). The surface is taken as the tangent plane for that direction, with inward normal

$$\hat{\mathbf n} = -\hat{\mathbf d}$$

(`GetCrystalFixedMasterEquivalentBeamAndSurface`, 260321Ch). `getVecK0(kvac, u0, beamDirection = d̂, surface = n̂)` solves

$$\mathbf K_0 = k_\mathrm{vac}\hat{\mathbf d} + x\hat{\mathbf n},\qquad x = \sqrt{b^2+U_0}-b,\quad b = \hat{\mathbf n}\cdot k_\mathrm{vac}\hat{\mathbf d} = -k_\mathrm{vac},$$

so for $\hat{\mathbf n}=-\hat{\mathbf d}$

$$\mathbf K_0 = -\sqrt{k_\mathrm{vac}^2+U_0}\;\hat{\mathbf d}$$

(of the two roots, the one pointing into the crystal is chosen; the header of `KikuchiProfileCalculator` states the same closed form). Depth $z$ is along $\hat{\mathbf n}$, i.e. the **path length along the exit direction**, not the depth below a flat sample surface (relating it to the real depth is the EBSD side's job; the note on `MeanAbsorptionLengthNm` also says "path length").

`GetCrystalFixedMasterEquivalentBeamAndSurface` builds the rotation $R$ that maps the local tangent basis $(\hat u_l,\hat v_l)$ onto the reference basis $(\hat u_r,\hat v_r)$ and brings the reference surface and beam back with $R^{\mathsf T}$. Since `GetSurfaceTangentialAxes` satisfies $\hat u\times\hat v=\hat{\mathbf n}$, the return value is always $(\hat{\mathbf d},\,-\hat{\mathbf d})$ up to rounding (§4, item 3).

### 2.3 Beam-set selection and thinning of directions

1. **Symmetry representatives** (`CreateMasterPatternSymmetryRepresentativeDirectionMapping`, 260327Ch / 260331Cl)  
   Each direction is mapped by the exact symmetry operations of the direction grid (`MasterPattern.GetMasterPatternSquareSymmetryOperations` / `GetMasterPatternHexSymmetryOperations`, implemented in `EBSD/EBSD.cs`); the smallest linear index in the orbit is the representative. Only representatives are computed dynamically; results are copied to the rest of the orbit at the end. On the hexagonal grid (`MasterPattern.ShouldUseHexGrid`) cells not satisfying $\lvert u\rvert,\lvert v\rvert,\lvert u+v\rvert\le N$ are marked `-1` and skipped.
2. **Sharing the beam set within a block** (`EbsdRepresentativeDirectionBlockSize` = $b$)  
   The grid coordinates $(w,h)$ of a representative are divided into $b\times b$ blocks and `Find_gVectors` is called once at the block centre $(\lfloor w/b\rfloor b + \lfloor b/2\rfloor,\ \ldots)$ (clamped to $G-1$). Shared: **beam set, potential matrix and atomic phase factors**. $\mathbf K_0$, $P_\mathbf g$, $Q_\mathbf g$ are recomputed per direction by `reset_gVectors` (bug fix 260321Cl: the old code reused the representative's $\mathbf K_0$ for all directions). $b=1$ means no approximation.
3. **Keep forward beams only**: per direction, only beams with $P_\mathbf g>0$ are kept (order preserved, original indices in `gMap`), and the potential matrix is sliced from the precomputed shared matrix (`getEigenMatrix(bLen, beams, gMap, baseLen, basePotential, ...)`, 260602Cl; the old path had a length mismatch and silently recomputed the local form, losing non-local absorption).
4. Incident condition $\boldsymbol\psi(0)=\mathbf e_0$. $\mathbf e_0$ is the first beam after `Find_gVectors` sorts by Rating $=\lvert\mathbf g\rvert Q_\mathbf g^2$ (= (000), since $Q_{000}=0$).

### 2.4 Eigenproblem and excitation amplitudes

- Native (`NativeWrapper.Enabled` → `EigenEnabled`) and `Solver.Eigen_Eigen`: `NativeWrapper.EigenSolver`, then $\alpha$ from `NativeWrapper.PartialPivLuSolve(C, e0)`.
- Otherwise: MathNet `DenseMatrix.Evd(Symmetricity.Asymmetric)` and `LU().Solve`.
- Work arrays are reused across directions with `ThreadLocal` + `ArrayPool` (260321Ch).

### 2.5 Backscattering source and intensity

#### (a) Local point source (default)

Using the wave field at atom $n$ (expanded to symmetry-equivalent positions; occupancy goes into $\sigma_n$), $\psi(\mathbf r_n,z) = \sum_\mathbf g \psi_\mathbf g(z)e^{2\pi i\mathbf g\cdot\mathbf r_n}$, the rate of backscattered electrons is

$$\frac{dI}{dz} = \sum_n \sigma_n\,\lvert\psi(\mathbf r_n,z)\rvert^2,\qquad I(t) = \int_0^t\frac{dI}{dz}dz.$$

Backscattering weights (made physical in 260919Cl; previously the dimensionless $Z^{1.7}\cdot\mathrm{Occ}$, off from the non-local term by 7–8 orders):

$$\sigma_n = \frac{2\pi}{k_\mathrm{vac}}\cdot\frac{\gamma_\mathrm{rel}}{\pi V}\cdot\mathrm{Occ}_n\cdot f'_{n,\mathrm{back}}(0)\quad[\mathrm{nm}^{-1}]$$

- $f'_{n,\mathrm{back}}(0)$ = `ElectronScatteringPeng[Z][sub].FactorImaginaryAnnular(kV, 0, m, π/2, π, nTheta = 30, nPhi = 12)` (TDS absorptive form factor integrated over the back hemisphere $\theta\in[\pi/2,\pi]$, $\mathbf g=\mathbf h=0$). `sub` is 0 when `ElasticIonModel == Neutral`.
- $m = \max(B_\mathrm{eff},\ 10^{-3}\ \mathrm{nm^2})$; $B_\mathrm{eff}$ = `Dsf.BisoEffective` ($B_\mathrm{iso}$ if isotropic, Hamilton's $B_\mathrm{eq}$ if anisotropic, 261105Cl). The floor `DiffuseScatteringFactor.BisoFloorNm2 = 1E-3` nm² (0.1 Å²) prevents a vanishing source for atoms without B (the code note says $1-e^{-2M}$ saturates to ~1 at backscattering angles, so the effect is small). **The absorption used in propagation (`getU`) does not get this floor.**
- With the same definition $\sum_n\sigma_n = (2\pi/k)\,U'_{0,\mathrm{back}}$ (when the floor does not act). In the plane-wave limit $I(t)=\mu_\mathrm{back}t$.

Inserting the Bloch expansion (`EBSDSolverManaged`):

$$\beta_{nj} = \alpha_j\sum_\mathbf g C_{\mathbf g j}e^{2\pi i\mathbf g\cdot\mathbf r_n},\qquad
S_{jj'} = \sum_n\sigma_n\,\beta_{nj}\,\beta_{nj'}^{*}\quad(S = B^{\mathsf T}\operatorname{diag}(\sigma)\,B^{*}),$$

$$I(t) = \operatorname{Re}\sum_{jj'}S_{jj'}F_{jj'}(t),\qquad I\leftarrow\max(0,I).$$

260922Cl: the old comment $S=B^\dagger\operatorname{diag}(\sigma)B$ is the complex conjugate of $S$, and the native version implemented that formula (found by the Temari audit; the managed implementation always computed the expression above). Commit 585e395 says the native side was fixed in Crystallography.Native (§5, to verify).

- Native: `NativeWrapper.EBSDSolver` (C++/Eigen, source not in this repository); otherwise `EBSDSolverManaged`.
- `EBSDSolverManaged` does not keep the full $B$; it `stackalloc`s $\beta_n$ per atom and accumulates into $S$.

#### (b) Non-local source (`IncludeTDSBackground = true`)

Instead of point sources, the absorptive potential matrix integrated over the back hemisphere is used as the source (same framework as STEM-HAADF $\int t_c^\dagger U'_\mathrm{det}t_c\,dz$, 260316Cl):

$$\frac{dI}{dz} = \frac{2\pi}{k_\mathrm{vac}}\,\boldsymbol\psi^\dagger(z)\,U'_\mathrm{back}\,\boldsymbol\psi(z),\qquad
(U'_\mathrm{back})_{\mathbf g\mathbf h} = \Bigl[U'^{\,\texttt{getU}}_{\mathbf g-\mathbf h}(\theta\in[\pi/2,\pi])\Bigr]^{*}$$

(`CreateMasterPatternMuBack`: `getU(voltage, beams[row] − beams[col], null, π/2, π, 30, 12).Imag.Conjugate()`. This is the **local** annular form (depends on $\mathbf g-\mathbf h$ only) and is different from the "non-local absorption" of §2.7; "non-local" here means that the source is not a point.)

$$M_{jj'} = \frac{2\pi}{k_\mathrm{vac}}\,\alpha_j^{*}\,(C^\dagger U'_\mathrm{back}C)_{jj'}\,\alpha_{j'},\qquad
I_\mathrm{TDS}(t) = \operatorname{Re}\sum_{jj'}M_{jj'}F_{j'j}(t).$$

The exponent of $\psi^\dagger U\psi$ is $2\pi i(\gamma_{j'}-\gamma_j^{*})=\lambda_{j'j}$, hence $F_{j'j}$ (fixed 260922Cl; the old managed code and native `_EBSDSolverWithTDS` paired it with $F_{jj'}$).

- Since 260919Cl this **replaces** the local source (in the "Disk への格納" block: `value = ebsdBackground?[r]?[t] ?? signal`; the old code added both, double-counting the same process).
- Native: `NativeWrapper.EBSDSolverWithTDS` computes $I$ and $I_\mathrm{TDS}$ with the same $F$; otherwise `EBSDSolverManaged` + `ComputeTDSMatrixBackground`.

#### (c) Smeared local source (finite `SourceSmearingSigmaNm`, trial)

`CreateSmearedLocalSourceMatrix` builds the source matrix

$$(U_\mathrm{src})_{\mathbf g\mathbf h} = \frac{1}{2\pi/k_\mathrm{vac}}\sum_n\sigma_n e^{2\pi i(\mathbf h-\mathbf g)\cdot\mathbf r_n}\;\exp\!\bigl(-2\pi^2\sigma^2\lvert\mathbf g-\mathbf h\rvert^2\bigr)$$

and evaluates it through the same matrix path $(2\pi/k)\psi^\dagger U_\mathrm{src}\psi$ as (b). Since $\lvert\mathbf g-\mathbf h\rvert$ has no 2π, the Fourier transform of an isotropic real-space Gaussian (standard deviation σ) is $\exp(-2\pi^2\sigma^2q^2)$ (same form as the DW factor $\exp(-Bs^2)$, $B=8\pi^2\langle u^2\rangle$, per the code note). $\sigma=0$ reproduces (a) (equivalence check). A finite σ takes precedence over the non-local source of `IncludeTDSBackground`.

Hypothesis (code note): most last events near the surface are small-angle scattering, so the source spreads over $\sim1/(2\pi s)$ and cannot resolve the 1s state (0.2 Å). If σ ≈ 0.5–1 Å darkens only the zone axes, that supports the hypothesis.

### 2.6 Re-injection of absorbed flux (`IncludeAbsorbedFluxBackground = true`, 260919Cl)

`ComputeAbsorbedFluxBackground`. With the cell-averaged density of the coherent wave

$$N(z) = \sum_\mathbf g\lvert\psi_\mathbf g(z)\rvert^2 = \operatorname{Re}\sum_{jj'}T_{jj'}e^{\lambda_{jj'}z},\qquad T_{jj'} = \alpha_j\alpha_{j'}^{*}\sum_\mathbf g C_{\mathbf g j}C^{*}_{\mathbf g j'}$$

($N(0)=1$), the lost part $1-N(z)$ is treated as diffuse electrons without directional information that backscatter with the average rate $\sum_n\sigma_n$:

$$D(t) = \sum_n\sigma_n\int_0^t\bigl(1-N(z)\bigr)dz = \sum_n\sigma_n\,\operatorname{Re}\sum_{jj'}T_{jj'}\bigl(t - F_{jj'}(t)\bigr),\qquad D\leftarrow\max(0,D),$$

added to whichever source is in use (the non-local source if present, otherwise the local one).

Numerical details:
- Uses $\sum T = N(0) = 1$ to evaluate $t-\sum TF$ as $\sum T(t-F)$, avoiding cancellation at weak absorption.
- For $\lvert\lambda t\rvert^2<10^{-6}$, the series $t-F = -(\lambda t^2/2)(1+\lambda t/3+(\lambda t)^2/12)$.
- Hermitian symmetry of $T$ and $t-F$: only diagonal + 2×lower triangle are summed.
- $e^{\lambda_{jj'}t} = E_{jt}E^{*}_{j't}$ with $E_{jt}=e^{2\pi i\gamma_j t}$, so only $b_\mathrm{Len}\cdot t_\mathrm{Len}$ exponentials.

Note: $N(z)$ is a sum of **densities**, not the conserved quantity without absorption, $\sum_\mathbf g P_\mathbf g\lvert\psi_\mathbf g\rvert^2$, mentioned in the 260829Cl note (§4, item 2).

### 2.7 Form of the absorptive potential (`UseNonLocalAbsorption`)

`RentMasterPatternPotentialMatrix`:

- false (default): `getPotentialMatrix(dim, b, ref A)`. $\tilde U_{\mathbf g\mathbf h}=U^{r}_{\mathbf h-\mathbf g}+iU'_{\mathbf h-\mathbf g}$ (diagonal $iU'_0$). $U'$ is the full-sphere `FactorImaginary`. Toeplitz structure: one evaluation per distinct $\mathbf g-\mathbf h$ (`uDictionary` cache).
- true: `getPotentialMatrix(dim, b, ref A, 0, π)`. The real part stays local; the imaginary part is `getU(kV, b[col], b[row], 0, π)` → `FactorImaginaryAnnular(kV, g, h, m, 0, π, 60, 20)`, depending on both $\mathbf g$ and $\mathbf h$ (2D Gauss–Legendre). The diagonal $U'(\mathbf g,\mathbf g)$ differs per $\mathbf g$. The difference is large for heavy elements with strongly forward-peaked factors (260316Cl note). Cost: one numerical integral per each of $N^2$ elements.

### 2.8 Thickness array, multiple voltages, output

- For each voltage `uDictionary.Clear()` is called before $U_0$ is read (fix 261014Cl: with the old order the second and later voltages used the previous voltage's $U_0$, because the cache key does not depend on voltage). $\sigma_n$ is also recomputed per voltage.
- Output: `Disks[vIndex][t]` is a `CBED_Disk`; `Amplitudes[r]` = $\sqrt{\max(0,\,I_r(t))}$ (the **square root of the intensity** of direction $r$, stored as a real number, no phase). $I_r(t)$ is **cumulative** over depth $0\ldots t$ (plane-wave limit $\mu_\mathrm{back}t$). Dimensionless ($\sigma$ [nm⁻¹] × depth [nm], per unit incident flux).
- Exchange with the EBSD side (definitions only):
  - Input: energies `Energies` [keV] (= `voltages` [kV]), depths `Depths` [nm] (= `thickness`, the path length above), directions (Rosca–Lambert grid, two hemispheres), grid width → `EbsdRepresentativeDirectionBlockSize`.
  - Auxiliary: `MeanAbsorptionLengthNm(double[] energiesKeV)` returns $\lambda_\mathrm{abs}(E) = k_\mathrm{vac}/(2\pi U'_0)$ [nm] (intensity attenuation length, along the path; 260923Cl), $+\infty$ if $U'_0\le0$. Used to remove double counting between MC survival and the mean absorption of the master pattern (`EbsdMonteCarloDistribution.FillPathLengthWeights`).
- Progress is reported per 50 directions of the full grid (representatives weighted by orbit size).

### 2.9 Caching caveat

`CreateMasterPatternMuBack` calls `getU(..., null, π/2, π, 30, 12)`, but the `uDictionary` key is (hkl of $\mathbf g$, `int.MaxValue`), **the same key** as the ordinary local $U$. Hence `uDictionary.Clear()` before and after the muBack stage. The potential matrices are finished before that stage, so results are correct, but any concurrent `getU` on the same instance would read back-hemisphere values (§4, item 9).

---

## 3. 1D dynamical Kikuchi band profiles (`Diffraction/Kikuchi/`, 260805Cl–)

According to the code notes, the design master is `ReciPro_菊池線動力学化設計.md` (outside this repository); § numbers below refer to it.

### 3.1 Scope and geometry

- For **transmission (forward) Kikuchi bands**. In the lab frame the beam travels along $\hat{\mathbf b}=(0,0,-1)$. The exit direction $\hat{\mathbf d}$ is near $\hat{\mathbf b}$ (by default $\angle(\hat{\mathbf d}(x{=}0),\hat{\mathbf b})\le0.35$ rad).
- Reversed geometry identical to the EBSD master pattern: $\hat{\mathbf n}=-\hat{\mathbf d}$, $\mathbf K_0=-\sqrt{k_\mathrm{vac}^2+U_0}\,\hat{\mathbf d}$.
- Momentum transfer at the source (reversed-wave beam $\mathbf h$ ↔ outgoing beam $-\mathbf h$): $\mathbf q_0 = k_\mathrm{vac}(\hat{\mathbf d}-\hat{\mathbf b})$, $\mathbf q_i = \mathbf q_0-\mathbf h_i$.
- `thickness` [nm] is the integration range $0\ldots t$ of a uniform source distribution.
- Provisional scheme "propagation = imaginary part of `getU` / source = independently built $Q$"; source and loss are not balanced (`DisplayNormalizedTag = "source-loss not balanced"`).

### 3.2 Snapshot (`KikuchiPotentialSnapshot`)

Freezes crystal, kV, kernel and absorption on/off. Holds its own `BetheMethod` instance and uses its `uDictionary` (never touches the shared `Crystal.Bethe`). $U_0 =$ `getU(kV).Real.Real`, $U'_0 =$ `getU(kV).Imag` (0 if `absorptionOff`, which also zeroes off-diagonal $U'$). Freshness `Matches`: crystal reference, kV and `ComputeCrystalHash` (cell, Z, `SubNumberElectron`, Occ, `Biso`/`Biso000`/`B11..B31`, positions, `ElasticIonModel`, `AbsorptionElasticModel`).

### 3.3 Band families (`KikuchiBandFamily`)

- `Pair`: merges $\pm\mathbf g$ into one family; the representative has its first non-zero index positive; input order preserved.
- `CollapseSystematicRows(families, xMax)` (260807Cl): families with the same direction $(h,k,l)/\gcd$ form one systematic row; only the generator (smallest $\lvert\mathbf g\rvert$) is kept, except higher orders with multiple $n>x_\mathrm{max}$, which lie outside the generator's profile range and stay as independent families. Applied automatically when `RowOrder > 0` (otherwise double counting; the code note reports MgO (Omoto Fig 6) corr 0.914 → 0.898 without collapsing and 0.927 with it).

### 3.4 Band coordinate

`BandGeometry`: $\mathbf g_\mathrm{lab}=R\,\mathbf g$, $\sin\theta_B = \lvert\mathbf g\rvert/(2k_\mathrm{vac})$, $\cos\varepsilon = \hat{\mathbf g}\cdot\hat{\mathbf b}$, $\sin\varepsilon = \sqrt{\hat g_x^2+\hat g_y^2}$ (avoids cancellation, 260806Cl). `Valid` = $\sin\varepsilon>10^{-6}$ and $\sin\theta_B<0.5$.

Profile coordinate $x = \sin\theta'/\sin\theta_B$ with $\sin\theta'\equiv-\hat{\mathbf g}\cdot\hat{\mathbf d}$ ($x\sin\theta_B$ clamped to $\pm0.99$). Of the two solutions $\hat{\mathbf d}=\cos\psi\,\hat{\mathbf b}+\sin\psi\,\hat{\mathbf g}_\perp$ the one closer to $\hat{\mathbf b}$ is taken. If that point is closer to $\hat{\mathbf b}$ than `SampleAngle` (default 0.025 rad), it is slid out of plane along the same cone $-\hat{\mathbf g}\cdot\hat{\mathbf d}=\sin\theta'$ to the point at `SampleAngle` from $\hat{\mathbf b}$ (prevents $\mathbf q_0\to0$, where the TDS source vanishes and $c$ diverges; seen on Si [111]).

### 3.5 Two-beam calculation (`ComputePoint`, `RowOrder = 0`)

Each family is computed as $\{0,+\mathbf g\}$ and $\{0,-\mathbf g\}$ (member $=\pm1$) separately with 2×2 matrices, and summed.

$$A = \begin{pmatrix} iU'_0/P_0 & U_{+\mathbf g}/P_0\\ U_{-\mathbf g}/P_1 & (iU'_0+Q_1)/P_1\end{pmatrix},\quad
P_0 = 2\hat{\mathbf n}\cdot\mathbf K_0,\ P_1 = 2\hat{\mathbf n}\cdot(\mathbf K_0+\mathbf g_m),\ Q_1 = \lvert\mathbf K_0\rvert^2-\lvert\mathbf K_0+\mathbf g_m\rvert^2$$

($U_{\pm\mathbf g}=U^r+iU'$ from `snap.GetU`; row division = 260829Cl fix; eigenvalues unchanged, eigenvectors and $\alpha$ corrected.) `DisableBraggCoupling` zeroes the off-diagonals. Points with $P_1\le10^{-3}P_0$ contribute 0.

Closed-form 2×2 eigenvalues $\gamma_{0,1} = \tfrac{a_{00}+a_{11}}{2}\pm\sqrt{(a_{00}-a_{11})^2/4+a_{01}a_{10}}$; eigenvectors $(a_{01},\gamma-a_{00})^\mathsf T$ if $\lvert a_{01}\rvert\ge\lvert a_{10}\rvert$, else $(\gamma-a_{11},a_{10})^\mathsf T$, column-normalised. If $\lvert\det C\rvert<10^{-14}$ (defective): pass-through $C=I,\ \alpha=(1,0)$. $\alpha = C^{-1}(1,0)^\mathsf T$.

### 3.6 Source density matrix $Q$ and kernels

$$Q_{ii'} = \sum_a \mathrm{Occ}_a\;\mathcal C_a(s_i^2,s_{i'}^2,s_{ii'}^2)\;e^{-2\pi i(\mathbf h_{i'}-\mathbf h_i)\cdot\mathbf r_a},\qquad s_i^2=\lvert\mathbf q_i\rvert^2/4,\ s_{ii'}^2=\lvert\mathbf h_{i'}-\mathbf h_i\rvert^2/4$$

(sum over sites expanded to symmetry-equivalent positions; phase convention as in `CreatePhaseFactors`, same structure as $\sum_n\sigma_n e^{2\pi i(\mathbf g-\mathbf g')\cdot\mathbf r_n}$ in the EBSD $S$).

`IKikuchiInelasticKernel`:
- `SourceAmplitude(a, s²)` = $\tau_a(q)$. The default `SourceCoherence` is the factorised approximation $\tau_i\tau_j$. Contract: $\mathcal C(s^2,s^2,0)=\tau^2$.
- `KikuchiTdsEinsteinKernel` (default): independent atoms, isotropic Einstein model:
  $$\tau_a^2 = f_a(s^2)^2\bigl[1-e^{-2B_a s^2}\bigr],\qquad
  \mathcal C_a = f_a(s_i^2)f_a(s_j^2)\bigl[e^{-B_a s_{ij}^2}-e^{-B_a(s_i^2+s_j^2)}\bigr]$$
  (Omoto–Tsuda–Tanaka 2002 eq. 38 = Hall–Hirsch, per the code note). $f_a$ = `ElectronScatteringPeng[Z][sub].Factor` (`ElasticIonModel` convention; the `IonFull` monopole is not included), $B_a$ = `Dsf.BisoEffective` [nm²]. If $B$ is NaN or $\le10^{-6}$, `DefaultBiso` = 0.005 nm² (0.5 Å²) is substituted and `UsedDefaultBiso = true` (response to a report that all bands vanished for B = 0).
- `KikuchiUniformSourceKernel`: $\tau=1$ (diagnostics).
- `DiagonalSourceOnly` zeroes the off-diagonals of $Q$; `SwapSourceWeights` swaps the $s^2$ of the two waves (diagnostics).

### 3.7 Intensity, reference intensity and contrast

$$S_{jj'} = \alpha_j\alpha_{j'}^{*}\sum_{ii'}C_{ij}Q_{ii'}C^{*}_{i'j'},\qquad I(t)=\operatorname{Re}\sum_{jj'}S_{jj'}F_{jj'}(t).$$

Reference intensity: the straight-through channel with Bragg coupling cut, $\gamma_\mathrm{ref}=iU'_0/P_0$, $I_\mathrm{ref}=Q_{00}F(\gamma_\mathrm{ref},\gamma_\mathrm{ref})$ (under the swap diagnostic the reference stays the real background).

$$c(x) = \frac{I(x) - n_\mathrm{ref}\,I_\mathrm{ref}(x)}{\max\bigl(I_\mathrm{ref}(x),\ 10^{-3}\max_x I_\mathrm{ref}\bigr)},\qquad n_\mathrm{ref}=\begin{cases}2&\text{two-beam (sum over }\pm\mathbf g)\\1&\text{systematic row}\end{cases}$$

(0 if $t\le0$; non-finite values set to 0). `DynamicalRank` = $\mathrm{RMS}_x(c)$.

### 3.8 Systematic row (`RowOrder = N > 0`, 260807Cl)

Beams $\{n\mathbf g: n=-N..N\}$ ($2N+1$ waves, $N$ clamped to ≤ 15) in one basis. `RowPotential` precomputes $U(m\mathbf g)$ and $e^{-2\pi i m\mathbf g\cdot\mathbf r_a}$ ($m=-2N..2N$) once per band. Matrix as in the two-beam case (diagonal $(iU'_0+Q_i)/P_i$, off-diagonal $U((c-r)\mathbf g)/P_r$). EVD via `NativeWrapper.EigenSolver` + `Inverse` if native, else MathNet. $\alpha_j=[C^{-1}]_{j,N}$. If any $P_i\le10^{-3}P_N$ (or $P_N\le0$), the point contributes 0. The $CQC^\dagger$-type contraction is two matrix products ($O(n_b^3)$).

### 3.9 `Options` defaults and selection

| Item | Default | Meaning |
|---|---|---|
| `SampleCount` | 129 | number of $x$ points (uniform over $-x_\mathrm{max}..x_\mathrm{max}$) |
| `XMax` | 2.5 | upper bound of $\lvert x\rvert$ (in Bragg half-widths) |
| `MaxScatteringAngle` | 0.35 rad | max angle between the exit direction at $x=0$ and $\hat{\mathbf b}$; beyond → `Valid=false` |
| `SampleAngle` | 0.025 rad | minimum sampling angle from the optic axis |
| `RowOrder` | 0 | 0 = two-beam path (bit-identical to before), max 15 |
| `SwapSourceWeights` / `DisableBraggCoupling` / `DiagonalSourceOnly` | false | diagnostics |

`ComputeProfiles(snap, candidates, rotation, thickness, topN, opt, previousSelection)`: computes candidates in parallel (`CollapseSystematicRows` first if `RowOrder>0`) and returns the top `topN` by rank among those with `Valid` and `DynamicalRank > 0`. Previously selected families get a 1.15× hysteresis.

---

## 4. Questions about the implementation

1. **Hexagonal grid + block sharing: the block centre can be an invalid cell** (`Diffraction/BetheMethod.cs:1600-1608, 1630`, `EBSD/EBSD.cs:235, 1191-1194`).  
   The block centre $(w_c,h_c)$ is computed in square-array coordinates without checking hexagonal validity ($\lvert u+v\rvert\le N$). The EBSD side uses $b=3$ or 5 for $G>512$, so in blocks along the rhombus boundary ($\lvert u+v\rvert\approx N$) the centre can fall on an invalid cell. That direction is `null`; `GetCrystalFixedMasterEquivalentBeamAndSurface(null, …)` returns the reference surface $(0,0,-1)$, so every direction in that block shares **the beam set selected for the [00(−1)] direction** ($P,Q$ are recomputed per direction, but the beam set is far from the Ewald sphere). Counterexample recipe: hexagonal crystal, `GridSize = 1025` ($b=3$); find a block where $(u,v)$ is valid but $(u_c,v_c)$ is not.
2. **$N(z)$ in absorbed-flux re-injection is a density sum** (`BetheMethod.cs:2108-2118, 2120-2168`).  
   After the 260829Cl row scaling, the quantity conserved without absorption is $J=\sum_\mathbf g (P_\mathbf g/P_0)\lvert\psi_\mathbf g\rvert^2$ (current along the normal), not $N=\sum_\mathbf g\lvert\psi_\mathbf g\rvert^2$. In the master pattern $\hat{\mathbf d}=-\hat{\mathbf s}$, so $1-P_\mathbf g/P_0=(Q_\mathbf g+g^2)/(2K^2)$, which is small and positive near the Bragg condition. Further investigation (2026-11-05, simplified Si, 20 kV calculation): without absorption $N\ge1$ (up to $1+6\times10^{-3}$), and `Max(0, …)` clips D to 0, so "$D\ne0$ even with $U'=0$" does not appear in practice. With absorption, D comes out 0.01–0.3% low. The summary's claim that N decreases monotonically is wrong. Details, and the question of whether $1-N$ or $1-J$ is intended, are in "A4 in detail" in [implementation-questions](implementation-questions.md).
3. **`GetCrystalFixedMasterEquivalentBeamAndSurface` reduces to the identity** (`BetheMethod.cs:815-825`).  
   Because `GetSurfaceTangentialAxes` satisfies $\hat u\times\hat v = \hat{\mathbf n}$, $R$ maps $\hat{\mathbf n}_l\mapsto\hat{\mathbf n}_r$, so $R^{\mathsf T}\hat{\mathbf n}_r = -\hat{\mathbf d}$ and $R^{\mathsf T}(-\hat{\mathbf n}_r)=\hat{\mathbf d}$. The return value is $(\hat{\mathbf d},-\hat{\mathbf d})$ regardless of the reference surface. Correct, but not what the comment ("push the inverse rotation onto beam/surface") suggests; the orientation of the tangent basis is not part of the return value.
4. **Sign label of the Kikuchi profile coordinate** (`Diffraction/Kikuchi/KikuchiProfileCalculator.cs:362`).  
   The comment says "x = +1 is the Bragg line of +g", but from $\mathbf K_0=-K\hat{\mathbf d}$ and $Q_{\mathbf g_m}=\lvert\mathbf K_0\rvert^2-\lvert\mathbf K_0+\mathbf g_m\rvert^2 = 2K\hat{\mathbf d}\cdot\mathbf g_m-g^2$, $Q=0$ for member $+\mathbf g$ is at $\hat{\mathbf d}\cdot\hat{\mathbf g}=g/2K$, i.e. $x = -k_\mathrm{vac}/K\approx-1$. It is correct if "+g" refers to the outgoing beam (reversed-wave $\mathbf h$ ↔ outgoing $-\mathbf h$). Which convention the ReciPro renderer assumes (direction of the E/D asymmetry) needs checking. Also, since $\sin\theta_B$ uses $k_\mathrm{vac}$, the Bragg lines are strictly at $\lvert x\rvert = k_\mathrm{vac}/K = 1-O(U_0/2k^2)$ (negligible).
5. **Old formula left in comments after the 260829Cl row-division fix** (`KikuchiProfileCalculator.cs:7, 592`): the header and the `ComputeRowPoint` summary write the off-diagonal as `/P_col` (the code uses `/P_row`). Comments only. Likewise `Diffraction/Alchemi/AlchemiReduction.cs:12` still writes $S=B^\dagger\operatorname{diag}(\sigma)B$ (the old notation corrected in 260922Cl; ALCHEMI document's scope).
6. **`ComputeCrystalHash` omits `Dsf.UseIso`** (`KikuchiProfileCalculator.cs:94-113`). The kernel's `BisoEffective` switches between `Biso` and `Biso000` by `UseIso`, so toggling only `UseIso` leaves a snapshot that is judged fresh.
7. **The two-beam path takes $I_\mathrm{ref}$ from member $+1$ only** (`KikuchiProfileCalculator.cs:261-267`). If the $+\mathbf g$ side returns `default` from the $P_1$ guard ($I=I_\mathrm{ref}=0$) while the $-\mathbf g$ side is valid, $I_\mathrm{ref}=0$, the floor is used and $c$ spikes. The comment says the `BandGeometry` invariants ($\sin\theta_B<0.5$ and the $\pm0.99$ clamp) guarantee $P_1>0.01P_0$, so this is latent.
8. **Atoms without B are treated differently by source and propagation** (design question). EBSD $\sigma_n$ floors $B$ at $10^{-3}$ nm², the Kikuchi kernel substitutes `DefaultBiso` = $5\times10^{-3}$ nm², while the propagation absorption `getU` has no floor ($U'=0$ for $B=0$). Atoms with B = 0 therefore "backscatter but do not absorb". Is this intended?
9. **Side effects of the shared `uDictionary`** (`BetheMethod.cs:860, 1641-1658, 3951-3963`). muBack caches annular values under the same keys as ordinary $U$ (cleared before and after). `MeanAbsorptionLengthNm` also clears the cache (two consecutive `Clear()` at the end are redundant). If these run concurrently with a calculation on the same instance, wrong $U$ values are read. Not a problem with the current call order.
10. **Cancellation in $F_{jj'}$ for small $\lambda$** (`BetheMethod.cs:2100, 2258`, `KikuchiProfileCalculator.cs:733`). The test is $\lvert\lambda\rvert^2<10^{-30}$ on $\lambda$ itself; for $10^{-15}\lesssim\lvert\lambda t\rvert\ll1$, $(e^{\lambda t}-1)/\lambda$ has relative error $\sim\epsilon/\lvert\lambda t\rvert$. Only `ComputeAbsorbedFluxBackground` uses a series. No practical effect for normal absorption ($\lvert\lambda t\rvert\gtrsim10^{-3}$).

---

## 5. Items to verify

1. Whether native `Crystallography.Native` `_EBSDSolver` now computes $S=B^{\mathsf T}\operatorname{diag}(\sigma)B^{*}$ (as managed) and `_EBSDSolverWithTDS` uses $F_{j'j}$, and whether the bundled DLL is that version. The source is not in this repository (commit 585e395 says "fixed in Crystallography.Native").
2. Physical reference axis of the back-hemisphere range $\theta\in[\pi/2,\pi]$ (reversed-wave direction or original incident beam). To be checked together with the angle definition of `FactorImaginaryAnnular` in [scattering-factors.md](scattering-factors.md).
3. The local point source uses $\lvert\psi(\mathbf r_n)\rvert^2$ at the mean atomic position; thermal spread of the source enters only through $\sigma_n$ (cross-section containing $1-e^{-2M}$). Is this the intended approximation (`SourceSmearingSigmaNm` being the trial alternative)?
4. The non-local source of `IncludeTDSBackground` uses the local annular $U'_\mathrm{back}(\mathbf g-\mathbf h)$. Is a non-local $U'_\mathrm{back}(\mathbf g,\mathbf h)$ planned, as for `UseNonLocalAbsorption`?
5. Whether the accuracy of `EbsdRepresentativeDirectionBlockSize > 1` (degraded beam selection at block edges) has been evaluated.
6. Bibliographic details of the design document `ReciPro_菊池線動力学化設計.md` and Omoto–Tsuda–Tanaka 2002 (the code gives only names and equation numbers).
7. Whether falling back to EVD in the EBSD path for `Solver.MtxExp_*` is intended.
8. Range of validity of solving $\pm\mathbf g$ as separate 2×2 problems and summing (ignoring simultaneous excitation of $+\mathbf g$ and $-\mathbf g$) in the two-beam Kikuchi path; is `RowOrder ≥ 1` the intended improvement?
9. Whether the denominator floor $10^{-3}\max I_\mathrm{ref}$ in `ComputeProfile` (implementation of ε in design §4) matches the design.
10. Whether the EBSD side (separate document) squares `Amplitudes` back to intensity and builds depth weights from differences of the cumulative values (exchange contract).
11. Whether it is intended that `beamDirections` must be two concatenated hemispheres (single-hemisphere detector-pattern calculation is currently impossible because the old path `ebsd_DoWork` was retired).

Items to verify: 11.
