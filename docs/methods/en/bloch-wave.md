Draft written by Claude from the code (2026-11-05). Not yet reviewed by the author.

# Bloch-wave (Bethe) method — potential, eigenproblem, CBED, diffracted amplitudes, PED

Scope: the `BetheMethod` class in `Diffraction/BetheMethod.cs`, excluding EBSD ([kikuchi-dynamical.md](kikuchi-dynamical.md)) and ALCHEMI ([alchemi-ionization.md](alchemi-ionization.md)).
HRTEM, STEM and potential images are in [bloch-wave-imaging.md](bloch-wave-imaging.md).
The scattering factors themselves (Peng coefficients in `AtomStatic.ES`, absorptive form factors `FactorImaginary*`) are in [scattering-factors.md](scattering-factors.md); this page only describes how `getU` uses them.

Line numbers refer to `BetheMethod.cs` as of 2026-11-05.

---

## 1. Conventions

| Item | Implementation | Source |
|---|---|---|
| Length unit | **nm** (reciprocal space nm⁻¹) | class summary "すべて、単位はnm" (all units nm) |
| Wave number | $k = 1/\lambda$ (no 2π) [nm⁻¹] | `UniversalConstants.Convert.EnergyToElectronWaveNumber` ("not divided by 2π") |
| Accelerating voltage | kV (`AccVoltage`) | |
| Plane wave | $\exp(+2\pi i\,\mathbf{k}\cdot\mathbf{r})$ | propagator `Exp(TwoPiI * γ * t)`, image synthesis `Exp(TwoPiI * g·r)` |
| Scattering angles (apertures, detectors, PED semi-angle) | **radians** | `Math.Sin(angle)/λ`, `Matrix3D.Rot(axis, θ)` |
| Frame | +X right, +Y up, +Z toward the viewer; incidence along $-Z$ by default | `Surface` summary, `getVecK0` |
| Surface normal `Surface` | unit vector pointing **into** the specimen, default $(0,0,-1)$ | `Surface` summary |
| Thickness $t$ | depth along the surface normal (not the path length along the beam) | 260829Cl comment in `cbed_DoWork` (~L503): multiplying by 1/cosτ would double-correct, so it is not applied |

### 1.1 Wavelength and relativistic correction

`UniversalConstants.Convert`:

$$
\lambda\,[\mathrm{nm}] = \frac{1.2264262862108}{\sqrt{V\,(1 + 0.97847537252\times10^{-6}\,V)}},\qquad V = 1000\times\text{kV}
$$

(the numerical form of the commented formula $\lambda = h/\sqrt{2m_0 eV(1+eV/2m_0c^2)}$). `EnergyToElectronWaveNumber` is its reciprocal; `WaveLengthToElectronEnergy` is the positive root of $aE^2 + bE - c^2/\lambda^2 = 0$ (E in keV).

The relativistic factor inside `getU` is

$$
\gamma = 1 + \frac{e_0 \cdot 10^3\,\text{kV}}{m_0 c^2}
$$

(`UniversalConstants.e0, m0, c2`). The constants are `m0 = 9.1093897E-31`, `e0 = 1.60217733E-19` (CODATA-1986-era values) and `h = 6.62606896E-34` (2006). The wavelength coefficients 1.22642629 and 0.978475e-6 are separate literals, not recomputed from the same constant set as `gamma`.

`WaveProperty` (`Diffraction/WaveProperty.cs`) is a plain container: radiation type, wavelength [nm], energy [eV], characteristic-line element and line. Conversions call `Convert` above. For electrons `Energy` is in eV and is divided by 1000 (→ kV) for the wavelength. 260718Cl fixed the Å→nm conversion and the wavelength-assigned-to-Energy bug in the characteristic-line constructor.

`DiffractionOptics` (`Diffraction/DiffractionOptics.cs`) only reads the detector-parameter XML used by IPAnalyzer/ReciPro (string fields plus renaming of legacy tags); no physics.

`DynamicTheory` (`Diffraction/DynamicTheory.cs`) is an abandoned multislice class whose implementation is entirely commented out; it is `[Obsolete]` (260718Cl, unused across repositories, kept for macro compatibility). All dynamical calculations in ReciPro go through `BetheMethod`.

---

## 2. Fourier coefficients of the potential, `getU`

`getU(in double kV, in Beam g, in Beam h = null, double inner = NaN, double outer = NaN, int nTheta = 60, int nPhi = 20)` (L3863–3934)

The return value `(Complex Real, Complex Imag)` is

- `Real` $= \hat U_\mathbf{g}$: Fourier coefficient of the elastic (real) potential $U(\mathbf r)$. **Complex** (carries a phase for non-centrosymmetric crystals).
- `Imag` $= U'_\mathbf{g}$: Fourier coefficient of the absorptive (imaginary) potential $U'(\mathbf r)$; also complex.

The total crystal potential is $U(\mathbf r) + iU'(\mathbf r)$; matrix elements are built as `Real + i·Imag` (§3).

### 2.1 Implemented formula

$s^2 = |\mathbf g|^2/4$ ($s = \sin\theta/\lambda$, nm⁻²); in the non-local form with $h$ given, $s^2 = |\mathbf g - \mathbf h|^2/4$ and the index is $\mathbf g - \mathbf h$.

$$
\texttt{Real} = \frac{\gamma}{\pi V_\mathrm{cell}} \sum_{\text{site } a}\; \sum_{j\in a} \mathrm{Occ}_a\, f^{e}_a(s)\, e^{-M_{aj} s^2}\, e^{+2\pi i\,\mathbf g\cdot \mathbf r_j}
$$

$$
\texttt{Imag} = \frac{\gamma}{\pi V_\mathrm{cell}} \sum_{a}\sum_{j\in a} \mathrm{Occ}_a\, f'_a(s; M_{aj})\, e^{-M_{aj} s^2}\, e^{+2\pi i\,\mathbf g\cdot \mathbf r_j}
$$

- $V_\mathrm{cell}$ = `Crystal.Volume` [nm³], $f^e$ [nm] → $U$ in **nm⁻²**. Code comment: "for the coefficients see Kirkland's textbook p.120".
- $U$ corresponds to $U = 2meV/h^2$ ($m=\gamma m_0$) when the Schrödinger equation is written $[\nabla^2 + 4\pi^2(k^2 + U)]\psi = 0$ (consistent with `getVecK0` using $|\mathbf K|^2 = k^2 + U_0$, §4).
- **The phase sign is $+2\pi i\,\mathbf g\cdot\mathbf r$** (20240524 comment: sign changed to plus "so that centrosymmetric crystals should now work — Mitsubishi / Nakamura"). This is opposite to the usual crystallographic $V_\mathbf g = V_\mathrm{cell}^{-1}\int V e^{-2\pi i\mathbf g\cdot\mathbf r}$, so **`getU(g)` equals the textbook $U_{-\mathbf g}$**. Matrix elements are built with `getU(b[col] - b[row])` (§3), so element $(g,h)$ is the textbook $U_{\mathbf g-\mathbf h}$. This "reversed-sign coefficient called with the reversed difference" structure is the background of the `.Imag.Conjugate()` in STEM ([bloch-wave-imaging.md](bloch-wave-imaging.md) §2.5) and of the direct negative-phase construction in STEM-EDX (260801Cl, `ComputeIonizationU`).
- Sum over `Crystal.Atoms` (independent sites) × `atoms.Atom` (all equivalent positions), times occupancy `Occ`.

### 2.2 Elastic factor and ions

`AtomStatic.ElectronScatteringPeng[Z][sub]` (Peng, 5 Gaussians). `es.Factor(s2)` **takes s² in nm⁻² and returns f in nm** (×0.01, ×0.1 internally). The 8-Gaussian table (`ElectronScatteringEightGaussian`) is commented out.

Static setting `BetheMethod.ElasticIonModel` (default `Neutral`, 260613Cl):

| Model | `sub` | Monopole |
|---|---|---|
| `Neutral` (default) | 0 (neutral atom; assigned ion ignored) | none |
| `IonScreened` | `atoms.SubNumberElectron` | none (finite part of Peng ion fit only; former default) |
| `IonFull` | `atoms.SubNumberElectron` | adds the Mott–Bethe monopole $0.2393366\,\Delta Z/s^2$ [nm] ($s^2$ in nm⁻²) for $g\neq0$ |

`IonFull` adds the monopole only when the cell is neutral (`|NetCellCharge()| ≤ 1e-6`, $\sum \text{Multiplicity}\cdot\text{Occ}\cdot\text{Valence}$); otherwise it silently falls back to `IonScreened` behaviour (260614Cl). No monopole at $s^2 = 0$ ($U_0$, mean inner potential).

### 2.3 Debye–Waller factor

`atoms.Dsf` (`DiffuseScatteringFactor`; details in [scattering-factors.md](scattering-factors.md)).

- Isotropic: $M = B_\mathrm{iso}$ [nm²], $e^{-B s^2}$ with $s^2 = g^2/4$. Conversion from $U$ ($B = 8\pi^2 U$) is in `Dsf.Biso`.
- Anisotropic, $g\neq0$: for each equivalent position $j$ the index is rotated by `atom.Operation.ConvertPlaneIndex(index)` and

  $$M_{aj} = \frac{\beta_{11}H^2 + \beta_{22}K^2 + \beta_{33}L^2 + 2\beta_{12}HK + 2\beta_{23}KL + 2\beta_{31}LH}{s^2}$$

  (L3912; `dsf.B11…` read as the dimensionless β form), i.e. converted to an "effective isotropic B along g" so that $e^{-Ms^2} = e^{-\beta_{ij}h_ih_j}$.
- Anisotropic, $g = 0$: `dsf.BisoEffective` (261105Cl; Hamilton's $B_\mathrm{eq}$ = `Biso000` for anisotropic).
- NaN `M` → 0. `dsf.IsZero` → no absorption.

### 2.4 Absorption (imaginary part)

`AbsorptionSource ImaginaryPotentialAbsorption => AbsorptionSource.TdsRedistributable` (L3848): **the only source of the imaginary part is currently TDS** (no mean absorption, no empirical damping; the note says to update this declaration if that changes).

$f'$ comes from these `AtomStatic.ES` methods; it is returned in nm like $f^e$ (the code comment says "dimensionless", but `FactorImaginary` returns $\gamma k_0/2 \times \int f f\,(\ldots)\,d\Omega \times 0.01$, i.e. f[Å]² × k₀[nm⁻¹] × 0.01 = nm):

| Call | Method | Meaning |
|---|---|---|
| `inner`/`outer` NaN (default) | `es.FactorImaginary(kV, s2, M)` | TDS absorption over the full solid angle (local form $U'(\mathbf g)$) |
| `inner`,`outer` given, `h = null` | `es.FactorImaginaryAnnular(kV, g.Vec, M, inner, outer, nTheta, nPhi)` | local form restricted to scattering angles [inner, outer] rad (STEM ADF) |
| `inner`,`outer` given, `h` given | `es.FactorImaginaryAnnular(kV, g.Vec, h.Vec, M, …)` | non-local $U'(\mathbf g,\mathbf h)$ (EBSD `UseNonLocalAbsorption`) |

The integrand is $f(|\mathbf k-\mathbf g/2|)f(|\mathbf k+\mathbf g/2|)\,[1-e^{M(s^2-s_1^2-s_2^2)}]$; `getU` multiplies by $e^{-Ms^2}$, giving the Hall–Hirsch-type $e^{-Ms^2} - e^{-M(s_1^2+s_2^2)}$. The $f$ inside the integrand is `ElasticFactorWithTail` (static `AbsorptionElasticModel`, default `TemariHybrid`, 260922Cl), using Peng → Mott–Bethe tail or Temari $f_e$ (this does not affect the real part `es.Factor`).

For isotropic sites $M$ is the same for all equivalent positions, so $f'$ is computed once per site. For anisotropic sites with $g\neq0$ the numerical integral is redone for every equivalent position (condition at L3902).

### 2.5 Cache `uDictionary`

`ConcurrentDictionary<(int Key1, int Key2), (Complex, Complex)>`. Key `(compose(g.Index), h is null ? int.MaxValue : compose(h.Index))`; it **does not include kV, inner/outer, nTheta/nPhi, ElasticIonModel or AbsorptionElasticModel**. `compose` packs each index +255 into 10 bits (range roughly ±255). Stored only when `kV > 0`.

Code paths that change voltage, detector angles or models are therefore expected to call `uDictionary.Clear()` / `ClearUCache()`: `cbed_DoWork` (L452, 462), `GetPrecessionElectronDiffraction` (L2360) and `StemDoWork` (L2523, 2535, 2877) clear unconditionally (261014Cl); `GetDifractedBeamAmpriltudes` clears only when the voltage changes (L2282); `MeanAbsorptionLengthNm` clears for each energy and once more at the end. → Question Q1.

Allocation-saving overloads (`getU(kV, index, vec)`, the STEM `getU(kV, index, vq, vi, vj, inner, outer, …)`) just return the cached value without creating a `Beam`; values are identical to the main method.

### 2.6 Mean absorption length `MeanAbsorptionLengthNm`

$$
\lambda_\mathrm{abs}(E) = \frac{k_\mathrm{vac}}{2\pi\,U'_0}\quad[\text{nm}]
$$

(intensity attenuation length; $+\infty$ if $U'_0 \le 0$; 260923Cl; used for EBSD depth distributions). Consistent with a plane wave $e^{2\pi i(K+i\kappa)z}$, $(K+i\kappa)^2 = K^2 + iU'$, $\kappa \approx U'/2K$, intensity $e^{-2\pi U' z/K}$. Vacuum $k$ is used instead of $K$ (difference $U_0/2K^2 \sim 10^{-4}$ per the comment).

---

## 3. Matrices

### 3.1 Potential matrix `getPotentialMatrix` (L4084–4136)

Column-major 1-D array `potentialMatrix[row + col*dim]`:

$$
\mathcal U_{gh} = \begin{cases} i\,U'_0 & g = h \\ \hat U_{\mathbf g-\mathbf h} + i\,U'_{\mathbf g-\mathbf h} & g\neq h\end{cases}
\qquad(\texttt{getU(b[col] - b[row])},\ g=\text{row},\ h=\text{col})
$$

- The diagonal carries only the absorption $U'_0$; $\hat U_0$ (mean inner potential) is moved into the refraction of $\mathbf K$ (§4).
- The non-local version (`getPotentialMatrix(b, nonLocalInner, nonLocalOuter)`) uses `getU(AccVoltage, b[col], b[row], inner, outer)`: same real part, imaginary part $U'(\mathbf g,\mathbf h)$.

### 3.2 Eigenproblem `getEigenMatrix` (L4147–4210)

With

$$
\psi(\mathbf r) = \sum_j \alpha_j \sum_\mathbf g C_\mathbf g^{(j)} \exp\!\big[2\pi i(\mathbf K + \mathbf g + \gamma^{(j)}\mathbf n)\cdot\mathbf r\big]
$$

and $|\mathbf K|^2 = k^2 + U_0$, dropping $\gamma^2$ (forward-scattering approximation), the linearised Bethe equation

$$
P_\mathbf g\,\gamma\,C_\mathbf g = Q_\mathbf g\,C_\mathbf g + \sum_\mathbf h \mathcal U_{\mathbf g\mathbf h}\,C_\mathbf h
$$

is solved with

$$
A = D_P^{-1}\big(\mathrm{diag}(Q) + \mathcal U\big),\qquad A_{gh} = \frac{\mathcal U_{gh} + Q_g\delta_{gh}}{P_g}
$$

(**row g divided by $P_g$**). 260829Cl restored this row division: the 2021 conversion to 1-D arrays had regressed to column division $\mathcal U D_P^{-1}$. Eigenvalues are invariant under the similarity transform but eigenvectors became $P_gC_g$, multiplying diffracted amplitudes by $P_g/P_0$ (the non-absorbing conserved quantity became $\sum|\psi_g|^2/P_g$ instead of $\sum P_g|\psi_g|^2$). The comment cites the check tool `tools/BetheDivCheck` and the public manual appendix a3-bloch-wave.

$Q, P$ (`Beam.Q`, `Beam.P`, L4639–4645):

$$
Q_\mathbf g = |\mathbf K|^2 - |\mathbf K + \mathbf g|^2,\qquad P_\mathbf g = 2\,\mathbf n\cdot(\mathbf K + \mathbf g)
$$

(`Beam.Q` summary: "$-g(2k_0+g)$; slightly different from the $2k_0 S$ excitation error"; `P` summary: "after Mr. Otsuka's notes".) The excitation error measured along the surface normal is

$$
S_\mathbf g = \sqrt{P_\mathbf g^2/4 + Q_\mathbf g} - P_\mathbf g/2
$$

(`Beam.S`), the root of $|\mathbf K + \mathbf g + S\mathbf n| = |\mathbf K|$; positive inside the Ewald sphere ($Q>0$).

The subset version `getEigenMatrix(bLen, beams, gMap, baseLen, basePotential, ref eigenMatrix)` (260602Cl) extracts rows/columns of a precomputed potential matrix through `gMap` and builds the same $A$ (for EBSD; fixes non-local absorption silently reverting to the local form).

### 3.3 Solution and boundary conditions

$\psi_\mathbf g(0) = \delta_{\mathbf g\mathbf 0}$ at the entrance surface, $\boldsymbol\alpha = C^{-1}\boldsymbol\psi_0$ (`LU().Solve` or explicit inverse). At depth $t$ (along the normal)

$$
\psi_\mathbf g(t) = e^{\pi i P_\mathbf g t}\sum_j C_\mathbf g^{(j)}\,e^{2\pi i\gamma^{(j)} t}\,\alpha_j
$$

The prefactor $e^{\pi iP_\mathbf g t} = e^{2\pi i\,\mathbf n\cdot(\mathbf K+\mathbf g)t}$ is the "phase that accounts for the exit-surface boundary condition" (GetDifractedBeamAmpriltudes: 20230827; CBED: 20220803). Differences by path:

| Path | Exit phase |
|---|---|
| `GetDifractedBeamAmpriltudes` | $e^{\pi i P_g t}$ |
| `cbed_DoWork` | $e^{\pi i (P_g - 2k_\mathrm{vac} n_z) t}$ (extra phase $e^{-2\pi i k n_z t}$ common to all beams; no effect on intensity) |
| `StemDoWork` | $e^{2\pi i (P_g/2) t}$ (same). Incident vector $\psi_0 = -1$ (global sign) |
| `GetPrecessionElectronDiffraction` | none (only intensities are used) |

With absorption $A$ is non-Hermitian, so MathNet `Evd(Symmetricity.Asymmetric)` or native Eigen is used. No Bloch-wave normalisation of $C$ is needed; $\alpha = C^{-1}\psi_0$ absorbs it.

---

## 4. Incident wave vector `getVecK0` (L4658–4679)

The crystal wave vector keeps the surface-parallel component and changes only the normal component:

$$
\mathbf K = \mathbf k_\mathrm{vac} + x\,\mathbf n,\quad |\mathbf K|^2 = k_\mathrm{vac}^2 + \hat U_0
\ \Rightarrow\ x = \sqrt{b^2 + \hat U_0} - b,\ \ b = \mathbf n\cdot\mathbf k_\mathrm{vac}
$$

$\mathbf k_\mathrm{vac} = (0,0,-k)$ if `beamDirection` is null, otherwise $k\,\hat{\mathbf b}$. $\hat U_0$ = `getU(voltage).Real.Real` (mean inner potential, no absorption). With `UseLocalSurfacePerBeamDirection = true` (EBSD master pattern) $\mathbf n = -\hat{\mathbf b}$ (`GetEbsdSurfaceNormal`).

---

## 5. Beam selection `Find_gVectors` (L4352–4552)

Two phases (260318Cl).

### 5.1 Phase 1: BFS cache `Build_gCache`

- BFS from the origin with a lattice-type-specific neighbour set (F: 8 ⟨111⟩, I: 12 ⟨110⟩, A/B/C, R (hexagonal setting), hexagonal/trigonal, P: 6 ⟨100⟩). **261030Cl** changed the sets from `FrozenSet` to arrays (the 12-element I set enumerated in a per-process hash order and results did not reproduce across processes).
- A priority queue keyed by $|\mathbf g|$ is expanded one shell at a time (everything below min $|\mathbf g|$ + 3·shift; shift = half of the longest neighbour vector) until `cacheLimit = max(64N, 50,000)` points.
- The result `(key, gX, gY, gZ, |g|)` (lab frame, nm⁻¹) and the covered radius `CoverRadius = (front min |g|) − 6·shift` form an immutable `GCacheSnapshot`, swapped once under a lock (261004Cl: PED parallel calls could see a new orientation with an old cache). Reused if the orientation matches and it is long enough.

### 5.2 Phase 2: Ewald-shell test and ranking

With $\mathbf v = \mathbf K + \mathbf g$, $q = |\mathbf K|^2 - |\mathbf v|^2$:

- Accept if $|q| < \text{maxQ} = \big||\mathbf K|^2 - (|\mathbf K| + \text{shift})^2\big|$ and $\mathbf n\cdot\mathbf v > 0$ (forward-going).
- Collect up to `limit = 8N` candidates from the head of the cache (≈ increasing $|\mathbf g|$) (default; with the diagnostic `BeamCandidateMultiplier` ≠ 0, exactly $N\times$multiplier, growing the cache if needed).
- **Rating = $|\mathbf g|\,Q^2$** (`Beam.Rating`): QuickSelect the $N+1$ smallest and sort ascending. $\mathbf 0$ has Rating 0 and comes first (all solvers treat index 0 as the incident beam).
- **Duplicate removal in surface-tangential coordinates**: beams whose projections on the tangent-plane basis (`GetSurfaceTangentialAxes`) round to the same key at $10^{-6}$ nm⁻¹ are reduced to the one with smaller Rating (260321Ch).
- **Drop the trailing tie group** (`BeamBoundaryTieMode.DropGroup`, default): beam $N+1$ is a sentinel. Find the last index $n$ where the Rating step exceeds $10^{-6}$ and return `beams[..n]`. The result therefore has **at most N beams**, and a symmetry-equivalent group straddling the boundary is dropped as a whole. `IncludeGroup` (diagnostic) includes the whole group.
- `getU(AccVoltage, hkl, g)` is stored on each beam (`Beam.Ureal/Uimag`).

Note (from comments): with the default limit 8N, for N ≳ 250 candidates can run out at the cache edge first (Ewald-shell candidates grow only as the 2/3 power of the cache size), so the candidate range is set by the cache radius (|g| ≈ 80 nm⁻¹ for MnAl). Exhaustion count: `BeamCandidateExhaustedCount`; last scan: `LastBeamCandidateScan` (ThreadStatic).

### 5.3 `reset_gVectors`

Recomputes only $Q, P$ for a new $\mathbf K$ (incident direction) and surface normal, keeping hkl, g and U. Used in the CBED and STEM direction loops; the beam set is chosen once at the central direction.

---

## 6. Solvers

`enum Solver { Eigen_MKL, Eigen_Eigen, MtxExp_MKL, MtxExp_Eigen, Auto }`

| Name | Content |
|---|---|
| `Eigen_Eigen` | native (`Crystallography.Native.dll`, C++ Eigen) EVD → $\psi(t)$; `NativeWrapper.CBEDSolver_Eigen` / `CBEDSolver2` (STEM; also returns eigenvalues, eigenvectors, α) |
| `Eigen_MKL` | MathNet `Evd(Asymmetric)`; MKL provider if loaded, otherwise managed (shown as `Eigen_Managed`) |
| `MtxExp_Eigen` | native $\exp(2\pi i\,\Delta t\,A)$ applied repeatedly (uniform thickness steps only) |
| `MtxExp_MKL` | MathNet `Matrix.Exponential()`, rebuilt whenever the step width changes (260711Cl; previously wrong for start≠step or non-uniform lists) |

- CBED `Auto`: `MtxExp_Eigen` with all logical cores if native Eigen is present, otherwise `Eigen_MKL` (cores/4 with MKL, all cores without).
- MtxExp requires uniform thicknesses (`NativeWrapper.IsUniformThickness`, relative 1e-12); otherwise EVD is used instead (260711Cl; also guarded natively).
- STEM needs eigenvalues/eigenvectors/α downstream, so MtxExp is normalised to EVD. `Auto` → `Eigen_Eigen` with native, except `Eigen_MKL` when MKL is active and bLen ≥ `MklStemBlenThreshold` (default 500; env `RECIPRO_STEM_MKL_BLEN`).
- `GetDifractedBeamAmpriltudes` and PED use native `EigenSolver` + `Inverse` when native is present and N < 400, otherwise MathNet EVD + `Inverse()`.
- The native library source is not in this repository (`NativeWrapper.cs` contains only the P/Invoke declarations).

---

## 7. CBED / LACBED `RunCBED` → `cbed_DoWork` (L430–735)

Inputs: `maxNumOfBloch` (N), voltage [kV], base orientation, thicknesses [nm], incident directions `beamDirections` (unit vectors on a square grid), `LACBED`, solver, threads.

1. `uDictionary.Clear()`, $\hat U_0$, $\mathbf K$ for the central direction, and one `Find_gVectors` call fix the beam set; the potential matrix is built once.
2. `beamDirections` is treated as a $w\times w$ grid, $w = \sqrt{\text{Length}}$; only directions inside the inscribed circle (tested at pixel centres) are computed.
3. Per direction: `getVecK0(kvac, u0, dir)` → `reset_gVectors` → `getEigenMatrix` → solver → exit phase (§3.3).
4. LACBED: only the transmitted beam (index 0) amplitudes go into `Disks[t][0]`.
5. CBED: amplitudes are arranged per disk with detector coordinates $\mathbf p = ((k\hat{\mathbf b}+\mathbf g)_x/(\ldots)_z,\ (\ldots)_y/(\ldots)_z)$ (camera length 1). Where disks overlap, **intensities $|\psi|^2$ of the other disks are added incoherently** and `Amplitudes = √(ΣI)` (real). Overlaps are found from a bounding-rectangle intersection graph (260610Cl). The nearest-pixel search `getIndex` marks disk-edge pixels as −1 (invalid).

---

## 8. Diffracted amplitudes `GetDifractedBeamAmpriltudes` (L2280–2341)

Single direction (`rotation`, incidence $-Z$), single thickness. If N, voltage and orientation are unchanged, the eigendecomposition is reused and only the thickness changes (used by the HRTEM thickness loop). `uDictionary` is cleared only when the voltage changes. Returns `Beam.Psi` = $\psi_g$ including the exit phase $e^{\pi i P_g t}$.

---

## 9. Precession electron diffraction `GetPrecessionElectronDiffraction` (L2355–2462)

- Semi-angle `semiangle` [rad], `step` divisions. Direction $k$ is $R(\hat{\mathbf u}_k, \alpha)\,(0,0,-1)$ with $\hat{\mathbf u}_k = (\cos 2\pi k/\text{step}, \sin 2\pi k/\text{step}, 0)$.
- Per direction `Find_gVectors` (sharing the BFS cache) → eigenproblem. Reused if N, voltage, orientation, semi-angle and step are unchanged.
- **Incoherent average** over directions: $I_\mathbf g = \frac1{\text{step}}\sum_k |\psi_\mathbf g^{(k)}(t)|^2$ (a g absent from a direction's beam set contributes 0 there).
- Finally $Q, P$ are recomputed at the base orientation, `Psi = √I` (real), sorted by Rating.
- No exit phase (not needed for intensities).

---

## 10. Main inputs and defaults (this page)

| Parameter | Default | Where |
|---|---|---|
| `Surface` | (0,0,−1) | property |
| `ElasticIonModel` | `Neutral` | static |
| `AbsorptionElasticModel` | `TemariHybrid` | static |
| annular quadrature `nTheta`, `nPhi` in `getU` | 60, 20 | argument defaults |
| θ segmentation / points in `FactorImaginary` | geometric in s: [0, 0.3, 1, 3, 9, 27, 81, 243, 729, ∞] nm⁻¹ / θ 12 points, φ 8 points (folded) | `AtomStatic.ES` ([scattering-factors.md](scattering-factors.md)) |
| BFS cache size | max(64N, 50,000) | `Find_gVectors` |
| Candidate limit | 8N (`BeamCandidateMultiplier = 0`) | same |
| Tie and tangential-duplicate tolerance | 1e-6 | same |
| `BeamBoundaryTieMode` | `DropGroup` | static |
| `MklStemBlenThreshold` | 500 | static (env var) |
| `EbsdRepresentativeDirectionBlockSize` | 3 | EBSD |
| `UseLocalSurfacePerBeamDirection` | false | EBSD master |

---

## Implementation questions

**Q1. The U cache key ignores detector angles; after STEM the elastic path may read the annular U′**
The `getU` key is `(compose(g), int.MaxValue)` and contains neither kV nor inner/outer (L3865–3867, 3976–3978, 3997–3999). `StemDoWork` clears at L2877 and then `fillUq` (L2895) stores **U′ restricted to the annulus [detAngleInner, detAngleOuter]** under the same key as the ordinary local U, and does not clear at the end. A subsequent `GetDifractedBeamAmpriltudes` on the same instance at the same voltage (L2282 clears only on a voltage change) lets `getU` inside `Find_gVectors` (L4466) hit the cache, so the absorption becomes the ADF-annulus TDS instead of the full-solid-angle TDS. `GetHRTEMImage` calls this function internally. CBED, PED, STEM and EBSD clear unconditionally at the start and are unaffected. Whether ReciPro shares one `BetheMethod` between STEM and HRTEM/diffraction cannot be verified from this repository.

**Q2. The temporal-coherence envelope differs by a factor 2 between STEM and HRTEM**
HRTEM (L3735 quasi-coherent, L3763 TCC) uses $\exp[-\tfrac12\pi^2\lambda^2\delta^2(q_1^2-q_2^2)^2]$; the STEM `Lenz` (L2741–2742) uses $\exp[-\pi^2\lambda^2\delta^2(k_q^2-k^2)^2]$ without the $\tfrac12$. If $\delta$ is the rms focal spread (standard deviation of a Gaussian), averaging gives the former. The same input "delta = Cc·ΔE/E" therefore acts as $\sqrt2\,\delta$ in STEM. One definition (rms or 1/e width) should be fixed. See [bloch-wave-imaging.md](bloch-wave-imaging.md).

**Q3. Swapped comments in the HRTEM quasi-coherent branch** (numbers are correct)
The β term at L3734 is labelled "Ec temporal incoherence" and the δ term at L3735 "Es spatial incoherence"; it is the other way round. The TCC branch (L3762–3763) is labelled correctly.

**Q4. Unit in the doc comment of `MottBetheMonopoleCoefficient`** (numbers are correct)
`AtomStatic.cs` L1210 says "s² = (sinθ/λ)² [Å⁻²], f in [nm]". Converting Mott–Bethe $f[\text{Å}] = 0.023934\,\Delta Z/s^2[\text{Å}^{-2}]$ to f[nm], s²[nm⁻²] gives 0.1×0.023934/0.01 = 0.23934, so 0.2393366 is the coefficient **for s² in nm⁻²**. `getU` (L3893) passes s2 in nm⁻², so the numbers are right, and the comment near `ElasticFactorWithTail` ("s² in nm⁻², f in nm") agrees. Only L1210 disagrees.

**Q5. The eigendecomposition reuse test ignores crystal and model changes**
`GetDifractedBeamAmpriltudes` (L2292) and PED (L2370) decide recomputation from N, voltage and orientation only (PED also semi-angle and step), not from atoms, DW factors, `Surface`, `ElasticIonModel` or `AbsorptionElasticModel`. PED clears the U cache at the start, but if the eigendecomposition is reused the new U is never used. This is harmless if callers recreate `BetheMethod` whenever the crystal is edited (to be confirmed).

---

## Items to confirm

1. Duplicate removal by tangential coordinates in `Find_gVectors` (L4475–4493) reduces beams that differ only along the surface normal (e.g. HOLZ rows where the difference of g is parallel to n) to one. Is this intended as the exit condition of the linearised Bethe method?
2. With anisotropic DW, the absorption integrand $1-e^{M(s^2-s_1^2-s_2^2)}$ uses the "effective isotropic B along g" (L3912 → `FactorImaginary(kV, s2, m)`). The intermediate TDS vectors $s_1, s_2$ are not along g, so this is an approximation. Acceptable?
3. The `getU` comment "inelastic factor, the answer is dimensionless" (L3918) vs. `FactorImaginary` returning nm (added with the same unit as the real part). Should the summaries be corrected?
4. Intent of the $-2k n_z$ term in the CBED exit phase $e^{\pi i(P - 2k n_z)t}$ (L589) (common to all beams, no effect on intensity). Does anything downstream use the phase of `Disks[].RawAmplitudes`?
5. If `UseLocalSurfacePerBeamDirection = true` remains set when CBED/STEM run, `getVecK0` (per-direction normal) and `reset_gVectors` (fixed `Surface`) use different normals. Can the EBSD flag remain on a shared instance?
6. Constants of mixed vintage (m0, e0 from 1986; h from 2006; λ coefficients separate literals). Intentional?
7. A beam with exactly $Q = 0$ (Rating 0) in `Find_gVectors` would tie with $\mathbf 0$ and an unstable sort could swap index 0 (exact zero is practically impossible in floating point).
8. `getIndex` (CBED disk overlap) invalidates disk-edge pixels, so intensity is not added on the 1-pixel rim of overlap regions. Intended?

(8 items to confirm, 5 implementation questions)
