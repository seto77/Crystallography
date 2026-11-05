This is a draft written by Claude from the code (2026-11-05). Not yet reviewed by the author.

# Inner-shell ionization (STEM-EDX) and ALCHEMI — computational methods

Scope: `Diffraction/IonizationChannel.cs`, `Diffraction/XrayLineSeries.cs`, `Diffraction/AbsorptionSource.cs`,
`Diffraction/Alchemi/*.cs`, and `ComputeIonizationU` / `FillIonizationUq` plus the STEM-EDX path in `Diffraction/BetheMethod.cs`.

- Provenance, licence, binary layout and version history of the data (`IonizationFsE0.bin`, `BoteSalvat.bin`) are documented in
  [`Diffraction/Ionization/README.md`](../../../Diffraction/Ionization/README.md) (Japanese, authoritative). Not repeated here.
- The Bloch-wave eigenproblem (`getEigenMatrix`, definitions of `P_g`, `Q_g`, `getU`): [bloch-wave.md](bloch-wave.md).
- The imaginary (absorptive) potential (TDS absorption form factors `FactorImaginary*`): [scattering-factors.md](scattering-factors.md).

---

## 1. Symbols and units

| Symbol | Meaning | Units / convention | Code |
|---|---|---|---|
| $k_\mathrm{vac}$ | vacuum electron wavenumber $=1/\lambda$ (no $2\pi$) | nm⁻¹, relativistic: $\lambda = 1.22643/\sqrt{V(1+0.97848\times10^{-3}V)}$ ($V$ in volts) | `UniversalConstants.Convert.EnergyToElectronWaveNumber` |
| $\mathbf G$ | momentum transfer (integer combination of reciprocal vectors) | nm⁻¹; **rebuilt canonically from the combined hkl** as `mat * index` (no floating-point vector sums) | `AlchemiMuBuilder` ctor, `ComputeIonizationU` |
| $s$ | $s=\lvert\mathbf G\rvert/2$ (sinθ/λ-type; **not** $q=4\pi\sin\theta/\lambda$) | internal API in **nm⁻¹**, F table in **Å⁻¹**; the only conversion is `sA = s * 0.1` in `IonizationTableShape.Evaluate` | `INormalizedIonizationShape` |
| $\sigma$ | total (sub)shell ionization cross section | nm² (`BoteSalvat.SigmaNm2` = cm² × 10¹⁴) | `IonizationData.TotalCrossSectionNm2` |
| $F(s)$ | normalized ionization form factor, $F(0)=1$, signed | dimensionless | `IonizationData.Shape` |
| $E_0$ | incident energy | keV, **30–400 keV only** (no extrapolation or clamping) | `IonizationDataProvider.MinE0KeV/MaxE0KeV` |
| $V_c$ | unit-cell volume | nm³ | `Crystal.Volume` |
| $M_a(\mathbf G)$ | Debye–Waller exponent | isotropic: $M=B_\mathrm{iso}s^2$ ($B$-type, nm²); anisotropic: $M=\beta_{11}H^2+\beta_{22}K^2+\beta_{33}L^2+2\beta_{12}HK+2\beta_{23}KL+2\beta_{31}LH$ ($\beta$ dimensionless; the index is transformed per equivalent position by `atom.Operation.ConvertPlaneIndex`) | `AlchemiMuBuilder.Build`, `ComputeIonizationU` |
| $\gamma_j$ | Bloch-wave eigenvalue (wavenumber correction along the surface normal) | nm⁻¹; $\mathrm{Im}\,\gamma_j>0$ decays ($e^{2\pi i\gamma_j z}$) | `AlchemiReduction` ctor docs |
| $z,t$ | depth / thickness along the surface normal | nm | |

The phase factor **$e^{-2\pi i\mathbf G\cdot\mathbf r_a}$ is generated directly** (the `getU` habit of building the positive phase and conjugating is not carried into the ionization path, 260801Cl).

---

## 2. Atomic data: σ and F(s)

### 2.1 Total cross section σ (Bote–Salvat)

`BoteSalvat.SigmaCm2(z, subshell, energyEv)`, subshell index 1=K, 2=L1, 3=L2, 4=L3, 5..9=M1..M5. Overvoltage $U=E/E_\mathrm{edge}$ (Bote edge unless `edgeEvOverride`).

- $U\le1$: $\sigma=0$.
- $1<U\le16$ (low-energy fit):
  $$\xi=(U-1)\left[\frac{a_1+a_2U+\frac{1}{1+U}\left(a_3+\frac{1}{(1+U)^2}\left(a_4+\frac{a_5}{(1+U)^2}\right)\right)}{U}\right]^2$$
- $U>16$: $\beta^2=E(E+2mc^2)/(E+mc^2)^2$, $x=\sqrt{E(E+2mc^2)}/mc^2$,
  $$\xi=\frac{A_{nlj}}{\beta^2}\frac{U}{U+B_e}\Big[(2\ln x-\beta^2)(1+g_1/x)+g_2+g_3\sqrt{mc^2/(E+mc^2)}+g_4/x\Big]$$
- $\sigma=4\pi a_0^2\,\xi$ with $a_0=5.291772108\times10^{-9}$ cm, $mc^2=5.10998918\times10^5$ eV (same as xion.f).

The operation order matches the Python reference `botesalvat.py` exactly (golden vectors). References: Bote & Salvat, PRA 77 (2008) 042701; xion.f (ADNDT 95 (2009) 871). The uncertainty of σ itself (K 10 % / L 15 % / M 24 %) is stated in the README.

### 2.2 Interpolating F(s, E0) (`IonizationChannelTable.GridAt`)

Each channel (shellCode, Z) stores E0 rows × 321 s nodes ($s=0..16$ Å⁻¹, step 0.05). Interpolation is two-stage, **E0 first, then s**:

1. **E0** (per s node $j$): abscissa $x=\ln(u-1)$, $u=E_0/E_\mathrm{th}$ (the row $u$ values are the serialized, 4-digit-rounded knots), PCHIP (scipy `PchipInterpolator`-compatible, class `Pchip`).
   - If all rows are positive at that node, interpolate $\ln F$ and exponentiate; otherwise interpolate the signed $F$ directly.
   - **Only rows with $s_\mathrm{cert}\ge s_j$ enter the basis** (260810Cl): low-E0 rows cannot reach $s_j$ kinematically and hold zero padding. If rows were dropped and the query E0 lies outside the remaining rows (= extrapolation), or fewer than 2 rows remain, the node is set to 0 (no exception).
   - $F(0)=1$ is imposed exactly.
2. **s** (`IonizationTableShape`): the 321 values are PCHIP-interpolated as signed $F$. When $s_\mathrm{cert}<16$ the **support of the interpolant is cut at $s_\mathrm{cert}$** so padding cannot contaminate the last interval.

Kinematic limit: $s\le k$ ($\lvert\mathbf G\rvert\le2k$), i.e. $1/\lambda$ in Å⁻¹; 14.33 Å⁻¹ at 30 kV (code comment).

### 2.3 Truncation and error bound (formatVersion ≥ 4, dataset v5 and later)

- $s\le s_\mathrm{cert}$: interpolated table value.
- $s>s_\mathrm{cert}$: **returns $F=0$ and declares a bound $\varepsilon$** (`TruncatedBeyondSMax=true`, `TruncationBound=ε`); no extrapolation.
- Between E0 rows: $\varepsilon=\max(\varepsilon_i,\varepsilon_{i+1})$, $s_\mathrm{cert}=\min(\cdot)$ — **ε is never interpolated** (an interpolated bound is not a bound; `TryGetTailModel`).
- Composite shells: bound $\sum_i w_i\varepsilon_i$ (triangle inequality, `IonizationLTotalShape.TruncationBound`).
- Code for the old exponential tail $F(s_\mathrm{max})e^{-\hat b(s-s_\mathrm{max})}$ remains but is unreachable with the bundled data (`IsBoundedTail`).

### 2.4 Composite shells (LTotal / MTotal) and line series (Kα, Kβ, Lα, Lβ, Mα)

`IonizationDataProvider.Resolve` → `IonizationLTotalShape` (named after L but used for every composite).

- Composite shells: $\sigma=\sum_i\sigma_i$, $F=\sum_iw_iF_i$, $w_i=\sigma_i/\sum\sigma$; each $\sigma_i$ is evaluated **at its own subshell edge** (closed subshells drop out as 0). LTotal = L1+L2+L3 (j-resolved dataset). MTotal uses **only the subshells present in the table** (Z = 30–32 lack M4/M5, so M1–M3).
- `EdgeEnergyKeV` = minimum of the **table** $E_\mathrm{th}$ of the constituents; overvoltage $U=E_0/E_\mathrm{edge}$; $U<1.2$ triggers a GUI caution (`LowOvervoltage`).
- Line series (260925Cl, `XrayLineSeries`): photon-production cross section $\sigma_s=\sum_iB_{si}\sigma_i$, shape $F_s=\sum_i(B_{si}\sigma_i/\sigma_s)F_i$,
  $$B_{si}=\sum_{\ell\in s}P_{k(\ell),i}\,\omega_{k(\ell)}R_\ell$$
  where $P_{ki}$ is the number of vacancies in shell $k$ produced by one primary vacancy in shell $i$ after the cascade, solved from inside out as $P_{ki}=\delta_{ki}+\sum_{j<k}t_{kj}P_{ji}$; transfers $t_{kj}$ = radiative $\omega_jR_{(k\leftarrow j)}$ + Coster–Kronig $f$ (L3←L1 uses $f_{13}+f'_{13}$). **No Auger cascade** (the bundled xraylib build has the Auger rates removed). $\omega$, $R$, $f$ from xraylib 4.2.1.
  - Vacancies moved by the cascade **keep the shape $F_i$ of the primary shell** (local approximation, per comment).
  - Feeding shells absent from the table (e.g. K for Z>50) are omitted (recorded in the provenance string).
  - Series: Kα=KA_LINE, Kβ=KB_LINE, Lα=LA_LINE (L3), Lβ = the 13 lines of xraylib LB_LINE_MACROS, Mα = M5N7+M5N6.
- Quantity: shell channels give `IonizationVacanciesGenerated`, line series `XrayPhotonsGenerated` (generated into 4π; no self-absorption, no detector efficiency).

---

## 3. Local ionization "potential" matrices

Both paths use the **local form-factor approximation** (single-variable $F(s)$, $s=\lvert\mathbf G\rvert/2$, $\mathbf G=\mathbf g_h-\mathbf g_g(+\mathbf q)$). A two-momentum MDFF kernel $K(Q_\mathrm{bra},Q_\mathrm{ket})$ is not supported (`AlchemiModelTier.FullMdff` → `NotSupportedException`). No dipole approximation is applied in this code — all shape information comes from the $F(s)$ table (generated by Temari, which is authoritative for its content).

The DWF is the **thermal average over the difference vector**, $\langle e^{2\pi i(\mathbf g-\mathbf h)\cdot\mathbf u}\rangle=e^{-M(\mathbf g-\mathbf h)}$. The product form $e^{-M(g)}e^{-M(h)}$ is forbidden (it makes the diagonal $\sigma e^{-2M}<\sigma$; header of AlchemiMuBuilder). At $\mathbf G=0$, $e^{-M}=1$ so **the diagonal is always $\sum_a\mathrm{Occ}_a\sigma$** (thermal vibration does not remove the total cross section). The $(1-e^{-M})$ structure of TDS absorption is not used.

### 3.1 STEM-EDX: `BetheMethod.ComputeIonizationU`

$$U_\mathrm{ion}(\mathbf G)=\frac{k_\mathrm{vac}}{2\pi V_c}\,\sigma\,F(\lvert\mathbf G\rvert/2)\sum_{a\in Z}\mathrm{Occ}_a\,e^{-M_a(\mathbf G)}\,e^{-2\pi i\mathbf G\cdot\mathbf r_a},\qquad\mathbf G=\mathbf q+\mathbf g_i-\mathbf g_j$$

- Sum over **all sites whose atomic number matches** (element filter); equivalent positions are already expanded in `atoms.Atom`.
- No relativistic $\gamma$ factor (the E0 dependence of σ and F already contains it, per comment).
- Meaning of the prefactor (Claude's addition, to be confirmed): with the $k=1/\lambda$ convention an imaginary potential $U'_0$ attenuates intensity with coefficient $\mu_\mathrm{abs}=2\pi U'_0/k$. Setting this to $n\sigma$ ($n$ atoms per volume) gives $U'_0=kn\sigma/(2\pi)$, i.e. $k\sigma/(2\pi V_c)$ per atom; the later factor $2\pi/k_\mathrm{vac}$ (§5) cancels it, leaving $\sigma/V_c$ (ionization probability per unit length).
- `FillIonizationUq`: for one $\mathbf q$, `dest[off + j*bLen + i] = U(q + g_i − g_j)` (column-major), with a per-channel combined-hkl cache (`ConcurrentDictionary`). Values are deterministic functions of the canonical $\mathbf G$, so the winning thread does not matter.

### 3.2 ALCHEMI: `AlchemiMuBuilder.Build`

$$\mu_{hg}=\sigma\,F(\lvert\mathbf G\rvert/2)\sum_{a\in\text{site}}\mathrm{Occ}_a\,e^{-M_a(\mathbf G)}\,e^{-2\pi i\mathbf G\cdot\mathbf r_a},\qquad\mathbf G=\mathbf g_h-\mathbf g_g\quad[\mathrm{nm}^2]$$

Deliberate differences from STEM-EDX (file header):
- **No element filter**: geometry of the site hypothesis (`AlchemiSiteBasis.AtomsIndices`: $\mathbf r_a$, ADP, Occ) × the dopant channel (σ, F) — the tracer combination.
- No $k_\mathrm{vac}/(2\pi V_c)$ ("cross-section matrix" in nm²); $V_c$ is divided out in the yield.
- Occupancy can be overridden with `OccupancyFraction` (the Tracer basis passes 1.0 and applies real occupancies by linear combination afterwards).
- Displacement via `SiteDisplacement`: `InheritHost` (default, host `Dsf`) / `IsotropicB` ($B_\mathrm{iso}$, nm²) / `AnisotropicB` ($\beta$, dimensionless).
- `AlchemiModelTier.StaticDelta`: $F\equiv1$, DWF$\equiv1$ internal oracle (for testing the reduction; not exposed in the GUI).

Numerics: $\Delta\mathbf G$ depends only on the basis, so it is made unique once with a pair→slot table (for $n=200$, 40 000 pairs map to a few thousand unique differences). $F(s)$ is batch-evaluated on the unique set and shared by all (site × channel). μ is orientation-independent (phases use integer hkl and fractional coordinates; $\lvert\mathbf G\rvert$ is rotation-invariant), so with a FixedUnion basis it is built **once** per (site × channel).

Diagonal: `Mu00` $=\mu_{00}=\sigma\sum_a\mathrm{Occ}_a\cdot m_a$ ($m_a$ = orbit multiplicity `atoms.Atom.Length`).

---

## 4. Ionization yield from Bloch waves (ALCHEMI, `AlchemiReduction`)

### 4.1 Wavefield and local ionization rate

$(\gamma_j, C^{(j)}_g)$ come from the eigendecomposition of `getEigenMatrix` (linearized Bethe equation with rows divided by $P_g$, after the 260829Cl fix). The boundary condition $\psi_g(0)=\delta_{g0}$ gives $\boldsymbol\alpha=C^{-1}\mathbf e_0$ (LU).

$$\psi_g(z)=\sum_j\alpha_jC^{(j)}_ge^{2\pi i\gamma_jz},\qquad\frac{dP}{dz}=\psi^\dagger(z)\,\mu\,\psi(z)$$

With $D=C\,\mathrm{diag}(\alpha)$ and $S=D^\dagger\mu D$ ($S_{j'j}=\sum_{h,g}D^*_{hj'}\mu_{hg}D_{gj}$) the thickness integral is analytic:

$$P(t)=\mathrm{Re}\sum_{j,j'}S_{j'j}F_{jj'}(t),\qquad F_{jj'}(t)=\frac{e^{\lambda_{jj'}t}-1}{\lambda_{jj'}},\quad\lambda_{jj'}=2\pi i(\gamma_j-\gamma^*_{j'})$$

- $\lvert\lambda\rvert^2<10^{-30}$ (nm⁻²) is treated as degenerate, $F=t$ (L'Hôpital).
- $e^{\lambda_{jj'}t}=u_j(t)u^*_{j'}(t)$ with $u_j=e^{2\pi i\gamma_jt}$, so only $b$ exponentials per thickness (not $b^2$).
- **Negative values are not clamped** (so truncation, sign and normalization errors stay visible).
- In the StaticDelta limit $\psi^\dagger\mu\psi=\sigma\sum_n\lvert\psi(\mathbf r_n)\rvert^2$, the same quantity as $S=B^\dagger\mathrm{diag}(\sigma)B$ in EBSD's `EBSDSolverManaged` (used as a regression oracle).

### 4.2 Dynamical and dechannelled terms

$$Y_\mathrm{dyn}(t)=\frac{P(t)}{V_c},\qquad Y_\mathrm{dech}(t)=\frac{\mu_{00}}{V_c}\big\{t-L_\mathrm{coh}(t)\big\},\qquad Y=Y_\mathrm{dyn}+Y_\mathrm{dech}$$

$$L_\mathrm{coh}(t)=\int_0^t\sum_g\frac{P_g}{P_0}\lvert\psi_g(z)\rvert^2dz,\qquad P_g=2\hat{\mathbf n}\cdot(\mathbf k_0+\mathbf g)$$

- $L_\mathrm{coh}$ uses the same reduction with $\mu=\mathrm{diag}(P_g/P_0)$, computed once per orientation and cached (`CoherentPathLengthNm`). The $P_g/P_0$ weight was added in 260829Cl (total normal current, exactly 1 without absorption); with `pOverP0=null` the old $\sum\lvert\psi_g\rvert^2$ is used.
- Meaning: electrons removed from the coherent field by TDS absorption stay in the specimen with essentially full energy, become randomly directed, and keep ionizing over the remaining thickness at the random-orientation rate $\mu_{00}/V_c$ (ICSC 2003 eq. (1), per comment).
- Units: $Y$ = **ionizations per incident electron (dimensionless)** (`SignalNormalization.PerIncidentElectron`).
- Guard: if $t-L_\mathrm{coh}<-10^{-6}\max(t,10^{-12})$ an `InvalidOperationException` is thrown (fail, no clamping); only smaller negative rounding is set to 0.
- **Re-injection contract (`AbsorptionSource`)**: $t-L_\mathrm{coh}$ picks up all attenuation produced by the imaginary part, whatever its origin. Re-injection is justified only when the imaginary part is pure TDS. `BetheMethod.ImaginaryPotentialAbsorption` declares `TdsRedistributable`, and `Yield` throws unless `IsFullyRedistributable()` (no bit other than TDS set). Mixing in `TrueLoss` / `Phenomenological` would require splitting the absorptive matrix by source and re-injecting only the TDS part. The declaration itself relies on discipline (update it whenever the imaginary part changes).
- Dynamic / Dechannelled / Total are **stored separately** ($t-L_\mathrm{coh}$ depends on thickness and orientation, so it cannot be absorbed into a constant nuisance parameter).

### 4.3 Run flow (`BetheMethod.RunAlchemi`, `BetheMethod.Alchemi.cs`)

1. Resolve channels (`IonizationDataProvider.Resolve`; immutable during the run). `ChannelDataOverride` is for sensitivity tests only (put δF on F; since the Bloch solution does not depend on F, the difference of two runs isolates the δF contribution exactly).
2. **FixedUnion basis**: call `Find_gVectors` at every scan orientation and take the true union of hkl (a single reference orientation drops g's with orientation-correlated errors that a fit would mistake for site-specific ICP). Ordered by $\lvert g\rvert^2\to h\to k\to l$, 000 first. `BasisOrientationIndices` can restrict the source orientations (to reproduce fixed literature bases; diagnostic only).
3. Build μ once per (site × channel).
4. Orientation loop (`Parallel.For`): `reset_gVectors` → `getEigenMatrix` (potential matrix `getPotentialMatrix(union)` shared across orientations, local form) → EVD (native Eigen, or MathNet `Evd(Asymmetric)`) → α (LU) → `AlchemiReduction.Yield`.
5. Expanded-basis check: re-solve the centre and both ends with a union built at `MaxNumOfBloch × ExpandedBasisFactor` (default 1.25); metric = max |ΔTotal| / max over the whole tensor; threshold `ExpandedBasisTolerance` default 3×10⁻³.

Defaults (`AlchemiRequest`): `Surface=(0,0,−1)`, `ModelTier=LocalFormFactor`, `MaxNumOfBloch=300`, `IncludeDechannelledComponent=true`, `ExpandedBasisFactor=1.25`, `ExpandedBasisTolerance=3e-3`, `UseNativeSolver=true`. Orientation is expressed as a beam direction (changing BaseRotation would rebuild gCache). `AlchemiScan.Rocking1D` rotates evenly about the tilt axis after removing its component parallel to $\hat d_0$.

Warnings (`BuildFixedUnion`): the basis needs $s_\mathrm{max}=\max\lvert g\rvert$ (Å⁻¹) above 16; tilt span > 30 mrad (outside the v1 FixedUnion guarantee) or 10–30 mrad (union + expanded-basis mandatory); union > 4× the centre-only basis; beams with $P\le0$.

`AlchemiBasisDiagnostic.Eligibility` **always returns `NotEvaluated`** (author decision 260809Cl). Reasons listed in code: ⑤ the denominator is the maximum of the whole tensor, diluting small-σ channels; ⑥ the numerator is an absolute yield and counts common scale changes that ICP cancels; ⑦ where Find_gVectors is limited by Ewald-sphere screening, 1.25× does not enlarge the basis, so the check compares the basis with itself (false positive).

### 4.4 Angular spread (`AlchemiAngularSpread.Gaussian`)

Post-processing along the orientation axis (the engine is untouched). The order is fixed: **convolve, then normalize**.
$$\tilde Y_i=\frac{\sum_jw_{ij}Y_j}{\sum_jw_{ij}},\quad w_{ij}=\exp\!\left[-4\ln2\,\frac{(\theta_i-\theta_j)^2}{\mathrm{FWHM}^2}\right],\quad\lvert\theta_i-\theta_j\rvert\le1.7\,\mathrm{FWHM}$$
The kernel is renormalized at the ends (treating outside points as 0 would create spurious dips); also valid for non-uniform spacing. FWHM ≤ 0 or fewer than 2 points returns the input.

---

## 5. STEM-EDX image (inside `BetheMethod.StemDoWork`)

The same slice integration as the STEM inelastic (TDS) path, with $U$ replaced by $U_\mathrm{ion}$ (shared in `AccumulateInelasticQ`, 260802Cl).

- q-space quantity, per reciprocal-lattice $\mathbf q$, thickness $t$, defocus $d$:
  $I(\mathbf q,t,d)=\sum_{\text{slices}\le t}\frac{2\pi}{k_\mathrm{vac}}\Delta z\cdot[\text{reduction in native }\texttt{STEM\_InelasticQ}]$.
  The reduction is "GEMM $W_q=U_q\times TC$, then 4 dots per entry (incident direction $k$, 4-point bilinear interpolation, lens function), accumulated over $d$" (comment). The factor $2\pi/k_\mathrm{vac}$ cancels the $k_\mathrm{vac}/2\pi V_c$ of $U_\mathrm{ion}$.
- ±q handling (order per design §3.4): independent computation → record the residual $\max\lvert I(\mathbf q)-I(-\mathbf q)^*\rvert/\max\lvert I\rvert$ and the $q=0$ imaginary part → if above tolerance, **fail without symmetrizing** → otherwise $\bar I(\mathbf q)=(I(\mathbf q)+I(-\mathbf q)^*)/2$, $q=0$ forced real. Fail if $\mathrm{Re}\,I(q=0)<-10^{-9}\max\lvert I\rvert$.
  - `HermitianTolerance` default 0.01, capped at 0.01 (tightening only). The residual scales roughly as $O(h^2)$ in the direction-grid division (measured 0.107/0.018/0.0093/0.0017 at div 10/20/32/48); `RecommendedProbeDivision=48`.
  - Environment variable `RECIPRO_STEM_EDX_HERMTOL` lifts the cap (escape hatch for aberration studies, 260805Cl).
- Real space: $\mathrm{image}(\mathbf R)=\mathrm{Re}[\sum_\mathbf qI(\mathbf q)e^{2\pi i\mathbf q\cdot\mathbf R}]/N$ with $N$ = **number of qEntries at $q=0$ = number of incident directions inside the aperture** (260802Cl fix; the old `radiusPix²` was ~2.7× too large and division-dependent). The real part, not `.Magnitude`, is used (keeps linearity and absolute scale). Source size is applied as a Gaussian blur. Finally negatives are clamped to 0 (pre-clamp minimum in `MinPixelBeforeClamp`).
- Constraints: native helpers required (no managed implementation → `NotSupportedException`); duplicate (Z, Shell) requests are a hard error.
- Diagnostics exposed in `StemSignalMap`: `HermitianResidualMax`, `QZeroImagMax`, `MinPixelBeforeClamp`, `TruncatedBeyondSMax`, `TruncationBound`, `Iq` (symmetrized), `IqBeforeSymmetrization` (with `CaptureRawIq`).

---

## 6. Approximations and range of validity (summary)

- Isolated-atom (DHFS/Dirac) F(s), local single-variable form factor; no bonding/solid-state effects.
- σ from Bote–Salvat (uncertainty in README); E0 = 30–400 keV.
- Cascade vacancies keep the primary-shell shape; no Auger cascade.
- ALCHEMI dechannelled electrons: energy conserved, direction fully randomized, ionize at the random-orientation cross section; non-TDS absorption cannot be re-injected.
- Self-absorption and detector efficiency (`XrayPhotonsSelfAbsorbed/Detected`) exist only as enum values; not implemented.

---

## Possible implementation issues

1. **Wrong comment on the truncation width of the angular spread** (`Diffraction/Alchemi/AlchemiAngularSpread.cs:48–50`). It says weights beyond "4σ (≈1.7 FWHM)" are < 1e-7, but $w(1.7\,\mathrm{FWHM})=\exp(-4\ln2\cdot1.7^2)=e^{-8.01}\approx3.3\times10^{-4}$; the discarded one-sided mass is ≈ 3×10⁻⁵ of the total (Claude's estimate). Small numerical effect, but reaching the stated 1e-7 would need ≈ 2.4 FWHM (≈ 5.7σ).
2. **$L_\mathrm{coh}$ counts current, $Y_\mathrm{dyn}$ counts density** (`AlchemiReduction.cs:146–171` vs `:194`, 260829Cl). For $\mu=\mu_{00}I$, $Y=\frac{\mu_{00}}{V_c}[t+\int\sum_g(1-P_g/P_0)\lvert\psi_g\rvert^2dz]\ne\mu_{00}t/V_c$. With the surface normal ⟂ $\mathbf g$ (ZOLZ) $P_g=P_0$ and the difference vanishes; it appears only for HOLZ beams or inclined surfaces. Which is physically right (oblique components have a longer path per unit depth) needs a decision.
3. **No $2K/P_0$ (≈ $1/\cos\theta_0$) factor in $Y_\mathrm{dyn}$** (`AlchemiReduction.Yield`, `:196`). The incident current is $\propto P_0\lvert\psi_0\rvert^2$, so ionizations per incident electron should be $(n\sigma)(2K/P_0)\int\lvert\psi\rvert^2dz$ (Claude's derivation, to be confirmed); $Y_\mathrm{dech}$ likewise uses the normal thickness $t$ as the path. About 4.5×10⁻⁴ over a 30 mrad scan, but it matters when the beam is strongly inclined to the surface normal.
4. **Diagnostic flags in `IonizationTableShape` are shared mutable state** (`IonizationChannel.cs`, `_truncated/_usedTail`). STEM evaluates the same shape concurrently from `Parallel.For` (writes only `true`, so benign), but if an `IonizationData` is reused across runs (e.g. via `ChannelDataOverride`), a flag set by an earlier run persists and `TruncatedBeyondSMax` can be reported spuriously.
5. **E0-direction PCHIP extrapolation at the table ends** (`IonizationChannelTable.GridAt`). At nodes where no row is dropped ($n=$rows), `Pchip.Evaluate` extrapolates with the end cubic if the query E0 is outside the first/last row; the only check is 30–400 keV. Harmless if the table rows always include 30 and 400 keV, which cannot be verified without reading the `.bin`.
6. **`RECIPRO_STEM_EDX_HERMTOL=NaN` disables the check** (`BetheMethod.cs:3450–3452, 3481`): `double.TryParse("NaN")` succeeds and `hermMax > NaN` is always false, whereas the API `HermitianTolerance` rejects NaN before the run.
7. **Stale numbers in messages/comments**: exception text "s > 4" and class summary "s_max (=8 Å⁻¹)" in `IonizationTableShape`, "s ≤ 8 Å⁻¹" in `AlchemiBasisDiagnostic`, "(s ≤ 8 Å⁻¹)" in `BuildFixedUnion:186–189`, "81" in the `IonizationFsTable.SCount` note. Current: 321 nodes, 16 Å⁻¹.
8. **Comment at `AlchemiMuBuilder.cs:206`** says "same rule as Dsf.BisoEffective" but the code uses `dwf.Biso000` (`ComputeIonizationU` uses `BisoEffective`). At $\mathbf G=0$, $s^2=0$, so the values coincide and there is no effect (261105Cl change).

## Open questions (to be confirmed)

1. The formula inside the native kernel `STEM_InelasticQ` (Crystallography.Native, outside this repository); this document reconstructs it from C# comments.
2. Physical derivation of the $k_\mathrm{vac}/(2\pi V_c)$ normalization of $U_\mathrm{ion}$ (the explanation in §3.1 is Claude's).
3. Origin of the 0.6 % excess in the 260802Cl measurement "$I(q=0)$/entries/$(\sigma tN/V)$ = 1.00586" (was it measured before the 260829Cl $P_g$ row-division fix?).
4. Whether adding $f_{13}+f'_{13}$ for L3←L1 in the Coster–Kronig transfers double-counts, given xraylib's definition of $f'_{13}$ (Claude's addition, to be confirmed).
5. Size of the difference between the table $E_\mathrm{th}$ used for `EdgeEnergyKeV` of LTotal/MTotal and the Bote edges used for σ.
6. Placement of E0 rows in the F table (do they include 30 and 400 keV as end points?); related to issue 5.
7. How the design document (`.project-guidance/ReciPro/ReciPro_ALCHEMI設計.md` §3.4, outside this repository) settles issues 2 and 3 (flux normalization and $1/\cos\theta_0$).
8. Whether `XrayLineSeries` Kβ = KB_LINE "= 1 − Kα" matches xraylib's RadRate definition.
9. That the TDS path's `getU(...).Imag.Conjugate()` ("taking the conjugate somehow works") and the ionization path's "direct negative phase" follow the same phase convention (to be checked together with [bloch-wave.md](bloch-wave.md)).
