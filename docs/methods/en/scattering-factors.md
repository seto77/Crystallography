This is a draft written by Claude from the code (2026-11-05). Not yet reviewed by the author.

# Atomic scattering factors, temperature factors and TDS absorptive form factors

Scope: `Atom/AtomStatic.cs` (`AtomStatic.ES`, Temari factors, neutron scattering lengths), `Atom/Atoms.cs` (`Atoms.GetAtomicScatteringFactorFor*`, `DiffuseScatteringFactor`), `Xraylib.cs` (f′/f″).
For how these enter the structure factor see [structure-factor.md](structure-factor.md); for the Bloch-wave potential $U_g$ see [bloch-wave.md](bloch-wave.md).

---

## 0. Conventions (most important)

| Quantity | Symbol | Unit / definition in the code |
|---|---|---|
| Scattering-vector magnitude | $s=\sin\theta/\lambda=|g|/2=1/(2d)$ | **$q=4\pi s$ is never used** (only the Kirkland table uses $q_K=2s$, see below) |
| Argument of `ES.Factor(S2)` | $s^2$ | **nm⁻²** (callers pass `vec.Length2/4`, `1/d²/4`) |
| Return value of `ES.Factor` | — | **electron tables: nm; X-ray tables: electrons × 0.1** (corrected 260818Cl) |
| `ES.ElasticFactorWithTail(s2)` | $s^2$ | **Å⁻²**, returns Å (a different convention inside the same class) |
| `TemariFactor.Fx(s)`, `Fe(s)` | $s$ | **Å⁻¹**, returns electrons / Å |
| Temperature factor | $B$ | **nm²**. `DiffuseScatteringFactor` converts between B and U types ($B=8\pi^2U$) |
| Anisotropic | $\beta_{ij}$ | `B11`…`B31` are the **dimensionless β** ($T=\exp(-\mathbf h^{\mathsf T}\beta\mathbf h)$) |
| γ of electron factors | — | **None of the tables (Peng / 8-Gaussian / Kirkland / Temari) contains γ.** γ is applied in `BetheMethod.getU` |
| Angles | — | radians internally (detector inner/outer are also radians) |

### Unit conversion in practice (`ES` constructors, AtomStatic.cs L2779–2822)
The Gaussian parameters are in Å units ($a_i$ [Å] or electrons, $b_i$ [Å²]). `Factor` is

$$
\texttt{Factor}(S_2)=0.1\times\Big[\sum_i a_i\exp\!\big(-b_i\cdot 0.01\,S_2\big)\;(+c)\Big],\qquad S_2\,[\mathrm{nm^{-2}}]
$$

where $0.01$ converts nm⁻² → Å⁻² and $0.1$ converts Å → nm. The same ×0.1 is applied to the X-ray tables (with constant $c$), so X-ray values come out in the artificial unit "electrons × 0.1". `FormBeamInteraction` multiplies by 10 for X-rays to get electrons back (comment at L2734–2737).

---

## 1. Elastic scattering-factor tables (`AtomStatic` static fields)

| Field | Form | Variable | Source (as cited in the code) | Use inside this library |
|---|---|---|---|---|
| `XrayScatteringWK[Z][sub]` | 5 Gaussians + c | $s$ [Å⁻¹] | Waasmaier & Kirfel (1995) Acta Cryst. A51, 416–431 | **kinematic X-ray F** (`Atoms.GetAtomicScatteringFactorForXray`), $f_x$ of the Mott–Bethe tail |
| `XrayScattering[Z][sub]` | 4 Gaussians + c | $s$ | Cromer & Mann (1968) Acta Cryst. A24, 321–324 | unused in this library (no references) |
| `ElectronScatteringPeng[Z][sub]` | 5 Gaussians | $s$ | Peng et al. (1996) Acta Cryst. A52, 257–276 (neutral) / Peng (1998) A54, 481–485 (ions) | **real part of the elastic $U_g$**, kinematic F, integrand of the TDS absorption |
| `ElectronScatteringEightGaussian[Z]` | 8 Gaussians | $s$ | private communication from Dr. Ohtsuka (unpublished) | unused (commented out in `getU`); only `AttachTail` is called on it |
| `ElectronScatteringKirkland[Z]` | 3 Lorentzians + 3 Gaussians | $q_K=2s$ [Å⁻¹] | Kirkland, *Advanced Computing in Electron Microscopy*, App. C (fparams.dat) | `ProjectedPotential` only, from the obsolete `DynamicTheory.MultiSlice` |

`[Z][0]` is always the neutral-atom entry (shared convention of the `AtomStatic` static constructor and the `ElasticIonModel.Neutral` branch of `getU`). `sub ≥ 1` are ions (`ES.Valence` holds the charge).

### 1.1 Peng 5-Gaussian (electrons)

$$
f_e(s)=\sum_{i=1}^{5}a_i\,e^{-b_i s^2}\quad[\text{Å}],\qquad s\ [\text{Å}^{-1}]
$$

- Note corrected at 260923Cl: the coefficients are ITC Vol. C Table 4.3.2.2 = Peng 1996 fit **for $s\le 2$ Å⁻¹** (the $s\le6$ fit is the separate Table 4.3.2.3). $s>2$ is an extrapolation.
- Measurement at 260811Cl (`AlchemiCheck absorb -join`): relative to Kirkland, Peng deviates by more than 5 % from $s\approx2.1$ Å⁻¹, is 0.09 at $s=4$ and $10^{-7}$ at $s=8$. The Gaussian sum goes exactly to zero at high $s$ and lacks the Mott–Bethe $1/s^2$ tail.
- **The real part of the elastic $U_g$ is still pure Peng** (no tail). The comment (L2918–2921) explains: practical reflections have $|g|\le 50$ nm⁻¹ ($s\approx0.25$ Å⁻¹), far from the joining region, and changing it would move all existing results.

### 1.2 Kirkland (electrons, 3 Lorentzians + 3 Gaussians)

$$
f_e(q_K)=\sum_{i=1}^{3}\frac{a_i}{q_K^2+b_i}+\sum_{i=1}^{3}c_i e^{-d_i q_K^2},\qquad q_K=2s=1/d\ [\text{Å}^{-1}]
$$

The implementation (L2830–2841) converts $q_K^2[\text{Å}^{-2}]=0.04\,S_2[\text{nm}^{-2}]$ and multiplies by 0.1 to get nm (fixed at 260606Cl; previously all conversions were missing). Projected potential (Kirkland eq. C.20):

$$
v_z(r)=2\pi^2 a_0 e\Big[\sum_i 2a_i K_0\!\big(2\pi r\sqrt{b_i}\big)+\sum_i\frac{c_i}{d_i}e^{-\pi^2r^2/d_i}\Big]\quad[\text{V·Å}],\ r\ [\text{Å}]
$$

$a_0e$ is `UniversalConstants.a0·1e10 × 2Ry·a0` (= 0.529 Å × 14.4 V·Å). Returns 0 for $r\le0$.

### 1.3 X-rays (Waasmaier–Kirfel / Cromer–Mann)

$$
f_x(s)=\sum_i a_i e^{-b_i s^2}+c\quad[\text{electrons}]
$$

Constructor arguments are interleaved $(a_1,b_1,a_2,b_2,\dots,c,\text{valence},\text{method})$; the original $(a_1..a_5,c,b_1..b_5)$ ordering of the paper has been rearranged on storage (checked with the H entry).

### 1.4 X-ray anomalous dispersion f′, f″ (`Xraylib.cs`)

- `Xraylib.Fprime(Z, E_keV)` = xraylib `Fi`; `Xraylib.Fdoubleprime(Z, E_keV)` = **−`Fii`** (xraylib's raw value has the opposite sign; this returns the crystallographic $f''>0$). E in keV, result in electrons.
- Self-test: Fe at 8.04778 keV gives $f'=-1.106$, $f''=+3.192$ (vs. ITC $-1.13$, $3.20$; Xraylib.cs L24).
- Activation is decided by `Crystal.IsXrayDispersionActive`: `anomalousDispersion && wave==Xray && E finite > 0 && Xraylib.Enabled`. **Depends only on Z and energy; ionic charge/chemical state is ignored.**
- How it enters the structure factor (including the sign convention): [structure-factor.md §2](structure-factor.md).

### 1.5 Neutron scattering lengths (`NeutronCoherentScattering[Z][i]`)

- Unit **fm**, complex. `[Z][0]` = natural abundance, `[Z][i≥1]` = isotopes (in the enumeration order of `IsotopeAbundance[Z]`). Source: periodictable 2.1.0 nsf (Sears 1992, Neutron Data Booklet; 260531Ch).
- The imaginary part is absorption at 1.798 Å: $\mathrm{Im}\,b=-\sigma_a/(2000\,\lambda_{th})$ (b [fm], σ [barn], λ [Å]). `NeutronAbsorptionCrossSection(Z, λ)` recovers $\sigma_a(\lambda)=-\mathrm{Im}\,b\cdot2000\,\lambda$ (1/v law). Cd/Sm/Eu/Gd break 1/v (`NeutronIsResonantAbsorber`).
- $\sigma_{coh}$ uses the tabulated value, never $4\pi|b|^2$ (they differ for resonant nuclei; Gd: 35.3 vs 29.3 barn).
- With isotopes specified (`Atoms.GetAtomicScatteringFactorForNeutron`): $b=\sum_i b_i\,x_i/100$ (`Isotope[i]` in %) × Occ. If the array length differs from `IsotopeAbundance[Z].Count`, the natural value is used. Independent of $s$.

---

## 2. Scattering factors from `Atoms` (occupancy included)

`Atoms.cs` L355–375:

| Method | Expression | Notes |
|---|---|---|
| `GetAtomicScatteringFactorForElectron(s2)` | `ElectronScatteringPeng[Z][SubNumberElectron].Factor(s2) * Occ` | s2 [nm⁻²], returns nm. **Uses the assigned ion (SubNumberElectron) and ignores `BetheMethod.ElasticIonModel`** |
| `GetAtomicScatteringFactorForXray(s2)` | `XrayScatteringWK[Z][SubNumberXray].Factor(s2) * Occ` | returns electrons×0.1 |
| `GetAtomicScatteringFactorForNeutron()` | §1.5 × Occ | fm |

**No temperature factor here** (applied in the structure factor / $U_g$).

---

## 3. Temperature factors (`DiffuseScatteringFactor`, Atoms.cs L442–)

Stored as `OriginalType ∈ {U, B}`, `Iso`, `Aniso11..31` (U type: nm²; B type: dimensionless β). Getters:

$$
B_{iso}=\begin{cases}\mathrm{Iso} & (\text{B type})\\ 8\pi^2\,\mathrm{Iso} & (\text{U type})\end{cases}\quad[\mathrm{nm^2}],\qquad
\beta_{ij}=\begin{cases}\mathrm{Aniso}_{ij} & (\text{B type})\\ 2\pi^2 a^*_i a^*_j\,U_{ij} & (\text{U type})\end{cases}
$$

- Fixed at 261105Cl: the off-diagonal coefficient is $2\pi^2a^*_ia^*_j$ (**no $\cos\gamma^*$ etc.**; CIF `_atom_site_aniso_U_ij`, Grosse-Kunstleve & Adams (2002) eqs. 4–6). Previously the cosines were included (so $U_{12}$ etc. vanished for orthogonal cells and $U_{12}$ was halved for hexagonal cells).
- Equivalent isotropic $B$ (Hamilton 1959 eq. 18; corrected to real-lattice dot products at 261105Cl):
  $$B_{eq}=\tfrac43\sum_{ij}\beta_{ij}\,(\mathbf a_i\cdot\mathbf a_j)$$
  (`Biso000`). For U-type input this equals $8\pi^2U_{eq}$.
- `BisoEffective` (for g = 0; changed 261105Cl): `IsZero`→0, `UseIso`→`Biso`, anisotropic→always `Biso000`. NaN and negative → 0.
- `BisoFloorNm2 = 1e-3` nm² (0.1 Å²): floor used for the EBSD backscatter source and MC. **Not applied to the absorption in `getU`.**
- Coefficients are built from `Dsf.Cell` (nm, rad). `Cell` is `[XmlIgnore]`, so `Crystal` sets it again after loading (Crystal.cs L711–719).
- For an anisotropic atom $j$ at an equivalent position $\mathbf r_j=R\mathbf r+\mathbf t$, the temperature factor uses $R^{\mathsf T}\mathbf h$ returned by `SymmetryOperation.ConvertPlaneIndex` (corrected comment 261105Cl, SymmetryOperation.cs L99):
  $$T_j(\mathbf h)=\exp\!\big[-(R^{\mathsf T}\mathbf h)^{\mathsf T}\beta\,(R^{\mathsf T}\mathbf h)\big]$$
  equivalent to $\beta'=R\beta R^{\mathsf T}$.
- `Transform` (β′ = MβMᵀ) and `Deform` (keeps the Cartesian ellipsoid under cell deformation, M = A⁻¹S⁻¹A) were added at 261105Cl. Errors are propagated treating the 6 components as independent (no covariances).

### How the isotropic factor is applied
Both the structure factor and $U_g$ use $T=\exp(-B s^2)$ ($B$ [nm²], $s^2$ [nm⁻²]). I checked that **there is no confusion with $\exp(-2B s^2)$ or $\exp(-8\pi^2U s^2)$** (Crystal.cs L1590, L1673; BetheMethod.cs L3921).

---

## 4. High-s tail of the electron factor (`ES.ElasticFactorWithTail`, L2957)

Elastic factor used **only inside the absorption (TDS) integrand**. Argument $s^2$ [Å⁻²], returns Å. Branches on the global setting `BetheMethod.AbsorptionElasticModel` (default `TemariHybrid`, author decision 260922Cl):

1. **TemariHybrid, neutral entry (Valence = 0), $s\le 6$ Å⁻¹, Z = 1–86** → `TemariScattering(Z).Fe(s)`.
2. Otherwise (ions, Z = 87–98, `PengTail`, $s>6$):
   - $s\le 1.5$: pure Peng (`PengSum`)
   - $s\ge 2.5$: Mott–Bethe
     $$f_e^{MB}(s)=0.023934\,\frac{Z-f_x(s)}{s^2}\ [\text{Å}],\qquad f_x\ \text{clamped to}\ [0,Z]$$
     $f_x$ is the WK factor of the **neutral atom of the same element** (also for ions: at high $s$ only the nuclear charge is seen). Coefficient `MottBetheCoefficientAngstrom = 0.2393366/10`.
   - $1.5<s<2.5$: smoothstep $w=t^2(3-2t)$, $t=(s-1.5)/1.0$, $f=(1-w)f_{Peng}+wf^{MB}$.
3. `useTail=false` is for verification only (pure Peng).

Rationale for the join (table in the 260811Cl comment): Mott–Bethe agrees with Kirkland to 0.97–1.00 everywhere, Peng deviates by 5 % at $s\approx2.1$. Rationale for the clamp: outside their fit range some WK fits diverge to $f_x<0$ (C: $a_4=-4.24$, N: $c=-11.80$). Rationale for Temari (260922Cl): the $s>6$ contribution is 1.5–4.8 % of $f''(0)$, the tail-model spread is 0.04–0.43 %, the effect on observables ≤ 0.45 %. The jump at $s=6$ is ≤ 0.041 % and is deliberately not matched.

Paths with the tail: `FactorImaginary`, `FactorImaginaryAnnular` (local and non-local). Paths without: `FactorImaginaryAnnularFlatEwald`, `FactorImaginaryAnnular2` (both unused in production, marked unfinished).

`AttachTail(z, neutralXray)` is called by the `AtomStatic` static constructor for all entries of the Peng and 8-Gaussian tables (L3335–3357), with `neutralXray = XrayScatteringWK[z][0]`. X-ray and Kirkland entries keep `AtomicNumber = 0` (tail disabled).

---

## 5. Temari factors (`AtomStatic.TemariFactor`, L7569–)

- Source (as cited): Temari dataset-factors v2.0.0 (Y. Seto), CC-BY-4.0, DOI 10.5281/zenodo.22820415. Full Dirac SCF + exact-exchange KLI (model_id `DHFS-KLI-DTM1-dt16-neutral-v1`). The dataset **makes no error-bound claim** (`certification_status = not_certified`).
- **Neutral atoms only, Z = 1–86, $0\le s\le 6$ Å⁻¹.** Out-of-range input throws (no extrapolation, no clamping; `Guard`).
- `Fx(s)` [electrons]: C² cubic spline in $s$, left end clamped $f_x'(0)=0$, right end not-a-knot.
- `Fe(s)` [Å]: spline **in $t=s^2$**, not-a-knot at both ends. **No γ** (`TemariGammaIncluded = false`). First-Born (Mott–Bethe) non-relativistic $f_e$.
- 7681 nodes ($s_i=6\,i/7680$). Stored with absolute quantisation ($f_x$ 1e-10 electrons, $f_e$ 1e-11 Å) + third differences + byte-plane shuffle + Brotli (`TemariFactors.bin`, magic "TMRF", format 1). Max deviation from the published values 5e-11 / 5e-12.
- `TemariSpline` is a verbatim port of the reference loader `tools/temari_factors_contract.py` (the comment says not to "simplify" the formulas). Interval search changed at 260922Cl to an equal-spacing guess plus correction (bit-identical values).
- **The only use inside this library is `Fe` in `ElasticFactorWithTail`.** `Fx` is not used for the elastic X-ray structure factor (that stays WK).

---

## 6. TDS absorptive form factor $f''$ (`ES.FactorImaginary*`)

Absorption enters as the **imaginary part of the potential**: $U_g=U^{(r)}_g+iU'_g$ (`BetheMethod.getU` returns the pair `(Real, Imag)`, each complex). Currently the only source of the imaginary part is TDS (`BetheMethod.ImaginaryPotentialAbsorption = AbsorptionSource.TdsRedistributable`).

### 6.1 Local form, full solid angle: `FactorImaginary(kV, s2, m)` (production, L3012)

Inputs: `s2` = $s_g^2$ [nm⁻²], `m` = $B$ [nm²] (returns 0 for NaN or 0). $\gamma=1+eV/(m_0c^2)$; $k_0$ = `EnergyToElectronWaveNumber` (**$k_0=1/\lambda$, relativistically corrected, nm⁻¹**).

$$
f''(\mathbf g)=\gamma\,\frac{k_0}{2}\int_{4\pi}\! d\Omega\;f_e(s_1)\,f_e(s_2)\Big[1-e^{-B\,(s_1^2+s_2^2-s_g^2)}\Big]\times 0.01
$$

$$
s_{1,2}=\tfrac12\big|\mathbf K\mp\tfrac{\mathbf g}{2}\big|,\qquad \mathbf K=k_0(\sin\theta\cos\phi,\ \sin\theta\sin\phi,\ \cos\theta-1)
$$

- $f_e$ is `ElasticFactorWithTail` (Å); ×0.01 converts Å² → nm². The result has the unit **nm** (same as $f_e$). The doc comment says "dimensionless", but `getU` multiplies it by the same $\gamma/(\pi V)$ as the real part, so effectively nm (→ implementation question Q2).
- **$\mathbf g$ is placed perpendicular to the beam (x axis)** — only $|\mathbf g|$ is passed. The old closed form was also a function of $|g|$ only.
- Relation to the optical theorem (added by Claude, to be confirmed): with $k=1/\lambda$, $\mathrm{Im}f(0)=k\sigma/2$. With the relativistic $f=\gamma f_e$ a factor $\gamma^2$ appears; the code puts one $\gamma$ in $f''$ and another in `getU`, giving $\gamma^2$ in total, which is consistent.
- Numerical integration:
  - θ is split at **geometrically spaced $s=k_0\sin(\theta/2)$** `AbsorptionSEdges = [0, 0.3, 1, 3, 9, 27, 81, 243, 729, ∞]` nm⁻¹ ($\theta=2\arcsin(s/k_0)$; segments with $s\ge k_0$ end at θ=π; the trailing ∞ guarantees full solid-angle coverage even at MV energies).
  - **12** Gauss–Legendre points per segment (`AbsorptionThetaPoints`); φ uses the **8-point** rule on [0, π] folded at φ=π/2 into 4 nodes with weights ×4 (`PhiFolded`; the weight sum 2π is asserted), using the even symmetries φ→−φ and φ→π−φ.
  - Accuracy (260811Cl): against a 48×48 reference of the same integrand, within 0.05 % in $f''(0)$ ratio for $|g|\le50$ nm⁻¹ (12×8: 0.023 %, 16×16: 0.0024 %).
  - The DW term $1-\exp[m(s_g^2-2q_0)]$ is φ-invariant and evaluated once per θ ($q_0=(|\mathbf K|^2+s_g^2)/4$).
- Verification reference: `FactorImaginaryGaussianAnalytic` (old closed form, Gaussian-only) should agree with `useTail:false`.

Closed form (`FactorImaginaryGaussianAnalytic`, L3129, not used in production):

$$
f''=\frac{2\pi\gamma}{k_0}\sum_{i,j}a_ia_j\left[\frac{e^{-s_g^2\,b_ib_j/(b_i+b_j)}}{b_i+b_j}-\frac{e^{-s_g^2\,(b_ib_j-B^2)/(b_i+b_j+2B)}}{b_i+b_j+2B}\right]
$$

($s_g^2$ and $B$ converted to Å units: `s2*=0.01`, `m*=100`; terms with $b_i+b_j=0$ are 0.)

### 6.2 Local form, annular detector: `FactorImaginaryAnnular(kV, g, m, inner, outer, nTheta=60, nPhi=20)` (L3221)

Same integrand, θ ∈ [inner, outer] (**radians**), φ ∈ [0, 2π], Gauss–Legendre (default 60×20). All three components of $\mathbf g$ are used ($g_z\ne0$ allowed; no φ folding). Used for STEM inelastic (e.g. HAADF) potentials.

### 6.3 Non-local form: `FactorImaginaryAnnular(kV, g, h, m, inner, outer, …)` (L3274)

$$
f''(\mathbf g,\mathbf h)=\gamma\frac{k_0}{2}\int d\Omega\ f_e\!\big(\tfrac12|\mathbf K-\mathbf g|\big)f_e\!\big(\tfrac12|\mathbf K-\mathbf h|\big)\Big[1-e^{-B(s_{Kg}^2+s_{Kh}^2-s_{g-h}^2)}\Big]\times0.01
$$

Called from `getU(kV, g, h, inner, outer)`.

### 6.4 Unused variants
- `FactorImaginaryAnnularFlatEwald` ("unfinished") and `FactorImaginaryAnnular2` ("slow, abandoned?"). Both are hand-written Gaussian sums without the tail. → implementation question Q6.

### 6.5 Use in `getU` (overview; details in [bloch-wave.md](bloch-wave.md))
`BetheMethod.getU` (BetheMethod.cs L3863):

$$
U_g=\frac{\gamma}{\pi V_c}\sum_{\text{sites}}\mathrm{Occ}\,f_e(s_g)\,e^{-B s_g^2}\,e^{+2\pi i\mathbf g\cdot\mathbf r},\qquad
U'_g=\frac{\gamma}{\pi V_c}\sum\mathrm{Occ}\,f''(\mathbf g)\,e^{-Bs_g^2}\,e^{+2\pi i\mathbf g\cdot\mathbf r}
$$

- $V_c$ [nm³], $U$ [nm⁻²]. The comment says the coefficient follows "Kirkland's textbook p.120".
- The phase is **$e^{+2\pi i\mathbf g\cdot\mathbf r}$** ("20240524 sign changed to plus, Mitsubishi/Nakamura") — opposite to the kinematic F ($e^{-2\pi i}$; [structure-factor.md §4](structure-factor.md)).
- For anisotropic atoms, the direction-dependent effective value $m=\mathbf h'^{\mathsf T}\beta\mathbf h'/s_g^2$ ($\mathbf h'=R^{\mathsf T}\mathbf h$) is passed to $f''$ **as an isotropic B**. For $g=0$, `BisoEffective` is used.
- Ions, `ElasticIonModel` (default `Neutral`): Neutral = `[Z][0]`; IonScreened = finite Peng part of the assigned ion; IonFull = finite part + Mott–Bethe monopole $0.2393366\,\Delta Z/s^2$ [nm] ($s^2$ [nm⁻²]) for $g\ne0$. If the cell is not neutral ($|\sum q|>10^{-6}$) the monopole is dropped (260614Cl). The monopole enters only the real part, not the absorption integral.
- The cache `uDictionary` is keyed by the hkl key only, not by energy (`MeanAbsorptionLengthNm` clears it per energy).

---

## 7. NIST elastic scattering (link only)

`Atom/NistElastic.cs` and [`Atom/NistElastic/README.md`](../../../Atom/NistElastic/README.md) serve only the elastic-angle sampling of the EBSD Monte Carlo (`EBSD/MonteCarlo.cs`). I confirmed by grep that **they are not referenced from the diffraction paths (Bloch waves, structure factors)**.

---

## Implementation questions

- **Q1** The doc comment of `MottBetheMonopoleCoefficient` (`AtomStatic.cs` L1209–1213) says "$s^2$ [Å⁻²], f in [nm]", but the value 0.2393366 corresponds to **$s^2$ [nm⁻²] and f [nm]** ($0.023934\,\text{Å}\cdot\text{Å}^{-2}$ → $f[\text{nm}]=0.1\times0.023934\,\Delta Z/(0.01\,s^2_{nm})=0.23934\,\Delta Z/s^2_{nm}$). `getU` passes nm⁻², so the computation is right; the note at L2936 of the same file correctly says nm⁻². Comment-only inconsistency.
- **Q2** The `<returns>` of `FactorImaginary*` ("dimensionless") and the comment in `getU` ("the answer is dimensionless", BetheMethod.cs L3918) disagree with the actual dimension ($\gamma k_0/2\cdot\int f_ef_e\,d\Omega\times0.01$ = nm, the same as the real $f_e$). `getU` applies the same $\gamma/(\pi V)$ as for the real part, so treating it as nm is consistent. Probably a comment-only issue.
- **Q3** The kinematic F (`Crystal.GetStructureFactor` → `Atoms.GetAtomicScatteringFactorForElectron`) uses Peng for `SubNumberElectron` (finite part of the assigned ion = IonScreened-like) and does not consult `BetheMethod.ElasticIonModel` (default Neutral). For a crystal with ions assigned, the kinematic $|F|$ and the dynamical $U_g$ of the same reflection use different models. Intended?
- **Q4** TDS absorption with anisotropic temperature factors (BetheMethod.cs L3900–3918): the DW exponent $B(s_1^2+s_2^2-s_g^2)$ inside the $f''$ integrand uses the effective isotropic value $m=\mathbf h'^{\mathsf T}\beta\mathbf h'/s_g^2$ defined along $\mathbf g$. Since $\mathbf K\pm\mathbf g/2$ point in other directions, β should strictly be evaluated per direction. Acceptable as an approximation? (Exact for isotropic atoms.)
- **Q5** Neutron isotope weighting (Atoms.cs L366–371) relies on the enumeration order of `Dictionary<int,double>` matching the order of `NeutronCoherentScattering[Z][i+1]`. .NET's Dictionary enumerates in insertion order when nothing has been removed, as an implementation behaviour, but this is not guaranteed by the specification (added by Claude, to be confirmed).
- **Q6** `FactorImaginaryAnnularFlatEwald` (L3174) lacks the Jacobian $R$ in the radial integral (no multiplication by `R`; the comment still says "sinTheta taken outside"), and its prefactor $\gamma/(\pi^2k_0)$ does not match the $\gamma k_0/2\cdot(1/k_0^2)=\gamma/(2k_0)$ of the other variants. Unused in production ("unfinished"), so no actual harm; needs fixing if revived.
- **Q7** Conflicting comments on the Peng fit range: the 260811Cl block (L2883–2884) says "the source claims a fit for s ≤ 6 Å⁻¹", while the 260923Cl table note (L1222–1223) says "for s ≤ 2 Å⁻¹ (ITC 4.3.2.2)". The latter reads as the newer correction; the former comment was not updated.

## Items to confirm

1. Error of the approximation in `FactorImaginary` that places $\mathbf g$ perpendicular to the beam, for HOLZ reflections and large EBSD $\mathbf g$ (large $g_z$).
2. Applying the symmetric local form ($\mathbf K\mp\mathbf g/2$) to an annular detector (`FactorImaginaryAnnular`, local) makes the integration domain centred on a point shifted by $\mathbf g/2$ rather than on the incident beam. Is this the intended definition for STEM HAADF?
3. X-ray white (`WaveColor.FlatWhite`) intensity ∝ $d^2$ while neutron white ∝ $d^4$ (Crystal.cs L1759, L1770). Possibly "flat per energy" vs "flat per wavelength" (added by Claude, to be confirmed); needs the author's intent. → also listed in [structure-factor.md](structure-factor.md).
4. Any plan to use Temari `Fx` for the X-ray structure factor (currently WK only).
5. `XrayScattering` (Cromer–Mann) and `ElectronScatteringEightGaussian` are unused inside this library. Are they used by upper applications (ReciPro)?
6. The 8-Gaussian table is sourced only from a private communication.
7. Edition of `IsotopeAbundance` (the code itself says "edition to be identified").
8. f′/f″ of the neutral element are used for ions too (chemical shift and charge dependence ignored).
9. The elastic real part (`es.Factor`) extrapolates Peng beyond $s>2$ Å⁻¹ also for high-angle HOLZ and EBSD calculations (left unchanged as outside the practical range per the comment).

Items to confirm: 9
