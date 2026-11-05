This is a draft written by Claude from the code (2026-11-05). Not yet reviewed by the author.

# Structure factor and kinematic intensities

Scope: `Crystal/Crystal.cs` (three versions of `GetStructureFactor`, `SetVectorOfG`, `SetVectorOfPlane`, `SetVectorOfG_KikuchiLine`, `SetPeakIntensity`, `GetDspacingList`).
Definitions and units of atomic scattering factors and temperature factors: [scattering-factors.md](scattering-factors.md). The dynamical potential $U_g$: [bloch-wave.md](bloch-wave.md).

---

## 1. Symbols and units

| Quantity | Definition in the code |
|---|---|
| $\mathbf g=h\mathbf a^*+k\mathbf b^*+l\mathbf c^*$ | `A_Star` etc. (nm⁻¹), $|\mathbf g|=1/d$ |
| $s^2$ | `vec.Length2/4.0` or `1/d²/4.0` = $(\sin\theta/\lambda)^2$ [**nm⁻²**] |
| Atomic coordinates | fractional `atom.X, atom.Y, atom.Z` (`Atoms.Atom[]` already expanded to equivalent positions) |
| Scattering-factor units | electrons: nm; X-rays: electrons×0.1; neutrons: fm (all including Occ) |
| Unit of $F$ | as above (not divided by the volume) |

## 2. Structure factor `Crystal.GetStructureFactor`

Three versions:

| Signature | Purpose |
|---|---|
| `public static Complex GetStructureFactor(WaveSource, Atoms[], (double h,k,l), double s2, double xrayEnergyKeV=NaN, bool anomalousDispersion=true)` (L1545) | non-integer indices, unoptimised |
| `private static Complex GetStructureFactor(WaveSource, Atoms[], (int h,k,l), double s2, Complex[] neutronFactors=null, double xrayEnergyKeV=NaN, (double fp,double fpp)[] xrayDispFactors=null, bool anomalousDispersion=true)` (L1620) | main implementation (real-amplitude fast path, 260605Cl) |
| `public Complex GetStructureFactor(WaveSource, (int h,k,l), double xrayEnergyKeV=NaN, bool anomalousDispersion=true)` (L1712) | instance version; computes $s^2=|\mathbf g|^2/4$ internally |

### 2.1 Formula (as implemented)

$$
F(\mathbf h)=\sum_{\text{Atoms}}\ \sum_{j\in\text{equiv. positions}} f_{\text{Atoms}}(s)\;T_j(\mathbf h)\;\exp\!\big(-2\pi i\,\mathbf h\cdot\mathbf r_j\big)
$$

- **Phase convention $e^{-2\pi i\mathbf h\cdot\mathbf r}$**: `AddStructureFactorContribution` (L1511) does `realSum += A.Re cosφ + A.Im sinφ`, `imagSum += A.Im cosφ − A.Re sinφ` ($\varphi=2\pi\mathbf h\cdot\mathbf r$), i.e. adds $A e^{-i\varphi}$.
- $f$ includes Occ ([scattering-factors.md §2](scattering-factors.md)). $f$ is reused for consecutive sites with the same element, sub-number **and Occ** (Occ condition added at 260715Ch; before, the first site's value was reused even when Occ differed).
- Temperature factor:
  - isotropic (`Dsf.UseIso`): $T=\exp(-B_{iso}s^2)$, $B$ [nm²]; NaN → 1.
  - anisotropic: $T_j=\exp[-(\mathbf h'^{\mathsf T}\beta\,\mathbf h')]$ with $\mathbf h'=$ `atom.Operation.ConvertPlaneIndex(h,k,l)` $=R^{\mathsf T}\mathbf h$ (β dimensionless, off-diagonal terms doubled); NaN → 1.
  - no special case for $\mathbf h=0$ when anisotropic (gives $T=1$, harmless).
- Neutrons: $f$ independent of $s$; `SetVectorOfG` passes `neutronFactors` precomputed once for all reflections.
- Electrons: real-amplitude path (`AddStructureFactorContributionReal`). **No γ**, and **`ElasticIonModel` is not consulted** (finite Peng part of the assigned ion).

### 2.2 X-ray anomalous dispersion (260606Cl)

$$
f=f_0+0.1\,\mathrm{Occ}\,\big(f'-i f''\big),\qquad f_0=\texttt{WK.Factor}(s^2)\cdot\mathrm{Occ}
$$

(`XrayDispersionFactor`, L1529; 0.1 converts electrons to the internal unit; a NaN f′/f″ term is set to 0.)

- The minus sign on $f''$ comes from the $e^{-2\pi i\mathbf h\cdot\mathbf r}$ convention (comment: "matches the + convention of International Tables, verified"). That is, $F_{\text{code}}(\mathbf h)=\overline{F_{\text{ITC}}(\mathbf h)}$ (ITC: $f=f_0+f'+if''$, $e^{+2\pi i}$), so $|F|^2$ and the sign relation of Bijvoet differences agree (added by Claude, to be confirmed: the conjugation relation is my derivation from the code).
- Activation `IsXrayDispersionActive`: `anomalousDispersion && wave==Xray && E>0 finite && Xraylib.Enabled`. The default `anomalousDispersion = true`, but `xrayEnergyKeV` defaults to NaN, so **calls without an energy get no dispersion**.
- `SetVectorOfG` and `SetPeakIntensity` precompute f′/f″ once per site (`xrayDispFactors`).
- With dispersion $F(-\mathbf h)\ne\overline{F(\mathbf h)}$, so `SetVectorOfG` computes both signs independently (the conjugate fill was removed at 260606Cl).

## 3. Reciprocal-lattice enumeration and kinematic intensities

### 3.1 `SetVectorOfG(dMin, dMax, wavesource, maxNum=25000, xrayEnergyKeV, anomalousDispersion)` (L1186)
- Grows shells from the origin, collecting up to `maxNum` points with $|\mathbf g|<1/d_{min}$, sorted by $|\mathbf g|$; removes those below $1/d_{max}$.
- Stores ± of each point in `VectorOfG`. Reflections with a non-null `Symmetry.GetFirstExtinctionRule(h,k,l)` get $F=0$ (not computed).
- `RawIntensity = |F|²`, `RelativeIntensity = RawIntensity / max` (no normalisation if all zero). **No Lorentz, polarisation or volume correction.**

### 3.2 `SetVectorOfPlane`, `SetVectorOfG_KikuchiLine`
Both give $|F|^2$ normalised to the maximum. The Kikuchi version sets $F=0$ for forbidden reflections and sorts by `RelativeIntensity` ascending. Dynamical Kikuchi band profiles: [kikuchi-dynamical.md](kikuchi-dynamical.md).

### 3.3 Powder intensities `SetPeakIntensity(waveSource, waveColor, λ, …)` (L1720)
$\sin\theta=\lambda/(2d)$, $2\theta=2\arcsin(\sin\theta)$ (**radians**). For each hkl in `Plane[i].strHKL` split at `+`, $F$ is computed, multiplied by the multiplicity `Multi[j]`, and summed. $V_c^2$ = `CellVolumeSquare`.

| Source | Colour | `eachIntensity` |
|---|---|---|
| X-ray | Monochrome | $|F|^2 m\,\dfrac{1+\cos^2 2\theta}{\sin 2\theta\,\sin\theta}\Big/V_c^2$ |
| X-ray | FlatWhite | $|F|^2 m\,d^2/V_c^2$ |
| Electron | — | $|F|^2 m\,\dfrac{1}{\sin2\theta\,\sin\theta}\Big/V_c^2$ |
| Neutron | Monochrome | $|F|^2 m\,\dfrac{1}{\sin2\theta\,\sin\theta}\Big/V_c^2$ |
| Neutron | white | $|F|^2 m\,d^4/V_c^2$ |

Finally normalised by the maximum over all peaks (`Intensity`).
- The monochromatic X-ray LP factor $\dfrac{1+\cos^22\theta}{\sin2\theta\sin\theta}=\dfrac{1+\cos^22\theta}{2\sin^2\theta\cos\theta}$ is 1/2 of the textbook form (added by Claude, to be confirmed: $\frac{1+\cos^22\theta}{\sin^2\theta\cos\theta}$). No effect after normalisation. No monochromator polarisation, absorption or extinction corrections.
- Friedel pairs overlapping at the same 2θ in a powder are approximated by $|F|^2$ of the representative hkl (no multiplicity averaging of Bijvoet asymmetry; 260606Cl comment).

### 3.4 `GetDspacingList(λ, maxNum=512, bestNum=16, maxTwoThetaDeg=150)` (L1343)
Returns the `bestNum` strongest d values by X-ray intensity with LP, $|F|^2 m\,(1+\cos^22\theta)/(\sin2\theta\sin\theta)$.
- 261105Cl: upper limit on $2\theta$ (**degrees**, default 150°), because the LP factor diverges as 2θ→180°.
- 261105Cl: equivalent planes counted with the Laue group (`IsRootPlane(..., addInversion: true)`), multiplicity also from the Laue group. d values within a relative $10^{-6}$ are merged and their intensities added.
- No anomalous dispersion or energy (calls `GetStructureFactor(WaveSource.Xray, …)` without dispersion).

## 4. Differences in convention from the dynamical calculation

| Location | Phase | Scattering factor | γ | Normalisation |
|---|---|---|---|---|
| `Crystal.GetStructureFactor` | $e^{-2\pi i\mathbf h\cdot\mathbf r}$ | Peng of the assigned ion (electrons) / WK (+f′,f″) (X-rays) | no | none ($F$ itself) |
| `BetheMethod.getU` | $e^{+2\pi i\mathbf g\cdot\mathbf r}$ (changed 20240524) | Peng per `ElasticIonModel` | $\gamma$ | $\gamma/(\pi V_c)$ → nm⁻² |
| `BetheMethod.ComputeIonizationU` | $e^{-2\pi i\mathbf G\cdot\mathbf r}$ (260801Cl "generate the negative phase directly") | inner-shell ionisation form factor | no | $k_{vac}/(2\pi V_c)$ |

Where the positive phase of `getU` is conjugated to be consistent overall is described in [bloch-wave.md](bloch-wave.md) (the comment of `ComputeIonizationU` calls it "the empirical convention of the existing getU: build the positive phase and conjugate"). For ALCHEMI/EDX see [alchemi-ionization.md](alchemi-ionization.md).

---

## Implementation questions

- **S1** Electron kinematic $F$ and dynamical $U_g$ treat ions differently ([scattering-factors.md](scattering-factors.md) Q3). `Crystal.GetStructureFactor` uses `SubNumberElectron` as is and does not follow `BetheMethod.ElasticIonModel` (default Neutral). For crystals with ions assigned, kinematic spot intensities and the dynamical $U_g$ use different models.
- **S2** The phase convention is opposite between `Crystal.GetStructureFactor` (−) and `BetheMethod.getU` (+). The `getU` comment gives an empirical reason ("centrosymmetric crystals should now work"); whether the overall convention is right for non-centrosymmetric crystals with absorption (complex $U'_g$) is for the Bloch-wave side ([bloch-wave.md](bloch-wave.md)) to confirm. Only the facts in the code are stated here.

## Items to confirm

1. The factor 1/2 in the LP factor (cancels after normalisation, but matters if absolute intensities are ever needed).
2. Intent of X-ray white ∝ $d^2$ vs neutron white ∝ $d^4$ ("flat per energy" vs "flat per wavelength"?) (added by Claude, to be confirmed). Absence of a polarisation factor for X-ray white (constant at fixed 2θ).
3. Source of the Lorentz factor $1/(\sin2\theta\sin\theta)$ used for electron powder intensities (none cited in the code).
4. `GetDspacingList` ignores anomalous dispersion; behaviour for wavelengths near absorption edges.
5. Whether falling back to $T=1$ when Occ/DW are unset (NaN) is the intended behaviour for $F$ in `SetVectorOfG` etc.
6. Check the relation $F_{\text{code}}=\overline{F_{\text{ITC}}}$ under anomalous dispersion (orientation of Bijvoet pairs) against the ReciPro display (which one is called $+\mathbf h$).

Items to confirm: 6
