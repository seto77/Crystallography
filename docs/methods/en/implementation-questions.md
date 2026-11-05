Draft written by Claude from the code (2026-11-05). Not yet reviewed by the author.

# Implementation questions (diffraction and dynamical theory)

Points where the implementation may be wrong, or where comments disagree with the code, collected while writing the method documents. **No code was changed.**
The full reasoning (equations, counterexamples) is in the "Implementation questions" section of each document. Line numbers refer to `main` at 9ffcf47 (2026-11-05).

Severity: **A** = may affect results / **B** = affects results only under some conditions (latent) / **C** = comment or doc mismatch only.

## A: may affect results

| # | Location | Issue | Details |
|---|---|---|---|
| A1 | `Diffraction/BetheMethod.cs:3865–3866, 2877, 2895, 2282–2283` | The U cache `uDictionary` is keyed only on `(compose(g), compose(h) or int.MaxValue)`, not on kV or detector angles (inner/outer). After clearing at line 2877, STEM writes U′ restricted to the ADF detector [inner, outer] under the same keys as the ordinary U, for every `q+g_i−g_j` (qList holds all differences $g_i-g_j$, so 000 is included). It does not clear afterwards; the only Clears in StemDoWork are at 2523, 2535 and 2877. The real part does not depend on inner/outer, so only the imaginary (absorptive) part is contaminated. `GetDifractedBeamAmpriltudes` clears only when `AccVoltage != voltage` (line 2282). 261014Cl added unconditional Clears to CBED, PED and STEM, but this one stayed conditional. Called after STEM on the same instance and voltage, with N or orientation changed so that it recomputes, it builds a **mixed absorption matrix**: elements whose indices STEM touched get the detector-restricted U′, the others the full-solid-angle U′. HRTEM (`GetHRTEMImage` → `GetDifractedBeamAmpriltudes`) is affected the same way. There is one shared instance per crystal, `Crystal.Bethe` (`Crystal/Crystal.cs:167, 477`), so this can happen if ReciPro runs STEM and the diffraction simulator on the same crystal (not checked on the ReciPro side). EBSD clears before and after muBack (`:1645, :1657`) and uses a separate master-pattern instance (`EBSD/EBSD.cs:233`), so it is not affected; neither are Kikuchi profiles (own instance, `KikuchiProfileCalculator.cs:64`). | [bloch-wave](bloch-wave.md) Q1, [kikuchi-dynamical](kikuchi-dynamical.md) §4-9 |
| A1b | `BetheMethod.cs:2292, 2528, 2337–2338` | **Separate issue found while investigating A1 (see "A1b in detail" below).** STEM overwrites the shared field `Beams` re-sorted by |g|². If `GetDifractedBeamAmpriltudes` then reuses the old eigenvectors, diffracted-beam amplitudes are assigned to the wrong reflections. | below |
| A2 | `BetheMethod.cs:2741–2742` vs `:3735, :3763` | The temporal-coherence envelope differs between STEM, $\exp[-\pi^2\lambda^2\delta^2(k_q^2-k^2)^2]$, and HRTEM (TCC and quasi-coherent), $\exp[-\tfrac12\pi^2\lambda^2\delta^2(\cdot)^2]$. The same `delta` acts √2 times larger in STEM. If δ is meant to be defined differently in the two paths, that should be stated. | [bloch-wave-imaging](bloch-wave-imaging.md) |
| A3 | `BetheMethod.cs:1600–1608, 1630`, `EBSD/EBSD.cs:235, 1191` | For hexagonal lattices with master-pattern grid width G>512 (block width 3 or 5), a block centre can fall on an invalid cell (direction `null`); the whole block then shares the beam set chosen for [00−1]. | [kikuchi-dynamical](kikuchi-dynamical.md) §4-1 |
| A4 | `BetheMethod.cs:2120–2168` `ComputeAbsorbedFluxBackground` | $N(z)=\sum_g|\psi_g|^2$ sums densities. After the 260829Cl row-scaling fix, the conserved quantity without absorption is $\sum_g P_g|\psi_g|^2$, so $D(t)\neq0$ is possible even with $U'=0$ (visible for tilted surfaces / HOLZ). | [kikuchi-dynamical](kikuchi-dynamical.md) §4-2 |
| A5 | `Diffraction/Alchemi/AlchemiReduction.cs:146–171, 194` | Same kind of mismatch: $L_\mathrm{coh}$ counts $P_g/P_0$-weighted current while $Y_\mathrm{dyn}$ counts density. In the limit $\mu=\mu_{00}I$, Total $\neq\mu_{00}t/V_c$. Differs only for HOLZ beams and tilted surfaces. | [alchemi-ionization](alchemi-ionization.md) |
| A6 | `AlchemiReduction.Yield :196` | $Y_\mathrm{dyn}$ may be missing a factor $2K/P_0\ (\approx1/\cos\theta_0)$ (derivation by Claude, needs checking). About 4.5e-4 at 30 mrad. | [alchemi-ionization](alchemi-ionization.md) |
| A7 | `BetheMethod.cs:2292, 2370` | Diffracted amplitudes and PED decide whether to reuse the eigendecomposition from N, voltage and orientation only; changes to atoms, DW factors, Surface or model settings are not detected. | [bloch-wave](bloch-wave.md) Q5 |
| A8 | `Atom/Atoms.cs:355` → `Crystal/Crystal.cs:1644` | The kinematical electron F uses the Peng factors of the assigned ion (SubNumberElectron) and ignores `BetheMethod.ElasticIonModel` (default Neutral). For crystals with assigned ions, kinematical and dynamical calculations use different $f_e$ models (state it if intended). | [structure-factor](structure-factor.md) S1 |

### A1b in detail (additional investigation, 2026-11-05)

**Why the reuse condition passes**
- Only `GetDifractedBeamAmpriltudes` (lines 2310–2318) writes the eigenvalues/eigenvectors (`EigenValues/EigenVectors/EigenVectorsInverse`). STEM, CBED, PED and ALCHEMI neither write nor null them.
- `RunSTEM` sets `MaxNumOfBloch/AccVoltage/BaseRotation` to the caller's values (lines 2496–2501). The worker then sets `Beams = [.. Find_gVectors(BaseRotation, vecK0).OrderBy(e => e.Vec.Length2)]` (line 2528).
- `Matrix3D.Equals` compares the nine components by value (`Mathematics/Matrix.cs:650`), so `rotation.Equals(BaseRotation)` is true for the same orientation.
- `Find_gVectors` returns the same set for the same inputs. Candidates are collected in a sequential loop (lines 4418–4425), and `QuickSelect` uses the last element as pivot, which is deterministic (`NthItemUtils.cs`). STEM computes vecK0 from a U0 read after a Clear, which gives the same value at the same voltage. Hence `Beams.Length == EigenValues.Length` also holds.
- So in `GetDifractedBeamAmpriltudes(N,V,R,t₁)` → `RunSTEM(N,V,R,…)` → `GetDifractedBeamAmpriltudes(N,V,R,t₂)`, the third call satisfies every condition at line 2292 and reuses the old eigendecomposition without recomputing.

**What happens**
- The eigendecomposition was built in `Find_gVectors` order (ascending Rating = |g|·Q², line 4472), whereas `Beams` is now in STEM's |g|² order.
- In `alpha = V⁻¹ e₀`, e₀ is 000 in both orders (both Rating and |g|² are 0), so the incident condition is correct.
- Lines 2337–2338 (`Beams[i].Psi = psi_atZ[i]`) assign the amplitude of the i-th reflection in the old order to the i-th reflection in the new order. The exit-surface phase `exp(iπ P_i t)` (line 2331) uses the new order's P. Since P is real, |ψ|² is **permuted** among reflections, and the phases are also wrong.
- Only the transmitted beam (000) is correct.

**Realistic sequences**
- HRTEM(N,V,R) → STEM(N,V,R) → HRTEM with only defocus or thickness changed. `GetHRTEMImage` calls `GetDifractedBeamAmpriltudes(BlochNum, AccVol, rot, t)` for each thickness (line 3723), so the reuse starts at its first call and the whole image is wrong.
- Diffracted amplitudes(N,V,R) → STEM(N,V,R) → amplitudes again with only the thickness changed.
- This assumes that the instance is shared per crystal as `Crystal.Bethe` (as in A1). How ReciPro calls these is not checked.

**How much the order differs (simplified calculation by Claude, approximate)**
Only the Rating formula was reproduced in Python (Si, 200 kV, near [001]; U0, tangent-plane duplicate removal and the BFS cache limit are omitted). The first N+1 beams in Rating order were compared with the same beams re-sorted by |g|².

| Tilt | N | Beams that change position | of which swapped with a reflection of different |g| |
|---|---|---|---|
| 0° | 50 | 0/51 | 0 |
| 0° | 200 | 200/201 | 196 |
| 0.5° | 50 | 45/51 | 33 |
| 2° | 200 | 200/201 | 198 |

- Exactly on a zone axis with only ZOLZ beams, Q = −g², so Rating = |g|⁵ and the order matches |g|² order. Swaps occur only within groups of equal |g|. If those are symmetry-equivalent, their intensities are equal and nothing is visible, except for accidentally equal non-equivalent pairs such as Si 333/511.
- When N is large enough to include HOLZ-ring reflections (Q≈0 with large |g|), or with any tilt, almost every beam is swapped with a different reflection.

**Related points**
- CBED (line 457) uses the same `Find_gVectors` call, so its order is unchanged. PED (line 2460) re-sorts its union by Rating, so even when its count equals N+1 the order mostly matches.
- The STEM worker re-reads the `Beams` field 23 times after line 2536 (e.g. lines 2647, 2759, 2765). If `GetDifractedBeamAmpriltudes` recomputes on the same instance while STEM is running, `Beams` is replaced mid-run. The potential matrix (built at line 2536) and `qList` (line 2759) then refer to different beam sets, and with a different N this can index out of range. `IsSTEM_Busy` (line 317) exists, but `GetDifractedBeamAmpriltudes` does not check it.
- Possible fixes (no code changed): keep STEM's beams in a local variable instead of the `Beams` field, or keep a reference to the beam array used for the eigendecomposition and add `ReferenceEquals(Beams, eigenBeams)` to the condition at line 2292.

## B: conditional / latent

| # | Location | Issue |
|---|---|---|
| B1 | `BetheMethod.cs:3900–3918` | With anisotropic DW, the effective value along g, $m=h'^\mathsf{T}\beta h'/s_g^2$, is used as an isotropic B inside the TDS integral, whose integrand involves $K\pm g/2$ in other directions — an approximation ([scattering-factors](scattering-factors.md) Q4). |
| B2 | `Diffraction/Kikuchi/KikuchiProfileCalculator.cs:94–113` | `ComputeCrystalHash` omits `Dsf.UseIso`; toggling only UseIso leaves a stale snapshot considered fresh. |
| B3 | `KikuchiProfileCalculator.cs:261–267` | The two-beam path takes $I_\mathrm{ref}$ from the +g side only; if that side returns 0 from a guard, c jumps (unreachable under current invariants). |
| B4 | `KikuchiProfileCalculator.cs:362` | Comment says "x=+1 is the Bragg line of +g", but Q=0 for member +g is at x≈−1. Correct if named after the outgoing beam; check what ReciPro's drawing assumes. |
| B5 | σₙ / Kikuchi kernel / `getU` | Atoms without B are handled inconsistently (σₙ floors at 1e-3 nm², the Kikuchi kernel substitutes 5e-3 nm², `getU` has no floor), giving atoms that backscatter but do not absorb. |
| B6 | `BetheMethod.cs:2100, 2258`, `KikuchiProfileCalculator.cs:733` | $F_{jj'}$ loses precision when $\lambda t$ is very small (the branch tests λ itself). Harmless for normal absorption levels. |
| B7 | `Diffraction/IonizationChannel.cs` `IonizationTableShape` | `_truncated/_usedTail` are mutable and shared; reusing IonizationData across runs can carry over flags from a previous run. |
| B8 | `IonizationChannelTable.GridAt` | For E0 beyond the outermost table rows, PCHIP may extrapolate silently (only 30–400 keV is range-checked). |
| B9 | `BetheMethod.cs:3450–3452, 3481` | Setting the environment variable `RECIPRO_STEM_EDX_HERMTOL=NaN` bypasses the Hermitian check. |
| B10 | `Atom/Atoms.cs:366–371` | Neutron isotope weighting depends on Dictionary enumeration order (not guaranteed by the spec). |
| B11 | `Atom/AtomStatic.cs:3174` `FactorImaginaryAnnularFlatEwald` | Missing the radial Jacobian R, and the prefactor $\gamma/(\pi^2k_0)$ does not match the other variants. **Not used in production.** |
| B13 | `Diffraction/Alchemi/BetheMethod.Alchemi.cs:73, 258` | ALCHEMI reads U0 with `getU(AccVoltage)` (line 73), but clears only at line 258, just before building the potential matrix. If a previous run left values at another voltage, the real part of U0 (which includes γ) is from that voltage, which shifts $K=\sqrt{k^2+U_0}$, Q and P. The U that `Find_gVectors` stores in each Beam inside `BuildFixedUnion` also comes from the cache, but the potential matrix is rebuilt after the Clear, so the only effect is through U0. Of the entry points, only ALCHEMI lacks the "Clear before reading U0" that 261014Cl added to CBED, PED and STEM. |
| B12 | Phase sign | `Crystal.GetStructureFactor` uses $e^{-2\pi i h\cdot r}$, `BetheMethod.getU` uses $e^{+2\pi i g\cdot r}$ (`BetheMethod.cs:3921`, i.e. the textbook $U_{-g}$). Each is internally consistent as checked in the documents, but any place comparing kinematical and dynamical results needs care ([structure-factor](structure-factor.md) S2, [bloch-wave](bloch-wave.md)). |

## C: comment/doc mismatches (numbers are correct)

| # | Location | Issue |
|---|---|---|
| C1 | `Atom/AtomStatic.cs:1209–1213` | `MottBetheMonopoleCoefficient` doc says "s² [Å⁻²], f [nm]", but 0.2393366 is the coefficient for s² [nm⁻²], f [nm] (0.023934 Å in Å units). Value checked by Claude. |
| C2 | `AtomStatic.cs` `FactorImaginary*` returns, `BetheMethod.cs:3918` | Documented as "dimensionless"; actually in nm, like $f_e$. |
| C3 | `BetheMethod.cs:3734–3735` | In the HRTEM quasi-coherent branch the Ec/Es (temporal/spatial) labels are swapped: the β term is spatial, the δ term temporal. Checked in the code by Claude. |
| C4 | `Diffraction/Alchemi/AlchemiAngularSpread.cs:48–50` | "Weight <1e-7 beyond 4σ (≈1.7 FWHM)"; the actual weight is $e^{-8}\approx3.3\times10^{-4}$. Checked by Claude. |
| C5 | `AtomStatic.cs:2883` vs `:1222` | The Peng fitting range is given as "s≤6" (260811Cl) and "for s≤2" (260923Cl). |
| C6 | `BetheMethod.cs:815–825` | `GetCrystalFixedMasterEquivalentBeamAndSurface` always returns (d, −d) up to rounding — an identity. Results are correct but the comment describes something else. |
| C7 | `KikuchiProfileCalculator.cs:7, 592`, `AlchemiReduction.cs:12` | The old `/P_col` form and the old notation $S=B^\dagger\mathrm{diag}(\sigma)B$ remain after the 260829Cl fix. |
| C8 | `IonizationTableShape`, `AlchemiBasisDiagnostic`, `BuildFixedUnion:186–189`, `IonizationFsTable.SCount` | Numbers from older dataset versions remain ("s>4", "s≤8 Å⁻¹", "81 points"; current: 321 points, 16 Å⁻¹). |
| C9 | `AlchemiMuBuilder.cs:206` | Comment says "same rule as BisoEffective", code uses `Biso000` (identical at G=0, no effect). |

## Other observations (not errors)

- `Atom/NistElastic.cs` is used only by the EBSD MonteCarlo; nothing on the diffraction path references it.
- Of the Temari factors, only `Fe` is used (inside the TDS-absorption integrand); `Fx` is unused and X-ray f remains WK.
- The Cromer–Mann and 8-Gaussian tables are not used inside the library.
- `Diffraction/DynamicTheory.cs` (multislice) is `[Obsolete]` and fully commented out.
- The native `_EBSDSolver` / `_EBSDSolverWithTDS` / `STEM_InelasticQ` sources are not in this repository, so their equations were not checked.
