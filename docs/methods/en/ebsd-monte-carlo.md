This is a draft written by Claude from the code (2026-11-05). Not yet reviewed by the author.

# EBSD Monte Carlo: electron transport and the (energy × depth) distribution

Target commit: `9ffcf47` (main). Files: `EBSD/MonteCarlo.cs`, `EBSD/EbsdBackscatterSimulator.cs`, `EBSD/EbsdMonteCarloDistribution.cs`. For coordinate frames, see [ebsd-geometry.md](ebsd-geometry.md) §1.

Overall roles:

1. `MonteCarlo` tracks one electron in the sample and returns the depth, direction, energy, etc. of an escaped electron.
2. `EbsdBackscatterSimulator.Run` runs this in large numbers in parallel and selects by an energy threshold.
3. `EbsdMonteCarloDistribution` sorts the electrons into bins of the exit hemisphere and fits, per bin, "the exit energy distribution $G(E)$", "the mean normal depth of the source $\lambda_d(E)$", and "the electron fraction $F$". The composer ([ebsd-pattern-composition.md](ebsd-pattern-composition.md)) interpolates these per pixel to build the (energy × depth) weights.

---

## 1. Material parameters (`MonteCarlo` ctor)

```
MonteCarlo(z, a, ρ, kev, tilt, thresholdKev = 2,
           stoppingPowerModel = JablonskiModified2008,
           elasticScatteringModel = MottNistSampler2023,
           inelasticScatteringModel = DiscreteBulkDiimfpApproximation,
           valenceElectronCount = null, bandGapEv = null, atoms = null, meanIonizationPotentialEv = null)
```

| Symbol | Argument | Unit | Convention |
|---|---|---|---|
| $Z$ | `z` | – | **Atom-number average** $\sum n_iZ_i/\sum n_i$ |
| $A$ | `a` | g/mol | **Atom-number average** (mass per atom) |
| $\rho$ | `_ρ` | g/cm³ | |
| $E_0$ | `kev` | keV | Incident energy |
| $\tau$ | `tilt` | rad | Sample tilt (about the $x$ axis) |
| $E_{\rm th}$ | `thresholdKev` | keV | Tracking cutoff (default 2) |
| $N_v$ | `valenceElectronCount` | – | Number of valence electrons for TPP-2M. If null, the atom-number average from `atoms`; if that is also absent, estimated from $Z$ |
| $E_g$ | `bandGapEv` | eV | Default 0 |
| $J$ | `meanIonizationPotentialEv` | eV | If null, from `atoms` by **Bragg's rule** (average of $\ln J_i$ weighted by electron number); otherwise $J(Z)$ |

Caution in the doc (260921Cl): $Z$, $A$, $N_v$ are used internally as the **ratios** $\rho Z/A$ (Bethe stopping power) and $N_v\rho/A$ (TPP-2M), so they must be averaged with the same weights (number of atoms). Previously $Z$ and $N_v$ were mass-weighted, and the ratios were increasingly wrong for compounds containing hydrogen (for Cu₂(OH)₃Cl, $Z/A$ was off by a factor of 1.93 and $N_v\rho/A$ by 1.63). `GetMeanAtomicParameters(atoms)` weights by Multiplicity × Occ and excludes atoms with $Z\notin[1,99]$.

Elemental $J$ (Joy & Luo 1989; `ElementIonizationPotentialEv`):

$$ J(Z) = \begin{cases} 11.5\,Z & (Z\le12)\\ 9.76\,Z + 58.5\,Z^{-0.19} & (Z>12)\end{cases}\quad[\mathrm{eV}]. $$

Other derived quantities:

- Composition-averaged Debye–Waller $\bar B$ [nm²] = $\sum(\text{Occ}\cdot\text{mult})\,B_{\rm iso,eff}/\sum(\text{Occ}\cdot\text{mult})$, with lower bound `DiffuseScatteringFactor.BisoFloorNm2` (1e-3 nm² according to the comment).
- TPP-2M material constants: $u = N_v\rho/A$, $E_p = 28.8\sqrt u$ [eV], $\beta = -0.10 + 0.944/\sqrt{E_p^2+E_g^2} + 0.069\rho^{0.1}$, $\gamma = 0.191\rho^{-1/2}$, $C = 1.97 - 0.91u$, $D = 53.4 - 20.8u$.
- Plasmon cutoff $q_c = \omega_p/v_F$: $k_F = (3\pi^2 n_v)^{1/3}$, $n_v = N_v\rho N_A/A$ [nm⁻³], $q_c = E_p/((\hbar^2/m_e)k_F)$ ($\hbar^2/m_e = 0.0761996$ eV·nm²), divided by $2\pi$ and stored in the $k=1/\lambda$ convention (`InelasticLocalizationQNm`; the comment gives ≈ 1.9 nm⁻¹ for Si).
- The loss regarded as localized, `InelasticLocalizedLossEv` $=\max\bigl(1.8\max(E_p,\,E_{\min}+0.5),\,30\bigr)$ eV, $E_{\min}=\max(E_g\text{ or }1,\,0.5)$.

---

## 2. Transport parameters

Precomputed in 1 eV steps (0 eV–$\lceil 1000E_0\rceil$ eV, with $E_0$ in keV; the 0 eV element is filled with the 1 eV value) (`BuildTransportParameterCache`), and looked up in the hot loop by the index `Math.Round(E·1000)` (`GetTransportParametersRef`; out-of-range is clamped to the end).

### 2.1 Elastic scattering

**Mott / NIST SRD 64** (`MottNistSampler2023`, default): code-generated PCHIP tables per element (`AtomStatic`). 111 energy points, logarithmically equally spaced from 50 eV to 36.4 keV (50 eV–20 keV from the NIST sampler, above 20 keV appended from the DCS; 260603Cl).

- Mixtures: from the number density $n_i$ [nm⁻³] and cross section $\sigma_i$ of element $i$, the macroscopic cross section $\Sigma = \sum_i n_i\sigma_i$, $\lambda_{\rm el} = 1/\Sigma$, and the effective microscopic cross section $\sigma = \Sigma/\sum n_i$. Linear interpolation between the 2 adjacent points on a logarithmic energy axis.
- Scattering angle: the scattering element is chosen from the CDF of $n_i\sigma_i/\Sigma$, and the PCHIP cumulative distribution of that element is evaluated by the inverse-transform method (nearest energy point).
- For out-of-range energies, or elements without generated data, falls back to screened Rutherford.

**Screened Rutherford** (`ScreenedRutherford`):

$$ \alpha = \frac{3.4\times10^{-3}\,Z^{2/3}}{E\,[\mathrm{keV}]},\qquad
\sigma = \pi\Bigl(\frac{Ze^2}{4\pi\varepsilon_0\,m v^2}\Bigr)^2\frac{1}{\alpha(1+\alpha)},\qquad
\lambda_{\rm el} = \frac{A}{N_A\rho\,\sigma}, $$

$m v^2$ is relativistic (`EnergyToElectronMass` × `EnergyToElectronVelositySquared`). Substituting $mv^2 = E(E+2m_0c^2)/(E+m_0c^2)$ gives agreement with Joy's formula left as a comment in the code, $5.21\times10^{-21}(Z/E)^2\,\tfrac{4\pi}{\alpha(1+\alpha)}\bigl(\tfrac{E+511}{E+1022}\bigr)^2$ [cm²] (checked by Claude). The scattering angle is

$$ \cos\theta = 1 - \frac{2\alpha R}{1+\alpha-R},\qquad R\sim U(0,1). $$

### 2.2 Stopping power

**Joy & Luo (1989)** (`JoyLuo1989`):

$$ \frac{dE}{ds} = \frac{c_3}{mv^2}\ln\Bigl(1.166\,k + \frac{0.583\,mv^2}{e\,J}\Bigr),\qquad k = 0.0299\ln Z + 0.7307, $$

$c_3 = -Z N_A\rho e^4/(4\pi\varepsilon_0^2 A)$ (including unit conversion, keV/nm).

**Jablonski (2008) modified formula** (`JablonskiModified2008`, default):

$$ S\,[\mathrm{eV/Å}] = \frac{D_1\,E^{D_2}\,\ln(D_3E)\,(\rho+D_4)^{D_5}}{\lambda_{\rm in}(E)\,[\mathrm{Å}]},\quad (D_1..D_5) = (7.89271,\ 0.0117088,\ 0.0545126,\ -0.0254488,\ 0.326907), $$

$E$ in eV. $dE/ds = -0.01\,S$ [keV/nm]. If $\lambda_{\rm in}$ is invalid, falls back to Joy–Luo.

### 2.3 Inelastic mean free path (TPP-2M)

$$ \lambda_{\rm in}(E)\,[\mathrm{Å}] = \frac{E}{E_p^2\bigl[\beta\ln(\gamma E) - C/E + D/E^2\bigr]},\qquad \text{NaN if } \gamma E\le1 \text{ or the denominator} \le0. $$

Mean loss per event $\langle\Delta E\rangle = |dE/ds|\,\lambda_{\rm in}$ ("on average $\langle\Delta E\rangle/\lambda_{\rm in}$ = stopping power"). Total scattering rate $\Sigma_{\rm tot} = 1/\lambda_{\rm el} + 1/\lambda_{\rm in}$, elastic probability $p_{\rm el} = (1/\lambda_{\rm el})/\Sigma_{\rm tot}$.

### 2.4 Loss distribution per event

| Model | Loss $\Delta E$ |
|---|---|
| `ContinuousSlowingDownApproximation` | No discrete events. On each flight $E \mathrel{+}= s\,dE/ds$ |
| `DiscreteMeanLoss` | Always $\langle\Delta E\rangle$ |
| `DiscreteBulkDiimfpApproximation` (default) | Mixture PDF of 3 components (low loss $\propto(\omega+0.35E_p)^{-1.6}$, plasmon = Gaussian at $E_p$ [$\sigma = \max(1.5, 0.18E_p)$ eV], high-loss tail $\propto\omega^{-2.2}$ [onset $\min(\omega_{\max},\max(1.8E_p,30))$ eV]) weighted so that the mean is $\langle\Delta E\rangle$. CDFs built at energy points in 10 eV steps; inverse-transform method with a guide table |
| `DiscreteDrudeValenceInnerShell` (260922Cl) | Event rate and mean loss same as above. Valence: DIIMFP from the $q$-integral of the extended Drude loss function $\mathrm{Im}[-1/\varepsilon] = \omega_p^2\gamma\omega/((\omega^2-\omega_q^2)^2+\gamma^2\omega^2)$, $\omega_q = \omega_p + q^2/2$ (atomic units, $\gamma = 4$ eV); inner shell: subshell chosen by the ratio of Bote–Salvat ionization cross sections, $\propto\omega^{-2}$ on $[B, E/2]$. The fraction of inner-shell selections is set so that the mean loss equals $S\lambda_{\rm in}$ (doc: 11.7 % for Si at 20 keV, in agreement with 12.0 % in Table 1 of Vos & Winkelmann 2019). Sampler with 16 energy points (logarithmically equally spaced) |

In all cases $\Delta E\le E$ is enforced by clamping.

---

## 3. Trajectory (`GetBackscatteredElectronDetail`)

Initial state: origin, $\mathbf v=(0,0,-1)$, $E=E_0$, $d=0$.

### 3.1 Discrete inelastic models (default path)

```
while E > E_th:
    s = −ln(max(R, ε)) / Σ_tot(E)               # distance to the next event
    d' = d + s (v_y sinτ − v_z cosτ)
    if d' < 0: break                            # crossed the surface → escape
    d = d'
    if R < p_el:  elastic: sample cosθ, φ = 2πR, rotate the direction
                  if the thermal diffuse test (§4.1) is true, "last decoherence" = d
    else:         inelastic: sample ΔE, "last inelastic" = (d, E, E−ΔE, v)
                  if the localization test (§4.2) is true, "last decoherence" = d
                  E ← max(E − ΔE, 0)
```

The direction is rotated by a closed form that rotates $(\sin\theta\cos\phi,\sin\theta\sin\phi,\cos\theta)$ in a local frame with the travel direction $\mathbf v$ as the polar axis (direct substitution at the singular point $v_z+1<10^{-7}$). $\mathbf v$ is renormalized once every 10 steps.

### 3.2 CSDA path

In each loop $s = -\lambda_{\rm el}\ln R$; before the second and later flights the direction is changed by elastic scattering, $d$ is updated, then $E \mathrel{+}= s\,dE/ds$. Since there are no discrete inelastic events, "last inelastic" is not defined (`HasLastInelasticEvent = false`, NaN).

### 3.3 Return value (`BackscatteredElectronDetail`)

| Field | Meaning |
|---|---|
| `Depth` | $d$ at the time the loop exits = **normal depth of the last event before escape (elastic or inelastic)** [nm]. For electrons that exit by the energy cutoff, the depth at that moment |
| `Direction` | Last flight direction (lab frame, unit vector) |
| `Energy`, `TotalEnergyLoss` | Energy at escape, $E_0 - E$ [keV] |
| `HasLastInelasticEvent`, `LastInelasticDepth`, `…EnergyBeforeLoss`, `…AfterLoss`, `LastInelasticDirection` | Last inelastic event |
| `HasLastDecoherenceEvent`, `LastDecoherenceDepth` | Last decoherence event (§4) |
| `CoherentInelasticCount`, `CoherentAngularVarianceRad2`, `CoherentPathLengthNm` | For investigation (the first 2 are 0 when `CollectCoherenceDiagnostics` = false) |

Escape is decided by `Energy > energyThreshold` (an argument of `EbsdBackscatterSimulator.Run`). `Run` runs `loop` electrons with `Parallel.For` and collects the escaped electrons into `EbsdBackscatteredElectron` (the above + `Position` = Schmidt projection `Stereonet.ConvertVectorToSchmidt(sampleRotation · Direction)`). Progress is added to a shared counter every 1024 electrons and reported roughly every 1 %.

> **Which depth is taken as the "source"** (which of `Depth`, `LastInelasticDepth`, `LastDecoherenceDepth` is passed to the MC distribution) is decided by the caller (`FormEBSD`, not in this repository). The doc of `EbsdMonteCarloDistribution` mentions two modes: "the default mode that takes the source at the 'last inelastic' event" and "the mode that takes it at the last decoherence event (dec)" (Item to confirm C-1).

---

## 4. Decoherence tests (260919Cl)

### 4.1 Thermal diffuse (elastic) — `IsThermalDiffuseElasticEvent`

$$ x = \bar B\,k^2(1-\cos\theta) = 2\bar B s^2,\quad s = k\sin(\theta/2),\ k = 1/\lambda\,[\mathrm{nm}^{-1}];\qquad p_{\rm TDS} = 1-e^{-x}. $$

Probability 1 if $x>20$, $x(1-x/2)$ if $x<0.03$. The coherent component $e^{-x}$ is handled on the Bloch-wave side (master pattern), so in the MC it is treated the same as "nothing happened" (the direction is still changed).

### 4.2 Inelastic — `InelasticDecoherenceProbability`

- 1 if $\Delta E \ge$ `InelasticLocalizedLossEv` (inner shell, single-particle excitation).
- Below that, the probability that the momentum transfer $q$ exceeds the cutoff $q_c$. With the angular distribution being a Lorentzian $dP\propto\theta\,d\theta/(\theta^2+\theta_E^2)$ truncated at the Bethe ridge $\theta_c^2=\Delta E/E$, $\theta_E = \Delta E/2E$, $q^2 = k^2(\theta^2+\theta_E^2)$:

$$ \theta_q^2 = q_c^2/k^2 - \theta_E^2,\qquad
p = \begin{cases} 1 & \theta_q^2\le0\\ 0 & \theta_q^2\ge\theta_c^2\\ 1 - \dfrac{\ln(1+\theta_q^2/\theta_E^2)}{\ln(1+\theta_c^2/\theta_E^2)} & \text{otherwise.}\end{cases} $$

(Collective excitations with $q<q_c$ are treated as preserving the Bloch state. The doc relates this to the fact that Kikuchi bands remain even for plasmon-loss electrons in energy-filtered experiments.)

---

## 5. MC distribution (`EbsdMonteCarloDistribution`)

### 5.1 Input

```
EbsdMonteCarloDistribution((Depth, Vec, Energy)[] bseList, beamEnergy, sampleTilt,
    energies, depths, binCount = 18, amorphousLayerNm = 0,
    energyWeightDeadKeV = NaN, energyFilterMinKeV = NaN,
    absorptionLengthNm = null, energyHistogram = false)
```

- `Depth`: **normal depth** of the source [nm] (which depth it is: see the note in §3.3). `Vec`: exit direction (lab). `Energy`: exit energy [keV].
- `sampleTilt`: the tilt **when the MC was run**. `energies`: the energy grid of the master pattern (descending). `depths`: the depth grid of the master pattern (path length).

### 5.2 Electron weight

$$ w(E) = \phi(E)\,\Theta(E - E_{\min}),\qquad \phi(E) = \begin{cases}\max(0,\,E - E_{\rm dead}) & E_{\rm dead} \text{ valid}\\ 1\end{cases} $$

When $E_{\rm dead}\ge E_0$ or $E_{\min}\ge E_0$, these are disabled (since all electrons would become 0). $\phi$ is the "phosphor response (trial)", and $\Theta$ is an ideal energy filter (260922Cl).

### 5.3 Binning of the exit hemisphere (260921Cl)

Exit direction in the sample frame $\mathbf s = R_x(\tau)\,\mathbf{Vec}$; $s_z\le0$ is discarded (but included in the weight denominator $\sum w$). **Lambert equal-area disk**:

$$ (X,Y) = \sqrt{\frac{2}{1+\hat s_z}}\,(\hat s_x,\hat s_y),\qquad X^2+Y^2\le2, $$

and the square $[-\sqrt2,\sqrt2]^2$ is cut into `binCount` × `binCount` (`DirectionToBinCoords`, bin centers at integer coordinates):

$$ b_x = (X+\sqrt2)\frac{n_b}{2\sqrt2} - \tfrac12,\qquad b_y = (\sqrt2 - Y)\frac{n_b}{2\sqrt2} - \tfrac12 . $$

Electrons are put into the nearest bin $(b_i,b_j)$ by `Round`. The solid angle of a bin inside the disk is $8/n_b^2$ sr (0.0247 sr at the default 18). The side that builds the distribution and the composing side use the same `DirectionToBinCoords` (emphasized in the doc). Reason for dropping the Rosca–Lambert square map (doc): it is only $C^0$ on the diagonals $|x|=|y|$, so background flattening (≈ ∇²) produced X-shaped dark lines.

Bin coverage `BinCoverage` (fraction of area inside the disk, sampled at 64×64 points). Bins with coverage ≥ `MinBinCoverage` = 0.5 use the measured values, and $F$ is divided by the coverage to make it "per full bin". Bins below that are extended from the inside by `ExtendFieldOutward` (weighted average of the 8 neighbors [edges 1, corners 1/2], one layer at a time).

### 5.4 Per-bin fits

- Electron fraction: $F_b = \sum_{b} w / \sum_{\rm all} w$ (including sources within the amorphous layer).
- Amorphous fraction: $f_{A,b} = \sum_{d<a} w / \sum_b w$ ($a$ = `amorphousLayerNm`). Sources within the layer are assigned to "no modulation", and depths in the crystal are counted as $d-a$ (from the bottom of the layer).
- **Energy distribution $G_b(E)$** (normalized to sum 1):
  - Default (`energyHistogram = false`): asymmetric Gaussian. Weighted mean $\bar E$, one-sided weighted variances $\sigma_L,\sigma_R$ (0.5 keV if $n\le1$ or $\sigma<0.01$), $G(E_i)\propto\exp(-(E_i-\bar E)^2/2\sigma_{L/R}^2)$. If a filter is present, multiplied by the fraction of each slice's interval $[E_i-h^-, E_i+h^+]$ (midpoints to the neighbors) that exceeds $E_{\min}$.
  - `energyHistogram = true` (260925Cl): raw histogram. Each electron's weight is distributed linearly to the 2 adjacent grid points (triangular kernel). Outside the grid, to the end point. Doc: the asymmetric Gaussian underestimates the sharp peak on the $E_0$ side by a factor of 2–3.
- **Depth $\lambda_{d,b}(E)$**: for each energy window $[E_i-h, E_i+h)$ ($h$ = half the mean step), the normal depths of the sources in the crystal are averaged (only windows with 4 or more electrons, lower bound 1 nm), and $\lambda(E) = l_a + l_bE + l_cE^2$ is fitted by least squares (a constant if there is 1 valid window, a line if 2). Clamped to $[1,10^6]$ nm on the grid. Interpreted as the mean of an exponential distribution = $\lambda$.
- Fallbacks: bins with fewer than 10 electrons use the $G$ and $\lambda$ summed over all bins (state `Global`); bins with few sources in the crystal or no valid windows use the summed $\lambda$ (`GlobalLambda`); if even the summed one is unavailable, $\lambda = 10^6$ nm (uniform in depth, `NoDepthData`).

### 5.5 Normal depth → path length (`FillPathLengthWeights`, depth mapping A2)

The depth grid of the master pattern is the **path length along the exit direction** $t$, while the MC source depth is the normal depth $d$. $t = d/\mu$ ($\mu=\cos\chi$ = sample-frame $z$ component of the exit direction). If $d$ is exponentially distributed with mean $\lambda_d$, $t$ is also exponential, with decay rate

$$ \alpha = \mu/\lambda_d(E)\quad[\mathrm{nm}^{-1}]. $$

When `absorptionLengthNm` $=\lambda_{\rm abs}(E)$ is passed (for dec mode, 260923Cl), $\alpha \leftarrow \max(\alpha - 1/\lambda_{\rm abs},\ 0.1\alpha)$ (removing the double counting of MC thermal-diffuse survival and the mean absorption of the master pattern).

- `sliceMass = true` (model 2): interval mass $m_j = e^{-\alpha t_{j-1}}(1-e^{-\alpha\Delta t_j})$. For $\alpha\Delta t\le10^{-4}$, $1-e^{-x}$ is evaluated by a 4th-order series.
- `sliceMass = false` (model 0/1): $m_j = e^{-\alpha(t_j-t_0)}\Delta t_j$ (right-endpoint rule of density × interval width; dividing by $t_0$ prevents underflow).
- In both cases, **conditional normalization per energy** is applied before multiplying by the weight $G(E)$ (or $F\,G(E)$): $w_{e,j} = G_e\,m_j/\sum_j m_j$. The tail beyond grid point $t_{n_D}$ is redistributed within the same energy, so the energy distribution is unchanged. If the sum is non-finite, the mass goes to the shallowest slice.

The compatibility arrays `BinWeights` (model 0/1, sum 1 within a bin) and `BinAbsoluteSliceWeights` (model 2, sum $F_b$ within a bin) are approximations built with the bin-center $\mu$ (`LambertBinCenterMu` $=\max(10^{-6}, 1-(X^2+Y^2)/2)$). Display composition rebuilds them with the pixel's $\mu$.

### 5.6 Other

- `GlobalDepthWeights` / `GlobalDepthSliceWeights` (only when there is an amorphous layer): $\sum_b F_b\,(\text{bin weights})$. Used as the reference intensity for sources within the amorphous layer.
- `ComputeElectronWiseCoherence(E0, E_c)` (260927Cl, only when `energyHistogram = true`): the histogram $G_A$ with $A(E)$ applied per electron is stored after dividing by the same normalization constant as $G$. $0\le G_A\le G$ is checked at all nodes; a violation throws an exception.
- `ComputeRangesFromMC`: the 95 % point of the loss and the 99 % point of the depth (QuickSelect).
- `ComputePathLengthUpperBound`: quantile (default 0.999) of the path length $t=(d-a)/\mu$. If a detector is passed, only electrons that hit its image; otherwise $\mu\ge0.2$. If none qualify, the 99 % point of the normal depth / 0.5.
- `ComputeGridFromRanges(E0, ΔE95, T, β = 2)`: energy in 16 steps, step $\max(0.1, \mathrm{round}(\Delta E_{95}/15, 0.1))$ keV, descending from $E_0$, lower bound 1 keV. Depth with 40 points; if $\beta>0$, `BuildGeometricDepthGrid`:

$$ t_i = T\,\frac{e^{\beta i/N}-1}{e^\beta-1}\quad(i=1..N),\qquad T \text{ rounded up to 2 significant digits; each point rounded to 0.001 nm, strictly monotonic}. $$

  $\beta\le0$ gives the conventional equal spacing.

---

## 6. Implementation questions

### C-I1 The CSDA path does not subtract the continuous loss of the escaping flight segment

The order in `GetBackscatteredElectronDetail` (CSDA) is "draw $s$ → scatter → $d' = d + s(\cdots)$ → break if $d'<0$ → $d=d'$ → $E \mathrel{+}= s\,dE/ds$". On the last flight that crosses the surface, the continuous loss for the part inside the sample is not subtracted either. Since this is the loss over one segment ($\sim\lambda_{\rm el}$), the effect should be small, but strictly speaking the escape energy is overestimated (in the discrete models the energy does not change during a segment, so there is no problem).

### C-I2 The meaning of `Depth`, representing "the depth of the last event", may disagree with the doc

The doc of `GetBackscatteredElectrons` only says "the return values are depth (in nm)…", and the doc on the MC distribution side calls `Depth` "the normal depth of the source". In the implementation, `Depth` is, as in §3.3, "the depth of the last event before escape", neither the maximum depth reached nor the depth of the last inelastic event. It needs to be confirmed which one the caller passes (Item to confirm C-1).

### C-I3 Small discrepancies between comments and implementation

- The comment in `ComputeRangesFromMC` (260602Cl) says "80%/99% percentiles", but the implementation uses 95 % / 99 % (changed in 260919Cl).
- The comment on `InelasticLocalizedLossEv` says "same definition as `tailOnsetEv` in `CreateBulkLossSamplerEntry`", but the former is $\max(1.8\max(E_p, E_{\min}+0.5), 30)$ and the latter $\min(\omega_{\max}, \max(1.8E_p, 30))$; they do not agree when $E_p < E_{\min}+0.5$ or when $\omega_{\max}$ is small.

---

## 7. Items to confirm

- **C-1** The choice and default of the source depth passed to the MC distribution (`Depth` / `LastInelasticDepth` / `LastDecoherenceDepth`). To be checked in the implementation of `FormEBSD`.
- **C-2** The exponent of the screening parameter is $Z^{2/3}$ (code), whereas Joy's formula (comment in the code) is often written with $Z^{0.67}$ (added by Claude; to be confirmed). The difference is numerically negligible, but which should be regarded as correct?
- **C-3** `IsThermalDiffuseElasticEvent` is also applied on the screened Rutherford path. Is it acceptable, with either the Mott cross section or the screened Rutherford cross section, to carve out the thermal diffuse fraction as $1-e^{-2Bs^2}$ (assuming the cross section is already the static-atom value)?
- **C-4** The Lorentzian and Bethe-ridge truncation used for the inelastic localization test are an approximation of "the angular distribution after the loss of one event has been determined", and the MC itself does not change direction at inelastic events. The range of applicability of this approximation (low energy, heavy elements).
- **C-5** The fit of mean depth per energy window → quadratic polynomial $\lambda(E)$ assumes that the depth distribution is exponential. Evaluation of the difference from distributions peaked on the shallow side (arising from penetration of the incident beam).
- **C-6** Are `MinBinCoverage = 0.5` and the weights of `ExtendFieldOutward` (edges 1, corners 1/2) empirical choices?
- **C-7** The electron weight $\phi(E)=\max(0,E-E_{\rm dead})$ (phosphor response) is described in the doc as a "trial". Its default value and range of applicability.
