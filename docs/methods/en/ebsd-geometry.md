This is a draft written by Claude from the code (2026-11-05). Not yet reviewed by the author.

# EBSD coordinate frames, detector geometry, and orientation representation

Target commit: `9ffcf47` (main). Files: `EBSD/EbsdDetectorGeometry.cs`, `EBSD/EbsdPatternComposer.cs` (`UpdateTiltCoefficients`, `BuildLookupTable`, `ExitRay`), `EBSD/MonteCarlo.cs` (definition of the sample surface), `EBSD/EbsdIndexer.cs` (rotation utilities), `EBSD/EbsdDictionaryIndexer.cs` (`GetProperRotations`), `Crystal/Crystal.cs` (crystal Cartesian frame).

All EBSD-related classes share the conventions in this document. The doc of `EbsdDetectorGeometry` says "the conventions are made to match the display forward model (`EbsdPatternComposer`) exactly", and states that the agreement of the two sets of formulas was checked with the `EbsdIndexCheck` harness (the harness is not in this repository).

---

## 1. Coordinate frames

### 1.1 Laboratory frame (lab)

The Monte Carlo (`MonteCarlo`) and the detector geometry use the same lab frame.

- Origin = the beam impact point (on the sample surface).
- The initial direction of the incident electron is $\mathbf{v}_0 = (0,0,-1)$ (pointing along $-z$).
- The sample is tilted about the $x$ axis by the tilt angle $\tau$ (`MonteCarlo.Tilt`, `EbsdDetectorGeometry.SampleTilt`, rad). A point $\mathbf{p}$ is inside the sample when

$$ p_y \tan\tau \ge p_z . $$

- Hence the outward normal of the sample surface (lab) is $\hat{\mathbf n}_s^{\rm lab} = (0,\,-\sin\tau,\,\cos\tau)$, and the normal depth below the surface is

$$ d = p_y\sin\tau - p_z\cos\tau \quad(\text{in MC, updated as } d \mathrel{+}= s\,(v_y\sin\tau - v_z\cos\tau)). $$

  The vacuum side is the $-y$, $+z$ side for $\tau>0$ and the $+y$, $+z$ side for $\tau<0$. **The UI (`FormEBSD`) sample tilt is in $[-90^\circ,0^\circ]$, default $-70^\circ$**, so in practice it is the latter (detector on the $+y$ side, default `DetY` = −35 mm).

### 1.2 Sample frame (sample)

lab → sample frame is the rotation $R_x(\tau)$ about the $x$ axis (`EbsdDetectorGeometry.LabToSample`; according to the comment, identical to the $M\cdot v$ action of OpenTK `CreateRotationX`):

$$
s_x = l_x,\qquad s_y = \cos\tau\, l_y + \sin\tau\, l_z,\qquad s_z = -\sin\tau\, l_y + \cos\tau\, l_z .
$$

The inverse is `SampleToLab`. This rotation maps $\hat{\mathbf n}_s^{\rm lab}\mapsto(0,0,1)$, so **$+z$ of the sample frame is the outward normal of the sample surface**. The incident beam in the sample frame is $(0,-\sin\tau,-\cos\tau)$.

The same rotation is also used in the ctor of `EbsdMonteCarloDistribution` and in `ComputePathLengthUpperBound` (which brings the electron exit direction `Vec` [lab] back to the sample frame). The doc states explicitly that the tilt passed to the MC distribution is "the value when the MC was run", not the current UI value.

### 1.3 Crystal Cartesian frame (crystal Cartesian)

`Crystal.cs` (around L684) places the lattice vectors as follows:

$$
\mathbf c = (0,0,c),\quad
\mathbf b = (0,\, b\sin\alpha,\, b\cos\alpha),\quad
\mathbf a = \Bigl(a\sqrt{1-\cos^2\beta-\tfrac{(\cos\gamma-\cos\alpha\cos\beta)^2}{\sin^2\alpha}},\ a\,\tfrac{\cos\gamma-\cos\alpha\cos\beta}{\sin\alpha},\ a\cos\beta\Bigr).
$$

That is, **$\mathbf c\parallel z$, $\mathbf b$ lies in the $yz$ plane, $x\parallel \mathbf b\times\mathbf c\ (\parallel\mathbf a^*)$**. Reciprocal lattice vectors, zone axes, and the directions of Kikuchi reflections (`VectorOfG_KikuchiLine`) are all handled as components in this Cartesian frame.

> Difference from other software: the choice of axes differs from conventions commonly used in EBSD (e.g. $x\parallel\mathbf a$, $z\parallel\mathbf c^*$) (added by Claude; to be confirmed). For hexagonal crystals, $\mathbf b\parallel y$, and $\mathbf a$ lies in the $xy$ plane at $-30^\circ$ from $x$ ($x\parallel\mathbf a^*$). When comparing with orientation data from other software, this difference in the choice of crystal frame must be included in the rotation.

---

## 2. Orientation representation

### 2.1 Rotation matrix

The EBSD code handles orientations **only as rotation matrices**. `EbsdOrientationCandidate.Rotation`, `EbsdZoneAxisSolution.Rotation`, and `EbsdMatchingContext.Rotation` are all "the rotation crystal → sample frame (same convention as `Crystal.RotationMatrix`)":

$$ \mathbf v^{\rm s} = R\,\mathbf v^{\rm c}. $$

- The columns of $R$ are the basis vectors of the crystal Cartesian frame expressed in the sample frame.
- Read passively, it is the coordinate transformation "crystal-frame components → sample-frame components"; read actively, it is "the rotation that turns the sample-frame axes onto the crystal axes".
- Projection uses the inverse: a sample-frame direction $\mathbf d^{\rm s}$ is brought back to the crystal frame $\mathbf d^{\rm c}=R^{-1}\mathbf d^{\rm s}$ before looking up the master pattern (both `EbsdPatternProjector.Project` and `EbsdPatternComposer.BuildLookupTable` use `rotation.Inverse()`).
- The 9-argument constructor of `Matrix3D` is column-major (comment in `EbsdIndexer.SolveWahba`).

### 2.2 Euler angles

The code under `EBSD/` does not use Euler angles (both input and output are rotation matrices). The repository has a Z-X-Z `FromMatrix` in `Mathematics/Euler.cs`, but it is not called from the EBSD code. The correspondence with Bunge angles $(\varphi_1,\Phi,\varphi_2)$ is not defined in the EBSD code (Item to confirm G-1).

For reference (added by Claude; to be confirmed): under the convention in which the Bunge orientation matrix $g$ is "sample frame → crystal frame" ($\mathbf v^{\rm c} = g\,\mathbf v^{\rm s}$), $R = g^{\mathsf T}$. However, because the choice of crystal Cartesian frame (§1.3) and of sample frame (§1.2) differ from other software, $R$ cannot be compared directly with another program's $g$.

### 2.3 Symmetry equivalence and misorientation

- The proper rotations of the point group (`EbsdDictionaryIndexer.GetProperRotations`) are obtained by transferring the linear part $W$ of the symmetry operations (in lattice coordinates, those with $\det W=+1$) to the Cartesian frame as $S = M\,W\,M^{-1}$ ($M$ = `MatrixReal`, columns $\mathbf a,\mathbf b,\mathbf c$), with duplicates and the identity removed. If there is no symmetry, as in P1, the result is `null`.
- Equivalent orientations are $R\,S$ (multiplied from the right).
- Misorientation angle (`EbsdIndexer.MisorientationDeg`):

$$ \theta(R_1,R_2) = \min_{S\in\{I\}\cup\mathcal S}\ \arccos\!\Bigl(\tfrac{\operatorname{tr}(R_1^{\mathsf T}R_2S)-1}{2}\Bigr). $$

  If `properSymmetries` is not passed, only $S=I$.

### 2.4 Orientation perturbation

The perturbation used in orientation refinement (Nelder–Mead) is a **left perturbation in the sample frame** (`EbsdIndexer.PerturbRotation`):

$$ R' = \mathrm{Rot}\bigl(\hat{\boldsymbol\omega},|\boldsymbol\omega|\bigr)\,R_0,\qquad \boldsymbol\omega=(\omega_x,\omega_y,\omega_z)\ [\deg]. $$

### 2.5 Orientation grid

- `EbsdIndexer.FibonacciSphereRotation(di, nSphere)`: the minimal rotation that turns $+z$ onto the Fibonacci sphere lattice point $\hat{\mathbf n}_i$ ($z_i = 1-2(i+\tfrac12)/N$, azimuth $i\cdot\pi(3-\sqrt5)$).
- The orientations for dictionary search and Radon coarse search are $R = R_{\rm sphere}(i)\,R_z(\phi_j)$. Since $R\,\hat{\mathbf z}^{\rm c} = \hat{\mathbf n}_i$, the sphere point is **the direction of the crystal $z$ axis ($\mathbf c$) in the sample frame**, and $\phi_j$ is the rotation about it. This is not a strictly uniform grid on SO(3) (Item to confirm G-2).

---

## 3. Detector geometry (`EbsdDetectorGeometry`)

### 3.1 Parameters

| Symbol | Property | Unit | Meaning |
|---|---|---|---|
| $\delta$ | `DetTilt` | rad | Detector tilt. Default in the doc 90° ($\pi/2$) |
| $(X_d,Y_d,Z_d)$ | `DetX`, `DetY`, `DetZ` | mm | Detector center parameters (values of the UI numeric fields) |
| $s$ | `PixelSize` | mm/px | Pixel size |
| $W,H$ | `WidthPx`, `HeightPx` | px | Number of pixels |
| $x_m$ | `XMirror` | ±1 | Left–right flip (same as `FormEBSD.DetectorXMirror`) |
| $\tau$ | `SampleTilt` | rad | Sample tilt |

The ctor rejects with an exception non-finite values, $s\le0$, $W,H\le0$, $x_m\ne\pm1$, and configurations in which the detector plane passes through the origin ($|\mathbf n\cdot\mathbf C|<10^{-12}$).

### 3.2 Derived quantities (lab frame)

$$
\mathbf C = (X_d,\,-Y_d,\,-Z_d),\quad
\mathbf e_x=(1,0,0),\quad
\mathbf e_y=(0,-\cos\delta,-\sin\delta),\quad
\mathbf n=(0,\sin\delta,-\cos\delta).
$$

$\{\mathbf e_x,\mathbf e_y,\mathbf n\}$ is orthonormal. The camera length (perpendicular distance from the impact point to the detector plane) is

$$ L = |\mathbf n\cdot\mathbf C| \quad(\texttt{CameraLength} = \texttt{FormEBSD.CameraLength2}). $$

### 3.3 Pixel → direction

A pixel $(\mathrm{col},\mathrm{row})$ follows the convention **corner origin, pixel center at $(\mathrm{col}+0.5,\mathrm{row}+0.5)$**. The display mm coordinates relative to the image center are

$$ u = (\mathrm{col}+0.5-W/2)\,s,\qquad v=(\mathrm{row}+0.5-H/2)\,s \quad(\texttt{PixelToMm}). $$

Point on the detector plane (lab):

$$ \mathbf P = \mathbf C + x_m\,u\,\mathbf e_x + v\,\mathbf e_y \quad(\texttt{PixelToLabPoint}). $$

For $\delta=90^\circ$, $\mathbf e_y=(0,0,-1)$, so increasing row moves toward lab $-z$ (down).

- **Line of sight (viewing direction)** (`PixelToSampleDirection`): $\hat{\mathbf d}^{\rm s} = R_x(\tau)\,(-\mathbf P/|\mathbf P|)$. The direction from the point on the detector toward the impact point.
- **Exit direction** (the direction in which the electron leaves the sample): $\hat{\mathbf k}_{\rm out}^{\rm s} = -\hat{\mathbf d}^{\rm s} = R_x(\tau)\,\hat{\mathbf P}$. Zone-axis indexing (`EbsdZoneAxisIndexer`) and lookups into the MC distribution (`EbsdPatternComposer.ExitRay`, `ComputeDetectorAverageSliceWeights`) use this one, flipping the sign explicitly.
- **The master pattern is looked up along the line of sight $\hat{\mathbf d}$ ($=-\hat{\mathbf k}_{\rm out}$)**: `EbsdPatternProjector` and `EbsdPatternComposer.BuildLookupTable` select the hemisphere by the sign of the crystal-frame $z$ component of $\hat{\mathbf d}^{\rm c}=R^{-1}\hat{\mathbf d}^{\rm s}$ and map to Rosca–Lambert coordinates ([ebsd-master-pattern.md](ebsd-master-pattern.md)). Which direction the master pattern is computed for (incident or exit direction) is a convention of the dynamical calculation; see the diffraction document (Item to confirm G-3).

Inverse map `SampleDirectionToPixel(k_out)`: $\mathbf v = R_x(\tau)^{-1}\mathbf k_{\rm out}$, $t = (\mathbf n\cdot\mathbf C)/(\mathbf n\cdot\mathbf v)$; `null` if $t\le0$ or $|\mathbf n\cdot\mathbf v|<10^{-12}$. From $\mathbf q=t\mathbf v-\mathbf C$: $u = x_m\,\mathbf q\cdot\mathbf e_x$, $v=\mathbf q\cdot\mathbf e_y$, $\mathrm{col}=u/s+W/2-0.5$, $\mathrm{row}=v/s+H/2-0.5$.

### 3.4 Pattern center (PC)

The PC is **the foot of the perpendicular dropped from the impact point onto the detector plane**, $\mathbf F=(\mathbf n\cdot\mathbf C)\,\mathbf n$. The physical in-plane coordinates relative to the detector center (`PatternCenterMm`) are

$$ u_F = -X_d,\qquad v_F = -(Y_d\cos\delta + Z_d\sin\delta). $$

In image pixel coordinates (`GetPatternCenter`), $x_m$ is applied to $x$ on the display side:

$$ \mathrm{col}_{\rm PC} = x_m u_F/s + W/2 - 0.5,\qquad \mathrm{row}_{\rm PC} = v_F/s + H/2 - 0.5,\qquad DD = L . $$

That is, the PC convention of this code is "**(PC pixel coordinates [corner origin, pixel center], DD [mm])**"; values made dimensionless by the detector width or DD (the PC conventions of Bruker, EDAX/TSL, Oxford, EMsoft, etc.) are not implemented (Item to confirm G-4).

The inverse `FromPatternCenter(u_F, v_F, DD, δ)` solves the following simultaneous equations (as written in the code):

$$ Y_d\cos\delta + Z_d\sin\delta = -v_F,\qquad Y_d\sin\delta - Z_d\cos\delta = -DD,\qquad X_d=-u_F, $$

$$ \Rightarrow\ Y_d = -(v_F\cos\delta + DD\sin\delta),\quad Z_d = -v_F\sin\delta + DD\cos\delta. $$

The second equation means $\mathbf n\cdot\mathbf C = +DD$. The $L$ in the doc's "standard configuration signed L = −DD" is `FormEBSD.SignedCameraLength` $=Y_d\sin\delta-Z_d\cos\delta=-\mathbf n\cdot\mathbf C$, so they agree (G-I1, G-I2 in §5).

### 3.5 Band center line ↔ plane normal

An image line $A\,\mathrm{col}+B\,\mathrm{row}+C=0$ (pixel-center convention) is converted to a line $a u + b v + c = 0$ in display mm coordinates:

$$ a = A/s,\quad b=B/s,\quad c = A c_x + B c_y + C,\qquad c_x = W/2-0.5,\ c_y=H/2-0.5 . $$

Normal of the plane containing the band center line and the origin (lab, normalized, sign undetermined; `LineToLabNormal`):

$$ \mathbf g \propto a\,\mathbf e_u + b\,\mathbf e_y + \gamma\,\mathbf n,\qquad \mathbf e_u = x_m\mathbf e_x,\qquad \gamma = \frac{c - a\,(\mathbf e_u\cdot\mathbf C) - b\,(\mathbf e_y\cdot\mathbf C)}{\mathbf n\cdot\mathbf C}. $$

(Chosen so that $\mathbf g\cdot\mathbf P = 0$ is equivalent to $au+bv+c=0$.) The inverse is `LabNormalToLine`: from $a_{\rm mm}=x_m g_x$, $b_{\rm mm}=\mathbf g\cdot\mathbf e_y$, $c_{\rm mm}=\mathbf g\cdot\mathbf C$, it gives $(a_{\rm mm}s,\ b_{\rm mm}s,\ c_{\rm mm}-a_{\rm mm}s\,c_x-b_{\rm mm}s\,c_y)$.

Bragg angle from a band edge point $\mathbf P_e$ (`SinBraggFromEdgePoint`): $\sin\theta_B = |\hat{\mathbf g}\cdot\mathbf P_e|/|\mathbf P_e|$.

### 3.6 Coefficient form for the display raster (`EbsdPatternComposer`)

Display composition does not use `EbsdDetectorGeometry`; it holds the same formulas as coefficients (`UpdateTiltCoefficients`):

$$
y_c = \cos\tau\cos\delta + \sin\tau\sin\delta,\quad z_c = -\sin\tau\cos\delta + \cos\tau\sin\delta,\quad
y_0 = \cos\tau\,Y_d + \sin\tau\,Z_d,\quad z_0 = -\sin\tau\,Y_d + \cos\tau\,Z_d .
$$

The display coordinates of a raster pixel $(w,h)$ are, in `EbsdRasterView`,

$$ p_x = (2w+1-\mathrm{width})\,S_W + O_X,\qquad p_y = (2h+1-\mathrm{height})\,S_H + O_Y $$

($O_X,O_Y$ are the display pan offsets; $S_W$ corresponds to half the pixel pitch). The line of sight (sample frame, unnormalized) is

$$ \mathbf d^{\rm s} = \bigl(-x_m p_x - X_d,\ \ y_c\,p_y + y_0,\ \ z_c\,p_y + z_0\bigr), $$

which is rotated into the crystal frame by $R^{-1}$. It agrees with $-R_x(\tau)\mathbf P$ of §3.3 (with $u=p_x$, $v=p_y$). The exit direction is `ExitRay`: $(x_m p_x + X_d,\ -(y_c p_y + y_0),\ -(z_c p_y+z_0))$. `ExitRay.PlaneDistance` $=|z_0 y_c - y_0 z_c|$ equals $L$.

---

## 4. Unit conventions

| Quantity | Unit | Notes |
|---|---|---|
| Detector coordinates, PC, DD | mm | |
| Angles (public API) | rad (geometry), deg (tolerances, step sizes, perturbations) | Arguments whose names contain `Deg` are in degrees |
| Wavelength (indexing) | nm | Doc of `EbsdIndexer.Index`: the default 0.00859 nm is 20 kV. Passing Å saturates the width likelihood for all bands |
| Wavenumber (MC) | nm⁻¹, $k=1/\lambda$ (without $2\pi$) | `UniversalConstants.Convert.EnergyToElectronWaveNumber` |
| Depth, path length | nm | |
| Energy | keV (MC, master pattern), eV (TPP-2M, losses) | |

---

## 5. Implementation questions

### G-I1 `FromPatternCenter` returns only the configuration with $\mathbf n\cdot\mathbf C=+DD$ (**checked: no problem within the UI range**)

- Formula (§3.4): the round trip `GetPatternCenter` → `FromPatternCenter` is the identity only when $\mathbf n\cdot\mathbf C>0$. A configuration with $\mathbf n\cdot\mathbf C<0$ (e.g. $\delta=90^\circ$, $(X_d,Y_d,Z_d)=(1,20,5)$) is mirrored to $(1,-20,5)$ by the round trip (checked by evaluating the formulas in Python).
- Conventions in `FormEBSD` (seto77/ReciPro, checked at `2add19b`, 2026-11-05):
  - Sample tilt `numericBoxSampleTilt`: range $[-90^\circ, 0^\circ]$, **default $-70^\circ$**. Passed with the same sign to MC (`new MonteCarlo(..., sampleTilt, ...)`, `sampleTilt = SmpTilt`) and to `BuildDetectorGeometry` (`SampleTilt = SmpTilt`).
  - Detector: default `DetTilt` = 90°, `DetY` = −35 mm, `DetZ` = 30 mm (`numericBoxYofDet` / `numericBoxZofDet`).
  - `FormEBSD.SignedCameraLength` $= Y_d\sin\delta - Z_d\cos\delta\ (=-\mathbf n\cdot\mathbf C)$. Comment: "in the standard configuration (DetTilt=90°, DetY<0) this is −CameraLength2 — the same convention as 'standard configuration signed L = −DD' in `EbsdDetectorGeometry.FromPatternCenter`".
- For the default configuration, $\mathbf n\cdot\mathbf C = +35 = +DD$, the round trip is the identity $(0,-35,30)$, and the sample-frame $z$ component of the exit direction at the detector center is $+22.6$ ($\tau=-70^\circ$, above the sample surface). For $\tau<0$ the vacuum side of §1.1 is the $+y$ side, and the detector sits on the $+y$ side ($Y_d<0$).
- Conclusion: within the UI operating range ($\tau\le0$, detector on the exit side) $\mathbf n\cdot\mathbf C>0$ always holds, so geometries built through `FromPatternCenter` in `EbsdGeometryCalibrator` and `EbsdOrientationSearch` are correct. The remaining caveat is only an `EbsdDetectorGeometry` constructed outside the UI with $\tau>0$ (detector at $Y_d>0$). Suggest stating "assumes $\mathbf n\cdot\mathbf C>0$" in `FromPatternCenter`, or passing $\operatorname{sgn}(\mathbf n\cdot\mathbf C)$ of the original geometry as an argument (suggestion; code unchanged).

### G-I2 The definition of "signed L" in the doc of `EbsdDetectorGeometry` is not in the same file (**checked**)

"signed L" is `FormEBSD.SignedCameraLength` $=Y_d\sin\delta - Z_d\cos\delta = -\mathbf n\cdot\mathbf C$. The second equation of `FromPatternCenter`, $Y_d\sin\delta - Z_d\cos\delta = -DD$, is exactly this, so the code and the doc agree. Reading only the Crystallography side one may misread it as $L=\mathbf n\cdot\mathbf C$; adding a one-line definition to the doc would help.

---

## 6. Items to confirm

- **G-1** Whether to formally define the correspondence with Bunge Euler angles (the relation between $R$ and $g$, including the differences in the choice of axes of the crystal Cartesian frame and the sample frame). At present `EBSD/` has no conversion.
- **G-2** The orientation grid $R_{\rm sphere}(i)R_z(\phi_j)$ (Fibonacci × equally spaced $\phi$) is not uniform on SO(3) (the $\phi$ direction is dense near the poles of $\hat{\mathbf n}$). Is this non-uniformity taken into account in evaluating the coverage of the coarse search?
- **G-3** The correspondence between looking up the master pattern along the line of sight $\hat{\mathbf d}=-\hat{\mathbf k}_{\rm out}$ and the meaning of the directions on the Bethe calculation side (`directions` passed to `BetheMethod.RunEBSDNew`) (incident or exit direction, application of the reciprocity theorem). To be checked in the diffraction document.
- **G-4** Conversion formulas to the PC conventions of other vendors and software (Bruker PCx/PCy/DD, EDAX/TSL $x^*,y^*,z^*$, Oxford, EMsoft $x_{pc},y_{pc},L$) are not implemented (added by Claude; to be confirmed).
- **G-5** The crystal Cartesian frame ($\mathbf c\parallel z$, $x\parallel\mathbf a^*$) differs from the choice commonly used in other EBSD software (added by Claude; to be confirmed). If data exchange is anticipated, document the conversion.
- **G-6** The physical meaning of `XMirror` (whether the left–right flip of the detector image is a flip on the camera side or on the display side). The code only specifies "multiply the $x$ component of the display mm coordinate $u$ by $x_m$".
