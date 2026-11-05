This is a draft written by Claude from the code (2026-11-05). Not yet reviewed by the author.

# EBSD master pattern: storage format, projection, file

Target commit: `9ffcf47` (main). Files: `EBSD/EBSD.cs` (`EBSD`, `MasterPattern`), `EBSD/MasterPatternFile.cs`, `EBSD/EbsdPatternProjector.cs`, `EBSD/EbsdPatternComposer.cs` (lookup table part).

For the dynamical calculation itself (the intensity for each direction, energy and thickness by the Bethe method), see the diffraction document (the Bethe method and Kikuchi chapters in `docs/methods/`). This document only describes how its output is arranged, stored, and looked up. For coordinate frames and the sense of rotations, see [ebsd-geometry.md](ebsd-geometry.md).

---

## 1. Construction (`EBSD.RunMasterPatternBuild`)

### 1.1 Input (`EBSD.MasterPatternBuildRequest`)

| Argument | Meaning | Default |
|---|---|---|
| `Crystal` | Crystal | Required |
| `MaxNumOfBloch` | Maximum number of Bloch waves for the Bethe calculation | Required |
| `Energies` | Energy grid [keV] | Required (descending order per the doc on the MC distribution side) |
| `Depths` | Depth grid [nm] (path length along the exit direction, monotonically increasing; §2.3) | Required |
| `GridSize` | Number of divisions per side ($2N+1$ for the hexagonal grid) | Required |
| `Solver` | `BetheMethod.Solver` | `Auto` |
| `Thread` | Number of threads | 1 |
| `UseNonLocalAbsorption` | Non-local absorption | false |
| `IncludeTDSBackground` | TDS background (non-local source) | false |
| `IncludeAbsorbedFluxBackground` | Re-injection of the absorbed flux (260919Cl) | false |

The grid type `GridType` is determined automatically from the crystal system: if `MasterPattern.ShouldUseHexGrid` = `UseHexGridEnabled && CrystalSystemNumber ∈ {5, 6}` (trigonal, hexagonal), the hexagonal grid is used; otherwise the square grid. `UseHexGridEnabled` is a static flag (default true; the comment says "for debugging").

### 1.2 Processing flow

1. `CreateDirectionsSquare` / `CreateDirectionsHex` create the cell-center directions of the $+Z$ hemisphere and the $-Z$ hemisphere, and these are **concatenated in the order $+Z$ → $-Z$** into `directions` (invalid cells of the hexagonal grid are `null`).
2. `BetheMethod.RunEBSDNew(MaxNumOfBloch, Energies, I, Depths, directions, Solver, Thread, UseNonLocalAbsorption, IncludeTDSBackground, IncludeAbsorbedFluxBackground)`. The crystal rotation is the identity matrix. The `BetheMethod` settings are
   - `EbsdRepresentativeDirectionBlockSize` = 1 (GridSize ≤ 512), 3 (≤ 1024), 5 (larger)
   - `UseLocalSurfacePerBeamDirection = false` (comment: in the new path the crystal is rotated, so the local surface approximation is not used)
3. After completion, a separate `BackgroundWorker` bakes `Disks[ei][di].Amplitudes` (complex, length 2 × number of hemisphere cells) into planes:

$$ \texttt{PositivePlanes}[k][i] = |A_i|^2,\qquad \texttt{NegativePlanes}[k][i] = |A_{i+N_{\rm hemi}}|^2,\qquad k = e\cdot n_D + d . $$

   Missing entries are 0. Planes are independent, so `Parallel.For` is used. Each baked `Disks[ei][di]` is set to `null` as it goes to lower peak memory (comment: complex amplitudes total ~2.6 GB at grid 512).
4. Progress is allocated 99 % to the Bethe stage and 1 % to the baking stage (`BetheStageWeight`, `CompileStageWeight`).

> The plane values are the squared magnitudes of `CBED_Disk.Amplitudes`. Whether this quantity is an "amplitude" or already corresponds to the square root of an intensity is a convention on the `BetheMethod` side, so check it in the diffraction document (Item to confirm M-1).

`BetheMethod` maps equivalent cells on the grid to each other using `MasterPattern.GetMasterPatternSquareSymmetryOperations` / `GetMasterPatternHexSymmetryOperations` (§3.3) (around L950–980 of `Diffraction/BetheMethod.cs`).

---

## 2. Storage format of `MasterPattern`

### 2.1 Arrays

- `PositivePlanes[k]`, `NegativePlanes[k]`: `float[GridSize²]`. $k = e\cdot n_D + d$ (**energy-major**). Missing planes have length 0 or are `null`.
- `GetPlane(hemisphere, e, d)`: `null` if out of range.
- `Energies`, `Depths`: the grid values (held by reference, not copied).
- `PlaneCount` $= n_E n_D$.

### 2.2 Hemispheres

`Hemisphere.PositiveZ` / `NegativeZ` are the hemispheres split by **the sign of $z$ in the crystal Cartesian frame** ($z\parallel\mathbf c$ of [ebsd-geometry.md](ebsd-geometry.md) §1.3). For both hemispheres, $(x,y)$ on the plane is mapped by the same formula (§3.1) in the orientation seen from the $+z$ side (the $-Z$ hemisphere is not mirrored).

### 2.3 Meaning of depth

Definition according to the doc of `MasterPattern.DepthIntervals`:

- The plane at `Depths[j]` $=t_j$ is the **cumulative response** up to thickness $t_j$, $M(t_j) = \int_0^{t_j} R(z)\,dz$.
- The mean response of interval $j$ is

$$ \bar R_j = \frac{M(t_j)-M(t_{j-1})}{\Delta t_j},\qquad \Delta t_0 = t_0,\ \ \Delta t_j = t_j - t_{j-1}\ (t_{-1}=0). $$

- `ComputeDepthIntervals`: returns 1 for intervals whose width is ≤ 0 or non-finite (no division). The construction of MC distribution weights uses the same function.
- Doc of `EbsdMonteCarloDistribution.Depths`: the depth grid is the **path length along the exit direction** ("because the master pattern uses the tangent-sphere approximation, and for each direction the equivalent surface normal = −exit direction").

---

## 3. Rosca–Lambert equal-area projection

### 3.1 Square grid

Half side $L=\sqrt{\pi/2}$ (`SquareLimit`), step $\Delta = 2L/N$. Center of cell $(w,h)$ ($w$ = column, $h$ = row, linear index $h N + w$):

$$ a = -L + (w+\tfrac12)\Delta,\qquad b = L - (h+\tfrac12)\Delta \quad(b \text{ decreases as the row goes down}). $$

**Square → sphere** (`RoscaLambertToSphereSquare`; the comment calls it Rosca's square → disk formula):

$$
|b|\le|a|:\ \ (A,B) = \tfrac{2a}{\sqrt\pi}\bigl(\cos\tfrac{\pi b}{4a},\ \sin\tfrac{\pi b}{4a}\bigr),\qquad
|a|<|b|:\ \ (A,B) = \tfrac{2b}{\sqrt\pi}\bigl(\sin\tfrac{\pi a}{4b},\ \cos\tfrac{\pi a}{4b}\bigr),
$$

$$ \rho^2 = A^2+B^2,\qquad (x,y,z) = \Bigl(\sqrt{1-\rho^2/4}\,A,\ \sqrt{1-\rho^2/4}\,B,\ \pm\bigl(1-\rho^2/2\bigr)\Bigr). $$

$(a,b)=(0,0)$ maps to $(0,0,\pm1)$. The corner $(L,L)$ maps to $(A,B)=(1,1)$, $z=0$ (equator). The square's area $2\pi$ = the solid angle of a hemisphere, so the map is equal-area.

**Sphere → square** (`SphereToRoscaLambertSquare`; the comment calls it the Shirley/Rosca inverse): for normalized $(x,y,z)$,

$$ r_s = \sqrt{(1+|z|)/2}\ \ (\text{branch on the sign; the implementation uses } (1+z)/2 \text{ if } z\ge0,\ \text{otherwise } (1-z)/2),\qquad (A,B) = (x,y)/r_s, $$

$$
|A|\ge|B|:\ a = \operatorname{sgn}(A)\,\sqrt{A^2+B^2}\,\tfrac{\sqrt\pi}{2},\ \ b = \tfrac{4a}{\pi}\arctan\tfrac{B}{A};\qquad
|A|<|B|:\ b = \operatorname{sgn}(B)\sqrt{A^2+B^2}\,\tfrac{\sqrt\pi}{2},\ \ a=\tfrac{4b}{\pi}\arctan\tfrac{A}{B}.
$$

Comment (260325Cl): `atan` is used, not `atan2` (`atan2` returns $\theta\pm\pi$ for $a<0$, which is wrong). The projection side always passes $|z|$, so only the $z\ge0$ branch is used.

**Interpolation** (`InterpolatePlaneSquare`): bilinear with $w_f = (a+L)/\Delta - 0.5$, $h_f = (L-b)/\Delta - 0.5$. Indices are clamped to $[0,N-1]$ (values of the opposite hemisphere beyond the equator are not used). The square-grid path of `EbsdPatternComposer.BuildLookupTable` differs only in clamping $(w_0,f_w)=(0,0)$ if $w_0<0$ and $(N-2,1)$ if $w_0\ge N-1$ (the values come out the same).

### 3.2 Hexagonal grid (trigonal, hexagonal)

- $N = (\texttt{GridSize}-1)/2$. Axial coordinates $(u,v)$; valid cells satisfy $|u|\le N,\ |v|\le N,\ |u+v|\le N$. Linear index $(v+N)(2N+1)+(u+N)$. Directions of invalid cells are `null`; plane values there are 0.
- Step $\Delta_h = P_g/N$, $P_g = 2\sqrt\pi/3^{3/4}$ (`Hex_preg` = circumradius of the hexagon; changed from the apothem to this in 260331Cl).
- Axial → Cartesian: $(x,y) = \bigl((u+v/2)\Delta_h,\ \tfrac{\sqrt3}{2}v\Delta_h\bigr)$.
- Hexagon → disk (`HexagonToDisk`; the comment says it follows EMsoft `Lambert2DHexForward`): swap the coordinates to $(X,Y)=(h_y,h_x)$, determine the sextant (`GetSextant`, following EMsoft), rotate into sextant 0,

$$ q = \frac{Y_0\,\pi/(2\sqrt3)}{X_0},\quad (XX_0, YY_0) = 3^{1/4}\sqrt{2/\pi}\;X_0\,(\cos q,\ \sin q), $$

  rotate back to the original sextant, and undo the swap. The inverse (`DiskToHexagon`) is $x_0 = r_0\,3^{-1/4}\sqrt{\pi/2}$, $y_0 = 3^{-1/4}\sqrt{6/\pi}\,r_0\,\operatorname{atan2}(YY_0,XX_0)$.
- Disk ↔ sphere uses the same formulas as the square grid ($\rho^2$, $\sqrt{1-\rho^2/4}$, $\pm(1-\rho^2/2)$). Sphere → disk is $(A,B)=\sqrt{2/(1+|z|)}\,(x,y)$.
- Interpolation (`InterpolatePlaneHex`, `GetHexBarycentricLookup`): the fractional axial coordinates are rounded to the nearest cell $(r_u,r_v)$ by cube rounding, one of the 6 triangles is chosen from the signs of the offset $(f_u,f_v)$, and the value is interpolated with barycentric coordinates over the 3 neighbors. `GetHexBarycentricLookup` redistributes the weight of an invalid cell half each to the remaining 2 points (`InterpolatePlaneHex` merely treats invalid cells as 0).
- For display and conversion: `RenderHexPlaneToImage` (to a square image while keeping hexagonal coordinates), `ResampleHexToSquarePlane` (resampling onto the square Rosca–Lambert grid; directions are created with the $+Z$ hemisphere formula).

### 3.3 Symmetry operations exact on the grid

Only operations that map indices without interpolation are listed.

- Square grid (4/mmm family, `SymmOper`, 16 kinds): identity, 90°/180°/270° about $z$, left–right and up–down flips, 2 diagonal mirrors, hemisphere inversion $m_z$, inversion, $\bar4$-type compositions, 2-fold rotations about the $x$/$y$/diagonal axes (accompanied by hemisphere inversion). `GetMasterPatternSquareSymmetryOperations` returns a subset for each point group number. Monoclinic (point groups 2, m, 2/m) determines the main axis with `MainAxis`, the 2mm family branches on `PointGroupHMStr`, and $\bar42m$/$\bar4m2$ on `StrSE2p`. For trigonal, hexagonal and cubic crystals, only the operations that map exactly on the square grid (the 3-fold axes of cubic crystals etc. are not included).
- Hexagonal grid (6/mmm family, `HexSymmOper`, 24 kinds): rotations about $z$ in 60° steps, 6 mirrors (`MirrorA0`–`A5`), compositions of these with hemisphere inversion, 6 $C_2'$. Branches for point groups 16–27 (32/3m/$\bar3$m on `StrSE2p`, $\bar6$m2/$\bar6$2m also on `StrSE2p`).

`TransformMasterPatternSquareIndex` / `TransformMasterPatternHexIndex` map the index over the two concatenated hemispheres ($\text{sphere}\cdot N_{\rm hemi} + \text{local}$).

---

## 4. Projection onto the detector

### 4.1 `EbsdPatternProjector` (for matching, single slice)

- Constructor / `Rebuild(geometry)`: for each pixel of the raster $W_r\times H_r$, caches the **line of sight (viewing direction)** (sample frame) `PixelToSampleDirection` at the detector pixel coordinates $\mathrm{col} = (c+0.5)/W_r\cdot W - 0.5$, $\mathrm{row} = (r+0.5)/H_r\cdot H - 0.5$. Independent of orientation. `Rebuild` reuses the buffers (because the geometry changes at every evaluation during geometry calibration; 260920Cl).
- `Project(mp, R, posPlane, negPlane, output)`: at each pixel $\mathbf d^{\rm c} = R^{-1}\mathbf d^{\rm s}$; `posPlane` if $d_z^{\rm c}\ge0$, otherwise `negPlane`. Square grid: `SphereToRoscaLambertSquare(d_x, d_y, |d_z|)` + `InterpolatePlaneSquare`; hexagonal grid: `SphereToRoscaLambertHex` + `GetHexBarycentricLookup`. If one of the planes is missing, that hemisphere is 0.
- `parallel=false` makes the row loop sequential (for callers that parallelize per orientation).
- Note in the doc: single slice only (speed first). MC-weighted composition is reduced to one plane with `EbsdMonteCarloDistribution.ComposeGlobalWeightedPattern` before being passed in ([ebsd-pattern-composition.md](ebsd-pattern-composition.md) §5).

### 4.2 In-plane rotation decomposition (for dictionary brute force, square grid only)

For an orientation $R = r_0 R_z(\phi)$, the crystal-frame line of sight is $\mathbf d^{\rm c} = R_z(-\phi)\,(r_0^{-1}\mathbf d^{\rm s})$. The polar coordinates of $\mathbf u = r_0^{-1}\mathbf d^{\rm s}$ on the Lambert disk do not depend on $\phi$, so once per sphere point (`PrepareSpherePoint`)

$$ q_0 = \tfrac{2}{\pi}\operatorname{atan2}(u_y,u_x),\quad r_a = \frac{\sqrt{u_x^2+u_y^2}/|\mathbf u|}{\sqrt{(1+|u_z|/|\mathbf u|)/2}}\cdot\frac{\sqrt\pi}{2},\quad r_b = \tfrac{4}{\pi} r_a,\quad \text{neg} = [u_z<0] $$

is computed, and for each $\phi$ (`ProjectInPlane`) $g = q_0 - 2\phi/\pi$, $s=\lfloor g+\tfrac12\rfloor$, local angle $t = (g-s)\pi/2\in[-\pi/4,\pi/4)$,

| $(s+8)\bmod 4$ | $a$ | $b$ |
|---|---|---|
| 0 | $r_a$ | $r_b t$ |
| 1 | $-r_b t$ | $r_a$ |
| 2 | $-r_a$ | $-r_b t$ |
| 3 | $r_b t$ | $-r_a$ |

after which only bilinear interpolation is needed. Doc: mathematically equivalent to the inverse transform of `SphereToRoscaLambertSquare` ($\arctan(\tan t) = t$). Numerically the difference is of order ULP ("equivalence group — verified by dict regression"). The comment notes that a 2D rotation in the $(a,b)$ plane after Lambert is not possible (because the Rosca map does not commute with rotations).

### 4.3 `EbsdPatternComposer.BuildLookupTable` (for display)

For the whole display raster, line of sight → Rosca–Lambert → interpolation coefficients are built only once (independent of energy and depth). For the square grid, "index of the top-left cell + $(f_w, f_h)$", and for the hexagonal grid, "3 indices + 3 weights" are stored per pixel, with the hemisphere flag held separately. The line-of-sight formula is in [ebsd-geometry.md](ebsd-geometry.md) §3.6. The composition model is in [ebsd-pattern-composition.md](ebsd-pattern-composition.md).

---

## 5. File format (`MasterPatternFile`, extension `.rpmp`, FormatVersion = 1)

Purpose (doc): building at grid 512 takes several minutes, and the (energy × depth) grid changes slightly every time due to the MC random numbers, so the build is not redone when search or display is repeated under the same conditions.

### 5.1 Layout (little-endian)

| Offset | Type | Content |
|---|---|---|
| 0 | 4 byte | Magic `"RPMP"` |
| 4 | int32 | Version (1) |
| 8 | int32 | Header length $n_h$ (0–$2^{24}$) |
| 12 | $n_h$ byte | Header = UTF-8 JSON (`MasterPatternFileInfo`) |
| Following | zlib | Body |

Body (`BinaryWriter`): `GridSize` (int32), `GridType` (int32), $n_E$ (int32), $n_D$ (int32), `Energies` (double × $n_E$), `Depths` (double × $n_D$), followed by **all planes of the $+Z$ hemisphere → all planes of the $-Z$ hemisphere** in the order $k=0..n_En_D-1$ as "length (int32, 0 = missing) + float32 sequence", and finally the number of electrons (int32) and the electrons (61 byte each: float32 × 15 + 1 flag byte).

Layout of one electron: Depth, Vec.X/Y/Z, Position.X/Y, Energy, TotalEnergyLoss, LastInelasticDepth, LastInelasticEnergyBeforeLoss, LastInelasticEnergyAfterLoss, LastInelasticDirection.X/Y/Z, LastDecoherenceDepth, flags (bit0 = HasLastInelasticEvent, bit1 = HasLastDecoherenceEvent).

### 5.2 Header (`MasterPatternFileInfo`)

`FormatVersion`, `CreatedUtc` (ISO 8601), `Creator`, `CrystalName`, `CrystalFingerprint`, `BeamEnergyKeV`, `MaxNumOfBloch`, `UseNonLocalAbsorption`, `IncludeTDSBackground`, `IncludeAbsorbedFluxBackground`, `Note`, `BackscatteredElectronCount`, `McSampleTiltRad` (default NaN), `McEnergyThresholdKeV` (default NaN), `McIncidentCount`. NaN and ±∞ are written as strings with `JsonNumberHandling.AllowNamedFloatingPointLiterals` (261028Cl; before that, saving with NaN threw an exception).

### 5.3 Crystal matching

`CrystalFingerprint`: `cell=a,b,c(nm, F6),α,β,γ(deg, F4);sg=<SymmetrySeriesNumber>;` followed, for each atom, by `Z:x,y,z(F5),Occ(F4),B_iso,eff(nm², F6)`, plus $\beta_{11},\beta_{22},\beta_{33},\beta_{12},\beta_{23},\beta_{31}$ (E6) for anisotropic atoms, then `;`. Names are not included. `CheckCompatibility` requires an exact fingerprint match and $|E_{\rm file}-E|\le10^{-6}$ keV, and returns a string with the reason if they do not match.

### 5.4 Read/write details

- Saving writes to `path.tmp` and then does `File.Move(overwrite)` (a failure midway does not corrupt an existing file). Compression defaults to `Fastest` (floats do not compress much).
- `Load(path, loadElectrons)`: with `false`, the electron section is not read (avoids a contiguous array of about 5.45 million electrons ≈ 0.7 GiB for MnAl with 1e7 incident electrons).
- `ReadElectronChunks(path, chunkSize=512)`: returns the electrons in arrays of 512 (69,632 B < LOH threshold 85,000 B). Planes are skipped by reading into a 64 KiB discard buffer.
- Reason for storing electrons (doc): changing the MC distribution (source depth mode, amorphous layer, phosphor weighting, grid) requires re-binning from the electrons; keeping the electrons avoids redoing the MC and gives the same distribution every time. The MC distribution itself is not stored.

---

## 6. Implementation questions

### M-I1 The names of the square-grid diagonal mirrors are swapped (no functional impact)

In `TransformMasterPatternSquareIndex`, setting $x = w - c$, $y = c - h$ ($c=(N-1)/2$, $b$ runs opposite to the row),

- `MirrorDiagXY` $(w,h)\mapsto(N-1-h,\,N-1-w)$ is $(x,y)\mapsto(y,x)$ = **mirror about $y=x$** (the comment says "$y=-x$").
- `MirrorDiagX_Y` $(w,h)\mapsto(h,w)$ is $(x,y)\mapsto(-y,-x)$ = **mirror about $y=-x$** (the comment says "$y=x$").

`GetMasterPatternSquareSymmetryOperations` always includes both as a pair for every point group, so the result does not change. The maps of `RotDiagXY180` (2-fold about $[110]$) and `RotDiagX_Y180` (2-fold about $[1\bar10]$) are as in the comments.

---

## 7. Items to confirm

- **M-1** Whether the physical meaning of the plane value $|A|^2$ (what `CBED_Disk.Amplitudes` is the amplitude of) and the interpretation as the cumulative $M(t)$ (doc of `DepthIntervals`) agree with the output definition on the Bethe side. To be cross-checked with the diffraction document.
- **M-2** Square-grid interpolation only clamps indices and does not use values of the opposite hemisphere near the equator. When the grid is coarse, a discontinuity of about half a cell can arise near the equator; is this accepted?
- **M-3** When the hexagonal grid is used for trigonal and hexagonal crystals, is the correspondence between the $x$ axis of the crystal Cartesian frame ($\parallel\mathbf a^*$, [ebsd-geometry.md](ebsd-geometry.md) §1.3) and the orientation of the hexagon's sextants (the coordinate swap in `HexagonToDisk`) consistent with the definition of the mirror axes A0–A5 of `HexSymmOper` (A0 = 0°, A1 = 30°, …)? (The choice of `MirrorA*` and `RotA*_180` is explained per point group in comments, but the basis for the reference angles of the axes is not written.)
- **M-4** Cubic crystals (point groups 28–32) use only operations that map exactly on the square grid. Is it intended that equivalence through the 3-fold axes is not exploited (the reduction in computation would be limited)?
- **M-5** With `MasterPatternFile.FormatVersion` left at 1, the fingerprint format changed in 261105Cl (addition of $\beta$ for anisotropic atoms). Isotropic crystals give the same string as before and are therefore compatible, but old files of anisotropic crystals will no longer load (stated in the doc). Is it right not to bump the version?
- **M-6** The `MasterPattern` constructor holds `Energies`/`Depths`/the plane arrays without copying them. `MasterPatternBuildRequest` copies them, but how is the case handled where the arrays are modified externally, e.g. via the file-loading path?
