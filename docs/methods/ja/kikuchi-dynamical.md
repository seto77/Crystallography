本稿はコードから Claude が起こした下書き（2026-11-05）。著者未査読。

# 菊池回折の動力学計算 — EBSD マスターパターンと TEM 菊池バンドプロファイル

本稿は Bloch 波（Bethe）法による菊池回折の 2 つの実装を扱う。

| 経路 | 用途 | 入口 |
|---|---|---|
| A. EBSD マスターパターン | 後方散乱電子の出射方向ごとの強度（全球・結晶固定） | `BetheMethod.RunEBSDNew` → `ebsdNew_DoWork`（`Diffraction/BetheMethod.cs`） |
| B. 菊池バンド 1D プロファイル | 透過（前方）菊池バンドの x 方向プロファイル（2 波 / 系統列） | `KikuchiProfileCalculator.ComputeProfile(s)`（`Diffraction/Kikuchi/`） |

- ポテンシャル行列 `getU`・固有値問題 `getEigenMatrix`・ビーム探索 `Find_gVectors` の一般形は [bloch-wave.md](bloch-wave.md) を参照。本稿は菊池・EBSD 固有の部分を書く。
- 散乱因子 $f^{e}$、TDS 吸収形状因子 $f'$（`FactorImaginary*`）、Debye–Waller 因子の扱いは [scattering-factors.md](scattering-factors.md) を参照。
- 単位胞・原子位置の展開や構造因子の規約は [structure-factor.md](structure-factor.md)。
- EBSD のモンテカルロ（深さ・エネルギー分布）、検出器への投影、マスターパターンの格子（Rosca–Lambert）と対称操作の実体は `EBSD/` 側（別文書）。本稿ではそれらとの受け渡し量だけを定義する。

---

## 1. 記号と規約（両経路に共通）

| 記号 | 意味 | 単位・規約 | 実装 |
|---|---|---|---|
| $\mathbf g,\mathbf h$ | 逆格子ベクトル | nm⁻¹、**2π を含まない**（$\lvert\mathbf g\rvert = 1/d$） | `Beam.Vec`。`Crystal.MatrixInverseTransposed` に回転 `BaseRotation` を掛けたもの |
| $s^2$ | 散乱因子の引数 | $s^2 = \lvert\mathbf q\rvert^2/4$ [nm⁻²]（$s=\sin\theta/\lambda$） | `getU` の `s2`、カーネルの `s2` |
| $k_\mathrm{vac}$ | 真空中波数 | nm⁻¹（相対論補正込み、2π なし） | `UniversalConstants.Convert.EnergyToElectronWaveNumber(kV)` |
| $\gamma_\mathrm{rel}$ | 相対論因子 | $1 + eE/(m_0c^2)$ | `getU` 内、`ebsdNew_DoWork` の σ 計算内 |
| $U_\mathbf g = U^{r}_\mathbf g + iU'_\mathbf g$ | 光学ポテンシャルの Fourier 係数 | nm⁻²。$U_\mathbf g = \gamma_\mathrm{rel}F_\mathbf g/(\pi V)$、$V$ [nm³] | `getU(...)` の `(Real, Imag)`。`Imag` が吸収（TDS） |
| $U_0$ | 平均内部ポテンシャル項 | $U^{r}_{000}$ | `u0 = getU(AccVoltage).Real.Real` |
| 位相規約 | 原子位置因子 | $e^{+2\pi i\,\mathbf g\cdot\mathbf r_n}$（`getU`、`CreatePhaseFactors`） | 20240524 に符号を + に変更（`getU` 内コメント） |
| $\hat{\mathbf n}$ | 表面法線（**試料内向き**） | 単位ベクトル。既定 `Surface = (0,0,−1)` | `BetheMethod.Surface` |
| $\mathbf K_0$ | 結晶内の入射波数ベクトル | $\lvert\mathbf K_0\rvert^2 = k_\mathrm{vac}^2 + U_0$、接線成分保存 | `getVecK0` |
| $Q_\mathbf g,\ P_\mathbf g$ | Bethe 方程式の係数 | $Q_\mathbf g = \lvert\mathbf K_0\rvert^2 - \lvert\mathbf K_0+\mathbf g\rvert^2$、$P_\mathbf g = 2\hat{\mathbf n}\cdot(\mathbf K_0+\mathbf g)$ | `getQ`, `getP`, `getQP` |
| $\gamma_j,\ C_{\mathbf g j}$ | 固有値・固有ベクトル | $\gamma$ [nm⁻¹]。波動場 $\psi_\mathbf g(z) = \sum_j C_{\mathbf g j}\alpha_j e^{2\pi i\gamma_j z}$。$\operatorname{Im}\gamma_j>0$ で減衰 | `eigenValues`, `eigenVectors`（column-major: `[j*bLen+g]` = $C_{\mathbf g j}$） |
| $\alpha_j$ | 励起振幅 | $\alpha = C^{-1}\boldsymbol\psi(0)$、$\boldsymbol\psi(0)=\mathbf e_0$（(000) のみ 1） | `PartialPivLuSolve` / `LU().Solve` |
| $z,\ t$ | 深さ | nm。$\hat{\mathbf n}$ 方向の距離 | `Thicknesses` |

固有値問題の行列（[bloch-wave.md](bloch-wave.md) 参照）は 260829Cl 以降、行 $\mathbf g$ を $P_\mathbf g$ で割る形

$$A_{\mathbf g\mathbf h} = \frac{1}{P_\mathbf g}\Bigl[Q_\mathbf g\,\delta_{\mathbf g\mathbf h} + \tilde U_{\mathbf g\mathbf h}\Bigr],\qquad
\tilde U_{\mathbf g\mathbf h} = \begin{cases} iU'_{0} & (\mathbf g=\mathbf h)\\ U_{\mathbf h-\mathbf g}^{\texttt{getU}} & (\mathbf g\neq\mathbf h)\end{cases}$$

である（`getEigenMatrix`。行 = $\mathbf g$、列 = $\mathbf h$、要素 `getU(b[col] − b[row])`。対角の実部 $U_0$ は $\mathbf K_0$ の屈折に畳み込まれていて行列には入らない）。`getU` は $e^{+2\pi i\mathbf q\cdot\mathbf r}$ で和を取るので、$U^{\texttt{getU}}_{\mathbf h-\mathbf g}$ は展開 $V(\mathbf r)=\sum_\mathbf q V_\mathbf q e^{2\pi i\mathbf q\cdot\mathbf r}$ の $V_{\mathbf g-\mathbf h}$ に当たる。本稿の行列（源行列・TDS 行列）もこの並び（行 $\mathbf g$・列 $\mathbf h$ が $V_{\mathbf g-\mathbf h}$）に揃っている。

厚さ積分の共通核（`EBSDSolverManaged`、`ComputeTDSMatrixBackground`、`KikuchiProfileCalculator.Fjj`）:

$$\lambda_{jj'} = 2\pi i\,(\gamma_j - \gamma_{j'}^{*}),\qquad
F_{jj'}(t) = \int_0^t e^{\lambda_{jj'}z}\,dz = \frac{e^{\lambda_{jj'}t}-1}{\lambda_{jj'}}\quad(\lvert\lambda_{jj'}\rvert^2<10^{-30}\ \text{なら}\ F=t).$$

$\lambda_{j'j}=\lambda_{jj'}^{*}$、$F_{j'j}=F_{jj'}^{*}$ が成り立つ。

---

## 2. EBSD マスターパターン（`BetheMethod.RunEBSDNew`）

### 2.1 入力と既定値

`RunEBSDNew(int maxNumOfBloch, double[] voltages, Matrix3D rotation, double[] thickness, Vector3DBase[] beamDirections, Solver solver = Solver.Auto, int thread = 1, bool useNonLocalAbsorption = false, bool includeTDSBackground = false, bool includeAbsorbedFluxBackground = false)`

| 引数・プロパティ | 意味 | 既定・備考 |
|---|---|---|
| `maxNumOfBloch` | ビーム（Bloch 波）数 $N$ | `Find_gVectors` は $N+1$ 本を選び末尾の同点組を落とす（[bloch-wave.md](bloch-wave.md)） |
| `voltages` | 加速電圧（= 電子のエネルギー）[kV] の配列 | 電圧ごとに全計算をやり直す |
| `rotation` | 結晶回転 `BaseRotation` | EBSD 側（`EBSD/EBSD.cs`）は `Matrix3D.IdentityMatrix` を渡す → 方向は結晶固定の直交座標系で与えられる |
| `thickness` | 深さの配列 $t_k$ [nm] | 各 $t_k$ で $0\ldots t_k$ の**累積**強度を返す（§2.8） |
| `beamDirections` | 出射方向 $\hat{\mathbf d}$（単位ベクトル）の配列 | **+Z 半球と −Z 半球の正方配列 2 枚を連結**した長さ $2G^2$ でなければ例外（`InvalidOperationException`）。六方格子の無効セルは `null` |
| `solver` | 固有値ソルバ | `Auto`: native 有効なら `Eigen_Eigen`、無効なら `Eigen_MKL`。どちらもスレッド数は `ProcessorCount` に上書き（260712Cl）。`MtxExp_*` は EBSD では MathNet の EVD 経路になる |
| `useNonLocalAbsorption` | 吸収ポテンシャルを非局所形 $U'(\mathbf g,\mathbf h)$ にする | false。§2.7 |
| `includeTDSBackground` | 局所点源を非局所源 $(2\pi/k)\psi^\dagger U'_\mathrm{back}\psi$ で**置換** | false。名前に反して「加算」ではない（260919Cl）。§2.5 |
| `includeAbsorbedFluxBackground` | 吸収で失われた流束を拡散背景として再注入 | false。§2.6 |
| `SourceSmearingSigmaNm`（static） | 局所源の Gaussian ぼかし幅 σ [nm] | `NaN`（無効）。有限値なら `IncludeTDSBackground` も強制的に true（試行機能、260919Cl）。§2.5(c) |
| `EbsdRepresentativeDirectionBlockSize` | 近傍方向でビーム集合を共有するブロック幅 | プロパティ既定 3。EBSD 側は格子幅 $G\le512$ で 1、$\le1024$ で 3、それ以上 5（`EBSD/EBSD.cs`）。§2.3 |
| `UseLocalSurfacePerBeamDirection` | 旧経路用の局所表面フラグ | 既定 false、EBSD 側も false。新経路は後述のとおり常に局所表面 $\hat{\mathbf n}=-\hat{\mathbf d}$ を使う |
| `ElasticIonModel` / `AbsorptionElasticModel`（static） | 弾性散乱因子のイオン模型 / 吸収の被積分関数に入る弾性散乱因子の模型 | `Neutral` / `TemariHybrid`。[scattering-factors.md](scattering-factors.md) |

### 2.2 幾何 — 相反定理による出射方向の扱い

出射方向 $\hat{\mathbf d}$ ごとに、検出器側から結晶へ入る「反転波」を解く（相反定理）。各方向で表面を接平面とみなし、内向き法線を

$$\hat{\mathbf n} = -\hat{\mathbf d}$$

とする（`GetCrystalFixedMasterEquivalentBeamAndSurface`、260321Ch）。`getVecK0(kvac, u0, beamDirection = d̂, surface = n̂)` は

$$\mathbf K_0 = k_\mathrm{vac}\hat{\mathbf d} + x\hat{\mathbf n},\qquad x = \sqrt{b^2+U_0}-b,\quad b = \hat{\mathbf n}\cdot k_\mathrm{vac}\hat{\mathbf d} = -k_\mathrm{vac}$$

を解くので、$\hat{\mathbf n}=-\hat{\mathbf d}$ では

$$\mathbf K_0 = -\sqrt{k_\mathrm{vac}^2+U_0}\;\hat{\mathbf d}$$

となる（二次方程式の 2 根のうち、結晶内へ向かう側が選ばれる。`KikuchiProfileCalculator` のヘッダ注記にも同じ閉形式が書かれている）。深さ $z$ は $\hat{\mathbf n}$ 方向、すなわち**出射方向に沿った経路長**であり、平らな試料表面からの深さではない（実試料の深さとの対応は EBSD 側の責務。`MeanAbsorptionLengthNm` の注記も「経路長」）。

`GetCrystalFixedMasterEquivalentBeamAndSurface` は「局所接平面基底 $(\hat u_l,\hat v_l)$ を基準基底 $(\hat u_r,\hat v_r)$ へ写す回転 $R$」を作り、$R^{\mathsf T}$ で基準の表面・ビーム方向を戻している。`GetSurfaceTangentialAxes` が $\hat u\times\hat v=\hat{\mathbf n}$ を満たすので、返り値は（丸め誤差を除き）常に $(\hat{\mathbf d},\,-\hat{\mathbf d})$ である（§4 疑問点 3）。

### 2.3 ビーム集合の選択と方向の間引き

1. **対称性による代表方向**（`CreateMasterPatternSymmetryRepresentativeDirectionMapping`、260327Ch / 260331Cl）  
   方向格子上の厳密な対称操作（`MasterPattern.GetMasterPatternSquareSymmetryOperations` / `GetMasterPatternHexSymmetryOperations`、実体は `EBSD/EBSD.cs`）で各方向を写し、軌道内で最小の線形インデックスを代表とする。代表だけを動力学計算し、最後に軌道内の他の方向へ結果をコピーする。六方格子（`MasterPattern.ShouldUseHexGrid`）では $\lvert u\rvert,\lvert v\rvert,\lvert u+v\rvert\le N$ を満たさないセルを `-1` としてスキップする。
2. **近傍ブロックでのビーム集合の共有**（`EbsdRepresentativeDirectionBlockSize` = $b$）  
   代表方向の格子座標 $(w,h)$ を $b\times b$ ブロックに割り、ブロック中心 $(\lfloor w/b\rfloor b + \lfloor b/2\rfloor,\ \ldots)$（端は $G-1$ で打ち切り）の方向で `Find_gVectors` を一度だけ呼ぶ。共有するのは **ビーム集合・ポテンシャル行列・原子位相因子** で、$\mathbf K_0$・$P_\mathbf g$・$Q_\mathbf g$ は各方向で `reset_gVectors` により計算し直す（260321Cl のバグ修正: 旧実装は代表の $\mathbf K_0$ を全方向に流用していた）。$b=1$ なら近似なし。
3. **前方波のみ残す**: 各方向で $P_\mathbf g>0$ のビームだけを元の順序で前に詰め（`gMap` に元のインデックスを保持）、ポテンシャル行列は事前計算した $b\times b$ 共有行列から部分抽出する（`getEigenMatrix(bLen, beams, gMap, baseLen, basePotential, ...)`、260602Cl。旧実装は長さ不一致で毎回局所形を再計算し、非局所吸収が静かに失われていた）。
4. 入射条件は $\boldsymbol\psi(0)=\mathbf e_0$。$\mathbf e_0$ は `Find_gVectors` が Rating $=\lvert\mathbf g\rvert Q_\mathbf g^2$ の昇順で並べた先頭（$Q_{000}=0$ なので (000)）。

### 2.4 固有値問題と励起振幅

- native（`NativeWrapper.Enabled` → `EigenEnabled`）かつ `Solver.Eigen_Eigen`: `NativeWrapper.EigenSolver` → `NativeWrapper.PartialPivLuSolve(C, e0)` で $\alpha$。
- それ以外: MathNet `DenseMatrix.Evd(Symmetricity.Asymmetric)` と `LU().Solve`。
- 作業配列は `ThreadLocal` + `ArrayPool` で方向をまたいで再利用（260321Ch）。

### 2.5 後方散乱源と強度

#### (a) 局所点源（既定）

原子 $n$（対称等価位置まで展開、占有率は $\sigma_n$ に含める）での波動場 $\psi(\mathbf r_n,z) = \sum_\mathbf g \psi_\mathbf g(z)e^{2\pi i\mathbf g\cdot\mathbf r_n}$ を用い、後方散乱電子の発生率を

$$\frac{dI}{dz} = \sum_n \sigma_n\,\lvert\psi(\mathbf r_n,z)\rvert^2,\qquad I(t) = \int_0^t\frac{dI}{dz}dz$$

とする。後方散乱の重み（260919Cl で物理スケール化。旧実装は無次元の $Z^{1.7}\cdot\mathrm{Occ}$ で、非局所項と 7〜8 桁ずれていた）:

$$\sigma_n = \frac{2\pi}{k_\mathrm{vac}}\cdot\frac{\gamma_\mathrm{rel}}{\pi V}\cdot\mathrm{Occ}_n\cdot f'_{n,\mathrm{back}}(0)\quad[\mathrm{nm}^{-1}]$$

- $f'_{n,\mathrm{back}}(0)$ = `ElectronScatteringPeng[Z][sub].FactorImaginaryAnnular(kV, 0, m, π/2, π, nTheta = 30, nPhi = 12)`（散乱角 $\theta\in[\pi/2,\pi]$ の後方半球で積分した TDS 吸収形状因子、$\mathbf g=\mathbf h=0$）。`sub` は `ElasticIonModel == Neutral` なら 0。
- $m = \max(B_\mathrm{eff},\ 10^{-3}\ \mathrm{nm^2})$。$B_\mathrm{eff}$ = `Dsf.BisoEffective`（等方なら $B_\mathrm{iso}$、非等方なら Hamilton 式の $B_\mathrm{eq}$、261105Cl）。下限 `DiffuseScatteringFactor.BisoFloorNm2 = 1E-3` nm²（0.1 Å²）は B 未設定で源が 0 になる非物理を避けるため（後方散乱角では $1-e^{-2M}$ がほぼ 1 に飽和するので結果への影響は小さい、とコード注記）。**伝播側の吸収（`getU`）にはこの下限は入らない。**
- 同じ定義で $\sum_n\sigma_n = (2\pi/k)\,U'_{0,\mathrm{back}}$（下限が効かない場合）。平面波極限で $I(t)=\mu_\mathrm{back}t$。

Bloch 展開を代入すると（`EBSDSolverManaged`）

$$\beta_{nj} = \alpha_j\sum_\mathbf g C_{\mathbf g j}e^{2\pi i\mathbf g\cdot\mathbf r_n},\qquad
S_{jj'} = \sum_n\sigma_n\,\beta_{nj}\,\beta_{nj'}^{*}\quad(S = B^{\mathsf T}\operatorname{diag}(\sigma)\,B^{*}),$$

$$I(t) = \operatorname{Re}\sum_{jj'}S_{jj'}F_{jj'}(t),\qquad I\leftarrow\max(0,I).$$

260922Cl: 旧コメントの $S=B^\dagger\operatorname{diag}(\sigma)B$ は $S$ の複素共役で、ネイティブ版がその式で実装されていた（Temari の監査で発覚。managed 版の実装自体は元から上式）。コミット 585e395 の記述では native 側は Crystallography.Native で修正済み（§5 要確認）。

- native 有効時は `NativeWrapper.EBSDSolver`（C++/Eigen、ソースはこのリポジトリに無い）、無効時は `EBSDSolverManaged`。
- `EBSDSolverManaged` は $B$ 全体を持たず、原子ごとに $\beta_n$ を `stackalloc` して $S$ に累積する。

#### (b) 非局所源（`IncludeTDSBackground = true`）

点源の代わりに、後方半球で積分した吸収ポテンシャル行列 $U'_\mathrm{back}$ を源として使う（STEM-HAADF の $\int t_c^\dagger U'_\mathrm{det}t_c\,dz$ と同じ枠組み、260316Cl）:

$$\frac{dI}{dz} = \frac{2\pi}{k_\mathrm{vac}}\,\boldsymbol\psi^\dagger(z)\,U'_\mathrm{back}\,\boldsymbol\psi(z),\qquad
(U'_\mathrm{back})_{\mathbf g\mathbf h} = \Bigl[U'^{\,\texttt{getU}}_{\mathbf g-\mathbf h}(\theta\in[\pi/2,\pi])\Bigr]^{*}$$

（`CreateMasterPatternMuBack`: `getU(voltage, beams[row] − beams[col], null, π/2, π, 30, 12).Imag.Conjugate()`。**局所形**（$\mathbf g-\mathbf h$ のみに依存）の annular 積分で、§2.7 の「非局所吸収」とは別物。ここでの「非局所」は源が点でないという意味）。

$$M_{jj'} = \frac{2\pi}{k_\mathrm{vac}}\,\alpha_j^{*}\,(C^\dagger U'_\mathrm{back}C)_{jj'}\,\alpha_{j'},\qquad
I_\mathrm{TDS}(t) = \operatorname{Re}\sum_{jj'}M_{jj'}F_{j'j}(t).$$

$\psi^\dagger U\psi$ の指数は $2\pi i(\gamma_{j'}-\gamma_j^{*})=\lambda_{j'j}$ なので $F_{j'j}$ と組む（260922Cl 修正。旧実装と native `_EBSDSolverWithTDS` は $F_{jj'}$ と組んでいた）。

- 260919Cl 以降、この値は局所源を**置換**する（`Disk への格納` で `value = ebsdBackground?[r]?[t] ?? signal`。旧実装は両者を加算しており同じ過程の二重計上だった）。
- native 有効時は `NativeWrapper.EBSDSolverWithTDS` が $I$ と $I_\mathrm{TDS}$ を同じ $F$ で一括計算、無効時は `EBSDSolverManaged` + `ComputeTDSMatrixBackground`。

#### (c) ぼかした局所源（`SourceSmearingSigmaNm` 有限、試行）

`CreateSmearedLocalSourceMatrix` が源行列

$$(U_\mathrm{src})_{\mathbf g\mathbf h} = \frac{1}{2\pi/k_\mathrm{vac}}\sum_n\sigma_n e^{2\pi i(\mathbf h-\mathbf g)\cdot\mathbf r_n}\;\exp\!\bigl(-2\pi^2\sigma^2\lvert\mathbf g-\mathbf h\rvert^2\bigr)$$

を作り、(b) と同じ行列経路 $(2\pi/k)\psi^\dagger U_\mathrm{src}\psi$ で評価する。$\lvert\mathbf g-\mathbf h\rvert$ は 2π なしの逆格子なので、実空間の等方 Gaussian（標準偏差 σ）の Fourier 変換は $\exp(-2\pi^2\sigma^2q^2)$（DW 因子 $\exp(-Bs^2)$、$B=8\pi^2\langle u^2\rangle$ と同形、とコード注記）。$\sigma=0$ で (a) と一致する（同値性の確認用）。$\sigma$ が有限のときは `IncludeTDSBackground` の非局所源より優先する。

仮説（コード注記）: 表面近くの最後のイベントの多くは小角散乱で源が $\sim1/(2\pi s)$ に広がり 1s 状態（0.2 Å）を分解できない。σ ≈ 0.5〜1 Å で晶帯軸だけが暗くなるなら仮説を支持する。

### 2.6 吸収フラックスの再注入（`IncludeAbsorbedFluxBackground = true`、260919Cl）

`ComputeAbsorbedFluxBackground`。干渉性の波の単位胞平均密度を

$$N(z) = \sum_\mathbf g\lvert\psi_\mathbf g(z)\rvert^2 = \operatorname{Re}\sum_{jj'}T_{jj'}e^{\lambda_{jj'}z},\qquad T_{jj'} = \alpha_j\alpha_{j'}^{*}\sum_\mathbf g C_{\mathbf g j}C^{*}_{\mathbf g j'}$$

とし（$N(0)=1$）、失われた分 $1-N(z)$ が方向情報を持たない拡散電子として平均の後方散乱率 $\sum_n\sigma_n$ を持つとして

$$D(t) = \sum_n\sigma_n\int_0^t\bigl(1-N(z)\bigr)dz = \sum_n\sigma_n\,\operatorname{Re}\sum_{jj'}T_{jj'}\bigl(t - F_{jj'}(t)\bigr),\qquad D\leftarrow\max(0,D)$$

を、使われている源（非局所源があればそれ、無ければ局所源）に加算する。

数値的工夫:
- $\sum T = N(0) = 1$ を使って $t-\sum TF$ を $\sum T(t-F)$ の形で評価し、弱吸収での桁落ちを避ける。
- $\lvert\lambda t\rvert^2<10^{-6}$ では $t-F = -(\lambda t^2/2)(1+\lambda t/3+(\lambda t)^2/12)$ の級数。
- $T$ と $t-F$ の Hermite 対称性から対角 + 2×下三角だけを足す。
- $e^{\lambda_{jj'}t} = E_{jt}E^{*}_{j't}$、$E_{jt}=e^{2\pi i\gamma_j t}$ として指数を $b_\mathrm{Len}\cdot t_\mathrm{Len}$ 回に減らす。

注: $N(z)$ は**密度**の和であり、260829Cl の注記にある無吸収時の保存量 $\sum_\mathbf g P_\mathbf g\lvert\psi_\mathbf g\rvert^2$ とは異なる（§4 疑問点 2）。

### 2.7 吸収ポテンシャルの形（`UseNonLocalAbsorption`）

`RentMasterPatternPotentialMatrix`:

- false（既定）: `getPotentialMatrix(dim, b, ref A)`。$\tilde U_{\mathbf g\mathbf h}=U^{r}_{\mathbf h-\mathbf g}+iU'_{\mathbf h-\mathbf g}$（対角は $iU'_0$）。$U'$ は全立体角の `FactorImaginary`。Toeplitz 構造なので一意な $\mathbf g-\mathbf h$ ごとに 1 回の計算で済む（`uDictionary` キャッシュ）。
- true: `getPotentialMatrix(dim, b, ref A, 0, π)`。実部は局所形のまま、虚部は `getU(kV, b[col], b[row], 0, π)` → `FactorImaginaryAnnular(kV, g, h, m, 0, π, 60, 20)` の $\mathbf g,\mathbf h$ 両依存（2 次元 Gauss–Legendre 積分）。対角 $U'(\mathbf g,\mathbf g)$ も $\mathbf g$ ごとに異なる。重い元素で前方ピークが強いと差が大きい（260316Cl 注記）。コストは $N^2$ 要素ごとの数値積分。

### 2.8 厚さ配列・多電圧・出力

- 各電圧で `uDictionary.Clear()` してから $U_0$ を読む（261014Cl 修正: 旧順序では 2 つ目以降の電圧が前の電圧の $U_0$ を使っていた。キャッシュ鍵が電圧に依らないため）。$\sigma_n$ も電圧ごとに計算。
- 出力: `Disks[vIndex][t]` は `CBED_Disk`。`Amplitudes[r]` = $\sqrt{\max(0,\,I_r(t))}$（方向 $r$ の**強度の平方根**を実数で格納、位相なし）。$I_r(t)$ は深さ $0\ldots t$ の**累積**値（平面波極限で $\mu_\mathrm{back}t$）。単位は無次元（$\sigma$ [nm⁻¹] × 深さ [nm]、入射流束 1 あたり）。
- EBSD 側との受け渡し（本稿では定義のみ）:
  - 入力: エネルギー列 `Energies` [keV]（= `voltages` [kV]）、深さ列 `Depths` [nm]（= `thickness`、上記の経路長）、方向列（Rosca–Lambert 格子 2 半球）、格子幅 → `EbsdRepresentativeDirectionBlockSize`。
  - 補助: `MeanAbsorptionLengthNm(double[] energiesKeV)` が $\lambda_\mathrm{abs}(E) = k_\mathrm{vac}/(2\pi U'_0)$ [nm]（強度の減衰長、経路長）を返す（260923Cl）。$U'_0\le0$ なら $+\infty$。MC 側の生存確率とマスターパターン平均吸収の二重計上を除くのに使う（`EbsdMonteCarloDistribution.FillPathLengthWeights`）。
- 進捗は全格子基準（代表方向は軌道サイズで重み付け）で 50 方向ごとに報告。

### 2.9 キャッシュ上の注意

`CreateMasterPatternMuBack` は `getU(..., null, π/2, π, 30, 12)` を呼ぶが、`uDictionary` の鍵は $(\mathbf g\text{ の hkl},\ \texttt{int.MaxValue})$ で、通常の局所 $U$ と**同じ鍵**になる。そのため muBack 段の前後で `uDictionary.Clear()` している。ポテンシャル行列はこの段より前に作り終えているので結果は正しいが、同じインスタンスで同時に別の `getU` を走らせると後方半球の値を読む（§4 疑問点 9）。

---

## 3. 菊池バンドの 1D 動力学プロファイル（`Diffraction/Kikuchi/`、260805Cl〜）

設計正本はコード注記によれば `ReciPro_菊池線動力学化設計.md`（リポジトリ外）。§番号はその設計書のもの。

### 3.1 位置づけと幾何

- **透過（前方）菊池バンド**用。実験室系でビームは $\hat{\mathbf b}=(0,0,-1)$ へ進む。出射方向 $\hat{\mathbf d}$ は $\hat{\mathbf b}$ の近く（既定で $\angle(\hat{\mathbf d}(x{=}0),\hat{\mathbf b})\le0.35$ rad）。
- 反転幾何は EBSD マスターパターンと同一: $\hat{\mathbf n}=-\hat{\mathbf d}$、$\mathbf K_0=-\sqrt{k_\mathrm{vac}^2+U_0}\,\hat{\mathbf d}$。
- 源の運動量移行（反転波のビーム $\mathbf h$ ↔ 出射ビーム $-\mathbf h$）: $\mathbf q_0 = k_\mathrm{vac}(\hat{\mathbf d}-\hat{\mathbf b})$、$\mathbf q_i = \mathbf q_0-\mathbf h_i$。
- 厚さ `thickness` [nm] は一様な源分布の積分範囲 $0\ldots t$。
- 「伝播 = `getU` 虚部 / 源 = 独立に構成した $Q$」の暫定運用で、源と損失の収支は取っていない（`DisplayNormalizedTag = "source-loss not balanced"`）。

### 3.2 スナップショット（`KikuchiPotentialSnapshot`）

結晶・kV・カーネル・吸収 on/off を凍結。菊池専用の `BetheMethod` インスタンスを内部に持ち、その `uDictionary` を使う（共有の `Crystal.Bethe` には触れない）。$U_0 = $ `getU(kV).Real.Real`、$U'_0 = $ `getU(kV).Imag`（`absorptionOff` なら 0、非対角の $U'$ も 0）。鮮度判定 `Matches` は結晶参照・kV・`ComputeCrystalHash`（格子定数、Z、`SubNumberElectron`、Occ、`Biso`/`Biso000`/`B11..B31`、原子位置、`ElasticIonModel`、`AbsorptionElasticModel`）。

### 3.3 バンド族（`KikuchiBandFamily`）

- `Pair`: $\pm\mathbf g$ を 1 族にまとめる。代表は先頭の非零指数が正の側。入力順を保存。
- `CollapseSystematicRows(families, xMax)`（260807Cl）: 方向 $(h,k,l)/\gcd$ が同じ族を 1 本の系統列とみなし、$\lvert\mathbf g\rvert$ 最小の生成元だけを残す。ただし倍率 $n>x_\mathrm{max}$ の高次は生成元のプロファイル範囲外なので独立な族として残す。`RowOrder > 0` で自動適用（畳まないと二重計上。MgO（Omoto Fig 6）で corr 0.914 → 0.898 に悪化、畳むと 0.927、とコード注記）。

### 3.4 バンド座標

`BandGeometry`: $\mathbf g_\mathrm{lab}=R\,\mathbf g$、$\sin\theta_B = \lvert\mathbf g\rvert/(2k_\mathrm{vac})$、$\cos\varepsilon = \hat{\mathbf g}\cdot\hat{\mathbf b}$、$\sin\varepsilon = \sqrt{\hat g_x^2+\hat g_y^2}$（桁落ち回避、260806Cl）。`Valid` = $\sin\varepsilon>10^{-6}$ かつ $\sin\theta_B<0.5$。

プロファイル座標 $x = \sin\theta'/\sin\theta_B$、$\sin\theta'\equiv-\hat{\mathbf g}\cdot\hat{\mathbf d}$（$x\sin\theta_B$ は $\pm0.99$ でクランプ）。$\hat{\mathbf d}=\cos\psi\,\hat{\mathbf b}+\sin\psi\,\hat{\mathbf g}_\perp$ の 2 解のうち $\hat{\mathbf b}$ に近い側を取る。その点が $\hat{\mathbf b}$ から `SampleAngle`（既定 0.025 rad）より近ければ、同じ円錐 $-\hat{\mathbf g}\cdot\hat{\mathbf d}=\sin\theta'$ 上を面外へ滑らせて $\hat{\mathbf b}$ から `SampleAngle` の点で評価する（$\mathbf q_0\to0$ で TDS 源が消え $c$ が発散するのを防ぐ。Si [111] で確認）。

### 3.5 2 波計算（`ComputePoint`、`RowOrder = 0`）

族ごとに $\{0,+\mathbf g\}$ と $\{0,-\mathbf g\}$（member $=\pm1$）を別々に 2×2 で解き、和を取る。

$$A = \begin{pmatrix} iU'_0/P_0 & U_{+\mathbf g}/P_0\\ U_{-\mathbf g}/P_1 & (iU'_0+Q_1)/P_1\end{pmatrix},\quad
P_0 = 2\hat{\mathbf n}\cdot\mathbf K_0,\ P_1 = 2\hat{\mathbf n}\cdot(\mathbf K_0+\mathbf g_m),\ Q_1 = \lvert\mathbf K_0\rvert^2-\lvert\mathbf K_0+\mathbf g_m\rvert^2$$

（$U_{\pm\mathbf g}=U^r+iU'$ は `snap.GetU`。行割り = 260829Cl 修正、固有値は不変で固有ベクトルと $\alpha$ が正しくなる）。`DisableBraggCoupling` で非対角 0。$P_1\le10^{-3}P_0$ の点は寄与 0。

2×2 固有値は閉形式 $\gamma_{0,1} = \tfrac{a_{00}+a_{11}}{2}\pm\sqrt{(a_{00}-a_{11})^2/4+a_{01}a_{10}}$、固有ベクトルは $\lvert a_{01}\rvert\ge\lvert a_{10}\rvert$ で $(a_{01},\gamma-a_{00})^\mathsf T$、それ以外 $(\gamma-a_{11},a_{10})^\mathsf T$、列正規化。$\lvert\det C\rvert<10^{-14}$（defective）は $C=I,\ \alpha=(1,0)$ の素通し。$\alpha = C^{-1}(1,0)^\mathsf T$。

### 3.6 源密度行列 $Q$ とカーネル

$$Q_{ii'} = \sum_a \mathrm{Occ}_a\;\mathcal C_a(s_i^2,s_{i'}^2,s_{ii'}^2)\;e^{-2\pi i(\mathbf h_{i'}-\mathbf h_i)\cdot\mathbf r_a},\qquad s_i^2=\lvert\mathbf q_i\rvert^2/4,\ s_{ii'}^2=\lvert\mathbf h_{i'}-\mathbf h_i\rvert^2/4$$

（和は対称等価位置まで展開したサイト。位相規約は `CreatePhaseFactors` と同じで、EBSD の $S$ の $\sum_n\sigma_n e^{2\pi i(\mathbf g-\mathbf g')\cdot\mathbf r_n}$ と同型）。

`IKikuchiInelasticKernel`:
- `SourceAmplitude(a, s²)` = $\tau_a(q)$。既定の `SourceCoherence` は因子化近似 $\tau_i\tau_j$。対角整合 $\mathcal C(s^2,s^2,0)=\tau^2$ が契約。
- `KikuchiTdsEinsteinKernel`（既定）: 独立原子・等方 Einstein。
  $$\tau_a^2 = f_a(s^2)^2\bigl[1-e^{-2B_a s^2}\bigr],\qquad
  \mathcal C_a = f_a(s_i^2)f_a(s_j^2)\bigl[e^{-B_a s_{ij}^2}-e^{-B_a(s_i^2+s_j^2)}\bigr]$$
  （Omoto–Tsuda–Tanaka 2002 eq. 38 = Hall–Hirsch、とコード注記）。$f_a$ = `ElectronScatteringPeng[Z][sub].Factor`（`ElasticIonModel` 規約。`IonFull` の単極子は入れない）、$B_a$ = `Dsf.BisoEffective` [nm²]。$B$ が NaN または $\le10^{-6}$ なら `DefaultBiso = 0.005` nm²（0.5 Å²）を代用し `UsedDefaultBiso = true`（B=0 で全バンドが消える実機報告への対処）。
- `KikuchiUniformSourceKernel`: $\tau=1$（診断用）。
- `DiagonalSourceOnly` で $Q$ の非対角を 0、`SwapSourceWeights` で 2 波の $s^2$ を入れ替える（診断）。

### 3.7 強度・参照強度・コントラスト

$$S_{jj'} = \alpha_j\alpha_{j'}^{*}\sum_{ii'}C_{ij}Q_{ii'}C^{*}_{i'j'},\qquad I(t)=\operatorname{Re}\sum_{jj'}S_{jj'}F_{jj'}(t).$$

参照強度は Bragg 結合を切った直進チャネル: $\gamma_\mathrm{ref}=iU'_0/P_0$、$I_\mathrm{ref}=Q_{00}F(\gamma_\mathrm{ref},\gamma_\mathrm{ref})$（swap 診断時も参照は実背景のまま）。

$$c(x) = \frac{I(x) - n_\mathrm{ref}\,I_\mathrm{ref}(x)}{\max\bigl(I_\mathrm{ref}(x),\ 10^{-3}\max_x I_\mathrm{ref}\bigr)},\qquad n_\mathrm{ref}=\begin{cases}2&\text{2 波（}\pm\mathbf g\text{ の和）}\\1&\text{系統列}\end{cases}$$

（$t\le0$ なら 0、非有限値は 0）。`DynamicalRank` = $\mathrm{RMS}_x(c)$。

### 3.8 系統列（`RowOrder = N > 0`、260807Cl）

ビーム $\{n\mathbf g: n=-N..N\}$（$2N+1$ 波、$N\le15$ にクランプ）を 1 つの基底で解く。`RowPotential` が $U(m\mathbf g)$ と $e^{-2\pi i m\mathbf g\cdot\mathbf r_a}$（$m=-2N..2N$）をバンドごとに 1 回前計算。行列は 2 波版と同じ構成（対角 $(iU'_0+Q_i)/P_i$、非対角 $U((c-r)\mathbf g)/P_r$）。EVD は native があれば `NativeWrapper.EigenSolver` + `Inverse`、無ければ MathNet。$\alpha_j=[C^{-1}]_{j,N}$。どれか 1 本でも $P_i\le10^{-3}P_N$（または $P_N\le0$）ならその点は寄与 0。$W=CQC^\dagger$ 型の縮約は 2 段の行列積（$O(n_b^3)$）。

### 3.9 `Options` の既定値と選定

| 項目 | 既定 | 意味 |
|---|---|---|
| `SampleCount` | 129 | $x$ 格子点数（$-x_\mathrm{max}..x_\mathrm{max}$ 等間隔） |
| `XMax` | 2.5 | $\lvert x\rvert$ 上限（Bragg 半幅の倍数） |
| `MaxScatteringAngle` | 0.35 rad | $x=0$ の出射方向と $\hat{\mathbf b}$ の角度上限。超えると `Valid=false` |
| `SampleAngle` | 0.025 rad | 視軸からの最小サンプル角 |
| `RowOrder` | 0 | 0 = 2 波経路（数値は従来とビット同一）、上限 15 |
| `SwapSourceWeights` / `DisableBraggCoupling` / `DiagonalSourceOnly` | false | 診断 |

`ComputeProfiles(snap, candidates, rotation, thickness, topN, opt, previousSelection)`: 候補を並列に計算し（`RowOrder>0` なら先に `CollapseSystematicRows`）、`Valid` かつ `DynamicalRank > 0` のうち順位上位 `topN` を返す。前回選ばれた族には 1.15 倍のヒステリシス。

---

## 4. 実装上の疑問点

1. **六方格子 + ブロック共有で、ブロック中心が無効セルになりうる**（`Diffraction/BetheMethod.cs:1600-1608, 1630`、`EBSD/EBSD.cs:235, 1191-1194`）。  
   代表方向のブロック中心 $(w_c,h_c)$ は正方配列の座標で計算され、六方格子の有効性（$\lvert u+v\rvert\le N$）を検査しない。EBSD 側は $G>512$ で $b=3$ または 5 を使うので、菱形配列の境界（$\lvert u+v\rvert=N$ 付近）のブロックでは中心が無効セルに落ちうる。その方向は `null` で、`GetCrystalFixedMasterEquivalentBeamAndSurface(null, …)` は基準表面 $(0,0,-1)$ を返すため、そのブロックの全方向が **[00(−1)] 方向で選んだビーム集合**を共有することになる（$P,Q$ は方向ごとに再計算されるが、ビーム集合はエワルド球から遠い）。反例の作り方: 六方晶で `GridSize = 1025`（$b=3$）とし、$(u,v)$ が有効で $(u_c,v_c)$ が無効になるブロックを探す。
2. **吸収フラックス再注入の $N(z)$ が密度和**（`BetheMethod.cs:2108-2118, 2120-2168`）。  
   260829Cl の注記どおり、行割り後の Bethe 方程式で無吸収時に保存するのは $\sum_\mathbf g P_\mathbf g\lvert\psi_\mathbf g\rvert^2$（法線方向の流束）であり、$\sum_\mathbf g\lvert\psi_\mathbf g\rvert^2$ ではない。$U'=0$ でも $P_\mathbf g\ne P_0$ のビームがあれば $N(z)$ は振動し、$D(t)\ne0$（`Max(0, …)` でクランプ）になる。「失われた流束」を測るなら $\sum_\mathbf g (P_\mathbf g/P_0)\lvert\psi_\mathbf g\rvert^2$ が自然と思われる。検証案: `UseNonLocalAbsorption=false` で全原子 $B=0$（$U'\equiv0$）にし、`IncludeAbsorbedFluxBackground=true` の $D(t)$ が 0 になるか見る（ただし $\sigma_n$ は下限 $B$ で非零なので $\sum\sigma_n>0$）。
3. **`GetCrystalFixedMasterEquivalentBeamAndSurface` は恒等変換に帰着する**（`BetheMethod.cs:815-825`）。  
   `GetSurfaceTangentialAxes` は $\hat u\times\hat v = \hat{\mathbf n}$ を満たすので、$R$ は $\hat{\mathbf n}_l\mapsto\hat{\mathbf n}_r$ を満たし、$R^{\mathsf T}\hat{\mathbf n}_r = -\hat{\mathbf d}$、$R^{\mathsf T}(-\hat{\mathbf n}_r)=\hat{\mathbf d}$。返り値は基準表面に依らず $(\hat{\mathbf d},-\hat{\mathbf d})$。結果は正しいが、コメント（「逆回転を beam / surface 側へ押し戻す」）から想像される処理とは異なる。接平面基底の向きは返り値に含まれない。
4. **菊池プロファイル座標の符号ラベル**（`Diffraction/Kikuchi/KikuchiProfileCalculator.cs:362`）。  
   コメントは「x = +1 が +g の Bragg 線」とするが、$\mathbf K_0=-K\hat{\mathbf d}$ と $Q_{\mathbf g_m}=\lvert\mathbf K_0\rvert^2-\lvert\mathbf K_0+\mathbf g_m\rvert^2 = 2K\hat{\mathbf d}\cdot\mathbf g_m-g^2$ から、member $+\mathbf g$ の $Q=0$ は $\hat{\mathbf d}\cdot\hat{\mathbf g}=g/2K$、すなわち $x = -k_\mathrm{vac}/K\approx-1$。「出射ビーム $-\mathbf h$ ↔ 反転波 $\mathbf h$」の読み替えで出射側の $+\mathbf g$ を指すなら正しい。ReciPro 側の描画（E/D 非対称の向き）がどちらを前提にしているか確認が要る。あわせて $\sin\theta_B$ は $k_\mathrm{vac}$ で定義されるため、Bragg 線は厳密には $\lvert x\rvert = k_\mathrm{vac}/K = 1-O(U_0/2k^2)$ にある（実害は小さい）。
5. **260829Cl の行割り修正後も古い式がコメントに残る**（`KikuchiProfileCalculator.cs:7, 592`）。ヘッダと `ComputeRowPoint` の summary が非対角を `/P_col` と書いている（実装は `/P_row`）。コメントのみ。同様に `Diffraction/Alchemi/AlchemiReduction.cs:12` は $S=B^\dagger\operatorname{diag}(\sigma)B$ と書いている（260922Cl で訂正された旧表記。分担 D の範囲）。
6. **`ComputeCrystalHash` が `Dsf.UseIso` を含まない**（`KikuchiProfileCalculator.cs:94-113`）。カーネルが使う `BisoEffective` は `UseIso` で `Biso` と `Biso000` を切り替えるので、値を変えずに `UseIso` だけ切り替えるとスナップショットが新鮮と判定される。
7. **2 波経路で $I_\mathrm{ref}$ を member $+1$ からしか取らない**（`KikuchiProfileCalculator.cs:261-267`）。$+\mathbf g$ 側が $P_1$ ガードで `default`（$I=I_\mathrm{ref}=0$）を返し $-\mathbf g$ 側が有効な場合、$I_\mathrm{ref}=0$ で床値割りになり $c$ が跳ねる。`BandGeometry` の不変条件（$\sin\theta_B<0.5$ と $\pm0.99$ クランプ）で $P_1>0.01P_0$ が保証され現状は到達しないとコメントにあり、潜在的な問題にとどまる。
8. **B 未設定原子の扱いが源と伝播で異なる**（設計上の確認事項）。EBSD の $\sigma_n$ は $B\ge10^{-3}$ nm² の下限、菊池カーネルは `DefaultBiso` = $5\times10^{-3}$ nm² を代用、伝播の吸収 `getU` は下限なし（$B=0$ なら $U'=0$）。結果として B=0 の原子は「吸収しないが後方散乱する」。意図どおりか確認したい。
9. **`uDictionary` を共有することによる副作用**（`BetheMethod.cs:860, 1641-1658, 3951-3963`）。muBack は通常の $U$ と同じ鍵で annular 値をキャッシュする（前後で Clear）。`MeanAbsorptionLengthNm` もキャッシュを Clear する（末尾で 2 回連続の Clear は冗長）。同じインスタンスで計算中にこれらが並行して呼ばれると誤った $U$ を読む。現状の呼び出し順では問題にならない。
10. **$F_{jj'}$ の小さい $\lambda$ での桁落ち**（`BetheMethod.cs:2100, 2258`、`KikuchiProfileCalculator.cs:733`）。判定が $\lvert\lambda\rvert^2<10^{-30}$ と $\lambda$ 自身に対するもので、$10^{-15}\lesssim\lvert\lambda t\rvert\ll1$ では $(e^{\lambda t}-1)/\lambda$ の相対誤差が $\sim\epsilon/\lvert\lambda t\rvert$ になる。`ComputeAbsorbedFluxBackground` だけは級数で対処済み。通常の吸収（$\lvert\lambda t\rvert\gtrsim10^{-3}$）では実害なし。

---

## 5. 要確認事項

1. native `Crystallography.Native` の `_EBSDSolver` が $S=B^{\mathsf T}\operatorname{diag}(\sigma)B^{*}$（managed と同じ）に、`_EBSDSolverWithTDS` が $F_{j'j}$ に修正済みで、同梱 DLL がその版か。ソースがこのリポジトリに無く確認できない（コミット 585e395 は「fixed in Crystallography.Native」と記載）。
2. 後方散乱の積分範囲 $\theta\in[\pi/2,\pi]$ の基準軸（反転波の進行方向か、元の入射ビームか）の物理的な意味づけ。`FactorImaginaryAnnular` の角度定義は [scattering-factors.md](scattering-factors.md) 側の記述と合わせて確認したい。
3. 局所点源は原子の平均位置での $\lvert\psi(\mathbf r_n)\rvert^2$ を使い、熱振動による源の広がりは $\sigma_n$（$1-e^{-2M}$ を含む断面積）にしか入らない。これが意図した近似か（`SourceSmearingSigmaNm` が試行中の代替）。
4. `IncludeTDSBackground` の非局所源は局所形の annular $U'_\mathrm{back}(\mathbf g-\mathbf h)$ を使う。`UseNonLocalAbsorption` と同様に $U'_\mathrm{back}(\mathbf g,\mathbf h)$ の非局所形にする予定があるか。
5. `EbsdRepresentativeDirectionBlockSize > 1` の精度（ブロック端でのビーム選択の劣化）を評価した記録の有無。
6. 菊池プロファイル（§3）の設計書 `ReciPro_菊池線動力学化設計.md` と Omoto–Tsuda–Tanaka 2002 の書誌（コード中には書名・式番号のみ）。
7. `Solver.MtxExp_*` を指定したとき EBSD 経路が EVD にフォールバックすることが意図どおりか。
8. 菊池 2 波経路で $\pm\mathbf g$ を別々の 2×2 で解いて和を取る近似（$+\mathbf g$ と $-\mathbf g$ の同時励起を無視）の適用範囲。`RowOrder ≥ 1` がその改善版という理解でよいか。
9. `ComputeProfile` の $c$ の分母の床 $10^{-3}\max I_\mathrm{ref}$（設計 §4 の ε の実装）の数値が設計どおりか。
10. EBSD 側（別文書）で `Amplitudes` を二乗して強度に戻しているか、深さ方向の重みを累積値の差分で作っているか（受け渡しの約束事）。
11. `beamDirections` が 2 半球の連結でなければ例外となる仕様（単一半球の検出器パターン計算は旧経路 `ebsd_DoWork` がお蔵入りしたため現状は不可）が意図どおりか。

要確認事項: 11 件。
