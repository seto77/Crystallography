本稿はコードから Claude が起こした下書き（2026-11-05）。著者未査読。

# Bloch 波（Bethe）法 — ポテンシャル・固有値問題・CBED・回折振幅・PED

対象: `Diffraction/BetheMethod.cs` の `BetheMethod` クラスのうち、EBSD（[kikuchi-dynamical.md](kikuchi-dynamical.md)）と ALCHEMI（[alchemi-ionization.md](alchemi-ionization.md)）以外。
HRTEM・STEM・ポテンシャル像は [bloch-wave-imaging.md](bloch-wave-imaging.md) に分けた。
散乱因子そのもの（`AtomStatic.ES` の Peng 係数、吸収形状因子 `FactorImaginary*`）は [scattering-factors.md](scattering-factors.md) を参照。ここでは `getU` の中での使い方だけを書く。

行番号は 2026-11-05 時点の `BetheMethod.cs`。

---

## 1. 規約（全体）

| 項目 | 実装 | 根拠 |
|---|---|---|
| 長さの単位 | **nm**（逆格子は nm⁻¹） | クラス summary「すべて、単位はnm」 |
| 波数 | $k = 1/\lambda$（2π を含まない）[nm⁻¹] | `UniversalConstants.Convert.EnergyToElectronWaveNumber`（summary「2πで割った数値では無い」） |
| 加速電圧 | kV（`AccVoltage`） | |
| 平面波 | $\exp(+2\pi i\,\mathbf{k}\cdot\mathbf{r})$ | 伝播因子 `Exp(TwoPiI * γ * t)`、像合成 `Exp(TwoPiI * g·r)` |
| 散乱角（絞り・検出器・PED 半角） | **ラジアン** | `Math.Sin(angle)/λ`、`Matrix3D.Rot(axis, θ)` |
| 座標系 | 画面右 +X、上 +Y、手前 +Z。入射は既定で $-Z$ | `Surface` の summary、`getVecK0` |
| 表面法線 `Surface` | 試料表面から**内部へ向かう**単位ベクトル、既定 $(0,0,-1)$ | `Surface` の summary |
| 厚さ $t$ | 表面法線方向の深さ（ビーム方向の経路長ではない） | 260829Cl コメント（`cbed_DoWork`、L503 付近）: 1/cosτ を掛けると二重補正になるので掛けない |

### 1.1 波長と相対論補正

`UniversalConstants.Convert`:

$$
\lambda\,[\mathrm{nm}] = \frac{1.2264262862108}{\sqrt{V\,(1 + 0.97847537252\times10^{-6}\,V)}},\qquad V = 1000\times\text{kV}
$$

（コメントの式 $\lambda = h/\sqrt{2m_0 eV(1+eV/2m_0c^2)}$ を数値化したもの。）`EnergyToElectronWaveNumber` はその逆数、`WaveLengthToElectronEnergy` は $aE^2 + bE - c^2/\lambda^2 = 0$（E は keV）の正根。

相対論因子は `getU` 内で

$$
\gamma = 1 + \frac{e_0 \cdot 10^3\,\text{kV}}{m_0 c^2}
$$

（`UniversalConstants.e0, m0, c2`）。定数は `m0 = 9.1093897E-31`, `e0 = 1.60217733E-19`（CODATA 1986 相当の値）、`h = 6.62606896E-34`（2006）。波長の係数 1.22642629 と 0.978475e-6 は別途の数値で、`gamma` と同じ定数組から毎回計算しているわけではない。

`WaveProperty`（`Diffraction/WaveProperty.cs`）は波の種類・波長 [nm]・エネルギー [eV]・特性 X 線の元素と線種を保持するだけの入れ物で、変換は上の `Convert` を呼ぶ。電子線では `Energy` は eV で、波長変換時に `/1000` して kV に直す。260718Cl で特性 X 線コンストラクタの Å→nm 換算漏れと Energy への波長代入のバグが直されている。

`DiffractionOptics`（`Diffraction/DiffractionOptics.cs`）は IPAnalyzer/ReciPro 用の検出器パラメータ XML（文字列の入れ物と旧タグ名の置換）を読むだけで、物理計算は含まない。

`DynamicTheory`（`Diffraction/DynamicTheory.cs`）は実装がすべてコメントアウトされた旧マルチスライスのクラスで、`[Obsolete]`（260718Cl、全リポで未使用、マクロ互換のため残置）。ReciPro の動力学計算はすべて `BetheMethod` が担う。

---

## 2. ポテンシャルのフーリエ係数 `getU`

`getU(in double kV, in Beam g, in Beam h = null, double inner = NaN, double outer = NaN, int nTheta = 60, int nPhi = 20)`（L3863–3934）

戻り値 `(Complex Real, Complex Imag)` は

- `Real` $= \hat U_\mathbf{g}$: 弾性（実）ポテンシャル $U(\mathbf r)$ のフーリエ係数。**複素数**（非対称心結晶で位相を持つ）。
- `Imag` $= U'_\mathbf{g}$: 吸収（虚）ポテンシャル $U'(\mathbf r)$ のフーリエ係数。これも複素数。

結晶内の全ポテンシャルは $U(\mathbf r) + iU'(\mathbf r)$ で、行列要素では `Real + i·Imag` として組む（§3）。

### 2.1 実装式

$s^2 = |\mathbf g|^2/4$（$s = \sin\theta/\lambda$、nm⁻²）、$h$ を与えた非局所形式では $s^2 = |\mathbf g - \mathbf h|^2/4$、指数は $\mathbf g - \mathbf h$。

$$
\texttt{Real} = \frac{\gamma}{\pi V_\mathrm{cell}} \sum_{\text{site } a}\; \sum_{j\in a} \mathrm{Occ}_a\, f^{e}_a(s)\, e^{-M_{aj} s^2}\, e^{+2\pi i\,\mathbf g\cdot \mathbf r_j}
$$

$$
\texttt{Imag} = \frac{\gamma}{\pi V_\mathrm{cell}} \sum_{a}\sum_{j\in a} \mathrm{Occ}_a\, f'_a(s; M_{aj})\, e^{-M_{aj} s^2}\, e^{+2\pi i\,\mathbf g\cdot \mathbf r_j}
$$

- $V_\mathrm{cell}$ = `Crystal.Volume` [nm³]、$f^e$ [nm] → $U$ は **nm⁻²**。コード中のコメント「係数については、Kirklandの教科書のp120参照」。
- $U$ は Schrödinger 方程式を $[\nabla^2 + 4\pi^2(k^2 + U)]\psi = 0$ と書いたときの $U = 2m e V/h^2$（$m=\gamma m_0$）に相当（`getVecK0` が $|\mathbf K|^2 = k^2 + U_0$ を使うことと整合。§4）。
- **位相の符号は $+2\pi i\,\mathbf g\cdot\mathbf r$**（20240524 のコメント「位相項の符号をプラスに変更（これで、対称心の結晶の計算が上手くいくはず 三菱・中村）」）。結晶学で普通の $V_\mathbf g = V_\mathrm{cell}^{-1}\int V e^{-2\pi i\mathbf g\cdot\mathbf r}$ とは逆符号なので、**`getU(g)` は教科書の $U_{-\mathbf g}$ に等しい**。行列要素は `getU(b[col] - b[row])` で組むので（§3）、要素 $(g,h)$ は教科書どおりの $U_{\mathbf g-\mathbf h}$ になる。この「符号を逆にした係数を逆向きの差で呼ぶ」構造が、STEM の `.Imag.Conjugate()`（[bloch-wave-imaging.md](bloch-wave-imaging.md) §2.5）と STEM-EDX の負位相の直接生成（260801Cl, `ComputeIonizationU`）の背景にある。
- 和は `Crystal.Atoms`（独立サイト）× `atoms.Atom`（等価位置すべて）。占有率 `Occ` を掛ける。

### 2.2 弾性散乱因子とイオンの扱い

`AtomStatic.ElectronScatteringPeng[Z][sub]`（Peng の 5 Gauss）。`es.Factor(s2)` は **s² を nm⁻² で受け、f を nm で返す**（内部で ×0.01, ×0.1）。8 Gauss（`ElectronScatteringEightGaussian`）はコメントアウトで不使用。

静的設定 `BetheMethod.ElasticIonModel`（既定 `Neutral`、260613Cl）:

| モデル | `sub` | 単極子 |
|---|---|---|
| `Neutral`（既定） | 0（中性原子、割当イオンを無視） | なし |
| `IonScreened` | `atoms.SubNumberElectron` | なし（Peng イオン係数の有限部分のみ。旧既定） |
| `IonFull` | `atoms.SubNumberElectron` | $g\neq0$ に Mott–Bethe 単極子 $0.2393366\,\Delta Z/s^2$ [nm]（$s^2$ は nm⁻²）を加える |

`IonFull` は単位胞が電荷中性（`|NetCellCharge()| ≤ 1e-6`、$\sum \text{Multiplicity}\cdot\text{Occ}\cdot\text{Valence}$）のときだけ単極子を足し、非中性なら自動的に `IonScreened` 相当に落ちる（260614Cl）。$s^2 = 0$（$U_0$, 平均内部ポテンシャル）には単極子を足さない。

### 2.3 Debye–Waller 因子

`atoms.Dsf`（`DiffuseScatteringFactor`、詳細は [scattering-factors.md](scattering-factors.md)）。

- 等方: $M = B_\mathrm{iso}$ [nm²]、$e^{-B s^2}$（$s^2 = g^2/4$）。$B = 8\pi^2 U$ への換算は `Dsf.Biso` 側。
- 非等方・$g\neq0$: 等価位置 $j$ ごとに指数を `atom.Operation.ConvertPlaneIndex(index)` で回し

  $$M_{aj} = \frac{\beta_{11}H^2 + \beta_{22}K^2 + \beta_{33}L^2 + 2\beta_{12}HK + 2\beta_{23}KL + 2\beta_{31}LH}{s^2}$$

  （L3912。`dsf.B11…` は無次元の β 形式と読める）。$e^{-Ms^2}$ が $e^{-\beta_{ij}h_ih_j}$ になるように「g 方向の実効等方 B」に直してから使う。
- 非等方・$g = 0$: `dsf.BisoEffective`（261105Cl。非等方では Hamilton の $B_\mathrm{eq}$ = `Biso000`）。
- `M` が NaN なら 0。`dsf.IsZero` なら吸収 0。

### 2.4 吸収（虚部）

`AbsorptionSource ImaginaryPotentialAbsorption => AbsorptionSource.TdsRedistributable`（L3848）: **虚部の出所は現状 TDS（熱散漫散乱）だけ**。平均吸収や経験的減衰は足していない（足したらここを更新せよ、という注記つき）。

$f'$ は `AtomStatic.ES` の次のメソッドで、戻り値は $f^e$ と同じ nm 単位（コメントには「答えは無次元」とあるが、`FactorImaginary` の末尾は $\gamma k_0/2 \times \int f f\,(\ldots)\,d\Omega \times 0.01$ で、f[Å]²×k₀[nm⁻¹]×0.01 = nm）:

| 呼び方 | 使うメソッド | 意味 |
|---|---|---|
| `inner`/`outer` が NaN（既定） | `es.FactorImaginary(kV, s2, M)` | 全立体角の TDS 吸収（局所形式 $U'(\mathbf g)$） |
| `inner`,`outer` あり, `h = null` | `es.FactorImaginaryAnnular(kV, g.Vec, M, inner, outer, nTheta, nPhi)` | 散乱角 [inner, outer] rad の環状域に限った局所形式（STEM ADF） |
| `inner`,`outer` あり, `h` あり | `es.FactorImaginaryAnnular(kV, g.Vec, h.Vec, M, …)` | 非局所形式 $U'(\mathbf g,\mathbf h)$（EBSD の `UseNonLocalAbsorption`） |

被積分関数は $f(|\mathbf k-\mathbf g/2|)f(|\mathbf k+\mathbf g/2|)\,[1-e^{M(s^2-s_1^2-s_2^2)}]$ で、`getU` 側でさらに $e^{-Ms^2}$ を掛けるので合わせて Hall–Hirsch 型の $e^{-Ms^2} - e^{-M(s_1^2+s_2^2)}$ になる。被積分関数の $f$ は `ElasticFactorWithTail`（静的設定 `AbsorptionElasticModel`、既定 `TemariHybrid`、260922Cl）で、Peng → Mott–Bethe テールや Temari の $f_e$ を使う（実部の `es.Factor` には影響しない）。

等方のときは、サイト内の等価位置で $M$ が同じなので $f'$ を 1 回だけ計算して使い回す。非等方・$g\neq0$ では等価位置ごとに数値積分をやり直す（L3902 の条件）。

### 2.5 キャッシュ `uDictionary`

`ConcurrentDictionary<(int Key1, int Key2), (Complex, Complex)>`。キーは `(compose(g.Index), h is null ? int.MaxValue : compose(h.Index))` で、**kV・inner/outer・nTheta/nPhi・ElasticIonModel・AbsorptionElasticModel を含まない**。`compose` は各成分を +255 して 10 ビットずつ詰める（$h$ は上位、範囲は概ね ±255）。`kV > 0` のときだけ格納する。

このため、電圧・検出器角・モデルを変える経路は呼び出し側で `uDictionary.Clear()` / `ClearUCache()` する約束になっている: `cbed_DoWork`（L452, 462）、`GetPrecessionElectronDiffraction`（L2360）、`StemDoWork`（L2523, 2535, 2877）は無条件に Clear（261014Cl）。`GetDifractedBeamAmpriltudes` は電圧が変わったときだけ Clear（L2282）。`MeanAbsorptionLengthNm` はエネルギーごとに Clear し、最後にも Clear する。→ 実装上の疑問点 Q1。

割り当て削減用のオーバーロード（`getU(kV, index, vec)`、STEM 用 `getU(kV, index, vq, vi, vj, inner, outer, …)`）はキャッシュヒット時に `Beam` を作らずに返すだけで、値は本体と同じ。

### 2.6 平均吸収長 `MeanAbsorptionLengthNm`

$$
\lambda_\mathrm{abs}(E) = \frac{k_\mathrm{vac}}{2\pi\,U'_0}\quad[\text{nm}]
$$

（強度の減衰長。$U'_0 \le 0$ なら $+\infty$。260923Cl。EBSD の深さ分布で使う。）$U'$ だけの平面波 $e^{2\pi i(K+i\kappa)z}$ で $(K+i\kappa)^2 = K^2 + iU'$ から $\kappa \approx U'/2K$、強度は $e^{-2\pi U' z/K}$ となることと一致する。$K$ の代わりに真空の $k$ を使う（差は $U_0/2K^2 \sim 10^{-4}$ とコメント）。

---

## 3. 行列

### 3.1 ポテンシャル行列 `getPotentialMatrix`（L4084–4136）

column-major の 1 次元配列 `potentialMatrix[row + col*dim]`:

$$
\mathcal U_{gh} = \begin{cases} i\,U'_0 & g = h \\ \hat U_{\mathbf g-\mathbf h} + i\,U'_{\mathbf g-\mathbf h} & g\neq h\end{cases}
\qquad(\texttt{getU(b[col] - b[row])}、\ g=\text{row},\ h=\text{col})
$$

- 対角には $U'_0$ の吸収だけを置き、$\hat U_0$（平均内部ポテンシャル）は入れない。$\hat U_0$ は $\mathbf K$ の屈折（§4）に回してある。
- 非局所版（`getPotentialMatrix(b, nonLocalInner, nonLocalOuter)`）は `getU(AccVoltage, b[col], b[row], inner, outer)` を使う。実部は局所版と同じ、虚部が $U'(\mathbf g,\mathbf h)$。

### 3.2 固有値問題 `getEigenMatrix`（L4147–4210）

Bloch 波を

$$
\psi(\mathbf r) = \sum_j \alpha_j \sum_\mathbf g C_\mathbf g^{(j)} \exp\!\big[2\pi i(\mathbf K + \mathbf g + \gamma^{(j)}\mathbf n)\cdot\mathbf r\big]
$$

と置き、$|\mathbf K|^2 = k^2 + U_0$ で $\gamma^2$ を落とした（前方散乱近似の）線形化 Bethe 方程式

$$
P_\mathbf g\,\gamma\,C_\mathbf g = Q_\mathbf g\,C_\mathbf g + \sum_\mathbf h \mathcal U_{\mathbf g\mathbf h}\,C_\mathbf h
$$

を解く。行列は

$$
A = D_P^{-1}\big(\mathrm{diag}(Q) + \mathcal U\big),\qquad A_{gh} = \frac{\mathcal U_{gh} + Q_g\delta_{gh}}{P_g}
$$

（**行 g を $P_g$ で割る**）。260829Cl でこの行割りに直された: 2021 年の 1 次元配列化で列割り $\mathcal U D_P^{-1}$ に退行しており、固有値は相似変換で不変だが固有ベクトルが $P_gC_g$ になって回折振幅に $P_g/P_0$ が乗っていた（無吸収での保存量が $\sum P_g|\psi_g|^2$ ではなく $\sum|\psi_g|^2/P_g$ になっていた）。検証ツール `tools/BetheDivCheck`、公開マニュアル appendix a3-bloch-wave と同一、とコメントにある。

$Q, P$（`Beam.Q`, `Beam.P`、L4639–4645）:

$$
Q_\mathbf g = |\mathbf K|^2 - |\mathbf K + \mathbf g|^2,\qquad P_\mathbf g = 2\,\mathbf n\cdot(\mathbf K + \mathbf g)
$$

（`Beam.Q` の summary「$-g(2k_0+g)$。$2k_0 S$ の励起誤差とわずかに定義が違う」、`P` の summary「大塚さんの資料参考」。）表面法線方向に測った励起誤差は

$$
S_\mathbf g = \sqrt{P_\mathbf g^2/4 + Q_\mathbf g} - P_\mathbf g/2
$$

（`Beam.S`）。$|\mathbf K + \mathbf g + S\mathbf n| = |\mathbf K|$ の解で、Ewald 球の内側（$Q>0$）で正。

部分集合版 `getEigenMatrix(bLen, beams, gMap, baseLen, basePotential, ref eigenMatrix)`（260602Cl）は、事前に作ったポテンシャル行列から `gMap` で抜き出して同じ $A$ を組む（EBSD 用。非局所吸収が黙って局所形式に戻る問題を直すため）。

### 3.3 解と境界条件

入射面で $\psi_\mathbf g(0) = \delta_{\mathbf g\mathbf 0}$ とし、$\boldsymbol\alpha = C^{-1}\boldsymbol\psi_0$（`LU().Solve` または逆行列）。深さ $t$（法線方向）で

$$
\psi_\mathbf g(t) = e^{\pi i P_\mathbf g t}\sum_j C_\mathbf g^{(j)}\,e^{2\pi i\gamma^{(j)} t}\,\alpha_j
$$

先頭の $e^{\pi iP_\mathbf g t} = e^{2\pi i\,\mathbf n\cdot(\mathbf K+\mathbf g)t}$ は「出射面での境界条件を考慮した位相」（GetDifractedBeamAmpriltudes: 20230827、CBED: 20220803）。経路ごとの違い:

| 経路 | 出射位相 |
|---|---|
| `GetDifractedBeamAmpriltudes` | $e^{\pi i P_g t}$ |
| `cbed_DoWork` | $e^{\pi i (P_g - 2k_\mathrm{vac} n_z) t}$（全ビーム共通の位相 $e^{-2\pi i k n_z t}$ が余分。強度には効かない） |
| `StemDoWork` | $e^{2\pi i (P_g/2) t}$（同じ式）。入射は $\psi_0 = -1$（符号は全体位相） |
| `GetPrecessionElectronDiffraction` | 掛けない（強度だけ使う） |

$A$ は吸収があると非エルミートなので、MathNet の `Evd(Symmetricity.Asymmetric)` または native の Eigen を使う。Bloch 波の規格化（$C$ の列の規格化）は不要で、$\alpha = C^{-1}\psi_0$ が吸収する。

---

## 4. 入射波数ベクトル `getVecK0`（L4658–4679）

結晶内波数ベクトルは表面に平行な成分を保存し、法線成分だけを変える:

$$
\mathbf K = \mathbf k_\mathrm{vac} + x\,\mathbf n,\quad |\mathbf K|^2 = k_\mathrm{vac}^2 + \hat U_0
\ \Rightarrow\ x = \sqrt{b^2 + \hat U_0} - b,\ \ b = \mathbf n\cdot\mathbf k_\mathrm{vac}
$$

$\mathbf k_\mathrm{vac}$ は `beamDirection` が null なら $(0,0,-k)$、それ以外は $k\,\hat{\mathbf b}$。$\hat U_0$ = `getU(voltage).Real.Real`（平均内部ポテンシャル、吸収なし）。`UseLocalSurfacePerBeamDirection = true`（EBSD マスターパターン用）のときは $\mathbf n = -\hat{\mathbf b}$（`GetEbsdSurfaceNormal`）。

---

## 5. ビームの選択 `Find_gVectors`（L4352–4552）

2 段構成（260318Cl）。

### 5.1 Phase 1: BFS キャッシュ `Build_gCache`

- 格子型ごとの隣接方向集合（F: ⟨111⟩ 8 本、I: ⟨110⟩ 12 本、A/B/C、R(H 設定)、六方/三方、P: ⟨100⟩ 6 本）で原点から BFS。**261030Cl** で `FrozenSet` → 配列に変更（I 格子の列挙順がプロセスごとに変わり、結果が再現しなかった）。
- 優先度付きキュー（キー $|\mathbf g|$）から、最小 $|\mathbf g|$ + 3·shift 未満を一巡ずつ展開（shift = 方向ベクトルの最大長/2）。点数 `cacheLimit = max(64N, 50,000)` に達するまで。
- 結果は `(key, gX, gY, gZ, |g|)`（実験室系、nm⁻¹）と、覆い切った半径 `CoverRadius = (前線の最小|g|) − 6·shift` を不変スナップショット `GCacheSnapshot` にまとめ、ロック下で 1 回だけ差し替える（261004Cl: PED の並列呼び出しで方位とキャッシュの不整合が起きたための修正）。同じ方位で長さが足りれば使い回す。

### 5.2 Phase 2: Ewald 殻の判定と順位付け

$\mathbf v = \mathbf K + \mathbf g$、$q = |\mathbf K|^2 - |\mathbf v|^2$ として

- 採用条件: $|q| < \text{maxQ} = \big||\mathbf K|^2 - (|\mathbf K| + \text{shift})^2\big|$ かつ $\mathbf n\cdot\mathbf v > 0$（前方へ進むもの）。
- キャッシュ先頭（概ね $|\mathbf g|$ の小さい順）から最大 `limit = 8N` 個集める（既定。診断用 `BeamCandidateMultiplier` ≠ 0 では $N\times$倍率、足りなければキャッシュを広げて探し直す）。
- 評価値 **Rating = $|\mathbf g|\,Q^2$**（`Beam.Rating`）の小さい $N+1$ 個を QuickSelect で選び、Rating 昇順にソート。$\mathbf 0$ は Rating 0 で先頭に来る（各ソルバは index 0 を入射ビームとして扱う）。
- **接平面内座標の重複除去**: 表面接平面の基底（`GetSurfaceTangentialAxes`）への射影を $10^{-6}$ nm⁻¹ で丸めたキーが同じビームは、Rating の小さい方だけ残す（260321Ch）。
- **末尾の同点組を落とす**（`BeamBoundaryTieMode.DropGroup`、既定）: $N+1$ 本目は番兵。末尾から Rating の差が $10^{-6}$ を超える最初の位置 $n$ を探し `beams[..n]` を返す。したがって結果は **N 本以下**で、対称等価な組が境界をまたぐときは組ごと外れる。`IncludeGroup`（診断）は組を丸ごと入れる。
- 各ビームに `getU(AccVoltage, hkl, g)` を格納（`Beam.Ureal/Uimag`）。

注意（コメントより）: 既定の候補上限 8N は、N ≳ 250 ではキャッシュの端で先に尽きることがあり（Ewald 殻の候補はキャッシュ点数の 2/3 乗でしか増えない）、候補範囲はキャッシュ半径で決まる（MnAl で |g| ≈ 80 nm⁻¹）。尽きた回数は `BeamCandidateExhaustedCount`、直近の状況は `LastBeamCandidateScan`（ThreadStatic）。

### 5.3 `reset_gVectors`

同じビーム集合（hkl・g・U）のまま、新しい $\mathbf K$（入射方向）と表面法線に対して $Q, P$ だけを再計算する。CBED・STEM の方向ループで使う。ビーム集合は中心方向で 1 回だけ選ぶ。

---

## 6. ソルバ

`enum Solver { Eigen_MKL, Eigen_Eigen, MtxExp_MKL, MtxExp_Eigen, Auto }`

| 名前 | 中身 |
|---|---|
| `Eigen_Eigen` | native（`Crystallography.Native.dll`、C++ Eigen）で EVD → $\psi(t)$。`NativeWrapper.CBEDSolver_Eigen` / `CBEDSolver2`（STEM 用、固有値・固有ベクトル・α も返す） |
| `Eigen_MKL` | MathNet `Evd(Asymmetric)`。MKL provider があれば MKL、無ければ managed（表示は `Eigen_Managed`） |
| `MtxExp_Eigen` | native で $\exp(2\pi i\,\Delta t\,A)$ を作って逐次掛ける（等間隔厚みのみ） |
| `MtxExp_MKL` | MathNet `Matrix.Exponential()`。厚みの区間幅が変わるたびに作り直す（260711Cl。以前は開始値≠間隔や不等間隔で誤った厚みを計算していた） |

- CBED の `Auto`: Eigen native があれば `MtxExp_Eigen`（スレッド = 論理コア数）、無ければ `Eigen_MKL`（MKL 有りで コア数/4、無しでコア数）。
- MtxExp 系は等間隔厚み前提（`NativeWrapper.IsUniformThickness`、相対 1e-12）で、そうでなければ EVD 系へ振り替える（260711Cl、native 側にも防御）。
- STEM は下流で固有値・固有ベクトル・α を使うので MtxExp 系は EVD に正規化。`Auto` は native 有りなら `Eigen_Eigen`、ただし MKL 有効かつ bLen ≥ `MklStemBlenThreshold`（既定 500、環境変数 `RECIPRO_STEM_MKL_BLEN`）なら `Eigen_MKL`。
- `GetDifractedBeamAmpriltudes` と PED は「native 有り かつ N < 400」で native `EigenSolver` + `Inverse`、それ以外は MathNet EVD + `Inverse()`。
- native ライブラリのソースはこのリポジトリには無い（`NativeWrapper.cs` は P/Invoke 宣言のみ）。

---

## 7. CBED / LACBED `RunCBED` → `cbed_DoWork`（L430–735）

入力: `maxNumOfBloch`（N）、電圧 [kV]、基準方位 `rotation`、厚み配列 [nm]、入射方向配列 `beamDirections`（単位ベクトル、正方格子に並べたもの）、`LACBED`、ソルバ、スレッド数。

1. `uDictionary.Clear()`、$\hat U_0$、中心方向の $\mathbf K$、`Find_gVectors` でビーム集合を 1 回だけ決める。ポテンシャル行列も 1 回だけ組む。
2. `beamDirections` を幅 $w = \sqrt{\text{Length}}$ の正方格子とみなし、内接円内（ピクセル中心で判定）の方向だけ計算。
3. 各方向: `getVecK0(kvac, u0, dir)` → `reset_gVectors` → `getEigenMatrix` → ソルバ → 出射位相（§3.3）。
4. LACBED: 透過波（index 0）の振幅だけを `Disks[t][0]` に。
5. CBED: ディスクごとに振幅を並べ、検出器座標 $\mathbf p = (\,(k\hat{\mathbf b}+\mathbf g)_x/(\ldots)_z,\ (\ldots)_y/(\ldots)_z\,)$（カメラ長 1）を計算。ディスクが重なる画素では、**他ディスクの強度 $|\psi|^2$ を非干渉に足し**、`Amplitudes = √(ΣI)`（実数）にする。重なり判定は外接矩形の交差グラフで近傍だけ見る（260610Cl）。最近接画素の探索 `getIndex` はディスク縁の画素を −1（無効）にする。

---

## 8. 回折振幅 `GetDifractedBeamAmpriltudes`（L2280–2341）

単一方向（`rotation`、入射 $-Z$）・単一厚み。N・電圧・方位が前回と同じなら固有値分解を使い回し、厚みだけ変える（HRTEM の厚みループが利用）。電圧が変わったときだけ `uDictionary` を Clear。出射位相 $e^{\pi i P_g t}$ を掛けた $\psi_g$ を `Beam.Psi` に入れて返す。

---

## 9. 歳差電子回折 PED `GetPrecessionElectronDiffraction`（L2355–2462）

- 半角 `semiangle` [rad]、分割数 `step`。$k$ 番目の入射方向は $R(\hat{\mathbf u}_k, \alpha)\,(0,0,-1)$、$\hat{\mathbf u}_k = (\cos 2\pi k/\text{step}, \sin 2\pi k/\text{step}, 0)$。
- 方向ごとに `Find_gVectors`（同じ BFS キャッシュを共有）→ 固有値問題。N・電圧・方位・半角・step が同じなら固有値分解を使い回す。
- 強度は方向について**非干渉平均** $I_\mathbf g = \frac1{\text{step}}\sum_k |\psi_\mathbf g^{(k)}(t)|^2$（ある方向でビーム集合に入らなかった g はその方向の寄与 0）。
- 最後に基準方位で $Q, P$ を再計算し、`Psi = √I`（実数）、Rating 昇順に並べて返す。
- 出射位相は掛けない（強度だけ使うため不要）。

---

## 10. 主な入力と既定値（このページの範囲）

| パラメータ | 既定 | 場所 |
|---|---|---|
| `Surface` | (0,0,−1) | プロパティ |
| `ElasticIonModel` | `Neutral` | static |
| `AbsorptionElasticModel` | `TemariHybrid` | static |
| `getU` の annular 求積 `nTheta`, `nPhi` | 60, 20 | 引数既定 |
| `FactorImaginary` の θ 区切り / 点数 | s の等比 [0, 0.3, 1, 3, 9, 27, 81, 243, 729, ∞] nm⁻¹ / θ 12 点・φ 8 点（折り畳み） | `AtomStatic.ES`（[scattering-factors.md](scattering-factors.md)） |
| BFS キャッシュ点数 | max(64N, 50,000) | `Find_gVectors` |
| 候補上限 | 8N（`BeamCandidateMultiplier = 0`） | 同 |
| 同点判定・接平面重複の許容 | 1e-6 | 同 |
| `BeamBoundaryTieMode` | `DropGroup` | static |
| `MklStemBlenThreshold` | 500 | static（環境変数） |
| `EbsdRepresentativeDirectionBlockSize` | 3 | EBSD 用 |
| `UseLocalSurfacePerBeamDirection` | false | EBSD マスター用 |

---

## 実装上の疑問点

**Q1. U キャッシュのキーが検出器角を含まず、STEM 後に弾性経路が環状 U′ を読む可能性**
`getU` のキーは `(compose(g), int.MaxValue)` で kV・inner/outer を含まない（L3865–3867, 3976–3978, 3997–3999）。`StemDoWork` は L2877 で Clear した後、`fillUq`（L2895）で **環状域 [detAngleInner, detAngleOuter] の U′** を通常の局所 U と同じキーで格納し、終了時に Clear しない。その後、同じインスタンス・同じ電圧で `GetDifractedBeamAmpriltudes` を呼ぶと（L2282 は電圧が変わったときしか Clear しない）、`Find_gVectors` 内の `getU`（L4466）がキャッシュを引き、吸収が全立体角 TDS ではなく ADF 検出器域の TDS になる。`GetHRTEMImage` は内部でこの関数を呼ぶ。CBED・PED・STEM・EBSD は冒頭で無条件 Clear するので影響しない。ReciPro が STEM と HRTEM/回折で同じ `BetheMethod` を共有しているかはこのリポジトリからは確認できない。

**Q2. STEM と HRTEM で時間的コヒーレンス包絡の係数が 2 倍違う**
HRTEM（L3735 準コヒーレント、L3763 TCC）は $\exp[-\tfrac12\pi^2\lambda^2\delta^2(q_1^2-q_2^2)^2]$、STEM の `Lenz`（L2741–2742）は $\exp[-\pi^2\lambda^2\delta^2(k_q^2-k^2)^2]$ で $\tfrac12$ が無い。$\delta$ を焦点ずれの rms（ガウス分布の標準偏差）とすると平均は前者になる。同じ「delta = Cc·ΔE/E」という入力が STEM では実効的に $\sqrt2$ 倍に効く。どちらかの定義（rms か 1/e 幅か）を明示する必要がある。詳細は [bloch-wave-imaging.md](bloch-wave-imaging.md)。

**Q3. HRTEM 準コヒーレント分岐のコメントの取り違え**（数値は正しい）
L3734 の β の項に「Ec 時間的インコヒーレンス」、L3735 の δ の項に「Es 空間的インコヒーレンス」とあるが逆。TCC 分岐（L3762–3763）は正しく書かれている。

**Q4. `MottBetheMonopoleCoefficient` の doc コメントの単位**（数値は正しい）
`AtomStatic.cs` L1210 は「s² = (sinθ/λ)² [Å⁻²]、f は [nm]」と書くが、Mott–Bethe の $f[\text{Å}] = 0.023934\,\Delta Z/s^2[\text{Å}^{-2}]$ を f[nm], s²[nm⁻²] に換算すると係数は 0.1×0.023934/0.01 = 0.23934 で、0.2393366 は **s² が nm⁻²** の係数。`getU`（L3893）は nm⁻² の s2 を渡しているので数値は正しく、`ElasticFactorWithTail` 付近のコメント（「s² が nm⁻²・f が nm の規約」）とも一致する。L1210 の記述だけが食い違う。

**Q5. 固有値分解の再利用条件が結晶・モデルの変更を見ない**
`GetDifractedBeamAmpriltudes`（L2292）と PED（L2370）は N・電圧・方位（PED は半角・step も）だけで再計算を判定し、原子・DW 因子・`Surface`・`ElasticIonModel`・`AbsorptionElasticModel` の変更を見ない。PED は冒頭で U キャッシュを Clear するが、固有値分解が再利用されれば新しい U は使われない。結晶を編集するたびに呼び出し側が `BetheMethod` を作り直す運用なら問題ない（要確認）。

---

## 要確認事項

1. `Find_gVectors` の接平面座標による重複除去（L4475–4493）は、表面法線方向の成分だけが違うビーム（HOLZ 列などで $\mathbf g$ が $\mathbf n$ に平行な差を持つ組）を 1 本に減らす。線形化 Bethe 法の出射条件として意図どおりか。
2. 非等方 DW のとき吸収の被積分関数 $1-e^{M(s^2-s_1^2-s_2^2)}$ に「g 方向の実効等方 B」を使っている（L3912 → `FactorImaginary(kV, s2, m)`）。TDS の中間散乱ベクトル $s_1, s_2$ は g 方向とは限らないので近似。許容範囲として扱ってよいか。
3. `getU` のコメント「非弾性散乱因子 答えは無次元」（L3918）と、`FactorImaginary` の戻り値（nm、実部と同じ単位で足し合わされる）の食い違い。summary 側の記述を直すべきか。
4. CBED の出射位相 $e^{\pi i(P - 2k n_z)t}$（L589）の $-2k n_z$ 項の意図（全ビーム共通の位相で強度には効かない）。`Disks[].RawAmplitudes` の位相を使う下流があるか。
5. `UseLocalSurfacePerBeamDirection = true` のまま CBED/STEM を走らせると、`getVecK0`（方向ごとの法線）と `reset_gVectors`（固定 `Surface`）で法線が食い違う。EBSD 用に立てたフラグが同じインスタンスに残る運用があるか。
6. 定数の年代が混在（m0, e0 は 1986 年値、h は 2006 年値、λ の係数は別の数値）。意図的か。
7. `Find_gVectors` で厳密に $Q = 0$（Rating = 0）のビームがあると $\mathbf 0$ と同点になり、不安定ソートで index 0 が入れ替わる可能性（浮動小数で厳密 0 はまず起きない）。
8. `getIndex`（CBED のディスク重なり）がディスク縁の画素を無効にするため、重なり領域の縁 1 画素で強度が足されない。仕様か。

（要確認 8 件、実装上の疑問点 5 件）
