本稿はコードから Claude が起こした下書き（2026-11-05）。著者未査読。

# EBSD の指数付け・方位探索・幾何較正

対象コミット: `9ffcf47`（main）。対象ファイル: `EBSD/EbsdBandDetector.cs`、`EBSD/EbsdIndexer.cs`、`EBSD/EbsdRadonIndexer.cs`、`EBSD/EbsdDictionaryIndexer.cs`、`EBSD/EbsdZoneAxisIndexer.cs`、`EBSD/EbsdOrientationSearch.cs`、`EBSD/EbsdGeometryCalibrator.cs`、`EBSD/EbsdPatternProjector.cs`（`EbsdPatternScorer`）。

方位は回転行列 $R$（結晶 → 試料系、$\mathbf v^{\rm s}=R\mathbf v^{\rm c}$）、幾何は `EbsdDetectorGeometry`（[ebsd-geometry.md](ebsd-geometry.md)）。反射の入力はどの経路も `Crystal.VectorOfG_KikuchiLine`（結晶直交系の $\mathbf g$、$d$ [nm]、`RelativeIntensity` = 運動学的強度）。

---

## 0. 全体像

| 経路 | 入力 | 評価関数 | 実装 |
|---|---|---|---|
| バンド検出 → pair-angle | 実測画像 | 割当の角度残差・幅尤度 | `EbsdBandDetector.Detect` → `EbsdIndexer.Index` |
| Radon テンプレート照合 | 実測画像 | Radon 証拠の SNR 正規化和 | `EbsdBandDetector.ComputeRadonMap` → `EbsdRadonIndexer.Index` |
| 辞書照合 | 実測画像 + マスターパターン | robust 前処理後の ZNCC | `EbsdDictionaryIndexer.Index` |
| 晶帯軸 | 手で拾ったバンド交点 | 角度残差 + 低指数の優遇 | `EbsdZoneAxisIndexer.Index` |
| 統合（GUI の Find） | 上記 + マスターパターン | 複合ランク（Radon z + ZNCC） | `EbsdOrientationSearch.Run` |
| 幾何較正 | 実測 + マスターパターン + 方位 | ZNCC | `EbsdGeometryCalibrator.Run` |

---

## 1. 共通: 反射カタログ

`EbsdIndexer.BuildCatalog` / `EbsdRadonIndexer.BuildCatalog`:

- $\hat{\mathbf g} = \mathbf g/|\mathbf g|$ の符号を辞書順で正に揃え（$z>0$、次に $y>0$、次に $x>0$）、$\hat{\mathbf g}\cdot\hat{\mathbf g}' > 0.99995$（約 0.6°）のものを 1 ノードにまとめる（$\pm\mathbf g$ と調和反射）。
- `EbsdIndexer`: ノードは反射の一覧・$d$ の一覧・最大強度を持ち、強度上位 `maxCatalogNodes`（既定 120）。
- `EbsdRadonIndexer`: ノードは代表（最強）反射と重み $w$。強度上位 `maxNodes`（`DefaultMaxNodes`）。重みは `weightExponent` = 0.5（既定）なら $w=\min(\sqrt I,\ 3\,\mathrm{median}\sqrt I)$、それ以外は $w=\mathrm{clip}((I/\mathrm{median}I)^{p},\,0.5,\,2)$。

---

## 2. バンド検出（`EbsdBandDetector.Detect`）

doc の手順: 長辺 256 px へ縮小 → ガウシアン背景除算 → 正規化 → Radon（sum/√N 正規化）→ 幅バンク butterfly 畳み込み → MAD 閾値 + 画像内線距離 NMS → θ/ρ 局所再探索 → 横断プロファイルからエッジ・幅。

| 段 | 内容 | 値 |
|---|---|---|
| 縮小 | 長辺 `WorkLongSide` | 256 px（位置の定数は `WS = 256/256` で換算） |
| 背景 | 広域ガウシアン（無効画素を重みから除く normalized convolution）で除算 | $\sigma = 0.20\times$短辺（260724Cl に 0.1 から変更） |
| 平滑・正規化 | ガウシアン → 平均 0・分散 1 → ±5σ クリップ | $\sigma = 0.8$ |
| Radon | $\theta\in[0,180°)$、$\rho$ は画像中心からの符号付き距離 [px]。線積分はバイリニア・1 px 刻みの和 / √N | $\Delta\theta = 0.5°$、$\Delta\rho = 1$ px |
| butterfly | $K_w(t) = G_\sigma(t) - \tfrac12G_\sigma(t-w/2) - \tfrac12G_\sigma(t+w/2)$、$\sigma=\max(0.75, 0.12w)$、平均 0・L2 ノルム 1。各 $(\theta,\rho)$ で $\lvert s\rvert$ 最大の幅の応答を**符号付き**で保持（明バンド正・暗バンド負） | $w\in\{3,5,7,10,14,20,28,38\}$ px |
| θ 平滑 | 中心と同符号の隣だけを足す 3 タップ（$\theta$ と $\theta+180°$ で $\rho$ 反転） | 1 パス |
| 閾値 | $\lvert\text{応答}\rvert > \mathrm{median} + 4.5\,\mathrm{MAD}$ かつ $\rho$ 方向の極大 | |
| NMS | 同符号: 画像内 RMS 距離 6 px。異符号かつほぼ平行（$\Delta\theta<10°$）: $\max(6,\ 0.75\max w)$（サイドローブ抑制）。異符号・近平行でスコア比 ≥ 0.7 は「明/暗ローブ対」としてタグ付け | |
| 局所再探索 | $\rho$ は $\pm2$ px を 0.2 px 刻み、$\theta$ も近傍で | |
| 重複除去 | 同符号 5 px、異符号近平行 $\max(5, 0.75\max w)$ | |
| 線検証ゲート | `CenterQuality` < (Score ≥ 17 なら 0.3、それ以外 0.5)、または（ローブ対でなく Score < 17 で `EdgeLinearity` ≥ 0.55）、または `WidthQuality` ≤ 0.05 のバンドを捨てる | |

出力 `EbsdBand`（native 画素座標）: 正規化した中心線 $A\,\mathrm{col}+B\,\mathrm{row}+C=0$（$A^2+B^2=1$）、中心線上のアンカー点、両縁の点（各アンカー ± 半幅）、幅、Score（MAD 単位）、CenterQuality（狭窓アンカー射影の中央値）、WidthQuality ほか。`maxBands` 既定 12。

`ComputeRadonMap` は同じ前処理〜θ 平滑までを共有し、離散検出をせずに $|\text{応答}|$ のマップ `EbsdRadonMap.Abs` と robust null 統計（$\mu_0$ = 中央値、$\sigma_0 = 1.4826\,\mathrm{MAD}$）を返す。

---

## 3. pair-angle 指数付け（`EbsdIndexer.Index`）

### 3.1 観測

各バンドについて面法線（lab）$\mathbf g^{\rm lab}$ = `LineToLabNormal`、試料系 $\hat{\mathbf m} = R_x(\tau)\mathbf g^{\rm lab}$（符号不定）。縁の点ごとに $\sin\theta_B = |\hat{\mathbf g}\cdot\mathbf P_e|/|\mathbf P_e|$ を求め、2 点以上かつ WidthQuality > 0.05 なら中央値 $\sin\theta_B$ と $\sigma = \max(1.4826\,\mathrm{MAD},\ 0.2\sin\theta_B)$。重み $w_o = 0.3 + 0.7\,\text{CenterQuality}$。

### 3.2 手順

1. 理論ペア表: カタログ方向の対の角 $\arccos|\hat{\mathbf g}_p\cdot\hat{\mathbf g}_q|\in(8°, 85°)$ を角度順に並べる。
2. 観測は重み上位 8 本。対 $(i,j)$ の角度 ± `pairTolDeg`（1.5°）に入る理論対について、対応の入れ替え 2 通り × 観測法線の符号 4 通りで Wahba 問題を解く。
3. **Wahba / Kabsch**（`SolveWahba`）: $H=\sum w\,\mathbf m\mathbf g^{\mathsf T}$、$H = U\Sigma V^{\mathsf T}$、

$$ R = U\,\mathrm{diag}\bigl(1,1,\det(UV^{\mathsf T})\bigr)\,V^{\mathsf T}\qquad(\sigma_2<10^{-10}\text{ なら不定}). $$

4. 採点（`ScoreAndAssign`）: 予測 $R\hat{\mathbf g}_n$ と観測 $\hat{\mathbf m}_o$ の角 $\psi = \arccos|\hat{\mathbf m}\cdot R\hat{\mathbf g}|$ が `assignTolDeg`（2.5°）未満の対を残差の昇順で**一対一**に貪欲割当。

$$ \text{Score} = \sum_{\rm assigned} w_o\,\Bigl(1-\frac{\psi}{\psi_{\rm tol}}\Bigr)\,f_w,\qquad
f_w = 0.4 + 0.6\exp\!\Bigl(-\tfrac12\min(r,4)^2\Bigr),\quad r = \min_m\frac{|\sin\theta_B - \lambda/(2d_m)|}{\sigma}. $$

   ($\sin\theta_B$ が無ければ $f_w=1$)。$\lambda$ は nm。割当は最も近い $d_m$ の調和反射の $(hkl)$。割当 3 本以上を候補にする。
5. 上位 60 件を「割当に合わせて観測法線の符号を揃えた重み付き Kabsch」で 2 回まで精密化（スコアが上がったときだけ更新）。
6. 方位差 3° 未満（対称考慮なし）の重複を除いて上位 `maxCandidates`（10）。`AngularRmsDeg` も返す。

---

## 4. Radon テンプレート照合（`EbsdRadonIndexer.Index`）

### 4.1 評価関数

方位 $R$（と幾何）の下で検出器と交わる予測バンド中心線 $k$ について、試料系 $R\hat{\mathbf g}_k$ → lab → `LabNormalToLine` と同じ係数 $(a,b,c)$ [mm] から $\theta_k=\operatorname{atan2}(b,a)$、$\rho_k = -c/(s\sqrt{a^2+b^2})\times\text{Scale}$（work px）。

$$ \text{Score}(R) = \frac{\sum_k w_k\,\psi\bigl(E(\theta_k,\rho_k)-\mu_0\bigr)}{\sqrt{\sum_k w_k^2}\;(\sigma_0\ \text{if cap}=0)},\qquad
\psi(e) = \begin{cases}\max(e,-\sigma_0) & \text{cap}=0\\ \text{cap}\cdot\tanh\bigl(z/\text{cap}\bigr)\ (z>0),\ \max(z,-1)\ (z\le0),\ z=e/\sigma_0 & \text{cap}>0\end{cases} $$

（`ScoreExactCore`）。$E$ は証拠マップのバイリニア値（マップ外の $\rho$ は 0）。近接した予測線（$\Delta\theta<2°$ かつ $\Delta\rho<5$ px、$\theta$ 循環と $\rho$ 反転を考慮）は強度順で先着の 1 本だけ得点する（二重得点の防止）。$(\sum w)^2/\sum w^2 < 4$ なら無効。

### 4.2 探索

- シード ①: 証拠マップの局所極大（$\mu_0+3\sigma_0$ 超、NMS）を擬似バンドとして pair-angle + Kabsch（§3）。コメント: グリッド + NM だけではリッジ（θ 幅 ~1°）に届かない例があったため。
- シード ②: SO(3) 粗格子 $R = R_{\rm sphere}(i)R_z(\phi)$（刻み `coarseStepDeg` = 3°）を、強度上位 `coarseNodes` = 20 ノードと**膨張マップ**（$\theta\pm$`dilTheta` 行 × $\rho\pm$`dilRho` 列の窓 max、最近傍参照）で採点。膨張幅は半刻みの回転で予測線が動く量: `dilRho` $=\max(3,\lceil 0.7\sqrt{W_w^2+H_w^2}\,\delta/2\rceil)$、`dilTheta` $=\max(2,\lceil(\delta/2)/\Delta\theta\rceil)$。上位 `maxGridSeeds`（52）。
- 精密化: 厳密スコアを Nelder–Mead で最大化。変数は $[\omega_x,\omega_y,\omega_z]$（試料系左摂動、deg）に、幾何探索が有効なら $(\Delta Z,\Delta X,\Delta Y)$ [mm] を足す。
- 等価の縮約: $Q=R_1^{\mathsf T}R_2$ がカタログ方向集合を（± 込みで）保存するならパターン上区別できないので 1 代表に縮約。幾何が違う候補は別に残す（260922Cl）。
- 幾何探索（260922Cl、既定は無効）: `zSearchRangeMm`（DetZ を格子で探す。刻みの既定は膨張幅を mm に戻した量）、`xSearchRangeMm`、`xyRefineRangeMm`（精密化だけで X・Y を動かす）、`coarseGeometryToleranceMm`（膨張幅に足す）。

---

## 5. 辞書照合（`EbsdDictionaryIndexer.Index`）

### 5.1 前処理と類似度（`EbsdPatternScorer`）

- `Downsample`: 長辺 `targetLongSide` への box 縮小（下限 8 px）。
- `RobustPreprocess`（実測とシミュレーションの**両方**に掛ける）: 広域背景 $\sigma = 0.1\times$短辺のガウシアン、$v = \ln\bigl(\max(x, 0.01f)/\max(\mathrm{bg}, f)\bigr)$（$f = \max(10^{-10}, 0.05\bar x)$）→ DoG（$\sigma_1=1.5$、$\sigma_2=6$）→ 標準化 → $\tanh(z/3)$ → 再標準化。doc: 動力学の単一スライスは晶帯軸の明点が分散を支配する heavy-tailed な分布で、生の ZNCC では正解が偽方位に負けるため。
- `RobustPreprocessFast`: ガウシアンを box 3 連（$r = \mathrm{round}((\sqrt{4\sigma^2+1}-1)/2)$、合成分散 $\sigma^2$）で近似し、SIMD（`Vector.Log`）・作業領域の再利用・逐次実行。
- ZNCC: 参照 $\mathbf a$ は平均 0・分散 1（母分散）、$\mathrm{ZNCC} = \sum a_i(b_i-\bar b)/(n\,\sigma_b)$、$\sigma_b<10^{-12}$ なら 0。`ZnccSimd` は加算順だけ違う（較正専用）。

### 5.2 3 段構成

1. **粗段**: Fibonacci 球 $n_s = \max(64,\lfloor4\pi/\delta^2\rfloor)$ 点 × 面内 $n_\phi=\max(16,\lfloor360/\delta\rfloor)$ 回転（$\delta$ = `coarseStepDeg`、GUI 経路は 3°）。正方格子では面内回転分解投影（[ebsd-master-pattern.md](ebsd-master-pattern.md) §4.2）、六方格子は通常の `Project`。GUI は常に `thoroughCoarse = true`（96 px・`RobustPreprocessFast`）。`false`（48 px・軽量前処理: 平均で正規化した log-ratio → box DoG 5×5 − 21×21 → tanh → 正規化）は検証ハーネスの対照用と doc にある。方位差 $1.5\delta$ 未満を捨てる NMS で上位 `coarseKeep`（64）盆地。
2. **中段**: 上位を完全 robust 前処理（96 px）で再採点して上位 `refineKeep`（12）。doc の注: `thoroughCoarse = true` では参照・解像度・投影器が粗段と同じなので「同一パイプでの再採点」になる。
3. **精段**: Nelder–Mead（試料系左摂動。初期刻み $\delta/2$ で 120 評価 → 0.5° で 80 評価）→ 方位差 NMS → 上位 `maxCandidates`。返す Score は NaN（Radon z の付与と複合ランクは統合層が行う）。

対称による削減（`properSymmetries`）: 同値類 $\{RS\}$ のうち $\operatorname{tr}$ 最大（= 回転角最小）の代表だけを評価する。$R=r_0R_z$ では $\operatorname{tr}(R_z M)$ を $M=r_0$、$M=Sr_0$ について閉形式で比べる。doc の訂正（260725Cl）: tie 帯 $10^{-12}$ は厳密な対称不変点しか救わず、被覆は「軌道の相手の近傍に別の格子点がある」機構で保たれている。そのため `EbsdOrientationSearch` は FZ 除外を proper 回転 1 個（C2）に限定している。

---

## 6. 晶帯軸からの指数付け（`EbsdZoneAxisIndexer.Index`）

- カタログ: 反射の整数外積から既約 $[uvw]$ を作り、`maxIndex`（4）以下に限る。重み = その晶帯に属する反射の相対強度の和。上位 `maxNodes`（300）を ± 両方向で。doc: 高指数を許すとどんな点でも説明できてしまう（Botallackite で $[10\,\bar2\,5]$ が割り当たった）。
- 観測: 拾った点 $(\mathrm{col},\mathrm{row})$ → **射出方向** $\hat{\mathbf o} = -$`PixelToSampleDirection`。
- 種: 2 点の対の角度がカタログの対と `toleranceDeg`（2°）以内で合うものを総当たりし Kabsch。手入力の指数があればその点の割当を $\pm[uvw]$ に限定、2 点以上手入力なら対をそれに絞る。
- 採点（`Evaluate`）: 各点を最も近い（内積最大の）回転後カタログ方向に割り当て、

$$ \text{Score} = \sum_{\rm assigned}\Bigl[\Bigl(1-\frac{\psi}{\psi_{\rm tol}}\Bigr) + 0.5\,\frac{W_c}{W_{\max}}\Bigr] + 2\,n_{\rm assigned}. $$

  異なる点が同じ晶帯軸に割り当たる解、手入力の点が割り当たらない解は棄却。`RmsDeg` と `RmsPx`（予測位置 `SampleDirectionToPixel` と拾った位置の画素 RMS）を返す。
- 任意（`refineGeometry`）: 割当を固定して方位 3 + 検出器中心 3（$\Delta X_d,\Delta Y_d,\Delta Z_d$）を Nelder–Mead で同時に最適化し、画素残差を最小化。4 点未満なら幾何は動かさない。

---

## 7. 統合探索（`EbsdOrientationSearch.Run`）

1. エンジン: `useDictionary` で辞書照合、そうでなければ Radon 照合。辞書の候補には Radon z を後付けする。
2. 複合ランク（`context` があるとき）: 各候補を投影し robust ZNCC を計算、**候補集合内で標準化**して

$$ \text{combo} = z_{\rm Radon} + 0.5\cdot\mathrm{clip}(z_{\rm ZNCC},\,\pm2). $$

   Radon の採点は証拠飽和 cap = `SaturateCap`（8）を使う（doc: 単独では劣化するので複合とセットでのみ）。
3. 精密化: 複合トップ 1 件だけ ZNCC で精密化（Nelder–Mead、初期刻み 0.25°、120 評価）。Radon z の低下が `GuardMaxZDrop`（0.2）を超えたら採用しない（誤収束ガード）。
4. 仕上げ: 最終方位を、幾何較正の最終段と**同じ目的関数・同じ刻み**（`EbsdGeometryCalibrator.CreateFinalScorer`。Nelder–Mead を初期刻み 0.7° で 150 評価 → `OrientationPolishStepDeg` = 0.1° で 100 評価）で仕上げ直す。doc: 目的関数が違うと Find と Calibrate の繰り返しで方位が約 1° 往復した。仕上げにも同じガード。
5. 幾何も探す構成（`EbsdSearchOptions.GeometryRobust`、260922Cl）: `ZSearchRangeMm = 12`、`XSearchRangeMm = 3`（刻み 3）、`CalibrateTopCandidates = 5`（上位 5 件を開始点 5 個・1 段で軽量較正し、較正後の ZNCC で並べ直す）、`CalibrateIncludePlain = true`（幾何固定の素の探索の上位も較正に回す）、`DictionaryZHypotheses = 1`、`DictionaryMergeRadon = true`、`SaturateCap = 4`、`RadonMaxGridSeeds = 150`。候補を較正したときは §7.3・7.4 の微調整を行わない（doc: 較正で 0.03° まで合った方位を 0.66° ずらしていた）。

`EbsdMatchingContext`: 幾何、マスターパターン、照合に使う平面（MC 重み合成または単一スライス）、縮小済みの参照、ラスター寸法、現在の方位、フル解像度の表示値 `DisplayReference` と、シミュレーション側に掛ける高域通過の半値幅 `SimFlattenFwhmPx`（実測画像 px）。

`EbsdSearchOptions.Default` の主な値: Radon 粗刻み 3°、粗ノード 20、`RadonMaxGridSeeds` 52、辞書粗刻み 3°、`DictionaryCoarseKeep` 64、`DictionaryRefineKeep` 12、`DictionaryThoroughCoarse = true`、`ZnccCoef` 0.5、`ZnccClip` 2、`GuardMaxZDrop` 0.2、幾何探索 0。

---

## 8. 幾何較正（`EbsdGeometryCalibrator.Run`）

### 8.1 変数と拘束

- 変数: PC の物理面内座標 $(u_F, v_F)$ [mm]、$\ln DD$、方位 3（試料系左摂動）。**$\delta$（DetTilt）は固定**（doc: 単一パターンでは DetTilt と方位の $x$ 回転がゲージ自由度になる）。
- 幾何への戻し: `FromPatternCenter(u_F, v_F, e^{\ln DD}, δ)`（[ebsd-geometry.md](ebsd-geometry.md) §3.4、疑問点 G-I1）。
- ソフト境界: 開始点から $|\Delta u_F| > 0.25\,W_{\rm phys}$、$|\Delta v_F| > 0.25\,H_{\rm phys}$、$|\Delta\ln DD| > 0.35$ のいずれかなら罰則 $10 + |\Delta u|/W + |\Delta v|/H + |\Delta\ln DD|$ を返す。

### 8.2 目的関数

$-\mathrm{ZNCC}$。参照は `DisplayReference`（表示中の値、ユーザーの「背景を平坦化」設定が反映されたもの）を各段の長辺へ box 縮小して正規化。シミュレーション側は `EbsdPatternProjector` で投影し、`SimFlattenFwhmPx` > 0 なら同じ半値幅（段の縮小率で換算）の高域通過（box 3 連で近似したガウシアン背景を差し引く `SubtractBoxBackground`）を掛ける。

### 8.3 手順

- 解像度 3 段: 長辺 160 → 480 → フル解像度（`StageLongSides` = [160, 480, 0]）。段の収束条件 `(Tol, XTol)` = (1e-6, 1e-2)、(1e-7, 3e-3)、(1e-9, 3e-4)。評価上限（方位 / 幾何 / 6 変数）= (150, 120, 600)、(200, 180, 900)、(300, 250, 1500)。
- 1 開始点の較正: 交互法（① 幾何固定で方位〔初期刻み 0.7°〕、② 方位固定で幾何〔刻み = 検出器幅・高さの 1 %、$\ln DD$ 0.02〕）を最大 `MaxRounds` = 20 ラウンド、1 ラウンドの改善が `ZnccTolerance` = 1e-4 未満で打ち切り → 方位の仕上げ → 6 変数同時の Nelder–Mead（最終段だけ `NelderMeadRestart` で刻みを 1/4 にして最大 3 回張り直す）。
- 多点開始: 段 0 で `StartCount` = 40 点。オフセットは $[-1,1]^3$ の Halton 列（[0] は現在の幾何）× (PC: 検出器幅・高さの 8 %、$\ln DD$: 0.08)。乱数を使わないので同じ入力なら同じ結果。段 1 へは上位 $\max(6,$ CPU 数$)$、段 2 へは 1 点を引き継ぐ。
- 診断値: 最良の開始点番号、最良と最悪の ZNCC 差、最良から 1e-3 以内の開始点数、その集団の PC・DD の広がり（`FlatU/V/Dd`）。

---

## 9. 実装上の疑問点

### I-I1 `SaturateCap` の doc と値が合わない

`EbsdOrientationSearch.SaturateCap` の doc は「係数スイープで決定（プラトー 0.4-1.0 の中央寄り）」だが、値は 8（`GeometryRobust` は 4）。プラトーの値域が別の正規化（例: $\sigma_0$ 単位でない ψ）を指しているのか、doc が古いのかが分からない。

### I-I2 pair-angle 指数付けの重複除去は対称を考慮しない

`EbsdIndexer.Index` の最終段は `MisorientationDeg(c.Rotation, current.Rotation)`（`properSymmetries` なし）で 3° 未満を重複とみなす。対称等価な方位は別の候補として残り、上位 `maxCandidates` 枠を消費しうる（doc に「対称等価は含まれ得る」と明記されているので仕様とも読める）。

### I-I3 幾何較正の開始点が `FromPatternCenter` の符号規約に依存する（**確認済み: UI の範囲では問題なし**）

§8.1 の通り、較正は $(u_F, v_F, \ln DD)$ から `FromPatternCenter` で幾何を作る。$\mathbf n\cdot\mathbf C<0$ の配置なら開始点 [0] が現在の幾何と一致しないが、`FormEBSD` の運用範囲（試料傾斜 ≤ 0、既定 −70°、`DetY` = −35 mm）では常に $\mathbf n\cdot\mathbf C>0$ で、開始点 [0] は現在の幾何に一致する（[ebsd-geometry.md](ebsd-geometry.md) の G-I1）。

---

## 10. 要確認事項

- **I-1** 反射カタログの強度 `RelativeIntensity` は運動学的強度（`VectorOfG_KikuchiLine` の定義）。Radon・pair-angle とも動力学の強度とは無関係に重みを決めている。重元素・低対称の結晶での妥当性。
- **I-2** バンド検出の定数（背景 $\sigma$ 0.2×短辺、閾値 4.5 MAD、NMS 6 px、線検証 0.5/0.3/17/0.55/0.05）は特定の実測画像群での定量評価で決めたとコメントにある。対象画像の条件（検出器・解像度・ビニング）を文書化するか。
- **I-3** 幅尤度の $\lambda/(2d)$ は運動学的な Bragg 角。バンド縁の位置（Kikuchi 線の位置）と、動力学的なバンドの見かけの幅との差（[ebsd-pattern-composition.md](ebsd-pattern-composition.md) §4.2 の $A(E)$ の議論とも関係）。
- **I-4** 辞書照合は MC 重み合成の 1 枚（または単一スライス）を投影する。検出器上の位置による重みの違い（[ebsd-pattern-composition.md](ebsd-pattern-composition.md) §5 の近似）が照合精度に与える影響。
- **I-5** 幾何較正のソフト境界・多点開始の幅・段ごとの収束条件は作者の実機での経験値（コメント）。他の検出器での妥当性。
- **I-6** 晶帯軸指数付けの採点式（$+0.5W/W_{\max}$、$+2n$）の係数の根拠。
- **I-7** Radon の等価判定（$Q$ がカタログ方向集合を保存）は、カタログが ± を同一視しているため、結晶の点群より大きい群（Laue 群 + 運動学的に区別できない操作）で縮約しうる。極性結晶で意図通りか。
