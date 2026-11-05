本稿はコードから Claude が起こした下書き（2026-11-05）。著者未査読。

# EBSD マスターパターン: 格納形式・投影・ファイル

対象コミット: `9ffcf47`（main）。対象ファイル: `EBSD/EBSD.cs`（`EBSD`、`MasterPattern`）、`EBSD/MasterPatternFile.cs`、`EBSD/EbsdPatternProjector.cs`、`EBSD/EbsdPatternComposer.cs`（参照テーブル部分）。

動力学計算そのもの（Bethe 法による各方向・各エネルギー・各厚さの強度）は回折担当の文書（`docs/methods/` の Bethe 法・Kikuchi の章）を参照。ここでは、その出力をどう並べて持ち、どう引くかだけを書く。座標系と回転の向きは [ebsd-geometry.md](ebsd-geometry.md)。

---

## 1. 構築（`EBSD.RunMasterPatternBuild`）

### 1.1 入力（`EBSD.MasterPatternBuildRequest`）

| 引数 | 意味 | 既定 |
|---|---|---|
| `Crystal` | 結晶 | 必須 |
| `MaxNumOfBloch` | Bethe 計算のブロッホ波の最大数 | 必須 |
| `Energies` | エネルギー格子 [keV] | 必須（MC 分布側の doc では降順） |
| `Depths` | 深さ格子 [nm]（出射方向に沿った経路長、単調増加。§2.3） | 必須 |
| `GridSize` | 1 辺の分割数（六方格子では $2N+1$） | 必須 |
| `Solver` | `BetheMethod.Solver` | `Auto` |
| `Thread` | スレッド数 | 1 |
| `UseNonLocalAbsorption` | 非局所吸収 | false |
| `IncludeTDSBackground` | TDS 背景（非局所源） | false |
| `IncludeAbsorbedFluxBackground` | 吸収フラックスの再注入（260919Cl） | false |

格子の種類 `GridType` は結晶系から自動で決まる: `MasterPattern.ShouldUseHexGrid` = `UseHexGridEnabled && CrystalSystemNumber ∈ {5, 6}`（三方晶・六方晶）なら六方格子、それ以外は正方格子。`UseHexGridEnabled` は static フラグ（既定 true、コメントでは「デバッグ用」）。

### 1.2 処理の流れ

1. `CreateDirectionsSquare` / `CreateDirectionsHex` で $+Z$ 半球と $-Z$ 半球のセル中心の方向を作り、**$+Z$ → $-Z$ の順に連結**して `directions` にする（六方格子の無効セルは `null`）。
2. `BetheMethod.RunEBSDNew(MaxNumOfBloch, Energies, I, Depths, directions, Solver, Thread, UseNonLocalAbsorption, IncludeTDSBackground, IncludeAbsorbedFluxBackground)`。結晶の回転は単位行列。`BetheMethod` の設定は
   - `EbsdRepresentativeDirectionBlockSize` = 1（GridSize ≤ 512）、3（≤ 1024）、5（それ以上）
   - `UseLocalSurfacePerBeamDirection = false`（コメント: 新経路では結晶を回すので局所表面近似は使わない）
3. 完了後、別の `BackgroundWorker` で `Disks[ei][di].Amplitudes`（複素数、長さ 2 × 半球セル数）を平面へ焼き込む:

$$ \texttt{PositivePlanes}[k][i] = |A_i|^2,\qquad \texttt{NegativePlanes}[k][i] = |A_{i+N_{\rm hemi}}|^2,\qquad k = e\cdot n_D + d . $$

   不足分は 0。平面ごとに独立なので `Parallel.For`。焼き込み済みの `Disks[ei][di]` は逐次 `null` にしてピークメモリを下げる（コメント: grid 512 で複素振幅の総計 ~2.6 GB）。
4. 進捗は Bethe 段 99 %、焼き込み段 1 % で配分（`BetheStageWeight`、`CompileStageWeight`）。

> 平面の値は `CBED_Disk.Amplitudes` の絶対値の 2 乗である。この量が「振幅」なのか、すでに強度の平方根に相当する量なのかは `BetheMethod` 側の約束事なので、回折担当の文書で確認すること（要確認事項 M-1）。

`BetheMethod` は、`MasterPattern.GetMasterPatternSquareSymmetryOperations` / `GetMasterPatternHexSymmetryOperations`（§3.3）で格子上の等価セルを対応づける（`Diffraction/BetheMethod.cs` L950–980 付近）。

---

## 2. `MasterPattern` の格納形式

### 2.1 配列

- `PositivePlanes[k]`、`NegativePlanes[k]`: `float[GridSize²]`。$k = e\cdot n_D + d$（**energy-major**）。欠落は長さ 0 か `null`。
- `GetPlane(hemisphere, e, d)`: 範囲外は `null`。
- `Energies`、`Depths`: 格子の値（コピーではなく参照で保持）。
- `PlaneCount` $= n_E n_D$。

### 2.2 半球

`Hemisphere.PositiveZ` / `NegativeZ` は**結晶直交系の $z$ の符号**で分けた半球（[ebsd-geometry.md](ebsd-geometry.md) §1.3 の $z\parallel\mathbf c$）。どちらの半球も、平面上の $(x,y)$ は $+z$ 側から見た向きで同じ式（§3.1）に写す（$-Z$ 半球を鏡映しない）。

### 2.3 深さの意味

`MasterPattern.DepthIntervals` の doc による定義:

- `Depths[j]` $=t_j$ の平面は、厚さ $t_j$ までの**累積応答** $M(t_j) = \int_0^{t_j} R(z)\,dz$。
- 区間 $j$ の平均応答は

$$ \bar R_j = \frac{M(t_j)-M(t_{j-1})}{\Delta t_j},\qquad \Delta t_0 = t_0,\ \ \Delta t_j = t_j - t_{j-1}\ (t_{-1}=0). $$

- `ComputeDepthIntervals`: 幅が 0 以下・非有限の区間は 1 を返す（割らない）。MC 分布の重み作りも同じ関数を使う。
- `EbsdMonteCarloDistribution.Depths` の doc: 深さ格子は**出射方向に沿った経路長**（「マスターパターンは接球面近似で、方向ごとに等価表面法線 = −出射方向なので」）。

---

## 3. Rosca–Lambert 等積投影

### 3.1 正方格子

半辺 $L=\sqrt{\pi/2}$（`SquareLimit`）、刻み $\Delta = 2L/N$。セル $(w,h)$（$w$ = 列、$h$ = 行、線形添字 $h N + w$）の中心:

$$ a = -L + (w+\tfrac12)\Delta,\qquad b = L - (h+\tfrac12)\Delta \quad(\text{行が下へ進むと } b \text{ が減る}). $$

**正方形 → 球面**（`RoscaLambertToSphereSquare`、コメントでは Rosca の square → disk 式）:

$$
|b|\le|a|:\ \ (A,B) = \tfrac{2a}{\sqrt\pi}\bigl(\cos\tfrac{\pi b}{4a},\ \sin\tfrac{\pi b}{4a}\bigr),\qquad
|a|<|b|:\ \ (A,B) = \tfrac{2b}{\sqrt\pi}\bigl(\sin\tfrac{\pi a}{4b},\ \cos\tfrac{\pi a}{4b}\bigr),
$$

$$ \rho^2 = A^2+B^2,\qquad (x,y,z) = \Bigl(\sqrt{1-\rho^2/4}\,A,\ \sqrt{1-\rho^2/4}\,B,\ \pm\bigl(1-\rho^2/2\bigr)\Bigr). $$

$(a,b)=(0,0)$ は $(0,0,\pm1)$。角 $(L,L)$ は $(A,B)=(1,1)$、$z=0$（赤道）。正方形の面積 $2\pi$ = 半球の立体角なので等積。

**球面 → 正方形**（`SphereToRoscaLambertSquare`、コメントでは Shirley/Rosca 逆変換）: 正規化した $(x,y,z)$ について

$$ r_s = \sqrt{(1+|z|)/2}\ \ (\text{符号で分岐。実装は } z\ge0 \text{ なら } (1+z)/2,\ \text{そうでなければ } (1-z)/2),\qquad (A,B) = (x,y)/r_s, $$

$$
|A|\ge|B|:\ a = \operatorname{sgn}(A)\,\sqrt{A^2+B^2}\,\tfrac{\sqrt\pi}{2},\ \ b = \tfrac{4a}{\pi}\arctan\tfrac{B}{A};\qquad
|A|<|B|:\ b = \operatorname{sgn}(B)\sqrt{A^2+B^2}\,\tfrac{\sqrt\pi}{2},\ \ a=\tfrac{4b}{\pi}\arctan\tfrac{A}{B}.
$$

コメント（260325Cl）: `atan2` ではなく `atan` を使う（`atan2` は $a<0$ で $\theta\pm\pi$ を返して誤る）。投影側は常に $|z|$ を渡すので $z\ge0$ の枝だけを使う。

**補間**（`InterpolatePlaneSquare`）: $w_f = (a+L)/\Delta - 0.5$、$h_f = (L-b)/\Delta - 0.5$ のバイリニア。添字を $[0,N-1]$ にクランプする（赤道を越えて反対半球の値は使わない）。`EbsdPatternComposer.BuildLookupTable` の正方格子経路は、$w_0<0$ なら $(w_0,f_w)=(0,0)$、$w_0\ge N-1$ なら $(N-2,1)$ とクランプする点だけが異なる（値は同じになる）。

### 3.2 六方格子（三方晶・六方晶）

- $N = (\texttt{GridSize}-1)/2$。軸座標 $(u,v)$、有効セルは $|u|\le N,\ |v|\le N,\ |u+v|\le N$。線形添字 $(v+N)(2N+1)+(u+N)$。無効セルの方向は `null`、平面値は 0。
- 刻み $\Delta_h = P_g/N$、$P_g = 2\sqrt\pi/3^{3/4}$（`Hex_preg` = 六角形の外接円半径。260331Cl に apothem からこちらへ修正）。
- 軸座標 → 直交: $(x,y) = \bigl((u+v/2)\Delta_h,\ \tfrac{\sqrt3}{2}v\Delta_h\bigr)$。
- 六角形 → 円板（`HexagonToDisk`、EMsoft `Lambert2DHexForward` に準拠とコメント）: 座標を $(X,Y)=(h_y,h_x)$ と入れ替え、sextant（`GetSextant`、EMsoft 準拠）を判定して sextant 0 へ回し、

$$ q = \frac{Y_0\,\pi/(2\sqrt3)}{X_0},\quad (XX_0, YY_0) = 3^{1/4}\sqrt{2/\pi}\;X_0\,(\cos q,\ \sin q), $$

  を元の sextant へ戻し、入れ替えを戻す。逆（`DiskToHexagon`）は $x_0 = r_0\,3^{-1/4}\sqrt{\pi/2}$、$y_0 = 3^{-1/4}\sqrt{6/\pi}\,r_0\,\operatorname{atan2}(YY_0,XX_0)$。
- 円板 ↔ 球面は正方格子と同じ式（$\rho^2$、$\sqrt{1-\rho^2/4}$、$\pm(1-\rho^2/2)$）。球面 → 円板は $(A,B)=\sqrt{2/(1+|z|)}\,(x,y)$。
- 補間（`InterpolatePlaneHex`、`GetHexBarycentricLookup`）: 分数軸座標を cube rounding で最近接セル $(r_u,r_v)$ に丸め、オフセット $(f_u,f_v)$ の符号から 6 つの三角形のどれかを決めて 3 近傍の重心座標で補間。`GetHexBarycentricLookup` は無効セルの重みを残り 2 点へ半分ずつ配り直す（`InterpolatePlaneHex` は無効セルを 0 として扱うだけ）。
- 表示・変換用: `RenderHexPlaneToImage`（六方座標のまま正方画像へ）、`ResampleHexToSquarePlane`（正方 Rosca–Lambert 格子へ再標本化。$+Z$ 半球の式で方向を作る）。

### 3.3 格子上で厳密な対称操作

補間なしに添字を写せる操作だけを列挙している。

- 正方格子（4/mmm 系、`SymmOper`、16 種）: 恒等、$z$ まわり 90°/180°/270°、左右・上下反転、対角鏡映 2 種、半球反転 $m_z$、反転、$\bar4$ 型の合成、$x$/$y$/対角軸まわりの 2 回回転（半球反転を伴う）。`GetMasterPatternSquareSymmetryOperations` は点群番号ごとの部分集合を返す。単斜晶（点群 2, m, 2/m）は `MainAxis` で主軸を判定、2mm 系は `PointGroupHMStr`、$\bar42m$/$\bar4m2$ は `StrSE2p` で分岐。三方・六方・立方晶では、正方格子で厳密に写せる操作だけ（立方晶の 3 回軸などは含まれない）。
- 六方格子（6/mmm 系、`HexSymmOper`、24 種）: $z$ まわり 60° 刻みの回転、鏡映 6 種（`MirrorA0`〜`A5`）、これらと半球反転の合成、$C_2'$ 6 種。点群 16–27 で分岐（32/3m/$\bar3$m は `StrSE2p`、$\bar6$m2/$\bar6$2m も `StrSE2p`）。

`TransformMasterPatternSquareIndex` / `TransformMasterPatternHexIndex` が 2 半球を連結した添字（$\text{sphere}\cdot N_{\rm hemi} + \text{local}$）を写す。

---

## 4. 検出器への投影

### 4.1 `EbsdPatternProjector`（照合用、単一スライス）

- コンストラクタ / `Rebuild(geometry)`: ラスター $W_r\times H_r$ の各画素について、検出器の画素座標 $\mathrm{col} = (c+0.5)/W_r\cdot W - 0.5$、$\mathrm{row} = (r+0.5)/H_r\cdot H - 0.5$ の**視線**（試料系）`PixelToSampleDirection` をキャッシュする。方位に依らない。`Rebuild` はバッファを再利用する（幾何較正で 1 評価ごとに幾何が変わるため。260920Cl）。
- `Project(mp, R, posPlane, negPlane, output)`: 各画素で $\mathbf d^{\rm c} = R^{-1}\mathbf d^{\rm s}$、$d_z^{\rm c}\ge0$ なら `posPlane`、そうでなければ `negPlane`。正方格子は `SphereToRoscaLambertSquare(d_x, d_y, |d_z|)` + `InterpolatePlaneSquare`、六方格子は `SphereToRoscaLambertHex` + `GetHexBarycentricLookup`。片方の平面が無ければその半球は 0。
- `parallel=false` で行ループを逐次にできる（方位単位で並列化する呼び出し向け）。
- doc の注記: 単一スライス専用（速度優先）。MC 重み付きの合成は `EbsdMonteCarloDistribution.ComposeGlobalWeightedPattern` で 1 枚の平面に縮約してから渡す（[ebsd-pattern-composition.md](ebsd-pattern-composition.md) §5）。

### 4.2 面内回転分解（辞書総当たり用、正方格子のみ）

方位 $R = r_0 R_z(\phi)$ について、結晶系の視線は $\mathbf d^{\rm c} = R_z(-\phi)\,(r_0^{-1}\mathbf d^{\rm s})$。$\mathbf u = r_0^{-1}\mathbf d^{\rm s}$ の Lambert 円板上の極座標は $\phi$ に依らないので、球点ごとに 1 回だけ（`PrepareSpherePoint`）

$$ q_0 = \tfrac{2}{\pi}\operatorname{atan2}(u_y,u_x),\quad r_a = \frac{\sqrt{u_x^2+u_y^2}/|\mathbf u|}{\sqrt{(1+|u_z|/|\mathbf u|)/2}}\cdot\frac{\sqrt\pi}{2},\quad r_b = \tfrac{4}{\pi} r_a,\quad \text{neg} = [u_z<0] $$

を計算し、各 $\phi$ では（`ProjectInPlane`）$g = q_0 - 2\phi/\pi$、$s=\lfloor g+\tfrac12\rfloor$、局所角 $t = (g-s)\pi/2\in[-\pi/4,\pi/4)$、

| $(s+8)\bmod 4$ | $a$ | $b$ |
|---|---|---|
| 0 | $r_a$ | $r_b t$ |
| 1 | $-r_b t$ | $r_a$ |
| 2 | $-r_a$ | $-r_b t$ |
| 3 | $r_b t$ | $-r_a$ |

としてバイリニア補間するだけにする。doc: `SphereToRoscaLambertSquare` の逆変換と数学的に等価（$\arctan(\tan t) = t$）。数値は ULP 程度の差（「等価群 — dict 回帰で検証」）。Lambert 後の $(a,b)$ 平面での 2 次元回転は不可（Rosca 写像が回転と可換でないため）とコメントにある。

### 4.3 `EbsdPatternComposer.BuildLookupTable`（表示用）

表示ラスター全体について、視線 → Rosca–Lambert → 補間係数を 1 回だけ作る（エネルギー・深さに依らない）。正方格子は「左上セルの添字 + $(f_w, f_h)$」、六方格子は「3 添字 + 3 重み」を画素ごとに保持し、半球フラグを別に持つ。視線の式は [ebsd-geometry.md](ebsd-geometry.md) §3.6。合成モデルは [ebsd-pattern-composition.md](ebsd-pattern-composition.md)。

---

## 5. ファイル形式（`MasterPatternFile`、拡張子 `.rpmp`、FormatVersion = 1）

目的（doc）: grid 512 の構築は数分かかり、MC の乱数で (energy × depth) 格子が毎回わずかに変わるので、同じ条件で探索や表示を繰り返すときに構築をやり直さない。

### 5.1 レイアウト（リトルエンディアン）

| 位置 | 型 | 内容 |
|---|---|---|
| 0 | 4 byte | マジック `"RPMP"` |
| 4 | int32 | 版（1） |
| 8 | int32 | ヘッダ長 $n_h$（0〜$2^{24}$） |
| 12 | $n_h$ byte | ヘッダ = UTF-8 JSON（`MasterPatternFileInfo`） |
| 以降 | zlib | 本体 |

本体（`BinaryWriter`）: `GridSize` (int32)、`GridType` (int32)、$n_E$ (int32)、$n_D$ (int32)、`Energies` (double × $n_E$)、`Depths` (double × $n_D$)、続いて **$+Z$ 半球の全平面 → $-Z$ 半球の全平面**を $k=0..n_En_D-1$ の順に「長さ (int32、0 は欠落) + float32 列」、最後に電子の本数 (int32) と電子（1 本 61 byte: float32 × 15 + フラグ 1 byte）。

電子 1 本の並び: Depth、Vec.X/Y/Z、Position.X/Y、Energy、TotalEnergyLoss、LastInelasticDepth、LastInelasticEnergyBeforeLoss、LastInelasticEnergyAfterLoss、LastInelasticDirection.X/Y/Z、LastDecoherenceDepth、フラグ（bit0 = HasLastInelasticEvent、bit1 = HasLastDecoherenceEvent）。

### 5.2 ヘッダ（`MasterPatternFileInfo`）

`FormatVersion`、`CreatedUtc`（ISO 8601）、`Creator`、`CrystalName`、`CrystalFingerprint`、`BeamEnergyKeV`、`MaxNumOfBloch`、`UseNonLocalAbsorption`、`IncludeTDSBackground`、`IncludeAbsorbedFluxBackground`、`Note`、`BackscatteredElectronCount`、`McSampleTiltRad`（既定 NaN）、`McEnergyThresholdKeV`（既定 NaN）、`McIncidentCount`。NaN・±∞ は `JsonNumberHandling.AllowNamedFloatingPointLiterals` で文字列として書く（261028Cl。それ以前は NaN で保存が例外になっていた）。

### 5.3 結晶の照合

`CrystalFingerprint`: `cell=a,b,c(nm, F6),α,β,γ(deg, F4);sg=<SymmetrySeriesNumber>;` に続けて原子ごとに `Z:x,y,z(F5),Occ(F4),B_iso,eff(nm², F6)` と、異方性原子なら $\beta_{11},\beta_{22},\beta_{33},\beta_{12},\beta_{23},\beta_{31}$（E6）を足して `;`。名前は含めない。`CheckCompatibility` は指紋の完全一致と $|E_{\rm file}-E|\le10^{-6}$ keV を要求し、不一致なら理由の文字列を返す。

### 5.4 読み書きの工夫

- 保存は `path.tmp` に書いてから `File.Move(overwrite)`（途中失敗で既存ファイルを壊さない）。圧縮は既定 `Fastest`（float はあまり縮まないため）。
- `Load(path, loadElectrons)`: `false` なら電子の区画を読まない（MnAl 1e7 入射で約 545 万本 ≈ 0.7 GiB の連続配列を避ける）。
- `ReadElectronChunks(path, chunkSize=512)`: 電子を 512 本（69,632 B < LOH 閾値 85,000 B）ずつの配列で返す。平面は 64 KiB の捨てバッファへ読み飛ばす。
- 電子を保存する理由（doc）: MC 分布（源の深さモード・非晶質層・蛍光体重み・格子）を変えると電子からの再ビニングが要るが、電子を持っておけば MC をやり直さずに済み、毎回同じ分布になる。MC 分布そのものは保存しない。

---

## 6. 実装上の疑問点

### M-I1 正方格子の対角鏡映の名前が入れ替わっている（機能への影響なし）

`TransformMasterPatternSquareIndex` で $x = w - c$、$y = c - h$（$c=(N-1)/2$、$b$ は行と逆向き）と置くと

- `MirrorDiagXY` $(w,h)\mapsto(N-1-h,\,N-1-w)$ は $(x,y)\mapsto(y,x)$ = **$y=x$ に関する鏡映**（コメントは「$y=-x$」）。
- `MirrorDiagX_Y` $(w,h)\mapsto(h,w)$ は $(x,y)\mapsto(-y,-x)$ = **$y=-x$ に関する鏡映**（コメントは「$y=x$」）。

`GetMasterPatternSquareSymmetryOperations` はどの点群でも両者を必ず対で入れているので、結果は変わらない。`RotDiagXY180`（$[110]$ 軸 2 回）と `RotDiagX_Y180`（$[1\bar10]$ 軸 2 回）の写像はコメント通り。

---

## 7. 要確認事項

- **M-1** 平面値 $|A|^2$ の物理的意味（`CBED_Disk.Amplitudes` が何の振幅か）と、累積 $M(t)$ という解釈（`DepthIntervals` の doc）が Bethe 側の出力定義と一致しているか。回折担当の文書と突き合わせること。
- **M-2** 正方格子の補間は添字をクランプするだけで、赤道付近では反対半球の値を使わない。grid が粗いときに赤道付近で半セル程度の不連続が生じうるが、許容としているか。
- **M-3** 三方晶・六方晶で六方格子を使うとき、結晶直交系の $x$ 軸（$\parallel\mathbf a^*$、[ebsd-geometry.md](ebsd-geometry.md) §1.3）と六角形の sextant の向き（`HexagonToDisk` の座標入れ替え）の対応が、`HexSymmOper` の鏡映軸 A0〜A5 の定義（A0 = 0°、A1 = 30°…）と整合しているか（`MirrorA*` と `RotA*_180` の選び方は点群ごとにコメントで説明されているが、軸の基準角の根拠は書かれていない）。
- **M-4** 立方晶（点群 28–32）は正方格子で厳密に写せる操作だけを使う。3 回軸による等価は利用していない（計算量の削減が限定的）ことを意図通りとしてよいか。
- **M-5** `MasterPatternFile.FormatVersion` = 1 のまま、261105Cl に指紋の書式（異方性原子の $\beta$ の追加）が変わった。等方の結晶は以前と同じ文字列なので互換だが、異方性の結晶の旧ファイルは読み込まれなくなる（doc に明記あり）。版を上げない判断でよいか。
- **M-6** `MasterPattern` のコンストラクタは `Energies`/`Depths`/平面配列を複製せずに保持する。`MasterPatternBuildRequest` 側は複製しているが、ファイル読み込み経路などで外部から配列を書き換えたときの扱い。
