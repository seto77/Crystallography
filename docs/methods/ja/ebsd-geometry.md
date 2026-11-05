本稿はコードから Claude が起こした下書き（2026-11-05）。著者未査読。

# EBSD の座標系・検出器幾何・方位の表し方

対象コミット: `9ffcf47`（main）。対象ファイル: `EBSD/EbsdDetectorGeometry.cs`、`EBSD/EbsdPatternComposer.cs`（`UpdateTiltCoefficients`・`BuildLookupTable`・`ExitRay`）、`EBSD/MonteCarlo.cs`（試料面の定義）、`EBSD/EbsdIndexer.cs`（回転のユーティリティ）、`EBSD/EbsdDictionaryIndexer.cs`（`GetProperRotations`）、`Crystal/Crystal.cs`（結晶直交座標系）。

EBSD 関連の全クラスはこの文書の約束事を共有している。`EbsdDetectorGeometry` の doc には「規約は表示フォワードモデル (`EbsdPatternComposer`) と厳密に一致させる」とあり、両者の式が一致することは `EbsdIndexCheck` ハーネスで確認したと書かれている（ハーネスはこのリポジトリには無い）。

---

## 1. 座標系

### 1.1 実験室系（lab）

モンテカルロ（`MonteCarlo`）と検出器幾何は同じ lab 系を使う。

- 原点 = 電子線の照射点（試料表面上）。
- 入射電子の初期方向は $\mathbf{v}_0 = (0,0,-1)$（$-z$ 向き）。
- 試料は $x$ 軸まわりに傾斜角 $\tau$（`MonteCarlo.Tilt`、`EbsdDetectorGeometry.SampleTilt`、rad）だけ傾いている。点 $\mathbf{p}$ が試料内にある条件は

$$ p_y \tan\tau \ge p_z . $$

- したがって試料表面の外向き法線（lab）は $\hat{\mathbf n}_s^{\rm lab} = (0,\,-\sin\tau,\,\cos\tau)$、表面からの垂直深さは

$$ d = p_y\sin\tau - p_z\cos\tau \quad(\text{MC では } d \mathrel{+}= s\,(v_y\sin\tau - v_z\cos\tau) \text{ で更新}). $$

  真空側は $-y$・$+z$ 側である（$\tau>0$ のとき）。

### 1.2 試料系（sample）

lab → 試料系は $x$ 軸まわりの回転 $R_x(\tau)$（`EbsdDetectorGeometry.LabToSample`。コメントによれば OpenTK `CreateRotationX` の $M\cdot v$ 作用と同一）:

$$
s_x = l_x,\qquad s_y = \cos\tau\, l_y + \sin\tau\, l_z,\qquad s_z = -\sin\tau\, l_y + \cos\tau\, l_z .
$$

逆変換は `SampleToLab`。この回転で $\hat{\mathbf n}_s^{\rm lab}\mapsto(0,0,1)$ になるので、**試料系の $+z$ は試料表面の外向き法線**である。入射ビームは試料系で $(0,-\sin\tau,-\cos\tau)$。

同じ回転は `EbsdMonteCarloDistribution` の ctor と `ComputePathLengthUpperBound`（電子の射出方向 `Vec`〔lab〕を試料系へ戻す）でも使われる。MC 分布に渡す傾斜は「MC を回したときの値」であって UI の現在値ではない、と doc に明記されている。

### 1.3 結晶直交座標系（crystal Cartesian）

`Crystal.cs`（L684 付近）で格子ベクトルを次のように置く:

$$
\mathbf c = (0,0,c),\quad
\mathbf b = (0,\, b\sin\alpha,\, b\cos\alpha),\quad
\mathbf a = \Bigl(a\sqrt{1-\cos^2\beta-\tfrac{(\cos\gamma-\cos\alpha\cos\beta)^2}{\sin^2\alpha}},\ a\,\tfrac{\cos\gamma-\cos\alpha\cos\beta}{\sin\alpha},\ a\cos\beta\Bigr).
$$

つまり **$\mathbf c\parallel z$、$\mathbf b$ は $yz$ 面内、$x\parallel \mathbf b\times\mathbf c\ (\parallel\mathbf a^*)$**。逆格子ベクトル・晶帯軸・Kikuchi 反射の方向（`VectorOfG_KikuchiLine`）はすべてこの直交系の成分で扱われる。

> 他ソフトとの違い: EBSD でよく使われる流儀（例: $x\parallel\mathbf a$、$z\parallel\mathbf c^*$）とは軸の取り方が異なる（Claude 追記・要確認）。六方晶では $\mathbf b\parallel y$、$\mathbf a$ は $xy$ 面内で $x$ から $-30^\circ$ の向き（$x\parallel\mathbf a^*$）になる。他ソフトの方位データと突き合わせるときは、この結晶系の取り方の差を回転に含める必要がある。

---

## 2. 方位の表し方

### 2.1 回転行列

EBSD のコードは方位を**回転行列だけ**で扱う。`EbsdOrientationCandidate.Rotation`・`EbsdZoneAxisSolution.Rotation`・`EbsdMatchingContext.Rotation` はいずれも「結晶 → 試料系の回転（`Crystal.RotationMatrix` と同じ規約）」:

$$ \mathbf v^{\rm s} = R\,\mathbf v^{\rm c}. $$

- $R$ の列は、結晶直交系の基底ベクトルを試料系で表したもの。
- 受動的に読めば「結晶系の成分 → 試料系の成分」への座標変換、能動的に読めば「試料系の軸を結晶の軸へ回す回転」。
- 投影では逆を使う: 試料系の方向 $\mathbf d^{\rm s}$ を結晶系 $\mathbf d^{\rm c}=R^{-1}\mathbf d^{\rm s}$ へ戻してからマスターパターンを引く（`EbsdPatternProjector.Project`、`EbsdPatternComposer.BuildLookupTable` とも `rotation.Inverse()`）。
- `Matrix3D` の 9 引数コンストラクタは列優先（`EbsdIndexer.SolveWahba` のコメント）。

### 2.2 オイラー角

`EBSD/` 以下のコードはオイラー角を使わない（入力も出力も回転行列）。リポジトリには `Mathematics/Euler.cs` に Z-X-Z の `FromMatrix` があるが、EBSD のコードからは呼ばれていない。Bunge 角 $(\varphi_1,\Phi,\varphi_2)$ との対応は EBSD のコードでは定義されていない（要確認事項 G-1）。

参考（Claude 追記・要確認）: Bunge の方位行列 $g$ を「試料系 → 結晶系」($\mathbf v^{\rm c} = g\,\mathbf v^{\rm s}$) とする流儀なら $R = g^{\mathsf T}$。ただし結晶直交系の取り方（§1.3）と試料系の取り方（§1.2）が他ソフトと違うので、$R$ と他ソフトの $g$ をそのまま比べることはできない。

### 2.3 対称等価・方位差

- 点群の proper 回転（`EbsdDictionaryIndexer.GetProperRotations`）は、対称操作の線形部 $W$（格子座標、$\det W=+1$ のもの）を $S = M\,W\,M^{-1}$（$M$ = `MatrixReal`、列が $\mathbf a,\mathbf b,\mathbf c$）で直交系へ移し、重複と単位元を除いたもの。P1 などで対称が無ければ `null`。
- 等価な方位は $R\,S$（右から掛ける）。
- 方位差角（`EbsdIndexer.MisorientationDeg`）:

$$ \theta(R_1,R_2) = \min_{S\in\{I\}\cup\mathcal S}\ \arccos\!\Bigl(\tfrac{\operatorname{tr}(R_1^{\mathsf T}R_2S)-1}{2}\Bigr). $$

  `properSymmetries` を渡さなければ $S=I$ のみ。

### 2.4 方位の摂動

方位の精密化（Nelder–Mead）で使う摂動は**試料系の左摂動**（`EbsdIndexer.PerturbRotation`）:

$$ R' = \mathrm{Rot}\bigl(\hat{\boldsymbol\omega},|\boldsymbol\omega|\bigr)\,R_0,\qquad \boldsymbol\omega=(\omega_x,\omega_y,\omega_z)\ [\deg]. $$

### 2.5 方位格子

- `EbsdIndexer.FibonacciSphereRotation(di, nSphere)`: Fibonacci 球面格子の点 $\hat{\mathbf n}_i$（$z_i = 1-2(i+\tfrac12)/N$、方位角 $i\cdot\pi(3-\sqrt5)$）に $+z$ を回す最小回転。
- 辞書探索・Radon 粗探索の方位は $R = R_{\rm sphere}(i)\,R_z(\phi_j)$。$R\,\hat{\mathbf z}^{\rm c} = \hat{\mathbf n}_i$ なので、球点は**結晶の $z$ 軸（$\mathbf c$）の試料系での向き**、$\phi_j$ はそのまわりの回転。SO(3) 上で厳密に一様な格子ではない（要確認事項 G-2）。

---

## 3. 検出器幾何（`EbsdDetectorGeometry`）

### 3.1 パラメータ

| 記号 | プロパティ | 単位 | 意味 |
|---|---|---|---|
| $\delta$ | `DetTilt` | rad | 検出器傾斜。doc の既定 90°（$\pi/2$） |
| $(X_d,Y_d,Z_d)$ | `DetX`, `DetY`, `DetZ` | mm | 検出器中心パラメータ（UI の数値欄の値） |
| $s$ | `PixelSize` | mm/px | ピクセルサイズ |
| $W,H$ | `WidthPx`, `HeightPx` | px | 画素数 |
| $x_m$ | `XMirror` | ±1 | 左右反転（`FormEBSD.DetectorXMirror` と同じ） |
| $\tau$ | `SampleTilt` | rad | 試料傾斜 |

ctor は非有限値、$s\le0$、$W,H\le0$、$x_m\ne\pm1$、検出器面が原点を通る配置（$|\mathbf n\cdot\mathbf C|<10^{-12}$）を例外で拒否する。

### 3.2 派生量（lab 系）

$$
\mathbf C = (X_d,\,-Y_d,\,-Z_d),\quad
\mathbf e_x=(1,0,0),\quad
\mathbf e_y=(0,-\cos\delta,-\sin\delta),\quad
\mathbf n=(0,\sin\delta,-\cos\delta).
$$

$\{\mathbf e_x,\mathbf e_y,\mathbf n\}$ は正規直交。カメラ長（照射点から検出器面への垂直距離）は

$$ L = |\mathbf n\cdot\mathbf C| \quad(\texttt{CameraLength} = \texttt{FormEBSD.CameraLength2}). $$

### 3.3 ピクセル → 方向

画素 $(\mathrm{col},\mathrm{row})$ は**コーナー原点・ピクセル中心 $(\mathrm{col}+0.5,\mathrm{row}+0.5)$** の規約。画像中心基準の表示 mm 座標は

$$ u = (\mathrm{col}+0.5-W/2)\,s,\qquad v=(\mathrm{row}+0.5-H/2)\,s \quad(\texttt{PixelToMm}). $$

検出器面上の点（lab）:

$$ \mathbf P = \mathbf C + x_m\,u\,\mathbf e_x + v\,\mathbf e_y \quad(\texttt{PixelToLabPoint}). $$

$\delta=90^\circ$ では $\mathbf e_y=(0,0,-1)$ なので、row が増えると lab の $-z$（下）へ進む。

- **視線**（`PixelToSampleDirection`）: $\hat{\mathbf d}^{\rm s} = R_x(\tau)\,(-\mathbf P/|\mathbf P|)$。検出器上の点から照射点へ向かう向き。
- **射出方向**（電子が試料から出ていく向き）: $\hat{\mathbf k}_{\rm out}^{\rm s} = -\hat{\mathbf d}^{\rm s} = R_x(\tau)\,\hat{\mathbf P}$。晶帯軸の指数付け（`EbsdZoneAxisIndexer`）、MC 分布の参照（`EbsdPatternComposer.ExitRay`、`ComputeDetectorAverageSliceWeights`）はこちらを使い、符号を明示的に反転している。
- **マスターパターンを引く方向は視線 $\hat{\mathbf d}$（$=-\hat{\mathbf k}_{\rm out}$）**: `EbsdPatternProjector` と `EbsdPatternComposer.BuildLookupTable` は $\hat{\mathbf d}^{\rm c}=R^{-1}\hat{\mathbf d}^{\rm s}$ の結晶系 $z$ 成分の符号で半球を選び、Rosca–Lambert 座標へ写す（[ebsd-master-pattern.md](ebsd-master-pattern.md)）。マスターパターンがどちらの向き（入射方向か射出方向か）で計算されているかは動力学計算側の約束事で、回折担当の文書を参照（要確認事項 G-3）。

逆写像 `SampleDirectionToPixel(k_out)`: $\mathbf v = R_x(\tau)^{-1}\mathbf k_{\rm out}$、$t = (\mathbf n\cdot\mathbf C)/(\mathbf n\cdot\mathbf v)$、$t\le0$ または $|\mathbf n\cdot\mathbf v|<10^{-12}$ なら `null`。$\mathbf q=t\mathbf v-\mathbf C$ から $u = x_m\,\mathbf q\cdot\mathbf e_x$、$v=\mathbf q\cdot\mathbf e_y$、$\mathrm{col}=u/s+W/2-0.5$、$\mathrm{row}=v/s+H/2-0.5$。

### 3.4 パターン中心（PC）

PC は**照射点から検出器面へ下ろした垂線の足** $\mathbf F=(\mathbf n\cdot\mathbf C)\,\mathbf n$。検出器中心基準の物理面内座標（`PatternCenterMm`）は

$$ u_F = -X_d,\qquad v_F = -(Y_d\cos\delta + Z_d\sin\delta). $$

画像ピクセル座標（`GetPatternCenter`）は表示側で $x$ に $x_m$ が掛かる:

$$ \mathrm{col}_{\rm PC} = x_m u_F/s + W/2 - 0.5,\qquad \mathrm{row}_{\rm PC} = v_F/s + H/2 - 0.5,\qquad DD = L . $$

つまりこのコードの PC 規約は「**(PC の画素座標〔コーナー原点・ピクセル中心〕, DD〔mm〕)**」で、検出器幅や DD で無次元化した値（Bruker・EDAX/TSL・Oxford・EMsoft などの PC 規約）は実装されていない（要確認事項 G-4）。

逆算 `FromPatternCenter(u_F, v_F, DD, δ)` は次の連立方程式を解く（コード中の式のまま）:

$$ Y_d\cos\delta + Z_d\sin\delta = -v_F,\qquad Y_d\sin\delta - Z_d\cos\delta = -DD,\qquad X_d=-u_F, $$

$$ \Rightarrow\ Y_d = -(v_F\cos\delta + DD\sin\delta),\quad Z_d = -v_F\sin\delta + DD\cos\delta. $$

第 2 式は $\mathbf n\cdot\mathbf C = +DD$ を意味する。doc の「標準配置 signed L = −DD を仮定」との関係は §5 の疑問点 G-I1 を参照。

### 3.5 バンド中心線 ↔ 面法線

画像の直線 $A\,\mathrm{col}+B\,\mathrm{row}+C=0$（ピクセル中心規約）を、表示 mm 座標の直線 $a u + b v + c = 0$ に直す:

$$ a = A/s,\quad b=B/s,\quad c = A c_x + B c_y + C,\qquad c_x = W/2-0.5,\ c_y=H/2-0.5 . $$

バンド中心線と原点を含む平面の法線（lab、正規化、符号不定。`LineToLabNormal`）:

$$ \mathbf g \propto a\,\mathbf e_u + b\,\mathbf e_y + \gamma\,\mathbf n,\qquad \mathbf e_u = x_m\mathbf e_x,\qquad \gamma = \frac{c - a\,(\mathbf e_u\cdot\mathbf C) - b\,(\mathbf e_y\cdot\mathbf C)}{\mathbf n\cdot\mathbf C}. $$

（$\mathbf g\cdot\mathbf P = 0$ が $au+bv+c=0$ と同値になるように選んだもの。）逆は `LabNormalToLine`: $a_{\rm mm}=x_m g_x$、$b_{\rm mm}=\mathbf g\cdot\mathbf e_y$、$c_{\rm mm}=\mathbf g\cdot\mathbf C$ から $(a_{\rm mm}s,\ b_{\rm mm}s,\ c_{\rm mm}-a_{\rm mm}s\,c_x-b_{\rm mm}s\,c_y)$。

バンド縁の点 $\mathbf P_e$ からの Bragg 角（`SinBraggFromEdgePoint`）: $\sin\theta_B = |\hat{\mathbf g}\cdot\mathbf P_e|/|\mathbf P_e|$。

### 3.6 表示ラスターの係数形（`EbsdPatternComposer`）

表示合成は `EbsdDetectorGeometry` を使わず、同じ式を係数で持つ（`UpdateTiltCoefficients`）:

$$
y_c = \cos\tau\cos\delta + \sin\tau\sin\delta,\quad z_c = -\sin\tau\cos\delta + \cos\tau\sin\delta,\quad
y_0 = \cos\tau\,Y_d + \sin\tau\,Z_d,\quad z_0 = -\sin\tau\,Y_d + \cos\tau\,Z_d .
$$

ラスター画素 $(w,h)$ の表示座標は `EbsdRasterView` で

$$ p_x = (2w+1-\mathrm{width})\,S_W + O_X,\qquad p_y = (2h+1-\mathrm{height})\,S_H + O_Y $$

（$O_X,O_Y$ は表示のパン量。$S_W$ は画素ピッチの半分に当たる）。視線（試料系、未正規化）は

$$ \mathbf d^{\rm s} = \bigl(-x_m p_x - X_d,\ \ y_c\,p_y + y_0,\ \ z_c\,p_y + z_0\bigr), $$

これを $R^{-1}$ で結晶系へ回す。§3.3 の $-R_x(\tau)\mathbf P$（$u=p_x$、$v=p_y$）と一致する。射出方向は `ExitRay`: $(x_m p_x + X_d,\ -(y_c p_y + y_0),\ -(z_c p_y+z_0))$。`ExitRay.PlaneDistance` $=|z_0 y_c - y_0 z_c|$ は $L$ に等しい。

---

## 4. 単位の約束

| 量 | 単位 | 備考 |
|---|---|---|
| 検出器座標・PC・DD | mm | |
| 角度（公開 API） | rad（幾何）、deg（許容角・刻み・摂動） | 引数名に `Deg` が付くものは度 |
| 波長（指数付け） | nm | `EbsdIndexer.Index` の doc: 既定 0.00859 nm は 20 kV。Å を渡すと幅尤度が全バンドで頭打ちになる |
| 波数（MC） | nm⁻¹、$k=1/\lambda$（$2\pi$ を含まない） | `UniversalConstants.Convert.EnergyToElectronWaveNumber` |
| 深さ・経路長 | nm | |
| エネルギー | keV（MC・マスターパターン）、eV（TPP-2M・損失） | |

---

## 5. 実装上の疑問点

### G-I1 `FromPatternCenter` が $\mathbf n\cdot\mathbf C=+DD$ の配置を返す

- 根拠: §3.4 の式。`GetPatternCenter` → `FromPatternCenter` の往復は $\mathbf n\cdot\mathbf C>0$ のときだけ恒等になる。
- 反例（`docs` 作成時に式をそのまま Python で評価）: $\delta=90^\circ$、$(X_d,Y_d,Z_d)=(1,20,5)$ では $\mathbf n\cdot\mathbf C=-20$、往復後は $(1,-20,5)$。
- 物理的な向き: §1.1 の MC の試料定義（$\tau=+70^\circ$）では、検出器中心の射出方向の試料系 $z$ 成分は $(1,20,5)$ で $+17.1$（試料表面より上 = 見える側）、$(1,-20,5)$ で $-20.5$（表面より下）。つまり $\tau>0$ で検出器が射出側にある配置は $\mathbf n\cdot\mathbf C<0$ になり、`FromPatternCenter` はそれを鏡像の側へ戻す。
- 影響しうる箇所: `EbsdGeometryCalibrator.Run` の `MakeGeom`（較正の開始点 $(u_{F0},v_{F0},\ln DD_0)$ を `FromPatternCenter` で幾何に戻す）、`EbsdOrientationSearch` の候補較正（L299）。
- 未確認の前提: UI（`FormEBSD`、このリポジトリに無い）での $\tau$ と $Y_d$ の符号の取り方。実運用の配置が $\mathbf n\cdot\mathbf C>0$（たとえば $\tau<0$）なら問題は無く、doc の「signed L = −DD」は $L \equiv Y_d\sin\delta - Z_d\cos\delta\,(=-\mathbf n\cdot\mathbf C)$ を指していることになる。どちらの意味かをコメントに明記するのがよい。

### G-I2 `EbsdDetectorGeometry` の doc 中の「signed L」が未定義

上記の通り、`FromPatternCenter` の doc の「標準配置 signed L = −DD」は $L$ の定義が書かれていないため、$\mathbf n\cdot\mathbf C$ と読むとコードと逆になる。

---

## 6. 要確認事項

- **G-1** Bunge オイラー角との対応（$R$ と $g$ の関係、結晶直交系・試料系の軸の取り方の差を含む）を公式に定義するか。現状 `EBSD/` には変換が無い。
- **G-2** 方位格子 $R_{\rm sphere}(i)R_z(\phi_j)$（Fibonacci × 等間隔 $\phi$）は SO(3) 上で一様ではない（$\hat{\mathbf n}$ の極付近で $\phi$ 方向が密）。粗探索の網羅性の評価にこの非一様性を考慮しているか。
- **G-3** マスターパターンを視線 $\hat{\mathbf d}=-\hat{\mathbf k}_{\rm out}$ で引くことと、Bethe 計算側（`BetheMethod.RunEBSDNew` に渡す `directions`）の方向の意味（入射方向か射出方向か、相反定理の適用）の対応。回折担当の文書で確認すること。
- **G-4** 他社・他ソフトの PC 規約（Bruker PCx/PCy/DD、EDAX/TSL $x^*,y^*,z^*$、Oxford、EMsoft の $x_{pc},y_{pc},L$）との変換式は実装されていない（Claude 追記・要確認）。
- **G-5** 結晶直交系（$\mathbf c\parallel z$、$x\parallel\mathbf a^*$）は他の EBSD ソフトでよく使われる取り方（Claude 追記・要確認）と異なる。データ交換を想定するなら変換を文書化すること。
- **G-6** `XMirror` の物理的な意味（検出器画像の左右反転がカメラ側の反転か表示側の反転か）。コードは「表示 mm 座標 $u$ の $x$ 成分に $x_m$ を掛ける」とだけ定めている。
