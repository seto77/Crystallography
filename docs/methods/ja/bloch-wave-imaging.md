本稿はコードから Claude が起こした下書き（2026-11-05）。著者未査読。

# Bloch 波法による像計算 — HRTEM・STEM（ADF/TDS）・ポテンシャル像

対象: `Diffraction/BetheMethod.cs` の `GetHRTEMImage`（L3709–3834）、`RunSTEM` / `StemDoWork`（L2478–3634）、`GetPotentialImage`（L3640–3690）。
固有値問題・ビーム選択・`getU` の規約は [bloch-wave.md](bloch-wave.md) を前提とする（nm 単位、$k = 1/\lambda$、平面波 $e^{+2\pi i\mathbf k\cdot\mathbf r}$、角度はラジアン）。
STEM-EDX（イオン化チャネル）の行列要素と規格化は [alchemi-ionization.md](alchemi-ionization.md) と `Diffraction/Ionization/README.md` を参照。ここでは EDX が TDS と共有する部分だけ触れる。

---

## 1. HRTEM `GetHRTEMImage`

### 1.1 入力

| 引数 | 単位・意味（コードの doc コメント） |
|---|---|
| `BlochNum` | Bloch 波（ビーム）数 N |
| `AccVol` | kV |
| `rot` | 結晶方位 |
| `aperture (R, X, Y)` | 対物絞りの半角とシフト。doc は「mrad」だが `ExtractInsideBeams` は `2 sin(R/2)/λ` で使うので**この関数に届く値はラジアン**（→要確認 1） |
| `size`, `res` | 画素数、nm/pixel |
| `cs` | 球面収差 $C_s$ [nm] |
| `beta` | 照射半角 $\beta$ [rad]（空間的コヒーレンス） |
| `delta` | 焦点ずれ $\delta = C_c\,\Delta V/V$ [nm]（時間的コヒーレンス） |
| `thicknesses`, `defocusses` | [nm] の配列 |
| `quasiMode` | true: 準コヒーレント（既定） / false: TCC |
| `native` | true かつ native 有効なら C++ の格子版ソルバ |

260610Cl: 波長はフィールド `AccVoltage` ではなく引数 `AccVol` から計算する（以前は新規インスタンスで λ = ∞ → 全画素 NaN）。

### 1.2 手順

1. 厚みごとに `GetDifractedBeamAmpriltudes(BlochNum, AccVol, rot, t)`（[bloch-wave.md](bloch-wave.md) §8。2 回目以降は固有値分解を使い回す）で出射振幅 $\psi_\mathbf g$。
2. `ExtractInsideBeams`: $|\mathbf g_{xy} - \mathbf c| < 2\sin(R/2)/\lambda$、$\mathbf c = (2\sin(X/2)/\lambda,\ 2\sin(Y/2)/\lambda)$。$R = \infty$ なら全ビーム。
3. 入射の傾き $\mathbf k_{xy}$ は 0 に固定（`var k = new PointD(0, 0)`）。$\mathbf q = \mathbf g_{xy}$。

### 1.3 収差関数と包絡関数

$$
\chi(q) = \pi\lambda q^2\Big(\Delta f + \tfrac12 C_s\lambda^2 q^2\Big),\qquad \text{伝達} = e^{-i\chi(q)}
$$

（$\Delta f$ = `defocus`。この符号では $C_s > 0$ に対する Scherzer 条件は $\Delta f < 0$ 側。→要確認 2）

**準コヒーレント**（各ビームに独立に包絡を掛ける）:

$$
L(q) = e^{-i\chi(q)}\;
\exp\!\big[-\pi^2\beta^2 q^2(\Delta f + \lambda^2 C_s q^2)^2\big]\;
\exp\!\big[-\tfrac12\pi^2\lambda^2\delta^2 q^4\big]
$$

$$
I(\mathbf r) = \Big|\sum_\mathbf g \psi_\mathbf g\,L(q_\mathbf g)\,e^{2\pi i\mathbf g\cdot\mathbf r}\Big|^2
$$

第 2 因子が空間的、第 3 因子が時間的コヒーレンス（コード L3734–3735 のコメントは逆に書いてある → [bloch-wave.md](bloch-wave.md) 疑問点 Q3）。

**TCC**（相互透過係数）: $g \le h$ の対について

$$
T(\mathbf q_1,\mathbf q_2) = e^{-i\chi(q_1)}e^{+i\chi(q_2)}
\exp\!\big[-\pi^2\beta^2\,|\Delta f(\mathbf q_1-\mathbf q_2) + \lambda^2C_s(q_1^2\mathbf q_1 - q_2^2\mathbf q_2)|^2\big]
\exp\!\big[-\tfrac12\pi^2\lambda^2\delta^2(q_1^2-q_2^2)^2\big]
$$

係数は $g = h$ で $|\psi_g|^2$、$g<h$ で $2\psi_g\psi_h^*$。$I(\mathbf r) = \big|\mathrm{Re}\sum_{g\le h} (\ldots)\,T\,e^{2\pi i(\mathbf g-\mathbf h)\cdot\mathbf r}\big|$。$T(\mathbf q_2,\mathbf q_1) = T(\mathbf q_1,\mathbf q_2)^*$ なので半分の対で実部を取るのは全対の和と等価。

### 1.4 画素座標

$\mathbf r = \big((x - c_X)\,\text{res} + s_X,\ -(y - c_Y)\,\text{res} + s_Y\big)$、$c = $ 画素数/2、$\mathbf s = R\,(\mathbf a+\mathbf b+\mathbf c)/2$（`Crystal.RotationMatrix` で回した単位胞の体対角の半分）。native 版は同じ格子を始点+刻みで渡し、行方向を位相回転の再帰で進め 256 画素ごとに `Exp` で再アンカーする（260610Cl、差 ~1e-13）。像は規格化しない（準コヒーレントで $|\cdot|^2$、TCC で $|\mathrm{Re}|$ の生値）。

---

## 2. STEM `RunSTEM` → `StemDoWork`

### 2.1 入力

| 引数 | 単位・意味 |
|---|---|
| `maxNumOfBloch` | N |
| `voltage` | kV |
| `cs` | $C_s$ [nm] |
| `delta` | 焦点ずれ $\delta$ [nm] |
| `sliceThickness` | 非弾性（TDS/EDX）積分の深さ刻みの目安 [nm] |
| `imageSize`, `resolution` | 画素数、nm/pixel |
| `sourceSize` | 光源像の大きさ。像に掛けるガウスぼかしの**半値半幅** [nm]（`GaussianBlurFast(…, sourceSize/resolution)`、引数名 `hwhm`） |
| `baseRotation`, `thicknesses`, `defocusses` | 方位、[nm]、[nm] |
| `beamDirections` | 入射方向の配列。**正方格子**（一辺 $\sqrt{\text{Length}}$）として扱い、内接円の内側だけ計算 |
| `convergenceAngle` | 収束半角 $\alpha$ [rad] |
| `detAngleInner`, `detAngleOuter` | ADF 検出器の内角・外角 [rad] |
| `solver`, `thread` | [bloch-wave.md](bloch-wave.md) §6 |
| `ionizations` | STEM-EDX のチャネル（空なら EDX なし） |

`RunSTEM` は busy なら何もしない（260711Cl）。配列は複製してから worker に渡す（260801Cl）。結果は成功時に一度だけ `ResultStem`（`StemSimulationResult`）として公開し、`ResultSTEM` はその互換ビュー。

### 2.2 Bloch 波の部分

- 中心方向で `Find_gVectors` → **$|\mathbf g|^2$ 昇順に並べ直す**（`OrderBy(e => e.Vec.Length2)`）。入射ベクトル $\psi_0 = (-1, 0, \ldots)$。ポテンシャル行列（全立体角 TDS 吸収の局所形式）は 1 回だけ組む。
- 各入射方向 $\mathbf k$（円内のみ）で `reset_gVectors` → `getEigenMatrix` → EVD。透過係数

  $$\mathrm{tc}_\mathbf g(\mathbf k, t) = e^{2\pi i (P_\mathbf g/2)\,t}\sum_j C^{(j)}_\mathbf g e^{2\pi i\gamma^{(j)}t}\alpha_j$$

  をユーザー厚み $t$ と全スライス深さ $z$ について持つ。
- 260711Cl の「全スライス物化」経路（native Eigen かつメモリ見積が有利なとき、環境変数 `RECIPRO_STEM_ALLSLICE` で強制）は EVD 直後に全深さの tc を 1 回の GEMM で作り、固有ベクトルを保持しない。そうでない経路は固有値・α を畳み込んだ固有ベクトルを保持し、スライスごとに `GenerateTC1B` で再生成する。どちらも同じ式。

### 2.3 レンズ・絞り・検出器関数（L2730–2756）

$$
L(\mathbf k, \mathbf k+\mathbf q; \Delta f) = \exp\!\Big\{-i\big[\chi(k) - \chi(|\mathbf k+\mathbf q|)\big]\Big\}\,
\exp\!\Big\{-\big[\pi\lambda\delta\,(|\mathbf k+\mathbf q|^2 - k^2)\big]^2\Big\}
$$

$\chi$ は §1.3 と同じ定義（コメント: 「大塚さんの資料『瀬戸先生・質問3への回答(瀬戸追記).docx』を参照」）。時間的コヒーレンスの指数に HRTEM の $\tfrac12$ が無い（→ [bloch-wave.md](bloch-wave.md) 疑問点 Q2）。空間的コヒーレンスはレンズ関数に入れず、最後に実空間ガウスぼかし（`sourceSize`）で入れる。

- 収束絞り: $A(\mathbf k) = [\,|\mathbf k_{xy}|^2 \le (\sin\alpha/\lambda)^2\,]$
- 検出器: $D(\mathbf k_f) = [\,(\sin\theta_\mathrm{in}/\lambda)^2 \le |\mathbf k_{f,xy}|^2 \le (\sin\theta_\mathrm{out}/\lambda)^2\,]$

$\mathbf k_{xy}$ は `getVecK0` の結果の XY 成分（表面平行成分は真空と同じ）。

### 2.4 弾性（コヒーレント ADF）

- $\mathbf q$ の集合 `qList`: ビーム差 $\mathbf g_i - \mathbf g_j$ のうち、ある円内方向 $\mathbf k$ で $A(\mathbf k)$ かつ $A(\mathbf k+\mathbf q)$ を満たすもの。$|\mathbf q|^2$ 昇順。
- $\mathbf k+\mathbf q$ の振幅は方向格子の 4 近傍から**双一次補間**（`list[kIndex]` に近傍添字 `n0..n3` と重み `r`）。格子座標への写像は $d = q_x\cdot(r_\mathrm{pix}-0.5)/\max_k k_x + (r_\mathrm{pix}-0.5)$（y は符号反転）で、`beamDirections` が $\mathbf k_{xy}$ について等間隔の正方格子であることを前提にしている。4 近傍すべてが円内でないと捨てる。

$$
I_\mathrm{el}(\mathbf q; t, \Delta f) = \sum_{\mathbf k} L(\mathbf k,\mathbf k+\mathbf q;\Delta f)
\sum_{\substack{\mathbf g:\ D(\mathbf k+\mathbf g)\\ \mathbf g-\mathbf q\in\text{Beams}}}
\mathrm{tc}_\mathbf g(\mathbf k,t)\,\Big[\sum_{m=0}^{3} r_m\,\mathrm{tc}_{\mathbf g-\mathbf q}(\mathbf k_m,t)\Big]^*
$$

（出射方向 $\mathbf k+\mathbf g = (\mathbf k+\mathbf q)+(\mathbf g-\mathbf q)$ で同じ検出器画素に入る 2 成分の干渉。）

### 2.5 TDS（インコヒーレント ADF）

吸収ポテンシャルで失われた流束のうち、検出器域に散乱される分を深さ方向に積分する:

$$
I_\mathrm{TDS}(\mathbf q; t, \Delta f) = \sum_{z_s \le t} \frac{2\pi}{k_\mathrm{vac}}\,\Delta z\;
\sum_\mathbf k L(\mathbf k,\mathbf k+\mathbf q;\Delta f)\;
\boldsymbol\psi(\mathbf k+\mathbf q, z_s)^\dagger\,\mathsf U'_\mathbf q\,\boldsymbol\psi(\mathbf k, z_s)
$$

$$
(\mathsf U'_\mathbf q)_{ij} = \overline{\texttt{getU}(\mathbf q + \mathbf g_i - \mathbf g_j;\ \theta_\mathrm{in},\theta_\mathrm{out}).\texttt{Imag}}
$$

- `getU` の環状版（局所形式。[bloch-wave.md](bloch-wave.md) §2.4）で、$f'$ を検出器角 [inner, outer] に限って積分したもの。非局所形式はコメントアウト。
- 複素共役（L2895「共役とると、なぜかいい感じ。」）: `getU` の位相は $e^{+2\pi i\mathbf G\cdot\mathbf r}$ なので `getU(G)` は教科書の $U'_{-\mathbf G}$。$f'$ が実数なので共役を取ると $U'_{\mathbf G}$ になる。$\psi^*_{\mathbf k+\mathbf q,i}\,\psi_{\mathbf k,j}$ を結ぶフーリエ成分は $\mathbf G = \mathbf q + \mathbf g_i - \mathbf g_j$ の $U'_\mathbf G$（教科書の符号）なので、共役は符号規約の補正として正しい（Claude の解析。経験則ではなく必然と考えられる）。
- 係数 $2\pi/k_\mathrm{vac}$ は吸収率 $d|\psi|^2/dz = -(2\pi/K)U'|\psi|^2$ に対応（[bloch-wave.md](bloch-wave.md) §2.6 と同じ）。
- **深さの刻み**: 区間 $[T_{t-1}, T_t]$（$T_{-1} = 0$）を $n = \max(1, \lfloor \Delta T/\text{sliceThickness}\rfloor)$ 等分し、評価点は各小区間の**右端** $z_s = T_{t-1} + e\,\Delta z$（$e = 1..n$）。右端の Riemann 和なので誤差は $O(\Delta z)$。厚み方向は累積（$I(T_t) = I(T_{t-1}) + \ldots$）。
- 計算経路（`STEM_InelasticQ` native、q ごとの U ストリーミング `RECIPRO_STEM_USTREAM`、managed fallback）は和の順序だけが違い、数式は同じ（260610Cl, 260711Cl のコメント）。

### 2.6 実空間像

$$
I_X(\mathbf R) = \frac{1}{r_\mathrm{pix}^2}\,\mathrm{Re}\sum_\mathbf q I_X(\mathbf q)\,e^{2\pi i\mathbf q\cdot\mathbf R},\qquad
\mathbf R = \big(\text{res}(x-c_X),\ -\text{res}(y-c_Y)\big) + \mathbf s_{xy}
$$

- $r_\mathrm{pix}$ = 方向格子の半径（画素）。この $r_\mathrm{pix}^2$ は「相対強度表示の既存規約」で、絶対規格化ではない（260802Cl。EDX 側は q = 0 のエントリ数 = 絞り内方向数で割るように直されたが、ADF 側は回帰互換のため変えていない）。
- **実部を取る**（260802Cl。以前は `.Magnitude`）。$I(-\mathbf q) = I(\mathbf q)^*$ が補間の非対称で $O(h^2)$ だけ破れる虚部は数値誤差なので直交射影で落とす、という理由がコメントにある。残差（`MeasureHermitian`）と捨てた虚部の最大値・負の最小画素値は診断として記録するだけで、ADF では例外にしない（EDX は閾値超えで例外、さらに ±q の対称化を行う）。
- その後 `sourceSize > 0` ならガウスぼかし（HWHM = sourceSize/resolution 画素）、負値を 0 にクランプ、`Both = Ela + TDS`。

### 2.7 STEM-EDX との共有

EDX は TDS と同じ透過係数・同じ q エントリ・同じレンズ関数・同じ累積式を使い、行列 $\mathsf U'_\mathbf q$ を `FillIonizationUq`（$U_\mathrm{ion}$、負位相を直接生成、`uDictionary` とは別キャッシュ）に替えるだけ（`AccumulateInelasticQ`、260802Cl）。EDX は native ヘルパ必須（managed 経路では `NotSupportedException`）。

---

## 3. ポテンシャル像 `GetPotentialImage`（L3640–3690）

与えたビーム列の `Ureal`（$\hat U_\mathbf g$）と `Uimag`（$U'_\mathbf g$）を実空間に合成する:

$$
S_\mathrm{cry}(\mathbf r) = \sum_\mathbf g \hat U_\mathbf g\,e^{2\pi i\mathbf g_{xy}\cdot\mathbf r},\qquad
S_\mathrm{th}(\mathbf r) = \sum_\mathbf g U'_\mathbf g\,e^{2\pi i\mathbf g_{xy}\cdot\mathbf r}
$$

出力 4 枚: `phase = true` なら $|S_\mathrm{cry}|, \arg S_\mathrm{cry}, |S_\mathrm{th}|, \arg S_\mathrm{th}$、false なら実部・虚部。座標は

$$
\mathbf r = \big(-(x-c_X)\,\text{res} + s_X,\ -(H-1-y-c_Y)\,\text{res} + s_Y\big)
$$

で、HRTEM/STEM の $\mathbf r$ に対して X・Y とも符号が反転している。`getU` の係数が教科書の $U_{-\mathbf g}$ であること（[bloch-wave.md](bloch-wave.md) §2.1）と合わせると $S(\mathbf r) = U(-\mathbf r)$ なので、反転した座標で正しい向きの $U$ が出る（Claude の解析）。ただし Y は $(H-1-y)$ を使うので HRTEM/STEM と 1 画素ずれる（→要確認 4）。単位は nm⁻²（$U$ そのもの）。

---

## 要確認事項

1. `GetHRTEMImage` の doc は絞り半角を「mrad」と書くが、`ExtractInsideBeams` は `2 sin(R/2)/λ` で使う（ラジアン）。呼び出し側（ReciPro）が mrad→rad 換算して渡しているか。
2. デフォーカスの符号: $\chi = \pi\lambda q^2(\Delta f + \tfrac12 C_s\lambda^2q^2)$ なので負の $\Delta f$ が不足焦点。UI 表示の符号と一致しているか。
3. STEM の時間的コヒーレンス包絡に $\tfrac12$ が無い件（[bloch-wave.md](bloch-wave.md) Q2）。`delta` の定義（rms か 1/e 半幅か）を決めたい。
4. `GetPotentialImage` の Y 座標が $(H-1-y)$ で、HRTEM/STEM（$y$）と 1 画素ずれる。重ね合わせ表示で問題になるか。
5. STEM の TDS 深さ積分が右端の Riemann 和（L2578）。中点則にすると $O(\Delta z^2)$ になるが、既定 `sliceThickness` との兼ね合いで十分か。
6. 方向格子の双一次補間は $\mathbf k_{xy}$ の等間隔正方格子を前提にしている（`maxK = max k_x`, `coeff2`）。`beamDirections` を作る側（ReciPro）がこの前提を満たしているか（球面上の等角度格子だと $\sin$ の分だけ非等間隔になる）。
7. ADF 像の規格化 $1/r_\mathrm{pix}^2$ は EDX（q = 0 のエントリ数）と比が $\pi/1.05^2$ 程度違う（260802Cl コメント）。ADF と EDX を同じ図で比べる用途があるか。
8. HRTEM は入射の傾き $\mathbf k_{xy} = 0$ 固定（L3716）。傾斜照明は扱わない仕様か。
9. STEM の方向格子の外縁は収束角より大きい（L2791 のコメント「収束角ではないことに注意(5%大きい)」）。格子座標への写像は `maxK`（計算した方向の最大 $k_x$）で、絞り判定 $A$ は $\sin\alpha/\lambda$ で行う。両者の関係が呼び出し側の格子生成と一致しているか。

（要確認 9 件。実装上の疑問点は [bloch-wave.md](bloch-wave.md) の Q1–Q3 に含めた）
