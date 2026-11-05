本稿はコードから Claude が起こした下書き（2026-11-05）。著者未査読。

# 実装上の疑問点（回折・動力学）

各手法文書を書く過程で見つけた、実装に誤りの疑いがある点・コメントと実装の食い違いを集めた。**コードは変更していない。**
詳しい根拠（式・反例）は各文書の「実装上の疑問点」節にある。行番号は 2026-11-05 時点の `main`（9ffcf47）。

重要度の目安: **A** = 計算結果に影響しうる / **B** = 条件次第で影響（潜在的） / **C** = コメント・doc の不整合のみ。

## A: 計算結果に影響しうるもの

| # | 場所 | 内容 | 詳細 |
|---|---|---|---|
| A1 | `Diffraction/BetheMethod.cs:3865–3866, 2877, 2895, 2282–2283` | U キャッシュ `uDictionary` のキーは `(compose(g), compose(h) or int.MaxValue)` だけで、kV と検出器角（inner/outer）を含まない。STEM は 2877 行で Clear したあと、ADF 検出器域 [inner, outer] に限った U′ を、通常の U と同じキーで `q+g_i−g_j` の全組（qList は $g_i-g_j$ の全差なので 000 も含む）に書き込み、終了時に消さない（StemDoWork 内の Clear は 2523・2535・2877 の 3 か所のみ）。実部は inner/outer に依らないので、汚染されるのは虚部（吸収）だけ。`GetDifractedBeamAmpriltudes` は `AccVoltage != voltage` のときしか Clear しない（2282 行。261014Cl で CBED・PED・STEM には無条件の Clear が入ったが、ここは条件付きのまま）。そのため同一インスタンス・同一電圧で STEM の後に呼び、N か方位が変わって再計算に入ると、`getPotentialMatrix` の要素のうち STEM が触れた指数は検出器域の U′、それ以外は全立体角の U′ という**混在した吸収行列**になる。HRTEM（`GetHRTEMImage` → `GetDifractedBeamAmpriltudes`）も同じ。インスタンスは `Crystal.Bethe`（`Crystal/Crystal.cs:167, 477`）として結晶ごとに 1 個で共有されるので、ReciPro で STEM と回折シミュレータが同じ結晶を使えば起こりうる（ReciPro 側は未確認）。EBSD は muBack の前後で Clear し（`:1645, :1657`）、マスターパターン用の別インスタンス（`EBSD/EBSD.cs:233`）を使うので影響しない。Kikuchi プロファイルも別インスタンス（`KikuchiProfileCalculator.cs:64`） | [bloch-wave](bloch-wave.md) Q1、[kikuchi-dynamical](kikuchi-dynamical.md) §4-9 |
| A1b | `BetheMethod.cs:2292, 2528, 2337–2338` | **A1 の調査で見つけた別件（詳細は下の「A1b の詳細」）**: STEM が共有フィールド `Beams` を |g|² 順に並べ替えて上書きする。その後 `GetDifractedBeamAmpriltudes` が古い固有ベクトルを再利用すると、回折波の振幅が別の反射に代入される | 下記 |
| A2 | `BetheMethod.cs:2741–2742` と `:3735, :3763` | 時間的コヒーレンス包絡の係数が STEM と HRTEM で異なる。STEM は $\exp[-\pi^2\lambda^2\delta^2(k_q^2-k^2)^2]$、HRTEM（TCC・準コヒーレントとも）は $\exp[-\tfrac12\pi^2\lambda^2\delta^2(\cdot)^2]$。同じ `delta` が STEM では √2 倍に効く。δ の定義（1σ か否か）が両者で違う意図なら、その旨の明記が必要 | [bloch-wave-imaging](bloch-wave-imaging.md) |
| A3 | `BetheMethod.cs:1600–1608, 1630`、`EBSD/EBSD.cs:235, 1191` | 六方格子でマスターパターン格子幅 G>512（ブロック幅 3 または 5）のとき、ブロック中心が無効セル（方向 `null`）に落ちうる。するとブロック全体が [00−1] で選んだビーム集合を共有する | [kikuchi-dynamical](kikuchi-dynamical.md) §4-1 |
| A4 | `BetheMethod.cs:2120–2168` `ComputeAbsorbedFluxBackground` | $N(z)=\sum_g|\psi_g|^2$ は密度の和。260829Cl の行割り修正後、吸収なしで保存されるのは $\sum_g P_g|\psi_g|^2$ なので、$U'=0$ でも $D(t)\neq0$ になりうる（傾いた表面・HOLZ で顕在化） | [kikuchi-dynamical](kikuchi-dynamical.md) §4-2 |
| A5 | `Diffraction/Alchemi/AlchemiReduction.cs:146–171, 194` | 同種の不整合: $L_\mathrm{coh}$ は $P_g/P_0$ 重みの電流、$Y_\mathrm{dyn}$ は密度で数える。$\mu=\mu_{00}I$ の極限で Total $\neq\mu_{00}t/V_c$。差は HOLZ ビームと傾いた表面のときのみ | [alchemi-ionization](alchemi-ionization.md) |
| A6 | `AlchemiReduction.Yield :196` | $Y_\mathrm{dyn}$ に $2K/P_0\ (\approx1/\cos\theta_0)$ の因子がない疑い（導出は Claude による・要検算）。30 mrad で 4.5e-4 程度 | [alchemi-ionization](alchemi-ionization.md) |
| A7 | `BetheMethod.cs:2292, 2370` | 回折振幅・PED が固有値分解を再利用するかの判定が N・電圧・方位だけで、原子・DW 因子・Surface・モデル設定の変更を見ない | [bloch-wave](bloch-wave.md) Q5 |
| A8 | `Atom/Atoms.cs:355` → `Crystal/Crystal.cs:1644` | 運動学的な電子の F は割当イオン（SubNumberElectron）の Peng 因子を使い、`BetheMethod.ElasticIonModel`（既定 Neutral）を見ない。イオンを割り当てた結晶で運動学と動力学の $f_e$ が別模型になる（設計意図なら明記を） | [structure-factor](structure-factor.md) S1 |

### A1b の詳細（2026-11-05 追加調査）

**再利用の判定が通ることの確認**
- 固有値・固有ベクトル（`EigenValues/EigenVectors/EigenVectorsInverse`）を書くのは `GetDifractedBeamAmpriltudes`（2310–2318 行）だけ。STEM・CBED・PED・ALCHEMI は書かず、null にもしない。
- `RunSTEM` は `MaxNumOfBloch/AccVoltage/BaseRotation` を呼び出し側と同じ値に設定し（2496–2501 行）、worker が `Beams = [.. Find_gVectors(BaseRotation, vecK0).OrderBy(e => e.Vec.Length2)]` とする（2528 行）。
- `Matrix3D.Equals` は 9 成分の値の比較（`Mathematics/Matrix.cs:650`）なので、同じ方位なら `rotation.Equals(BaseRotation)` は真。
- `Find_gVectors` は同じ入力なら同じ集合を返す。候補は逐次ループで集め（4418–4425 行）、`QuickSelect` のピボットは末尾要素で決定的（`NthItemUtils.cs`）。STEM の vecK0 も U0 を Clear 後に計算し、同じ電圧なら同じ値になる。よって `Beams.Length == EigenValues.Length` も真。
- 以上から、`GetDifractedBeamAmpriltudes(N,V,R,t₁)` → `RunSTEM(N,V,R,…)` → `GetDifractedBeamAmpriltudes(N,V,R,t₂)` の 3 回目は 2292 行の条件をすべて満たし、再計算せずに古い固有分解を使う。

**何が起こるか**
- 固有分解は `Find_gVectors` の順（Rating = |g|·Q² の昇順、4472 行）で作られている。一方 `Beams` は STEM の |g|² 順。
- `alpha = V⁻¹ e₀` の e₀ は両方の順で 000（Rating も |g|² も 0）なので、入射条件は正しい。
- 2337–2338 行の `Beams[i].Psi = psi_atZ[i]` で、旧順の i 番目の反射の振幅が新順の i 番目の反射に代入される。出射面の位相 `exp(iπ P_i t)`（2331 行）は新順の P を使う。P は実数なので |ψ|² は反射間の**置換**になり、位相も合わない。
- 透過波（000）だけは正しい。

**現実に起きる手順**
- HRTEM(N,V,R) → STEM(N,V,R) → HRTEM（デフォーカスや厚さだけ変更）。`GetHRTEMImage` は厚さごとに `GetDifractedBeamAmpriltudes(BlochNum, AccVol, rot, t)` を呼ぶ（3723 行）ので、最初の呼び出しから再利用に入り、像全体が誤る。
- 回折振幅(N,V,R) → STEM(N,V,R) → 厚さだけ変えて回折振幅。
- インスタンスが `Crystal.Bethe` として結晶ごとに共有されること（A1 と同じ）が前提。ReciPro 側の呼び方は未確認。

**並びがどれだけ違うか（Claude による簡易計算・近似）**
Rating の式だけを Python で再現し（Si、200 kV、[001] 近傍。U0・接平面での重複除去・BFS キャッシュの上限は省略）、Rating 順の先頭 N+1 本と、それを |g|² 順に並べ直したものを比べた。

| 傾き | N | 位置が変わる本数 | うち |g| の違う反射と入れ替わる本数 |
|---|---|---|---|
| 0° | 50 | 0/51 | 0 |
| 0° | 200 | 200/201 | 196 |
| 0.5° | 50 | 45/51 | 33 |
| 2° | 200 | 200/201 | 198 |

- 晶帯軸ちょうどで ZOLZ だけのときは Q = −g² なので Rating = |g|⁵ となり、|g|² 順と一致する。入れ替わるのは |g| の等しい組の中だけで、対称等価なら強度も等しいので目に見えない（ただし Si の 333/511 のように |g| が偶然等しい非等価な組は除く）。
- HOLZ リングの反射（Q≈0 で |g| が大きい）が選ばれるほど N が大きいとき、または少しでも傾いているときは、ほぼ全部が別の反射と入れ替わる。

**関連する点**
- CBED（457 行）は同じ `Find_gVectors` 呼び出しなので並びは変わらない。PED（2460 行）は和集合を Rating 順に並べ直すので、本数が N+1 と一致しても大半は同じ並びになる。
- STEM の worker は `Beams` フィールドを 2536 行以降も 23 回読み直す（2647, 2759, 2765 行など）。STEM の実行中に同じインスタンスで `GetDifractedBeamAmpriltudes` が再計算に入ると、STEM の途中で `Beams` が別の配列に差し替わる。ポテンシャル行列（2536 行で作成）と `qList`（2759 行）が別のビーム集合になり、N が違えば範囲外アクセスにもなりうる。`IsSTEM_Busy`（317 行）はあるが、`GetDifractedBeamAmpriltudes` 側では見ていない。
- 対処の方向（コードは変更していない）: STEM は `Beams` フィールドではなく局所変数でビームを持つ。あるいは、固有分解を作ったときのビーム配列の参照を保持し、2292 行の条件に `ReferenceEquals(Beams, eigenBeams)` を加える。

## B: 条件次第・潜在的なもの

| # | 場所 | 内容 |
|---|---|---|
| B1 | `BetheMethod.cs:3900–3918` | 異方性 DW のとき、g 方向の実効値 $m=h'^\mathsf{T}\beta h'/s_g^2$ を等方 B として TDS 積分に入れる。被積分の $K\pm g/2$ は別方向なので近似（[scattering-factors](scattering-factors.md) Q4） |
| B2 | `Diffraction/Kikuchi/KikuchiProfileCalculator.cs:94–113` | `ComputeCrystalHash` が `Dsf.UseIso` を含まず、UseIso だけ切り替えるとスナップショットを新鮮と誤判定する |
| B3 | `KikuchiProfileCalculator.cs:261–267` | 2 波経路の $I_\mathrm{ref}$ を +g 側からしか取らず、+g 側がガードで 0 を返すと c が跳ねる（現状の不変条件では未到達） |
| B4 | `KikuchiProfileCalculator.cs:362` | 「x=+1 が +g の Bragg 線」とあるが、member +g の Q=0 は x≈−1。出射ビーム側の呼び方なら正しい。ReciPro の描画の前提を要確認 |
| B5 | σₙ / 菊池カーネル / `getU` | B 未設定原子の扱いが不統一（σₙ は下限 1e-3 nm²、菊池カーネルは 5e-3 nm² を代用、`getU` は下限なし）。「吸収しないが後方散乱する」原子になる |
| B6 | `BetheMethod.cs:2100, 2258`、`KikuchiProfileCalculator.cs:733` | $F_{jj'}$ は $\lambda t$ が非常に小さいとき桁落ちする（判定が λ 自身に対して行われる）。通常の吸収の大きさでは実害なし |
| B7 | `Diffraction/IonizationChannel.cs` `IonizationTableShape` | `_truncated/_usedTail` が可変のまま共有され、IonizationData を run 間で使い回すと前回のフラグが残りうる |
| B8 | `IonizationChannelTable.GridAt` | E0 が表の端の行より外だと PCHIP が黙って外挿する可能性（範囲チェックは 30–400 keV のみ） |
| B9 | `BetheMethod.cs:3450–3452, 3481` | 環境変数 `RECIPRO_STEM_EDX_HERMTOL=NaN` で Hermitian 検査が素通りする |
| B10 | `Atom/Atoms.cs:366–371` | 中性子の同位体加重が Dictionary の列挙順に依存（仕様上の保証なし） |
| B11 | `Atom/AtomStatic.cs:3174` `FactorImaginaryAnnularFlatEwald` | 動径積分の Jacobian R が欠け、前因子 $\gamma/(\pi^2k_0)$ も他版と合わない。**本番では未使用** |
| B13 | `Diffraction/Alchemi/BetheMethod.Alchemi.cs:73, 258` | ALCHEMI は U0 を `getU(AccVoltage)` で読む（73 行）が、Clear は 258 行（ポテンシャル行列の直前）にしかない。前の run が別電圧で残した値があると、U0 の実部（γ を含む）が古い電圧のものになり、$K=\sqrt{k^2+U_0}$ と Q・P がずれる。`BuildFixedUnion` 内の `Find_gVectors` が Beam に入れる U もキャッシュ由来だが、ポテンシャル行列は Clear 後に作り直すので、影響は U0 経由だけ。261014Cl で CBED・PED・STEM の入口に入れた「U0 を読む前に Clear」が ALCHEMI にだけない |
| B12 | 位相の符号 | `Crystal.GetStructureFactor` は $e^{-2\pi i h\cdot r}$、`BetheMethod.getU` は $e^{+2\pi i g\cdot r}$（`BetheMethod.cs:3921`。教科書の $U_{-g}$ に相当）。それぞれ内部で一貫しているかは各文書で確認したが、運動学と動力学の結果を突き合わせる箇所があれば要注意（[structure-factor](structure-factor.md) S2、[bloch-wave](bloch-wave.md)） |

## C: コメント・doc と実装の不一致（数値は正しい）

| # | 場所 | 内容 |
|---|---|---|
| C1 | `Atom/AtomStatic.cs:1209–1213` | `MottBetheMonopoleCoefficient` の doc は「s² [Å⁻²]・f [nm]」だが、0.2393366 は s² [nm⁻²]・f [nm] の係数（Å 系では 0.023934 Å）。Claude が値を検算済み |
| C2 | `AtomStatic.cs` `FactorImaginary*` の returns、`BetheMethod.cs:3918` | 「無次元」とあるが実際は $f_e$ と同じ nm 単位 |
| C3 | `BetheMethod.cs:3734–3735` | HRTEM 準コヒーレント分岐で Ec/Es（時間的/空間的）のラベルが逆。β の項が空間的、δ の項が時間的。Claude がコードで確認済み |
| C4 | `Diffraction/Alchemi/AlchemiAngularSpread.cs:48–50` | 「4σ（≈1.7 FWHM）より外は重み <1e-7」とあるが、実際は $e^{-8}\approx3.3\times10^{-4}$。Claude が検算済み |
| C5 | `AtomStatic.cs:2883` と `:1222` | Peng の当てはめ範囲が「s≤6」（260811Cl）と「s≤2 用」（260923Cl）で食い違う |
| C6 | `BetheMethod.cs:815–825` | `GetCrystalFixedMasterEquivalentBeamAndSurface` は丸め誤差を除けば常に (d, −d) を返す恒等変換。結果は正しいが、コメントの説明と処理が合わない |
| C7 | `KikuchiProfileCalculator.cs:7, 592`、`AlchemiReduction.cs:12` | 260829Cl 修正後も旧式 `/P_col`、旧表記 $S=B^\dagger\mathrm{diag}(\sigma)B$ が残る |
| C8 | `IonizationTableShape`、`AlchemiBasisDiagnostic`、`BuildFixedUnion:186–189`、`IonizationFsTable.SCount` | 「s>4」「s≤8 Å⁻¹」「81 点」など旧データ版の数値が残る（現行は 321 点・16 Å⁻¹） |
| C9 | `AlchemiMuBuilder.cs:206` | コメントは「BisoEffective と同じ規則」、実装は `Biso000`（G=0 では一致し実害なし） |

## その他の観察（誤りではない）

- `Atom/NistElastic.cs` は EBSD の MonteCarlo 専用で、回折の経路からは参照されていない。
- Temari 因子は `Fe` が TDS 吸収の被積分関数で使われるだけで、`Fx` は未使用。X 線の f は WK のまま。
- Cromer–Mann 表と 8-Gaussian 表はライブラリ内で使われていない。
- `Diffraction/DynamicTheory.cs`（マルチスライス）は `[Obsolete]` で、中身はすべてコメントアウトされている。
- native の `_EBSDSolver` / `_EBSDSolverWithTDS` / `STEM_InelasticQ` はソースがこのリポジトリにないため、式を確認できていない。
