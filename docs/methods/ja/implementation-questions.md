本稿はコードから Claude が起こした下書き（2026-11-05）。著者未査読。

# 実装上の疑問点（回折・動力学）

各手法文書を書く過程で見つけた、実装に誤りの疑いがある点・コメントと実装の食い違いを集めた。**コードは変更していない。**
詳しい根拠（式・反例）は各文書の「実装上の疑問点」節にある。行番号は 2026-11-05 時点の `main`（9ffcf47）。

重要度の目安: **A** = 計算結果に影響しうる / **B** = 条件次第で影響（潜在的） / **C** = コメント・doc の不整合のみ。

## A: 計算結果に影響しうるもの

| # | 場所 | 内容 | 詳細 |
|---|---|---|---|
| A1 | `Diffraction/BetheMethod.cs:3865, 2895, 2282` | U キャッシュ `uDictionary` のキーが反射指数だけで、kV と検出器角（inner/outer）を含まない。STEM は ADF 検出器域に限った U′ を通常の U と同じキーで保存し、終了時に消さない。同一インスタンス・同一電圧で後から `GetDifractedBeamAmpriltudes`（HRTEM も経由）を呼ぶと、全立体角の TDS 吸収ではなく検出器域の吸収を読む。EBSD の muBack（`:860`）も同じキーを使う（現状の呼び出し順では起きない） | [bloch-wave](bloch-wave.md) Q1、[kikuchi-dynamical](kikuchi-dynamical.md) §4-9 |
| A2 | `BetheMethod.cs:2741–2742` と `:3735, :3763` | 時間的コヒーレンス包絡の係数が STEM と HRTEM で異なる。STEM は $\exp[-\pi^2\lambda^2\delta^2(k_q^2-k^2)^2]$、HRTEM（TCC・準コヒーレントとも）は $\exp[-\tfrac12\pi^2\lambda^2\delta^2(\cdot)^2]$。同じ `delta` が STEM では √2 倍に効く。δ の定義（1σ か否か）が両者で違う意図なら、その旨の明記が必要 | [bloch-wave-imaging](bloch-wave-imaging.md) |
| A3 | `BetheMethod.cs:1600–1608, 1630`、`EBSD/EBSD.cs:235, 1191` | 六方格子でマスターパターン格子幅 G>512（ブロック幅 3 または 5）のとき、ブロック中心が無効セル（方向 `null`）に落ちうる。するとブロック全体が [00−1] で選んだビーム集合を共有する | [kikuchi-dynamical](kikuchi-dynamical.md) §4-1 |
| A4 | `BetheMethod.cs:2120–2168` `ComputeAbsorbedFluxBackground` | $N(z)=\sum_g|\psi_g|^2$ は密度の和。260829Cl の行割り修正後、吸収なしで保存されるのは $\sum_g P_g|\psi_g|^2$ なので、$U'=0$ でも $D(t)\neq0$ になりうる（傾いた表面・HOLZ で顕在化） | [kikuchi-dynamical](kikuchi-dynamical.md) §4-2 |
| A5 | `Diffraction/Alchemi/AlchemiReduction.cs:146–171, 194` | 同種の不整合: $L_\mathrm{coh}$ は $P_g/P_0$ 重みの電流、$Y_\mathrm{dyn}$ は密度で数える。$\mu=\mu_{00}I$ の極限で Total $\neq\mu_{00}t/V_c$。差は HOLZ ビームと傾いた表面のときのみ | [alchemi-ionization](alchemi-ionization.md) |
| A6 | `AlchemiReduction.Yield :196` | $Y_\mathrm{dyn}$ に $2K/P_0\ (\approx1/\cos\theta_0)$ の因子がない疑い（導出は Claude による・要検算）。30 mrad で 4.5e-4 程度 | [alchemi-ionization](alchemi-ionization.md) |
| A7 | `BetheMethod.cs:2292, 2370` | 回折振幅・PED が固有値分解を再利用するかの判定が N・電圧・方位だけで、原子・DW 因子・Surface・モデル設定の変更を見ない | [bloch-wave](bloch-wave.md) Q5 |
| A8 | `Atom/Atoms.cs:355` → `Crystal/Crystal.cs:1644` | 運動学的な電子の F は割当イオン（SubNumberElectron）の Peng 因子を使い、`BetheMethod.ElasticIonModel`（既定 Neutral）を見ない。イオンを割り当てた結晶で運動学と動力学の $f_e$ が別模型になる（設計意図なら明記を） | [structure-factor](structure-factor.md) S1 |

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
