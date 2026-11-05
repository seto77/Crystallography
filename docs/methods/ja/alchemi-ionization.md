本稿はコードから Claude が起こした下書き（2026-11-05）。著者未査読。

# 内殻イオン化（STEM-EDX）と ALCHEMI の計算手法

対象: `Diffraction/IonizationChannel.cs`、`Diffraction/XrayLineSeries.cs`、`Diffraction/AbsorptionSource.cs`、
`Diffraction/Alchemi/*.cs`、`Diffraction/BetheMethod.cs` の `ComputeIonizationU` / `FillIonizationUq` と STEM-EDX 経路。

- データ（`IonizationFsE0.bin` / `BoteSalvat.bin`）の出所・ライセンス・バイナリ形式・版の履歴は
  [`Diffraction/Ionization/README.md`](../../../Diffraction/Ionization/README.md) が正本。本稿では繰り返さない。
- Bloch 波の固有値問題（`getEigenMatrix` の行列、`P_g`・`Q_g` の定義、`getU`）は [bloch-wave.md](bloch-wave.md)。
- 吸収ポテンシャルの虚部（TDS 吸収形状因子 `FactorImaginary*`）は [scattering-factors.md](scattering-factors.md)。

---

## 1. 記号と単位の約束

| 記号 | 意味 | 単位・規約 | 実装 |
|---|---|---|---|
| $k_\mathrm{vac}$ | 真空中の電子波数 $=1/\lambda$（$2\pi$ を含まない） | nm⁻¹。相対論補正込み $\lambda = 1.22643/\sqrt{V(1+0.97848\times10^{-3}V)}$（$V$ は V） | `UniversalConstants.Convert.EnergyToElectronWaveNumber` |
| $\mathbf G$ | 運動量移送ベクトル（逆格子ベクトルの整数結合） | nm⁻¹。**結合 hkl から `mat * index` で正準再構築**（浮動小数のベクトル和は使わない） | `AlchemiMuBuilder` ctor、`ComputeIonizationU` |
| $s$ | $s=\lvert\mathbf G\rvert/2$（= sinθ/λ 相当。$q=4\pi\sin\theta/\lambda$ ではない） | 内部 API は **nm⁻¹**。F テーブルは **Å⁻¹**。換算は `IonizationTableShape.Evaluate` 内の `sA = s * 0.1` の 1 か所だけ | `INormalizedIonizationShape` |
| $\sigma$ | 内殻（副殻）イオン化の全断面積 | nm²（`BoteSalvat.SigmaNm2` = cm² × 10¹⁴） | `IonizationData.TotalCrossSectionNm2` |
| $F(s)$ | 正規化イオン化形状因子、$F(0)=1$（符号付き） | 無次元 | `IonizationData.Shape` |
| $E_0$ | 入射電子エネルギー | keV。**30–400 keV のみ**（外挿・clamp なし） | `IonizationDataProvider.MinE0KeV/MaxE0KeV` |
| $V_c$ | 単位胞体積 | nm³ | `Crystal.Volume` |
| $M_a(\mathbf G)$ | Debye–Waller 指数 | 等方: $M=B_\mathrm{iso}s^2$（$B$ 型、nm²）。非等方: $M=\beta_{11}H^2+\beta_{22}K^2+\beta_{33}L^2+2\beta_{12}HK+2\beta_{23}KL+2\beta_{31}LH$（$\beta$ は無次元、等価位置ごとに `atom.Operation.ConvertPlaneIndex` で指数を変換） | `AlchemiMuBuilder.Build`、`ComputeIonizationU` |
| $\gamma_j$ | Bloch 波の固有値（表面法線方向の波数補正） | nm⁻¹。$\mathrm{Im}\,\gamma_j>0$ が減衰（$e^{2\pi i\gamma_j z}$） | `AlchemiReduction` ctor 引数注記 |
| $z,t$ | 表面法線方向の深さ・厚み | nm | |

位相は**負位相 $e^{-2\pi i\mathbf G\cdot\mathbf r_a}$ を直接生成**する（`getU` の「正位相を作って共役を取る」規約はイオン化経路に持ち込まない、260801Cl）。

---

## 2. 原子データ: σ と F(s)

### 2.1 全断面積 σ（Bote–Salvat）

`BoteSalvat.SigmaCm2(z, subshell, energyEv)`。副殻番号は 1=K, 2=L1, 3=L2, 4=L3, 5..9=M1..M5。過電圧 $U=E/E_\mathrm{edge}$（`edge` は Bote 表の端、`edgeEvOverride` で差し替え可）。

- $U\le 1$: $\sigma=0$。
- $1<U\le16$（低エネルギー fit）:
  $$\xi=(U-1)\left[\frac{a_1+a_2U+\frac{1}{1+U}\left(a_3+\frac{1}{(1+U)^2}\left(a_4+\frac{a_5}{(1+U)^2}\right)\right)}{U}\right]^2$$
- $U>16$（高エネルギー形）: $\beta^2=E(E+2m c^2)/(E+mc^2)^2$、$x=\sqrt{E(E+2mc^2)}/mc^2$、
  $$\xi=\frac{A_{nlj}}{\beta^2}\frac{U}{U+B_e}\Big[(2\ln x-\beta^2)(1+g_1/x)+g_2+g_3\sqrt{mc^2/(E+mc^2)}+g_4/x\Big]$$
- $\sigma=4\pi a_0^2\,\xi$（$a_0=5.291772108\times10^{-9}$ cm、$mc^2=5.10998918\times10^5$ eV、いずれも xion.f と同値）。

演算順は Python 参照実装 `botesalvat.py` と厳密一致させてある（golden vector 照合、コメント）。文献: Bote & Salvat, PRA 77 (2008) 042701 / xion.f（ADNDT 95 (2009) 871）。σ 自体の不確かさ（K 10 %/L 15 %/M 24 %）は README に記載。

### 2.2 形状因子 F(s, E0) の補間（`IonizationChannelTable.GridAt`）

テーブルは各チャネル（shellCode, Z）について E0 行 × s 節点 321 点（$s=0..16$ Å⁻¹、刻み 0.05）。補間は **E0 方向 → s 方向**の 2 段:

1. **E0 方向**（各 s 節点 $j$ ごと）: 横軸 $x=\ln(u-1)$、$u=E_0/E_\mathrm{th}$（行の $u$ はシリアライズ値＝4 桁丸め）で PCHIP（scipy `PchipInterpolator` 互換、`Pchip` クラス）。
   - その節点で全行が正なら $\ln F$ を補間して $\exp$、非正を含めば符号付き $F$ を直接補間。
   - **$s_\mathrm{cert}\ge s_j$ の行だけを基底に使う**（260810Cl）。低い E0 の行は運動学的に $s_j$ に届かず 0 の埋め草なので混ぜない。行を落とした節点で問い合わせ E0 が残った行の範囲外（＝外挿）になる、または残り 2 行未満なら、その節点は 0（throw しない）。
   - $s=0$ は厳密に $F=1$ を代入。
2. **s 方向**（`IonizationTableShape`）: 上で得た 321 点を符号付き $F$ のまま PCHIP。$s_\mathrm{cert}<16$ のときは**補間の台を $s_\mathrm{cert}$ までに切る**（埋め草で直下の区間が汚染されないように）。

運動学的上限: $s\le k$（$\lvert\mathbf G\rvert\le 2k$）、Å⁻¹ で $1/\lambda$。30 kV で 14.33 Å⁻¹（コメント値）。

### 2.3 打ち切りと誤差上界（formatVersion ≥ 4 = dataset v5 以降）

- $s\le s_\mathrm{cert}$: 表の補間値。
- $s> s_\mathrm{cert}$: **$F=0$ を返し、上界 $\varepsilon$ を宣言**（`TruncatedBeyondSMax=true`、`TruncationBound=ε`）。外挿はしない。
- 問い合わせ E0 が表の行の間にあるとき、$\varepsilon=\max(\varepsilon_i,\varepsilon_{i+1})$、$s_\mathrm{cert}=\min(\cdot)$（**ε は内挿しない**。内挿した上界は上界でない、`TryGetTailModel`）。
- 合成殻の上界は $\sum_i w_i\varepsilon_i$（三角不等式、`IonizationLTotalShape.TruncationBound`）。
- 旧 formatVersion の指数 tail $F(s_\mathrm{max})e^{-\hat b(s-s_\mathrm{max})}$ のコードは残っているが、現在の同梱データでは通らない（`IsBoundedTail`）。

### 2.4 合成殻（LTotal / MTotal）と線の系列（Kα, Kβ, Lα, Lβ, Mα）

`IonizationDataProvider.Resolve` → `IonizationLTotalShape`（名前は L だが合成全般に使う）。

- 合成殻: $\sigma=\sum_i\sigma_i$、$F=\sum_i w_iF_i$、$w_i=\sigma_i/\sum\sigma$。各 $\sigma_i$ は**各副殻自身の端**で計算（閉じた副殻は 0 で落ちる）。LTotal = L1+L2+L3（j 分離 dataset）。MTotal は**表にある副殻だけ**（Z=30–32 は M4/M5 が無いので M1–M3）。
- 吸収端 `EdgeEnergyKeV` は構成副殻の**テーブルの** $E_\mathrm{th}$ の最小値、過電圧 $U=E_0/E_\mathrm{edge}$。$U<1.2$ は GUI 警告（`LowOvervoltage`）。
- 線の系列（260925Cl、`XrayLineSeries`）: 光子生成断面積 $\sigma_s=\sum_i B_{si}\sigma_i$、形 $F_s=\sum_i(B_{si}\sigma_i/\sigma_s)F_i$。
  $$B_{si}=\sum_{\ell\in s}P_{k(\ell),i}\,\omega_{k(\ell)}R_\ell$$
  $P_{ki}$ は殻 $i$ の一次空孔 1 個が連鎖後に殻 $k$ に作る空孔数で、内側から順に $P_{ki}=\delta_{ki}+\sum_{j<k}t_{kj}P_{ji}$。移送 $t_{kj}$ = 放射遷移 $\omega_jR_{(k\leftarrow j)}$ + Coster–Kronig $f$（L3←L1 は $f_{13}+f'_{13}$）。**Auger 連鎖は含まない**（同梱 xraylib が Auger rates を削った build）。$\omega$, $R$, $f$ は xraylib 4.2.1。
  - 連鎖で移った空孔は**一次の殻の形 $F_i$ を持つ**とする（局所近似、コメント）。
  - 表に無い供給元（Z>50 の K など）は省く（provenance 文字列に記録）。
  - 系列の定義: Kα=KA_LINE、Kβ=KB_LINE、Lα=LA_LINE（L3）、Lβ=xraylib LB_LINE_MACROS の 13 本、Mα=M5N7+M5N6。
- 物理量: 殻チャネルは `IonizationVacanciesGenerated`、線の系列は `XrayPhotonsGenerated`（全方位の生成数。自己吸収・検出効率なし）。

---

## 3. 局所イオン化「ポテンシャル」行列

両経路とも **局所形状因子近似**（一変数 $F(s)$、$s=\lvert\mathbf G\rvert/2$、$\mathbf G=\mathbf g_h-\mathbf g_g(+\mathbf q)$）。bra/ket で別の運動量移送を持つ二変数 MDFF $K(Q_\mathrm{bra},Q_\mathrm{ket})$ は扱わない（`AlchemiModelTier.FullMdff` は `NotSupportedException`）。双極子近似は使っていない — 形の情報は $F(s)$ テーブル側（Temari 生成）にあり、その内訳は Temari 側が正本。

DWF は**差ベクトルの熱平均** $\langle e^{2\pi i(\mathbf g-\mathbf h)\cdot\mathbf u}\rangle=e^{-M(\mathbf g-\mathbf h)}$ を使う。積近似 $e^{-M(g)}e^{-M(h)}$ は禁止（対角が $\sigma e^{-2M}<\sigma$ になり総断面積が壊れる、AlchemiMuBuilder 冒頭）。$\mathbf G=0$ では $e^{-M}=1$ で、**対角は必ず $\sum_a\mathrm{Occ}_a\sigma$**（熱振動で総断面積は消えない）。TDS 吸収の $(1-e^{-M})$ 構造は持ち込まない。

### 3.1 STEM-EDX: `BetheMethod.ComputeIonizationU`

$$U_\mathrm{ion}(\mathbf G)=\frac{k_\mathrm{vac}}{2\pi V_c}\,\sigma\,F(\lvert\mathbf G\rvert/2)\sum_{a\in Z}\mathrm{Occ}_a\,e^{-M_a(\mathbf G)}\,e^{-2\pi i\mathbf G\cdot\mathbf r_a},\qquad \mathbf G=\mathbf q+\mathbf g_i-\mathbf g_j$$

- 和は**原子番号が一致する全サイト**（元素フィルタ）。等価位置は `atoms.Atom` で展開済み。
- 相対論の $\gamma$ 補正は掛けない（σ・F の E0 依存に含まれる、とコメント）。
- 規格化の意味（Claude 追記・要確認）: $k=1/\lambda$ 規約で虚ポテンシャル $U'_0$ による強度の減衰係数は $\mu_\mathrm{abs}=2\pi U'_0/k$。これを $n\sigma$（$n$ = 単位体積あたりの原子数）と置くと $U'_0=k\,n\sigma/(2\pi)$ で、1 原子あたり $k\sigma/(2\pi V_c)$ になる。後段の係数 $2\pi/k_\mathrm{vac}$（§5）と相殺して $\sigma/V_c$（単位長あたりのイオン化確率密度）になる。
- `FillIonizationUq`: 1 つの $\mathbf q$ について `dest[off + j*bLen + i] = U(q + g_i − g_j)`（column-major）。チャネル内限定の combined-hkl キャッシュ（`ConcurrentDictionary`）。値は正準 $\mathbf G$ から決定的に計算されるので勝者スレッドに依らない。

### 3.2 ALCHEMI: `AlchemiMuBuilder.Build`

$$\mu_{hg}=\sigma\,F(\lvert\mathbf G\rvert/2)\sum_{a\in\text{site}}\mathrm{Occ}_a\,e^{-M_a(\mathbf G)}\,e^{-2\pi i\mathbf G\cdot\mathbf r_a},\qquad\mathbf G=\mathbf g_h-\mathbf g_g\quad[\mathrm{nm}^2]$$

STEM-EDX との意図的な違い（ファイル冒頭コメント）:
- **元素フィルタを掛けない**: サイト仮説（`AlchemiSiteBasis.AtomsIndices`）の原子軌道の幾何（$\mathbf r_a$、ADP、Occ）× ドーパントのチャネル（σ, F）。tracer 仮説の組み合わせ。
- $k_\mathrm{vac}/(2\pi V_c)$ を掛けない（「断面積行列」、単位 nm²）。$V_c$ は収率側で割る。
- 占有率は `OccupancyFraction` で上書き可（Tracer 基底では 1.0 を渡し、実占有率は縮約後に線形結合で掛ける）。
- 変位は `SiteDisplacement`: `InheritHost`（既定、ホストの `Dsf`）/ `IsotropicB`（$B_\mathrm{iso}$ nm²）/ `AnisotropicB`（$\beta$ 無次元）。
- `AlchemiModelTier.StaticDelta`: $F\equiv1$, DWF$\equiv1$ の内部 oracle（縮約経路の検証用、GUI 非公開）。

数値上の工夫: $\Delta\mathbf G$ は基底だけで決まるので、構築時に unique 化して pair→slot 表を持つ（$n=200$ で $n^2=40000$ ペアに対し unique は数千）。$F(s)$ は unique 側でバッチ評価、全（サイト×チャネル）で表を共有。μ は方位に依存しない（位相は整数 hkl と分率座標、$\lvert\mathbf G\rvert$ は回転不変）ので、FixedUnion 基底の run では（サイト×チャネル）ごとに **1 回だけ**組む。

対角値: `Mu00` $=\mu_{00}=\sigma\sum_{a}\mathrm{Occ}_a\cdot m_a$（$m_a$ は軌道の多重度 `atoms.Atom.Length`）。

---

## 4. Bloch 波からのイオン化収率（ALCHEMI、`AlchemiReduction`）

### 4.1 波動場と局所イオン化率

固有系 $(\gamma_j, C^{(j)}_g)$ は `getEigenMatrix`（行を $P_g$ で割った線形化 Bethe 方程式、260829Cl 修正後）の固有分解。境界条件 $\psi_g(0)=\delta_{g0}$ から $\boldsymbol\alpha=C^{-1}\mathbf e_0$（LU）。

$$\psi_g(z)=\sum_j\alpha_jC^{(j)}_g e^{2\pi i\gamma_j z},\qquad \frac{dP}{dz}=\psi^\dagger(z)\,\mu\,\psi(z)$$

$D=C\,\mathrm{diag}(\alpha)$、$S=D^\dagger\mu D$（$S_{j'j}=\sum_{h,g}D^*_{hj'}\mu_{hg}D_{gj}$）とおくと厚み積分は解析的:

$$P(t)=\mathrm{Re}\sum_{j,j'}S_{j'j}F_{jj'}(t),\qquad F_{jj'}(t)=\frac{e^{\lambda_{jj'}t}-1}{\lambda_{jj'}},\quad\lambda_{jj'}=2\pi i(\gamma_j-\gamma^*_{j'})$$

- $\lvert\lambda\rvert^2<10^{-30}$（nm⁻²）は縮退とみなし $F=t$（ロピタル）。
- $e^{\lambda_{jj'}t}=u_j(t)u^*_{j'}(t)$、$u_j=e^{2\pi i\gamma_jt}$ と分解し、指数関数は厚みあたり $b$ 回（$b^2$ 回でなく）。
- **負値を clamp しない**（打切り・符号・規格化の誤りを潰さない）。
- StaticDelta 極限では $\psi^\dagger\mu\psi=\sigma\sum_n\lvert\psi(\mathbf r_n)\rvert^2$ となり、EBSD の `EBSDSolverManaged` の $S=B^\dagger\mathrm{diag}(\sigma)B$ と同じ量（回帰 oracle）。

### 4.2 動力学項と非チャネリング項

$$Y_\mathrm{dyn}(t)=\frac{P(t)}{V_c},\qquad Y_\mathrm{dech}(t)=\frac{\mu_{00}}{V_c}\big\{t-L_\mathrm{coh}(t)\big\},\qquad Y=Y_\mathrm{dyn}+Y_\mathrm{dech}$$

$$L_\mathrm{coh}(t)=\int_0^t\sum_g\frac{P_g}{P_0}\lvert\psi_g(z)\rvert^2dz,\qquad P_g=2\hat{\mathbf n}\cdot(\mathbf k_0+\mathbf g)$$

- $L_\mathrm{coh}$ は $\mu=\mathrm{diag}(P_g/P_0)$ の同じ縮約で計算し、方位ごとに 1 回だけ作ってキャッシュ（`CoherentPathLengthNm`）。重み $P_g/P_0$ は 260829Cl 追加（法線方向全電流＝無吸収で厳密に 1）。`pOverP0=null` なら旧来の $\sum\lvert\psi_g\rvert^2$。
- 意味: TDS 吸収で coherent 場から抜けた電子はエネルギーを保ったまま試料内に残り、方向がランダム化された電子として残り厚みで random 方位の断面積 $\mu_{00}/V_c$ でイオン化を続ける（ICSC 2003 の式(1)、コメント）。
- 単位: $Y$ は**入射電子 1 個あたりの発生イオン化数（無次元）**（`SignalNormalization.PerIncidentElectron`）。
- ガード: $t-L_\mathrm{coh}<-10^{-6}\max(t,10^{-12})$ なら `InvalidOperationException`（clamp せず fail）。それより小さい負の丸めだけ 0 にする。
- **再注入の契約（`AbsorptionSource`）**: $t-L_\mathrm{coh}$ は虚部が作る減衰を出所を問わず拾う。再注入が正当なのは虚部が TDS のみのとき。`BetheMethod.ImaginaryPotentialAbsorption` が `TdsRedistributable` を宣言し、`Yield` は `IsFullyRedistributable()`（TDS 以外のビットが立っていない）でなければ例外。`TrueLoss` / `Phenomenological` を混ぜたら、出所ごとに分けて TDS 部分だけで再注入する設計を要する。宣言自体は人間の規律（虚部の作り方を変えたら宣言を更新）。
- 結果は Dynamic / Dechannelled / Total を**分離して保存**（$t-L_\mathrm{coh}$ は厚み・方位依存なので fit の定数 nuisance では吸収できない）。

### 4.3 run の流れ（`BetheMethod.RunAlchemi`、`BetheMethod.Alchemi.cs`）

1. チャネル解決（`IonizationDataProvider.Resolve`、run 中は immutable）。`ChannelDataOverride` は感度試験専用（F に δF を載せて差を取ると、Bloch 解が F に依らないので δF の寄与だけが厳密に出る）。
2. **FixedUnion 基底**: 走査の全方位で `Find_gVectors` を呼び hkl の真の union を取る（基準方位 1 回だけでは落ちた g の誤差が方位相関を持ち、サイト固有 ICP と誤認される）。順序は $\lvert g\rvert^2\to h\to k\to l$ で 000 が先頭。`BasisOrientationIndices` で基底方位を限定可能（文献の固定基底の再現用、診断専用）。
3. μ を（サイト×チャネル）ぶん 1 回だけ構築。
4. 方位ループ（`Parallel.For`）: `reset_gVectors` → `getEigenMatrix`（ポテンシャル行列 `getPotentialMatrix(union)` は全方位共通、局所形式）→ EVD（native Eigen、または MathNet `Evd(Asymmetric)`）→ α（LU）→ `AlchemiReduction.Yield`。
5. expanded-basis 診断: 中心・両端の 3 方位を `MaxNumOfBloch × ExpandedBasisFactor`（既定 1.25）の union で解き直し、Total の最大差 / テンソル全体の最大値。閾値 `ExpandedBasisTolerance` 既定 3×10⁻³。

既定値（`AlchemiRequest`）: `Surface=(0,0,−1)`、`ModelTier=LocalFormFactor`、`MaxNumOfBloch=300`、`IncludeDechannelledComponent=true`、`ExpandedBasisFactor=1.25`、`ExpandedBasisTolerance=3e-3`、`UseNativeSolver=true`。方位はビーム方向で表す（BaseRotation を変えると gCache 再構築で遅い）。`AlchemiScan.Rocking1D` は傾斜軸の $\hat d_0$ 平行成分を除いた軸まわりに等間隔。

警告（`BuildFixedUnion`）: 基底が要求する $s_\mathrm{max}=\max\lvert g\rvert$（Å⁻¹）が 16 を超える / 傾斜幅 >30 mrad（v1 FixedUnion 保証外）、10–30 mrad（union + expanded-basis 必須）/ union が中心のみ基底の 4 倍超 / $P\le0$ のビームあり。

`AlchemiBasisDiagnostic.Eligibility` は**常に `NotEvaluated`**（260809Cl 作者決定）。理由はコード中に列挙: ⑤ 分母がテンソル全体の最大値で σ の小さいチャネルが薄まる、⑥ 分子が絶対収率で ICP では落ちる共通スケールも数える、⑦ Find_gVectors がエワルド球スクリーニング律速で 1.25 倍でも基底が増えない領域では自己比較になり偽陽性。

### 4.4 角度広がり（`AlchemiAngularSpread.Gaussian`）

方位軸上の後処理（engine は触らない）。**畳み込み → 規格化**の順に固定。
$$\tilde Y_i=\frac{\sum_j w_{ij}Y_j}{\sum_j w_{ij}},\quad w_{ij}=\exp\!\left[-4\ln2\,\frac{(\theta_i-\theta_j)^2}{\mathrm{FWHM}^2}\right],\quad \lvert\theta_i-\theta_j\rvert\le1.7\,\mathrm{FWHM}$$
端部はカーネルを再規格化（走査外を 0 とみなすと端に偽の谷が出る）。非等間隔でも動く。FWHM ≤ 0 または点数 < 2 は入力をそのまま返す。

---

## 5. STEM-EDX の像（`BetheMethod.StemDoWork` 内）

STEM の非弾性（TDS）経路と**同じ**スライス積分を、$U$ を $U_\mathrm{ion}$ に替えて回す（`AccumulateInelasticQ` に集約、260802Cl）。

- q 空間量: 各 $\mathbf q$（逆格子点）・厚み $t$・デフォーカス $d$ について、
  $I(\mathbf q,t,d)=\sum_{\text{slice}\le t}\frac{2\pi}{k_\mathrm{vac}}\Delta z\cdot\big[\text{native }\texttt{STEM\_InelasticQ}\text{ の縮約}\big]$。
  縮約の中身は「$W_q=U_q\times TC$ を GEMM → エントリ（入射方向 $k$、近傍 4 点の bilinear 補間、レンズ関数）ごとに 4 dot → $\sum_d$」（コメント）。係数 $2\pi/k_\mathrm{vac}$ が $U_\mathrm{ion}$ の $k_\mathrm{vac}/2\pi V_c$ と相殺し、$\sigma/V_c$ の単位長あたり確率になる。
- ±q 処理（設計書 §3.4 の順序）: 独立計算 → 残差 $\max\lvert I(\mathbf q)-I(-\mathbf q)^*\rvert/\max\lvert I\rvert$ と $q=0$ 虚部を記録 → 許容超過は**対称化せず fail** → 合格時のみ $\bar I(\mathbf q)=(I(\mathbf q)+I(-\mathbf q)^*)/2$、$q=0$ は実部固定。$I(q=0)$ の実部が $-10^{-9}\max\lvert I\rvert$ 未満なら fail。
  - 許容値 `HermitianTolerance` 既定 0.01、上限 0.01（厳しくする方向のみ）。残差は方向グリッド div に対し概ね $O(h^2)$（div 10/20/32/48 で 0.107/0.018/0.0093/0.0017 実測）、推奨 `RecommendedProbeDivision=48`。
  - 環境変数 `RECIPRO_STEM_EDX_HERMTOL` で上限を外せる（収差検討用の escape hatch、260805Cl）。
- 実空間: $\mathrm{image}(\mathbf R)=\mathrm{Re}\big[\sum_\mathbf q I(\mathbf q)e^{2\pi i\mathbf q\cdot\mathbf R}\big]/N$、$N$ は **$q=0$ の qEntry 数＝絞り内の入射方向数**（260802Cl 修正。旧 `radiusPix²` は約 2.7 倍過大かつ div 依存だった）。`.Magnitude` でなく実部（線形性と絶対規格化を保つ）。ソースサイズはガウスぼかし。最後に負値を 0 に clamp（clamp 前最小値は `MinPixelBeforeClamp`）。
- 制約: native ヘルパ必須（managed 経路は未実装 → `NotSupportedException`）。同一 (Z, Shell) の重複要求は hard error。
- 診断の公開: `StemSignalMap.HermitianResidualMax / QZeroImagMax / MinPixelBeforeClamp / TruncatedBeyondSMax / TruncationBound / Iq`（対称化後）/ `IqBeforeSymmetrization`（`CaptureRawIq` 時）。

---

## 6. 近似と適用範囲（まとめ）

- 孤立原子（DHFS/Dirac）の F(s)、局所（一変数）形状因子近似。結合・固体効果なし。
- σ は Bote–Salvat（不確かさは README）。E0 = 30–400 keV。
- 連鎖空孔は一次の殻の形を持つ。Auger 連鎖なし。
- ALCHEMI の非チャネリング電子は「エネルギー保存・方向完全ランダム・random 方位断面積でイオン化」。TDS 以外の吸収は再注入不可。
- 自己吸収・検出器効率（`XrayPhotonsSelfAbsorbed/Detected`）は enum だけ用意され未実装。

---

## 実装上の疑問点

1. **角度広がりの打ち切り幅のコメント誤り**（`Diffraction/Alchemi/AlchemiAngularSpread.cs:48–50`）。「4σ (≈1.7 FWHM) より外は重み <1e-7」とあるが、$w(1.7\,\mathrm{FWHM})=\exp(-4\ln2\cdot1.7^2)=e^{-8.01}\approx3.3\times10^{-4}$。打ち切られる片側の質量は約 $3\times10^{-5}$（全体比、Claude 計算）。数値への影響は小さいが、コメントの根拠どおりの精度（1e-7）を狙うなら $\approx2.4$ FWHM（$\approx5.7\sigma$）が必要。
2. **$L_\mathrm{coh}$ は電流重み、$Y_\mathrm{dyn}$ は密度で数えている**（`AlchemiReduction.cs:146–171` と `:194`、260829Cl）。$\mu=\mu_{00}I$ の極限で $Y=\frac{\mu_{00}}{V_c}\big[t+\int\sum_g(1-P_g/P_0)\lvert\psi_g\rvert^2dz\big]\ne\mu_{00}t/V_c$。表面法線 ⊥ $\mathbf g$（ZOLZ）なら $P_g=P_0$ で差は 0、HOLZ や傾斜表面でのみ現れる。物理的にどちらが正しいか（斜行成分は深さあたりの行路が長い）を含めて要判断。
3. **$Y_\mathrm{dyn}$ に $2K/P_0$（≈ $1/\cos\theta_0$）の因子が無い**（`AlchemiReduction.Yield`、`:196`）。入射電流は $\propto P_0\lvert\psi_0\rvert^2$ なので、入射電子 1 個あたりの深さ積分イオン化は $(n\sigma)(2K/P_0)\int\lvert\psi\rvert^2dz$ になるはず（Claude の導出・要確認）。$Y_\mathrm{dech}$ も行路長 $t$（法線方向）で数えている。走査幅 30 mrad では $\sim4.5\times10^{-4}$ と小さいが、表面法線とビームが大きく傾く条件では効く。
4. **`IonizationTableShape` の診断フラグが共有の可変状態**（`IonizationChannel.cs` `IonizationTableShape._truncated/_usedTail`）。STEM の `Parallel.For` から同じ shape を並行評価して書く（true を書くだけなので結果は良性）が、`IonizationData` を run 間で使い回す（例: `ChannelDataOverride`）と前の run のフラグが残り、`TruncatedBeyondSMax` が誤って true になり得る。
5. **E0 方向 PCHIP の端の外挿**（`IonizationChannelTable.GridAt`）。行を落とさない節点（$n=$rows）では、問い合わせ E0 が表の最初/最後の行の外でも `Pchip.Evaluate` は端区間の 3 次式で外挿する（範囲チェックは 30–400 keV のみ）。表の E0 行が 30 と 400 keV を必ず含むなら問題なし — 含むかは `.bin` を読まないと確認できない。
6. **`RECIPRO_STEM_EDX_HERMTOL` に `NaN` を与えると検査が素通り**（`BetheMethod.cs:3450–3452, 3481`）。`double.TryParse("NaN")` は成功し、`hermMax > NaN` は常に false。API 側の `HermitianTolerance` は NaN を run 前に拒否しているので、それと揃っていない。
7. **古い数値が残るメッセージ/コメント**: `IonizationTableShape.Evaluate` の例外文「s > 4」、クラス要約「s_max (=8 Å⁻¹)」、`AlchemiBasisDiagnostic` の「s ≤ 8 Å⁻¹」、`BuildFixedUnion:186–189` の「テーブル上限 (s ≤ 8 Å⁻¹)」、`IonizationFsTable.SCount` の注記「81」。現行は 321 点・16 Å⁻¹。
8. **`AlchemiMuBuilder.cs:206` のコメント**は「Dsf.BisoEffective と同じ規則」とするが実装は `dwf.Biso000`（`ComputeIonizationU` は `BisoEffective`）。$\mathbf G=0$ では $s^2=0$ なので値は一致し実害は無い（261105Cl の変更）。

## 要確認事項

1. STEM-EDX の native カーネル `STEM_InelasticQ`（Crystallography.Native、本リポジトリ外）の式。本稿は C# 側のコメントからの再構成。
2. $U_\mathrm{ion}$ の規格化 $k_\mathrm{vac}/(2\pi V_c)$ の物理的導出（§3.1 の説明は Claude 追記）。
3. 260802Cl コメントの実測「$I(q=0)$/entries/$(\sigma tN/V)$ = 1.00586」の 0.6 % 余剰の由来（260829Cl の $P_g$ 行割り修正以前の測定か）。
4. Coster–Kronig で L3←L1 に $f_{13}+f'_{13}$ を足す扱いが、xraylib の $f'_{13}$ の定義（Claude 追記・要確認: Krause の放射を含む L1→L3 遷移か）と二重計上にならないか。
5. LTotal/MTotal の `EdgeEnergyKeV` はテーブルの $E_\mathrm{th}$、σ は Bote の端を使う。両者の値の差の大きさ。
6. F テーブルの E0 行の配置（30 / 400 keV を端点に含むか、`IonizationFsE0.bin` の行数）。疑問点 5 と関連。
7. 疑問点 2・3（流束規格化と $1/\cos\theta_0$）の扱いを設計書（`.project-guidance/ReciPro/ReciPro_ALCHEMI設計.md` §3.4、本リポジトリ外）がどう定めているか。
8. `XrayLineSeries` の Kβ=KB_LINE「= 1 − Kα」が xraylib の RadRate の定義どおりか。
9. TDS 経路の `getU(...).Imag.Conjugate()`（「共役とると、なぜかいい感じ」）と、イオン化経路の「負位相直生成」が同じ位相規約になっていることの確認（[bloch-wave.md](bloch-wave.md) 側の検証と合わせる）。
