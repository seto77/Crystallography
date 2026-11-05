本稿はコードから Claude が起こした下書き（2026-11-05）。著者未査読。

# 原子散乱因子・温度因子・TDS 吸収形状因子

対象: `Atom/AtomStatic.cs`（`AtomStatic.ES`, Temari 因子, 中性子散乱長）、`Atom/Atoms.cs`（`Atoms.GetAtomicScatteringFactorFor*`, `DiffuseScatteringFactor`）、`Xraylib.cs`（f′/f″）。
構造因子への組み込みは [structure-factor.md](structure-factor.md)、Bloch 波ポテンシャル $U_g$ への組み込みは [bloch-wave.md](bloch-wave.md) を参照。

---

## 0. 共通の約束事（最重要）

| 量 | 記号 | 実装上の単位・定義 |
|---|---|---|
| 散乱ベクトルの大きさ | $s=\sin\theta/\lambda=|g|/2=1/(2d)$ | **$q=4\pi s$ は使わない**（Kirkland 表だけ $q_K=2s$。後述）|
| `ES.Factor(S2)` の引数 | $s^2$ | **nm⁻²**（呼び出し側は `vec.Length2/4`, `1/d²/4`）|
| `ES.Factor` の戻り値 | — | **電子用 = nm、X 線用 = electrons × 0.1**（260818Cl 訂正）|
| `ES.ElasticFactorWithTail(s2)` | $s^2$ | **Å⁻²**、戻り値 Å（同じクラス内で規約が異なる）|
| `TemariFactor.Fx(s)`, `Fe(s)` | $s$ | **Å⁻¹**、戻り値 electrons / Å |
| 温度因子 | $B$ | **nm²**。`DiffuseScatteringFactor` 内部は B 型と U 型を相互変換（$B=8\pi^2U$）|
| 異方性 | $\beta_{ij}$ | `B11`…`B31` は**無次元の β**（$T=\exp(-\mathbf h^{\mathsf T}\beta\mathbf h)$）|
| 電子散乱因子の γ | — | **すべての表（Peng / 8-Gaussian / Kirkland / Temari）は γ を含まない**。γ は `BetheMethod.getU` で掛ける |
| 角度 | — | 内部はラジアン（検出器 inner/outer もラジアン）|

### 換算の実例（`ES` コンストラクタ, AtomStatic.cs L2779–2822）
Gauss 和の表パラメータは Å 系（$a_i$ [Å] または electrons, $b_i$ [Å²]）。`Factor` は

$$
\texttt{Factor}(S_2)=0.1\times\Big[\sum_i a_i\exp\!\big(-b_i\cdot 0.01\,S_2\big)\;(+c)\Big],\qquad S_2\,[\mathrm{nm^{-2}}]
$$

で、$0.01$ は nm⁻² → Å⁻²、$0.1$ は Å → nm。X 線（定数項 $c$ あり）にも同じ ×0.1 が掛かるため、X 線の戻り値は「electrons×0.1」という人工的な単位になる。`FormBeamInteraction` は X 線のとき ×10 して electrons に戻す（コメント L2734–2737）。

---

## 1. 弾性散乱因子の表（`AtomStatic` 静的フィールド）

| フィールド | 形 | 変数 | 出典（コード内記載） | ライブラリ内の利用 |
|---|---|---|---|---|
| `XrayScatteringWK[Z][sub]` | 5 Gauss + c | $s$ [Å⁻¹] | Waasmaier & Kirfel (1995) Acta Cryst. A51, 416–431 | **X 線の運動学 F**（`Atoms.GetAtomicScatteringFactorForXray`）、Mott–Bethe テールの $f_x$ |
| `XrayScattering[Z][sub]` | 4 Gauss + c | $s$ | Cromer & Mann (1968) Acta Cryst. A24, 321–324 | 本ライブラリ内では未使用（参照箇所なし）|
| `ElectronScatteringPeng[Z][sub]` | 5 Gauss | $s$ | Peng et al. (1996) Acta Cryst. A52, 257–276（中性）/ Peng (1998) A54, 481–485（イオン）| **電子の弾性 $U_g$ の実部**、運動学 F、TDS 吸収の被積分関数 |
| `ElectronScatteringEightGaussian[Z]` | 8 Gauss | $s$ | Ohtsuka 氏の私信（公表文献なし）| 未使用（`getU` 内にコメントアウト）。`AttachTail` だけは呼ばれる |
| `ElectronScatteringKirkland[Z]` | 3 Lorentz + 3 Gauss | $q_K=2s$ [Å⁻¹] | Kirkland, *Advanced Computing in Electron Microscopy*, App. C（fparams.dat）| `ProjectedPotential` のみ、廃止予定の `DynamicTheory.MultiSlice` から |

`[Z][0]` は常に中性原子のエントリ（`AtomStatic` 静的コンストラクタと `getU` の `ElasticIonModel.Neutral` 分岐の共通規約）。`sub ≥ 1` はイオン（`ES.Valence` に価数）。

### 1.1 Peng 5-Gaussian（電子）

$$
f_e(s)=\sum_{i=1}^{5}a_i\,e^{-b_i s^2}\quad[\text{Å}],\qquad s\ [\text{Å}^{-1}]
$$

- 260923Cl の注記訂正: 収録係数は ITC Vol. C Table 4.3.2.2 = Peng 1996 の **$s\le 2$ Å⁻¹ 用**の当てはめ（$s\le 6$ 用は別表 4.3.2.3）。$s>2$ は外挿。
- 260811Cl の実測（`AlchemiCheck absorb -join`）: Kirkland 比で $s\approx 2.1$ Å⁻¹ から 5 % 以上外れ、$s=4$ で 0.09、$s=8$ で $10^{-7}$。Gauss 和は高 $s$ で厳密に 0 へ落ち、Mott–Bethe の $1/s^2$ テールを持たない。
- **弾性 $U_g$ の実部は純 Peng のまま**（テールなし）。理由はコメントに「実用反射は $|g|\le 50$ nm⁻¹（$s\approx0.25$ Å⁻¹）で接続域に届かず、変えると全結果が動く」と記載（L2918–2921）。

### 1.2 Kirkland（電子、3 Lorentz + 3 Gauss）

$$
f_e(q_K)=\sum_{i=1}^{3}\frac{a_i}{q_K^2+b_i}+\sum_{i=1}^{3}c_i e^{-d_i q_K^2},\qquad q_K=2s=1/d\ [\text{Å}^{-1}]
$$

実装（L2830–2841）は $q_K^2[\text{Å}^{-2}]=0.04\,S_2[\text{nm}^{-2}]$ と換算し、×0.1 で nm に直す（260606Cl 修正。以前は換算をすべて欠いていた）。投影ポテンシャル（Kirkland 式 C.20）:

$$
v_z(r)=2\pi^2 a_0 e\Big[\sum_i 2a_i K_0\!\big(2\pi r\sqrt{b_i}\big)+\sum_i\frac{c_i}{d_i}e^{-\pi^2r^2/d_i}\Big]\quad[\text{V·Å}],\ r\ [\text{Å}]
$$

$a_0e$ は `UniversalConstants.a0·1e10 × 2Ry·a0`（= 0.529 Å × 14.4 V·Å）。$r\le0$ は 0 を返す。

### 1.3 X 線（Waasmaier–Kirfel / Cromer–Mann）

$$
f_x(s)=\sum_i a_i e^{-b_i s^2}+c\quad[\text{electrons}]
$$

コンストラクタの引数並びは $(a_1,b_1,a_2,b_2,\dots,c,\text{valence},\text{method})$ の交互形。元論文の $(a_1..a_5,c,b_1..b_5)$ 並びから並べ替えて格納されている（H の値で確認）。

### 1.4 X 線異常分散 f′, f″（`Xraylib.cs`）

- `Xraylib.Fprime(Z, E_keV)` = xraylib `Fi`、`Xraylib.Fdoubleprime(Z, E_keV)` = **−`Fii`**（xraylib の生値は符号が逆。結晶学慣用の $f''>0$ に揃える）。E は keV、戻り値 electrons。
- 自己テスト: Fe @ 8.04778 keV で $f'=-1.106$, $f''=+3.192$（ITC の $-1.13$, $3.20$ と照合、Xraylib.cs L24）。
- 有効条件は `Crystal.IsXrayDispersionActive`: `anomalousDispersion && wave==Xray && E 有限 > 0 && Xraylib.Enabled`。**元素（Z）とエネルギーだけに依存し、イオン価数・化学状態は無視**。
- 構造因子への入れ方（符号規約を含む）は [structure-factor.md §2](structure-factor.md) 参照。

### 1.5 中性子散乱長（`NeutronCoherentScattering[Z][i]`）

- 単位 **fm**、複素数。`[Z][0]` = 天然組成、`[Z][i≥1]` = 同位体（`IsotopeAbundance[Z]` の列挙順）。出典: periodictable 2.1.0 nsf（Sears 1992, Neutron Data Booklet 系、260531Ch）。
- 虚部は 1.798 Å 基準の吸収: $\mathrm{Im}\,b=-\sigma_a/(2000\,\lambda_{th})$（b [fm], σ [barn], λ [Å]）。`NeutronAbsorptionCrossSection(Z, λ)` は $\sigma_a(\lambda)=-\mathrm{Im}\,b\cdot2000\,\lambda$（1/v 則）で復元。Cd/Sm/Eu/Gd は 1/v が破れる（`NeutronIsResonantAbsorber`）。
- $\sigma_{coh}$ は表値を使い、$4\pi|b|^2$ から計算しない（共鳴核で乖離、Gd: 35.3 vs 29.3 barn）。
- 同位体指定時（`Atoms.GetAtomicScatteringFactorForNeutron`）: $b=\sum_i b_i\,x_i/100$（$x_i$ は % 単位の `Isotope[i]`）× Occ。配列長が `IsotopeAbundance[Z].Count` と一致しないときは天然値。$s$ 非依存。

---

## 2. `Atoms` からの散乱因子（占有率込み）

`Atoms.cs` L355–375:

| メソッド | 式 | 備考 |
|---|---|---|
| `GetAtomicScatteringFactorForElectron(s2)` | `ElectronScatteringPeng[Z][SubNumberElectron].Factor(s2) * Occ` | s2 [nm⁻²]、戻り nm。**割当イオン（SubNumberElectron）をそのまま使い `BetheMethod.ElasticIonModel` は見ない** |
| `GetAtomicScatteringFactorForXray(s2)` | `XrayScatteringWK[Z][SubNumberXray].Factor(s2) * Occ` | 戻り electrons×0.1 |
| `GetAtomicScatteringFactorForNeutron()` | 上記 §1.5 × Occ | fm |

**温度因子はここでは掛けない**（構造因子・$U_g$ 側で掛ける）。

---

## 3. 温度因子（`DiffuseScatteringFactor`, Atoms.cs L442–）

内部保持は `OriginalType ∈ {U, B}`、`Iso`、`Aniso11..31`（U 型: nm²、B 型: 無次元 β）。取得側:

$$
B_{iso}=\begin{cases}\mathrm{Iso} & (B\ 型)\\ 8\pi^2\,\mathrm{Iso} & (U\ 型)\end{cases}\quad[\mathrm{nm^2}],\qquad
\beta_{ij}=\begin{cases}\mathrm{Aniso}_{ij} & (B\ 型)\\ 2\pi^2 a^*_i a^*_j\,U_{ij} & (U\ 型)\end{cases}
$$

- 261105Cl 修正: 非対角の係数は $2\pi^2a^*_ia^*_j$（**$\cos\gamma^*$ 等を掛けない**。CIF `_atom_site_aniso_U_ij`、Grosse-Kunstleve & Adams (2002) 式 4–6）。以前は cos を掛けていた（直交格子で $U_{12}$ 等が消え、六方晶で $U_{12}$ が半分）。
- 等価等方 $B$（Hamilton 1959 式 18、261105Cl で実格子の内積に修正）:
  $$B_{eq}=\tfrac43\sum_{ij}\beta_{ij}\,(\mathbf a_i\cdot\mathbf a_j)$$
  `Biso000`。U 型からは $B_{eq}=8\pi^2U_{eq}$ に一致。
- `BisoEffective`（g = 0 用, 261105Cl 変更）: `IsZero`→0、`UseIso`→`Biso`、非等方→常に `Biso000`。NaN・負は 0。
- `BisoFloorNm2 = 1e-3` nm²（0.1 Å²）: EBSD 後方散乱源と MC 用の下限。**`getU` の吸収には適用しない**。
- 係数は `Dsf.Cell`（nm, rad）から作る。`Cell` は `[XmlIgnore]` なので読込後に `Crystal` が設定し直す（Crystal.cs L711–719）。
- 非等方の原子 $j$（等価位置 $\mathbf r_j=R\mathbf r+\mathbf t$）での温度因子は、`SymmetryOperation.ConvertPlaneIndex` が返す $R^{\mathsf T}\mathbf h$ を代入する（261105Cl の訂正コメント, SymmetryOperation.cs L99）:
  $$T_j(\mathbf h)=\exp\!\big[-(R^{\mathsf T}\mathbf h)^{\mathsf T}\beta\,(R^{\mathsf T}\mathbf h)\big]$$
  これは $\beta'=R\beta R^{\mathsf T}$ と等価。
- `Transform`（β′ = MβMᵀ）と `Deform`（格子変形で Cartesian 楕円体を保つ、M = A⁻¹S⁻¹A）は 261105Cl 追加。誤差は 6 成分独立として伝播（共分散なし）。

### 等方温度因子の掛け方
構造因子・$U_g$ とも $T=\exp(-B s^2)$（$B$ [nm²], $s^2$ [nm⁻²]）。**$\exp(-2B s^2)$ や $\exp(-8\pi^2U s^2)$ の取り違えはない**ことを確認した（Crystal.cs L1590・L1673, BetheMethod.cs L3921）。

---

## 4. 電子散乱因子の高 s テール（`ES.ElasticFactorWithTail`, L2957）

**吸収（TDS）の被積分関数だけ**に使う弾性散乱因子。引数 $s^2$ [Å⁻²]、戻り Å。グローバル設定 `BetheMethod.AbsorptionElasticModel`（既定 `TemariHybrid`, 260922Cl 作者決定）で分岐:

1. **TemariHybrid かつ中性（Valence = 0）かつ $s\le 6$ Å⁻¹ かつ Z = 1–86** → `TemariScattering(Z).Fe(s)`。
2. それ以外（イオン、Z = 87–98、`PengTail`、$s>6$）:
   - $s\le 1.5$: 純 Peng（`PengSum`）
   - $s\ge 2.5$: Mott–Bethe
     $$f_e^{MB}(s)=0.023934\,\frac{Z-f_x(s)}{s^2}\ [\text{Å}],\qquad f_x\in[0,Z]\ \text{へクランプ}$$
     $f_x$ は**同元素の中性原子**の WK（イオンでも中性を使う。高 $s$ では核電荷だけが見えるため）。係数 `MottBetheCoefficientAngstrom = 0.2393366/10`。
   - $1.5<s<2.5$: smoothstep $w=t^2(3-2t)$, $t=(s-1.5)/1.0$ で $f=(1-w)f_{Peng}+wf^{MB}$。
3. `useTail=false` は検証専用（純 Peng）。

接続の根拠（260811Cl コメントの表）: Mott–Bethe は Kirkland と全域で 0.97–1.00 の比、Peng は $s\approx2.1$ で 5 % 外れる。クランプの根拠: WK の当てはめ外で $f_x<0$ に発散する元素（C: $a_4=-4.24$, N: $c=-11.80$）。Temari 採用の根拠（260922Cl）: $s>6$ の寄与は $f''(0)$ の 1.5–4.8 %、テール模型差は 0.04–0.43 %、観測量への効き ≤ 0.45 %。$s=6$ での値の跳びは ≤ 0.041 % で、合わせていない。

テールを入れた経路: `FactorImaginary`, `FactorImaginaryAnnular`（局所・非局所）。入れていない経路: `FactorImaginaryAnnularFlatEwald`, `FactorImaginaryAnnular2`（いずれも本番未使用・未完成扱い）。

`AttachTail(z, neutralXray)` は `AtomStatic` 静的コンストラクタが Peng 表と 8-Gaussian 表の全エントリに対して呼ぶ（L3335–3357）。`neutralXray = XrayScatteringWK[z][0]`。X 線表・Kirkland 表では `AtomicNumber = 0` のまま（テール無効）。

---

## 5. Temari 因子（`AtomStatic.TemariFactor`, L7569–）

- 出典（コード記載）: Temari dataset-factors v2.0.0（Y. Seto）, CC-BY-4.0, DOI 10.5281/zenodo.22820415。完全 Dirac SCF + 厳密交換 KLI（model_id `DHFS-KLI-DTM1-dt16-neutral-v1`）。データセットは**誤差上界を主張しない**（`certification_status = not_certified`）。
- **中性原子のみ、Z = 1–86、$0\le s\le 6$ Å⁻¹**。範囲外は例外（外挿もクランプもしない、`Guard`）。
- `Fx(s)` [electrons]: $s$ 上の C² 3 次スプライン、左端 clamped $f_x'(0)=0$、右端 not-a-knot。
- `Fe(s)` [Å]: **$t=s^2$ 上**のスプライン、両端 not-a-knot。**γ を含まない**（`TemariGammaIncluded = false`）。first-Born（Mott–Bethe）の非相対論的 $f_e$。
- 節点 7681 点（$s_i=6\,i/7680$）。格納は絶対量子化（$f_x$ 1e-10 electrons, $f_e$ 1e-11 Å）+ 3 階差分 + byte-plane shuffle + Brotli（`TemariFactors.bin`, magic "TMRF", format 1）。公開値との差は最大 5e-11 / 5e-12。
- `TemariSpline` は参照ローダ `tools/temari_factors_contract.py` の消去順をそのまま移植（式を整理しないこと、とコメント）。区間探索は 260922Cl で等間隔見当 + 補正に変更（値はビット同一）。
- **ライブラリ内での利用は `ElasticFactorWithTail` の `Fe` のみ**。`Fx` は弾性 X 線構造因子には使われていない（X 線 F は WK のまま）。

---

## 6. TDS 吸収形状因子 $f''$（`ES.FactorImaginary*`）

吸収は**ポテンシャルの虚部**として入る。$U_g=U^{(r)}_g+iU'_g$（`BetheMethod.getU` は `(Real, Imag)` の組で返し、それぞれ複素数）。虚部の出所は現状 TDS のみ（`BetheMethod.ImaginaryPotentialAbsorption = AbsorptionSource.TdsRedistributable`）。

### 6.1 局所形式・全立体角 `FactorImaginary(kV, s2, m)`（本番, L3012）

入力: `s2` = $s_g^2$ [nm⁻²]、`m` = $B$ [nm²]（NaN・0 なら 0 を返す）。$\gamma=1+eV/(m_0c^2)$、$k_0$ = `EnergyToElectronWaveNumber`（**$k_0=1/\lambda$、相対論補正込み、nm⁻¹**）。

$$
f''(\mathbf g)=\gamma\,\frac{k_0}{2}\int_{4\pi}\! d\Omega\;f_e(s_1)\,f_e(s_2)\Big[1-e^{-B\,(s_1^2+s_2^2-s_g^2)}\Big]\times 0.01
$$

$$
s_{1,2}=\tfrac12\big|\mathbf K\mp\tfrac{\mathbf g}{2}\big|,\qquad \mathbf K=k_0(\sin\theta\cos\phi,\ \sin\theta\sin\phi,\ \cos\theta-1)
$$

- $f_e$ は `ElasticFactorWithTail`（Å）、×0.01 で Å² → nm²。戻り値の単位は **nm**（$f_e$ と同じ）。※ doc コメントは「無次元」と書いているが、`getU` は実部と同じく $\gamma/(\pi V)$ を掛けているので実質 nm（→ 実装上の疑問点 Q2）。
- **$\mathbf g$ はビームに垂直（x 軸）に置く**近似（$|\mathbf g|$ しか渡らない）。旧閉形式も $|g|$ だけの関数だった。
- 光学定理との対応（Claude 追記・要確認）: $k=1/\lambda$ 規約で $\mathrm{Im}f(0)=k\sigma/2$。相対論的 $f=\gamma f_e$ とすると $\gamma^2$ が出るが、実装は $f''$ に $\gamma$ を 1 つ、`getU` でさらに $\gamma$ を 1 つ掛けて合計 $\gamma^2$ となり整合する。
- 数値積分:
  - θ を **$s=k_0\sin(\theta/2)$ の等比的区切り** `AbsorptionSEdges = [0, 0.3, 1, 3, 9, 27, 81, 243, 729, ∞]` nm⁻¹ で分割（$\theta=2\arcsin(s/k_0)$、$s\ge k_0$ の区間は θ=π で打ち切り。末尾 ∞ で数 MV でも全立体角を覆う）。
  - 各区間 Gauss–Legendre **12 点**（`AbsorptionThetaPoints`）、φ は [0, π] の **8 点**則を φ=π/2 で折り畳んだ 4 節点・重み ×4（`PhiFolded`, 総和 2π を Assert）。φ→−φ, φ→π−φ の偶対称を利用。
  - 精度根拠（260811Cl）: 同じ被積分関数の 48×48 参照に対し $|g|\le50$ nm⁻¹ で $f''(0)$ 比 0.05 % 以内（12×8 で 0.023 %、16×16 で 0.0024 %）。
  - DW 項 $1-\exp[m(s_g^2-2q_0)]$ は φ 不変なので θ ごとに 1 回評価（$q_0=(|\mathbf K|^2+s_g^2)/4$）。
- 検証基準: `FactorImaginaryGaussianAnalytic`（旧閉形式、Gauss 型専用）と `useTail:false` で一致するはず。

閉形式（`FactorImaginaryGaussianAnalytic`, L3129, 本番未使用）:

$$
f''=\frac{2\pi\gamma}{k_0}\sum_{i,j}a_ia_j\left[\frac{e^{-s_g^2\,b_ib_j/(b_i+b_j)}}{b_i+b_j}-\frac{e^{-s_g^2\,(b_ib_j-B^2)/(b_i+b_j+2B)}}{b_i+b_j+2B}\right]
$$

（$s_g^2$, $B$ は Å 系へ換算: `s2*=0.01`, `m*=100`。$b_i+b_j=0$ の項は 0。）

### 6.2 局所形式・円環検出器 `FactorImaginaryAnnular(kV, g, m, inner, outer, nTheta=60, nPhi=20)`（L3221）

同じ被積分関数を θ ∈ [inner, outer]（**ラジアン**）、φ ∈ [0, 2π] で Gauss–Legendre（既定 60×20）。$\mathbf g$ は 3 成分すべて使う（$g_z\ne0$ 可、φ 折り畳みは不適用）。STEM の HAADF 等の非弾性ポテンシャル用。

### 6.3 非局所形式 `FactorImaginaryAnnular(kV, g, h, m, inner, outer, …)`（L3274）

$$
f''(\mathbf g,\mathbf h)=\gamma\frac{k_0}{2}\int d\Omega\ f_e\!\big(\tfrac12|\mathbf K-\mathbf g|\big)f_e\!\big(\tfrac12|\mathbf K-\mathbf h|\big)\Big[1-e^{-B(s_{Kg}^2+s_{Kh}^2-s_{g-h}^2)}\Big]\times0.01
$$

`getU(kV, g, h, inner, outer)` から使われる。

### 6.4 本番未使用のもの
- `FactorImaginaryAnnularFlatEwald`（「未完成」）、`FactorImaginaryAnnular2`（「遅いので没?」）。どちらもテール未導入の Gauss 直書き。→ 実装上の疑問点 Q6。

### 6.5 `getU` での使い方（概要。詳細は [bloch-wave.md](bloch-wave.md)）
`BetheMethod.getU`（BetheMethod.cs L3863）:

$$
U_g=\frac{\gamma}{\pi V_c}\sum_{\text{sites}}\mathrm{Occ}\,f_e(s_g)\,e^{-B s_g^2}\,e^{+2\pi i\mathbf g\cdot\mathbf r},\qquad
U'_g=\frac{\gamma}{\pi V_c}\sum\mathrm{Occ}\,f''(\mathbf g)\,e^{-Bs_g^2}\,e^{+2\pi i\mathbf g\cdot\mathbf r}
$$

- $V_c$ [nm³]、$U$ [nm⁻²]。係数は「Kirkland の教科書 p120 参照」とコメント。
- 位相は **$e^{+2\pi i\mathbf g\cdot\mathbf r}$**（「20240524 符号をプラスに変更 三菱・中村」）。運動学 F（$e^{-2\pi i}$）と逆（[structure-factor.md §4](structure-factor.md)）。
- 非等方のとき $B$ は方向依存の実効値 $m=\mathbf h'^{\mathsf T}\beta\mathbf h'/s_g^2$（$\mathbf h'=R^{\mathsf T}\mathbf h$）を**等方 B として** $f''$ にも渡す。$g=0$ は `BisoEffective`。
- イオンの扱い `ElasticIonModel`（既定 `Neutral`）: Neutral = `[Z][0]`、IonScreened = 割当イオンの Peng 有限部分、IonFull = 有限部分 + $g\ne0$ に Mott–Bethe 単極子 $0.2393366\,\Delta Z/s^2$ [nm]（$s^2$ [nm⁻²]）。単位胞が非中性（$|\sum q|>10^{-6}$）なら単極子を落とす（260614Cl）。単極子は実部にだけ入り、吸収積分には入らない。
- キャッシュ `uDictionary` は hkl 合成キーのみでエネルギーを含まない（`MeanAbsorptionLengthNm` はエネルギーごとに Clear）。

---

## 7. NIST 弾性散乱（参考リンク）

`Atom/NistElastic.cs` と [`Atom/NistElastic/README.md`](../../../Atom/NistElastic/README.md) は EBSD モンテカルロ（`EBSD/MonteCarlo.cs`）の弾性散乱角サンプリング専用。**回折（Bloch 波・構造因子）の経路からは参照されていない**ことを grep で確認した。

---

## 実装上の疑問点

- **Q1** `AtomStatic.cs` L1209–1213 `MottBetheMonopoleCoefficient` の doc コメントが「$s^2$ [Å⁻²]、f は [nm]」と書くが、値 0.2393366 は **$s^2$ [nm⁻²]・f [nm]** に対応する（$0.023934\,\text{Å}\cdot\text{Å}^{-2}$ → $f[\text{nm}]=0.1\times0.023934\,\Delta Z/(0.01\,s^2_{nm})=0.23934\,\Delta Z/s^2_{nm}$）。`getU` は nm⁻² を渡しており計算は正しい。同ファイル L2936 の注記は正しく nm⁻² と書いている。コメントのみの不整合。
- **Q2** `FactorImaginary*` の `<returns>`「戻り値が無次元量」と `getU` 内コメント「答えは無次元」（BetheMethod.cs L3918）は、実際の次元（$\gamma k_0/2\cdot\int f_ef_e\,d\Omega\times0.01$ = nm、実部 $f_e$ と同次元）と食い違う。`getU` は実部と同じ $\gamma/(\pi V)$ を掛けており、nm として扱うのが整合的。コメントのみの問題と思われる。
- **Q3** 運動学 F（`Crystal.GetStructureFactor` → `Atoms.GetAtomicScatteringFactorForElectron`）は `SubNumberElectron` の Peng（割当イオンの有限部分 = IonScreened 相当）を使い、`BetheMethod.ElasticIonModel`（既定 Neutral）を参照しない。イオン指定した結晶では、同じ反射の運動学 $|F|$ と動力学 $U_g$ が別のモデルになる。意図か要確認。
- **Q4** 非等方温度因子のときの TDS 吸収（BetheMethod.cs L3900–3918）: $f''$ の被積分関数の DW 指数 $B(s_1^2+s_2^2-s_g^2)$ に、$\mathbf g$ 方向で定義した実効等方値 $m=\mathbf h'^{\mathsf T}\beta\mathbf h'/s_g^2$ を使っている。$\mathbf K\pm\mathbf g/2$ は $\mathbf g$ と別方向を向くので、厳密には方向ごとに $\beta$ を評価すべき量。近似として許容か要確認（等方なら厳密）。
- **Q5** 中性子の同位体加重（Atoms.cs L366–371）は `Dictionary<int,double>` の列挙順と `NeutronCoherentScattering[Z][i+1]` の並びが一致することに依存する。.NET の Dictionary は削除がなければ挿入順に列挙されるのが実装上の挙動だが、仕様上の保証はない（Claude 追記・要確認）。
- **Q6** `FactorImaginaryAnnularFlatEwald`（L3174）は動径 $R$ の積分に Jacobian $R$ が掛かっていない（被積分関数に `R` の乗算がなく、コメントは「sinTheta を外に出して」のまま）うえ、前因子 $\gamma/(\pi^2k_0)$ が他の版の $\gamma k_0/2\cdot(1/k_0^2)=\gamma/(2k_0)$ と合わない。本番未使用（「未完成」）なので実害なし。復活させる際は要修正。
- **Q7** 弾性散乱因子に関するコメントの食い違い: 260811Cl のブロック（L2883–2884）は「表の出典は s ≤ 6 Å⁻¹ の当てはめと称している」、260923Cl の表注記（L1222–1223）は「s ≤ 2 Å⁻¹ 用（ITC 4.3.2.2）」。後者が新しい訂正と読めるが、前者のコメントが未更新。

## 要確認事項

1. TDS 吸収の局所形式で $\mathbf g$ をビームに垂直に置く近似（`FactorImaginary`）の誤差評価が、HOLZ 反射・EBSD の大角 $\mathbf g$（$g_z$ が大きい）でどの程度か。
2. 局所形式の対称化（$\mathbf K\mp\mathbf g/2$）を円環検出器（`FactorImaginaryAnnular` 局所版）に適用すると、検出器の積分域が「入射ビーム基準」ではなく「$\mathbf g/2$ だけずれた点基準」になる。STEM の HAADF 用としてこの定義で意図どおりか。
3. X 線白色（`WaveColor.FlatWhite`）の強度が $d^2$、中性子白色が $d^4$（Crystal.cs L1759, L1770）。「平坦」がエネルギー当たりか波長当たりかの違いと解釈できる（Claude 追記・要確認）が、意図の確認を要する。→ [structure-factor.md](structure-factor.md) 要確認事項にも記載。
4. Temari `Fx` をX 線構造因子に使う予定の有無（現状 WK のみ）。
5. `XrayScattering`（Cromer–Mann）・`ElectronScatteringEightGaussian` がライブラリ内で未使用。上位アプリ（ReciPro）での利用の有無。
6. 8-Gaussian 表の出典が私信のみ。
7. `IsotopeAbundance` の収録版（コード内に「※収録版の特定は要確認」）。
8. 異常分散 f′/f″ をイオンでも中性元素の値で代用していること（化学シフト・価数依存は無視）。
9. Peng 表を $s>2$ Å⁻¹ で外挿している弾性実部（`es.Factor`）を、大角の HOLZ や EBSD 計算でも使っていること（コメントでは実用範囲外として据え置き）。

要確認事項: 9 件
