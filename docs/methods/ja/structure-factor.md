本稿はコードから Claude が起こした下書き（2026-11-05）。著者未査読。

# 構造因子と運動学的強度

対象: `Crystal/Crystal.cs`（`GetStructureFactor` 3 版、`SetVectorOfG`、`SetVectorOfPlane`、`SetVectorOfG_KikuchiLine`、`SetPeakIntensity`、`GetDspacingList`）。
原子散乱因子・温度因子の定義と単位は [scattering-factors.md](scattering-factors.md)、動力学計算のポテンシャル $U_g$ は [bloch-wave.md](bloch-wave.md) を参照。

---

## 1. 記号と単位

| 量 | 実装上の定義 |
|---|---|
| $\mathbf g=h\mathbf a^*+k\mathbf b^*+l\mathbf c^*$ | `A_Star` 等（nm⁻¹）。$|\mathbf g|=1/d$ |
| $s^2$ | `vec.Length2/4.0` または `1/d²/4.0` = $(\sin\theta/\lambda)^2$ [**nm⁻²**] |
| 原子座標 | 分率座標 `atom.X, atom.Y, atom.Z`（`Atoms.Atom[]` は等価位置展開済み）|
| 散乱因子の単位 | 電子: nm、X 線: electrons×0.1、中性子: fm（いずれも Occ 込み）|
| $F$ の単位 | 上記と同じ（体積で割らない）|

## 2. 構造因子 `Crystal.GetStructureFactor`

実装は 3 版:

| シグネチャ | 用途 |
|---|---|
| `public static Complex GetStructureFactor(WaveSource, Atoms[], (double h,k,l), double s2, double xrayEnergyKeV=NaN, bool anomalousDispersion=true)` (L1545) | 非整数指数対応の最適化なし版 |
| `private static Complex GetStructureFactor(WaveSource, Atoms[], (int h,k,l), double s2, Complex[] neutronFactors=null, double xrayEnergyKeV=NaN, (double fp,double fpp)[] xrayDispFactors=null, bool anomalousDispersion=true)` (L1620) | 本体（260605Cl 実数振幅高速パス）|
| `public Complex GetStructureFactor(WaveSource, (int h,k,l), double xrayEnergyKeV=NaN, bool anomalousDispersion=true)` (L1712) | インスタンス版、$s^2=|\mathbf g|^2/4$ を内部で計算 |

### 2.1 式（実装どおり）

$$
F(\mathbf h)=\sum_{\text{Atoms}}\ \sum_{j\in\text{等価位置}} f_{\text{Atoms}}(s)\;T_j(\mathbf h)\;\exp\!\big(-2\pi i\,\mathbf h\cdot\mathbf r_j\big)
$$

- **位相規約は $e^{-2\pi i\mathbf h\cdot\mathbf r}$**: `AddStructureFactorContribution` (L1511) は `realSum += A.Re cosφ + A.Im sinφ`, `imagSum += A.Im cosφ − A.Re sinφ`（$\varphi=2\pi\mathbf h\cdot\mathbf r$）で、$A e^{-i\varphi}$ を足している。
- $f$ には Occ が含まれる（[scattering-factors.md §2](scattering-factors.md)）。連続する同元素・同 sub・**同 Occ** のサイトでは $f$ を再利用（260715Ch で Occ 条件を追加。以前は Occ の違うサイトで先頭の値を使い回していた）。
- 温度因子:
  - 等方（`Dsf.UseIso`）: $T=\exp(-B_{iso}s^2)$、$B$ [nm²]。NaN なら 1。
  - 非等方: $T_j=\exp[-(\mathbf h'^{\mathsf T}\beta\,\mathbf h')]$、$\mathbf h'=$`atom.Operation.ConvertPlaneIndex(h,k,l)`$=R^{\mathsf T}\mathbf h$（β は無次元、非対角は 2 倍で足す）。NaN なら 1。
  - 非等方でも $\mathbf h=0$ の特別扱いはない（$T=1$ になるので問題なし）。
- 中性子: $f$ は $s$ 非依存。`SetVectorOfG` は全反射共通で事前計算した `neutronFactors` を渡す。
- 電子: 実数振幅パス（`AddStructureFactorContributionReal`）。**γ は掛けない**、**`ElasticIonModel` は参照しない**（割当イオンの Peng 有限部分）。

### 2.2 X 線異常分散（260606Cl）

$$
f=f_0+0.1\,\mathrm{Occ}\,\big(f'-i f''\big),\qquad f_0=\texttt{WK.Factor}(s^2)\cdot\mathrm{Occ}
$$

（`XrayDispersionFactor`, L1529。0.1 は electrons→内部単位。f′/f″ が NaN なら該当項 0。）

- $f''$ に負号が付くのは $e^{-2\pi i\mathbf h\cdot\mathbf r}$ 規約のため（コメント「International Tables の + 規約と一致させる、検証済」）。すなわち $F_{\text{code}}(\mathbf h)=\overline{F_{\text{ITC}}(\mathbf h)}$（ITC: $f=f_0+f'+if''$, $e^{+2\pi i}$）で、$|F|^2$ と Bijvoet 差の符号関係は一致する（Claude 追記・要確認: 共役関係はコードからの導出）。
- 適用条件 `IsXrayDispersionActive`: `anomalousDispersion && wave==Xray && E>0 有限 && Xraylib.Enabled`。既定 `anomalousDispersion = true` だが `xrayEnergyKeV` 既定 NaN なので、**エネルギーを渡さない呼び出しでは分散なし**。
- `SetVectorOfG`, `SetPeakIntensity` はサイトごとに f′/f″ を 1 回だけ事前計算（`xrayDispFactors`）。
- 異常分散があると $F(-\mathbf h)\ne\overline{F(\mathbf h)}$ なので、`SetVectorOfG` は ± 両側を独立計算する（260606Cl で共役埋めを撤廃）。

## 3. 逆格子点の列挙と運動学的強度

### 3.1 `SetVectorOfG(dMin, dMax, wavesource, maxNum=25000, xrayEnergyKeV, anomalousDispersion)` (L1186)
- 原点から殻状に拡張して $|\mathbf g|<1/d_{min}$ の点を最大 `maxNum` 個集め、$|\mathbf g|$ 昇順。$1/d_{max}$ 未満を除く。
- 各点の ± を `VectorOfG` に格納。消滅則 `Symmetry.GetFirstExtinctionRule(h,k,l)` が非 null の反射は $F=0$（計算しない）。
- `RawIntensity = |F|²`、`RelativeIntensity = RawIntensity / max`（全 0 なら正規化しない）。**ローレンツ・偏光・体積の補正なし**。

### 3.2 `SetVectorOfPlane`, `SetVectorOfG_KikuchiLine`
いずれも $|F|^2$ と最大値規格化のみ。菊池線用は消滅則で $F=0$、`RelativeIntensity` 昇順にソート。動力学的な菊池バンドプロファイルは [kikuchi-dynamical.md](kikuchi-dynamical.md)。

### 3.3 粉末強度 `SetPeakIntensity(waveSource, waveColor, λ, …)` (L1720)
$\sin\theta=\lambda/(2d)$, $2\theta=2\arcsin(\sin\theta)$（**ラジアン**）。`Plane[i].strHKL` を `+` で分割した各 hkl について $F$ を計算し、多重度 `Multi[j]` を掛けて和をとる。$V_c^2$ = `CellVolumeSquare`。

| 線源 | 色 | `eachIntensity` |
|---|---|---|
| X 線 | Monochrome | $|F|^2 m\,\dfrac{1+\cos^2 2\theta}{\sin 2\theta\,\sin\theta}\Big/V_c^2$ |
| X 線 | FlatWhite | $|F|^2 m\,d^2/V_c^2$ |
| 電子 | — | $|F|^2 m\,\dfrac{1}{\sin2\theta\,\sin\theta}\Big/V_c^2$ |
| 中性子 | Monochrome | $|F|^2 m\,\dfrac{1}{\sin2\theta\,\sin\theta}\Big/V_c^2$ |
| 中性子 | 白色 | $|F|^2 m\,d^4/V_c^2$ |

最後に全ピークの最大値で規格化（`Intensity`）。
- X 線単色の LP 因子は $\dfrac{1+\cos^22\theta}{\sin2\theta\sin\theta}=\dfrac{1+\cos^22\theta}{2\sin^2\theta\cos\theta}$ で、教科書形（Claude 追記・要確認: $\frac{1+\cos^22\theta}{\sin^2\theta\cos\theta}$）の 1/2。規格化するので無影響。モノクロメータ偏光・吸収・消衰の補正はない。
- 粉末で Friedel 対が同じ 2θ に重なることは、代表 hkl の $|F|^2$ を使う近似（Bijvoet 非対称の多重度平均はしない、260606Cl コメント）。

### 3.4 `GetDspacingList(λ, maxNum=512, bestNum=16, maxTwoThetaDeg=150)` (L1343)
X 線・LP 因子込み強度 $|F|^2 m\,(1+\cos^22\theta)/(\sin2\theta\sin\theta)$ の上位 `bestNum` 個の d を返す。
- 261105Cl: $2\theta$ 上限（**度**、既定 150°）を追加（2θ→180° で LP が発散するため）。
- 261105Cl: 同値面は Laue 群（`IsRootPlane(..., addInversion: true)`）で数え、多重度も Laue 群で。相対差 $10^{-6}$ 未満の d は 1 本にまとめて強度を加算。
- 異常分散・エネルギーは使わない（`GetStructureFactor(WaveSource.Xray, …)` を分散なしで呼ぶ）。

## 4. 動力学計算との規約の違い

| 箇所 | 位相 | 散乱因子 | γ | 規格化 |
|---|---|---|---|---|
| `Crystal.GetStructureFactor` | $e^{-2\pi i\mathbf h\cdot\mathbf r}$ | 割当イオンの Peng（電子）/ WK(+f′,f″)（X 線）| なし | なし（$F$ そのもの）|
| `BetheMethod.getU` | $e^{+2\pi i\mathbf g\cdot\mathbf r}$（20240524 変更）| `ElasticIonModel` に従う Peng | $\gamma$ | $\gamma/(\pi V_c)$ → nm⁻² |
| `BetheMethod.ComputeIonizationU` | $e^{-2\pi i\mathbf G\cdot\mathbf r}$（260801Cl「負位相を直接生成」）| 内殻イオン化形状因子 | なし | $k_{vac}/(2\pi V_c)$ |

`getU` の正位相がどこで共役されて全体の規約と整合しているかは [bloch-wave.md](bloch-wave.md) 側の記述を参照（`ComputeIonizationU` のコメントは「既存 getU の『正位相を作って共役』という経験的規約」と書いている）。ALCHEMI/EDX 側は [alchemi-ionization.md](alchemi-ionization.md)。

---

## 実装上の疑問点

- **S1** 電子線の運動学 $F$ と動力学 $U_g$ でイオンの扱いが異なる（[scattering-factors.md](scattering-factors.md) Q3）。`Crystal.GetStructureFactor` は `SubNumberElectron` をそのまま使い、`BetheMethod.ElasticIonModel`（既定 Neutral）に従わない。イオンを割り当てた結晶では、運動学スポット強度と動力学計算の $U_g$ が別の模型になる。
- **S2** 位相規約が `Crystal.GetStructureFactor`（$-$）と `BetheMethod.getU`（$+$）で逆。`getU` のコメントは「対称心の結晶の計算が上手くいくはず」と経験的な変更理由を書いており、非中心対称結晶・吸収あり（$U'_g$ が複素）の場合に全体として正しい規約になっているかは B 側（[bloch-wave.md](bloch-wave.md)）での確認事項。ここではコードの事実のみ記す。

## 要確認事項

1. LP 因子の 1/2 の差（規格化で消えるが、絶対強度を出す用途が将来あるなら要注意）。
2. X 線白色 $d^2$・中性子白色 $d^4$ の意図（「平坦」がエネルギー当たり／波長当たりの違いか）（Claude 追記・要確認）。X 線白色に偏光因子がない（固定 2θ なら定数）こと。
3. 電子の粉末強度に用いた $1/(\sin2\theta\sin\theta)$ のローレンツ因子の出典（コード内に記載なし）。
4. `GetDspacingList` の強度は異常分散を使わない。波長が吸収端近傍の場合の扱い。
5. `SetVectorOfG` などの $F$ は温度因子込みだが、Occ・DW が未設定（NaN）のとき $T=1$ に落ちる仕様で良いか。
6. 異常分散時の $F_{\text{code}}=\overline{F_{\text{ITC}}}$ の関係（Bijvoet 対の向き）を、ReciPro 側の表示（どちらを $+\mathbf h$ と呼ぶか）と照合すること。

要確認事項: 6 件
