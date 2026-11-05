本稿はコードから Claude が起こした下書き（2026-11-05）。著者未査読。

# EBSD のモンテカルロ: 電子輸送と (エネルギー × 深さ) 分布

対象コミット: `9ffcf47`（main）。対象ファイル: `EBSD/MonteCarlo.cs`、`EBSD/EbsdBackscatterSimulator.cs`、`EBSD/EbsdMonteCarloDistribution.cs`。座標系は [ebsd-geometry.md](ebsd-geometry.md) §1。

全体の役割:

1. `MonteCarlo` が電子 1 本を試料内で追跡し、脱出した電子の深さ・方向・エネルギーなどを返す。
2. `EbsdBackscatterSimulator.Run` がそれを大量に並列に回し、エネルギーしきい値で選別する。
3. `EbsdMonteCarloDistribution` が電子を射出半球のビンに分け、ビンごとに「射出エネルギー分布 $G(E)$」「源の垂直深さの平均 $\lambda_d(E)$」「電子の割合 $F$」を当てはめる。合成器（[ebsd-pattern-composition.md](ebsd-pattern-composition.md)）はこれを画素ごとに内挿して (エネルギー × 深さ) の重みを作る。

---

## 1. 物質パラメータ（`MonteCarlo` ctor）

```
MonteCarlo(z, a, ρ, kev, tilt, thresholdKev = 2,
           stoppingPowerModel = JablonskiModified2008,
           elasticScatteringModel = MottNistSampler2023,
           inelasticScatteringModel = DiscreteBulkDiimfpApproximation,
           valenceElectronCount = null, bandGapEv = null, atoms = null, meanIonizationPotentialEv = null)
```

| 記号 | 引数 | 単位 | 規約 |
|---|---|---|---|
| $Z$ | `z` | – | **原子数平均** $\sum n_iZ_i/\sum n_i$ |
| $A$ | `a` | g/mol | **原子数平均**（原子 1 個あたりの質量） |
| $\rho$ | `_ρ` | g/cm³ | |
| $E_0$ | `kev` | keV | 入射エネルギー |
| $\tau$ | `tilt` | rad | 試料傾斜（$x$ 軸まわり） |
| $E_{\rm th}$ | `thresholdKev` | keV | 追跡打ち切り（既定 2） |
| $N_v$ | `valenceElectronCount` | – | TPP-2M の価電子数。null なら `atoms` から原子数平均、それも無ければ $Z$ から推定 |
| $E_g$ | `bandGapEv` | eV | 既定 0 |
| $J$ | `meanIonizationPotentialEv` | eV | null なら `atoms` から **Bragg 則**（電子数重みで $\ln J_i$ を平均）、無ければ $J(Z)$ |

doc の注意（260921Cl）: $Z$・$A$・$N_v$ は内部で $\rho Z/A$（Bethe 阻止能）と $N_v\rho/A$（TPP-2M）という**比**として使うので、同じ重み（原子数）で揃えること。以前は $Z$ と $N_v$ が質量加重で、水素を含む化合物ほど比が狂っていた（Cu₂(OH)₃Cl で $Z/A$ が 1.93 倍、$N_v\rho/A$ が 1.63 倍）。`GetMeanAtomicParameters(atoms)` は Multiplicity × Occ で重み付けし、$Z\notin[1,99]$ の原子を除く。

元素の $J$（Joy & Luo 1989。`ElementIonizationPotentialEv`）:

$$ J(Z) = \begin{cases} 11.5\,Z & (Z\le12)\\ 9.76\,Z + 58.5\,Z^{-0.19} & (Z>12)\end{cases}\quad[\mathrm{eV}]. $$

その他の派生量:

- Debye–Waller の組成平均 $\bar B$ [nm²] = $\sum(\text{Occ}\cdot\text{mult})\,B_{\rm iso,eff}/\sum(\text{Occ}\cdot\text{mult})$、下限 `DiffuseScatteringFactor.BisoFloorNm2`（コメントでは 1e-3 nm²）。
- TPP-2M の材料定数: $u = N_v\rho/A$、$E_p = 28.8\sqrt u$ [eV]、$\beta = -0.10 + 0.944/\sqrt{E_p^2+E_g^2} + 0.069\rho^{0.1}$、$\gamma = 0.191\rho^{-1/2}$、$C = 1.97 - 0.91u$、$D = 53.4 - 20.8u$。
- プラズモンカットオフ $q_c = \omega_p/v_F$: $k_F = (3\pi^2 n_v)^{1/3}$、$n_v = N_v\rho N_A/A$ [nm⁻³]、$q_c = E_p/((\hbar^2/m_e)k_F)$（$\hbar^2/m_e = 0.0761996$ eV·nm²）を $2\pi$ で割って $k=1/\lambda$ 系で保持（`InelasticLocalizationQNm`。Si で ≈ 1.9 nm⁻¹ とコメント）。
- 局在とみなす損失 `InelasticLocalizedLossEv` $=\max\bigl(1.8\max(E_p,\,E_{\min}+0.5),\,30\bigr)$ eV、$E_{\min}=\max(E_g\text{ or }1,\,0.5)$。

---

## 2. 輸送パラメータ

1 eV 刻み（0 eV〜$\lceil 1000E_0\rceil$ eV、$E_0$ は keV。0 eV の要素には 1 eV の値を入れる）で事前計算し（`BuildTransportParameterCache`）、ホットループでは `Math.Round(E·1000)` の添字で引く（`GetTransportParametersRef`。範囲外は末尾にクランプ）。

### 2.1 弾性散乱

**Mott / NIST SRD 64**（`MottNistSampler2023`、既定）: 元素ごとのコード生成済み PCHIP テーブル（`AtomStatic`）。エネルギー点は 111 点、50 eV〜36.4 keV の対数等間隔（50 eV〜20 keV は NIST の sampler 由来、20 keV 超は DCS から継ぎ足し。260603Cl）。

- 混合物: 元素 $i$ の数密度 $n_i$ [nm⁻³] と断面積 $\sigma_i$ から巨視的断面積 $\Sigma = \sum_i n_i\sigma_i$、$\lambda_{\rm el} = 1/\Sigma$、有効微視断面積 $\sigma = \Sigma/\sum n_i$。対数エネルギー軸で隣接 2 点を線形補間。
- 散乱角: 散乱元素を $n_i\sigma_i/\Sigma$ の CDF で選び、その元素の PCHIP 累積分布を逆関数法で評価（エネルギー点は最近傍）。
- 範囲外のエネルギー、または生成データが無い元素は遮蔽 Rutherford へフォールバック。

**遮蔽 Rutherford**（`ScreenedRutherford`）:

$$ \alpha = \frac{3.4\times10^{-3}\,Z^{2/3}}{E\,[\mathrm{keV}]},\qquad
\sigma = \pi\Bigl(\frac{Ze^2}{4\pi\varepsilon_0\,m v^2}\Bigr)^2\frac{1}{\alpha(1+\alpha)},\qquad
\lambda_{\rm el} = \frac{A}{N_A\rho\,\sigma}, $$

$m v^2$ は相対論的（`EnergyToElectronMass` × `EnergyToElectronVelositySquared`）。$mv^2 = E(E+2m_0c^2)/(E+m_0c^2)$ を代入すると、コード中にコメントで残っている Joy の式 $5.21\times10^{-21}(Z/E)^2\,\tfrac{4\pi}{\alpha(1+\alpha)}\bigl(\tfrac{E+511}{E+1022}\bigr)^2$ [cm²] と一致する（Claude 確認）。散乱角は

$$ \cos\theta = 1 - \frac{2\alpha R}{1+\alpha-R},\qquad R\sim U(0,1). $$

### 2.2 阻止能

**Joy & Luo (1989)**（`JoyLuo1989`）:

$$ \frac{dE}{ds} = \frac{c_3}{mv^2}\ln\Bigl(1.166\,k + \frac{0.583\,mv^2}{e\,J}\Bigr),\qquad k = 0.0299\ln Z + 0.7307, $$

$c_3 = -Z N_A\rho e^4/(4\pi\varepsilon_0^2 A)$（単位換算込み、keV/nm）。

**Jablonski (2008) 修正式**（`JablonskiModified2008`、既定）:

$$ S\,[\mathrm{eV/Å}] = \frac{D_1\,E^{D_2}\,\ln(D_3E)\,(\rho+D_4)^{D_5}}{\lambda_{\rm in}(E)\,[\mathrm{Å}]},\quad (D_1..D_5) = (7.89271,\ 0.0117088,\ 0.0545126,\ -0.0254488,\ 0.326907), $$

$E$ は eV。$dE/ds = -0.01\,S$ [keV/nm]。$\lambda_{\rm in}$ が無効なら Joy–Luo に戻す。

### 2.3 非弾性平均自由行程（TPP-2M）

$$ \lambda_{\rm in}(E)\,[\mathrm{Å}] = \frac{E}{E_p^2\bigl[\beta\ln(\gamma E) - C/E + D/E^2\bigr]},\qquad \gamma E\le1 \text{ または分母} \le0 \text{ なら NaN}. $$

1 事象の平均損失 $\langle\Delta E\rangle = |dE/ds|\,\lambda_{\rm in}$（「平均的には $\langle\Delta E\rangle/\lambda_{\rm in}$ = 阻止能」）。全散乱率 $\Sigma_{\rm tot} = 1/\lambda_{\rm el} + 1/\lambda_{\rm in}$、弾性の確率 $p_{\rm el} = (1/\lambda_{\rm el})/\Sigma_{\rm tot}$。

### 2.4 1 事象の損失分布

| モデル | 損失 $\Delta E$ |
|---|---|
| `ContinuousSlowingDownApproximation` | 離散事象なし。各飛行で $E \mathrel{+}= s\,dE/ds$ |
| `DiscreteMeanLoss` | 常に $\langle\Delta E\rangle$ |
| `DiscreteBulkDiimfpApproximation`（既定） | 3 成分の混合 PDF（低損失 $\propto(\omega+0.35E_p)^{-1.6}$、プラズモン = $E_p$ のガウス〔$\sigma = \max(1.5, 0.18E_p)$ eV〕、高損失テール $\propto\omega^{-2.2}$〔開始 $\min(\omega_{\max},\max(1.8E_p,30))$ eV〕）を平均が $\langle\Delta E\rangle$ になるよう重み付け。10 eV 刻みのエネルギー点で CDF を作り、ガイドテーブル付き逆関数法 |
| `DiscreteDrudeValenceInnerShell`（260922Cl） | 事象率と平均損失は上と同じ。価電子は拡張 Drude 損失関数 $\mathrm{Im}[-1/\varepsilon] = \omega_p^2\gamma\omega/((\omega^2-\omega_q^2)^2+\gamma^2\omega^2)$、$\omega_q = \omega_p + q^2/2$（原子単位、$\gamma = 4$ eV）を $q$ 積分した DIIMFP、内殻は Bote–Salvat の電離断面積比で副殻を選び $[B, E/2]$ で $\propto\omega^{-2}$。内殻を選ぶ割合は平均損失が $S\lambda_{\rm in}$ になるように決める（doc: Si 20 keV で 11.7 %、Vos & Winkelmann 2019 表 1 の 12.0 % と一致）。16 エネルギー点（対数等間隔）のサンプラー |

いずれも $\Delta E\le E$ にクランプ。

---

## 3. 軌跡（`GetBackscatteredElectronDetail`）

初期状態: 原点、$\mathbf v=(0,0,-1)$、$E=E_0$、$d=0$。

### 3.1 離散非弾性モデル（既定の経路）

```
while E > E_th:
    s = −ln(max(R, ε)) / Σ_tot(E)               # 次事象までの距離
    d' = d + s (v_y sinτ − v_z cosτ)
    if d' < 0: break                            # 表面を抜けた → 脱出
    d = d'
    if R < p_el:  弾性: cosθ をサンプル、φ = 2πR、方向を回転
                  熱散漫判定 (§4.1) が真なら「最後のコヒーレンス破壊」= d
    else:         非弾性: ΔE をサンプル、「最後の非弾性」= (d, E, E−ΔE, v)
                  局在判定 (§4.2) が真なら「最後のコヒーレンス破壊」= d
                  E ← max(E − ΔE, 0)
```

方向の回転は進行方向 $\mathbf v$ を極軸とする局所座標で $(\sin\theta\cos\phi,\sin\theta\sin\phi,\cos\theta)$ を回す閉形式（$v_z+1<10^{-7}$ の特異点では直接代入）。10 回に 1 回 $\mathbf v$ を再正規化。

### 3.2 CSDA 経路

各ループで $s = -\lambda_{\rm el}\ln R$、2 回目以降の飛行の前に弾性散乱で方向を変え、$d$ を更新してから $E \mathrel{+}= s\,dE/ds$。離散非弾性が無いので「最後の非弾性」は定義されない（`HasLastInelasticEvent = false`、NaN）。

### 3.3 返り値（`BackscatteredElectronDetail`）

| 項目 | 意味 |
|---|---|
| `Depth` | ループを抜けた時点の $d$ = **脱出直前の最後の事象（弾性・非弾性を問わない）の垂直深さ** [nm]。エネルギー打ち切りで抜けた電子ではその時点の深さ |
| `Direction` | 最後の飛行方向（lab 系、単位ベクトル） |
| `Energy`, `TotalEnergyLoss` | 脱出時のエネルギー、$E_0 - E$ [keV] |
| `HasLastInelasticEvent`, `LastInelasticDepth`, `…EnergyBeforeLoss`, `…AfterLoss`, `LastInelasticDirection` | 最後の非弾性事象 |
| `HasLastDecoherenceEvent`, `LastDecoherenceDepth` | 最後のコヒーレンス破壊事象（§4） |
| `CoherentInelasticCount`, `CoherentAngularVarianceRad2`, `CoherentPathLengthNm` | 調査用（`CollectCoherenceDiagnostics` = false では前 2 つは 0） |

脱出の判定は `Energy > energyThreshold`（`EbsdBackscatterSimulator.Run` の引数）で行う。`Run` は `Parallel.For` で `loop` 本を回し、脱出電子を `EbsdBackscatteredElectron`（上記 + `Position` = Schmidt 投影 `Stereonet.ConvertVectorToSchmidt(sampleRotation · Direction)`）にまとめる。進捗は 1024 電子ごとに共有カウンタへ加算し、約 1 % ごとに通知。

> **どの深さを「源」とするか**（`Depth`、`LastInelasticDepth`、`LastDecoherenceDepth` のどれを MC 分布に渡すか）は呼び出し側（`FormEBSD`、このリポジトリに無い）が決める。`EbsdMonteCarloDistribution` の doc には「源を『最後の非弾性』で取る既定モード」と「最後のコヒーレンス破壊事象で取るモード (dec)」の 2 つが出てくる（要確認事項 C-1）。

---

## 4. コヒーレンス破壊の判定（260919Cl）

### 4.1 熱散漫（弾性）— `IsThermalDiffuseElasticEvent`

$$ x = \bar B\,k^2(1-\cos\theta) = 2\bar B s^2,\quad s = k\sin(\theta/2),\ k = 1/\lambda\,[\mathrm{nm}^{-1}];\qquad p_{\rm TDS} = 1-e^{-x}. $$

$x>20$ なら確率 1、$x<0.03$ なら $x(1-x/2)$。干渉性の成分 $e^{-x}$ は Bloch 波側（マスターパターン）が扱うので、MC では「何も起きなかった」と同じ扱い（方向は変える）。

### 4.2 非弾性 — `InelasticDecoherenceProbability`

- $\Delta E \ge$ `InelasticLocalizedLossEv` なら 1（内殻・単一粒子励起）。
- それ未満は、運動量移行 $q$ がカットオフ $q_c$ を超える確率。角度分布を Bethe ridge $\theta_c^2=\Delta E/E$ で打ち切った Lorentzian $dP\propto\theta\,d\theta/(\theta^2+\theta_E^2)$、$\theta_E = \Delta E/2E$、$q^2 = k^2(\theta^2+\theta_E^2)$:

$$ \theta_q^2 = q_c^2/k^2 - \theta_E^2,\qquad
p = \begin{cases} 1 & \theta_q^2\le0\\ 0 & \theta_q^2\ge\theta_c^2\\ 1 - \dfrac{\ln(1+\theta_q^2/\theta_E^2)}{\ln(1+\theta_c^2/\theta_E^2)} & \text{otherwise.}\end{cases} $$

（$q<q_c$ の集団励起は Bloch 状態を保つ、という扱い。doc はエネルギーフィルター実験でプラズモン損失電子にも菊池バンドが残る事実に対応させている。）

---

## 5. MC 分布（`EbsdMonteCarloDistribution`）

### 5.1 入力

```
EbsdMonteCarloDistribution((Depth, Vec, Energy)[] bseList, beamEnergy, sampleTilt,
    energies, depths, binCount = 18, amorphousLayerNm = 0,
    energyWeightDeadKeV = NaN, energyFilterMinKeV = NaN,
    absorptionLengthNm = null, energyHistogram = false)
```

- `Depth`: 源の**垂直深さ** [nm]（どの深さかは §3.3 の注）。`Vec`: 射出方向（lab）。`Energy`: 射出エネルギー [keV]。
- `sampleTilt`: **MC を回したときの**傾斜。`energies`: マスターパターンのエネルギー格子（降順）。`depths`: マスターパターンの深さ格子（経路長）。

### 5.2 電子の重み

$$ w(E) = \phi(E)\,\Theta(E - E_{\min}),\qquad \phi(E) = \begin{cases}\max(0,\,E - E_{\rm dead}) & E_{\rm dead} \text{ 有効}\\ 1\end{cases} $$

$E_{\rm dead}\ge E_0$ や $E_{\min}\ge E_0$ のときは（全電子が 0 になるので）無効に落とす。$\phi$ は「蛍光体応答（試行）」、$\Theta$ は理想的なエネルギーフィルター（260922Cl）。

### 5.3 射出半球のビニング（260921Cl）

試料系の射出方向 $\mathbf s = R_x(\tau)\,\mathbf{Vec}$、$s_z\le0$ は捨てる（ただし重みの母数 $\sum w$ には入れる）。**Lambert 等積ディスク**:

$$ (X,Y) = \sqrt{\frac{2}{1+\hat s_z}}\,(\hat s_x,\hat s_y),\qquad X^2+Y^2\le2, $$

正方形 $[-\sqrt2,\sqrt2]^2$ を `binCount` × `binCount` に切る（`DirectionToBinCoords`、ビン中心が整数座標）:

$$ b_x = (X+\sqrt2)\frac{n_b}{2\sqrt2} - \tfrac12,\qquad b_y = (\sqrt2 - Y)\frac{n_b}{2\sqrt2} - \tfrac12 . $$

電子は `Round` で最寄りのビン $(b_i,b_j)$ に入れる。円板内側のビンの立体角は $8/n_b^2$ sr（既定 18 で 0.0247 sr）。分布を作る側と合成側は同じ `DirectionToBinCoords` を使う（doc で強調）。Rosca–Lambert 正方写像をやめた理由（doc）: 対角線 $|x|=|y|$ で $C^0$ なので、背景平坦化（≈ ∇²）が X 字の暗線を作った。

ビンの被覆率 `BinCoverage`（円板に入る面積の割合、64×64 点で標本化）。被覆率 ≥ `MinBinCoverage` = 0.5 のビンは測った値を使い、$F$ は被覆率で割って「ビン全面あたり」にする。それ未満のビンは `ExtendFieldOutward`（8 近傍の重み付き平均〔辺 1、角 1/2〕を 1 層ずつ）で内側から延長する。

### 5.4 ビンごとの当てはめ

- 電子の割合: $F_b = \sum_{b} w / \sum_{\rm all} w$（非晶質層内の源も含む）。
- 非晶質割合: $f_{A,b} = \sum_{d<a} w / \sum_b w$（$a$ = `amorphousLayerNm`）。層内の源は「変調なし」に回し、結晶内の深さは $d-a$（層の底から）で数える。
- **エネルギー分布 $G_b(E)$**（総和 1 に正規化）:
  - 既定（`energyHistogram = false`）: 左右非対称ガウス。重み付き平均 $\bar E$、片側ずつの重み付き分散 $\sigma_L,\sigma_R$（$n\le1$ または $\sigma<0.01$ なら 0.5 keV）、$G(E_i)\propto\exp(-(E_i-\bar E)^2/2\sigma_{L/R}^2)$。フィルターがあれば各スライスの区間 $[E_i-h^-, E_i+h^+]$（隣との中点）のうち $E_{\min}$ を超える割合を掛ける。
  - `energyHistogram = true`（260925Cl）: 生のヒストグラム。電子の重みを隣り合う 2 枚の格子点へ線形に配る（三角核）。格子の外は端の 1 枚へ。doc: 非対称ガウスは $E_0$ 側の鋭い山を 2〜3 分の 1 に過小評価する。
- **深さ $\lambda_{d,b}(E)$**: エネルギー窓 $[E_i-h, E_i+h)$（$h$ = 平均刻みの半分）ごとに結晶内の源の垂直深さを平均し（4 本以上の窓だけ、下限 1 nm）、$\lambda(E) = l_a + l_bE + l_cE^2$ を最小二乗（有効窓 1 点なら定数、2 点なら直線）。格子上で $[1,10^6]$ nm にクランプ。指数分布の平均 = $\lambda$ という解釈。
- 退避: 電子 10 本未満のビンは全ビン合算の $G$ と $\lambda$（状態 `Global`）、結晶内の源が少ない・有効窓が無いビンは合算の $\lambda$（`GlobalLambda`）、合算でも無ければ $\lambda = 10^6$ nm（深さ一様、`NoDepthData`）。

### 5.5 垂直深さ → 経路長（`FillPathLengthWeights`、深さ写像 A2）

マスターパターンの深さ格子は**出射方向に沿った経路長** $t$ で、MC の源深さは垂直深さ $d$。$t = d/\mu$（$\mu=\cos\chi$ = 射出方向の試料系 $z$ 成分）。$d$ を平均 $\lambda_d$ の指数分布とすると $t$ も指数分布で、減衰率

$$ \alpha = \mu/\lambda_d(E)\quad[\mathrm{nm}^{-1}]. $$

`absorptionLengthNm` $=\lambda_{\rm abs}(E)$ を渡したとき（dec モード用、260923Cl）は $\alpha \leftarrow \max(\alpha - 1/\lambda_{\rm abs},\ 0.1\alpha)$（MC の熱散漫の生存とマスターパターンの平均吸収の二重計上を除く）。

- `sliceMass = true`（model 2）: 区間質量 $m_j = e^{-\alpha t_{j-1}}(1-e^{-\alpha\Delta t_j})$。$\alpha\Delta t\le10^{-4}$ では $1-e^{-x}$ を 4 次の級数で評価。
- `sliceMass = false`（model 0/1）: $m_j = e^{-\alpha(t_j-t_0)}\Delta t_j$（密度 × 区間幅の右端則。$t_0$ で割ってアンダーフローを防ぐ）。
- どちらも**エネルギーごとに条件付き正規化**してから重み $G(E)$（または $F\,G(E)$）を掛ける: $w_{e,j} = G_e\,m_j/\sum_j m_j$。格子 $t_{n_D}$ を超える尾部は同じエネルギーの中で配り直され、エネルギー分布は変わらない。合計が非有限なら質量を最浅スライスへ。

互換用の配列 `BinWeights`（model 0/1、ビン内総和 1）と `BinAbsoluteSliceWeights`（model 2、ビン内総和 $F_b$）はビン中心の $\mu$（`LambertBinCenterMu` $=\max(10^{-6}, 1-(X^2+Y^2)/2)$）で作った近似。表示合成は画素の $\mu$ で作り直す。

### 5.6 その他

- `GlobalDepthWeights` / `GlobalDepthSliceWeights`（非晶質層があるときだけ）: $\sum_b F_b\,(\text{ビンの重み})$。非晶質層内の源の基準強度に使う。
- `ComputeElectronWiseCoherence(E0, E_c)`（260927Cl、`energyHistogram = true` のときだけ）: 電子ごとに $A(E)$ を掛けたヒストグラム $G_A$ を $G$ と同じ規格化定数で割って保持。$0\le G_A\le G$ を全節点で検査し、違反は例外。
- `ComputeRangesFromMC`: 損失の 95 % 点と深さの 99 % 点（QuickSelect）。
- `ComputePathLengthUpperBound`: 経路長 $t=(d-a)/\mu$ の分位点（既定 0.999）。検出器を渡せばその画像に当たる電子だけ、渡さなければ $\mu\ge0.2$。該当なしなら垂直深さの 99 % 点 / 0.5。
- `ComputeGridFromRanges(E0, ΔE95, T, β = 2)`: エネルギーは 16 段、刻み $\max(0.1, \mathrm{round}(\Delta E_{95}/15, 0.1))$ keV、$E_0$ から降順、下限 1 keV。深さは 40 点、$\beta>0$ なら `BuildGeometricDepthGrid`:

$$ t_i = T\,\frac{e^{\beta i/N}-1}{e^\beta-1}\quad(i=1..N),\qquad T \text{ は有効数字 2 桁へ切り上げ、各点 0.001 nm に丸めて狭義単調}. $$

  $\beta\le0$ は従来の等間隔。

---

## 6. 実装上の疑問点

### C-I1 CSDA 経路で、脱出する飛行区間の連続損失を引かない

`GetBackscatteredElectronDetail`（CSDA）の順序は「$s$ を引く → 散乱 → $d' = d + s(\cdots)$ → $d'<0$ なら break → $d=d'$ → $E \mathrel{+}= s\,dE/ds$」。表面を抜ける最後の飛行では、試料内にある部分の連続損失も引かれない。1 区間（$\sim\lambda_{\rm el}$）ぶんの損失なので影響は小さいはずだが、厳密には脱出エネルギーを過大評価する（離散モデルでは区間中にエネルギーは変わらないので問題ない）。

### C-I2 「最後の事象の深さ」を表す `Depth` の意味が doc と食い違う可能性

`GetBackscatteredElectrons` の doc は「返り値は、深さ (nm単位)…」とだけあり、MC 分布側の doc は `Depth` を「源の垂直深さ」と呼ぶ。実装上の `Depth` は §3.3 の通り「脱出直前の最後の事象の深さ」で、最大到達深さでも最後の非弾性の深さでもない。呼び出し側がどれを渡しているかの確認が要る（要確認事項 C-1）。

### C-I3 コメントと実装の小さなずれ

- `ComputeRangesFromMC` 内のコメント（260602Cl）は「80%/99% パーセンタイル」だが、実装は 95 % / 99 %（260919Cl に変更済み）。
- `InelasticLocalizedLossEv` のコメントは「`CreateBulkLossSamplerEntry` の `tailOnsetEv` と同じ定義」だが、前者は $\max(1.8\max(E_p, E_{\min}+0.5), 30)$、後者は $\min(\omega_{\max}, \max(1.8E_p, 30))$ で、$E_p < E_{\min}+0.5$ や $\omega_{\max}$ が小さいときに一致しない。

---

## 7. 要確認事項

- **C-1** MC 分布へ渡す源の深さ（`Depth` / `LastInelasticDepth` / `LastDecoherenceDepth`）の選択と既定値。`FormEBSD` の実装で確認すること。
- **C-2** 遮蔽パラメータの指数は $Z^{2/3}$（コード）で、Joy の式（コード中のコメント）は $Z^{0.67}$ と書かれることが多い（Claude 追記・要確認）。差は数値的に無視できるが、どちらを正とするか。
- **C-3** `IsThermalDiffuseElasticEvent` は遮蔽 Rutherford の経路でも適用される。Mott 断面積・遮蔽 Rutherford 断面積のどちらを使う場合も、熱散漫の割合を $1-e^{-2Bs^2}$ で切り出すという扱いでよいか（断面積がすでに静的原子の値であることを前提にしている）。
- **C-4** 非弾性の局在判定に使う Lorentzian と Bethe ridge 打ち切りは「1 事象の損失が決まった後の角度分布」の近似で、MC 自体は非弾性事象で方向を変えない。この近似の適用範囲（低エネルギー・重元素）。
- **C-5** エネルギー窓ごとの平均深さ → $\lambda(E)$ の 2 次多項式という当てはめは、深さ分布が指数分布であることを前提にしている。浅い側に山を持つ分布（入射ビームの侵入で生じる）との差の評価。
- **C-6** `MinBinCoverage = 0.5` と `ExtendFieldOutward` の重み（辺 1・角 1/2）は経験的な選択か。
- **C-7** 電子の重み $\phi(E)=\max(0,E-E_{\rm dead})$（蛍光体応答）は doc で「試行」とされている。既定値と適用範囲。
