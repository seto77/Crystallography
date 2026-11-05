本稿はコードから Claude が起こした下書き（2026-11-05）。著者未査読。

# EBSD パターンの合成: マスターパターン × MC 重み

対象コミット: `9ffcf47`（main）。対象ファイル: `EBSD/EbsdPatternComposer.cs`、`EBSD/EbsdMonteCarloDistribution.cs`（`ComposeGlobalWeightedPattern`）。前提: [ebsd-geometry.md](ebsd-geometry.md)（画素 → 方向）、[ebsd-master-pattern.md](ebsd-master-pattern.md)（平面 $M_{e,d}$ の意味）、[ebsd-monte-carlo.md](ebsd-monte-carlo.md)（MC 分布）。

`EbsdPatternComposer` は `FormEBSD` から移した表示用の合成器で、GUI に依存しない。`EbsdPatternProjector`（照合用、単一スライス）との違いは doc の通り: こちらは視野全体のラスター・MC 重み付きの全スライス合成・方位は結晶の現在値。

---

## 1. 記号

- $M_{e,d}(\hat{\mathbf d})$: エネルギー $E_e$、深さ（経路長）$t_d$ の平面を、画素の視線 $\hat{\mathbf d}$（結晶系）で補間した値。累積応答 $M(t_d)=\int_0^{t_d}R\,dz$ と解釈される（`DepthIntervals` の doc）。
- $\Delta M_{e,d} = \max\bigl(0,\ M_{e,d}-M_{e,d-1}\bigr)$（$M_{e,-1}=0$）。$\Delta t_d$ は `ComputeDepthIntervals`。
- $W_{e,d}(i)$: 画素 $i$ の (エネルギー × 深さ) 重み（§3）。
- $\bar M_{e,d}$: 平面の**全球**方向平均（$+Z$・$-Z$ 半球の平均の平均。§4.3）。
- $A_e$: 損失依存のコントラスト係数（§4.1）。$f_A(i)$: 非晶質層に源を持つ電子の割合（§4.2）。

---

## 2. 単一スライスのモデル

参照テーブル（`BuildLookupTable`、エネルギー・深さに依らない）を作ってから、指定スライスを書き込む。

| モデル | メソッド | 画素値 |
|---|---|---|
| 0 current | `ApplySingleSliceModel0` | $M_{e,d}$（= model 1 で係数 1） |
| 1 globally normalized | `ApplySingleSliceModel1` | $c_k\,M_{e,d}$、$c_k = 1/\sum_{\rm cells}(M^{+}_k + M^{-}_k)$（両半球の全セルの和。0 以下なら 0） |
| 2 differential | `ApplySingleSliceModel2` | $\max(0,\,M_{e,d}-M_{e,d-1})$（**$\Delta t$ では割らない**） |

正方格子の Rosca–Lambert 格子は等積なので、$c_k$ のセル和は全球積分に比例する。どちらの平面も `GridSize²` に満たなければ全画素 0。

---

## 3. 画素ごとの重み $W_{e,d}(i)$（`EvaluatePathLengthWeights`、深さ写像 A2）

1. 画素の**射出方向**（試料系）を `ExitRay` で作る（[ebsd-geometry.md](ebsd-geometry.md) §3.6）。$o_z\le0$（試料表面より下）の画素は 0。
2. $\mu = o_z/|\mathbf o|$、ビン座標 $(b_x,b_y)$ = `DirectionToBinCoords(o)`（MC 分布と同じ写像）。
3. **3 次 B スプライン**の 4 × 4 タップ（`BinSplineTaps`、$t = b-\lfloor b\rfloor$）:

$$ q(t) = \tfrac16\bigl[(1-t)^3,\ 3t^3-6t^2+4,\ -3t^3+3t^2+3t+1,\ t^3\bigr],\quad \text{タップ } \lfloor b\rfloor-1,\dots,\lfloor b\rfloor+2 \text{（添字だけ端でクランプ）}. $$

   $\sum q = 1$、$\sum q\cdot(\text{タップ位置}) = b$（1 次を厳密に再現）、$C^2$。近似であって補間ではない（ビン中心の値は $(v_{i-1}+4v_i+v_{i+1})/6$）。$b$ は $[-3, n_b+2]$ に丸めるだけで、クランプしない。
4. タップ $t$ の重み $\omega_t = c_t\,s_b$、$c_t = q_x q_y$。model 2 は $s_b = F_b$（電子の割合）、model 0/1 は $s_b = [F_b>0]$。$\omega_t\le0$ のタップは飛ばす。

$$ N(E_e) = \sum_t \omega_t G_{b_t}(E_e),\qquad \lambda(E_e) = \frac{\sum_t\omega_t G_{b_t}(E_e)\,\lambda_{b_t}(E_e)}{N(E_e)}\ \ (N=0 \text{ なら } 10^6\ \mathrm{nm}). $$

5. `FillPathLengthWeights(N, λ, μ, depths, Δt, sliceMass)`（[ebsd-monte-carlo.md](ebsd-monte-carlo.md) §5.5）で $W_{e,d}$。model 2 は `sliceMass = true`（区間質量、総和 $\sum_e N_e$）、model 0/1 は `false`（密度 × 区間幅）。
6. 非晶質割合 $f_A = \sum_t\omega_t f_{A,b_t}/\sum_t\omega_t$。

doc の注: 旧版は「ビン中心の $\mu$ で作った重み配列」を内挿していたが、16×16 の射出半球ビン 1 個の中でも $\chi$ は約 11° 動くので、**パラメータ（$G$、$\lambda$、$F$）を内挿してから画素の $\mu$ で換算**する。$\mu=1$ でも旧版とは一致しない（非線形）。

---

## 4. MC 重み付きのモデル（`ApplyWeightedModel0/1/2`）

前提検査: `dist.MatchesGridOf(mp)`（エネルギー・深さ格子の値が相対 $10^{-9}$ で一致）。一致しなければ例外（260727Cl。格子がずれると別スライスの重みを黙って引くため）。

### 4.1 A(E) 無効・非晶質層なしの基本形

$$
\begin{aligned}
\text{model 0:}\quad V_i &= \sum_{e,d} W_{e,d}(i)\,M_{e,d}(\hat{\mathbf d}_i),\\
\text{model 1:}\quad V_i &= \sum_{e,d} W_{e,d}(i)\,c_k\,M_{e,d}(\hat{\mathbf d}_i),\\
\text{model 2:}\quad V_i &= \Omega_i\sum_{e,d} W_{e,d}(i)\,\frac{\Delta M_{e,d}(\hat{\mathbf d}_i)}{\Delta t_d},\qquad \Omega_i = \cos^3\alpha_i .
\end{aligned}
$$

- model 2 の $\Delta M/\Delta t$ は区間平均応答 $\bar R$（MC 側の重みは区間質量なので割らない。260726Cl）。
- $\Omega_i$ は画素の立体角 $d\Omega/dA = D/|\mathbf r|^3 = \cos^3\alpha/D^2$ の $\cos^3\alpha$ 部分（$\cos\alpha = D/|\mathbf r|$、$D$ = `ExitRay.PlaneDistance`）。射出半球のビンは等立体角なので「立体角あたり」の重みに画素の立体角を掛ける。定数 $1/D^2$ は落とす（表示は自動伸張、ZNCC は定数倍に不変）。**model 0/1 には掛けない**。

### 4.2 損失依存のコントラスト係数 $A(E)$（260920Cl〜）

$$ A_e = \exp\!\Bigl(-\frac{\max(0, E_0-E_e)}{E_c}\Bigr),\qquad E_0 = \texttt{BeamEnergyKeV}\ (\text{NaN なら } \max_e E_e),\quad E_c = \texttt{CoherenceLossDecayKeV}. $$

$E_c$ が NaN・0 以下なら無効（全スライスが満額のコントラスト）。全部が $1-10^{-12}$ 以上なら `null`（無効と同じ）。式は `CoherenceFactors` の 1 箇所だけにある（表示合成と ZNCC の目的関数が同じ重みを使う前提のため）。

有効時の合成（`CohAccum`、260921Cl）:

$$
V^{\rm coh}_i = \underbrace{\sum W\bar M}_{D_c}
\;+\; \bar A\,\Bigl(\sum W\Bigr)\,\frac{\sum W A\,I - \sum W A\bar M}{\sum W A},
$$

ここで $I$ は model に応じた $M$ / $c_kM$ / $\Delta M/\Delta t$、$W$ は $W_{e,d}(i)$（model 1 は $W c_k$）、和は $(e,d)$ について。すなわち $V = \sum W\bar M + \bar A\,(\sum W)\,(\langle I\rangle_A - \langle\bar M\rangle_A)$、$\langle\cdot\rangle_A$ は $WA$ で重み付けした正規化平均。第 1 項は $A$ 無効時の明るさ（変調なし成分）で、第 2 項はコントラストを全面で $\bar A$ 倍にする。

$\bar A$（`MeanCoherentFraction`）は画素に依らない 1 個のスカラーで、**検出器面**を 48 × 48 点で標本化し、各点の重みに立体角 $\cos^3\alpha$ を掛けて

$$ \bar A = \frac{\sum_{\rm pts}\sum_{e,d}\cos^3\alpha\;W_{e,d}\,A_e}{\sum_{\rm pts}\sum_{e,d}\cos^3\alpha\;W_{e,d}} $$

（model 1 は $W c_k$）。行ごとに足してから行の順に足すので、並列の順序に依らない。範囲は `EbsdRasterView.HalfWidth/HalfHeight`（検出器の物理寸法）で、表示のパン・ズームに依らない。

`energyHistogram = true` の分布では、節点の $A_e$ の代わりに画素ごとの実効 $A$（$G_A/G$ を同じ 4 × 4 タップで内挿、`ComputeElectronWiseCoherence`。260927Cl）を使う。

doc の経緯（要約）: 旧版は $A$ で失った分を 1 個のスカラーの台座として足していたが、$\sum W A$ が画素ごとに 8×8 の痕跡を持つため、背景平坦化の後にビン間隔のブロードな暗線が出た。上の式で明るさの痕跡を消している。コントラストの面内変化（取り出し角による低損失電子の割合の違い）は意図的に捨てている。

$E_c$ の既定値についての doc の記述: Si 20 kV の実測 1 枚で校正した経験値（0.7 → 0.8 keV）で、物質・電圧に依らない定数ではない。2 物質目（Botallackite）の結果は誤差の中で物質依存を決められない、とある（数値の詳細は `CoherenceLossDecayKeV` の doc）。

### 4.3 非晶質層（260919Cl）

層内に源を持つ電子は菊池変調を持たない一様成分として寄与する:

$$ V_i = \Omega_i\Bigl[(1-f_A)\,V^{\rm coh}_i + f_A\sum_{e,d}W_{e,d}\,\bar M^{\rm ref}_e\Bigr] \quad(\text{model 0/1 は }\Omega_i=1\text{、model 1 は } W \text{ に } c_k \text{ を掛ける}). $$

- $\bar M^{\rm ref}_e$（`sphereMeans`）: 平面ごとの方向平均（値が 0 のセルは除く。model 2 は $\Delta M/\Delta t$ の平均）を、全ビン合算の深さ重み（`GlobalDepthWeights` / `GlobalDepthSliceWeights`）で深さ方向に平均した「エネルギーごとの基準値」を、$+Z$・$-Z$ で平均したもの（`CollapseToEnergyReference`）。半球ごとだと赤道で段差が出るので全球平均にする。
- これに対し §4.2 の $\bar M_{e,d}$（`planeMeans`）は、エネルギーへ畳まない平面ごとの全球平均（$A$ の台座はコヒーレント項と $(e,d)$ で 1 対 1 に打ち消す必要があるため）。
- 平面平均は `(MasterPattern, Distribution, differential)` が同じならキャッシュする。

### 4.4 背景（model 2 の任意出力）

`ApplyWeightedModel2(..., background)` を渡すと $B_i = \Omega_i[(1-f_A)\sum W\bar M + f_A\sum W\bar M^{\rm ref}]$（菊池の変調を除いたもの）も書く。

### 4.5 診断用の出力

`ComputeEffectiveCoherenceModel2`: model 2 の合成が画素 $i$・節点 $e$ で掛ける $A$ の係数と $\bar A$ を返す（外部の参照実装 `tools/EbsdProfileFit` で同じ離散化を再現するため。ツールはこのリポジトリに無い）。

---

## 5. 照合用の 1 枚への縮約（`EbsdMonteCarloDistribution.ComposeGlobalWeightedPattern`）

ZNCC 照合・幾何較正では画素ごとに重みを作らず、重みを 1 本に平均してから平面を 1 枚に縮約し、`EbsdPatternProjector` で投影する:

$$ \tilde w_{e,d} = \frac{A_e\,\bar W_{e,d}}{\sum A\bar W},\qquad
\tilde M^{\pm}(\text{cell}) = \sum_{e,d}\frac{\tilde w_{e,d}}{\Delta t_d}\,\max\bigl(0, M^{\pm}_{e,d}-M^{\pm}_{e,d-1}\bigr)\quad(\tilde w<10^{-15} \text{ は飛ばす}). $$

- $\bar W$ は引数 `sliceWeights`（推奨: `EbsdPatternComposer.ComputeDetectorAverageSliceWeights` = 検出器画像の画素平均、`stride` 4 で間引き、各画素に $\cos^3\alpha$、地平線より下の画素は 0 として分母に数える）。null なら全ビンの `BinAbsoluteSliceWeights` の和（射出半球全体の平均）。
- $(1-A)$ の台座は足さない（ZNCC は定数加算・定数倍に不変）。
- doc の注: 「重みを平均してから 1 枚に縮約して投影」する近似で、$E_c$ の最終較正は画素ごとの合成で行うこと。

---

## 6. 実装上の疑問点

### P-I1 $\bar A$ の重みと画素値の重みが model 0/1 で揃っていない

`MeanCoherentFraction` は model 0/1 でも各標本点の重みに $\cos^3\alpha$ を掛けるが、`ApplyWeightedModel0/1` の画素値には $\Omega_i$ を掛けない。doc は「W は各モデルの画素の重み × 画素の立体角」と定義しているので意図的とも読めるが、model 2 では「表示合成の重みの画素平均」と一致するのに対し、model 0/1 では一致しない。

### P-I2 単一スライスの model 2 は $\Delta t$ で割らない

`ApplySingleSliceModel2` は $\max(0, M_d - M_{d-1})$ をそのまま書く。重み付きの model 2 と `ComposeGlobalWeightedPattern` は $/\Delta t_d$ する（260726Cl）。単一スライス表示はスライス間の比較をしないので見かけ上は問題にならないが、不等間隔の深さ格子でスライスを切り替えて見比べると、明るさが区間幅に比例して変わる。

---

## 7. 要確認事項

- **P-1** model 0（累積 $M$ を密度 × 区間幅の重みで足す）と model 1 の物理的な位置付け。doc は名前（current / globally normalized master）だけで、どの近似に対応するかの説明が無い。
- **P-2** model 2 の $\Omega_i=\cos^3\alpha$ は「検出器が画素あたりに受ける電子数」に対応するが、蛍光体の発光の角度依存・カメラの感度は含まない。実測との比較で背景平坦化を前提にしているか。
- **P-3** $A(E)$ の既定 $E_c$（UI の欄で 0.8 keV とコメント）と、`energyHistogram = true` のときの再較正の要否（`EbsdMonteCarloDistribution` の doc に「true では E_c の較正をやり直す必要がある」）。
- **P-4** `MeanCoherentFraction` は 48 × 48 点の固定標本化。検出器のアスペクト比が大きいときの十分性。
- **P-5** `BinSplineTaps` は近似（ビン値を 1 ビンぶん平滑化）なので、射出分布の鋭い角度依存（例: チャネリングによる BSE 収率の異方性）は表現できない。MC 分布は結晶性を含まない（非晶質 MC）ので問題にならない、という理解でよいか。
