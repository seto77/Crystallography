# EBSD 計算手法 — 目次 / EBSD methods — index

本稿はコードから Claude が起こした下書き（2026-11-05）。著者未査読。
This is a draft written by Claude from the code (2026-11-05). Not yet reviewed by the author.

対象コミット / Target commit: `9ffcf47` (main). 対象 / Scope: `EBSD/`（`MonteCarlo.cs` を含む / including `MonteCarlo.cs`）。

> この目次は EBSD 担当の分だけを載せる。回折・動力学（Bethe 法、Kikuchi の動力学計算）の文書は別ブランチで作成中のため、`docs/methods/README.md` への統合は両方のマージ後に行う。
> This index covers the EBSD part only. The diffraction/dynamical documents (Bethe method, dynamical Kikuchi calculation) are being written on a separate branch; merging into `docs/methods/README.md` is left until both are merged.

## 文書 / Documents

| # | 日本語 | English | 内容 / Contents |
|---|---|---|---|
| 1 | [ja/ebsd-geometry.md](ja/ebsd-geometry.md) | [en/ebsd-geometry.md](en/ebsd-geometry.md) | lab・試料・結晶の座標系、回転行列の向き、検出器幾何、パターン中心、バンド線 ↔ 面法線 / frames, rotation convention, detector geometry, pattern centre, band line ↔ plane normal |
| 2 | [ja/ebsd-master-pattern.md](ja/ebsd-master-pattern.md) | [en/ebsd-master-pattern.md](en/ebsd-master-pattern.md) | マスターパターンの構築手順・格納形式・Rosca–Lambert（正方/六方）・格子上の対称操作・`.rpmp` ファイル・検出器への投影 / build, storage, Rosca–Lambert grids, symmetry, file format, projection |
| 3 | [ja/ebsd-monte-carlo.md](ja/ebsd-monte-carlo.md) | [en/ebsd-monte-carlo.md](en/ebsd-monte-carlo.md) | MC の断面積・阻止能・IMFP・損失分布・軌跡・コヒーレンス破壊、MC 分布（射出半球ビン、エネルギー分布、垂直深さ → 経路長） / MC transport and the (energy × depth) distribution |
| 4 | [ja/ebsd-pattern-composition.md](ja/ebsd-pattern-composition.md) | [en/ebsd-pattern-composition.md](en/ebsd-pattern-composition.md) | 合成モデル 0/1/2、画素ごとの重み、$A(E)$、非晶質層、照合用の 1 枚への縮約 / composition models, per-pixel weights, $A(E)$, amorphous layer |
| 5 | [ja/ebsd-indexing.md](ja/ebsd-indexing.md) | [en/ebsd-indexing.md](en/ebsd-indexing.md) | バンド検出、pair-angle、Radon テンプレート、辞書、晶帯軸、統合探索、幾何較正 / band detection, pair-angle, Radon template, dictionary, zone-axis, combined search, geometry calibration |

読む順 / Suggested order: 1 → 2 → 3 → 4 → 5.

## 実装上の疑問点（一覧） / Implementation questions (summary)

コードは変更していない。各項目の根拠は本文。 / No code was changed; see each document for the evidence.

| ID | 要旨 / Summary | 文書 / Doc |
|---|---|---|
| G-I1 | `FromPatternCenter` は $\mathbf n\cdot\mathbf C=+DD$ の配置を返す。$\tau>0$ で検出器が射出側にある配置（$\mathbf n\cdot\mathbf C<0$）では PC 経由の往復で DetY の符号が反転する / returns the $\mathbf n\cdot\mathbf C=+DD$ branch; round trip flips DetY when $\mathbf n\cdot\mathbf C<0$ | 1 |
| G-I2 | doc の「signed L = −DD」の $L$ が未定義 / "signed L" undefined in the doc | 1 |
| M-I1 | 正方格子の対角鏡映 2 種の名前（コメント）が入れ替わっている（常に対で使うので結果に影響なし） / diagonal-mirror labels swapped (no functional effect) | 2 |
| C-I1 | CSDA 経路で脱出区間の連続損失を引かない / CSDA omits the loss on the escaping flight | 3 |
| C-I2 | `Depth` は「脱出直前の最後の事象の深さ」で、doc の「源の深さ」との対応は呼び出し側次第 / `Depth` is the last-event depth | 3 |
| C-I3 | コメントと実装の小さなずれ（80 %/95 %、`tailOnsetEv` の定義） / small comment–code mismatches | 3 |
| P-I1 | model 0/1 で $\bar A$ の重み（$\cos^3\alpha$ あり）と画素値（なし）が揃っていない / weighting mismatch in $\bar A$ for models 0/1 | 4 |
| P-I2 | 単一スライスの model 2 は $\Delta t$ で割らない / single-slice model 2 is not divided by $\Delta t$ | 4 |
| I-I1 | `SaturateCap` の doc（プラトー 0.4–1.0）と値（8）が合わない / doc vs value mismatch | 5 |
| I-I2 | pair-angle の重複除去が対称を考慮しない / duplicate removal ignores symmetry | 5 |
| I-I3 | 幾何較正の開始点が G-I1 の符号規約に依存する / calibration start depends on G-I1 | 5 |

## 要確認事項 / Items to confirm

合計 31 件（G: 6、M: 6、C: 7、P: 5、I: 7）。代表例:
31 items in total (G 6, M 6, C 7, P 5, I 7). Examples:

- G-1: Bunge オイラー角との対応が `EBSD/` に無い（回転行列のみ）。結晶直交系は $\mathbf c\parallel z$、$x\parallel\mathbf a^*$。 / No Bunge-angle mapping in `EBSD/`; crystal frame is $\mathbf c\parallel z$, $x\parallel\mathbf a^*$.
- G-3: マスターパターンは視線 $-\hat{\mathbf k}_{\rm out}$ で引く。Bethe 側の方向の意味（相反定理）との対応は回折担当の文書で確認。 / Master pattern is sampled along $-\hat{\mathbf k}_{\rm out}$; check against the Bethe-side convention.
- G-4: 他社の PC 規約（Bruker, EDAX/TSL, Oxford, EMsoft）との変換は未実装。 / No conversion to vendor PC conventions.
- M-1: 平面値 $|A|^2$ と「累積 $M(t)$」という解釈が Bethe 側の出力定義と一致するか。 / Meaning of $|A|^2$ and the cumulative-$M(t)$ interpretation.
- C-1: MC 分布へ渡す源の深さ（最後の事象 / 最後の非弾性 / 最後のコヒーレンス破壊）の選択は `FormEBSD`（このリポジトリに無い）にある。 / The choice of source depth lives in `FormEBSD`, outside this repository.
- P-3: $A(E)$ の $E_c$ 既定値（0.8 keV）は Si 20 kV の実測 1 枚での経験値。 / Default $E_c$ is an empirical value from one Si 20 kV pattern.
- I-7: Radon の等価判定が点群より大きい群で縮約しうる（極性結晶）。 / Radon equivalence may over-reduce for polar crystals.

## 書き方の約束 / Conventions used in these documents

- すべての記述は実装箇所（ファイル・クラス・メソッド）で確認したもの。確信が持てない点は各文書末尾の「要確認事項」。 / Every statement was checked against the code; uncertain points are listed under "Items to confirm".
- 文献はコード中のコメントに出てくるものだけを名前で挙げた（Joy & Luo 1989、Jablonski 2008、TPP-2M、NIST SRD 64、Bote–Salvat、Vos & Winkelmann 2019、Rosca、Shirley/Rosca、EMsoft `Lambert.f90`）。それ以外は「（Claude 追記・要確認）」。 / Only references named in code comments are cited; anything else is marked "(added by Claude; to be confirmed)".
- コメント中の `tools/…`・`.project-guidance/…`・`FormEBSD` はこのリポジトリに無い。 / `tools/…`, `.project-guidance/…` and `FormEBSD` referenced in comments are not in this repository.
