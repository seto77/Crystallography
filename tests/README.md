# Crystallography の回帰テスト

`Crystallography.Tests` は、空間群まわりのデータとコードを独立した基準と照合する xUnit テストです (261105Cl 新規)。

## 実行方法

```bash
# Windows
dotnet test tests/Crystallography.Tests -c Release

# Linux / macOS (net10.0-windows をビルドするため EnableWindowsTargeting が必要)
dotnet test tests/Crystallography.Tests -c Release -p:EnableWindowsTargeting=true
```

Linux でも動きます。WinForms に触れる経路はテストしていません。テストプロジェクトはソリューション (`Crystallography.sln`) には入れていないので、通常のビルドには影響しません。本体の `Crystallography.csproj` は `tests\**` をコンパイル対象から外しています。

## 基準データ

`Crystallography.Tests/Reference/symmetry_reference.json` には次のデータが入っています。

- **spglib**：Hall 番号 1〜530 の対称操作（センタリング込み）、Hall 記号、H-M 記号（短縮・完全）、Schoenflies 記号、点群
- **cctbx**：各設定の Wyckoff 位置（文字、多重度、サイト対称の型、Wyckoff 位置への射影演算子）。生成時に、cctbx の対称操作が spglib と全設定で一致することを確認しています。
- ReciPro 独自の通し番号 531〜540（センタリング付き三斜晶 C-1, I-1, A-1, B-1, F-1, A1, B1, C1, F1, I1）は、cctbx で Hall 記号から生成しています。

ReciPro の通し番号 1〜530 は、spglib の Hall 番号と同じ順番です。再生成は次のとおりです（cctbx は spglib より先に import する必要があり、スクリプトはそうなっています）。

```bash
pip install spglib cctbx-base
python tests/Crystallography.Tests/Reference/generate_symmetry_reference.py
```

## 照合している内容

| テスト | 内容 |
|---|---|
| `GeneralPositions_MatchReference` | 一般位置の座標式（センタリング込み）が、基準の対称操作の集合と一致する（540 設定） |
| `SeitzOperations_MatchGeneralPositionsOneToOne` | Seitz 形式の対称操作 (`PositionOperations`) が、一般位置と添字ごとに一致する。原子に持たせる対称操作（異方性 DW に使う）がこれに依存する |
| `ExtinctionRules_MatchOperations` | 手書きの消滅則 (`HasExtinction`) が、対称操作から導いた系統的消滅と一致する（\|h\|,\|k\|,\|l\| ≤ 10） |
| `EquivalentPlanes_/EquivalentAxes_MatchPointGroupOrbit` | 等価な面・軸が、点群の回転部分による軌道と一致する（対称心の付加あり・なし） |
| `Symbols_MatchSpglib` | 標準 Hall 記号、H-M 記号、Schoenflies 記号、点群、晶系 |
| `PointGroupNumber_IsConsistentWithPointGroupSymbol` | 点群番号と点群記号が 1 対 1 に対応する |
| `OrientedPointGroupSymbol_Mm2_MatchesTwofoldAxis` | mm2 型の向き付き点群記号 (2mm/m2m/mm2) が、2 回軸の向きと合う |
| `HallSymbols_IdentifyGroupsUniquely` | 同じ標準 Hall 記号を持つ設定は、同じ群である |
| `WyckoffPositions_MatchCctbx` | Wyckoff 位置の数、多重度、サイト対称。本番の判定 (`GetEquivalentAtomsPosition`) に Wyckoff 位置上の点を動かしてから渡し、多重度・サイト対称・文字が正しく返る |
| `CifSpaceGroupSearchTests` | CIF 読み込みの空間群判定 (`SearchSGseriesNumberForCIF`) が、標準 Hall 記号、旧版の空白なし Hall 記号、H-M 記号のみ、のどれからでも正しい設定を返す。261105Cl 追加：H-M 記号が無く、空白なしの Hall 記号が複数の設定に当たる場合（`P32`、`P62`、`C22-1bc` など 8 通り × 2 設定）に、例外にならず、同じ群なら採り、違う群なら CIF の対称操作か番号で決まる（どちらも無ければ -1）。`ConvertFromCifText` で CIF 全体も読む |
| `BetheEigenReuseTests` | 261105Cl 追加：同じ `BetheMethod` で Beams が置き換えられたあと（STEM・PED など）に平行ビームの振幅を求めても、新しいインスタンスと一致する（固有解を古い並び順のまま使い回さない）。実際に STEM を走らせた後の振幅と HRTEM 像も一致する（STEM が検出器の範囲に限った U' を U のキャッシュに残さない）。STEM の実行中に同じインスタンスで別の条件の回折を計算しても、STEM の結果も回折の結果も変わらない (261105Cl。STEM は専用の複製で走る) |
| `MonoclinicWyckoffItaTests` | 261105Cl 追加：単斜晶 (No. 3–15) の全 57 設定で、Wyckoff 文字が ITA (2016) Vol. A と一致する。b 軸・c 軸の設定は ITA に印刷された特殊位置の座標 (頁つき) と、a 軸の設定と c-ba・ba-c・-acb の設定は ITA の変換 (Table 1.5.1.1、Table 1.5.4.4) で印刷された設定から写したものと比べる。写した設定の対称操作が写す先と一致すること、印刷された c 軸の設定が b 軸の設定を写したものと同じ点の集合であること (書き写しの確認) も確かめる |
| `AffineGenParseTests` | 261105Cl 追加：`AffineGen.Parse` (対称操作の x,y,z 表記) が、読める式 (係数、分数、空白) を正しく読み、読めない式 (小数、大文字、`*`、括弧、分母 0、項の無い符号、空、桁あふれ) は `FormatException` ですぐに止まる (以前は無限ループや黙った誤読) |

## 既知の例外と要確認事項

- **Wyckoff 文字**：cctbx の参照設定（基底変換が恒等）でだけ比べています。非標準設定の文字は基底変換の取り方に依存し、cctbx と spglib でも一致しない場合があるためです。
- **522 (Pn-3m 原点選択 2)**：cctbx の 24i と 24j が ITA と逆です。ReciPro、spglib、pyxtal（Bilbao WYCKPOS 由来）は ITA どおりなので、この設定では文字を比べていません。
- **331 (No. 68, 2bca)**：spglib 2.x は短縮記号を旧記号 `Bbcb` で返します（同じ番号のほかの設定は `Bbeb`）。テストでは `Bbeb` を期待値にしています。
- **H-M 記号だけでは区別できない設定の組**：Cmme (abc) と Cmme (ba-c) などです。CIF の H-M 記号だけのテストでは、組のどれかに割り当てられれば合格としています。
- ~~**要確認（テスト対象外）**：単斜晶の非標準設定 7 つで、ReciPro の Wyckoff 文字が spglib・cctbx の両方と異なります（I12/m1、I112/m、I2/m11、P12/n1、P112/n、P2/n11 の e↔f、I12/c1 の c↔d）。ITA のセル選択の表で確認が必要です。~~ (261105Cl 解決: 下の項)
- **単斜晶の非標準設定の Wyckoff 文字は ITA に合わせています** (261105Cl)。ITA (2016) Vol. A はセル選択 1→2→3 を基底の変換 (a,b,c)→(−a−c, b, a) で結び (Table 1.5.1.1)、b 軸・c 軸の設定のセル選択 2・3 の座標と文字を印刷しています (2002 年版も同じ)。印刷されていない設定 (a 軸の設定と、No. 15 の c-ba・ba-c・-acb の設定) は、ITA の変換 (b 軸→a 軸 (a,b,c)→(b,c,a)、Table 1.5.4.4 の設定記号) で印刷された設定から写しています。同じ多重度・同じサイト対称の 2 つの位置のどちらにどの文字を付けるかは、この基底の変換の取り方で決まります。
  - 261105Cl に ReciPro の表 (`Symmetry/PositionsDictionary.csv` の行の順) を ITA に合わせた 21 設定：No. 12 の A12/m1・I12/m1 (p. 239)、B112/m・I112/m (p. 243)、C2/m11・I2/m11 (写したもの) の c↔d、No. 14 の P12₁/n1・P12₁/a1 (p. 255)、P112₁/n・P112₁/b (p. 259)、P2₁/n11・P2₁/c11 (写したもの) の b↔d、No. 15 の A12/n1 (p. 263)、B112/n (p. 267) と、写した A12/a1・C12/n1・B112/b・A112/n・C2/n11・C2/c11・B2/n11 の c↔d。例えば P12₁/n1 の (0,0,½) は ITA では 2b (以前の ReciPro と cctbx では 2d) です。
  - 基準データ (cctbx) の文字は、27 設定で ITA と異なります：No. 12 の A12/m1・B112/m・C2/m11 (c↔d)、I12/m1・I112/m・I2/m11 (c↔d と e↔f)、No. 13 の P12/n1 (p. 247)・P112/n (p. 251)・P2/n11 (e↔f)、No. 14 の上の 6 設定 (b↔d)、No. 15 の A12/n1・I12/a1 (p. 263)・I12/c1・B112/n・I112/b (p. 267)・A112/n・I112/a・C2/n11・I2/c11・B2/n11・I2/b11 (c↔d)、C12/n1 (a〜d がすべて違う)。以前の確認で spglib も cctbx と同じだったのは I12/m1・I112/m・I2/m11・P12/n1・P112/n・P2/n11 (e↔f) と I12/c1 (c↔d) です。
  - `WyckoffPositions_MatchCctbx` は cctbx の参照設定でだけ文字を比べるので、これらの設定では多重度・サイト対称だけを比べ、文字は `MonoclinicWyckoffItaTests` で ITA と比べます。
