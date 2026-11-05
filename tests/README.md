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

## 既知の例外と要確認事項

- **Wyckoff 文字**：cctbx の参照設定（基底変換が恒等）でだけ比べています。非標準設定の文字は基底変換の取り方に依存し、cctbx と spglib でも一致しない場合があるためです。
- **522 (Pn-3m 原点選択 2)**：cctbx の 24i と 24j が ITA と逆です。ReciPro、spglib、pyxtal（Bilbao WYCKPOS 由来）は ITA どおりなので、この設定では文字を比べていません。
- **331 (No. 68, 2bca)**：spglib 2.x は短縮記号を旧記号 `Bbcb` で返します（同じ番号のほかの設定は `Bbeb`）。テストでは `Bbeb` を期待値にしています。
- **H-M 記号だけでは区別できない設定の組**：Cmme (abc) と Cmme (ba-c) などです。CIF の H-M 記号だけのテストでは、組のどれかに割り当てられれば合格としています。
- **要確認（テスト対象外）**：単斜晶の非標準設定 7 つで、ReciPro の Wyckoff 文字が spglib・cctbx の両方と異なります（I12/m1、I112/m、I2/m11、P12/n1、P112/n、P2/n11 の e↔f、I12/c1 の c↔d）。ITA のセル選択の表で確認が必要です。
