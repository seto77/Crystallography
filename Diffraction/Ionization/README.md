# STEM-EDX イオン化データリソース (260801Cl 作成 / 260805Cl v3 更新 / 260809Cl v4 更新 / 260811Cl v5 反映)

| ファイル | 内容 | リーダー |
|---|---|---|
| `IonizationFsE0.bin` | 内殻イオン化形状因子 F(s,E0) 本番テーブル dataset 5.0.0 (K: Z=6–50 / L1・L2・L3: Z=20–86 / M1–M5: Z=30–86、計525ch、**s=0..16Å⁻¹ 321点**、符号付きF) | `IonizationFsTable` (Diffraction/IonizationChannel.cs) |
| `BoteSalvat.bin` | Bote–Salvat 2008 電子衝撃イオン化断面積係数 (Z=1–99, K/L/M subshell) | `BoteSalvat` (同上) |

## 出所とライセンス

- `IonizationFsE0.bin` は**完全自前計算** (生成器 = Temari `src/gen_production.jl`、
  モデルID `DHFS-KS23-DiracB-KDIRAC2C-jsplit-fullrange-sym-v4-DSCF` = **κ 分解 Dirac
  連続状態 + 完全 Dirac SCF 原子場** + HIGH 求積 + E0 倍密度グリッド)。
  OA2000 表・µSTEM データは一切含まれない。
- ⚠ **自前計算だが MIT ではない。**Temari は**ソフトを MIT、生成データを CC BY 4.0** と
  意図的に分けている (MIT はソフトウェア向けの文言で EU のデータベース権に触れないため。
  根拠は Temari `licenses/README.md`)。ここに置いてあるのは**その改変版** (下記の再パック) なので、
  **CC BY 4.0 の帰属表示が要る**。表示の実体はリポジトリルートの `THIRD-PARTY-NOTICES.md`。
  - 原典: dataset **v5.0.0** (2026)、DOI [10.5281/zenodo.21872050](https://doi.org/10.5281/zenodo.21872050)
    (版非依存 DOI は 10.5281/zenodo.21872049)
  - **改変あり**: 公開テーブルを `pack_resource.py` (method 2 = 1e-6 量子化 + s 方向 delta +
    byte-plane shuffle) で 1 本のバイナリに詰め直してある。**値の再計算はしていない**
- `BoteSalvat.bin` は usnistgov/BoteSalvatICX.jl (Unlicense) の xione.jl 係数を機械抽出した
  `tools/IonizationGen/bote_full.json` 由来。原典: D. Bote & F. Salvat, PRA 77 (2008) 042701 /
  xion.f (ADNDT 95 (2009) 871)。

## v3.0.0 → v4.0.0 (260809Cl)

1. **M 殻 (M1–M5) を追加**して 246 → 525 チャネル。`shellCode` に **M1..M5 = 5..9**
   が増え、**formatVersion は 3** へ。L23=2 と同じ規律で**番号は再利用しない**
   (2 は欠番のまま)
2. **連続状態を κ 分解 Dirac へ差し替え。**v3 のスカラー相対論連続状態 (SRC) は
   真の相対論効果の 5–20 倍の偽項を持つことが判明した (Temari
   `docs/src_defect_2026-08-07.md`)。あわせて原子場も完全 Dirac SCF になっている
3. **s グリッド・E0 ノード規則・blob 構造は v3 と同一** — リーダー側の追随は
   formatVersion の受け入れと shellCode の上限、そして M 殻の合成だけで済む

⚠ **Z = 30–32 (Zn/Ga/Ge) には M4/M5 の収録が無い** (M1–M3 が 57 本、M4/M5 が 54 本)。
Bote–Salvat の係数表が Z ≤ 32 では 7 副殻 (K, L1–L3, M1–M3) までしか持たないためで、
F テーブルの生成側も同じ条件で弾いている。σ が無い = 重みが作れないので、
`MTotal` は**表にある副殻だけ**を σ 重みで束ねる契約 (`IonizationDataProvider.CodesFor`)。
`LTotal` のように全副殻を要求すると、これらの元素の M 線がまるごと落ちてしまう。

⚠⚠ **σ (Bote–Salvat) 自体の不確かさは、実験との RMS で K 10 % / L 15 % / M 24 %**
(NSRDS 164 / Llovet ら 2014)。**M 殻を出す以上、UI かドキュメントで利用者に伝えること。**
形状 F(s) の精度とは別の話で、こちらは絶対値に効く。

## v4.0.0 → v5.0.0 (260810Cl・**現在の同梱版**)

1. **s グリッドを 161 → 321 点へ延長** (s ≤ 8 → **16 Å⁻¹**、刻みは 0.05 Å⁻¹ のまま)。
   **formatVersion は 4** へ。ALCHEMI は基底が大きいと s = 10.5 Å⁻¹ を要求する
   (β-AlCo・1600 Bloch の実測。**E0 にはほぼ依存しない**) ので、8 Å⁻¹ では足りなかった
2. ★**tail の意味論が変わった。**旧 formatVersion の指数 tail `a·exp(−b·s)` は**撤回**。
   上界でも近似でもなく、43 % の行で hard fail し、高 l では**符号の逆の値**を返していた
   (Au L3 @300 kV, s=12 で +4.15e-5 vs 真値 −1.62e-3)。代わりに行ごとに
   **s_cert (運動学的保証上限)** と **ε (s > s_cert の実測上界)** を持つ
3. ⚠ **formatVersion 1/2/3 の .bin はこのリーダーでは読めない** (s グリッド検査で拒否)。
   `.bin` と `Crystallography.dll` は**必ず同時に差し替える**こと

## フォーマット・契約の正本

- バイナリレイアウト: `tools/IonizationGen/pack_resource.py` 冒頭コメント
- 補間・範囲・tail の C# 契約 + 生成の QC: `tools/IonizationGen/handout/prod_v5_jl/MANIFEST.md`
- 設計: `.project-guidance/ReciPro/ReciPro_STEM-EDX設計.md`

`IonizationFsE0.bin` の method フィールド: 1 = F を float32 可逆格納 (v4 実測 8.25MB) /
2 = 1e-6 量子化 + s方向delta + byte-plane shuffle (v4 実測 2,636,031 bytes、
m2 復元 vs フル精度の最大差 **1.03e-6** = EdxCheck golden 実測。v3 では 7.5e-7 だった —
量子化幅は同じで、M 殻ぶんチャネルが増えて最悪値を引く機会が増えただけ)。
同梱は **method 2 で確定** (golden-vector 実測 → codex 17巡同意 → 作者確定 2026-08-01。
v3/v4 でも同方式を継続。方式選定の根拠は tools/IonizationGen/prod/MANIFEST.md「同梱リソース方式」)。

⚠ **旧 dataset の .bin とは相互に非互換** — リーダーの s グリッド検査
(**SCount 321 / SStep 0.05**) で必ず拒否される。意図した非互換で、混用事故は起きない。
v1/v2 は 81 点 (s ≤ 4 Å⁻¹)、v3/v4 は 161 点 (s ≤ 8 Å⁻¹) なので、**いずれもここで落ちる**
(260810Cl 以前は「v3 は s グリッドが同じなので読めてしまう」状態だったが、
v5 で s グリッド自体が変わったのでその抜け道も閉じた)。

## 再生成手順

```powershell
cd tools\IonizationGen
python -X utf8 pack_resource.py     # build/ に m1/m2/Bote の 3 バイナリ
python -X utf8 gen_golden.py        # golden (EdxCheck 用) を作り直す
# 採用する方式の .bin を本フォルダへコピー (IonizationFsE0.m2.bin → IonizationFsE0.bin)
```

⚠ `BoteSalvat.bin` は**係数が変わらない限り差し替えない**。パッカーを回すと
`packer` 文字列 (git HEAD) だけが変わって別バイトになるが、payload SHA-256 は同一で、
中身は 1 ビットも変わらない。無意味なバイナリ差分を履歴に残さないこと。

検証は `tools/EdxCheck` (golden 照合・破損拒否・恒等式テスト)。
prod テーブル自体を作り直す場合は MANIFEST の dataset_version を上げること (指示書 §3)。
