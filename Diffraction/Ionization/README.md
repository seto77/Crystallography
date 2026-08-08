# STEM-EDX イオン化データリソース (260801Cl 作成 / 260805Cl v3 更新 / 260809Cl v4 更新)

| ファイル | 内容 | リーダー |
|---|---|---|
| `IonizationFsE0.bin` | 内殻イオン化形状因子 F(s,E0) 本番テーブル dataset 4.0.0 (K: Z=6–50 / L1・L2・L3: Z=20–86 / M1–M5: Z=30–86、計525ch、s=0..8Å⁻¹ 161点、符号付きF) | `IonizationFsTable` (Diffraction/IonizationChannel.cs) |
| `BoteSalvat.bin` | Bote–Salvat 2008 電子衝撃イオン化断面積係数 (Z=1–99, K/L/M subshell) | `BoteSalvat` (同上) |

## 出所とライセンス

- `IonizationFsE0.bin` は**完全自前計算** (生成器 = Temari `src/gen_production.jl`、
  モデルID `DHFS-KS23-DiracB-KDIRAC2C-jsplit-fullrange-sym-v4-DSCF` = **κ 分解 Dirac
  連続状態 + 完全 Dirac SCF 原子場** + HIGH 求積 + E0 倍密度グリッド)。
  OA2000 表・µSTEM データは一切含まれない。
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

## フォーマット・契約の正本

- バイナリレイアウト: `tools/IonizationGen/pack_resource.py` 冒頭コメント
- 補間・範囲・tail の C# 契約 + 生成の QC: `tools/IonizationGen/handout/prod_v4_jl/MANIFEST.md`
- 設計: `.project-guidance/ReciPro/ReciPro_STEM-EDX設計.md`

`IonizationFsE0.bin` の method フィールド: 1 = F を float32 可逆格納 (v4 実測 8.25MB) /
2 = 1e-6 量子化 + s方向delta + byte-plane shuffle (v4 実測 2,636,031 bytes、
m2 復元 vs フル精度の最大差 **1.03e-6** = EdxCheck golden 実測。v3 では 7.5e-7 だった —
量子化幅は同じで、M 殻ぶんチャネルが増えて最悪値を引く機会が増えただけ)。
同梱は **method 2 で確定** (golden-vector 実測 → codex 17巡同意 → 作者確定 2026-08-01。
v3/v4 でも同方式を継続。方式選定の根拠は tools/IonizationGen/prod/MANIFEST.md「同梱リソース方式」)。

⚠ **v1/v2 (s=0..4Å⁻¹ 81点) の .bin とは相互に非互換** — リーダーの s グリッド検査
(SCount 161 / SStep 0.05) で必ず拒否される。意図した非互換で、混用事故は起きない。
v3 (formatVersion 2) は s グリッドが同じなので**読めてしまう**が、M 殻が索引に無いので
`HasMShell=false` になり MTotal が UnsupportedShell に落ちる — これも意図した挙動
(旧 dataset で A/B を取るため)。

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
