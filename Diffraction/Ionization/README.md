# STEM-EDX イオン化データリソース (260801Cl)

| ファイル | 内容 | リーダー |
|---|---|---|
| `IonizationFsE0.bin` | 内殻イオン化形状因子 F(s,E0) 本番テーブル dataset 1.0.0 (K: Z=6–50 / L1・L23: Z=20–60、計127ch、s=0..4Å⁻¹ 81点、符号付きF) | `IonizationFsTable` (Diffraction/IonizationChannel.cs) |
| `BoteSalvat.bin` | Bote–Salvat 2008 電子衝撃イオン化断面積係数 (Z=1–99, K/L/M subshell) | `BoteSalvat` (同上) |

## 出所とライセンス

- `IonizationFsE0.bin` は**完全自前計算** (DHFS-SCF 生成器 `tools/IonizationGen`、モデルID
  `DHFS-KS23-semi-rel-fullrange-sym-v1`)。OA2000 表・µSTEM データは一切含まれない。
- `BoteSalvat.bin` は usnistgov/BoteSalvatICX.jl (Unlicense) の xione.jl 係数を機械抽出した
  `tools/IonizationGen/bote_full.json` 由来。原典: D. Bote & F. Salvat, PRA 77 (2008) 042701 /
  xion.f (ADNDT 95 (2009) 871)。

## フォーマット・契約の正本

- バイナリレイアウト: `tools/IonizationGen/pack_resource.py` 冒頭コメント
- 補間・範囲・tail の C# 契約: `tools/IonizationGen/prod/MANIFEST.md` (codex 15–16巡)
- 設計: `.project-guidance/ReciPro/ReciPro_STEM-EDX設計.md`

`IonizationFsE0.bin` の method フィールド: 1 = F を float32 可逆格納 (683KB) /
2 = 1e-6 量子化 + s方向delta + byte-plane shuffle (309KB、最大誤差 5.0e-7)。
同梱は **method 2 で確定** (golden-vector 実測 → codex 17巡同意 → 作者確定 2026-08-01。
根拠と既知の符号消失 1 節点は tools/IonizationGen/prod/MANIFEST.md「同梱リソース方式」)。

## 再生成手順

```powershell
cd tools\IonizationGen
python -X utf8 pack_resource.py     # build/ に m1/m2/Bote の 3 バイナリ
python -X utf8 gen_golden.py        # build/golden_v1.json + method_report.md
# 採用する方式の .bin を本フォルダへコピー (IonizationFsE0.m2.bin → IonizationFsE0.bin)
```

検証は `tools/EdxCheck` (golden 照合・破損拒否・恒等式テスト)。
prod テーブル自体を作り直す場合は MANIFEST の dataset_version を上げること (指示書 §3)。
