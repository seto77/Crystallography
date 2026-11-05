# 計算手法ドキュメント / Computational methods

本稿はコードから Claude が起こした下書き（2026-11-05）。著者未査読。
*Draft written by Claude from the source code (2026-11-05). Not yet reviewed by the author.*

各文書は日本語版（`ja/`）と英語版（`en/`）が同名ファイルで対応する。
Each document exists in Japanese (`ja/`) and English (`en/`) under the same file name.

## 回折・動力学 / Diffraction and dynamical theory

| 文書 / Document | 内容 / Scope | 主な実装 / Main code |
|---|---|---|
| scattering-factors ([ja](ja/scattering-factors.md) / [en](en/scattering-factors.md)) | X 線・電子・中性子の原子散乱因子、Temari 因子、Debye–Waller 因子、TDS 吸収形状因子 / atomic scattering factors (X-ray, electron, neutron), Temari factors, Debye–Waller factor, TDS absorptive form factors | `Atom/AtomStatic.cs`, `Atom/Atoms.cs` |
| structure-factor ([ja](ja/structure-factor.md) / [en](en/structure-factor.md)) | 構造因子・運動学的強度・異常分散 / structure factor, kinematical intensity, anomalous dispersion | `Crystal/Crystal.cs` |
| bloch-wave ([ja](ja/bloch-wave.md) / [en](en/bloch-wave.md)) | Bethe（Bloch 波）法: U_g、ビーム選択、固有値問題、CBED・PED・HRTEM・STEM / Bethe (Bloch-wave) method: U_g, beam selection, eigenproblem, CBED, PED, HRTEM, STEM | `Diffraction/BetheMethod.cs`, `WaveProperty.cs`, `DiffractionOptics.cs` |
| bloch-wave-imaging ([ja](ja/bloch-wave-imaging.md) / [en](en/bloch-wave-imaging.md)) | Bloch 波法による像計算: HRTEM・STEM（ADF/TDS）・ポテンシャル像 / image simulation: HRTEM, STEM (ADF/TDS), potential images | `Diffraction/BetheMethod.cs` |
| kikuchi-dynamical ([ja](ja/kikuchi-dynamical.md) / [en](en/kikuchi-dynamical.md)) | Kikuchi 回折の動力学計算、EBSD マスターパターン / dynamical Kikuchi diffraction, EBSD master pattern | `Diffraction/BetheMethod.cs` (EBSD), `Diffraction/Kikuchi/` |
| alchemi-ionization ([ja](ja/alchemi-ionization.md) / [en](en/alchemi-ionization.md)) | 内殻イオン化（STEM-EDX）と ALCHEMI / inner-shell ionization (STEM-EDX) and ALCHEMI | `Diffraction/IonizationChannel.cs`, `Diffraction/Alchemi/`, `AbsorptionSource.cs`, `XrayLineSeries.cs` |
| implementation-questions ([ja](ja/implementation-questions.md) / [en](en/implementation-questions.md)) | 文書作成中に見つけた実装上の疑問点の一覧 / implementation questions found while writing | — |

## 既存の資料 / Existing references

- [Atom/NistElastic/README.md](../../Atom/NistElastic/README.md) — NIST SRD 64 弾性散乱データの出所と再生成（EBSD モンテカルロ用）/ provenance of the NIST SRD 64 elastic-scattering data (EBSD Monte Carlo)
- [Diffraction/Ionization/README.md](../../Diffraction/Ionization/README.md) — イオン化データ（F(s,E0)、Bote–Salvat）の出所・形式・ライセンス / provenance, format and license of the ionization data

## 対象外 / Out of scope

- `Symmetry/`（空間群処理）と `EBSD/`（モンテカルロ・パターン処理）は別文書。EBSD 側の文書は Kikuchi の動力学計算について本ディレクトリの kikuchi-dynamical を参照する。
- `Symmetry/` (space-group handling) and `EBSD/` (Monte Carlo, pattern processing) are documented separately; the EBSD documents refer to kikuchi-dynamical here for the dynamical Kikuchi calculation.
