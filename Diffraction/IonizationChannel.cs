#region using
using System;
using System.Collections.Generic;
using System.IO;
using System.IO.Compression;
using System.Linq;
using System.Numerics;
using System.Reflection;
using System.Security.Cryptography;
using System.Text;
#endregion

namespace Crystallography;

//260801Cl 追加: STEM-EDX (内殻イオン化) チャネルの型群・F(s,E0) テーブル・補間器。
//設計正本 = .project-guidance/ReciPro/ReciPro_STEM-EDX設計.md §5.1、データ契約 = tools/IonizationGen/prod/MANIFEST.md (codex 15-16巡)。
//テーブルは完全自前計算 (DHFS-KS23-semi-rel-fullrange-sym-v1)。OA2000/µSTEM のデータは一切含まれない。

#region 公開 enum / record (設計書 §5.1)

/// <summary>260801Cl 追加: イオン化殻。v1 のプロバイダが返すのは K / LTotal のみ (L1/L2/L3 は v2 で分離)。
/// 260802Cl 変更: **値を明示** (旧: 暗黙の 0,1,2,3,4)。ReciPro のプリセットが (Z, Shell) の組を
/// この基底値のまま永続化する (ImageSimulatorSetting.EdxChannels、設計書 §5.9.1-6) ので、
/// 既存の値は変更・並べ替え禁止。新しい殻は末尾に追加すること。</summary>
public enum IonizationShell { K = 0, LTotal = 1, L1 = 2, L2 = 3, L3 = 4 }

/// <summary>260801Cl 追加: 元素×殻のチャネル指定。</summary>
public record IonizationChannelSpec(int Z, IonizationShell Shell)
{
    /// <summary>260801Cl 追加: 短い表示名 (例 "Fe-K" / "Sr-L")。要約・凡例・チャネル選択 UI で共通に使う
    /// (呼び出し側で殻の三項演算子を書き散らさない)。</summary>
    public string ShortLabel => $"{AtomStatic.AtomicName(Z)}-{(Shell == IonizationShell.LTotal ? "L" : Shell.ToString())}";

    /// <summary>260801Cl 追加: 原子番号と殻を明示した表示名 (例 "Fe (26) K" / "Sr (38) L (total)")。</summary>
    public string Label => $"{AtomStatic.AtomicName(Z)} ({Z}) {(Shell == IonizationShell.LTotal ? "L (total)" : Shell.ToString())}";
}

/// <summary>260801Cl 追加: データ出所 (σ と形状で分離して保持する)。</summary>
public sealed record IonizationDataProvenance(string ModelId, string DatasetVersion, string Detail);

/// <summary>260801Cl 追加: 正規化イオン化形状 F(s)。F(0)=1。s の単位は nm⁻¹ (s=|G|/2)。batch 評価 (N² 内の virtual call 回避)。</summary>
public interface INormalizedIonizationShape
{
    void Evaluate(ReadOnlySpan<double> sPerNm, Span<double> values);
}

/// <summary>260801Cl 追加: run 開始時に immutable へ解決されたチャネルデータ (プロバイダ選択と範囲判定を実行中に持ち込まない)。</summary>
public sealed record IonizationData(
    IonizationChannelSpec Target,
    double EdgeEnergyKeV,                 // Bote/xion edge (LTotal は最小 subshell edge)
    double TotalCrossSectionNm2,          // Bote–Salvat (LTotal は開いている subshell の合算)
    INormalizedIonizationShape Shape,     // F(0)=1
    IonizationDataProvenance CrossSectionSource,
    IonizationDataProvenance ShapeSource);

/// <summary>260801Cl 追加: 物理シグナル量。深さ分解自己吸収 (v3) の伏線として v1 から区別する (設計書 §5.5)。</summary>
public enum SignalQuantity { IonizationVacanciesGenerated, XrayPhotonsGenerated, XrayPhotonsSelfAbsorbed, XrayPhotonsDetected }

/// <summary>260801Cl 追加: モデル上の規格化状態 (表示正規化と混同しない)。</summary>
public enum SignalNormalization { ModelAbsoluteNotAudited, PerIncidentElectron }

/// <summary>260801Cl 追加: 表示正規化 (GUI 専用)。</summary>
public enum DisplayNormalization { PerMaximum, Absolute }

/// <summary>260801Cl 追加: RunSTEM への EDX 要求 (多チャネルは配列で渡す。同一 (Z,Shell) の重複は hard error、codex 20巡)。
/// HermitianTolerance は ±q 非 Hermitian 残差 (相対) の許容値。超過時は対称化せず hard fail (設計書 §3.4)。
/// 0.01 を上限として「厳しくする方向のみ」有効 (それ以上は 0.01 に clamp、codex 17巡)。非有限・負値は run 前に hard error。
/// 残差は方向グリッド div に対しほぼ O(h²) (bilinear 補間誤差由来): div=10 で ~0.11 / 32 で ~0.009 / 48 で ~0.0017 実測。
/// CaptureRawIq=true で ±q 対称化前の I(q,t,d) を StemSignalMap.IqBeforeSymmetrization に保持 (fixture 凍結用、設計書 §6.3。通常 run は false)。</summary>
public sealed record StemIonizationRequest(IonizationChannelSpec Channel, double HermitianTolerance = 0.01, bool CaptureRawIq = false)//260801Cl CaptureRawIq 追加 (旧: (Channel, HermitianTolerance) の 2 引数)
{
    /// <summary>260801Cl 追加: ±q Hermitian 残差が許容 0.01 に十分収まる方向グリッド分割数の推奨下限。
    /// 残差は div に対しほぼ O(h²) で、実測 div=32 → 0.0093 (許容ぎりぎり) / 48 → 0.0017。
    /// GUI の事前警告もこの値を参照する (数値の根拠がドメイン側にあるため、しきい値も UI に置かない)。</summary>
    public const int RecommendedProbeDivision = 48;
}

//260801Cl 削除: StemEdxResult (v0a の 1 チャネル内部形) は StemSimulationResult/StemSignalMap (§5.5、codex 20巡) に置換。旧定義は git 履歴 (Crystallography 4e0f39e) 参照。

/// <summary>260801Cl 追加: STEM 画像スタック (storage-neutral な公開型、codex 20巡)。
/// v0b の backing は既存 jagged [t][d][pix] (byte-exact 接続コスト優先)。将来連続 double[] へ変えても公開 API は不変。</summary>
public sealed class StemImageStack
{
    private readonly double[][][] _planes; // [thickness][defocus][pixel]
    public System.Drawing.Size Size { get; }
    public int ThicknessCount => _planes.Length;
    public int DefocusCount => _planes.Length > 0 ? _planes[0].Length : 0;

    internal StemImageStack(System.Drawing.Size size, double[][][] planes) { Size = size; _planes = planes; }

    public ReadOnlyMemory<double> GetPlane(int thicknessIndex, int defocusIndex) => _planes[thicknessIndex][defocusIndex];

    /// <summary>260801Cl 変更: 内部 backing を公開 (旧: internal Backing)。legacy ResultSTEM タプルが同じ配列を
    /// 既に GUI へ渡しているので新たな露出は無く、これが internal だと消費側が全員 GetPlane→ToArray の
    /// 全画素コピーを書く羽目になっていた (GUI は再描画ごと)。**公開後は変更禁止**の契約は ResultSTEM と同じ。</summary>
    public double[][][] Planes => _planes;
}

/// <summary>260801Cl 追加: 1 チャネル分の STEM-EDX 信号マップ (設計書 §5.5)。公開後は不変。</summary>
public sealed class StemSignalMap
{
    public IonizationData Data { get; init; }
    /// <summary>チャネル指定 (Data.Target から導出、二重保存を避ける)</summary>
    public IonizationChannelSpec Channel => Data.Target;
    public SignalQuantity Quantity { get; init; }
    public SignalNormalization Normalization { get; init; }
    /// <summary>実空間マップ (±q 対称化・実部合成・負値 clamp 済み)</summary>
    public StemImageStack Image { get; init; }
    /// <summary>±q 対称化前の非 Hermitian 残差最大値 (相対)</summary>
    public double HermitianResidualMax { get; init; }
    /// <summary>この run に適用された許容値 (clamp 後)</summary>
    public double HermitianToleranceApplied { get; init; }
    /// <summary>q=0 の虚部残差最大値 (相対)</summary>
    public double QZeroImagMax { get; init; }
    /// <summary>clamp 前の最小画素値 (負値診断)</summary>
    public double MinPixelBeforeClamp { get; init; }
    /// <summary>形状評価が s>4 Å⁻¹ の tail 外挿を使ったか (診断フラグ、silent extrapolation 禁止の契約)</summary>
    public bool UsedTailExtrapolation { get; init; }
    /// <summary>対称化後の I(q,t,d) (検証・回帰用の生値。設計書 §6.2「保存対象は最終画像だけでなく」)</summary>
    public Complex[,,] Iq { get; init; }
    /// <summary>±q 対称化前の I(q,t,d) (CaptureRawIq=true の run のみ、通常 run は null。fixture 用 §6.3 — 対称化後だけでは位相符号・共役ミスが隠れるため)</summary>
    public Complex[,,] IqBeforeSymmetrization { get; init; }
}

/// <summary>260801Cl 追加: STEM run の primary 結果 (設計書 §5.5、codex 20巡)。worker 終端で一度だけ Volatile.Write で公開され、以後不変。
/// 失敗・cancel した run は公開されない (以前の成功結果が残る)。legacy ResultSTEM タプルはこの型からの互換 view。</summary>
public sealed class StemSimulationResult
{
    public long RunId { get; init; }
    public System.Drawing.Size Size { get; init; }
    public double Resolution { get; init; }
    public double[] Thicknesses { get; init; }
    public double[] Defocusses { get; init; }
    public Matrix3D Rotation { get; init; }
    public StemImageStack ImageEla { get; init; }
    public StemImageStack ImageTDS { get; init; }
    public StemImageStack ImageBoth { get; init; }
    /// <summary>EDX 信号 (要求順を保持。EDX off の run では空配列、null にはならない)</summary>
    public StemSignalMap[] EdxSignals { get; init; }
    /// <summary>EDX の q 次元に対応する combined hkl (全チャネル共通。EDX off では空配列)</summary>
    public (int H, int K, int L)[] QIndices { get; init; }
    /// <summary>各 q のエントリ数 (aperture 重なりで有効だった方向数。0 = 未計算。EDX off では空配列)</summary>
    public int[] QEntryCounts { get; init; }

    /// <summary>260802Cl 追加: 参照像 (弾性・TDS) の数値品質。設計書 §3.4 の規律を EDX 以外にも広げた際の実測値。
    /// 像 I(r) は実数なので I(−q)=I(q)* が厳密に成り立つべきで、破れは方向グリッド上の bilinear 補間による
    /// O(h²) の数値誤差にすぎない。**大きいほど角度分解能 (probe division) が粗い**ことを意味する。
    /// EDX と違い run は止めない (欠陥ではなくユーザーの設定なので、報告して判断に委ねる)。</summary>
    public StemReferenceQuality ReferenceQuality { get; init; }
}

/// <summary>260802Cl 追加: 参照像の数値品質診断 (値そのものには影響しない観測量)。</summary>
/// <param name="ElasticHermitianMax">弾性 I(q) の非 Hermitian 残差 (相対)</param>
/// <param name="ElasticQZeroImagMax">弾性 q=0 の虚部残差 (相対)</param>
/// <param name="TdsHermitianMax">TDS I(q) の非 Hermitian 残差 (相対)</param>
/// <param name="TdsQZeroImagMax">TDS q=0 の虚部残差 (相対)</param>
/// <param name="ElasticMinPixelBeforeClamp">弾性像の clamp 前最小画素 (負なら打切り誤差の目安)</param>
/// <param name="TdsMinPixelBeforeClamp">TDS 像の clamp 前最小画素</param>
public sealed record StemReferenceQuality(
    double ElasticHermitianMax, double ElasticQZeroImagMax,
    double TdsHermitianMax, double TdsQZeroImagMax,
    double ElasticMinPixelBeforeClamp, double TdsMinPixelBeforeClamp,
    double ElasticImagOverRealMax, double TdsImagOverRealMax);

/// <summary>260802Cl 追加: STEM run の進捗ステージ (設計書 §5.9-8)。値は GUI の "Stage n" 表示と一致させてある。</summary>
public enum StemStage { EigenSolve = 1, ElasticQ = 2, PotentialMatrix = 3, InelasticQ = 4, IonizationQ = 5 }

/// <summary>260802Cl 追加: STEM の進捗通知 (<see cref="System.ComponentModel.ProgressChangedEventArgs.UserState"/> に載せる。設計書 §5.9-8)。
/// v0b までは "Calculating I_EDX(Q) (ch 2/3)" のような文字列前方一致プロトコルだったが、ステージが増えるたびに
/// 前方一致の分岐が増え、チャネル番号を文字列からパースする必要もあった (負長 Substring でクラッシュし得た) ので型に置き換えた。
/// 表示文言は GUI 側で組む — backend は表示都合の文字列を持たない。</summary>
/// <param name="Stage">ステージ</param>
/// <param name="Fraction">当該ステージ内の進捗 0–1。<c>ProgressPercentage</c> は同じ値の 1E6 倍なので両者は必ず一致する</param>
/// <param name="SolverLabel">Stage1 のみ: 実際に使った solver とスレッド数 (例 "Eigen8")。それ以外は null</param>
/// <param name="ChannelIndex">Stage5 のみ: 計算中の EDX チャネル (0 始まり)。それ以外は -1</param>
/// <param name="ChannelCount">この run が要求された EDX チャネル数 (EDX なしは 0)。**全ステージで有効** = GUI は Stage4 の
/// 進捗配分を Stage5 に入る前から決められる (旧実装が GUI 側のチェック状態を見に行かずに済む)</param>
/// <param name="Channel">Stage5 のみ: 計算中のチャネル指定。それ以外は null</param>
public sealed record StemProgressInfo(StemStage Stage, double Fraction, string SolverLabel = null,
    int ChannelIndex = -1, int ChannelCount = 0, IonizationChannelSpec Channel = null);

/// <summary>260801Cl 追加: チャネル利用可否 (GUI は本 enum を 11 言語の表示文へ変換する。例外文字列を UI に直接出さない。設計書 §5.9-3)</summary>
public enum IonizationAvailability { Available, BelowEdge, UnsupportedShell, UnsupportedElement, E0OutOfRange }

/// <summary>260801Cl 追加: GUI 向けチャネル照会結果 (設計書 §5.9-3。GUI 側に Z 範囲をハードコードしない)。
/// Status==Available は「同じ引数で Resolve が成功する」と同値 (codex 20巡の契約)。未定義値は NaN。</summary>
public sealed record IonizationChannelInfo
{
    /// <summary>260801Cl 追加: 照会したチャネル (呼び出し側が spec と info を対で持ち回らなくて済む)</summary>
    public IonizationChannelSpec Channel { get; init; }
    public IonizationAvailability Status { get; init; }
    /// <summary>吸収端 [keV] (テーブル収録があれば E0 範囲外でも返す。UnsupportedElement/UnsupportedShell では NaN)</summary>
    public double EdgeEnergyKeV { get; init; } = double.NaN;
    /// <summary>過電圧 U = E0/E_edge (edge が取れれば返す)</summary>
    public double Overvoltage { get; init; } = double.NaN;
    /// <summary>総イオン化断面積 [nm²] (Available のみ。それ以外は NaN)</summary>
    public double SigmaNm2 { get; init; } = double.NaN;
    /// <summary>σ の出所 (provider 品質タグ、§5.6)</summary>
    public IonizationDataProvenance CrossSectionSource { get; init; }
    /// <summary>形状の出所 (provider 品質タグ、§5.6)</summary>
    public IonizationDataProvenance ShapeSource { get; init; }

    /// <summary>260801Cl 追加: 過電圧 U がこの値を下回ると断面積の信頼度が落ちる (選択は可能だが警告する)。</summary>
    public const double LowOvervoltage = 1.2;

    /// <summary>選択に注意が要るか (Available だが U が低い、または選択不可)。</summary>
    public bool HasCaution => Status != IonizationAvailability.Available || Overvoltage < LowOvervoltage;

    //260801Cl 追加 (作者指示 2026-08-01: GUI 側に置いていた 11 言語の状態文を、状態 enum と同じ場所へ移す)。
    //§5.6 の「変換が物理的に許されるかの判定・provenance・単位は Crystallography.dll が持つ」に沿う。
    //短文 (一覧の末尾に括弧書き) と長文 (ToolTip) を 1 つの switch で対にし、しきい値も 1 か所しか書かない。
    private (string Short, string Long) StatusText() => Status switch
    {
        IonizationAvailability.Available when Overvoltage < LowOvervoltage => (
            Localization.Loc(en: "low U", ja: "U 小", de: "U klein", fr: "U faible", es: "U baja", pt: "U baixa",
                it: "U bassa", ru: "малое U", zhHans: "U 偏低", zhHant: "U 偏低", ko: "U 낮음"),
            Localization.Loc(
                en: "Available, but the overvoltage U = E0/E_edge is below 1.2 — the cross section is less reliable there.",
                ja: "利用可能ですが過電圧 U = E0/E_edge が 1.2 未満です。この領域は断面積の信頼度が下がります。",
                de: "Verfügbar, aber die Überspannung U = E0/E_Kante liegt unter 1,2 — der Wirkungsquerschnitt ist dort weniger zuverlässig.",
                fr: "Disponible, mais le survoltage U = E0/E_seuil est inférieur à 1,2 : la section efficace y est moins fiable.",
                es: "Disponible, pero la sobretensión U = E0/E_borde es menor que 1,2: la sección eficaz es menos fiable ahí.",
                pt: "Disponível, mas a sobretensão U = E0/E_borda é inferior a 1,2: a secção eficaz é menos fiável aí.",
                it: "Disponibile, ma la sovratensione U = E0/E_soglia è inferiore a 1,2: la sezione d'urto è meno affidabile.",
                ru: "Доступно, но перенапряжение U = E0/E_края меньше 1,2 — сечение там менее надёжно.",
                zhHans: "可用，但过电压 U = E0/E_边 低于 1.2，该区域的截面可靠性较低。",
                zhHant: "可用，但過電壓 U = E0/E_邊 低於 1.2，該區域的截面可靠性較低。",
                ko: "사용 가능하지만 과전압 U = E0/E_단 이 1.2 미만입니다. 이 영역은 단면적 신뢰도가 낮습니다.")),
        IonizationAvailability.Available => ("",
            Localization.Loc(en: "Available.", ja: "利用可能です。", de: "Verfügbar.", fr: "Disponible.", es: "Disponible.",
                pt: "Disponível.", it: "Disponibile.", ru: "Доступно.", zhHans: "可用。", zhHant: "可用。", ko: "사용 가능합니다.")),
        IonizationAvailability.BelowEdge => (
            Localization.Loc(en: "below edge", ja: "端以下", de: "unter Kante", fr: "sous seuil", es: "bajo borde",
                pt: "sob borda", it: "sotto soglia", ru: "ниже края", zhHans: "低于边", zhHant: "低於邊", ko: "단 이하"),
            Localization.Loc(
                en: "The incident energy is below the absorption edge, so this shell cannot be ionized.",
                ja: "入射エネルギーが吸収端より低いため、この殻はイオン化されません。",
                de: "Die Primärenergie liegt unter der Absorptionskante, diese Schale kann nicht ionisiert werden.",
                fr: "L'énergie incidente est inférieure au seuil d'absorption : cette couche ne peut pas être ionisée.",
                es: "La energía incidente está por debajo del borde de absorción: esta capa no puede ionizarse.",
                pt: "A energia incidente está abaixo da borda de absorção: esta camada não pode ser ionizada.",
                it: "L'energia incidente è sotto la soglia di assorbimento: questo guscio non può essere ionizzato.",
                ru: "Энергия пучка ниже края поглощения, поэтому эта оболочка не ионизируется.",
                zhHans: "入射能量低于吸收边，该壳层无法被电离。",
                zhHant: "入射能量低於吸收邊，該殼層無法被游離。",
                ko: "입사 에너지가 흡수단보다 낮아 이 껍질은 이온화되지 않습니다.")),
        IonizationAvailability.E0OutOfRange => (
            Localization.Loc(en: "E0 range", ja: "E0 範囲外", de: "E0-Bereich", fr: "plage E0", es: "rango E0", pt: "faixa E0",
                it: "intervallo E0", ru: "диапазон E0", zhHans: "E0 范围", zhHant: "E0 範圍", ko: "E0 범위"),
            Localization.Loc(
                en: "STEM-EDX supports 30-400 kV only (the ionization form-factor table is not extrapolated).",
                ja: "STEM-EDX は 30-400 kV のみ対応です (イオン化形状因子テーブルを外挿しないため)。",
                de: "STEM-EDX unterstützt nur 30-400 kV (die Tabelle der Ionisationsformfaktoren wird nicht extrapoliert).",
                fr: "STEM-EDX ne prend en charge que 30-400 kV (la table des facteurs de forme d'ionisation n'est pas extrapolée).",
                es: "STEM-EDX solo admite 30-400 kV (la tabla de factores de forma de ionización no se extrapola).",
                pt: "O STEM-EDX suporta apenas 30-400 kV (a tabela de fatores de forma de ionização não é extrapolada).",
                it: "STEM-EDX supporta solo 30-400 kV (la tabella dei fattori di forma di ionizzazione non viene estrapolata).",
                ru: "STEM-EDX поддерживает только 30-400 кВ (таблица форм-факторов ионизации не экстраполируется).",
                zhHans: "STEM-EDX 仅支持 30-400 kV（电离形状因子表不做外推）。",
                zhHant: "STEM-EDX 僅支援 30-400 kV（游離形狀因子表不做外推）。",
                ko: "STEM-EDX 는 30-400 kV 만 지원합니다 (이온화 형상 인자 표를 외삽하지 않음).")),
        _ => ("", "")
    };

    /// <summary>一覧の 1 行分の表示文 (例 <c>O (8) K   0.537 keV   U = 372</c>)。注意が要る状態は末尾に括弧書き。</summary>
    public string ToListItemText()
    {
        var text = Channel.Label;
        if (!double.IsNaN(EdgeEnergyKeV)) text += $"   {EdgeEnergyKeV:f3} keV";
        if (!double.IsNaN(Overvoltage)) text += $"   U = {Overvoltage.ToString(Overvoltage < 100 ? "f2" : "f0")}";
        var status = StatusText().Short;
        return status.Length == 0 ? text : $"{text}   ({status})";
    }

    /// <summary>状態の完全な説明 + provider 品質タグ (ToolTip 用)。</summary>
    public string ToDescription()
    {
        var head = StatusText().Long;
        return ShapeSource is null ? head : $"{head}\r\nσ: {CrossSectionSource.ModelId} / F(s): {ShapeSource.ModelId} {ShapeSource.DatasetVersion}";
    }
}

#endregion

#region scipy 互換 PCHIP (260801Cl 追加)

/// <summary>260801Cl 追加: scipy.interpolate.PchipInterpolator 互換の単調 3 次 Hermite 補間。
/// 導関数は Fritsch–Carlson 加重調和平均 + scipy 流エッジ処理 (_edge_case)。評価は PPoly と同じ
/// 局所 power 基底 Horner。Python golden vector (tools/IonizationGen/build/golden_v1.json) とロックする。</summary>
public static class Pchip
{
    /// <summary>ノード導関数 (scipy _find_derivatives 互換)。x は狭義単調増加、n≥2。</summary>
    public static double[] Derivatives(double[] x, double[] y)
    {
        int n = x.Length;
        var d = new double[n];
        if (n == 2)
        {
            d[0] = d[1] = (y[1] - y[0]) / (x[1] - x[0]);
            return d;
        }
        var h = new double[n - 1];
        var m = new double[n - 1];
        for (int k = 0; k < n - 1; k++)
        {
            h[k] = x[k + 1] - x[k];
            m[k] = (y[k + 1] - y[k]) / h[k];
        }
        for (int k = 1; k < n - 1; k++)
        {
            if (Math.Sign(m[k]) != Math.Sign(m[k - 1]) || m[k] == 0 || m[k - 1] == 0)
                d[k] = 0;
            else
            {
                double w1 = 2 * h[k] + h[k - 1], w2 = h[k] + 2 * h[k - 1];
                var whmean = (w1 / (w1 + w2)) / m[k - 1] + (w2 / (w1 + w2)) / m[k];
                d[k] = 1.0 / whmean;
            }
        }
        d[0] = EdgeCase(h[0], h[1], m[0], m[1]);
        d[n - 1] = EdgeCase(h[n - 2], h[n - 3], m[n - 2], m[n - 3]);
        return d;
    }

    private static double EdgeCase(double h0, double h1, double m0, double m1)
    {
        var d = ((2 * h0 + h1) * m0 - h0 * m1) / (h0 + h1);
        if (Math.Sign(d) != Math.Sign(m0)) return 0.0;
        if (Math.Sign(m0) != Math.Sign(m1) && Math.Abs(d) > 3.0 * Math.Abs(m0)) return 3.0 * m0;
        return d;
    }

    /// <summary>1 点評価。範囲外は端区間の 3 次式で外挿 (scipy extrapolate=True 相当)。</summary>
    public static double Evaluate(double[] x, double[] y, double[] d, double xq)
    {
        int n = x.Length;
        int i = Array.BinarySearch(x, xq);
        if (i < 0) i = ~i - 1;
        i = Math.Clamp(i, 0, n - 2);
        double hh = x[i + 1] - x[i], slope = (y[i + 1] - y[i]) / hh;
        double c0 = (d[i] + d[i + 1] - 2 * slope) / (hh * hh);
        double c1 = (3 * slope - 2 * d[i] - d[i + 1]) / hh;
        double s = xq - x[i];
        return ((c0 * s + c1) * s + d[i]) * s + y[i];
    }
}

#endregion

#region Bote–Salvat 断面積 (260801Cl 追加)

/// <summary>260801Cl 追加: Bote–Salvat 2008 電子衝撃イオン化総断面積 (K/L/M subshell, Z=1–99)。
/// 移植元 = usnistgov/BoteSalvatICX.jl (Unlicense) / xion.f (Bote, Salvat, Jablonski, Powell, ADNDT 95 (2009) 871)。
/// 係数は埋め込みリソース Crystallography.BoteSalvat.bin (tools/IonizationGen/pack_resource.py が bote_full.json から生成)。
/// Python 参照実装 = tools/IonizationGen/botesalvat.py (golden vector で照合)。</summary>
public static class BoteSalvat
{
    private const string ResourceName = "Crystallography.BoteSalvat.bin"; // csproj の LogicalName と一致させること
    private const int Magic = 0x45544F42; // "BOTE"
    private const double Rev = 5.10998918e5;   // 電子静止エネルギー [eV] (xion.f と同値)
    private const double A0Cm = 5.291772108e-9; // Bohr 半径 [cm]

    private sealed class Element
    {
        public double[] Be, Anlj, G, EdgeEv, A; // G は [nss*4]、A は [nss*5] row-major
    }

    private static volatile Element[] _elements; // [z-1]、非 null になった時点で全構築済み (NistElastic と同じ volatile+lock 公開)
    private static readonly object _sync = new();

    private static Element[] Load()
    {
        var el = _elements;
        if (el is not null) return el;
        lock (_sync)
        {
            if (_elements is not null) return _elements;
            using var stream = Assembly.GetExecutingAssembly().GetManifestResourceStream(ResourceName)
                ?? throw new InvalidOperationException($"Embedded resource not found: {ResourceName}");
            using var reader = new BinaryReader(stream);
            if (reader.ReadInt32() != Magic) throw new InvalidDataException("BoteSalvat.bin: bad magic");
            reader.ReadInt32(); // formatVersion
            var codec = reader.ReadInt32();
            if (codec != 1) throw new InvalidDataException($"BoteSalvat.bin: unknown codec {codec}");
            ReadString(reader); // source_ref
            ReadString(reader); // packer
            var sha = reader.ReadBytes(32);
            var compLen = reader.ReadInt32();
            var comp = reader.ReadBytes(compLen);
            using var ms = new MemoryStream(comp, writable: false);
            using var br = new BrotliStream(ms, CompressionMode.Decompress);
            using var payload = new MemoryStream();
            br.CopyTo(payload);
            var raw = payload.ToArray();
            if (!SHA256.HashData(raw).AsSpan().SequenceEqual(sha))
                throw new InvalidDataException("BoteSalvat.bin: payload SHA-256 mismatch");
            using var pr = new BinaryReader(new MemoryStream(raw, writable: false));
            var zCount = pr.ReadInt32();
            var arr = new Element[100];
            for (int i = 0; i < zCount; i++)
            {
                int z = pr.ReadInt32(), nss = pr.ReadInt32();
                var e = new Element
                {
                    Be = ReadDoubles(pr, nss),
                    Anlj = ReadDoubles(pr, nss),
                    G = ReadDoubles(pr, nss * 4),
                    EdgeEv = ReadDoubles(pr, nss),
                    A = ReadDoubles(pr, nss * 5),
                };
                arr[z] = e;
            }
            _elements = arr;
            return arr;
        }
    }

    internal static string ReadString(BinaryReader r)
    {
        var len = r.ReadInt32();
        return Encoding.UTF8.GetString(r.ReadBytes(len));
    }

    private static double[] ReadDoubles(BinaryReader r, int count)
    {
        var a = new double[count];
        for (int i = 0; i < count; i++) a[i] = r.ReadDouble();
        return a;
    }

    private static Element Get(int z)
        => (uint)z <= 99 && Load()[z] is Element e ? e
           : throw new ArgumentOutOfRangeException(nameof(z), z, "Bote–Salvat coefficients cover Z=1–99");

    /// <summary>Z の収録サブシェル数 (1〜9)。index 順 = K, L1, L2, L3, M1..M5。</summary>
    public static int SubshellCount(int z) => Get(z).EdgeEv.Length;

    /// <summary>吸収端エネルギー [eV]。subshell は 1 始まり (1=K, 2=L1, 3=L2, 4=L3, 5..9=M1..M5)。</summary>
    public static double EdgeEv(int z, int subshell) => Get(z).EdgeEv[CheckSubshell(z, subshell)];

    private static int CheckSubshell(int z, int subshell)
        => subshell >= 1 && subshell <= SubshellCount(z) ? subshell - 1
           : throw new ArgumentOutOfRangeException(nameof(subshell), subshell, $"Z={z} has {SubshellCount(z)} subshells");

    /// <summary>イオン化断面積 [cm²]。演算順は botesalvat.py sigma_cm2 と厳密一致 (golden vector 照合)。</summary>
    public static double SigmaCm2(int z, int subshell, double energyEv, double edgeEvOverride = double.NaN)
    {
        var el = Get(z);
        int ss = CheckSubshell(z, subshell);
        var edge = double.IsNaN(edgeEvOverride) ? el.EdgeEv[ss] : edgeEvOverride;
        var overv = energyEv / edge;
        if (overv <= 1.0) return 0.0;
        double xione;
        if (overv <= 16.0)
        {
            double a1 = el.A[ss * 5], a2 = el.A[ss * 5 + 1], a3 = el.A[ss * 5 + 2], a4 = el.A[ss * 5 + 3], a5 = el.A[ss * 5 + 4];
            var opu = 1.0 / (1.0 + overv);
            var ffitlo = a1 + a2 * overv + opu * (a3 + opu * opu * (a4 + opu * opu * a5));
            var r = ffitlo / overv;
            xione = (overv - 1.0) * (r * r);
        }
        else
        {
            var beta2 = (energyEv * (energyEv + 2.0 * Rev)) / ((energyEv + Rev) * (energyEv + Rev));
            var x = Math.Sqrt(energyEv * (energyEv + 2.0 * Rev)) / Rev;
            double g1 = el.G[ss * 4], g2 = el.G[ss * 4 + 1], g3 = el.G[ss * 4 + 2], g4 = el.G[ss * 4 + 3];
            var ffitup = (2.0 * Math.Log(x) - beta2) * (1.0 + g1 / x) + g2
                + g3 * Math.Sqrt(Rev / (energyEv + Rev)) + g4 / x;
            var factr = el.Anlj[ss] / beta2;
            xione = ((factr * overv) / (overv + el.Be[ss])) * ffitup;
        }
        return 4.0 * Math.PI * (A0Cm * A0Cm) * xione;
    }

    /// <summary>イオン化断面積 [nm²]。</summary>
    public static double SigmaNm2(int z, int subshell, double energyEv, double edgeEvOverride = double.NaN)
        => SigmaCm2(z, subshell, energyEv, edgeEvOverride) * 1e14;
}

#endregion

#region F(s,E0) テーブル (260801Cl 追加)

/// <summary>260801Cl 追加: 本番 F(s,E0) テーブルのリーダー。
/// フォーマット・契約 = tools/IonizationGen/pack_resource.py ヘッダコメント + prod/MANIFEST.md。
/// NistElasticPchipResource と同じ「blob 常駐 + チャネル単位 lazy Brotli decode + volatile 公開」。
/// 260802Cl: formatVersion 2 (dataset 2.0.0) で 2p が j 分離され L2/L3 の shellCode が増えた。
/// v1 (formatVersion 1, L23 1 本) の .bin も読める — 生成側の A/B 比較に使うため。</summary>
public sealed class IonizationFsTable
{
    private const string ResourceName = "Crystallography.IonizationFsE0.bin"; // csproj の LogicalName と一致させること
    private const int Magic = 0x31534649; // "IFS1"
    //260802Cl: L2=3 / L3=4 を追加。**L23=2 は欠番として予約**し番号を再利用しない
    //(v1 の .bin を読んだときに意味が入れ替わらないようにするため)。
    public const int ShellCodeK = 0, ShellCodeL1 = 1, ShellCodeL23 = 2, ShellCodeL2 = 3, ShellCodeL3 = 4;
    //260805Cl 変更: 4.0 → 8.0 (dataset v3.0.0 で s グリッドを 81 → 161 点へ延長)。
    //正典の SrTiO₃ 条件 (a=0.3905nm, 125 beams, 200kV) が実際に要求する s は
    //max|q+g_i−g_j|/2 = 5.56 Å⁻¹ で、s≤4 では全行列要素の 5.5 % が tail 外挿頼みだった。
    //v1/v2 の .bin (SCount=81) はこのリーダーでは読めない (下の s グリッド検査で拒否)。
    public const double SMaxAngstromInv = 8.0;
    //public const double SMaxAngstromInv = 4.0;   //260804Cl まで (dataset v1/v2)

    /// <summary>260802Cl 追加: 2p が j 分離された dataset か (L2/L3 が索引にある = v2)。
    /// LTotal の合成を「L1+L2+L3」にするか「L1+L23」にするかの分岐に使う。</summary>
    public bool HasJResolvedL { get; }

    public int Method { get; }             // 1=float32 / 2=1e-6 量子化+delta+shuffle
    public int SCount { get; }             // 81
    public double SStep { get; }           // 0.05 Å⁻¹
    public string DatasetVersion { get; }
    public string ModelId { get; }
    public string BoteRef { get; }
    public string Packer { get; }
    public byte[] PayloadSha256 { get; }

    private readonly byte[] _blob;
    private readonly int _payloadStart;
    private readonly Dictionary<(int ShellCode, int Z), (int Offset, int Length)> _index = [];
    private readonly Dictionary<(int ShellCode, int Z), IonizationChannelTable> _cache = [];
    private readonly object _sync = new();
    internal readonly double[] SGrid;

    private static volatile IonizationFsTable _default;
    private static readonly object _defaultSync = new();

    /// <summary>埋め込みリソースから構築される既定インスタンス。</summary>
    public static IonizationFsTable Default
    {
        get
        {
            var t = _default;
            if (t is not null) return t;
            lock (_defaultSync)
            {
                if (_default is null)
                {
                    using var s = Assembly.GetExecutingAssembly().GetManifestResourceStream(ResourceName)
                        ?? throw new InvalidOperationException($"Embedded resource not found: {ResourceName}");
                    _default = new IonizationFsTable(s);
                }
                return _default;
            }
        }
    }

    /// <summary>任意ストリームから構築 (方式比較・破損テスト用)。ヘッダは厳格検査し、unknown codec/method は拒否する。</summary>
    public IonizationFsTable(Stream stream)
    {
        var blob = new byte[stream.Length];
        stream.ReadExactly(blob);
        _blob = blob;
        using var reader = new BinaryReader(new MemoryStream(blob, writable: false));
        if (reader.ReadInt32() != Magic) throw new InvalidDataException("IonizationFsE0.bin: bad magic");
        var formatVersion = reader.ReadInt32();
        //260802Cl: 1 (v1.0.0, L23 1 本) と 2 (v2.0.0, L2/L3 j 分離) を受け入れる。旧: != 1 で拒否
        if (formatVersion is not (1 or 2)) throw new InvalidDataException($"IonizationFsE0.bin: unknown format version {formatVersion}");
        var codec = reader.ReadInt32();
        if (codec != 1) throw new InvalidDataException($"IonizationFsE0.bin: unknown codec {codec}");
        Method = reader.ReadInt32();
        if (Method is not (1 or 2)) throw new InvalidDataException($"IonizationFsE0.bin: unknown method {Method}");
        SCount = reader.ReadInt32();
        SStep = reader.ReadDouble();
        //260805Cl 変更: 81 → 161 (s≤4 → s≤8)。SMaxAngstromInv と整合を取る
        if (SCount != 161 || SStep != 0.05) throw new InvalidDataException("IonizationFsE0.bin: unexpected s grid");
        //if (SCount != 81 || SStep != 0.05) throw ... //260804Cl まで (dataset v1/v2)
        reader.ReadInt32(); // schemaVersion
        DatasetVersion = BoteSalvat.ReadString(reader);
        ModelId = BoteSalvat.ReadString(reader);
        BoteRef = BoteSalvat.ReadString(reader);
        Packer = BoteSalvat.ReadString(reader);
        PayloadSha256 = reader.ReadBytes(32);
        reader.ReadBytes(32); // sourceSha256 (記録用)
        var channelCount = reader.ReadInt32();
        if (channelCount is <= 0 or > 10000) throw new InvalidDataException("IonizationFsE0.bin: bad channel count");
        long payloadLen = 0;
        for (int i = 0; i < channelCount; i++)
        {
            int shellCode = reader.ReadInt32(), z = reader.ReadInt32(), offset = reader.ReadInt32(), length = reader.ReadInt32();
            //260802Cl: 上限を L23(2) → L3(4) へ。旧: shellCode is < ShellCodeK or > ShellCodeL23
            if (shellCode is < ShellCodeK or > ShellCodeL3 || length <= 0 || offset != payloadLen)
                throw new InvalidDataException("IonizationFsE0.bin: bad index entry"); // offset 連続 = 重複/オーバーラップ拒否
            if (!_index.TryAdd((shellCode, z), (offset, length)))
                throw new InvalidDataException($"IonizationFsE0.bin: duplicate channel ({shellCode},{z})");
            payloadLen += length;
        }
        //j 分離の有無は索引の実体で判定する (formatVersion は「読める形式か」だけを表す)
        HasJResolvedL = _index.Keys.Any(k => k.ShellCode == ShellCodeL2);
        _payloadStart = (int)reader.BaseStream.Position;
        if (_payloadStart + payloadLen != blob.Length) throw new InvalidDataException("IonizationFsE0.bin: payload length mismatch");
        SGrid = new double[SCount];
        for (int i = 0; i < SCount; i++) SGrid[i] = i * SStep;
    }

    public bool Contains(int shellCode, int z) => _index.ContainsKey((shellCode, z));

    public IonizationChannelTable GetChannel(int shellCode, int z)
    {
        lock (_sync)
        {
            if (_cache.TryGetValue((shellCode, z), out var hit)) return hit;
            if (!_index.TryGetValue((shellCode, z), out var entry))
                throw new NotSupportedException($"Ionization table has no channel shellCode={shellCode}, Z={z} (K: Z=6–50, L1/L23: Z=20–60)");
            var table = Decode(entry, shellCode, z);
            _cache.Add((shellCode, z), table); // lock 内構築 = ExecutionAndPublication (半初期化を見せない)
            return table;
        }
    }

    /// <summary>全チャネルを展開して canonical payload の SHA-256 を検証 (ハーネス用。runtime 起動時には呼ばない)。</summary>
    public bool VerifyPayloadHash()
    {
        using var sha = IncrementalHash.CreateHash(HashAlgorithmName.SHA256);
        foreach (var kv in _index.OrderBy(e => e.Value.Offset))
            sha.AppendData(DecompressBlob(kv.Value));
        return sha.GetHashAndReset().AsSpan().SequenceEqual(PayloadSha256);
    }

    private byte[] DecompressBlob((int Offset, int Length) entry)
    {
        using var ms = new MemoryStream(_blob, _payloadStart + entry.Offset, entry.Length, writable: false);
        using var br = new BrotliStream(ms, CompressionMode.Decompress);
        using var outMs = new MemoryStream();
        br.CopyTo(outMs);
        return outMs.ToArray();
    }

    private IonizationChannelTable Decode((int Offset, int Length) entry, int expectShell, int expectZ)
    {
        var raw = DecompressBlob(entry);
        using var r = new BinaryReader(new MemoryStream(raw, writable: false));
        int z = r.ReadInt32(), shellCode = r.ReadInt32();
        if (z != expectZ || shellCode != expectShell) throw new InvalidDataException("IonizationFsE0.bin: channel blob/index mismatch");
        var eth = r.ReadDouble();
        var rowCount = r.ReadInt32();
        if (rowCount is < 2 or > 1000) throw new InvalidDataException("IonizationFsE0.bin: bad row count");
        var e0 = new double[rowCount];
        for (int i = 0; i < rowCount; i++) e0[i] = r.ReadDouble();
        var u = new double[rowCount];
        for (int i = 0; i < rowCount; i++) u[i] = r.ReadDouble();
        var tailFlag = r.ReadBytes(rowCount);
        var tailA = new double[rowCount];
        for (int i = 0; i < rowCount; i++) tailA[i] = r.ReadDouble();
        var tailB = new double[rowCount];
        for (int i = 0; i < rowCount; i++) tailB[i] = r.ReadDouble();
        var f = new double[rowCount][];
        if (Method == 1)
        {
            for (int i = 0; i < rowCount; i++)
            {
                var row = new double[SCount];
                for (int j = 0; j < SCount; j++) row[j] = r.ReadSingle();
                f[i] = row;
            }
        }
        else // method 2: int32 量子化 + 行内 s 方向 delta + byte-plane shuffle (out[p*n+i] = raw[i*4+p])
        {
            int n = rowCount * SCount;
            var shuffled = r.ReadBytes(n * 4);
            if (shuffled.Length != n * 4) throw new InvalidDataException("IonizationFsE0.bin: truncated F block");
            var q = new int[n];
            for (int i = 0; i < n; i++)
                q[i] = shuffled[i] | (shuffled[n + i] << 8) | (shuffled[2 * n + i] << 16) | (shuffled[3 * n + i] << 24);
            for (int i = 0; i < rowCount; i++)
            {
                var row = new double[SCount];
                long acc = 0;
                for (int j = 0; j < SCount; j++)
                {
                    acc += q[i * SCount + j];
                    row[j] = acc * 1e-6;
                }
                f[i] = row;
            }
        }
        if (r.BaseStream.Position != raw.Length) throw new InvalidDataException("IonizationFsE0.bin: blob length mismatch");
        // 構造検証: F(0)=1 (method2 は量子化後も 1e6*1e-6)、tail ノード整合 a=F(4)·e^{4b} (量子化誤差ぶんの許容)
        var aTol = Method == 2 ? 1.1e-6 : 6e-8;
        for (int i = 0; i < rowCount; i++)
        {
            if (Math.Abs(f[i][0] - 1.0) > 1e-12) throw new InvalidDataException($"IonizationFsE0.bin: F(0)≠1 (Z={z}, row {i})");
            //260805Cl 変更: 4.0 決め打ち → SMaxAngstromInv (s グリッド上限の延長に追随)
            if (tailFlag[i] != 0 && Math.Abs(tailA[i] * Math.Exp(-SMaxAngstromInv * tailB[i]) - f[i][SCount - 1]) > aTol)
                throw new InvalidDataException($"IonizationFsE0.bin: tail/F(s_max) inconsistent (Z={z}, row {i})");
        }
        return new IonizationChannelTable(this, z, shellCode, eth, e0, u, tailFlag, tailA, tailB, f);
    }
}

/// <summary>260801Cl 追加: 1 チャネル分の展開済みテーブルと E0 補間 (契約 = prod/MANIFEST.md)。</summary>
public sealed class IonizationChannelTable
{
    private readonly IonizationFsTable _owner;
    public readonly int Z;
    public readonly int ShellCode;
    public readonly double EthKeV;
    public readonly double[] E0KeV;   // 厳密昇順
    public readonly double[] U;       // serialized row.u (4桁丸め) = 補間ノット契約値
    private readonly double[] _x;     // ln(u-1)
    private readonly byte[] _tailFlag;
    private readonly double[] _tailA, _tailB;
    private readonly double[][] _f;   // [row][81]
    private readonly List<(int Lo, int Hi)> _tailRuns = [];

    internal IonizationChannelTable(IonizationFsTable owner, int z, int shellCode, double eth,
        double[] e0, double[] u, byte[] tailFlag, double[] tailA, double[] tailB, double[][] f)
    {
        _owner = owner; Z = z; ShellCode = shellCode; EthKeV = eth;
        E0KeV = e0; U = u; _tailFlag = tailFlag; _tailA = tailA; _tailB = tailB; _f = f;
        _x = new double[u.Length];
        for (int i = 0; i < u.Length; i++) _x[i] = Math.Log(u[i] - 1.0);
        for (int i = 0; i < tailFlag.Length;)
        {
            if (tailFlag[i] != 0)
            {
                int k = i;
                while (k + 1 < tailFlag.Length && tailFlag[k + 1] != 0) k++;
                _tailRuns.Add((i, k));
                i = k + 1;
            }
            else i++;
        }
    }

    public int RowCount => E0KeV.Length;

    /// <summary>全 (row, s節点) の F 総和 (golden との構造照合用の診断値)。</summary>
    public double SumF()
    {
        var sum = 0.0;
        foreach (var row in _f)
            foreach (var v in row) sum += v;
        return sum;
    }

    /// <summary>E0 [keV] における 81 節点グリッドを契約どおり補間 (各 s 節点で x=ln(u−1) PCHIP、
    /// 全行正なら lnF・非正含みは符号付き F 直接)。E0 は 30–400 keV 限定 (外挿・clamp 禁止)。</summary>
    public double[] GridAt(double e0KeV)
    {
        if (!(e0KeV >= 30.0 && e0KeV <= 400.0))
            throw new ArgumentOutOfRangeException(nameof(e0KeV), e0KeV, "F(s,E0) table covers E0 = 30–400 keV only (no extrapolation)");
        var xq = Math.Log(e0KeV / EthKeV - 1.0);
        int rows = RowCount, sCount = _owner.SCount;
        var grid = new double[sCount];
        var col = new double[rows];
        for (int j = 0; j < sCount; j++)
        {
            var allPositive = true;
            for (int i = 0; i < rows; i++)
            {
                col[i] = _f[i][j];
                if (col[i] <= 0) allPositive = false;
            }
            if (allPositive)
            {
                for (int i = 0; i < rows; i++) col[i] = Math.Log(col[i]);
                grid[j] = Math.Exp(Pchip.Evaluate(_x, col, Pchip.Derivatives(_x, col), xq));
            }
            else
                grid[j] = Pchip.Evaluate(_x, col, Pchip.Derivatives(_x, col), xq);
        }
        grid[0] = 1.0; // s=0 は厳密 1 (契約)
        return grid;
    }

    /// <summary>s>4 tail の減衰係数 b̂(E0) (連続性アンカー方式、codex 16巡)。取得不能なら false。
    /// tail≠null の連続 E0 区間内のみで PCHIP。E0 を挟む行に null が絡む場合は不可。exact node はその行の b。</summary>
    public bool TryGetTailB(double e0KeV, out double bHat)
    {
        var hit = Array.BinarySearch(E0KeV, e0KeV);
        if (hit >= 0 && _tailFlag[hit] != 0) { bHat = _tailB[hit]; return true; }
        int i = hit >= 0 ? hit : ~hit - 1;
        i = Math.Clamp(i, 0, RowCount - 2);
        foreach (var (lo, hi) in _tailRuns)
            if (lo <= i && i + 1 <= hi && hi > lo)
            {
                var xs = _x[lo..(hi + 1)];
                var bs = _tailB[lo..(hi + 1)];
                bHat = Pchip.Evaluate(xs, bs, Pchip.Derivatives(xs, bs), Math.Log(e0KeV / EthKeV - 1.0));
                return true;
            }
        bHat = double.NaN;
        return false;
    }

    /// <summary>E0 を固定して解決した形状 (run-scoped)。</summary>
    public IonizationTableShape BuildShape(double e0KeV) => new(this, _owner.SGrid, GridAt(e0KeV), e0KeV);
}

/// <summary>260801Cl 追加: E0 解決済みの単一殻形状。s 方向は符号付き F 直接 PCHIP、
/// s>s_max (=8 Å⁻¹, 260805Cl に 4→8 へ延長) は連続性アンカー tail F(s_max)e^{−b̂(s−s_max)}
/// (b̂ 不能なら hard fail)。入力 s の単位は nm⁻¹。</summary>
public sealed class IonizationTableShape : INormalizedIonizationShape
{
    private readonly double[] _sGrid, _grid, _deriv;
    private readonly double _f4, _bHat;   //_f4 = F(s_max)。名前は s_max=4 時代の名残
    private readonly bool _tailAvailable;
    private readonly int _z, _shellCode;
    private bool _usedTail;

    internal IonizationTableShape(IonizationChannelTable table, double[] sGrid, double[] grid, double e0KeV)
    {
        _sGrid = sGrid; _grid = grid;
        _deriv = Pchip.Derivatives(sGrid, grid);
        _f4 = grid[^1];
        _tailAvailable = table.TryGetTailB(e0KeV, out _bHat);
        _z = table.Z; _shellCode = table.ShellCode;
    }

    /// <summary>s>4 Å⁻¹ の tail 外挿を使ったか (診断)。</summary>
    public bool UsedTailExtrapolation => _usedTail;

    public void Evaluate(ReadOnlySpan<double> sPerNm, Span<double> values)
    {
        for (int k = 0; k < sPerNm.Length; k++)
        {
            var sA = sPerNm[k] * 0.1; // nm⁻¹ → Å⁻¹
            if (sA == 0.0)
                values[k] = 1.0;
            else if (sA <= IonizationFsTable.SMaxAngstromInv)
                values[k] = Pchip.Evaluate(_sGrid, _grid, _deriv, sA);
            else if (_tailAvailable)
            {
                //260805Cl 変更: 4.0 決め打ち → SMaxAngstromInv (s グリッド上限の延長に追随)
                values[k] = _f4 * Math.Exp(-_bHat * (sA - IonizationFsTable.SMaxAngstromInv));
                _usedTail = true;
            }
            else
                throw new InvalidOperationException(
                    $"F(s) tail unavailable for s={sA:f3} Å⁻¹ > 4 (Z={_z}, shellCode={_shellCode}: null-tail rows bracket this E0). Reduce gMax or refuse the channel.");
        }
    }
}

/// <summary>260801Cl 追加: LTotal 合成形状 F_L = [σ_L1·F_L1 + (σ_L2+σ_L3)·F_L23]/Σσ (実行時 Bote 重み、MANIFEST 契約)。</summary>
public sealed class IonizationLTotalShape : INormalizedIonizationShape
{
    //260802Cl 変更: 副殻 2 本 (L1 + L23) 固定だったのを可変本数へ。v2 dataset では
    //L1 + L2 + L3 の 3 本になる。σ=0 の副殻は shape=null + 重み 0 で渡ってくる契約。
    //旧: private readonly IonizationTableShape _l1, _l23;  private readonly double _w1, _w23;
    private readonly IonizationTableShape[] _shapes;
    private readonly double[] _weights;   // σ 重み (総和 1 に正規化済み)

    internal IonizationLTotalShape(IonizationTableShape[] shapes, double[] sigmas)
    {
        var total = sigmas.Sum();
        _shapes = shapes;
        _weights = [.. sigmas.Select(s => s / total)];
    }

    public bool UsedTailExtrapolation => _shapes.Any(s => s?.UsedTailExtrapolation ?? false);

    public void Evaluate(ReadOnlySpan<double> sPerNm, Span<double> values)
    {
        Span<double> tmp = sPerNm.Length <= 256 ? stackalloc double[sPerNm.Length] : new double[sPerNm.Length];
        values.Clear();
        for (int i = 0; i < _shapes.Length; i++)
        {
            if (!(_weights[i] > 0)) continue;
            _shapes[i].Evaluate(sPerNm, tmp);
            for (int k = 0; k < values.Length; k++) values[k] += _weights[i] * tmp[k];
        }
    }
}

#endregion

#region プロバイダ (260801Cl 追加)

/// <summary>260801Cl 追加: チャネル指定 → 解決済み IonizationData。run 開始時に 1 回だけ呼び、実行中は immutable を使う (設計書 §5.1)。</summary>
public static class IonizationDataProvider
{
    /// <summary>サポート対象の E0 範囲 [keV] (F テーブルの収録域。外挿・clamp はしない)。</summary>
    public const double MinE0KeV = 30.0, MaxE0KeV = 400.0;

    //260801Cl 変更 (codex 21巡): Resolve と Inspect が E0 範囲・対応殻・edge 規則・σ 規則・provenance を
    //それぞれ独立に書いていた (同値性はハーネスのテストだけが担保。tools はリモート無しのローカルリポで CI にも乗らない)。
    //判定を Describe() 1 か所へ集約し、Resolve = Describe + 例外化 + shape 構築、Inspect = Describe そのもの、とする。
    //これで「Available ⇔ Resolve 成功」がテストではなく構造で保証される。
    //260802Cl 追加: 殻 → その殻を構成するテーブル shellCode 列。空 = その dataset では扱えない殻。
    //LTotal だけが複数本になり、v2 (j 分離) では L1+L2+L3、v1 では L1+L23 を束ねる。
    //ここが「dataset の版差を吸収する唯一の場所」で、Describe / Resolve の分岐はこの 1 本に集約する。
    private static int[] ShellCodesOf(IonizationShell shell, IonizationFsTable table) => shell switch
    {
        IonizationShell.K => [IonizationFsTable.ShellCodeK],
        IonizationShell.L1 => [IonizationFsTable.ShellCodeL1],
        IonizationShell.L2 when table.HasJResolvedL => [IonizationFsTable.ShellCodeL2],
        IonizationShell.L3 when table.HasJResolvedL => [IonizationFsTable.ShellCodeL3],
        IonizationShell.LTotal when table.HasJResolvedL =>
            [IonizationFsTable.ShellCodeL1, IonizationFsTable.ShellCodeL2, IonizationFsTable.ShellCodeL3],
        IonizationShell.LTotal => [IonizationFsTable.ShellCodeL1, IonizationFsTable.ShellCodeL23],
        _ => [],
    };

    //260802Cl 追加: shellCode → Bote の副殻番号 (1=K, 2=L1, 3=L2, 4=L3)。
    //L23 (v1 の 2p 平均) だけは L2+L3 の合算なので σ を 2 本足す必要があり、ここでは扱えない。
    private static double SigmaOf(int shellCode, int z, double eV) => shellCode switch
    {
        IonizationFsTable.ShellCodeK => BoteSalvat.SigmaNm2(z, 1, eV),
        IonizationFsTable.ShellCodeL1 => BoteSalvat.SigmaNm2(z, 2, eV),
        IonizationFsTable.ShellCodeL2 => BoteSalvat.SigmaNm2(z, 3, eV),
        IonizationFsTable.ShellCodeL3 => BoteSalvat.SigmaNm2(z, 4, eV),
        _ => BoteSalvat.SigmaNm2(z, 3, eV) + BoteSalvat.SigmaNm2(z, 4, eV),   // ShellCodeL23
    };

    /// <summary>チャネルの状態・edge・過電圧・σ・provenance を算出する共通コア (shape は作らない = Inspect が安く済む)。</summary>
    private static IonizationChannelInfo Describe(IonizationChannelSpec spec, double e0KeV, IonizationFsTable table)
    {
        //対応殻・収録有無の判定は provenance を作る前に済ませる (早期 return パスで捨てる record を作らない)
        //260802Cl: v2 dataset では L1/L2/L3 も単独で解決できる (旧: K と LTotal だけ)。
        var codes = ShellCodesOf(spec.Shell, table);
        if (codes.Length == 0)
            return new IonizationChannelInfo { Channel = spec, Status = IonizationAvailability.UnsupportedShell };
        if (codes.Any(c => !table.Contains(c, spec.Z)))
            return new IonizationChannelInfo { Channel = spec, Status = IonizationAvailability.UnsupportedElement };
        //LTotal は開いている副殻のうち最小の端 (どれか 1 本でも励起できれば信号は出る)
        var edge = codes.Min(c => table.GetChannel(c, spec.Z).EthKeV);
        var partial = new IonizationChannelInfo
        {
            Channel = spec,
            EdgeEnergyKeV = edge,
            Overvoltage = e0KeV / edge,
            CrossSectionSource = new IonizationDataProvenance("Bote-Salvat-2008", "xion.f/ADNDT95", table.BoteRef),
            ShapeSource = new IonizationDataProvenance(table.ModelId, table.DatasetVersion, "self-generated DHFS tables (tools/IonizationGen prod)")
        };
        if (!(e0KeV >= MinE0KeV && e0KeV <= MaxE0KeV))
            return partial with { Status = IonizationAvailability.E0OutOfRange };

        // σ は各サブシェル自身の edge で計算 (MANIFEST 契約)。閉じている subshell は 0 で自然に落ちる
        //260802Cl: 殻ごとの場合分けをやめ、構成 shellCode の σ を足す形に統一した
        //(旧: K なら subshell 1、それ以外は 2+3+4 決め打ち)。
        var eV = e0KeV * 1e3;
        var sigma = codes.Sum(c => SigmaOf(c, spec.Z, eV));
        //現行 dataset (K: Z=6-50 / L: Z=20-86、E0≥30 keV) では全収録チャネルが励起可能なため BelowEdge は実データでは到達しない
        //(全収録 edge < 30 keV)。synthetic table でのみテスト可能 (codex 20巡)
        return sigma <= 0
            ? partial with { Status = IonizationAvailability.BelowEdge }
            : partial with { Status = IonizationAvailability.Available, SigmaNm2 = sigma };
    }

    /// <summary>解決。E0 範囲外 (30–400 keV 以外) は ArgumentOutOfRangeException、
    /// 未収録 Z/殻・below-edge は NotSupportedException。
    /// 260802Cl: v2 dataset では K / LTotal に加えて L1 / L2 / L3 も単独で解決できる
    /// (v1 dataset を読んでいるときは L2 / L3 が UnsupportedShell になる)。</summary>
    public static IonizationData Resolve(IonizationChannelSpec spec, double e0KeV, IonizationFsTable table = null)
    {
        ArgumentNullException.ThrowIfNull(spec);
        table ??= IonizationFsTable.Default;
        var info = Describe(spec, e0KeV, table);
        switch (info.Status)
        {
            case IonizationAvailability.E0OutOfRange:
                throw new ArgumentOutOfRangeException(nameof(e0KeV), e0KeV, $"STEM-EDX supports E0 = {MinE0KeV}–{MaxE0KeV} keV only (F table range, no extrapolation)");
            case IonizationAvailability.UnsupportedShell:
                throw new NotSupportedException($"IonizationShell.{spec.Shell} is not available in this dataset (j-resolved L: {table.HasJResolvedL})");
            case IonizationAvailability.UnsupportedElement:
                throw new NotSupportedException($"Ionization table has no channel for Z={spec.Z} {spec.Shell} (K: Z=6–50, L: Z=20–86)");
            case IonizationAvailability.BelowEdge:
                throw new NotSupportedException($"Z={spec.Z} {spec.Shell}: below edge at E0={e0KeV} keV (σ=0)");
        }
        //ここから先は Available 確定。shape は σ>0 の成分だけ構築する (σ=0 成分は null + 重み 0 で合成。Evaluate は w>0 の成分しか触らない契約)
        //260802Cl: 単一副殻も複数副殻も同じ経路で組む (旧: K を特別扱いし、L は L1+L23 決め打ち)。
        var codes = ShellCodesOf(spec.Shell, table);
        var eV = e0KeV * 1e3;
        var sigmas = codes.Select(c => SigmaOf(c, spec.Z, eV)).ToArray();
        if (codes.Length == 1)
        {
            var ch = table.GetChannel(codes[0], spec.Z);
            return new IonizationData(spec, info.EdgeEnergyKeV, info.SigmaNm2, ch.BuildShape(e0KeV), info.CrossSectionSource, info.ShapeSource);
        }
        var shapes = codes.Select((c, i) => sigmas[i] > 0 ? table.GetChannel(c, spec.Z).BuildShape(e0KeV) : null).ToArray();
        return new IonizationData(spec, info.EdgeEnergyKeV, info.SigmaNm2,
            new IonizationLTotalShape(shapes, sigmas), info.CrossSectionSource, info.ShapeSource);
    }

    /// <summary>260801Cl 追加: GUI 向け照会 (設計書 §5.9-3)。throw せず状態 enum を返す。
    /// Resolve と同じ <see cref="Describe"/> を使うので <c>Status==Available ⇔ 同じ引数で Resolve が成功する</c> は構造的に保証される。
    /// 判定順: UnsupportedShell → UnsupportedElement → E0OutOfRange (edge/U は返す) → BelowEdge (edge/U は返す) → Available。</summary>
    public static IonizationChannelInfo Inspect(IonizationChannelSpec spec, double e0KeV, IonizationFsTable table = null)
    {
        ArgumentNullException.ThrowIfNull(spec);
        return Describe(spec, e0KeV, table ?? IonizationFsTable.Default);
    }

    /// <summary>260801Cl 追加: 結晶の構成元素から STEM-EDX 候補チャネルを列挙する (設計書 §5.9-3)。
    /// 収録外の元素・殻は返さない。below-edge / E0 範囲外は「理由付きで選べない候補」として返す
    /// (GUI 側にデータ収録範囲 (K: Z=6–50 等) をハードコードさせないための入口)。</summary>
    public static IonizationChannelInfo[] EnumerateChannels(Crystal crystal, double e0KeV, IonizationFsTable table = null)
    {
        if (crystal?.Atoms is null || crystal.Atoms.Length == 0) return [];
        table ??= IonizationFsTable.Default;
        var list = new List<IonizationChannelInfo>();
        foreach (var z in crystal.Atoms.Select(a => a.AtomicNumber).Distinct().OrderBy(z => z))
            //260802Cl: 列挙は K と LTotal のまま据え置く。v2 で L1/L2/L3 も解決可能になったが、
            //EDX は「列挙された全チャネルを計算する」仕様 (9daab2f4) なので、ここに副殻を足すと
            //L 元素の計算量が 3 倍になるうえ、EDS 検出器が見るのは Lα/Lβ という線であって
            //副殻ごとの空孔マップではない。副殻を分けて使うのは蛍光収率・線分岐を入れる発光層の仕事。
            foreach (var shell in new[] { IonizationShell.K, IonizationShell.LTotal })
            {
                var info = Describe(new IonizationChannelSpec(z, shell), e0KeV, table);
                if (info.Status is not (IonizationAvailability.UnsupportedElement or IonizationAvailability.UnsupportedShell))
                    list.Add(info);
            }
        return [.. list];
    }
}

#endregion
