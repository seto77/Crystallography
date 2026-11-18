// 261118Cl 新規作成 (codex2 の単位 6・7 の第 1 段): 動力学計算の数値の状態と生値の集計。正常な値は変えず、記録を増やし、失敗は明示する。
//   設計は .project-guidance/reports/codex2_261118/wf_native_failure_261118.md の statsDesign と第 6 回 (round6_answer.md §3)・第 7 回 (round7_answer.md §1)。
//   集計は worker ごとの NumericStatsAccumulator (ロックなし) に貯め、Parallel.For が戻った後に NumericRunStatsCollector.Merge で合わせる。
//   最悪の条件・最初の失敗は順序に依らない全順序で選ぶので、スレッドへの割り振りが変わっても同じ記録になる。
//   261118Cl 改訂 (第 7 回): 段の失敗・観測した値の非有限は呼び出し側が NumericFailureException で止める (次の段へ渡さない)。
//   要約は「未検査」「評価なし」「観測しなかった源」を成功と区別し、run の状態 (完了・取消・失敗) を持つ。公開する集計は変更できない写し。
using System;
using System.Collections.Frozen;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Crystallography;

/// <summary>261118Cl 追加: 集計する物理量 (源)。KikuchiI・KikuchiIRef は菊池のプロファイルの I と参照強度 I_ref (261118Cl 追加)</summary>
public enum NumericQuantity { EbsdLocal, EbsdTds, EbsdReinjection, EbsdTotal, KikuchiI, KikuchiIRef }

/// <summary>261118Cl 追加: 段。並びは同じ方向の中の先後 (EVD → 境界 → 積分 → 値の検査)。ValueCheck は観測した値の有限性 (solver の状態ではない)</summary>
public enum NumericStage { Evd, Boundary, Integrator, ValueCheck }

/// <summary>261118Cl 追加: 源の値をどの段階で観測したか (Raw = 丸めの前の生値、AfterClamp = 負・NaN を 0 にする丸めの後)</summary>
[Flags]
public enum NumericObservation { None = 0, Raw = 1, AfterClamp = 2 }

/// <summary>261118Cl 追加: run の終わり方</summary>
public enum NumericRunState { Completed, Cancelled, Failed }

/// <summary>261118Cl 追加: 最悪の値が出た場所 (電圧の番号・方向の番号・深さの番号) と値・尺度</summary>
public readonly record struct NumericWorstCase(int VoltageIndex, int Direction, int DepthIndex, double Value, double Scale);

/// <summary>261118Cl 追加 (codex2 第 7 回): 失敗の場所と種類。値の検査 (<see cref="NumericStage.ValueCheck"/>) なら源と深さの番号も持つ (それ以外は Quantity = null・DepthIndex = −1)</summary>
public readonly record struct NumericFailure(int VoltageIndex, int Direction, NumericStage Stage, NativeWrapper.NativeSolveStatus Status, NumericQuantity? Quantity = null, int DepthIndex = -1)
{
    /// <summary>順序に依らない全順序: 電圧 → 方向 → 段 → 源 → 深さ (同じ方向の EVD の失敗が、後の段の失敗に上書きされない)</summary>
    internal static int Compare(in NumericFailure a, in NumericFailure b)
    {
        int c = a.VoltageIndex.CompareTo(b.VoltageIndex); if (c != 0) return c;
        c = a.Direction.CompareTo(b.Direction); if (c != 0) return c;
        c = a.Stage.CompareTo(b.Stage); if (c != 0) return c;
        c = ((int?)a.Quantity ?? -1).CompareTo((int?)b.Quantity ?? -1); if (c != 0) return c;
        c = a.DepthIndex.CompareTo(b.DepthIndex); if (c != 0) return c;
        return a.Status.CompareTo(b.Status);
    }

    public override string ToString() =>
        $"voltage #{VoltageIndex}, direction #{Direction}, {Stage} {Status}" + (Quantity is { } q ? $" ({q}, depth #{DepthIndex})" : "");
}

/// <summary>261118Cl 追加 (codex2 第 7 回): 動力学計算の数値の失敗 (solver の非成功の状態・観測した値の NaN・±∞)。
/// BackgroundWorker の完了の通知の Error にそのまま入る</summary>
public sealed class NumericFailureException(NumericFailure failure, string message) : Exception(message)
{
    public NumericFailure Failure { get; } = failure;
}

/// <summary>261118Cl 追加: 1 つの源の集計 (変更できない写し)。Evaluated は計算した値の数 (計算済みの印の 0 などは Excluded)。
/// Negative は有限で x &lt; 0 (−0 は数えない)、NonFinite は NaN・±∞。NegativeDecades[k] は |x|/尺度 が 10^−k 以上 10^−(k−1) 未満の負値の数
/// (尺度 = 同じ方向の同じ源の正の最大。k = 0 は 1 以上、19 は 10^−18 未満)</summary>
public sealed class NumericQuantityStats
{
    public NumericObservation Observation { get; }
    public long Evaluated { get; }
    public long Negative { get; }
    public long NonFinite { get; }
    public long Excluded { get; }
    /// <summary>有限の値の最小 (有限の値が 1 つも無ければ null)</summary>
    public double? Min { get; }
    public double MaxPositive { get; }
    public IReadOnlyList<long> NegativeDecades { get; }
    public NumericWorstCase? WorstNegative { get; }
    public NumericWorstCase? FirstNonFinite { get; }
    /// <summary>最悪の負値の、同じ方向の尺度に対する比 (負値が無ければ 0、尺度が 0 なら NaN)</summary>
    public double WorstNegativeRelative => WorstNegative is { } w ? (w.Scale > 0 ? -w.Value / w.Scale : double.NaN) : 0;

    internal NumericQuantityStats(NumericQuantityAccumulator a)
    {
        Observation = a.Observation; Evaluated = a.Evaluated; Negative = a.Negative; NonFinite = a.NonFinite; Excluded = a.Excluded;
        Min = a.Evaluated > a.NonFinite ? a.Min : null; MaxPositive = a.MaxPositive;
        NegativeDecades = Array.AsReadOnly((long[])a.NegativeDecades.Clone()); // 写し (集め役の配列を公開しない)
        WorstNegative = a.WorstNegative; FirstNonFinite = a.FirstNonFinite;
    }
}

/// <summary>261118Cl 追加: 1 つの源の集計の途中 (1 つのスレッドだけが書く)</summary>
internal sealed class NumericQuantityAccumulator
{
    internal NumericObservation Observation;
    internal long Evaluated, Negative, NonFinite, Excluded;
    internal double Min = double.PositiveInfinity, MaxPositive = 0;
    internal readonly long[] NegativeDecades = new long[20];
    internal NumericWorstCase? WorstNegative, FirstNonFinite;

    internal void Merge(NumericQuantityAccumulator o)
    {
        Observation |= o.Observation;
        Evaluated += o.Evaluated; Negative += o.Negative; NonFinite += o.NonFinite; Excluded += o.Excluded;
        if (o.Min < Min) Min = o.Min;
        if (o.MaxPositive > MaxPositive) MaxPositive = o.MaxPositive;
        for (int k = 0; k < NegativeDecades.Length; k++) NegativeDecades[k] += o.NegativeDecades[k];
        WorstNegative = PickWorstNegative(WorstNegative, o.WorstNegative);
        FirstNonFinite = PickFirst(FirstNonFinite, o.FirstNonFinite);
    }

    /// <summary>|値|/尺度 の大きい方。同じなら (電圧, 方向, 深さ) の辞書順の小さい方 (順序に依らない)</summary>
    internal static NumericWorstCase? PickWorstNegative(NumericWorstCase? a, NumericWorstCase? b)
    {
        if (a is null) return b; if (b is null) return a;
        double ra = a.Value.Scale > 0 ? -a.Value.Value / a.Value.Scale : double.PositiveInfinity, rb = b.Value.Scale > 0 ? -b.Value.Value / b.Value.Scale : double.PositiveInfinity;
        if (ra != rb) return ra > rb ? a : b;
        return Before(a.Value, b.Value) ? a : b;
    }
    internal static NumericWorstCase? PickFirst(NumericWorstCase? a, NumericWorstCase? b) => a is null ? b : b is null ? a : Before(a.Value, b.Value) ? a : b;
    /// <summary>(電圧, 方向, 深さ) の辞書順で x が y より前か同じ (同じ場所は同じ値なので、どちらを採っても同じ)</summary>
    internal static bool Before(in NumericWorstCase x, in NumericWorstCase y) =>
        x.VoltageIndex != y.VoltageIndex ? x.VoltageIndex < y.VoltageIndex : x.Direction != y.Direction ? x.Direction < y.Direction : x.DepthIndex <= y.DepthIndex;
}

/// <summary>261118Cl 追加: worker ごとの集計 (ロックなし。1 つのスレッドだけが書く)</summary>
internal sealed class NumericStatsAccumulator
{
    internal readonly NumericQuantityAccumulator[] Quantities = Enumerable.Range(0, Enum.GetValues<NumericQuantity>().Length).Select(_ => new NumericQuantityAccumulator()).ToArray();
    internal readonly Dictionary<(NumericStage, NativeWrapper.NativeSolveStatus), long> Statuses = [];
    /// <summary>解きに入った (電圧, 代表方向) の数と、そのうちビームが無くて解かずに戻った数</summary>
    internal long Directions, NoBeamDirections;
    internal NumericFailure? FirstFailure;

    /// <summary>値の列 (深さごと) を 1 つの方向として数える。尺度は同じ列の正の最大 (scales を渡せば深さごとにその値。真値が 0 になりうる再注入は Σσ·t)。
    /// NaN・±∞ が有れば列の最初の場所を失敗として記録して返す (呼び出し側はその値を使わずに止める)。全部有限なら null</summary>
    internal NumericFailure? Observe(NumericQuantity q, NumericObservation kind, ReadOnlySpan<double> values, int voltageIndex, int direction, ReadOnlySpan<double> scales = default)
    {
        var s = Quantities[(int)q];
        s.Observation |= kind;
        double columnScale = 0;
        foreach (var x in values) if (double.IsFinite(x) && x > columnScale) columnScale = x;
        NumericFailure? failure = null;
        for (int t = 0; t < values.Length; t++)
        {
            double x = values[t], scale = scales.IsEmpty ? columnScale : scales[t];
            s.Evaluated++;
            if (!double.IsFinite(x))
            {
                s.NonFinite++; s.FirstNonFinite = NumericQuantityAccumulator.PickFirst(s.FirstNonFinite, new(voltageIndex, direction, t, x, scale));
                failure ??= new NumericFailure(voltageIndex, direction, NumericStage.ValueCheck, NativeWrapper.NativeSolveStatus.NonFiniteOutput, q, t);
                continue;
            }
            if (x < s.Min) s.Min = x;
            if (x > s.MaxPositive) s.MaxPositive = x;
            if (x < 0)
            {
                s.Negative++;
                double rel = scale > 0 ? -x / scale : double.PositiveInfinity;
                int k = rel >= 1 ? 0 : rel > 0 ? Math.Min(19, (int)Math.Ceiling(-Math.Log10(rel))) : 19;
                s.NegativeDecades[k]++;
                s.WorstNegative = NumericQuantityAccumulator.PickWorstNegative(s.WorstNegative, new(voltageIndex, direction, t, x, scale));
            }
        }
        if (failure is { } f) RecordFailure(f);
        return failure;
    }

    /// <summary>計算しない (計算済みの印の 0 など) 値の数</summary>
    internal void Exclude(NumericQuantity q, int count) => Quantities[(int)q].Excluded += count;

    /// <summary>段の状態を数える。Ok と Unchecked なら null、それ以外は失敗として記録して返す (呼び出し側はその方向の次の段へ進まない)</summary>
    internal NumericFailure? Status(NumericStage stage, NativeWrapper.NativeSolveStatus status, int voltageIndex, int direction)
    {
        var key = (stage, status);
        Statuses[key] = Statuses.TryGetValue(key, out var c) ? c + 1 : 1;
        if (status is NativeWrapper.NativeSolveStatus.Ok or NativeWrapper.NativeSolveStatus.Unchecked) return null;
        var f = new NumericFailure(voltageIndex, direction, stage, status);
        RecordFailure(f);
        return f;
    }

    internal void RecordFailure(in NumericFailure f)
    {
        if (FirstFailure is null || NumericFailure.Compare(f, FirstFailure.Value) < 0) FirstFailure = f;
    }

    internal void Merge(NumericStatsAccumulator o)
    {
        for (int q = 0; q < Quantities.Length; q++) Quantities[q].Merge(o.Quantities[q]);
        foreach (var (k, v) in o.Statuses) Statuses[k] = Statuses.TryGetValue(k, out var c) ? c + v : v;
        Directions += o.Directions; NoBeamDirections += o.NoBeamDirections;
        if (o.FirstFailure is { } f) RecordFailure(f);
    }
}

/// <summary>261118Cl 追加: 1 回の計算 (run) の集計 (変更できない写し)。<see cref="Summary"/> は記録用の 1 行</summary>
public sealed class NumericRunStats
{
    public string Kind { get; }
    public NumericRunState RunState { get; }
    /// <summary>積分器: native-raw (生値を見た) / native-legacy (旧い DLL、丸め済み) / managed / dressed-tail (従来の集計を通らない)</summary>
    public string IntegratorBackend { get; }
    /// <summary>固有値: native-checked / native-legacy-unchecked / managed-MathNet / research-hook / dressed-tail</summary>
    public string EvdBackend { get; }
    public string NativeLibrary { get; }
    /// <summary>局所・TDS の源を生値で見たか (旧い DLL・managed・着衣の裾では見ない)</summary>
    public bool RawObserved { get; }
    /// <summary>解きに入った (電圧, 代表方向) の数と、そのうちビームが無くて解かずに戻った数 (実験の分母の確認用)</summary>
    public long Directions { get; }
    public long NoBeamDirections { get; }
    public IReadOnlyDictionary<NumericQuantity, NumericQuantityStats> Quantities { get; }
    public IReadOnlyDictionary<(NumericStage Stage, NativeWrapper.NativeSolveStatus Status), long> Statuses { get; }
    public NumericFailure? FirstFailure { get; }

    internal NumericRunStats(string kind, NumericRunState state, string integratorBackend, string evdBackend, bool rawObserved, NumericStatsAccumulator total)
    {
        Kind = kind; RunState = state; IntegratorBackend = integratorBackend; EvdBackend = evdBackend; RawObserved = rawObserved;
        NativeLibrary = NativeWrapper.LoadedNativeLibrary ?? "none";
        Directions = total.Directions; NoBeamDirections = total.NoBeamDirections;
        Quantities = Enum.GetValues<NumericQuantity>().ToFrozenDictionary(q => q, q => new NumericQuantityStats(total.Quantities[(int)q]));
        Statuses = total.Statuses.ToFrozenDictionary(kv => kv.Key, kv => kv.Value);
        FirstFailure = total.FirstFailure;
    }

    /// <summary>記録用の 1 行。段は「unchecked (状態を見ていない)」と状態ごとの数、源は「unobserved (見ていない)」・観測の段階 (raw / after clamp)・
    /// 「no evaluated values (評価なし)」を区別する</summary>
    public string Summary()
    {
        var sb = new StringBuilder($"{Kind} {RunState.ToString().ToLowerInvariant()}: directions {Directions} (no beams {NoBeamDirections}); integrator {IntegratorBackend}, EVD {EvdBackend}, native {NativeLibrary}, raw {(RawObserved ? "observed" : "not observed")}; stages");
        foreach (var stage in new[] { NumericStage.Evd, NumericStage.Boundary, NumericStage.Integrator })
        {
            var counts = Statuses.Where(kv => kv.Key.Stage == stage).OrderBy(kv => (int)kv.Key.Status).ToArray();
            sb.Append(stage == NumericStage.Evd ? " " : ", ").Append(stage).Append(' ')
              .Append(counts.Length == 0 ? "unchecked" : string.Join("/", counts.Select(kv => $"{kv.Key.Status} {kv.Value}")));
        }
        foreach (var (q, s) in Quantities.Where(kv => kv.Key.ToString().StartsWith(Kind, StringComparison.OrdinalIgnoreCase)).OrderBy(kv => (int)kv.Key)) // その計算の種類の源だけ (EBSD → Ebsd*、Kikuchi → Kikuchi*)
        {
            sb.Append($"; {q} ");
            if (s.Observation == NumericObservation.None && s.Evaluated == 0 && s.Excluded == 0) { sb.Append("unobserved"); continue; }
            sb.Append(s.Observation switch { NumericObservation.Raw => "raw", NumericObservation.AfterClamp => "after clamp", NumericObservation.None => "not evaluated", _ => "raw+after clamp" });
            sb.Append($": n {s.Evaluated} (excluded {s.Excluded}), ");
            sb.Append(s.Evaluated == 0 ? "no evaluated values"
                : $"negative {s.Negative}, non-finite {s.NonFinite}, min {(s.Min is { } m ? m.ToString("g4") : "n/a")}, max {s.MaxPositive:g4}, worst negative/scale {s.WorstNegativeRelative:e2}");
            if (s.Negative > 0) // 負値の |x|/尺度 の桁の度数 (k: 10^−k 以上 10^−(k−1) 未満)
                sb.Append(" (negative/scale decades " + string.Join(" ", s.NegativeDecades.Select((c, k) => (c, k)).Where(x => x.c > 0).Select(x => $"1e-{x.k}:{x.c}")) + ")");
        }
        sb.Append("; first failure ").Append(FirstFailure is { } f ? f.ToString() : "none");
        return sb.ToString();
    }
}

/// <summary>261118Cl 追加: run ごとの集め役 (スレッドに安全)。worker の集計を Merge し、Build で <see cref="NumericRunStats"/> (写し) にする。
/// backend の欄は構築の開始時に決め、着衣の裾を採ったら書き換える</summary>
internal sealed class NumericRunStatsCollector(string kind)
{
    internal string IntegratorBackend = "", EvdBackend = "";
    internal bool RawObserved;
    private readonly object gate = new();
    private readonly NumericStatsAccumulator total = new();

    internal void Merge(NumericStatsAccumulator acc)
    {
        if (acc is null) return;
        lock (gate) total.Merge(acc);
    }

    internal NumericRunStats Build(NumericRunState state)
    {
        lock (gate) return new NumericRunStats(kind, state, IntegratorBackend, EvdBackend, RawObserved, total);
    }
}
