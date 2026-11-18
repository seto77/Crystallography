// 261118Cl 新規作成 (codex2 の単位 6・7 の第 1 段): 動力学計算の数値の状態と生値の集計。利用者に見える値は変えず、記録だけを増やす。
//   設計は .project-guidance/reports/codex2_261118/wf_native_failure_261118.md の statsDesign と第 6 回 (round6_answer.md §3)。
//   集計は worker ごとの NumericStatsAccumulator (ロックなし) に貯め、Parallel.For が戻った後に NumericRunStatsCollector.Merge で合わせる。
//   最悪の条件は順序に依らない全順序で選ぶので、スレッドへの割り振りが変わっても同じ記録になる。
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;

namespace Crystallography;

/// <summary>261118Cl 追加: 集計する物理量 (源)</summary>
public enum NumericQuantity { EbsdLocal, EbsdTds, EbsdReinjection, EbsdTotal }

/// <summary>261118Cl 追加: 状態を数える段</summary>
public enum NumericStage { Evd, Boundary, Integrator }

/// <summary>261118Cl 追加: 最悪の値が出た場所 (電圧の番号・方向の番号・深さの番号) と値・尺度</summary>
public readonly record struct NumericWorstCase(int VoltageIndex, int Direction, int DepthIndex, double Value, double Scale);

/// <summary>261118Cl 追加: 1 つの源の集計。Evaluated は計算した生値の数 (計算済みの印の 0 などは Excluded)。
/// Negative は有限で x &lt; 0 (−0 は数えない)、NonFinite は NaN・±∞。NegativeDecades[k] は |x|/尺度 が 10^−k 以上 10^−(k−1) 未満の負値の数
/// (尺度 = 同じ方向の同じ源の正の最大。k = 0 は 1 以上、19 は 10^−18 未満)</summary>
public sealed class NumericQuantitySummary
{
    public long Evaluated, Negative, NonFinite, Excluded;
    public double Min = double.PositiveInfinity, MaxPositive = 0;
    public readonly long[] NegativeDecades = new long[20];
    public NumericWorstCase? WorstNegative, FirstNonFinite;
    /// <summary>最悪の負値の、同じ方向の尺度に対する比 (負値が無ければ 0、尺度が 0 なら NaN)</summary>
    public double WorstNegativeRelative => WorstNegative is { } w ? (w.Scale > 0 ? -w.Value / w.Scale : double.NaN) : 0;

    internal void Merge(NumericQuantitySummary o)
    {
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
    internal static bool Before(in NumericWorstCase x, in NumericWorstCase y) =>
        x.VoltageIndex != y.VoltageIndex ? x.VoltageIndex < y.VoltageIndex : x.Direction != y.Direction ? x.Direction < y.Direction : x.DepthIndex <= y.DepthIndex;
}

/// <summary>261118Cl 追加: worker ごとの集計 (ロックなし。1 つのスレッドだけが書く)</summary>
internal sealed class NumericStatsAccumulator
{
    internal readonly NumericQuantitySummary[] Quantities = Enumerable.Range(0, Enum.GetValues<NumericQuantity>().Length).Select(_ => new NumericQuantitySummary()).ToArray();
    internal readonly Dictionary<(NumericStage, NativeWrapper.NativeSolveStatus), long> Statuses = [];
    internal NumericWorstCase? FirstFailure; internal (NumericStage, NativeWrapper.NativeSolveStatus)? FirstFailureKind;

    /// <summary>生値の列 (深さごと) を 1 つの方向として数える。尺度は同じ列の正の最大</summary>
    internal void Observe(NumericQuantity q, ReadOnlySpan<double> raw, int voltageIndex, int direction)
    {
        var s = Quantities[(int)q];
        double scale = 0;
        foreach (var x in raw) if (double.IsFinite(x) && x > scale) scale = x;
        for (int t = 0; t < raw.Length; t++)
        {
            double x = raw[t];
            s.Evaluated++;
            if (!double.IsFinite(x)) { s.NonFinite++; s.FirstNonFinite = NumericQuantitySummary.PickFirst(s.FirstNonFinite, new(voltageIndex, direction, t, x, scale)); continue; }
            if (x < s.Min) s.Min = x;
            if (x > s.MaxPositive) s.MaxPositive = x;
            if (x < 0)
            {
                s.Negative++;
                double rel = scale > 0 ? -x / scale : double.PositiveInfinity;
                int k = rel >= 1 ? 0 : rel > 0 ? Math.Min(19, (int)Math.Ceiling(-Math.Log10(rel))) : 19;
                s.NegativeDecades[k]++;
                s.WorstNegative = NumericQuantitySummary.PickWorstNegative(s.WorstNegative, new(voltageIndex, direction, t, x, scale));
            }
        }
    }

    /// <summary>計算しない (計算済みの印の 0 など) 値の数</summary>
    internal void Exclude(NumericQuantity q, int count) => Quantities[(int)q].Excluded += count;

    internal void Status(NumericStage stage, NativeWrapper.NativeSolveStatus status, int voltageIndex, int direction)
    {
        var key = (stage, status);
        Statuses[key] = Statuses.TryGetValue(key, out var c) ? c + 1 : 1;
        if (status is not (NativeWrapper.NativeSolveStatus.Ok or NativeWrapper.NativeSolveStatus.Unchecked))
        {
            var here = new NumericWorstCase(voltageIndex, direction, -1, double.NaN, double.NaN);
            if (FirstFailure is null || NumericQuantitySummary.Before(here, FirstFailure.Value)) { FirstFailure = here; FirstFailureKind = key; }
        }
    }
}

/// <summary>261118Cl 追加: 1 回の計算 (run) の集計 (不変)。<see cref="Summary"/> は記録用の 1 行</summary>
public sealed class NumericRunStats
{
    public string Kind { get; init; } = "";
    /// <summary>積分器: native-raw (生値を見た) / native-legacy (旧い DLL、丸め済み) / managed</summary>
    public string IntegratorBackend { get; init; } = "";
    /// <summary>固有値: native-checked / native-legacy-unchecked / managed-MathNet / research-hook</summary>
    public string EvdBackend { get; init; } = "";
    public string NativeLibrary { get; init; } = "";
    /// <summary>生値を見たか (旧い DLL の丸め済みの値では負値は見えない)</summary>
    public bool RawObserved { get; init; }
    public IReadOnlyDictionary<NumericQuantity, NumericQuantitySummary> Quantities { get; init; }
    public IReadOnlyDictionary<(NumericStage Stage, NativeWrapper.NativeSolveStatus Status), long> Statuses { get; init; }
    public NumericWorstCase? FirstFailure { get; init; }
    public (NumericStage Stage, NativeWrapper.NativeSolveStatus Status)? FirstFailureKind { get; init; }

    public string Summary()
    {
        var sb = new StringBuilder($"{Kind}: integrator {IntegratorBackend}, EVD {EvdBackend}, native {NativeLibrary}, raw {(RawObserved ? "observed" : "not observed")}");
        foreach (var (q, s) in Quantities.Where(kv => kv.Value.Evaluated + kv.Value.Excluded > 0))
            sb.Append($"; {q}: n {s.Evaluated} (excluded {s.Excluded}), negative {s.Negative}, non-finite {s.NonFinite}, min {s.Min:g4}, max {s.MaxPositive:g4}, worst negative/scale {s.WorstNegativeRelative:e2}");
        var bad = Statuses.Where(kv => kv.Key.Status is not (NativeWrapper.NativeSolveStatus.Ok or NativeWrapper.NativeSolveStatus.Unchecked)).ToArray();
        sb.Append(bad.Length == 0 ? "; status all ok" : "; status " + string.Join(", ", bad.Select(kv => $"{kv.Key.Stage}/{kv.Key.Status} {kv.Value}")));
        return sb.ToString();
    }
}

/// <summary>261118Cl 追加: run ごとの集め役 (スレッドに安全)。worker の集計を Merge し、Build で <see cref="NumericRunStats"/> にする</summary>
internal sealed class NumericRunStatsCollector(string kind, string integratorBackend, string evdBackend, bool rawObserved)
{
    private readonly object gate = new();
    private readonly NumericStatsAccumulator total = new();

    internal void Merge(NumericStatsAccumulator acc)
    {
        if (acc is null) return;
        lock (gate)
        {
            for (int q = 0; q < total.Quantities.Length; q++) total.Quantities[q].Merge(acc.Quantities[q]);
            foreach (var (k, v) in acc.Statuses) total.Statuses[k] = total.Statuses.TryGetValue(k, out var c) ? c + v : v;
            if (acc.FirstFailure is { } f && (total.FirstFailure is null || NumericQuantitySummary.Before(f, total.FirstFailure.Value)))
            { total.FirstFailure = f; total.FirstFailureKind = acc.FirstFailureKind; }
        }
    }

    internal NumericRunStats Build()
    {
        lock (gate)
            return new NumericRunStats
            {
                Kind = kind, IntegratorBackend = integratorBackend, EvdBackend = evdBackend, NativeLibrary = NativeWrapper.LoadedNativeLibrary ?? "none", RawObserved = rawObserved,
                Quantities = Enum.GetValues<NumericQuantity>().ToDictionary(q => q, q => total.Quantities[(int)q]),
                Statuses = new Dictionary<(NumericStage, NativeWrapper.NativeSolveStatus), long>(total.Statuses),
                FirstFailure = total.FirstFailure, FirstFailureKind = total.FirstFailureKind,
            };
    }
}
