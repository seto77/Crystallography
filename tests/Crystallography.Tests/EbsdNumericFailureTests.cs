using Xunit;

namespace Crystallography.Tests;

/// <summary>261118Cl 新規 (codex2 第 7 回): EBSD の数値の検査の失敗の伝わり方と、集計 (<see cref="NumericRunStats"/>) の区別の回帰。
/// ・検査つきの EVD が失敗したら、その方向の境界・積分を呼ばず、構築は <see cref="NumericFailureException"/> で失敗する (完了の通知の Error)。
/// ・生値の積分器が失敗したら、その出力を評価しない (評価数 0)。
/// ・集計は「未検査」「観測しなかった源」「評価なし」を成功と区別し、最初の失敗は (電圧, 方向, 段) の順で決まり、公開した集計は後から変わらない。
/// 集計の型の試験は native が無くても走る。EBSD の構築の試験は検査つきの native の入口が要る (無ければ Skip)</summary>
[Collection("Bethe global settings")]
public class EbsdNumericFailureTests
{
    private static Vector3DBase[] Directions(int grid) =>
        [.. MasterPattern.CreateDirectionsSquare(grid, MasterPattern.Hemisphere.PositiveZ), .. MasterPattern.CreateDirectionsSquare(grid, MasterPattern.Hemisphere.NegativeZ)];

    private readonly record struct EbsdOutcome(Exception Error, bool Cancelled, NumericRunStats Stats, BetheMethod.CBED_Disk[][] Disks);

    /// <summary>EBSD (Si、6×6×2 方向、N16、20 kV、深さ 10・20 nm) を走らせて完了を待つ。完了の通知は SynchronizationContext の無いテストではスレッドプールで届く</summary>
    private static EbsdOutcome Run(bool tds, Action<BetheMethod> configure = null)
    {
        var (crystal, _) = BetheEigenReuseTests.TiltedSi();
        var bethe = new BetheMethod(crystal) { EbsdRepresentativeDirectionBlockSize = 1, UseLocalSurfacePerBeamDirection = false };
        configure?.Invoke(bethe);
        using var done = new ManualResetEventSlim();
        EbsdOutcome outcome = default;
        void C(object s, System.ComponentModel.RunWorkerCompletedEventArgs e)
        {
            outcome = new(e.Error, e.Cancelled, bethe.EbsdNumericStats, bethe.Disks);
            done.Set();
        }
        bethe.EBSD_Completed += C;
        try
        {
            bethe.RunEBSDNew(16, [20.0], Matrix3D.IdentityMatrix, [10.0, 20.0], Directions(6), BetheMethod.Solver.Eigen_Eigen, 4, includeTDSBackground: tds);
            Assert.True(done.Wait(TimeSpan.FromMinutes(5)), "EBSD did not finish");
        }
        finally { bethe.EBSD_Completed -= C; }
        return outcome;
    }

    private static long Count(NumericRunStats s, NumericStage stage) => s.Statuses.Where(kv => kv.Key.Stage == stage).Sum(kv => kv.Value);

    [NativeNumericFact]
    public void NormalRun_IsCompleted_AndEveryStageIsCounted()
    {
        foreach (var tds in new[] { false, true })
        {
            var r = Run(tds);
            Assert.Null(r.Error); Assert.False(r.Cancelled);
            var s = r.Stats;
            Assert.Equal(NumericRunState.Completed, s.RunState);
            Assert.Equal("native-raw", s.IntegratorBackend); Assert.Equal("native-checked", s.EvdBackend); Assert.True(s.RawObserved);
            long solved = s.Directions - s.NoBeamDirections;
            Assert.True(solved > 0);
            foreach (var stage in new[] { NumericStage.Evd, NumericStage.Boundary, NumericStage.Integrator })
            {
                Assert.Equal(solved, Count(s, stage));
                Assert.Equal(solved, s.Statuses.GetValueOrDefault((stage, NativeWrapper.NativeSolveStatus.Ok)));
            }
            var source = s.Quantities[tds ? NumericQuantity.EbsdTds : NumericQuantity.EbsdLocal];
            Assert.Equal(NumericObservation.Raw, source.Observation);
            Assert.Equal(solved * 2, source.Evaluated); // 深さ 2 点
            Assert.Equal(0, source.NonFinite);
            Assert.Equal(NumericObservation.AfterClamp, s.Quantities[NumericQuantity.EbsdTotal].Observation);
            Assert.Null(s.FirstFailure);
            Assert.StartsWith("EBSD completed:", s.Summary());
            Assert.Contains("first failure none", s.Summary());
            Assert.NotNull(r.Disks[0]);
        }
    }

    /// <summary>native の試験の旗 (setMaxIterations(0)) で EVD を未収束にする。境界・積分は 1 度も呼ばれず、源は評価されず、マスターの段 (Disks) も作られない</summary>
    [NativeNumericFact]
    public void EvdFailure_StopsBeforeBoundaryAndIntegrator_AndFailsTheRun()
    {
        var r = Run(tds: false, b => b.EbsdEvdCheckFlagsForTest = NativeWrapper.NumericCheckFlags.TestForceNoConvergence);
        var ex = Assert.IsType<NumericFailureException>(r.Error);
        Assert.Equal(NumericStage.Evd, ex.Failure.Stage);
        Assert.Equal(NativeWrapper.NativeSolveStatus.EvdNoConvergence, ex.Failure.Status);
        Assert.Contains("20 kV", ex.Message);
        var s = r.Stats;
        Assert.Equal(NumericRunState.Failed, s.RunState);
        Assert.True(s.Statuses.GetValueOrDefault((NumericStage.Evd, NativeWrapper.NativeSolveStatus.EvdNoConvergence)) > 0);
        // 次の段は前の段が成功した方向でだけ呼ばれる (未収束の旗でも 1×1 などは収束しうるので、0 とは決めない)
        Assert.Equal(s.Statuses.GetValueOrDefault((NumericStage.Evd, NativeWrapper.NativeSolveStatus.Ok)), Count(s, NumericStage.Boundary));
        Assert.Equal(s.Statuses.GetValueOrDefault((NumericStage.Boundary, NativeWrapper.NativeSolveStatus.Ok)), Count(s, NumericStage.Integrator));
        Assert.Equal(s.Statuses.GetValueOrDefault((NumericStage.Integrator, NativeWrapper.NativeSolveStatus.Ok)) * 2, s.Quantities[NumericQuantity.EbsdLocal].Evaluated);
        Assert.Equal(NumericStage.Evd, s.FirstFailure?.Stage);
        Assert.Contains("EBSD failed:", s.Summary());
        Assert.Null(r.Disks[0]);
    }

    /// <summary>生値の積分器の状態を失敗に差し替える。出力 (確保直後の 0) は評価・公開されない</summary>
    [NativeNumericFact]
    public void IntegratorFailure_IsNotEvaluated_AndFailsTheRun()
    {
        foreach (var tds in new[] { false, true })
        {
            var r = Run(tds, b => b.EbsdNumericStatusForTest = (stage, _, status) => stage == NumericStage.Integrator ? NativeWrapper.NativeSolveStatus.OutOfMemory : status);
            var ex = Assert.IsType<NumericFailureException>(r.Error);
            Assert.Equal(NumericStage.Integrator, ex.Failure.Stage);
            Assert.Equal(NativeWrapper.NativeSolveStatus.OutOfMemory, ex.Failure.Status);
            var s = r.Stats;
            Assert.Equal(NumericRunState.Failed, s.RunState);
            Assert.All(s.Quantities.Values, q => Assert.Equal(0, q.Evaluated));
            Assert.Equal(0, s.Statuses.GetValueOrDefault((NumericStage.Integrator, NativeWrapper.NativeSolveStatus.Ok)));
            Assert.Null(r.Disks[0]);
        }
    }

    // ---- 集計の型 (native なしで走る) ----

    [Fact]
    public void EmptyStats_AreUncheckedAndUnobserved_NotOk()
    {
        var s = new NumericRunStatsCollector("EBSD") { IntegratorBackend = "managed", EvdBackend = "managed-MathNet" }.Build(NumericRunState.Completed);
        var text = s.Summary();
        Assert.Contains("Evd unchecked", text); Assert.Contains("Boundary unchecked", text); Assert.Contains("Integrator unchecked", text);
        Assert.Contains("EbsdLocal unobserved", text);
        Assert.DoesNotContain("all ok", text);
        Assert.DoesNotContain("∞", text);
    }

    [Fact]
    public void ExcludedOnly_ReadsAsNoEvaluatedValues()
    {
        var acc = new NumericStatsAccumulator();
        acc.Exclude(NumericQuantity.EbsdLocal, 4);
        var c = new NumericRunStatsCollector("EBSD"); c.Merge(acc);
        var s = c.Build(NumericRunState.Completed);
        Assert.Null(s.Quantities[NumericQuantity.EbsdLocal].Min);
        Assert.Contains("EbsdLocal not evaluated: n 0 (excluded 4), no evaluated values", s.Summary());
    }

    [Fact]
    public void FirstFailure_PrefersTheEarlierStageOfTheSameDirection_RegardlessOfOrder()
    {
        var late = new NumericStatsAccumulator();
        Assert.NotNull(late.Status(NumericStage.Boundary, NativeWrapper.NativeSolveStatus.Singular, 0, 5));
        Assert.NotNull(late.Status(NumericStage.Evd, NativeWrapper.NativeSolveStatus.EvdNoConvergence, 0, 5));
        Assert.Equal(NumericStage.Evd, late.FirstFailure?.Stage);

        var other = new NumericStatsAccumulator();
        Assert.NotNull(other.Status(NumericStage.Integrator, NativeWrapper.NativeSolveStatus.Args, 0, 3));
        foreach (var order in new[] { new[] { late, other }, new[] { other, late } })
        {
            var c = new NumericRunStatsCollector("EBSD");
            foreach (var a in order) c.Merge(a);
            Assert.Equal(new NumericFailure(0, 3, NumericStage.Integrator, NativeWrapper.NativeSolveStatus.Args), c.Build(NumericRunState.Failed).FirstFailure);
        }
        Assert.Null(new NumericStatsAccumulator().Status(NumericStage.Evd, NativeWrapper.NativeSolveStatus.Ok, 0, 0));
    }

    [Fact]
    public void Observe_ReportsTheFirstNonFiniteValueAsAFailure()
    {
        var acc = new NumericStatsAccumulator();
        Assert.Null(acc.Observe(NumericQuantity.EbsdLocal, NumericObservation.Raw, [1.0, -1e-14, 2.0], 0, 7));
        var f = acc.Observe(NumericQuantity.EbsdTds, NumericObservation.Raw, [1.0, double.NaN, double.PositiveInfinity], 1, 2);
        Assert.Equal(new NumericFailure(1, 2, NumericStage.ValueCheck, NativeWrapper.NativeSolveStatus.NonFiniteOutput, NumericQuantity.EbsdTds, 1), f);
        var c = new NumericRunStatsCollector("EBSD"); c.Merge(acc);
        var s = c.Build(NumericRunState.Failed);
        Assert.Equal(1, s.Quantities[NumericQuantity.EbsdLocal].Negative);
        Assert.Equal(2, s.Quantities[NumericQuantity.EbsdTds].NonFinite);
        Assert.Equal(1.0, s.Quantities[NumericQuantity.EbsdTds].Min);
        Assert.Equal(f, s.FirstFailure);
    }

    [Fact]
    public void PublishedStats_DoNotChangeAfterLaterMerges()
    {
        var c = new NumericRunStatsCollector("EBSD");
        var a = new NumericStatsAccumulator();
        a.Observe(NumericQuantity.EbsdLocal, NumericObservation.Raw, [-1e-3, 1.0], 0, 0);
        a.Status(NumericStage.Evd, NativeWrapper.NativeSolveStatus.Ok, 0, 0);
        c.Merge(a);
        var first = c.Build(NumericRunState.Completed);
        var b = new NumericStatsAccumulator();
        b.Observe(NumericQuantity.EbsdLocal, NumericObservation.Raw, [-0.5, 1.0], 0, 1);
        b.Status(NumericStage.Evd, NativeWrapper.NativeSolveStatus.Ok, 0, 1);
        c.Merge(b);
        var local = first.Quantities[NumericQuantity.EbsdLocal];
        Assert.Equal(1, local.Negative); Assert.Equal(2, local.Evaluated); Assert.Equal(1, local.NegativeDecades.Sum());
        Assert.Equal(1, first.Statuses[(NumericStage.Evd, NativeWrapper.NativeSolveStatus.Ok)]);
        Assert.IsNotType<long[]>(local.NegativeDecades); // 読み取り専用の写し (配列のまま渡さない)
        Assert.Equal(2, c.Build(NumericRunState.Completed).Quantities[NumericQuantity.EbsdLocal].Negative);
    }
}
