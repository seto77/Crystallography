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
    private static EbsdOutcome Run(bool tds, Action<BetheMethod> configure = null, bool reinject = false)
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
            bethe.RunEBSDNew(16, [20.0], Matrix3D.IdentityMatrix, [10.0, 20.0], Directions(6), BetheMethod.Solver.Eigen_Eigen, 4, includeTDSBackground: tds, includeAbsorbedFluxBackground: reinject);
            Assert.True(done.Wait(TimeSpan.FromMinutes(5)), "EBSD did not finish");
        }
        finally { bethe.EBSD_Completed -= C; }
        return outcome;
    }

    private static long Count(NumericRunStats s, NumericStage stage) => s.Statuses.Where(kv => kv.Key.Stage == stage).Sum(kv => kv.Value);

    /// <summary>managed の経路 (native を切る) でも局所・TDS・再注入を生値で観測する。段は状態を見ないので unchecked。native が無い CI でも走る</summary>
    [Fact]
    public void ManagedRun_ObservesRawSources_AndStagesAreUnchecked()
    {
        bool saved = BetheMethod.EigenEnabled;
        try
        {
            BetheMethod.EigenEnabled = false;
            foreach (var tds in new[] { false, true })
            {
                var r = Run(tds, reinject: true);
                Assert.Null(r.Error);
                var s = r.Stats;
                Assert.Equal(NumericRunState.Completed, s.RunState);
                Assert.Equal("managed-raw", s.IntegratorBackend); Assert.Equal("managed-MathNet", s.EvdBackend); Assert.True(s.RawCapable);
                Assert.Empty(s.Statuses);
                Assert.Contains("Evd unchecked", s.Summary());
                long solved = s.Directions - s.NoBeamDirections;
                Assert.True(solved > 0);
                var source = s.Quantities[tds ? NumericQuantity.EbsdTds : NumericQuantity.EbsdLocal];
                Assert.Equal(NumericObservation.Raw, source.Observation);
                Assert.Equal(solved * 2, source.Evaluated);
                Assert.Equal(NumericObservation.Raw, s.Quantities[NumericQuantity.EbsdReinjection].Observation);
                Assert.Equal(solved * 2, s.Quantities[NumericQuantity.EbsdReinjection].Evaluated);
                Assert.Equal(0, s.Quantities[NumericQuantity.EbsdTotal].NonFinite);
            }
        }
        finally { BetheMethod.EigenEnabled = saved; }
    }

    /// <summary>再注入の重み Σσ が NaN なら、方向を解く前にその電圧で止める (旧は 0 の源として黙って進めた)</summary>
    [Fact]
    public void InvalidReinjectionWeight_FailsTheRun()
    {
        var saved = BetheMethod.BackscatterSourceScaleForResearch;
        try
        {
            BetheMethod.BackscatterSourceScaleForResearch = (_, _) => double.NaN;
            var r = Run(tds: false, reinject: true);
            var ex = Assert.IsType<NumericFailureException>(r.Error);
            Assert.Equal(new NumericFailure(0, -1, NumericStage.ValueCheck, NativeWrapper.NativeSolveStatus.Args, NumericQuantity.EbsdReinjection), ex.Failure);
            Assert.Equal(NumericRunState.Failed, r.Stats.RunState);
            Assert.Equal(0, r.Stats.Directions);
            Assert.Equal(ex.Failure, r.Stats.FirstFailure); // 261118Cl (codex2 第 8 回): 例外と同じ失敗が集計にも残る
        }
        finally { BetheMethod.BackscatterSourceScaleForResearch = saved; }
    }

    /// <summary>261118Cl 追加 (codex2 第 8 回): 局所の源を使う構築では、各席の σ_n の NaN・負は電圧の入口で失敗 (源の非負性の前提)</summary>
    [Fact]
    public void InvalidLocalSigma_FailsTheRun()
    {
        var saved = BetheMethod.BackscatterSourceScaleForResearch;
        try
        {
            BetheMethod.BackscatterSourceScaleForResearch = (z, _) => z == 14 ? -1.0 : 1.0;
            var r = Run(tds: false);
            var ex = Assert.IsType<NumericFailureException>(r.Error);
            Assert.Equal(new NumericFailure(0, -1, NumericStage.ValueCheck, NativeWrapper.NativeSolveStatus.Args, NumericQuantity.EbsdLocal), ex.Failure);
            Assert.Equal(ex.Failure, r.Stats.FirstFailure);
            Assert.Equal(NumericRunState.Failed, r.Stats.RunState);
        }
        finally { BetheMethod.BackscatterSourceScaleForResearch = saved; }
    }

    /// <summary>真値がゼロと分かる条件 (codex2 第 7 回 §2-3): 二波・無吸収・非等価な P の系は流束 Σ(P_g/P_0)|ψ_g|² が保存し、再注入の真値は厳密に 0。
    /// 生値は丸めの範囲の大きさ (|x| ≤ 1e-12·Σσ·t) に収まり、旧い入口 (包み) は生値に Math.Max(0, ·) を当てたものに等しい</summary>
    [Theory]
    [InlineData(2.0)]
    [InlineData(0.5)]
    public void ReinjectionRaw_IsRoundoffSmall_WhenTheTrueValueIsZero(double relativeP)
    {
        double rootP = Math.Sqrt(relativeP), halfRoot = 1 / Math.Sqrt(2), frequency = 0.17 / rootP;
        System.Numerics.Complex[] eigen = [new(frequency, 0), new(-frequency, 0)];
        System.Numerics.Complex[] vectors = [halfRoot, halfRoot / rootP, halfRoot, -halfRoot / rootP];
        System.Numerics.Complex[] alpha = [halfRoot, halfRoot];
        double[] depths = [0, 0.1, 0.7, 1.9, 3.0, 40, 400], weights = [1, relativeP];
        const double sigma = 0.3;
        var raw = BetheMethod.ComputeAbsorbedFluxBackgroundRaw(2, eigen, vectors, alpha, sigma, depths, weights);
        var legacy = (double[])typeof(BetheMethod).GetMethod("ComputeAbsorbedFluxBackground", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Static)!
            .Invoke(null, [2, eigen, vectors, alpha, sigma, depths, weights]);
        for (int t = 0; t < depths.Length; t++)
        {
            Assert.True(double.IsFinite(raw[t]));
            // 261118Cl 変更 (codex2 第 8 回): 深さごとの Σσ·t で判定し、t = 0 は厳密に 0 (旧は max(t, 1) で薄い深さの幅を見ていなかった)
            if (depths[t] == 0) Assert.Equal(0.0, raw[t]);
            else Assert.True(Math.Abs(raw[t]) <= BetheMethod.EbsdReinjectionNegativeTolerance * sigma * depths[t], $"depth {depths[t]}: raw {raw[t]:e3} vs {BetheMethod.EbsdReinjectionNegativeTolerance * sigma * depths[t]:e3}");
            Assert.Equal(BitConverter.DoubleToInt64Bits(Math.Max(0.0, raw[t])), BitConverter.DoubleToInt64Bits(legacy[t]));
        }
        Assert.All(BetheMethod.ComputeAbsorbedFluxBackgroundRaw(2, eigen, vectors, alpha, 0.0, depths, weights), x => Assert.Equal(0.0, x)); // Σσ = 0 は正当なゼロの源
    }

    [NativeNumericFact]
    public void NormalRun_IsCompleted_AndEveryStageIsCounted()
    {
        foreach (var tds in new[] { false, true })
        {
            var r = Run(tds);
            Assert.Null(r.Error); Assert.False(r.Cancelled);
            var s = r.Stats;
            Assert.Equal(NumericRunState.Completed, s.RunState);
            Assert.Equal("native-raw", s.IntegratorBackend); Assert.Equal("native-checked", s.EvdBackend); Assert.True(s.RawCapable);
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

    // ---- 261118Cl 追加 (codex2 第 8 回): 解析値が決まる負の源の対照 (native と managed の生値、受け取り側の記録) ----

    /// <summary>一波・一席・単位振幅・吸収なし (γ 実数で λ = 0、F = t) で σ &lt; 0 なら、局所の源の生値は σ·t (負)。native と managed の両方</summary>
    [Fact]
    public void LocalSource_NegativeSigma_GivesSigmaTimesDepth()
    {
        System.Numerics.Complex[] eigen = [new(0.04, 0)], vectors = [1], alpha = [1], phase = [1];
        double[] sigma = [-0.5], depths = [0, 1.5, 4];
        var managed = BetheMethod.EBSDSolverManagedRaw(1, 1, eigen, vectors, alpha, phase, sigma, depths);
        for (int t = 0; t < depths.Length; t++) Assert.Equal(-0.5 * depths[t], managed[t], 12);
        if (NativeWrapper.NumericChecksAvailable)
        {
            var native = new double[depths.Length];
            Assert.Equal(NativeWrapper.NativeSolveStatus.Ok, NativeWrapper.EBSDSolverRaw(1, eigen, vectors, alpha, phase, sigma, depths, native));
            for (int t = 0; t < depths.Length; t++) Assert.Equal(-0.5 * depths[t], native[t], 12);
        }
        // 受け取り側: 局所の源は判定しない (診断だけ) が、負値の数と最悪の場所 (源・深さ) は残る
        var acc = new NumericStatsAccumulator();
        double[] scale = [0, 0.5 * 1.5, 0.5 * 4]; // |Σσ|·t
        Assert.Null(acc.Observe(NumericQuantity.EbsdLocal, NumericObservation.Raw, managed, 0, 3, scale));
        var c = new NumericRunStatsCollector("EBSD"); c.Merge(acc);
        var local = c.Build(NumericRunState.Completed).Quantities[NumericQuantity.EbsdLocal];
        Assert.Equal(2, local.Negative); Assert.Equal(0, local.SignificantNegative);
        Assert.Equal(3, local.WorstNegative?.Direction); Assert.Equal(1.0, local.WorstNegativeRelative, 12);
    }

    /// <summary>二波・伝播一定 (γ = 0)・固有ベクトル単位・入射 (1, −1)/√2 で、源 μ = [[u, b], [b, u]] の TDS の生値は tdsCoeff·t·(u − b)。
    /// b &gt; u &gt; 0 は対角も尺度も正のまま負 (エルミートでも半正定値でない)、b = u は半正定値・暗い状態で真値 0</summary>
    [Theory]
    [InlineData(0.3, 0.7)]
    [InlineData(0.5, 0.5)]
    public void TdsSource_TwoBeamControl_GivesAnalyticValue(double u, double b)
    {
        double h = 1 / Math.Sqrt(2), coeff = 0.9;
        System.Numerics.Complex[] eigen = [0, 0], vectors = [1, 0, 0, 1], alpha = [h, -h], mu = [u, b, b, u];
        double[] depths = [0, 2, 7];
        var managed = BetheMethod.ComputeTDSMatrixBackgroundRaw(2, eigen, vectors, alpha, mu, coeff, depths);
        for (int t = 0; t < depths.Length; t++) Assert.Equal(coeff * depths[t] * (u - b), managed[t], 12);
        if (NativeWrapper.NumericChecksAvailable)
        {
            double[] local = new double[depths.Length], tds = new double[depths.Length];
            Assert.Equal(NativeWrapper.NativeSolveStatus.Ok, NativeWrapper.EBSDSolverWithTDSRaw(2, eigen, vectors, alpha, [], [], mu, coeff, depths, local, tds));
            for (int t = 0; t < depths.Length; t++) Assert.Equal(coeff * depths[t] * (u - b), tds[t], 12);
        }
    }

    /// <summary>判定つきの観測 (再注入) では、尺度に対して許容を越える負が SignificantNegative の失敗になり、源と深さが残る。許容以下の負は数えるだけ</summary>
    [Fact]
    public void Observe_WithTolerance_FailsOnlySignificantNegatives()
    {
        var acc = new NumericStatsAccumulator();
        double[] scale = [0, 1, 1, 1];
        Assert.Null(acc.Observe(NumericQuantity.EbsdReinjection, NumericObservation.Raw, [-1e-3, -5e-13, 2, 1], 0, 4, scale, 1e-12)); // 尺度 0 の深さは判定しない
        var f = acc.Observe(NumericQuantity.EbsdReinjection, NumericObservation.Raw, [0, -5e-13, -2e-12, 1], 1, 5, scale, 1e-12);
        Assert.Equal(new NumericFailure(1, 5, NumericStage.ValueCheck, NativeWrapper.NativeSolveStatus.SignificantNegative, NumericQuantity.EbsdReinjection, 2), f);
        var c = new NumericRunStatsCollector("EBSD"); c.Merge(acc);
        var s = c.Build(NumericRunState.Failed);
        Assert.Equal(1, s.Quantities[NumericQuantity.EbsdReinjection].SignificantNegative);
        Assert.Equal(4, s.Quantities[NumericQuantity.EbsdReinjection].Negative);
        Assert.Contains("significant negative 1", s.Summary());
        Assert.Contains("(value ", s.Summary());
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
