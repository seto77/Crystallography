using Xunit;

namespace Crystallography.Tests;

/// <summary>261118Cl 新規 (codex2 第 7・8 回の順の 4): CBED の検査つきの伝播と最終の振幅の有限性の検査。
/// Eigen_Eigen は検査つきの入口 (EVD → α → 伝播の状態) を通り、段ごとの数が結果の集計に残ること、EVD の未収束を注入すると
/// その方向で NumericFailureException になり結果を公開しないこと、行列指数 (Auto の既定) でも完了して経路が記録されることを確かめる</summary>
[Collection("Bethe global settings")]
public class CbedNumericCheckTests
{
    private static Vector3DBase[] Grid(int division, double sin)
    {
        var dirs = new Vector3DBase[division * division]; double radius = division / 2.0;
        for (int h = 0; h < division; h++)
            for (int w = 0; w < division; w++)
            {
                double x = (w - radius + 0.5) / (radius - 0.5) * sin, y = -(h - radius + 0.5) / (radius - 0.5) * sin;
                dirs[h * division + w] = new Vector3DBase(x, y, -Math.Sqrt(1 - x * x - y * y));
            }
        return dirs;
    }

    private static (Exception Error, BetheMethod.CbedSimulationResult Result, BetheMethod Owner) Run(BetheMethod.Solver solver, Action<BetheMethod> configure = null, double[] thickness = null)
    {
        var (crystal, rot) = BetheEigenReuseTests.TiltedSi();
        var bethe = new BetheMethod(crystal);
        configure?.Invoke(bethe);
        using var done = new ManualResetEventSlim();
        Exception error = null; BetheMethod.CbedSimulationResult result = null;
        void C(object s, System.ComponentModel.RunWorkerCompletedEventArgs e) { error = e.Error; result = e.Error is null && !e.Cancelled ? e.Result as BetheMethod.CbedSimulationResult : null; done.Set(); }
        bethe.CBED_Completed += C;
        try
        {
            bethe.RunCBED(24, 200, rot, thickness ?? [30, 60], Grid(8, 2e-3), false, solver, 4);
            Assert.True(done.Wait(TimeSpan.FromMinutes(5)), "CBED did not finish");
        }
        finally { bethe.CBED_Completed -= C; }
        return (error, result, bethe);
    }

    [NativeNumericFact]
    public void EigenEigen_UsesTheCheckedPropagation_AndCountsEveryStage()
    {
        var (error, result, _) = Run(BetheMethod.Solver.Eigen_Eigen);
        Assert.Null(error);
        var s = result.NumericStats;
        Assert.NotNull(s);
        Assert.Equal("CBED", s.Kind); Assert.Equal("native-checked", s.IntegratorBackend);
        Assert.True(s.Directions > 0);
        foreach (var stage in new[] { NumericStage.Evd, NumericStage.Boundary, NumericStage.Integrator })
            Assert.Equal(s.Directions, s.Statuses.GetValueOrDefault((stage, NativeWrapper.NativeSolveStatus.Ok)));
        Assert.Null(s.FirstFailure);
    }

    [NativeNumericFact]
    public void EvdFailure_FailsTheRun_AndPublishesNothing()
    {
        var (error, result, owner) = Run(BetheMethod.Solver.Eigen_Eigen, b => b.CbedCheckFlagsForTest = NativeWrapper.NumericCheckFlags.TestForceNoConvergence);
        var ex = Assert.IsType<NumericFailureException>(error);
        Assert.Equal(NumericStage.Evd, ex.Failure.Stage);
        Assert.Equal(NativeWrapper.NativeSolveStatus.EvdNoConvergence, ex.Failure.Status);
        Assert.Contains("CBED numeric failure at 200 kV", ex.Message);
        Assert.Null(result);
        Assert.Null(owner.ResultCbed);
    }

    [Fact]
    public void Auto_Completes_AndRecordsThePath()
    {
        var (error, result, _) = Run(BetheMethod.Solver.Auto);
        Assert.Null(error);
        var s = result.NumericStats;
        Assert.Equal(BetheMethod.EigenEnabled ? "native-matexp" : "managed-MathNet-evd", s.IntegratorBackend);
        if (BetheMethod.EigenEnabled) Assert.Equal("none (matrix exponential, output checked)", s.EvdBackend); // 261118Cl (codex2 第 9 回): 行列指数は EVD を使わない
        Assert.True(s.Directions > 0);
        Assert.Null(s.FirstFailure);
    }

    /// <summary>261118Cl 追加 (codex2 第 9 回): 有限でも値が違う場合を見分ける。Auto (行列指数) と Eigen_Eigen (検査つきの EVD) の円盤の振幅が一致し、
    /// 厚さが不等間隔なら Auto は EVD へ振り替わる (経路の記録も)</summary>
    [NativeNumericFact]
    public void Auto_AgreesWithEigen_AndFallsBackForNonUniformThickness()
    {
        var auto = Run(BetheMethod.Solver.Auto).Result;
        var eigen = Run(BetheMethod.Solver.Eigen_Eigen).Result;
        double max = 0, diff = 0;
        for (int t = 0; t < auto.Disks.Length; t++)
            for (int g = 0; g < auto.Disks[t].Length; g++)
            {
                var a = auto.Disks[t][g].RawAmplitudes; var b = eigen.Disks[t][g].RawAmplitudes;
                for (int r = 0; r < a.Length; r++) { max = Math.Max(max, b[r].Magnitude); diff = Math.Max(diff, (a[r] - b[r]).Magnitude); }
            }
        Assert.True(diff <= 1e-8 * max, $"matrix exponential vs EVD: max |Δψ| {diff:e2} of max |ψ| {max:e2}");
        var nonUniform = Run(BetheMethod.Solver.Auto, thickness: [20, 45, 60]).Result;
        Assert.Equal("native-checked", nonUniform.NumericStats.IntegratorBackend);
    }

    /// <summary>261118Cl 追加 (codex2 第 10 回): 解いた振幅の検査の後に作る円盤の値 (対称の写し・重なりの二乗和の平方根) の異常は、
    /// 公開の前の検査で画素・厚さの場所つきの失敗になり、結果を公開しない</summary>
    [Fact]
    public void NonFiniteCompiledDisk_FailsWithPixelAndThickness_AndPublishesNothing()
    {
        int pixel = -1;
        var (error, result, owner) = Run(BetheMethod.Solver.Auto, b => b.CbedCompiledDisksHookForTest = disks =>
        {
            var amps = disks[1][0].Amplitudes;
            pixel = Array.FindIndex(amps, a => a != System.Numerics.Complex.Zero);
            amps[pixel] = new System.Numerics.Complex(double.NaN, 0);
        });
        var ex = Assert.IsType<NumericFailureException>(error);
        Assert.Equal(NumericStage.ValueCheck, ex.Failure.Stage);
        Assert.Equal(NativeWrapper.NativeSolveStatus.NonFiniteOutput, ex.Failure.Status);
        Assert.Equal(pixel, ex.Failure.Direction);
        Assert.Equal(1, ex.Failure.DepthIndex);
        Assert.Contains("after the overlap", ex.Message);
        Assert.Null(result);
        Assert.Null(owner.ResultCbed);
    }

    /// <summary>261118Cl 追加 (codex2 第 11 回): 解いた振幅の写し (RawAmplitudes、対称の写しの位相の後) だけの異常も検査で止まる (重なりの後の値の検査とは別に)</summary>
    [Fact]
    public void NonFiniteRawAmplitude_FailsWithPixelAndThickness()
    {
        int pixel = -1;
        var (error, result, owner) = Run(BetheMethod.Solver.Auto, b => b.CbedCompiledDisksHookForTest = disks =>
        {
            var raw = disks[1][0].RawAmplitudes;
            pixel = Array.FindIndex(raw, a => a != System.Numerics.Complex.Zero);
            raw[pixel] = new System.Numerics.Complex(0, double.PositiveInfinity);
        });
        var ex = Assert.IsType<NumericFailureException>(error);
        Assert.Equal(pixel, ex.Failure.Direction);
        Assert.Equal(1, ex.Failure.DepthIndex);
        Assert.Contains("non-finite amplitude of beam #0", ex.Message);
        Assert.DoesNotContain("after the overlap", ex.Message);
        Assert.Null(result);
        Assert.Null(owner.ResultCbed);
    }
}
