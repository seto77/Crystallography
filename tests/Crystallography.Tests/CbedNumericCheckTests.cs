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

    private static (Exception Error, BetheMethod.CbedSimulationResult Result, BetheMethod Owner) Run(BetheMethod.Solver solver, Action<BetheMethod> configure = null)
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
            bethe.RunCBED(24, 200, rot, [30, 60], Grid(8, 2e-3), false, solver, 4);
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
        Assert.True(s.Directions > 0);
        Assert.Null(s.FirstFailure);
    }
}
