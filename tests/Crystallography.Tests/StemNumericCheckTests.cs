using Xunit;

namespace Crystallography.Tests;

/// <summary>261118Cl 新規 (codex2 第 7・8 回の順の 4): STEM の検査つきの伝播 (全スライスの CBEDSolver_EigenChecked と Eigen2Checked) の失敗の伝わり方。
/// EVD の未収束を注入すると、その方向で NumericFailureException (段 = EVD) になって完了の Error に届き、正常な run は完了すること</summary>
[Collection("Bethe global settings")]
public class StemNumericCheckTests
{
    private static Exception Run(NativeWrapper.NumericCheckFlags flags)
    {
        var (crystal, rot) = BetheEigenReuseTests.TiltedSi();
        var bethe = new BetheMethod(crystal) { StemCheckFlagsForTest = flags };
        const int div = 3;
        var sin = Math.Sin(0.010);
        var directions = new Vector3DBase[div * div];
        for (int h = 0; h < div; h++)
            for (int w = 0; w < div; w++)
            {
                double x = (w - div / 2.0 + 0.5) / (div / 2.0 - 0.5) * sin, y = -(h - div / 2.0 + 0.5) / (div / 2.0 - 0.5) * sin;
                directions[h * div + w] = new Vector3DBase(x, y, -Math.Sqrt(1 - x * x - y * y));
            }
        using var done = new ManualResetEventSlim();
        Exception error = null;
        bethe.StemCompleted += (s, e) => { error = e.Error; done.Set(); };
        var context = SynchronizationContext.Current;
        SynchronizationContext.SetSynchronizationContext(null); // 完了通知をスレッドプールへ
        try { bethe.RunSTEM(24, 200, 1.0e6, 3, 2, new System.Drawing.Size(4, 4), 0.05, 0, rot, [20.0], [0.0], directions, 0.010, 0.040, 0.150, BetheMethod.Solver.Eigen_Eigen, 2); }
        finally { SynchronizationContext.SetSynchronizationContext(context); }
        Assert.True(done.Wait(TimeSpan.FromMinutes(5)), "STEM did not finish");
        return error;
    }

    [NativeNumericFact]
    public void NormalRun_Completes()
    {
        Assert.Null(Run(NativeWrapper.NumericCheckFlags.None));
    }

    [NativeNumericFact]
    public void EvdFailure_FailsTheRun_WithTheDirection()
    {
        var ex = Assert.IsType<NumericFailureException>(Run(NativeWrapper.NumericCheckFlags.TestForceNoConvergence));
        Assert.Equal(NumericStage.Evd, ex.Failure.Stage);
        Assert.Equal(NativeWrapper.NativeSolveStatus.EvdNoConvergence, ex.Failure.Status);
        Assert.Contains("STEM numeric failure at 200 kV", ex.Message);
    }
}
