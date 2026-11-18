using Xunit;

namespace Crystallography.Tests;

/// <summary>261118Cl 追加: 検査つきの native の入口が要る Theory (無ければ Skip。CRYSTALLOGRAPHY_REQUIRE_NUMERIC_CHECKS=1 なら走らせて落とす)</summary>
public sealed class NativeNumericTheoryAttribute : TheoryAttribute
{
    public NativeNumericTheoryAttribute()
    {
        if (!NativeWrapper.NumericChecksAvailable && Environment.GetEnvironmentVariable("CRYSTALLOGRAPHY_REQUIRE_NUMERIC_CHECKS") != "1")
            Skip = "Crystallography.Native with the numeric-check entry points is required: " + NativeWrapper.NumericChecksUnavailableReason;
    }
}

/// <summary>261118Cl 新規 (codex2 第 7・8 回の順の 4): STEM の検査つきの伝播 (全スライスの CBEDSolver_EigenChecked と Eigen2Checked) の失敗の伝わり方。
/// EVD の未収束を注入すると、その方向で NumericFailureException (段 = EVD) になって完了の Error に届き、正常な run は完了すること (261118Cl 第 9 回: 全スライスの方式と Eigen2 の方式の両方)</summary>
[Collection("Bethe global settings")]
public class StemNumericCheckTests
{
    // private static Exception Run(NativeWrapper.NumericCheckFlags flags, string allSlice = null) // 261118Cl 変更前
    private static Exception Run(NativeWrapper.NumericCheckFlags flags, string allSlice = null, Action<BetheMethod> configure = null)
    {
        var savedEnv = Environment.GetEnvironmentVariable("RECIPRO_STEM_ALLSLICE"); // 261118Cl (codex2 第 9 回): 全スライスの方式と Eigen2 の方式を明示して両方を通す
        Environment.SetEnvironmentVariable("RECIPRO_STEM_ALLSLICE", allSlice);
        try
        {
        var (crystal, rot) = BetheEigenReuseTests.TiltedSi();
        var bethe = new BetheMethod(crystal) { StemCheckFlagsForTest = flags };
        configure?.Invoke(bethe); // 261118Cl 追加
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
        finally { Environment.SetEnvironmentVariable("RECIPRO_STEM_ALLSLICE", savedEnv); }
    }

    [NativeNumericTheory]
    [InlineData("1")]
    [InlineData("0")]
    public void NormalRun_Completes(string allSlice)
    {
        Assert.Null(Run(NativeWrapper.NumericCheckFlags.None, allSlice));
    }

    [NativeNumericTheory]
    [InlineData("1")]
    [InlineData("0")]
    public void EvdFailure_FailsTheRun_WithTheDirection(string allSlice)
    {
        var ex = Assert.IsType<NumericFailureException>(Run(NativeWrapper.NumericCheckFlags.TestForceNoConvergence, allSlice));
        Assert.Equal(NumericStage.Evd, ex.Failure.Stage);
        Assert.Equal(NativeWrapper.NativeSolveStatus.EvdNoConvergence, ex.Failure.Status);
        Assert.Contains("STEM numeric failure at 200 kV", ex.Message);
    }

    /// <summary>261118Cl 追加 (codex2 第 10 回): 旧経路 (全スライスを物化しない) で作り直したスライスの透過係数の異常は、非弾性の計算に渡す前の検査で
    /// 方向・厚さの場所つきの失敗になる (中心の方向 #4 の最初の厚さ・最初のスライスに注入)</summary>
    [NativeNumericFact]
    public void NonFiniteRegeneratedSlice_FailsWithDirectionAndThickness()
    {
        var ex = Assert.IsType<NumericFailureException>(Run(NativeWrapper.NumericCheckFlags.None, "0", b => b.StemRegeneratedSliceHookForTest = (tc, t, si) =>
        {
            if (t == 0 && si == 0) tc[4 * (tc.Length / 9)] = new System.Numerics.Complex(double.NaN, 0);
        }));
        Assert.Equal(NumericStage.ValueCheck, ex.Failure.Stage);
        Assert.Equal(NativeWrapper.NativeSolveStatus.NonFiniteOutput, ex.Failure.Status);
        Assert.Equal(4, ex.Failure.Direction);
        Assert.Equal(0, ex.Failure.DepthIndex);
        Assert.Contains("regenerated transmission coefficient", ex.Message);
    }
}
