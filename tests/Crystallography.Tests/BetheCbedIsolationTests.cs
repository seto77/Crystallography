using Xunit;

namespace Crystallography.Tests;

/// <summary>261118Cl 新規 (codex2 第 3・4 回の単位 5): CBED を専用の BetheMethod で計算すること (<see cref="BetheMethod.RunCBED"/>) の回帰。
/// 以前の RunCBED は受理の前に所有元の MaxNumOfBloch・AccVoltage・BaseRotation・BeamDirections・Thicknesses を書き、計算中も所有元の値を読んだので、
/// 同じ結晶の SAED が割り込むと結果が混ざりえた。二重起動の判定も無かった。完了の通知は SynchronizationContext の無いテストではスレッドプールで届く</summary>
[Collection("Bethe global settings")]
public class BetheCbedIsolationTests
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

    private static (Crystal Crystal, Matrix3D Rotation) Si() => BetheEigenReuseTests.TiltedSi();

    /// <summary>CBED を走らせて完了を待つ。onFirstProgress は最初の進捗の通知で 1 度だけ (計算のスレッドとは別のスレッドで) 呼ぶ</summary>
    private static (BetheMethod.CbedSimulationResult Result, RunWorkerCompletedEventArgsLite Done) Run(BetheMethod bethe, Matrix3D rot, Vector3DBase[] dirs, double[] thick,
        Action onFirstProgress = null, bool cancelAtFirstProgress = false)
    {
        using var done = new ManualResetEventSlim();
        RunWorkerCompletedEventArgsLite lite = default;
        int fired = 0;
        void P(object s, System.ComponentModel.ProgressChangedEventArgs e)
        {
            if (Interlocked.Exchange(ref fired, 1) != 0) return;
            if (cancelAtFirstProgress) bethe.CancelCBED();
            onFirstProgress?.Invoke();
        }
        void C(object s, System.ComponentModel.RunWorkerCompletedEventArgs e)
        {
            lite = new(e.Error, e.Cancelled, e.Error is null && !e.Cancelled ? e.Result as BetheMethod.CbedSimulationResult : null,
                bethe.IsCBED_Busy, ReferenceEquals(bethe.ResultCbed, e.Error is null && !e.Cancelled ? e.Result : null));
            done.Set();
        }
        bethe.CBED_ProgressChanged += P; bethe.CBED_Completed += C;
        try
        {
            bethe.RunCBED(40, 200, rot, thick, dirs, false, BetheMethod.Solver.Eigen_Eigen, 4);
            Assert.True(done.Wait(TimeSpan.FromMinutes(5)), "CBED did not finish");
        }
        finally { bethe.CBED_ProgressChanged -= P; bethe.CBED_Completed -= C; }
        return (lite.Result, lite);
    }

    internal readonly record struct RunWorkerCompletedEventArgsLite(Exception Error, bool Cancelled, BetheMethod.CbedSimulationResult Result, bool BusyAtNotify, bool PublishedBeforeNotify);

    private static void AssertSameCbed(BetheMethod.CbedSimulationResult expected, BetheMethod.CbedSimulationResult actual, string what)
    {
        Assert.NotNull(expected); Assert.NotNull(actual);
        Assert.Equal(expected.Beams.Select(b => b.Index), actual.Beams.Select(b => b.Index));
        Assert.Equal(expected.Disks.Length, actual.Disks.Length);
        for (int t = 0; t < expected.Disks.Length; t++)
            for (int g = 0; g < expected.Disks[t].Length; g++)
            {
                var a = expected.Disks[t][g].RawAmplitudes; var b = actual.Disks[t][g].RawAmplitudes;
                Assert.Equal(a.Length, b.Length);
                for (int r = 0; r < a.Length; r++) Assert.True((a[r] - b[r]).Magnitude <= 1e-12, $"{what}: beam {expected.Beams[g].Index}, pixel {r}, t {t}: {a[r]} vs {b[r]}");
            }
    }

    [Fact]
    public void SaedDuringCbed_LeavesBothUnaffected_AndKeepsOwnerBeams()
    {
        var (crystal, rot) = Si();
        var dirs = Grid(12, 2e-3); double[] thick = [30, 60];
        var saedRot = Matrix3D.Rot(new Vector3DBase(0.3, 1, 0.1), 0.04);
        var expectedCbed = Run(new BetheMethod(crystal), rot, dirs, thick).Result;
        var expectedSaed = new BetheMethod(crystal).GetDifractedBeamAmpriltudes(48, 120, saedRot, 25);

        var owner = new BetheMethod(crystal);
        var ownerBeamsBefore = owner.GetDifractedBeamAmpriltudes(48, 120, saedRot, 25); // 所有元の SAED の状態
        BetheMethod.Beam[] interleaved = null;
        var (actualCbed, done) = Run(owner, rot, dirs, thick, onFirstProgress: () => interleaved = owner.GetDifractedBeamAmpriltudes(48, 120, saedRot, 25));

        Assert.Null(done.Error);
        AssertSameCbed(expectedCbed, actualCbed, "CBED with SAED interleaved");
        Assert.NotNull(interleaved);
        Assert.Equal(expectedSaed.Select(b => b.Index), interleaved.Select(b => b.Index));
        for (int i = 0; i < expectedSaed.Length; i++) Assert.True((expectedSaed[i].Psi - interleaved[i].Psi).Magnitude <= 1e-12);
        Assert.Same(interleaved, owner.Beams); // CBED は所有元の Beams (SAED の結果) を書き換えない
        Assert.Equal(ownerBeamsBefore.Select(b => b.Index), owner.Beams.Select(b => b.Index));
        Assert.True(done.PublishedBeforeNotify, "結果は通知の前に公開する");
        Assert.False(done.BusyAtNotify, "門は通知の前に開ける (完了で次の run を始められる)");
        Assert.Same(actualCbed.Disks, owner.Disks); // 互換の写し
    }

    [Fact]
    public void SecondStart_IsRejected_WithoutTouchingTheRunningOne()
    {
        var (crystal, rot) = Si();
        var dirs = Grid(12, 2e-3); double[] thick = [30, 60];
        var expected = Run(new BetheMethod(crystal), rot, dirs, thick).Result;
        var owner = new BetheMethod(crystal);
        Exception second = null;
        var (actual, done) = Run(owner, rot, dirs, thick, onFirstProgress: () =>
        {
            try { owner.RunCBED(16, 300, Matrix3D.Rot(new Vector3DBase(1, 0, 0), 0.1), [5], Grid(4, 1e-3), false); }
            catch (Exception ex) { second = ex; }
        });
        Assert.IsType<InvalidOperationException>(second);
        Assert.Null(done.Error);
        AssertSameCbed(expected, actual, "CBED after a rejected second start");
    }

    [Fact]
    public void CancelledRun_KeepsPreviousResult_AndRerunMatches()
    {
        var (crystal, rot) = Si();
        var dirs = Grid(12, 2e-3); double[] thick = [30, 60];
        var owner = new BetheMethod(crystal);
        var first = Run(owner, rot, dirs, thick).Result;
        var (_, cancelled) = Run(owner, Matrix3D.Rot(new Vector3DBase(1, 1, 0), 0.01) * rot, dirs, thick, cancelAtFirstProgress: true);
        Assert.True(cancelled.Cancelled);
        Assert.Same(first, owner.ResultCbed); // 取消では前回の結果のまま
        Assert.Same(first.Disks, owner.Disks);
        var rerun = Run(owner, rot, dirs, thick).Result; // 直後の再実行
        AssertSameCbed(first, rerun, "rerun after cancel");
    }

    [Fact]
    public void InvalidInputs_AreRejectedBeforeTheWorkerStarts()
    {
        var (crystal, rot) = Si();
        var owner = new BetheMethod(crystal);
        Assert.Throws<ArgumentException>(() => owner.RunCBED(16, 200, rot, [10], new Vector3DBase[3], false)); // 平方数でない
        Assert.Throws<ArgumentException>(() => owner.RunCBED(16, 200, rot, [10], [new Vector3DBase(double.NaN, 0, -1)], false));
        Assert.Throws<ArgumentException>(() => owner.RunCBED(16, 200, rot, [], Grid(4, 1e-3), false));
        Assert.Throws<ArgumentOutOfRangeException>(() => owner.RunCBED(16, 0, rot, [10], Grid(4, 1e-3), false));
        Assert.False(owner.IsCBED_Busy);
        Assert.Null(owner.ResultCbed);
    }
}
