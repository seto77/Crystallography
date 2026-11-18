using Xunit;

namespace Crystallography.Tests;

/// <summary>261118Cl 新規 (codex2 第 9 回): 平行ビーム (SAED) と PED の検査つきの EVD・逆行列。失敗は NumericFailureException で、
/// 失敗した後に同じ条件でやり直すと、新しく作った BetheMethod と同じ値になる (途中の固有解・振幅を有効なものとして使い回さない)</summary>
[Collection("Bethe global settings")]
public class DiffractionNumericCheckTests
{
    private static Dictionary<(int, int, int), System.Numerics.Complex> Psi(BetheMethod.Beam[] beams) => beams.ToDictionary(b => b.Index, b => b.Psi);

    [NativeNumericFact]
    public void Saed_FailureThenRetry_EqualsFresh()
    {
        var (crystal, rot) = BetheEigenReuseTests.TiltedSi();
        var expected = Psi(new BetheMethod(crystal).GetDifractedBeamAmpriltudes(48, 200, rot, 30));
        var bethe = new BetheMethod(crystal) { DiffractionCheckFlagsForTest = _ => NativeWrapper.NumericCheckFlags.TestForceNoConvergence };
        var ex = Assert.Throws<NumericFailureException>(() => bethe.GetDifractedBeamAmpriltudes(48, 200, rot, 30));
        Assert.Equal(NumericStage.Evd, ex.Failure.Stage);
        Assert.Contains("SAED numeric failure", ex.Message);
        bethe.DiffractionCheckFlagsForTest = null;
        BetheEigenReuseTests.AssertSameAmplitudes(expected, Psi(bethe.GetDifractedBeamAmpriltudes(48, 200, rot, 30)), "SAED retry after a failure");
    }

    [NativeNumericFact]
    public void Ped_OneDirectionFailure_ThenRetry_EqualsFresh()
    {
        var (crystal, rot) = BetheEigenReuseTests.TiltedSi();
        var expected = new BetheMethod(crystal).GetPrecessionElectronDiffraction(32, 200, rot, 30, 0.02, 8).ToDictionary(b => b.Index, b => b.Psi);
        var bethe = new BetheMethod(crystal);
        bethe.GetPrecessionElectronDiffraction(32, 200, rot, 30, 0.02, 8); // 成功 (固有解が作られる)
        bethe.DiffractionCheckFlagsForTest = k => k == 3 ? NativeWrapper.NumericCheckFlags.TestForceNoConvergence : NativeWrapper.NumericCheckFlags.None;
        var ex = Assert.Throws<NumericFailureException>(() => bethe.GetPrecessionElectronDiffraction(32, 200, rot, 30, 0.025, 8)); // 条件を替えて固有解を作り直させる
        Assert.Equal(3, ex.Failure.Direction);
        Assert.Contains("PED numeric failure", ex.Message);
        bethe.DiffractionCheckFlagsForTest = null;
        var retry = bethe.GetPrecessionElectronDiffraction(32, 200, rot, 30, 0.02, 8).ToDictionary(b => b.Index, b => b.Psi);
        BetheEigenReuseTests.AssertSameAmplitudes(expected, retry, "PED retry after a one-direction failure");
    }
}
