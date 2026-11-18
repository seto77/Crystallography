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

    /// <summary>261118Cl 追加 (codex2 第 10 回): 失敗しても前に公開した Beams・固有解・Psi はそのまま (旧は EVD の前に Beams を差し替えていた)。
    /// 失敗した条件そのもののやり直しは新しく作ったものと同じ</summary>
    [NativeNumericFact]
    public void Saed_Failure_KeepsPublishedResult_AndFailedConditionRetryEqualsFresh()
    {
        var (crystal, rot) = BetheEigenReuseTests.TiltedSi();
        var rot2 = Matrix3D.Rot(new Vector3DBase(0, 1, 0), 0.004) * rot;
        var expected2 = Psi(new BetheMethod(crystal).GetDifractedBeamAmpriltudes(48, 200, rot2, 30));
        var bethe = new BetheMethod(crystal);
        var published = bethe.GetDifractedBeamAmpriltudes(48, 200, rot, 30);
        var publishedPsi = Psi(published);
        var publishedValues = bethe.EigenValues;
        bethe.DiffractionCheckFlagsForTest = _ => NativeWrapper.NumericCheckFlags.TestForceNoConvergence;
        Assert.Throws<NumericFailureException>(() => bethe.GetDifractedBeamAmpriltudes(48, 200, rot2, 30));
        Assert.Same(published, bethe.Beams);
        Assert.Same(publishedValues, bethe.EigenValues);
        BetheEigenReuseTests.AssertSameAmplitudes(publishedPsi, Psi(bethe.Beams), "SAED published amplitudes after a failure");
        bethe.DiffractionCheckFlagsForTest = null;
        BetheEigenReuseTests.AssertSameAmplitudes(expected2, Psi(bethe.GetDifractedBeamAmpriltudes(48, 200, rot2, 30)), "SAED retry of the failed condition");
    }

    /// <summary>261118Cl 追加 (codex2 第 10 回): PED の 1 方位の失敗の後も、前に公開した Beams・BeamsPED・固有解はそのまま。失敗した半角 0.025 そのもののやり直しは新しく作ったものと同じ</summary>
    [NativeNumericFact]
    public void Ped_Failure_KeepsPublishedResult_AndFailedConditionRetryEqualsFresh()
    {
        var (crystal, rot) = BetheEigenReuseTests.TiltedSi();
        var expected025 = new BetheMethod(crystal).GetPrecessionElectronDiffraction(32, 200, rot, 30, 0.025, 8).ToDictionary(b => b.Index, b => b.Psi);
        var bethe = new BetheMethod(crystal);
        var published = bethe.GetPrecessionElectronDiffraction(32, 200, rot, 30, 0.02, 8);
        var publishedPsi = Psi(published);
        var (beamsPED, valuesPED) = (bethe.BeamsPED, bethe.EigenValuesPED);
        var directionPsi = bethe.BeamsPED[5].Select(b => b.Psi).ToArray();
        bethe.DiffractionCheckFlagsForTest = k => k == 3 ? NativeWrapper.NumericCheckFlags.TestForceNoConvergence : NativeWrapper.NumericCheckFlags.None;
        Assert.Throws<NumericFailureException>(() => bethe.GetPrecessionElectronDiffraction(32, 200, rot, 30, 0.025, 8));
        Assert.Same(published, bethe.Beams);
        Assert.Same(beamsPED, bethe.BeamsPED);
        Assert.Same(valuesPED, bethe.EigenValuesPED);
        BetheEigenReuseTests.AssertSameAmplitudes(publishedPsi, Psi(bethe.Beams), "PED published amplitudes after a failure");
        Assert.Equal(directionPsi, bethe.BeamsPED[5].Select(b => b.Psi).ToArray());
        bethe.DiffractionCheckFlagsForTest = null;
        var retry = bethe.GetPrecessionElectronDiffraction(32, 200, rot, 30, 0.025, 8).ToDictionary(b => b.Index, b => b.Psi);
        BetheEigenReuseTests.AssertSameAmplitudes(expected025, retry, "PED retry of the failed condition");
    }
}
