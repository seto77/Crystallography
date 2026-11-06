using Xunit;

namespace Crystallography.Tests;

/// <summary>261106Ch 追加: PED と平行ビーム回折の切替・表面変更で、以前の条件の固有解が混ざらないことを確かめる。</summary>
[Collection("Bethe global settings")]
public class BethePedReuseTests
{
    /// <summary>261106Ch 追加: 方位・電圧・ビーム数ごとに、実際のモード切替を再現する。</summary>
    [Theory]
    [InlineData(0)]
    [InlineData(1)]
    [InlineData(2)]
    public void Ped_AfterParallelDiffractionChangesCondition_MatchesFreshInstance(int changed)
    {
        var (crystal, rotation) = BetheEigenReuseTests.TiltedSi();
        var bethe = new BetheMethod(crystal);
        const int initialN = 32, steps = 4;
        const double initialVoltage = 200, semiangle = 0.015, thickness = 35;
        bethe.GetPrecessionElectronDiffraction(initialN, initialVoltage, rotation, thickness, semiangle, steps);

        var targetRotation = changed == 0 ? Matrix3D.Rot(new Vector3DBase(0.2, 1, 0), 0.035) : rotation;
        var targetVoltage = changed == 1 ? 120 : initialVoltage;
        var targetN = changed == 2 ? 48 : initialN;
        bethe.GetDifractedBeamAmpriltudes(targetN, targetVoltage, targetRotation, thickness);
        var actual = BetheEigenReuseTests.Amplitudes(bethe.GetPrecessionElectronDiffraction(targetN, targetVoltage, targetRotation, thickness, semiangle, steps));
        var expected = BetheEigenReuseTests.Amplitudes(new BetheMethod(crystal).GetPrecessionElectronDiffraction(targetN, targetVoltage, targetRotation, thickness, semiangle, steps));
        BetheEigenReuseTests.AssertSameAmplitudes(expected, actual, $"PED after parallel diffraction, setting {changed}");
    }

    /// <summary>261106Ch 追加: 同じ Surface オブジェクトの成分を変更した場合も、新しい表面で解く。</summary>
    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void Diffraction_AfterSurfaceChangedInPlace_MatchesFreshInstance(bool ped)
    {
        var (crystal, rotation) = BetheEigenReuseTests.TiltedSi();
        var bethe = new BetheMethod(crystal);
        BetheMethod.Beam[] calculate(BetheMethod b) => ped
            ? b.GetPrecessionElectronDiffraction(32, 200, rotation, 35, 0.015, 4)
            : b.GetDifractedBeamAmpriltudes(32, 200, rotation, 35);
        calculate(bethe);
        bethe.Surface.X = Math.Sin(0.4);
        bethe.Surface.Z = -Math.Cos(0.4);
        var actual = BetheEigenReuseTests.Amplitudes(calculate(bethe));
        var fresh = new BetheMethod(crystal) { Surface = new Vector3D(Math.Sin(0.4), 0, -Math.Cos(0.4)) };
        BetheEigenReuseTests.AssertSameAmplitudes(BetheEigenReuseTests.Amplitudes(calculate(fresh)), actual, $"surface change, PED={ped}");
    }

    /// <summary>261106Ch 追加: 厚さだけを変える高速な再利用は維持する。</summary>
    [Fact]
    public void Ped_ThicknessOnlyChange_ReusesEigenSolutionAndMatchesFreshInstance()
    {
        var (crystal, rotation) = BetheEigenReuseTests.TiltedSi();
        var bethe = new BetheMethod(crystal);
        bethe.GetPrecessionElectronDiffraction(32, 200, rotation, 35, 0.015, 4);
        var eigenvalues = bethe.EigenValuesPED;
        var actual = BetheEigenReuseTests.Amplitudes(bethe.GetPrecessionElectronDiffraction(32, 200, rotation, 70, 0.015, 4));
        Assert.Same(eigenvalues, bethe.EigenValuesPED);
        var expected = BetheEigenReuseTests.Amplitudes(new BetheMethod(crystal).GetPrecessionElectronDiffraction(32, 200, rotation, 70, 0.015, 4));
        BetheEigenReuseTests.AssertSameAmplitudes(expected, actual, "PED thickness series");
    }

    /// <summary>261106Ch 追加: 歳差角ゼロの方位平均が平行ビームの強度へ戻る、物理的な極限を確認する。</summary>
    [Fact]
    public void Ped_ZeroSemiangle_MatchesParallelIntensities()
    {
        var (crystal, rotation) = BetheEigenReuseTests.TiltedSi();
        var expected = BetheEigenReuseTests.Amplitudes(new BetheMethod(crystal).GetDifractedBeamAmpriltudes(32, 200, rotation, 35));
        var actual = BetheEigenReuseTests.Amplitudes(new BetheMethod(crystal).GetPrecessionElectronDiffraction(32, 200, rotation, 35, 0, 4));
        Assert.Equal(expected.Count, actual.Count);
        foreach (var (index, psi) in expected)
            Assert.True(Math.Abs(psi.Magnitude * psi.Magnitude - actual[index].Magnitude * actual[index].Magnitude) < 1e-10, $"zero-precession intensity {index}");
    }

    /// <summary>261106Ch 追加: 吸収散乱因子の模型を切り替えたとき、U と固有解の両方を更新する。</summary>
    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void Diffraction_AfterAbsorptionModelChanged_MatchesFreshInstance(bool ped)
    {
        var original = BetheMethod.AbsorptionElasticModel;
        try
        {
            var (crystal, rotation) = BetheEigenReuseTests.TiltedSi();
            var bethe = new BetheMethod(crystal);
            BetheMethod.Beam[] calculate(BetheMethod b) => ped
                ? b.GetPrecessionElectronDiffraction(32, 200, rotation, 75, 0.015, 4)
                : b.GetDifractedBeamAmpriltudes(32, 200, rotation, 75);
            BetheMethod.AbsorptionElasticModel = AbsorptionElasticModel.PengTail;
            var before = BetheEigenReuseTests.Amplitudes(calculate(bethe));
            BetheMethod.AbsorptionElasticModel = AbsorptionElasticModel.TemariHybrid;
            var actual = BetheEigenReuseTests.Amplitudes(calculate(bethe));
            var expected = BetheEigenReuseTests.Amplitudes(calculate(new BetheMethod(crystal)));
            Assert.True(expected.Any(pair => (pair.Value - before[pair.Key]).Magnitude > 1e-8), "模型の変更で数値が変わる条件で確かめる");
            BetheEigenReuseTests.AssertSameAmplitudes(expected, actual, $"absorption model change, PED={ped}");
        }
        finally { BetheMethod.AbsorptionElasticModel = original; }
    }

    /// <summary>261106Ch 追加: 公開 API で ADP を変更した後に明示的に cache を空にすれば、固有解も更新される。</summary>
    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void Diffraction_AfterAdpMutationAndClearUCache_MatchesFreshInstance(bool ped)
    {
        var (crystal, rotation) = BetheEigenReuseTests.TiltedSi();
        var bethe = new BetheMethod(crystal);
        BetheMethod.Beam[] calculate(BetheMethod b) => ped
            ? b.GetPrecessionElectronDiffraction(32, 200, rotation, 75, 0.015, 4)
            : b.GetDifractedBeamAmpriltudes(32, 200, rotation, 75);
        var before = BetheEigenReuseTests.Amplitudes(calculate(bethe));
        crystal.Atoms[0].Dsf.Iso *= 2;
        bethe.ClearUCache();
        var actual = BetheEigenReuseTests.Amplitudes(calculate(bethe));
        var expected = BetheEigenReuseTests.Amplitudes(calculate(new BetheMethod(crystal)));
        Assert.True(expected.Any(pair => (pair.Value - before[pair.Key]).Magnitude > 1e-5), "ADP の変更で強度が変わる条件で確かめる");
        BetheEigenReuseTests.AssertSameAmplitudes(expected, actual, $"ADP mutation + ClearUCache, PED={ped}");
    }
}

/// <summary>261106Ch 追加: プロセス全体の模型を変更するテストを他の計算と同時に走らせない。</summary>
[CollectionDefinition("Bethe global settings", DisableParallelization = true)]
public class BetheGlobalSettingsCollection { }
