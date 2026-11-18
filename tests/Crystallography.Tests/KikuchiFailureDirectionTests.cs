using Xunit;

namespace Crystallography.Tests;

/// <summary>261118Cl 新規 (codex2 第 16・17 回): <see cref="KikuchiProfileCalculator.ComputeProfiles"/> の外側は、方向 −1 の数値の失敗 (2 ビームの伝播などが投げる)
/// に候補の族の番号を入れて投げ直す。製品に口を足さず、既存の差し替えの口 <see cref="IKikuchiInelasticKernel"/> に方向 −1 の失敗を投げる試験の実装を渡す</summary>
public class KikuchiFailureDirectionTests
{
    private sealed class FailingKernel : IKikuchiInelasticKernel
    {
        public double SourceAmplitude(int atomsIndex, double s2) =>
            throw new NumericFailureException(new NumericFailure(0, -1, NumericStage.Integrator, NativeWrapper.NativeSolveStatus.NonFiniteOutput), "injected two-beam failure");
    }

    [Fact]
    public void FailureWithoutDirection_GetsTheFamilyIndex()
    {
        var crystal = EbsdBackSourceAxisTests.GaAs(0.6);
        var snap = new KikuchiPotentialSnapshot(crystal, 20, new FailingKernel());
        var rot = Matrix3D.Rot(new Vector3DBase(1, 0, 0), 0.35) * Matrix3D.Rot(new Vector3DBase(0, 1, 0), 0.22);
        var family = KikuchiRawStatsTests.Family(crystal, 2, 2, 0);
        var ex = Assert.Throws<NumericFailureException>(() => KikuchiProfileCalculator.ComputeProfiles(snap, [family], rot, 30, 5));
        Assert.Equal(0, ex.Failure.Direction); // 旧は −1 のまま
        Assert.Equal(NumericStage.Integrator, ex.Failure.Stage);
        Assert.Contains("injected two-beam failure", ex.Message);
    }
}
