using Xunit;

namespace Crystallography.Tests;

/// <summary>261118Cl 新規 (codex2 第 7 回・実験 (b)): 菊池のプロファイルの I・I_ref の生値 (非有限・負を 0 にする前) の観測。
/// 見える値は変えず、各サンプル点の生値を集計すること (2 ビーム版と systematic row 版)、点の診断の I = max(0, 生値) であることを確かめる</summary>
public class KikuchiRawStatsTests
{
    private static KikuchiBandFamily Family(Crystal c, int h, int k, int l)
    {
        var v = c.MatrixInverseTransposed * (h, k, l);
        return new KikuchiBandFamily { Index = (h, k, l), Vec = new Vector3DBase(v.X, v.Y, v.Z), Text = $"{h} {k} {l}", RelativeIntensity = 1 };
    }

    [Theory]
    [InlineData(0)]
    [InlineData(2)]
    public void Profile_ObservesRawIntensities_WithoutChangingTheContrast(int rowOrder)
    {
        var crystal = EbsdBackSourceAxisTests.GaAs(0.6);
        var snap = new KikuchiPotentialSnapshot(crystal, 20);
        var rot = Matrix3D.Rot(new Vector3DBase(1, 0, 0), 0.35) * Matrix3D.Rot(new Vector3DBase(0, 1, 0), 0.22);
        var opt = new KikuchiProfileCalculator.Options { RowOrder = rowOrder };
        var families = new[] { Family(crystal, 2, 2, 0), Family(crystal, 1, 1, 1), Family(crystal, 4, 0, 0) };
        var profiles = families.Select(f => KikuchiProfileCalculator.ComputeProfile(snap, f, rot, 30, opt)).ToArray();
        var valid = profiles.Where(p => p.Valid).ToArray();
        Assert.NotEmpty(valid);
        foreach (var p in valid)
        {
            var s = p.NumericStats;
            Assert.NotNull(s);
            Assert.Equal("Kikuchi", s.Kind);
            Assert.Equal(rowOrder > 0 ? "managed-row" : "managed-2beam", s.IntegratorBackend);
            Assert.Equal(NumericObservation.Raw, s.Quantities[NumericQuantity.KikuchiI].Observation);
            Assert.Equal(p.X.Length, s.Quantities[NumericQuantity.KikuchiI].Evaluated);
            Assert.Equal(p.X.Length, s.Quantities[NumericQuantity.KikuchiIRef].Evaluated);
            Assert.Contains("KikuchiI raw", s.Summary());
            Assert.DoesNotContain("Ebsd", s.Summary()); // 種類の違う源は出さない
        }
        Assert.Null(profiles.FirstOrDefault(p => !p.Valid)?.NumericStats);
        var combined = KikuchiProfileCalculator.CombineNumericStats(profiles);
        Assert.Equal(valid.Length, combined.Directions);
        Assert.Equal(valid.Sum(p => p.X.Length), combined.Quantities[NumericQuantity.KikuchiI].Evaluated);
        if (rowOrder == 0)
            foreach (var x in new[] { -1.5, -1.0, -0.3, 0.0, 0.7, 1.0 })
            {
                var d = KikuchiProfileCalculator.ComputePointDiagnostics(snap, families[0], rot, 30, x, +1, opt);
                Assert.Equal(double.IsFinite(d.RawI) ? Math.Max(0, d.RawI) : 0, d.I);
                Assert.Equal(double.IsFinite(d.RawIRef) ? Math.Max(0, d.RawIRef) : 0, d.IRef);
            }
    }
}
