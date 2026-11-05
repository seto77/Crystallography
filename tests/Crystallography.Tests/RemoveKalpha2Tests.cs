using Xunit;

namespace Crystallography.Tests;

/// <summary>261105Cl 新規: DiffractionProfile2.RemoveKalpha2 が誤差 Err を保つこと (旧は捨てていた)。
/// 強度の結果は従来と同じで、誤差は元の誤差以上 (引いた量の誤差を二乗和で加える) になることを確かめる。</summary>
public class RemoveKalpha2Tests
{
    private const double Ka1 = 0.1540562, Ka2 = 0.1544390, Ratio = 0.5; // nm (Cu)

    private static Profile Sample(bool withErr)
    {
        var p = new Profile();
        for (int i = 0; i < 400; i++)
        {
            double x = 25 + i * 0.01, y = 100 + 1000 * Math.Exp(-Math.Pow((x - 28.44) / 0.05, 2));
            p.Pt.Add(new PointD(x, y));
            if (withErr)
                p.Err.Add(new PointD(x, Math.Sqrt(y)));
        }
        return p;
    }

    [Fact]
    public void KeepsErrorsAndIntensities()
    {
        var withErr = DiffractionProfile2.RemoveKalpha2(Sample(true), Ka1, Ka2, Ratio);
        var noErr = DiffractionProfile2.RemoveKalpha2(Sample(false), Ka1, Ka2, Ratio);
        var src = Sample(true);

        Assert.Equal(src.Pt.Count, withErr.Err.Count);
        Assert.Empty(noErr.Err);
        for (int i = 0; i < src.Pt.Count; i++)
        {
            Assert.Equal(noErr.Pt[i].Y, withErr.Pt[i].Y); // 強度は誤差の有無に依らない (従来と同じ計算)
            Assert.Equal(src.Err[i].X, withErr.Err[i].X);
            Assert.True(withErr.Err[i].Y >= src.Err[i].Y - 1e-12);
        }
    }
}
