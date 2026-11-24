using System.Drawing;
using System.Drawing.Imaging;
using Xunit;

namespace Crystallography.Tests;

/// <summary>261125Cl 追加 (IPAnalyzer の依頼 2026-11-25): GetTiltAndOffset の u の標本標準偏差と、PseudoBitmap の呼び出し側所有の描画 (RenderInto / CreateImage)</summary>
public class IPAnalyzerRequestTests
{
    private static readonly PointD[] centers = [new(0.012, 0.021), new(0.031, -0.008), new(-0.018, 0.016), new(0.004, 0.001)];
    private static readonly double[] radii = [10, 20, 30, 40];
    private const double cameraLength = 100;

    /// <summary>組ごとのオフセットを 6 引数版で求め直して u = X − Y·tanξ の標本標準偏差 (平均を引く二乗和) を作る</summary>
    private static double ExpectedUDeviation(double ksi) // 261125Cl 追加
    {
        var u = new List<double>();
        for (int i = 0; i < centers.Length; i++)
            for (int j = i + 1; j < centers.Length; j++)
            {
                PointD o = new(0, 0);
                double t = 0, p = 0;
                Geometry.GetTiltAndOffset([centers[i], centers[j]], [radii[i], radii[j]], cameraLength, ref o, ref t, ref p);
                u.Add(o.X - o.Y * Math.Tan(ksi));
            }
        var mean = u.Average();
        return Math.Sqrt(u.Sum(v => (v - mean) * (v - mean)) / (u.Count - 1));
    }

    [Theory]
    [InlineData(0.0)]
    [InlineData(0.05)]
    [InlineData(-0.2)]
    public void GetTiltAndOffset_ReturnsSampleDeviationOfU(double ksi) // 261125Cl 追加
    {
        PointD offset = new(0, 0), offsetDev = new();
        double tau = 0, tauDev = 0, phi = 0, phiDev = 0;
        var sU = Geometry.GetTiltAndOffset(centers, radii, cameraLength, ref offset, ref offsetDev, ref tau, ref tauDev, ref phi, ref phiDev, ksi);
        var expected = ExpectedUDeviation(ksi);
        Assert.True(double.IsFinite(sU));
        Assert.Equal(expected, sU, 1e-12 * Math.Max(1, Math.Abs(expected)));
        if (ksi == 0) // ξ = 0 なら X の標本標準偏差 (offsetDev.X は従来の式なので丸め誤差の範囲で一致)
            Assert.Equal(offsetDev.X, sU, 1e-9 * Math.Max(1e-12, Math.Abs(sU)) + 1e-15);
    }

    [Fact]
    public void GetTiltAndOffset_FewerThanTwoPairsIsNaN() // 261125Cl 追加
    {
        PointD offset = new(0, 0), offsetDev = new();
        double tau = 0, tauDev = 0, phi = 0, phiDev = 0;
        var sU = Geometry.GetTiltAndOffset([centers[0], centers[1]], [radii[0], radii[1]], cameraLength, ref offset, ref offsetDev, ref tau, ref tauDev, ref phi, ref phiDev, 0.1);
        Assert.True(double.IsNaN(sU)); // 組は 1 つだけ
    }

    [Fact]
    public void GetTiltAndOffset_PairsWithInvalidRadiusAreNotCounted() // 261125Cl 追加 (codex2 の査読): 半径が NaN の組は下位の探索が (0, 0) を返すが、u には数えない
    {
        PointD offset = new(0, 0), offsetDev = new();
        double tau = 0, tauDev = 0, phi = 0, phiDev = 0;
        var sU = Geometry.GetTiltAndOffset(centers[..3], [10, 20, double.NaN], cameraLength, ref offset, ref offsetDev, ref tau, ref tauDev, ref phi, ref phiDev, 0.0);
        Assert.True(double.IsNaN(sU)); // 有効なのは (10, 20) の 1 組だけ
        sU = Geometry.GetTiltAndOffset(centers[..3], [double.NaN, double.NaN, double.NaN], cameraLength, ref offset, ref offsetDev, ref tau, ref tauDev, ref phi, ref phiDev, 0.0);
        Assert.True(double.IsNaN(sU)); // 旧の数え方では u = [0, 0, 0] で 0 になっていた
    }

    private static PseudoBitmap MakeBitmap() // 261125Cl 追加
    {
        var values = Enumerable.Range(0, 64 * 48).Select(i => (double)(i % 97)).ToArray();
        return new PseudoBitmap(values, 64) { MinValue = 0, MaxValue = 96 };
    }

    private static void AssertSamePixels(Bitmap a, Bitmap b) // 261125Cl 追加
    {
        Assert.Equal(a.Size, b.Size);
        for (int y = 0; y < a.Height; y++)
            for (int x = 0; x < a.Width; x++)
                Assert.Equal(a.GetPixel(x, y).ToArgb(), b.GetPixel(x, y).ToArgb());
    }

    [Fact]
    public void CreateImage_MatchesGetImageAndSurvivesCacheResize() // 261125Cl 追加
    {
        var pb = MakeBitmap();
        var rect = new RectangleD(5, 3, 40, 30);
        var size = new Size(37, 23);
        using var created = pb.CreateImage(rect, size);
        Assert.NotNull(created);
        AssertSamePixels(created, pb.GetImage(rect, size));

        pb.GetImage(rect, new Size(10, 10)); // 内部のキャッシュは作り直される (前のは破棄)
        using var again = pb.CreateImage(rect, size);
        AssertSamePixels(created, again); // 呼び出し側の画像は生きている
        pb.Dispose();
        Assert.Equal(size, created.Size);
    }

    [Fact]
    public void RenderInto_DrawsIntoCallerBitmapAndChecksFormat() // 261125Cl 追加
    {
        var pb = MakeBitmap();
        var rect = new RectangleD(0, 0, 64, 48);
        using var target = new Bitmap(20, 15, pb.ImagePixelFormat);
        Assert.True(pb.RenderInto(rect, target));
        using var reference = pb.CreateImage(rect, new Size(20, 15));
        AssertSamePixels(reference, target);

        using var wrong = new Bitmap(20, 15, PixelFormat.Format32bppArgb);
        Assert.Throws<ArgumentException>(() => pb.RenderInto(rect, wrong));
        Assert.False(pb.RenderInto(new RectangleD(0, 0, 0, 10), target));
        Assert.Null(pb.CreateImage(rect, new Size(0, 5)));
    }

    [Fact]
    public void RenderInto_KnownPixels() // 261125Cl 追加 (codex2 の査読): 3 つの入口が同じ RenderInto を通るので、期待値を外から与えて確かめる
    {
        var pb = new PseudoBitmap([0.0, 96.0], 2) { MinValue = 0, MaxValue = 96 }; // 既定のグレー (線形)
        using var bmp = pb.CreateImage(new RectangleD(0, 0, 2, 1), new Size(2, 1));
        Assert.Equal(Color.FromArgb(0, 0, 0).ToArgb(), bmp.GetPixel(0, 0).ToArgb()); // 最小値 = 黒
        var last = PseudoBitmap.ColorScaleGrayLiner[^1];
        Assert.Equal(Color.FromArgb(last.R, last.G, last.B).ToArgb(), bmp.GetPixel(1, 0).ToArgb()); // 最大値 = 色表の最後
    }

    [Fact]
    public void GetImage_AfterDisposeRebuildsCache() // 261125Cl 追加 (codex2 の査読): 共有する画面の片方が破棄しても、残った側の次の GetImage は落ちない
    {
        var pb = MakeBitmap();
        var rect = new RectangleD(0, 0, 64, 48);
        Assert.NotNull(pb.GetImage(rect, new Size(16, 12)));
        pb.Dispose();
        var again = pb.GetImage(rect, new Size(16, 12));
        Assert.NotNull(again);
        Assert.Equal(16, again.Width);
    }
}
