using System;
using System.Drawing;
using System.Linq;
using System.Threading.Tasks;
using Xunit;
using Xunit.Abstractions;

namespace Crystallography.Tests;

/// <summary>261124Cl 追加 (IPAnalyzer の codex2 査読の所見 1〜5 の検証、作者決定 11-24): Ring の同心円・放射・扇形の積分を、
/// 合成した平面検出器の画像で確かめる。Ring は static の状態を持つので、このクラスの試験は並列に走らせない (Collection)</summary>
[Collection("Ring static state")]
public class RingIntegrationTests(ITestOutputHelper output)
{
    const int N = 201;                 // 画像の一辺 (画素)
    const double Lambda = 0.1;          // 波長 (nm)
    const double FD = 100.0;            // カメラ長 (mm)
    const double Pix = 1.0;             // 画素の大きさ (mm)

    /// <summary>中心 (cx, cy) の平面検出器に、半径 r0 (mm) の細い輪 (幅 sigma mm のガウス) と、弱い一様な背景を描く</summary>
    static IntegralProperty Setup(double cx, double cy, Func<double, double, double> f)
    {
        var ip = new IntegralProperty
        {
            WaveLength = Lambda, FilmDistance = FD, Camera = IntegralProperty.CameraEnum.FlatPanel,
            SrcWidth = N, SrcHeight = N, CenterX = cx, CenterY = cy, PixSizeX = Pix, PixSizeY = Pix,
            IsFull = true, IsRectangle = true, ConcentricMode = true,
        };
        Ring.SrcImgSize = new Size(N, N);
        var img = new double[N * N];
        const int sub = 6;
        for (int y = 0; y < N; y++)
            for (int x = 0; x < N; x++)
            {
                double s = 0;
                for (int j = 0; j < sub; j++)
                    for (int i = 0; i < sub; i++)
                        s += f(x - 0.5 + (i + 0.5) / sub - cx, y - 0.5 + (j + 0.5) / sub - cy);
                img[x + y * N] = s / sub / sub;
            }
        Ring.Intensity = img;
        Ring.IsSpots = new bool[N * N];
        Ring.IsThresholdOver = new bool[N * N];
        Ring.IsThresholdUnder = new bool[N * N];
        Ring.IsValid = new bool[N * N]; // SetMask は IsValid の長さが合わないと何もしない
        Ring.IP = ip; // SetInsideArea の中の SetTiltParameter が static の IP を読む
        Ring.SetInsideArea(ip);
        Ring.SetMask(false, false, false);
        return ip;
    }

    static double TwoThetaOfD(double d) => 2 * Math.Asin(Lambda / 2 / d);
    static double RadiusOfD(double d) => FD * Math.Tan(TwoThetaOfD(d)) / Pix;

    /// <summary>所見 1: d 値モードのビンの付け方。d0 の細い輪のピークが d0 のビンに出ること</summary>
    [Fact]
    public void DspacingModePeakIsAtTheRingD()
    {
        double d0 = 0.20; // nm
        double r0 = RadiusOfD(d0);
        var ip = Setup(100.0, 100.0, (dx, dy) => 1 + 1000 * Math.Exp(-Math.Pow(Math.Sqrt(dx * dx + dy * dy) - r0, 2) / (2 * 0.3 * 0.3)));
        ip.Mode = HorizontalAxis.d; ip.StartDspacing = 0.15; ip.EndDspacing = 0.40; ip.StepDspacing = 0.01;
        Ring.SetInsideArea(ip); Ring.SetMask(false, false, false);
        var p = Ring.GetProfile(ip);
        var peak = p.Pt.MaxBy(q => q.Y);
        output.WriteLine($"r0={r0:f2} px; peak at d={peak.X:f3} Å (expected {d0 * 10:f3})");
        foreach (var q in p.Pt.Where(q => Math.Abs(q.X - d0 * 10) < 0.35)) output.WriteLine($"  d={q.X:f3} I={q.Y:f2}");
        Assert.Equal(d0 * 10, peak.X, 3);
    }

    /// <summary>所見 1 の対照: 角度モードでは同じ輪のピークが 2θ(d0) のビンに出る</summary>
    [Fact]
    public void AngleModePeakIsAtTheRingTwoTheta()
    {
        double d0 = 0.20;
        double r0 = RadiusOfD(d0);
        var ip = Setup(100.0, 100.0, (dx, dy) => 1 + 1000 * Math.Exp(-Math.Pow(Math.Sqrt(dx * dx + dy * dy) - r0, 2) / (2 * 0.3 * 0.3)));
        ip.Mode = HorizontalAxis.Angle; ip.StartAngle = 10 / 180.0 * Math.PI; ip.EndAngle = 40 / 180.0 * Math.PI; ip.StepAngle = 0.5 / 180.0 * Math.PI;
        Ring.SetInsideArea(ip); Ring.SetMask(false, false, false);
        var p = Ring.GetProfile(ip);
        var peak = p.Pt.MaxBy(q => q.Y);
        var expected = TwoThetaOfD(d0) / Math.PI * 180;
        output.WriteLine($"peak at 2θ={peak.X:f3}° (expected {expected:f3})");
        Assert.True(Math.Abs(peak.X - expected) <= 0.25 + 1e-9, $"peak {peak.X} vs {expected}");
    }

    /// <summary>所見 2: 直射の中心が整数画素に重なり、開始 2θ = 0 のとき、扇形に分けた積分が終わること</summary>
    [Fact]
    public async Task SectorIntegrationTerminatesWhenTheCenterIsOnAPixel()
    {
        var ip = Setup(100.0, 100.0, (dx, dy) => 1.0);
        ip.Mode = HorizontalAxis.Angle; ip.StartAngle = 0; ip.EndAngle = 20 / 180.0 * Math.PI; ip.StepAngle = 1 / 180.0 * Math.PI;
        Ring.SetInsideArea(ip); Ring.SetMask(false, false, false);
        // 終わらないときは 60 s で TimeoutException (旧は終わらなかった)
        var result = await Task.Run(() => Ring.GetConcenrticProfilesBySector(ip, 4)).WaitAsync(TimeSpan.FromSeconds(60));
        Assert.Equal(4, result.Length);
    }

    /// <summary>所見 2 の続き: 中心の画素を外したうえで、扇形ごとの強度が有限 (NaN が 0 に化けていない) であること</summary>
    [Fact]
    public void SectorProfilesAreFiniteWhenTheCenterIsOnAPixel()
    {
        var ip = Setup(100.0, 100.0, (dx, dy) => 1.0);
        ip.Mode = HorizontalAxis.Angle; ip.StartAngle = 0; ip.EndAngle = 20 / 180.0 * Math.PI; ip.StepAngle = 1 / 180.0 * Math.PI;
        Ring.SetInsideArea(ip); Ring.SetMask(false, false, false);
        var full = Ring.GetProfile(ip);
        // 出力は立体角の補正を含む (一様な画像で ≈ 1/cos³2θ)。2θ ≤ 0.5° のビン 0 はほぼ 1 (旧は中心画素の NaN でビン 0 が 0)
        Assert.Equal(1.0, full.Pt[0].Y, 2);
        var sectors = Ring.GetConcenrticProfilesBySector(ip, 4);
        // 扇形ごとの値は、全周の値と画素の離散化の範囲 (5 %) で一致する
        foreach (var s in sectors)
            for (int i = 1; i < s.Pt.Count; i++)
                Assert.True(Math.Abs(s.Pt[i].Y / full.Pt[i].Y - 1) < 0.05, $"bin {i}: sector {s.Pt[i].Y} vs full {full.Pt[i].Y}");
    }

    /// <summary>所見 1 の端: d ≤ λ/2 の境界 (π) を含む範囲でも例外・NaN にならない</summary>
    [Fact]
    public void DspacingModeNearHalfWavelengthIsFinite()
    {
        var ip = Setup(100.0, 100.0, (dx, dy) => 1.0);
        ip.Mode = HorizontalAxis.d; ip.StartDspacing = 0.052; ip.EndDspacing = 0.40; ip.StepDspacing = 0.01; // 最初の境界 0.047 nm < λ/2
        Ring.SetInsideArea(ip); Ring.SetMask(false, false, false);
        var p = Ring.GetProfile(ip);
        Assert.All(p.Pt, q => Assert.True(double.IsFinite(q.Y)));
        // 検出器の内側 (d ≥ 0.13 nm) は、立体角の補正 1/cos³2θ に従う
        Assert.All(p.Pt.Where(q => q.X >= 1.3 && q.X <= 3.9), q =>
        {
            double c = Math.Cos(TwoThetaOfD(q.X / 10));
            Assert.True(Math.Abs(q.Y * c * c * c - 1) < 0.03, $"d={q.X}: {q.Y}");
        });
    }

    /// <summary>所見 3: 放射 (Radial) の d 値モードで、中心 d から半幅の 3/4 だけ離れた輪も積分される (旧は画素の除外が半幅の半分だった)</summary>
    [Fact]
    public void RadialDspacingUsesTheFullHalfWidth()
    {
        double Total(double dRing)
        {
            double r0 = RadiusOfD(dRing);
            var ip = Setup(100.0, 100.0, (dx, dy) => 1000 * Math.Exp(-Math.Pow(Math.Sqrt(dx * dx + dy * dy) - r0, 2) / (2 * 0.3 * 0.3)));
            ip.ConcentricMode = false; ip.Mode = HorizontalAxis.d; ip.RadialSectorAngle = 10;
            ip.RadialRadiusDspacing = 0.20; ip.RadialRadiusDspacingRange = 0.02;
            Ring.SetInsideArea(ip); Ring.SetMask(false, false, false);
            return Ring.GetProfile(ip).Pt.Sum(q => q.Y);
        }
        double inner = Total(0.20), outer = Total(0.215);
        output.WriteLine($"ring at d=0.200: {inner:f1}, ring at d=0.215: {outer:f1}");
        Assert.True(outer > 0.5 * inner, $"ring at 0.215 nm: {outer} vs {inner}");
    }

    /// <summary>所見 4: 反時計回りの χ で、+x 方向 (χ = 0°) の斑点が 0° に、+y 方向の斑点が時計回りと逆の角に出る</summary>
    [Fact]
    public void RadialCounterclockwiseChiIsTheMirrorOfClockwise()
    {
        double Peak(double sx, double sy, Ring.Rotation rot)
        {
            var ip = Setup(100.0, 100.0, (dx, dy) => 1000 * Math.Exp(-((dx - sx) * (dx - sx) + (dy - sy) * (dy - sy)) / (2 * 1.5 * 1.5)));
            ip.ConcentricMode = false; ip.Mode = HorizontalAxis.Angle; ip.RadialSectorAngle = 1;
            ip.RadialRadiusAngle = Math.Atan(50 / FD); ip.RadialRadiusAngleRange = 5 / 180.0 * Math.PI;
            Ring.SetInsideArea(ip); Ring.SetMask(false, false, false);
            var old = Ring.ChiRotation;
            try { Ring.ChiRotation = rot; return Ring.GetProfile(ip).Pt.MaxBy(q => q.Y).X; }
            finally { Ring.ChiRotation = old; }
        }
        double cw0 = Peak(50, 0, Ring.Rotation.Clockwise), ccw0 = Peak(50, 0, Ring.Rotation.Counterclockwise);
        double cw90 = Peak(0, 50, Ring.Rotation.Clockwise), ccw90 = Peak(0, 50, Ring.Rotation.Counterclockwise);
        output.WriteLine($"+x: CW {cw0}, CCW {ccw0};  +y: CW {cw90}, CCW {ccw90}");
        Assert.Equal(0, cw0, 6);
        Assert.Equal(0, ccw0, 6); // 旧は 359°
        Assert.Equal(360, cw90 + ccw90, 6);
    }

    /// <summary>所見 5: Gandolfi の最終ビンにも強度が入る (旧は常に 0)</summary>
    [Fact]
    public void GandolfiLastBinIsFilled()
    {
        var ip = Setup(100.0, 100.0, (dx, dy) => 1.0);
        ip.Camera = IntegralProperty.CameraEnum.Gandolfi; ip.GandolfiRadius = FD;
        ip.Mode = HorizontalAxis.Angle; ip.StartAngle = 10 / 180.0 * Math.PI; ip.EndAngle = 30 / 180.0 * Math.PI; ip.StepAngle = 1 / 180.0 * Math.PI;
        Ring.SetInsideArea(ip); Ring.SetMask(false, false, false);
        var p = Ring.GetProfile(ip);
        output.WriteLine(string.Join(" ", p.Pt.Select(q => $"{q.X:f0}:{q.Y:f3}")));
        // ビン数は Ring と同じ切り捨ての式 (20°/1° はラジアンでは 19.999… になり 20 点)
        Assert.Equal((int)((ip.EndAngle - ip.StartAngle) / ip.StepAngle) + 1, p.Pt.Count);
        Assert.Equal(1.0, p.Pt[^1].Y, 3); // 一様な画像なので 1 (旧は 0)
    }
}
