using System.Drawing;
using System.Numerics;
using Xunit;
using Xunit.Abstractions;

namespace Crystallography.Tests;

/// <summary>261106Ch 追加: HRTEM の完全コヒーレンス極限、TCC の全ペア和、空の対物絞りを確認する。</summary>
public class HrtemTccTests
{
    /// <summary>261106Ch 追加: native 指定が managed fallback で通っただけかをテスト出力で区別する。</summary>
    public HrtemTccTests(ITestOutputHelper output)
    {
        output.WriteLine($"NativeWrapper.Enabled = {NativeWrapper.Enabled}; EigenEnabled = {BetheMethod.EigenEnabled}");
    }

    // 261106Ch 追加: q の集約を使わず全順序対 (g,h) を直接足す独立な oracle。
    private static double[] DirectTcc(Crystal crystal, BetheMethod.Beam[] beams, Size size, double resolution,
        double voltage, double cs, double beta, double delta, double defocus)
    {
        double lambda = UniversalConstants.Convert.EnergyToElectronWaveLength(voltage);
        var shift = crystal.RotationMatrix * (crystal.A_Axis + crystal.B_Axis + crystal.C_Axis) / 2;
        var image = new double[size.Width * size.Height];
        for (int i = 0; i < image.Length; i++)
        {
            double x = (i % size.Width - size.Width / 2.0) * resolution + shift.X;
            double y = -(i / size.Width - size.Height / 2.0) * resolution + shift.Y;
            Complex intensity = 0;
            foreach (var g in beams)
                foreach (var h in beams)
                {
                    double gx = g.Vec.X, gy = g.Vec.Y, hx = h.Vec.X, hy = h.Vec.Y;
                    double g2 = gx * gx + gy * gy, h2 = hx * hx + hy * hy;
                    double phase = -Math.PI * lambda * (defocus * (g2 - h2) + cs * lambda * lambda * (g2 * g2 - h2 * h2) / 2);
                    double dx = defocus * (gx - hx) + cs * lambda * lambda * (g2 * gx - h2 * hx);
                    double dy = defocus * (gy - hy) + cs * lambda * lambda * (g2 * gy - h2 * hy);
                    double temporal = Math.PI * lambda * delta * (g2 - h2);
                    double envelope = Math.Exp(-Math.PI * Math.PI * beta * beta * (dx * dx + dy * dy) - temporal * temporal / 2);
                    intensity += g.Psi * Complex.Conjugate(h.Psi) * envelope *
                        Complex.FromPolarCoordinates(1, phase + 2 * Math.PI * ((gx - hx) * x + (gy - hy) * y));
                }
            Assert.True(Math.Abs(intensity.Imaginary) < 2e-12, "全順序対の TCC 強度は実数");
            Assert.True(intensity.Real > -2e-12, "TCC 強度は非負");
            image[i] = intensity.Real;
        }
        return image;
    }

    /// <summary>261106Ch 追加: 実結晶・傾斜方位・複数厚み/焦点で、差反射の集約と全順序対の和を比較。</summary>
    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void Tcc_MatchesDirectOrderedPairs(bool native)
    {
        var (crystal, rot) = BetheEigenReuseTests.TiltedSi();
        var bethe = new BetheMethod(crystal);
        var size = new Size(9, 7);
        double[] thicknesses = [0, 17, 43], defocusses = [-45, 0, 20];
        const double resolution = 0.027, cs = 0.8e6, beta = 0.00015, delta = 3;
        bethe.GetHRTEMImage(50, 200, rot, (double.PositiveInfinity, 0, 0), size, resolution,
            cs, beta, delta, thicknesses, defocusses, quasiMode: false, native: native);
        var actual = bethe.ResultHRTEM.Image;
        for (int t = 0; t < thicknesses.Length; t++)
        {
            var beams = bethe.GetDifractedBeamAmpriltudes(50, 200, rot, thicknesses[t]);
            for (int d = 0; d < defocusses.Length; d++)
            {
                var expected = DirectTcc(crystal, beams, size, resolution, 200, cs, beta, delta, defocusses[d]);
                double scale = Math.Max(1, expected.Max(Math.Abs));
                for (int i = 0; i < expected.Length; i++)
                    Assert.True(Math.Abs(expected[i] - actual[t][d][i]) < 2e-11 * scale, $"TCC t={t}, d={d}, i={i}");
            }
        }
    }

    /// <summary>261106Ch 追加: 偏心絞り、強い/負の球面収差、複数厚さと焦点でも光学位相の再利用を独立な全順序対で検証する。</summary>
    [Theory]
    [InlineData(false, -0.6e6)]
    [InlineData(true, -0.6e6)]
    [InlineData(false, 3.0e7)]
    [InlineData(true, 3.0e7)]
    public void Tcc_OffAxisApertureAndLargePhase_MatchesDirectOrderedPairs(bool native, double cs)
    {
        var (crystal, rot) = BetheEigenReuseTests.TiltedSi();
        var bethe = new BetheMethod(crystal);
        var aperture = (R: 0.014, X: 0.007, Y: -0.004);
        var size = new Size(9, 7);
        double[] thicknesses = [0, 31, 83], defocusses = [-1200, 35, 2300];
        const double resolution = 0.031, beta = 0.00003, delta = 1.7;
        bethe.GetHRTEMImage(80, 200, rot, aperture, size, resolution,
            cs, beta, delta, thicknesses, defocusses, quasiMode: false, native: native);
        var actual = bethe.ResultHRTEM.Image;
        for (int t = 0; t < thicknesses.Length; t++)
        {
            var allBeams = bethe.GetDifractedBeamAmpriltudes(80, 200, rot, thicknesses[t]);
            var beams = BetheMethod.ExtractInsideBeams(allBeams, 200, aperture.R, aperture.X, aperture.Y);
            Assert.True(beams.Length > 1 && beams.Length < allBeams.Length, "偏心絞りで反射の一部だけを選ぶ");
            for (int d = 0; d < defocusses.Length; d++)
            {
                var expected = DirectTcc(crystal, beams, size, resolution, 200, cs, beta, delta, defocusses[d]);
                double scale = Math.Max(1, expected.Max(Math.Abs));
                for (int p = 0; p < expected.Length; p++)
                    Assert.True(Math.Abs(expected[p] - actual[t][d][p]) < 2e-11 * scale, $"off-axis TCC t={t}, d={d}, p={p}");
            }
        }
    }

    /// <summary>261106Ch 追加: β=δ=0 で TCC は波動関数の絶対値二乗と一致する。</summary>
    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void Tcc_CoherentLimit_MatchesWaveIntensity(bool native)
    {
        var (crystal, rot) = BetheEigenReuseTests.TiltedSi();
        var bethe = new BetheMethod(crystal);
        bethe.GetHRTEMImage(50, 200, rot, (0.04, 0, 0), new Size(11, 9), 0.025, 0.6e6, 0, 0, [30], [-30, 10], true, native);
        var coherent = bethe.ResultHRTEM.Image;
        bethe.GetHRTEMImage(50, 200, rot, (0.04, 0, 0), new Size(11, 9), 0.025, 0.6e6, 0, 0, [30], [-30, 10], false, native);
        var tcc = bethe.ResultHRTEM.Image;
        for (int d = 0; d < coherent[0].Length; d++)
            for (int i = 0; i < coherent[0][d].Length; i++)
                Assert.True(Math.Abs(coherent[0][d][i] - tcc[0][d][i]) < 2e-11);
    }

    /// <summary>261106Ch 追加: 反射が一つも通らない絞りは、両モデル/両描画経路でゼロ像になる。</summary>
    [Theory]
    [InlineData(false, false)]
    [InlineData(false, true)]
    [InlineData(true, false)]
    [InlineData(true, true)]
    public void EmptyAperture_ReturnsZeroImage(bool quasi, bool native)
    {
        var (crystal, rot) = BetheEigenReuseTests.TiltedSi();
        var bethe = new BetheMethod(crystal);
        bethe.GetHRTEMImage(20, 200, rot, (0, 0.1, 0.1), new Size(7, 5), 0.03, 1e6, 0.0001, 3, [0, 20], [-30, 0], quasi, native);
        Assert.All(bethe.ResultHRTEM.Image.SelectMany(t => t).SelectMany(d => d), value => Assert.Equal(0, value));
    }
}
