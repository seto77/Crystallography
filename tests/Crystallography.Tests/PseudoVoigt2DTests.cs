using System;
using Xunit;

namespace Crystallography.Tests;

public sealed class PseudoVoigt2DTests // 261126Cl 追加
{
    private const double Eps = 2.2204460492503131e-16;
    private static readonly double C = Math.Pow(2, 2.0 / 3) - 1;
    private static readonly double Ln2 = Math.Log(2);

    private static void Close(double expected, double actual, double quadratureError = 0) // 261126Cl 追加
    {
        double tol = 64 * Eps * Math.Abs(expected) + quadratureError;
        Assert.True(Math.Abs(actual - expected) <= tol,
            $"expected={expected:R}, actual={actual:R}, tolerance={tol:R}");
    }

    [Theory]
    [InlineData(2.0, 2.0)]
    [InlineData(2.0, 0.75)]
    public void ValuesHalfHeightAndLegacyCurve(double hx, double hy) // 261126Cl 追加
    {
        const double x0 = 0.25, y0 = -0.5, t = 0.375, area = 12;
        double cs = Math.Cos(t), sn = Math.Sin(t);
        foreach (var (u, v) in new[] { (0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (-1.5, 0.75), (3.0, -2.0) })
        {
            double[] xy = [x0 + hx * u * cs - hy * v * sn, y0 + hx * u * sn + hy * v * cs];
            double q = u * u + v * v, d = 1 / (Math.PI * hx * hy);
            double g = d * Ln2 * Math.Pow(2, -q);
            double oldL = 0.5 * d * Math.Pow(1 + C * q, -1.5);
            double l = C * oldL;
            Close(area * l, Marquardt.Lorenzian(xy[0], xy[1], x0, y0, hx, hy, t, area));
            if (hx == hy) Close(area * l, Marquardt.Lorenzian(xy[0], xy[1], x0, y0, hx, area));

            foreach (double eta in new[] { 0.0, 0.375, 1.0 })
            {
                double expected = area * ((1 - eta) * g + eta * l);
                double actual = Marquardt.PseudoVoigt(xy, x0, y0, hx, hy, t, eta, area);
                Close(expected, actual);
                if (hx == hy)
                    Close(expected, Marquardt.PseudoVoigt(xy[0], xy[1], x0, y0, hx, eta, area));
                if (q == 1)
                    Close(Marquardt.PseudoVoigt(new[] { x0, y0 }, x0, y0, hx, hy, t, eta, area) / 2, actual);

                double k = 1 - eta + eta / C;
                Close(area * ((1 - eta) * g + eta * oldL),
                    Marquardt.PseudoVoigt(xy, x0, y0, hx, hy, t, eta / (C * k), area * k));
            }
        }
    }

    [Theory]
    [InlineData(0.0)]
    [InlineData(0.375)]
    [InlineData(1.0)]
    public void IntegratedAreaIncludingTail(double eta) // 261126Cl 追加
    {
        const double hx = 2, hy = 0.75, area = 12, qMax = 8;
        const int n = 16384;
        const double step = qMax / n;
        // 楕円座標で角度積分を済ませると面積要素は π hx hy dq。
        double[] xy = new double[2];
        double sum = 0, compensation = 0;
        for (int i = 0; i <= n; i++)
        {
            xy[0] = hx * Math.Sqrt(i * step);
            int weight = i == 0 || i == n ? 1 : (i % 2 == 0 ? 2 : 4);
            double term = weight * Marquardt.PseudoVoigt(xy, 0, 0, hx, hy, 0, eta, area) - compensation;
            double next = sum + term;
            compensation = (next - sum) - term;
            sum = next;
        }
        double inside = Math.PI * hx * hy * step / 3 * sum;
        double tail = area * ((1 - eta) * Math.Pow(2, -qMax) + eta / Math.Sqrt(1 + C * qMax));
        double maxFourth = area * ((1 - eta) * Math.Pow(Ln2, 5) + eta * 945.0 / 32 * Math.Pow(C, 5));
        Close(area, inside + tail, qMax * Math.Pow(step, 4) * maxFourth / 180);
    }

    [Fact]
    public void AllSevenDerivatives() // 261126Cl 追加
    {
        // 正規化式を x0,y0,hx,hy,t,eta,A で記号微分し、80 桁評価後に double に丸めた値。
        foreach (var (eta, expected) in new (double, double[])[]
        {
            (0, [0.5001723196618556, -0.6051676894056406, 0.05236609322932006,
                 0.7358019351272665, -0.8089572966289432, -0.0635265934495308, 0.022811867205068236]),
            (0.375, [0.3774605813944862, -0.4566965001954532, 0.010627189201295318,
                     0.4974981545355265, -0.6104885846448, -0.0635265934495308, 0.020826661159770398]),
            (1, [0.17294101761553707, -0.2092445181784742, -0.05893765084541259,
                 0.10032518688262636, -0.2797073980045613, -0.0635265934495308, 0.017517984417607338])
        })
        {
            double[] actual = Marquardt.PseudoVoigtDiff(new[] { 2.0, -1.0 }, 0.25, -0.5, 1.5, 0.75, 0.375, eta, 12);
            Assert.Equal(7, actual.Length);
            for (int i = 0; i < expected.Length; i++) Close(expected[i], actual[i]);
        }
    }

    [Theory]
    [InlineData(-0.25, 0.0)]
    [InlineData(1.25, 1.0)]
    public void TwoDimensionalFitLimits(double input, double expected) // 261126Cl 追加
    {
        var circle = new Marquardt.Function(Marquardt.FuncType.PV2, 0, 0, 2, input, 12);
        var ellipse = new Marquardt.Function(Marquardt.FuncType.PV2E, 0, 0, 2, 1, 0, input, 12);
        circle.ApplyConstraints();
        ellipse.ApplyConstraints();
        Assert.Equal(expected, circle.Prms[3]);
        Assert.Equal(expected, ellipse.Prms[5]);
    }

    /// <summary>261126Cl 追加: IPAnalyzer の Find Center が使う中心検出 (FitPeakAsPseudoVoigtByMarcal2D) も同じ正規化した式。
    /// その式で作った楕円のピーク (サブピクセルの中心・背景なし) から中心を取り戻す。片側をマスク (NaN) しても取り戻す</summary>
    [Theory]
    [InlineData(0.0, false)]
    [InlineData(0.6, false)]
    [InlineData(0.6, true)]
    public void Marcal2DRecoversTheCenterOfItsOwnModel(double eta, bool maskLeftColumns)
    {
        const int n = 25;
        const double x0 = 12.3, y0 = 11.6, hk = 6, cc = 1.4, area = 50;
        double hx = hk / 2, hy = hk / (2 * Math.Sqrt(cc));
        var data = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                data[i, j] = maskLeftColumns && i < 4 ? double.NaN
                    : Marquardt.PseudoVoigt(new double[] { i, j }, x0, y0, hx, hy, 0, eta, area);
        var center = FittingPeak.FitPeakAsPseudoVoigtByMarcal2D(data);
        Assert.True(Math.Abs(center.X - x0) < 1e-4 && Math.Abs(center.Y - y0) < 1e-4, $"center = ({center.X}, {center.Y})");
    }

    [Fact]
    public void Marcal2DRejectsAllMaskedOrNonPositiveSums()
    {
        var masked = new double[5, 5];
        for (int i = 0; i < 5; i++) for (int j = 0; j < 5; j++) masked[i, j] = double.NaN;
        Assert.True(double.IsNaN(FittingPeak.FitPeakAsPseudoVoigtByMarcal2D(masked).X));
        Assert.True(double.IsNaN(FittingPeak.FitPeakAsPseudoVoigtByMarcal2D(new double[5, 5]).X));
    }

    /// <summary>261126Cl 追加 (codex2 査読 s13 中 1): 一様な背景と片側のマスクがあっても、中心検出は背景の列を解かずに (残差と同じく固定して) 中心を大きく外さない。
    /// 以前の 9 列の正規方程式では X が約 0.015 ずれた (codex2 の独立計算)</summary>
    [Fact]
    public void Marcal2DWithFlatBackgroundAndMask()
    {
        const int n = 25;
        const double x0 = 12.3, y0 = 11.6, hk = 6, cc = 1.4, area = 50;
        double hx = hk / 2, hy = hk / (2 * Math.Sqrt(cc));
        var data = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                data[i, j] = i < 7 ? double.NaN : Marquardt.PseudoVoigt(new double[] { i, j }, x0, y0, hx, hy, 0, 0.6, area) + 0.03; // 一様な背景 0.03
        var center = FittingPeak.FitPeakAsPseudoVoigtByMarcal2D(data);
        Assert.True(Math.Abs(center.X - x0) < 0.01 && Math.Abs(center.Y - y0) < 0.01, $"center = ({center.X}, {center.Y})");
    }

    /// <summary>261126Cl 追加 (codex2 査読 s13 中 2): 範囲外の初期値 (η = 1.25) から始めても、Solve は η を [0,1] の中で返す</summary>
    [Fact]
    public void SolveKeepsEtaInsideTheConstraint()
    {
        var obs = new System.Collections.Generic.List<(double[] x, double y, double w)>();
        for (int i = -6; i <= 6; i++)
            for (int j = -6; j <= 6; j++)
            {
                var xy = new double[] { i, j };
                obs.Add((xy, Marquardt.PseudoVoigt(xy, 0, 0, 2, 1, 0, 1.25, 12), 1));
            }
        var f = new Marquardt.Function(Marquardt.FuncType.PV2E, 0, 0, 2, 1, 0, 1.25, 12);
        var result = Marquardt.Solve(obs.ToArray(), [f], Marquardt.Precision.Low);
        Assert.InRange(result.Prms[0][5], 0.0, 1.0);
    }
}
