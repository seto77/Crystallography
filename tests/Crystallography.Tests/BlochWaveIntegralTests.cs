using System.Numerics;
using System.Reflection;
using Xunit;

namespace Crystallography.Tests;

/// <summary>261106Ch 追加: 弱吸収の厚さ積分を、固有モード積分式を使わない Simpson 求積で検証する。</summary>
public class BlochWaveIntegralTests
{
    [Theory]
    [InlineData(0.0, 0.1)]
    [InlineData(1e-14, 0.001)]
    [InlineData(1e-12, 0.1)]
    [InlineData(1e-8, 0.1)]
    [InlineData(1e-4, 0.1)]
    [InlineData(0.1, 0.1)]
    [InlineData(1.0, 0.1)]
    [InlineData(1e-16, 1e16)]
    public void OneAbsorbingMode_AllManagedIntegralsMatchQuadrature(double rate, double thickness)
    {
        // |exp(2πiγz)|² = exp(−rate z)。source=1 なので独立参照はこの実関数の積分。
        Complex[] eigen = [new(0, rate / (4 * Math.PI))];
        Complex[] unit = [Complex.One];
        double[] depths = [0, thickness, 2 * thickness, 4 * thickness];
        var local = (double[])typeof(BetheMethod).GetMethod("EBSDSolverManaged", BindingFlags.NonPublic | BindingFlags.Static)!
            .Invoke(null, [1, 1, eigen, unit, unit, unit, new double[] { 1 }, depths]);
        var tds = (double[])typeof(BetheMethod).GetMethod("ComputeTDSMatrixBackground", BindingFlags.NonPublic | BindingFlags.Static)!
            .Invoke(null, [1, eigen, unit, unit, unit, 1.0, depths]);
        var reduction = new AlchemiReduction(1, eigen, unit, unit, depths);
        var ionization = reduction.Contract(unit);
        var coherent = reduction.CoherentPathLengthNm();
        var kikuchiType = typeof(BetheMethod).Assembly.GetType("Crystallography.KikuchiProfileCalculator", throwOnError: true)!;
        var kikuchi = kikuchiType.GetMethod("Fjj", BindingFlags.NonPublic | BindingFlags.Static)!;
        for (int i = 0; i < depths.Length; i++)
        {
            double expected = Integrate(z => Math.Exp(-rate * z), depths[i]);
            var k = (Complex)kikuchi.Invoke(null, [eigen[0], eigen[0], depths[i]]);
            foreach (double actual in new[] { local[i], tds[i], ionization[i], coherent[i], k.Real })
                Assert.True(Math.Abs(actual - expected) <= 5e-13 * Math.Max(depths[i], 1e-300),
                    $"rate={rate:g17}, t={depths[i]:g17}: expected {expected:g17}, got {actual:g17}");
        }
    }

    [Theory]
    [InlineData(1e-14)]
    [InlineData(1e-8)]
    [InlineData(0.1)]
    public void TwoModes_ComplexInterferenceMatchesDirectWaveQuadrature(double splitting)
    {
        Complex[] eigen = [new(0, 1e-14), new(splitting, 2e-14)];
        Complex[] vectors = [Complex.One, Complex.Zero, Complex.Zero, Complex.One];
        Complex[] alpha = [new(0.6, 0.2), new(-0.1, 0.7)];
        Complex[] phase = [Complex.One, Complex.One];
        double[] depths = [0.001, 0.1, 0.2, 0.7];
        var local = (double[])typeof(BetheMethod).GetMethod("EBSDSolverManaged", BindingFlags.NonPublic | BindingFlags.Static)!
            .Invoke(null, [2, 1, eigen, vectors, alpha, phase, new double[] { 1 }, depths]);
        var source = new Complex[] { Complex.One, Complex.One, Complex.One, Complex.One };
        var ionization = new AlchemiReduction(2, eigen, vectors, alpha, depths).Contract(source);
        for (int i = 0; i < depths.Length; i++)
        {
            double expected = Integrate(z =>
            {
                var psi = alpha[0] * Complex.Exp(2 * Math.PI * Complex.ImaginaryOne * eigen[0] * z)
                    + alpha[1] * Complex.Exp(2 * Math.PI * Complex.ImaginaryOne * eigen[1] * z);
                return psi.Real * psi.Real + psi.Imaginary * psi.Imaginary;
            }, depths[i]);
            Assert.True(Math.Abs(local[i] - expected) <= 1e-12 * depths[i]);
            Assert.True(Math.Abs(ionization[i] - expected) <= 1e-12 * depths[i]);
        }
    }

    /// <summary>261106Ch 追加: [0,t] の直接求積。打切り誤差が丸めより小さくなる滑らかな fixture に使う。</summary>
    private static double Integrate(Func<double, double> value, double thickness)
    {
        const int steps = 4096;
        double sum = value(0) + value(thickness);
        for (int i = 1; i < steps; i++) sum += (i % 2 == 0 ? 2 : 4) * value(thickness * i / steps);
        return thickness / (3 * steps) * sum;
    }
}
