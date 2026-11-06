using System.Numerics;
using System.Reflection;
using Xunit;

namespace Crystallography.Tests;

/// <summary>261106Ch 追加: 非等Pの2波系で、密度の変化を吸収フラックスと取り違えないこと。</summary>
public class EbsdAbsorbedFluxTests
{
    [Theory]
    [InlineData(2.0, 0.0)]
    [InlineData(2.0, 0.2)]
    [InlineData(0.5, 0.2)]
    public void ReinjectsLostFluxInsteadOfUnweightedDensity(double relativeP, double absorptionRate)
    {
        // 物理振幅は ψ0=cos(2πνz)e^(−az/2), ψ1=i sin(2πνz)/sqrt(P1/P0)e^(−az/2)。
        // P=diag(1,p), A=P^-1 H+i a/(4π) I, H=[[0,c],[c,0]], ν=c/sqrt(p)。
        // 直接の流束は e^(−az)。固有ベクトルはEuclidean直交しないがPで直交する。
        double rootP = Math.Sqrt(relativeP), halfRoot = 1 / Math.Sqrt(2), frequency = 0.17 / rootP;
        Complex[] eigen = [new(frequency, absorptionRate / (4 * Math.PI)), new(-frequency, absorptionRate / (4 * Math.PI))];
        Complex[] vectors = [halfRoot, halfRoot / rootP, halfRoot, -halfRoot / rootP];
        Complex[] alpha = [halfRoot, halfRoot];
        double[] depths = [0, 0.1, 0.7, 1.9, 3.0], weights = [1, relativeP];
        const double sigma = 0.3;
        var method = typeof(BetheMethod).GetMethod("ComputeAbsorbedFluxBackground", BindingFlags.NonPublic | BindingFlags.Static)!;
        var actual = (double[])method.Invoke(null, [2, eigen, vectors, alpha, sigma, depths, weights]);
        var legacy = (double[])method.Invoke(null, [2, eigen, vectors, alpha, sigma, depths, null]);
        for (int t = 0; t < depths.Length; t++)
        {
            // 参照は明示的な2波振幅からの直接求積。modalのFやT行列を使わない。
            const int steps = 4096;
            double sum = 0;
            for (int i = 0; i <= steps; i++)
            {
                double z = depths[t] * i / steps, envelope = Math.Exp(-absorptionRate * z / 2);
                double psi0 = Math.Cos(2 * Math.PI * frequency * z) * envelope;
                double psi1 = Math.Sin(2 * Math.PI * frequency * z) * envelope / rootP;
                double lost = 1 - (psi0 * psi0 + relativeP * psi1 * psi1);
                sum += (i == 0 || i == steps ? 1 : i % 2 == 0 ? 2 : 4) * lost;
            }
            double expected = sigma * depths[t] / (3 * steps) * sum;
            Assert.True(actual[t] >= 0, "Lost flux must be nonnegative.");
            Assert.True(Math.Abs(actual[t] - expected) < 2e-12 * Math.Max(1, depths[t]),
                $"p={relativeP}, absorption={absorptionRate}, t={depths[t]}: {actual[t]:g17} vs {expected:g17}");
        }
        if (absorptionRate == 0)
        {
            Assert.All(actual, value => Assert.True(value < 1e-12));
            Assert.True(legacy[^1] > 0.01, "The fixture must expose the old density/flux error.");
        }
    }
}
