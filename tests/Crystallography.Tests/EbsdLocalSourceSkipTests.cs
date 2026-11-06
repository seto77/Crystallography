using System.Numerics;
using Xunit;
using Xunit.Abstractions;

namespace Crystallography.Tests;

/// <summary>261106Ch 追加: native DLL が無い環境で、実行されなかった検証を成功と数えない。</summary>
public sealed class NativeEbsdTheoryAttribute : TheoryAttribute
{
    /// <summary>261106Ch 追加: DLL の有無はテスト発見時に判定し、無ければ明示的なSkipにする。</summary>
    public NativeEbsdTheoryAttribute()
    {
        if (!NativeWrapper.Enabled)
            Skip = "Crystallography.Native is required to verify the zero-local-source native ABI contract.";
    }
}

/// <summary>261106Ch 追加: 既存APIのsigma空配列が、局所出力0・TDS不変を意味することを検証する。</summary>
public class EbsdLocalSourceSkipTests
{
    /// <summary>261106Ch 追加: 実際にnativeがロードされたことをTRXのテスト出力にも残す。</summary>
    public EbsdLocalSourceSkipTests(ITestOutputHelper output)
    {
        output.WriteLine($"NativeWrapper.Enabled = {NativeWrapper.Enabled}; {NativeWrapper.LastLoadError}");
    }

    /// <summary>261106Ch 追加: 両overload、等間隔/非等間隔、弱吸収で、TDSを直接波のSimpson積分と比較。</summary>
    [NativeEbsdTheory]
    [InlineData(false, false, false)]
    [InlineData(false, false, true)]
    [InlineData(false, true, false)]
    [InlineData(false, true, true)]
    [InlineData(true, false, false)]
    [InlineData(true, false, true)]
    [InlineData(true, true, false)]
    [InlineData(true, true, true)]
    public void EmptyAtoms_RemoveOnlyLocalSignal(bool explicitDimension, bool weakAbsorption, bool uniformDepths)
    {
        Assert.True(NativeWrapper.Enabled);
        double absorption = weakAbsorption ? 1e-14 : 0.007;
        Complex[] eigen = [new(0.03, absorption), new(-0.05, 2 * absorption)];
        Complex[] vectors = [Complex.One, new(0.2, 0.1), new(-0.15, 0.05), Complex.One];
        Complex[] alpha = [new(0.7, 0.1), new(0.2, -0.3)];
        // column-major: phase[n,g] と Hermitian source[row,column]。
        Complex[] phase = [Complex.One, Complex.One, new(0.6, 0.8), Complex.FromPolarCoordinates(1, 1.9)];
        double[] sigma = [1.1, 0.6];
        Complex[] source = [new(1.2, 0), new(0.15, -0.22), new(0.15, 0.22), new(0.8, 0)];
        double[] depths = uniformDepths ? [0.2, 0.4, 0.6, 0.8] : [0, 0.003, 0.1, 1.7];
        const double coefficient = 0.7;

        var both = explicitDimension
            ? NativeWrapper.EBSDSolverWithTDS(2, eigen, vectors, alpha, phase, sigma, source, coefficient, depths)
            : NativeWrapper.EBSDSolverWithTDS(eigen, vectors, alpha, phase, sigma, source, coefficient, depths);
        var noLocal = explicitDimension
            ? NativeWrapper.EBSDSolverWithTDS(2, eigen, vectors, alpha, [], [], source, coefficient, depths)
            : NativeWrapper.EBSDSolverWithTDS(eigen, vectors, alpha, [], [], source, coefficient, depths);
        var local = explicitDimension
            ? NativeWrapper.EBSDSolver(2, eigen, vectors, alpha, phase, sigma, depths)
            : NativeWrapper.EBSDSolver(eigen, vectors, alpha, phase, sigma, depths);

        Assert.Equal(depths.Length, noLocal.intensity.Length);
        Assert.Equal(depths.Length, noLocal.tdsIntensity.Length);
        for (int t = 0; t < depths.Length; t++)
        {
            Assert.Equal(0.0, noLocal.intensity[t]);
            Assert.True(double.IsFinite(noLocal.tdsIntensity[t]));
            double scale = Math.Max(depths[t], 1e-300);
            Assert.True(Math.Abs(both.intensity[t] - local[t]) <= 2e-12 * scale,
                "The original nonempty-atom call must still return the local signal.");
            Assert.True(Math.Abs(both.tdsIntensity[t] - noLocal.tdsIntensity[t]) <= 2e-12 * scale,
                "Removing local atoms must leave the independent TDS source unchanged.");
            double reference = DirectSourceIntegral(eigen, vectors, alpha, source, coefficient, depths[t]);
            Assert.True(Math.Abs(noLocal.tdsIntensity[t] - reference) <= 2e-11 * scale,
                $"TDS direct reference t={depths[t]:g17}: actual={noLocal.tdsIntensity[t]:g17}, expected={reference:g17}");
            if (depths[t] > 0) Assert.True(both.intensity[t] > 0);
        }
    }

    /// <summary>261106Ch 追加: 0厚さ・空厚さの両方で配列契約を保つ。既存native DLLでも同じ契約。</summary>
    [NativeEbsdTheory]
    [InlineData(false)]
    [InlineData(true)]
    public void EmptyAtoms_PreserveEmptyAndZeroDepthOutputs(bool explicitDimension)
    {
        Assert.True(NativeWrapper.Enabled);
        Complex[] eigen = [new(0, 1e-14)], unit = [Complex.One];
        foreach (double[] depths in new double[][] { [], [0], [0, 0] })
        {
            var result = explicitDimension
                ? NativeWrapper.EBSDSolverWithTDS(1, eigen, unit, unit, [], [], unit, 1.0, depths)
                : NativeWrapper.EBSDSolverWithTDS(eigen, unit, unit, [], [], unit, 1.0, depths);
            Assert.Equal(depths.Length, result.intensity.Length);
            Assert.Equal(depths.Length, result.tdsIntensity.Length);
            Assert.All(result.intensity, value => Assert.Equal(0.0, value));
            Assert.All(result.tdsIntensity, value => Assert.Equal(0.0, value));
        }
    }

    /// <summary>261106Ch 追加: λペアの解析積分を使わず、ψ(z)†Uψ(z)を直接Simpson積分する。</summary>
    private static double DirectSourceIntegral(Complex[] eigen, Complex[] vectors, Complex[] alpha,
        Complex[] source, double coefficient, double thickness)
    {
        const int steps = 4096;
        double sum = 0;
        for (int step = 0; step <= steps; step++)
        {
            double z = thickness * step / steps;
            var f0 = alpha[0] * Complex.Exp(2 * Math.PI * Complex.ImaginaryOne * eigen[0] * z);
            var f1 = alpha[1] * Complex.Exp(2 * Math.PI * Complex.ImaginaryOne * eigen[1] * z);
            var psi0 = vectors[0] * f0 + vectors[2] * f1;
            var psi1 = vectors[1] * f0 + vectors[3] * f1;
            var u0 = source[0] * psi0 + source[2] * psi1;
            var u1 = source[1] * psi0 + source[3] * psi1;
            double value = (Complex.Conjugate(psi0) * u0 + Complex.Conjugate(psi1) * u1).Real;
            int weight = step == 0 || step == steps ? 1 : step % 2 == 0 ? 2 : 4;
            sum += weight * value;
        }
        return coefficient * thickness * sum / (3 * steps);
    }
}
