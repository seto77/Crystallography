using System.Numerics;
using Xunit;
using Xunit.Abstractions;

namespace Crystallography.Tests;

/// <summary>261118Cl 新規 (codex2 第 10 回): 菊池の 2 ビームの例外点 (固有ベクトルがほぼ平行、|det C| が小さい)。
/// 旧は単位の基底と対角の伝播を返し、下三角の Jordan 型 (J e₀ ≠ 0) で 2πi z J e₀ の項を落とした。
/// 例外点は対角化を使わない 2×2 の行列指数の伝播で I を求め、その両側の対角化の値とつながることを確かめる</summary>
public class KikuchiTwoBeamExceptionalTests(ITestOutputHelper output)
{
    private static readonly Complex Lambda = new(0.012, 0.0015); // 対角 (吸収あり)
    private const double Q00 = 1, Q11 = 0.6, T = 40;
    private static readonly Complex Q01 = new(0.3, 0.1), Q10 = Complex.Conjugate(new Complex(0.3, 0.1));

    private static KikuchiTwoBeamSolution Solve(Complex a00, Complex a01, Complex a10, Complex a11)
        => KikuchiProfileCalculator.SolveTwoBeam(a00, a01, a10, a11, Q00, Q01, Q10, Q11, T);

    /// <summary>独立の参照: ψ(z) を関数で与え、Simpson の公式 (細かい刻み) で I = Re ∫ Σ ψ_i ψ_i'* Q_ii' dz</summary>
    private static double Reference(Func<double, (Complex, Complex)> psi, int steps = 20000)
    {
        double h = T / steps, sum = 0;
        for (int k = 0; k <= steps; k++)
        {
            var (p0, p1) = psi(k * h);
            var f = (p0 * Complex.Conjugate(p0) * Q00 + p0 * Complex.Conjugate(p1) * Q01 + p1 * Complex.Conjugate(p0) * Q10 + p1 * Complex.Conjugate(p1) * Q11).Real;
            sum += f * (k == 0 || k == steps ? 1 : k % 2 == 1 ? 4 : 2);
        }
        return sum * h / 3;
    }

    private static Complex E(Complex gamma, double z) => Complex.Exp(2 * Math.PI * Complex.ImaginaryOne * gamma * z);

    [Fact]
    public void Generic_PropagationAgreesWithTheModalSum()
    {
        Complex a00 = Lambda, a11 = Lambda + 0.02, a01 = new(0.004, 0.0005), a10 = new(0.003, -0.0002);
        var sol = Solve(a00, a01, a10, a11);
        Assert.False(sol.Exceptional);
        var prop = KikuchiProfileCalculator.TwoBeamIntensityByPropagation(a00, a01, a10, a11, Q00, Q01, Q10, Q11, T);
        Assert.True(Math.Abs(sol.Intensity - prop) <= 1e-10 * Math.Abs(prop), $"modal {sol.Intensity:R} vs propagation {prop:R}");
    }

    [Fact]
    public void LowerJordan_IncludesTheSecularTerm()
    {
        Complex a10 = 0.005;
        var sol = Solve(Lambda, 0, a10, Lambda);
        Assert.True(sol.Exceptional);
        // ψ(z) = e^{2πiλz} (e₀ + 2πi z a10 e₁) (J² = 0 なので厳密)
        var reference = Reference(z => (E(Lambda, z), E(Lambda, z) * 2 * Math.PI * Complex.ImaginaryOne * z * a10));
        Assert.True(Math.Abs(sol.Intensity - reference) <= 1e-9 * Math.Abs(reference), $"Jordan {sol.Intensity:R} vs reference {reference:R}");
        var dropped = Reference(z => (E(Lambda, z), Complex.Zero)); // 旧の素通し (対角の伝播) の値
        Assert.True(Math.Abs(dropped - reference) > 1e-3 * Math.Abs(reference), "the secular term must matter in this test");
    }

    [Fact]
    public void UpperJordan_HasNoSecularTermForTheIncidentBeam()
    {
        var sol = Solve(Lambda, 0.005, 0, Lambda);
        Assert.True(sol.Exceptional);
        var reference = Reference(z => (E(Lambda, z), Complex.Zero)); // J e₀ = 0
        Assert.True(Math.Abs(sol.Intensity - reference) <= 1e-9 * Math.Abs(reference), $"{sol.Intensity:R} vs {reference:R}");
    }

    /// <summary>261118Cl 追加 (codex2 第 11 回): 複素の対角差を持つ下三角 (例外点でない) の厳密な振幅 ψ₀ = e^{2πiλ₁z}、ψ₁ = a10 (e^{2πiλ₁z} − e^{2πiλ₂z})/(λ₁ − λ₂) と、
    /// 細かい Simpson を独立の参照にする。厚い強吸収 (3,000 nm で e^{iθ} は 0 に、cos(x) は ∞ に溢れる) でも伝播は有限で参照と一致する</summary>
    [Theory]
    [InlineData(40.0)]
    [InlineData(3000.0)]
    public void LowerTriangular_StrongAbsorption_MatchesTheExactAmplitudes(double thickness)
    {
        Complex l1 = new(0.01, 0.11), l2 = new(0.02, 0.01), a10 = 0.005;
        const int steps = 400000;
        double h = thickness / steps, reference = 0;
        for (int k = 0; k <= steps; k++)
        {
            double z = k * h;
            var p0 = E(l1, z); var p1 = a10 * (E(l1, z) - E(l2, z)) / (l1 - l2);
            var f = (p0 * Complex.Conjugate(p0) * Q00 + p0 * Complex.Conjugate(p1) * Q01 + p1 * Complex.Conjugate(p0) * Q10 + p1 * Complex.Conjugate(p1) * Q11).Real;
            reference += f * (k == 0 || k == steps ? 1 : k % 2 == 1 ? 4 : 2);
        }
        reference *= h / 3;
        var prop = KikuchiProfileCalculator.TwoBeamIntensityByPropagation(l1, 0, a10, l2, Q00, Q01, Q10, Q11, thickness);
        var modal = KikuchiProfileCalculator.SolveTwoBeam(l1, 0, a10, l2, Q00, Q01, Q10, Q11, thickness);
        output.WriteLine($"t {thickness}: reference {reference:R}, propagation {prop:R}, modal {modal.Intensity:R} (exceptional {modal.Exceptional})");
        Assert.True(double.IsFinite(prop));
        Assert.True(Math.Abs(prop - reference) <= 1e-9 * Math.Abs(reference), $"propagation {prop:R} vs reference {reference:R}");
        Assert.True(Math.Abs(modal.Intensity - reference) <= 1e-9 * Math.Abs(reference), $"modal {modal.Intensity:R} vs reference {reference:R}");
    }

    /// <summary>261118Cl 追加 (codex2 第 12 回): 伝播の入力と止め方。厚さ 0 は正当なゼロ、NaN の厚さは入力の失敗、
    /// 区画の上限まで倍にしても収束しない (細かく振動する) ときは成功の値を返さず失敗</summary>
    [Fact]
    public void Propagation_RejectsInvalidInput_AndUnconvergedIntegrals()
    {
        Assert.Equal(0, KikuchiProfileCalculator.TwoBeamIntensityByPropagation(Lambda, 0, 0.005, Lambda, Q00, Q01, Q10, Q11, 0));
        var nan = Assert.Throws<NumericFailureException>(() => KikuchiProfileCalculator.TwoBeamIntensityByPropagation(Lambda, 0, 0.005, Lambda, Q00, Q01, Q10, Q11, double.NaN));
        Assert.Equal(NumericStage.Input, nan.Failure.Stage);
        var ex = Assert.Throws<NumericFailureException>(() => KikuchiProfileCalculator.TwoBeamIntensityByPropagation(1000, 0.3, 0.3, -1000, Q00, Q01, Q10, Q11, 1000));
        Assert.Equal(NumericStage.Integrator, ex.Failure.Stage);
        Assert.Contains("not converged", ex.Message);
    }

    /// <summary>例外点の両側 (a00 − a11 = ±2ε) の値は例外点の値につながる (差は ε に比例して小さくなる)。
    /// 近い点で対角化の和の桁落ちが残らないことも、伝播の値と比べて確かめる</summary>
    [Fact]
    public void BothSidesOfTheExceptionalPoint_AreContinuous()
    {
        Complex a10 = 0.005;
        double atPoint = Solve(Lambda, 0, a10, Lambda).Intensity;
        foreach (var eps in new[] { 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12 })
            foreach (var sign in new[] { +1, -1 })
            {
                Complex a00 = Lambda + sign * eps, a11 = Lambda - sign * eps;
                var sol = Solve(a00, 0, a10, a11);
                var prop = KikuchiProfileCalculator.TwoBeamIntensityByPropagation(a00, 0, a10, a11, Q00, Q01, Q10, Q11, T);
                double det = Math.Abs(sol.C00.Real * sol.C11.Real - sol.C01.Real * sol.C10.Real);
                output.WriteLine($"eps {sign * eps:e0}: exceptional {sol.Exceptional}, |det| {det:e2}, I {sol.Intensity:R}, rel to propagation {(sol.Intensity - prop) / prop:e2}, rel to the point {(sol.Intensity - atPoint) / atPoint:e2}");
                // 対角化の和の誤差は約 1e-16/|det|² なので、|det| < 1e-3 では伝播へ切り替わる (切り替えの後は伝播と同じ値)
                Assert.Equal(det < KikuchiProfileCalculator.ExceptionalDetThreshold || sol.C01 == 0 && sol.C10 == 0 && sol.C00 == 1, sol.Exceptional);
                Assert.True(Math.Abs(sol.Intensity - prop) <= 1e-9 * Math.Abs(prop), $"eps {sign * eps:e0}: {sol.Intensity:R} vs propagation {prop:R}");
                Assert.True(Math.Abs(sol.Intensity - atPoint) <= 100 * eps * T * Math.Abs(atPoint), $"eps {sign * eps:e0}: {sol.Intensity:R} vs the point {atPoint:R}");
            }
    }
}
