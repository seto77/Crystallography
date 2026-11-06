using System.Numerics;

namespace Crystallography;

/// <summary>261106Ch 追加: Bloch 波の厚さ積分 ∫₀ᵗ exp(λz) dz。弱吸収・近縮退で exp(λt)−1 が桁落ちしないようにする。</summary>
internal static class BlochWaveIntegral
{
    // 判定は次元を持つ λ ではなく無次元の λt。|λt|≤0.01 なら6次の残余は 2.5e-19 以下。
    internal const double SeriesLimitSquared = 1e-4; // 261106Ch 追加

    /// <summary>261106Ch 追加: t·exprel(λt)。λ=0、t=0 も連続な極限を返す。</summary>
    internal static Complex Evaluate(Complex lambda, double thickness)
    {
        var z = lambda * thickness;
        return (z.Real * z.Real + z.Imaginary * z.Imaginary) <= SeriesLimitSquared
            ? Series(z, thickness)
            : (Complex.Exp(z) - Complex.One) / lambda;
    }

    /// <summary>261106Ch 追加: |z|≤0.01 専用の6次級数。指数・逆数を再利用する縮約経路からも呼ぶ。</summary>
    internal static Complex Series(Complex z, double thickness)
    {
        return thickness * (Complex.One + z * (0.5 + z * (1.0 / 6 + z * (1.0 / 24
            + z * (1.0 / 120 + z * (1.0 / 720 + z / 5040))))));
    }
}
