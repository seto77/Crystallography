using Xunit;
using static Crystallography.Tests.SymmetryReference;

namespace Crystallography.Tests;

/// <summary>261124Cl 追加: SymmetryOperation.ConvertPlaneIndex が、各操作の線形部 R (ApplyMatrix) の転置 Rᵀ·h を返すことを、
/// 操作 1 個ずつ確かめる。異方性 DW の hᵀ(RβRᵀ)h = (Rᵀh)ᵀβ(Rᵀh) はこの対応に依存する。
/// 等価な面の「集合」の比較では、二つの操作の結果が入れ替わっていても集合は一致してしまうので見逃す
/// (六方軸の [210] と [120] の二回軸が入れ替わっていた: codex2 第 1 回 C1-01)。</summary>
public class PlaneIndexTransposeTests
{
    private static readonly (int H, int K, int L)[] Probes = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 2, 3), (-2, 1, 4)];

    [Theory, MemberData(nameof(AllSeries), MemberType = typeof(SymmetryReference))]
    public void ConvertPlaneIndex_IsTransposeOfApplyMatrix(int series)
    {
        foreach (var so0 in SymmetryStatic.WyckoffPositions[series][0].PositionOperations)
        {
            var so = new SymmetryOperation(so0, series);
            var (m, _) = Decompose(so.ApplyMatrix);
            foreach (var (h, k, l) in Probes)
            {
                // (Rᵀh)_j = Σ_i R_ij h_i  (m は行優先: m[3i+j] = R_ij)
                var expected = (H: (int)Math.Round(m[0] * h + m[3] * k + m[6] * l),
                                K: (int)Math.Round(m[1] * h + m[4] * k + m[7] * l),
                                L: (int)Math.Round(m[2] * h + m[5] * k + m[8] * l));
                var actual = so.ConvertPlaneIndex(h, k, l);
                Assert.True(expected == actual,
                    $"series {series}, order {so.Order}, axis {so.Direction}: Rᵀ({h},{k},{l}) = {expected} だが ConvertPlaneIndex = {actual}");

                var actualD = so.ConvertPlaneIndex((double)h, k, l);
                Assert.True(Math.Abs(actualD.H - expected.H) < 1e-12 && Math.Abs(actualD.K - expected.K) < 1e-12 && Math.Abs(actualD.L - expected.L) < 1e-12,
                    $"series {series}, order {so.Order}, axis {so.Direction} (実数版): Rᵀ({h},{k},{l}) = {expected} だが {actualD}");
            }
        }
    }
}
