using System.Drawing;
using Xunit;

namespace Crystallography.Tests;

/// <summary>261126Cl 追加 (±255 キーの段 C・比較版): Bethe 用の全列挙 (Crystal.EnumerateReciprocalPoints の betheAllowed つき) が、
/// 斜交・心出し・長い軸の格子で、半径 R の球の中の許される格子点を総当たりと過不足なく返し、座標が hkl と合い、(|g|, h, k, l) の順に並ぶこと。
/// BetheMethod の候補の列挙を C に替える比較 (環境変数 RECIPRO_BEAM_ENUMERATION=C) の土台</summary>
public class QuadraticBeamEnumerationTests
{
    private static Crystal Make(int series, double a, double b, double c, double alpha, double beta, double gamma)
    {
        var cell = (a, b, c, alpha, beta, gamma);
        Atoms[] atoms = [new("C1", 6, 0, 0, null, series, new Vector3DBase(0, 0, 0), 1,
            new DiffuseScatteringFactor(DiffuseScatteringFactor.Type.B, true, 0.5 / 100.0, 0, new double[6], new double[6], cell))];
        var crystal = new Crystal(cell, null, series, "QuadraticBeamEnumerationTests", Color.Gray, new Matrix3D(), atoms, ("", "", "", ""), []);
        crystal.SetAxis();
        return crystal;
    }

    private static int Series(int spaceGroup)
    {
        for (int i = 1; i < SymmetryStatic.TotalSpaceGroupNumber; i++)
            if (SymmetryStatic.Symmetries[i].SpaceGroupNumber == spaceGroup) return i;
        throw new ArgumentException($"space group {spaceGroup}");
    }

    public static TheoryData<string, int, double, double, double, double, double, double, double> Cases => new()
    {
        { "P", 1, 1, 1, 0.5, 90, 90, 30, 1.5 },        // 斜交 P1 (γ = 30°)
        { "P", 1, 0.31, 0.47, 300, 90, 90, 90, 1.0 },  // 長い軸 (l が 255 を越える)
        { "F", 225, 0.4, 0.4, 0.4, 90, 90, 90, 12 },   // F
        { "I", 229, 0.3, 0.3, 0.3, 90, 90, 90, 12 },   // I
        { "RH", 166, 0.5, 0.5, 1.7, 90, 90, 120, 8 },  // R (六方の軸)
        { "C", 15, 0.9, 0.8, 0.5, 90, 105, 90, 6 },    // C (単斜)
    };

    [Theory]
    [MemberData(nameof(Cases))]
    public void MatchesBruteForce(string centering, int spaceGroup, double a, double b, double c, double alphaDeg, double betaDeg, double gammaDeg, double gMax)
    {
        const double d2r = Math.PI / 180;
        var crystal = Make(Series(spaceGroup), a, b, c, alphaDeg * d2r, betaDeg * d2r, gammaDeg * d2r);
        Func<int, int, int, bool> allowed = centering switch
        {
            "F" => (h, k, l) => ((h ^ k) & 1) == 0 && ((k ^ l) & 1) == 0,
            "I" => (h, k, l) => ((h + k + l) & 1) == 0,
            "C" => (h, k, l) => ((h + k) & 1) == 0,
            "RH" => (h, k, l) => ((-h + k + l) % 3 + 3) % 3 == 0,
            _ => (_, _, _) => true,
        };
        var mat = crystal.MatrixInverseTransposed; // g = M (h, k, l)ᵀ (BetheMethod の Build_gCache と同じ、回転なし)
        var inverse = mat.Inverse();
        double ra = inverse.Row1.Length, rb = inverse.Row2.Length, rc = inverse.Row3.Length;
        var got = Crystal.EnumerateReciprocalPoints(mat.E11, mat.E21, mat.E31, mat.E12, mat.E22, mat.E32, mat.E13, mat.E23, mat.E33,
            ra, rb, rc, gMax, halfSpace: false, betheAllowed: allowed);

        // 総当たり (|h| ≤ gMax·|a_i| + 1)
        var expected = new HashSet<(int, int, int)>();
        int hb = (int)Math.Ceiling(gMax * ra) + 1, kb = (int)Math.Ceiling(gMax * rb) + 1, lb = (int)Math.Ceiling(gMax * rc) + 1;
        for (int h = -hb; h <= hb; h++)
            for (int k = -kb; k <= kb; k++)
                for (int l = -lb; l <= lb; l++)
                {
                    if ((h, k, l) == (0, 0, 0) || !allowed(h, k, l)) continue;
                    var g = mat * new Vector3DBase(h, k, l);
                    if (g.Length < gMax * (1 - 1e-12)) expected.Add((h, k, l));
                }
        var gotSet = got.Select(p => (p.h, p.k, p.l)).ToHashSet();
        Assert.Equal(got.Count, gotSet.Count); // 重複なし
        Assert.True(expected.IsSubsetOf(gotSet), $"抜け: {string.Join(" ", expected.Except(gotSet).Take(5))}");
        foreach (var p in got)
        {
            Assert.True(p.len < gMax && allowed(p.h, p.k, p.l));
            var g = mat * new Vector3DBase(p.h, p.k, p.l);
            Assert.True((g - new Vector3DBase(p.x, p.y, p.z)).Length < 1e-9 * Math.Max(1, g.Length), $"{(p.h, p.k, p.l)} の座標");
        }
        for (int i = 1; i < got.Count; i++) // (|g|, h, k, l) の順
            Assert.True((got[i - 1].len, got[i - 1].h, got[i - 1].k, got[i - 1].l).CompareTo((got[i].len, got[i].h, got[i].k, got[i].l)) < 0);
    }
}
