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

    private static int SeriesWithSetting(int spaceGroup, string sub)
    {
        for (int i = 1; i < SymmetryStatic.TotalSpaceGroupNumber; i++)
            if (SymmetryStatic.Symmetries[i].SpaceGroupNumber == spaceGroup && SymmetryStatic.Symmetries[i].SpaceGroupHMsubStr == sub) return i;
        throw new ArgumentException($"space group {spaceGroup} {sub}");
    }

    /// <summary>261126Cl 追加 (codex2 の議論 e2): 一般の回転をかけた実際の C のキャッシュ (BuildQuadraticGCache) を、独立に書いた心出しの判定と
    /// 製品と同じ丸め (BetheMethod.Dot3) の総当たりで照合する。覆い切った半径 Cover 未満の許される点が全部あり、余分も重複も無く、
    /// (|g|, h, k, l) の順で、点の数が K 以上。R は六方の軸 (−h+k+l = 3n) と菱面体の軸 (全整数) を区別する</summary>
    [Theory]
    [InlineData("P", 1, "", 1.0, 1.0, 0.5, 90, 90, 30)]
    [InlineData("F", 227, "", 0.5431, 0.5431, 0.5431, 90, 90, 90)]
    [InlineData("I", 229, "", 0.3, 0.3, 0.3, 90, 90, 90)]
    [InlineData("RH", 166, "H", 0.5, 0.5, 1.7, 90, 90, 120)]
    [InlineData("P", 166, "R", 0.6, 0.6, 0.6, 50, 50, 50)]
    [InlineData("C", 15, "", 0.9, 0.8, 0.5, 90, 105, 90)]
    public void RotatedSnapshotMatchesBruteForce(string centering, int spaceGroup, string setting, double a, double b, double c, double alphaDeg, double betaDeg, double gammaDeg)
    {
        const double d2r = Math.PI / 180;
        int series = setting == "" ? Series(spaceGroup) : SeriesWithSetting(spaceGroup, setting);
        var crystal = Make(series, a, b, c, alphaDeg * d2r, betaDeg * d2r, gammaDeg * d2r);
        Func<int, int, int, bool> allowed = centering switch
        {
            "F" => (h, k, l) => ((h & 1) == (k & 1)) && ((k & 1) == (l & 1)),
            "I" => (h, k, l) => ((h + k + l) % 2 + 2) % 2 == 0,
            "C" => (h, k, l) => ((h + k) % 2 + 2) % 2 == 0,
            "RH" => (h, k, l) => ((-h + k + l) % 3 + 3) % 3 == 0,
            _ => (_, _, _) => true,
        };
        var bm = new BetheMethod(crystal);
        var flags = System.Reflection.BindingFlags.Instance | System.Reflection.BindingFlags.NonPublic;
        var direction = typeof(BetheMethod).GetMethod("GetDirection", flags).Invoke(bm, null);
        var rot = Matrix3D.Rot(new Vector3DBase(0.3, -0.5, 0.81), 0.7);
        var mat = rot * crystal.MatrixInverseTransposed;
        const int K = 3000;
        var snap = typeof(BetheMethod).GetMethod("BuildQuadraticGCache", flags).Invoke(bm, [mat, direction, K]);
        var items = (ValueTuple<long, double, double, double, double>[])snap.GetType().GetProperty("Items").GetValue(snap);
        double cover = (double)snap.GetType().GetProperty("CoverRadius").GetValue(snap);
        Assert.True(items.Length >= K);
        Assert.Equal((0, 0, 0), Decompose(items[0].Item1));

        var got = new HashSet<(int, int, int)>();
        for (int i = 0; i < items.Length; i++)
        {
            var hkl = Decompose(items[i].Item1);
            Assert.True(got.Add(hkl), $"重複 {hkl}");
            Assert.True(allowed(hkl.h, hkl.k, hkl.l), $"心出しで消える点 {hkl}");
            // 261126Cl (codex2 の議論 e3): 座標と長さを試験の中で計算し直して照合する (保存された値に頼らない)
            double gx = BetheMethod.Dot3(mat.E11, hkl.h, mat.E12, hkl.k, mat.E13, hkl.l), gy = BetheMethod.Dot3(mat.E21, hkl.h, mat.E22, hkl.k, mat.E23, hkl.l), gz = BetheMethod.Dot3(mat.E31, hkl.h, mat.E32, hkl.k, mat.E33, hkl.l);
            Assert.Equal((gx, gy, gz, Math.Sqrt(gx * gx + gy * gy + gz * gz)), (items[i].Item2, items[i].Item3, items[i].Item4, items[i].Item5));
            Assert.True(items[i].Item5 < cover, $"{hkl} は Cover の外");
            if (i > 0)
                Assert.True((items[i - 1].Item5, Decompose(items[i - 1].Item1)).CompareTo((items[i].Item5, hkl)) < 0, "順");
        }
        // K 番目を含む殻で止めている (その手前の殻の終わりは K 未満): 殻 = 先頭の長さ s から s(1+1e-9) まで
        int lastShellStart = items.Length - 1;
        while (lastShellStart > 1 && items[lastShellStart - 1].Item5 >= items[lastShellStart].Item5 / (1 + 1E-9) * (1 - 1E-15)) lastShellStart--;
        Assert.True(lastShellStart < K, $"K 番目より先の殻まで採っている (最後の殻の始まり {lastShellStart}、K {K})");
        var inverse = mat.Inverse();
        var expected = new HashSet<(int, int, int)>(); // 261126Cl 追加
        int hb = (int)Math.Ceiling(cover * inverse.Row1.Length) + 1, kb = (int)Math.Ceiling(cover * inverse.Row2.Length) + 1, lb = (int)Math.Ceiling(cover * inverse.Row3.Length) + 1;
        for (int h = -hb; h <= hb; h++)
            for (int k = -kb; k <= kb; k++)
                for (int l = -lb; l <= lb; l++)
                {
                    if (!allowed(h, k, l)) continue;
                    double gx = BetheMethod.Dot3(mat.E11, h, mat.E12, k, mat.E13, l), gy = BetheMethod.Dot3(mat.E21, h, mat.E22, k, mat.E23, l), gz = BetheMethod.Dot3(mat.E31, h, mat.E32, k, mat.E33, l);
                    if (Math.Sqrt(gx * gx + gy * gy + gz * gz) < cover)
                        Assert.True(expected.Add((h, k, l)) && got.Contains((h, k, l)), $"抜け {(h, k, l)}");
                }
        Assert.True(expected.SetEquals(got), $"余分 {string.Join(" ", got.Except(expected).Take(5))}"); // 261126Cl: 余分も独立に
        if (spaceGroup == 166) // 六方の軸では (1,0,0) は消え、菱面体の軸では残る
            Assert.Equal(setting == "R", got.Contains((1, 0, 0)));
    }

    private static (int h, int k, int l) Decompose(long key)
    {
        const int bias = 1 << 20; const long mask = (1L << 21) - 1;
        return ((int)(key >> 42) - bias, (int)((key >> 21) & mask) - bias, (int)(key & mask) - bias);
    }
}
