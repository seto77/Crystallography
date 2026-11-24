using System.Drawing;
using Xunit;

namespace Crystallography.Tests;

/// <summary>261124Cl 追加 (codex2 査読 C1-03・C1-11、作者の決定 = 二次式で範囲を解く): 逆格子の列挙 (SetPlanes / SetVectorOfG) が、
/// 斜交格子でも範囲内の点を落とさず (総当たりと一致)、高指数でも指数が化けないこと</summary>
public class ReciprocalEnumerationTests
{
    private static Crystal Make(double a, double b, double c, double alpha, double beta, double gamma)
    {
        int series = 1; // P1
        var cell = (a, b, c, alpha, beta, gamma);
        Atoms[] atoms = [new("C1", 6, 0, 0, null, series, new Vector3DBase(0, 0, 0), 1,
            new DiffuseScatteringFactor(DiffuseScatteringFactor.Type.B, true, 0.5 / 100.0, 0, new double[6], new double[6], cell))];
        var crystal = new Crystal(cell, null, series, "P1 (ReciprocalEnumerationTests)", Color.Gray, new Matrix3D(), atoms, ("", "", "", ""), []);
        crystal.SetAxis();
        return crystal;
    }

    /// <summary>総当たり (指数 −n…n) で |G| &lt; gMax の点を数える</summary>
    private static HashSet<(int, int, int)> BruteForce(Crystal c, double gMax, int n)
    {
        var set = new HashSet<(int, int, int)>();
        for (int h = -n; h <= n; h++)
            for (int k = -n; k <= n; k++)
                for (int l = -n; l <= n; l++)
                {
                    if (h == 0 && k == 0 && l == 0) continue;
                    var g = h * c.A_Star + k * c.B_Star + l * c.C_Star;
                    if (g.Length < gMax) set.Add((h, k, l));
                }
        return set;
    }

    [Fact]
    public void SetPlanes_ObliqueLatticeMatchesBruteForce()
    {
        var c = Make(1, 1, 0.5, Math.PI / 2, Math.PI / 2, Math.PI / 6); // γ = 30°
        double gMax = 1.5;
        c.SetPlanes(double.PositiveInfinity, 1 / gMax, false, false, false, false, HorizontalAxis.Angle, 0, 0.1);
        var got = c.Plane.Select(p => (p.h, p.k, p.l)).ToHashSet();
        var expected = BruteForce(c, gMax, 12);
        Assert.True(expected.SetEquals(got), $"列挙 {got.Count} 点、総当たり {expected.Count} 点。抜け: {string.Join(" ", expected.Except(got).Take(5))}");
        Assert.True(got.Contains((1, 1, 0)) || got.Contains((1, -1, 0))); // 旧は短い (110) 系が抜けた
    }

    [Fact]
    public void SetVectorOfG_ObliqueLatticeMatchesBruteForce()
    {
        var c = Make(1, 1, 0.5, Math.PI / 2, Math.PI / 2, Math.PI / 6);
        double gMax = 1.5;
        c.SetVectorOfG(1 / gMax, WaveSource.None);
        var got = c.VectorOfG.Select(v => v.Index).ToHashSet();
        var expected = BruteForce(c, gMax, 12);
        Assert.Equal(expected.Count, c.VectorOfG.Length);
        Assert.True(expected.SetEquals(got));
    }

    [Fact]
    public void SetVectorOfG_HighIndexIsNotCorrupted()
    {
        var c = Make(0.3, 0.3, 300, Math.PI / 2, Math.PI / 2, Math.PI / 2); // l が 255 を越える
        c.SetVectorOfG(1 / 1.2, WaveSource.None, maxNum: 100000);
        Assert.Contains(c.VectorOfG, v => Math.Abs(v.Index.l) > 255);
        foreach (var v in c.VectorOfG)
        {
            var g = v.Index.h * c.A_Star + v.Index.k * c.B_Star + v.Index.l * c.C_Star;
            Assert.True((g - new Vector3DBase(v.X, v.Y, v.Z)).Length < 1e-9, $"指数 {v.Index} と座標が合わない");
        }
    }
}
