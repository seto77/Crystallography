using System.Drawing;
using Xunit;

namespace Crystallography.Tests;

/// <summary>261126Cl 追加 (作者の決定 11-25「±255 キー」・codex2 の 11-25 の議論 観点 F §3): BetheMethod の hkl の鍵を 64 bit にした。
/// 旧の 32 bit の鍵では k・l の欄が 10 bit しか無く、(0,0,257) が (0,0,−767) に復号され、(0,0,512) と (0,1,−512) が同じ鍵になった。
/// 長い軸の P1 で、(1) getU のキャッシュが別の指数の値を返さないこと、(2) Find_gVectors が返すビームの指数と座標が合うこと、を確かめる。
/// 旧い鍵が正しく働く範囲 (|h|,|k|,|l| ≤ 255) では、10,800 条件の選択 (hkl の列・Q・P・U) が旧と同じことを別に確かめた (.project-guidance/reports/codex2_261126_keys)</summary>
public class HklKeyTests
{
    private static Crystal Make(double a, double b, double c)
    {
        int series = 1; // P1
        var cell = (a, b, c, Math.PI / 2, Math.PI / 2, Math.PI / 2);
        Atoms[] atoms = [new("Si1", 14, 0, 0, null, series, new Vector3DBase(0.1, 0.2, 0.3), 1,
            new DiffuseScatteringFactor(DiffuseScatteringFactor.Type.B, true, 0.5 / 100.0, 0, new double[6], new double[6], cell))];
        var crystal = new Crystal(cell, null, series, "P1 (HklKeyTests)", Color.Gray, new Matrix3D(), atoms, ("", "", "", ""), []);
        crystal.SetAxis();
        return crystal;
    }

    private static Vector3DBase G(Crystal c, (int h, int k, int l) i) => i.h * c.A_Star + i.k * c.B_Star + i.l * c.C_Star;

    [Theory]
    [InlineData(0, 0, 512, 0, 1, -512)]   // 旧い鍵で同じ値
    [InlineData(0, 0, 257, 0, 0, -767)]   // 旧い鍵でも鍵は別 (257 は −767 に復号されたが、getU は復号しない)。対照
    [InlineData(0, 1, -512, 0, 0, 512)]   // 逆の順にキャッシュへ入れる
    public void GetUCacheKeepsHighIndicesApart(int h1, int k1, int l1, int h2, int k2, int l2)
    {
        var c = Make(0.31, 0.47, 300);
        var bm = new BetheMethod(c);
        const double kV = 200;
        var first = new BetheMethod.Beam((h1, k1, l1), G(c, (h1, k1, l1)));
        var second = new BetheMethod.Beam((h2, k2, l2), G(c, (h2, k2, l2)));
        var cached1 = bm.getU(kV, first);   // キャッシュへ入る
        var cached2 = bm.getU(kV, second);  // 旧は first の値を返した
        var direct1 = bm.getU(kV, first, cache: false);
        var direct2 = bm.getU(kV, second, cache: false);
        Assert.NotEqual(direct1.Real, direct2.Real); // この 2 つは別の値 (試験の前提)
        Assert.Equal(direct1.Real, cached1.Real);
        Assert.Equal(direct1.Imag, cached1.Imag);
        Assert.Equal(direct2.Real, cached2.Real);
        Assert.Equal(direct2.Imag, cached2.Imag);
    }

    [Fact]
    public void FindGVectorsIndicesMatchTheirVectorsBeyond255()
    {
        var c = Make(0.31, 0.47, 300); // c* が短いので、晶帯軸 [100] の近くで l が 255 を越えるビームが選ばれる
        var bm = new BetheMethod(c);
        const double kV = 200;
        double kvac = UniversalConstants.Convert.EnergyToElectronWaveNumber(kV), u0 = bm.getU(kV).Real.Real;
        var d = c.A_Axis / c.A_Axis.Length;
        var s = new Vector3DBase(-d.X, -d.Y, -d.Z);
        var beams = bm.Find_gVectors(Matrix3D.IdentityMatrix, bm.getVecK0(kvac, u0, d, s), s, 800);
        Assert.Contains(beams, b => b.Index.L >= 257); // 256 は旧い鍵でも正しく復号できたので、257 以上を求める (codex2 の査読 k1)
        Assert.Contains(beams, b => b.Index.L <= -257);
        foreach (var b in beams)
        {
            var g = G(c, (b.Index.H, b.Index.K, b.Index.L));
            Assert.True((g - b.Vec).Length < 1e-9 * Math.Max(1, g.Length), $"指数 {b.Index} と座標が合わない");
        }
    }
}
