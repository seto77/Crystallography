using System.Drawing;
using Xunit;

namespace Crystallography.Tests;

/// <summary>261118Cl 新規 (codex2 第 15 回 §4 の 2・第 16 回で写し on を追加、後方源の検証 wf_backsource_remaining_261118 の付随の所見): 六方格子 (三方・六方の結晶) の
/// マスターパターンでブロック共有 (EbsdRepresentativeDirectionBlockSize &gt; 1、製品では格子 512 より上で 3) を使うと、ブロックの中心が
/// 六角形の外の無効なセル (方向 null) に落ちることがある。そのときは画素そのものを代表にして解き、構築が失敗しないこと</summary>
[Collection("Bethe global settings")]
public class EbsdHexBlockShareTests
{
    /// <summary>P3 (空間群 143) の仮の結晶 (Mg を一般の位置に。六方の格子定数)。写しの軌道が小さく (点群 3)、写しの代表の中にもブロックの中心が無効なものが両半球にある
    /// (6/mmm の Mg では、最小の番号の代表は中心が無効な位置に来ない)</summary>
    private static Crystal TrigonalP3()
    {
        int series = 0;
        for (var i = 1; i < SymmetryStatic.TotalSpaceGroupNumber; i++)
            if (SymmetryStatic.Symmetries[i].SpaceGroupNumber == 143) { series = i; break; }
        var cell = (0.32094, 0.32094, 0.52108, Math.PI / 2, Math.PI / 2, 2 * Math.PI / 3);
        Atoms[] atoms = [new("Mg1", 12, 0, 0, null, series, new Vector3DBase(0.1, 0.2, 0.3), 1,
            new DiffuseScatteringFactor(DiffuseScatteringFactor.Type.B, true, 0.005, 0, new double[6], new double[6], cell))];
        return new Crystal(cell, null, series, "P3 (EbsdHexBlockShareTests)", Color.Gray, new Matrix3D(), atoms, ("", "", "", ""), []);
    }

    private static double[] Run(Crystal crystal, Vector3DBase[] dirs, int block, bool symmetry = false)
    {
        // 写し (対称) を切って全部の画素を解く (写しの代表の選び方に依らず、中心が無効な画素そのものを比べるため)。無効なセル (null) も基準の向きで解かれるが比べない
        var bethe = new BetheMethod(crystal) { EbsdRepresentativeDirectionBlockSize = block, UseLocalSurfacePerBeamDirection = false, EbsdUseMasterPatternSymmetry = symmetry };
        using var done = new ManualResetEventSlim();
        Exception error = null;
        bethe.EBSD_Completed += (_, e) => { error = e.Error; done.Set(); };
        bethe.RunEBSDNew(16, [20.0], Matrix3D.IdentityMatrix, [10.0], dirs, BetheMethod.Solver.Eigen_Eigen, 4, false, true, false);
        Assert.True(done.Wait(TimeSpan.FromMinutes(10)), "EBSD did not finish");
        Assert.Null(error);
        return [.. bethe.Disks[0][0].Amplitudes.Take(dirs.Length).Select(a => a.Real * a.Real + a.Imaginary * a.Imaginary)];
    }

    [Fact]
    public void BlockCentreOutsideTheHexagon_UsesThePixelItself()
    {
        var crystal = TrigonalP3();
        Assert.True(MasterPattern.ShouldUseHexGrid(crystal.Symmetry));
        const int gs = 13, block = 3; // N = 6。ブロック 3 の中心 (u, v) = (−2, −5) などが六角形の外 (|u + v| > N) で、その中に有効な画素 (−1, −5) がある
        var dirs = MasterPattern.CreateDirectionsHex(gs, MasterPattern.Hemisphere.PositiveZ)
            .Concat(MasterPattern.CreateDirectionsHex(gs, MasterPattern.Hemisphere.NegativeZ)).ToArray();
        int n = (gs - 1) / 2, half = gs * gs;
        var targets = new List<int>(); // 有効な画素で、ブロックの中心が無効なもの (両半球)
        for (int v = -n; v <= n; v++)
            for (int u = -n; u <= n; u++)
            {
                if (!MasterPattern.IsValidHexCell(u, v, n)) continue;
                int w = u + n, h = v + n;
                int cu = Math.Min(w / block * block + block / 2, gs - 1) - n, cv = Math.Min(h / block * block + block / 2, gs - 1) - n;
                if (!MasterPattern.IsValidHexCell(cu, cv, n)) { int idx = MasterPattern.HexLinearIndex(u, v, n); targets.Add(idx); targets.Add(half + idx); }
            }
        Assert.NotEmpty(targets);

        var shared = Run(crystal, dirs, block);
        var single = Run(crystal, dirs, 1);
        Assert.All(Enumerable.Range(0, dirs.Length).Where(i => dirs[i] != null), i => Assert.True(double.IsFinite(shared[i]) && shared[i] > 0, $"direction #{i}: {shared[i]}"));
        // 中心が無効な画素は自分で解く = ブロック 1 と同じ計算 (旧は基準の向き (表面の法線) で選んだビームの組で解き、大きく違った)
        var worst = targets.Max(i => Math.Abs(shared[i] / single[i] - 1));
        Assert.True(worst < 1e-9, $"{targets.Count} pixels with an invalid block centre; worst relative difference from block 1 = {worst:g3}");

        // 261118Cl 追加 (codex2 第 16 回): 製品の経路 (写し on)。写しの代表 r の中心が無効なら r は自分で解かれ、軌道の全部の像に写る。
        //   両半球にそういう代表があることを確かめ、軌道の全部の像がブロック 1・写し off の r の値と一致することを見る
        var mapping = new BetheMethod(crystal).CreateMasterPatternSymmetryRepresentativeDirectionMapping(dirs.Length, half, gs);
        var reps = targets.Where(i => mapping[i] == i).ToArray();
        Assert.Contains(reps, i => i < half);
        Assert.Contains(reps, i => i >= half);
        var symmetric = Run(crystal, dirs, block, symmetry: true);
        var orbitWorst = Enumerable.Range(0, dirs.Length).Where(i => mapping[i] >= 0 && reps.Contains(mapping[i])).Max(i => Math.Abs(symmetric[i] / single[mapping[i]] - 1));
        Assert.True(orbitWorst < 1e-9, $"{reps.Length} representatives with an invalid block centre (symmetry on); worst relative difference = {orbitWorst:g3}");
    }
}
