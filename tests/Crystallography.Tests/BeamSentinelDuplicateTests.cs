using System.Drawing;
using Xunit;

namespace Crystallography.Tests;

/// <summary>261118Cl 新規 (codex2 第 11 回 §2): Find_gVectors の番兵 (N+1 本目) は「接平面の重複を除いた一意な反射の N+1 本目」。
/// 旧は生の候補の N+1 本を取ってから重複を除いたので、晶帯軸ちょうど (表面の法線が格子の方向に平行) で重複が出ると、残った本物の反射が番兵として落ちた
/// (opx の [100] で 448 条件中 2 回)。opx と同じ胞の単純な直方晶 (Pmmm) の [100]・0.05 kV では、候補は殻を覆い切って N に依らず同じ 364 本
/// (一意な反射は 67 本) なので、大きい N の結果 (番兵なしで一意な反射を全部) を参照にできる。旧は N64・N128・N256 で 53・57・61 本だった</summary>
[Collection("Bethe global settings")]
public class BeamSentinelDuplicateTests
{
    private static Crystal Orthorhombic()
    {
        int series = 0;
        for (var i = 1; i < SymmetryStatic.TotalSpaceGroupNumber; i++)
            if (SymmetryStatic.Symmetries[i].SpaceGroupNumber == 47) { series = i; break; }
        var cell = (1.8228, 0.8805, 0.5185, Math.PI / 2, Math.PI / 2, Math.PI / 2);
        Atoms[] atoms = [new("Mg1", 12, 0, 0, null, series, new Vector3DBase(0.1, 0.2, 0.3), 1,
            new DiffuseScatteringFactor(DiffuseScatteringFactor.Type.B, true, 0.5 / 100.0, 0, new double[6], new double[6], cell))];
        return new Crystal(cell, null, series, "Pmmm (BeamSentinelDuplicateTests)", Color.Gray, new Matrix3D(), atoms, ("", "", "", ""), []);
    }

    private static BetheMethod.Beam[] Select(Crystal crystal, int n)
    {
        const double kv = 0.05;
        var bm = new BetheMethod(crystal);
        var d = crystal.A_Axis / crystal.A_Axis.Length; // [100] (表面の法線が a に平行)
        double kvac = UniversalConstants.Convert.EnergyToElectronWaveNumber(kv), u0 = bm.getU(kv).Real.Real;
        var s = new Vector3DBase(-d.X, -d.Y, -d.Z);
        var beams = bm.Find_gVectors(Matrix3D.IdentityMatrix, bm.getVecK0(kvac, u0, d, s), s, n);
        Assert.True(BetheMethod.LastBeamCandidateScan.ShellCovered && BetheMethod.LastBeamCandidateScan.Collected < 8 * n, "the candidates must cover the Ewald shell (independent of N)");
        return beams;
    }

    [Theory]
    [InlineData(64)]
    [InlineData(128)]
    [InlineData(256)]
    public void ZoneAxisDuplicates_DoNotDropARealReflection(int n)
    {
        var saved = BetheMethod.BeamCandidateCut;
        try
        {
            BetheMethod.BeamCandidateCut = BeamCandidateCutMode.Legacy;
            var crystal = Orthorhombic();
            var all = Select(crystal, 512); // 一意な反射の全部 (番兵なし)、評価値の順
            long before = BetheMethod.SentinelDuplicateCount;
            var beams = Select(crystal, n);
            Assert.True(BetheMethod.SentinelDuplicateCount > before, "the test needs duplicates among the N+1 candidates");
            var got = beams.Select(b => b.Index).ToHashSet();
            var expected = all.Take(n).Select(b => b.Index).ToArray(); // 一意な反射の評価値の下から N 本 (N 以上なら全部)
            var missing = expected.Where(i => !got.Contains(i)).ToArray();
            Assert.True(missing.Length == 0, $"N {n}: {beams.Length} beams, missing {string.Join(" ", missing)}");
            Assert.True(beams.Length >= Math.Min(n, all.Length), $"N {n}: {beams.Length} beams of {all.Length} unique");
            if (n >= all.Length) Assert.Equal(all.Select(b => b.Index).ToHashSet(), got); // 一意な反射が N 以下なら全部
        }
        finally { BetheMethod.BeamCandidateCut = saved; }
    }
}
