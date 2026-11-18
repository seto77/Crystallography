using System.Drawing;
using System.Numerics;
using Xunit;

namespace Crystallography.Tests;

/// <summary>261118Cl 新規: EBSD マスターの後方 TDS 源の軸の直し (<see cref="BetheMethod.MasterBackSourceAxisCorrection"/>) の部品の検査。
/// hkl ごとの 3 点 (入射の軸と G の角 α = 0°・45°・90°) から cos²α の 2 次式で組んだ値が、任意の α で直接評価した値 (同じ求積) と合うこと。
/// 局所の後方源は G → −G で不変なので、α と 180° − α は同じ値になること (軸の正負が効かない)。</summary>
public class EbsdBackSourceAxisTests
{
    [Theory]
    [InlineData(10.0)]
    [InlineData(20.0)]
    [InlineData(30.0)]
    public void ThreeNodeCombination_MatchesDirectEvaluation(double kV)
    {
        var crystal = GaAs(0.6);
        var bethe = new BetheMethod(crystal);
        foreach (var hkl in new[] { (2, 2, 0), (4, 0, 0), (3, 3, 1), (4, 4, 4), (6, 2, 0) })
        {
            double len = (crystal.MatrixInverseTransposed * hkl).Length;
            var nodes = bethe.BackSourceAxisNodes(kV, hkl, len);
            foreach (var deg in new[] { 10.0, 30, 60, 75, 89, 120, 170 })
            {
                double a = deg * Math.PI / 180;
                var direct = Complex.Conjugate(bethe.getU(kV, new BetheMethod.Beam(hkl, new Vector3DBase(len * Math.Sin(a) * Math.Cos(0.4), len * Math.Sin(a) * Math.Sin(0.4), len * Math.Cos(a))),
                    null, Math.PI / 2, Math.PI, 30, 12, cache: false).Imag);
                var combined = BetheMethod.CombineBackSourceAxisNodes(nodes, Math.Cos(a) * Math.Cos(a));
                double rel = (combined - direct).Magnitude / direct.Magnitude;
                // 2026-11-18 の実測 (原子の因子、Si・Ga・As・Au、|G| ≤ 30 nm⁻¹) で ≤ 1e-5 (求積の φ の向きの雑音を含む)。余裕を取って 5e-5
                Assert.True(rel < 5e-5, $"{kV} kV {hkl} α {deg}°: combined vs direct {rel:e2}");
            }
        }
    }

    internal static Crystal GaAs(double bisoAng2) // 261118Cl: EbsdBackSourceBlockAxisTests からも使う (旧 private)
    {
        int series = 0;
        for (var i = 1; i < SymmetryStatic.TotalSpaceGroupNumber; i++)
            if (SymmetryStatic.Symmetries[i].SpaceGroupNumber == 216) { series = i; break; }
        var cell = (0.56533, 0.56533, 0.56533, Math.PI / 2, Math.PI / 2, Math.PI / 2);
        Atoms Make(string label, int z, double x) => new(label, z, 0, 0, null, series, new Vector3DBase(x, x, x), 1,
            new DiffuseScatteringFactor(DiffuseScatteringFactor.Type.B, true, bisoAng2 / 100.0, 0, new double[6], new double[6], cell));
        Atoms[] atoms = [Make("Ga1", 31, 0), Make("As1", 33, 0.25)];
        return new Crystal(cell, null, series, "GaAs (EbsdBackSourceAxisTests)", Color.Gray, new Matrix3D(), atoms, ("", "", "", ""), []);
    }
}
