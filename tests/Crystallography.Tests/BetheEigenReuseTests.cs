using System.Numerics;
using Xunit;

namespace Crystallography.Tests;

/// <summary>261105Cl 新規: 同じ BetheMethod (ReciPro では結晶ごとに 1 つの Crystal.Bethe を回折・像シミュレータで共有) で、
/// Beams が置き換えられたあと (STEM・PED など) に平行ビームの振幅を求めても、新しいインスタンスで求めたものと一致すること。
/// 以前は N・電圧・方位が同じだと古い固有解を使い回し、振幅を並び順の違う Beams に入れていた (implementation-questions A1b)。</summary>
public class BetheEigenReuseTests
{
    private const string SiCif = """
data_Si
_chemical_name_common Si
_cell_length_a 5.4309
_cell_length_b 5.4309
_cell_length_c 5.4309
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'F d -3 m :2'
_symmetry_Int_Tables_number 227
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_B_iso_or_equiv
Si1 Si 0.125 0.125 0.125 0.46
""";

    private const int N = 50;
    private const double KV = 200, Thickness = 30;

    /// <summary>Si、[001] から 0.5° 傾けた方位 (晶帯軸から外すと、Find_gVectors の順と |g|² 順がほとんど一致しない)</summary>
    internal static (Crystal Crystal, Matrix3D Rotation) TiltedSi()
    {
        var crystal = ConvertCrystalData.ConvertFromCifText(SiCif).ToCrystal();
        var rot = Matrix3D.Rot(new Vector3DBase(1, 0.3, 0), 0.5 * Math.PI / 180);
        crystal.RotationMatrix = rot;
        return (crystal, rot);
    }

    internal static Dictionary<(int H, int K, int L), Complex> Amplitudes(BetheMethod.Beam[] beams) => beams.ToDictionary(b => b.Index, b => b.Psi);

    internal static void AssertSameAmplitudes(Dictionary<(int H, int K, int L), Complex> expected, Dictionary<(int H, int K, int L), Complex> actual, string label)
    {
        Assert.Equal(expected.Count, actual.Count);
        foreach (var (index, psi) in expected)
            Assert.True((psi - actual[index]).Magnitude <= 1e-10 * Math.Max(1, psi.Magnitude),
                $"{label}: {index} |ψ|² {actual[index].Magnitude * actual[index].Magnitude:g6}, 正しくは {psi.Magnitude * psi.Magnitude:g6}");
    }

    [Fact]
    public void Amplitudes_AfterBeamsReplacedByReorderedArray_MatchFreshInstance()
    {
        var (crystal, rot) = TiltedSi();
        var bethe = new BetheMethod(crystal);
        bethe.GetDifractedBeamAmpriltudes(N, KV, rot, Thickness);
        var before = bethe.Beams.Select(b => b.Index).ToArray();

        // StemDoWork と同じ置き換え: 同じ反射を |g|² 順に並べた別の配列 (固有解は消さない)
        bethe.Beams = [.. bethe.Beams.OrderBy(b => b.Vec.Length2)];
        Assert.True(before.Zip(bethe.Beams).Count(p => p.First != p.Second.Index) > N / 2, "並び順がほとんど変わる方位で確かめる");

        var after = Amplitudes(bethe.GetDifractedBeamAmpriltudes(N, KV, rot, Thickness));
        var fresh = Amplitudes(new BetheMethod(crystal).GetDifractedBeamAmpriltudes(N, KV, rot, Thickness));
        AssertSameAmplitudes(fresh, after, "Beams を並べ替えた後");
    }
}
