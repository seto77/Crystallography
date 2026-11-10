using System.Threading;
using Xunit;

namespace Crystallography.Tests;

/// <summary>261111Cl 追加: CBED の晶帯軸の対称の利用 (BetheMethod.CbedLegacyUseSymmetry) で、代表の方向から写した振幅が、全方向を解いた振幅と丸めの範囲で一致することを確かめる。</summary>
[Collection("Bethe global settings")]
public class CbedSymmetryTests
{
    private const string GaAsCif = """
data_GaAs
_chemical_name_common GaAs
_cell_length_a 5.6533
_cell_length_b 5.6533
_cell_length_c 5.6533
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'F -4 3 m'
_symmetry_Int_Tables_number 216
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_B_iso_or_equiv
Ga1 Ga 0 0 0 0.6
As1 As 0.25 0.25 0.25 0.5
""";

    /// <summary>Si (Fd−3m の原点の選択 1、原子が 0,0,0)。[110] の操作の並進の位相が ±i になり、写しの位相の符号の誤りが振幅に出る (原点の選択 2 では ±1 で見えない)</summary>
    private const string SiOrigin1Cif = """
data_Si
_chemical_name_common Si
_cell_length_a 5.4309
_cell_length_b 5.4309
_cell_length_c 5.4309
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'F d -3 m :1'
_symmetry_Int_Tables_number 227
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_B_iso_or_equiv
Si1 Si 0 0 0 0.46
""";

    /// <summary>晶帯軸 [uvw] を実験室の +z へ向ける回転 (CBED の入射は −z)</summary>
    private static Matrix3D ZoneToZ(Crystal crystal, int u, int v, int w)
    {
        var axis = u * crystal.A_Axis + v * crystal.B_Axis + w * crystal.C_Axis; axis /= axis.Length;
        var z = new Vector3DBase(0, 0, 1);
        var rot = Vector3DBase.VectorProduct(axis, z);
        double th = Math.Acos(Math.Clamp(axis * z, -1, 1));
        var r1 = Matrix3D.Rot(rot / rot.Length, th);
        return ((r1 * axis) - z).Length < 1e-9 ? r1 : Matrix3D.Rot(rot / rot.Length, -th);
    }

    /// <summary>画面と同じ方向の正方格子 (一辺 division、半角 alpha [rad])</summary>
    private static Vector3DBase[] Grid(int division, double alpha)
    {
        var dirs = new Vector3DBase[division * division];
        double radius = division / 2.0, sin = Math.Sin(alpha);
        for (int h = 0; h < division; h++)
            for (int w = 0; w < division; w++)
            {
                double x = (w - radius + 0.5) / (radius - 0.5) * sin, y = -(h - radius + 0.5) / (radius - 0.5) * sin;
                dirs[h * division + w] = new Vector3DBase(x, y, -Math.Sqrt(1 - x * x - y * y));
            }
        return dirs;
    }

    private static BetheMethod RunLegacy(Crystal crystal, Matrix3D rot, Vector3DBase[] dirs, double[] thick, bool symmetry, int n)
    {
        var bethe = new BetheMethod(crystal) { CbedLegacyUseSymmetry = symmetry };
        using var done = new ManualResetEventSlim();
        Exception error = null;
        bethe.CBED_Completed += (s, e) => { error = e.Error; done.Set(); };
        bethe.RunCBED(n, 200, rot, thick, dirs, false, BetheMethod.Solver.Auto, 4);
        Assert.True(done.Wait(TimeSpan.FromMinutes(5)), "CBED did not finish");
        Assert.Null(error);
        return bethe;
    }

    /// <summary>Si [110] 200 kV (2mm、鏡映は画面の斜め): 対称を使うと 4 分の 1 近くの方向だけを解き、振幅 (位相込み) は全方向を解いた値と 1e−10 で一致する。
    /// N200 は FOLZ のビーム (|g| ≈ 44 nm⁻¹) を含み、並進のビーム軸方向の成分による HOLZ の層の位相を確かめる (照合の組がビームの組の ZOLZ だけだった頃は e^{−i 105°} ずれた)</summary>
    [Theory]
    [InlineData(32, 10)]
    [InlineData(200, 6)]
    public void LegacyCbed_WithSymmetry_MatchesFullSolve(int n, int division)
    {
        var (crystal, _) = BetheEigenReuseTests.TiltedSi(); // 結晶だけ使う (向きは下の rot を RunCBED に渡す)
        var rot = ZoneToZ(crystal, 1, 1, 0);
        var dirs = Grid(division, 1e-3);
        double[] thick = [30, 60];
        var full = RunLegacy(crystal, rot, dirs, thick, false, n);
        var sym = RunLegacy(crystal, rot, dirs, thick, true, n);
        if (n >= 200) Assert.Contains(full.Beams, b => Math.Abs(b.Vec.Z) > 1); // HOLZ のビームが組にある
        Assert.Equal("", full.CbedSymmetryDetail);
        Assert.StartsWith("symmetry 2 ", sym.CbedSymmetryDetail);
        Assert.Equal(full.Beams.Select(b => b.Index), sym.Beams.Select(b => b.Index));
        for (int t = 0; t < thick.Length; t++)
            for (int g = 0; g < full.Beams.Length; g++)
            {
                var a = full.Disks[t][g].RawAmplitudes; var b = sym.Disks[t][g].RawAmplitudes;
                for (int r = 0; r < a.Length; r++)
                    Assert.True((a[r] - b[r]).Magnitude <= 1e-10, $"beam {full.Beams[g].Index}, pixel {r}, t {thick[t]}: {a[r]} vs {b[r]}");
            }
    }

    /// <summary>対称心のない GaAs [110] (−43m): ビーム軸まわりの 2 回回転は ZOLZ の |U| が Friedel で等しくても対称ではないので使わず、鏡映 (1−10) だけを使う。写した振幅は全方向を解いた値と一致する</summary>
    [Fact]
    public void LegacyCbed_GaAs110_UsesOnlyTheMirror()
    {
        var crystal = ConvertCrystalData.ConvertFromCifText(GaAsCif).ToCrystal();
        var rot = ZoneToZ(crystal, 1, 1, 0);
        var dirs = Grid(8, 1e-3);
        double[] thick = [40];
        var full = RunLegacy(crystal, rot, dirs, thick, false, 32);
        var sym = RunLegacy(crystal, rot, dirs, thick, true, 32);
        Assert.StartsWith("symmetry m(", sym.CbedSymmetryDetail);
        Assert.DoesNotContain("symmetry 2", sym.CbedSymmetryDetail);
        Assert.DoesNotContain(" 2 ", sym.CbedSymmetryDetail);
        for (int g = 0; g < full.Beams.Length; g++)
        {
            var a = full.Disks[0][g].RawAmplitudes; var b = sym.Disks[0][g].RawAmplitudes;
            for (int r = 0; r < a.Length; r++) Assert.True((a[r] - b[r]).Magnitude <= 1e-10, $"beam {full.Beams[g].Index}, pixel {r}");
        }
    }

    /// <summary>晶帯軸から外れた向き (Si、[001] から 0.5° 傾けた) では対称を使わない</summary>
    [Fact]
    public void LegacyCbed_OffZone_DoesNotUseSymmetry()
    {
        var (crystal, rot) = BetheEigenReuseTests.TiltedSi();
        var sym = RunLegacy(crystal, rot, Grid(6, 1e-3), [30], true, 32);
        Assert.Equal("", sym.CbedSymmetryDetail);
    }

    /// <summary>新方式 (着衣の裾) の対称の写し: 各方向を自分の区画で解く設定 (区画 1e−6 mrad) なら、写した振幅は全方向を解いた振幅と単精度の裾の積の丸め (1e−6) の範囲で一致する。
    /// 軽い設定 (N60・候補の半径 20 nm⁻¹・裾は全部) で、Si [110] 200 kV は 2 回回転と 2 つの鏡映を使う</summary>
    [Fact]
    public void DressedTailCbed_WithSymmetry_MatchesFullSolve()
    {
        var crystal = ConvertCrystalData.ConvertFromCifText(SiOrigin1Cif).ToCrystal();
        Assert.Equal(8, crystal.Atoms.Sum(a => a.Atom.Count())); // 単位胞の 8 原子 (原点の選択 1 が読めている)
        var rot = ZoneToZ(crystal, 1, 1, 0);
        var dirs = Grid(6, 0.3e-3);
        double[] thick = [40, 80];
        BetheMethod Run(bool sym)
        {
            var bethe = new BetheMethod(crystal) { CbedBeamSelection = BetheMethod.CbedBeamSelectionMethod.DressedTail };
            bethe.CbedDressedTail = BetheMethod.CbedDressedTailDefaults with { N = 60, PoolRadiusNm = 20, CoreRadiusNm = double.NaN, NeighborRadiusNm = double.NaN, PatchRadiusMrad = 1e-6, UseSymmetry = sym };
            using var done = new ManualResetEventSlim();
            Exception error = null;
            bethe.CBED_Completed += (s, e) => { error = e.Error; done.Set(); };
            bethe.RunCBED(64, 200, rot, thick, dirs, false, BetheMethod.Solver.Auto, 4);
            Assert.True(done.Wait(TimeSpan.FromMinutes(10)), "CBED did not finish");
            Assert.Null(error);
            Assert.Equal(BetheMethod.CbedBeamSelectionMethod.DressedTail, bethe.CbedBeamSelectionUsed);
            return bethe;
        }
        var full = Run(false);
        var sym = Run(true);
        Assert.Equal("", full.CbedSymmetryDetail);
        Assert.StartsWith("symmetry 2 ", sym.CbedSymmetryDetail);
        var idx = sym.Beams.Select((b, i) => (b.Index, i)).ToDictionary(x => x.Index, x => x.i);
        for (int t = 0; t < thick.Length; t++)
            for (int g = 0; g < full.Beams.Length; g++)
            {
                Assert.True(idx.TryGetValue(full.Beams[g].Index, out var g1), $"beam {full.Beams[g].Index} missing");
                var a = full.Disks[t][g].RawAmplitudes; var b = sym.Disks[t][g1].RawAmplitudes;
                for (int r = 0; r < a.Length; r++)
                    Assert.True((a[r] - b[r]).Magnitude <= 1e-5, $"beam {full.Beams[g].Index}, pixel {r}, t {thick[t]}: {a[r]} vs {b[r]}");
            }
    }
}
