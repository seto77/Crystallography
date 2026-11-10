using System.Threading;
using Xunit;

namespace Crystallography.Tests;

/// <summary>261111Cl 追加: CBED の晶帯軸の対称の利用 (BetheMethod.CbedLegacyUseSymmetry) で、代表の方向から写した振幅が、全方向を解いた振幅と丸めの範囲で一致することを確かめる。</summary>
[Collection("Bethe global settings")]
public class CbedSymmetryTests
{
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
        Assert.Equal("", full.CbedLegacySymmetryDetail);
        Assert.StartsWith("symmetry 2 ", sym.CbedLegacySymmetryDetail);
        Assert.Equal(full.Beams.Select(b => b.Index), sym.Beams.Select(b => b.Index));
        for (int t = 0; t < thick.Length; t++)
            for (int g = 0; g < full.Beams.Length; g++)
            {
                var a = full.Disks[t][g].RawAmplitudes; var b = sym.Disks[t][g].RawAmplitudes;
                for (int r = 0; r < a.Length; r++)
                    Assert.True((a[r] - b[r]).Magnitude <= 1e-10, $"beam {full.Beams[g].Index}, pixel {r}, t {thick[t]}: {a[r]} vs {b[r]}");
            }
    }
}
