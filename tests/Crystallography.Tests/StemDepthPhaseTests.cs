using Xunit;

namespace Crystallography.Tests;

/// <summary>261111Cl 新規: STEM の TDS (と EDX) の深さ方向の位相 (StemDoWork の depthPhase、codex2 第 40 回、作者の判断 2026-11-11)。
/// スライスの振幅から各ビームの搬送波を外し、q ごとに exp(+2πi q_n z)·sinc(π q_n Δz) を掛ける。以前は局所の吸収ポテンシャルの z の位相を欠き、
/// 晶帯軸から傾けた結晶で HAADF が数十 % 外れ、HOLZ のビームが入ると N で収束しなかった。</summary>
public class StemDepthPhaseTests
{
    private const double KV = 200, Conv = 0.015;

    /// <summary>BetheBench と同じ形の入射方向 (収束角の 1.05 倍の正方格子) で STEM を走らせ、TDS の像を返す</summary>
    private static double[][][] RunTds(Crystal crystal, Matrix3D rot, int n, double slice, Vector3D surface = null, double maxThickness = 30, double[] thicknesses = null)
    {
        crystal.RotationMatrix = rot;
        var bethe = new BetheMethod(crystal);
        if (surface != null) bethe.Surface = surface;
        const int div = 8;
        var sin = Math.Sin(Conv * 1.05);
        var directions = new Vector3DBase[div * div];
        for (int h = 0; h < div; h++)
            for (int w = 0; w < div; w++)
            {
                double x = (w - div / 2.0 + 0.5) / (div / 2.0 - 0.5) * sin, y = -(h - div / 2.0 + 0.5) / (div / 2.0 - 0.5) * sin;
                directions[h * div + w] = new Vector3DBase(x, y, -Math.Sqrt(1 - x * x - y * y));
            }
        using var done = new ManualResetEventSlim();
        Exception error = null;
        bethe.StemCompleted += (s, e) => { error = e.Error; done.Set(); };
        var context = SynchronizationContext.Current;
        SynchronizationContext.SetSynchronizationContext(null); // 完了通知をスレッドプールへ
        try { bethe.RunSTEM(n, KV, 5e5, 3, slice, new System.Drawing.Size(12, 12), 0.032, 0, rot, thicknesses ?? [10.0, maxThickness], [0.0], directions, Conv, 0.040, 0.120); }
        finally { SynchronizationContext.SetSynchronizationContext(context); }
        Assert.True(done.Wait(TimeSpan.FromMinutes(5)), "STEM が終わらない");
        Assert.Null(error);
        return bethe.ResultSTEM.ImageTDS;
    }

    /// <summary>max|a−b| / max(|a|, |b|) (BetheBench compare の rel と同じ)</summary>
    private static double MaxRelDiff(double[][][] a, double[][][] b)
    {
        double diff = 0, scale = 0;
        for (int t = 0; t < a.Length; t++)
            for (int d = 0; d < a[t].Length; d++)
                for (int i = 0; i < a[t][d].Length; i++)
                {
                    diff = Math.Max(diff, Math.Abs(a[t][d][i] - b[t][d][i]));
                    scale = Math.Max(scale, Math.Max(Math.Abs(a[t][d][i]), Math.Abs(b[t][d][i])));
                }
        return diff / scale;
    }

    /// <summary>[001] から 5 mrad 傾けた Si: 表面の法線を入射方向のまま (ZOLZ のビームにも g_n・q_n が付く) と、結晶と一緒に傾けた場合 (ZOLZ の g_n = q_n = 0、
    /// 位相の扱いに依らない) の TDS の像が合う。以前の式では約 50 % ずれた (BetheBench、N200)。残りの差は表面の向きの違い (弾性の像も約 2 % 動く)</summary>
    [Fact]
    public void TiltedCrystal_FixedSurface_MatchesSurfaceTiltedWithCrystal()
    {
        var crystal = ConvertCrystalData.ConvertFromCifText(SiCifText).ToCrystal();
        var tilt = Matrix3D.Rot(new Vector3DBase(1, 0, 0), 5e-3);
        var fixedSurface = RunTds(crystal, tilt, 80, 1.0);
        var tiltedSurface = RunTds(crystal, tilt, 80, 1.0, new Vector3D(tilt * new Vector3DBase(0, 0, -1)));
        var rel = MaxRelDiff(fixedSurface, tiltedSurface);
        Assert.True(rel < 0.04, $"表面の向きだけが違う 2 つの計算の TDS の差 {rel:e2} (以前の式では約 0.5)");
    }

    /// <summary>Si [001] の N300 (FOLZ のビームが入り、q_n = 1/a の q の写しが出る) で、スライス 1 nm と 0.5 nm の TDS の像が合う
    /// (q の位相はスライスの中で sinc の重みで積分するので、格子の周期を刻まなくてよい。BetheBench の N500 で 1.4e−4)</summary>
    [Fact]
    public void HolzReplicas_SliceThickness_Converged()
    {
        var crystal = ConvertCrystalData.ConvertFromCifText(SiCifText).ToCrystal();
        var coarse = RunTds(crystal, new Matrix3D(), 300, 1.0, maxThickness: 20);
        var fine = RunTds(crystal, new Matrix3D(), 300, 0.5, maxThickness: 20);
        var rel = MaxRelDiff(coarse, fine);
        Assert.True(rel < 2e-3, $"スライス 1 nm と 0.5 nm の TDS の差 {rel:e2}");
    }

    /// <summary>厚さの列が重なって幅 0 の区間がある (10, 10, 20 nm) 傾けた Si (全部の q に q_n ≠ 0): sinc(0)/0 で NaN にならず、
    /// 2 つ目の厚さの像は 1 つ目と同じ (codex2 第 41 回)</summary>
    [Fact]
    public void ZeroWidthThicknessSegment_NoNaN()
    {
        var crystal = ConvertCrystalData.ConvertFromCifText(SiCifText).ToCrystal();
        var img = RunTds(crystal, Matrix3D.Rot(new Vector3DBase(1, 0, 0), 5e-3), 40, 1.0, thicknesses: [10.0, 10.0, 20.0]);
        Assert.All(img.SelectMany(t => t).SelectMany(d => d), v => Assert.False(double.IsNaN(v), "TDS の像に NaN"));
        var scale = img[0][0].Max(Math.Abs);
        for (int i = 0; i < img[0][0].Length; i++)
            Assert.True(Math.Abs(img[1][0][i] - img[0][0][i]) <= 1e-12 * scale, $"幅 0 の区間で像が変わった (画素 {i})");
    }

    private const string SiCifText = """
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
}
