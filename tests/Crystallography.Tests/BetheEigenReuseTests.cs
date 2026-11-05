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

    /// <summary>261105Cl 追加: 実際に STEM を (小さな像・少ない入射方向で) 走らせ、同じインスタンスで振幅と HRTEM 像を求め直す。
    /// 固有解の使い回し (A1b) に加え、STEM が検出器の角度範囲に限った U' を通常の U と同じ鍵で uDictionary に残す件 (A1) も確かめる。</summary>
    [Fact]
    public void AmplitudesAndHrtem_AfterStemOnSameInstance_MatchFreshInstance()
    {
        var (crystal, rot) = TiltedSi();
        var bethe = new BetheMethod(crystal);
        double[] thicknesses = [Thickness], defocusses = [0, -40];
        var aperture = (R: 0.040, X: 0.0, Y: 0.0);
        var size = new System.Drawing.Size(16, 16);
        double[][][] hrtem(BetheMethod b)
        {
            b.GetHRTEMImage(N, KV, rot, aperture, size, 0.02, 1.0e6, 0.0001, 3, thicknesses, defocusses);
            return b.ResultHRTEM.Image;
        }
        bethe.GetDifractedBeamAmpriltudes(N, KV, rot, Thickness);
        hrtem(bethe);

        // STEM: 3×3 の入射方向 (収束半角 10 mrad)、検出器 40–150 mrad、4×4 画素
        const int div = 3;
        var sin = Math.Sin(0.010);
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
        SynchronizationContext.SetSynchronizationContext(null); // 完了通知をスレッドプールへ (テストのスレッドを待たせても詰まらない)
        try { bethe.RunSTEM(N, KV, 1.0e6, 3, 2, new System.Drawing.Size(4, 4), 0.05, 0, rot, thicknesses, [0.0], directions, 0.010, 0.040, 0.150); }
        finally { SynchronizationContext.SetSynchronizationContext(context); }
        Assert.True(done.Wait(TimeSpan.FromMinutes(5)), "STEM が終わらない");
        Assert.Null(error);

        var after = Amplitudes(bethe.GetDifractedBeamAmpriltudes(N, KV, rot, Thickness));
        var imageAfter = hrtem(bethe);

        var fresh = new BetheMethod(crystal);
        AssertSameAmplitudes(Amplitudes(fresh.GetDifractedBeamAmpriltudes(N, KV, rot, Thickness)), after, "STEM の後");
        var imageFresh = hrtem(fresh);
        var scale = imageFresh.SelectMany(t => t).SelectMany(d => d).Max(Math.Abs);
        for (int t = 0; t < imageFresh.Length; t++)
            for (int d = 0; d < imageFresh[t].Length; d++)
                for (int i = 0; i < imageFresh[t][d].Length; i++)
                    Assert.True(Math.Abs(imageAfter[t][d][i] - imageFresh[t][d][i]) <= 1e-9 * scale, $"HRTEM 像 (厚さ {t}, デフォーカス {d}, 画素 {i})");
    }

    /// <summary>261105Cl 追加: 小さな STEM (3×3 方向、4×4 画素) を走らせ、終わるまで whileRunning を繰り返し呼ぶ。返り値は結果と、実行中に呼べた回数</summary>
    private static (StemSimulationResult Result, int Calls) RunStem(BetheMethod bethe, Matrix3D rot, Action whileRunning = null)
    {
        const int div = 3;
        var sin = Math.Sin(0.010);
        var directions = new Vector3DBase[div * div];
        for (int h = 0; h < div; h++)
            for (int w = 0; w < div; w++)
            {
                double x = (w - div / 2.0 + 0.5) / (div / 2.0 - 0.5) * sin, y = -(h - div / 2.0 + 0.5) / (div / 2.0 - 0.5) * sin;
                directions[h * div + w] = new Vector3DBase(x, y, -Math.Sqrt(1 - x * x - y * y));
            }
        using var done = new ManualResetEventSlim();
        Exception error = null;
        void completed(object s, System.ComponentModel.RunWorkerCompletedEventArgs e) { error = e.Error; done.Set(); }
        bethe.StemCompleted += completed;
        var context = SynchronizationContext.Current;
        SynchronizationContext.SetSynchronizationContext(null); // 完了通知をスレッドプールへ
        int calls = 0;
        try
        {
            bethe.RunSTEM(N, KV, 1.0e6, 3, 2, new System.Drawing.Size(4, 4), 0.05, 0, rot, [Thickness], [0.0], directions, 0.010, 0.040, 0.150);
            while (whileRunning != null && !done.IsSet)
            {
                whileRunning();
                calls++;
            }
            Assert.True(done.Wait(TimeSpan.FromMinutes(5)), "STEM が終わらない");
        }
        finally
        {
            SynchronizationContext.SetSynchronizationContext(context);
            bethe.StemCompleted -= completed;
        }
        Assert.Null(error);
        return (bethe.ResultStem, calls);
    }

    /// <summary>261105Cl 追加: STEM の実行中に、同じインスタンスで別の条件 (電圧・N・方位) の平行ビームの振幅を何度も求めても、
    /// STEM の結果は邪魔の無いときと同じで、実行中に求めた振幅も新しいインスタンスと同じ (STEM は RunSTEM で作った専用の複製で走る)</summary>
    [Fact]
    public void StemAndDiffraction_OnSameInstanceConcurrently_DoNotInterfere()
    {
        var (crystal, rot) = TiltedSi();
        var referenceBethe = new BetheMethod(crystal); // 261105Cl 追加: StemBeams を比べるため
        var reference = RunStem(referenceBethe, rot).Result;

        var bethe = new BetheMethod(crystal);
        var otherRot = Matrix3D.Rot(new Vector3DBase(0.2, 1, 0), 1.5 * Math.PI / 180);
        var settings = new (int N, double KV, Matrix3D Rot)[] { (30, 300, otherRot), (N, KV, rot), (40, 120, rot) };
        var during = new List<((int N, double KV, Matrix3D Rot) Setting, Dictionary<(int H, int K, int L), System.Numerics.Complex> Amp)>();
        var (result, calls) = RunStem(bethe, rot, () =>
        {
            var s = settings[during.Count % settings.Length];
            during.Add((s, Amplitudes(bethe.GetDifractedBeamAmpriltudes(s.N, s.KV, s.Rot, Thickness))));
        });
        Assert.True(calls > 0, "STEM の実行中に回折を計算できていない");

        // STEM の結果は、邪魔の無い実行と同じ
        foreach (var (got, exp, name) in new[] { (result.ImageBoth, reference.ImageBoth, "Both"), (result.ImageEla, reference.ImageEla, "Ela"), (result.ImageTDS, reference.ImageTDS, "TDS") })
        {
            var scale = exp.Planes.SelectMany(t => t).SelectMany(d => d).Max(Math.Abs);
            for (int t = 0; t < exp.Planes.Length; t++)
                for (int d = 0; d < exp.Planes[t].Length; d++)
                    for (int i = 0; i < exp.Planes[t][d].Length; i++)
                        Assert.True(Math.Abs(got.Planes[t][d][i] - exp.Planes[t][d][i]) <= 1e-9 * scale, $"STEM {name} (画素 {i}) が実行中の回折の計算で変わった ({calls} 回)");
        }
        // 261105Cl 追加: STEM で使ったビームは StemBeams で受け取れ、実行中の回折の計算に左右されない
        Assert.NotNull(referenceBethe.StemBeams);
        Assert.Equal(referenceBethe.StemBeams.Select(b => b.Index), bethe.StemBeams.Select(b => b.Index));
        // 実行中に求めた振幅は、新しいインスタンスで求めたものと同じ
        foreach (var (s, amp) in during.Take(settings.Length))
            AssertSameAmplitudes(Amplitudes(new BetheMethod(crystal).GetDifractedBeamAmpriltudes(s.N, s.KV, s.Rot, Thickness)), amp, $"STEM 実行中の回折 (N={s.N}, {s.KV} kV)");
    }
}
