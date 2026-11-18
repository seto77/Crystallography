using System.Drawing;
using System.Numerics;
using Xunit;

namespace Crystallography.Tests;

/// <summary>261118Cl 新規: 動力学計算 (Bloch 波法) の最小の独立検査。作者の決定 (2026-11-18) で公開の tests に置き、CI で走らせる。
/// codex2 との相談 (第 2 回の単位 1) で、過去の不具合 (P の列割り・EBSD の S の共役・STEM の深さ位相) を見つけたのは、
/// 経路どうしの照合ではなく独立の参照だったことから、それを自動の経路に載せる。
/// <list type="number">
/// <item>流束の保存: 吸収なし (B=0)・表面を 20° 傾けた Si と GaAs で、(1/P₀) Σ_g P_g |ψ_g|² ≡ 1。行の P_g で割る規約の検査
/// (列で割ると、傾いた表面で崩れる。tools/BetheDivCheck と同じ形)。晶帯軸に厳密に平行な表面では両方の割り方が同じ行列になるので、20° 傾ける。</item>
/// <item>RK4 参照: SAED の振幅 (EVD → 入射の境界条件 → 厚さの伝播) を、同じ行列 A から dψ/dz = 2πi A ψ を RK4 で直接積分した強度と比べる。
/// EVD を通らない独立の計算。出射面の位相 exp(πi P_g z) は大きさ 1 なので、強度で比べる。</item>
/// </list></summary>
public class BetheIndependentReferenceTests(Xunit.Abstractions.ITestOutputHelper output)
{
    private const double KV = 200;

    public static TheoryData<string> Crystals => new() { "Si", "GaAs" };

    [Theory]
    [MemberData(nameof(Crystals))]
    public void Flux_IsConserved_WithoutAbsorption_OnTiltedSurface(string name)
    {
        var crystal = name == "Si" ? Silicon(0) : GaAs(0);
        double maxFluxDev = 0, maxAntiDev = 0;
        foreach (var t in new[] { 10.0, 30, 80, 200 })
        {
            var beams = Solve(crystal, 32, t, tiltDeg: 20);
            double p0 = beams[0].P, flux = 0, anti = 0;
            foreach (var b in beams)
            {
                double i = b.Psi.Real * b.Psi.Real + b.Psi.Imaginary * b.Psi.Imaginary;
                Assert.True(double.IsFinite(i), $"非有限の強度 ({b.Index})");
                flux += b.P * i;
                anti += i / b.P;
            }
            maxFluxDev = Math.Max(maxFluxDev, Math.Abs(flux / p0 - 1));
            maxAntiDev = Math.Max(maxAntiDev, Math.Abs(anti * p0 - 1));
        }
        output.WriteLine($"{name}: max|flux−1| = {maxFluxDev:e3}, max|anti−1| = {maxAntiDev:e3}");
        Assert.True(maxFluxDev < 1e-9, $"{name}: max|flux−1| = {maxFluxDev:e2} (行割りなら ≡1)");
        // 検査に判別力があること: 列割りが保存してしまう量 (anti) は、この傾きでは 1 から十分に離れる
        Assert.True(maxAntiDev > 1e-6, $"{name}: max|anti−1| = {maxAntiDev:e2} — 傾きが足りず、列割りと行割りを区別できない");
    }

    [Theory]
    [MemberData(nameof(Crystals))]
    public void Amplitudes_MatchRk4OfTheSameMatrix(string name)
    {
        var crystal = name == "Si" ? Silicon(0.46) : GaAs(0.6); // 吸収あり (U′≠0)
        const double t = 30;
        var bethe = new BetheMethod(crystal);
        bethe.Surface = TiltedSurface(20);
        var beams = bethe.GetDifractedBeamAmpriltudes(24, KV, new Matrix3D(), t);
        int n = beams.Length;
        Assert.True(n >= 10, $"ビームが少なすぎる ({n})");

        // A を規約どおりに組み直す: 非対角 A_gh = [U(h−g) + iU′(h−g)] / P_g、対角 A_gg = (Q_g + iU′₀) / P_g
        var u0Imag = bethe.getU(KV, (0, 0, 0), new Vector3DBase(0, 0, 0)).Imag;
        var a = new Complex[n, n];
        for (int g = 0; g < n; g++)
            for (int h = 0; h < n; h++)
            {
                Complex v;
                if (g == h)
                    v = beams[g].Q + Complex.ImaginaryOne * u0Imag;
                else
                {
                    var (hi, gi) = (beams[h].Index, beams[g].Index);
                    var u = bethe.getU(KV, (hi.H - gi.H, hi.K - gi.K, hi.L - gi.L), beams[h].Vec - beams[g].Vec);
                    v = u.Real + Complex.ImaginaryOne * u.Imag;
                }
                a[g, h] = v / beams[g].P;
            }

        // dψ/dz = 2πi A ψ、ψ(0) = e₀ を RK4 で積分。刻みを半分にして参照自身の収束を確かめる
        var coarse = Rk4(a, t, 4000);
        var fine = Rk4(a, t, 8000);
        double maxI = beams.Max(b => b.Psi.Magnitude * b.Psi.Magnitude);
        double refConv = Enumerable.Range(0, n).Max(i => Math.Abs(coarse[i].Magnitude * coarse[i].Magnitude - fine[i].Magnitude * fine[i].Magnitude)) / maxI;
        Assert.True(refConv < 1e-11, $"{name}: RK4 参照が収束していない (刻み半分との差 {refConv:e2})");

        double maxDiff = Enumerable.Range(0, n).Max(i => Math.Abs(beams[i].Psi.Magnitude * beams[i].Psi.Magnitude - fine[i].Magnitude * fine[i].Magnitude)) / maxI;
        output.WriteLine($"{name}: beams {n}, RK4 self-convergence {refConv:e3}, SAED vs RK4 max rel diff {maxDiff:e3}");
        // 2026-11-18 の実測は Si 9.4e-15・GaAs 4.0e-14 (RK4 参照の自己収束 1e-13)。4 桁の余裕を取って 1e-10
        Assert.True(maxDiff < 1e-10, $"{name}: SAED の強度が RK4 参照と合わない (最大相対差 {maxDiff:e2}、ビーム {n})");
    }

    //---- 補助 ----

    private static BetheMethod.Beam[] Solve(Crystal crystal, int maxBloch, double thickness, double tiltDeg)
    {
        var bethe = new BetheMethod(crystal) { Surface = TiltedSurface(tiltDeg) };
        return bethe.GetDifractedBeamAmpriltudes(maxBloch, KV, new Matrix3D(), thickness);
    }

    private static Vector3D TiltedSurface(double deg)
    {
        double r = deg * Math.PI / 180;
        return new Vector3D(Math.Sin(r), 0, -Math.Cos(r));
    }

    private static Complex[] Rk4(Complex[,] a, double t, int steps)
    {
        int n = a.GetLength(0);
        var psi = new Complex[n];
        psi[0] = Complex.One;
        double h = t / steps;
        var twoPiI = new Complex(0, 2 * Math.PI);
        Complex[] F(Complex[] y)
        {
            var r = new Complex[n];
            for (int i = 0; i < n; i++)
            {
                Complex s = 0;
                for (int j = 0; j < n; j++) s += a[i, j] * y[j];
                r[i] = twoPiI * s;
            }
            return r;
        }
        Complex[] Add(Complex[] y, Complex[] k, double c) => [.. y.Select((v, i) => v + c * k[i])];
        for (int s = 0; s < steps; s++)
        {
            var k1 = F(psi);
            var k2 = F(Add(psi, k1, h / 2));
            var k3 = F(Add(psi, k2, h / 2));
            var k4 = F(Add(psi, k3, h));
            for (int i = 0; i < n; i++) psi[i] += h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
        }
        return psi;
    }

    /// <summary>Si ダイヤモンド構造 (Fd-3m, IT #227, a=5.4309 Å)。B は Å² (0 なら U′=0)。tools/BetheDivCheck と同じ作り</summary>
    private static Crystal Silicon(double bisoAng2)
    {
        var series = SeriesOf(227);
        var cell = (0.54309, 0.54309, 0.54309, Math.PI / 2, Math.PI / 2, Math.PI / 2);
        Atoms[] atoms = [CreateAtom("Si1", 14, 0, 0, 0, bisoAng2, series, cell)];
        return new Crystal(cell, null, series, "Si (BetheIndependentReferenceTests)", Color.Gray, new Matrix3D(),
            atoms, ("", "", "", ""), Bonds.GetVestaBonds(atoms.Select(x => x.AtomicNumber)));
    }

    /// <summary>GaAs 閃亜鉛鉱 (F-43m, IT #216, a=5.6533 Å)。非中心対称の代表</summary>
    private static Crystal GaAs(double bisoAng2)
    {
        var series = SeriesOf(216);
        var cell = (0.56533, 0.56533, 0.56533, Math.PI / 2, Math.PI / 2, Math.PI / 2);
        Atoms[] atoms =
        [
            CreateAtom("Ga1", 31, 0, 0, 0, bisoAng2, series, cell),
            CreateAtom("As1", 33, 0.25, 0.25, 0.25, bisoAng2, series, cell),
        ];
        return new Crystal(cell, null, series, "GaAs (BetheIndependentReferenceTests)", Color.Gray, new Matrix3D(),
            atoms, ("", "", "", ""), Bonds.GetVestaBonds(atoms.Select(x => x.AtomicNumber)));
    }

    private static Atoms CreateAtom(string label, int atomicNumber, double x, double y, double z, double bisoAng2, int series,
        (double A, double B, double C, double Alpha, double Beta, double Gamma) cellNmRad)
        => new(label, atomicNumber, 0, 0, null, series, new Vector3DBase(x, y, z), 1,
            new DiffuseScatteringFactor(DiffuseScatteringFactor.Type.B, true, bisoAng2 / 100.0, 0, new double[6], new double[6], cellNmRad));

    private static int SeriesOf(int itNumber)
    {
        for (var i = 1; i < SymmetryStatic.TotalSpaceGroupNumber; i++)
            if (SymmetryStatic.Symmetries[i].SpaceGroupNumber == itNumber)
                return i;
        throw new ArgumentException($"No space group with IT number {itNumber}.");
    }
}
