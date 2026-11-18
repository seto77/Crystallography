using System.Numerics;
using Xunit;
using Xunit.Abstractions;

namespace Crystallography.Tests;

/// <summary>261118Cl 追加: 検査つき・生値の native の入口 (<see cref="NativeWrapper.NumericChecksAvailable"/>) が要る試験。
/// 読み込んだ DLL に入口が無ければ明示的に Skip する。ただし環境変数 CRYSTALLOGRAPHY_REQUIRE_NUMERIC_CHECKS=1 のときは Skip せずに走らせ、
/// 入口が無ければ失敗にする (Skip を合格に数えない。codex2 の単位 8)</summary>
public sealed class NativeNumericFactAttribute : FactAttribute
{
    public NativeNumericFactAttribute()
    {
        if (!NativeWrapper.NumericChecksAvailable && Environment.GetEnvironmentVariable("CRYSTALLOGRAPHY_REQUIRE_NUMERIC_CHECKS") != "1")
            Skip = "Crystallography.Native with the numeric-check entry points is required: " + NativeWrapper.NumericChecksUnavailableReason;
    }
}

/// <summary>261118Cl 新規 (codex2 の単位 6・7 の第 1 段): native の検査つきの EVD・逆行列・LU・CBED と、生値の EBSD の入口の契約。
/// 正常な行列で旧の入口と同じ量を返すこと、非有限の入力・未収束・特異・あふれを状態で返して出力を NaN で埋めること、
/// 旧の EBSD の入口の丸めが生値に旧の意味 (NaN・負・−0・−∞ → +0) を当てたものに等しいことを確かめる</summary>
public class NativeNumericChecksTests
{
    private readonly ITestOutputHelper output;
    public NativeNumericChecksTests(ITestOutputHelper output)
    {
        this.output = output;
        output.WriteLine($"NativeWrapper.Enabled = {NativeWrapper.Enabled} ({NativeWrapper.LoadedNativeLibrary}); NumericChecksAvailable = {NativeWrapper.NumericChecksAvailable}; {NativeWrapper.NumericChecksUnavailableReason}");
    }

    private static Complex[] RandomMatrix(int n, int seed, double imagScale = 0.05)
    {
        var rng = new Random(seed);
        var a = new Complex[n * n];
        for (int i = 0; i < a.Length; i++) a[i] = new Complex(rng.NextDouble() - 0.5, (rng.NextDouble() - 0.5) * imagScale);
        return a;
    }
    private static bool AllNaN(Complex[] a, int n) { for (int i = 0; i < n; i++) if (!double.IsNaN(a[i].Real) || !double.IsNaN(a[i].Imaginary)) return false; return true; }

    /// <summary>managed だけ: probe の判定の本体 (旧い DLL の場合を作って確かめる)</summary>
    [Fact]
    public void Probe_ReportsMissingEntryPointsAndAbiMismatch()
    {
        var notEnabled = NativeWrapper.EvaluateNumericChecksProbe(false, "x", IntPtr.Zero, null, _ => true, () => 1);
        Assert.False(notEnabled.Ok);
        var old = NativeWrapper.EvaluateNumericChecksProbe(true, null, new IntPtr(1), "old.dll", name => !name.EndsWith("Raw"), () => 1);
        Assert.False(old.Ok);
        Assert.Contains("_EBSDSolverRaw", old.Reason);
        var abi = NativeWrapper.EvaluateNumericChecksProbe(true, null, new IntPtr(1), "new.dll", _ => true, () => 2);
        Assert.False(abi.Ok);
        Assert.True(NativeWrapper.EvaluateNumericChecksProbe(true, null, new IntPtr(1), "new.dll", _ => true, () => 1).Ok);
    }

    /// <summary>managed だけ: 旧の丸めの意味</summary>
    [Fact]
    public void LegacyEbsdClamp_MapsNaNNegativeAndMinusZeroToPlusZero()
    {
        Assert.Equal(0.0, NativeWrapper.LegacyEbsdClamp(double.NaN));
        Assert.Equal(0.0, NativeWrapper.LegacyEbsdClamp(-1e-300));
        Assert.Equal(0.0, NativeWrapper.LegacyEbsdClamp(double.NegativeInfinity));
        Assert.False(double.IsNegative(NativeWrapper.LegacyEbsdClamp(-0.0)));
        Assert.Equal(double.PositiveInfinity, NativeWrapper.LegacyEbsdClamp(double.PositiveInfinity));
        Assert.Equal(2.5, NativeWrapper.LegacyEbsdClamp(2.5));
    }

    [NativeNumericFact]
    public void EigenSolverChecked_NormalMatrix_MatchesLegacyAndHasSmallResidual()
    {
        Assert.True(NativeWrapper.NumericChecksAvailable, NativeWrapper.NumericChecksUnavailableReason);
        const int n = 12;
        var a = RandomMatrix(n, 1);
        var (legacyValues, _) = NativeWrapper.EigenSolver(n, a);
        var values = new Complex[n]; var vectors = new Complex[n * n]; var diag = new double[NativeWrapper.NumericDiagLength];
        var st = NativeWrapper.EigenSolverChecked(n, a, values, vectors, NativeWrapper.NumericCheckFlags.VerifyResidual, diag);
        Assert.Equal(NativeWrapper.NativeSolveStatus.Ok, st);
        output.WriteLine($"EVD residual {diag[2]:e2}");
        Assert.True(diag[2] < 1e-13, $"EVD の残差 {diag[2]:e2}");
        var s1 = legacyValues.OrderBy(z => z.Real).ThenBy(z => z.Imaginary).ToArray();
        var s2 = values.OrderBy(z => z.Real).ThenBy(z => z.Imaginary).ToArray();
        for (int i = 0; i < n; i++) Assert.True((s1[i] - s2[i]).Magnitude < 1e-12, $"固有値 {i}: {s1[i]} vs {s2[i]}");
    }

    [NativeNumericFact]
    public void EigenSolverChecked_NonFiniteInputAndNoConvergence_ReturnStatusAndNaN()
    {
        Assert.True(NativeWrapper.NumericChecksAvailable, NativeWrapper.NumericChecksUnavailableReason);
        const int n = 4;
        var a = RandomMatrix(n, 2); a[5] = new Complex(double.NaN, 0);
        var values = new Complex[n]; var vectors = new Complex[n * n];
        Assert.Equal(NativeWrapper.NativeSolveStatus.NonFiniteInput, NativeWrapper.EigenSolverChecked(n, a, values, vectors));
        Assert.True(AllNaN(values, n) && AllNaN(vectors, n * n));

        // 未収束 (試験の旗で反復を 0 回にする): 三角でない 2×2 は未収束、三角 (列優先 [1,0,3,4] = [[1,3],[0,4]]) は反復なしで Success
        var v2 = new Complex[2]; var w2 = new Complex[4];
        Assert.Equal(NativeWrapper.NativeSolveStatus.EvdNoConvergence, NativeWrapper.EigenSolverChecked(2, [1, 2, 3, 4], v2, w2, NativeWrapper.NumericCheckFlags.TestForceNoConvergence));
        Assert.True(AllNaN(v2, 2) && AllNaN(w2, 4));
        Assert.Equal(NativeWrapper.NativeSolveStatus.Ok, NativeWrapper.EigenSolverChecked(2, [1, 0, 3, 4], v2, w2, NativeWrapper.NumericCheckFlags.TestForceNoConvergence));
    }

    [NativeNumericFact]
    public void InverseAndLuChecked_SingularMatrices_ReturnSingular()
    {
        Assert.True(NativeWrapper.NumericChecksAvailable, NativeWrapper.NumericChecksUnavailableReason);
        var inv = new Complex[4]; var x = new Complex[2]; var diag = new double[NativeWrapper.NumericDiagLength];
        foreach (var m in new[] { new Complex[] { 1, 1, 1, 1 }, new Complex[] { new(-0.0, 0), new(-0.0, -0.0), 1, 1 } })
        {
            Assert.Equal(NativeWrapper.NativeSolveStatus.Singular, NativeWrapper.InverseChecked(2, m, inv, diag: diag));
            Assert.True(AllNaN(inv, 4));
            Assert.Equal(NativeWrapper.NativeSolveStatus.Singular, NativeWrapper.PartialPivLuSolveChecked(2, m, [1, 0], x));
            Assert.True(AllNaN(x, 2));
        }
        // 正常な行列は旧の逆行列と同じ (丸めの範囲)
        var a = RandomMatrix(8, 3, 0.3);
        var legacy = NativeWrapper.Inverse(8, a);
        var inv8 = new Complex[64];
        Assert.Equal(NativeWrapper.NativeSolveStatus.Ok, NativeWrapper.InverseChecked(8, a, inv8, NativeWrapper.NumericCheckFlags.VerifyRcond, diag));
        output.WriteLine($"pivot ratio {diag[0]:e2}, rcond {diag[1]:e2}");
        Assert.True(diag[0] > 0 && diag[1] > 0);
        for (int i = 0; i < 64; i++) Assert.True((legacy[i] - inv8[i]).Magnitude <= 1e-12 * Math.Max(1, legacy[i].Magnitude));
    }

    /// <summary>欠陥 (Jordan) の行列: EVD は Ok (縮退や gap では判定しない)、V の LU のピボットの比が小さく記録される</summary>
    [NativeNumericFact]
    public void DefectiveMatrix_EvdOk_PivotRatioRecorded()
    {
        Assert.True(NativeWrapper.NumericChecksAvailable, NativeWrapper.NumericChecksUnavailableReason);
        Complex[] jordan = [1, 0, 1, 1]; // 列優先 [[1,1],[0,1]]
        var v = new Complex[2]; var w = new Complex[4];
        Assert.Equal(NativeWrapper.NativeSolveStatus.Ok, NativeWrapper.EigenSolverChecked(2, jordan, v, w));
        var result = new Complex[2 * 2]; var diag = new double[NativeWrapper.NumericDiagLength];
        var st = NativeWrapper.CBEDSolver_EigenChecked(2, jordan, [1, 0], [1.0, 2.0], result, diag: diag);
        output.WriteLine($"Jordan: status {st}, pivot ratio {diag[0]:e2}");
        Assert.True(st is NativeWrapper.NativeSolveStatus.Ok or NativeWrapper.NativeSolveStatus.Singular, $"{st}");
        if (st == NativeWrapper.NativeSolveStatus.Ok) Assert.True(diag[0] < 1e-6, $"pivot ratio {diag[0]:e2}");
    }

    [NativeNumericFact]
    public void CbedChecked_MatchesLegacyAndFlagsOverflow()
    {
        Assert.True(NativeWrapper.NumericChecksAvailable, NativeWrapper.NumericChecksUnavailableReason);
        const int n = 10;
        var a = RandomMatrix(n, 4, 0.02);
        var psi0 = new Complex[n]; psi0[0] = 1;
        double[] t = [5, 12.5, 40];
        var legacy = NativeWrapper.CBEDSolver_Eigen(a, psi0, t);
        var result = new Complex[n * t.Length];
        Assert.Equal(NativeWrapper.NativeSolveStatus.Ok, NativeWrapper.CBEDSolver_EigenChecked(n, a, psi0, t, result));
        double max = legacy.Max(z => z.Magnitude);
        for (int i = 0; i < result.Length; i++) Assert.True((legacy[i] - result[i]).Magnitude <= 1e-11 * max, $"{i}: {legacy[i]} vs {result[i]}");

        // あふれ: 虚部が大きく負の固有値と厚い t で exp(2πiλt) が +∞ → NonFiniteOutput (計算した値は残す)
        Complex[] grow = [new(1, -50), 0, 0, 2];
        var r2 = new Complex[2];
        Assert.Equal(NativeWrapper.NativeSolveStatus.NonFiniteOutput, NativeWrapper.CBEDSolver_EigenChecked(2, grow, [1, 0], [1000.0], r2));
        Assert.False(double.IsFinite(r2[0].Real) && double.IsFinite(r2[0].Imaginary)); // NaN で埋めずに、計算した非有限の値が残る (行列の積の 0×∞ で他の成分も NaN になりうる)
    }

    /// <summary>261118Cl 追加 (codex2 第 6 回): _CBEDSolver_Eigen2Checked 自体。正常なら Vα = ψ₀ (境界の残差) と厚さ 0 の Tg = ψ₀、強制の未収束なら 4 つの出力が全部 NaN</summary>
    [NativeNumericFact]
    public void Eigen2Checked_BoundaryAndForcedFailure()
    {
        Assert.True(NativeWrapper.NumericChecksAvailable, NativeWrapper.NumericChecksUnavailableReason);
        const int n = 8;
        var a = RandomMatrix(n, 9, 0.03);
        var psi0 = new Complex[n]; psi0[0] = 1;
        double[] t = [0, 15];
        var v = new Complex[n]; var w = new Complex[n * n]; var al = new Complex[n]; var tg = new Complex[n * t.Length]; var diag = new double[NativeWrapper.NumericDiagLength];
        Assert.Equal(NativeWrapper.NativeSolveStatus.Ok, NativeWrapper.CBEDSolver_Eigen2Checked(n, a, psi0, t, v, w, al, tg, NativeWrapper.NumericCheckFlags.VerifyResidual, diag));
        output.WriteLine($"Eigen2: EVD residual {diag[2]:e2}, boundary residual {diag[3]:e2}, pivot ratio {diag[0]:e2}");
        Assert.True(diag[3] < 1e-13, $"境界の残差 {diag[3]:e2}");
        for (int g = 0; g < n; g++) Assert.True((tg[g] - psi0[g]).Magnitude < 1e-12, $"Tg(0)[{g}] = {tg[g]}");
        Assert.Equal(NativeWrapper.NativeSolveStatus.EvdNoConvergence,
            NativeWrapper.CBEDSolver_Eigen2Checked(n, a, psi0, t, v, w, al, tg, NativeWrapper.NumericCheckFlags.TestForceNoConvergence));
        Assert.True(AllNaN(v, n) && AllNaN(w, n * n) && AllNaN(al, n) && AllNaN(tg, n * t.Length));
    }

    /// <summary>261118Cl 追加 (codex2 第 6 回): 生値の入口は native が読む入力の長さを先に確かめる。出力の配列の論理範囲 (tLen) の外の番兵は書き換えない。
    /// diag は書かなかった欄を NaN にする</summary>
    [NativeNumericFact]
    public void RawInputsAreLengthChecked_AndDiagUnwrittenSlotsAreNaN()
    {
        Assert.True(NativeWrapper.NumericChecksAvailable, NativeWrapper.NumericChecksUnavailableReason);
        const int b = 5, nAtoms = 2;
        var values = new Complex[b]; var vectors = new Complex[b * b]; var alpha = new Complex[b]; var phase = new Complex[nAtoms * b];
        for (int i = 0; i < b; i++) { values[i] = new Complex(0.1 * i, -0.001); vectors[i * b + i] = 1; alpha[i] = 0.2; }
        for (int i = 0; i < phase.Length; i++) phase[i] = 1;
        double[] sigma = [0.3, 0.2], thick = [5.0, 10, 15];
        Assert.Throws<ArgumentException>(() => NativeWrapper.EBSDSolverRaw(b, values, new Complex[b * b - 1], alpha, phase, sigma, thick, new double[3]));
        Assert.Throws<ArgumentException>(() => NativeWrapper.EBSDSolverRaw(b, values, vectors, alpha, new Complex[nAtoms * b - 1], sigma, thick, new double[3]));
        Assert.Throws<ArgumentException>(() => NativeWrapper.EBSDSolverWithTDSRaw(b, values, vectors, alpha, phase, sigma, new Complex[b * b - 1], 0.5, thick, new double[3], new double[3]));
        var pooled = new double[8]; Array.Fill(pooled, 12345.0); // プールの配列 (論理の長さ 3 より長い)
        Assert.Equal(NativeWrapper.NativeSolveStatus.Ok, NativeWrapper.EBSDSolverRaw(b, values, vectors, alpha, phase, sigma, thick, pooled));
        for (int i = thick.Length; i < pooled.Length; i++) Assert.Equal(12345.0, pooled[i]);
        // TDS つき・局所の原子 0 (製品で使う形) は通る
        var tds = new double[3]; var loc = new double[3];
        Assert.Equal(NativeWrapper.NativeSolveStatus.Ok, NativeWrapper.EBSDSolverWithTDSRaw(b, values, vectors, alpha, [], [], RandomMatrix(b, 3, 0), 0.5, thick, loc, tds));
        Assert.All(loc, x => Assert.Equal(0.0, x));
        var diag = new double[NativeWrapper.NumericDiagLength]; Array.Fill(diag, 5.0);
        Assert.Equal(NativeWrapper.NativeSolveStatus.Ok, NativeWrapper.InverseChecked(b, vectors, new Complex[b * b], diag: diag));
        Assert.Equal(1.0, diag[0]); // 単位行列のピボットの比
        Assert.True(double.IsNaN(diag[1]) && double.IsNaN(diag[2]) && double.IsNaN(diag[3]), "書かなかった欄は NaN");
    }

    /// <summary>旧の EBSD の入口 = 生値に旧の意味の丸めを当てたもの (ビット単位)。負の σ と NaN の σ で丸めを通す</summary>
    [NativeNumericFact]
    public void EbsdLegacyExport_EqualsLegacyClampOfRaw()
    {
        Assert.True(NativeWrapper.NumericChecksAvailable, NativeWrapper.NumericChecksUnavailableReason);
        const int b = 6, nAtoms = 3;
        var a = RandomMatrix(b, 5, 0.1);
        var values = new Complex[b]; var vectors = new Complex[b * b];
        Assert.Equal(NativeWrapper.NativeSolveStatus.Ok, NativeWrapper.EigenSolverChecked(b, a, values, vectors));
        var psi0 = new Complex[b]; psi0[0] = 1;
        var alpha = new Complex[b];
        Assert.Equal(NativeWrapper.NativeSolveStatus.Ok, NativeWrapper.PartialPivLuSolveChecked(b, vectors, psi0, alpha));
        var rng = new Random(7);
        var phase = new Complex[nAtoms * b];
        for (int i = 0; i < phase.Length; i++) phase[i] = Complex.FromPolarCoordinates(1, rng.NextDouble() * 6.28);
        var mu = RandomMatrix(b, 8, 0.0);
        foreach (var sigma in new[] { new[] { 0.3, 0.2, 0.1 }, new[] { -0.3, -0.2, -0.1 }, new[] { double.NaN, 0.2, 0.1 } })
            foreach (var thick in new[] { new[] { 3.0, 6, 9, 12 }, new[] { 2.0, 7, 20 } }) // 等間隔と非等間隔
            {
                var legacy = NativeWrapper.EBSDSolver(b, values, vectors, alpha, phase, sigma, thick);
                var raw = new double[thick.Length];
                Assert.Equal(NativeWrapper.NativeSolveStatus.Ok, NativeWrapper.EBSDSolverRaw(b, values, vectors, alpha, phase, sigma, thick, raw));
                var (li, lt) = NativeWrapper.EBSDSolverWithTDS(values, vectors, alpha, phase, sigma, mu, 0.7, thick);
                var ri = new double[thick.Length]; var rt = new double[thick.Length];
                Assert.Equal(NativeWrapper.NativeSolveStatus.Ok, NativeWrapper.EBSDSolverWithTDSRaw(b, values, vectors, alpha, phase, sigma, mu, 0.7, thick, ri, rt));
                for (int i = 0; i < thick.Length; i++)
                {
                    Assert.Equal(BitConverter.DoubleToInt64Bits(NativeWrapper.LegacyEbsdClamp(raw[i])), BitConverter.DoubleToInt64Bits(legacy[i]));
                    Assert.Equal(BitConverter.DoubleToInt64Bits(NativeWrapper.LegacyEbsdClamp(ri[i])), BitConverter.DoubleToInt64Bits(li[i]));
                    Assert.Equal(BitConverter.DoubleToInt64Bits(NativeWrapper.LegacyEbsdClamp(rt[i])), BitConverter.DoubleToInt64Bits(lt[i]));
                }
                if (sigma[0] < 0) Assert.True(raw.All(x => x < 0), "負の σ では生値が負 (丸めの前の値が見える)");
                if (double.IsNaN(sigma[0])) Assert.True(raw.All(double.IsNaN), "NaN の σ では生値が NaN");
            }
    }
}
