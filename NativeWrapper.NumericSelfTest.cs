// 261118Cl 新規 (codex2 第 15 回 §4 の 1「配布物の native の数値検査」): 読み込んだ native DLL そのもので、検査つき・生値の入口が
//   正しく働くかを確かめる自己検査。配布物の起動の検査 (ReciPro.exe --smoke) が呼び、CI の smoke (ARM64・x64 の各 flavor) が結果を要求する。
//   内容は公開の試験 (tests/NativeNumericChecksTests) の要点を小さくしたもの: 入口と ABI、正常な EVD (残差・旧の入口との一致)、
//   生値の EBSD の局所の積分 (閉じた式)、CBED の伝播 (旧の入口との一致)、非有限の入力・強制の未収束・特異・出力のあふれを状態で返すこと
//   (clang-cl の /fp:fast が有限性の判定を畳むと、NaN の入力やあふれを Ok と返す。ARM64 の配布物ではここでしか確かめられない)
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;

namespace Crystallography;

public static partial class NativeWrapper
{
    /// <summary>261118Cl 追加: 読み込んだ native DLL の数値の自己検査。Ok = 全部の項目が合格。Detail は合格なら要約、不合格なら失敗した項目</summary>
    public static (bool Ok, string Detail) NumericSelfTest()
    {
        if (!Enabled) return (false, "native library is not loaded: " + LastLoadError);
        if (!NumericChecksAvailable) return (false, "numeric-check entry points are not available: " + NumericChecksUnavailableReason);
        var failures = new List<string>();
        void Check(bool ok, string what) { if (!ok) failures.Add(what); }
        static bool AllNaN(Complex[] a, int n) { for (int i = 0; i < n; i++) if (!double.IsNaN(a[i].Real) || !double.IsNaN(a[i].Imaginary)) return false; return true; }
        double evdResidual = double.NaN, ebsdRel = double.NaN, cbedRel = double.NaN;
        try
        {
            // 1. 正常な EVD: エルミートに吸収 (対角の正の虚部) を足した 6×6 (深さで減衰する固有値)
            const int b = 6, nAtoms = 3;
            var a = new Complex[b * b];
            for (int i = 0; i < b; i++)
                for (int j = 0; j < b; j++)
                    a[j * b + i] = i == j ? new Complex(0.05 * i, 0.004 * (1 + i)) : new Complex(0.01 / (1 + i + j), 0);
            var values = new Complex[b]; var vectors = new Complex[b * b]; var diag = new double[NumericDiagLength];
            var st = EigenSolverChecked(b, a, values, vectors, NumericCheckFlags.VerifyResidual, diag);
            Check(st == NativeSolveStatus.Ok, $"EVD status {st}");
            evdResidual = diag[2];
            Check(evdResidual < 1e-12, $"EVD residual {evdResidual:e2}");
            var (legacyValues, _) = EigenSolver(b, a);
            var s1 = legacyValues.OrderBy(z => z.Real).ThenBy(z => z.Imaginary).ToArray();
            var s2 = values.OrderBy(z => z.Real).ThenBy(z => z.Imaginary).ToArray();
            Check(Enumerable.Range(0, b).All(i => (s1[i] - s2[i]).Magnitude < 1e-12), "EVD eigenvalues differ from the legacy entry point");

            // 2. 生値の EBSD の局所の源: I(t) = Σ_n σ_n ∫₀ᵗ |Σ_j B_nj e^{2πiγ_j z}|² dz (B = P·C·diag(α)) の閉じた式
            //    ∫₀ᵗ e^{2πi(γ_j − γ_k*) z} dz = (e^{2πi d t} − 1) / (2πi d)、d = γ_j − γ_k*
            var psi0 = new Complex[b]; psi0[0] = 1; var alpha = new Complex[b];
            st = PartialPivLuSolveChecked(b, vectors, psi0, alpha);
            Check(st == NativeSolveStatus.Ok, $"LU status {st}");
            var phase = new Complex[nAtoms * b];
            // for (int i = 0; i < phase.Length; i++) phase[i] = Complex.FromPolarCoordinates(1, 0.7 * i); // 261118Cl 変更前 (codex2 第 16 回: 原子と反射に分離でき、原子の位相は強度で消えた)
            for (int g = 0; g < b; g++) for (int n = 0; n < nAtoms; n++) phase[g * nAtoms + n] = Complex.FromPolarCoordinates(1, 0.7 * g * (n + 1) + 0.3 * n); // 原子ごとに反射に依る位相 (原子の添字の誤りにも感度)
            double[] sigma = [0.4, 0.25, 0.1], thick = [3.0, 8.0];
            var raw = new double[thick.Length];
            st = EBSDSolverRaw(b, values, vectors, alpha, phase, sigma, thick, raw);
            Check(st == NativeSolveStatus.Ok, $"EBSD raw status {st}");
            var B = new Complex[nAtoms, b];
            for (int n = 0; n < nAtoms; n++)
                for (int j = 0; j < b; j++)
                {
                    Complex s = 0;
                    for (int g = 0; g < b; g++) s += phase[g * nAtoms + n] * vectors[j * b + g];
                    B[n, j] = s * alpha[j];
                }
            ebsdRel = 0;
            for (int ti = 0; ti < thick.Length; ti++)
            {
                Complex sum = 0;
                for (int n = 0; n < nAtoms; n++)
                    for (int j = 0; j < b; j++)
                        for (int k = 0; k < b; k++)
                        {
                            var d = values[j] - Complex.Conjugate(values[k]);
                            var f = d.Magnitude < 1e-14 ? thick[ti] : (Complex.Exp(new Complex(0, 2 * Math.PI) * d * thick[ti]) - 1) / (new Complex(0, 2 * Math.PI) * d);
                            sum += sigma[n] * B[n, j] * Complex.Conjugate(B[n, k]) * f;
                        }
                ebsdRel = Math.Max(ebsdRel, Math.Abs(raw[ti] - sum.Real) / Math.Abs(sum.Real));
            }
            Check(ebsdRel < 1e-9, $"EBSD local source differs from the closed form (rel {ebsdRel:e2})");

            // 2b. 旧の EBSD の丸め (NbLegacyClamp): 負・NaN の σ で生値は負・NaN、旧の入口はビットで +0 (fast-math で比較が畳まれると負や NaN が漏れる)
            foreach (var badSigma in new[] { new[] { -0.4, -0.25, -0.1 }, new[] { double.NaN, 0.25, 0.1 } })
            {
                var r = new double[thick.Length];
                st = EBSDSolverRaw(b, values, vectors, alpha, phase, badSigma, thick, r);
                // Check(st == NativeSolveStatus.Ok && r.All(x => double.IsNaN(x) || x < 0), ...); // 261118Cl 変更前 (codex2 第 17 回: 負と NaN を取り違えても、−∞ でも合格した)
                bool rawOk = double.IsNaN(badSigma[0]) ? r.All(double.IsNaN) : r.All(x => double.IsFinite(x) && x < 0); // 負の σ = 全部が有限の負、NaN の σ = 全部が NaN
                Check(st == NativeSolveStatus.Ok && rawOk, $"EBSD raw with σ {badSigma[0]}: status {st}, values {string.Join(", ", r)}");
                var legacyEbsd = EBSDSolver(b, values, vectors, alpha, phase, badSigma, thick);
                Check(legacyEbsd.All(x => BitConverter.DoubleToInt64Bits(x) == 0), $"legacy EBSD clamp with σ {badSigma[0]}: {string.Join(", ", legacyEbsd)} (expected +0)");
            }

            // 3. CBED の伝播: 検査つきの入口が旧の入口と同じ値
            double[] t = [5, 12.5, 40];
            var legacy = CBEDSolver_Eigen(a, psi0, t);
            var result = new Complex[b * t.Length];
            st = CBEDSolver_EigenChecked(b, a, psi0, t, result);
            Check(st == NativeSolveStatus.Ok, $"CBED status {st}");
            Check(legacy.All(z => double.IsFinite(z.Real) && double.IsFinite(z.Imaginary)) && result.All(z => double.IsFinite(z.Real) && double.IsFinite(z.Imaginary)), "CBED values are not finite"); // 261118Cl 追加 (codex2 第 16 回)
            double max = legacy.Max(z => z.Magnitude);
            cbedRel = Enumerable.Range(0, result.Length).Max(i => (legacy[i] - result[i]).Magnitude) / max;
            Check(cbedRel < 1e-11, $"CBED differs from the legacy entry point (rel {cbedRel:e2})");

            // 4. 非有限の入力: 状態 NonFiniteInput、出力は NaN
            var bad = (Complex[])a.Clone(); bad[7] = new Complex(double.NaN, 0);
            st = EigenSolverChecked(b, bad, values, vectors);
            Check(st == NativeSolveStatus.NonFiniteInput && AllNaN(values, b) && AllNaN(vectors, b * b), $"NaN input: status {st}");

            // 5. 強制の未収束 (三角でない 2×2): 状態 EvdNoConvergence
            var v2 = new Complex[2]; var w2 = new Complex[4];
            st = EigenSolverChecked(2, [1, 2, 3, 4], v2, w2, NumericCheckFlags.TestForceNoConvergence);
            Check(st == NativeSolveStatus.EvdNoConvergence && AllNaN(v2, 2), $"forced no-convergence: status {st}");

            // 6. 特異な行列: 状態 Singular
            st = InverseChecked(2, [1, 1, 1, 1], new Complex[4]);
            Check(st == NativeSolveStatus.Singular, $"singular inverse: status {st}");

            // 7. 出力のあふれ: exp(2πiλt) が +∞ → NonFiniteOutput
            st = CBEDSolver_EigenChecked(2, [new(1, -50), 0, 0, 2], [1, 0], [1000.0], new Complex[2]);
            Check(st == NativeSolveStatus.NonFiniteOutput, $"overflow: status {st}");
        }
        catch (Exception ex) { failures.Add($"{ex.GetType().Name}: {ex.Message}"); }
        return failures.Count == 0
            ? (true, $"ok (EVD residual {evdResidual:e1}, EBSD rel {ebsdRel:e1}, CBED rel {cbedRel:e1})")
            : (false, string.Join("; ", failures));
    }
}
