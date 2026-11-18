// 261118Cl 新規作成 (codex2 の単位 6・7 の第 1 段): 数値の失敗を C# へ伝える native の入口 (Crystallography.Native/EigenFuncs.cpp の
//   _NumericChecksAbiVersion・…Checked・…Raw、ABI の版 1)。既存の入口 (_EigenSolver・_Inverse・_PartialPivLuSolve・_CBEDSolver_Eigen・
//   _CBEDSolver_Eigen2・_EBSDSolver・_EBSDSolverWithTDS) は変えずに残り、旧い DLL でもそちらは動く。新しい入口は、読み込んだ DLL に
//   8 本すべてがあり ABI の版が合うときだけ呼ぶ (<see cref="NumericChecksAvailable"/>。LibraryImport は最初の呼び出しで結び付くので、
//   確かめてから呼べば旧い DLL で EntryPointNotFoundException は出ない)。
//   状態の意味 (261118Cl codex2 第 6 回で明記): Checked の数値の失敗 (NonFiniteInput・EvdNoConvergence・NonFiniteEigen・Singular・NonFiniteSolve・
//   OutOfMemory・Unknown) では、native は有効な出力の領域を quiet NaN で埋める (前の方向の解やごみが、もっともらしい値として残らない)。
//   NonFiniteOutput では調べるために計算した値を残す。Args は検査の直後に戻るので出力は未規定。Raw は値を検査せず、OutOfMemory・Unknown でも埋めない。
//   呼び出し側は状態を最初に読み、Ok でない出力を次の段へ渡さないこと。diag は書いた欄だけが有効 (包みは呼ぶ前に NaN で初期化する)。
//   Checked の EVD は出力に触る前に ComplexEigenSolver.info() を見る (旧の _EigenSolver は未収束のとき未初期化の領域を写していた)。
//   LU は同梱の Eigen の PartialPivLU.info() が常に Success なので、U の対角が厳密に 0 か非有限なら Singular を返す。
using System;
using System.Numerics;
using System.Runtime.InteropServices;

namespace Crystallography;

public static partial class NativeWrapper
{
    /// <summary>261118Cl 追加: 検査つきの native の入口の戻り値 (EigenFuncs.h の EigenFuncsStatus。5・6・9 は DressedTailStatus と同じ値)。
    /// <see cref="Unchecked"/> は managed だけの値 (旧い DLL などで検査せずに通った)</summary>
    public enum NativeSolveStatus
    {
        Ok = 0, NonFiniteInput = 1, EvdNoConvergence = 2, NonFiniteEigen = 3, Singular = 4, Args = 5, OutOfMemory = 6,
        NonFiniteSolve = 7, NonFiniteOutput = 8, Unknown = 9, Unchecked = 100,
    }

    /// <summary>261118Cl 追加: 検査の旗 (EigenFuncs.h の EigenFuncsFlags)</summary>
    [Flags]
    public enum NumericCheckFlags
    {
        None = 0,
        /// <summary>diag[1] = rcond (O(N²))</summary>
        VerifyRcond = 1,
        /// <summary>diag[2] = EVD の残差 ‖AV−VΛ‖_F/(‖A‖_F‖V‖_F)、diag[3] = 境界の残差 ‖Vα−ψ₀‖/‖ψ₀‖ (O(N³)、検証モードだけ)</summary>
        VerifyResidual = 2,
        /// <summary>試験専用: setMaxIterations(0) で未収束を起こす</summary>
        TestForceNoConvergence = 4,
    }

    /// <summary>261118Cl 追加: diag の長さ ([0] ピボットの比・[1] rcond・[2] EVD の残差・[3] 境界の残差)</summary>
    public const int NumericDiagLength = 4;
    internal const int NumericChecksAbiVersionExpected = 1;

    [LibraryImport("Crystallography.Native.dll")] private static partial int _NumericChecksAbiVersion();
    [LibraryImport("Crystallography.Native.dll")]
    private static unsafe partial int _EigenSolverChecked(int dim, double* mat, double* eigenValues, double* eigenVectors, int flags, double* diag, int diagLength);
    [LibraryImport("Crystallography.Native.dll")]
    private static unsafe partial int _InverseChecked(int dim, double* mat, double* inverse, int flags, double* diag, int diagLength);
    [LibraryImport("Crystallography.Native.dll")]
    private static unsafe partial int _PartialPivLuSolveChecked(int dim, double* mat, double* vec, double* result, int flags, double* diag, int diagLength);
    [LibraryImport("Crystallography.Native.dll")]
    private static unsafe partial int _CBEDSolver_EigenChecked(int dim, double* potential, double* psi0, int tDim, double* thickness, double* result, int flags, double* diag, int diagLength);
    [LibraryImport("Crystallography.Native.dll")]
    private static unsafe partial int _CBEDSolver_Eigen2Checked(int dim, double* potential, double* psi0, int tDim, double* thickness,
        double* values, double* vectors, double* alphas, double* tg, int flags, double* diag, int diagLength);
    [LibraryImport("Crystallography.Native.dll")]
    private static unsafe partial int _EBSDSolverRaw(int bLen, int nAtoms, int tLen, double* eigenValues, double* eigenVectors, double* alpha,
        double* phaseNG, double* sigma, double* thicknesses, double* intensity);
    [LibraryImport("Crystallography.Native.dll")]
    private static unsafe partial int _EBSDSolverWithTDSRaw(int bLen, int nAtoms, int tLen, double* eigenValues, double* eigenVectors, double* alpha,
        double* phaseNG, double* sigma, double* muBack, double tdsCoeff, double* thicknesses, double* intensity, double* tdsIntensity);

    /// <summary>261118Cl 追加: 検査つき・生値の入口の名前 (probe で全部の有無を確かめる)</summary>
    internal static readonly string[] NumericCheckExports =
        ["_NumericChecksAbiVersion", "_EigenSolverChecked", "_InverseChecked", "_PartialPivLuSolveChecked",
         "_CBEDSolver_EigenChecked", "_CBEDSolver_Eigen2Checked", "_EBSDSolverRaw", "_EBSDSolverWithTDSRaw"];

    /// <summary>261118Cl 追加: probe の判定の本体 (純関数。CI の試験から直接呼べる)。missing = 欠けた入口の名前、abi = 読めた ABI の版 (読めなければ null)</summary>
    internal static (bool Ok, string Reason) EvaluateNumericChecksProbe(bool enabled, string lastLoadError, IntPtr handle, string library,
        Func<string, bool> hasExport, Func<int> abiVersion)
    {
        if (!enabled) return (false, "native library is not enabled" + (lastLoadError != null ? $" ({lastLoadError})" : ""));
        if (handle == IntPtr.Zero) return (false, $"the loaded native library ({library}) has no handle");
        var missing = Array.FindAll(NumericCheckExports, name => !hasExport(name));
        if (missing.Length > 0) return (false, $"the loaded native library ({library}) has no numeric-check entry points ({string.Join(", ", missing)}; older build)");
        int v = abiVersion();
        if (v != NumericChecksAbiVersionExpected) return (false, $"numeric-check ABI version {v} (expected {NumericChecksAbiVersionExpected}) in {library}");
        return (true, null);
    }

    private static readonly Lazy<(bool Ok, string Reason)> numericChecksProbe = new(() =>
        EvaluateNumericChecksProbe(Enabled, LastLoadError, LoadedNativeHandle, LoadedNativeLibrary,
            name => NativeLibrary.TryGetExport(LoadedNativeHandle, name, out _), () => _NumericChecksAbiVersion()));

    /// <summary>261118Cl 追加: 検査つき・生値の native の入口が使えるか (読み込んだ DLL に 8 本すべてがあり ABI の版が合う)</summary>
    public static bool NumericChecksAvailable => numericChecksProbe.Value.Ok;
    /// <summary>261118Cl 追加: 使えないときの理由 (使えるなら null)</summary>
    public static string NumericChecksUnavailableReason => numericChecksProbe.Value.Reason;

    private static void RequireNumericChecks()
    {
        if (!NumericChecksAvailable) throw new InvalidOperationException("numeric-check entry points are not available: " + NumericChecksUnavailableReason);
    }
    private static void RequireLength<T>(T[] a, long n, string name)
    {
        if (a is null || a.LongLength < n) throw new ArgumentException($"{name} must have at least {n} elements (has {a?.LongLength ?? 0})", name);
    }
    private static void RequireDiag(double[] diag)
    {
        if (diag is not null && diag.Length < NumericDiagLength) throw new ArgumentException($"diag must have at least {NumericDiagLength} elements", nameof(diag));
        if (diag is not null) Array.Fill(diag, double.NaN); // 261118Cl 追加 (codex2 第 6 回): 書かれなかった欄を前の呼び出しの値と取り違えない
    }

    /// <summary>261118Cl 追加: 検査つきの EVD (<see cref="EigenSolver(int, Complex[], ref Complex[], ref Complex[])"/> と同じ量)。
    /// 出力の配列は呼び出し側が用意する (長さ dim・dim²以上)。失敗のときは出力の [0,dim)・[0,dim²) が NaN になる</summary>
    public static unsafe NativeSolveStatus EigenSolverChecked(int dim, Complex[] mat, Complex[] eigenValues, Complex[] eigenVectors,
        NumericCheckFlags flags = NumericCheckFlags.None, double[] diag = null)
    {
        RequireNumericChecks();
        RequireLength(mat, (long)dim * dim, nameof(mat)); RequireLength(eigenValues, dim, nameof(eigenValues)); RequireLength(eigenVectors, (long)dim * dim, nameof(eigenVectors)); RequireDiag(diag);
        fixed (Complex* m = mat, v = eigenValues, w = eigenVectors)
        fixed (double* d = diag)
            return (NativeSolveStatus)_EigenSolverChecked(dim, (double*)m, (double*)v, (double*)w, (int)flags, d, diag?.Length ?? 0);
    }

    /// <summary>261118Cl 追加: 検査つきの逆行列 (<see cref="Inverse(int, Complex[])"/> と同じ量)。diag[0] にピボットの比</summary>
    public static unsafe NativeSolveStatus InverseChecked(int dim, Complex[] mat, Complex[] inverse, NumericCheckFlags flags = NumericCheckFlags.None, double[] diag = null)
    {
        RequireNumericChecks();
        RequireLength(mat, (long)dim * dim, nameof(mat)); RequireLength(inverse, (long)dim * dim, nameof(inverse)); RequireDiag(diag);
        fixed (Complex* m = mat, r = inverse)
        fixed (double* d = diag)
            return (NativeSolveStatus)_InverseChecked(dim, (double*)m, (double*)r, (int)flags, d, diag?.Length ?? 0);
    }

    /// <summary>261118Cl 追加: 検査つきの LU 解 (<see cref="PartialPivLuSolve(in int, Complex[], Complex[], ref Complex[])"/> と同じ量)</summary>
    public static unsafe NativeSolveStatus PartialPivLuSolveChecked(int dim, Complex[] mat, Complex[] vec, Complex[] result,
        NumericCheckFlags flags = NumericCheckFlags.None, double[] diag = null)
    {
        RequireNumericChecks();
        RequireLength(mat, (long)dim * dim, nameof(mat)); RequireLength(vec, dim, nameof(vec)); RequireLength(result, dim, nameof(result)); RequireDiag(diag);
        fixed (Complex* m = mat, v = vec, r = result)
        fixed (double* d = diag)
            return (NativeSolveStatus)_PartialPivLuSolveChecked(dim, (double*)m, (double*)v, (double*)r, (int)flags, d, diag?.Length ?? 0);
    }

    /// <summary>261118Cl 追加: 検査つきの CBED/STEM の伝播 (_CBEDSolver_Eigen と同じ量。result は dim × thickness.Length の列優先)。
    /// <see cref="NativeSolveStatus.NonFiniteOutput"/> のときは計算した値が残る</summary>
    public static unsafe NativeSolveStatus CBEDSolver_EigenChecked(int dim, Complex[] potential, Complex[] psi0, double[] thickness, Complex[] result,
        NumericCheckFlags flags = NumericCheckFlags.None, double[] diag = null)
    {
        RequireNumericChecks();
        int tDim = thickness?.Length ?? 0;
        RequireLength(potential, (long)dim * dim, nameof(potential)); RequireLength(psi0, dim, nameof(psi0)); if (tDim > 0) RequireLength(result, (long)dim * tDim, nameof(result)); RequireDiag(diag); // 261118Cl: 厚さが空なら result は読まない (native と同じ約束)
        fixed (Complex* p = potential, s = psi0, r = result)
        fixed (double* t = thickness, d = diag)
            return (NativeSolveStatus)_CBEDSolver_EigenChecked(dim, (double*)p, (double*)s, tDim, t, (double*)r, (int)flags, d, diag?.Length ?? 0);
    }

    /// <summary>261118Cl 追加: 検査つきの _CBEDSolver_Eigen2 (λ・V・α・Tg)</summary>
    public static unsafe NativeSolveStatus CBEDSolver_Eigen2Checked(int dim, Complex[] potential, Complex[] psi0, double[] thickness,
        Complex[] values, Complex[] vectors, Complex[] alphas, Complex[] tg, NumericCheckFlags flags = NumericCheckFlags.None, double[] diag = null)
    {
        RequireNumericChecks();
        int tDim = thickness?.Length ?? 0;
        RequireLength(potential, (long)dim * dim, nameof(potential)); RequireLength(psi0, dim, nameof(psi0)); RequireLength(values, dim, nameof(values));
        RequireLength(vectors, (long)dim * dim, nameof(vectors)); RequireLength(alphas, dim, nameof(alphas)); if (tDim > 0) RequireLength(tg, (long)dim * tDim, nameof(tg)); RequireDiag(diag); // 261118Cl: 厚さが空なら tg は読まない
        fixed (Complex* p = potential, s = psi0, v = values, w = vectors, a = alphas, g = tg)
        fixed (double* t = thickness, d = diag)
            return (NativeSolveStatus)_CBEDSolver_Eigen2Checked(dim, (double*)p, (double*)s, tDim, t, (double*)v, (double*)w, (double*)a, (double*)g, (int)flags, d, diag?.Length ?? 0);
    }

    /// <summary>261118Cl 追加: 生値の EBSD (<see cref="EBSDSolver(int, Complex[], Complex[], Complex[], Complex[], double[], double[])"/> の clamp の前の和)。
    /// 値の検査はしない (呼び出し側が源ごとに集計・判定する)。旧の意味の丸めは <see cref="LegacyEbsdClamp"/></summary>
    internal static unsafe NativeSolveStatus EBSDSolverRaw(int bLen, Complex[] eigenValues, Complex[] eigenVectors, Complex[] alpha,
        Complex[] phaseNG, double[] sigma, double[] thicknesses, double[] intensity)
    {
        RequireNumericChecks();
        ArgumentNullException.ThrowIfNull(sigma); ArgumentNullException.ThrowIfNull(thicknesses);
        int nAtoms = sigma.Length, tLen = thicknesses.Length;
        RequireLength(intensity, tLen, nameof(intensity));
        RequireEbsdInputs(bLen, nAtoms, tLen, eigenValues, eigenVectors, alpha, phaseNG, null, withTds: false); // 261118Cl 追加 (codex2 第 6 回): native が読む入力の長さ
        fixed (Complex* v = eigenValues, w = eigenVectors, a = alpha, p = phaseNG)
        fixed (double* s = sigma, t = thicknesses, r = intensity)
            return (NativeSolveStatus)_EBSDSolverRaw(bLen, nAtoms, tLen, (double*)v, (double*)w, (double*)a, (double*)p, s, t, r);
    }

    /// <summary>261118Cl 追加: 生値の EBSD + TDS (<see cref="EBSDSolverWithTDS(Complex[], Complex[], Complex[], Complex[], double[], Complex[], double, double[])"/> の clamp の前の和)</summary>
    internal static unsafe NativeSolveStatus EBSDSolverWithTDSRaw(int bLen, Complex[] eigenValues, Complex[] eigenVectors, Complex[] alpha,
        Complex[] phaseNG, double[] sigma, Complex[] muBack, double tdsCoeff, double[] thicknesses, double[] intensity, double[] tdsIntensity)
    {
        RequireNumericChecks();
        ArgumentNullException.ThrowIfNull(sigma); ArgumentNullException.ThrowIfNull(thicknesses);
        int nAtoms = sigma.Length, tLen = thicknesses.Length;
        RequireLength(intensity, tLen, nameof(intensity)); RequireLength(tdsIntensity, tLen, nameof(tdsIntensity));
        RequireEbsdInputs(bLen, nAtoms, tLen, eigenValues, eigenVectors, alpha, phaseNG, muBack, withTds: true); // 261118Cl 追加 (codex2 第 6 回)
        fixed (Complex* v = eigenValues, w = eigenVectors, a = alpha, p = phaseNG, u = muBack)
        fixed (double* s = sigma, t = thicknesses, r = intensity, q = tdsIntensity)
            return (NativeSolveStatus)_EBSDSolverWithTDSRaw(bLen, nAtoms, tLen, (double*)v, (double*)w, (double*)a, (double*)p, s, (double*)u, tdsCoeff, t, r, q);
    }

    /// <summary>261118Cl 追加 (codex2 第 6 回): 生値の EBSD の入口が読む入力の長さを、native を呼ぶ前に確かめる (短い非 null の配列は native の null の検査を通ってしまう)。
    /// 固有値・α は bLen、固有ベクトル・muBack は bLen²、位相は nAtoms × bLen (局所の原子が 0 のときは要らない)。tLen == 0 は何も読まない</summary>
    private static void RequireEbsdInputs(int bLen, int nAtoms, int tLen, Complex[] eigenValues, Complex[] eigenVectors, Complex[] alpha, Complex[] phaseNG, Complex[] muBack, bool withTds)
    {
        if (bLen < 0) throw new ArgumentOutOfRangeException(nameof(bLen), bLen, "bLen must be non-negative");
        if (tLen == 0 || (!withTds && nAtoms == 0)) return; // 何も読まない (局所の源だけで原子 0 は 0 を書いて戻る)
        long b2 = (long)bLen * bLen;
        RequireLength(eigenValues, bLen, nameof(eigenValues)); RequireLength(alpha, bLen, nameof(alpha)); RequireLength(eigenVectors, b2, nameof(eigenVectors));
        if (nAtoms > 0) RequireLength(phaseNG, (long)nAtoms * bLen, nameof(phaseNG));
        if (withTds) RequireLength(muBack, b2, nameof(muBack)); // TDS つきの入口では muBack を必ず読む
    }

    /// <summary>261118Cl 追加: 旧の native の EBSD の丸め max(0.0, sum) の意味 (11-13 版の DLL の実測): NaN・負・−0・−∞ は +0、正・+∞ はそのまま。
    /// C# には fast-math が無いので、この式どおりになる (managed の積分器の Math.Max(0, x) は NaN を通すので意味が違う)</summary>
    public static double LegacyEbsdClamp(double raw) => raw > 0 ? raw : 0.0;
}
