// 261110Cl 新規作成: EBSD のビーム選択の新方式 (着衣つきの裾) の native の入口 (Crystallography.Native/DressedTail.cpp、ABI の版 1)。
//   表 (エネルギーごと) と作業領域 (作業ごと) は明示の手綱 (SafeHandle) で持ち、使い終わったら Dispose で返す (codex2 第 31 回)。
using System;
using System.Numerics;
using System.Runtime.InteropServices;

namespace Crystallography;

public static partial class NativeWrapper
{
    /// <summary>261110Cl 追加: _DressedTailSolve などに渡す設定 (DressedTail.h の DressedTailSettings と同じ並び、ABI の版 1)</summary>
    [StructLayout(LayoutKind.Sequential)]
    internal struct DressedTailNativeSettings
    {
        public int AbiVersion, AllowUBf16;
        public double ModeTau, DegTol, EdeltaThreshold, RegDepthNm;
        public double Reserved0, Reserved1, Reserved2, Reserved3;
    }

    /// <summary>261110Cl 追加: 着衣つきの裾の native の戻り値 (DressedTail.h の DressedTailStatus)</summary>
    internal enum DressedTailStatus { Ok = 0, Stage1Evd = 1, FinalEvd = 2, Gram = 3, Capacity = 4, Args = 5, OutOfMemory = 6, NonFinite = 7, Unknown = 9 }

    internal const int DressedTailAbiVersionExpected = 1;

    [LibraryImport("Crystallography.Native.dll")] private static partial int _DressedTailAbiVersion();
    [LibraryImport("Crystallography.Native.dll")] private static partial int _DressedTailCapabilities();
    [LibraryImport("Crystallography.Native.dll")]
    private static unsafe partial IntPtr _DressedTailCreateTables(double* diffU, double* diffF, int diffLength, int DH, int DK, int DL,
        int* dHkl, float* dRe, float* dIm, int nD, int H, int K, int L, double* phase, int nPhaseCols, int nSites, double* wsite, int* status);
    [LibraryImport("Crystallography.Native.dll")] private static partial void _DressedTailFreeTables(IntPtr tables);
    [LibraryImport("Crystallography.Native.dll")] private static partial long _DressedTailTablesBytes(IntPtr tables);
    [LibraryImport("Crystallography.Native.dll")] private static unsafe partial IntPtr _DressedTailCreateWorkspace(int* status);
    [LibraryImport("Crystallography.Native.dll")] private static partial void _DressedTailFreeWorkspace(IntPtr workspace);
    [LibraryImport("Crystallography.Native.dll")] private static partial long _DressedTailWorkspaceBytes(IntPtr workspace);
    [LibraryImport("Crystallography.Native.dll")]
    private static unsafe partial int _DressedTailSolve(IntPtr tables, IntPtr workspace, DressedTailNativeSettings* settings,
        int nS, int nQ, int* sHkl, int* qHkl, double* Ds, double* Dq, double* hss, double* lamQ, int* qPhaseCol, double sqrtP0,
        double* gammaOut, double* sourceOut, int capacity, int* nM, double* seconds);

    private static readonly Lazy<(bool ok, string reason, int caps)> dressedTailProbe = new(() =>
    {
        if (!Enabled) return (false, "native library is not enabled" + (LastLoadError != null ? $" ({LastLoadError})" : ""), 0);
        if (LoadedNativeHandle == IntPtr.Zero || !NativeLibrary.TryGetExport(LoadedNativeHandle, "_DressedTailSolve", out _))
            return (false, $"the loaded native library ({LoadedNativeLibrary}) has no dressed-tail entry points (older build)", 0);
        int v = _DressedTailAbiVersion();
        if (v != DressedTailAbiVersionExpected) return (false, $"dressed-tail ABI version {v} (expected {DressedTailAbiVersionExpected})", 0);
        return (true, null, _DressedTailCapabilities());
    });

    /// <summary>261110Cl 追加: 着衣つきの裾の native の入口が使えるか (ロードした DLL に入口があり ABI の版が合う)</summary>
    public static bool DressedTailAvailable => dressedTailProbe.Value.ok;
    /// <summary>261110Cl 追加: 使えないときの理由 (使えるなら null)</summary>
    public static string DressedTailUnavailableReason => dressedTailProbe.Value.reason;
    /// <summary>261110Cl 追加: bit 0 = AVX-512 の版、bit 1 = U の bf16 が使える、bit 2 = AVX2 の版</summary>
    public static int DressedTailCapabilities => dressedTailProbe.Value.caps;

    /// <summary>261110Cl 追加: エネルギーごとの共有の表の手綱</summary>
    internal sealed class DressedTailTablesHandle : SafeHandle
    {
        public DressedTailTablesHandle() : base(IntPtr.Zero, true) { }
        public override bool IsInvalid => handle == IntPtr.Zero;
        protected override bool ReleaseHandle() { _DressedTailFreeTables(handle); return true; }
        public long Bytes => IsInvalid || IsClosed ? 0 : _DressedTailTablesBytes(handle);
        internal void Set(IntPtr h) => SetHandle(h);
    }

    /// <summary>261110Cl 追加: 1 つの作業が専有する作業領域の手綱</summary>
    internal sealed class DressedTailWorkspaceHandle : SafeHandle
    {
        public DressedTailWorkspaceHandle() : base(IntPtr.Zero, true) { }
        public override bool IsInvalid => handle == IntPtr.Zero;
        protected override bool ReleaseHandle() { _DressedTailFreeWorkspace(handle); return true; }
        public long Bytes => IsInvalid || IsClosed ? 0 : _DressedTailWorkspaceBytes(handle);
        internal void Set(IntPtr h) => SetHandle(h);
    }

    internal static unsafe DressedTailTablesHandle DressedTailCreateTables(Complex[] diffU, Complex[] diffF, int dh, int dk, int dl,
        int[] dHkl, float[] dRe, float[] dIm, int h, int k, int l, Complex[] phase, int nSites, double[] wsite)
    {
        if (diffU.Length != diffF.Length) throw new ArgumentException("diffU and diffF must have the same length");
        int nPhaseCols = nSites > 0 ? phase.Length / nSites : 0;
        var t = new DressedTailTablesHandle();
        int status;
        IntPtr p;
        fixed (Complex* u = diffU, f = diffF, ph = phase)
        fixed (int* d = dHkl)
        fixed (float* re = dRe, im = dIm)
        fixed (double* w = wsite)
            p = _DressedTailCreateTables((double*)u, (double*)f, diffU.Length, dh, dk, dl, d, re, im, dRe.Length, h, k, l, (double*)ph, nPhaseCols, nSites, w, &status);
        if (p == IntPtr.Zero) throw status == (int)DressedTailStatus.OutOfMemory ? new OutOfMemoryException("dressed-tail native tables") : new InvalidOperationException($"dressed-tail native tables failed ({(DressedTailStatus)status})");
        t.Set(p);
        return t;
    }

    internal static unsafe DressedTailWorkspaceHandle DressedTailCreateWorkspace()
    {
        var w = new DressedTailWorkspaceHandle();
        int status;
        var p = _DressedTailCreateWorkspace(&status);
        if (p == IntPtr.Zero) throw new OutOfMemoryException($"dressed-tail native workspace ({(DressedTailStatus)status})");
        w.Set(p);
        return w;
    }

    /// <summary>261110Cl 追加: 1 方向を解く。成功なら (γ, B = Z†F_W Z (nM × nM、列優先))、失敗なら例外 (理由つき)</summary>
    internal static unsafe (Complex[] gamma, Complex[] source) DressedTailSolve(DressedTailTablesHandle tables, DressedTailWorkspaceHandle workspace, in DressedTailNativeSettings settings,
        int[] sHkl, int[] qHkl, double[] ds, double[] dq, Complex[] hss, Complex[] lamQ, int[] qPhaseCol, double sqrtP0, double[] seconds = null)
    {
        int nS = ds.Length, nQ = dq.Length;
        // 261110Cl 追加 (codex2 第 32 回): 配列の長さを検査する (native は長さを知らない)
        if (sHkl.Length != 3 * nS || qHkl.Length != 3 * nQ || hss.Length != nS * nS || lamQ.Length != nQ || qPhaseCol.Length != nQ || (seconds != null && seconds.Length < 9))
            throw new ArgumentException("dressed-tail native solve: inconsistent array lengths");
        var gamma = new Complex[nS]; var src = new Complex[nS * nS];
        int nM = 0, rc;
        bool addT = false, addW = false;
        try
        {
            tables.DangerousAddRef(ref addT); workspace.DangerousAddRef(ref addW);
            var st = settings;
            fixed (int* s = sHkl, q = qHkl, pc = qPhaseCol)
            fixed (double* pDs = ds, pDq = dq, sec = seconds)
            fixed (Complex* h = hss, lq = lamQ, g = gamma, so = src)
                rc = _DressedTailSolve(tables.DangerousGetHandle(), workspace.DangerousGetHandle(), &st, nS, nQ, s, q, pDs, pDq, (double*)h, (double*)lq, pc, sqrtP0,
                    (double*)g, (double*)so, nS, &nM, sec);
        }
        finally
        {
            if (addW) workspace.DangerousRelease();
            if (addT) tables.DangerousRelease();
        }
        if (rc != 0 || nM <= 0)
            throw rc == (int)DressedTailStatus.OutOfMemory ? new OutOfMemoryException("dressed-tail native solve") : new InvalidOperationException($"dressed-tail native solve failed ({(DressedTailStatus)rc})");
        if (nM == nS) return (gamma, src);
        var g2 = new Complex[nM]; Array.Copy(gamma, g2, nM);
        var s2 = new Complex[nM * nM]; Array.Copy(src, s2, nM * nM);
        return (g2, s2);
    }
}
