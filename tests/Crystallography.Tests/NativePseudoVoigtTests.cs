using System;
using System.Runtime.InteropServices;
using Xunit;

namespace Crystallography.Tests;

/// <summary>261126Cl 追加 (作者の決定 11-26・codex2 の相談 n2 §7): native の 2D 擬 Voigt (_PseudoVoigt・_PseudoVoigtDiff) が、
/// C# の Marquardt.PseudoVoigt / PseudoVoigtDiff と同じ正規化した式であること。製品の経路からは呼んでいない (Marquardt の呼び出しはコメント内) ので、
/// 読み込んだ DLL の export を直接呼ぶ。入口が無い・native が読めないときは Skip (CRYSTALLOGRAPHY_REQUIRE_NUMERIC_CHECKS=1 なら失敗)</summary>
public sealed class NativePseudoVoigtFactAttribute : FactAttribute
{
    public NativePseudoVoigtFactAttribute()
    {
        if (NativePseudoVoigtTests.Load() == null && Environment.GetEnvironmentVariable("CRYSTALLOGRAPHY_REQUIRE_NUMERIC_CHECKS") != "1")
            Skip = $"Crystallography.Native with _PseudoVoigt / _PseudoVoigtDiff is required ({NativeWrapper.LoadedNativeLibrary ?? "not loaded"}).";
    }
}

public class NativePseudoVoigtTests
{
    // C ABI: x と results はポインタ、7 個のスカラーは値渡し (PseudoVoigt.h)
    [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
    internal delegate void PseudoVoigtFn(double[] x, double x0, double y0, double hx, double hy, double t, double eta, double a, double[] results);

    internal static (PseudoVoigtFn Value, PseudoVoigtFn Diff)? Load()
    {
        if (!NativeWrapper.Enabled || NativeWrapper.LoadedNativeHandle == IntPtr.Zero
            || !NativeLibrary.TryGetExport(NativeWrapper.LoadedNativeHandle, "_PseudoVoigt", out var v)
            || !NativeLibrary.TryGetExport(NativeWrapper.LoadedNativeHandle, "_PseudoVoigtDiff", out var d))
            return null;
        return (Marshal.GetDelegateForFunctionPointer<PseudoVoigtFn>(v), Marshal.GetDelegateForFunctionPointer<PseudoVoigtFn>(d));
    }

    [NativePseudoVoigtFact]
    public void NativeMatchesManaged()
    {
        var fns = Load();
        Assert.True(fns != null, $"Crystallography.Native with _PseudoVoigt is required ({NativeWrapper.LoadedNativeLibrary ?? "not loaded"}).");
        var (value, diff) = fns.Value;
        const double tol = 64 * 2.2204460492503131e-16;
        foreach (var xy in new[] { new[] { 2.0, -1.0 }, new[] { 0.25, -0.5 }, new[] { -3.5, 4.0 } })
            foreach (double eta in new[] { 0.0, 0.375, 1.0 })
            {
                var r = new double[1];
                value(xy, 0.25, -0.5, 1.5, 0.75, 0.375, eta, 12, r);
                double m = Marquardt.PseudoVoigt(xy, 0.25, -0.5, 1.5, 0.75, 0.375, eta, 12);
                Assert.True(Math.Abs(r[0] - m) <= tol * Math.Abs(m), $"value: native {r[0]:R}, managed {m:R}");
                var rd = new double[7];
                diff(xy, 0.25, -0.5, 1.5, 0.75, 0.375, eta, 12, rd);
                var md = Marquardt.PseudoVoigtDiff(xy, 0.25, -0.5, 1.5, 0.75, 0.375, eta, 12);
                double scale = 0; foreach (var e in md) scale = Math.Max(scale, Math.Abs(e)); // 中心では x0・y0 の微分が 0 になるので、成分の最大で尺度を取る
                for (int k = 0; k < 7; k++)
                    Assert.True(Math.Abs(rd[k] - md[k]) <= tol * scale, $"diff[{k}]: native {rd[k]:R}, managed {md[k]:R}");
            }
    }
}
