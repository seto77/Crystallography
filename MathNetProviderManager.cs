using System;

namespace Crystallography;

/// <summary>
/// 260711Cl 追加: MathNet の LinearAlgebra provider (MKL native / managed) をプロセス起動時に一元初期化する。
///
/// 旧構成の問題 (codex 調査で確認):
/// (1) BetheMethod / Geometry / PeakSearch の static ctor が無条件に TryUseNativeMKL() を呼んでおり、
///     GUI の「Use MKL」メニュー状態と実際の provider が食い違い得た (未チェックでも DLL があれば MKL が有効)。
/// (2) MathNet 6.0.0-beta2 の MKL provider は「ロード時の Control.MaxDegreeOfParallelism を set_max_threads に渡し、
///     以後の変更は反映しない」ため、既定 (=論理コア数) でロードされると EVD の外側並列 (方向並列) と
///     オーバーサブクリプションを起こす。実測 (9950X) では 外32×内1 が最速で、外8×内4 は半分以下の性能。
///
/// 本クラスは「MKL native DLL が配置済みなら内側スレッド数 1 でロードし、無ければ managed を明示する」を
/// プロセスで一度だけ行う。GUI (Program.cs) と BetheBench の両方がここを通る。
/// </summary>
public static class MathNetProviderManager
{
    /// <summary>MKL native ライブラリのファイル名 (実行時ダウンロード方式。GUI の Download メニューが取得する)</summary>
    public static readonly string[] MklFileNames = ["libMathNetNumericsMKL.dll", "libiomp5md.dll"];

    /// <summary>現在の LinearAlgebra provider が MKL native かどうか (Initialize 後に確定)</summary>
    public static bool MklActive { get; private set; } = false;

    /// <summary>260712Cl 追加: 環境変数 RECIPRO_DISABLE_MKL=1 による強制無効化が指定されているか (GUI のメニュー可視性判定用)</summary>
    public static bool MklDisabledByEnvironment => Environment.GetEnvironmentVariable("RECIPRO_DISABLE_MKL") == "1";

    private static readonly object initLock = new();//260712Cl 追加 (codex 指摘: 多重初期化の競合防止)

    /// <summary>MKL native DLL が実行フォルダに配置済みかどうか</summary>
    public static bool MklNativeFilesExist()
    {
        var dir = AppDomain.CurrentDomain.BaseDirectory;
        foreach (var f in MklFileNames)
            if (!System.IO.File.Exists(System.IO.Path.Combine(dir, f)))
                return false;
        return true;
    }

    /// <summary>
    /// プロセス起動時 (MathNet を使う計算の前) に呼ぶ。useMkl=true なら MKL native provider を内側スレッド数
    /// mklInnerThreads でロードし、失敗時および useMkl=false は managed provider を明示する。
    /// MKL は x64 専用 (ARM64 では ロードせず managed へ)。
    /// 環境変数 RECIPRO_DISABLE_MKL=1 で強制無効化 (DL 済みでも使いたくない場合の非常口)。
    ///
    /// 260712Cl 状態遷移の規約 (codex レビュー反映): **Managed→MKL の一方向のみ**。
    /// - 既に MklActive の場合は何もしない。MKL はロード後のスレッド数変更を受け付けないため、
    ///   異なる mklInnerThreads での再呼び出しも無視される (内側スレッド数は初回ロード時のみ有効)
    /// - MKL→Managed へ戻す手段は提供しない (provider は process-global で、計算中の切替は危険)
    /// </summary>
    public static void Initialize(bool useMkl, int mklInnerThreads = 1)
    {
      lock (initLock)//260712Cl 追加 (codex 指摘: 多重初期化の競合防止)
      {
        if (MklActive) return;

        if (mklInnerThreads < 1) mklInnerThreads = 1;//260712Cl 追加 (codex 指摘: 0/負数の防御)
        if (MklDisabledByEnvironment)
            useMkl = false;
        if (System.Runtime.InteropServices.RuntimeInformation.ProcessArchitecture != System.Runtime.InteropServices.Architecture.X64)
            useMkl = false;

        if (useMkl)
        {
            //MKL provider は初回ロード時に set_max_threads(Control.MaxDegreeOfParallelism) を呼ぶため、
            //ロードの瞬間だけ MaxDoP を内側スレッド数に落とし、直後に復元する (managed 並列は実行時参照なので影響なし)
            var saved = MathNet.Numerics.Control.MaxDegreeOfParallelism;
            try
            {
                MathNet.Numerics.Control.MaxDegreeOfParallelism = mklInnerThreads;
                Environment.SetEnvironmentVariable("MKL_NUM_THREADS", mklInnerThreads.ToString());
                Environment.SetEnvironmentVariable("OMP_NUM_THREADS", mklInnerThreads.ToString());
                //TryUseNativeMKL の戻り値だけでは実 provider を保証しないため、実型も確認する (codex 指摘)
                MklActive = MathNet.Numerics.Control.TryUseNativeMKL() &&
                    MathNet.Numerics.Providers.LinearAlgebra.LinearAlgebraControl.Provider.GetType().Name.Contains("Mkl");
            }
            catch { MklActive = false; }
            finally { MathNet.Numerics.Control.MaxDegreeOfParallelism = saved; }
        }

        if (!MklActive)
            MathNet.Numerics.Control.UseManaged();
      }
    }
}
