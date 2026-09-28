#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;

namespace Crystallography;

/// <summary>260925Cl 追加: EDX の線の系列 (Kα・Kβ・Lα・Lβ・Mα) の「一次の空孔 → 光子」の係数。
/// <para>
/// 線の系列 s の光子の生成断面積 σ_s = Σ_i B_si σ_i、実効の形 F_s = Σ_i (B_si σ_i / σ_s) F_i
/// (i = 一次に電離された殻、σ_i・F_i = 既存の Bote–Salvat の σ と Temari の F)。
/// B_si = Σ_{系列の線 ℓ} P_k(i) · ω_k · R_ℓ (k = 線 ℓ の始状態の殻、P_k(i) = 殻 i の空孔 1 個が連鎖の後に殻 k に作る空孔の数)。
/// 連鎖で移った空孔は最初に電離された殻の形 F_i を持つ (局所の近似。EDX は高速電子の終状態を全部積分する)。
/// </para>
/// <para>
/// 連鎖の模型 = xraylib の Kissel の連鎖の <b>radiative</b> (Coster–Kronig + 内殻からの放射遷移)。
/// **Auger の連鎖は含まない** — 同梱の <c>libxrl-11.dll</c> は Auger rates を削った build (作者決定 2026-09-25、Temari の I86)。
/// Temari の <c>temari_engine.xray_lines</c> (cascade="radiative") と同じ組み方・同じ系列の定義で、そこでの影響の実測は
/// K 線は full と同じ、Lα (Z ≤ 50) は 0.4〜1.7 % 小さい、Mα は約 5 % 小さい (L → M の Auger の分)。
/// </para>
/// <para>
/// 系列の定義は xraylib の合成の線と同じ: Kα = KA_LINE、Kβ = KB_LINE (= 1 − Kα)、Lα = LA_LINE (L3)、
/// Lβ = xraylib の LB_LINE_MACROS の 13 本の和 (L1・L2・L3 の混合)、Mα = MA1 + MA2 (M5N7 + M5N6、xraylib に合成は無い)。
/// </para>
/// </summary>
public static class XrayLineSeries
{
    // 殻の添字 = xraylib の殻の番号 (K=0, L1=1, L2=2, L3=3, M1=4 … M5=8)
    private const int NShell = 9;
    private const int K = 0, L1 = 1, L2 = 2, L3 = 3, M1 = 4, M2 = 5, M3 = 6, M4 = 7, M5 = 8;

    /// <summary>殻の添字 → F テーブルの shellCode (<see cref="IonizationFsTable"/>)。L23 (v1) は使わない。</summary>
    internal static readonly int[] ShellCodeOf =
    [
        IonizationFsTable.ShellCodeK, IonizationFsTable.ShellCodeL1, IonizationFsTable.ShellCodeL2, IonizationFsTable.ShellCodeL3,
        IonizationFsTable.ShellCodeM1, IonizationFsTable.ShellCodeM2, IonizationFsTable.ShellCodeM3, IonizationFsTable.ShellCodeM4,
        IonizationFsTable.ShellCodeM5,
    ];

    // 放射遷移で内殻 j の空孔が殻 k に移る線 (xraylib の符号): (k, j, line)。xraylib の P*_rad_cascade_kissel と同じ組
    private static readonly (int k, int j, int line)[] RadiativeTransfers =
    [
        (L1, K, -1), (L2, K, -2), (L3, K, -3),                                    // KL1 KL2 KL3
        (M1, K, -4), (M1, L1, -32), (M1, L2, -60), (M1, L3, -86),                 // KM1 L1M1 L2M1 L3M1
        (M2, K, -5), (M2, L1, -33), (M2, L2, -61), (M2, L3, -87),
        (M3, K, -6), (M3, L1, -34), (M3, L2, -62), (M3, L3, -88),
        (M4, K, -7), (M4, L1, -35), (M4, L2, -63), (M4, L3, -89),
        (M5, K, -8), (M5, L1, -36), (M5, L2, -64), (M5, L3, -90),
    ];

    // Coster–Kronig (同じ主量子数の殻の中): (k, j, 遷移)。L3 ← L1 は f13 + f′13
    private static readonly (int k, int j, Xraylib.XrlTrans t)[] CosterKronigTransfers =
    [
        (L2, L1, Xraylib.XrlTrans.FL12), (L3, L1, Xraylib.XrlTrans.FL13), (L3, L1, Xraylib.XrlTrans.FLP13), (L3, L2, Xraylib.XrlTrans.FL23),
        (M2, M1, Xraylib.XrlTrans.FM12), (M3, M1, Xraylib.XrlTrans.FM13), (M3, M2, Xraylib.XrlTrans.FM23),
        (M4, M1, Xraylib.XrlTrans.FM14), (M4, M2, Xraylib.XrlTrans.FM24), (M4, M3, Xraylib.XrlTrans.FM34),
        (M5, M1, Xraylib.XrlTrans.FM15), (M5, M2, Xraylib.XrlTrans.FM25), (M5, M3, Xraylib.XrlTrans.FM35), (M5, M4, Xraylib.XrlTrans.FM45),
    ];

    // 系列の線: (xraylib の符号, 始状態の殻)。KA/KB/LA は xraylib の合成の符号 (RadRate が殻の中の合計を返す)
    private static readonly (int line, int shell)[] KAlphaLines = [(0, K)];
    private static readonly (int line, int shell)[] KBetaLines = [(1, K)];
    private static readonly (int line, int shell)[] LAlphaLines = [(2, L3)];
    private static readonly (int line, int shell)[] LBetaLines =        // xraylib の LB_LINE_MACROS (kissel_pe.c)
    [
        (-63, L2), (-95, L3), (-34, L1), (-33, L1), (-102, L3), (-91, L3), (-98, L3),   // LB1 LB2 LB3 LB4 LB5 LB6 LB7
        (-36, L1), (-35, L1), (-94, L3), (-62, L2), (-96, L3), (-97, L3),               // LB9 LB10 LB15 LB17 L3N6 L3N7
    ];
    private static readonly (int line, int shell)[] MAlphaLines = [(-207, M5), (-206, M5)];   // MA1 = M5N7、MA2 = M5N6

    /// <summary>線の系列か (殻ではなく X 線の線を数えるチャネル)。</summary>
    public static bool IsLineSeries(IonizationShell shell)
        => shell is IonizationShell.KAlpha or IonizationShell.KBeta or IonizationShell.LAlpha or IonizationShell.LBeta or IonizationShell.MAlpha;

    private static (int line, int shell)[] LinesOf(IonizationShell series) => series switch
    {
        IonizationShell.KAlpha => KAlphaLines,
        IonizationShell.KBeta => KBetaLines,
        IonizationShell.LAlpha => LAlphaLines,
        IonizationShell.LBeta => LBetaLines,
        IonizationShell.MAlpha => MAlphaLines,
        _ => throw new ArgumentOutOfRangeException(nameof(series), series, "not an X-ray line series"),
    };

    /// <summary>系列の始状態の殻の添字 (この殻が F テーブルに無ければ系列は出せない)。</summary>
    internal static int[] InitialShellsOf(IonizationShell series) => [.. LinesOf(series).Select(l => l.shell).Distinct()];

    /// <summary>xraylib の値の NaN (「データが無い」も「0」も含む) を 0 に畳む。xraylib 自身の連鎖も 0 として扱う。</summary>
    private static double Z0(double v) => double.IsFinite(v) && v > 0 ? v : 0.0;

    /// <summary>空孔の行列 P[k, i] = 殻 i の一次の空孔 1 個が連鎖の後に殻 k に作る空孔の数 (期待値。規格化しない)。
    /// <paramref name="radiative"/> = false なら Coster–Kronig だけ (xraylib の pure)。</summary>
    public static double[,] VacancyMatrix(int z, bool radiative = true)
    {
        var omega = new double[NShell];
        for (int s = 0; s < NShell; s++) omega[s] = Z0(Xraylib.FluorescenceYield(z, (Xraylib.XrlShell)s));
        var t = new double[NShell, NShell];   // 直接の移送 t[k, j]
        if (radiative)
            foreach (var (k, j, line) in RadiativeTransfers)
                t[k, j] += omega[j] * Z0(Xraylib.LineRadRate(z, (Xraylib.XrlLine)line));
        foreach (var (k, j, tr) in CosterKronigTransfers)
            t[k, j] += Z0(Xraylib.CosterKronig(z, tr));
        // 殻の順 (K → M5) に解く: P[k, i] = δ_ki + Σ_j t[k, j] P[j, i] (j は k より内側)
        var p = new double[NShell, NShell];
        for (int k = 0; k < NShell; k++)
            for (int i = 0; i < NShell; i++)
            {
                double v = k == i ? 1.0 : 0.0;
                for (int j = 0; j < k; j++) v += t[k, j] * p[j, i];
                p[k, i] = v;
            }
        return p;
    }

    /// <summary>系列 <paramref name="series"/> の B_si (殻 i = 0..8 の一次の空孔 1 個あたりに出る系列の光子の数、全方位)。
    /// 系列の線が 1 本も無い (R = 0)、または R > 0 の線の始状態の殻に xraylib の蛍光収率が無いとき (例: M5 の ω は Z 48–57 で 0) は
    /// null と理由を返す — 0 と区別する (データの穴を黙って 0 にしない)。</summary>
    public static double[]? Yields(int z, IonizationShell series, out string? reason, bool radiative = true)
    {
        reason = null;
        if (!Xraylib.Enabled || !Xraylib.CosterKronigEnabled)
        {
            reason = "xraylib (with Coster-Kronig) is not available: " + (Xraylib.LastLoadError ?? Xraylib.CosterKronigLoadError);
            return null;
        }
        var lines = LinesOf(series).Select(l => (l.line, l.shell, r: Z0(Xraylib.LineRadRate(z, (Xraylib.XrlLine)l.line)))).Where(l => l.r > 0).ToArray();
        if (lines.Length == 0)
        {
            reason = $"xraylib has no {series} line for Z={z}";
            return null;
        }
        var noYield = lines.Select(l => l.shell).Distinct().Where(s => Z0(Xraylib.FluorescenceYield(z, (Xraylib.XrlShell)s)) == 0).ToArray();
        if (noYield.Length > 0)
        {
            reason = $"xraylib has no fluorescence yield for shell {string.Join(",", noYield.Select(s => (Xraylib.XrlShell)s))} (Z={z})";
            return null;
        }
        var p = VacancyMatrix(z, radiative);
        var b = new double[NShell];
        foreach (var (line, k, r) in lines)
        {
            double omegaR = Xraylib.FluorescenceYield(z, (Xraylib.XrlShell)k) * r;
            for (int i = 0; i < NShell; i++) b[i] += p[k, i] * omegaR;
        }
        return b;
    }
}
