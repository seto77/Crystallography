// 261110Cl 新規作成: EBSD のビーム選択の新方式 (着衣つきの裾) の本番の経路。作者の判断 (2026-11-10): 移植の候補 = N300・τ3e−3・一体化した native の解、時間 3.7 倍を許す。
//   設計 = .project-guidance/ReciPro/ReciPro_ビーム選択_新方式_移植計画_261110.md (codex2 第 31 回の点検を反映)。
//   RunEBSDNew の代表方向ごとの解に差し込む (研究用の口 EbsdDressedTailResearch と同じ所・同じ外枠 SolveEbsdDressedTailVoltage)。
//   1 方向の解き手だけを本番の native (Crystallography.Native/DressedTail.cpp) にし、選択の設定は引数で渡す (研究用の static を使わない)。
//   従来の経路へ戻す判断は開始前・マスター全体で (native の入口が無い・非対応の設定・メモリの予算が足りない)。固有分解の失敗やキャンセルは失敗として止める。
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Runtime.InteropServices;

namespace Crystallography;

public partial class BetheMethod
{
    /// <summary>261110Cl 追加: EBSD のマスターパターンのビーム選択の方式</summary>
    public enum EbsdBeamSelectionMethod
    {
        /// <summary>従来 (Find_gVectors で N 本、弱いビームは捨てる)</summary>
        Legacy = 0,
        /// <summary>着衣つきの裾 (E1 の選択 + 遅い反射の繰り上げ + 残りを静的な応答として着せる。後方 TDS の源だけ)</summary>
        DressedTail = 1,
    }

    /// <summary>261110Cl 追加: 要求する方式 (既定 = 従来)。新方式が使えない条件なら構築の開始前に従来へ戻し、理由を <see cref="EbsdBeamSelectionFallbackReason"/> に残す</summary>
    public EbsdBeamSelectionMethod EbsdBeamSelection { get; set; } = EbsdBeamSelectionMethod.Legacy;
    /// <summary>261110Cl 追加: 直前の構築で実際に使った方式</summary>
    public EbsdBeamSelectionMethod EbsdBeamSelectionUsed { get; private set; } = EbsdBeamSelectionMethod.Legacy;
    /// <summary>261110Cl 追加: 新方式を要求したのに従来へ戻したときの理由 (戻さなければ null)</summary>
    public string EbsdBeamSelectionFallbackReason { get; private set; }
    /// <summary>261110Cl 追加: 直前の構築で新方式が使った作業の数 (並列の数)</summary>
    public int EbsdBeamSelectionWorkers { get; private set; }
    /// <summary>261110Cl 追加: 直前の構築の新方式の資源の見積もり (共有・1 作業・予算・予算で許す作業の数) の文 (表示と記録用)</summary>
    public string EbsdBeamSelectionDetail { get; private set; }
    /// <summary>261110Cl 追加 (codex2 第 32 回): 新方式を必須にする (研究・q の検証用)。使えない条件なら従来へ戻さず、構築の開始前に例外で止める</summary>
    public bool EbsdBeamSelectionRequireDressedTail { get; set; }
    /// <summary>261110Cl 追加 (codex2 第 32 回): 直前の構築で新方式が使った native の版 (DLL の名前・機能の印・U の bf16 を使ったか)</summary>
    public string EbsdBeamSelectionNative { get; private set; } = "";
    /// <summary>261110Cl 追加: 直前の構築の新方式の設定の版 (既定と違う設定なら +custom)</summary>
    public string EbsdBeamSelectionSettingsVersion { get; private set; } = "";

    /// <summary>261110Cl 追加: 新方式の本番の設定 (研究で確定した 1 つの組。仕様書 §12 末尾 (aaa)〜(hhh)、作者の判断 2026-11-10)。値で比べる (既定と違えば記録の版に +custom)</summary>
    public sealed record DressedTailProductionSettings
    {
        /// <summary>設定の版 (保存情報とキャッシュの区別に使う)</summary>
        public const string Version = "DT1-261110";
        public double PoolRadiusNm { get; init; } = 64;
        public double TailRadiusNm { get; init; } = 64;
        public double CoreRadiusNm { get; init; } = 44;
        public double NeighborRadiusNm { get; init; } = 10;
        public double NeighborPs { get; init; } = 0.3;
        public double Eta { get; init; } = 0.03;
        /// <summary>E1 の基準の本数 (000 を含む。遅い反射の繰り上げで明示の本数は増える)</summary>
        public int N { get; init; } = 300;
        public double Ps { get; init; } = 0.1;
        public double ModeTau { get; init; } = 3e-3;
        public double DegTol { get; init; } = 1e-4;
        public double EdeltaThreshold { get; init; } = 1e-2;
        public double RegDepthNm { get; init; } = 45;
        public int SourceQuadrature { get; init; } = 77;
        public bool AllowUBf16 { get; init; } = true;
        // /// <summary>作業の数の上限 (研究の経路は論理 CPU 24 の機械で 24 並列が最速、32 では遅い。他の CPU で 24 が最適とは言わない)</summary> // 261110Cl 変更前 (機械の記述の誤り)
        /// <summary>作業の数の上限 (研究の経路は 16 コア・32 論理 CPU の機械 (Ryzen 9 9950X) で 24 並列が最速、32 では SMT の取り合いで遅い。他の CPU で 24 が最適とは言わない)</summary>
        public int MaxWorkers { get; init; } = 24;

        /// <summary>261110Cl 追加 (codex2 第 34 回): 数値に効く設定の全部の指紋 (文化に依らない文字列の SHA-256 の先頭 8 桁)。既定と違う設定の版を区別する (N300 と N600 など)。作業の数の上限は結果を変えないので含めない</summary>
        public string Fingerprint()
        {
            string text = string.Create(System.Globalization.CultureInfo.InvariantCulture,
                $"pool={PoolRadiusNm:R};tail={TailRadiusNm:R};core={CoreRadiusNm:R};rn={NeighborRadiusNm:R};p2={NeighborPs:R};eta={Eta:R};N={N};ps={Ps:R};tau={ModeTau:R};deg={DegTol:R};ed={EdeltaThreshold:R};reg={RegDepthNm:R};quad={SourceQuadrature};ubf16={AllowUBf16}");
            return Convert.ToHexString(System.Security.Cryptography.SHA256.HashData(System.Text.Encoding.UTF8.GetBytes(text)))[..8].ToLowerInvariant();
        }
    }
    /// <summary>261110Cl 追加: 新方式の本番の設定 (既定)</summary>
    public static DressedTailProductionSettings DressedTailDefaults { get; } = new();
    /// <summary>261110Cl 追加 (作者の判断 2026-11-10「EBSD は 3 択」): 新方式の軽量の設定 (手早く確かめる用途)。E1 の基準の本数 32、候補と裾の半径 40 nm⁻¹、芯 31 nm⁻¹。
    /// Si 20 keV の本番の重みで最大 2.8 % (従来 N32 は 27 %)、格子 256 の 1 エネルギー 約 26 s (従来 N32 は 1 s)。MnAl 14 keV では最大 13 % が残る (従来 28 %) ので、定量の当てはめの検証済みの設定ではない。
    /// 版は DT1-261110-light</summary>
    public static DressedTailProductionSettings DressedTailLight { get; } = new() { N = 32, PoolRadiusNm = 40, TailRadiusNm = 40, CoreRadiusNm = 31 };
    /// <summary>261110Cl 追加: この BetheMethod の構築で使う新方式の設定 (既定 = <see cref="DressedTailDefaults"/>。構築の要求から渡す。研究の照合では値を変えた複製を渡せる)</summary>
    public DressedTailProductionSettings EbsdDressedTailSettings { get; set; } = DressedTailDefaults;

    /// <summary>261110Cl 追加: このプロセスで走っている新方式の構築が見込んだメモリの合計 [byte] (同時に走る構築が予算を二重に使わないため)</summary>
    static long s_dressedTailReserved;
    static readonly object s_dressedTailBudgetGate = new();

    /// <summary>261110Cl 追加: 構築の開始前に、新方式を使えるか・作業の数を決める。使えなければ (null, 理由)</summary>
    private (EbsdDressedTailOptions opt, DressedTailNativeRun run, string reason) TryCreateDressedTailProduction(int threads, bool tds)
    {
        var set = EbsdDressedTailSettings ?? DressedTailDefaults;
        if (!NativeWrapper.DressedTailAvailable) return (null, null, "native: " + NativeWrapper.DressedTailUnavailableReason);
        if (!tds) return (null, null, "the dressed-tail method supports only the TDS backscatter source (IncludeTDSBackground)");
        if (UseNonLocalAbsorption) return (null, null, "non-local absorption is not supported by the dressed-tail method");
        if (IncludeAbsorbedFluxBackground) return (null, null, "absorbed-flux reinjection is not supported by the dressed-tail method");
        if (double.IsFinite(SourceSmearingSigmaNm)) return (null, null, "source smearing is not supported by the dressed-tail method");
        var (shared, perWorker) = EstimateDressedTailMemory(set);
        // var info = GC.GetGCMemoryInfo(); // 261110Cl 変更前 (GC の値はヒープの上限 (GCHeapHardLimit) を反映し、機械のメモリではない)
        // long total = info.TotalAvailableMemoryBytes, free = Math.Max(0, total - info.MemoryLoadBytes);
        DressedTailNativeRun run; int workers;
        lock (s_dressedTailBudgetGate) // 261110Cl 追加: 予算の計算と見込みの登録を一度に (同時に始めた構築が同じ空きを二重に数えない)
        {
        var (total, free) = PhysicalMemory(); // 261110Cl: OS の物理メモリと使えるメモリ (GlobalMemoryStatusEx)
        long reserve = Math.Max(2L << 30, (long)(0.15 * total));
        // long budget = Math.Min(total / 2, free - reserve); // 261110Cl 変更前
        long others = s_dressedTailReserved; // 261110Cl 追加: 走っているほかの構築の見込み (もう確保した分は free にも出るので保守的に二重に数える)
        long budget = Math.Min(total / 2, free - reserve) - others;
        int byBudget = perWorker > 0 ? (int)Math.Min(int.MaxValue, Math.Max(0, (budget - shared) / perWorker)) : 0;
        workers = Math.Min(Math.Min(Math.Max(1, threads), set.MaxWorkers), Math.Min(Environment.ProcessorCount, byBudget));
        EbsdBeamSelectionDetail = $"memory estimate: shared {shared / 1e9:f2} GB, {perWorker / 1e9:f2} GB per worker; budget {budget / 1e9:f1} GB (total {total / 1e9:f1}, free {free / 1e9:f1}, reserve {reserve / 1e9:f1}, other builds {others / 1e9:f1}); workers {workers} (requested {threads}, cap {set.MaxWorkers}, logical CPUs {Environment.ProcessorCount}, by budget {byBudget})";
        if (workers < 1)
            return (null, null, $"estimated memory {(shared + perWorker) / 1e9:f1} GB (shared {shared / 1e9:f1} GB + {perWorker / 1e9:f1} GB per worker) exceeds the budget {budget / 1e9:f1} GB");
        // var run = new DressedTailNativeRun(set, workers); // 261110Cl 変更前
        // 261110Cl 変更 (codex2 第 34 回): 1 作業に許す量と登録する見込みを一致させる (許す量 = 見積もりの 2 倍まで、予算の残りの等分まで)。前は登録が見積もりで許す量がそれより大きく、同時の構築の安全の枠に抜けがあった
        long allowedPerWorker = Math.Max(perWorker, Math.Min(2 * perWorker, (budget - shared) / workers));
        run = new DressedTailNativeRun(set, workers, allowedPerWorker, shared + workers * allowedPerWorker); // 261110Cl (codex2 第 32 回): 方向ごとに実際の大きさで予算を確かめる。見込みを登録
        }
        int caps = NativeWrapper.DressedTailCapabilities;
        EbsdBeamSelectionNative = $"{NativeWrapper.LoadedNativeLibrary} caps {caps}{(set.AllowUBf16 && (caps & 2) != 0 ? " U-bf16" : "")}";
        // EbsdBeamSelectionSettingsVersion = DressedTailProductionSettings.Version + (set == DressedTailDefaults ? "" : "+custom"); // 261110Cl 変更前
        // EbsdBeamSelectionSettingsVersion = DressedTailProductionSettings.Version + (set == DressedTailDefaults ? "" : "+custom-" + set.Fingerprint()); // 261110Cl (codex2 第 34 回): 既定と違う設定どうしも区別する // 261110Cl 変更前
        EbsdBeamSelectionSettingsVersion = DressedTailProductionSettings.Version + (set == DressedTailDefaults ? "" : set == DressedTailLight ? "-light" : "+custom-" + set.Fingerprint()); // 261110Cl: 軽量の設定には名前
        var opt = new EbsdDressedTailOptions
        {
            TailRadiusNm = set.TailRadiusNm, PoolRadiusNm = set.PoolRadiusNm, Eta = set.Eta, N = set.N, Ps = set.Ps, RegDepthNm = set.RegDepthNm,
            OneProduct = true, BlockSize = 1, SourceTailQuadrature = set.SourceQuadrature, MaxParallelBlocks = workers, BuildCsr = false,
            SingleSolver = run.Solve, VoltageFinished = run.VoltageFinished,
        };
        return (opt, run, null);
    }

    [StructLayout(LayoutKind.Sequential)]
    private struct MemoryStatusEx { public uint Length, MemoryLoad; public ulong TotalPhys, AvailPhys, TotalPageFile, AvailPageFile, TotalVirtual, AvailVirtual, AvailExtendedVirtual; }
    [DllImport("kernel32.dll", SetLastError = true)]
    [return: MarshalAs(UnmanagedType.Bool)]
    private static extern bool GlobalMemoryStatusEx(ref MemoryStatusEx status);

    /// <summary>261110Cl 追加: 物理メモリの全体と、今使える量 (待機中のキャッシュを含む) [byte]。取れなければ GC の値</summary>
    private static (long total, long free) PhysicalMemory()
    {
        try
        {
            var m = new MemoryStatusEx { Length = (uint)Marshal.SizeOf<MemoryStatusEx>() };
            if (GlobalMemoryStatusEx(ref m)) return ((long)m.TotalPhys, (long)m.AvailPhys);
        }
        catch (Exception) { }
        var info = GC.GetGCMemoryInfo();
        return (info.TotalAvailableMemoryBytes, Math.Max(0, info.TotalAvailableMemoryBytes - info.MemoryLoadBytes));
    }

    /// <summary>261110Cl 追加 (codex2 第 31 回の 4): 新方式のメモリの見積もり (共有の表、1 作業あたり) [byte]。逆格子の点の密度 ρ = V_cell / (格子の心の多重度) から、
    /// 候補 (半径 R_pool の球)・裾 (芯の球の 1.4 倍)・繰り上げ (Ewald 球の厚さ 2ps の殻、π R_pool² × 2ps) を数える。消滅則は心の分だけ (保守的)</summary>
    private (long shared, long perWorker) EstimateDressedTailMemory(DressedTailProductionSettings set)
    {
        double v = Crystal.Volume; // [nm³]
        string lt = Crystal.Symmetry.LatticeTypeStr ?? "P";
        double mult = lt.StartsWith('F') ? 4 : lt.StartsWith('I') || lt.StartsWith('A') || lt.StartsWith('B') || lt.StartsWith('C') ? 2
            : lt.StartsWith('R') && Math.Abs(Crystal.Gamma - 2 * Math.PI / 3) < 0.02 ? 3 : 1;
        double rho = v / mult;
        int atoms = Crystal.Atoms.Sum(a => a.Atom.Count());
        double pool = 4.0 / 3 * Math.PI * Math.Pow(set.PoolRadiusNm, 3) * rho;
        double sites = Math.Max(1, atoms) * set.SourceQuadrature;
        double box = (4 * set.PoolRadiusNm * Crystal.A + 1) * (4 * set.PoolRadiusNm * Crystal.B + 1) * (4 * set.PoolRadiusNm * Crystal.C + 1);
        double shared = pool * sites * 24 + box * 48 + pool * 200; // 位相の表 (倍精度 + 単精度の写し)、差の表 (U・F、倍精度 + 単精度)、候補の配列
        // double nS = set.N + 2 * Math.PI * set.PoolRadiusNm * set.PoolRadiusNm * set.Ps * rho; // 261110Cl 変更前
        // double nQ = 1.4 * 4.0 / 3 * Math.PI * Math.Pow(set.CoreRadiusNm, 3) * rho; // 261110Cl 変更前
        // 261110Cl 変更: 明示・裾の数は候補の数 (前方の全部 ≤ pool) を超えない (繰り上げの閾値を大きくした照合の設定で桁があふれた)
        double nS = Math.Min(pool, set.N + 2 * Math.PI * set.PoolRadiusNm * set.PoolRadiusNm * set.Ps * rho);
        double nQ = Math.Min(pool, 1.4 * 4.0 / 3 * Math.PI * Math.Pow(set.CoreRadiusNm, 3) * rho);
        // double nM = 0.5 * nS; // 261110Cl 変更前
        double nM = nS; // 261110Cl 変更 (codex2 第 33 回): モードの数の上限 (τ = 0 なら nM = nS)。方向ごとの検査と同じ式
        double perWorker = 2 * nS * nQ * 8 + 3 * nQ * nM * 8 + nQ * nM * 4 + 10 * nS * nS * 16 + nQ * 64 + nS * nQ * 16; // H_SQ・X₀、T・EX₀・Y、X の bf16、密な小行列、裾の配列、C# のブロック
        return ((long)Math.Min(long.MaxValue / 4, shared), (long)Math.Min(long.MaxValue / 4, perWorker));
    }

    /// <summary>261110Cl 追加: 新方式の 1 回の構築の native の資源 (エネルギーごとの表・作業ごとの作業領域) の持ち主。構築の終わりに Dispose で全部返す</summary>
    private sealed class DressedTailNativeRun : IDisposable
    {
        readonly DressedTailProductionSettings set;
        readonly NativeWrapper.DressedTailNativeSettings nset;
        readonly DressedTailSelectionSettings sel;
        readonly ConcurrentBag<NativeWrapper.DressedTailWorkspaceHandle> idle = new();
        readonly List<NativeWrapper.DressedTailWorkspaceHandle> all = [];
        readonly object gate = new();
        readonly long perWorkerBudget; // 261110Cl 追加 (codex2 第 32 回): 1 作業が使ってよい量 [byte]
        readonly long reservedBytes; // 261110Cl 追加: この構築が s_dressedTailReserved に登録した見込み [byte] (Dispose で戻す)
        readonly int workerCount; // 261110Cl 追加: 作業の数 (方向ごとの検査で、そのときの空きを作業で割るため)
        DressedTailTable curTable; NativeWrapper.DressedTailTablesHandle curNative; int[] curPoolCol;
        bool disposed;

        // public DressedTailNativeRun(DressedTailProductionSettings s, int workers) // 261110Cl 旧シグネチャ
        // public DressedTailNativeRun(DressedTailProductionSettings s, int workers, long perWorkerBudget) // 261110Cl 旧シグネチャ
        public DressedTailNativeRun(DressedTailProductionSettings s, int workers, long perWorkerBudget, long reservedBytes)
        {
            set = s; this.perWorkerBudget = perWorkerBudget; workerCount = Math.Max(1, workers);
            lock (s_dressedTailBudgetGate) { s_dressedTailReserved += reservedBytes; this.reservedBytes = reservedBytes; } // 261110Cl 追加
            nset = new NativeWrapper.DressedTailNativeSettings
            {
                AbiVersion = NativeWrapper.DressedTailAbiVersionExpected, AllowUBf16 = s.AllowUBf16 ? 1 : 0,
                ModeTau = s.ModeTau, DegTol = s.DegTol, EdeltaThreshold = s.EdeltaThreshold, RegDepthNm = s.RegDepthNm,
            };
            sel = new DressedTailSelectionSettings { CoreRadiusNm = s.CoreRadiusNm, NeighborRadiusNm = s.NeighborRadiusNm, NeighborPs = s.NeighborPs, SkipQBlocks = true, RecordDiagnostics = false };
            try
            {
                for (int i = 0; i < workers; i++) { var w = NativeWrapper.DressedTailCreateWorkspace(); all.Add(w); idle.Add(w); }
            }
            catch { Dispose(); throw; }
        }

        (NativeWrapper.DressedTailTablesHandle tables, int[] poolCol) GetTables(DressedTailTable table, (double x, double y, double z)[] atomSites, (double weight, double sigmaNm)[] quadAtoms)
        {
            lock (gate)
            {
                ObjectDisposedException.ThrowIf(disposed, this);
                if (ReferenceEquals(curTable, table) && curNative != null) return (curNative, curPoolCol);
                curNative?.Dispose(); curNative = null; curTable = null;
                var (dU, dF, dh, dk, dl) = DressedTailDiffTablesResearch(table);
                var (dHkl, re, im) = DressedTailStencilResearch(table);
                int h = 0, k = 0, l = 0;
                foreach (var (a, b, c) in table.Pool.Hkl) { h = Math.Max(h, Math.Abs(a)); k = Math.Max(k, Math.Abs(b)); l = Math.Max(l, Math.Abs(c)); }
                var (phase, poolCol, w, nSites) = DressedTailQuadTableResearch(table, atomSites, quadAtoms, set.SourceQuadrature);
                curNative = NativeWrapper.DressedTailCreateTables(dU, dF, dh, dk, dl, dHkl, re, im, h, k, l, phase, nSites, w);
                curTable = table; curPoolCol = poolCol;
                return (curNative, curPoolCol);
            }
        }

        /// <summary>1 方向を解く (EbsdDressedTailOptions.SingleSolver)。並列で呼ばれる</summary>
        public DressedTailSolution Solve(BetheMethod bm, DressedTailTable table, Vector3DBase k0, Vector3DBase n, (double x, double y, double z)[] atomSites, (double weight, double sigmaNm)[] quadAtoms, bool tds)
        {
            if (!tds || quadAtoms == null) throw new NotSupportedException("The dressed-tail native solve supports only the TDS backscatter source");
            var (tables, poolCol) = GetTables(table, atomSites, quadAtoms);
            var b = bm.BuildDressedTailBlocksResearch(table, k0, n, set.N, set.Ps, sel);
            int nS = b.SArr.Length, nQ = b.QArr.Length; var hkl = table.Pool.Hkl;
            // 261110Cl 追加 (codex2 第 32 回): 実際の大きさで大きな確保の前に予算を確かめる (式は EstimateDressedTailMemory と同じ・モードの数は上限の nM = nS で見込む)。超えたら失敗として止める (途中で従来へ戻さない)
            double nMe = nS, need = 2.0 * nS * nQ * 8 + 3.0 * nQ * nMe * 8 + nQ * nMe * 4 + 10.0 * nS * nS * 16 + nQ * 64.0 + (double)nS * nQ * 16;
            // if (need > perWorkerBudget) // 261110Cl 変更前
            // 261110Cl 変更: 見積もりは裾を芯の球で数えるので、繰り上げた高次の明示のまわりの近傍が多い設定 (芯を小さくした軽い設定など) では実際の裾が数倍になる
            //   (Si・N64・芯 25 で見積もりの約 4 倍)。許す量を超えても、そのときの物理メモリの空き (予約を引いて作業で割った分) に収まれば進める。収まらなければ止める
            if (need > perWorkerBudget && need > (double)(PhysicalMemory().free - Math.Max(2L << 30, (long)(0.15 * PhysicalMemory().total))) / workerCount)
                throw new InvalidOperationException($"Dressed-tail beam selection: a direction needs {need / 1e9:f2} GB (explicit {nS}, tail {nQ}), more than the per-worker budget {perWorkerBudget / 1e9:f2} GB. Build again with fewer threads or with the conventional beam selection.");
            var sHkl = new int[3 * nS]; var qHkl = new int[3 * nQ]; var ds = new double[nS]; var dq = new double[nQ]; var qCol = new int[nQ];
            for (int i = 0; i < nS; i++) { var (h, k, l) = hkl[b.SArr[i]]; sHkl[3 * i] = h; sHkl[3 * i + 1] = k; sHkl[3 * i + 2] = l; ds[i] = b.D[b.SArr[i]]; }
            for (int r = 0; r < nQ; r++) { var (h, k, l) = hkl[b.QArr[r]]; qHkl[3 * r] = h; qHkl[3 * r + 1] = k; qHkl[3 * r + 2] = l; dq[r] = b.D[b.QArr[r]]; qCol[r] = poolCol[b.QArr[r]]; }
            // 261110Cl 変更 (codex2 第 32 回): 作業の数より多く並列に呼ばれたら、予算を破って足さずに止める (不変条件の違反)
            if (!idle.TryTake(out var ws))
                throw new InvalidOperationException("Dressed-tail beam selection: no idle workspace (more parallel solves than workers)");
            try
            {
                var (gamma, source) = NativeWrapper.DressedTailSolve(tables, ws, nset, sHkl, qHkl, ds, dq, b.HSS.Values, b.LamQ, qCol, Math.Sqrt(b.P[0]));
                return new DressedTailSolution { Pool = table.Pool, Gamma = gamma, SourceMatrix = source, NExplicit = nS, NSelected = b.NSelected, PoolIndices = new int[nS + nQ] };
            }
            finally { idle.Add(ws); }
        }

        /// <summary>1 エネルギーの解が終わったら native の表を返す (EbsdDressedTailOptions.VoltageFinished)</summary>
        public void VoltageFinished(DressedTailTable table)
        {
            lock (gate)
                if (ReferenceEquals(curTable, table)) { curNative?.Dispose(); curNative = null; curTable = null; curPoolCol = null; }
        }

        public void Dispose()
        {
            lock (gate)
            {
                if (disposed) return;
                disposed = true;
                curNative?.Dispose(); curNative = null; curTable = null;
                foreach (var w in all) w.Dispose();
                all.Clear();
                lock (s_dressedTailBudgetGate) s_dressedTailReserved -= reservedBytes; // 261110Cl 追加: 見込みを戻す (disposed で 1 回だけ)
            }
        }
    }
}
