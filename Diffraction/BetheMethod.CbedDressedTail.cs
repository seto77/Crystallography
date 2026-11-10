// 261110Cl 新規作成: CBED のビーム選択の新方式 (着衣つきの裾)。作者の判断 (2026-11-10): 選択式で本体へ入れる、既定は従来、区画の半径は 1 mrad から。
//   研究 = .project-guidance/ReciPro/ReciPro_着衣の裾_高速化_研究計画_261110.md §5.4・§5.9・§5.15〜 (今の CBED N64 より 100〜600 倍正確)。
//   円盤の方向を区画 (正方形、中心から角まで PatchRadiusMrad) に分け、区画の中心の方向で 1 回だけ選択と着衣 (DressCore) と基底の行列 (Cm・G₀・G_a) を作り、
//   区画の各方向は K = Cm + iU'₀G₀ − 2Σ k0_a G_a、M = 2(n·k0)G₀ + 2Σ n_a G_a (Q・P が k0・n の一次式なので厳密) の A = M⁻¹K を厚さへ伝播する
//   (SolveDressedTailBlockResearch と同じ式 = 10 月の WCphys)。
//   表面の法線 n は円盤の全方向で同じなので、M の方向による違いはスカラー n·k0 だけ。区画ごとに一般化固有分解 G_n v = μ G₀ v (V†G₀V = I、V†G_nV = diag μ) を 1 回解けば
//   M⁻¹ = ½ V diag(1/(n·k0 + μ)) V† となり、方向ごとの Cholesky と解き直しが要らない (基底 V の上で A' = ½ diag(1/(n·k0 + μ)) V†KV、厳密に同じ)。
//   出射の振幅は物理振幅 ψ = W c (W = D_c [I; X]) を、従来と同じ出射面の位相 exp(πi(P_g − 2k Surface.Z)t) で返す (位相の規約は従来と 8e−12 で一致、day11a)。
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Numerics;
using System.Threading;
using System.Threading.Tasks;
using MathNet.Numerics.LinearAlgebra;
using CDense = MathNet.Numerics.LinearAlgebra.Complex.DenseMatrix;

namespace Crystallography;

public partial class BetheMethod
{
    /// <summary>261110Cl 追加: CBED のビーム選択の方式</summary>
    public enum CbedBeamSelectionMethod
    {
        /// <summary>従来 (中心で Find_gVectors の N 本を選んで固定、弱いビームは捨てる)</summary>
        Legacy = 0,
        /// <summary>着衣つきの裾 (区画ごとに中心の着衣を使い回す)</summary>
        DressedTail = 1,
    }

    /// <summary>261110Cl 追加: CBED で要求する方式 (既定 = 従来)。新方式が使えない条件なら従来へ戻し、理由を <see cref="CbedBeamSelectionFallbackReason"/> に残す</summary>
    public CbedBeamSelectionMethod CbedBeamSelection { get; set; } = CbedBeamSelectionMethod.Legacy;
    /// <summary>261110Cl 追加: 直前の CBED で実際に使った方式</summary>
    public CbedBeamSelectionMethod CbedBeamSelectionUsed { get; private set; } = CbedBeamSelectionMethod.Legacy;
    /// <summary>261110Cl 追加: 新方式を要求したのに従来へ戻したときの理由 (戻さなければ null)</summary>
    public string CbedBeamSelectionFallbackReason { get; private set; }
    /// <summary>261110Cl 追加: 直前の新方式の CBED の内訳 (候補・区画・明示と裾の本数・時間) の文 (表示と記録用)</summary>
    public string CbedBeamSelectionDetail { get; private set; } = "";
    /// <summary>261110Cl 追加: 直前の新方式の CBED の段ごとの時間 [s] (0 表、1 出射のビームの選択、2 着衣 (CPU 秒)、3 区画の基底 (CPU 秒)、
    /// 4 方向ごとの行列の組み立て (CPU 秒)、5 厚さへの伝播 (固有分解など、CPU 秒)、6 出射の振幅と位相 (CPU 秒)、7 全体 (実時間))</summary>
    public double[] CbedDressedTailStageSeconds { get; } = new double[8];

    /// <summary>261110Cl 追加 (研究用): メモリの予算を上書きする [byte] (0 以下 = 物理メモリの空きから決める。込み合った機械での検査用)</summary>
    public long CbedMemoryBudgetOverride { get; set; }

    /// <summary>261110Cl 追加: CBED の新方式の設定。値で比べる</summary>
    public sealed record CbedDressedTailSettings
    {
        /// <summary>設定の版</summary>
        public const string Version = "CDT1-261110";
        /// <summary>E1 の基準の本数 (000 を含む。遅い反射の繰り上げで明示の本数は増える)</summary>
        public int N { get; init; } = 300;
        /// <summary>候補と裾の半径 [nm⁻¹] (研究の Si [110] 200 kV は FOLZ を含めて 60)</summary>
        public double PoolRadiusNm { get; init; } = 60;
        /// <summary>裾どうしの結合の閾値 (|U| ≥ η max)</summary>
        public double Eta { get; init; } = 0.05;
        /// <summary>遅い反射の繰り上げの閾値 |Q/P| [nm⁻¹]</summary>
        public double Ps { get; init; } = 0.1;
        /// <summary>着衣の正則化 δ = 1/(πT) の T [nm]</summary>
        public double RegDepthNm { get; init; } = 45;
        /// <summary>区画の半径 (中心から角まで) [mrad]。研究 (§5.4): 1 mrad で Si [111] 100 kV の局所解との差 1.5e−4、[110] 200 kV で 4e−6</summary>
        public double PatchRadiusMrad { get; init; } = 1;
        /// <summary>裾の芯の半径 [nm⁻¹] (NaN = 候補の全部を裾に)。近傍の半径と組で使う (EBSD の本番と同じ形)</summary>
        public double CoreRadiusNm { get; init; } = double.NaN;
        /// <summary>芯の外の裾 = 明示 (と |Q/P| &lt; NeighborPs の候補) からこの距離 [nm⁻¹] の中 (NaN = 芯の外も全部)</summary>
        public double NeighborRadiusNm { get; init; } = double.NaN;
        /// <summary>近傍の中心に |Q/P| &lt; これの候補も (NaN = 明示だけ)</summary>
        public double NeighborPs { get; init; } = double.NaN;
        /// <summary>射影の E·X を着衣前の E·X₀ で近似し、疎な積を 1 回にする (EBSD の本番と同じ)</summary>
        public bool OneProduct { get; init; }
        /// <summary>着衣のモードの組の打ち切り τ (区画の中心の入射の重みの最大に対する比。0 = 打ち切らない)。EBSD の本番は 3e−3</summary>
        public double ModeTau { get; init; }
        /// <summary>適応的な区画の最大の半径 [mrad] (NaN = 使わない = 一様な区画 PatchRadiusMrad)。使うときは、この半径の区画から始め、
        /// 目安 η (<see cref="CbedPatchIndicator"/>) が AdaptiveIndicatorLimit を超えたら 4 分割する (最小は PatchRadiusMrad)</summary>
        public double AdaptiveMaxRadiusMrad { get; init; } = double.NaN;
        /// <summary>適応的な区画の目安 η の上限</summary>
        public double AdaptiveIndicatorLimit { get; init; } = double.NaN;
        /// <summary>区画の中のモードの打ち切り (NaN = 使わない)。区画の中心の A' の固有モード (右 R・左 L = R⁻¹) のうち、区画の中心と端の探りの方向 (最大 8) で
        /// 厚さのどこかの振幅 |(L y(t))_j|·‖R_j‖ が最大のこの比以上のモード (近い固有値の組は一緒に) だけを残し、各方向は L_k A'(d) R_k (次元 nM) で解く</summary>
        public double ModeReduceTau { get; init; } = double.NaN;
        /// <summary>区画の基底の行列のうち裾 (A_Q) の積を単精度 (complex64) で積む (明示の部分は倍精度)。裾の寄与は小さいので丸めの影響は小さい見込み</summary>
        public bool SingleTailProducts { get; init; }
    }
    /// <summary>261110Cl 追加: CBED の新方式の既定の設定</summary>
    public static CbedDressedTailSettings CbedDressedTailDefaults { get; } = new();
    /// <summary>261110Cl 追加: この BetheMethod の CBED で使う新方式の設定</summary>
    public CbedDressedTailSettings CbedDressedTail { get; set; } = CbedDressedTailDefaults;

    /// <summary>261110Cl 追加: CBED の新方式で円盤の全方向を解く。返り値は出射のビーム (000 が先頭) と、方向ごとの振幅 [方向][t × nOut + b] (従来の cbed_DoWork と同じ並び・位相の規約)。
    /// 使えない条件 (非局所の吸収・メモリの予算) なら beams = null を返し、理由を <see cref="CbedBeamSelectionFallbackReason"/> に入れる。
    /// 出射のビーム = 各区画の中心の明示 (E1 + 繰り上げ) の和集合。区画の中でそのビームが明示でも裾でもなければ (後ろ向きなど) 振幅 0。</summary>
    /// <param name="kvac">真空の波数 [nm⁻¹]</param>
    /// <param name="u0">U₀ の実部 (屈折)</param>
    /// <param name="directions">入射の方向 (単位ベクトル、実験室の座標)</param>
    /// <param name="thicknesses">厚さ [nm]</param>
    /// <param name="threads">並列の数</param>
    /// <param name="solver">厚さへの伝播の解き手 (従来の CBED と同じ選択。MtxExp は等間隔の厚さのときだけ、Eigen の native が無ければ MathNet の固有分解)</param>
    /// <param name="cancelled">中止の問い合わせ (true で途中の方向を null のまま返す)</param>
    /// <param name="progressed">終わった方向の累計を知らせる</param>
    public (Beam[] beams, Complex[][] amplitudes) SolveCbedDressedTail(double kvac, double u0, IReadOnlyList<Vector3DBase> directions, double[] thicknesses, int threads,
        Solver solver = Solver.Eigen_Eigen, Func<bool> cancelled = null, Action<int> progressed = null)
    {
        var set = CbedDressedTail ?? CbedDressedTailDefaults;
        CbedBeamSelectionFallbackReason = null;
        Array.Clear(CbedDressedTailStageSeconds);
        threads = Math.Max(1, threads);
        if (UseNonLocalAbsorption) { CbedBeamSelectionFallbackReason = "non-local absorption is not supported by the dressed-tail method"; return (null, null); }
        var swAll = Stopwatch.StartNew();
        // 厚さへの伝播の形: 0 = native の固有分解、1 = native の行列指数 (等間隔の厚さ)、2 = MathNet の固有分解
        int prop = !EigenEnabled ? 2 : solver == Solver.Eigen_Eigen ? 0 : solver == Solver.MtxExp_Eigen && NativeWrapper.IsUniformThickness(thicknesses) ? 1
            : solver == Solver.MtxExp_Eigen ? 0 : 2;

        // メモリの見積もり (1 区画): 着衣と基底の nQ × nS の密な行列 約 12 本。候補の数は逆格子の点の密度から (EstimateDressedTailMemory と同じ数え方)。
        //   CBED は k0 が大きく候補はほぼ全部が前方 (P > 0) なので、裾 ≈ 候補 (芯と近傍を使うときは芯の球の 1.4 倍)
        string lt = Crystal.Symmetry.LatticeTypeStr ?? "P";
        double mult = lt.StartsWith('F') ? 4 : lt.StartsWith('I') || lt.StartsWith('A') || lt.StartsWith('B') || lt.StartsWith('C') ? 2
            : lt.StartsWith('R') && Math.Abs(Crystal.Gamma - 2 * Math.PI / 3) < 0.02 ? 3 : 1;
        double poolEst = 4.0 / 3 * Math.PI * Math.Pow(set.PoolRadiusNm, 3) * Crystal.Volume / mult;
        double nSEst = Math.Min(poolEst, set.N + 150);
        double nQEst = double.IsNaN(set.CoreRadiusNm) || double.IsNaN(set.NeighborRadiusNm) ? poolEst : Math.Min(poolEst, 1.4 * 4.0 / 3 * Math.PI * Math.Pow(set.CoreRadiusNm, 3) * Crystal.Volume / mult);
        long perPatch = (long)Math.Min(long.MaxValue / 4, 12 * nQEst * nSEst * 16);
        var (total, free) = PhysicalMemory();
        long reserve = Math.Max(2L << 30, (long)(0.15 * total));
        long budget = CbedMemoryBudgetOverride > 0 ? CbedMemoryBudgetOverride : Math.Min(total / 2, free - reserve);
        int byBudget = (int)Math.Min(int.MaxValue, Math.Max(0, budget / Math.Max(1, perPatch)));
        if (byBudget < 1)
        {
            CbedBeamSelectionFallbackReason = $"estimated memory {perPatch / 1e9:f1} GB per patch (about {poolEst:f0} candidates within {set.PoolRadiusNm} nm⁻¹) exceeds the budget {budget / 1e9:f1} GB";
            return (null, null);
        }

        // 候補のプールと結合の表 (方向に依らない)
        DressedTailPrecomputeOptions = new ParallelOptions { MaxDegreeOfParallelism = threads };
        var pool = CreateBeamSelectionPool(set.PoolRadiusNm, buildE0: false);
        var table = CreateDressedTailTable(pool, set.Eta, buildCsr: true);
        CbedDressedTailStageSeconds[0] = swAll.Elapsed.TotalSeconds;
        var sel = new DressedTailSelectionSettings
        {
            CoreRadiusNm = double.IsNaN(set.CoreRadiusNm) ? set.PoolRadiusNm : set.CoreRadiusNm, NeighborRadiusNm = set.NeighborRadiusNm, NeighborPs = set.NeighborPs, RecordDiagnostics = false
        };
        var n = Vector3DBase.Normarize(Surface);
        var (nx, ny, nz) = n.Tuple;
        double surfZ = Surface.Z; // 従来の出射面の位相と同じ値 (正規化しない)

        // 区画: 円盤の中心の方向 c の接平面の座標 (u, v) = (d·e1, d·e2)/(d·c) を一辺 √2 r の升に分ける (c を升の中心に)
        int nDir = directions.Count;
        var sum = new Vector3DBase(0, 0, 0);
        foreach (var d in directions) sum += d;
        var c = Vector3DBase.Normarize(sum);
        var e1 = Vector3DBase.Normarize(Vector3DBase.VectorProduct(c, Math.Abs(c.Z) < 0.9 ? new Vector3DBase(0, 0, 1) : new Vector3DBase(1, 0, 0)));
        var e2 = Vector3DBase.Normarize(Vector3DBase.VectorProduct(c, e1));
        var uu = new double[nDir]; var vv = new double[nDir];
        for (int i = 0; i < nDir; i++) { var d = directions[i]; double dc = d * c; uu[i] = (d * e1) / dc; vv[i] = (d * e2) / dc; }
        // 升 (中心 (uc, vc)・一辺 side・方向の番号)。原点を升の中心にした格子で分ける
        List<(double uc, double vc, double side, int[] members)> Squares(IEnumerable<int> idx, double sd)
        {
            var g = new Dictionary<(long, long), List<int>>();
            foreach (var i in idx)
            {
                var key = ((long)Math.Floor(uu[i] / sd + 0.5), (long)Math.Floor(vv[i] / sd + 0.5));
                if (!g.TryGetValue(key, out var list)) g[key] = list = [];
                list.Add(i);
            }
            return g.OrderBy(kv => kv.Key).Select(kv => ((double)kv.Key.Item1 * sd, (double)kv.Key.Item2 * sd, sd, kv.Value.ToArray())).ToList();
        }
        // 升を 2^L × 2^L に分ける
        List<(double uc, double vc, double side, int[] members)> Split((double uc, double vc, double side, int[] members) sq, int L)
        {
            int m = 1 << L; double cs = sq.side / m, u0c = sq.uc - sq.side / 2, v0c = sq.vc - sq.side / 2;
            var g = new Dictionary<(int, int), List<int>>();
            foreach (var i in sq.members)
            {
                var key = (Math.Clamp((int)Math.Floor((uu[i] - u0c) / cs), 0, m - 1), Math.Clamp((int)Math.Floor((vv[i] - v0c) / cs), 0, m - 1));
                if (!g.TryGetValue(key, out var list)) g[key] = list = [];
                list.Add(i);
            }
            return g.OrderBy(kv => kv.Key).Select(kv => (u0c + (kv.Key.Item1 + 0.5) * cs, v0c + (kv.Key.Item2 + 0.5) * cs, cs, kv.Value.ToArray())).ToList();
        }
        Vector3DBase CenterOf(int[] members) { var s = new Vector3DBase(0, 0, 0); foreach (var i in members) s += directions[i]; return Vector3DBase.Normarize(s); }
        double sideMin = Math.Sqrt(2) * set.PatchRadiusMrad * 1e-3;
        bool adaptive = double.IsFinite(set.AdaptiveMaxRadiusMrad) && double.IsFinite(set.AdaptiveIndicatorLimit) && set.AdaptiveMaxRadiusMrad > set.PatchRadiusMrad;
        int maxLevel = adaptive ? (int)Math.Max(0, Math.Round(Math.Log2(set.AdaptiveMaxRadiusMrad / set.PatchRadiusMrad))) : 0;
        double side0 = sideMin * (1 << maxLevel);
        // 出射のビームの選択は最小の升の中心で (適応的な区画でも、細かく分けた所の明示を落とさない)
        var patches = Squares(Enumerable.Range(0, nDir), sideMin).Select(q => (center: CenterOf(q.members), q.members)).ToArray();

        // 出射のビーム = 各区画の中心の明示の和集合 (選択だけ。000 を先頭、円盤の中心に近い区画の順)
        var sw = Stopwatch.StartNew();
        var outPool = new List<int> { 0 };
        var inOut = new bool[pool.Count]; inOut[0] = true;
        var countSel = new DressedTailSelectionSettings { CoreRadiusNm = sel.CoreRadiusNm, NeighborRadiusNm = sel.NeighborRadiusNm, NeighborPs = sel.NeighborPs, RecordDiagnostics = false, CountOnly = true };
        var patchOrder = Enumerable.Range(0, patches.Length).OrderBy(p => -(patches[p].center * c)).ToArray();
        var selections = new int[patches.Length][];
        Parallel.For(0, patches.Length, new ParallelOptions { MaxDegreeOfParallelism = threads }, p =>
        {
            var ck = getVecK0(kvac, u0, patches[p].center, n);
            selections[p] = DressCore(table, ck, n, set.N, set.Ps, set.RegDepthNm, set.OneProduct, false, null, selection: countSel).SArr;
        });
        foreach (var p in patchOrder)
            foreach (var i in selections[p])
                if (!inOut[i]) { inOut[i] = true; outPool.Add(i); }
        int nOut = outPool.Count, tLen = thicknesses.Length;
        var gOut = outPool.Select(i => (pool.GX[i], pool.GY[i], pool.GZ[i])).ToArray();
        var cK0 = getVecK0(kvac, u0, c, n);
        var beams = new Beam[nOut];
        for (int o = 0; o < nOut; o++)
        {
            var g = new Vector3DBase(gOut[o].Item1, gOut[o].Item2, gOut[o].Item3);
            beams[o] = new Beam(pool.Hkl[outPool[o]], g, getU(AccVoltage, pool.Hkl[outPool[o]], g), getQP(g, cK0, n));
        }
        CbedDressedTailStageSeconds[1] = sw.Elapsed.TotalSeconds;

        // 並列: 区画を外側で (メモリの予算の内)、区画の中の方向を内側で
        int outer = Math.Max(1, Math.Min(Math.Min(threads, patches.Length), byBudget));
        int inner = Math.Max(1, threads / outer);
        var amplitudes = new Complex[nDir][];
        var iu0 = Complex.ImaginaryOne * getU(AccVoltage, new Beam((0, 0, 0), new Vector3DBase(0, 0, 0)), cache: false).Imag;
        int done = 0, nPatchesDone = 0;
        long sumS = 0, sumQ = 0, sumB = 0; int minS = int.MaxValue, maxS = 0, nBackward = 0; double etaMin = double.PositiveInfinity, etaMax = 0;
        var patchEta = new List<(Vector3DBase center, double eta)>();
        int nSplit = 0, nSolvedPatches = 0;
        var gate = new object();
        var backward = new System.Collections.Concurrent.ConcurrentBag<int>();
        var stage = CbedDressedTailStageSeconds;
        // 出射面の位相の厚さの因子 (方向に依らない部分): exp(πi(P_g − 2k Surface.Z)t) を方向ごとに作る
        var twoPiI = 2 * Math.PI * Complex.ImaginaryOne;

        // MathNet の固有分解で y(t) を作る (native の無いとき)
        Complex[] ManagedPropagate(Complex[] a, Complex[] y0, int dim)
        {
            var evd = new CDense(dim, dim, a).Evd(MathNet.Numerics.LinearAlgebra.Symmetricity.Asymmetric);
            var alpha = evd.EigenVectors.Solve(MathNet.Numerics.LinearAlgebra.Vector<Complex>.Build.Dense(y0));
            var ph = new CDense(dim, tLen);
            for (int t = 0; t < tLen; t++)
                for (int j = 0; j < dim; j++)
                    ph[j, t] = Complex.Exp(twoPiI * evd.EigenValues[j] * thicknesses[t]) * alpha[j];
            return (evd.EigenVectors * ph).ToColumnMajorArray();
        }

        // 1 区画を解く (members の方向を center の着衣で)。後ろ向きのビームがある方向は backward へ返す。
        //   適応的な区画では、目安 η が上限を超え、まだ分けられるなら、分けた升を next へ返して解かない (η は区画の半径にほぼ比例するので、必要な段数を一度に分ける)
        void SolvePatch((double uc, double vc, double side, int[] members) sq, bool parallelInside, System.Collections.Concurrent.ConcurrentBag<(double, double, double, int[])> next)
        {
            if (cancelled?.Invoke() == true) return;
            var members = sq.members;
            var center = CenterOf(members);
            double rEff = 0; foreach (var i in members) rEff = Math.Max(rEff, Math.Acos(Math.Min(1, directions[i] * center))); // 中心から最も遠い方向までの角度 [rad]
            var swp = Stopwatch.StartNew();
            var ck = getVecK0(kvac, u0, center, n);
            // 全次元 (nB) の A'(d) と y(0) (基底の行列 bc, bx, by, bz, mu, v0 を作った後で使う)
            Complex[] bc = null, bx = null, by = null, bz = null, v0 = null; double[] mu = null; int nB = 0;
            (Complex[] a, Complex[] y0, double nk, Vector3DBase k0) FullSystem(int di)
            {
                var k0 = getVecK0(kvac, u0, directions[di], n);
                var (kx, ky, kz) = k0.Tuple;
                double nk = Dot3(nx, kx, ny, ky, nz, kz);
                var a = new Complex[nB * nB]; var y0 = new Complex[nB]; var sc = new double[nB];
                for (int i = 0; i < nB; i++) { sc[i] = 0.5 / (nk + mu[i]); y0[i] = 2 * nk * sc[i] * v0[i]; }
                for (int col = 0; col < nB; col++)
                {
                    int o = col * nB;
                    for (int i = 0; i < nB; i++) a[o + i] = (bc[o + i] - 2 * (kx * bx[o + i] + ky * by[o + i] + kz * bz[o + i])) * sc[i];
                }
                return (a, y0, nk, k0);
            }
            // 着衣はモードの基底 (A_S = V_b、A_Q = X̃、張る空間は [I; X] と同じ) のまま使う: X = X̃V_b⁻¹ の積と逆行列を省き、τ > 0 ならモードを打ち切る
            var cr = DressCore(table, ck, n, set.N, set.Ps, set.RegDepthNm, set.OneProduct, parallelInside, null, modeBasis: true, selection: sel, modeTau: set.ModeTau);
            double tDress = swp.Elapsed.TotalSeconds; swp.Restart();
            var (sArr, qArr, Dc) = (cr.SArr, cr.QArr, cr.D);
            int nS = sArr.Length, nQ = qArr.Length; nB = cr.Vb.ColumnCount;
            double eta = CbedPatchIndicator(cr, pool, n, kvac * rEff, set.RegDepthNm);
            if (next != null && adaptive && eta > set.AdaptiveIndicatorLimit && sq.side > sideMin * 1.001 && members.Length > 1)
            {
                int levels = Math.Clamp((int)Math.Ceiling(Math.Log2(eta / set.AdaptiveIndicatorLimit)), 1, (int)Math.Round(Math.Log2(sq.side / sideMin)));
                foreach (var ch in Split(sq, levels)) next.Add(ch);
                lock (gate) { stage[2] += tDress; nSplit++; }
                return;
            }
            // 出射のビームの行 W_out (nOut × nB): 明示なら D_c A_S の行、裾なら D_c A_Q の行
            var posS = new Dictionary<int, int>(nS); for (int i = 0; i < nS; i++) posS[sArr[i]] = i;
            var posQ = new Dictionary<int, int>(nQ); for (int r = 0; r < nQ; r++) posQ[qArr[r]] = r;
            var wout = new CDense(nOut, nB);
            for (int o = 0; o < nOut; o++)
            {
                int pi = outPool[o];
                if (posS.TryGetValue(pi, out var si)) for (int j = 0; j < nB; j++) wout[o, j] = Dc[pi] * cr.Vb[si, j];
                else if (posQ.TryGetValue(pi, out var r)) for (int j = 0; j < nB; j++) wout[o, j] = Dc[pi] * cr.Xt[r, j];
            }
            var row0 = cr.Vb.Row(0); // 000 (明示の先頭) の行
            var all = sArr.Concat(qArr).ToArray();
            var (Cm, G0, Gx, Gy, Gz) = CbedSharedBasis(cr, pool, set.OneProduct, set.SingleTailProducts);
            cr = null; // 大きな行列を手放す
            // 一般化固有分解 G_n v = μ G₀ v (G₀ = LL†、C = L⁻¹G_nL⁻† のエルミートの固有分解、V = L⁻†U)
            var Gn = nx * Gx + ny * Gy + nz * Gz;
            var Linv = G0.Cholesky().Factor.Inverse();
            var Ch = Linv * Gn * Linv.ConjugateTranspose();
            Ch = 0.5 * (Ch + Ch.ConjugateTranspose());
            var evdG = Ch.Evd(MathNet.Numerics.LinearAlgebra.Symmetricity.Hermitian);
            var V = Linv.ConjugateTranspose() * evdG.EigenVectors;
            mu = evdG.EigenValues.Select(z => z.Real).ToArray();
            var Vh = V.ConjugateTranspose();
            bc = ((CDense)(Vh * (Cm + iu0 * G0) * V)).Values;
            bx = ((CDense)(Vh * Gx * V)).Values; by = ((CDense)(Vh * Gy * V)).Values; bz = ((CDense)(Vh * Gz * V)).Values;
            var wp = (CDense)(wout * V); // 出射のビームの行を基底 V で
            var r0v = row0 * V; // W†e₀ = D_c,000 conj(A_S[0, :]) → (V†W†e₀)_j = D_c,000 conj((A_S V)[0, j])
            v0 = new Complex[nB]; for (int j = 0; j < nB; j++) v0[j] = Complex.Conjugate(r0v[j]) * Dc[0];
            Cm = G0 = Gx = Gy = Gz = Gn = Linv = Ch = V = Vh = null;
            // 区画の中のモードの打ち切り (2 側の射影)。A'(d) = S(d)B(d)、S(d)_i = ½/(n·k0 + μ_i) = S_c,i (1 − Δnk/(n·k0 + μ_i)) ≈ S_c,i (1 − Δnk w_i) (w_i = 1/(n·k_c + μ_i)、Δnk² の項は 1e−11 程度)
            //   → L_k A'(d) R_k = Σ_M c_M (P_M − Δnk Q_M)、P_M = L_k S_c M R_k、Q_M = L_k S_c W M R_k (M = B_c, B_x, B_y, B_z、c = 1, −2k_x, −2k_y, −2k_z)
            int nR = nB; CDense[] redP = null, redQ = null; Matrix<Complex> Lk = null; CDense wpR = null; double nkc = Dot3(nx, ck.X, ny, ck.Y, nz, ck.Z);
            var probeDone = new Dictionary<int, Complex[]>(); // 探りの方向は全次元で解いた結果をそのまま使う
            if (double.IsFinite(set.ModeReduceTau) && set.ModeReduceTau > 0 && members.Length > 9)
            {
                var (ckx, cky, ckz) = ck.Tuple;
                var scc = new double[nB]; for (int i = 0; i < nB; i++) scc[i] = 0.5 / (nkc + mu[i]);
                var acv = new Complex[nB * nB];
                for (int col = 0; col < nB; col++) { int o = col * nB; for (int i = 0; i < nB; i++) acv[o + i] = (bc[o + i] - 2 * (ckx * bx[o + i] + cky * by[o + i] + ckz * bz[o + i])) * scc[i]; }
                var evc = new CDense(nB, nB, acv).Evd(MathNet.Numerics.LinearAlgebra.Symmetricity.Asymmetric);
                var Rc = evc.EigenVectors; var gc = evc.EigenValues;
                var Lc = Rc.Inverse();
                var rnorm = new double[nB]; for (int j = 0; j < nB; j++) rnorm[j] = Rc.Column(j).L2Norm();
                // 探り: 区画の中心から 8 方位で最も遠い方向
                var probes = new HashSet<int>();
                for (int m8 = 0; m8 < 8; m8++)
                {
                    double ca = Math.Cos(m8 * Math.PI / 4), sa = Math.Sin(m8 * Math.PI / 4), best = double.NegativeInfinity; int bi = members[0];
                    foreach (var i in members) { double pr = (uu[i] - sq.uc) * ca + (vv[i] - sq.vc) * sa; if (pr > best) { best = pr; bi = i; } }
                    probes.Add(bi);
                }
                var wgt = new double[nB];
                foreach (var di in probes)
                {
                    var (a, y0, _, _) = FullSystem(di);
                    var yt = prop == 1 ? NativeWrapper.CBEDSolver_MatExp(a, y0, thicknesses) : NativeWrapper.CBEDSolver_Eigen(a, y0, thicknesses);
                    probeDone[di] = yt;
                    var z = Lc * new CDense(nB, tLen, yt);
                    for (int j = 0; j < nB; j++) for (int t = 0; t < tLen; t++) wgt[j] = Math.Max(wgt[j], z[j, t].Magnitude * rnorm[j]);
                }
                // 区画の中心 (方向の平均) の励起も入れる: z(t) = e^{2πiγt} L_c y0_c
                {
                    var y0c = new Complex[nB]; for (int i = 0; i < nB; i++) y0c[i] = 2 * nkc * scc[i] * v0[i];
                    var zc = Lc * MathNet.Numerics.LinearAlgebra.Vector<Complex>.Build.Dense(y0c);
                    for (int j = 0; j < nB; j++) wgt[j] = Math.Max(wgt[j], zc[j].Magnitude * rnorm[j] * Math.Exp(-2 * Math.PI * Math.Min(0, gc[j].Imaginary) * thicknesses.Max()));
                }
                double wmax = wgt.Max();
                // 近い固有値の組 (|Δγ| < 1e−3 nm⁻¹) は一緒に残す・落とす
                var order = Enumerable.Range(0, nB).OrderBy(j => gc[j].Real).ToArray();
                var keepJ = new bool[nB];
                for (int a0 = 0; a0 < nB;)
                {
                    int b0 = a0 + 1;
                    while (b0 < nB && (gc[order[b0]] - gc[order[b0 - 1]]).Magnitude < 1e-3) b0++;
                    bool k = false; for (int q = a0; q < b0; q++) if (wgt[order[q]] >= set.ModeReduceTau * wmax) k = true;
                    for (int q = a0; q < b0; q++) keepJ[order[q]] = k;
                    a0 = b0;
                }
                var keep = Enumerable.Range(0, nB).Where(j => keepJ[j]).ToArray();
                nR = keep.Length;
                var Rk = CDense.Create(nB, nR, (i, j) => Rc[i, keep[j]]);
                Lk = CDense.Create(nR, nB, (i, j) => Lc[keep[i], j]);
                var wv = new double[nB]; for (int i = 0; i < nB; i++) wv[i] = 1 / (nkc + mu[i]);
                redP = new CDense[4]; redQ = new CDense[4];
                var mats = new[] { bc, bx, by, bz };
                for (int mm = 0; mm < 4; mm++)
                {
                    var Ms = new CDense(nB, nB); var Mw = new CDense(nB, nB); var msv = Ms.Values; var mwv = Mw.Values; var src = mats[mm];
                    for (int col = 0; col < nB; col++) { int o = col * nB; for (int i = 0; i < nB; i++) { msv[o + i] = src[o + i] * scc[i]; mwv[o + i] = src[o + i] * (scc[i] * wv[i]); } }
                    redP[mm] = (CDense)(Lk * (Ms * Rk)); redQ[mm] = (CDense)(Lk * (Mw * Rk));
                }
                wpR = (CDense)(wp * Rk);
            }
            double tBasis = swp.Elapsed.TotalSeconds;
            double tForm = 0, tProp = 0, tOut = 0;
            Parallel.For(0, members.Length, new ParallelOptions { MaxDegreeOfParallelism = parallelInside ? inner : 1 }, () => (0.0, 0.0, 0.0), (m, _, acc) =>
            {
                if (cancelled?.Invoke() == true) return acc;
                var swd = Stopwatch.StartNew();
                int di = members[m];
                var k0 = getVecK0(kvac, u0, directions[di], n);
                var (kx, ky, kz) = k0.Tuple;
                foreach (var i in all)
                    if (2 * Dot3(nx, pool.GX[i] + kx, ny, pool.GY[i] + ky, nz, pool.GZ[i] + kz) <= 0) { backward.Add(di); return acc; }
                double t1;
                Complex[] yt; // y(t) (列優先)
                CDense wUse = wp; int nUse = nB;
                if (probeDone.TryGetValue(di, out var ypr)) { yt = ypr; t1 = swd.Elapsed.TotalSeconds; }
                else if (redP != null)
                {
                    // 縮約: A_red = Σ c_M (P_M − Δnk Q_M)、y0_red = L_k y0(d)
                    double nk = Dot3(nx, kx, ny, ky, nz, kz), dnk = nk - nkc;
                    var cM = new[] { 1.0, -2 * kx, -2 * ky, -2 * kz };
                    var ar = new Complex[nR * nR];
                    for (int mm = 0; mm < 4; mm++)
                    {
                        var pv = redP[mm].Values; var qv = redQ[mm].Values; double cm = cM[mm];
                        for (int i = 0; i < ar.Length; i++) ar[i] += cm * (pv[i] - dnk * qv[i]);
                    }
                    var y0f = new Complex[nB]; for (int i = 0; i < nB; i++) y0f[i] = nk / (nk + mu[i]) * v0[i];
                    var y0r = (Lk * MathNet.Numerics.LinearAlgebra.Vector<Complex>.Build.Dense(y0f)).ToArray();
                    t1 = swd.Elapsed.TotalSeconds;
                    yt = prop == 1 ? NativeWrapper.CBEDSolver_MatExp(ar, y0r, thicknesses) : prop == 0 ? NativeWrapper.CBEDSolver_Eigen(ar, y0r, thicknesses) : ManagedPropagate(ar, y0r, nR);
                    wUse = wpR; nUse = nR;
                }
                else
                {
                    var (a, y0, _, _) = FullSystem(di);
                    t1 = swd.Elapsed.TotalSeconds;
                    yt = prop == 0 ? NativeWrapper.CBEDSolver_Eigen(a, y0, thicknesses) : prop == 1 ? NativeWrapper.CBEDSolver_MatExp(a, y0, thicknesses) : ManagedPropagate(a, y0, nB);
                }
                double t2 = swd.Elapsed.TotalSeconds;
                var psi = (wUse * new CDense(nUse, tLen, yt)).ToColumnMajorArray(); // nOut × tLen
                // 出射面の位相 (従来と同じ: exp(πi(P_g − 2k Surface.Z)t)、P_g = 2n·(k0 + g))
                for (int o = 0; o < nOut; o++)
                {
                    var (gx, gy, gz) = gOut[o];
                    double p = 2 * Dot3(nx, gx + kx, ny, gy + ky, nz, gz + kz);
                    for (int t = 0; t < tLen; t++)
                        psi[o + t * nOut] *= Complex.Exp(Math.PI * Complex.ImaginaryOne * (p - 2 * kvac * surfZ) * thicknesses[t]);
                }
                amplitudes[di] = psi;
                int dn = Interlocked.Increment(ref done);
                progressed?.Invoke(dn);
                double t3 = swd.Elapsed.TotalSeconds;
                return (acc.Item1 + t1, acc.Item2 + t2 - t1, acc.Item3 + t3 - t2);
            }, acc => { lock (gate) { tForm += acc.Item1; tProp += acc.Item2; tOut += acc.Item3; } });
            lock (gate)
            {
                stage[2] += tDress; stage[3] += tBasis; stage[4] += tForm; stage[5] += tProp; stage[6] += tOut; nPatchesDone++; nSolvedPatches++;
                sumS += nS; sumQ += nQ; sumB += nR; minS = Math.Min(minS, nS); maxS = Math.Max(maxS, nS); etaMin = Math.Min(etaMin, eta); etaMax = Math.Max(etaMax, eta); patchEta.Add((center, eta));
            }
        }

        var work = adaptive ? Squares(Enumerable.Range(0, nDir), side0) : Squares(Enumerable.Range(0, nDir), sideMin);
        int nLevels = 0;
        while (work.Count > 0 && cancelled?.Invoke() != true)
        {
            var next = new System.Collections.Concurrent.ConcurrentBag<(double, double, double, int[])>();
            int o2 = Math.Max(1, Math.Min(outer, work.Count)), i2 = Math.Max(1, threads / o2);
            inner = i2;
            var cur = work;
            Parallel.For(0, cur.Count, new ParallelOptions { MaxDegreeOfParallelism = o2 }, p => SolvePatch(cur[p], o2 == 1 || i2 > 1, next));
            work = next.OrderBy(q => q.Item1).ThenBy(q => q.Item2).ToList();
            nLevels++;
        }
        // 区画の中心で後ろ向きのビームがあった方向は、自分を中心に解き直す (CBED の通常の幾何ではまず起きない)
        var back = backward.Distinct().ToArray();
        nBackward = back.Length;
        Parallel.For(0, back.Length, new ParallelOptions { MaxDegreeOfParallelism = threads }, b => SolvePatch((uu[back[b]], vv[back[b]], 0, new[] { back[b] }), false, null));
        stage[7] = swAll.Elapsed.TotalSeconds;

        string[] propName = ["native eigen", "native matrix exponential", "MathNet eigen"];
        CbedBeamSelectionDetail = $"dressed tail {CbedDressedTailSettings.Version}: N {set.N}, pool {pool.Count} within {set.PoolRadiusNm} nm⁻¹, eta {set.Eta}" +
            (double.IsNaN(set.NeighborRadiusNm) ? "" : $", core {set.CoreRadiusNm} + neighbours {set.NeighborRadiusNm} nm⁻¹") + (set.OneProduct ? ", one product" : "") + "; " +
            (adaptive ? $"adaptive patches (radius {set.AdaptiveMaxRadiusMrad} down to {set.PatchRadiusMrad} mrad, indicator limit {set.AdaptiveIndicatorLimit}): {nSolvedPatches} solved, {nSplit} split, {nLevels} rounds" : $"{nSolvedPatches} patches (radius {set.PatchRadiusMrad} mrad)") + $" for {nDir} directions, {nBackward} re-solved alone; explicit {minS}-{maxS} (mean {(nPatchesDone > 0 ? sumS / nPatchesDone : 0)}), tail mean {(nPatchesDone > 0 ? sumQ / nPatchesDone : 0)}, basis mean {(nPatchesDone > 0 ? sumB / nPatchesDone : 0)} (mode tau {set.ModeTau}{(double.IsFinite(set.ModeReduceTau) ? $", reduced by {set.ModeReduceTau}" : "")}), indicator {etaMin:e2}-{etaMax:e2}; " +
            $"output beams {nOut}; propagation {propName[prop]}; time: tables {stage[0]:f1} s, output-beam selection {stage[1]:f1} s, dressing {stage[2]:f1} CPU-s, basis {stage[3]:f1} CPU-s, " +
            $"per direction form {stage[4]:f1} / propagate {stage[5]:f1} / output {stage[6]:f1} CPU-s, total {stage[7]:f1} s; parallel {outer} patches × {inner}";
        return (beams, amplitudes);
    }

    /// <summary>261110Cl 追加: 区画の基底の行列。物理振幅の基底 W = D_c [A_S; A_Q] (A_S = V_b、A_Q = X̃ = 着衣のモードの基底。張る空間は [I; X] と同じ) の上で
    /// G₀ = W†W、G_a = W† diag(g_a) W、Cm = W†(V − diag|g|²)W (V = 流束座標の非対角の結合: A_S†H'_SS A_S + A_S†H_SQ A_Q + A_Q†(H_QS A_S + E A_Q)、
    /// 疎な積 1 回なら H_QS A_S + E A_Q ≈ R = (H_QS + E X₀)V_b)。SolveDressedTailBlockResearch と同じ式を一般の基底で。
    /// A† diag(w) A (w = D², D² g_x, D² g_y, D² g_z, D² |g|²) は行の塊ごとに積んで足す (5 倍幅の一時行列を全部の行で持たない)</summary>
    private static (Matrix<Complex> Cm, Matrix<Complex> G0, Matrix<Complex> Gx, Matrix<Complex> Gy, Matrix<Complex> Gz) CbedSharedBasis(DressCoreResult cr, BeamSelectionPool pool, bool oneProduct, bool singleTail = false)
    {
        var (sArr, qArr, AS, AQ, HSS, HQS, HSQ) = (cr.SArr, cr.QArr, cr.Vb, cr.Xt, cr.HSS, cr.HQS, cr.HSQ);
        double[] Dc = cr.D;
        int nS = sArr.Length, nQ = qArr.Length, nB = AS.ColumnCount;
        var hssOff = HSS.Clone();
        for (int i = 0; i < nS; i++) hssOff[i, i] = 0;
        var rq = oneProduct ? cr.R : HQS * AS + cr.EXt; // H_QS A_S + E A_Q
        // var Koff = AS.ConjugateTransposeThisAndMultiply(hssOff * AS + HSQ * AQ) + AQ.ConjugateTransposeThisAndMultiply(rq);
        Matrix<Complex> Koff;
        MathNet.Numerics.LinearAlgebra.Complex32.DenseMatrix aq32 = null;
        if (singleTail)
        {
            aq32 = ToC32(AQ);
            Koff = AS.ConjugateTransposeThisAndMultiply(hssOff * AS + ToC64(ToC32(HSQ) * aq32)) + ToC64(aq32.ConjugateTransposeThisAndMultiply(ToC32(rq)));
        }
        else Koff = AS.ConjugateTransposeThisAndMultiply(hssOff * AS + HSQ * AQ) + AQ.ConjugateTransposeThisAndMultiply(rq);
        rq = null;
        var T = new CDense(nB, 5 * nB);
        void Accumulate(Matrix<Complex> A, int[] idx, int r0, int rows, bool single)
        {
            var av = ((CDense)A).Values; int ld = A.RowCount;
            if (single) // 裾の塊: 単精度で積む
            {
                var ac32 = new MathNet.Numerics.Complex32[rows * nB]; var y32 = new MathNet.Numerics.Complex32[rows * 5 * nB];
                for (int col = 0; col < nB; col++)
                    for (int rr = 0; rr < rows; rr++)
                    {
                        int r = r0 + rr, q = idx[r]; double d2 = Dc[q] * Dc[q];
                        double gx = pool.GX[q], gy = pool.GY[q], gz = pool.GZ[q], g2 = Dot3(gx, gx, gy, gy, gz, gz);
                        var x = av[r + col * ld];
                        ac32[rr + col * rows] = new MathNet.Numerics.Complex32((float)x.Real, (float)x.Imaginary);
                        void Put(int blk, double w) { var z = x * w; y32[rr + (col + blk * nB) * rows] = new MathNet.Numerics.Complex32((float)z.Real, (float)z.Imaginary); }
                        Put(0, d2); Put(1, d2 * gx); Put(2, d2 * gy); Put(3, d2 * gz); Put(4, d2 * g2);
                    }
                var p32 = new MathNet.Numerics.LinearAlgebra.Complex32.DenseMatrix(rows, nB, ac32).ConjugateTransposeThisAndMultiply(new MathNet.Numerics.LinearAlgebra.Complex32.DenseMatrix(rows, 5 * nB, y32));
                T.Add(ToC64(p32), T);
                return;
            }
            var Ac = new CDense(rows, nB); var ac = Ac.Values;
            var Y = new CDense(rows, 5 * nB); var yv = Y.Values;
            for (int col = 0; col < nB; col++)
                for (int rr = 0; rr < rows; rr++)
                {
                    int r = r0 + rr, q = idx[r]; double d2 = Dc[q] * Dc[q];
                    double gx = pool.GX[q], gy = pool.GY[q], gz = pool.GZ[q], g2 = Dot3(gx, gx, gy, gy, gz, gz);
                    var x = av[r + col * ld];
                    ac[rr + col * rows] = x;
                    yv[rr + col * rows] = x * d2;
                    yv[rr + (col + nB) * rows] = x * (d2 * gx);
                    yv[rr + (col + 2 * nB) * rows] = x * (d2 * gy);
                    yv[rr + (col + 3 * nB) * rows] = x * (d2 * gz);
                    yv[rr + (col + 4 * nB) * rows] = x * (d2 * g2);
                }
            T.Add(Ac.ConjugateTransposeThisAndMultiply(Y), T);
        }
        Accumulate(AS, sArr, 0, nS, false);
        const int chunk = 2048;
        for (int r0 = 0; r0 < nQ; r0 += chunk) Accumulate(AQ, qArr, r0, Math.Min(chunk, nQ - r0), singleTail);
        Matrix<Complex> Part(int m) => T.SubMatrix(0, nB, m * nB, nB);
        return (Koff - Part(4), Part(0), Part(1), Part(2), Part(3));
    }

    /// <summary>261110Cl 追加: 単精度の複素の密な行列を倍精度へ</summary>
    private static CDense ToC64(Matrix<MathNet.Numerics.Complex32> m)
    {
        var src = ((MathNet.Numerics.LinearAlgebra.Complex32.DenseMatrix)m).Values;
        var dst = new Complex[src.Length];
        for (int i = 0; i < src.Length; i++) dst[i] = new Complex(src[i].Real, src[i].Imaginary);
        return new CDense(m.RowCount, m.ColumnCount, dst);
    }

    /// <summary>261110Cl 追加: 区画の目安 η (区画の端までの裾の応答の相対の変化の見積もり)。着衣の 1 段 X̃_qj ∝ 1/(λ_q − γ_j) で、方向が δk (|δk| ≤ k·r) 動くと
    /// λ_q は Δλ_q ≈ 2|δk||g_q,t|/P_q、γ_j は Δγ_j ≈ 2|δk|⟨|g_t|⟩_j/P_000 (⟨⟩ = モードの明示の成分の重み) だけ動く。η_j = ‖D_Q X̃_j ∘ (Δλ + Δγ_j)/√(|λ − γ_j|² + δ²)‖ / ‖D [A_S; A_Q]_j‖、
    /// η = max_j η_j (区画の中心で励起されないモードも含める: 晶帯軸では対称で励起されないモードが軸から外れると励起される)</summary>
    private static double CbedPatchIndicator(DressCoreResult cr, BeamSelectionPool pool, Vector3DBase n, double dkMax, double regDepthNm)
    {
        var (sArr, qArr, AS, AQ, lam, gb, Dc) = (cr.SArr, cr.QArr, cr.Vb, cr.Xt, cr.LamQ, cr.Gb, cr.D);
        int nS = sArr.Length, nQ = qArr.Length, nB = AS.ColumnCount;
        double d2 = 1 / (Math.PI * regDepthNm); d2 *= d2;
        var (nx, ny, nz) = n.Tuple;
        double Gt(int i) { double gx = pool.GX[i], gy = pool.GY[i], gz = pool.GZ[i], gn = gx * nx + gy * ny + gz * nz; return Math.Sqrt(Math.Max(0, gx * gx + gy * gy + gz * gz - gn * gn)); }
        double p0 = cr.P[0];
        var dl = new double[nQ];
        for (int r = 0; r < nQ; r++) dl[r] = 2 * dkMax * Gt(qArr[r]) / cr.P[qArr[r]];
        var asv = ((CDense)AS).Values; var aqv = ((CDense)AQ).Values;
        double worst = 0;
        for (int j = 0; j < nB; j++)
        {
            double ws = 0, wsg = 0, num = 0, den;
            for (int i = 0; i < nS; i++) { double a2 = (asv[i + j * nS] * Dc[sArr[i]]).Magnitude; a2 *= a2; ws += a2; wsg += a2 * Gt(sArr[i]); }
            double dg = ws > 0 ? 2 * dkMax * (wsg / ws) / p0 : 0;
            den = ws;
            for (int r = 0; r < nQ; r++)
            {
                var x = aqv[r + j * nQ]; double a2 = (x * Dc[qArr[r]]).Magnitude; a2 *= a2;
                den += a2;
                var dd = lam[r] - gb[j];
                num += a2 * (dl[r] + dg) * (dl[r] + dg) / (dd.Real * dd.Real + dd.Imaginary * dd.Imaginary + d2);
            }
            if (den > 0) worst = Math.Max(worst, Math.Sqrt(num / den));
        }
        return worst;
    }
}
