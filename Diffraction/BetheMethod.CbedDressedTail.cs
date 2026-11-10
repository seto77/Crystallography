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
    /// 4 方向ごとの行列の組み立て (CPU 秒)、5 厚さへの伝播 (固有分解など、CPU 秒)、6 出射の振幅と位相 (CPU 秒)、7 全体 (実時間)、
    /// 3 の内訳: 8 Gram の積 (CbedSharedBasis)、9 一般化固有分解と基底 V への変換、10 モードの打ち切り (中心の固有分解・探り・縮約の行列)。
    /// 2 の内訳 (DressCore の段): 11 選択、12 密なブロック、13 X₀ と H_SQ X₀、14 H_b の固有分解、15 着衣の残り (選び分け・T・倍率)、16 疎な積 1、17 疎な積 2)</summary>
    public double[] CbedDressedTailStageSeconds { get; } = new double[18];

    /// <summary>261111Cl 追加: 従来の CBED (N 本の固定のビーム) でも晶帯軸の対称を使う (既定は使わない)。表面の法線・方向の格子・ビームの組・ポテンシャル (並進の位相込み) を
    /// 正確に移す操作 (<see cref="CbedSymmetryOps"/>) があれば、各軌道の代表の方向だけを解き、残りは振幅を位相込みで写す (各方向は独立に解くので結果は丸めの範囲で同じ)。
    /// ビームの組がその操作で閉じていなければ (末尾の同点の組の扱いなど) その操作は使わない</summary>
    public bool CbedLegacyUseSymmetry { get; set; }
    /// <summary>261111Cl 追加: 直前の従来の CBED で使った対称 (使わなければ空)</summary>
    public string CbedLegacySymmetryDetail { get; private set; } = "";

    /// <summary>261110Cl 追加 (研究用): メモリの予算を上書きする [byte] (0 以下 = 物理メモリの空きから決める。込み合った機械での検査用)</summary>
    public long CbedMemoryBudgetOverride { get; set; }

    /// <summary>261110Cl 追加: CBED の新方式の設定。値で比べる</summary>
    public sealed record CbedDressedTailSettings
    {
        /// <summary>設定の版</summary>
        // public const string Version = "CDT1-261110"; // 261111Cl 変更前
        public const string Version = "CDT2-261111"; // 261111Cl 変更: 晶帯軸の対称の利用 (UseSymmetry) を足した
        /// <summary>E1 の基準の本数 (000 を含む。遅い反射の繰り上げで明示の本数は増える)</summary>
        public int N { get; init; } = 300;
        /// <summary>候補の半径 [nm⁻¹] (明示の繰り上げはこの中の全部から)。R40 では Si [111] 100 kV の 40〜60 nm⁻¹ の HOLZ を落として誤差 6.5e−3 (R60 は 2.6e−4、§5.19)</summary>
        public double PoolRadiusNm { get; init; } = 60;
        /// <summary>裾どうしの結合の閾値 (|U| ≥ η max)</summary>
        public double Eta { get; init; } = 0.05;
        /// <summary>遅い反射の繰り上げの閾値 |Q/P| [nm⁻¹]</summary>
        public double Ps { get; init; } = 0.1;
        /// <summary>着衣の正則化 δ = 1/(πT) の T [nm]</summary>
        public double RegDepthNm { get; init; } = 45;
        /// <summary>区画の半径 (中心から角まで) [mrad]。区画の端でも繰り上げる (PromoteOverPatch) と、四隅の局所解との差は Si [111] 100 kV で r 1 / 2 / 3 / 4 mrad: 5.7e−5 / 4.3e−5 / 9.4e−5 / 1.3e−4 (§5.21)。
        /// 261110Cl: 既定 2 mrad (作者の初期案は 1 mrad。繰り上げで 2 mrad でも 1 mrad (繰り上げなし、1.25e−4) より正確になったので変更、作者に諮る)</summary>
        public double PatchRadiusMrad { get; init; } = 2;
        /// <summary>裾の芯の半径 [nm⁻¹] (NaN = 候補の全部を裾に)。近傍の半径と組で使う (EBSD の本番と同じ形)。既定 40 + 近傍 10 で裾 3.6 万 → 1.4 万本、精度はほぼ同じ (芯 30 は 2.7 倍悪い)</summary>
        public double CoreRadiusNm { get; init; } = 40;
        /// <summary>芯の外の裾 = 明示 (と |Q/P| &lt; NeighborPs の候補) からこの距離 [nm⁻¹] の中 (NaN = 芯の外も全部)</summary>
        public double NeighborRadiusNm { get; init; } = 10;
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
        public double ModeReduceTau { get; init; } = 3e-3; // 261110Cl: 既定 3e−3 (伝播 15〜20 倍速く、強度の誤差 4e−5)
        /// <summary>着衣の直後のモードの選び分け (NaN = 使わない)。明示の着衣のモード V_b (H_b の固有ベクトル) のうち、区画の中心と探りの方向 (8 方位で最も遠い方向) の励起を
        /// H_b + diag(Δλ_s(d)) (明示の励起誤差の変化だけ、裾の変化は無視) の伝播で見積もり、最大のこの比以上のモード (近い固有値の組は一緒に) だけを着衣・基底に使う。
        /// Gram の積と着衣の後半が (nM/nS)² に減る。後段の ModeReduceTau (厳密な 2 側の射影) と併用する前提で、それより小さく取る</summary>
        public double EarlyModeTau { get; init; } = 1e-3; // 261110Cl: 既定 1e−3 (結果を変えず着衣のモード 303 → 180)
        /// <summary>区画の端 (8 方位で最も遠い方向、方向が 16 以下なら全部) のどれかで遅い (|Q/P| &lt; ps) 反射も明示へ繰り上げる (区画の端で Ewald 球に近づく HOLZ を裾にしない)</summary>
        public bool PromoteOverPatch { get; init; } = true; // 261110Cl: 既定で繰り上げる (区画 2 mrad の使い回しの誤差 5.5e−4 → 4.3e−5)
        /// <summary>探りの方向で、縮約の解と全次元の解の出射の強度の差の最大がこれを超えたら、その区画は縮約しない (事後の検査)</summary>
        public double ReductionCheckLimit { get; init; } = 3e-4;
        /// <summary>区画の中の交差の検出 (落とすモードと残すモードの対角が区画の中で交差する組を残す)</summary>
        public bool ReductionCrossingCheck { get; init; } = false; // 261110Cl: 既定で無効 (残したモードで縮約の解が発散した区画があった。一般の位置の探りで足りるかを確かめる)
        /// <summary>裾の積を単精度 (complex64) で積む: 区画の基底の行列のうち裾 (A_Q) の積と、着衣の密な積 (H_SQ X₀・H_QS V_b) と疎な積 (E·X)。明示の部分は倍精度。
        /// 裾の寄与は小さいので丸めの影響は小さい (基底の行列の単精度で強度の差 5e−8、day11r)</summary>
        public bool SingleTailProducts { get; init; } = true; // 261110Cl: 既定で単精度 (強度の差 5e−8)
        /// <summary>261111Cl 追加: 晶帯軸の対称を使う。ビーム軸まわりの 2 回・4 回の回転と縦・横・斜めの鏡映のうち、表面の法線・方向の格子・候補のビーム・
        /// ポテンシャル (U と U'、並進の位相込み) を正確に移すものがあれば、代表の方向だけを解き、残りの方向の振幅は位相込みで写す (<see cref="CbedSymmetryOps"/>)</summary>
        public bool UseSymmetry { get; init; } = true;
        /// <summary>261111Cl 追加: 後段の打ち切り (ModeReduceTau) でも、区画のどこか (中心・8 方位の端) で遅い (|Q/P| &lt; ps) か区画の中で Q/P の符号が変わる明示の出射のビームに
        /// 重み (出射のビームの上の振幅の 2 乗の 5 % 以上) のあるモードは残す (早い選び分けと同じ規則)。探り (4 隅 + 中の 4 点) から離れた所だけを通る HOLZ の線のモードを落とさない
        /// (Si [111] 100 kV・区画 2 mrad の中ほどの線で 1.5e−3、day13f)</summary>
        public bool LateSlowKeep { get; init; } = true;
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
        if (cancelled?.Invoke() == true) return (null, null); // 261110Cl (codex2 第 37 回): 取り消し (呼び出し側は CancellationPending で区別する)
        if (thicknesses == null || thicknesses.Length == 0 || directions == null || directions.Count == 0 || !(set.PatchRadiusMrad > 0) || !(set.RegDepthNm > 0))
        { CbedBeamSelectionFallbackReason = "invalid input for the dressed-tail method (thicknesses, directions, patch radius or regularization depth)"; return (null, null); }
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
        // long perPatch = (long)Math.Min(long.MaxValue / 4, 12 * nQEst * nSEst * 16); // 261110Cl 変更前
        // 261110Cl 変更: 着衣と基底で同時に生きる nQ 行の行列 ≈ 3 nS (H_QS・H_SQ・X₀) + 10 nM (T・X₀V・疎な積の作業・R・X̃・E X̃・基底の一時)。早い選び分けなら nM ≈ 0.65 nS
        double nMEst = double.IsFinite(set.EarlyModeTau) && set.EarlyModeTau > 0 ? 0.65 * nSEst : nSEst;
        long perPatch = (long)Math.Min(long.MaxValue / 4, (3 * nSEst + 10 * nMEst) * nQEst * 16);
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
        if (cancelled?.Invoke() == true) return (null, null);
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
        // 261111Cl 追加: 晶帯軸の対称。使える操作があれば、各軌道の代表 (区画の数が最も少なくなる向きと升の格子のずれを選ぶ) だけを解き、最後に残りの方向へ写す
        int nAll = nDir; var allDirections = directions;
        double symOffU = 0, symOffV = 0, symCos = 1, symSin = 0, sideMinSym = Math.Sqrt(2) * set.PatchRadiusMrad * 1e-3; // 升の格子のずれと回転 (斜めの鏡映の境に辺を合わせる。対称を使わなければ 0)
        // var symOps = set.UseSymmetry ? CbedSymmetryOps(pool, directions, n) : []; // 261111Cl 変更前
        var symOps = set.UseSymmetry ? CbedSymmetryOps(pool.Hkl, pool.GX, pool.GY, pool.GZ, pool.IndexOf, pool.Mat, Math.Min(25, 0.999 * pool.RadiusNm), directions, n) : [];
        int[] repOrig = null; (int rep, int op)[] mapFrom = null; (double, double) bestOff = (0, 0); (double, double) bestRot = (1, 0);
        bool adaptiveSym = double.IsFinite(set.AdaptiveMaxRadiusMrad) && double.IsFinite(set.AdaptiveIndicatorLimit) && set.AdaptiveMaxRadiusMrad > set.PatchRadiusMrad; // 適応的な区画 (Split は軸にそろった升) では回さない
        if (symOps.Length > 0)
        {
            int bestCount = int.MaxValue;
            foreach (var (px, py) in new[] { (1.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, -1.0), (1.0, 0.31), (0.31, 1.0) })
            {
                // 軌道の代表 = (p·d, p⊥·d) が辞書順で最大の方向 (同じなら番号の小さい方)
                (double, double) KeyOf(int i) { var d = directions[i]; return (Math.Round((px * d.X + py * d.Y) * 1e12), Math.Round((-py * d.X + px * d.Y) * 1e12)); }
                var reps = new List<int>();
                for (int i = 0; i < nDir; i++)
                {
                    var ki = KeyOf(i); bool isRep = true;
                    foreach (var op in symOps) { int j = op.DirMap[i]; if (j == i) continue; var kj = KeyOf(j); if (kj.Item1 > ki.Item1 || (kj.Item1 == ki.Item1 && (kj.Item2 > ki.Item2 || (kj.Item2 == ki.Item2 && j < i)))) { isRep = false; break; } }
                    if (isRep) reps.Add(i);
                }
                foreach (var rot in adaptiveSym ? [0.0] : new[] { 0.0, Math.PI / 4 })
                    foreach (var (ou, ov) in new[] { (0.0, 0.0), (0.5, 0.0), (0.0, 0.5), (0.5, 0.5) })
                    {
                        (symOffU, symOffV, symCos, symSin) = (ou * sideMinSym, ov * sideMinSym, Math.Cos(rot), Math.Sin(rot));
                        int cnt = Squares(reps, sideMinSym).Count;
                        if (cnt < bestCount) { bestCount = cnt; repOrig = [.. reps]; bestOff = (symOffU, symOffV); bestRot = (symCos, symSin); }
                    }
            }
            (symOffU, symOffV) = bestOff; (symCos, symSin) = bestRot;
            // 代表でない方向 i = T_k (代表 r) の (r の代表の中の番号, k)
            var repPos = new int[nDir]; Array.Fill(repPos, -1);
            for (int r = 0; r < repOrig.Length; r++) repPos[repOrig[r]] = r;
            mapFrom = new (int, int)[nDir]; Array.Fill(mapFrom, (-1, -1));
            for (int r = 0; r < repOrig.Length; r++)
                for (int k = 0; k < symOps.Length; k++) { int i = symOps[k].DirMap[repOrig[r]]; if (repPos[i] < 0 && mapFrom[i].Item1 < 0) mapFrom[i] = (r, k); }
            for (int i = 0; i < nDir; i++) if (repPos[i] < 0 && mapFrom[i].Item1 < 0) throw new InvalidOperationException("CBED symmetry: a direction is not reached from its representative");
            directions = repOrig.Select(i => allDirections[i]).ToArray();
            uu = repOrig.Select(i => uu[i]).ToArray(); vv = repOrig.Select(i => vv[i]).ToArray();
            nDir = repOrig.Length;
        }
        // 升 (中心 (uc, vc)・一辺 side・方向の番号)。原点を升の中心にした格子で分ける
        List<(double uc, double vc, double side, int[] members)> Squares(IEnumerable<int> idx, double sd)
        {
            var g = new Dictionary<(long, long), List<int>>();
            foreach (var i in idx)
            {
                // var key = ((long)Math.Floor(uu[i] / sd + 0.5), (long)Math.Floor(vv[i] / sd + 0.5)); // 261111Cl 変更前
                double pu = symCos * uu[i] + symSin * vv[i], pv = -symSin * uu[i] + symCos * vv[i]; // 261111Cl: 回した格子の座標 (対称を使わなければ uu, vv のまま)
                var key = ((long)Math.Floor((pu - symOffU) / sd + 0.5), (long)Math.Floor((pv - symOffV) / sd + 0.5)); // 261111Cl 変更: 対称の代表の境を升の辺に合わせるずれと回転
                if (!g.TryGetValue(key, out var list)) g[key] = list = [];
                list.Add(i);
            }
            // return g.OrderBy(kv => kv.Key).Select(kv => ((double)kv.Key.Item1 * sd, (double)kv.Key.Item2 * sd, sd, kv.Value.ToArray())).ToList(); // 261111Cl 変更前
            return g.OrderBy(kv => kv.Key).Select(kv => { double pc = kv.Key.Item1 * sd + symOffU, qc = kv.Key.Item2 * sd + symOffV; return (symCos * pc - symSin * qc, symSin * pc + symCos * qc, sd, kv.Value.ToArray()); }).ToList(); // 261111Cl 変更 (中心は uu, vv の座標へ戻す)
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
        // 升の中心 (uc, vc) から 8 方位で最も遠い方向 (重複は 1 つに)
        int[] Extremes(int[] members, double uc, double vc)
        {
            var ex = new HashSet<int>();
            for (int m8 = 0; m8 < 8; m8++)
            {
                double ca = Math.Cos(m8 * Math.PI / 4), sa = Math.Sin(m8 * Math.PI / 4), best = double.NegativeInfinity; int bi = members[0];
                foreach (var i in members) { double pr = (uu[i] - uc) * ca + (vv[i] - vc) * sa; if (pr > best) { best = pr; bi = i; } }
                ex.Add(bi);
            }
            return [.. ex];
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
        var selNQ = new int[patches.Length];
        Parallel.For(0, patches.Length, new ParallelOptions { MaxDegreeOfParallelism = threads }, p =>
        {
            var ck = getVecK0(kvac, u0, patches[p].center, n);
            // 261110Cl: 区画の端でも繰り上げるなら、出射のビームの選択でも升の端 (8 方位で最も遠い方向) で繰り上げる (端で繰り上げた HOLZ の円盤を落とさない)
            (Vector3DBase, Vector3DBase)[] promoteSel = null;
            if (set.PromoteOverPatch && patches[p].members.Length > 1)
            {
                var mem = patches[p].members; double mu0 = mem.Average(i => uu[i]), mv0 = mem.Average(i => vv[i]);
                promoteSel = Extremes(mem, mu0, mv0).Select(di => (getVecK0(kvac, u0, directions[di], n), n)).ToArray();
            }
            var crSel = DressCore(table, ck, n, set.N, set.Ps, set.RegDepthNm, set.OneProduct, false, null, promoteAlso: promoteSel, selection: countSel);
            selections[p] = crSel.SArr; selNQ[p] = crSel.QArr.Length;
        });
        if (cancelled?.Invoke() == true) return (null, null);
        foreach (var p in patchOrder)
            foreach (var i in selections[p])
                if (!inOut[i]) { inOut[i] = true; outPool.Add(i); }
        // 261111Cl 追加: 対称で写す先のビームも出射に含める (出射のビームの組を対称の操作で閉じる)
        for (bool grown = symOps.Length > 0; grown;)
        {
            grown = false;
            foreach (var op in symOps)
                for (int o = 0; o < outPool.Count; o++) { int j = op.BeamMap[outPool[o]]; if (!inOut[j]) { inOut[j] = true; outPool.Add(j); grown = true; } }
        }
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

        // 261110Cl (codex2 第 37 回): 選択の段の実際の明示・裾の数でメモリを見積もり直す。表を作った後の空きから、出力の振幅 (方向 × 出射のビーム × 厚さ) を引く
        {
            int maxNS = selections.Max(a => a.Length), maxNQ = selNQ.Max();
            double nMReal = double.IsFinite(set.EarlyModeTau) && set.EarlyModeTau > 0 ? 0.75 * maxNS : maxNS;
            perPatch = (long)Math.Min(long.MaxValue / 4, (3.0 * maxNS + 10.0 * nMReal) * maxNQ * 16);
            // long outMem = (long)nDir * nOut * tLen * 16; // 261111Cl 変更前
            long outMem = (long)nAll * nOut * tLen * 16; // 261111Cl 変更: 対称で写した後の全方向
            var (total2, free2) = PhysicalMemory();
            long budget2 = (CbedMemoryBudgetOverride > 0 ? CbedMemoryBudgetOverride : Math.Min(total2 / 2, free2 - reserve)) - outMem;
            byBudget = (int)Math.Min(int.MaxValue, Math.Max(0, budget2 / Math.Max(1, perPatch)));
            if (byBudget < 1)
            {
                CbedBeamSelectionFallbackReason = $"memory {perPatch / 1e9:f1} GB per patch (explicit {maxNS}, tail {maxNQ}) and {outMem / 1e9:f1} GB of output exceed the budget {(budget2 + outMem) / 1e9:f1} GB";
                return (null, null);
            }
        }

        // 並列: 区画を外側で (メモリの予算の内)、区画の中の方向を内側で
        int outer = Math.Max(1, Math.Min(Math.Min(threads, patches.Length), byBudget));
        int inner = 1; // 区画の中で使ったスレッドの数の最大 (記録用)
        var amplitudes = new Complex[nDir][];
        var iu0 = Complex.ImaginaryOne * getU(AccVoltage, new Beam((0, 0, 0), new Vector3DBase(0, 0, 0)), cache: false).Imag;
        int done = 0, nPatchesDone = 0;
        long sumS = 0, sumQ = 0, sumB = 0, sumNb = 0; int minS = int.MaxValue, maxS = 0, nBackward = 0; double etaMin = double.PositiveInfinity, etaMax = 0;
        var patchEta = new List<(Vector3DBase center, double eta)>();
        int nSplit = 0, nSolvedPatches = 0;
        var gate = new object();
        var backward = new System.Collections.Concurrent.ConcurrentBag<int>();
        var stage = CbedDressedTailStageSeconds;
        // 出射面の位相の厚さの因子 (方向に依らない部分): exp(πi(P_g − 2k Surface.Z)t) を方向ごとに作る
        var twoPiI = 2 * Math.PI * Complex.ImaginaryOne;

        // 方向 di で、候補の番号 idx のどれかが後ろ向き (P ≤ 0) か
        bool BackwardAt(int di, int[] idx)
        {
            var kd = getVecK0(kvac, u0, directions[di], n); var (kx, ky, kz) = kd.Tuple;
            foreach (var i in idx) if (2 * Dot3(nx, pool.GX[i] + kx, ny, pool.GY[i] + ky, nz, pool.GZ[i] + kz) <= 0) return true;
            return false;
        }
        int nReduceSkipped = 0, nCrossKept = 0, nGuardResolved = 0; double maxProbeDiff = 0;
        long nSlowKept = 0; // 261111Cl 追加: 後段の打ち切りで遅いビームのために残したモードの数 (区画の和)

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
        // 261110Cl (codex2 第 37 回): innerLocal = この区画の中で使うスレッドの数 (外側で同時に動く区画の数で割った分)。DressCore の上限なしの並列は区画が 1 つだけのときに限る
        void SolvePatch((double uc, double vc, double side, int[] members) sq, int innerLocal, System.Collections.Concurrent.ConcurrentBag<(double, double, double, int[])> next)
        {
            if (cancelled?.Invoke() == true) return;
            bool parallelInside = innerLocal >= threads && threads > 1;
            var members = sq.members;
            var center = CenterOf(members);
            double rEff = 0; foreach (var i in members) rEff = Math.Max(rEff, Math.Acos(Math.Min(1, directions[i] * center))); // 中心から最も遠い方向までの角度 [rad]
            var swp = Stopwatch.StartNew();
            var ck = getVecK0(kvac, u0, center, n);
            // 探り: 区画の中心から 8 方位で最も遠い方向
            var probes = new HashSet<int>(); var probesDiag = new HashSet<int>(); // 8 方位 (後段の打ち切り) と斜めの 4 つ = 四隅 (着衣の直後の選び分け)
            if (members.Length > 9)
                for (int m8 = 0; m8 < 8; m8++)
                {
                    double ca = Math.Cos(m8 * Math.PI / 4), sa = Math.Sin(m8 * Math.PI / 4), best = double.NegativeInfinity; int bi = members[0];
                    foreach (var i in members) { double pr = (uu[i] - sq.uc) * ca + (vv[i] - sq.vc) * sa; if (pr > best) { best = pr; bi = i; } }
                    probes.Add(bi);
                    if (m8 % 2 == 1) probesDiag.Add(bi);
                }
            // 261110Cl 追加: 励起の見積もりの探り (早い選び分け・後段の打ち切り) は対称の線から外れた一般の位置に置く: 方位を 0.37 rad ずらした 8 方向の端と、区画の中の非対称な 4 点。
            //   晶帯軸の区画では x・y 軸と対角線が鏡映の線に乗り、そこでは鏡映で反対称なモードの励起が 0 に見える (Si [110] 200 kV・区画 2 mrad で 4.4e−3 を見落とした)。
            //   繰り上げ (λ は方向の一次なので四隅で足りる) には今までの 8 方位 (probes) を使う
            var eprobes = new HashSet<int>(); var eprobesEarly = new HashSet<int>(); // 後段 (12 点) と、早い選び分け (6 点: ずらした方位の斜めの 4 つ + 中の 2 点。近似なので少なく)
            if (members.Length > 9)
            {
                for (int m8 = 0; m8 < 8; m8++)
                {
                    double ang = m8 * Math.PI / 4 + 0.37, ca = Math.Cos(ang), sa = Math.Sin(ang), best = double.NegativeInfinity; int bi = members[0];
                    foreach (var i in members) { double pr = (uu[i] - sq.uc) * ca + (vv[i] - sq.vc) * sa; if (pr > best) { best = pr; bi = i; } }
                    eprobes.Add(bi); if (m8 % 2 == 1) eprobesEarly.Add(bi);
                }
                int qn = 0;
                foreach (var (qu, qv) in new[] { (0.23, 0.11), (-0.13, 0.27), (-0.27, -0.09), (0.07, -0.29) })
                {
                    double tu = sq.uc + qu * sq.side, tv = sq.vc + qv * sq.side; int bi = -1; double bd = double.PositiveInfinity;
                    foreach (var i in members) { double d2u = (uu[i] - tu) * (uu[i] - tu) + (vv[i] - tv) * (vv[i] - tv); if (d2u < bd) { bd = d2u; bi = i; } }
                    if (bi >= 0) { eprobes.Add(bi); if (qn++ % 2 == 0) eprobesEarly.Add(bi); }
                }
            }
            // 着衣の直後のモードの選び分け (EarlyModeTau): 探りの方向 (四隅) の励起を H_b + diag(Δλ_s) で見積もる
            Func<Matrix<Complex>, MathNet.Numerics.LinearAlgebra.Vector<Complex>, Matrix<Complex>, int[], bool[]> selector = null;
            if (double.IsFinite(set.EarlyModeTau) && set.EarlyModeTau > 0 && probes.Count > 0)
                selector = (Vb, gbv, hb, sArrSel) =>
                {
                    int ns = sArrSel.Length;
                    var vbInv = Vb.Inverse();
                    if (InverseResidual(vbInv, Vb) > 1e-6) return Enumerable.Repeat(true, ns).ToArray(); // 261110Cl (codex2 第 37 回): V_b が悪条件なら選び分けない
                    var vnorm = new double[ns]; for (int j = 0; j < ns; j++) vnorm[j] = Vb.Column(j).L2Norm();
                    var w = new double[ns];
                    var e0 = MathNet.Numerics.LinearAlgebra.Vector<Complex>.Build.Dense(ns, i => i == 0 ? Complex.One : Complex.Zero);
                    var a0 = vbInv * e0; // 中心
                    for (int j = 0; j < ns; j++) w[j] = a0[j].Magnitude * vnorm[j];
                    var (cx, cy, cz) = ck.Tuple;
                    double Lam(int s0, double kx, double ky, double kz) { double vx = pool.GX[s0] + kx, vy = pool.GY[s0] + ky, vz = pool.GZ[s0] + kz; double k2 = kx * kx + ky * ky + kz * kz; return (k2 - (vx * vx + vy * vy + vz * vz)) / (2 * (nx * vx + ny * vy + nz * vz)); }
                    // 探りの方向 (四隅 = 8 方位の斜めの 4 つ) の励起を H_d = H_b + diag(Δλ_s(d)) の固有分解で厳密に伝播し、中心のモードへ射影する
                    //   (一次の摂動の上界は 2 段の経由の励起を落とし (誤差 3.9e−4)、くり返すと保守的すぎてモードが減らない (261110Cl の試し) ので使わない)
                    // var probeList = probesDiag.ToArray(); // 261110Cl 変更前
                    var probeList = eprobesEarly.Where(di => !BackwardAt(di, sArrSel)).ToArray(); // 261110Cl (codex2 第 37 回): 明示に後ろ向きのビームがある探りは使わない。一般の位置の探り (eprobes) を使う
                    var zs = new Matrix<Complex>[probeList.Length]; // 探りごとの中心のモードの係数 (ns × tLen)
                    Parallel.For(0, probeList.Length, new ParallelOptions { MaxDegreeOfParallelism = Math.Max(1, innerLocal) }, pi =>
                    {
                        var kd = getVecK0(kvac, u0, directions[probeList[pi]], n); var (dx, dy, dz) = kd.Tuple;
                        var hd = hb.Clone();
                        for (int i = 0; i < ns; i++) hd[i, i] += Lam(sArrSel[i], dx, dy, dz) - Lam(sArrSel[i], cx, cy, cz);
                        var evd = hd.Evd(MathNet.Numerics.LinearAlgebra.Symmetricity.Asymmetric);
                        var ad = evd.EigenVectors.Solve(e0);
                        var tr = vbInv * evd.EigenVectors;
                        var ph = new CDense(ns, tLen);
                        for (int t = 0; t < tLen; t++) for (int j = 0; j < ns; j++) ph[j, t] = Complex.Exp(twoPiI * evd.EigenValues[j] * thicknesses[t]) * ad[j];
                        zs[pi] = tr * ph;
                    });
                    // 中心: z(t) = e^{2πiγt} a0
                    var zc = new CDense(ns, tLen);
                    for (int t = 0; t < tLen; t++) for (int j = 0; j < ns; j++) zc[j, t] = Complex.Exp(twoPiI * gbv[j] * thicknesses[t]) * a0[j];
                    var keepE = KeepClusters(gbv, Vb, zs.Append(zc), set.EarlyModeTau, perMode: true);
                    // 261110Cl: 区画のどこか (中心か探り) で遅い (|Q/P| < ps) 明示の反射 (端で繰り上げた HOLZ など) に重み (5 % 以上) のあるモードは必ず残す
                    //   (H_b + diag(Δλ_s) の見積もりはこれらの励起を約 3 倍小さく見て、繰り上げと併用すると 1.5e−4 に悪化した)
                    var slow = new bool[ns];
                    for (int i = 0; i < ns; i++)
                    {
                        double l0 = Lam(sArrSel[i], cx, cy, cz), lo = l0, hi = l0;
                        if (Math.Abs(l0) < set.Ps) { slow[i] = true; continue; }
                        foreach (var di in probes) { var kd = getVecK0(kvac, u0, directions[di], n); double ld = Lam(sArrSel[i], kd.X, kd.Y, kd.Z); lo = Math.Min(lo, ld); hi = Math.Max(hi, ld); if (Math.Abs(ld) < set.Ps) { slow[i] = true; break; } }
                        if (lo < 0 && hi > 0) slow[i] = true; // 261110Cl: 区画の中で符号が変わる (Bragg の線が横切る)
                    }
                    for (int j = 0; j < ns; j++)
                    {
                        if (keepE[j]) continue;
                        double tot = 0, sl = 0;
                        for (int i = 0; i < ns; i++) { double a2 = Vb[i, j].Magnitude; a2 *= a2; tot += a2; if (slow[i]) sl += a2; }
                        if (tot > 0 && sl >= 0.05 * tot) keepE[j] = true;
                    }
                    return keepE; // 261110Cl: 早い選び分けはモードごとの最大 (組の寄与 ‖R_G z_G‖ では 1.5e−4 に悪化した。近い組の打ち消しで小さく見積もる)
                };
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
            // 縮約の方向の解: A_red = Σ c_M (P_M − Δnk Q_M)、y0_red = L_k y0(d) (redP・redQ・Lk・nR・nkc を作った後で使う)
            int nR = 0; CDense[] redP = null, redQ = null; Matrix<Complex> Lk = null; CDense wpR = null; double nkc = Dot3(nx, ck.X, ny, ck.Y, nz, ck.Z);
            Complex[] ReducedY(int di)
            {
                var k0 = getVecK0(kvac, u0, directions[di], n);
                var (kx, ky, kz) = k0.Tuple;
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
                return prop == 1 ? NativeWrapper.CBEDSolver_MatExp(ar, y0r, thicknesses) : prop == 0 ? NativeWrapper.CBEDSolver_Eigen(ar, y0r, thicknesses) : ManagedPropagate(ar, y0r, nR);
            }
            // 着衣はモードの基底 (A_S = V_b、A_Q = X̃、張る空間は [I; X] と同じ) のまま使う: X = X̃V_b⁻¹ の積と逆行列を省き、τ > 0 ならモードを打ち切る
            var dtm = new DressedTailTiming();
            (Vector3DBase k0, Vector3DBase n)[] promote = null;
            if (set.PromoteOverPatch && members.Length > 1)
                promote = (members.Length <= 16 ? members : probes.ToArray()).Select(di => (getVecK0(kvac, u0, directions[di], n), n)).ToArray();
            var cr = DressCore(table, ck, n, set.N, set.Ps, set.RegDepthNm, set.OneProduct, parallelInside, dtm, promoteAlso: promote, modeBasis: true, selection: sel, modeTau: set.ModeTau, modeSelector: selector,
                denseSingle: set.SingleTailProducts, spmmSingle: set.SingleTailProducts);
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
            double tb0 = swp.Elapsed.TotalSeconds;
            var (Cm, G0, Gx, Gy, Gz) = CbedSharedBasis(cr, pool, set.OneProduct, set.SingleTailProducts);
            double tb1 = swp.Elapsed.TotalSeconds;
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
            double tb2 = swp.Elapsed.TotalSeconds;
            nR = nB;
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
                // 探り (前方の条件を満たすもの) を全次元で解く
                var zsL = new List<Matrix<Complex>>();
                foreach (var di in eprobes.Where(di => !BackwardAt(di, all)))
                {
                    var (a, y0, _, _) = FullSystem(di);
                    var yt = prop == 0 ? NativeWrapper.CBEDSolver_Eigen(a, y0, thicknesses) : prop == 1 ? NativeWrapper.CBEDSolver_MatExp(a, y0, thicknesses) : ManagedPropagate(a, y0, nB); // 261110Cl (codex2 第 37 回): 通常の方向と同じ 3 分岐
                    probeDone[di] = yt;
                    zsL.Add(Lc * new CDense(nB, tLen, yt));
                }
                // 区画の中心 (方向の平均) の励起: z(t) = e^{2πiγt} L_c y0_c
                {
                    var y0c = new Complex[nB]; for (int i = 0; i < nB; i++) y0c[i] = 2 * nkc * scc[i] * v0[i];
                    var zc0 = Lc * MathNet.Numerics.LinearAlgebra.Vector<Complex>.Build.Dense(y0c);
                    var zc = new CDense(nB, tLen);
                    for (int t = 0; t < tLen; t++) for (int j = 0; j < nB; j++) zc[j, t] = Complex.Exp(twoPiI * gc[j] * thicknesses[t]) * zc0[j];
                    zsL.Add(zc);
                }
                // 261110Cl (codex2 第 37 回): 近い固有値の組は全組の連結成分、組の寄与 ‖R_G z_G‖ で判定 (実部の隣どうしでは複素で近い組を落とす)
                var keepJ = KeepClusters(gc, Rc, zsL, set.ModeReduceTau);
                // 261111Cl 追加: 区画のどこかで遅い明示の出射のビームに重みのあるモードは残す (LateSlowKeep)
                if (set.LateSlowKeep && probes.Count > 0)
                {
                    var (sx, sy, sz) = ck.Tuple;
                    double LamS(int s0, double kx, double ky, double kz) { double vx = pool.GX[s0] + kx, vy = pool.GY[s0] + ky, vz = pool.GZ[s0] + kz; double k2 = kx * kx + ky * ky + kz * kz; return (k2 - (vx * vx + vy * vy + vz * vz)) / (2 * (nx * vx + ny * vy + nz * vz)); }
                    var probeK = probes.Select(di => getVecK0(kvac, u0, directions[di], n)).ToArray();
                    var slowO = new List<int>();
                    for (int o = 1; o < nOut; o++)
                    {
                        int pi = outPool[o];
                        if (!posS.ContainsKey(pi)) continue;
                        double l0 = LamS(pi, sx, sy, sz), lo = l0, hi = l0; bool sl = Math.Abs(l0) < set.Ps;
                        foreach (var kd in probeK) { double ld = LamS(pi, kd.X, kd.Y, kd.Z); lo = Math.Min(lo, ld); hi = Math.Max(hi, ld); if (Math.Abs(ld) < set.Ps) sl = true; }
                        if (sl || (lo < 0 && hi > 0)) slowO.Add(o);
                    }
                    if (slowO.Count > 0)
                    {
                        var Mo = wp * Rc; // 出射のビームの上のモードの物理振幅 (nOut × nB)
                        int added = 0;
                        for (int j = 0; j < nB; j++)
                        {
                            if (keepJ[j]) continue;
                            double tot = 0, sw2 = 0;
                            for (int o = 0; o < nOut; o++) { double a2 = Mo[o, j].Magnitude; tot += a2 * a2; }
                            foreach (var o in slowO) { double a2 = Mo[o, j].Magnitude; sw2 += a2 * a2; }
                            if (tot > 0 && sw2 >= 0.05 * tot) { keepJ[j] = true; added++; }
                        }
                        lock (gate) nSlowKept += added;
                    }
                }
                if (set.ReductionCrossingCheck) // 261110Cl 追加: 区画の中の交差の検出。探りは区画の端だけなので、中で残すモード j と落とすモード k が交差 (縮退) して振幅が移る所を見落とす
                //   (Si [110] 200 kV・区画 2 mrad で FOLZ のモードが交差し、端の探りの照合 6.9e−5 のまま中で 4.4e−3・区画の境界の段差)。
                //   中心の固有モードの基底で A'(d) ≈ A'_c + Σ δk_a E_a (E_a = L_c S_c (−2B_a) R_c) とし、区画の全方向で対角の差 Δ_kj(d) = (γ_k + E_kk(d)) − (γ_j + E_jj(d)) を見る。
                //   移る振幅 ≈ max|E_kj| |z_j| min(2/min|Δ|, 2π t_max) (Δ の実部が符号を変えれば交差 = min|Δ| 0) が τ 以上なら k も残す。経由を拾うため 2 回くり返す
                {
                    var Ea = new Matrix<Complex>[3]; var bsrc = new[] { bx, by, bz };
                    for (int a = 0; a < 3; a++)
                    {
                        var Ma = new CDense(nB, nB); var mv = Ma.Values; var src = bsrc[a];
                        for (int col = 0; col < nB; col++) { int o = col * nB; for (int i = 0; i < nB; i++) mv[o + i] = -2 * scc[i] * src[o + i]; }
                        Ea[a] = Lc * (Ma * Rc);
                    }
                    var (ckx2, cky2, ckz2) = ck.Tuple;
                    var dks = members.Select(di => { var kd = getVecK0(kvac, u0, directions[di], n); return (kd.X - ckx2, kd.Y - cky2, kd.Z - ckz2); }).ToArray();
                    double dxm = dks.Max(d => Math.Abs(d.Item1)), dym = dks.Max(d => Math.Abs(d.Item2)), dzm = dks.Max(d => Math.Abs(d.Item3));
                    var rnormC = new double[nB]; for (int j = 0; j < nB; j++) rnormC[j] = Rc.Column(j).L2Norm();
                    var zmag = new double[nB];
                    foreach (var z in zsL) for (int j = 0; j < nB; j++) for (int t = 0; t < tLen; t++) zmag[j] = Math.Max(zmag[j], z[j, t].Magnitude);
                    double wAll = 0; for (int j = 0; j < nB; j++) wAll = Math.Max(wAll, zmag[j] * rnormC[j]);
                    double tmax = thicknesses.Max(), cap = 2 * Math.PI * tmax;
                    var dg = new Complex[3][]; for (int a = 0; a < 3; a++) { dg[a] = new Complex[nB]; for (int j = 0; j < nB; j++) dg[a][j] = Ea[a][j, j]; }
                    for (int iter = 0; iter < 2; iter++)
                    {
                        var addW = new double[nB];
                        for (int k = 0; k < nB; k++)
                        {
                            if (keepJ[k]) continue;
                            double best = 0;
                            for (int j = 0; j < nB; j++)
                            {
                                if (!keepJ[j] || zmag[j] == 0) continue;
                                double ekj = dxm * Ea[0][k, j].Magnitude + dym * Ea[1][k, j].Magnitude + dzm * Ea[2][k, j].Magnitude;
                                if (ekj * zmag[j] * cap < set.ModeReduceTau * wAll / Math.Max(rnormC[k], 1e-300)) continue; // 交差しても届かない
                                var d0 = gc[k] - gc[j]; double gmin = double.PositiveInfinity; bool neg = false, pos = false;
                                foreach (var (ex, ey, ez) in dks)
                                {
                                    var dd = d0 + ex * (dg[0][k] - dg[0][j]) + ey * (dg[1][k] - dg[1][j]) + ez * (dg[2][k] - dg[2][j]);
                                    gmin = Math.Min(gmin, dd.Magnitude); if (dd.Real < 0) neg = true; else pos = true;
                                }
                                if (neg && pos) gmin = 0;
                                best = Math.Max(best, ekj * zmag[j] * Math.Min(gmin > 0 ? 2 / gmin : cap, cap));
                            }
                            addW[k] = best;
                        }
                        bool changed = false;
                        for (int k = 0; k < nB; k++)
                            if (!keepJ[k] && addW[k] * rnormC[k] >= set.ModeReduceTau * wAll) { keepJ[k] = true; zmag[k] = Math.Max(zmag[k], addW[k]); changed = true; lock (gate) nCrossKept++; }
                        if (!changed) break;
                    }
                }
                // 健全性: L_c R_c ≈ I、S(d) の一次近似の誤差 max x² (x = Δnk/(n·k_c + μ))
                double xmax = 0; double wmaxMu = 0; for (int i = 0; i < nB; i++) wmaxMu = Math.Max(wmaxMu, 1 / (nkc + mu[i]));
                foreach (var di in members) { var kd = getVecK0(kvac, u0, directions[di], n); xmax = Math.Max(xmax, Math.Abs(Dot3(nx, kd.X, ny, kd.Y, nz, kd.Z) - nkc) * wmaxMu); }
                bool healthy = InverseResidual(Lc, Rc) <= 1e-6 && xmax * xmax <= 1e-9;
                if (!healthy) { keepJ = Enumerable.Repeat(true, nB).ToArray(); lock (gate) nReduceSkipped++; }
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
                // 261110Cl (codex2 第 37 回): 事後の検査。探りの方向で縮約の解と全次元の解の出射の強度を比べ、上限を超えたら縮約しない
                //   261110Cl 追加: 区画の内側の 4 方向 (升の中心から ±side/4) も全次元で解いて比べる (端の探りだけでは中の交差を見落とした)
                foreach (var (qu, qv) in new[] { (-0.31, 0.17), (0.19, -0.13), (0.33, 0.29), (-0.09, -0.34) }) // 261110Cl: 対称の線から外れた点 (探りとも別)
                {
                    double tu = sq.uc + qu * sq.side, tv = sq.vc + qv * sq.side; int bi = -1; double bd = double.PositiveInfinity;
                    foreach (var i in members) { double d2u = (uu[i] - tu) * (uu[i] - tu) + (vv[i] - tv) * (vv[i] - tv); if (d2u < bd) { bd = d2u; bi = i; } }
                    if (bi < 0 || probeDone.ContainsKey(bi) || BackwardAt(bi, all)) continue;
                    var (a, y0, _, _) = FullSystem(bi);
                    probeDone[bi] = prop == 0 ? NativeWrapper.CBEDSolver_Eigen(a, y0, thicknesses) : prop == 1 ? NativeWrapper.CBEDSolver_MatExp(a, y0, thicknesses) : ManagedPropagate(a, y0, nB);
                }
                double probeDiff = 0;
                foreach (var kv in probeDone)
                {
                    var yr = ReducedY(kv.Key);
                    var pr = (wpR * new CDense(nR, tLen, yr)).Values; var pf = (wp * new CDense(nB, tLen, kv.Value)).Values;
                    for (int i = 0; i < pr.Length; i++) probeDiff = Math.Max(probeDiff, Math.Abs(pr[i].Magnitude * pr[i].Magnitude - pf[i].Magnitude * pf[i].Magnitude));
                }
                lock (gate) maxProbeDiff = Math.Max(maxProbeDiff, probeDiff);
                if (probeDiff > set.ReductionCheckLimit || nR == nB) { redP = null; nR = nB; if (nR != nB || probeDiff > set.ReductionCheckLimit) lock (gate) nReduceSkipped++; }
            }
            double tBasis = swp.Elapsed.TotalSeconds;
            if (cancelled?.Invoke() == true) return;
            double tForm = 0, tProp = 0, tOut = 0;
            Parallel.For(0, members.Length, new ParallelOptions { MaxDegreeOfParallelism = Math.Max(1, innerLocal) }, () => (0.0, 0.0, 0.0), (m, _, acc) =>
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
                    t1 = swd.Elapsed.TotalSeconds;
                    yt = ReducedY(di);
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
                // 261110Cl 追加: 縮約の解の安全装置。斜めの射影は非正規な行列で成長するにせのモードを作りうる → 出射の強度の和が 1.1 を超えたら (吸収のある結晶では起きない) 全次元で解き直す
                if (nUse != nB)
                {
                    bool bad = false;
                    for (int t = 0; t < tLen && !bad; t++) { double sI = 0; for (int o = 0; o < nOut; o++) { var z = psi[o + t * nOut]; sI += z.Real * z.Real + z.Imaginary * z.Imaginary; } if (!(sI <= 1.1)) bad = true; }
                    if (bad)
                    {
                        var (a, y0, _, _) = FullSystem(di);
                        yt = prop == 0 ? NativeWrapper.CBEDSolver_Eigen(a, y0, thicknesses) : prop == 1 ? NativeWrapper.CBEDSolver_MatExp(a, y0, thicknesses) : ManagedPropagate(a, y0, nB);
                        psi = (wp * new CDense(nB, tLen, yt)).ToColumnMajorArray();
                        Interlocked.Increment(ref nGuardResolved);
                    }
                }
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
                // progressed?.Invoke(dn); // 261111Cl 変更前
                progressed?.Invoke(nAll == nDir ? dn : (int)((long)dn * nAll / nDir)); // 261111Cl 変更: 対称の代表だけを解くときは全方向に換算
                double t3 = swd.Elapsed.TotalSeconds;
                return (acc.Item1 + t1, acc.Item2 + t2 - t1, acc.Item3 + t3 - t2);
            }, acc => { lock (gate) { tForm += acc.Item1; tProp += acc.Item2; tOut += acc.Item3; } });
            lock (gate)
            {
                stage[2] += tDress; stage[3] += tBasis; stage[8] += tb1 - tb0; stage[9] += tb2 - tb1; stage[10] += tBasis - tb2;
                stage[11] += dtm.Select; stage[12] += dtm.Blocks; stage[13] += dtm.Project; stage[14] += dtm.Eig; stage[15] += dtm.Dress - dtm.Project - dtm.Eig; stage[16] += dtm.SpMM1; stage[17] += dtm.SpMM2; stage[4] += tForm; stage[5] += tProp; stage[6] += tOut; nPatchesDone++; nSolvedPatches++;
                sumS += nS; sumQ += nQ; sumB += nR; sumNb += nB; minS = Math.Min(minS, nS); maxS = Math.Max(maxS, nS); etaMin = Math.Min(etaMin, eta); etaMax = Math.Max(etaMax, eta); patchEta.Add((center, eta));
            }
        }

        var work = adaptive ? Squares(Enumerable.Range(0, nDir), side0) : Squares(Enumerable.Range(0, nDir), sideMin);
        int nLevels = 0;
        while (work.Count > 0 && cancelled?.Invoke() != true)
        {
            var next = new System.Collections.Concurrent.ConcurrentBag<(double, double, double, int[])>();
            // int o2 = Math.Max(1, Math.Min(outer, work.Count)), i2 = Math.Max(1, threads / o2); // 261110Cl 変更前 (区画の中の並列を threads / 外側 に固定 → 最後の巡で区画が少ないとコアが余った)
            // 261110Cl 変更: 区画の中の方向の並列はスレッド数まで許す (入れ子の並列はスレッドプールが分け合う)。メモリは同時に着衣する区画の数 (外側) で抑える。区画は大きい順に 1 つずつ取る
            // 261110Cl (codex2 第 37 回): 区画の中のスレッドは「スレッド数 / 同時に動く区画の数 (残りの区画の数まで)」(入れ子で指定のスレッド数を大きく超えない。始めた時点の値なので近似)
            int o2 = Math.Max(1, Math.Min(outer, work.Count));
            var cur = work.OrderByDescending(q => q.members.Length).ToList();
            int remaining = cur.Count;
            Parallel.ForEach(System.Collections.Concurrent.Partitioner.Create(cur, System.Collections.Concurrent.EnumerablePartitionerOptions.NoBuffering),
                new ParallelOptions { MaxDegreeOfParallelism = o2 }, q =>
                {
                    int rem = Interlocked.Decrement(ref remaining) + 1;
                    int il = Math.Max(1, threads / Math.Max(1, Math.Min(o2, rem)));
                    inner = Math.Max(inner, il);
                    SolvePatch(q, il, next);
                });
            work = next.OrderBy(q => q.Item1).ThenBy(q => q.Item2).ToList();
            nLevels++;
        }
        // 区画の中心で後ろ向きのビームがあった方向は、自分を中心に解き直す (CBED の通常の幾何ではまず起きない)
        var back = backward.Distinct().ToArray();
        nBackward = back.Length;
        // Parallel.For(0, back.Length, new ParallelOptions { MaxDegreeOfParallelism = threads }, b => SolvePatch((uu[back[b]], vv[back[b]], 0, new[] { back[b] }), false, null)); // 261110Cl 変更前 (予算を迂回していた)
        Parallel.For(0, back.Length, new ParallelOptions { MaxDegreeOfParallelism = Math.Max(1, Math.Min(threads, byBudget)) }, b => SolvePatch((uu[back[b]], vv[back[b]], 0, new[] { back[b] }), 1, null));
        stage[7] = swAll.Elapsed.TotalSeconds;

        // 261111Cl 追加: 対称で残りの方向へ写す。ψ_{Tg}(Td) = ψ_g(d) exp(+2πi h_{Tg}·w) (U_{Tg} = U_g exp(−2πi h_{Tg}·w) の共役。出射面の位相 exp(πi(P_g − 2k Surface.Z)t) は T で変わらない)
        if (symOps.Length > 0)
        {
            var outPos = new Dictionary<int, int>(nOut); for (int o = 0; o < nOut; o++) outPos[outPool[o]] = o;
            var oMap = new int[symOps.Length][]; var oPh = new Complex[symOps.Length][];
            for (int k = 0; k < symOps.Length; k++)
            {
                oMap[k] = new int[nOut]; oPh[k] = new Complex[nOut]; var (w1, w2, w3) = symOps[k].W;
                for (int o = 0; o < nOut; o++)
                {
                    int j = outPos[symOps[k].BeamMap[outPool[o]]]; oMap[k][o] = j;
                    var (h, kk, l) = pool.Hkl[outPool[j]];
                    oPh[k][o] = Complex.Exp(new Complex(0, 2 * Math.PI * (h * w1 + kk * w2 + l * w3))); // 振幅は U の比の共役で移る (行列の (g, h) が U_{h−g} の規約。各方向を自分の区画で解いた対称なし/ありの差 6e−7、day13s)
                }
            }
            var full = new Complex[nAll][];
            for (int r = 0; r < nDir; r++) full[repOrig[r]] = amplitudes[r];
            for (int i = 0; i < nAll; i++)
            {
                var (r, k) = mapFrom[i];
                if (r < 0 || amplitudes[r] == null) continue;
                var src = amplitudes[r]; var dst = new Complex[src.Length]; var om = oMap[k]; var ph = oPh[k];
                for (int t = 0; t < tLen; t++) for (int o = 0; o < nOut; o++) dst[t * nOut + om[o]] = src[t * nOut + o] * ph[o];
                full[i] = dst;
            }
            amplitudes = full;
        }
        string[] propName = ["native eigen", "native matrix exponential", "MathNet eigen"];
        CbedBeamSelectionDetail = $"dressed tail {CbedDressedTailSettings.Version}: N {set.N}, pool {pool.Count} within {set.PoolRadiusNm} nm⁻¹, eta {set.Eta}" +
            (double.IsNaN(set.NeighborRadiusNm) ? "" : $", core {set.CoreRadiusNm} + neighbours {set.NeighborRadiusNm} nm⁻¹") + (set.OneProduct ? ", one product" : "") + "; " +
            (adaptive ? $"adaptive patches (radius {set.AdaptiveMaxRadiusMrad} down to {set.PatchRadiusMrad} mrad, indicator limit {set.AdaptiveIndicatorLimit}): {nSolvedPatches} solved, {nSplit} split, {nLevels} rounds" : $"{nSolvedPatches} patches (radius {set.PatchRadiusMrad} mrad)") + $" for {nDir} directions, {nBackward} re-solved alone; explicit {minS}-{maxS} (mean {(nPatchesDone > 0 ? sumS / nPatchesDone : 0)}), tail mean {(nPatchesDone > 0 ? sumQ / nPatchesDone : 0)}, dressed modes mean {(nPatchesDone > 0 ? sumNb / nPatchesDone : 0)}, basis mean {(nPatchesDone > 0 ? sumB / nPatchesDone : 0)} (mode tau {set.ModeTau}{(double.IsFinite(set.ModeReduceTau) ? $", reduced by {set.ModeReduceTau}" : "")}), indicator {etaMin:e2}-{etaMax:e2}{(set.PromoteOverPatch ? ", promoted over patch edges" : "")}, reduction skipped in {nReduceSkipped} patches, modes kept for crossings {nCrossKept}, modes kept for slow beams {nSlowKept}, probe check max |dI| {maxProbeDiff:e1}, directions re-solved by the guard {nGuardResolved}; " +
            $"output beams {nOut}; propagation {propName[prop]}, dense products {(MathNetProviderManager.MklActive ? "MKL" : "managed (MKL not loaded: about 4 times slower)")}; time (patch-s = sum of patch wall times, includes inner parallel work): tables {stage[0]:f1} s, output-beam selection {stage[1]:f1} s, dressing {stage[2]:f1} patch-s (select {stage[11]:f1}, blocks {stage[12]:f1}, X0 {stage[13]:f1}, eig {stage[14]:f1}, rest {stage[15]:f1}, spmm {stage[16]:f1}+{stage[17]:f1}), basis {stage[3]:f1} patch-s (Gram {stage[8]:f1}, pencil {stage[9]:f1}, reduction {stage[10]:f1}), " +
            $"per direction form {stage[4]:f1} / propagate {stage[5]:f1} / output {stage[6]:f1} patch-s, total {stage[7]:f1} s; parallel {outer} patches × {inner}" +
            (symOps.Length > 0 ? $"; symmetry {string.Join(" ", symOps.Select(o => o.Name))}: {nDir} of {nAll} directions solved" : "; symmetry not used"); // 261111Cl 追加
        return (beams, amplitudes);
    }

    /// <summary>261111Cl 追加: CBED の円盤の対称の操作 (晶帯軸)。ビーム軸 (実験室の z) まわりの 2 回・4 回の回転と、x・y・対角の鏡映 (画面の方向の正方格子を正確に移す 7 つ) のうち、
    /// 表面の法線を保ち、方向の組と候補のビームの組を移し (結晶の向きのずれ 1e−9 以内)、ポテンシャルが並進の位相つきで移る
    /// (U_{Tg} = U_g exp(−2πi h_{Tg}·w)、U' も同じ。w (分率) を 1/48 刻みで探し、|g| ≤ 25 nm⁻¹ の全部で 1e−7 の精度で照合) ものを返す。
    /// 振幅は共役の位相で移る: ψ_{Tg}(Td) = ψ_g(d) exp(+2πi h_{Tg}·w) (Bloch の方程式が対角の位相の相似変換で移るため。行列の (g, h) 成分が U_{h−g} の規約で共役になる)。群として閉じていなければ何も返さない。</summary>
    // private (string Name, int[] DirMap, int[] BeamMap, (double, double, double) W)[] CbedSymmetryOps(BeamSelectionPool pool, IReadOnlyList<Vector3DBase> directions, Vector3DBase n) // 261111Cl 変更前
    // 261111Cl 変更: 候補のプールに依らない形 (従来の CBED のビームの組でも使う)。hkl[0] は 000、g = mat · hkl (実験室)、indexOf は hkl の番号 (無ければ −1)、checkRadius 以下の g でポテンシャルを照合
    private (string Name, int[] DirMap, int[] BeamMap, (double, double, double) W)[] CbedSymmetryOps(IReadOnlyList<(int h, int k, int l)> hkl, double[] gx, double[] gy, double[] gz,
        Func<(int, int, int), int> indexOf, (double, double, double, double, double, double, double, double, double) mat, double checkRadius, IReadOnlyList<Vector3DBase> directions, Vector3DBase n)
    {
        var cands = new (string name, int a, int b, int c, int d)[] { ("2", -1, 0, 0, -1), ("4", 0, -1, 1, 0), ("4^3", 0, 1, -1, 0), ("m(x)", -1, 0, 0, 1), ("m(y)", 1, 0, 0, -1), ("m(x=y)", 0, 1, 1, 0), ("m(x=-y)", 0, -1, -1, 0) };
        var found = new List<(string, int[], int[], (double, double, double), (int, int, int, int))>();
        int nDir = directions.Count;
        static long Key(double v) => (long)Math.Round(v * 1e12);
        var dirKey = new Dictionary<(long, long), int>(nDir);
        for (int i = 0; i < nDir; i++) dirKey.TryAdd((Key(directions[i].X), Key(directions[i].Y)), i);
        var (m11, m12, m13, m21, m22, m23, m31, m32, m33) = mat; // 実験室の g = mat · hkl
        var M = new Matrix3D(m11, m21, m31, m12, m22, m32, m13, m23, m33);
        var Minv = M.Inverse();
        int count = hkl.Count;
        // var chk = Enumerable.Range(1, count - 1).Where(i => gx[i] * gx[i] + gy[i] * gy[i] + gz[i] * gz[i] <= checkRadius * checkRadius).ToArray(); // 261111Cl 変更前 (ビームの組の中で照合)
        // var uCache = new Dictionary<int, (Complex Real, Complex Imag)>();
        // 261111Cl 変更: ポテンシャルの照合の組は、ビームの組に依らない |g| ≤ checkRadius の 3 次元の逆格子の点 (多くても約 6000 点に半径を抑える)。
        //   従来の N 本のビームの組は |g| ≤ 25 が ZOLZ ばかりで、並進 w のビーム軸方向の成分が決まらず、HOLZ の層ごとに振幅の位相がずれた (Si [110] 200 kV N200 の FOLZ で e^{−i 105°}、強度は一致)
        var mi = Minv.Tuple;
        double vol = 1 / Math.Abs(m11 * (m22 * m33 - m23 * m32) - m12 * (m21 * m33 - m23 * m31) + m13 * (m21 * m32 - m22 * m31)); // 単位胞の体積 [nm³] (逆格子の胞の体積の逆数)
        double rChk = Math.Min(checkRadius, Math.Cbrt(6000 / (4.19 * vol)));
        int hb = (int)Math.Ceiling(rChk * Math.Sqrt(mi.E11 * mi.E11 + mi.E12 * mi.E12 + mi.E13 * mi.E13)), kb = (int)Math.Ceiling(rChk * Math.Sqrt(mi.E21 * mi.E21 + mi.E22 * mi.E22 + mi.E23 * mi.E23)),
            lb = (int)Math.Ceiling(rChk * Math.Sqrt(mi.E31 * mi.E31 + mi.E32 * mi.E32 + mi.E33 * mi.E33));
        var chk = new List<(int, int, int)>();
        for (int h = -hb; h <= hb; h++) for (int k = -kb; k <= kb; k++) for (int l = -lb; l <= lb; l++)
        {
            if (h == 0 && k == 0 && l == 0) continue;
            double x = m11 * h + m12 * k + m13 * l, y = m21 * h + m22 * k + m23 * l, z = m31 * h + m32 * k + m33 * l;
            if (x * x + y * y + z * z <= rChk * rChk) chk.Add((h, k, l));
        }
        var uCache = new Dictionary<(int, int, int), (Complex Real, Complex Imag)>();
        (Complex Real, Complex Imag) U((int h, int k, int l) q)
        {
            if (!uCache.TryGetValue(q, out var u)) uCache[q] = u = getU(AccVoltage, q, new Vector3DBase(m11 * q.h + m12 * q.k + m13 * q.l, m21 * q.h + m22 * q.k + m23 * q.l, m31 * q.h + m32 * q.k + m33 * q.l));
            return u;
        }
        foreach (var (name, a, b, c, d) in cands)
        {
            if (Math.Abs(a * n.X + b * n.Y - n.X) + Math.Abs(c * n.X + d * n.Y - n.Y) > 1e-12) continue; // 表面の法線
            // 方向の組
            var dm = new int[nDir]; bool ok = true;
            for (int i = 0; i < nDir && ok; i++)
            {
                var v = directions[i]; double x = a * v.X + b * v.Y, y = c * v.X + d * v.Y;
                if (!dirKey.TryGetValue((Key(x), Key(y)), out var j) || Math.Abs(directions[j].X - x) + Math.Abs(directions[j].Y - y) + Math.Abs(directions[j].Z - v.Z) > 1e-12) ok = false;
                else dm[i] = j;
            }
            if (!ok) continue;
            // 候補のビーム (hkl の整数の行列 K = Mat⁻¹ T Mat)
            var K = Minv * new Matrix3D(a, c, 0, b, d, 0, 0, 0, 1) * M;
            var kt = K.Tuple; var ke = new[] { kt.E11, kt.E12, kt.E13, kt.E21, kt.E22, kt.E23, kt.E31, kt.E32, kt.E33 };
            if (ke.Any(e => Math.Abs(e - Math.Round(e)) > 1e-6)) continue;
            var ki = ke.Select(e => (int)Math.Round(e)).ToArray();
            var bm = new int[count];
            for (int i = 0; i < count && ok; i++)
            {
                var (h, k, l) = hkl[i];
                int j = indexOf((ki[0] * h + ki[1] * k + ki[2] * l, ki[3] * h + ki[4] * k + ki[5] * l, ki[6] * h + ki[7] * k + ki[8] * l));
                double tx = a * gx[i] + b * gy[i], ty = c * gx[i] + d * gy[i], tz = gz[i];
                if (j < 0 || Math.Abs(gx[j] - tx) + Math.Abs(gy[j] - ty) + Math.Abs(gz[j] - tz) > 1e-9 * (1 + Math.Abs(tx) + Math.Abs(ty) + Math.Abs(tz))) ok = false;
                else bm[i] = j;
            }
            if (!ok || chk.Count < 3) continue;
            // ポテンシャル: 並進 w の位相を、|U| の大きい 16 本の比 U_{Tg}/U_g から 1/48 刻みで探す
            (int, int, int) Kq((int h, int k, int l) q) => (ki[0] * q.h + ki[1] * q.k + ki[2] * q.l, ki[3] * q.h + ki[4] * q.k + ki[5] * q.l, ki[6] * q.h + ki[7] * q.k + ki[8] * q.l);
            double uMax = chk.Max(q => U(q).Real.Magnitude), upMax = chk.Max(q => U(q).Imag.Magnitude);
            if (!(uMax > 0)) continue;
            var fit = chk.Where(q => U(q).Real.Magnitude >= 1e-3 * uMax).OrderByDescending(q => U(q).Real.Magnitude).Take(16).ToArray();
            const int G = 48;
            var rr = fit.Select(q => U(Kq(q)).Real / U(q).Real).ToArray();
            var eh = new Complex[fit.Length, 3, G];
            for (int f = 0; f < fit.Length; f++)
            {
                var (h, k, l) = Kq(fit[f]); var hk = new[] { h, k, l };
                for (int ax = 0; ax < 3; ax++) for (int q = 0; q < G; q++) eh[f, ax, q] = Complex.Exp(new Complex(0, -2 * Math.PI * hk[ax] * q / (double)G));
            }
            double bestCost = double.PositiveInfinity; (int, int, int) bq = (0, 0, 0);
            for (int q1 = 0; q1 < G; q1++) for (int q2 = 0; q2 < G; q2++) for (int q3 = 0; q3 < G; q3++)
            {
                double cost = 0;
                for (int f = 0; f < fit.Length && cost < bestCost; f++) cost += (rr[f] - eh[f, 0, q1] * eh[f, 1, q2] * eh[f, 2, q3]).Magnitude;
                if (cost < bestCost) { bestCost = cost; bq = (q1, q2, q3); }
            }
            var w = (bq.Item1 / (double)G, bq.Item2 / (double)G, bq.Item3 / (double)G);
            foreach (var q in chk)
            {
                var tq = Kq(q); var (h, k, l) = tq;
                var ph = Complex.Exp(new Complex(0, -2 * Math.PI * (h * w.Item1 + k * w.Item2 + l * w.Item3)));
                if ((U(tq).Real - U(q).Real * ph).Magnitude > 1e-7 * uMax || (U(tq).Imag - U(q).Imag * ph).Magnitude > 1e-7 * upMax) { ok = false; break; }
            }
            if (ok) found.Add((name, dm, bm, w, (a, b, c, d)));
        }
        // 群として閉じているか (2 つの合成が恒等か、見つけた操作のどれか)
        foreach (var x in found)
            foreach (var y in found)
            {
                var (a1, b1, c1, d1) = x.Item5; var (a2, b2, c2, d2) = y.Item5;
                var z = (a1 * a2 + b1 * c2, a1 * b2 + b1 * d2, c1 * a2 + d1 * c2, c1 * b2 + d1 * d2);
                if (z != (1, 0, 0, 1) && !found.Any(f => f.Item5 == z)) return [];
            }
        return [.. found.Select(f => (f.Item1, f.Item2, f.Item3, f.Item4))];
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
        // var rq = oneProduct ? cr.R : HQS * AS + cr.EXt; // H_QS A_S + E A_Q // 261110Cl 変更前
        var rq = oneProduct ? cr.R : (cr.HqsVb ?? HQS * AS) + cr.EXt; // H_QS A_S + E A_Q (261110Cl: 打ち切りの経路では DressCore の T = H_QS V_k を使う)
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

    /// <summary>261110Cl 追加 (codex2 第 37 回): ‖A B − I‖_max (逆行列の健全性)</summary>
    private static double InverseResidual(Matrix<Complex> A, Matrix<Complex> B)
    {
        var p = A * B; double r = 0;
        for (int j = 0; j < p.ColumnCount; j++) for (int i = 0; i < p.RowCount; i++) r = Math.Max(r, (p[i, j] - (i == j ? Complex.One : Complex.Zero)).Magnitude);
        return r;
    }

    /// <summary>261110Cl 追加 (codex2 第 37 回): モードの組の選び分け。固有値が 1e−3 nm⁻¹ より近い組を全組の連結成分でまとめ、組 G の寄与 = 探り・厚さごとの ‖R_G z_G‖ の最大
    /// (R = モードの列、z = モードの係数 (モード × 厚さ)) が全体の最大の τ 倍以上の組を残す</summary>
    private static bool[] KeepClusters(MathNet.Numerics.LinearAlgebra.Vector<Complex> gam, Matrix<Complex> R, IEnumerable<Matrix<Complex>> zs, double tau,
        bool perMode = false) // perMode: 組の寄与の代わりにモードごとの |z_j|·‖R_j‖ の最大 (組の誰かが τ 以上なら組ごと残す、保守的)
    {
        int n = gam.Count;
        var par = Enumerable.Range(0, n).ToArray();
        int Find(int x) { while (par[x] != x) x = par[x] = par[par[x]]; return x; }
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
                if ((gam[i] - gam[j]).Magnitude < 1e-3) { int ri = Find(i), rj = Find(j); if (ri != rj) par[ri] = rj; }
        var groups = Enumerable.Range(0, n).GroupBy(Find).Select(g => g.ToArray()).ToArray();
        var w = new double[groups.Length];
        if (perMode)
        {
            var rn = new double[n]; for (int j = 0; j < n; j++) rn[j] = R.Column(j).L2Norm();
            foreach (var z in zs)
                for (int gi = 0; gi < groups.Length; gi++)
                    foreach (var j in groups[gi])
                        for (int t = 0; t < z.ColumnCount; t++) w[gi] = Math.Max(w[gi], z[j, t].Magnitude * rn[j]);
        }
        else
        foreach (var z in zs)
            for (int gi = 0; gi < groups.Length; gi++)
            {
                var G = groups[gi];
                for (int t = 0; t < z.ColumnCount; t++)
                {
                    double s2 = 0;
                    for (int r = 0; r < R.RowCount; r++)
                    {
                        Complex v = 0;
                        foreach (var j in G) v += R[r, j] * z[j, t];
                        s2 += v.Real * v.Real + v.Imaginary * v.Imaginary;
                    }
                    w[gi] = Math.Max(w[gi], Math.Sqrt(s2));
                }
            }
        double wmax = w.Length > 0 ? w.Max() : 0;
        var keep = new bool[n];
        for (int gi = 0; gi < groups.Length; gi++) if (w[gi] >= tau * wmax) foreach (var j in groups[gi]) keep[j] = true;
        return keep;
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
