// 261109Cl 新規作成 (研究用): 新しいビーム選択の「方式の核」(仕様書 ReciPro_ビーム選択_新方式_仕様書_261108.md §11) の C# 版の試作。
//   作者 11-09「速さを先に」→ 本体へ入れる前に、方式の核を C# で実装し、Python (tools/BeamSelStudy/py/e1_stage_b_261109.py の wc_variant) と照合して
//   1 方向・1 エネルギーあたりの時間を測る。既定の経路 (Find_gVectors・RunEBSDNew・他の手法) は変えない。研究用の道具だけが呼ぶ。
//   方式の核:
//     (1) 明示 S = E1 の選択 (SelectIndicesResearch) + 遅い反射 (|Q/P| < ps、裾の半径の中) の繰り上げ。000 が先頭。
//     (2) 裾 Q = 候補 (|g| ≤ 裾の半径・P > 0) のうち S 以外の全部 (窓なし)。
//     (3) 着衣: 流束座標 (u = P^{1/2}ψ) で X₀ = −Λ⁻¹H_QS (正則化 δ = 1/(πT))、H_b = H_SS + H_SQ X₀ の固有分解 (γ_j, v_j)、
//         X = X̃ V⁻¹、X̃_wj = −[(H_QS + E X₀) V]_wj (λ_w − γ_j)* / (|λ_w − γ_j|² + δ²) (Sylvester の 1 段。E = 裾どうしの疎な結合 |U| ≥ η max)。
//     (4) Galerkin 射影: R = [I; X]、K = R†HR、G = R†R = I + X†X、G⁻¹K の固有分解 (γ, C)、入射の初期値 a = C⁻¹G⁻¹e₀ (= 直交射影)。
//     (5) 出力: 明示 + 裾の物理振幅 ψ_g(z) = Σ_j M_gj e^{2πiγ_j z} (ψ_000(0) ≈ 1)。
//   結合の表 (DressedTailTable) は方向に依らない: 候補全体の CSR (行 g・列 h = U(h − g)、|U| ≥ η max、U = 弾性 + i·吸収)。方向ごとは P^{−1/2} の対角の倍率だけ。
//   規約は本体の行列 (getPotentialMatrix: 要素 (行 g, 列 h) = getU(h − g)、対角 = i·Imag(000)、Q/P を対角に) と同じ。
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Numerics;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics;
using System.Runtime.Intrinsics.X86;
using System.Threading; // 261109Cl 追加 (Interlocked)
using System.Threading.Tasks;
using MathNet.Numerics.LinearAlgebra;
using CDense = MathNet.Numerics.LinearAlgebra.Complex.DenseMatrix;

namespace Crystallography;

public partial class BetheMethod
{
    /// <summary>261109Cl 追加 (研究用): 着衣つきの解の不変データ (結晶・電圧・半径・η ごと。方向に依らない)。</summary>
    public sealed class DressedTailTable
    {
        internal BetheMethod Owner;
        /// <summary>候補のプール (E1 の選択と共有)</summary>
        public BeamSelectionPool Pool { get; internal set; }
        /// <summary>電圧 [kV]</summary>
        public double KV { get; internal set; }
        /// <summary>裾どうしの結合の閾値 (|U| ≥ η·max|U_{d≠0}|)</summary>
        public double Eta { get; internal set; }
        /// <summary>裾の半径 [nm⁻¹] (既定 = プールの半径)</summary>
        public double TailRadiusNm { get; internal set; }
        internal Complex[] diffUfull; // プールの差の表と同じ箱: U_d = Real + i·Imag (弾性 + i·吸収) [nm⁻²]
        /// <summary>000 の吸収 (getU(000).Imag、対角に i 倍で入る)</summary>
        public Complex U000Imag { get; internal set; }
        /// <summary>差の |U| の最大 (d ≠ 0)</summary>
        public double MaxAbsU { get; internal set; }
        internal int[] rowPtr, colIdx; internal double[] vRe, vIm; // CSR (候補の番号): 要素 (i, j) = U(g_j − g_i)、i ≠ j、|U| ≥ η max
        internal float[] vReF, vImF; // 261109Cl 追加: 単精度の疎な積の値 (初めて使うときに作る)
        /// <summary>結合の表の非零の数</summary>
        public long Nonzeros => rowPtr[^1];
        public double SecondsUfull, SecondsCsr;
        // 261109Cl 追加: 後方散乱の源の行列 F(行 g, 列 h) = conj(U'_back(g − h)) (CreateMasterPatternMuBack と同じ値、後方半球 π/2〜π の吸収の積分)。
        //   差の箱 diffF (d = 0 を含む) と、候補全体の CSR (要素 (i, j) = F(g_i − g_j)、|F| ≥ SourceEta·|F(0)|、対角を含む)。CreateDressedTailTable に sourceEta を渡したときだけ作る
        internal Complex[] diffF;
        internal int[] fRowPtr, fColIdx; internal double[] fRe, fIm; internal float[] fReF, fImF;
        /// <summary>261109Cl 追加: 源の行列の疎化の閾値 (|F(d)| ≥ SourceEta·|F(0)|。0 = 差の箱の全部、NaN = 源の表なし)</summary>
        public double SourceEta { get; internal set; } = double.NaN;
        /// <summary>261109Cl 追加: 源の行列の非零の数 (表なしなら 0)</summary>
        public long SourceNonzeros => fRowPtr?[^1] ?? 0;
        public double SecondsSource;
    }

    /// <summary>261109Cl 追加 (研究用): 着衣つきの解 1 方向分。</summary>
    public sealed class DressedTailSolution
    {
        /// <summary>候補のプール (PoolIndices の hkl を引く)</summary>
        public BeamSelectionPool Pool;
        /// <summary>明示 S と裾 Q の候補の番号 (先頭 NExplicit 本が S、先頭は 000)</summary>
        public int[] PoolIndices;
        /// <summary>明示の本数 (繰り上げを含む)</summary>
        public int NExplicit;
        /// <summary>E1 で選んだ本数 (000 を含む、繰り上げの前)</summary>
        public int NSelected;
        /// <summary>縮約した固有値 γ_j (N_explicit 個)</summary>
        public Complex[] Gamma;
        /// <summary>物理振幅 M (行 = PoolIndices の順、列 = モード)、列優先 (行数 = PoolIndices.Length)。原子位置の振幅だけを求めたときは null</summary>
        public Complex[] Amp;
        /// <summary>261109Cl 追加: 原子位置の振幅 Φ_aj = Σ_g e^{2πi g·r_a} M_gj (行 = 原子、列 = モード、列優先)。共有の経路で atomSites を渡したときだけ</summary>
        public Complex[] AtomAmp;
        /// <summary>261109Cl 追加: AtomAmp の原子の数</summary>
        public int AtomCount;
        /// <summary>261109Cl 追加: モードの基底の源の行列 B = Z† F_W Z (F_W = W†FW、モード × モード、列優先)。共有の経路で backscatterSource を指定したときだけ。
        /// native の EBSDSolverWithTDS に固有ベクトル = 単位行列・α = 1 と一緒に、muBack (CreateMasterPatternMuBack と同じ並び) として渡す</summary>
        public Complex[] SourceMatrix;
    }

    /// <summary>261109Cl 追加 (研究用): 段ごとの時間 [s]</summary>
    public sealed class DressedTailTiming
    {
        public double Select, Blocks, Dress, SpMM1, SpMM2, Project, Eig;
        public double Total => Select + Blocks + Dress + SpMM1 + SpMM2 + Project + Eig;
    }

    /// <summary>261109Cl 追加 (研究用): 方向に依らない結合の表を作る (結晶・電圧・半径・η ごとに 1 回)。電圧は今の AccVoltage (0 以下なら 20 kV)。</summary>
    /// <param name="pool"><see cref="CreateBeamSelectionPool"/> のプール (候補と差の箱を共有)</param>
    /// <param name="eta">裾どうしの結合の閾値 (研究の既定 0.02、通常用の候補 0.05)</param>
    /// <param name="sourceEta">261109Cl 追加: 後方散乱の源の行列の表も作るときの閾値 (|F(d)| ≥ sourceEta·|F(0)|、0 = 差の箱の全部。NaN = 作らない)</param>
    // public DressedTailTable CreateDressedTailTable(BeamSelectionPool pool, double eta) // 261109Cl 旧シグネチャ
    public DressedTailTable CreateDressedTailTable(BeamSelectionPool pool, double eta, double sourceEta = double.NaN)
    {
        if (pool.Owner != this) throw new ArgumentException("The pool was built for another BetheMethod", nameof(pool));
        var sw = Stopwatch.StartNew();
        double kV = AccVoltage > 0 ? AccVoltage : 20;
        var (m11, m12, m13, m21, m22, m23, m31, m32, m33) = pool.Mat;
        var (allowed, _) = CenteringExtinctionRule(Crystal, true);
        int dh = pool.DiffBoundH, dk = pool.DiffBoundK, dl = pool.DiffBoundL, dnk = 2 * dk + 1, dnl = 2 * dl + 1;
        var full = new Complex[(2 * dh + 1) * dnk * dnl];
        double r2 = pool.RadiusNm * pool.RadiusNm * (1 + 1e-9), d2max = 4 * r2;
        // 全ての U (弾性 + i·吸収) を差の箱で求める。getU (cache: false) をそのまま使う (式を複製しない)
        Parallel.For(-dh, dh + 1, h =>
        {
            for (int k = -dk; k <= dk; k++)
                for (int l = -dl; l <= dl; l++)
                {
                    if ((h, k, l) == (0, 0, 0) || !allowed(h, k, l)) continue;
                    double x = Dot3(m11, h, m12, k, m13, l), y = Dot3(m21, h, m22, k, m23, l), z = Dot3(m31, h, m32, k, m33, l);
                    if (Dot3(x, x, y, y, z, z) > d2max) continue;
                    var (re, im) = getU(kV, new Beam((h, k, l), new Vector3DBase(x, y, z)), cache: false);
                    full[((h + dh) * dnk + k + dk) * dnl + l + dl] = re + Complex.ImaginaryOne * im;
                }
        });
        double maxAbs = 0;
        for (int i = 0; i < full.Length; i++) maxAbs = Math.Max(maxAbs, full[i].Magnitude);
        var table = new DressedTailTable
        {
            Owner = this, Pool = pool, KV = kV, Eta = eta, TailRadiusNm = pool.RadiusNm, diffUfull = full, MaxAbsU = maxAbs,
            U000Imag = getU(kV, new Beam((0, 0, 0), new Vector3DBase(0, 0, 0)), cache: false).Imag,
        };
        table.SecondsUfull = sw.Elapsed.TotalSeconds;
        // CSR: 行 g_i、列 g_j = g_i + d (|U_d| ≥ η max、d ≠ 0)。d の候補を先に絞ってから行ごとに引く
        sw.Restart();
        double thr = eta * maxAbs;
        var ds = new List<(int h, int k, int l, Complex u)>();
        for (int h = -dh; h <= dh; h++)
            for (int k = -dk; k <= dk; k++)
                for (int l = -dl; l <= dl; l++)
                {
                    var u = full[((h + dh) * dnk + k + dk) * dnl + l + dl];
                    if ((h, k, l) != (0, 0, 0) && u.Magnitude >= thr) ds.Add((h, k, l, u));
                }
        int n = pool.Count;
        var rows = new (int j, Complex u)[n][];
        Parallel.For(0, n, i =>
        {
            var (h0, k0, l0) = pool.Hkl[i];
            var list = new List<(int, Complex)>();
            foreach (var d in ds)
            {
                int j = pool.IndexOf((h0 + d.h, k0 + d.k, l0 + d.l));
                if (j >= 0) list.Add((j, d.u)); // 要素 (i, j) = U(g_j − g_i) = U_d
            }
            list.Sort((a, b) => a.Item1.CompareTo(b.Item1));
            rows[i] = [.. list];
        });
        var rowPtr = new int[n + 1];
        for (int i = 0; i < n; i++) rowPtr[i + 1] = rowPtr[i] + rows[i].Length;
        var col = new int[rowPtr[n]]; var vr = new double[rowPtr[n]]; var vi = new double[rowPtr[n]];
        Parallel.For(0, n, i =>
        {
            int p = rowPtr[i];
            foreach (var (j, u) in rows[i]) { col[p] = j; vr[p] = u.Real; vi[p] = u.Imaginary; p++; }
        });
        table.rowPtr = rowPtr; table.colIdx = col; table.vRe = vr; table.vIm = vi;
        table.SecondsCsr = sw.Elapsed.TotalSeconds;
        if (!double.IsNaN(sourceEta)) // 261109Cl 追加: 後方散乱の源の行列 (CreateMasterPatternMuBack と同じ getU の積分、30 × 12 点)
        {
            sw.Restart();
            var fFull = new Complex[full.Length];
            Parallel.For(-dh, dh + 1, h =>
            {
                for (int k = -dk; k <= dk; k++)
                    for (int l = -dl; l <= dl; l++)
                    {
                        if (!allowed(h, k, l)) continue;
                        double x = Dot3(m11, h, m12, k, m13, l), y = Dot3(m21, h, m22, k, m23, l), z = Dot3(m31, h, m32, k, m33, l);
                        if (Dot3(x, x, y, y, z, z) > d2max) continue;
                        fFull[((h + dh) * dnk + k + dk) * dnl + l + dl] = Complex.Conjugate(getU(kV, new Beam((h, k, l), new Vector3DBase(x, y, z)), null, Math.PI / 2, Math.PI, 30, 12, cache: false).Imag);
                    }
            });
            int Box(int h, int k, int l) => ((h + dh) * dnk + k + dk) * dnl + l + dl;
            double fThr = sourceEta * fFull[Box(0, 0, 0)].Magnitude;
            // 差 d (列 g_j = g_i + d) の値は F(g_i − g_j) = F(−d)。d と −d を対で残す (|F(d)| と |F(−d)| の大きい方で判定)
            var fds = new List<(int h, int k, int l, Complex f)>();
            for (int h = -dh; h <= dh; h++)
                for (int k = -dk; k <= dk; k++)
                    for (int l = -dl; l <= dl; l++)
                    {
                        if (!allowed(h, k, l)) continue;
                        double x = Dot3(m11, h, m12, k, m13, l), y = Dot3(m21, h, m22, k, m23, l), z = Dot3(m31, h, m32, k, m33, l);
                        if (Dot3(x, x, y, y, z, z) > d2max) continue;
                        var fm = fFull[Box(-h, -k, -l)];
                        if (Math.Max(fm.Magnitude, fFull[Box(h, k, l)].Magnitude) >= fThr) fds.Add((h, k, l, fm));
                    }
            var fRows = new (int j, Complex f)[n][];
            Parallel.For(0, n, i =>
            {
                var (h0, k0, l0) = pool.Hkl[i];
                var list = new List<(int, Complex)>();
                foreach (var d in fds)
                {
                    int j = pool.IndexOf((h0 + d.h, k0 + d.k, l0 + d.l));
                    if (j >= 0) list.Add((j, d.f));
                }
                list.Sort((a, b) => a.Item1.CompareTo(b.Item1));
                fRows[i] = [.. list];
            });
            var fPtr = new int[n + 1];
            long fTotal = 0;
            for (int i = 0; i < n; i++) { fTotal += fRows[i].Length; if (fTotal > int.MaxValue) throw new OutOfMemoryException("The backscatter source table exceeds int.MaxValue nonzeros; raise sourceEta"); fPtr[i + 1] = (int)fTotal; }
            var fCol = new int[fTotal]; var fr = new double[fTotal]; var fi = new double[fTotal];
            Parallel.For(0, n, i =>
            {
                int p = fPtr[i];
                foreach (var (j, f) in fRows[i]) { fCol[p] = j; fr[p] = f.Real; fi[p] = f.Imaginary; p++; }
            });
            table.diffF = fFull; table.fRowPtr = fPtr; table.fColIdx = fCol; table.fRe = fr; table.fIm = fi; table.SourceEta = sourceEta;
            table.SecondsSource = sw.Elapsed.TotalSeconds;
        }
        return table;
    }

    /// <summary>261109Cl 追加 (研究用): 1 方向の着衣つきの解 (方式の核)。</summary>
    /// <param name="table"><see cref="CreateDressedTailTable"/> の表</param>
    /// <param name="vecK0">屈折後の k0</param>
    /// <param name="surface">表面の法線</param>
    /// <param name="n">E1 で選ぶ本数 (000 を含む)</param>
    /// <param name="ps">遅い反射の繰り上げの閾値 |Q/P| [nm⁻¹]</param>
    /// <param name="regDepthNm">着衣の正則化 δ = 1/(πT) の T [nm] (研究の既定 45)</param>
    /// <param name="oneProduct">射影の E·X を着衣前の E·X₀ で近似し、疎な積を 1 回にする (精度は DW で約 0.1 % 落ちる)</param>
    /// <param name="parallel">疎な積を並列にする (時間の計測では false = 1 スレッド)</param>
    public DressedTailSolution SolveDressedTailResearch(DressedTailTable table, Vector3DBase vecK0, Vector3DBase surface, int n, double ps = 0.1,
        double regDepthNm = 45, bool oneProduct = false, bool parallel = false, DressedTailTiming timing = null)
    {
        var cr = DressCore(table, vecK0, surface, n, ps, regDepthNm, oneProduct, parallel, timing);
        var (sArr, qArr, P, D, X, HSS, HQS, HSQ, EX, lamQ) = (cr.SArr, cr.QArr, cr.P, cr.D, cr.X, cr.HSS, cr.HQS, cr.HSQ, cr.EX, cr.LamQ);
        int nS = sArr.Length, nQ = qArr.Length;
        var sw = Stopwatch.StartNew();
        // (4) Galerkin 射影: K = H_SS + H_SQ X + X†(H_QS + Λ X + E X)、G = I + X†X
        var hrq = (HQS + EX).ToColumnMajorArray();
        var xv = ((CDense)X).Values;
        for (int c = 0; c < nS; c++)
            for (int r = 0; r < nQ; r++)
                hrq[r + c * nQ] += lamQ[r] * xv[r + c * nQ];
        var HRQ = new CDense(nQ, nS, hrq);
        var K = HSS + HSQ * X + X.ConjugateTransposeThisAndMultiply(HRQ);
        var G = X.ConjugateTransposeThisAndMultiply(X);
        for (int i = 0; i < nS; i++) G[i, i] += 1;
        var ch = G.Cholesky();
        var W = ch.Solve(K); // G⁻¹K (G はエルミート正定値)
        if (timing != null) timing.Project = sw.Elapsed.TotalSeconds;
        sw.Restart();
        var evd = W.Evd(Symmetricity.Asymmetric);
        var gam = evd.EigenValues; var C = evd.EigenVectors;
        var e0 = MathNet.Numerics.LinearAlgebra.Vector<Complex>.Build.Dense(nS); e0[0] = 1;
        var a = C.Solve(ch.Solve(e0)); // C a = G⁻¹e₀ (直交射影の初期値)
        if (timing != null) timing.Eig = sw.Elapsed.TotalSeconds;
        // (5) 物理振幅: [C; X C] diag(a) × √(P₀/P_g)
        var Z = C * Matrix<Complex>.Build.DiagonalOfDiagonalVector(a);
        return MakeDressedSolution(table.Pool, sArr, qArr, X, Z, D, Math.Sqrt(P[0]), gam.ToArray(), cr.NSelected);
    }

    /// <summary>261109Cl 追加: 明示と裾の物理振幅 ψ = [D_S ∘ Z; D_Q ∘ (X Z)] × scale (D = 中心の P^{−1/2}) を解の形にする。</summary>
    private static DressedTailSolution MakeDressedSolution(BeamSelectionPool pool, int[] sArr, int[] qArr, Matrix<Complex> X, Matrix<Complex> Z, double[] D, double scale, Complex[] gam, int nSelected)
    {
        int nS = sArr.Length, nQ = qArr.Length, nAll = nS + nQ, nm = Z.ColumnCount;
        var UQ = X * Z;
        var amp = new Complex[nAll * nm];
        var usv = Z.ToColumnMajorArray(); var uqv = UQ.ToColumnMajorArray();
        for (int c = 0; c < nm; c++)
        {
            for (int r = 0; r < nS; r++) amp[r + c * nAll] = usv[r + c * nS] * (scale * D[sArr[r]]);
            for (int r = 0; r < nQ; r++) amp[nS + r + c * nAll] = uqv[r + c * nQ] * (scale * D[qArr[r]]);
        }
        var idxAll = new int[nAll];
        Array.Copy(sArr, idxAll, nS); Array.Copy(qArr, 0, idxAll, nS, nQ);
        return new DressedTailSolution { Pool = pool, PoolIndices = idxAll, NExplicit = nS, NSelected = nSelected, Gamma = gam, Amp = amp };
    }

    /// <summary>261109Cl 追加: 着衣の中間結果 (中心の方向の流束座標)</summary>
    private sealed class DressCoreResult
    {
        public int[] SArr, QArr; public double[] P, D; public int NSelected;
        public Matrix<Complex> X, HSS, HQS, HSQ, EX; public Complex[] LamQ;
    }

    /// <summary>261109Cl 追加 (SolveDressedTailResearch から切り出し): 選択・繰り上げ・密なブロック・着衣 (Sylvester の 1 段) と、射影に使う E·X。</summary>
    private DressCoreResult DressCore(DressedTailTable table, Vector3DBase vecK0, Vector3DBase surface, int n, double ps, double regDepthNm, bool oneProduct, bool parallel, DressedTailTiming timing)
    {
        if (table.Owner != this) throw new ArgumentException("The table was built for another BetheMethod", nameof(table));
        if (table.KV != (AccVoltage > 0 ? AccVoltage : 20)) // 261109Cl 追加 (codex2 第18回): 表の電圧と今の電圧の混在を禁じる
            throw new ArgumentException($"The table was built for {table.KV} kV but the current voltage is {AccVoltage} kV", nameof(table));
        var pool = table.Pool;
        var sw = Stopwatch.StartNew();
        // (1) 選択と繰り上げ
        var e1 = SelectIndicesResearch(pool, vecK0, surface, n, BeamPartnerSet.E1StrongSet, 170, out _, out _);
        var (kX, kY, kZ) = vecK0.Tuple;
        var (sX, sY, sZ) = surface.Tuple;
        double k02 = vecK0.Length2;
        int cnt = pool.Count;
        var Q = new double[cnt]; var P = new double[cnt];
        double rt2 = table.TailRadiusNm * table.TailRadiusNm * (1 + 1e-9);
        var inTail = new bool[cnt];
        for (int i = 0; i < cnt; i++)
        {
            double vX = pool.GX[i] + kX, vY = pool.GY[i] + kY, vZ = pool.GZ[i] + kZ;
            Q[i] = k02 - Dot3(vX, vX, vY, vY, vZ, vZ);
            P[i] = 2 * Dot3(sX, vX, sY, vY, sZ, vZ);
            inTail[i] = P[i] > 0 && Dot3(pool.GX[i], pool.GX[i], pool.GY[i], pool.GY[i], pool.GZ[i], pool.GZ[i]) <= rt2;
        }
        var inS = new bool[cnt];
        var S = new List<int>(e1.Length + 256) { 0 };
        inS[0] = true;
        foreach (var i in e1) { S.Add(i); inS[i] = true; }
        int nSelected = S.Count;
        for (int i = 1; i < cnt; i++)
            if (inTail[i] && !inS[i] && Math.Abs(Q[i] / P[i]) < ps) { S.Add(i); inS[i] = true; }
        var Qi = new List<int>(cnt);
        for (int i = 1; i < cnt; i++)
            if (inTail[i] && !inS[i]) Qi.Add(i);
        int nS = S.Count, nQ = Qi.Count;
        var sArr = S.ToArray(); var qArr = Qi.ToArray();
        var map = new int[cnt]; Array.Fill(map, -1);
        for (int r = 0; r < nQ; r++) map[qArr[r]] = r;
        var D = new double[cnt];
        for (int i = 0; i < cnt; i++) D[i] = P[i] > 0 ? 1 / Math.Sqrt(P[i]) : 0;
        var u0i = Complex.ImaginaryOne * table.U000Imag;
        Complex Lam(int i) => (Q[i] + u0i) / P[i];
        if (timing != null) timing.Select = sw.Elapsed.TotalSeconds;

        // (2) 密なブロック (流束座標): H(g, h) = U(h − g) D_g D_h、H_SS の対角 = λ
        sw.Restart();
        Complex Ufull(int ig, int ih)
        {
            var (h1, k1, l1) = pool.Hkl[ig]; var (h2, k2, l2) = pool.Hkl[ih];
            return table.diffUfull[pool.DiffIndex(h2 - h1, k2 - k1, l2 - l1)];
        }
        var hss = new Complex[nS * nS];
        for (int c = 0; c < nS; c++)
            for (int r = 0; r < nS; r++)
                hss[r + c * nS] = r == c ? Lam(sArr[r]) : Ufull(sArr[r], sArr[c]) * (D[sArr[r]] * D[sArr[c]]);
        var hqs = new Complex[nQ * nS]; var hsq = new Complex[nS * nQ];
        Parallel.For(0, nS, new ParallelOptions { MaxDegreeOfParallelism = parallel ? -1 : 1 }, c =>
        {
            int sc = sArr[c]; double dc = D[sc];
            for (int r = 0; r < nQ; r++)
            {
                int qr = qArr[r]; double dd = D[qr] * dc;
                hqs[r + c * nQ] = Ufull(qr, sc) * dd;      // (行 Q_r, 列 S_c)
                hsq[c + r * nS] = Ufull(sc, qr) * dd;      // (行 S_c, 列 Q_r)
            }
        });
        var HSS = new CDense(nS, nS, hss); var HQS = new CDense(nQ, nS, hqs); var HSQ = new CDense(nS, nQ, hsq);
        var lamQ = new Complex[nQ];
        for (int r = 0; r < nQ; r++) lamQ[r] = Lam(qArr[r]);
        if (timing != null) timing.Blocks = sw.Elapsed.TotalSeconds;

        // (3) 着衣 (Sylvester の 1 段)
        sw.Restart();
        double delta = 1 / (Math.PI * regDepthNm), d2 = delta * delta;
        var x0 = new Complex[nQ * nS];
        for (int c = 0; c < nS; c++)
            for (int r = 0; r < nQ; r++)
            {
                var l = lamQ[r];
                x0[r + c * nQ] = -hqs[r + c * nQ] * Complex.Conjugate(l) / (l.Real * l.Real + l.Imaginary * l.Imaginary + d2);
            }
        var X0 = new CDense(nQ, nS, x0);
        var evdB = (HSS + HSQ * X0).Evd(Symmetricity.Asymmetric);
        var gb = evdB.EigenValues; var Vb = evdB.EigenVectors;
        if (timing != null) timing.Dress = sw.Elapsed.TotalSeconds;
        sw.Restart();
        var EX0 = SpMMTail(table, qArr, map, D, X0, parallel);
        if (timing != null) timing.SpMM1 = sw.Elapsed.TotalSeconds;
        sw.Restart();
        var rhs = (HQS + EX0) * Vb;
        var xt = ((CDense)rhs).Values;
        for (int c = 0; c < nS; c++)
        {
            var g = gb[c];
            for (int r = 0; r < nQ; r++)
            {
                var dd = lamQ[r] - g;
                xt[r + c * nQ] = -xt[r + c * nQ] * Complex.Conjugate(dd) / (dd.Real * dd.Real + dd.Imaginary * dd.Imaginary + d2);
            }
        }
        var X = rhs * Vb.Inverse();
        if (timing != null) timing.Dress += sw.Elapsed.TotalSeconds;
        // 射影に使う E·X (oneProduct なら着衣前の E·X₀ で近似)
        sw.Restart();
        var EX = oneProduct ? EX0 : SpMMTail(table, qArr, map, D, X, parallel);
        if (timing != null) timing.SpMM2 = sw.Elapsed.TotalSeconds;
        return new DressCoreResult { SArr = sArr, QArr = qArr, P = P, D = D, NSelected = nSelected, X = X, HSS = HSS, HQS = HQS, HSQ = HSQ, EX = EX, LamQ = lamQ };
    }

    /// <summary>261109Cl 追加 (研究用): 方向のブロックで着衣を共有する解 (10 月の WCphys と同じ式)。中心の方向で S・Q・着衣 X を作り、
    /// 物理振幅の基底 W = D_c [I; X] (D_c = 中心の P^{−1/2}) の上で G₀ = W†W、G_a = W† diag(g_a) W、C = W†(V − diag|g|²)W を 1 回だけ作る。
    /// 各方向では K = C + iU'₀G₀ − 2Σ k0_a G_a、M = 2(n·k0)G₀ + 2Σ n_a G_a (Q・P が k0・n の一次式なので厳密) とし、M⁻¹K の固有分解だけを解く。
    /// 入射の初期値は P で重みづけた最小二乗: M a(0) = W†P e₀ = (P_000 / √P_c,000) e₀。共有の上乗せの精度は ±2.45 mrad まで確かめた (仕様書 §11 の追記)。</summary>
    /// <param name="dirs">ブロックの各方向の (屈折後の k0、表面の法線)。どれも S∪Q の全ビームが P &gt; 0 であること</param>
    public DressedTailSolution[] SolveDressedTailBlockResearch(DressedTailTable table, Vector3DBase centerK0, Vector3DBase centerSurface, (Vector3DBase k0, Vector3DBase n)[] dirs,
        int n, double ps = 0.1, double regDepthNm = 45, bool oneProduct = false, bool parallel = false, DressedTailTiming timing = null, double[] blockSeconds = null,
        (double x, double y, double z)[] atomSites = null,   // 261109Cl atomSites 追加: 渡すと全振幅の代わりに原子位置の振幅だけを求める (方向ごとの X·Z を中心の 1 回の積に回す)
        bool backscatterSource = false,   // 261109Cl backscatterSource 追加: 後方散乱の源の行列をモードの基底へ射影する (表に源が要る)。旧シグネチャは atomSites までで、これを除いた形
        (double weight, double sigmaNm)[] tailQuadratureAtoms = null, // 261109Cl 追加: 源の裾どうしの項を求積で積むときの原子ごとの (F(0) の分 A_n、Gauss の幅 σ_n = √(B/8π²) [nm])。atomSites と同じ並び
        int tailQuadratureOrder = 0)      // 261109Cl 追加: その Gauss–Hermite の 1 軸の点の数 (0 = 求積せず源の表の疎な積)
    {
        var cr = DressCore(table, centerK0, centerSurface, n, ps, regDepthNm, oneProduct, parallel, timing);
        var pool = table.Pool;
        var (sArr, qArr, X, HSS, HQS, HSQ, EX) = (cr.SArr, cr.QArr, cr.X, cr.HSS, cr.HQS, cr.HSQ, cr.EX);
        double[] Dc = cr.D;
        int nS = sArr.Length, nQ = qArr.Length;
        var sw = Stopwatch.StartNew();
        // Koff = W†VW (流束座標の非対角: H_SS の対角を除く + H_SQ X + X†(H_QS + E X))
        var hssOff = HSS.Clone();
        for (int i = 0; i < nS; i++) hssOff[i, i] = 0;
        var Koff = hssOff + HSQ * X + X.ConjugateTransposeThisAndMultiply(HQS + EX);
        // X† diag(w) X を 5 種 (w = 1, g_x, g_y, g_z, |g|²、どれも × D_c² = 1/P_c) まとめて 1 回の積に
        var xv = ((CDense)X).Values;
        var Y = new CDense(nQ, 5 * nS);
        var yv = Y.Values;
        for (int c = 0; c < nS; c++)
            for (int r = 0; r < nQ; r++)
            {
                int q = qArr[r]; double d2 = Dc[q] * Dc[q];
                double gx = pool.GX[q], gy = pool.GY[q], gz = pool.GZ[q], g2 = Dot3(gx, gx, gy, gy, gz, gz);
                var x = xv[r + c * nQ];
                yv[r + c * nQ] = x * d2;
                yv[r + (c + nS) * nQ] = x * (d2 * gx);
                yv[r + (c + 2 * nS) * nQ] = x * (d2 * gy);
                yv[r + (c + 3 * nS) * nQ] = x * (d2 * gz);
                yv[r + (c + 4 * nS) * nQ] = x * (d2 * g2);
            }
        var T = X.ConjugateTransposeThisAndMultiply(Y);
        Matrix<Complex> Part(int m) => T.SubMatrix(0, nS, m * nS, nS);
        var G0 = Part(0); var Gx = Part(1); var Gy = Part(2); var Gz = Part(3); var Gg2 = Part(4);
        for (int i = 0; i < nS; i++)
        {
            int s = sArr[i]; double d2 = Dc[s] * Dc[s];
            double gx = pool.GX[s], gy = pool.GY[s], gz = pool.GZ[s];
            G0[i, i] += d2; Gx[i, i] += d2 * gx; Gy[i, i] += d2 * gy; Gz[i, i] += d2 * gz; Gg2[i, i] += d2 * Dot3(gx, gx, gy, gy, gz, gz);
        }
        var Cm = Koff - Gg2;
        var iu0 = Complex.ImaginaryOne * table.U000Imag;
        // 原子位置の振幅の前計算: Φ = (Ph_S D_S + Ph_Q D_Q X) Z (Ph = e^{2πi hkl·r}、方向に依らない)。方向ごとは (原子 × 明示) × (明示 × モード) だけ
        Matrix<Complex> phiBase = null;
        if (atomSites != null)
        {
            int na = atomSites.Length;
            var phS = new CDense(na, nS); var phQ = new CDense(na, nQ);
            for (int a = 0; a < na; a++)
            {
                var (ax, ay, az) = atomSites[a];
                for (int i = 0; i < nS; i++) { var (h, k, l) = pool.Hkl[sArr[i]]; phS[a, i] = Complex.FromPolarCoordinates(Dc[sArr[i]], 2 * Math.PI * (h * ax + k * ay + l * az)); }
                for (int i = 0; i < nQ; i++) { var (h, k, l) = pool.Hkl[qArr[i]]; phQ[a, i] = Complex.FromPolarCoordinates(Dc[qArr[i]], 2 * Math.PI * (h * ax + k * ay + l * az)); }
            }
            phiBase = phS + phQ * X;
        }
        // 261109Cl 追加: 後方散乱の源の行列の射影 F_W = W†FW = F'_SS + F'_SQ X + X†(F'_QS + F'_QQ X) (' = 中心の D を両側に掛けたもの。F(行 g, 列 h) = F(g − h))。
        //   方向に依らないので中心で 1 回。各方向は B = Z† F_W Z (明示 × 明示) だけ
        Matrix<Complex> FW = null;
        if (backscatterSource)
        {
            if (table.fRowPtr == null) throw new ArgumentException("The table has no backscatter source; create it with sourceEta", nameof(table));
            Complex Ff(int ig, int ih)
            {
                var (h1, k1, l1) = pool.Hkl[ig]; var (h2, k2, l2) = pool.Hkl[ih];
                return table.diffF[pool.DiffIndex(h1 - h2, k1 - k2, l1 - l2)];
            }
            var fss = new CDense(nS, nS); var fsq = new CDense(nS, nQ); var fqs = new CDense(nQ, nS);
            for (int c = 0; c < nS; c++)
                for (int r = 0; r < nS; r++)
                    fss[r, c] = Ff(sArr[r], sArr[c]) * (Dc[sArr[r]] * Dc[sArr[c]]);
            var fsqv = fsq.Values; var fqsv = fqs.Values;
            Parallel.For(0, nS, new ParallelOptions { MaxDegreeOfParallelism = parallel ? -1 : 1 }, c =>
            {
                int sc = sArr[c]; double dc = Dc[sc];
                for (int r = 0; r < nQ; r++)
                {
                    int qr = qArr[r]; double dd = Dc[qr] * dc;
                    fqsv[r + c * nQ] = Ff(qr, sc) * dd;   // (行 Q_r, 列 S_c)
                    fsqv[c + r * nS] = Ff(sc, qr) * dd;   // (行 S_c, 列 Q_r)
                }
            });
            // 裾どうしの項 X†F'_QQ X。F はほとんど減衰しない (R40 で非零 5 千万〜1 億、η_F 0.1 で 1 方向 +15 s) ので、既定は求積:
            //   F の包絡は DW の Gauss exp(−B|d|²/4) とほぼ同じ (Si 20 kV で比 0.98〜1.00、day5e) → F(g − h) ≈ Σ_n A_n e^{−2πi(g−h)·r_n} ∫G_σn(x) e^{−2πi(g−h)·x} d³x
            //   (σ_n² = B_n/8π²、A_n = 原子 n の F(0) の分)。よって X†F'_QQ X ≈ Σ_n A_n Σ_q w_q χ_n,q† χ_n,q (χ = 裾の場 Σ_g D_g X_g e^{2πi g·(r_n + σ_n t_q)}、t_q = Gauss–Hermite の積の点)
            Matrix<Complex> qq;
            if (tailQuadratureOrder > 0)
            {
                if (atomSites == null || tailQuadratureAtoms == null || tailQuadratureAtoms.Length != atomSites.Length)
                    throw new ArgumentException("The tail quadrature needs atomSites and one (weight, sigma) per site", nameof(tailQuadratureAtoms));
                var (tq, wq) = GaussHermiteProbabilists(tailQuadratureOrder);
                int nq = tq.Length, nPer = nq * nq * nq, na = atomSites.Length, nSites = na * nPer;
                var ph = new CDense(nSites, nQ); var phv = ph.Values;
                Parallel.For(0, nQ, new ParallelOptions { MaxDegreeOfParallelism = parallel ? -1 : 1 }, r =>
                {
                    int q = qArr[r]; var (h, k, l) = pool.Hkl[q];
                    double gx = pool.GX[q], gy = pool.GY[q], gz = pool.GZ[q];
                    Span<Complex> ex = stackalloc Complex[nq], ey = stackalloc Complex[nq], ez = stackalloc Complex[nq];
                    for (int a = 0; a < na; a++)
                    {
                        var (ax, ay, az) = atomSites[a]; double sig = tailQuadratureAtoms[a].sigmaNm;
                        var b0 = Complex.FromPolarCoordinates(Dc[q], 2 * Math.PI * (h * ax + k * ay + l * az));
                        for (int i = 0; i < nq; i++)
                        {
                            ex[i] = Complex.FromPolarCoordinates(1, 2 * Math.PI * sig * gx * tq[i]);
                            ey[i] = Complex.FromPolarCoordinates(1, 2 * Math.PI * sig * gy * tq[i]);
                            ez[i] = Complex.FromPolarCoordinates(1, 2 * Math.PI * sig * gz * tq[i]);
                        }
                        int o = a * nPer;
                        for (int i = 0; i < nq; i++)
                            for (int j = 0; j < nq; j++)
                            {
                                var bij = b0 * ex[i] * ey[j];
                                for (int m = 0; m < nq; m++) phv[o + (i * nq + j) * nq + m + (long)r * nSites] = bij * ez[m];
                            }
                    }
                });
                var chi = ph * X; // (点 × 明示)
                var chiW = chi.Clone();
                for (int a = 0; a < na; a++)
                    for (int i = 0; i < nq; i++)
                        for (int j = 0; j < nq; j++)
                            for (int m = 0; m < nq; m++)
                            {
                                int row = a * nPer + (i * nq + j) * nq + m;
                                double w = tailQuadratureAtoms[a].weight * wq[i] * wq[j] * wq[m];
                                for (int c = 0; c < nS; c++) chiW[row, c] *= w;
                            }
                qq = chi.ConjugateTransposeThisAndMultiply(chiW);
            }
            else
            {
                var mapQ = new int[pool.Count]; Array.Fill(mapQ, -1);
                for (int r = 0; r < nQ; r++) mapQ[qArr[r]] = r;
                var FQX = SpMMTail(table, qArr, mapQ, Dc, X, parallel, source: true);
                qq = X.ConjugateTransposeThisAndMultiply(FQX);
            }
            FW = fss + fsq * X + X.ConjugateTransposeThisAndMultiply(fqs) + qq;
        }
        if (blockSeconds != null) blockSeconds[0] += sw.Elapsed.TotalSeconds;
        var result = new DressedTailSolution[dirs.Length];
        for (int t = 0; t < dirs.Length; t++)
        {
            sw.Restart();
            var (k0, nv) = dirs[t];
            var (kx, ky, kz) = k0.Tuple; var (nx, ny, nz) = nv.Tuple;
            double nk = Dot3(nx, kx, ny, ky, nz, kz);
            // 全ビームが前方 (P > 0) であることを確かめる
            foreach (var i in sArr.Concat(qArr))
                if (2 * Dot3(nx, pool.GX[i] + kx, ny, pool.GY[i] + ky, nz, pool.GZ[i] + kz) <= 0)
                    // throw new InvalidOperationException("A shared beam is not forward (P <= 0) in this direction; shrink the block"); // 261109Cl 変更前
                    throw new DressedTailBackwardBeamException(); // 261109Cl 変更: 呼び出し側がブロックを分ける印 (InvalidOperationException の派生)
            var K = Cm + iu0 * G0 - 2 * (kx * Gx + ky * Gy + kz * Gz);
            var M = 2 * nk * G0 + 2 * (nx * Gx + ny * Gy + nz * Gz);
            var ch = M.Cholesky();
            var evd = ch.Solve(K).Evd(Symmetricity.Asymmetric);
            var Cv = evd.EigenVectors;
            var e0 = MathNet.Numerics.LinearAlgebra.Vector<Complex>.Build.Dense(nS); e0[0] = 2 * nk * Dc[0]; // W†P e₀ = P_000 · D_c,000 e₀ (P_000 = 2 n·k0)
            var a = Cv.Solve(ch.Solve(e0));
            var Z = Cv * Matrix<Complex>.Build.DiagonalOfDiagonalVector(a);
            if (phiBase != null)
            {
                var idxAll = new int[nS + nQ];
                Array.Copy(sArr, idxAll, nS); Array.Copy(qArr, 0, idxAll, nS, nQ);
                result[t] = new DressedTailSolution { Pool = pool, PoolIndices = idxAll, NExplicit = nS, NSelected = cr.NSelected, Gamma = evd.EigenValues.ToArray(),
                    AtomAmp = (phiBase * Z).ToColumnMajorArray(), AtomCount = atomSites.Length };
            }
            else
                result[t] = MakeDressedSolution(pool, sArr, qArr, X, Z, Dc, 1.0, evd.EigenValues.ToArray(), cr.NSelected);
            if (FW != null) result[t].SourceMatrix = Z.ConjugateTransposeThisAndMultiply(FW * Z).ToColumnMajorArray(); // 261109Cl 追加
            if (blockSeconds != null) blockSeconds[1] += sw.Elapsed.TotalSeconds;
        }
        return result;
    }

    /// <summary>261109Cl 追加 (研究用): 裾どうしの疎な結合と密な行列の積 E·X (流束座標: E_wv = U(g_v − g_w) D_w D_v、w, v ∈ Q、w ≠ v)。
    /// X を行優先へ並べ替え、行ごとに CSR の非零 × X の行を AVX の複素の積和で足す。S の列 (map &lt; 0) は飛ばす。</summary>
    // private static Matrix<Complex> SpMMTail(DressedTailTable table, int[] qArr, int[] map, double[] D, Matrix<Complex> X, bool parallel) // 261109Cl 旧シグネチャ
    private static Matrix<Complex> SpMMTail(DressedTailTable table, int[] qArr, int[] map, double[] D, Matrix<Complex> X, bool parallel, bool source = false) // 261109Cl source 追加: 後方散乱の源の表 (F) で積む
    {
        // if (DressedSpMMSingle) return SpMMTailSingle(table, qArr, map, D, X, parallel); // 261109Cl 変更前
        if (DressedSpMMSingle) return SpMMTailSingle(table, qArr, map, D, X, parallel, source); // 261109Cl
        var sw = Stopwatch.StartNew();
        int nQ = X.RowCount, k = X.ColumnCount;
        var xcm = ((CDense)X).Values;
        var xr = new Complex[nQ * k]; // 行優先、行に D_v を掛けておく (32×32 のタイルで並べ替える)
        for (int r0 = 0; r0 < nQ; r0 += 32)
            for (int c0 = 0; c0 < k; c0 += 32)
                for (int r = r0, rEnd = Math.Min(r0 + 32, nQ); r < rEnd; r++)
                {
                    double dv = D[qArr[r]];
                    for (int c = c0, cEnd = Math.Min(c0 + 32, k); c < cEnd; c++)
                        xr[r * k + c] = xcm[r + c * nQ] * dv;
                }
        var outr = new Complex[nQ * k];
        DressedSpMMSeconds[0] += sw.Elapsed.TotalSeconds; sw.Restart();
        // var rowPtr = table.rowPtr; var col = table.colIdx; var vr = table.vRe; var vi = table.vIm; // 261109Cl 変更前
        var (rowPtr, col, vr, vi) = source ? (table.fRowPtr, table.fColIdx, table.fRe, table.fIm) : (table.rowPtr, table.colIdx, table.vRe, table.vIm); // 261109Cl
        // 行ごとに、列の塊 (AVX 4 本 = 複素 8 個) の積算先をレジスタに置いたまま行の非零を全部足す (積算先のメモリの読み書きを非零ごとから塊ごとへ減らす)。
        // Σ a·x と Σ b·swap(x) を別々に積み、最後に AddSubtract で (Σ a·re − Σ b·im, Σ a·im + Σ b·re) にする
        unsafe void Row(int r)
        {
            int pw = qArr[r], p0 = rowPtr[pw], p1 = rowPtr[pw + 1], k2 = 2 * k;
            double dw = D[pw];
            fixed (Complex* xp0 = xr, yp0 = outr)
            fixed (int* colp = col, mapp = map)
            fixed (double* vrp = vr, vip = vi)
            {
                double* X = (double*)xp0, Y = (double*)yp0 + (long)r * k2;
                int c = 0;
                if (Fma.IsSupported && Avx.IsSupported)
                {
                    var vdw = Vector256.Create(dw);
                    for (; c + 16 <= k2; c += 16)
                    {
                        Vector256<double> a0 = default, a1 = default, a2 = default, a3 = default, b0 = default, b1 = default, b2 = default, b3 = default;
                        for (int p = p0; p < p1; p++)
                        {
                            int qv = mapp[colp[p]];
                            if (qv < 0) continue;
                            var va = Vector256.Create(vrp[p]); var vb = Vector256.Create(vip[p]);
                            double* xs = X + (long)qv * k2 + c;
                            var x0 = Avx.LoadVector256(xs); var x1 = Avx.LoadVector256(xs + 4); var x2 = Avx.LoadVector256(xs + 8); var x3 = Avx.LoadVector256(xs + 12);
                            a0 = Fma.MultiplyAdd(va, x0, a0); b0 = Fma.MultiplyAdd(vb, Avx.Permute(x0, 0b0101), b0);
                            a1 = Fma.MultiplyAdd(va, x1, a1); b1 = Fma.MultiplyAdd(vb, Avx.Permute(x1, 0b0101), b1);
                            a2 = Fma.MultiplyAdd(va, x2, a2); b2 = Fma.MultiplyAdd(vb, Avx.Permute(x2, 0b0101), b2);
                            a3 = Fma.MultiplyAdd(va, x3, a3); b3 = Fma.MultiplyAdd(vb, Avx.Permute(x3, 0b0101), b3);
                        }
                        Avx.Store(Y + c, Avx.Multiply(vdw, Avx.AddSubtract(a0, b0)));
                        Avx.Store(Y + c + 4, Avx.Multiply(vdw, Avx.AddSubtract(a1, b1)));
                        Avx.Store(Y + c + 8, Avx.Multiply(vdw, Avx.AddSubtract(a2, b2)));
                        Avx.Store(Y + c + 12, Avx.Multiply(vdw, Avx.AddSubtract(a3, b3)));
                    }
                    for (; c + 4 <= k2; c += 4)
                    {
                        Vector256<double> a0 = default, b0 = default;
                        for (int p = p0; p < p1; p++)
                        {
                            int qv = mapp[colp[p]];
                            if (qv < 0) continue;
                            var x0 = Avx.LoadVector256(X + (long)qv * k2 + c);
                            a0 = Fma.MultiplyAdd(Vector256.Create(vrp[p]), x0, a0); b0 = Fma.MultiplyAdd(Vector256.Create(vip[p]), Avx.Permute(x0, 0b0101), b0);
                        }
                        Avx.Store(Y + c, Avx.Multiply(vdw, Avx.AddSubtract(a0, b0)));
                    }
                }
                for (; c < k2; c += 2)
                {
                    double sre = 0, sim = 0;
                    for (int p = p0; p < p1; p++)
                    {
                        int qv = mapp[colp[p]];
                        if (qv < 0) continue;
                        double a = vrp[p], b = vip[p], xre = X[(long)qv * k2 + c], xim = X[(long)qv * k2 + c + 1];
                        sre += a * xre - b * xim; sim += a * xim + b * xre;
                    }
                    Y[c] = dw * sre; Y[c + 1] = dw * sim;
                }
            }
        }
        // 261109Cl 追加: 非零ごとに X の行を頭から流して読む形 (numba と同じ順。ハードウェアの先読みに合う)。積算先の行 (複素 k 個) は L1 に載る
        unsafe void RowStream(int r)
        {
            int pw = qArr[r], p0 = rowPtr[pw], p1 = rowPtr[pw + 1], k2 = 2 * k;
            double dw = D[pw];
            fixed (Complex* xp0 = xr, yp0 = outr)
            fixed (int* colp = col, mapp = map)
            fixed (double* vrp = vr, vip = vi)
            {
                double* X = (double*)xp0, Y = (double*)yp0 + (long)r * k2;
                for (int p = p0; p < p1; p++)
                {
                    int qv = mapp[colp[p]];
                    if (qv < 0) continue;
                    double a = vrp[p], b = vip[p];
                    double* xs = X + (long)qv * k2;
                    int c = 0;
                    if (Fma.IsSupported && Avx.IsSupported)
                    {
                        var va = Vector256.Create(a); var vb = Vector256.Create(b);
                        for (; c + 16 <= k2; c += 16)
                        {
                            var x0 = Avx.LoadVector256(xs + c); var x1 = Avx.LoadVector256(xs + c + 4); var x2 = Avx.LoadVector256(xs + c + 8); var x3 = Avx.LoadVector256(xs + c + 12);
                            Avx.Store(Y + c, Avx.AddSubtract(Fma.MultiplyAdd(va, x0, Avx.LoadVector256(Y + c)), Avx.Multiply(vb, Avx.Permute(x0, 0b0101))));
                            Avx.Store(Y + c + 4, Avx.AddSubtract(Fma.MultiplyAdd(va, x1, Avx.LoadVector256(Y + c + 4)), Avx.Multiply(vb, Avx.Permute(x1, 0b0101))));
                            Avx.Store(Y + c + 8, Avx.AddSubtract(Fma.MultiplyAdd(va, x2, Avx.LoadVector256(Y + c + 8)), Avx.Multiply(vb, Avx.Permute(x2, 0b0101))));
                            Avx.Store(Y + c + 12, Avx.AddSubtract(Fma.MultiplyAdd(va, x3, Avx.LoadVector256(Y + c + 12)), Avx.Multiply(vb, Avx.Permute(x3, 0b0101))));
                        }
                        for (; c + 4 <= k2; c += 4)
                        {
                            var x0 = Avx.LoadVector256(xs + c);
                            Avx.Store(Y + c, Avx.AddSubtract(Fma.MultiplyAdd(va, x0, Avx.LoadVector256(Y + c)), Avx.Multiply(vb, Avx.Permute(x0, 0b0101))));
                        }
                    }
                    for (; c < k2; c += 2)
                    {
                        double xre = xs[c], xim = xs[c + 1];
                        Y[c] += a * xre - b * xim; Y[c + 1] += a * xim + b * xre;
                    }
                }
                for (int c = 0; c < k2; c++) Y[c] *= dw;
            }
        }
        Action<int> kernel = DressedSpMMKernel == 1 ? RowStream : Row;
        if (parallel) Parallel.For(0, nQ, kernel);
        else for (int r = 0; r < nQ; r++) kernel(r);
        DressedSpMMSeconds[1] += sw.Elapsed.TotalSeconds; sw.Restart();
        var res = new CDense(nQ, k);
        var rv = res.Values;
        for (int c0 = 0; c0 < k; c0 += 32)
            for (int r0 = 0; r0 < nQ; r0 += 32)
                for (int c = c0, cEnd = Math.Min(c0 + 32, k); c < cEnd; c++)
                    for (int r = r0, rEnd = Math.Min(r0 + 32, nQ); r < rEnd; r++)
                        rv[r + c * nQ] = outr[r * k + c];
        DressedSpMMSeconds[2] += sw.Elapsed.TotalSeconds;
        return res;
    }

    /// <summary>261109Cl 追加 (研究用): 単精度 (複素 float) の疎な積。X と結合の値を float にし、積和も float (E·X は着衣の補正なので丸め 1e−6 程度は効かない見込み、照合で確かめる)。
    /// 核は RowStream と同じ順 (非零ごとに X の行を流す)、AVX 1 本に複素 4 個。</summary>
    // private static Matrix<Complex> SpMMTailSingle(DressedTailTable table, int[] qArr, int[] map, double[] D, Matrix<Complex> X, bool parallel) // 261109Cl 旧シグネチャ
    private static Matrix<Complex> SpMMTailSingle(DressedTailTable table, int[] qArr, int[] map, double[] D, Matrix<Complex> X, bool parallel, bool source = false) // 261109Cl source 追加
    {
        var sw = Stopwatch.StartNew();
        lock (table)
        {
            if (!source && table.vReF == null)
            {
                table.vReF = table.vRe.Select(v => (float)v).ToArray();
                table.vImF = table.vIm.Select(v => (float)v).ToArray();
            }
            if (source && table.fReF == null) // 261109Cl 追加
            {
                table.fReF = table.fRe.Select(v => (float)v).ToArray();
                table.fImF = table.fIm.Select(v => (float)v).ToArray();
            }
        }
        int nQ = X.RowCount, k = X.ColumnCount, k2 = 2 * k;
        var xcm = ((CDense)X).Values;
        var xr = new float[nQ * k2];
        for (int r0 = 0; r0 < nQ; r0 += 32)
            for (int c0 = 0; c0 < k; c0 += 32)
                for (int r = r0, rEnd = Math.Min(r0 + 32, nQ); r < rEnd; r++)
                {
                    double dv = D[qArr[r]];
                    for (int c = c0, cEnd = Math.Min(c0 + 32, k); c < cEnd; c++)
                    {
                        var x = xcm[r + c * nQ];
                        xr[r * k2 + 2 * c] = (float)(x.Real * dv); xr[r * k2 + 2 * c + 1] = (float)(x.Imaginary * dv);
                    }
                }
        var outr = new float[nQ * k2];
        DressedSpMMSeconds[0] += sw.Elapsed.TotalSeconds; sw.Restart();
        // var rowPtr = table.rowPtr; var col = table.colIdx; var vr = table.vReF; var vi = table.vImF; // 261109Cl 変更前
        var (rowPtr, col, vr, vi) = source ? (table.fRowPtr, table.fColIdx, table.fReF, table.fImF) : (table.rowPtr, table.colIdx, table.vReF, table.vImF); // 261109Cl
        unsafe void RowF(int r)
        {
            int pw = qArr[r], p0 = rowPtr[pw], p1 = rowPtr[pw + 1];
            fixed (float* xp0 = xr, yp0 = outr, vrp = vr, vip = vi)
            fixed (int* colp = col, mapp = map)
            {
                float* Y = yp0 + (long)r * k2;
                for (int p = p0; p < p1; p++)
                {
                    int qv = mapp[colp[p]];
                    if (qv < 0) continue;
                    float a = vrp[p], b = vip[p];
                    float* xs = xp0 + (long)qv * k2;
                    int c = 0;
                    if (Fma.IsSupported && Avx.IsSupported)
                    {
                        var va = Vector256.Create(a); var vb = Vector256.Create(b);
                        for (; c + 32 <= k2; c += 32)
                        {
                            var x0 = Avx.LoadVector256(xs + c); var x1 = Avx.LoadVector256(xs + c + 8); var x2 = Avx.LoadVector256(xs + c + 16); var x3 = Avx.LoadVector256(xs + c + 24);
                            Avx.Store(Y + c, Avx.AddSubtract(Fma.MultiplyAdd(va, x0, Avx.LoadVector256(Y + c)), Avx.Multiply(vb, Avx.Permute(x0, 0b10110001))));
                            Avx.Store(Y + c + 8, Avx.AddSubtract(Fma.MultiplyAdd(va, x1, Avx.LoadVector256(Y + c + 8)), Avx.Multiply(vb, Avx.Permute(x1, 0b10110001))));
                            Avx.Store(Y + c + 16, Avx.AddSubtract(Fma.MultiplyAdd(va, x2, Avx.LoadVector256(Y + c + 16)), Avx.Multiply(vb, Avx.Permute(x2, 0b10110001))));
                            Avx.Store(Y + c + 24, Avx.AddSubtract(Fma.MultiplyAdd(va, x3, Avx.LoadVector256(Y + c + 24)), Avx.Multiply(vb, Avx.Permute(x3, 0b10110001))));
                        }
                        for (; c + 8 <= k2; c += 8)
                        {
                            var x0 = Avx.LoadVector256(xs + c);
                            Avx.Store(Y + c, Avx.AddSubtract(Fma.MultiplyAdd(va, x0, Avx.LoadVector256(Y + c)), Avx.Multiply(vb, Avx.Permute(x0, 0b10110001))));
                        }
                    }
                    for (; c < k2; c += 2)
                    {
                        float xre = xs[c], xim = xs[c + 1];
                        Y[c] += a * xre - b * xim; Y[c + 1] += a * xim + b * xre;
                    }
                }
            }
        }
        if (parallel) Parallel.For(0, nQ, RowF);
        else for (int r = 0; r < nQ; r++) RowF(r);
        DressedSpMMSeconds[1] += sw.Elapsed.TotalSeconds; sw.Restart();
        var res = new CDense(nQ, k);
        var rv = res.Values;
        for (int c0 = 0; c0 < k; c0 += 32)
            for (int r0 = 0; r0 < nQ; r0 += 32)
                for (int c = c0, cEnd = Math.Min(c0 + 32, k); c < cEnd; c++)
                    for (int r = r0, rEnd = Math.Min(r0 + 32, nQ); r < rEnd; r++)
                    {
                        double dw = D[qArr[r]];
                        rv[r + c * nQ] = new Complex(outr[r * k2 + 2 * c] * dw, outr[r * k2 + 2 * c + 1] * dw);
                    }
        DressedSpMMSeconds[2] += sw.Elapsed.TotalSeconds;
        return res;
    }

    /// <summary>261109Cl 追加 (研究用): 疎な積を単精度にする (既定 false)</summary>
    public static bool DressedSpMMSingle { get; set; } = false;
    /// <summary>261109Cl 追加 (研究用): 疎な積の核の選択 (0 = 列の塊をレジスタに置く、1 = 非零ごとに行を流す)</summary>
    public static int DressedSpMMKernel { get; set; } = 1;
    /// <summary>261109Cl 追加 (研究用): 疎な積の内訳の積算 [s] (並べ替え・核・戻し)</summary>
    public static double[] DressedSpMMSeconds { get; } = new double[3];

    /// <summary>261109Cl 追加 (研究用、照合用の観測量): 点源の深さ積分 Σ_a w_a ∫₀^Z κe^{−κz} |Σ_g e^{2πi g·r_a} ψ_g(z)|² dz (tools の EbsdObs.point と同じ式)。</summary>
    /// <param name="sites">原子の分率座標と重み (占有率 × Z²)</param>
    public double PointSourceIntensityResearch(DressedTailSolution sol, (double x, double y, double z, double w)[] sites, double kappa = 1.0 / 15, double zMax = 75)
    {
        int nAll = sol.PoolIndices.Length, nm = sol.Gamma.Length;
        var F = new Complex[nm * nm];
        for (int j = 0; j < nm; j++)
            for (int kk = 0; kk < nm; kk++)
            {
                var dlt = sol.Gamma[j] - Complex.Conjugate(sol.Gamma[kk]);
                var aa = kappa - 2 * Math.PI * Complex.ImaginaryOne * dlt;
                F[j + kk * nm] = kappa * (1 - Complex.Exp(-aa * zMax)) / aa;
            }
        return PointSourceCore(sol, sites, F, nm, nAll);
    }

    private static double PointSourceCore(DressedTailSolution sol, (double x, double y, double z, double w)[] sites, Complex[] F, int nm, int nAll)
    {
        double total = 0;
        var phi = new Complex[nm];
        for (int si = 0; si < sites.Length; si++)
        {
            var s = sites[si];
            Array.Clear(phi);
            if (sol.AtomAmp != null) // 261109Cl: 共有の経路で求めた原子位置の振幅 (sites の順は atomSites と同じであること)
                for (int j = 0; j < nm; j++) phi[j] = sol.AtomAmp[si + j * sol.AtomCount];
            else
            for (int r = 0; r < nAll; r++)
            {
                var (h, k, l) = sol.Pool.Hkl[sol.PoolIndices[r]];
                var ph = Complex.FromPolarCoordinates(1, 2 * Math.PI * (h * s.x + k * s.y + l * s.z));
                for (int j = 0; j < nm; j++) phi[j] += ph * sol.Amp[r + j * nAll];
            }
            Complex acc = 0;
            for (int j = 0; j < nm; j++)
            {
                Complex t = 0;
                for (int kk = 0; kk < nm; kk++) t += F[j + kk * nm] * Complex.Conjugate(phi[kk]);
                acc += phi[j] * t;
            }
            total += s.w * acc.Real;
        }
        return total;
    }

    /// <summary>261109Cl 追加 (研究用): 源の裾どうしの項の求積に使う、原子ごとの F(0) の分 A_n (getU の後方半球の積分と同じ式) と Gauss の幅 σ_n = √(B/8π²) [nm]
    /// (等方の B_eq で近似)。並びは RunEBSDNew の atomArray と同じ (Crystal.Atoms → 各 Atom)。電圧は今の AccVoltage。表の F(0) と和が合うことを確かめる。</summary>
    public (double weight, double sigmaNm)[] DressedTailSourceAtomsResearch(DressedTailTable table)
    {
        if (table.diffF == null) throw new ArgumentException("The table has no backscatter source", nameof(table));
        var gammaRel = 1 + UniversalConstants.e0 * AccVoltage * 1E3 / UniversalConstants.m0 / UniversalConstants.c2;
        var list = new List<(double, double)>();
        foreach (var atoms in Crystal.Atoms)
        {
            double m = atoms.Dsf?.BisoEffective ?? 0;
            if (!(m > 0)) m = 0;
            int sub = _elasticIonModel == ElasticIonModel.Neutral ? 0 : atoms.SubNumberElectron;
            double f = m == 0 ? 0 : AtomStatic.ElectronScatteringPeng[atoms.AtomicNumber][sub].FactorImaginaryAnnular(AccVoltage, new Vector3DBase(0, 0, 0), m, Math.PI / 2, Math.PI, 30, 12);
            if (BackscatterSourceScaleForResearch != null) f *= BackscatterSourceScaleForResearch(atoms.AtomicNumber, AccVoltage);
            double a = f * atoms.Occ * gammaRel / Math.PI / Crystal.Volume, sigma = Math.Sqrt(m / (8 * Math.PI * Math.PI));
            foreach (var _ in atoms.Atom) list.Add((a, sigma));
        }
        var quadAtoms = list.ToArray();
        double sumA = quadAtoms.Sum(x => x.Item1), f0 = table.diffF[table.Pool.DiffIndex(0, 0, 0)].Real;
        if (Math.Abs(sumA - f0) > 1e-9 * Math.Max(Math.Abs(f0), 1e-300)) // 原子の分の和 = F(0) (規格化の取り違えを止める)
            throw new InvalidOperationException($"Tail quadrature weights {sumA} do not add up to F(0) = {f0}");
        return quadAtoms;
    }

    /// <summary>261109Cl 追加: 標準正規分布 N(0, 1) の Gauss–Hermite 求積の点と重み (Golub–Welsch、重みの和 = 1)。</summary>
    private static (double[] t, double[] w) GaussHermiteProbabilists(int n)
    {
        var J = Matrix<double>.Build.Dense(n, n);
        for (int k = 1; k < n; k++) J[k - 1, k] = J[k, k - 1] = Math.Sqrt(k);
        var evd = J.Evd(Symmetricity.Symmetric);
        var t = new double[n]; var w = new double[n];
        for (int i = 0; i < n; i++) { t[i] = evd.EigenValues[i].Real; w[i] = evd.EigenVectors[0, i] * evd.EigenVectors[0, i]; }
        double sw = w.Sum();
        for (int i = 0; i < n; i++) w[i] /= sw;
        return (t, w);
    }

    /// <summary>261109Cl 追加 (研究用、参照用): 候補の前方 (P &gt; 0) の全ビームを明示にした本体と同じ規約の行列 A (行 g・列 h、列優先、000 が先頭):
    /// A_gh = U(h − g)/P_g (g ≠ h)、A_gg = (Q_g + iU'₀)/P_g。固有分解 A C = C diag(γ) の C が物理振幅 (従来の getEigenMatrix と同じ)。</summary>
    public static Complex[] DressedTailExplicitMatrixForResearch(DressedTailTable table, Vector3DBase vecK0, Vector3DBase surface, out int[] poolIndices)
    {
        var pool = table.Pool;
        var (kX, kY, kZ) = vecK0.Tuple; var (sX, sY, sZ) = surface.Tuple;
        double k02 = vecK0.Length2;
        var idx = new List<int>();
        var Q = new List<double>(); var P = new List<double>();
        for (int i = 0; i < pool.Count; i++)
        {
            double vX = pool.GX[i] + kX, vY = pool.GY[i] + kY, vZ = pool.GZ[i] + kZ;
            double p = 2 * Dot3(sX, vX, sY, vY, sZ, vZ);
            if (i != 0 && p <= 0) continue;
            idx.Add(i); Q.Add(k02 - Dot3(vX, vX, vY, vY, vZ, vZ)); P.Add(p);
        }
        int n = idx.Count;
        var a = new Complex[(long)n * n];
        var u0i = Complex.ImaginaryOne * table.U000Imag;
        Parallel.For(0, n, c =>
        {
            var (h2, k2, l2) = pool.Hkl[idx[c]];
            for (int r = 0; r < n; r++)
            {
                var (h1, k1, l1) = pool.Hkl[idx[r]];
                a[r + (long)c * n] = (r == c ? Q[r] + u0i : table.diffUfull[pool.DiffIndex(h2 - h1, k2 - k1, l2 - l1)]) / P[r];
            }
        });
        poolIndices = [.. idx];
        return a;
    }

    /// <summary>261109Cl 追加 (研究用、照合用): 源の行列の密な形 F(行 r, 列 c) = F(g_r − g_c) (候補の番号の並び、列優先、差の箱の全部 = 疎化なし)。</summary>
    public static Complex[] DressedTailSourceForResearch(DressedTailTable table, int[] poolIndices)
    {
        if (table.diffF == null) throw new ArgumentException("The table has no backscatter source", nameof(table));
        var pool = table.Pool;
        int n = poolIndices.Length;
        var f = new Complex[n * n];
        Parallel.For(0, n, c =>
        {
            var (h2, k2, l2) = pool.Hkl[poolIndices[c]];
            for (int r = 0; r < n; r++)
            {
                var (h1, k1, l1) = pool.Hkl[poolIndices[r]];
                f[r + c * n] = table.diffF[pool.DiffIndex(h1 - h2, k1 - k2, l1 - l2)];
            }
        });
        return f;
    }

    /// <summary>261109Cl 追加 (研究用): ブロックの中に、中心で選んだビームが後ろ向き (P ≤ 0) になる方向がある。呼び出し側はブロックを分けて解き直す。</summary>
    public sealed class DressedTailBackwardBeamException : InvalidOperationException
    {
        public DressedTailBackwardBeamException() : base("A shared beam is not forward (P <= 0) in this direction; shrink the block") { }
    }

    /// <summary>261109Cl 追加 (研究用): MasterPattern の経路 (RunEBSDNew) を、着衣つきの裾の方式 (方式の核 + 方向のブロックの共有) で解く設定。
    /// <see cref="EbsdDressedTailResearch"/> に入れると、各エネルギーで「プール → 結合の表 → ブロックごとの解 → native の厚さ積分 (C = I、α = 1、phaseNG = 原子位置の振幅)」を回す。
    /// 局所源と後方 TDS 行列 (IncludeTDSBackground、F_W = W†FW を中心で 1 回、各方向 B = Z†F_W Z) に対応。非局所吸収・吸収フラックスの再注入・源のぼかしは非対応で、指定すると例外。単精度の疎な積は <see cref="DressedSpMMSingle"/> (静的) で選ぶ。</summary>
    public sealed class EbsdDressedTailOptions
    {
        /// <summary>裾の半径 R [nm⁻¹] (研究の既定 40)</summary>
        public double TailRadiusNm { get; set; } = 40;
        /// <summary>裾どうしの結合の閾値 η (通常用の候補 0.05、研究 0.02)</summary>
        public double Eta { get; set; } = 0.05;
        /// <summary>E1 で選ぶ本数 (000 を含む)</summary>
        public int N { get; set; } = 300;
        /// <summary>遅い反射の繰り上げの閾値 |Q/P| [nm⁻¹]</summary>
        public double Ps { get; set; } = 0.1;
        /// <summary>着衣の正則化の T [nm]</summary>
        public double RegDepthNm { get; set; } = 45;
        /// <summary>射影の E·X を着衣前の E·X₀ で近似する (疎な積 1 回。Galerkin の行列そのものを変える近似なので、既定は false = 2 回)</summary>
        public bool OneProduct { get; set; } = false;
        /// <summary>方向のブロックの一辺 (格子の画素の数。1 = 方向ごとに着衣、2 = 2×2 で共有し中心は画素の境目)</summary>
        public int BlockSize { get; set; } = 2;
        /// <summary>後方 TDS 行列 (IncludeTDSBackground) の源の表の疎化の閾値 (|F(d)| ≥ SourceEta·|F(0)|、0 = 差の箱の全部)</summary>
        public double SourceEta { get; set; } = 1e-3;
        /// <summary>源の裾どうしの項を Gauss–Hermite の求積で積むときの 1 軸の点の数 (0 = 源の表の疎な積。求積なら源の表は作らない)。
        /// Si 20 kV・G16・36 方向で、全部の F に対する差の最大: 1 点 0.96 %、3 点 0.18 %、5 点 0.022 % (η_F 0.1 の疎な積 0.037 %・1 方向 +15 s に対し +1 s)。既定 5</summary>
        public int SourceTailQuadrature { get; set; } = 5;
        /// <summary>同時に解くブロックの数の上限 (1 ブロックの作業領域は裾 1 万本・明示 300 本で約 0.5 GB)</summary>
        public int MaxParallelBlocks { get; set; } = 4;
        /// <summary>走行の記録 (エネルギーごとに 1 行追加)</summary>
        public List<EbsdDressedTailRecord> Records { get; } = [];
    }

    /// <summary>261109Cl 追加 (研究用): <see cref="EbsdDressedTailOptions"/> の 1 エネルギー分の記録。</summary>
    public sealed class EbsdDressedTailRecord
    {
        public double KV;
        public int Blocks, Directions, SplitBlocks, MinExplicit = int.MaxValue, MaxExplicit, MaxTail;
        /// <summary>中心とメンバーの方向の角度の最大 [mrad]</summary>
        public double MaxMemberAngleMrad;
        public double SecondsPool, SecondsTable, SecondsSource, SecondsSolve, SecondsIntegrate;
        public long TableNonzeros, SourceNonzeros;
    }

    /// <summary>261109Cl 追加 (研究用の口。既定 null = 従来の経路): MasterPattern の経路 (RunEBSDNew) を着衣つきの裾の方式で解く。研究用の道具だけが設定する。</summary>
    public EbsdDressedTailOptions EbsdDressedTailResearch { get; set; } = null;

    /// <summary>261109Cl 追加 (研究用): RunEBSDNew の 1 エネルギー分を着衣つきの裾の方式で解き、代表方向ごとの厚さの強度を返す (従来の ebsdIntensity と同じ形)。
    /// 代表方向を元の半球の格子の BlockSize × BlockSize の升目でまとめ、升目の中心 (有効な画素の方向の平均を正規化) で着衣を作って共有する。
    /// 各方向の k0・表面は従来と同じ式 (GetCrystalFixedMasterEquivalentBeamAndSurface → getVecK0)。ブロックの中に後ろ向きのビームがあれば方向ごとに解き直す。</summary>
    private double[][] SolveEbsdDressedTailVoltage(EbsdDressedTailOptions opt, ref BeamSelectionPool pool, double kvac, double u0,
        Vector3DBase referenceSurface, Vector3DBase referenceBeamDirection, Vector3DBase referenceAxisU, Vector3DBase referenceAxisV,
        List<int> representatives, int[] symmetryMapping, int hemisphereLength, int gridSize,
        (double x, double y, double z)[] atomArray, double[] sigmaArray, double tdsCoeff, double[][] background, int threads, Action<int> completed)
    {
        var rec = new EbsdDressedTailRecord { KV = AccVoltage };
        var sw = Stopwatch.StartNew();
        var mat = ((BaseRotation ?? new Matrix3D()) * Crystal.MatrixInverseTransposed).Tuple;
        if (pool == null || pool.RadiusNm != opt.TailRadiusNm || pool.Mat != mat)
            pool = CreateBeamSelectionPool(opt.TailRadiusNm, buildE0: false); // E1 の尺度は U/γ (電圧に依らない) なので、エネルギーの間で共有する
        rec.SecondsPool = sw.Elapsed.TotalSeconds; sw.Restart();
        bool tds = background != null;
        var table = CreateDressedTailTable(pool, opt.Eta, !tds ? double.NaN : opt.SourceTailQuadrature > 0 ? double.PositiveInfinity : opt.SourceEta); // 吸収は電圧で変わるので、表はエネルギーごとに作り直す
        (double weight, double sigmaNm)[] quadAtoms = tds && opt.SourceTailQuadrature > 0 ? DressedTailSourceAtomsResearch(table) : null;
        rec.SecondsTable = sw.Elapsed.TotalSeconds - table.SecondsSource; rec.SecondsSource = table.SecondsSource;
        rec.TableNonzeros = table.Nonzeros; rec.SourceNonzeros = table.SourceNonzeros;

        int b = Math.Max(1, opt.BlockSize), cellsPerRow = (gridSize + b - 1) / b;
        var cells = new Dictionary<long, List<int>>();
        foreach (var i in representatives)
        {
            int hemi = i / hemisphereLength, local = i - hemi * hemisphereLength, w = local % gridSize, h = local / gridSize;
            long key = ((long)hemi * cellsPerRow + h / b) * cellsPerRow + w / b;
            if (!cells.TryGetValue(key, out var list)) cells.Add(key, list = []);
            list.Add(i);
        }
        var cellList = cells.OrderBy(c => c.Key).ToArray();
        rec.Blocks = cellList.Length; rec.Directions = representatives.Count;

        (Vector3DBase k0, Vector3DBase n) Dir(Vector3DBase d)
        {
            var (eb, es) = GetCrystalFixedMasterEquivalentBeamAndSurface(d, referenceSurface, referenceBeamDirection, referenceAxisU, referenceAxisV);
            return (getVecK0(kvac, u0, eb, es), es);
        }
        var result = new double[BeamDirections.Length][];
        int nAtoms = atomArray.Length, split = 0;
        long solveTicks = 0, integrateTicks = 0;
        var gate = new object();
        Parallel.For(0, cellList.Length, new ParallelOptions { MaxDegreeOfParallelism = Math.Max(1, Math.Min(threads, opt.MaxParallelBlocks)) }, (c, state) =>
        {
            if (bwEBSDNew.CancellationPending) { state.Stop(); return; }
            var (key, members) = cellList[c];
            int hemi = (int)(key / ((long)cellsPerRow * cellsPerRow)), ch = (int)(key / cellsPerRow % cellsPerRow), cw = (int)(key % cellsPerRow);
            // 升目の中心: 有効な画素 (六方の無効セルを除く) の方向の平均。画素が 1 つなら、その方向そのもの (従来の経路とビット単位で同じ k0)
            double sx = 0, sy = 0, sz = 0; int valid = 0, last = -1;
            for (int h = ch * b; h < Math.Min(ch * b + b, gridSize); h++)
                for (int w = cw * b; w < Math.Min(cw * b + b, gridSize); w++)
                {
                    int p = hemi * hemisphereLength + h * gridSize + w;
                    if (symmetryMapping[p] < 0) continue;
                    var d = BeamDirections[p]; sx += d.X; sy += d.Y; sz += d.Z; valid++; last = p;
                }
            var center = valid == 1 ? BeamDirections[last] : Vector3DBase.Normarize(new Vector3DBase(sx, sy, sz));
            var (ck0, cn) = Dir(center);
            var dirs = members.Select(i => Dir(BeamDirections[i])).ToArray();
            double maxAngle = 0;
            foreach (var i in members)
                maxAngle = Math.Max(maxAngle, Math.Acos(Math.Clamp(Vector3DBase.Normarize(BeamDirections[i]) * center, -1, 1)) * 1e3);
            var swc = Stopwatch.StartNew();
            DressedTailSolution[] sols;
            try
            {
                sols = SolveDressedTailBlockResearch(table, ck0, cn, dirs, opt.N, opt.Ps, opt.RegDepthNm, opt.OneProduct, false, null, null, atomArray, tds, quadAtoms, opt.SourceTailQuadrature);
            }
            catch (DressedTailBackwardBeamException) when (members.Count > 1)
            {
                // P のクリップやビームの削除はせず、方向ごとに (自分を中心に) 解き直す
                Interlocked.Increment(ref split);
                sols = new DressedTailSolution[members.Count];
                for (int m = 0; m < members.Count; m++)
                    sols[m] = SolveDressedTailBlockResearch(table, dirs[m].k0, dirs[m].n, [dirs[m]], opt.N, opt.Ps, opt.RegDepthNm, opt.OneProduct, false, null, null, atomArray, tds, quadAtoms, opt.SourceTailQuadrature)[0];
            }
            Interlocked.Add(ref solveTicks, swc.Elapsed.Ticks);
            swc.Restart();
            int minS = int.MaxValue, maxS = 0, maxQ = 0;
            for (int m = 0; m < members.Count; m++)
            {
                var s = sols[m];
                int nm = s.Gamma.Length;
                // native の B = P·C·diag(α) に C = I・α = 1・P = Φ (原子 × モード、Φ が初期係数と固有ベクトルを含む) を渡す
                var ident = new Complex[nm * nm];
                for (int j = 0; j < nm; j++) ident[j * nm + j] = Complex.One;
                var ones = new Complex[nm];
                Array.Fill(ones, Complex.One);
                if (tds)
                {
                    // 従来の経路と同じ: 非局所源があれば局所値は使わない (原子 0 個として省き、計算済みの印の 0 配列を残す)。B は CreateMasterPatternMuBack と同じ並び (列優先、行 = 第 1 添字) で渡す
                    if (EigenEnabled)
                    {
                        var r2 = NativeWrapper.EBSDSolverWithTDS(nm, s.Gamma, ident, ones, Array.Empty<Complex>(), Array.Empty<double>(), s.SourceMatrix, tdsCoeff, Thicknesses);
                        result[members[m]] = r2.intensity; background[members[m]] = r2.tdsIntensity;
                    }
                    else
                    {
                        result[members[m]] = new double[Thicknesses.Length];
                        background[members[m]] = ComputeTDSMatrixBackground(nm, s.Gamma, ident, ones, s.SourceMatrix, tdsCoeff, Thicknesses);
                    }
                }
                else
                    result[members[m]] = EigenEnabled
                        ? NativeWrapper.EBSDSolver(nm, s.Gamma, ident, ones, s.AtomAmp, sigmaArray, Thicknesses)
                        : EBSDSolverManaged(nm, nAtoms, s.Gamma, ident, ones, s.AtomAmp, sigmaArray, Thicknesses);
                minS = Math.Min(minS, s.NExplicit); maxS = Math.Max(maxS, s.NExplicit); maxQ = Math.Max(maxQ, s.PoolIndices.Length - s.NExplicit);
                completed(members[m]);
            }
            Interlocked.Add(ref integrateTicks, swc.Elapsed.Ticks);
            lock (gate)
            {
                rec.MinExplicit = Math.Min(rec.MinExplicit, minS); rec.MaxExplicit = Math.Max(rec.MaxExplicit, maxS); rec.MaxTail = Math.Max(rec.MaxTail, maxQ);
                rec.MaxMemberAngleMrad = Math.Max(rec.MaxMemberAngleMrad, maxAngle);
            }
        });
        rec.SplitBlocks = split;
        rec.SecondsSolve = TimeSpan.FromTicks(solveTicks).TotalSeconds; rec.SecondsIntegrate = TimeSpan.FromTicks(integrateTicks).TotalSeconds;
        lock (opt.Records) opt.Records.Add(rec);
        return result;
    }
}
