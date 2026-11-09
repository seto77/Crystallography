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
        /// <summary>結合の表の非零の数</summary>
        public long Nonzeros => rowPtr[^1];
        public double SecondsUfull, SecondsCsr;
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
        /// <summary>物理振幅 M (行 = PoolIndices の順、列 = モード)、列優先 (行数 = PoolIndices.Length)</summary>
        public Complex[] Amp;
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
    public DressedTailTable CreateDressedTailTable(BeamSelectionPool pool, double eta)
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
        int n, double ps = 0.1, double regDepthNm = 45, bool oneProduct = false, bool parallel = false, DressedTailTiming timing = null, double[] blockSeconds = null)
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
                    throw new InvalidOperationException("A shared beam is not forward (P <= 0) in this direction; shrink the block");
            var K = Cm + iu0 * G0 - 2 * (kx * Gx + ky * Gy + kz * Gz);
            var M = 2 * nk * G0 + 2 * (nx * Gx + ny * Gy + nz * Gz);
            var ch = M.Cholesky();
            var evd = ch.Solve(K).Evd(Symmetricity.Asymmetric);
            var Cv = evd.EigenVectors;
            var e0 = MathNet.Numerics.LinearAlgebra.Vector<Complex>.Build.Dense(nS); e0[0] = 2 * nk * Dc[0]; // W†P e₀ = P_000 · D_c,000 e₀ (P_000 = 2 n·k0)
            var a = Cv.Solve(ch.Solve(e0));
            var Z = Cv * Matrix<Complex>.Build.DiagonalOfDiagonalVector(a);
            result[t] = MakeDressedSolution(pool, sArr, qArr, X, Z, Dc, 1.0, evd.EigenValues.ToArray(), cr.NSelected);
            if (blockSeconds != null) blockSeconds[1] += sw.Elapsed.TotalSeconds;
        }
        return result;
    }

    /// <summary>261109Cl 追加 (研究用): 裾どうしの疎な結合と密な行列の積 E·X (流束座標: E_wv = U(g_v − g_w) D_w D_v、w, v ∈ Q、w ≠ v)。
    /// X を行優先へ並べ替え、行ごとに CSR の非零 × X の行を AVX の複素の積和で足す。S の列 (map &lt; 0) は飛ばす。</summary>
    private static Matrix<Complex> SpMMTail(DressedTailTable table, int[] qArr, int[] map, double[] D, Matrix<Complex> X, bool parallel)
    {
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
        var rowPtr = table.rowPtr; var col = table.colIdx; var vr = table.vRe; var vi = table.vIm;
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
        foreach (var s in sites)
        {
            Array.Clear(phi);
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
}
