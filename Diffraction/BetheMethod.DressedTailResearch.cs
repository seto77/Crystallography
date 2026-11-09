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
        /// <summary>261109Cl 追加: F(−d) = conj F(d) が差の箱の全部で成り立つ (相対 1e−12)。成り立てば交差項を T + T† (T = F'_SQ X) で 1 回の積にする</summary>
        public bool SourceIsHermitian { get; internal set; }
        internal DressedQuadCache quadCache; // 261109Cl 追加: 源の求積の幾何の位相の表 (速い形が初めて使うときに作る)
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
    // public DressedTailTable CreateDressedTailTable(BeamSelectionPool pool, double eta, double sourceEta = double.NaN) // 261109Cl 旧シグネチャ
    public DressedTailTable CreateDressedTailTable(BeamSelectionPool pool, double eta, double sourceEta = double.NaN,
        double tailRadiusNm = double.NaN) // 261109Cl tailRadiusNm 追加: 裾の半径 (プールの半径以下、NaN = プールの半径)。遅い反射の繰り上げと E1 はプールの半径の全部から
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
            // Owner = this, Pool = pool, KV = kV, Eta = eta, TailRadiusNm = pool.RadiusNm, diffUfull = full, MaxAbsU = maxAbs, // 261109Cl 変更前
            Owner = this, Pool = pool, KV = kV, Eta = eta, TailRadiusNm = double.IsNaN(tailRadiusNm) ? pool.RadiusNm : Math.Min(tailRadiusNm, pool.RadiusNm), diffUfull = full, MaxAbsU = maxAbs, // 261109Cl
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
        if (!DressedTailBuildCsr) // 261110Cl 追加 (研究用、codex2 第 27 回): stencil の疎な積だけ使うなら CSR を作らない (表のメモリと構築の時間を省く)
        {
            table.rowPtr = new int[n + 1]; table.colIdx = []; table.vRe = []; table.vIm = [];
            table.SecondsCsr = sw.Elapsed.TotalSeconds;
            goto afterCsr;
        }
        // 261109Cl 追加 (codex2 第 21 回): 疎な表は裾の半径の中だけ (行・列とも)。積で使うのは裾どうしの結合だけなので値は変わらず、表と積の走査が軽くなる
        double rtt = table.TailRadiusNm * table.TailRadiusNm * (1 + 1e-9);
        var inT = new bool[n];
        for (int i = 0; i < n; i++) inT[i] = Dot3(pool.GX[i], pool.GX[i], pool.GY[i], pool.GY[i], pool.GZ[i], pool.GZ[i]) <= rtt;
        var rows = new (int j, Complex u)[n][];
        Parallel.For(0, n, i =>
        {
            if (!inT[i]) { rows[i] = []; return; } // 261109Cl
            var (h0, k0, l0) = pool.Hkl[i];
            var list = new List<(int, Complex)>();
            foreach (var d in ds)
            {
                int j = pool.IndexOf((h0 + d.h, k0 + d.k, l0 + d.l));
                // if (j >= 0) list.Add((j, d.u)); // 要素 (i, j) = U(g_j − g_i) = U_d // 261109Cl 変更前
                if (j >= 0 && inT[j]) list.Add((j, d.u)); // 要素 (i, j) = U(g_j − g_i) = U_d (261109Cl: 列も裾の半径の中だけ)
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
        afterCsr: // 261110Cl 追加
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
            double herm = 0, f0abs = fFull[Box(0, 0, 0)].Magnitude; // 261109Cl 追加: エルミート性の確認
            for (int h = -dh; h <= dh; h++)
                for (int k = -dk; k <= dk; k++)
                    for (int l = -dl; l <= dl; l++)
                        herm = Math.Max(herm, (fFull[Box(h, k, l)] - Complex.Conjugate(fFull[Box(-h, -k, -l)])).Magnitude);
            table.SourceIsHermitian = herm <= 1e-12 * f0abs;
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
        public Complex[] Gb; // 261109Cl 追加: 着衣のモードの固有値 γ_j (modeBasis のとき)
        public int[] Map; // 261109Cl 追加: 候補の番号 → 裾の行 (blocksOnly のとき)
        public Matrix<Complex> Xt, Vb, R, EXt; // 261109Cl 追加: モードの基底 (modeBasis) のとき = X̃ (X = X̃V_b⁻¹)、V_b、R = (H_QS + E X₀)V_b (着衣の倍率の前)、E·X̃ (疎な積 2 回のときだけ)
    }

    /// <summary>261109Cl 追加 (SolveDressedTailResearch から切り出し): 選択・繰り上げ・密なブロック・着衣 (Sylvester の 1 段) と、射影に使う E·X。</summary>
    // private DressCoreResult DressCore(DressedTailTable table, Vector3DBase vecK0, Vector3DBase surface, int n, double ps, double regDepthNm, bool oneProduct, bool parallel, DressedTailTiming timing) // 261109Cl 旧シグネチャ
    private DressCoreResult DressCore(DressedTailTable table, Vector3DBase vecK0, Vector3DBase surface, int n, double ps, double regDepthNm, bool oneProduct, bool parallel, DressedTailTiming timing,
        (Vector3DBase k0, Vector3DBase n)[] promoteAlso = null, // 261109Cl promoteAlso 追加: これらの方向のどれかで遅い (|Q/P| < ps) 裾の反射も明示へ繰り上げる (ブロックの共有用)
        int[] explicitSet = null, // 261109Cl explicitSet 追加 (研究用、照合用): 明示の集合をこの候補の番号に固定する (000 が先頭。選択と繰り上げをしない)
        bool modeBasis = false, // 261109Cl modeBasis 追加: X = X̃V_b⁻¹ を作らず、モードの基底 (V_b, X̃) と R を返す (span([I; X]) = span([V_b; X̃]) なので Galerkin は同じ)
        bool blocksOnly = false) // 261109Cl blocksOnly 追加 (研究用、native の試作): 選択と密なブロックまでで返す
    {
        if (table.Owner != this) throw new ArgumentException("The table was built for another BetheMethod", nameof(table));
        if (table.KV != (AccVoltage > 0 ? AccVoltage : 20)) // 261109Cl 追加 (codex2 第18回): 表の電圧と今の電圧の混在を禁じる
            throw new ArgumentException($"The table was built for {table.KV} kV but the current voltage is {AccVoltage} kV", nameof(table));
        var pool = table.Pool;
        var sw = Stopwatch.StartNew();
        // (1) 選択と繰り上げ
        // var e1 = SelectIndicesResearch(pool, vecK0, surface, n, BeamPartnerSet.E1StrongSet, 170, out _, out _); // 261109Cl 変更前
        var e1 = explicitSet != null ? [] : SelectIndicesResearch(pool, vecK0, surface, n, BeamPartnerSet.E1StrongSet, 170, out _, out _); // 261109Cl
        var (kX, kY, kZ) = vecK0.Tuple;
        var (sX, sY, sZ) = surface.Tuple;
        double k02 = vecK0.Length2;
        int cnt = pool.Count;
        var Q = new double[cnt]; var P = new double[cnt];
        double rt2 = table.TailRadiusNm * table.TailRadiusNm * (1 + 1e-9);
        var inTail = new bool[cnt];
        double rcore2 = DressedTailCoreRadiusNm * DressedTailCoreRadiusNm * (1 + 1e-9); // 261109Cl 追加 (研究用)
        for (int i = 0; i < cnt; i++)
        {
            double vX = pool.GX[i] + kX, vY = pool.GY[i] + kY, vZ = pool.GZ[i] + kZ;
            Q[i] = k02 - Dot3(vX, vX, vY, vY, vZ, vZ);
            P[i] = 2 * Dot3(sX, vX, sY, vY, sZ, vZ);
            // inTail[i] = P[i] > 0 && Dot3(pool.GX[i], pool.GX[i], pool.GY[i], pool.GY[i], pool.GZ[i], pool.GZ[i]) <= rt2; // 261109Cl 変更前
            double g2i = Dot3(pool.GX[i], pool.GX[i], pool.GY[i], pool.GY[i], pool.GZ[i], pool.GZ[i]);
            inTail[i] = P[i] > 0 && (double.IsNaN(DressedTailShellLambdaCut) ? g2i <= rt2
                : g2i <= rcore2 || (g2i <= rt2 && Math.Abs(Q[i] / P[i]) < DressedTailShellLambdaCut)); // 261109Cl 追加 (研究用): 芯の半径の外は Ewald 球に近い裾だけ
        }
        bool SlowInAny(int i) // 261109Cl 追加: ブロックの方向のどれかで遅い反射か
        {
            foreach (var (mk, mn) in promoteAlso)
            {
                var (ax, ay, az) = mk.Tuple; var (bx, by, bz) = mn.Tuple;
                double vX = pool.GX[i] + ax, vY = pool.GY[i] + ay, vZ = pool.GZ[i] + az;
                double q = mk.Length2 - Dot3(vX, vX, vY, vY, vZ, vZ), p = 2 * Dot3(bx, vX, by, vY, bz, vZ);
                if (p > 0 && Math.Abs(q / p) < ps) return true;
            }
            return false;
        }
        var inS = new bool[cnt];
        var S = new List<int>(e1.Length + 256) { 0 };
        inS[0] = true;
        foreach (var i in e1) { S.Add(i); inS[i] = true; }
        if (explicitSet != null) // 261109Cl 追加: 固定の明示の集合
        {
            if (explicitSet.Length == 0 || explicitSet[0] != 0) throw new ArgumentException("explicitSet must start with 000", nameof(explicitSet));
            foreach (var i in explicitSet.Skip(1))
            {
                // if (!inTail[i]) throw new DressedTailBackwardBeamException(); // 261109Cl 変更前
                if (!(P[i] > 0)) throw new DressedTailBackwardBeamException(); // 261109Cl: 明示はプールの半径の全部から (裾の半径の外でもよい)
                if (!inS[i]) { S.Add(i); inS[i] = true; }
            }
        }
        if (explicitSet == null && promoteAlso != null && DressedBlockUnionE1) // 261109Cl 追加 (研究用): ブロックの各方向の E1 の和集合
            foreach (var (mk, mn) in promoteAlso)
                foreach (var i in SelectIndicesResearch(pool, mk, mn, n, BeamPartnerSet.E1StrongSet, 170, out _, out _))
                    // if (!inS[i] && inTail[i]) { S.Add(i); inS[i] = true; } // 261109Cl 変更前
                    if (!inS[i] && P[i] > 0) { S.Add(i); inS[i] = true; } // 261109Cl
        int nSelected = S.Count;
        for (int i = 1; i < cnt && explicitSet == null; i++) // 261109Cl: 固定の集合なら繰り上げない
            // if (inTail[i] && !inS[i] && Math.Abs(Q[i] / P[i]) < ps) { S.Add(i); inS[i] = true; } // 261109Cl 変更前
            // if (inTail[i] && !inS[i] && (Math.Abs(Q[i] / P[i]) < ps || (promoteAlso != null && SlowInAny(i)))) { S.Add(i); inS[i] = true; } // 261109Cl 変更前
            if (P[i] > 0 && !inS[i] && (Math.Abs(Q[i] / P[i]) < ps || (promoteAlso != null && SlowInAny(i)))) { S.Add(i); inS[i] = true; } // 261109Cl: 繰り上げはプールの半径の全部から (裾の半径を小さくしても高次の遅い反射を明示に残す)
        if (!double.IsNaN(DressedTailNeighborRadiusNm)) // 261109Cl 追加 (研究用): 裾 = 芯の半径の中の全部 + 明示のビームから DressedTailNeighborRadiusNm 以内 (表の裾の半径まで)。繰り上げた高次の明示の近傍を着衣に入れる
        {
            double rn2 = DressedTailNeighborRadiusNm * DressedTailNeighborRadiusNm * (1 + 1e-9);
            // 261109Cl: 近傍の中心 = 明示 (DressedTailNeighborPs が NaN) か、候補のうち |Q/P| < DressedTailNeighborPs の反射 (繰り上げより広い閾値で、繰り上げの前から近傍を裾に入れて方向の段差を和らげる)
            IEnumerable<int> centers = double.IsNaN(DressedTailNeighborPs) ? S
                : S.Concat(Enumerable.Range(1, cnt - 1).Where(i => P[i] > 0 && Math.Abs(Q[i] / P[i]) < DressedTailNeighborPs));
            double rSkip = Math.Max(0, DressedTailCoreRadiusNm - DressedTailNeighborRadiusNm); // codex2 第 24 回: |g_s| + r_n ≤ R_c の中心は近傍が芯に入るので調べない
            var outerS = centers.Distinct().Where(s => Dot3(pool.GX[s], pool.GX[s], pool.GY[s], pool.GY[s], pool.GZ[s], pool.GZ[s]) > rSkip * rSkip).ToArray();
            if (DressedTailNeighborUseCsr) // 261110Cl 追加 (codex2 第 24 回の比較の候補): 近傍 = 中心から強い結合 (|U(d)| ≥ η max、表の CSR の列) で一段つながる裾。距離の球の代わり
            {
                var nearC = new bool[cnt];
                foreach (var c in centers.Distinct()) { nearC[c] = true; for (int pp = table.rowPtr[c]; pp < table.rowPtr[c + 1]; pp++) nearC[table.colIdx[pp]] = true; }
                for (int i = 1; i < cnt; i++)
                    if (inTail[i] && Dot3(pool.GX[i], pool.GX[i], pool.GY[i], pool.GY[i], pool.GZ[i], pool.GZ[i]) > rcore2) inTail[i] = nearC[i];
            }
            else
            for (int i = 1; i < cnt; i++)
            {
                if (!inTail[i]) continue;
                double gx = pool.GX[i], gy = pool.GY[i], gz = pool.GZ[i];
                if (Dot3(gx, gx, gy, gy, gz, gz) <= rcore2) continue;
                bool near = false;
                foreach (var s in outerS)
                {
                    double dx = gx - pool.GX[s], dy = gy - pool.GY[s], dz = gz - pool.GZ[s];
                    if (Dot3(dx, dx, dy, dy, dz, dz) <= rn2) { near = true; break; }
                }
                inTail[i] = near;
            }
        }
        var Qi = new List<int>(cnt);
        for (int i = 1; i < cnt; i++)
            if (inTail[i] && !inS[i]) Qi.Add(i);
        if (DressedTailMortonOrder) // 261110Cl 追加 (研究用、速さ): 裾の行を (h, k, l) の Morton 順に並べ替える (疎な積で隣の行の近傍が近くに並び、X の行のキャッシュの再利用が増える見込み。値は丸めの範囲で同じ)
        {
            static ulong Spread(uint v) { ulong x = v & 0x1FFFFF; x = (x | x << 32) & 0x1F00000000FFFF; x = (x | x << 16) & 0x1F0000FF0000FF; x = (x | x << 8) & 0x100F00F00F00F00F; x = (x | x << 4) & 0x10C30C30C30C30C3; x = (x | x << 2) & 0x1249249249249249; return x; }
            ulong Key(int i) { var (h, k, l) = pool.Hkl[i]; return Spread((uint)(h + 1024)) | Spread((uint)(k + 1024)) << 1 | Spread((uint)(l + 1024)) << 2; }
            Qi.Sort((a, c) => Key(a).CompareTo(Key(c)));
        }
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
        if (blocksOnly) // 261109Cl 追加
            return new DressCoreResult { SArr = sArr, QArr = qArr, P = P, D = D, NSelected = nSelected, HSS = HSS, HQS = HQS, HSQ = HSQ, LamQ = lamQ, Map = map };

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
        // var evdB = (HSS + HSQ * X0).Evd(Symmetricity.Asymmetric); // 261109Cl 変更前
        var hbm = HSS + MulMaybeSingle(HSQ, X0); // 261109Cl: 着衣の積は単精度を選べる (DressedDenseSingle)
        if (timing != null) timing.Project = sw.Elapsed.TotalSeconds; // 261109Cl 追加: X₀ と H_SQ X₀ まで
        var evdB = hbm.Evd(Symmetricity.Asymmetric);
        if (timing != null) timing.Eig = sw.Elapsed.TotalSeconds - timing.Project; // 261109Cl 追加: H_b の固有分解
        var gb = evdB.EigenValues; var Vb = evdB.EigenVectors;
        if (timing != null) timing.Dress = sw.Elapsed.TotalSeconds;
        sw.Restart();
        // var EX0 = SpMMTail(table, qArr, map, D, X0, parallel); // 261109Cl 変更前
        bool[] rowKeep = null; // 261109Cl 追加 (研究用、速さ): |λ_q| が DressedSpMMLambdaCut 以上の裾の行は E·X₀ を積まない (着衣で 1/λ_q 倍になる寄与)
        if (double.IsFinite(DressedSpMMLambdaCut)) { rowKeep = new bool[nQ]; for (int r = 0; r < nQ; r++) rowKeep[r] = lamQ[r].Magnitude < DressedSpMMLambdaCut; }
        Matrix<Complex> rhs, EX0 = null; // 261109Cl: EX0 は打ち切りでないときだけ (着衣前の E·X₀)
        int nM = nS; // 261109Cl 追加: 着衣のモードの数 (打ち切りなら nS より少ない)
        if (modeBasis && DressedModeTau > 0) // 261109Cl 追加 (研究用、native の試作と同じ規則): モードの組ごとの打ち切り。T = H_QS V_k、X₀V_k = −diag(f)T、R = T + E(X₀V_k)
        {
            var aIn = Vb.Solve(MathNet.Numerics.LinearAlgebra.Vector<Complex>.Build.Dense(nS, i => i == 0 ? Complex.One : Complex.Zero));
            var par = Enumerable.Range(0, nS).ToArray();
            int Find(int x) { while (par[x] != x) x = par[x] = par[par[x]]; return x; }
            if (DressedModeDegTol > 0)
                for (int i = 0; i < nS; i++)
                    for (int j = i + 1; j < nS; j++)
                        if ((gb[i] - gb[j]).Magnitude < DressedModeDegTol) { int ri = Find(i), rj = Find(j); if (ri != rj) par[ri] = rj; }
            var gsum = new Complex[nS * nS]; var used = new bool[nS]; var vbv = ((CDense)Vb).Values;
            for (int j = 0; j < nS; j++)
            {
                int r = Find(j); used[r] = true; var aj = aIn[j];
                for (int i = 0; i < nS; i++) gsum[i + r * nS] += vbv[i + j * nS] * aj;
            }
            var gw = new double[nS]; double wmax = 0;
            for (int r = 0; r < nS; r++)
                if (used[r]) { double s = 0; for (int i = 0; i < nS; i++) { var z = gsum[i + r * nS]; s += z.Real * z.Real + z.Imaginary * z.Imaginary; } gw[r] = Math.Sqrt(s); wmax = Math.Max(wmax, gw[r]); }
            var keep = Enumerable.Range(0, nS).Where(j => gw[Find(j)] >= DressedModeTau * wmax).ToArray();
            nM = keep.Length;
            var vk = new Complex[nS * nM];
            for (int t = 0; t < nM; t++) Array.Copy(vbv, keep[t] * nS, vk, t * nS, nS);
            Vb = new CDense(nS, nM, vk);
            gb = MathNet.Numerics.LinearAlgebra.Vector<Complex>.Build.Dense(nM, t => gb[keep[t]]);
            var T = MulMaybeSingle(HQS, Vb);
            var x0v = new CDense(nQ, nM); var x0vv = x0v.Values; var tv = ((CDense)T).Values;
            for (int r = 0; r < nQ; r++)
            {
                var l = lamQ[r]; var f = Complex.Conjugate(l) / (l.Real * l.Real + l.Imaginary * l.Imaginary + d2);
                for (int c = 0; c < nM; c++) x0vv[r + c * nQ] = -f * tv[r + c * nQ];
            }
            var EX0V = SpMMTail(table, qArr, map, D, x0v, parallel, false, rowKeep);
            if (timing != null) timing.SpMM1 = sw.Elapsed.TotalSeconds;
            sw.Restart();
            rhs = T + EX0V;
        }
        else
        {
            EX0 = SpMMTail(table, qArr, map, D, X0, parallel, false, rowKeep);
            if (timing != null) timing.SpMM1 = sw.Elapsed.TotalSeconds;
            sw.Restart();
            // var rhs = (HQS + EX0) * Vb; // 261109Cl 変更前
            rhs = MulMaybeSingle(HQS + EX0, Vb); // 261109Cl
        }
        var Rkeep = modeBasis ? rhs.Clone() : null; // 261109Cl 追加: 倍率の前の R
        var xt = ((CDense)rhs).Values;
        // for (int c = 0; c < nS; c++) // 261109Cl 変更前
        for (int c = 0; c < nM; c++) // 261109Cl: 打ち切りならモードの数
        {
            var g = gb[c];
            for (int r = 0; r < nQ; r++)
            {
                var dd = lamQ[r] - g;
                xt[r + c * nQ] = -xt[r + c * nQ] * Complex.Conjugate(dd) / (dd.Real * dd.Real + dd.Imaginary * dd.Imaginary + d2);
            }
        }
        if (modeBasis) // 261109Cl 追加: V_b⁻¹ の逆行列と積を省く
        {
            if (timing != null) timing.Dress += sw.Elapsed.TotalSeconds;
            sw.Restart();
            var EXt = oneProduct ? null : SpMMTail(table, qArr, map, D, rhs, parallel);
            if (timing != null) timing.SpMM2 = sw.Elapsed.TotalSeconds;
            return new DressCoreResult { SArr = sArr, QArr = qArr, P = P, D = D, NSelected = nSelected, HSS = HSS, HQS = HQS, HSQ = HSQ, LamQ = lamQ, Xt = rhs, Vb = Vb, R = Rkeep, EXt = EXt, Gb = gb.ToArray() };
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
        int tailQuadratureOrder = 0,      // 261109Cl 追加: その Gauss–Hermite の 1 軸の点の数 (0 = 求積せず源の表の疎な積)
        bool promoteOverBlock = false,    // 261109Cl 追加: ブロックのどれかの方向で遅い裾の反射も明示へ繰り上げる
        int[] explicitSet = null)         // 261109Cl 追加 (研究用、照合用): 明示の集合を固定する (候補の番号、000 が先頭)
    {
        // var cr = DressCore(table, centerK0, centerSurface, n, ps, regDepthNm, oneProduct, parallel, timing); // 261109Cl 変更前
        var cr = DressCore(table, centerK0, centerSurface, n, ps, regDepthNm, oneProduct, parallel, timing, promoteOverBlock ? dirs : null, explicitSet);
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
            bool herm = table.SourceIsHermitian; // 261109Cl: エルミートなら F'_QS = F'_SQ† なので作らない
            var fss = new CDense(nS, nS); var fsq = new CDense(nS, nQ); var fqs = new CDense(herm ? 0 : nQ, nS);
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
                    if (!herm) fqsv[r + c * nQ] = Ff(qr, sc) * dd;   // (行 Q_r, 列 S_c)
                    fsqv[c + r * nS] = Ff(sc, qr) * dd;   // (行 S_c, 列 Q_r)
                }
            });
            // 裾どうしの項 X†F'_QQ X。F はほとんど減衰しない (R40 で非零 5 千万〜1 億、η_F 0.1 で 1 方向 +15 s) ので、既定は求積:
            //   F の包絡は DW の Gauss exp(−B|d|²/4) とほぼ同じ (Si 20 kV で比 0.98〜1.00、day5e) → F(g − h) ≈ Σ_n A_n e^{−2πi(g−h)·r_n} ∫G_σn(x) e^{−2πi(g−h)·x} d³x
            //   (σ_n² = B_n/8π²、A_n = 原子 n の F(0) の分)。よって X†F'_QQ X ≈ Σ_n A_n Σ_q w_q χ_n,q† χ_n,q (χ = 裾の場 Σ_g D_g X_g e^{2πi g·(r_n + σ_n t_q)}、t_q = Gauss–Hermite の積の点)
            Matrix<Complex> qq;
            if (tailQuadratureOrder == 52 || tailQuadratureOrder == 77) // 261109Cl 追加 (codex2 第 19 回): 球面 (Lebedev) × 半径の求積 (次数 7 = 52 点、次数 9 = 77 点)
            {
                if (atomSites == null || tailQuadratureAtoms == null || tailQuadratureAtoms.Length != atomSites.Length)
                    throw new ArgumentException("The tail quadrature needs atomSites and one (weight, sigma) per site", nameof(tailQuadratureAtoms));
                var (pts, wpt) = GaussianCubatureSpherical(tailQuadratureOrder);
                int np = pts.Length, na = atomSites.Length, nSites = na * np;
                var ph = new CDense(nSites, nQ); var phv = ph.Values;
                Parallel.For(0, nQ, new ParallelOptions { MaxDegreeOfParallelism = parallel ? -1 : 1 }, r =>
                {
                    int q = qArr[r]; var (h, k, l) = pool.Hkl[q];
                    double gx = pool.GX[q], gy = pool.GY[q], gz = pool.GZ[q];
                    for (int a = 0; a < na; a++)
                    {
                        var (ax, ay, az) = atomSites[a]; double s2 = 2 * Math.PI * tailQuadratureAtoms[a].sigmaNm;
                        var b0 = Complex.FromPolarCoordinates(Dc[q], 2 * Math.PI * (h * ax + k * ay + l * az));
                        for (int t = 0; t < np; t++)
                        {
                            var (tx, ty, tz) = pts[t];
                            phv[a * np + t + (long)r * nSites] = b0 * Complex.FromPolarCoordinates(1, s2 * (gx * tx + gy * ty + gz * tz));
                        }
                    }
                });
                var chi = ph * X;
                var chiW = chi.Clone();
                for (int a = 0; a < na; a++)
                    for (int t = 0; t < np; t++)
                    {
                        double w = tailQuadratureAtoms[a].weight * wpt[t];
                        for (int c = 0; c < nS; c++) chiW[a * np + t, c] *= w;
                    }
                qq = chi.ConjugateTransposeThisAndMultiply(chiW);
            }
            else if (tailQuadratureOrder > 0)
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
            if (herm) { var cross = fsq * X; FW = fss + cross + cross.ConjugateTranspose() + qq; } // 261109Cl 追加 (codex2 第 19 回): 交差項を 1 回の積で
            else FW = fss + fsq * X + X.ConjugateTransposeThisAndMultiply(fqs) + qq;
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

    /// <summary>261109Cl 追加 (研究用、codex2 第 20 回): 一般の基底 W = D_c [A_S; A_Q] (D_c = 中心の P^{−1/2}、列 nB 本) での Galerkin と後方 TDS の源。
    /// 各方向で K = A†V'A − G_|g|² + iU'₀G₀ − 2Σ k0_a G_a、M = 2(n·k0)G₀ + 2Σ n_a G_a (G_w = A_S† diag(D_S² w) A_S + A_Q† diag(D_Q² w) A_Q、V' = 流束座標の非対角の結合)、
    /// 初期値 M a(0) = W†P e₀ = conj(A_S[0, :])·2(n·k0)·D_c,000、源 B = Z†F_W Z (裾どうしは球面の求積)。EAQ = E·A_Q (中心の流束座標の疎な積)。返り値は方向ごとの (γ, B)。</summary>
    private (Complex[] gamma, Complex[] source)[] GalerkinBasisResearch(DressedTailTable table, DressCoreResult cr, Matrix<Complex> AS, Matrix<Complex> AQ, Matrix<Complex> EAQ,
        (Vector3DBase k0, Vector3DBase n)[] dirs, (double x, double y, double z)[] atomSites, (double weight, double sigmaNm)[] quadAtoms, int quadOrder)
    {
        if (!table.SourceIsHermitian) throw new NotSupportedException("The general-basis research path assumes a Hermitian source");
        var pool = table.Pool;
        var (sArr, qArr, Dc) = (cr.SArr, cr.QArr, cr.D);
        int nS = sArr.Length, nQ = qArr.Length, nB = AS.ColumnCount;
        var hssOff = cr.HSS.Clone();
        for (int i = 0; i < nS; i++) hssOff[i, i] = 0;
        var Koff = AS.ConjugateTransposeThisAndMultiply(hssOff * AS + cr.HSQ * AQ) + AQ.ConjugateTransposeThisAndMultiply(cr.HQS * AS + EAQ);
        Matrix<Complex> Gram(Func<int, double> w)
        {
            var aS = Matrix<Complex>.Build.DenseOfMatrix(AS); var aQ = Matrix<Complex>.Build.DenseOfMatrix(AQ);
            var vs = ((CDense)aS).Values; var vq = ((CDense)aQ).Values;
            for (int c = 0; c < nB; c++)
            {
                for (int r = 0; r < nS; r++) { int s = sArr[r]; vs[r + c * nS] *= Dc[s] * Dc[s] * w(s); }
                for (int r = 0; r < nQ; r++) { int q = qArr[r]; vq[r + c * nQ] *= Dc[q] * Dc[q] * w(q); }
            }
            return AS.ConjugateTransposeThisAndMultiply(aS) + AQ.ConjugateTransposeThisAndMultiply(aQ);
        }
        var G0 = Gram(_ => 1); var Gx = Gram(i => pool.GX[i]); var Gy = Gram(i => pool.GY[i]); var Gz = Gram(i => pool.GZ[i]);
        var Gg2 = Gram(i => Dot3(pool.GX[i], pool.GX[i], pool.GY[i], pool.GY[i], pool.GZ[i], pool.GZ[i]));
        var Cm = Koff - Gg2;
        var iu0 = Complex.ImaginaryOne * table.U000Imag;
        // 源 F_W = A_S†F'_SS A_S + (A_S†F'_SQ A_Q) + (…)† + 裾どうし (求積)
        Complex Ff(int ig, int ih)
        {
            var (h1, k1, l1) = pool.Hkl[ig]; var (h2, k2, l2) = pool.Hkl[ih];
            return table.diffF[pool.DiffIndex(h1 - h2, k1 - k2, l1 - l2)];
        }
        var fss = new CDense(nS, nS);
        for (int c = 0; c < nS; c++)
            for (int r = 0; r < nS; r++) fss[r, c] = Ff(sArr[r], sArr[c]) * (Dc[sArr[r]] * Dc[sArr[c]]);
        var fsq = new CDense(nS, nQ); var fsqv = fsq.Values;
        Parallel.For(0, nQ, r =>
        {
            int qr = qArr[r];
            for (int c = 0; c < nS; c++) fsqv[c + r * nS] = Ff(sArr[c], qr) * (Dc[sArr[c]] * Dc[qr]);
        });
        var cross = AS.ConjugateTransposeThisAndMultiply(fsq * AQ);
        var (pts, wpt) = GaussianCubatureSpherical(quadOrder);
        int np = pts.Length, na = atomSites.Length, nSites = na * np;
        var ph = new CDense(nSites, nQ); var phv = ph.Values;
        Parallel.For(0, nQ, r =>
        {
            int q = qArr[r]; var (h, k, l) = pool.Hkl[q];
            double gx = pool.GX[q], gy = pool.GY[q], gz = pool.GZ[q];
            for (int a = 0; a < na; a++)
            {
                var (ax, ay, az) = atomSites[a]; double s2 = 2 * Math.PI * quadAtoms[a].sigmaNm;
                var b0 = Complex.FromPolarCoordinates(Dc[q], 2 * Math.PI * (h * ax + k * ay + l * az));
                for (int t = 0; t < np; t++)
                {
                    var (tx, ty, tz) = pts[t];
                    phv[a * np + t + (long)r * nSites] = b0 * Complex.FromPolarCoordinates(1, s2 * (gx * tx + gy * ty + gz * tz));
                }
            }
        });
        var chi = ph * AQ;
        var chiW = Matrix<Complex>.Build.DenseOfMatrix(chi);
        for (int a = 0; a < na; a++)
            for (int t = 0; t < np; t++)
            {
                double w = quadAtoms[a].weight * wpt[t];
                for (int c = 0; c < nB; c++) chiW[a * np + t, c] *= w;
            }
        var FW = AS.ConjugateTransposeThisAndMultiply(fss * AS) + cross + cross.ConjugateTranspose() + chi.ConjugateTransposeThisAndMultiply(chiW);
        var result = new (Complex[] gamma, Complex[] source)[dirs.Length];
        for (int t = 0; t < dirs.Length; t++)
        {
            var (k0, nv) = dirs[t];
            var (kx, ky, kz) = k0.Tuple; var (nx, ny, nz) = nv.Tuple;
            double nk = Dot3(nx, kx, ny, ky, nz, kz);
            var K = Cm + iu0 * G0 - 2 * (kx * Gx + ky * Gy + kz * Gz);
            var M = 2 * nk * G0 + 2 * (nx * Gx + ny * Gy + nz * Gz);
            var ch = M.Cholesky();
            var evd = ch.Solve(K).Evd(Symmetricity.Asymmetric);
            var Cv = evd.EigenVectors;
            var rhs = MathNet.Numerics.LinearAlgebra.Vector<Complex>.Build.Dense(nB, j => Complex.Conjugate(AS[0, j]) * (2 * nk * Dc[0]));
            var a = Cv.Solve(ch.Solve(rhs));
            var Z = Cv * Matrix<Complex>.Build.DiagonalOfDiagonalVector(a);
            result[t] = (evd.EigenValues.ToArray(), Z.ConjugateTransposeThisAndMultiply(FW * Z).ToColumnMajorArray());
        }
        return result;
    }

    /// <summary>261109Cl 追加 (研究用、codex2 第 20 回): 基底を足した共有の可能性の検査の結果。Variants の鍵 = "sharedPlain"・"localPlain"・"sharedAug{r}"・"localAug{r}"。
    /// 値はメンバーごとの (γ, B) (local は各メンバーを自分の方向で解いたもの)。</summary>
    public sealed class DressedTailAugmentedResult
    {
        public int NExplicit, NTail;
        public double[] SingularValues;
        public Dictionary<string, (Complex[] gamma, Complex[] source)[]> Variants = [];
    }

    /// <summary>261109Cl 追加 (研究用、codex2 第 20 回の可能性の検査): 中心の着衣 X_c (ブロック全体の繰り上げ + E1 の和集合で決めた S) に、
    /// メンバーの実際の着衣 (同じ S) との差 Δ_m = A_m X_m B_m − X_c (A_m = D_Q,c⁻¹D_Q,m、B_m = D_S,m⁻¹D_S,c) を、span([I; X_c]) に直交化してから
    /// 特異値分解で r 本に圧縮した追加の基底 T を足して共有する。比べる組: 共有 [W_c] / 方向ごと [W_m B_m] (= 方向ごとの解と同じ部分空間)、
    /// 共有 [W_c, T] / 方向ごと [W_m B_m, T]。源は後方 TDS (求積 quadOrder = 52 / 77)。</summary>
    public DressedTailAugmentedResult SolveDressedTailAugmentedResearch(DressedTailTable table, Vector3DBase centerK0, Vector3DBase centerSurface, (Vector3DBase k0, Vector3DBase n)[] dirs,
        int n, double ps, double regDepthNm, (double x, double y, double z)[] atomSites, (double weight, double sigmaNm)[] quadAtoms, int quadOrder, int[] ranks)
    {
        bool keep = DressedBlockUnionE1;
        DressedBlockUnionE1 = true;
        var crc = DressCore(table, centerK0, centerSurface, n, ps, regDepthNm, false, false, null, dirs);
        DressedBlockUnionE1 = keep;
        var (sArr, qArr) = (crc.SArr, crc.QArr);
        int nS = sArr.Length, nQ = qArr.Length, nm = dirs.Length;
        var Xc = Matrix<Complex>.Build.DenseOfMatrix(crc.X);
        var map = new int[table.Pool.Count]; Array.Fill(map, -1);
        for (int r = 0; r < nQ; r++) map[qArr[r]] = r;
        var deltas = new Matrix<Complex>[nm];
        for (int m = 0; m < nm; m++)
        {
            var crm = DressCore(table, dirs[m].k0, dirs[m].n, n, ps, regDepthNm, false, false, null, null, sArr);
            if (!crm.QArr.SequenceEqual(qArr) || !crm.SArr.SequenceEqual(sArr)) throw new InvalidOperationException("The member's tail differs from the centre's");
            var xm = Matrix<Complex>.Build.DenseOfMatrix(crm.X); var v = ((CDense)xm).Values;
            for (int c = 0; c < nS; c++)
            {
                double bc = crc.D[sArr[c]] / crm.D[sArr[c]];
                for (int r = 0; r < nQ; r++) v[r + c * nQ] *= crm.D[qArr[r]] / crc.D[qArr[r]] * bc;
            }
            deltas[m] = xm - Xc;
        }
        // 候補 [0; Δ_m] を span([I; X_c]) (流束の計量 = 係数のユークリッド) に直交化: R = C − Y (Y†Y)⁻¹ Y†C、Y†C = X_c†Δ
        var Ginv = (Xc.ConjugateTransposeThisAndMultiply(Xc) + Matrix<Complex>.Build.DenseIdentity(nS)).Inverse();
        var Rs = new CDense(nS, nm * nS); var Rq = new CDense(nQ, nm * nS);
        for (int m = 0; m < nm; m++)
        {
            var coef = Ginv * Xc.ConjugateTransposeThisAndMultiply(deltas[m]); // (nS × nS)
            Rs.SetSubMatrix(0, m * nS, -coef);
            Rq.SetSubMatrix(0, m * nS, deltas[m] - Xc * coef);
        }
        var gram = Rs.ConjugateTransposeThisAndMultiply(Rs) + Rq.ConjugateTransposeThisAndMultiply(Rq);
        var eg = gram.Evd(Symmetricity.Hermitian);
        var order = Enumerable.Range(0, gram.RowCount).OrderByDescending(i => eg.EigenValues[i].Real).ToArray();
        var res = new DressedTailAugmentedResult { NExplicit = nS, NTail = nQ, SingularValues = order.Take(80).Select(i => Math.Sqrt(Math.Max(0, eg.EigenValues[i].Real))).ToArray() };
        int rMax = ranks.Length == 0 ? 0 : ranks.Max();
        CDense TS = null, TQ = null;
        if (rMax > 0)
        {
            var Vr = new CDense(gram.RowCount, rMax);
            for (int j = 0; j < rMax; j++)
            {
                double sv = Math.Sqrt(Math.Max(eg.EigenValues[order[j]].Real, 1e-300));
                for (int i = 0; i < gram.RowCount; i++) Vr[i, j] = eg.EigenVectors[i, order[j]] / sv;
            }
            TS = (CDense)(Rs * Vr); TQ = (CDense)(Rq * Vr); // 正規直交 (流束の計量)
        }
        var ETQ = rMax > 0 ? SpMMTail(table, qArr, map, crc.D, TQ, false) : null;
        var EDel = deltas.Select(d => SpMMTail(table, qArr, map, crc.D, d, false)).ToArray();
        var I = Matrix<Complex>.Build.DenseIdentity(nS);
        Matrix<Complex> H(Matrix<Complex> a, Matrix<Complex> b) => a.Append(b);
        res.Variants["sharedPlain"] = GalerkinBasisResearch(table, crc, I, Xc, crc.EX, dirs, atomSites, quadAtoms, quadOrder);
        res.Variants["localPlain"] = Enumerable.Range(0, nm).Select(m => GalerkinBasisResearch(table, crc, I, Xc + deltas[m], crc.EX + EDel[m], [dirs[m]], atomSites, quadAtoms, quadOrder)[0]).ToArray();
        foreach (var r in ranks.Where(r => r > 0))
        {
            var ts = TS.SubMatrix(0, nS, 0, r); var tq = TQ.SubMatrix(0, nQ, 0, r); var etq = ETQ.SubMatrix(0, nQ, 0, r);
            res.Variants[$"sharedAug{r}"] = GalerkinBasisResearch(table, crc, H(I, ts), H(Xc, tq), H(crc.EX, etq), dirs, atomSites, quadAtoms, quadOrder);
            res.Variants[$"localAug{r}"] = Enumerable.Range(0, nm).Select(m => GalerkinBasisResearch(table, crc, H(I, ts), H(Xc + deltas[m], tq), H(crc.EX + EDel[m], etq), [dirs[m]], atomSites, quadAtoms, quadOrder)[0]).ToArray();
        }
        return res;
    }

    /// <summary>261109Cl 追加: 3 次元の標準正規分布の求積の点と重み。order = 52 / 77 は球面 (Lebedev) × 半径、それ以外は Gauss–Hermite の積 (order³ 点)。</summary>
    private static ((double x, double y, double z)[] pts, double[] w) GaussianCubature3D(int order)
    {
        if (order == 52 || order == 77) return GaussianCubatureSpherical(order);
        var (t, w1) = GaussHermiteProbabilists(order);
        int nq = t.Length;
        var pts = new (double, double, double)[nq * nq * nq]; var w = new double[nq * nq * nq];
        for (int i = 0; i < nq; i++)
            for (int j = 0; j < nq; j++)
                for (int m = 0; m < nq; m++) { int p = (i * nq + j) * nq + m; pts[p] = (t[i], t[j], t[m]); w[p] = w1[i] * w1[j] * w1[m]; }
        return (pts, w);
    }

    /// <summary>261109Cl 追加 (研究用、codex2 第 20 回の案 3): 1 方向を自分の着衣で解く軽い形 (共有の G 群を作らない)。
    /// 流束座標 (自分の D) で K = H_SS + H_SQ X + X†(H_QS + Λ X + E X)、G = I + X†X、G⁻¹K の固有分解、初期値 G c = √P₀ e₀。
    /// 原子位置の振幅 Φ = (Ph_S D_S + Ph_Q D_Q X) Z と、後方 TDS の源 B = Z†F_W Z (F_W の裾どうしは求積) を返す。
    /// 中心 = 自分の方向の共有の解 (SolveDressedTailBlockResearch の 1 方向) と同じ値になる。</summary>
    public DressedTailSolution SolveDressedTailSingleResearch(DressedTailTable table, Vector3DBase vecK0, Vector3DBase surface, int n, double ps, double regDepthNm, bool oneProduct,
        (double x, double y, double z)[] atomSites, bool backscatterSource, (double weight, double sigmaNm)[] quadAtoms, int quadOrder, int[] explicitSet = null,
        DressedTailTiming timing = null, double[] seconds = null)
    {
        var cr = DressCore(table, vecK0, surface, n, ps, regDepthNm, oneProduct, false, timing, null, explicitSet);
        var pool = table.Pool;
        var (sArr, qArr, P, D, X, HSS, HQS, HSQ, EX, lamQ) = (cr.SArr, cr.QArr, cr.P, cr.D, cr.X, cr.HSS, cr.HQS, cr.HSQ, cr.EX, cr.LamQ);
        int nS = sArr.Length, nQ = qArr.Length;
        var sw = Stopwatch.StartNew();
        var hrq = (HQS + EX).ToColumnMajorArray();
        var xv = ((CDense)X).Values;
        for (int c = 0; c < nS; c++)
            for (int r = 0; r < nQ; r++) hrq[r + c * nQ] += lamQ[r] * xv[r + c * nQ];
        var K = HSS + HSQ * X + X.ConjugateTransposeThisAndMultiply(new CDense(nQ, nS, hrq));
        var G = X.ConjugateTransposeThisAndMultiply(X);
        for (int i = 0; i < nS; i++) G[i, i] += 1;
        var ch = G.Cholesky();
        var evd = ch.Solve(K).Evd(Symmetricity.Asymmetric);
        var Cv = evd.EigenVectors;
        var e0 = MathNet.Numerics.LinearAlgebra.Vector<Complex>.Build.Dense(nS); e0[0] = Math.Sqrt(P[0]);
        var a = Cv.Solve(ch.Solve(e0));
        var Z = Cv * Matrix<Complex>.Build.DiagonalOfDiagonalVector(a);
        if (seconds != null) seconds[0] += sw.Elapsed.TotalSeconds;
        sw.Restart();
        var sol = new DressedTailSolution { Pool = pool, NExplicit = nS, NSelected = cr.NSelected, Gamma = evd.EigenValues.ToArray() };
        var idxAll = new int[nS + nQ];
        Array.Copy(sArr, idxAll, nS); Array.Copy(qArr, 0, idxAll, nS, nQ);
        sol.PoolIndices = idxAll;
        if (atomSites != null)
        {
            int na = atomSites.Length;
            var phS = new CDense(na, nS); var phQ = new CDense(na, nQ);
            for (int at = 0; at < na; at++)
            {
                var (ax, ay, az) = atomSites[at];
                for (int i = 0; i < nS; i++) { var (h, k, l) = pool.Hkl[sArr[i]]; phS[at, i] = Complex.FromPolarCoordinates(D[sArr[i]], 2 * Math.PI * (h * ax + k * ay + l * az)); }
                for (int i = 0; i < nQ; i++) { var (h, k, l) = pool.Hkl[qArr[i]]; phQ[at, i] = Complex.FromPolarCoordinates(D[qArr[i]], 2 * Math.PI * (h * ax + k * ay + l * az)); }
            }
            sol.AtomAmp = ((phS + phQ * X) * Z).ToColumnMajorArray();
            sol.AtomCount = na;
        }
        if (backscatterSource)
        {
            if (table.diffF == null) throw new ArgumentException("The table has no backscatter source; create it with sourceEta", nameof(table));
            if (!table.SourceIsHermitian) throw new NotSupportedException("The single-direction path assumes a Hermitian source");
            if (atomSites == null || quadAtoms == null || quadAtoms.Length != atomSites.Length || quadOrder <= 0)
                throw new ArgumentException("The single-direction path needs atomSites and the tail quadrature", nameof(quadAtoms));
            Complex Ff(int ig, int ih)
            {
                var (h1, k1, l1) = pool.Hkl[ig]; var (h2, k2, l2) = pool.Hkl[ih];
                return table.diffF[pool.DiffIndex(h1 - h2, k1 - k2, l1 - l2)];
            }
            var fss = new CDense(nS, nS);
            for (int c = 0; c < nS; c++)
                for (int r = 0; r < nS; r++) fss[r, c] = Ff(sArr[r], sArr[c]) * (D[sArr[r]] * D[sArr[c]]);
            var fsq = new CDense(nS, nQ); var fsqv = fsq.Values;
            for (int r = 0; r < nQ; r++)
            {
                int qr = qArr[r];
                for (int c = 0; c < nS; c++) fsqv[c + r * nS] = Ff(sArr[c], qr) * (D[sArr[c]] * D[qr]);
            }
            var cross = fsq * X;
            var (pts, wpt) = GaussianCubature3D(quadOrder);
            int np = pts.Length, na = atomSites.Length, nSites = na * np;
            var ph = new CDense(nSites, nQ); var phv = ph.Values;
            for (int r = 0; r < nQ; r++)
            {
                int q = qArr[r]; var (h, k, l) = pool.Hkl[q];
                double gx = pool.GX[q], gy = pool.GY[q], gz = pool.GZ[q];
                for (int at = 0; at < na; at++)
                {
                    var (ax, ay, az) = atomSites[at]; double s2 = 2 * Math.PI * quadAtoms[at].sigmaNm;
                    var b0 = Complex.FromPolarCoordinates(D[q], 2 * Math.PI * (h * ax + k * ay + l * az));
                    for (int t = 0; t < np; t++)
                    {
                        var (tx, ty, tz) = pts[t];
                        phv[at * np + t + (long)r * nSites] = b0 * Complex.FromPolarCoordinates(1, s2 * (gx * tx + gy * ty + gz * tz));
                    }
                }
            }
            var chi = ph * X;
            var chiW = Matrix<Complex>.Build.DenseOfMatrix(chi);
            var cwv = ((CDense)chiW).Values;
            for (int c = 0; c < nS; c++)
                for (int at = 0; at < na; at++)
                    for (int t = 0; t < np; t++) cwv[at * np + t + c * nSites] *= quadAtoms[at].weight * wpt[t];
            var FW = fss + cross + cross.ConjugateTranspose() + chi.ConjugateTransposeThisAndMultiply(chiW);
            sol.SourceMatrix = Z.ConjugateTransposeThisAndMultiply(FW * Z).ToColumnMajorArray();
        }
        if (seconds != null) seconds[1] += sw.Elapsed.TotalSeconds;
        return sol;
    }

    /// <summary>261109Cl 追加 (研究用、速さ): 1 方向の軽い解の速い形。モードの基底 A = [V_b; X̃] (span = [I; X]) で Galerkin を組み、X = X̃V_b⁻¹ の積と逆行列を省く。
    /// K' = V_b†H_SS V_b + V_b†(H_SQ X̃) + X̃†(H_QS V_b + ΛX̃ + E X̃) (疎な積 1 回なら H_QS V_b + E X̃ ≈ R = (H_QS + E X₀)V_b、着衣で作った積をそのまま使う)、
    /// G' = V_b†V_b + X̃†X̃、初期値 G'c = A†(√P₀ e₀) = conj(V_b[0, :])·√P₀。後方 TDS の源の求積は位相を (原子) × (点) に因数分解して作る。
    /// 後方 TDS のとき原子位置の振幅は作らない (本体の経路は非局所源があれば局所値を使わない)。値は SolveDressedTailSingleResearch と丸めの範囲で同じ。
    /// seconds (長さ 7 以上) に段ごとの時間を足す: 0 着衣 (選択・ブロック・疎な積を含む)、1 Galerkin、2 原子の振幅、3 源の密なブロック、4 求積の位相、5 求積の積、6 B = Z†F_W Z。</summary>
    public DressedTailSolution SolveDressedTailFastResearch(DressedTailTable table, Vector3DBase vecK0, Vector3DBase surface, int n, double ps, double regDepthNm, bool oneProduct,
        (double x, double y, double z)[] atomSites, bool backscatterSource, (double weight, double sigmaNm)[] quadAtoms, int quadOrder, int[] explicitSet = null, double[] seconds = null)
    {
        var sw = Stopwatch.StartNew();
        var tmc = seconds != null && seconds.Length >= 13 ? new DressedTailTiming() : null; // 261109Cl: 着衣の内訳 (seconds[7..12])
        var cr = DressCore(table, vecK0, surface, n, ps, regDepthNm, oneProduct, false, tmc, null, explicitSet, modeBasis: true);
        if (tmc != null) { seconds[7] += tmc.Select; seconds[8] += tmc.Blocks; seconds[9] += tmc.Dress; seconds[10] += tmc.SpMM1; seconds[11] += tmc.SpMM2; if (seconds.Length >= 15) { seconds[13] += tmc.Project; seconds[14] += tmc.Eig; } }
        var pool = table.Pool;
        var (sArr, qArr, P, D, HSS, HQS, HSQ, lamQ, Xt, Vb) = (cr.SArr, cr.QArr, cr.P, cr.D, cr.HSS, cr.HQS, cr.HSQ, cr.LamQ, cr.Xt, cr.Vb);
        int nS = sArr.Length, nQ = qArr.Length, nM = Vb.ColumnCount; // 261109Cl nM 追加: 着衣のモードの数 (打ち切りなら nS より少ない)
        void Lap(int slot) { if (seconds != null) seconds[slot] += sw.Elapsed.TotalSeconds; sw.Restart(); }
        Lap(0);
        // (1) Galerkin (モードの基底)
        Matrix<Complex> K, G;
        if (oneProduct && !double.IsNaN(DressedEdeltaThreshold))
        {
            // 261109Cl 追加 (codex2 第 21 回): 疎な積 1 回では ΛX̃ + R = X̃Γ + E_δ、(E_δ)_qj = R_qj δ²/(|λ_q − γ_j|² + δ²) なので
            //   X̃†(R + ΛX̃) = (X̃†X̃)Γ + X̃†E_δ。X̃†X̃ は G にも要るので使い回し、E_δ は共鳴の近く (倍率 > 閾値) だけの疎な積にする
            var YY = Xt.ConjugateTransposeThisAndMultiply(Xt);
            var gb = cr.Gb;
            double delta = 1 / (Math.PI * regDepthNm), d2 = delta * delta, thr = DressedEdeltaThreshold;
            var Rv = ((CDense)cr.R).Values; var xv = ((CDense)Xt).Values;
            var YE = new Complex[nM * nM];
            long nnz = 0;
            for (int j = 0; j < nM; j++)
            {
                var g = gb[j];
                for (int q = 0; q < nQ; q++)
                {
                    var dd = lamQ[q] - g;
                    double f = d2 / (dd.Real * dd.Real + dd.Imaginary * dd.Imaginary + d2);
                    if (f <= thr) continue;
                    var e = Rv[q + j * nQ] * f;
                    nnz++;
                    for (int i = 0; i < nM; i++) YE[i + j * nM] += Complex.Conjugate(xv[q + i * nQ]) * e; // (X̃†E_δ)_ij
                }
            }
            Interlocked.Add(ref DressedEdeltaNonzeros, nnz);
            var YYG = YY.Clone();
            for (int j = 0; j < nM; j++) for (int i = 0; i < nM; i++) YYG[i, j] *= gb[j];
            K = Vb.ConjugateTransposeThisAndMultiply(HSS * Vb + HSQ * Xt) + YYG + new CDense(nM, nM, YE);
            G = Vb.ConjugateTransposeThisAndMultiply(Vb) + YY;
        }
        else
        {
            var right = oneProduct ? cr.R : HQS * Vb + cr.EXt; // H_QS V_b + E X̃
            var rv = ((CDense)right).Values; var xv = ((CDense)Xt).Values;
            for (int c = 0; c < nM; c++)
                for (int r = 0; r < nQ; r++) rv[r + c * nQ] += lamQ[r] * xv[r + c * nQ];
            K = Vb.ConjugateTransposeThisAndMultiply(HSS * Vb + HSQ * Xt) + Xt.ConjugateTransposeThisAndMultiply(right);
            G = Vb.ConjugateTransposeThisAndMultiply(Vb) + Xt.ConjugateTransposeThisAndMultiply(Xt);
        }
        if (seconds != null && seconds.Length >= 13) { seconds[12] -= sw.Elapsed.TotalSeconds; } // 261109Cl: 固有分解だけの時間 (下で足す)
        var ch = G.Cholesky();
        var evd = ch.Solve(K).Evd(Symmetricity.Asymmetric);
        if (seconds != null && seconds.Length >= 13) { seconds[12] += sw.Elapsed.TotalSeconds; }
        var Cv = evd.EigenVectors;
        double sp0 = Math.Sqrt(P[0]);
        var e0 = MathNet.Numerics.LinearAlgebra.Vector<Complex>.Build.Dense(nM, j => Complex.Conjugate(Vb[0, j]) * sp0);
        var a = Cv.Solve(ch.Solve(e0));
        var Z = Cv * Matrix<Complex>.Build.DiagonalOfDiagonalVector(a);
        Lap(1);
        var sol = new DressedTailSolution { Pool = pool, NExplicit = nS, NSelected = cr.NSelected, Gamma = evd.EigenValues.ToArray() };
        var idxAll = new int[nS + nQ];
        Array.Copy(sArr, idxAll, nS); Array.Copy(qArr, 0, idxAll, nS, nQ);
        sol.PoolIndices = idxAll;
        // (2) 原子位置の振幅 (局所源のときだけ): Φ = (Ph_S D_S V_b + Ph_Q D_Q X̃) Z
        if (atomSites != null && !backscatterSource)
        {
            int na = atomSites.Length;
            var phS = new CDense(na, nS); var phQ = new CDense(na, nQ);
            for (int at = 0; at < na; at++)
            {
                var (ax, ay, az) = atomSites[at];
                for (int i = 0; i < nS; i++) { var (h, k, l) = pool.Hkl[sArr[i]]; phS[at, i] = Complex.FromPolarCoordinates(D[sArr[i]], 2 * Math.PI * (h * ax + k * ay + l * az)); }
                for (int i = 0; i < nQ; i++) { var (h, k, l) = pool.Hkl[qArr[i]]; phQ[at, i] = Complex.FromPolarCoordinates(D[qArr[i]], 2 * Math.PI * (h * ax + k * ay + l * az)); }
            }
            sol.AtomAmp = ((phS * Vb + phQ * Xt) * Z).ToColumnMajorArray();
            sol.AtomCount = na;
        }
        Lap(2);
        // (3) 後方 TDS の源: F_W = V_b†F'_SS V_b + V_b†(F'_SQ X̃) + (…)† + X̃†F'_QQ X̃ (求積)、B = Z†F_W Z
        if (backscatterSource)
        {
            if (table.diffF == null) throw new ArgumentException("The table has no backscatter source; create it with sourceEta", nameof(table));
            if (!table.SourceIsHermitian) throw new NotSupportedException("The fast path assumes a Hermitian source");
            if (atomSites == null || quadAtoms == null || quadAtoms.Length != atomSites.Length || quadOrder <= 0)
                throw new ArgumentException("The fast path needs atomSites and the tail quadrature", nameof(quadAtoms));
            // 261109Cl 追加 (研究用): 源だけの外殻。S ∪ Q の外で半径 DressedSourceShellRadiusNm 以内のビーム A の振幅を、明示との結合の一次 X_A = −f_A(γ_j) (H_AS V_k)_j で作り、
            //   源の行列にだけ足す (Galerkin には入れない)。多数の弱い高次の振幅が原子の位置で揃って効く分 (源の交差項) を安く取り込む試し
            var qArrSrc = qArr; int nQs = nQ; var XtAll = Xt;
            var qcEarly = !double.IsNaN(DressedSourceShellRadiusNm) ? GetQuadCache(table, atomSites, quadAtoms, quadOrder) : null;
            if (qcEarly != null)
            {
                double ro2 = DressedSourceShellRadiusNm * DressedSourceShellRadiusNm * (1 + 1e-9);
                var used = new bool[pool.Count]; foreach (var i in sArr) used[i] = true; foreach (var i in qArr) used[i] = true;
                var (kX, kY, kZ) = vecK0.Tuple; double k02 = vecK0.Length2;
                var aList = new List<int>();
                for (int i = 1; i < pool.Count; i++)
                    if (!used[i] && P[i] > 0 && qcEarly.PoolCol[i] >= 0 && Dot3(pool.GX[i], pool.GX[i], pool.GY[i], pool.GY[i], pool.GZ[i], pool.GZ[i]) <= ro2) aList.Add(i);
                int nA = aList.Count;
                var has = new CDense(nA, nS); var hasv = has.Values;
                for (int c = 0; c < nS; c++)
                {
                    var (h2, k2, l2) = pool.Hkl[sArr[c]]; double dc = D[sArr[c]];
                    for (int r = 0; r < nA; r++) { var (h1, k1, l1) = pool.Hkl[aList[r]]; hasv[r + c * nA] = table.diffUfull[pool.DiffIndex(h2 - h1, k2 - k1, l2 - l1)] * (D[aList[r]] * dc); }
                }
                var TA = MulMaybeSingle(has, Vb); var tav = ((CDense)TA).Values;
                double dl = 1 / (Math.PI * regDepthNm), dl2 = dl * dl; var gbA = cr.Gb;
                for (int r = 0; r < nA; r++)
                {
                    int i = aList[r];
                    double vX = pool.GX[i] + kX, vY = pool.GY[i] + kY, vZ = pool.GZ[i] + kZ;
                    var lam = (k02 - Dot3(vX, vX, vY, vY, vZ, vZ) + Complex.ImaginaryOne * table.U000Imag) / P[i]; // DressCore の Lam と同じ
                    for (int j = 0; j < nM; j++) { var dd = lam - gbA[j]; tav[r + j * nA] = -tav[r + j * nA] * Complex.Conjugate(dd) / (dd.Real * dd.Real + dd.Imaginary * dd.Imaginary + dl2); }
                }
                qArrSrc = [.. qArr, .. aList]; nQs = nQ + nA;
                XtAll = Xt.Stack(TA);
                Interlocked.Add(ref DressedSourceShellCount, nA);
            }
            // var fss = new CDense(nS, nS); var fsq = new CDense(nS, nQ); var fssv = fss.Values; var fsqv = fsq.Values; // 261109Cl 変更前
            var fss = new CDense(nS, nS); var fsq = new CDense(nS, nQs); var fssv = fss.Values; var fsqv = fsq.Values;
            var dF = table.diffF;
            int dkn = 2 * pool.DiffBoundK + 1, dln = 2 * pool.DiffBoundL + 1;
            var sh = new (int h, int k, int l)[nS];
            for (int i = 0; i < nS; i++) sh[i] = pool.Hkl[sArr[i]];
            for (int c = 0; c < nS; c++)
            {
                var (h2, k2, l2) = sh[c]; double dc = D[sArr[c]];
                for (int r = 0; r < nS; r++) { var (h1, k1, l1) = sh[r]; fssv[r + c * nS] = dF[pool.DiffIndex(h1 - h2, k1 - k2, l1 - l2)] * (D[sArr[r]] * dc); }
            }
            // for (int r = 0; r < nQ; r++) // 261109Cl 変更前
            for (int r = 0; r < nQs; r++) // 261109Cl: 源だけの外殻も
            {
                // int qr = qArr[r]; var (h2, k2, l2) = pool.Hkl[qr]; double dq = D[qr]; // 261109Cl 変更前
                int qr = qArrSrc[r]; var (h2, k2, l2) = pool.Hkl[qr]; double dq = D[qr];
                for (int c = 0; c < nS; c++) { var (h1, k1, l1) = sh[c]; fsqv[c + r * nS] = dF[pool.DiffIndex(h1 - h2, k1 - k2, l1 - l2)] * (D[sArr[c]] * dq); }
            }
            // 261109Cl 追加 (研究用、codex2 第 23 回の診断): 裾を半径 DressedSourceSplitRadiusNm で B (内、明示を含む) と A (外) に分け、源の B–B (part 1) か A–A (part 2) だけを積む
            var XtSrc = XtAll; // 261109Cl: 源だけの外殻があればその行も
            int splitPart = double.IsNaN(DressedSourceSplitRadiusNm) || qcEarly != null ? 0 : DressedSourceSplitPart;
            if (splitPart != 0)
            {
                double rs2 = DressedSourceSplitRadiusNm * DressedSourceSplitRadiusNm * (1 + 1e-9);
                XtSrc = Xt.Clone(); var xsv = ((CDense)XtSrc).Values;
                for (int r = 0; r < nQ; r++)
                {
                    int q = qArr[r]; bool outer = Dot3(pool.GX[q], pool.GX[q], pool.GY[q], pool.GY[q], pool.GZ[q], pool.GZ[q]) > rs2;
                    if (outer == (splitPart == 1)) for (int c = 0; c < nM; c++) xsv[r + c * nQ] = 0;
                }
            }
            // var xt32 = DressedDenseSingle ? ToC32(Xt) : null; // 261109Cl 変更前
            // var cross = Vb.ConjugateTransposeThisAndMultiply(MulMaybeSingle(fsq, Xt, xt32)); // 261109Cl 変更前
            // var FW = Vb.ConjugateTransposeThisAndMultiply(fss * Vb) + cross + cross.ConjugateTranspose(); // 261109Cl 変更前
            var xt32 = DressedDenseSingle ? ToC32(XtSrc) : null; // 単精度なら X̃ の変換は 1 回
            var cross = Vb.ConjugateTransposeThisAndMultiply(MulMaybeSingle(fsq, XtSrc, xt32));
            var FW = splitPart == 2 ? new CDense(nM, nM) : Vb.ConjugateTransposeThisAndMultiply(fss * Vb) + cross + cross.ConjugateTranspose(); // 261109Cl: A–A だけなら明示の項と交差を除く
            Lap(3);
            // 求積の位相 (261109Cl 変更、codex2 第 21 回): 方向・電圧に依らない幾何の位相 e^{2πi(hkl·r_a + σ_a g·t_p)} を表に保存して使い回し、
            //   心格子の並進で重なる原子 (Si は 8 → 2 組) をまとめる。方向ごとは列を集めて D_q を掛けるだけ
            var qc = GetQuadCache(table, atomSites, quadAtoms, quadOrder);
            int nSites = qc.NSites;
            // var ph = new CDense(nSites, nQ); var phv = ph.Values; // 261109Cl 変更前
            var ph = new CDense(nSites, nQs); var phv = ph.Values; // 261109Cl: 源だけの外殻も
            for (int r = 0; r < nQs; r++)
            {
                int q = qArrSrc[r]; long src = (long)qc.PoolCol[q] * nSites, dst = (long)r * nSites; double dq = D[q];
                for (int i = 0; i < nSites; i++) phv[dst + i] = qc.Phase[src + i] * dq;
            }
            Lap(4);
            // var chi = MulMaybeSingle(ph, Xt, xt32); // 261109Cl 変更前
            var chi = MulMaybeSingle(ph, XtSrc, xt32);
            var chiW = Matrix<Complex>.Build.DenseOfMatrix(chi);
            var cwv = ((CDense)chiW).Values;
            for (int c = 0; c < nM; c++)
                for (int i = 0; i < nSites; i++) cwv[i + c * nSites] *= qc.Weights[i];
            FW += chi.ConjugateTransposeThisAndMultiply(chiW);
            Lap(5);
            sol.SourceMatrix = Z.ConjugateTransposeThisAndMultiply(FW * Z).ToColumnMajorArray();
            Lap(6);
        }
        return sol;
    }

    /// <summary>261109Cl 追加 (研究用): 裾どうしの疎な結合と密な行列の積 E·X (流束座標: E_wv = U(g_v − g_w) D_w D_v、w, v ∈ Q、w ≠ v)。
    /// X を行優先へ並べ替え、行ごとに CSR の非零 × X の行を AVX の複素の積和で足す。S の列 (map &lt; 0) は飛ばす。</summary>
    // private static Matrix<Complex> SpMMTail(DressedTailTable table, int[] qArr, int[] map, double[] D, Matrix<Complex> X, bool parallel) // 261109Cl 旧シグネチャ
    // private static Matrix<Complex> SpMMTail(DressedTailTable table, int[] qArr, int[] map, double[] D, Matrix<Complex> X, bool parallel, bool source = false) // 261109Cl 旧シグネチャ
    private static Matrix<Complex> SpMMTail(DressedTailTable table, int[] qArr, int[] map, double[] D, Matrix<Complex> X, bool parallel, bool source = false, // 261109Cl source 追加: 後方散乱の源の表 (F) で積む
        bool[] rowKeep = null) // 261109Cl rowKeep 追加: false の行は積まない (結果は 0)
    {
        // if (DressedSpMMSingle) return SpMMTailSingle(table, qArr, map, D, X, parallel); // 261109Cl 変更前
        // if (DressedSpMMSingle) return SpMMTailSingle(table, qArr, map, D, X, parallel, source); // 261109Cl 変更前
        if (DressedSpMMSingle && DressedSpMMKernel == 2 && !source && !parallel) return SpMMTailSingleSoA(table, qArr, map, D, X, rowKeep); // 261109Cl 追加: SoA の核
        if (DressedSpMMSingle) return SpMMTailSingle(table, qArr, map, D, X, parallel, source, rowKeep); // 261109Cl
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
        // if (parallel) Parallel.For(0, nQ, kernel); // 261109Cl 変更前
        // else for (int r = 0; r < nQ; r++) kernel(r);
        if (parallel) Parallel.For(0, nQ, r => { if (rowKeep == null || rowKeep[r]) kernel(r); }); // 261109Cl
        else for (int r = 0; r < nQ; r++) if (rowKeep == null || rowKeep[r]) kernel(r);
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
    // private static Matrix<Complex> SpMMTailSingle(DressedTailTable table, int[] qArr, int[] map, double[] D, Matrix<Complex> X, bool parallel, bool source = false) // 261109Cl 旧シグネチャ
    /// <summary>261109Cl 追加 (研究用、速さ): 単精度の疎な積で、実部と虚部を別の並び (SoA) にした形。非零ごとに X の行を流し、複素 8 個を FMA 4 回で積む
    /// (AoS の並べ替え・AddSubtract が要らない)。DressedSpMMKernel = 2 で選ぶ (源の表には使わない)。値は SpMMTailSingle と丸めの範囲で同じ。</summary>
    private static Matrix<Complex> SpMMTailSingleSoA(DressedTailTable table, int[] qArr, int[] map, double[] D, Matrix<Complex> X, bool[] rowKeep)
    {
        lock (table)
            if (table.vReF == null)
            {
                table.vReF = table.vRe.Select(v => (float)v).ToArray();
                table.vImF = table.vIm.Select(v => (float)v).ToArray();
            }
        int nQ = X.RowCount, k = X.ColumnCount, kp = (k + 15) & ~15;
        var xcm = ((CDense)X).Values;
        var xre = new float[(long)nQ * kp]; var xim = new float[(long)nQ * kp];
        for (int r0 = 0; r0 < nQ; r0 += 32)
            for (int c0 = 0; c0 < k; c0 += 32)
                for (int r = r0, rEnd = Math.Min(r0 + 32, nQ); r < rEnd; r++)
                {
                    double dv = D[qArr[r]];
                    for (int c = c0, cEnd = Math.Min(c0 + 32, k); c < cEnd; c++)
                    {
                        var x = xcm[r + c * nQ];
                        xre[(long)r * kp + c] = (float)(x.Real * dv); xim[(long)r * kp + c] = (float)(x.Imaginary * dv);
                    }
                }
        var yre = new float[(long)nQ * kp]; var yim = new float[(long)nQ * kp];
        var rowPtr = table.rowPtr; var col = table.colIdx; var vr = table.vReF; var vi = table.vImF;
        unsafe
        {
            fixed (float* xr0 = xre, xi0 = xim, yr0 = yre, yi0 = yim, vrp = vr, vip = vi)
            fixed (int* colp = col, mapp = map)
            {
                for (int r = 0; r < nQ; r++)
                {
                    if (rowKeep != null && !rowKeep[r]) continue;
                    int pw = qArr[r], p0 = rowPtr[pw], p1 = rowPtr[pw + 1];
                    float* Yr = yr0 + (long)r * kp, Yi = yi0 + (long)r * kp;
                    for (int p = p0; p < p1; p++)
                    {
                        int qv = mapp[colp[p]];
                        if (qv < 0) continue;
                        float a = vrp[p], b = vip[p];
                        float* Xr = xr0 + (long)qv * kp, Xi = xi0 + (long)qv * kp;
                        if (Fma.IsSupported && Avx.IsSupported)
                        {
                            var va = Vector256.Create(a); var vb = Vector256.Create(b);
                            for (int c = 0; c < kp; c += 16)
                            {
                                var xr = Avx.LoadVector256(Xr + c); var xi = Avx.LoadVector256(Xi + c);
                                var xr2 = Avx.LoadVector256(Xr + c + 8); var xi2 = Avx.LoadVector256(Xi + c + 8);
                                var yr = Avx.LoadVector256(Yr + c); var yi = Avx.LoadVector256(Yi + c);
                                var yr2 = Avx.LoadVector256(Yr + c + 8); var yi2 = Avx.LoadVector256(Yi + c + 8);
                                yr = Fma.MultiplyAddNegated(vb, xi, Fma.MultiplyAdd(va, xr, yr));
                                yi = Fma.MultiplyAdd(vb, xr, Fma.MultiplyAdd(va, xi, yi));
                                yr2 = Fma.MultiplyAddNegated(vb, xi2, Fma.MultiplyAdd(va, xr2, yr2));
                                yi2 = Fma.MultiplyAdd(vb, xr2, Fma.MultiplyAdd(va, xi2, yi2));
                                Avx.Store(Yr + c, yr); Avx.Store(Yi + c, yi); Avx.Store(Yr + c + 8, yr2); Avx.Store(Yi + c + 8, yi2);
                            }
                        }
                        else
                            for (int c = 0; c < k; c++) { float xr = Xr[c], xi = Xi[c]; Yr[c] += a * xr - b * xi; Yi[c] += a * xi + b * xr; }
                    }
                }
            }
        }
        var res = new CDense(nQ, k);
        var rv = res.Values;
        for (int c0 = 0; c0 < k; c0 += 32)
            for (int r0 = 0; r0 < nQ; r0 += 32)
                for (int c = c0, cEnd = Math.Min(c0 + 32, k); c < cEnd; c++)
                    for (int r = r0, rEnd = Math.Min(r0 + 32, nQ); r < rEnd; r++)
                    {
                        double dw = D[qArr[r]];
                        rv[r + c * nQ] = new Complex(yre[(long)r * kp + c] * dw, yim[(long)r * kp + c] * dw);
                    }
        return res;
    }

    private static Matrix<Complex> SpMMTailSingle(DressedTailTable table, int[] qArr, int[] map, double[] D, Matrix<Complex> X, bool parallel, bool source = false, bool[] rowKeep = null) // 261109Cl source・rowKeep 追加
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
        // 261109Cl 追加 (研究用、速さ): 列の塊ごとに全行を回す (X の列の塊 = 裾 × CB 個の float がキャッシュに載る。非零ごとに X の行の全部を読み直す帯域の律速を避ける)
        unsafe void RowFB(int r, int c0, int c1)
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
                    int c = c0;
                    if (Fma.IsSupported && Avx.IsSupported)
                    {
                        var va = Vector256.Create(a); var vb = Vector256.Create(b);
                        for (; c + 32 <= c1; c += 32)
                        {
                            var x0 = Avx.LoadVector256(xs + c); var x1 = Avx.LoadVector256(xs + c + 8); var x2 = Avx.LoadVector256(xs + c + 16); var x3 = Avx.LoadVector256(xs + c + 24);
                            Avx.Store(Y + c, Avx.AddSubtract(Fma.MultiplyAdd(va, x0, Avx.LoadVector256(Y + c)), Avx.Multiply(vb, Avx.Permute(x0, 0b10110001))));
                            Avx.Store(Y + c + 8, Avx.AddSubtract(Fma.MultiplyAdd(va, x1, Avx.LoadVector256(Y + c + 8)), Avx.Multiply(vb, Avx.Permute(x1, 0b10110001))));
                            Avx.Store(Y + c + 16, Avx.AddSubtract(Fma.MultiplyAdd(va, x2, Avx.LoadVector256(Y + c + 16)), Avx.Multiply(vb, Avx.Permute(x2, 0b10110001))));
                            Avx.Store(Y + c + 24, Avx.AddSubtract(Fma.MultiplyAdd(va, x3, Avx.LoadVector256(Y + c + 24)), Avx.Multiply(vb, Avx.Permute(x3, 0b10110001))));
                        }
                        for (; c + 8 <= c1; c += 8)
                        {
                            var x0 = Avx.LoadVector256(xs + c);
                            Avx.Store(Y + c, Avx.AddSubtract(Fma.MultiplyAdd(va, x0, Avx.LoadVector256(Y + c)), Avx.Multiply(vb, Avx.Permute(x0, 0b10110001))));
                        }
                    }
                    for (; c < c1; c += 2)
                    {
                        float xre = xs[c], xim = xs[c + 1];
                        Y[c] += a * xre - b * xim; Y[c + 1] += a * xim + b * xre;
                    }
                }
            }
        }
        int cbF = DressedSpMMColumnBlock * 2; // 複素の個数 → float の個数
        if (cbF > 0 && cbF < k2)
        {
            for (int c0 = 0; c0 < k2; c0 += cbF)
            {
                int c1 = Math.Min(c0 + cbF, k2);
                if (parallel) Parallel.For(0, nQ, r => { if (rowKeep == null || rowKeep[r]) RowFB(r, c0, c1); });
                else for (int r = 0; r < nQ; r++) if (rowKeep == null || rowKeep[r]) RowFB(r, c0, c1);
            }
        }
        else
        {
        // if (parallel) Parallel.For(0, nQ, RowF); // 261109Cl 変更前
        // else for (int r = 0; r < nQ; r++) RowF(r);
        if (parallel) Parallel.For(0, nQ, r => { if (rowKeep == null || rowKeep[r]) RowF(r); }); // 261109Cl
        else for (int r = 0; r < nQ; r++) if (rowKeep == null || rowKeep[r]) RowF(r);
        }
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

    /// <summary>261109Cl 追加 (研究用): ブロックの全方向で繰り上げるとき、各方向の E1 の和集合も明示にする (既定 false、試験用)</summary>
    public static bool DressedBlockUnionE1 { get; set; } = false;
    /// <summary>261109Cl 追加 (研究用、native の試作): 1 方向の選択と密なブロック (流束座標、H_SS の対角 = λ)。速い形の部品を外 (研究用の native) で組むために返す。</summary>
    public sealed class DressedTailBlocks
    {
        public int[] SArr, QArr, Map; public double[] P, D; public int NSelected;
        public CDense HSS, HQS, HSQ; public Complex[] LamQ;
    }

    /// <summary>261109Cl 追加 (研究用、native の試作): 選択と密なブロックだけを作る</summary>
    public DressedTailBlocks BuildDressedTailBlocksResearch(DressedTailTable table, Vector3DBase vecK0, Vector3DBase surface, int n, double ps)
    {
        var cr = DressCore(table, vecK0, surface, n, ps, 45, true, false, null, null, null, false, true);
        return new DressedTailBlocks { SArr = cr.SArr, QArr = cr.QArr, Map = cr.Map, P = cr.P, D = cr.D, NSelected = cr.NSelected,
            HSS = (CDense)cr.HSS, HQS = (CDense)cr.HQS, HSQ = (CDense)cr.HSQ, LamQ = cr.LamQ };
    }

    /// <summary>261109Cl 追加 (研究用、native の試作): 裾どうしの疎な積 E·X (流束座標、DressedSpMMSingle などの静的な選択に従う)</summary>
    public static CDense DressedTailSpMMResearch(DressedTailTable table, DressedTailBlocks b, CDense X)
        => (CDense)SpMMTail(table, b.QArr, b.Map, b.D, X, false);

    /// <summary>261109Cl 追加 (研究用、native の試作): 裾どうしの結合の表 (CSR、単精度の値)。行・列は候補の番号、要素 (i, j) = U(g_j − g_i)</summary>
    /// <summary>261110Cl 追加 (研究用、配布の設計の試作): 裾どうしの結合は差 d = g_v − g_w だけの関数 U(d) (流束座標の D は前後で掛ける) なので、E·X は逆格子の上の畳み込み (stencil)。
    /// |U(d)| ≥ η max|U| (d ≠ 0、表の CSR と同じ条件) の差の一覧 (h, k, l の 3 つ組) と U(d) (単精度) を返す。CSR を持たずに疎な積ができる (表のメモリがほぼ要らない)</summary>
    public static (int[] dHkl, float[] re, float[] im) DressedTailStencilResearch(DressedTailTable table)
    {
        var pool = table.Pool; double thr = table.Eta * table.MaxAbsU;
        var dh = new List<int>(); var re = new List<float>(); var im = new List<float>();
        for (int h = -pool.DiffBoundH; h <= pool.DiffBoundH; h++)
            for (int k = -pool.DiffBoundK; k <= pool.DiffBoundK; k++)
                for (int l = -pool.DiffBoundL; l <= pool.DiffBoundL; l++)
                {
                    if (h == 0 && k == 0 && l == 0) continue;
                    var u = table.diffUfull[pool.DiffIndex(h, k, l)];
                    if (!(u.Magnitude >= thr)) continue;
                    dh.Add(h); dh.Add(k); dh.Add(l); re.Add((float)u.Real); im.Add((float)u.Imaginary);
                }
        return ([.. dh], [.. re], [.. im]);
    }

    /// <summary>261110Cl 追加 (研究用): 上の stencil の疎な積の方向ごとの入力。裾の行の (h, k, l)、(h, k, l) → 裾の行 (無ければ −1) の 3 次元の表 (範囲 ±H, ±K, ±L)、裾の行の D</summary>
    public static (int[] rowHkl, int[] map3D, int H, int K, int L, double[] Dq) DressedTailStencilMapResearch(DressedTailTable table, DressedTailBlocks b)
    {
        var pool = table.Pool; int H = 0, K = 0, L = 0;
        for (int i = 0; i < pool.Count; i++) { var (h, k, l) = pool.Hkl[i]; H = Math.Max(H, Math.Abs(h)); K = Math.Max(K, Math.Abs(k)); L = Math.Max(L, Math.Abs(l)); }
        var map3D = new int[(2 * H + 1) * (2 * K + 1) * (2 * L + 1)]; Array.Fill(map3D, -1);
        int nQ = b.QArr.Length; var rowHkl = new int[3 * nQ]; var Dq = new double[nQ];
        for (int r = 0; r < nQ; r++)
        {
            var (h, k, l) = pool.Hkl[b.QArr[r]];
            rowHkl[3 * r] = h; rowHkl[3 * r + 1] = k; rowHkl[3 * r + 2] = l; Dq[r] = b.D[b.QArr[r]];
            map3D[((h + H) * (2 * K + 1) + k + K) * (2 * L + 1) + l + L] = r;
        }
        return (rowHkl, map3D, H, K, L, Dq);
    }

    public static (int[] rowPtr, int[] col, float[] re, float[] im) DressedTailCsrResearch(DressedTailTable table)
    {
        lock (table)
            if (table.vReF == null)
            {
                table.vReF = table.vRe.Select(v => (float)v).ToArray();
                table.vImF = table.vIm.Select(v => (float)v).ToArray();
            }
        return (table.rowPtr, table.colIdx, table.vReF, table.vImF);
    }

    /// <summary>261109Cl 追加 (研究用、native の試作): 後方 TDS の源の密なブロック F'_SS (明示 × 明示) と F'_SQ (明示 × 裾) (流束座標)</summary>
    public static (CDense fss, CDense fsq) DressedTailSourceBlocksResearch(DressedTailTable table, DressedTailBlocks b)
    {
        var pool = table.Pool; var sArr = b.SArr; var qArr = b.QArr; var D = b.D; int nS = sArr.Length, nQ = qArr.Length;
        var fss = new CDense(nS, nS); var fsq = new CDense(nS, nQ); var fssv = fss.Values; var fsqv = fsq.Values; var dF = table.diffF;
        var sh = new (int h, int k, int l)[nS];
        for (int i = 0; i < nS; i++) sh[i] = pool.Hkl[sArr[i]];
        for (int c = 0; c < nS; c++)
        {
            var (h2, k2, l2) = sh[c]; double dc = D[sArr[c]];
            for (int r = 0; r < nS; r++) { var (h1, k1, l1) = sh[r]; fssv[r + c * nS] = dF[pool.DiffIndex(h1 - h2, k1 - k2, l1 - l2)] * (D[sArr[r]] * dc); }
        }
        for (int r = 0; r < nQ; r++)
        {
            int qr = qArr[r]; var (h2, k2, l2) = pool.Hkl[qr]; double dq = D[qr];
            for (int c = 0; c < nS; c++) { var (h1, k1, l1) = sh[c]; fsqv[c + r * nS] = dF[pool.DiffIndex(h1 - h2, k1 - k2, l1 - l2)] * (D[sArr[c]] * dq); }
        }
        return (fss, fsq);
    }

    /// <summary>261109Cl 追加 (研究用、native の試作): 源の求積の位相 (点 × 裾、D_q を掛けたもの) と点の重み</summary>
    public static (Complex[] ph, double[] w, int nSites) DressedTailCubatureResearch(DressedTailTable table, DressedTailBlocks b, (double x, double y, double z)[] atomSites, (double weight, double sigmaNm)[] quadAtoms, int order)
    {
        var qc = GetQuadCache(table, atomSites, quadAtoms, order);
        int nSites = qc.NSites, nQ = b.QArr.Length;
        var ph = new Complex[(long)nSites * nQ];
        for (int r = 0; r < nQ; r++)
        {
            int q = b.QArr[r]; long src = (long)qc.PoolCol[q] * nSites, dst = (long)r * nSites; double dq = b.D[q];
            for (int i = 0; i < nSites; i++) ph[dst + i] = qc.Phase[src + i] * dq;
        }
        return (ph, qc.Weights, nSites);
    }

    /// <summary>261109Cl 追加 (研究用、速さ、codex2 第 21 回): 速い形で疎な積 1 回のとき、X̃†(R + ΛX̃) = (X̃†X̃)Γ + X̃†E_δ を使う閾値 (E_δ の倍率 δ²/(|λ − γ|² + δ²) がこれを超える要素だけ。NaN = 使わない、0 = 全部)</summary>
    public static double DressedEdeltaThreshold { get; set; } = double.NaN;
    /// <summary>261109Cl 追加 (研究用、速さ): 着衣の疎な積 E·X₀ を |λ_q| がこれより小さい裾の行だけで積む (+∞ = 全部)</summary>
    public static double DressedSpMMLambdaCut { get; set; } = double.PositiveInfinity;
    /// <summary>261109Cl 追加 (研究用): 裾を「芯の半径 DressedTailCoreRadiusNm の中の全部」+「その外 (表の裾の半径まで) の |Q/P| がこれより小さいもの」にする (NaN = 使わない)</summary>
    public static double DressedTailShellLambdaCut { get; set; } = double.NaN;
    /// <summary>261109Cl 追加 (研究用): 上の芯の半径 [nm⁻¹]</summary>
    public static double DressedTailCoreRadiusNm { get; set; } = 40;
    /// <summary>261109Cl 追加 (研究用、速さ、native の試作と同じ規則): 速い形の着衣のモードの打ち切り。近縮退 (|γ_i − γ_j| &lt; DressedModeDegTol、推移的に) の組 C の
    /// 重み ‖V_C a_C‖ (a = V_b⁻¹e₀) が最大のこの割合より小さい組を丸ごと捨てる (0 = 打ち切らない)</summary>
    public static double DressedModeTau { get; set; } = 0;
    /// <summary>261109Cl 追加 (研究用): 速い形の後方 TDS の源に、S ∪ Q の外でこの半径以内 (表の裾の半径まで) のビームの一次の振幅を足す [nm⁻¹] (NaN = 使わない)</summary>
    public static double DressedSourceShellRadiusNm { get; set; } = double.NaN;
    /// <summary>261109Cl 追加 (研究用): 上の外殻のビームの数の積算</summary>
    public static long DressedSourceShellCount;
    /// <summary>261109Cl 追加 (研究用): 裾を「芯の半径 DressedTailCoreRadiusNm の中の全部」+「明示のビームからこの距離以内 (表の裾の半径まで)」にする [nm⁻¹] (NaN = 使わない)</summary>
    public static double DressedTailNeighborRadiusNm { get; set; } = double.NaN;
    /// <summary>261109Cl 追加 (研究用): 上の近傍の中心に、明示に加えて候補のうち |Q/P| がこれより小さい反射も含める (NaN = 明示だけ)</summary>
    public static double DressedTailNeighborPs { get; set; } = double.NaN;
    /// <summary>261110Cl 追加 (研究用、速さ): 裾の行を (h, k, l) の Morton 順に並べる (false = 候補の番号の順)</summary>
    public static bool DressedTailMortonOrder { get; set; } = false;
    /// <summary>261110Cl 追加 (研究用): 疎な表 (CSR) を作るか (false = stencil の疎な積だけ使う研究の経路。C# の疎な積・CSR の近傍は使えない)</summary>
    public static bool DressedTailBuildCsr { get; set; } = true;
    /// <summary>261110Cl 追加 (研究用): 近傍を距離の球でなく、中心から表の CSR の列 (|U(d)| ≥ η max) で一段つながる裾にする (DressedTailNeighborRadiusNm は NaN 以外で有効化だけに使う)</summary>
    public static bool DressedTailNeighborUseCsr { get; set; } = false;
    /// <summary>261109Cl 追加 (研究用): 上の近縮退の判定の幅 (0 = 組にしない)</summary>
    public static double DressedModeDegTol { get; set; } = 1e-4;
    /// <summary>261109Cl 追加 (研究用、診断): 速い形の後方 TDS の源を裾の半径で分けるときの境界 [nm⁻¹] (NaN = 分けない)</summary>
    public static double DressedSourceSplitRadiusNm { get; set; } = double.NaN;
    /// <summary>261109Cl 追加 (研究用、診断): 0 = 全部、1 = 内側 (明示と境界の内の裾) どうし、2 = 外側の裾どうしだけ。交差 = 全部 − 1 − 2</summary>
    public static int DressedSourceSplitPart { get; set; } = 0;
    /// <summary>261109Cl 追加 (研究用、速さ): 単精度の疎な積の列の塊の大きさ (複素の個数、0 = 塊にしない)</summary>
    public static int DressedSpMMColumnBlock { get; set; } = 0;
    /// <summary>261109Cl 追加 (研究用): E_δ の疎な積で使った要素の数の積算</summary>
    public static long DressedEdeltaNonzeros;

    /// <summary>261109Cl 追加: 源の求積の幾何の位相の表 (表 1 つにつき 1 つ)</summary>
    internal sealed class DressedQuadCache
    {
        public int Order, NSites; public (double x, double y, double z)[] AtomSites; public (double weight, double sigmaNm)[] QuadAtoms;
        public int[] PoolCol; public Complex[] Phase; public double[] Weights; public int Groups;
    }

    /// <summary>261109Cl 追加 (codex2 第 21 回): 源の求積の点と幾何の位相の表を作る (裾の半径の中の候補ごとに e^{2πi(hkl·r_a + σ_a g·t_p)})。
    /// 同じ σ・同じ重みで、裾の全候補について hkl·(r_a − r_b) が整数になる原子 (心格子の並進で重なる) は 1 組にまとめ重みを足す (厳密)。</summary>
    private static DressedQuadCache GetQuadCache(DressedTailTable table, (double x, double y, double z)[] atomSites, (double weight, double sigmaNm)[] quadAtoms, int order)
    {
        lock (table)
        {
            if (table.quadCache is { } c0 && c0.Order == order && c0.AtomSites.SequenceEqual(atomSites) && c0.QuadAtoms.SequenceEqual(quadAtoms)) return c0;
            var pool = table.Pool; int n = pool.Count;
            double rtt = table.TailRadiusNm * table.TailRadiusNm * (1 + 1e-9);
            var col = new int[n]; Array.Fill(col, -1); var tail = new List<int>();
            for (int i = 0; i < n; i++)
                if (Dot3(pool.GX[i], pool.GX[i], pool.GY[i], pool.GY[i], pool.GZ[i], pool.GZ[i]) <= rtt) { col[i] = tail.Count; tail.Add(i); }
            int na = atomSites.Length;
            var reps = new List<int>(); var gw = new List<double>();
            for (int a = 0; a < na; a++)
            {
                int found = -1;
                for (int g = 0; g < reps.Count && found < 0; g++)
                {
                    int b = reps[g];
                    if (quadAtoms[b].sigmaNm != quadAtoms[a].sigmaNm || quadAtoms[b].weight != quadAtoms[a].weight) continue;
                    double dx = atomSites[a].x - atomSites[b].x, dy = atomSites[a].y - atomSites[b].y, dz = atomSites[a].z - atomSites[b].z;
                    bool same = true;
                    foreach (var i in tail)
                    {
                        var (h, k, l) = pool.Hkl[i]; double v = h * dx + k * dy + l * dz;
                        if (Math.Abs(v - Math.Round(v)) > 1e-9) { same = false; break; }
                    }
                    if (same) found = g;
                }
                if (found >= 0) gw[found] += quadAtoms[a].weight;
                else { reps.Add(a); gw.Add(quadAtoms[a].weight); }
            }
            var (pts, wpt) = GaussianCubature3D(order);
            int np = pts.Length, ng = reps.Count, nSites = ng * np;
            var w = new double[nSites];
            for (int g = 0; g < ng; g++) for (int t = 0; t < np; t++) w[g * np + t] = gw[g] * wpt[t];
            var phase = new Complex[(long)tail.Count * nSites];
            Parallel.For(0, tail.Count, ci =>
            {
                int i = tail[ci]; var (h, k, l) = pool.Hkl[i]; double gx = pool.GX[i], gy = pool.GY[i], gz = pool.GZ[i];
                long o = (long)ci * nSites;
                for (int g = 0; g < ng; g++)
                {
                    var (ax, ay, az) = atomSites[reps[g]]; double s2 = 2 * Math.PI * quadAtoms[reps[g]].sigmaNm, bph = 2 * Math.PI * (h * ax + k * ay + l * az);
                    for (int t = 0; t < np; t++) { var (tx, ty, tz) = pts[t]; phase[o + g * np + t] = Complex.FromPolarCoordinates(1, bph + s2 * (gx * tx + gy * ty + gz * tz)); }
                }
            });
            var cache = new DressedQuadCache { Order = order, NSites = nSites, AtomSites = [.. atomSites], QuadAtoms = [.. quadAtoms], PoolCol = col, Phase = phase, Weights = w, Groups = ng };
            table.quadCache = cache;
            return cache;
        }
    }

    /// <summary>261109Cl 追加 (研究用、速さ): 着衣 (H_SQ X₀、(H_QS + E X₀)V_b) と源 (F'_SQ X̃、求積の積) の密な積を単精度 (complex64) にする (既定 false)。
    /// Galerkin の K・G の積は倍精度のまま (Im γ の誤差は深さに比例して効くため)。</summary>
    public static bool DressedDenseSingle { get; set; } = false;

    /// <summary>261109Cl 追加: 複素の密な行列を単精度へ</summary>
    private static MathNet.Numerics.LinearAlgebra.Complex32.DenseMatrix ToC32(Matrix<Complex> m)
    {
        var src = m is CDense d ? d.Values : ((CDense)Matrix<Complex>.Build.DenseOfMatrix(m)).Values;
        var dst = new MathNet.Numerics.Complex32[src.Length];
        for (int i = 0; i < src.Length; i++) dst[i] = new MathNet.Numerics.Complex32((float)src[i].Real, (float)src[i].Imaginary);
        return new MathNet.Numerics.LinearAlgebra.Complex32.DenseMatrix(m.RowCount, m.ColumnCount, dst);
    }

    /// <summary>261109Cl 追加: a·b を、DressedDenseSingle なら単精度で積んで倍精度へ戻す (b32 = b の単精度を渡すと変換を省く)</summary>
    private static Matrix<Complex> MulMaybeSingle(Matrix<Complex> a, Matrix<Complex> b, MathNet.Numerics.LinearAlgebra.Complex32.DenseMatrix b32 = null)
    {
        if (!DressedDenseSingle) return a * b;
        var p = ToC32(a) * (b32 ?? ToC32(b));
        var pv = ((MathNet.Numerics.LinearAlgebra.Complex32.DenseMatrix)p).Values;
        var r = new Complex[pv.Length];
        for (int i = 0; i < pv.Length; i++) r[i] = new Complex(pv[i].Real, pv[i].Imaginary);
        return new CDense(p.RowCount, p.ColumnCount, r);
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

    /// <summary>261109Cl 追加 (codex2 第 19 回): 3 次元の標準正規分布の正の重みの求積 (球面の Lebedev × 半径の Gauss)。
    /// 52 点 = 次数 7 (Lebedev 26 方向 × r² = 5 ∓ √10)、77 点 = 次数 9 (Lebedev 38 方向 × r² = 7 ∓ √14 + 原点)。重みの和 = 1。作るときに 2 次・4 次のモーメントを確かめる。</summary>
    private static ((double x, double y, double z)[] pts, double[] w) GaussianCubatureSpherical(int points)
    {
        var dirs = new List<((double, double, double) d, double w)>();
        void AddPerm(double a, double b, double c, double w) // 符号と並べ替えの全部 (重複は除く)
        {
            var set = new HashSet<(double, double, double)>();
            foreach (var (x, y, z) in new[] { (a, b, c), (a, c, b), (b, a, c), (b, c, a), (c, a, b), (c, b, a) })
                foreach (var sx in new[] { 1, -1 }) foreach (var sy in new[] { 1, -1 }) foreach (var sz in new[] { 1, -1 })
                            set.Add((sx * x + 0.0, sy * y + 0.0, sz * z + 0.0));
            foreach (var d in set) dirs.Add((d, w));
        }
        double r3 = 1 / Math.Sqrt(3), r2 = 1 / Math.Sqrt(2);
        double w0 = 0; double[] u, wu;
        if (points == 52)
        {
            AddPerm(1, 0, 0, 1.0 / 21); AddPerm(0, r2, r2, 4.0 / 105); AddPerm(r3, r3, r3, 9.0 / 280);
            double u1 = 5 - Math.Sqrt(10), u2 = 5 + Math.Sqrt(10);
            u = [u1, u2]; wu = [(u2 - 3) / (u2 - u1), (3 - u1) / (u2 - u1)];
        }
        else if (points == 77)
        {
            AddPerm(1, 0, 0, 1.0 / 105); AddPerm(r3, r3, r3, 9.0 / 280); AddPerm(0.4597008433809831, 0.8880738339771153, 0, 1.0 / 35);
            u = [7 - Math.Sqrt(14), 7 + Math.Sqrt(14)]; wu = [0.706427776009195, 0.065000795419376]; w0 = 0.228571428571429;
        }
        else throw new ArgumentException("Only 52 or 77 points", nameof(points));
        var pts = new List<(double, double, double)>(); var w = new List<double>();
        if (w0 > 0) { pts.Add((0, 0, 0)); w.Add(w0); }
        for (int i = 0; i < u.Length; i++)
        {
            double r = Math.Sqrt(u[i]);
            foreach (var ((dx, dy, dz), wd) in dirs) { pts.Add((r * dx, r * dy, r * dz)); w.Add(wu[i] * wd); }
        }
        double sw = w.Sum();
        for (int i = 0; i < w.Count; i++) w[i] /= sw;
        double m2 = 0, m4 = 0, m22 = 0; // E[x²] = 1、E[x⁴] = 3、E[x²y²] = 1
        for (int i = 0; i < w.Count; i++) { var (x, y, _) = pts[i]; m2 += w[i] * x * x; m4 += w[i] * x * x * x * x; m22 += w[i] * x * x * y * y; }
        if (pts.Count != points || Math.Abs(sw - 1) > 1e-12 || Math.Abs(m2 - 1) > 1e-12 || Math.Abs(m4 - 3) > 1e-11 || Math.Abs(m22 - 1) > 1e-11)
            throw new InvalidOperationException($"Cubature check failed: {pts.Count} points, Σw {sw}, E[x²] {m2}, E[x⁴] {m4}, E[x²y²] {m22}");
        return ([.. pts], [.. w]);
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
        /// <summary>261109Cl 追加: 候補のプールの半径 [nm⁻¹] (E1 の選択と遅い反射の繰り上げはこの半径の全部から、裾は TailRadiusNm まで)。
        /// Si 20 kV の 2 方向で、80:80 (候補:裾) に対し 40:40 は晶帯軸で −2.15 %、64:40 は −0.31 % / −0.17 % (1 方向 2.2 s)、64:64 は +0.03 % (15 s)。既定 64</summary>
        public double PoolRadiusNm { get; set; } = 64;
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
        /// <summary>方向のブロックの一辺 (格子の画素の数。1 = 方向ごとに着衣、2 = 2×2 で共有し中心は画素の境目)。
        /// 261109Cl: 既定を 1 に (G512 の 2×2 は中心から最大 4.5 mrad で、集合を固定した純粋な共有の上乗せが 32 升すべてで 0.1 % 超 (0.12〜0.31 %)。
        /// 追加の基底 (着衣の差 Δ の圧縮、r ≤ 40) も特異値がほとんど減らず効かない。codex2 第 20 回の案 3 = 着衣は方向ごと)。2 は速さ優先の選択肢</summary>
        // public int BlockSize { get; set; } = 2; // 261109Cl 変更前
        public int BlockSize { get; set; } = 1;
        /// <summary>後方 TDS 行列 (IncludeTDSBackground) の源の表の疎化の閾値 (|F(d)| ≥ SourceEta·|F(0)|、0 = 差の箱の全部)</summary>
        public double SourceEta { get; set; } = 1e-3;
        /// <summary>源の裾どうしの項を Gauss–Hermite の求積で積むときの 1 軸の点の数 (0 = 源の表の疎な積。求積なら源の表は作らない)。
        /// Si 20 kV・G16・36 方向で、全部の F に対する差の最大: 1 点 0.96 %、3 点 0.18 %、5 点 0.022 % (η_F 0.1 の疎な積 0.037 %・1 方向 +15 s に対し +1 s)。
        /// 52 / 77 は球面 (Lebedev) × 半径の規則 (次数 7 / 9): 0.077 % / 0.020 %。7 点 (343) は 0.013 % (残りは Gauss の包絡と正確な F の差)。既定 77 (125 点と同じ精度で点は 62 %)</summary>
        // public int SourceTailQuadrature { get; set; } = 5; // 261109Cl 変更前
        public int SourceTailQuadrature { get; set; } = 77;
        /// <summary>ブロックのどれかの方向で遅い (|Q/P| &lt; Ps) 裾の反射も明示へ繰り上げる (G512 の 2×2 は中心から最大 4.45 mrad あり、中心だけで判定すると共有の上乗せが 1 % 級。
        /// 繰り上げると 0.4 %、各方向の E1 の和集合 (<see cref="DressedBlockUnionE1"/>) も足すと 0.3 %、day5p)</summary>
        public bool PromoteOverBlock { get; set; } = true;
        /// <summary>261109Cl 追加 (研究用): 1 方向を解く外の解き手 (null = SolveDressedTailSingleResearch)。引数 = (BetheMethod、表、k0、表面の法線、原子位置、源の求積の原子、後方 TDS か)。
        /// 返す解は Gamma と、後方 TDS なら SourceMatrix、局所源なら AtomAmp を持つこと。並列で呼ばれる</summary>
        public Func<BetheMethod, DressedTailTable, Vector3DBase, Vector3DBase, (double x, double y, double z)[], (double weight, double sigmaNm)[], bool, DressedTailSolution> SingleSolver { get; set; }
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
        // if (pool == null || pool.RadiusNm != opt.TailRadiusNm || pool.Mat != mat) // 261109Cl 変更前
        //     pool = CreateBeamSelectionPool(opt.TailRadiusNm, buildE0: false);
        double poolRadius = Math.Max(opt.PoolRadiusNm, opt.TailRadiusNm); // 261109Cl: 候補の半径と裾の半径を分ける
        if (pool == null || pool.RadiusNm != poolRadius || pool.Mat != mat)
            pool = CreateBeamSelectionPool(poolRadius, buildE0: false); // E1 の尺度は U/γ (電圧に依らない) なので、エネルギーの間で共有する
        rec.SecondsPool = sw.Elapsed.TotalSeconds; sw.Restart();
        bool tds = background != null;
        // var table = CreateDressedTailTable(pool, opt.Eta, !tds ? double.NaN : opt.SourceTailQuadrature > 0 ? double.PositiveInfinity : opt.SourceEta); // 261109Cl 変更前
        var table = CreateDressedTailTable(pool, opt.Eta, !tds ? double.NaN : opt.SourceTailQuadrature > 0 ? double.PositiveInfinity : opt.SourceEta, opt.TailRadiusNm); // 吸収は電圧で変わるので、表はエネルギーごとに作り直す
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
            // 261109Cl 追加 (codex2 第 20 回の案 3): 自分を中心にする 1 方向は軽い形 (共有の G 群を作らない。同じ値で約 2 割速い)
            DressedTailSolution Single((Vector3DBase k0, Vector3DBase n) d)
                // => SolveDressedTailSingleResearch(table, d.k0, d.n, opt.N, opt.Ps, opt.RegDepthNm, opt.OneProduct, atomArray, tds, quadAtoms, opt.SourceTailQuadrature); // 261109Cl 変更前
                => opt.SingleSolver != null ? opt.SingleSolver(this, table, d.k0, d.n, atomArray, quadAtoms, tds) // 261109Cl 追加 (研究用): 外の解き手 (native の試作など)
                : SolveDressedTailSingleResearch(table, d.k0, d.n, opt.N, opt.Ps, opt.RegDepthNm, opt.OneProduct, atomArray, tds, quadAtoms, opt.SourceTailQuadrature);
            if (members.Count == 1 && valid == 1)
                sols = [Single(dirs[0])];
            else
            try
            {
                sols = SolveDressedTailBlockResearch(table, ck0, cn, dirs, opt.N, opt.Ps, opt.RegDepthNm, opt.OneProduct, false, null, null, atomArray, tds, quadAtoms, opt.SourceTailQuadrature, opt.PromoteOverBlock);
            }
            catch (DressedTailBackwardBeamException) when (members.Count > 1)
            {
                // P のクリップやビームの削除はせず、方向ごとに (自分を中心に) 解き直す
                Interlocked.Increment(ref split);
                sols = new DressedTailSolution[members.Count];
                for (int m = 0; m < members.Count; m++)
                    // sols[m] = SolveDressedTailBlockResearch(table, dirs[m].k0, dirs[m].n, [dirs[m]], opt.N, opt.Ps, opt.RegDepthNm, opt.OneProduct, false, null, null, atomArray, tds, quadAtoms, opt.SourceTailQuadrature)[0]; // 261109Cl 変更前
                    sols[m] = Single(dirs[m]); // 261109Cl
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
