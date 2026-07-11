// 260709Cl 新規 (Phase N1): affine normalizer エンジン。
// 目的: isomorphic (IIc) 部分群の「ITA A1 流の系列表示」に向けた基盤 — 現在の G-共役類別では P1 の
// index≤8 で 108 行に爆発する (codex R7 で保留)。N_Aff(G) の軌道で束ねるための normalizer データ
// (離散線形部の生成集合 + 純並進核) を、既存の型同定機構と同じ「self-identification」方式で計算する。
// 設計は codex 相談 R9 で確定 (.project-guidance/ReciPro_FormGroupRelations改修計画.md):
//   - n=(U,a) ∈ N_Aff(G) ⟺ U·A_i·U⁻¹ = A_{π(i)} (置換 π) かつ (I−A_{π(i)})·a ≡ t_{π(i)} − U·t_i (mod Z³) ∀i。
//     この rectangular な連立合同式を Smith 標準形 (SNF) で解くのが中核 (可解性・particular 解・解核が同時に出る)。
//   - 純並進核 K = {a | (I−A_i)a ∈ Z³ ∀i} ≅ (R/Z)^{3−r} × Π Z/d_i は同じソルバの斉次版。
//     連続部分 (polar 方向) は G の全元を点ごとに centralize するため軌道計算では恒等 — メタデータとして保持。
//   - 線形部候補は SmallUnimodular(k) ∩ 正規化フィルタの有界探索。「k=1 で完全」の一般保証は無い
//     (codex R9 が centralizer の合同条件による反例を提示) ため、Completeness を BoundedVerified(k) と
//     明示し、完全とは主張しない。exact 生成系 (Schreier 法) は将来の拡張。
// 座標系: 内部正本は primitive 座標 (KSubgroupFinder.PointGroupData の A/T0 と直結)。conventional 表現は
// 派生値 (中心化格子では有理行列になるため int では持てない、codex R9)。
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;

namespace Crystallography;

/// <summary>N_Aff(G) の純並進部分 K = {a mod Z³ | (I−R)a ∈ Z³ ∀R ∈ P_G} (primitive 座標)。
/// K ≅ (R/Z)^(連続次元) × Π Z/d_i (離散)。連続方向は G を点ごとに centralize するため部分群の
/// 類別には作用しない (polar 方向のメタデータ)。260709Cl 追加。</summary>
public sealed class TranslationNormalizerKernel
{
    /// <summary>離散生成元 (primitive 座標、mod 1 代表)。位数は対応する <see cref="DiscreteInvariantFactors"/>。</summary>
    public Fraction[][] DiscreteGenerators { get; init; }
    /// <summary>離散生成元の位数 d_i (>1 のもののみ)。</summary>
    public int[] DiscreteInvariantFactors { get; init; }
    /// <summary>連続並進方向の基底 (primitive 座標)。polar 群でのみ非空。</summary>
    public Fraction[][] ContinuousBasis { get; init; }
}

/// <summary>N_Aff(G) の (G と純並進核を法として) 非自明な線形部を持つ生成元 1 つ。260709Cl 追加。</summary>
public sealed class NormalizerGenerator
{
    /// <summary>線形部 U ∈ GL(3,Z) (primitive 基底、row-major 9 要素)。</summary>
    public int[] LinearPrimitive { get; init; }
    /// <summary>particular shift a (primitive 座標、mod 1 代表)。他の有効 shift は a + K で尽くされる。</summary>
    public Fraction[] ShiftPrimitive { get; init; }
    /// <summary>det(U) の符号。full affine (±) を保持し、proper (向き保存) への制限は利用側で行う (codex R9)。</summary>
    public int DetSign { get; init; }
}

/// <summary>normalizer 生成集合の完全性の主張レベル。有界探索 (SmallUnimodular(k)) は
/// 「その範囲で見つかったものが全て」であり、N_Aff の生成系である一般保証は無い (codex R9)。260709Cl 追加。</summary>
public enum NormalizerCompleteness
{
    /// <summary>成分 [-1,1] の unimodular 行列の範囲で検証済み。</summary>
    BoundedVerified1 = 1,
    /// <summary>成分 [-2,2] まで拡張して検証済み。</summary>
    BoundedVerified2 = 2,
    /// <summary>厳密な生成系 (将来の拡張: exact linear-normalizer + extension-class stabilizer)。</summary>
    Exact = 100,
}

/// <summary>空間群設定 1 つの affine normalizer データ。260709Cl 追加。</summary>
public sealed class AffineNormalizerData
{
    public int SeriesNumber { get; init; }
    /// <summary>conventional → primitive の基底行列 B (行列 P: (a′b′c′)=(abc)B、KSubgroupFinder.GetPrimitiveBasis)。
    /// LinearPrimitive の conventional 表現は B·U·B⁻¹ (有理) で得る。</summary>
    public Fraction[] PrimitiveBasis { get; init; }
    /// <summary>非自明線形部の生成元 (P_G 左剰余で重複除去済み)。空 = 線形部は P_G で尽きる (例: Pm-3m)。</summary>
    public NormalizerGenerator[] Generators { get; init; }
    public TranslationNormalizerKernel TranslationKernel { get; init; }
    public NormalizerCompleteness Completeness { get; init; }
}

/// <summary>affine normalizer N_Aff(G) の計算 (260709Cl 追加、Phase N1)。結果は series ごとにキャッシュ。
/// ⚠ キャッシュは IT 番号でなく series 単位 — 単斜設定・origin choice・R hex/rhombo で行列は変わる (codex R9)。</summary>
public static class NormalizerFinder
{
    private static readonly ConcurrentDictionary<int, AffineNormalizerData> _cache = new();

    public static AffineNormalizerData Get(int seriesNumber) => _cache.GetOrAdd(seriesNumber, Compute);

    private static AffineNormalizerData Compute(int sn)
    {
        var pg = KSubgroupFinder.BuildPointGroupData(sn);
        int m = pg.LinKeys.Length;

        // ---- 純並進核 K: 斉次系 stack(I − A_i)·a ≡ 0 (mod Z³) ----
        var (rows, _) = BuildCongruenceRows(pg, null, IdentityPermutation(m));
        var zero = new Fraction[3 * m];
        for (int i = 0; i < zero.Length; i++) zero[i] = Fraction.Zero;
        var kernelSol = SolveCongruence(rows, zero);
        if (kernelSol == null)
            throw new InvalidOperationException($"translation kernel must be solvable (sn={sn})"); // 斉次系は a=0 が常に解
        var kernel = new TranslationNormalizerKernel
        {
            DiscreteGenerators = kernelSol.DiscreteKernelGenerators,
            DiscreteInvariantFactors = kernelSol.DiscreteKernelOrders,
            ContinuousBasis = kernelSol.ContinuousKernelBasis,
        };

        // ---- 非自明線形部の有界探索: SmallUnimodular(k) ∩ {U | U·{A}·U⁻¹ = {A}} → shift 解決 ----
        var gens = new List<NormalizerGenerator>();
        var cosetSeen = new List<int[]>(); // P_G 左剰余 (A_j·U) の正準代表で重複除去

        // 探索 1 バウンド分。戻り値 = 「点群は正規化するが lift 不能 (extension class を保存する shift 無し)」だった剰余の数。
        int ScanBound(int k)
        {
            int liftRejected = 0;
            foreach (var u in KSubgroupFinder.SmallUnimodular(k))
            {
                // U が P_G の元そのもの (⇒ G·T を法として自明) はスキップ
                if (pg.A.Any(a => KSubgroupFinder.SameIntVec(a, u)))
                    continue;
                var perm = FindConjugationPermutation(pg, u);
                if (perm == null)
                    continue; // U は点群を正規化しない
                // P_G 左剰余の重複除去: {A_j·U} の中の辞書順最小を代表に
                var rep = CosetCanonical(pg, u);
                if (cosetSeen.Any(r => KSubgroupFinder.SameIntVec(r, rep)))
                    continue;
                var (rows2, rhs) = BuildCongruenceRows(pg, u, perm);
                var sol = SolveCongruence(rows2, rhs);
                if (sol == null)
                {
                    liftRejected++; // nonsymmorphic の lift 制約が生成集合を痩せさせるシグナル (下記 k=2 拡張の判定に使う)
                    continue;
                }
                cosetSeen.Add(rep);
                int det = u[0] * (u[4] * u[8] - u[5] * u[7]) - u[1] * (u[3] * u[8] - u[5] * u[6]) + u[2] * (u[3] * u[7] - u[4] * u[6]);
                gens.Add(new NormalizerGenerator { LinearPrimitive = u, ShiftPrimitive = sol.Particular, DetSign = Math.Sign(det) });
            }
            return liftRejected;
        }

        // 260709Cl (codex R10): k=1 で「正規化するが lift 不能」な U が 1 つでもあった群は k=2 へ拡張する。
        // 実証された失敗モード — P2₁/c では lift 条件が q ≡ 0 (mod 2) を課すため、成分 [-1,1] の lift 可能元
        // (q=0) の積からは q=2 の shear U₀=[[1,0,2],[0,1,0],[0,0,1]] (shift 0 で lift 可能、ITA A1 の
        // IIc 2 軌道化に必須) が生成できず、index3 が 3 軌道 (正 2)・index2 が 2 軌道 (正 1) に割れた。
        // symmorphic 群 (lift 拒否 0) は拡張不要。nonsymmorphic 群は点群正規化フィルタが効いて
        // SmallUnimodular(2) (~10⁵ 元) でも候補が激減するため、コストは許容範囲。
        var completeness = NormalizerCompleteness.BoundedVerified1;
        if (ScanBound(1) > 0)
        {
            ScanBound(2); // k=1 の元は cosetSeen 済みで自然にスキップされる
            completeness = NormalizerCompleteness.BoundedVerified2;
        }

        return new AffineNormalizerData
        {
            SeriesNumber = sn,
            PrimitiveBasis = KSubgroupFinder.GetPrimitiveBasis(sn),
            Generators = [.. gens],
            TranslationKernel = kernel,
            Completeness = completeness,
        };
    }

    #region 共役置換・剰余代表
    private static int[] IdentityPermutation(int m)
    {
        var p = new int[m];
        for (int i = 0; i < m; i++) p[i] = i;
        return p;
    }

    /// <summary>U·A_i·U⁻¹ = A_{π(i)} となる置換 π を求める (無ければ null = U は P_G を正規化しない)。
    /// 260709Cl: Phase 2 (KSubgroupFinder.GetNormalizerOrbits の軌道作用) からも使うため private → internal。</summary>
    internal static int[] FindConjugationPermutation(KSubgroupFinder.PointGroupData pg, int[] u)
    {
        int det = u[0] * (u[4] * u[8] - u[5] * u[7]) - u[1] * (u[3] * u[8] - u[5] * u[6]) + u[2] * (u[3] * u[7] - u[4] * u[6]);
        var uInv = KSubgroupFinder.AdjTimesDet(u, det);
        int m = pg.A.Length;
        var perm = new int[m];
        for (int i = 0; i < m; i++)
        {
            var conj = KSubgroupFinder.MatMulInt(KSubgroupFinder.MatMulInt(u, pg.A[i]), uInv);
            int j = Array.FindIndex(pg.A, a => KSubgroupFinder.SameIntVec(a, conj));
            if (j < 0) return null;
            perm[i] = j;
        }
        return perm;
    }

    /// <summary>P_G 左剰余 {A_j·U | j} の辞書順最小行列 (剰余の正準代表)。</summary>
    private static int[] CosetCanonical(KSubgroupFinder.PointGroupData pg, int[] u)
    {
        int[] best = null;
        foreach (var a in pg.A)
        {
            var cand = KSubgroupFinder.MatMulInt(a, u);
            if (best == null || CompareIntVec(cand, best) < 0) best = cand;
        }
        return best;
    }

    private static int CompareIntVec(int[] a, int[] b)
    {
        for (int i = 0; i < a.Length; i++)
            if (a[i] != b[i]) return a[i].CompareTo(b[i]);
        return 0;
    }
    #endregion

    #region 連立合同式 (SNF ソルバ)
    /// <summary>合同式系 stack(I − A_{π(i)})·a ≡ t_{π(i)} − U·t_i (mod Z³) の係数行列 (3m×3, long) と右辺を組む。
    /// u = null は斉次系 (核の計算、右辺 0 用に rhs は捨てる)。</summary>
    private static (long[,] Rows, Fraction[] Rhs) BuildCongruenceRows(KSubgroupFinder.PointGroupData pg, int[] u, int[] perm)
    {
        int m = pg.A.Length;
        var rows = new long[3 * m, 3];
        var rhs = new Fraction[3 * m];
        for (int i = 0; i < m; i++)
        {
            var aPi = pg.A[perm[i]];
            for (int r = 0; r < 3; r++)
                for (int c = 0; c < 3; c++)
                    rows[3 * i + r, c] = (r == c ? 1 : 0) - aPi[3 * r + c];
            if (u == null)
            {
                rhs[3 * i] = Fraction.Zero; rhs[3 * i + 1] = Fraction.Zero; rhs[3 * i + 2] = Fraction.Zero;
            }
            else
            {
                var tPi = pg.T0[perm[i]];
                var ti = pg.T0[i];
                for (int r = 0; r < 3; r++)
                {
                    Fraction ut = new Fraction(u[3 * r]) * ti[0] + new Fraction(u[3 * r + 1]) * ti[1] + new Fraction(u[3 * r + 2]) * ti[2];
                    rhs[3 * i + r] = tPi[r] - ut;
                }
            }
        }
        return (rows, rhs);
    }

    private sealed class CongruenceSolution
    {
        public Fraction[] Particular;             // mod 1 代表
        public Fraction[][] DiscreteKernelGenerators; // 位数 >1 の離散核生成元 (mod 1)
        public int[] DiscreteKernelOrders;
        public Fraction[][] ContinuousKernelBasis;    // 連続方向
    }

    /// <summary>M·a ≡ b (mod Z^rows) を Smith 標準形で解く。不能なら null。
    /// 行変形は (M | b) に、列変形は M と追跡行列 V に適用し、対角化後
    ///   d_j·y_j ≡ c_j (d_j≠0: 常に可解、y_j = c_j/d_j、核 (1/d_j)Z/Z ≅ Z/d_j)
    ///   0·y_j ≡ c_j  (c_j ∈ Z が可解条件。j&lt;3 なら y_j 自由 = 連続方向)
    /// a = V·y。成分は I−R 由来で小さく (±2)、行数 ≤ 144 のため long で桁あふれしない。260709Cl 追加。</summary>
    private static CongruenceSolution SolveCongruence(long[,] mIn, Fraction[] bIn)
    {
        int rows = mIn.GetLength(0);
        var M = (long[,])mIn.Clone();
        var b = (Fraction[])bIn.Clone();
        // V: 3×3 列変換の追跡 (a = V·y)
        var V = new long[3, 3];
        for (int i = 0; i < 3; i++) V[i, i] = 1;

        void SwapRows(int r1, int r2)
        {
            for (int c = 0; c < 3; c++) (M[r1, c], M[r2, c]) = (M[r2, c], M[r1, c]);
            (b[r1], b[r2]) = (b[r2], b[r1]);
        }
        void AddRow(int dst, int src, long f) // row_dst += f·row_src
        {
            for (int c = 0; c < 3; c++) M[dst, c] += f * M[src, c];
            b[dst] += new Fraction(f) * b[src];
        }
        void SwapCols(int c1, int c2)
        {
            for (int r = 0; r < rows; r++) (M[r, c1], M[r, c2]) = (M[r, c2], M[r, c1]);
            for (int r = 0; r < 3; r++) (V[r, c1], V[r, c2]) = (V[r, c2], V[r, c1]);
        }
        void AddCol(int dst, int src, long f) // col_dst += f·col_src
        {
            for (int r = 0; r < rows; r++) M[r, dst] += f * M[r, src];
            for (int r = 0; r < 3; r++) V[r, dst] += f * V[r, src];
        }
        void NegateCol(int c)
        {
            for (int r = 0; r < rows; r++) M[r, c] = -M[r, c];
            for (int r = 0; r < 3; r++) V[r, c] = -V[r, c];
        }

        // SNF (対角化のみ。d_1 | d_2 | d_3 の整除連鎖は解法に不要なので省略)
        for (int k = 0; k < 3; k++)
        {
            while (true)
            {
                // ピボット選択: 残り小行列で絶対値最小の非ゼロ
                long best = 0; int br = -1, bc = -1;
                for (int r = k; r < rows; r++)
                    for (int c = k; c < 3; c++)
                        if (M[r, c] != 0 && (best == 0 || Math.Abs(M[r, c]) < Math.Abs(best)))
                        { best = M[r, c]; br = r; bc = c; }
                if (br < 0) break; // 残りは全零
                if (br != k) SwapRows(k, br);
                if (bc != k) SwapCols(k, bc);
                if (M[k, k] < 0) NegateCol(k);

                bool reduced = true;
                for (int r = k + 1; r < rows; r++)
                    if (M[r, k] != 0)
                    {
                        AddRow(r, k, -Div(M[r, k], M[k, k]));
                        if (M[r, k] != 0) reduced = false; // 剰余が残った → 次周でより小さいピボット
                    }
                for (int c = k + 1; c < 3; c++)
                    if (M[k, c] != 0)
                    {
                        AddCol(c, k, -Div(M[k, c], M[k, k]));
                        if (M[k, c] != 0) reduced = false;
                    }
                if (reduced) break;
            }
        }

        // 可解性と解の構成
        var y = new Fraction[3];
        var discGens = new List<Fraction[]>();
        var discOrders = new List<int>();
        var contBasis = new List<Fraction[]>();
        for (int j = 0; j < 3; j++)
        {
            long d = M[j, j];
            if (d != 0)
            {
                y[j] = b[j] / new Fraction(d);
                if (d > 1 || d < -1)
                {
                    long dAbs = Math.Abs(d);
                    var g = MulV(V, [new Fraction(1, dAbs), Fraction.Zero, Fraction.Zero], j);
                    discGens.Add(RationalMatrix.ModVec1(g));
                    discOrders.Add((int)dAbs);
                }
            }
            else
            {
                if (!b[j].IsInteger) return null; // 0·y_j ≡ 非整数 → 不能
                y[j] = Fraction.Zero;
                contBasis.Add(MulV(V, [Fraction.One, Fraction.Zero, Fraction.Zero], j));
            }
        }
        // 対角より下の行 (rank 超過分): 0 ≡ c_j が可解条件
        for (int r = 3; r < rows; r++)
            if (!b[r].IsInteger) return null;

        var particular = RationalMatrix.ModVec1(MulVFull(V, y));
        return new CongruenceSolution
        {
            Particular = particular,
            DiscreteKernelGenerators = [.. discGens],
            DiscreteKernelOrders = [.. discOrders],
            ContinuousKernelBasis = [.. contBasis],
        };
    }

    /// <summary>floor 除算 (C# の / はゼロ方向切り捨てなので負値で剰余が負にならないように)。</summary>
    private static long Div(long a, long b) => (long)Math.Floor((double)a / b);

    /// <summary>V の第 j 列方向の寄与: V·(s·e_j) (s = v[0] に格納された係数)。</summary>
    private static Fraction[] MulV(long[,] V, Fraction[] v, int j)
        => [new Fraction(V[0, j]) * v[0], new Fraction(V[1, j]) * v[0], new Fraction(V[2, j]) * v[0]];

    private static Fraction[] MulVFull(long[,] V, Fraction[] y)
        =>
        [
            new Fraction(V[0, 0]) * y[0] + new Fraction(V[0, 1]) * y[1] + new Fraction(V[0, 2]) * y[2],
            new Fraction(V[1, 0]) * y[0] + new Fraction(V[1, 1]) * y[1] + new Fraction(V[1, 2]) * y[2],
            new Fraction(V[2, 0]) * y[0] + new Fraction(V[2, 1]) * y[1] + new Fraction(V[2, 2]) * y[2],
        ];
    #endregion
}
