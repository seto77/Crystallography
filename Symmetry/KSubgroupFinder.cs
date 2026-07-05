// 260705Cl 新規: klassengleiche (k-) 部分群の実行時計算エンジン (Phase 2c)。
// TSubgroupFinder (t-) とは別クラス (計画書 §4.2: 列挙対象が異なるため共有しない。
// 型同定・操作集合完全一致検証・格子ユーティリティは今後の共通化候補)。
//
// 現時点では Step 1 (部分格子 T′ の列挙) のみ実装。complement 列挙 (Step 2) 以降は後続コミットで追加する。
// 詳細ロードマップ: .project-guidance/ReciPro_FormGroupRelations改修計画.md §4。
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;

namespace Crystallography;

/// <summary>k-最大部分群の列挙 (実装中、Phase 2c)。</summary>
public static class KSubgroupFinder
{
    #region 実並進格子の primitive 基底
    /// <summary>親空間群 (通し番号) の実並進格子 T の primitive 基底を、慣用胞座標系で返す (row-major 3 行)。
    /// 中心化ベクトル (F/I/A/B/C/R 等) を恒等線形部の操作から実データとして集め、Z³ と合わせた生成元集合から
    /// 整数格子基底を抽出する (codex 提案: 手書き分岐でなく、実データ駆動で中心化の型を問わず統一的に扱う)。
    /// 260705Cl 追加。</summary>
    public static Fraction[] GetPrimitiveBasis(int seriesNumber)
    {
        var ops = TSubgroupFinder.GetExpandedOps(seriesNumber);
        var centering = ops.Where(o => IsIdentityLinear(o))
                            .Select(o => { var t = o.SeitzTranslation; return new[] { Fraction.FromDouble(t.U), Fraction.FromDouble(t.V), Fraction.FromDouble(t.W) }; })
                            .Where(c => !(c[0].IsZero && c[1].IsZero && c[2].IsZero))
                            .ToList();
        if (centering.Count == 0)
            return RationalMatrix.FromInt([1, 0, 0, 0, 1, 0, 0, 0, 1]); // P: 慣用胞 = primitive

        // 共通分母 d で整数化し、d·Z³ と合わせて整数格子基底を抽出、最後に d で割って戻す。
        BigInteger dBig = 1;
        foreach (var c in centering)
            foreach (var f in c)
                dBig = Lcm(dBig, f.Den);
        long d = (long)dBig;

        List<long[]> generators = [[d, 0, 0], [0, d, 0], [0, 0, d]];
        foreach (var c in centering)
            generators.Add([(long)(c[0].Num * (d / c[0].Den)), (long)(c[1].Num * (d / c[1].Den)), (long)(c[2].Num * (d / c[2].Den))]);

        var basisInt = IntegerLattice.BasisFromGenerators(generators);
        if (basisInt == null)
            throw new InvalidOperationException($"primitive basis extraction failed (sn={seriesNumber}): generators did not span rank 3");

        // 260705Cl 修正: 既存コードベースの規約 ((a′,b′,c′)=(a,b,c)·P、P の列=新基底) に合わせ、
        // 返す基底行列は「列=基底ベクトル」にする (IntegerLattice.BasisFromGenerators は行=基底ベクトルで
        // 返すため転記時に添字を入れ替える)。これを row=基底ベクトルのまま返すと、FilterPointGroupInvariant
        // の A_R=B⁻¹RB (R は列ベクトル左作用) が数学的に誤った式になる実バグがあった。
        var basis = new Fraction[9];
        for (int i = 0; i < 3; i++)      // i = 成分 (行)
            for (int j = 0; j < 3; j++)  // j = 基底ベクトル番号 (列)
                basis[i * 3 + j] = new Fraction(basisInt[j][i], d);
        return basis;
    }

    private static bool IsIdentityLinear(in SymmetryOperation op)
    {
        var m = SeitzNotation.LinearMatrix(op);
        return m[0, 0] == 1 && m[1, 1] == 1 && m[2, 2] == 1 &&
               m[0, 1] == 0 && m[0, 2] == 0 && m[1, 0] == 0 && m[1, 2] == 0 && m[2, 0] == 0 && m[2, 1] == 0;
    }

    private static BigInteger Lcm(BigInteger a, BigInteger b) => a / BigInteger.GreatestCommonDivisor(a, b) * b;
    #endregion

    #region 部分格子 T′ の列挙 (HNF)
    /// <summary>指数 n の Z³ 部分格子を表す HNF (row-major 9 要素、列ベクトル規約: 各列が部分格子の基底ベクトル、
    /// R が列ベクトルへ左作用する既存規約と整合)。H = [[a,0,0],[x,b,0],[y,z,c]] (列は (a,x,y),(0,b,z),(0,0,c))、
    /// abc=n、0≤x&lt;b、0≤y,z&lt;c で完全列挙する (n=2: 7 個 / 3: 13 個 / 4: 35 個)。
    /// 260705Cl 修正: 初版は codex 提案の行/列を取り違えて転置してしまい (列が (a,0,0),(x,b,0),(y,z,c) になる
    /// 誤った形で実装していた)、既知の I-centering 部分格子 (index4, {(2,0,0),(0,2,0),(1,1,1)} が張る格子)
    /// が列として表現できず点群不変フィルタが常に空を返す実バグがあった。IntegerLattice.BasisFromGenerators
    /// で実際に求まる基底と照合して修正 (下記コメント参照)。260705Cl 追加。</summary>
    public static List<int[]> EnumerateHnf(int index)
    {
        var list = new List<int[]>();
        for (int a = 1; a <= index; a++)
        {
            if (index % a != 0) continue;
            int rest = index / a;
            for (int b = 1; b <= rest; b++)
            {
                if (rest % b != 0) continue;
                int c = rest / b;
                for (int x = 0; x < b; x++)
                    for (int y = 0; y < c; y++)
                        for (int z = 0; z < c; z++)
                            list.Add([a, 0, 0, x, b, 0, y, z, c]);
            }
        }
        return list;
    }

    /// <summary>親空間群の点群 (線形部) すべてで不変な部分格子 T′ (primitive 基底 B 上、HNF H で指定) だけを残す。
    /// 判定: A_R = B⁻¹RB が整数行列であること (実並進格子が R で不変なのは自明、assert 用) を確認したうえで、
    /// A_R·H の各列が H の列の整数combinationで表せるか (= H⁻¹·A_R·H が整数行列か) を確認する。260705Cl 追加。</summary>
    public static List<int[]> FilterPointGroupInvariant(int seriesNumber, IEnumerable<int[]> candidates)
    {
        var ops = TSubgroupFinder.GetExpandedOps(seriesNumber);
        var linKeys = new List<int[]>();
        foreach (var op in ops)
        {
            var key = LinKeyOf(op);
            if (!linKeys.Any(k => SameIntVec(k, key))) linKeys.Add(key);
        }
        var basis = GetPrimitiveBasis(seriesNumber);
        var basisInv = RationalMatrix.Invert3(basis);
        if (basisInv == null) throw new InvalidOperationException("primitive basis is singular");

        var result = new List<int[]>();
        foreach (var h in candidates)
        {
            bool invariantForAll = true;
            foreach (var lin in linKeys)
            {
                var rRational = RationalMatrix.FromInt(lin);
                var aR = RationalMatrix.Mul(RationalMatrix.Mul(basisInv, rRational), basis);
                var aRInt = RationalMatrix.ToIntOrNull(aR);
                if (aRInt == null)
                    throw new InvalidOperationException("B^-1 R B is not integral — primitive basis is wrong");

                if (!IsLatticeInvariant(aRInt, h)) { invariantForAll = false; break; }
            }
            if (invariantForAll) result.Add(h);
        }
        return result;
    }

    /// <summary>整数行列 aR による HNF h の像が h と同じ部分格子に留まるか (aR·h の各列が h の整数列combinationか) を、
    /// 有理逆行列を使わず整数演算のクラメル解 (adj/det) で判定する (codex 提案: 高速・頑丈)。</summary>
    private static bool IsLatticeInvariant(int[] aR, int[] h)
    {
        // aR * h (3x3 整数行列積)
        var m = new int[9];
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                m[i * 3 + j] = aR[i * 3] * h[j] + aR[i * 3 + 1] * h[3 + j] + aR[i * 3 + 2] * h[6 + j];

        int detH = h[0] * (h[4] * h[8] - h[5] * h[7]) - h[1] * (h[3] * h[8] - h[5] * h[6]) + h[2] * (h[3] * h[7] - h[4] * h[6]);
        // adj(h) (余因子転置): h * adj(h) = det(h) * I
        int[] adj =
        [
            h[4] * h[8] - h[5] * h[7], h[2] * h[7] - h[1] * h[8], h[1] * h[5] - h[2] * h[4],
            h[5] * h[6] - h[3] * h[8], h[0] * h[8] - h[2] * h[6], h[2] * h[3] - h[0] * h[5],
            h[3] * h[7] - h[4] * h[6], h[1] * h[6] - h[0] * h[7], h[0] * h[4] - h[1] * h[3],
        ];
        // h^-1 * m = adj(h) * m / det(h) が整数行列かどうか
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                int num = adj[i * 3] * m[j] + adj[i * 3 + 1] * m[3 + j] + adj[i * 3 + 2] * m[6 + j];
                if (num % detH != 0) return false;
            }
        return true;
    }

    private static int[] LinKeyOf(in SymmetryOperation op)
    {
        var m = SeitzNotation.LinearMatrix(op);
        return [m[0, 0], m[0, 1], m[0, 2], m[1, 0], m[1, 1], m[1, 2], m[2, 0], m[2, 1], m[2, 2]];
    }

    private static bool SameIntVec(int[] a, int[] b)
    {
        for (int i = 0; i < a.Length; i++) if (a[i] != b[i]) return false;
        return true;
    }
    #endregion

    #region complement 列挙 (Step 2, Q=G/T′ 内の section)
    // 260705Cl 追加 (Phase 2c Step2)。設計は codex との3回目の相談で確定
    // (.project-guidance/ReciPro_FormGroupRelations改修計画.md §4.1 item2)。
    //
    // 数学的な骨格: T′ (index n) を法とした有限商 Q=G/T′ (|Q|=m·n、m=|P_G|) の中で、
    // 並進部分群 T/T′ と交わらず点群へ全射する complement (H の実体、位数 m) を全数列挙する。
    // 各線形部 R_i の実並進 t_i (親の実データ由来、primitive 座標) は mod T′ で n 通りの coset の
    // どれかに属する。空間群の合成則 (R_i,t_i)(R_j,t_j)=(R_iR_j, t_i+R_i·t_j) から、
    // cocycle f[i,j] = t_i^0 + A_i·t_j^0 − t_{mul[i,j]}^0 (実データの group law が保証するので
    // 必ず整数、primitive 座標) が直接計算でき、これと coset オフセットの組合せで Q の乗積表が
    // 具体的に (代表ベクトルの加算だけで) 書ける。
    //
    // 生成集合は点群の小さい生成元を1組だけ固定して構わない (codex 確認済み): 任意の complement C は
    // 射影 C→P_G が全単射なので、生成元の像 (lift) で一意に決まる。したがって、固定した生成元への
    // n^k 通りの coset 割当を全数試せば、必ずどの complement もそのうちの1通りと一致し、取り逃さない。
    // 複数の生成集合を試す必要はない (dedupe の保険にしかならない)。
    //
    // 恒等線形部の基準並進 t_e^0 は必ず 0 に固定する。GetExpandedOps の並び順で中心化コピーが
    // 最初に来ると Q の単位元が (E,0) にならず cocycle 計算全体が破綻するため。

    /// <summary>点群 (線形部) の乗積表・逆元・primitive 座標での A_R=B⁻¹RB・基準並進 t_i^0 をまとめた
    /// 内部データ (T′ に依存しない部分、Step2 で繰り返し使う)。260705Cl 追加。</summary>
    private sealed class PointGroupData
    {
        public int[][] LinKeys { get; init; }
        public int[,] Mul { get; init; }
        public int[] Inv { get; init; }
        public int E { get; init; }
        public int[][] A { get; init; }
        public Fraction[][] T0 { get; init; }
    }

    /// <summary>指定 T′ (HNF hnf) に対する有限商 Q=G/T′ の代数データ (coset 代表・正準ラベル・乗積表)。
    /// 260705Cl 追加。</summary>
    private sealed class QuotientData
    {
        public int N { get; init; }
        public long[][] Reps { get; init; }
        public Fraction[][] Labels { get; init; }
        /// <summary>フラット化した乗積表 [qa*(m*N)+qb]。qa/qb = 線形部index*N + coset index。</summary>
        public int[] MulQ { get; init; }
    }

    private static PointGroupData BuildPointGroupData(int sn)
    {
        var ops = TSubgroupFinder.GetExpandedOps(sn);
        var basis = GetPrimitiveBasis(sn);
        var basisInv = RationalMatrix.Invert3(basis) ?? throw new InvalidOperationException("primitive basis is singular");

        var linKeys = new List<int[]>();
        var t0List = new List<Fraction[]>();
        foreach (var op in ops)
        {
            var key = LinKeyOf(op);
            if (linKeys.Any(k => SameIntVec(k, key))) continue;
            linKeys.Add(key);
            var t = op.SeitzTranslation;
            var tConv = new Fraction[] { Fraction.FromDouble(t.U), Fraction.FromDouble(t.V), Fraction.FromDouble(t.W) };
            t0List.Add(RationalMatrix.MulVec(basisInv, tConv));
        }
        int m = linKeys.Count;
        int e = Enumerable.Range(0, m).First(i => IsIdentityLin(linKeys[i]));
        // 260705Cl: 恒等線形部の基準並進は必ず 0 に固定 (上記コメント参照)。
        t0List[e] = [Fraction.Zero, Fraction.Zero, Fraction.Zero];

        var mul = new int[m, m];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < m; j++)
                mul[i, j] = FindIntKey(linKeys, MatMulInt(linKeys[i], linKeys[j]));

        var inv = new int[m];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < m; j++)
                if (mul[i, j] == e) { inv[i] = j; break; }

        var A = new int[m][];
        for (int i = 0; i < m; i++)
        {
            var af = RationalMatrix.Mul(RationalMatrix.Mul(basisInv, RationalMatrix.FromInt(linKeys[i])), basis);
            A[i] = RationalMatrix.ToIntOrNull(af) ?? throw new InvalidOperationException($"B^-1 R B is not integral (sn={sn}, lin={i})");
        }

        return new PointGroupData { LinKeys = [.. linKeys], Mul = mul, Inv = inv, E = e, A = A, T0 = [.. t0List] };
    }

    private static QuotientData BuildQuotient(PointGroupData pg, int[] hnf)
    {
        int m = pg.LinKeys.Length;
        var hInv = RationalMatrix.Invert3(RationalMatrix.FromInt(hnf)) ?? throw new InvalidOperationException("HNF is singular");
        int a = hnf[0], b = hnf[4], c = hnf[8];
        int n = a * b * c;

        // T′ の T における coset 代表: 標準基底の (i,j,k), 0≤i<a,j<b,k<c がちょうど n 個の代表系になる
        // (列 HNF が (a,x,y),(0,b,z),(0,0,c) の形であることに由来する代数的事実、codex 確認済み)。
        var reps = new List<long[]>();
        var labels = new List<Fraction[]>();
        for (int i = 0; i < a; i++)
            for (int j = 0; j < b; j++)
                for (int k = 0; k < c; k++)
                {
                    reps.Add([i, j, k]);
                    labels.Add(RationalMatrix.ModVec1(RationalMatrix.MulVec(hInv, [new Fraction(i), new Fraction(j), new Fraction(k)])));
                }

        int CosetIndexOf(long[] v)
        {
            var lab = RationalMatrix.ModVec1(RationalMatrix.MulVec(hInv, [new Fraction(v[0]), new Fraction(v[1]), new Fraction(v[2])]));
            for (int t = 0; t < labels.Count; t++)
                if (RationalMatrix.VecEquals(labels[t], lab)) return t;
            throw new InvalidOperationException("coset representative not found (HNF label construction bug)");
        }

        // cocycle f[i,j] = t0[i] + A_i·t0[j] − t0[mul[i,j]] (primitive 座標)。
        // 実データ (親空間群の実並進) の group law から必ず整数になる — 整数でなければ基準並進・A 行列の規約が壊れている。
        var f = new int[m * m][];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < m; j++)
            {
                var aiT0j = RationalMatrix.MulVec(RationalMatrix.FromInt(pg.A[i]), pg.T0[j]);
                var sum = RationalMatrix.SubVec(RationalMatrix.AddVec(pg.T0[i], aiT0j), pg.T0[pg.Mul[i, j]]);
                var fij = new int[3];
                for (int c2 = 0; c2 < 3; c2++)
                {
                    if (!sum[c2].IsInteger)
                        throw new InvalidOperationException($"cocycle f[{i},{j}] is not integer ({sum[c2]}) — base translation convention bug");
                    fij[c2] = (int)sum[c2].Num;
                }
                f[i * m + j] = fij;
            }

        var mulQ = new int[m * n * m * n];
        for (int i = 0; i < m; i++)
            for (int aRep = 0; aRep < n; aRep++)
                for (int j = 0; j < m; j++)
                    for (int bRep = 0; bRep < n; bRep++)
                    {
                        var aj = MatMulIntVec(pg.A[i], reps[bRep]);
                        var fij = f[i * m + j];
                        var target = new long[] { fij[0] + reps[aRep][0] + aj[0], fij[1] + reps[aRep][1] + aj[1], fij[2] + reps[aRep][2] + aj[2] };
                        int qa = i * n + aRep, qb = j * n + bRep;
                        mulQ[qa * (m * n) + qb] = pg.Mul[i, j] * n + CosetIndexOf(target);
                    }

        return new QuotientData { N = n, Reps = [.. reps], Labels = [.. labels], MulQ = mulQ };
    }

    /// <summary>親空間群 sn の点群を保つ部分格子 T′ (HNF hnf) について、Q=G/T′ 内の complement を全数列挙する。
    /// 返り値は各 complement を表す σ (m 個、線形部 index → T′ の coset index)。260705Cl 追加。</summary>
    public static List<int[]> EnumerateComplements(int seriesNumber, int[] hnf)
    {
        var pg = BuildPointGroupData(seriesNumber);
        var q = BuildQuotient(pg, hnf);
        int m = pg.LinKeys.Length, n = q.N, size = m * n;

        var gens = ChooseGenerators(pg);
        int k = gens.Count;
        long total = 1;
        for (int t = 0; t < k; t++) total *= n;

        var found = new List<int[]>();
        var seen = new HashSet<string>();
        var offsets = new int[k];
        for (long combo = 0; combo < total; combo++)
        {
            long rem = combo;
            for (int t = 0; t < k; t++) { offsets[t] = (int)(rem % n); rem /= n; }

            var seed = new List<int> { pg.E * n };
            for (int t = 0; t < k; t++) seed.Add(gens[t] * n + offsets[t]);

            var closure = ClosureQ(seed, q.MulQ, size, m);
            if (closure.Count != m) continue;

            var sigma = new int[m];
            var seenLin = new bool[m];
            bool ok = true;
            foreach (var elem in closure)
            {
                int lin = elem / n, c = elem % n;
                if (seenLin[lin]) { ok = false; break; }
                seenLin[lin] = true;
                sigma[lin] = c;
            }
            if (!ok || seenLin.Any(x => !x)) continue; // 射影が全単射でない (理論上起きないはずの防御的チェック)

            var key = string.Join(",", sigma);
            if (seen.Add(key)) found.Add(sigma);
        }
        return found;
    }

    /// <summary>見つかった complement (σ 配列) を Q 内の共役 (q·H·q⁻¹, q は Q の全元) で類別する。
    /// 計画書 §4.4「複数 complement の分類」対応。260705Cl 追加。</summary>
    public static List<List<int[]>> GroupComplementsByConjugacy(int seriesNumber, int[] hnf, List<int[]> complements)
    {
        var pg = BuildPointGroupData(seriesNumber);
        var q = BuildQuotient(pg, hnf);
        int m = pg.LinKeys.Length, n = q.N, size = m * n;
        int qIdentity = pg.E * n;

        var invQ = new int[size];
        for (int aElem = 0; aElem < size; aElem++)
            for (int bElem = 0; bElem < size; bElem++)
                if (q.MulQ[aElem * size + bElem] == qIdentity) { invQ[aElem] = bElem; break; }

        int[] Conjugate(int[] sigma, int qElem)
        {
            int qInv = invQ[qElem];
            var result = new int[m];
            for (int i = 0; i < m; i++)
            {
                int h = i * n + sigma[i];
                int step1 = q.MulQ[qElem * size + h];
                int step2 = q.MulQ[step1 * size + qInv];
                result[step2 / n] = step2 % n;
            }
            return result;
        }

        var keys = complements.Select(s => string.Join(",", s)).ToList();
        var classes = new List<List<int[]>>();
        var assigned = new bool[complements.Count];
        for (int i = 0; i < complements.Count; i++)
        {
            if (assigned[i]) continue;
            var cls = new List<int[]> { complements[i] };
            assigned[i] = true;
            for (int qElem = 0; qElem < size; qElem++)
            {
                var key = string.Join(",", Conjugate(complements[i], qElem));
                for (int j = 0; j < complements.Count; j++)
                    if (!assigned[j] && keys[j] == key) { cls.Add(complements[j]); assigned[j] = true; }
            }
            classes.Add(cls);
        }
        return classes;
    }

    private static List<int> ChooseGenerators(PointGroupData pg)
    {
        int m = pg.LinKeys.Length;
        var spanned = ClosureLin([pg.E], pg.Mul);
        var gens = new List<int>();
        for (int i = 0; i < m && spanned.Count < m; i++)
        {
            if (spanned.Contains(i)) continue;
            gens.Add(i);
            spanned = ClosureLin([.. spanned, i], pg.Mul);
        }
        return gens;
    }

    private static HashSet<int> ClosureLin(IEnumerable<int> seed, int[,] mul)
    {
        var s = new HashSet<int>(seed);
        var queue = new Queue<int>(s);
        while (queue.Count > 0)
        {
            int a = queue.Dequeue();
            foreach (var b in s.ToArray())
            {
                if (s.Add(mul[a, b])) queue.Enqueue(mul[a, b]);
                if (s.Add(mul[b, a])) queue.Enqueue(mul[b, a]);
            }
        }
        return s;
    }

    /// <summary>Q の乗積表 mulQ 上で closure を取る。要素数が maxSize を超えたら即座に打ち切る
    /// (complement 候補は必ず位数 m なので、それを超えたら不採用確定。総当たりの高速化)。260705Cl 追加。</summary>
    private static HashSet<int> ClosureQ(IEnumerable<int> seed, int[] mulQ, int size, int maxSize)
    {
        var s = new HashSet<int>(seed);
        var queue = new Queue<int>(s);
        while (queue.Count > 0)
        {
            int a = queue.Dequeue();
            foreach (var b in s.ToArray())
            {
                int ab = mulQ[a * size + b];
                if (s.Add(ab)) { if (s.Count > maxSize) return s; queue.Enqueue(ab); }
                int ba = mulQ[b * size + a];
                if (s.Add(ba)) { if (s.Count > maxSize) return s; queue.Enqueue(ba); }
            }
        }
        return s;
    }

    private static long[] MatMulIntVec(int[] mat, long[] vec)
    {
        var r = new long[3];
        for (int i = 0; i < 3; i++)
            r[i] = mat[i * 3] * vec[0] + mat[i * 3 + 1] * vec[1] + mat[i * 3 + 2] * vec[2];
        return r;
    }

    private static int[] MatMulInt(int[] a, int[] b)
    {
        var c = new int[9];
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                c[i * 3 + j] = a[i * 3] * b[j] + a[i * 3 + 1] * b[3 + j] + a[i * 3 + 2] * b[6 + j];
        return c;
    }

    private static int FindIntKey(List<int[]> list, int[] key)
    {
        for (int i = 0; i < list.Count; i++)
            if (SameIntVec(list[i], key)) return i;
        throw new InvalidOperationException("linear part not found in point-group closure (multiplication table bug)");
    }

    private static bool IsIdentityLin(int[] m)
        => m[0] == 1 && m[4] == 1 && m[8] == 1 && m[1] == 0 && m[2] == 0 && m[3] == 0 && m[5] == 0 && m[6] == 0 && m[7] == 0;
    #endregion

    #region 型同定 (Step 3, IdentifyK)
    // 260705Cl 追加 (Phase 2c Step3)。設計は codex との4回目の相談で確定
    // (.project-guidance/ReciPro_FormGroupRelations改修計画.md §4.1 item3)。
    //
    // P = S·U·C⁻¹ (x_parent = P·x_child + p の規約、既存 TSubgroupFinder と同じ)。
    //   S = B·H   … T′ の primitive 基底 (親慣用胞座標)
    //   C = GetPrimitiveBasis(candSn) … 候補設定の primitive 基底 (候補慣用胞座標)
    //   U … T′ の primitive 基底の取り替えを表す unimodular 整数行列 (det=±1)
    // U は「小さい成分の総当たり」で探す (K=1 で見つからなければ K=2 にフォールバック)。
    // 誤同定は最終的な操作集合完全一致検証 (SolveOriginShiftK) で原理的に排除されるため、
    // U の探索が理論上完全でなくても安全側 (見つからなければ「未同定」、t-エンジンと同じ正直な方針)。
    // フィルタとして、A_H[i]=H⁻¹A[i]H (T′ 自身の座標系での親点群表現) を U で共役したものが、
    // 候補の A_cand[i]=C⁻¹R_i C と整数行列の集合として一致するかを先に見る (安価、格子一致は構成上自動)。

    private static readonly Dictionary<int, List<int[]>> _unimodularCache = [];

    /// <summary>成分が [-k,k] の unimodular (det=±1) 整数 3×3 行列を総当たりで列挙する (キャッシュ済み)。260705Cl 追加。</summary>
    private static List<int[]> SmallUnimodular(int k)
    {
        if (_unimodularCache.TryGetValue(k, out var cached)) return cached;
        var result = new List<int[]>();
        var m = new int[9];
        void Rec(int idx)
        {
            if (idx == 9)
            {
                int det = m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6]) + m[2] * (m[3] * m[7] - m[4] * m[6]);
                if (det == 1 || det == -1) result.Add((int[])m.Clone());
                return;
            }
            for (int v = -k; v <= k; v++) { m[idx] = v; Rec(idx + 1); }
        }
        Rec(0);
        _unimodularCache[k] = result;
        return result;
    }

    /// <summary>候補設定 candSn の点群線形部・primitive 基底・その逆・primitive 座標での線形部・
    /// 恒等線形部の中心化並進 (centering cosets) をまとめたキャッシュ。260705Cl 追加。</summary>
    private sealed class CandidateData
    {
        public int[][] LinKeys { get; init; }
        public Fraction[] CInv { get; init; }
        public int[][] ACand { get; init; }         // C⁻¹ R C (整数、candidate primitive 座標)
        public Fraction[][] Rt { get; init; }        // LinKeys[i] に対応する実際の並進 (candidate 慣用胞座標、代表1つ、無還元)
        public Fraction[][] Centering { get; init; } // 恒等線形部の中心化並進 (candidate 慣用胞座標、mod1、重複無し)
    }

    private static readonly Dictionary<int, CandidateData> _candidateCache = [];
    private static Dictionary<int, List<int>> _candidatesByOrder;

    private static CandidateData BuildCandidateData(int sn)
    {
        if (_candidateCache.TryGetValue(sn, out var cached)) return cached;
        var ops = TSubgroupFinder.GetExpandedOps(sn);
        var c = GetPrimitiveBasis(sn);
        var cInv = RationalMatrix.Invert3(c) ?? throw new InvalidOperationException("candidate primitive basis is singular");

        var linKeys = new List<int[]>();
        var rt = new List<Fraction[]>();
        foreach (var op in ops)
        {
            var key = LinKeyOf(op);
            if (linKeys.Any(k => SameIntVec(k, key))) continue;
            linKeys.Add(key);
            var t = op.SeitzTranslation;
            rt.Add([Fraction.FromDouble(t.U), Fraction.FromDouble(t.V), Fraction.FromDouble(t.W)]);
        }
        int m = linKeys.Count;
        var aCand = new int[m][];
        for (int i = 0; i < m; i++)
        {
            var af = RationalMatrix.Mul(RationalMatrix.Mul(cInv, RationalMatrix.FromInt(linKeys[i])), c);
            aCand[i] = RationalMatrix.ToIntOrNull(af) ?? throw new InvalidOperationException($"candidate C^-1 R C not integral (sn={sn})");
        }

        var centering = new List<Fraction[]>();
        foreach (var op in ops)
        {
            if (!IsIdentityLin(LinKeyOf(op))) continue;
            var t = op.SeitzTranslation;
            var v = RationalMatrix.ModVec1([Fraction.FromDouble(t.U), Fraction.FromDouble(t.V), Fraction.FromDouble(t.W)]);
            if (!centering.Any(x => RationalMatrix.VecEquals(x, v))) centering.Add(v);
        }

        var data = new CandidateData { LinKeys = [.. linKeys], CInv = cInv, ACand = aCand, Rt = [.. rt], Centering = [.. centering] };
        _candidateCache[sn] = data;
        return data;
    }

    /// <summary>点群位数 (相異なる線形部数) ごとに候補設定の通し番号をまとめた索引 (初回のみ全 530 設定を走査)。
    /// U 探索のたびに全設定を試すコストを避けるための絞り込み。260705Cl 追加。</summary>
    private static Dictionary<int, List<int>> CandidatesByOrder()
    {
        if (_candidatesByOrder != null) return _candidatesByOrder;
        var map = new Dictionary<int, List<int>>();
        for (int sn = 1; sn < SymmetryStatic.TotalSpaceGroupNumber; sn++)
        {
            if (SymmetryStatic.Symmetries[sn].SpaceGroupNumber == 0) continue;
            int order = BuildCandidateData(sn).LinKeys.Length;
            if (!map.TryGetValue(order, out var list)) map[order] = list = [];
            list.Add(sn);
        }
        _candidatesByOrder = map;
        return map;
    }

    /// <summary>親空間群 parentSn の complement (T′=hnf, σ=sigma) を型同定する。
    /// 成功時 (childSn, P, p)（x_parent = P·x_child + p、親慣用胞座標、Fraction）、失敗時 (-1, null, null)。
    /// 260705Cl 追加 (Phase 2c Step3)。</summary>
    public static (int Child, Fraction[] P, Fraction[] Shift) IdentifyK(int parentSn, int[] hnf, int[] sigma)
    {
        var pg = BuildPointGroupData(parentSn);
        var q = BuildQuotient(pg, hnf);
        int m = pg.LinKeys.Length;

        var basis = GetPrimitiveBasis(parentSn);
        var hFrac = RationalMatrix.FromInt(hnf);
        var hInv = RationalMatrix.Invert3(hFrac) ?? throw new InvalidOperationException("HNF is singular");
        var s = RationalMatrix.Mul(basis, hFrac); // S = B·H (親慣用胞座標、T′ の primitive 基底)

        // A_H[i] = H⁻¹·A[i]·H (T′ 自身の座標系での親点群の表現。Step1 の point-group-invariant 検証により必ず整数)
        var aH = new int[m][];
        for (int i = 0; i < m; i++)
        {
            var af = RationalMatrix.Mul(RationalMatrix.Mul(hInv, RationalMatrix.FromInt(pg.A[i])), hFrac);
            aH[i] = RationalMatrix.ToIntOrNull(af) ?? throw new InvalidOperationException($"H^-1 A H not integral (sn={parentSn}) — sublattice is not point-group-invariant");
        }

        // H の実際の並進 (親慣用胞座標、無還元): t_parent[i] = B·(t0[i]+reps[σ(i)])
        var tParent = new Fraction[m][];
        for (int i = 0; i < m; i++)
        {
            var repI = q.Reps[sigma[i]];
            var sum = RationalMatrix.AddVec(pg.T0[i], [new Fraction(repI[0]), new Fraction(repI[1]), new Fraction(repI[2])]);
            tParent[i] = RationalMatrix.MulVec(basis, sum);
        }

        var byOrder = CandidatesByOrder();
        if (!byOrder.TryGetValue(m, out var candList)) return (-1, null, null);

        foreach (int k in new[] { 1, 2 })
        {
            foreach (var u in SmallUnimodular(k))
            {
                var uFrac = RationalMatrix.FromInt(u);
                var uInv = RationalMatrix.Invert3(uFrac); // det=±1 なので必ず存在

                var conjugated = new int[m][];
                bool intAll = true;
                for (int i = 0; i < m && intAll; i++)
                {
                    var cf = RationalMatrix.Mul(RationalMatrix.Mul(uInv, RationalMatrix.FromInt(aH[i])), uFrac);
                    var ci = RationalMatrix.ToIntOrNull(cf);
                    if (ci == null) { intAll = false; break; } // U が unimodular なら理論上必ず整数 (防御的チェック)
                    conjugated[i] = ci;
                }
                if (!intAll) continue;

                foreach (var candSn in candList)
                {
                    var cand = BuildCandidateData(candSn);
                    if (!SetEqualsIntMats(conjugated, cand.ACand)) continue;

                    var p = RationalMatrix.Mul(RationalMatrix.Mul(s, uFrac), cand.CInv);
                    var pInv = RationalMatrix.Invert3(p);
                    if (pInv == null) continue;

                    var rChild = new int[m][];
                    var tChild = new Fraction[m][];
                    bool ok = true;
                    for (int i = 0; i < m && ok; i++)
                    {
                        var rf = RationalMatrix.Mul(RationalMatrix.Mul(pInv, RationalMatrix.FromInt(pg.LinKeys[i])), p);
                        var ri = RationalMatrix.ToIntOrNull(rf);
                        if (ri == null) { ok = false; break; }
                        rChild[i] = ri;
                        tChild[i] = RationalMatrix.MulVec(pInv, tParent[i]);
                    }
                    if (!ok) continue;

                    var origin = SolveOriginShiftK(rChild, tChild, cand);
                    if (origin == null) continue;

                    var pShift = RationalMatrix.MulVec(p, origin);
                    return (candSn, p, pShift);
                }
            }
        }
        return (-1, null, null);
    }

    private static bool SetEqualsIntMats(int[][] a, int[][] b)
    {
        if (a.Length != b.Length) return false;
        var used = new bool[b.Length];
        foreach (var m in a)
        {
            int idx = -1;
            for (int j = 0; j < b.Length; j++)
                if (!used[j] && SameIntVec(m, b[j])) { idx = j; break; }
            if (idx < 0) return false;
            used[idx] = true;
        }
        return true;
    }

    /// <summary>子基準系の操作集合 (rChild,tChild、無還元) が候補設定と原点シフト q で一致するかを検証し、
    /// 一致すれば q (candidate 座標系) を返す。既存 TSubgroupFinder.SolveOriginShift の Fraction 厳密版
    /// (k- は index 2/3/4 由来の分数が 1/24 格子スナップに乗りにくいため厳密演算が必要)。260705Cl 追加。</summary>
    private static Fraction[] SolveOriginShiftK(int[][] rChild, Fraction[][] tChild, CandidateData cand)
    {
        int m = rChild.Length;
        var setB = new HashSet<string>();
        for (int i = 0; i < cand.LinKeys.Length; i++)
            foreach (var c in cand.Centering)
                setB.Add(KeyOfK(cand.LinKeys[i], RationalMatrix.ModVec1(RationalMatrix.AddVec(cand.Rt[i], c))));

        int pivot = -1;
        Fraction[] rmiInv = null;
        for (int i = 0; i < m; i++)
        {
            Fraction[] rmi = [rChild[i][0] - 1, rChild[i][1], rChild[i][2], rChild[i][3], rChild[i][4] - 1, rChild[i][5], rChild[i][6], rChild[i][7], rChild[i][8] - 1];
            var invRmi = RationalMatrix.Invert3(rmi);
            if (invRmi != null) { pivot = i; rmiInv = invRmi; break; }
        }

        var qCands = new List<Fraction[]>();
        if (pivot >= 0)
        {
            int candIdx = -1;
            for (int j = 0; j < cand.LinKeys.Length; j++)
                if (SameIntVec(cand.LinKeys[j], rChild[pivot])) { candIdx = j; break; }
            if (candIdx < 0) return null; // 呼び出し元で線形部集合の一致は確認済みのため理論上起きない

            foreach (var cc in cand.Centering)
            {
                var ts = RationalMatrix.AddVec(cand.Rt[candIdx], cc);
                for (int nx = -1; nx <= 1; nx++)
                    for (int ny = -1; ny <= 1; ny++)
                        for (int nz = -1; nz <= 1; nz++)
                        {
                            Fraction[] d = [ts[0] - tChild[pivot][0] + nx, ts[1] - tChild[pivot][1] + ny, ts[2] - tChild[pivot][2] + nz];
                            var qc = RationalMatrix.ModVec1(RationalMatrix.MulVec(rmiInv, d));
                            if (!qCands.Any(x => RationalMatrix.VecEquals(x, qc))) qCands.Add(qc);
                        }
            }
        }
        else
        {
            for (int i = 0; i < 24; i++)
                for (int j = 0; j < 24; j++)
                    for (int k2 = 0; k2 < 24; k2++)
                        qCands.Add([new Fraction(i, 24), new Fraction(j, 24), new Fraction(k2, 24)]);
        }

        foreach (var qc in qCands)
        {
            var setA = new HashSet<string>();
            bool ok = true;
            for (int i = 0; i < m && ok; i++)
            {
                var rq = RationalMatrix.MulVec(RationalMatrix.FromInt(rChild[i]), qc);
                var shift = RationalMatrix.SubVec(rq, qc); // (R-I)q
                var t2 = RationalMatrix.AddVec(tChild[i], shift);
                foreach (var cc in cand.Centering)
                {
                    var key = KeyOfK(rChild[i], RationalMatrix.ModVec1(RationalMatrix.AddVec(t2, cc)));
                    if (!setB.Contains(key)) { ok = false; break; }
                    setA.Add(key);
                }
            }
            if (ok && setA.Count == setB.Count) return qc;
        }
        return null;
    }

    private static string KeyOfK(int[] r, Fraction[] t) => $"{string.Join(" ", r)}|{t[0]}/{t[1]}/{t[2]}";
    #endregion
}
