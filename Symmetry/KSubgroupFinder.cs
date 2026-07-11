// 260705Cl 新規: klassengleiche (k-) 部分群の実行時計算エンジン (Phase 2c)。
// TSubgroupFinder (t-) とは別クラス (計画書 §4.2: 列挙対象が異なるため共有しない。
// 型同定・操作集合完全一致検証・格子ユーティリティは今後の共通化候補)。
//
// Step 1 (部分格子 T′ の列挙) → Step 2 (complement 列挙・共役類分け) → Step 3 (型同定 IdentifyK) →
// Step 4 (極大性判定 + GroupRelation(Kind=K) への配線、GetMaximalKSubgroups) まで実装済み。
// 未対応: k- 専用の軌道分裂・New reflections (親胞 mod1 前提の t- 用ロジックは流用不可、
// FormGroupRelations 側でガード済み)、isomorphic 系列の UI 分離表示 (Phase 2d)。
// 詳細ロードマップ: .project-guidance/ReciPro_FormGroupRelations改修計画.md §4。
using System;
using System.Collections.Concurrent; // 260708Cl 追加: 並列化に伴うキャッシュのスレッド安全化
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Threading; // 260708Cl 追加: 並列化 (LazyThreadSafetyMode)
using System.Threading.Tasks; // 260708Cl 追加: 並列化 (Parallel/Task)

namespace Crystallography;

/// <summary>k-最大部分群の列挙・型同定 (Phase 2c)。公開エントリポイントは <see cref="GetMaximalKSubgroups"/>。</summary>
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

    internal static bool SameIntVec(int[] a, int[] b)
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
    internal sealed class PointGroupData // 260709Cl: NormalizerFinder と共有するため private → internal
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

    // 260708Cl 追加: PointGroupData は純関数の結果で不変。IdentifyK がクラスごとに呼ぶため
    // (立方晶で 1 タイプに複数クラス)、キャッシュ化して再構築を排除する (並列呼び出しにも安全)。
    private static readonly ConcurrentDictionary<int, PointGroupData> _pgCache = new();

    internal static PointGroupData BuildPointGroupData(int sn) => _pgCache.GetOrAdd(sn, BuildPointGroupDataCore); // 260709Cl: private → internal

    private static PointGroupData BuildPointGroupDataCore(int sn)
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
        return EnumerateComplementsCore(pg, BuildQuotient(pg, hnf)); // 260708Cl: core 分離 (ComputeMaximalK が memo 化した q を再利用するため)
    }

    /// <summary>260708Cl 追加: 構築済みの (pg, q) を受け取る本体 (旧 EnumerateComplements の中身)。</summary>
    private static List<int[]> EnumerateComplementsCore(PointGroupData pg, QuotientData q)
    {
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

    internal static int[] MatMulInt(int[] a, int[] b) // 260709Cl: private → internal
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

    //private static readonly Dictionary<int, List<int[]>> _unimodularCache = []; // 260708Cl: 並列化に伴い ConcurrentDictionary+Lazy へ
    private static readonly ConcurrentDictionary<int, Lazy<List<int[]>>> _unimodularCache = new();

    /// <summary>成分が [-k,k] の unimodular (det=±1) 整数 3×3 行列を総当たりで列挙する (キャッシュ済み)。260705Cl 追加。
    /// 260708Cl: 複数スレッドから呼ばれるため Lazy (ExecutionAndPublication) で二重構築を防ぐ (k=3 は 7^9 ≈ 4千万の再帰列挙で高価)。</summary>
    internal static List<int[]> SmallUnimodular(int k) // 260709Cl: private → internal
        => _unimodularCache.GetOrAdd(k, kk => new Lazy<List<int[]>>(() =>
        {
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
                for (int v = -kk; v <= kk; v++) { m[idx] = v; Rec(idx + 1); }
            }
            Rec(0);
            return result;
        }, LazyThreadSafetyMode.ExecutionAndPublication)).Value;

    /// <summary>unimodular 整数行列 (row-major 9 要素) の逆行列 U⁻¹ = det(U)·adj(U) を整数で返す (det=±1 前提)。260708Cl 追加。</summary>
    internal static int[] AdjTimesDet(int[] m, int det) => // 260709Cl: private → internal
    [
        det * (m[4] * m[8] - m[5] * m[7]), det * (m[2] * m[7] - m[1] * m[8]), det * (m[1] * m[5] - m[2] * m[4]),
        det * (m[5] * m[6] - m[3] * m[8]), det * (m[0] * m[8] - m[2] * m[6]), det * (m[2] * m[3] - m[0] * m[5]),
        det * (m[3] * m[7] - m[4] * m[6]), det * (m[1] * m[6] - m[0] * m[7]), det * (m[0] * m[4] - m[1] * m[3]),
    ];

    /// <summary>候補設定 candSn の点群線形部・primitive 基底・その逆・primitive 座標での線形部・
    /// 恒等線形部の中心化並進 (centering cosets) をまとめたキャッシュ。260705Cl 追加。</summary>
    private sealed class CandidateData
    {
        public int[][] LinKeys { get; init; }
        public Fraction[] CInv { get; init; }
        public int[][] ACand { get; init; }         // C⁻¹ R C (整数、candidate primitive 座標)
        public Fraction[][] Rt { get; init; }        // LinKeys[i] に対応する実際の並進 (candidate 慣用胞座標、代表1つ、無還元)
        public Fraction[][] Centering { get; init; } // 恒等線形部の中心化並進 (candidate 慣用胞座標、mod1、重複無し)
        public int CDetSign { get; init; }           // 260708Cl 追加: sign(det C)。det(P)=det(S)·det(U)·det(C⁻¹) の符号事前判定用
    }

    //private static readonly Dictionary<int, CandidateData> _candidateCache = []; // 260708Cl: 並列化に伴い ConcurrentDictionary へ
    //private static Dictionary<int, List<int>> _candidatesByOrder;
    private static readonly ConcurrentDictionary<int, CandidateData> _candidateCache = new();
    private static readonly Lazy<Dictionary<int, List<int>>> _candidatesByOrder = new(BuildCandidatesByOrder, LazyThreadSafetyMode.ExecutionAndPublication); // 260708Cl

    // 260708Cl: 並列呼び出しで同じ sn を同時構築しても純関数なので無害 (最初の格納が勝つ)。
    private static CandidateData BuildCandidateData(int sn) => _candidateCache.GetOrAdd(sn, BuildCandidateDataCore);

    private static CandidateData BuildCandidateDataCore(int sn)
    {
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

        return new CandidateData { LinKeys = [.. linKeys], CInv = cInv, ACand = aCand, Rt = [.. rt], Centering = [.. centering], CDetSign = RationalMatrix.Det3(c).Sign }; // 260708Cl: CDetSign 追加。格納は GetOrAdd 側
        //_candidateCache[sn] = data;
        //return data;
    }

    /// <summary>点群位数 (相異なる線形部数) ごとに候補設定の通し番号をまとめた索引 (初回のみ全 530 設定を走査)。
    /// U 探索のたびに全設定を試すコストを避けるための絞り込み。260705Cl 追加。
    /// 260708Cl: Lazy 化 (_candidatesByOrder) + 全 530 設定の候補データ前計算を並列化。</summary>
    private static Dictionary<int, List<int>> BuildCandidatesByOrder()
    {
        Parallel.For(1, SymmetryStatic.TotalSpaceGroupNumber, sn =>
        {
            if (SymmetryStatic.Symmetries[sn].SpaceGroupNumber != 0)
                BuildCandidateData(sn); // キャッシュ warm (純関数、順序不問)
        });
        var map = new Dictionary<int, List<int>>();
        for (int sn = 1; sn < SymmetryStatic.TotalSpaceGroupNumber; sn++)
        {
            if (SymmetryStatic.Symmetries[sn].SpaceGroupNumber == 0) continue;
            int order = BuildCandidateData(sn).LinKeys.Length;
            if (!map.TryGetValue(order, out var list)) map[order] = list = [];
            list.Add(sn);
        }
        return map;
    }

    /// <summary>親空間群 parentSn の complement (T′=hnf, σ=sigma) を型同定する。
    /// 成功時 (childSn, P, p)（x_parent = P·x_child + p、親慣用胞座標、Fraction）、失敗時 (-1, null, null)。
    /// 260705Cl 追加 (Phase 2c Step3)。</summary>
    public static (int Child, Fraction[] P, Fraction[] Shift) IdentifyK(int parentSn, int[] hnf, int[] sigma)
        => IdentifyK(parentSn, hnf, sigma, null);

    /// <summary>260708Cl: 呼び出し元 (ComputeMaximalK) が memo 化した QuotientData を渡せる内部版。
    /// QuotientData が private 型のため public 側にデフォルト引数では追加できず、オーバーロードに分離。</summary>
    private static (int Child, Fraction[] P, Fraction[] Shift) IdentifyK(int parentSn, int[] hnf, int[] sigma, QuotientData q)
    {
        var pg = BuildPointGroupData(parentSn);
        q ??= BuildQuotient(pg, hnf); // 260708Cl (旧: var q = BuildQuotient(pg, hnf);)
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

        var byOrder = _candidatesByOrder.Value; // 260708Cl (旧: CandidatesByOrder())
        if (!byOrder.TryGetValue(m, out var candList)) return (-1, null, null);

        // 260708Cl: ITA 慣行 (右手系を保つ det(P)>0 の胞変換) を優先する 2 パス探索 (pass=0 は det(P)>0 のみ)。
        // 実 GUI 目視で Pm-3m→Fm-3m が P=[[0,-2,0],[-2,0,0],[0,0,2]] (det=-8、左手系) と同定され、変換行列と
        // 親分数指数の表示が非慣行になっていた。任意の有効同定には det>0 の変換が必ず存在する
        // (右手系どうしの基底取り替え) ため、パス 1 (det<0 のみ) は防御的フォールバック。
        // sign(det P) = sign(det S)·sign(det U)·sign(det C⁻¹) (行列式の乗法性) なので、各パスで試すべき
        // det(U) の符号は (S, 候補 C) から事前に決まる。これで 2 パス化しても U ごとの共役フィルタは
        // どちらか一方のパスでしか走らず、総探索量は 1 パス時と同じに保たれる。
        int sSign = RationalMatrix.Det3(s).Sign;
        // 260708Cl (/simplify): 候補データを事前に配列へホイスト。旧実装は最内ループ (u × 候補) で毎回
        // BuildCandidateData (ConcurrentDictionary.GetOrAdd) を引いていたが、candList は呼び出し中不変。
        var cands = new CandidateData[candList.Count];
        for (int ci = 0; ci < candList.Count; ci++) cands[ci] = BuildCandidateData(candList[ci]);
        //var candDetSigns = candList.Select(c => BuildCandidateData(c).CDetSign).Distinct().ToArray();
        var candDetSigns = cands.Select(c => c.CDetSign).Distinct().ToArray();
        foreach (var (pass, k) in new[] { (0, 1), (0, 2), (0, 3), (1, 1), (1, 2), (1, 3) })
        {
            int wanted = pass == 0 ? 1 : -1; // sign(det P) の目標
            // 260709Cl 追加: 同じ (pass, k) バケット内では最初のヒットで確定せず全候補を評価し、最も簡明な P
            // (①非対角の非ゼロ成分数 ②負成分数 ③成分絶対値和 の辞書式最小、同点は列挙順の先勝ちで決定的) を返す。
            // 旧実装は SmallUnimodular の列挙順 (成分 -k..k の辞書順) で最初に通った U を採用したため、
            // 対角 P=2I で書ける Pm-3m→Fm-3m が [[0,2,0],[2,0,0],[0,0,2]] のような置換込みの P で表示され得た
            // (det>0 は 260708Cl の 2 パス化で保証済みだが最簡形は未保証だった)。バケット単位なので
            // det(P)>0 優先・k 小優先のフォールバック構造は変わらない。
            (int Child, Fraction[] P, Fraction[] Shift) best = (-1, null, null);
            var bestScore = (int.MaxValue, int.MaxValue, double.MaxValue);
            foreach (var u in SmallUnimodular(k))
            {
                int uDet = u[0] * (u[4] * u[8] - u[5] * u[7]) - u[1] * (u[3] * u[8] - u[5] * u[6]) + u[2] * (u[3] * u[7] - u[4] * u[6]);
                // 全候補の det(C) 符号が同一 (通常ケース) なら、共役フィルタの前に u ごと枝刈りできる。
                if (candDetSigns.Length == 1 && sSign * uDet * candDetSigns[0] != wanted) continue;
                //var uFrac = RationalMatrix.FromInt(u);
                //var uInv = RationalMatrix.Invert3(uFrac); // det=±1 なので必ず存在
                //
                //var conjugated = new int[m][];
                //bool intAll = true;
                //for (int i = 0; i < m && intAll; i++)
                //{
                //    var cf = RationalMatrix.Mul(RationalMatrix.Mul(uInv, RationalMatrix.FromInt(aH[i])), uFrac);
                //    var ci = RationalMatrix.ToIntOrNull(cf);
                //    if (ci == null) { intAll = false; break; } // U が unimodular なら理論上必ず整数 (防御的チェック)
                //    conjugated[i] = ci;
                //}
                //if (!intAll) continue;
                // 260708Cl: U は unimodular 整数行列なので U⁻¹ = det(U)·adj(U) も厳密に整数。最内ループの共役
                // U⁻¹·A_H·U を Fraction (BigInteger) から純整数演算へ置換 (数学的に同値のまま桁違いに速い。
                // この共役フィルタが k-エンジン全体の支配的ホットループだった。整数なので ToIntOrNull 検査も不要)。
                var uInvInt = AdjTimesDet(u, uDet);
                var conjugated = new int[m][];
                for (int i = 0; i < m; i++)
                    conjugated[i] = MatMulInt(MatMulInt(uInvInt, aH[i]), u);

                //foreach (var candSn in candList) // 260708Cl (/simplify): ホイスト済み cands を index で走査 (辞書引き排除)
                //{
                //    var cand = BuildCandidateData(candSn);
                for (int ci = 0; ci < cands.Length; ci++)
                {
                    var cand = cands[ci];
                    if (sSign * uDet * cand.CDetSign != wanted) continue; // 260708Cl: このパスの det(P) 符号と不一致
                    if (!SetEqualsIntMats(conjugated, cand.ACand)) continue;

                    var p = RationalMatrix.Mul(RationalMatrix.Mul(s, RationalMatrix.FromInt(u)), cand.CInv); // 260708Cl: uFrac は共役フィルタ整数化で不要になったためここで直接変換
                    //if (pass == 0 && RationalMatrix.Det3(p).Sign < 0) continue; // 260708Cl: 上の符号事前判定に置換 (乗法性より等価)
                    // 260709Cl: 即 return せずバケット内の最簡 P を選ぶ (①非対角の非ゼロ ②負成分 ③絶対値和の
                    // 辞書式最小、同点は先勝ち = 決定的)。スコア判定は重い検証 (rChild 構築 + SolveOriginShiftK、
                    // Fraction/BigInteger 演算) の前に行い、既知 best を改善しない候補はここで棄却する
                    // (同一操作集合には対称性の分だけ多数の U が通るため、全候補で origin 解決すると
                    // 2/m 逆引きが 1.4 s → 416 s になった)。ただし最初の成功 (best 未確定) までは棄却しない。
                    // フォールバックバケット (k ≥ 2、希少経路) は従来どおり首ヒット確定で、最簡化は k=1 のみ。
                    var score = ScoreP(p);
                    bool perfect = score.OffDiagNonZero == 0 && score.Negatives == 0; // 対角・全非負 = それ以上本質的に改善しない
                    if (k == 1 && best.Child >= 0 && !perfect && score.CompareTo(bestScore) >= 0)
                        continue;
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
                    //return (candSn, p, pShift);
                    //return (candList[ci], p, pShift); // 260708Cl (/simplify): candSn → candList[ci]
                    if (k > 1 || perfect)
                        return (candList[ci], p, pShift);
                    if (best.Child < 0 || score.CompareTo(bestScore) < 0)
                    {
                        best = (candList[ci], p, pShift);
                        bestScore = score;
                    }
                }
            }
            if (best.Child >= 0) // 260709Cl: バケット内で見つかったらフォールバック (次の k / det<0 パス) へは進まない
                return best;
        }
        return (-1, null, null);
    }

    /// <summary>260709Cl 追加: 同定候補の基底変換 P の「簡明さ」スコア (小さいほど良い)。
    /// ITA の表記慣行に合わせ、①非対角の非ゼロ成分が少ない (対角形優先) ②負成分が少ない
    /// ③成分絶対値和が小さい、の辞書式で比較する。</summary>
    private static (int OffDiagNonZero, int Negatives, double AbsSum) ScoreP(Fraction[] p)
    {
        int offDiag = 0, negatives = 0;
        double absSum = 0;
        for (int i = 0; i < 9; i++)
        {
            int sign = p[i].Sign;
            if (sign != 0 && i % 4 != 0) offDiag++; // i=0,4,8 が対角
            if (sign < 0) negatives++;
            absSum += Math.Abs((double)p[i].Num / (double)p[i].Den);
        }
        return (offDiag, negatives, absSum);
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

        // 260705Cl 修正 (codex R5 指摘): 旧実装は候補 q を nx,ny,nz∈{-1,0,1} の 27 通りに固定していたが、
        // これは Z³/(R-I)Z³ の完全代表系である保証がない (例: det(R-I)=-4 の roto-inversion では
        // mod4 の代表 {-1,0,1}≡{3,0,1} が剰余 "2" を取り逃がす)。Z³/(R-I)Z³ (位数|det(R-I)|) を
        // Step2 の coset 代表構築と同じ手法 (0..|det|-1 の箱を総当たり→正準ラベルで重複除去) で
        // 安全に列挙する。最初に見つかった pivot で q が見つからなければ、他の full-rank pivot も試す
        // (実装途中の中心化展開・候補対応のズレを拾いやすくするための保険、codex 提案)。
        for (int pivot = 0; pivot < m; pivot++)
        {
            var rmiInt = new[] { rChild[pivot][0] - 1, rChild[pivot][1], rChild[pivot][2], rChild[pivot][3], rChild[pivot][4] - 1, rChild[pivot][5], rChild[pivot][6], rChild[pivot][7], rChild[pivot][8] - 1 };
            var rmi = RationalMatrix.FromInt(rmiInt);
            var rmiInv = RationalMatrix.Invert3(rmi);
            if (rmiInv == null) continue; // det(R-I)=0: この操作は origin shift を決められない (screw/glide 等)

            int candIdx = -1;
            for (int j = 0; j < cand.LinKeys.Length; j++)
                if (SameIntVec(cand.LinKeys[j], rChild[pivot])) { candIdx = j; break; }
            if (candIdx < 0) continue; // 呼び出し元で線形部集合の一致は確認済みのため理論上起きない

            var zReps = QuotientRepsForFullRankM(rmiInt);
            var qCands = new List<Fraction[]>();
            foreach (var cc in cand.Centering)
            {
                var diff = RationalMatrix.SubVec(RationalMatrix.AddVec(cand.Rt[candIdx], cc), tChild[pivot]);
                foreach (var z in zReps)
                {
                    var qc = RationalMatrix.ModVec1(RationalMatrix.MulVec(rmiInv, RationalMatrix.AddVec(diff, z)));
                    if (!qCands.Any(x => RationalMatrix.VecEquals(x, qc))) qCands.Add(qc);
                }
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
        }

        // 全操作が並進のみを固定する特殊ケース (全線形部で det(R-I)=0、通常起きない) — 1/24 格子総当たりへフォールバック
        for (int i = 0; i < 24; i++)
            for (int j = 0; j < 24; j++)
                for (int k2 = 0; k2 < 24; k2++)
                {
                    var qc = new Fraction[] { new(i, 24), new(j, 24), new(k2, 24) };
                    var setA = new HashSet<string>();
                    bool ok = true;
                    for (int ii = 0; ii < m && ok; ii++)
                    {
                        var rq = RationalMatrix.MulVec(RationalMatrix.FromInt(rChild[ii]), qc);
                        var shift = RationalMatrix.SubVec(rq, qc);
                        var t2 = RationalMatrix.AddVec(tChild[ii], shift);
                        foreach (var cc in cand.Centering)
                        {
                            var key = KeyOfK(rChild[ii], RationalMatrix.ModVec1(RationalMatrix.AddVec(t2, cc)));
                            if (!setB.Contains(key)) { ok = false; break; }
                            setA.Add(key);
                        }
                    }
                    if (ok && setA.Count == setB.Count) return qc;
                }
        return null;
    }

    //private static readonly Dictionary<string, List<Fraction[]>> _quotientRepsCache = []; // 260708Cl: 並列化に伴い ConcurrentDictionary へ
    private static readonly ConcurrentDictionary<string, List<Fraction[]>> _quotientRepsCache = new();

    /// <summary>整数行列 M (det≠0) について Z³/MZ³ (位数|det(M)|) の完全代表系を安全に構築する
    /// (0..|det|-1 の箱を総当たりし、M⁻¹v mod1 を正準ラベルとして重複除去、|det| 個集まったら確定)。
    /// M·adj(M)=det(M)·I なので、この箱に完全代表系が含まれることは保証される。260705Cl 追加 (Phase 2c Step3)。
    /// 260708Cl: 並列呼び出しで同じ M を同時構築しても純関数なので無害 (最初の格納が勝つ)。</summary>
    private static List<Fraction[]> QuotientRepsForFullRankM(int[] m)
        => _quotientRepsCache.GetOrAdd(string.Join(",", m), _ =>
        {
            var mInv = RationalMatrix.Invert3(RationalMatrix.FromInt(m)) ?? throw new InvalidOperationException("M is singular");
            int d = Math.Abs((int)RationalMatrix.Det3(RationalMatrix.FromInt(m)).Num);

            var reps = new List<Fraction[]>();
            var seen = new HashSet<string>();
            for (int x = 0; x < d && reps.Count < d; x++)
                for (int y = 0; y < d && reps.Count < d; y++)
                    for (int z = 0; z < d && reps.Count < d; z++)
                    {
                        Fraction[] v = [x, y, z];
                        var label = RationalMatrix.ModVec1(RationalMatrix.MulVec(mInv, v));
                        if (seen.Add($"{label[0]}/{label[1]}/{label[2]}")) reps.Add(v);
                    }
            if (reps.Count != d)
                throw new InvalidOperationException($"failed to enumerate Z^3/MZ^3 (got {reps.Count}, expect {d})");
            return reps;
        });

    private static string KeyOfK(int[] r, Fraction[] t) => $"{string.Join(" ", r)}|{t[0]}/{t[1]}/{t[2]}";
    #endregion

    #region 極大 k-部分群 (Step 4, GroupRelation(Kind=K) への配線)
    // 260705Cl 追加 (Phase 2c Step4)。設計は codex との5回目の相談で確定
    // (.project-guidance/ReciPro_FormGroupRelations改修計画.md §4.1 item5)。
    //
    // 極大性判定: index n の complement H (T′=hnf) が非極大なのは、T′⊊T″⊊T な point-group-invariant
    // 中間格子 T″ (index n の真の約数、index2 のみ既存 index=2/3/4 列挙内で該当し得る) が存在し、
    // かつその T″ 上の何らかの complement H″ が H を包含する (各線形部で H の並進が H″ の並進と
    // T″ を法として一致する) ときに限る。index2/3 は中間指数が無いので自動的に極大。
    //private static readonly Dictionary<int, GroupRelation[]> _kCache = []; // 260708Cl: 並列化に伴い per-type Lazy へ
    //private static readonly object _kLock = new();
    private static readonly ConcurrentDictionary<int, Lazy<GroupRelation[]>> _kCache = new();

    /// <summary>親空間群 (通し番号) の maximal k-部分群を共役類単位で返す (index 2,3,4)。計算は初回のみ (キャッシュ)。
    /// 260705Cl 追加 (Phase 2c Step4)。260708Cl: グローバルロックを per-type Lazy に置換 — 異なるタイプは
    /// 並列に計算でき、同一タイプは ExecutionAndPublication で一度だけ計算される (二重計算・ブロッキング最小)。</summary>
    public static GroupRelation[] GetMaximalKSubgroups(int seriesNumber)
        => _kCache.GetOrAdd(seriesNumber, sn => new Lazy<GroupRelation[]>(() => ComputeMaximalK(sn), LazyThreadSafetyMode.ExecutionAndPublication)).Value;

    // 260708Cl 追加 (Phase 2d 後段): k-超群 (minimal k-supergroup) の逆引き。klassengleiche は幾何結晶類
    // (Schoenflies 点群) を変えないため、候補親は同じ PointGroupSFStr のタイプに限られる
    // (t-超群索引のような全 230 型走査は不要。この不変量は tools/SymmetryPropsCheck の広域 sweep で検証)。
    // 候補親の GetMaximalKSubgroups (キャッシュ共有) から子タイプ == itNumber の関係を拾う。
    // 同型 (Kind=Isomorphic) も klassengleiche なので含まれる (自分自身が超群になる場合を含む)。
    //private static readonly Dictionary<int, GroupRelation[]> _kSupergroupCache = []; // 260708Cl: ConcurrentDictionary へ
    private static readonly ConcurrentDictionary<int, GroupRelation[]> _kSupergroupCache = new();

    /// <summary>itNumber (IT 番号) の k-超群逆引きが計算済みか。初回計算は同じ結晶類の全タイプの k-部分群計算を
    /// 伴い重い場合があるため、GUI はこれを見てバックグラウンド構築 + 「計算中…」表示を選べる。260708Cl 追加。</summary>
    public static bool KSupergroupsReady(int itNumber) => _kSupergroupCache.ContainsKey(itNumber);

    /// <summary>指定タイプ (IT 番号) を maximal k-部分群 (同型含む) に持つ関係 (= minimal k-supergroup) を返す。
    /// ParentSeriesNumber = 超群 (第 1 設定)、ChildSeriesNumber = itNumber 側の設定。260708Cl 追加 (Phase 2d 後段)。
    /// 260708Cl 並列化: 候補タイプ (同一結晶類) ごとの k-部分群計算は独立なので Parallel 実行する
    /// (per-type Lazy キャッシュにより同一タイプは一度だけ計算・他呼び出しと共有)。結果の組み立ては
    /// IT 番号昇順で逐次行い、表示順の決定性を保つ。同一 itNumber の並行呼び出しは二重集計になり得るが、
    /// type 結果はキャッシュ済みで 2 回目は安価、格納結果も同一なので無害。</summary>
    public static GroupRelation[] GetMinimalKSupergroups(int itNumber)
    {
        if (_kSupergroupCache.TryGetValue(itNumber, out var cached)) return cached;
        string sfTarget = SymmetryStatic.Symmetries[SymmetryStatic.GetSeriesNumber(itNumber, 1)].PointGroupSFStr;
        var candSns = new List<int>();
        for (int it = 1; it <= 230; it++)
        {
            int sn = SymmetryStatic.GetSeriesNumber(it, 1);
            if (sn >= 0 && SymmetryStatic.Symmetries[sn].PointGroupSFStr == sfTarget)
                candSns.Add(sn);
        }
        var perType = new GroupRelation[candSns.Count][];
        Parallel.For(0, candSns.Count, i => perType[i] = GetMaximalKSubgroups(candSns[i]));
        var list = new List<GroupRelation>();
        foreach (var subs in perType)
            foreach (var sub in subs)
                if (sub.ChildSeriesNumber >= 0 && SymmetryStatic.Symmetries[sub.ChildSeriesNumber].SpaceGroupNumber == itNumber)
                    list.Add(sub);
        return _kSupergroupCache.GetOrAdd(itNumber, [.. list]);
    }

    #region normalizer 軌道 (Phase 2, 260709Cl)
    // G-共役類 (GetMaximalKSubgroups の ConjugacyClassId) を affine normalizer N_Aff(G) の軌道に束ねる。
    // 「固定 index ごとの normalizer 軌道」が ITA A1 の同型 (IIc) 系列表示の分類粒度に対応する (codex R9。
    // 異なる素数 index を 1 本の系列式にまとめる記号処理は Phase 3 の後段)。
    // 作用: n=(U,a) は点群不変な部分格子 T′ を U·T′ へ、complement の操作 g_i=(A_i, t_i) を
    // (A_{π(i)}, U·t_i+(I−A_{π(i)})·a) へ写す。N は G を正規化し「index n の極大部分群」の集合を保つため、
    // 作用先は必ず列挙済みの極大クラスのどれかに一致する (見つからなければ hard fail)。
    // union-find で辺を張るだけで連結成分 = 生成群の軌道になる (作用は全単射、無向成分に逆元も含まれる)。
    // ⚠ 生成集合は NormalizerFinder の BoundedVerified(k=1) — 軌道が「粗すぎる」ことは原理的に無いが
    // 「細かすぎる」(本来 1 軌道が複数に見える) 可能性は残る。既知ケース照合 (PART 11) で監視する。

    /// <summary>ComputeMaximalK が保存する classId 順の生データ (HNF と共役類の全メンバー σ)。260709Cl 追加。</summary>
    private static readonly ConcurrentDictionary<int, (int[] Hnf, List<int[]> Members)[]> _rawClassesCache = new();

    private static readonly ConcurrentDictionary<int, Lazy<int[]>> _orbitCache = new();

    /// <summary>GetMaximalKSubgroups(seriesNumber) の各共役類を N_Aff(G) の軌道へ束ねた軌道 ID
    /// (0 始まり、最小 classId 順の連番) を返す。orbits[rel.ConjugacyClassId] が rel の軌道 ID。
    /// K と Isomorphic の両方の類を含む (軌道は Kind・index・子タイプを保つ)。260709Cl 追加 (Phase 2)。</summary>
    public static int[] GetNormalizerOrbits(int seriesNumber)
        => _orbitCache.GetOrAdd(seriesNumber, snn => new Lazy<int[]>(() => ComputeNormalizerOrbits(snn), LazyThreadSafetyMode.ExecutionAndPublication)).Value;

    private static int[] ComputeNormalizerOrbits(int sn)
    {
        _ = GetMaximalKSubgroups(sn); // per-series Lazy: _rawClassesCache を確実に埋める
        return ComputeNormalizerOrbitsCore(sn, _rawClassesCache[sn]);
    }

    /// <summary>260709Cl (Phase 3): 軌道計算の本体 (rawClasses をパラメータ化し、index≤4 の既定リストと
    /// 拡張列挙 (GetNormalizerOrbitsAt) で共有)。</summary>
    private static int[] ComputeNormalizerOrbitsCore(int sn, (int[] Hnf, List<int[]> Members)[] rawClasses)
    {
        var pg = BuildPointGroupData(sn);
        int m = pg.LinKeys.Length;
        var nd = NormalizerFinder.Get(sn);

        // hnf ごとの Quotient / その hnf 上の classId 一覧。260709Cl (codex R11): キーは格子の一意な
        // canonical HNF 文字列 — 像格子 U·T′ の照合を線形探索 (IsSameLattice 全走査) から辞書 O(1) へ
        // (P1 の高指数では 類数 × 生成元数 × HNF 数 の積が爆発していた)。
        var hnfKeys = new string[rawClasses.Length];
        var qByKey = new Dictionary<string, QuotientData>();
        var classesByKey = new Dictionary<string, List<int>>();
        var hnfByKey = new Dictionary<string, int[]>();
        for (int c = 0; c < rawClasses.Length; c++)
        {
            var key = string.Join(",", CanonicalHnf(rawClasses[c].Hnf));
            hnfKeys[c] = key;
            if (!classesByKey.TryGetValue(key, out var list))
            {
                classesByKey[key] = list = [];
                hnfByKey[key] = rawClasses[c].Hnf;
                qByKey[key] = BuildQuotient(pg, rawClasses[c].Hnf);
            }
            list.Add(c);
        }

        // 生成元 = 非自明線形部 + 純並進核の離散生成元。連続方向 (polar) は G を点ごとに centralize する
        // ため部分群への作用が恒等であり、含めない (codex R9)。
        int[] idU = [1, 0, 0, 0, 1, 0, 0, 0, 1];
        var gens = new List<(int[] U, Fraction[] Shift)>();
        foreach (var g in nd.Generators)
            gens.Add((g.LinearPrimitive, g.ShiftPrimitive));
        foreach (var d in nd.TranslationKernel.DiscreteGenerators)
            gens.Add((idU, d));

        // union-find (260709Cl codex R11: 連結成分数を追跡し、1 まで潰れたら残りの生成元をスキップ —
        // P1 のように生成元が数千あっても最初の数個で 1 軌道に潰れるケースの実効コストを抑える)
        var parent = new int[rawClasses.Length];
        for (int i = 0; i < parent.Length; i++) parent[i] = i;
        int components = rawClasses.Length;
        int Find(int x) { while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; } return x; }
        void Union(int x, int y) { x = Find(x); y = Find(y); if (x != y) { parent[Math.Max(x, y)] = Math.Min(x, y); components--; } }

        foreach (var (u, shift) in gens)
        {
            if (components <= 1) break; // 260709Cl: これ以上束ねようがない
            var perm = NormalizerFinder.FindConjugationPermutation(pg, u)
                ?? throw new InvalidOperationException($"normalizer generator does not normalize the point group (sn={sn})");
            for (int c = 0; c < rawClasses.Length; c++)
            {
                var (hnf, members) = rawClasses[c];
                var q = qByKey[hnfKeys[c]];
                // T″ = U·T′。同 index の点群不変 HNF (クラスを持つもの) の中に必ずある。
                // 260709Cl (codex R11): canonical HNF 化して辞書 O(1) 照合 (旧: IsSameLattice の線形探索)。
                var uh = MatMulInt(u, hnf);
                string key2 = string.Join(",", CanonicalHnf(uh));
                var hnf2 = classesByKey.ContainsKey(key2) ? hnfByKey[key2]
                    : throw new InvalidOperationException($"image sublattice not found among enumerated HNFs (sn={sn})");
                var q2 = qByKey[key2];

                // σ″: 代表 complement の各操作を共役して T″ 上の coset index として読み取る
                var sigma = members[0];
                var sigma2 = new int[m];
                for (int i = 0; i < m; i++)
                {
                    var rep = q.Reps[sigma[i]];
                    Fraction[] t = [pg.T0[i][0] + new Fraction(rep[0]), pg.T0[i][1] + new Fraction(rep[1]), pg.T0[i][2] + new Fraction(rep[2])];
                    var ut = MulIntMatVec(u, t);
                    var ra = MulIntMatVec(pg.A[perm[i]], shift);
                    // δ = U·t + (I−A_π)·a − T0[π(i)]。n が normalizer なら g′∈G ゆえ必ず整数ベクトル。
                    Fraction[] delta =
                    [
                        ut[0] + shift[0] - ra[0] - pg.T0[perm[i]][0],
                        ut[1] + shift[1] - ra[1] - pg.T0[perm[i]][1],
                        ut[2] + shift[2] - ra[2] - pg.T0[perm[i]][2],
                    ];
                    if (!delta[0].IsInteger || !delta[1].IsInteger || !delta[2].IsInteger)
                        throw new InvalidOperationException($"conjugated translation is not integral (sn={sn}) — generator is not in the normalizer");
                    sigma2[perm[i]] = CosetIndexOf(q2, hnf2, [(long)delta[0].Num, (long)delta[1].Num, (long)delta[2].Num]);
                }

                // 作用先の共役類を同定 (共役類メンバーは完全列挙済みなので int[] 一致照合で必ず見つかる)
                int target = -1;
                foreach (var c2 in classesByKey[key2])
                    if (rawClasses[c2].Members.Any(s => SameSigma(s, sigma2))) { target = c2; break; }
                if (target < 0)
                    throw new InvalidOperationException($"image complement not found among maximal classes (sn={sn})");
                Union(c, target);
            }
        }

        // 軌道 ID を最小 classId 順の連番へ正規化 (決定的)
        var orbitId = new Dictionary<int, int>();
        var result = new int[rawClasses.Length];
        for (int c = 0; c < rawClasses.Length; c++)
        {
            int root = Find(c);
            if (!orbitId.TryGetValue(root, out var id))
                orbitId[root] = id = orbitId.Count;
            result[c] = id;
        }
        return result;
    }

    private static bool SameSigma(int[] a, int[] b)
    {
        for (int i = 0; i < a.Length; i++)
            if (a[i] != b[i]) return false;
        return true;
    }

    // 260709Cl (codex R11): 像格子の照合を CanonicalHnf 辞書キーへ置換したため Det3Int/IsSameLattice は不要に。
    //private static long Det3Int(int[] h)
    //    => (long)h[0] * (h[4] * h[8] - h[5] * h[7]) - (long)h[1] * (h[3] * h[8] - h[5] * h[6]) + (long)h[2] * (h[3] * h[7] - h[4] * h[6]);
    ///// <summary>2 つの部分格子基底 (整数、同 det) が同一格子か: H1⁻¹·H2 が整数 (det 同一なので unimodular)。260709Cl 追加。</summary>
    //private static bool IsSameLattice(int[] h1, int[] h2)
    //{
    //    var inv = RationalMatrix.Invert3(RationalMatrix.FromInt(h1));
    //    if (inv == null) return false;
    //    return RationalMatrix.ToIntOrNull(RationalMatrix.Mul(inv, RationalMatrix.FromInt(h2))) != null;
    //}

    /// <summary>整数 3×3 基底 (列ベクトルが格子基底) の張る格子の一意な HNF 正規形を列演算で求める
    /// (EnumerateHnf と同じ規約: 下三角 [a,0,0, x,b,0, y,z,c]、対角正、0 ≤ x &lt; b・0 ≤ y,z &lt; c)。
    /// 格子が同一 ⟺ 正規形が一致するため、辞書キーに使える。260709Cl 追加 (codex R11)。</summary>
    private static int[] CanonicalHnf(int[] mIn)
    {
        var h = (int[])mIn.Clone();
        void SwapCol(int i, int j) { for (int r = 0; r < 3; r++) (h[r * 3 + i], h[r * 3 + j]) = (h[r * 3 + j], h[r * 3 + i]); }
        void AddCol(int dst, int src, int f) { for (int r = 0; r < 3; r++) h[r * 3 + dst] += f * h[r * 3 + src]; }
        void NegCol(int j) { for (int r = 0; r < 3; r++) h[r * 3 + j] = -h[r * 3 + j]; }
        static int FloorDiv(int a, int b) => (int)Math.Floor((double)a / b);

        for (int row = 0; row < 3; row++) // 行 row の右側列 (row+1..2) をユークリッドで 0 化し、対角を正へ
        {
            while (true)
            {
                int best = 0, bc = -1;
                for (int j = row; j < 3; j++)
                    if (h[row * 3 + j] != 0 && (best == 0 || Math.Abs(h[row * 3 + j]) < Math.Abs(best))) { best = h[row * 3 + j]; bc = j; }
                if (bc < 0) throw new InvalidOperationException("singular lattice basis in CanonicalHnf");
                if (bc != row) SwapCol(row, bc);
                if (h[row * 3 + row] < 0) NegCol(row);
                bool done = true;
                for (int j = row + 1; j < 3; j++)
                {
                    int f = FloorDiv(h[row * 3 + j], h[row * 3 + row]);
                    if (f != 0) AddCol(j, row, -f);
                    if (h[row * 3 + j] != 0) done = false;
                }
                if (done) break;
            }
        }
        // reduce: 対角の左側成分を [0, 対角) へ
        if (h[4] > 0) { int f = FloorDiv(h[3], h[4]); if (f != 0) AddCol(0, 1, -f); }
        if (h[8] > 0)
        {
            int f = FloorDiv(h[6], h[8]); if (f != 0) AddCol(0, 2, -f);
            f = FloorDiv(h[7], h[8]); if (f != 0) AddCol(1, 2, -f);
        }
        return h;
    }

    private static Fraction[] MulIntMatVec(int[] mm, Fraction[] v)
        =>
        [
            new Fraction(mm[0]) * v[0] + new Fraction(mm[1]) * v[1] + new Fraction(mm[2]) * v[2],
            new Fraction(mm[3]) * v[0] + new Fraction(mm[4]) * v[1] + new Fraction(mm[5]) * v[2],
            new Fraction(mm[6]) * v[0] + new Fraction(mm[7]) * v[1] + new Fraction(mm[8]) * v[2],
        ];
    #endregion

    private sealed class RawComplement
    {
        public int[] Hnf;
        public int Index;
        public int[] Sigma;
        public bool IsMaximal = true;
    }

    private static GroupRelation[] ComputeMaximalK(int sn)
    {
        // 260709Cl (Phase 3): 本体を ComputeMaximalKCore へ一般化 (index リスト指定)。既存挙動は不変。
        var (rels, rawClasses) = ComputeMaximalKCore(sn, [2, 3, 4], targetIndices: null);
        _rawClassesCache[sn] = rawClasses;
        return rels;
    }

    /// <summary>260709Cl 追加 (Phase 3): 指定 index の極大 k-部分群を列挙する (index ≥ 5 で極大に残るのは
    /// 理論上 isomorphic のみ — 非同型 klassengleiche 極大の index は 2,3,4 [ITA A1]。スピナーによる
    /// 同型系列の拡張列挙に使う)。極大性判定のため index の真の約数の格子・complement も内部で列挙する。
    /// (sn, index) 単位でキャッシュ。</summary>
    public static GroupRelation[] GetMaximalKSubgroupsAt(int seriesNumber, int index)
    {
        // 260709Cl (codex R11): 極大 k-部分群の index は p^r (r ≤ 3) に限る — T/T′ が単純有限 Z[P]-加群
        // ⟺ 極大、単純加群の加法群は (F_p)^r、T の階数 3 より r ≤ 3。非素数冪は Sylow 分解の
        // characteristic な中間格子で必ず落ちるため、列挙せずに即時空を返す。
        if (!IsPrimePowerAtMostCube(index))
            return [];
        return _isoAtCache.GetOrAdd((seriesNumber, index), key => new Lazy<(GroupRelation[] Rels, (int[] Hnf, List<int[]> Members)[] Raw)>(() =>
        {
            var divisors = new List<int>();
            for (int d = 2; d <= key.Item2; d++)
                if (key.Item2 % d == 0) divisors.Add(d); // 1 < d ≤ index (自身含む、真の約数は極大性判定の coarse)
            return ComputeMaximalKCore(key.Item1, [.. divisors], targetIndices: [key.Item2]);
        }, LazyThreadSafetyMode.ExecutionAndPublication)).Value.Rels;
    }

    /// <summary>index が p^r (p 素数、1 ≤ r ≤ 3) か。260709Cl 追加 (codex R11)。</summary>
    private static bool IsPrimePowerAtMostCube(int n)
    {
        if (n < 2) return false;
        int p = 2;
        while (p * p <= n && n % p != 0) p++;
        if (n % p != 0) p = n; // n 自身が素数
        int r = 0;
        while (n > 1) { if (n % p != 0) return false; n /= p; r++; }
        return r <= 3;
    }

    private static readonly ConcurrentDictionary<(int Sn, int Index), Lazy<(GroupRelation[] Rels, (int[] Hnf, List<int[]> Members)[] Raw)>> _isoAtCache = new();
    private static readonly ConcurrentDictionary<(int Sn, int Index), Lazy<int[]>> _orbitAtCache = new(); // 260709Cl (codex R11): 毎回再計算を排除

    /// <summary>GetMaximalKSubgroupsAt(sn, index) の各共役類に対する normalizer 軌道 ID。260709Cl 追加 (Phase 3)。</summary>
    public static int[] GetNormalizerOrbitsAt(int seriesNumber, int index)
        => !IsPrimePowerAtMostCube(index)
            ? []
            : _orbitAtCache.GetOrAdd((seriesNumber, index), key => new Lazy<int[]>(() =>
            {
                _ = GetMaximalKSubgroupsAt(key.Sn, key.Index); // Lazy を確実に評価
                return ComputeNormalizerOrbitsCore(key.Sn, _isoAtCache[key].Value.Raw);
            }, LazyThreadSafetyMode.ExecutionAndPublication)).Value;

    /// <summary>k-部分群列挙の本体 (260709Cl: ComputeMaximalK から一般化)。
    /// enumIndices = HNF/complement を列挙する index の集合 (極大性判定の coarse も含めること)。
    /// targetIndices = GroupRelation を構築する index (null = enumIndices 全部)。</summary>
    private static (GroupRelation[] Rels, (int[] Hnf, List<int[]> Members)[] RawClasses) ComputeMaximalKCore(int sn, int[] enumIndices, int[] targetIndices)
    {
        var pg = BuildPointGroupData(sn);
        int m = pg.LinKeys.Length;

        // 260708Cl 並列化: (index, hnf) ごとの Quotient 構築 + complement 列挙は独立なので並列実行し、
        // 結果は元の列挙順で組み立てる (classId・表示順の決定性維持)。あわせて hnf ごとの
        // BuildQuotient / complement 列挙を memo 化 (旧実装は極大性判定の (fine×coarse) 二重ループ内で
        // 毎回再構築・再列挙していた)。
        var byIndex = new Dictionary<int, List<int[]>>();
        var hnfItems = new List<(int Index, int[] Hnf)>();
        foreach (int index in enumIndices) // 260709Cl (旧: new[] { 2, 3, 4 })
        {
            var inv = FilterPointGroupInvariant(sn, EnumerateHnf(index));
            byIndex[index] = inv;
            foreach (var hnf in inv)
                hnfItems.Add((index, hnf));
        }
        var quotients = new QuotientData[hnfItems.Count];
        var sigmasPerHnf = new List<int[]>[hnfItems.Count];
        Parallel.For(0, hnfItems.Count, i =>
        {
            quotients[i] = BuildQuotient(pg, hnfItems[i].Hnf);
            sigmasPerHnf[i] = EnumerateComplementsCore(pg, quotients[i]);
        });
        var qByHnf = new Dictionary<string, (QuotientData Q, List<int[]> Sigmas)>();
        var raws = new List<RawComplement>();
        for (int i = 0; i < hnfItems.Count; i++)
        {
            qByHnf[string.Join(",", hnfItems[i].Hnf)] = (quotients[i], sigmasPerHnf[i]);
            foreach (var sigma in sigmasPerHnf[i])
                raws.Add(new RawComplement { Hnf = hnfItems[i].Hnf, Index = hnfItems[i].Index, Sigma = sigma });
        }

        // 中間格子を経由できる complement は非極大。260708Cl: fine ごとに独立 (自身の IsMaximal のみ書く) → 並列。
        // 260709Cl (Phase 3): 「index4 → coarse=index2」固定から「fine.Index の真の約数 d (1<d<Index) の
        // 全 coarse」へ一般化 (index 8/9/25/27 等の素数冪で必要。素数 index は約数が無く自動極大)。
        Parallel.ForEach(raws.Where(r => byIndex.Keys.Any(d => d > 1 && d < r.Index && r.Index % d == 0)), fine =>
        {
            var fineQ = qByHnf[string.Join(",", fine.Hnf)].Q;
            foreach (int d in byIndex.Keys.Where(d => d > 1 && d < fine.Index && fine.Index % d == 0))
            {
                foreach (var coarseHnf in byIndex[d])
                {
                    if (!IsLatticeSubset(coarseHnf, fine.Hnf)) continue; // T′(fine) ⊂ T″(coarse) か
                    var (coarseQ, coarseSigmas) = qByHnf[string.Join(",", coarseHnf)];
                    foreach (var coarseSigma in coarseSigmas)
                    {
                        bool contained = true;
                        for (int i = 0; i < m; i++)
                        {
                            int mapped = CosetIndexOf(coarseQ, coarseHnf, fineQ.Reps[fine.Sigma[i]]);
                            if (mapped != coarseSigma[i]) { contained = false; break; }
                        }
                        if (contained) { fine.IsMaximal = false; break; }
                    }
                    if (!fine.IsMaximal) break;
                }
                if (!fine.IsMaximal) break;
            }
        });

        // 共役類分け (hnf ごと、順序決定的に逐次)。型同定 (BuildGroupRelation 内の IdentifyK が支配的コスト) は
        // クラスごとに独立なので並列。classId = 収集順 index で旧逐次カウンタと同一。260708Cl。
        var items = new List<(int[] Hnf, QuotientData Q, List<int[]> Cls)>();
        // 260709Cl (Phase 3): targetIndices 指定時は対象 index の complement だけを類分け・構築する
        // (真の約数の complement は極大性判定にのみ使う)。
        foreach (var grp in raws.Where(r => r.IsMaximal && (targetIndices == null || targetIndices.Contains(r.Index))).GroupBy(r => string.Join(",", r.Hnf)))
        {
            var hnf = grp.First().Hnf;
            var sigmas = grp.Select(r => r.Sigma).ToList();
            foreach (var cls in GroupComplementsByConjugacy(sn, hnf, sigmas))
                items.Add((hnf, qByHnf[string.Join(",", hnf)].Q, cls));
        }
        var rels = new GroupRelation[items.Count];
        Parallel.For(0, items.Count, i =>
            rels[i] = BuildGroupRelation(sn, pg, items[i].Hnf, items[i].Cls[0], i, items[i].Cls.Count, items[i].Q));
        // 260709Cl (Phase 2/3): normalizer 軌道計算用の classId 順生データ (HNF と共役類の全メンバー σ) を
        // 戻り値で返す (classId = items の収集順 = GroupRelation.ConjugacyClassId。保存は呼び出し元)。
        return ([.. rels.OrderBy(r => r.Index).ThenBy(r => r.ChildSeriesNumber < 0 ? 1 : 0).ThenBy(r => r.ChildSeriesNumber)],
                [.. items.Select(it => (it.Hnf, it.Cls))]);
    }

    /// <summary>T′(fineHnf) ⊆ T″(coarseHnf) か (fineHnf の列が coarseHnf の列の整数combinationで表せるか)。260705Cl 追加。</summary>
    private static bool IsLatticeSubset(int[] coarseHnf, int[] fineHnf)
    {
        var coarseInv = RationalMatrix.Invert3(RationalMatrix.FromInt(coarseHnf));
        if (coarseInv == null) return false;
        return RationalMatrix.ToIntOrNull(RationalMatrix.Mul(coarseInv, RationalMatrix.FromInt(fineHnf))) != null;
    }

    /// <summary>整数ベクトル v (primitive 座標) が Q=G/hnf のどの coset に属するかを返す。260705Cl 追加。</summary>
    private static int CosetIndexOf(QuotientData q, int[] hnf, long[] v)
    {
        var hInv = RationalMatrix.Invert3(RationalMatrix.FromInt(hnf)) ?? throw new InvalidOperationException("HNF is singular");
        var label = RationalMatrix.ModVec1(RationalMatrix.MulVec(hInv, [new Fraction(v[0]), new Fraction(v[1]), new Fraction(v[2])]));
        for (int t = 0; t < q.Labels.Length; t++)
            if (RationalMatrix.VecEquals(q.Labels[t], label)) return t;
        throw new InvalidOperationException("coset representative not found");
    }

    // 260708Cl: enantiomorphic 対テーブルの複製を廃止し、SymmetryProperties.GetEnantiomorphPartnerNumber
    // (既存の掌性対 11 組テーブル) へ一本化 (二重保守防止、/simplify レビュー指摘)。
    ///// <summary>ITA の enantiomorphic 対 (11 対、双方向)。同型 (IIc) 判定は「同一タイプまたは enantiomorphic 対」。260708Cl 追加。</summary>
    //private static readonly Dictionary<int, int> _enantiomorphicPair = new()
    //{
    //    { 76, 78 }, { 78, 76 }, { 91, 95 }, { 95, 91 }, { 92, 96 }, { 96, 92 }, { 144, 145 }, { 145, 144 },
    //    { 151, 153 }, { 153, 151 }, { 152, 154 }, { 154, 152 }, { 169, 170 }, { 170, 169 }, { 171, 172 }, { 172, 171 },
    //    { 178, 179 }, { 179, 178 }, { 180, 181 }, { 181, 180 }, { 212, 213 }, { 213, 212 },
    //};

    /// <summary>1 つの complement (T′=hnf, σ=sigma) から GroupRelation (Kind=K または Isomorphic) を構築する。260705Cl 追加。
    /// 260708Cl: 呼び出し元 (ComputeMaximalK) が memo 化した QuotientData を q で渡せる (null なら構築)。</summary>
    private static GroupRelation BuildGroupRelation(int sn, PointGroupData pg, int[] hnf, int[] sigma, int classId, int conjugateCount, QuotientData q = null)
    {
        q ??= BuildQuotient(pg, hnf); // 260708Cl (旧: var q = BuildQuotient(pg, hnf);)
        int m = pg.LinKeys.Length;
        var basis = GetPrimitiveBasis(sn);
        var (child, pFrac, shiftFrac) = IdentifyK(sn, hnf, sigma, q);

        var baseOpByLin = new SymmetryOperation?[m];
        foreach (var op in TSubgroupFinder.GetExpandedOps(sn))
        {
            var key = LinKeyOf(op);
            int idx = -1;
            for (int i = 0; i < m; i++) if (SameIntVec(pg.LinKeys[i], key)) { idx = i; break; }
            if (idx >= 0 && baseOpByLin[idx] == null) baseOpByLin[idx] = op;
        }

        var reps = new SymmetryOperation[m];
        for (int i = 0; i < m; i++)
        {
            var repI = q.Reps[sigma[i]];
            var targetPrim = RationalMatrix.AddVec(pg.T0[i], [new Fraction(repI[0]), new Fraction(repI[1]), new Fraction(repI[2])]);
            var targetConv = RationalMatrix.MulVec(basis, targetPrim);
            var baseOp = baseOpByLin[i].Value;
            var baseT = baseOp.SeitzTranslation;
            reps[i] = new SymmetryOperation(baseOp, sn, ToDouble(targetConv[0]) - baseT.U, ToDouble(targetConv[1]) - baseT.V, ToDouble(targetConv[2]) - baseT.W);
        }

        // CosetRepresentatives: T/T′ の非自明な coset (恒等 coset=T′自身を除く) を代表する純並進操作。
        var identityOp = baseOpByLin[pg.E].Value;
        var identityT = identityOp.SeitzTranslation;
        var cosetReps = new List<SymmetryOperation>();
        for (int c = 1; c < q.N; c++)
        {
            var targetConv = RationalMatrix.MulVec(basis, [new Fraction(q.Reps[c][0]), new Fraction(q.Reps[c][1]), new Fraction(q.Reps[c][2])]);
            cosetReps.Add(new SymmetryOperation(identityOp, sn, ToDouble(targetConv[0]) - identityT.U, ToDouble(targetConv[1]) - identityT.V, ToDouble(targetConv[2]) - identityT.W));
        }

        double[] pDouble = null, shiftDouble = null;
        if (pFrac != null)
        {
            pDouble = [.. pFrac.Select(ToDouble)];
            shiftDouble = [.. shiftFrac.Select(ToDouble)];
        }
        string pointGroupHM = SymmetryStatic.Symmetries[sn].PointGroupHMStr switch { "2mm" or "m2m" => "mm2", var t2 => t2 };

        // 260708Cl: 同型 (ITA IIc) 判定 — 子タイプが親と同一または enantiomorphic 対なら Kind=Isomorphic。
        // 同型は klassengleiche の特殊例なので、データ構造とタブ表示ロジックは K と共通 (codex R7 合意)。
        int parentNo = SymmetryStatic.Symmetries[sn].SpaceGroupNumber;
        int childNo = child >= 0 ? SymmetryStatic.Symmetries[child].SpaceGroupNumber : -1;
        //bool isIso = childNo == parentNo || (_enantiomorphicPair.TryGetValue(parentNo, out int enPair) && enPair == childNo); // 260708Cl: テーブル一本化
        bool isIso = childNo == parentNo || (childNo > 0 && SymmetryProperties.GetEnantiomorphPartnerNumber(parentNo) == childNo);

        return new GroupRelation
        {
            Kind = isIso ? GroupRelationKind.Isomorphic : GroupRelationKind.K, // 260708Cl (旧: 常に K)
            ParentSeriesNumber = sn,
            Index = q.N,
            ConjugacyClassId = classId,
            ConjugateCount = conjugateCount,
            PointGroupHM = pointGroupHM,
            // 260705Cl: Operations は Representatives と同一 (親胞 mod1 の軌道生成・消滅則判定には未対応。
            // k- 専用の Orbit splitting / New reflections ロジックは Phase 2d 以降の課題、UI 側でガードする)。
            Operations = reps,
            Representatives = reps,
            CosetRepresentatives = [.. cosetReps],
            ChildSeriesNumber = child,
            TransformP = pDouble,
            TransformShift = shiftDouble,
            SublatticeBasis = [.. RationalMatrix.Mul(basis, RationalMatrix.FromInt(hnf)).Select(ToDouble)],
        };
    }

    private static double ToDouble(Fraction f) => (double)f.Num / (double)f.Den;
    #endregion
}
