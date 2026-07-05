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
}
