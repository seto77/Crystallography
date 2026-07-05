// 260705Cl 新規: k-部分群エンジン (KSubgroupFinder, Phase 2c) の有理数演算基盤。
//
// t-エンジン (TSubgroupFinder) は double + 1/24 格子スナップで足りたが、k- は指数 2/3/4 の分数
// (1/2,1/3,1/4 等の組み合わせ) がスナップ許容に乗りにくく、HNF・格子一致検証・q_parent 判定を
// 有理数で厳密に行う必要がある (計画書 §3.2「有理数演算」、.project-guidance/ReciPro_FormGroupRelations改修計画.md §4.2)。
//
// 既約分数は BigInteger 分子分母 (long は HNF/Smith 標準形・isomorphic 系列展開で溢れ得るとの codex 指摘)。
// 3×3 行列・3 ベクトルは既存 TSubgroupFinder の int[9]/double[9] 配列スタイルに合わせ、Fraction[9]/Fraction[3]
// の軽量配列 + static helper で統一する (構造体化した汎用行列クラスは導入しない)。
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;

namespace Crystallography;

/// <summary>既約分数 (分子・分母は BigInteger、分母は常に正)。260705Cl 追加 (Phase 2c)。</summary>
public readonly struct Fraction : IEquatable<Fraction>
{
    public BigInteger Num { get; }
    public BigInteger Den { get; }

    public static readonly Fraction Zero = new(0);
    public static readonly Fraction One = new(1);

    public Fraction(BigInteger num, BigInteger den)
    {
        if (den.IsZero)
            throw new DivideByZeroException("Fraction denominator is zero");
        if (den.Sign < 0) { num = -num; den = -den; }
        var g = BigInteger.GreatestCommonDivisor(BigInteger.Abs(num), den);
        if (g > BigInteger.One) { num /= g; den /= g; }
        Num = num; Den = den;
    }

    /// <summary>整数 (分母 1)。default(Fraction) もこれ経由で 0/1 になる。</summary>
    public Fraction(BigInteger num) : this(num, BigInteger.One) { }

    public bool IsZero => Num.IsZero;
    public int Sign => Num.Sign;

    public static implicit operator Fraction(int n) => new(n);
    public static implicit operator Fraction(BigInteger n) => new(n);

    public static Fraction operator +(Fraction a, Fraction b) => new(a.Num * b.Den + b.Num * a.Den, a.Den * b.Den);
    public static Fraction operator -(Fraction a, Fraction b) => new(a.Num * b.Den - b.Num * a.Den, a.Den * b.Den);
    public static Fraction operator -(Fraction a) => new(-a.Num, a.Den);
    public static Fraction operator *(Fraction a, Fraction b) => new(a.Num * b.Num, a.Den * b.Den);
    public static Fraction operator /(Fraction a, Fraction b)
    {
        if (b.IsZero) throw new DivideByZeroException();
        return new Fraction(a.Num * b.Den, a.Den * b.Num);
    }
    public static bool operator ==(Fraction a, Fraction b) => a.Num == b.Num && a.Den == b.Den;
    public static bool operator !=(Fraction a, Fraction b) => !(a == b);
    public static bool operator <(Fraction a, Fraction b) => (a - b).Sign < 0;
    public static bool operator >(Fraction a, Fraction b) => (a - b).Sign > 0;
    public static bool operator <=(Fraction a, Fraction b) => (a - b).Sign <= 0;
    public static bool operator >=(Fraction a, Fraction b) => (a - b).Sign >= 0;

    /// <summary>[0,1) への正規化 (floor 除算。負数も exact、double を経由しない)。</summary>
    public Fraction Mod1()
    {
        var r = BigInteger.Remainder(Num, Den);
        if (r.Sign < 0) r += Den;
        return new Fraction(r, Den);
    }

    /// <summary>整数なら true (Mod1 したうえで判定するのが通例の呼び方)。</summary>
    public bool IsInteger => Den.IsOne;

    public bool Equals(Fraction other) => this == other;
    public override bool Equals(object obj) => obj is Fraction f && Equals(f);
    public override int GetHashCode() => HashCode.Combine(Num, Den);
    public override string ToString() => Den.IsOne ? Num.ToString() : $"{Num}/{Den}";

    /// <summary>既存の SymmetryOperation 由来 double を厳密有理数へ変換する唯一の入口。
    /// maxDen 以下の分母で一致しなければ例外 (内部計算では二度と double へ戻さない・再スナップしない)。</summary>
    public static Fraction FromDouble(double d, int maxDen = 96)
    {
        for (int den = 1; den <= maxDen; den++)
        {
            double x = d * den;
            var r = (BigInteger)Math.Round(x);
            if (Math.Abs(x - (double)r) < 1e-6)
                return new Fraction(r, den);
        }
        throw new InvalidOperationException($"cannot rationalize {d} with denominator <= {maxDen}");
    }
}

/// <summary>有理数の 3×3 行列 (row-major 9 要素) / 3-ベクトルの static helper。260705Cl 追加 (Phase 2c)。
/// 既存 TSubgroupFinder の int[9] row-major 整数行列と同じ添字規約 (m[r*3+c])。</summary>
public static class RationalMatrix
{
    public static Fraction[] Mul(Fraction[] a, Fraction[] b)
    {
        var c = new Fraction[9];
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                c[i * 3 + j] = a[i * 3] * b[j] + a[i * 3 + 1] * b[3 + j] + a[i * 3 + 2] * b[6 + j];
        return c;
    }

    public static Fraction[] MulVec(Fraction[] m, Fraction[] v)
    {
        var r = new Fraction[3];
        for (int i = 0; i < 3; i++)
            r[i] = m[i * 3] * v[0] + m[i * 3 + 1] * v[1] + m[i * 3 + 2] * v[2];
        return r;
    }

    public static Fraction[] AddVec(Fraction[] a, Fraction[] b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
    public static Fraction[] SubVec(Fraction[] a, Fraction[] b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];

    public static Fraction[] Transpose(Fraction[] m) => [m[0], m[3], m[6], m[1], m[4], m[7], m[2], m[5], m[8]];

    public static Fraction Det3(Fraction[] m)
        => m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6]) + m[2] * (m[3] * m[7] - m[4] * m[6]);

    /// <summary>逆行列 (余因子 / det)。特異なら null。</summary>
    public static Fraction[] Invert3(Fraction[] m)
    {
        var det = Det3(m);
        if (det.IsZero) return null;
        Fraction[] adjT =
        [
            m[4] * m[8] - m[5] * m[7], m[2] * m[7] - m[1] * m[8], m[1] * m[5] - m[2] * m[4],
            m[5] * m[6] - m[3] * m[8], m[0] * m[8] - m[2] * m[6], m[2] * m[3] - m[0] * m[5],
            m[3] * m[7] - m[4] * m[6], m[1] * m[6] - m[0] * m[7], m[0] * m[4] - m[1] * m[3],
        ];
        var inv = new Fraction[9];
        for (int i = 0; i < 9; i++) inv[i] = adjT[i] / det;
        return inv;
    }

    public static Fraction[] FromInt(int[] m)
    {
        var r = new Fraction[9];
        for (int i = 0; i < 9; i++) r[i] = m[i];
        return r;
    }

    /// <summary>各成分が整数ならその int[9] を返す。1 つでも非整数なら null。</summary>
    public static int[] ToIntOrNull(Fraction[] m)
    {
        var r = new int[9];
        for (int i = 0; i < 9; i++)
        {
            if (!m[i].IsInteger) return null;
            r[i] = (int)m[i].Num;
        }
        return r;
    }

    public static Fraction[] ModVec1(Fraction[] v) => [v[0].Mod1(), v[1].Mod1(), v[2].Mod1()];

    public static bool VecEquals(Fraction[] a, Fraction[] b) => a[0] == b[0] && a[1] == b[1] && a[2] == b[2];
}

/// <summary>整数格子の小規模ユーティリティ (行ベクトル規約)。260705Cl 追加 (Phase 2c、k-エンジンの primitive 基底構築用)。
/// 値は常に小さい (中心化ベクトルの分母を掛けた程度) ため long で十分、BigInteger は使わない。</summary>
public static class IntegerLattice
{
    /// <summary>拡張ユークリッド互除法。a*x + b*y = gcd(|a|,|b|) を満たす (gcd, x, y) を返す (gcd ≥ 0)。
    /// 260705Cl 修正: 初期値 (x0,x1)/(y0,y1) が標準アルゴリズムと入れ替わっており a*x+b*y=gcd を満たさない
    /// 実バグがあった (例: ExtGcd(2,1) が (1,1,0) を返し 2*1+1*0=2≠1 になっていた。primitive 基底構築で
    /// 中心化ベクトルが実質無視され F/I 中心化を取り逃がす原因)。標準の old_r/old_s/old_t 初期値に修正。</summary>
    public static (long Gcd, long X, long Y) ExtGcd(long a, long b)
    {
        long r0 = a, r1 = b;
        long x0 = 1, x1 = 0;
        long y0 = 0, y1 = 1;
        while (r1 != 0)
        {
            long q = r0 / r1, r2 = r0 % r1;
            (r0, r1) = (r1, r2);
            (x0, x1) = (x1, x0 - q * x1);
            (y0, y1) = (y1, y0 - q * y1);
        }
        if (r0 < 0) { r0 = -r0; x0 = -x0; y0 = -y0; }
        return (r0, x0, y0);
    }

    /// <summary>整数生成元集合 (行ベクトル、3 要素) が張る格子の基底 (3 本、行ベクトル) を抽出する。
    /// 列ごとに gcd 結合で 1 本に絞り込む標準的な手法 (行 HNF の簡易版、正準形までは求めない)。
    /// 生成元がランク 3 を張らない (退化) 場合は null。</summary>
    public static long[][] BasisFromGenerators(IReadOnlyList<long[]> generators)
    {
        var pool = generators.Select(g => (long[])g.Clone()).ToList();
        var basis = new long[3][];
        for (int col = 0; col < 3; col++)
        {
            while (true)
            {
                int i = pool.FindIndex(r => r[col] != 0);
                if (i < 0) break;
                int j = pool.FindIndex(i + 1, r => r[col] != 0);
                if (j < 0) break; // ちょうど 1 本だけ非零 = この列の pivot 候補
                var (g, x, y) = ExtGcd(pool[i][col], pool[j][col]);
                long p = pool[i][col] / g, q = pool[j][col] / g;
                var newI = new long[3];
                var newJ = new long[3];
                for (int c = 0; c < 3; c++)
                {
                    newI[c] = x * pool[i][c] + y * pool[j][c];
                    newJ[c] = q * pool[i][c] - p * pool[j][c];
                }
                pool[i] = newI; pool[j] = newJ;
            }
            int pivot = pool.FindIndex(r => r[col] != 0);
            if (pivot < 0) return null; // ランク不足 (この使い方では発生しない想定)
            basis[col] = pool[pivot];
            pool.RemoveAt(pivot);
        }
        return basis;
    }
}
