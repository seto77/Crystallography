// 260805Cl 新規作成: 菊池線動力学化 Phase 0 (設計正本 = ReciPro_菊池線動力学化設計.md §2, §5, §7)。
// g/−g 面族ペアリングと帯幾何。設計の3分離 (幾何 = バンド幅 2θ_B / 動力学強度プロファイル / 表示線幅) のうち
// 「幾何」の静的部分を担う。強度は KikuchiProfileCalculator、表示は ReciPro 側 renderer。
// ⚠ この Diffraction/Kikuchi/ フォルダは WinForms / System.Drawing 非依存を規律で維持する (設計 §5)。

using System;
using System.Collections.Generic;

namespace Crystallography;

/// <summary>菊池バンドの g/−g 面族 (canonical 代表 = 先頭非零指数が正の側)。260805Cl 追加</summary>
public sealed class KikuchiBandFamily
{
    /// <summary>canonical 側の指数</summary>
    public (int H, int K, int L) Index { get; init; }

    /// <summary>canonical 側の逆格子ベクトル (結晶固定系, nm⁻¹)。回転は使用時に掛ける</summary>
    public Vector3DBase Vec { get; init; }

    /// <summary>|g| [nm⁻¹]。260806Cl /simplify: 導出値 (Vec と食い違い得る init フィールドを廃止)</summary>
    public double GLength => Vec.Length;

    /// <summary>ラベル (例 "1 1 1")</summary>
    public string Text { get; init; }

    /// <summary>運動学相対強度 (静的候補順位用)</summary>
    public double RelativeIntensity { get; init; }

    /// <summary>
    /// Vector3D 集合 (VectorOfG_KikuchiLine 相当) を g/−g 面族へペアリングする。
    /// canonical 代表は (h,k,l) の先頭非零成分が正の側。片符号しか含まれない入力でも canonical 形で1族にする。
    /// 返り値は入力の初出順を保存する (既存の静的選定順位を壊さない)。
    /// 注 (260806Cl): 現行の唯一の生産者 Crystal.SetVectorOfG_KikuchiLine は最初から canonical 半空間しか
    /// 列挙しないため、flip / dedup は防御的正規化 (他の入力源に備えたもの) であり通常は素通りする。
    /// </summary>
    public static List<KikuchiBandFamily> Pair(IEnumerable<Vector3D> vectors)
    {
        var seen = new HashSet<(int, int, int)>();
        var list = new List<KikuchiBandFamily>();
        foreach (var v in vectors)
        {
            var (h, k, l) = v.Index;
            if (h == 0 && k == 0 && l == 0)
                continue;
            bool flip = h < 0 || (h == 0 && (k < 0 || (k == 0 && l < 0)));
            var idx = flip ? (-h, -k, -l) : (h, k, l);
            if (!seen.Add(idx))
                continue;
            list.Add(new KikuchiBandFamily
            {
                Index = idx,
                Vec = flip ? -(Vector3DBase)v : new Vector3DBase(v.X, v.Y, v.Z), //260806Cl /simplify: 既存の単項マイナス演算子を使用
                Text = $"{idx.Item1} {idx.Item2} {idx.Item3}",
                RelativeIntensity = v.RelativeIntensity,
            });
        }
        return list;
    }

    /// <summary>
    /// 260807Cl 追加 (設計 Phase 2.5): systematic row 計算と組み合わせるときの候補整理。
    /// 方向 (h,k,l)/gcd が同じ族を 1 本の row と見なし、その中で |g| 最小の族 (= row の生成元) だけを残す。
    /// ただし生成元に対する倍率 n が xMax を超える高次は、生成元のプロファイル範囲 |x| ≤ xMax の外にあって
    /// 描かれないので**独立な族として残す**。
    ///
    /// ⚠これをしないと**同一 row の高次を二重に計上する**: 候補には {020} と {040} が別々に並ぶのに、
    /// {020} の row 計算 (N ≥ 2) は {040} を beam として既に含んでいる。
    /// 260807Cl 実測 (MgO Omoto Fig 6): 畳まずに row を入れると corr が 0.914 → 0.898 と**悪化**し、
    /// 畳むと 0.927 まで改善した。RowOrder &gt; 0 のときは ComputeProfiles が自動で適用する。
    /// </summary>
    public static List<KikuchiBandFamily> CollapseSystematicRows(IReadOnlyList<KikuchiBandFamily> families, double xMax)
    {
        static int Gcd(int a, int b) => b == 0 ? a : Gcd(b, a % b);
        static (int, int, int) Dir(in (int H, int K, int L) idx)
        {
            var g = Gcd(Gcd(Math.Abs(idx.H), Math.Abs(idx.K)), Math.Abs(idx.L));
            return g == 0 ? idx : (idx.H / g, idx.K / g, idx.L / g);
        }

        var generator = new Dictionary<(int, int, int), KikuchiBandFamily>();
        foreach (var f in families)
        {
            var d = Dir(f.Index);
            if (!generator.TryGetValue(d, out var cur) || f.GLength < cur.GLength)
                generator[d] = f;
        }
        var kept = new List<KikuchiBandFamily>(families.Count);
        foreach (var f in families)
        {
            var gen = generator[Dir(f.Index)];
            // 倍率 1 (生成元自身) か、生成元の描画域外に落ちる高次だけ残す
            if (ReferenceEquals(f, gen) || f.GLength > xMax * gen.GLength + 1e-9)
                kept.Add(f);
        }
        return kept;
    }
}
