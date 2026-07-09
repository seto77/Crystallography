// 260704Cl 追加: 対称操作 (SymmetryOperation) を「座標トリプレット / Seitz 記号 / 幾何的解釈」の
// 3 表現へ文字列化するヘルパ。FormSymmetryInformation の Operations タブ (Phase 1-D) と CIF コピーで使用。
// 全て操作の線形部 (Order/Sense/Direction) と SeitzTranslation から算出し、Cartesian 変換に依存しない。
// 参照: ReciPro_SymmetryInformation拡張計画.md §4.1。
using System;
using System.Collections.Concurrent; // 260709Cl 追加: series ごとの中心化並進キャッシュ
using System.Collections.Generic; // 260709Cl 追加
using System.Globalization;
using System.Linq; // 260709Cl 追加
using System.Text;

namespace Crystallography;

/// <summary>対称操作を ITA / CIF 流の 3 表現へ整形する静的ヘルパ (260704Cl 追加)。</summary>
public static class SeitzNotation
{
    // 並進成分をスナップする分母候補 (螺旋/映進で現れるのは 1/2, 1/3, 1/4, 1/6, 1/8 など)。
    // 260705Cl: 先頭にあった 1 は SnapFraction の条件 (0 < n < den) を満たせない死にエントリだったため除去。
    private static readonly int[] Denominators = [2, 3, 4, 6, 8, 12];
    private const double Tol = 1e-6;

    #region 線形部行列の取り出し
    /// <summary>操作の線形部 3×3 整数行列 R を取り出す。R[i,j] は「出力成分 i に対する入力 j の係数」。
    /// 六方晶系では係数が -2..2 になり得る (x−y など) が全て整数。</summary>
    public static int[,] LinearMatrix(in SymmetryOperation op)
    {
        var cx = op.ApplyMatrix(1, 0, 0); // = R の第 0 列 (x 基底の像)
        var cy = op.ApplyMatrix(0, 1, 0);
        var cz = op.ApplyMatrix(0, 0, 1);
        return new int[,]
        {
            { Ri(cx.X), Ri(cy.X), Ri(cz.X) },
            { Ri(cx.Y), Ri(cy.Y), Ri(cz.Y) },
            { Ri(cx.Z), Ri(cy.Z), Ri(cz.Z) },
        };
    }

    private static int Ri(double d) => (int)Math.Round(d);
    #endregion

    #region 座標トリプレット (例 "-y, x-y, z+1/3")
    /// <summary>操作を座標トリプレット文字列 "…, …, …" に整形する。</summary>
    // 260705Cl: 全呼び出しで未使用だった投機的パラメータ extraTranslation を削除 (ops は中心化展開済みで不要)。
    //public static string Triplet(in SymmetryOperation op, (double U, double V, double W) extraTranslation = default)
    public static string Triplet(in SymmetryOperation op)
    {
        var R = LinearMatrix(op);
        var t = op.SeitzTranslation;
        double[] tr = [t.U, t.V, t.W];
        var sb = new StringBuilder(24);
        for (int i = 0; i < 3; i++)
        {
            if (i > 0) sb.Append(", ");
            sb.Append(RowExpression(R[i, 0], R[i, 1], R[i, 2], tr[i]));
        }
        return sb.ToString();
    }

    /// <summary>1 出力成分ぶんの式 (係数×x,y,z + 並進) を作る。空なら "0"。</summary>
    private static string RowExpression(int cx, int cy, int cz, double t)
    {
        var sb = new StringBuilder(8);
        AppendVar(sb, cx, 'x');
        AppendVar(sb, cy, 'y');
        AppendVar(sb, cz, 'z');
        AppendFraction(sb, t);
        return sb.Length == 0 ? "0" : sb.ToString();
    }

    private static void AppendVar(StringBuilder sb, int coef, char v)
    {
        if (coef == 0) return;
        bool first = sb.Length == 0;
        if (coef > 0 && !first) sb.Append('+');
        if (coef == -1) sb.Append('-');
        else if (coef == 1) { /* 係数 1 は省略 */ }
        else sb.Append(coef.ToString(CultureInfo.InvariantCulture)); // -2x など
        sb.Append(v);
    }

    private static void AppendFraction(StringBuilder sb, double t)
    {
        var (num, den) = SnapFraction(t);
        if (num == 0) return;
        if (num > 0 && sb.Length > 0) sb.Append('+');
        else if (num < 0) sb.Append('-');
        int a = Math.Abs(num);
        sb.Append(den == 1 ? a.ToString(CultureInfo.InvariantCulture)
                           : $"{a}/{den}"); // 260705Cl: string への no-op ToString を除去
    }

    /// <summary>並進 t を [0,1) に落として最も近い num/den (den は既定候補) に丸める。</summary>
    private static (int Num, int Den) SnapFraction(double t)
    {
        // [0,1) へ正規化
        t -= Math.Floor(t + Tol);
        if (t < Tol || t > 1 - Tol) return (0, 1);
        foreach (var den in Denominators)
        {
            double x = t * den;
            int n = (int)Math.Round(x);
            if (Math.Abs(x - n) < 1e-4 && n > 0 && n < den)
            {
                int g = (int)GammaFunction.Gcd(n, den); // 260705Cl: 既存 GammaFunction.Gcd に一本化
                return (n / g, den / g);
            }
        }
        // フォールバック: 12 分割で最近似
        int nn = (int)Math.Round(t * 12);
        int gg = (int)GammaFunction.Gcd(nn, 12); // 260705Cl: 同上 (nn=0 でも Gcd(0,12)=12 で従来と同値)
        return (nn / gg, 12 / gg);
    }

    // 260705Cl: 私製 Gcd を削除し既存 GammaFunction.Gcd (Mathematics/GammaFunction.cs) に一本化。
    //private static int Gcd(int a, int b) { a = Math.Abs(a); b = Math.Abs(b); while (b != 0) (a, b) = (b, a % b); return a == 0 ? 1 : a; }
    #endregion

    #region Seitz 記号 (例 "3+ [111]", "m [1-10]", "-1")
    /// <summary>ITA 流の簡易 Seitz 記号 (回転記号 + 向き [uvw] + 並進があれば付記)。</summary>
    public static string Seitz(in SymmetryOperation op)
    {
        int order = op.Order;
        var t = op.SeitzTranslation; // 260708Ch: 中心化並進を 1/-1 でも落とさない
        string trans = HasTranslation(t) ? $" {FractionTriplet(t)}" : ""; // 260708Ch
        //if (order == 1) return "1"; // 260708Ch: F/I/A/B/C/R 格子の中心化並進が消えていた
        //if (order == -1) return "-1"; // 260708Ch
        if (order == 1) return $"1{trans}"; // 260708Ch
        if (order == -1) return $"-1{trans}"; // 260708Ch

        string rot = RotationSymbol(order, op.Sense); // 例 "3+", "m", "-4+"
        string dir = DirectionStr(op.Direction);
        return $"{rot} {dir}{trans}".TrimEnd();
    }

    /// <summary>回転/回反の記号。order は ±2,±3,±4,±6 (負=回反, -2 は鏡映 m)。</summary>
    private static string RotationSymbol(int order, bool sense)
    {
        if (order == -2) return "m";
        int n = Math.Abs(order);
        string sign = order < 0 ? "-" : "";
        string s = n >= 3 ? (sense ? "+" : "-") : ""; // 2 回は向きの区別なし
        return $"{sign}{n}{s}";
    }

    /// <summary>ITA 流の Seitz 記号を LaTeX の {R|t} 記法へ整形する (FormSymmetryInformation の Operations タブ用)。
    /// 260708Ch 追加: 旧実装 (FormSymmetryInformation.SeitzToLatex) は Seitz() が返す文字列を正規表現で再パースして
    /// LaTeX へ組み直していたが、Seitz() と同じ構造データ (Order/Sense/Direction/SeitzTranslation) から直接組み立てる形に刷新。</summary>
    public static string SeitzLatex(in SymmetryOperation op)
    {
        int order = op.Order;
        var t = op.SeitzTranslation;
        string trans = (HasTranslation(t) ? FractionTriplet(t) : "0,0,0").Replace(",", @",\,");

        if (order == 1) return $@"\{{\,1\mid {trans}\,\}}";
        if (order == -1) return $@"\{{\,\bar{{1}}\mid {trans}\,\}}";

        string rot = RotationSymbolLatex(order, op.Sense);
        string dir = $"_{{{DirectionLatex(op.Direction)}}}";
        return $@"\{{\,{rot}{dir}\mid {trans}\,\}}";
    }

    /// <summary>RotationSymbol の LaTeX 版。回反 (負の order) は \bar{n}、n≧3 の向き (+/-) は上付きにする。260708Ch 追加。</summary>
    private static string RotationSymbolLatex(int order, bool sense)
    {
        if (order == -2) return "m";
        int n = Math.Abs(order);
        string body = order < 0 ? $@"\bar{{{n}}}" : I(n);
        string s = n >= 3 ? (sense ? "^{+}" : "^{-}") : "";
        return $"{body}{s}";
    }

    /// <summary>DirectionStr の LaTeX 版。負の成分を \bar{} で包む。文字列を正規表現で再分割せず (int U,V,W)
    /// から直接組み立てるため、DirectionStr の桁連結表記 (例 "[1-10]") が抱える複数桁成分の曖昧さを引きずらない。260708Ch 追加。</summary>
    private static string DirectionLatex((int U, int V, int W) d)
        => $"{Bar(d.U)}{Bar(d.V)}{Bar(d.W)}";

    private static string Bar(int v) => v < 0 ? $@"\bar{{{I(-v)}}}" : I(v);
    #endregion

    #region 幾何的解釈 (例 "3-fold rotation [111]", "c-glide ⊥[001]")
    /// <summary>操作の幾何的な種類 (回転/螺旋/鏡映/映進/回反/反転) を短い英語句で返す。多言語化は呼び出し側で不要
    /// (記号中心の表示のため据置)。位置情報は簡潔さのため向き [uvw] のみ添える。</summary>
    public static string GeometricType(in SymmetryOperation op)
    {
        int order = op.Order;
        if (order == 1) // 260708Ch: 中心化並進つき {1|t} は Identity ではなく純並進として表示
        {
            var t = op.SeitzTranslation; // 260708Ch
            return HasTranslation(t) ? $"Translation {FractionTriplet(t)}" : "Identity"; // 260708Ch
        }
        //if (order == -1) return $"Inversion centre at {PointStr(op.Position)}"; // 260708Ch: {-1|t} の中心は t/2
        if (order == -1) // 260708Ch
        {
            var t = op.SeitzTranslation; // 260708Ch
            return $"Inversion centre at {PointStr((t.U / 2, t.V / 2, t.W / 2))}"; // 260708Ch
        }

        var dir = op.Direction;
        int n = Math.Abs(order);

        if (order == -2)
        {
            // 鏡映または映進。IntrinsicTranslation が面内並進成分。
            var g = op.IntrinsicTranslation;
            string letter = GlideLetter(g, dir);
            string kind = letter == "m" ? "Mirror plane m" : $"{letter}-glide plane";
            return $"{kind} ⊥{DirectionStr(dir)}";
        }

        if (order < 0) // -3, -4, -6 回反
            return $"{n}-fold rotoinversion (-{n}{(op.Sense ? "+" : "-")}) {DirectionStr(dir)}";

        // 正の回転 or 螺旋
        //int pitch = ScrewPitch(op.IntrinsicTranslation, dir, n);
        int pitch = ScrewPitch(op.IntrinsicTranslation, dir, n, op.SeriesNumber); // 260709Cl: 中心化補正 (下記)
        if (pitch > 0)
            return $"{ScrewLabel(n, pitch)} screw axis {DirectionStr(dir)}";
        return $"{n}-fold rotation ({n}{(op.Sense ? "+" : "-")}) {DirectionStr(dir)}";
    }

    /// <summary>螺旋のピッチ p (n_p の p)。intrinsic 並進のうち軸方向成分 ≈ p/n。無ければ 0。
    /// 260709Cl シグネチャ変更 (旧: ScrewPitch(it, dir, n)): 中心化格子 (I/F/R) の対角軸では軸方向の最小
    /// 格子並進 primitive_along_d が方向ベクトル d より短い (I-立方の体対角 = d/2 等) ため、ITA 規約
    /// (axial = (p/n)·primitive_along_d) に合わせ seriesNumber の中心化並進で補正する。全 530 設定の検証
    /// (SymmetryPropsCheck PART 9) で I 格子立方晶の [111] 3 回軸 200 操作が誤添字だった実バグの修正
    /// (例: I23 の {3⁺|½½½} を 3₂ 螺旋と表示 — ½½½ は I 中心化並進そのもので、実際は純回転)。
    /// seriesNumber が無効 (0) のときは中心化なし (P 格子) として従来と同じ値を返す。</summary>
    private static int ScrewPitch((double U, double V, double W) it, (int U, int V, int W) dir, int n, int seriesNumber)
    {
        double len2 = (double)dir.U * dir.U + (double)dir.V * dir.V + (double)dir.W * dir.W;
        if (len2 < Tol) return 0;
        double proj = (it.U * dir.U + it.V * dir.V + it.W * dir.W) / len2; // 軸方向の並進 (格子単位)
        //int p = (int)Math.Round(proj * n);
        //p = ((p % n) + n) % n;
        //return p;
        double primitive = SymmetryElementsTable.PrimitiveAlongDirectionInDUnits(dir, CenteringsOf(seriesNumber)); // 260709Cl
        double alongPrim = proj / primitive; // 260709Cl: d 単位 → primitive_along_d 単位
        alongPrim -= Math.Floor(alongPrim);
        if (alongPrim < 1e-3 || alongPrim > 1 - 1e-3) return 0; // 格子並進と同値 = 純回転
        int p = ((int)Math.Round(alongPrim * n)) % n;
        return p < 0 ? p + n : p;
    }

    /// <summary>260709Cl 追加: series の中心化並進 (恒等線形部の非ゼロ mod1 並進)。ScrewPitch の
    /// primitive_along_d 補正用 (series ごとにキャッシュ、複数スレッド安全)。無効 series は空 (P 格子扱い)。</summary>
    private static readonly ConcurrentDictionary<int, (double U, double V, double W)[]> _centeringCache = new();
    private static (double U, double V, double W)[] CenteringsOf(int seriesNumber)
    {
        if (seriesNumber <= 0 || seriesNumber >= SymmetryStatic.TotalSpaceGroupNumber)
            return [];
        return _centeringCache.GetOrAdd(seriesNumber, sn =>
        {
            var list = new List<(double U, double V, double W)>();
            foreach (var op in TSubgroupFinder.GetExpandedOps(sn))
            {
                if (op.Order != 1) continue;
                var t = op.SeitzTranslation;
                double cu = Frac(t.U), cv = Frac(t.V), cw = Frac(t.W);
                if (cu + cv + cw < Tol) continue;
                if (!list.Any(c => Math.Abs(c.U - cu) + Math.Abs(c.V - cv) + Math.Abs(c.W - cw) < 1e-6))
                    list.Add((cu, cv, cw));
            }
            return [.. list];
        });
    }

    private static string ScrewLabel(int n, int p) => $"{n}{Subscript(p)}";

    private static string Subscript(int p) => p switch
    {
        0 => "₀", 1 => "₁", 2 => "₂", 3 => "₃",
        4 => "₄", 5 => "₅", _ => p.ToString(CultureInfo.InvariantCulture)
    };

    /// <summary>映進面の並進成分から glide 記号 (m/a/b/c/n/d) を推定する。</summary>
    private static string GlideLetter((double U, double V, double W) g, (int U, int V, int W) normal)
    {
        double gu = Frac(g.U), gv = Frac(g.V), gw = Frac(g.W);
        int nnz = (gu > Tol ? 1 : 0) + (gv > Tol ? 1 : 0) + (gw > Tol ? 1 : 0);
        if (nnz == 0) return "m";
        // 1/4 成分を含めば d、単一 1/2 成分なら軸名、2 成分なら n。
        bool anyQuarter = IsNear(gu, 0.25) || IsNear(gu, 0.75) || IsNear(gv, 0.25) || IsNear(gv, 0.75) || IsNear(gw, 0.25) || IsNear(gw, 0.75);
        if (anyQuarter) return "d";
        if (nnz >= 2) return "n";
        if (IsNear(gu, 0.5)) return "a";
        if (IsNear(gv, 0.5)) return "b";
        if (IsNear(gw, 0.5)) return "c";
        return "g"; // その他 (e-glide などは簡易化のため g)
    }

    private static double Frac(double d) { d -= Math.Floor(d); return d > 1 - Tol ? 0 : d; }
    private static bool IsNear(double d, double v) => Math.Abs(Frac(d) - v) < 1e-3;
    #endregion

    #region CIF ループ
    /// <summary>全対称操作 (中心化展開済み ops) を CIF の _space_group_symop_operation_xyz ループ文字列にする。</summary>
    public static string ToCifSymopLoop(System.Collections.Generic.IReadOnlyList<SymmetryOperation> ops)
    {
        var sb = new StringBuilder(64 + ops.Count * 24);
        sb.AppendLine("loop_");
        sb.AppendLine("_space_group_symop_id");
        sb.AppendLine("_space_group_symop_operation_xyz");
        for (int i = 0; i < ops.Count; i++)
        {
            var op = ops[i];
            sb.Append((i + 1).ToString(CultureInfo.InvariantCulture)).Append(" '")
              .Append(Triplet(op)).Append('\'').Append('\n');
        }
        return sb.ToString();
    }
    #endregion

    #region 共通の小整形
    private static bool HasTranslation((double U, double V, double W) t)
        => Frac(t.U) > Tol || Frac(t.V) > Tol || Frac(t.W) > Tol;

    /// <summary>並進ベクトルを "1/2,0,1/2" 形式に。</summary>
    private static string FractionTriplet((double U, double V, double W) t)
        => $"{Frac1(t.U)},{Frac1(t.V)},{Frac1(t.W)}";

    private static string Frac1(double d)
    {
        var (num, den) = SnapFraction(d);
        return num == 0 ? "0" : (den == 1 ? num.ToString(CultureInfo.InvariantCulture) : $"{num}/{den}");
    }

    /// <summary>点 (fractional) を "1/4,0,1/4" 形式の分数表記で簡潔に。(260705Cl: 実装と食い違っていた説明を修正)</summary>
    private static string PointStr((double U, double V, double W) p)
        => $"{Frac1(p.U)},{Frac1(p.V)},{Frac1(p.W)}";

    /// <summary>整数方向 (uvw) を "[uvw]" 形式に (負は連結: [1-10])。</summary>
    public static string DirectionStr((int U, int V, int W) d) => $"[{I(d.U)}{I(d.V)}{I(d.W)}]";

    private static string I(int v) => v.ToString(CultureInfo.InvariantCulture);

    /// <summary>HM 文字列内の螺旋添字 "sub1".."sub6" を Unicode 下付き文字へ (表セル等の簡易整形)。
    /// 260705Cl 追加: FormSymmetryInformation / FormGroupRelations に重複していた同一実装をここへ集約。</summary>
    public static string PrettyHM(string hm)
    {
        if (string.IsNullOrEmpty(hm)) return hm;
        var sb = new StringBuilder(hm);
        sb.Replace("sub1", "₁").Replace("sub2", "₂").Replace("sub3", "₃")
          .Replace("sub4", "₄").Replace("sub5", "₅").Replace("sub6", "₆");
        return sb.ToString();
    }
    #endregion
}
