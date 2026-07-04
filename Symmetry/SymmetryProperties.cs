// 260704Cl 追加: 空間群 1 つの群論的性質を算出する。FormSymmetryInformation の Properties タブ (Phase 1-E/G)。
// 反転心・掌性(Sohncke)・極性(方向)・物性(焦電/圧電/SHG/旋光)の許容を「点群操作の整数線形部」から計算し、
// symmorphic / enantiomorphic / arithmetic crystal class / Patterson 等は小テーブルで補う。Cartesian 変換に非依存。
// 物性の許容則は Nye "Physical Properties of Crystals" 準拠 (10 polar / 20 piezoelectric / 15 gyrotropic)。
// 参照: ReciPro_SymmetryInformation拡張計画.md §4.2, §4.3。
using System;
using System.Collections.Generic;
using System.Linq;

namespace Crystallography;

/// <summary>空間群の群論的性質・物性許容を保持する (260704Cl 追加)。1 結晶切替につき 1 度生成する。</summary>
public sealed class SymmetryProperties
{
    #region 公開プロパティ
    /// <summary>一般位置の多重度 (= 空間群の対称操作数, 中心化を含む)。</summary>
    public int GeneralMultiplicity { get; }
    /// <summary>点群の位数 (相異なる線形部の数)。</summary>
    public int PointGroupOrder { get; }

    /// <summary>反転心を持つか (中心対称)。</summary>
    public bool IsCentrosymmetric { get; }
    /// <summary>Sohncke 群か (回反・鏡映・反転を含まず、真回転のみ = キラル環境で許される)。</summary>
    public bool IsSohncke { get; }
    /// <summary>Symmorphic 群か (73 群)。</summary>
    public bool IsSymmorphic { get; }

    /// <summary>極性群か (全操作で不変な方向が存在)。</summary>
    public bool IsPolar { get; }
    /// <summary>極性方向の説明 ("[001]" / "any" / "any in (001)" / "none")。</summary>
    public string PolarDirectionStr { get; }

    /// <summary>掌性対 (エナンチオモルフ) を持つか。</summary>
    public bool HasEnantiomorph { get; }
    /// <summary>掌性対の相手の空間群 IT 番号 (無ければ 0)。</summary>
    public int EnantiomorphPartnerNumber { get; }

    public string CrystalFamilyStr { get; }
    public string LatticeSystemStr { get; }
    public string BravaisTypeStr { get; }
    public string ArithmeticCrystalClassStr { get; }
    public string PattersonSymmetryStr { get; }

    // --- 物性の対称性許容 (点群で許容されるか。存在の主張ではない) ---
    /// <summary>焦電性 (rank-1 極性ベクトル) が許容されるか = 極性群。</summary>
    public bool PyroelectricAllowed { get; }
    /// <summary>圧電性 (rank-3 極性テンソル) が許容されるか = 非中心対称かつ点群 ≠ 432。</summary>
    public bool PiezoelectricAllowed { get; }
    /// <summary>第二高調波発生 (SHG, χ⁽²⁾ rank-3) が許容されるか = 圧電性と同条件。</summary>
    public bool SHGAllowed { get; }
    /// <summary>旋光性 (自然光学活性, rank-2 軸性テンソル) が許容されるか (15 gyrotropic 点群)。</summary>
    public bool OpticalActivityAllowed { get; }
    #endregion

    #region コンストラクタ (点群操作から算出)
    public SymmetryProperties(in Symmetry sym)
    {
        int sn = sym.SeriesNumber;

        // 一般位置の対称操作 (中心化展開済み)。SeriesNumber を付け直して hex 系の ApplyMatrix を正しく効かせる。
        // 260705Cl: 展開処理を TSubgroupFinder.GetExpandedOps に一本化 (4 箇所に散在していた同型ブロックの解消)。
        //var wyck = SymmetryStatic.WyckoffPositions[sn][0];
        //var ops = wyck.PositionOperations;
        //GeneralMultiplicity = wyck.Multiplicity;
        GeneralMultiplicity = SymmetryStatic.WyckoffPositions[sn][0].Multiplicity;
        var ops = TSubgroupFinder.GetExpandedOps(sn);

        // 相異なる線形部 R (= 点群)。
        var distinct = new List<int[,]>();
        foreach (var op in ops)
        {
            var R = SeitzNotation.LinearMatrix(op);
            if (!distinct.Any(d => Same(d, R)))
                distinct.Add(R);
        }
        PointGroupOrder = distinct.Count;

        IsCentrosymmetric = distinct.Any(IsMinusIdentity);
        IsSohncke = distinct.All(d => Det(d) == 1); // 真回転のみ (det=+1)

        // --- 極性方向: 全操作 R に対し R·v = v となる v の共通部分空間 (整数計算) ---
        var (dim, n1, n2) = FixedSubspace(distinct);
        IsPolar = dim >= 1;
        PolarDirectionStr = dim switch
        {
            >= 3 => "any",
            2 => $"any in ({Miller(Cross(n1, n2))})",  // 鏡映面など 2 次元 → 面内自由
            1 => $"[{Uvw(n1)}]",
            _ => "none",
        };

        // --- 物性の対称性許容 ---
        // 260705Cl 修正: StrArray の点群表記には mm2 の設定バリアント "2mm"/"m2m" が存在し (計 67 設定)、
        // 生文字列比較ではそれらの設定で旋光性が誤って「禁止」になっていた。TSubgroupFinder と同じ正規化を適用。
        // is432 も導出ロジックをやめ既存データ (PointGroupHMStr) の直参照に単純化。
        //bool is432 = sym.CrystalSystemNumber == 7 && IsSohncke && PointGroupOrder == 24; // 立方の回転群 O
        //var pg = sym.PointGroupHMStr;
        var pg = sym.PointGroupHMStr switch { "2mm" or "m2m" => "mm2", var t => t };
        bool is432 = pg == "432";
        PyroelectricAllowed = IsPolar;
        PiezoelectricAllowed = !IsCentrosymmetric && !is432;
        SHGAllowed = PiezoelectricAllowed;
        OpticalActivityAllowed = IsSohncke || pg is "m" or "mm2" or "-4" or "-42m";

        // --- 小テーブル / 合成 ---
        int itno = sym.SpaceGroupNumber;
        IsSymmorphic = SymmorphicNumbers.Contains(itno);
        EnantiomorphPartnerNumber = EnantiomorphPairs.TryGetValue(itno, out var partner) ? partner : 0;
        HasEnantiomorph = EnantiomorphPartnerNumber != 0;

        string lat = sym.LatticeTypeStr;
        (CrystalFamilyStr, LatticeSystemStr, BravaisTypeStr) = LatticeDescriptors(sym.CrystalSystemNumber, lat);
        ArithmeticCrystalClassStr = pg + lat;
        PattersonSymmetryStr = Patterson(sym, lat);
    }
    #endregion

    #region 整数線形代数
    private static bool Same(int[,] a, int[,] b)
    {
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                if (a[i, j] != b[i, j]) return false;
        return true;
    }

    private static bool IsMinusIdentity(int[,] m)
    {
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                if (m[i, j] != (i == j ? -1 : 0)) return false;
        return true;
    }

    private static int Det(int[,] m) =>
        m[0, 0] * (m[1, 1] * m[2, 2] - m[1, 2] * m[2, 1])
      - m[0, 1] * (m[1, 0] * m[2, 2] - m[1, 2] * m[2, 0])
      + m[0, 2] * (m[1, 0] * m[2, 1] - m[1, 1] * m[2, 0]);

    /// <summary>全 R について R·v = v を満たす v の部分空間 (∩ ker(R−I)) の次元と基底を返す。
    /// 各 (R−I) を行として積み上げ、rank r を求め dim = 3 − r。dim=1 は null ベクトル、dim=2 は 2 基底。</summary>
    private static (int Dim, double[] N1, double[] N2) FixedSubspace(List<int[,]> mats)
    {
        var rows = new List<double[]>();
        foreach (var m in mats)
            for (int i = 0; i < 3; i++)
                rows.Add([m[i, 0] - (i == 0 ? 1 : 0), m[i, 1] - (i == 1 ? 1 : 0), m[i, 2] - (i == 2 ? 1 : 0)]);
        return NullSpace(rows);
    }

    /// <summary>m×3 行列 (double) の零空間を求める。返り値は (次元, 基底1, 基底2)。次元 0..3。</summary>
    private static (int Dim, double[] N1, double[] N2) NullSpace(List<double[]> rows)
    {
        const double eps = 1e-9;
        var a = rows.Select(r => (double[])r.Clone()).ToList();
        int nrow = a.Count, ncol = 3;
        var pivotCol = new List<int>();
        int r = 0;
        for (int c = 0; c < ncol && r < nrow; c++)
        {
            int sel = -1; double best = eps;
            for (int i = r; i < nrow; i++)
                if (Math.Abs(a[i][c]) > best) { best = Math.Abs(a[i][c]); sel = i; }
            if (sel < 0) continue;
            (a[r], a[sel]) = (a[sel], a[r]);
            double pv = a[r][c];
            for (int j = 0; j < ncol; j++) a[r][j] /= pv;
            for (int i = 0; i < nrow; i++)
                if (i != r && Math.Abs(a[i][c]) > eps)
                {
                    double f = a[i][c];
                    for (int j = 0; j < ncol; j++) a[i][j] -= f * a[r][j];
                }
            pivotCol.Add(c);
            r++;
        }
        var freeCols = Enumerable.Range(0, ncol).Where(c => !pivotCol.Contains(c)).ToList();
        var basis = new List<double[]>();
        foreach (var fc in freeCols)
        {
            var v = new double[ncol];
            v[fc] = 1;
            for (int i = 0; i < pivotCol.Count; i++)
                v[pivotCol[i]] = -a[i][fc]; // RREF から従属成分を復元
            basis.Add(v);
        }
        return (basis.Count,
                basis.Count > 0 ? basis[0] : null,
                basis.Count > 1 ? basis[1] : null);
    }

    private static double[] Cross(double[] a, double[] b) =>
        [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];

    /// <summary>実数方向を最小整数比の "uvw" 文字列へ (負は連結: 1-10)。</summary>
    private static string Uvw(double[] v)
    {
        var (u, w, t) = Integerize(v);
        return $"{u}{w}{t}";
    }

    private static string Miller(double[] v)
    {
        var (h, k, l) = Integerize(v);
        return $"{h}{k}{l}";
    }

    /// <summary>実数 3 ベクトルを最小整数比に整形 (最小非零成分で正規化して丸め)。符号は最初の非零を正に。(260705Cl: 実装と食い違っていた説明を修正)</summary>
    private static (int, int, int) Integerize(double[] v)
    {
        const double eps = 1e-6;
        // 最初の非零を正にする
        double sign = 0;
        foreach (var x in v) if (Math.Abs(x) > eps) { sign = Math.Sign(x); break; }
        if (sign == 0) sign = 1;
        double[] w = [v[0] / sign, v[1] / sign, v[2] / sign];
        // 最小の非零絶対値で正規化 → 整数化
        double min = double.MaxValue;
        foreach (var x in w) if (Math.Abs(x) > eps && Math.Abs(x) < min) min = Math.Abs(x);
        if (min == double.MaxValue) return (0, 0, 0);
        int[] o = new int[3];
        for (int i = 0; i < 3; i++) o[i] = (int)Math.Round(w[i] / min);
        int g = (int)GammaFunction.Gcd(GammaFunction.Gcd(Math.Abs(o[0]), Math.Abs(o[1])), Math.Abs(o[2])); // 260705Cl: 既存 GammaFunction.Gcd に一本化
        if (g > 1) for (int i = 0; i < 3; i++) o[i] /= g;
        return (o[0], o[1], o[2]);
    }

    // 260705Cl: 私製 Gcd を削除し既存 GammaFunction.Gcd (Mathematics/GammaFunction.cs) に一本化。
    //private static int Gcd(int a, int b) { while (b != 0) (a, b) = (b, a % b); return a == 0 ? 1 : a; }
    #endregion

    #region 格子記述子 / Patterson
    private static (string Family, string LatticeSystem, string Bravais) LatticeDescriptors(int crystalSystem, string lat)
    {
        // crystalSystem: 1 triclinic .. 7 cubic。trigonal は格子で hexagonal/rhombohedral が分かれる。
        return crystalSystem switch
        {
            1 => ("triclinic", "triclinic", "a" + lat),
            2 => ("monoclinic", "monoclinic", "m" + lat),
            3 => ("orthorhombic", "orthorhombic", "o" + lat),
            4 => ("tetragonal", "tetragonal", "t" + lat),
            5 => lat == "R"
                    ? ("hexagonal", "rhombohedral", "hR")
                    : ("hexagonal", "hexagonal", "h" + lat),
            6 => ("hexagonal", "hexagonal", "h" + lat),
            7 => ("cubic", "cubic", "c" + lat),
            _ => ("unknown", "unknown", lat),
        };
    }

    /// <summary>Patterson 対称性 = 格子 + Laue 類 HM (中心対称・symmorphic)。三方 -3m は空間群 HM の向きで -31m/-3m1 を判別。</summary>
    private static string Patterson(in Symmetry sym, string lat)
    {
        string laue = sym.LaueGroupStr; // 例 "-1","2/m","mmm","4/m","4/mmm","-3","-32m","6/m","6/mmm","m-3","m-3m"
        // 三方 -3m 系 (データ表記 "-32m") は向きで分岐。
        if (laue == "-32m")
        {
            if (lat == "R") return "R-3m";
            // 260705Cl 修正: 旧 HM 文字列パース (TrigonalMinus3mOrientation) は螺旋添字 (1,2) をスキップできず、
            // P3sub121/P3sub212 (No.152/153, 水晶) で -31m/-3m1 が逆転していた。構造化データの第 2 方向
            // 回転要素 (StrSE2p: 312/31m/31c 系="1"、321/3m1/3c1 系="2" または空) で判別する。
            //laue = TrigonalMinus3mOrientation(sym.SpaceGroupHMStr);
            laue = sym.StrSE2p == "1" ? "-31m" : "-3m1";
        }
        else if (laue == "-3" && lat == "R")
            return "R-3";
        return lat + laue;
    }

    // 260705Cl: 上記の修正に伴い廃止 (螺旋添字を含む HM で誤判別するバグがあった)。
    ///// <summary>P-31m 系か P-3m1 系か: HM の '3' 直後の文字が '1' なら -31m、それ以外なら -3m1。</summary>
    //private static string TrigonalMinus3mOrientation(string hm)
    //{
    //    var s = hm.Replace("sub", "").Replace("-", "").Replace("Hex", "").Replace("Rho", "");
    //    int idx = s.IndexOf('3');
    //    if (idx >= 0 && idx + 1 < s.Length)
    //    {
    //        // '3' の次の桁 (螺旋番号) を飛ばして最初の英字/桁を見る
    //        int j = idx + 1;
    //        while (j < s.Length && char.IsDigit(s[j]) && s[j] != '1' && s[j] != '2') j++;
    //        if (j < s.Length && s[j] == '1') return "-31m";
    //    }
    //    return "-3m1";
    //}
    #endregion

    #region 静的テーブル
    /// <summary>Symmorphic 空間群の IT 番号 (73 群)。</summary>
    private static readonly HashSet<int> SymmorphicNumbers =
    [
        1, 2, 3, 5, 6, 8, 10, 12, 16, 21, 22, 23, 25, 35, 38, 42, 44, 47, 65, 69, 71,
        75, 79, 81, 82, 83, 87, 89, 97, 99, 107, 111, 115, 119, 121, 123, 139,
        143, 146, 147, 148, 149, 150, 155, 156, 157, 160, 162, 164, 166,
        168, 174, 175, 177, 183, 187, 189, 191,
        195, 196, 197, 200, 202, 204, 207, 209, 211, 215, 216, 217, 221, 225, 229,
    ];

    /// <summary>エナンチオモルフ (掌性) 対 11 組 (双方向)。IT 番号 → 相手の IT 番号。</summary>
    private static readonly Dictionary<int, int> EnantiomorphPairs = new()
    {
        [76] = 78, [78] = 76,     // P4₁ / P4₃
        [91] = 95, [95] = 91,     // P4₁22 / P4₃22
        [92] = 96, [96] = 92,     // P4₁2₁2 / P4₃2₁2
        [144] = 145, [145] = 144, // P3₁ / P3₂
        [151] = 153, [153] = 151, // P3₁12 / P3₂12
        [152] = 154, [154] = 152, // P3₁21 / P3₂21
        [169] = 170, [170] = 169, // P6₁ / P6₅
        [171] = 172, [172] = 171, // P6₂ / P6₄
        [178] = 179, [179] = 178, // P6₁22 / P6₅22
        [180] = 181, [181] = 180, // P6₂22 / P6₄22
        [212] = 213, [213] = 212, // P4₃32 / P4₁32
    };
    #endregion
}
