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

    /// <summary>結晶族 (英語の小文字。triclinic … cubic)。画面の訳は表示側 (FormSymmetryInformation) で付ける (261105Cl)。</summary>
    public string CrystalFamilyStr { get; }
    /// <summary>格子系 (英語の小文字。triclinic … rhombohedral, hexagonal, cubic)。画面の訳は表示側で付ける (261105Cl)。</summary>
    public string LatticeSystemStr { get; }
    /// <summary>Bravais 型。ITA の正式な 14 型で、設定の格子記号が違うときは括弧で併記 (例 "oS (oA)")。261105Cl 変更。</summary>
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
        var pg = sym.PointGroupHMStr switch { "2mm" or "m2m" => "mm2", var t => t };
        bool is432 = pg == "432";
        PyroelectricAllowed = IsPolar;
        PiezoelectricAllowed = !IsCentrosymmetric && !is432;
        SHGAllowed = PiezoelectricAllowed;
        OpticalActivityAllowed = IsSohncke || pg is "m" or "mm2" or "-4" or "-42m";

        // --- 小テーブル / 合成 ---
        int itno = sym.SpaceGroupNumber;
        IsSymmorphic = SymmorphicNumbers.Contains(itno);
        EnantiomorphPartnerNumber = GetEnantiomorphPartnerNumber(itno);
        HasEnantiomorph = EnantiomorphPartnerNumber != 0;

        string lat = sym.LatticeTypeStr;
        (CrystalFamilyStr, LatticeSystemStr, BravaisTypeStr) = LatticeDescriptors(sym.CrystalSystemNumber, lat);
        //ArithmeticCrystalClassStr = pg + lat; // 261105Cl 変更前: 点群の記号に向きが無いので 312P/321P が両方「32P」など 6 組が同じ記号になり、
        //                                       非標準設定では中心化の記号が設定のもの (A2/m の「2/mA」、I2/a の「2/mI」) になって 73 類より多くなった (検査で 91 種)
        ArithmeticCrystalClassStr = ArithmeticCrystalClass(itno) ?? pg + lat; // 261105Cl 変更: IT 番号から決める (設定に依らない)
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

    //260725Cl 変更: private → internal (EbsdDictionaryIndexer.GetProperRotations と共有 — int 3×3 行列式の重複実装解消)
    internal static int Det(int[,] m) =>
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
    #endregion

    #region 格子記述子 / Patterson
    /// <summary>261105Cl 追加: IT 番号 (1..230) から算術結晶類 (73 類) の記号を返す。範囲外は null。
    /// 記号は標準設定の点群 + 格子 (ITA の表記)。向きで分かれる 6 組 (312P/321P, 3m1P/31mP, -31mP/-3m1P, -6m2P/-62mP, -42mP/-4m2P, -42mI/-4m2I)
    /// を区別し、非標準設定 (A2/m, I2/a, Bmm2 など) も標準設定の記号にする。Amm2 型 (38–41) は mm2A (ITA の 2mmC と同じ類)。</summary>
    private static string ArithmeticCrystalClass(int itno) => itno switch
    {
        1 => "1P", 2 => "-1P",
        3 or 4 => "2P", 5 => "2C", 6 or 7 => "mP", 8 or 9 => "mC", 10 or 11 or 13 or 14 => "2/mP", 12 or 15 => "2/mC",
        >= 16 and <= 19 => "222P", 20 or 21 => "222C", 22 => "222F", 23 or 24 => "222I",
        >= 25 and <= 34 => "mm2P", >= 35 and <= 37 => "mm2C", >= 38 and <= 41 => "mm2A", 42 or 43 => "mm2F", >= 44 and <= 46 => "mm2I",
        >= 47 and <= 62 => "mmmP", >= 63 and <= 68 => "mmmC", 69 or 70 => "mmmF", >= 71 and <= 74 => "mmmI",
        >= 75 and <= 78 => "4P", 79 or 80 => "4I", 81 => "-4P", 82 => "-4I", >= 83 and <= 86 => "4/mP", 87 or 88 => "4/mI",
        >= 89 and <= 96 => "422P", 97 or 98 => "422I", >= 99 and <= 106 => "4mmP", >= 107 and <= 110 => "4mmI",
        >= 111 and <= 114 => "-42mP", >= 115 and <= 118 => "-4m2P", 119 or 120 => "-4m2I", 121 or 122 => "-42mI",
        >= 123 and <= 138 => "4/mmmP", >= 139 and <= 142 => "4/mmmI",
        >= 143 and <= 145 => "3P", 146 => "3R", 147 => "-3P", 148 => "-3R",
        149 or 151 or 153 => "312P", 150 or 152 or 154 => "321P", 155 => "32R",
        156 or 158 => "3m1P", 157 or 159 => "31mP", 160 or 161 => "3mR",
        162 or 163 => "-31mP", 164 or 165 => "-3m1P", 166 or 167 => "-3mR",
        >= 168 and <= 173 => "6P", 174 => "-6P", 175 or 176 => "6/mP", >= 177 and <= 182 => "622P", >= 183 and <= 186 => "6mmP",
        187 or 188 => "-6m2P", 189 or 190 => "-62mP", >= 191 and <= 194 => "6/mmmP",
        195 or 198 => "23P", 196 => "23F", 197 or 199 => "23I",
        200 or 201 or 205 => "m-3P", 202 or 203 => "m-3F", 204 or 206 => "m-3I",
        207 or 208 or 212 or 213 => "432P", 209 or 210 => "432F", 211 or 214 => "432I",
        215 or 218 => "-43mP", 216 or 219 => "-43mF", 217 or 220 => "-43mI",
        >= 221 and <= 224 => "m-3mP", >= 225 and <= 228 => "m-3mF", 229 or 230 => "m-3mI",
        _ => null,
    };

    private static (string Family, string LatticeSystem, string Bravais) LatticeDescriptors(int crystalSystem, string lat)
    {
        // crystalSystem: 1 triclinic .. 7 cubic。trigonal は格子で hexagonal/rhombohedral が分かれる。
        // 261105Cl 変更: Bravais 型は ITA の正式な 14 型 (aP, mP, mS, oP, oS, oI, oF, tP, tI, hP, hR, cP, cI, cF) で示し、設定の格子記号が
        //   正式な型と違うときは括弧で併記する (例 oS (oA)、mS (mI)、aP (aC))。旧は「結晶族の文字 + 設定の格子記号」で、oA・mI・aC など
        //   設定の記号が正式な型のように出ていた。中心化の対応: 三斜晶は全て aP、単斜晶の A・B・C・I は mS、直方晶の A・B・C は oS、
        //   正方晶の C は tP・F は tI (今の 530 余の設定には無いが ITA に沿う)、六方晶の H 胞は hP、R は六方軸・菱面体軸とも hR。
        string Bravais(string formal) => formal[1..] == lat ? formal : $"{formal} ({formal[0]}{lat})"; // 261105Cl 追加
        return crystalSystem switch
        {
            //1 => ("triclinic", "triclinic", "a" + lat), // 261105Cl 変更前 (以下 4 行も)
            //2 => ("monoclinic", "monoclinic", "m" + lat),
            //3 => ("orthorhombic", "orthorhombic", "o" + lat),
            //4 => ("tetragonal", "tetragonal", "t" + lat),
            1 => ("triclinic", "triclinic", Bravais("aP")), // 261105Cl 変更
            2 => ("monoclinic", "monoclinic", Bravais(lat == "P" ? "mP" : "mS")), // 261105Cl 変更
            3 => ("orthorhombic", "orthorhombic", Bravais(lat is "A" or "B" or "C" ? "oS" : "o" + lat)), // 261105Cl 変更
            4 => ("tetragonal", "tetragonal", Bravais(lat is "P" or "C" ? "tP" : "tI")), // 261105Cl 変更
            5 => lat == "R"
                    ? ("hexagonal", "rhombohedral", "hR")
                    //: ("hexagonal", "hexagonal", "h" + lat), // 261105Cl 変更前
                    : ("hexagonal", "hexagonal", Bravais("hP")), // 261105Cl 変更
            //6 => ("hexagonal", "hexagonal", "h" + lat), // 261105Cl 変更前
            6 => ("hexagonal", "hexagonal", Bravais("hP")), // 261105Cl 変更
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

    /// <summary>掌性対の相手の IT 番号を返す (無ければ 0)。260708Cl 追加: KSubgroupFinder の同型 (IIc) 判定と
    /// テーブルを共有するための静的入口 (二重保守防止)。</summary>
    public static int GetEnantiomorphPartnerNumber(int itNumber) => EnantiomorphPairs.GetValueOrDefault(itNumber);
    #endregion
}
