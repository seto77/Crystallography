using System.Reflection;
using Xunit;
using static Crystallography.Tests.SymmetryReference;

namespace Crystallography.Tests;

/// <summary>261105Cl 新規: CIF 読み込みで空間群の設定を決める ConvertCrystalData.SearchSGseriesNumberForCIF (private) を、
/// 530 設定それぞれの「よくある CIF の書き方」で呼び、正しい設定 (対称操作が一致する通し番号) が返るかを確かめる。
/// 呼び出し時の引数 (セルの形の判定) は ConvertFromCIF と同じ規則で作る。</summary>
public class CifSpaceGroupSearchTests
{
    private static readonly MethodInfo Search = typeof(ConvertCrystalData).GetMethod("SearchSGseriesNumberForCIF", BindingFlags.NonPublic | BindingFlags.Static)
        ?? throw new MissingMethodException(nameof(ConvertCrystalData), "SearchSGseriesNumberForCIF");

    /// <summary>CIF の H-M 記号 (例 'P 1 21/n 1', 'F 41/d -3 2/m :2', 'R -3 2/m :H')。spglib の完全記号から '_' を除き、
    /// 原点選択 (1/2, 軸の取り方が違う '2ba-c' などを含む) と六方晶軸/菱面体軸の接尾辞を付ける。</summary>
    private static string CifHm(Group r)
        => r.HmFull.Replace("_", "") + r.Choice switch
        {
            "H" => " :H",
            "R" => " :R",
            var c when c.StartsWith('1') => " :1",
            var c when c.StartsWith('2') => " :2",
            _ => "",
        };

    /// <summary>同じ CIF H-M 記号になる設定 (自分を含む)。例: Cmme (abc) と Cmme (ba-c) は完全記号が同じで、H-M 記号だけでは区別できない</summary>
    private static int[] SameHmSettings(int series)
        => [.. Enumerable.Range(1, 530).Where(i => CifHm(Groups[i]) == CifHm(Groups[series]))];

    /// <summary>ConvertFromCIF が渡すセルの形の判定 (a=b=c かつ α=β=γ, α=90°, β=90°, γ=90°) を設定から作る</summary>
    private static (bool rho, bool a90, bool b90, bool g90) CellFlags(int series)
    {
        var sym = SymmetryStatic.Symmetries[series];
        var choice = Groups[series].Choice;
        return sym.CrystalSystemNumber switch
        {
            1 => (false, false, false, false),
            2 => (false, sym.MainAxis != "a", sym.MainAxis != "b", sym.MainAxis != "c"),
            5 or 6 when choice == "R" => (true, false, false, false),
            5 or 6 => (false, true, true, false),
            7 => (true, true, true, true),
            _ => (false, true, true, true),
        };
    }

    //private static int Invoke(string hm, string hall, int number, int series) // 261105Cl 変更前
    private static int Invoke(string hm, string hall, int number, int series, string[] ops = null) // 261105Cl ops (CIF の対称操作) を追加
    {
        var (rho, a90, b90, g90) = CellFlags(series);
        //return (int)Search.Invoke(null, [hm, hall, number, rho, a90, b90, g90]); // 261105Cl 変更前
        return (int)Search.Invoke(null, [hm, hall, number, rho, a90, b90, g90, ops]); // 261105Cl: 引数 operations (既定値 null) もリフレクションでは明示する
    }

    private static void AssertSameGroup(int expected, int got, string input)
    {
        Assert.True(got > 0, $"{input} → {got} (見つからない)");
        Assert.True(Groups[expected].OpKeys.SetEquals(GeneralPositionKeys(got)),
            $"{input} → {got} ({SymmetryStatic.Symmetries[got].SpaceGroupHMStr}), 正しくは {expected} ({SymmetryStatic.Symmetries[expected].SpaceGroupHMStr})");
    }

    /// <summary>標準の Hall 記号 (spglib / ITA の表記, 空白あり) と H-M 記号・番号を持つ CIF</summary>
    [Theory, MemberData(nameof(ItaSeries), MemberType = typeof(SymmetryReference))]
    public void StandardHallSymbol_IsAssignedToCorrectSetting(int series)
    {
        var r = Groups[series];
        AssertSameGroup(series, Invoke(CifHm(r), r.HallSymbol, r.Number, series), $"Hall '{r.HallSymbol}', H-M '{CifHm(r)}'");
    }

    /// <summary>旧版の ReciPro が書き出した CIF (Hall 記号が空白なしの従来表記 SpaceGroupHallStr) も正しい設定に戻る。
    /// 'P32' (P3_2 と P312) や 'P62' (P6_2 と P622) のように 1 つに決まらないものは H-M 記号で決まる。</summary>
    [Theory, MemberData(nameof(ItaSeries), MemberType = typeof(SymmetryReference))]
    public void LegacyReciProHallSymbol_IsAssignedToCorrectSetting(int series)
    {
        var r = Groups[series];
        var legacy = SymmetryStatic.Symmetries[series].SpaceGroupHallStr;
        AssertSameGroup(series, Invoke(CifHm(r), legacy, r.Number, series), $"旧 Hall '{legacy}', H-M '{CifHm(r)}'");
    }

    /// <summary>Hall 記号が無く、H-M 記号 (完全記号) と番号だけを持つ CIF。
    /// H-M 記号だけでは区別できない設定の組 (SameHmSettings が複数の群を含む) では、その組のどれかに割り当てられれば良しとする。</summary>
    [Theory, MemberData(nameof(ItaSeries), MemberType = typeof(SymmetryReference))]
    public void HermannMauguinOnly_IsAssignedToCorrectSetting(int series)
    {
        var r = Groups[series];
        var got = Invoke(CifHm(r), "", r.Number, series);
        var input = $"H-M '{CifHm(r)}', No. {r.Number}";
        var candidates = SameHmSettings(series);
        if (candidates.Select(i => string.Join("|", Groups[i].OpKeys.Order())).Distinct().Count() == 1)
            AssertSameGroup(series, got, input);
        else
            Assert.True(got > 0 && candidates.Any(i => Groups[i].OpKeys.SetEquals(GeneralPositionKeys(got))),
                $"{input} → {got}, 候補 {string.Join(",", candidates)} のどれでもない");
    }

    #region 261105Cl 追加: H-M 記号が無く、空白なしの Hall 記号が複数の設定に当たる CIF
    // 'P32' (P3_2 と P312)、'P62' (P6_2 と P622)、'C22-1bc' (Ccce(1) の abc と ba-c) など、標準表記・従来表記で計 10 通り。
    // 以前は H-M 記号の経路に空文字列のまま進み、SgNameHM[0] で IndexOutOfRangeException になって CIF が開けなかった。

    /// <summary>261105Cl 追加: 空白を除くと複数の設定に当たる Hall 記号と、その設定 (標準表記・従来表記の両方。'P32' と 'P62' は両方で同じ文字列なので 1 回だけ)。
    /// 8 通りの記号 × 2 設定 = 16 組</summary>
    public static IEnumerable<object[]> AmbiguousSpacelessHall()
        => (from standard in new[] { true, false }
            from g in Enumerable.Range(1, 530).GroupBy(i => (standard ? SymmetryStatic.Symmetries[i].SpaceGroupHallStandardStr : SymmetryStatic.Symmetries[i].SpaceGroupHallStr).Replace(" ", ""))
            where g.Count() > 1
            from series in g
            select (g.Key, series)).Distinct().Select(t => new object[] { t.Key, t.series });

    /// <summary>261105Cl 追加: 同じ空白なし Hall 記号になる設定がすべて同じ群か</summary>
    private static bool CandidatesAreSameGroup(string hall)
        => AmbiguousSpacelessHall().Where(o => (string)o[0] == hall).Select(o => string.Join("|", Groups[(int)o[1]].OpKeys.Order())).Distinct().Count() == 1;

    /// <summary>261105Cl 追加: 基準の対称操作 [r11..r33, t1..t3 (1/24 単位)] を CIF の書き方 ('-x+y, -x, z+2/3') にする</summary>
    private static string Xyz(int[] op)
    {
        static int gcd(int a, int b) => b == 0 ? a : gcd(b, a % b);
        string row(int k)
        {
            var s = "";
            for (int j = 0; j < 3; j++)
                if (op[3 * k + j] != 0)
                    s += (op[3 * k + j] switch { 1 => "+", -1 => "-", var c => c > 0 ? $"+{c}" : $"{c}" }) + "xyz"[j];
            var t = ((op[9 + k] % Den) + Den) % Den;
            if (t != 0)
                s += $"+{t / gcd(t, Den)}/{Den / gcd(t, Den)}";
            return s.TrimStart('+');
        }
        return $"{row(0)}, {row(1)}, {row(2)}";
    }

    private static string[] CifOps(int series) => [.. Groups[series].Ops.Select(Xyz)];

    /// <summary>261105Cl 追加: Hall 記号だけ (H-M 記号・番号・対称操作なし)。候補が同じ群なら採り、違う群なら決まらない (-1) と返す。例外にならない</summary>
    [Theory, MemberData(nameof(AmbiguousSpacelessHall))]
    public void AmbiguousSpacelessHallOnly_DoesNotThrow(string hall, int series)
    {
        var got = Invoke("", hall, -1, series);
        if (CandidatesAreSameGroup(hall))
            AssertSameGroup(series, got, $"Hall '{hall}' のみ");
        else
            Assert.Equal(-1, got);
    }

    /// <summary>261105Cl 追加: Hall 記号と CIF の対称操作だけ (H-M 記号・番号なし)。対称操作の一致する設定が選ばれる</summary>
    [Theory, MemberData(nameof(AmbiguousSpacelessHall))]
    public void AmbiguousSpacelessHall_WithOperations_IsAssignedToCorrectSetting(string hall, int series)
        => AssertSameGroup(series, Invoke("", hall, -1, series, CifOps(series)), $"Hall '{hall}' + 対称操作 ({Groups[series].HmShort})");

    /// <summary>261105Cl 追加: Hall 記号と番号だけ (H-M 記号・対称操作なし)。番号で決まる</summary>
    [Theory, MemberData(nameof(AmbiguousSpacelessHall))]
    public void AmbiguousSpacelessHall_WithNumber_IsAssignedToCorrectSetting(string hall, int series)
        => AssertSameGroup(series, Invoke("", hall, Groups[series].Number, series), $"Hall '{hall}' + No. {Groups[series].Number}");

    /// <summary>261105Cl 追加: 対称操作が読めない書き方 (小数) なら使わず、番号で決まる (読めない式で止まったり例外になったりしない)</summary>
    [Fact]
    public void AmbiguousSpacelessHall_UnreadableOperations_FallBackToNumber()
    {
        string[] ops = ["x,y,z", "-y,x-y,z+0.6666667", "-x+y,-x,z+0.3333333"];
        AssertSameGroup(432, Invoke("", "P32", 145, 432, ops), "Hall 'P32' + 小数の対称操作 + No. 145");
        Assert.Equal(-1, Invoke("", "P32", -1, 432, ops));
        Assert.Equal(-1, Invoke("", "P32", -1, 432, ["X,Y,Z", "foo", "x,y"]));
    }

    /// <summary>261105Cl 追加: CIF 全体 (ConvertFromCifText) を読む。以前は例外を内部で捕まえて null (開けない) になっていた</summary>
    [Theory]
    [InlineData("P32", 438, true, -1)]   // P312 の対称操作 → P312
    [InlineData("P32", 432, true, -1)]   // P3_2 の対称操作 → P3_2
    [InlineData("P32", 432, false, 145)] // 番号 145 → P3_2
    [InlineData("P62", 471, true, -1)]   // P622 の対称操作 → P622
    [InlineData("C22-1bc", 322, false, -1)] // Ccce(1): abc と ba-c は同じ群
    [InlineData("A22-1ac", 326, false, -1)]
    [InlineData("B22-1bc", 330, false, -1)]
    public void ConvertFromCifText_AmbiguousSpacelessHallWithoutHm_Opens(string hall, int series, bool withOps, int number)
    {
        var hex = SymmetryStatic.Symmetries[series].CrystalSystemNumber is 5 or 6;
        var cif = new System.Text.StringBuilder();
        cif.Append("data_test\n_chemical_name_common test\n");
        cif.Append($"_cell_length_a 5.0\n_cell_length_b {(hex ? "5.0" : "6.0")}\n_cell_length_c 7.0\n");
        cif.Append($"_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma {(hex ? "120" : "90")}\n");
        cif.Append($"_symmetry_space_group_name_Hall '{hall}'\n");
        if (number > 0)
            cif.Append($"_symmetry_Int_Tables_number {number}\n");
        if (withOps)
        {
            cif.Append("loop_\n_symmetry_equiv_pos_as_xyz\n");
            foreach (var op in CifOps(series))
                cif.Append($"'{op}'\n");
        }
        cif.Append("loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\nSi1 Si 0.1 0.2 0.3\n");

        var c2 = ConvertCrystalData.ConvertFromCifText(cif.ToString());
        Assert.NotNull(c2);
        AssertSameGroup(series, c2.sym, $"CIF (Hall '{hall}', 対称操作 {withOps}, No. {number})");
    }
    #endregion
}
