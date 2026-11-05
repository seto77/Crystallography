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

    private static int Invoke(string hm, string hall, int number, int series)
    {
        var (rho, a90, b90, g90) = CellFlags(series);
        return (int)Search.Invoke(null, [hm, hall, number, rho, a90, b90, g90]);
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
}
