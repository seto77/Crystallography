using Xunit;

namespace Crystallography.Tests;

/// <summary>261105Cl 新規: AffineGen.Parse (対称操作の x,y,z 表記の読み取り)。読める式は正しく、読めない式は FormatException ですぐに止まる
/// (以前は小数点や大文字の X などで先へ進まず無限ループになっていた。無限ループに戻ってもテストが止まらないよう、別スレッドで時間を区切る)。</summary>
public class AffineGenParseTests
{
    private static Exception ParseWithin(string a, string b, string c, out AffineGen result)
    {
        AffineGen r = default;
        var task = Task.Run(() => Record.Exception(() => { r = AffineGen.Parse(a, b, c); }));
        Assert.True(task.Wait(TimeSpan.FromSeconds(10)), $"AffineGen.Parse(\"{a}\", \"{b}\", \"{c}\") が終わらない (無限ループ)");
        result = r;
        return task.Result;
    }

    [Theory]
    [InlineData("x", "y", "z", 0.1, 0.2, 0.3)]
    [InlineData("-x+y", "-x", "z+2/3", 0.1, -0.1, 0.3 + 2.0 / 3)]
    [InlineData("2x", "x+1/2", " -z + 1/4 ", 0.2, 0.6, -0.3 + 0.25)]
    [InlineData("1/2-x", "+y", "-1/4", 0.4, 0.2, -0.25)]
    [InlineData("x-y", "1", "0", -0.1, 1, 0)]
    [InlineData("x - 1/4", "- y", "z", -0.15, -0.2, 0.3)] // 符号の後の空白 (以前は符号が捨てられていた)
    public void ValidExpressions_AreParsed(string a, string b, string c, double ex, double ey, double ez)
    {
        Assert.Null(ParseWithin(a, b, c, out var g));
        var (x, y, z) = g.Apply(0.1, 0.2, 0.3);
        Assert.Equal(ex, x, 12);
        Assert.Equal(ey, y, 12);
        Assert.Equal(ez, z, 12);
    }

    [Theory]
    [InlineData("x+0.5")]      // 小数 (以前は無限ループ)
    [InlineData("X")]          // 大文字 (以前は無限ループ)
    [InlineData("2*x")]        // 掛け算の記号
    [InlineData("x+(1/2)")]    // 括弧
    [InlineData("x+1/")]       // 分母が無い
    [InlineData("x+/2")]       // 分子が無い
    [InlineData("x+1/0")]      // 分母が 0
    [InlineData("x+")]         // 符号の後に項が無い
    [InlineData("x--y")]       // 符号が続く
    [InlineData("")]           // 空
    [InlineData("99999999999x")] // 桁あふれ
    public void InvalidExpressions_ThrowFormatException(string bad)
    {
        var e = ParseWithin("x", "y", bad, out _);
        Assert.IsType<FormatException>(e);
        Assert.Contains("AffineGen.Parse", e.Message);
    }

    /// <summary>空間群の表 (PositionStringList.csv) の全行が読めている (型の初期化で Parse を通る)</summary>
    [Fact]
    public void SymmetryTables_AreParsed()
        => Assert.All(SymmetryStatic.WyckoffPositions[1..], w => Assert.NotEmpty(w[0].PositionGenerator));
}
