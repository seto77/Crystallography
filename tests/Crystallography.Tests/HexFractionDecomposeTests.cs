using Xunit;

namespace Crystallography.Tests;

/// <summary>261125Cl 追加 (codex2 の 11-25 の議論 観点 I §1-2): 六方晶の座標の近似小数を分数に寄せるとき、符号・整数部・指数を保つこと
/// (旧は "-0.3333" を +1/3、"1.3333" を 1/3、"0.333E-1" を 1/3 にしていた)</summary>
public class HexFractionDecomposeTests
{
    [Theory] // 261125Cl 追加
    [InlineData("0.3333", 1.0 / 3.0)]
    [InlineData("-0.3333", -1.0 / 3.0)]
    [InlineData("0.6667", 2.0 / 3.0)]
    [InlineData("-0.6667", -2.0 / 3.0)]
    [InlineData("1.3333", 4.0 / 3.0)]
    [InlineData("-0.0833", -1.0 / 12.0)]
    [InlineData("0.3333|2|", 1.0 / 3.0)]
    [InlineData("0.333E-1", 0.0333)]
    [InlineData("0.25", 0.25)]
    public void HexValuesKeepSignIntegerPartAndExponent(string text, double expected)
    {
        var (value, _) = Crystal2.Decompose(text, IsHex: true);
        Assert.Equal(expected, value, 1e-12);
    }

    [Theory] // 261125Cl 追加 (codex2 の査読 L §2・§3): Compose が書く指数の形 ("3.333|20|E-1") でも 1/3 に戻り、誤差にも指数が掛かる
    [InlineData("3.333|20|E-1", 1.0 / 3.0, 0.002)]
    [InlineData("-6.667|3|E-1", -2.0 / 3.0, 0.0003)]
    [InlineData("0.3333|2|", 1.0 / 3.0, 0.0002)]
    public void HexExponentFormRoundTrips(string text, double expectedValue, double expectedError)
    {
        var (value, error) = Crystal2.Decompose(text, IsHex: true);
        Assert.Equal(expectedValue, value, 1e-12);
        Assert.Equal(expectedError, error, 1e-12);
    }

    [Fact] // 261125Cl 追加: 六方晶でなければ小数のまま
    public void NonHexIsUnchanged()
    {
        Assert.Equal(-0.3333, Crystal2.Decompose("-0.3333", IsHex: false).Value, 1e-15);
    }
}
