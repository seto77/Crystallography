using System.Globalization;
using System.Text;
using System.Xml;
using Xunit;

namespace Crystallography.Tests;

/// <summary>261105Cl 新規: DiffractionOptics.ParseNumber (.prm の数値の読み込み) が地域設定に依らず、
/// 新形式 (InvariantCulture, "0.1") と旧形式 (書いた PC の地域設定, "0,1" や U+2212 の負号) の両方を正しく読み、
/// 桁区切りと区別できない文字列は推測せず例外にすることを確かめる。</summary>
public class PrmNumberTests
{
    private static readonly string[] Cultures = ["en-US", "ja-JP", "de-DE", "fr-FR", "sv-SE"];

    public static TheoryData<string, double> Valid => new()
    {
        { "0.1", 0.1 },
        { "0,1", 0.1 },
        { " 12.5 ", 12.5 },
        { "-1.5", -1.5 },
        { "\u22121,5", -1.5 },          // sv-SE などの旧形式 (U+2212 の負号)
        { "\u22121,5E\u221205", -1.5e-5 },
        { "1.234E-05", 1.234e-5 },
        { "1500", 1500 },
        { "-1", -1 },
        { "NaN", double.NaN },
    };

    public static TheoryData<string> Invalid => new() { "1,234.5", "1.234,5", "1,234,567", "1.234.5", "abc", "" };

    [Theory]
    [MemberData(nameof(Valid))]
    public void ParseNumber_ReadsBothFormats_InAnyCulture(string s, double expected)
    {
        foreach (var c in Cultures)
            WithCulture(c, () => Assert.Equal(expected, DiffractionOptics.ParseNumber(s)));
    }

    [Theory]
    [MemberData(nameof(Invalid))]
    public void ParseNumber_RejectsAmbiguousOrInvalidText(string s)
    {
        foreach (var c in Cultures)
            WithCulture(c, () => Assert.Throws<FormatException>(() => DiffractionOptics.ParseNumber(s)));
    }

    [Fact]
    public void ParseNumber_Null_Throws() => Assert.Throws<ArgumentNullException>(() => DiffractionOptics.ParseNumber(null));

    /// <summary>旧形式 (書いた PC の CurrentCulture) と新形式 (InvariantCulture) の .prm を、5 つの地域設定で書いて 5 つの地域設定で読む。
    /// 書き方は IPAnalyzer の SaveParameter と同じ (XmlSerializer + StreamWriter(UTF8) + XmlWriter Indent)。</summary>
    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void PrmFile_RoundTrips_AcrossCultures(bool invariantWriter)
    {
        (string field, double value)[] values =
        [
            ("DirectSpotX", 1024.5), ("DirectSpotY", -3.25e-5), ("CameraLength1", 123.456789), ("CameraLength2", 98.7654321),
            ("FootX", 1500.25), ("FootY", -0.125), ("waveLength", 0.4132), ("pixSizeX", 0.1), ("pixSizeY", 0.100000001),
            ("pixKsi", -0.0123), ("tiltPhi", -1.5), ("tiltTau", 179.99), ("ConcentricStep", 0.01), ("RectangleDirection", -1),
        ];
        var t = typeof(DiffractionOptics.Parameter);
        foreach (var writer in Cultures)
        {
            var fn = Path.Combine(Path.GetTempPath(), $"PrmNumberTests_{Guid.NewGuid():N}.prm");
            try
            {
                WithCulture(writer, () =>
                {
                    var prm = new DiffractionOptics.Parameter();
                    foreach (var (f, v) in values)
                        t.GetField(f).SetValue(prm, invariantWriter ? v.ToString(CultureInfo.InvariantCulture) : v.ToString());
                    using var sw = new StreamWriter(fn, false, new UTF8Encoding());
                    using var xw = XmlWriter.Create(sw, new XmlWriterSettings { Indent = true });
                    new System.Xml.Serialization.XmlSerializer(typeof(DiffractionOptics.Parameter)).Serialize(xw, prm);
                });
                foreach (var reader in Cultures)
                    WithCulture(reader, () =>
                    {
                        var prm = DiffractionOptics.Read(fn);
                        Assert.NotNull(prm);
                        foreach (var (f, v) in values)
                            Assert.True(v.Equals(DiffractionOptics.ParseNumber((string)t.GetField(f).GetValue(prm))),
                                $"{f}: written in {writer} ({(invariantWriter ? "invariant" : "current culture")}), read in {reader}");
                    });
            }
            finally { File.Delete(fn); }
        }
    }

    private static void WithCulture(string name, Action action)
    {
        var saved = CultureInfo.CurrentCulture;
        CultureInfo.CurrentCulture = new CultureInfo(name);
        try { action(); }
        finally { CultureInfo.CurrentCulture = saved; }
    }
}
