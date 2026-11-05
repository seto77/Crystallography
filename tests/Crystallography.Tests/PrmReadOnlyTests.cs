using Xunit;

namespace Crystallography.Tests;

/// <summary>261105Cl 新規: DiffractionOptics.Read が .prm の元ファイルを書き換えないこと (旧実装は旧版の要素名を置き換えた本文を
/// 同じファイルへ書き戻し、時刻だけ元に戻していた)。旧版の要素名 (centerX / cameraLength / SACLA_EH5_FootX) は従来どおり読めること。</summary>
public class PrmReadOnlyTests
{
    [Fact]
    public void ReadDoesNotRewriteTheFileAndStillReadsOldNames()
    {
        const string oldFormat =
            "<?xml version=\"1.0\" encoding=\"utf-8\"?>\r\n" +
            "<parameter xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xmlns:xsd=\"http://www.w3.org/2001/XMLSchema\">\r\n" +
            "  <centerX>1500.25</centerX>\r\n" +
            "  <centerY>-0.125</centerY>\r\n" +
            "  <cameraLength>455.5</cameraLength>\r\n" +
            "  <SACLA_EH5_FootX>12.5</SACLA_EH5_FootX>\r\n" +
            "</parameter>\r\n";
        var fn = Path.Combine(Path.GetTempPath(), $"PrmReadOnlyTests_{Guid.NewGuid():N}.prm");
        try
        {
            File.WriteAllText(fn, oldFormat);
            var before = File.ReadAllBytes(fn);
            var prm = DiffractionOptics.Read(fn);
            Assert.NotNull(prm);
            Assert.Equal(before, File.ReadAllBytes(fn));
            Assert.Equal(1500.25, DiffractionOptics.ParseNumber(prm.DirectSpotX));
            Assert.Equal(-0.125, DiffractionOptics.ParseNumber(prm.DirectSpotY));
            Assert.Equal(455.5, DiffractionOptics.ParseNumber(prm.CameraLength1));
            Assert.Equal(12.5, DiffractionOptics.ParseNumber(prm.FootX));
        }
        finally { File.Delete(fn); }
    }
}
