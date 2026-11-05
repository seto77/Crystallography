using Xunit;

namespace Crystallography.Tests;

/// <summary>261105Cl 新規: XYFile.ReadPdi2File (PDIndexer の .pdi / .pdi2 読み込み) が元ファイルを書き換えないこと、
/// 名前やコメント中の "pt" を変えないこと、旧版の要素名 &lt;pt&gt; (現 Profile.Pt) の互換読み込みは残っていることを確かめる。
/// (旧実装は全行に Replace("pt", "Pt") を掛けて同じファイルへ上書き保存していた)</summary>
public class PdiFileTests
{
    private const int NumPoints = 50;

    private static DiffractionProfile2 Sample()
    {
        var dp = new DiffractionProfile2 { Name = "Sample opt pt", Comment = "script: adapt, Pt foil" };
        for (int i = 0; i < NumPoints; i++)
            dp.SourceProfile.Pt.Add(new PointD(10 + i * 0.02, 100 + i));
        return dp;
    }

    private static string TempFile() => Path.Combine(Path.GetTempPath(), $"PdiFileTests_{Guid.NewGuid():N}.pdi2");

    [Fact]
    public void Read_DoesNotModifyFile_AndKeepsTextContainingPt()
    {
        var fn = TempFile();
        try
        {
            Assert.True(XYFile.SavePdi2File([Sample()], fn));
            var before = File.ReadAllBytes(fn);

            var d = Assert.Single(XYFile.ReadPdi2File(fn));

            Assert.Equal(before, File.ReadAllBytes(fn));
            Assert.Equal("Sample opt pt", d.Name);
            Assert.Equal("script: adapt, Pt foil", d.Comment);
            Assert.Equal(NumPoints, d.SourceProfile.Pt.Count);
            Assert.Equal(10 + 49 * 0.02, d.SourceProfile.Pt[^1].X, 12);
        }
        finally { File.Delete(fn); }
    }

    [Fact]
    public void Read_LegacyLowercasePtElement_IsStillLoaded_WithoutRewritingFile()
    {
        var fn = TempFile();
        try
        {
            Assert.True(XYFile.SavePdi2File([Sample()], fn));
            var xml = File.ReadAllText(fn);
            Assert.Contains("<Pt>", xml);
            File.WriteAllText(fn, xml.Replace("<Pt>", "<pt>").Replace("</Pt>", "</pt>")); // 旧版の要素名 (小文字の pt) を再現
            var before = File.ReadAllBytes(fn);

            var d = Assert.Single(XYFile.ReadPdi2File(fn));

            Assert.Equal(before, File.ReadAllBytes(fn));
            Assert.Equal(NumPoints, d.SourceProfile.Pt.Count);
            Assert.Equal("Sample opt pt", d.Name);
        }
        finally { File.Delete(fn); }
    }
}
