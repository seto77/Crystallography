using System.Drawing;
using System.IO;
using Xunit;

namespace Crystallography.Tests;

/// <summary>261124Cl 追加 (IPAnalyzer の codex2 査読の所見 9・10、作者決定 11-24): ImageIO.ReadImage が失敗したときに、読む前の状態 (Ring の画像と連番データ) が残ること。
/// Ring は static の状態を持つので、Ring の試験と同じ Collection で並列に走らせない</summary>
[Collection("Ring static state")]
public class ImageIOReadStateTests
{
    /// <summary>所見 9: 形式が分からないファイルを読んでも、表示中の画像・連番データ・原画像は変わらない</summary>
    [Fact]
    public void UnknownFormatKeepsTheCurrentImage()
    {
        var path = Path.Combine(Path.GetTempPath(), $"recipro_imageio_{System.Guid.NewGuid():N}.unknownext");
        File.WriteAllBytes(path, [1, 2, 3, 4]);
        var oldIntensity = Ring.Intensity; var oldSize = Ring.SrcImgSize; var oldOriginal = Ring.IntensityOriginal;
        var oldNames = Ring.SequentialImageNames; var oldSeq = Ring.SequentialImageIntensities; var oldComments = Ring.Comments;
        try
        {
            double[] img = [1, 2, 3];
            double[] original = [7, 8, 9];
            Ring.Intensity = img; Ring.SrcImgSize = new Size(3, 1); Ring.IntensityOriginal = original;
            Ring.SequentialImageNames = ["a"]; Ring.SequentialImageIntensities = [img]; Ring.Comments = "keep";

            Assert.False(ImageIO.ReadImage(path));

            Assert.Same(img, Ring.Intensity);
            Assert.Equal(new double[] { 1, 2, 3 }, Ring.Intensity);
            Assert.Equal(new Size(3, 1), Ring.SrcImgSize);
            Assert.Same(original, Ring.IntensityOriginal); // 旧は (読みかけの) Ring.Intensity で上書きされることがあった
            Assert.Equal(["a"], Ring.SequentialImageNames); // 旧は形式の判定より前に null にしていた
            Assert.Single(Ring.SequentialImageIntensities);
            Assert.Equal("keep", Ring.Comments);
        }
        finally
        {
            File.Delete(path);
            (Ring.Intensity, Ring.SrcImgSize, Ring.IntensityOriginal, Ring.SequentialImageNames, Ring.SequentialImageIntensities, Ring.Comments)
                = (oldIntensity, oldSize, oldOriginal, oldNames, oldSeq, oldComments);
        }
    }

    /// <summary>所見 10: 読込本体の無い HDF5 (.h5) は読める形式に数えない。NeXus (.nxs) は読める</summary>
    [Fact]
    public void Hdf5IsNotListedAsReadable()
    {
        Assert.False(ImageIO.IsReadable(".h5"));
        Assert.DoesNotContain("*.h5", ImageIO.FilterString);
        Assert.True(ImageIO.IsReadable(".nxs"));
        Assert.True(ImageIO.IsReadable("tif"));
    }
}
