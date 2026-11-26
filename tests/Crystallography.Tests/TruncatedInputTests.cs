using System;
using System.IO;
using System.IO.Compression;
using System.Linq;
using System.Text;
using Xunit;

namespace Crystallography.Tests;

/// <summary>261126Cl 追加 (codex2 の 11-25 の議論 観点 I §2・§3): 途中で切れた TIFF を成功にしない (以前は末尾の画素を 0 で埋めていた)、
/// .rpmp の深さの軸 (有限・厳密に増加) を書く前と読むときに確かめる</summary>
public sealed class TruncatedInputTests : IDisposable
{
    private readonly DirectoryInfo dir = Directory.CreateTempSubdirectory("recipro_truncated_");
    public void Dispose() => dir.Delete(recursive: true);

    /// <summary>II (little endian) の 32 bit 浮動小数点・幅 2 の TIFF を書く。IFD を画素の前に置くので、画素だけを切ったファイルを作れる。
    /// stripTableBeyondEof なら、2 本のストリップの位置の表をファイルの外に置く</summary>
    string WriteTiff(int compression, int declared, int rawCount, int height, bool stripTableBeyondEof = false, int writtenTags = 8)
    {
        string path = Path.Combine(dir.FullName, Guid.NewGuid().ToString("N") + ".tif");
        byte[] raw = [0, 0, 128, 63, 0, 0, 0, 64]; // float の 1 と 2
        using var w = new BinaryWriter(File.Create(path));
        void Tag(ushort id, ushort type, int count, int value) { w.Write(id); w.Write(type); w.Write(count); w.Write(value); }
        const ushort tags = 8;
        int pixelOffset = 8 + 2 + tags * 12 + 4;
        w.Write((ushort)0x4949); w.Write((ushort)42); w.Write(8); w.Write(tags);
        var entries = new Action[]
        {
            () => Tag(256, 4, 1, 2), () => Tag(257, 4, 1, height), () => Tag(258, 3, 1, 32), () => Tag(259, 3, 1, compression),
            () => { if (stripTableBeyondEof) Tag(273, 4, 2, 1_000_000); else Tag(273, 4, 1, pixelOffset); },
            () => Tag(278, 4, 1, 1), () => Tag(279, 4, 1, declared), () => Tag(339, 3, 1, 3),
        };
        foreach (var e in entries.Take(writtenTags)) e();
        if (writtenTags < tags) return path; // IFD の表の途中で切れたファイル
        w.Write(0); // 次の IFD はない
        if (compression == 32773) w.Write((byte)(rawCount - 1)); // PackBits のリテラル列
        w.Write(raw, 0, rawCount);
        return path;
    }

    [Theory]
    [InlineData(1, 8, 8, 1)]     // 非圧縮の対照
    [InlineData(32773, 9, 8, 1)] // PackBits の対照
    public void CompleteTiffIsRead(int compression, int declared, int rawCount, int height)
        => Assert.Equal(new double[] { 1, 2 }, new Tiff.Loader(WriteTiff(compression, declared, rawCount, height)).Images.Single().Value);

    [Theory]
    [InlineData(1, 8, 5, 1)]     // 2 つ目の float の途中でファイルが尽きる (以前は 0 で埋めて成功)
    [InlineData(1, 4, 4, 1)]     // 宣言どおりに読めても画素数に足りない (以前は残りが 0 のまま成功)
    [InlineData(1, 8, 8, 2)]     // 2 行目のストリップが無い
    [InlineData(32773, 5, 4, 1)] // 圧縮の列は完結しているが、展開すると 1 画素だけ
    [InlineData(1, 9, 8, 1)]     // 宣言の範囲がファイルの外へはみ出す (読める分だけで画素はそろう)
    [InlineData(32773, 10, 8, 1)] // 同上 (PackBits)
    public void TruncatedTiffFails(int compression, int declared, int rawCount, int height)
        => Assert.Throws<EndOfStreamException>(() => new Tiff.Loader(WriteTiff(compression, declared, rawCount, height)));

    [Fact]
    public void TruncatedStripTableFails()
        => Assert.Throws<EndOfStreamException>(() => new Tiff.Loader(WriteTiff(1, 8, 8, 2, stripTableBeyondEof: true)));

    /// <summary>2 本の帯 (高さ 2・1 行ずつ) の位置と長さの表を、画素の後ろ (ファイルの末尾) に置く。dropLast なら最後の 1 バイト (長さ 8 の上位の 0) を落とす。
    /// 以前は外部の表を 0 で埋めて読んだので、落としても長さが 8 のままで成功していた</summary>
    string WriteTwoStripTiffWithTrailingTables(bool dropLast)
    {
        string path = Path.Combine(dir.FullName, Guid.NewGuid().ToString("N") + ".tif");
        using (var w = new BinaryWriter(File.Create(path)))
        {
            void Tag(ushort id, ushort type, int count, int value) { w.Write(id); w.Write(type); w.Write(count); w.Write(value); }
            const ushort tags = 8;
            int pixelOffset = 8 + 2 + tags * 12 + 4, offsetsTable = pixelOffset + 16, countsTable = offsetsTable + 8;
            w.Write((ushort)0x4949); w.Write((ushort)42); w.Write(8); w.Write(tags);
            Tag(256, 4, 1, 2); Tag(257, 4, 1, 2); Tag(258, 3, 1, 32); Tag(259, 3, 1, 1);
            Tag(273, 4, 2, offsetsTable); Tag(278, 4, 1, 1); Tag(279, 4, 2, countsTable); Tag(339, 3, 1, 3);
            w.Write(0);
            foreach (var v in new float[] { 1, 2, 3, 4 }) w.Write(v);
            w.Write(pixelOffset); w.Write(pixelOffset + 8);
            w.Write(8); w.Write(8);
        }
        if (dropLast) using (var f = new FileStream(path, FileMode.Open)) f.SetLength(f.Length - 1);
        return path;
    }

    [Fact]
    public void TruncatedTrailingStripTableFails()
    {
        Assert.Equal(new double[] { 1, 2, 3, 4 }, new Tiff.Loader(WriteTwoStripTiffWithTrailingTables(false)).Images.Single().Value); // 対照
        Assert.Throws<EndOfStreamException>(() => new Tiff.Loader(WriteTwoStripTiffWithTrailingTables(true)));
    }

    /// <summary>MM (big endian) の 16 bit・幅 2 の TIFF。以前は MM の整数だけ末尾を 0 で埋めて成功にしていた</summary>
    [Theory]
    [InlineData(4, true)]
    [InlineData(3, false)]
    public void BigEndianTiffRequiresEveryPixel(int rawCount, bool valid)
    {
        string path = Path.Combine(dir.FullName, Guid.NewGuid().ToString("N") + ".tif");
        using (var w = new BinaryWriter(File.Create(path)))
        {
            void U16(int v) { w.Write((byte)(v >> 8)); w.Write((byte)v); }
            void U32(int v) { U16(v >> 16); U16(v & 0xFFFF); }
            void Tag(ushort id, ushort type, int value) { U16(id); U16(type); U32(1); if (type == 3) { U16(value); U16(0); } else U32(value); }
            const ushort tags = 7;
            w.Write((byte)'M'); w.Write((byte)'M'); U16(42); U32(8); U16(tags);
            Tag(256, 3, 2); Tag(257, 3, 1); Tag(258, 3, 16); Tag(259, 3, 1);
            Tag(273, 4, 8 + 2 + tags * 12 + 4); Tag(278, 3, 1); Tag(279, 4, 4);
            U32(0);
            w.Write(new byte[] { 0x12, 0x34, 0x56, 0x78 }, 0, rawCount);
        }
        if (valid) Assert.Equal(new double[] { 0x1234, 0x5678 }, new Tiff.Loader(path).Images.Single().Value);
        else Assert.Throws<EndOfStreamException>(() => new Tiff.Loader(path));
    }

    [Fact]
    public void TruncatedImageFileDirectoryFails()
        => Assert.Throws<EndOfStreamException>(() => new Tiff.Loader(WriteTiff(1, 8, 8, 1, writtenTags: 3)));

    [Fact]
    public void PartialBigEndianSampleIsNotZeroPadded()
    {
        using var input = new MemoryStream([0x12, 0x34, 0x56]);
        double[] pixels = [-1, -1];
        Assert.Throws<EndOfStreamException>(() => ImageIO.ReadPixels<ushort>(input, pixels, bigEndian: true));
        Assert.Equal(new double[] { 0x1234, -1 }, pixels);
    }

    static MasterPattern Pattern(double[] depths) => new(2, [20.0], depths, new float[depths.Length][], new float[depths.Length][], MasterPattern.Types.Square);

    [Theory]
    [InlineData(1, 1)]
    [InlineData(2, 1)]
    [InlineData(double.NaN, 1)]
    [InlineData(1, double.PositiveInfinity)]
    [InlineData(double.NegativeInfinity, 1)]
    public void InvalidDepthsAreNotSavedAndKeepTheExistingFile(double first, double second)
    {
        string path = Path.Combine(dir.FullName, "keep.rpmp");
        byte[] original = [1, 2, 3];
        File.WriteAllBytes(path, original);
        Assert.Throws<ArgumentException>(() => MasterPatternFile.Save(Pattern([first, second]), path, null));
        Assert.Equal(original, File.ReadAllBytes(path));
        Assert.Throws<ArgumentException>(() => MasterPatternFile.SaveVersion1(Pattern([first, second]), path, null));
        Assert.Equal(original, File.ReadAllBytes(path));
        Assert.Single(dir.GetFiles()); // 一時ファイルも残さない
    }

    [Theory]
    [InlineData(new double[0])]
    [InlineData(new double[] { 0 })]
    [InlineData(new double[] { 0, 2, 5 })]
    public void ValidDepthAxesRoundTrip(double[] depths)
    {
        string path = Path.Combine(dir.FullName, "ok.rpmp");
        MasterPatternFile.Save(Pattern(depths), path, null);
        Assert.Equal(depths, MasterPatternFile.Load(path).MasterPattern.Depths);
        MasterPatternFile.SaveVersion1(Pattern(depths), path, null);
        Assert.Equal(depths, MasterPatternFile.Load(path).MasterPattern.Depths);
    }

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    public void ReadersRejectDuplicateDepthInAnOtherwiseValidFile(int version)
    {
        string path = Path.Combine(dir.FullName, "axis.rpmp");
        var mp = Pattern([0, 2]);
        if (version == 1) MasterPatternFile.SaveVersion1(mp, path, null);
        else MasterPatternFile.Save(mp, path, null);
        Assert.Equal(mp.Depths, MasterPatternFile.Load(path).MasterPattern.Depths); // 書き換える前は読める
        Assert.Empty(MasterPatternFile.ReadElectronChunks(path));
        byte[] bytes = File.ReadAllBytes(path);
        int headerLength = BitConverter.ToInt32(bytes, 8), bodyOffset = 12 + headerLength;
        if (version == 2)
        {
            var header = bytes.AsSpan(12, headerLength);
            ReadOnlySpan<byte> needle = "\"Depths\":[0,2]"u8;
            int pos = header.IndexOf(needle);
            Assert.True(pos >= 0);
            header[pos + needle.Length - 2] = (byte)'0'; // [0,2] → [0,0]。見出しの長さ・区画の位置・CRC は保つ
            File.WriteAllBytes(path, bytes);
            Assert.Contains("Depths[1]", Assert.Throws<InvalidDataException>(() => MasterPatternFile.ReadInfo(path)).Message);
        }
        else
        {
            using var body = new MemoryStream();
            using (var encoded = new MemoryStream(bytes, bodyOffset, bytes.Length - bodyOffset, false))
            using (var z = new ZLibStream(encoded, CompressionMode.Decompress)) z.CopyTo(body);
            body.Position = 4 * sizeof(int) + 2 * sizeof(double); // 格子の 4 整数、エネルギー 1 個、深さ 1 個の次
            using (var bw = new BinaryWriter(body, Encoding.UTF8, leaveOpen: true)) bw.Write(0.0);
            using var rewritten = new MemoryStream();
            rewritten.Write(bytes, 0, bodyOffset);
            body.Position = 0;
            using (var z = new ZLibStream(rewritten, CompressionLevel.Fastest, leaveOpen: true)) body.CopyTo(z);
            File.WriteAllBytes(path, rewritten.ToArray());
            Assert.NotNull(MasterPatternFile.ReadInfo(path)); // v1 の見出しは深さの列を持たない
        }
        Assert.Contains("Depths[1]", Assert.Throws<InvalidDataException>(() => MasterPatternFile.Load(path, false)).Message);
        Assert.Contains("Depths[1]", Assert.Throws<InvalidDataException>(() => MasterPatternFile.ReadElectronChunks(path).ToArray()).Message);
    }
}
