using System;
using System.IO;
using System.IO.Compression;
using System.Linq;
using Xunit;

namespace Crystallography.Tests;

/// <summary>261124Cl 追加 (作者決定 11-24): マスターパターンのファイルの形式 v2 (区画・並べ替え + Brotli・CRC32・結晶の同梱・物理の版) の試験。
/// v2 の往復がビットで同じこと、v1 のファイルを読んだ結果と v2 の結果がビットで同じこと、電子の塊の読み方 (ReadElectronChunks) が
/// Load と同じこと、区画の化けを CRC32 で見つけること、同梱の結晶が同じ指紋に戻ることを確かめる</summary>
public class MasterPatternFileV2Tests
{
    static MasterPattern Pattern(int seed)
    {
        var rnd = new Random(seed);
        int g = 16;
        double[] es = [20.0, 15.0], ds = [1.0, 5.0, 10.0];
        float[][] Planes() => Enumerable.Range(0, es.Length * ds.Length)
            .Select(k => k == 3 ? null : Enumerable.Range(0, g * g).Select(i => (float)(Math.Sin(i * 0.1 + k) + 1 + rnd.NextDouble() * 1e-3)).ToArray()).ToArray(); // 平面 3 は欠落
        return new MasterPattern(g, es, ds, Planes(), Planes(), MasterPattern.Types.Square) { BackSourceAxis = BetheMethod.BackSourceAxisIncident };
    }

    static EbsdBackscatteredElectron[] Electrons(int n, int seed)
    {
        var rnd = new Random(seed);
        double R() => rnd.NextDouble() * 20 - 10;
        // 非弾性の事象が無い電子は MC と同じく NaN を入れる (Vector3d の == は NaN を等しいとしないので、比べるのは Bits で)
        return Enumerable.Range(0, n).Select(i => i % 3 == 0
            ? new EbsdBackscatteredElectron(R(), new OpenTK.Mathematics.Vector3d(R(), R(), R()), new PointD(R(), R()), R(), R(),
                true, R(), R(), R(), new OpenTK.Mathematics.Vector3d(R(), R(), R()), i % 5 == 0, R())
            : new EbsdBackscatteredElectron(R(), new OpenTK.Mathematics.Vector3d(R(), R(), R()), new PointD(R(), R()), R(), R(),
                false, double.NaN, double.NaN, double.NaN, new OpenTK.Mathematics.Vector3d(double.NaN, double.NaN, double.NaN), i % 5 == 0, R())).ToArray();
    }

    /// <summary>電子の全成分のビット列 (NaN も同じビットなら同じ)</summary>
    static long[] Bits(EbsdBackscatteredElectron[] es) => es.SelectMany(e => new[] {
        e.Depth, e.Vec.X, e.Vec.Y, e.Vec.Z, e.Position.X, e.Position.Y, e.Energy, e.TotalEnergyLoss, e.LastInelasticDepth, e.LastInelasticEnergyBeforeLoss,
        e.LastInelasticEnergyAfterLoss, e.LastInelasticDirection.X, e.LastInelasticDirection.Y, e.LastInelasticDirection.Z, e.LastDecoherenceDepth }
        .Select(BitConverter.DoubleToInt64Bits).Append((e.HasLastInelasticEvent ? 1L : 0) | (e.HasLastDecoherenceEvent ? 2L : 0))).ToArray();

    static string TempPath() => Path.Combine(Path.GetTempPath(), $"rpmp_v2_{Guid.NewGuid():N}.rpmp");

    static void AssertSamePlanes(MasterPattern a, MasterPattern b)
    {
        Assert.Equal(a.GridSize, b.GridSize);
        Assert.Equal(a.GridType, b.GridType);
        Assert.Equal(a.Energies, b.Energies);
        Assert.Equal(a.Depths, b.Depths);
        foreach (var (pa, pb) in new[] { (a.PositivePlanes, b.PositivePlanes), (a.NegativePlanes, b.NegativePlanes) })
            for (int k = 0; k < pa.Length; k++)
            {
                if (pa[k] is not { Length: > 0 }) { Assert.True(pb[k] is not { Length: > 0 }); continue; }
                Assert.True(pa[k].AsSpan().SequenceEqual(pb[k]), $"plane {k}");
            }
    }

    /// <summary>v2 の往復: 平面はビットで同じ、電子は float32 に丸めた値 (v1 と同じ) で同じ。電子は 3 区画にまたがる</summary>
    [Theory]
    [InlineData(CompressionLevel.Fastest)]
    [InlineData(CompressionLevel.NoCompression)]
    [InlineData(CompressionLevel.SmallestSize)]
    public void Version2_RoundTrip(CompressionLevel level)
    {
        var mp = Pattern(1);
        var es = Electrons(2 * MasterPatternFile.ElectronsPerBlock + 1234, 2);
        string p2 = TempPath(), p1 = TempPath();
        try
        {
            MasterPatternFile.Save(mp, p2, new MasterPatternFileInfo { McSampleTiltRad = 1.2, McEnergyThresholdKeV = 5, McIncidentCount = 100 }, es, level);
            var info = MasterPatternFile.ReadInfo(p2);
            Assert.Equal(2, info.FormatVersion);
            Assert.Equal(mp.GridSize, info.GridSize);
            Assert.Equal(mp.Energies, info.Energies);
            Assert.Equal(es.Length, info.BackscatteredElectronCount);
            var v2 = MasterPatternFile.Load(p2);
            AssertSamePlanes(mp, v2.MasterPattern);
            Assert.Equal(BetheMethod.BackSourceAxisIncident, v2.MasterPattern.BackSourceAxis);

            // 同じ中身を v1 で書いて読んだ結果と、ビットで同じ (電子の丸めも同じ)
            MasterPatternFile.SaveVersion1(mp, p1, new MasterPatternFileInfo { McSampleTiltRad = 1.2, McEnergyThresholdKeV = 5 }, es);
            Assert.Equal(1, MasterPatternFile.ReadInfo(p1).FormatVersion);
            var v1 = MasterPatternFile.Load(p1);
            AssertSamePlanes(v1.MasterPattern, v2.MasterPattern);
            Assert.Equal(Bits(v1.BackscatteredElectrons), Bits(v2.BackscatteredElectrons));
            Assert.Equal(es.Length, v2.BackscatteredElectrons.Length);

            // 塊の読み方は Load と同じ (区画の境をまたぐ大きさで)
            Assert.Equal(Bits(v2.BackscatteredElectrons), Bits(MasterPatternFile.ReadElectronChunks(p2, 1000).SelectMany(c => c).ToArray()));
            Assert.Equal(Bits(v1.BackscatteredElectrons), Bits(MasterPatternFile.ReadElectronChunks(p1, 1000).SelectMany(c => c).ToArray()));

            // 電子を読まない読み方
            var noE = MasterPatternFile.Load(p2, loadElectrons: false);
            Assert.Empty(noE.BackscatteredElectrons);
            AssertSamePlanes(mp, noE.MasterPattern);
        }
        finally { File.Delete(p2); File.Delete(p1); }
    }

    /// <summary>なめらかな平面は v2 で v1 より小さくなる (並べ替え + Brotli の効き)</summary>
    [Fact]
    public void Version2_IsSmallerThanVersion1_ForSmoothPlanes()
    {
        var mp = Pattern(3);
        string p2 = TempPath(), p1 = TempPath();
        try
        {
            MasterPatternFile.Save(mp, p2, null);
            MasterPatternFile.SaveVersion1(mp, p1, null);
            Assert.True(new FileInfo(p2).Length < new FileInfo(p1).Length, $"v2 {new FileInfo(p2).Length} B, v1 {new FileInfo(p1).Length} B");
        }
        finally { File.Delete(p2); File.Delete(p1); }
    }

    /// <summary>区画の中身が化けたら、CRC32 か展開の失敗で InvalidDataException になる (黙って違う値を返さない)</summary>
    [Fact]
    public void Version2_DetectsCorruptBlock()
    {
        var mp = Pattern(4);
        var path = TempPath();
        try
        {
            MasterPatternFile.Save(mp, path, null, compression: CompressionLevel.NoCompression);
            var bytes = File.ReadAllBytes(path);
            bytes[^100] ^= 0x40; // 最後の区画 (−Z の平面) の中ほど
            File.WriteAllBytes(path, bytes);
            Assert.Throws<InvalidDataException>(() => MasterPatternFile.Load(path));
        }
        finally { File.Delete(path); }
    }

    /// <summary>同梱の結晶が同じ指紋に戻り、構築の条件からは今の物理の版と結晶が入る</summary>
    [Fact]
    public void EmbeddedCrystal_RestoresTheSameFingerprint()
    {
        var si = ConvertCrystalData.ConvertFromCifText(StemDepthPhaseTests.SiCifText).ToCrystal();
        var info = MasterPatternFileInfo.FromBuildRequest(new EBSD.MasterPatternBuildRequest(si, 32, [20.0], [10.0], 16), 20);
        Assert.Equal(MasterPatternFile.PhysicsVersion, info.PhysicsVersion);
        var restored = info.GetCrystal();
        Assert.NotNull(restored);
        Assert.Equal(MasterPatternFile.CrystalFingerprint(si), MasterPatternFile.CrystalFingerprint(restored));
        Assert.Null(MasterPatternFile.CheckCompatibility(info, restored, 20));

        // 見出しを通しても同じ (保存と読み込み)
        var path = TempPath();
        try
        {
            MasterPatternFile.Save(Pattern(5), path, info);
            var back = MasterPatternFile.ReadInfo(path).GetCrystal();
            Assert.Equal(MasterPatternFile.CrystalFingerprint(si), MasterPatternFile.CrystalFingerprint(back));
        }
        finally { File.Delete(path); }

        // 結晶の無い見出し (v1 のファイル) は null
        Assert.Null(new MasterPatternFileInfo().GetCrystal());
    }
}
