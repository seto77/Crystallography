using Xunit;

namespace Crystallography.Tests;

/// <summary>261124Cl 追加: PDIndexer の査読から引き継いだ共有ライブラリの修正 (C-02, C-14, C-16, C-21) の回帰テスト。</summary>
public class SharedLibReviewTests
{
    /// <summary>C-02: 短い試料の後ろの行 (,,,4,400,,) の値が、先頭の試料に入らない</summary>
    [Fact]
    public void CombinedCsv_ShortSampleKeepsItsColumns()
    {
        var path = Path.Combine(Path.GetTempPath(), $"recipro_c02_{Guid.NewGuid():N}.csv");
        File.WriteAllLines(path, ["A,,,B,,,", "X,Y,,X,Y,,", "1,10,,1,100,,", "2,20,,2,200,,", "3,30,,3,300,,", ",,,4,400,,"]);
        try
        {
            var dp = XYFile.ReadCSVFile(path);
            Assert.Equal(2, dp.Length);
            Assert.Equal("A", dp[0].Name);
            Assert.Equal(3, dp[0].SourceProfile.Pt.Count);
            Assert.Equal(4, dp[1].SourceProfile.Pt.Count);
            Assert.Equal(400, dp[1].SourceProfile.Pt[^1].Y);
        }
        finally { File.Delete(path); }
    }

    /// <summary>C-02: データが 1 行だけのファイルと、小文字の見出し (IPAnalyzer の結合 CSV) も読める</summary>
    [Fact]
    public void CombinedCsv_OneRowAndLowercaseHeader()
    {
        var path = Path.Combine(Path.GetTempPath(), $"recipro_c02b_{Guid.NewGuid():N}.csv");
        File.WriteAllLines(path, ["A,,,,,", "x,y,,x,y,", "1,10,,2,20,"]);
        try
        {
            var dp = XYFile.ReadCSVFile(path);
            Assert.Equal(2, dp.Length);
            Assert.Single(dp[0].SourceProfile.Pt);
            Assert.Single(dp[1].SourceProfile.Pt);
            Assert.False(string.IsNullOrEmpty(dp[1].Name)); // 空の名前は補う
        }
        finally { File.Delete(path); }
    }

    /// <summary>C-14: 両側から点を集めるとき、点が足りなくても終わる。0 点は NaN、1 点は定数</summary>
    [Fact]
    public void GetValues_TerminatesWithFewPoints()
    {
        var p = new Profile();
        foreach (var x in new[] { 0.0, 1, 3, 4 })
            p.Pt.Add(new PointD(x, x * x));
        var v = p.GetValues([2.0], 6, 2, true, true); // 両側で 6 点を求めるが 4 点しかない
        Assert.True(double.IsFinite(v[0]));

        Assert.True(double.IsNaN(new Profile().GetValues([1.0], 4, 2, true)[0]));
        var one = new Profile();
        one.Pt.Add(new PointD(1, 5));
        Assert.Equal(5, one.GetValues([3.0], 4, 2, true)[0]);
    }

    /// <summary>C-16: q = 0 で NaN にならず、q → 0 の極限と連続</summary>
    [Fact]
    public void MieGruneisen_QZeroIsTheLimit()
    {
        double p0 = EOS.MieGruneisen(4, 1, 170, 2.9, 0, 300, 67.85, 1000, 60);
        double pSmall = EOS.MieGruneisen(4, 1, 170, 2.9, 1e-7, 300, 67.85, 1000, 60);
        Assert.True(double.IsFinite(p0));
        Assert.Equal(pSmall, p0, 4);
    }

    /// <summary>C-21: 強度と誤差を一組として変換する。対数では σ/(I ln10)、I ≤ 0 の行は両方から除く</summary>
    [Fact]
    public void ConvertedProfile_LogKeepsIntensityErrorPairs()
    {
        var dp = new DiffractionProfile2();
        dp.SourceProfile.Pt.AddRange([new PointD(1, 0), new PointD(2, 100)]);
        dp.SourceProfile.Err.AddRange([new PointD(1, 1), new PointD(2, 10)]);
        dp.IsCPS = false;
        dp.IsLogIntensity = true;
        dp.SetConvertedProfile(dp.SrcProperty);
        Assert.Single(dp.ConvertedProfile.Pt);
        Assert.Single(dp.ConvertedProfile.Err);
        Assert.Equal(2, dp.ConvertedProfile.Pt[0].X);
        Assert.Equal(2, dp.ConvertedProfile.Err[0].X);
        Assert.Equal(2, dp.ConvertedProfile.Pt[0].Y, 12);
        Assert.Equal(10 / (100 * Math.Log(10)), dp.ConvertedProfile.Err[0].Y, 12);
    }
}
