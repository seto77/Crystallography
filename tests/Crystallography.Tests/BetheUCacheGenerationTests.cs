using Xunit;

namespace Crystallography.Tests;

/// <summary>261118Cl 新規 (codex2 第 4 回): getU のキャッシュ (<see cref="BetheMethod.UCache"/>) の世代の切り離しを、決まった順に逐次で検査する。
/// Clear の前に取った入れ物へ後から格納しても (Clear と重なった評価が終わったときに起きる)、今の世代からは見えないこと。
/// 以前の形 (一覧を空にして直前を null に・格納のときに入れ物を取り直す) では、Clear の後に古い入れ物が直前として戻りえた。</summary>
public class BetheUCacheGenerationTests
{
    private static BetheMethod.UContext Ctx(double kV) =>
        new(kV, double.NaN, double.NaN, 0, 0, ElasticIonModel.Neutral, AbsorptionElasticModel.TemariHybrid);

    [Fact]
    public void BucketTakenBeforeClear_IsNotVisibleAfterClear()
    {
        var cache = new BetheMethod.UCache();
        var oldBucket = cache.Get(Ctx(200));          // 評価の入口で取った入れ物
        cache.Clear();                                // 別のスレッドが Clear
        var newBucket = cache.Get(Ctx(200));          // Clear の後の評価
        oldBucket.Map.TryAdd((1, 2), (1, 1));         // Clear の前に始まった評価が、取ってあった入れ物へ格納
        Assert.NotSame(oldBucket, newBucket);
        Assert.False(newBucket.Map.ContainsKey((1, 2)));
        Assert.False(cache.Get(Ctx(200)).Map.ContainsKey((1, 2))); // 「直前の入れ物」を経由しても見えない
        Assert.Same(newBucket, cache.Get(Ctx(200)));
    }

    [Fact]
    public void DifferentContexts_GetDifferentBuckets_AndTheLastBucketHintIsChecked()
    {
        var cache = new BetheMethod.UCache();
        var a = cache.Get(Ctx(200));
        a.Map.TryAdd((1, 2), (1, 1));
        var b = cache.Get(Ctx(300));                  // 直前の入れ物と条件が違う
        Assert.NotSame(a, b);
        Assert.False(b.Map.ContainsKey((1, 2)));
        Assert.Same(a, cache.Get(Ctx(200)));          // 一覧から同じ入れ物が戻る
    }
}
