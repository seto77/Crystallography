using Xunit;

namespace Crystallography.Tests;

/// <summary>261118Cl 新規 (codex2 第 15・16 回、番兵の正準化): <see cref="BetheMethod.SelectCanonicalBeams"/> の合成の候補の試験。
/// 候補の順を全部入れ替えても、独立の参照 (全部を並べる → 鍵ごとに最初 → 境の同点の組) と返る hkl の列が完全に一致すること。
/// 旧の Find_gVectors は N+1 本の部分選択の後に重複除去したので、同じ評価値・同じ接平面の鍵の候補のどれが残るかが部分選択の並びに依存した</summary>
public class BeamCanonicalSelectionTests
{
    private sealed record Candidate((int H, int K, int L) Hkl, double Rating, long Key);

    private static readonly (int, int, int) Zero = (0, 0, 0);

    /// <summary>独立の参照: 000 を先頭に (評価値, 指数) で並べ、鍵ごとに最初を残し、N+1 本目の番兵の組を扱う</summary>
    private static List<(int, int, int)> Reference(IReadOnlyList<Candidate> c, int n, BeamBoundaryTieMode mode)
    {
        var u = c.OrderBy(x => x.Hkl == Zero ? 0 : 1).ThenBy(x => x.Rating).ThenBy(x => x.Hkl.H).ThenBy(x => x.Hkl.K).ThenBy(x => x.Hkl.L)
            .DistinctBy(x => x.Key).ToList();
        if (u.Count <= n) return [.. u.Select(x => x.Hkl)];
        static double tol(double r) => 1E-6 + 1E-11 * Math.Abs(r);
        int j = n;
        while (j >= 1 && Math.Abs(u[j].Rating - u[j - 1].Rating) <= tol(u[j].Rating)) j--;
        if (j == n) return [.. u.Take(n).Select(x => x.Hkl)];
        if (mode == BeamBoundaryTieMode.IncludeGroup)
        {
            double hi = u[n].Rating; int e = n + 1;
            while (e < u.Count && u[e].Rating <= hi + tol(hi)) { hi = Math.Max(hi, u[e].Rating); e++; }
            return [.. u.Take(e).Select(x => x.Hkl)];
        }
        if (j == 0) return u[0].Hkl == Zero ? [Zero] : [];
        return [.. u.Take(j).Select(x => x.Hkl)];
    }

    private static List<(int, int, int)> Select(IReadOnlyList<Candidate> c, int n, BeamBoundaryTieMode mode)
    {
        var (selected, _) = BetheMethod.SelectCanonicalBeams([.. c.Select(x => x.Rating)], [.. c.Select(x => x.Hkl)], [.. c.Select(x => x.Key)], n, mode);
        return [.. selected.Select(i => c[i].Hkl)];
    }

    private static IEnumerable<List<T>> Permutations<T>(List<T> items)
    {
        if (items.Count <= 1) { yield return [.. items]; yield break; }
        for (int i = 0; i < items.Count; i++)
        {
            var rest = items.Where((_, k) => k != i).ToList();
            foreach (var p in Permutations(rest)) { p.Insert(0, items[i]); yield return p; }
        }
    }

    /// <summary>全部の順列で参照と完全に一致すること。期待の列が与えられればそれとも一致すること</summary>
    private static void AssertAllOrders(List<Candidate> c, int n, BeamBoundaryTieMode mode, (int, int, int)[] expected = null)
    {
        var reference = Reference(c, n, mode);
        if (expected != null) Assert.Equal(expected, reference);
        foreach (var p in Permutations(c))
            Assert.Equal(reference, Select(p, n, mode));
    }

    /// <summary>codex2 第 16 回の反例: 000 と同じ正の評価値の X₁ &lt; X₂、Y。X どうしだけ同じ接平面の鍵。N = 2。
    /// 旧は部分選択の先頭に X₂, Y が来ると、残りの X₁ を「既出の鍵」として捨て、X₂ を返した (補充が起きなくても正準でない)</summary>
    [Theory]
    [InlineData(BeamBoundaryTieMode.IncludeGroup)]
    [InlineData(BeamBoundaryTieMode.DropGroup)]
    public void DuplicateKeyAtTheBoundary_KeepsTheCanonicalRepresentative(BeamBoundaryTieMode mode)
    {
        List<Candidate> c = [new(Zero, 0, 0), new((1, 0, 0), 2.0, 10), new((2, 0, 0), 2.0, 10), new((0, 1, 0), 2.0, 20)];
        AssertAllOrders(c, 2, mode, mode == BeamBoundaryTieMode.IncludeGroup ? [Zero, (0, 1, 0), (1, 0, 0)] : [Zero]);
    }

    /// <summary>000 と負の指数の反射がどちらも評価値 0 (ブラッグ条件ちょうど) でも 000 が先頭</summary>
    [Fact]
    public void ZeroBeamIsFirst_EvenWhenAnotherReflectionHasRatingZero()
    {
        List<Candidate> c = [new((-1, 1, 1), 0, 1), new(Zero, 0, 0), new((1, -1, -1), 0.5, 2), new((2, 0, 0), 3.0, 3)];
        AssertAllOrders(c, 3, BeamBoundaryTieMode.IncludeGroup, [Zero, (-1, 1, 1), (1, -1, -1)]);
        var (sel, _) = BetheMethod.SelectCanonicalBeams([0, 0], [(-1, 1, 1), Zero], [1, 0], 5, BeamBoundaryTieMode.IncludeGroup);
        Assert.Equal(1, sel[0]); // 000 の番号
    }

    /// <summary>全部が同点 (000 も評価値 0 の組につながる): IncludeGroup は全部、DropGroup は 000 だけ (旧は先頭の N 本を返した)</summary>
    [Theory]
    [InlineData(BeamBoundaryTieMode.IncludeGroup)]
    [InlineData(BeamBoundaryTieMode.DropGroup)]
    public void AllTied(BeamBoundaryTieMode mode)
    {
        List<Candidate> c = [new(Zero, 0, 0), new((1, 0, 0), 0, 1), new((0, 1, 0), 0, 2), new((0, 0, 1), 0, 3)];
        AssertAllOrders(c, 2, mode, mode == BeamBoundaryTieMode.IncludeGroup ? [Zero, (0, 0, 1), (0, 1, 0), (1, 0, 0)] : [Zero]);
    }

    /// <summary>重複を除くと N 本以下 (補充しても番兵に届かない): 一意なものを全部</summary>
    [Fact]
    public void DuplicatesLeaveNoSentinel_ReturnsAllUnique()
    {
        List<Candidate> c = [new(Zero, 0, 0), new((1, 0, 0), 1, 5), new((-1, 0, 0), 2, 5), new((0, 1, 0), 3, 6), new((0, -1, 0), 4, 6)];
        AssertAllOrders(c, 3, BeamBoundaryTieMode.DropGroup, [Zero, (1, 0, 0), (0, 1, 0)]); // 候補 5 本 > N = 3 だが、一意なのは 3 本 (番兵なし)
        var (_, dup) = BetheMethod.SelectCanonicalBeams([.. c.Select(x => x.Rating)], [.. c.Select(x => x.Hkl)], [.. c.Select(x => x.Key)], 3, BeamBoundaryTieMode.DropGroup);
        Assert.True(dup); // 番兵までの候補で重複除去が効いた (診断の数え方)
    }

    /// <summary>許容の内側でつながる連鎖 (隣どうしは同点、端どうしは同点でない) と、外側で切れる場合</summary>
    [Theory]
    [InlineData(BeamBoundaryTieMode.IncludeGroup)]
    [InlineData(BeamBoundaryTieMode.DropGroup)]
    public void ToleranceChain(BeamBoundaryTieMode mode)
    {
        List<Candidate> c = [new(Zero, 0, 0), new((1, 0, 0), 1.0, 1), new((2, 0, 0), 5.0, 2), new((3, 0, 0), 5.0 + 0.9e-6, 3),
            new((4, 0, 0), 5.0 + 1.8e-6, 4), new((5, 0, 0), 5.0 + 4e-6, 5)];
        AssertAllOrders(c, 3, mode, mode == BeamBoundaryTieMode.IncludeGroup
            ? [Zero, (1, 0, 0), (2, 0, 0), (3, 0, 0), (4, 0, 0)]   // 番兵 (3,0,0) の組は前へ (2,0,0)、後ろへ (4,0,0) までつながり、(5,0,0) で切れる
            : [Zero, (1, 0, 0)]);
    }

    /// <summary>番兵が組をつくらない (N 本目と同点でない): 先頭の N 本</summary>
    [Fact]
    public void SentinelAlone_DropsOnlyTheSentinel()
    {
        List<Candidate> c = [new(Zero, 0, 0), new((1, 0, 0), 1, 1), new((0, 1, 0), 2, 2), new((0, 0, 1), 3, 3)];
        AssertAllOrders(c, 2, BeamBoundaryTieMode.IncludeGroup, [Zero, (1, 0, 0)]);
        AssertAllOrders(c, 2, BeamBoundaryTieMode.DropGroup, [Zero, (1, 0, 0)]);
    }

    [Fact]
    public void EmptyAndZeroOnly()
    {
        Assert.Empty(BetheMethod.SelectCanonicalBeams([], [], [], 3, BeamBoundaryTieMode.IncludeGroup).Selected);
        Assert.Equal([0], BetheMethod.SelectCanonicalBeams([0], [Zero], [0], 3, BeamBoundaryTieMode.DropGroup).Selected);
    }

    /// <summary>無作為の合成の候補 (同点・重複の鍵を多く含む) で参照と一致すること (順は無作為に入れ替える)</summary>
    [Fact]
    public void RandomCandidates_MatchTheReference()
    {
        var rng = new Random(261118);
        for (int trial = 0; trial < 400; trial++)
        {
            int m = rng.Next(1, 30);
            var hkls = new HashSet<(int, int, int)> { Zero };
            while (hkls.Count < m) hkls.Add((rng.Next(-3, 4), rng.Next(-3, 4), rng.Next(-3, 4)));
            var c = hkls.Select(h => new Candidate(h, h == Zero ? 0 : rng.Next(0, 6) * 1.0 + (rng.Next(3) == 0 ? 0.5e-6 : 0), h == Zero ? 0 : rng.Next(1, 12))).ToList();
            c = [.. c.OrderBy(_ => rng.Next())];
            int n = rng.Next(1, m + 2);
            foreach (var mode in new[] { BeamBoundaryTieMode.IncludeGroup, BeamBoundaryTieMode.DropGroup })
                Assert.Equal(Reference(c, n, mode), Select(c, n, mode));
        }
    }
}
