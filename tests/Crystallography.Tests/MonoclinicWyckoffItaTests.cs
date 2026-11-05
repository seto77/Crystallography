using Xunit;
using static Crystallography.Tests.SymmetryReference;

namespace Crystallography.Tests;

/// <summary>261105Cl 新規: 単斜晶 (No. 3–15) の全設定で、Wyckoff 文字が ITA (International Tables for Crystallography Vol. A, 2016) と一致すること。
/// b 軸・c 軸の設定 (セル選択 1–3) は ITA に印刷された座標 (下の Printed、頁は本の頁番号) と、a 軸の設定と c-ba・ba-c・-acb の設定は
/// ITA の変換 (Table 1.5.1.1 の b 軸→a 軸、Table 1.5.4.4 の設定記号) で印刷された設定から写したものと比べる。
/// cctbx (基準データ) はセル選択の変換を ITA と別に取るため、非標準設定 27 個で文字が ITA と食い違う (README の「既知の例外」)。
/// ITA に印刷された特殊位置の代表座標で、自由な座標は一般の値にして GetEquivalentAtomsPosition に渡し、返る文字を確かめる。</summary>
public class MonoclinicWyckoffItaTests
{
    /// <summary>ITA に印刷された特殊位置の最初の座標 (文字 座標; …)。通し番号 = spglib の Hall 番号。頁は ITA (2016) Vol. A の本の頁</summary>
    private static readonly Dictionary<int, (int Page, string Positions)> Printed = new()
    {
        [3] = (196, "d 1/2,y,1/2; c 1/2,y,0; b 0,y,1/2; a 0,y,0"),                        // P 1 2 1
        [4] = (197, "d 1/2,1/2,z; c 0,1/2,z; b 1/2,0,z; a 0,0,z"),                        // P 1 1 2
        [18] = (208, "b x,1/2,z; a x,0,z"),                                               // P 1 m 1
        [19] = (209, "b x,y,1/2; a x,y,0"),                                               // P 1 1 m
        [57] = (230, "n x,1/2,z; m x,0,z; l 1/2,y,1/2; k 0,y,1/2; j 1/2,y,0; i 0,y,0; h 1/2,1/2,1/2; g 1/2,0,1/2; f 0,1/2,1/2; e 1/2,1/2,0; d 1/2,0,0; c 0,0,1/2; b 0,1/2,0; a 0,0,0"),
        [58] = (232, "n x,y,1/2; m x,y,0; l 1/2,1/2,z; k 1/2,0,z; j 0,1/2,z; i 0,0,z; h 1/2,1/2,1/2; g 1/2,1/2,0; f 1/2,0,1/2; e 0,1/2,1/2; d 0,1/2,0; c 1/2,0,0; b 0,0,1/2; a 0,0,0"),
        [60] = (234, "e x,1/4,z; d 1/2,0,1/2; c 0,0,1/2; b 1/2,0,0; a 0,0,0"),          // P 1 21/m 1
        [61] = (235, "e x,y,1/4; d 1/2,1/2,0; c 1/2,0,0; b 0,1/2,0; a 0,0,0"),          // P 1 1 21/m
        [63] = (238, "i x,0,z; h 0,y,1/2; g 0,y,0; f 1/4,1/4,1/2; e 1/4,1/4,0; d 0,1/2,1/2; c 0,0,1/2; b 0,1/2,0; a 0,0,0"),        // C 1 2/m 1
        [64] = (239, "i x,0,z; h 1/2,y,1/2; g 0,y,0; f 1/2,1/4,3/4; e 0,1/4,1/4; d 1/2,1/2,1/2; c 1/2,0,1/2; b 0,1/2,0; a 0,0,0"),  // A 1 2/m 1
        [65] = (239, "i x,0,z; h 1/2,y,0; g 0,y,0; f 1/4,1/4,3/4; e 3/4,1/4,3/4; d 1/2,1/2,0; c 1/2,0,0; b 0,1/2,0; a 0,0,0"),      // I 1 2/m 1
        [66] = (242, "i x,y,0; h 1/2,0,z; g 0,0,z; f 1/2,1/4,1/4; e 0,1/4,1/4; d 1/2,0,1/2; c 1/2,0,0; b 0,0,1/2; a 0,0,0"),        // A 1 1 2/m
        [67] = (243, "i x,y,0; h 1/2,1/2,z; g 0,0,z; f 3/4,1/2,1/4; e 1/4,0,1/4; d 1/2,1/2,1/2; c 1/2,1/2,0; b 0,0,1/2; a 0,0,0"),  // B 1 1 2/m
        [68] = (243, "i x,y,0; h 0,1/2,z; g 0,0,z; f 3/4,1/4,1/4; e 3/4,3/4,1/4; d 0,1/2,1/2; c 0,1/2,0; b 0,0,1/2; a 0,0,0"),      // I 1 1 2/m
        [72] = (246, "f 1/2,y,1/4; e 0,y,1/4; d 1/2,0,0; c 0,1/2,0; b 1/2,1/2,0; a 0,0,0"),     // P 1 2/c 1
        [73] = (247, "f 3/4,y,1/4; e 3/4,y,3/4; d 0,0,1/2; c 0,1/2,0; b 0,1/2,1/2; a 0,0,0"),   // P 1 2/n 1
        [74] = (247, "f 3/4,y,1/2; e 1/4,y,0; d 1/2,0,1/2; c 0,1/2,0; b 1/2,1/2,1/2; a 0,0,0"), // P 1 2/a 1
        [75] = (250, "f 1/4,1/2,z; e 1/4,0,z; d 0,1/2,0; c 0,0,1/2; b 0,1/2,1/2; a 0,0,0"),     // P 1 1 2/a
        [76] = (251, "f 1/4,3/4,z; e 3/4,3/4,z; d 1/2,0,0; c 0,0,1/2; b 1/2,0,1/2; a 0,0,0"),   // P 1 1 2/n
        [77] = (251, "f 1/2,3/4,z; e 0,1/4,z; d 1/2,1/2,0; c 0,0,1/2; b 1/2,1/2,1/2; a 0,0,0"), // P 1 1 2/b
        [81] = (254, "d 1/2,0,1/2; c 0,0,1/2; b 1/2,0,0; a 0,0,0"),   // P 1 21/c 1
        [82] = (255, "d 1/2,0,0; c 1/2,0,1/2; b 0,0,1/2; a 0,0,0"),   // P 1 21/n 1
        [83] = (255, "d 0,0,1/2; c 1/2,0,0; b 1/2,0,1/2; a 0,0,0"),   // P 1 21/a 1
        [84] = (258, "d 1/2,1/2,0; c 1/2,0,0; b 0,1/2,0; a 0,0,0"),   // P 1 1 21/a
        [85] = (259, "d 0,1/2,0; c 1/2,1/2,0; b 1/2,0,0; a 0,0,0"),   // P 1 1 21/n
        [86] = (259, "d 1/2,0,0; c 0,1/2,0; b 1/2,1/2,0; a 0,0,0"),   // P 1 1 21/b
        [90] = (262, "e 0,y,1/4; d 1/4,1/4,1/2; c 1/4,1/4,0; b 0,1/2,0; a 0,0,0"),     // C 1 2/c 1
        [91] = (263, "e 3/4,y,3/4; d 1/2,1/4,3/4; c 0,1/4,1/4; b 0,1/2,0; a 0,0,0"),   // A 1 2/n 1
        [92] = (263, "e 1/4,y,0; d 1/4,1/4,3/4; c 3/4,1/4,3/4; b 0,1/2,0; a 0,0,0"),   // I 1 2/a 1
        [96] = (266, "e 1/4,0,z; d 1/2,1/4,1/4; c 0,1/4,1/4; b 0,0,1/2; a 0,0,0"),     // A 1 1 2/a
        [97] = (267, "e 3/4,3/4,z; d 3/4,1/2,1/4; c 1/4,0,1/4; b 0,0,1/2; a 0,0,0"),   // B 1 1 2/n
        [98] = (267, "e 0,1/4,z; d 3/4,1/4,1/4; c 3/4,3/4,1/4; b 0,0,1/2; a 0,0,0"),   // I 1 1 2/b
    };

    /// <summary>基底の変換 P: 新しい基底 (a′,b′,c′) を元の基底 (a,b,c) で表したものを列に持つ (下では列を並べて転置する)。x′ = P⁻¹ x</summary>
    private static readonly int[][] BtoC = Transpose([[0, 0, 1], [1, 0, 0], [0, 1, 0]]);  // (a′,b′,c′) = (c, a, b)  ITA Table 1.5.1.1 b 軸→c 軸
    private static readonly int[][] BtoA = Transpose([[0, 1, 0], [0, 0, 1], [1, 0, 0]]);  // (a′,b′,c′) = (b, c, a)  ITA Table 1.5.1.1 b 軸→a 軸
    private static readonly int[][] CmBA = Transpose([[0, 0, 1], [0, -1, 0], [1, 0, 0]]);  // 設定 c-ba: (c, -b, a)    ITA Table 1.5.4.4
    private static readonly int[][] BAmC = Transpose([[0, 1, 0], [1, 0, 0], [0, 0, -1]]);  // 設定 ba-c: (b, a, -c)
    private static readonly int[][] mACB = Transpose([[-1, 0, 0], [0, 0, 1], [0, 1, 0]]);  // 設定 -acb: (-a, c, b)

    /// <summary>印刷されていない設定: (写す元の通し番号, 基底の列)。a 軸の設定は b 軸の同じセル選択から、c-ba・ba-c・-acb は No. 15 だけ
    /// (ほかの群では同じ記号・同じ群の設定と同一なので spglib の通し番号を持たない)</summary>
    private static readonly Dictionary<int, (int Source, int[][] Basis)> Derived = new()
    {
        [5] = (3, BtoA), [20] = (18, BtoA), [59] = (57, BtoA), [62] = (60, BtoA),
        [69] = (63, BtoA), [70] = (64, BtoA), [71] = (65, BtoA),
        [78] = (72, BtoA), [79] = (73, BtoA), [80] = (74, BtoA),
        [87] = (81, BtoA), [88] = (82, BtoA), [89] = (83, BtoA),
        [102] = (90, BtoA), [103] = (91, BtoA), [104] = (92, BtoA),
        [93] = (90, CmBA), [94] = (91, CmBA), [95] = (92, CmBA),
        [99] = (96, BAmC), [100] = (97, BAmC), [101] = (98, BAmC),
        [105] = (102, mACB), [106] = (103, mACB), [107] = (104, mACB),
    };

    private static readonly Dictionary<char, double> Free = new() { ['x'] = 0.1371, ['y'] = 0.2293, ['z'] = 0.3187 };

    /// <summary>"1/2,y,0" → 成分ごとに (固定値, 自由な変数 or null)</summary>
    private static (double Value, char? Var)[] Parse(string rep)
        => [.. rep.Split(',').Select(c => c.Trim() is var t && t.Length == 1 && Free.ContainsKey(t[0]) ? (0.0, (char?)t[0])
            : (t.Contains('/') ? double.Parse(t.Split('/')[0]) / double.Parse(t.Split('/')[1]) : double.Parse(t), (char?)null))];

    private static int[][] Transpose(int[][] m) => [.. Enumerable.Range(0, 3).Select(i => Enumerable.Range(0, 3).Select(j => m[j][i]).ToArray())];

    /// <summary>符号付き置換の基底の列 P で座標を写す (x′ = P⁻¹ x = Pᵀ x)。自由な成分は新しい軸の変数名にする</summary>
    private static (double, char?)[] Transform((double Value, char? Var)[] rep, int[][] basis)
    {
        var q = Transpose(basis);
        return [.. Enumerable.Range(0, 3).Select(i =>
        {
            var j = Array.FindIndex(q[i], v => v != 0);
            return rep[j].Var is null ? ((q[i][j] * rep[j].Value % 1 + 1) % 1, (char?)null) : (0.0, (char?)"xyz"[i]);
        })];
    }

    /// <summary>通し番号 series の、ITA の文字と代表座標</summary>
    private static Dictionary<string, (double Value, char? Var)[]> ItaLetters(int series)
    {
        if (Printed.TryGetValue(series, out var p))
            return p.Positions.Split(';').Select(s => s.Trim().Split(' ', 2)).ToDictionary(s => s[0], s => Parse(s[1]));
        var (source, basis) = Derived[series];
        return ItaLetters(source).ToDictionary(kv => kv.Key, kv => Transform(kv.Value, basis));
    }

    /// <summary>基準の対称操作の集合を、基底の列 P で写したキーの集合 (R′ = Pᵀ R P, t′ = Pᵀ t)</summary>
    private static HashSet<string> TransformedOpKeys(int series, int[][] basis)
    {
        var q = Transpose(basis);
        int[] mul(int[][] a, int[] r9, int[][] b)
        {
            var r = new int[9];
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    for (int k = 0; k < 3; k++)
                        for (int l = 0; l < 3; l++)
                            r[i * 3 + l] += a[i][j] * r9[j * 3 + k] * b[k][l];
            return r;
        }
        return [.. Groups[series].Ops.Select(o =>
        {
            var t = Enumerable.Range(0, 3).Select(i => q[i][0] * o[9] + q[i][1] * o[10] + q[i][2] * o[11]).ToArray();
            return Key(mul(q, o[..9], basis), t[0], t[1], t[2]);
        })];
    }

    public static IEnumerable<object[]> MonoclinicSeries() => Printed.Keys.Concat(Derived.Keys).Order().Select(s => new object[] { s });

    [Theory, MemberData(nameof(MonoclinicSeries))]
    public void WyckoffLetters_MatchIta(int series)
    {
        var wrong = new List<string>();
        foreach (var (letter, rep) in ItaLetters(series))
        {
            var p = rep.Select(c => c.Var is char v ? Free[v] : c.Value).ToArray();
            var a = WyckoffPosition.GetEquivalentAtomsPosition((p[0], p[1], p[2]), series);
            if (a.WyckoffLeter != letter)
                wrong.Add($"{letter} ({string.Join(",", p.Select(v => v.ToString("0.###")))}) → '{a.WyckoffLeter}'");
        }
        Assert.True(wrong.Count == 0, $"series {series} ({SymmetryStatic.Symmetries[series].SpaceGroupHMfullStr}): {string.Join(", ", wrong)}");
    }

    /// <summary>写した設定の対称操作が、写す先の設定の対称操作と一致する (変換の取り方の確認)</summary>
    [Theory, MemberData(nameof(DerivedSeries))]
    public void DerivedSettings_TransformOntoTargetGroup(int series)
    {
        var (source, basis) = Derived[series];
        Assert.True(TransformedOpKeys(source, basis).SetEquals(Groups[series].OpKeys), $"series {series} ← {source}");
    }

    public static IEnumerable<object[]> DerivedSeries() => Derived.Keys.Order().Select(s => new object[] { s });

    /// <summary>書き写しの確認: 印刷された c 軸の設定は、印刷された b 軸の設定 (同じセル選択) を b 軸→c 軸の変換で写したものと同じ点の集合</summary>
    [Theory]
    [InlineData(3, 4)] [InlineData(18, 19)] [InlineData(57, 58)] [InlineData(60, 61)]
    [InlineData(63, 66)] [InlineData(64, 67)] [InlineData(65, 68)]
    [InlineData(72, 75)] [InlineData(73, 76)] [InlineData(74, 77)]
    [InlineData(81, 84)] [InlineData(82, 85)] [InlineData(83, 86)]
    [InlineData(90, 96)] [InlineData(91, 97)] [InlineData(92, 98)]
    public void PrintedCAxisSettings_AreBAxisSettingsTransformed(int bSeries, int cSeries)
    {
        Assert.True(TransformedOpKeys(bSeries, BtoC).SetEquals(Groups[cSeries].OpKeys));
        // ReciPro の表を使わずに、基準の対称操作だけで比べる: 写した代表点の像のどれかが、印刷された同じ文字の座標 (固定された成分) に乗る
        var printedC = ItaLetters(cSeries);
        foreach (var (letter, rep) in ItaLetters(bSeries))
        {
            var p = Transform(rep, BtoC).Select(c => c.Item2 is char v ? Free[v] : c.Item1).ToArray();
            var target = printedC[letter];
            bool onLocus(double[] q) => Enumerable.Range(0, 3).All(i => target[i].Var is not null || Math.Abs(Math.IEEERemainder(q[i] - target[i].Value, 1)) < 1e-9);
            var hit = Groups[cSeries].Ops.Any(o => onLocus([.. Enumerable.Range(0, 3).Select(i => o[3 * i] * p[0] + o[3 * i + 1] * p[1] + o[3 * i + 2] * p[2] + (double)o[9 + i] / Den)]));
            Assert.True(hit, $"{letter}: series {bSeries} → {cSeries}");
        }
    }
}
