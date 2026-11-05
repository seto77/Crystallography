using Xunit;
using static Crystallography.Tests.SymmetryReference;

namespace Crystallography.Tests;

/// <summary>261105Cl 新規: 空間群データ (SymmetryStatic / Symmetry / WyckoffPosition) を spglib・cctbx 由来の基準と照合する。
/// 通し番号 1..530 は ITA の 530 設定 (spglib の Hall 番号と一致)、531..540 は ReciPro 独自のセンタリング付き三斜晶。</summary>
public class SpaceGroupTests
{
    #region 対称操作
    /// <summary>一般位置の座標式 (センタリング込み) が基準の対称操作の集合と一致する</summary>
    [Theory, MemberData(nameof(AllSeries), MemberType = typeof(SymmetryReference))]
    public void GeneralPositions_MatchReference(int series)
    {
        var keys = GeneralPositionKeys(series);
        Assert.Equal(keys.Count, keys.Distinct().Count());
        Assert.True(Groups[series].OpKeys.SetEquals(keys),
            $"series {series}: 基準にない操作 [{string.Join(" ", keys.Except(Groups[series].OpKeys).Take(4))}] / 欠けている操作 [{string.Join(" ", Groups[series].OpKeys.Except(keys).Take(4))}]");
    }

    /// <summary>Seitz 形式の対称操作 (PositionOperations) が一般位置の座標式と添字ごとに一致する。
    /// GeneratePositions は両者を同じ添字で対応づけて原子に操作を持たせ、異方性 DW の計算がこれに依存する。</summary>
    [Theory, MemberData(nameof(AllSeries), MemberType = typeof(SymmetryReference))]
    public void SeitzOperations_MatchGeneralPositionsOneToOne(int series)
        => Assert.Equal(GeneralPositionKeys(series), SeitzKeys(series));
    #endregion

    #region 消滅則
    /// <summary>手書きの消滅則 (Symmetry.HasExtinction) が、対称操作から導いた系統的消滅と一致する (|h|,|k|,|l| ≤ 10)。
    /// 反射 h は、h·R = h かつ h·t ∉ Z となる操作 (R, t) があるとき系統的に消滅する。</summary>
    [Theory, MemberData(nameof(AllSeries), MemberType = typeof(SymmetryReference))]
    public void ExtinctionRules_MatchOperations(int series)
    {
        const int n = 10;
        var sym = SymmetryStatic.Symmetries[series];
        var ops = Groups[series].Ops;
        var wrong = new List<string>();
        for (int h = -n; h <= n; h++)
            for (int k = -n; k <= n; k++)
                for (int l = -n; l <= n; l++)
                {
                    if (h == 0 && k == 0 && l == 0) continue;
                    bool absent = false;
                    foreach (var o in ops)
                    {
                        // h·R (行ベクトル × R)
                        int h2 = h * o[0] + k * o[3] + l * o[6], k2 = h * o[1] + k * o[4] + l * o[7], l2 = h * o[2] + k * o[5] + l * o[8];
                        if (h2 == h && k2 == k && l2 == l && (h * o[9] + k * o[10] + l * o[11]) % Den != 0) { absent = true; break; }
                    }
                    if (absent != sym.HasExtinction(h, k, l) && wrong.Count < 10)
                        wrong.Add($"({h} {k} {l}) 期待={(absent ? "消滅" : "出現")}");
                }
        Assert.True(wrong.Count == 0, $"series {series} ({sym.SpaceGroupHMStr}): {string.Join(", ", wrong)}");
    }
    #endregion

    #region 等価な面・軸
    private static readonly (int, int, int)[] TestVectors =
    [
        (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 0, 1), (0, 1, 1), (1, 1, 1), (1, -1, 0), (2, 1, 0), (1, 2, 0),
        (2, -1, 0), (1, -2, 0), (1, 2, 3), (3, 1, 2), (2, 3, 1), (1, -2, 3), (-1, 2, 3), (1, 1, 2), (2, 1, 1), (1, 2, 1),
        (1, 1, -2), (2, -1, -1), (3, 0, 1), (1, 0, 3), (0, 3, 1), (4, 2, 1), (1, 1, 3), (2, 2, 1), (0, 1, 2), (2, 0, 1),
    ];

    /// <summary>GenerateEquivalentPlanes が点群の回転部分による軌道 {h·R} と一致する (対称心の付加あり・なし)</summary>
    [Theory, MemberData(nameof(AllSeries), MemberType = typeof(SymmetryReference))]
    public void EquivalentPlanes_MatchPointGroupOrbit(int series)
    {
        var sym = SymmetryStatic.Symmetries[series];
        var rots = Groups[series].Rotations;
        var wrong = new List<string>();
        foreach (var (h, k, l) in TestVectors)
        {
            var exp = rots.Select(r => (h * r[0] + k * r[3] + l * r[6], h * r[1] + k * r[4] + l * r[7], h * r[2] + k * r[5] + l * r[8])).ToHashSet();
            var expInv = exp.Concat(exp.Select(v => (-v.Item1, -v.Item2, -v.Item3))).ToHashSet();
            foreach (var (inv, e) in new[] { (false, exp), (true, expInv) })
            {
                var got = SymmetryStatic.GenerateEquivalentPlanes((h, k, l), sym, inv).Select(v => (v.H, v.K, v.L)).ToList();
                if (got.Count != got.Distinct().Count() || !e.SetEquals(got))
                    wrong.Add($"({h} {k} {l}) inv={inv}: {got.Count} 個 (期待 {e.Count} 個)");
            }
        }
        Assert.True(wrong.Count == 0, $"series {series} ({sym.SpaceGroupHMStr}): {string.Join(", ", wrong.Take(6))}");
    }

    /// <summary>GenerateEquivalentAxes が点群の回転部分による軌道 {R·u} と一致する (対称心の付加あり・なし)</summary>
    [Theory, MemberData(nameof(AllSeries), MemberType = typeof(SymmetryReference))]
    public void EquivalentAxes_MatchPointGroupOrbit(int series)
    {
        var sym = SymmetryStatic.Symmetries[series];
        var rots = Groups[series].Rotations;
        var wrong = new List<string>();
        foreach (var (u, v, w) in TestVectors)
        {
            var exp = rots.Select(r => (r[0] * u + r[1] * v + r[2] * w, r[3] * u + r[4] * v + r[5] * w, r[6] * u + r[7] * v + r[8] * w)).ToHashSet();
            var expInv = exp.Concat(exp.Select(x => (-x.Item1, -x.Item2, -x.Item3))).ToHashSet();
            foreach (var (inv, e) in new[] { (false, exp), (true, expInv) })
            {
                var got = SymmetryStatic.GenerateEquivalentAxes((u, v, w), sym, inv).Select(x => (x.U, x.V, x.W)).ToList();
                if (got.Count != got.Distinct().Count() || !e.SetEquals(got))
                    wrong.Add($"[{u} {v} {w}] inv={inv}: {got.Count} 個 (期待 {e.Count} 個)");
            }
        }
        Assert.True(wrong.Count == 0, $"series {series} ({sym.SpaceGroupHMStr}): {string.Join(", ", wrong.Take(6))}");
    }
    #endregion

    #region 記号・番号
    private static string NormHmShort(string s)
    {
        s = s.Split('=')[0];
        foreach (var suf in new[] { "(1)", "(2)", "Hex", "Rho" }) s = s.Replace(suf, "");
        return s.Replace("sub", "_").Replace(" ", "");
    }

    private static string NormHmFull(string s) => s.Replace("Hex", "").Replace("Rho", "").Replace("sub", "_").Replace(" ", "");

    private static int CrystalSystemOf(int number)
        => number <= 2 ? 1 : number <= 15 ? 2 : number <= 74 ? 3 : number <= 142 ? 4 : number <= 167 ? 5 : number <= 194 ? 6 : 7;

    private static readonly string[] CrystalSystemNames = ["", "triclinic", "monoclinic", "orthorhombic", "tetragonal", "trigonal", "hexagonal", "cubic"];

    /// <summary>Hall 記号・H-M 記号・Schoenflies 記号・点群・晶系が基準 (spglib) と一致する</summary>
    [Theory, MemberData(nameof(ItaSeries), MemberType = typeof(SymmetryReference))]
    public void Symbols_MatchSpglib(int series)
    {
        var sym = SymmetryStatic.Symmetries[series];
        var r = Groups[series];
        Assert.Equal(r.Number, sym.SpaceGroupNumber);
        Assert.Equal(r.HallSymbol, sym.SpaceGroupHallStandardStr);
        // spglib 2.x は Hall 331 (No. 68, 2bca) の短縮記号だけを旧記号 'Bbcb' で返す (同じ番号の 330/332/333 は 'Bbeb')。
        var hmShort = series == 331 ? "Bbeb" : r.HmShort;
        Assert.Equal(hmShort.Replace(" ", ""), NormHmShort(sym.SpaceGroupHMStr));
        Assert.Equal(r.HmFull.Replace(" ", ""), NormHmFull(sym.SpaceGroupHMfullStr));
        Assert.Equal(r.Schoenflies, sym.SpaceGroupSFStr);
        Assert.Equal(r.PgSchoenflies, sym.PointGroupSFStr);
        Assert.Equal(CrystalSystemOf(r.Number), sym.CrystalSystemNumber);
        Assert.Equal(CrystalSystemNames[CrystalSystemOf(r.Number)], sym.CrystalSystemStr);
    }

    /// <summary>点群番号 (PointGroupNumber) と点群記号 (PointGroupSFStr) が全設定で 1 対 1 に対応する。
    /// GenerateEquivalentPlanes/Axes と EBSD のマスターパターン対称化は点群番号で分岐する。</summary>
    [Fact]
    public void PointGroupNumber_IsConsistentWithPointGroupSymbol()
    {
        var bad = Enumerable.Range(1, 540).Select(i => SymmetryStatic.Symmetries[i])
            .GroupBy(s => s.PointGroupNumber)
            .Where(g => g.Select(s => s.PointGroupSFStr).Distinct().Count() > 1)
            .Select(g => $"点群番号 {g.Key}: " + string.Join(", ", g.GroupBy(s => s.PointGroupSFStr).Select(x => $"{x.Key}×{x.Count()} (例 {x.First().SeriesNumber} {x.First().SpaceGroupHMStr})")))
            .ToList();
        Assert.True(bad.Count == 0, string.Join(" / ", bad));
    }

    /// <summary>mm2 型の向き付き点群記号 (2mm / m2m / mm2) が、対称操作の 2 回軸の向きと一致する</summary>
    [Theory, MemberData(nameof(ItaSeries), MemberType = typeof(SymmetryReference))]
    public void OrientedPointGroupSymbol_Mm2_MatchesTwofoldAxis(int series)
    {
        var sym = SymmetryStatic.Symmetries[series];
        if (Groups[series].PgSchoenflies != "C2v") return;
        // 回転部分の中から真性 2 回回転 (対角成分の 1 つが +1, 他が -1) を探す
        var axis = Groups[series].Rotations
            .Where(r => r[0] + r[4] + r[8] == -1 && r[0] * (r[4] * r[8] - r[5] * r[7]) - r[1] * (r[3] * r[8] - r[5] * r[6]) + r[2] * (r[3] * r[7] - r[4] * r[6]) == 1)
            .Select(r => r[0] == 1 ? 0 : r[4] == 1 ? 1 : 2).Single();
        Assert.Equal(new[] { "2mm", "m2m", "mm2" }[axis], sym.PointGroupHMStr);
    }

    /// <summary>同じ標準 Hall 記号を持つ設定は同じ群でなければならない (CIF の Hall 記号から設定を引くため)。
    /// 従来の空白なし表記 (SpaceGroupHallStr) は 'P32' (P3_2 と P312) などが重なるので、CIF には標準表記を使う。</summary>
    [Fact]
    public void HallSymbols_IdentifyGroupsUniquely()
    {
        static string Norm(string s) => string.Join(" ", s.Split(' ', StringSplitOptions.RemoveEmptyEntries));
        var bad = Enumerable.Range(1, 540)
            .GroupBy(i => Norm(SymmetryStatic.Symmetries[i].SpaceGroupHallStandardStr))
            .Where(g => g.Count() > 1 && g.Select(i => string.Join("|", Groups[i].OpKeys.Order())).Distinct().Count() > 1)
            .Select(g => $"'{g.Key}': " + string.Join(", ", g.Select(i => $"{i} {SymmetryStatic.Symmetries[i].SpaceGroupHMStr}")))
            .ToList();
        Assert.True(bad.Count == 0, string.Join(" / ", bad));
    }

    /// <summary>261105Cl 追加: 算術結晶類 (SymmetryProperties.ArithmeticCrystalClassStr) は設定に依らず、IT 番号 1..230 で 73 類。
    /// 向きで分かれる 312P/321P などを区別する (以前は点群 + 格子の記号で 91 種になり、312P と 321P が同じ「32P」だった)。</summary>
    [Fact]
    public void ArithmeticCrystalClass_IsSettingIndependentAnd73()
    {
        var byIt = Enumerable.Range(1, 540).GroupBy(i => (int)SymmetryStatic.Symmetries[i].SpaceGroupNumber)
            .ToDictionary(g => g.Key, g => g.Select(i => new SymmetryProperties(SymmetryStatic.Symmetries[i]).ArithmeticCrystalClassStr).Distinct().ToArray());
        var multi = byIt.Where(kv => kv.Value.Length != 1).Select(kv => $"{kv.Key}: {string.Join("/", kv.Value)}").ToList();
        Assert.True(multi.Count == 0, string.Join(", ", multi));
        Assert.Equal(73, byIt.Values.Select(v => v[0]).Distinct().Count());
        Assert.Equal(new[] { "312P", "321P", "3m1P", "31mP", "-31mP", "-3m1P", "-6m2P", "-62mP", "-42mP", "-4m2P", "-4m2I", "-42mI" },
            new[] { 149, 150, 156, 157, 162, 164, 187, 189, 111, 115, 119, 121 }.Select(it => byIt[it][0]));
    }
    #endregion

    #region Wyckoff 位置
    private static string NormSite(string s)
    {
        s = s.Replace(".", "");
        return s switch { "2mm" or "m2m" => "mm2", "-4m2" => "-42m", "-6m2" => "-62m", _ => s };
    }

    /// <summary>Wyckoff 位置の表が基準 (cctbx) と一致する: 位置の数、各位置の多重度とサイト対称の型。
    /// さらに本番の判定 (WyckoffPosition.GetEquivalentAtomsPosition) に、各 Wyckoff 位置上の点を
    /// 対称操作と格子並進で動かしてから渡し、多重度・サイト対称・文字が正しく返るかを確かめる。
    /// 文字は cctbx の参照設定 (基底変換が恒等) のときだけ比べる。非標準設定の文字は基底変換の取り方に依存し、
    /// cctbx と spglib でも一致しない場合があるため。cctbx の '@' (ITA の α) は ReciPro では 'A'。</summary>
    /// <summary>cctbx の Wyckoff 文字が ITA と異なる設定。522 (Pn-3m 原点選択 2) は 24i (1/2,y,y+1/2) と 24j (1/2,y,-y) が逆
    /// (ReciPro・spglib・pyxtal (Bilbao WYCKPOS 由来) は ITA どおり)。これらの設定では文字を比べない。</summary>
    private static readonly HashSet<int> CctbxLetterErrata = [522];

    [Theory, MemberData(nameof(AllSeries), MemberType = typeof(SymmetryReference))]
    public void WyckoffPositions_MatchCctbx(int series)
    {
        var r = Groups[series];
        bool checkLetter = r.CctbxReferenceSetting && !CctbxLetterErrata.Contains(series);
        var table = SymmetryStatic.WyckoffPositions[series];
        Assert.True(r.Wyckoff.Length == table.Length, $"series {series}: Wyckoff 位置の数 {table.Length} (期待 {r.Wyckoff.Length})");
        Assert.Equal(r.Wyckoff.Select(w => (w.Multiplicity, w.Site)).Order(), table.Select(w => ((int)w.Multiplicity, NormSite(w.SiteSymmetry))).Order());

        double[] free = [0.1234567, 0.2718282, 0.3819660]; // 有理数から十分離れた値 (偶然に特殊位置へ乗らないように)
        var ops = r.Ops;
        var wrong = new List<string>();
        foreach (var w in r.Wyckoff)
        {
            var s = w.Special;
            double px = s[0] * free[0] + s[1] * free[1] + s[2] * free[2] + s[9];
            double py = s[3] * free[0] + s[4] * free[1] + s[5] * free[2] + s[10];
            double pz = s[6] * free[0] + s[7] * free[1] + s[8] * free[2] + s[11];
            foreach (var o in new[] { ops[0], ops[ops.Length / 2], ops[^1] })
            {
                double x = o[0] * px + o[1] * py + o[2] * pz + (double)o[9] / Den + 1;
                double y = o[3] * px + o[4] * py + o[5] * pz + (double)o[10] / Den - 1;
                double z = o[6] * px + o[7] * py + o[8] * pz + (double)o[11] / Den + 1;
                var a = WyckoffPosition.GetEquivalentAtomsPosition((x, y, z), series);
                var expLetter = w.Letter == "@" ? "A" : w.Letter;
                if (a.Multiplicity != w.Multiplicity || NormSite(a.SiteSymmetry) != w.Site || (checkLetter && a.WyckoffLeter != expLetter))
                {
                    wrong.Add($"{w.Letter}({w.Multiplicity},{w.Site}) → '{a.WyckoffLeter}'({a.Multiplicity},{a.SiteSymmetry})");
                    break;
                }
            }
        }
        Assert.True(wrong.Count == 0, $"series {series} ({SymmetryStatic.Symmetries[series].SpaceGroupHMStr}): {string.Join(", ", wrong.Take(6))}");
    }
    #endregion
}
