// 260712Cl 新規: 32 の幾何結晶類 (結晶学的点群型) の包含 poset — FormGroupRelations の点群 Hasse 図タブの
// データ層。設計は codex R12 で確定 (.project-guidance/ReciPro_FormGroupRelations改修計画.md):
//   - 型 B ≤ A ⟺ 「A 型のある代表点群が B 型に共役な部分群を含む」(型レベル包含は代表の取り方に依存しない)。
//   - 32 ノード全体は束 (lattice) ではなく poset — 互いに比較不能な最大元が 2 つある (m-3m と 6/mmm)。
//   - 被覆辺 (Hasse 辺) は 80 本。「具体的極大部分群を型へ畳んだ辺」と「全包含関係の transitive reduction」が
//     32 型では一致する (codex R12 が独立に列挙・確認。本クラスは両方を計算して不一致なら throw する自己検証付き)。
//   - 計算は手書きテーブルなしの実行時自己構築: 代表の線形部 → 乗積表 → 全部分群列挙 (TSubgroupFinder の
//     既存機構) → 回転型シグネチャで 32 型へ分類 → 極大化 / transitive reduction。
using System;
using System.Collections.Generic;
using System.Linq;

namespace Crystallography;

/// <summary>幾何結晶類 (点群型) 1 つ分の情報。260712Cl 追加。</summary>
public sealed class PointGroupTypeInfo
{
    /// <summary>正規化 HM 記号 (2mm/m2m → mm2)。</summary>
    public string Name { get; init; }
    /// <summary>Schoenflies 記号。</summary>
    public string Schoenflies { get; init; }
    /// <summary>位数 (点群の元の数)。</summary>
    public int Order { get; init; }
    /// <summary>この点群を持つ空間群型の数 (230 型中、SpaceGroupNumber の distinct)。</summary>
    public int SpaceGroupTypeCount { get; init; }
}

/// <summary>32 点群型の包含 poset (Hasse 図データ)。初回アクセス時に全 530 設定から自己構築。260712Cl 追加 (codex R12)。</summary>
public static class PointGroupCatalog
{
    /// <summary>32 型 (位数降順 → 名前昇順の決定的順序)。</summary>
    public static IReadOnlyList<PointGroupTypeInfo> Types => _data.Value.Types;

    /// <summary>Hasse 図の被覆辺 (80 本)。Index = 親位数 / 子位数。</summary>
    public static IReadOnlyList<(string Parent, string Child, int Index)> CoverEdges => _data.Value.Edges;

    /// <summary>型 → 極大部分型のリスト (被覆辺の親側索引)。</summary>
    public static IReadOnlyDictionary<string, string[]> MaximalSubtypes => _data.Value.MaxSub;

    /// <summary>全 32 型の極大部分群の親内共役類数の合計 (= 95、検証用)。</summary>
    public static int MaximalConjugacyClassTotal => _data.Value.MaxClassTotal;

    /// <summary>設定 seriesNumber の点群型名 (正規化 HM)。</summary>
    public static string NormalizedName(int seriesNumber)
        => Normalize(SymmetryStatic.Symmetries[seriesNumber].PointGroupHMStr);

    private static string Normalize(string hm) => hm switch { "2mm" or "m2m" => "mm2", var t => t };

    private sealed class Data
    {
        public PointGroupTypeInfo[] Types;
        public (string Parent, string Child, int Index)[] Edges;
        public Dictionary<string, string[]> MaxSub;
        public int MaxClassTotal;
    }

    private static readonly Lazy<Data> _data = new(Compute, System.Threading.LazyThreadSafetyMode.ExecutionAndPublication);

    private static Data Compute()
    {
        // ---- 代表設定の選定と空間群型数の集計 ----
        var repSn = new Dictionary<string, int>();       // 型名 → 代表 series
        var sgTypes = new Dictionary<string, HashSet<int>>(); // 型名 → 空間群番号集合
        for (int sn = 1; sn < SymmetryStatic.TotalSpaceGroupNumber; sn++)
        {
            var sym = SymmetryStatic.Symmetries[sn];
            if (sym.SpaceGroupNumber == 0) continue;
            var name = Normalize(sym.PointGroupHMStr);
            if (string.IsNullOrEmpty(name)) continue;
            if (!repSn.ContainsKey(name)) repSn[name] = sn;
            if (!sgTypes.TryGetValue(name, out var set)) sgTypes[name] = set = [];
            set.Add(sym.SpaceGroupNumber);
        }
        if (repSn.Count != 32)
            throw new InvalidOperationException($"expected 32 point-group types, got {repSn.Count}");

        // ---- 各型: 部分群列挙 → 型分類 → 包含集合・極大部分群 ----
        var order = new Dictionary<string, int>();
        var containsTypes = new Dictionary<string, HashSet<string>>(); // 真部分群として現れる型
        var maxSub = new Dictionary<string, string[]>();
        int maxClassTotal = 0;
        foreach (var (name, sn) in repSn)
        {
            // 線形部の抽出 (conventional 整数行列、重複除去)
            var linKeys = new List<int[]>();
            foreach (var op in TSubgroupFinder.GetExpandedOps(sn))
            {
                var key = TSubgroupFinder.LinKey(op);
                if (TSubgroupFinder.FindKey(linKeys, key) < 0) linKeys.Add(key);
            }
            int n = linKeys.Count;
            order[name] = n;
            var mul = new int[n, n];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                {
                    mul[i, j] = TSubgroupFinder.FindKey(linKeys, KSubgroupFinder.MatMulInt(linKeys[i], linKeys[j]));
                    if (mul[i, j] < 0) throw new InvalidOperationException("point group not closed under multiplication");
                }
            //int e = Enumerable.Range(0, n).First(i => linKeys[i][0] == 1 && linKeys[i][4] == 1 && linKeys[i][8] == 1
            //    && linKeys[i][1] == 0 && linKeys[i][2] == 0 && linKeys[i][3] == 0 && linKeys[i][5] == 0 && linKeys[i][6] == 0 && linKeys[i][7] == 0);
            int e = Enumerable.Range(0, n).First(i => TSubgroupFinder.IsIdentity(linKeys[i])); // 260717Cl: インライン判定を既存 IsIdentity へ集約

            var subs = TSubgroupFinder.EnumerateSubgroups(n, mul, e);
            string TypeOf(SortedSet<int> s) => TSubgroupFinder.SignatureNameMap.Value[TSubgroupFinder.Signature(s.Select(i => linKeys[i]))];

            var proper = subs.Where(s => s.Count < n).ToList();
            containsTypes[name] = [.. proper.Select(TypeOf)];

            // 極大部分群 = 他の真部分群に真に含まれない真部分群
            var maximal = proper.Where(s => !proper.Any(t => t.Count > s.Count && s.IsSubsetOf(t))).ToList();
            maxSub[name] = [.. maximal.Select(TypeOf).Distinct().OrderBy(t => t)];
            maxClassTotal += TSubgroupFinder.GroupByConjugacy(maximal, n, mul, linKeys).Count;
        }

        // ---- transitive reduction (全包含) と極大型辺の一致検証 (codex R12: 32 型では一致するはず) ----
        var edges = new List<(string Parent, string Child, int Index)>();
        foreach (var (parent, kids) in containsTypes)
        {
            foreach (var child in kids)
            {
                // 被覆 ⟺ parent ⊃ child で、中間 c (parent ⊃ c ⊃ child) が無い
                bool covered = kids.Any(c => c != child && containsTypes[c].Contains(child));
                if (!covered)
                    edges.Add((parent, child, order[parent] / order[child]));
            }
        }
        foreach (var (parent, kids) in maxSub)
        {
            var reduced = edges.Where(ed => ed.Parent == parent).Select(ed => ed.Child).OrderBy(t => t).ToArray();
            if (!reduced.SequenceEqual(kids))
                throw new InvalidOperationException($"maximal-subtype edges disagree with transitive reduction for {parent}");
        }

        var types = repSn.Keys
            .Select(nm => new PointGroupTypeInfo
            {
                Name = nm,
                Schoenflies = SymmetryStatic.Symmetries[repSn[nm]].PointGroupSFStr,
                Order = order[nm],
                SpaceGroupTypeCount = sgTypes[nm].Count,
            })
            .OrderByDescending(t => t.Order).ThenBy(t => t.Name, StringComparer.Ordinal)
            .ToArray();

        return new Data
        {
            Types = types,
            Edges = [.. edges.OrderBy(ed => ed.Parent, StringComparer.Ordinal).ThenBy(ed => ed.Child, StringComparer.Ordinal)],
            MaxSub = maxSub,
            MaxClassTotal = maxClassTotal,
        };
    }
}
