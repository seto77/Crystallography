using System.Text.Json;
using System.Text.Json.Serialization;

namespace Crystallography.Tests;

/// <summary>261105Cl 新規: 基準データ (Reference/symmetry_reference.json, spglib / cctbx 由来) の読み込みと、
/// ReciPro 側の対称操作を比較用のキーに直す補助関数。
/// 対称操作のキーは "r11,…,r33;t1,t2,t3" (t は 1/24 単位の整数, mod 24)。x' = R x + t。</summary>
public static class SymmetryReference
{
    public const int Den = 24;

    public sealed class Group
    {
        [JsonPropertyName("series")] public int Series { get; set; }
        [JsonPropertyName("hall_symbol")] public string HallSymbol { get; set; }
        [JsonPropertyName("ops")] public int[][] Ops { get; set; }
        [JsonPropertyName("cctbx_reference_setting")] public bool CctbxReferenceSetting { get; set; }
        [JsonPropertyName("wyckoff")] public Wyckoff[] Wyckoff { get; set; }
        [JsonPropertyName("number")] public int Number { get; set; }
        [JsonPropertyName("hm_short")] public string HmShort { get; set; }
        [JsonPropertyName("hm_full")] public string HmFull { get; set; }
        [JsonPropertyName("schoenflies")] public string Schoenflies { get; set; }
        [JsonPropertyName("choice")] public string Choice { get; set; }
        [JsonPropertyName("pg_hm")] public string PgHm { get; set; }
        [JsonPropertyName("pg_schoenflies")] public string PgSchoenflies { get; set; }

        /// <summary>対称操作キーの集合 (センタリング込み)</summary>
        public HashSet<string> OpKeys => _opKeys ??= [.. Ops.Select(o => Key(o[..9], o[9], o[10], o[11]))];
        private HashSet<string> _opKeys;

        /// <summary>回転部分 (重複なし)。等価な面・軸の基準に使う</summary>
        public int[][] Rotations => _rot ??= [.. Ops.Select(o => o[..9]).DistinctBy(r => string.Join(",", r))];
        private int[][] _rot;
    }

    public sealed class Wyckoff
    {
        [JsonPropertyName("letter")] public string Letter { get; set; }
        [JsonPropertyName("mult")] public int Multiplicity { get; set; }
        [JsonPropertyName("site")] public string Site { get; set; }
        /// <summary>Wyckoff 位置への射影演算子 [r11..r33, t1, t2, t3] (cctbx の special_op)</summary>
        [JsonPropertyName("special")] public double[] Special { get; set; }
    }

    private sealed class Root
    {
        [JsonPropertyName("groups")] public Group[] Groups { get; set; }
    }

    public static readonly IReadOnlyDictionary<int, Group> Groups = Load();

    private static Dictionary<int, Group> Load()
    {
        var path = Path.Combine(AppContext.BaseDirectory, "Reference", "symmetry_reference.json");
        var root = JsonSerializer.Deserialize<Root>(File.ReadAllText(path));
        return root.Groups.ToDictionary(g => g.Series);
    }

    /// <summary>ITA の 530 設定 (ReciPro 通し番号 1..530 = spglib の Hall 番号)</summary>
    public static IEnumerable<object[]> ItaSeries() => Enumerable.Range(1, 530).Select(i => new object[] { i });

    /// <summary>ReciPro 独自の設定 (531..540: センタリング付き三斜晶) を含む全設定</summary>
    public static IEnumerable<object[]> AllSeries() => Enumerable.Range(1, 540).Select(i => new object[] { i });

    public static string Key(IReadOnlyList<int> r, int t1, int t2, int t3)
        => $"{string.Join(",", r)};{Mod(t1)},{Mod(t2)},{Mod(t3)}";

    private static int Mod(int t) => ((t % Den) + Den) % Den;

    /// <summary>実数の (R, t) をキーに直す。R は整数、t は 1/24 の倍数でなければ例外。</summary>
    public static string Key(double[] r9, double t1, double t2, double t3)
    {
        var r = r9.Select(v =>
        {
            var n = Math.Round(v);
            if (Math.Abs(v - n) > 1e-9) throw new InvalidOperationException($"non-integer rotation element {v}");
            return (int)n;
        }).ToArray();
        return Key(r, T24(t1), T24(t2), T24(t3));
    }

    private static int T24(double t)
    {
        var n = Math.Round(t * Den);
        if (Math.Abs(t * Den - n) > 1e-6) throw new InvalidOperationException($"translation {t} is not a multiple of 1/{Den}");
        return (int)n;
    }

    /// <summary>アフィン変換 f(x) = M x + t を、原点と基底ベクトルの像から (M, t) に分解する</summary>
    public static (double[] M, double[] T) Decompose(Func<double, double, double, (double X, double Y, double Z)> f)
    {
        var o = f(0, 0, 0); var a = f(1, 0, 0); var b = f(0, 1, 0); var c = f(0, 0, 1);
        return ([a.X - o.X, b.X - o.X, c.X - o.X, a.Y - o.Y, b.Y - o.Y, c.Y - o.Y, a.Z - o.Z, b.Z - o.Z, c.Z - o.Z], [o.X, o.Y, o.Z]);
    }

    /// <summary>ReciPro の一般位置 (Wyckoff[0] の生成器, センタリング込み) のキー列</summary>
    public static List<string> GeneralPositionKeys(int series)
        => [.. SymmetryStatic.WyckoffPositions[series][0].PositionGenerator.Select(g =>
        {
            var (m, t) = Decompose(g.Apply);
            return Key(m, t[0], t[1], t[2]);
        })];

    /// <summary>ReciPro の Seitz 形式の対称操作 (Wyckoff[0].PositionOperations) のキー列。
    /// 本番 (GeneratePositions / GetExpandedOps) と同じく通し番号で包み直してから行列を取り出す。</summary>
    public static List<string> SeitzKeys(int series)
        => [.. SymmetryStatic.WyckoffPositions[series][0].PositionOperations.Select(so0 =>
        {
            var so = new SymmetryOperation(so0, series);
            var (m, _) = Decompose(so.ApplyMatrix);
            var t = so.SeitzTranslation;
            return Key(m, t.U, t.V, t.W);
        })];
}
