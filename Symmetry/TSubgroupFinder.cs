// 260704Cl 追加: translationengleiche (t-) 最大部分群の実行時厳密計算エンジン。FormGroupRelations (Phase 2) のデータ源。
//
// 原理: t-部分群は並進格子を保つため、空間群 G の t-部分群は点群 P_G の部分群と 1:1 対応する。
//   1. G の対称操作の線形部 (整数行列) から点群 P_G を取り出し、全部分群を閉包列挙 → 極大部分群 → P_G 共役で類別。
//   2. 各類 M について H = { g∈G : 線形部(g) ∈ M } (中心化並進込み) を構成。coset 代表 = ドメイン/双晶操作。
//   3. 型同定: 候補基底変換 P (カタログ) × 候補原点シフト p を試し、H を子基準系へ写して
//      既存 530 設定の操作集合と「格子法で完全一致」するものを探す。一致 = 同定 (誤同定は原理的に不可能)。
//      一致が見つからない場合は「未同定」として点群名のみ提示する (正直なフォールバック)。
//   4. 軌道分裂 (親 Wyckoff → H 軌道への分割) と、消滅則差分 (親で消滅・子で許容になる反射) も操作集合から直接計算。
//
// この計算は将来のオフライン生成 CSV (GAP/spglib パイプライン) の独立差分検証器を兼ねる
// (ReciPro_SymmetryInformation拡張計画.md §5.1「自作 enumerator」)。k-部分群・超群はここでは扱わない (Phase 2 データ待ち)。
//
// 規約: 基底 (a′,b′,c′) = (a,b,c)·P、座標 x_parent = P·x_child + p、演算子は x′ = R·x + t (列ベクトル左作用)。
// 子基準系での操作は R_c = P⁻¹RP、t_c = P⁻¹((R−I)p + t)。det(P) > 0 に限定 (エナンチオモルフ対の混同防止)。
using System;
using System.Collections.Generic;
using System.Linq;

namespace Crystallography;

/// <summary>群-部分群関係の種別。260705Cl 追加 (Phase 2e)。k/Isomorphic は将来の KSubgroupFinder 用の予約 (未使用)。</summary>
public enum GroupRelationKind { T, K, Isomorphic }

/// <summary>群-部分群関係の 1 共役類を表す共通 DTO (260705Cl 追加, Phase 2e。旧 TSubgroup/TSubgroupFinder.TSupergroup を統合)。
/// <see cref="TSubgroupFinder.GetMaximalTSubgroups"/> / <see cref="TSubgroupFinder.GetMinimalTSupergroups"/> が返す。
/// UI (FormGroupRelations) はこの DTO のみを読み、TSubgroupFinder / 将来の KSubgroupFinder を直接知らない。
/// Parent は常に群として大きい側 (index の分母側)、Child は小さい側。P,p は常に「子基準系 → 親基準系」の向き
/// (x_parent = P·x_child + p) で格納する。Minimal supergroups 一覧に載る関係は、逆引き元である Parent 自身の
/// 部分群表から取った値をそのまま使う (Parent 側から見た Child への変換)。現在の閲覧対象が Parent でなく Child 側
/// (supergroup を見ている) のときは <see cref="GetInverseTransform"/> で (P,p)⁻¹ を求めて表示する。</summary>
public sealed class GroupRelation
{
    /// <summary>関係の種別 (現在は T のみ実データ)。</summary>
    public GroupRelationKind Kind { get; init; } = GroupRelationKind.T;
    /// <summary>親空間群の通し番号。</summary>
    public int ParentSeriesNumber { get; init; }
    /// <summary>指数 [G:H] = |P_G| / |M|。t-部分群では並進指数 1 なので点群指数に等しい。</summary>
    public int Index { get; init; }
    /// <summary>親の中での共役類 ID (0 始まり、Compute() 呼び出し内で安定)。260705Cl 追加 (Phase 2e、将来の k- 共役類区別に使う予約)。</summary>
    public int ConjugacyClassId { get; init; }
    /// <summary>この共役類に属する部分群の個数 (方位バリアント数)。</summary>
    public int ConjugateCount { get; init; }
    /// <summary>H の点群の HM 記号 (正規化済み: 2mm/m2m→mm2 等)。</summary>
    public string PointGroupHM { get; init; }
    /// <summary>H の全対称操作 (親基準系、中心化展開済み)。</summary>
    public SymmetryOperation[] Operations { get; init; }
    /// <summary>H の coset 代表 (線形部が相異なる代表操作、恒等を含む)。</summary>
    public SymmetryOperation[] Representatives { get; init; }
    /// <summary>G の H に対する coset 代表 (恒等 coset を除く)。ドメインを結ぶ操作 = 双晶則。</summary>
    public SymmetryOperation[] CosetRepresentatives { get; init; }
    /// <summary>同定された標準設定の通し番号。未同定なら -1。</summary>
    public int ChildSeriesNumber { get; init; } = -1;
    /// <summary>基底変換 P (row-major 9 要素、(a′,b′,c′)=(a,b,c)·P)。未同定なら null。</summary>
    public double[] TransformP { get; init; }
    /// <summary>原点シフト p (親座標系)。未同定なら null。</summary>
    public double[] TransformShift { get; init; }
    /// <summary>部分格子基底 T′ (k- のみ、t- では null)。260705Cl 追加 (Phase 2c KSubgroupFinder 用の予約)。</summary>
    public double[] SublatticeBasis { get; init; }

    /// <summary>同定済みなら子設定の HM 記号 (sub 表記含む生文字列)、未同定なら点群 HM。</summary>
    public string ChildLabel => ChildSeriesNumber >= 0 ? SymmetryStatic.Symmetries[ChildSeriesNumber].SpaceGroupHMStr : PointGroupHM;

    /// <summary>(P,p) の逆変換 (P⁻¹, −P⁻¹·p) を返す。Minimal supergroups 側 (Child から Parent を見る向き) の
    /// Matrix タブ表示用 (260705Cl 追加, Phase 2e)。未同定 (TransformP=null) なら (null, null)。</summary>
    public (double[] P, double[] Shift) GetInverseTransform()
    {
        if (TransformP == null) return (null, null);
        var pinv = TSubgroupFinder.Invert3(TransformP);
        if (pinv == null) return (null, null); // 特異行列は理論上発生しない (P は基底変換で正則) が防御的に
        var shift = new double[]
        {
            -(pinv[0] * TransformShift[0] + pinv[1] * TransformShift[1] + pinv[2] * TransformShift[2]),
            -(pinv[3] * TransformShift[0] + pinv[4] * TransformShift[1] + pinv[5] * TransformShift[2]),
            -(pinv[6] * TransformShift[0] + pinv[7] * TransformShift[1] + pinv[8] * TransformShift[2]),
        };
        return (pinv, shift);
    }
}

/// <summary>親 Wyckoff 軌道が t-部分群でどう分裂するかの 1 成分 (260704Cl 追加)。</summary>
public readonly struct OrbitPart
{
    /// <summary>親の慣用胞内での軌道点数。</summary>
    public int CountInParentCell { get; init; }
    /// <summary>子設定での Wyckoff 文字 (同定済みのときのみ、それ以外 null)。</summary>
    public string ChildWyckoffLetter { get; init; }
    /// <summary>子設定での多重度 (子の慣用胞基準)。未同定なら 0。</summary>
    public int ChildMultiplicity { get; init; }
    /// <summary>子設定でのサイト対称性。未同定なら null。</summary>
    public string ChildSiteSymmetry { get; init; }
}

/// <summary>t-最大部分群の列挙・型同定・軌道分裂・消滅則差分 (260704Cl 追加)。結果は空間群ごとにキャッシュされる。</summary>
public static class TSubgroupFinder
{
    private const double Tol = 1e-6;
    private static readonly Dictionary<int, GroupRelation[]> _cache = [];
    private static readonly object _lock = new();

    #region 公開 API
    /// <summary>親空間群 (通し番号) の maximal t-部分群を共役類単位で返す。計算は初回のみ (キャッシュ)。</summary>
    public static GroupRelation[] GetMaximalTSubgroups(int seriesNumber)
    {
        lock (_lock)
        {
            if (_cache.TryGetValue(seriesNumber, out var cached))
                return cached;
            var result = Compute(seriesNumber);
            _cache[seriesNumber] = result;
            return result;
        }
    }

    // 260705Cl: 専用の軽量 TSupergroup 型 (SupergroupSeriesNumber/Index/ConjugateCount のみ) を廃し、
    // GroupRelation をそのまま逆引き索引に格納する (Phase 2e DTO 統合)。ParentSeriesNumber が supergroup の
    // 通し番号、ChildSeriesNumber が引数 itNumber 側の設定。Operations/TransformP 等の全データが引き続き手に入るため、
    // Minimal supergroups 側でも Matrix/Orbit/Reflections タブが (P,p)⁻¹ 経由で正しく表示できる。

    private static Dictionary<int, List<GroupRelation>> _supergroupIndex;

    /// <summary>逆引き索引が構築済みか。260705Cl 追加: 初回構築は全 230 タイプの部分群計算 (数秒) を伴うため、
    /// GUI 側はこれを見て「未構築ならバックグラウンドで構築 → 完了時に表示を差し替える」を選べる。</summary>
    public static bool SupergroupIndexReady { get { lock (_lock) return _supergroupIndex != null; } }

    /// <summary>指定空間群タイプ (IT 番号) を maximal t-部分群に持つ空間群 (= minimal t-supergroup) を返す。
    /// 全 230 タイプの第 1 設定を 1 度だけ走査して逆引き索引を構築する (translationengleiche のみ)。</summary>
    public static IReadOnlyList<GroupRelation> GetMinimalTSupergroups(int itNumber)
    {
        lock (_lock)
        {
            if (_supergroupIndex == null)
            {
                _supergroupIndex = [];
                for (int it = 1; it <= 230; it++)
                {
                    //int sn = FirstSeriesOf(it);
                    int sn = SymmetryStatic.GetSeriesNumber(it, 1); // 260705Cl: 既存 API を再利用 (全 IT 番号で第 1 設定は sub=1)
                    if (sn < 0) continue;
                    foreach (var sub in GetMaximalTSubgroups(sn))
                    {
                        if (sub.ChildSeriesNumber < 0) continue;
                        int childIt = SymmetryStatic.Symmetries[sub.ChildSeriesNumber].SpaceGroupNumber;
                        if (!_supergroupIndex.TryGetValue(childIt, out var list))
                            _supergroupIndex[childIt] = list = [];
                        list.Add(sub); // 260705Cl: sub (GroupRelation) をそのまま格納 (旧: 3 フィールドだけの TSupergroup を新規生成)
                    }
                }
            }
            return _supergroupIndex.TryGetValue(itNumber, out var result) ? result : [];
        }
    }

    // 260705Cl: 既存 SymmetryStatic.GetSeriesNumber(number, sub:1) の近似重複だったため削除。
    //private static int FirstSeriesOf(int itNumber)
    //{
    //    for (int sn = 0; sn < SymmetryStatic.TotalSpaceGroupNumber; sn++)
    //        if (SymmetryStatic.Symmetries[sn].SpaceGroupNumber == itNumber)
    //            return sn;
    //    return -1;
    //}

    /// <summary>操作集合ベースの系統的消滅判定: ∃(R,t): h·R = h かつ h·t ∉ Z。</summary>
    public static bool IsExtinct(IReadOnlyList<SymmetryOperation> ops, int h, int k, int l)
    {
        foreach (var op in ops)
        {
            var (h2, k2, l2) = op.ConvertPlaneIndex(h, k, l);
            if (h2 == h && k2 == k && l2 == l)
            {
                var t = op.SeitzTranslation;
                double phase = h * t.U + k * t.V + l * t.W;
                if (Math.Abs(phase - Math.Round(phase)) > 1e-4)
                    return true;
            }
        }
        return false;
    }

    /// <summary>親の各 Wyckoff 位置 (index 順) の H による軌道分裂を返す。generic 代表点によるサンプル計算。</summary>
    public static OrbitPart[][] GetOrbitSplitting(int parentSeries, GroupRelation sub)
    {
        if (sub.Kind == GroupRelationKind.K) return GetOrbitSplittingK(parentSeries, sub); // 260708Cl (Phase 2d)
        // 特殊関係 (x=y, 2x=z 等) を偶然踏まない generic 値
        const double gx = 0.127743, gy = 0.291317, gz = 0.437129;
        var wycks = SymmetryStatic.WyckoffPositions[parentSeries];
        var parentOps = GetExpandedOps(parentSeries);
        var result = new OrbitPart[wycks.Length][];
        var inv = sub.ChildSeriesNumber >= 0 ? Invert3(sub.TransformP) : null; // 260705Cl: 不変値を軌道ループ外へ

        for (int w = 0; w < wycks.Length; w++)
        {
            var (rx, ry, rz) = wycks[w].PositionGenerator[0].Apply(gx, gy, gz);
            var orbit = GenerateOrbit(parentOps, rx, ry, rz);
            var parts = new List<OrbitPart>();
            var used = new bool[orbit.Count];
            for (int i = 0; i < orbit.Count; i++)
            {
                if (used[i]) continue;
                var subOrbit = GenerateOrbit(sub.Operations, orbit[i].X, orbit[i].Y, orbit[i].Z);
                foreach (var p in subOrbit)
                    for (int j = 0; j < orbit.Count; j++)
                        if (!used[j] && Near(orbit[j], p))
                            used[j] = true;

                string letter = null, siteSym = null;
                int childMult = 0;
                if (sub.ChildSeriesNumber >= 0)
                {
                    // x_child = P⁻¹(x_parent − p)
                    //var inv = Invert3(sub.TransformP); // 260705Cl: ループ外へ hoist
                    double px = orbit[i].X - sub.TransformShift[0], py = orbit[i].Y - sub.TransformShift[1], pz = orbit[i].Z - sub.TransformShift[2];
                    var cx = inv[0] * px + inv[1] * py + inv[2] * pz;
                    var cy = inv[3] * px + inv[4] * py + inv[5] * pz;
                    var cz = inv[6] * px + inv[7] * py + inv[8] * pz;
                    var atoms = WyckoffPosition.GetEquivalentAtomsPosition((cx, cy, cz), sub.ChildSeriesNumber);
                    letter = atoms.WyckoffLeter;
                    childMult = atoms.Multiplicity;
                    siteSym = atoms.SiteSymmetry;
                }
                parts.Add(new OrbitPart { CountInParentCell = subOrbit.Count, ChildWyckoffLetter = letter, ChildMultiplicity = childMult, ChildSiteSymmetry = siteSym });
            }
            result[w] = [.. parts];
        }
        return result;
    }

    /// <summary>親で系統的消滅・部分群 H で許容になる反射 (超構造反射の t-版) を列挙する。
    /// 子の等価反射 (Friedel 込み) で代表 1 つに集約し、(代表 hkl, 等価数, 親の消滅則) を返す。</summary>
    public static (int H, int K, int L, int EquivCount, string ParentRule)[] GetNewReflections(int parentSeries, GroupRelation sub, int maxIndex = 4)
    {
        if (sub.Kind == GroupRelationKind.K) return GetNewReflectionsK(parentSeries, sub, maxIndex); // 260708Cl (Phase 2d)
        var parentOps = GetExpandedOps(parentSeries);
        var parentSym = SymmetryStatic.Symmetries[parentSeries];
        var seen = new HashSet<(int, int, int)>();
        var list = new List<(int, int, int, int, string)>();

        for (int h = maxIndex; h >= -maxIndex; h--)
            for (int k = maxIndex; k >= -maxIndex; k--)
                for (int l = maxIndex; l >= -maxIndex; l--)
                {
                    if ((h == 0 && k == 0 && l == 0) || seen.Contains((h, k, l)))
                        continue;
                    if (!IsExtinct(parentOps, h, k, l) || IsExtinct(sub.Operations, h, k, l))
                        continue;

                    // 子の等価反射 + Friedel 対で軌道を張り、代表 (辞書順最大) のみ採用
                    var orbit = new HashSet<(int, int, int)> { (h, k, l), (-h, -k, -l) };
                    bool grown = true;
                    while (grown)
                    {
                        grown = false;
                        foreach (var q in orbit.ToArray())
                            foreach (var rep in sub.Representatives)
                            {
                                var r = rep.ConvertPlaneIndex(q.Item1, q.Item2, q.Item3);
                                if (orbit.Add(r) | orbit.Add((-r.H, -r.K, -r.L)))
                                    grown = true;
                            }
                    }
                    foreach (var q in orbit)
                        seen.Add(q);
                    var canon = orbit.OrderByDescending(q => q).First();
                    list.Add((canon.Item1, canon.Item2, canon.Item3, orbit.Count, parentSym.GetFirstExtinctionRule(h, k, l) ?? ""));
                }

        list.Sort((a, b) => (a.Item1 * a.Item1 + a.Item2 * a.Item2 + a.Item3 * a.Item3).CompareTo(b.Item1 * b.Item1 + b.Item2 * b.Item2 + b.Item3 * b.Item3));
        return [.. list];
    }

    /// <summary>260708Cl 追加 (Phase 2d): k-部分群 (klassengleiche) の軌道分裂。t- と異なり並進格子が粗くなる
    /// (子胞が大きくなる) ため、親の G-軌道を子座標系へ写して子慣用胞で mod1 集約し、子の全操作 (中心化込み) で
    /// H-軌道へ分割する。これにより中心化変化 (IIa/IIb) を含め子の慣用胞多重度が正しく数えられる。
    /// 子が未同定 (k ではほぼ発生しない) の場合は各 Wyckoff を空成分で返す。</summary>
    private static OrbitPart[][] GetOrbitSplittingK(int parentSeries, GroupRelation sub)
    {
        const double gx = 0.127743, gy = 0.291317, gz = 0.437129; // 特殊関係を偶然踏まない generic 値
        var wycks = SymmetryStatic.WyckoffPositions[parentSeries];
        var result = new OrbitPart[wycks.Length][];
        if (sub.ChildSeriesNumber < 0 || sub.TransformP == null)
        {
            for (int w = 0; w < wycks.Length; w++) result[w] = [];
            return result;
        }
        var parentOps = GetExpandedOps(parentSeries);
        var childOps = GetExpandedOps(sub.ChildSeriesNumber);
        var P = sub.TransformP;
        var inv = Invert3(P); // x_child = P⁻¹(x_parent − p)
        var p = sub.TransformShift;

        // 子胞は親胞の |det P| 倍の体積なので、親操作を generic 点へ作用させただけでは 1 親胞分しか埋まらない。
        // 子胞を満たすには親格子並進 n のコピーが要る。子座標での並進コセット Frac(P⁻¹·n) を列挙する
        // (親格子/子格子の剰余類群、位数 |det P|)。子胞が親より小さい (中心化喪失で det P<1) 場合は親格子点が
        // 子格子に含まれるため Frac(P⁻¹·n) が全て 0 となり自動的に {0} に縮退する。
        double detP = P[0] * (P[4] * P[8] - P[5] * P[7]) - P[1] * (P[3] * P[8] - P[5] * P[6]) + P[2] * (P[3] * P[7] - P[4] * P[6]);
        int C = Math.Max(2, (int)Math.Ceiling(Math.Abs(detP)) + 1); // コセット群を確実に列挙する n の範囲 [0,C)
        var fills = new List<(double X, double Y, double Z)>();
        for (int a = 0; a < C; a++)
            for (int b = 0; b < C; b++)
                for (int c = 0; c < C; c++)
                {
                    var f = (X: Frac(inv[0] * a + inv[1] * b + inv[2] * c), Y: Frac(inv[3] * a + inv[4] * b + inv[5] * c), Z: Frac(inv[6] * a + inv[7] * b + inv[8] * c));
                    if (!fills.Any(u => Near(u, f))) fills.Add(f);
                }

        for (int w = 0; w < wycks.Length; w++)
        {
            var (rx, ry, rz) = wycks[w].PositionGenerator[0].Apply(gx, gy, gz);
            // 親 G-軌道を子座標系へ写し、格子並進コピーで子胞を満たして mod1 集約 (子胞内の全 G-軌道点)。
            // dedup は許容誤差ベースの Near で行う (座標が 1e-4 グリッド境界にちょうど乗ると整数キー方式は
            // 同一点でも演算経路差で丸めが割れて別点扱いになり、分割が過剰計上する不具合があったため)。
            var big = new List<(double X, double Y, double Z)>();
            foreach (var op in parentOps)
            {
                var (mx, my, mz) = op.ApplyMatrix(rx, ry, rz);
                var t = op.SeitzTranslation;
                double qx = mx + t.U - p[0], qy = my + t.V - p[1], qz = mz + t.W - p[2];
                double bx = inv[0] * qx + inv[1] * qy + inv[2] * qz;
                double by = inv[3] * qx + inv[4] * qy + inv[5] * qz;
                double bz = inv[6] * qx + inv[7] * qy + inv[8] * qz;
                foreach (var f in fills)
                {
                    var v = (X: Frac(bx + f.X), Y: Frac(by + f.Y), Z: Frac(bz + f.Z));
                    if (!big.Any(u => Near(u, v))) big.Add(v);
                }
            }
            // 子の全操作で H-軌道へ分割 (mod1 子胞)。子軌道は G-軌道の部分集合なので必ず big に含まれる。
            var parts = new List<OrbitPart>();
            var used = new bool[big.Count];
            for (int i = 0; i < big.Count; i++)
            {
                if (used[i]) continue;
                var subOrbit = GenerateOrbit(childOps, big[i].X, big[i].Y, big[i].Z);
                foreach (var q in subOrbit)
                    for (int j = 0; j < big.Count; j++)
                        if (!used[j] && Near(big[j], q))
                            used[j] = true;
                var atoms = WyckoffPosition.GetEquivalentAtomsPosition((big[i].X, big[i].Y, big[i].Z), sub.ChildSeriesNumber);
                parts.Add(new OrbitPart { CountInParentCell = subOrbit.Count, ChildWyckoffLetter = atoms.WyckoffLeter, ChildMultiplicity = atoms.Multiplicity, ChildSiteSymmetry = atoms.SiteSymmetry });
            }
            result[w] = [.. parts];
        }
        return result;
    }

    /// <summary>260708Cl 追加 (Phase 2d): k-部分群で新たに現れる反射 (超格子反射) を子の指数で列挙する。
    /// 子指数 (h',k',l') を親指数 (h,k,l)=(h',k',l')·P⁻¹ に写して 2 分類する:
    /// ① fractional-index (超格子): 親で非整数指数 → 胞拡大で新規出現。ParentRule に親分数指数を "(…)" で格納。
    /// ② released: 親で整数だが系統消滅・子で許容 → 消滅則解除。ParentRule に親の消滅則を格納。
    /// いずれも子で許容 (子の消滅則で消えない) 反射のみ。子の対称等価 + Friedel 対で代表 1 つに集約する。
    /// 子が未同定なら空 (子の消滅則が判定できないため)。</summary>
    private static (int H, int K, int L, int EquivCount, string ParentRule)[] GetNewReflectionsK(int parentSeries, GroupRelation sub, int maxIndex)
    {
        if (sub.ChildSeriesNumber < 0 || sub.TransformP == null)
            return [];
        var parentOps = GetExpandedOps(parentSeries);
        var parentSym = SymmetryStatic.Symmetries[parentSeries];
        var childOps = GetExpandedOps(sub.ChildSeriesNumber);
        var inv = Invert3(sub.TransformP); // (h,k,l) = (h',k',l')·P⁻¹
        var seen = new HashSet<(int, int, int)>();
        var list = new List<(int, int, int, int, string)>();

        for (int h = maxIndex; h >= -maxIndex; h--)
            for (int k = maxIndex; k >= -maxIndex; k--)
                for (int l = maxIndex; l >= -maxIndex; l--)
                {
                    if ((h == 0 && k == 0 && l == 0) || seen.Contains((h, k, l)))
                        continue;
                    if (IsExtinct(childOps, h, k, l)) // 子で消える反射は観測されない
                        continue;
                    // 子指数 → 親指数 (行ベクトル×P⁻¹)。分類 (整数=消滅則解除 / 非整数=超格子) は子の対称軌道
                    // (= 親点群軌道、k- では点群不変) で不変なので、まず (h,k,l) で基本反射 (retained) を除外する。
                    double ph = h * inv[0] + k * inv[3] + l * inv[6];
                    double pk = h * inv[1] + k * inv[4] + l * inv[7];
                    double pl = h * inv[2] + k * inv[5] + l * inv[8];
                    if (IsInt(ph) && IsInt(pk) && IsInt(pl) && !IsExtinct(parentOps, (int)Math.Round(ph), (int)Math.Round(pk), (int)Math.Round(pl)))
                        continue; // 親でも許容 = 基本反射 (新規でない)

                    // 子の対称等価 + Friedel 対で軌道を張り、代表 (辞書順最大) のみ採用。
                    var orbit = new HashSet<(int, int, int)> { (h, k, l), (-h, -k, -l) };
                    bool grown = true;
                    while (grown)
                    {
                        grown = false;
                        foreach (var q in orbit.ToArray())
                            foreach (var op in childOps)
                            {
                                var r = op.ConvertPlaneIndex(q.Item1, q.Item2, q.Item3);
                                if (orbit.Add(r) | orbit.Add((-r.H, -r.K, -r.L)))
                                    grown = true;
                            }
                    }
                    foreach (var q in orbit)
                        seen.Add(q);
                    var canon = orbit.OrderByDescending(q => q).First();
                    // 表示する代表 canon の親指数で分類ラベルを作る (表示子指数と親指数の対応を整合させる)。
                    double ch = canon.Item1 * inv[0] + canon.Item2 * inv[3] + canon.Item3 * inv[6];
                    double ck = canon.Item1 * inv[1] + canon.Item2 * inv[4] + canon.Item3 * inv[7];
                    double cl = canon.Item1 * inv[2] + canon.Item2 * inv[5] + canon.Item3 * inv[8];
                    string rule = IsInt(ch) && IsInt(ck) && IsInt(cl)
                        ? parentSym.GetFirstExtinctionRule((int)Math.Round(ch), (int)Math.Round(ck), (int)Math.Round(cl)) ?? "" // 消滅則解除
                        : $"({FracStr(ch)} {FracStr(ck)} {FracStr(cl)})"; // 親分数指数 = 超格子反射
                    list.Add((canon.Item1, canon.Item2, canon.Item3, orbit.Count, rule));
                }

        list.Sort((a, b) => (a.Item1 * a.Item1 + a.Item2 * a.Item2 + a.Item3 * a.Item3).CompareTo(b.Item1 * b.Item1 + b.Item2 * b.Item2 + b.Item3 * b.Item3));
        return [.. list];
    }

    private static bool IsInt(double d) => Math.Abs(d - Math.Round(d)) < 1e-4; // 260708Cl 追加
    /// <summary>260708Cl 追加: 親分数指数の 1 成分を短い分数/整数文字列へ (超格子反射の親指数表示用)。</summary>
    private static string FracStr(double d)
    {
        d = Math.Round(d, 6);
        if (Math.Abs(d - Math.Round(d)) < 1e-4) return ((int)Math.Round(d)).ToString(System.Globalization.CultureInfo.InvariantCulture);
        foreach (int den in new[] { 2, 3, 4, 6 })
        {
            double x = d * den;
            if (Math.Abs(x - Math.Round(x)) < 1e-4)
                return $"{(int)Math.Round(x)}/{den}";
        }
        return d.ToString("0.##", System.Globalization.CultureInfo.InvariantCulture);
    }
    #endregion

    #region 本体計算
    private static GroupRelation[] Compute(int sn)
    {
        var ops = GetExpandedOps(sn);
        if (ops.Length == 0)
            return [];

        // --- 1. 線形部 (整数行列) の抽出と点群の乗積表 ---
        var linKeys = new List<int[]>();               // 相異なる線形部 (int[9], row-major)
        var opsByLin = new List<List<SymmetryOperation>>(); // 線形部ごとの操作 (中心化コピー含む)
        foreach (var op in ops)
        {
            var key = LinKey(op);
            int idx = FindKey(linKeys, key);
            if (idx < 0) { linKeys.Add(key); opsByLin.Add([]); idx = linKeys.Count - 1; }
            opsByLin[idx].Add(op);
        }
        int n = linKeys.Count; // |P_G|
        var mul = new int[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                mul[i, j] = FindKey(linKeys, MatMul(linKeys[i], linKeys[j]));
        int e = Enumerable.Range(0, n).First(i => IsIdentity(linKeys[i]));

        // --- 2. 全部分群の閉包列挙 → 極大 → 共役類 ---
        var subgroups = EnumerateSubgroups(n, mul, e);
        var maximal = subgroups.Where(s => s.Count < n && !subgroups.Any(t2 => t2.Count < n && t2.Count > s.Count && s.IsSubsetOf(t2))).ToList();
        var classes = GroupByConjugacy(maximal, n, mul, linKeys);

        // --- 3. 各共役類の代表について H を構成・型同定 ---
        var sigName = SignatureNameMap.Value;
        var result = new List<GroupRelation>();
        for (int classIdx = 0; classIdx < classes.Count; classIdx++)
        {
            var cls = classes[classIdx];
            // 同定は類内のどの共役でも等価 (型は共通)。カタログの向きに合う共役があれば拾えるよう全共役を試す。
            int child = -1;
            double[] bestP = null, bestShift = null;
            SortedSet<int> reprM = cls[0];
            foreach (var m in cls)
            {
                var hOps = m.SelectMany(i => opsByLin[i]).ToArray();
                (child, bestP, bestShift) = Identify(sn, hOps, m.Select(i => linKeys[i]).ToList());
                if (child >= 0) { reprM = m; break; }
            }

            var mSet = reprM;
            var subOps = mSet.SelectMany(i => opsByLin[i]).ToArray();
            var reps = mSet.Select(i => opsByLin[i][0]).ToArray();

            // G の coset 代表 (恒等 coset 以外): 双晶操作
            var visited = new HashSet<int>(mSet);
            var cosetReps = new List<SymmetryOperation>();
            for (int g = 0; g < n; g++)
            {
                if (visited.Contains(g)) continue;
                foreach (var m2 in mSet)
                    visited.Add(mul[m2, g]);
                cosetReps.Add(opsByLin[g][0]);
            }

            result.Add(new GroupRelation
            {
                ParentSeriesNumber = sn,
                Index = n / mSet.Count,
                ConjugacyClassId = classIdx,
                ConjugateCount = cls.Count,
                PointGroupHM = sigName.TryGetValue(Signature(mSet.Select(i => linKeys[i])), out var nm) ? nm : "?",
                Operations = subOps,
                Representatives = reps,
                CosetRepresentatives = [.. cosetReps],
                ChildSeriesNumber = child,
                TransformP = bestP,
                TransformShift = bestShift,
            });
        }
        // index 昇順 → 点群位数降順で安定表示
        return [.. result.OrderBy(r => r.Index).ThenBy(r => r.ChildSeriesNumber < 0 ? 1 : 0).ThenBy(r => r.ChildSeriesNumber)];
    }

    /// <summary>一般位置の全対称操作 (中心化展開済み・seriesNumber 付替え済み) を返す。
    /// 260705Cl: SymmetryProperties / FormSymmetryInformation と共用するため public 化 (4 箇所の同型展開を一本化)。</summary>
    //private static SymmetryOperation[] GetExpandedOps(int sn)
    public static SymmetryOperation[] GetExpandedOps(int sn)
    {
        if (sn <= 0 || sn >= SymmetryStatic.TotalSpaceGroupNumber)
            return [];
        var raw = SymmetryStatic.WyckoffPositions[sn][0].PositionOperations;
        if (raw == null || raw.Length == 0)
            return [];
        return [.. raw.Select(o => new SymmetryOperation(o, sn))];
    }
    #endregion

    #region 部分群列挙・共役類
    /// <summary>乗積表 mul 上で全部分群 (要素 index の集合) を閉包列挙する。</summary>
    private static List<SortedSet<int>> EnumerateSubgroups(int n, int[,] mul, int e)
    {
        var found = new List<SortedSet<int>>();
        var keys = new HashSet<string>();
        void Add(SortedSet<int> s)
        {
            var key = string.Join(",", s);
            if (keys.Add(key)) found.Add(s);
        }
        Add([e]);
        // 巡回部分群から開始
        for (int g = 0; g < n; g++)
            Add(Closure([e, g], mul));
        // 既知の部分群 + 1 要素で拡張し飽和するまで
        for (int pass = 0; pass < found.Count; pass++)
        {
            var s = found[pass];
            if (s.Count == n) continue;
            for (int g = 0; g < n; g++)
            {
                if (s.Contains(g)) continue;
                var ext = new SortedSet<int>(s) { g };
                Add(Closure(ext, mul));
            }
        }
        return found;
    }

    private static SortedSet<int> Closure(IEnumerable<int> seed, int[,] mul)
    {
        var s = new SortedSet<int>(seed);
        var queue = new Queue<int>(s);
        while (queue.Count > 0)
        {
            int a = queue.Dequeue();
            foreach (var b in s.ToArray())
            {
                foreach (var c in new[] { mul[a, b], mul[b, a] })
                    if (s.Add(c)) queue.Enqueue(c);
            }
        }
        return s;
    }

    /// <summary>部分群集合を P_G 共役 (gMg⁻¹) で類別する。</summary>
    private static List<List<SortedSet<int>>> GroupByConjugacy(List<SortedSet<int>> subs, int n, int[,] mul, List<int[]> linKeys)
    {
        // 逆元表
        int e = Enumerable.Range(0, n).First(i => IsIdentity(linKeys[i]));
        var inv = new int[n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                if (mul[i, j] == e) { inv[i] = j; break; }

        var classes = new List<List<SortedSet<int>>>();
        var assigned = new bool[subs.Count];
        for (int i = 0; i < subs.Count; i++)
        {
            if (assigned[i]) continue;
            var cls = new List<SortedSet<int>> { subs[i] };
            assigned[i] = true;
            for (int g = 0; g < n; g++)
            {
                var conj = new SortedSet<int>(subs[i].Select(m => mul[mul[g, m], inv[g]]));
                for (int j = 0; j < subs.Count; j++)
                    if (!assigned[j] && conj.SetEquals(subs[j]))
                    {
                        cls.Add(subs[j]);
                        assigned[j] = true;
                    }
            }
            classes.Add(cls);
        }
        return classes;
    }
    #endregion

    #region 点群の命名 (回転型シグネチャ)
    /// <summary>回転型 (±1,±2,±3,±4,±6) の多重集合 → 点群 HM 名。全 530 設定から自己構築 (手書きテーブル不使用)。</summary>
    private static readonly Lazy<Dictionary<string, string>> SignatureNameMap = new(() =>
    {
        var map = new Dictionary<string, string>();
        for (int s = 1; s < SymmetryStatic.TotalSpaceGroupNumber; s++)
        {
            // 260705Cl: 線形部抽出〜シグネチャ生成の重複実装を SettingSignature に一本化 (キャッシュも共用)。
            // 点群名も生配列添字 (StrArray[s][13]) から既存プロパティ参照へ。
            //var ops = GetExpandedOps(s);
            //if (ops.Length == 0) continue;
            //var lin = new List<int[]>();
            //foreach (var op in ops)
            //{
            //    var key = LinKey(op);
            //    if (FindKey(lin, key) < 0) lin.Add(key);
            //}
            //var sig = Signature(lin);
            //var name = SymmetryStatic.StrArray[s][13] switch { "2mm" or "m2m" => "mm2", var t2 => t2 };
            var sig = SettingSignature(s);
            if (sig == "") continue;
            var name = SymmetryStatic.Symmetries[s].PointGroupHMStr switch { "2mm" or "m2m" => "mm2", var t2 => t2 };
            if (map.TryGetValue(sig, out var prev))
            {
                if (prev != name)
                    throw new InvalidOperationException($"point-group signature collision: {prev} vs {name}");
            }
            else
                map[sig] = name;
        }
        return map;
    });

    private static string Signature(IEnumerable<int[]> lins)
        => string.Join(",", lins.Select(RotationType).OrderBy(v => v));

    /// <summary>整数行列の回転型: det=+1 → trace 3,-1,0,1,2 = 1,2,3,4,6 / det=−1 → -1,-2(m),-3,-4,-6。trace は基底不変。</summary>
    private static int RotationType(int[] m)
    {
        int det = Det3(m);
        int tr = m[0] + m[4] + m[8];
        return det > 0
            ? tr switch { 3 => 1, -1 => 2, 0 => 3, 1 => 4, 2 => 6, _ => throw new InvalidOperationException("bad rotation") }
            : tr switch { -3 => -1, 1 => -2, 0 => -3, -1 => -4, -2 => -6, _ => throw new InvalidOperationException("bad rotoinversion") };
    }
    #endregion

    #region 型同定 (候補 (P,p) の全数検証)
    // 候補基底変換カタログ (列 = 子基底ベクトルの親座標成分)。24 の proper 符号付置換と右から合成する。
    private static readonly double[][] CellChanges =
    [
        [1, 0, 0, 0, 1, 0, 0, 0, 1],                                     // 恒等
        [0.5, 0.5, 0, -0.5, 0.5, 0, 0, 0, 1],                            // a′=(a−b)/2, b′=(a+b)/2 (F 立方 → I 正方 等)
        [1, 1, 0, -1, 1, 0, 0, 0, 1],                                    // a′=a−b, b′=a+b (P 正方 → C 直方 等)
        [1, 1, 0, 0, 2, 0, 0, 0, 1],                                     // a′=a, b′=a+2b (六方 → 直方 C, orthohexagonal)
        [0, 0.5, 0.5, 0.5, 0, 0.5, 0.5, 0.5, 0],                         // F 立方 → 菱面体 primitive (÷2 面心ベクトル)
        [-0.5, 0.5, 0.5, 0.5, -0.5, 0.5, 0.5, 0.5, -0.5],                // I 立方 → 菱面体 primitive
        [1, 0, 1, -1, 1, 1, 0, -1, 1],                                   // 立方 → 三方 Hex 軸 (a′=a−b, b′=b−c, c′=a+b+c)
        // 260704Cl 追加: F/A/B/C 底心直方晶 → 単斜 C (unique 軸 3 通り)。子の C centering (a′+b′)/2 が親格子ベクトルになる組。
        // 誤った候補は操作集合の完全一致検証で必ず落ちるため、カタログ追加は安全 (同定率が上がるだけ)。
        [1, 0, 0, 0, 0, 0.5, 0, 1, 0.5],                                 // a′=a, b′=c, c′=(b+c)/2
        [0, 1, 0.5, 1, 0, 0, 0, 0, 0.5],                                 // a′=b, b′=a, c′=(a+c)/2
        [0, 0, 0.5, 0, 1, 0.5, 1, 0, 0],                                 // a′=c, b′=b, c′=(a+b)/2
        // 260704Cl 追加: R (Hex 軸) → 三斜 primitive (obverse 菱面体基底: (2/3,1/3,1/3) 系)
        [2.0 / 3, -1.0 / 3, -1.0 / 3, 1.0 / 3, 1.0 / 3, -2.0 / 3, 1.0 / 3, 1.0 / 3, 1.0 / 3],
        // 260704Cl 追加: R (Hex 軸) → 単斜 C (面内 2 回軸が unique 軸。方位の違いは共役ループが吸収)
        [1, 1, -1.0 / 3, 2, 0, -2.0 / 3, 0, 0, 1.0 / 3],
    ];

    private static readonly Lazy<double[][]> ProperSignedPermutations = new(() =>
    {
        var list = new List<double[]>();
        int[][] perms = [[0, 1, 2], [1, 2, 0], [2, 0, 1], [0, 2, 1], [1, 0, 2], [2, 1, 0]];
        foreach (var pm in perms)
            for (int sMask = 0; sMask < 8; sMask++)
            {
                var m = new double[9];
                for (int c = 0; c < 3; c++)
                    m[pm[c] * 3 + c] = ((sMask >> c) & 1) == 0 ? 1 : -1;
                if (Det3(m) > 0.5)
                    list.Add(m);
            }
        return [.. list];
    });

    /// <summary>H (親基準系の操作集合) の空間群型を同定する。成功時 (childSeries, P, p)、失敗時 (-1, null, null)。</summary>
    private static (int Child, double[] P, double[] Shift) Identify(int parentSn, SymmetryOperation[] hOps, List<int[]> mLin)
    {
        //var sigName = SignatureNameMap.Value; // 260705Cl: 未使用の死にローカル (Lazy の全設定走査を無駄に誘発) を削除
        string sig = Signature(mLin);

        // 親格子の生成元 (整数基底 + 中心化) — 恒等線形部の並進から抽出
        var parentCentering = hOps.Where(o => IsIdentity(LinKey(o)))
                                  .Select(o => o.SeitzTranslation)
                                  .Select(t => new[] { Frac(t.U), Frac(t.V), Frac(t.W) })
                                  .ToList();

        // 候補設定: 点群シグネチャが一致する全設定。
        // 260705Cl 修正: cand のみに依存するデータ (中心化格子・線形部キー・(R,t) 対) を基底候補ループ
        // (12 cell × 24 rot) の内側で毎回再構築していたため、候補選定時に 1 度だけ前計算する
        // (FormGroupRelations 初回表示・全 230 群索引構築の応答性改善)。
        var candidates = new List<(int Sn, List<double[]> Lattice, List<int[]> Lin, List<(int[] R, double[] T)> RT)>();
        for (int s = 1; s < SymmetryStatic.TotalSpaceGroupNumber; s++)
        {
            if (SettingSignature(s) != sig) continue;
            var candOps = GetExpandedOps(s);
            var candLattice = candOps.Where(o => IsIdentity(LinKey(o)))
                                     .Select(o => { var t = o.SeitzTranslation; return new[] { Frac(t.U), Frac(t.V), Frac(t.W) }; }).ToList();
            var candLin = new List<int[]>();
            foreach (var op in candOps)
            {
                var key = LinKey(op);
                if (FindKey(candLin, key) < 0) candLin.Add(key);
            }
            candidates.Add((s, candLattice, candLin, [.. candOps.Select(op => (R: LinKey(op), T: TVec(op)))]));
        }

        foreach (var cell in CellChanges)
            foreach (var rot in ProperSignedPermutations.Value)
            {
                var P = MatMulD(cell, rot);
                var Pinv = Invert3(P);
                if (Pinv == null) continue;

                // 260704Cl 追加検証: P の各列 (子基底ベクトルの親座標) は親格子のベクトルでなければならない。
                // これを怠ると det(P)<1 の胞変換 (F→I 等) が非中心化格子の親にも適用され、
                // 「子の整数並進 ⊄ 実格子」のまま mod-1 比較して偽陽性同定を起こす (P6_3/mmc の mmm 部分群が
                // Cmcm でなく Pcmm と誤同定された実バグの原因)。
                bool columnsInLattice = true;
                for (int c = 0; c < 3 && columnsInLattice; c++)
                {
                    var col = new[] { Frac(P[c]), Frac(P[3 + c]), Frac(P[6 + c]) };
                    if (col[0] < Tol && col[1] < Tol && col[2] < Tol) continue; // 整数列は常に格子ベクトル
                    columnsInLattice = parentCentering.Any(pc => NearVec(pc, col));
                }
                if (!columnsInLattice) continue;

                // H の線形部を子基準系へ: 整数行列にならない P は棄却
                var hLinChild = new List<int[]>(mLin.Count);
                bool ok = true;
                foreach (var lm in mLin)
                {
                    var rc = ConjugateInt(Pinv, lm, P);
                    if (rc == null) { ok = false; break; }
                    hLinChild.Add(rc);
                }
                if (!ok) continue;

                // 子基準系の格子 (親格子の像): {0}∪中心化集合
                var childLattice = LatticeCosets(Pinv, parentCentering);
                if (childLattice == null) continue;

                // 260705Cl: 前計算済み候補データを参照するだけにする (旧: ここで GetExpandedOps/格子/線形部を毎回再構築)
                foreach (var (candSn, candLattice, candLin, candRT) in candidates)
                {
                    // 設定側の中心化集合と一致するか (格子の同一性)
                    if (!SameVecSet(childLattice, candLattice)) continue;

                    // 線形部集合の一致 (順不同)
                    if (candLin.Count != hLinChild.Count || candLin.Any(cl => FindKey(hLinChild, cl) < 0)) continue;

                    // 原点シフト q (子基準系) を解いて完全一致検証
                    var q = SolveOriginShift(hOps, P, Pinv, candRT, childLattice);
                    if (q != null)
                    {
                        // p (親座標系) = P·q
                        var p = new[]
                        {
                            P[0] * q[0] + P[1] * q[1] + P[2] * q[2],
                            P[3] * q[0] + P[4] * q[1] + P[5] * q[2],
                            P[6] * q[0] + P[7] * q[1] + P[8] * q[2],
                        };
                        return (candSn, P, p);
                    }
                }
            }
        return (-1, null, null);
    }

    private static readonly Dictionary<int, string> _settingSig = [];
    private static string SettingSignature(int s)
    {
        if (_settingSig.TryGetValue(s, out var sig)) return sig;
        var lin = new List<int[]>();
        foreach (var op in GetExpandedOps(s))
        {
            var key = LinKey(op);
            if (FindKey(lin, key) < 0) lin.Add(key);
        }
        sig = lin.Count == 0 ? "" : Signature(lin);
        _settingSig[s] = sig;
        return sig;
    }

    /// <summary>親格子 (整数基底+中心化) を P⁻¹ で子基準系へ写し、[0,1)³ の剰余類集合 (加法閉) を返す。桁が合わないときは null。</summary>
    private static List<double[]> LatticeCosets(double[] Pinv, List<double[]> parentCentering)
    {
        var gens = new List<double[]>();
        for (int i = 0; i < 3; i++)
            gens.Add([Pinv[i], Pinv[3 + i], Pinv[6 + i]]); // P⁻¹ の列 = 親整数基底の像
        foreach (var c in parentCentering)
            gens.Add([
                Pinv[0] * c[0] + Pinv[1] * c[1] + Pinv[2] * c[2],
                Pinv[3] * c[0] + Pinv[4] * c[1] + Pinv[5] * c[2],
                Pinv[6] * c[0] + Pinv[7] * c[1] + Pinv[8] * c[2]]);

        var set = new List<double[]> { new double[3] };
        var queue = new Queue<double[]>(set);
        while (queue.Count > 0)
        {
            var a = queue.Dequeue();
            foreach (var g in gens)
            {
                var v = new[] { Frac(a[0] + g[0]), Frac(a[1] + g[1]), Frac(a[2] + g[2]) };
                if (!set.Any(x => NearVec(x, v)))
                {
                    if (set.Count > 16) return null; // 想定外の桁 (P が格子を保っていない)
                    set.Add(v);
                    queue.Enqueue(v);
                }
            }
        }
        return set;
    }

    /// <summary>子基準系での原点シフト q を求め、H と候補設定の操作集合が (子格子法で) 完全一致するか検証する。</summary>
    //private static double[] SolveOriginShift(SymmetryOperation[] hOps, double[] P, double[] Pinv, SymmetryOperation[] candOps, List<double[]> lattice) // 260705Cl 旧シグネチャ: candOps から (R,T) を毎回構築していた
    private static double[] SolveOriginShift(SymmetryOperation[] hOps, double[] P, double[] Pinv, List<(int[] R, double[] T)> cand, List<double[]> lattice)
    {
        // H を子基準系へ (q=0): (R_c, t_c0)
        var hChild = new List<(int[] R, double[] T)>();
        foreach (var op in hOps)
        {
            var R = LinKey(op);
            var rc = ConjugateInt(Pinv, R, P);
            var t = op.SeitzTranslation;
            hChild.Add((rc, new[]
            {
                Pinv[0] * t.U + Pinv[1] * t.V + Pinv[2] * t.W,
                Pinv[3] * t.U + Pinv[4] * t.V + Pinv[5] * t.W,
                Pinv[6] * t.U + Pinv[7] * t.V + Pinv[8] * t.W,
            }));
        }
        //var cand = candOps.Select(op => (R: LinKey(op), T: TVec(op))).ToList(); // 260705Cl: 呼び出し元で前計算済み

        // 候補 q の生成: det(R−I) ≠ 0 の R があれば解析解、無ければ 1/24 格子の総当たり
        var qCands = new List<double[]>();
        var pivot = hChild.Select(x => x.R).FirstOrDefault(r => Math.Abs(DetRmI(r)) > 0.5);
        if (pivot != null)
        {
            double det = DetRmI(pivot);
            var adj = AdjRmI(pivot);
            var hT = hChild.Where(x => SameKey(x.R, pivot)).Select(x => x.T).ToList();
            var cT = cand.Where(x => SameKey(x.R, pivot)).Select(x => x.T).ToList();
            foreach (var tc in hT)
                foreach (var ts in cT)
                    for (int nx = -1; nx <= 1; nx++)
                        for (int ny = -1; ny <= 1; ny++)
                            for (int nz = -1; nz <= 1; nz++)
                            {
                                double dx = ts[0] - tc[0] + nx, dy = ts[1] - tc[1] + ny, dz = ts[2] - tc[2] + nz;
                                var q = new[]
                                {
                                    Frac((adj[0] * dx + adj[1] * dy + adj[2] * dz) / det),
                                    Frac((adj[3] * dx + adj[4] * dy + adj[5] * dz) / det),
                                    Frac((adj[6] * dx + adj[7] * dy + adj[8] * dz) / det),
                                };
                                if (!qCands.Any(x => NearVec(x, q)))
                                    qCands.Add(q);
                            }
        }
        else
        {
            for (int i = 0; i < 24; i++)
                for (int j = 0; j < 24; j++)
                    for (int k = 0; k < 24; k++)
                        qCands.Add([i / 24.0, j / 24.0, k / 24.0]);
        }

        // 260705Cl 修正: q に依存しない候補側キー集合 setB を q ループの外で 1 度だけ構築する
        // (旧実装は q ごとに再構築しており、特異ケースの 24³=13,824 総当たりで支配的なコストだった)。
        var setB = BuildKeySet(cand, lattice);
        if (setB == null) return null; // 1/24 格子に乗らない候補は q によらず不一致
        foreach (var q in qCands)
            if (VerifySetEqual(hChild, setB, q, lattice))
                return q;
        return null;
    }

    /// <summary>(R, t) を 1/24 格子スナップの文字列キーへ。スナップに乗らなければ null (= 不一致扱い)。
    /// 260705Cl: VerifySetEqual のローカル関数から昇格 (BuildKeySet と共用)。</summary>
    private static string OpKey(int[] R, double x, double y, double z)
    {
        int qx = (int)Math.Round(Frac(x) * 24), qy = (int)Math.Round(Frac(y) * 24), qz = (int)Math.Round(Frac(z) * 24);
        if (Math.Abs(Frac(x) * 24 - qx) > 1e-3 || Math.Abs(Frac(y) * 24 - qy) > 1e-3 || Math.Abs(Frac(z) * 24 - qz) > 1e-3)
            return null;
        return $"{string.Join(" ", R)}|{qx % 24},{qy % 24},{qz % 24}";
    }

    /// <summary>操作集合を格子剰余類で展開したキー集合を作る。スナップに乗らない操作があれば null。260705Cl 追加。</summary>
    private static HashSet<string> BuildKeySet(List<(int[] R, double[] T)> ops, List<double[]> lattice)
    {
        var set = new HashSet<string>();
        foreach (var (R, T) in ops)
            foreach (var g in lattice)
            {
                var key = OpKey(R, T[0] + g[0], T[1] + g[1], T[2] + g[2]);
                if (key == null) return null;
                set.Add(key);
            }
        return set;
    }

    /// <summary>操作集合の完全一致検証: H 側を原点シフト q で写し、全キーが setB に属し要素数も一致するか
    /// (A⊆B かつ |A|=|B| ⇔ A=B)。260705Cl 修正: 最初の不一致キーで即 false を返す早期棄却方式に。</summary>
    private static bool VerifySetEqual(List<(int[] R, double[] T)> hChild, HashSet<string> setB, double[] q, List<double[]> lattice)
    {
        var setA = new HashSet<string>();
        foreach (var (R, T) in hChild)
        {
            // t′ = t + (R−I)q
            double tx = T[0] + (R[0] - 1) * q[0] + R[1] * q[1] + R[2] * q[2];
            double ty = T[1] + R[3] * q[0] + (R[4] - 1) * q[1] + R[5] * q[2];
            double tz = T[2] + R[6] * q[0] + R[7] * q[1] + (R[8] - 1) * q[2];
            foreach (var g in lattice)
            {
                var key = OpKey(R, tx + g[0], ty + g[1], tz + g[2]);
                if (key == null || !setB.Contains(key)) return false;
                setA.Add(key);
            }
        }
        return setA.Count == setB.Count;
    }
    #endregion

    #region 行列・ベクトル小物 (整数 3×3 は row-major int[9])
    private static int[] LinKey(in SymmetryOperation op)
    {
        var m = SeitzNotation.LinearMatrix(op);
        return [m[0, 0], m[0, 1], m[0, 2], m[1, 0], m[1, 1], m[1, 2], m[2, 0], m[2, 1], m[2, 2]];
    }

    private static double[] TVec(in SymmetryOperation op)
    {
        var t = op.SeitzTranslation;
        return [t.U, t.V, t.W];
    }

    private static int FindKey(List<int[]> list, int[] key)
    {
        for (int i = 0; i < list.Count; i++)
            if (SameKey(list[i], key)) return i;
        return -1;
    }

    private static bool SameKey(int[] a, int[] b)
    {
        for (int i = 0; i < 9; i++)
            if (a[i] != b[i]) return false;
        return true;
    }

    private static bool IsIdentity(int[] m)
        => m[0] == 1 && m[4] == 1 && m[8] == 1 && m[1] == 0 && m[2] == 0 && m[3] == 0 && m[5] == 0 && m[6] == 0 && m[7] == 0;

    private static int[] MatMul(int[] a, int[] b)
    {
        var c = new int[9];
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                c[i * 3 + j] = a[i * 3] * b[j] + a[i * 3 + 1] * b[3 + j] + a[i * 3 + 2] * b[6 + j];
        return c;
    }

    private static double[] MatMulD(double[] a, double[] b)
    {
        var c = new double[9];
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                c[i * 3 + j] = a[i * 3] * b[j] + a[i * 3 + 1] * b[3 + j] + a[i * 3 + 2] * b[6 + j];
        return c;
    }

    private static int Det3(int[] m)
        => m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6]) + m[2] * (m[3] * m[7] - m[4] * m[6]);

    private static double Det3(double[] m)
        => m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6]) + m[2] * (m[3] * m[7] - m[4] * m[6]);

    /// <summary>3×3 (row-major) の逆行列。特異なら null。260705Cl: GroupRelation.GetInverseTransform から使うため internal 化。</summary>
    internal static double[] Invert3(double[] m)
    {
        double det = Det3(m);
        if (Math.Abs(det) < 1e-12) return null;
        return
        [
            (m[4] * m[8] - m[5] * m[7]) / det, (m[2] * m[7] - m[1] * m[8]) / det, (m[1] * m[5] - m[2] * m[4]) / det,
            (m[5] * m[6] - m[3] * m[8]) / det, (m[0] * m[8] - m[2] * m[6]) / det, (m[2] * m[3] - m[0] * m[5]) / det,
            (m[3] * m[7] - m[4] * m[6]) / det, (m[1] * m[6] - m[0] * m[7]) / det, (m[0] * m[4] - m[1] * m[3]) / det,
        ];
    }

    /// <summary>P⁻¹ · R · P。結果が整数行列でなければ null。</summary>
    private static int[] ConjugateInt(double[] Pinv, int[] R, double[] P)
    {
        var Rd = new double[9];
        for (int i = 0; i < 9; i++) Rd[i] = R[i];
        var c = MatMulD(MatMulD(Pinv, Rd), P);
        var r = new int[9];
        for (int i = 0; i < 9; i++)
        {
            r[i] = (int)Math.Round(c[i]);
            if (Math.Abs(c[i] - r[i]) > 1e-6) return null;
        }
        return r;
    }

    /// <summary>det(R − I)。</summary>
    private static double DetRmI(int[] r)
    {
        var m = new double[] { r[0] - 1, r[1], r[2], r[3], r[4] - 1, r[5], r[6], r[7], r[8] - 1 };
        return Det3(m);
    }

    /// <summary>adj(R − I) (余因子転置)。(R−I)⁻¹ = adj/det。</summary>
    private static double[] AdjRmI(int[] r)
    {
        double a = r[0] - 1, b = r[1], c = r[2], d = r[3], e2 = r[4] - 1, f = r[5], g = r[6], h = r[7], i2 = r[8] - 1;
        return
        [
            e2 * i2 - f * h, c * h - b * i2, b * f - c * e2,
            f * g - d * i2, a * i2 - c * g, c * d - a * f,
            d * h - e2 * g, b * g - a * h, a * e2 - b * d,
        ];
    }

    private static double Frac(double d)
    {
        d -= Math.Floor(d);
        return d > 1 - Tol ? 0 : d;
    }

    private static bool NearVec(double[] a, double[] b)
        => FracDist(a[0], b[0]) < 1e-4 && FracDist(a[1], b[1]) < 1e-4 && FracDist(a[2], b[2]) < 1e-4;

    private static double FracDist(double a, double b)
    {
        double d = Math.Abs(Frac(a) - Frac(b));
        return Math.Min(d, 1 - d);
    }

    private static bool SameVecSet(List<double[]> a, List<double[]> b)
        => a.Count == b.Count && a.All(x => b.Any(y => NearVec(x, y)));

    private static List<(double X, double Y, double Z)> GenerateOrbit(IReadOnlyList<SymmetryOperation> ops, double x, double y, double z)
    {
        var pts = new List<(double X, double Y, double Z)>();
        foreach (var op in ops)
        {
            var (px, py, pz) = op.ApplyMatrix(x, y, z);
            var t = op.SeitzTranslation;
            var v = (X: Frac(px + t.U), Y: Frac(py + t.V), Z: Frac(pz + t.W));
            if (!pts.Any(q => Near(q, v)))
                pts.Add(v);
        }
        return pts;
    }

    private static bool Near((double X, double Y, double Z) a, (double X, double Y, double Z) b)
        => FracDist(a.X, b.X) < 1e-4 && FracDist(a.Y, b.Y) < 1e-4 && FracDist(a.Z, b.Z) < 1e-4;
    #endregion
}
