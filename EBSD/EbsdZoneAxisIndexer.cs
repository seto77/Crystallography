#region using
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using V3 = OpenTK.Mathematics.Vector3d;
#endregion

namespace Crystallography;

/// <summary>手動で拾ったバンド交点 (晶帯軸) 1 点。260921Cl 追加</summary>
/// <param name="Col">実測画像のピクセル座標 (コーナー原点)</param>
/// <param name="FixedIndex">作者が指数を手入力した場合の [uvw]。null なら自動で割り当てる</param>
public readonly record struct EbsdZoneAxisPick(double Col, double Row, (int U, int V, int W)? FixedIndex = null);

/// <summary>晶帯軸からの方位 (と検出器幾何) の解。260921Cl 追加</summary>
public sealed class EbsdZoneAxisSolution
{
    /// <summary>結晶 → 試料系の回転 (Crystal.RotationMatrix と同じ規約)</summary>
    public Matrix3D Rotation;
    /// <summary>幾何も最適化した場合は更新後の幾何。していなければ入力と同じ</summary>
    public EbsdDetectorGeometry Geometry;
    /// <summary>拾った点ごとの割当 [uvw] (未割当は null)</summary>
    public (int U, int V, int W)?[] Assignments = [];
    /// <summary>割り当てられた点の数</summary>
    public int Assigned;
    /// <summary>割当点の残差 (画像ピクセル rms と角度 rms)</summary>
    public double RmsPx, RmsDeg;
    /// <summary>順位付け用スコア (大きいほど良い)</summary>
    public double Score;
    /// <summary>幾何を動かした量 (検出器中心 mm。動かしていなければ 0)</summary>
    public double GeometryShiftMm;

    /// <summary>260921Cl 追加: 割当の表示文字列。「点番号 (1 始まり):[u v w]」を空白区切りで並べ、未割当の点は - になる</summary>
    public string AssignmentText => string.Join(" ", Assignments.Select((a, i) =>
        a is null ? $"{i + 1}:-" : $"{i + 1}:[{a.Value.U} {a.Value.V} {a.Value.W}]"));
}

/// <summary>
/// 手動で拾ったバンド交点 (晶帯軸) から結晶方位を決める。260921Cl 追加。GUI 非依存。
///
/// 【なぜ要るか】ノイズの多い実測パターンでは Radon も Dictionary も動かないが、バンドの交点は目視で拾える。
/// 晶帯軸は 1 点で 2 自由度を拘束するので、バンド (1 点あたり 1 自由度) より情報量が多く、3 点あれば方位が決まる。
///
/// 【手順】
///   ① 晶帯軸カタログを作る (反射の整数外積 → 既約 [uvw]。重みはその晶帯に属する反射の相対強度の和 = 交点の目立ちやすさ)
///   ② 拾った点を試料系の「出ていく方向」へ変換する
///   ③ 2 点の対の角度がカタログ側の対の角度と合う組み合わせを総当たりし、Wahba/Kabsch で回転を作る
///   ④ 全点を採点 → 上位を反復精緻化 → 重複除去
///   ⑤ (任意) 方位と検出器幾何 (検出器中心 3 成分) を同時に Nelder-Mead で最適化し、予測位置と拾った位置の
///      ピクセル残差を最小化する
///
/// 指数を手入力した点があれば、その点の割当をその [uvw] (と符号反転) に限定する。2 点以上が手入力なら
/// 探索そのものを飛ばして ③ の組み合わせをその 2 点に絞る。
/// </summary>
public static class EbsdZoneAxisIndexer
{
    /// <summary>対の角度がこれより小さい/大きい組は方位を決める力が弱いので使わない</summary>
    static readonly double MinPairAngle = 8 * Math.PI / 180, MaxPairAngle = 172 * Math.PI / 180;

    #region 晶帯軸カタログ

    sealed class ZoneNode
    {
        public V3 Dir;                    //結晶直交座標系の単位ベクトル
        public (int U, int V, int W) Index;
        public double Weight;             //この晶帯に属する反射の相対強度の和 (交点の目立ちやすさ)
    }

    /// <summary>反射リスト (VectorOfG_KikuchiLine) から晶帯軸カタログを作る。重み上位 maxNodes 本を ± 両方向で返す。
    /// ⚠ maxIndex で指数の上限を切るのが要点。EBSD で交点として見えるのは低指数の晶帯軸だけで、
    /// 高指数まで許すと**どんな点でも説明できてしまい**、残差の小さい偽解が上位に来る
    /// (260921Cl に Botallackite の実測で実際に起きた: [10 -2 5] のような軸が割り当たった)。</summary>
    static List<ZoneNode> BuildCatalog(Crystal crystal, int maxNodes, int maxIndex)
    {
        var gList = crystal.VectorOfG_KikuchiLine;
        if (gList == null || gList.Count < 2) return [];

        //① 反射の整数外積で晶帯軸を列挙 (FormEBSD の晶帯軸ラベルと同じ作り方)
        var seen = new HashSet<(int, int, int)>();
        for (int i = 0; i < gList.Count - 1; i++)
        {
            var (h1, k1, l1) = gList[i].Index;
            for (int j = i + 1; j < gList.Count; j++)
            {
                var (h2, k2, l2) = gList[j].Index;
                int u = k1 * l2 - l1 * k2, v = l1 * h2 - h1 * l2, w = h1 * k2 - k1 * h2;
                if (u == 0 && v == 0 && w == 0) continue;
                int gcd = Algebra.Irreducible(u, v, w);
                if (gcd > 1) { u /= gcd; v /= gcd; w /= gcd; }
                if (u < 0 || (u == 0 && v < 0) || (u == 0 && v == 0 && w < 0)) { u = -u; v = -v; w = -w; }
                if (Math.Abs(u) > maxIndex || Math.Abs(v) > maxIndex || Math.Abs(w) > maxIndex) continue;
                seen.Add((u, v, w));
            }
        }
        if (seen.Count == 0) return [];

        //② 重み = その晶帯に属する反射 (g・[uvw] = 0) の相対強度の和
        var nodes = new List<ZoneNode>(seen.Count);
        foreach (var (u, v, w) in seen)
        {
            double weight = 0;
            foreach (var g in gList)
            {
                var (h, k, l) = g.Index;
                if (h * u + k * v + l * w == 0) weight += g.RelativeIntensity;
            }
            var real = u * crystal.A_Axis + v * crystal.B_Axis + w * crystal.C_Axis;
            double len = Math.Sqrt(real.X * real.X + real.Y * real.Y + real.Z * real.Z);
            if (len < 1E-12) continue;
            nodes.Add(new ZoneNode { Dir = new V3(real.X / len, real.Y / len, real.Z / len), Index = (u, v, w), Weight = weight });
        }

        //③ 重み上位を ± 両方向に展開する。観測方向には実際の符号があるので、カタログ側も符号付きで持つ
        var top = nodes.OrderByDescending(n => n.Weight).Take(Math.Max(2, maxNodes)).ToList();
        var result = new List<ZoneNode>(top.Count * 2);
        foreach (var n in top)
        {
            result.Add(n);
            result.Add(new ZoneNode { Dir = -n.Dir, Index = (-n.Index.U, -n.Index.V, -n.Index.W), Weight = n.Weight });
        }
        return result;
    }

    #endregion

    /// <summary>拾った晶帯軸から方位候補を返す (スコア降順)。</summary>
    /// <param name="picks">拾った点。2 点以上必要 (幾何も最適化するなら 4 点以上を推奨)</param>
    /// <param name="geometry">実測画像のピクセルグリッドを基準にした検出器幾何</param>
    /// <param name="crystal">VectorOfG_KikuchiLine を設定済みの結晶</param>
    /// <param name="toleranceDeg">角度の許容差。拾いの誤差 5 px は約 0.4° (カメラ長 35 mm・0.05 mm/px) なので既定 2°</param>
    /// <param name="refineGeometry">true で検出器中心 (DetX/DetY/DetZ) も同時に最適化する</param>
    /// <param name="maxNodes">カタログに載せる晶帯軸の本数 (重み上位から)。± 展開するので実際の要素数は 2 倍</param>
    /// <param name="maxCandidates">返す候補の最大数</param>
    /// <param name="properSymmetries">結晶点群の proper 回転。重複候補の除去に使う。null なら素の misorientation で比べる</param>
    /// <param name="maxIndex">カタログに載せる晶帯軸指数の上限。
    ///   ⚠ <b>ここを緩めると偽解が上位に来る。</b> 高指数の軸はどんな点でも説明できてしまうため
    ///   (実測で [10 -2 5] が割り当たった)。既定 4 より上げるときは必ずオーバーレイで目視確認すること</param>
    /// <param name="cancel">中止トークン</param>
    public static List<EbsdZoneAxisSolution> Index(
        IReadOnlyList<EbsdZoneAxisPick> picks, EbsdDetectorGeometry geometry, Crystal crystal,
        double toleranceDeg = 2.0, int maxNodes = 300, int maxCandidates = 10,
        bool refineGeometry = false, Matrix3D[] properSymmetries = null,
        int maxIndex = 4, CancellationToken cancel = default)
    {
        ArgumentNullException.ThrowIfNull(picks);
        ArgumentNullException.ThrowIfNull(geometry);
        ArgumentNullException.ThrowIfNull(crystal);
        if (picks.Count < 2) return [];
        if (!(toleranceDeg > 0)) throw new ArgumentOutOfRangeException(nameof(toleranceDeg));

        var catalog = BuildCatalog(crystal, maxNodes, maxIndex);
        if (catalog.Count < 2) return [];
        double tol = toleranceDeg * Math.PI / 180;
        double maxWeight = catalog.Max(n => n.Weight);
        if (!(maxWeight > 0)) maxWeight = 1;

        //拾った点 → 試料系の「出ていく方向」。PixelToSampleDirection は視線 (-P̂) なので符号を反転する
        var obs = new V3[picks.Count];
        for (int i = 0; i < picks.Count; i++)
            obs[i] = -geometry.PixelToSampleDirection(picks[i].Col, picks[i].Row);

        //手入力された指数は、その点で許すカタログ添字を限定する
        var allowed = new List<int>[picks.Count];
        for (int i = 0; i < picks.Count; i++)
        {
            if (picks[i].FixedIndex is not { } fx) { allowed[i] = null; continue; }
            allowed[i] = [];
            for (int c = 0; c < catalog.Count; c++)
                if (SameAxis(catalog[c].Index, fx)) allowed[i].Add(c);
            //カタログに無い指数を入力された場合はその場で作る (作者が高指数を知っている場合)
            if (allowed[i].Count == 0)
            {
                var real = fx.U * crystal.A_Axis + fx.V * crystal.B_Axis + fx.W * crystal.C_Axis;
                double len = Math.Sqrt(real.X * real.X + real.Y * real.Y + real.Z * real.Z);
                if (len > 1E-12)
                {
                    var d = new V3(real.X / len, real.Y / len, real.Z / len);
                    catalog.Add(new ZoneNode { Dir = d, Index = fx, Weight = maxWeight });
                    catalog.Add(new ZoneNode { Dir = -d, Index = (-fx.U, -fx.V, -fx.W), Weight = maxWeight });
                    allowed[i] = [catalog.Count - 2, catalog.Count - 1];
                }
                else allowed[i] = null;
            }
        }

        //③ 対の角度で総当たり。手入力が 2 点以上あればその組だけを種にする
        var fixedIdx = Enumerable.Range(0, picks.Count).Where(i => allowed[i] != null).ToArray();
        var seedPairs = new List<(int I, int J)>();
        for (int i = 0; i < picks.Count - 1; i++)
            for (int j = i + 1; j < picks.Count; j++)
            {
                if (fixedIdx.Length >= 2 && !(allowed[i] != null && allowed[j] != null)) continue;
                double a = Angle(obs[i], obs[j]);
                if (a < MinPairAngle || a > MaxPairAngle) continue;
                seedPairs.Add((i, j));
            }
        if (seedPairs.Count == 0) return [];

        var solutions = new List<EbsdZoneAxisSolution>();
        var rotated = new V3[catalog.Count]; //260921Cl: Evaluate が使う作業バッファ (逐次呼び出しなので 1 本で足りる)
        foreach (var (i, j) in seedPairs)
        {
            cancel.ThrowIfCancellationRequested();
            double target = Angle(obs[i], obs[j]);
            //260921Cl 変更: |角度 − target| ≤ tol は、cos が [0,π] で単調減少なので内積の窓で厳密に同値。
            //  旧は Acos + sqrt×2 + 除算をカタログ対の総当たり (種ペア × catalog²) で回していた
            //  (点 10 個・カタログ 600 本で 1.6e7 回)。Dir も obs も単位ベクトルなので長さの除算も要らない
            double dotLo = Math.Cos(target + tol), dotHi = Math.Cos(target - tol);
            //260921Cl 変更: 手入力が無いときは全カタログ。種ペアごとに List を作らず添字で回す
            var listI = allowed[i]; var listJ = allowed[j];
            int nI = listI?.Count ?? catalog.Count, nJ = listJ?.Count ?? catalog.Count;
            for (int ti = 0; ti < nI; ti++)
                for (int tj = 0; tj < nJ; tj++)
                {
                    int ca = listI is null ? ti : listI[ti], cb = listJ is null ? tj : listJ[tj];
                    if (ca == cb) continue;
                    double dt = V3.Dot(catalog[ca].Dir, catalog[cb].Dir);
                    if (dt < dotLo || dt > dotHi) continue;
                    var r = EbsdIndexer.SolveWahba([(obs[i], catalog[ca].Dir, 1.0), (obs[j], catalog[cb].Dir, 1.0)]);
                    if (r == null) continue;
                    var sol = Evaluate(r, obs, catalog, allowed, tol, maxWeight, geometry, picks, rotated);
                    if (sol != null) solutions.Add(sol);
                }
        }
        if (solutions.Count == 0) return [];

        //④ 上位を反復精緻化 (割当 → Wahba → 再割当)
        var refined = new List<EbsdZoneAxisSolution>();
        foreach (var s in solutions.OrderByDescending(x => x.Score).Take(80))
        {
            cancel.ThrowIfCancellationRequested();
            var cur = s;
            for (int iter = 0; iter < 3; iter++)
            {
                var eqs = new List<(V3 m, V3 g, double w)>();
                for (int k = 0; k < picks.Count; k++)
                {
                    if (cur.Assignments[k] is not { } uvw) continue;
                    var node = catalog.FirstOrDefault(n => n.Index == uvw);
                    if (node == null) continue;
                    eqs.Add((obs[k], node.Dir, 1.0));
                }
                if (eqs.Count < 2) break;
                var r = EbsdIndexer.SolveWahba([.. eqs]);
                if (r == null) break;
                var next = Evaluate(r, obs, catalog, allowed, tol, maxWeight, geometry, picks, rotated);
                if (next == null || next.Score <= cur.Score) break;
                cur = next;
            }
            refined.Add(cur);
        }

        //重複除去 (対称等価も含めて) → 上位だけ返す
        var unique = new List<EbsdZoneAxisSolution>();
        foreach (var s in refined.OrderByDescending(x => x.Score))
        {
            if (unique.Any(u => Misorientation(u.Rotation, s.Rotation, properSymmetries) < 1.0)) continue;
            unique.Add(s);
            if (unique.Count >= maxCandidates) break;
        }

        //⑤ 幾何も最適化する (方位 3 + 検出器中心 3)
        if (refineGeometry)
            for (int k = 0; k < unique.Count; k++)
            {
                cancel.ThrowIfCancellationRequested();
                var better = RefineWithGeometry(unique[k], picks, catalog, geometry);
                if (better != null) unique[k] = better;
            }

        return [.. unique.OrderByDescending(s => s.Score)];
    }

    /// <summary>回転 r について全点を採点し、割当と残差を詰めた解を返す。手入力の指数に反すれば null。</summary>
    /// <param name="rotated">回転済みカタログを書き込む作業バッファ (長さ = catalog.Count)。呼び出し側で 1 本確保して使い回す</param>
    static EbsdZoneAxisSolution Evaluate(Matrix3D r, V3[] obs, List<ZoneNode> catalog, List<int>[] allowed,
        double tol, double maxWeight, EbsdDetectorGeometry geometry, IReadOnlyList<EbsdZoneAxisPick> picks, V3[] rotated)
    {
        var assign = new (int U, int V, int W)?[obs.Length];
        double score = 0, sumSqDeg = 0;
        int assigned = 0;
        //260921Cl 追加: 回転済みカタログを 1 回だけ作る。旧は拾った点の数だけ r * Dir を作り直していた
        //  (点 6〜10 × カタログ 600 で、1 回の Evaluate に matvec が数千回。Evaluate 自体が種ペアの数だけ回る)
        //  ⚠ バッファは呼び出し側から借りる。Evaluate は種ペア × 合致対の数だけ呼ばれるので、
        //  ここで new すると 600 要素 × 数万回 = GB 級の Gen0 割当になる (260921Cl の /simplify2 で判明)
        for (int c = 0; c < catalog.Count; c++) rotated[c] = r * catalog[c].Dir;
        //最近傍は「角度が最小」= 「内積が最大」。obs もカタログ Dir も単位ベクトルで回転も長さを変えないので、
        //  Angle() の Acos と 2 回の sqrt は選ぶだけなら不要。採用した 1 本にだけ Angle を掛けて従来と同じ値を使う
        double cosTol = Math.Cos(tol);
        for (int k = 0; k < obs.Length; k++)
        {
            var candidates = allowed[k];
            double bestDot = double.MinValue; int bestC = -1;
            int n = candidates?.Count ?? catalog.Count;
            for (int t = 0; t < n; t++)
            {
                int c = candidates is null ? t : candidates[t];
                double d = V3.Dot(obs[k], rotated[c]);
                if (d > bestDot) { bestDot = d; bestC = c; }
            }
            if (bestC < 0 || bestDot < cosTol) continue;
            double bestAng = Angle(obs[k], rotated[bestC]);
            if (bestAng > tol) continue;
            assign[k] = catalog[bestC].Index;
            assigned++;
            double deg = bestAng * 180 / Math.PI;
            sumSqDeg += deg * deg;
            //低指数 (= 目立つ) 晶帯軸への割当を優遇する。高指数の軸はどんな点でも説明できてしまうため
            score += (1 - bestAng / tol) + 0.5 * catalog[bestC].Weight / maxWeight;
        }
        //手入力した点が割り当たらなかった解は棄却する
        for (int k = 0; k < obs.Length; k++)
            if (allowed[k] != null && assign[k] is null) return null;
        if (assigned < 2) return null;
        //⚠ 画面上の別々の点が同じ晶帯軸に割り当たることはあり得ない。そういう解は棄却する
        //  (260921Cl: 近接した 2 点が同じ軸を取り合って偽解が上位に来た)
        for (int a = 0; a < assign.Length; a++)
            for (int b = a + 1; b < assign.Length; b++)
                if (assign[a] is { } x && assign[b] is { } y && SameAxis(x, y)) return null;

        var sol = new EbsdZoneAxisSolution
        {
            Rotation = r,
            Geometry = geometry,
            Assignments = assign,
            Assigned = assigned,
            RmsDeg = Math.Sqrt(sumSqDeg / assigned),
            Score = score + 2.0 * assigned, //まず割当数、次に質
        };
        sol.RmsPx = PixelRms(r, geometry, picks, assign, catalog);
        return sol;
    }

    /// <summary>割当済みの点について、予測位置と拾った位置のピクセル rms を返す。</summary>
    static double PixelRms(Matrix3D r, EbsdDetectorGeometry geom, IReadOnlyList<EbsdZoneAxisPick> picks,
        (int U, int V, int W)?[] assign, List<ZoneNode> catalog)
    {
        double sum = 0; int n = 0;
        for (int k = 0; k < picks.Count; k++)
        {
            if (assign[k] is not { } uvw) continue;
            var node = catalog.FirstOrDefault(x => x.Index == uvw);
            if (node == null) continue;
            var p = geom.SampleDirectionToPixel(r * node.Dir);
            if (p is not { } q) return double.NaN;
            double dx = q.Col - picks[k].Col, dy = q.Row - picks[k].Row;
            sum += dx * dx + dy * dy; n++;
        }
        return n == 0 ? double.NaN : Math.Sqrt(sum / n);
    }

    /// <summary>方位 3 + 検出器中心 3 を Nelder-Mead で同時に最適化し、ピクセル残差を最小にする。
    /// 割当は固定したまま動かす。点が 4 点未満なら幾何は動かさない (自由度が足りない)。</summary>
    static EbsdZoneAxisSolution RefineWithGeometry(EbsdZoneAxisSolution s, IReadOnlyList<EbsdZoneAxisPick> picks,
        List<ZoneNode> catalog, EbsdDetectorGeometry geom0)
    {
        if (s.Assigned < 4) return null;
        var pairs = new List<(V3 g, double Col, double Row)>();
        for (int k = 0; k < picks.Count; k++)
        {
            if (s.Assignments[k] is not { } uvw) continue;
            var node = catalog.FirstOrDefault(x => x.Index == uvw);
            if (node != null) pairs.Add((node.Dir, picks[k].Col, picks[k].Row));
        }
        if (pairs.Count < 4) return null;

        double Objective(double[] p)
        {
            var r = EbsdIndexer.PerturbRotation(s.Rotation, p[0], p[1], p[2]);
            EbsdDetectorGeometry g;
            try
            {
                g = new EbsdDetectorGeometry(geom0.DetTilt, geom0.DetX + p[3], geom0.DetY + p[4], geom0.DetZ + p[5],
                    geom0.PixelSize, geom0.WidthPx, geom0.HeightPx, geom0.XMirror, geom0.SampleTilt);
            }
            catch { return 1E9; }
            double sum = 0;
            foreach (var (gv, col, row) in pairs)
            {
                var q = g.SampleDirectionToPixel(r * gv);
                if (q is not { } t) return 1E9;
                double dx = t.Col - col, dy = t.Row - row;
                sum += dx * dx + dy * dy;
            }
            return Math.Sqrt(sum / pairs.Count);
        }

        var (best, value, _) = EbsdPatternScorer.NelderMeadRestart(Objective,
            [0, 0, 0, 0, 0, 0], [0.5, 0.5, 0.5, 0.3, 0.3, 0.3], maxEvalPerRun: 600, restarts: 2);
        if (!(value < s.RmsPx) || !double.IsFinite(value)) return null;

        var rot = EbsdIndexer.PerturbRotation(s.Rotation, best[0], best[1], best[2]);
        EbsdDetectorGeometry geom;
        try
        {
            geom = new EbsdDetectorGeometry(geom0.DetTilt, geom0.DetX + best[3], geom0.DetY + best[4], geom0.DetZ + best[5],
                geom0.PixelSize, geom0.WidthPx, geom0.HeightPx, geom0.XMirror, geom0.SampleTilt);
        }
        catch { return null; }

        double sumSqDeg = 0; int n = 0;
        for (int k = 0; k < picks.Count; k++)
        {
            if (s.Assignments[k] is not { } uvw) continue;
            var node = catalog.FirstOrDefault(x => x.Index == uvw);
            if (node == null) continue;
            double deg = Angle(-geom.PixelToSampleDirection(picks[k].Col, picks[k].Row), rot * node.Dir) * 180 / Math.PI;
            sumSqDeg += deg * deg; n++;
        }
        return new EbsdZoneAxisSolution
        {
            Rotation = rot,
            Geometry = geom,
            Assignments = s.Assignments,
            Assigned = s.Assigned,
            RmsPx = value,
            RmsDeg = n > 0 ? Math.Sqrt(sumSqDeg / n) : s.RmsDeg,
            Score = s.Score,
            GeometryShiftMm = Math.Sqrt(best[3] * best[3] + best[4] * best[4] + best[5] * best[5]),
        };
    }

    #region 小物

    static double Angle(in V3 a, in V3 b)
        => Math.Acos(Math.Clamp(V3.Dot(a, b) / (a.Length * b.Length), -1, 1));

    static bool SameAxis((int U, int V, int W) a, (int U, int V, int W) b)
        => (a.U == b.U && a.V == b.V && a.W == b.W) || (a.U == -b.U && a.V == -b.V && a.W == -b.W);

    /// <summary>対称操作を考慮した misorientation [deg]。symmetries が null なら素の値</summary>
    static double Misorientation(Matrix3D r1, Matrix3D r2, Matrix3D[] symmetries)
    {
        double best = EbsdIndexer.MisorientationDeg(r1, r2);
        if (symmetries == null) return best;
        foreach (var s in symmetries)
            best = Math.Min(best, EbsdIndexer.MisorientationDeg(r1, r2 * s));
        return best;
    }

    #endregion
}
