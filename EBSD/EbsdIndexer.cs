#region using
using System;
using System.Collections.Generic;
using System.Linq;
using V3 = OpenTK.Mathematics.Vector3d;
#endregion

namespace Crystallography;

/// <summary>指数付け候補 (結晶方位)。260724Cl 追加</summary>
public sealed class EbsdOrientationCandidate
{
    /// <summary>結晶→試料系の回転行列 (Crystal.RotationMatrix と同じ規約。FormMain.SetRotation へ渡せる)</summary>
    public Matrix3D Rotation;

    /// <summary>総合スコア (大きいほど良い)</summary>
    public double Score;

    /// <summary>割り当てられたバンド数 / 検出バンド数</summary>
    public int AssignedBands, TotalBands;

    /// <summary>割当バンドの角度残差 RMS (deg)</summary>
    public double AngularRmsDeg;

    /// <summary>バンド index → 割当反射 (h,k,l)。未割当は含まない</summary>
    public Dictionary<int, (int H, int K, int L)> Assignments = [];

    /// <summary>ZNCC (Phase 3 で充填。未計算は NaN)</summary>
    public double Zncc = double.NaN;

    public string AssignmentText => string.Join(", ", Assignments.OrderBy(p => p.Key).Select(p => $"{p.Key}:({p.Value.H} {p.Value.K} {p.Value.L})"));

    /// <summary>hkl のみの表示 (Radon 方位探索用: キーは強度順ランクでバンド番号でないため)。260724Cl 追加</summary>
    public string HklText => string.Join(" ", Assignments.OrderBy(p => p.Key).Select(p => $"({p.Value.H} {p.Value.K} {p.Value.L})"));
}

/// <summary>
/// 検出バンドからの結晶方位候補探索 (pair-angle マッチ + Kabsch + 反復割当)。260724Cl 追加。GUI 非依存。
/// 参照カタログは VectorOfG_KikuchiLine から ±g・同方向調和反射をまとめた方向ノード (PlaneCatalog) を構築して使う。
/// </summary>
public static class EbsdIndexer
{
    /// <summary>pair-angle 対応に使う面間角の範囲 (rad)。260724Cl (/simplify): 理論側・観測側の 2 箇所に散っていた 8°/85° を一元化</summary>
    static readonly double MinPairAngle = 8 * Math.PI / 180, MaxPairAngle = 85 * Math.PI / 180;

    #region PlaneCatalog
    sealed class PlaneNode
    {
        public V3 Direction;                 //結晶直交座標系の単位ベクトル (±は辞書順で正規化)
        public (int H, int K, int L)[] Indices = [];  //この方向の反射 (低次から)
        public double[] DValues = [];
        public double BestIntensity;
    }

    /// <summary>±g・同方向調和反射を 1 ノードにまとめ、強度上位 maxNodes 方向に制限したカタログを作る</summary>
    static List<PlaneNode> BuildCatalog(IEnumerable<Vector3D> kikuchiReflections, int maxNodes)
    {
        var nodes = new List<PlaneNode>();
        foreach (var g in kikuchiReflections)
        {
            var len = Math.Sqrt(g.X * g.X + g.Y * g.Y + g.Z * g.Z);
            if (len < 1E-12) continue;
            var dir = new V3(g.X / len, g.Y / len, g.Z / len);
            //±を統一 (辞書順で正の側へ)
            if (dir.Z < 0 || (dir.Z == 0 && dir.Y < 0) || (dir.Z == 0 && dir.Y == 0 && dir.X < 0)) dir = -dir;

            var node = nodes.FirstOrDefault(n => V3.Dot(n.Direction, dir) > 0.99995); //~0.6°
            if (node == null)
            {
                node = new PlaneNode { Direction = dir };
                nodes.Add(node);
            }
            node.Indices = [.. node.Indices, (g.Index.h, g.Index.k, g.Index.l)];
            node.DValues = [.. node.DValues, g.d];
            node.BestIntensity = Math.Max(node.BestIntensity, g.RelativeIntensity);
        }
        return [.. nodes.OrderByDescending(n => n.BestIntensity).Take(maxNodes)];
    }
    #endregion

    #region バンド観測 (試料系法線 + sinθB 統計)
    sealed class BandObservation
    {
        public int BandIndex;
        public V3 NormalSample;   //符号不定・正規化
        public double SinBragg;    //エッジ点群の median (エッジ不明瞭なら NaN)
        public double SinBraggSigma;
        public double Weight;      //センター品質由来
    }

    static List<BandObservation> BuildObservations(IReadOnlyList<EbsdBand> bands, EbsdDetectorGeometry geometry)
    {
        var obs = new List<BandObservation>();
        for (int i = 0; i < bands.Count; i++)
        {
            var b = bands[i];
            var gLab = geometry.LineToLabNormal(b.LineA, b.LineB, b.LineC);
            var sinArr = b.EdgePoints.Select(e => geometry.SinBraggFromEdgePoint(gLab, e.Col, e.Row)).OrderBy(v => v).ToArray();
            double sinB = double.NaN, sigma = double.NaN;
            if (sinArr.Length >= 2 && b.WidthQuality > 0.05)
            {
                sinB = sinArr[sinArr.Length / 2];
                var mad = sinArr.Select(v => Math.Abs(v - sinB)).OrderBy(v => v).ElementAt(sinArr.Length / 2);
                sigma = Math.Max(1.4826 * mad, sinB * 0.2); //最低 20% 相対
            }
            obs.Add(new BandObservation
            {
                BandIndex = i,
                NormalSample = geometry.LabToSample(gLab).Normalized(),
                SinBragg = sinB,
                SinBraggSigma = sigma,
                Weight = 0.3 + 0.7 * b.CenterQuality,
            });
        }
        return obs;
    }
    #endregion

    /// <summary>
    /// バンド群から方位候補を探索する。
    /// kikuchiReflections = Crystal.VectorOfG_KikuchiLine、waveLength = 電子波長 (Å)。
    /// 返り値はスコア降順 (misorientation 3° 以内の重複は除去済み。対称等価は含まれ得る)。
    /// </summary>
    public static List<EbsdOrientationCandidate> Index(IReadOnlyList<EbsdBand> bands, EbsdDetectorGeometry geometry,
        IEnumerable<Vector3D> kikuchiReflections, double waveLength,
        int maxCandidates = 10, int maxCatalogNodes = 120, double pairTolDeg = 1.5, double assignTolDeg = 2.5,
        System.Threading.CancellationToken cancel = default)
    {
        var catalog = BuildCatalog(kikuchiReflections, maxCatalogNodes);
        var obs = BuildObservations(bands, geometry);
        if (catalog.Count < 2 || obs.Count < 2) return [];

        //260724Cl (/simplify): ScoreAndAssign は seed 毎に数千回呼ばれるため、作業バッファを 1 回確保して使い回す (GC churn 削減。呼び出しは逐次)
        var predictedBuf = new V3[catalog.Count];
        var usedObsBuf = new bool[obs.Count];
        var usedNodeBuf = new bool[catalog.Count];

        #region 理論 pair-angle テーブル (角度昇順ソート)
        var pairs = new List<(int P, int Q, double Angle)>();
        for (int p = 0; p < catalog.Count; p++)
            for (int q = p + 1; q < catalog.Count; q++)
            {
                double ang = Math.Acos(Math.Clamp(Math.Abs(V3.Dot(catalog[p].Direction, catalog[q].Direction)), 0, 1));
                if (ang > MinPairAngle && ang < MaxPairAngle)
                    pairs.Add((p, q, ang));
            }
        var pairAngles = pairs.OrderBy(p => p.Angle).ToArray();
        var angleKeys = pairAngles.Select(p => p.Angle).ToArray();
        #endregion

        #region 観測ペア → 理論ペア対応 → seed 回転 → 全バンド割当スコア
        double pairTol = pairTolDeg * Math.PI / 180, assignTol = assignTolDeg * Math.PI / 180;
        var candidates = new List<EbsdOrientationCandidate>();
        //観測ペアは品質上位から (組合せ爆発防止のため上位 8 バンドまで)
        var topObs = obs.OrderByDescending(o => o.Weight).Take(8).ToList();

        for (int i = 0; i < topObs.Count; i++)
            for (int j = i + 1; j < topObs.Count; j++)
            {
                cancel.ThrowIfCancellationRequested();
                var oi = topObs[i]; var oj = topObs[j];
                double angObs = Math.Acos(Math.Clamp(Math.Abs(V3.Dot(oi.NormalSample, oj.NormalSample)), 0, 1));
                if (angObs < MinPairAngle || angObs > MaxPairAngle) continue;

                int lo = LowerBound(angleKeys, angObs - pairTol), hi = UpperBound(angleKeys, angObs + pairTol);
                for (int k = lo; k < hi; k++)
                {
                    var (p, q, _) = pairAngles[k];
                    //(i→p, j→q) と (i→q, j→p) の両対応 × ± 符号 4 通り
                    foreach (var (np, nq) in new[] { (p, q), (q, p) })
                        foreach (var si in new[] { 1.0, -1.0 })
                            foreach (var sj in new[] { 1.0, -1.0 })
                            {
                                var r = SolveWahba([(si * oi.NormalSample, catalog[np].Direction, oi.Weight), (sj * oj.NormalSample, catalog[nq].Direction, oj.Weight)]);
                                if (r == null) continue;
                                var cand = ScoreAndAssign(r, obs, catalog, waveLength, assignTol, predictedBuf, usedObsBuf, usedNodeBuf);
                                if (cand != null && cand.AssignedBands >= 3)
                                    candidates.Add(cand);
                            }
                }
            }
        #endregion

        #region 上位 seed の反復精緻化 (weighted Kabsch × 2) → 重複除去 → 上位返却
        var refinedList = new List<EbsdOrientationCandidate>();
        foreach (var cand in candidates.OrderByDescending(c => c.Score).Take(60))
        {
            cancel.ThrowIfCancellationRequested();
            var current = cand;
            for (int iter = 0; iter < 2; iter++)
            {
                var eqs = new List<(V3 m, V3 g, double w)>();
                foreach (var (bandIdx, hkl) in current.Assignments)
                {
                    var o = obs.First(x => x.BandIndex == bandIdx);
                    var node = catalog.First(n => n.Indices.Contains(hkl));
                    //観測法線の符号を現在回転に合わせる
                    var predicted = current.Rotation * node.Direction; //260724Cl (/simplify): 手書き MultiplySample → 既存の Matrix3D×Vector3d 演算子
                    var m = V3.Dot(o.NormalSample, predicted) >= 0 ? o.NormalSample : -o.NormalSample;
                    eqs.Add((m, node.Direction, o.Weight));
                }
                if (eqs.Count < 2) break;
                var r = SolveWahba([.. eqs]);
                if (r == null) break;
                var next = ScoreAndAssign(r, obs, catalog, waveLength, assignTol, predictedBuf, usedObsBuf, usedNodeBuf);
                if (next == null || next.Score <= current.Score) break;
                current = next;
            }
            //重複除去 (misorientation < 3°)
            if (!refinedList.Any(c => MisorientationDeg(c.Rotation, current.Rotation) < 3))
                refinedList.Add(current);
        }
        return [.. refinedList.OrderByDescending(c => c.Score).Take(maxCandidates)];
        #endregion
    }

    /// <summary>回転 R (crystal→sample) の下で全バンドをカタログ方向へ一対一割当し、スコアを付ける。
    /// 260724Cl: 残差昇順の greedy 一対一 (複数バンドが同一理論面へ重複割当されて AssignedBands が水増しされるのを防ぐ。Codex 指摘)
    /// 260724Cl (/simplify) シグネチャ変更: predicted/usedObs/usedNode を呼び出し元バッファ受け取りに (旧: 呼び出し毎に new。数千回×catalog 長の GC churn)</summary>
    static EbsdOrientationCandidate ScoreAndAssign(Matrix3D rot, List<BandObservation> obs, List<PlaneNode> catalog, double waveLength, double assignTol,
        V3[] predicted, bool[] usedObs, bool[] usedNode)
    {
        //角度残差 < assignTol の全 (バンド, ノード) ペアを列挙
        //260724Cl (/simplify): cos は単調なので dot 比較でゲートし、Acos は通過ペア (通常数個) のみに遅延 (旧: obs×catalog 全ペアで Acos → 指数付け時間の主因)
        double cosAssignTol = Math.Cos(assignTol);
        for (int ni = 0; ni < catalog.Count; ni++) predicted[ni] = rot * catalog[ni].Direction; //既存の Matrix3D×Vector3d 演算子
        var pairs = new List<(int ObsIdx, int NodeIdx, double Ang)>();
        for (int oi = 0; oi < obs.Count; oi++)
            for (int ni = 0; ni < catalog.Count; ni++)
            {
                double d = Math.Abs(V3.Dot(obs[oi].NormalSample, predicted[ni]));
                if (d > cosAssignTol) pairs.Add((oi, ni, Math.Acos(Math.Min(d, 1))));
            }

        var cand = new EbsdOrientationCandidate { Rotation = rot, TotalBands = obs.Count };
        Array.Clear(usedObs); Array.Clear(usedNode);
        double score = 0, sqSum = 0; int n = 0;
        foreach (var (oi, ni, ang) in pairs.OrderBy(p => p.Ang))
        {
            if (usedObs[oi] || usedNode[ni]) continue;
            usedObs[oi] = true; usedNode[ni] = true;
            var o = obs[oi]; var node = catalog[ni];

            //d (バンド幅) の soft likelihood: 調和反射のうち最も近い d を採用
            double widthFactor = 1;
            var hkl = node.Indices[0];
            if (!double.IsNaN(o.SinBragg))
            {
                double bestR = double.MaxValue;
                for (int m = 0; m < node.DValues.Length; m++)
                {
                    double r = Math.Abs(o.SinBragg - waveLength / (2 * node.DValues[m])) / o.SinBraggSigma;
                    if (r < bestR) { bestR = r; hkl = node.Indices[m]; }
                }
                widthFactor = Math.Exp(-0.5 * Math.Min(bestR, 4) * Math.Min(bestR, 4)) * 0.6 + 0.4; //幅不一致でも 0.4 は残す (角度主体)
            }

            double angFactor = 1 - ang / assignTol; //0..1
            score += o.Weight * angFactor * widthFactor;
            sqSum += ang * ang; n++;
            cand.Assignments[o.BandIndex] = hkl;
        }
        if (n == 0) return null;
        cand.AssignedBands = n;
        cand.AngularRmsDeg = Math.Sqrt(sqSum / n) * 180 / Math.PI;
        cand.Score = score;
        return cand;
    }

    #region 数学ユーティリティ

    /// <summary>Wahba/Kabsch: Σ w·m·gᵀ の SVD から R (g→m) を求める。m=試料系観測、g=結晶系。縮退時は null</summary>
    static Matrix3D SolveWahba((V3 m, V3 g, double w)[] eqs)
    {
        var h = MathNet.Numerics.LinearAlgebra.Double.DenseMatrix.Create(3, 3, 0);
        foreach (var (m, g, w) in eqs)
        {
            h[0, 0] += w * m.X * g.X; h[0, 1] += w * m.X * g.Y; h[0, 2] += w * m.X * g.Z;
            h[1, 0] += w * m.Y * g.X; h[1, 1] += w * m.Y * g.Y; h[1, 2] += w * m.Y * g.Z;
            h[2, 0] += w * m.Z * g.X; h[2, 1] += w * m.Z * g.Y; h[2, 2] += w * m.Z * g.Z;
        }
        var svd = h.Svd(true);
        //2 対応では H の rank=2 だが SVD は直交補空間を埋めるので回転は一意に決まる (符号は det 補正)
        if (svd.S[1] < 1E-10) return null; //rank<2 (平行ベクトル) は不定
        var u = svd.U; var vt = svd.VT;
        double det = (u * vt).Determinant();
        var d = MathNet.Numerics.LinearAlgebra.Double.DenseMatrix.CreateDiagonal(3, 3, i => i == 2 ? det : 1);
        var r = u * d * vt;
        //Matrix3D の 9 引数コンストラクタは column-major (e11,e21,e31, e12,...) なので明示代入で行優先を避ける
        return new Matrix3D(
            r[0, 0], r[1, 0], r[2, 0],
            r[0, 1], r[1, 1], r[2, 1],
            r[0, 2], r[1, 2], r[2, 2]);
    }

    //260724Cl (/simplify): MultiplySample (行列×ベクトルの手書き展開) は既存の Matrix3D×Vector3d 演算子 (Matrix.cs、FMA 使用) の再実装だったため削除

    /// <summary>2 つの回転行列間の misorientation 角 (deg、対称考慮なし)</summary>
    public static double MisorientationDeg(Matrix3D r1, Matrix3D r2)
    {
        //trace(R1ᵀ·R2)
        double tr = r1.E11 * r2.E11 + r1.E21 * r2.E21 + r1.E31 * r2.E31
                  + r1.E12 * r2.E12 + r1.E22 * r2.E22 + r1.E32 * r2.E32
                  + r1.E13 * r2.E13 + r1.E23 * r2.E23 + r1.E33 * r2.E33;
        return Math.Acos(Math.Clamp((tr - 1) / 2, -1, 1)) * 180 / Math.PI;
    }

    static int LowerBound(double[] arr, double v)
    {
        int lo = 0, hi = arr.Length;
        while (lo < hi) { int mid = (lo + hi) / 2; if (arr[mid] < v) lo = mid + 1; else hi = mid; }
        return lo;
    }

    static int UpperBound(double[] arr, double v)
    {
        int lo = 0, hi = arr.Length;
        while (lo < hi) { int mid = (lo + hi) / 2; if (arr[mid] <= v) lo = mid + 1; else hi = mid; }
        return lo;
    }

    #endregion
}
