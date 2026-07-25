#region using
using System;
using System.Collections.Generic;
using System.Linq;
using V3 = OpenTK.Mathematics.Vector3d;
#endregion

namespace Crystallography;

/// <summary>
/// 実測 EBSD パターンの Radon/butterfly 応答マップ (方位テンプレート照合の証拠)。260724Cl 追加。
/// Abs = |butterfly 平滑応答| (運動学では excess/deficiency の符号を予測できないため絶対値)。
/// 座標系: θ = 中心線法線の角度 (0.5° 刻み、[0,180) 循環で θ+180 は ρ 反転)、ρ = work 画像中心からの符号付き距離 (px)。
/// </summary>
public sealed class EbsdRadonMap
{
    public double[] Abs;
    public int NTheta, NRho, RhoOffset;
    public double ThetaStepDeg;
    public int WorkW, WorkH;
    /// <summary>native px → work px の縮小率</summary>
    public double Scale;
    /// <summary>Abs の robust null 統計 (median と 1.4826×MAD)。スコアの SNR 正規化に使う</summary>
    public double Mu0, Sigma0;

    double[] dilated; //粗探索用の異方膨張マップ (遅延構築)

    /// <summary>θ を [0,180) に折り返す (θ+180 は ρ 反転)。戻り値 = (θdeg∈[0,180), ρ)</summary>
    static (double thetaDeg, double rho) Fold(double thetaDeg, double rho)
    {
        thetaDeg %= 360;
        if (thetaDeg < 0) thetaDeg += 360;
        if (thetaDeg >= 180) { thetaDeg -= 180; rho = -rho; }
        return (thetaDeg, rho);
    }

    /// <summary>証拠マップの bilinear 補間値。マップ外の ρ は 0 (=無証拠)</summary>
    public double Sample(double thetaDeg, double rhoWork)
    {
        (thetaDeg, rhoWork) = Fold(thetaDeg, rhoWork);
        double tb = thetaDeg / ThetaStepDeg, rb = rhoWork + RhoOffset;
        int t0 = (int)Math.Floor(tb);
        double ft = tb - t0;
        if (t0 >= NTheta) { t0 = 0; } //数値誤差ガード (θ≈180)
        int t1 = t0 + 1;
        bool wrap = t1 >= NTheta; //θ=180 の継ぎ目: 行 0 を ρ 反転で参照
        if (wrap) t1 = 0;

        double RowSample(int t, double r)
        {
            int r0 = (int)Math.Floor(r);
            if (r0 < 0 || r0 >= NRho - 1) return 0;
            double fr = r - r0;
            return Abs[t * NRho + r0] * (1 - fr) + Abs[t * NRho + r0 + 1] * fr;
        }
        double v0 = RowSample(t0, rb);
        double v1 = RowSample(t1, wrap ? NRho - 1 - rb : rb);
        return v0 * (1 - ft) + v1 * ft;
    }

    /// <summary>
    /// 粗探索用の異方膨張マップを構築する (θ±thetaBins 行 × ρ±rhoPx 列の窓 max)。
    /// 粗方位グリッドの離散化誤差で予測線が真のピークから外れても盆地を取りこぼさないため (Codex 裁定: θ 方向の膨張も必須)
    /// </summary>
    public void BuildDilated(int thetaBins, int rhoPx)
    {
        //① ρ 方向の窓 max (行内)
        var tmp = new double[Abs.Length];
        System.Threading.Tasks.Parallel.For(0, NTheta, t =>
        {
            for (int r = 0; r < NRho; r++)
            {
                double m = 0;
                int lo = Math.Max(0, r - rhoPx), hi = Math.Min(NRho - 1, r + rhoPx);
                for (int j = lo; j <= hi; j++) m = Math.Max(m, Abs[t * NRho + j]);
                tmp[t * NRho + r] = m;
            }
        });
        //② θ 方向の窓 max (循環: 継ぎ目は ρ 反転)
        var dst = new double[Abs.Length];
        System.Threading.Tasks.Parallel.For(0, NTheta, t =>
        {
            for (int r = 0; r < NRho; r++)
            {
                double m = 0;
                for (int dt = -thetaBins; dt <= thetaBins; dt++)
                {
                    int tt = t + dt; int rr = r;
                    if (tt < 0) { tt += NTheta; rr = NRho - 1 - r; }
                    else if (tt >= NTheta) { tt -= NTheta; rr = NRho - 1 - r; }
                    m = Math.Max(m, tmp[tt * NRho + rr]);
                }
                dst[t * NRho + r] = m;
            }
        });
        dilated = dst;
    }

    /// <summary>膨張マップの最近傍値 (粗探索用。bilinear 不要)</summary>
    public double SampleDilatedNearest(double thetaDeg, double rhoWork)
    {
        (thetaDeg, rhoWork) = Fold(thetaDeg, rhoWork);
        int t = (int)Math.Round(thetaDeg / ThetaStepDeg);
        if (t >= NTheta) t = 0;
        int r = (int)Math.Round(rhoWork) + RhoOffset;
        return (uint)r < (uint)NRho ? dilated[t * NRho + r] : 0;
    }
}

/// <summary>
/// Radon 証拠マップに対する運動学的テンプレート照合で結晶方位候補を直接探索する。260724Cl 追加。GUI 非依存。
/// バンドの離散検出 (検出/非検出の二値判断) を経由せず、SO(3) を階層探索して
/// Score(R) = Σ_k w_k·(E(θ_k,ρ_k) − μ₀) / √(Σ_k w_k²·σ₀²) (SNR 正規化、Codex 裁定) を最大化する。
/// k は回転 R の下で検出器と交差する予測 Kikuchi バンド中心線、w_k = √(運動学的強度) (上限クリップ)。
/// </summary>
public static class EbsdRadonIndexer
{
    #region カタログ (方向ノード)
    sealed class Node
    {
        public V3 Dir;         //結晶直交系の単位ベクトル (±は正規化)
        public double Weight;  //√強度 (クリップ済)
        public (int H, int K, int L) Hkl; //代表 (最強) 反射
    }

    /// <summary>±g・同方向調和反射を 1 ノードに統合し、強度上位 maxNodes に制限。
    /// weightExponent=0.5 (既定): w=√I を median の 3 倍でクリップ (Codex 裁定 Q5)。
    /// それ以外 (例 0.25): w=clip((I/median I)^exp, 0.5, 2) — 強度重みの弱化 (実マップ採点頑健化、Codex 裁定 260724)。260724Cl シグネチャ変更 (weightExponent 追加)</summary>
    static List<Node> BuildCatalog(IEnumerable<Vector3D> kikuchiReflections, int maxNodes, double weightExponent = 0.5)
    {
        var nodes = new List<(V3 Dir, double I, (int, int, int) Hkl)>();
        foreach (var g in kikuchiReflections)
        {
            var len = Math.Sqrt(g.X * g.X + g.Y * g.Y + g.Z * g.Z);
            if (len < 1E-12) continue;
            var dir = new V3(g.X / len, g.Y / len, g.Z / len);
            if (dir.Z < 0 || (dir.Z == 0 && dir.Y < 0) || (dir.Z == 0 && dir.Y == 0 && dir.X < 0)) dir = -dir;

            int found = nodes.FindIndex(n => V3.Dot(n.Dir, dir) > 0.99995); //~0.6°
            if (found < 0) nodes.Add((dir, g.RelativeIntensity, (g.Index.h, g.Index.k, g.Index.l)));
            else if (g.RelativeIntensity > nodes[found].I) nodes[found] = (nodes[found].Dir, g.RelativeIntensity, (g.Index.h, g.Index.k, g.Index.l));
        }
        var top = nodes.OrderByDescending(n => n.I).Take(maxNodes).ToList();
        if (Math.Abs(weightExponent - 0.5) > 1E-9) //260724Cl: 弱化重み (I/medI)^exp を [0.5, 2] クリップ
        {
            double medI = Math.Max(top.Select(n => n.I).OrderBy(v => v).ElementAt(top.Count / 2), 1E-12);
            return [.. top.Select(n => new Node { Dir = n.Dir, Weight = Math.Clamp(Math.Pow(Math.Max(n.I, 0) / medI, weightExponent), 0.5, 2), Hkl = n.Hkl })];
        }
        var ws = top.Select(n => Math.Sqrt(Math.Max(n.I, 0))).ToArray();
        double wClip = 3 * ws.OrderBy(v => v).ElementAt(ws.Length / 2);
        return [.. top.Select((n, i) => new Node { Dir = n.Dir, Weight = Math.Min(ws[i], Math.Max(wClip, 1E-12)), Hkl = n.Hkl })];
    }
    #endregion

    /// <summary>単一方位の Radon 証拠 z スコア (Index の厳密スコアと同一評価)。ZNCC 精密化のガードや複合ランクの再評価用。260724Cl 追加</summary>
    public static double ScoreOrientation(EbsdRadonMap map, EbsdDetectorGeometry geometry, IEnumerable<Vector3D> kikuchiReflections,
        Matrix3D rotation, double saturateCap = 0, double weightExponent = 0.5, int maxNodes = 90)
        => ScoreExactCore(map, geometry, BuildCatalog(kikuchiReflections, maxNodes, weightExponent), rotation, saturateCap);

    /// <summary>厳密スコア本体 (bilinear・全ノード・非膨張マップ・近接予測線の排他込み)。
    /// Index 内クロージャからインライン実装を抽出 (公開 ScoreOrientation と共用)。260724Cl 追加</summary>
    static double ScoreExactCore(EbsdRadonMap map, EbsdDetectorGeometry geometry, List<Node> catalog, Matrix3D rot, double saturateCap)
    {
        double xm = geometry.XMirror, pix = geometry.PixelSize;
        var ey = geometry.Ey; var center = geometry.Center;
        double mu0 = map.Mu0, sigma0 = map.Sigma0;
        double rhoLimit = map.RhoOffset - 2;
        var exTh = new double[catalog.Count]; var exRho = new double[catalog.Count];
        double num = 0, wSum = 0, w2Sum = 0;
        int nAcc = 0;
        foreach (var node in catalog)
        {
            var gs = rot * node.Dir;
            var gl = geometry.SampleToLab(gs);
            double aMm = xm * gl.X, bMm = gl.Y * ey.Y + gl.Z * ey.Z, cMm = gl.X * center.X + gl.Y * center.Y + gl.Z * center.Z;
            double norm = Math.Sqrt(aMm * aMm + bMm * bMm);
            if (norm < 1E-9) continue;
            double rhoWork = -cMm / (pix * norm) * map.Scale;
            if (Math.Abs(rhoWork) > rhoLimit) continue;
            var (thF, rhoF) = FoldLine(Math.Atan2(bMm, aMm) * 180 / Math.PI, rhoWork);
            bool dup = false;
            for (int j = 0; j < nAcc; j++)
                if (SameLine(thF, rhoF, exTh[j], exRho[j], 2, 5)) { dup = true; break; }
            if (dup) continue; //二重得点防止 (強度降順 → 先着優先)
            exTh[nAcc] = thF; exRho[nAcc] = rhoF; nAcc++;
            double e = map.Sample(thF, rhoF) - mu0;
            //証拠飽和 (粗探索と同一の ψ)。cap=0 で旧動作
            if (saturateCap > 0)
            {
                double zk = e / sigma0;
                num += node.Weight * (zk > 0 ? saturateCap * Math.Tanh(zk / saturateCap) : Math.Max(zk, -1));
            }
            else
                num += node.Weight * Math.Max(e, -sigma0);
            wSum += node.Weight; w2Sum += node.Weight * node.Weight;
        }
        if (w2Sum <= 0 || wSum * wSum / w2Sum < 4) return double.MinValue;
        return num / (Math.Sqrt(w2Sum) * (saturateCap > 0 ? 1 : sigma0)); //飽和時は num が z 単位
    }

    /// <summary>2 つの (θ[deg]∈[0,180), ρ) 線が実質同一か (θ 循環と ρ 反転を考慮)。260724Cl 追加。
    /// 同一方位評価内で近接する複数の予測線 (調和・近縁反射) が同じ観測リッジから二重に得点するのを防ぐ
    /// (旧 pair-angle 指数付けの greedy 1対1 割当に相当する排他。これがないと二重得点方位が正解を上回る)</summary>
    static bool SameLine(double t1, double r1, double t2, double r2, double tolDeg, double tolRho)
    {
        double d = Math.Abs(t1 - t2);
        if (d < tolDeg && Math.Abs(r1 - r2) < tolRho) return true;
        return d > 180 - tolDeg && Math.Abs(r1 + r2) < tolRho;
    }

    /// <summary>θ を [0,180) に折り返す (θ+180 は ρ 反転)。260724Cl 追加</summary>
    static (double thetaDeg, double rho) FoldLine(double thetaDeg, double rho)
    {
        thetaDeg %= 360;
        if (thetaDeg < 0) thetaDeg += 360;
        if (thetaDeg >= 180) { thetaDeg -= 180; rho = -rho; }
        return (thetaDeg, rho);
    }

    /// <summary>回転 rot の下で検出器と交差する予測バンド中心線 (native px 係数 A·col+B·row+C=0、正規化済) を返す。検証・診断用。260724Cl 追加</summary>
    public static List<(double A, double B, double C, (int H, int K, int L) Hkl, double Weight)> PredictLines(
        Matrix3D rot, EbsdDetectorGeometry geometry, IEnumerable<Vector3D> kikuchiReflections, int maxNodes = 90)
    {
        var lines = new List<(double, double, double, (int, int, int), double)>();
        foreach (var node in BuildCatalog(kikuchiReflections, maxNodes))
        {
            var gl = geometry.SampleToLab(rot * node.Dir);
            var (a, b, c) = geometry.LabNormalToLine(gl);
            double norm = Math.Sqrt(a * a + b * b);
            if (norm < 1E-12) continue;
            //中心線が画像矩形と交差するか (画像中心からの距離 < 半対角)
            double cx = geometry.WidthPx / 2.0 - 0.5, cy = geometry.HeightPx / 2.0 - 0.5;
            double rho = Math.Abs(a * cx + b * cy + c) / norm;
            if (rho > Math.Sqrt(cx * cx + cy * cy)) continue;
            lines.Add((a / norm, b / norm, c / norm, node.Hkl, node.Weight));
        }
        return lines;
    }

    /// <summary>証拠マップの局所極大 (μ₀+3σ₀ 超) を NMS 付きで抽出し、native 線係数で返す。260724Cl 追加。
    /// pair-angle シード用の内部ピーク (ユーザー向けバンド検出ではないため検証ゲート等は不要)</summary>
    static List<(double A, double B, double C, double Score)> ExtractPeakLines(EbsdRadonMap map, int topK)
    {
        double th0 = map.Mu0 + 3 * map.Sigma0;
        var peaks = new List<(int T, int R, double V)>();
        for (int t = 0; t < map.NTheta; t++)
            for (int r = 1; r < map.NRho - 1; r++)
            {
                double v = map.Abs[t * map.NRho + r];
                if (v > th0 && v >= map.Abs[t * map.NRho + r - 1] && v >= map.Abs[t * map.NRho + r + 1])
                    peaks.Add((t, r, v));
            }
        var picked = new List<(double Th, double Rho, double V)>();
        foreach (var (t, r, v) in peaks.OrderByDescending(p => p.V))
        {
            double th = t * map.ThetaStepDeg, rho = r - map.RhoOffset;
            if (picked.Any(p => SameLine(th, rho, p.Th, p.Rho, 3, 6))) continue;
            picked.Add((th, rho, v));
            if (picked.Count >= topK) break;
        }
        //(θ, ρ_work) → native 線係数 A·col+B·row+C=0 (A=cosθ, B=sinθ, C=−(A·cx+B·cy+ρ_native))
        double cx = map.WorkW / map.Scale / 2.0 - 0.5, cy = map.WorkH / map.Scale / 2.0 - 0.5;
        return [.. picked.Select(p =>
        {
            var (sin, cos) = Math.SinCos(p.Th * Math.PI / 180);
            return (cos, sin, -(cos * cx + sin * cy + p.Rho / map.Scale), p.V);
        })];
    }

    /// <summary>
    /// 方位候補を探索する。kikuchiReflections = Crystal.VectorOfG_KikuchiLine、waveLength = 電子波長 (nm)。
    /// シードは 2 系統 (①証拠マップ内部ピークの pair-angle+Kabsch = サブ度精度、②SO(3) 粗グリッド = ピーク抽出漏れの保険)、
    /// 採点は密な Radon 証拠 z 値に一本化。返り値はスコア降順。
    /// 対称等価 (カタログ方向集合を保存する回転で結ばれる方位=パターン上区別不能) は 1 代表に縮約。
    /// saturateCap: 0 で旧動作 (floor のみ)。>0 で証拠 z の正側を ψ(z)=cap·tanh(z/cap) で飽和・負側を max(z,−1) — 少数強リッジの支配を抑える (推奨 4)。
    /// weightExponent: 0.5 (既定) = w=√I 3median クリップ。0.25 = w=clip((I/medI)^0.25, 0.5, 2) の弱化重み。(実マップ採点頑健化、Codex 裁定 260724)
    /// </summary>
    //260724Cl シグネチャ変更 (実験パラメータ saturateCap / weightExponent 追加)。旧:
    //public static List<EbsdOrientationCandidate> Index(EbsdRadonMap map, EbsdDetectorGeometry geometry,
    //    IEnumerable<Vector3D> kikuchiReflections, double waveLength = 0.00859,
    //    int maxCandidates = 10, double coarseStepDeg = 3, int maxNodes = 90, int coarseNodes = 20,
    //    System.Threading.CancellationToken cancel = default)
    public static List<EbsdOrientationCandidate> Index(EbsdRadonMap map, EbsdDetectorGeometry geometry,
        IEnumerable<Vector3D> kikuchiReflections, double waveLength = 0.00859,
        int maxCandidates = 10, double coarseStepDeg = 3, int maxNodes = 90, int coarseNodes = 20,
        double saturateCap = 0, double weightExponent = 0.5,
        System.Threading.CancellationToken cancel = default)
    {
        var refl = kikuchiReflections as IReadOnlyList<Vector3D> ?? [.. kikuchiReflections]; //260724Cl: 多重列挙防止
        var catalog = BuildCatalog(refl, maxNodes, weightExponent);
        if (catalog.Count < 4) return [];
        var coarseCatalog = catalog.Take(coarseNodes).ToList(); //強度上位のみで粗探索 (Codex 裁定 Q6)

        //膨張窓: 粗グリッド半刻み (coarseStep/2) の回転で予測線が動く最大量をカバーする。
        //ρ 移動 ≈ (DD + |PC からの距離|)·δ ≲ 対角長·δ/2 → work px 換算。θ 移動 ≈ δ
        double halfStepRad = coarseStepDeg / 2 * Math.PI / 180;
        int dilRho = Math.Max(3, (int)Math.Ceiling(Math.Sqrt(map.WorkW * map.WorkW + map.WorkH * map.WorkH) * 0.7 * halfStepRad));
        int dilTheta = Math.Max(2, (int)Math.Ceiling(coarseStepDeg / 2 / map.ThetaStepDeg));
        map.BuildDilated(dilTheta, dilRho);

        #region 粗探索: Fibonacci 球面 × 面内回転のグリッドを膨張マップで採点
        int nSphere = Math.Max(64, (int)(4 * Math.PI / (coarseStepDeg * Math.PI / 180 * (coarseStepDeg * Math.PI / 180))));
        int nPhi = Math.Max(16, (int)(360 / coarseStepDeg));
        double golden = Math.PI * (3 - Math.Sqrt(5));

        //ノード方向を配列化 (SoA)
        int nc = coarseCatalog.Count;
        var ndx = new double[nc]; var ndy = new double[nc]; var ndz = new double[nc]; var nw = new double[nc];
        for (int i = 0; i < nc; i++) { ndx[i] = coarseCatalog[i].Dir.X; ndy[i] = coarseCatalog[i].Dir.Y; ndz[i] = coarseCatalog[i].Dir.Z; nw[i] = coarseCatalog[i].Weight; }

        double mu0 = map.Mu0, sigma0 = map.Sigma0;
        double xm = geometry.XMirror, pix = geometry.PixelSize;
        var ey = geometry.Ey; var center = geometry.Center;
        double sinSmp = Math.Sin(geometry.SampleTilt), cosSmp = Math.Cos(geometry.SampleTilt);
        double rhoLimit = map.RhoOffset - 2;

        var survivors = new List<(double S, int Di, int Pi)>();
        var lockObj = new object();
        const int keepPerMerge = 600;

        System.Threading.Tasks.Parallel.For(0, nSphere,
            () => new List<(double S, int Di, int Pi)>(),
            (di, _, local) =>
            {
                cancel.ThrowIfCancellationRequested();
                //Fibonacci 球面点 n̂ = R·ẑ (結晶 ẑ 軸の行き先)
                double z = 1 - 2.0 * (di + 0.5) / nSphere;
                double rxy = Math.Sqrt(Math.Max(0, 1 - z * z));
                double az = di * golden;
                double nxs = rxy * Math.Cos(az), nys = rxy * Math.Sin(az), nzs = z;
                //R0: ẑ→n̂ の最短回転 (Rodrigues)。R = R0·Rz(φ)
                double r00, r01, r02, r10, r11, r12, r20, r21, r22;
                {
                    double ax = -nys, ay = nxs; //軸 = ẑ×n̂ (正規化前)、sinA=|axis|, cosA=nz
                    double sinA = Math.Sqrt(ax * ax + ay * ay);
                    if (sinA < 1E-9)
                    {
                        if (nzs > 0) { r00 = 1; r01 = 0; r02 = 0; r10 = 0; r11 = 1; r12 = 0; r20 = 0; r21 = 0; r22 = 1; }
                        else { r00 = 1; r01 = 0; r02 = 0; r10 = 0; r11 = -1; r12 = 0; r20 = 0; r21 = 0; r22 = -1; } //x 軸 180°
                    }
                    else
                    {
                        double ux = ax / sinA, uy = ay / sinA, c = nzs, s = sinA, t = 1 - c;
                        r00 = c + ux * ux * t; r01 = ux * uy * t; r02 = uy * s;
                        r10 = ux * uy * t; r11 = c + uy * uy * t; r12 = -ux * s;
                        r20 = -uy * s; r21 = ux * s; r22 = c;
                    }
                }
                var accTh = new double[nc]; var accRho = new double[nc]; //260724Cl: 同一評価内の予測線排他バッファ (強度降順に採用)
                for (int pi = 0; pi < nPhi; pi++)
                {
                    var (sinP, cosP) = Math.SinCos(pi * 2 * Math.PI / nPhi);
                    double num = 0, den = 0, wSum = 0, w2Sum = 0;
                    int nAcc = 0;
                    for (int k = 0; k < nc; k++)
                    {
                        //g_s = R0·Rz(φ)·d
                        double dx = ndx[k] * cosP - ndy[k] * sinP, dy = ndx[k] * sinP + ndy[k] * cosP, dz = ndz[k];
                        double gx = r00 * dx + r01 * dy + r02 * dz;
                        double gy = r10 * dx + r11 * dy + r12 * dz;
                        double gz = r20 * dx + r21 * dy + r22 * dz;
                        //sample → lab (Rx(SmpTilt) の逆 = SampleToLab)
                        double lx = gx, ly = cosSmp * gy - sinSmp * gz, lz = sinSmp * gy + cosSmp * gz;
                        //lab 法線 → 検出器線 (θ, ρ_work)
                        double aMm = xm * lx, bMm = ly * ey.Y + lz * ey.Z, cMm = lx * center.X + ly * center.Y + lz * center.Z;
                        double norm = Math.Sqrt(aMm * aMm + bMm * bMm);
                        if (norm < 1E-9) continue; //バンド面がほぼ視軸に垂直 (線が無限遠)
                        double rhoWork = -cMm / (pix * norm) * map.Scale;
                        if (Math.Abs(rhoWork) > rhoLimit) continue; //検出器と交差しない
                        //260724Cl: 近接予測線の排他 (カタログは強度降順 → 先着=強い方を採用)。二重得点防止
                        var (thF, rhoF) = FoldLine(Math.Atan2(bMm, aMm) * 180 / Math.PI, rhoWork);
                        bool dup = false;
                        for (int j = 0; j < nAcc; j++)
                            if (SameLine(thF, rhoF, accTh[j], accRho[j], 2, 5)) { dup = true; break; }
                        if (dup) continue;
                        accTh[nAcc] = thF; accRho[nAcc] = rhoF; nAcc++;
                        double e = map.SampleDilatedNearest(thF, rhoF) - mu0;
                        double w = nw[k];
                        //260724Cl: 証拠飽和 (正側 ψ(z)=cap·tanh(z/cap)、負側 max(z,−1)) — 少数強リッジの支配抑制 (Codex 裁定 260724)。cap=0 で旧動作
                        if (saturateCap > 0)
                        {
                            double zk = e / sigma0;
                            num += w * (zk > 0 ? saturateCap * Math.Tanh(zk / saturateCap) : Math.Max(zk, -1));
                        }
                        else
                            num += w * Math.Max(e, -sigma0); //視野隅の希薄線の過剰ペナルティを floor
                        wSum += w; w2Sum += w * w;
                    }
                    if (w2Sum <= 0) continue;
                    double nEff = wSum * wSum / w2Sum;
                    if (nEff < 4) continue; //有効バンド数下限 (少数バンド方位の上振れ防止、Codex 裁定 Q1)
                    den = Math.Sqrt(w2Sum) * (saturateCap > 0 ? 1 : sigma0); //260724Cl: 飽和時は num が z 単位 (σ₀ 正規化済)。旧: den = Math.Sqrt(w2Sum) * sigma0
                    double score = num / den;
                    local.Add((score, di, pi));
                }
                if (local.Count > keepPerMerge * 4)
                {
                    local.Sort((a, b) => b.S.CompareTo(a.S));
                    local.RemoveRange(keepPerMerge, local.Count - keepPerMerge);
                }
                return local;
            },
            local => { lock (lockObj) survivors.AddRange(local); });

        survivors.Sort((a, b) => b.S.CompareTo(a.S));
        if (survivors.Count > keepPerMerge) survivors.RemoveRange(keepPerMerge, survivors.Count - keepPerMerge);
        #endregion

        #region 粗 seed の重複除去 → NelderMead 局所精密化 (非膨張マップ・全ノード)
        //260725Cl (/simplify): 球点部を EbsdIndexer.FibonacciSphereRotation へ統合 (EbsdDictionaryIndexer と重複していた。同一演算・ビット一致)
        //旧: Matrix3D SeedRotation(int di, int pi) { ...(z/rxy/az/nHat/axis から r0 生成)... return r0 * Matrix3D.Rot(new V3(0, 0, 1), pi * 2 * Math.PI / nPhi); }
        Matrix3D SeedRotation(int di, int pi)
            => EbsdIndexer.FibonacciSphereRotation(di, nSphere) * Matrix3D.Rot(new V3(0, 0, 1), pi * 2 * Math.PI / nPhi);

        var seeds = new List<(double S, Matrix3D R)>();

        //260724Cl: ① pair-angle シード — 証拠マップの内部ピーク線を擬似バンドとして旧 pair-angle+Kabsch 指数付けに掛ける。
        //グリッド+NM だけではリッジ (θ 幅 ~1°) への到達精度が不足し、真の方位が z 最適値まで到達できないことを
        //診断で確認 (4-1_33: z(真)=21.2 > 探索トップ 17.3 = カバレッジ欠落、5-2_22: 残差 ~1° の減点で真が負ける)。
        //Kabsch は割当ペアの角度残差を直接最小化するためサブ度精度のシードが得られる
        var peakLines = ExtractPeakLines(map, 24);
        if (peakLines.Count >= 3)
        {
            var pseudo = peakLines.Select(p => new EbsdBand
            {
                LineA = p.A, LineB = p.B, LineC = p.C,
                CenterAnchors = [(geometry.WidthPx / 2.0, geometry.HeightPx / 2.0)],
                EdgePoints = [], CenterQuality = 1, WidthQuality = 0,
            }).ToList();
            foreach (var c in EbsdIndexer.Index(pseudo, geometry, refl, waveLength, maxCandidates: 12, cancel: cancel))
                seeds.Add((double.MaxValue, c.Rotation));
        }

        //② SO(3) 粗グリッドの生存者 (ピーク抽出漏れ・擬似バンド不足時の保険)
        foreach (var s in survivors)
        {
            var r = SeedRotation(s.Di, s.Pi);
            if (seeds.All(x => EbsdIndexer.MisorientationDeg(x.R, r) > coarseStepDeg * 0.8))
                seeds.Add((s.S, r));
            if (seeds.Count >= 52) break;
        }

        //厳密スコア (bilinear・全ノード・非膨張・近接予測線の排他込み)。Parallel から呼ばれるため排他バッファはローカル確保
        //260724Cl: 本体を ScoreExactCore へ抽出 (公開 ScoreOrientation と共用のため。旧インライン実装は ScoreExactCore に移動)
        double ScoreExact(Matrix3D rot) => ScoreExactCore(map, geometry, catalog, rot, saturateCap);

        static Matrix3D Perturb(Matrix3D r0, double wxDeg, double wyDeg, double wzDeg)
        {
            double wx = wxDeg * Math.PI / 180, wy = wyDeg * Math.PI / 180, wz = wzDeg * Math.PI / 180;
            double len = Math.Sqrt(wx * wx + wy * wy + wz * wz);
            return len < 1E-12 ? r0 : Matrix3D.Rot(new V3(wx / len, wy / len, wz / len), len) * r0;
        }

        var refined = new (double S, Matrix3D R)[seeds.Count];
        System.Threading.Tasks.Parallel.For(0, seeds.Count, si =>
        {
            cancel.ThrowIfCancellationRequested();
            var r0 = seeds[si].R;
            double Obj(double[] v) => -ScoreExact(Perturb(r0, v[0], v[1], v[2]));
            var (b1, _, _) = EbsdPatternScorer.NelderMead(Obj, [0, 0, 0], [coarseStepDeg * 0.5, coarseStepDeg * 0.5, coarseStepDeg * 0.5], 120);
            var (b2, v2, _) = EbsdPatternScorer.NelderMead(Obj, b1, [0.4, 0.4, 0.4], 80);
            refined[si] = (-v2, Perturb(r0, b2[0], b2[1], b2[2]));
        });
        #endregion

        #region 対称等価の縮約 → 候補構築
        //等価判定: Q = R1ᵀ·R2 がカタログ方向集合を (± 込みで) 保存するなら、2 方位の予測バンド線集合は同一 = パターン上区別不能
        bool Equivalent(Matrix3D r1, Matrix3D r2)
        {
            if (EbsdIndexer.MisorientationDeg(r1, r2) < 2) return true;
            var q = new Matrix3D(r1.E11, r1.E12, r1.E13, r1.E21, r1.E22, r1.E23, r1.E31, r1.E32, r1.E33); //R1ᵀ (column-major ctor に行順で渡す = 転置)
            foreach (var node in catalog)
            {
                var mapped = q * (r2 * node.Dir); //R1ᵀ·R2·d
                bool hit = false;
                foreach (var other in catalog)
                    if (Math.Abs(V3.Dot(mapped, other.Dir)) > 0.99996) { hit = true; break; } //~0.5°
                if (!hit) return false;
            }
            return true;
        }

        var result = new List<EbsdOrientationCandidate>();
        foreach (var (s, r) in refined.OrderByDescending(x => x.S))
        {
            cancel.ThrowIfCancellationRequested();
            if (s == double.MinValue || result.Any(c => Equivalent(c.Rotation, r))) continue;

            //強い証拠を持つノード (z ≥ 3) を情報として列挙 (スコアと同じ近接排他を適用)
            var cand = new EbsdOrientationCandidate { Rotation = r, Score = s, AngularRmsDeg = double.NaN };
            int inView = 0, nAcc2 = 0;
            var acc2Th = new double[catalog.Count]; var acc2Rho = new double[catalog.Count];
            var strong = new List<(double Z, (int H, int K, int L) Hkl)>();
            foreach (var node in catalog)
            {
                var gl = geometry.SampleToLab(r * node.Dir);
                double aMm = xm * gl.X, bMm = gl.Y * ey.Y + gl.Z * ey.Z, cMm = gl.X * center.X + gl.Y * center.Y + gl.Z * center.Z;
                double norm = Math.Sqrt(aMm * aMm + bMm * bMm);
                if (norm < 1E-9) continue;
                double rhoWork = -cMm / (pix * norm) * map.Scale;
                if (Math.Abs(rhoWork) > rhoLimit) continue;
                var (thF, rhoF) = FoldLine(Math.Atan2(bMm, aMm) * 180 / Math.PI, rhoWork);
                bool dup = false;
                for (int j = 0; j < nAcc2; j++)
                    if (SameLine(thF, rhoF, acc2Th[j], acc2Rho[j], 2, 5)) { dup = true; break; }
                if (dup) continue;
                acc2Th[nAcc2] = thF; acc2Rho[nAcc2] = rhoF; nAcc2++;
                inView++;
                double zVal = (map.Sample(thF, rhoF) - mu0) / sigma0;
                if (zVal >= 3) strong.Add((zVal, node.Hkl));
            }
            cand.TotalBands = inView;
            cand.AssignedBands = strong.Count;
            foreach (var (x, i) in strong.OrderByDescending(x => x.Z).Take(12).Select((x, i) => (x, i)))
                cand.Assignments[i] = x.Hkl;
            result.Add(cand);
            if (result.Count >= maxCandidates) break;
        }
        return result;
        #endregion
    }
}
