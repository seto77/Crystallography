#region using
using System;
using System.Collections.Generic;
using System.Linq;
using V3 = OpenTK.Mathematics.Vector3d;
#endregion

namespace Crystallography;

/// <summary>
/// Primary indexing: MasterPattern から方位空間を網羅する辞書パターンを on-the-fly 生成し、実測との ZNCC 総当たりで方位候補を得る。260724Cl 追加。
/// 3 段構成 (Codex 裁定 260724):
///   ①粗段 = coarseStepDeg (既定 5°) の Fibonacci 球 × 面内回転を 48px・軽量前処理 (log → box DoG 5/21 → tanh → 正規化、両側同一) で総当たり
///     → 方位距離 NMS で上位 coarseKeep 盆地を抽出
///   ②中段 = 上位を 96px・完全 RobustPreprocess で再スコア → 上位 refineKeep
///   ③精段 = NelderMead 精密化 (±coarseStep/2 → 0.5°) → misor NMS → 上位 maxCandidates
/// 辞書は検出器幾何 (PC/DD) に依存するため事前保存せず毎回生成する (実行 ~数秒)。
/// 返り値の Score は NaN — Radon z の付与と複合ランクへの接続は呼び出し側の統合層 (FormEBSD.Indexing) が行う。
/// </summary>
public static class EbsdDictionaryIndexer
{
    /// <summary>
    /// 方位候補を辞書総当たりで探索する。expValues = 実測生強度 (expWidth×expHeight、前処理は内部で解像度別に整合させる)。
    /// posPlane/negPlane = MasterPattern.GetPlane の単一スライス。
    /// </summary>
    public static List<EbsdOrientationCandidate> Index(
        MasterPattern mp, float[] posPlane, float[] negPlane, EbsdDetectorGeometry geometry,
        double[] expValues, int expWidth, int expHeight,
        double coarseStepDeg = 5, int maxCandidates = 10, int coarseKeep = 64, int refineKeep = 12,
        System.Threading.CancellationToken cancel = default)
    {
        //実測参照を 2 解像度で準備 (粗段 = 軽量前処理 48px / 精密段 = 完全 robust 96px)
        var (refCoarse, cw, ch) = PrepareLight(expValues, expWidth, expHeight, 48);
        var (refFine, fw, fh) = EbsdPatternScorer.PrepareReferenceRobust(expValues, expWidth, expHeight, 96);
        var projCoarse = new EbsdPatternProjector(geometry, cw, ch);
        var projFine = new EbsdPatternProjector(geometry, fw, fh);

        #region ①粗段: Fibonacci 球 × 面内回転の総当たり (方位単位並列・バッファ再利用)
        int nSphere = Math.Max(64, (int)(4 * Math.PI / (coarseStepDeg * Math.PI / 180 * (coarseStepDeg * Math.PI / 180))));
        int nPhi = Math.Max(16, (int)(360 / coarseStepDeg));
        double golden = Math.PI * (3 - Math.Sqrt(5));

        Matrix3D GridRotation(int di, int pi)
        {
            double z = 1 - 2.0 * (di + 0.5) / nSphere;
            double rxy = Math.Sqrt(Math.Max(0, 1 - z * z));
            double az = di * golden;
            var nHat = new V3(rxy * Math.Cos(az), rxy * Math.Sin(az), z);
            var axis = V3.Cross(V3.UnitZ, nHat);
            var r0 = axis.Length < 1E-9
                ? (nHat.Z > 0 ? Matrix3D.IdentityMatrix : Matrix3D.Rot(new V3(1, 0, 0), Math.PI))
                : Matrix3D.Rot(axis.Normalized(), Math.Acos(Math.Clamp(nHat.Z, -1, 1)));
            return r0 * Matrix3D.Rot(new V3(0, 0, 1), pi * 2 * Math.PI / nPhi);
        }

        int keepPerThread = Math.Max(64, coarseKeep * 4);
        var survivors = new List<(double S, int Di, int Pi)>();
        var lockObj = new object();
        System.Threading.Tasks.Parallel.For(0, nSphere,
            () => (Local: new List<(double S, int Di, int Pi)>(), Buf: new double[cw * ch]),
            (di, _, state) =>
            {
                cancel.ThrowIfCancellationRequested();
                for (int pi = 0; pi < nPhi; pi++)
                {
                    projCoarse.Project(mp, GridRotation(di, pi), posPlane, negPlane, state.Buf, parallel: false);
                    ApplyLight(state.Buf, cw, ch);
                    state.Local.Add((EbsdPatternScorer.Zncc(refCoarse, state.Buf), di, pi));
                }
                if (state.Local.Count > keepPerThread * 4)
                {
                    state.Local.Sort((a, b) => b.S.CompareTo(a.S));
                    state.Local.RemoveRange(keepPerThread, state.Local.Count - keepPerThread);
                }
                return state;
            },
            state => { lock (lockObj) survivors.AddRange(state.Local); });
        survivors.Sort((a, b) => b.S.CompareTo(a.S));

        //方位距離 NMS: 上位から順に、既採用と近い (misor < 1.5×粗刻み) 方位を捨てて盆地の代表だけ残す
        var basins = new List<(double S, Matrix3D R)>();
        foreach (var s in survivors)
        {
            var r = GridRotation(s.Di, s.Pi);
            if (basins.All(b => EbsdIndexer.MisorientationDeg(b.R, r) > coarseStepDeg * 1.5))
                basins.Add((s.S, r));
            if (basins.Count >= coarseKeep) break;
        }
        #endregion

        #region ②中段: 96px・完全 RobustPreprocess で再スコア
        var rescored = new (double S, Matrix3D R)[basins.Count];
        System.Threading.Tasks.Parallel.For(0, basins.Count,
            () => new double[fw * fh],
            (bi, _, buf) =>
            {
                cancel.ThrowIfCancellationRequested();
                projFine.Project(mp, basins[bi].R, posPlane, negPlane, buf, parallel: false);
                rescored[bi] = (EbsdPatternScorer.Zncc(refFine, EbsdPatternScorer.RobustPreprocess(buf, fw, fh)), basins[bi].R);
                return buf;
            },
            _ => { });
        var top = rescored.OrderByDescending(t => t.S).Take(refineKeep).ToList();
        #endregion

        #region ③精段: NelderMead 精密化 → misor NMS → 候補構築
        var refined = new (double S, Matrix3D R)[top.Count];
        System.Threading.Tasks.Parallel.For(0, top.Count,
            () => new double[fw * fh],
            (ti, _, buf) =>
            {
                cancel.ThrowIfCancellationRequested();
                var r0 = top[ti].R;
                double Obj(double[] v)
                {
                    projFine.Project(mp, Perturb(r0, v[0], v[1], v[2]), posPlane, negPlane, buf, parallel: false);
                    return -EbsdPatternScorer.Zncc(refFine, EbsdPatternScorer.RobustPreprocess(buf, fw, fh));
                }
                var (b1, _, _) = EbsdPatternScorer.NelderMead(Obj, [0, 0, 0], [coarseStepDeg * 0.5, coarseStepDeg * 0.5, coarseStepDeg * 0.5], 120);
                var (b2, v2, _) = EbsdPatternScorer.NelderMead(Obj, b1, [0.5, 0.5, 0.5], 80);
                refined[ti] = (-v2, Perturb(r0, b2[0], b2[1], b2[2]));
                return buf;
            },
            _ => { });

        var result = new List<EbsdOrientationCandidate>();
        foreach (var (s, r) in refined.OrderByDescending(t => t.S))
        {
            if (result.Any(c => EbsdIndexer.MisorientationDeg(c.Rotation, r) < 2)) continue;
            result.Add(new EbsdOrientationCandidate { Rotation = r, Score = double.NaN, Zncc = s, AngularRmsDeg = double.NaN });
            if (result.Count >= maxCandidates) break;
        }
        return result;
        #endregion
    }

    /// <summary>粗段の軽量前処理 (Codex 裁定: 背景除算を省き log → box DoG (5×5−21×21) → tanh(z/3) → 正規化)。in place。
    /// log は平均で正規化した log-ratio (floor も平均比例) — 実測 (輝度 0-255) とシミュレーション (強度 ~1E-3) のレンジ差で
    /// 前処理が非対称にならないようにする (絶対値 floor は不可)。260724Cl</summary>
    static void ApplyLight(double[] v, int w, int h)
    {
        double mean = 0;
        foreach (var x in v) mean += x;
        mean = Math.Max(mean / v.Length, 1E-30);
        double floor = mean * 0.01;
        for (int i = 0; i < v.Length; i++) v[i] = Math.Log(Math.Max(v[i], floor) / mean);
        var b1 = BoxBlur(v, w, h, 2);  //5×5
        var b2 = BoxBlur(v, w, h, 10); //21×21
        for (int i = 0; i < v.Length; i++) v[i] = b1[i] - b2[i];
        EbsdPatternScorer.NormalizeInPlace(v);
        for (int i = 0; i < v.Length; i++) v[i] = Math.Tanh(v[i] / 3);
        EbsdPatternScorer.NormalizeInPlace(v);
    }

    /// <summary>実測を box 縮小して軽量前処理を掛けた参照を返す (粗段用、1 回きり)</summary>
    static (double[] Data, int W, int H) PrepareLight(double[] values, int width, int height, int targetLongSide)
    {
        var (dst, w, h) = EbsdPatternScorer.Downsample(values, width, height, targetLongSide);
        ApplyLight(dst, w, h);
        return (dst, w, h);
    }

    /// <summary>running box 平均 (半径 radius、境界は有効画素数で正規化)。分離 2 パス O(N)</summary>
    static double[] BoxBlur(double[] src, int w, int h, int radius)
    {
        var tmp = new double[w * h];
        for (int y = 0; y < h; y++) //横パス
        {
            int row = y * w;
            double sum = 0; int n = 0;
            for (int x = 0; x <= Math.Min(radius, w - 1); x++) { sum += src[row + x]; n++; }
            for (int x = 0; x < w; x++)
            {
                tmp[row + x] = sum / n;
                int add = x + radius + 1, rem = x - radius;
                if (add < w) { sum += src[row + add]; n++; }
                if (rem >= 0) { sum -= src[row + rem]; n--; }
            }
        }
        var dst = new double[w * h];
        for (int x = 0; x < w; x++) //縦パス
        {
            double sum = 0; int n = 0;
            for (int y = 0; y <= Math.Min(radius, h - 1); y++) { sum += tmp[y * w + x]; n++; }
            for (int y = 0; y < h; y++)
            {
                dst[y * w + x] = sum / n;
                int add = y + radius + 1, rem = y - radius;
                if (add < h) { sum += tmp[add * w + x]; n++; }
                if (rem >= 0) { sum -= tmp[rem * w + x]; n--; }
            }
        }
        return dst;
    }

    /// <summary>方位摂動 R'=Rot(ω̂,|ω|)·R0 (deg)。EbsdRadonIndexer.Perturb と同一規約</summary>
    static Matrix3D Perturb(Matrix3D r0, double wxDeg, double wyDeg, double wzDeg)
    {
        double wx = wxDeg * Math.PI / 180, wy = wyDeg * Math.PI / 180, wz = wzDeg * Math.PI / 180;
        double len = Math.Sqrt(wx * wx + wy * wy + wz * wz);
        return len < 1E-12 ? r0 : Matrix3D.Rot(new V3(wx / len, wy / len, wz / len), len) * r0;
    }
}
