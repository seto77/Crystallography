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
    //260724Cl シグネチャ変更 (thoroughCoarse 追加): true で粗段を 96px 完全 RobustPreprocess の総当たりにする (数倍遅いが判別力最大。
    //48px 軽量前処理の粗段が正解盆地を落とす画像 (コントラスト弱め) への対策。作者指示「パワープレーで良い」)。
    //旧: (..., int refineKeep = 12, CancellationToken cancel = default)
    //260725Cl シグネチャ変更 (properSymmetries 追加): 結晶の proper 対称回転 (単位元は含めない)。指定時は R·S 同値集合の
    //fundamental-zone 代表のみ粗段評価する (探索 ~1/(1+個数)。monoclinic C2 で総当たり半減 — Codex 裁定 260725)
    //旧: (..., bool thoroughCoarse = false, CancellationToken cancel = default)
    public static List<EbsdOrientationCandidate> Index(
        MasterPattern mp, float[] posPlane, float[] negPlane, EbsdDetectorGeometry geometry,
        double[] expValues, int expWidth, int expHeight,
        double coarseStepDeg = 5, int maxCandidates = 10, int coarseKeep = 64, int refineKeep = 12,
        bool thoroughCoarse = false,
        Matrix3D[] properSymmetries = null,
        System.Threading.CancellationToken cancel = default)
    {
        //実測参照を 2 解像度で準備 (粗段 = 軽量前処理 48px または完全 robust 96px / 精密段 = 完全 robust 96px)。
        //260724Cl 高速化: thorough では全段 RobustPreprocessFast (box3 近似・scratch 再利用・逐次 — 入れ子 Parallel 競合と GC 圧を解消)。
        //両側同一の原則に従い参照側も Fast パイプで生成する
        double[] refCoarse; int cw, ch;
        if (thoroughCoarse)
        {
            var (dRef, w96, h96) = EbsdPatternScorer.Downsample(expValues, expWidth, expHeight, 96);
            refCoarse = new double[w96 * h96]; cw = w96; ch = h96;
            EbsdPatternScorer.RobustPreprocessFast(dRef, w96, h96, refCoarse, new double[w96 * h96], new double[w96 * h96]);
        }
        else
            (refCoarse, cw, ch) = PrepareLight(expValues, expWidth, expHeight, 48);
        double[] refFine; int fw, fh;
        if (thoroughCoarse) { refFine = refCoarse; fw = cw; fh = ch; } //thorough は全段 96px Fast で参照共有
        else (refFine, fw, fh) = EbsdPatternScorer.PrepareReferenceRobust(expValues, expWidth, expHeight, 96);
        var projCoarse = new EbsdPatternProjector(geometry, cw, ch);
        var projFine = new EbsdPatternProjector(geometry, fw, fh);

        #region ①粗段: Fibonacci 球 × 面内回転の総当たり (方位単位並列・バッファ再利用)
        int nSphere = Math.Max(64, (int)(4 * Math.PI / (coarseStepDeg * Math.PI / 180 * (coarseStepDeg * Math.PI / 180))));
        int nPhi = Math.Max(16, (int)(360 / coarseStepDeg));
        double golden = Math.PI * (3 - Math.Sqrt(5));

        //260725Cl 変更: 旧 GridRotation(di,pi) を球点部 (SphereRotation) と面内部 (rzTable) に分離 — pi 毎の r0 再計算と
        //Rot(z,φ) 再計算を除去 (SphereRotation(di) * rzTable[pi] は旧式と同一演算・ビット一致)。
        //旧: Matrix3D GridRotation(int di, int pi) { ...(r0 計算)... return r0 * Matrix3D.Rot(new V3(0, 0, 1), pi * 2 * Math.PI / nPhi); }
        var rzTable = new Matrix3D[nPhi];
        for (int pi0 = 0; pi0 < nPhi; pi0++) rzTable[pi0] = Matrix3D.Rot(new V3(0, 0, 1), pi0 * 2 * Math.PI / nPhi);

        Matrix3D SphereRotation(int di)
        {
            double z = 1 - 2.0 * (di + 0.5) / nSphere;
            double rxy = Math.Sqrt(Math.Max(0, 1 - z * z));
            double az = di * golden;
            var nHat = new V3(rxy * Math.Cos(az), rxy * Math.Sin(az), z);
            var axis = V3.Cross(V3.UnitZ, nHat);
            return axis.Length < 1E-9
                ? (nHat.Z > 0 ? Matrix3D.IdentityMatrix : Matrix3D.Rot(new V3(1, 0, 0), Math.PI))
                : Matrix3D.Rot(axis.Normalized(), Math.Acos(Math.Clamp(nHat.Z, -1, 1)));
        }

        //260725Cl 追加: proper 対称回転による fundamental-zone 判定 — 同値集合 {R·S} の trace 最大代表のみ評価する
        //(最小回転角代表。monoclinic C2b では E11+E33>=0 と等価)。trace 同値の境界では両側 true (重複評価する安全側 —
        //Fibonacci 格子は R·S を厳密には含まないため片側破棄は被覆穴の危険。Codex 裁定 260725)
        bool IsFundamental(Matrix3D r)
        {
            if (properSymmetries == null) return true;
            double t = r.E11 + r.E22 + r.E33;
            foreach (var s in properSymmetries)
            {
                var rs = r * s;
                if (rs.E11 + rs.E22 + rs.E33 > t + 1E-12) return false;
            }
            return true;
        }

        int keepPerThread = Math.Max(64, coarseKeep * 4);
        var survivors = new List<(double S, int Di, int Pi)>();
        var lockObj = new object();
        //260725Cl: square 格子は面内回転分解プロジェクション (球点毎に Lambert 極座標を 1 回計算、面内 120 回は sector 折り返しのみ —
        //3×3 積・sqrt・atan 全除去。Codex 裁定 260725)。hex 格子は従来 Project へフォールバック
        bool inPlaneFast = mp.GridType != MasterPattern.Types.Hexagonal;
        System.Threading.Tasks.Parallel.For(0, nSphere,
            () => (Local: new List<(double S, int Di, int Pi)>(), Buf: new double[cw * ch], Dst: new double[cw * ch], T1: new double[cw * ch], T2: new double[cw * ch],
                   Th0: inPlaneFast ? new double[cw * ch] : null, Ra: inPlaneFast ? new double[cw * ch] : null, Rb: inPlaneFast ? new double[cw * ch] : null, Neg: inPlaneFast ? new bool[cw * ch] : null), //260725Cl
            (di, _, state) =>
            {
                cancel.ThrowIfCancellationRequested();
                var r0 = SphereRotation(di); //260725Cl: 球点回転は pi ループ外で 1 回 (旧: GridRotation が pi 毎に再計算)
                bool prepared = false; //260725Cl: PrepareSpherePoint は FZ を通る pi が現れた時のみ (全スキップ球点では省略)
                for (int pi = 0; pi < nPhi; pi++)
                {
                    var rot = r0 * rzTable[pi];
                    if (!IsFundamental(rot)) continue; //260725Cl: 対称同値の非代表側は Project 前にスキップ (C2 で総当たり半減)
                    //projCoarse.Project(mp, GridRotation(di, pi), posPlane, negPlane, state.Buf, parallel: false); //260725Cl 変更前
                    if (inPlaneFast) //260725Cl
                    {
                        if (!prepared) { projCoarse.PrepareSpherePoint(r0, state.Th0, state.Ra, state.Rb, state.Neg); prepared = true; }
                        projCoarse.ProjectInPlane(mp, pi * 2 * Math.PI / nPhi, posPlane, negPlane, state.Th0, state.Ra, state.Rb, state.Neg, state.Buf);
                    }
                    else
                        projCoarse.Project(mp, rot, posPlane, negPlane, state.Buf, parallel: false);
                    if (thoroughCoarse) //260724Cl: 完全 robust 前処理の総当たり (Fast 版 = box3 近似+バッファ再利用+逐次)
                    {
                        EbsdPatternScorer.RobustPreprocessFast(state.Buf, cw, ch, state.Dst, state.T1, state.T2);
                        state.Local.Add((EbsdPatternScorer.Zncc(refCoarse, state.Dst), di, pi));
                    }
                    else
                    {
                        ApplyLight(state.Buf, cw, ch);
                        state.Local.Add((EbsdPatternScorer.Zncc(refCoarse, state.Buf), di, pi));
                    }
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
            //var r = GridRotation(s.Di, s.Pi); //260725Cl 変更前
            var r = SphereRotation(s.Di) * rzTable[s.Pi]; //260725Cl: 旧 GridRotation と同一演算 (survivors は少数なので r0 再計算で可)
            if (basins.All(b => EbsdIndexer.MisorientationDeg(b.R, r) > coarseStepDeg * 1.5))
                basins.Add((s.S, r));
            if (basins.Count >= coarseKeep) break;
        }
        #endregion

        #region ②中段: 96px・完全 RobustPreprocess で再スコア (thorough 時は Fast 版)
        double ScoreFine(double[] buf, double[][] sc) //260724Cl: sc = [dst, t1, t2] (thorough 用 scratch)
        {
            if (thoroughCoarse)
            {
                EbsdPatternScorer.RobustPreprocessFast(buf, fw, fh, sc[0], sc[1], sc[2]);
                return EbsdPatternScorer.Zncc(refFine, sc[0]);
            }
            return EbsdPatternScorer.Zncc(refFine, EbsdPatternScorer.RobustPreprocess(buf, fw, fh));
        }
        var rescored = new (double S, Matrix3D R)[basins.Count];
        System.Threading.Tasks.Parallel.For(0, basins.Count,
            () => (Buf: new double[fw * fh], Sc: new[] { new double[fw * fh], new double[fw * fh], new double[fw * fh] }),
            (bi, _, st) =>
            {
                cancel.ThrowIfCancellationRequested();
                projFine.Project(mp, basins[bi].R, posPlane, negPlane, st.Buf, parallel: false);
                rescored[bi] = (ScoreFine(st.Buf, st.Sc), basins[bi].R);
                return st;
            },
            _ => { });
        var top = rescored.OrderByDescending(t => t.S).Take(refineKeep).ToList();
        #endregion

        #region ③精段: NelderMead 精密化 → misor NMS → 候補構築
        var refined = new (double S, Matrix3D R)[top.Count];
        System.Threading.Tasks.Parallel.For(0, top.Count,
            () => (Buf: new double[fw * fh], Sc: new[] { new double[fw * fh], new double[fw * fh], new double[fw * fh] }),
            (ti, _, st) =>
            {
                cancel.ThrowIfCancellationRequested();
                var r0 = top[ti].R;
                double Obj(double[] v)
                {
                    projFine.Project(mp, Perturb(r0, v[0], v[1], v[2]), posPlane, negPlane, st.Buf, parallel: false);
                    return -ScoreFine(st.Buf, st.Sc);
                }
                var (b1, _, _) = EbsdPatternScorer.NelderMead(Obj, [0, 0, 0], [coarseStepDeg * 0.5, coarseStepDeg * 0.5, coarseStepDeg * 0.5], 120);
                var (b2, v2, _) = EbsdPatternScorer.NelderMead(Obj, b1, [0.5, 0.5, 0.5], 80);
                refined[ti] = (-v2, Perturb(r0, b2[0], b2[1], b2[2]));
                return st;
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

    /// <summary>結晶の点群 proper 回転 (単位元を除く、crystal Cartesian 系) を返す。Index の properSymmetries 用。
    /// 対称操作の線形部 W (格子座標系、det=+1 のみ) を MatrixReal·W·MatrixReal⁻¹ で Cartesian 化し重複除去する。
    /// 対称が無い (P1 等) 場合は null。260725Cl 追加 (Codex 裁定 260725: 対称削減はグリッド生成後の FZ 判定で)</summary>
    public static Matrix3D[] GetProperRotations(Crystal crystal)
    {
        var ops = TSubgroupFinder.GetExpandedOps(crystal.SymmetrySeriesNumber);
        var a = crystal.MatrixReal;
        var ai = a.Inverse();
        var result = new List<Matrix3D>();
        foreach (var op in ops)
        {
            var w = SeitzNotation.LinearMatrix(op);
            int det = w[0, 0] * (w[1, 1] * w[2, 2] - w[1, 2] * w[2, 1]) - w[0, 1] * (w[1, 0] * w[2, 2] - w[1, 2] * w[2, 0]) + w[0, 2] * (w[1, 0] * w[2, 1] - w[1, 1] * w[2, 0]);
            if (det != 1) continue; //improper (回反・鏡映) は除外 — 検出器像はキラリティを区別する
            //Matrix3D 9 引数 ctor は column-major (第 1 列, 第 2 列, 第 3 列)
            var r = a * new Matrix3D(w[0, 0], w[1, 0], w[2, 0], w[0, 1], w[1, 1], w[2, 1], w[0, 2], w[1, 2], w[2, 2]) * ai;
            if (AbsDiff(r, Matrix3D.IdentityMatrix) < 1E-9) continue; //単位元
            if (result.Any(x => AbsDiff(x, r) < 1E-9)) continue;      //中心化等による重複
            result.Add(r);
        }
        return result.Count == 0 ? null : [.. result];

        static double AbsDiff(in Matrix3D x, in Matrix3D y)
            => Math.Abs(x.E11 - y.E11) + Math.Abs(x.E12 - y.E12) + Math.Abs(x.E13 - y.E13)
             + Math.Abs(x.E21 - y.E21) + Math.Abs(x.E22 - y.E22) + Math.Abs(x.E23 - y.E23)
             + Math.Abs(x.E31 - y.E31) + Math.Abs(x.E32 - y.E32) + Math.Abs(x.E33 - y.E33);
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
