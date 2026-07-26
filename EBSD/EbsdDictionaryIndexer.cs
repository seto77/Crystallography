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
///   ①粗段 = coarseStepDeg の Fibonacci 球 × 面内回転を総当たり → 方位距離 NMS で上位 coarseKeep 盆地を抽出
///   ②中段 = 上位を 96px・完全 RobustPreprocess で再スコア → 上位 refineKeep
///   ③精段 = NelderMead 精密化 (±coarseStep/2 → 0.5°) → misor NMS → 上位 maxCandidates
/// 粗段の前処理は thoroughCoarse で切り替わる。**GUI は常に thoroughCoarse=true (96px 完全 RobustPreprocessFast、刻み 3°)** —
/// 既定の false (48px・軽量前処理 log → box DoG 5/21 → tanh → 正規化) はコントラストの弱い画像で正解盆地を落とすため実運用では使わず、
/// 検証ハーネスの対照条件としてのみ残している (260725Cl 追記: 出荷経路と既定値が食い違う点を明示)。
/// ⚠thoroughCoarse=true では refFine/fw/fh/projFine が①粗段と同一 (下の参照共有・projector 共有) なので、
/// ②中段は「別解像度での再スコア」ではなく **同一パイプでの再採点** になる (実質は面内分解 ProjectInPlane と
/// 通常 Project の投影経路差の再確認)。①②を別段の解像度として読むと誤解するので注記する。260725Cl
/// 辞書は検出器幾何 (PC/DD) に依存するため事前保存せず毎回生成する (実行 ~数秒)。
/// 返り値の Score は NaN — Radon z の付与と複合ランクへの接続は呼び出し側の統合層 (EbsdOrientationSearch) が行う。
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
        System.Threading.CancellationToken cancel = default,
        //260725Cl 追加: 粗段の進捗 (0-1) を通知する。null で従来どおり無通知。ワーカースレッドから呼ばれるので受け手側でマーシャリングすること
        Action<double> progress = null)
    {
        //260725Ch: ホットループへ入る前に公開APIの寸法・plane境界・探索個数を一度だけ検証
        ArgumentNullException.ThrowIfNull(mp);
        ArgumentNullException.ThrowIfNull(geometry);
        ArgumentNullException.ThrowIfNull(expValues);
        if (expWidth <= 0) throw new ArgumentOutOfRangeException(nameof(expWidth));
        if (expHeight <= 0) throw new ArgumentOutOfRangeException(nameof(expHeight));
        if (expValues.Length != checked(expWidth * expHeight)) throw new ArgumentException("expValues.Length must equal expWidth * expHeight.", nameof(expValues));
        if (!(coarseStepDeg > 0) || !double.IsFinite(coarseStepDeg)) throw new ArgumentOutOfRangeException(nameof(coarseStepDeg));
        if (maxCandidates <= 0) throw new ArgumentOutOfRangeException(nameof(maxCandidates));
        if (coarseKeep <= 0) throw new ArgumentOutOfRangeException(nameof(coarseKeep));
        if (refineKeep <= 0) throw new ArgumentOutOfRangeException(nameof(refineKeep));
        if (mp.GridSize < 2) throw new ArgumentException("MasterPattern.GridSize must be at least 2.", nameof(mp));
        int requiredPlaneLength = checked(mp.GridSize * mp.GridSize);
        if (posPlane == null && negPlane == null) throw new ArgumentException("At least one master-pattern hemisphere is required.");
        if (posPlane != null && posPlane.Length < requiredPlaneLength) throw new ArgumentException("The positive master-pattern plane is too short.", nameof(posPlane));
        if (negPlane != null && negPlane.Length < requiredPlaneLength) throw new ArgumentException("The negative master-pattern plane is too short.", nameof(negPlane));

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
        //260725Cl (/simplify): thorough では粗段と精密段が同一寸法 (上の refFine 共有と同条件) なので projector も共有する
        //(旧: 常に new。96px の ray キャッシュ (Vector3d × w·h) を二重に構築・保持していた)
        var projFine = (fw == cw && fh == ch) ? projCoarse : new EbsdPatternProjector(geometry, fw, fh);

        #region ①粗段: Fibonacci 球 × 面内回転の総当たり (方位単位並列・バッファ再利用)
        int nSphere = Math.Max(64, (int)(4 * Math.PI / (coarseStepDeg * Math.PI / 180 * (coarseStepDeg * Math.PI / 180))));
        int nPhi = Math.Max(16, (int)(360 / coarseStepDeg));

        //260725Cl 変更: 旧 GridRotation(di,pi) を球点部 (EbsdIndexer.FibonacciSphereRotation — Radon 側 SeedRotation と共通化) と
        //面内部 (rzTable) に分離 — pi 毎の r0 再計算と Rot(z,φ) 再計算を除去 (積は旧式と同一演算・ビット一致)。
        //旧: Matrix3D GridRotation(int di, int pi) { ...(r0 計算)... return r0 * Matrix3D.Rot(new V3(0, 0, 1), pi * 2 * Math.PI / nPhi); }
        var rzTable = new Matrix3D[nPhi];
        for (int pi0 = 0; pi0 < nPhi; pi0++) rzTable[pi0] = Matrix3D.Rot(new V3(0, 0, 1), pi0 * 2 * Math.PI / nPhi);

        int keepPerThread = Math.Max(64, coarseKeep * 4);
        var survivors = new List<(double S, int Di, int Pi)>();
        var lockObj = new object();
        int coarseDone = 0; //260725Cl: 進捗通知用 (粗段が総時間の大半)
        var parallelOptions = new System.Threading.Tasks.ParallelOptions { CancellationToken = cancel }; //260725Ch: OCEをAggregateException化せず、TPL全ワーカーへ中止を伝播
        //260725Cl: square 格子は面内回転分解プロジェクション (球点毎に Lambert 極座標を 1 回計算、面内 120 回は sector 折り返しのみ —
        //3×3 積・sqrt・atan 全除去。Codex 裁定 260725)。hex 格子は従来 Project へフォールバック
        bool inPlaneFast = mp.GridType != MasterPattern.Types.Hexagonal;
        int symCount = properSymmetries?.Length ?? 0;
        //System.Threading.Tasks.Parallel.For(0, nSphere, //260725Ch 変更前
        System.Threading.Tasks.Parallel.For(0, nSphere, parallelOptions, //260725Ch
            () => new CoarseScratch(cw * ch, inPlaneFast, symCount), //260725Cl (/simplify): 旧 9 要素 ValueTuple の名前付き化
            (di, _, st) =>
            {
                cancel.ThrowIfCancellationRequested();
                var r0 = EbsdIndexer.FibonacciSphereRotation(di, nSphere); //260725Cl: 球点回転は pi ループ外で 1 回
                //260725Cl (/simplify): fundamental-zone 判定を trace 閉形式化 — Rz は z 回転 (E13=E23=E31=E32=0, E33=1) なので
                //trace(Rz·M) = Rz.E11·M.E11 + Rz.E12·M.E21 + Rz.E21·M.E12 + Rz.E22·M.E22 + M.E33。R=r0·Rz の trace は循環置換で M=r0、
                //R·S の trace は trace(Rz·S·r0) で M=S·r0 (di 毎に 1 回だけ行列積)。同値集合 {R·S} の trace 最大代表のみ評価する
                //(最小回転角代表。monoclinic C2b では E11+E33>=0 と等価)。
                //⚠260725Cl 訂正: 旧コメントは「±1E-12 の tie 帯が被覆穴の保険」としていたが、これは実態と合わない。
                //trace = 1+2cos θ なので粗刻み 3° の半刻みで trace は最大 2·sin θ·0.026 ≈ 0.05 動く = 境界のどちら側かは
                //1E-12 より 10 桁大きい誤差で決まる。1E-12 は厳密な対称不変点しか救わない。実際に被覆を保っているのは
                //「軌道の相手 R·S の近傍にも別の格子点があり、そちらが代表として残る」機構で、粗段の方位誤差が最悪 1 刻み分
                //余分に乗り得る (NM の初期シンプレックス ±1.5° が通常は吸収)。これが EbsdOrientationSearch で FZ 除外を
                //proper 回転 1 個 (C2) に限定している理由の数値的裏付け。恒久対策は tie 帯を格子由来の変動幅へ広げること (未実施 —
                //C2 は現行 1E-12 で候補一致を A/B 検証済みなので、検証済み構成を崩さない)。
                //旧: pi 毎に rot=r0*rzTable[pi] と IsFundamental 内の r*s を Matrix3D 生成 (~55万×2 個の Gen0 圧) → 全廃
                double r0A = r0.E11, r0B = r0.E21, r0C = r0.E12, r0D = r0.E22, r0E = r0.E33;
                for (int si = 0; si < symCount; si++)
                {
                    var m = properSymmetries[si] * r0;
                    st.SymTr[si * 5] = m.E11; st.SymTr[si * 5 + 1] = m.E21; st.SymTr[si * 5 + 2] = m.E12; st.SymTr[si * 5 + 3] = m.E22; st.SymTr[si * 5 + 4] = m.E33;
                }
                bool prepared = false; //260725Cl: PrepareSpherePoint は FZ を通る pi が現れた時のみ (全スキップ球点では省略)
                for (int pi = 0; pi < nPhi; pi++)
                {
                    if (symCount > 0) //260725Cl: 対称同値の非代表側は Project 前にスキップ (C2 で総当たり半減)
                    {
                        var rz = rzTable[pi];
                        double tR = rz.E11 * r0A + rz.E12 * r0B + rz.E21 * r0C + rz.E22 * r0D + r0E;
                        bool fundamental = true;
                        for (int si = 0; si < symCount && fundamental; si++)
                            if (rz.E11 * st.SymTr[si * 5] + rz.E12 * st.SymTr[si * 5 + 1] + rz.E21 * st.SymTr[si * 5 + 2] + rz.E22 * st.SymTr[si * 5 + 3] + st.SymTr[si * 5 + 4] > tR + 1E-12)
                                fundamental = false;
                        if (!fundamental) continue;
                    }
                    if (inPlaneFast) //260725Cl
                    {
                        if (!prepared) { projCoarse.PrepareSpherePoint(r0, st.Q0, st.Ra, st.Rb, st.Neg); prepared = true; }
                        projCoarse.ProjectInPlane(mp, pi * 2 * Math.PI / nPhi, posPlane, negPlane, st.Q0, st.Ra, st.Rb, st.Neg, st.Buf);
                    }
                    else
                        projCoarse.Project(mp, r0 * rzTable[pi], posPlane, negPlane, st.Buf, parallel: false);
                    if (thoroughCoarse) //260724Cl: 完全 robust 前処理の総当たり (Fast 版 = box3 近似+バッファ再利用+逐次)
                    {
                        EbsdPatternScorer.RobustPreprocessFast(st.Buf, cw, ch, st.Dst, st.T1, st.T2);
                        st.Local.Add((EbsdPatternScorer.Zncc(refCoarse, st.Dst), di, pi));
                    }
                    else
                    {
                        ApplyLight(st.Buf, cw, ch);
                        st.Local.Add((EbsdPatternScorer.Zncc(refCoarse, st.Buf), di, pi));
                    }
                }
                if (st.Local.Count > keepPerThread * 4)
                {
                    st.Local.Sort((a, b) => b.S.CompareTo(a.S));
                    st.Local.RemoveRange(keepPerThread, st.Local.Count - keepPerThread);
                }
                //260725Cl 追加: 進捗通知 (64 球点ごと。粗段が総時間の大半なので 0-0.9 を割り当て、残りは再スコア+NM)
                if (progress != null && (System.Threading.Interlocked.Increment(ref coarseDone) & 63) == 0)
                    progress(0.9 * coarseDone / nSphere);
                return st;
            },
            st => { lock (lockObj) survivors.AddRange(st.Local); });
        survivors.Sort((a, b) => b.S.CompareTo(a.S));

        //方位距離 NMS: 上位から順に、既採用と近い (misor < 1.5×粗刻み) 方位を捨てて盆地の代表だけ残す
        var basins = new List<(double S, Matrix3D R)>();
        foreach (var s in survivors)
        {
            //var r = GridRotation(s.Di, s.Pi); //260725Cl 変更前
            var r = EbsdIndexer.FibonacciSphereRotation(s.Di, nSphere) * rzTable[s.Pi]; //260725Cl: 旧 GridRotation と同一演算 (survivors は少数なので r0 再計算で可)
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
        //System.Threading.Tasks.Parallel.For(0, basins.Count, //260725Ch 変更前
        System.Threading.Tasks.Parallel.For(0, basins.Count, parallelOptions, //260725Ch
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
        progress?.Invoke(0.95); //260725Cl: 粗段 0.9 → 盆地再スコア完了
        #endregion

        #region ③精段: NelderMead 精密化 → misor NMS → 候補構築
        var refined = new (double S, Matrix3D R)[top.Count];
        //System.Threading.Tasks.Parallel.For(0, top.Count, //260725Ch 変更前
        System.Threading.Tasks.Parallel.For(0, top.Count, parallelOptions, //260725Ch
            () => (Buf: new double[fw * fh], Sc: new[] { new double[fw * fh], new double[fw * fh], new double[fw * fh] }),
            (ti, _, st) =>
            {
                cancel.ThrowIfCancellationRequested();
                var r0 = top[ti].R;
                double Obj(double[] v)
                {
                    projFine.Project(mp, EbsdIndexer.PerturbRotation(r0, v[0], v[1], v[2]), posPlane, negPlane, st.Buf, parallel: false);
                    return -ScoreFine(st.Buf, st.Sc);
                }
                var (b1, _, _) = EbsdPatternScorer.NelderMead(Obj, [0, 0, 0], [coarseStepDeg * 0.5, coarseStepDeg * 0.5, coarseStepDeg * 0.5], 120);
                var (b2, v2, _) = EbsdPatternScorer.NelderMead(Obj, b1, [0.5, 0.5, 0.5], 80);
                refined[ti] = (-v2, EbsdIndexer.PerturbRotation(r0, b2[0], b2[1], b2[2]));
                return st;
            },
            _ => { });

        progress?.Invoke(1.0); //260725Cl: NM 精密化まで完了 (以降の候補構築は瞬時)
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

    /// <summary>粗段 Parallel.For のスレッドローカル作業領域。260725Cl 追加 (/simplify: 旧 9 要素 ValueTuple の名前付き化)</summary>
    sealed class CoarseScratch
    {
        public readonly List<(double S, int Di, int Pi)> Local = [];
        public readonly double[] Buf, Dst, T1, T2;
        public readonly double[] Q0, Ra, Rb; //面内分解プロジェクション用 (square 格子のみ、hex では null)
        public readonly bool[] Neg;
        public readonly double[] SymTr;      //FZ 判定の trace 係数 (対称 1 つあたり 5 値、di 毎に更新)
        public CoarseScratch(int n, bool inPlane, int symCount)
        {
            Buf = new double[n]; Dst = new double[n]; T1 = new double[n]; T2 = new double[n];
            if (inPlane) { Q0 = new double[n]; Ra = new double[n]; Rb = new double[n]; Neg = new bool[n]; }
            SymTr = new double[symCount * 5];
        }
    }

    /// <summary>結晶の点群 proper 回転 (単位元を除く、crystal Cartesian 系) を返す。Index の properSymmetries 用。
    /// 対称操作の線形部 W (格子座標系、det=+1 のみ) を MatrixReal·W·MatrixReal⁻¹ で Cartesian 化し重複除去する。
    /// 対称が無い (P1 等) 場合は null。260725Cl 追加 (Codex 裁定 260725: 対称削減はグリッド生成後の FZ 判定で)</summary>
    public static Matrix3D[] GetProperRotations(Crystal crystal)
    {
        var ops = TSubgroupFinder.GetExpandedOps(crystal.SymmetrySeriesNumber);
        var a = crystal.MatrixReal;
        // var ai = a.Inverse(); //260725Cl 変更前: 同一値の再計算
        var ai = crystal.MatrixInverse; //260725Cl 変更 (/simplify): SetAxis がキャッシュ済みの実格子逆行列を使う
        var result = new List<Matrix3D>();
        foreach (var op in ops)
        {
            var w = SeitzNotation.LinearMatrix(op);
            //int det = w[0, 0] * (w[1, 1] * w[2, 2] - ...); //260725Cl 変更前 (/simplify): SymmetryProperties.Det (既存の int 3×3 行列式) を internal 昇格して共有
            if (SymmetryProperties.Det(w) != 1) continue; //improper (回反・鏡映) は除外 — 検出器像はキラリティを区別する
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

    //260725Cl (/simplify): Perturb は EbsdIndexer.PerturbRotation へ統合 (EbsdRadonIndexer・FormEBSD 側と 3 重複していた。式・演算順は同一)
    //旧: static Matrix3D Perturb(Matrix3D r0, double wxDeg, double wyDeg, double wzDeg) { ...(ω を rad 化し Rot(ω̂,|ω|)·r0)... }
}
