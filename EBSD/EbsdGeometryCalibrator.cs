using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;

namespace Crystallography;

/// <summary>
/// 幾何較正の結果。260726Cl 追加 (FormEBSD.cs の指数付け region — 統合前は FormEBSD.Indexing.cs — の匿名タプルを名前付きにしたもの)。 //260727Cl: 移設元ファイル名を実在するものへ訂正
/// PatternCenterU/V と CameraLength は <see cref="EbsdDetectorGeometry.FromPatternCenter"/> で DetX/DetY/DetZ へ戻せる。
/// BestIndex / Spread / NearBest / Flat* は「単一パターンで幾何がどこまで決まっているか」を利用側が表示するための診断値。
/// </summary>
/// <param name="Rotation">較正後の方位</param>
/// <param name="PatternCenterU">PC (垂線の足) の物理面内 mm 座標 (検出器中心基準)</param>
/// <param name="PatternCenterV">同上 (面内 V)</param>
/// <param name="CameraLength">検出器距離 DD (mm)</param>
/// <param name="Zncc">到達 ZNCC</param>
/// <param name="ZnccStart">較正前 ZNCC</param>
/// <param name="Evaluations">目的関数の総評価回数</param>
/// <param name="Rounds">最良解が使った交互最適化のラウンド数</param>
/// <param name="Converged">収束で打ち切ったか (false = 上限到達)</param>
/// <param name="JointGain">6 変数同時最適化が交互法の到達点から伸ばした量</param>
/// <param name="Starts">多点開始の点数</param>
/// <param name="BestIndex">最良だった開始点の番号</param>
/// <param name="Spread">最良と最悪の ZNCC 差</param>
/// <param name="NearBest">最良から 1E-3 以内に入った開始点の数 (= 最良解の basin の広さ)</param>
/// <param name="FlatU">その集団における PC (U) の広がり (半値幅、mm)</param>
/// <param name="FlatV">同上 (PC の V)</param>
/// <param name="FlatDd">同上 (DD)</param>
public sealed record EbsdCalibrationResult(
    Matrix3D Rotation,
    double PatternCenterU, double PatternCenterV,
    double CameraLength,
    double Zncc, double ZnccStart,
    int Evaluations,
    int Rounds, bool Converged,
    double JointGain,
    int Starts, int BestIndex, double Spread,
    int NearBest,
    double FlatU, double FlatV, double FlatDd);

/// <summary>
/// 検出器のパターンセンター (PC) と検出器距離 (DD) の較正 (方位も交互に微調整)。DetTilt は固定。
/// 260726Cl 追加: FormEBSD.cs の buttonCalibrateGeometry_Click 内 Task.Run 本体をそのまま移設したもの (GUI 非依存)。 //260727Cl: 旧 FormEBSD.Indexing.cs は同日 FormEBSD.cs へ統合済みなのでファイル名を訂正
/// 単一パターンでは DetTilt と方位 X 回転がゲージ自由度になるため Tilt は較正しない (設計正本 §7.2 / Codex 裁定)。
/// 最適化は (PC_u, PC_v, ln DD) と方位 3 変数の alternating fit → 方位仕上げ → 6 変数同時最適化を、多点開始で繰り返す。
/// </summary>
public static class EbsdGeometryCalibrator
{
    /// <summary>幾何較正の交互最適化 (方位 ⇄ PC/DD) の最大ラウンド数。260725Cl: 2 固定 → 10 → 20 (作者指示: 10 でも十分速い)。
    /// PC・DD・方位は単一パターンで強く相関しており、交互法は谷底でジグザグするため 2 ラウンドでは収束の保証が無かった。
    /// 実際には <see cref="ZnccTolerance"/> で早期終了するので、上限まで回るのは収束が遅い配置のときだけ</summary>
    public const int MaxRounds = 20;

    /// <summary>1 ラウンドの ZNCC 改善がこれ未満なら収束とみなして較正を打ち切る。260725Cl 追加</summary>
    const double ZnccTolerance = 1E-4;

    /// <summary>較正の多点開始の点数。260726Cl: 10 → 200 → 40 (作者指示)。
    /// 1 点あたり 0.2 秒程度なので全体で 8 秒前後。200 点で ±8% を探しても最良は現在の幾何のままだったので、日常はこの点数で足りる</summary>
    const int StartCount = 40;

    /// <summary>多点開始の振れ幅。PC は検出器幅・高さに対する割合、DD は lnDD の絶対値 (0.08 ≈ 8%)。260726Cl 追加。
    /// 当初の PC ±1%・lnDD ±0.02 では実機で 200 点すべてが同じ谷に落ち (best #0、200 within 1E-3、spread 0.0007)、
    /// 多点開始が機能していなかった。作者が観測した別の谷はもっと離れているので広げる。
    /// 較正のソフト境界 (<see cref="SoftBoundPcFraction"/> / <see cref="SoftBoundLnDd"/>) の内側に収めること</summary>
    const double StartSpreadPc = 0.08, StartSpreadLnDd = 0.08;

    /// <summary>較正のソフト境界。較正開始位置からの PC のずれ上限 (検出器 W/H に対する割合) と lnDD のずれ上限
    /// (0.35 ≈ DD ±40%)。単一パターンでは PC・DD・方位が縮退するため、これを超える解は非物理として
    /// <see cref="SoftBoundPenaltyBase"/> 起点の罰則値を返し、Nelder-Mead に採らせない。
    /// 260727Cl 追加: 交互法②と 6 変数同時仕上げに裸のリテラルで 2 重に書かれていたので命名した</summary>
    const double SoftBoundPcFraction = 0.25, SoftBoundLnDd = 0.35;

    /// <summary>ソフト境界の外で返す罰則値の下駄。目的関数は -ZNCC (= 高々 1 程度) なので、10 なら確実に棄却される。260727Cl 追加</summary>
    const double SoftBoundPenaltyBase = 10;

    //260920Cl (/simplify) 削除: const int JointPolishMaxEval = 600; — 多段化で StageMaxEval[].Joint へ吸収され、
    //  doc も「1 評価ごとに projector を作り直す」という撤去済みの設計を説明したままだった

    /// <summary>較正の多点開始オフセット。値は無次元 [-1,1]³ で、消費側で <see cref="StartSpreadPc"/> (検出器幅・高さ比) と
    /// <see cref="StartSpreadLnDd"/> (lnDD) を掛けてスケールする。260726Cl 追加 (作者要望)。 //260727Cl: doc が旧値 (1%・0.02) のままで実装 (8%・0.08) と食い違っていたので訂正
    /// 乱数を使わず決定的にする (同じ入力なら同じ結果)。[0] は現在の幾何そのもの、以降は Halton 列で [-1,1]³ を準一様に埋める。
    /// 局所解が多く (初期 DetX/Y/Z で最終スコアが 0.3 程度ばらつく)、同時最適化でも壁は越えられないので、開始点を変えて拾う。
    /// 260726Cl 変更: 旧は軸方向 6 点+対角 3 点の手書き 10 点。点数を増やすには系統的な列が要る</summary>
    static readonly (double U, double V, double D)[] StartOffsets = BuildStartOffsets(StartCount);

    static (double U, double V, double D)[] BuildStartOffsets(int count)
    {
        //Halton 列 (基数 2,3,5) を [0,1) → [-1,1] へ。低食い違い列なので、点数を増やすほど隙間なく埋まる
        static double Halton(int index, int b)
        {
            double f = 1, r = 0;
            for (int i = index; i > 0; i /= b) { f /= b; r += f * (i % b); }
            return r;
        }
        var offsets = new (double U, double V, double D)[count];
        offsets[0] = (0, 0, 0); //現在の幾何そのもの
        for (int i = 1; i < count; i++)
            offsets[i] = (2 * Halton(i, 2) - 1, 2 * Halton(i, 3) - 1, 2 * Halton(i, 5) - 1);
        return offsets;
    }

    /// <summary>260920Cl 追加: 較正の比較解像度 (長辺 px)。粗 → 中 → フル解像度の 3 段で、段ごとに前段の解を引き継ぐ。
    /// 旧実装は全段 160 px 固定で、1 画素が検出器上 0.4 mm 相当 (1344 px の検出器) と最終の詰めには粗すぎた。
    /// 最終段は 0 = 実測画像のフル解像度。作者指示「計算時間は長くかかっても構わない。なるべく元の解像度で」</summary>
    static readonly int[] StageLongSides = [160, 480, 0];

    /// <summary>260920Cl 追加: 2 段目・3 段目へ持ち上げる上位解の数。1 段目 (粗) の谷が浅いと最良が入れ替わるため 1 点に絞らない</summary>
    //260920Cl: 中段の点数はコア数に合わせる (6 点だと多コア機で遊ぶ)。最終段は 1 点にして評価の内部を並列にする
    static readonly int[] StageKeep = [Math.Max(6, Environment.ProcessorCount), 1];

    /// <summary>260920Cl 追加: 各段の Nelder-Mead の収束条件。段が細かくなるほど厳しくする
    /// (関数値の幅 tol と、初期刻みに対するシンプレックスの広がり xtolRel の両方)</summary>
    static readonly (double Tol, double XTol)[] StageTolerance = [(1E-6, 1E-2), (1E-7, 3E-3), (1E-9, 3E-4)];

    /// <summary>260920Cl 追加: 各段の評価上限 (方位段 / 幾何段 / 6 変数同時段)。フル解像度は 1 評価が重いので点数を絞る代わりに上限を上げる</summary>
    static readonly (int Ori, int Geo, int Joint)[] StageMaxEval = [(150, 120, 600), (200, 180, 900), (300, 250, 1500)];

    /// <summary>260920Cl 追加: 最終段で 6 変数同時最適化を刻みを縮めて張り直す回数。素の Nelder-Mead は谷で停滞する</summary>
    const int JointRestarts = 3;

    /// <summary>260920Cl 追加: 1 つの比較解像度。Reference は正規化済み (ZNCC の参照)、FlattenFwhm はシミュレーション側に掛ける高域通過の半値幅 [この解像度の px]</summary>
    sealed class Scale { public int W, H; public double[] Reference; public double FlattenFwhm; }

    /// <summary>260920Cl 追加 (作者指示): 幾何を固定して方位だけ動かすときの採点器。較正の**最終段とまったく同じ**
    /// 参照 (表示中の実測値をフル解像度で正規化)・同じ高域通過・同じ ZNCC を使う。
    /// 方位探索の仕上げ段がこれを通ることで Find と Calibrate の目的関数が一致する。
    /// 両者が違うと「Find → トップ選択 → Calibrate → 再び Find」で方位が 2 値を約 1° で往復して収束しない (作者の実機報告)。
    /// DisplayReference が無い旧経路では context.Reference (縮小 + 強制背景除算) の 1 段になるので、従来の採点と同じになる</summary>
    public sealed class FixedGeometryScorer
    {
        readonly EbsdMatchingContext context;
        readonly Scale scale;
        readonly EbsdPatternProjector projector;
        readonly double[] buf, work1, work2;

        /// <summary>採点に使う解像度 (較正の最終段と同じ)</summary>
        public int Width => scale.W;
        /// <summary>同上</summary>
        public int Height => scale.H;

        internal FixedGeometryScorer(EbsdMatchingContext context)
        {
            this.context = context;
            scale = BuildScales(context)[^1]; //最終段 = 最も細かい解像度
            projector = new EbsdPatternProjector(context.Geometry, scale.W, scale.H);
            int n = scale.W * scale.H;
            buf = new double[n]; work1 = new double[n]; work2 = new double[n];
        }

        /// <summary>与えた方位での ZNCC。値が大きいほど実測に近い (較正は −ZNCC を最小化している)</summary>
        public double Zncc(Matrix3D rotation)
        {
            projector.Project(context.MasterPattern, rotation, context.PositivePlane, context.NegativePlane, buf);
            if (scale.FlattenFwhm >= 1) EbsdPatternScorer.SubtractBoxBackground(buf, work1, work2, scale.W, scale.H, scale.FlattenFwhm, true);
            return EbsdPatternScorer.ZnccSimd(scale.Reference, buf, true);
        }
    }

    /// <summary>260920Cl 追加 (作者指示): 較正の最終段と同じ採点器を作る。方位探索の仕上げ段と目的関数を揃えるための入口</summary>
    public static FixedGeometryScorer CreateFinalScorer(EbsdMatchingContext context)
    {
        ArgumentNullException.ThrowIfNull(context);
        return new FixedGeometryScorer(context);
    }

    /// <summary>260920Cl 追加: スレッドごとの作業バッファ (投影先・box blur の作業用・使い回すプロジェクタ)。
    /// Proj は幾何が変わるたび Rebuild するだけで、配列の確保はこの 1 回きり</summary>
    sealed class Work
    {
        public double[] Buf, W1, W2;
        public EbsdPatternProjector Proj;
        public Work(Scale sc, EbsdDetectorGeometry geom)
        {
            int n = sc.W * sc.H;
            Buf = new double[n]; W1 = new double[n]; W2 = new double[n];
            Proj = new EbsdPatternProjector(geom, sc.W, sc.H);
        }
    }

    /// <summary>PC/DD と方位を較正する。結果は <see cref="EbsdDetectorGeometry.FromPatternCenter"/> で DetX/DetY/DetZ へ戻す。
    /// 260920Cl 全面改修 (作者指示): 粗 → 中 → フル解像度の多段、多点開始の並列化、シミュレーション側にも実測と同じ平坦化、収束判定の強化。
    /// 旧実装 (全段 160 px・逐次・シミュレーション側は生値・関数値だけの収束判定) は git 履歴 260727Cl 版を参照</summary>
    /// <param name="context">実測パターン・MasterPattern・現在の幾何と方位のスナップショット</param>
    /// <param name="detectorWidthMm">検出器の物理幅 (mm)。ソフト境界と Nelder-Mead の初期ステップに使う</param>
    /// <param name="detectorHeightMm">検出器の物理高さ (mm)</param>
    /// <param name="cancel">中止トークン</param>
    /// <param name="progress">進捗 (0-1) と段の名前。ワーカースレッドから呼ばれるので受け手側でマーシャリングすること</param>
    public static EbsdCalibrationResult Run(EbsdMatchingContext context, double detectorWidthMm, double detectorHeightMm,
        CancellationToken cancel = default, Action<double, string> progress = null)
    {
        ArgumentNullException.ThrowIfNull(context);
        if (!(detectorWidthMm > 0) || !double.IsFinite(detectorWidthMm)) throw new ArgumentOutOfRangeException(nameof(detectorWidthMm));
        if (!(detectorHeightMm > 0) || !double.IsFinite(detectorHeightMm)) throw new ArgumentOutOfRangeException(nameof(detectorHeightMm));

        var geom0 = context.Geometry;
        double detTilt = geom0.DetTilt, smpTilt = geom0.SampleTilt, xm = geom0.XMirror, pixelSize = geom0.PixelSize;
        int imgW = geom0.WidthPx, imgH = geom0.HeightPx;
        var (footU0, footV0) = geom0.PatternCenterMm;
        double dd0 = geom0.CameraLength, lnDd0 = Math.Log(dd0);
        double physW = detectorWidthMm, physH = detectorHeightMm;

        //--- 比較スケールを作る。DisplayReference (フル解像度の表示値) があればそれを段ごとに縮小 + 正規化する。
        //    無ければ旧来どおり context.Reference (縮小 + 強制背景除算済み) の 1 段だけで動く
        var scales = BuildScales(context);

        EbsdDetectorGeometry MakeGeom(double u, double v, double ld)
        {
            var (dx, dy, dz) = EbsdDetectorGeometry.FromPatternCenter(u, v, Math.Exp(ld), detTilt);
            return new EbsdDetectorGeometry(detTilt, dx, dy, dz, pixelSize, imgW, imgH, xm, smpTilt);
        }
        bool OutOfSoftBounds(double u, double v, double lnDd, out double penalty)
        {
            double du = u - footU0, dv = v - footV0, dlnDd = lnDd - lnDd0;
            penalty = SoftBoundPenaltyBase + Math.Abs(du) / physW + Math.Abs(dv) / physH + Math.Abs(dlnDd);
            return Math.Abs(du) > physW * SoftBoundPcFraction || Math.Abs(dv) > physH * SoftBoundPcFraction
                || Math.Abs(dlnDd) > SoftBoundLnDd;
        }

        int evalTotal = 0; //260920Cl (/simplify): 進捗は完了した開始点の数で出すので、評価回数の逐次カウンタ (旧 evalsDone) は不要
        double ScoreWith(Scale sc, Work w, EbsdPatternProjector proj, Matrix3D rot, bool innerParallel)
        {
            cancel.ThrowIfCancellationRequested();
            proj.Project(context.MasterPattern, rot, context.PositivePlane, context.NegativePlane, w.Buf, innerParallel);
            //260920Cl: 実測側が平坦化されているならシミュレーション側にも同じ高域通過を掛ける (これを欠くと ZNCC がモデル由来の背景勾配に引かれる)
            if (sc.FlattenFwhm >= 1) EbsdPatternScorer.SubtractBoxBackground(w.Buf, w.W1, w.W2, sc.W, sc.H, sc.FlattenFwhm, innerParallel); //260920Cl: 最終段は内部も並列
            return -EbsdPatternScorer.ZnccSimd(sc.Reference, w.Buf, innerParallel); //260920Cl: SIMD 版。最終段は並列
        }
        //260920Cl: 旧実装は評価ごとに new EbsdPatternProjector していた (フル解像度で 1 評価あたり 33 MB の確保)。バッファを使い回す
        //旧: => ScoreWith(sc, w, new EbsdPatternProjector(MakeGeom(fu, fv, lnDd), sc.W, sc.H), rot, innerParallel);
        double ScoreAt(Scale sc, Work w, double fu, double fv, double lnDd, Matrix3D rot, bool innerParallel)
        {
            w.Proj.Rebuild(MakeGeom(fu, fv, lnDd), innerParallel);
            return ScoreWith(sc, w, w.Proj, rot, innerParallel);
        }

        var coarse = scales[0];
        //260920Cl: 較正前後の ZNCC は**同じ解像度**で測る。段ごとに ZNCC の絶対値は変わる (細かいほど下がる) ので、
        //  旧: 開始値を粗い段、結果を最終段で測っており、改善しても悪化したように見えていた
        var finest = scales[^1];
        double startZncc = -ScoreAt(finest, new Work(finest, geom0), footU0, footV0, lnDd0, context.Rotation, true);

        //--- 1 開始点ぶんの較正 (交互法 → 方位仕上げ → 6 変数同時 (最終段は再起動付き))
        (double Zncc, double Fu, double Fv, double LnDd, Matrix3D Rot, int Rounds, bool Converged, double JointGain)
            RunFrom(int stage, double fu, double fv, double lnDd, Matrix3D rot, bool innerParallel)
        {
            var sc = scales[Math.Min(stage, scales.Length - 1)];
            var (tol, xtol) = StageTolerance[Math.Min(stage, StageTolerance.Length - 1)];
            var (maxOri, maxGeo, maxJoint) = StageMaxEval[Math.Min(stage, StageMaxEval.Length - 1)];
            var w = new Work(sc, geom0);
            var r0 = rot;
            int roundsUsed = 0; bool converged = false;
            double prevZncc = -ScoreAt(sc, w, fu, fv, lnDd, r0, innerParallel);
            for (int round = 0; round < MaxRounds; round++)
            {
                cancel.ThrowIfCancellationRequested();
                //① 幾何固定で方位 (粗 0.7°)
                //260920Cl: 方位段は幾何が固定なので視線を 1 回だけ作り直して使い回す (旧: 毎ラウンド new = フル解像度で 33 MB)
                w.Proj.Rebuild(MakeGeom(fu, fv, lnDd), innerParallel);
                var projFixed = w.Proj;
                var (bo, _, eo) = EbsdPatternScorer.NelderMead(v => ScoreWith(sc, w, projFixed, EbsdIndexer.PerturbRotation(r0, v[0], v[1], v[2]), innerParallel),
                    [0, 0, 0], [0.7, 0.7, 0.7], maxOri, tol, xtol);
                r0 = EbsdIndexer.PerturbRotation(r0, bo[0], bo[1], bo[2]); Interlocked.Add(ref evalTotal, eo);

                //② 方位固定で幾何 (dU, dV [mm], dlnDD)。ステップ = 検出器幅/高の 1%、lnDD 0.02
                var rFixed = r0;
                var (bg, vg, eg) = EbsdPatternScorer.NelderMead(
                    v => OutOfSoftBounds(fu + v[0], fv + v[1], lnDd + v[2], out var pen) ? pen
                        : ScoreAt(sc, w, fu + v[0], fv + v[1], lnDd + v[2], rFixed, innerParallel),
                    [0, 0, 0], [physW * 0.01, physH * 0.01, 0.02], maxGeo, tol, xtol);
                fu += bg[0]; fv += bg[1]; lnDd += bg[2]; Interlocked.Add(ref evalTotal, eg);
                roundsUsed = round + 1;

                double zncc = -vg;
                if (zncc - prevZncc < ZnccTolerance) { converged = true; break; }
                prevZncc = zncc;
            }
            //仕上げの方位微調整
            const double polishStep = EbsdOrientationSearch.OrientationPolishStepDeg;
            w.Proj.Rebuild(MakeGeom(fu, fv, lnDd), innerParallel); //260920Cl: 同上
            var projFinal = w.Proj;
            var (bf, vf, ef) = EbsdPatternScorer.NelderMead(v => ScoreWith(sc, w, projFinal, EbsdIndexer.PerturbRotation(r0, v[0], v[1], v[2]), innerParallel),
                [0, 0, 0], [polishStep, polishStep, polishStep], maxOri, tol, xtol);
            r0 = EbsdIndexer.PerturbRotation(r0, bf[0], bf[1], bf[2]); Interlocked.Add(ref evalTotal, ef);

            //6 変数 (PC_u, PC_v, lnDD, 方位 3) 同時最適化。交互法は相関のある谷を斜めに下れないのでここで下る
            var rBase = r0;
            double fuBase = fu, fvBase = fv, lnDdBase = lnDd;
            double ScoreJoint(double[] v)
            {
                if (OutOfSoftBounds(fuBase + v[0], fvBase + v[1], lnDdBase + v[2], out var pen)) return pen;
                return ScoreAt(sc, w, fuBase + v[0], fvBase + v[1], lnDdBase + v[2], EbsdIndexer.PerturbRotation(rBase, v[3], v[4], v[5]), innerParallel);
            }
            double[] jointSteps = [physW * 0.005, physH * 0.005, 0.01, polishStep, polishStep, polishStep];
            //260920Cl: 最終段だけ刻みを縮めて張り直す (再起動)。粗い段でそこまで詰めても次段で作り直すので無駄
            var (bj, vj, ej) = stage >= scales.Length - 1
                ? EbsdPatternScorer.NelderMeadRestart(ScoreJoint, [0, 0, 0, 0, 0, 0], jointSteps, maxJoint, tol, xtol, JointRestarts)
                : EbsdPatternScorer.NelderMead(ScoreJoint, [0, 0, 0, 0, 0, 0], jointSteps, maxJoint, tol, xtol);
            fu = fuBase + bj[0]; fv = fvBase + bj[1]; lnDd = lnDdBase + bj[2];
            r0 = EbsdIndexer.PerturbRotation(rBase, bj[3], bj[4], bj[5]); Interlocked.Add(ref evalTotal, ej);

            return (Zncc: -vj, Fu: fu, Fv: fv, LnDd: lnDd, Rot: r0, Rounds: roundsUsed, Converged: converged, JointGain: -vj - -vf);
        }

        //--- 段 0: 粗い解像度で多点開始 (並列)。1 点が独立なので Parallel.For し、内側の投影は逐次にする
        var runs = new (double Zncc, double Fu, double Fv, double LnDd, Matrix3D Rot, int Rounds, bool Converged, double JointGain)[StartOffsets.Length];
        int completed = 0;
        progress?.Invoke(0, $"stage 1/{scales.Length}: {StartOffsets.Length} starts at {coarse.W}x{coarse.H}");
        Parallel.For(0, StartOffsets.Length, new ParallelOptions { CancellationToken = cancel }, i =>
        {
            var (ou, ov, od) = StartOffsets[i];
            runs[i] = RunFrom(0, footU0 + ou * physW * StartSpreadPc, footV0 + ov * physH * StartSpreadPc, lnDd0 + od * StartSpreadLnDd, context.Rotation, false);
            int done = Interlocked.Increment(ref completed);
            progress?.Invoke(0.5 * done / StartOffsets.Length, $"stage 1/{scales.Length}: {done}/{StartOffsets.Length} starts");
        });

        //--- 診断値 (幾何がどこまで決まっているか) は多点開始を行ったこの段で測る
        var order = Enumerable.Range(0, runs.Length).OrderByDescending(i => runs[i].Zncc).ToArray();
        int bestIndex = order[0];
        double worstZncc = runs.Min(r => r.Zncc), coarseBest = runs[bestIndex].Zncc;
        var near = runs.Where(r => r.Zncc >= coarseBest - 1E-3).ToArray();
        double flatU = (near.Max(r => r.Fu) - near.Min(r => r.Fu)) / 2;
        double flatV = (near.Max(r => r.Fv) - near.Min(r => r.Fv)) / 2;
        double flatDd = (near.Max(r => Math.Exp(r.LnDd)) - near.Min(r => Math.Exp(r.LnDd))) / 2;

        //--- 段 1 以降: 上位だけを細かい解像度で引き継ぐ。最終段は 1 点なので内側 (投影) を並列にする
        var carried = order.ToArray();
        for (int stage = 1; stage < scales.Length; stage++)
        {
            int keep = Math.Min(StageKeep[Math.Min(stage - 1, StageKeep.Length - 1)], carried.Length);
            var sc = scales[stage];
            var stageRuns = new (double Zncc, double Fu, double Fv, double LnDd, Matrix3D Rot, int Rounds, bool Converged, double JointGain)[keep];
            int doneStage = 0;
            progress?.Invoke(0.5 + 0.5 * (stage - 1) / Math.Max(1, scales.Length - 1), $"stage {stage + 1}/{scales.Length}: {keep} start(s) at {sc.W}x{sc.H}");
            if (keep > 1)
                Parallel.For(0, keep, new ParallelOptions { CancellationToken = cancel }, k =>
                {
                    var src = runs[carried[k]];
                    stageRuns[k] = RunFrom(stage, src.Fu, src.Fv, src.LnDd, src.Rot, false);
                    int d = Interlocked.Increment(ref doneStage);
                    progress?.Invoke(0.5 + 0.5 * ((stage - 1) + (double)d / keep) / Math.Max(1, scales.Length - 1), null);
                });
            else
            {
                var src = runs[carried[0]];
                stageRuns[0] = RunFrom(stage, src.Fu, src.Fv, src.LnDd, src.Rot, true);
            }
            //次段へ渡すため runs/carried を差し替える
            runs = stageRuns;
            carried = Enumerable.Range(0, keep).OrderByDescending(i => stageRuns[i].Zncc).ToArray();
        }
        var bestRun = runs[carried[0]];

        return new EbsdCalibrationResult(bestRun.Rot, bestRun.Fu, bestRun.Fv, Math.Exp(bestRun.LnDd), bestRun.Zncc, startZncc,
            evalTotal, bestRun.Rounds, bestRun.Converged, bestRun.JointGain,
            StartOffsets.Length, bestIndex, coarseBest - worstZncc, near.Length,
            flatU, flatV, flatDd);
    }

    /// <summary>260920Cl 追加: 比較スケールを作る。DisplayReference (表示中の実測値、フル解像度) を各段の長辺へ box 縮小し、
    /// ZNCC 用に正規化する。シミュレーション側へ掛ける高域通過の半値幅も同じ縮小率で換算する。
    /// DisplayReference が無い場合 (旧経路) は context.Reference の 1 段だけを返す</summary>
    static Scale[] BuildScales(EbsdMatchingContext context)
    {
        if (context.DisplayReference == null || context.DisplayWidth <= 0 || context.DisplayHeight <= 0
            || context.DisplayReference.Length != context.DisplayWidth * context.DisplayHeight)
            return [new Scale { W = context.RasterWidth, H = context.RasterHeight, Reference = context.Reference, FlattenFwhm = 0 }];

        int fullW = context.DisplayWidth, fullH = context.DisplayHeight;
        var list = new List<Scale>();
        foreach (var target in StageLongSides)
        {
            int longSide = target > 0 ? Math.Min(target, Math.Max(fullW, fullH)) : Math.Max(fullW, fullH);
            var (data, w, h) = EbsdPatternScorer.Downsample(context.DisplayReference, fullW, fullH, longSide);
            if (list.Count > 0 && list[^1].W == w && list[^1].H == h) continue; //同じ解像度が続いたら 1 段にまとめる
            EbsdPatternScorer.NormalizeInPlace(data);
            list.Add(new Scale { W = w, H = h, Reference = data, FlattenFwhm = context.SimFlattenFwhmPx * w / fullW });
        }
        return [.. list];
    }
}
