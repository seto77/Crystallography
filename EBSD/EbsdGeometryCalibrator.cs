#region using
using System;
using System.Linq;
using System.Threading;
#endregion

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
    /// 較正のソフト境界 (初期値から W/H の 25%、lnDD 0.35) の内側に収めること</summary>
    const double StartSpreadPc = 0.08, StartSpreadLnDd = 0.08;

    /// <summary>較正の最後に行う 6 変数 (PC_u, PC_v, lnDD, 方位 3) 同時最適化の評価上限。260726Cl 追加。
    /// 6 次元なので交互法の 3 変数段 (120-150) より多く要る。1 評価ごとに projector を作り直す重い段だが、
    /// 交互法では下れない斜めの谷をここで下る</summary>
    const int JointPolishMaxEval = 600;

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

    /// <summary>PC/DD と方位を較正する。結果は <see cref="EbsdDetectorGeometry.FromPatternCenter"/> で DetX/DetY/DetZ へ戻す。</summary>
    /// <param name="context">実測パターン・MasterPattern・現在の幾何と方位のスナップショット</param>
    /// <param name="detectorWidthMm">検出器の物理幅 (mm)。ソフト境界と Nelder-Mead の初期ステップに使う</param>
    /// <param name="detectorHeightMm">検出器の物理高さ (mm)</param>
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
        var (footU0, footV0) = geom0.PatternCenterMm; //260724Cl (/simplify): PC 式の手書き重複 (-DetX, -(DetY cosδ+DetZ sinδ)) を幾何オブジェクトへ一元化
        double dd0 = geom0.CameraLength;
        double physW = detectorWidthMm, physH = detectorHeightMm;

        var buf = new double[context.RasterWidth * context.RasterHeight];
        int evalTotal = 0;
        //260726Cl 変更 (作者報告「プログレスバーの挙動がおかしい」): 旧実装は「評価回数 / 静的な予算」で進捗を出していたが、
        //予算は最大ラウンド (20) を使い切る前提なのに実際は 1-2 ラウンドで収束するため、バーは 3 割ほどで止まって最後に 100% へ飛んでいた
        //(実測 334,551 評価 / 予算 1,220,000)。**完了した開始点の数**を主軸にし、実行中の開始点の内側だけを
        //「これまでの 1 点あたり実測平均」で按分する。1 点目だけは実測が無いので静的な予算で見積もる。
        //旧: int evalBudget = StartOffsets.Length * PerStartBudget; ratio = evalsDone / evalBudget
        int evalsDone = 0;
        const int PerStartBudget = MaxRounds * (150 + 120) + 100 + JointPolishMaxEval;
        int completedStarts = 0, evalsAtStartBegin = 0;
        double avgEvalsPerStart = PerStartBudget;

        EbsdDetectorGeometry MakeGeom(double u, double v, double ld)
        {
            var (dx, dy, dz) = EbsdDetectorGeometry.FromPatternCenter(u, v, Math.Exp(ld), detTilt);
            return new EbsdDetectorGeometry(detTilt, dx, dy, dz, pixelSize, imgW, imgH, xm, smpTilt);
        }
        //260727Cl (/simplify): soft bounds の判定とペナルティ式が交互法② と 6 変数同時仕上げの 2 箇所に同じ形で書かれ、
        //  閾値 (W/H の 25%・lnDD 0.35) とペナルティ基底 10 も 2 重にハードコードされていたので 1 本にまとめた。
        //  ペナルティ値は常に計算するが副作用が無いので、境界内で使われないだけ (式・戻り値は旧実装と同一)。
        bool OutOfSoftBounds(double du, double dv, double dlnDd, out double penalty)
        {
            penalty = 10 + Math.Abs(du) / physW + Math.Abs(dv) / physH + Math.Abs(dlnDd);
            return Math.Abs(du) > physW * 0.25 || Math.Abs(dv) > physH * 0.25 || Math.Abs(dlnDd) > 0.35;
        }
        double ScoreWith(EbsdPatternProjector proj, Matrix3D rot)
        {
            cancel.ThrowIfCancellationRequested(); //260725Ch: 各評価の投影前に中止を反映
            //260726Cl: 完了した開始点 + 実行中の開始点の按分。NM は逐次なので単純加算で足りる
            evalsDone++;
            double inCurrentStart = Math.Min(0.99, (evalsDone - evalsAtStartBegin) / Math.Max(1, avgEvalsPerStart));
            progress?.Invoke(Math.Min(0.99, (completedStarts + inCurrentStart) / StartOffsets.Length), null);
            proj.Project(context.MasterPattern, rot, context.PositivePlane, context.NegativePlane, buf);
            return -EbsdPatternScorer.Zncc(context.Reference, buf);
        }
        double startZncc = -ScoreWith(new EbsdPatternProjector(MakeGeom(footU0, footV0, Math.Log(dd0)), context.RasterWidth, context.RasterHeight), context.Rotation);

        //260726Cl 追加 (作者要望): 1 開始点ぶんの較正 (交互法 → 方位仕上げ → 6 変数同時) を関数化し、多点開始から呼ぶ
        (double Zncc, double Fu, double Fv, double LnDd, Matrix3D Rot, int Rounds, bool Converged, double JointGain) RunFrom(double fu, double fv, double lnDd)
        {
            var r0 = context.Rotation;
            //260725Cl 変更 (作者指示): 旧 for (int round = 0; round < 2; round++) — 2 ラウンド固定で収束判定なし。
            //PC・DD・方位の相関で交互法はジグザグするため、改善が止まるまで最大 MaxRounds 回まわす
            int roundsUsed = 0;
            bool converged = false;
            double prevZncc = -ScoreWith(new EbsdPatternProjector(MakeGeom(fu, fv, lnDd), context.RasterWidth, context.RasterHeight), r0);
            for (int round = 0; round < MaxRounds; round++)
            {
                cancel.ThrowIfCancellationRequested(); //260725Ch
                //① 幾何固定で方位 (粗 0.7°)
                var projFixed = new EbsdPatternProjector(MakeGeom(fu, fv, lnDd), context.RasterWidth, context.RasterHeight);
                var (bo, _, eo) = EbsdPatternScorer.NelderMead(v => ScoreWith(projFixed, EbsdIndexer.PerturbRotation(r0, v[0], v[1], v[2])), [0, 0, 0], [0.7, 0.7, 0.7], 150);
                r0 = EbsdIndexer.PerturbRotation(r0, bo[0], bo[1], bo[2]); evalTotal += eo;

                //② 方位固定で幾何 (dU, dV [mm], dlnDD)。ステップ = 検出器幅/高の 1%、lnDD 0.02
                //260724Cl: 単一パターンの PC-DD-方位縮退で非物理領域へ流れないよう soft bounds (初期値から W/H の 25%・DD ±40% でペナルティ)
                var rFixed = r0;
                var (bg, vg, eg) = EbsdPatternScorer.NelderMead(
                    v => OutOfSoftBounds(v[0], v[1], v[2], out var pen) ? pen //260727Cl: 判定+罰則式を OutOfSoftBounds へ集約
                        : ScoreWith(new EbsdPatternProjector(MakeGeom(fu + v[0], fv + v[1], lnDd + v[2]), context.RasterWidth, context.RasterHeight), rFixed),
                    [0, 0, 0], [physW * 0.01, physH * 0.01, 0.02], 120);
                fu += bg[0]; fv += bg[1]; lnDd += bg[2]; evalTotal += eg;
                roundsUsed = round + 1;

                //260725Cl: このラウンドの ZNCC 到達点で収束判定 (soft bounds のペナルティ値が返った場合は改善なしとして扱われる)
                double zncc = -vg;
                if (zncc - prevZncc < ZnccTolerance) { converged = true; break; }
                prevZncc = zncc;
            }
            //仕上げの方位微調整。260725Cl 変更: 0.2° → OrientationPolishStepDeg (0.1°、作者指示)。Find の仕上げ段と同じ値
            const double polishStep = EbsdOrientationSearch.OrientationPolishStepDeg;
            var projFinal = new EbsdPatternProjector(MakeGeom(fu, fv, lnDd), context.RasterWidth, context.RasterHeight);
            var (bf, vf, ef) = EbsdPatternScorer.NelderMead(v => ScoreWith(projFinal, EbsdIndexer.PerturbRotation(r0, v[0], v[1], v[2])),
                [0, 0, 0], [polishStep, polishStep, polishStep], 100);
            r0 = EbsdIndexer.PerturbRotation(r0, bf[0], bf[1], bf[2]); evalTotal += ef;

            //260726Cl 追加 (作者要望): 6 変数 (PC_u, PC_v, lnDD, 方位 3) の同時最適化を仕上げに 1 段。
            //交互法は変数を片方ずつしか動かせないので、相関のある谷では斜め方向に下れずジグザグして止まる。
            //実機報告でも初期 DetX/Y/Z を変えると最終スコアが 20.0〜20.3 程度ばらついていた。
            //開始点 (増分ゼロ) が初期シンプレックスの頂点 0 で、NelderMead は最良頂点を返すので、この段で悪化することはない。
            //ソフト境界は交互法の②と同じ判定を増分に対して掛ける (この段の増分は小さいので通常は発火しない)。
            var rBase = r0;
            double fuBase = fu, fvBase = fv, lnDdBase = lnDd;
            double ScoreJoint(double[] v)
            {
                if (OutOfSoftBounds(v[0], v[1], v[2], out var pen)) return pen; //260727Cl: 交互法②と同じ判定を共通関数へ
                return ScoreWith(new EbsdPatternProjector(MakeGeom(fuBase + v[0], fvBase + v[1], lnDdBase + v[2]), context.RasterWidth, context.RasterHeight),
                    EbsdIndexer.PerturbRotation(rBase, v[3], v[4], v[5]));
            }
            //幾何側は交互法②の半分のステップ (もう最適点の近くにいる)、方位側は仕上げと同じ 0.1°
            var (bj, vj, ej) = EbsdPatternScorer.NelderMead(ScoreJoint, [0, 0, 0, 0, 0, 0],
                [physW * 0.005, physH * 0.005, 0.01, polishStep, polishStep, polishStep], JointPolishMaxEval);
            fu = fuBase + bj[0]; fv = fvBase + bj[1]; lnDd = lnDdBase + bj[2];
            r0 = EbsdIndexer.PerturbRotation(rBase, bj[3], bj[4], bj[5]); evalTotal += ej;

            return (Zncc: -vj, Fu: fu, Fv: fv, LnDd: lnDd, Rot: r0, Rounds: roundsUsed, Converged: converged,
                JointGain: -vj - -vf); //260726Cl: 同時最適化が交互法の到達点からどれだけ伸ばしたか
        }

        //260726Cl 追加 (作者要望): 多点開始。局所解が多く、初期 DetX/Y/Z を変えると最終スコアが 0.3 程度ばらつくため、
        //現在の幾何と、そこから決定的に振った開始点から同じ較正を走らせ、最も ZNCC の高い解を採る。
        //同時最適化は交互法の停滞は解消するが局所解の壁は越えないので、壁の向こう側は開始点を変えて拾うしかない
        (double Zncc, double Fu, double Fv, double LnDd, Matrix3D Rot, int Rounds, bool Converged, double JointGain) bestRun = default;
        int bestIndex = -1;
        double worstZncc = double.MaxValue;
        var runs = new (double Zncc, double Fu, double Fv, double Dd)[StartOffsets.Length]; //260726Cl: 最良解へ到達した点の数と、その幾何の広がりを見るため
        for (int s = 0; s < StartOffsets.Length; s++)
        {
            cancel.ThrowIfCancellationRequested();
            progress?.Invoke(Math.Min(0.99, (double)s / StartOffsets.Length), $"start {s + 1}/{StartOffsets.Length}"); //260726Cl
            var (ou, ov, od) = StartOffsets[s];
            //260726Cl 変更: 振れ幅を定数化 (旧 physW*0.01 / physH*0.01 / 0.02 は狭すぎて全点が同じ谷に落ちていた)
            var run = RunFrom(footU0 + ou * physW * StartSpreadPc, footV0 + ov * physH * StartSpreadPc,
                Math.Log(dd0) + od * StartSpreadLnDd);
            runs[s] = (run.Zncc, run.Fu, run.Fv, Math.Exp(run.LnDd));
            worstZncc = Math.Min(worstZncc, run.Zncc);
            if (bestIndex < 0 || run.Zncc > bestRun.Zncc) { bestRun = run; bestIndex = s; }
            //260726Cl: 進捗の按分に使う「1 点あたりの実測評価数」を更新する
            completedStarts = s + 1;
            avgEvalsPerStart = (double)evalsDone / completedStarts;
            evalsAtStartBegin = evalsDone;
        }
        //最良から 1E-3 以内に入った開始点の数 = 最良解の basin の広さ。spread (最良−最悪) だけだと外れ値に引きずられる
        var near = runs.Where(r => r.Zncc >= bestRun.Zncc - 1E-3).ToArray();
        //260726Cl 追加 (作者要望): その集団の PC・DD の広がり (半値幅) = ZNCC で幾何がどこまで決まっているか。
        //ZNCC 1E-3 以内で PC が数 mm 動くなら、単一パターンでは幾何がその精度までしか決まっていない (正本 §2.4)
        double flatU = (near.Max(r => r.Fu) - near.Min(r => r.Fu)) / 2;
        double flatV = (near.Max(r => r.Fv) - near.Min(r => r.Fv)) / 2;
        double flatDd = (near.Max(r => r.Dd) - near.Min(r => r.Dd)) / 2;

        return new EbsdCalibrationResult(bestRun.Rot, bestRun.Fu, bestRun.Fv, Math.Exp(bestRun.LnDd), bestRun.Zncc, startZncc,
            evalTotal, bestRun.Rounds, bestRun.Converged, bestRun.JointGain,
            StartOffsets.Length, bestIndex, bestRun.Zncc - worstZncc, near.Length, //260726Cl: 局所解のばらつきを可視化
            flatU, flatV, flatDd); //260726Cl: ZNCC が同等な解の集団における PC・DD の広がり (半値幅、mm)
    }
}
