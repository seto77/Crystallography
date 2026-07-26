#region using
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
#endregion

namespace Crystallography;

/// <summary>
/// ZNCC 系 (辞書照合・方位仕上げ・幾何較正) が必要とする状態のスナップショット。260726Cl 追加。
/// GUI 側 (FormEBSD.SnapshotMatchingContext) が UI スレッド上で 1 回だけ作り、以降はワーカースレッドから読み取り専用で使う。
/// Positive/NegativePlane は「MC 重み合成パターン」または単一 energy/depth スライス、Reference は前処理済み実測パターン
/// (RasterWidth×RasterHeight に縮小済み)、Rotation は現在の結晶方位。
/// </summary>
public sealed record EbsdMatchingContext(
    EbsdDetectorGeometry Geometry,
    MasterPattern MasterPattern,
    float[] PositivePlane,
    float[] NegativePlane,
    double[] Reference,
    int RasterWidth,
    int RasterHeight,
    Matrix3D Rotation);

/// <summary>
/// 実測 EBSD パターンからの方位候補探索 (Radon テンプレート照合 or MasterPattern 辞書照合 + ZNCC 複合ランク + 仕上げ)。
/// 260726Cl 追加: FormEBSD.cs の buttonFindOrientation_Click 内 Task.Run 本体をそのまま移設したもの (GUI 非依存)。 //260727Cl: 旧 FormEBSD.Indexing.cs は同日 FormEBSD.cs へ統合済みなのでファイル名を訂正
/// UI 側に残るのは前提チェック・スナップショット作成・進捗表示・結果の適用だけ。
///
/// 260724Cl 改訂 (ベンチ+Codex 裁定、設計正本 §2.1): 生 ZNCC の再ランクは有害 (シミュレーションの heavy-tailed 生強度が支配し
/// 正解方位が偽方位に負ける) と実測で判明。現在の構成は
///   ① Radon 採点は複合前提のとき証拠飽和 cap=<see cref="SaturateCap"/> (少数強リッジ支配の抑制。単独では 5-2_22 のトップが劣化するため複合とセットでのみ使う)
///   ② 実測・シミュレーション両方に RobustPreprocess を掛けた ZNCC を候補集合内で標準化し、combo = zRadon + 0.5·clip(z,±2) で再ランク
///   ③ ZNCC 精密化は複合トップ 1 件のみ ±0.25° (ガード: Radon z 低下が <see cref="GuardMaxZDrop"/> 超で棄却)。ベンチ 3 画像で複合トップ全勝 (12/20, 5/15, 11/14) //260727Cl: 裸の 0.2 を定数参照へ
/// </summary>
public static class EbsdOrientationSearch
{
    /// <summary>Radon 採点の証拠飽和 cap。260724Cl: EbsdIndexCheck ハーネスの係数スイープで決定 (プラトー 0.4-1.0 の中央寄り)</summary>
    public const double SaturateCap = 8;

    /// <summary>複合ランクにおける標準化 ZNCC の係数</summary>
    const double ZnccCoef = 0.5;

    /// <summary>方位仕上げ (Find のトップ候補・較正の最終段) の Nelder-Mead 初期ステップ [°]。260725Cl 追加: 0.2 → 0.1 (作者指示)。
    /// <see cref="EbsdGeometryCalibrator"/> と同じ値を使う — 目的関数だけでなくステップも揃えないと、Find と Calibrate を
    /// 繰り返したときに方位が微妙に往復する</summary>
    public const double OrientationPolishStepDeg = 0.1;

    /// <summary>誤収束ガードで許容する Radon z の劣化量。260727Cl (/simplify): 精密化と仕上げの 2 箇所に
    /// 裸のリテラル 0.2 で書かれていたので命名した (SaturateCap/ZnccCoef と同じ扱いに揃える)</summary>
    const double GuardMaxZDrop = 0.2;

    /// <summary>方位候補を探索する。</summary>
    /// <param name="image">実測パターンの生強度 (width×height)</param>
    /// <param name="geometry">実測画像のピクセルグリッドを基準にした検出器幾何</param>
    /// <param name="reflections">指数付け用の反射リスト (VectorOfG_KikuchiLine)</param>
    /// <param name="waveLength">波長 (nm)。pair-angle シードの幅尤度に使う</param>
    /// <param name="useDictionary">true = MasterPattern 辞書の総当たり ZNCC (Primary indexing)、false = Radon テンプレート照合</param>
    /// <param name="context">動力学 MasterPattern 由来のスナップショット。null なら Radon 単独 (ZNCC 複合ランクと仕上げを行わない)</param>
    /// <param name="properSymmetries">点群 proper 回転 (辞書探索の fundamental-zone 除外用)。null で無効</param>
    public static List<EbsdOrientationCandidate> Run(
        double[] image, int width, int height,
        EbsdDetectorGeometry geometry, Vector3D[] reflections, double waveLength,
        bool useDictionary, EbsdMatchingContext context,
        Matrix3D[] properSymmetries = null, int maxCandidates = 10,
        CancellationToken cancel = default, Action<double> progress = null)
    {
        ArgumentNullException.ThrowIfNull(image);
        ArgumentNullException.ThrowIfNull(geometry);
        ArgumentNullException.ThrowIfNull(reflections);
        if (useDictionary && context == null) throw new ArgumentNullException(nameof(context), "Dictionary search requires the dynamical master pattern.");

        bool refineByZncc = context != null;
        var map = EbsdBandDetector.ComputeRadonMap(image, width, height);
        //260724Cl 追加: 探索エンジン切替 (ラジオボタン、作者指示)。Dictionary = MasterPattern 辞書の総当たり ZNCC (Primary indexing)。
        //候補には Radon z を後付けし、以降の複合ランク+ガード付きトップ精密化は両エンジン共通
        List<EbsdOrientationCandidate> cands;
        if (useDictionary)
        {
            //260724Cl: thoroughCoarse=true (粗段も 96px 完全 robust 総当たり)。作者方針=辞書はパワープレーで精度優先。
            //ベンチ (正しい共通幾何+MC 合成): 3 画像とも辞書トップ=正解系 (14/20・13/15・11/14、5-2_22 では Radon 経路を上回る)
            //260725Cl: properSymmetries (点群 proper 回転の FZ 除外) + 面内分解プロジェクション + SIMD 前処理で
            //12.5s→**2.4〜2.8s/画像** (結果は同一、C2 重複候補も解消)。260725Cl 訂正: 旧コメントの「→4.3s」は中間段階の値
            cands = EbsdDictionaryIndexer.Index(context.MasterPattern, context.PositivePlane, context.NegativePlane, context.Geometry,
                image, width, height, coarseStepDeg: 3, maxCandidates: maxCandidates, thoroughCoarse: true,
                properSymmetries: properSymmetries, cancel: cancel, progress: progress); //260725Ch (progress は 260725Cl)
            //260725Cl (/simplify): 候補ごとの ScoreOrientation はカタログを毎回組み直していた → 一括版で 1 回に (スコアは同一)
            var radonZ = EbsdRadonIndexer.ScoreOrientations(map, geometry, reflections, [.. cands.Select(c => c.Rotation)], SaturateCap);
            for (int i = 0; i < cands.Count; i++)
                cands[i].Score = radonZ[i];
        }
        else
            cands = EbsdRadonIndexer.Index(map, geometry, reflections, waveLength, maxCandidates: maxCandidates,
                saturateCap: refineByZncc ? SaturateCap : 0, cancel: cancel, progress: progress); //260725Ch (progress は 260725Cl)

        if (refineByZncc && cands.Count > 0)
        {
            var projector = new EbsdPatternProjector(context.Geometry, context.RasterWidth, context.RasterHeight);
            var buf = new double[context.RasterWidth * context.RasterHeight];
            var (refRobust, _, _) = EbsdPatternScorer.PrepareReferenceRobust(image, width, height, 160);
            //260727Cl (/simplify): 「投影 → robust ZNCC」の 2 行組が 3 箇所に散っており、前処理を変えるたび 3 箇所同時修正が要る形だったので 1 本にまとめた。
            //  キャンセル判定は元の 3 箇所で意味が違う (候補ループと Nelder-Mead 評価境界だけに置く) ので、ここには入れず呼び出し側に残す。
            double RobustZncc(Matrix3D rot)
            {
                projector.Project(context.MasterPattern, rot, context.PositivePlane, context.NegativePlane, buf);
                return EbsdPatternScorer.Zncc(refRobust, EbsdPatternScorer.RobustPreprocess(buf, context.RasterWidth, context.RasterHeight));
            }
            foreach (var c in cands) //全候補の robust ZNCC (未精密化 — 精密化はどの方位でも ZNCC を伸ばすため判別には使えない)
            {
                cancel.ThrowIfCancellationRequested(); //260725Ch
                c.Zncc = RobustZncc(c.Rotation);
            }
            //候補集合内で ZNCC を標準化 → 複合ランク (Radon の幾何証拠を主、ZNCC は ±2σ クリップの補助)
            double mZ = cands.Average(c => c.Zncc);
            double sZ = Math.Sqrt(Math.Max(cands.Average(c => (c.Zncc - mZ) * (c.Zncc - mZ)), 1E-12));
            cands = [.. cands.OrderByDescending(c => c.Score + ZnccCoef * Math.Clamp((c.Zncc - mZ) / sZ, -2, 2))];
            //複合トップのみ ZNCC 精密化 (±0.25°)。Radon z が GuardMaxZDrop 超劣化する精密化は棄却 (誤収束ガード) //260727Cl: 裸の 0.2 を定数名へ
            var top = cands[0];
            double Score(double[] v)
            {
                cancel.ThrowIfCancellationRequested(); //260725Ch: Nelder-Mead の評価境界で停止
                return -RobustZncc(EbsdIndexer.PerturbRotation(top.Rotation, v[0], v[1], v[2])); //260727Cl
            }
            var (b2, v2, _) = EbsdPatternScorer.NelderMead(Score, [0, 0, 0], [0.25, 0.25, 0.25], 120);
            var rFin = EbsdIndexer.PerturbRotation(top.Rotation, b2[0], b2[1], b2[2]);
            //260725Cl (/simplify): ガードの 2 回採点も一括版へ (旧: ScoreOrientation ×2 でカタログを 2 回構築)
            var guard = EbsdRadonIndexer.ScoreOrientations(map, geometry, reflections, [rFin, top.Rotation], SaturateCap);
            if (guard[0] >= guard[1] - GuardMaxZDrop) //260727Cl: 旧 `- 0.2` (定数化)
            { top.Rotation = rFin; top.Zncc = -v2; }

            //260725Cl 追加 (作者指示): ここまでは「候補の順位付け」のための保守的な微調整 (robust ZNCC を ±0.25° だけ)。
            //順位が確定したあとの最終方位は、Calibrate geometry と同じ目的関数 (素の前処理での ZNCC = context.Reference) と
            //同じステップ (0.7°→0.2°) で仕上げ直す。両者の目的関数が違うと「Find→トップ選択→Calibrate→再び Find」を
            //繰り返したときに方位が 2 値を約 1° で往復して収束しない (作者の実機報告)。順位付けの保護 (±0.25°+ガード) は
            //誤候補を ZNCC で押し上げないためのもので、最終方位の精度のためのものではない、という切り分け。
            double ScoreRaw(double[] v)
            {
                cancel.ThrowIfCancellationRequested();
                projector.Project(context.MasterPattern, EbsdIndexer.PerturbRotation(top.Rotation, v[0], v[1], v[2]), context.PositivePlane, context.NegativePlane, buf);
                return -EbsdPatternScorer.Zncc(context.Reference, buf);
            }
            var (p1, _, _) = EbsdPatternScorer.NelderMead(ScoreRaw, [0, 0, 0], [0.7, 0.7, 0.7], 150);
            //260725Cl 変更: 仕上げステップ 0.2 → OrientationPolishStepDeg (0.1、作者指示)。較正の最終段と同じ値を使う
            var (p2, _, _) = EbsdPatternScorer.NelderMead(ScoreRaw, p1, [OrientationPolishStepDeg, OrientationPolishStepDeg, OrientationPolishStepDeg], 100);
            var rPolished = EbsdIndexer.PerturbRotation(top.Rotation, p2[0], p2[1], p2[2]);
            //仕上げでも同じ誤収束ガード (Radon の幾何証拠を GuardMaxZDrop 超失うなら採用しない) //260727Cl: 裸の 0.2 を定数名へ
            var guardPolish = EbsdRadonIndexer.ScoreOrientations(map, geometry, reflections, [rPolished, top.Rotation], SaturateCap);
            if (guardPolish[0] >= guardPolish[1] - GuardMaxZDrop) //260727Cl: 旧 `- 0.2` (定数化)
            {
                top.Rotation = rPolished;
                //表示中の ZNCC 列は順位付けに使った robust 値なので、仕上げ後の方位で取り直して列と方位の意味を一致させる
                top.Zncc = RobustZncc(top.Rotation); //260727Cl
            }
        }
        #region お蔵入り //260724Cl: 旧 ZNCC 連結 (上位 5 候補を ±1° 精密化して ZNCC 降順に再ランク)。精密化 ZNCC は誤方位ほど伸び正解を落とすため廃止
        //if (refineByZncc && cands.Count > 0)
        //{
        //    var projector = new EbsdPatternProjector(ctx.Geom, ctx.Rw, ctx.Rh);
        //    var buf = new double[ctx.Rw * ctx.Rh];
        //    foreach (var c in cands.Take(5)) //ZNCC は上位 5 候補のみ (1 候補 ~250 評価)
        //    {
        //        double Score(double[] v)
        //        {
        //            projector.Project(ctx.Mp, PerturbRotation(c.Rotation, v[0], v[1], v[2]), ctx.Pos, ctx.Neg, buf);
        //            return -EbsdPatternScorer.Zncc(ctx.Ref, buf);
        //        }
        //        var (b1, _, _) = EbsdPatternScorer.NelderMead(Score, [0, 0, 0], [1.0, 1.0, 1.0], 150);
        //        var (b2, v2, _) = EbsdPatternScorer.NelderMead(Score, b1, [0.25, 0.25, 0.25], 100);
        //        c.Rotation = PerturbRotation(c.Rotation, b2[0], b2[1], b2[2]);
        //        c.Zncc = -v2;
        //    }
        //    cands = [.. cands.OrderByDescending(c => double.IsNaN(c.Zncc) ? double.MinValue : c.Zncc)];
        //}
        #endregion
        return cands;
    }
}
