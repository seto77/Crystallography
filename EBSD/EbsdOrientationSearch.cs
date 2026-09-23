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
    Matrix3D Rotation,
    //260920Cl 追加 (幾何較正の高解像度化・前処理の対称化): 表示中の実測値をフル解像度のまま持つ。
    //  Reference は Find 用に縮小 + 強制背景除算済みで、較正には粗すぎる & シミュレーション側と非対称だった。
    //  DisplayReference はユーザーの「背景を平坦化」設定がそのまま反映された値 (生 or 平坦化後)。
    //  SimFlattenFwhmPx > 0 のとき、較正はシミュレーション側にも同じ半値幅の高域通過を掛けて比べる。
    //  ⚠単位は **DisplayWidth と同じ実測画像 px** (呼び出し側で検出器 px から換算して渡すこと)。
    double[] DisplayReference = null,
    int DisplayWidth = 0,
    int DisplayHeight = 0,
    double SimFlattenFwhmPx = 0);

/// <summary>
/// 260921Cl 追加: 方位探索 (<see cref="EbsdOrientationSearch.Run"/>) の調整パラメータ。既定値は従来のハードコード値と同じ
/// (既定のまま渡せば結果は変わらない)。実測パターンでの調整・比較ハーネス (tools/EbsdIndexTune) と、将来の GUI 設定のため。
/// </summary>
public sealed record EbsdSearchOptions
{
    /// <summary>既定値 (= 従来の動作)</summary>
    public static EbsdSearchOptions Default { get; } = new();

    /// <summary>260922Cl 追加: 検出器幾何が DetZ ±10 mm・DetX/DetY ±3 mm 程度ずれていても方位を見つける構成 (作者指示の条件)。
    /// 検証 (tools/EbsdIndexTune): 実測 Forsterite Ol002 で幾何を最大 (±3, ±3, ±10) mm ずらして全条件 1 位 (0.03°、幾何誤差 ≤ 0.05 mm)。
    /// 合成データ 25 方位 (幾何ずれありランダム) で Si・Botallackite・Al2O3 は 25/25 (従来 2〜5/25)、Forsterite は 21/25。
    /// 時間は Radon 10〜17 s、Dictionary 18〜45 s (従来 0.4 s / 3〜14 s)</summary>
    public static EbsdSearchOptions GeometryRobust { get; } = new()
    {
        ZSearchRangeMm = 12, ZSearchStepMm = 0, XSearchRangeMm = 3, XSearchStepMm = 3,
        CalibrateTopCandidates = 5, CalibrateStages = 1, CalibrateStarts = 5, CalibrateIncludePlain = true,
        DictionaryZHypotheses = 1, DictionaryMergeRadon = true,
        SaturateCap = 4, RadonMaxGridSeeds = 150,
    };

    #region Radon
    /// <summary>Radon 粗探索の方位刻み [°]</summary>
    public double RadonCoarseStepDeg { get; init; } = 3;
    /// <summary>Radon の厳密スコアに使う反射ノード数 (運動学的強度の上位)</summary>
    public int RadonMaxNodes { get; init; } = EbsdRadonIndexer.DefaultMaxNodes;
    /// <summary>Radon 粗探索に使う反射ノード数</summary>
    public int RadonCoarseNodes { get; init; } = 20;
    /// <summary>反射の重み w = I^exp (0.5 は √I を中央値の 3 倍で頭打ち、それ以外は (I/medI)^exp を [0.5, 2] に制限)</summary>
    public double RadonWeightExponent { get; init; } = 0.5;
    /// <summary>260922Cl 追加: Radon 粗探索の格子から精密化へ回すシードの上限 (旧: 52 固定)</summary>
    public int RadonMaxGridSeeds { get; init; } = 52;
    /// <summary>証拠の飽和 cap (ZNCC 複合ランクを使うとき。0 で飽和なし)</summary>
    public double SaturateCap { get; init; } = EbsdOrientationSearch.SaturateCap;
    #endregion

    #region Dictionary
    /// <summary>辞書粗探索の方位刻み [°]</summary>
    public double DictionaryCoarseStepDeg { get; init; } = 3;
    /// <summary>辞書粗探索から精密化へ残す数</summary>
    public int DictionaryCoarseKeep { get; init; } = 64;
    /// <summary>辞書精密化から候補へ残す数</summary>
    public int DictionaryRefineKeep { get; init; } = 12;
    /// <summary>粗段も完全 robust 前処理で総当たりする (遅いが精度優先)</summary>
    public bool DictionaryThoroughCoarse { get; init; } = true;
    #endregion

    #region 複合ランク・仕上げ
    /// <summary>複合ランク combo = Radon z + ZnccCoef · clip(標準化 ZNCC, ±ZnccClip)</summary>
    public double ZnccCoef { get; init; } = 0.5;
    /// <summary>標準化 ZNCC のクリップ幅</summary>
    public double ZnccClip { get; init; } = 2;
    /// <summary>精密化・仕上げの誤収束ガード: Radon z の劣化がこれを超えたら採用しない</summary>
    public double GuardMaxZDrop { get; init; } = 0.2;
    #endregion

    #region 検出器 Z の同時探索 (260922Cl 追加)
    /// <summary>検出器 Z (DetZ) も探す範囲 [mm] (±)。0 で幾何固定 (従来の動作)。
    /// 試料をポールピースぎりぎりまで寄せる EBSD では DetZ が観察ごとに ±10 mm 程度動く (作者談)</summary>
    public double ZSearchRangeMm { get; init; } = 0;
    /// <summary>Z の粗探索の刻み [mm]。0 なら Radon 粗探索の ρ 膨張幅から決める</summary>
    public double ZSearchStepMm { get; init; } = 0;
    /// <summary>精密化で検出器の X (横の平行移動) と Y (カメラ長) も動かす範囲 [mm] (±)。0 で動かさない。
    /// X・Y は粗探索の格子には入れず、Radon の精密化 (Nelder-Mead) の自由度に足すだけ</summary>
    public double XYRefineRangeMm { get; init; } = 0;
    /// <summary>Radon 粗探索の ρ 膨張幅に足す、検出器幾何 (X, Y) の不確かさ [mm]。0 で従来どおり。
    /// 精密化は膨張させない元のマップで行うので、最終的な精度は落ちない (粗探索の生き残りに正解の谷を残すためのもの)</summary>
    public double RadonCoarseGeometryToleranceMm { get; init; } = 0;
    /// <summary>検出器 X (DetX) も Radon 粗探索の格子で探す範囲 [mm] (±)。0 で探さない。Z 探索 (<see cref="ZSearchRangeMm"/>) が有効なときだけ効く。
    /// 格子の刻みは Z と同じで、計算量は X の格子点数倍になる</summary>
    public double XSearchRangeMm { get; init; } = 0;
    /// <summary>X の格子の刻み [mm]。0 なら Z と同じ。膨張マップが格子点の前後を吸収するので、Z より粗くてよい場合がある</summary>
    public double XSearchStepMm { get; init; } = 0;
    /// <summary>複合ランクの上位この数の候補を、それぞれの方位と幾何から軽量に較正 (開始点 1 つ・<see cref="CalibrateStages"/> 段) し、
    /// 較正後の ZNCC で並べ直す。0 で行わない (従来の動作)。
    /// Radon の採点はバンドの中心線の位置しか見ないので、幾何 (特にカメラ長) のずれを方位の回転で補った偽の解と正解が競る。
    /// パターン全体の ZNCC で幾何ごと合わせ直して比べると区別しやすい (実測 Ol002 で X・Y・Z が同時にずれた条件)</summary>
    public int CalibrateTopCandidates { get; init; } = 0;
    /// <summary>候補の較正に使う解像度の段数 (<see cref="EbsdGeometryCalibrator.Run"/> の maxStages)。1 = 160 px、2 = 480 px まで</summary>
    public int CalibrateStages { get; init; } = 1;
    /// <summary>候補の較正の多点開始の点数 (1 = 候補の幾何と方位からだけ)。幾何が X・Y・Z とも数 mm ずれていると 1 点では谷に届かないことがある</summary>
    public int CalibrateStarts { get; init; } = 1;
    /// <summary>較正後の並べ直しに掛ける幾何の事前分布の強さ (0 = なし)。順位 = 較正 ZNCC − w/2·[(ΔZ/σz)² + (ΔX/σxy)² + (ΔY/σxy)²]。
    /// 中心は入力された幾何 (全候補共通)、σz = ZSearchRangeMm、σxy = <see cref="CalibratePriorSigmaXYMm"/></summary>
    public double CalibratePriorWeight { get; init; } = 0;
    /// <summary>事前分布の X・Y の幅 [mm]</summary>
    public double CalibratePriorSigmaXYMm { get; init; } = 3;
    /// <summary>辞書探索で試す Z の仮説の数 (Radon の Z 探索の上位から。Z = 0 は別に必ず試す)</summary>
    public int DictionaryZHypotheses { get; init; } = 3;
    /// <summary>Z の仮説どうし・候補の重複判定で「別の Z」とみなす最小の差 [mm]</summary>
    public double DictionaryZSeparationMm { get; init; } = 1.5;
    /// <summary>辞書探索で Z を探すとき、Z の仮説を立てるのに走らせた Radon の候補も候補の集合に合流させる</summary>
    public bool DictionaryMergeRadon { get; init; } = false;
    /// <summary>幾何を探すとき、幾何を固定した素の Radon 探索も走らせ、その上位も較正に回す (<see cref="CalibrateTopCandidates"/> 件ずつ)。
    /// 幾何の自由度を足すと偽の解が高い点数を取る機会が増え (多重比較)、幾何のずれが無くても素の探索より悪くなる方位があった
    /// (合成 25 方位: 素の探索 23/25 に対し幾何探索のみ 18/25)。素の探索は 0.4 秒程度</summary>
    public bool CalibrateIncludePlain { get; init; } = false;
    #endregion
}

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
    /// <param name="options">260921Cl 追加: 調整パラメータ。null なら <see cref="EbsdSearchOptions.Default"/> (= 従来の動作)</param>
    //旧シグネチャ: public static List<EbsdOrientationCandidate> Run(double[] image, int width, int height, EbsdDetectorGeometry geometry, Vector3D[] reflections, double waveLength,
    //    bool useDictionary, EbsdMatchingContext context, Matrix3D[] properSymmetries = null, int maxCandidates = 10, CancellationToken cancel = default, Action<double> progress = null)
    public static List<EbsdOrientationCandidate> Run(
        double[] image, int width, int height,
        EbsdDetectorGeometry geometry, Vector3D[] reflections, double waveLength,
        bool useDictionary, EbsdMatchingContext context,
        Matrix3D[] properSymmetries = null, int maxCandidates = 10,
        CancellationToken cancel = default, Action<double> progress = null,
        EbsdSearchOptions options = null)
    {
        var o = options ?? EbsdSearchOptions.Default; //260921Cl 追加
        ArgumentNullException.ThrowIfNull(image);
        ArgumentNullException.ThrowIfNull(geometry);
        ArgumentNullException.ThrowIfNull(reflections);
        if (useDictionary && context == null) throw new ArgumentNullException(nameof(context), "Dictionary search requires the dynamical master pattern.");

        //260922Cl 変更 (作者指示「DetZ は観察ごとに ±10 mm 動く」、Codex 相談 260922): 検出器 Z も探せるようにした。
        //  候補ごとに幾何 (EbsdOrientationCandidate.Geometry) を持たせ、Radon 採点・投影・ZNCC・仕上げ・ガードは全部その幾何で行う。
        //  ZSearchRangeMm = 0 (既定) では幾何は 1 つで、処理と結果は旧実装と同じ (下にコメントで残す)。
        //⚠ context.Geometry は geometry と同じ幾何であること (FormEBSD は両方とも BuildDetectorGeometry で作る)。
        //  幾何を持たない候補は、Radon 側は geometry、投影・較正側は context.Geometry で扱う
        bool refineByZncc = context != null;
        bool searchZ = o.ZSearchRangeMm > 0;
        var map = EbsdBandDetector.ComputeRadonMap(image, width, height);
        EbsdDetectorGeometry GeomOf(EbsdOrientationCandidate c) => c.Geometry ?? geometry;
        EbsdMatchingContext CtxOf(EbsdOrientationCandidate c) => c.Geometry == null ? context : context with { Geometry = c.Geometry };
        double[] RadonZ(IReadOnlyList<Matrix3D> rots, EbsdDetectorGeometry g)
            => EbsdRadonIndexer.ScoreOrientations(map, g, reflections, rots, o.SaturateCap, o.RadonWeightExponent, o.RadonMaxNodes);
        List<EbsdOrientationCandidate> RadonIndex(bool withZ, Action<double> prog)
            => EbsdRadonIndexer.Index(map, geometry, reflections, waveLength, maxCandidates: maxCandidates,
                coarseStepDeg: o.RadonCoarseStepDeg, maxNodes: o.RadonMaxNodes, coarseNodes: o.RadonCoarseNodes,
                saturateCap: refineByZncc ? o.SaturateCap : 0, weightExponent: o.RadonWeightExponent, cancel: cancel, progress: prog,
                zSearchRangeMm: withZ ? o.ZSearchRangeMm : 0, zStepMm: o.ZSearchStepMm, xyRefineRangeMm: o.XYRefineRangeMm,
                coarseGeometryToleranceMm: o.RadonCoarseGeometryToleranceMm, xSearchRangeMm: withZ ? o.XSearchRangeMm : 0, xStepMm: o.XSearchStepMm, maxGridSeeds: o.RadonMaxGridSeeds);

        //260724Cl 追加: 探索エンジン切替 (ラジオボタン、作者指示)。Dictionary = MasterPattern 辞書の総当たり ZNCC (Primary indexing)。
        //候補には Radon z を後付けし、以降の複合ランク+ガード付きトップ精密化は両エンジン共通
        List<EbsdOrientationCandidate> cands;
        if (useDictionary)
        {
            //Z を探すときは、先に Radon で Z の仮説を立てる (互いに DictionaryZSeparationMm 以上離れた Z を上位から。0 も保険として必ず含める)。
            //  Radon 首位の Z 1 つに固定すると、Radon が外れたとき辞書も巻き添えになる (Codex 助言)
            //260922Cl: 仮説は Radon 候補の幾何そのもの (X, Y も精密化していればそれも含む)。元の幾何 (ずれ 0) は保険として必ず含める
            var gHyp = new List<EbsdDetectorGeometry> { context.Geometry };
            List<EbsdOrientationCandidate> radonCands = searchZ ? RadonIndex(true, null) : [];
            if (searchZ)
                foreach (var c in radonCands)
                {
                    var gc = GeomOf(c);
                    if (gHyp.All(h => Math.Abs(h.DetZ - gc.DetZ) >= o.DictionaryZSeparationMm)) gHyp.Add(gc);
                    if (gHyp.Count > o.DictionaryZHypotheses) break;
                }
            cands = [];
            foreach (var g in gHyp)
            {
                //260724Cl: thoroughCoarse=true (粗段も 96px 完全 robust 総当たり)。作者方針=辞書はパワープレーで精度優先。
                //260725Cl: properSymmetries (点群 proper 回転の FZ 除外) + 面内分解プロジェクション + SIMD 前処理で 12.5s→2.4〜2.8s/画像
                var list = EbsdDictionaryIndexer.Index(context.MasterPattern, context.PositivePlane, context.NegativePlane, g,
                    image, width, height, coarseStepDeg: o.DictionaryCoarseStepDeg, maxCandidates: maxCandidates,
                    coarseKeep: o.DictionaryCoarseKeep, refineKeep: o.DictionaryRefineKeep, thoroughCoarse: o.DictionaryThoroughCoarse,
                    properSymmetries: properSymmetries, cancel: cancel, progress: progress);
                if (searchZ) foreach (var c in list) { c.Geometry = g; c.GeometryShiftMm = Math.Sqrt(Math.Pow(g.DetX - geometry.DetX, 2) + Math.Pow(g.DetY - geometry.DetY, 2) + Math.Pow(g.DetZ - geometry.DetZ, 2)); }
                var radonZ = RadonZ([.. list.Select(c => c.Rotation)], g);
                for (int i = 0; i < list.Count; i++) list[i].Score = radonZ[i];
                cands.AddRange(list);
            }
            //260922Cl 追加: Z の仮説を立てるために走らせた Radon の候補も合流させる (候補ごとの幾何付き)。辞書の 10 件に正解が入らないと、
            //  Radon が見つけていた正解まで捨ててしまっていた (合成データ: Radon は当たるのに辞書が外す方位が多数)。
            //  以降の複合ランク・較正で辞書の候補と同じ土俵で並べ直す
            if (o.DictionaryMergeRadon) cands.AddRange(radonCands);
        }
        else
            cands = RadonIndex(searchZ, progress); //X・Y の精密化だけの指定 (Z 探索なし) も効く
        //260922Cl 追加: 幾何を固定した素の探索の候補も加える (理由は CalibrateIncludePlain の doc)。較正はグループごとに上位を回す
        var plainSet = new HashSet<EbsdOrientationCandidate>(ReferenceEqualityComparer.Instance);
        if (o.CalibrateIncludePlain && searchZ && refineByZncc)
            foreach (var c in RadonIndex(false, null)) { plainSet.Add(c); cands.Add(c); }

        if (refineByZncc && cands.Count > 0)
        {
            var buf = new double[context.RasterWidth * context.RasterHeight];
            var (refRobust, _, _) = EbsdPatternScorer.PrepareReferenceRobust(image, width, height, 160);
            var projectors = new Dictionary<EbsdDetectorGeometry, EbsdPatternProjector>(ReferenceEqualityComparer.Instance);
            //260727Cl (/simplify): 「投影 → robust ZNCC」を 1 本にまとめたもの。260922Cl: 幾何を引数にとる (候補ごとの幾何で投影する)
            double RobustZncc(Matrix3D rot, EbsdDetectorGeometry g)
            {
                if (!projectors.TryGetValue(g, out var projector))
                    projectors[g] = projector = new EbsdPatternProjector(g, context.RasterWidth, context.RasterHeight);
                projector.Project(context.MasterPattern, rot, context.PositivePlane, context.NegativePlane, buf);
                return EbsdPatternScorer.Zncc(refRobust, EbsdPatternScorer.RobustPreprocess(buf, context.RasterWidth, context.RasterHeight));
            }
            foreach (var c in cands) //全候補の robust ZNCC (未精密化 — 精密化はどの方位でも ZNCC を伸ばすため判別には使えない)
            {
                cancel.ThrowIfCancellationRequested(); //260725Ch
                c.Zncc = RobustZncc(c.Rotation, GeomOf(c));
            }
            //候補集合内で ZNCC を標準化 → 複合ランク (Radon の幾何証拠を主、ZNCC は ±ZnccClip のクリップの補助)。
            //  260922Cl: Z の仮説ごとではなく、統合した全候補で 1 回だけ標準化する (Codex 助言)
            double mZ = cands.Average(c => c.Zncc);
            double sZ = Math.Sqrt(Math.Max(cands.Average(c => (c.Zncc - mZ) * (c.Zncc - mZ)), 1E-12));
            cands = [.. cands.OrderByDescending(c => c.Score + o.ZnccCoef * Math.Clamp((c.Zncc - mZ) / sZ, -o.ZnccClip, o.ZnccClip))];
            if (cands.Count > maxCandidates && o.CalibrateTopCandidates <= 0) //260922Cl: Z の仮説を統合したとき (辞書) は、近い重複を除いてから上限で切る (較正するときは較正のあとで切る)
            {
                var kept = new List<EbsdOrientationCandidate>();
                foreach (var c in cands)
                {
                    if (kept.Any(k => Math.Abs(GeomOf(k).DetZ - GeomOf(c).DetZ) < o.DictionaryZSeparationMm && EbsdIndexer.MisorientationDeg(k.Rotation, c.Rotation, properSymmetries) < 2)) continue;
                    kept.Add(c);
                    if (kept.Count >= maxCandidates) break;
                }
                cands = kept;
            }
            //260922Cl 追加: 上位の候補をそれぞれ軽量に較正 (方位 + 検出器幾何) し、較正後の ZNCC で並べ直す (理由は CalibrateTopCandidates の doc)
            if (o.CalibrateTopCandidates > 0)
            {
                //較正に回す候補: 複合ランクの上位 K 件。素の探索の候補があれば、それとは別に素の探索の上位 K 件も (グループごとに枠を取る)
                var toCalibrate = cands.Where(c => !plainSet.Contains(c)).Take(o.CalibrateTopCandidates)
                    .Concat(cands.Where(plainSet.Contains).Take(o.CalibrateTopCandidates)).ToList();
                int k = toCalibrate.Count;
                double detW = geometry.WidthPx * geometry.PixelSize, detH = geometry.HeightPx * geometry.PixelSize;
                var calibrated = new EbsdOrientationCandidate[k];
                System.Threading.Tasks.Parallel.For(0, k, new System.Threading.Tasks.ParallelOptions { CancellationToken = cancel }, i =>
                {
                    var c = toCalibrate[i];
                    var g0 = GeomOf(c);
                    var res = EbsdGeometryCalibrator.Run(context with { Geometry = g0, Rotation = c.Rotation }, detW, detH, cancel, null,
                        startCount: o.CalibrateStarts, maxStages: o.CalibrateStages);
                    var (gx, gy, gz) = EbsdDetectorGeometry.FromPatternCenter(res.PatternCenterU, res.PatternCenterV, res.CameraLength, g0.DetTilt);
                    var gNew = new EbsdDetectorGeometry(g0.DetTilt, gx, gy, gz, g0.PixelSize, g0.WidthPx, g0.HeightPx, g0.XMirror, g0.SampleTilt);
                    c.Rotation = res.Rotation;
                    c.Geometry = gNew;
                    c.GeometryShiftMm = Math.Sqrt(Math.Pow(gx - geometry.DetX, 2) + Math.Pow(gy - geometry.DetY, 2) + Math.Pow(gz - geometry.DetZ, 2));
                    c.CalibratedZncc = res.Zncc;
                    calibrated[i] = c;
                });
                //幾何の事前分布 (既定オフ): 入力された幾何 (全候補共通) を中心に、Z は広く・X と Y は狭く。候補自身を中心にすると誤シードを守ってしまう (Codex 助言)
                double RankScore(EbsdOrientationCandidate c)
                {
                    if (!(o.CalibratePriorWeight > 0)) return c.CalibratedZncc;
                    var gc = GeomOf(c);
                    double sz = Math.Max(1, o.ZSearchRangeMm), sxy = Math.Max(0.1, o.CalibratePriorSigmaXYMm);
                    double q = Math.Pow((gc.DetZ - geometry.DetZ) / sz, 2) + Math.Pow((gc.DetX - geometry.DetX) / sxy, 2) + Math.Pow((gc.DetY - geometry.DetY) / sxy, 2);
                    return c.CalibratedZncc - o.CalibratePriorWeight * q / 2;
                }
                var ordered = calibrated.OrderByDescending(RankScore).ToList();
                //較正で同じ解に収束した候補は 1 つにまとめる (Codex 指摘)
                var uniq = new List<EbsdOrientationCandidate>();
                foreach (var c in ordered)
                    if (!uniq.Any(u => EbsdIndexer.MisorientationDeg(u.Rotation, c.Rotation, properSymmetries) < 1
                        && Math.Abs(GeomOf(u).DetX - GeomOf(c).DetX) < 0.5 && Math.Abs(GeomOf(u).DetY - GeomOf(c).DetY) < 0.5 && Math.Abs(GeomOf(u).DetZ - GeomOf(c).DetZ) < 0.5))
                        uniq.Add(c);
                //表示する列 (Radon z・robust ZNCC) を較正後の方位と幾何で取り直す (旧値のままだと選んだ行の画像と数値が食い違う。Codex 指摘)
                foreach (var c in uniq)
                {
                    c.Score = RadonZ([c.Rotation], GeomOf(c))[0];
                    c.Zncc = RobustZncc(c.Rotation, GeomOf(c));
                }
                //較正した候補を前に、残りを複合ランクの順で後ろに付け、近い重複を除いて上限で切る
                var rest = cands.Where(c => !toCalibrate.Contains(c));
                var merged = new List<EbsdOrientationCandidate>(uniq);
                foreach (var c in rest)
                {
                    if (merged.Count >= maxCandidates) break;
                    if (merged.Any(m => EbsdIndexer.MisorientationDeg(m.Rotation, c.Rotation, properSymmetries) < 2
                        && Math.Abs(GeomOf(m).DetZ - GeomOf(c).DetZ) < o.DictionaryZSeparationMm)) continue;
                    merged.Add(c);
                }
                cands = merged.Count > maxCandidates ? merged.GetRange(0, maxCandidates) : merged;
            }

            //複合トップのみ ZNCC 精密化 (±0.25°)。Radon z が GuardMaxZDrop 超劣化する精密化は棄却 (誤収束ガード)
            //260922Cl: 候補を較正したときは、方位と幾何を ZNCC で同時に合わせ終えているので、以下の微調整と仕上げは行わない
            //  (実測 Ol002: 較正で正解から 0.03° まで合った方位を、固定幾何の微調整と仕上げが 0.66° ずらしていた)
            if (o.CalibrateTopCandidates > 0) return cands;
            var top = cands[0];
            var topGeom = GeomOf(top);
            double Score(double[] v)
            {
                cancel.ThrowIfCancellationRequested(); //260725Ch: Nelder-Mead の評価境界で停止
                return -RobustZncc(EbsdIndexer.PerturbRotation(top.Rotation, v[0], v[1], v[2]), topGeom);
            }
            var (b2, v2, _) = EbsdPatternScorer.NelderMead(Score, [0, 0, 0], [0.25, 0.25, 0.25], 120);
            var rFin = EbsdIndexer.PerturbRotation(top.Rotation, b2[0], b2[1], b2[2]);
            var guard = RadonZ([rFin, top.Rotation], topGeom);
            if (guard[0] >= guard[1] - o.GuardMaxZDrop)
            { top.Rotation = rFin; top.Zncc = -v2; }

            //260725Cl 追加 (作者指示): 順位が確定したあとの最終方位は、Calibrate geometry と同じ目的関数・同じステップ (0.7°→0.1°) で仕上げ直す
            //  (目的関数が違うと Find と Calibrate の繰り返しで方位が約 1° 往復する)。260920Cl: 較正の最終段 (CreateFinalScorer) を共有。
            //  260922Cl: 仕上げもトップ候補自身の幾何で行う
            var finalScorer = EbsdGeometryCalibrator.CreateFinalScorer(CtxOf(top));
            double ScoreRaw(double[] v)
            {
                cancel.ThrowIfCancellationRequested();
                return -finalScorer.Zncc(EbsdIndexer.PerturbRotation(top.Rotation, v[0], v[1], v[2]));
            }
            var (p1, _, _) = EbsdPatternScorer.NelderMead(ScoreRaw, [0, 0, 0], [0.7, 0.7, 0.7], 150);
            var (p2, _, _) = EbsdPatternScorer.NelderMead(ScoreRaw, p1, [OrientationPolishStepDeg, OrientationPolishStepDeg, OrientationPolishStepDeg], 100);
            var rPolished = EbsdIndexer.PerturbRotation(top.Rotation, p2[0], p2[1], p2[2]);
            //仕上げでも同じ誤収束ガード (Radon の幾何証拠を GuardMaxZDrop 超失うなら採用しない)
            var guardPolish = RadonZ([rPolished, top.Rotation], topGeom);
            if (guardPolish[0] >= guardPolish[1] - o.GuardMaxZDrop)
            {
                top.Rotation = rPolished;
                //表示中の ZNCC 列は順位付けに使った robust 値なので、仕上げ後の方位で取り直して列と方位の意味を一致させる
                top.Zncc = RobustZncc(top.Rotation, topGeom); //260727Cl
            }
        }
        else if (cands.Count > maxCandidates) cands = [.. cands.Take(maxCandidates)];
        #region 旧実装 (260922Cl 変更前: 幾何は 1 つに固定)
        //旧: bool refineByZncc = context != null;
        //旧: var map = EbsdBandDetector.ComputeRadonMap(image, width, height);
        //旧: //260724Cl 追加: 探索エンジン切替 (ラジオボタン、作者指示)。Dictionary = MasterPattern 辞書の総当たり ZNCC (Primary indexing)。
        //旧: //候補には Radon z を後付けし、以降の複合ランク+ガード付きトップ精密化は両エンジン共通
        //旧: List<EbsdOrientationCandidate> cands;
        //旧: if (useDictionary)
        //旧: {
        //旧:     //260724Cl: thoroughCoarse=true (粗段も 96px 完全 robust 総当たり)。作者方針=辞書はパワープレーで精度優先。
        //旧:     //ベンチ (正しい共通幾何+MC 合成): 3 画像とも辞書トップ=正解系 (14/20・13/15・11/14、5-2_22 では Radon 経路を上回る)
        //旧:     //260725Cl: properSymmetries (点群 proper 回転の FZ 除外) + 面内分解プロジェクション + SIMD 前処理で
        //旧:     //12.5s→**2.4〜2.8s/画像** (結果は同一、C2 重複候補も解消)。260725Cl 訂正: 旧コメントの「→4.3s」は中間段階の値
        //旧:     cands = EbsdDictionaryIndexer.Index(context.MasterPattern, context.PositivePlane, context.NegativePlane, context.Geometry,
        //旧:         //image, width, height, coarseStepDeg: 3, maxCandidates: maxCandidates, thoroughCoarse: true, //260921Cl 変更前 (オプション化)
        //旧:         image, width, height, coarseStepDeg: o.DictionaryCoarseStepDeg, maxCandidates: maxCandidates,
        //旧:         coarseKeep: o.DictionaryCoarseKeep, refineKeep: o.DictionaryRefineKeep, thoroughCoarse: o.DictionaryThoroughCoarse,
        //旧:         properSymmetries: properSymmetries, cancel: cancel, progress: progress); //260725Ch (progress は 260725Cl)
        //旧:     //260725Cl (/simplify): 候補ごとの ScoreOrientation はカタログを毎回組み直していた → 一括版で 1 回に (スコアは同一)
        //旧:     //var radonZ = EbsdRadonIndexer.ScoreOrientations(map, geometry, reflections, [.. cands.Select(c => c.Rotation)], SaturateCap); //260921Cl 変更前
        //旧:     var radonZ = EbsdRadonIndexer.ScoreOrientations(map, geometry, reflections, [.. cands.Select(c => c.Rotation)], o.SaturateCap, o.RadonWeightExponent, o.RadonMaxNodes);
        //旧:     for (int i = 0; i < cands.Count; i++)
        //旧:         cands[i].Score = radonZ[i];
        //旧: }
        //旧: else
        //旧:     //cands = EbsdRadonIndexer.Index(map, geometry, reflections, waveLength, maxCandidates: maxCandidates,
        //旧:     //    saturateCap: refineByZncc ? SaturateCap : 0, cancel: cancel, progress: progress); //260725Ch (progress は 260725Cl) //260921Cl 変更前 (オプション化)
        //旧:     cands = EbsdRadonIndexer.Index(map, geometry, reflections, waveLength, maxCandidates: maxCandidates,
        //旧:         coarseStepDeg: o.RadonCoarseStepDeg, maxNodes: o.RadonMaxNodes, coarseNodes: o.RadonCoarseNodes,
        //旧:         saturateCap: refineByZncc ? o.SaturateCap : 0, weightExponent: o.RadonWeightExponent, cancel: cancel, progress: progress);
        //旧:
        //旧: if (refineByZncc && cands.Count > 0)
        //旧: {
        //旧:     var projector = new EbsdPatternProjector(context.Geometry, context.RasterWidth, context.RasterHeight);
        //旧:     var buf = new double[context.RasterWidth * context.RasterHeight];
        //旧:     var (refRobust, _, _) = EbsdPatternScorer.PrepareReferenceRobust(image, width, height, 160);
        //旧:     //260727Cl (/simplify): 「投影 → robust ZNCC」の 2 行組が 3 箇所に散っており、前処理を変えるたび 3 箇所同時修正が要る形だったので 1 本にまとめた。
        //旧:     //  キャンセル判定は元の 3 箇所で意味が違う (候補ループと Nelder-Mead 評価境界だけに置く) ので、ここには入れず呼び出し側に残す。
        //旧:     double RobustZncc(Matrix3D rot)
        //旧:     {
        //旧:         projector.Project(context.MasterPattern, rot, context.PositivePlane, context.NegativePlane, buf);
        //旧:         return EbsdPatternScorer.Zncc(refRobust, EbsdPatternScorer.RobustPreprocess(buf, context.RasterWidth, context.RasterHeight));
        //旧:     }
        //旧:     foreach (var c in cands) //全候補の robust ZNCC (未精密化 — 精密化はどの方位でも ZNCC を伸ばすため判別には使えない)
        //旧:     {
        //旧:         cancel.ThrowIfCancellationRequested(); //260725Ch
        //旧:         c.Zncc = RobustZncc(c.Rotation);
        //旧:     }
        //旧:     //候補集合内で ZNCC を標準化 → 複合ランク (Radon の幾何証拠を主、ZNCC は ±2σ クリップの補助)
        //旧:     double mZ = cands.Average(c => c.Zncc);
        //旧:     double sZ = Math.Sqrt(Math.Max(cands.Average(c => (c.Zncc - mZ) * (c.Zncc - mZ)), 1E-12));
        //旧:     //cands = [.. cands.OrderByDescending(c => c.Score + ZnccCoef * Math.Clamp((c.Zncc - mZ) / sZ, -2, 2))]; //260921Cl 変更前 (オプション化)
        //旧:     cands = [.. cands.OrderByDescending(c => c.Score + o.ZnccCoef * Math.Clamp((c.Zncc - mZ) / sZ, -o.ZnccClip, o.ZnccClip))];
        //旧:     //複合トップのみ ZNCC 精密化 (±0.25°)。Radon z が GuardMaxZDrop 超劣化する精密化は棄却 (誤収束ガード) //260727Cl: 裸の 0.2 を定数名へ
        //旧:     var top = cands[0];
        //旧:     double Score(double[] v)
        //旧:     {
        //旧:         cancel.ThrowIfCancellationRequested(); //260725Ch: Nelder-Mead の評価境界で停止
        //旧:         return -RobustZncc(EbsdIndexer.PerturbRotation(top.Rotation, v[0], v[1], v[2])); //260727Cl
        //旧:     }
        //旧:     var (b2, v2, _) = EbsdPatternScorer.NelderMead(Score, [0, 0, 0], [0.25, 0.25, 0.25], 120);
        //旧:     var rFin = EbsdIndexer.PerturbRotation(top.Rotation, b2[0], b2[1], b2[2]);
        //旧:     //260725Cl (/simplify): ガードの 2 回採点も一括版へ (旧: ScoreOrientation ×2 でカタログを 2 回構築)
        //旧:     //var guard = EbsdRadonIndexer.ScoreOrientations(map, geometry, reflections, [rFin, top.Rotation], SaturateCap); //260921Cl 変更前 (オプション化)
        //旧:     //if (guard[0] >= guard[1] - GuardMaxZDrop) //260727Cl: 旧 `- 0.2` (定数化)
        //旧:     var guard = EbsdRadonIndexer.ScoreOrientations(map, geometry, reflections, [rFin, top.Rotation], o.SaturateCap, o.RadonWeightExponent, o.RadonMaxNodes);
        //旧:     if (guard[0] >= guard[1] - o.GuardMaxZDrop)
        //旧:     { top.Rotation = rFin; top.Zncc = -v2; }
        //旧:
        //旧:     //260725Cl 追加 (作者指示): ここまでは「候補の順位付け」のための保守的な微調整 (robust ZNCC を ±0.25° だけ)。
        //旧:     //順位が確定したあとの最終方位は、Calibrate geometry と同じ目的関数・同じステップ (0.7°→0.1°) で仕上げ直す。
        //旧:     //両者の目的関数が違うと「Find→トップ選択→Calibrate→再び Find」を繰り返したときに方位が 2 値を約 1° で
        //旧:     //往復して収束しない (作者の実機報告)。順位付けの保護 (±0.25°+ガード) は誤候補を ZNCC で押し上げないための
        //旧:     //もので、最終方位の精度のためのものではない、という切り分け。
        //旧:     //260920Cl 変更 (作者指示「Calibrate は ZNCC。Find もこれに統一」): 較正が多段解像度 + ユーザー指定の平坦化 +
        //旧:     //  ZnccSimd へ移ったので、160 px + 強制背景除算の素の Zncc を見ていたこちらとは目的関数がずれていた。
        //旧:     //  較正の最終段そのものを作る EbsdGeometryCalibrator.CreateFinalScorer を通して、定義を 1 つに戻す。
        //旧:     //旧: projector.Project(...160 px...); return -EbsdPatternScorer.Zncc(context.Reference, buf);
        //旧:     var finalScorer = EbsdGeometryCalibrator.CreateFinalScorer(context);
        //旧:     double ScoreRaw(double[] v)
        //旧:     {
        //旧:         cancel.ThrowIfCancellationRequested();
        //旧:         return -finalScorer.Zncc(EbsdIndexer.PerturbRotation(top.Rotation, v[0], v[1], v[2]));
        //旧:     }
        //旧:     var (p1, _, _) = EbsdPatternScorer.NelderMead(ScoreRaw, [0, 0, 0], [0.7, 0.7, 0.7], 150);
        //旧:     //260725Cl 変更: 仕上げステップ 0.2 → OrientationPolishStepDeg (0.1、作者指示)。較正の最終段と同じ値を使う
        //旧:     var (p2, _, _) = EbsdPatternScorer.NelderMead(ScoreRaw, p1, [OrientationPolishStepDeg, OrientationPolishStepDeg, OrientationPolishStepDeg], 100);
        //旧:     var rPolished = EbsdIndexer.PerturbRotation(top.Rotation, p2[0], p2[1], p2[2]);
        //旧:     //仕上げでも同じ誤収束ガード (Radon の幾何証拠を GuardMaxZDrop 超失うなら採用しない) //260727Cl: 裸の 0.2 を定数名へ
        //旧:     //var guardPolish = EbsdRadonIndexer.ScoreOrientations(map, geometry, reflections, [rPolished, top.Rotation], SaturateCap); //260921Cl 変更前 (オプション化)
        //旧:     //if (guardPolish[0] >= guardPolish[1] - GuardMaxZDrop) //260727Cl: 旧 `- 0.2` (定数化)
        //旧:     var guardPolish = EbsdRadonIndexer.ScoreOrientations(map, geometry, reflections, [rPolished, top.Rotation], o.SaturateCap, o.RadonWeightExponent, o.RadonMaxNodes);
        //旧:     if (guardPolish[0] >= guardPolish[1] - o.GuardMaxZDrop)
        //旧:     {
        //旧:         top.Rotation = rPolished;
        //旧:         //表示中の ZNCC 列は順位付けに使った robust 値なので、仕上げ後の方位で取り直して列と方位の意味を一致させる
        //旧:         top.Zncc = RobustZncc(top.Rotation); //260727Cl
        //旧:     }
        //旧: }
        #endregion

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
