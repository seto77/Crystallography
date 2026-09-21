using System;
using System.Collections.Generic;
using System.Threading.Tasks; // 260602Cl: Parallel.For (64 ビンフィットの並列化)
using V3 = OpenTK.Mathematics.Vector3d;

namespace Crystallography;

/// <summary>
/// モンテカルロ BSE シミュレーション結果から、検出器の 8×8 ビン毎に
/// エネルギー・深さの重み分布 w(E, z) = g(E) × h(z|E) をフィッティングし、
/// 任意のピクセル位置で重みを返すクラス。260325Cl 追加
/// </summary>
public sealed class EbsdMonteCarloDistribution
{
    public int BinCount { get; }
    /// <summary>260919Cl 追加 (試行): 電子ごとの蛍光体応答重み φ(E) = max(0, E − E_dead) [keV] の E_dead。NaN なら重み無し (従来 = 1 本 1 票)</summary>
    public double EnergyWeightDeadKeV { get; }

    /// <summary>ビンごとの正規化重み。BinWeights[binI, binJ] は double[energyCount * depthCount]。 // (260327Ch)
    /// <para>⚠ 260921Cl: 深さ写像の改修 (A2) 以降は、<b>ビン中心の射出方向 <see cref="BinCenterMu"/> で経路長へ換算した近似</b>。
    /// 表示合成 (<see cref="EbsdPatternComposer"/>) は画素ごとの μ で <see cref="FillPathLengthWeights"/> を呼ぶのでこの配列を使わない。
    /// 残しているのは全ビン平均を取る大域的な消費者 (A(E) の Ā・ZNCC 用の大域合成・診断ツール) のため。
    /// 密度標本 × 区間幅 (右端則)、各エネルギーで条件付き正規化、ビン内総和 1 (電子の無いビンは 0)。</para></summary>
    public double[,][] BinWeights { get; }

    /// <summary>model 2 用。detector bin の absolute 強度を保った depth-slice 重み。260325Ch 追加
    /// <para>⚠ 260921Cl: <see cref="BinWeights"/> と同じくビン中心の μ で換算した近似。区間質量、各エネルギーで条件付き正規化、ビン内総和 = <see cref="BinFraction"/>。</para></summary>
    public double[,][] BinAbsoluteSliceWeights { get; }

    /// <summary>260921Cl 追加 (深さ写像 A2): ビンごとの当てはめの状態。退避の分岐が 10 nm / 10⁶ nm / 一様重み と散らばっていたのを 1 つにまとめた。</summary>
    public enum BinFitState
    {
        /// <summary>ビンの電子からエネルギー分布と λ(E) を当てはめた</summary>
        Fitted,
        /// <summary>エネルギー分布はビンから、λ(E) は全ビン合算から (結晶内の源が 10 本未満、または λ の有効標本が無い)</summary>
        GlobalLambda,
        /// <summary>ビンの電子が 10 本未満なので、エネルギー分布も λ(E) も全ビン合算から</summary>
        Global,
        /// <summary>全ビン合算でも λ の有効標本が無い。深さ一様 (λ = <see cref="UniformDepthLambdaNm"/>) で代用</summary>
        NoDepthData,
    }

    /// <summary>260921Cl 追加 (深さ写像 A2): 深さの情報が無いときの代用 λ [nm]。事実上の深さ一様 (旧コードの 1E6 と同じ値)。</summary>
    public const double UniformDepthLambdaNm = 1E6;

    /// <summary>260921Cl 追加 (深さ写像 A2): ビンごとの当てはめの状態。</summary>
    public BinFitState[,] BinFitStates { get; }

    /// <summary>260921Cl 追加 (深さ写像 A2): ビンごとの λ_d(E) [nm]。源の<b>垂直深さ</b>を指数分布とみなしたときの平均で、
    /// エネルギー格子ごとに正値化済み ([1, <see cref="UniformDepthLambdaNm"/>] にクランプ)。
    /// 経路長への換算 λ_t = λ_d / μ は使う側が画素ごとに行う (<see cref="FillPathLengthWeights"/>)。</summary>
    public double[,][] BinLambdaNm { get; }

    /// <summary>260921Cl 追加 (深さ写像 A2): ビンごとの射出エネルギー分布 G(E) (エネルギー格子上、総和 1)。電子の無いビンにも全体の分布を入れてある。</summary>
    public double[,][] BinEnergyDistribution { get; }

    /// <summary>260921Cl 追加 (深さ写像 A2): ビンの重み割合 F = ビンの Σφ / 全電子の Σφ (model 2 の絶対スケール)。</summary>
    public double[,] BinFraction { get; }

    /// <summary>260921Cl 追加 (深さ写像 A2): ビン中心の射出方向の μ = cos χ (試料系の z 成分、0 &lt; μ ≤ 1)。</summary>
    public double[,] BinCenterMu { get; }

    /// <summary>260921Cl 追加 (深さ写像 A2): この分布を作ったエネルギー格子 [keV] (コピー)。<see cref="MatchesGridOf"/> が値まで比べるのに使う。</summary>
    public double[] Energies { get; }

    /// <summary>260921Cl 追加 (深さ写像 A2): この分布を作った深さ格子 [nm] (コピー)。<b>出射方向に沿った経路長</b>の格子
    /// (マスターパターンは接球面近似で、方向ごとに等価表面法線 = −出射方向なので)。</summary>
    public double[] Depths { get; }

    //260921Cl 追加 (深さ写像 A2): 合成器のホットループ用の平坦配列 (添字 b = bi·BinCount + bj、エネルギーは b·eLen + ei)
    internal double[] FlatEnergyDistribution { get; }
    internal double[] FlatLambdaNm { get; }
    internal double[] FlatFraction { get; }

    /// <summary>260919Cl 追加: 表面非晶質層の厚さ [nm] (0 = 無し)。輸送は結晶と同じ組成・密度で計算し、層内に源を持つ電子だけを「変調なし」に振り分ける。</summary>
    public double AmorphousLayerNm { get; }
    /// <summary>260919Cl 追加: ビンごとの、干渉性後方散乱源が表面非晶質層内にある電子の割合 (0..1)。合成時に (1−fA)·変調あり + fA·方向平均 として使う。</summary>
    public double[,] BinAmorphousFraction { get; }
    /// <summary>260919Cl 追加: 非晶質層が有効 (厚さ &gt; 0) か。false のとき合成器は平面平均の計算を省く。</summary>
    public bool HasAmorphousLayer => AmorphousLayerNm > 0; // 260919Cl (/simplify) 派生値にする
    /// <summary>260919Cl 追加: 全ビン合算の λ(E) による深さ重み (エネルギー因子は 1、添字 ei*DepthCount+di、model 0/1 用の累積形)。
    /// 非晶質層内の源の基準強度 ⟨M⟩(e) をエネルギーごとに作るのに使い、ビンごとの深さフィットの揺らぎが一様成分のムラにならないようにする。非晶質層が無い (厚さ 0) ときは null。</summary>
    public double[] GlobalDepthWeights { get; }
    /// <summary>260919Cl 追加: 同上の depth-slice 区間質量形 (model 2 用)。非晶質層が無いときは null。</summary>
    public double[] GlobalDepthSliceWeights { get; }

    /// <summary>この分布を作った MasterPattern の energy 数。260727Cl 追加</summary>
    public int EnergyCount { get; }

    /// <summary>この分布を作った MasterPattern の depth 数。260727Cl 追加:
    /// 重み配列の添字は wIdx = ei * DepthCount + di なので、消費側の MasterPattern が別の格子だと
    /// 例外を出さずに別スライスの重みを引く。消費側はこの 2 値で格子一致を検査すること
    /// (<see cref="MatchesGridOf"/>)。</summary>
    public int DepthCount { get; }

    /// <summary>重み配列の (energy × depth) 格子が <paramref name="mp"/> と一致するか。260727Cl 追加
    /// <para>260921Cl 変更 (深さ写像 A2): 長さだけでなく<b>値</b>まで比べる。重みは格子の値 (深さ t と区間幅) から作るので、
    /// 同じ点数でも値の違う格子を通すと、例外も警告も無く別の深さの重みで合成してしまう。
    /// 許容差は相対 1e-9 (UI の開始・刻みから作り直した格子の足し算誤差は 1e-14 程度)。</para></summary>
    //旧: public bool MatchesGridOf(MasterPattern mp)
    //旧:     => mp != null && mp.Energies.Length == EnergyCount && mp.Depths.Length == DepthCount;
    public bool MatchesGridOf(MasterPattern mp)
        => mp != null && SameGrid(mp.Energies, Energies) && SameGrid(mp.Depths, Depths);

    /// <summary>260921Cl 追加: 2 つの格子が同じ長さで、各値が相対 1e-9 以内で一致するか。</summary>
    static bool SameGrid(double[] a, double[] b)
    {
        if (a == null || b == null || a.Length != b.Length) return false;
        for (int i = 0; i < a.Length; i++)
            if (!(Math.Abs(a[i] - b[i]) <= 1E-9 * Math.Max(1.0, Math.Abs(b[i])))) return false; //NaN もここで落ちる
        return true;
    }

    /// <summary>
    /// MasterPattern の全 (energy, depth) スライスを、この MC 分布の全ビン平均重みで微分合成した 1 枚 (pos/neg 半球) を返す。260724Cl 追加。
    /// EbsdPatternComposer の表示合成 (ApplyWeightedModel2 = absolute MC × differential master) のグローバル近似で、位置依存重みを全ビンで平均している。 //260727Cl: 移設に伴い参照先を訂正
    /// ⚠ 260921Cl: ビニングが検出器から射出半球へ変わったので、この平均は「検出器に当たる電子」ではなく「射出半球全体」の平均になった (ctor の doc)。
    /// ZNCC 方位照合 (複合ランク・幾何較正) 用の実測に忠実なシミュレーションパターン。単一スライスより実測との相関が上がることをベンチで確認済み。
    /// </summary>
    //260920Cl シグネチャ変更 (作者指示): 損失依存のコントラスト係数 A(E) = exp(−(E0 − E)/E_c) を掛けられるようにした。
    //  表示合成 (EbsdPatternComposer.CoherenceLossDecayKeV) と同じ重みにしないと、ZNCC の目的関数だけ別のパターンを見ることになる。
    //  ⚠ここでは (1 − A) の平坦な台座を足していない。ZNCC は画像全体への定数加算・定数倍に不変なので、方向に依らない成分は結果を変えないため。
    //旧: public (float[] Pos, float[] Neg) ComposeGlobalWeightedPattern(MasterPattern mp)
    //260921Cl シグネチャ変更 (深さ写像 A2、段階 4): 大域の重みを外から渡せるようにした。
    //  EbsdPatternComposer.ComputeDetectorAverageSliceWeights (検出器画像の画素平均) を渡すのが本命。
    //  null なら従来どおり全ビンの BinAbsoluteSliceWeights (ビン中心の μ で換算) を射出半球全体で足したもの
    //  (A2 以降は検出器が見ない表面すれすれの方向の長い経路まで混ざるので、検出器があるときは渡すこと)。
    //旧: public (float[] Pos, float[] Neg) ComposeGlobalWeightedPattern(MasterPattern mp, double beamEnergyKeV = double.NaN, double coherenceLossDecayKeV = double.NaN)
    public (float[] Pos, float[] Neg) ComposeGlobalWeightedPattern(MasterPattern mp, double beamEnergyKeV = double.NaN, double coherenceLossDecayKeV = double.NaN,
        double[] sliceWeights = null)
    {
        ArgumentNullException.ThrowIfNull(mp); //260725Ch: null を後段の不明瞭な参照例外にしない
        if (mp.GridSize < 2) throw new ArgumentException("MasterPattern.GridSize must be at least 2.", nameof(mp)); //260725Ch
        //260727Cl 追加: 重み配列は (EnergyCount × DepthCount) 格子で添字 wIdx = ei*DepthCount + di を前提にしている。
        //  別格子の mp を渡されたとき、下の `k < bw.Length` クランプは範囲外アクセスこそ防ぐが、ei/di の対応が
        //  ずれた重みを黙って積むだけだった (ZNCC 照合・幾何較正が静かに劣化する)。入口で弾く。
        if (!MatchesGridOf(mp))
            throw new ArgumentException(
                $"The MC distribution grid ({EnergyCount} energies x {DepthCount} depths) does not match the MasterPattern ({mp.Energies.Length} x {mp.Depths.Length}).", nameof(mp));
        int eLen = mp.Energies.Length, dLen = mp.Depths.Length;
        var wG = new double[eLen * dLen];
        if (sliceWeights != null) //260921Cl 追加: 外から渡された大域の重み (検出器平均)
        {
            if (sliceWeights.Length != wG.Length) throw new ArgumentException($"sliceWeights must have {wG.Length} elements (energies x depths).", nameof(sliceWeights));
            Array.Copy(sliceWeights, wG, wG.Length);
        }
        else
        for (int bi = 0; bi < BinCount; bi++)
            for (int bj = 0; bj < BinCount; bj++)
            {
                var bw = BinAbsoluteSliceWeights[bi, bj];
                if (bw == null) continue;
                for (int k = 0; k < wG.Length && k < bw.Length; k++) wG[k] += bw[k];
            }
        //260920Cl 追加: エネルギースライスごとに A(E) を掛ける。直後に総和で正規化するので、重みの規約 (総和 1) は保たれる
        //260921Cl 変更: 式を EbsdPatternComposer.CoherenceFactors に一本化 (表示合成と同じ重みであることが機能の前提なので、
        //  同じ式を 2 クラスに置かない)。全部が 1 のときは null が返るが、直後に正規化するので掛けても掛けなくても同じ
        var cohA = EbsdPatternComposer.CoherenceFactors(mp.Energies, beamEnergyKeV, coherenceLossDecayKeV);
        if (cohA != null)
            for (int ei = 0; ei < eLen; ei++)
            {
                double a = cohA[ei];
                for (int di = 0; di < dLen; di++) wG[ei * dLen + di] *= a;
            }
        double wSum = 0;
        foreach (var v in wG) wSum += v;
        if (wSum > 0) for (int k = 0; k < wG.Length; k++) wG[k] /= wSum;

        //int gs2 = mp.GridSize * mp.GridSize; //260725Ch 変更前
        int gs2 = checked(mp.GridSize * mp.GridSize); //260725Ch
        var pos = new float[gs2];
        var neg = new float[gs2];
        /*260725Ch 変更前: energy×depth×hemisphere の各スライスで Parallel.ForEach を起動して同期していた
        //260725Cl (/simplify): 画素ループを並列化 (UI スレッド上で走るため。grid 512×320 スライスで 5000 万反復あり体感フリーズになる)。
        //並列軸は画素 i のみ — スライス順の加算順序は不変なので結果はビット一致。旧: 逐次 for + 半球 2 要素の配列を毎スライス確保
        void Accumulate(float[] dst, MasterPattern.Hemisphere hemi, int ei, int di, double wgt)
        {
            var p = mp.GetPlane(hemi, ei, di);
            if (p == null || p.Length == 0) return;
            var pPrev = di > 0 ? mp.GetPlane(hemi, ei, di - 1) : null;
            bool hasPrev = pPrev is { Length: > 0 };
            System.Threading.Tasks.Parallel.ForEach(System.Collections.Concurrent.Partitioner.Create(0, gs2), range =>
            {
                for (int i = range.Item1; i < range.Item2; i++)
                    dst[i] += (float)(wgt * Math.Max(0, p[i] - (hasPrev ? pPrev[i] : 0))); //model 2 の differential 合成
            });
        }
        for (int ei = 0; ei < eLen; ei++)
            for (int di = 0; di < dLen; di++)
            {
                double wgt = wG[ei * dLen + di];
                if (wgt < 1E-15) continue;
                Accumulate(pos, MasterPattern.Hemisphere.PositiveZ, ei, di, wgt);
                Accumulate(neg, MasterPattern.Hemisphere.NegativeZ, ei, di, wgt);
            }
        */
        //260725Ch: 有効スライス参照を先に集め、半球ごとに Parallel.ForEach を 1 回だけ起動する。
        //各画素内のスライス加算順と float 丸めは旧実装と同じなので結果はビット一致する。
        var posSlices = new List<(float[] Plane, float[] Previous, double Weight)>();
        var negSlices = new List<(float[] Plane, float[] Previous, double Weight)>();
        //260726Cl 追加 (正本 §1.4): plane は厚さ t までの累積 M(t) なので、隣接差は区間積分。
        //区間平均 R̄ = ΔM/Δt にするため区間幅で割る。MC 側の重みは区間の「質量」なので割らない (質量 × 平均応答が正しい寄与)。
        //等間隔グリッドでは全体が定数倍になるだけだが、不等間隔では区間ごとの重み比が変わる
        var depthWidths = mp.DepthIntervals;
        for (int ei = 0; ei < eLen; ei++)
            for (int di = 0; di < dLen; di++)
            {
                double wgt = wG[ei * dLen + di];
                if (wgt < 1E-15) continue;
                wgt /= depthWidths[di]; //260726Cl
                var p = mp.GetPlane(MasterPattern.Hemisphere.PositiveZ, ei, di);
                //if (p is { Length: > 0 }) //260725Ch 変更前: 短い非空 plane は並列画素ループ内で範囲外になった
                if (p != null && p.Length >= gs2) //260725Ch
                {
                    var previous = di > 0 ? mp.GetPlane(MasterPattern.Hemisphere.PositiveZ, ei, di - 1) : null;
                    posSlices.Add((p, previous != null && previous.Length >= gs2 ? previous : null, wgt));
                }
                p = mp.GetPlane(MasterPattern.Hemisphere.NegativeZ, ei, di);
                if (p != null && p.Length >= gs2) //260725Ch
                {
                    var previous = di > 0 ? mp.GetPlane(MasterPattern.Hemisphere.NegativeZ, ei, di - 1) : null;
                    negSlices.Add((p, previous != null && previous.Length >= gs2 ? previous : null, wgt));
                }
            }

        Accumulate(pos, posSlices);
        Accumulate(neg, negSlices);
        return (pos, neg);

        static void Accumulate(float[] destination, List<(float[] Plane, float[] Previous, double Weight)> sliceList)
        {
            if (sliceList.Count == 0) return;
            var slices = sliceList.ToArray();
            System.Threading.Tasks.Parallel.ForEach(System.Collections.Concurrent.Partitioner.Create(0, destination.Length), range =>
            {
                for (int i = range.Item1; i < range.Item2; i++)
                {
                    float value = destination[i];
                    for (int si = 0; si < slices.Length; si++)
                    {
                        var (plane, previous, weight) = slices[si];
                        value += (float)(weight * Math.Max(0, plane[i] - (previous == null ? 0 : previous[i])));
                    }
                    destination[i] = value;
                }
            });
        }
    }

    // 260718Cl: smpTilt 引数を削除。BSE の Vec は MC 内で試料傾斜を織り込んで lab 座標系で追跡され、検出器 (detY/detZ/detTilt も lab 座標系) への投影に試料傾斜は不要 (Codex 検証済: 適用すると二重計上)。
    // 260723Cl: 円形検出器 (半径 detR) → 矩形検出器 (半幅 halfWidth × 半高 halfHeight) + 中心 X オフセット (detX) へ変更。
    // 旧シグネチャ (〜260921Cl): EbsdMonteCarloDistribution(bseList, beamEnergy, detTilt, detX, detY, detZ, halfWidth, halfHeight, energies, depths, binCount = 8, ...)
    /// <summary>260921Cl シグネチャ変更 (作者指示): <b>検出器ではなく、試料表面側の半球全体をビニングする。</b>
    /// <para>【なぜ】旧版は検出器に当たった電子だけを検出器面上の 8×8 に切り、外れた電子は捨てていた。そのため
    /// 検出器より広い視野を表示すると、そこは全部「端のビンの外挿」で物理的な裏付けが無かった
    /// (実測: 視野 683 mm のとき検出器はその 1/10。正本 §1.5.4)。半球を切れば視野をどれだけ広げても実データで埋まる。</para>
    /// <para>【格子】MasterPattern と同じ <b>Rosca-Lambert 等積正方格子</b>。等積なので全ビンの立体角が等しく、
    /// 極にも縁にも特異点が無い。既定 16×16 は、検出器 8×8 (立体角 ≈ 1.7 sr を 64 分割) と同程度の角度分解能で
    /// 半球 (2π sr) を覆う数 (2π ÷ (1.7/64) ≈ 236 ≈ 15.4²)。</para>
    /// <para>【座標系】bseList.Vec は lab 系なので X 軸まわりに −sampleTilt 回して試料系へ戻し、射出半球 (z &gt; 0) を切る。
    /// 合成側 (<see cref="EbsdPatternComposer"/>) はまったく同じ写像を使う。⚠ 合成側の画素方向は射出方向の
    /// **逆向き**なので、あちらでは符号を反転してからこの写像へ渡す (両者の式を突き合わせて数値で確認済)。</para>
    /// <para>⚠ 検出器に依存しなくなったので、検出器幾何を変えても MC 分布を作り直す必要はない。</para>
    /// <para>⚠ <see cref="ComposeGlobalWeightedPattern"/> の全ビン平均は「検出器に当たる電子の平均」から
    /// 「射出半球全体の平均」へ意味が変わる。ZNCC 目的関数と E_c 校正 (正本 §2.7) はこの重みを使うので、
    /// 旧値へ厳密に戻したいときは検出器の立体角で重み付けし直す必要がある。</para></summary>
    public EbsdMonteCarloDistribution(
        (double Depth, V3 Vec, double Energy)[] bseList,
        double beamEnergy,
        double sampleTilt, //260921Cl 変更: detTilt/detX/detY/detZ/halfWidth/halfHeight を置き換え
        double[] energies, double[] depths,
        int binCount = 16, //260921Cl 変更 (旧 8): 半球を検出器 8×8 と同程度の角度分解能で覆う数
        double amorphousLayerNm = 0, // 260919Cl 追加: 表面非晶質層の厚さ [nm]。層内の源は変調なし成分、結晶内の深さは層厚を引く
        double energyWeightDeadKeV = double.NaN) // 260919Cl 追加 (試行): 蛍光体応答 φ(E)=max(0,E−E_dead) で電子を重み付け。NaN = 重み無し (従来)
    {
        //260725Ch: 下流のビン補間は非空の energy/depth 軸を前提とする (260921Cl: 補間は 4×4 の 3 次 B スプライン。添字は端でクランプするので binCount ≥ 2 でよい)
        ArgumentNullException.ThrowIfNull(bseList);
        ArgumentNullException.ThrowIfNull(energies);
        ArgumentNullException.ThrowIfNull(depths);
        if (binCount < 2) throw new ArgumentOutOfRangeException(nameof(binCount), "binCount must be at least 2.");
        if (energies.Length == 0) throw new ArgumentException("At least one energy is required.", nameof(energies));
        if (depths.Length == 0) throw new ArgumentException("At least one depth is required.", nameof(depths));
        if (!double.IsFinite(sampleTilt)) throw new ArgumentOutOfRangeException(nameof(sampleTilt)); //260921Cl 変更 (旧: halfWidth/halfHeight の検証)

        BinCount = binCount;
        //260920Cl (/simplify2): E_dead がビームエネルギー以上だと全電子の φ が 0 になり、パターンが例外も警告も無く恒等的にゼロになる。
        //  低加速 (5〜10 keV) の EBSD で E_dead を上限 10 keV にすると実際に起こるので、重み無しへ落とす
        bool useEnergyWeight = double.IsFinite(energyWeightDeadKeV) && energyWeightDeadKeV < beamEnergy; // 260919Cl 追加 / 260920Cl ガード追加
        EnergyWeightDeadKeV = useEnergyWeight ? energyWeightDeadKeV : double.NaN;
        double totalWeight = 0; // 260919Cl 追加: Σφ(E) (重み無しなら電子数)

        //260921Cl 変更: 検出器面との交点ではなく、試料系での射出方向そのものでビニングする
        //  旧: var (sinDet, cosDet) = Math.SinCos(detTilt); double lamDenom = detY * sinDet - detZ * cosDet;
        var (sinSmp, cosSmp) = Math.SinCos(sampleTilt);

        var bins = new List<(double depth, double energy)>[binCount, binCount];
        for (int i = 0; i < binCount; i++)
            for (int j = 0; j < binCount; j++)
                bins[i, j] = new List<(double, double)>();

        // var binTotals = new int[binCount, binCount]; // 260919Cl 変更前 (エネルギー重み試行で double 化)
        // var amorphousCounts = new int[binCount, binCount]; // 260919Cl 変更前
        var binTotals = new double[binCount, binCount]; // 260919Cl 変更: 非晶質層内の源も含めたビンの重み和 Σφ (重み無しなら電子数)
        var binCounts = new int[binCount, binCount]; // 260919Cl 追加: ビンの電子数 (一様フォールバックの判定用)
        var amorphousCounts = new double[binCount, binCount]; // 260919Cl 変更: 層内の源の重み和
        var binEnergies = new List<double>[binCount, binCount]; // 260919Cl 追加: 非晶質層内の源も含む全電子のエネルギー (エネルギー分布のフィットに使う)
        var binEnergyWeights = new List<double>[binCount, binCount]; // 260919Cl 追加: 同じ並びの重み φ(E)
        for (int i = 0; i < binCount; i++) for (int j = 0; j < binCount; j++) { binEnergies[i, j] = new List<double>(); binEnergyWeights[i, j] = new List<double>(); }
        var allCrystalline = new List<(double depth, double energy)>(); // 260919Cl 追加: 結晶内の源が少ないビンの λ(E) フォールバック用 (全ビン合算)
        if (!(amorphousLayerNm > 0) || !double.IsFinite(amorphousLayerNm)) amorphousLayerNm = 0; // 260919Cl 追加
        foreach (var (depth, vec, energy) in bseList)
        {
            //260920Cl (/simplify2): 重みの母数は全電子。検出器を外れる continue より前で積む (旧 binFraction の分母 bseList.Length と同義)
            totalWeight += useEnergyWeight ? Math.Max(0, energy - energyWeightDeadKeV) : 1.0;
            //260921Cl 変更: 検出器面との交点 (px, py) で 8×8 に切っていたのをやめ、試料系の射出方向を
            //  Rosca-Lambert 等積正方格子へ写して半球を切る。旧コードは git 履歴 (260723Cl 版) を参照。
            //  lab → 試料系: X 軸まわりに −sampleTilt。RotationX(−θ): y′ = y cosθ + z sinθ、z′ = −y sinθ + z cosθ
            //lab → 試料系 = EbsdDetectorGeometry.LabToSample と同じ回転 (X 軸まわり、y′ = c·y + s·z、z′ = −s·y + c·z)。
            //  (/simplify: 一時的に公開ヘルパーにしていたが、呼び出しがここ 1 か所だけで同名の既存メソッドと紛らわしいので戻した)
            double sy = vec.Y * cosSmp + vec.Z * sinSmp, sz = -vec.Y * sinSmp + vec.Z * cosSmp;
            if (!(sz > 0)) continue; //試料内部へ向かう方向 (物理的に射出しない)。NaN もここで落ちる
            //⚠ 合成側とまったく同じ写像を使う (DirectionToBinCoords の doc)
            var (fbx, fby) = DirectionToBinCoords(vec.X, sy, sz, binCount);
            int bi = Math.Clamp((int)Math.Round(fbx), 0, binCount - 1);
            int bj = Math.Clamp((int)Math.Round(fby), 0, binCount - 1);
            // binTotals[bi, bj]++; binEnergies[bi, bj].Add(energy); // 260919Cl 変更前
            double phi = useEnergyWeight ? Math.Max(0, energy - energyWeightDeadKeV) : 1.0; // 260919Cl 追加 (試行): 蛍光体の発光量 ∝ E − E_dead
            binTotals[bi, bj] += phi; binCounts[bi, bj]++; // 260919Cl 変更
            //260920Cl (/simplify2): totalWeight はループ先頭で全電子ぶん積む。ここで積むと母数が「検出器に当たった電子」になり、
            //  重み OFF (φ=1) でも旧 bseList.Length と一致しなくなっていた
            binEnergies[bi, bj].Add(energy); binEnergyWeights[bi, bj].Add(phi); // 260919Cl 追加
            // if (depth < amorphousLayerNm) { amorphousCounts[bi, bj]++; continue; } // 260919Cl 変更前
            if (depth < amorphousLayerNm) { amorphousCounts[bi, bj] += phi; continue; } // 260919Cl 追加: 層内の源は菊池変調を持たない一様成分として数える (重み付き)
            // bins[bi, bj].Add((depth, energy)); // 260919Cl 変更前
            bins[bi, bj].Add((depth - amorphousLayerNm, energy)); // 260919Cl 変更: 結晶内の深さは結晶表面 (層の底) から測る
            allCrystalline.Add((depth - amorphousLayerNm, energy)); // 260919Cl 追加
        }

        BinWeights = new double[binCount, binCount][];
        BinAbsoluteSliceWeights = new double[binCount, binCount][]; // (260325Ch) model 2 用
        AmorphousLayerNm = amorphousLayerNm; BinAmorphousFraction = new double[binCount, binCount]; // 260919Cl 追加
        int eLen = energies.Length, dLen = depths.Length;
        EnergyCount = eLen; DepthCount = dLen; // 260727Cl 追加: 消費側が格子一致を検査できるようにする

        //260921Cl 追加 (深さ写像 A2): ビンの当てはめは「パラメータ (G(E), λ_d(E), F)」として持ち、重みは経路長へ換算してから作る。
        //  【なぜ】MC の源深さ d は傾斜試料表面からの**垂直深さ**だが、マスターパターンの深さ格子 t は接球面近似の
        //  **出射方向に沿った経路長**で、t = d/μ (μ = cos χ)。旧版は d をそのまま t として使い、経路長を
        //  最大 4 倍 (検出器下端 χ ≈ 75°) 過小評価していた (正本: .project-guidance/ReciPro/ReciPro_EBSD深さ写像_改修計画.md)。
        //  λ はビンごとに垂直深さで当てはめ、μ での換算は使う側が画素ごとに行う (A2)。ここではビン中心の μ で換算した
        //  互換用の重み配列 (BinWeights / BinAbsoluteSliceWeights) も作る。
        Energies = [.. energies]; Depths = [.. depths]; //呼び出し側が配列を書き換えても分布の意味が変わらないようにコピー
        var depthWidths = MasterPattern.ComputeDepthIntervals(depths);
        BinFitStates = new BinFitState[binCount, binCount];
        BinLambdaNm = new double[binCount, binCount][];
        BinEnergyDistribution = new double[binCount, binCount][];
        BinFraction = new double[binCount, binCount];
        BinCenterMu = new double[binCount, binCount];
        int nBins = binCount * binCount;
        FlatEnergyDistribution = new double[nBins * eLen];
        FlatLambdaNm = new double[nBins * eLen];
        FlatFraction = new double[nBins];

        // 260602Cl 変更: 64 ビンは互いに独立 (distinct な BinWeights[bi,bj]/BinAbsoluteSliceWeights[bi,bj] へ書く) なので
        //   Parallel.For 化。bins への集約 (上の foreach) は逐次のまま、ここはフィット段だけ並列化する。
        // 260919Cl 追加: 全ビン合算の λ(E) 多項式。結晶内の源が 10 本未満のビン (厚い非晶質層・検出器端) で使う
        //260921Cl 変更 (深さ写像 A2): 合算の当てはめが「有効な標本ゼロ」で失敗したときは 10 nm を黙って使わず、深さ一様 (NoDepthData) と明示する
        //旧: bool hasGlobalLambda = allCrystalline.Count >= 10;
        //旧: double gla = 1E6, glb = 0, glc = 0; // λ→∞ = 深さ一様 (合算でも足りないときのフォールバック)
        //旧: if (hasGlobalLambda) FitLambdaFromData(allCrystalline, energies, out gla, out glb, out glc);
        double gla = UniformDepthLambdaNm, glb = 0, glc = 0; // λ→∞ = 深さ一様 (合算でも足りないときのフォールバック)
        bool hasGlobalLambda = allCrystalline.Count >= 10 && FitLambdaFromData(allCrystalline, energies, out gla, out glb, out glc);
        if (!hasGlobalLambda) (gla, glb, glc) = (UniformDepthLambdaNm, 0, 0);
        var globalLambda = EvaluateLambda(energies, gla, glb, glc);
        //260921Cl 追加 (深さ写像 A2): 全電子のエネルギー分布 (電子 10 本未満のビンの退避先。旧版はそのビンを (E, 深さ) 一様にしていた)
        var energyGroups = new List<double>[nBins]; var weightGroups = useEnergyWeight ? new List<double>[nBins] : null;
        for (int b = 0; b < nBins; b++) { energyGroups[b] = binEnergies[b / binCount, b % binCount]; if (weightGroups != null) weightGroups[b] = binEnergyWeights[b / binCount, b % binCount]; }
        var globalEnergy = NormalizeToUnitSum(ComputeEnergyGaussian(energyGroups, energies, weightGroups));
        // 260919Cl 追加: 合算 λ(E) の深さ重み (エネルギー因子 1)。合成器で非晶質源の基準強度 ⟨M⟩(e) を作るのに使う
        //260921Cl 変更 (深さ写像 A2): 全ビン合算の λ ではなく、各ビンの (μ_b で換算した) 条件付き深さ分布を F_b で混ぜたものにする (下の Parallel.For の後)。
        //  全 (E, 深さ) を一括正規化してから混ぜると、有限の深さ範囲で捕まえられる割合がビンごとに違うせいでビン間の比が変わる (Codex 指摘)
        //旧: if (amorphousLayerNm > 0) // (/simplify) 非晶質層が無ければ使われないので作らない
        //旧: {
        //旧:     var ones = new double[eLen]; Array.Fill(ones, 1.0);
        //旧:     GlobalDepthWeights = new double[eLen * dLen]; GlobalDepthSliceWeights = new double[eLen * dLen];
        //旧:     FillBinWeights(GlobalDepthWeights, energies, depths, ones, gla, glb, glc, useSliceMass: false, totalScale: 1.0);
        //旧:     FillBinWeights(GlobalDepthSliceWeights, energies, depths, ones, gla, glb, glc, useSliceMass: true, totalScale: 1.0);
        //旧: }
        Parallel.For(0, binCount * binCount, (int idx) =>
        {
            int bi = idx / binCount, bj = idx % binCount;
            var binData = bins[bi, bj];
            var weights = new double[eLen * dLen];
            var absoluteSliceWeights = new double[eLen * dLen]; // (260325Ch) model 2 用
            // double binFraction = bseList.Length > 0 ? (double)binData.Count / bseList.Length : 0.0; // (260325Ch) 260919Cl 変更前
            // int binTotal = binTotals[bi, bj]; double binFraction = bseList.Length > 0 ? (double)binTotal / bseList.Length : 0.0; // 260919Cl 変更前
            double binTotal = binTotals[bi, bj]; // 260919Cl 変更: 重み和
            double binFraction = totalWeight > 0 ? binTotal / totalWeight : 0.0; // 260919Cl 変更: 重みの総量はビンの全電子の Σφ (非晶質分は合成時に fA で振り分ける)
            BinAmorphousFraction[bi, bj] = binTotal > 0 ? amorphousCounts[bi, bj] / binTotal : 0.0; // 260919Cl 追加

            //260921Cl 変更 (深さ写像 A2): 一様フォールバックをやめ、パラメータを当てはめる (足りなければ全ビン合算へ退避) → ビン中心の μ で重みにする
            //旧: if (binCounts[bi, bj] < 10) // 260919Cl 変更: ビンの全電子が 10 本未満のときだけ一様フォールバック (判定は本数、重みではない)
            //旧: {
            //旧:     double uniform = binTotal > 0 ? 1.0 / (eLen * dLen) : 0.0; // 260919Cl 変更: binData.Count → binTotal (全電子が層内でも重みを持つ)
            //旧:     Array.Fill(weights, uniform);
            //旧:     double absoluteUniform = binTotal > 0 ? binFraction / (eLen * dLen) : 0.0; // (260325Ch) 260919Cl binData.Count → binTotal
            //旧:     Array.Fill(absoluteSliceWeights, absoluteUniform);
            //旧: }
            //旧: else
            //旧: {
            //旧:     FitBinDistribution(binData, binEnergies[bi, bj], useEnergyWeight ? binEnergyWeights[bi, bj] : null, energies, depths, weights, absoluteSliceWeights, binFraction, hasGlobalLambda, gla, glb, glc); // 260919Cl 変更: エネルギーは全電子 (重み付き)、λ は結晶内の源 (少なければ全ビン合算)
            //旧: }
            double[] energyDistribution, lambda;
            BinFitState state;
            if (binCounts[bi, bj] < 10)
            {
                energyDistribution = globalEnergy; lambda = globalLambda;
                state = hasGlobalLambda ? BinFitState.Global : BinFitState.NoDepthData;
            }
            else
            {
                energyDistribution = NormalizeToUnitSum(ComputeEnergyGaussian([binEnergies[bi, bj]], energies, useEnergyWeight ? [binEnergyWeights[bi, bj]] : null)); // 260919Cl: エネルギーは全電子 (重み付き)
                if (binData.Count >= 10 && FitLambdaFromData(binData, energies, out double la, out double lb, out double lc)) // 260919Cl: λ は結晶内の源
                {
                    lambda = EvaluateLambda(energies, la, lb, lc); state = BinFitState.Fitted;
                }
                else
                {
                    lambda = globalLambda; state = hasGlobalLambda ? BinFitState.GlobalLambda : BinFitState.NoDepthData;
                }
            }
            double mu = LambertBinCenterMu(bi, bj, binCount);
            BinFitStates[bi, bj] = state; BinEnergyDistribution[bi, bj] = energyDistribution; BinLambdaNm[bi, bj] = lambda;
            BinFraction[bi, bj] = binFraction; BinCenterMu[bi, bj] = mu;
            Array.Copy(energyDistribution, 0, FlatEnergyDistribution, idx * eLen, eLen);
            Array.Copy(lambda, 0, FlatLambdaNm, idx * eLen, eLen);
            FlatFraction[idx] = binFraction;

            //互換用の重み (ビン中心 μ の近似)。電子の無いビンは旧版どおり 0
            if (binTotal > 0)
            {
                var scaled = new double[eLen];
                FillPathLengthWeights(weights, energyDistribution, lambda, mu, depths, depthWidths, sliceMass: false); //総和 1
                for (int ei = 0; ei < eLen; ei++) scaled[ei] = binFraction * energyDistribution[ei];
                FillPathLengthWeights(absoluteSliceWeights, scaled, lambda, mu, depths, depthWidths, sliceMass: true); //総和 F
            }

            BinWeights[bi, bj] = weights;
            BinAbsoluteSliceWeights[bi, bj] = absoluteSliceWeights; // (260325Ch)
        });

        //260921Cl 追加 (深さ写像 A2): 非晶質層の基準強度用の深さ重み = Σ_b F_b·G_b(E)·p_b(t|E) (ビン中心 μ で換算)。
        //  消費側 (EbsdPatternComposer.CollapseToEnergyReference) がエネルギーごとに正規化するので、結果は
        //  「エネルギー E の電子全体 (全ビンを電子数で混ぜたもの) の、経路長の条件付き分布」になる。
        //  G_b(E) がビンごとに違うので、F_b だけで混ぜたものとは違う (Codex 指摘でコメントを訂正)
        if (amorphousLayerNm > 0) // (/simplify) 非晶質層が無ければ使われないので作らない
        {
            GlobalDepthWeights = new double[eLen * dLen]; GlobalDepthSliceWeights = new double[eLen * dLen];
            for (int bi = 0; bi < binCount; bi++)
                for (int bj = 0; bj < binCount; bj++)
                {
                    double f = BinFraction[bi, bj];
                    if (!(f > 0)) continue;
                    var bw = BinWeights[bi, bj]; var bs = BinAbsoluteSliceWeights[bi, bj];
                    for (int k = 0; k < GlobalDepthWeights.Length; k++) { GlobalDepthWeights[k] += f * bw[k]; GlobalDepthSliceWeights[k] += bs[k]; }
                }
        }
    }

    /// <summary>260921Cl 追加 (深さ写像 A2): 源の深さ分布を<b>出射方向の経路長</b>へ換算して、(エネルギー × 深さ) の重みを w へ書く。
    /// <para>垂直深さ d を平均 λ_d(E) の指数分布とみなすと、経路長 t = d/μ も指数分布で平均 λ_t = λ_d/μ。
    /// 減衰率 α = μ/λ_d なので、<b>指数の引数だけでなく密度の係数も μ 倍になる</b> (ここでは各エネルギーで正規化するので係数は消える)。</para>
    /// <para>・<paramref name="sliceMass"/> = true (model 2): 区間質量 P(t_{d−1} &lt; t &lt; t_d) = e^{−α t_{d−1}}·(1 − e^{−α Δt_d})。
    ///   Δt が小さい・μ → 0 のときの桁落ちを避けるため、α Δt が小さいときは (1 − e^{−x}) を級数で評価する。</para>
    /// <para>・false (model 0/1): 密度標本 × 区間幅 α e^{−α t_d}·Δt_d (右端則。不等間隔の格子でも区間の重みが正しくなる)。</para>
    /// <para>どちらも<b>エネルギーごとに条件付き正規化</b>してから <paramref name="energyWeight"/>[e] を掛ける。
    /// つまり有限の深さ範囲 (最後の格子点 T) を超える尾部は同じエネルギーの中で配り直され、エネルギー分布は変えない
    /// (旧版は全 (E, 深さ) を一括で正規化していたので、尾部の長い低エネルギーの比が格子の打ち切りで削れていた)。
    /// 結果の総和は Σ_e energyWeight[e]。</para>
    /// <para>表示合成 (<see cref="EbsdPatternComposer"/>) は画素の μ で、互換用の <see cref="BinWeights"/> 等はビン中心の μ でこれを呼ぶ。</para></summary>
    /// <param name="w">出力 (長さ eLen·dLen、添字 ei·dLen + di)</param>
    /// <param name="energyWeight">エネルギーごとの重み (長さ eLen、非負)</param>
    /// <param name="lambdaNm">垂直深さの平均 λ_d(E) [nm] (長さ eLen、正)</param>
    /// <param name="mu">出射方向の cos χ (0 &lt; μ ≤ 1)</param>
    /// <param name="depths">経路長の格子 t_d [nm] (単調増加、t₀ = 0 は暗黙)</param>
    /// <param name="depthWidths">区間幅 Δt_d (<see cref="MasterPattern.ComputeDepthIntervals"/>)</param>
    /// <param name="sliceMass">true = 区間質量 (model 2)、false = 密度 × 区間幅 (model 0/1)</param>
    public static void FillPathLengthWeights(Span<double> w, ReadOnlySpan<double> energyWeight, ReadOnlySpan<double> lambdaNm,
        double mu, ReadOnlySpan<double> depths, ReadOnlySpan<double> depthWidths, bool sliceMass)
    {
        int eLen = energyWeight.Length, dLen = depths.Length;
        for (int ei = 0; ei < eLen; ei++)
        {
            var row = w.Slice(ei * dLen, dLen);
            double ew = energyWeight[ei];
            if (!(ew > 0)) { row.Clear(); continue; }
            double alpha = mu / lambdaNm[ei]; //経路長の減衰率 [1/nm] (λ_t = λ_d/μ)
            double sum = 0;
            if (sliceMass)
            {
                double tPrev = 0, ePrev = 1; //t₀ = 0、e^{−α·0} = 1
                for (int di = 0; di < dLen; di++)
                {
                    double t = depths[di], x = alpha * (t - tPrev);
                    double eCur = Math.Exp(-alpha * t);
                    //1 − e^{−x} の桁落ち対策: x が小さいときは級数 x − x²/2 + x³/6 − x⁴/24 (打ち切り誤差は相対 x⁴/120、x = 1e-4 で 1e-18)。
                    //  x > 1e-4 の差分 e_prev − e_cur は相対 1e-16/x ≲ 1e-12 の丸めで済む
                    double m = x > 1E-4 ? ePrev - eCur : ePrev * x * (1 - x * (0.5 - x * (1.0 / 6 - x / 24)));
                    row[di] = m; sum += m;
                    tPrev = t; ePrev = eCur;
                }
            }
            else
            {
                //e^{−α t} を先頭の点で割って、α t が大きいときのアンダーフローで全部 0 になるのを防ぐ (正規化で消える)
                double t0 = depths[0];
                for (int di = 0; di < dLen; di++)
                {
                    double m = Math.Exp(-alpha * (depths[di] - t0)) * depthWidths[di];
                    row[di] = m; sum += m;
                }
            }
            if (sum > 0 && double.IsFinite(sum))
            {
                double s = ew / sum;
                for (int di = 0; di < dLen; di++) row[di] *= s;
            }
            else { row.Clear(); row[0] = ew; } //到達しないはず (α が NaN など)。質量は捨てずに最浅スライスへ
        }
    }

    /// <summary>260921Cl 追加 (深さ写像 A2): λ(E) = la + lb·E + lc·E² をエネルギー格子で評価し、[1, UniformDepthLambdaNm] nm に正値化する (旧 FillBinWeights の下限 1 nm と同じ規約)。</summary>
    static double[] EvaluateLambda(double[] energies, double la, double lb, double lc)
    {
        var l = new double[energies.Length];
        for (int ei = 0; ei < l.Length; ei++)
        {
            double e = energies[ei], v = la + lb * e + lc * e * e;
            l[ei] = double.IsFinite(v) ? Math.Clamp(v, 1.0, UniformDepthLambdaNm) : UniformDepthLambdaNm;
        }
        return l;
    }

    /// <summary>260921Cl 追加: 総和 1 に正規化した新しい配列を返す (総和が 0 なら一様)。</summary>
    static double[] NormalizeToUnitSum(double[] v)
    {
        double s = 0; foreach (var x in v) s += x;
        var o = new double[v.Length];
        if (s > 0 && double.IsFinite(s)) for (int i = 0; i < v.Length; i++) o[i] = v[i] / s;
        else Array.Fill(o, 1.0 / Math.Max(1, v.Length));
        return o;
    }

    /// <summary>260921Cl 追加 (深さ写像 A2): ビン (bi, bj) の中心の射出方向の μ = cos χ。<see cref="DirectionToBinCoords"/> の逆 (ビン中心が整数座標)。</summary>
    public static double LambertBinCenterMu(int bi, int bj, int binCount)
    {
        double scale = binCount / (2.0 * MasterPattern.SquareLimit);
        double la = (bi + 0.5) / scale - MasterPattern.SquareLimit, lb = MasterPattern.SquareLimit - (bj + 0.5) / scale;
        var v = MasterPattern.RoscaLambertToSphereSquare(la, lb, MasterPattern.Hemisphere.PositiveZ);
        double len = Math.Sqrt(v.X * v.X + v.Y * v.Y + v.Z * v.Z);
        return len > 0 ? Math.Clamp(v.Z / len, 1E-6, 1.0) : 1.0;
    }

    // 260602Cl 変更: weights と absoluteSliceWeights を 1 回の共通 fit から両方埋める形へ統合。
    //   旧実装は同一 binData に対して本メソッドを 2 回呼び (useSliceMass=false/true)、
    //   meanE/sigmaL,R/gE[]/lambdaPerEnergy/FitLambdaPolynomial を二重計算していた。
    //   これらは useSliceMass/totalScale に依存しない共通部 (Stage1) なので 1 回だけ計算し、
    //   weights の書き込み (Stage2) のみ 2 種類行う。物理的な値は旧 2 回呼びと厳密一致。
    // 旧シグネチャ: FitBinDistribution(data, beamEnergy, energies, depths, weights, bool useSliceMass=false, double totalScale=1.0)
    // private static void FitBinDistribution(List<(double depth, double energy)> data, double beamEnergy, double[] energies, double[] depths, double[] weights, double[] absoluteSliceWeights, double binFraction) // 260919Cl 変更前のシグネチャ
    //260921Cl 削除 (深さ写像 A2): ctor が当てはめ (ComputeEnergyGaussian / FitLambdaFromData) と重み作り (FillPathLengthWeights) を直接呼ぶ形にしたので未使用。
    ///// <summary>260919Cl 変更: エネルギー分布は非晶質層内の源も含む全電子 (allEnergies) から、λ(E) は結晶内の源 (data) から作る。
    ///// 結晶内の源が 10 本未満なら全ビン合算の λ(E) (hasGlobalLambda) を使い、それも無ければ深さ一様 (λ→∞)。
    ///// 旧 Stage1 の本体は ComputeEnergyGaussian / FitLambdaFromData へ切り出した (中身は同じ)。</summary>
    //private static void FitBinDistribution(
    //    List<(double depth, double energy)> data, List<double> allEnergies, List<double> allWeights, // 260919Cl 追加: allWeights (null なら等重み)
    //    double[] energies, double[] depths, // 260919Cl (/simplify2) 未使用だった beamEnergy を除去
    //    double[] weights,             // useSliceMass=false, totalScale=1.0 相当
    //    double[] absoluteSliceWeights, // useSliceMass=true,  totalScale=binFraction 相当
    //    double binFraction, bool hasGlobalLambda, double gla, double glb, double glc)
    //{
    //    if (allEnergies.Count == 0) return;
    //    var gE = ComputeEnergyGaussian(allEnergies, energies, allWeights); // 260919Cl 変更: 重み付き
    //    double la = gla, lb = glb, lc = glc;
    //    if (data.Count >= 10) FitLambdaFromData(data, energies, out la, out lb, out lc);
    //    else if (!hasGlobalLambda) { la = 1E6; lb = 0; lc = 0; }
    //    FillBinWeights(weights, energies, depths, gE, la, lb, lc, useSliceMass: false, totalScale: 1.0);
    //    FillBinWeights(absoluteSliceWeights, energies, depths, gE, la, lb, lc, useSliceMass: true, totalScale: binFraction);
    //}

    /// <summary>260919Cl 追加 (旧 FitBinDistribution Stage1 前半): 電子のエネルギー分布を左右非対称ガウシアンで近似し、energies 格子上の重み gE を返す。
    /// <para>260921Cl 変更 (深さ写像 A2): 複数の電子リスト (グループ) をまとめて当てはめられるようにした。全ビン合算の分布 (疎なビンの退避先) を
    /// 1 本のリストへ連結せずに作るため (1000 万本の MC で 160 MB の一時配列を避ける)。ビン単体は 1 要素のグループで呼ぶ。</para></summary>
    // private static double[] ComputeEnergyGaussian(List<double> electronEnergies, double[] energies) // 260919Cl 変更前のシグネチャ
    //旧: private static double[] ComputeEnergyGaussian(List<double> electronEnergies, double[] energies, List<double> electronWeights = null) // 260919Cl 変更: 重み付き平均・片側分散 (weights=null なら従来と同値)
    private static double[] ComputeEnergyGaussian(IReadOnlyList<List<double>> energyGroups, double[] energies, IReadOnlyList<List<double>> weightGroups = null) // 260921Cl 変更: グループ化
    {
        int eLen = energies.Length;
        double sumE = 0, sumW = 0, firstE = double.NaN;
        for (int g = 0; g < energyGroups.Count; g++)
        {
            var es = energyGroups[g]; var ws = weightGroups?[g];
            for (int i = 0; i < es.Count; i++) { double w = ws == null ? 1.0 : ws[i]; sumE += w * es[i]; sumW += w; if (double.IsNaN(firstE)) firstE = es[i]; }
        }
        if (double.IsNaN(firstE)) return new double[eLen]; //電子が 1 本も無い (呼び出し側は総和 0 を一様へ正規化する)
        double meanE = sumW > 0 ? sumE / sumW : firstE;

        double varL = 0, varR = 0, wL = 0, wR = 0;
        int nL = 0, nR = 0;
        for (int g = 0; g < energyGroups.Count; g++)
        {
            var es = energyGroups[g]; var ws = weightGroups?[g];
            for (int i = 0; i < es.Count; i++)
            {
                double e = es[i], w = ws == null ? 1.0 : ws[i];
                if (e < meanE) { varL += w * (e - meanE) * (e - meanE); wL += w; nL++; }
                else { varR += w * (e - meanE) * (e - meanE); wR += w; nR++; }
            }
        }
        //旧 (単一リスト版): int count = electronEnergies.Count; double sumE = 0, sumW = 0;
        //旧: for (int i = 0; i < count; i++) { double w = electronWeights == null ? 1.0 : electronWeights[i]; sumE += w * electronEnergies[i]; sumW += w; }
        //旧: double meanE = sumW > 0 ? sumE / sumW : electronEnergies[0];
        //旧: double varL = 0, varR = 0, wL = 0, wR = 0; int nL = 0, nR = 0;
        //旧: for (int i = 0; i < count; i++)
        //旧: {
        //旧:     double e = electronEnergies[i], w = electronWeights == null ? 1.0 : electronWeights[i];
        //旧:     if (e < meanE) { varL += w * (e - meanE) * (e - meanE); wL += w; nL++; }
        //旧:     else { varR += w * (e - meanE) * (e - meanE); wR += w; nR++; }
        //旧: }
        double sigmaL = nL > 1 && wL > 0 ? Math.Sqrt(varL / wL) : 0.5;
        double sigmaR = nR > 1 && wR > 0 ? Math.Sqrt(varR / wR) : 0.5;
        double Ep = meanE;
        if (sigmaL < 0.01) sigmaL = 0.5;
        if (sigmaR < 0.01) sigmaR = 0.5;

        // 260327Cl: 2*sigma*sigma を事前計算して Exp 内の除算を軽減
        double inv2SigmaSqL = 1.0 / (2 * sigmaL * sigmaL);
        double inv2SigmaSqR = 1.0 / (2 * sigmaR * sigmaR);

        var gE = new double[eLen];
        for (int ei = 0; ei < eLen; ei++)
        {
            double dE = energies[ei] - Ep;
            double inv2SigmaSq = dE < 0 ? inv2SigmaSqL : inv2SigmaSqR;
            gE[ei] = Math.Exp(-dE * dE * inv2SigmaSq);
        }

        return gE;
    }

    /// <summary>260919Cl 追加 (旧 FitBinDistribution Stage1 後半): エネルギー格子ごとの平均深さから λ(E) の 2 次多項式 la + lb·E + lc·E² を作る。
    /// <para>260921Cl 変更 (深さ写像 A2): 当てはめの成否を返す。有効なエネルギー窓 (電子 4 本以上) が 1 つも無いとき false
    /// (旧版はこのとき黙って λ = 10 nm を返していた)。深さはすべて<b>垂直深さ</b>のまま扱う (経路長への換算は使う側)。</para></summary>
    //旧: private static void FitLambdaFromData(List<(double depth, double energy)> data, double[] energies, out double la, out double lb, out double lc)
    private static bool FitLambdaFromData(List<(double depth, double energy)> data, double[] energies, out double la, out double lb, out double lc)
    {
        int eLen = energies.Length;
        double eStep = eLen > 1 ? Math.Abs(energies[0] - energies[^1]) / (eLen - 1) : 1.0;
        double halfStep = eStep * 0.5, e0 = energies[0];

        // 260602Cl 変更: 旧実装は ei ごとに binData 全体を energy 窓でフルスキャン (O(eLen × N_bin))。
        //   energies が ComputeGridFromRanges 由来 (beamEnergy から降順・ほぼ等間隔・0.1 丸め) の場合は、
        //   各電子を 1 パスで energy 窓へ割り当てる高速経路 (O(N_bin)) を使う。候補 index cand を算出し、
        //   現行と同じ半開窓 [energies[ei]-eStep/2, energies[ei]+eStep/2) を満たす ei (丸めズレ・窓の
        //   微小な重なりを cand±2 で吸収) すべてに加算 → 旧ループと同一集合・同一加算順で bit 一致。
        //   constructor は energies を外部受け取りなので、降順・ほぼ等間隔 (|δ'|<=eStep/2) でない、または
        //   eStep<=0 のときは安全のため旧 eLen フルスキャンへ fallback する (任意 energies で常に正しい)。
        var depthSumPerEnergy = new double[eLen];
        var depthCountPerEnergy = new int[eLen];

        bool uniformDescending = eStep > 0;
        if (uniformDescending)
            for (int ei = 0; ei < eLen; ei++)
                if (Math.Abs(energies[ei] - (e0 - ei * eStep)) > halfStep) { uniformDescending = false; break; } // cand±2 が全マッチ窓を覆う前提を保証

        if (uniformDescending)
        {
            double invEStep = 1.0 / eStep;
            foreach (var (depth, energy) in data)
            {
                int cand = (int)((e0 - energy) * invEStep); // 降順 energies に対する floor 候補 (energy<=e0)
                for (int ei = cand - 2; ei <= cand + 2; ei++)
                    if ((uint)ei < (uint)eLen && energy >= energies[ei] - halfStep && energy < energies[ei] + halfStep)
                    {
                        depthSumPerEnergy[ei] += depth;
                        depthCountPerEnergy[ei]++;
                    }
            }
        }
        else
        {
            // fallback: 旧 eLen フルスキャン (任意 energies に対し常に正しい)
            for (int ei = 0; ei < eLen; ei++)
            {
                double eLow = energies[ei] - halfStep, eHigh = energies[ei] + halfStep;
                foreach (var (depth, energy) in data)
                    if (energy >= eLow && energy < eHigh)
                    {
                        depthSumPerEnergy[ei] += depth;
                        depthCountPerEnergy[ei]++;
                    }
            }
        }

        var lambdaPerEnergy = new double[eLen];
        for (int ei = 0; ei < eLen; ei++)
            lambdaPerEnergy[ei] = depthCountPerEnergy[ei] > 3 ? Math.Max(1.0, depthSumPerEnergy[ei] / depthCountPerEnergy[ei]) : -1;

        //旧: FitLambdaPolynomial(energies, lambdaPerEnergy, out la, out lb, out lc);
        return FitLambdaPolynomial(energies, lambdaPerEnergy, out la, out lb, out lc); //260921Cl 変更: 成否を返す
    }

    //260921Cl 削除 (深さ写像 A2): 重みは FillPathLengthWeights (経路長へ換算・エネルギーごとの条件付き正規化) で作るので未使用
    ///// <summary>
    ///// 260602Cl 追加: 共通 fit パラメータ (gE, lambda 多項式 la/lb/lc) から 1 つの weights 配列を埋める。
    ///// useSliceMass=false: 連続深さ重み g(E)·exp(-z/λ)/λ。
    ///// useSliceMass=true : depth slice 区間質量 g(E)·(exp(-z_prev/λ) - exp(-z/λ))。
    ///// 旧 FitBinDistribution の weights 書き込み部 (Stage2) をそのまま切り出したもの。(260919Cl: 分割時に誤って消した summary を復元)
    ///// </summary>
    //private static void FillBinWeights(
    //    double[] weights, double[] energies, double[] depths, double[] gE,
    //    double la, double lb, double lc, bool useSliceMass, double totalScale)
    //{
    //    int eLen = energies.Length, dLen = depths.Length;
    //    double totalWeight = 0;
    //    for (int ei = 0; ei < eLen; ei++)
    //    {
    //        double lambda = la + lb * energies[ei] + lc * energies[ei] * energies[ei];
    //        if (lambda < 1.0) lambda = 1.0;
    //        double invLambda = 1.0 / lambda; // 260327Cl: 除算を事前計算
    //
    //        // 260327Cl: useSliceMass 時、隣接スライスで Exp 値を再利用
    //        if (useSliceMass)
    //        {
    //            double expPrev = 1.0; // di==0 の lowerDepth=0 → exp(0)=1
    //            for (int di = 0; di < dLen; di++)
    //            {
    //                double expCur = Math.Exp(-depths[di] * invLambda);
    //                double w = gE[ei] * (expPrev - expCur); // (260325Ch) model 2 は depth slice 区間質量
    //                weights[ei * dLen + di] = w;
    //                totalWeight += w;
    //                expPrev = expCur; // 260327Cl: 次スライスの lowerDepth 用にキャッシュ
    //            }
    //        }
    //        else
    //        {
    //            for (int di = 0; di < dLen; di++)
    //            {
    //                double w = gE[ei] * Math.Exp(-depths[di] * invLambda) * invLambda;
    //                weights[ei * dLen + di] = w;
    //                totalWeight += w;
    //            }
    //        }
    //    }
    //
    //    if (totalWeight > 0)
    //        for (int k = 0; k < weights.Length; k++)
    //            weights[k] = weights[k] / totalWeight * totalScale; // (260325Ch)
    //}

    // 260327Cl: List<double> を固定長配列に置き換えて GC 圧力を軽減
    //260921Cl 変更 (深さ写像 A2): 成否を返す。有効な点が 1 つも無いときは false (旧版は λ = 10 nm という根拠の無い値を返していた)
    //旧: private static void FitLambdaPolynomial(double[] energies, double[] lambdaValues, out double a, out double b, out double c)
    private static bool FitLambdaPolynomial(double[] energies, double[] lambdaValues, out double a, out double b, out double c)
    {
        int validCount = 0;
        for (int i = 0; i < energies.Length; i++)
            if (lambdaValues[i] > 0) validCount++;

        if (validCount == 0)
        {
            //旧: a = 10; b = 0; c = 0;
            //旧: return;
            a = UniformDepthLambdaNm; b = 0; c = 0;
            return false;
        }

        var validE = new double[validCount];
        var validL = new double[validCount];
        int idx = 0;
        for (int i = 0; i < energies.Length; i++)
            if (lambdaValues[i] > 0)
            {
                validE[idx] = energies[i];
                validL[idx] = lambdaValues[i];
                idx++;
            }

        if (validCount == 1)
        {
            a = validL[0]; b = 0; c = 0;
            return true; //260921Cl 旧: return;
        }
        if (validCount == 2)
        {
            double e0 = validE[0], e1 = validE[1], l0 = validL[0], l1 = validL[1];
            b = (l1 - l0) / (e1 - e0);
            a = l0 - b * e0;
            c = 0;
            return true; //260921Cl 旧: return;
        }

        int n = validCount;
        double s0 = n, s1 = 0, s2 = 0, s3 = 0, s4 = 0;
        double r0 = 0, r1 = 0, r2 = 0;
        for (int k = 0; k < n; k++)
        {
            double e = validE[k], l = validL[k];
            double e2 = e * e;
            s1 += e; s2 += e2; s3 += e * e2; s4 += e2 * e2;
            r0 += l; r1 += l * e; r2 += l * e2;
        }

        double det = s0 * (s2 * s4 - s3 * s3) - s1 * (s1 * s4 - s3 * s2) + s2 * (s1 * s3 - s2 * s2);
        if (Math.Abs(det) < 1e-30)
        {
            a = r0 / s0; b = 0; c = 0;
            return true; //260921Cl 旧: return;
        }

        double invDet = 1.0 / det;
        a = ((s2 * s4 - s3 * s3) * r0 + (s2 * s3 - s1 * s4) * r1 + (s1 * s3 - s2 * s2) * r2) * invDet;
        b = ((s2 * s3 - s1 * s4) * r0 + (s0 * s4 - s2 * s2) * r1 + (s1 * s2 - s0 * s3) * r2) * invDet;
        c = ((s1 * s3 - s2 * s2) * r0 + (s1 * s2 - s0 * s3) * r1 + (s0 * s2 - s1 * s1) * r2) * invDet;
        return true; //260921Cl 追加
    }

    // public static (double energyLoss80, double depth99) ComputeRangesFromMC( // 260919Cl 変更前: エネルギー損失は 80 パーセンタイル
    public static (double energyLoss95, double depth99) ComputeRangesFromMC( // 260919Cl 変更: エネルギー損失は 95 パーセンタイルまで
        (double Depth, V3 Vec, double Energy)[] bseList,
        double beamEnergy)
    {
        if (bseList == null || bseList.Length == 0)
            return (5.0, 50.0);

        int n = bseList.Length;
        var losses = new double[n];
        var depths = new double[n];
        for (int i = 0; i < n; i++)
        {
            losses[i] = beamEnergy - bseList[i].Energy;
            depths[i] = bseList[i].Depth;
        }

        // int idxLoss80 = Math.Min((int)(n * 0.8), n - 1); // 260919Cl 変更前
        int idxLoss95 = Math.Min((int)(n * 0.95), n - 1); // 260919Cl 変更: 95 パーセンタイル
        int idxDepth99 = Math.Min((int)(n * 0.99), n - 1); // (260326Ch)

        // 260602Cl 変更: 80%/99% パーセンタイルのためのフルソート 2 回を QuickSelect (nth_element) に置換。
        //   QuickSelect.Execute(span, k, cmp) は span[k] に sorted index k の値を置く (前は <=, 後は >=) ので、
        //   Array.Sort 後の losses[idxLoss80]/depths[idxDepth99] と同一値。O(n log n) → 平均 O(n)。
        // Array.Sort(losses); Array.Sort(depths); // 260602Cl 変更前
        QuickSelect.Execute(losses.AsSpan(), idxLoss95, static (a, b) => a.CompareTo(b));
        QuickSelect.Execute(depths.AsSpan(), idxDepth99, static (a, b) => a.CompareTo(b));

        return (losses[idxLoss95], depths[idxDepth99]); // 260919Cl idxLoss80 → idxLoss95
    }

    /// <summary>260921Cl 追加 (作者指示): <b>試料系の射出方向 → ビン座標</b>。ビン中心が整数、範囲は [−0.5, binCount−0.5]。
    /// <para>⚠ <b>分布を作る側 (このクラスの ctor) と使う側 (<see cref="EbsdPatternComposer"/>) は必ずこれを共有する。</b>
    /// 同じ式を 2 か所に書くと、片方だけ直したときに「絵はそれらしいが微妙にずれている」という見つけにくい壊れ方をする。</para>
    /// <para>格子は MasterPattern と同じ Rosca-Lambert 等積正方格子。b (縦) は MasterPattern の行の向きに合わせて反転する。</para></summary>
    public static (double bx, double by) DirectionToBinCoords(double sx, double sy, double sz, int binCount)
    {
        var (la, lb) = MasterPattern.SphereToRoscaLambertSquare(sx, sy, sz);
        double scale = binCount / (2.0 * MasterPattern.SquareLimit);
        return ((la + MasterPattern.SquareLimit) * scale - 0.5, (MasterPattern.SquareLimit - lb) * scale - 0.5);
    }

    /// <summary>260921Cl 追加 (深さ写像 A2): 不等間隔の深さ格子の既定の形状パラメータ β (<see cref="ComputeGridFromRanges"/>)。
    /// Si 20 kV (Si004 の幾何) で T ≈ 170 nm のとき、最浅の区間 ≈ 1.4 nm (旧・等間隔の 1.56 nm と同程度)、最深 ≈ 10 nm。</summary>
    public const double DefaultDepthGridBeta = 2.0;

    /// <summary>260921Cl 追加 (深さ写像 A2): 検出器を渡さないとき、経路長の上限 T の統計に入れる射出方向の下限 μ (= cos 78.5°)。</summary>
    public const double MinMuWithoutDetector = 0.2;

    /// <summary>260921Cl 追加 (深さ写像 A2): マスターパターンの深さ格子の上限 T [nm] = <b>経路長</b> t = (d − a)/μ の <paramref name="quantile"/> 分位点。
    /// <para>マスターパターンは接球面近似なので深さ格子は出射方向に沿った経路長。MC の源深さ d は垂直深さなので、
    /// 旧版のように垂直深さの 99 % 点を上限にすると、斜めに出る電子の経路が格子に収まらない
    /// (Si 20 kV・Si004 の幾何で、垂直深さの 99 % = 62 nm に対し、検出器に当たる電子の経路長の 99 % = 105 nm、99.9 % = 163 nm。
    /// 画像下端は 1/μ = 4.7)。上限を超えた分は画素ごとに同じエネルギーの中で配り直される (FillPathLengthWeights) ので、
    /// T は十分な余裕を持って取る。</para>
    /// <para>統計に入れるのは、結晶内に源を持つ (d ≥ a) 射出電子のうち、<paramref name="detector"/> を渡したときは<b>検出器画像に当たるもの</b>、
    /// 渡さないときは μ ≥ <see cref="MinMuWithoutDetector"/> のもの。該当が無いときは垂直深さの 99 % 点 / 0.5。</para>
    /// <para>⚠ 表示の視野を検出器より広げると、T が足りない方向 (μ の小さい方向) が出る。そこは近似 (尾部の配り直し) になる。</para></summary>
    public static double ComputePathLengthUpperBound((double Depth, V3 Vec, double Energy)[] bseList, double sampleTilt,
        EbsdDetectorGeometry detector = null, double amorphousLayerNm = 0, double quantile = 0.999)
    {
        if (bseList == null || bseList.Length == 0) return 50.0;
        double a = amorphousLayerNm > 0 && double.IsFinite(amorphousLayerNm) ? amorphousLayerNm : 0;
        var (sinSmp, cosSmp) = Math.SinCos(sampleTilt);
        var paths = new List<double>(bseList.Length / 2);
        var normals = new List<double>(bseList.Length);
        foreach (var (depth, vec, _) in bseList)
        {
            if (!(depth >= a)) continue; //非晶質層内の源は変調なし成分なので格子に関係しない (NaN もここで落ちる)
            double d = depth - a;
            normals.Add(d);
            //lab → 試料系 (ctor と同じ回転)
            double sy = vec.Y * cosSmp + vec.Z * sinSmp, sz = -vec.Y * sinSmp + vec.Z * cosSmp;
            double len = Math.Sqrt(vec.X * vec.X + sy * sy + sz * sz);
            if (!(sz > 0) || !(len > 0)) continue;
            double mu = sz / len;
            if (detector != null)
            {
                var pix = detector.SampleDirectionToPixel(new V3(vec.X, sy, sz));
                if (pix is not { } p || p.Col < -0.5 || p.Col > detector.WidthPx - 0.5 || p.Row < -0.5 || p.Row > detector.HeightPx - 0.5) continue; //画像の外
            }
            else if (mu < MinMuWithoutDetector) continue;
            paths.Add(d / mu);
        }
        if (paths.Count == 0)
        {
            if (normals.Count == 0) return 50.0;
            var nArr = normals.ToArray();
            int k = Math.Min((int)(nArr.Length * 0.99), nArr.Length - 1);
            QuickSelect.Execute(nArr.AsSpan(), k, static (x, y) => x.CompareTo(y));
            return Math.Max(1.0, nArr[k] / 0.5);
        }
        var arr = paths.ToArray();
        int idx = Math.Clamp((int)(arr.Length * quantile), 0, arr.Length - 1);
        QuickSelect.Execute(arr.AsSpan(), idx, static (x, y) => x.CompareTo(y));
        return Math.Max(1.0, arr[idx]);
    }

    public static (double[] energies, double energyStart, double energyEnd, double energyStep,
                    double[] depths, double depthStart, double depthEnd, double depthStep)
        // ComputeGridFromRanges(double beamEnergy, double energyLoss80, double depth99) // 260919Cl 変更前
        //260921Cl シグネチャ変更 (深さ写像 A2): 深さ格子の形状パラメータ β を追加。第 3 引数は経路長の上限 T として使う (ComputePathLengthUpperBound)
        //旧: ComputeGridFromRanges(double beamEnergy, double energyLoss95, double depth99) // 260919Cl 変更: 第 2 引数は 95 パーセンタイルのエネルギー損失
        ComputeGridFromRanges(double beamEnergy, double energyLoss95, double depth99, double depthGridBeta = DefaultDepthGridBeta)
    {
        // int numEnergyLevels = 8; // 260919Cl 変更前
        int numEnergyLevels = 16; // 260919Cl 変更: 8 → 16 段
        double rawEnergyStep = energyLoss95 / (numEnergyLevels - 1); // 260919Cl energyLoss80 → energyLoss95
        double energyStep = Math.Max(0.1, Math.Round(rawEnergyStep / 0.1) * 0.1);
        double energyStart = beamEnergy;
        double energyEnd = Math.Round((beamEnergy - energyStep * (numEnergyLevels - 1)) / 0.1) * 0.1;
        if (energyEnd < 1) energyEnd = 1;

        var energyList = new List<double>();
        for (double e = energyStart; e >= energyEnd - 0.001; e -= energyStep)
            energyList.Add(Math.Round(e * 10) / 10.0);
        var energies = energyList.ToArray();

        int maxDepthDivisions = 40;
        //260921Cl 追加 (深さ写像 A2): β > 0 なら不等間隔 (浅い側を細かく、深い側を幾何級数的に粗く) の 40 点。β ≤ 0 は旧来の等間隔 (以下のコードそのまま)
        if (depthGridBeta > 0 && double.IsFinite(depthGridBeta))
        {
            var g = BuildGeometricDepthGrid(depth99, maxDepthDivisions, depthGridBeta);
            return (energies, energyStart, energyEnd, energyStep,
                    g, g[0], g[^1], g[0]); //depthStep は最初の区間幅 (UI の目安表示用。格子そのものは配列で持ち回すこと)
        }
        double depthStep = Math.Max(0.01, Math.Round(depth99 / maxDepthDivisions * 100) / 100.0); // (260326Ch)
        double depthStart = depthStep;
        double depthEnd = Math.Max(depthStep, Math.Round(depth99 * 100) / 100.0); // (260326Ch)
        var depthList = new List<double>();
        for (double d = depthStart; d <= depthEnd + 0.0001; d += depthStep)
            depthList.Add(Math.Round(d * 100) / 100.0);
        if (depthList.Count > maxDepthDivisions)
            depthList.RemoveRange(maxDepthDivisions, depthList.Count - maxDepthDivisions);
        var depths = depthList.ToArray();

        return (energies, energyStart, energyEnd, energyStep,
                depths, depthStart, depthEnd, depthStep);
    }

    /// <summary>260921Cl 追加 (深さ写像 A2): 経路長の不等間隔格子 t_i = T·expm1(β i/N)/expm1(β)、i = 1…N (t₀ = 0 は暗黙)。
    /// <para>・β → 0 の極限は等間隔 T·i/N。β が大きいほど浅い側が細かい。区間幅は i とともに単調に増える (ゼロ幅は作らない)。</para>
    /// <para>・T は有効数字 2 桁へ切り上げる (MC の乱数で分位点が揺れても格子が変わりにくくし、マスターパターンのキャッシュを効かせるため。
    ///   切り上げなので経路長の被覆は減らない)。各点は 0.001 nm に丸め、狭義単調増加を保証する。</para>
    /// <para>【なぜ不等間隔か】経路長の分布は浅い側に集中する (Si 20 kV で中央値 16 nm、99.9 % 点 163 nm) が、
    /// 下端の画素 (1/μ ≈ 4.7) では 100 nm を超える経路も効く。40 点のまま上限を 3 倍に広げると等間隔では浅い側が 3 倍粗くなる。
    /// マスターパターンは深さ点数に比例してメモリを食う (grid 512 で 40 点 ≈ 1.3 GB) ので点数は増やさない (Codex 推奨 (b))。</para></summary>
    public static double[] BuildGeometricDepthGrid(double upperBoundNm, int count, double beta)
    {
        if (count < 1) throw new ArgumentOutOfRangeException(nameof(count));
        double T = RoundUpToTwoSignificantDigits(Math.Max(upperBoundNm, 0.01));
        var g = new double[count];
        double denom = beta > 1E-9 ? Math.Exp(beta) - 1 : 0;
        for (int i = 1; i <= count; i++)
        {
            double f = beta > 1E-9 ? (Math.Exp(beta * i / count) - 1) / denom : (double)i / count;
            double v = Math.Round(T * f * 1000) / 1000.0;
            if (i > 1 && !(v > g[i - 2])) v = g[i - 2] + 0.001; //丸めで重なったら 0.001 nm だけずらす (ゼロ幅区間を作らない)
            g[i - 1] = v;
        }
        g[count - 1] = Math.Max(g[count - 1], T); //最後の点はちょうど T
        return g;
    }

    /// <summary>260921Cl 追加: 正の値を有効数字 2 桁へ切り上げる (例 162.6 → 170、58.6 → 59、0.123 → 0.13)。</summary>
    static double RoundUpToTwoSignificantDigits(double v)
    {
        if (!(v > 0) || !double.IsFinite(v)) return v;
        double p = Math.Pow(10, Math.Floor(Math.Log10(v)) - 1);
        return Math.Ceiling(v / p - 1E-9) * p;
    }
}
