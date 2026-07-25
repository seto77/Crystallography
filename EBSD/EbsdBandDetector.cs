#region using
using System;
using System.Buffers; //260725Ch: GaussianBlur の大きな一時配列を例外安全に再利用
using System.Collections.Generic;
using System.Linq;
#endregion

namespace Crystallography;

/// <summary>検出された菊池バンド 1 本。座標はすべて native 画像ピクセル (コーナー原点・ピクセル中心規約)。260724Cl 追加</summary>
public sealed class EbsdBand
{
    /// <summary>正規化同次中心線 A·col + B·row + C = 0 (A²+B²=1)</summary>
    public double LineA, LineB, LineC;

    /// <summary>中心線上のアンカー点 (画像内区間を等分)</summary>
    public (double Col, double Row)[] CenterAnchors = [];

    /// <summary>バンド両縁の点 (各アンカー ± 半幅。3D へ戻して sinθB を求めるための観測点)</summary>
    public (double Col, double Row)[] EdgePoints = [];

    /// <summary>バンド幅 (native px、中央アンカーでの推定。表示・初期値用の派生値)</summary>
    public double WidthPx;

    /// <summary>Radon/butterfly ピークスコア (画像の MAD 単位)</summary>
    public double Score;

    /// <summary>中心線の信頼度 (0-1 目安)。260724Cl: 狭窓アンカー射影の中央値 (線検証)</summary>
    public double CenterQuality;

    /// <summary>線検証の診断値。260724Cl 追加: 射影が有意なアンカー割合 / 平均形状の線形フィット R² (edge-runner 判別) / ローブ対バンドか</summary>
    public double VerifyFrac, EdgeLinearity;
    public bool IsLobePair;

    /// <summary>幅推定の信頼度 (0-1 目安。エッジ不明瞭なら低)</summary>
    public double WidthQuality;

    /// <summary>表示用: work 座標系の (θ[deg], ρ[px])</summary>
    public double ThetaDeg, RhoWorkPx;
}

/// <summary>バンド検出の所要時間内訳 (ベンチマーク用)。260724Cl 追加</summary>
public sealed class EbsdBandDetectionTiming
{
    public double PreprocessMs, RadonMs, ButterflyMs, PeakMs, RefineMs, WidthMs;
    public override string ToString()
        => $"pre {PreprocessMs:f0} / radon {RadonMs:f0} / butterfly {ButterflyMs:f0} / peak {PeakMs:f0} / refine {RefineMs:f0} / width {WidthMs:f0} ms";
}

/// <summary>
/// 実測 EBSD パターンからの菊池バンド検出 (Radon 変換 + マルチスケール butterfly カーネル)。260724Cl 追加。GUI 非依存。
/// 手順: 長辺 256px へ縮小 → ガウシアン背景除算 → 正規化 → Radon (sum/√N 正規化) → 幅バンク butterfly 畳み込み →
/// MAD 閾値 + 画像内線距離 NMS → θ/ρ 局所再探索 → 横断プロファイルからエッジ・幅推定。
/// </summary>
public static class EbsdBandDetector
{
    const int WorkLongSide = 256;
    //260724Cl: 384 への引き上げを等価スケール定数付きで実測したが撤回 — ダウンサンプル平均が浅くなる分ピクセルノイズが増え
    //線積分の総情報量は不変 (SNR √1.5 改善の仮定は誤り)。むしろ ρ ビン細分化で幅広バンドの応答が分散し、強バンド (4-1_33 θ84) が
    //pair/検証相互作用の変動で落ちるなど全画像で劣化 (recall 30→28/49、precision 88.2→77.8%)。再チューニングは過適合リスクのため断念
    //260724Cl: 位置次元の定数は 256px 基準の値を WS で等価スケールする (解像度変更時の再チューニング回避、384 実験の遺産)
    const double WS = WorkLongSide / 256.0;
    //static readonly double[] WidthBank = [3, 5, 7, 10, 14, 20, 28, 38]; //work px。260724Cl: 太バンド (低指数) 対応で 38 を追加
    static readonly double[] WidthBank = [.. new double[] { 3, 5, 7, 10, 14, 20, 28, 38 }.Select(v => v * WS)]; //260724Cl: WS スケール

    const double RadonCoreThetaStepDeg = 0.5; //260724Cl: Detect と ComputeRadonMap で共有する θ 刻み

    /// <summary>前処理〜butterfly 平滑応答までの共有コア成果物。260724Cl 追加 (Detect と ComputeRadonMap で共有)</summary>
    sealed class RadonCore
    {
        public double[] Work; public bool[] WorkValid;
        public int W, H; public double Scale, CxW, CyW;
        public double[] Smoothed, BestWidth;
        public int NTheta, NRho, RhoOffset;
    }

    /// <summary>前処理 (縮小→背景除算→正規化) → Radon → butterfly バンク → θ 平滑までを計算する。260724Cl 追加 (Detect 前段の切り出し)</summary>
    static RadonCore ComputeCore(double[] values, int width, int height, bool[] valid, EbsdBandDetectionTiming timing)
    {
        //260725Ch: 公開 Detect/ComputeRadonMap から不正な寸法が入ったとき、並列ループ内の分かりにくい範囲外例外や NaN 連鎖にしない
        ArgumentNullException.ThrowIfNull(values);
        if (width <= 0) throw new ArgumentOutOfRangeException(nameof(width));
        if (height <= 0) throw new ArgumentOutOfRangeException(nameof(height));
        if (values.Length != checked(width * height))
            throw new ArgumentException("values.Length must equal width * height.", nameof(values));
        if (valid != null && valid.Length != values.Length)
            throw new ArgumentException("valid.Length must equal values.Length.", nameof(valid));

        var sw = System.Diagnostics.Stopwatch.StartNew();

        #region 前処理: 縮小 → 背景除算 → 正規化
        double scale = (double)WorkLongSide / Math.Max(width, height);
        if (scale > 1) scale = 1;
        int w = Math.Max(16, (int)Math.Round(width * scale)), h = Math.Max(16, (int)Math.Round(height * scale));
        var work = Downsample(values, width, height, valid, w, h, out var workValid);

        //広域ガウシアン (σ=0.2×短辺) で背景推定し除算 (蛍光体照明の乗算性)
        //260724Cl: 0.1→0.2×短辺。σ がバンド幅と同オーダーだと太バンド自体が背景として除去され縁だけが残り (ハイパス効果)、
        //中心線検出が縁へ引かれる (Cu2(OH)3Cl 4-1_33 の太バンドで顕在化)
        var bg = GaussianBlur(work, workValid, w, h, 0.20 * Math.Min(w, h));
        double floor = Math.Max(1E-10, work.Where((_, i) => workValid[i]).DefaultIfEmpty(1).Average() * 0.05);
        for (int i = 0; i < work.Length; i++)
            work[i] = workValid[i] ? work[i] / Math.Max(bg[i], floor) : 0;

        //軽い平滑 (σ≈0.8) → zero-mean/unit-variance → ±5σ clip (検出専用)
        work = GaussianBlur(work, workValid, w, h, 0.8 * WS); //260724Cl: WS スケール (native 換算の平滑量を維持)
        Normalize(work, workValid, clipSigma: 5);
        if (timing != null) timing.PreprocessMs = sw.Elapsed.TotalMilliseconds;
        #endregion

        #region Radon 変換 (θ 0.5° 刻み、ρ 1 work-px 刻み、sum/√N 正規化)
        sw.Restart();
        const double thetaStepDeg = RadonCoreThetaStepDeg;
        int nTheta = (int)(180 / thetaStepDeg);
        double cxw = w / 2.0 - 0.5, cyw = h / 2.0 - 0.5;
        double rhoMax = Math.Sqrt(w * w + h * h) / 2;
        int nRho = 2 * (int)Math.Ceiling(rhoMax) + 1;
        int rhoOffset = nRho / 2;
        var radon = new double[nTheta * nRho];

        System.Threading.Tasks.Parallel.For(0, nTheta, it =>
        {
            var (sinT, cosT) = Math.SinCos(it * thetaStepDeg * Math.PI / 180);
            for (int ir = 0; ir < nRho; ir++)
            {
                double rho = ir - rhoOffset;
                radon[it * nRho + ir] = LineIntegral(work, workValid, w, h, cxw, cyw, cosT, sinT, rho);
            }
        });
        if (timing != null) timing.RadonMs = sw.Elapsed.TotalMilliseconds;
        #endregion

        #region butterfly バンク畳み込み (ρ 方向 1D、幅横断 max) + θ 方向 3 タップ平滑
        sw.Restart();
        var response = new double[nTheta * nRho];
        var bestWidth = new double[nTheta * nRho];
        //260725Ch: ButterflyKernel のキャッシュ lock を θ 行×幅バンク回数だけ踏まないよう、並列領域の前で一度だけ解決する
        //var kernels = new double[WidthBank.Length][]; //260725Ch 変更前: ComputeCoreごとに参照配列を再確保
        //for (int wi = 0; wi < WidthBank.Length; wi++) kernels[wi] = ButterflyKernel(WidthBank[wi]);
        var kernels = WidthKernels; //260725Ch: 固定WidthBankのカーネル参照も型初期化時に一度だけ構築
        System.Threading.Tasks.Parallel.For(0, nTheta, it =>
        {
            //var row = new double[nRho]; //260725Ch 変更前: θ 行ごとに配列を確保して radon から全コピー
            //Array.Copy(radon, it * nRho, row, 0, nRho);
            int rowOffset = it * nRho; //260725Ch: 元配列の行を直接参照して 360 回の確保・コピーを除去
            //foreach (var bw in WidthBank) //260725Ch 変更前: 各反復で ButterflyKernel の lock 付きキャッシュ検索
            for (int wi = 0; wi < WidthBank.Length; wi++)
            {
                double bw = WidthBank[wi];
                var kernel = kernels[wi]; //260725Ch
                int half = kernel.Length / 2;
                for (int ir = 0; ir < nRho; ir++)
                {
                    double s = 0;
                    for (int k = -half; k <= half; k++)
                    {
                        int j = ir + k;
                        //if ((uint)j < (uint)nRho) s += row[j] * kernel[k + half]; //260725Ch 変更前
                        if ((uint)j < (uint)nRho) s += radon[rowOffset + j] * kernel[k + half]; //260725Ch
                    }
                    //260724Cl: |s| 最大の応答を「符号付き」で保持。明バンド=正応答、暗バンド (deficiency) =負応答。
                    //符号を捨てる (旧 |s| 化) と、zero-mean カーネルが明バンド両脇に作る負のサイドローブ (振幅~50%) が
                    //独立ピークに昇格し「1 バンドに 2-3 本」の偽線が出る。符号はサイドローブ抑制 NMS (異符号判定) に使う。
                    //if (Math.Abs(s) > Math.Abs(response[it * nRho + ir])) { response[it * nRho + ir] = s; bestWidth[it * nRho + ir] = bw; } //260725Ch 変更前
                    if (Math.Abs(s) > Math.Abs(response[rowOffset + ir])) { response[rowOffset + ir] = s; bestWidth[rowOffset + ir] = bw; } //260725Ch
                }
            }
        });
        //θ 方向 3 タップ (循環: θ と θ+180 で ρ 符号反転)
        //260724Cl: 中心と同符号の隣接のみ加算 (交差・重畳領域で異符号応答が相殺し、正解バンドの応答が沈むのを防ぐ。Codex 指摘)
        //260724Cl: 多段化実験 — 同符号 3 タップを ThetaSmoothPasses 回反復 (2 回で三角 5 タップ相当)。
        //2〜3 パスは定量評価で完全中立 (recall/precision 不変)、かつ多段でも候補プール拡大の precision 崩壊 (88→64%) は救えなかった
        //= 微弱バンドの Radon 証拠は θ 統合で浮上する水準になく、これ以上の回収は結晶学的事前知識 (指数付けフィードバック探索) が必要。
        //既定は従来同等の 1 パスに固定
        const int ThetaSmoothPasses = 1;
        var smoothed = response;
        for (int pass = 0; pass < ThetaSmoothPasses; pass++)
        {
            var src = smoothed;
            var dst = new double[nTheta * nRho];
            System.Threading.Tasks.Parallel.For(0, nTheta, it =>
            {
                for (int ir = 0; ir < nRho; ir++)
                {
                    double s0 = src[it * nRho + ir];
                    double s1 = NeighborTheta(src, nTheta, nRho, it - 1, ir);
                    double s2 = NeighborTheta(src, nTheta, nRho, it + 1, ir);
                    double s = 2 * s0; int n = 2;
                    if (s1 * s0 > 0) { s += s1; n++; }
                    if (s2 * s0 > 0) { s += s2; n++; }
                    dst[it * nRho + ir] = s / n;
                }
            });
            smoothed = dst;
        }
        if (timing != null) timing.ButterflyMs = sw.Elapsed.TotalMilliseconds;
        #endregion

        return new RadonCore
        {
            Work = work, WorkValid = workValid, W = w, H = h, Scale = scale, CxW = cxw, CyW = cyw,
            Smoothed = smoothed, BestWidth = bestWidth, NTheta = nTheta, NRho = nRho, RhoOffset = rhoOffset,
        };
    }

    /// <summary>
    /// 実測パターンから Radon/butterfly 応答マップ (|応答| + robust null 統計) を計算する。260724Cl 追加。
    /// バンドを離散検出せず、EbsdRadonIndexer の方位テンプレート照合の証拠マップとして使う。
    /// </summary>
    public static EbsdRadonMap ComputeRadonMap(double[] values, int width, int height, bool[] valid = null, EbsdBandDetectionTiming timing = null)
    {
        var core = ComputeCore(values, width, height, valid, timing);
        var abs = new double[core.Smoothed.Length];
        for (int i = 0; i < abs.Length; i++) abs[i] = Math.Abs(core.Smoothed[i]);
        //var sorted = (double[])abs.Clone(); Array.Sort(sorted); double median = sorted[sorted.Length / 2]; //260725Ch 変更前
        //double mad = sorted.Select(v => Math.Abs(v - median)).OrderBy(v => v).ElementAt(sorted.Length / 2);
        var (median, mad) = MedianAndMad(abs); //260725Ch: 同じ scratch を median と MAD の2回のソートに再利用
        return new EbsdRadonMap
        {
            Abs = abs, NTheta = core.NTheta, NRho = core.NRho, RhoOffset = core.RhoOffset,
            ThetaStepDeg = RadonCoreThetaStepDeg, WorkW = core.W, WorkH = core.H, Scale = core.Scale,
            Mu0 = median, Sigma0 = Math.Max(1.4826 * mad, 1E-12),
        };
    }

    /// <summary>バンドを検出する。values は native 画像の生強度 (row-major)、valid は有効画素マスク (null=全有効)、debugLog は診断ログ出力 (検証ハーネス用、null=無効)</summary>
    //260724Cl (/simplify) シグネチャ変更: 診断フックを internal static フィールド DebugLog (reflection 設定) から引数へ (グローバル可変状態の排除)
    //旧: public static List<EbsdBand> Detect(double[] values, int width, int height, bool[] valid = null, int maxBands = 12, EbsdBandDetectionTiming timing = null)
    public static List<EbsdBand> Detect(double[] values, int width, int height, bool[] valid = null, int maxBands = 12, EbsdBandDetectionTiming timing = null, Action<string> debugLog = null)
    {
        //260724Cl: 前処理〜butterfly 平滑応答までを ComputeCore へ切り出し (Radon 方位テンプレート照合 ComputeRadonMap と共有)
        var core = ComputeCore(values, width, height, valid, timing);
        double scale = core.Scale, cxw = core.CxW, cyw = core.CyW;
        int w = core.W, h = core.H, nTheta = core.NTheta, nRho = core.NRho, rhoOffset = core.RhoOffset;
        var work = core.Work; var workValid = core.WorkValid;
        var smoothed = core.Smoothed; var bestWidth = core.BestWidth;
        const double thetaStepDeg = RadonCoreThetaStepDeg;
        var sw = System.Diagnostics.Stopwatch.StartNew();

        #region ピーク検出 (|応答| の median+MAD 閾値、サイドローブ抑制 NMS)
        sw.Restart();
        //260724Cl: 閾値・極大判定は |smoothed| で行い (明暗両極性)、符号はピーク属性として保持する
        var absSmoothed = new double[smoothed.Length];
        for (int i = 0; i < smoothed.Length; i++) absSmoothed[i] = Math.Abs(smoothed[i]);
        //var sorted = (double[])absSmoothed.Clone(); Array.Sort(sorted); double median = sorted[sorted.Length / 2]; //260725Ch 変更前
        //double mad = sorted.Select(v => Math.Abs(v - median)).OrderBy(v => v).ElementAt(sorted.Length / 2);
        var (median, mad) = MedianAndMad(absSmoothed); //260725Ch
        if (mad < 1E-12) mad = 1E-12;
        double threshold = median + 4.5 * mad;

        var peakIdx = new List<int>();
        for (int it = 0; it < nTheta; it++)
            for (int ir = 1; ir < nRho - 1; ir++)
            {
                int i = it * nRho + ir;
                if (absSmoothed[i] > threshold && absSmoothed[i] >= absSmoothed[i - 1] && absSmoothed[i] >= absSmoothed[i + 1]
                    && absSmoothed[i] >= Math.Abs(NeighborTheta(smoothed, nTheta, nRho, it - 1, ir)) && absSmoothed[i] >= Math.Abs(NeighborTheta(smoothed, nTheta, nRho, it + 1, ir)))
                    peakIdx.Add(i);
            }

        var candidates = peakIdx.OrderByDescending(i => absSmoothed[i]).ToList();
        if (debugLog != null) //260724Cl: 診断 — NMS 前の候補上位 40 (work 座標)
            foreach (var (ci, rank) in candidates.Take(40).Select((c, r) => (c, r)))
                debugLog($"cand#{rank}: θ={(ci / nRho) * thetaStepDeg,6:f1} ρ={ci % nRho - rhoOffset,5} score={(absSmoothed[ci] - median) / mad,5:f1} sign={Math.Sign(smoothed[ci])} w={bestWidth[ci]:f0}");
        //260724Cl: pairDelta = 異符号ローブ対 (excess/deficiency バンド) の相方への符号付き距離 (0=対なし)。
        //サイドローブ (振幅~50%) と違い、本物の明/暗ローブ対はスコア比が高い (実測 0.85-0.98) — 比 0.7 以上を対と判定し、
        //棄却はそのまま (二重線は出さない) だが survivor にタグ付けして BuildBand の反対称センタリングに使う
        var picked = new List<(double theta, double rho, double score, double bw, int sign, double pairDelta, double pairScore)>();
        foreach (var i in candidates)
        {
            double th = (i / nRho) * thetaStepDeg * Math.PI / 180, rho = i % nRho - rhoOffset;
            double bwCand = bestWidth[i];
            int sgn = Math.Sign(smoothed[i]);
            double scoreCand = (absSmoothed[i] - median) / mad;
            //260724Cl: NMS を 2 本立てに —
            //  同符号: 6px (実在する近接平行バンドは残す)
            //  異符号: max(6, 0.75×max(幅)) (zero-mean カーネルのサイドローブは主ピークの ±(w/2+σ) に必ず「異符号」で立つため、
            //          先に採用された強いピークの幅レンジ内にある異符号ピークはサイドローブとして棄却する)
            bool suppressed = false;
            //foreach (var p in picked)
            for (int pi = 0; pi < picked.Count; pi++) //260724Cl: 対タグ更新のため index ループへ
            {
                var p = picked[pi];
                double dist = LineRmsDistance(p.theta, p.rho, th, rho, w, h, cxw, cyw);
                //260724Cl: 異符号の幅依存排他は「ほぼ平行 (θ差<10°)」ペアに限定 (Codex 指摘: サイドローブは主ピークとほぼ同 θ に立つ。
                //交差・斜交する実在の明/暗バンドまで広半径で消さない)
                bool nearParallel = ThetaDiffDeg(p.theta, th) < 10;
                //double limit = p.sign == sgn || !nearParallel ? 6 : Math.Max(6, 0.75 * Math.Max(p.bw, bwCand));
                double limit = p.sign == sgn || !nearParallel ? 6 * WS : Math.Max(6 * WS, 0.75 * Math.Max(p.bw, bwCand)); //260724Cl: WS スケール
                if (dist < limit)
                {
                    suppressed = true;
                    //260724Cl: 異符号・近平行かつスコア比 ≥0.7 なら本物のローブ対 → survivor に相方位置を記録 (最初=最強の相方のみ)
                    if (p.sign != sgn && nearParallel && p.pairDelta == 0 && scoreCand >= 0.7 * p.score)
                    {
                        //相方線の中心最近点から survivor 線への符号付き距離 = ρc·cos(θc−θp) − ρp (θ 180° 折返しも吸収)
                        double delta = rho * Math.Cos(th - p.theta) - p.rho;
                        picked[pi] = (p.theta, p.rho, p.score, p.bw, p.sign, delta, scoreCand);
                    }
                    break;
                }
            }
            if (suppressed) continue;
            picked.Add((th, rho, scoreCand, bwCand, sgn, 0, 0));
            if (picked.Count >= maxBands) break;
            //260724Cl: 候補プール 2 倍化 (maxBands*2) を試したが撤回 — スコア 7-10 帯の追加候補は微弱実バンドでなくほぼ偽線で、
            //低振幅域では線検証 (med≥0.5) も擦り抜ける (precision 85.7→64.4% に悪化)。微弱バンドは Radon 段階の証拠が弱く、
            //拾いに行くと「誤配置は見逃しより重罪」の方針に反する
        }
        if (timing != null) timing.PeakMs = sw.Elapsed.TotalMilliseconds;
        #endregion

        #region 局所再探索 (θ ±1°を0.1°、ρ ±2px を 0.2px) — butterfly 応答をピーク符号の方向に最大化
        sw.Restart();
        //260724Cl (/simplify): バンド毎に独立なので並列化。さらに ρ 側は θ 毎に 0.2px 刻みの 1D プロファイルを 1 回だけ算出して
        //カーネルと畳み込む (旧: drho×カーネルタップの二重ループで、0.2 グリッド上の同一 ρ の LineIntegral を約 4 回重複評価していた)
        var refinedArr = new (double theta, double rho, double score, double bw, int sign, double pairDelta, double pairScore)[picked.Count]; //260724Cl: 対タグを持ち回り
        System.Threading.Tasks.Parallel.For(0, picked.Count, pi =>
        {
            var p = picked[pi];
            var kernel = ButterflyKernel(p.bw);
            int half = kernel.Length / 2;
            int nProf = 10 * half + 21; //ρ オフセット -2-half .. 2+half (0.2 刻み)。drho + k = -2 - half + 0.2m と一対一
            var prof = new double[nProf];
            double bestTh = p.theta, bestRho = p.rho, bestVal = double.MinValue;
            for (double dth = -1; dth <= 1.0001; dth += 0.1)
            {
                double th = p.theta + dth * Math.PI / 180;
                var (sinT, cosT) = Math.SinCos(th);
                for (int m = 0; m < nProf; m++)
                    prof[m] = LineIntegral(work, workValid, w, h, cxw, cyw, cosT, sinT, p.rho - 2 - half + 0.2 * m);
                for (int j = 0; j <= 20; j++) //drho = -2 + 0.2j
                {
                    double s = 0;
                    for (int k = 0; k < kernel.Length; k++)
                        s += kernel[k] * prof[j + 5 * k];
                    s *= p.sign; //260724Cl: 暗バンド (負応答) はより負へ = |応答| 最大化
                    if (s > bestVal) { bestVal = s; bestTh = th; bestRho = p.rho - 2 + 0.2 * j; }
                }
            }
            //refinedArr[pi] = (bestTh, bestRho, p.score, p.bw, p.sign);
            refinedArr[pi] = (bestTh, bestRho, p.score, p.bw, p.sign, p.pairDelta, p.pairScore); //260724Cl
        });
        var refined = refinedArr.ToList();

        //260724Cl 追加: 局所再探索で同一線へ収束した重複を除去 (スコア降順で高スコア側を残す)。
        //異符号ペアはサイドローブ抑制と同じ幅依存距離で除去 (再探索でピーク段 NMS を通過した残党が近づくことがある)
        refined = [.. refined.OrderByDescending(r => r.score)];
        for (int i = 0; i < refined.Count; i++)
            for (int j = refined.Count - 1; j > i; j--)
            {
                double dist = LineRmsDistance(refined[i].theta, refined[i].rho, refined[j].theta, refined[j].rho, w, h, cxw, cyw);
                bool nearParallel = ThetaDiffDeg(refined[i].theta, refined[j].theta) < 10; //260724Cl: ピーク段 NMS と同じ平行限定
                //double limit = refined[i].sign == refined[j].sign || !nearParallel ? 5 : Math.Max(5, 0.75 * Math.Max(refined[i].bw, refined[j].bw));
                double limit = refined[i].sign == refined[j].sign || !nearParallel ? 5 * WS : Math.Max(5 * WS, 0.75 * Math.Max(refined[i].bw, refined[j].bw)); //260724Cl: WS スケール
                if (dist < limit) refined.RemoveAt(j);
            }

        //260724Cl: 旧「縁ペア統合」は撤去。サイドローブ抑制 NMS (符号ベース) が原因 (1 バンド 2 本) を上流で断つため不要になり、
        //実在する明バンド・暗バンドの隣接ペアを誤って 1 本に潰すリスクの方が大きい。
        if (timing != null) timing.RefineMs = sw.Elapsed.TotalMilliseconds;
        #endregion

        #region エッジ・幅推定 (横断平均プロファイル → 1次微分極値) + native 座標へ変換
        sw.Restart();
        var bands = new List<EbsdBand>();
        foreach (var p in refined)
        {
            //var band = BuildBand(p.theta, p.rho, p.score, p.bw, work, workValid, w, h, cxw, cyw, scale, width, height);
            var band = BuildBand(p.theta, p.rho, p.score, p.bw, p.sign, p.pairDelta, work, workValid, w, h, cxw, cyw, scale, width, height, debugLog); //260724Cl: 対タグ+符号+診断ログを追加
            if (band != null) bands.Add(band);
        }

        //260724Cl 追加: 線検証ゲート — 「間違った場所に中心線を引く」ことを「見逃し」より重い罪として濾過する (作者方針)。
        //CenterQuality (狭窓アンカー射影の中央値) が低い = 線に沿って共有される横断構造がない (交差点連結・幽霊)。
        //Radon+butterfly スコアが非常に強い (≥17 MAD) 場合は大域証拠を尊重してバーを 0.3 へ緩和
        //(交差密集地帯を貫く幅広の主要バンドはアンカー汚染で med が沈む — 4-1_33 中央バンドの救済。定量評価で決定)。
        //EdgeLinearity ≥ 0.55 = 平均形状が単調勾配 = バンド縁を走る線 (ローブ対バンドは S 字形状が本質なので免除)
        bands.RemoveAll(b => b.CenterQuality < (b.Score >= 17 ? 0.3 : 0.5)
            || (!b.IsLobePair && b.Score < 17 && b.EdgeLinearity >= 0.55) //260724Cl: r2 kill は弱スコア線に限定 (強証拠の実バンドまで殺さない)
            || b.WidthQuality <= 0.05); //260724Cl: 自窓内でエッジが見つからない線は誤配置の疑い (正解一致線の widthQ 最小実測 0.18、偽線に 0.00 が頻出)
        if (timing != null) timing.WidthMs = sw.Elapsed.TotalMilliseconds;
        #endregion

        return bands;
    }

    //260724Cl (/simplify): Let 拡張メソッド (timing?.Let(t => ...) 用) は間接化の価値がなく削除。直接 if (timing != null) 代入へ

    /// <summary>2 直線の法線角差 (deg、180° 循環考慮)。260724Cl 追加</summary>
    static double ThetaDiffDeg(double th1, double th2)
    {
        double d = Math.Abs(th1 - th2) * 180 / Math.PI % 180;
        return Math.Min(d, 180 - d);
    }

    static double NeighborTheta(double[] map, int nTheta, int nRho, int it, int ir)
    {
        //θ の循環: θ+180° は ρ 符号反転
        if (it < 0) { it += nTheta; ir = nRho - 1 - ir; }
        else if (it >= nTheta) { it -= nTheta; ir = nRho - 1 - ir; }
        return (uint)ir < (uint)nRho ? map[it * nRho + ir] : 0;
    }

    /// <summary>値の median と MAD を、1 本の scratch 配列を再利用して求める。260725Ch 追加</summary>
    static (double Median, double Mad) MedianAndMad(double[] values)
    {
        var scratch = (double[])values.Clone();
        Array.Sort(scratch);
        double median = scratch[scratch.Length / 2];
        for (int i = 0; i < scratch.Length; i++) scratch[i] = Math.Abs(values[i] - median);
        Array.Sort(scratch);
        return (median, scratch[scratch.Length / 2]);
    }

    /// <summary>直線 (法線角 θ、距離 ρ、画像中心基準) に沿った線積分 sum/√N (bilinear、1px ステップ)</summary>
    static double LineIntegral(double[] img, bool[] valid, int w, int h, double cx, double cy, double cosT, double sinT, double rho)
    {
        //線: (x-cx)·cosT + (y-cy)·sinT = ρ。方向 = (-sinT, cosT)
        double px = cx + rho * cosT, py = cy + rho * sinT;
        double dx = -sinT, dy = cosT;
        //画像矩形との交差区間 [t0,t1]
        if (!ClipLine(px, py, dx, dy, w, h, out double t0, out double t1)) return 0;
        double sum = 0; int n = 0;
        for (double t = t0; t <= t1; t += 1.0)
        {
            double x = px + t * dx, y = py + t * dy;
            int x0 = (int)x, y0 = (int)y;
            if ((uint)x0 >= (uint)(w - 1) || (uint)y0 >= (uint)(h - 1)) continue;
            int i = y0 * w + x0;
            if (!valid[i]) continue;
            double fx = x - x0, fy = y - y0;
            sum += (img[i] * (1 - fx) + img[i + 1] * fx) * (1 - fy) + (img[i + w] * (1 - fx) + img[i + w + 1] * fx) * fy;
            n++;
        }
        //return n >= 8 ? sum / Math.Sqrt(n) : 0;
        return n >= 8 * WS ? sum / Math.Sqrt(n) : 0; //260724Cl: 最小サンプル数を WS スケール
    }

    static bool ClipLine(double px, double py, double dx, double dy, int w, int h, out double t0, out double t1)
    {
        t0 = double.MinValue; t1 = double.MaxValue;
        //Liang-Barsky (境界 0..w-1, 0..h-1)
        /*260725Ch 変更前: LineIntegral ごと (通常 13 万回以上) に 4 要素タプル配列を確保していた
        foreach (var (p, q) in new[] { (-dx, px - 0.0), (dx, w - 1.0 - px), (-dy, py - 0.0), (dy, h - 1.0 - py) })
        {
            if (Math.Abs(p) < 1E-12) { if (q < 0) return false; }
            else
            {
                double r = q / p;
                if (p < 0) { if (r > t1) return false; if (r > t0) t0 = r; }
                else { if (r < t0) return false; if (r < t1) t1 = r; }
            }
        }
        */
        //260725Ch: 4 境界を直接評価し、ホットループのヒープ確保をゼロにする
        if (!ClipBoundary(-dx, px, ref t0, ref t1)
            || !ClipBoundary(dx, w - 1.0 - px, ref t0, ref t1)
            || !ClipBoundary(-dy, py, ref t0, ref t1)
            || !ClipBoundary(dy, h - 1.0 - py, ref t0, ref t1))
            return false;
        return t1 > t0;

        static bool ClipBoundary(double p, double q, ref double lower, ref double upper)
        {
            if (Math.Abs(p) < 1E-12) return q >= 0;
            double r = q / p;
            if (p < 0)
            {
                if (r > upper) return false;
                if (r > lower) lower = r;
            }
            else
            {
                if (r < lower) return false;
                if (r < upper) upper = r;
            }
            return true;
        }
    }

    /// <summary>butterfly 1D カーネル K_w(t)=G_σ(t) − ½G_σ(t−w/2) − ½G_σ(t+w/2)、σ=max(0.75, 0.12w)。平均0・L2=1 に正規化</summary>
    static readonly Dictionary<double, double[]> kernelCache = [];
    static readonly double[][] WidthKernels = [.. WidthBank.Select(ButterflyKernel)]; //260725Ch: kernelCache初期化後に固定バンクを一括解決
    static double[] ButterflyKernel(double w)
    {
        lock (kernelCache)
        {
            if (kernelCache.TryGetValue(w, out var cached)) return cached;
            double sigma = Math.Max(0.75, 0.12 * w);
            int half = (int)Math.Ceiling(w / 2 + 3 * sigma);
            var k = new double[2 * half + 1];
            for (int i = -half; i <= half; i++)
                k[i + half] = Gauss(i, sigma) - 0.5 * Gauss(i - w / 2, sigma) - 0.5 * Gauss(i + w / 2, sigma);
            double mean = k.Average();
            for (int i = 0; i < k.Length; i++) k[i] -= mean;
            double norm = Math.Sqrt(k.Sum(v => v * v));
            for (int i = 0; i < k.Length; i++) k[i] /= norm;
            kernelCache[w] = k;
            return k;
        }
    }

    static double Gauss(double x, double sigma) => Math.Exp(-x * x / (2 * sigma * sigma));

    /// <summary>2 直線の画像内 RMS 距離 (同一バンド判定用の簡易指標: 画像内 5 点で相互距離)</summary>
    static double LineRmsDistance(double th1, double rho1, double th2, double rho2, int w, int h, double cx, double cy)
    {
        var (sin1, cos1) = Math.SinCos(th1);
        double px = cx + rho1 * cos1, py = cy + rho1 * sin1, dx = -sin1, dy = cos1;
        if (!ClipLine(px, py, dx, dy, w, h, out double t0, out double t1)) return double.MaxValue;
        var (sin2, cos2) = Math.SinCos(th2);
        double sum = 0;
        for (int i = 0; i < 5; i++)
        {
            double t = t0 + (t1 - t0) * i / 4.0;
            double x = px + t * dx, y = py + t * dy;
            double d = (x - cx) * cos2 + (y - cy) * sin2 - rho2;
            sum += d * d;
        }
        return Math.Sqrt(sum / 5);
    }

    /// <summary>精緻化済み (θ,ρ) からエッジ・幅を推定し、native ピクセル座標の EbsdBand を構築する</summary>
    //260724Cl シグネチャ変更: sign (ローブ極性)・pairDelta (異符号ローブ対の相方への符号付き距離、0=対なし)・debugLog を追加
    //旧: static EbsdBand BuildBand(double theta, double rho, double score, double initialWidth,
    //    double[] work, bool[] valid, int w, int h, double cx, double cy, double scale, int nativeW, int nativeH)
    static EbsdBand BuildBand(double theta, double rho, double score, double initialWidth, int sign, double pairDelta,
        double[] work, bool[] valid, int w, int h, double cx, double cy, double scale, int nativeW, int nativeH, Action<string> debugLog = null)
    {
        var (sinT, cosT) = Math.SinCos(theta);
        double px = cx + rho * cosT, py = cy + rho * sinT, dx = -sinT, dy = cosT;
        if (!ClipLine(px, py, dx, dy, w, h, out double t0, out double t1)) return null;
        if (t1 - t0 < Math.Min(w, h) * 0.25) return null; //画像内の支持長が短すぎる線は棄却

        #region 横断平均プロファイル (アンカー 5 点 × ±1.6×初期幅) → 1 次微分極値で両縁
        //int halfRange = (int)Math.Ceiling(Math.Max(initialWidth * 1.3, 8)); //260724Cl: 1.3×幅 (1.6 は隣接バンド混入で悪化することを定量評価で確認)
        //260724Cl: ローブ対バンドは相方まで窓に収める (|Δ|+6)。定数は WS スケール
        int halfRange = (int)Math.Ceiling(Math.Max(Math.Max(initialWidth * 1.3, 8 * WS), pairDelta == 0 ? 0 : Math.Abs(pairDelta) + 6 * WS));
        //260724Cl: 中心決定用プロファイルは従来どおり 5 アンカー平均 (7 に増やすと平均形状が変わり
        //対称相関ゲートの判定が揺れて recall が劣化することを定量評価で確認 — 検証は後段の VerifyLine が別サンプリングで行う)
        var profile = new double[2 * halfRange + 1];
        var count = new int[2 * halfRange + 1];
        const int nAnchor = 5;
        var anchorsT = Enumerable.Range(0, nAnchor).Select(i => t0 + (t1 - t0) * (i + 0.5) / nAnchor).ToArray();
        foreach (var t in anchorsT)
        {
            double ax = px + t * dx, ay = py + t * dy;
            for (int s = -halfRange; s <= halfRange; s++)
            {
                double x = ax + s * cosT, y = ay + s * sinT;
                int x0 = (int)x, y0 = (int)y;
                if ((uint)x0 >= (uint)(w - 1) || (uint)y0 >= (uint)(h - 1)) continue;
                int idx = y0 * w + x0;
                if (!valid[idx]) continue;
                double fx = x - x0, fy = y - y0;
                profile[s + halfRange] += (work[idx] * (1 - fx) + work[idx + 1] * fx) * (1 - fy) + (work[idx + w] * (1 - fx) + work[idx + w + 1] * fx) * fy;
                count[s + halfRange]++;
            }
        }
        for (int i = 0; i < profile.Length; i++)
            profile[i] = count[i] > 0 ? profile[i] / count[i] : 0;

        //3 タップ平滑 → 1 次微分
        var smooth = new double[profile.Length];
        for (int i = 1; i < profile.Length - 1; i++) smooth[i] = (profile[i - 1] + profile[i] + profile[i + 1]) / 3;
        smooth[0] = profile[0]; smooth[^1] = profile[^1];
        var deriv = new double[profile.Length];
        for (int i = 1; i < profile.Length - 1; i++) deriv[i] = (smooth[i + 1] - smooth[i - 1]) / 2;

        //260724Cl 変更 (作者指摘: 急峻な線で赤線が excess 側へ系統的にずれる):
        //中心補正を「エッジ中点」から「プロファイルの対称相関」へ。菊池バンドは excess/deficiency で
        //PC 基準の同じ側が常に明るいため、明るさ重心・片極性エッジ検出は全バンドで同方向へバイアスする。
        //対称相関 S(δ)=Σ_t p̃(c+δ+t)·p̃(c+δ−t) は反対称成分 (excess/deficiency) に不感で、対称成分 (バンド台形) の中心にロックする。
        double meanP = smooth.Average();
        var sym = new double[smooth.Length];
        for (int i = 0; i < sym.Length; i++) sym[i] = smooth[i] - meanP;
        double Sample(double x) //bilinear
        {
            int i0 = (int)Math.Floor(x);
            if (i0 < 0 || i0 >= sym.Length - 1) return 0;
            return sym[i0] * (1 - (x - i0)) + sym[i0 + 1] * (x - i0);
        }
        double SymCorr(double delta) //260724Cl: 対称相関をローブ対センタリングと共用するため関数化 (数式は従来と同一)
        {
            int tMax = (int)(halfRange - Math.Abs(delta) - 1);
            if (tMax < 3) return double.NaN;
            double s = 0;
            for (int t = 1; t <= tMax; t++)
                s += Sample(halfRange + delta + t) * Sample(halfRange + delta - t);
            return s / tMax; //項数で正規化 (δ 依存の項数差を補正)
        }
        double maxShift = Math.Min(halfRange * 0.5, Math.Max(initialWidth * 0.5, 3 * WS)); //260724Cl: 幅比例に厳格化 (細バンドが窓内の隣接構造へ大きく誤移動するのを防ぐ。定量評価で決定)。下限は WS スケール
        double bestShift = 0, bestCorr = double.MinValue, corrAtZero = 0;
        for (double delta = -maxShift; delta <= maxShift + 1E-9; delta += 0.25)
        {
            double s = SymCorr(delta);
            if (double.IsNaN(s)) continue;
            if (Math.Abs(delta) < 0.125) corrAtZero = s;
            if (s > bestCorr) { bestCorr = s; bestShift = delta; }
        }

        //260724Cl: 対称相関のゲート判定を先に行う (pair 経路は「sym が何も見つけられなかった」場合のみのフォールバック)
        {
            //260724Cl: 探索境界に張り付いた解は「真の対称中心が範囲外」のサインなので採用しない (初期線がバンド縁にある場合の誤ロック防止)
            if (Math.Abs(bestShift) > maxShift * 0.9) bestShift = 0;
            //260724Cl: シフトは対称相関が「有意に」改善する場合のみ採用 (S(δ*)>0 かつ S(0) 比 +15% 以上)。
            //微差で動かすと、Radon ピーク (正解位置) にあった線が隣接構造の対称中心へ引っ張られる (定量評価: ref2/ref6 が 7 work px 誤移動して正解から外れていた)
            if (!(bestCorr > 0 && bestCorr > corrAtZero + 0.15 * Math.Abs(corrAtZero))) bestShift = 0;
            //260724Cl: シフト先が「同じバンドの内部」に留まることも要求 (プロファイル値が同符号かつ |強度| 50% 以上維持)。
            //細いバンドでは対称相関だけだと窓内の隣接構造の対称中心が勝つことがある (定量評価: ref7 の 5 work px 誤移動)
            if (bestShift != 0)
            {
                double v0 = Sample(halfRange), v1 = Sample(halfRange + bestShift);
                if (v0 * v1 <= 0 || Math.Abs(v1) < 0.5 * Math.Abs(v0)) bestShift = 0;
            }
        }

        //260724Cl 追加: ローブ対バンド (pairDelta≠0 = NMS で異符号・近平行・スコア比 0.7 以上の相方を棄却した survivor)。
        //excess/deficiency 型の非対称プロファイル (明ピーク+暗谷が隣接) では、正解中心 ≈ 両ローブの中点 = 反対称中心。
        //そこでは S(δ) が強い負になる (p̃(c+t)≈−p̃(c−t)) ので、期待中点 Δ/2 の近傍で S を最小化して中心を求める。
        //ただし「幅広の明バンド+片側暗フリンジ」(正解=バンド本体中心、sym 相関が正しく見つける) と紛らわしいため、
        //sym 相関がシフトを見つけた場合はそちらを優先し、pair 経路は sym 不発時のフォールバックに限定する (3 画像の定量評価で決定)。
        bool pairApplied = false;
        int lobeA = 0, lobeB = 0; double lobeAAmp = 0, lobeBAmp = 0;
        if (pairDelta != 0 && bestShift == 0)
        {
            double expected = pairDelta / 2, span = Math.Max(2 * WS, Math.Abs(pairDelta) * 0.35); //260724Cl: 下限 WS スケール
            double pairShift = 0, pairMin = double.MaxValue;
            for (double delta = expected - span; delta <= expected + span + 1E-9; delta += 0.25)
            {
                double s = SymCorr(delta);
                if (!double.IsNaN(s) && s < pairMin) { pairMin = s; pairShift = delta; }
            }
            //260724Cl: エネルギー正規化反相関 N = S(δ*)/E(δ*) ∈ [-1,1]。真の e/d 対はローブ振幅が拮抗し N→−1、
            //対称バンド+サイドローブ (振幅比~3:1) は N≈−0.6 に留まる (N=−2ab/(a²+b²))
            double energy = 0; int tMaxP = (int)(halfRange - Math.Abs(pairShift) - 1);
            for (int t = 1; t <= tMaxP; t++)
            {
                double p1 = Sample(halfRange + pairShift + t), p2 = Sample(halfRange + pairShift - t);
                energy += 0.5 * (p1 * p1 + p2 * p2);
            }
            energy /= Math.Max(1, tMaxP);
            double normAnti = energy > 1E-12 ? pairMin / energy : 0;
            //両ローブ位置と振幅 (=Kikuchi 線対) をプロファイル極値から実測 (振幅比もゲート候補のため常時計測)
            int c0 = halfRange, cB = (int)Math.Round(halfRange + pairDelta);
            int r = Math.Max((int)Math.Round(2 * WS), (int)(Math.Abs(pairDelta) * 0.4)); //260724Cl: 下限 WS スケール
            (int idx, double amp) FindLobe(int center, int lobeSign)
            {
                int best = center; double bestV = double.MinValue;
                for (int i = Math.Max(1, center - r); i <= Math.Min(sym.Length - 2, center + r); i++)
                    if (lobeSign * sym[i] > bestV) { bestV = lobeSign * sym[i]; best = i; }
                return (best, bestV);
            }
            (lobeA, lobeAAmp) = FindLobe(c0, sign);
            (lobeB, lobeBAmp) = FindLobe(cB, -sign);
            double ampRatio = lobeAAmp > 0 && lobeBAmp > 0 ? Math.Min(lobeAAmp, lobeBAmp) / Math.Max(lobeAAmp, lobeBAmp) : 0;
            //260724Cl: ローブ半値幅 (極値から両側へ |sym|≥½振幅 が続く範囲)。真の e/d 対は両ローブとも細い「線」(≲0.45Δ)、
            //「幅広の明バンド+片側暗フリンジ」はバンド本体が幅広 (実測 0.64-0.8Δ) — この場合の作者正解はバンド本体中心なので pair 化しない
            double LobeHalfWidth(int idx, int lobeSign, double amp)
            {
                double halfAmp = 0.5 * amp;
                int lo = idx, hi = idx;
                while (lo > 0 && lobeSign * sym[lo - 1] >= halfAmp) lo--;
                while (hi < sym.Length - 1 && lobeSign * sym[hi + 1] >= halfAmp) hi++;
                return hi - lo + 1;
            }
            double maxLobeW = Math.Max(LobeHalfWidth(lobeA, sign, lobeAAmp), LobeHalfWidth(lobeB, -sign, lobeBAmp));
            //260724Cl ゲート (3 画像の定量評価+±30-50% プラトー確認で決定):
            //  ① normAnti < −0.55: 真の反対称構造の存在 (−0.4〜−0.7 でスコア不変のプラトー)
            //  ② ampRatio > 0.5: ローブ振幅の拮抗 (サイドローブは~50% 以下)
            //  ③ maxLobeW < 1.0·|Δ|: ローブが「線」であること (実測: 採用すべき対 ≤0.8Δ / 棄却すべき幅広本体 ≥1.2Δ の中点)
            //  ④ |δ*−Δ/2| ≤ 0.5·span: S 最小が期待中点近傍にあること (境界張り付き=別構造への誤ロック防止。0.8 で隣接バンド誤シフトの崖)
            bool ok = pairMin < 0 && normAnti < -0.55 && ampRatio > 0.5
                && maxLobeW < 1.0 * Math.Abs(pairDelta) && Math.Abs(pairShift - expected) <= 0.5 * span;
            debugLog?.Invoke($"pair θ={theta * 180 / Math.PI,6:f1} ρ={rho,6:f1} sign={sign,2} Δ={pairDelta,5:f1}: symMax={bestCorr,6:f2} pairMin={pairMin,6:f2}@{pairShift,5:f2} N={normAnti,6:f2} r={ampRatio:f2} lw={maxLobeW:f0}/{1.0 * Math.Abs(pairDelta):f1} dev={Math.Abs(pairShift - expected):f2}/{0.5 * span:f2} -> {(ok ? "PAIR" : "reject")}");
            if (ok)
            {
                pairApplied = true;
                bestShift = pairShift;
            }
        }
        px += bestShift * cosT; py += bestShift * sinT; //中心線を対称 (または反対称) 中心へ平行移動

        //エッジ: 補正後中心から外側へ走査し、最初の有意な |勾配| 極大を両縁とする。260724Cl 変更
        //(明暗極性に依存せず、かつ範囲内の隣接バンドの強エッジへ飛びついて幅が過大になるのを防ぐ)
        int center = Math.Clamp((int)Math.Round(halfRange + bestShift), 1, profile.Length - 2);
        int leftEdge = -1, rightEdge = -1;
        double leftVal = 0, rightVal = 0;
        double derivMax = 0;
        for (int i = 1; i < profile.Length - 1; i++) derivMax = Math.Max(derivMax, Math.Abs(deriv[i]));
        double significant = derivMax * 0.35;
        for (int i = center - 1; i >= 1; i--)
            if (Math.Abs(deriv[i]) > significant && Math.Abs(deriv[i]) >= Math.Abs(deriv[i - 1]) && Math.Abs(deriv[i]) >= Math.Abs(deriv[i + 1]))
            { leftEdge = i; leftVal = Math.Abs(deriv[i]); break; }
        for (int i = center + 1; i < profile.Length - 1; i++)
            if (Math.Abs(deriv[i]) > significant && Math.Abs(deriv[i]) >= Math.Abs(deriv[i - 1]) && Math.Abs(deriv[i]) >= Math.Abs(deriv[i + 1]))
            { rightEdge = i; rightVal = Math.Abs(deriv[i]); break; }
        //260724Cl: 左右エッジの実測オフセット (補正後中心線基準、非対称のまま保持)
        double leftOffset, rightOffset, widthWork, widthQuality;
        if (pairApplied) //260724Cl: ローブ対バンドはローブ極値 (=Kikuchi 線対) を両縁とする (勾配走査は反対称中心の急勾配で幅過小になる)
        {
            double oa = lobeA - (halfRange + bestShift), ob = lobeB - (halfRange + bestShift);
            leftOffset = Math.Min(oa, ob); rightOffset = Math.Max(oa, ob);
            widthWork = rightOffset - leftOffset;
            widthQuality = lobeAAmp > 0 && lobeBAmp > 0 ? Math.Min(lobeAAmp, lobeBAmp) / Math.Max(lobeAAmp, lobeBAmp) : 0;
        }
        else if (leftEdge >= 0 && rightEdge >= 0)
        {
            leftOffset = leftEdge - (halfRange + bestShift);   //負値
            rightOffset = rightEdge - (halfRange + bestShift); //正値
            widthWork = rightOffset - leftOffset;
            //勾配の鋭さ (プロファイル振幅比) を品質に
            var (smoothMin, smoothMax) = smooth.MinMax(); //260724Cl (/simplify): Max()+Min() の 2 走査 → SimdLinq MinMax 1 走査
            double amp = smoothMax - smoothMin;
            widthQuality = amp > 1E-9 ? Math.Clamp((leftVal + rightVal) / amp, 0, 1) : 0;
        }
        else { leftOffset = -initialWidth / 2; rightOffset = initialWidth / 2; widthWork = initialWidth; widthQuality = 0; }
        #endregion

        #region 線検証 (260724Cl 追加): 狭窓アンカー多数決 — 偽陽性 (誤配置の中心線) の抑制
        //シフト後中心 ±0.75×実測幅 の狭い台座で、7 アンカー各々の横断プロファイルが平均形状 (leave-one-out) と相関するかを見る。
        //Radon/butterfly は「複数バンドの明るい交差点を通るだけの線」(アンカー間で形状が無相関) や
        //「強いバンドの平行オフセット幽霊」(狭窓中心が無構造) にも高スコアを与えるため、独立した検証が要る。
        //窓を広げると幽霊線でも窓内の実バンドで相関が出てしまう — 狭窓 (実測幅比例) が本質
        const int nVerify = 7;
        //260724Cl: 上限 12 の絶対クランプが本質 — 幅推定が壊れた幽霊線 (エッジ不検出で initialWidth 28 等へフォールバック) は
        //vr が大きいと窓内に隣の実バンドの裾が入り「一貫した勾配」として射影 ≈1 になってしまう
        //int vr = Math.Clamp((int)Math.Ceiling(widthWork * 0.75), 4, halfRange);
        //260724Cl: 下限 6 — 極小窓 (±4) は「どんな局所勾配でも一貫」になり判別力を失う。細バンドも背景を含めた形状で照合する。WS スケール
        int vr = Math.Clamp((int)Math.Ceiling(widthWork * 0.75), (int)Math.Round(6 * WS), Math.Min((int)Math.Round(12 * WS), halfRange));
        var vProf = new double[nVerify, 2 * vr + 1];
        for (int a = 0; a < nVerify; a++)
        {
            double t = t0 + (t1 - t0) * (a + 0.5) / nVerify;
            double ax = px + t * dx, ay = py + t * dy;
            for (int s = -vr; s <= vr; s++)
            {
                vProf[a, s + vr] = double.NaN;
                double x = ax + s * cosT, y = ay + s * sinT;
                int x0 = (int)x, y0 = (int)y;
                if ((uint)x0 >= (uint)(w - 1) || (uint)y0 >= (uint)(h - 1)) continue;
                int idx = y0 * w + x0;
                if (!valid[idx]) continue;
                double fx = x - x0, fy = y - y0;
                vProf[a, s + vr] = (work[idx] * (1 - fx) + work[idx + 1] * fx) * (1 - fy) + (work[idx + w] * (1 - fx) + work[idx + w + 1] * fx) * fy;
            }
        }
        var vMean = new double[2 * vr + 1]; var vCount = new int[2 * vr + 1];
        for (int i = 0; i < vMean.Length; i++)
        {
            double s2 = 0; int n = 0;
            for (int a = 0; a < nVerify; a++) if (!double.IsNaN(vProf[a, i])) { s2 += vProf[a, i]; n++; }
            vMean[i] = n > 0 ? s2 / n : 0; vCount[i] = n;
        }
        //260724Cl: 指標は Pearson 相関でなく「射影係数」 c_a = ⟨P̃_a, M̃⟩/⟨M̃, M̃⟩ (両者 zero-mean、M̃ は LOO 平均形状)。
        //Pearson は交差バンドの混入がアンカーの分散を膨らませて実バンドの値まで潰す (busy 領域を通る正解バンドが 0 に落ちた)。
        //射影は M̃ に直交する混入に不感 — 「平均形状がそのアンカーにどれだけ含まれるか」(実バンド ≈1、幽霊/ジャンク ≈0±)
        var corrs = new List<double>();
        double meanRms = 0; int meanRmsN = 0;
        for (int a = 0; a < nVerify; a++)
        {
            double spm = 0, smm = 0, sp = 0, sm = 0; int n = 0;
            for (int i = 0; i < vMean.Length; i++)
            {
                double v = vProf[a, i];
                if (double.IsNaN(v) || vCount[i] < 2) continue;
                double m = (vCount[i] * vMean[i] - v) / (vCount[i] - 1); //自分を除いた他アンカー平均 (LOO)
                sp += v; sm += m; spm += v * m; smm += m * m; n++;
            }
            if (n < vMean.Length * 0.7) continue; //画像外にはみ出したアンカーは統計から除外
            double cov = spm - sp * sm / n, vm = smm - sm * sm / n;
            if (vm / n > 1E-12) { corrs.Add(cov / vm); meanRms += Math.Sqrt(vm / n); meanRmsN++; }
        }
        double medianCorr = 0, fracGood = 0;
        if (corrs.Count >= 4)
        {
            var srt = corrs.OrderBy(c => c).ToList();
            medianCorr = srt[srt.Count / 2];
            fracGood = (double)corrs.Count(c => c > 0.25) / corrs.Count;
        }
        meanRms = meanRmsN > 0 ? meanRms / meanRmsN : 0; //LOO 平均形状の振幅 (無構造の幽霊は小)

        //260724Cl: 平均形状の線形フィット R² — 「バンドの縁に沿って走る幽霊」(edge-runner) は窓内の M が単調勾配 ≈ 直線で R²→1。
        //真のバンド中心の M は極値 (対称バンド) か S 字 (e/d 対 — こちらは pair ゲート承認済みなので除外判定しない) を持つ
        double linR2 = 0;
        {
            double sxs = 0, sys = 0, sxx2 = 0, sxy2 = 0, syy2 = 0; int n2 = 0;
            for (int i = 0; i < vMean.Length; i++)
            {
                if (vCount[i] < 2) continue;
                sxs += i; sys += vMean[i]; sxx2 += (double)i * i; sxy2 += i * vMean[i]; syy2 += vMean[i] * vMean[i]; n2++;
            }
            if (n2 >= 5)
            {
                double cov = sxy2 - sxs * sys / n2, vx = sxx2 - sxs * sxs / n2, vy = syy2 - sys * sys / n2;
                if (vx > 1E-12 && vy > 1E-12) linR2 = cov * cov / (vx * vy);
            }
        }
        debugLog?.Invoke($"verify θ={theta * 180 / Math.PI,6:f1} ρ={rho + bestShift,6:f1} vr={vr,2} nA={corrs.Count} med={medianCorr,5:f2} frac={fracGood:f2} rms={meanRms:f2} r2={linR2:f2} c=[{string.Join(",", corrs.Select(c => c.ToString("f2")))}]");
        #endregion

        #region native 座標へ変換
        //work ピクセル中心 (x,y) → native (col,row): col = (x - cx)/scale·? — work は縮小画像なので native = work/scale (ピクセル中心規約で近似)
        double toNative = 1.0 / scale;
        (double Col, double Row) ToNative(double x, double y) => ((x + 0.5) * toNative - 0.5, (y + 0.5) * toNative - 0.5);

        var centerAnchors = anchorsT.Select(t => ToNative(px + t * dx, py + t * dy)).ToArray();
        var edges = new List<(double Col, double Row)>();
        foreach (var t in anchorsT)
        {
            //260724Cl: 実測した左右オフセットを別々に適用 (非対称保持)
            edges.Add(ToNative(px + t * dx + leftOffset * cosT, py + t * dy + leftOffset * sinT));
            edges.Add(ToNative(px + t * dx + rightOffset * cosT, py + t * dy + rightOffset * sinT));
        }

        //中心線係数 (native): 法線 (cosT,sinT) は等方スケールで不変。ρ を native へ
        //work: (x-cx)cosT + (y-cy)sinT = ρ。x_native=(x+0.5)/scale-0.5 → 代入で ρ_native = ρ/scale + (中心シフト)
        var a0 = centerAnchors[0]; var a1 = centerAnchors[^1];
        double lineA = -(a1.Row - a0.Row), lineB = a1.Col - a0.Col;
        double lineNorm = Math.Sqrt(lineA * lineA + lineB * lineB);
        lineA /= lineNorm; lineB /= lineNorm;
        double lineC = -(lineA * a0.Col + lineB * a0.Row);

        return new EbsdBand
        {
            LineA = lineA, LineB = lineB, LineC = lineC,
            CenterAnchors = centerAnchors,
            EdgePoints = [.. edges],
            WidthPx = widthWork * toNative,
            Score = score,
            //CenterQuality = Math.Clamp(score / 20.0, 0, 1),
            CenterQuality = Math.Clamp(medianCorr, 0, 1), //260724Cl: スコア比例の仮値 → 狭窓アンカー射影の中央値へ。偽線ほど低い
            VerifyFrac = fracGood, EdgeLinearity = linR2, IsLobePair = pairApplied, //260724Cl
            WidthQuality = widthQuality,
            ThetaDeg = theta * 180 / Math.PI,
            RhoWorkPx = rho + bestShift, //260724Cl: 対称相関の中心補正を表示用派生値にも反映 (Codex 指摘)
        };
        #endregion
    }

    #region 画像ユーティリティ (縮小・ガウシアン・正規化)

    static double[] Downsample(double[] src, int sw, int sh, bool[] srcValid, int dw, int dh, out bool[] dstValid)
    {
        var dst = new double[dw * dh];
        dstValid = new bool[dw * dh];
        double sx = (double)sw / dw, sy = (double)sh / dh;
        for (int y = 0; y < dh; y++)
        {
            int y0 = (int)(y * sy), y1 = Math.Min(sh, (int)Math.Ceiling((y + 1) * sy));
            for (int x = 0; x < dw; x++)
            {
                int x0 = (int)(x * sx), x1 = Math.Min(sw, (int)Math.Ceiling((x + 1) * sx));
                double sum = 0; int n = 0;
                for (int yy = y0; yy < y1; yy++)
                    for (int xx = x0; xx < x1; xx++)
                    {
                        int i = yy * sw + xx;
                        if (srcValid == null || srcValid[i]) { sum += src[i]; n++; }
                    }
                if (n > 0) { dst[y * dw + x] = sum / n; dstValid[y * dw + x] = true; }
            }
        }
        return dst;
    }

    /// <summary>分離ガウシアン (EbsdPatternScorer の背景除去からも利用)。260724Cl internal 公開</summary>
    internal static double[] GaussianBlurGrid(double[] src, bool[] valid, int w, int h, double sigma) => GaussianBlur(src, valid, w, h, sigma);

    /// <summary>分離ガウシアン (invalid 画素は重みから除外する normalized convolution)</summary>
    static double[] GaussianBlur(double[] src, bool[] valid, int w, int h, double sigma)
    {
        int half = Math.Max(1, (int)Math.Ceiling(3 * sigma));
        var kernel = new double[2 * half + 1];
        for (int i = -half; i <= half; i++) kernel[i + half] = Gauss(i, sigma);

        //var tmp = new double[w * h]; var tmpW = new double[w * h]; //260725Ch 変更前: 呼び出しごとに大配列 2 本を確保
        int length = checked(w * h);
        var tmp = ArrayPool<double>.Shared.Rent(length); //260725Ch
        double[] tmpW = null; //260725Ch: 2本目のRent自体が失敗しても1本目をfinallyで返す
        try
        {
            tmpW = ArrayPool<double>.Shared.Rent(length); //260725Ch
            System.Threading.Tasks.Parallel.For(0, h, y =>
            {
                for (int x = 0; x < w; x++)
                {
                    double s = 0, wsum = 0;
                    for (int k = -half; k <= half; k++)
                    {
                        int xx = x + k;
                        if ((uint)xx >= (uint)w) continue;
                        int i = y * w + xx;
                        if (!valid[i]) continue;
                        s += src[i] * kernel[k + half]; wsum += kernel[k + half];
                    }
                    tmp[y * w + x] = s; tmpW[y * w + x] = wsum;
                }
            });
            var dst = new double[length];
            System.Threading.Tasks.Parallel.For(0, h, y =>
            {
                for (int x = 0; x < w; x++)
                {
                    double s = 0, wsum = 0;
                    for (int k = -half; k <= half; k++)
                    {
                        int yy = y + k;
                        if ((uint)yy >= (uint)h) continue;
                        s += tmp[yy * w + x] * kernel[k + half]; wsum += tmpW[yy * w + x] * kernel[k + half];
                    }
                    dst[y * w + x] = wsum > 1E-12 ? s / wsum : 0;
                }
            });
            return dst;
        }
        finally
        {
            //260725Ch: Parallel.For が例外終了しても、全ワーカー停止後に必ずプールへ返す
            ArrayPool<double>.Shared.Return(tmp);
            if (tmpW != null) ArrayPool<double>.Shared.Return(tmpW); //260725Ch: 2本目のRent失敗時はnull
        }
    }

    static void Normalize(double[] img, bool[] valid, double clipSigma)
    {
        double mean = 0; int n = 0;
        for (int i = 0; i < img.Length; i++) if (valid[i]) { mean += img[i]; n++; }
        if (n == 0) return;
        mean /= n;
        double var = 0;
        for (int i = 0; i < img.Length; i++) if (valid[i]) { double d = img[i] - mean; var += d * d; }
        double std = Math.Sqrt(var / n);
        if (std < 1E-12) std = 1;
        for (int i = 0; i < img.Length; i++)
            img[i] = valid[i] ? Math.Clamp((img[i] - mean) / std, -clipSigma, clipSigma) : 0;
    }

    #endregion
}
