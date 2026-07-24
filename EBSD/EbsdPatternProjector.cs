#region using
using System;
using System.Linq;
using System.Threading.Tasks;
using V3 = OpenTK.Mathematics.Vector3d;
#endregion

namespace Crystallography;

/// <summary>
/// MasterPattern から指定方位のシミュレーション EBSD パターンを検出器 native グリッド (縮小可) に投影する純関数プロジェクター。260724Cl 追加。
/// 2 段構成: コンストラクタで検出器グリッドの視線 (試料系、方位非依存) をキャッシュし、Project で回転適用+Rosca-Lambert+補間のみ行う。
/// 座標規約 (視線 = -P 方向) は EbsdDetectorGeometry.PixelToSampleDirection = FormEBSD.BuildEbsdLookupTable と同一。single slice 専用 (速度優先、§6.15)。
/// </summary>
public sealed class EbsdPatternProjector
{
    readonly V3[] raysSample;
    public int Width { get; }
    public int Height { get; }

    /// <summary>geometry のピクセルグリッドを rasterW×rasterH に縮小したグリッドで視線をキャッシュする</summary>
    public EbsdPatternProjector(EbsdDetectorGeometry geometry, int rasterW, int rasterH)
    {
        Width = rasterW; Height = rasterH;
        raysSample = new V3[rasterW * rasterH];
        //260724Cl (/simplify): 各ピクセル独立なので並列化 (幾何較正では評価毎に本コンストラクタが再構築されるため、逐次だと 1 Project 相当の逐次コストが毎評価に乗っていた)
        Parallel.For(0, rasterH, r =>
        {
            for (int c = 0; c < rasterW; c++)
            {
                double col = (c + 0.5) / rasterW * geometry.WidthPx - 0.5;
                double row = (r + 0.5) / rasterH * geometry.HeightPx - 0.5;
                raysSample[r * rasterW + c] = geometry.PixelToSampleDirection(col, row);
            }
        });
    }

    #region 面内回転分解プロジェクション (辞書総当たり用、square 格子専用) 260725Cl 追加
    //R(di,φ)=r0·Rz(φ) の構造を利用: 結晶系視線 v = Rz(−φ)·(r0⁻¹·d) なので、u=r0⁻¹·d の Lambert ディスク極座標
    //(方位角 θ0・半径由来の ra/rb) と半球フラグを球点毎に 1 回だけ計算し、面内 φ 毎は θ=θ0−φ の sector 折り返し+
    //バイリニア補間のみにする (3×3 回転積・sqrt・atan を全除去。Codex 裁定 260725: Lambert 後の (a,b) 2D 回転は不可、
    //ディスク極座標の再利用が正解)。SphereToRoscaLambertSquare の Shirley 逆変換と数学的に等価 (atan(tanθ')=局所角 t)。

    /// <summary>球点回転 r0 の面内共通量を計算する。theta0/ra/rb/neg は呼び出し側が確保する長さ Width×Height のバッファ。
    /// ra = ディスク半径×√π/2 (sector 支配軸の座標)、rb = 4·ra/π (直交軸の係数)。260725Cl 追加</summary>
    public void PrepareSpherePoint(Matrix3D r0, double[] theta0, double[] ra, double[] rb, bool[] neg)
    {
        var ri = r0.Inverse();
        double sqrtPiHalf = Math.Sqrt(Math.PI) / 2;
        for (int i = 0; i < raysSample.Length; i++)
        {
            var d = raysSample[i];
            double ux = ri.E11 * d.X + ri.E12 * d.Y + ri.E13 * d.Z;
            double uy = ri.E21 * d.X + ri.E22 * d.Y + ri.E23 * d.Z;
            double uz = ri.E31 * d.X + ri.E32 * d.Y + ri.E33 * d.Z;
            double len = Math.Sqrt(ux * ux + uy * uy + uz * uz);
            neg[i] = uz < 0;
            double az = len < 1E-15 ? 1 : Math.Abs(uz) / len; //len=0 は中心特異点 — ra=0 で (0,0) に落ちる
            double radialScale = Math.Sqrt(Math.Max(0, (1 + az) / 2)); //SphereToRoscaLambertSquare の z>=0 枝 (|z| を渡すため常にこちら)
            double rxy = len < 1E-15 ? 0 : Math.Sqrt(ux * ux + uy * uy) / len;
            double rDisk = radialScale < 1E-15 ? 0 : rxy / radialScale;
            ra[i] = rDisk * sqrtPiHalf;
            rb[i] = 4 * ra[i] / Math.PI;
            theta0[i] = Math.Atan2(uy, ux);
        }
    }

    /// <summary>PrepareSpherePoint 済みの球点について面内回転角 phi のパターンを output へ書き込む (完全逐次、square 格子専用)。260725Cl 追加</summary>
    public void ProjectInPlane(MasterPattern mp, double phi, float[] posPlane, float[] negPlane, double[] theta0, double[] ra, double[] rb, bool[] neg, double[] output)
    {
        int gs = mp.GridSize;
        bool hasPos = posPlane is { Length: > 0 }, hasNeg = negPlane is { Length: > 0 };
        const double halfPi = Math.PI / 2, quarterPi = Math.PI / 4;
        for (int i = 0; i < output.Length; i++)
        {
            bool n = neg[i];
            if (n ? !hasNeg : !hasPos) { output[i] = 0; continue; }
            double th = theta0[i] - phi;
            int s = (int)Math.Floor((th + quarterPi) / halfPi); //支配軸 sector (0=+a,1=+b,2=−a,3=−b)、t = sector 内局所角 ∈ [−π/4, π/4)
            double t = th - s * halfPi;
            double a, b;
            switch (((s % 4) + 4) % 4)
            {
                case 0: a = ra[i]; b = rb[i] * t; break;
                case 1: b = ra[i]; a = -rb[i] * t; break;
                case 2: a = -ra[i]; b = -rb[i] * t; break;
                default: b = -ra[i]; a = rb[i] * t; break;
            }
            output[i] = MasterPattern.InterpolatePlaneSquare(n ? negPlane : posPlane, gs, a, b);
        }
    }
    #endregion

    /// <summary>回転 rotation (crystal→sample) のパターンを output (Width×Height) へ書き込む。posPlane/negPlane = MasterPattern.GetPlane の単一スライス。
    /// parallel=false で行ループを逐次実行 (辞書総当たりのような方位単位で並列化する呼び出し向け。小ラスターでは行並列のオーバーヘッドが支配的)。260724Cl シグネチャ変更 (parallel 追加)</summary>
    //260724Cl 旧: public void Project(MasterPattern mp, Matrix3D rotation, float[] posPlane, float[] negPlane, double[] output)
    public void Project(MasterPattern mp, Matrix3D rotation, float[] posPlane, float[] negPlane, double[] output, bool parallel = true)
    {
        var ri = rotation.Inverse();
        int gs = mp.GridSize;
        bool isHex = mp.GridType == MasterPattern.Types.Hexagonal;
        bool hasPos = posPlane is { Length: > 0 }, hasNeg = negPlane is { Length: > 0 };

        if (parallel)
            Parallel.For(0, Height, r => ProjectRow(r));
        else
            for (int r = 0; r < Height; r++) ProjectRow(r);

        void ProjectRow(int r)
        {
            for (int c = 0; c < Width; c++)
            {
                int i = r * Width + c;
                var d = raysSample[i];
                //結晶系へ (Ri・d)
                double dx = ri.E11 * d.X + ri.E12 * d.Y + ri.E13 * d.Z;
                double dy = ri.E21 * d.X + ri.E22 * d.Y + ri.E23 * d.Z;
                double dz = ri.E31 * d.X + ri.E32 * d.Y + ri.E33 * d.Z;

                bool posZ = dz >= 0;
                float[] plane = posZ ? posPlane : negPlane;
                if (posZ ? !hasPos : !hasNeg) { output[i] = 0; continue; }

                if (isHex) //六方格子 (FormEBSD.BuildEbsdLookupTable と同式)
                {
                    double invLen = 1.0 / Math.Sqrt(dx * dx + dy * dy + dz * dz);
                    var (hx, hy) = MasterPattern.SphereToRoscaLambertHex(dx * invLen, dy * invLen, Math.Abs(dz) * invLen);
                    MasterPattern.GetHexBarycentricLookup(hx, hy, gs, out int idx0, out int idx1, out int idx2, out float bw0, out float bw1, out float bw2);
                    output[i] = bw0 * plane[idx0] + bw1 * plane[idx1] + bw2 * plane[idx2];
                }
                else //正方格子 (同式)
                {
                    //260724Cl (/simplify): インライン展開していた Lambert 逆写像+バイリニアを hex 分岐と同方針で共有ヘルパへ委譲
                    //(数式は同一: 旧 edgeRadius=√(π/2·(1−|z|)) は SphereToRoscaLambertSquare の radialScale 経由の式と等価。ヘルパ内で正規化される)
                    var (a, b) = MasterPattern.SphereToRoscaLambertSquare(dx, dy, Math.Abs(dz));
                    output[i] = MasterPattern.InterpolatePlaneSquare(plane, gs, a, b);
                }
            }
        }
    }
}

/// <summary>実測/シミュレーションパターンの比較スコア (masked ZNCC) と最適化ユーティリティ。260724Cl 追加</summary>
public static class EbsdPatternScorer
{
    /// <summary>実測画像を targetLongSide へ box 縮小 → 広域ガウシアン背景で除算 → zero-mean/unit-variance 化した参照配列を返す (最適化前に 1 回だけ準備)。
    /// 260724Cl: 背景除算 (蛍光体照明の乗算勾配の除去) が無いと ZNCC が方位でなく照明勾配に引かれ、特に PC/DD 較正で致命的 (Codex 指摘)</summary>
    public static (double[] Data, int W, int H) PrepareReference(double[] values, int width, int height, int targetLongSide = 160)
    {
        var (dst, w, h) = Downsample(values, width, height, targetLongSide); //260724Cl: 縮小部を Downsample へ抽出 (PrepareReferenceRobust と共用)
        //広域ガウシアン (σ=0.1×短辺) 背景で除算 (EbsdBandDetector の前処理と同方式)
        var validAll = new bool[w * h];
        Array.Fill(validAll, true);
        var bg = EbsdBandDetector.GaussianBlurGrid(dst, validAll, w, h, 0.10 * Math.Min(w, h));
        double floor = Math.Max(1E-10, dst.Average() * 0.05);
        for (int i = 0; i < dst.Length; i++)
            dst[i] /= Math.Max(bg[i], floor);
        NormalizeInPlace(dst);
        return (dst, w, h);
    }

    /// <summary>実測画像を box 縮小し robust 前処理 (RobustPreprocess) を掛けた参照配列を返す。
    /// 方位候補の複合ランク (Radon z + ZNCC) 用 — シミュレーション側にも同じ RobustPreprocess を掛けて比較する (Codex 裁定 260724)。260724Cl 追加</summary>
    //260725Cl シグネチャ変更 (dogSigma1/dogSigma2 追加): フル解像度の公正比較 (σ を解像度比例スケール) 用。既定値は従来と同一
    //旧: public static (double[] Data, int W, int H) PrepareReferenceRobust(double[] values, int width, int height, int targetLongSide = 160)
    public static (double[] Data, int W, int H) PrepareReferenceRobust(double[] values, int width, int height, int targetLongSide = 160, double dogSigma1 = 1.5, double dogSigma2 = 6.0)
    {
        var (dst, w, h) = Downsample(values, width, height, targetLongSide);
        return (RobustPreprocess(dst, w, h, dogSigma1, dogSigma2), w, h);
    }

    /// <summary>box 縮小 (targetLongSide = 長辺の目標 px)。260724Cl 追加 (PrepareReference からの抽出。EbsdDictionaryIndexer と共用のため internal)</summary>
    internal static (double[] Data, int W, int H) Downsample(double[] values, int width, int height, int targetLongSide)
    {
        double scale = Math.Min(1.0, (double)targetLongSide / Math.Max(width, height));
        int w = Math.Max(8, (int)Math.Round(width * scale)), h = Math.Max(8, (int)Math.Round(height * scale));
        var dst = new double[w * h];
        double sx = (double)width / w, sy = (double)height / h;
        for (int y = 0; y < h; y++)
        {
            int y0 = (int)(y * sy), y1 = Math.Min(height, (int)Math.Ceiling((y + 1) * sy));
            for (int x = 0; x < w; x++)
            {
                int x0 = (int)(x * sx), x1 = Math.Min(width, (int)Math.Ceiling((x + 1) * sx));
                double sum = 0; int n = 0;
                for (int yy = y0; yy < y1; yy++)
                    for (int xx = x0; xx < x1; xx++) { sum += values[yy * width + xx]; n++; }
                dst[y * w + x] = n > 0 ? sum / n : 0;
            }
        }
        return (dst, w, h);
    }

    /// <summary>
    /// robust ZNCC 前処理 (Codex 裁定 260724): 広域 bg 除算 (log-ratio, σ=0.1×短辺) → DoG σ1=1.5/σ2=6 → 標準化 → tanh(z/3) ソフトクリップ → 再標準化。
    /// 実測とシミュレーション投影の両方に同一処理を掛けることで、動力学単一スライスの heavy-tailed な生強度分布 (zone axis 明点が分散を支配し
    /// ZNCC が正解方位で偽方位に負ける) を等質化する。src は非破壊。260724Cl 追加
    /// </summary>
    //260725Cl シグネチャ変更 (dogSigma1/dogSigma2 追加): フル解像度の公正比較 (σ を解像度比例スケール) 用。既定値は従来と同一
    //旧: public static double[] RobustPreprocess(double[] src, int w, int h)
    public static double[] RobustPreprocess(double[] src, int w, int h, double dogSigma1 = 1.5, double dogSigma2 = 6.0)
    {
        var validAll = new bool[w * h];
        Array.Fill(validAll, true);
        var bg = EbsdBandDetector.GaussianBlurGrid(src, validAll, w, h, 0.10 * Math.Min(w, h));
        double floor = Math.Max(1E-10, src.Average() * 0.05);
        var v = new double[w * h];
        for (int i = 0; i < v.Length; i++) v[i] = Math.Log(Math.Max(src[i], floor * 0.01) / Math.Max(bg[i], floor));
        var g1 = EbsdBandDetector.GaussianBlurGrid(v, validAll, w, h, dogSigma1);
        var g2 = EbsdBandDetector.GaussianBlurGrid(v, validAll, w, h, dogSigma2);
        for (int i = 0; i < v.Length; i++) v[i] = g1[i] - g2[i];
        NormalizeInPlace(v);
        for (int i = 0; i < v.Length; i++) v[i] = Math.Tanh(v[i] / 3);
        NormalizeInPlace(v);
        return v;
    }

    /// <summary>
    /// RobustPreprocess の高速版 (260724Cl 追加、辞書総当たりの高速化)。ガウシアンを running-box 3 連 (分散一致) で近似し、
    /// 呼び出し側が渡す scratch (長さ w·h × 3 本) を再利用してアロケーションゼロ・完全逐次で実行する
    /// (呼び出し側が方位単位で並列化する前提。GaussianBlurGrid 内蔵の Parallel.For との入れ子競合を解消)。
    /// 結果は dst に書く。数値はガウシアン版と厳密一致しないが帯域特性は同等 (採否はベンチの refHit/順位一致で判定)。
    /// </summary>
    public static void RobustPreprocessFast(double[] src, int w, int h, double[] dst, double[] tmp1, double[] tmp2)
    {
        if (box3Scratch == null || box3Scratch.Length < w * h) box3Scratch = new double[w * h];
        var tmp3 = box3Scratch;
        //①広域背景 box3(σ=0.1×短辺): src → tmp1 (tmp2 作業)
        Box3Seq(src, tmp1, tmp2, w, h, 0.10 * Math.Min(w, h));
        double mean = 0;
        foreach (var x in src) mean += x;
        double floor = Math.Max(1E-10, mean / src.Length * 0.05);
        //②log-ratio → dst
        for (int i = 0; i < dst.Length; i++) dst[i] = Math.Log(Math.Max(src[i], floor * 0.01) / Math.Max(tmp1[i], floor));
        //③DoG: g1(σ1.5) = dst→tmp1、g2(σ6) = dst→tmp2 (dst は両方の入力なので温存)
        Box3Seq(dst, tmp1, tmp2, w, h, 1.5);
        Box3Seq(dst, tmp2, tmp3, w, h, 6.0);
        for (int i = 0; i < dst.Length; i++) dst[i] = tmp1[i] - tmp2[i];
        //④標準化 → tanh(z/3) → 再標準化
        NormalizeInPlace(dst);
        for (int i = 0; i < dst.Length; i++) dst[i] = Math.Tanh(dst[i] / 3);
        NormalizeInPlace(dst);
    }

    [ThreadStatic] static double[] box3Scratch; //260724Cl: RobustPreprocessFast の第 4 バッファ (スレッドローカル再利用)

    /// <summary>box blur 3 連 (ガウシアン分散一致近似、完全逐次)。src → dst (src は不変)。work は作業バッファ。260724Cl 追加</summary>
    static void Box3Seq(double[] src, double[] dst, double[] work, int w, int h, double sigma)
    {
        int r = Math.Max(1, (int)Math.Round((Math.Sqrt(4 * sigma * sigma + 1) - 1) / 2)); //3 連の合成分散 3(w²−1)/12 = σ² となる box 幅
        //260725Cl: 境界正規化 1/n を位置別に事前計算し BoxPassSeq 内の毎画素除算 (~18 回/画素) を乗算化 (prof: box が前処理の 52%)
        Span<double> invX = w <= 2048 ? stackalloc double[w] : new double[w];
        Span<double> invY = h <= 2048 ? stackalloc double[h] : new double[h];
        for (int x = 0; x < w; x++) invX[x] = 1.0 / (Math.Min(x + r, w - 1) - Math.Max(x - r, 0) + 1);
        for (int y = 0; y < h; y++) invY[y] = 1.0 / (Math.Min(y + r, h - 1) - Math.Max(y - r, 0) + 1);
        BoxPassSeq(src, dst, w, h, r, invX, invY);
        BoxPassSeq(dst, work, w, h, r, invX, invY);
        BoxPassSeq(work, dst, w, h, r, invX, invY);
    }

    /// <summary>running box 平均 1 回分 (横+縦、境界は有効画素数で正規化、完全逐次)。src → dst (src 不変)。260724Cl 追加</summary>
    //260725Cl シグネチャ変更 (invX/invY 追加): sum/n 除算 → sum·(1/n) 乗算 (テーブルは Box3Seq が半径別に事前計算)
    //旧: static void BoxPassSeq(double[] src, double[] dst, int w, int h, int radius)
    static void BoxPassSeq(double[] src, double[] dst, int w, int h, int radius, ReadOnlySpan<double> invX, ReadOnlySpan<double> invY)
    {
        for (int y = 0; y < h; y++) //横パス: src → dst
        {
            int row = y * w;
            double sum = 0;
            for (int x = 0; x <= Math.Min(radius, w - 1); x++) sum += src[row + x];
            for (int x = 0; x < w; x++)
            {
                //dst[row + x] = sum / n; //260725Cl 変更前 (n は running カウント)
                dst[row + x] = sum * invX[x];
                int add = x + radius + 1, rem = x - radius;
                if (add < w) sum += src[row + add];
                if (rem >= 0) sum -= src[row + rem];
            }
        }
        Span<double> col = h <= 2048 ? stackalloc double[2048] : new double[h]; //縦パス用の一時列
        for (int x = 0; x < w; x++) //縦パス: dst → dst
        {
            double sum = 0;
            for (int y = 0; y <= Math.Min(radius, h - 1); y++) sum += dst[y * w + x];
            for (int y = 0; y < h; y++)
            {
                //col[y] = sum / n; //260725Cl 変更前
                col[y] = sum * invY[y];
                int add = y + radius + 1, rem = y - radius;
                if (add < h) sum += dst[add * w + x];
                if (rem >= 0) sum -= dst[rem * w + x];
            }
            for (int y = 0; y < h; y++) dst[y * w + x] = col[y];
        }
    }

    /// <summary>zero-mean/unit-variance 化 (in place)</summary>
    public static void NormalizeInPlace(double[] data)
    {
        double mean = 0;
        foreach (var v in data) mean += v;
        mean /= data.Length;
        double var = 0;
        foreach (var v in data) { double d = v - mean; var += d * d; }
        double std = Math.Sqrt(var / data.Length);
        if (std < 1E-12) std = 1;
        for (int i = 0; i < data.Length; i++) data[i] = (data[i] - mean) / std;
    }

    /// <summary>zero-normalized cross correlation。a は正規化済み参照、b は生 (内部で正規化)</summary>
    public static double Zncc(double[] normalizedRef, double[] pattern)
    {
        double mean = 0;
        foreach (var v in pattern) mean += v;
        mean /= pattern.Length;
        double var = 0, dot = 0;
        for (int i = 0; i < pattern.Length; i++)
        {
            double d = pattern[i] - mean;
            var += d * d;
            dot += normalizedRef[i] * d;
        }
        double std = Math.Sqrt(var / pattern.Length);
        return std < 1E-12 ? 0 : dot / (pattern.Length * std);
    }

    /// <summary>
    /// 簡易 Nelder-Mead (初期ステップ明示・下降単体法)。objective を最小化する。260724Cl 追加
    /// MathNet の NelderMeadSimplex は初期シンプレックス制御が弱いため自前実装 (数変数・数百評価の用途限定)。
    /// </summary>
    public static (double[] Best, double Value, int Evaluations) NelderMead(Func<double[], double> objective, double[] start, double[] step, int maxEval = 400, double tol = 1E-5)
    {
        int n = start.Length;
        var simplex = new double[n + 1][];
        var values = new double[n + 1];
        int eval = 0;
        for (int i = 0; i <= n; i++)
        {
            simplex[i] = (double[])start.Clone();
            if (i > 0) simplex[i][i - 1] += step[i - 1];
            values[i] = objective(simplex[i]); eval++;
        }
        while (eval < maxEval)
        {
            Array.Sort(values, simplex);
            if (Math.Abs(values[n] - values[0]) < tol) break;

            var centroid = new double[n];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++) centroid[j] += simplex[i][j] / n;

            double[] Blend(double coeff)
            {
                var p = new double[n];
                for (int j = 0; j < n; j++) p[j] = centroid[j] + coeff * (simplex[n][j] - centroid[j]);
                return p;
            }
            double Eval(double[] p) { double v = objective(p); eval++; return double.IsNaN(v) ? double.PositiveInfinity : v; } //260724Cl: NaN は +∞ 扱い

            var reflected = Blend(-1);
            double fr = Eval(reflected);
            if (fr < values[0])
            {
                var expanded = Blend(-2);
                double fe = Eval(expanded);
                if (fe < fr) { simplex[n] = expanded; values[n] = fe; }
                else { simplex[n] = reflected; values[n] = fr; }
            }
            else if (fr < values[n - 1]) { simplex[n] = reflected; values[n] = fr; }
            else if (fr < values[n]) //260724Cl: outside contraction (標準 Nelder-Mead。旧実装はこの場合も inside contraction だった)
            {
                var contracted = Blend(-0.5);
                double fc = Eval(contracted);
                if (fc <= fr) { simplex[n] = contracted; values[n] = fc; }
                else
                    Shrink();
            }
            else
            {
                var contracted = Blend(0.5); //inside contraction
                double fc = Eval(contracted);
                if (fc < values[n]) { simplex[n] = contracted; values[n] = fc; }
                else
                    Shrink();
            }

            void Shrink()
            {
                for (int i = 1; i <= n; i++)
                {
                    for (int j = 0; j < n; j++) simplex[i][j] = simplex[0][j] + 0.5 * (simplex[i][j] - simplex[0][j]);
                    values[i] = Eval(simplex[i]);
                }
            }
        }
        Array.Sort(values, simplex);
        return (simplex[0], values[0], eval);
    }
}
