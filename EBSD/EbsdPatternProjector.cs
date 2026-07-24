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

    /// <summary>回転 rotation (crystal→sample) のパターンを output (Width×Height) へ書き込む。posPlane/negPlane = MasterPattern.GetPlane の単一スライス</summary>
    public void Project(MasterPattern mp, Matrix3D rotation, float[] posPlane, float[] negPlane, double[] output)
    {
        var ri = rotation.Inverse();
        int gs = mp.GridSize;
        bool isHex = mp.GridType == MasterPattern.Types.Hexagonal;
        bool hasPos = posPlane is { Length: > 0 }, hasNeg = negPlane is { Length: > 0 };

        Parallel.For(0, Height, r =>
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
        });
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
    public static (double[] Data, int W, int H) PrepareReferenceRobust(double[] values, int width, int height, int targetLongSide = 160)
    {
        var (dst, w, h) = Downsample(values, width, height, targetLongSide);
        return (RobustPreprocess(dst, w, h), w, h);
    }

    /// <summary>box 縮小 (targetLongSide = 長辺の目標 px)。260724Cl 追加 (PrepareReference からの抽出)</summary>
    static (double[] Data, int W, int H) Downsample(double[] values, int width, int height, int targetLongSide)
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
    public static double[] RobustPreprocess(double[] src, int w, int h)
    {
        var validAll = new bool[w * h];
        Array.Fill(validAll, true);
        var bg = EbsdBandDetector.GaussianBlurGrid(src, validAll, w, h, 0.10 * Math.Min(w, h));
        double floor = Math.Max(1E-10, src.Average() * 0.05);
        var v = new double[w * h];
        for (int i = 0; i < v.Length; i++) v[i] = Math.Log(Math.Max(src[i], floor * 0.01) / Math.Max(bg[i], floor));
        var g1 = EbsdBandDetector.GaussianBlurGrid(v, validAll, w, h, 1.5);
        var g2 = EbsdBandDetector.GaussianBlurGrid(v, validAll, w, h, 6.0);
        for (int i = 0; i < v.Length; i++) v[i] = g1[i] - g2[i];
        NormalizeInPlace(v);
        for (int i = 0; i < v.Length; i++) v[i] = Math.Tanh(v[i] / 3);
        NormalizeInPlace(v);
        return v;
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
