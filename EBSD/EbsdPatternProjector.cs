#region using
using System;
using System.Linq;
using System.Numerics; //260725Cl: Vector<double> SIMD (RobustPreprocessFast の log/tanh/正規化)
using System.Threading.Tasks;
using V3 = OpenTK.Mathematics.Vector3d;
#endregion

namespace Crystallography;

/// <summary>
/// MasterPattern から指定方位のシミュレーション EBSD パターンを検出器 native グリッド (縮小可) に投影する純関数プロジェクター。260724Cl 追加。
/// 2 段構成: コンストラクタで検出器グリッドの視線 (試料系、方位非依存) をキャッシュし、Project で回転適用+Rosca-Lambert+補間のみ行う。
/// 座標規約 (視線 = -P 方向) は EbsdDetectorGeometry.PixelToSampleDirection = EbsdPatternComposer.BuildLookupTable と同一。single slice 専用 (速度優先、§6.15)。
/// </summary>
public sealed class EbsdPatternProjector
{
    readonly V3[] raysSample;
    public int Width { get; }
    public int Height { get; }

    /// <summary>geometry のピクセルグリッドを rasterW×rasterH に縮小したグリッドで視線をキャッシュする</summary>
    public EbsdPatternProjector(EbsdDetectorGeometry geometry, int rasterW, int rasterH)
    {
        ArgumentNullException.ThrowIfNull(geometry); //260725Ch: 負寸法・積overflowを配列確保より前に明瞭化
        if (rasterW <= 0) throw new ArgumentOutOfRangeException(nameof(rasterW));
        if (rasterH <= 0) throw new ArgumentOutOfRangeException(nameof(rasterH));
        Width = rasterW; Height = rasterH;
        //raysSample = new V3[rasterW * rasterH]; //260725Ch 変更前
        raysSample = new V3[checked(rasterW * rasterH)]; //260725Ch
        Rebuild(geometry); //260920Cl: 視線の計算は Rebuild と共通化 (幾何較正が同じバッファを使い回せるようにするため)
    }

    /// <summary>260920Cl 追加: 同じラスターのまま検出器幾何だけ差し替えて視線を計算し直す。
    /// 幾何較正は 1 評価ごとに幾何が変わるが、旧実装はそのたびにコンストラクタで V3[W×H] を確保していた
    /// (フル解像度 1344×1024 では 1 評価あたり 33 MB。数千評価ぶんの GC 負荷が投影本体より重かった)。
    /// バッファを使い回すことで確保をゼロにする。数値は従来と同一 (同じ式・同じ順序)</summary>
    public void Rebuild(EbsdDetectorGeometry geometry, bool parallel = true)
    {
        ArgumentNullException.ThrowIfNull(geometry);
        int rasterW = Width, rasterH = Height;
        //260724Cl (/simplify): 各ピクセル独立なので並列化。260920Cl: 多点開始を並列に走らせるときは入れ子を避けるため parallel: false で呼ぶ
        if (parallel) Parallel.For(0, rasterH, BuildRow);
        else for (int r = 0; r < rasterH; r++) BuildRow(r);

        void BuildRow(int r)
        {
            for (int c = 0; c < rasterW; c++)
            {
                double col = (c + 0.5) / rasterW * geometry.WidthPx - 0.5;
                double row = (r + 0.5) / rasterH * geometry.HeightPx - 0.5;
                raysSample[r * rasterW + c] = geometry.PixelToSampleDirection(col, row);
            }
        }
    }

    /// <summary>公開投影APIの共通バッファ前提をホット画素ループの外で一度だけ検証する。260725Ch 追加</summary>
    int ValidateProjectionBuffers(MasterPattern mp, float[] posPlane, float[] negPlane, double[] output)
    {
        ArgumentNullException.ThrowIfNull(mp);
        ArgumentNullException.ThrowIfNull(output);
        if (output.Length != raysSample.Length) throw new ArgumentException("output.Length must equal Width * Height.", nameof(output));
        if (mp.GridSize < 2) throw new ArgumentException("MasterPattern.GridSize must be at least 2.", nameof(mp));
        int requiredPlaneLength = checked(mp.GridSize * mp.GridSize);
        if (posPlane != null && posPlane.Length < requiredPlaneLength) throw new ArgumentException("The positive master-pattern plane is too short.", nameof(posPlane));
        if (negPlane != null && negPlane.Length < requiredPlaneLength) throw new ArgumentException("The negative master-pattern plane is too short.", nameof(negPlane));
        return mp.GridSize;
    }

    #region 面内回転分解プロジェクション (辞書総当たり用、square 格子専用) 260725Cl 追加
    //R(di,φ)=r0·Rz(φ) の構造を利用: 結晶系視線 v = Rz(−φ)·(r0⁻¹·d) なので、u=r0⁻¹·d の Lambert ディスク極座標
    //(方位角 θ0・半径由来の ra/rb) と半球フラグを球点毎に 1 回だけ計算し、面内 φ 毎は θ=θ0−φ の sector 折り返し+
    //バイリニア補間のみにする (3×3 回転積・sqrt・atan を全除去。Codex 裁定 260725: Lambert 後の (a,b) 2D 回転は不可、
    //ディスク極座標の再利用が正解)。SphereToRoscaLambertSquare の Shirley 逆変換と数学的に等価 (atan(tanθ')=局所角 t)。

    /// <summary>球点回転 r0 の面内共通量を計算する。q0/ra/rb/neg は呼び出し側が確保する長さ Width×Height のバッファ。
    /// q0 = ディスク方位角×(2/π) (sector 判定を除算なしで行う正規化角)、ra = ディスク半径×√π/2 (sector 支配軸の座標)、rb = 4·ra/π (直交軸の係数)。260725Cl 追加</summary>
    //260725Cl (/simplify) 変更: theta0 (方位角そのもの) → q0=theta0·2/π を保存し、ProjectInPlane の毎画素 /halfPi 除算を除去
    public void PrepareSpherePoint(Matrix3D r0, double[] q0, double[] ra, double[] rb, bool[] neg)
    {
        ArgumentNullException.ThrowIfNull(q0); ArgumentNullException.ThrowIfNull(ra); //260725Ch: 公開作業バッファの短配列を画素ループ前に拒否
        ArgumentNullException.ThrowIfNull(rb); ArgumentNullException.ThrowIfNull(neg);
        int requiredLength = raysSample.Length;
        if (q0.Length < requiredLength || ra.Length < requiredLength || rb.Length < requiredLength || neg.Length < requiredLength)
            throw new ArgumentException("Projection scratch buffers must be at least Width * Height elements long.");
        var ri = r0.Inverse();
        double sqrtPiHalf = Math.Sqrt(Math.PI) / 2, twoOverPi = 2 / Math.PI;
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
            q0[i] = Math.Atan2(uy, ux) * twoOverPi;
        }
    }

    /// <summary>PrepareSpherePoint 済みの球点について面内回転角 phi のパターンを output へ書き込む (完全逐次、square 格子専用)。260725Cl 追加</summary>
    //260725Cl (/simplify) 変更: ①sector 判定を q0−pc の正規化角で除算なし化 (s∈[−6,2] なので (s+8)&3 = 整数一致の mod4)
    //②両半球あり (通常ケース) は半球 active チェックと InterpolatePlaneSquare のループ不変ガード・呼び出しを除去した特化ループ
    //(バイリニアは同一演算のインライン、step 除算は invStep 乗算化)。片半球欠けのみ従来形。数値は ULP 差 (等価群 — dict 回帰で検証)
    public void ProjectInPlane(MasterPattern mp, double phi, float[] posPlane, float[] negPlane, double[] q0, double[] ra, double[] rb, bool[] neg, double[] output)
    {
        int gs = ValidateProjectionBuffers(mp, posPlane, negPlane, output); //260725Ch
        ArgumentNullException.ThrowIfNull(q0); ArgumentNullException.ThrowIfNull(ra); //260725Ch
        ArgumentNullException.ThrowIfNull(rb); ArgumentNullException.ThrowIfNull(neg);
        if (q0.Length < output.Length || ra.Length < output.Length || rb.Length < output.Length || neg.Length < output.Length)
            throw new ArgumentException("Projection scratch buffers must be at least Width * Height elements long.");
        bool hasPos = posPlane is { Length: > 0 }, hasNeg = negPlane is { Length: > 0 };
        const double halfPi = Math.PI / 2;
        double pc = phi * (2 / Math.PI);
        double lim = MasterPattern.SquareLimit, invStep = gs / (2.0 * lim);
        if (hasPos && hasNeg) //通常ケース (実稼働合成は両半球あり)
        {
            for (int i = 0; i < output.Length; i++)
            {
                double g = q0[i] - pc;
                int s = (int)Math.Floor(g + 0.5); //支配軸 sector (0=+a,1=+b,2=−a,3=−b)、t = sector 内局所角 ∈ [−π/4, π/4)
                double t = (g - s) * halfPi;
                double a, b;
                switch ((s + 8) & 3)
                {
                    case 0: a = ra[i]; b = rb[i] * t; break;
                    case 1: b = ra[i]; a = -rb[i] * t; break;
                    case 2: a = -ra[i]; b = -rb[i] * t; break;
                    default: b = -ra[i]; a = rb[i] * t; break;
                }
                var plane = neg[i] ? negPlane : posPlane;
                double wf = (a + lim) * invStep - 0.5, hf = (lim - b) * invStep - 0.5;
                int w0 = (int)Math.Floor(wf), h0 = (int)Math.Floor(hf);
                double fw = wf - w0, fh = hf - h0;
                int w1 = Math.Clamp(w0 + 1, 0, gs - 1), h1 = Math.Clamp(h0 + 1, 0, gs - 1);
                w0 = Math.Clamp(w0, 0, gs - 1); h0 = Math.Clamp(h0, 0, gs - 1);
                output[i] = (float)((1 - fw) * (1 - fh) * plane[h0 * gs + w0] + fw * (1 - fh) * plane[h0 * gs + w1]
                                  + (1 - fw) * fh * plane[h1 * gs + w0] + fw * fh * plane[h1 * gs + w1]);
            }
        }
        else //片半球欠け (稀ケース): active チェック付きの従来形
        {
            for (int i = 0; i < output.Length; i++)
            {
                bool n = neg[i];
                if (n ? !hasNeg : !hasPos) { output[i] = 0; continue; }
                double g = q0[i] - pc;
                int s = (int)Math.Floor(g + 0.5);
                double t = (g - s) * halfPi;
                double a, b;
                switch ((s + 8) & 3)
                {
                    case 0: a = ra[i]; b = rb[i] * t; break;
                    case 1: b = ra[i]; a = -rb[i] * t; break;
                    case 2: a = -ra[i]; b = -rb[i] * t; break;
                    default: b = -ra[i]; a = rb[i] * t; break;
                }
                output[i] = MasterPattern.InterpolatePlaneSquare(n ? negPlane : posPlane, gs, a, b);
            }
        }
    }
    #endregion

    /// <summary>回転 rotation (crystal→sample) のパターンを output (Width×Height) へ書き込む。posPlane/negPlane = MasterPattern.GetPlane の単一スライス。
    /// parallel=false で行ループを逐次実行 (辞書総当たりのような方位単位で並列化する呼び出し向け。小ラスターでは行並列のオーバーヘッドが支配的)。260724Cl シグネチャ変更 (parallel 追加)</summary>
    //260724Cl 旧: public void Project(MasterPattern mp, Matrix3D rotation, float[] posPlane, float[] negPlane, double[] output)
    //260920Cl: 視線の計算をこのループへ融合する版も試したが、フル解像度の較正で 86.6 s → 93.5 s と遅くなったので採らない
    //  (raysSample 33 MB の往復を省くより、ホットループに分岐と PixelToSampleDirection が入る分が上回る)。視線は Rebuild で別に作る
    public void Project(MasterPattern mp, Matrix3D rotation, float[] posPlane, float[] negPlane, double[] output, bool parallel = true)
    {
        //int gs = mp.GridSize; //260725Ch 変更前
        int gs = ValidateProjectionBuffers(mp, posPlane, negPlane, output); //260725Ch
        var ri = rotation.Inverse();
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

                if (isHex) //六方格子 (EbsdPatternComposer.BuildLookupTable と同式)
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
        //double mean = 0; foreach (var x in src) mean += x; //260725Cl 変更前 (/simplify): SIMD 和 (SumSimd) へ
        double floor = Math.Max(1E-10, SumSimd(src) / src.Length * 0.05);
        //②log-ratio → dst。260725Cl: Vector.Log (net10 標準 SIMD) 化 — prof で log が前処理の 27% (Math.Log と数 ULP 差、
        //辞書ベンチの候補/refHit/misor 不変を確認して採用)
        //旧: for (int i = 0; i < dst.Length; i++) dst[i] = Math.Log(Math.Max(src[i], floor * 0.01) / Math.Max(tmp1[i], floor));
        int n = dst.Length, vc = Vector<double>.Count, i0 = 0;
        var vF2 = new Vector<double>(floor * 0.01);
        var vF = new Vector<double>(floor);
        for (; i0 <= n - vc; i0 += vc)
            Vector.Log(Vector.Max(new Vector<double>(src, i0), vF2) / Vector.Max(new Vector<double>(tmp1, i0), vF)).CopyTo(dst, i0);
        for (; i0 < n; i0++) dst[i0] = Math.Log(Math.Max(src[i0], floor * 0.01) / Math.Max(tmp1[i0], floor));
        //③DoG: g1(σ1.5) = dst→tmp1、g2(σ6) = dst→tmp2 (dst は両方の入力なので温存)
        Box3Seq(dst, tmp1, tmp2, w, h, 1.5);
        Box3Seq(dst, tmp2, tmp3, w, h, 6.0);
        //for (int i = 0; i < dst.Length; i++) dst[i] = tmp1[i] - tmp2[i]; //260725Cl 変更前 (/simplify): SIMD 差へ
        i0 = 0;
        for (; i0 <= n - vc; i0 += vc) (new Vector<double>(tmp1, i0) - new Vector<double>(tmp2, i0)).CopyTo(dst, i0);
        for (; i0 < n; i0++) dst[i0] = tmp1[i0] - tmp2[i0];
        //④標準化+tanh(z/3) の 1 パス融合 → 再標準化 (NormalizeInPlace は SIMD 統一済み)。
        //260725Cl (/simplify): 旧 NormalizeSimd(dst)→TanhOver3Simd(dst)→NormalizeSimd(dst) の中間 1 往復 (書いて読み直す全画素パス) を融合で削減
        NormalizeTanhSimd(dst);
        NormalizeInPlace(dst);
    }

    /// <summary>SIMD 水平和。260725Cl 追加 (/simplify: NormalizeInPlace/NormalizeTanhSimd/RobustPreprocessFast の 3 箇所で共用)</summary>
    static double SumSimd(double[] data)
    {
        int n = data.Length, vc = Vector<double>.Count, i = 0;
        var vs = Vector<double>.Zero;
        for (; i <= n - vc; i += vc) vs += new Vector<double>(data, i);
        double sum = Vector.Sum(vs);
        for (; i < n; i++) sum += data[i];
        return sum;
    }

    /// <summary>標準化 → tanh(x/3) の 1 パス融合 (in place)。tanh(x/3) = (e^{2x/3}−1)/(e^{2x/3}+1)、標準化後の値は ±60 に
    /// クランプ (tanh(20)=1−4E-18 で数値同一、Vector.Exp のオーバーフロー→NaN を防止)。260725Cl 追加</summary>
    static void NormalizeTanhSimd(double[] data)
    {
        int n = data.Length, vc = Vector<double>.Count, i = 0;
        double mean = SumSimd(data) / n;
        var vm = new Vector<double>(mean);
        var vv = Vector<double>.Zero;
        for (; i <= n - vc; i += vc) { var d = new Vector<double>(data, i) - vm; vv += d * d; }
        double var = Vector.Sum(vv);
        for (; i < n; i++) { double d = data[i] - mean; var += d * d; }
        double std = Math.Sqrt(var / n);
        if (std < 1E-12) std = 1;
        var vi = new Vector<double>(1 / std);
        var one = Vector<double>.One;
        var c = new Vector<double>(2.0 / 3.0);
        var lim = new Vector<double>(60.0);
        i = 0;
        for (; i <= n - vc; i += vc)
        {
            var x = Vector.Max(Vector.Min((new Vector<double>(data, i) - vm) * vi, lim), -lim);
            var e = Vector.Exp(x * c);
            ((e - one) / (e + one)).CopyTo(data, i);
        }
        for (; i < n; i++) data[i] = Math.Tanh((data[i] - mean) / std / 3);
    }

    [ThreadStatic] static double[] box3Scratch; //260724Cl: RobustPreprocessFast の第 4 バッファ (スレッドローカル再利用)

    /// <summary>260920Cl 追加: data から自身の広域背景 (box blur 3 連 = ガウシアン近似、半値全幅 fwhm) を差し引く。
    /// 表示側の <see cref="ImageProcess.SubtractGaussianBackground"/> と同じ高域通過だが、直接畳み込みでなく O(1)/画素 の box 近似なので
    /// 1 評価ごとに掛ける較正でも実用になる (fwhm 100 px のガウシアンは 600 演算/画素、box3 は 25 演算/画素)。work1/work2 は data と同じ長さ</summary>
    //260920Cl: parallel=true で box パスを分割して並列化する。較正の最終段 (開始点 1 つ・フル解像度) では投影は並列なのに
    //  背景差し引きだけ逐次だと、ここが Amdahl のボトルネックになる。辞書探索のように呼び出し側が既に並列なときは false (入れ子回避)
    internal static void SubtractBoxBackground(double[] data, double[] work1, double[] work2, int w, int h, double fwhm, bool parallel = false)
    {
        if (!(fwhm >= 1)) { Array.Clear(data); return; } //ぼかし幅が 1 px 未満なら 原画像 − 原画像 = 0
        Box3Seq(data, work1, work2, w, h, fwhm / 2.354820045, parallel);
        if (parallel) Parallel.For(0, h, y => { for (int i = y * w, e = i + w; i < e; i++) data[i] -= work1[i]; });
        else for (int i = 0; i < data.Length; i++) data[i] -= work1[i];
    }

    /// <summary>box blur 3 連 (ガウシアン分散一致近似)。src → dst (src は不変)。work は作業バッファ。260724Cl 追加。
    /// 260920Cl: parallel=false (既定) なら完全逐次 — 呼び出し側が既に並列な辞書探索用。true は行/列帯で並列化する (較正の最終段用)</summary>
    internal static void Box3Seq(double[] src, double[] dst, double[] work, int w, int h, double sigma, bool parallel = false) // 260920Cl: internal 化 + 並列フラグ
    {
        int r = Math.Max(1, (int)Math.Round((Math.Sqrt(4 * sigma * sigma + 1) - 1) / 2)); //3 連の合成分散 3(w²−1)/12 = σ² となる box 幅
        //260725Cl: 境界正規化 1/n を位置別に事前計算し BoxPassSeq 内の毎画素除算 (~18 回/画素) を乗算化 (prof: box が前処理の 52%)
        //260920Cl (/simplify2): parallel のときは後段の Parallel.For デリゲートへ配列で渡すので、stackalloc に作ってから
        //  ToArray で複製し直すのは二重手間 (w+h≈2400 要素 × 評価回数)。最初から配列に作る
        Span<double> invX = !parallel && w <= 2048 ? stackalloc double[w] : new double[w];
        Span<double> invY = !parallel && h <= 2048 ? stackalloc double[h] : new double[h];
        //260725Cl (/simplify): フル窓の内部は 1/(2r+1) 定数を共有 — 除算は両端 ~2r 個のみ (旧: 全 w+h 要素で除算。値は同一 = ビット一致)
        double invMid = 1.0 / (2 * r + 1);
        for (int x = 0; x < w; x++) invX[x] = x >= r && x < w - r ? invMid : 1.0 / (Math.Min(x + r, w - 1) - Math.Max(x - r, 0) + 1);
        for (int y = 0; y < h; y++) invY[y] = y >= r && y < h - r ? invMid : 1.0 / (Math.Min(y + r, h - 1) - Math.Max(y - r, 0) + 1);
        //260725Cl: 縦パスの行アキュムレータ化で横パス出力の一時バッファが必要 (in-place 縦パス廃止)
        if (parallel)
        {   //260920Cl (/simplify): tmpH は呼び出しスレッドで確保してワーカーへ渡すだけ (各ワーカーは重ならない範囲しか触らない) なので、
            //  逐次版と同じ ThreadStatic スクラッチが使える。旧実装はここで毎回 new double[w*h] = フル解像度で 1 評価 11 MB を捨てていた
            if (box3ScratchH == null || box3ScratchH.Length < w * h) box3ScratchH = new double[w * h];
            var tmpH = box3ScratchH;
            var ix = invX.ToArray(); var iy = invY.ToArray(); //Parallel.For のデリゲートへ Span は渡せないので配列へ (w+h 要素のみ)
            BoxPassPar(src, dst, tmpH, w, h, r, ix, iy);
            BoxPassPar(dst, work, tmpH, w, h, r, ix, iy);
            BoxPassPar(work, dst, tmpH, w, h, r, ix, iy);
            return;
        }
        if (box3ScratchH == null || box3ScratchH.Length < w * h) box3ScratchH = new double[w * h];
        BoxPassSeq(src, dst, box3ScratchH, w, h, r, invX, invY);
        BoxPassSeq(dst, work, box3ScratchH, w, h, r, invX, invY);
        BoxPassSeq(work, dst, box3ScratchH, w, h, r, invX, invY);
    }

    /// <summary>260920Cl 追加: BoxPassSeq の並列版。横パスは行ごと、縦パスは列帯ごとに分ける
    /// (縦パスは行アキュムレータを上から走らせるので y では割れないが、列で割れば各帯が独立)</summary>
    static void BoxPassPar(double[] src, double[] dst, double[] tmpH, int w, int h, int radius, double[] invX, double[] invY)
    {
        Parallel.For(0, h, y => //横パス: src → tmpH (行内 running sum)
        {
            int row = y * w;
            double sum = 0;
            for (int x = 0; x <= Math.Min(radius, w - 1); x++) sum += src[row + x];
            for (int x = 0; x < w; x++)
            {
                tmpH[row + x] = sum * invX[x];
                int add = x + radius + 1, rem = x - radius;
                if (add < w) sum += src[row + add];
                if (rem >= 0) sum -= src[row + rem];
            }
        });
        //260920Cl (/simplify2): 帯の境界を double 8 個 = 64 B (キャッシュライン) に揃える。揃えないと隣接する 2 帯が
        //  同じラインへ dst を書き、行ごとに無効化が往復する (w=1344・32 コアだと帯幅 42 で毎回ずれていた)。結果はビット一致
        int bands = Math.Min(Environment.ProcessorCount, Math.Max(1, w / 32));
        int bandWidth = ((w + bands - 1) / bands + 7) & ~7;
        Parallel.For(0, bands, b => //縦パス: tmpH → dst (列帯ごとに独立なアキュムレータ)
        {
            int x0 = b * bandWidth, x1 = Math.Min(w, x0 + bandWidth);
            if (x0 >= x1) return;
            var acc = new double[x1 - x0];
            for (int y = 0; y <= Math.Min(radius, h - 1); y++)
            {
                int row = y * w;
                for (int x = x0; x < x1; x++) acc[x - x0] += tmpH[row + x];
            }
            for (int y = 0; y < h; y++)
            {
                int row = y * w;
                int addRow = (y + radius + 1) < h ? (y + radius + 1) * w : -1, remRow = (y - radius) >= 0 ? (y - radius) * w : -1;
                double iv = invY[y];
                for (int x = x0; x < x1; x++)
                {
                    double a = acc[x - x0];
                    dst[row + x] = a * iv;
                    if (addRow >= 0) a += tmpH[addRow + x];
                    if (remRow >= 0) a -= tmpH[remRow + x];
                    acc[x - x0] = a;
                }
            }
        });
    }

    [ThreadStatic] static double[] box3ScratchH; //260725Cl: BoxPassSeq 横パス出力の一時 (スレッドローカル再利用)

    /// <summary>running box 平均 1 回分 (横+縦、境界は有効画素数で正規化、完全逐次)。src → dst (src 不変)。260724Cl 追加</summary>
    //260725Cl シグネチャ変更 (invX/invY/tmpH 追加): ①sum/n 除算 → sum·(1/n) 乗算 (テーブルは Box3Seq が事前計算)
    //②縦パスを列毎 running sum (stride アクセス) から行アキュムレータ+Vector<double> 融合 (行順アクセス) へ変更
    //旧: static void BoxPassSeq(double[] src, double[] dst, int w, int h, int radius)
    static void BoxPassSeq(double[] src, double[] dst, double[] tmpH, int w, int h, int radius, ReadOnlySpan<double> invX, ReadOnlySpan<double> invY)
    {
        for (int y = 0; y < h; y++) //横パス: src → tmpH (行内 running sum)
        {
            int row = y * w;
            double sum = 0;
            for (int x = 0; x <= Math.Min(radius, w - 1); x++) sum += src[row + x];
            for (int x = 0; x < w; x++)
            {
                tmpH[row + x] = sum * invX[x];
                int add = x + radius + 1, rem = x - radius;
                if (add < w) sum += src[row + add];
                if (rem >= 0) sum -= src[row + rem];
            }
        }
        //縦パス: tmpH → dst。行アキュムレータ acc[w] を上から走査 (dst 書き込み・acc 更新とも行順、SIMD 融合)
        Span<double> acc = w <= 2048 ? stackalloc double[w] : new double[w];
        acc.Clear();
        for (int y = 0; y <= Math.Min(radius, h - 1); y++)
        {
            int row = y * w;
            for (int x = 0; x < w; x++) acc[x] += tmpH[row + x];
        }
        int vc = Vector<double>.Count;
        for (int y = 0; y < h; y++)
        {
            int row = y * w;
            int addRow = (y + radius + 1) < h ? (y + radius + 1) * w : -1, remRow = (y - radius) >= 0 ? (y - radius) * w : -1;
            var vIv = new Vector<double>(invY[y]);
            int x = 0;
            for (; x <= w - vc; x += vc) //融合: dst 行 = acc·(1/n) → acc += 次行 → acc −= 抜け行
            {
                var a = new Vector<double>(acc.Slice(x, vc));
                (a * vIv).CopyTo(dst, row + x);
                if (addRow >= 0) a += new Vector<double>(tmpH, addRow + x);
                if (remRow >= 0) a -= new Vector<double>(tmpH, remRow + x);
                a.CopyTo(acc.Slice(x, vc));
            }
            for (; x < w; x++)
            {
                double a = acc[x];
                dst[row + x] = a * invY[y];
                if (addRow >= 0) a += tmpH[addRow + x];
                if (remRow >= 0) a -= tmpH[remRow + x];
                acc[x] = a;
            }
        }
    }

    /// <summary>zero-mean/unit-variance 化 (in place)</summary>
    //260725Cl 変更 (/simplify): スカラー実装 → SIMD 化し、旧 RobustPreprocessFast 専用 NormalizeSimd を本体へ統一 (同一カーネル 2 実装の解消)。
    //旧スカラー版とは加算順・逆数乗算の ULP 差 — dict/combo ベンチで候補・refHit・misor 不変を確認して採用。
    //旧: double mean = 0; foreach (var v in data) mean += v; mean /= data.Length;
    //    double var = 0; foreach (var v in data) { double d = v - mean; var += d * d; }
    //    double std = Math.Sqrt(var / data.Length); if (std < 1E-12) std = 1;
    //    for (int i = 0; i < data.Length; i++) data[i] = (data[i] - mean) / std;
    public static void NormalizeInPlace(double[] data)
    {
        int n = data.Length, vc = Vector<double>.Count, i = 0;
        double mean = SumSimd(data) / n;
        var vm = new Vector<double>(mean);
        var vv = Vector<double>.Zero;
        for (; i <= n - vc; i += vc) { var d = new Vector<double>(data, i) - vm; vv += d * d; }
        double var = Vector.Sum(vv);
        for (; i < n; i++) { double d = data[i] - mean; var += d * d; }
        double std = Math.Sqrt(var / n);
        if (std < 1E-12) std = 1;
        double inv = 1 / std;
        var vi = new Vector<double>(inv);
        i = 0;
        for (; i <= n - vc; i += vc) ((new Vector<double>(data, i) - vm) * vi).CopyTo(data, i);
        for (; i < n; i++) data[i] = (data[i] - mean) * inv;
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

    /// <summary>260920Cl 追加: <see cref="Zncc"/> の SIMD 版。フル解像度 (1.4 M 画素) を 1 評価ごとに 2 周する幾何較正のための高速版。
    /// 加算順が変わるため結果は ULP 程度ずれる。辞書探索の golden 順位を動かさないよう、既存の <see cref="Zncc"/> は残して
    /// 較正だけがこちらを使う (同じ数式・同じ規格化)</summary>
    //260920Cl: parallel=true で固定ブロック分割して並列化する (ブロック数は固定なので同じ入力なら常に同じ値 = 決定的)
    public static double ZnccSimd(double[] normalizedRef, double[] pattern, bool parallel = false)
    {
        int n = pattern.Length, vc = Vector<double>.Count, i = 0;
        if (parallel && n >= 1 << 16) return ZnccSimdParallel(normalizedRef, pattern);
        double mean = SumSimd(pattern) / n;
        var vm = new Vector<double>(mean);
        Vector<double> vvar = Vector<double>.Zero, vdot = Vector<double>.Zero;
        for (; i <= n - vc; i += vc)
        {
            var d = new Vector<double>(pattern, i) - vm;
            vvar += d * d;
            vdot += new Vector<double>(normalizedRef, i) * d;
        }
        double var = Vector.Sum(vvar), dot = Vector.Sum(vdot);
        for (; i < n; i++) { double d = pattern[i] - mean; var += d * d; dot += normalizedRef[i] * d; }
        double std = Math.Sqrt(var / n);
        return std < 1E-12 ? 0 : dot / (n * std);
    }

    /// <summary>260920Cl 追加: ZnccSimd の並列版。ブロック数を固定して部分和を決まった順に畳み込むので、同じ入力なら常に同じ値</summary>
    static double ZnccSimdParallel(double[] normalizedRef, double[] pattern)
    {
        int n = pattern.Length;
        //260920Cl (/simplify2): ブロック境界も 8 要素に揃える (未整列開始とスカラー末尾を減らす)。分割規則は固定なので決定的なまま
        int blocks = Math.Min(Environment.ProcessorCount, Math.Max(1, n / 8192));
        int size = ((n + blocks - 1) / blocks + 7) & ~7;
        var sums = new double[blocks];
        Parallel.For(0, blocks, b =>
        {
            int lo = b * size, hi = Math.Min(n, lo + size);
            //260920Cl (/simplify): ブロック本体もベクトル化する (旧実装はスカラーで、SIMD が効くのは逐次経路だけだった)
            int vc = Vector<double>.Count, k = lo;
            var vs = Vector<double>.Zero;
            for (; k <= hi - vc; k += vc) vs += new Vector<double>(pattern, k);
            double s = Vector.Sum(vs);
            for (; k < hi; k++) s += pattern[k];
            sums[b] = s;
        });
        double total = 0;
        foreach (var v in sums) total += v;
        double mean = total / n;
        var vars = new double[blocks]; var dots = new double[blocks];
        Parallel.For(0, blocks, b =>
        {
            int lo = b * size, hi = Math.Min(n, lo + size);
            int vc = Vector<double>.Count, k = lo;
            var vm = new Vector<double>(mean);
            Vector<double> vv = Vector<double>.Zero, vd = Vector<double>.Zero;
            for (; k <= hi - vc; k += vc)
            {
                var d = new Vector<double>(pattern, k) - vm;
                vv += d * d;
                vd += new Vector<double>(normalizedRef, k) * d;
            }
            double va = Vector.Sum(vv), dt = Vector.Sum(vd);
            for (; k < hi; k++) { double d = pattern[k] - mean; va += d * d; dt += normalizedRef[k] * d; }
            vars[b] = va; dots[b] = dt;
        });
        double var2 = 0, dot2 = 0;
        for (int b = 0; b < blocks; b++) { var2 += vars[b]; dot2 += dots[b]; }
        double std = Math.Sqrt(var2 / n);
        return std < 1E-12 ? 0 : dot2 / (n * std);
    }

    /// <summary>
    /// 簡易 Nelder-Mead (初期ステップ明示・下降単体法)。objective を最小化する。260724Cl 追加
    /// MathNet の NelderMeadSimplex は初期シンプレックス制御が弱いため自前実装 (数変数・数百評価の用途限定)。
    /// </summary>
    /// <summary>260920Cl 追加: <see cref="NelderMead"/> を最良点から刻みを縮めて繰り返す (再起動付き)。素の Nelder-Mead は谷で停滞するので、
    /// 収まった点を中心に初期シンプレックスを張り直すと残差がさらに下がる。改善が tol 未満になるか restarts 回で打ち切る。
    /// ⚠260920Cl (/simplify2): 直上の「簡易 Nelder-Mead…」の doc は下の <see cref="NelderMead"/> のものだったが、
    ///   本メソッドを挿し込んだときに帰属がずれていた。下の本体側へ doc を戻してある</summary>
    /// <param name="restarts">刻みを縮めて張り直す最大回数</param>
    /// <param name="shrink">1 回の再起動で刻みに掛ける係数</param>
    /// <param name="xtolRel">初期刻みに対するシンプレックスの広がりの打ち切り閾値 (0 = 関数値の幅だけで打ち切る)</param>
    public static (double[] Best, double Value, int Evaluations) NelderMeadRestart(Func<double[], double> objective, double[] start, double[] step,
        int maxEvalPerRun = 400, double tol = 1E-5, double xtolRel = 0, int restarts = 3, double shrink = 0.25)
    {
        var best = (double[])start.Clone();
        var scale = (double[])step.Clone();
        double bestValue = double.PositiveInfinity;
        int evalTotal = 0;
        for (int k = 0; k <= restarts; k++)
        {
            var (b, v, e) = NelderMead(objective, best, scale, maxEvalPerRun, tol, xtolRel);
            evalTotal += e;
            //260920Cl (/simplify2): 1 回目は無条件に採用する。旧実装は v が NaN だと 2 つの比較が両方 false になり、
            //  bestValue が +∞ のまま返って呼び出し側の順位付けとステータス表示が壊れていた
            if (k == 0) { bestValue = v; best = b; if (double.IsNaN(v)) break; continue; }
            if (!(v < bestValue - tol)) { if (v < bestValue) { bestValue = v; best = b; } break; } //改善が止まったら終了
            bestValue = v; best = b;
            for (int i = 0; i < scale.Length; i++) scale[i] *= shrink;
        }
        return (best, bestValue, evalTotal);
    }

    /// <summary>簡易 Nelder-Mead (初期ステップ明示・下降単体法)。objective を最小化する。260724Cl 追加。
    /// MathNet の NelderMeadSimplex は初期シンプレックス制御が弱いため自前実装 (数変数・数百評価の用途限定)</summary>
    /// <param name="xtolRel">初期刻みに対するシンプレックスの広がりの打ち切り閾値。0 (既定) なら関数値の幅だけで打ち切る = 従来動作</param>
    //260920Cl シグネチャ変更: xtolRel を追加。既存の呼び出し側 (方位探索) は既定値のままなので挙動不変
    //旧: public static (double[] Best, double Value, int Evaluations) NelderMead(Func<double[], double> objective, double[] start, double[] step, int maxEval = 400, double tol = 1E-5)
    public static (double[] Best, double Value, int Evaluations) NelderMead(Func<double[], double> objective, double[] start, double[] step, int maxEval = 400, double tol = 1E-5, double xtolRel = 0)
    {
        int n = start.Length;
        var simplex = new double[n + 1][];
        var values = new double[n + 1];
        int eval = 0;
        for (int i = 0; i <= n; i++)
        {
            simplex[i] = (double[])start.Clone();
            if (i > 0) simplex[i][i - 1] += step[i - 1];
            //values[i] = objective(simplex[i]); eval++; //260725Ch 変更前: 初期 simplex だけ NaN を +∞ へ正規化していなかった
            double initialValue = objective(simplex[i]); eval++; //260725Ch
            values[i] = double.IsNaN(initialValue) ? double.PositiveInfinity : initialValue;
        }
        while (eval < maxEval)
        {
            Array.Sort(values, simplex);
            // if (Math.Abs(values[n] - values[0]) < tol) break; // 260920Cl 変更前: 関数値の幅だけで打ち切っていた。
            //   PC-DD が縮退した平らな谷では ZNCC が tol 以下しか動かなくてもパラメータは大きく動けるので、
            //   xtolRel > 0 のときはシンプレックスの広がり (初期刻み比) も同時に小さいことを要求する
            if (Math.Abs(values[n] - values[0]) < tol)
            {
                if (xtolRel <= 0) break;
                double extent = 0;
                for (int i = 1; i <= n; i++)
                    for (int j = 0; j < n; j++)
                        extent = Math.Max(extent, Math.Abs(simplex[i][j] - simplex[0][j]) / Math.Max(1E-300, Math.Abs(step[j])));
                if (extent < xtolRel) break;
            }

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
