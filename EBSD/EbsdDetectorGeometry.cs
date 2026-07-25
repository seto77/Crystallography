#region using
using System;
using V3 = OpenTK.Mathematics.Vector3d;
#endregion

namespace Crystallography;

/// <summary>
/// EBSD 検出器幾何の不変スナップショット。260724Cl 追加 (実測パターン指数付け・幾何較正の共有幾何モデル)。
/// 規約は FormEBSD (フォワードモデル) と厳密に一致させる:
///   検出器中心 C=(DetX,-DetY,-DetZ)、面内基底 ex=(1,0,0)・ey=(0,-cosδ,-sinδ)、法線 n=(0,sinδ,-cosδ) (δ=DetTilt)。
///   画像ピクセル (col,row) はコーナー原点・ピクセル中心 (col+0.5,row+0.5)。画像中心からの表示 mm 座標 (u,v) に対し
///   検出器面上の点 P = C + xm·u·ex + v·ey (xm=左右反転 ±1)。視線 (lab) = -P/|P|、試料系 = Rx(SmpTilt)·(lab)。
///   (FormEBSD.BuildEbsdLookupTable の labDir と数式一致することを EbsdIndexCheck ハーネスで検証)
/// </summary>
public sealed class EbsdDetectorGeometry
{
    /// <summary>検出器傾斜 δ (radian)。既定 90°=π/2</summary>
    public double DetTilt { get; }
    /// <summary>検出器中心座標パラメータ (numericBox 値。lab 中心は (DetX,-DetY,-DetZ))</summary>
    public double DetX { get; }
    public double DetY { get; }
    public double DetZ { get; }
    /// <summary>ピクセルサイズ (mm/px)</summary>
    public double PixelSize { get; }
    /// <summary>検出器のピクセル数</summary>
    public int WidthPx { get; }
    public int HeightPx { get; }
    /// <summary>左右反転 (+1 or -1)。FormEBSD.DetectorXMirror と同じ</summary>
    public double XMirror { get; }
    /// <summary>試料傾斜 (radian)</summary>
    public double SampleTilt { get; }

    //派生量 (lab 系)
    public V3 Center { get; }
    public V3 Ex { get; }
    public V3 Ey { get; }
    public V3 Normal { get; }
    /// <summary>照射点から検出器面への垂直距離 |n・C| (= FormEBSD.CameraLength2)</summary>
    public double CameraLength { get; }

    readonly double sinSmp, cosSmp;

    public EbsdDetectorGeometry(double detTilt, double detX, double detY, double detZ, double pixelSize, int widthPx, int heightPx, double xMirror, double sampleTilt)
    {
        //260725Ch: 公開スナップショットの除算・正規化前提を入口で保証し、後段のNaN連鎖を明瞭な引数例外にする
        if (!double.IsFinite(detTilt)) throw new ArgumentOutOfRangeException(nameof(detTilt));
        if (!double.IsFinite(detX)) throw new ArgumentOutOfRangeException(nameof(detX));
        if (!double.IsFinite(detY)) throw new ArgumentOutOfRangeException(nameof(detY));
        if (!double.IsFinite(detZ)) throw new ArgumentOutOfRangeException(nameof(detZ));
        if (!(pixelSize > 0) || !double.IsFinite(pixelSize)) throw new ArgumentOutOfRangeException(nameof(pixelSize));
        if (widthPx <= 0) throw new ArgumentOutOfRangeException(nameof(widthPx));
        if (heightPx <= 0) throw new ArgumentOutOfRangeException(nameof(heightPx));
        if (xMirror is not (1.0 or -1.0)) throw new ArgumentOutOfRangeException(nameof(xMirror), "xMirror must be +1 or -1.");
        if (!double.IsFinite(sampleTilt)) throw new ArgumentOutOfRangeException(nameof(sampleTilt));

        DetTilt = detTilt; DetX = detX; DetY = detY; DetZ = detZ;
        PixelSize = pixelSize; WidthPx = widthPx; HeightPx = heightPx;
        XMirror = xMirror; SampleTilt = sampleTilt;

        var (sinDet, cosDet) = Math.SinCos(detTilt);
        Center = new V3(detX, -detY, -detZ);
        Ex = new V3(1, 0, 0);
        Ey = new V3(0, -cosDet, -sinDet);
        Normal = new V3(0, sinDet, -cosDet);
        CameraLength = Math.Abs(V3.Dot(Normal, Center));
        if (CameraLength < 1E-12) throw new ArgumentException("The detector plane must not pass through the sample origin."); //260725Ch: LineToLabNormal等の0除算を防ぐ
        (sinSmp, cosSmp) = Math.SinCos(sampleTilt);
    }

    /// <summary>ピクセル (col,row) の中心 (col+0.5, row+0.5) → 画像中心基準の表示 mm 座標 (u,v)</summary>
    public (double u, double v) PixelToMm(double col, double row)
        => ((col + 0.5 - WidthPx / 2.0) * PixelSize, (row + 0.5 - HeightPx / 2.0) * PixelSize);

    /// <summary>ピクセル中心 → 検出器面上の点 P (lab)</summary>
    public V3 PixelToLabPoint(double col, double row)
    {
        var (u, v) = PixelToMm(col, row);
        return Center + XMirror * u * Ex + v * Ey;
    }

    /// <summary>lab ベクトル → 試料系 (Rx(SampleTilt): y'=c·y+s·z, z'=-s·y+c·z。OpenTK CreateRotationX の M・v 作用と同一)</summary>
    public V3 LabToSample(in V3 lab)
        => new(lab.X, cosSmp * lab.Y + sinSmp * lab.Z, -sinSmp * lab.Y + cosSmp * lab.Z);

    /// <summary>ピクセル中心 → 視線方向 (試料系、正規化)。FormEBSD の labDir=-P 規約と同一</summary>
    public V3 PixelToSampleDirection(double col, double row)
    {
        var p = PixelToLabPoint(col, row);
        return LabToSample(-p.Normalized());
    }

    /// <summary>
    /// 画像ピクセル座標系の直線 A·col + B·row + C = 0 (ピクセル中心規約) から、
    /// そのバンド中心線と試料原点を含む面の法線 (lab、正規化、符号不定) を返す。
    /// </summary>
    public V3 LineToLabNormal(double lineA, double lineB, double lineC)
    {
        //ピクセル → 表示 mm へ係数変換: col = u/s + cx, row = v/s + cy (s=PixelSize, cx=W/2-0.5, cy=H/2-0.5)
        double cx = WidthPx / 2.0 - 0.5, cy = HeightPx / 2.0 - 0.5;
        double a = lineA / PixelSize, b = lineB / PixelSize;
        double c = lineA * cx + lineB * cy + lineC;

        //g ∝ a·eu + b·ey + γ·n、eu = xm·ex、γ = (c - a·(eu・C) - b·(ey・C)) / (n・C)
        var eu = XMirror * Ex;
        double nDotC = V3.Dot(Normal, Center);
        double gamma = (c - a * V3.Dot(eu, Center) - b * V3.Dot(Ey, Center)) / nDotC;
        return (a * eu + b * Ey + gamma * Normal).Normalized();
    }

    //260724Cl (/simplify): LineToSampleNormal (LabToSample∘LineToLabNormal の合成) は未使用のため削除 (利用側は gLab を挟んで個別に呼ぶ)

    /// <summary>バンド縁の 1 点 (ピクセル) から sinθB = |ĝ・P_e|/|P_e| を返す (ĝ=バンド中心面法線 lab)</summary>
    public double SinBraggFromEdgePoint(in V3 gLabNormalized, double col, double row)
    {
        var p = PixelToLabPoint(col, row);
        return Math.Abs(V3.Dot(gLabNormalized, p)) / p.Length;
    }

    /// <summary>試料系ベクトル → lab (LabToSample の逆変換)</summary>
    public V3 SampleToLab(in V3 sample)
        => new(sample.X, cosSmp * sample.Y - sinSmp * sample.Z, sinSmp * sample.Y + cosSmp * sample.Z);

    /// <summary>面法線 (lab) → 画像ピクセル座標系のバンド中心線係数 A·col + B·row + C = 0 (LineToLabNormal の逆投影)</summary>
    public (double A, double B, double C) LabNormalToLine(in V3 gLab)
    {
        double aMm = XMirror * gLab.X, bMm = V3.Dot(gLab, Ey), cMm = V3.Dot(gLab, Center);
        double cx = WidthPx / 2.0 - 0.5, cy = HeightPx / 2.0 - 0.5;
        return (aMm * PixelSize, bMm * PixelSize, cMm - aMm * PixelSize * cx - bMm * PixelSize * cy);
    }

    /// <summary>PC (垂線の足) の物理面内 mm 座標 (検出器中心基準)。260724Cl (/simplify) 追加: GetPatternCenter・幾何較正に重複していた式を一元化</summary>
    public (double footU, double footV) PatternCenterMm => (-DetX, -(DetY * Math.Cos(DetTilt) + DetZ * Math.Sin(DetTilt)));

    /// <summary>パターンセンター (PC、垂線の足) の画像ピクセル座標と検出器距離 DD (mm) を返す。
    /// PC の物理面内座標 (検出器中心基準) は (-DetX, -(DetY cosδ + DetZ sinδ))、表示側では X に xm が掛かる。</summary>
    public (double pcCol, double pcRow, double dd) GetPatternCenter()
    {
        var (footU, footV) = PatternCenterMm; //物理面内 (検出器中心基準)
        double uView = XMirror * footU; //表示 mm
        return (uView / PixelSize + WidthPx / 2.0 - 0.5, footV / PixelSize + HeightPx / 2.0 - 0.5, CameraLength);
    }

    /// <summary>PC (物理面内 mm、検出器中心基準) + DD + δ から DetX/DetY/DetZ を逆算 (標準配置 signed L = -DD を仮定)</summary>
    public static (double detX, double detY, double detZ) FromPatternCenter(double footU, double footV, double dd, double detTilt)
    {
        var (sinDet, cosDet) = Math.SinCos(detTilt);
        // DetY·cosδ + DetZ·sinδ = -footV, DetY·sinδ - DetZ·cosδ = -DD
        double detY = -(footV * cosDet + dd * sinDet);
        double detZ = -footV * sinDet + dd * cosDet;
        return (-footU, detY, detZ);
    }
}
