#region using 
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Drawing.Imaging;
using System.Linq;
using System.Numerics;
using System.Runtime.InteropServices;
using System.Threading.Tasks;
#endregion

namespace Crystallography;

[Serializable]
public class PseudoBitmap : IDisposable
{
    public void Dispose()
    {
#pragma warning disable CS0618 // 261111Cl 追加: SrcBitmap は [Obsolete] (4 アプリで未使用)。削除までは従来どおり解放する
        SrcBitmap?.Dispose();
#pragma warning restore CS0618
        destBmp?.Dispose();
    }

    #region static なコンストラクタ、フィールド

    public static (byte R, byte G, byte B)[] ColorScaleGrayLiner = new (byte R, byte G, byte B)[65536];
    public static (byte R, byte G, byte B)[] ColorScaleGrayLog = new (byte R, byte G, byte B)[65536];


    public static (byte R, byte G, byte B)[] ColorScaleColdWarmLiner = new (byte R, byte G, byte B)[65536];
    public static (byte R, byte G, byte B)[] ColorScaleColdWarmLog = new (byte R, byte G, byte B)[65536];

    public static (byte R, byte G, byte B)[] ColorScaleFireLiner = new (byte R, byte G, byte B)[65536];
    public static (byte R, byte G, byte B)[] ColorScaleFireLog = new (byte R, byte G, byte B)[65536];

    public static (byte R, byte G, byte B)[] ColorScaleSpectrumLiner = new (byte R, byte G, byte B)[65536];
    public static (byte R, byte G, byte B)[] ColorScaleSpectrumLog = new (byte R, byte G, byte B)[65536];

    public static (byte R, byte G, byte B)[] ColorScaleRotationLiner = new (byte R, byte G, byte B)[65536];
    public static (byte R, byte G, byte B)[] ColorScaleRotationLog = new (byte R, byte G, byte B)[65536];

    static PseudoBitmap()
    {
        #region スケールの初期化

        var gray = new (int I, byte R, byte G, byte B)[]{
            (0,0,0,0),
            (65536,255,255,255)
        };

        var cold_warm = new (int I, byte R, byte G, byte B)[]{
            (0,0,0,0),
            (5376,0,0,255),
            (16256,0,191,191),
            (27136,0,255,0),
            (38016,255,255,0),
            (48896,255,0,0),
            (59904,255,0,255),
            (65536,255,255,255)
        };

        var fire = new (int I, byte R, byte G, byte B)[]{
            (0,0,0,0),
            (3840, 0,0,91),
            (7936,46,0,160),
            (12032,94,0,216),
            (14336,122,0,227),
            (16128,143,0, 212),
            (20224,171,0,154),
            (24320,193,0,96),
            (24576,195,0,93),
            (24832,196,1,89),
            (27392,210,21,53),
            (32512,238,76,0),
            (36608,254,115,0),
            (36864,255,117,0),
            (40704,255, 145,0),
            (45312,255,176, 0 ),
            (48896,255,203, 0 ),
            (53248,255,234, 0 ),
            (57088,255,254,90 ),
            (57344,255,255,98 ),
            (61184,255,255, 215 ),
            (65536,255,255,255)
        };

        var spectrum = new (int I, byte R, byte G, byte B)[]{
            (0,255,0,0),
            (10752,255,252,0),
            (11008,252,255,0),
            (21760,0,255,0),
            (32512,0,255,252),
            (32768,0,252,255),
            (43520,0,0,255),
            (54272,252,0,255),
            (54528,255,0,252),
            (65536,255, 0, 0),
        };

        var rotation = new (int I, byte R, byte G, byte B)[]{
            (0,     80, 80,    255),       //0
            (10922, 80,    255,    255),           //60°
            (21845, 80,    255,    80),            //120°
            (32767, 255, 255,    80),     //180°
            (43690, 255, 80,    80),         //240°
            (54613, 255, 80,    255),         //300°
            (65536, 80, 80,    255),               //360°

        };

        (byte R, byte G, byte B)[] setScaleLinear((int I, byte R, byte G, byte B)[] src)
        {
            var result = new (byte R, byte G, byte B)[65536];
            for (int j = 0; j < src.Length - 1; j++)
            {
                (int I, byte R, byte G, byte B) a = src[j], b = src[j + 1];
                for (int k = src[j].I; k < src[j + 1].I; k++)
                {

                    var x = (double)(b.I - k) / (b.I - a.I);
                    result[k] = ((byte)(a.R * x + b.R * (1 - x)), (byte)(a.G * x + b.G * (1 - x)), (byte)(a.B * x + b.B * (1 - x)));
                }
            }
            return result;
        }

        (byte R, byte G, byte B)[] setScaleLog((int I, byte R, byte G, byte B)[] src)
        {
            var result = new (byte R, byte G, byte B)[65536];
            for (int i = 0; i < 65536; i++)
            {
                var k = (int)(Math.Log(i, 65536) * 65536);
                for (int j = 0; j < src.Length - 1; j++)
                    if (src[j].I < k && src[j + 1].I >= k)
                    {
                        (int I, byte R, byte G, byte B) a = src[j], b = src[j + 1];
                        double x = (double)(b.I - k) / (b.I - a.I);
                        result[i] = ((byte)(a.R * x + b.R * (1 - x)), (byte)(a.G * x + b.G * (1 - x)), (byte)(a.B * x + b.B * (1 - x)));
                        break;
                    }
            }
            return result;
        }

        ColorScaleGrayLiner = setScaleLinear(gray);
        ColorScaleColdWarmLiner = setScaleLinear(cold_warm);
        ColorScaleFireLiner = setScaleLinear(fire);
        ColorScaleRotationLiner = setScaleLinear(rotation);
        ColorScaleSpectrumLiner = setScaleLinear(spectrum);

        ColorScaleGrayLog = setScaleLog(gray);
        ColorScaleColdWarmLog = setScaleLog(cold_warm);
        ColorScaleFireLog = setScaleLog(fire);
        ColorScaleRotationLog = setScaleLog(rotation);
        ColorScaleSpectrumLog = setScaleLog(spectrum);

        #endregion スケールの初期化
    }
    #endregion 

    #region コンストラクタ

    public PseudoBitmap()
    {
        ColorScale = ColorScaleGrayLiner;
    }

    /// <summary>List<uint>型がソースデータであるグレー実画像のコンストラクタ</summary>
    /// <param name="valueGray"></param>
    /// <param name="bitsPerPixels"></param>
    /// <param name="width"></param>
    /// <param name="scaleR"></param>
    /// <param name="scaleG"></param>
    /// <param name="scaleB"></param>
    /// <param name="normarize"></param>
    /// <param name="realImage"></param>
    public PseudoBitmap(List<uint> valueGray, int width, byte[] scaleR = null, byte[] scaleG = null, byte[] scaleB = null, bool realImage = true)
        : this(valueGray.Select(a => (double)a).ToArray(), width, scaleR, scaleG, scaleB, realImage)
    { }

    /// <summary>List<uint>型がソースデータであるグレー実画像のコンストラクタ</summary>
    /// <param name="valueGray"></param>
    /// <param name="scale"></param>
    /// <param name="width"></param>
    /// <param name="scaleR"></param>
    /// <param name="scaleG"></param>
    /// <param name="scaleB"></param>
    /// <param name="realImage"></param>
    public PseudoBitmap(List<uint> valueGray, double scale, int width, byte[] scaleR = null, byte[] scaleG = null, byte[] scaleB = null, bool realImage = true)
      : this(valueGray.Select(a => a * scale).ToArray(), width, scaleR, scaleG, scaleB, realImage)
    {
    }

    /// <summary>グレーの実画像を描画する.</summary>
    /// <param name="valueR"></param>
    /// <param name="valueG"></param>
    /// <param name="valueB"></param>
    /// <param name="width"></param>
    public PseudoBitmap(uint[] value, int width, byte[] scaleR = null, byte[] scaleG = null, byte[] scaleB = null, bool realImage = true)
        : this(value.Select(a => (double)a).ToArray(), width, scaleR, scaleG, scaleB, realImage)
    { }

    /// <summary>double型がソースデータであるグレー実画像のコンストラクタ</summary>
    /// <param name="values"></param>
    /// <param name="width"></param>
    /// <param name="scaleR"></param>
    /// <param name="scaleG"></param>
    /// <param name="scaleB"></param>
    /// <param name="realImage"></param>
    public PseudoBitmap(double[] values, int width, byte[] scaleR, byte[] scaleG, byte[] scaleB, bool realImage = true)
    {
        if (scaleR == null || scaleG == null || scaleB == null || scaleR.Length != scaleG.Length || scaleG.Length != scaleB.Length || scaleR.Length < 2)
            constructor(values, width, null, realImage);
        else
        {
            var scale = new (byte R, byte G, byte B)[scaleR.Length];
            for (int i = 0; i < scale.Length; i++)
                scale[i] = (scaleR[i], scaleG[i], scaleB[i]);
            constructor(values, width, scale, realImage);
        }
    }

    /// <summary>double型がソースデータであるグレー実画像のコンストラクタ</summary>
    /// <param name="values"></param>
    /// <param name="width"></param>
    /// <param name="scaleR"></param>
    /// <param name="scaleG"></param>
    /// <param name="scaleB"></param>
    /// <param name="realImage"></param>
    public PseudoBitmap(double[] values, int width, (byte R, byte G, byte B)[] scale = null, bool realImage = true)
    {
        constructor(values, width, scale, realImage);
    }

    private void constructor(double[] values, int width, (byte R, byte G, byte B)[] scale = null, bool realImage = true)
    {
        if (values.Length == 0 || values.Length % width != 0) throw new ArgumentNullException("Invalid input");
        Width = width;
        Height = values.Length / width;

        GrayScale = true;
        RealImage = true;

        SrcValuesGrayOriginal = SrcValuesGray = values;

        IsSrcGray = true;

#pragma warning disable CS0618 // 261111Cl 追加: SrcValuesR/G/B は [Obsolete] (4 アプリで未使用)
        SrcValuesR.Clear();
        SrcValuesG = SrcValuesB = SrcValuesR;
#pragma warning restore CS0618

        if (scale == null || scale.Length < 2)
            ColorScale = ColorScaleGrayLiner;
        else
        {
            ColorScale = scale;
        }


        initFilter();
    }


    #endregion コンストラクタ

    #region 初期化関連

    [Obsolete("261111Cl: 4 アプリで未使用。将来削除予定")]
    private static Complex[][] initComplex(int width, int height)
    {
        Complex[][] complex = new Complex[height][];
        for (int i = 0; i < height; i++)
        {
            complex[i] = new Complex[width];
            for (int j = 0; j < width; j++)
                complex[i][j] = new Complex(0, 0);
        }
        return complex;
    }

    public void initFilter()
    {
        Filter1.Clear();
        Filter1.AddRange(new bool[Height * Width]);

        Filter2.Clear();
        Filter2.AddRange(new bool[Height * Width]);

        Filter3.Clear();
        Filter3.AddRange(new bool[Height * Width]);

        Filter4.Clear();
        Filter4.AddRange(new bool[Height * Width]);

        Filter5.Clear();
        Filter5.AddRange(new bool[Height * Width]);

        FilterTemporary.Clear();
        FilterTemporary.AddRange(new bool[Height * Width]);

#pragma warning disable CS0618 // 261111Cl 追加: FFT_Filter は [Obsolete] (4 アプリで未使用)
        if (Height * Width < 3001 * 3001)
            FFT_Filter = Enumerable.Repeat(0f, Height * Width).ToList();
#pragma warning restore CS0618

    }

    #endregion 初期化関連

    #region プロパティ

    /// <summary>グレイスケールモードの時、SrcValuesGrayOriginalをSrcValuesGrayとは異なる配列として確保する。newの直後に設定しないといけない。</summary>
    public bool ReserveSrcValuesGrayOriginal
    {
        set
        {
            if (value)
            {
                SrcValuesGrayOriginal = new double[SrcValuesGray.Length];
                Array.Copy(SrcValuesGray, SrcValuesGrayOriginal, SrcValuesGray.Length);
            }
            else
            {
                SrcValuesGrayOriginal = SrcValuesGray;
            }
        }
    }

    /// <summary>アルファチャンネルが有効かどうか. FilterAlphaもセットしないと、機能しない.</summary>
    public bool AlphaEnabled { set; get; } = false;

    /// <summary>上下方向の反転をするかどうか</summary>
    public bool VerticalFlip { set; get; } = false;

    /// <summary>左右方向の反転をするかどうか</summary>
    public bool HorizontalFlip { set; get; } = false;

    public enum Scales { GrayLinear, ColdWarmLinear, RotationLinear, SpectrumLinear, FireLinear, GrayLog, ColdWarmLog, RotationLog, SpectrumLog, FireLog, }

    private Scales scale = Scales.GrayLinear;
    public Scales Scale
    {
        get => scale;
        set
        {
            scale = value;
            if (scale == Scales.GrayLinear)
                SetScaleGray(true);
            else if (scale == Scales.GrayLog)
                SetScaleGray(false);
            else if (scale == Scales.ColdWarmLinear)
                SetScaleColdWarm(true);
            else if (scale == Scales.ColdWarmLog)
                SetScaleColdWarm(false);
            else if (scale == Scales.RotationLinear)
                SetScaleRotation();
            else if (scale == Scales.SpectrumLinear)
                SetScaleSpectrum(true);
            else if (scale == Scales.SpectrumLog)
                SetScaleSpectrum(false);
            else if (scale == Scales.FireLinear)
                SetScaleFire(true);
            else if (scale == Scales.FireLog)
                SetScaleFire(false);
        }
    }




    public bool RealImage = true;

    public PointF Center;//虚画像のときの中心ピクセル位置

    /// <summary>Imageの幅</summary>
    public int Width { get; set; }

    public int Height { get; set; }//高さ

    [Obsolete("261111Cl: 4 アプリで未使用。将来削除予定")]
    public Complex[][] ComplexGray;//元画像の複素数;
    [Obsolete("261111Cl: 4 アプリで未使用。将来削除予定")]
    public Complex[][] ComplexR;//元画像Rの複素数;
    [Obsolete("261111Cl: 4 アプリで未使用。将来削除予定")]
    public Complex[][] ComplexG;//元画像Gの複素数;
    [Obsolete("261111Cl: 4 アプリで未使用。将来削除予定")]
    public Complex[][] ComplexB;//元画像Bの複素数;

    [Obsolete("261111Cl: 4 アプリで未使用。将来削除予定")]
    public List<float> FFT_Filter = [];
    public List<bool> FilterTemporary = [];
    public List<bool> Filter1 = [], Filter2 = [], Filter3 = [], Filter4 = [], Filter5 = [];
    public bool Filter1Visible = true, Filter2Visible = true, Filter3Visible = true, Filter4Visible = true, Filter5Visible = true;

    /// <summary>透明度フィルター. 0が透明、255が不透明. AlphaEnabledがTrueの時だけ使われる.</summary>
    public List<byte> FilterAlfha { get; set; } = [];

    public (byte R, byte G, byte B)[] ColorScale;

    [Obsolete("261111Cl: 4 アプリで未使用。将来削除予定")]
    public List<uint> SrcValuesR = [], SrcValuesG = [], SrcValuesB = [];

    /// <summary>グレースケールの画像データ (BlurModeが有効の時は、フィルター後のデータ)</summary>
    public double[] SrcValuesGray;

    /// <summary>グレースケールのもとの画像データ</summary>
    public double[] SrcValuesGrayOriginal;

    /// <summary>Blurモードの列挙体</summary>
    public enum BlurModeEnum { Gaussian, Lorentzian, None }

    private enum ImageSouceEnum { GrayDouble, GrayUint, ColorDouble, ColorUint }

    /// <summary>グレイスケールとして表示するかどうか　(ソースがグレースケールかどうかは無関係)</summary>
    public bool GrayScale;

    /// <summary>ソースがグレースケールかどうか (グレースケールとして表示するかどうかは無関係)</summary>
    public bool IsSrcGray = false;

    public UniversalConstants.LengthUnit Unit = UniversalConstants.LengthUnit.mm;
    public double PixelSizeX = 1;//単位はunit
    public double PixelSizeY = 1;

    /// <summary>表示する強度の上限 (画像中の最大強度という意味ではない）</summary>
    public double MaxValue = 255;

    /// <summary>表示する強度の下限 (画像中の最大強度という意味ではない）</summary>
    public double MinValue = 0;

    //public int BitsPerPixels = 8;

    public bool IsNegative = false;

    [Obsolete("261111Cl: 4 アプリで未使用。将来削除予定")]
    public Bitmap SrcBitmap;//実画像のときの元画像

    /// <summary>画像に関連付けた情報 (なんでも)</summary>
    public object Tag { get; set; }

    public Profile FrequencyProfile { get; set; }

    /// <summary>1ピクセルが相当する数値 (単位は PixelUnitで指定)</summary>
    public double PixelScale;

    /// <summary>ピクセルの単位</summary>
    public LengthUnitEnum PixelUnit;

    #endregion プロパティ

    #region 操作

    /// <summary>オリジナルのピクセル強度データにリセットする。Grayスケールモードの時だけ有効。</summary>
    public void SetOriginalGray()
    {
        if (IsSrcGray)
            SrcValuesGray = SrcValuesGrayOriginal;
    }

    /// <summary>ダスト＆スクラッチを施す。Grayスケールモードの時だけ有効。</summary>
    /// <param name="radius"></param>
    /// <param name="threshold"></param>
    /// <param name="applyToOriginal"></param>
    public void SetDustAndScratches(double radius, double threshold, bool applyToOriginal = true)
    {
        if (IsSrcGray)
        {
            if (applyToOriginal)
                SrcValuesGray = ImageProcess.DustAndScratches(SrcValuesGrayOriginal, Width, radius, threshold);
            else
                SrcValuesGray = ImageProcess.DustAndScratches(SrcValuesGray, Width, radius, threshold);
        }
    }

    /// <summary>ガウシアンBlurを施す。Grayスケールモードの時だけ有効。</summary>
    /// <param name="hwhm">半値半幅</param>
    /// <param name="mode"></param>
    /// <param name="applyToOriginal"></param>
    public void SetBlurImage(double hwhm, BlurModeEnum mode, bool applyToOriginal = true)
    {
        if (IsSrcGray)
        {
            if (mode == BlurModeEnum.Gaussian)
            {
                if (applyToOriginal)
                    SrcValuesGray = ImageProcess.GaussianBlurFast(SrcValuesGrayOriginal, Width, hwhm);
                else
                    SrcValuesGray = ImageProcess.GaussianBlurFast(SrcValuesGray, Width, hwhm);
            }
            else
                SrcValuesGray = SrcValuesGrayOriginal;
        }
    }

    /// <summary>指定されたピクセルのカラー値を返す</summary>
    /// <param name="x">xピクセル値</param>
    /// <param name="y">yピクセル値</param>
    /// <returns></returns>
    [Obsolete("261111Cl: 4 アプリで未使用。将来削除予定")]
    public Color GetPixelColor(int x, int y)
    {
        if (x < 0 || Width <= x || y < 0 || Height <= y || this.RealImage == false) return Color.FromArgb(0, 0, 0);

        var offset = x + (y * Width);
        var index = Math.Max(Math.Min((int)((double)(SrcValuesR[offset] - MinValue) / (MaxValue - MinValue) * ColorScale.Length + 0.5), ColorScale.Length - 1), 0);
        return Color.FromArgb(ColorScale[index].R, ColorScale[index].G, ColorScale[index].B);
    }

    [Obsolete("261111Cl: 4 アプリで未使用。将来削除予定")]
    public Color GetPixelColor(double x, double y) => GetPixelColor((int)(x + 0.5), (int)(y + 0.5));

    public double GetPixelRawValue(PointD pt) => GetPixelRawValue(pt.X, pt.Y);

    public void SetScaleColdWarm(bool linear = true)
    {
        if (linear)
            ColorScale = ColorScaleColdWarmLiner;
        else
            ColorScale = ColorScaleColdWarmLog;

        GrayScale = false;
    }

    public void SetScaleSpectrum(bool linear = true)
    {
        if (linear)
            ColorScale = ColorScaleSpectrumLiner;
        else
            ColorScale = ColorScaleSpectrumLog;

        GrayScale = false;
    }
    public void SetScaleFire(bool linear = true)
    {
        if (linear)
            ColorScale = ColorScaleFireLiner;
        else
            ColorScale = ColorScaleFireLog;

        GrayScale = false;
    }

    public void SetScaleGray(bool linear = true)
    {
        if (linear)
            ColorScale = ColorScaleGrayLiner;
        else
            ColorScale = ColorScaleGrayLog;
        GrayScale = false;
    }

    public void SetScaleRotation()
    {
        ColorScale = ColorScaleRotationLiner;
        GrayScale = false;
    }

    /// <summary>グレースケールの時のx, y位置の生の強度を返す. グレースケールではなかったら常に0を返す.</summary>
    /// <param name="x"></param>
    /// <param name="y"></param>
    /// <returns></returns>
    public double GetPixelRawValue(double x, double y) => GetPixelRawValue((int)(x + 0.5), (int)(y + 0.5));

    /// <summary>グレースケールの時のx, y位置の生の強度を返す グレースケールではなかったら常に0を返す.</summary>
    /// <param name="x"></param>
    /// <param name="y"></param>
    /// <returns></returns>
    public double GetPixelRawValue(int x, int y)
    {
        int offset = x + y * Width;
        if (IsSrcGray == true && SrcValuesGray != null && offset < SrcValuesGray.Length && 0 <= offset)
            return SrcValuesGray[offset];
        else
            return 0;
    }

    /// <summary>指定されたピクセルの複素数値を返す</summary>
    /// <param name="x"></param>
    /// <param name="y"></param>
    /// <returns></returns>
    public Complex[] GetPixelComplex(int x, int y)
    {
        if (x < 0 || Width <= x || y < 0 || Height <= y || this.RealImage == true)
            return [new(), new(), new()];
        else
#pragma warning disable CS0618 // 261111Cl 追加: ComplexGray/R/G/B は [Obsolete] (4 アプリで未使用)
            return this.GrayScale ? ([ComplexGray[y][x], ComplexGray[y][x], ComplexGray[y][x]]) : ([ComplexR[y][x], ComplexG[y][x], ComplexB[y][x]]);
#pragma warning restore CS0618
    }

    public Complex[] GetPixelComplex(double x, double y) => GetPixelComplex((int)Math.Round(x), (int)Math.Round(y));

    /// <summary>等倍のイメージを取得する</summary>
    /// <param name="srcCenter"></param>
    /// <param name="zoom"></param>
    /// <param name="destSize"></param>
    /// <returns></returns>
    public Bitmap GetImage()
    {
        var bmp = GetImage(new RectangleD(0, 0, Width, Height), new Size(Width, Height));
        return bmp;
    }

    /// <summary>指定された範囲のイメージを取得する</summary>
    /// <param name="srcCenter"></param>
    /// <param name="zoom"></param>
    /// <param name="destSize"></param>
    /// <returns></returns>
    public Bitmap GetImage(PointD srcCenter, double zoom, Size destSize) => srcCenter.IsNaN
        ? null
        : GetImage(GetDrawingArea(srcCenter, zoom, destSize), destSize);

    /// <summary>描画範囲を得る</summary>
    /// <param name="srcCenter"></param>
    /// <param name="zoom"></param>
    /// <param name="destSize"></param>
    /// <returns></returns>
    public static RectangleD GetDrawingArea(PointD srcCenter, double zoom, Size destSize)
        => new(srcCenter.X - destSize.Width / zoom / 2.0, srcCenter.Y - destSize.Height / zoom / 2.0, destSize.Width / zoom, destSize.Height / zoom);

    private Bitmap destBmp;

    public delegate Point GetSrcPosition(int x, int y);

    public delegate byte GetValue(double value, byte[] scale);

    /// <summary>指定された範囲のイメージを取得する</summary>
    /// <param name="srcRect"></param>
    /// <param name="destSize"></param>
    /// <returns></returns>
    public unsafe Bitmap GetImage(RectangleD srcRect, Size destSize)
    {
        //if (Height == 0 || Width == 0 || srcRect.Width == 0 || srcRect.Height == 0 || destSize.Width == 0 || destSize.Height == 0) return null; //260731Cl 変更前
        if (Height <= 0 || Width <= 0 || srcRect.Width <= 0 || srcRect.Height <= 0 || destSize.Width <= 0 || destSize.Height <= 0) return null; //260731Cl 変更: レイアウト過渡で ClientSize 由来の destSize が負になり得る (new Bitmap の ArgumentException 防止)

        int width = destSize.Width, height = destSize.Height;

        //まずbmpが前回作られていなかったら作成
        //if (destBmp == null || destBmp.Width != width || destBmp.Height != height || destBmp.PixelFormat == (AlphaEnabled ? PixelFormat.Format32bppPArgb : PixelFormat.Format24bppRgb)) // 旧: 同じ形式のたびに再生成していた
        var pixelFormat = AlphaEnabled ? PixelFormat.Format32bppPArgb : PixelFormat.Format24bppRgb; // (260715Ch) キャッシュの期待形式を一度だけ決定
        if (destBmp == null || destBmp.Width != width || destBmp.Height != height || destBmp.PixelFormat != pixelFormat) // (260715Ch) サイズまたは形式が変わった場合だけ再生成
        {
            //destBmp = new Bitmap(width, height, AlphaEnabled ? PixelFormat.Format32bppPArgb : PixelFormat.Format24bppRgb); // 旧: 置換前の Bitmap が未解放
            //destBmp?.Dispose(); // (260715Ch) GDI ハンドルを失わないよう、キャッシュ差し替え前に解放 //260802Cl 変更前
            //destBmp = new Bitmap(width, height, pixelFormat); // (260715Ch) //260802Cl 変更前
            //260802Cl 変更: 新しいキャッシュを先に作ってから古い方を解放する。旧順序では「破棄済みなのに
            //まだ呼び出し側 (PictureBox.Image) から参照されている」区間が長く、new Bitmap が失敗した場合も
            //destBmp が破棄済みのまま残ってしまう
            var previous = destBmp;
            destBmp = new Bitmap(width, height, pixelFormat); // (260715Ch)
            previous?.Dispose();
        }
        //bmpをロック
        BitmapData bmpData;
        try
        {
            bmpData = destBmp.LockBits(new Rectangle(0, 0, destBmp.Width, destBmp.Height), ImageLockMode.ReadWrite, destBmp.PixelFormat);
        }
        catch
        {
            return null;
        }
        //261112Cl 追加: 画素ループなどで例外が出ても必ず UnlockBits する (旧はロックが残り、以後の LockBits が失敗して GetImage が null を返し続けた)。
        //  try の中身は字下げを変えずに残す (差分を小さくするため)。finally は下の UnlockBits の位置
        try
        {
        double zoom = (width / srcRect.Width + height / srcRect.Height) / 2.0;

        # region 描画する画素位置をソース画像位置に変換するローカル関数
        double srcX = srcRect.X + 0.5, srcY = srcRect.Y + 0.5, w = srcRect.Width / width, h = srcRect.Height / height;
        #region 261111Cl 変更前: 画素ごとに Func を呼んでいた
        //Func<int, int, (int X, int Y)> getSrcPosition;
        //if (HorizontalFlip & VerticalFlip)
        //    getSrcPosition = (x, y) => ((int)(srcX + (width - x) * w), (int)(srcY + (height - y) * h));
        //else if (!HorizontalFlip & VerticalFlip)
        //    getSrcPosition = (x, y) => ((int)(srcX + x * w), (int)(srcY + (height - y) * h));
        //else if (HorizontalFlip & !VerticalFlip)
        //    getSrcPosition = (x, y) => ((int)(srcX + (width - x) * w), (int)(srcY + y * h));
        //else
        //    getSrcPosition = (x, y) => ((int)(srcX + x * w), (int)(srcY + y * h));
        #endregion
        //261111Cl 変更: X の元座標は表 (width 個) を 1 回だけ作り、Y は行ごとに計算する (スレッド本体の中)。式は従来と同じ
        var srcXs = new int[width];
        if (HorizontalFlip)
            for (int x = 0; x < width; x++) srcXs[x] = (int)(srcX + (width - x) * w);
        else
            for (int x = 0; x < width; x++) srcXs[x] = (int)(srcX + x * w);
        bool verticalFlip = VerticalFlip;
        #endregion

        #region  入力値doubleを表示する値に変換するローカル関数
        //261111Cl 変更 (作者承認・.NET 8 以前の見た目に戻す): (int)x は .NET 8 以前の x64 では cvttsd2si で、NaN と int の範囲外 (正負とも) を int.MinValue にしていた
        //  (下の丸めで添字 0 = 最暗色)。.NET 9 以降は飽和変換になり、正の範囲外 (MaxValue == MinValue のときの +∞ など) が int.MaxValue (最明色) に変わった。
        //  旧の変換を明示的に書いて、NaN と範囲外は int.MinValue (添字 0) に戻す。範囲内 (-2^31 ≦ x < 2^31) の値は (int)x のまま。
        var length = ColorScale.Length;
        var coeff = length / (MaxValue - MinValue);
        byte getValueR(double rawValue)
        {
            //var rawIndex = (int)((rawValue - MinValue) * coeff + 0.5);
            var x = (rawValue - MinValue) * coeff + 0.5;// 261111Cl 変更
            var rawIndex = x is >= -2147483648.0 and < 2147483648.0 ? (int)x : int.MinValue;// 261111Cl 変更: .NET 8 以前の (int)x と同じ
            return ColorScale[Math.Min(Math.Max(0, rawIndex), length - 1)].R;
        }

        byte getValueG(double rawValue)
        {
            //var rawIndex = (int)((rawValue - MinValue) * coeff + 0.5);
            var x = (rawValue - MinValue) * coeff + 0.5;// 261111Cl 変更
            var rawIndex = x is >= -2147483648.0 and < 2147483648.0 ? (int)x : int.MinValue;// 261111Cl 変更: .NET 8 以前の (int)x と同じ
            return ColorScale[Math.Min(Math.Max(0, rawIndex), length - 1)].G;
        }

        byte getValueB(double rawValue)
        {
            //var rawIndex = (int)((rawValue - MinValue) * coeff + 0.5);
            var x = (rawValue - MinValue) * coeff + 0.5;// 261111Cl 変更
            var rawIndex = x is >= -2147483648.0 and < 2147483648.0 ? (int)x : int.MinValue;// 261111Cl 変更: .NET 8 以前の (int)x と同じ
            return ColorScale[Math.Min(Math.Max(0, rawIndex), length - 1)].B;
        }

        //261111Cl 追加: 色表の添字 (getValueR/G/B と同じ式)。グレーの源では R/G/B で 3 回計算していたのを 1 回にする
        [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.AggressiveInlining)]
        static int toIndex(double rawValue, double minValue, double coeff, int length)
        {
            var x = (rawValue - minValue) * coeff + 0.5;
            var rawIndex = x is >= -2147483648.0 and < 2147483648.0 ? (int)x : int.MinValue;// .NET 8 以前の (int)x と同じ (NaN・範囲外は添字 0)
            return Math.Min(Math.Max(0, rawIndex), length - 1);
        }
        #endregion

        //int range = (int)(1 / zoom / 2 + 0.5);// 261112Cl 変更前
        // 261112Cl 変更: zoom がほぼ 0 (描画範囲が出力の約 40 億倍) だと range が int.MaxValue に張り付き、下の sy + range があふれて標本数が 0 になり、
        //  前の画素の色が残っていた。あふれない大きさ (int.MaxValue / 4) で頭打ちにする。それ未満では旧と同じ値
        int range = (int)Math.Min(1 / zoom / 2 + 0.5, int.MaxValue / 4);
        int increase = (zoom >= 1) ? 0 : Math.Max(range / 3, 1);

        int step = AlphaEnabled ? 4 : 3;
        int nResidual = bmpData.Stride - destBmp.Width * step;
        //double filter; //261111Cl 変更前: 全スレッドで 1 個を共有していた (FFT_Filter が 0 でないと競合する)。スレッド本体の中へ移した

        int thread = Environment.ProcessorCount;


#if DEBUG
        //thread = 1;
#endif
        int yStep = height / thread + 1;
        #region 261111Cl 変更前 (画素ループ)
        //Parallel.For(0, thread, t =>
        //{
        //    byte* p = (byte*)(void*)bmpData.Scan0;
        //    p += t * yStep * bmpData.Stride;
        //
        //    byte r = 0, g = 0, b = 0;
        //
        //    (int X, int Y) beforePt = (int.MinValue, int.MinValue);
        //    for (int y = t * yStep; y < Math.Min((t + 1) * yStep, height); y++)
        //    {
        //        for (int x = 0; x < width; x++)
        //        {
        //            var pt = getSrcPosition(x, y);
        //            if (pt.X == beforePt.X && pt.Y == beforePt.Y)
        //            {
        //                for (int i = 0; i < step; i++)
        //                    p[i] = p[i - step];
        //            }
        //            else if ((uint)pt.X < (uint)Width && (uint)pt.Y < (uint)Height)//描画域内のとき
        //            {
        //                var srcPosition = pt.X + pt.Y * Width;
        //                if (increase > 0)//補完モード
        //                {
        //                    double rSum = 0, gSum = 0, bSum = 0, count = 0;
        //
        //                    for (int j = Math.Max(0, pt.Y - range); j <= Math.Min(pt.Y + range, Height - 1); j += increase)
        //                        for (int i = Math.Max(0, pt.X - range); i <= Math.Min(pt.X + range, Width - 1); i += increase)
        //                        {
        //                            var tempPosition = i + j * Width;
        //
        //                            if (IsSrcGray)
        //                                bSum += SrcValuesGray[tempPosition];
        //                            else
        //                            {
//#pragma warning disable CS0618 // 261111Cl 追加: SrcValuesR/G/B は [Obsolete] (4 アプリで未使用)
        //                                bSum += SrcValuesB[tempPosition];
        //                                if (!GrayScale)
        //                                {
        //                                    gSum += SrcValuesG[tempPosition];
        //                                    rSum += SrcValuesR[tempPosition];
        //                                }
//#pragma warning restore CS0618
        //                            }
        //                            count++;
        //                        }
        //                    if (count > 0)
        //                    {
        //                        b = getValueB(bSum / count);
        //                        if (GrayScale)
        //                            g = r = b;
        //                        else
        //                        {
        //                            if (IsSrcGray)
        //                                gSum = rSum = bSum;
        //                            g = getValueG(gSum / count);
        //                            r = getValueR(rSum / count);
        //                        }
        //                    }
        //                }
        //                else
        //                {
        //                    if (IsSrcGray)
        //                    {
        //                        b = getValueB(SrcValuesGray[srcPosition]);
        //                        if (GrayScale)
        //                            g = r = b;
        //                        else
        //                        {
        //                            g = getValueG(SrcValuesGray[srcPosition]);
        //                            r = getValueR(SrcValuesGray[srcPosition]);
        //                        }
        //                    }
        //                    else
        //                    {
//#pragma warning disable CS0618 // 261111Cl 追加: SrcValuesR/G/B は [Obsolete] (4 アプリで未使用)
        //                        b = getValueB(SrcValuesB[srcPosition]);
        //                        if (GrayScale)
        //                            g = r = b;
        //                        else
        //                        {
        //                            g = getValueG(SrcValuesG[srcPosition]);
        //                            r = getValueR(SrcValuesR[srcPosition]);
        //                        }
//#pragma warning restore CS0618
        //                    }
        //                }
        //
        //                if (IsNegative)
        //                {
        //                    b = (byte)(byte.MaxValue - b); g = (byte)(byte.MaxValue - g); r = (byte)(byte.MaxValue - r);
        //                }
        //
        //                p[0] = b;
        //                p[1] = g;
        //                p[2] = r;
        //
        //                if (AlphaEnabled)
        //                    p[3] = (srcPosition < FilterAlfha.Count) ? FilterAlfha[srcPosition] : (byte)255;
        //
        //                #region 各種フィルター
//#pragma warning disable CS0618 // 261111Cl 追加: FFT_Filter は [Obsolete] (4 アプリで未使用)
        //                if (FFT_Filter != null && FFT_Filter.Count > 0 && (filter = FFT_Filter[srcPosition]) != 0)//もしフィルターされていたら全体を1/2にして(0,127,127)をたす
//#pragma warning restore CS0618
        //                { //BG
        //                    p[0] = (byte)(p[0] * (2 - filter) * 0.5 + 127 * filter);
        //                    p[1] = (byte)(p[1] * (2 - filter) * 0.5 + 127 * filter);
        //                    p[2] = (byte)(p[2] * (2 - filter) * 0.5);
        //                }
        //
        //                if (FilterTemporary != null && FilterTemporary.Count > 0 && FilterTemporary[srcPosition])//もし一時フィルターされていたら全体を1/2にして(127,127,0)をたす
        //                {//GR
        //                    p[0] = (byte)(p[0] * 0.8);
        //                    p[1] = (byte)(p[1] * 0.8 + 51);
        //                    p[2] = (byte)(p[2] * 0.8 + 51);
        //                }
        //                if (Filter1Visible && Filter1 != null && Filter1[srcPosition])//もし飽和点にいろをつけるなら
        //                {//R
        //                    p[0] = (byte)(p[0] * 0.8);
        //                    p[1] = (byte)(p[1] * 0.8);
        //                    p[2] = (byte)(p[2] * 0.8 + 51);
        //                }
        //                if (Filter2Visible && Filter2 != null && Filter2[srcPosition])
        //                {//B
        //                    p[0] = (byte)(p[0] * 0.8 + 51);
        //                    p[1] = (byte)(p[1] * 0.8);
        //                    p[2] = (byte)(p[2] * 0.8);
        //                }
        //                if (Filter3Visible && Filter3 != null && Filter3[srcPosition])//もしスポットだったら全体を1/2にして(0,127,127)をたす
        //                {//BG
        //                    p[0] = (byte)(p[0] * 0.8 + 51);
        //                    p[1] = (byte)(p[1] * 0.8 + 51);
        //                    p[2] = (byte)(p[2] * 0.8);
        //                }
        //                if (Filter4Visible && Filter4 != null && Filter4[srcPosition])//もし範囲外だったら全体を2/3にして(0,85,0)をたす
        //                {//G
        //                    p[0] = (byte)(p[0] * 0.8);
        //                    p[1] = (byte)(p[1] * 0.8 + 51);
        //                    p[2] = (byte)(p[2] * 0.8);
        //                }
        //                if (Filter5Visible && Filter5 != null && Filter5[srcPosition])//もし範囲外だったら全体を2/3にして(0,85,0)をたす
        //                {//BR
        //                    // 260718Cl: p[1]=p[2]*0.8; p[2]=(書換後の)p[1]*0.8+51 という自己二重依存のチャネル取り違えを、他フィルタと同じ自チャネル参照に修正
        //                    p[0] = (byte)(p[0] * 0.8 + 51);
        //                    p[1] = (byte)(p[1] * 0.8);
        //                    p[2] = (byte)(p[2] * 0.8 + 51);
        //                }
        //                #endregion
        //
        //            }
        //            else//描画域がはみ出したとき
        //            {
        //                p[0] = 0; p[1] = 127; p[2] = 0;
        //                if (AlphaEnabled)
        //                    p[3] = (byte)255;
        //            }
        //
        //            p += step;
        //
        //            beforePt = pt;
        //        }
        //        p += nResidual;
        //    }
        //});
        #endregion

        //261111Cl 変更: 画素ループの改修 (出力は従来と同じ。ただし下の不具合修正の分だけ変わる)
        //  - 元座標は X の表 (srcXs) と行ごとの Y で求める (画素ごとの Func 呼び出しをやめた)
        //  - グレーの源では色表の添字を 1 回だけ計算して色を 1 回読む (旧は R/G/B で 3 回)
        //  - Filter1〜5・FilterTemporary・FilterAlfha は List<T> のインデクサでなく Span で読む
        //  - 不具合修正 (挙動が変わる): 「前の画素と同じ元座標なら複写」の比較 (旧 beforePt) を行ごとにリセットする。旧は行をまたいで
        //    持ち越したため、1 行全体が元画像の 1 列に入るほど拡大すると、行頭で p[i - step] が前の行の stride の詰め物を読んでいた
        //  - 直接モード (縮小なし) で元画像の行が前の行と同じなら、出力も同じなので行ごと複写する
        var colorScale = ColorScale;
        var srcGray = SrcValuesGray;
        double minValue = MinValue;
        int srcWidth = Width, srcHeight = Height, stride = bmpData.Stride, rowBytes = width * step;
        bool isSrcGray = IsSrcGray, grayScale = GrayScale, isNegative = IsNegative, alphaEnabled = AlphaEnabled;
        Parallel.For(0, thread, t =>
        {
            // Span はラムダに取り込めないので、スレッド本体の中で取る
            ReadOnlySpan<int> xs = srcXs;
            ReadOnlySpan<bool> filterTemporary = CollectionsMarshal.AsSpan(FilterTemporary);
            ReadOnlySpan<bool> filter1 = CollectionsMarshal.AsSpan(Filter1), filter2 = CollectionsMarshal.AsSpan(Filter2), filter3 = CollectionsMarshal.AsSpan(Filter3),
                filter4 = CollectionsMarshal.AsSpan(Filter4), filter5 = CollectionsMarshal.AsSpan(Filter5);
            bool useFilter1 = Filter1Visible && Filter1 != null, useFilter2 = Filter2Visible && Filter2 != null, useFilter3 = Filter3Visible && Filter3 != null,
                useFilter4 = Filter4Visible && Filter4 != null, useFilter5 = Filter5Visible && Filter5 != null;
            ReadOnlySpan<byte> filterAlpha = CollectionsMarshal.AsSpan(FilterAlfha);
            double filter;

            byte* p = (byte*)(void*)bmpData.Scan0;
            p += t * yStep * stride;

            byte r = 0, g = 0, b = 0;

            int yStart = t * yStep, yEnd = Math.Min((t + 1) * yStep, height), beforeSy = 0;
            for (int y = yStart; y < yEnd; y++)
            {
                int sy = verticalFlip ? (int)(srcY + (height - y) * h) : (int)(srcY + y * h);
                if (increase == 0 && y > yStart && sy == beforeSy)
                {// 直接モードで元画像の行が前の行と同じ → 前の行をそのまま複写
                    Buffer.MemoryCopy(p - stride, p, rowBytes, rowBytes);
                    p += stride;
                    continue;
                }
                beforeSy = sy;
                bool rowInside = (uint)sy < (uint)srcHeight;
                int rowOffset = rowInside ? sy * srcWidth : 0;

                int beforeSx = 0;
                for (int x = 0; x < xs.Length; x++)
                {
                    int sx = xs[x];
                    if (x > 0 && sx == beforeSx)// 前の画素と同じ元座標なら複写 (比較は行ごとにリセット)
                    {
                        for (int i = 0; i < step; i++)
                            p[i] = p[i - step];
                    }
                    else if (rowInside && (uint)sx < (uint)srcWidth)//描画域内のとき
                    {
                        var srcPosition = sx + rowOffset;
                        if (increase > 0)//補完モード
                        {
                            double rSum = 0, gSum = 0, bSum = 0, count = 0;
                            int jEnd = Math.Min(sy + range, srcHeight - 1), iStart = Math.Max(0, sx - range), iEnd = Math.Min(sx + range, srcWidth - 1);
                            for (int j = Math.Max(0, sy - range); j <= jEnd; j += increase)
                                for (int i = iStart; i <= iEnd; i += increase)
                                {
                                    var tempPosition = i + j * srcWidth;

                                    if (isSrcGray)
                                        bSum += srcGray[tempPosition];
                                    else
                                    {
#pragma warning disable CS0618 // SrcValuesR/G/B は [Obsolete] (4 アプリで未使用)。扱いは従来どおり
                                        bSum += SrcValuesB[tempPosition];
                                        if (!grayScale)
                                        {
                                            gSum += SrcValuesG[tempPosition];
                                            rSum += SrcValuesR[tempPosition];
                                        }
#pragma warning restore CS0618
                                    }
                                    count++;
                                }
                            if (count > 0)
                            {
                                if (isSrcGray)
                                {// 旧は gSum = rSum = bSum として R/G/B で同じ添字を 3 回計算していた
                                    var c = colorScale[toIndex(bSum / count, minValue, coeff, length)];
                                    if (grayScale)
                                        g = r = b = c.B;
                                    else
                                    {
                                        b = c.B; g = c.G; r = c.R;
                                    }
                                }
                                else
                                {
                                    b = getValueB(bSum / count);
                                    if (grayScale)
                                        g = r = b;
                                    else
                                    {
                                        g = getValueG(gSum / count);
                                        r = getValueR(rSum / count);
                                    }
                                }
                            }
                        }
                        else
                        {
                            if (isSrcGray)
                            {
                                var c = colorScale[toIndex(srcGray[srcPosition], minValue, coeff, length)];
                                if (grayScale)
                                    g = r = b = c.B;
                                else
                                {
                                    b = c.B; g = c.G; r = c.R;
                                }
                            }
                            else
                            {
#pragma warning disable CS0618 // SrcValuesR/G/B は [Obsolete] (4 アプリで未使用)。扱いは従来どおり
                                b = getValueB(SrcValuesB[srcPosition]);
                                if (grayScale)
                                    g = r = b;
                                else
                                {
                                    g = getValueG(SrcValuesG[srcPosition]);
                                    r = getValueR(SrcValuesR[srcPosition]);
                                }
#pragma warning restore CS0618
                            }
                        }

                        if (isNegative)
                        {
                            b = (byte)(byte.MaxValue - b); g = (byte)(byte.MaxValue - g); r = (byte)(byte.MaxValue - r);
                        }

                        p[0] = b;
                        p[1] = g;
                        p[2] = r;

                        if (alphaEnabled)
                            p[3] = (srcPosition < filterAlpha.Length) ? filterAlpha[srcPosition] : (byte)255;

                        #region 各種フィルター
#pragma warning disable CS0618 // FFT_Filter は [Obsolete] (4 アプリで未使用)。扱いは従来どおり
                        if (FFT_Filter != null && FFT_Filter.Count > 0 && (filter = FFT_Filter[srcPosition]) != 0)//もしフィルターされていたら全体を1/2にして(0,127,127)をたす
#pragma warning restore CS0618
                        { //BG
                            p[0] = (byte)(p[0] * (2 - filter) * 0.5 + 127 * filter);
                            p[1] = (byte)(p[1] * (2 - filter) * 0.5 + 127 * filter);
                            p[2] = (byte)(p[2] * (2 - filter) * 0.5);
                        }

                        if (filterTemporary.Length > 0 && filterTemporary[srcPosition])//もし一時フィルターされていたら全体を1/2にして(127,127,0)をたす
                        {//GR
                            p[0] = (byte)(p[0] * 0.8);
                            p[1] = (byte)(p[1] * 0.8 + 51);
                            p[2] = (byte)(p[2] * 0.8 + 51);
                        }
                        if (useFilter1 && filter1[srcPosition])//もし飽和点にいろをつけるなら
                        {//R
                            p[0] = (byte)(p[0] * 0.8);
                            p[1] = (byte)(p[1] * 0.8);
                            p[2] = (byte)(p[2] * 0.8 + 51);
                        }
                        if (useFilter2 && filter2[srcPosition])
                        {//B
                            p[0] = (byte)(p[0] * 0.8 + 51);
                            p[1] = (byte)(p[1] * 0.8);
                            p[2] = (byte)(p[2] * 0.8);
                        }
                        if (useFilter3 && filter3[srcPosition])//もしスポットだったら全体を1/2にして(0,127,127)をたす
                        {//BG
                            p[0] = (byte)(p[0] * 0.8 + 51);
                            p[1] = (byte)(p[1] * 0.8 + 51);
                            p[2] = (byte)(p[2] * 0.8);
                        }
                        if (useFilter4 && filter4[srcPosition])//もし範囲外だったら全体を2/3にして(0,85,0)をたす
                        {//G
                            p[0] = (byte)(p[0] * 0.8);
                            p[1] = (byte)(p[1] * 0.8 + 51);
                            p[2] = (byte)(p[2] * 0.8);
                        }
                        if (useFilter5 && filter5[srcPosition])//もし範囲外だったら全体を2/3にして(0,85,0)をたす
                        {//BR
                            // 260718Cl: p[1]=p[2]*0.8; p[2]=(書換後の)p[1]*0.8+51 という自己二重依存のチャネル取り違えを、他フィルタと同じ自チャネル参照に修正
                            p[0] = (byte)(p[0] * 0.8 + 51);
                            p[1] = (byte)(p[1] * 0.8);
                            p[2] = (byte)(p[2] * 0.8 + 51);
                        }
                        #endregion

                    }
                    else//描画域がはみ出したとき
                    {
                        p[0] = 0; p[1] = 127; p[2] = 0;
                        if (alphaEnabled)
                            p[3] = (byte)255;
                    }

                    p += step;

                    beforeSx = sx;
                }
                p += nResidual;
            }
        });

        //FilterTemporary.Clear();
        //FilterTemporary.AddRange(new bool[Filter1.Count]);
        //261111Cl 変更: 再描画のたびに Filter1.Count 個の bool[] (数千×数千画素で 16 MB) を確保していたのをやめ、要素数をそろえてその場で false に戻す (結果は同じ: 全要素 false・要素数 Filter1.Count)
        FilterTemporary ??= [];// 261112Cl 追加: public フィールドが null だと SetCount で例外になっていた
        //CollectionsMarshal.SetCount(FilterTemporary, Filter1.Count);// 261112Cl 変更前
        CollectionsMarshal.SetCount(FilterTemporary, Filter1?.Count ?? 0);// 261112Cl 変更: Filter1 が null でも例外にしない (画素ループは Filter1 != null を確かめている)
        CollectionsMarshal.AsSpan(FilterTemporary).Clear();

        }
        finally
        {
        //destBmp.UnlockBits(bmpData); //261112Cl 変更前
        destBmp.UnlockBits(bmpData); //261112Cl 変更: finally の中へ (上の try を参照)
        }

        return destBmp;
    }

    /// <summary>実画像モードのとき回転操作を施す</summary>
    /// <param name="bmp"></param>
    /// <param name="angle"></param>
    /// <returns></returns>
    [Obsolete("261111Cl: 4 アプリで未使用。将来削除予定")]
    public void Rotate(double angle)
    {
        //var destBitmap = new Bitmap(SrcBitmap.Width, SrcBitmap.Height, PixelFormat.Format24bppRgb); // 旧: ローカル Bitmap が未解放
        using var destBitmap = new Bitmap(SrcBitmap.Width, SrcBitmap.Height, PixelFormat.Format24bppRgb); // (260715Ch)
        //var g = Graphics.FromImage(destBitmap); // 旧: Graphics が未解放
        using (var g = Graphics.FromImage(destBitmap)) // (260715Ch) LockBits 前に描画を完了して解放
        {
            g.SmoothingMode = System.Drawing.Drawing2D.SmoothingMode.HighQuality;
            g.TranslateTransform((float)SrcBitmap.Width / 2, (float)SrcBitmap.Height / 2);
            g.RotateTransform((float)(angle * 180 / Math.PI));
            g.TranslateTransform(-(float)SrcBitmap.Width / 2, -(float)SrcBitmap.Height / 2);
            g.Clear(Color.White);
            g.DrawImage(SrcBitmap, new RectangleF(0, 0, SrcBitmap.Width, SrcBitmap.Height), new RectangleF(0, 0, SrcBitmap.Width, SrcBitmap.Height), GraphicsUnit.Pixel);
        }

        //SourceBitmapをロックする
        var bitmapData = destBitmap.LockBits(new Rectangle(0, 0, Width, Height), ImageLockMode.ReadWrite, SrcBitmap.PixelFormat);
        //SouceStrideを設定しておく
        int BitmapStride = bitmapData.Stride;
        //次にその情報をargbValuesに転送
        byte[] bitmapByteArray = GC.AllocateUninitializedArray<byte>(bitmapData.Stride * SrcBitmap.Height);
        Marshal.Copy(bitmapData.Scan0, bitmapByteArray, 0, bitmapByteArray.Length);
        //SourceBitmapをアンロックして終了
        destBitmap.UnlockBits(bitmapData);

        if (GrayScale)
        {
            SrcValuesR.Clear();

            for (int yInit = 0; yInit < Height * bitmapData.Stride; yInit += BitmapStride)
                for (int x = 0; x < Width; x++)
                    SrcValuesR.Add(bitmapByteArray[yInit + x * 3]);
            SrcValuesG = SrcValuesR;
            SrcValuesB = SrcValuesR;
        }
        else
        {
            SrcValuesR.Clear();
            SrcValuesG.Clear();
            SrcValuesB.Clear();
            for (int yInit = 0; yInit < Height * bitmapData.Stride; yInit += BitmapStride)
                for (int x = 0; x < Width; x++)
                {
                    SrcValuesR.Add(bitmapByteArray[yInit + x * 3 + 2]);
                    SrcValuesG.Add(bitmapByteArray[yInit + x * 3 + 1]);
                    SrcValuesB.Add(bitmapByteArray[yInit + x * 3 + 0]);
                }
        }
        bitmapByteArray = null;
    }



    #endregion 操作
}
