using System;
using System.Buffers;// 261111Cl 追加
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;// 261111Cl 追加
using System.Text;

namespace Crystallography;

public class Tiff
{
    /// <summary>マルチバイトデータのバイト順 正順の場合はMotorola, 逆順はIntel</summary>
    public enum TiffByteOrder
    {
        Motorola, Intel
    }

    /// <summary>データタイプ</summary>
    public enum TiffDataType
    {
        /// <summary>1バイト整数</summary>
        ByteType,

        AsciiType, //１バイトのASCII文字
        ShortType,//２バイト短整数
        LongType,//(４バイト長整数)
        RationalType,//(８バイト分数、４バイトの分子とそれに続く４バイトの分母)
        Sbyte,//１バイト符号付き整数
        UndefinedType,//(あらゆる１バイトデータ)
        SshortType,//２バイト符号付き短整数)
        SlongType,//(４バイト符号付き長整数)
        Srational,//(８バイト符号付き分数、４バイトの分子とそれに続く４バイトの分母)
        FloatType,//(４バイト実数、IEEE浮動小数点形式)
        DoubleType//(８バイト倍精度実数、IEEE倍精度浮動小数点形式)
    }

    public struct IFD(int tag, Type dataType, object[] data)
    {
        public int Tag = tag;
        public Type DataType = dataType;
        public object[] Data = data;

        public static void Write(BinaryWriter bw, ushort tag, byte value)
        {
            bw.Write(BitConverter.GetBytes(tag));
            bw.Write(BitConverter.GetBytes((ushort)1));//type
            bw.Write(BitConverter.GetBytes((uint)1));//count
            bw.Write(value); bw.Write((byte)0x00); bw.Write((byte)0x00); bw.Write((byte)0x00);
        }

        public static void Write(BinaryWriter bw, ushort tag, char value)
        {
            bw.Write(BitConverter.GetBytes(tag));
            bw.Write(BitConverter.GetBytes((ushort)2));//type
            bw.Write(BitConverter.GetBytes((uint)1));//count
            bw.Write(value); bw.Write((byte)0x00); bw.Write((byte)0x00); bw.Write((byte)0x00);
        }

        public static void Write(BinaryWriter bw, ushort tag, ushort value)
        {
            bw.Write(BitConverter.GetBytes(tag));
            bw.Write(BitConverter.GetBytes((ushort)3));//type
            bw.Write(BitConverter.GetBytes((uint)1));//count
            bw.Write(BitConverter.GetBytes(value)); bw.Write((byte)0x00); bw.Write((byte)0x00);
        }

        public static void Write(BinaryWriter bw, ushort tag, uint value)
        {
            bw.Write(BitConverter.GetBytes(tag));
            bw.Write(BitConverter.GetBytes((ushort)4));//type
            bw.Write(BitConverter.GetBytes((uint)1));//count
            bw.Write(BitConverter.GetBytes(value));
        }

        public static void Write(BinaryWriter bw, ushort tag, float value)
        {
            bw.Write(BitConverter.GetBytes(tag));
            bw.Write(BitConverter.GetBytes((ushort)11));//type
            bw.Write(BitConverter.GetBytes((uint)1));//count
            bw.Write(BitConverter.GetBytes(value));
        }

        public static uint Write(BinaryWriter bw, ushort tag, double value, uint offset)
        {
            bw.Write(BitConverter.GetBytes(tag));
            bw.Write(BitConverter.GetBytes((ushort)12));//type
            bw.Write(BitConverter.GetBytes((uint)1));//count
            bw.Write(BitConverter.GetBytes(offset));
            long pos = bw.BaseStream.Position;
            bw.BaseStream.Position = offset;
            bw.Write(BitConverter.GetBytes(value));
            bw.BaseStream.Position = pos;
            return offset + 8;
        }

        public static uint Write(BinaryWriter bw, ushort tag, uint value1, uint value2, uint offset)
        {
            bw.Write(BitConverter.GetBytes(tag));
            bw.Write(BitConverter.GetBytes((ushort)5));//type
            bw.Write(BitConverter.GetBytes((uint)1));//count
            bw.Write(BitConverter.GetBytes(offset));
            long pos = bw.BaseStream.Position;
            bw.BaseStream.Position = offset;
            bw.Write(BitConverter.GetBytes(value1));
            bw.Write(BitConverter.GetBytes(value2));
            bw.BaseStream.Position = pos;
            return offset + 8;
        }

        public static uint Write(BinaryWriter bw, ushort tag, string str, uint offset)
        {
            bw.Write(BitConverter.GetBytes(tag));
            bw.Write(BitConverter.GetBytes((ushort)2));//type
            bw.Write(BitConverter.GetBytes((uint)str.Length + 1));//count
            // 261124Cl 追加 (codex2 の査読): 終端込み 4 バイト以下の文字列は、TIFF の決まりどおり値の欄に直接書く (読み手はそう読む。旧は常にオフセットを書いていた)
            if (str.Length + 1 <= 4)
            {
                var inline = new byte[4];
                for (int i = 0; i < str.Length; i++) inline[i] = (byte)str[i];
                bw.Write(inline);
                return offset;
            }
            bw.Write(BitConverter.GetBytes(offset));
            long pos = bw.BaseStream.Position;
            bw.BaseStream.Position = offset;
            for (int i = 0; i < str.Length; i++)
            {
                bw.Write((byte)str[i]);
                //bw.BaseStream.Position -= 1;
            }

            bw.Write((byte)0);

            bw.BaseStream.Position = pos;
            return offset + (uint)str.Length + 1;
        }
    }

    public static void Writer(string filename, double[] srcData, int sampleFormat, int imageWidth, IFD[] additionalIFD = null)
    {
        if (additionalIFD != null)
            //Writer(filename, [srcData], sampleFormat, imageWidth, [additionalIFD]); // 261124Cl 変更前 (codex2 の査読): タグの配列をフレームの配列として包んでいたので、タグが 1 個以外だとフレーム数の検査で黙って保存しなかった
            Writer(filename, [srcData], sampleFormat, imageWidth, [.. additionalIFD.Select(f => new[] { f })]); // 261124Cl 変更: タグごとに 1 フレームの配列
        else
            Writer(filename, [srcData], sampleFormat, imageWidth);
    }

    /// <summary>TIFFデータを書き込む. SampeleFormatが1のときは整数(byte, ushort, uint)、3の時はfloat, 1と3以外: 全て整数の場合は整数, そうでなければfloatに変換</summary>
    /// <typeparam name="T"></typeparam>
    /// <param name="srcData"></param>
    /// <param name="imageWidth"></param>
    /// <param name="filename"></param>
    /// <param name="sampleFormat">1のときは整数(byte, ushort, uint)、3の時はfloat, 1と3以外: 全て整数の場合は整数, そうでなければfloatに変換</param>
    public static void Writer(string filename, double[][] srcData, int sampleFormat, int imageWidth, IFD[][] additionalIFD = null)
    {
        double pixelSizeX = 1;
        double pixelSizeY = 1;
        //string unit = "";

        if (srcData == null || srcData.Length == 0 || srcData[0].Length == 0 || imageWidth == 0 || srcData[0].Length % imageWidth != 0) return;
        //additionalIFDの正当性チェック
        if (additionalIFD != null)
            for (int i = 0; i < additionalIFD.Length; i++)
                if (additionalIFD[i].Length != srcData.Length)
                    return;

        if (sampleFormat != 3 && sampleFormat != 1)
        {
            //if (!srcData[0].Any(d => d % 1 != 0))
            //261111Cl 変更: 「全部が整数か」を double.IsInteger で判定する (NaN と ±∞ は d % 1 = NaN ≠ 0 で「整数でない」、IsInteger も false なので判定は同じ。-0.0 はどちらも整数)
            if (Array.TrueForAll(srcData[0], double.IsInteger))
                sampleFormat = 1;
            // else if (srcData[0].Max() < float.MaxValue && srcData[0].Min() > float.MinValue) // 261111Cl 変更前 (SimdLinq。NaN の扱いが位置で変わっていた)
            else if (srcData[0].MinMaxIgnoreNaN() is var (srcMin, srcMax) && srcMax < float.MaxValue && srcMin > float.MinValue) // 261111Cl 変更: NaN を除いた最小・最大 (作者判断。NaN は float でも表せる)
                sampleFormat = 3;
            else
                return;
        }
        //var bw = new BinaryWriter(new FileStream(filename, FileMode.Create, FileAccess.ReadWrite));
        using var bw = new BinaryWriter(new FileStream(filename, FileMode.Create, FileAccess.ReadWrite));// 261111Cl 変更: 例外時もファイルを閉じる

        bw.Write((byte)0x49); bw.Write((byte)0x49);//バイトオーダー
        bw.Write((byte)0x2A); bw.Write((byte)0x00);//バージョン
        uint offset = (uint)bw.BaseStream.Position + 4;

        for (int n = 0; n < srcData.Length; n++)
        {
            int bytePerPixel = 4;
            long stripOffsets = offset;
            if (sampleFormat == 1)
            {
                var data = srcData[n].Select(d => d < uint.MaxValue && d > uint.MinValue ? (uint)d : (uint)0).ToArray();
                uint max = data.Max();
                if (max < 256)
                    bytePerPixel = 1;
                else if (max < 65536)
                    bytePerPixel = 2;
                else
                    bytePerPixel = 4;
                //IFDのポインタを書き込み
                bw.Write(BitConverter.GetBytes((uint)(offset + data.Length * bytePerPixel)));

                bw.BaseStream.Position = stripOffsets;
                //if (bytePerPixel == 1)
                //    for (int i = 0; i < data.Length; i++)
                //        bw.Write((byte)data[i]);
                //else if (bytePerPixel == 2)
                //    for (int i = 0; i < data.Length; i++)
                //        bw.Write(BitConverter.GetBytes((ushort)data[i]));
                //else
                //    for (int i = 0; i < data.Length; i++)
                //        bw.Write(BitConverter.GetBytes(data[i]));
                //261111Cl 変更: 画素ごとの BitConverter.GetBytes をやめ、型付きの配列にしてからバイト列として一度に書く (リトルエンディアン、出力はバイト単位で同一)
                if (bytePerPixel == 1)
                {
                    var buf = GC.AllocateUninitializedArray<byte>(data.Length);
                    for (int i = 0; i < data.Length; i++)
                        buf[i] = (byte)data[i];
                    bw.Write(buf);
                }
                else if (bytePerPixel == 2)
                {
                    var buf = GC.AllocateUninitializedArray<ushort>(data.Length);
                    for (int i = 0; i < data.Length; i++)
                        buf[i] = (ushort)data[i];
                    bw.Write(MemoryMarshal.AsBytes(buf.AsSpan()));
                }
                else
                    bw.Write(MemoryMarshal.AsBytes(data.AsSpan()));
            }
            else if (sampleFormat == 3)
            {
                var data = srcData[n].Select(d => d < float.MaxValue && d > float.MinValue ? (float)d : 0).ToArray();
                //IFDのポインタを書き込み
                bw.Write(BitConverter.GetBytes((uint)(offset + data.Length * bytePerPixel)));
                bw.BaseStream.Position = stripOffsets;
                //for (int i = 0; i < data.Length; i++)
                //    bw.Write(BitConverter.GetBytes(data[i]));
                bw.Write(MemoryMarshal.AsBytes(data.AsSpan()));// 261111Cl 変更: float 配列をバイト列として一度に書く
            }

            //IFDの総数
            short totalIFD = 15;
            //additionalIFD = null;
            //if (additionalIFD != null) // 261124Cl 変更前: 1 フレームあたりの追加タグの数を固定 (additionalIFD.Length) で数えていたので、
            //    totalIFD += (short)additionalIFD.Length; //  あるフレームだけタグを省く (PulsePower が NaN のフレーム) ことができなかった (IPAnalyzer 所見 8 の残り)
            // 261124Cl 変更: このフレームで実際に書くタグだけを数える (既知のタグで、値が空でも NaN でもないもの)。全フレームに値がある従来の呼び方は同じ数
            // 261124Cl (codex2 の査読): 値の型を確かめる (float・int なども double に変換できるもの。名前は空でない文字列)
            static bool writable(IFD f) => f.Data is { Length: > 0 } && f.Tag switch
            {
                60000 or 60002 => f.Data[0] is IConvertible c && c is not string && double.IsFinite(Convert.ToDouble(c, System.Globalization.CultureInfo.InvariantCulture)),
                60001 => f.Data[0] is string { Length: > 0 },
                _ => false
            };
            if (additionalIFD != null)
                for (int i = 0; i < additionalIFD.Length; i++)
                    if (writable(additionalIFD[i][n])) totalIFD++;

            bw.Write(BitConverter.GetBytes(totalIFD));
            offset = (uint)(bw.BaseStream.Position + totalIFD * 12 + 4);

            IFD.Write(bw, 254, (ushort)0);                              //NewSubfileType
            IFD.Write(bw, 256, (uint)imageWidth);                       //ImageWidth
            IFD.Write(bw, 257, (uint)(srcData[n].Length / imageWidth));    //ImageLength
            IFD.Write(bw, 258, (ushort)(bytePerPixel * 8));             //BitsPerSample
            IFD.Write(bw, 259, (ushort)1);                              //Compression
            IFD.Write(bw, 262, (ushort)1);                              //PhotometricInterpretation
                                                                        //offset = IFD.Write(bw, (ushort)270, "ImageJ=1.44p\nunit=" + unit, offset);         //ImageDescription
            IFD.Write(bw, 273, (uint)stripOffsets);                      //StripOffsets
            IFD.Write(bw, 278, (uint)(srcData[n].Length / imageWidth));    //RowsPerStrip
            IFD.Write(bw, 279, (uint)(srcData[n].Length * bytePerPixel));  //StripByteCounts
            IFD.Write(bw, 280, (ushort)0);                              //MinSampleValue
            // IFD.Write(bw, 281, (ushort)255);                         //MaxSampleValue // 260907Cl 変更前: bit 深度に関わらず 255 固定だった
            //260907Cl 変更: 16/32bit 整数で書いたときも 255 と名乗っていたため、MaxSampleValue を実際の bit 深度に合わせる。
            //  タグ 281 の型は SHORT なので 65535 が上限 (32bit 整数のときはそこで頭打ちにする)。
            //  値そのものは advisory (TIFF 6.0 で deprecated) で画素データの解釈には使われないが、
            //  これを表示レンジの初期値に使う読み手 (ImageJ 等) では 16bit 画像が飽和して見える。
            //  8bit (bytePerPixel==1) と float (sampleFormat==3) は従来どおり 255 のまま。
            IFD.Write(bw, 281, (ushort)(sampleFormat == 1 && bytePerPixel >= 2 ? 65535 : 255)); //MaxSampleValue
            offset = IFD.Write(bw, 282, (uint)1, (uint)pixelSizeX, offset); //XResolution
            offset = IFD.Write(bw, 283, (uint)1, (uint)pixelSizeY, offset); //YResolution
            IFD.Write(bw, 296, (ushort)1);                              //ResolutionUnit
            IFD.Write(bw, 339, (ushort)sampleFormat);                    //SampleFormat, 1: 整数, 3: float

            if (additionalIFD != null)
            {
                for (int i = 0; i < additionalIFD.Length; i++)
                {
                    if (!writable(additionalIFD[i][n])) continue; // 261124Cl 追加: 数えなかったタグは書かない
                    if (additionalIFD[i][n].Tag == 60000)//XrayEnergy
                        //offset = IFD.Write(bw, (ushort)additionalIFD[i][n].Tag, (double)additionalIFD[i][n].Data[0], offset); // 261124Cl 変更前: double 以外の数値で例外
                        offset = IFD.Write(bw, (ushort)additionalIFD[i][n].Tag, Convert.ToDouble(additionalIFD[i][n].Data[0], System.Globalization.CultureInfo.InvariantCulture), offset); // 261124Cl
                    if (additionalIFD[i][n].Tag == 60001)//Name
                        offset = IFD.Write(bw, (ushort)additionalIFD[i][n].Tag, (string)additionalIFD[i][n].Data[0], offset);
                    if (additionalIFD[i][n].Tag == 60002)//PulsePower
                        //offset = IFD.Write(bw, (ushort)additionalIFD[i][n].Tag, (double)additionalIFD[i][n].Data[0], offset); // 261124Cl 変更前: double 以外の数値で例外
                        offset = IFD.Write(bw, (ushort)additionalIFD[i][n].Tag, Convert.ToDouble(additionalIFD[i][n].Data[0], System.Globalization.CultureInfo.InvariantCulture), offset); // 261124Cl
                }
            }
        }

        bw.Close();
    }

    public class Loader
    {
        /// <summary>バイトオーダー</summary>
        private TiffByteOrder byteOrder;

        /// <summary>バージョン番号</summary>
        public int Version;

        /// <summary>格納しているイメージの数</summary>
        public int NumberOfFrames => Images.Count;

        public TiffByteOrder ByteOrder { get => byteOrder; set => byteOrder = value; }

        public List<imageProperty> Images = [];

        /// <summary>イメージの幅(ピクセル)</summary>
        public int ImageWidth;

        /// <summary>イメージの高さ(ピクセル)</summary>
        public int ImageLength;

        public bool IsGray;

        public int BitsPerSampleGray;

        public class imageProperty
        {
            public int BitsPerSampleGray;
            public int BitsPerSampleRed;
            public int BitsPerSampleGreen;
            public int BitsPerSampleBlue;

            /// <summary>1: 非圧縮</summary>
            public int Compression;

            /// <summary>261124Cl 追加: 圧縮前の差分予測 (タグ 317)。1: なし、2: 水平差分、3: 浮動小数点の差分</summary>
            public int Predictor = 1;

            /// <summary>0: 黒モードモノクロ (0が白い),   1: 白モードモノクロ (0が黒い),  2: RGBダイレクトカラー (0が黒い)</summary>
            public int PhotometricInterpretation;

            /// <summary>イメージの幅(ピクセル)</summary>
            public int ImageWidth;

            /// <summary>イメージの高さ(ピクセル)</summary>
            public int ImageLength;

            public double ResolutionX;
            public double ResolutionY;

            /// <summary>1: 整数, 3:float</summary>
            public int SampleFormat = 1;

            /// <summary>1: 単位なし, 2:inch,  3:cm</summary>
            public int ResolutionUnit = 2;

            public int[] ColorMapRed;
            public int[] ColorMapGreen;
            public int[] ColorMapBlue;

            public bool IsGray;

            public int[] StripOffsets;
            public int[] RowsPerStrip;
            public int[] StripByteCounts;

            //public uint[] ValueGray;
            public uint[] ValueRed;

            public uint[] ValueGreen;
            public uint[] ValueBlue;
            public double[] Value;

            public string ImageDescription = "";

            //gelファイルのタグ
            public double MDScaleFactor = 1;//強度に乗ずるスケール

            public int MDFileTag = 128;

            //h5ファイルのタグ
            public double XrayEnergy = double.NaN;

            public string Name = "";
            public double PulsePower = double.NaN;

            //TIAというソフトからエクスポートされるTEM写真(emiファイル)のピクセルサイズ
            public double TIA_PixelSizeX;//65450
            public double TIA_PixelSizeY;//65451
            public int TIA_PixelSizeUnit;//65452
        }

        public static byte[] Read(BinaryReader br, long position, int length, TiffByteOrder byteOrder)
        {
            br.BaseStream.Position = position;
            var buffer = new byte[length];
            _ = br.Read(buffer, 0, buffer.Length);
            if (byteOrder == TiffByteOrder.Intel) buffer = buffer.Reverse().ToArray();
            return buffer;
        }

        public Loader(string fileName)
        {
            //MDFileTag = 128;
            //MDScaleFactor = 1;
            bool originalEndian = BitConverter.IsLittleEndian;

            //var br = new BinaryReader(new FileStream(fileName, FileMode.Open, FileAccess.Read));// (260715Ch) 旧: 不正ヘッダーの早期 return / 例外時にファイルを解放しない
            using var br = new BinaryReader(new FileStream(fileName, FileMode.Open, FileAccess.Read));// (260715Ch)

            //まず最初の2Byteを読み込んでバイトオーダーを決める
            var temp = new byte[2];
            br.Read(temp, 0, 2);
            if (temp[0] == 0x49 && temp[1] == 0x49)
                ByteOrder = TiffByteOrder.Motorola;
            else if (temp[0] == 0x4D && temp[1] == 0x4D)
                ByteOrder = TiffByteOrder.Intel;
            else
                return;

            //次に2Byteよんでバージョン
            Version = BitConverter.ToInt16(Read(br, 2, 2, ByteOrder), 0);

            var ifdPointa = long.MaxValue;

            while (ifdPointa != 0)
            {
                var image = new imageProperty();
                //ifdPointaがlong.MaxValueの時(何も読み込んでいないとき) 4Byte読んでIFDのポインタ
                if (ifdPointa == long.MaxValue)
                    ifdPointa = BitConverter.ToUInt32(Read(br, 4, 4, ByteOrder), 0);

                //ポインタ位置から2ByteよんでIFDの総数
                int totalIFD = BitConverter.ToInt16(Read(br, ifdPointa, 2, ByteOrder), 0);

                //IFDの読み込み開始
                ifdPointa += 2;
                var iFD = new IFD[totalIFD];
                for (int i = 0; i < totalIFD; i++)
                {
                    try { iFD[i] = ReadIFD(br, (int)ifdPointa, ByteOrder); }
                    catch { }
                    ifdPointa += 12;
                }
                //ここまででIFDデータを読み込み完了

                //次のifdPointaを読み込む

                ifdPointa = BitConverter.ToUInt32(Read(br, ifdPointa, 4, ByteOrder), 0);

                int bitsPerSampleLength = 1;

                #region  IFDのデータ解析

                for (int i = 0; i < totalIFD; i++)
                    switch (iFD[i].Tag)
                    {
                        case 256:
                            image.ImageWidth = (int)iFD[i].Data[0]; break;
                        case 257:
                            image.ImageLength = (int)iFD[i].Data[0]; break;
                        case 258:

                            if (iFD[i].Data.Length == 1)
                            {
                                image.IsGray = true;
                                image.BitsPerSampleGray = (int)iFD[i].Data[0];
                            }
                            else if (iFD[i].Data.Length == 3 || iFD[i].Data.Length == 4)
                            {
                                bitsPerSampleLength = iFD[i].Data.Length;
                                image.IsGray = false;
                                image.BitsPerSampleRed = (int)iFD[i].Data[0];
                                image.BitsPerSampleGreen = (int)iFD[i].Data[1];
                                image.BitsPerSampleBlue = (int)iFD[i].Data[2];
                            }
                            break;

                        case 259:
                            image.Compression = (int)iFD[i].Data[0]; break;
                        case 262:
                            // 261112Cl 注記: 0 (WhiteIsZero) でも値は反転しない (作者判断 261112)。FLA/Typhoon の MD gel・Rad-icon・PF の IP などの検出器の画像が
                            //  この印付きで、保存された値がそのまま強度 (大きいほど強い) になっている。反転すると強度が逆になる。
                            image.PhotometricInterpretation = (int)iFD[i].Data[0]; break;
                        case 270:
                            break;
                        case 273:
                            image.StripOffsets = new int[iFD[i].Data.Length];
                            for (int j = 0; j < image.StripOffsets.Length; j++)
                                image.StripOffsets[j] = (int)iFD[i].Data[j];
                            break;

                        case 278:
                            image.RowsPerStrip = new int[iFD[i].Data.Length];
                            for (int j = 0; j < image.RowsPerStrip.Length; j++)
                                image.RowsPerStrip[j] = (int)iFD[i].Data[j];
                            break;

                        case 279:
                            image.StripByteCounts = new int[iFD[i].Data.Length];
                            for (int j = 0; j < image.StripByteCounts.Length; j++)
                                image.StripByteCounts[j] = (int)iFD[i].Data[j];
                            break;

                        case 282:
                            if (iFD[i].Data[0] is int n)
                                image.ResolutionX = n;
                            else if (iFD[i].Data[0] is float f)
                                image.ResolutionX = f;
                            else if (iFD[i].Data[0] is double d)
                                image.ResolutionX = d;
                            break;
                        case 283:
                            if (iFD[i].Data[0] is int n2)
                                image.ResolutionY = n2;
                            else if (iFD[i].Data[0] is float f)
                                image.ResolutionY = f;
                            else if (iFD[i].Data[0] is double d)
                                image.ResolutionY = d;
                            break;
                        case 284:
                            image.ResolutionUnit = (int)iFD[i].Data[0]; break;
                        case 320:
                            // 261112Cl 追加: 色表の配列を確保する (以前は確保していなかったので、下の if が常に偽で色表を読んでいなかった)
                            image.ColorMapRed = new int[iFD[i].Data.Length / 3];
                            image.ColorMapGreen = new int[iFD[i].Data.Length / 3];
                            image.ColorMapBlue = new int[iFD[i].Data.Length / 3];
                            for (int j = 0; j < iFD[i].Data.Length / 3; j++)
                            {
                                if (image.ColorMapRed != null)
                                {
                                    image.ColorMapRed[j] = (int)iFD[i].Data[j];
                                    image.ColorMapGreen[j] = (int)iFD[i].Data[j + iFD[i].Data.Length / 3];
                                    image.ColorMapBlue[j] = (int)iFD[i].Data[j + 2 * iFD[i].Data.Length / 3];
                                }
                            }
                            break;

                        case 317:
                            image.Predictor = (int)iFD[i].Data[0]; break; // 261124Cl 追加
                        case 339:
                            image.SampleFormat = (int)iFD[i].Data[0]; break;  //SampleFormat, 1: 整数, 3: float

                        case 33449:
                            {
                                var sb = new StringBuilder();
                                for (int j = 0; j < iFD[i].Data.Length; j++)
                                    //if ((char)iFD[i].Data[j] != '\r\n')
                                    sb.Append((char)iFD[i].Data[j]);
                                image.ImageDescription = sb.ToString();
                            }
                            break;

                        case 33446:
                            image.MDScaleFactor = (double)iFD[i].Data[0]; break;
                        case 33445:
                            image.MDFileTag = (int)iFD[i].Data[0]; break;

                        case 60000:
                            image.XrayEnergy = (double)iFD[i].Data[0]; break;
                        case 60001:
                            image.Name = new string(iFD[i].Data.Cast<char>().ToArray());
                            if (image.Name.EndsWith('\0'))
                                image.Name = image.Name.TrimEnd('\0');
                            break;

                        case 60002:
                            image.PulsePower = (double)iFD[i].Data[0]; break;

                        case 65450:
                            image.TIA_PixelSizeX = (double)iFD[i].Data[0]; break;
                        case 65451:
                            image.TIA_PixelSizeY = (double)iFD[i].Data[0]; break;
                        case 65452:
                            image.TIA_PixelSizeUnit = (int)iFD[i].Data[0]; break;
                    }
                #endregion

                //画像のサイズやカラーモードが違ったら、読み込み停止
                if (Images.Count == 0 || (ImageLength == image.ImageLength && ImageWidth == image.ImageWidth && image.IsGray == IsGray && BitsPerSampleGray == image.BitsPerSampleGray))
                {
                    ImageLength = image.ImageLength;
                    ImageWidth = image.ImageWidth;
                    IsGray = image.IsGray;
                    BitsPerSampleGray = image.BitsPerSampleGray;
                }
                else
                    break;

                int bytePerPixel = 1;
                if (image.IsGray || image.SampleFormat == 3)//grayか、浮動小数点の時
                {
                    bytePerPixel = image.BitsPerSampleGray / 8;
                    if (bytePerPixel == 0) // 261124Cl 追加 (査読): 1・4 bit の画像は読めない。以前は下の割り算で 0 除算になっていた
                        throw new NotSupportedException($"{image.BitsPerSampleGray}-bit TIFF images are not supported.");
                    image.Value = new double[image.ImageLength * image.ImageWidth];
                }
                else
                {
                    bytePerPixel = image.BitsPerSampleRed / 8;
                    image.ValueRed = new uint[image.ImageLength * image.ImageWidth];
                    image.ValueGreen = new uint[image.ImageLength * image.ImageWidth];
                    image.ValueBlue = new uint[image.ImageLength * image.ImageWidth];
                }

                //ここからデータを読み込む

                #region PFで導入された新フォーマットへの対処
                if (iFD.Any(i => i.Tag == 270) && image.StripByteCounts == null)
                {
                    var ifd = iFD.First(i => i.Tag == 270);
                    var comment = new string(ifd.Data.Cast<char>().ToArray());

                    //software
                    var software = comment[..16];
                    //Sensor
                    string sensor = comment.Substring(16, 20);
                    //Image Info
                    var correction = comment[38];
                    var type = comment[39];
                    var binning_x = comment[40];
                    var binning_y = comment[41];
                    var direction = comment[42];
                    var bitLength = comment[43];

                    //Accumulation number　画像取得時の積算回数	
                    var accumulation_number = comment.Substring(44, 4).ToInt();

                    //ExposureTime
                    var exposure_time = comment.Substring(48, 8).ToInt();

                    var frame = comment.Substring(56, 4).ToInt();

                    //Number of stored images
                    var num_of_stored_images = comment.Substring(60, 4).ToInt();

                    var total = image.ImageWidth * image.ImageLength;
                    image.StripByteCounts = [image.ImageWidth * image.ImageLength * bytePerPixel];
                }
                #endregion


                //261111Cl 追加 (T33): 圧縮された画素の展開には対応していない。誤読した画像を黙って返さず、明確なエラーにする (タグが無い = 0 は非圧縮として扱う)
                //if (image.Compression != 0 && image.Compression != 1) // 261124Cl 変更前
                //    throw new NotSupportedException($"Compressed TIFF (Compression = {image.Compression}) is not supported. Please save the image without compression.");
                // 261124Cl 変更 (作者指示: Spot ID v2 で圧縮 TIFF を読む): LZW (5)・Deflate (8, 32946)・PackBits (32773) は展開して読む。
                //  差分予測は水平差分 (Predictor 2) まで。それ以外の圧縮 (JPEG など) と浮動小数点の差分 (Predictor 3) は従来どおり明確なエラーにする
                bool compressed = image.Compression is not (0 or 1);
                if (compressed && image.Compression is not (5 or 8 or 32946 or 32773))
                    throw new NotSupportedException($"This TIFF compression (Compression = {image.Compression}) is not supported. Uncompressed, LZW, Deflate (zip) and PackBits TIFFs can be read.");
                if (image.Predictor is not (1 or 2))
                    throw new NotSupportedException($"This TIFF predictor (Predictor = {image.Predictor}) is not supported. Save without a predictor or with the horizontal predictor.");

                #region 261111Cl 変更前: 画素ごとに Func<double> を呼ぶ (MM の 16/32bit と浮動小数点は、画素ごとに byte[] を確保して反転していた)
                //Func<double> toInt = () => 0;
                //if (bytePerPixel == 1)
                //    toInt = () => br.ReadByte();
                //else if (byteOrder == TiffByteOrder.Motorola)
                //{
                //
                //    if (bytePerPixel == 2 && image.SampleFormat == 1)
                //        toInt = () => br.ReadUInt16();
                //    else if (bytePerPixel == 2 && image.SampleFormat == 2)
                //        toInt = () => br.ReadInt16();
                //    else if (bytePerPixel == 4 && image.SampleFormat == 1)
                //        toInt = () => br.ReadUInt32();
                //    else if (bytePerPixel == 4 && image.SampleFormat == 2)
                //        toInt = () => br.ReadInt32();
                //}
                //else
                //{
                //    if (bytePerPixel == 2)
                //    {
                //        if (image.SampleFormat == 1)
                //            toInt = () =>
                //            {
                //                var temp = new byte[2];
                //                br.Read(temp, 0, 2);
                //                temp = temp.Reverse().ToArray();
                //                return BitConverter.ToUInt16(temp, 0);
                //            };
                //
                //        else
                //            toInt = () =>
                //            {
                //                var temp = new byte[2];
                //                br.Read(temp, 0, 2);
                //                temp = temp.Reverse().ToArray();
                //                return BitConverter.ToInt16(temp, 0);
                //            };
                //    }
                //    else if (bytePerPixel == 4)
                //    {
                //        if (image.SampleFormat == 1)
                //            toInt = () =>
                //            {
                //                var temp = new byte[4];
                //                br.Read(temp, 0, 4);
                //                temp = temp.Reverse().ToArray();
                //                return BitConverter.ToUInt32(temp, 0);
                //            };
                //
                //        else
                //            toInt = () =>
                //            {
                //                var temp = new byte[4];
                //                br.Read(temp, 0, 4);
                //                temp = temp.Reverse().ToArray();
                //                return BitConverter.ToInt32(temp, 0);
                //            };
                //    }
                //}
                //
                //
                //
                //for (int i = 0, n = 0; i < image.StripByteCounts.Length; i++)
                //{
                //    br.BaseStream.Position = image.StripOffsets[i];
                //    if (image.SampleFormat == 3)//浮動小数点データの時
                //    {
                //        for (int j = 0; j < image.StripByteCounts[i] / bytePerPixel; j++)
                //            image.Value[n++] = toFloat(br, bytePerPixel, ByteOrder);
                //    }
                //    else//整数データの時 1の時は符号なし。2の時は符号アリ
                //    {
                //        var sign = image.SampleFormat == 2;
                //        if (image.IsGray)
                //        {
                //
                //            if (image.MDFileTag == 2)//gelファイルの時
                //            {
                //                for (int j = 0; j < image.StripByteCounts[i] / bytePerPixel; j++)
                //                {
                //                    var intensity = toInt();
                //                    image.Value[n++] = intensity * intensity;
                //                }
                //            }
                //            else
                //            {
                //                for (int j = 0; j < image.StripByteCounts[i] / bytePerPixel; j++)
                //                    image.Value[n++] = toInt();
                //            }
                //        }
                //        else
                //        {
                //            for (int j = 0; j < image.StripByteCounts[i] / bytePerPixel / bitsPerSampleLength; j++)
                //            {
                //                image.ValueRed[n] = (uint)toInt();
                //                image.ValueGreen[n] = (uint)toInt();
                //                image.ValueBlue[n] = (uint)toInt();
                //                if (bitsPerSampleLength == 4)
                //                    toInt();
                //                n++;
                //            }
                //        }
                //    }
                //}
                #endregion

                //261111Cl 変更 (T1): ストリップ単位で一括で読む (ImageIO.ReadPixels: MemoryMarshal で型を変え、MM ならバイト順を反転してから double に広げる)。
                //  型の選び方 (bytePerPixel と SampleFormat の組) は旧式と同じ。ファイル末尾での振る舞いも旧式と同じにする:
                //  BinaryReader.ReadXxx で読んでいた 8bit と II の 16/32bit 整数は末尾で例外、BinaryReader.Read(byte[]) で読んでいた MM の 16/32bit 整数と浮動小数点は足りない分を 0。
                //  F2 (挙動が変わる): 浮動小数点は旧式では常に 4 バイトずつ読んでいたため、64bit (double) と 16bit (half) の TIFF を誤読していた。幅に合わせて読む。
                var bigEndian = byteOrder == TiffByteOrder.Intel;//列挙子の名前が逆: Intel = MM (ビッグエンディアン)、Motorola = II
                //var kind = image.SampleFormat == 3 ? bytePerPixel switch { 8 => SampleKind.F64, 2 => SampleKind.F16, _ => SampleKind.F32 } // 261112Cl 変更前
                //    : bytePerPixel == 1 ? SampleKind.U8
                //    : !bigEndian ? (bytePerPixel, image.SampleFormat) switch { (2, 1) => SampleKind.U16, (2, 2) => SampleKind.S16, (4, 1) => SampleKind.U32, (4, 2) => SampleKind.S32, _ => SampleKind.None }
                //    : bytePerPixel == 2 ? (image.SampleFormat == 1 ? SampleKind.U16 : SampleKind.S16)
                //    : bytePerPixel == 4 ? (image.SampleFormat == 1 ? SampleKind.U32 : SampleKind.S32)
                //    : SampleKind.None;
                // 261112Cl 変更: 64 bit 整数 (符号なし・符号あり) も読む (以前は None で全画素 0 になっていた)。ほかの組は変更前と同じ
                var kind = image.SampleFormat == 3 ? bytePerPixel switch { 8 => SampleKind.F64, 2 => SampleKind.F16, _ => SampleKind.F32 }
                    : bytePerPixel == 1 ? SampleKind.U8
                    : !bigEndian ? (bytePerPixel, image.SampleFormat) switch { (2, 1) => SampleKind.U16, (2, 2) => SampleKind.S16, (4, 1) => SampleKind.U32, (4, 2) => SampleKind.S32, (8, 1) => SampleKind.U64, (8, 2) => SampleKind.S64, _ => SampleKind.None }
                    : bytePerPixel == 2 ? (image.SampleFormat == 1 ? SampleKind.U16 : SampleKind.S16)
                    : bytePerPixel == 4 ? (image.SampleFormat == 1 ? SampleKind.U32 : SampleKind.S32)
                    : bytePerPixel == 8 ? (image.SampleFormat == 1 ? SampleKind.U64 : SampleKind.S64)
                    : SampleKind.None;

                // 261124Cl 追加: 圧縮 (または差分予測) のあるストリップは、生のバイト列に展開してから同じ読み方をする
                bool buffered = compressed || image.Predictor == 2;
                int samplesPerPixel = image.SampleFormat == 3 || image.IsGray ? 1 : bitsPerSampleLength;
                int rowBytes = image.ImageWidth * samplesPerPixel * bytePerPixel;
                int stripRows = image.RowsPerStrip is { Length: > 0 } rps && rps[0] > 0 ? Math.Min(rps[0], image.ImageLength) : image.ImageLength;
                for (int i = 0, n = 0; i < image.StripByteCounts.Length; i++)
                {
                    br.BaseStream.Position = image.StripOffsets[i];
                    Stream src = br.BaseStream; // 261124Cl 追加
                    int stripBytes = image.StripByteCounts[i]; // 261124Cl 追加
                    if (buffered) // 261124Cl 追加
                    {
                        // 展開後の大きさ = この帯の行数 × 1 行のバイト数 (最後の帯は残りの行数)。int を超える大きさは扱えない
                        long expected = (long)rowBytes * Math.Max(1, Math.Min(stripRows, image.ImageLength - i * stripRows));
                        if (expected > Array.MaxLength - 4096)
                            throw new NotSupportedException("This compressed TIFF strip is too large to decode.");
                        var dec = DecodeStrip(br.ReadBytes(stripBytes), image.Compression, (int)expected);
                        if (image.Predictor == 2)
                            UndoHorizontalPredictor(dec, image.ImageWidth, samplesPerPixel, bytePerPixel, bigEndian);
                        src = new MemoryStream(dec, false);
                        stripBytes = dec.Length;
                    }
                    if (image.SampleFormat == 3 || image.IsGray)//浮動小数点データ、またはグレーの整数データの時 (整数は 1 の時は符号なし。2 の時は符号アリ)
                    {
                        //int count = Math.Max(0, image.StripByteCounts[i] / bytePerPixel); // 261124Cl 変更前
                        int count = Math.Max(0, stripBytes / bytePerPixel); // 261124Cl 変更: 展開した後のバイト数で数える
                        int fit = Math.Min(count, image.Value.Length - n);
                        var dst = image.Value.AsSpan(n, fit);
                        ReadSamples(src, dst, kind, bigEndian); // 261124Cl 変更 (旧 br.BaseStream)
                        if (image.SampleFormat != 3 && image.MDFileTag == 2)//gelファイルの時
                            for (int j = 0; j < dst.Length; j++)
                                dst[j] *= dst[j];
                        n += fit;
                        if (fit < count)//旧式と同じく、配列を超える分があれば IndexOutOfRangeException
                            throw new IndexOutOfRangeException();
                    }
                    else//RGB(A) の整数データ: 画素ごとに bitsPerSampleLength 個の標本。一時配列に読んで振り分ける
                    {
                        //int count = Math.Max(0, image.StripByteCounts[i] / bytePerPixel / bitsPerSampleLength); // 261124Cl 変更前
                        int count = Math.Max(0, stripBytes / bytePerPixel / bitsPerSampleLength); // 261124Cl 変更
                        int fit = Math.Min(count, image.ValueRed.Length - n);
                        const int chunk = 1 << 14;
                        var tmp = ArrayPool<double>.Shared.Rent(Math.Min(fit, chunk) * bitsPerSampleLength);
                        try
                        {
                            for (int p = 0; p < fit; p += chunk)
                            {
                                int m = Math.Min(chunk, fit - p);
                                var t = tmp.AsSpan(0, m * bitsPerSampleLength);
                                ReadSamples(src, t, kind, bigEndian); // 261124Cl 変更 (旧 br.BaseStream)
                                for (int k = 0, q = n + p; k < m; k++, q++)
                                {
                                    image.ValueRed[q] = (uint)t[k * bitsPerSampleLength];
                                    image.ValueGreen[q] = (uint)t[k * bitsPerSampleLength + 1];
                                    image.ValueBlue[q] = (uint)t[k * bitsPerSampleLength + 2];
                                }
                            }
                        }
                        finally { ArrayPool<double>.Shared.Return(tmp); }
                        n += fit;
                        if (fit < count)
                            throw new IndexOutOfRangeException();
                    }
                }

                // 261112Cl 追加: パレット画像 (PhotometricInterpretation = 3) は、画素値 (色表の添字) を色表の色の輝度に置き換える (以前は添字をそのまま強度にしていた)。
                //  輝度は BT.601 の重み (0.299, 0.587, 0.114) で、目盛りは 0〜255 (色表の値 / 257)。灰色の色表 (i 番目が i×257) なら値は変わらない。
                if (image.PhotometricInterpretation == 3 && image.ColorMapRed is { Length: > 0 } && image.Value != null)
                {
                    var lum = new double[image.ColorMapRed.Length];
                    for (int k = 0; k < lum.Length; k++)
                    {
                        int r = image.ColorMapRed[k], g = image.ColorMapGreen[k], b = image.ColorMapBlue[k];
                        lum[k] = (r == g && g == b ? r : 0.299 * r + 0.587 * g + 0.114 * b) / 257.0;
                    }
                    var v = image.Value;
                    for (int k = 0; k < v.Length; k++)
                        if (v[k] >= 0 && v[k] < lum.Length)
                            v[k] = lum[(int)v[k]];
                }
                Images.Add(image);
            }
            br.Close();
        }


        //261111Cl 変更前: 浮動小数点の 1 画素読み (常に 4 バイト・画素ごとに byte[4] を確保)。ReadSamples に置き換えたので使っていない
        //private static float toFloat(BinaryReader br, int byteCount, TiffByteOrder byteOrder)
        //{
        //    float intensity = 0;
        //    var temp = new byte[4];
        //    br.Read(temp, 0, 4);
        //
        //    if (byteOrder == TiffByteOrder.Intel) temp = temp.Reverse().ToArray();
        //
        //    intensity = BitConverter.ToSingle(temp, 0);
        //
        //    return intensity;
        //}

        /// <summary>標本の型 (一括読み込み用)。261111Cl 追加</summary>
        //private enum SampleKind { None, U8, U16, S16, U32, S32, F16, F32, F64 } // 261112Cl 変更前
        private enum SampleKind { None, U8, U16, S16, U32, S32, F16, F32, F64, U64, S64 } // 261112Cl 変更: 64 bit 整数を追加

        /// <summary>261124Cl 追加: 圧縮されたストリップを展開する。5: LZW、8 / 32946: Deflate (zlib)、32773: PackBits。
        /// expected はこの帯の展開後の大きさ。結果はそこで打ち切る (余分な展開で次の帯の画素がずれる・際限なくメモリを使うのを防ぐ)</summary>
        private static byte[] DecodeStrip(byte[] src, int compression, int expected) => compression switch
        {
            5 => DecodeLzw(src, expected),
            8 or 32946 => DecodeDeflate(src, expected),
            32773 => DecodePackBits(src, expected),
            _ => src,
        };

        /// <summary>261124Cl 追加: TIFF の LZW (MSB から詰めた 9〜12 ビットの符号、符号の幅は 1 つ早く増える "early change")</summary>
        private static byte[] DecodeLzw(byte[] src, int expected)
        {
            const int clear = 256, eoi = 257;
            var prefix = new short[4096];
            var suffix = new byte[4096];
            var first = new byte[4096];
            var length = new int[4096];
            for (int k = 0; k < 256; k++) { prefix[k] = -1; suffix[k] = first[k] = (byte)k; length[k] = 1; }
            var dst = new byte[expected + 4096]; // 1 つの符号の展開は最長 4096 バイト。expected を超えた分は最後に捨てる
            int o = 0, next = 258, width = 9, prev = -1, bits = 0, pos = 0;
            uint buffer = 0;

            void Write(int code)
            {
                int len = length[code];
                for (int c = code, q = o + len - 1; q >= o; c = prefix[c], q--)
                    dst[q] = suffix[c];
                o += len;
            }
            void Add(int pre, byte last)
            {
                if (next >= 4096) return;
                prefix[next] = (short)pre; suffix[next] = last; first[next] = first[pre]; length[next] = length[pre] + 1;
                next++;
                if (next + 1 == 1 << width && width < 12) width++; // early change: 511・1023・2047 で幅を増やす
            }

            while (o < expected)
            {
                while (bits < width && pos < src.Length) { buffer = (buffer << 8) | src[pos++]; bits += 8; }
                if (bits < width) break;
                int code = (int)((buffer >> (bits - width)) & ((1u << width) - 1));
                bits -= width;
                if (code == eoi) break;
                if (code == clear) { next = 258; width = 9; prev = -1; continue; }
                if (prev == -1) { if (code >= 256) break; Write(code); prev = code; continue; }
                if (code < next) { Write(code); Add(prev, first[code]); }
                else if (code == next) { Add(prev, first[prev]); Write(code); }
                else break; // 壊れたデータ
                prev = code;
            }
            Array.Resize(ref dst, Math.Min(o, expected));
            return dst;
        }

        /// <summary>261124Cl 追加: Deflate (zlib の枠つき)。帯の大きさ (expected) の配列へ直接読み、それを超える分は読まない</summary>
        private static byte[] DecodeDeflate(byte[] src, int expected)
        {
            using var z = new System.IO.Compression.ZLibStream(new MemoryStream(src, false), System.IO.Compression.CompressionMode.Decompress);
            var dst = new byte[expected];
            int o = 0, n;
            while (o < dst.Length && (n = z.Read(dst, o, dst.Length - o)) > 0)
                o += n;
            if (o != dst.Length) Array.Resize(ref dst, o);
            return dst;
        }

        /// <summary>261124Cl 追加: PackBits (0〜127: 続く n+1 バイトをそのまま、-1〜-127: 次の 1 バイトを 1−n 回、-128: 何もしない)</summary>
        private static byte[] DecodePackBits(byte[] src, int expected)
        {
            var dst = new byte[expected];
            int o = 0;
            for (int p = 0; p < src.Length && o < dst.Length;)
            {
                int n = (sbyte)src[p++], m;
                if (n >= 0) m = Math.Min(n + 1, src.Length - p);
                else if (n != -128 && p < src.Length) m = 1 - n;
                else continue;
                int w = Math.Min(m, dst.Length - o); // 帯の大きさで打ち切る
                if (n >= 0) { Array.Copy(src, p, dst, o, w); p += m; }
                else dst.AsSpan(o, w).Fill(src[p++]);
                o += w;
            }
            if (o != dst.Length) Array.Resize(ref dst, o);
            return dst;
        }

        /// <summary>261124Cl 追加: 水平差分の予測 (Predictor 2) を戻す。行ごとに、各標本へ同じ成分の左隣の値を足す。
        /// 標本の幅 (1・2・4・8 バイト) によらず、下位のバイトから繰り上がりつきで 1 バイトずつ足す (バイト順はファイルのとおり)</summary>
        private static void UndoHorizontalPredictor(byte[] d, int width, int samplesPerPixel, int bytesPerSample, bool bigEndian)
        {
            int rowSamples = width * samplesPerPixel, rowBytes = rowSamples * bytesPerSample, step = samplesPerPixel * bytesPerSample;
            if (rowBytes <= 0) return;
            for (int r = 0; r + rowBytes <= d.Length; r += rowBytes)
                for (int a = r + step; a < r + rowBytes; a += bytesPerSample)
                {
                    int carry = 0;
                    for (int j = 0; j < bytesPerSample; j++)
                    {
                        int k = bigEndian ? bytesPerSample - 1 - j : j;
                        int sum = d[a + k] + d[a - step + k] + carry;
                        d[a + k] = (byte)sum;
                        carry = sum >> 8;
                    }
                }
        }

        /// <summary>stream の現在位置から dst.Length 個の標本を読み、double にして dst に書く。261111Cl 追加
        /// 8bit と II の整数はファイル末尾で EndOfStreamException (旧 BinaryReader.ReadXxx と同じ)、
        /// MM の整数と浮動小数点は足りない分を 0 とする (旧 BinaryReader.Read(byte[], ..) と同じ)。None は読まずに 0 (旧 toInt の既定 () => 0 と同じ)。</summary>
        private static void ReadSamples(Stream s, Span<double> dst, SampleKind kind, bool bigEndian)
        {
            switch (kind)
            {
                case SampleKind.U8: ImageIO.ReadPixels<byte>(s, dst); break;
                case SampleKind.U16: ImageIO.ReadPixels<ushort>(s, dst, bigEndian, zeroFillAtEnd: bigEndian); break;
                case SampleKind.S16: ImageIO.ReadPixels<short>(s, dst, bigEndian, zeroFillAtEnd: bigEndian); break;
                case SampleKind.U32: ImageIO.ReadPixels<uint>(s, dst, bigEndian, zeroFillAtEnd: bigEndian); break;
                case SampleKind.S32: ImageIO.ReadPixels<int>(s, dst, bigEndian, zeroFillAtEnd: bigEndian); break;
                case SampleKind.F16: ImageIO.ReadPixels<Half>(s, dst, bigEndian, zeroFillAtEnd: true); break;
                case SampleKind.F32: ImageIO.ReadPixels<float>(s, dst, bigEndian, zeroFillAtEnd: true); break;
                case SampleKind.F64: ImageIO.ReadPixels<double>(s, dst, bigEndian, zeroFillAtEnd: true); break;
                case SampleKind.U64: ImageIO.ReadPixels<ulong>(s, dst, bigEndian, zeroFillAtEnd: bigEndian); break; // 261112Cl 追加
                case SampleKind.S64: ImageIO.ReadPixels<long>(s, dst, bigEndian, zeroFillAtEnd: bigEndian); break; // 261112Cl 追加
                default: dst.Clear(); break;
            }
        }

        public static IFD ReadIFD(BinaryReader br, int index, TiffByteOrder byteOrder)
        {
            br.BaseStream.Position = index;

            var ifd = new IFD();

            byte[] temp;

            //まずタグを読み込み
            temp = new byte[2];
            br.Read(temp, 0, 2);
            if (byteOrder == TiffByteOrder.Intel)
                temp = temp.Reverse().ToArray();
            ifd.Tag = BitConverter.ToUInt16(temp, 0);

            if (ifd.Tag == 0)
                return ifd;

            //次にデータタイプを読み込み
            temp = new byte[2];
            br.Read(temp, 0, 2);
            if (byteOrder == TiffByteOrder.Intel)
                temp = temp.Reverse().ToArray();
            int datatype = BitConverter.ToUInt16(temp, 0);

            int dataLengh = 1;
            switch (datatype)//データのバイト数を読み込み
            {
                case 1: dataLengh = 1; ifd.DataType = typeof(byte); break;
                case 2: dataLengh = 1; ifd.DataType = typeof(char); break;
                case 3: dataLengh = 2; ifd.DataType = typeof(ushort); break;
                case 4: dataLengh = 4; ifd.DataType = typeof(uint); break;
                case 5: dataLengh = 8; ifd.DataType = typeof(double); break;
                case 6: dataLengh = 1; ifd.DataType = typeof(sbyte); break;
                case 7: dataLengh = 1; ifd.DataType = typeof(byte); break;
                case 8: dataLengh = 2; ifd.DataType = typeof(short); break;
                case 9: dataLengh = 4; ifd.DataType = typeof(int); break;
                case 10: dataLengh = 8; ifd.DataType = typeof(double); break;
                case 11: dataLengh = 4; ifd.DataType = typeof(float); break;
                case 12: dataLengh = 8; ifd.DataType = typeof(double); break;
            }

            //次にデータカウントを読み込み
            temp = new byte[4];
            br.Read(temp, 0, 4);
            if (byteOrder == TiffByteOrder.Intel) temp = temp.Reverse().ToArray();
            int datacount = BitConverter.ToInt32(temp, 0);//データカウント読み込み

            temp = new byte[4];
            br.Read(temp, 0, 4);

            if (dataLengh * datacount > 4)//データがポインタを示しているときtempを設定しなおす
            {
                if (byteOrder == TiffByteOrder.Intel)
                    temp = temp.Reverse().ToArray();
                int pointa = BitConverter.ToInt32(temp, 0);
                temp = new byte[dataLengh * datacount];
                br.BaseStream.Position = pointa;
                br.Read(temp, 0, dataLengh * datacount);
            }

            //次にTempをDataCount数に分割する

            byte[][] data = new byte[datacount][];
            for (int i = 0; i < datacount; i++)
            {
                data[i] = new byte[dataLengh];
                Array.Copy(temp, i * dataLengh, data[i], 0, dataLengh);
                if (byteOrder == TiffByteOrder.Intel)
                    data[i] = data[i].Reverse().ToArray();
            }

            //次にデータフィールドあるいはデータポインタを読み込み

            ifd.Data = new object[datacount];

            for (int i = 0; i < datacount; i++)
            {
                switch (datatype)
                {
                    case 1: ifd.Data[i] = data[i][0]; break;
                    case 2: ifd.Data[i] = (char)data[i][0]; break;
                    case 3: ifd.Data[i] = (int)BitConverter.ToUInt16(data[i], 0); break;
                    case 4: ifd.Data[i] = (int)BitConverter.ToUInt32(data[i], 0); break;
                    case 5:
                        //double a1 = (double)(BitConverter.ToUInt32(data[i], 0));
                        //double a2 = (double)(BitConverter.ToUInt32(data[i], 4));
                        //261111Cl 変更 (F1, 挙動が変わる): MM では上で 8 バイトをまとめて反転しているので、分子は後半 4 バイト・分母は前半 4 バイトに来る。
                        //  旧式は MM で分母/分子 (逆数) を返していた。II は従来どおり
                        double a1 = (double)(BitConverter.ToUInt32(data[i], byteOrder == TiffByteOrder.Intel ? 4 : 0));
                        double a2 = (double)(BitConverter.ToUInt32(data[i], byteOrder == TiffByteOrder.Intel ? 0 : 4));
                        ifd.Data[i] = a1 / a2; break;
                    case 6: ifd.Data[i] = (sbyte)data[i][0]; break;
                    case 7: ifd.Data[i] = data[i][0]; break;
                    case 8: ifd.Data[i] = (int)BitConverter.ToInt16(data[i], 0); break;
                    case 9: ifd.Data[i] = (int)BitConverter.ToInt32(data[i], 0); break;
                    //case 10: ifd.Data[i] = (double)(BitConverter.ToInt32(data[i], 0) / BitConverter.ToInt32(data[i], 4)); break;// (260715Ch) 旧: SRATIONAL を整数除算して小数部を失っていた
                    //case 10: ifd.Data[i] = (double)BitConverter.ToInt32(data[i], 0) / BitConverter.ToInt32(data[i], 4); break;// (260715Ch)
                    case 10: ifd.Data[i] = (double)BitConverter.ToInt32(data[i], byteOrder == TiffByteOrder.Intel ? 4 : 0) / BitConverter.ToInt32(data[i], byteOrder == TiffByteOrder.Intel ? 0 : 4); break;// 261111Cl 変更 (F1): SRATIONAL も RATIONAL と同じく MM で分子と分母が入れ替わっていた
                    case 11: ifd.Data[i] = (double)BitConverter.ToSingle(data[i], 0); break;
                    case 12: ifd.Data[i] = (double)BitConverter.ToDouble(data[i], 0); break;
                }
            }
            return ifd;
        }
    }
}
