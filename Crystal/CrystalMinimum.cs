using MemoryPack;
using MemoryPack.Compression;
using System;
using System.Collections.Frozen;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
// using ZLinq; // 261111Cl 削除: ZLinq を外したため
using System.Text;
using System.Windows.Forms;

namespace Crystallography;

//必要最小限の情報だけを保存するクラス
[MemoryPackable]
public partial class Crystal2
{
    #region フィールド プライベートメンバーの場合[MemoryPackInclude]が必要
    public enum DataType { None = 0, AMCSD = 1, COD = 2, }


    [MemoryPackInclude]
    private byte[][] cellBytes;

    public int argb;

    public float density;

    public string name;

    public string note;

    public string jour;

    public string auth;

    public string sect;

    public string formula;//計算可能な場合は。

    public short sym;

    public List<Atoms2> atoms;

    public float[] d;//強度8位までのd値

    public string fileName;

    public byte datatype; //DatabaseType (AMCSD, CODなど) の値を格納する

    #endregion

    #region プロパティ
    [MemoryPackIgnore]
    public static byte ID => 0;

    /// <summary>a,b,c,alpha,beta,gammaの順番. 単位はAと度. エラーは 「9.726|5|」のような形式で表現</summary>
    [MemoryPackIgnore]
    public string[] CellTexts
    {
        get => cellBytes == null ? null : Array.ConvertAll(cellBytes, ToString);
        // (260320Ch) null 指定時も内部状態を同期して古い値を残さない
        set => cellBytes = value == null ? null : Array.ConvertAll(value, ToBytes);
    }


    /// <summary>a,b,c,α,β,γ の順番. Getのみ. 長さはA, 角度は度単位.</summary>
    [MemoryPackIgnore]
    public ((double A, double B, double C, double Alpha, double Beta, double Gamma) Values, (double A, double B, double C, double Alpha, double Beta, double Gamma) Errors) Cell
    {
        get
        {
            var c = CellTexts.Select(t => Decompose(t)).ToArray();
            return ((c[0].Value, c[1].Value, c[2].Value, c[3].Value, c[4].Value, c[5].Value),
                     (c[0].Error, c[1].Error, c[2].Error, c[3].Error, c[4].Error, c[5].Error));
        }
    }

    /// <summary>a,b,c,α,β,γ の順番. Getのみ. 長さはA, 角度は度単位. エラーの値は含まない.</summary>
    [MemoryPackIgnore]
    public (double A, double B, double C, double Alpha, double Beta, double Gamma) CellOnlyValue
    {
        get
        {
            //var cellTexts = CellTexts; // ローカルに持つことで呼び出しコスト削減 // 261111Cl 変更前 (3 行)
            //return (
            //    DecomposeOnlyValue(cellTexts[0]), DecomposeOnlyValue(cellTexts[1]), DecomposeOnlyValue(cellTexts[2]),
            //    DecomposeOnlyValue(cellTexts[3]), DecomposeOnlyValue(cellTexts[4]), DecomposeOnlyValue(cellTexts[5]));
            // 261111Cl 変更: 6 本の文字列を作らず、4bit 符号を stackalloc の文字列へ直接展開して読む (値はビット一致)
            var cb = cellBytes;
            return (
                ValueFromBytes(cb[0]), ValueFromBytes(cb[1]), ValueFromBytes(cb[2]),
                ValueFromBytes(cb[3]), ValueFromBytes(cb[4]), ValueFromBytes(cb[5]));
        }
    }

    /// <summary>a,b,c,α,β,γ の順番. Getのみ. 長さはA, 角度は度単位. エラーの値は含まない.</summary>
    [MemoryPackIgnore]
    public (float A, float B, float C, float Alpha, float Beta, float Gamma) CellOnlyValueFloat
    {
        get
        {
            //var cellTexts = CellTexts; // ローカルに持つことで呼び出しコスト削減 // 261111Cl 変更前 (3 行)
            //return (
            //    (float)DecomposeOnlyValue(cellTexts[0]), (float)DecomposeOnlyValue(cellTexts[1]), (float)DecomposeOnlyValue(cellTexts[2]),
            //    (float)DecomposeOnlyValue(cellTexts[3]), (float)DecomposeOnlyValue(cellTexts[4]), (float)DecomposeOnlyValue(cellTexts[5]));
            // 261111Cl 変更: 上の CellOnlyValue と同じく文字列を作らない
            var cb = cellBytes;
            return (
                (float)ValueFromBytes(cb[0]), (float)ValueFromBytes(cb[1]), (float)ValueFromBytes(cb[2]),
                (float)ValueFromBytes(cb[3]), (float)ValueFromBytes(cb[4]), (float)ValueFromBytes(cb[5]));
        }
    }


    /// <summary>a,b,c,α,β,γ の順番. Getのみ. 長さはnm, 角度はradian.</summary>
    [MemoryPackIgnore]
    public ((double A, double B, double C, double Alpha, double Beta, double Gamma) Values, (double A, double B, double C, double Alpha, double Beta, double Gamma) Errors) Cell_nm_radian
    {
        get
        {
            var (Values, Errors) = Cell;
            return ((Values.A / 10, Values.B / 10, Values.C / 10, Values.Alpha / 180 * Math.PI, Values.Beta / 180 * Math.PI, Values.Gamma / 180 * Math.PI),
                    (Errors.A / 10, Errors.B / 10, Errors.C / 10, Errors.Alpha / 180 * Math.PI, Errors.Beta / 180 * Math.PI, Errors.Gamma / 180 * Math.PI));
        }
    }

    /// <summary>a,b,c,α,β,γ の順番. Getのみ. 長さはnm, 角度はradian. エラーの値は含まない.</summary>
    [MemoryPackIgnore]
    public (double A, double B, double C, double Alpha, double Beta, double Gamma) CellOnlyValue_nm_radian
    {
        get
        {
            var (A, B, C, Alpha, Beta, Gamma) = CellOnlyValue;
            return ((A / 10, B / 10, C / 10, Alpha / 180 * Math.PI, Beta / 180 * Math.PI, Gamma / 180 * Math.PI));
        }
    }
    #endregion

    public Crystal2()
    {
    }

    public Crystal ToCrystal() => GetCrystal(this);

    public static Crystal GetCrystal(Crystal2 c)
    {
        if (c == null)
            return null;

        var cell = c.Cell_nm_radian;

        //var atom = new List<Atoms>();
        // (260320Ch) 要素数が分かっているので事前確保して再割り当てを減らす
        var atom = new List<Atoms>(c.atoms?.Count ?? 0);
        if (c.atoms != null)
            foreach (Atoms2 a in c.atoms)
            {
                var pos = a.PositionTexts.Select(x => Decompose(x, c.sym)).ToArray();
                var occ = Decompose(a.OccText);
                var iso = Decompose(a.IsoText);
                //Atoms2の単位はA^2なので、nm^2に変換
                iso = (iso.Value / 100, double.IsNaN(iso.Error) ? iso.Error : iso.Error / 100);

                (double Value, double Error)[] aniso = a.AnisoTexts != null ? a.AnisoTexts.Select(x => Decompose(x)).ToArray() :
                     [(0.0, double.NaN), (0.0, double.NaN), (0.0, double.NaN), (0.0, double.NaN), (0.0, double.NaN), (0.0, double.NaN)];

                var anisoValues = a.IsU ? aniso.Select(an => an.Value / 100).ToArray() : aniso.Select(an => an.Value).ToArray();
                var anisoErrors = a.IsU ? aniso.Select(an => an.Error / 100).ToArray() : aniso.Select(an => an.Error).ToArray();
                var _atom = new Atoms(
                        a.Label, a.AtomNo, a.SubXray, a.SubElectron, null, c.sym,
                        new Vector3DBase(pos[0].Value, pos[1].Value, pos[2].Value),
                        new Vector3DBase(pos[0].Error, pos[1].Error, pos[2].Error),
                        occ.Value, occ.Error,
                        new DiffuseScatteringFactor(a.IsU ? DiffuseScatteringFactor.Type.U : DiffuseScatteringFactor.Type.B, a.IsIso,
                            iso.Value, iso.Error, anisoValues, anisoErrors, cell.Values)
                        );
                atom.Add(_atom);
                atom[^1].ResetVesta();

                //AtomNoが255(重水素D)だった時の処理
                if (atom[^1].AtomicNumber == 255)
                {
                    atom[^1].AtomicNumber = 1;
                    atom[^1].Isotope = [0.0, 100.0, 0.0];
                }
            }

        var bonds = Bonds.GetVestaBonds(atom.Select(a => a.AtomicNumber));

        return new Crystal(cell.Values, cell.Errors, c.sym, c.name, System.Drawing.Color.FromArgb(c.argb), new Matrix3D(), [.. atom], (c.note, c.auth, c.jour, c.sect), bonds);
    }

    public static Crystal2 FromCrystal(Crystal c)
    {
        if (c == null) return null;
        var c2 = new Crystal2
        {
            sym = (short)c.SymmetrySeriesNumber,
            name = c.Name,
            note = c.Note,
            argb = c.Argb,
            auth = c.PublAuthorName,
            sect = c.PublSectionTitle,
            jour = c.Journal,
            formula = c.ChemicalFormulaSum,
            density = (float)c.Density,
            CellTexts = [
                Compose(c.A * 10, c.A_err * 10), Compose(c.B * 10, c.B_err * 10), Compose(c.C * 10, c.C_err * 10),
                Compose(c.Alpha /Math.PI*180, c.Alpha_err/Math.PI*180), Compose(c.Beta /Math.PI*180, c.Beta_err/Math.PI*180), Compose(c.Gamma /Math.PI*180, c.Gamma_err/Math.PI*180) ],
            atoms = []
        };

        foreach (Atoms a in c.Atoms)
        {
            var atom2 = new Atoms2
            {
                Label = a.Label,
                AtomNo = (byte)a.AtomicNumber,
                SubXray = (byte)a.SubNumberXray,
                SubElectron = (byte)a.SubNumberElectron,
                PositionTexts = [Compose(a.X, a.X_err), Compose(a.Y, a.Y_err), Compose(a.Z, a.Z_err)],
                OccText = Compose(a.Occ, a.Occ_err),
                IsIso = a.Dsf.UseIso,
                IsU = a.Dsf.OriginalType == DiffuseScatteringFactor.Type.U,
                IsoText = Compose(a.Dsf.Iso * 100, a.Dsf.Iso_err * 100),
            };

            atom2.AnisoTexts = atom2.IsU ?
                [ Compose(a.Dsf.Aniso11*100, a.Dsf.Aniso11_err*100),
                            Compose(a.Dsf.Aniso22*100, a.Dsf.Aniso22_err*100),
                            Compose(a.Dsf.Aniso33*100, a.Dsf.Aniso33_err*100),
                            Compose(a.Dsf.Aniso12*100, a.Dsf.Aniso12_err*100),
                            Compose(a.Dsf.Aniso23*100, a.Dsf.Aniso23_err*100),
                            Compose(a.Dsf.Aniso31*100, a.Dsf.Aniso31_err*100)] :
                [ Compose(a.Dsf.Aniso11, a.Dsf.Aniso11_err),
                            Compose(a.Dsf.Aniso22, a.Dsf.Aniso22_err),
                            Compose(a.Dsf.Aniso33, a.Dsf.Aniso33_err),
                            Compose(a.Dsf.Aniso12, a.Dsf.Aniso12_err),
                            Compose(a.Dsf.Aniso23, a.Dsf.Aniso23_err),
                            Compose(a.Dsf.Aniso31, a.Dsf.Aniso31_err)];

            c2.atoms.Add(atom2);
        }
        return c2;
    }

    /// <summary>MemoryPackでシリアライズして byte[] 配列を返す</summary>
    /// <returns></returns>
    public byte[] Serialize() => Serialize(this);

    /// <summary>静的メソッド　MemoryPackでシリアライズしてbyte[]配列を返す</summary>
    /// <param name="c"></param>
    /// <returns></returns>
    public static byte[] Serialize(Crystal2 c )
    {
        using var compressor = new BrotliCompressor(System.IO.Compression.CompressionLevel.Optimal, 22);
        MemoryPackSerializer.Serialize(compressor, c);
        return compressor.ToArray();
    }

    /// <summary>MemoryPackでシリアライズされた byte[] 配列からCrystal2を返す. 不適切な入力値だった場合は null 返し。</summary>
    /// <param name="bytes"></param>
    /// <returns></returns>
    public static Crystal2 Deserialize(byte[] bytes)
    {
        try
        {
            using var decompressor = new BrotliDecompressor();// Decompression(require using)
            return MemoryPackSerializer.Deserialize<Crystal2>(decompressor.Decompress(bytes));
        }
        catch { return null; }
    }

    [MemoryPackIgnore]
    private static readonly CultureInfo culture = CultureInfo.InvariantCulture;
    [MemoryPackIgnore]
    private static readonly NumberStyles style = NumberStyles.Number;
    [MemoryPackIgnore]
    private static readonly StringComparison Ord = StringComparison.Ordinal;
    [MemoryPackIgnore]
    // 261111Cl 注: ToString(byte[]) は nibbleChars から直接展開するようになり、この表はもう参照しない (4bit 符号と文字の対応の記録として残す)
    static readonly string[] toStringDic =
        [
            #region 
            "00"            ,"10"            ,"20"            ,"30"            ,"40"            ,"50"            ,"60"            ,"70"            ,"80"            ,"90"            ,".0"            ,"/0"            ,"-0"            ,"|0"            ,"E0"            ,"0"
            ,"01"            ,"11"            ,"21"            ,"31"            ,"41"            ,"51"            ,"61"            ,"71"            ,"81"            ,"91"            ,".1"            ,"/1"            ,"-1"            ,"|1"            ,"E1"            ,"1"
            ,"02"            ,"12"            ,"22"            ,"32"            ,"42"            ,"52"            ,"62"            ,"72"            ,"82"            ,"92"            ,".2"            ,"/2"            ,"-2"            ,"|2"            ,"E2"            ,"2"
            ,"03"            ,"13"            ,"23"            ,"33"            ,"43"            ,"53"            ,"63"            ,"73"            ,"83"            ,"93"            ,".3"            ,"/3"            ,"-3"            ,"|3"            ,"E3"            ,"3"
            ,"04"            ,"14"            ,"24"            ,"34"            ,"44"            ,"54"            ,"64"            ,"74"            ,"84"            ,"94"            ,".4"            ,"/4"            ,"-4"            ,"|4"            ,"E4"            ,"4"
            ,"05"            ,"15"            ,"25"            ,"35"            ,"45"            ,"55"            ,"65"            ,"75"            ,"85"            ,"95"            ,".5"            ,"/5"            ,"-5"            ,"|5"            ,"E5"            ,"5"
            ,"06"            ,"16"            ,"26"            ,"36"            ,"46"            ,"56"            ,"66"            ,"76"            ,"86"            ,"96"            ,".6"            ,"/6"            ,"-6"            ,"|6"            ,"E6"            ,"6"
            ,"07"            ,"17"            ,"27"            ,"37"            ,"47"            ,"57"            ,"67"            ,"77"            ,"87"            ,"97"            ,".7"            ,"/7"            ,"-7"            ,"|7"            ,"E7"            ,"7"
            ,"08"            ,"18"            ,"28"            ,"38"            ,"48"            ,"58"            ,"68"            ,"78"            ,"88"            ,"98"            ,".8"            ,"/8"            ,"-8"            ,"|8"            ,"E8"            ,"8"
            ,"09"            ,"19"            ,"29"            ,"39"            ,"49"            ,"59"            ,"69"            ,"79"            ,"89"            ,"99"            ,".9"            ,"/9"            ,"-9"            ,"|9"            ,"E9"            ,"9"
            ,"0."            ,"1."            ,"2."            ,"3."            ,"4."            ,"5."            ,"6."            ,"7."            ,"8."            ,"9."            ,".."            ,"/."            ,"-."            ,"|."            ,"E."            ,"."            
            ,"0/"            ,"1/"            ,"2/"            ,"3/"            ,"4/"            ,"5/"            ,"6/"            ,"7/"            ,"8/"            ,"9/"            ,"./"            ,"//"            ,"-/"            ,"|/"            ,"E/"            ,"/"            
            ,"0-"            ,"1-"            ,"2-"            ,"3-"            ,"4-"            ,"5-"            ,"6-"            ,"7-"            ,"8-"            ,"9-"            ,".-"            ,"/-"            ,"--"            ,"|-"            ,"E-"            ,"-"            
            ,"0|"            ,"1|"            ,"2|"            ,"3|"            ,"4|"            ,"5|"            ,"6|"            ,"7|"            ,"8|"            ,"9|"            ,".|"            ,"/|"            ,"-|"            ,"||"            ,"E|"            ,"|"            
            ,"0E"            ,"1E"            ,"2E"            ,"3E"            ,"4E"            ,"5E"            ,"6E"            ,"7E"            ,"8E"            ,"9E"            ,".E"            ,"/E"            ,"-E"            ,"|E"            ,"EE"            ,"E"
            ,"0"            ,"1"            ,"2"            ,"3"            ,"4"            ,"5"            ,"6"            ,"7"            ,"8"            ,"9"            ,"."            ,"/"            ,"-"            ,"|"            ,"E"            ,""
            #endregion
        ];

    //[MemoryPackIgnore] // 261111Cl 変更前 (5 行): 16 文字の対応を FrozenDictionary で引いていた。下の ToNibble (switch 式) にした
    //static readonly FrozenDictionary<char, byte> toByteDic = new Dictionary<char, byte>()
    //{
    //    { '0', 0 }, { '1', 1 }, { '2', 2 }, { '3', 3 },  { '4', 4 }, { '5', 5 }, { '6', 6 }, { '7', 7 }, { '8', 8 }, { '9', 9 }, { '.', 10 }, { '/', 11 }, { '-', 12 }, { '(', 13 }, { ')', 13 }, { 'E', 14 },
    //}.ToFrozenDictionary();

    /// <summary>261111Cl 追加: 文字を 4bit 符号に直す (旧 toByteDic と同じ対応)。対応の無い文字は -1</summary>
    private static int ToNibble(char c) => c switch
    {
        >= '0' and <= '9' => c - '0',
        '.' => 10,
        '/' => 11,
        '-' => 12,
        '(' or ')' => 13,
        'E' => 14,
        _ => -1,
    };

    /// <summary></summary>
    /// <param name="s"></param>
    /// <returns></returns>
    public static byte[] ToBytes(string s)
    {
        s = s.Trim().TrimEnd().Replace('e', 'E');
        if (s.Length == 0 || s == "?" || s == "NaN")
            return [255];
        else if (s == "0")
            return [(240 + 0)];
        else
        {
            if (s.StartsWith("0.", Ord))
            {
                if (s == "0.")
                    s = "0";
                else
                    s = s[1..];
            }
            else if (s.StartsWith("-0.", Ord))
                s = s.Replace("-0.", "-.");
            try
            {
                var result = new byte[(s.Length + 1) / 2];
                for (int i = 0; i < s.Length; i += 2)
                {
                    if (i + 1 < s.Length)
                    {
                        //if (toByteDic.TryGetValue(s[i], out var b1) && toByteDic.TryGetValue(s[i + 1], out var b2)) // 261111Cl 変更前
                        if (ToNibble(s[i]) is var b1 and >= 0 && ToNibble(s[i + 1]) is var b2 and >= 0) // 261111Cl 変更
                            result[i / 2] = (byte)(b1 + (b2 << 4));
                    }
                    else
                    {
                        //if (toByteDic.TryGetValue(s[i], out var b1)) // 261111Cl 変更前
                        if (ToNibble(s[i]) is var b1 and >= 0) // 261111Cl 変更
                            result[i / 2] = (byte)(b1 + 240);
                    }
                }
                return result;
            }
            catch (Exception e)
            {
                if (AssemblyState.IsDebug)
                    MessageBox.Show(e.ToString());
                return [255];
            }
        }
    }

    //public static string ToString(byte[] bytes) // 261111Cl 変更前
    //{
    //    var sb = new StringBuilder(bytes.Length * 2);
    //    foreach (var b in bytes)
    //        sb.Append(toStringDic[b]);
    //    return sb.ToString();
    //}
    // 261111Cl 変更: 長さを数えてから string.Create で 1 回だけ確保する。並びは toStringDic[b] と同じ (下位 4bit → 上位 4bit、15 は文字なし)
    public static string ToString(byte[] bytes)
    {
        int len = 0;
        foreach (var b in bytes)
            len += ((b & 15) != 15 ? 1 : 0) + ((b >> 4) != 15 ? 1 : 0);
        return string.Create(len, bytes, static (dst, bs) => DecodeChars(bs, dst));
    }

    /// <summary>261111Cl 追加: 4bit 符号の文字 (toStringDic の 1 文字目・2 文字目と同じ対応)。15 は「文字なし」</summary>
    [MemoryPackIgnore]
    private const string nibbleChars = "0123456789./-|E";

    /// <summary>261111Cl 追加: 4bit 符号のバイト列を dst へ文字として展開し、書いた文字数を返す (dst は bytes.Length * 2 以上)</summary>
    private static int DecodeChars(ReadOnlySpan<byte> bytes, Span<char> dst)
    {
        int n = 0;
        foreach (var b in bytes)
        {
            int lo = b & 15, hi = b >> 4;
            if (lo != 15) dst[n++] = nibbleChars[lo];
            if (hi != 15) dst[n++] = nibbleChars[hi];
        }
        return n;
    }

    /// <summary>261111Cl 追加: 4bit 符号のバイト列のまま DecomposeOnlyValue(ToString(bytes)) と同じ値を返す (文字列を作らない)</summary>
    private static double ValueFromBytes(byte[] bytes)
    {
        Span<char> buf = bytes.Length <= 64 ? stackalloc char[128] : new char[bytes.Length * 2];
        return DecomposeOnlyValue((ReadOnlySpan<char>)buf[..DecodeChars(bytes, buf)]);
    }

    private static (double Value, double Error) Decompose(string str) => Decompose(str, false);
    public static (double Value, double Error) Decompose(string str, int sgnum) => Decompose(str, sgnum >= 430 && sgnum <= 488);


    /// <summary>
    /// 9.726|5|, 1.234|12|E-6 のような文字列を、ValueとErrorに分解してタプルで返す. 
    /// 例外の場合は(double.NaN,double.NaN). 括弧が存在しない場合、Errorは double.NaN. 
    /// </summary>
    /// <param name="str"></param>
    /// <param name="IsHex"></param>
    /// <returns></returns>
    public static (double Value, double Error) Decompose(string str, bool IsHex = false)
    {
        var expValue = 1.0;
        int i;
        if ((i = str.IndexOf("E", Ord)) > 0)
        {
            //_ = double.TryParse("1" + str[i..], style, culture, out expValue); // 261105Cl 変更前 (2 行)
            //str = str[..i];
            // 261105Cl 変更: 指数は NumberStyles.Float で読む。以前の style (= NumberStyles.Number) は指数を読めず、TryParse の失敗で expValue が 0 になり、
            // 値も誤差も 0 になっていた。Compose は誤差つきの値を "2.500(10)E-1" の形で書くので、クリップボードで送った結晶の格子定数・座標・
            // 温度因子がほぼ全部 0 になっていた (tools/AdpCheck r3)。CIF 1.1 の正式な形 "7.896E-1(13)" (内部形式 7.896E-1|13|) のように
            // 指数の後ろに誤差があるときは、誤差を仮数の後ろへ移す。読めない指数は NaN を返す (0 にしない)
            var exp = str[(i + 1)..];
            var bar = exp.IndexOf('|');
            str = bar >= 0 ? str[..i] + exp[bar..] : str[..i];
            if (!double.TryParse("1E" + (bar >= 0 ? exp[..bar] : exp), NumberStyles.Float, culture, out expValue))
                return (double.NaN, double.NaN);
        }
        string valStr;
        double err;
        if ((i = str.IndexOf("|", Ord)) > 0)
        {
            valStr = str[0..i];
            if (str.Length - 1 > i && double.TryParse(str[(i + 1)..^1], style, culture, out err))
            {
                var j = valStr.IndexOf(".", Ord);
                if (j >= 0 && valStr.Length - j - 1 > 0)
                    err *= Math.Pow(10, -valStr.Length + j + 1);
            }
            else
                err = double.NaN;
        }
        else
        {
            valStr = str;
            err = double.NaN;
        }

        if (valStr == "0")
            return (0, err * expValue);

        if (IsHex)
        {
            if (str.Contains(".0833"))
                return (1.0 / 12.0, err);
            else if (str.Contains(".167") || str.Contains(".1667") || str.Contains(".16667") || str.Contains(".166667"))
                return (1.0 / 6.0, err);
            else if (str.Contains(".333"))
                return (1.0 / 3.0, err);
            else if (str.Contains(".667") || str.Contains(".6667") || str.Contains(".66667") || str.Contains(".666667"))
                return (2.0 / 3.0, err);
            else if (str.Contains(".4167") || str.Contains(".41667") || str.Contains(".416667"))
                return (5.0 / 12.0, err);
            else if (str.Contains(".5833"))
                return (7.0 / 12.0, err);
            else if (str.Contains(".8333")|| str.Contains(".83333"))
                return (5.0 / 6.0, err);
            else if (str.Contains(".9167") || str.Contains(".91667") || str.Contains(".916667"))
                return (11.0 / 12.0, err);
        }

        if ((i = valStr.IndexOf("/", Ord)) >= 0)
        {
            if (double.TryParse(valStr.AsSpan()[0..i], style, culture, out var temp0)
                && double.TryParse(valStr.AsSpan()[(i + 1)..^0], style, culture, out var temp1))
                return (temp0 / temp1, double.NaN);
            else
                return (double.NaN, double.NaN);
        }

        return double.TryParse(valStr, style, culture, out var val) ? (val * expValue, err * expValue) : (double.NaN, double.NaN);
    }

    /// <summary>
    /// 9.726|5|, 1.234|12|E-6 のような文字列を、Valueの部分だけ返す. 
    /// 格子定数を変換するときだけに呼ばれる.
    /// </summary>
    /// <param name="str"></param>
    /// <returns></returns>
    public static double DecomposeOnlyValue(string str) => DecomposeOnlyValue(str.AsSpan()); // 261111Cl 変更: 本体は下の span 版へ移した (旧本体はその中にコメントで残す)

    /// <summary>261111Cl 追加: DecomposeOnlyValue(string) の span 版。指数の "1E…" も stackalloc で組み立てて読む (値はビット一致)</summary>
    public static double DecomposeOnlyValue(ReadOnlySpan<char> s)
    {
        //var s = str.AsSpan(); // 261111Cl 変更前 (引数で受ける)
        var expValue = 1.0;
        if (s.Contains('E'))
        {
            int i = s.IndexOf('E');
            //_ = double.TryParse("1" + s[i..].ToString(), style, culture, out expValue); // 261105Cl 変更前
            // 261105Cl 変更: Decompose と同じく指数は NumberStyles.Float で読む (以前は 10 Å 以上の格子定数などが 0 になっていた)。指数の後ろの |誤差| は捨てる
            var exp = s[(i + 1)..];
            if (exp.IndexOf('|') is var bar and >= 0)
                exp = exp[..bar];
            //if (!double.TryParse(string.Concat("1E", exp), NumberStyles.Float, culture, out expValue)) // 261111Cl 変更前
            //    return double.NaN;
            // 261111Cl 変更: "1E" + exp を stackalloc に組み立てる (TryParse(string) は中で同じ span 版を呼ぶので結果は同じ)
            Span<char> e = exp.Length <= 62 ? stackalloc char[64] : new char[exp.Length + 2];
            e[0] = '1';
            e[1] = 'E';
            exp.CopyTo(e[2..]);
            if (!double.TryParse(e[..(exp.Length + 2)], NumberStyles.Float, culture, out expValue))
                return double.NaN;
            s = s[..i];
        }

        if (s.Contains('|'))
            s = s[0..s.IndexOf('|')];

        return double.TryParse(s, style, culture, out var val) ? val * expValue : double.NaN;
    }

    public static string Compose(double val, double err = double.NaN)
    {
        if (double.IsNaN(val))
            return "";
        else if (double.IsNaN(err) || err <= 0)
            //return val.ToString(); // 261105Cl 変更前
            return val.ToString(culture); // 261105Cl 変更: 小数点がカンマの環境でも '.' で書く (ToBytes はカンマを表せず、値が壊れていた)。以下の ToString・Parse も同じ
        else
        {
            //まず、誤差を 23E-6 みたいな形にする
            //var temp = err.ToString("E1"); // 261105Cl 変更前 (この下の 4 行は culture の指定を足しただけ)
            //var errLog = Convert.ToInt32(temp[(temp.IndexOf("E", Ord) + 1)..]) - 1;
            //var valStr = val.ToString("E15");
            //var valLog = Convert.ToInt32(valStr[(valStr.IndexOf("E", Ord) + 1)..]);
            var temp = err.ToString("E1", culture);
            var errStr = temp[..temp.IndexOf("E", Ord)].Replace(".", "");
            var errLog = int.Parse(temp[(temp.IndexOf("E", Ord) + 1)..], NumberStyles.AllowLeadingSign, culture) - 1;

            //vを取りあえず十分な精度で出力する
            var valStr = val.ToString("E15", culture);
            var valLog = int.Parse(valStr[(valStr.IndexOf("E", Ord) + 1)..], NumberStyles.AllowLeadingSign, culture);

            //var result = valLog >= errLog ? valStr[..(valLog - errLog + 2)] : val.ToString("E0")[..1]; // 261105Cl 変更前
            // 261105Cl 変更: 負の値は先頭の '-' の分だけ 1 桁短く切られ、誤差の桁が 1 つずれて 10 倍の誤差として読まれていた (tools/AdpCheck r3)
            var sign = val < 0 ? 1 : 0;
            var result = valLog >= errLog ? valStr[..(valLog - errLog + 2 + sign)] : val.ToString("E0", culture)[..(1 + sign)];

            if (valLog < errLog)
                for (int i = 0; i < errLog - valLog; i++)
                    errStr += "0";

            result += "(" + errStr + ")";

            if (valLog != 0)
                //result += $"E{valLog}"; // 261105Cl 変更前
                result += "E" + valLog.ToString(culture); // 261105Cl 変更: 負号が U+2212 の環境でも '-' で書く

            return result;
        }
    }
}
