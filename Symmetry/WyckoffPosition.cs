using System;
using System.Collections.Generic;
using System.Linq;
using System.Xml.Serialization;

namespace Crystallography;

/// <summary>ワイコフポジション構造体. SymmetryStaticで初期化される.</summary>
[Serializable()]
public readonly struct WyckoffPosition
{
    #region 定数
    const double th = 0.00015;
    const double one_th = 1 - th * 2;
    #endregion

    #region フィールド、プロパティ
    /// <summary>空間群の番号</summary>
    public ushort SymmetrySeriesNumber { get; }

    /// <summary>格子のタイプ</summary>
    public string LatticeType { get; }

    /// <summary>多重度 (整数)</summary>
    public byte Multiplicity { get; }

    /// <summary>ワイコフ文字</summary>
    public string WyckoffLetter { get; }

    /// <summary>ワイコフナンバー (一般位置が0, 特殊になるほど数字が大)</summary>
    public byte WyckoffNumber { get; }

    /// <summary>サイトシンメトリ</summary>
    public string SiteSymmetry { get; }

    /// <summary>等価位置を生成するアフィン対称操作の配列 (260604Cl: Func[] → AffineGen[]。デリゲート除去で高速・データ駆動)</summary>
    [XmlIgnore]
    public AffineGen[] PositionGenerator { get; }

    /// <summary>等価位置の文字列(x,y,zなど)の配列</summary>
    public string[] PositionStr { get; }

    /// <summary>等価位置の対称操作をSymmetryOperationクラスとして格納したもの</summary>
    public SymmetryOperation[] PositionOperations { get; }

    /// <summary>自由度 (このワイコフ位置がx,y,zなどの変数を含む場合はtrue, 含まない場合はfalse)</summary>
    public (bool X, bool Y, bool Z) Free { get; }

    private static readonly char[] separator = [','];

    #endregion

    #region コンストラクタ
    public WyckoffPosition(int symSeries, string latticeType, string wykLetter, int wykNum, string siteSym,
        string[] coordStr,
        AffineGen[] generators,
        SymmetryOperation[] operations = null)
    {
        SymmetrySeriesNumber = (ushort)symSeries;
        LatticeType = latticeType;

        WyckoffLetter = wykLetter;
        WyckoffNumber = (byte)wykNum;
        SiteSymmetry = siteSym;

        PositionStr = coordStr;
        PositionGenerator = generators;
        Multiplicity = (byte)generators.Length;

        PositionOperations = operations;

        // (260320Ch) パターンマッチで null/空配列判定を明示する
        if (PositionStr is not { Length: > 0 })
            Free = (true, true, true);
        else
        {
            var tempStr = PositionStr[0].Split(separator, StringSplitOptions.RemoveEmptyEntries);
            if (tempStr.Length == 3)
                Free = (tempStr[0].Contains('x'), tempStr[1].Contains('y'), tempStr[2].Contains('z'));
            else
                Free = (true, true, true);
        }
    }
    #endregion

    #region メソッド
    /// <summary>与えられたx,y,zで、このワイコフ位置を再生</summary>
    /// <param name="x"></param>
    /// <param name="y"></param>
    /// <param name="z"></param>
    /// <returns></returns>
    public readonly Vector3D[] GeneratePositions(in double x, in double y, in double z)
    {
        var pos = new List<Vector3D>(PositionGenerator.Length);

        for (int i = 0; i < PositionGenerator.Length; i++)
        {
            var (X, Y, Z) = PositionGenerator[i].Apply(x, y, z);

            //当たり判定
            // (260320Ch) LINQ All をループへ置き換えてホットパスの割り当てを減らす
            if (ContainsPosition(pos, X, Y, Z))
                continue;

            var v = new Vector3D(X, Y, Z, false);
            //0~1の範囲に収まるかどうかチェックし、修正
            v.InnerLatticeThis();
            if (PositionOperations != null)
                v.Operation = new SymmetryOperation(PositionOperations[i], SymmetrySeriesNumber);//PositionOperatorsを格納
            pos.Add(v);
        }
        return [.. pos];
    }
    /// <summary>与えられたposがこのWykoffPositionかどうかを判定する</summary>
    /// <param name="pos"></param>
    /// <returns></returns>
    public readonly bool CheckPosition(double x, double y, double z)
    {
        // (260320Ch) Any のラムダ割り当てを避ける
        foreach (var g in PositionGenerator)
        {
            var (X, Y, Z) = g.Apply(x, y, z);
            if (chk(X, x) && chk(Y, y) && chk(Z, z))
                return true;
        }
        return false;
    }

    // (260320Ch) 重複位置の判定をループ化してクロージャ生成を避ける
    private static bool ContainsPosition(List<Vector3D> positions, double x, double y, double z)
    {
        foreach (var p in positions)
            if (chk(z, p.Z) && chk(x, p.X) && chk(y, p.Y))
                return true;

        return false;
    }

    static bool chk(in double d1, in double d2)
    {
        var d = Math.Abs(d1 - d2);
        if (d >= 1 || d < 0) d -= Math.Floor(d);
        if (d > one_th) d -= 1;
        return d <= th;
    }
    #endregion

    #region static メソッド
    /// <summary>引数の空間群による対称操作で映る原子位置(pos)の等価な原子位置をクラスAtomsでかえす</summary>
    /// <param name="Pos"></param>
    /// <param name="SymmetrySeriesNumber"></param>
    /// <returns></returns>
    public static Atoms GetEquivalentAtomsPosition(in (double X, double Y, double Z) Pos, in int SymmetrySeriesNumber)
    {
        if (double.IsNaN(Pos.X) || double.IsNaN(Pos.Y) || double.IsNaN(Pos.Z))
            return new Atoms();

        var atoms = new Atoms();
        var wyck = SymmetryStatic.WyckoffPositions[SymmetrySeriesNumber];
        //まず、もっとも対称性の低いワイコフ位置で原子位置を再生
        atoms.Atom = wyck[0].GeneratePositions(Pos.X, Pos.Y, Pos.Z);

        //ワイコフ位置判定
        atoms.WyckoffLeter = "{";
        atoms.SiteSymmetry = "";
        atoms.Multiplicity = 0;
        string wyckLet = "";
        string siteSym = "";
        int multi = 0;
        int wyckNum = 0;

        var atomsTemp = atoms.Atom.Select(a => (a.X, a.Y, a.Z)).ToArray();

        for (int j = wyck.Length - 1; j >= 0; j--)
        {

            //2020/05/15 一部のtrigonal, hexagonalで正しい判定が出来ない問題の対応
            // (0.2, 0.1, 0)という座標は (x, -x, 0), (x, 2x, 0), (-2x, -x, 0)というワイコフ位置 (P321, 3j)
            // にもかかわらず、正しく判定することが出来ない.
            //2022/06/16 服部さんから、P4_2/nmc(2)でも、(3/4,y,z)が上手く判定できていないとの指摘を受ける

            // そのため、もっとも低対称性で再生した位置全てに対して判定をおこない、一回でもOKだったらそのワイコフ位置だと判定する
            if (wyck[j].Multiplicity == atoms.Atom.Length)
                foreach (var (X, Y, Z) in atomsTemp)
                    if (wyck[j].CheckPosition(X, Y, Z))
                    {
                        multi = wyck[j].Multiplicity;
                        wyckLet = wyck[j].WyckoffLetter;
                        siteSym = wyck[j].SiteSymmetry;
                        wyckNum = j;
                        break;
                    }

            if (multi != 0)
                break;
        }

        atoms.WyckoffLeter = wyckLet;
        atoms.SiteSymmetry = siteSym;
        atoms.Multiplicity = atoms.Atom.Length;
        atoms.WyckoffNumber = wyckNum;

        return atoms;
    }

    #endregion
}

/// <summary>260604Cl アフィン対称操作 (X,Y,Z) = M·(x,y,z) + t + post。
/// 旧 SymmetryStatic.PositionGeneratorList の 1,969 個のラムダと _PositionStringList の文字列を統合し、
/// 可読 CSV(x,y,z 標準表記)から構築する。post は格子併進(A/B/C/I/R1/R2)を base 適用後に最後に加算する後置併進
/// (旧コードの結合順と一致させ bit-exact を保つ)。struct のため JIT がインライン化でき、旧 Func[] より高速。</summary>
public readonly struct AffineGen
{
    // 出力行 r = cxr*x + cyr*y + czr*z + tr + (後置併進)
    readonly double cx0, cy0, cz0, t0;
    readonly double cx1, cy1, cz1, t1;
    readonly double cx2, cy2, cz2, t2;
    readonly double px, py, pz; // 後置併進(格子センタリング)。base 生成器は 0。

    private AffineGen(double cx0, double cy0, double cz0, double t0,
                      double cx1, double cy1, double cz1, double t1,
                      double cx2, double cy2, double cz2, double t2,
                      double px, double py, double pz)
    {
        this.cx0 = cx0; this.cy0 = cy0; this.cz0 = cz0; this.t0 = t0;
        this.cx1 = cx1; this.cy1 = cy1; this.cz1 = cz1; this.t1 = t1;
        this.cx2 = cx2; this.cy2 = cy2; this.cz2 = cz2; this.t2 = t2;
        this.px = px; this.py = py; this.pz = pz;
    }

    /// <summary>(X,Y,Z) = M·(x,y,z) + t + post。post は最後に加算(旧 base 適用→併進 の結合順に一致)。</summary>
    public (double X, double Y, double Z) Apply(double x, double y, double z)
        => (cx0 * x + cy0 * y + cz0 * z + t0 + px,
            cx1 * x + cy1 * y + cz1 * z + t1 + py,
            cx2 * x + cy2 * y + cz2 * z + t2 + pz);

    /// <summary>格子併進を後置で加えた新生成器(A/B/C/I/R1/R2 用)。base 適用後に dx,dy,dz を最後に加算する。</summary>
    public AffineGen AddTranslation(double dx, double dy, double dz)
        => new(cx0, cy0, cz0, t0, cx1, cy1, cz1, t1, cx2, cy2, cz2, t2, px + dx, py + dy, pz + dz);

    /// <summary>"a,b,c" の 3 成分(x,y,z 標準表記)からアフィン変換を構築。整数係数(2x 等)と分数併進(1/3 等)に対応。
    /// 261105Cl 追記: 読めない成分 (小数点、大文字の X、* や括弧、空の成分、分子・分母の無い分数、分母 0、桁あふれ) は FormatException を投げる
    /// (以前は小数点や大文字などで先へ進まず無限ループになっていた)。</summary>
    /// <exception cref="FormatException">261105Cl 追加: 成分を読めないとき</exception>
    /// <param name="allowDecimal">261111Cl 追加: true なら小数 ("x+0.5"、"z-.25"、"1/2.0") も読む (CIF の対称操作の確認用)。
    /// 既定の false では従来どおり小数点を FormatException にする (空間群の表は整数と分数だけ)。</param>
    //public static AffineGen Parse(string a, string b, string c) // 261111Cl 変更前
    public static AffineGen Parse(string a, string b, string c, bool allowDecimal = false) // 261111Cl allowDecimal を追加
    {
        //var (cx0, cy0, cz0, t0) = ParseComponent(a); // 261111Cl 変更前 (3 行)
        //var (cx1, cy1, cz1, t1) = ParseComponent(b);
        //var (cx2, cy2, cz2, t2) = ParseComponent(c);
        var (cx0, cy0, cz0, t0) = ParseComponent(a, allowDecimal); // 261111Cl
        var (cx1, cy1, cz1, t1) = ParseComponent(b, allowDecimal); // 261111Cl
        var (cx2, cy2, cz2, t2) = ParseComponent(c, allowDecimal); // 261111Cl
        return new(cx0, cy0, cz0, t0, cx1, cy1, cz1, t1, cx2, cy2, cz2, t2, 0, 0, 0);
    }

    //private static (double cx, double cy, double cz, double t) ParseComponent(string s) // 261111Cl 変更前
    private static (double cx, double cy, double cz, double t) ParseComponent(string s, bool allowDecimal = false) // 261111Cl allowDecimal を追加
    {
        // 261105Cl 追加: 読めない成分は FormatException。整数は桁あふれも FormatException にそろえる
        FormatException bad(string why) => new($"AffineGen.Parse: cannot read \"{s}\" ({why})"); // 261105Cl 追加
        int toInt(string digits) => int.TryParse(digits, System.Globalization.NumberStyles.None, System.Globalization.CultureInfo.InvariantCulture, out var v) ? v : throw bad($"integer '{digits}'"); // 261105Cl 追加
        // 261111Cl 追加: 数の値。小数点を含むもの (allowDecimal のときだけ現れる) は double で読み、整数は従来どおり int で読む (整数だけの式の結果は従来と同じ)
        double toNum(string digits) => !digits.Contains('.') ? toInt(digits)
            : double.TryParse(digits, System.Globalization.NumberStyles.AllowDecimalPoint, System.Globalization.CultureInfo.InvariantCulture, out var v) ? v : throw bad($"number '{digits}'");
        // 261111Cl 追加: 数字の並びを読み進める。allowDecimal なら小数点 1 つと続く数字も含める
        void skipNumber(ref int k)
        {
            while (k < s.Length && char.IsDigit(s[k])) k++;
            if (allowDecimal && k < s.Length && s[k] == '.')
            {
                k++;
                while (k < s.Length && char.IsDigit(s[k])) k++;
            }
        }
        if (string.IsNullOrWhiteSpace(s)) throw bad("empty"); // 261105Cl 追加
        double cx = 0, cy = 0, cz = 0, t = 0; int i = 0;
        while (i < s.Length)
        {
            int sign = 1;
            if (s[i] == '+') { i++; }
            else if (s[i] == '-') { sign = -1; i++; }
            else if (char.IsWhiteSpace(s[i])) { i++; continue; }
            while (i < s.Length && char.IsWhiteSpace(s[i])) i++; // 261105Cl 追加: 符号の後の空白 ('z - 1/4')。以前は符号が捨てられ +1/4 と読んでいた
            int ns = i;
            //while (i < s.Length && char.IsDigit(s[i])) i++; // 261111Cl 変更前
            skipNumber(ref i); // 261111Cl
            string num = s[ns..i];
            if (i < s.Length && (s[i] == 'x' || s[i] == 'y' || s[i] == 'z'))
            {
                //double coef = sign * (num.Length > 0 ? int.Parse(num) : 1); // 261105Cl 変更前
                //double coef = sign * (num.Length > 0 ? toInt(num) : 1); // 261105Cl // 261111Cl 変更前
                double coef = sign * (num.Length > 0 ? toNum(num) : 1); // 261111Cl
                if (s[i] == 'x') cx += coef; else if (s[i] == 'y') cy += coef; else cz += coef;
                i++;
            }
            else if (i < s.Length && s[i] == '/')
            {
                //i++; int ds = i; while (i < s.Length && char.IsDigit(s[i])) i++; // 261111Cl 変更前
                i++; int ds = i; skipNumber(ref i); // 261111Cl
                //t += sign * ((double)int.Parse(num) / (double)int.Parse(s[ds..i])); // 261105Cl 変更前
                if (num.Length == 0 || ds == i) throw bad("fraction without numerator or denominator"); // 261105Cl 追加
                //var den = toInt(s[ds..i]); // 261105Cl 追加 // 261111Cl 変更前
                var den = toNum(s[ds..i]); // 261111Cl
                if (den == 0) throw bad("zero denominator"); // 261105Cl 追加
                //t += sign * ((double)toInt(num) / den); // 261105Cl // 261111Cl 変更前
                t += sign * (toNum(num) / den); // 261111Cl (整数どうしなら (double)分子 / 分母 で従来と同じ)
            }
            else
            {
                // 261105Cl 追加: ここで止まる文字 (符号・空白以外) は次の周回でも読めず、以前は i が進まないまま無限ループになっていた。
                //  数字が無い (符号の後に項が無い 'x+'、'x--y' など) のも読めない
                if (num.Length == 0) throw bad(i < s.Length ? $"unexpected '{s[i]}' at {i}" : "sign without a term"); // 261105Cl 追加
                if (i < s.Length && s[i] != '+' && s[i] != '-' && !char.IsWhiteSpace(s[i])) throw bad($"unexpected '{s[i]}' at {i}"); // 261105Cl 追加
                //t += sign * (num.Length > 0 ? int.Parse(num) : 0); // 261105Cl 変更前
                //t += sign * (num.Length > 0 ? toInt(num) : 0); // 261105Cl // 261111Cl 変更前
                t += sign * (num.Length > 0 ? toNum(num) : 0); // 261111Cl
            }
        }
        return (cx, cy, cz, t);
    }
}
