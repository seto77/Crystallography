#region using
using System;
using System.Collections.Frozen;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Diagnostics;
using System.Drawing;
using System.IO;
using System.Linq;
// using System.Linq.Dynamic.Core;       // 260427Cl DynamicExpresso へ置換
// using System.Linq.Expressions;         // 260427Cl Expression.Parameter 不要になったため削除
using DynamicExpresso;
using System.Text;
using System.Text.RegularExpressions; // 260427Cl 追加: 対称操作式の暗黙乗算変換に使用
using System.Threading;
using ZLinq;
using V3 = OpenTK.Mathematics.Vector3d;
#endregion

namespace Crystallography;

//public class ConvertCrystalData // 261111Cl 変更前
public partial class ConvertCrystalData // 261111Cl 変更: [GeneratedRegex] のため partial (このクラスの宣言はこのファイルだけ)
{
    private static readonly StringComparison Ord = StringComparison.Ordinal;
    /// <summary>261111Cl 追加: ファイル形式の数値の読み書きに使うカルチャ (小数点 '.'。日本語・英語の環境と同じ書式)</summary>
    private static readonly System.Globalization.CultureInfo Inv = System.Globalization.CultureInfo.InvariantCulture;

    // 260427Cl 追加: CIF/特殊位置の "2x" のような暗黙乗算を "2*x" に変換 (DynamicExpresso は暗黙乗算非対応)
    //private static readonly Regex ImplicitMulRegex = new(@"([0-9])([xyz])", RegexOptions.Compiled); // 261111Cl 変更前
    // 261111Cl 変更: 実行時にコンパイルする Regex の代わりにソース生成 (パターン・オプションは同じ。呼び出し側はそのまま ImplicitMulRegex.Replace)
    private static Regex ImplicitMulRegex => ImplicitMulRegexGenerated();
    [GeneratedRegex(@"([0-9])([xyz])")]
    private static partial Regex ImplicitMulRegexGenerated();

    /// <summary>261111Cl 追加: scaleValueText が読む CIF の数値 ("0.071(1)"、"1.2E-3(3)"、"1.2(3)E-3" など)。旧は呼ぶたびに静的 Regex.Match へパターン文字列を渡していた (パターンは同じ)</summary>
    [GeneratedRegex(@"^([+-]?(?:\d+\.?\d*|\.\d+))(?:\((\d+)\))?(?:[eE]([+-]?\d+))?(?:\((\d+)\))?$")]
    private static partial Regex CifNumberWithUncertaintyRegex();

    #region 文字列操作の高速化を狙った関数群

    static bool Is90(string s) => s == "90" || s == "90.0" || s == "90.00" || s == "90.000" || s == "90.0000";

    static void Replace(ref string str, string oldValue, string newValue)
    {
        if (str.Contains(oldValue))
            str = str.Replace(oldValue, newValue);
    }

    static void TrimStartEnd(ref string str)
    {
        if (str.StartsWith(' '))
            str = str.TrimStart(' ');
        if (str.EndsWith(' '))
            str = str.TrimEnd(' ');
    }
    public static void EscapeString(ref string s)
    {
        Replace(ref s, "&amp;", "&");
        Replace(ref s, "&lt;", "<");
        Replace(ref s, "&gt;", ">");
        Replace(ref s, "&quot;", "\"");
        Replace(ref s, "&nbsp;", " ");
        Replace(ref s, "&#39;", "'");
    }

    // 261111Cl 変更前: 呼ぶたびに配列を作り直す (O(N)) ので、CIF の前処理で行を消すたびに全体を複写していた (O(N²))。
    // 前処理は List<string> の RemoveRange / Insert と 1 パスのフィルタにしたので使わない
    //static void RemoveAt(ref string[] str, int index) => str = [..str[..index] , ..str[(index + 1)..]];
    //static void RemoveRange(ref string[] str, int start, int count) => str = [..str[..start] ,.. str[(start + count)..]];
    //
    //static void Insert(ref string[] str, int index, string value) => str = [.. str[..index], value, .. str[index..]];

    #endregion

    #region CrystalList(xml形式)の読み込み/書き込み
    public static bool SaveCrystalListXml(Crystal[] crystals, string filename)
    {
        if (crystals != null && crystals.Length != 0 && crystals[0].FlexibleMode)
            crystals[0].FlexiblePlane = crystals[0].Plane;

        try
        {
            var serializer = new System.Xml.Serialization.XmlSerializer(typeof(Crystal[]));
            //var fs = new FileStream(filename, FileMode.Create);
            //serializer.Serialize(fs, crystals);
            //fs.Close();
            using var fs = new FileStream(filename, FileMode.Create); // 260718Cl: Serialize 例外時に FileStream がリークしファイルロックが残るため using 化
            serializer.Serialize(fs, crystals);
            return true;
        }
        catch { return false; }
    }


    //CrystalListを読み込むとき
    public static Crystal[] ConvertToCrystalList(string filename)
    {
        var cry = Array.Empty<Crystal>();
        //if (filename.ToLower().EndsWith("xml", Ord))//XML形式のリストを読み込んだとき // 261111Cl 変更前
        if (filename.EndsWith("xml", StringComparison.OrdinalIgnoreCase))//XML形式のリストを読み込んだとき // 261111Cl 変更: 一時文字列と現在のカルチャへの依存をなくす
        {
            if (new FileInfo(filename).Length > 10000000)//なぜかファイルが3GBとかになったことが有ったので、それに対する対処. 10MB以上だったらスキップすることにした.
                return cry;

            #region old code
            //プロパティ文字列が変更にたいする対処
            // 261105Cl 修正: 開くだけで元のファイルを書き換えていた (全行を置換して同じファイルへ上書き。読み取り専用なら読込そのものが失敗し、
            //  書き込みの途中で失敗すると原本が壊れる)。旧書式との互換の置換はメモリ上だけで行い、ファイルへは書き戻さない (PdiFile の 5c4c99a と同じ)
            string xmlText; // 261105Cl 追加: 互換の置換を済ませた XML 本文
            try
            {
                //var reader = new StreamReader(filename, Encoding.GetEncoding("UTF-8")); // 261111Cl 変更前
                var reader = new StreamReader(filename, Encoding.UTF8); // 261111Cl 変更: 同じ UTF-8 (読み込みは BOM の有無を問わない)
                var strList = new List<string>();
                string tempstr;
                while ((tempstr = reader.ReadLine()) != null)
                {
                    tempstr = tempstr.Replace("Kprime0", "Kp0");
                    tempstr = tempstr.Replace("Birch_Murnaghan", "BM3");
                    strList.Add(tempstr);
                }

                reader.Close();

                //var writer = new StreamWriter(filename, false, Encoding.GetEncoding("UTF-8")); // 261105Cl 変更前: 元のファイルへ書き戻していた
                //for (int i = 0; i < strList.Count; i++)
                //    writer.WriteLine(strList[i]);
                //writer.Flush();
                //writer.Close();
                xmlText = string.Join(Environment.NewLine, strList); // 261105Cl 追加
            }
            catch { return null; };
            //プロパティ文字列が変更にたいする対処　ここまで
            #endregion old code
            try
            {
                //using var fs = new FileStream(filename, FileMode.Open); // 261105Cl 変更前
                using var fs = new StringReader(xmlText); // 261105Cl 変更: 置換済みの本文から読む
                var serializer = new System.Xml.Serialization.XmlSerializer(typeof(Crystal[]));


                #region 261111Cl 変更前: BackgroundWorker + UI スレッドで Thread.Sleep(1) を最大 1 万回ポーリング。例外や 0 件のときは約 2.5 分固まってから空を返していた
                //var worker = new BackgroundWorker();
                //worker.DoWork += (sender, e) => cry = (Crystal[])serializer.Deserialize(fs);
                //worker.RunWorkerAsync();
                //
                //var timeout = 10000;
                //var sw = new Stopwatch();
                //for (int i = 0; i < timeout; i++)
                //{
                //    if (cry.Length != 0)
                //        break;
                //    if (i == timeout - 1)
                //    {
                //        throw new Exception();
                //    }
                //    Thread.Sleep(1);
                //}
                #endregion
                // 261111Cl 変更: Task で読み、終わるまで (上限 30 秒) 待つ。読込の例外はその場で Wait から上がり、下の catch で従来どおり空の配列を返す。
                // 0 件のリストもすぐ返る (従来は 0 件を「未完了」と区別できず、上限まで待っていた)。30 秒を過ぎたら TimeoutException → 空の配列
                var task = System.Threading.Tasks.Task.Run(() => (Crystal[])serializer.Deserialize(fs));
                if (!task.Wait(TimeSpan.FromSeconds(30)))
                    throw new TimeoutException();
                cry = task.Result;

                foreach (var c in cry.Where(e => e.FlexiblePlane != null && e.FlexiblePlane.Count > 0))
                    c.Plane = c.FlexiblePlane;
            }
            catch { }
        }
        //else if (filename.EndsWith("out", Ord))//SMAP形式を読み込んだとき // 261105Cl 変更前
        else if (filename.EndsWith("out", StringComparison.OrdinalIgnoreCase))//SMAP形式を読み込んだとき // 261105Cl 変更 (監査 A21): 大文字の ".OUT" も読む (xml の判定と同じく大文字小文字を区別しない)
        {
            var stringList = new List<string>();
            string strTemp;
            var reader = new StreamReader(filename);
            while ((strTemp = reader.ReadLine()) != null)
                stringList.Add(strTemp);
            reader.Close();
            cry = ConvertFromSMAP([.. stringList]);
        }

        for (int i = 0; i < cry.Length; i++)
        {
            cry[i].Reset();
            cry[i].SaveInitialCellConstants();
        }

        return cry;
    }

    #endregion

    #region SMAPの出力ファイル(*.out)読込

    /// <summary>SMAPのoutファイル読み込み</summary>
    /// <param name="str"></param>
    /// <returns></returns>
    public static Crystal[] ConvertFromSMAP(string[] str)
    {
        string wavelength = "";
        string ChemicalFormula = "";
        double a = 0, b = 0, c = 0, alpha = 0, beta = 0, gamma = 0;
        int spaceGroupSeriesNum = 0;
        //空間群
        int line = 0;

        for (; line < str.Length; line++)
            if (str[line].StartsWith(" wave length", Ord))
            {
                wavelength = str[line];
                break;
            }
        for (; line < str.Length; line++)
            if (str[line].StartsWith(" Chemical Formula", Ord))
            {
                ChemicalFormula = str[line];
                break;
            }
        for (; line < str.Length; line++)
            if (str[line].StartsWith(" Cell Constants", Ord))
            {
                string[] cellconstants = str[line + 1].Split(' ', true);
                //a = Convert.ToDouble(cellconstants[0]) / 10.0; // 261111Cl 変更前 (6 行)
                //b = Convert.ToDouble(cellconstants[1]) / 10.0;
                //c = Convert.ToDouble(cellconstants[2]) / 10.0;
                //alpha = Convert.ToDouble(cellconstants[3]) / 180 * Math.PI;
                //beta = Convert.ToDouble(cellconstants[4]) / 180 * Math.PI;
                //gamma = Convert.ToDouble(cellconstants[5]) / 180 * Math.PI;
                // 261111Cl 変更: ファイルの数値は現在のカルチャでなく Invariant で読む (小数点がカンマの環境で "5.43" が 543 になっていた)
                a = Convert.ToDouble(cellconstants[0], Inv) / 10.0;
                b = Convert.ToDouble(cellconstants[1], Inv) / 10.0;
                c = Convert.ToDouble(cellconstants[2], Inv) / 10.0;
                alpha = Convert.ToDouble(cellconstants[3], Inv) / 180 * Math.PI;
                beta = Convert.ToDouble(cellconstants[4], Inv) / 180 * Math.PI;
                gamma = Convert.ToDouble(cellconstants[5], Inv) / 180 * Math.PI;
                break;
            }

        for (; line < str.Length; line++)
            if (str[line].StartsWith(" Space Group Number", Ord))
            {
                string s = (str[line].Split('=', true)[1]).Trim();
                int num = Convert.ToInt32(s.Substring(1, 3));
                int sub = Convert.ToInt32(s.Substring(4, 1));
                spaceGroupSeriesNum = SymmetryStatic.GetSeriesNumber(num, sub);
                break;
            }

        int n = 0;
        //ここから原始座標読み取り
        var cry = new List<Crystal>();
        while (line < str.Length)
        {
            for (; line < str.Length; line++)
                if (str[line].StartsWith("No =", Ord))
                {
                    var tempCrystal = new Crystal((a, b, c, alpha, beta, gamma), null, spaceGroupSeriesNum,
                        ChemicalFormula.Split('=', true)[1].Trim() + "-" + (n++).ToString(), Color.Blue)
                    {
                        Note = $"{wavelength}  {ChemicalFormula}\r\n{str[line]}"
                    };
                    line++;
                    for (; line < str.Length; line++)
                    {
                        if (str[line].Length == 0)
                            break;
                        string[] s = str[line].Split(' ', true);
                        string label = s[0];
                        //ラベル名から元素を決める
                        string temp;

                        int atomicNumber = 0;
                        //for (int q = label.Length; q > 0; q--)
                        for (int q = label.Length; q > 0 && atomicNumber == 0; q--) // 260718Cl: 一致後も外側ループが回り続け、より短いプレフィックス一致 (例: "Si1"→Si の後に S) が上書きするバグを修正 (AMC 経路と同形の最長一致ガード)
                        {
                            temp = label[..q];
                            for (int k = 0; k <= 96; k++)
                            {
                                if (AtomStatic.AtomicName(k) == temp)
                                {
                                    atomicNumber = k;
                                    break;
                                }
                            }
                        }
                        //double x = ZeroToOne(Convert.ToDouble(s[2])); // 261111Cl 変更前 (3 行)
                        //double y = ZeroToOne(Convert.ToDouble(s[3]));
                        //double z = ZeroToOne(Convert.ToDouble(s[4]));
                        double x = ZeroToOne(Convert.ToDouble(s[2], Inv)); // 261111Cl 変更: Invariant で読む
                        double y = ZeroToOne(Convert.ToDouble(s[3], Inv));
                        double z = ZeroToOne(Convert.ToDouble(s[4], Inv));

                        tempCrystal.AddAtoms(label, atomicNumber, 0, 0, null, x, y, z, 1, new DiffuseScatteringFactor(), true);
                    }
                    //SetOpenGL_property(tempCrystal);
                    cry.Add(tempCrystal);
                }
        }
        return [.. cry];
    }

    /// <summary>ConvertFromSMAPから呼び出される</summary>
    /// <param name="x"></param>
    /// <returns></returns>
    public static double ZeroToOne(double x)
    {
        if (x < 0)
            while (x < 0)
                x++;
        else if (x > 1)
            while (x > 1)
                x--;
        return x;
    }

    #endregion SMAPの出力ファイル(*.out)読込

    #region AMCとCIFの読み込みインターフェース ConvertCrystal(filename)  

    public static Crystal2 ConvertToCrystal2(string fileName)
    {
        try
        {
            Crystal2 c2 = null;
            //if (fileName.EndsWith("amc", Ord)) // 261105Cl 変更前 (4 行)
            //    c2 = ConvertFromAmc(fileName);
            //else if (fileName.EndsWith("cif", Ord))
            //    c2 = ConvertFromCIF(fileName);
            // 261105Cl 変更: 拡張子の大文字小文字を区別しない (以前は ".CIF" などを読めなかった。監査 #28)
            if (fileName.EndsWith("amc", StringComparison.OrdinalIgnoreCase))
                c2 = ConvertFromAmc(fileName);
            else if (fileName.EndsWith("cif", StringComparison.OrdinalIgnoreCase))
                c2 = ConvertFromCIF(fileName);

            if (c2 != null)
            {
                c2.fileName = Path.GetFileNameWithoutExtension(fileName);
                return c2;
            }
            else
                return null;
        }
        catch (Exception e)
        {
            if (AssemblyState.IsDebug)
                _ = System.Windows.Forms.MessageBox.Show(fileName + " " + e.Message);
            return null;
        }
    }

    public static Crystal ConvertToCrystal(string fileName)
    {
        try
        {
            //if (fileName.EndsWith("amc", Ord)) // 261105Cl 変更前 (拡張子の判定 2 か所)
            //else if (fileName.EndsWith("cif", Ord))
            if (fileName.EndsWith("amc", StringComparison.OrdinalIgnoreCase)) // 261105Cl 変更: 大文字小文字を区別しない (監査 #28)
            {
                var c = ConvertFromAmc(fileName);
                return c?.ToCrystal();
            }
            else if (fileName.EndsWith("cif", StringComparison.OrdinalIgnoreCase))
            {
                var c = ConvertFromCIF(fileName);
                return c?.ToCrystal();
            }
            else
                return null;
        }
        catch (Exception e)
        {
            if (AssemblyState.IsDebug)
                System.Windows.Forms.MessageBox.Show(e.Message);
            return null;
        }

    }

    /// <summary>260805Cl 追加: CIF テキスト (ファイル内容そのもの) から Crystal を生成する。失敗時は null。ReciPro マクロ Crystal.LoadCifText 用。</summary>
    public static Crystal ConvertToCrystalFromCifText(string cifText)
    {
        try { return ConvertFromCifText(cifText)?.ToCrystal(); }
        catch (Exception e)
        {
            if (AssemblyState.IsDebug)
                System.Windows.Forms.MessageBox.Show(e.Message);
            return null;
        }
    }
    #endregion

    #region amcファイルの読み込み
    private static Crystal2 ConvertFromAmc(string fileName)
    {
        using var reader = new StreamReader(fileName);
        var stringList = new List<string>();
        string strTemp;
        while ((strTemp = reader.ReadLine()) != null)
            if (strTemp != "")
                stringList.Add(strTemp);
        return ConvertFromAmc([.. stringList]);
    }


    /// <summary>amcファイルの読み込み</summary>
    /// <param name="str"></param>
    /// <returns></returns>
    public static Crystal2 ConvertFromAmc(string[] str)
    {
        var n = 0;
        if (str[n].Length == 0)
            n++;

        var Name = str[n];//結晶の名前

        n++; if (str.Length <= n) return null;

        var AuthorName = str[n];//著者の名前
        n++; if (str.Length <= n) return null;
        while ((str[n][^1] < '.' || str[n][^1] > '9') && !str[n].Contains("doi.org") && !str[n].Contains("DOI: ")
            && !str[n].Contains("Mineralogy") )//行の最後の文字が数字ではないとき　(ただし、"doi.org"の文字列が存在するときは除外) 
        {
            AuthorName += ", " + str[n];
            n++; if (str.Length <= n) return null;
        }

        var Reference = str[n];//引用文献
        n++; if (str.Length <= n) return null;

        var Title = str[n];//文献タイトル
        n++; if (str.Length <= n) return null;

        //ここで格子定数、対称性とタイトル
        Crystal2 crystal;
        string extra;
        while (((crystal,extra) = CellParamForAmc(str[n])) == (null,null))
        {
            if (str[n].Length > 0 && char.IsLower(str[n][0]))
                Title += " " + str[n];
            else
                Title += "\r\n" + str[n];
            n++;
            if (str.Length <= n) return null;
        }

        if (Title.Contains("_cod_database_code"))
            Title = Title.Replace("_cod_database_code", "\r\nCOD:");
        else if (Title.Contains("_database_code_amcsd"))
            Title = Title.Replace("_database_code_amcsd", "\r\nAMCSD:");

        n++; if (str.Length <= n) return null;

        double xShift = 0, yShift = 0, zShift = 0;
        if (!str[n].StartsWith("atom", Ord) && str[n].Split(" ", true).Length == 3)
        {
            xShift = ConvertToDouble(str[n].Split(" ", true)[0]);
            yShift = ConvertToDouble(str[n].Split(" ", true)[1]);
            zShift = ConvertToDouble(str[n].Split(" ", true)[2]);
            n++; if (str.Length <= n) return null;
        }

        //2ndセッティングにもかかわらず、shift量がゼロでないときは、1stセッティングに戻す。
        if ((xShift != 0 || yShift != 0 || zShift != 0) && SymmetryStatic.SpaceGroupListWithoutSpace[crystal.sym].EndsWith("(2)", Ord))
            crystal.sym--;


        //ここから原子座標の読み取り
        bool IsOcc = false, IsisoUsed = false, IsanisoUsed = false, IsUtypeUsed = false;
        if (str[n].Contains("occ"))
            IsOcc = true;


        if (str[n].Contains("Uiso") || str[n].Contains("uiso"))
        {
            IsisoUsed = true;
            IsUtypeUsed = true;
        }
        else if (str[n].Contains("Biso") || str[n].Contains("biso"))
        {
            IsisoUsed = true;
            IsUtypeUsed = false;
        }

        if (str[n].Contains("U(1,1)") || str[n].Contains("u(1,1)"))
        {
            IsanisoUsed = true;
            IsUtypeUsed = true;
        }
        else if (str[n].Contains("B(1,1)") || str[n].Contains("b(1,1)"))
        {
            IsanisoUsed = true;
            IsUtypeUsed = false;
        }

        var tempStr = str[n].Split(" ", true);
        var l = new int[tempStr.Length + 1];
        l[0] = 0;

        for (int i = 0; i < tempStr.Length; i++)//各入力値の文字位置を決める。
            l[i + 1] = str[n].IndexOf(tempStr[i]) + tempStr[i].Length;
        for (int i = n + 1; i < str.Length; i++)//最初のatomラベルだけ例外があるようなのでそれに対処
            if (l[1] < str[i].Split(" ", true)[0].Length)
                l[1] = str[i].Split(" ", true)[0].Length;

        // 261105Cl 追加: AMC の異方性の列は U(1,1) U(2,2) U(3,3) U(1,2) U(1,3) U(2,3) (B(i,j) も同じ) の順だが、Atoms2 と ReciPro の並びは
        // (11, 22, 33, 12, 23, 31) (CIF の読込も [u11, u22, u33, u12, u23, u13] に並べ替えている)。以前は列の順のまま入れていたため、
        // AMC 由来の結晶では 23 と 31 (13) が入れ替わっていた (同梱の AMCSD.cdb3 で、単斜晶 b 軸の 2 回軸上の原子が U23≠0・U13=0 になる型が 14,841 原子)
        int column13 = str[n].IndexOf("1,3)", Ord), column23 = str[n].IndexOf("2,3)", Ord);
        bool anisoColumn13Before23 = column13 >= 0 && column23 >= 0 && column13 < column23;

        var atoms = new List<Atoms2>();
        for (int i = n + 1; i < str.Length; i++)
        {//原子座標読み取りループ開始
            str[i] = str[i].PadRight(str[n].Length, ' ');

            var item = str[i].Split(" ", true);
            if (item.Length != l.Length - 1)
            {
                item = new string[l.Length - 1];
                for (int k = 0; k < item.Length; k++)
                {
                    item[k] = str[i][l[k]..l[k + 1]].Trim().TrimEnd();
                    item[k] = item[k].Replace(',', '.');//たまにカンマとピリオドが間違えられている
                    if (item[k].Length > 3 && item[k][1] == ' ')
                    {
                        //二文字目がスペースの場合は、数字がずれている可能性が考えられる。(.123 .456 => ".12", "3 .456" )　
                        //この場合は、二文字目までを削除して対応する
                        item[k] = item[k][2..];
                    }
                }
            }

            int j = 0;
            var label = item[j++];
            var x = item[j++];
            var y = item[j++];
            var z = item[j++];

            if (xShift != 0 || yShift != 0 || zShift != 0)
            {
                //x = (x.ToDouble() + xShift).ToString("f8").TrimEnd(['0']); // 261111Cl 変更前 (3 行)
                //y = (y.ToDouble() + yShift).ToString("f8").TrimEnd(['0']);
                //z = (z.ToDouble() + zShift).ToString("f8").TrimEnd(['0']);
                // 261111Cl 変更: Invariant で書く (小数点がカンマの環境では "0,375" になり、内部形式 (ToBytes) がカンマを表せず座標が壊れていた)
                x = (x.ToDouble() + xShift).ToString("f8", Inv).TrimEnd(['0']);
                y = (y.ToDouble() + yShift).ToString("f8", Inv).TrimEnd(['0']);
                z = (z.ToDouble() + zShift).ToString("f8", Inv).TrimEnd(['0']);
            }

            var occ = "1";
            if (IsOcc)
            {
                occ = item[j++];
                if (occ.Length == 0) occ = "1";
            }

            var iso = "";
            if (IsisoUsed)
                iso = item[j++];

            var aniso = new string[6];
            if (IsanisoUsed)
            {
                for (int k = 0; k < 6; k++)
                    aniso[k] = item[j++];
                if (anisoColumn13Before23) // 261105Cl 追加: (12, 13, 23) → (12, 23, 31)
                    (aniso[4], aniso[5]) = (aniso[5], aniso[4]);
            }
            else
                aniso = null;

            //ラベル名から元素を決める
            var atomicNumber = 0;
            for (int q = label.Length; q > 0 && atomicNumber == 0; q--)
            {
                //var temp = label[..q]; // 261111Cl 変更前 (11 行)
                //for (int k = 0; k <= 96 && atomicNumber == 0; k++)
                //{
                //    if (AtomStatic.AtomicName(k) == temp)
                //    {
                //        atomicNumber = k;
                //        q = -1;
                //        break;
                //    }
                //}
                //if (temp == "OH") //OH基のとき
                //    atomicNumber = -1;
                // 261111Cl 変更: 部分文字列を作らず span 版で 1〜96 を大文字小文字を区別して探す (AtomicName(0) = "" は一致しないので同じ。
                // 一致したら atomicNumber ≠ 0 で外側のループを抜けるので q = -1 は不要)。OH/D/Wat の特例を後で見る順は同じ
                var temp = label.AsSpan(0, q);
                atomicNumber = AtomStatic.AtomicNumber(temp, caseSensitive: true);
                if (temp is "OH") //OH基のとき
                    atomicNumber = -1;
                else if (temp is "D") //重水素のとき
                    atomicNumber = 255;
                else if (temp is "Wat" or "WAT" or "wat") //水のとき
                    atomicNumber = -2;
            }

            var IsIso = aniso == null || aniso.All(e => e.Length == 0);

            if (atomicNumber > 0)
            {
                atoms.Add(new Atoms2(label, (byte)atomicNumber, 0, 0, [x, y, z], occ, IsIso, IsUtypeUsed, iso, aniso));
            }
            else if (atomicNumber == -1)//"OH"のときの対処
            {
                atoms.Add(new Atoms2(label, 1, 0, 0, [x, y, z], occ, IsIso, IsUtypeUsed, iso, aniso));
                atoms.Add(new Atoms2(label, 8, 0, 0, [x, y, z], occ, IsIso, IsUtypeUsed, iso, aniso));
            }
            else if (atomicNumber == -2)//"Wat"水のときの対処
            {
                atoms.Add(new Atoms2(label, 1, 0, 0, [x, y, z], occ, IsIso, IsUtypeUsed, iso, aniso));
                atoms.Add(new Atoms2(label, 1, 0, 0, [x, y, z], occ, IsIso, IsUtypeUsed, iso, aniso));
                atoms.Add(new Atoms2(label, 8, 0, 0, [x, y, z], occ, IsIso, IsUtypeUsed, iso, aniso));
            }
        }
        crystal.name = Name;
        crystal.sect = Title;
        crystal.auth = AuthorName;
        crystal.jour = Reference;
        crystal.atoms = atoms;

        //最後に,extra (ファイル中の空間格子と、ソフト内の空間格子の記号が違う場合のみ)を処理
        if (extra != null)
        {
            var c = crystal.ToCrystal();
            (double X, double Y, double Z)[] translation = extra switch
            {
                "A" => [(0, 0, 0), (0.0, 0.5, 0.5)],
                "B" => [(0, 0, 0), (0.5, 0.0, 0.5)],
                "C" => [(0, 0, 0), (0.5, 0.5, 0)],
                "I" => [(0, 0, 0), (0.5, 0.5, 0.5)],
                "F" => [(0, 0, 0), (0.0, 0.5, 0.5), (0.5, 0.0, 0.5), (0.5, 0.5, 0)],
                _ => [(0, 0, 0)],
            };

            var atomsEx = new List<Atoms>();
            foreach (var a in c.Atoms)
            {
                foreach (var (X, Y, Z) in translation) 
                    atomsEx.Add(new Atoms(a.Label,a.AtomicNumber,a.SubNumberXray,a.SubNumberElectron,a.Isotope,a.SymmetrySeriesNumber,
                        new Vector3DBase( a.X+X,a.Y+Y,a.Z+Z ), a.PositionError,a.Occ,a.Occ_err,a.Dsf,a.Material,a.Radius,true,false));
            }
            c.Atoms = [.. atomsEx];
            crystal = c.ToCrystal2();
        }



        return crystal;
    }

    private static (Crystal2, string) CellParamForAmc(string str)
    {
        string extra = null;

        double A, B, C, Alpha, Beta, Gamma;
        int symmetrySeriesNumber = -1;
        var s = str.Split(" ", true);
        if (s.Length != 7)
            return (null,null);
        try
        {
            for (int i = 0; i < 6; i++) if (s[i].EndsWith(',')) s[i] = s[i].TrimEnd(','); //最後に','が入っているときは削除

            A = s[0].ToDouble();
            B = s[1].ToDouble();
            C = s[2].ToDouble();
            Alpha = s[3].ToDouble();
            Beta = s[4].ToDouble();
            Gamma = s[5].ToDouble();
            if (double.IsNaN(A) || double.IsNaN(B) || double.IsNaN(C) || double.IsNaN(Alpha) || double.IsNaN(Beta) || double.IsNaN(Gamma))
                return (null, null);
        }
        catch { return (null, null); }
        string SgName = s[6];

        SgName = SgName.Replace("_", "sub");
        bool isAsterisk = SgName.Contains('*');
        SgName = SgName.Replace("*", "");

        if (SgName.Contains(":2"))
            isAsterisk = true;
        SgName = SgName.Replace(":1", "");
        SgName = SgName.Replace(":2", "");

        #region 空間群の場合分け

        if (SgName == "Pm3") SgName = "Pm-3";
        else if (SgName == "Pn3") SgName = "Pn-3";
        else if (SgName == "Fm3") SgName = "Fm-3";
        else if (SgName == "Fd3") SgName = "Fd-3";
        else if (SgName == "Im3") SgName = "Im-3";
        else if (SgName == "Pa3") SgName = "Pa-3";
        else if (SgName == "Ia3") SgName = "Ia-3";

        else if (SgName == "Pm3m") SgName = "Pm-3m";
        //else if (SgName == "Pn3n") SgName = "Pm-3m";
        else if (SgName == "Pn3n") SgName = "Pn-3n"; // 260718Cl: Pn3n (SG222) が Pm-3m (SG221) に誤マップされていたコピペバグを修正 (前後の行は全て文字保存の変換)
        else if (SgName == "Pm3n") SgName = "Pm-3n";
        else if (SgName == "Pn3m") SgName = "Pn-3m";
        else if (SgName == "Fm3m") SgName = "Fm-3m";
        else if (SgName == "Fm3c") SgName = "Fm-3c";
        else if (SgName == "Fd3m") SgName = "Fd-3m";
        else if (SgName == "Fd3c") SgName = "Fd-3c";
        else if (SgName == "Im3m") SgName = "Im-3m";
        else if (SgName == "Ia3d") SgName = "Ia-3d";
        else if (SgName == "P4/m-32/m") SgName = "Pm-3m";

        else if (SgName == "I2sub1/a-3") SgName = "Ia-3";

        else if (SgName == "R-32/c") SgName = "R-3c";

        else if (SgName == "P2" && Alpha == 90.0 && Gamma == 90.0) SgName = "P121";
        else if (SgName == "P2" && Alpha == 90.0 && Beta == 90.0) SgName = "P112";
        else if (SgName == "P2" && Beta == 90.0 && Gamma == 90.0) SgName = "P211";
        else if (SgName == "P2sub1" && Alpha == 90.0 && Gamma == 90.0) SgName = "P12sub11";
        else if (SgName == "P2sub1" && Alpha == 90.0 && Beta == 90.0) SgName = "P112sub1";
        else if (SgName == "P2sub1" && Beta == 90.0 && Gamma == 90.0) SgName = "P2sub111";
        else if (SgName == "C2" && Alpha == 90.0 && Gamma == 90.0) SgName = "C121";
        else if (SgName == "A2" && Alpha == 90.0 && Gamma == 90.0) SgName = "A121";
        else if (SgName == "I2" && Alpha == 90.0 && Gamma == 90.0) SgName = "I121";
        else if (SgName == "A2" && Alpha == 90.0 && Beta == 90.0) SgName = "A112";
        else if (SgName == "B2" && Alpha == 90.0 && Beta == 90.0) SgName = "B112";
        else if (SgName == "I2" && Alpha == 90.0 && Beta == 90.0) SgName = "I112";
        else if (SgName == "B2" && Beta == 90.0 && Gamma == 90.0) SgName = "B211";
        else if (SgName == "C2" && Beta == 90.0 && Gamma == 90.0) SgName = "C211";
        else if (SgName == "I2" && Beta == 90.0 && Gamma == 90.0) SgName = "I211";
        else if (SgName == "Pm" && Alpha == 90.0 && Gamma == 90.0) SgName = "P1m1";
        else if (SgName == "Pm" && Alpha == 90.0 && Beta == 90.0) SgName = "P11m";
        else if (SgName == "Pm" && Beta == 90.0 && Gamma == 90.0) SgName = "Pm11";
        else if (SgName == "Pc" && Alpha == 90.0 && Gamma == 90.0) SgName = "P1c1";
        else if (SgName == "Pn" && Alpha == 90.0 && Gamma == 90.0) SgName = "P1n1";
        else if (SgName == "Pa" && Alpha == 90.0 && Gamma == 90.0) SgName = "P1a1";
        else if (SgName == "Pa" && Alpha == 90.0 && Beta == 90.0) SgName = "P11a";
        else if (SgName == "Pn" && Alpha == 90.0 && Beta == 90.0) SgName = "P11n";
        else if (SgName == "Pb" && Alpha == 90.0 && Beta == 90.0) SgName = "P11b";
        else if (SgName == "Pb" && Beta == 90.0 && Gamma == 90.0) SgName = "Pb11";
        else if (SgName == "Pn" && Beta == 90.0 && Gamma == 90.0) SgName = "Pn11";
        else if (SgName == "Pc" && Beta == 90.0 && Gamma == 90.0) SgName = "Pc11";
        else if (SgName == "Cm" && Alpha == 90.0 && Gamma == 90.0) SgName = "C1m1";
        else if (SgName == "Am" && Alpha == 90.0 && Gamma == 90.0) SgName = "A1m1";
        else if (SgName == "Im" && Alpha == 90.0 && Gamma == 90.0) SgName = "I1m1";
        else if (SgName == "Am" && Alpha == 90.0 && Beta == 90.0) SgName = "A11m";
        else if (SgName == "Bm" && Alpha == 90.0 && Beta == 90.0) SgName = "B11m";
        else if (SgName == "Im" && Alpha == 90.0 && Beta == 90.0) SgName = "I11m";
        else if (SgName == "Bm" && Beta == 90.0 && Gamma == 90.0) SgName = "Bm11";
        else if (SgName == "Cm" && Beta == 90.0 && Gamma == 90.0) SgName = "Cm11";
        else if (SgName == "Im" && Beta == 90.0 && Gamma == 90.0) SgName = "Im11";
        else if (SgName == "Cc" && Alpha == 90.0 && Gamma == 90.0) SgName = "C1c1";
        else if (SgName == "An" && Alpha == 90.0 && Gamma == 90.0) SgName = "A1n1";
        else if (SgName == "Ia" && Alpha == 90.0 && Gamma == 90.0) SgName = "I1a1";
        else if (SgName == "Aa" && Alpha == 90.0 && Gamma == 90.0) SgName = "A1a1";
        else if (SgName == "Cn" && Alpha == 90.0 && Gamma == 90.0) SgName = "C1n1";
        else if (SgName == "Ic" && Alpha == 90.0 && Gamma == 90.0) SgName = "I1c1";
        else if (SgName == "Aa" && Alpha == 90.0 && Beta == 90.0) SgName = "A11a";
        else if (SgName == "Bn" && Alpha == 90.0 && Beta == 90.0) SgName = "B11n";
        else if (SgName == "Ib" && Alpha == 90.0 && Beta == 90.0) SgName = "I11b";
        else if (SgName == "Bb" && Alpha == 90.0 && Beta == 90.0) SgName = "B11b";
        else if (SgName == "An" && Alpha == 90.0 && Beta == 90.0) SgName = "A11n";
        else if (SgName == "Ia" && Alpha == 90.0 && Beta == 90.0) SgName = "I11a";
        else if (SgName == "Bb" && Beta == 90.0 && Gamma == 90.0) SgName = "Bb11";
        else if (SgName == "Cn" && Beta == 90.0 && Gamma == 90.0) SgName = "Cn11";
        else if (SgName == "Ic" && Beta == 90.0 && Gamma == 90.0) SgName = "Ic11";
        else if (SgName == "Cc" && Beta == 90.0 && Gamma == 90.0) SgName = "Cc11";
        else if (SgName == "Bn" && Beta == 90.0 && Gamma == 90.0) SgName = "Bn11";
        else if (SgName == "Ib" && Beta == 90.0 && Gamma == 90.0) SgName = "Ib11";
        else if (SgName == "P2/m" && Alpha == 90.0 && Gamma == 90.0) SgName = "P12/m1";
        else if (SgName == "P2/m" && Alpha == 90.0 && Beta == 90.0) SgName = "P112/m";
        else if (SgName == "P2/m" && Beta == 90.0 && Gamma == 90.0) SgName = "P2/m11";
        else if (SgName == "P2sub/m" && Alpha == 90.0 && Gamma == 90.0) SgName = "P12sub1/m1";
        else if (SgName == "P2sub/m" && Alpha == 90.0 && Beta == 90.0) SgName = "P112sub1/m";
        else if (SgName == "P2sub/m" && Beta == 90.0 && Gamma == 90.0) SgName = "P2sub1/m11";
        else if (SgName == "C2/m" && Alpha == 90.0 && Gamma == 90.0) SgName = "C12/m1";
        else if (SgName == "A2/m" && Alpha == 90.0 && Gamma == 90.0) SgName = "A12/m1";
        else if (SgName == "I2/m" && Alpha == 90.0 && Gamma == 90.0) SgName = "I12/m1";
        else if (SgName == "A2/m" && Alpha == 90.0 && Beta == 90.0) SgName = "A112/m";
        else if (SgName == "B2/m" && Alpha == 90.0 && Beta == 90.0) SgName = "B112/m";
        else if (SgName == "I2/m" && Alpha == 90.0 && Beta == 90.0) SgName = "I112/m";
        else if (SgName == "B2/m" && Beta == 90.0 && Gamma == 90.0) SgName = "B2/m11";
        else if (SgName == "C2/m" && Beta == 90.0 && Gamma == 90.0) SgName = "C2/m11";
        else if (SgName == "I2/m" && Beta == 90.0 && Gamma == 90.0) SgName = "I2/m11";
        else if (SgName == "P2/c" && Alpha == 90.0 && Gamma == 90.0) SgName = "P12/c1";
        else if (SgName == "P2/n" && Alpha == 90.0 && Gamma == 90.0) SgName = "P12/n1";
        else if (SgName == "P2/a" && Alpha == 90.0 && Gamma == 90.0) SgName = "P12/a1";
        else if (SgName == "P2/a" && Alpha == 90.0 && Beta == 90.0) SgName = "P112/a";
        else if (SgName == "P2/n" && Alpha == 90.0 && Beta == 90.0) SgName = "P112/n";
        else if (SgName == "P2/b" && Alpha == 90.0 && Beta == 90.0) SgName = "P112/b";
        else if (SgName == "P2/b" && Beta == 90.0 && Gamma == 90.0) SgName = "P2/b11";
        else if (SgName == "P2/n" && Beta == 90.0 && Gamma == 90.0) SgName = "P2/n11";
        else if (SgName == "P2/c" && Beta == 90.0 && Gamma == 90.0) SgName = "P2/c11";
        else if (SgName == "P2sub1/c" && Alpha == 90.0 && Gamma == 90.0) SgName = "P12sub1/c1";
        else if (SgName == "P2sub1/n" && Alpha == 90.0 && Gamma == 90.0) SgName = "P12sub1/n1";
        else if (SgName == "P2sub1/a" && Alpha == 90.0 && Gamma == 90.0) SgName = "P12sub1/a1";
        else if (SgName == "P2sub1/a" && Alpha == 90.0 && Beta == 90.0) SgName = "P112sub1/a";
        else if (SgName == "P2sub1/n" && Alpha == 90.0 && Beta == 90.0) SgName = "P112sub1/n";
        else if (SgName == "P2sub1/b" && Alpha == 90.0 && Beta == 90.0) SgName = "P112sub1/b";
        else if (SgName == "P2sub1/b" && Beta == 90.0 && Gamma == 90.0) SgName = "P2sub1/b11";
        else if (SgName == "P2sub1/n" && Beta == 90.0 && Gamma == 90.0) SgName = "P2sub1/n11";
        else if (SgName == "P2sub1/c" && Beta == 90.0 && Gamma == 90.0) SgName = "P2sub1/c11";
        else if (SgName == "C2/c" && Alpha == 90.0 && Gamma == 90.0) SgName = "C12/c1";
        else if (SgName == "A2/n" && Alpha == 90.0 && Gamma == 90.0) SgName = "A12/n1";
        else if (SgName == "I2/a" && Alpha == 90.0 && Gamma == 90.0) SgName = "I12/a1";
        else if (SgName == "A2/a" && Alpha == 90.0 && Gamma == 90.0) SgName = "A12/a1";
        else if (SgName == "C2/n" && Alpha == 90.0 && Gamma == 90.0) SgName = "C12/n1";
        else if (SgName == "I2/c" && Alpha == 90.0 && Gamma == 90.0) SgName = "I12/c1";
        else if (SgName == "A2/a" && Alpha == 90.0 && Beta == 90.0) SgName = "A112/a";
        else if (SgName == "B2/n" && Alpha == 90.0 && Beta == 90.0) SgName = "B112/n";
        else if (SgName == "I2/b" && Alpha == 90.0 && Beta == 90.0) SgName = "I112/b";
        else if (SgName == "B2/b" && Alpha == 90.0 && Beta == 90.0) SgName = "B112/b";
        else if (SgName == "A2/n" && Alpha == 90.0 && Beta == 90.0) SgName = "A112/n";
        else if (SgName == "I2/a" && Alpha == 90.0 && Beta == 90.0) SgName = "I112/a";
        else if (SgName == "B2/b" && Beta == 90.0 && Gamma == 90.0) SgName = "B2/b11";
        else if (SgName == "C2/n" && Beta == 90.0 && Gamma == 90.0) SgName = "C2/n11";
        else if (SgName == "I2/c" && Beta == 90.0 && Gamma == 90.0) SgName = "I2/c11";
        else if (SgName == "C2/c" && Beta == 90.0 && Gamma == 90.0) SgName = "C2/c11";
        else if (SgName == "B2/n" && Beta == 90.0 && Gamma == 90.0) SgName = "B2/n11";
        else if (SgName == "I2/b" && Beta == 90.0 && Gamma == 90.0) SgName = "I2/b11";

        else if (SgName == "C2sub1" && Alpha == 90.0 && Gamma == 90.0) { SgName = "P12sub11"; extra = "C"; }
        else if (SgName == "C2sub1" && Alpha == 90.0 && Beta == 90.0) { SgName = "P112sub1"; extra = "C"; }
        else if (SgName == "C2sub1" && Beta == 90.0 && Gamma == 90.0) { SgName = "P2sub111"; extra = "C"; }
        else if (SgName == "C2/a" && Alpha == 90.0 && Gamma == 90.0) { SgName = "P12/a1"; extra = "C"; }
        else if (SgName == "C2/a" && Alpha == 90.0 && Beta == 90.0) { SgName = "P112/a"; extra = "C"; }
        else if (SgName == "C2sub1/a" && Alpha == 90.0 && Gamma == 90.0) { SgName = "P12sub1/a1"; extra = "C"; }
        else if (SgName == "C2sub1/a" && Alpha == 90.0 && Beta == 90.0) { SgName = "P112sub1/a"; extra = "C"; }

        else if (SgName == "F2/m" && Alpha == 90.0 && Gamma == 90.0) { SgName = "P12/m1"; extra = "F"; }
        else if (SgName == "F2/m" && Beta == 90.0 && Gamma == 90.0) { SgName = "P2/m11"; extra = "F"; }
        else if (SgName == "F2/m" && Alpha == 90.0 && Beta == 90.0) { SgName = "P112/m"; extra = "F"; }

        else if (SgName == "B2sub1" && Alpha == 90.0 && Gamma == 90.0) { SgName = "P12sub11"; extra = "B"; }
        else if (SgName == "B2sub1" && Alpha == 90.0 && Beta == 90.0) { SgName = "P112sub1"; extra = "B"; }
        else if (SgName == "B2sub1" && Beta == 90.0 && Gamma == 90.0) { SgName = "P2sub111"; extra = "B"; }
        else if (SgName == "B2sub1/m" && Alpha == 90.0 && Gamma == 90.0) { SgName = "P12sub1/m1"; extra = "B"; }
        else if (SgName == "B2sub1/m" && Alpha == 90.0 && Beta == 90.0) { SgName = "P112sub1/m"; extra = "B"; }
        else if (SgName == "B2sub1/m" && Beta == 90.0 && Gamma == 90.0) { SgName = "P2sub1/m11"; extra = "B"; }
        
        else if (SgName == "Imab" && Beta == 90.0 && Gamma == 90.0) { SgName = "Pmab"; extra = "I"; }

        else if (SgName == "F4/mmm") { SgName = "P4/mmm"; extra = "F"; }

        #endregion amcファイルの読み込み

        //文字列を含んでいて、かつ、文字数が少ない空間群を選択する (C1とC121などを見分けるため)
        int length = int.MaxValue;

        var sg = SymmetryStatic.SpaceGroupListWithoutSpace;
        for (int i = 0; i < SymmetryStatic.TotalSpaceGroupNumber; i++)
        {
            if (sg[i].Length < length)
            {
                var sg_list = new List<string>([sg[i]]);

                if (sg[i].Contains('e'))
                {
                    sg_list.Add(sg[i].Replace('e', 'a'));
                    sg_list.Add(sg[i].Replace('e', 'b'));
                    sg_list.Add(sg[i].Replace('e', 'c'));
                }

                if (sg_list.Any(s => s.IndexOf(SgName, 0) != -1))
                {
                    symmetrySeriesNumber = i;
                    length = sg[i].Length;
                }
            }
        }
        if (symmetrySeriesNumber == -1)
            return (null, null);
        //Rhombohedoralのときの処置
        if (A == B && B == C && Alpha == Beta && Beta == Gamma && SymmetryStatic.Symmetries[symmetrySeriesNumber].SpaceGroupHMStr.Contains("Hex", StringComparison.CurrentCulture))
            symmetrySeriesNumber++;

        //Asteriskの時(2nd setting)の処理
        if (isAsterisk && sg[symmetrySeriesNumber].EndsWith("(1)", Ord))
            symmetrySeriesNumber++;


        if (symmetrySeriesNumber >= 0)
        {
            var r = new Random();
            return (new Crystal2
            {

                CellTexts = [s[0], s[1], s[2], s[3], s[4], s[5]],
                sym = (short)symmetrySeriesNumber,
                argb = Color.FromArgb(r.Next(255), r.Next(255), r.Next(255)).ToArgb()
            },extra);
        }
        else
            return (null, null);
    }

    private static double ConvertToDouble(string str, bool IsHex = false)
    {
        try
        {
            if (str.Trim().Length == 0)
                return 0;
            else if (str.IndexOf('/') > 0)
            {
                var temp = str.Split('/', true);
                //return temp.Length == 2 ? Convert.ToDouble(temp[0]) / Convert.ToDouble(temp[1]) : 0; // 261111Cl 変更前
                return temp.Length == 2 ? Convert.ToDouble(temp[0], Inv) / Convert.ToDouble(temp[1], Inv) : 0; // 261111Cl 変更: Invariant で読む
            }
            else if (IsHex)
            {
                if (str.Contains(".0833")) return 1.0 / 12.0;
                else if (str.Contains(".1667")) return 1.0 / 6.0;
                else if (str.Contains(".16667")) return 1.0 / 6.0;
                else if (str.Contains(".166667")) return 1.0 / 6.0;
                else if (str.Contains(".333333")) return 1.0 / 3.0;
                else if (str.Contains(".33333")) return 1.0 / 3.0;
                else if (str.Contains(".3333")) return 1.0 / 3.0;
                else if (str.Contains(".4167")) return 5.0 / 12.0;
                else if (str.Contains(".5833")) return 7.0 / 12.0;
                else if (str.Contains(".666667")) return 2.0 / 3.0;
                else if (str.Contains(".66667")) return 2.0 / 3.0;
                else if (str.Contains(".6667")) return 2.0 / 3.0;
                else if (str.Contains(".8333")) return 5.0 / 6.0;
                else if (str.Contains(".9167")) return 11.0 / 12.0;

                if (str.Contains(",0833")) return 1.0 / 12.0;
                else if (str.Contains(",1667")) return 1.0 / 6.0;
                else if (str.Contains(",16667")) return 1.0 / 6.0;
                else if (str.Contains(",166667")) return 1.0 / 6.0;
                else if (str.Contains(",333333")) return 1.0 / 3.0;
                else if (str.Contains(",33333")) return 1.0 / 3.0;
                else if (str.Contains(",3333")) return 1.0 / 3.0;
                else if (str.Contains(",4167")) return 5.0 / 12.0;
                else if (str.Contains(",5833")) return 7.0 / 12.0;
                else if (str.Contains(",666667")) return 2.0 / 3.0;
                else if (str.Contains(",66667")) return 2.0 / 3.0;
                else if (str.Contains(",6667")) return 2.0 / 3.0;
                else if (str.Contains(",8333")) return 5.0 / 6.0;
                else if (str.Contains(",9167")) return 11.0 / 12.0;
            }

            //return Convert.ToDouble(str); // 261111Cl 変更前
            return Convert.ToDouble(str, Inv); // 261111Cl 変更: Invariant で読む (小数点がカンマの環境で AMC の原点シフト "0.25" が 25 になっていた)
        }
        catch
        {
            return 0;
        }
    }

    #endregion

    #region CIFファイルの読み込み
    static readonly Random rnd = Random.Shared;

    static readonly FrozenSet<string> ignoreWords1 = ["_shelx_hkl_file", "_shelxl_hkl_file", "_shelx_fab_file", "_shelx_res_file", "_shelxl_res_file", "_iucr_refine_reflections_details"];//語尾に何もつかない

    static readonly FrozenSet<string> ignoreWords2 = ["_shelx_hkl_", "_shelx_fab_", "_shelx_res_"];//語尾に fileかchecksumがつく

    static readonly FrozenSet<string> ignoreWords3 = ["_refln", "_geom", "_platon", "_diffrn_refln"];
    private static Crystal2 ConvertFromCIF(string fileName)
    {
        //260805Cl 変更: 本文をテキスト版 ConvertFromCifText へ抽出 (ReciPro マクロ Crystal.LoadCifText 用)。処理内容は不変
        using var reader = new StreamReader(fileName);
        return ConvertFromCifText(reader.ReadToEnd());
    }

    /// <summary>260805Cl 追加: CIF ファイルの内容 (テキストそのもの) から Crystal2 を生成する。ReciPro マクロ Crystal.LoadCifText 用の入口。</summary>
    public static Crystal2 ConvertFromCifText(string cifText)
    {
        string[] str;
        {
            var strTemp = cifText;
            Replace(ref strTemp, "\r\n", "\n");
            Replace(ref strTemp, "\r", "\n");

            EscapeString(ref strTemp);

            foreach (var word in ignoreWords1 )
            {
                if (strTemp.Contains(word + "\n;"))
                {
                    int start = strTemp.IndexOf(word + "\n;");
                    int end = strTemp.IndexOf("\n;", start + word.Length + 3);
                    if (end > start && end + 3 < strTemp.Length)
                        strTemp = string.Concat(strTemp[..start], strTemp[(end + 3)..]);
                }
            }

            str = strTemp.Split('\n', true);
        }

        // 261111Cl 変更: 行の削除は List<string> の RemoveRange (その場で詰める) で行う。旧は消すたびに配列全体を作り直していた (O(N²))。
        // 探す条件・消す範囲・順序は旧と同じ。最後に配列へ戻す
        var lines = new List<string>(str); // 261111Cl 追加
        foreach (var word in ignoreWords2)
        {
            int start = -1, end = -1;
            //while ((start = Array.IndexOf( str, word + "file")) > -1 && // 261111Cl 変更前 (3 行)
            //    (end =Array.FindIndex(str, s => s.StartsWith(word + "checksum", Ord))) > -1)
            //    RemoveRange(ref str, Math.Min(start, end), Math.Abs(start - end) + 1);
            while ((start = lines.IndexOf(word + "file")) > -1 &&
                (end = lines.FindIndex(s => s.StartsWith(word + "checksum", Ord))) > -1)
                lines.RemoveRange(Math.Min(start, end), Math.Abs(start - end) + 1);
        }


        foreach (var word in ignoreWords3)
        {
            int start = -1;
            int from = 0; // 261111Cl 追加: 次に探し始める行。前回見つけた行の 1 つ前 (start) より前の行は一致せず、消されもしないので、そこから探しても旧 (毎回先頭から) と同じ
            //while ((start = Array.FindIndex(str, s => s.Trim().StartsWith(word, Ord)) - 1) > -1) // 261111Cl 変更前
            while ((start = lines.FindIndex(from, s => s.Trim().StartsWith(word, Ord)) - 1) > -1)
            {
                from = start; // 261111Cl 追加
                //if (str[start] == "loop_") // 261111Cl 変更前 (この中の str[…] を lines[…] に、RemoveRange を lines.RemoveRange にした)
                if (lines[start] == "loop_")
                {
                    var end = start + 1;
                    //while (str[end].StartsWith(word, Ord)) // 261111Cl 変更前
                    while ((end < lines.Count ? lines[end] : throw new IndexOutOfRangeException()).StartsWith(word, Ord)) // 261111Cl 変更: 末尾を越えたときの例外も旧 (配列) と同じ型にする
                        end++;

                    //for (; end < str.Length; end++) // 261111Cl 変更前 (4 行)
                    //    if (str[end] == "loop_" || str[end].StartsWith("_", Ord) || str[end].StartsWith("#", Ord))
                    //        break;
                    //RemoveRange(ref str, start, end - start);
                    for (; end < lines.Count; end++)
                        if (lines[end] == "loop_" || lines[end].StartsWith("_", Ord) || lines[end].StartsWith("#", Ord))
                            break;
                    lines.RemoveRange(start, end - start);
                }
                else
                {
                    //RemoveAt(ref str, start + 1); // 261105Cl 変更前
                    // 261105Cl 変更: loop_ の外の単独のタグは、値もいっしょに消す。値が ; で囲まれた複数行のときは閉じの ; の行まで。
                    // 以前はタグの行だけを消して ; の文章を残したため、文章が直前のループの行として読まれ、ループ全体が読めなくなっていた。
                    // SHELXL の CIF では _geom_special_details が _atom_site_aniso のループの直後に来るので、異方性の温度因子が全部捨てられていた
                    // (tools/EbsdCheck/testdata/koseki_cif/888769.cif など。tools/AdpCheck realcif)
                    var tag = start + 1;
                    var last = tag;
                    // 値はタグの次の「コメントでない行」にある。タグと値の間のコメント行もいっしょに消す (Codex2 r3。以前はタグとコメントだけを消して値を残していた)
                    var next = tag + 1;
                    //while (next < str.Length && str[next].TrimStart().StartsWith('#')) // 261111Cl 変更前 (12 行。str → lines、RemoveRange → lines.RemoveRange)
                    //    next++;
                    //if (next < str.Length && str[next].StartsWith(';'))
                    //{
                    //    last = next + 1;
                    //    while (last < str.Length - 1 && !str[last].StartsWith(';'))
                    //        last++;
                    //}
                    //else if (next < str.Length && str[tag].Trim().Split((char[])null, StringSplitOptions.RemoveEmptyEntries).Length == 1
                    //    && !str[next].TrimStart().StartsWith('_') && !str[next].TrimStart().StartsWith("loop_", Ord) && !str[next].StartsWith("data_", Ord))
                    //    last = next; // タグだけの行で、値が次の (コメントでない) 行にあるとき
                    //RemoveRange(ref str, tag, last - tag + 1);
                    while (next < lines.Count && lines[next].TrimStart().StartsWith('#'))
                        next++;
                    if (next < lines.Count && lines[next].StartsWith(';'))
                    {
                        last = next + 1;
                        while (last < lines.Count - 1 && !lines[last].StartsWith(';'))
                            last++;
                    }
                    else if (next < lines.Count && lines[tag].Trim().Split((char[])null, StringSplitOptions.RemoveEmptyEntries).Length == 1
                        && !lines[next].TrimStart().StartsWith('_') && !lines[next].TrimStart().StartsWith("loop_", Ord) && !lines[next].StartsWith("data_", Ord))
                        last = next; // タグだけの行で、値が次の (コメントでない) 行にあるとき
                    lines.RemoveRange(tag, last - tag + 1);
                }
            }
        }

        //return ConvertFromCIF(str); // 261111Cl 変更前
        return ConvertFromCIF([.. lines]); // 261111Cl 変更
    }

    /// <summary>CIFファイルを読み込む</summary>
    /// <param name="str"></param>
    /// <returns></returns>
    private static Crystal2 ConvertFromCIF(string[] str)
    {
        //まず ;と; で囲まれている複数にわたる行を一行にする
        //var str = new List<string>();
        var note = "";
        if (str[0].StartsWith("data", Ord))
            note = str[0];
        //for (int n = 0; n < str.Length; n++) // 261105Cl 変更前 (5 行)
        //{
        //    //全てのコメントを削除
        //    while (str[n].StartsWith("#", Ord) && n < str.Length - 1)
        //        RemoveAt(ref str, n);
        //}
        // 261105Cl 変更: 全てのコメントを削除 (最後の行も)。以前は n < str.Length - 1 の条件で最後の行のコメントだけ残り、
        // 前処理で末尾の _geom のループが除かれたあと最後に残る「# 10. MOLECULAR GEOMETRY」のような区切り (CCDC の CIF に多い) が、
        // 直前のループ (_atom_site_aniso) のデータとして読まれ、異方性の温度因子が全部捨てられていた (tools/EbsdCheck/testdata/koseki_cif/818186_OAc.cif)
        str = [.. str.Where(s => !s.StartsWith('#'))];
        #region 261111Cl 変更前: ; で囲まれた複数行を 1 行にするたびに RemoveRange で配列全体を作り直していた
        //for (int n = 0; n < str.Length; n++)
        //{
        //    if (str[n].StartsWith(";", Ord))//;で始まる行を見つけたら
        //    {
        //        //次の;の位置を見つけて mに格納
        //        int m = n + 1;
        //        for (; m < str.Length; m++)
        //            if (str[m].StartsWith(';'))
        //                break;
        //        //直後の行が「loop_」の場合や「_」で始まる場合は、;を一文字だけ消す
        //        if (n + 1 < str.Length && (str[n + 1].StartsWith("loop_", Ord) || str[n + 1].StartsWith("_", Ord)))
        //            str[n] = "";
        //        //次の;が見つからない場合も、;を一文字だけ消す
        //        else if (m == str.Length)
        //            str[n] = "";
        //        else
        //        {
        //            var temp = new StringBuilder();
        //            //次に;が出てくるところまですすめてまとめて一行にする
        //            for (int i = n + 1; i < m; i++)
        //                temp.Append(str[i] + " ");
        //            str[n] = $"'{temp.ToString().Trim().TrimEnd()}'";
        //            RemoveRange(ref str, n + 1, m - n);
        //        }
        //    }
        //}
        #endregion
        // 261111Cl 変更: 1 パスで新しい List に書き出す。まとめた行の中身 (n+1〜m 行目) は読み飛ばす (旧の RemoveRange(n + 1, m - n) と同じ行)。
        // 旧で消した行はどれも処理済みの行なので、後ろの行から見た「次の行」「次の ;」は同じ
        {
            var merged = new List<string>(str.Length);
            for (int n = 0; n < str.Length; n++)
            {
                if (str[n].StartsWith(";", Ord))//;で始まる行を見つけたら
                {
                    //次の;の位置を見つけて mに格納
                    int m = n + 1;
                    for (; m < str.Length; m++)
                        if (str[m].StartsWith(';'))
                            break;
                    //直後の行が「loop_」の場合や「_」で始まる場合は、;を一文字だけ消す
                    if (n + 1 < str.Length && (str[n + 1].StartsWith("loop_", Ord) || str[n + 1].StartsWith("_", Ord)))
                        merged.Add("");
                    //次の;が見つからない場合も、;を一文字だけ消す
                    else if (m == str.Length)
                        merged.Add("");
                    else
                    {
                        var temp = new StringBuilder();
                        //次に;が出てくるところまですすめてまとめて一行にする
                        for (int i = n + 1; i < m; i++)
                            temp.Append(str[i] + " ");
                        merged.Add($"'{temp.ToString().Trim().TrimEnd()}'");
                        n = m;
                    }
                }
                else
                    merged.Add(str[n]);
            }
            str = [.. merged];
        }

        if (str[^1]=="")
            str[^1] = "#End of data";

        //for (int n = 0; n <str.Length; n++) // 261111Cl 変更前 (7 行): 空行を 1 行消すたびに配列全体を作り直していた
        //{
        //    Replace(ref str[n], "\t", " ");
        //    TrimStartEnd(ref str[n]);
        //    if (str[n] == "")
        //        RemoveAt(ref str, n--);
        //}
        // 261111Cl 変更: 1 パスのフィルタ (タブを空白に → 前後の空白を除く → 空なら捨てる。順は同じ)
        {
            var kept = new List<string>(str.Length);
            foreach (var line in str)
            {
                var s = line;
                Replace(ref s, "\t", " ");
                TrimStartEnd(ref s);
                if (s != "")
                    kept.Add(s);
            }
            str = [.. kept];
        }

        //次に'あるいは"で囲まれている文字列中の空白を偶然出てこないような文字列に変換する
        for (int n = 0; n < str.Length; n++)
        {
            if (str[n].Contains("''"))//''という文字列が含まれていたら
            {
                //var temp = str[n][str[n].IndexOf("'")..]; // 261111Cl 変更前
                var temp = str[n][str[n].IndexOf('\'')..]; // 261111Cl 変更: 文字で探す (Ordinal。旧の文字列版は現在のカルチャで比べていた)
                Replace(ref temp, "''", "薔");
                //str[n] = $"{str[n][..str[n].IndexOf("'")]}{temp}"; // 261111Cl 変更前
                str[n] = $"{str[n][..str[n].IndexOf('\'')]}{temp}"; // 261111Cl 変更
            }

            if (str[n].Contains('\''))//'が含まれていたら
            {
                while (str[n].Contains('\''))
                {
                    var firstIndex = str[n].IndexOf('\'');
                    var next = str[n].IndexOf('\'', firstIndex + 1);
                    if (next == -1)
                        break;
                    var substring = str[n].Substring(firstIndex, next - firstIndex + 1);
                    Replace(ref substring, " ", "薔");
                    Replace(ref substring, "'", "");

                    str[n] = $"{str[n][..firstIndex]}{substring}{str[n][(next + 1)..]}";
                }
            }

            if (str[n].Contains('"'))//\"が含まれていたら
            {
                while (str[n].Contains('"'))
                {
                    var firstIndex = str[n].IndexOf('"');
                    var next = str[n].IndexOf('"', firstIndex + 1);
                    if (next == -1)
                        break;
                    var substring = str[n].Substring(firstIndex, next - firstIndex + 1);
                    Replace(ref substring, " ", "薔");
                    Replace(ref substring, "\"", "");
                    str[n] = $"{str[n][..firstIndex]}{substring}{str[n][(next + 1)..]}";
                }
            }
        }

        //次にloop_に続く行が　"label  data"だった時に対応
        //for (int n = 0; n < str.Length - 1; n++) // 261111Cl 変更前 (9 行)
        //{
        //    if (str[n].StartsWith("loop_", Ord) && str[n + 1].Contains(' '))
        //    {
        //        var temp = str[n + 1].Split(' ', true);
        //        str[n + 1] = temp[0];
        //        Insert(ref str, n + 2, temp[1]);
        //    }
        //}
        // 261111Cl 変更: 挿入は List<string>.Insert で行う (旧は挿入のたびに配列全体を作り直していた)。挿入した行も続けて調べるのは同じ
        {
            var list = new List<string>(str);
            for (int n = 0; n < list.Count - 1; n++)
            {
                if (list[n].StartsWith("loop_", Ord) && list[n + 1].Contains(' '))
                {
                    var temp = list[n + 1].Split(' ', true);
                    list[n + 1] = temp[0];
                    list.Insert(n + 2, temp[1]);
                }
            }
            str = [.. list];
        }

        var CIF = new List<List<(string Label, string Data)>>();
        for (int n = 0; n < str.Length; n++)
        {
            if (str[n].Trim().StartsWith("_", Ord))
            {//単体アイテムのとき
                var temp = str[n].Split(' ', true);
                string tempLabel = temp[0], tempData;
                if (temp.Length == 2)
                    tempData = temp[1];
                else if (temp.Length == 1)
                {
                    n++;
                    tempData = str[n];
                }
                else
                    tempData = "";
                CIF.Add(new List<(string Label, string Data)>([(tempLabel, tempData.Replace("薔", " "))]));
            }
            else if (str[n].Trim().StartsWith("loop_", Ord))
            {//ループのとき
                var tempLoopLabels = new List<string>();
                var tempLoopDatas = new List<string>();
                n++;
                //"_"で始まるラベルを数える
                while (n < str.Length && str[n].Trim().StartsWith("_", Ord))
                    tempLoopLabels.Add(str[n++].Trim());

                //次に"_"か"loop_"か"#End of"で始まる行が出てくるまでループ
                while (n < str.Length && !str[n].Trim().StartsWith("_", Ord) && !str[n].Trim().StartsWith("loop_", Ord) && !str[n].Trim().StartsWith("#End of", Ord))
                {
                    var temp = str[n].Split(' ', true);
                    for (int i = 0; i < temp.Length; i++)
                        tempLoopDatas.Add(temp[i]);
                    n++;
                }
                n--;
                if (tempLoopDatas.Count % tempLoopLabels.Count == 0)
                {
                    for (int i = 0; i < tempLoopDatas.Count / tempLoopLabels.Count; i++)
                    {
                        var cif_Group = new List<(string Label, string Data)>();
                        for (int j = 0; j < tempLoopLabels.Count; j++)
                            cif_Group.Add((tempLoopLabels[j], tempLoopDatas[i * tempLoopLabels.Count + j].Replace("薔", " ")));
                        CIF.Add(cif_Group);
                    }
                }
            }
        }
        //ここまででCIF_Groupクラスのリストが完成

        //格子定数は、CIFファイル中に何回も記載されている場合があるため、リストにする。
        List<(int index, string value)> aList = [], bList = [], cList = [], alphaList = [], betaList = [], gammaList = [];

        string name = "", sectionTitle = "", journalNameFull = "", journalCodenASTM = "";
        string volume = "", year = "", pageFirst = "", pageLast = "", issue = "";
        var journal = new StringBuilder();
        List<string> spaceGroupNameHM = [], spaceGroupNameHall = [];
        string chemical_formula_sum = "", chemical_formula_structural = "";
        var symmetry_Int_Tables_number = -1;
        var author = new List<string>();
        var operations = new List<string>();

        for (int i = 0; i < CIF.Count; i++)
            for (int j = 0; j < CIF[i].Count; j++)
            {
                var label = CIF[i][j].Label;
                var data = CIF[i][j].Data;

                if (label.StartsWith("_chemical_name", Ord)) name += data + " ";

                //ここから格子定数
                if (label == "_cell_length_a") aList.Add((i, data));
                else if (label == "_cell_length_b") bList.Add((i, data));
                else if (label == "_cell_length_c") cList.Add((i, data));
                else if (label == "_cell_angle_alpha") alphaList.Add((i, data));
                else if (label == "_cell_angle_beta") betaList.Add((i, data));
                else if (label == "_cell_angle_gamma") gammaList.Add((i, data));

                //ここからジャーナル情報
                else if (label == "_publ_author_name") author.Add(data);
                else if (label == "_publ_section_title") sectionTitle = data;
                else if (label == "_journal_name_full") journalNameFull = data;
                else if (label == "_journal_coden_ASTM") journalCodenASTM = data;
                else if (label == "_journal_year") year = data;
                else if (label == "_journal_volume") volume = data;
                else if (label == "_journal_page_first") pageFirst = data;
                else if (label == "_journal_page_last") pageLast = data;
                else if (label == "_journal_issue") issue = data;
                //ここから対称性
                else if (label.Contains("_space_group_name_H-M")) spaceGroupNameHM.Add(data);
                else if (label.Contains("_space_group_name_Hall")) spaceGroupNameHall.Add(data);
                else if (label.Contains("_Int_Tables_number")) int.TryParse(data, out symmetry_Int_Tables_number);
                else if (label == "_chemical_formula_sum") chemical_formula_sum = data;
                else if (label == "_chemical_formula_structural") chemical_formula_structural = data;
                else if (label == "_space_group_symop_operation_xyz") operations.Add(data);
                else if (label == "_symmetry_equiv_pos_as_xyz") operations.Add(data);
            }

        if (aList.Count == 0 || bList.Count == 0 || cList.Count == 0 || alphaList.Count == 0 || betaList.Count == 0 || gammaList.Count == 0) return null;

        if (name.Length == 0 || name == "?" || name == "? ?" || name.Trim().Length == 0)
            name = chemical_formula_sum;

        //"_atom_site_label"というラベルを含むCIF番号をリストする。 (最初の連続した番号のみを使う)
        //さらに、格子定数Listの中で、どれが正しい値なのかを判別する
        string a = "", b = "", c = "", alpha = "", beta = "", gamma = "";
        var atomCIF = new List<List<(string Label, string Data)>>();
        bool flag = false;
        for (int i = 0; i < CIF.Count; i++)
            if (flag || atomCIF.Count == 0)
            {
                if (CIF[i].Exists(item => item.Label == "_atom_site_label"))
                {
                    atomCIF.Add(CIF[i]);
                    flag = true;
                }
                else if (atomCIF.Count != 0)
                {
                    flag = false;
                    a = aList.Count == 1 ? aList[0].value : aList.Last(e => e.index < i).value;
                    b = bList.Count == 1 ? bList[0].value : bList.Last(e => e.index < i).value;
                    c = cList.Count == 1 ? cList[0].value : cList.Last(e => e.index < i).value;
                    alpha = alphaList.Count == 1 ? alphaList[0].value : alphaList.Last(e => e.index < i).value;
                    beta = betaList.Count == 1 ? betaList[0].value : betaList.Last(e => e.index < i).value;
                    gamma = gammaList.Count == 1 ? gammaList[0].value : gammaList.Last(e => e.index < i).value;
                }
            }
        if (a == "") a = aList[0].value;
        if (b == "") b = bList[0].value;
        if (c == "") c = cList[0].value;
        if (alpha == "") alpha = alphaList[0].value;
        if (beta == "") beta = betaList[0].value;
        if (gamma == "") gamma = gammaList[0].value;

        #region 空間群を調べる部分
        //空間群を検索
        int sgnum = 0;
        if (spaceGroupNameHM.Count == 0 && spaceGroupNameHall.Count == 0)
            sgnum = 0;
        else
        {
            for (int i = 0; i < Math.Max(spaceGroupNameHall.Count, spaceGroupNameHM.Count); i++)
            {
                var HM = i < spaceGroupNameHM.Count ? spaceGroupNameHM[i] : "";
                var Hall = i < spaceGroupNameHall.Count ? spaceGroupNameHall[i] : "";
                //sgnum = SearchSGseriesNumberForCIF(HM, Hall, symmetry_Int_Tables_number,  // 261105Cl 変更前
                //    a == b && b == c && alpha == beta && beta == gamma, Is90(alpha), Is90(beta), Is90(gamma));
                sgnum = SearchSGseriesNumberForCIF(HM, Hall, symmetry_Int_Tables_number, // 261105Cl: 対称操作も渡す (曖昧な空白なし Hall 記号の判定用)
                    a == b && b == c && alpha == beta && beta == gamma, Is90(alpha), Is90(beta), Is90(gamma), operations); // 261105Cl
                if (sgnum != -1)
                    break;
            }
        }

        if (sgnum == -1)
            return null;

        #region 対象操作がCIFファイル中に記載されている場合は、本当に現在の空間群でよいかどうかをチェック

        var shift = new V3(0, 0, 0);
        var p = new V3(0.111, 0.234, 0.457);//適当な一般位置
        var tempAtom = WyckoffPosition.GetEquivalentAtomsPosition((p.X, p.Y, p.Z), sgnum).Atom;

        if (operations.Count != 0 && operations.Count == tempAtom.Length)
        {
            var th = 0.0000001;
            // 260427Cl DynamicExpresso へ置換。"x, y+1/2, z" のような対称操作 3 成分文字列を、カンマで分割して各成分を個別に Func<double,double,double,double> へコンパイルする。
            // var prms = new[] { "x", "y", "z" }.Select(s => Expression.Parameter(typeof(double), s)).ToArray();
            var interpreter = new Interpreter();
            //文字列からラムダ式を返すローカル関数
            Func<double, double, double, V3> func(string sExpr)
            {
                try
                {
                    sExpr = sExpr.Replace(" ", "").Replace(",+", ",").TrimStart(['+']);
                    sExpr = sExpr.Replace("/", ".0/").Replace(".0.0", ".0");//分子に小数点を加える
                    // (260427Ch) DynamicExpresso は "2x" のような暗黙の掛け算を受け付けないため、CIF/特殊位置表記を明示的な掛け算に直す。
                    // foreach (var v in new[] { 'x', 'y', 'z' })                                  // 260427Cl 旧: 30回の Replace ループ → 単一 Regex に統合
                    //     for (var coefficient = 0; coefficient <= 9; coefficient++)
                    //         sExpr = sExpr.Replace($"{coefficient}{v}", $"{coefficient}*{v}");
                    sExpr = ImplicitMulRegex.Replace(sExpr, "$1*$2"); // 260427Cl
                    // sExpr = "new [] {" + sExpr.Replace("/", ".0/").Replace(".0.0", ".0") + "}"; // 260427Cl 旧: System.Linq.Dynamic.Core 用の配列リテラル形式
                    // var f = DynamicExpressionParser.ParseLambda(prms, typeof(double[]), sExpr).Compile() as Func<double, double, double, double[]>;
                    // return (x, y, z) => { var d = f(x, y, z); return new V3(d[0], d[1], d[2]); };
                    var parts = sExpr.Split(',');
                    if (parts.Length != 3) return null;
                    var fx = interpreter.ParseAsDelegate<Func<double, double, double, double>>(parts[0], "x", "y", "z");
                    var fy = interpreter.ParseAsDelegate<Func<double, double, double, double>>(parts[1], "x", "y", "z");
                    var fz = interpreter.ParseAsDelegate<Func<double, double, double, double>>(parts[2], "x", "y", "z");
                    return (x, y, z) => new V3(fx(x, y, z), fy(x, y, z), fz(x, y, z));
                }
                catch (Exception e)
                {
                    if (AssemblyState.IsDebug)
                        System.Windows.Forms.MessageBox.Show(e.Message);
                    return null;
                }
            }

            var funcs = operations.Select(s => func(s)).ToArray();
            if (funcs.All(f => f != null))
            {
                var shiftCandidates = new[] { 0, 0.125, 0.25, -0.125, -0.25 };

                var temp = tempAtom.Select(a => norm(new V3(a.X, a.Y, a.Z))).ToList();
                temp.Sort((o1, o2) => Math.Abs(o1.X - o2.X) > th ? o1.X.CompareTo(o2.X) : (Math.Abs(o1.Y - o2.Y) > th ? o1.Y.CompareTo(o2.Y) : o1.Z.CompareTo(o2.Z)));
                var source = temp.Select((pos, i) => (pos, i)).ToList();

                var match = false;
                foreach (var sX in shiftCandidates)
                    foreach (var sY in shiftCandidates)
                        foreach (var sZ in shiftCandidates)
                            if (!match)
                            {
                                // var target = funcs.Select(f => f(p.X + sX, p.Y + sY, p.Z + sZ)).Select(r => norm(new V3(r.X - sX, r.Y - sY, r.Z - sY))).ToList(); // 260427Ch 旧: Z 成分から sY を引いていた
                                var target = funcs.Select(f => f(p.X + sX, p.Y + sY, p.Z + sZ)).Select(r => norm(new V3(r.X - sX, r.Y - sY, r.Z - sZ))).ToList(); // (260427Ch)
                                target.Sort((o1, o2) => Math.Abs(o1.X - o2.X) > th ? o1.X.CompareTo(o2.X) : (Math.Abs(o1.Y - o2.Y) > th ? o1.Y.CompareTo(o2.Y) : o1.Z.CompareTo(o2.Z)));
                                if (source.All(o => (o.pos - target[o.i]).Length < th))
                                {
                                    match = true;
                                    shift = new V3(sX, sY, sZ);
                                }

                            }
            }
        }
        #endregion

        #endregion

        var atoms = new List<Atoms2>();
        foreach (var cif in atomCIF)
        {
            //原子の情報
            string atomLabel = "", atomSymbol = "", thermalDisplaceType = "";
            string x = "", y = "", z = "", occ = "1";
            string uIso = "", u11 = "", u22 = "", u33 = "", u12 = "", u13 = "", u23 = "";
            string bIso = "", b11 = "", b22 = "", b33 = "", b12 = "", b13 = "", b23 = "";

            //まず基本的な原子位置や占有率などの情報を探す
            for (int j = 0; j < cif.Count; j++)
            {
                var label = cif[j].Label;
                var data = cif[j].Data;
                if (label == "_atom_site_type_symbol") atomSymbol = data;
                else if (label == "_atom_site_label") atomLabel = data;
                else if (label == "_atom_site_thermal_displace_type") thermalDisplaceType = data;
                else if (data != "?" && data != ".")
                {
                    if (label == "_atom_site_fract_x") x = data;// + shift.X;
                    else if (label == "_atom_site_fract_y") y = data;// + shift.X;
                    else if (label == "_atom_site_fract_z") z = data;// + shift.X;
                    else if (label == "_atom_site_occupancy") occ = data;
                    else if (label == "_atom_site_U_iso_or_equiv") uIso = data;
                    else if (label == "_atom_site_B_iso_or_equiv") bIso = data;
                }
            }

            if (shift.X != 0 || shift.Y != 0 || shift.Z != 0)
            {
                //var _x = Crystal2.Decompose(x, sgnum); // 261105Cl 変更前 (6 行)
                //var _y = Crystal2.Decompose(y, sgnum);
                //var _z = Crystal2.Decompose(z, sgnum);
                //x = Crystal2.Compose(_x.Value + shift.X, _x.Error);
                //y = Crystal2.Compose(_y.Value + shift.Y, _y.Error);
                //z = Crystal2.Compose(_z.Value + shift.Z, _z.Error);
                // 261105Cl 変更: x などは CIF の生の値 ("0.4689(2)")。Crystal2.Decompose は内部形式 (|誤差|) 用で括弧を読めず NaN になり、
                // Compose が空文字を返すので、原点をずらした結晶の誤差つきの座標が全部消えていた (COD 2104403 など。tools/AdpCheck r3)。
                // CIF の数値として読んでずらす。CIF の数値でないもの (分数 "1/4" など) は従来どおり Decompose で読む
                string shiftText(string s, double d) => scaleValueText(s, 1, d) is { Length: > 0 } t ? t
                    : Crystal2.Decompose(s, sgnum) is var (v, e) && !double.IsNaN(v) ? Crystal2.Compose(v + d, e) : s;
                (x, y, z) = (shiftText(x, shift.X), shiftText(y, shift.Y), shiftText(z, shift.Z));
            }

            //次に異方性の温度散乱因子をさがす (等方性の温度因子は既に上のループで読み込まれている)
            for (int k = 0; k < CIF.Count; k++)
            {
                if (CIF[k].Exists(item => item.Label == "_atom_site_aniso_label" && item.Data == atomLabel))
                {
                    for (int l = 0; l < CIF[k].Count; l++)
                    {
                        var label = CIF[k][l].Label;
                        var data = CIF[k][l].Data;
                        if (data != "?" && data != ".")
                        {
                            if (label == "_atom_site_aniso_U_11") u11 = data;
                            else if (label == "_atom_site_aniso_U_22") u22 = data;
                            else if (label == "_atom_site_aniso_U_33") u33 = data;
                            else if (label == "_atom_site_aniso_U_12") u12 = data;
                            else if (label == "_atom_site_aniso_U_13") u13 = data;
                            else if (label == "_atom_site_aniso_U_23") u23 = data;
                            else if (label == "_atom_site_aniso_B_11") b11 = data;
                            else if (label == "_atom_site_aniso_B_22") b22 = data;
                            else if (label == "_atom_site_aniso_B_33") b33 = data;
                            else if (label == "_atom_site_aniso_B_12") b12 = data;
                            else if (label == "_atom_site_aniso_B_13") b13 = data;
                            else if (label == "_atom_site_aniso_B_23") b23 = data;
                        }
                    }
                    break;
                }
            }
            //ラベル名から元素を探す
            int atomicNumber = 0;
            var atomName = atomSymbol.Length == 0 ? atomLabel : atomSymbol;

            for (int q = atomName.Length; q > 0 && atomicNumber == 0; q--)
            {
                //var temp = atomName[..q]; // 261111Cl 変更前 (6 行)
                //for (int k = 0; k <= 96 && atomicNumber == 0; k++)
                //{
                //    if (AtomStatic.AtomicName(k).ToLower() == temp.ToLower())
                //        atomicNumber = k;
                //}
                // 261111Cl 変更: 1 文字ずつ ToLower() した文字列を q×97 回作っていたのを、span 版 (ASCII の大文字小文字を区別しない Ordinal 比較) にする。
                // 現在のカルチャに依らなくなる (tr-TR で "SI1" が S になっていた)。日本語・英語の環境では結果は同じ
                var temp = atomName.AsSpan(0, q);
                atomicNumber = AtomStatic.AtomicNumber(temp, caseSensitive: false);

                if (temp is "OH")
                {
                    atomicNumber = -1;
                    break;
                }
                else if (temp is "D")
                {
                    atomicNumber = 255;
                    break;
                }
            }

            ////Bタイプが全て ”” だったら、Uタイプと判定 // 261105Cl 変更前 (14 行)
            //var isU = bIso.Length == 0 && b11.Length == 0 && b12.Length == 0 && b13.Length == 0 && b22.Length == 0 && b23.Length == 0 && b33.Length == 0;
            ////非等方性が全て ”” だったら、等方性と判断
            //var isIso = isU ?
            //    u11.Length == 0 && u12.Length == 0 && u13.Length == 0 && u22.Length == 0 && u23.Length == 0 && u33.Length == 0 :
            //    b11.Length == 0 && b12.Length == 0 && b13.Length == 0 && b22.Length == 0 && b23.Length == 0 && b33.Length == 0;
            //
            //var iso = isU ? uIso : bIso;
            //
            //if (iso.Length == 0)
            //    iso = "0";
            //
            //string[] aniso = isU ? //11, 22, 33, 12, 23, 31の順番
            //    [u11, u22, u33, u12, u23, u13] :
            //    [b11, b22, b33, b12, b23, b13];
            // 261105Cl 変更: 異方性があるときは U 型にそろえる。
            // (1) CIF の _atom_site_aniso_B_ij は Å² の B (= 8π²U、CIF core 辞書) だが、ReciPro の B 型の異方性は無次元の β (exp(−hᵀβh)) なので、
            //     以前のように B_ij をそのまま B 型として入れると DW が桁違いになっていた (格子の大きさで変わる。tools/AdpCheck cifb で 260 倍)。B_ij / 8π² = U_ij に直す。
            // (2) 以前は _atom_site_B_iso_or_equiv があるだけで B 型と判定したため、B_iso と _atom_site_aniso_U_ij が並ぶ CIF では U_ij を捨てていた
            bool anyAnisoU = u11.Length + u22.Length + u33.Length + u12.Length + u13.Length + u23.Length > 0;
            bool anyAnisoB = b11.Length + b22.Length + b33.Length + b12.Length + b13.Length + b23.Length > 0;
            const double bToU = 1 / (8 * Math.PI * Math.PI);
            bool isU, isIso;
            string iso;
            string[] aniso; //11, 22, 33, 12, 23, 31の順番
            if (anyAnisoU || anyAnisoB)
            {
                isU = true;
                isIso = false;
                aniso = anyAnisoU ? [u11, u22, u33, u12, u23, u13] : [.. new[] { b11, b22, b33, b12, b23, b13 }.Select(s => scaleValueText(s, bToU))];
                iso = uIso.Length > 0 ? uIso : bIso.Length > 0 ? scaleValueText(bIso, bToU) : "0";
            }
            else
            {
                isU = bIso.Length == 0;
                isIso = true;
                iso = isU ? uIso : bIso;
                aniso = isU ? [u11, u22, u33, u12, u23, u13] : [b11, b22, b33, b12, b23, b13];
            }
            if (iso.Length == 0)
                iso = "0";

            if (atomicNumber > 0)
                atoms.Add(new Atoms2(atomLabel, atomicNumber, 0, 0, [x, y, z], occ, isIso, isU, iso, aniso));
            else if (atomicNumber == -1)//"OH"のときの対処
            {
                atoms.Add(new Atoms2(atomLabel, 1, 0, 0, [x, y, z], occ, isIso, isU, iso, aniso));
                atoms.Add(new Atoms2(atomLabel, 8, 0, 0, [x, y, z], occ, isIso, isU, iso, aniso));
            }
        }

        if (journalNameFull != "" || journalCodenASTM != "") journal.Append(journalNameFull).Append(' ').Append(journalCodenASTM);
        if (issue != "") journal.Append(", ").Append(issue);
        if (volume != "") journal.Append(", ").Append(volume);
        if (year != "") journal.Append('(').Append(year).Append(')');
        if (pageFirst != "") journal.Append(", ").Append(pageFirst);
        if (pageLast != "") journal.Append('-').Append(pageLast);

        var authors = new StringBuilder();
        for(int i = 0; i< author.Count; i++)
            if(i== author.Count -1)
                authors.Append(author[i].Replace(",", ""));
            else
                authors.Append(author[i].Replace(",", "") + ", \r\n");

        return new Crystal2
        {
            CellTexts = [a, b, c, alpha, beta, gamma],
            sym = (short)sgnum,
            name = name,
            atoms = atoms,
            auth = authors.ToString().TrimEnd([' ', ',']).Replace(".", ""),
            jour = journal.ToString(),
            sect = sectionTitle,
            argb = Color.FromArgb(rnd.Next(255), rnd.Next(255), rnd.Next(255)).ToArgb()
        };
    }

 

    private static V3 norm(V3 v)
    {
        var d = new[] { v.X, v.Y, v.Z };
        for (int i = 0; i < 3; i++)
        {
            while (d[i] > 0.9999999)
                d[i]--;
            while (d[i] < -0.0000001)
                d[i]++;
        }
        return new V3(d[0], d[1], d[2]);
    }

    /// <summary>261111Cl 追加: 530 設定の空間群記号を、SearchSGseriesNumberForCIF で比べる形に直した表。CIF 1 件ごとに 530 設定ぶんの
    /// Split/Join と Replace (計 530×4 回) をしていたのを、最初に使うときの 1 回だけにする。値は旧来の式と同じ</summary>
    private static class SgNameTables
    {
        /// <summary>Hall 記号の連続する空白・タブを空白 1 つにまとめる (旧 SearchSGseriesNumberForCIF 内の normalizeHall と同じ式)</summary>
        public static string NormalizeHall(string s) => string.Join(" ", s.Split([' ', '\t'], StringSplitOptions.RemoveEmptyEntries));
        /// <summary>標準 Hall 記号 (SpaceGroupHallStandardStr) を NormalizeHall したもの</summary>
        public static readonly string[] HallStandardNormalized;
        /// <summary>標準 Hall 記号から空白を除いたもの</summary>
        public static readonly string[] HallStandardCompact;
        /// <summary>従来の Hall 記号 (SpaceGroupHallStr) から空白を除いたもの</summary>
        public static readonly string[] HallCompact;
        /// <summary>完全 H-M 記号 (SpaceGroupHMfullStr) から空白を除いたもの</summary>
        public static readonly string[] HMfullCompact;
        static SgNameTables()
        {
            var syms = SymmetryStatic.Symmetries;
            int n = SymmetryStatic.TotalSpaceGroupNumber;
            HallStandardNormalized = new string[n];
            HallStandardCompact = new string[n];
            HallCompact = new string[n];
            HMfullCompact = new string[n];
            for (int i = 0; i < n; i++)
            {
                HallStandardNormalized[i] = NormalizeHall(syms[i].SpaceGroupHallStandardStr);
                HallStandardCompact[i] = syms[i].SpaceGroupHallStandardStr.Replace(" ", "");
                HallCompact[i] = syms[i].SpaceGroupHallStr.Replace(" ", "");
                HMfullCompact[i] = syms[i].SpaceGroupHMfullStr.Replace(" ", "");
            }
        }
    }

    //private static int SearchSGseriesNumberForCIF(string SgNameHM, string SgNameHall, int SpaceGroupNumber, bool isRhomboShape, bool isAlpha90, bool isBeta90, bool isGamma90) // 261105Cl 変更前
    /// <param name="operations">261105Cl 追加: CIF の対称操作 (_space_group_symop_operation_xyz / _symmetry_equiv_pos_as_xyz)。
    /// 空白なしの Hall 記号が違う群の複数の設定に当たるとき、どれを採るかの判定に使う。null や空なら使わない。</param>
    private static int SearchSGseriesNumberForCIF(string SgNameHM, string SgNameHall, int SpaceGroupNumber, bool isRhomboShape, bool isAlpha90, bool isBeta90, bool isGamma90, IEnumerable<string> operations = null) // 261105Cl operations を追加
    {
        int symmetrySeriesNumber = -1;
        if (SgNameHall != "")
        {
            // 261105Cl 修正: Hall 記号では空白が区切りの意味を持つ (例 'P 3 2' = P312 と 'P 32' = P3_2)。旧実装は CIF 側の空白を除き、
            //   空白なしの表 (SpaceGroupHallStr) と比べていたため、P312 → P3_2、P622 → P6_2 に誤って割り当てていた。また表の表記が標準
            //   (ITA/spglib) と異なる 57 設定では一致せず H-M 経路に落ちていた。標準表記 (SpaceGroupHallStandardStr, 空白あり) と、
            //   連続する空白を 1 つにまとめた文字列どうしで比べる。空白なしの表記 (旧版が書き出した CIF など) は、空白を除いた比較で
            //   1 つの設定に決まる場合だけ受け付ける。
            //var hall = SgNameHall.Replace(" ", ""); // 261105Cl 変更前
            //for (int i = 0; i < SymmetryStatic.TotalSpaceGroupNumber; i++)
            //    if (hall == SymmetryStatic.Symmetries[i].SpaceGroupHallStr)
            //    {
            //        symmetrySeriesNumber = i;
            //        break;
            //    }
            //if (symmetrySeriesNumber != -1)
            //    return symmetrySeriesNumber;
            //static string normalizeHall(string s) => string.Join(" ", s.Split([' ', '\t'], StringSplitOptions.RemoveEmptyEntries)); // 261111Cl 変更前 (SgNameTables.NormalizeHall へ)
            //var hall = normalizeHall(SgNameHall);
            //for (int i = 0; i < SymmetryStatic.TotalSpaceGroupNumber; i++)
            //    if (hall == normalizeHall(SymmetryStatic.Symmetries[i].SpaceGroupHallStandardStr))
            //        return i;
            // 261111Cl 変更: 530 設定側の正規化は SgNameTables に 1 回だけ作った表を使う (最初に一致した設定を採るのは同じ)
            var hall = SgNameTables.NormalizeHall(SgNameHall);
            for (int i = 0; i < SymmetryStatic.TotalSpaceGroupNumber; i++)
                if (hall == SgNameTables.HallStandardNormalized[i])
                    return i;

            // 空白なしの比較は、標準表記と従来表記 (SpaceGroupHallStr。旧版が CIF に書き出していた表記) の両方で行う。
            // 'P32' (P3_2 と P312) や 'P62' (P6_2 と P622) のように複数の設定に当たるときは採らず、H-M 記号で決める。
            // 261105Cl 追記: 複数の設定に当たっても、それらがすべて同じ群 (対称操作の集合が同じ。'C22-1bc' の Ccce(1) abc と ba-c など) なら最初の設定を採る。
            //   違う群のときは、CIF の対称操作と一致する設定がちょうど 1 つならそれを採る。決まらなければ H-M 記号で、H-M 記号も無ければ番号で決める。
            var compact = hall.Replace(" ", "");
            // 261105Cl 追加: 群どうしの比較は、一般位置 p の像の集合で行う (座標の差が格子併進なら同じ点)。
            //   CIF の対称操作は、x, y, z と整数・分数だけの式 ('-x+y, -x, z+2/3' など) なら読む。それ以外 (小数など) を含むときは使わない。
            var p = (X: 0.111, Y: 0.234, Z: 0.457); // 261105Cl 追加 (下の対称操作の確認と同じ一般位置)
            static bool covers(IEnumerable<(double X, double Y, double Z)> a, IEnumerable<(double X, double Y, double Z)> b) // 261105Cl 追加: a の各点が b のどれかと格子併進を除いて一致
                => a.All(u => b.Any(v => Math.Abs(Math.IEEERemainder(u.X - v.X, 1)) < 1e-6 && Math.Abs(Math.IEEERemainder(u.Y - v.Y, 1)) < 1e-6 && Math.Abs(Math.IEEERemainder(u.Z - v.Z, 1)) < 1e-6));
            (double X, double Y, double Z)[] cifOrbit = null; // 261105Cl 追加
            try // 261105Cl 追加
            { // 261105Cl 追加
                var ops = operations?.Select(o => o.ToLowerInvariant().Replace(" ", "").Split(',')).ToArray(); // 261105Cl 追加
                if (ops?.Length > 0 && ops.All(c => c.Length == 3 && c.All(t => t.Length > 0 && t.All("xyz0123456789+-/".Contains)))) // 261105Cl 追加 (AffineGen.Parse は想定外の文字 (小数点など) があると先へ進まず無限ループになるので、先に弾く)
                    cifOrbit = [.. ops.Select(c => AffineGen.Parse(c[0], c[1], c[2]).Apply(p.X, p.Y, p.Z))]; // 261105Cl 追加
            } // 261105Cl 追加
            catch (Exception e) when (e is FormatException or OverflowException) { cifOrbit = null; } // 261105Cl 追加 ('1/' など)
            // 261105Cl 追記: その後 AffineGen.Parse 自身も読めない式で FormatException を投げるようにした (無限ループしない)。上の文字の確認は、大文字を小文字にしてから読むために残す
            foreach (var standard in new[] { true, false })
            {
                //int hit = -1, hitCount = 0; // 261105Cl 変更前
                var hits = new List<int>(); // 261105Cl
                var table = standard ? SgNameTables.HallStandardCompact : SgNameTables.HallCompact; // 261111Cl 追加
                for (int i = 0; i < SymmetryStatic.TotalSpaceGroupNumber; i++)
                {
                    //var s = SymmetryStatic.Symmetries[i]; // 261111Cl 変更前 (2 行)
                    //if (compact == (standard ? s.SpaceGroupHallStandardStr : s.SpaceGroupHallStr).Replace(" ", ""))
                    if (compact == table[i]) // 261111Cl 変更: 空白を除いた表を 1 回だけ作って比べる (一致の順・件数は同じ)
                        hits.Add(i); // 261105Cl
                    //{ // 261105Cl 変更前
                    //    hit = i;
                    //    hitCount++;
                    //}
                }
                //if (hitCount == 1) // 261105Cl 変更前
                //    return hit;
                if (hits.Count == 1) // 261105Cl
                    return hits[0]; // 261105Cl
                if (hits.Count > 1) // 261105Cl 追加
                {
                    var orbits = hits.Select(i => SymmetryStatic.WyckoffPositions[i][0].GeneratePositions(p.X, p.Y, p.Z).Select(v => (v.X, v.Y, v.Z)).ToArray()).ToArray(); // 261105Cl 追加
                    if (orbits.All(o => covers(o, orbits[0]) && covers(orbits[0], o))) // 261105Cl 追加: すべて同じ群
                        return hits[0]; // 261105Cl 追加
                    var matched = hits.Where((_, k) => cifOrbit != null && covers(orbits[k], cifOrbit) && covers(cifOrbit, orbits[k])).ToArray(); // 261105Cl 追加
                    if (matched.Length == 1) // 261105Cl 追加
                        return matched[0]; // 261105Cl 追加
                }
            }
        }

        SgNameHM = SgNameHM.TrimStart(' ').TrimEnd(' ');


        if (!SgNameHM.Contains(' ') && SgNameHM.Contains('_'))
            for (int i = 1; i < SgNameHM.Length - 1; i++)
                if (SgNameHM[i] == '_' && '0' < SgNameHM[i - 1] && '9' > SgNameHM[i - 1] && '0' < SgNameHM[i + 1] && '9' > SgNameHM[i + 1])
                {
                    SgNameHM = SgNameHM.Remove(i, 1);
                    SgNameHM = SgNameHM.Insert(i, "sub");
                    i += 3;
                }

        Replace(ref SgNameHM, "_", " ");

        Replace(ref SgNameHM, "{hexagonal axes}", " ");
        Replace(ref SgNameHM, "{rhombohedral axes}", " ");

        SgNameHM = SgNameHM.TrimStart(' ').TrimEnd(' ');

        if (SgNameHM.EndsWith("RS", Ord) || SgNameHM.EndsWith("HR", Ord))
            SgNameHM = SgNameHM[..^2].TrimEnd(' ');

        if (SgNameHM.EndsWith("H", Ord) || SgNameHM.EndsWith("h", Ord) || SgNameHM.EndsWith("R", Ord) || SgNameHM.EndsWith("r", Ord))
            SgNameHM = SgNameHM[..^1].TrimEnd(' ');

        if (SgNameHM.EndsWith(":", Ord))
            SgNameHM = SgNameHM[..^1].TrimEnd(' ');

        bool IsOrigineChoice2 = false;
        if (SgNameHM.EndsWith("Z", Ord))//最後にZがついていたらOriginChoice2
        {
            IsOrigineChoice2 = true;
            SgNameHM = SgNameHM.TrimEnd('Z').TrimEnd();
        }
        if (SgNameHM.EndsWith("S", Ord))//最後がSだったらOriginChoice1
            SgNameHM = SgNameHM.TrimEnd('S').TrimEnd();

        if (SgNameHM.EndsWith("S1", Ord))//最後がS1だったらOriginChoice1
            SgNameHM = SgNameHM.Replace("S1", "").TrimEnd();
        if (SgNameHM.EndsWith("Z1", Ord))//最後がZ1だったらOriginChoice1
            SgNameHM = SgNameHM.Replace("Z1", "").TrimEnd();
        if (SgNameHM.EndsWith(":S1", Ord))//最後が:S1だったらOriginChoice1
            SgNameHM = SgNameHM.Replace(":S1", "").TrimEnd();
        if (SgNameHM.EndsWith(":1", Ord))//最後が:1だったらOriginChoice1
            SgNameHM = SgNameHM.Replace(":1", "").TrimEnd();

        if (SgNameHM.EndsWith("S2", Ord))//最後がS2だったらOriginChoice2
        {
            SgNameHM = SgNameHM.Replace("S2", "").TrimEnd();
            IsOrigineChoice2 = true;
        }
        if (SgNameHM.EndsWith("Z2", Ord))//最後がZ2だったらOriginChoice2
        {
            SgNameHM = SgNameHM.Replace("Z2", "").TrimEnd();
            IsOrigineChoice2 = true;
        }
        if (SgNameHM.EndsWith(":S2", Ord))//最後が:S2だったらOriginChoice2
        {
            SgNameHM = SgNameHM.Replace(":S2", "").TrimEnd();
            IsOrigineChoice2 = true;
        }
        if (SgNameHM.EndsWith(":2", Ord))//最後が:S2だったらOriginChoice2
        {
            SgNameHM = SgNameHM.Replace(":2", "").TrimEnd();
            IsOrigineChoice2 = true;
        }
        if (SgNameHM.EndsWith("O2", Ord))//最後が:S2だったらOriginChoice2
        {
            SgNameHM = SgNameHM.Replace("O2", "").TrimEnd();
            IsOrigineChoice2 = true;
        }

     

        Replace(ref SgNameHM, "~", "");

        // 261105Cl 追加: H-M 記号が無い (空白なしの Hall 記号が違う群の複数の設定に当たり、対称操作でも決まらなかった場合など) と、
        //   下の SgNameHM[0] で IndexOutOfRangeException になり、呼び出し側の catch で CIF が開けなかった。H-M 記号では探せないので、
        //   末尾の救済と同じく番号の最初の設定を採る。番号も無ければ -1 (決まらない)。
        if (SgNameHM.Length == 0) // 261105Cl 追加
            return SpaceGroupNumber is >= 1 and <= 230 ? Array.FindIndex(SymmetryStatic.Symmetries, s => s.SpaceGroupNumber == SpaceGroupNumber) : -1; // 261105Cl 追加

        //一文字目以降の英字は全て小文字に
        SgNameHM = SgNameHM[0] + SgNameHM[1..].ToLower();

        Replace(ref SgNameHM,"P(-1)", "P-1");

        Replace(ref SgNameHM,"R 32", "R32");

        Replace(ref SgNameHM,"P3(1)21", "P 3sub1 2 1");
        Replace(ref SgNameHM,"P3(2)21", "P 3sub2 2 1");
        Replace(ref SgNameHM,"P3121", "P 3sub1 2 1");
        Replace(ref SgNameHM,"P3221", "P 3sub2 2 1");
        Replace(ref SgNameHM,"P 31 21", "P 3sub1 2 1");
        Replace(ref SgNameHM,"P 32 21", "P 3sub2 2 1");
        Replace(ref SgNameHM,"P321", "P 3 2 1");

        Replace(ref SgNameHM,"P61", "P 6sub1");
        Replace(ref SgNameHM,"P62", "P 6sub2");
        Replace(ref SgNameHM,"P63", "P 6sub3");
        Replace(ref SgNameHM,"P64", "P 6sub4");
        Replace(ref SgNameHM,"P65", "P 6sub5");

        Replace(ref SgNameHM,"I41", "I 4sub1");
        Replace(ref SgNameHM,"P42", "P 4sub2");

        Replace(ref SgNameHM,"P42/m b c", "P 4sub2/m b c");

        Replace(ref SgNameHM,"P4322", "P 4sub3 2 2");

        Replace(ref SgNameHM," 61", " 6sub1");
        Replace(ref SgNameHM," 62", " 6sub2");
        Replace(ref SgNameHM," 63", " 6sub3");
        Replace(ref SgNameHM," 64", " 6sub4");
        Replace(ref SgNameHM," 65", " 6sub5");
        Replace(ref SgNameHM," 41", " 4sub1");
        Replace(ref SgNameHM," 42", " 4sub2");
        Replace(ref SgNameHM," 43", " 4sub3");
        Replace(ref SgNameHM," 31", " 3sub1");
        Replace(ref SgNameHM," 32", " 3sub2");
        Replace(ref SgNameHM," 21", " 2sub1");

        Replace(ref SgNameHM,"2(1)", " 2sub1");
        Replace(ref SgNameHM,"3(1)", " 3sub1");
        Replace(ref SgNameHM,"3(2)", " 3sub2");
        Replace(ref SgNameHM,"4(1)", " 4sub1");
        Replace(ref SgNameHM,"4(2)", " 4sub2");
        Replace(ref SgNameHM,"4(3)", " 4sub3");
        Replace(ref SgNameHM,"6(1)", " 6sub1");
        Replace(ref SgNameHM,"6(2)", " 6sub2");
        Replace(ref SgNameHM,"6(3)", " 6sub3");
        Replace(ref SgNameHM,"6(4)", " 6sub4");
        Replace(ref SgNameHM,"6(5)", " 6sub5");

        Replace(ref SgNameHM,"21", " 2sub1");

        Replace(ref SgNameHM,"  ", " ");

        string temp = SgNameHM.Replace(" ", "");
        #region
        if (temp == "Pm3") SgNameHM = "P m -3";
        else if (temp == "Pn3") SgNameHM = "P n -3";
        else if (temp == "Fm3") SgNameHM = "F m -3";
        else if (temp == "Fd3") SgNameHM = "F d -3";
        else if (temp == "Im3") SgNameHM = "I m -3";
        else if (temp == "Pa3") SgNameHM = "P a -3";
        else if (temp == "Ia3") SgNameHM = "I a -3";

        else if (temp == "Pm3m") SgNameHM = "P m -3 m";
        //else if (temp == "Pn3n") SgNameHM = "P m -3 m";
        else if (temp == "Pn3n") SgNameHM = "P n -3 n"; // 260718Cl: AMC 経路 (SearchSGseriesNumberForAmc) と同じ Pn3n→Pm-3m 誤マップを修正
        else if (temp == "Pm3n") SgNameHM = "P m -3 n";
        else if (temp == "Pn3m") SgNameHM = "P n -3 m";
        else if (temp == "Fm3m") SgNameHM = "F m -3 m";
        else if (temp == "Fm3c") SgNameHM = "F m -3 c";
        else if (temp == "Fd3m") SgNameHM = "F d -3 m";
        else if (temp == "Fd3c") SgNameHM = "F d -3 c";
        else if (temp == "Im3m") SgNameHM = "I m -3 m";
        else if (temp == "Ia3d") SgNameHM = "I a -3 d";
        //else if (temp == "I2sub1/a-3") SgNameHM = "I a 3";
        else if (temp == "I2sub1/a-3") SgNameHM = "I a -3"; // 260718Cl: マイナス欠落。空白除去後の検索語 "Ia3" はテーブル名 "Ia-3" と絶対に一致しない (AMC 経路は正しく Ia-3)


        else if (temp == "R-32/c") SgNameHM = "R -3 c";

        else if (temp == "P2" && isAlpha90 && isGamma90) SgNameHM = "P 1 2 1";
        else if (temp == "P2" && isAlpha90 && isBeta90) SgNameHM = "P 1 1 2";
        else if (temp == "P2" && isBeta90 && isGamma90) SgNameHM = "P 2 1 1";
        else if (temp == "P2sub1" && isAlpha90 && isGamma90) SgNameHM = "P 1 2sub1 1";
        else if (temp == "P2sub1" && isAlpha90 && isBeta90) SgNameHM = "P 1 1 2sub1";
        else if (temp == "P2sub1" && isBeta90 && isGamma90) SgNameHM = "P 2sub1 1 1";
        else if (temp == "C2" && isAlpha90 && isGamma90) SgNameHM = "C 1 2 1";
        else if (temp == "A2" && isAlpha90 && isGamma90) SgNameHM = "A 1 2 1";
        else if (temp == "I2" && isAlpha90 && isGamma90) SgNameHM = "I 1 2 1";
        else if (temp == "A2" && isAlpha90 && isBeta90) SgNameHM = "A 1 1 2";
        else if (temp == "B2" && isAlpha90 && isBeta90) SgNameHM = "B 1 1 2";
        else if (temp == "I2" && isAlpha90 && isBeta90) SgNameHM = "I 1 1 2";
        else if (temp == "B2" && isBeta90 && isGamma90) SgNameHM = "B 2 1 1";
        else if (temp == "C2" && isBeta90 && isGamma90) SgNameHM = "C 2 1 1";
        else if (temp == "I2" && isBeta90 && isGamma90) SgNameHM = "I 2 1 1";
        else if (temp == "Pm" && isAlpha90 && isGamma90) SgNameHM = "P 1 m 1";
        else if (temp == "Pm" && isAlpha90 && isBeta90) SgNameHM = "P 1 1 m";
        else if (temp == "Pm" && isBeta90 && isGamma90) SgNameHM = "P m 1 1";
        else if (temp == "Pc" && isAlpha90 && isGamma90) SgNameHM = "P 1 c 1";
        else if (temp == "Pn" && isAlpha90 && isGamma90) SgNameHM = "P 1 n 1";
        else if (temp == "Pa" && isAlpha90 && isGamma90) SgNameHM = "P 1 a 1";
        else if (temp == "Pa" && isAlpha90 && isBeta90) SgNameHM = "P 1 1 a";
        else if (temp == "Pn" && isAlpha90 && isBeta90) SgNameHM = "P 1 1 n";
        else if (temp == "Pb" && isAlpha90 && isBeta90) SgNameHM = "P 1 1 b";
        else if (temp == "Pb" && isBeta90 && isGamma90) SgNameHM = "P b 1 1";
        else if (temp == "Pn" && isBeta90 && isGamma90) SgNameHM = "P n 1 1";
        else if (temp == "Pc" && isBeta90 && isGamma90) SgNameHM = "P c 1 1";
        else if (temp == "Cm" && isAlpha90 && isGamma90) SgNameHM = "C 1 m 1";
        else if (temp == "Am" && isAlpha90 && isGamma90) SgNameHM = "A 1 m 1";
        else if (temp == "Im" && isAlpha90 && isGamma90) SgNameHM = "I 1 m 1";
        else if (temp == "Am" && isAlpha90 && isBeta90) SgNameHM = "A 1 1 m";
        else if (temp == "Bm" && isAlpha90 && isBeta90) SgNameHM = "B 1 1 m";
        else if (temp == "Im" && isAlpha90 && isBeta90) SgNameHM = "I 1 1 m";
        else if (temp == "Bm" && isBeta90 && isGamma90) SgNameHM = "B m 1 1";
        else if (temp == "Cm" && isBeta90 && isGamma90) SgNameHM = "C m 1 1";
        else if (temp == "Im" && isBeta90 && isGamma90) SgNameHM = "I m 1 1";
        else if (temp == "Cc" && isAlpha90 && isGamma90) SgNameHM = "C 1 c 1";
        else if (temp == "An" && isAlpha90 && isGamma90) SgNameHM = "A 1 n 1";
        else if (temp == "Ia" && isAlpha90 && isGamma90) SgNameHM = "I 1 a 1";
        else if (temp == "Aa" && isAlpha90 && isGamma90) SgNameHM = "A 1 a 1";
        else if (temp == "Cn" && isAlpha90 && isGamma90) SgNameHM = "C 1 n 1";
        else if (temp == "Ic" && isAlpha90 && isGamma90) SgNameHM = "I 1 c 1";
        else if (temp == "Aa" && isAlpha90 && isBeta90) SgNameHM = "A 1 1 a";
        else if (temp == "Bn" && isAlpha90 && isBeta90) SgNameHM = "B 1 1 n";
        else if (temp == "Ib" && isAlpha90 && isBeta90) SgNameHM = "I 1 1 b";
        else if (temp == "Bb" && isAlpha90 && isBeta90) SgNameHM = "B 1 1 b";
        else if (temp == "An" && isAlpha90 && isBeta90) SgNameHM = "A 1 1 n";
        else if (temp == "Ia" && isAlpha90 && isBeta90) SgNameHM = "I 1 1 a";
        else if (temp == "Bb" && isBeta90 && isGamma90) SgNameHM = "B b 1 1";
        else if (temp == "Cn" && isBeta90 && isGamma90) SgNameHM = "C n 1 1";
        else if (temp == "Ic" && isBeta90 && isGamma90) SgNameHM = "I c 1 1";
        else if (temp == "Cc" && isBeta90 && isGamma90) SgNameHM = "C c 1 1";
        else if (temp == "Bn" && isBeta90 && isGamma90) SgNameHM = "B n 1 1";
        else if (temp == "Ib" && isBeta90 && isGamma90) SgNameHM = "I b 1 1";
        else if (temp == "P2/m" && isAlpha90 && isGamma90) SgNameHM = "P 1 2/m 1";
        else if (temp == "P2/m" && isAlpha90 && isBeta90) SgNameHM = "P 1 1 2/m";
        else if (temp == "P2/m" && isBeta90 && isGamma90) SgNameHM = "P 2/m 1 1";
        else if (temp == "P2sub/m" && isAlpha90 && isGamma90) SgNameHM = "P 1 2sub1/m 1";
        else if (temp == "P2sub/m" && isAlpha90 && isBeta90) SgNameHM = "P 1 1 2sub1/m";
        else if (temp == "P2sub/m" && isBeta90 && isGamma90) SgNameHM = "P 2sub1/m 1 1";
        else if (temp == "C2/m" && isAlpha90 && isGamma90) SgNameHM = "C 1 2/m 1";
        else if (temp == "A2/m" && isAlpha90 && isGamma90) SgNameHM = "A 1 2/m 1";
        else if (temp == "I2/m" && isAlpha90 && isGamma90) SgNameHM = "I 1 2/m 1";
        else if (temp == "A2/m" && isAlpha90 && isBeta90) SgNameHM = "A 1 1 2/m";
        else if (temp == "B2/m" && isAlpha90 && isBeta90) SgNameHM = "B 1 1 2/m";
        else if (temp == "I2/m" && isAlpha90 && isBeta90) SgNameHM = "I 1 1 2/m";
        else if (temp == "B2/m" && isBeta90 && isGamma90) SgNameHM = "B 2/m 1 1";
        else if (temp == "C2/m" && isBeta90 && isGamma90) SgNameHM = "C 2/m 1 1";
        else if (temp == "I2/m" && isBeta90 && isGamma90) SgNameHM = "I 2/m 1 1";
        else if (temp == "P2/c" && isAlpha90 && isGamma90) SgNameHM = "P 1 2/c 1";
        else if (temp == "P2/n" && isAlpha90 && isGamma90) SgNameHM = "P 1 2/n 1";
        else if (temp == "P2/a" && isAlpha90 && isGamma90) SgNameHM = "P 1 2/a 1";
        else if (temp == "P2/a" && isAlpha90 && isBeta90) SgNameHM = "P 1 1 2/a";
        else if (temp == "P2/n" && isAlpha90 && isBeta90) SgNameHM = "P 1 1 2/n";
        else if (temp == "P2/b" && isAlpha90 && isBeta90) SgNameHM = "P 1 1 2/b";
        else if (temp == "P2/b" && isBeta90 && isGamma90) SgNameHM = "P 2/b 1 1";
        else if (temp == "P2/n" && isBeta90 && isGamma90) SgNameHM = "P 2/n 1 1";
        else if (temp == "P2/c" && isBeta90 && isGamma90) SgNameHM = "P 2/c 1 1";
        else if (temp == "P2sub1/c" && isAlpha90 && isGamma90) SgNameHM = "P 1 2sub1/c 1";
        else if (temp == "P2sub1/n" && isAlpha90 && isGamma90) SgNameHM = "P 1 2sub1/n 1";
        else if (temp == "P2sub1/a" && isAlpha90 && isGamma90) SgNameHM = "P 1 2sub1/a 1";
        else if (temp == "P2sub1/a" && isAlpha90 && isBeta90) SgNameHM = "P 1 1 2sub1/a";
        else if (temp == "P2sub1/n" && isAlpha90 && isBeta90) SgNameHM = "P 1 1 2sub1/n";
        else if (temp == "P2sub1/b" && isAlpha90 && isBeta90) SgNameHM = "P 1 1 2sub1/b";
        else if (temp == "P2sub1/b" && isBeta90 && isGamma90) SgNameHM = "P 2sub1/b 1 1";
        else if (temp == "P2sub1/n" && isBeta90 && isGamma90) SgNameHM = "P 2sub1/n 1 1";
        else if (temp == "P2sub1/c" && isBeta90 && isGamma90) SgNameHM = "P 2sub1/c 1 1";
        else if (temp == "C2/c" && isAlpha90 && isGamma90) SgNameHM = "C 1 2/c 1";
        else if (temp == "A2/n" && isAlpha90 && isGamma90) SgNameHM = "A 1 2/n 1";
        else if (temp == "I2/a" && isAlpha90 && isGamma90) SgNameHM = "I 1 2/a 1";
        else if (temp == "A2/a" && isAlpha90 && isGamma90) SgNameHM = "A 1 2/a 1";
        else if (temp == "C2/n" && isAlpha90 && isGamma90) SgNameHM = "C 1 2/n 1";
        else if (temp == "I2/c" && isAlpha90 && isGamma90) SgNameHM = "I 1 2/c 1";
        else if (temp == "A2/a" && isAlpha90 && isBeta90) SgNameHM = "A 1 1 2/a";
        else if (temp == "B2/n" && isAlpha90 && isBeta90) SgNameHM = "B 1 1 2/n";
        else if (temp == "I2/b" && isAlpha90 && isBeta90) SgNameHM = "I 1 1 2/b";
        else if (temp == "B2/b" && isAlpha90 && isBeta90) SgNameHM = "B 1 1 2/b";
        else if (temp == "A2/n" && isAlpha90 && isBeta90) SgNameHM = "A 1 1 2/n";
        else if (temp == "I2/a" && isAlpha90 && isBeta90) SgNameHM = "I 1 1 2/a";
        else if (temp == "B2/b" && isBeta90 && isGamma90) SgNameHM = "B 2/b 1 1";
        else if (temp == "C2/n" && isBeta90 && isGamma90) SgNameHM = "C 2/n 1 1";
        else if (temp == "I2/c" && isBeta90 && isGamma90) SgNameHM = "I 2/c 1 1";
        else if (temp == "C2/c" && isBeta90 && isGamma90) SgNameHM = "C 2/c 1 1";
        else if (temp == "B2/n" && isBeta90 && isGamma90) SgNameHM = "B 2/n 1 1";
        else if (temp == "I2/b" && isBeta90 && isGamma90) SgNameHM = "I 2/b 1 1";
        #endregion

        // 261105Cl 追加: 完全記号 (例 'P 21/b 21/n 21/m', 'F 41/d -3 2/m :2', 'R -3 2/m :R') と完全に一致する設定を先に探す。
        //   下の部分一致は短縮記号の表 (SpaceGroupListWithoutSpace) だけを見るため、直方晶・正方晶・立方晶の非標準設定や原点選択 2 を
        //   完全記号で書いた CIF は一致せず、番号による救済で標準設定 (原点選択 1) になっていた。
        //   同じ完全記号の設定が複数あるときは、菱面体軸/六方晶軸 (セルの形) と原点選択で絞る。それでも残る場合 (Cmme とその ba-c 設定の
        //   ように H-M 記号だけでは区別できない組) は最初の設定を採る。一致しなければ従来どおり下の部分一致に進む。
        var fullName = SgNameHM.Replace(" ", "");
        var fullHits = new List<int>();
        for (int i = 1; i < SymmetryStatic.TotalSpaceGroupNumber; i++)
            //if (SymmetryStatic.Symmetries[i].SpaceGroupHMfullStr.Replace(" ", "") == fullName) // 261111Cl 変更前
            if (SgNameTables.HMfullCompact[i] == fullName) // 261111Cl 変更: 空白を除いた表を 1 回だけ作って比べる
                fullHits.Add(i);
        if (fullHits.Count > 0)
        {
            static string hmOf(int i) => SymmetryStatic.Symmetries[i].SpaceGroupHMStr;
            IEnumerable<int> narrowed = fullHits;
            if (narrowed.Any(i => hmOf(i).Contains("Rho", Ord)))
                narrowed = narrowed.Where(i => hmOf(i).Contains(isRhomboShape ? "Rho" : "Hex", Ord));
            if (narrowed.Any(i => hmOf(i).Contains("(1)", Ord) || hmOf(i).Contains("(2)", Ord)))
                narrowed = narrowed.Where(i => hmOf(i).Contains(IsOrigineChoice2 ? "(2)" : "(1)", Ord));
            symmetrySeriesNumber = narrowed.DefaultIfEmpty(fullHits[0]).First();
        }

        //文字列を含んでいて、かつ、文字数が少ない空間群を選択する (C1とC121などを見分けるため)

        int length = int.MaxValue;
        var sg = SymmetryStatic.SpaceGroupListWithoutSpace;
        Replace(ref SgNameHM, " ", "");
        for (int i = 0; fullHits.Count == 0 && i < SymmetryStatic.TotalSpaceGroupNumber; i++) // 261105Cl: 完全記号で決まったときは部分一致を行わない
        {
            if (sg[i].Length < length)
            {
                var sg_list = new List<string>([sg[i]]);

                if (sg[i].Contains('e'))
                {
                    sg_list.Add(sg[i].Replace('e', 'a'));
                    sg_list.Add(sg[i].Replace('e', 'b'));
                    sg_list.Add(sg[i].Replace('e', 'c'));
                }

                if (sg_list.Any(s => s.IndexOf(SgNameHM, 0) != -1))
                {
                    symmetrySeriesNumber = i;
                    length = sg[i].Length;
                }
            }
        }

        //Rhombohedoralのときの処置
        //if (isRhomboShape && SymmetryStatic.Symmetries[symmetrySeriesNumber].SpaceGroupHMStr.Contains("Hex", StringComparison.Ordinal))
        if (isRhomboShape && symmetrySeriesNumber >= 0 && SymmetryStatic.Symmetries[symmetrySeriesNumber].SpaceGroupHMStr.Contains("Hex", StringComparison.Ordinal)) // 260718Cl: HM名不一致 (-1) のとき Symmetries[-1] で例外になり、番号による救済 (下の SpaceGroupNumber ループ) が死んでいたため -1 ガードを追加 (3箇所)
            symmetrySeriesNumber++;

        //originChoiceが2のときの対処
        //if (IsOrigineChoice2 && SymmetryStatic.Symmetries[symmetrySeriesNumber].SpaceGroupHMStr.Contains("(1)"))
        if (IsOrigineChoice2 && symmetrySeriesNumber >= 0 && SymmetryStatic.Symmetries[symmetrySeriesNumber].SpaceGroupHMStr.Contains("(1)"))
            symmetrySeriesNumber++;

        if (SpaceGroupNumber >= 1 && SpaceGroupNumber <= 230)
            //if (symmetrySeriesNumber >= 1 || SpaceGroupNumber != SymmetryStatic.Symmetries[symmetrySeriesNumber].SpaceGroupNumber)
            //if (symmetrySeriesNumber == -1 || symmetrySeriesNumber >= 1 || SpaceGroupNumber != SymmetryStatic.Symmetries[symmetrySeriesNumber].SpaceGroupNumber) // 261105Cl 変更前
            // 261105Cl 修正: 条件の `symmetrySeriesNumber >= 1 ||` により、H-M 記号で設定が見つかっても常に「その番号の最初の設定」で
            //   上書きしていた (Hall 記号の無い CIF で P 1 21/n 1 → P 1 21/c 1、P b n m → Pnma など、530 設定中 299 設定)。
            //   CIF の対称操作による後段の確認は原点の移動しか探さないため、誤った設定のまま原子が展開されていた。
            //   H-M 記号で見つからない (-1, 0) か、見つかった設定の番号が CIF の番号と食い違うときだけ、番号で救済する。
            if (symmetrySeriesNumber <= 0 || SpaceGroupNumber != SymmetryStatic.Symmetries[symmetrySeriesNumber].SpaceGroupNumber)
                for (int i = 0; i < SymmetryStatic.TotalSpaceGroupNumber; i++)
                    if (SymmetryStatic.Symmetries[i].SpaceGroupNumber == SpaceGroupNumber)
                    {
                        symmetrySeriesNumber = i;
                        break;
                    }

        return symmetrySeriesNumber;
    }
    #endregion

    #region CIFファイルへの変換

    /// <summary>261105Cl 追加: CIF の数値の文字列 ("0.071(1)"、"-.0030(9)"、"1.2E-3(3)"、"1.2(3)E-3" など。括弧は仮数の最後の桁の標準不確かさ) に係数を掛けた、
    /// 同じ形式の文字列を返す (不確かさも同じ係数)。Atoms2 のコンストラクタに CIF の生の値として渡す。
    /// Crystal2.Decompose は内部形式 (|誤差|) 用で CIF の括弧を読めないので使わない。指数を使わない InvariantCulture の書式で書く。
    /// 読めない文字列は換算できないので空 (= 未設定) を返す (未換算のまま別の型として採用しない)。
    /// 値は、不確かさの有効数字 3〜4 桁目の位まで出す (不確かさの桁で丸めて値の精度を落とさないため)。
    /// offset (261105Cl 追加) を足すと、原点の移動にも使える。係数 1 で指数が無いときは、元の小数の桁と不確かさの表記をそのまま使う。</summary>
    //private static string scaleValueText(string s, double factor) // 261105Cl 変更前
    private static string scaleValueText(string s, double factor, double offset = 0)
    {
        if (string.IsNullOrWhiteSpace(s))
            return s;
        // 不確かさは指数の前 ("7.896(13)E-1") と後ろ (CIF 1.1 の文法どおりの "7.896E-1(13)") のどちらでも読む。どちらも仮数の最後の桁に対する値
        //var m = Regex.Match(s.Trim(), @"^([+-]?(?:\d+\.?\d*|\.\d+))(?:\((\d+)\))?(?:[eE]([+-]?\d+))?(?:\((\d+)\))?$"); // 261111Cl 変更前
        var m = CifNumberWithUncertaintyRegex().Match(s.Trim()); // 261111Cl 変更: ソース生成の Regex
        if (!m.Success || (m.Groups[2].Success && m.Groups[4].Success))
            return "";
        var ci = System.Globalization.CultureInfo.InvariantCulture;
        var mantissa = m.Groups[1].Value;
        var v = double.Parse(mantissa, System.Globalization.NumberStyles.Float, ci);
        var decimals = mantissa.Contains('.') ? mantissa.Length - mantissa.IndexOf('.') - 1 : 0;
        var su = m.Groups[2].Success ? m.Groups[2] : m.Groups[4];
        var e = su.Success ? long.Parse(su.Value, ci) * Math.Pow(10, -decimals) : double.NaN;
        if (m.Groups[3].Success)
        {
            var p = Math.Pow(10, int.Parse(m.Groups[3].Value, ci));
            v *= p;
            e *= p;
        }
        //v *= factor; // 261105Cl 変更前
        v = v * factor + offset; // 261105Cl 変更: 原点の移動 (offset) も
        e *= Math.Abs(factor);
        if (factor == 1 && !m.Groups[3].Success) // 261105Cl 追加: 平行移動だけなら元の桁で書く (0.336 + 0.25 を 0.5860000000000001 にしない)
        {
            var outDecimals = Math.Max(decimals, offset == 0 ? 0 : 3); // 原点の移動は 1/8 刻みなので小数 3 桁あれば表せる ("0" + 0.25 を "0" にしない)
            var suText = su.Success ? "(" + long.Parse(su.Value, ci) * (long)Math.Pow(10, outDecimals - decimals) + ")" : "";
            return v.ToString("F" + outDecimals, ci) + suText;
        }
        if (double.IsNaN(e) || e <= 0)
            return v.ToString("0.####################", ci);
        var digits = Math.Clamp((int)Math.Ceiling(-Math.Log10(e)) + 3, 0, 20); // 不確かさの有効数字 3〜4 桁まで出す (値を不確かさの桁で丸めて精度を落とさないため)
        var errDigits = (long)Math.Round(e * Math.Pow(10, digits));
        return v.ToString("F" + digits, ci) + "(" + errDigits + ")";
    }

    public static string ConvertToCIF(Crystal crystal)
    {
        var sb = new StringBuilder();
        //sb.AppendLine("# This file is exported from \"" + Process.GetCurrentProcess().ProcessName + "\""); // 261111Cl 変更前: Process を Dispose していなかった
        sb.AppendLine("# This file is exported from \"" + Path.GetFileNameWithoutExtension(Environment.ProcessPath) + "\""); // 261111Cl 変更: 実行ファイル名 (拡張子なし) = 旧の ProcessName と同じ
        sb.AppendLine("# https://github.com/seto77/");

        sb.AppendLine("data_global");
        sb.AppendLine("_chemical_name '" + crystal.Name + "'");

        if (crystal.PublAuthorName != string.Empty)
        {
            sb.AppendLine("loop_");
            sb.AppendLine("_publ_author_name");
            foreach (string str in crystal.PublAuthorName.Split(',', true))
                sb.AppendLine("'" + str.Trim() + "'");
        }

        if(crystal.Journal != string.Empty)
            sb.AppendLine("_journal_name '" + crystal.Journal + "'");

        #region 論文タイトル
        if (crystal.PublSectionTitle != "")
        {
            sb.AppendLine("_publ_section_title");
            sb.AppendLine(";");
            string title = "";
            foreach (string t in crystal.PublSectionTitle.Split(' ', true))
            {
                if ((title + " " + t).Length > 80)
                {
                    sb.AppendLine(title);
                    title = "";
                }
                title += " " + t;
            }
            if (title != "")
                sb.AppendLine(title);
            sb.AppendLine(";");
        }
        #endregion

        #region 格子定数、対称性
        sb.AppendLine("_chemical_formula_sum '" + crystal.ChemicalFormulaSum + "'");
        //sb.AppendLine("_cell_length_a " + (crystal.A * 10).ToString("f6")); // 261105Cl 変更前 (8 行)
        //sb.AppendLine("_cell_length_b " + (crystal.B * 10).ToString("f6"));
        //sb.AppendLine("_cell_length_c " + (crystal.C * 10).ToString("f6"));
        //sb.AppendLine("_cell_angle_alpha " + (crystal.Alpha / Math.PI * 180).ToString("f6"));
        //sb.AppendLine("_cell_angle_beta " + (crystal.Beta / Math.PI * 180).ToString("f6"));
        //sb.AppendLine("_cell_angle_gamma " + (crystal.Gamma / Math.PI * 180).ToString("f6"));
        //sb.AppendLine("_cell_volume " + (crystal.Volume * 1000).ToString("f6"));
        //sb.AppendLine("_exptl_crystal_density_diffrn " + crystal.Density.ToString("f6"));
        // 261105Cl 変更: 数値は InvariantCulture で書く (小数点がカンマの環境では、正当な CIF の数値にならなかった)。原子の 2 つのループも同じ
        var ci = System.Globalization.CultureInfo.InvariantCulture;
        sb.AppendLine(ci, $"_cell_length_a {crystal.A * 10:f6}");
        sb.AppendLine(ci, $"_cell_length_b {crystal.B * 10:f6}");
        sb.AppendLine(ci, $"_cell_length_c {crystal.C * 10:f6}");
        sb.AppendLine(ci, $"_cell_angle_alpha {crystal.Alpha / Math.PI * 180:f6}");
        sb.AppendLine(ci, $"_cell_angle_beta {crystal.Beta / Math.PI * 180:f6}");
        sb.AppendLine(ci, $"_cell_angle_gamma {crystal.Gamma / Math.PI * 180:f6}");
        sb.AppendLine(ci, $"_cell_volume {crystal.Volume * 1000:f6}");
        sb.AppendLine(ci, $"_exptl_crystal_density_diffrn {crystal.Density:f6}");

        var sym = crystal.Symmetry;
        sb.AppendLine("_space_group_IT_number " + sym.SpaceGroupNumber);
        sb.AppendLine("_symmetry_Int_Tables_number " + sym.SpaceGroupNumber);
        sb.AppendLine("_symmetry_cell_setting '" + sym.CrystalSystemStr + "'");
        var hm = sym.SpaceGroupHMStr;
        hm = hm.Replace("Hex", "");
        hm = hm.Replace("Rho", "");
        sb.AppendLine("_symmetry_space_group_name_H-M '" + hm + "'");
        //sb.AppendLine("_symmetry_space_group_name_Hall '" + sym.SpaceGroupHallStr + "'"); // 261105Cl 変更前
        // 261105Cl 修正: 空白なしの Hall 記号は標準外で、P312 が 'P32' (= P3_2) と読まれるなど他のソフトで誤解釈される。標準表記を書き出す。
        sb.AppendLine("_symmetry_space_group_name_Hall '" + sym.SpaceGroupHallStandardStr + "'");
        #endregion

        #region 原子の等価位置
        sb.AppendLine("loop_");
        sb.AppendLine("_symmetry_equiv_pos_as_xyz");
        bool[][] flag = [];
        if (sym.LatticeTypeStr == "P") flag = [[false, false, false]];
        else if (sym.LatticeTypeStr == "A") flag = [[false, false, false], [false, true, true]];
        else if (sym.LatticeTypeStr == "B") flag = [[false, false, false], [true, false, true]];
        else if (sym.LatticeTypeStr == "C") flag = [[false, false, false], [true, true, false]];
        else if (sym.LatticeTypeStr == "I") flag = [[false, false, false], [true, true, true]];
        else if (sym.LatticeTypeStr == "F") flag = [[false, false, false], [false, true, true], [true, false, true], [true, true, false]];

        foreach (string wp in SymmetryStatic.WyckoffPositions[crystal.SymmetrySeriesNumber][0].PositionStr)
        {
            if (sym.SpaceGroupHMsubStr != "H")
            {
                for (int i = 0; i < flag.Length; i++)
                {
                    string[] xyz = wp.Split(',', true);
                    for (int j = 0; j < flag[i].Length; j++)
                    {
                        if (flag[i][j])
                        {
                            if (xyz[j].EndsWith("+1/2", Ord)) xyz[j] = xyz[j].Replace("+1/2", "");
                            else if (xyz[j].EndsWith("+1/4", Ord)) xyz[j] = xyz[j].Replace("+1/4", "+3/4");
                            else if (xyz[j].EndsWith("+3/4", Ord)) xyz[j] = xyz[j].Replace("+3/4", "+1/4");
                            else xyz[j] += "+1/2";
                        }
                    }
                    sb.AppendLine($"  '{xyz[0]},{xyz[1]},{xyz[2]}'");
                }
            }
            else//R格子のHexセッティングのとき
            {
                sb.AppendLine("  '" + wp + "'");//(0,0,0)
                                                //(1/3,2/3,2/3)
                string[] xyz = wp.Split(',', true);
                xyz[0] += "+1/3";
                xyz[1] += "+2/3";
                if (xyz[2].EndsWith("+1/2", Ord)) xyz[2] = xyz[2].Replace("+1/2", "+1/6");
                else xyz[2] += "+2/3";
                sb.AppendLine($"  '{xyz[0]},{xyz[1]},{xyz[2]}'");
                //(2/3,1/3,1/3)
                xyz = wp.Split(',', true);
                xyz[0] += "+2/3";
                xyz[1] += "+1/3";
                if (xyz[2].EndsWith("+1/2", Ord)) xyz[2] = xyz[2].Replace("+1/2", "+5/6");
                else xyz[2] += "+1/3";
                sb.AppendLine($"  '{xyz[0]},{xyz[1]},{xyz[2]}'");
            }
        }
        #endregion

        #region 各原子の情報
        sb.AppendLine("loop_");
        sb.AppendLine("_atom_site_label");
        sb.AppendLine("_atom_site_type_symbol");
        sb.AppendLine("_atom_site_fract_x");
        sb.AppendLine("_atom_site_fract_y");
        sb.AppendLine("_atom_site_fract_z");
        sb.AppendLine("_atom_site_occupancy");
        sb.AppendLine("_atom_site_U_iso_or_equiv");

        foreach (var a in crystal.Atoms)
        {
            //var u = double.IsNaN(a.Dsf.Uiso) ? 0 : a.Dsf.Uiso * 100; // 261105Cl 変更前
            // 261105Cl 変更: 異方性の原子は、テンソルから求めた等価な等方 U (U_eq = B_eq/8π²、Hamilton / Fischer & Tillmanns) を書く。
            // 以前は保存されている Iso をそのまま書いたため、異方性だけ入力した原子 (Iso=0) では 0 になっていた
            var uRaw = a.Dsf.UseIso ? a.Dsf.Uiso : a.Dsf.BisoEffective / (8 * Math.PI * Math.PI);
            var u = double.IsNaN(uRaw) ? 0 : uRaw * 100;
            //sb.AppendLine($" {a.Label} {AtomStatic.AtomicName(a.AtomicNumber)} {a.X:f6} {a.Y:f6} {a.Z:f6} {a.Occ:f6} {u:f6}"); // 261105Cl 変更前
            sb.AppendLine(ci, $" {a.Label} {AtomStatic.AtomicName(a.AtomicNumber)} {a.X:f6} {a.Y:f6} {a.Z:f6} {a.Occ:f6} {u:f6}"); // 261105Cl 変更: InvariantCulture
        }


        if (crystal.Atoms.Any(a => !a.Dsf.UseIso))
        {
            sb.AppendLine("loop_");
            sb.AppendLine("_atom_site_aniso_label");
            sb.AppendLine("_atom_site_aniso_U_11");
            sb.AppendLine("_atom_site_aniso_U_22");
            sb.AppendLine("_atom_site_aniso_U_33");
            sb.AppendLine("_atom_site_aniso_U_23");
            sb.AppendLine("_atom_site_aniso_U_13");
            sb.AppendLine("_atom_site_aniso_U_12");
            foreach (var a in crystal.Atoms.Where(e => !e.Dsf.UseIso))
            {
                var u11 = double.IsNaN(a.Dsf.U11) ? 0 : a.Dsf.U11 * 100;
                var u22 = double.IsNaN(a.Dsf.U22) ? 0 : a.Dsf.U22 * 100;
                var u33 = double.IsNaN(a.Dsf.U33) ? 0 : a.Dsf.U33 * 100;
                var u23 = double.IsNaN(a.Dsf.U23) ? 0 : a.Dsf.U23 * 100;
                var u31 = double.IsNaN(a.Dsf.U31) ? 0 : a.Dsf.U31 * 100;
                var u12 = double.IsNaN(a.Dsf.U12) ? 0 : a.Dsf.U12 * 100;
                //sb.AppendLine($" {a.Label} {u11:f6} {u22:f6} {u33:f6} {u23:f6} {u31:f6} {u12:f6}"); // 261105Cl 変更前
                sb.AppendLine(ci, $" {a.Label} {u11:f6} {u22:f6} {u33:f6} {u23:f6} {u31:f6} {u12:f6}"); // 261105Cl 変更: InvariantCulture
            }
        }
        #endregion

        return sb.ToString();
    }
    #endregion

    #region OpenGLのためのBond設定
    public static void SetOpenGL_property(Crystal c)
    {
        foreach (Atoms a in c.Atoms)
            a.ResetVesta();
        c.Bonds = Bonds.GetVestaBonds(c.Atoms.Select(a => a.ElementName).Distinct());
    }
    #endregion
}
