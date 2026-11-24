using MathNet.Numerics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace Crystallography;

//public static class XYFile // 261111Cl 変更前
public static partial class XYFile // 261111Cl 変更: [GeneratedRegex] のため partial
{
    /// <summary>261111Cl 追加: 旧版の要素名 &lt;pt&gt; / &lt;/pt&gt; / &lt;pt /&gt; (現 Profile.Pt) を読むための置換の対象。旧は 1 行ずつ Regex.Replace(tempstr, @"(&lt;/?)pt(?=[\s/&gt;])", "$1Pt") していた。
    /// 全文に 1 回だけ掛けるので、後ろを見る文字から改行 (\r, \n) を除く (1 行ずつのときは行末の後ろに文字が無く、一致しなかった)</summary>
    [System.Text.RegularExpressions.GeneratedRegex(@"(</?)pt(?=[\s/>-[\r\n]])")]
    private static partial System.Text.RegularExpressions.Regex OldPtTagRegex();

    /// <summary>261124Cl 追加: 保存先と同じフォルダの一時ファイルへ書き切ってから、保存先と置き換える。
    /// 書き込みの途中で失敗しても (ディスクが満杯・直列化の例外など) 元のファイルは残る。失敗したら一時ファイルを消して false を返す。
    /// (PDIndexer の査読 C-03: 旧実装は FileMode.Create で先に切り詰めていたので、失敗すると元の内容が消えていた)</summary>
    public static bool SaveReplacing(string fileName, Action<Stream> write)
    {
        string tmp = null;
        bool created = false;
        try
        {
            // 261124Cl (codex2 の査読): 保存先の絶対パスを 1 回だけ決め、一時ファイルと置換の両方に使う (作業フォルダが途中で変わっても別のファイルを置き換えない)
            var full = Path.GetFullPath(fileName);
            tmp = Path.Combine(Path.GetDirectoryName(full), $".{Guid.NewGuid():N}.tmp");
            using (var fs = new FileStream(tmp, FileMode.CreateNew))
            {
                created = true;
                write(fs);
            }
            File.Move(tmp, full, overwrite: true);
            return true;
        }
        catch
        {
            try { if (created && File.Exists(tmp)) File.Delete(tmp); } catch { } // 自分で作った一時ファイルだけ消す
            return false;
        }
    }

    //public static bool SavePdi2File(DiffractionProfile2[] dp, string fileName) // 261124Cl 変更前: FileMode.Create で先に切り詰めてから直列化していた (C-03)
    //{
    //    FileStream fs = null;
    //    try
    //    {
    //        System.Xml.Serialization.XmlSerializer serializer = new System.Xml.Serialization.XmlSerializer(typeof(DiffractionProfile2[]));
    //        fs = new FileStream(fileName, FileMode.Create);
    //        serializer.Serialize(fs, dp);
    //        fs.Close();
    //        return true;
    //    }
    //    catch
    //    {
    //        fs?.Close();
    //        return false;
    //    }
    //}
    public static bool SavePdi2File(DiffractionProfile2[] dp, string fileName) // 261124Cl 変更: 書き切ってから置き換える
        => SaveReplacing(fileName, fs => new System.Xml.Serialization.XmlSerializer(typeof(DiffractionProfile2[])).Serialize(fs, dp));

    public static DiffractionProfile2[] ReadPdi2File(string fileName, int version = 2)
    {
        //OriginalFormatType  -> SrcAxisMode
        //OriginalWaveLength  -> SrcWaveLength
        //OriginalTakeoffAngle　-> SrcTakeoffAngle
        // 261105Cl 修正: 開くだけで元ファイルを書き換えていた。全行に Replace("pt", "Pt") を掛けて同じファイルへ上書き保存するため、
        //  名前やコメント中の "pt" まで "Pt" に変わっていた (例: "Sample opt" → "Sample oPt")。
        //  旧書式との互換のための置換はメモリ上だけで行い、ファイルへは書き戻さない。
        //  "pt" → "Pt" は旧版の要素名 <pt> (現 Profile.Pt) を読むためのものなので、要素タグ <pt>, </pt>, <pt /> だけに限る。
        string xmlText = null; // 261105Cl 追加: 互換の置換を済ませた XML 本文 (ファイルは書き換えない)
        try//まずXMLのタグを変更
        {
            using var reader = new StreamReader(fileName, Encoding.GetEncoding("UTF-8"));
            var strList = new List<string>();
            string tempstr;
            while ((tempstr = reader.ReadLine()) != null)
            {
                //tempstr = tempstr.Replace("OriginalFormatType", "SrcAxisMode"); // 261111Cl 変更前 (5 行): 1 行ごとに置換していた。下で全文に 1 回だけ掛ける
                //tempstr = tempstr.Replace("OriginalWaveLength", "SrcWaveLength");
                //tempstr = tempstr.Replace("OriginalTakeoffAngle", "SrcTakeoffAngle");
                ////tempstr = tempstr.Replace("pt", "Pt"); // 261105Cl 変更前
                //tempstr = System.Text.RegularExpressions.Regex.Replace(tempstr, @"(</?)pt(?=[\s/>])", "$1Pt"); // 261105Cl
                strList.Add(tempstr);
            }

            reader.Close();

            //using var writer = new StreamWriter(fileName, false, Encoding.GetEncoding("UTF-8")); // 261105Cl 変更前: 元ファイルへ書き戻していた
            //for (int i = 0; i < strList.Count; i++)
            //    writer.WriteLine(strList[i]);
            //writer.Flush();
            //writer.Close();
            xmlText = string.Join(Environment.NewLine, strList); // 261105Cl 追加
            // 261111Cl 追加: 互換の置換を全文に 1 回だけ掛ける (順も同じ)。どの置換も改行を含まない語なので行をまたいで一致せず、1 行ずつ掛けたのと同じ結果になる
            xmlText = xmlText.Replace("OriginalFormatType", "SrcAxisMode").Replace("OriginalWaveLength", "SrcWaveLength").Replace("OriginalTakeoffAngle", "SrcTakeoffAngle");
            xmlText = OldPtTagRegex().Replace(xmlText, "$1Pt");
        }
        catch { };


        System.IO.FileStream fs = null;
        try
        {
            if (version == 2)
            {
                var serializer = new System.Xml.Serialization.XmlSerializer(typeof(DiffractionProfile2[]));
                //fs = new FileStream(fileName, FileMode.Open); // 261105Cl 変更前
                //DiffractionProfile2[] dp = (DiffractionProfile2[])serializer.Deserialize(fs);
                //fs.Close();
                DiffractionProfile2[] dp = (DiffractionProfile2[])serializer.Deserialize(new StringReader(xmlText)); // 261105Cl: 置換済みの本文から読む
                if (dp.Length > 0)
                    return dp;
                else
                    return [];
            }
            else //バージョン1
            {
                var serializer = new System.Xml.Serialization.XmlSerializer(typeof(DiffractionProfile[]));
                //fs = new FileStream(fileName, FileMode.Open); // 261105Cl 変更前
                //DiffractionProfile[] dp = (DiffractionProfile[])serializer.Deserialize(fs);
                //fs.Close();
                DiffractionProfile[] dp = (DiffractionProfile[])serializer.Deserialize(new StringReader(xmlText)); // 261105Cl: 置換済みの本文から読む

                if (dp.Length > 0)
                    return dp.Select(e => e.ConvertToDiffractionProfile2()).ToArray();
                else
                    return Array.Empty<DiffractionProfile2>();
            }

        }
        catch//もしシリアライズできなかったら、name部分に間違った日本語が書かれている可能性あり。
        {
            fs?.Close();
            try
            {
                var serializer = new System.Xml.Serialization.XmlSerializer(typeof(DiffractionProfile[]));
                // 261105Cl 変更前: ここでも行末の '\0' を除いた内容で元ファイルを上書きしていた
                //var reader = new StreamReader(fileName, Encoding.UTF8);
                //var strList = new List<string>();
                //string tempstr;
                //while ((tempstr = reader.ReadLine()) != null)
                //    strList.Add(tempstr);
                //reader.Close();
                //StreamWriter writer = new StreamWriter(fileName, false, Encoding.UTF8);
                //for (int i = 0; i < strList.Count; i++)
                //{
                //    if (strList[i].EndsWith("\0"))
                //        writer.WriteLine(strList[i].TrimEnd('\0'));
                //    else
                //        writer.WriteLine(strList[i]);
                //}
                //
                //writer.Flush();
                //writer.Close();
                //fs = new System.IO.FileStream(fileName, System.IO.FileMode.Open);
                //var dp = (DiffractionProfile[])serializer.Deserialize(fs);
                //fs.Close();
                // 261105Cl 変更: 同じ処理 (各行の行末の '\0' を除く) をメモリ上で行い、そこから読む。元ファイルは書き換えない。
                //  旧コードは先頭の置換を書き戻した後のファイルを読んでいたので、置換済みの xmlText を使えば同じ内容になる。
                var trimmed = string.Join(Environment.NewLine, xmlText.Split(Environment.NewLine).Select(s => s.TrimEnd('\0'))); // 261105Cl
                var dp = (DiffractionProfile[])serializer.Deserialize(new StringReader(trimmed)); // 261105Cl
                if (dp.Length > 0)
                    return dp.Select(e => e.ConvertToDiffractionProfile2()).ToArray();
                else return Array.Empty<DiffractionProfile2>();
            }
            catch
            {
                try
                {
                    if (fs != null)
                        fs.Close();
                    //var reader = new StreamReader(fileName, Encoding.GetEncoding("UTF-8")); // 261111Cl 変更前
                    using var reader = new StreamReader(fileName, Encoding.GetEncoding("UTF-8")); // 261111Cl 変更: 例外でも閉じる
                    var strList = new List<string>();
                    string tempstr;
                    while ((tempstr = reader.ReadLine()) != null)
                        strList.Add(tempstr);
                    reader.Close();
                    if (strList.Count <= 3)
                        return Array.Empty<DiffractionProfile2>();

                    var diffProf = new DiffractionProfile2();

                    //古いヘッダの書式
                    //Wave Length (0.1nm):0.4176811455              0
                    //Camera Length (mm):445.8055943768             1
                    //Pixel Size Horizontal(mm):0.10080533          2
                    //AspectRatio(Vertical/Horizontal):1            3
                    //Mode:Angle                                    4
                    //4,733.102005852614                            5

                    if (strList[0].Contains("Wave Length", StringComparison.Ordinal))
                    {
                        if (strList[0].Contains("(nm)", StringComparison.Ordinal))
                            //diffProf.SrcProperty.WaveLength = Convert.ToDouble((strList[0].Split(':', true))[1]);
                            diffProf.SrcProperty.WaveLength = (strList[0].Split(':', true))[1].ToDouble();
                        else if (strList[0].Contains("(0.1nm)", StringComparison.Ordinal))
                            //diffProf.SrcProperty.WaveLength = Convert.ToDouble((strList[0].Split(':', true))[1]) / 10.0;
                            diffProf.SrcProperty.WaveLength = (strList[0].Split(':', true))[1].ToDouble() / 10.0;
                    }

                    if ((strList[4].Split(':', true))[1] == "Angle")
                        diffProf.SrcProperty.AxisMode = HorizontalAxis.Angle;
                    else if ((strList[4].Split(':', true))[1] == "d-spacing")
                        diffProf.SrcProperty.AxisMode = HorizontalAxis.d;
                    else if ((strList[4].Split(':', true))[1] == "Energy")
                        diffProf.SrcProperty.AxisMode = HorizontalAxis.EnergyXray;
                    else
                        return Array.Empty<DiffractionProfile2>();

                    for (int i = 5; i < strList.Count; i++)
                    {
                        string[] str = strList[i].Split(',', true);
                        //diffProf.SourceProfile.Pt.Add(new PointD(Convert.ToDouble(str[0]), Convert.ToDouble(str[1])));
                        diffProf.SourceProfile.Pt.Add(new PointD(str[0].ToDouble(), str[1].ToDouble()));
                    }
                    diffProf.Name = fileName.Remove(0, fileName.LastIndexOf('\\') + 1);
                    return [diffProf];
                }
                catch { return []; }
            }
        }
    }

    public static DiffractionProfile2[] ReadRasFile(string fileName)
    {
        var strArray = new List<string>();
        //var reader = new StreamReader(fileName, Encoding.GetEncoding("UTF-8")); // 261111Cl 変更前
        using var reader = new StreamReader(fileName, Encoding.GetEncoding("UTF-8")); // 261111Cl 変更: 例外でも閉じる
        string tempstr;
        while ((tempstr = reader.ReadLine()) != null)
            strArray.Add(tempstr);
        reader.Close();
        if (strArray.Count <= 3)
            return [];

        var dp = new List<DiffractionProfile2>();

        for (int i = 0; i < strArray.Count; i++)
        {
            if (strArray[i] == "*RAS_INT_START")
            {
                i++;
                dp.Add(new DiffractionProfile2());
                for (; i < strArray.Count; i++)
                {
                    if (strArray[i] != "*RAS_INT_END")
                    {
                        string[] tempStr = strArray[i].Split(new[] { ' ' });
                        //double x = Convert.ToDouble(tempStr[0]);
                        //double y = Convert.ToDouble(tempStr[1]);
                        double x = tempStr[0].ToDouble();
                        double y = tempStr[1].ToDouble();
                        dp[^1].SourceProfile.Pt.Add(new PointD(x, y));
                    }
                    else
                        break;
                }
            }
        }

        for (int i = 0; i < dp.Count; i++)
            dp[i].Name = $"{Path.GetFileName(fileName)}{(dp.Count > 1 ? $" -{i}" : "")}";

        if (dp.Count > 0)
            return [.. dp];
        else return [];
    }

    /// <summary>複数のプロファイルを含むCSV形式(PDI独自形式)</summary>
    /// <param name="fileName"></param>
    /// <returns></returns>
    public static DiffractionProfile2[] ReadCSVFile(string fileName)
    {
        var strArray = new List<string>();
        using (var reader = new StreamReader(fileName, Encoding.GetEncoding("UTF-8")))
        {
            string tempstr;
            while ((tempstr = reader.ReadLine()) != null)
                strArray.Add(tempstr);
        }

        // 261124Cl 変更 (PDIndexer の査読 C-02): 空欄を消して詰めてから 2 列ずつ読んでいたので、短い試料の後ろの行 (,,,4,40,,) の値が
        //  先頭の試料に入った。見出し行の X・Y の列の位置で各試料を読み、その試料の 2 欄が両方空の行は飛ばす。
        //  あわせて、データが 1 行だけのファイル、小文字の見出し (IPAnalyzer の結合 CSV の x,y)、空の試料名を受ける
        List<string> axisHeader = strArray.Count >= 3 ? SplitCsvLine(strArray[1]) : null; // 261124Cl (codex2 の査読): 見出しも引用符を外してから判定する ("X","Y",...)
        if (axisHeader != null && axisHeader.Count >= 2 && axisHeader[0].Trim().Equals("X", StringComparison.OrdinalIgnoreCase) && axisHeader[1].Trim().Equals("Y", StringComparison.OrdinalIgnoreCase))
        {
            try
            {
                var title = SplitCsvLine(strArray[0]);
                var axis = axisHeader;
                var cols = new List<int>();
                for (int c = 0; c + 1 < axis.Count; c++)
                    if (axis[c].Trim().Equals("X", StringComparison.OrdinalIgnoreCase) && axis[c + 1].Trim().Equals("Y", StringComparison.OrdinalIgnoreCase))
                        cols.Add(c++);
                if (cols.Count == 0) return [];
                var dp2 = new DiffractionProfile2[cols.Count];
                for (int k = 0; k < cols.Count; k++)
                {
                    var name = cols[k] < title.Count ? title[cols[k]].Trim() : "";
                    dp2[k] = new DiffractionProfile2 { Name = name.Length > 0 ? name : $"{Path.GetFileName(fileName)}{(cols.Count > 1 ? $" -{k}" : "")}" };
                }
                for (int i = 2; i < strArray.Count; i++)
                {
                    var v = SplitCsvLine(strArray[i]);
                    for (int k = 0; k < cols.Count; k++)
                    {
                        int c = cols[k];
                        var xs = c < v.Count ? v[c].Trim() : "";
                        var ys = c + 1 < v.Count ? v[c + 1].Trim() : "";
                        if (xs.Length == 0 || ys.Length == 0) continue; // この試料はこの行で欠測 (短い試料の後ろ)
                        dp2[k].SourceProfile.Pt.Add(new PointD(xs.ToDouble(), ys.ToDouble()));
                    }
                }
                return dp2;
            }
            catch
            {
                return [];
            }
        }
        return [];

        #region 261124Cl 変更前 (空欄を消して詰めて読んでいた)
        //if (strArray.Count <= 3)
        //    return [];
//
        //if (!strArray[1].StartsWith("X,Y,"))
        //    return [];
//
        //try
        //{
        //    var title = strArray[0].Split([','], StringSplitOptions.RemoveEmptyEntries);
        //    var axis = strArray[1].Split([','], StringSplitOptions.RemoveEmptyEntries);
        //    var value = strArray[2].Split([','], StringSplitOptions.RemoveEmptyEntries);
        //    if (title.Length * 2 != axis.Length || axis.Length != value.Length)
        //        return [];
//
        //    var dp = new DiffractionProfile2[title.Length];
        //    for (int i = 0; i < dp.Length; i++)
        //    {
        //        dp[i] = new DiffractionProfile2();
        //        dp[i].Name = title[i];
        //    }
        //    for (int i = 2; i < strArray.Count; i++)
        //    {
        //        value = strArray[i].Split(new[] { ',' }, StringSplitOptions.RemoveEmptyEntries);
        //        for (int j = 0; j < value.Length / 2; j++)
        //        {
        //            //double x = Convert.ToDouble(value[j * 2]);
        //            //double y = Convert.ToDouble(value[j * 2 + 1]);
        //            var x = value[j * 2].ToDouble();
        //            var y = value[j * 2 + 1].ToDouble();
        //            dp[j].SourceProfile.Pt.Add(new PointD(x, y));
        //        }
        //    }
//
        //    for (int i = 0; i < dp.Length; i++)
        //    {
        //        if (title[i].Length == 0)
        //            dp[i].Name = $"{Path.GetFileName(fileName)}{(dp.Length > 1 ? $" -{i}" : "")}";
        //        else
        //            dp[i].Name = title[i];
        //    }
        //    return dp;
        //}
        //catch
        //{
        //    return [];
        //}
        #endregion
    }

    /// <summary>261124Cl 追加: CSV の 1 行を欄に分ける (空欄も残す)。二重引用符で囲んだ欄の中の , と "" (= ") を扱う</summary>
    private static List<string> SplitCsvLine(string line)
    {
        var result = new List<string>();
        var sb = new StringBuilder();
        bool quoted = false;
        for (int i = 0; i < line.Length; i++)
        {
            char ch = line[i];
            if (quoted)
            {
                if (ch == '"')
                {
                    if (i + 1 < line.Length && line[i + 1] == '"') { sb.Append('"'); i++; }
                    else quoted = false;
                }
                else sb.Append(ch);
            }
            else if (ch == '"') quoted = true;
            else if (ch == ',') { result.Add(sb.ToString()); sb.Clear(); }
            else sb.Append(ch);
        }
        result.Add(sb.ToString());
        return result;
    }



    /// <summary>よく分からないファイルを読み込む</summary>
    /// <param name="fileName"></param>
    /// <param name="separater"></param>
    /// <returns></returns>
    public static DiffractionProfile2 ConvertUnknownFileToProfileData(string fileName, char separater)
    {
        var strArray = new List<string>();
        //var reader = new StreamReader(fileName, Encoding.GetEncoding("UTF-8")); // 261111Cl 変更前
        using var reader = new StreamReader(fileName, Encoding.GetEncoding("UTF-8")); // 261111Cl 変更: 例外でも閉じる
        string tempstr;
        while ((tempstr = reader.ReadLine()) != null)
            strArray.Add(tempstr);
        reader.Close();
        if (strArray.Count <= 3)
            return null;

        var stringList = new List<string[]>();
        //まず指定されたセパレータで全てを区切る
        for (int i = 0; i < strArray.Count; i++)
            stringList.Add(strArray[i].Split([separater], StringSplitOptions.RemoveEmptyEntries));
        //その全てを数値に変換する
        var doubleList = new List<double[]>();
        for (int i = 0; i < stringList.Count; i++)
        {
            var doubleTemp = new List<double>();
            for (int j = 0; j < stringList[i].Length; j++)
            {
                //var str = Miscellaneous.IsDecimalPointComma ? stringList[i][j].Replace('.', ',') : stringList[i][j].Replace(',', '.');
                //if (double.TryParse(str, out double result))
                //    doubleTemp.Add(result);
                var val = stringList[i][j].ToDouble();
                if(val.IsFinite())
                    doubleTemp.Add(stringList[i][j].ToDouble());
            }
            doubleList.Add([.. doubleTemp]);
        }
        int count = 1;
        int beforeLength = 0;
        int countMax = int.MinValue;
        int startColumn = 1;
        int endColumn = 0;

        //doubleList中で2つ以上数値を含み、10個以上連続しているものばしょをえらぶ
        for (int i = 0; i < doubleList.Count; i++)
        {
            count = 0;
            beforeLength = doubleList[i].Length;
            if (beforeLength >= 2)
            {
                for (int j = i; j < doubleList.Count && beforeLength == doubleList[j].Length; j++)
                    count++;
                if (countMax < count)
                {
                    countMax = count;
                    startColumn = i;
                    endColumn = i + count - 1;
                }
                i += count - 1;
            }
        }
        if (countMax < 10) return null;//10以下だったらNullをかえす

        if (endColumn + 1 < doubleList.Count)
            doubleList.RemoveRange(endColumn + 1, doubleList.Count - endColumn - 1);
        doubleList.RemoveRange(0, startColumn);

        //X軸をきめる 0.00000001以上のステップで100個以上連続しているものをさがす
        int xRow = -1;
        double tempStep;
        countMax = int.MinValue;
        startColumn = 0;
        endColumn = 0;
        for (int i = 0; i < doubleList[0].Length; i++)
        {
            for (int j = 0; j < doubleList.Count - 1; j++)
            {
                count = 0;
                tempStep = doubleList[j + 1][i] - doubleList[j][i];
                for (int k = j + 1; k < doubleList.Count - 1 && Math.Abs(doubleList[k][i] + tempStep - doubleList[k + 1][i]) < 10 && tempStep > 0.00001; k++)
                    count++;
                if (countMax < count)
                {
                    countMax = count + 2;
                    startColumn = j;
                    endColumn = startColumn + count + 1;
                    xRow = i;
                }
                j += count;
            }
        }
        if (countMax < 10 || xRow == -1) return null;//10以下かXを見つけられなかったらNullをかえす

        if (endColumn + 1 < doubleList.Count)
            doubleList.RemoveRange(endColumn + 1, doubleList.Count - endColumn - 1);
        doubleList.RemoveRange(0, startColumn);

        //y軸を決める 標準偏差が一番大きい数が格納されている
        //↑のやり方だと、間違ってしまうことがあるみたい。
        //単にxRowの次(+1)のインデックスにしてみる。もし+1がダメだったら、-1にする。
        int yRow;

        if (xRow + 1 < doubleList[0].Length)
            yRow = xRow + 1;
        else if (xRow - 1 > -1)
            yRow = xRow - 1;
        else
            return null;
        //旧: 標準偏差が最大の列を y 軸に選ぶ方式 (上記理由で xRow±1 に変更)

        //最後に値を代入
        DiffractionProfile2 dif = new();
        for (int i = 0; i < doubleList.Count; i++)
            dif.SourceProfile.Pt.Add(new PointD(doubleList[i][xRow], doubleList[i][yRow]));

        return dif;
    }
}