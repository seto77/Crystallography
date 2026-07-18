using System;
using System.Collections.Generic;

namespace Crystallography;

public class NumericalFormula
{
    public static double GetNumetricValue(string[] str)
    {
        try
        {
            //'=' で定義された定数を計算し、それ以降の式へ展開する
            for (int i = 0; i < str.Length; i++)
                if (str[i].Contains('=')) // '=' の文字列がみつかったら
                {
                    // 260718Cl: Split を1回化。
                    var parts = str[i].Split(["="], StringSplitOptions.RemoveEmptyEntries);
                    string leftString = parts[0].Replace(" ", "");
                    string rightString = parts[1].Replace(" ", "");

                    for (int j = i + 1; j < str.Length; j++)
                        // 260718Cl: 旧実装は k=IndexOf(leftString) を走査位置 l を無視して常に先頭から取り直すため、識別子内に埋め込まれた先頭出現が棄却されると
                        //           以降が永久に置換されなかった。走査位置 l から検索し、置換時は挿入長ぶん、非置換時は k+1 へ進める。
                        for (int l = 0; l < str[j].Length;)
                        {
                            int k = str[j].IndexOf(leftString, l, StringComparison.Ordinal);
                            if (k < 0) break;

                            // 直前・直後が英字なら大きな識別子の一部とみなし置換しない。
                            // 260718Cl: 旧 'A'<=c && 'z'>=c は 'Z'(90) と 'a'(97) の間の [ \ ] ^ _ ` (91-96) を英字と誤判定し、x^2 等で置換に失敗して NaN になっていた → char.IsAsciiLetter に修正。
                            bool embedded =
                                (k + leftString.Length < str[j].Length && char.IsAsciiLetter(str[j][k + leftString.Length])) ||
                                (k > 0 && char.IsAsciiLetter(str[j][k - 1]));

                            if (embedded)
                                l = k + 1;
                            else
                            {
                                str[j] = str[j].Remove(k, leftString.Length).Insert(k, "(" + rightString + ")");
                                l = k + rightString.Length + 2; // 挿入した "(rightString)" の直後へ
                            }
                        }
                }
            return NumericalValue(str[^1]);
        }
        catch
        {
            return double.NaN;
        }
    }

    private static double NumericalValue(string str)
    {
        str = str.Replace(" ", "");

        //"e" があって直前直後が数値のときは "*10^"に変換しておく
        for (int i = 0; i < str.Length - 1; i++)
            if (str.Substring(i, 1).ToLower() == "e" && (str[i + 1] == '-' || str[i + 1] == '+' || (str[i + 1] >= '0' && str[i + 1] <= '9')))
            {
                str = str.Remove(i, 1);
                if (i == 0)
                    str = str.Insert(i, "10^");
                else
                    str = str.Insert(i, "*10^");
            }

        var list = new List<object>();

        for (int i = 0; i < str.Length; i++)
        {
            if (i + 1 < str.Length)//2文字の関数、定数を検索
            {
                string func = str.Substring(i, 2).ToLower();
                if (func == "pi" || func == "ln")
                {
                    list.Add(func);
                    str = str.Remove(i, 2);
                    i = 0;
                }
            }
            if (i + 2 < str.Length)//3文字の関数、定数を検索
            {
                string func = str.Substring(i, 3).ToLower();
                if (func == "sin" || func == "cos" || func == "tan" || func == "exp" || func == "log" || func == "abs")
                {
                    list.Add(func);
                    str = str.Remove(i, 3);
                    i = 0;
                }
            }
            if (i + 3 < str.Length)//4文字の関数、定数を検索
            {
                string func = str.Substring(i, 4).ToLower();
                if (func == "asin" || func == "acos" || func == "atan" || func == "sqrt")
                {
                    list.Add(func);
                    str = str.Remove(i, 4);
                    i = 0;
                }
            }

            if (i < str.Length && str[i] == '(')  //かっこの始まりが現れたら
            {
                int count = 1;
                for (int j = 1; j < str.Length; j++)//対応するかっこの終りをみつける
                {
                    if (str[j] == '(')
                        count++;
                    else if (str[j] == ')')
                        count--;

                    if (count == 0)//見つかったら
                    {
                        list.Add(NumericalValue(str[1..j]));
                        str = str.Remove(0, j + 1);
                        i = 0; //次を0に戻す
                        break;
                    }
                }
            }

            if (i < str.Length && (str[i] == '+' || str[i] == '-' || str[i] == '*' || str[i] == '/' || str[i] == '^')) //演算子が現れたら
            {
                if (i != 0)
                    list.Add(Convert.ToDouble(str[..i]));//演算子の直前までを数値に変換し格納

                list.Add(str[i]);//演算子を格納

                str = str.Remove(0, i + 1);
                i = -1;//次を0に戻す
            }
        }
        //最後にstrに余りがあればそれを数値に変換
        if (str.Length > 0)
            list.Add(Convert.ToDouble(str));

        if (list.Count == 0)
            return 0;

        //最初に関数をチェック
        for (int i = 0; i < list.Count; i++)
            if (list[i].GetType() == typeof(string))
            {
                if ((string)list[i] == "pi")//定数の場合
                    list[i] = Math.PI;
                else if (i + 1 < list.Count && list[i + 1].GetType() == typeof(double))
                {
                    switch (list[i])
                    {
                        case "ln": list[i] = Math.Log((double)list[i + 1]); break;
                        case "sin": list[i] = Math.Sin((double)list[i + 1] / 180 * Math.PI); break;
                        case "cos": list[i] = Math.Cos((double)list[i + 1] / 180 * Math.PI); break;
                        case "tan": list[i] = Math.Tan((double)list[i + 1] / 180 * Math.PI); break;
                        case "exp": list[i] = Math.Exp((double)list[i + 1]); break;
                        case "log": list[i] = Math.Log10((double)list[i + 1]); break;
                        case "abs": list[i] = Math.Abs((double)list[i + 1]); break;
                        case "asin": list[i] = Math.Asin((double)list[i + 1]) / Math.PI * 180; break;
                        case "acos": list[i] = Math.Acos((double)list[i + 1]) / Math.PI * 180; break;
                        case "atan": list[i] = Math.Atan((double)list[i + 1]) / Math.PI * 180; break;
                        case "sqrt": list[i] = Math.Sqrt((double)list[i + 1]); break;
                    }

                    list.RemoveAt(i + 1);
                }
                else
                    return double.NaN;
            }

        //先頭が'-'あるいは'+'で始まる場合を考慮する
        if (list.Count > 1)
            if (list[0].GetType() == typeof(char) && ((char)list[0] == '-' || (char)list[0] == '+'))
            {
                if (list[1].GetType() == typeof(double))
                {
                    if ((char)list[0] == '-')
                        list[1] = -(double)list[1];
                    list.RemoveAt(0);
                }
                else
                    return double.NaN;//先頭が'-'なのに次が数値でないときはNaNを返す
            }

        //次に 演算子の後に"-"あるいは'+'が来る場合に対処
        if (list.Count > 2)
            for (int i = 0; i < list.Count - 2; i++)
            {
                if (list[i].GetType() == typeof(char) && list[i + 1].GetType() == typeof(char))//2個演算子が続いているところ
                {
                    if (list[i + 2].GetType() != typeof(double))
                        return double.NaN;//2個続いたあとのやつが数値でないときはNaNを返す

                    if ((char)list[i + 1] == '-')// '-'だったら符号を変える
                        list[i + 2] = -(double)list[i + 2];
                    else if ((char)list[i + 1] != '+')// '-'でも'+'でもないときはNaNをかえす
                        return double.NaN;
                    list.RemoveAt(i + 1);
                }
            }

        //ここまでで数値はかならず演算子一つに挟まれるようになっているはずなのでチェック
        //listのカウントが奇数個で、かつ偶数番目が数値、奇数が演算子になるはず
        if (list.Count % 2 != 1)
            return double.NaN;

        // 260712Cl 修正: 条件が && だと i%2==0 と i%2==1 が同時成立せず恒偽で、型検証が死んでいた。|| が正しく、
        // 不正なトークン列 (偶数位置が非double / 奇数位置が非char) を下流のキャスト例外でなく NaN で穏当に弾く。正常入力は不変。
        // if ((i % 2 == 0 && list[i].GetType() != typeof(double)) && (i % 2 == 1 && list[i].GetType() != typeof(char))) // 260712Cl 変更前
        for (int i = 0; i < list.Count; i++)
            if ((i % 2 == 0 && list[i].GetType() != typeof(double)) || (i % 2 == 1 && list[i].GetType() != typeof(char)))
                return double.NaN;

        //まず^を後ろから探す
        for (int i = list.Count - 1; i >= 0; i--)
            if (list[i].GetType() == typeof(char) && (char)list[i] == '^')
            {
                double v = Math.Pow((double)list[i - 1], (double)list[i + 1]);
                list.RemoveRange(i, 2);
                list[i - 1] = v;
            }

        //次に*と/を探す
        for (int i = 1; i < list.Count; i++)
        {
            if (list[i].GetType() == typeof(char) && (char)list[i] == '*')
            {
                double v = (double)list[i - 1] * (double)list[i + 1];
                list.RemoveRange(i, 2);
                list[i - 1] = v;
                i = 0;
            }
            else if (list[i].GetType() == typeof(char) && (char)list[i] == '/')
            {
                double v = (double)list[i - 1] / (double)list[i + 1];
                list.RemoveRange(i, 2);
                list[i - 1] = v;
                i = 0;
            }
        }

        //最後に+と-を探す
        for (int i = 1; i < list.Count; i++)
        {
            if (list[i].GetType() == typeof(char) && (char)list[i] == '+')
            {
                double v = (double)list[i - 1] + (double)list[i + 1];
                list.RemoveRange(i, 2);
                list[i - 1] = v;
                i = 0;
            }
            else if (list[i].GetType() == typeof(char) && (char)list[i] == '-')
            {
                double v = (double)list[i - 1] - (double)list[i + 1];
                list.RemoveRange(i, 2);
                list[i - 1] = v;
                i = 0;
            }
        }

        if (list.Count == 1)
            return (double)list[0];
        else
            return double.NaN;
    }
}
