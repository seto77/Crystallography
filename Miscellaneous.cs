using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Drawing;
using System.IO;
using System.Reflection;
using System.Runtime.CompilerServices;
using System.Threading;

namespace Crystallography;

public struct FastSpinLock
{
    private const int SYNC_ENTER = 1;
    private const int SYNC_EXIT = 0;
    private int _syncFlag;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void Enter()
    {
        if (Interlocked.CompareExchange(ref _syncFlag, SYNC_ENTER, SYNC_EXIT) == SYNC_ENTER)
        {
            Spin();
        }
        return;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void Exit() => Volatile.Write(ref _syncFlag, SYNC_EXIT);

    [MethodImpl(MethodImplOptions.NoInlining)]
    private void Spin()
    {
        var spinner = new SpinWait();
        spinner.SpinOnce();
        while (Interlocked.CompareExchange(ref _syncFlag, SYNC_ENTER, SYNC_EXIT) == SYNC_ENTER)
        {
            spinner.SpinOnce();
        }
    }
}

public static class Miscellaneous
{
    public static Color BlendColor(Color col1, Color col2, double ratio)
    {
        //var r = ratio > 1 ? 1 : ratio < 0 ? 0 : ratio;
        var r = Math.Clamp(ratio, 0, 1); // 261111Cl 変更: 手書きのクランプ → Math.Clamp (NaN・-0 も旧と同じ値を返す)
        return Color.FromArgb(
            (byte)(col1.A * r + col2.A * (1 - r)),
            (byte)(col1.R * r + col2.R * (1 - r)),
            (byte)(col1.G * r + col2.G * (1 - r)),
            (byte)(col1.B * r + col2.B * (1 - r)));
    }


    /// <summary>数字に stとかthみたいな文字を追加した文字列で返す</summary>
    /// <param name="number"></param>
    /// <returns></returns>
    public static string Ordinal(this int number)
    {
        var ones = number % 10;
        var tens = Math.Floor(number / 10f) % 10;
        if (tens == 1)
            return number + "th";
        return ones switch
        {
            1 => number + "st",
            2 => number + "nd",
            3 => number + "rd",
            _ => number + "th",
        };
    }

    /// <summary>有限の数字かどうかを判定する</summary>
    /// <param name="d"></param>
    /// <returns></returns>
    //public static bool IsFiniteNumber(params double[] d)
    //{
    //    foreach (var value in d)
    //    {
    //        if (double.IsNaN(value) || double.IsInfinity(value))
    //            return false;
    //    }
    //    return true;
    //}
    // 261111Cl 変更: params double[] → params ReadOnlySpan<double> (可変長の呼び出しで配列を作らない)、判定は double.IsFinite (NaN・±∞ を偽にするのは同じ)。
    // 呼び出しは 4 アプリとも GraphControl の 4 引数の 1 か所だけで、ソースはそのまま通る。null を渡すと旧は例外・新は空 (true) になる (null を渡す呼び出しは無い)
    public static bool IsFiniteNumber(params ReadOnlySpan<double> d)
    {
        foreach (var value in d)
        {
            if (!double.IsFinite(value))
                return false;
        }
        return true;
    }

    private static bool isDecimalPointCommaFlag = true;
    private static bool isDecimalPointComma = false;

    /// <summary>小数点がカンマかどうかを判定する</summary>
    public static bool IsDecimalPointComma
    {
        get
        {
            if (isDecimalPointCommaFlag)
            {
                isDecimalPointComma = (1.01).ToString().Contains(',');
                isDecimalPointCommaFlag = false;
            }
            return isDecimalPointComma;
        }
    }

    //public static (int Division, int Modulus) DivMod(int n, int m) => (n / m, n % m);
    public static (int Division, int Modulus) DivMod(int n, int m) => Math.DivRem(n, m); // 261111Cl 変更: Math.DivRem (同じ値・同じ例外)

    /// <summary>ファイルが使用中かどうかをチェック</summary>
    /// <param name="path"></param>
    /// <returns></returns>
    public static bool isFileExistsAndLocked(string path)
    {
       
        if (File.Exists(path))
        {
            FileStream stream = null;

            try
            {
                stream = new FileStream(path, FileMode.Open, FileAccess.ReadWrite, FileShare.ReadWrite);
            }
            catch (DirectoryNotFoundException) //260317Cl 未使用変数e削除
            {
                return false;
            }
            catch (FileNotFoundException) //260317Cl 未使用変数e削除
            {
                return false;
            }
            catch (IOException) //260317Cl 未使用変数e削除
            {
                if (File.Exists(path))
                {
                    return true;
                }
            }
            catch (Exception) //260317Cl 未使用変数e削除
            {
                return false;
            }
            finally
            {
                stream?.Close();
            }

            return false;
        }

        return false;
    }
}



/// <summary>プロパティグリッドのプロパティの並び順をソート</summary>
public class DefinitionOrderTypeConverter : TypeConverter
{
    public override PropertyDescriptorCollection GetProperties(ITypeDescriptorContext context, object value, Attribute[] attributes)
    {
        // TypeDescriptorを使用してプロパティ一覧を取得する
        var pdc = TypeDescriptor.GetProperties(value, attributes);

        // プロパティ一覧をリフレクションから取得
        var type = value.GetType();
        var list = new List<string>();
        foreach (PropertyInfo propertyInfo in type.GetProperties())
            list.Add(propertyInfo.Name);
        // リフレクションから取得した順でソート
        return pdc.Sort([.. list]);
    }

    /// <summary>GetPropertiesをサポートしていることを表明する。</summary>
    public override bool GetPropertiesSupported(ITypeDescriptorContext context) => true;
}
