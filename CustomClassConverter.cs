using System;
using System.ComponentModel;

namespace Crystallography
{
    /// <summary>変換を行う</summary>
    public class ConverterBase : ExpandableObjectConverter
    {
        public override bool CanConvertFrom(ITypeDescriptorContext context, Type sourceType)
        {
            return sourceType == typeof(string) || base.CanConvertFrom(context, sourceType);// stringからなら変換可能
        }

        public override bool CanConvertTo(ITypeDescriptorContext context, Type destinationType)
        {
            return destinationType == typeof(string) || base.CanConvertTo(context, destinationType);// stringへなら変換可能
        }

        /// <summary>261111Cl 追加: 値を、引数の culture (null なら現在のカルチャ) の区切り (TextInfo.ListSeparator + 空白) と数値の書式で並べる
        /// (System.Drawing.PointConverter と同じ方式)。旧来は culture を無視して現在のカルチャの数値を ", " で並べていた
        /// (日本語・英語の環境では区切りが ',' なので結果は同じ)</summary>
        protected static string JoinValues(System.Globalization.CultureInfo culture, params double[] values)
        {
            culture ??= System.Globalization.CultureInfo.CurrentCulture;
            var texts = new string[values.Length];
            for (int i = 0; i < values.Length; i++)
                texts[i] = values[i].ToString(culture);
            return string.Join(culture.TextInfo.ListSeparator + " ", texts);
        }

        /// <summary>261111Cl 追加: JoinValues の逆。culture (null なら現在のカルチャ) の区切りの最初の文字で最大 maxParts 個に切り、先頭の count 個を
        /// その culture で読む (旧来は ',' で切って現在のカルチャで読んでいたので、小数点がカンマの環境では読めなかった)</summary>
        protected static double[] SplitValues(string text, System.Globalization.CultureInfo culture, int count, int maxParts)
        {
            culture ??= System.Globalization.CultureInfo.CurrentCulture;
            var ss = text.Split(culture.TextInfo.ListSeparator[0], maxParts);
            var values = new double[count];
            for (int i = 0; i < count; i++)
                values[i] = double.Parse(ss[i], culture);
            return values;
        }
    }

    public class RectangleDConverter : ConverterBase
    {
        public override PropertyDescriptorCollection GetProperties(ITypeDescriptorContext context, object value, Attribute[] filter)
        {
            return TypeDescriptor.GetProperties(value, filter).Sort(new[] { "X", "Y", "Width", "Height" });
        }

        public override object ConvertTo(ITypeDescriptorContext context, System.Globalization.CultureInfo culture, object value, Type destinationType)
        {
            if (destinationType == typeof(string) && value is RectangleD d)
            {
                //return $"{d.X}, {d.Y}, {d.Width}, {d.Height}"; //260317Cl string.Format → 文字列補間 // 261111Cl 変更前
                return JoinValues(culture, d.X, d.Y, d.Width, d.Height); // 261111Cl 変更: 引数の culture の区切りと書式で
            }
            return base.ConvertTo(context, culture, value, destinationType);
        }

        public override object ConvertFrom(ITypeDescriptorContext context, System.Globalization.CultureInfo culture, object value)
        {
            if (value is string @string)
            {
                //string[] ss = @string.Split(',', 4); //260317Cl new char[] { ',' } → char リテラル // 261111Cl 変更前 (3 行)
                ////260317Cl 変更: Convert.ToDouble → double.Parse
                //return new RectangleD(double.Parse(ss[0]), double.Parse(ss[1]), double.Parse(ss[2]), double.Parse(ss[3]));
                var v = SplitValues(@string, culture, 4, 4); // 261111Cl 変更: 引数の culture の区切りと書式で読む
                return new RectangleD(v[0], v[1], v[2], v[3]);
            }
            return base.ConvertFrom(context, culture, value);
        }
    }

    public class PointDConverter : ConverterBase
    {
        public override PropertyDescriptorCollection GetProperties(ITypeDescriptorContext context, object value, Attribute[] filter)
        {
            return TypeDescriptor.GetProperties(value, filter).Sort(new[] { "X", "Y" });
        }

        public override object ConvertTo(ITypeDescriptorContext context, System.Globalization.CultureInfo culture, object value, Type destinationType)
        {
            if (destinationType == typeof(string) && value is PointD g)
                //return $"{g.X}, {g.Y}"; //260317Cl string.Format → 文字列補間 // 261111Cl 変更前
                return JoinValues(culture, g.X, g.Y); // 261111Cl 変更: 引数の culture の区切りと書式で
            return base.ConvertTo(context, culture, value, destinationType);
        }

        public override object ConvertFrom(ITypeDescriptorContext context, System.Globalization.CultureInfo culture, object value)
        {
            if (value is string s) //260317Cl パターンマッチング + char リテラル
            {
                //string[] ss = s.Split(',', 2); // 261111Cl 変更前 (2 行)
                //return new PointD(double.Parse(ss[0]), double.Parse(ss[1]));
                var v = SplitValues(s, culture, 2, 2); // 261111Cl 変更: 引数の culture の区切りと書式で読む
                return new PointD(v[0], v[1]);
            }
            return base.ConvertFrom(context, culture, value);
        }
    }

    public class SizeDConverter : ConverterBase
    {
        public override PropertyDescriptorCollection GetProperties(ITypeDescriptorContext context, object value, Attribute[] filter)
        {
            return TypeDescriptor.GetProperties(value, filter).Sort(new[] { "Width", "Height" });
        }

        public override object ConvertTo(ITypeDescriptorContext context, System.Globalization.CultureInfo culture, object value, Type destinationType)
        {
            if (destinationType == typeof(string) && value is SizeD d)
                //return $"{d.Width}, {d.Height}"; //260317Cl string.Format → 文字列補間 // 261111Cl 変更前
                return JoinValues(culture, d.Width, d.Height); // 261111Cl 変更: 引数の culture の区切りと書式で
            return base.ConvertTo(context, culture, value, destinationType);
        }

        public override object ConvertFrom(ITypeDescriptorContext context, System.Globalization.CultureInfo culture, object value)
        {
            if (value is string s) //260317Cl パターンマッチング + char リテラル
            {
                //string[] ss = s.Split(',', 2); // 261111Cl 変更前 (2 行)
                //return new SizeD(double.Parse(ss[0]), double.Parse(ss[1]));
                var v = SplitValues(s, culture, 2, 2); // 261111Cl 変更: 引数の culture の区切りと書式で読む
                return new SizeD(v[0], v[1]);
            }
            return base.ConvertFrom(context, culture, value);
        }
    }

    public class Vector3DConverter : ConverterBase
    {
        public override PropertyDescriptorCollection GetProperties(ITypeDescriptorContext context, object value, Attribute[] filter)
        {
            return TypeDescriptor.GetProperties(value, filter).Sort(new[] { "X", "Y", "Z" });
        }

        public override object ConvertTo(ITypeDescriptorContext context, System.Globalization.CultureInfo culture, object value, Type destinationType)
        {
            if (destinationType != typeof(string) || value is not Vector3DBase g) //260317Cl パターンマッチング
                return base.ConvertTo(context, culture, value, destinationType);
            // 260517Cl Vector3D で Text が設定されていればそれを優先（CheckedListBox 等で FormattingEnabled=true のとき ToString ではなくこちらが呼ばれるため）
            if (value is Vector3D v && !string.IsNullOrEmpty(v.Text))
                return v.Text;
            // ConvertFrom が "x,y,z" を Split(',') で読むため、ここは敢えて括弧なしで返す（Vector3DBase.ToString の "(x, y, z)" とは別書式）
            //return $"{g.X}, {g.Y}, {g.Z}"; //260317Cl string.Format → 文字列補間 // 261111Cl 変更前
            return JoinValues(culture, g.X, g.Y, g.Z); // 261111Cl 変更: 引数の culture の区切りと書式で
        }

        public override object ConvertFrom(ITypeDescriptorContext context, System.Globalization.CultureInfo culture, object value)
        {
            if (value is string s) //260317Cl パターンマッチング + char リテラル
            {
                //string[] ss = s.Split(',', 4); // 261111Cl 変更前 (2 行)
                //return new Vector3DBase(double.Parse(ss[0]), double.Parse(ss[1]), double.Parse(ss[2]));
                var v = SplitValues(s, culture, 3, 4); // 261111Cl 変更: 引数の culture の区切りと書式で読む (旧と同じく 4 つまでに切って先頭の 3 つ)
                return new Vector3DBase(v[0], v[1], v[2]);
            }
            return base.ConvertFrom(context, culture, value);
        }
    }
}