using System;
using System.IO;

namespace Crystallography;

public static class DiffractionOptics
{
    public static Parameter Read(string filename)
    {
        try
        {
            var fi = new FileInfo(filename);
            var ct = fi.CreationTime;
            var lwt = fi.LastWriteTime;

            //var sr = new StreamReader(filename);
            //var str = sr.ReadToEnd();
            //sr.Close();
            var str = File.ReadAllText(filename); // 260718Cl: 例外時に StreamReader がリークしファイルロックが残るため置換 (既定エンコーディング挙動は同一)

            var replace = new Func<string, string, string, string>((string srcString, string oldValue, string newValue)
                => srcString.Replace("<" + oldValue, "<" + newValue).Replace("</" + oldValue, "</" + newValue));

            str = replace(str, "parameter", "Parameter");

            str = replace(str, "SACLA_EH5>", "FootMode>");

            str = replace(str, "centerX", "DirectSpotX");
            str = replace(str, "centerY", "DirectSpotY");
            str = replace(str, "cameraLength", "CameraLength1");

            str = replace(str, "SACLA_EH5_FootX", "FootX");
            str = replace(str, "SACLA_EH5_FootY", "FootY");
            str = replace(str, "SACLA_EH5_CameraLength2", "CameraLength2");

            //var sw = new StreamWriter(filename);
            //sw.Write(str);
            //sw.Close();
            //File.WriteAllText(filename, str); // 260718Cl: 同上 (StreamWriter リーク対策) // 261105Cl 変更前 (下の 5 行も)
            //
            //fi = new FileInfo(filename)
            //{
            //    CreationTime = ct,
            //    LastWriteTime = lwt
            //};
            // 261105Cl 変更: 旧版の要素名の置き換えはメモリ上だけで行い、元のファイルには書き戻さない (開くだけで .prm が書き換わり、時刻だけ元に戻していた。
            //  PdiFile・結晶リスト XML と同じ直し方)。読むのは置き換えた本文から

            System.Xml.Serialization.XmlSerializer serializer = new(typeof(Parameter));
            //var fs = new FileStream(filename, FileMode.Open);
            //var prm = (Parameter)serializer.Deserialize(fs);
            //fs.Close();
            //using var fs = new FileStream(filename, FileMode.Open); // 260718Cl: Deserialize 例外時 (不正XML) に FileStream がリークし、ユーザーがファイル修正後に再読込できなくなるため using 化 // 261105Cl 変更前
            using var fs = new StringReader(str); // 261105Cl 変更: 置き換えた本文から読む
            var prm = (Parameter)serializer.Deserialize(fs);

            return prm;
        }
        catch
        {
            return null;
        }
    }

    /// <summary>261105Cl 追加: .prm (<see cref="Parameter"/>) の数値文字列を地域設定に依らず解釈する。
    /// <para>IPAnalyzer は 261105Cl から数値を InvariantCulture ("0.1") で書くが、それより前は地域設定 (CurrentCulture) で書いていたので
    /// "0,1" (de, fr など) や U+2212 の負号を使う "−0,1" (sv, nb, fi など) も来る。CurrentCulture の Convert.ToDouble では、
    /// de-DE が "0.1" を桁区切りとみなして黙って 1 と読み、fr-FR / sv-SE は例外になっていた。</para>
    /// <para>.prm の各要素は単一の数値で、桁区切りは一度も書いていない (どれも double/decimal/int の ToString() そのもの)。そこで
    /// "," が 1 個だけなら小数点とみなす。"," が 2 個以上、または "," と "." の混在は、どちらの書き方でも生じないので
    /// 推測せず <see cref="FormatException"/> にする (黙って誤読しない)。null は <see cref="ArgumentNullException"/>。</para></summary>
    public static double ParseNumber(string s) // 261105Cl 追加
    {
        ArgumentNullException.ThrowIfNull(s); // 261105Cl 追加
        s = s.Trim().Replace('\u2212', '-'); // 261105Cl 追加
        if (s.Contains(',')) // 261105Cl 追加
        {
            if (s.Contains('.') || s.IndexOf(',') != s.LastIndexOf(',')) // 261105Cl 追加
                throw new FormatException($"Ambiguous number in the parameter file: {s}"); // 261105Cl 追加
            s = s.Replace(',', '.'); // 261105Cl 追加
        }
        return double.Parse(s, System.Globalization.NumberStyles.Float, System.Globalization.CultureInfo.InvariantCulture); // 261105Cl 追加
    }

    /// <summary>回折光学系のパラメータ. 主にIPAnzlyzerとReciProから利用。</summary>
    [Serializable()]
    public class Parameter
    {
        /// <summary>"FlatPanel" か "Gandolfi"</summary>
        public string cameraMode;

        public string FootMode = "False";// 以前は SACLA_EH5 TrueであればFootモード、FalseであればDirectSpotモード

        public string DirectSpotX;//以前は centerX
        public string DirectSpotY;//以前は centerY
        public string CameraLength1;//以前は cameraLength

        public string FootX;//以前は SACLA_EH5_FootX
        public string FootY;//以前は SACLA_EH5_FootY
        public string CameraLength2;//以前は SACLA_EH5_CameraLength2

        public string SACLA_EH5_PixelWidth;//廃止予定
        public string SACLA_EH5_PixelHeight;//廃止予定
        public string SACLA_EH5_PixleSize;//廃止予定
        public string SACLA_EH5_TwoTheta;//=Tau 廃止予定
        public string SACLA_EH5_Tau;//=Tau 廃止予定
        public string SACLA_EH5_Phi;//=Phi 廃止予定
        public string SACLA_EH5_Distance;//＝CameraLength2 廃止予定
        public string pixSize;//廃止予定
        public string aspectRatio;//廃止予定

        public string waveSource;
        public string electronEnergy;

        public string xRayElement;
        public string xRayLine;
        public string waveLength;

        public string pixSizeX;
        public string pixSizeY;
        public string pixKsi;

        public string tiltPhi;
        public string tiltTau;

        public string sphericalRadiusInverse;
        public string GandolfiRadius;

        /// <summary>Concentric か Radial</summary>
        public string IntegrationMode;

        public string ConcentricStart;
        public string ConcentricEnd;
        public string ConcentricStep;

        /// <summary>Angle, Length, d-spacing</summary>
        public string ConcentricUnit;

        public string RadialRadius;
        public string RadialWidth;
        public string RadialStep;

        /// <summary>Angle, d-spacing</summary>
        public string RadialUnit;

        //Rectangle or Sector
        public string RegionMode;

        public string RectangleDirection;
        public string RectangleBandWidth;
        public string RectangleBothSide;
        public string RectangleAngle;

        public string SectorStart;
        public string SectorEnd;

        public string ExceptMaskedSpots;
        public string ExceptEdges;
        public string ExceptOver;
        public string ExceptUnder;

        public Parameter()
        {
        }
    }
}
