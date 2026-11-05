using System;
using System.Numerics;
using System.Xml.Serialization;

namespace Crystallography;

[Serializable()]
public class Atoms : System.IEquatable<Atoms>, ICloneable
{
    #region 基本メソッド
    public object Clone()
    {
        // (260320Ch) MemberwiseClone 後に配列自体を複製してから各要素を複製する
        Atoms atoms = (Atoms)this.MemberwiseClone();
        atoms.Atom = Atom is null ? [] : new Vector3D[Atom.Length];
        for (int i = 0; i < atoms.Atom.Length; i++)
            atoms.Atom[i] = Atom[i] is null ? null : (Vector3D)Atom[i].Clone();
        // 261105Cl 追加: 温度因子と同位体組成も複製する (以前は MemberwiseClone のまま共有していたため、複製した結晶の格子に合わせて
        // Dsf.Cell を設定すると、元の結晶の温度因子の係数まで変わった)
        atoms.Dsf = Dsf?.Clone();
        atoms.Isotope = (double[])Isotope?.Clone();
        return atoms;
    }

    public override string ToString()
    {
        return $"{Label}\t{ElementName}\t{GetStringFromDouble(X)}\t{GetStringFromDouble(Y)}\t{GetStringFromDouble(Z)}\t{GetStringFromDouble(Occ)}\t{Multiplicity}\t{WyckoffLeter}\t{SiteSymmetry}";
    }

    public bool Equals(Atoms obj)
    {
        // (260320Ch) null 安全化して比較条件を明示する
        return ReferenceEquals(this, obj) || obj is not null && obj.Label == Label && obj.X == X && obj.Y == Y && obj.Z == Z && obj.Occ == Occ;
    }

    public override bool Equals(object obj)
    {
        // (260320Ch) as による中間代入を避けて型判定を明示する
        return obj is Atoms atoms && Equals(atoms);
    }

    public override int GetHashCode()
    {
        // (260320Ch) 匿名型依存をやめ、配列要素も含めて明示的にハッシュ化する
        var hash = new HashCode();
        hash.Add(X);
        hash.Add(Y);
        hash.Add(Z);
        hash.Add(Occ);
        hash.Add(Label);
        hash.Add(WyckoffNumber);
        hash.Add(SymmetrySeriesNumber);
        hash.Add(Dsf);
        hash.Add(Texture);

        if (Atom != null)
            foreach (var atom in Atom)
                hash.Add(atom);

        return hash.ToHashCode();
    }
    #endregion

    #region フィールド、プロパティ

    public int ID;

    [XmlIgnore]
    public Vector3D[] Atom = [];

    public double X, Y, Z;
    public double X_err, Y_err, Z_err;


    public double Occ, Occ_err;
    public int AtomicNumber;

    public int SubNumberXray = 0;
    public int SubNumberElectron = 0;
    public double[] Isotope = null;

    public string Label;

    public string ElementName;

    [XmlIgnore]
    public string WyckoffLeter, SiteSymmetry;

    [XmlIgnore]
    public int WyckoffNumber;

    [XmlIgnore]
    public int Multiplicity;

    public DiffuseScatteringFactor Dsf;

    [XmlIgnore]
    public int SymmetrySeriesNumber;

    public int Argb;

    public float Radius = 0.6f;

    public float Ambient = Material.DefaultTexture.Ambient;//環境光
    public float Diffusion = Material.DefaultTexture.Diffuse;//拡散光
    public float Emission = Material.DefaultTexture.Emission;//自己証明
    public float Shininess = Material.DefaultTexture.SpecularPow;//反射光の強度
    public float Specular = Material.DefaultTexture.Specular;//反射光

    [XmlIgnore]
    public Vector3DBase PositionError => new(X_err, Y_err, Z_err);
    [XmlIgnore]
    public Vector3DBase Position => new(X, Y, Z);

    [XmlIgnore]
    public (float Ambient, float Diffusion, float Specular, float Shininess, float Emission) Texture
    {
        get => (Ambient, Diffusion, Specular, Shininess, Emission);
        set
        {
            Ambient = value.Ambient;
            Diffusion = value.Diffusion;
            Specular = value.Specular;
            Shininess = value.Shininess;
            Emission = value.Emission;
        }
    }

    public Material Material => new(Argb, (Ambient, Diffusion, Specular, Shininess, Emission));

    /// <summary>OpenGL描画時に、ラベルを表示するか</summary>
    public bool ShowLabel = false;

    /// <summary>OpenGLの描画時に有効にするかどうか</summary>
    public bool GLEnabled = true;

    //260604Cl Reduce XML size: omit default-valued render members on serialize
    public bool ShouldSerializeArgb() => Argb != 0;
    public bool ShouldSerializeRadius() => Radius != 0.6f;
    public bool ShouldSerializeAmbient() => Ambient != Material.DefaultTexture.Ambient;
    public bool ShouldSerializeDiffusion() => Diffusion != Material.DefaultTexture.Diffuse;
    public bool ShouldSerializeEmission() => Emission != Material.DefaultTexture.Emission;
    public bool ShouldSerializeShininess() => Shininess != Material.DefaultTexture.SpecularPow;
    public bool ShouldSerializeSpecular() => Specular != Material.DefaultTexture.Specular;
    public bool ShouldSerializeShowLabel() => ShowLabel;
    public bool ShouldSerializeGLEnabled() => !GLEnabled;
    public bool ShouldSerializeSubNumberXray() => SubNumberXray != 0;
    public bool ShouldSerializeSubNumberElectron() => SubNumberElectron != 0;
    public bool ShouldSerializeID() => ID != 0;
    public bool ShouldSerializeX_err() => X_err != 0;
    public bool ShouldSerializeY_err() => Y_err != 0;
    public bool ShouldSerializeZ_err() => Z_err != 0;
    public bool ShouldSerializeOcc_err() => Occ_err != 0;
    #endregion

    #region コンストラクタ

    public Atoms()
    {
    }


    /// <summary>ワイコフポジションだけ指定して、原子位置は実態のないコンストラクタ</summary>
    /// <param name="wyk"></param>
    public Atoms(WyckoffPosition wyk, string label, int atomicNumber, int subXray, int subElectron, double[] isotope,
         double occ, DiffuseScatteringFactor dsf)
    {
        SymmetrySeriesNumber = wyk.SymmetrySeriesNumber;

        Label = label;

        this.Occ = 1;

        WyckoffLeter = wyk.WyckoffLetter;
        SiteSymmetry = wyk.SiteSymmetry;
        Multiplicity = wyk.Multiplicity;
        WyckoffNumber = wyk.WyckoffNumber;

        Atom = new Vector3D[Multiplicity];
        for (int i = 0; i < Multiplicity; i++)
            Atom[i] = new Vector3D(0, 0, 0);

        Dsf = dsf;

        SubNumberXray = subXray;
        SubNumberElectron = subElectron;
        AtomicNumber = atomicNumber;
        Isotope = isotope;
        ElementName = AtomicNumber.ToString() + ": " + AtomStatic.AtomicName(atomicNumber);

        Isotope = isotope;
    }

    /// <summary>基本コンストラクタ</summary>
    /// <param name="label"></param>
    /// <param name="atomicNumber"></param>
    /// <param name="subXray"></param>
    /// <param name="subElectron"></param>
    /// <param name="isotope"></param>
    /// <param name="symmetrySeriesNumber"></param>
    /// <param name="pos"></param>
    /// <param name="occ"></param>
    /// <param name="dsf"></param>
    public Atoms(string label, int atomicNumber, int subXray, int subElectron, double[] isotope, int symmetrySeriesNumber,
        Vector3DBase pos, double occ, DiffuseScatteringFactor dsf)
    {
        SymmetrySeriesNumber = symmetrySeriesNumber;

        Label = label;

        Occ = occ;
        X = pos.X;
        Y = pos.Y;
        Z = pos.Z;

        var temp = WyckoffPosition.GetEquivalentAtomsPosition((X, Y, Z), symmetrySeriesNumber);
        WyckoffLeter = temp.WyckoffLeter;
        SiteSymmetry = temp.SiteSymmetry;
        Multiplicity = temp.Multiplicity;
        WyckoffNumber = temp.WyckoffNumber;

        Atom = temp.Atom;
        Dsf = dsf;

        SubNumberXray = subXray;
        SubNumberElectron = subElectron;
        Isotope = isotope ?? [];
        AtomicNumber = atomicNumber;
        ElementName = AtomicNumber.ToString() + ": " + AtomStatic.AtomicName(atomicNumber);
    }

    /// <summary>基本コンストラクタ + エラー</summary>
    /// <param name="label"></param>
    /// <param name="atomicNumber"></param>
    /// <param name="subXray"></param>
    /// <param name="subElectron"></param>
    /// <param name="isotope"></param>
    /// <param name="symmetrySeriesNumber"></param>
    /// <param name="pos"></param>
    /// <param name="pos_err"></param>
    /// <param name="occ"></param>
    /// <param name="occ_err"></param>
    /// <param name="dsf"></param>
    public Atoms(string label, int atomicNumber, int subXray, int subElectron, double[] isotope, int symmetrySeriesNumber,
       Vector3DBase pos, Vector3DBase pos_err, double occ, double occ_err, DiffuseScatteringFactor dsf)
        : this(label, atomicNumber, subXray, subElectron, isotope, symmetrySeriesNumber, pos, occ, dsf)
    {

        X_err = pos_err.X;
        Y_err = pos_err.Y;
        Z_err = pos_err.Z;
        Occ_err = occ_err;
    }

    /// <summary>基本コンストラクタ + エラー + Material</summary>
    /// <param name="label"></param>
    /// <param name="atomicNumber"></param>
    /// <param name="subXray"></param>
    /// <param name="subElectron"></param>
    /// <param name="isotope"></param>
    /// <param name="symmetrySeriesNumber"></param>
    /// <param name="pos"></param>
    /// <param name="pos_err"></param>
    /// <param name="occ"></param>
    /// <param name="occ_err"></param>
    /// <param name="dsf"></param>
    /// <param name="mat"></param>
    /// <param name="radius"></param>
    public Atoms(string label, int atomicNumber, int subXray, int subElectron, double[] isotope, int symmetrySeriesNumber,
       Vector3DBase pos, Vector3DBase pos_err, double occ, double occ_err,
        DiffuseScatteringFactor dsf, Material mat, float radius, bool glEnabled = true, bool showLabel = false)
        : this(label, atomicNumber, subXray, subElectron, isotope, symmetrySeriesNumber, pos, pos_err, occ, occ_err, dsf)
    {
        Radius = radius;
        if (mat != null)
        {
            Argb = mat.Argb;
            Texture = mat.Texture;
        }
        GLEnabled = glEnabled;
        ShowLabel = showLabel;

    }
    #endregion

    #region ResetSymmetry  対称性をリセットして、等価な原子位置や、ワイコフ位置などを再設定する。
    /// <summary>対称性をリセットして、等価な原子位置や、ワイコフ位置などを再設定する。</summary>
    /// <param name="symmetrySeriesNumber"></param>
    public void ResetSymmetry(int symmetrySeriesNumber)
    {
        SymmetrySeriesNumber = symmetrySeriesNumber;

        var temp = WyckoffPosition.GetEquivalentAtomsPosition((X, Y, Z), symmetrySeriesNumber);
        WyckoffLeter = temp.WyckoffLeter;
        SiteSymmetry = temp.SiteSymmetry;
        Multiplicity = temp.Multiplicity;
        WyckoffNumber = temp.WyckoffNumber;
        ElementName = $"{AtomicNumber}: {AtomStatic.AtomicName(AtomicNumber)}";

        Atom = temp.Atom;
    }
    #endregion

    #region ResetVesta()  Vestaの色、材質にセットする。
    /// <summary>Vestaの色、材質にセットする。</summary>
    public void ResetVesta()
    {
        Texture = Material.DefaultTexture;
        (Radius, Argb) =AtomStatic.GetVesta(AtomicNumber);
    }
    #endregion

    #region 原子位置の乱数変化、微小変化      テストコード
    /// <summary>多重度を保ったまま、原子位置を乱数的的に変化させる。ワイコフ位置も多重度をもとにきまる</summary>
    /// <param name="r"></param>
    public static void RandomizeKeepintMultiplicity(Random r)
    {
    }

    /// <summary>ワイコフ位置を保ったまま、原子位置を乱数的に変化させる</summary>
    /// <param name="seed"></param>
    public void RandomizeKeepingWykoff(Random r)
    {
        var wyk = SymmetryStatic.WyckoffPositions[SymmetrySeriesNumber][WyckoffNumber];
        if (wyk.Free.X) X = r.NextDouble();
        if (wyk.Free.Y) Y = r.NextDouble();
        if (wyk.Free.Z) Z = r.NextDouble();

        Atom = SymmetryStatic.WyckoffPositions[SymmetrySeriesNumber][WyckoffNumber].GeneratePositions(X, Y, Z);
        X = Atom[0].X;
        Y = Atom[0].Y;
        Z = Atom[0].Z;
    }

    /// <summary>ワイコフ位置を保ったまま、原子位置を現在の位置からわずかに変化させる</summary>
    /// <param name="threshold">動かす最大値 (相対位置) </param>
    public void ShakeKeepingWykoff(double threshold, Random r)
    {
        var wyk = SymmetryStatic.WyckoffPositions[SymmetrySeriesNumber][WyckoffNumber];
        if (wyk.Free.X) X += (r.NextDouble() - 1) * (r.NextDouble() - 1) * 4 * threshold;
        if (wyk.Free.Y) Y += (r.NextDouble() - 1) * (r.NextDouble() - 1) * 4 * threshold;
        if (wyk.Free.Z) Z += (r.NextDouble() - 1) * (r.NextDouble() - 1) * 4 * threshold;

        Atom = SymmetryStatic.WyckoffPositions[SymmetrySeriesNumber][WyckoffNumber].GeneratePositions(X, Y, Z);
        X = Atom[0].X;
        Y = Atom[0].Y;
        Z = Atom[0].Z;
    }
    #endregion

    #region 電子線、X線、中性子の、原子散乱因子の計算 
    /// <summary>電子線の原子散乱因子を計算 s2の単位はnm^-2</summary>
    /// <param name="S2">S2: (sin(theta)/ramda)^2, unit is nm^-2</param>
    /// <returns></returns>
    public double GetAtomicScatteringFactorForElectron(double s2)
        => AtomStatic.ElectronScatteringPeng[AtomicNumber][SubNumberElectron].Factor(s2) * Occ;

    /// <summary>X線の原子散乱因子を計算 s2の単位はnm^-2</summary>
    /// <param name="s2"> unit is nm^-2</param>
    /// <returns></returns>
    public double GetAtomicScatteringFactorForXray(double s2)
        => AtomStatic.XrayScatteringWK[AtomicNumber][SubNumberXray].Factor(s2) * Occ;

    public Complex GetAtomicScatteringFactorForNeutron()
    {
        if (Isotope != null && Isotope.Length == AtomStatic.IsotopeAbundance[AtomicNumber].Count)
        {
            var f = new Complex();
            for (int i = 0; i < AtomStatic.IsotopeAbundance[AtomicNumber].Count; i++)
                f += AtomStatic.NeutronCoherentScattering[AtomicNumber][i + 1] * Isotope[i] / 100.0;
            return f * Occ;
        }
        else
            return AtomStatic.NeutronCoherentScattering[AtomicNumber][0] * Occ;
    }
    #endregion

    #region 静的メソッド
    public static string GetStringFromDouble(double d)
    {
        if (Math.Abs(d - 0.125) < 0.000000001) return "1/8";
        else if (Math.Abs(d - 0.375) < 0.000000001) return "3/8";
        else if (Math.Abs(d - 0.625) < 0.000000001) return "5/8";
        else if (Math.Abs(d - 0.875) < 0.000000001) return "7/8";
        else if (Math.Abs(d - 0.25) < 0.000000001) return "1/4";
        else if (Math.Abs(d - 0.75) < 0.000000001) return "3/4";
        else if (Math.Abs(d - 0.5) < 0.000000001) return "1/2";
        else if (Math.Abs(d - 1.0 / 3.0) < 0.000000001) return "1/3";
        else if (Math.Abs(d - 2.0 / 3.0) < 0.000000001) return "2/3";
        else if (Math.Abs(d - 1.0 / 6.0) < 0.000000001) return "1/6";
        else if (Math.Abs(d - 5.0 / 6.0) < 0.000000001) return "5/6";
        else if (Math.Abs(d - 1.0 / 12.0) < 0.000000001) return "1/12";
        else if (Math.Abs(d - 5.0 / 12.0) < 0.000000001) return "5/12";
        else if (Math.Abs(d - 7.0 / 12.0) < 0.000000001) return "7/12";
        else if (Math.Abs(d - 11.0 / 12.0) < 0.000000001) return "11/12";
        else if (Math.Abs(d - 1.0 / 24.0) < 0.000000001) return "1/24";
        else if (Math.Abs(d - 5.0 / 24.0) < 0.000000001) return "5/24";
        else if (Math.Abs(d - 7.0 / 24.0) < 0.000000001) return "7/24";
        else if (Math.Abs(d - 11.0 / 24.0) < 0.000000001) return "11/24";
        else if (Math.Abs(d - 13.0 / 24.0) < 0.000000001) return "13/24";
        else if (Math.Abs(d - 17.0 / 24.0) < 0.000000001) return "17/24";
        else if (Math.Abs(d - 19.0 / 24.0) < 0.000000001) return "19/24";
        else if (Math.Abs(d - 23.0 / 24.0) < 0.000000001) return "23/24";
        else return d.ToString("g6");
    }

    public static double GetDoubleFromString(string s)
    {
        if (s == "1/8" || s == "1.0/8.0") return 1.0 / 8.0;
        else if (s == "3/8" || s == "3.0/8.0") return 3.0 / 8.0;
        else if (s == "5/8" || s == "5.0/8.0") return 5.0 / 8.0;
        else if (s == "7/8" || s == "7.0/8.0") return 7.0 / 8.0;
        else if (s == "1/4" || s == "1.0/4.0") return 1.0 / 4.0;
        else if (s == "3/4" || s == "3.0/4.0") return 3.0 / 4.0;
        else if (s == "1/2" || s == "1.0/2.0") return 1.0 / 2.0;
        else if (s == "1/3" || s == "1.0/3.0") return 1.0 / 3.0;
        else if (s == "2/3" || s == "2.0/3.0") return 2.0 / 3.0;
        else if (s == "1/6" || s == "1.0/6.0") return 1.0 / 6.0;
        else if (s == "5/6" || s == "5.0/6.0") return 5.0 / 6.0;
        else if (s == "1/12" || s == "1.0/12.0") return 1.0 / 12.0;
        else if (s == "5/12" || s == "5.0/12.0") return 5.0 / 12.0;
        else if (s == "7/12" || s == "7.0/12.0") return 7.0 / 12.0;
        else if (s == "11/12" || s == "11.0/12.0") return 11.0 / 12.0;
        else if (s == "1/12" || s == "1.0/12.0") return 1.0 / 12.0;
        else if (s == "5/12" || s == "5.0/12.0") return 5.0 / 12.0;
        else if (s == "7/12" || s == "7.0/12.0") return 7.0 / 12.0;
        else if (s == "11/12" || s == "11.0/12.0") return 11.0 / 12.0;
        else if (s == "13/12" || s == "13.0/12.0") return 13.0 / 12.0;
        else if (s == "17/12" || s == "17.0/12.0") return 17.0 / 12.0;
        else if (s == "19/12" || s == "19.0/12.0") return 19.0 / 12.0;
        else if (s == "23/12" || s == "23.0/12.0") return 23.0 / 12.0;
        else if (s == "0.3333") return 1.0 / 3.0;
        else if (s == "0.6667") return 2.0 / 3.0;
        else
            try { return Convert.ToDouble(s); }
            catch { System.Windows.Forms.MessageBox.Show("Please input a valid value"); return 0; }
    }
    #endregion
}

[Serializable()]
public class DiffuseScatteringFactor
{
    #region プロパティ、フィールド
    const double PI2 = Math.PI * Math.PI;
    public enum Type { U, B }
    //Biomolecular Crystallography: Principles, Practice, and Application to Structural Biology
    //641ページ

    #region B type. Getのみ
    /// <summary>unit: nm^2</summary>
    public double Biso => OriginalType == Type.B ? Iso : Iso * PI2 * 8;
    /// <summary>260919Cl 追加: g=0 の実効等方 B [nm²]。IsZero なら 0、等方なら Biso、非等方なら (Biso が NaN のとき) Biso000。getU の g=0 経路と同じ選択規則を 1 か所にまとめたもの。NaN は 0 にする。</summary>
    //public double BisoEffective { get { double b = IsZero ? 0 : UseIso ? Biso : double.IsNaN(Biso) ? Biso000 : Biso; return double.IsNaN(b) || b < 0 ? 0 : b; } } // 261105Cl 変更前
    // 261105Cl 変更: 非等方なら常にテンソルから求めた等価な等方 B (Biso000、Hamilton の式) を使う。以前は Iso が NaN のときだけ Biso000 に落ちたため、
    // (1) GUI で異方性だけ入力した Iso=0 のサイトは g=0 の B が 0 になり、吸収の U'0 が 0 になっていた (v.4.950 の持ち越し)。
    // (2) 異方性を編集しても、隠れた Iso (以前の値や CIF の U_eq) が残って使われ、g≠0 のテンソルと食い違った。
    // CIF の U_eq とは U_ij の丸めの範囲で一致する (tools/AdpCheck realcif)
    public double BisoEffective { get { double b = IsZero ? 0 : UseIso ? Biso : Biso000; return double.IsNaN(b) || b < 0 ? 0 : b; } }
    /// <summary>260919Cl 追加: B 未設定 (0) の原子に与える零点振動相当の下限 [nm²] (= 0.1 Å²)。後方散乱源 σ_n と MC のコヒーレンス破壊判定で使う (getU の吸収ポテンシャルには適用しない)。</summary>
    public const double BisoFloorNm2 = 1E-3;

    /// <summary>unit: nm^2. g=000の時のBiso. Acta Cryst. (1959). 12, 609 , Hamilton の式に従って、Bisoを計算</summary>
    //public double Biso000 => (B11 * a2 + B22 * b2 + B33 * c2 + 2 * B12 * ab + 2 * B23 * bc + 2 * B31 * ca) * 4.0 / 3.0; // 261105Cl 変更前
    // 261105Cl 変更: Hamilton (1959) 式 (18) B = (4/3) Σ_ij β_ij (a_i·a_j) は実格子の内積 (cos を含む) を使う。以前は |a||b| などの積だけで、
    // 直交でない格子 (六方・単斜・三斜) の B_eq が誤っていた (六方晶で 16 % 程度。tools/AdpCheck beq)
    public double Biso000 => (B11 * a2 + B22 * b2 + B33 * c2 + 2 * B12 * gab + 2 * B23 * gbc + 2 * B31 * gca) * 4.0 / 3.0;

    /// <summary>温度因子がゼロの場合はtrue</summary>
    public bool IsZero => UseIso ? Biso == 0 : B11 == 0 && B22 == 0 && B33 == 0 && B12 == 0 && B23 == 0 && B31 == 0;

    /// <summary>unit: nm^2</summary>
    public double Biso_err => OriginalType == Type.B ? Iso_err : Iso_err * PI2 * 8;
    /// <summary>unit: none</summary>
    public double B11 => OriginalType == Type.B ? Aniso11 : Aniso11 * coeff11;
    /// <summary>unit: none</summary>
    public double B22 => OriginalType == Type.B ? Aniso22 : Aniso22 * coeff22;
    /// <summary>unit: none</summary>
    public double B33 => OriginalType == Type.B ? Aniso33 : Aniso33 * coeff33;
    /// <summary>unit: none</summary>
    public double B12 => OriginalType == Type.B ? Aniso12 : Aniso12 * coeff12;
    /// <summary>unit: none</summary>
    public double B23 => OriginalType == Type.B ? Aniso23 : Aniso23 * coeff23;
    /// <summary>unit: none</summary>
    public double B31 => OriginalType == Type.B ? Aniso31 : Aniso31 * coeff31;
    /// <summary>unit: none</summary>
    public double B11_err => OriginalType == Type.B ? Aniso11_err : Aniso11_err * coeff11;
    /// <summary>unit: none</summary>
    public double B22_err => OriginalType == Type.B ? Aniso22_err : Aniso22_err * coeff22;
    /// <summary>unit: none</summary>
    public double B33_err => OriginalType == Type.B ? Aniso33_err : Aniso33_err * coeff33;
    /// <summary>unit: none</summary>
    public double B12_err => OriginalType == Type.B ? Aniso12_err : Aniso12_err * coeff12;
    /// <summary>unit: none</summary>
    public double B23_err => OriginalType == Type.B ? Aniso23_err : Aniso23_err * coeff23;
    /// <summary>unit: none</summary>
    public double B31_err => OriginalType == Type.B ? Aniso31_err : Aniso31_err * coeff31;

    #endregion

    #region U type. Getのみ
    /// <summary>単位は nm^2</summary>
    public double Uiso => OriginalType == Type.U ? Iso : Iso / PI2 / 8;
    /// <summary>単位は nm^2</summary>
    public double Uiso_err => OriginalType == Type.U ? Iso_err : Iso_err / PI2 / 8;
    /// <summary>単位は nm^2</summary>
    public double U11 => OriginalType == Type.U ? Aniso11 : Aniso11 / coeff11;
    /// <summary>単位は nm^2</summary>
    public double U22 => OriginalType == Type.U ? Aniso22 : Aniso22 / coeff22;
    /// <summary>単位は nm^2</summary>
    public double U33 => OriginalType == Type.U ? Aniso33 : Aniso33 / coeff33;
    /// <summary>単位は nm^2</summary>
    public double U12 => OriginalType == Type.U ? Aniso12 : Aniso12 / coeff12;
    /// <summary>単位は nm^2</summary>
    public double U23 => OriginalType == Type.U ? Aniso23 : Aniso23 / coeff23;
    /// <summary>単位は nm^2</summary>
    public double U31 => OriginalType == Type.U ? Aniso31 : Aniso31 / coeff31;
    /// <summary>単位は nm^2</summary>
    public double U11_err => OriginalType == Type.U ? Aniso11_err : Aniso11_err / coeff11;
    /// <summary>単位は nm^2</summary>
    public double U22_err => OriginalType == Type.U ? Aniso22_err : Aniso22_err / coeff22;
    /// <summary>単位は nm^2</summary>
    public double U33_err => OriginalType == Type.U ? Aniso33_err : Aniso33_err / coeff33;
    /// <summary>単位は nm^2</summary>
    public double U12_err => OriginalType == Type.U ? Aniso12_err : Aniso12_err / coeff12;
    /// <summary>単位は nm^2</summary>
    public double U23_err => OriginalType == Type.U ? Aniso23_err : Aniso23_err / coeff23;
    /// <summary>単位は nm^2</summary>
    public double U31_err => OriginalType == Type.U ? Aniso31_err : Aniso31_err / coeff31;
    #endregion

    #region オリジナルの値 
    /// <summary>単位は nm^2</summary>
    public double Iso { get; set; }
    /// <summary>単位は nm^2</summary>
    public double Iso_err { get; set; }
    /// <summary>単位は Uの場合 nm^2,  Bの場合　無次元</summary>
    public double Aniso11 { get; set; }
    /// <summary>単位は Uの場合 nm^2,  Bの場合　無次元</summary>
    public double Aniso22 { get; set; }
    /// <summary>単位は Uの場合 nm^2,  Bの場合　無次元</summary>
    public double Aniso33 { get; set; }
    /// <summary>単位は Uの場合 nm^2,  Bの場合　無次元</summary>
    public double Aniso12 { get; set; }
    /// <summary>単位は Uの場合 nm^2,  Bの場合　無次元</summary>
    public double Aniso23 { get; set; }
    /// <summary>単位は Uの場合 nm^2,  Bの場合　無次元</summary>
    public double Aniso31 { get; set; }

    /// <summary>単位は Uの場合 nm^2,  Bの場合　無次元</summary>
    public double Aniso11_err { get; set; }
    /// <summary>単位は Uの場合 nm^2,  Bの場合　無次元</summary>
    public double Aniso22_err { get; set; }
    /// <summary>単位は Uの場合 nm^2,  Bの場合　無次元</summary>
    public double Aniso33_err { get; set; }
    /// <summary>単位は Uの場合 nm^2,  Bの場合　無次元</summary>
    public double Aniso12_err { get; set; }
    /// <summary>単位は Uの場合 nm^2,  Bの場合　無次元</summary>
    public double Aniso23_err { get; set; }
    /// <summary>単位は Uの場合 nm^2,  Bの場合　無次元</summary>
    public double Aniso31_err { get; set; }

    #endregion

    public bool UseIso { get; set; }
    public Type OriginalType { get; set; } = Type.B;

    //260604Cl Reduce XML size: omit default-valued Dsf members on serialize
    public bool ShouldSerializeIso() => Iso != 0;
    public bool ShouldSerializeIso_err() => Iso_err != 0;
    public bool ShouldSerializeAniso11() => Aniso11 != 0;
    public bool ShouldSerializeAniso22() => Aniso22 != 0;
    public bool ShouldSerializeAniso33() => Aniso33 != 0;
    public bool ShouldSerializeAniso12() => Aniso12 != 0;
    public bool ShouldSerializeAniso23() => Aniso23 != 0;
    public bool ShouldSerializeAniso31() => Aniso31 != 0;
    public bool ShouldSerializeAniso11_err() => Aniso11_err != 0;
    public bool ShouldSerializeAniso22_err() => Aniso22_err != 0;
    public bool ShouldSerializeAniso33_err() => Aniso33_err != 0;
    public bool ShouldSerializeAniso12_err() => Aniso12_err != 0;
    public bool ShouldSerializeAniso23_err() => Aniso23_err != 0;
    public bool ShouldSerializeAniso31_err() => Aniso31_err != 0;
    public bool ShouldSerializeUseIso() => UseIso;
    public bool ShouldSerializeOriginalType() => OriginalType != Type.B;

    [XmlIgnore]
    public (double A, double B, double C, double Alpha, double Beta, double Gamma) Cell
    {
        get => cell;
        set
        {
            cell = value;
            a2 = cell.A * cell.A;
            b2 = cell.B * cell.B;
            c2 = cell.C * cell.C;
            ab = cell.A * cell.B;
            bc = cell.B * cell.C;
            ca = cell.C * cell.A;
            var (sinAlpha, cosAlpha) = Math.SinCos(cell.Alpha);
            var (sinBeta,cosBeta) = Math.SinCos(cell.Beta);
            var (sinGamma, cosGamma) = Math.SinCos(cell.Gamma);
            var v = cell.A * cell.B * cell.C * Math.Sqrt(1 - cosAlpha * cosAlpha - cosBeta * cosBeta - cosGamma * cosGamma + 2 * cosAlpha * cosBeta * cosGamma);
            var aStar = bc * sinAlpha / v;
            var bStar = ca * sinBeta / v;
            var cStar = ab * sinGamma / v;
            var cosAlphaStar = (cosBeta * cosGamma - cosAlpha) / sinBeta / sinGamma;
            var cosBetaStar = (cosGamma * cosAlpha - cosBeta) / sinGamma / sinAlpha;
            var cosGammaStar = (cosAlpha * cosBeta - cosGamma) / sinAlpha / sinBeta;
            coeff11 = PI2 * 2 * aStar * aStar;
            coeff22 = PI2 * 2 * bStar * bStar;
            coeff33 = PI2 * 2 * cStar * cStar;
            //coeff12 = PI2 * 2 * aStar * bStar * cosGammaStar; // 261105Cl 変更前 (3 行)
            //coeff23 = PI2 * 2 * bStar * cStar * cosAlphaStar;
            //coeff31 = PI2 * 2 * cStar * aStar * cosBetaStar;
            // 261105Cl 変更: CIF の U_ij は β_ij = 2π² a*_i a*_j U_ij で、非対角にも逆格子の角の cos は付かない
            // (CIF core 辞書 _atom_site_aniso_U_ij、Grosse-Kunstleve & Adams (2002) J. Appl. Cryst. 35, 477 式 4-6)。
            // 以前は cos γ* などを掛けていたため、直交の格子では U12・U23・U31 が消え、六方晶では U12 が半分になっていた
            coeff12 = PI2 * 2 * aStar * bStar;
            coeff23 = PI2 * 2 * bStar * cStar;
            coeff31 = PI2 * 2 * cStar * aStar;
            // 261105Cl 追加: Biso000 (Hamilton の式) で使う実格子の内積 a_i·a_j
            gab = ab * cosGamma;
            gbc = bc * cosAlpha;
            gca = ca * cosBeta;
        }
    }
    private (double A, double B, double C, double Alpha, double Beta, double Gamma) cell = (0, 0, 0, 0, 0, 0);

    //private double a2, b2, c2, ab, bc, ca, coeff11, coeff22, coeff33, coeff12, coeff23, coeff31; // 261105Cl 変更前
    private double a2, b2, c2, ab, bc, ca, gab, gbc, gca, coeff11, coeff22, coeff33, coeff12, coeff23, coeff31; // 261105Cl 変更: gab, gbc, gca (実格子の内積) を追加
    #endregion

    #region コンストラクタ
    public DiffuseScatteringFactor()
    {
    }

    /// <summary>コンストラクタ. B##は無次元, 他はnm^2. Cellの 単位は nm & radians.</summary>
    /// <param name="t"></param>
    /// <param name="useIso"></param>
    /// <param name="iso"></param>
    /// <param name="iso_err"></param>
    /// <param name="aniso"></param>
    /// <param name="aniso_err"></param>
    /// <param name="cell"></param>
    public DiffuseScatteringFactor(Type t, bool useIso, double iso, double iso_err, double[] aniso, double[] aniso_err,
        (double A, double B, double C, double Alpha, double Beta, double Gamma) cell)
    {
        OriginalType = t;
        UseIso = useIso;
        Iso = iso;
        Iso_err = iso_err;
        Aniso11 = aniso[0];
        Aniso22 = aniso[1];
        Aniso33 = aniso[2];
        Aniso12 = aniso[3];
        Aniso23 = aniso[4];
        Aniso31 = aniso[5];
        Aniso11_err = aniso_err[0];
        Aniso22_err = aniso_err[1];
        Aniso33_err = aniso_err[2];
        Aniso12_err = aniso_err[3];
        Aniso23_err = aniso_err[4];
        Aniso31_err = aniso_err[5];

        Cell = cell;
    }
    #endregion

    /// <summary>261105Cl 追加: 値と格子 (係数) を含めた複製</summary>
    public DiffuseScatteringFactor Clone()
        => new(OriginalType, UseIso, Iso, Iso_err,
            [Aniso11, Aniso22, Aniso33, Aniso12, Aniso23, Aniso31],
            [Aniso11_err, Aniso22_err, Aniso33_err, Aniso12_err, Aniso23_err, Aniso31_err], Cell);

    #region 座標変換 261105Cl 追加
    /// <summary>261105Cl 追加: 分率座標の線形変換 x' = M x (+ 並進) に合わせて非等方テンソルを変換した、新しい DiffuseScatteringFactor を返す。
    /// β' = M β Mᵀ (Grosse-Kunstleve &amp; Adams (2002) J. Appl. Cryst. 35, 477 式 10 の U* = R U* Rᵀ と同じ。β = 2π² U*)。等方の値 (Iso) は変わらない。
    /// 等価位置への展開 (M = 対称操作の回転部分 R)、超構造 (M = diag(1/u, 1/v, 1/w))、軸の取り替え (M = 新しい分率座標 ← 古い分率座標) に使う。
    /// 元の型を保つ: B 型の異方性 (= 無次元の β) はそのまま、U 型は新しい格子の a* で U'_ij = β'_ij / (2π² a*'_i a*'_j)。
    /// 誤差は 6 成分を独立とみなして伝播する (共分散は持っていないので近似)。</summary>
    /// <param name="m">分率座標の変換行列。x'_i = Σ_j m[i, j] x_j</param>
    /// <param name="oldCell">変換前の格子 (nm, rad)。U から β を求めるのに使う (Cell が未設定のオブジェクトでも正しく変換するため明示する)</param>
    /// <param name="newCell">変換後の格子 (nm, rad)</param>
    public DiffuseScatteringFactor Transform(double[,] m, (double A, double B, double C, double Alpha, double Beta, double Gamma) oldCell,
        (double A, double B, double C, double Alpha, double Beta, double Gamma) newCell)
    {
        var src = new DiffuseScatteringFactor(OriginalType, UseIso, Iso, Iso_err,
            [Aniso11, Aniso22, Aniso33, Aniso12, Aniso23, Aniso31],
            [Aniso11_err, Aniso22_err, Aniso33_err, Aniso12_err, Aniso23_err, Aniso31_err], oldCell);
        double[] p = [src.B11, src.B22, src.B33, src.B12, src.B23, src.B31];
        double[] pe = [src.B11_err, src.B22_err, src.B33_err, src.B12_err, src.B23_err, src.B31_err];
        (int I, int J)[] index = [(0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (2, 0)]; // 独立な 6 成分の (行, 列)。並びは Aniso と同じ
        var beta = new double[6];
        var betaErr = new double[6];
        for (int n = 0; n < 6; n++)
        {
            var (i, j) = index[n];
            double v = 0, e2 = 0;
            bool errUnknown = false;
            for (int q = 0; q < 6; q++)
            {
                var (k, l) = index[q];
                // β'_ij = Σ_kl m_ik β_kl m_jl。非対角の β_kl (= β_lk) は (k,l) と (l,k) の 2 通りで現れる
                var c = k == l ? m[i, k] * m[j, k] : m[i, k] * m[j, l] + m[i, l] * m[j, k];
                if (c == 0)
                    continue; // 効かない成分は値も誤差も見ない (0 × NaN で結果を不明にしない)
                v += c * p[q];
                if (double.IsNaN(pe[q])) errUnknown = true; // 効く成分の誤差が不明なら、結果の誤差も不明
                else e2 += c * c * pe[q] * pe[q];
            }
            beta[n] = v;
            betaErr[n] = errUnknown ? double.NaN : Math.Sqrt(e2);
        }
        var dst = new DiffuseScatteringFactor(OriginalType, UseIso, Iso, Iso_err, new double[6], new double[6], newCell);
        double[] coeff = OriginalType == Type.B ? [1, 1, 1, 1, 1, 1] : [dst.coeff11, dst.coeff22, dst.coeff33, dst.coeff12, dst.coeff23, dst.coeff31];
        (dst.Aniso11, dst.Aniso22, dst.Aniso33, dst.Aniso12, dst.Aniso23, dst.Aniso31)
            = (beta[0] / coeff[0], beta[1] / coeff[1], beta[2] / coeff[2], beta[3] / coeff[3], beta[4] / coeff[4], beta[5] / coeff[5]);
        (dst.Aniso11_err, dst.Aniso22_err, dst.Aniso33_err, dst.Aniso12_err, dst.Aniso23_err, dst.Aniso31_err)
            = (betaErr[0] / coeff[0], betaErr[1] / coeff[1], betaErr[2] / coeff[2], betaErr[3] / coeff[3], betaErr[4] / coeff[4], betaErr[5] / coeff[5]);
        return dst;
    }

    /// <summary>261105Cl 追加: 格子が変形したとき (歪み・圧縮・状態方程式)、原子の揺れの大きさ (Cartesian の変位の分散 U_cart, Å²) を保つ温度因子を返す
    /// (作者の方針 2026-11-05: B 型・U 型のどちらで持っていても同じ結果にする)。原子は分率座標のまま格子と一緒に動き、揺れの楕円体は伸び縮みしない、とする。
    /// 変形前の実格子の行列を A (列 = a, b, c)、変形を F = R S (R は回転、S = √(FᵀF) は純粋な伸縮) とすると、楕円体は R で回るだけなので
    /// U* = A⁻¹ U_cart A⁻ᵀ (Grosse-Kunstleve &amp; Adams 2002) から β' = M β Mᵀ、M = A⁻¹ S⁻¹ A。FᵀF = A⁻ᵀ G' A⁻¹ (G' は変形後の計量) なので、
    /// M は変形前後の格子定数だけで決まる (A の向きの取り方にも、変形後の格子の向きにも依らない)。等方の値 (Iso) は変えない。
    /// 直交格子の軸方向の伸縮では M = diag(a/a', b/b', c/c') となり、U 型の U_ij を固定するのと同じ。B 型 (β) を固定すると揺れが格子と一緒に伸びる。</summary>
    /// <param name="oldCell">変形前の格子 (nm, rad)</param>
    /// <param name="newCell">変形後の格子 (nm, rad)</param>
    public DiffuseScatteringFactor Deform((double A, double B, double C, double Alpha, double Beta, double Gamma) oldCell,
        (double A, double B, double C, double Alpha, double Beta, double Gamma) newCell)
    {
        var (a, b, c, al, be, ga) = oldCell;
        // 変形前の実格子 (a を x 軸、b を xy 面に置く。M は向きの取り方に依らない)
        var cy = (Math.Cos(al) - Math.Cos(be) * Math.Cos(ga)) / Math.Sin(ga);
        var axes = MathNet.Numerics.LinearAlgebra.Double.DenseMatrix.OfArray(new double[,] {
            { a, b * Math.Cos(ga), c * Math.Cos(be) },
            { 0, b * Math.Sin(ga), c * cy },
            { 0, 0, c * Math.Sqrt(1 - Math.Cos(be) * Math.Cos(be) - cy * cy) } });
        var (a2, b2, c2, al2, be2, ga2) = newCell;
        var g2 = MathNet.Numerics.LinearAlgebra.Double.DenseMatrix.OfArray(new double[,] {
            { a2 * a2, a2 * b2 * Math.Cos(ga2), a2 * c2 * Math.Cos(be2) },
            { a2 * b2 * Math.Cos(ga2), b2 * b2, b2 * c2 * Math.Cos(al2) },
            { a2 * c2 * Math.Cos(be2), b2 * c2 * Math.Cos(al2), c2 * c2 } });
        var inv = axes.Inverse();
        var ftf = inv.Transpose() * g2 * inv; // = FᵀF (対称・正定値)
        var evd = ftf.Evd(MathNet.Numerics.LinearAlgebra.Symmetricity.Symmetric);
        var d = new double[3];
        for (int i = 0; i < 3; i++)
            d[i] = 1 / Math.Sqrt(evd.EigenValues[i].Real);
        var sInv = evd.EigenVectors * MathNet.Numerics.LinearAlgebra.Double.DenseMatrix.OfDiagonalArray(d) * evd.EigenVectors.Transpose(); // S⁻¹ = (FᵀF)^(-1/2)
        return Transform((inv * sInv * axes).ToArray(), oldCell, newCell);
    }
    #endregion
}
