using System;
using System.Collections.Generic;
using System.Globalization; // 261026Cl 追加: 小数の添字は文化圏によらず "." で書く
using System.Linq;

namespace Crystallography;

/// <summary>
/// 261025Cl 追加: 元素ごとの組成を、基準陰イオン B (O / ハロゲン / S) による化合物の形式的換算 (酸化物 wt% など) に直す。
/// <para>設定した価数による<b>形式的な換算</b>であり、実在する化学種や相の量を推定するものではない (作者決定 261025)。
/// 例えば NH4Cl を N=−3, H=+1 で換算すると HCl 272.6%, N 26.2%, −Cl≡N −198.8%, 合計 100% になる。</para>
/// <para>計算の定義 (B の価数を −k、全元素の質量の和を M とする):
/// 正の価数 v の元素 E は成分 E_aB_b (質量 n_E (A_E + v A_B / k))、B 以外の負の価数の元素 X は元素の行と
/// 「−B≡X」の補正 (n_X |v_X| A_B / k)、価数 0 の元素は質量だけの行。B の過不足 ΔN_B = n_B − (Σ⁺nv − Σ⁻n|v|)/k は参考行で、
/// 補正後合計 = 100 − 100 A_B ΔN_B / M、電荷合計 = Σ n v (B を含む全元素) = −k ΔN_B。</para>
/// </summary>
public static class CompositionConverter
{
    /// <summary>表示モード (261026Cl 各値の説明を追加)</summary>
    public enum Mode
    {
        /// <summary>元素ごとの重量%・原子%・原子数 (換算しない)</summary>
        Elements,
        /// <summary>O を基準陰イオンとする酸化物への換算</summary>
        Oxides,
        /// <summary>原子番号が最小のハロゲンを基準陰イオンとするハロゲン化物への換算 (O を含まない結晶)</summary>
        Halides,
        /// <summary>S を基準陰イオンとする硫化物への換算 (O を含まない結晶)</summary>
        Sulfides,
    }

    /// <summary>表の行の種類</summary>
    public enum RowKind
    {
        /// <summary>元素モードの行</summary>
        Element,
        /// <summary>正の価数の元素を基準陰イオンとの化合物にした成分 (Fe2O3 など)</summary>
        Component,
        /// <summary>価数 0 の元素 (質量だけ数える)</summary>
        MassOnly,
        /// <summary>価数が未設定の元素 (質量だけ数え、電荷合計は未確定になる)。
        /// 基準陰イオンの価数が負でなく換算できない時は、全元素をこの種類で出す (261026Cl 説明を追加)</summary>
        Unset,
        /// <summary>基準以外の負の価数の元素 (F, Cl など。元素のまま出す)</summary>
        Anion,
        /// <summary>「−B≡X」の補正</summary>
        Correction,
        /// <summary>補正後の合計</summary>
        Total,
        /// <summary>基準陰イオンの過不足 (参考。合計には加えない)</summary>
        Residual,
    }

    /// <summary>表の 1 行。Atoms は規格化の倍率を掛けた原子数 (補正・過不足の行は基準陰イオンの当量)。</summary>
    /// <param name="Kind">行の種類</param>
    /// <param name="Z">元素の原子番号 (補正・過不足の行は基準陰イオン、合計の行は 0)</param>
    /// <param name="Valence">適用した価数 (元素モード・補正・合計・過不足の行は NaN)</param>
    /// <param name="Formula">成分の組成式 (Component の行だけ。他は null)</param>
    /// <param name="WeightPercent">結晶の質量に対する重量% (再規格化しない)</param>
    /// <param name="AtomPercent">原子% (元素モードだけ。換算モードは NaN)</param>
    /// <param name="Atoms">規格化した原子数 (倍率が NaN なら NaN)</param>
    public sealed record Row(RowKind Kind, int Z, double Valence, string Formula, double WeightPercent, double AtomPercent, double Atoms);

    /// <summary>換算の結果</summary>
    public sealed class Result
    {
        /// <summary>実際に使ったモード (指定したモードが使えない結晶では Elements に戻る)。261026Cl 説明を追加</summary>
        public Mode Mode { get; init; }
        /// <summary>基準陰イオンの原子番号 (元素モードは 0)</summary>
        public int BaseZ { get; init; }
        /// <summary>補正に入った陰イオンの原子番号 (表示順)</summary>
        public int[] CorrectionAnions { get; init; } = [];
        /// <summary>表の行 (並び順は Convert の説明を参照)。261026Cl 説明を追加</summary>
        public List<Row> Rows { get; init; } = [];
        /// <summary>基準陰イオンの過不足 (単位胞あたり。正なら過剰)</summary>
        public double ResidualAtomsPerCell { get; init; }
        /// <summary>電荷合計 Σ n v (単位胞あたり。B を含む全元素)。原子数 &gt; 0 で価数が未設定の元素があれば NaN。元素モードと、質量が 0 以下の時も NaN。
        /// 261026Cl 変更: 規格化の倍率を掛けない (全体の電荷は単位胞あたり、作者決定)。フォームはこれをそのまま表示する</summary>
        public double ChargeSum { get; init; }
        /// <summary>原子数 &gt; 0 で価数が未設定の元素があるか。基準陰イオンの価数が負でなく換算できない時も true (261026Cl 説明を追加)</summary>
        public bool HasUnset { get; init; }
    }

    /// <summary>ハロゲン (F, Cl, Br, I)。At は作者の指定 (F, Cl, Br, I) に含めない</summary>
    public static readonly int[] Halogens = [9, 17, 35, 53];

    private static bool IsAnionFamily(int z) => z == 16 || Array.IndexOf(Halogens, z) >= 0;

    /// <summary>そのモードが使えるか (= ラジオを表示するか)。原子数は占有率込み (Crystal.ElementNum)。
    /// <para>酸化物: O &gt; 0 かつ O 以外の元素 &gt; 0。ハロゲン化物: O が無く、F/Cl/Br/I のどれか &gt; 0、かつ陰イオン族 (F, Cl, Br, I, S) 以外の元素 &gt; 0。
    /// 硫化物: O が無く、S &gt; 0、かつ陰イオン族以外の元素 &gt; 0。O がわずかでもあれば酸化物だけ (閾値は設けない、作者決定 261025)。</para></summary>
    public static bool IsAvailable(Mode mode, IReadOnlyList<(int Z, double N)> elements)
    {
        bool has(int z) => elements.Any(e => e.Z == z && e.N > 0);
        var hasO = has(8);
        return mode switch
        {
            Mode.Elements => elements.Count > 0,
            Mode.Oxides => hasO && elements.Any(e => e.Z != 8 && e.N > 0),
            Mode.Halides => !hasO && Halogens.Any(has) && elements.Any(e => !IsAnionFamily(e.Z) && e.N > 0),
            Mode.Sulfides => !hasO && has(16) && elements.Any(e => !IsAnionFamily(e.Z) && e.N > 0),
            _ => false,
        };
    }

    /// <summary>基準陰イオンの原子番号。ハロゲン化物は原子数 &gt; 0 のうち原子番号が最小のハロゲン (作者決定 261025)。使えなければ 0。</summary>
    public static int BaseZ(Mode mode, IReadOnlyList<(int Z, double N)> elements) => mode switch
    {
        Mode.Oxides => 8,
        Mode.Halides => Halogens.FirstOrDefault(z => elements.Any(e => e.Z == z && e.N > 0)),
        Mode.Sulfides => 16,
        _ => 0,
    };

    // 261026Cl 削除: 固定する価数は無くなった (O・S・ハロゲンも手入力で変えられる。作者指示 261025) ので FixedValence 自体を廃止。
    // O を含まない 2 モードの S −2・ハロゲン −1 は DefaultValence の既定として持つ。
    ///// <summary>そのモードで固定される価数 (固定でなければ NaN)。O を含まない 2 モードでは、
    ///// どちらが基準でもハロゲン −1・S −2 に固定する (作者決定 261025)。
    ///// 261025Cl 変更: 酸化物モードの O は固定しない (作者指示。既定 −2 のまま手入力で変えられる)。</summary>
    //public static double FixedValence(Mode mode, int z) => mode switch // 261025Cl 変更前
    //{
    //    //Mode.Oxides => z == 8 ? -2 : double.NaN, // 261025Cl 変更前
    //    Mode.Halides or Mode.Sulfides => z == 16 ? -2 : Array.IndexOf(Halogens, z) >= 0 ? -1 : double.NaN,
    //    _ => double.NaN,
    //};
    //public static double FixedValence(Mode mode, int z) => double.NaN; // 261025Cl (S・ハロゲンの固定もやめた時点で常に NaN だった)

    /// <summary>モード別の既定価数 (鉱物でよく見る酸化数。決まらなければ NaN = 未設定)。
    /// 既存の ChemicalFormulaInputControl の表 (Ag=4 など) とは別に持つ (作者決定 261025)。</summary>
    public static double DefaultValence(Mode mode, int z)
    {
        if (mode is Mode.Halides or Mode.Sulfides)
            switch (z) // O を含まないモードでだけ違う既定
            {
                case 6: return double.NaN; // C (炭化物など)。+4 を引き継がない
                case 7: case 15: return -3; // N, P (窒化物・リン化物)
                case 33: case 51: case 83: return 3; // As, Sb, Bi (硫塩鉱物の陽イオン)
                case 34: case 52: return -2; // Se, Te
                case 16: return -2; // S (硫化物)。261025Cl 追加: 固定をやめたので既定として持つ (O を含むモードの +6 を引き継がない)
            }
        return z switch
        {
            1 => 1,                                              // H → H2O
            2 or 10 or 18 or 36 or 54 or 86 => 0,                // 希ガス (質量だけ)
            3 or 11 or 19 or 37 or 55 or 87 => 1,                // アルカリ金属
            4 or 12 or 20 or 38 or 56 or 88 => 2,                // アルカリ土類金属
            5 => 3, 6 => 4, 7 => 5, 8 => -2,                     // B2O3, CO2, N2O5, O
            9 or 17 or 35 or 53 or 85 => -1,                     // ハロゲン
            13 => 3, 14 => 4, 15 => 5, 16 => 6,                  // Al2O3, SiO2, P2O5, SO3
            21 => 3, 22 => 4, 23 => 3, 24 => 3, 25 => 2,         // Sc2O3, TiO2, V2O3, Cr2O3, MnO
            26 => 2, 27 => 2, 28 => 2, 29 => 2, 30 => 2,         // FeO (作者決定: Fe=2), CoO, NiO, CuO, ZnO
            31 => 3, 32 => 4, 33 => 5, 34 => 4,                  // Ga2O3, GeO2, As2O5, SeO2
            39 => 3, 40 => 4, 41 => 5, 42 => 6, 43 => 7,         // Y2O3, ZrO2, Nb2O5, MoO3, Tc2O7
            44 => 4, 45 => 3, 46 => 2, 47 => 1, 48 => 2,         // RuO2, Rh2O3, PdO, Ag2O, CdO
            49 => 3, 50 => 4, 51 => 5, 52 => 4,                  // In2O3, SnO2, Sb2O5, TeO2
            >= 57 and <= 71 => 3,                                // ランタノイド (Ce も 3)
            72 => 4, 73 => 5, 74 => 6, 75 => 7, 76 => 4,         // HfO2, Ta2O5, WO3, Re2O7, OsO2
            77 => 4, 78 => 4, 79 => 3, 80 => 2,                  // IrO2, PtO2, Au2O3, HgO
            81 => 1, 82 => 2, 83 => 3, 84 => 4,                  // Tl2O, PbO, Bi2O3, PoO2
            89 => 3, 90 => 4, 91 => 5, 92 => 4, 93 => 5, 94 => 4, // Ac2O3, ThO2, Pa2O5, UO2, Np2O5, PuO2
            >= 95 and <= 103 => 3,                               // Am 以降
            _ => double.NaN,
        };
    }

    /// <summary>その元素の全サイト (原子数 &gt; 0) でイオンの X 線散乱因子が選ばれ、価数が一致していればその価数 (無ければ NaN)。</summary>
    public static double IonicValence(Crystal crystal, int z)
    {
        double v = double.NaN;
        foreach (var a in crystal?.Atoms ?? [])
        {
            if (a.AtomicNumber != z || !(a.Multiplicity * a.Occ > 0))
                continue;
            var table = AtomStatic.XrayScatteringWK;
            if (z < 0 || z >= table.Length || a.SubNumberXray < 0 || a.SubNumberXray >= table[z].Length)
                return double.NaN;
            var sv = table[z][a.SubNumberXray].Valence;
            if (sv == 0 || (!double.IsNaN(v) && v != sv))
                return double.NaN;
            v = sv;
        }
        return v;
    }

    /// <summary>適用する価数を決める。優先順位は 保存した手入力 &gt; Fe なら 2 &gt; 全サイトで一致するイオン散乱因子 &gt; モード別の既定 &gt; 未設定 (NaN)。
    /// (261026Cl 変更: 先頭にあった「固定」は廃止。硫化物モードの S の自動の釣り合わせは ResolveValences が行う)</summary>
    public static double ResolveValence(Mode mode, int z, double saved, double ionic)
    {
        //var fixedValence = FixedValence(mode, z); // 261026Cl 削除 (固定の価数は無い)
        //if (!double.IsNaN(fixedValence)) return fixedValence;
        if (!double.IsNaN(saved)) return saved;
        if (z == 26) return 2; // 作者決定 261025: Fe の既定は 2 に固定 (散乱因子で変えない)
        if (!double.IsNaN(ionic)) return ionic;
        return DefaultValence(mode, z);
    }

    /// <summary>全元素の価数を 1 度にまとめて決める (表示 1 回につき 1 回呼ぶ)。各元素は ResolveValence の優先順位で決め、
    /// 硫化物モードの S は、手入力が無ければ他の元素の価数から電荷が釣り合う値にする (FeS2 なら Fe 2 → S −1 で FeS2 100%。作者指示 261025)。
    /// 釣り合う値が求まらない (他に未設定の元素がある) か負にならない時は既定値のまま。
    /// saved は手入力の価数 (無ければ NaN)、ionic はイオン散乱因子の価数 (要る時だけ呼ぶ)。
    /// 261026Cl 追加: フォームにあった S の釣り合わせと、元素ごとに何度も価数を求め直していた処理をここへまとめた</summary>
    public static Dictionary<int, double> ResolveValences(Mode mode, IReadOnlyList<(int Z, double N)> elements, Func<int, double> saved, Func<int, double> ionic)
    {
        var valences = new Dictionary<int, double>();
        foreach (var (z, _) in elements)
        {
            if (valences.ContainsKey(z))
                continue;
            var s = saved(z);
            valences[z] = ResolveValence(mode, z, s, double.IsNaN(s) && z != 26 ? ionic(z) : double.NaN); // 手入力か Fe なら散乱因子は見ない
        }
        if (mode == Mode.Sulfides && valences.ContainsKey(16) && double.IsNaN(saved(16)))
        {
            var balanced = BalancingValence(elements, valences, 16);
            if (balanced < 0)
                valences[16] = balanced;
        }
        return valences;
    }

    /// <summary>他の元素の価数はそのままで、電荷合計 Σ 原子数 × 価数 が 0 になる元素 z の価数 (= −Σ_他 n v / n_z)。
    /// 他の元素 (原子数 ≠ 0) に未設定の価数がある時や、z の原子数が 0 以下の時は NaN。
    /// 261026Cl 追加: フォームの「Σ=0」ボタンと、硫化物モードの S の釣り合わせで使う (フォームから移した)</summary>
    public static double BalancingValence(IReadOnlyList<(int Z, double N)> elements, IReadOnlyDictionary<int, double> valences, int z)
    {
        var n = elements.Where(e => e.Z == z).Sum(e => e.N);
        if (!(n > 0))
            return double.NaN;
        double others = 0, magnitude = 0;
        foreach (var (zz, nn) in elements)
        {
            if (zz == z || nn == 0)
                continue;
            var v = valences.TryGetValue(zz, out var vv) ? vv : double.NaN;
            if (double.IsNaN(v))
                return double.NaN;
            others += nn * v;
            magnitude += Math.Abs(nn * v);
        }
        var balanced = -others / n;
        // 0 になるはずの値に浮動小数の残り (占有率 0.2 × 3 など、1e-16 程度) が乗ると、S の自動の釣り合わせで価数 ≈ −1e-16 となり、
        // 異常な成分名・重量% になるので 0 にそろえる (−0 も 0 にして「-0」と見せない)
        return Math.Abs(balanced) <= 1e-9 * Math.Max(1, magnitude / n) ? 0 : balanced;
    }

    /// <summary>
    /// 正の価数 v (分数可) の元素 e と基準陰イオン b (価数 −k) の成分の組成式 (ASCII)。v = p/q (q ≦ 12) で
    /// E_aB_b, a = kq/g, b = p/g, g = gcd(p, kq) (例: Fe 2 → FeO, 3 → Fe2O3, 8/3 → Fe3O4, Na Cl → NaCl, Ca F → CaF2)。
    /// 小さい整数比にならない (q &gt; 12、または添字 &gt; 24) ときは null (261026Cl 以降、Convert が「MnO1.107」のような小数の添字の式にする)。
    /// 共有の Molecule.CombineCationAndAnion は小数で壊れた式を返すので使わない。
    /// <para>261025Cl 変更: O の価数も手入力で変えられるようにしたので、k も分数 r/s (s ≦ 12) を許す。
    /// 電荷の釣り合い a·v = b·k から a = rq/g, b = sp/g, g = gcd(rq, sp) (k が整数なら従来どおり)。</para>
    /// </summary>
    //public static string ComponentFormula(string e, double v, string b, int k) // 261025Cl 変更前 (k は整数だった)
    public static string ComponentFormula(string e, double v, string b, double k)
    {
        //if (!(v > 0) || double.IsInfinity(v) || k <= 0) // 261025Cl 変更前
        //    return null;
        //static long gcd(long x, long y) { while (y != 0) (x, y) = (y, x % y); return Math.Abs(x); }
        //for (int q = 1; q <= 12; q++)
        //{
        //    var p = Math.Round(v * q);
        //    if (p < 1 || p > 1e6 || Math.Abs(v - p / q) > 1e-9 * Math.Max(1, v))
        //        continue;
        //    long kq = (long)k * q, g = gcd((long)p, kq), na = kq / g, nb = (long)p / g;
        //    if (na > 24 || nb > 24)
        //        return null;
        //    return $"{e}{(na > 1 ? na.ToString() : "")}{b}{(nb > 1 ? nb.ToString() : "")}";
        //}
        //return null;
        if (!(v > 0) || double.IsInfinity(v) || !(k > 0) || double.IsInfinity(k))
            return null;
        //static long gcd(long x, long y) { while (y != 0) (x, y) = (y, x % y); return Math.Abs(x); } // 261026Cl 変更前 (既存の GammaFunction.Gcd を使う)
        // x を分母 12 以下の分数 (分子, 分母) に直す。直せなければ null
        static (long N, long D)? fraction(double x)
        {
            for (int d = 1; d <= 12; d++)
            {
                var n = Math.Round(x * d);
                if (n >= 1 && n <= 1e6 && Math.Abs(x - n / d) <= 1e-9 * Math.Max(1, x))
                    return ((long)n, d);
            }
            return null;
        }
        if (fraction(v) is not var (p, q) || fraction(k) is not var (r, s))
            return null;
        //long rq = r * q, sp = s * p, g = gcd(rq, sp), na = rq / g, nb = sp / g; // 261026Cl 変更前
        long rq = r * q, sp = s * p, g = GammaFunction.Gcd(rq, sp), na = rq / g, nb = sp / g;
        if (na > 24 || nb > 24)
            return null;
        return $"{e}{(na > 1 ? na.ToString() : "")}{b}{(nb > 1 ? nb.ToString() : "")}";
    }

    /// <summary>
    /// 換算する。elements は表示順の (原子番号, 単位胞あたりの原子数)。valence は適用する価数 (NaN = 未設定)。scale は規格化の倍率 (原子数列に掛ける。NaN なら原子数列は NaN。
    /// 261026Cl 変更: 電荷合計には掛けない = 単位胞あたり)。
    /// <para>行の順: 元素モードは元素の順。換算モードは 正の価数の成分 → 価数 0・未設定の元素 → 基準以外の陰イオン → 補正 → 補正後合計 → 過不足 (参考)。</para>
    /// </summary>
    public static Result Convert(Mode mode, IReadOnlyList<(int Z, double N)> elements, Func<int, double> valence, double scale)
    {
        double massTotal = elements.Sum(e => e.N * AtomStatic.AtomicWeight(e.Z));
        double countTotal = elements.Sum(e => e.N);
        if (!(massTotal > 0))
            return new Result { Mode = mode, ChargeSum = double.NaN };

        if (mode == Mode.Elements || !IsAvailable(mode, elements))
            return new Result
            {
                Mode = Mode.Elements,
                ChargeSum = double.NaN,
                Rows = [.. elements.Select(e => new Row(RowKind.Element, e.Z, double.NaN, null,
                    100 * e.N * AtomStatic.AtomicWeight(e.Z) / massTotal, countTotal > 0 ? 100 * e.N / countTotal : double.NaN, e.N * scale))],
            };

        //int baseZ = BaseZ(mode, elements), k = (int)-FixedValence(mode, baseZ); // 261025Cl 変更前
        // 261025Cl 変更: O の価数も手入力で変えられるようにしたので、基準陰イオンの価数 (−k) も valence から取る (分数も可)
        int baseZ = BaseZ(mode, elements);
        double k = -valence(baseZ);
        if (!(k > 0) || double.IsInfinity(k))
        {
            // 基準陰イオンの価数が負でない (未設定・0・正) と換算できない。元素ごとの質量だけ「未設定」の行で出す
            List<Row> unsetRows = [.. elements.Select(e => new Row(RowKind.Unset, e.Z, valence(e.Z), null, 100 * e.N * AtomStatic.AtomicWeight(e.Z) / massTotal, double.NaN, e.N * scale))];
            unsetRows.Add(new Row(RowKind.Total, 0, double.NaN, null, unsetRows.Sum(r => r.WeightPercent), double.NaN, double.NaN));
            var anyUnset = elements.Any(e => e.N > 0 && double.IsNaN(valence(e.Z)));
            return new Result
            {
                Mode = mode,
                BaseZ = baseZ,
                Rows = unsetRows,
                HasUnset = true,
                //ChargeSum = anyUnset ? double.NaN : elements.Where(e => e.N != 0).Sum(e => e.N * valence(e.Z)) * scale, // 261026Cl 変更前
                ChargeSum = anyUnset ? double.NaN : elements.Where(e => e.N != 0).Sum(e => e.N * valence(e.Z)), // 261026Cl 変更: 単位胞あたり (倍率を掛けない)
            };
        }
        double weightB = AtomStatic.AtomicWeight(baseZ);
        double nB = elements.Where(e => e.Z == baseZ).Sum(e => e.N);
        double sumPositive = 0, sumNegative = 0;
        bool hasUnset = false;
        var components = new List<Row>();
        var massOnly = new List<Row>();
        var anions = new List<Row>();
        var correctionZ = new List<int>(); // 261025Cl 追加: 補正 (−B≡X) に入れる陰イオン (原子数 ≠ 0)
        foreach (var (z, n) in elements)
        {
            if (z == baseZ)
                continue;
            var v = valence(z);
            var massPercent = 100 * n * AtomStatic.AtomicWeight(z) / massTotal;
            if (double.IsNaN(v))
            {
                hasUnset |= n > 0;
                massOnly.Add(new Row(RowKind.Unset, z, v, null, massPercent, double.NaN, n * scale));
            }
            else if (v > 0)
            {
                sumPositive += n * v;
                //components.Add(new Row(RowKind.Component, z, v, ComponentFormula(AtomStatic.AtomicName(z), v, AtomStatic.AtomicName(baseZ), k), // 261025Cl 変更前
                //    100 * n * (AtomStatic.AtomicWeight(z) + v * weightB / k) / massTotal, double.NaN, n * scale));
                // 261026Cl 変更 (作者指示): 小さい整数比で書けない価数 (Mn 2.2143 など) は、陽イオン 1 個あたりの陰イオンの数を小数 3 桁で書く (MnO1.107)。
                // 以前は null を返して、表示側で「Mn–O 化合物 (価数 2.214)」としていた
                var formula = ComponentFormula(AtomStatic.AtomicName(z), v, AtomStatic.AtomicName(baseZ), k)
                    ?? $"{AtomStatic.AtomicName(z)}{AtomStatic.AtomicName(baseZ)}{(v / k).ToString("0.000", CultureInfo.InvariantCulture)}";
                components.Add(new Row(RowKind.Component, z, v, formula,
                    100 * n * (AtomStatic.AtomicWeight(z) + v * weightB / k) / massTotal, double.NaN, n * scale));
            }
            else if (v < 0)
            {
                sumNegative += n * -v;
                anions.Add(new Row(RowKind.Anion, z, v, null, massPercent, double.NaN, n * scale));
                if (n != 0)
                    correctionZ.Add(z); // 261025Cl 追加: 原子数 0 (占有率 0 のサイトだけ) の陰イオンは補正に寄与しないので「−B≡X」に入れない
            }
            else
                massOnly.Add(new Row(RowKind.MassOnly, z, v, null, massPercent, double.NaN, n * scale));
        }

        var rows = new List<Row>();
        rows.AddRange(components);
        rows.AddRange(massOnly);
        rows.AddRange(anions);
        //if (anions.Count > 0) // 261025Cl 変更前 (原子数 0 の陰イオンだけでも 0 の補正行が出ていた)
        if (correctionZ.Count > 0)
        {
            var equivalent = sumNegative / k; // 補正に当たる基準陰イオンの原子数
            rows.Add(new Row(RowKind.Correction, baseZ, double.NaN, null, -100 * weightB * equivalent / massTotal, double.NaN, -equivalent * scale));
        }
        rows.Add(new Row(RowKind.Total, 0, double.NaN, null, rows.Sum(r => r.WeightPercent), double.NaN, double.NaN));

        var residual = nB - (sumPositive - sumNegative) / k; // 正なら過剰
        rows.Add(new Row(RowKind.Residual, baseZ, double.NaN, null, 100 * weightB * residual / massTotal, double.NaN, residual * scale));

        return new Result
        {
            Mode = mode,
            BaseZ = baseZ,
            //CorrectionAnions = [.. anions.Select(r => r.Z)], // 261025Cl 変更前
            CorrectionAnions = [.. correctionZ],
            Rows = rows,
            ResidualAtomsPerCell = residual,
            HasUnset = hasUnset,
            //ChargeSum = hasUnset ? double.NaN : (sumPositive - sumNegative - k * nB) * scale, // = −k ΔN_B × 倍率 // 261026Cl 変更前
            ChargeSum = hasUnset ? double.NaN : sumPositive - sumNegative - k * nB, // = −k ΔN_B (単位胞あたり。261026Cl 変更: 倍率を掛けない)
        };
    }
}
