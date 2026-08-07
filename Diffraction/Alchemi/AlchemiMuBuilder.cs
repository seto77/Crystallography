// 260807Cl 新規作成: ALCHEMI 局所イオン化行列 μ の **production builder** (A2′)。
// 設計正本 = .project-guidance/ReciPro/ReciPro_ALCHEMI設計.md §3.1 (LocalFormFactor 階層) / §3.3 (DWF) / §5.3 (計算パス)。
//
//   μ_hg = Σ_a Occ_a · e^{−M_a(G)} · σ_c · F_c(|G|/2) · e^{−2πi G·r_a},   G = g_h − g_g
//
// 添字は h = row (bra) / g = col (ket)、格納は column-major (mu[row + col*n])。既存の
// getPotentialMatrix / EBSDSolver と同じレイアウト契約なので、C†μC の縮約へそのまま渡せる。
//
// ⚠ この Diffraction/Alchemi/ フォルダは WinForms / System.Drawing 非依存を規律で維持する
//   (ロードマップのファイル所有権: solver 本体は新規ファイル、BetheMethod.cs は薄い結線のみ)。
// ⚠ IonizationChannel.cs (共通型) は **消費のみ**。ALCHEMI 側で型を追加・変更しない。
//
// STEM-EDX の BetheMethod.ComputeIonizationU との関係 (意図的に別実装):
//   ・EDX は「元素 Z が一致する全原子」を 1 本に潰し、k_vac/(2πV_cell) を掛けた**ポテンシャル** U_ion を返す。
//   ・ALCHEMI は「サイト仮説ごと」に分けた**断面積行列** μ を返す (単位 nm²)。V_cell・k_vac の因子は
//     yield 側の規格化 (設計 §3.4 の μ00 = m_s·Occ_s·σ_s/V_c) が持つ。μ の対角が Σ_a Occ_a σ に
//     なるという契約を保つため、ここでは掛けない。
//   ・**元素フィルタを掛けない**のが最大の違い: tracer 仮説では「ホストサイトの幾何 (r_a, ADP, Occ) ×
//     ドーパントのチャネル (σ, F)」という組み合わせが正しい (設計 §3.5)。
//
// ⚠設計 §3.3 の禁止事項: **積近似 e^{−M(g)}·e^{−M(h)} を使ってはならない**。正しい熱平均は
//   ⟨e^{2πi(g−h)·u}⟩ = e^{−M(g−h)} で、積形にすると対角が σ·e^{−2M} < σ になり総断面積と
//   1 ビーム極限の規格化が壊れる。AlchemiCheck の「対角 = μ00」検査がこれを捕まえる。

using System;
using System.Collections.Generic;
using System.Numerics;

namespace Crystallography;

/// <summary>260807Cl 追加: 物理モデル階層 (設計 §3.1)。</summary>
public enum AlchemiModelTier
{
    /// <summary>δ 局在・静止原子 (F≡1, DWF≡1)。**GUI 非公開の内部 oracle** — 縮約経路の検証専用。</summary>
    StaticDelta,
    /// <summary>σ×F(s) + DWF の局所形状因子。v1 の公開モデル。</summary>
    LocalFormFactor,
    /// <summary>二変数 (Q_bra, Q_ket) カーネル。v2・別プロジェクト (一変数 F(s) からは生成不能)。</summary>
    FullMdff,
}

/// <summary>260807Cl 追加: サイト仮説の変位 (DWF) の与え方 (設計 §5.2 の SiteDisplacementSpec)。</summary>
public enum SiteDisplacementKind
{
    /// <summary>ホスト原子の Dsf をそのまま使う (tracer の既定。ドーパントは母格子の振幅で振れると仮定)。</summary>
    InheritHost,
    /// <summary>等方 Biso [nm²] で上書き。</summary>
    IsotropicB,
    /// <summary>非等方 B11..B31 (getU / Crystal.cs と同じ「B 型」無次元係数) で上書き。</summary>
    AnisotropicB,
}

/// <summary>260807Cl 追加: サイト仮説の変位指定。<see cref="InheritHost"/> 以外は
/// ホストの Dsf を無視してこの値を使う (設計 §5.2)。</summary>
public sealed record SiteDisplacement(SiteDisplacementKind Kind,
    double BisoNm2 = 0, double B11 = 0, double B22 = 0, double B33 = 0, double B12 = 0, double B23 = 0, double B31 = 0)
{
    /// <summary>ホスト原子の Dsf を継承する (既定)。</summary>
    public static readonly SiteDisplacement InheritHost = new(SiteDisplacementKind.InheritHost);

    /// <summary>等方 Biso [nm²] で上書き。</summary>
    public static SiteDisplacement Isotropic(double bisoNm2) => new(SiteDisplacementKind.IsotropicB, bisoNm2);

    /// <summary>非等方 B (無次元、getU:3167-3186 と同じ規約) で上書き。</summary>
    public static SiteDisplacement Anisotropic(double b11, double b22, double b33, double b12, double b23, double b31)
        => new(SiteDisplacementKind.AnisotropicB, 0, b11, b22, b33, b12, b23, b31);
}

/// <summary>260807Cl 追加: μ を組む対象の原子軌道集合 = サイト仮説の幾何部分 (設計 §5.2 の SiteHypothesis から
/// μ 構築に必要な部分だけを取り出したもの)。**Crystal.Atoms を書き換えない**のが規律で、空孔・アンチサイト・
/// 格子間もすべて「候補サイト上の成分分率」として表す。</summary>
/// <param name="Label">表示名 (凡例・診断・export に載る)</param>
/// <param name="AtomsIndices">この基底に属する <see cref="Crystal.Atoms"/> のインデックス列 (軌道単位)</param>
/// <param name="OccupancyFraction">占有率の上書き。null = 結晶側の Occ をそのまま使う。
/// Tracer 基底 (設計 §3.5) では 1.0 を渡し、実際の占有率は縮約後に線形結合で掛ける</param>
/// <param name="Displacement">変位 (DWF) の与え方。null = <see cref="SiteDisplacement.InheritHost"/></param>
public sealed record AlchemiSiteBasis(string Label, int[] AtomsIndices,
    double? OccupancyFraction = null, SiteDisplacement Displacement = null);

/// <summary>
/// 260807Cl 追加: 固定した反射集合 (FixedUnion 基底) に対する μ_hg の builder。
///
/// **μ は方位に依存しない** (結晶固定量。位相は整数 hkl と分率座標だけで決まり、|G| は回転不変)。
/// したがって FixedUnion 基底の run では (サイト × チャネル) ごとに **1 回だけ**組めばよく、
/// 方位ループの中では C†μC の縮約しか回らない (設計 §5.3)。
///
/// ΔG = g_h − g_g の重複は基底だけで決まるので、構築時に一度だけ unique 化して pair→slot 表を持つ。
/// n=200 本の基底で n² = 40000 ペアに対し unique ΔG は数千 — F(s) の batch 評価も原子和も
/// unique 側の回数で済み、全 (サイト×チャネル) で表を共有できる。
/// </summary>
public sealed class AlchemiMuBuilder
{
    private const double TwoPi = 2 * Math.PI;

    private readonly Crystal _crystal;
    private readonly (int H, int K, int L)[] _diff;   // unique ΔG (正準整数指数)
    private readonly int[] _slot;                     // [row + col*n] → _diff の添字
    private readonly double[] _s;                     // |ΔG|/2 [nm⁻¹] (unique 側)
    private readonly double[] _s2;                    // s² [nm⁻²]

    /// <summary>基底の反射本数 (μ は n×n)。</summary>
    public int BeamCount { get; }

    /// <summary>unique な ΔG の本数 (診断値。n² に対する圧縮率がそのまま構築コストの比になる)。</summary>
    public int UniqueDifferenceCount => _diff.Length;

    /// <summary>結晶と固定基底から ΔG 表を作る。gIndices は run 中不変であること
    /// (基底が変わったら builder を作り直す = FixedUnion の規律をコード側で強制する)。</summary>
    public AlchemiMuBuilder(Crystal crystal, (int H, int K, int L)[] gIndices)
    {
        ArgumentNullException.ThrowIfNull(crystal);
        ArgumentNullException.ThrowIfNull(gIndices);
        _crystal = crystal;
        int n = BeamCount = gIndices.Length;

        _slot = new int[n * n];
        var slotOf = new Dictionary<(int H, int K, int L), int>(n * 4);
        var diff = new List<(int H, int K, int L)>(n * 4);
        for (int col = 0; col < n; col++)
            for (int row = 0; row < n; row++)
            {
                var d = (gIndices[row].H - gIndices[col].H, gIndices[row].K - gIndices[col].K, gIndices[row].L - gIndices[col].L);
                if (!slotOf.TryGetValue(d, out var slot))
                {
                    slot = diff.Count;
                    diff.Add(d);
                    slotOf.Add(d, slot);
                }
                _slot[row + col * n] = slot;
            }
        _diff = [.. diff];

        //G は combined hkl から正準再構築する (Cartesian 差の浮動小数和はビーム集合の並び順で
        //最下位ビットが揺れるため使わない。STEM-EDX §5.7-4 と同じ規律)。
        var mat = crystal.MatrixInverseTransposed;
        _s = new double[_diff.Length];
        _s2 = new double[_diff.Length];
        for (int i = 0; i < _diff.Length; i++)
        {
            var v = mat * _diff[i];
            _s2[i] = v.Length2 / 4;      //(|G|/2)² [nm⁻²] (getU:3143・ComputeIonizationU と同一定義)
            _s[i] = Math.Sqrt(_s2[i]);
        }
    }

    /// <summary>μ の対角値 μ00 = Σ_a Occ_a σ [nm²] (等価位置の多重度を含む)。
    /// dechannelling 項 (設計 §3.4) の μ00 = m_s·Occ_s·σ_s と同じ量 (V_c は yield 側)。</summary>
    public double Mu00(IonizationData channel, AlchemiSiteBasis site)
    {
        ArgumentNullException.ThrowIfNull(channel);
        ArgumentNullException.ThrowIfNull(site);
        double occSum = 0;
        foreach (var index in site.AtomsIndices)
        {
            var atoms = _crystal.Atoms[index];
            occSum += (site.OccupancyFraction ?? atoms.Occ) * atoms.Atom.Length;
        }
        return occSum * channel.TotalCrossSectionNm2;
    }

    /// <summary>1 つの (サイト, チャネル) に対する μ 行列 [nm²] を column-major (mu[row + col*n]) で返す。
    /// row = h (bra) / col = g (ket)、G = g_h − g_g。</summary>
    /// <param name="channel">run 開始時に解決済みの IonizationData (σ と F(s) の出所)</param>
    /// <param name="site">サイト仮説の幾何部分</param>
    /// <param name="tier">モデル階層。<see cref="AlchemiModelTier.StaticDelta"/> は F≡1・DWF≡1 の内部 oracle</param>
    public Complex[] Build(IonizationData channel, AlchemiSiteBasis site, AlchemiModelTier tier = AlchemiModelTier.LocalFormFactor)
    {
        ArgumentNullException.ThrowIfNull(channel);
        ArgumentNullException.ThrowIfNull(site);
        if (tier == AlchemiModelTier.FullMdff)
            throw new NotSupportedException("AlchemiModelTier.FullMdff needs a two-momentum kernel K(Q_bra,Q_ket); it cannot be generated from the one-variable F(s) table (設計 §3.1)");

        int u = _diff.Length;
        //F(s) は batch API で評価する契約 (N² 内の virtual call 回避)。unique ΔG の本数だけで済む
        var f = new double[u];
        if (tier == AlchemiModelTier.StaticDelta)
            Array.Fill(f, 1.0);
        else
            channel.Shape.Evaluate(_s, f);

        //unique ΔG ごとに「サイト内の原子和」を作る。σ·F は原子に依らないので最後に 1 回だけ掛ける
        //(ComputeIonizationU の structureSum と同じ組み立て順)。
        //軌道 (Atoms エントリ) を外側に置き、占有率と DWF 係数の解決を ΔG ループの外へ出す
        var value = new Complex[u];
        foreach (var atomsIndex in site.AtomsIndices)
        {
            var atoms = _crystal.Atoms[atomsIndex];
            var occ = site.OccupancyFraction ?? atoms.Occ;
            var dwf = Resolve(atoms, site.Displacement);
            //DWF exp(−M_a(G)) は getU:3167-3186 / ComputeIonizationU と同一の m 算出 (等方 Biso / 非等方 B11..B31)。
            //TDS の (1−exp) 構造や m==0 → imag=0 の分岐は持ち込まない (§3.2: G=0 で 1、総断面積は熱振動で消えない)
            var zero = dwf.IsZero || tier == AlchemiModelTier.StaticDelta;
            for (int i = 0; i < u; i++)
            {
                var index = _diff[i];
                double s2 = _s2[i];
                double m = zero ? 0 : double.NaN;
                foreach (var atom in atoms.Atom)
                {
                    if (!zero && ((!dwf.UseIso && index != (0, 0, 0)) || double.IsNaN(m)))//非等方でg≠0の時、あるいは初めての時
                    {
                        if (dwf.UseIso)
                            m = dwf.Biso;
                        else if (index == (0, 0, 0))
                            m = double.IsNaN(dwf.Biso) ? dwf.Biso000 : dwf.Biso;
                        else
                        {
                            //等価位置ごとに対称操作で指数を変換する = 非等方テンソルを軌道全体へ正しく展開する
                            var (H, K, L) = atom.Operation.ConvertPlaneIndex(index);
                            m = (dwf.B11 * H * H + dwf.B22 * K * K + dwf.B33 * L * L
                                + 2 * dwf.B12 * H * K + 2 * dwf.B23 * K * L + 2 * dwf.B31 * L * H) / s2;
                        }
                        if (double.IsNaN(m))
                            m = 0;
                    }
                    //負位相 e^{−2πiG·r_a} を直接生成 (既存 getU の「正位相を作って共役」規約を持ち込まない)
                    var (sin, cos) = Math.SinCos(TwoPi * (atom * index));
                    var t = m == 0 ? occ : occ * Math.Exp(-m * s2);
                    value[i] += new Complex(t * cos, -t * sin);
                }
            }
        }
        var sigma = channel.TotalCrossSectionNm2;
        for (int i = 0; i < u; i++)
            value[i] *= sigma * f[i];

        int n = BeamCount;
        var mu = GC.AllocateUninitializedArray<Complex>(n * n);
        for (int i = 0; i < mu.Length; i++)
            mu[i] = value[_slot[i]];
        return mu;
    }

    /// <summary>DWF の実効係数 (ホストの Dsf か、サイト仮説による上書き)。
    /// <see cref="Atoms.Dsf"/> と同じフィールド意味を持つ値型に落として、以降の分岐を 1 本にする。</summary>
    private readonly record struct Dwf(bool IsZero, bool UseIso, double Biso, double Biso000,
        double B11, double B22, double B33, double B12, double B23, double B31);

    private static Dwf Resolve(Atoms atoms, SiteDisplacement displacement)
    {
        switch (displacement?.Kind ?? SiteDisplacementKind.InheritHost)
        {
            case SiteDisplacementKind.IsotropicB:
                return new Dwf(displacement.BisoNm2 == 0, true, displacement.BisoNm2, displacement.BisoNm2, 0, 0, 0, 0, 0, 0);
            case SiteDisplacementKind.AnisotropicB:
                var d = displacement;
                var zero = d.B11 == 0 && d.B22 == 0 && d.B33 == 0 && d.B12 == 0 && d.B23 == 0 && d.B31 == 0;
                //Biso/Biso000 は G=0 の分岐でしか読まれず、そこは s²=0 なので e^{−m·s²}=1 になる
                //(= μ の対角は必ず Σ Occ σ。§3.2「総断面積は熱振動で消えない」がここで担保される)。
                //上書き指定には格子定数が無く Hamilton 式 (Acta Cryst. 12, 609) の Biso000 を作れないが、値に影響しない
                return new Dwf(zero, false, 0, 0, d.B11, d.B22, d.B33, d.B12, d.B23, d.B31);
            default:
                var dsf = atoms.Dsf;
                return new Dwf(dsf.IsZero, dsf.UseIso, dsf.Biso, dsf.Biso000,
                    dsf.B11, dsf.B22, dsf.B33, dsf.B12, dsf.B23, dsf.B31);
        }
    }
}
