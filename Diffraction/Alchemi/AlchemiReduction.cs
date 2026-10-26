// 260807Cl 新規作成: ALCHEMI の **modal 縮約カーネル** (A2′、設計 §3.2 / §3.4)。
// 1 方位ぶんの固有系 (γ_j, C_g^{(j)}, α_j) を受け取り、(サイト × チャネル) の μ 行列を
// 解析的厚み積分まで畳み込む。方位ループの内側で回る唯一の O(b³) がここ。
//
//   ψ_g(z) = Σ_j α_j C_g^{(j)} e^{2πi γ_j z}                     … Bloch 波展開
//   dP/dz  = ψ†(z) μ ψ(z)                                         … 局所イオン化率 [nm²]
//   P(t)   = Σ_{j,j'} S_{j'j} F_{jj'}(t),  S = D†μD, D = C·diag(α)
//   F_{jj'}(t) = [e^{λ_{jj'}t} − 1] / λ_{jj'},  λ_{jj'} = 2πi(γ_j − γ_{j'}*)
//
// StaticDelta 極限 (μ_hg = σ Σ_n e^{−2πi(g_h−g_g)·r_n}) では
//   ψ†μψ = σ Σ_n |ψ(r_n)|²
// となり、EBSDSolverManaged (BetheMethod.cs:1745-1868) の S = B†diag(σ)B と同じ量になる
// = 既存経路が回帰 oracle として使える (設計 §6.1-1、AlchemiCheck reduce が実測)。
//
// dechannelling (設計 §3.4、ICSC 式(1)):
//   Y_dech(t) = (μ00/V_c)·{t − L_coh(t)},  L_coh(t) = ∫₀ᵗ Σ_g |ψ_g(z)|² dz = **μ = I の同じ縮約**
// TDS 吸収で coherent Bloch 部分空間から失われた電子を、方向ランダム化された電子として
// 残り厚みぶん random 方位相当の断面積で発生させる。吸収ゼロ → Σ_g|ψ_g|² ≡ 1 → L_coh = t → Y_dech = 0。
// **t − L_coh が有意に負なら clamp せず hard fail** (基底・符号・規格化エラーの検出器、§3.4)。
//
// ⚠ 260811Cl: t − L_coh は虚部が作る減衰を**出所を問わず全部**拾う。だから再注入が正当なのは
// 虚部が TDS だけのとき (= 電子がエネルギーを保ったまま試料内に残るとき) に限られる。
// この前提を型で持たせたのが <see cref="AbsorptionSource"/> で、コンストラクタで宣言し
// <see cref="AlchemiReduction.Yield"/> が検査する。現行の虚部は TDS のみ (ICSC 2003 との照合で二重計上が
// 無いことを確認済み) なので、既定値のままでは挙動は変わらない — 効くのは将来 mean absorption や
// 経験的 damping を混ぜたときで、そのとき黙って過大評価にならずに落ちる。
//
// ⚠ この Diffraction/Alchemi/ フォルダは WinForms / System.Drawing 非依存を規律で維持する。
// ⚠ 既存 CBED / STEM / EBSD の worker・バッファ・総和順序には触らない (別パス方式、設計 §8)。

using System;
using System.Numerics;

namespace Crystallography;

/// <summary>260807Cl 追加: 1 サイト × 1 チャネルぶんの yield [入射電子あたりの発生イオン化数、無次元]。
/// Dynamic / Dechannelled / Total を分離保存する (設計 §3.4: 単一の nuisance 定数では吸収できない)。</summary>
public sealed class AlchemiSiteYield
{
    /// <summary>動力学項 (Bloch 波の定在波が作るチャネリング成分)</summary>
    public double[] Dynamic { get; }
    /// <summary>非チャネリング項 (吸収で coherent 部分空間から抜けた電子の寄与)</summary>
    public double[] Dechannelled { get; }
    /// <summary>Dynamic + Dechannelled。表示既定だが、**分離した 2 本も必ず保存する** (設計 §3.4:
    /// t−L_coh は厚み・方位・結晶依存なので、fit の定数 nuisance では吸収できない)</summary>
    public double[] Total { get; }
    /// <summary>L_coh(t) [nm] — 診断・回帰用 (吸収ゼロなら t に一致)</summary>
    public double[] CoherentPathLengthNm { get; }

    internal AlchemiSiteYield(double[] dynamic, double[] dechannelled, double[] coherentPathLengthNm)
    {
        Dynamic = dynamic;
        Dechannelled = dechannelled;
        CoherentPathLengthNm = coherentPathLengthNm;
        Total = new double[dynamic.Length];
        for (int i = 0; i < Total.Length; i++) Total[i] = dynamic[i] + dechannelled[i];
    }
}

/// <summary>
/// 260807Cl 追加: 1 方位ぶんの縮約器。固有系と厚み一覧で構築し、(サイト × チャネル) の μ を次々に畳む。
/// D = C·diag(α) と λ⁻¹ は構築時に 1 回だけ作るので、チャネルが増えても再計算しない
/// (設計 §3.2 の「B を方位ごとに 1 回作り、サイト群へ振り分けて複数 S へ蓄積する」batch 設計)。
/// </summary>
public sealed class AlchemiReduction
{
    private static readonly Complex TwoPiI = new(0, 2 * Math.PI);

    private readonly int _b;
    private readonly Complex[] _d;          // D = C·diag(α)、column-major D[g + j*b]
    private readonly Complex[] _invLambda;  // 1/λ_{jj'}、[j' + j*b] (縮退は 0)
    private readonly bool[] _degenerate;    // λ ≈ 0 (ロピタルで F → t)
    private readonly Complex[] _eigen;
    private readonly double[] _t;
    private readonly AbsorptionSource _absorption;// 虚部の出所 (再注入してよいのは TDS だけ)
    private readonly double[] _pOverP0;     // 260829Cl 追加: L_coh の流束重み P_g/P_0 (null = 重みなし)
    private double[] _lcoh;                 // L_coh(t) の遅延キャッシュ (方位ごとに 1 回で足りる)

    /// <summary>厚み一覧 [nm] (構築時に渡したもの)。</summary>
    public double[] ThicknessesNm => _t;

    /// <summary>ビーム数 b。</summary>
    public int BeamCount => _b;

    /// <param name="bLen">ビーム数 b</param>
    /// <param name="eigenValues">固有値 γ_j (b 個)。Im(γ_j) &gt; 0 が減衰</param>
    /// <param name="eigenVectors">固有ベクトル C_g^{(j)} (b×b, column-major: eigenVectors[j*b + g])</param>
    /// <param name="alpha">励起振幅 α_j (b 個)</param>
    /// <param name="thicknessesNm">厚み [nm]</param>
    /// <param name="absorption">260811Cl 追加: この固有系の**虚部に何が入っているか** (<see cref="AbsorptionSource"/>)。
    /// 非チャネリング項は失われた流束を「試料内に残るランダム方位の電子」として**まるごと**再注入するので、
    /// TDS 以外が混ざっていたら <see cref="Yield"/> が落ちる。既定は v1 の実態
    /// (<see cref="BetheMethod.ImaginaryPotentialAbsorption"/> = TDS のみ) で、既存の呼び出しは 1 ビットも変わらない</param>
    /// <param name="pOverP0">260829Cl 追加: 各ビームの P_g/P_0 (P_g = 2n̂·(k₀+g))。L_coh の生存流束を
    /// Σ_g (P_g/P_0)|ψ_g|² (表面法線方向の全電流。無吸収で厳密に 1) で数えるための重み。
    /// null なら従来どおり Σ_g|ψ_g|² (P_g ≈ P_0 の近似。系統列 ALCHEMI では差は ~10⁻³ 以下)</param>
    public AlchemiReduction(int bLen, Complex[] eigenValues, Complex[] eigenVectors, Complex[] alpha, double[] thicknessesNm,
        AbsorptionSource absorption = AbsorptionSource.TdsRedistributable, double[] pOverP0 = null)
    {
        ArgumentNullException.ThrowIfNull(eigenValues);
        ArgumentNullException.ThrowIfNull(eigenVectors);
        ArgumentNullException.ThrowIfNull(alpha);
        ArgumentNullException.ThrowIfNull(thicknessesNm);
        if (eigenValues.Length < bLen || alpha.Length < bLen || eigenVectors.Length < bLen * bLen)
            throw new ArgumentException($"AlchemiReduction: eigen system is smaller than bLen = {bLen}");
        _b = bLen;
        _eigen = eigenValues;
        _t = thicknessesNm;
        _absorption = absorption;
        if (pOverP0 is not null && pOverP0.Length < bLen)
            throw new ArgumentException($"AlchemiReduction: pOverP0 is smaller than bLen = {bLen}");
        _pOverP0 = pOverP0; // 260829Cl 追加

        _d = GC.AllocateUninitializedArray<Complex>(bLen * bLen);
        for (int j = 0; j < bLen; j++)
        {
            var a = alpha[j];
            for (int g = 0; g < bLen; g++)
                _d[g + j * bLen] = eigenVectors[j * bLen + g] * a;
        }

        //λ_{jj'} = 2πi(γ_j − γ_{j'}*) を逆数で持つ (縮約の内側では除算をしない)。
        //λ ≈ 0 (縮退 = 吸収ゼロの対角など) は F → t (ロピタル)。判定は EBSDSolverManaged:1876 と同じ基準
        _invLambda = GC.AllocateUninitializedArray<Complex>(bLen * bLen);
        _degenerate = new bool[bLen * bLen];
        for (int j = 0; j < bLen; j++)
            for (int jp = 0; jp < bLen; jp++)
            {
                var lam = TwoPiI * (eigenValues[j] - Complex.Conjugate(eigenValues[jp]));
                int i = jp + j * bLen;
                _degenerate[i] = lam.Real * lam.Real + lam.Imaginary * lam.Imaginary < 1e-30;
                _invLambda[i] = _degenerate[i] ? Complex.Zero : Complex.One / lam;
            }
    }

    /// <summary>μ 行列 (column-major、row = bra) を厚みごとに畳んだ Re Σ_{jj'} S_{j'j} F_{jj'}(t) を返す。
    /// 単位は μ のそれ × nm (μ が nm² なら nm³)。物理量への換算 (1/V_cell) は呼び出し側。
    /// **負値を clamp しない** — 打切り・符号・規格化の誤りを潰さないため (設計 §3.4)。</summary>
    public double[] Contract(Complex[] mu)
    {
        ArgumentNullException.ThrowIfNull(mu);
        if (mu.Length < _b * _b) throw new ArgumentException($"AlchemiReduction.Contract: μ is smaller than {_b}×{_b}");
        return ContractCore(BuildS(mu));
    }

    /// <summary>L_coh(t) = ∫₀ᵗ Σ_g (P_g/P_0)|ψ_g(z)|² dz [nm] (μ = diag(P_g/P_0) の同じ縮約)。
    /// 260829Cl 変更: 生存流束を表面法線方向の全電流 Σ (P_g/P_0)|ψ_g|² で数える (コンストラクタに pOverP0 を渡した場合)。
    /// 吸収ゼロならこの量が厳密に 1 で L_coh = t になる (行割り修正後の正しい保存量。
    /// 旧記述「C はユニタリで Σ_g|ψ_g|² ≡ 1」は P_g ≈ P_0 の近似でのみ成立)。pOverP0 = null なら従来の Σ|ψ_g|²。
    /// 方位だけで決まる量なので初回で作って以後使い回す (チャネル・サイトが増えても O(b³) は 1 回)。
    /// 返す配列は共有 = **書き換えないこと**。</summary>
    public double[] CoherentPathLengthNm()
    {
        //競合しても同じ値を 2 度作るだけ (決定的計算 + 参照代入) なので lock を張らない
        if (_lcoh is not null) return _lcoh;
        //μ = diag(P_g/P_0) (null なら I) なので S = D†μD。μ を実体化せずに直接組む (b² の一時行列を作らない)
        var s = new Complex[_b * _b];
        for (int j = 0; j < _b; j++)
            for (int jp = 0; jp < _b; jp++)
            {
                Complex acc = 0;
                if (_pOverP0 is null)
                    for (int g = 0; g < _b; g++)
                        acc += Complex.Conjugate(_d[g + jp * _b]) * _d[g + j * _b];
                else
                    for (int g = 0; g < _b; g++)
                        acc += Complex.Conjugate(_d[g + jp * _b]) * _pOverP0[g] * _d[g + j * _b];
                s[jp + j * _b] = acc;
            }
        return _lcoh = ContractCore(s);
    }

    /// <summary>1 サイト × 1 チャネルの yield (設計 §3.2 / §3.4)。
    /// Y_dyn = Contract(μ)/V_c、Y_dech = (μ00/V_c)·{t − L_coh(t)}。</summary>
    /// <param name="mu">サイト × チャネルの μ 行列 [nm²]</param>
    /// <param name="mu00Nm2">μ の対角 Σ_a Occ_a σ [nm²] (<see cref="AlchemiMuBuilder.Mu00"/>)</param>
    /// <param name="unitCellVolumeNm3">単位胞体積 [nm³] (<see cref="Crystal.Volume"/>)</param>
    /// <param name="includeDechannelled">false なら Dechannelled を 0 で返す (診断用。既定は v1 から必須)</param>
    public AlchemiSiteYield Yield(Complex[] mu, double mu00Nm2, double unitCellVolumeNm3, bool includeDechannelled = true)
    {
        if (!(unitCellVolumeNm3 > 0))
            throw new ArgumentOutOfRangeException(nameof(unitCellVolumeNm3), unitCellVolumeNm3, "unit cell volume must be positive");
        //260811Cl: 再注入してよいのは TDS だけ。t − L_coh は虚部が作る減衰を**出所を問わず全部**拾うので、
        //TDS 以外が混じったまま再注入すると、真の非弾性損失や経験的 damping まで
        //「まだ試料内にいてイオン化する電子」として数えてしまう (静かに、必ず過大の側へ)。
        //ここで落とすのは「吸収を混ぜるな」ではなく「混ぜたなら出所ごとに分けて再注入せよ」の意味。
        if (includeDechannelled && !_absorption.IsFullyRedistributable())
            throw new InvalidOperationException(
                $"ALCHEMI dechannelling: the absorptive potential declares {_absorption}, which is not purely "
                + $"{AbsorptionSource.TdsRedistributable}. The dechannelled term re-injects the entire flux lost from the "
                + "coherent field as randomly directed electrons that keep ionizing, and that is only defensible for thermal "
                + "diffuse scattering. Split the absorptive matrix by source and drive the re-injection from the TDS part "
                + "alone, or set includeDechannelled = false.");
        var raw = Contract(mu);
        var dyn = new double[raw.Length];
        for (int i = 0; i < raw.Length; i++) dyn[i] = raw[i] / unitCellVolumeNm3;

        var lcoh = CoherentPathLengthNm();//方位ごとに 1 回だけ作られる (以後キャッシュ)
        var dech = new double[raw.Length];
        if (includeDechannelled)
            for (int i = 0; i < raw.Length; i++)
            {
                var lost = _t[i] - lcoh[i];
                //設計 §3.4: clamp せず fail。t−L_coh は「吸収で coherent 部分空間から抜けた電子の走行距離和」で
                //定義上非負 — 有意に負なら基底・符号・規格化のどれかが壊れている。
                //しきい値 1e-6·t の根拠: 符号規約や規格化を取り違えると t と同じオーダーで負に振れる (実測は §6 の
                //Im(γ)<0 fixture) 一方、固有分解の丸めは条件数込みでもこの水準に届かない。
                //残った微小負値だけを 0 にする (欠陥は落とし、丸めは通す)
                if (lost < -1e-6 * Math.Max(_t[i], 1e-12))
                    throw new InvalidOperationException(
                        $"ALCHEMI dechannelling: t − L_coh = {lost:e3} nm is significantly negative at t = {_t[i]:f4} nm "
                        + "(L_coh must not exceed the thickness). Refusing to clamp — check the basis, the eigenvalue sign convention and the excitation normalization (設計 §3.4).");
                dech[i] = mu00Nm2 / unitCellVolumeNm3 * Math.Max(0, lost);
            }
        return new AlchemiSiteYield(dyn, dech, lcoh);
    }

    /// <summary>S = D†μD (column-major、S[j' + j*b] = Σ_{h,g} conj(D_{h,j'}) μ_{hg} D_{g,j})。</summary>
    private Complex[] BuildS(Complex[] mu)
    {
        int b = _b;
        //T = μD : T[h + j*b] = Σ_g μ[h + g*b] · D[g + j*b]
        var t = new Complex[b * b];
        for (int j = 0; j < b; j++)
            for (int g = 0; g < b; g++)
            {
                var d = _d[g + j * b];
                if (d.Real == 0 && d.Imaginary == 0) continue;
                int muCol = g * b, tCol = j * b;
                for (int h = 0; h < b; h++)
                    t[h + tCol] += mu[h + muCol] * d;
            }
        //S = D†T
        var s = new Complex[b * b];
        for (int j = 0; j < b; j++)
            for (int jp = 0; jp < b; jp++)
            {
                Complex acc = 0;
                for (int h = 0; h < b; h++)
                    acc += Complex.Conjugate(_d[h + jp * b]) * t[h + j * b];
                s[jp + j * b] = acc;
            }
        return s;
    }

    /// <summary>Re Σ_{j,j'} S_{j'j} F_{jj'}(t) を厚みごとに評価する。
    /// e^{λ_{jj'}t} = u_j(t)·conj(u_{j'}(t)) と分解するので、指数関数は厚みあたり b 回で済む
    /// (素直に組むと b² 回)。数学的には同一で、丸めの出方だけが違う。</summary>
    private double[] ContractCore(Complex[] s)
    {
        int b = _b;
        var u = new Complex[b];//インスタンスに持たせない = 同一方位のチャネルを並列に畳んでも安全
        var result = new double[_t.Length];
        for (int ti = 0; ti < _t.Length; ti++)
        {
            double thick = _t[ti];
            for (int j = 0; j < b; j++)
                u[j] = Complex.Exp(TwoPiI * _eigen[j] * thick);
            Complex sum = 0;
            for (int j = 0; j < b; j++)
            {
                var uj = u[j];
                int col = j * b;
                for (int jp = 0; jp < b; jp++)
                {
                    int i = jp + col;
                    var f = _degenerate[i] ? thick : (uj * Complex.Conjugate(u[jp]) - Complex.One) * _invLambda[i];
                    sum += s[i] * f;
                }
            }
            result[ti] = sum.Real;
        }
        return result;
    }
}
