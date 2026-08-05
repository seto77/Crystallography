// 260805Cl 新規作成: 菊池線動力学化 Phase 0-1 (設計正本 = ReciPro_菊池線動力学化設計.md §3)。
// TDS 源カーネルの抽象境界。v1 は「伝播 = getU 虚部 / 源 = 独立構成の Q」の暫定運用 (設計 §3「吸収と源の整合」)。
// 将来の整合カーネル (EvaluateSource/EvaluateLoss を同居させ ∫Q dΩdE ↔ −2 Im H の収支を取る) は
// このインターフェースを拡張して追加する。
// ⚠ Diffraction/Kikuchi/ フォルダは WinForms / System.Drawing 非依存を規律で維持 (設計 §5)。

using System;

namespace Crystallography;

/// <summary>菊池バンド計算の非弾性「源」カーネル。260805Cl 追加</summary>
public interface IKikuchiInelasticKernel
{
    /// <summary>
    /// 原子種 atomsIndex (Crystal.Atoms の添字) の TDS 源振幅 τ_a(q) を返す。
    /// s2 = |q|²/4 [nm⁻²] (散乱因子表の引数と同じ規約)。
    /// </summary>
    double SourceAmplitude(int atomsIndex, double s2);

    /// <summary>
    /// 260805Cl 追加: 源密度行列 Q の非対角 (原子局在由来の源コヒーレンス) の実振幅。
    /// s2i, s2j = 2 波の |q_i|²/4, |q_j|²/4、s2ij = |q_i − q_j|²/4 (2 ビームでは = |g|²/4)。
    /// 既定実装は因子化近似 τ_i·τ_j (設計 §3 の最小モデル)。厳密な Einstein 混合動的形状因子は
    /// KikuchiTdsEinsteinKernel が override する (Omoto 2002 eq. 38 = Hall–Hirsch)。
    /// 対角との整合: SourceCoherence(s2, s2, 0) = SourceAmplitude(s2)² を満たすこと。
    /// </summary>
    double SourceCoherence(int atomsIndex, double s2i, double s2j, double s2ij)
        => SourceAmplitude(atomsIndex, s2i) * SourceAmplitude(atomsIndex, s2j);
}

/// <summary>
/// 独立原子・等方 Einstein の TDS 源 (設計 §3 の最小モデル)。
/// τ_a(q) = √T_a(q), T_a(q) = |f_a(q)|²[1 − e^{−2M_a(q)}], M_a = B_iso·s² (s² = |q|²/4)。
/// 260805Cl 追加
/// </summary>
public sealed class KikuchiTdsEinsteinKernel : IKikuchiInelasticKernel
{
    private readonly Func<double, double>[] _factor;
    private readonly double[] _biso;

    public KikuchiTdsEinsteinKernel(Crystal crystal)
    {
        int n = crystal.Atoms.Length;
        _factor = new Func<double, double>[n];
        _biso = new double[n];
        for (int i = 0; i < n; i++)
        {
            var atoms = crystal.Atoms[i];
            // 弾性散乱因子は getU と同じ ElasticIonModel 規約 (Neutral は中性エントリを強制)。
            // IonFull の g≠0 単極子は連続 q の TDS 源には持ち込まない (有限部分のみ)
            int sub = BetheMethod.ElasticIonModel == ElasticIonModel.Neutral ? 0 : atoms.SubNumberElectron;
            _factor[i] = AtomStatic.ElectronScatteringPeng[atoms.AtomicNumber][sub].Factor;
            var dsf = atoms.Dsf;
            var b = double.IsNaN(dsf.Biso) ? dsf.Biso000 : dsf.Biso; // v1 は等方 ADP のみ (設計 §3)
            _biso[i] = double.IsNaN(b) ? 0 : b;
        }
    }

    public double SourceAmplitude(int atomsIndex, double s2)
    {
        var f = _factor[atomsIndex](s2);
        var t = f * f * (1 - Math.Exp(-2 * _biso[atomsIndex] * s2));
        return t > 0 ? Math.Sqrt(t) : 0;
    }

    /// <summary>
    /// 260805Cl 追加: 厳密な Einstein 混合動的形状因子 (Omoto-Tsuda-Tanaka 2002 eq. 38 = Hall–Hirsch):
    /// f(q_i) f(q_j) [e^{−W(q_i−q_j)} − e^{−W(q_i)−W(q_j)}], W(q) = B s² (s² = |q|²/4)。
    /// 対角 (s2i=s2j=s², s2ij=0) は f²(1−e^{−2Bs²}) = SourceAmplitude² に一致。
    /// 因子化近似 τ_i τ_j との差は q_i ∦ q_j 成分 (2W(g) 相当) の分だけ。
    /// </summary>
    public double SourceCoherence(int atomsIndex, double s2i, double s2j, double s2ij)
    {
        var b = _biso[atomsIndex];
        return _factor[atomsIndex](s2i) * _factor[atomsIndex](s2j)
             * (Math.Exp(-b * s2ij) - Math.Exp(-b * (s2i + s2j)));
    }
}

/// <summary>診断用の一様源 (τ = 1)。E/D 消失テスト専用 (設計 §3, §6)。260805Cl 追加</summary>
public sealed class KikuchiUniformSourceKernel : IKikuchiInelasticKernel
{
    public double SourceAmplitude(int atomsIndex, double s2) => 1.0;
}
