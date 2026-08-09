// 260809Cl 新規作成: 角度広がりの畳み込み (設計 §3.6、指示書 §1-1 の ①)。
//
// engine は触らない。畳み込みは**方位軸上の後処理**なので、RunAlchemi が返した曲線に対して
// ここで掛ける。⚠処理順は **「畳み込み → 規格化」に固定** (作者決定 260809Cl)。
// 規格化された曲線を畳み込むと、規格化の基準 (走査平均) 自体が畳み込み前の値のままになって
// ICP の意味が変わるため。
//
// v1 は None / Gaussian1D(FWHM) のみ。離散カーネル CSV は公開後 (設計 §3.6)。
//
// 実装は素直な直接畳み込み。走査点は数百なので FFT は不要 (設計 §3.6「数百点なら直接畳み込み」)。
// **端部はカーネルを再規格化する** = 出力 i の重みを Σ_j w_ij で割る。走査の外側を 0 と見なすと
// 端が谷になって偽の構造が出るので、これは必須。副作用として、非等間隔の方位列でも
// (Rocking1D は等間隔だが) そのまま正しく動く。

using System;

namespace Crystallography;

/// <summary>260809Cl 追加: 角度広がりカーネルの種類 (設計 §3.6)。</summary>
public enum AlchemiSpreadKernel
{
    /// <summary>畳み込まない (生の forward 計算)。</summary>
    None,
    /// <summary>1 次元ガウシアン。幅は FWHM [mrad] で与える。</summary>
    Gaussian1D,
}

/// <summary>260809Cl 追加: ロッキング曲線への角度広がり畳み込み。WinForms 非依存。</summary>
public static class AlchemiAngularSpread
{
    /// <summary>FWHM → ガウシアン指数の係数。exp(−4 ln2 · Δ²/FWHM²) が FWHM で 1/2 になる。</summary>
    private const double FourLn2 = 4 * 0.693147180559945309;

    /// <summary>1 本の曲線に角度広がりを掛ける。<paramref name="fwhmRad"/> ≤ 0 なら入力をそのまま返す。</summary>
    /// <param name="curve">方位軸上の値 (長さ = 走査点数)</param>
    /// <param name="tiltRad">各点の傾斜角 [rad] (curve と同じ長さ・同じ順)</param>
    /// <param name="fwhmRad">ガウシアンの半値全幅 [rad]</param>
    public static double[] Gaussian(double[] curve, double[] tiltRad, double fwhmRad)
    {
        ArgumentNullException.ThrowIfNull(curve);
        ArgumentNullException.ThrowIfNull(tiltRad);
        if (curve.Length != tiltRad.Length)
            throw new ArgumentException("curve and tiltRad must have the same length", nameof(tiltRad));
        if (!double.IsFinite(fwhmRad)) throw new ArgumentException("FWHM must be finite", nameof(fwhmRad));
        if (fwhmRad <= 0 || curve.Length < 2) return curve;

        var a = FourLn2 / (fwhmRad * fwhmRad);
        //4σ (≈1.7 FWHM) より外は重み <1e-7 なので切る。O(N²) でも数百点なら十分速いが、
        //1001 点 × 全厚み × 全サイト × 全チャネルを毎回描き直すので効かせておく
        var cutoff = 1.7 * fwhmRad;
        var result = new double[curve.Length];
        for (int i = 0; i < curve.Length; i++)
        {
            double sum = 0, weight = 0;
            for (int j = 0; j < curve.Length; j++)
            {
                var d = tiltRad[i] - tiltRad[j];
                if (Math.Abs(d) > cutoff) continue;
                var w = Math.Exp(-a * d * d);
                sum += w * curve[j];
                weight += w;
            }
            //weight は自分自身 (w=1) を必ず含むので 0 にはならない
            result[i] = sum / weight;
        }
        return result;
    }
}
