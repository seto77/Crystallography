using System;
using System.Collections.Generic;

namespace Crystallography;

/// <summary>
/// 261106Cl 追加: 非弾性散乱の GOS 模型 (PENELOPE 2018 §3.2 に準拠: Liljequist の δ 振動子、Sternheimer の共鳴エネルギー、
/// 束縛殻の遠い衝突の三角分布、遠い縦・横の相互作用、交換を含む近い衝突 (Møller)、密度効果 δ_F)。
/// 損失 W と曲がり (運動量移行 Q からの散乱角) を同じ事象から引く。
/// <para>ReciPro との差分 (独自の調整係数は持たない):
/// (1) 価電子は 1 つの振動子 (PENELOPE の伝導帯の扱い: U=0、W_cb = √(f_cb/Z)·Ω_p、Q_cb = W_cb) にまとめる。f_cb は TPP-2M と同じ価電子数。
/// (2) 束縛殻はすべて三角分布と式 3.78・3.80 の閾値の修正で扱う。
/// (3) 呼び出し側 (<see cref="MonteCarlo"/>) が、内殻 (U ≥ <see cref="InnerShellThresholdEv"/>) の率を保ち外側の率を一律に倍して、
/// 既定の阻止能 S (相対論 Bethe + Joy–Luo) に規格化する。λ_in はその結果として出る (TPP-2M は比較値)。</para>
/// <para>261106Cl codex2 第 5 回の査読で直した点: 近い衝突の下限は束縛殻で U_k (修正後の Q'_k ではない)、価電子で W_cb /
/// 遠い衝突の損失に交換の上限 (E+U_k)/2 を入れ、率と阻止能は切った三角分布の 0 次・1 次モーメントで出す /
/// Møller の係数 a は (E/(E+mc²))² (E' = E+U_k は損失の比と上限にだけ使う。PenRed と同じ) / 元素ごとに集約してから振動子を作る /
/// Q₋ は桁落ちしない形で計算 / 棄却が尽きたら例外。</para>
/// <para>殻のデータ = <see cref="CarlsonAtomicShells"/> (Carlson 1975)。平均励起エネルギー I は呼び出し側の J (単体は ReciPro の式、化合物は Bragg 則)。</para>
/// </summary>
internal sealed class InelasticGos
{
    /// <summary>電子の静止エネルギー [eV]</summary>
    internal const double Mc2 = 510998.95;
    /// <summary>e²/(4πε₀) [eV·nm]</summary>
    private const double E2 = 1.43996448;
    /// <summary>ħc [eV·nm]</summary>
    private const double HbarC = 197.3269804;
    /// <summary>この値以上の束縛殻を「内殻」(率を保つ側) とする [eV]。追跡の下限 50 eV と揃える</summary>
    internal const double InnerShellThresholdEv = 50.0;

    /// <summary>事象の経路。振動子 k ごとに 3 つ (遠い縦・遠い横・近い衝突) が並ぶ</summary>
    internal enum Channel { DistantLongitudinal = 0, DistantTransverse = 1, Close = 2 }

    /// <summary>振動子の強さ (平均原子 1 個あたりの電子数)</summary>
    internal readonly double[] F;
    /// <summary>イオン化エネルギー U_k [eV] (価電子の振動子は 0)</summary>
    internal readonly double[] U;
    /// <summary>共鳴エネルギー W_k [eV] (修正前)</summary>
    internal readonly double[] W;
    /// <summary>内殻 (率を保つ側) か</summary>
    internal readonly bool[] IsInner;
    /// <summary>振動子の元素 (価電子の振動子は 0)</summary>
    internal readonly int[] ElementZ;
    /// <summary>振動子の殻の番号 (<see cref="CarlsonAtomicShells"/> の並び = K, L1, L2, L3, M1, …。価電子の振動子は −1)</summary>
    internal readonly int[] ShellIndex;
    /// <summary>振動子の数 (添字 0 が価電子の振動子)</summary>
    internal int Count => F.Length;
    /// <summary>平均原子 1 個あたりの電子数</summary>
    internal readonly double Z;
    /// <summary>平均原子の数密度 [1/nm³]</summary>
    internal readonly double AtomsPerNm3;
    /// <summary>全電子のプラズマエネルギー Ω_p [eV]</summary>
    internal readonly double OmegaP;
    /// <summary>平均励起エネルギー I [eV]</summary>
    internal readonly double MeanExcitationEv;
    /// <summary>Sternheimer の調整係数 a (式 3.64 の解)</summary>
    internal readonly double SternheimerA;

    /// <param name="composition">元素と原子数の重み (占有率×多重度)。Z は 1–99。同じ元素が複数あってよい (先に元素ごとに集約する)</param>
    /// <param name="electronsPerNm3">電子の数密度 [1/nm³] (呼び出し側の ρ·N_A·Z/A と揃える)</param>
    /// <param name="meanExcitationEv">平均励起エネルギー I [eV]</param>
    /// <param name="valenceElectrons">元素ごとの価電子数 (価電子の振動子へ回す電子数)</param>
    internal InelasticGos(IReadOnlyList<(int Z, double Count)> composition, double electronsPerNm3, double meanExcitationEv, Func<int, double> valenceElectrons)
    {
        // 元素ごとに集約 (Lorentz–Lorenz 項が f_k に依るので、同じ元素を分けて渡すと模型が変わるのを防ぐ)
        var byElement = new SortedDictionary<int, double>();
        foreach (var (zi, ni) in composition)
            if (ni > 0 && zi >= 1 && zi <= CarlsonAtomicShells.Electrons.Length)
                byElement[zi] = (byElement.TryGetValue(zi, out var v) ? v : 0) + ni;
        double sumN = 0, sumNZ = 0;
        foreach (var (zi, ni) in byElement) { sumN += ni; sumNZ += ni * zi; }
        if (!(sumN > 0)) throw new ArgumentException("InelasticGos: 組成が空");

        double fcb = 0;
        var bound = new List<(double f, double u, int z, int shell)>();
        foreach (var (zi, ni) in byElement)
        {
            double x = ni / sumN;
            var ue = CarlsonAtomicShells.BindingEnergiesEv[zi - 1];
            var ne = CarlsonAtomicShells.Electrons[zi - 1];
            // 束縛の弱い殻から価電子数ぶんを価電子の振動子へ回す (殻を割ってもよい)
            double nv = Math.Clamp(valenceElectrons(zi), 0, zi);
            fcb += x * nv;
            var order = new int[ue.Length];
            for (int i = 0; i < order.Length; i++) order[i] = i;
            Array.Sort(order, (p, q) => ue[p].CompareTo(ue[q]));
            double remaining = nv;
            var left = new double[ue.Length];
            for (int i = 0; i < ue.Length; i++) left[i] = ne[i];
            foreach (var i in order)
            {
                if (!(remaining > 0)) break;
                double take = Math.Min(left[i], remaining);
                left[i] -= take; remaining -= take;
            }
            for (int i = 0; i < ue.Length; i++)
                if (left[i] > 1e-12) bound.Add((x * left[i], ue[i], zi, i));
        }
        Z = sumNZ / sumN;
        AtomsPerNm3 = electronsPerNm3 / Z;
        OmegaP = Math.Sqrt(4 * Math.PI * electronsPerNm3 * E2 * HbarC * HbarC / Mc2);
        MeanExcitationEv = meanExcitationEv;

        int n = bound.Count + 1;
        F = new double[n]; U = new double[n]; W = new double[n]; IsInner = new bool[n]; ElementZ = new int[n]; ShellIndex = new int[n];
        F[0] = fcb; U[0] = 0; W[0] = fcb > 0 ? Math.Sqrt(fcb / Z) * OmegaP : 0; IsInner[0] = false; ShellIndex[0] = -1;
        for (int k = 1; k < n; k++)
        {
            (F[k], U[k], ElementZ[k], ShellIndex[k]) = bound[k - 1];
            IsInner[k] = U[k] >= InnerShellThresholdEv;
        }

        // 式 (3.64): Z ln I = f_cb ln W_cb + Σ f_k ln √((aU_k)² + (2/3)(f_k/Z)Ω_p²) を a について解く (左辺は a の増加関数)
        double target = Z * Math.Log(meanExcitationEv);
        double Lhs(double a)
        {
            double s = F[0] > 0 ? F[0] * Math.Log(W[0]) : 0;
            for (int k = 1; k < n; k++) s += F[k] * 0.5 * Math.Log(a * a * U[k] * U[k] + 2.0 / 3.0 * F[k] / Z * OmegaP * OmegaP);
            return s;
        }
        double lo = 1e-6, hi = 1e6;
        if (n == 1 || Lhs(lo) > target || Lhs(hi) < target)
            throw new InvalidOperationException($"InelasticGos: Sternheimer の係数が解けない (I = {meanExcitationEv:g4} eV)。I と殻の表が両立しない");
        for (int it = 0; it < 200; it++)
        {
            double mid = Math.Sqrt(lo * hi);
            if (Lhs(mid) < target) lo = mid; else hi = mid;
        }
        SternheimerA = Math.Sqrt(lo * hi);
        for (int k = 1; k < n; k++)
            W[k] = Math.Sqrt(SternheimerA * SternheimerA * U[k] * U[k] + 2.0 / 3.0 * F[k] / Z * OmegaP * OmegaP);
    }

    /// <summary>密度効果の補正 δ_F (式 3.70–3.72)。1−β² が F(0) 以上なら 0</summary>
    internal double DensityEffect(double beta2)
    {
        double rhs = 1 - beta2, op2 = OmegaP * OmegaP;
        double Fl(double l2) { double s = 0; for (int k = 0; k < Count; k++) if (W[k] > 0) s += F[k] / (W[k] * W[k] + l2); return op2 / Z * s; }
        if (!(Fl(0) > rhs)) return 0;
        double lo = 0, hi = op2;
        while (Fl(hi) > rhs) hi *= 4;
        for (int it = 0; it < 200 && hi - lo > 1e-13 * hi; it++) { double mid = 0.5 * (lo + hi); if (Fl(mid) > rhs) lo = mid; else hi = mid; }
        double L2 = 0.5 * (lo + hi), d = 0;
        for (int k = 0; k < Count; k++) if (W[k] > 0) d += F[k] * Math.Log(1 + L2 / (W[k] * W[k]));
        return d / Z - L2 * rhs / op2;
    }

    /// <summary>最小の反跳エネルギー Q₋ (式 A.31) [eV]。P−P' = W(2E−W+2mc²)/(P+P') で桁落ちを避ける</summary>
    internal static double QMinus(double e, double w)
    {
        double p = Math.Sqrt(e * (e + 2 * Mc2)), pp = Math.Sqrt(Math.Max(0, (e - w) * (e - w + 2 * Mc2)));
        double d = w * (2 * e - w + 2 * Mc2) / (p + pp);
        return d * d / (Math.Sqrt(Mc2 * Mc2 + d * d) + Mc2);
    }

    /// <summary>損失 W・反跳エネルギー Q から入射電子の散乱角の余弦 (運動量の保存)</summary>
    internal static double CosTheta(double e, double w, double q)
    {
        double p2 = e * (e + 2 * Mc2), pp2 = (e - w) * (e - w + 2 * Mc2), q2 = q * (q + 2 * Mc2);
        if (!(pp2 > 0)) return 1;
        return Math.Clamp((p2 + pp2 - q2) / (2 * Math.Sqrt(p2 * pp2)), -1.0, 1.0);
    }

    /// <summary>
    /// 振動子 k の遠い衝突のエネルギー e での量: 修正後の共鳴エネルギー W'_k と打ち切り Q'_k (式 3.78・3.80)、
    /// 損失の区間 [U_k, B] (B = min(3W'_k − 2U_k, (e+U_k)/2)) と、切った三角分布の 0 次・1 次モーメント M0・M1 (切らなければ 1・W'_k)。
    /// 価電子は δ(W − W_cb) で、W_cb ≤ e/2 のときだけ開く。閉じていれば false。
    /// </summary>
    internal bool Distant(int k, double e, out double wk, out double qk, out double bMax, out double m0, out double m1)
    {
        wk = W[k]; qk = k == 0 ? W[0] : U[k]; bMax = wk; m0 = 1; m1 = wk;
        if (!(F[k] > 0) || !(wk > 0)) return false;
        if (k == 0) return wk <= 0.5 * e;
        if (!(e > U[k])) return false;
        double wdis = 3 * W[k] - 2 * U[k];
        if (!(e > wdis))
        {
            wk = (e + 2 * U[k]) / 3;
            qk = U[k] * e / wdis;
        }
        double d = 3 * wk - 2 * U[k], u = U[k], b = Math.Min(d, 0.5 * (e + u));
        if (!(b > u)) return false;
        double dd = d - u;
        bMax = b;
        m0 = 1 - (d - b) * (d - b) / (dd * dd);
        m1 = 2 / (dd * dd) * (d * (b * b - u * u) / 2 - (b * b * b - u * u * u) / 3);
        return true;
    }

    /// <summary>近い衝突の損失の区間 [W_min, W_max] (束縛殻は [U_k, (e+U_k)/2]、価電子は [W_cb, e/2])。閉じていれば false</summary>
    internal bool Close(int k, double e, out double wMin, out double wMax)
    {
        wMin = k == 0 ? W[0] : U[k];
        wMax = 0.5 * (e + U[k]);
        return F[k] > 0 && wMin > 0 && wMin < wMax && (k == 0 || e > U[k]);
    }

    /// <summary>近い衝突の損失分布の不定積分 (式 3.86–3.87)。E' = e + U_k は損失の比と W²/E'² に、a = (e/(e+mc²))² は入射のエネルギーで (PenRed と同じ)。
    /// G0 = ∫F/W² dW、G1 = ∫F/W dW</summary>
    private static (double G0, double G1) CloseIntegrals(double ep, double a, double w)
    {
        double r = ep - w;
        double g0 = -1 / w + 1 / r + (a - 1) / ep * Math.Log(w / r) + a * w / (ep * ep);
        double g1 = Math.Log(w) + ep / r + (2 - a) * Math.Log(r) + a * w * w / (2 * ep * ep);
        return (g0, g1);
    }

    /// <summary>Møller の係数 a = (e/(e+mc²))²</summary>
    private static double MollerA(double e) { double t = e / (e + Mc2); return t * t; }

    /// <summary>
    /// エネルギー e [eV] での経路ごとの逆平均自由行程 [1/nm] と阻止能 [eV/nm] (規格化前)。配列の長さは 3·Count、添字 = 3k + (int)Channel。
    /// </summary>
    internal void Rates(double e, double[] rate, double[] stopping)
    {
        double g = 1 + e / Mc2, beta2 = 1 - 1 / (g * g);
        double c = 2 * Math.PI * E2 * E2 / (beta2 * Mc2) * AtomsPerNm3; // [eV·nm²]·[1/nm³] = eV/nm
        double trans = Math.Max(0, Math.Log(1 / (1 - beta2)) - beta2 - DensityEffect(beta2));
        double am = MollerA(e);
        for (int k = 0; k < Count; k++)
        {
            int i = 3 * k;
            rate[i] = rate[i + 1] = rate[i + 2] = stopping[i] = stopping[i + 1] = stopping[i + 2] = 0;
            double cf = c * F[k];
            if (Distant(k, e, out double wk, out double qk, out _, out double m0, out double m1))
            {
                double qm = QMinus(e, wk);
                if (qm < qk)
                {
                    double lnl = Math.Log(qk * (qm + 2 * Mc2) / (qm * (qk + 2 * Mc2)));
                    rate[i] = cf / wk * lnl * m0; stopping[i] = cf / wk * lnl * m1;
                    rate[i + 1] = cf / wk * trans * m0; stopping[i + 1] = cf / wk * trans * m1;
                }
            }
            if (Close(k, e, out double w1, out double w2))
            {
                double ep = e + U[k];
                var (a0, a1) = CloseIntegrals(ep, am, w1);
                var (b0, b1) = CloseIntegrals(ep, am, w2);
                rate[i + 2] = cf * (b0 - a0); stopping[i + 2] = cf * (b1 - a1);
            }
        }
    }

    /// <summary>経路 (振動子 k・種類) を決めたあとの損失 W [eV] と散乱角の余弦を引く (エネルギー e [eV]。率を出したのと同じ e を渡すこと)</summary>
    internal (double W, double CosTheta) Sample(int k, Channel ch, double e, Random rnd)
    {
        if (ch == Channel.Close)
        {
            if (!Close(k, e, out double w1, out double w2)) return (0, 1);
            double ep = e + U[k], a = MollerA(e);
            double fmax = 1 + 1.25 * a; // r = W/(E'−W) ∈ [0,1] で r²−r ≤ 0、r + (W/E')² ≤ 1.25
            for (int it = 0; it < 100000; it++)
            {
                double w = w1 * w2 / (w2 - rnd.NextDouble() * (w2 - w1)); // ∝ 1/W² で引く
                double r = w / (ep - w);
                double fm = 1 + r * r - r + a * (r + w * w / (ep * ep)); // 式 3.87
                if (rnd.NextDouble() * fmax <= fm)
                    return (w, CosTheta(e, w, w)); // 近い衝突は Q = W (二体の極限)
            }
            throw new InvalidOperationException("InelasticGos: 近い衝突の棄却が尽きた");
        }
        if (!Distant(k, e, out double wk, out double qk, out double b, out _, out _)) return (0, 1);
        // 損失: 価電子 = δ(W−W_cb)、束縛殻 = 三角分布 p ∝ (D − W) を [U_k, B] で切ったもの。x = D − W の密度は ∝ x なので x² を一様に
        double wloss = wk;
        if (k > 0)
        {
            double d = 3 * wk - 2 * U[k], x0 = d - b, x1 = d - U[k];
            wloss = d - Math.Sqrt(x0 * x0 + rnd.NextDouble() * (x1 * x1 - x0 * x0));
        }
        if (ch == Channel.DistantTransverse) return (wloss, 1); // 横の相互作用では曲げない (§3.2.2)
        // 縦: Q を ∝ 1/(Q(Q+2mc²)) で [Q₋, Q'_k] から引く。角度は修正後の共鳴エネルギー W'_k で決める (式 3.83)
        double qm = QMinus(e, wk);
        if (!(qm < qk)) return (wloss, 1);
        double y0 = Math.Log(qm / (qm + 2 * Mc2)), y1 = Math.Log(qk / (qk + 2 * Mc2));
        double ey = Math.Exp(y0 + rnd.NextDouble() * (y1 - y0));
        double q = 2 * Mc2 * ey / (1 - ey);
        return (wloss, CosTheta(e, wk, q));
    }
}
