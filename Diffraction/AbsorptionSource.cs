// 260811Cl 新規作成: 吸収 (複素ポテンシャルの虚部) の**出所を型で持つ**。
//
// 動機は予防であって、現状の不具合の修正ではない。2026-08-11 の文献精査 (ICSC 2003 の照合) で
// **現在の吸収経路は TDS だけであり、二重計上は無い**ことを確認している。危ないのはこの先で:
//
//   ALCHEMI の非チャネリング項 Y_dech = (μ00/V_c)(t − L_coh) は、
//   「**コヒーレント場から失われた電子は全部まだ試料内にいて、ランダム方位でイオン化を続ける**」
//   という仮定でできている。t − L_coh は虚部が作る減衰**すべて**を拾うので、
//   将来ここへ mean absorption / 真の非弾性損失 / 経験的 damping を混ぜると、
//   **それらまでイオン化電子として再注入されて破綻する** (静かに、しかも符号は必ず過大の側)。
//
// そこで「虚部に何が入っているか」を <see cref="AbsorptionSource"/> で宣言させ、
// 再注入する側 (AlchemiReduction.Yield) は **TDS だけのときしか走らない**ようにした。
// 型が無ければこの前提はコメントの中にしか無く、次に虚部を触る人には見えない。
//
// ⚠ 虚部の作り方を変えたら <see cref="BetheMethod.ImaginaryPotentialAbsorption"/> を必ず更新すること。
//   更新を忘れても ALCHEMI の数値が静かにずれることはない (Yield が例外で落ちる) が、
//   逆に「TDS のままだ」と嘘を書けば落ちない。ここは人間の規律に頼る唯一の点。

using System;

namespace Crystallography;

/// <summary>260811Cl 追加: 複素ポテンシャルの虚部に入っている減衰の**出所**。
/// 「どれだけ減衰するか」ではなく「その減衰した電子がその後どうなるか」の分類である。</summary>
[Flags]
public enum AbsorptionSource
{
    /// <summary>吸収なし (虚部が恒等的に 0)。</summary>
    None = 0,

    /// <summary>熱散漫散乱 (TDS)。電子はほぼ全エネルギーを保ったまま試料内に残り、方向だけを失う。
    /// **非チャネリング電子として再注入してよい唯一の出所**。
    /// 現行 ReciPro の虚部はここだけ (<see cref="AtomStatic.ES.FactorImaginary(double, double, double)"/> 系)。</summary>
    TdsRedistributable = 1,

    /// <summary>真の非弾性損失 (プラズモン励起・内殻イオン化に伴う損失など)。
    /// 電子はエネルギーを失う = その後のイオン化断面積が入射時と同じではないので、
    /// **弾性チャネルと同じ μ で再注入してはいけない**。</summary>
    TrueLoss = 2,

    /// <summary>経験的な damping / 実測に合わせた mean absorption。物理的な出所が混ざっているため、
    /// **どの割合を再注入してよいか原理的に決められない**。</summary>
    Phenomenological = 4,
}

/// <summary>260811Cl 追加: <see cref="AbsorptionSource"/> の判定をまとめた拡張。</summary>
public static class AbsorptionSourceExtensions
{
    /// <summary>失われた流束を**まるごと**「試料内に残る非チャネリング電子」として再注入してよいか。
    /// TDS だけのとき (および吸収が無いとき) に限り true。
    /// ⚠ <c>HasFlag(TdsRedistributable)</c> ではない — TDS **以外が混ざっていない**ことが条件。</summary>
    public static bool IsFullyRedistributable(this AbsorptionSource sources)
        => (sources & ~AbsorptionSource.TdsRedistributable) == 0;
}
