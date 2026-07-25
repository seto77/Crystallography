#region using
using System;
using System.Threading;
#endregion

namespace Crystallography;

/// <summary>
/// 実測 EBSD パターンの探索・較正を 1 回走らせる間の状態 (実行中フラグ・世代番号・キャンセル要求)。260726Cl 追加。
/// FormEBSD の指数付け UI から失効ロジックだけを抜き出したもので、WinForms へは依存しないので単体で検証できる (正本 §6 P1)。
///
/// 解く問題は 2 つ:
///   ①実行中に入力 (実測画像・検出器幾何・結晶・波長・MasterPattern など) が変わったら、走っている計算を止める。
///   ②止まりきる前に完了した結果を、変わったあとの入力へ適用しない。
/// ②のために <see cref="TryBegin"/> が返す世代番号を呼び出し側が保持し、await から戻ったところで
/// <see cref="IsCurrent"/> に問い合わせる。<see cref="Invalidate"/> が呼ばれていれば false になる。
///
/// すべて UI スレッドから呼ばれる前提 (トークンだけがワーカースレッドから読まれる)。
/// </summary>
public sealed class EbsdIndexingSession : IDisposable
{
    CancellationTokenSource cts;

    /// <summary>入力が変わるたびに進む。await を跨いだ結果が古くないかの判定に使う</summary>
    public int Generation { get; private set; }

    /// <summary>探索・較正を実行中か (ボタンの相互排他に使う)</summary>
    public bool Busy { get; private set; }

    /// <summary>実行中でなければ開始する。二重起動を防ぐため、実行中は false を返して何もしない</summary>
    /// <param name="token">この実行のキャンセルトークン。Invalidate / Cancel / Dispose で発火する</param>
    /// <param name="generation">await から戻ったあとに IsCurrent へ渡す世代番号</param>
    public bool TryBegin(out CancellationToken token, out int generation)
    {
        if (Busy) { token = CancellationToken.None; generation = Generation; return false; }
        Busy = true;
        cts?.Dispose(); //例外的な経路で End を通らなかった場合でも古い CTS を残さない
        cts = new CancellationTokenSource();
        token = cts.Token;
        generation = Generation;
        return true;
    }

    /// <summary>実行を終える (成功・中止・失敗のいずれでも必ず呼ぶ)</summary>
    public void End()
    {
        cts?.Dispose();
        cts = null;
        Busy = false;
    }

    /// <summary>
    /// 入力が変わったので、実行中の計算を止め、これまでの結果を失効させる。
    /// 実行中でなくても世代は進む (直前に完了した結果も古くなるため)。
    /// </summary>
    public void Invalidate()
    {
        Generation++;
        Cancel();
    }

    /// <summary>実行中ならキャンセルを要求する (世代は進めない)。フォームの Dispose など、結果を捨てるだけの場面で使う</summary>
    public void Cancel()
    {
        //End 済みなら cts は null。Cancel 済みの CTS へ再度 Cancel しても副作用はない
        if (cts is { IsCancellationRequested: false }) cts.Cancel();
    }

    /// <summary>TryBegin で受け取った世代番号が今も有効か (false なら結果を捨てる)</summary>
    public bool IsCurrent(int generation) => generation == Generation;

    public void Dispose()
    {
        Cancel();
        cts?.Dispose();
        cts = null;
        Busy = false;
    }
}
