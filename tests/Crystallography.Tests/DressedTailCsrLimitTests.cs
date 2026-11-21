using Xunit;

namespace Crystallography.Tests;

/// <summary>261121Cl 新規 (codex2 第 30 回 §5): 着衣の裾の疎な表 (CSR) の行の先頭は int。旧は累積も int で、非零が 2³¹ を超えると黙ってあふれた
/// (opx 3° の R 60 では非零 3.48 × 10⁹ と数えた。今は候補の上限 15 万が先に効くので届かない)。累積を long で数え、配列の上限
/// (Array.MaxLength) を超えたら資源の上限の例外にする。境の上下は合成の行の長さで試す (巨大な表は確保しない)</summary>
public class DressedTailCsrLimitTests
{
    [Fact]
    public void RowPointers_AreCumulativeLengths()
    {
        int[] lengths = [3, 0, 5, 1];
        var ptr = BetheMethod.CsrRowPointersChecked(lengths.Length, i => lengths[i], "test table");
        Assert.Equal([0, 3, 3, 8, 9], ptr);
    }

    [Fact]
    public void RowPointers_AcceptExactlyTheArrayLimit()
    {
        int[] lengths = [Array.MaxLength - 10, 10];
        var ptr = BetheMethod.CsrRowPointersChecked(lengths.Length, i => lengths[i], "test table");
        Assert.Equal(Array.MaxLength, ptr[^1]);
    }

    [Fact]
    public void RowPointers_ThrowOneAboveTheArrayLimit()
    {
        int[] lengths = [Array.MaxLength - 10, 11];
        var ex = Assert.Throws<BetheMethod.DressedTailTableLimitException>(() => BetheMethod.CsrRowPointersChecked(lengths.Length, i => lengths[i], "test table"));
        Assert.Equal((long)Array.MaxLength + 1, ex.Nonzeros);
        Assert.Contains("test table", ex.Message);
    }

    [Fact]
    public void RowPointers_ThrowInsteadOfWrappingPastInt32()
    {
        // 旧の int の累積なら int.MaxValue + int.MaxValue = −2 に化けて、負の長さの配列で別の例外になるか、黙って壊れた表になった
        int[] lengths = [int.MaxValue / 2 + 1, int.MaxValue / 2 + 1, 7];
        var ex = Assert.Throws<BetheMethod.DressedTailTableLimitException>(() => BetheMethod.CsrRowPointersChecked(lengths.Length, i => lengths[i], "test table", "hint text"));
        Assert.True(ex.Nonzeros > Array.MaxLength);
        Assert.Contains("hint text", ex.Message);
        Assert.IsAssignableFrom<InvalidOperationException>(ex); // 本物のメモリ不足 (OutOfMemoryException) とは分ける
    }
}
