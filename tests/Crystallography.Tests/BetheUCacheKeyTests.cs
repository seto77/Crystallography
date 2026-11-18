using System.Drawing;
using System.Numerics;
using Xunit;

namespace Crystallography.Tests;

/// <summary>261118Cl 新規: getU のキャッシュ (uDictionary) の鍵 (codex2 第 2 回の単位 3)。
/// 以前の鍵は index の組だけで、電圧・積分角度・模型を含まず、呼び出し側の Clear の時機に正しさを頼っていた
/// (検出器範囲の U′ が全立体角の値として残った不具合 A1、前の電圧の U0 を使った不具合 2f804a1・RunAlchemi)。
/// ここでは、別の条件で読んだ後でも、汚れていない新しいインスタンスと同じ値が返ることを確かめる (Clear は呼ばない)。</summary>
public class BetheUCacheKeyTests
{
    private static readonly BetheMethod.Beam G = new((2, 2, 0), new Vector3DBase(2 / 0.54309, 2 / 0.54309, 0));

    private static BetheMethod Fresh() => new(Silicon(0.5));

    [Fact]
    public void Voltage_IsPartOfTheKey()
    {
        var bethe = Fresh();
        var at200 = bethe.getU(200, G);
        var at300 = bethe.getU(300, G); // 同じ index を別の電圧で (Clear しない)
        Assert.Equal(Fresh().getU(300, G), at300);
        Assert.NotEqual(at200.Real, at300.Real); // 相対論の補正で電圧により変わる
    }

    [Fact]
    public void AnnularRange_IsPartOfTheKey()
    {
        var bethe = Fresh();
        var annular = bethe.getU(200, G, null, 0.04, 0.12);
        var full = bethe.getU(200, G); // 環状の後に全立体角
        Assert.Equal(Fresh().getU(200, G), full);
        Assert.NotEqual(annular.Imag, full.Imag);

        var bethe2 = Fresh();
        bethe2.getU(200, G); // 全立体角の後に
        var annular2 = bethe2.getU(200, G.Index, G.Vec, new Vector3DBase(0, 0, 0), new Vector3DBase(0, 0, 0), 0.04, 0.12); // STEM 用の環状の入口
        Assert.Equal(Fresh().getU(200, G, null, 0.04, 0.12), annular2);
    }

    [Fact]
    public void AllocationFreeOverload_ReturnsTheSameAsTheFullOverload()
    {
        var bethe = Fresh();
        var viaBeam = bethe.getU(200, G);
        var viaIndex = bethe.getU(200, G.Index, G.Vec); // キャッシュにある値を割り当てなしで返す入口
        Assert.Equal(viaBeam, viaIndex);
        var other = Fresh();
        Assert.Equal(viaBeam, other.getU(200, G.Index, G.Vec)); // キャッシュが空のときは本体へ委ねる
    }

    private static Crystal Silicon(double bisoAng2)
    {
        int series = 0;
        for (var i = 1; i < SymmetryStatic.TotalSpaceGroupNumber; i++)
            if (SymmetryStatic.Symmetries[i].SpaceGroupNumber == 227) { series = i; break; }
        var cell = (0.54309, 0.54309, 0.54309, Math.PI / 2, Math.PI / 2, Math.PI / 2);
        Atoms[] atoms = [new("Si1", 14, 0, 0, null, series, new Vector3DBase(0, 0, 0), 1,
            new DiffuseScatteringFactor(DiffuseScatteringFactor.Type.B, true, bisoAng2 / 100.0, 0, new double[6], new double[6], cell))];
        return new Crystal(cell, null, series, "Si (BetheUCacheKeyTests)", Color.Gray, new Matrix3D(),
            atoms, ("", "", "", ""), Bonds.GetVestaBonds(atoms.Select(x => x.AtomicNumber)));
    }
}
