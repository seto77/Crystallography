using Xunit;

namespace Crystallography.Tests;

/// <summary>261111Cl 新規: Find_gVectors の境の同値の組 (評価値の小さい順に N+1 本を選び、末尾の同点の組を丸ごと扱う。
/// 261112Cl: 既定は丸ごと採る (IncludeGroup、作者の判断)。それまでは丸ごと捨てる (DropGroup))。
/// 共有の QuickSelect の同値の扱いが変わったとき (Crystallography 6b0d049)、STEM-EDX のビームの並びが変わった。
/// 晶帯軸 [001] の Si で、選んだ集合が 4mm (入射方向を保つ点群の操作) で閉じている (同値の組を途中で切らない) ことと、
/// 並びが (評価値, 指数) の順で決まる (候補の並びに依らない) ことを、多くの N で確かめる</summary>
public class BeamBoundaryTieTests
{
    [Theory]
    [InlineData(200.0)]
    [InlineData(20.0)]
    public void FindGVectors_ZoneAxis_ClosedUnder4mm_AndOrderDeterministic(double kv)
    {
        var crystal = ConvertCrystalData.ConvertFromCifText(StemDepthPhaseTests.SiCifText).ToCrystal();
        var rot = Matrix3D.Rot(new Vector3DBase(0, 0, 1), 0);
        crystal.RotationMatrix = rot;
        var bethe = new BetheMethod(crystal);
        bethe.GetDifractedBeamAmpriltudes(10, kv, rot, 10); // 電圧 (getU が使う) を設定する
        var vecK0 = bethe.getVecK0(UniversalConstants.Convert.EnergyToElectronWaveNumber(kv), 0);
        var broken = new List<string>();
        var unordered = new List<string>();
        for (int n = 8; n <= 700; n += 3)
        {
            var beams = bethe.Find_gVectors(rot, vecK0, n);
            var set = beams.Select(b => b.Index).ToHashSet();
            foreach (var (h, k, l) in set)
            {
                var missing = new[] { (k, -h, l), (-h, -k, l), (-k, h, l), (h, -k, l), (-h, k, l), (k, h, l), (-k, -h, l) }.Where(e => !set.Contains(e)).ToArray();
                if (missing.Length > 0) { broken.Add($"N{n}: ({h},{k},{l}) はあるが {string.Join(" ", missing)} が無い"); break; }
            }
            for (int i = 1; i < beams.Length; i++)
            {
                int c = beams[i - 1].Rating.CompareTo(beams[i].Rating);
                if (c > 0 || (c == 0 && beams[i - 1].Index.CompareTo(beams[i].Index) >= 0)) { unordered.Add($"N{n}: {i - 1} と {i}"); break; }
            }
        }
        Assert.True(broken.Count == 0, $"{kv} kV: 対称で閉じていない N が {broken.Count} 個\n{string.Join("\n", broken.Take(10))}");
        Assert.True(unordered.Count == 0, $"{kv} kV: (評価値, 指数) の順でない N が {unordered.Count} 個\n{string.Join("\n", unordered.Take(10))}");
    }
}
