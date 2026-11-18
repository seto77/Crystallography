using Xunit;

namespace Crystallography.Tests;

/// <summary>261118Cl 新規 (後方源の検証 wf_26a4f176-803 の批判の段の検算 hook を公開の回帰に): ブロック共有 (EbsdRepresentativeDirectionBlockSize &gt; 1) の画素は
/// 「代表のビームの組 (代表の K0 で選んだもの) + 自分の K0」で解かれ、後方源の軸の直し (<see cref="BetheMethod.MasterBackSourceAxisCorrection"/>) では
/// 後方源も自分の入射の軸で組まれること。ブロック 1 で研究用の口から代表の組を返して解いた値と、製品のブロック 3 の値が丸めの範囲で一致すれば、
/// 割り当ての再現と画素ごとの軸の両方が確かめられる (代表の軸で組んでいれば、画素と代表の向きの差の分だけ、直しのときに差が出る)。
/// 直しの効き (ブロック 1 の on 対 off) が丸めより十分大きいことも確かめて、試験が軸の違いを見分けられることを保証する</summary>
[Collection("Bethe global settings")]
public class EbsdBackSourceBlockAxisTests
{
    private static double[][] Run(Crystal crystal, Vector3DBase[] dirs, int block, bool correction, Func<BetheMethod, Vector3DBase, Vector3DBase, int, BetheMethod.Beam[]> hook)
    {
        bool savedAxis = BetheMethod.MasterBackSourceAxisCorrection; var savedHook = BetheMethod.EbsdBeamSelectionResearchHook; // 261118Cl (codex2 第 8 回): 元の値に戻す (既定が変わっても他の試験を壊さない)
        BetheMethod.MasterBackSourceAxisCorrection = correction;
        BetheMethod.EbsdBeamSelectionResearchHook = hook;
        try
        {
            var bethe = new BetheMethod(crystal) { EbsdRepresentativeDirectionBlockSize = block, UseLocalSurfacePerBeamDirection = false, EbsdUseMasterPatternSymmetry = false };
            using var done = new ManualResetEventSlim();
            Exception error = null;
            bethe.EBSD_Completed += (_, e) => { error = e.Error; done.Set(); };
            bethe.RunEBSDNew(32, [20.0], Matrix3D.IdentityMatrix, [10.0, 30.0], dirs, BetheMethod.Solver.Eigen_Eigen, 4, false, true, false);
            Assert.True(done.Wait(TimeSpan.FromMinutes(10)), "EBSD did not finish");
            Assert.Null(error);
            Assert.Equal(correction ? BetheMethod.BackSourceAxisIncident : BetheMethod.BackSourceAxisCrystalZ, bethe.EbsdBackSourceAxisUsed);
            // Assert.Equal(BetheMethod.BeamCandidateCut.ToString(), bethe.EbsdBeamCandidateCutUsed); // 261118Cl 変更前 (codex2 第 15 回: 研究用の口では設定値を書かない)
            Assert.Equal(hook != null ? BetheMethod.BeamCandidateCutResearchHook : BetheMethod.BeamCandidateCut.ToString(), bethe.EbsdBeamCandidateCutUsed); // 261118Cl 追加: 実際に通った分岐が記録される
            return [.. Enumerable.Range(0, 2).Select(t => bethe.Disks[0][t].Amplitudes.Take(dirs.Length).Select(a => a.Real * a.Real + a.Imaginary * a.Imaginary).ToArray())];
        }
        finally { BetheMethod.MasterBackSourceAxisCorrection = savedAxis; BetheMethod.EbsdBeamSelectionResearchHook = savedHook; }
    }

    private static double MaxRel(double[][] a, double[][] b, IEnumerable<int> idx) =>
        idx.SelectMany(i => Enumerable.Range(0, a.Length).Where(t => b[t][i] > 0).Select(t => Math.Abs(a[t][i] / b[t][i] - 1))).DefaultIfEmpty(double.NaN).Max();

    [Fact]
    public void BlockSharedPixels_UseTheirOwnIncidentAxis()
    {
        var crystal = EbsdBackSourceAxisTests.GaAs(0.6); // B ≠ 0 (B = 0 では後方 TDS の源が 0 で、軸の違いが見えない)
        const int g = 12, bs = 3;
        Vector3DBase[] dirs = [.. MasterPattern.CreateDirectionsSquare(g, MasterPattern.Hemisphere.PositiveZ), .. MasterPattern.CreateDirectionsSquare(g, MasterPattern.Hemisphere.NegativeZ)];
        int n = dirs.Length, hl = g * g;
        int Nearest(Vector3DBase s) // 口に渡る表面 (= −方向) から方向の番号を引く
        {
            double l = s.Length, best = -2; int arg = -1;
            for (int i = 0; i < n; i++) { double d = -(s.X * dirs[i].X + s.Y * dirs[i].Y + s.Z * dirs[i].Z) / l; if (d > best) { best = d; arg = i; } }
            return arg;
        }
        // ブロック 1・従来の構築で、方向ごとに製品が渡す K0 と表面を記録する (口の中は製品と同じ Find_gVectors)
        var k0Rec = new Vector3DBase[n]; var sRec = new Vector3DBase[n];
        var r1Legacy = Run(crystal, dirs, 1, false, (b, k0, s, nmax) => { int p = Nearest(s); lock (k0Rec) { k0Rec[p] = k0; sRec[p] = s; } return b.Find_gVectors(Matrix3D.IdentityMatrix, k0, s, nmax); });
        Assert.All(k0Rec, k => Assert.NotNull(k));
        var r1Correction = Run(crystal, dirs, 1, true, null);
        // 製品のブロック 3 と同じ割り当て (BetheMethod の ebsdNew_DoWork の代表の式)
        int Rep(int p) { int off = p / hl * hl, local = p - off, w = local % g, h = local / g; return off + Math.Min(w / bs * bs + bs / 2, g - 1) + Math.Min(h / bs * bs + bs / 2, g - 1) * g; }
        var nonCenters = Enumerable.Range(0, n).Where(p => Rep(p) != p).ToArray();
        Assert.NotEmpty(nonCenters);
        // 試験の判別力: 直しの効きが丸めより十分大きい
        Assert.True(MaxRel(r1Correction, r1Legacy, Enumerable.Range(0, n)) > 1e-5, "the axis correction must change the intensity measurably");
        foreach (var correction in new[] { false, true })
        {
            var product = Run(crystal, dirs, bs, correction, null);
            var reference = Run(crystal, dirs, 1, correction, (b, k0, s, nmax) => { int r = Rep(Nearest(s)); return b.Find_gVectors(Matrix3D.IdentityMatrix, k0Rec[r], sRec[r], nmax); });
            double rel = MaxRel(product, reference, nonCenters);
            Assert.True(rel < 1e-10, $"correction={correction}: block {bs} vs block 1 with the representative beams, non-center pixels max rel {rel:e2}");
        }
    }

    /// <summary>261118Cl 追加 (codex2 第 8 回 §5): 後方源の軸の方式の記録が、マスターパターンのファイルの保存と読込で往復し、
    /// 記録の無い旧いヘッダは "" (方式の記録なし) として読めること</summary>
    [Fact]
    public void BackSourceAxisRecord_RoundTripsThroughTheFile()
    {
        int g = 4;
        float[][] Planes() => [Enumerable.Range(0, g * g).Select(i => (float)i).ToArray()];
        var mp = new MasterPattern(g, [20.0], [10.0], Planes(), Planes(), MasterPattern.Types.Square) { BackSourceAxis = BetheMethod.BackSourceAxisIncident, BeamCandidateCut = nameof(BeamCandidateCutMode.Shell) }; // 261118Cl 打ち切り方も
        var path = Path.Combine(Path.GetTempPath(), $"rpmp_axis_{Guid.NewGuid():N}.rpmp");
        try
        {
            MasterPatternFile.Save(mp, path, new MasterPatternFileInfo());
            var loaded = MasterPatternFile.Load(path, loadElectrons: false);
            Assert.Equal(BetheMethod.BackSourceAxisIncident, loaded.Info.BackSourceAxis);
            Assert.Equal(BetheMethod.BackSourceAxisIncident, loaded.MasterPattern.BackSourceAxis);
            Assert.Equal(nameof(BeamCandidateCutMode.Shell), loaded.Info.BeamCandidateCut); // 261118Cl 追加
            Assert.Equal(nameof(BeamCandidateCutMode.Shell), loaded.MasterPattern.BeamCandidateCut);
        }
        finally { File.Delete(path); }
        var old = System.Text.Json.JsonSerializer.Deserialize<MasterPatternFileInfo>("{}");
        Assert.Equal("", old.BackSourceAxis);
        Assert.Equal("", old.BeamCandidateCut); // 261118Cl 追加: 古いファイル = 従来の打ち切り方
        // 261118Cl 追加 (codex2 第 16 回): 打ち切り方を求めたら、記録の無いファイルは一致としない。空の要求は引数の誤り
        var si = EbsdBackSourceAxisTests.GaAs(0.6);
        var info = new MasterPatternFileInfo { CrystalFingerprint = MasterPatternFile.CrystalFingerprint(si), BeamEnergyKeV = 20 };
        Assert.Null(MasterPatternFile.CheckCompatibility(info, si, 20));
        Assert.Contains("does not record", MasterPatternFile.CheckCompatibility(info, si, 20, beamCandidateCut: "Legacy"));
        info.BeamCandidateCut = "Shell";
        Assert.Contains("Shell", MasterPatternFile.CheckCompatibility(info, si, 20, beamCandidateCut: "Legacy"));
        Assert.Null(MasterPatternFile.CheckCompatibility(info, si, 20, beamCandidateCut: "Shell"));
        Assert.Throws<ArgumentException>(() => MasterPatternFile.CheckCompatibility(info, si, 20, beamCandidateCut: ""));
    }
}
