using System.Drawing;
using Xunit;

namespace Crystallography.Tests;

/// <summary>261118Cl 新規 (codex2 第 7・8 回、N64 の非対称の直しの試作): Find_gVectors の候補の打ち切り方 <see cref="BetheMethod.BeamCandidateCut"/>。
/// 従来 (Legacy) は BFS の発見の順の先頭 8N 本で切るので、Si・20 kV・N64 の一般の方向とその m-3m の 48 像のうち 12 像で境の 1 本が入れ替わる
/// (段ごとの測定 wf_n64_stage_261118: (−12 −10 12) が通過のちょうど 512 番目か 513〜519 番目か)。殻の単位 (Shell) では 48 像が同じ組を選ぶ。
/// 従来の側の非対称も確かめて、試験が違いを見分けられることを保証する</summary>
[Collection("Bethe global settings")]
public class BeamCandidateShellCutTests
{
    private static Crystal Silicon()
    {
        int series = 0;
        for (var i = 1; i < SymmetryStatic.TotalSpaceGroupNumber; i++)
            if (SymmetryStatic.Symmetries[i].SpaceGroupNumber == 227) { series = i; break; }
        var cell = (0.54309, 0.54309, 0.54309, Math.PI / 2, Math.PI / 2, Math.PI / 2);
        Atoms[] atoms = [new("Si1", 14, 0, 0, null, series, new Vector3DBase(0, 0, 0), 1,
            new DiffuseScatteringFactor(DiffuseScatteringFactor.Type.B, true, 0.5 / 100.0, 0, new double[6], new double[6], cell))];
        return new Crystal(cell, null, series, "Si (BeamCandidateShellCutTests)", Color.Gray, new Matrix3D(), atoms, ("", "", "", ""), []);
    }

    /// <summary>48 の像それぞれで選んだ組を基準の方向へ戻し、像 0 と違う像の数と像 0 の組を返す</summary>
    private static (int Differ, HashSet<(int, int, int)> Image0) Images(BeamCandidateCutMode mode, int nb, double kv = 20)
    {
        var saved = BetheMethod.BeamCandidateCut;
        try
        {
            BetheMethod.BeamCandidateCut = mode;
            var bm = new BetheMethod(Silicon());
            double l = Math.Sqrt(0.21 * 0.21 + 0.47 * 0.47 + 0.857 * 0.857);
            double[] c = [0.21 / l, 0.47 / l, 0.857 / l];
            int[][] perms = [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]];
            double kvac = UniversalConstants.Convert.EnergyToElectronWaveNumber(kv), u0 = bm.getU(kv).Real.Real;
            HashSet<(int, int, int)> image0 = null; int differ = 0;
            for (int i = 0; i < 48; i++)
            {
                var pp = perms[i / 8]; int si = i % 8; int[] sg = [(si & 1) == 0 ? 1 : -1, (si & 2) == 0 ? 1 : -1, (si & 4) == 0 ? 1 : -1];
                var d = new Vector3DBase(sg[0] * c[pp[0]], sg[1] * c[pp[1]], sg[2] * c[pp[2]]);
                var s = new Vector3DBase(-d.X, -d.Y, -d.Z);
                var beams = bm.Find_gVectors(Matrix3D.IdentityMatrix, bm.getVecK0(kvac, u0, d, s), s, nb);
                var mapped = beams.Select(b => { int[] img = [b.H, b.K, b.L]; var bs = new int[3]; for (int j = 0; j < 3; j++) bs[pp[j]] = sg[j] * img[j]; return (bs[0], bs[1], bs[2]); }).ToHashSet();
                if (image0 == null) image0 = mapped;
                else if (!mapped.SetEquals(image0)) differ++;
            }
            return (differ, image0);
        }
        finally { BetheMethod.BeamCandidateCut = saved; }
    }

    [Fact]
    public void ShellCut_SelectsTheSameSetOnAllSymmetryImages()
    {
        var legacy = Images(BeamCandidateCutMode.Legacy, 64);
        Assert.True(legacy.Differ > 0, "the legacy BFS-order cut must show the asymmetry this test is about"); // 測定では 12/47
        var shell = Images(BeamCandidateCutMode.Shell, 64);
        Assert.Equal(0, shell.Differ);
        Assert.True(shell.Image0.SetEquals(legacy.Image0), "N64: the shell cut keeps the image-0 set of the legacy cut (measured)");
        foreach (var nb in new[] { 56, 63 }) // 従来で 4/47・24/47 が入れ替わる N
            Assert.Equal(0, Images(BeamCandidateCutMode.Shell, nb).Differ);
    }

    private static (int H, int K, int L)[] Select(BetheMethod bm, double kv, int nb, Vector3DBase d)
    {
        double kvac = UniversalConstants.Convert.EnergyToElectronWaveNumber(kv), u0 = bm.getU(kv).Real.Real;
        var s = new Vector3DBase(-d.X, -d.Y, -d.Z);
        return [.. bm.Find_gVectors(Matrix3D.IdentityMatrix, bm.getVecK0(kvac, u0, d, s), s, nb).Select(b => b.Index)];
    }

    /// <summary>261118Cl 追加 (codex2 第 9 回): 同じ入力なら結果はキャッシュの履歴に依らない。大きい N の呼び出しが作った大きいキャッシュの後でも、
    /// 小さい N は新しく作ったときと同じ組を返す (従来・殻の両方。旧は大きいキャッシュの全部を使い、候補が尽きる条件で結果が変わりえた)</summary>
    [Fact]
    public void Selection_DoesNotDependOnTheCacheHistory()
    {
        var saved = BetheMethod.BeamCandidateCut;
        try
        {
            var d = new Vector3DBase(0.21, 0.47, 0.857); d = d / d.Length;
            foreach (var mode in new[] { BeamCandidateCutMode.Legacy, BeamCandidateCutMode.Shell })
            {
                BetheMethod.BeamCandidateCut = mode;
                foreach (var nb in new[] { 64, 400 })
                {
                    var fresh = Select(new BetheMethod(Silicon()), 20, nb, d);
                    var used = new BetheMethod(Silicon());
                    Select(used, 20, 1600, d); // 大きいキャッシュを先に作らせる
                    Assert.Equal(fresh, Select(used, 20, nb, d));
                }
            }
        }
        finally { BetheMethod.BeamCandidateCut = saved; }
    }

    /// <summary>261118Cl 追加 (codex2 第 9 回): 候補が N 本以下なら番兵は無いので、本物のビームを落とさない (旧は最後の 1 本を落とし、000 だけなら空を返した)。
    /// 低い電圧ではエワルドの条件を通る候補が少ない</summary>
    [Fact]
    public void FewCandidates_AreAllReturned()
    {
        foreach (var mode in new[] { BeamCandidateCutMode.Legacy, BeamCandidateCutMode.Shell })
        {
            var saved = BetheMethod.BeamCandidateCut;
            try
            {
                BetheMethod.BeamCandidateCut = mode;
                var bm = new BetheMethod(Silicon());
                var d = new Vector3DBase(0.21, 0.47, 0.857); d = d / d.Length;
                foreach (var kv in new[] { 0.02, 0.05, 0.2 })
                {
                    var beams = Select(bm, kv, 400, d);
                    var collected = BetheMethod.LastBeamCandidateScan.Collected;
                    Assert.True(collected <= 400, $"{kv} kV: the test needs fewer candidates than N ({collected})");
                    Assert.Contains((0, 0, 0), beams);
                    Assert.Equal(collected, beams.Length); // 重複除去が無ければ全部
                }
            }
            finally { BetheMethod.BeamCandidateCut = saved; }
        }
    }

    /// <summary>261118Cl 追加 (codex2 第 10 回): 殻の外端 (殻の最初の |g| × (1 + 1e-9)) を覆い切った半径と殻ごとに 1 回だけ比べる。
    /// 同じ丸めの幅の 2 点 r1 &lt; r2 が R/(1 + 1e-9) の両側にあっても、殻を丸ごと採る (旧の点ごとの比較は r1 だけを採った)。外端が R を越える殻は丸ごと外す</summary>
    [Fact]
    public void SelectShellCut_ComparesTheCoverRadiusOncePerShell()
    {
        const double R = 10, eps = 1e-9, delta = 1e-12;
        double edge = R / (1 + eps);
        double[] lengths = [0, 1, 1 * (1 + 1e-10), 5, edge - delta, edge + delta, 12];
        Assert.True((edge - delta) * (1 + eps) < R && (edge + delta) * (1 + eps) > R, "the two points must straddle R/(1+eps)");
        Assert.Equal([0, 1, 2, 3, 4, 5], Run(lengths, R, 100));
        double[] crossing = [0, 1, edge + delta, edge + 2 * delta, 12]; // 外端が R を越える殻は丸ごと外す
        Assert.Equal([0, 1], Run(crossing, R, 100));
    }

    /// <summary>261118Cl 追加 (codex2 第 10 回): 上限の本数目を含む殻は、その殻の固定した外端まで採り、次の殻で止める (上限の点の |g| から測らない)</summary>
    [Fact]
    public void SelectShellCut_TakesTheLimitShellToItsFixedEnd()
    {
        double[] lengths = [0, 1, 1 * (1 + 0.9e-9), 1 * (1 + 1.5e-9), 2];
        Assert.Equal([0, 1, 2], Run(lengths, 100, 3)); // 1·(1 + 1.5e-9) は外端 1·(1 + 1e-9) の外 = 次の殻
        Assert.Equal([0, 1, 2, 3], Run(lengths, 100, 4));
    }

    private static int[] Run(double[] lengths, double cover, int limit)
    {
        var order = Enumerable.Range(0, lengths.Length).OrderBy(i => lengths[i]).ToArray();
        var picked = new List<int>();
        BetheMethod.SelectShellCut(order, i => lengths[i], _ => true, cover, limit, i => { picked.Add(i); return picked.Count; });
        return [.. picked];
    }

    /// <summary>261118Cl 追加 (codex2 第 10 回): エワルドの条件を通る候補が 0 本 (表面の向きを逆にして 000 も通らない) なら、従来・殻とも例外でなく空を返す
    /// (第 9 回の前は範囲外の例外)</summary>
    [Fact]
    public void NoCandidate_ReturnsEmpty()
    {
        foreach (var mode in new[] { BeamCandidateCutMode.Legacy, BeamCandidateCutMode.Shell })
        {
            var saved = BetheMethod.BeamCandidateCut;
            try
            {
                BetheMethod.BeamCandidateCut = mode;
                var bm = new BetheMethod(Silicon());
                var d = new Vector3DBase(0.21, 0.47, 0.857); d = d / d.Length;
                double kvac = UniversalConstants.Convert.EnergyToElectronWaveNumber(200), u0 = bm.getU(200).Real.Real;
                var s = new Vector3DBase(-d.X, -d.Y, -d.Z);
                var beams = bm.Find_gVectors(Matrix3D.IdentityMatrix, bm.getVecK0(kvac, u0, d, s), d, 64); // 表面の向きを逆に
                Assert.Equal(0, BetheMethod.LastBeamCandidateScan.Collected);
                Assert.Empty(beams);
            }
            finally { BetheMethod.BeamCandidateCut = saved; }
        }
    }
}
