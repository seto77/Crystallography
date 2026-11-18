using System.Drawing;
using System.Numerics;
using Xunit;

namespace Crystallography.Tests;

/// <summary>261118Cl 新規 (codex2 第 3・4 回の単位 4): ALCHEMI を専用の BetheMethod で計算すること (<see cref="BetheMethod.PrepareAlchemi"/>) の回帰。
/// 以前の RunAlchemi は所有元の BetheMethod の共有の状態 (AccVoltage・BaseRotation・MaxNumOfBloch・Surface・Thicknesses・getU のキャッシュ) を
/// 書き換えて計算し、要求の配列も計算の途中で読み続けていた。ここでは、進捗の通知 (計算のスレッドで同期に呼ばれる) の中で割り込みを起こし、
/// 単独の実行と一致することを確かめる。全体の設定を書き換えるので、並列化を止めた集まりに入れる。</summary>
[Collection("Bethe global settings")]
public class BetheAlchemiIsolationTests
{
    private const double KV = 200;

    /// <summary>同じ所有元で、ALCHEMI の途中に電圧・方位・本数の違う SAED を計算しても、どちらも単独の実行と一致する</summary>
    [Fact]
    public void SaedInterleavedDuringAlchemi_BothMatchStandalone()
    {
        var crystal = SrTiO3();
        var saedRotation = Matrix3D.Rot(new Vector3DBase(1, 0.3, 0.1), 0.02);
        var expectedAlchemi = new BetheMethod(crystal).RunAlchemi(Request());
        var expectedSaed = new BetheMethod(crystal).GetDifractedBeamAmpriltudes(48, 300, saedRotation, 30);

        var owner = new BetheMethod(crystal);
        BetheMethod.Beam[] interleavedSaed = null;
        var actualAlchemi = owner.RunAlchemi(Request(), p =>
        {
            // 基底を組み終えてポテンシャル行列を組む前 (以前の形では、ここで所有元の AccVoltage が 300 kV に替わった)
            if (interleavedSaed == null && p.Stage == AlchemiStage.BuildingUnionBasis && p.Fraction == 1)
                interleavedSaed = owner.GetDifractedBeamAmpriltudes(48, 300, saedRotation, 30);
        });

        Assert.NotNull(interleavedSaed);
        AssertSameAlchemi(expectedAlchemi, actualAlchemi, "ALCHEMI with SAED interleaved");
        AssertSameBeams(expectedSaed, interleavedSaed, "SAED interleaved in ALCHEMI");
        AssertSameBeams(expectedSaed, owner.GetDifractedBeamAmpriltudes(48, 300, saedRotation, 30), "SAED after ALCHEMI");
    }

    /// <summary>入力を確定した後に、呼び出し側の要求の配列・元の結晶の温度因子・全体の設定を変えても、その run は確定時の入力で計算する。
    /// 次に確定した run は新しい条件に応じる</summary>
    [Fact]
    public void PreparedRun_IsUnaffectedByLaterEdits()
    {
        var originalAbsorption = BetheMethod.AbsorptionElasticModel;
        var originalTie = BetheMethod.BeamBoundaryTieMode;
        try
        {
            var crystal = SrTiO3();
            var expected = new BetheMethod(SrTiO3()).RunAlchemi(Request());

            var owner = new BetheMethod(crystal);
            var request = Request();
            var run = owner.PrepareAlchemi(request);

            // 確定の後に呼び出し側の入力を書き換える
            request.Orientations[0].BeamDirection.X += 0.01;
            request.Sites[0].AtomsIndices[0] = 1;
            request.ThicknessesNm[1] = 99;
            request.BaseRotation.E11 = 0.5;
            crystal.Atoms[0].Dsf.Iso *= 3;
            BetheMethod.AbsorptionElasticModel = AbsorptionElasticModel.PengTail;
            BetheMethod.BeamBoundaryTieMode = BeamBoundaryTieMode.DropGroup;

            var actual = run.Run();
            AssertSameAlchemi(expected, actual, "prepared run after the caller edited its inputs");
            Assert.Equal(20.0, actual.ThicknessesNm[1]);
            Assert.Equal(0, actual.Sites[0].AtomsIndices[0]);

            // 次の run は新しい条件で計算する (ここでは厚さだけ確かめる。ほかの書き換えも効くので値は別物になる)
            BetheMethod.AbsorptionElasticModel = originalAbsorption;
            BetheMethod.BeamBoundaryTieMode = originalTie;
            var changed = Request() with { ThicknessesNm = [0, 99, 50] };
            var next = owner.PrepareAlchemi(changed).Run();
            Assert.Equal(99.0, next.ThicknessesNm[1]);
            Assert.True(MaxRel(expected.Dynamic, next.Dynamic) > 1e-6, "温度因子を変えた結晶の次の run は値が変わる");
        }
        finally
        {
            BetheMethod.AbsorptionElasticModel = originalAbsorption;
            BetheMethod.BeamBoundaryTieMode = originalTie;
        }
    }

    /// <summary>途中で取り消した後、同じ所有元の SAED と ALCHEMI の再実行が単独の実行と一致する</summary>
    [Fact]
    public void CancelledRun_LeavesOwnerClean()
    {
        var crystal = SrTiO3();
        var saedRotation = Matrix3D.Rot(new Vector3DBase(0.2, 1, 0), 0.03);
        var expectedAlchemi = new BetheMethod(crystal).RunAlchemi(Request());
        var expectedSaed = new BetheMethod(crystal).GetDifractedBeamAmpriltudes(40, 120, saedRotation, 25);

        var owner = new BetheMethod(crystal);
        using var cts = new CancellationTokenSource();
        Assert.ThrowsAny<OperationCanceledException>(() => owner.RunAlchemi(Request(), p =>
        {
            if (p.Stage == AlchemiStage.SolvingOrientations) cts.Cancel();
        }, cts.Token));

        AssertSameBeams(expectedSaed, owner.GetDifractedBeamAmpriltudes(40, 120, saedRotation, 25), "SAED after a cancelled ALCHEMI");
        AssertSameAlchemi(expectedAlchemi, owner.RunAlchemi(Request()), "ALCHEMI rerun after cancel");
    }

    /// <summary>結晶の写し (<c>new Crystal(source)</c>、PrepareAlchemi が使う) は、非直交の格子・異方性の温度因子・対称操作のある結晶でも、
    /// 逆格子・体積・等価位置・U を元と同じに保つ。所有元の上で本体を直接走らせた結果とも一致する</summary>
    [Fact]
    public void CrystalCopy_KeepsGeometryPotentialsAndAlchemi()
    {
        var source = MonoclinicAniso();
        var copy = new Crystal(source);
        Assert.Equal(source.Volume, copy.Volume, 15);
        foreach (var (s, c) in new[] { (source.A_Star, copy.A_Star), (source.B_Star, copy.B_Star), (source.C_Star, copy.C_Star) })
            Assert.True((s - c).Length < 1e-12 * s.Length, $"逆格子が変わった: {s} vs {c}");
        Assert.Equal(source.Atoms.Length, copy.Atoms.Length);
        for (int i = 0; i < source.Atoms.Length; i++)
        {
            Assert.Equal(source.Atoms[i].Atom.Length, copy.Atoms[i].Atom.Length);
            for (int j = 0; j < source.Atoms[i].Atom.Length; j++)
                Assert.True((source.Atoms[i].Atom[j] - copy.Atoms[i].Atom[j]).Length < 1e-14, $"等価位置が変わった ({i},{j})");
        }
        var bs = new BetheMethod(source); var bc = new BetheMethod(copy);
        foreach (var hkl in new[] { (1, 0, 0), (0, 1, 1), (2, -1, 1), (1, 2, -3), (0, 0, 4) })
        {
            var us = bs.getU(KV, new BetheMethod.Beam(hkl, source.MatrixInverseTransposed * hkl), cache: false);
            var uc = bc.getU(KV, new BetheMethod.Beam(hkl, copy.MatrixInverseTransposed * hkl), cache: false);
            Assert.True((us.Real - uc.Real).Magnitude <= 1e-14 * Math.Max(1, us.Real.Magnitude) && (us.Imag - uc.Imag).Magnitude <= 1e-14 * Math.Max(1, us.Imag.Magnitude), $"U が変わった {hkl}");
        }

        var request = new AlchemiRequest(KV, new Matrix3D(), AlchemiScan.Rocking1D(new Vector3DBase(0, 0, -1), new Vector3DBase(1, 0, 0), -6e-3, 6e-3, 5), [0, 30],
            [new AlchemiSiteBasis("Ti", [0]), new AlchemiSiteBasis("O", [1])], [new IonizationChannelSpec(22, IonizationShell.K)]) { MaxNumOfBloch = 40 };
        var owner = new BetheMethod(source);
        var run = owner.PrepareAlchemi(request);
        var viaCopy = run.Run();
        var viaOwner = owner.RunAlchemiCore(run.Request, null, default); // 所有元の上で本体を直接 (テストの照合だけ)
        AssertSameAlchemi(viaOwner, viaCopy, "ALCHEMI on the crystal copy vs on the owner");
    }

    //---- 補助 ----

    private static AlchemiRequest Request() => new(KV, new Matrix3D(),
        AlchemiScan.Rocking1D(new Vector3DBase(0, 0, -1), new Vector3DBase(1, 0, 0), -8e-3, 8e-3, 9), [0, 20, 50],
        [new AlchemiSiteBasis("Sr", [0]), new AlchemiSiteBasis("Ti", [1]), new AlchemiSiteBasis("O", [2, 3, 4])],
        [new IonizationChannelSpec(22, IonizationShell.K), new IonizationChannelSpec(38, IonizationShell.LTotal)])
    { MaxNumOfBloch = 40 };

    private static double MaxRel(double[] a, double[] b)
    {
        Assert.Equal(a.Length, b.Length);
        double scale = a.Select(Math.Abs).DefaultIfEmpty(0).Max();
        return a.Zip(b, (x, y) => Math.Abs(x - y)).DefaultIfEmpty(0).Max() / (scale > 0 ? scale : 1);
    }

    /// <summary>形・有限性・基底・成分ごとの値・拡張基底の判定を比べる (fresh 同士の再現差は 0。丸めの範囲として 1e-12)</summary>
    private static void AssertSameAlchemi(AlchemiResult expected, AlchemiResult actual, string what)
    {
        Assert.Equal(expected.Shape, actual.Shape);
        foreach (var arr in new[] { actual.Dynamic, actual.Dechannelled, actual.CoherentPathLengthNm })
            Assert.True(arr.All(double.IsFinite), $"{what}: 非有限の値");
        Assert.Equal(expected.Basis.BasisHash, actual.Basis.BasisHash);
        Assert.Equal(expected.Basis.AcceptedForFit, actual.Basis.AcceptedForFit);
        double d = MaxRel(expected.Dynamic, actual.Dynamic), e = MaxRel(expected.Dechannelled, actual.Dechannelled),
            l = MaxRel(expected.CoherentPathLengthNm, actual.CoherentPathLengthNm),
            x = Math.Abs(expected.Basis.ExpandedBasisMaxRelDiff - actual.Basis.ExpandedBasisMaxRelDiff);
        Assert.True(d <= 1e-12 && e <= 1e-12 && l <= 1e-12 && x <= 1e-12, $"{what}: dyn {d:e2} dech {e:e2} lcoh {l:e2} expanded {x:e2}");
    }

    private static void AssertSameBeams(BetheMethod.Beam[] expected, BetheMethod.Beam[] actual, string what)
    {
        Assert.Equal(expected.Length, actual.Length);
        double max = expected.Max(b => b.Psi.Magnitude);
        for (int i = 0; i < expected.Length; i++)
        {
            Assert.Equal(expected[i].Index, actual[i].Index);
            Assert.True((expected[i].Psi - actual[i].Psi).Magnitude <= 1e-12 * max, $"{what}: beam {expected[i].Index} {expected[i].Psi} vs {actual[i].Psi}");
        }
    }

    /// <summary>SrTiO₃ を P1 で 5 サイト明示 (tools/AlchemiCheck の Fixtures.MakeSrTiO3 と同じ)</summary>
    private static Crystal SrTiO3()
    {
        const double a = 0.3905;
        var cell = (a, a, a, Math.PI / 2, Math.PI / 2, Math.PI / 2);
        var series = SeriesOf(1);
        Atoms Make(string name, int z, double biso, double x, double y, double zz) => new(name, z, 0, 0, null, series, new Vector3DBase(x, y, zz), 1,
            new DiffuseScatteringFactor(DiffuseScatteringFactor.Type.B, true, biso, 0, new double[6], new double[6], cell));
        Atoms[] atoms = [Make("Sr", 38, 0.008, 0, 0, 0), Make("Ti", 22, 0.006, 0.5, 0.5, 0.5), Make("O1", 8, 0.010, 0.5, 0.5, 0), Make("O2", 8, 0.010, 0.5, 0, 0.5), Make("O3", 8, 0.010, 0, 0.5, 0.5)];
        return new Crystal(cell, null, series, "SrTiO3 (BetheAlchemiIsolationTests)", Color.Gray, new Matrix3D(), atoms, ("", "", "", ""), []);
    }

    /// <summary>単斜晶 P2₁/c、一般位置の Ti と O に異方性の温度因子 (B 型の無次元係数)</summary>
    private static Crystal MonoclinicAniso()
    {
        var cell = (0.512, 0.634, 0.719, Math.PI / 2, 103.0 * Math.PI / 180, Math.PI / 2);
        var series = SeriesOf(14);
        Atoms Make(string name, int z, double x, double y, double zz, double[] aniso) => new(name, z, 0, 0, null, series, new Vector3DBase(x, y, zz), 1,
            new DiffuseScatteringFactor(DiffuseScatteringFactor.Type.B, false, 0, 0, aniso, new double[6], cell));
        Atoms[] atoms = [Make("Ti1", 22, 0.13, 0.27, 0.41, [0.004, 0.003, 0.002, 0.0005, -0.0003, 0.0008]),
            Make("O1", 8, 0.31, 0.07, 0.19, [0.006, 0.005, 0.004, -0.0007, 0.0004, 0.0010])];
        return new Crystal(cell, null, series, "monoclinic aniso (BetheAlchemiIsolationTests)", Color.Gray, new Matrix3D(), atoms, ("", "", "", ""), []);
    }

    private static int SeriesOf(int itNumber)
    {
        for (var i = 1; i < SymmetryStatic.TotalSpaceGroupNumber; i++)
            if (SymmetryStatic.Symmetries[i].SpaceGroupNumber == itNumber)
                return i;
        throw new ArgumentException($"No space group with IT number {itNumber}.");
    }
}
