using System.Drawing;
using Xunit;

namespace Crystallography.Tests;

/// <summary>261118Cl 新規 (codex2 第 3 回の単位 3 の指摘): 吸収の模型 (<see cref="BetheMethod.AbsorptionElasticModel"/>、プロセス全体の設定) の読み方。
/// 以前は getU が鍵を作るときに模型を読み、積分の求積点ごとに (<c>ElasticFactorWithTail</c> の中で) 読み直していたので、
/// 評価の途中で設定が変わると、鍵の模型と値の模型が食い違いえた (プロセス全体で共有する FactorImaginary の記憶にも入る)。
/// 今は入口で 1 度だけ読み、積分の全部へ渡す。設定を書き換えるので、並列化を止めた集まりに入れる。</summary>
[Collection("Bethe global settings")]
public class BetheModelCaptureTests
{
    private const double KV = 200, S2 = 36, M = 0.005; // s² [nm⁻²]・m [nm²] (TdsHybridCheck と同じ点)

    /// <summary>渡した模型が、そのときの全体の設定より優先される (局所・環状・非局所の 3 つの積分)</summary>
    [Fact]
    public void AbsorptionIntegrals_UseThePassedModel_NotTheGlobalSetting()
    {
        var original = BetheMethod.AbsorptionElasticModel;
        try
        {
            var es = AtomStatic.ElectronScatteringPeng[14][0];
            var g = new Vector3DBase(3.2, 1.1, 0.4);
            var h = new Vector3DBase(-1.0, 2.3, 0.2);
            (double Local, double Annular, double NonLocal) eval(AbsorptionElasticModel? model) => (
                es.FactorImaginary(KV, S2, M, model: model),
                es.FactorImaginaryAnnular(KV, g, M, 0.04, 0.12, model: model),
                es.FactorImaginaryAnnular(KV, g, h, M, 0.04, 0.12, model: model));

            BetheMethod.AbsorptionElasticModel = AbsorptionElasticModel.TemariHybrid;
            var temari = eval(null);
            BetheMethod.AbsorptionElasticModel = AbsorptionElasticModel.PengTail;
            var peng = eval(null);
            Assert.NotEqual(temari.Local, peng.Local); // 模型で値が変わる点で確かめる
            Assert.NotEqual(temari.Annular, peng.Annular);
            Assert.NotEqual(temari.NonLocal, peng.NonLocal);

            Assert.Equal(temari, eval(AbsorptionElasticModel.TemariHybrid)); // 全体は PengTail のまま
            BetheMethod.AbsorptionElasticModel = AbsorptionElasticModel.TemariHybrid;
            Assert.Equal(peng, eval(AbsorptionElasticModel.PengTail));
        }
        finally { BetheMethod.AbsorptionElasticModel = original; }
    }

    /// <summary>模型を A → B → A と切り替えても、getU は Clear なしでそれぞれの模型の値を返す (模型は鍵に入っている)。
    /// 参照は、その模型に切り替えてから作った新しいインスタンス</summary>
    [Fact]
    public void GetU_ModelSwitchedAtoBtoA_ReturnsEachModelsValue()
    {
        var original = BetheMethod.AbsorptionElasticModel;
        try
        {
            var crystal = BetheUCacheKeyTestsSilicon(0.5);
            var beam = new BetheMethod.Beam((2, 2, 0), new Vector3DBase(2 / 0.54309, 2 / 0.54309, 0));
            var bethe = new BetheMethod(crystal);

            BetheMethod.AbsorptionElasticModel = AbsorptionElasticModel.TemariHybrid;
            var a1 = bethe.getU(KV, beam);
            var refA = new BetheMethod(crystal).getU(KV, beam);
            BetheMethod.AbsorptionElasticModel = AbsorptionElasticModel.PengTail;
            var b = bethe.getU(KV, beam);
            var refB = new BetheMethod(crystal).getU(KV, beam);
            BetheMethod.AbsorptionElasticModel = AbsorptionElasticModel.TemariHybrid;
            var a2 = bethe.getU(KV, beam);

            Assert.NotEqual(refA.Imag, refB.Imag); // 模型で虚部が変わる条件で確かめる
            Assert.Equal(refA, a1);
            Assert.Equal(refB, b);
            Assert.Equal(refA, a2);
        }
        finally { BetheMethod.AbsorptionElasticModel = original; }
    }

    private static Crystal BetheUCacheKeyTestsSilicon(double bisoAng2)
    {
        int series = 0;
        for (var i = 1; i < SymmetryStatic.TotalSpaceGroupNumber; i++)
            if (SymmetryStatic.Symmetries[i].SpaceGroupNumber == 227) { series = i; break; }
        var cell = (0.54309, 0.54309, 0.54309, Math.PI / 2, Math.PI / 2, Math.PI / 2);
        Atoms[] atoms = [new("Si1", 14, 0, 0, null, series, new Vector3DBase(0, 0, 0), 1,
            new DiffuseScatteringFactor(DiffuseScatteringFactor.Type.B, true, bisoAng2 / 100.0, 0, new double[6], new double[6], cell))];
        return new Crystal(cell, null, series, "Si (BetheModelCaptureTests)", Color.Gray, new Matrix3D(),
            atoms, ("", "", "", ""), Bonds.GetVestaBonds(atoms.Select(x => x.AtomicNumber)));
    }
}
