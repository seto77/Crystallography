using Xunit;

namespace Crystallography.Tests;

/// <summary>261118Cl 新規 (codex2 第 5 回「符号の公開回帰は今入れる価値がある」): 非局所の吸収の積分 (<see cref="AtomStatic.ES.FactorImaginaryAnnular(double, Vector3DBase, Vector3DBase, double, double, double, int, int, bool, AbsorptionElasticModel?)"/>)
/// の座標の約束。積分器は入射を K₀ = +k₀ẑ と決め打ちし、ビーム K₀ + g から終状態への散乱を組む。殻の上のビーム (|K₀ + g| = k₀) では、
/// 平面波 K₀ + g の全散乱断面積は向きに依らないので、全球の対角 U′(g,g) は U′(0,0) に等しい。これが成り立つのは g の z 成分を
/// 入射の向き (+z) で測ったときだけで、z の符号を逆にすると外れる。EBSD のマスターパターン (結晶の座標系を固定し、方向ごとに入射を変える) で
/// 非局所の吸収を直すときは、この約束に合わせて −Rg, −Rh (= 入射 K̂₀ を +ẑ へ送る変換) を渡す (ワークフロー wf_axis_261118 の結論)</summary>
public class NonLocalAbsorptionAxisTests
{
    [Theory]
    [InlineData(20.0, 10.0)]
    [InlineData(20.0, 20.0)]
    [InlineData(10.0, 15.0)]
    public void OnShellDiagonal_EqualsForwardValue_OnlyWhenZIsTheIncidentDirection(double kV, double gt)
    {
        var es = AtomStatic.ElectronScatteringPeng[14][0];
        const double m = 0.005; // nm² (B 0.5 Å²)
        var model = AbsorptionElasticModel.TemariHybrid;
        double k0 = UniversalConstants.Convert.EnergyToElectronWaveNumber(kV);
        double gz = Math.Sqrt(k0 * k0 - gt * gt) - k0; // |k₀ẑ + g| = k₀ (殻の上)
        var onShell = new Vector3DBase(gt * Math.Cos(0.3), gt * Math.Sin(0.3), gz);
        var mirrored = new Vector3DBase(onShell.X, onShell.Y, -gz); // z の符号を逆にした (入射と逆向きの軸で測った) ベクトル
        var zero = new Vector3DBase(0, 0, 0);
        // 山 (K_f = K₀ + g) が球の上に乗るので、製品の (60, 20) では φ の求積の誤差が 1e-3 を超える (反証の担当の測定)。細かい点数で比べる
        const int nTheta = 256, nPhi = 128;
        double u00 = es.FactorImaginaryAnnular(kV, zero, zero, m, 0, Math.PI, nTheta, nPhi, model: model);
        double ugg = es.FactorImaginaryAnnular(kV, onShell, onShell, m, 0, Math.PI, nTheta, nPhi, model: model);
        double umm = es.FactorImaginaryAnnular(kV, mirrored, mirrored, m, 0, Math.PI, nTheta, nPhi, model: model);
        double relOn = ugg / u00 - 1, relMirror = umm / u00 - 1;
        Assert.True(Math.Abs(relOn) < 1e-5, $"{kV} kV |g_t| {gt}: on-shell U'(g,g)/U'(0,0) − 1 = {relOn:e2}");
        Assert.True(Math.Abs(relMirror) > 1e-3, $"{kV} kV |g_t| {gt}: mirrored U'(g,g)/U'(0,0) − 1 = {relMirror:e2} (z の符号で区別できない)");
    }
}
