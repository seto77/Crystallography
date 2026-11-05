using System.Drawing;
using MathNet.Numerics.LinearAlgebra.Double;
using Xunit;
using static Crystallography.Tests.SymmetryReference;
using SC = Crystallography.Elasticity.SymmetryConstraint;

namespace Crystallography.Tests;

/// <summary>261105Cl 新規: 弾性定数の拘束 (<see cref="Elasticity.ApplySymmetryConstraint"/>) と応力→歪み (<see cref="Elasticity.GetStrain(Symmetry, Matrix3D)"/>) の検査。
/// 拘束が、ReciPro の結晶の直交座標 (<see cref="Crystal.SetAxis"/>) で点群の回転に不変な 4 階テンソルの部分空間と一致することを、
/// spglib・cctbx 由来の基準 (symmetry_reference.json) の回転部分で確かめる。典拠は ITC Vol. D (2003) §1.1.4.10.6。</summary>
public class ElasticityTests
{
    #region Voigt 表記 ↔ 4 階テンソル
    private static readonly int[,] V = { { 0, 5, 4 }, { 5, 1, 3 }, { 4, 3, 2 } };
    private static double F(int m) => m < 3 ? 1 : 2;

    /// <summary>Voigt の 6×6 を回転 R (v' = R v) で回す。コンプライアンスは Voigt の係数 (1, 2, 4) を外してから回す</summary>
    private static double[,] Rotate(double[,] c, double[,] R, bool compliance)
    {
        var t = new double[3, 3, 3, 3];
        for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) for (int k = 0; k < 3; k++) for (int l = 0; l < 3; l++)
                        t[i, j, k, l] = c[V[i, j], V[k, l]] / (compliance ? F(V[i, j]) * F(V[k, l]) : 1);
        var r = new double[6, 6];
        int[][] idx = [[0, 0], [1, 1], [2, 2], [1, 2], [0, 2], [0, 1]];
        for (int m = 0; m < 6; m++)
            for (int n = 0; n < 6; n++)
            {
                int i = idx[m][0], j = idx[m][1], k = idx[n][0], l = idx[n][1];
                double s = 0;
                for (int p = 0; p < 3; p++) for (int q = 0; q < 3; q++) for (int u = 0; u < 3; u++) for (int w = 0; w < 3; w++)
                                s += R[i, p] * R[j, q] * R[k, u] * R[l, w] * t[p, q, u, w];
                r[m, n] = s * (compliance ? F(m) * F(n) : 1);
            }
        return r;
    }

    private static double[,] Average(double[,] c, List<double[,]> rotations, bool compliance)
    {
        var sum = new double[6, 6];
        foreach (var R in rotations)
        {
            var r = Rotate(c, R, compliance);
            for (int i = 0; i < 6; i++) for (int j = 0; j < 6; j++) sum[i, j] += r[i, j] / rotations.Count;
        }
        return sum;
    }

    private static double RelDiff(double[,] a, double[,] b)
    {
        double d = 0, m = 0;
        for (int i = 0; i < a.GetLength(0); i++)
            for (int j = 0; j < a.GetLength(1); j++)
            { d = Math.Max(d, Math.Abs(a[i, j] - b[i, j])); m = Math.Max(m, Math.Max(Math.Abs(a[i, j]), Math.Abs(b[i, j]))); }
        return m == 0 ? d : d / m;
    }

    private static double[,] Full(params (int ij, double v)[] comps)
    {
        var c = new double[6, 6];
        foreach (var (ij, v) in comps) { int i = ij / 10 - 1, j = ij % 10 - 1; c[i, j] = c[j, i] = v; }
        return c;
    }

    private static double[,] OnlyIndependent(double[,] m, SC constraint, Elasticity.Mode mode)
    {
        var (_, isFixed) = Elasticity.ApplySymmetryConstraint(m, constraint, mode);
        var r = new double[6, 6];
        for (int i = 0; i < 6; i++) for (int j = i; j < 6; j++) if (!isFixed[i, j]) r[i, j] = m[i, j];
        return r;
    }

    private static int Series(string hm) => Array.FindIndex(SymmetryStatic.Symmetries, s => s.SpaceGroupHMStr == hm);
    #endregion

    /// <summary>直交座標の点群の回転 (基準の回転部分 W を M W M⁻¹。M の列 = 格子ベクトル a, b, c)</summary>
    private static (Crystal Crystal, List<double[,]> Rotations) CartesianRotations(int series, (double, double, double, double, double, double)? cell = null)
    {
        var crystal = new Crystal(cell ?? (5.1, 6.3, 7.7, 81 * Math.PI / 180, 97 * Math.PI / 180, 108 * Math.PI / 180), series, "x", Color.White);
        var M = crystal.MatrixReal;
        var m = DenseMatrix.OfArray(new double[,] { { M.E11, M.E12, M.E13 }, { M.E21, M.E22, M.E23 }, { M.E31, M.E32, M.E33 } });
        var mi = m.Inverse();
        var list = new List<double[,]>();
        foreach (var w in Groups[series].Rotations)
        {
            var W = DenseMatrix.OfArray(new double[,] { { w[0], w[1], w[2] }, { w[3], w[4], w[5] }, { w[6], w[7], w[8] } });
            var R = m * W * mi;
            Assert.True((R.Transpose() * R - DenseMatrix.CreateIdentity(3)).InfinityNorm() < 1e-9, $"series {series}: Cartesian rotation is not orthogonal");
            list.Add(R.ToArray());
        }
        return (crystal, list);
    }

    private static double[,] RandomPositiveDefinite(Random rnd)
    {
        var a = new double[6, 6];
        for (int i = 0; i < 6; i++) for (int j = 0; j < 6; j++) a[i, j] = rnd.NextDouble() * 2 - 1;
        var c = new double[6, 6];
        for (int i = 0; i < 6; i++) for (int j = 0; j < 6; j++) { for (int k = 0; k < 6; k++) c[i, j] += a[i, k] * a[j, k]; if (i == j) c[i, j] += 3; }
        return c;
    }

    #region 拘束 = 点群で不変な部分空間
    /// <summary>全設定で、(1) 拘束は冪等、(2) 点群で平均した (不変な) 剛性・コンプライアンスは拘束を満たす、
    /// (3) 拘束を掛けた行列は点群で不変 (拘束しない三斜晶・菱面体の軸を除く)、(4) 不変な剛性の逆行列はコンプライアンスの拘束 (係数込み) を満たす。
    /// (2) と (3) で、拘束の像 = 不変な部分空間。</summary>
    [Theory, MemberData(nameof(AllSeries), MemberType = typeof(SymmetryReference))]
    public void Constraint_EqualsPointGroupInvariantSubspace(int series)
    {
        var (_, rotations) = CartesianRotations(series);
        var constraint = Elasticity.GetSymmetryConstraint(SymmetryStatic.Symmetries[series]);
        var rnd = new Random(series);
        foreach (var mode in new[] { Elasticity.Mode.Stiffness, Elasticity.Mode.Compliance })
        {
            bool compliance = mode == Elasticity.Mode.Compliance;
            var random = RandomPositiveDefinite(rnd);
            var (constrained, _) = Elasticity.ApplySymmetryConstraint(random, constraint, mode);
            Assert.True(RelDiff(Elasticity.ApplySymmetryConstraint(constrained, constraint, mode).Matrix, constrained) < 1e-14, $"{mode}: not idempotent");

            var invariant = Average(random, rotations, compliance);
            Assert.True(RelDiff(Elasticity.ApplySymmetryConstraint(invariant, constraint, mode).Matrix, invariant) < 1e-10,
                $"series {series} ({constraint}, {mode}): a tensor invariant under the point group does not satisfy the constraint");

            if (constraint is not (SC.None or SC.Rhombohedral))
                Assert.True(RelDiff(Average(constrained, rotations, compliance), constrained) < 1e-10,
                    $"series {series} ({constraint}, {mode}): the constrained tensor is not invariant under the point group");

            if (!compliance)
            {
                var s = DenseMatrix.OfArray(invariant).Inverse().ToArray();
                Assert.True(RelDiff(Elasticity.ApplySymmetryConstraint(s, constraint, Elasticity.Mode.Compliance).Matrix, s) < 1e-10,
                    $"series {series} ({constraint}): the inverse of an invariant stiffness does not satisfy the compliance constraint");
            }
        }
    }

    /// <summary>独立な成分の数 (ITC Vol. D の表と同じ)</summary>
    [Theory]
    [InlineData(SC.None, 21), InlineData(SC.Rhombohedral, 21), InlineData(SC.MonoclinicX, 13), InlineData(SC.MonoclinicY, 13), InlineData(SC.MonoclinicZ, 13),
     InlineData(SC.Orthorhombic, 9), InlineData(SC.Tetragonal4m, 7), InlineData(SC.Tetragonal4mmm, 6), InlineData(SC.Trigonal3, 7),
     InlineData(SC.Trigonal3m2x, 6), InlineData(SC.Trigonal3m2y, 6), InlineData(SC.Hexagonal, 5), InlineData(SC.Cubic, 3)]
    public void Constraint_IndependentComponentCount(SC constraint, int expected)
    {
        var (_, isFixed) = Elasticity.ApplySymmetryConstraint(new double[6, 6], constraint, Elasticity.Mode.Stiffness);
        int n = 0;
        for (int i = 0; i < 6; i++) for (int j = i; j < 6; j++) if (!isFixed[i, j]) n++;
        Assert.Equal(expected, n);
    }

    /// <summary>三方晶 -3m の 2 回軸の向き: ReciPro の座標 (b ∥ y) では 321・3m1・-3m1 型 (と六方晶の軸の R 格子) は 2 回軸 ∥ y、312・31m・-31m 型は ∥ x</summary>
    [Theory]
    [InlineData("P321", SC.Trigonal3m2y), InlineData("P3sub121", SC.Trigonal3m2y), InlineData("P3m1", SC.Trigonal3m2y), InlineData("P-3m1", SC.Trigonal3m2y),
     InlineData("R32Hex", SC.Trigonal3m2y), InlineData("R3mHex", SC.Trigonal3m2y), InlineData("R-3cHex", SC.Trigonal3m2y),
     InlineData("P312", SC.Trigonal3m2x), InlineData("P31m", SC.Trigonal3m2x), InlineData("P-31c", SC.Trigonal3m2x),
     InlineData("R-3mRho", SC.Rhombohedral), InlineData("R3Rho", SC.Rhombohedral), InlineData("R-3Hex", SC.Trigonal3)]
    public void Constraint_TrigonalOrientation(string hm, SC expected)
        => Assert.Equal(expected, Elasticity.GetSymmetryConstraint(SymmetryStatic.Symmetries[Series(hm)]));
    #endregion

    #region 文献の弾性定数
    // α-石英: Ogi et al. (2006) J. Appl. Phys. 100, 053511 (30 °C)。IEEE 流の x ∥ a (2 回軸 ∥ x) の c14 型
    private const double Q11 = 86.76, Q12 = 6.868, Q13 = 11.85, Q14 = -18.02, Q33 = 105.5, Q44 = 58.14;
    private static double[,] QuartzIeee => Full((11, Q11), (22, Q11), (33, Q33), (12, Q12), (13, Q13), (23, Q13), (14, Q14), (24, -Q14), (56, Q14),
        (44, Q44), (55, Q44), (66, (Q11 - Q12) / 2));

    /// <summary>IEEE の座標 (x ∥ a) から ReciPro の座標への回転 (v_ReciPro = R v_IEEE)。z ∥ c は共通</summary>
    private static double[,] IeeeToReciPro(Crystal c)
    {
        double n = Math.Sqrt(c.A_Axis.X * c.A_Axis.X + c.A_Axis.Y * c.A_Axis.Y), ax = c.A_Axis.X / n, ay = c.A_Axis.Y / n;
        return new double[,] { { ax, -ay, 0 }, { ay, ax, 0 }, { 0, 0, 1 } };
    }

    /// <summary>石英 (P3₂21) の文献値 (IEEE の c14 型) は ReciPro の座標では c15 型 (c15 = c14(IEEE)、c14 = 0) になり、点群で不変。IEEE のままでは不変でない</summary>
    [Fact]
    public void Quartz_IeeeC14Type_IsC15TypeInReciProAxes()
    {
        int series = Series("P3sub221");
        var (crystal, rotations) = CartesianRotations(series, (4.913, 4.913, 5.405, Math.PI / 2, Math.PI / 2, 2 * Math.PI / 3));
        Assert.True(RelDiff(Average(QuartzIeee, rotations, false), QuartzIeee) > 1e-3);
        var c = Rotate(QuartzIeee, IeeeToReciPro(crystal), false);
        Assert.True(Math.Abs(c[0, 3]) < 1e-12 && Math.Abs(c[0, 4] - Q14) < 1e-12 && Math.Abs(c[1, 4] + Q14) < 1e-12 && Math.Abs(c[3, 5] + Q14) < 1e-12);
        Assert.True(RelDiff(Average(c, rotations, false), c) < 1e-12);
        Assert.Equal(SC.Trigonal3m2y, Elasticity.GetSymmetryConstraint(crystal.Symmetry));
    }

    /// <summary>文献の剛性 (と、その逆行列のコンプライアンス) を独立な成分だけから補完すると元に戻る (コンプライアンスの係数込み)</summary>
    [Theory]
    [InlineData("Si")]      // Hall (1967): c11 165.7, c12 63.9, c44 79.6 GPa
    [InlineData("Mg")]      // 代表値: c11 59.4, c12 25.6, c13 21.4, c33 61.6, c44 16.4
    [InlineData("rutile")]  // 代表値: c11 268, c12 175, c13 147, c33 484, c44 124, c66 190
    [InlineData("quartz")]  // Ogi et al. (2006)、ReciPro の座標へ回したもの
    public void Literature_CompletionFromIndependentComponents(string name)
    {
        var (hm, c) = name switch
        {
            "Si" => ("Fd-3m(1)", Full((11, 165.7), (22, 165.7), (33, 165.7), (12, 63.9), (13, 63.9), (23, 63.9), (44, 79.6), (55, 79.6), (66, 79.6))),
            "Mg" => ("P6sub3/mmc", Full((11, 59.4), (22, 59.4), (33, 61.6), (12, 25.6), (13, 21.4), (23, 21.4), (44, 16.4), (55, 16.4), (66, (59.4 - 25.6) / 2))),
            "rutile" => ("P4sub2/mnm", Full((11, 268), (22, 268), (33, 484), (12, 175), (13, 147), (23, 147), (44, 124), (55, 124), (66, 190))),
            _ => ("P3sub221", Rotate(QuartzIeee, IeeeToReciPro(new Crystal((4.913, 4.913, 5.405, Math.PI / 2, Math.PI / 2, 2 * Math.PI / 3), Series("P3sub221"), "q", Color.White)), false)),
        };
        var constraint = Elasticity.GetSymmetryConstraint(SymmetryStatic.Symmetries[Series(hm)]);
        var s = DenseMatrix.OfArray(c).Inverse().ToArray();
        Assert.True(RelDiff(Elasticity.ApplySymmetryConstraint(OnlyIndependent(c, constraint, Elasticity.Mode.Stiffness), constraint, Elasticity.Mode.Stiffness).Matrix, c) < 1e-12);
        Assert.True(RelDiff(Elasticity.ApplySymmetryConstraint(OnlyIndependent(s, constraint, Elasticity.Mode.Compliance), constraint, Elasticity.Mode.Compliance).Matrix, s) < 1e-12);
    }
    #endregion

    #region 応力 → 歪み
    /// <summary>Si の一軸応力 σ11 と剪断応力 σ12 (手計算: s11 = (c11+c12)/((c11−c12)(c11+2c12))、s12 = −c12/(…)、ε12 = s44/2 = 1/(2c44))</summary>
    [Fact]
    public void Strain_SiliconUniaxialAndShear()
    {
        double c11 = 165.7, c12 = 63.9, c44 = 79.6, den = (c11 - c12) * (c11 + 2 * c12);
        var el = new Elasticity(DenseMatrix.OfArray(Full((11, c11), (22, c11), (33, c11), (12, c12), (13, c12), (23, c12), (44, c44), (55, c44), (66, c44))), Elasticity.Mode.Stiffness);
        var sym = SymmetryStatic.Symmetries[Series("Fd-3m(1)")];
        var e = el.GetStrain(sym, new Matrix3D(1, 0, 0, 0, 0, 0, 0, 0, 0));
        Assert.Equal((c11 + c12) / den, e.E11, 12);
        Assert.Equal(-c12 / den, e.E22, 12);
        Assert.Equal(-c12 / den, e.E33, 12);
        var g = el.GetStrain(sym, new Matrix3D(0, 1, 0, 1, 0, 0, 0, 0, 0));
        Assert.Equal(1 / (2 * c44), g.E12, 12);
        Assert.Equal(1 / (2 * c44), g.E21, 12);
        Assert.Equal(0, g.E11, 15);
    }

    /// <summary>六方晶の剪断応力 σ23: ε23 = s44/2 = 1/(2 c44) (旧の式は 2 s44 で 4 倍)</summary>
    [Fact]
    public void Strain_HexagonalShear()
    {
        var el = new Elasticity(DenseMatrix.OfArray(Full((11, 59.4), (22, 59.4), (33, 61.6), (12, 25.6), (13, 21.4), (23, 21.4), (44, 16.4), (55, 16.4), (66, (59.4 - 25.6) / 2))), Elasticity.Mode.Stiffness);
        var e = el.GetStrain(SymmetryStatic.Symmetries[Series("P6sub3/mmc")], new Matrix3D(0, 0, 0, 0, 0, 1, 0, 1, 0));
        Assert.Equal(1 / (2 * 16.4), e.E23, 12);
    }

    /// <summary>回転した結晶子に一般の応力: GetStrain(sym, R, σ) = R (S : (Rᵀ σ R)) Rᵀ を 4 階テンソルで独立に計算したものと一致</summary>
    [Fact]
    public void Strain_RotatedCrystallite_MatchesTensorContraction()
    {
        var crystal = new Crystal((4.913, 4.913, 5.405, Math.PI / 2, Math.PI / 2, 2 * Math.PI / 3), Series("P3sub221"), "q", Color.White);
        var c = Rotate(QuartzIeee, IeeeToReciPro(crystal), false);
        var s = DenseMatrix.OfArray(c).Inverse().ToArray();
        var el = new Elasticity(DenseMatrix.OfArray(c), Elasticity.Mode.Stiffness);
        var R = Matrix3D.Rot(new Vector3DBase(0.3, -0.5, 0.8), 0.7);
        var sigma = new Matrix3D(0.3, -0.2, 0.5, -0.2, 1.1, 0.4, 0.5, 0.4, -0.7);
        var e = el.GetStrain(crystal.Symmetry, R, sigma);

        var r = new double[,] { { R.E11, R.E12, R.E13 }, { R.E21, R.E22, R.E23 }, { R.E31, R.E32, R.E33 } };
        var sg = new double[,] { { sigma.E11, sigma.E12, sigma.E13 }, { sigma.E21, sigma.E22, sigma.E23 }, { sigma.E31, sigma.E32, sigma.E33 } };
        var rm = DenseMatrix.OfArray(r);
        var sc = (rm.Transpose() * DenseMatrix.OfArray(sg) * rm).ToArray();
        var ec = new double[3, 3];
        for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) for (int k = 0; k < 3; k++) for (int l = 0; l < 3; l++)
                        ec[i, j] += s[V[i, j], V[k, l]] / (F(V[i, j]) * F(V[k, l])) * sc[k, l];
        var es = (rm * DenseMatrix.OfArray(ec) * rm.Transpose()).ToArray();
        var actual = new double[,] { { e.E11, e.E12, e.E13 }, { e.E21, e.E22, e.E23 }, { e.E31, e.E32, e.E33 } };
        Assert.True(RelDiff(actual, es) < 1e-12);
    }

    /// <summary>弾性定数が全部 0 (逆行列が無い) か NaN を含むとき、応力からの歪みは 0 (NaN・例外を出さない)</summary>
    [Fact]
    public void Strain_WithoutElasticConstants_IsZero()
    {
        var el = new Elasticity(new DenseMatrix(6, 6), Elasticity.Mode.Stiffness);
        Assert.Null(el.Compliance);
        Assert.True(el.GetStrainByHill(SymmetryStatic.Symmetries[1], new Matrix3D(), new Matrix3D(1, 0, 0, 0, 0, 0, 0, 0, 0), new Matrix3D(0, 0, 0, 0, 0, 0, 0, 0, 0), 1).IsZero());
        var nan = DenseMatrix.CreateIdentity(6); nan[0, 0] = double.NaN;
        Assert.True(new Elasticity(nan, Elasticity.Mode.Compliance).GetStrain(SymmetryStatic.Symmetries[1], new Matrix3D(1, 0, 0, 0, 0, 0, 0, 0, 0)).IsZero());
    }
    #endregion

    #region 261106Cl 追加: Hill 係数の既定 (作者の決定 C6) と、入力したコンプライアンスの保存 (C8-3)
    private static readonly (double, double, double, double, double, double) CubicCell = (0.5431, 0.5431, 0.5431, Math.PI / 2, Math.PI / 2, Math.PI / 2);

    /// <summary>結晶リスト (XML) に書いて読み戻す (アプリと同じ ConvertCrystalData の経路)。<paramref name="edit"/> で保存後の XML 本文を書き換えられる</summary>
    private static Crystal[] RoundTrip(Crystal[] crystals, Func<string, string> edit = null)
    {
        var path = Path.Combine(Path.GetTempPath(), $"recipro_test_{Guid.NewGuid():N}.xml");
        try
        {
            Assert.True(ConvertCrystalData.SaveCrystalListXml(crystals, path));
            if (edit != null) File.WriteAllText(path, edit(File.ReadAllText(path)));
            return ConvertCrystalData.ConvertToCrystalList(path);
        }
        finally { File.Delete(path); }
    }

    /// <summary>新しく作る結晶 (格子定数を受け取るコンストラクタ・データベースやアプリ間の受け渡し (Crystal2)) の Hill 係数は 0.5。
    /// 引数なしのコンストラクタ (XML の読み込みが使う) は 0 のまま、複製は元の値を引き継ぐ</summary>
    [Fact]
    public void HillCoefficient_NewCrystalIsHalf()
    {
        var c = new Crystal(CubicCell, Series("Fd-3m(1)"), "Si", Color.White);
        Assert.Equal(0.5, c.HillCoefficient);
        Assert.Equal(0.5, c.ToCrystal2().ToCrystal().HillCoefficient);
        Assert.Equal(0.0, new Crystal().HillCoefficient);
        c.HillCoefficient = 0.25;
        Assert.Equal(0.25, new Crystal(c).HillCoefficient);
        Assert.Equal(0.25, ((Crystal)c.Clone()).HillCoefficient);
    }

    /// <summary>保存済みの結晶は保存された値のまま読む。0 は XML に書かれない (要素が無い) ので、要素の無い古いデータは 0 として読まれる (既定の 0.5 にならない)</summary>
    [Fact]
    public void HillCoefficient_SavedValueIsKept()
    {
        double[] values = [0, 0.5, 1, 0.3];
        var crystals = values.Select((h, i) => new Crystal(CubicCell, Series("Fd-3m(1)"), $"c{i}", Color.White) { HillCoefficient = h }).ToArray();
        string xml = null;
        var read = RoundTrip(crystals, s => xml = s);
        Assert.Equal(3, System.Text.RegularExpressions.Regex.Matches(xml, "<HillCoefficient>").Count); // 0 だけ書かれない
        Assert.Equal(values, read.Select(c => c.HillCoefficient).ToArray());
    }

    /// <summary>コンプライアンスで入力した値 (逆行列の無いものも) を保存して読み戻せる。剛性を設定すると消える。要素の無い古いデータは null</summary>
    [Fact]
    public void ElasticCompliance_IsStoredAndRead()
    {
        // 逆行列の無いコンプライアンス (立方晶で s44 を入れていない)
        var s = Full((11, 0.0077), (22, 0.0077), (33, 0.0077), (12, -0.0021), (13, -0.0021), (23, -0.0021));
        Assert.Null(DenseMatrix.OfArray(s).TryInverse());
        var c = new Crystal(CubicCell, Series("Fd-3m(1)"), "singular", Color.White) { ElasticStiffness = new double[6, 6] };
        c.ElasticCompliance = s;
        var plain = new Crystal(CubicCell, Series("Fd-3m(1)"), "plain", Color.White) { ElasticStiffness = Full((11, 165.7), (22, 165.7), (33, 165.7), (12, 63.9), (13, 63.9), (23, 63.9), (44, 79.6), (55, 79.6), (66, 79.6)) };
        var read = RoundTrip([c, plain]);
        Assert.Equal(0, RelDiff(read[0].ElasticCompliance, s));
        Assert.All(read[0].ElasticStiffnessArray, v => Assert.Equal(0, v));
        Assert.Null(read[1].ElasticCompliance);
        Assert.Equal(0, RelDiff(read[1].ElasticStiffness, plain.ElasticStiffness));
        // 複製は引き継ぐ (別の配列)
        var copy = new Crystal(c);
        Assert.Equal(0, RelDiff(copy.ElasticCompliance, s));
        Assert.NotSame(copy.ElasticComplianceArray, c.ElasticComplianceArray);
        // 剛性を設定すると消える (食い違いを残さない)
        c.ElasticStiffness = plain.ElasticStiffness;
        Assert.Null(c.ElasticCompliance);
        // 逆行列の無いコンプライアンスでも応力から歪みを出せる (ε = S σ。σ11 = 1 GPa → ε11 = s11, ε22 = ε33 = s12)
        var e = new Elasticity(DenseMatrix.OfArray(s), Elasticity.Mode.Compliance).GetStrain(SymmetryStatic.Symmetries[1], new Matrix3D(1, 0, 0, 0, 0, 0, 0, 0, 0));
        Assert.Equal(0.0077, e.E11, 15);
        Assert.Equal(-0.0021, e.E22, 15);
        Assert.Equal(-0.0021, e.E33, 15);
    }

    /// <summary>新しい要素 (ElasticComplianceArray など) を知らない古い版は、その要素を読み飛ばす (XmlSerializer は知らない要素で失敗しない)。
    /// ここでは逆に、今の版が知らない要素を足した XML を読んで確かめる</summary>
    [Fact]
    public void CrystalListXml_UnknownElementIsIgnored()
    {
        var c = new Crystal(CubicCell, Series("Fd-3m(1)"), "x", Color.White) { HillCoefficient = 0.5 };
        var read = RoundTrip([c], s => s.Replace("<HillCoefficient>", "<SomeFutureArray><double>1</double></SomeFutureArray><HillCoefficient>"));
        Assert.Single(read);
        Assert.Equal(0.5, read[0].HillCoefficient);
    }
    #endregion
}
