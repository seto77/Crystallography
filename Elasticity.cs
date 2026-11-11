using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace Crystallography
{
    public class Elasticity
    {
        public enum Mode
        {
            Stiffness, Compliance
        }

        public Matrix<double> Compliance = new DenseMatrix(6, 6);
        public Matrix<double> Stiffness = new DenseMatrix(6, 6);

        public Elasticity()
        {
        }

        public Elasticity(Matrix<double> m, Mode mode)
        {
            if (mode == Mode.Compliance)
            {
                Compliance = m;
                Stiffness = m.TryInverse();
            }
            else
            {
                Stiffness = m;
                Compliance = m.TryInverse();
            }
        }

        #region 261105Cl 追加: 結晶系 (Laue 群) による弾性定数の拘束
        /// <summary>261105Cl 追加: 弾性定数の 6×6 行列 (Voigt 表記) に掛ける対称性の拘束の種類。
        /// 行列は ReciPro の結晶の直交座標 (<see cref="Crystal.SetAxis"/>: z ∥ c、b は yz 面内、x ∥ a*。六方晶・三方晶では y ∥ b) で表す。
        /// 典拠: International Tables for Crystallography Vol. D (2003) §1.1.4.10.6 (A. Authier) の行列形。同書の座標は
        /// 単斜晶が 2 回軸 ∥ Ox2、三方晶 32・3m・-3m が 2 回軸 ∥ Ox1 (鏡映面 ⊥ Ox1)、つまり IEEE 流の x ∥ a。
        /// ReciPro の座標では a は x から −30° 回った向き (b ∥ y) なので、321・3m1・-3m1 型 (2 回軸 ∥ a, b) は c15 型 (<see cref="SymmetryConstraint.Trigonal3m2y"/>)、
        /// 312・31m・-31m 型 (2 回軸 ∥ a−b, a+2b, −2a−b のうち −2a−b ∥ −x) は同書と同じ c14 型 (<see cref="SymmetryConstraint.Trigonal3m2x"/>) になる。</summary>
        public enum SymmetryConstraint
        {
            /// <summary>拘束なし (三斜晶・不明)</summary>
            None,
            /// <summary>拘束なし (三方晶の菱面体の軸。3 回軸が z に向かないので、直交座標で単純な形にならない)</summary>
            Rhombohedral,
            /// <summary>単斜晶 a 軸ユニーク (2 回軸 ∥ a ∥ x)。c14・c24・c34・c56 が残る</summary>
            MonoclinicX,
            /// <summary>単斜晶 b 軸ユニーク (2 回軸 ∥ b ∥ y)。c15・c25・c35・c46 が残る</summary>
            MonoclinicY,
            /// <summary>単斜晶 c 軸ユニーク (2 回軸 ∥ c ∥ z)。c16・c26・c36・c45 が残る</summary>
            MonoclinicZ,
            /// <summary>直方晶 (mmm)。独立 9</summary>
            Orthorhombic,
            /// <summary>正方晶 4/m (4, -4, 4/m)。独立 7 (c26 = −c16)</summary>
            Tetragonal4m,
            /// <summary>正方晶 4/mmm (422, 4mm, -42m, 4/mmm)。独立 6</summary>
            Tetragonal4mmm,
            /// <summary>三方晶 -3 (3, -3。六方晶の軸)。独立 7</summary>
            Trigonal3,
            /// <summary>三方晶 -3m (32, 3m, -3m) で 2 回軸 ∥ x (312・31m・-31m 型)。c14 型、独立 6</summary>
            Trigonal3m2x,
            /// <summary>三方晶 -3m (32, 3m, -3m) で 2 回軸 ∥ y (321・3m1・-3m1 型、六方晶の軸の R32・R3m・R-3m)。c15 型、独立 6</summary>
            Trigonal3m2y,
            /// <summary>六方晶 (6/m, 6/mmm)。独立 5</summary>
            Hexagonal,
            /// <summary>立方晶 (m-3, m-3m)。独立 3</summary>
            Cubic,
        }

        /// <summary>261105Cl 追加: 空間群から弾性定数の拘束の種類を決める。</summary>
        public static SymmetryConstraint GetSymmetryConstraint(Symmetry symmetry) => symmetry.CrystalSystemStr switch
        {
            "monoclinic" => symmetry.MainAxis switch
            {
                "a" => SymmetryConstraint.MonoclinicX,
                "b" => SymmetryConstraint.MonoclinicY,
                "c" => SymmetryConstraint.MonoclinicZ,
                _ => SymmetryConstraint.None
            },
            "orthorhombic" => SymmetryConstraint.Orthorhombic,
            "tetragonal" => symmetry.LaueGroupStr == "4/m" ? SymmetryConstraint.Tetragonal4m : SymmetryConstraint.Tetragonal4mmm,
            // 菱面体の軸の判定は Crystal.SetAxis と同じ。-3m の向きは 2 番目の対称方向 (a, b, a+b) に 2 回軸か鏡映面があるか (P321, P3m1, R32 など) で決める
            "trigonal" => symmetry.SpaceGroupHMStr.Contains("Rho") ? SymmetryConstraint.Rhombohedral
                : symmetry.LaueGroupStr == "-3" ? SymmetryConstraint.Trigonal3
                : symmetry.StrSE2p != "1" ? SymmetryConstraint.Trigonal3m2y : SymmetryConstraint.Trigonal3m2x,
            "hexagonal" => SymmetryConstraint.Hexagonal,
            "cubic" => SymmetryConstraint.Cubic,
            _ => SymmetryConstraint.None,
        };

        /// <summary>261105Cl 追加: 6×6 の弾性定数 <paramref name="m"/> (上三角を使う) に結晶系の拘束を掛けた対称行列と、
        /// 拘束で決まる成分 (0 になるか、他の成分から決まる) の印を返す。拘束の無い設定では印は全部 false で、行列は上三角を写しただけ。
        /// コンプライアンスでは Voigt 表記の係数 (剪断の添字 4〜6 を 1 つ含めば 2 倍、2 つ含めば 4 倍。ITC Vol. D 式 1.1.4.12) のため、
        /// 剛性と関係式が違う: c66 = (c11 − c12)/2 に対し s66 = 2(s11 − s12)、c56 = c14 に対し s56 = 2s14、c46 = −c15 に対し s46 = −2s15。</summary>
        public static (double[,] Matrix, bool[,] Fixed) ApplySymmetryConstraint(double[,] m, SymmetryConstraint constraint, Mode mode)
        {
            var r = new double[6, 6];
            for (int i = 0; i < 6; i++)
                for (int j = i; j < 6; j++)
                    r[i, j] = r[j, i] = m[i, j];
            var fix = new bool[6, 6];

            // ij は 1 始まりの 2 桁 (例 14 = c14)
            void set(int ij, double v) { int i = ij / 10 - 1, j = ij % 10 - 1; r[i, j] = r[j, i] = v; fix[i, j] = fix[j, i] = true; }
            double v(int ij) => r[ij / 10 - 1, ij % 10 - 1];
            // void zero(params int[] ijs) { foreach (var ij in ijs) set(ij, 0); } // 261106Ch 変更前
            void zero(params System.ReadOnlySpan<int> ijs) { foreach (var ij in ijs) set(ij, 0); } // 261106Ch: 固定引数の一時配列を不要にする。
            double k = mode == Mode.Compliance ? 2 : 1;          // c56 = c14 ↔ s56 = 2s14、c46 = −c15 ↔ s46 = −2s15
            double f66 = mode == Mode.Compliance ? 2 : 0.5;      // c66 = (c11 − c12)/2 ↔ s66 = 2(s11 − s12)

            switch (constraint)
            {
                case SymmetryConstraint.MonoclinicX: zero(15, 16, 25, 26, 35, 36, 45, 46); break;
                case SymmetryConstraint.MonoclinicY: zero(14, 16, 24, 26, 34, 36, 45, 56); break;
                case SymmetryConstraint.MonoclinicZ: zero(14, 15, 24, 25, 34, 35, 46, 56); break;
                case SymmetryConstraint.Orthorhombic: zero(14, 15, 16, 24, 25, 26, 34, 35, 36, 45, 46, 56); break;
                case SymmetryConstraint.Tetragonal4m:
                    zero(14, 15, 24, 25, 34, 35, 36, 45, 46, 56);
                    set(22, v(11)); set(23, v(13)); set(55, v(44)); set(26, -v(16));
                    break;
                case SymmetryConstraint.Tetragonal4mmm:
                    zero(14, 15, 16, 24, 25, 26, 34, 35, 36, 45, 46, 56);
                    set(22, v(11)); set(23, v(13)); set(55, v(44));
                    break;
                case SymmetryConstraint.Trigonal3 or SymmetryConstraint.Trigonal3m2x or SymmetryConstraint.Trigonal3m2y:
                    if (constraint == SymmetryConstraint.Trigonal3m2x) zero(15); // 2 回軸 ∥ x: c15 = 0 (c14 型)
                    if (constraint == SymmetryConstraint.Trigonal3m2y) zero(14); // 2 回軸 ∥ y: c14 = 0 (c15 型)
                    zero(16, 26, 34, 35, 36, 45);
                    set(22, v(11)); set(23, v(13)); set(24, -v(14)); set(25, -v(15)); set(55, v(44));
                    set(56, k * v(14)); set(46, -k * v(15)); set(66, f66 * (v(11) - v(12)));
                    break;
                case SymmetryConstraint.Hexagonal:
                    zero(14, 15, 16, 24, 25, 26, 34, 35, 36, 45, 46, 56);
                    set(22, v(11)); set(23, v(13)); set(55, v(44)); set(66, f66 * (v(11) - v(12)));
                    break;
                case SymmetryConstraint.Cubic:
                    zero(14, 15, 16, 24, 25, 26, 34, 35, 36, 45, 46, 56);
                    set(22, v(11)); set(33, v(11)); set(13, v(12)); set(23, v(12)); set(55, v(44)); set(66, v(44));
                    break;
            }
            return (r, fix);
        }
        #endregion

        /// <summary>応力 (結晶の直交座標) から歪みテンソルを返す。</summary>
        /// <param name="symmetry">261105Cl: 使わない (互換のため残す)。結晶系ごとの拘束は <see cref="Compliance"/> の 0 と等しい成分に入っている</param>
        public Matrix3D GetStrain(Symmetry symmetry, Matrix3D stress)
        {
            //return symmetry.CrystalSystemNumber switch // 261105Cl 変更前
            //{
            //    0 => GetStrainTriclinic(stress),
            //    1 => GetStrainTriclinic(stress),//tricrinic
            //    2 => GetStrainOrthorhombic(stress),//mono
            //    3 => GetStrainOrthorhombic(stress),//ortho
            //    4 => GetStrainTriclinic(stress),//tetra
            //    5 => GetStrainTrigonal(stress),//trigonal
            //    6 => GetStrainTriclinic(stress),//hexa
            //    7 => GetStrainCubic(stress),//cubic
            //    _ => new Matrix3D(),
            //};
            // 261105Cl 変更: 結晶系によらず 6×6 のコンプライアンス全体で計算する。旧は単斜晶を直方晶の式 (s15・s25・s35・s46 を無視)、
            // 三方晶を IEEE 流の c14 型専用の式 (ReciPro の座標の 321 型は c15 型。下の GetStrainTrigonal)、正方晶・六方晶・三斜晶を
            // 剪断の係数を誤った式 (GetStrainTriclinic の旧版) で計算していた。また回転した結晶子では結晶の座標の応力に全成分が出るので、一般の式が要る
            return GetStrainTriclinic(stress);
        }

        /// <summary>261105Cl 変更: Voigt 表記の一般の式 ε_α = Σ s_αβ σ_β (σ4 = σ23, σ5 = σ13, σ6 = σ12; ε4 = 2ε23, ε5 = 2ε13, ε6 = 2ε12)。
        /// 単位はコンプライアンスと応力に合わせる (GPa⁻¹ × GPa)。
        /// コンプライアンスが無い (弾性定数が未設定で剛性行列の逆行列が無い) か、結果に NaN・∞ が出るときは歪み 0 を返す (NaN を計算に流さない)。</summary>
        public Matrix3D GetStrainTriclinic(Matrix3D stress)
        {
            //double sigma11 = Compliance[0, 0] * stress.E11 + Compliance[0, 1] * stress.E22 + Compliance[0, 2] * stress.E33 + 2 * (Compliance[0, 3] * stress.E23 + Compliance[0, 4] * stress.E31 + Compliance[0, 5] * stress.E12); // 261105Cl 変更前 (剪断応力に余分な 2、剪断歪みの 1/2 が無い)
            //double sigma22 = Compliance[1, 0] * stress.E11 + Compliance[1, 1] * stress.E22 + Compliance[1, 2] * stress.E33 + 2 * (Compliance[1, 3] * stress.E23 + Compliance[1, 4] * stress.E31 + Compliance[1, 5] * stress.E12);
            //double sigma33 = Compliance[2, 0] * stress.E11 + Compliance[2, 1] * stress.E22 + Compliance[2, 2] * stress.E33 + 2 * (Compliance[2, 3] * stress.E23 + Compliance[2, 4] * stress.E31 + Compliance[2, 5] * stress.E12);
            //double sigma23 = Compliance[3, 0] * stress.E11 + Compliance[3, 1] * stress.E22 + Compliance[3, 2] * stress.E33 + 2 * (Compliance[3, 3] * stress.E23 + Compliance[3, 4] * stress.E31 + Compliance[3, 5] * stress.E12);
            //double sigma31 = Compliance[4, 0] * stress.E11 + Compliance[4, 1] * stress.E22 + Compliance[4, 2] * stress.E33 + 2 * (Compliance[4, 3] * stress.E23 + Compliance[4, 4] * stress.E31 + Compliance[4, 5] * stress.E12);
            //double sigma12 = Compliance[5, 0] * stress.E11 + Compliance[5, 1] * stress.E22 + Compliance[5, 2] * stress.E33 + 2 * (Compliance[5, 3] * stress.E23 + Compliance[5, 4] * stress.E31 + Compliance[5, 5] * stress.E12);
            //return new Matrix3D(sigma11, sigma12, sigma31, sigma12, sigma22, sigma23, sigma31, sigma23, sigma33);

            var zero = new Matrix3D(0, 0, 0, 0, 0, 0, 0, 0, 0); // 261105Cl 追加
            if (Compliance == null || Compliance.RowCount != 6 || Compliance.ColumnCount != 6 || stress == null) return zero;
            //double[] t = [stress.E11, stress.E22, stress.E33, (stress.E23 + stress.E32) / 2, (stress.E13 + stress.E31) / 2, (stress.E12 + stress.E21) / 2];
            //var e = new double[6];
            //for (int i = 0; i < 6; i++)
            //    for (int j = 0; j < 6; j++)
            //        e[i] += Compliance[i, j] * t[j];
            // 261111Cl 変更: 配列 2 本を stackalloc にし、密行列なら MathNet のインデクサ (36 回) を通さず列優先の内部配列を直接読む。加算の順序は旧と同じ (結果はビット一致)
            System.Span<double> t = [stress.E11, stress.E22, stress.E33, (stress.E23 + stress.E32) / 2, (stress.E13 + stress.E31) / 2, (stress.E12 + stress.E21) / 2];
            System.Span<double> e = stackalloc double[6];
            var c = Compliance.AsColumnMajorArray(); // 密行列の内部配列 (それ以外は null)
            for (int i = 0; i < 6; i++)
                for (int j = 0; j < 6; j++)
                    e[i] += (c != null ? c[i + 6 * j] : Compliance[i, j]) * t[j];
            foreach (var v in e)
                if (!double.IsFinite(v)) return zero;
            return new Matrix3D(e[0], e[5] / 2, e[4] / 2, e[5] / 2, e[1], e[3] / 2, e[4] / 2, e[3] / 2, e[2]);
        }

        /// <summary>261105Cl 注記: IEEE 流 (2 回軸 ∥ x) の c14 型専用の式。ReciPro の座標の 321・3m1・-3m1 型 (c15 型) や
        /// 回転した結晶子の応力では誤るので <see cref="GetStrain(Symmetry, Matrix3D)"/> からは使わない。</summary>
        public Matrix3D GetStrainTrigonal(Matrix3D stress)
        {
            double sigma11 = Compliance[0, 0] * stress.E11 + Compliance[0, 1] * stress.E22 + Compliance[0, 2] * stress.E33 + Compliance[0, 3] * stress.E23;
            double sigma22 = Compliance[1, 0] * stress.E11 + Compliance[1, 1] * stress.E22 + Compliance[1, 2] * stress.E33 + Compliance[1, 3] * stress.E23;
            double sigma33 = Compliance[2, 0] * (stress.E11 + stress.E22) + Compliance[2, 2] * stress.E33;
            double sigma23 = Compliance[3, 0] * (stress.E11 - stress.E22) + Compliance[3, 3] * stress.E23;
            double sigma31 = Compliance[4, 4] * stress.E31 + Compliance[4, 5] * stress.E12;
            double sigma12 = Compliance[4, 5] * stress.E31 + Compliance[5, 5] * stress.E12;

            return new Matrix3D(sigma11, sigma12 / 2.0, sigma31 / 2.0, sigma12 / 2.0, sigma22, sigma23 / 2.0, sigma31 / 2.0, sigma23 / 2.0, sigma33);
        }

        public Matrix3D GetStrainOrthorhombic(Matrix3D stress)
        {
            double sigma11 = Compliance[0, 0] * stress.E11 + Compliance[0, 1] * stress.E22 + Compliance[0, 2] * stress.E33;
            double sigma22 = Compliance[1, 0] * stress.E11 + Compliance[1, 1] * stress.E22 + Compliance[1, 2] * stress.E33;
            double sigma33 = Compliance[2, 0] * stress.E11 + Compliance[2, 1] * stress.E22 + Compliance[2, 2] * stress.E33;
            double sigma23 = Compliance[3, 3] * stress.E23 / 2;
            double sigma31 = Compliance[4, 4] * stress.E31 / 2;
            double sigma12 = Compliance[5, 5] * stress.E12 / 2;

            return new Matrix3D(sigma11, sigma12, sigma31, sigma12, sigma22, sigma23, sigma31, sigma23, sigma33);
        }

        public Matrix3D GetStrainCubic(Matrix3D stress)
        {
            double sigma11 = Compliance[0, 0] * stress.E11 + Compliance[0, 1] * (stress.E22 + stress.E33);
            double sigma22 = Compliance[0, 0] * stress.E22 + Compliance[0, 1] * (stress.E33 + stress.E11);
            double sigma33 = Compliance[0, 0] * stress.E33 + Compliance[0, 1] * (stress.E11 + stress.E22);
            double sigma23 = Compliance[3, 3] * stress.E23 / 2;
            double sigma31 = Compliance[3, 3] * stress.E31 / 2;
            double sigma12 = Compliance[3, 3] * stress.E12 / 2;

            return new Matrix3D(sigma11, sigma12, sigma31, sigma12, sigma22, sigma23, sigma31, sigma23, sigma33);
        }

        public Matrix3D GetStrain(Symmetry symmetry, Matrix3D rotation, Matrix3D stress)
        {
            Matrix3D m1 = rotation.Transpose() * stress * rotation;
            Matrix3D m2 = GetStrain(symmetry, m1);
            Matrix3D m = rotation * m2 * rotation.Transpose();
            return m;
        }

        public Matrix3D GetStrainByHill(Symmetry symmetry, Matrix3D rotation, Matrix3D stress, Matrix3D strain, double hillFactor)
        {
            return hillFactor * GetStrain(symmetry, rotation, stress) + (1 - hillFactor) * strain;
        }

        /*

        public double[] GetStrain(Symmetry symmetry, double e11, double e12, double e13, double e22, double e23, double e33)
        {
            switch (symmetry.CrystalSystemNumber)
            {
                case 0: return GetStrainTriclinic(e11,e12,e13,e22,e23,e33); break;
                case 1: return GetStrainTriclinic(e11, e12, e13, e22, e23, e33); break;//tricrinic
                case 2: return GetStrainOrthorhombic(e11, e12, e13, e22, e23, e33); break;//mono
                case 3: return GetStrainOrthorhombic(e11, e12, e13, e22, e23, e33); break;//ortho
                case 4: return GetStrainTriclinic(e11, e12, e13, e22, e23, e33); break;//tetra
                case 5: return GetStrainTrigonal(e11, e12, e13, e22, e23, e33); break;//trigonal
                case 6: return GetStrainTriclinic(e11, e12, e13, e22, e23, e33); break;//hexa
                case 7: return GetStrainCubic(e11, e12, e13, e22, e23, e33); break;//cubic
            }
            return new double[] { 0, 0, 0, 0, 0, 0 };
        }

        public double[] GetStrainTriclinic(double e11,double e12,double e13,double e22,double e23,double e33)
        {
            double sigma11 = Compliance.E[0, 0] * e11 + Compliance.E[0, 1] * e22 + Compliance.E[0, 2] * e33 + 2 * (Compliance.E[0, 3] * e23 + Compliance.E[0, 4] * e13 + Compliance.E[0, 5] * e12);
            double sigma22 = Compliance.E[1, 0] * e11 + Compliance.E[1, 1] * e22 + Compliance.E[1, 2] * e33 + 2 * (Compliance.E[1, 3] * e23 + Compliance.E[1, 4] * e13 + Compliance.E[1, 5] * e12);
            double sigma33 = Compliance.E[2, 0] * e11 + Compliance.E[2, 1] * e22 + Compliance.E[2, 2] * e33 + 2 * (Compliance.E[2, 3] * e23 + Compliance.E[2, 4] * e13 + Compliance.E[2, 5] * e12);
            double sigma23 = Compliance.E[3, 0] * e11 + Compliance.E[3, 1] * e22 + Compliance.E[3, 2] * e33 + 2 * (Compliance.E[3, 3] * e23 + Compliance.E[3, 4] * e13 + Compliance.E[3, 5] * e12);
            double sigma31 = Compliance.E[4, 0] * e11 + Compliance.E[4, 1] * e22 + Compliance.E[4, 2] * e33 + 2 * (Compliance.E[4, 3] * e23 + Compliance.E[4, 4] * e13 + Compliance.E[4, 5] * e12);
            double sigma12 = Compliance.E[5, 0] * e11 + Compliance.E[5, 1] * e22 + Compliance.E[5, 2] * e33 + 2 * (Compliance.E[5, 3] * e23 + Compliance.E[5, 4] * e13 + Compliance.E[5, 5] * e12);

            return new double[]{sigma11, sigma12, sigma31, sigma22, sigma23, sigma33};
        }

        public double[] GetStrainTrigonal(double e11, double e12, double e13, double e22, double e23, double e33)
        {
            double sigma11 = Compliance.E[0, 0] * e11 + Compliance.E[0, 1] * e22 + Compliance.E[0, 2] * e33 +  Compliance.E[0, 3] * e23;
            double sigma22 = Compliance.E[1, 0] * e11 + Compliance.E[1, 1] * e22 + Compliance.E[1, 2] * e33 +  Compliance.E[1, 3] * e23;
            double sigma33 = Compliance.E[2, 0] * (e11 + e22) + Compliance.E[2, 2] * e33;
            double sigma23 = Compliance.E[3, 0] * (e11 - e22) + Compliance.E[3, 3] * e23;
            double sigma31 = Compliance.E[4, 4] * e13 + Compliance.E[4, 5] * e12;
            double sigma12 = Compliance.E[4, 5] * e13 + Compliance.E[5, 5] * e12;

            return new double[] {sigma11, sigma12 / 2.0, sigma31 / 2.0, sigma22, sigma23 / 2.0,  sigma33};
        }

        public double[] GetStrainOrthorhombic(double e11, double e12, double e13, double e22, double e23, double e33)
        {
            double sigma11 = Compliance.E[0, 0] * e11 + Compliance.E[0, 1] * e22 + Compliance.E[0, 2] * e33;
            double sigma22 = Compliance.E[1, 0] * e11 + Compliance.E[1, 1] * e22 + Compliance.E[1, 2] * e33;
            double sigma33 = Compliance.E[2, 0] * e11 + Compliance.E[2, 1] * e22 + Compliance.E[2, 2] * e33;
            double sigma23 = Compliance.E[3, 3] * e23 / 2;
            double sigma31 = Compliance.E[4, 4] * e13 / 2;
            double sigma12 = Compliance.E[5, 5] * e12 / 2;

            return new double[] { sigma11, sigma12, sigma31, sigma22, sigma23, sigma33 };
        }

        public double[] GetStrainCubic(double e11, double e12, double e13, double e22, double e23, double e33)
        {
            double sigma11 = Compliance.E[0, 0] * e11 + Compliance.E[0, 1] * (e22 + e33);
            double sigma22 = Compliance.E[0, 0] * e22 + Compliance.E[0, 1] * (e33 + e11);
            double sigma33 = Compliance.E[0, 0] * e33 + Compliance.E[0, 1] * (e11 + e22);
            double sigma23 = Compliance.E[3, 3] * e23 / 2;
            double sigma31 = Compliance.E[3, 3] * e13 / 2;
            double sigma12 = Compliance.E[3, 3] * e12 / 2;

            return new double[] { sigma11, sigma12, sigma31, sigma22, sigma23, sigma33 };
        }

        //s: stress, r: rotation
        public double[] GetStrain(Symmetry symmetry,
            double r11, double r12, double r13, double r21, double r22, double r23, double r31, double r32, double r33,
            double s11, double s12, double s13, double s22, double s23, double s33)
        {
            //Matrix m1 = rotation.Transpose() * stress * rotation;
             //Matrix3D m2 = GetStrain(symmetry, m1);
            double[] m = GetStrain(symmetry,
                r11*r11*s11 + 2*r11*(r21*s12 + r31*s13) + r21*r21*s22 + 2*r21*r31*s23 + r31*r31*s33,
                r12*r21*s12 + r12*r31*s13 + r11*(r12*s11 + r22*s12 + r32*s13) + r21*r22*s22 + r22*r31*s23 + r21*r32*s23 + r31*r32*s33,
                r13*r21*s12 + r13*r31*s13 + r11*(r13*s11 + r23*s12 + r33*s13) + r21*r23*s22 + r23*r31*s23 + r21*r33*s23 + r31*r33*s33,
                r12*r12*s11 + 2*r12*(r22*s12 + r32*s13) + r22*r22*s22 + 2*r22*r32*s23 + r32*r32*s33,
                r13*r22*s12 + r13*r32*s13 + r12*(r13*s11 + r23*s12 + r33*s13) + r22*r23*s22 + r23*r32*s23 + r22*r33*s23 + r32*r33*s33,
                r13*r13*s11 + 2*r13*(r23*s12 + r33*s13) + r23*r23*s22 + 2*r23*r33*s23 + r33*r33*s33
                );

            //Matrix3D m = rotation * m2 * rotation.Transpose();
            //return m;
            return new double[]
            {
                r11*r11*m[0] + 2*r11*(r21*m[1] + r31*m[2]) + r21*r21*m[3] + 2*r21*r31*m[4] + r31*r31*m[5],
               r11*(r21*m[0] + r22*m[1] + r23*m[2]) + r12*(r21*m[1] + r22*m[3] + r23*m[4]) + r13*(r21*m[2] + r22*m[4] + r23*m[5]),
               r11*(r31*m[0] + r32*m[1] + r33*m[2]) + r12*(r31*m[1] + r32*m[3] + r33*m[4]) + r13*(r31*m[2] + r32*m[4] + r33*m[5]),
               r21*r21*m[0] + 2*r21*(r22*m[1] + r23*m[2]) + r22*r22*m[3] + 2*r22*r23*m[4] + r23*r23*m[5],
               r21*(r31*m[0] + r32*m[1] + r33*m[2]) + r22*(r31*m[1] + r32*m[3] + r33*m[4]) + r23*(r31*m[2] + r32*m[4] + r33*m[5]),
               r31*r31*m[0] + 2*r31*(r32*m[1] + r33*m[2]) + r32*r32*m[3] + 2*r32*r33*m[4] + r33*r33*m[5]
            };
        }

        public Matrix3D GetStrainByHill(Symmetry symmetry, Matrix3D rotation, Matrix3D stress, Matrix3D strain, double hillFactor)
        {
            double[] s=GetStrain(
                symmetry,
                rotation.E11,rotation.E12,rotation.E13,rotation.E21,rotation.E22,rotation.E23,rotation.E31,rotation.E32,rotation.E33,
                stress.E11,stress.E12,stress.E13,stress.E22,stress.E23,stress.E33);
            Matrix3D elasticStrain = new Matrix3D(s[0], s[1], s[2], s[1], s[3], s[4], s[2], s[3], s[5]);

            return hillFactor * elasticStrain + (1 - hillFactor) * strain;
        }
        */
    }
}