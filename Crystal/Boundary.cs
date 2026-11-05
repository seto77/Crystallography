using System;
using System.Drawing;
using System.Linq;
using System.Xml.Serialization;

namespace Crystallography
{
    /// <summary>境界を表現するクラス. 主にReciProのStructure Viewerから呼び出される.</summary>
    [Serializable()]
    public class Bound
    {
        public (int H, int K, int L) Index { get; set; }

        [XmlIgnore]
        public Color Color { get => Color.FromArgb(ColorArgb); }
        public int ColorArgb { get; set; } = Color.Gray.ToArgb();

        public bool Equivalency { get; set; } = true;
        /// <summary>261106Cl 追加 (作者決定 D5): 境界の面。(X, Y, Z) は面の単位法線 n、D は投影中心から面までの距離 (nm)。
        /// 面は投影中心から n の向きの距離 D にあり、残る (描かれる) 側は n·r ≤ D (r は投影中心からの座標)。
        /// ReciPro の Structure Viewer は内側を「n·r + D ≥ 0」で扱うので、(−n, D) にして使う</summary>
        public (double X, double Y, double Z, double D)[] PlaneParams { get; set; }
        public double Distance { get; set; }
        public double MultipleOfD { get; set; }

        /// <summary>平行移動量(nm)。261106Cl: 正の値で、面で囲んだ領域が指定した (hkl) の法線の向き (+ĝ) へ動く</summary>
        public double Translation { get; set; }

        public bool Enabled { get; set; } = true;

        public Bound()
        { }

        /// <summary>Boundクラスを生成. distanceはnm単位であることに注意</summary>
        /// <param name="enabled"></param>
        /// <param name="crystal"></param>
        /// <param name="h"></param>
        /// <param name="k"></param>
        /// <param name="l"></param>
        /// <param name="equivalency"></param>
        /// <param name="distance"></param>
        /// <param name="translation"></param>
        /// <param name="argb"></param>
        public Bound(bool enabled, Crystal crystal, int h, int k, int l, bool equivalency, double distance, double translation, int argb) : this()
        {
            Enabled = enabled;
            ColorArgb = argb;
            Equivalency = equivalency;
            Distance = distance;
            Translation = translation;
            Index = (h, k, l);

            Reset(crystal);
        }

        public void Reset(Crystal crystal)
        {
            var gBase = Index.H * crystal.A_Star + Index.K * crystal.B_Star + Index.L * crystal.C_Star;
            MultipleOfD = Distance * gBase.Length;

            (int H, int K, int L)[] planes = Equivalency ? SymmetryStatic.GenerateEquivalentPlanes(Index, crystal.Symmetry, false) : [Index];

            // 261105Cl 変更前: Translation が効くのは「等価な面を考慮」かつ等価な面がちょうど 2 枚のときだけで、3 枚以上・1 枚では無視されていた。
            //   また 2 枚のときは並び替えた先頭の面 (基底の指数) に +T を付けたので、指定した (hkl) が先頭でないと向きが逆になった
            //if (Equivalency && planes.Length == 2 && Translation != 0)
            //{
            //    var g1 = planes[0].H * crystal.A_Star + planes[0].K * crystal.B_Star + planes[0].L * crystal.C_Star;
            //    var g2 = planes[1].H * crystal.A_Star + planes[1].K * crystal.B_Star + planes[1].L * crystal.C_Star;
            //    g1.NormarizeThis();
            //    g2.NormarizeThis();
            //    PlaneParams = [(g1.X, g1.Y, g1.Z, Distance + Translation), (g2.X, g2.Y, g2.Z, Distance - Translation)];
            //}
            //else
            //{
            //    PlaneParams = planes.Select(p =>
            //    {
            //        var g = p.H * crystal.A_Star + p.K * crystal.B_Star + p.L * crystal.C_Star;
            //        g.NormarizeThis();
            //        return (g.X, g.Y, g.Z, Distance);
            //    }).ToArray();
            //}
            // 261105Cl 変更: 面で囲んだ領域全体を、指定した (hkl) の法線 n₀ の向きに Translation だけ平行移動する。各面 i (単位法線 nᵢ) の
            //   条件 nᵢ·r ≤ D は、移動後 nᵢ·(r − T n₀) ≤ D、すなわち nᵢ·r ≤ D + T (nᵢ·n₀) になる。どの面数でも効き、2 枚 (n₀, −n₀) なら
            //   D+T / D−T (旧と同じ)、等価な面を考えない 1 枚なら D+T (面を法線の向きに T だけずらす) になる
            // 261106Cl 変更 (作者決定 D5): 面の約束を、上の想定どおり「各面 i について nᵢ·r ≤ D」(r は投影中心からの座標。面は中心から nᵢ の向きの距離 D にあり、
            //   中心の側が残る) にそろえた (<see cref="PlaneParams"/>)。式 D + T (nᵢ·n₀) がそのまま正しく、正の T で領域が +ĝ 側へ動く (作者決定 C5 を保つ)。
            //   面が 1 枚だけ (等価な面を考えない・反転中心の無い結晶で 1 枚) なら、(hkl) の面が +ĝ 側の距離 D (+T) に来て、投影中心の側が残る。
            //   描画側 (Structure Viewer・Geometry.GetClippedPolygon) は内側を n·r + D ≥ 0 で扱うので、Structure Viewer が (−nᵢ, D) にして渡す。
            //   (C5 の時点では描画側の約束のまま D − T (nᵢ·n₀) にしていた。等価な面が ±n の組でそろう閉じた多面体は同じ領域だが、
            //    1 枚だけの面は −ĝ 側の距離 D に置かれて +ĝ 側が残り、反転中心の無い結晶の等価な面 (四面体など) は反転した形になっていた)
            var n0 = gBase.Length > 0 ? gBase / gBase.Length : gBase; // 261105Cl 追加 (Index が (0,0,0) なら 0 ベクトル。描画側で除外される)
            PlaneParams = planes.Select(p =>
            {
                var g = p.H * crystal.A_Star + p.K * crystal.B_Star + p.L * crystal.C_Star;
                g.NormarizeThis();
                return (g.X, g.Y, g.Z, Distance + Translation * (g.X * n0.X + g.Y * n0.Y + g.Z * n0.Z)); // 261106Cl 変更 (作者決定 D5): 約束 nᵢ·r ≤ D で、正の T で +ĝ 側へ
            }).ToArray();
        }
    }

    /// <summary>境界を表現するクラス. 主にReciProのStructure Viewerから呼び出される.</summary>
    [Serializable()]
    public class LatticePlane
    {
        public (int H, int K, int L) Index { get; set; }
        /// <summary>面間隔に対する割合で表した平行移動量。261106Cl: 正の値で面が (hkl) の法線の向き (+ĝ) へ動く (ReciPro の Structure Viewer が ĝ·r = (m + Translation) d に描く)</summary>
        public double Translation { get; set; } = 0;
        public double D { get; set; } = 0;
        public bool Enabled { get; set; } = true;

        public (double X, double Y, double Z, double D) PlaneParam { get; set; }

        [XmlIgnore]
        public Color Color => Color.FromArgb(ColorArgb);

        public int ColorArgb { get; set; } = Color.Gray.ToArgb();

        public LatticePlane()
        { }

        public LatticePlane(bool enabled, Crystal crystal, int h, int k, int l, double translation, int argb) : this()
        {
            Enabled = enabled;
            ColorArgb = argb;
            Translation = translation;
            Index = (h, k, l);

            Reset(crystal);
        }

        public void Reset(Crystal crystal)
        {
            var g = Index.H * crystal.A_Star + Index.K * crystal.B_Star + Index.L * crystal.C_Star;
            var d = 1 / g.Length;
            g.NormarizeThis();
            PlaneParam = (g.X, g.Y, g.Z, d);
        }
    }
}