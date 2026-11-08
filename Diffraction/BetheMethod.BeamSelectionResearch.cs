// 261108Cl 新規作成 (研究用): ebsdsim 型のビーム選択 (作者 11-08「ebsdsim の方法を軸にして、他の方法のいいところを取り入れつつ、
//   高精度で高速なビーム選択アルゴリズムを構築しましょう」)。設計 = .project-guidance/ReciPro/ReciPro_ビーム選択の改善_研究計画_261108.md §4.1
//   (codex2 第 14 回)、仕様 = ReciPro_ビーム選択_新方式_仕様書_261108.md。既定の経路 (Find_gVectors、CBED・STEM・PED・HRTEM・ALCHEMI・EBSD で共有) は変えない。
//   研究用の道具だけが EbsdBeamSelectionResearchHook を設定して EBSD の MasterPattern 経路 (RunEBSDNew) の選択を差し替える。
//
//   a0 候補: 計量 |g| ≤ R の逆格子点を完全に列挙する (BFS の個数制限・8N の打ち切り・|s| の窓・接平面の重複除去を持ち込まない)。
//      除くのは格子の centering (純粋な並進) で消える点だけ (映進・らせんの禁制反射は残す)。
//   差の表: U^el_d / γ (d = g − h、|d| ≤ 2R) を整数 3 組の密な配列で持つ (圧縮キーを使わない = 衝突しない)。
//      弾性の U は γ(E)·F/(πV) なので、エネルギーの間は γ で掛け直すだけで共通に使える。
//   結合の尺度 M_g = max_{h ∈ 相手の集合, h ≠ g} |U^el_{g−h}|:
//      E0 (対照) = 相手 = 半径の内側の全プール (ebsdsim の原法)、E1 (主候補) = 000 と |U^el| の上位 50 反射 (境界の同点の組は丸ごと)。
//   点数 r_g = max(|Q_g|, P_g/(πT)) / (2 M_g)、Q_g = k0² − |k0 + g|²、P_g = 2 n·(k0 + g) > 0 (n = 表面の法線、T = 伝播座標の最大)。
//      M = 0 は順位 ∞ (選ばない)。000 は必ず採る。
//   本数: 000 を含む N 本を最低の目安にし、N 本目の点数と同点 (相対 1e-9) の組を全部採る (実本数 ≥ N)。
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Numerics;
using System.Threading.Tasks;

namespace Crystallography;

public partial class BetheMethod
{
    /// <summary>261108Cl 追加 (研究用): ebsdsim 型のビーム選択の、結合の相手の集合。</summary>
    public enum BeamPartnerSet
    {
        /// <summary>E0 (対照): 半径の内側の全プール (ebsdsim の原法。広いプールの内側では尺度がほぼ定数になり、順位が励起誤差だけに退化しうる)</summary>
        E0AllPool,
        /// <summary>E1 (主候補): 000 と、弾性の |U| の大きい上位 K 反射 (境界の同点の組を丸ごと)</summary>
        E1StrongSet,
    }

    /// <summary>261108Cl 追加 (研究用の口。既定 null = 従来の Find_gVectors): EBSD の MasterPattern 経路 (RunEBSDNew) のビーム選択だけを差し替える。
    /// 引数 = (この BetheMethod、屈折後の k0、局所の表面の法線、MaxNumOfBloch)。複数のスレッドから同時に呼ばれる。
    /// プロセス全体に効く静的な値なので、研究用の道具だけが設定する。本体の GUI は使わない。</summary>
    public static Func<BetheMethod, Vector3DBase, Vector3DBase, int, Beam[]> EbsdBeamSelectionResearchHook { get; set; } = null;

    /// <summary>261108Cl 追加 (研究用): 選択 1 回の記録 (<see cref="SelectBeamsResearch"/> に渡すと書き込む)。</summary>
    public sealed class BeamSelectionTrace
    {
        /// <summary>要求の本数 (000 を含む)</summary>
        public int RequestedN;
        /// <summary>実際に採った本数 (000 と境界の同点の組を含む)</summary>
        public int SelectedCount;
        /// <summary>点数が有限の候補 (表面の前方で M &gt; 0) の数</summary>
        public int ScoredCandidates;
        /// <summary>採った最後の点数 (境界)</summary>
        public double ThresholdScore;
        /// <summary>採った反射の最大の |g| [nm⁻¹]</summary>
        public double MaxG;
        /// <summary>採った反射のうち |g| &gt; 0.9 R の数 (半径の打ち切りの診断)</summary>
        public int NearBoundary;
    }

    /// <summary>261108Cl 追加 (研究用): ビーム選択の不変データ (候補の完全なプール・差の表・結合の尺度)。
    /// 材料・占有率・温度因子・模型・半径・相手の数・方位 (BaseRotation) ごとに作る。方向やエネルギー、呼び出しの順・並列度に依らない。</summary>
    public sealed class BeamSelectionPool
    {
        internal BetheMethod Owner;
        internal (double, double, double, double, double, double, double, double, double) Mat;
        /// <summary>候補の半径 R [nm⁻¹]</summary>
        public double RadiusNm { get; internal set; }
        /// <summary>E1 の相手の数 K (同点の組で増えうる)</summary>
        public int StrongPartnerCount { get; internal set; }
        /// <summary>差の表を作った電圧 [kV] (表は U/γ で持つのでエネルギーに依らない)</summary>
        public double ReferenceKV { get; internal set; }
        /// <summary>centering の消滅に使った格子の型 (隠れた centering を足したものは + で続ける)</summary>
        public string CenteringRule { get; internal set; }

        /// <summary>候補の hkl (先頭は 000)</summary>
        public (int h, int k, int l)[] Hkl { get; internal set; }
        /// <summary>候補の逆格子ベクトル (BaseRotation·Crystal.MatrixInverseTransposed·hkl) [nm⁻¹]</summary>
        public double[] GX { get; internal set; }
        public double[] GY { get; internal set; }
        public double[] GZ { get; internal set; }
        public int Count => Hkl.Length;

        /// <summary>候補の hkl の範囲 (|h| ≤ BoundH など)</summary>
        public int BoundH { get; internal set; }
        public int BoundK { get; internal set; }
        public int BoundL { get; internal set; }
        internal int[] poolIndexBox; // (h + BoundH, k + BoundK, l + BoundL) → 候補の番号 (無ければ −1)

        /// <summary>差の表の hkl の範囲 (= 2 × 候補の範囲)</summary>
        public int DiffBoundH { get; internal set; }
        public int DiffBoundK { get; internal set; }
        public int DiffBoundL { get; internal set; }
        internal Complex[] diffU; // (d + DiffBound) → U^el_d / γ [nm⁻²] (|d| > 2R・centering で消える d は 0)
        /// <summary>差の表で値を求めた d の数</summary>
        public int DiffEntries { get; internal set; }
        /// <summary>数値の 0 とみなす |U| の閾値 (d ≠ 0 の最大の |U| × 1e-12)。これ未満は 0 に置いた</summary>
        public double ZeroThreshold { get; internal set; }

        /// <summary>E1 の相手の集合 (000 を含む)</summary>
        public (int h, int k, int l)[] StrongSet { get; internal set; }
        /// <summary>結合の尺度 M_g / γ (E0・E1)。000 は 0</summary>
        public double[] CouplingE0 { get; internal set; }
        public double[] CouplingE1 { get; internal set; }
        /// <summary>M_g を与えた相手の候補の番号 (E0・E1、無ければ −1)</summary>
        public int[] PartnerE0 { get; internal set; }
        public int[] PartnerE1 { get; internal set; }

        /// <summary>準備の時間 [s] (候補・差の表・E1・E0)</summary>
        public double SecondsEnumerate, SecondsDiffTable, SecondsE1, SecondsE0;

        /// <summary>hkl の候補の番号 (無ければ −1)</summary>
        public int IndexOf((int h, int k, int l) x)
        {
            if (Math.Abs(x.h) > BoundH || Math.Abs(x.k) > BoundK || Math.Abs(x.l) > BoundL) return -1;
            return poolIndexBox[((x.h + BoundH) * (2 * BoundK + 1) + x.k + BoundK) * (2 * BoundL + 1) + x.l + BoundL];
        }

        /// <summary>差の表の U^el_d / γ (表の範囲の外は例外)</summary>
        public Complex DifferenceU((int h, int k, int l) d)
        {
            if (Math.Abs(d.h) > DiffBoundH || Math.Abs(d.k) > DiffBoundK || Math.Abs(d.l) > DiffBoundL)
                throw new ArgumentOutOfRangeException(nameof(d), $"d {d} is outside the difference table");
            return diffU[DiffIndex(d.h, d.k, d.l)];
        }

        internal int DiffIndex(int h, int k, int l) => ((h + DiffBoundH) * (2 * DiffBoundK + 1) + k + DiffBoundK) * (2 * DiffBoundL + 1) + l + DiffBoundL;
    }

    /// <summary>261108Cl 追加 (研究用): 弾性の U の相対論の因子 γ (getU と同じ式)。</summary>
    static double ElasticGamma(double kV) => 1 + UniversalConstants.e0 * kV * 1E3 / UniversalConstants.m0 / UniversalConstants.c2;

    /// <summary>261108Cl 追加 (研究用): 格子の centering (純粋な並進) で消える hkl を除く規則。
    /// 結晶の格子の型に加え、includeHidden なら単純格子として与えた結晶の隠れた centering (<see cref="FindHiddenCenteringTranslations"/>) も除く。</summary>
    static (Func<int, int, int, bool> Allowed, string Name) CenteringExtinctionRule(Crystal crystal, bool includeHidden)
    {
        var rules = new List<Func<int, int, int, bool>>();
        var names = new List<string>();
        var type = crystal.Symmetry.LatticeTypeStr;
        void Add(string name)
        {
            names.Add(name);
            rules.Add(name switch
            {
                "F" => static (h, k, l) => ((h ^ k) & 1) == 0 && ((k ^ l) & 1) == 0,
                "I" => static (h, k, l) => ((h + k + l) & 1) == 0,
                "A" => static (h, k, l) => ((k + l) & 1) == 0,
                "B" => static (h, k, l) => ((h + l) & 1) == 0,
                "C" => static (h, k, l) => ((h + k) & 1) == 0,
                "RH" => static (h, k, l) => ((-h + k + l) % 3 + 3) % 3 == 0, // 六方の軸の R (obverse): −h + k + l = 3n (Find_gVectors の directionRH と同じ格子)
                _ => throw new NotSupportedException(name),
            });
        }
        if (type == "F" || type == "I" || type == "A" || type == "B" || type == "C") Add(type);
        else if (type == "R" && crystal.Symmetry.SpaceGroupHMsubStr == "H") Add("RH");
        if (includeHidden)
            foreach (var t in FindHiddenCenteringTranslations(crystal))
                Add(t);
        var array = rules.ToArray();
        return (array.Length == 0 ? static (_, _, _) => true : (h, k, l) => array.All(r => r(h, k, l)),
            names.Count == 0 ? "P" : string.Join("+", names));
    }

    /// <summary>261108Cl 追加 (研究用): ビーム選択の不変データを作る。電圧は今の AccVoltage (0 以下なら 20 kV) で差の表を求め、U/γ で持つ。
    /// 候補 = |g| ≤ radiusNm の逆格子点すべて (centering で消える点を除く)。差の表 = |d| ≤ 2·radiusNm。</summary>
    /// <param name="radiusNm">候補の半径 R [nm⁻¹]</param>
    /// <param name="strongPartnerCount">E1 の相手の数 K (既定 50 = 10 月の K2sum に合わせた研究の定数。物質ごとに調整しない)</param>
    /// <param name="includeHiddenCentering">単純格子として与えた結晶の隠れた centering で消える点も除く</param>
    /// <param name="buildE0">E0 (全プールの max) の尺度も作る</param>
    public BeamSelectionPool CreateBeamSelectionPool(double radiusNm, int strongPartnerCount = 50, bool includeHiddenCentering = true, bool buildE0 = true)
    {
        var sw = Stopwatch.StartNew();
        var rotation = BaseRotation ?? new Matrix3D();
        var mat = rotation * Crystal.MatrixInverseTransposed;
        var (m11, m12, m13, m21, m22, m23, m31, m32, m33) = mat.Tuple;
        var (allowed, ruleName) = CenteringExtinctionRule(Crystal, includeHiddenCentering);
        double kV = AccVoltage > 0 ? AccVoltage : 20;

        // a0: |g| ≤ R の完全な列挙。g = A⁻ᵀh なので h_i = a_i·g (a_i = 実格子の基底) → |h_i| ≤ |a_i|·R (斜交格子でも成り立つ。逆格子の軸の長さで割った上限は斜交格子で取りこぼす)
        int bh = (int)Math.Ceiling(radiusNm * Crystal.A), bk = (int)Math.Ceiling(radiusNm * Crystal.B), bl = (int)Math.Ceiling(radiusNm * Crystal.C);
        double r2 = radiusNm * radiusNm * (1 + 1e-9); // 対称で等価な点の丸めの差で境界の点が割れないよう、相対 1e-9 の余裕 (事前に固定)
        var hkl = new List<(int, int, int)> { (0, 0, 0) };
        var gx = new List<double> { 0 }; var gy = new List<double> { 0 }; var gz = new List<double> { 0 };
        int nk = 2 * bk + 1, nl = 2 * bl + 1;
        var box = new int[(2 * bh + 1) * nk * nl];
        Array.Fill(box, -1);
        box[(bh * nk + bk) * nl + bl] = 0;
        for (int h = -bh; h <= bh; h++)
            for (int k = -bk; k <= bk; k++)
                for (int l = -bl; l <= bl; l++)
                {
                    if ((h, k, l) == (0, 0, 0) || !allowed(h, k, l)) continue;
                    double x = Dot3(m11, h, m12, k, m13, l), y = Dot3(m21, h, m22, k, m23, l), z = Dot3(m31, h, m32, k, m33, l); // Find_gVectors (Build_gCache) と同じ式 = g がビット単位で同じ
                    if (Dot3(x, x, y, y, z, z) > r2) continue;
                    box[((h + bh) * nk + k + bk) * nl + l + bl] = hkl.Count;
                    hkl.Add((h, k, l)); gx.Add(x); gy.Add(y); gz.Add(z);
                }
        var pool = new BeamSelectionPool
        {
            Owner = this, Mat = mat.Tuple, RadiusNm = radiusNm, ReferenceKV = kV, CenteringRule = ruleName,
            Hkl = [.. hkl], GX = [.. gx], GY = [.. gy], GZ = [.. gz], BoundH = bh, BoundK = bk, BoundL = bl, poolIndexBox = box,
            DiffBoundH = 2 * bh, DiffBoundK = 2 * bk, DiffBoundL = 2 * bl,
        };
        pool.SecondsEnumerate = sw.Elapsed.TotalSeconds;

        // 差の表: |d| ≤ 2R の d の U^el / γ。getU (cache: false) をそのまま使う (式を複製しない・uDictionary を膨らませない)
        sw.Restart();
        int dh = pool.DiffBoundH, dk = pool.DiffBoundK, dl = pool.DiffBoundL, dnk = 2 * dk + 1, dnl = 2 * dl + 1;
        var diff = new Complex[(2 * dh + 1) * dnk * dnl];
        double d2max = 4 * r2, gammaRef = ElasticGamma(kV);
        int entries = 0;
        Parallel.For(-dh, dh + 1, () => 0, (h, _, local) =>
        {
            for (int k = -dk; k <= dk; k++)
                for (int l = -dl; l <= dl; l++)
                {
                    if (!allowed(h, k, l)) continue;
                    double x = Dot3(m11, h, m12, k, m13, l), y = Dot3(m21, h, m22, k, m23, l), z = Dot3(m31, h, m32, k, m33, l);
                    if (Dot3(x, x, y, y, z, z) > d2max) continue;
                    diff[((h + dh) * dnk + k + dk) * dnl + l + dl] = getU(kV, new Beam((h, k, l), new Vector3DBase(x, y, z)), cache: false).Real / gammaRef;
                    local++;
                }
            return local;
        }, local => System.Threading.Interlocked.Add(ref entries, local));
        // 数値の 0: d ≠ 0 の最大の |U| の 1e-12 倍未満 (消滅する反射の丸めの残り。事前に固定) を 0 に置く
        int zeroIndex = pool.DiffIndex(0, 0, 0);
        double maxAbs = 0;
        for (int i = 0; i < diff.Length; i++)
            if (i != zeroIndex) maxAbs = Math.Max(maxAbs, diff[i].Magnitude);
        double zeroTol = maxAbs * 1e-12;
        for (int i = 0; i < diff.Length; i++)
            if (diff[i] != Complex.Zero && diff[i].Magnitude < zeroTol) diff[i] = Complex.Zero;
        pool.diffU = diff; pool.DiffEntries = entries; pool.ZeroThreshold = zeroTol;
        pool.SecondsDiffTable = sw.Elapsed.TotalSeconds;

        int n = pool.Count;
        double AbsU(int i, int j) // |U^el_{g_i − g_j}| / γ
        {
            var (h1, k1, l1) = pool.Hkl[i]; var (h2, k2, l2) = pool.Hkl[j];
            return diff[pool.DiffIndex(h1 - h2, k1 - k2, l1 - l2)].Magnitude;
        }

        // E1: 相手 = 000 + |U_g| の上位 K (|U| の降順、同点は hkl の辞書順で並べるだけで、K 番目と同点 (相対 1e-9) の組は丸ごと入れる)
        sw.Restart();
        var order = Enumerable.Range(1, n - 1).Where(i => AbsU(i, 0) > 0)
            .OrderByDescending(i => AbsU(i, 0)).ThenBy(i => pool.Hkl[i]).ToArray();
        int take = Math.Min(strongPartnerCount, order.Length);
        if (take > 0)
        {
            double edge = AbsU(order[take - 1], 0) * (1 - 1e-9);
            while (take < order.Length && AbsU(order[take], 0) >= edge) take++;
        }
        var strong = new int[take + 1];
        strong[0] = 0;
        Array.Copy(order, 0, strong, 1, take);
        pool.StrongSet = [.. strong.Select(i => pool.Hkl[i])];
        pool.StrongPartnerCount = take;
        var c1 = new double[n]; var p1 = new int[n];
        Parallel.For(1, n, i =>
        {
            double best = 0; int partner = -1;
            foreach (var j in strong)
            {
                if (j == i) continue;
                var u = AbsU(i, j);
                if (u > best) { best = u; partner = j; }
            }
            c1[i] = best; p1[i] = partner;
        });
        p1[0] = -1;
        pool.CouplingE1 = c1; pool.PartnerE1 = p1;
        pool.SecondsE1 = sw.Elapsed.TotalSeconds;

        // E0: 相手 = 全プール。差 d (≠ 0) を |U| の降順に調べ、g − d がプールにある最初の d で止める (= max_h |U_{g−h}|、h = g − d ≠ g)
        if (buildE0)
        {
            sw.Restart();
            var ds = new List<(double absU, int h, int k, int l)>();
            for (int h = -dh; h <= dh; h++)
                for (int k = -dk; k <= dk; k++)
                    for (int l = -dl; l <= dl; l++)
                    {
                        if ((h, k, l) == (0, 0, 0)) continue;
                        var u = diff[((h + dh) * dnk + k + dk) * dnl + l + dl].Magnitude;
                        if (u > 0) ds.Add((u, h, k, l));
                    }
            var dsSorted = ds.OrderByDescending(x => x.absU).ThenBy(x => (x.h, x.k, x.l)).ToArray();
            var c0 = new double[n]; var p0 = new int[n];
            Parallel.For(1, n, i =>
            {
                var (h, k, l) = pool.Hkl[i];
                c0[i] = 0; p0[i] = -1;
                foreach (var d in dsSorted)
                {
                    var j = pool.IndexOf((h - d.h, k - d.k, l - d.l));
                    if (j < 0) continue;
                    c0[i] = d.absU; p0[i] = j;
                    break;
                }
            });
            p0[0] = -1;
            pool.CouplingE0 = c0; pool.PartnerE0 = p0;
            pool.SecondsE0 = sw.Elapsed.TotalSeconds;
        }
        return pool;
    }

    /// <summary>261108Cl 追加 (研究用): ebsdsim 型の点数でビームを選ぶ。返り値の先頭は 000、以下は点数の昇順 (同点は hkl の辞書順)。
    /// Beam の U・Q・P は Find_gVectors と同じ作り方 (getU(AccVoltage, hkl, g)、Q = k0² − |k0+g|²、P = 2 n·(k0+g))。</summary>
    /// <param name="pool"><see cref="CreateBeamSelectionPool"/> で作った不変データ (この BetheMethod・同じ方位のもの)</param>
    /// <param name="vecK0">屈折後の k0</param>
    /// <param name="surface">表面の法線 (EBSD のマスターは −出射方向)</param>
    /// <param name="n">000 を含む最低の本数。N 本目と同点の組は全部採る</param>
    /// <param name="partners">結合の相手の集合 (E0 対照 / E1 主候補)</param>
    /// <param name="depthNm">有限深さの正則化の T [nm] (ソルバーで評価する伝播座標の最大。EBSD 170 nm)</param>
    /// <param name="trace">記録 (null なら書かない)</param>
    public Beam[] SelectBeamsResearch(BeamSelectionPool pool, Vector3DBase vecK0, Vector3DBase surface, int n, BeamPartnerSet partners = BeamPartnerSet.E1StrongSet, double depthNm = 170, BeamSelectionTrace trace = null)
    {
        if (pool.Owner != this)
            throw new ArgumentException("The pool was built for another BetheMethod", nameof(pool));
        if (pool.Mat != ((BaseRotation ?? new Matrix3D()) * Crystal.MatrixInverseTransposed).Tuple)
            throw new ArgumentException("BaseRotation changed after the pool was built", nameof(pool));
        var coupling = partners == BeamPartnerSet.E0AllPool
            ? pool.CouplingE0 ?? throw new InvalidOperationException("E0 coupling was not built")
            : pool.CouplingE1;
        double gamma = ElasticGamma(AccVoltage), reg = 1 / (Math.PI * depthNm);
        var (kX, kY, kZ) = vecK0.Tuple;
        var (sX, sY, sZ) = surface.Tuple;
        double k02 = vecK0.Length2; // 以下 Q・P も Find_gVectors と同じ式 (Dot3) = 同じ hkl なら Q・P がビット単位で同じ
        int count = pool.Count;
        var scored = new (double r, int i)[count];
        int m = 0;
        for (int i = 1; i < count; i++)
        {
            double vX = pool.GX[i] + kX, vY = pool.GY[i] + kY, vZ = pool.GZ[i] + kZ;
            double p = 2 * Dot3(sX, vX, sY, vY, sZ, vZ);
            if (p <= 0) continue; // 表面の前方だけ (本体と同じ P > 0)
            double mg = coupling[i] * gamma;
            if (mg <= 0) continue; // M = 0 は順位 ∞
            double q = k02 - Dot3(vX, vX, vY, vY, vZ, vZ);
            scored[m++] = (Math.Max(Math.Abs(q), p * reg) / (2 * mg), i);
        }
        var span = scored.AsSpan(0, m);
        int need = Math.Max(0, n - 1);
        double threshold;
        int selected;
        if (need >= m)
        {
            threshold = double.PositiveInfinity;
            selected = m;
        }
        else if (need == 0)
        {
            threshold = double.NegativeInfinity;
            selected = 0;
        }
        else
        {
            QuickSelect.Execute(span, need, static (a, b) => a.r.CompareTo(b.r));
            threshold = double.NegativeInfinity;
            for (int j = 0; j < need; j++) threshold = Math.Max(threshold, span[j].r);
            // 境界の同点の組 (相対 1e-9、事前に固定) を丸ごと採る: 先頭 need 個の外にある組の残りを前へ寄せる
            double edge = threshold * (1 + 1e-9);
            selected = need;
            for (int j = need; j < m; j++)
                if (span[j].r <= edge) (span[selected], span[j]) = (span[j], span[selected++]);
        }
        var chosen = span[..selected].ToArray();
        Array.Sort(chosen, (a, b) => a.r != b.r ? a.r.CompareTo(b.r) : pool.Hkl[a.i].CompareTo(pool.Hkl[b.i]));

        var beams = new Beam[selected + 1];
        double p0 = 2 * Dot3(sX, kX, sY, kY, sZ, kZ);
        beams[0] = new Beam((0, 0, 0), new Vector3DBase(0, 0, 0), getU(AccVoltage, (0, 0, 0), new Vector3DBase(0, 0, 0)), (k02 - Dot3(kX, kX, kY, kY, kZ, kZ), p0)); // 000 の Q も Find_gVectors と同じ式 (丸めで 0 とは限らない)
        double maxG = 0;
        int nearBoundary = 0;
        for (int j = 0; j < selected; j++)
        {
            int i = chosen[j].i;
            var g = new Vector3DBase(pool.GX[i], pool.GY[i], pool.GZ[i]);
            double vX = g.X + kX, vY = g.Y + kY, vZ = g.Z + kZ;
            double q = k02 - Dot3(vX, vX, vY, vY, vZ, vZ), p = 2 * Dot3(sX, vX, sY, vY, sZ, vZ);
            beams[j + 1] = new Beam(pool.Hkl[i], g, getU(AccVoltage, pool.Hkl[i], g), (q, p));
            var len = g.Length;
            maxG = Math.Max(maxG, len);
            if (len > 0.9 * pool.RadiusNm) nearBoundary++;
        }
        if (trace != null)
        {
            trace.RequestedN = n; trace.SelectedCount = beams.Length; trace.ScoredCandidates = m;
            trace.ThresholdScore = threshold; trace.MaxG = maxG; trace.NearBoundary = nearBoundary;
        }
        return beams;
    }
}
