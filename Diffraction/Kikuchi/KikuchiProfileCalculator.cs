// 260805Cl 新規作成: 菊池線動力学化 Phase 0-1 (設計正本 = ReciPro_菊池線動力学化設計.md §3-§5)。
// 2ビーム (2×2) 複素 Bethe 行列 + 完全な TDS 源密度行列 Q + 既存 EBSD と同型の modal 厚み積分で
// 1D 動力学バンドプロファイル c_g(x) を計算する。
//
// 規約はすべて既存実装から流用する (取り違え防止のため新設しない。設計 §3):
//   - 固有値問題行列: BetheMethod.getEigenMatrix と同一の構成
//       A[row + col*2] = U(g_col − g_row)/P_col (非対角) / (i·U'(0) + Q_col)/P_col (対角)
//   - 厚み積分: EBSDSolverManaged と同一の F_{jj'}(t) = [e^{λt}−1]/λ, λ_{jj'} = 2πi(γ_j − γ_j'*)
//   - 原子位置位相: CreatePhaseFactors と同一の e^{+2πi(h x + k y + l z)}
//   - 反転幾何: EBSD master-pattern と同一 (beamDirection = 検出器方向 d̂, 内向き法線 = −d̂
//       → k0 = −√(k_vac²+u0)·d̂ が自動的に得られる)
//
// ⚠ WinForms / System.Drawing 非依存 (設計 §5)。既存 CBED/STEM/EBSD ソルバと共有 uDictionary には触れない。
//   getU の流用は「菊池専用の BetheMethod インスタンス」経由で行い、その per-instance uDictionary が
//   設計のいう「PotentialSnapshot 内のローカル辞書 (キー = Δg hkl)」に相当する。

using MathNet.Numerics; // Complex.MagnitudeSquared() 拡張 (BetheMethod.cs と同じ)
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;

namespace Crystallography;

#region KikuchiPotentialSnapshot

/// <summary>
/// 結晶・加速電圧・吸収/イオンモデルを固定した不変スナップショット。260805Cl 追加。
/// 共有 Crystal.Bethe とその uDictionary には一切触れない (設計 §5)。
/// </summary>
public sealed class KikuchiPotentialSnapshot
{
    public Crystal Crystal { get; }
    public double KV { get; }
    /// <summary>真空中波数 k_vac [nm⁻¹] (相対論補正込み)</summary>
    public double KVac { get; }
    /// <summary>平均内部ポテンシャル u0 (U(0) 実部)</summary>
    public double U0 { get; }
    /// <summary>吸収対角 U'(0)。AbsorptionOff のとき 0</summary>
    public Complex UPrime0 { get; }
    /// <summary>診断: 吸収 (U' 全成分) を落とす (設計 §6 の E/D 消失テスト用)</summary>
    public bool AbsorptionOff { get; }
    public IKikuchiInelasticKernel Kernel { get; }
    /// <summary>生成時の結晶内容ハッシュ (格子・原子・Occ・ADP・イオンモデル)。鮮度判定用</summary>
    public long ContentHash { get; }

    private readonly BetheMethod _bethe; // 菊池専用インスタンス (ローカル U 辞書)

    /// <summary>対称等価位置展開済みの原子サイト (分率座標, Crystal.Atoms 添字, 占有率)</summary>
    internal readonly (double X, double Y, double Z, int AtomsIndex, double Occ)[] Sites;

    public KikuchiPotentialSnapshot(Crystal crystal, double kV, IKikuchiInelasticKernel kernel = null, bool absorptionOff = false)
    {
        Crystal = crystal;
        KV = kV;
        AbsorptionOff = absorptionOff;
        _bethe = new BetheMethod(crystal);
        KVac = UniversalConstants.Convert.EnergyToElectronWaveNumber(kV);
        var u = _bethe.getU(kV);
        U0 = u.Real.Real;
        UPrime0 = absorptionOff ? Complex.Zero : u.Imag;
        Kernel = kernel ?? new KikuchiTdsEinsteinKernel(crystal);
        ContentHash = ComputeCrystalHash(crystal);

        var sites = new List<(double, double, double, int, double)>();
        for (int i = 0; i < crystal.Atoms.Length; i++)
            foreach (var atom in crystal.Atoms[i].Atom)
                sites.Add((atom.X, atom.Y, atom.Z, i, crystal.Atoms[i].Occ));
        Sites = [.. sites];
    }

    /// <summary>U(hkl)。菊池専用インスタンスのローカル辞書経由 (共有 uDictionary 非使用)</summary>
    public (Complex Real, Complex Imag) GetU(in (int h, int k, int l) index, Vector3DBase vec)
    {
        var u = _bethe.getU(KV, index, vec);
        return AbsorptionOff ? (u.Real, Complex.Zero) : u;
    }

    /// <summary>snapshot が現在の条件に対して新鮮か (結晶内容・kV・イオンモデル)</summary>
    public bool Matches(Crystal crystal, double kV)
        => ReferenceEquals(Crystal, crystal) && KV == kV && ContentHash == ComputeCrystalHash(crystal);

    /// <summary>格子・原子 (Z, イオン, Occ, 位置, Biso)・ElasticIonModel を畳んだ簡易内容ハッシュ</summary>
    public static long ComputeCrystalHash(Crystal c)
    {
        var h = new HashCode();
        h.Add(c.A); h.Add(c.B); h.Add(c.C); h.Add(c.Alpha); h.Add(c.Beta); h.Add(c.Gamma);
        h.Add((int)BetheMethod.ElasticIonModel);
        foreach (var atoms in c.Atoms)
        {
            h.Add(atoms.AtomicNumber); h.Add(atoms.SubNumberElectron); h.Add(atoms.Occ);
            h.Add(atoms.Dsf.Biso);
            foreach (var atom in atoms.Atom)
            { h.Add(atom.X); h.Add(atom.Y); h.Add(atom.Z); }
        }
        return ((long)h.ToHashCode() << 32) ^ (uint)c.Atoms.Length;
    }
}

#endregion

#region 結果型

/// <summary>1 バンド分の 1D 動力学プロファイル。260805Cl 追加</summary>
public sealed class KikuchiBandProfile
{
    public KikuchiBandFamily Family { get; init; }
    /// <summary>sinθ_B (このバンド・この kV)</summary>
    public double SinThetaB { get; init; }
    /// <summary>無次元プロファイル座標 x = sinθ'/sinθ_B (x = ±1 が両 Bragg 線, x = 0 がバンド中心)</summary>
    public double[] X { get; init; }
    /// <summary>無次元コントラスト c_g(x) = (I − I_ref)/max(I_ref, ε) (設計 §4)</summary>
    public double[] C { get; init; }
    /// <summary>動力学順位 C_g = RMS_x(c) (設計 §4 の二段階選定の2段目)</summary>
    public double DynamicalRank { get; init; }
    /// <summary>幾何退化 (バンド法線が視軸とほぼ平行) のとき false</summary>
    public bool Valid { get; init; }

    /// <summary>x での線形補間 (範囲外は 0)</summary>
    public double Interpolate(double x)
    {
        var xs = X;
        if (!Valid || xs.Length < 2 || x <= xs[0] || x >= xs[^1])
            return 0;
        var t = (x - xs[0]) / (xs[1] - xs[0]); // 等間隔格子
        int i = (int)t;
        if (i >= xs.Length - 1)
            return 0;
        var f = t - i;
        return C[i] * (1 - f) + C[i + 1] * f;
    }
}

/// <summary>1 サンプル点の内部量 (単体テスト・golden 保存用の診断出力。設計 §6)。260805Cl 追加</summary>
public struct KikuchiPointDiagnostics
{
    /// <summary>源密度行列 Q (2×2, Hermitian 半正定値)。[i + j*2] = Q_ij</summary>
    public Complex Q00, Q01, Q10, Q11;
    /// <summary>固有値 γ_j</summary>
    public Complex Gamma0, Gamma1;
    /// <summary>励起振幅 α_j = [C⁻¹]_{j,0}</summary>
    public Complex Alpha0, Alpha1;
    /// <summary>固有ベクトル C[g, j] (列 = Bloch 状態)</summary>
    public Complex C00, C10, C01, C11;
    /// <summary>縮約行列 S_{jj'} = α_j α_j'* Σ_{ii'} C_{i,j} C*_{i',j'} Q_{ii'}</summary>
    public Complex S00, S01, S10, S11;
    /// <summary>参照計算 (Bragg 結合 off) の対角 γ と源強度</summary>
    public Complex GammaRef;
    public double SRef;
    public double I, IRef;
}

#endregion

/// <summary>菊池 1D 動力学プロファイル計算 (2ビーム版)。260805Cl 追加</summary>
public static class KikuchiProfileCalculator
{
    /// <summary>v1 の非整合暫定運用タグ (設計 §3「吸収と源の整合」)</summary>
    public const string DisplayNormalizedTag = "Display-normalized / source-loss not balanced";

    public sealed class Options
    {
        /// <summary>x 格子点数 (奇数推奨: 中心 x=0 を含む)</summary>
        public int SampleCount { get; init; } = 129;
        /// <summary>|x| の上限 (Bragg 半幅の倍数)</summary>
        public double XMax { get; init; } = 2.5;
        /// <summary>診断: 2 波の源振幅 τ を入れ替える (E/D 交換テスト。設計 §6)</summary>
        public bool SwapSourceWeights { get; init; }
        /// <summary>診断: Bragg 結合 U_g を 0 にする (平坦テスト。設計 §6)</summary>
        public bool DisableBraggCoupling { get; init; }
        /// <summary>診断: 源コヒーレンス (Q の非対角) を 0 にする対角近似 (設計 §3: 診断・単体テスト用)</summary>
        public bool DiagonalSourceOnly { get; init; }
        /// <summary>
        /// バンド最近接点 (x=0 の出射方向) の視軸からの最大散乱角 [rad]。これを超えるバンドは Valid=false。
        /// 前方ピーク源に対し I_ref が潰れる遠方バンドの発散と、表示域外バンドの無駄な計算を防ぐ
        /// </summary>
        public double MaxScatteringAngle { get; init; } = 0.35;
        /// <summary>
        /// 260805Cl 追加: プロファイル代表点の視軸からの最小サンプル角 [rad]。
        /// 最近接点をそのまま使うと視軸をかすめるバンドで q_0 → 0 (TDS 源消失) となり
        /// c = I/I_ref が発散して帯全長を飽和させる (Si [111] 実写で確認)。代表点をバンドに沿って
        /// この角度まで滑らせて評価する (Omoto 比較も中心を外した行で行う流儀と同じ)。
        /// </summary>
        public double SampleAngle { get; init; } = 0.025;
    }

    #region 公開 API

    /// <summary>
    /// 1 バンドの 1D プロファイルを計算する。
    /// rotation = 結晶回転 (lab: ビームは −z へ進む)、thickness [nm]。
    /// </summary>
    public static KikuchiBandProfile ComputeProfile(KikuchiPotentialSnapshot snap, KikuchiBandFamily family,
        Matrix3D rotation, double thickness, Options opt = null)
    {
        opt ??= new Options();
        int n = Math.Max(3, opt.SampleCount);
        var xs = new double[n];
        var cs = new double[n];
        for (int i = 0; i < n; i++)
            xs[i] = -opt.XMax + 2 * opt.XMax * i / (n - 1);

        var geo = new BandGeometry(snap, family, rotation);
        // バンド最近接点が視軸から遠すぎる族は除外 (d̂(0)·b̂ = 視軸との余弦)
        if (geo.Valid && -geo.Direction(0).Z < Math.Cos(opt.MaxScatteringAngle))
            return new KikuchiBandProfile { Family = family, SinThetaB = geo.SinThetaB, X = xs, C = cs, Valid = false };
        if (!geo.Valid)
            return new KikuchiBandProfile { Family = family, SinThetaB = geo.SinThetaB, X = xs, C = cs, Valid = false };

        // 1 族 = ±g 両メンバーの 2×2 の和。検出方向が x を掃くとき、+g 系と −g 系の共鳴が
        // それぞれ片側の Bragg 線に現れるため、両方を足して初めてバンド全体になる。
        // E/D 非対称は member 間の鏡映を破る「前方ピーク源」(τ の q 依存) だけから生じ、
        // 一様源+吸収なしでは c(x) = c(−x) が厳密に成り立つ (設計 §1-3 と整合)。
        //
        // 260805Cl 正規化の見直し: 晶帯軸至近では q_0 → 0 で TDS 源が物理的に消え I_ref → 0 になる
        // (比 c が発散し float Inf → 描画 NaN 全透明、の実機事故)。設計 §4 の ε は
        // 「プロファイル内の max I_ref の 1e-3」という相対フロアとして実装する。
        var iArr = new double[n];
        var irefArr = new double[n];
        double maxIref = 0;
        for (int i = 0; i < n; i++)
        {
            var dP = ComputePoint(snap, geo, xs[i], thickness, opt, +1);
            var dM = ComputePoint(snap, geo, xs[i], thickness, opt, -1);
            iArr[i] = dP.I + dM.I;
            irefArr[i] = dP.IRef; // I_ref は直進チャネルのみで member 非依存 (dM.IRef と同値)
            maxIref = Math.Max(maxIref, irefArr[i]);
        }
        double sumSq = 0;
        var floor = 1e-3 * maxIref;
        for (int i = 0; i < n; i++)
        {
            var c = thickness <= 0 || floor <= 0 ? 0 : (iArr[i] - 2 * irefArr[i]) / Math.Max(irefArr[i], floor);
            if (!double.IsFinite(c))
                c = 0; // 非有限値は描画層に流さない
            cs[i] = c;
            sumSq += c * c;
        }
        return new KikuchiBandProfile
        {
            Family = family,
            SinThetaB = geo.SinThetaB,
            X = xs,
            C = cs,
            DynamicalRank = Math.Sqrt(sumSq / n),
            Valid = true,
        };
    }

    /// <summary>1 サンプル点・1 メンバー (member = ±1 = ±g 系) の内部量を返す (単体テスト・golden 用)</summary>
    public static KikuchiPointDiagnostics ComputePointDiagnostics(KikuchiPotentialSnapshot snap, KikuchiBandFamily family,
        Matrix3D rotation, double thickness, double x, int member = +1, Options opt = null)
    {
        opt ??= new Options();
        var geo = new BandGeometry(snap, family, rotation);
        if (!geo.Valid)
            throw new InvalidOperationException("band geometry degenerate (g // beam axis)");
        return ComputePoint(snap, geo, x, thickness, opt, member);
    }

    /// <summary>
    /// 候補族すべてのプロファイルを計算し、動力学順位 C_g = RMS(c) の上位 topN を返す (設計 §4 の二段階選定)。
    /// previousSelection には前回選ばれた族の Index を渡すと順位に 15% のヒステリシスが掛かる
    /// (厚み変更などでの順位境界の点滅防止。設計 §4: 10–20%)。
    /// </summary>
    public static List<KikuchiBandProfile> ComputeProfiles(KikuchiPotentialSnapshot snap, IReadOnlyList<KikuchiBandFamily> candidates,
        Matrix3D rotation, double thickness, int topN, Options opt = null, ISet<(int, int, int)> previousSelection = null)
    {
        opt ??= new Options();
        var profiles = new KikuchiBandProfile[candidates.Count];
        System.Threading.Tasks.Parallel.For(0, candidates.Count, i =>
            profiles[i] = ComputeProfile(snap, candidates[i], rotation, thickness, opt));

        double Rank(KikuchiBandProfile p) =>
            p.DynamicalRank * (previousSelection != null && previousSelection.Contains(p.Family.Index) ? 1.15 : 1.0);

        return [.. profiles.Where(p => p.Valid && p.DynamicalRank > 0)
                           .OrderByDescending(Rank)
                           .Take(topN)];
    }

    #endregion

    #region バンド幾何

    /// <summary>1 バンドの回転済み幾何 (プロファイル計算の共有前計算)</summary>
    private readonly struct BandGeometry
    {
        public readonly KikuchiBandFamily Family;
        public readonly bool Valid;
        public readonly double SinThetaB;
        public readonly Vector3DBase GLab;      // 回転済み g (nm⁻¹)
        public readonly Vector3DBase GHat;      // 単位 g
        public readonly Vector3DBase GPerp;     // ĝ のビーム垂直成分 (単位)
        public readonly double CosEps;          // ĝ·b̂ (b̂ = (0,0,−1))
        public readonly double SinEps;          // √(1−CosEps²)

        public BandGeometry(KikuchiPotentialSnapshot snap, KikuchiBandFamily family, Matrix3D rotation)
        {
            Family = family;
            GLab = rotation * family.Vec;
            var len = GLab.Length;
            SinThetaB = len / (2 * snap.KVac);
            GHat = GLab / len;
            CosEps = -GHat.Z; // b̂ = (0,0,−1) との内積
            var sp2 = 1 - CosEps * CosEps;
            SinEps = sp2 > 0 ? Math.Sqrt(sp2) : 0;
            Valid = SinEps > 1e-6 && SinThetaB < 0.5; // 退化 (バンド法線 ∥ 視軸) と非物理な広角を除外
            GPerp = Valid ? (GHat - CosEps * new Vector3DBase(0, 0, -1)) / SinEps : new Vector3DBase(1, 0, 0);
        }

        /// <summary>
        /// プロファイル座標 x (= sinθ'/sinθ_B) に対応する出射方向 d̂ を返す。sinθ' ≡ −ĝ·d̂ (x = +1 が +g の Bragg 線)。
        /// 既定はバンド上で視軸に最近接の点 (b̂ と ĝ の張る面内の解) だが、最近接点が sampleAngle より
        /// 視軸に近い場合は、同じ円錐 (−ĝ·d̂ = sinθ') 上を面外へ滑らせて視軸から sampleAngle の点で評価する
        /// (q_0 → 0 の発散防止。260805Cl 追加)。
        /// </summary>
        public Vector3DBase Direction(double x, double sampleAngle = 0)
        {
            var sinTp = Math.Clamp(x * SinThetaB, -0.99, 0.99);
            var bHat = new Vector3DBase(0, 0, -1);
            // d̂ = cosψ·b̂ + sinψ·ĝp,  −ĝ·d̂ = sinθ' → cos(ψ − ψ0) = sinθ', cosψ0 = −CosEps, sinψ0 = −SinEps
            var psi0 = Math.Atan2(-SinEps, -CosEps);
            var dPsi = Math.Acos(sinTp);
            double p1 = psi0 + dPsi, p2 = psi0 - dPsi;
            var psi = Math.Cos(p1) >= Math.Cos(p2) ? p1 : p2; // 視軸 (cosψ 最大) に近い側
            var (sinTs, cosTs) = Math.SinCos(sampleAngle);
            if (sampleAngle > 0 && Math.Cos(psi) > cosTs)
            {
                // 円錐上で視軸から sampleAngle の点: d̂ = cosθs·b̂ + sinθs·(cosφ·ĝp + sinφ·(b̂×ĝp))
                // −ĝ·d̂ = sinθ' → cosφ = −(sinθ' + CosEps·cosθs)/(SinEps·sinθs)
                var cosPhi = -(sinTp + CosEps * cosTs) / (SinEps * sinTs);
                if (Math.Abs(cosPhi) <= 1)
                {
                    var sinPhi = Math.Sqrt(1 - cosPhi * cosPhi); // 面外側は対称なので正側を取る
                    var bxg = new Vector3DBase(bHat.Y * GPerp.Z - bHat.Z * GPerp.Y, bHat.Z * GPerp.X - bHat.X * GPerp.Z, bHat.X * GPerp.Y - bHat.Y * GPerp.X);
                    return cosTs * bHat + sinTs * (cosPhi * GPerp + sinPhi * bxg);
                }
            }
            var (sinPsi, cosPsi) = Math.SinCos(psi);
            return cosPsi * bHat + sinPsi * GPerp;
        }
    }

    #endregion

    #region 1 点計算 (2×2 EVD + 完全 Q + modal 積分)

    /// <summary>member = +1 で {0, +g} 系、−1 で {0, −g} 系の 2 ビーム計算 (バンドは両者の和)</summary>
    private static KikuchiPointDiagnostics ComputePoint(KikuchiPotentialSnapshot snap, in BandGeometry geo,
        double x, double thickness, Options opt, int member)
    {
        var kVac = snap.KVac;
        var dHat = geo.Direction(x, opt.SampleAngle); //260805Cl: 代表点は視軸から SampleAngle 以上 (q_0→0 の発散防止)
        var bHat = new Vector3DBase(0, 0, -1);
        var gm = member * geo.GLab; // member 側の g
        var (mh, mk, ml) = (member * geo.Family.Index.H, member * geo.Family.Index.K, member * geo.Family.Index.L);

        // --- 源の運動量移行。反転定理の対応「反転波の beam h ↔ 出射 beam −h」より、
        //     反転波 beam 0 には q_0 = k_vac(d̂ − b̂)、beam g_m には q_1 = q_0 − g_m が対になる ---
        var q0 = kVac * (dHat - bHat);
        var q1 = q0 - gm;
        double s20 = q0.Length2 / 4, s21 = q1.Length2 / 4;

        // --- 反転幾何 (EBSD master-pattern と同一): beamDirection = d̂, 内向き法線 n̂ = −d̂ ---
        var nHat = -dHat;
        var b = nHat * (kVac * dHat); // = −k_vac
        var xr = Math.Sqrt(b * b + snap.U0) - b;
        var k0 = kVac * dHat + xr * nHat; // = −√(k_vac²+u0)·d̂

        // --- 2 ビーム {0, g_m} の Q, P (getQ / getP と同一式) ---
        var k0g = k0 + gm;
        double bigQ1 = k0.Length2 - k0g.Length2;
        double p0 = 2 * (nHat * k0);
        double p1 = 2 * (nHat * k0g);

        // --- 固有値問題行列 (getEigenMatrix と同一の構成) ---
        var uG = snap.GetU((mh, mk, ml), gm);
        var uMg = snap.GetU((-mh, -mk, -ml), -gm);
        var diagU = Complex.ImaginaryOne * snap.UPrime0;
        Complex a00 = diagU / p0;
        Complex a11 = (diagU + bigQ1) / p1;
        Complex a01 = opt.DisableBraggCoupling ? Complex.Zero : (uG.Real + Complex.ImaginaryOne * uG.Imag) / p1;   // [row0, col1] = U(g_1−g_0) = U(+g)
        Complex a10 = opt.DisableBraggCoupling ? Complex.Zero : (uMg.Real + Complex.ImaginaryOne * uMg.Imag) / p0; // [row1, col0] = U(g_0−g_1) = U(−g)

        // --- 2×2 EVD (閉形式) と α = C⁻¹ψ0, ψ0 = (1,0)ᵀ ---
        var half = (a00 + a11) / 2;
        var disc = Complex.Sqrt((a00 - a11) * (a00 - a11) / 4 + a01 * a10);
        Complex g0 = half + disc, g1 = half - disc;
        Complex c00, c10, c01, c11;
        if (a01.Magnitude < 1e-300 && a10.Magnitude < 1e-300)
        { c00 = 1; c10 = 0; c01 = 0; c11 = 1; g0 = a00; g1 = a11; } // 既に対角
        else if (a01.Magnitude >= a10.Magnitude)
        { c00 = a01; c10 = g0 - a00; c01 = a01; c11 = g1 - a00; }
        else
        { c00 = g0 - a11; c10 = a10; c01 = g1 - a11; c11 = a10; }
        // 列正規化 (α と S では列スケールは相殺するが、数値安定のため)
        var n0 = Math.Sqrt(c00.MagnitudeSquared() + c10.MagnitudeSquared());
        var n1 = Math.Sqrt(c01.MagnitudeSquared() + c11.MagnitudeSquared());
        if (n0 > 0) { c00 /= n0; c10 /= n0; }
        if (n1 > 0) { c01 /= n1; c11 /= n1; }
        var det = c00 * c11 - c01 * c10;
        Complex al0, al1;
        if (det.Magnitude < 1e-14)
        { c00 = 1; c10 = 0; c01 = 0; c11 = 1; al0 = 1; al1 = 0; } // 縮退 (defective) は素通し扱い
        else
        { al0 = c11 / det; al1 = -c10 / det; }

        // --- 源振幅 τ とコヒーレンス (原子種ごと・2 波ぶん)。診断 swap は s² の入れ替えとして両方に効かせる ---
        var atoms = snap.Crystal.Atoms;
        int nSpec = atoms.Length;
        var (sA, sB) = opt.SwapSourceWeights ? (s21, s20) : (s20, s21);
        var s2g = gm.Length2 / 4; // |q_0 − q_1|²/4 = |g|²/4 (2 ビームでは固定)
        Span<double> tau0 = stackalloc double[nSpec];
        Span<double> tau1 = stackalloc double[nSpec];
        Span<double> coh = stackalloc double[nSpec];
        for (int i = 0; i < nSpec; i++)
        {
            tau0[i] = snap.Kernel.SourceAmplitude(i, sA);
            tau1[i] = snap.Kernel.SourceAmplitude(i, sB);
            // 260805Cl: 非対角は厳密な Einstein 混合動的形状因子 (Omoto 2002 eq. 38)。因子化近似は既定実装側
            coh[i] = snap.Kernel.SourceCoherence(i, sA, sB, s2g);
        }

        // --- 完全な源密度行列 Q_ij = Σ_a (設計 §3。v1 既定は完全 Q) ---
        //     Q_00 = Σ Occ·τ0², Q_11 = Σ Occ·τ1², Q_01 = Σ Occ·coh·e^{−2πi g_m·r_a} (CreatePhaseFactors と同じ符号規約)
        double q00 = 0, q11 = 0;
        Complex q01 = Complex.Zero;
        foreach (var site in snap.Sites)
        {
            var occ = site.Occ;
            double t0 = tau0[site.AtomsIndex], t1 = tau1[site.AtomsIndex];
            q00 += occ * t0 * t0;
            q11 += occ * t1 * t1;
            var (sin, cos) = Math.SinCos(2 * Math.PI * (mh * site.X + mk * site.Y + ml * site.Z));
            q01 += occ * coh[site.AtomsIndex] * new Complex(cos, -sin);
        }
        if (opt.DiagonalSourceOnly)
            q01 = Complex.Zero;
        var q10 = Complex.Conjugate(q01);

        // --- S_{jj'} = α_j α_j'* Σ_{ii'} C_{i,j} C*_{i',j'} Q_{ii'} ---
        Complex W(Complex ci0, Complex ci1, Complex cj0, Complex cj1) // Σ_{ii'} C_{i,j} C*_{i',j'} Q_{ii'}
            => ci0 * Complex.Conjugate(cj0) * q00 + ci0 * Complex.Conjugate(cj1) * q01
             + ci1 * Complex.Conjugate(cj0) * q10 + ci1 * Complex.Conjugate(cj1) * q11;
        var s00 = al0 * Complex.Conjugate(al0) * W(c00, c10, c00, c10);
        var s01 = al0 * Complex.Conjugate(al1) * W(c00, c10, c01, c11);
        var s10 = al1 * Complex.Conjugate(al0) * W(c01, c11, c00, c10);
        var s11 = al1 * Complex.Conjugate(al1) * W(c01, c11, c01, c11);

        // --- modal 厚み積分 I(t) = Re Σ S_{jj'} F_{jj'}(t) (EBSDSolverManaged と同一) ---
        var intensity = (Fjj(g0, g0, thickness) * s00 + Fjj(g0, g1, thickness) * s01
                       + Fjj(g1, g0, thickness) * s10 + Fjj(g1, g1, thickness) * s11).Real;

        // --- 参照強度 I_ref: 同一源で Bragg 結合を切った計算 (設計 §4)。結合 0 では検出器へ届くのは
        //     直接波チャネルのみ → I_ref = Q_00 · F_00(γ_ref), γ_ref = A00 (吸収対角は残す)。
        //     260805Cl: swap 診断はバンド側 (Q) だけを反転し、参照 = 実背景は入れ替えない
        //     (swap 時は τ(q_0) が tau1 側に入っているので q11 がそのまま実背景の Q_00 になる) ---
        var q00ref = opt.SwapSourceWeights ? q11 : q00;
        var gRef = diagU / p0;
        var iRef = (Fjj(gRef, gRef, thickness) * q00ref).Real;

        return new KikuchiPointDiagnostics
        {
            Q00 = q00, Q01 = q01, Q10 = q10, Q11 = q11,
            Gamma0 = g0, Gamma1 = g1, Alpha0 = al0, Alpha1 = al1,
            C00 = c00, C10 = c10, C01 = c01, C11 = c11,
            S00 = s00, S01 = s01, S10 = s10, S11 = s11,
            GammaRef = gRef, SRef = q00ref,
            I = Math.Max(0, intensity), IRef = Math.Max(0, iRef),
        };
    }

    /// <summary>F_{jj'}(t) = [e^{λt}−1]/λ, λ = 2πi(γ_j − γ_j'*)。λ≈0 は t (EBSDSolverManaged と同一)</summary>
    private static Complex Fjj(in Complex gammaJ, in Complex gammaJp, double t)
    {
        var lam = 2 * Math.PI * Complex.ImaginaryOne * (gammaJ - Complex.Conjugate(gammaJp));
        return lam.MagnitudeSquared() < 1e-30 ? t : (Complex.Exp(lam * t) - Complex.One) / lam;
    }

    #endregion
}

#region 物理プロファイルキャッシュ (Phase 0 骨格)

/// <summary>
/// 物理プロファイルの二層キャッシュのうち物理側 (設計 §5)。
/// キー: 結晶内容ハッシュ / kV / 厚み / 族 / バンド傾き cosε (1e-3 量子化。プロファイルは方位のうち
/// 傾きのみに依存し、視軸まわりの方位角には不変) / 格子仕様 / 診断フラグ。
/// 260805Cl 追加
/// </summary>
public sealed class KikuchiProfileCache
{
    private readonly Dictionary<(long Hash, double KV, double T, (int, int, int) Index, int CosEpsQ, int Samples, double XMax, bool Abs), KikuchiBandProfile> _dic = [];
    private const int Capacity = 1024;

    public bool TryGet(KikuchiPotentialSnapshot snap, KikuchiBandFamily family, Matrix3D rotation, double thickness,
        KikuchiProfileCalculator.Options opt, out KikuchiBandProfile profile)
        => _dic.TryGetValue(Key(snap, family, rotation, thickness, opt), out profile);

    public void Add(KikuchiPotentialSnapshot snap, KikuchiBandFamily family, Matrix3D rotation, double thickness,
        KikuchiProfileCalculator.Options opt, KikuchiBandProfile profile)
    {
        if (_dic.Count >= Capacity)
            _dic.Clear(); // 素朴な全消し (LRU は将来必要になったら)
        _dic[Key(snap, family, rotation, thickness, opt)] = profile;
    }

    public void Clear() => _dic.Clear();

    private static (long, double, double, (int, int, int), int, int, double, bool) Key(
        KikuchiPotentialSnapshot snap, KikuchiBandFamily family, Matrix3D rotation, double thickness, KikuchiProfileCalculator.Options opt)
    {
        var gLab = rotation * family.Vec;
        var cosEps = -gLab.Z / gLab.Length;
        return (snap.ContentHash, snap.KV, thickness, family.Index, (int)Math.Round(cosEps * 1000), opt.SampleCount, opt.XMax, snap.AbsorptionOff);
    }
}

#endregion
