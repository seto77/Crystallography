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
    /// <summary>対象の結晶 (参照のみ保持。数値は生成時に凍結済み)</summary>
    public Crystal Crystal { get; }
    /// <summary>加速電圧 [kV]</summary>
    public double KV { get; }
    /// <summary>真空中波数 k_vac [nm⁻¹] (相対論補正込み)</summary>
    public double KVac { get; }
    /// <summary>平均内部ポテンシャル u0 (U(0) 実部)</summary>
    public double U0 { get; }
    /// <summary>吸収対角 U'(0)。AbsorptionOff のとき 0</summary>
    public Complex UPrime0 { get; }
    /// <summary>診断: 吸収 (U' 全成分) を落とす (設計 §6 の E/D 消失テスト用)</summary>
    public bool AbsorptionOff { get; }
    /// <summary>TDS 源カーネル (既定 = 厳密 Einstein 混合形)</summary>
    public IKikuchiInelasticKernel Kernel { get; }
    /// <summary>生成時の結晶内容ハッシュ (格子・原子・Occ・ADP・イオンモデル)。鮮度判定用</summary>
    public long ContentHash { get; }
    /// <summary>生成時の Crystal.Atoms.Length (260806Cl /simplify2: hot loop から live crystal を読まないための凍結値。
    /// Sites の AtomsIndex はこの値を上限とするため、計算中に UI が原子を削除しても配列境界が壊れない)</summary>
    public int SpeciesCount { get; }

    private readonly BetheMethod _bethe; // 菊池専用インスタンス (ローカル U 辞書)

    /// <summary>対称等価位置展開済みの原子サイト (分率座標, Crystal.Atoms 添字, 占有率)</summary>
    internal readonly (double X, double Y, double Z, int AtomsIndex, double Occ)[] Sites;

    /// <summary>snapshot を生成し、原子サイト・カーネル・U(0)・内容ハッシュを凍結する (UI スレッドで呼ぶこと。260806Cl /simplify2)</summary>
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
        SpeciesCount = crystal.Atoms.Length;

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

    /// <summary>
    /// snapshot が「標準用途」の現在条件に対して新鮮か (結晶内容・kV・イオンモデル)。
    /// 260806Cl /simplify: 診断用 snapshot (AbsorptionOff) は常に false (GUI が診断 snapshot を誤って再利用しない保険)
    /// </summary>
    public bool Matches(Crystal crystal, double kV)
        => !AbsorptionOff && ReferenceEquals(Crystal, crystal) && KV == kV && ContentHash == ComputeCrystalHash(crystal);

    /// <summary>格子・原子 (Z, イオン, Occ, 位置, ADP)・ElasticIonModel を畳んだ簡易内容ハッシュ</summary>
    public static long ComputeCrystalHash(Crystal c)
    {
        var h = new HashCode();
        h.Add(c.A); h.Add(c.B); h.Add(c.C); h.Add(c.Alpha); h.Add(c.Beta); h.Add(c.Gamma);
        h.Add((int)BetheMethod.ElasticIonModel);
        foreach (var atoms in c.Atoms)
        {
            h.Add(atoms.AtomicNumber); h.Add(atoms.SubNumberElectron); h.Add(atoms.Occ);
            var dsf = atoms.Dsf;
            h.Add(dsf.Biso);
            //260806Cl /simplify2: 異方 ADP 編集でも snapshot が失効するようフィールドを網羅 (kernel は Biso000 へ fallback するため)
            h.Add(dsf.Biso000); h.Add(dsf.B11); h.Add(dsf.B22); h.Add(dsf.B33); h.Add(dsf.B12); h.Add(dsf.B23); h.Add(dsf.B31);
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
    /// <summary>このプロファイルが属する g/−g 面族</summary>
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

    //260806Cl /simplify2 削除: Interpolate(double)。renderer が補間をインライン化して以降呼び出しゼロで、
    //正しさに敏感なインデックス計算の複製 (ドリフト危険) になっていた。必要になったら renderer 側の実装を移す。
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
    /// <summary>参照計算 (Bragg 結合 off) の対角 γ</summary>
    public Complex GammaRef;
    /// <summary>参照計算の源強度 (実背景の Q_00)</summary>
    public double SRef;
    /// <summary>この点の強度 I と参照強度 I_ref (非有限値は 0 に正規化済み。260806Cl /simplify2)</summary>
    public double I, IRef;
}

#endregion

/// <summary>菊池 1D 動力学プロファイル計算 (2ビーム版)。260805Cl 追加</summary>
public static class KikuchiProfileCalculator
{
    /// <summary>
    /// v1 の非整合暫定運用タグ (設計 §3「吸収と源の整合」の物理部分)。
    /// 260806Cl: 表示スケール (linear/log/tanh) が選択式になったため、スケール名は GUI 側で前置し
    /// 本定数は物理の注記のみに戻した (作者指摘による開示の実装)
    /// </summary>
    public const string DisplayNormalizedTag = "source-loss not balanced";

    /// <summary>プロファイル計算の格子仕様・診断フラグ (作者調整枠は設計 §9 参照)</summary>
    public sealed class Options
    {
        /// <summary>SampleCount の既定値 (GUI の品質プリセット "Standard" と共有。260806Cl /simplify)</summary>
        public const int DefaultSampleCount = 129;

        /// <summary>x 格子点数 (奇数推奨: 中心 x=0 を含む)</summary>
        public int SampleCount { get; init; } = DefaultSampleCount;
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
        /// 260807Cl 追加 (設計 Phase 2.5): systematic row の片側次数 N。
        /// ビーム = {n·g : n = −N..+N} の **2N+1 波**を 1 つの基底で解く。
        /// **0 (既定) = 従来の 2 ビーム経路** — {0,+g} と {0,−g} を別々に 2×2 で解いて和を取る。
        /// 既定のままなら数値は 1 ビットも変わらない (新経路は完全に分岐している)。
        /// ⚠ 2 ビーム経路が「1 族 = ±g 2 メンバーの和」で直進チャネルを 2 本数えるのに対し、
        /// row 経路は ±g が同じ基底に入るので直進チャネルは 1 本。c の参照強度の係数が違う
        /// (ComputeProfile の refCount)。
        /// </summary>
        public int RowOrder { get; init; }

        /// <summary>RowOrder の上限 (2N+1 = 31 波)。これを超える値は丸める</summary>
        public const int MaxRowOrder = 15;

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
        // 幾何退化、またはバンド最近接点が視軸から遠すぎる族は除外 (d̂(0)·b̂ = 視軸との余弦)。
        // 260806Cl /simplify: 同一の invalid プロファイルを返す 2 つの if を 1 つに統合 (|| は短絡するので退化時 Direction は呼ばれない)
        if (!geo.Valid || -geo.Direction(0).Z < Math.Cos(opt.MaxScatteringAngle))
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
        //260807Cl: row 経路は ±g が同じ基底に入るので直進チャネルが 1 本。2 ビーム経路は
        //±g 2 メンバーの和なので 2 本。差し引く参照強度の本数がここで決まる
        int refCount = opt.RowOrder > 0 ? 1 : 2;
        var row = opt.RowOrder > 0 ? new RowPotential(snap, geo, opt.RowOrder) : null;
        for (int i = 0; i < n; i++)
        {
            if (row != null)
            {
                var d = ComputeRowPoint(snap, geo, row, xs[i], thickness, opt);
                iArr[i] = d.I;
                irefArr[i] = d.IRef;
            }
            else
            {
                var dP = ComputePoint(snap, geo, xs[i], thickness, opt, +1);
                var dM = ComputePoint(snap, geo, xs[i], thickness, opt, -1);
                iArr[i] = dP.I + dM.I;
                irefArr[i] = dP.IRef; // I_ref は直進チャネルのみで member 非依存 (dM.IRef と同値)
            }
            maxIref = Math.Max(maxIref, irefArr[i]);
        }
        double sumSq = 0;
        var floor = 1e-3 * maxIref;
        for (int i = 0; i < n; i++)
        {
            var c = thickness <= 0 || floor <= 0 ? 0 : (iArr[i] - refCount * irefArr[i]) / Math.Max(irefArr[i], floor);
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
        //260807Cl (設計 Phase 2.5): row 計算を使うときは候補を row の生成元へ畳む。
        //畳まないと {020} の row が beam として含む {040} を族としても足してしまい二重計上になる
        //(MgO 比較で corr 0.914 → 0.898 と悪化するのを実測済み)。忘れると静かに悪化するので自動適用する
        if (opt.RowOrder > 0)
            candidates = KikuchiBandFamily.CollapseSystematicRows(candidates, opt.XMax);
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
            //260806Cl /simplify2 (数値レビュー F-3): √(1−cos²) は cosε≈±1 で桁落ちする (guard 境界で相対誤差 ~1e-4)。
            //面内成分 √(Gx²+Gy²) は同じ量を桁落ちなしで与え、GPerp = (Gx, Gy, 0)/SinEps が 1 ulp 単位で正規化される
            SinEps = Math.Sqrt(GHat.X * GHat.X + GHat.Y * GHat.Y);
            //⚠不変条件 (F-2): SinThetaB < 0.5 と Direction の Clamp(±0.99) の組で p1 = 2K + 2m|g|sinθ' > 0.02·k_vac が保証される。
            //どちらかを緩めるときは ComputePoint の p1 ガードが効くことを確認すること
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
                    var bxg = Vector3DBase.VectorProduct(bHat, GPerp); //260806Cl /simplify: 手書き外積を既存ヘルパーへ
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

        // --- 反転幾何 (EBSD master-pattern と同一): beamDirection = d̂, 内向き法線 n̂ = −d̂。
        //     getVecK0 の二次方程式は n̂ = −d̂ (単位) では閉形式 k0 = −√(k_vac²+u0)·d̂ に潰れる (260806Cl /simplify) ---
        var nHat = -dHat;
        var k0 = -Math.Sqrt(kVac * kVac + snap.U0) * dHat;

        // --- 2 ビーム {0, g_m} の Q, P (getQ / getP と同一式) ---
        var k0g = k0 + gm;
        double bigQ1 = k0.Length2 - k0g.Length2;
        double p0 = 2 * (nHat * k0);
        double p1 = 2 * (nHat * k0g);
        //260806Cl /simplify2 (数値レビュー F-2): p1 → 0 は現状 SinThetaB<0.5 + Clamp±0.99 の組合せで排除されているが、
        //その保証を暗黙にしない。ほぼ掠め角のビームは物理的にも 2 波近似の外なので、その点は寄与 0 で返す
        if (p1 <= 1e-3 * p0)
            return default;

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
        //260806Cl /simplify2 (C-1): 種数は snapshot の凍結値を使う (live crystal を並列 hot loop から読まない。
        //計算中の原子削除で Sites.AtomsIndex が配列境界を越える事故の防止)。
        //(F-4): stackalloc は上限付き (巨大セルで pool スレッドの StackOverflow を防ぐ)
        int nSpec = snap.SpeciesCount;
        var (sA, sB) = opt.SwapSourceWeights ? (s21, s20) : (s20, s21);
        var s2g = gm.Length2 / 4; // |q_0 − q_1|²/4 = |g|²/4 (2 ビームでは固定)
        Span<double> tau0 = nSpec <= 64 ? stackalloc double[64] : new double[nSpec];
        Span<double> tau1 = nSpec <= 64 ? stackalloc double[64] : new double[nSpec];
        Span<double> coh = nSpec <= 64 ? stackalloc double[64] : new double[nSpec];
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
            //260806Cl /simplify2 (F-1): Math.Max(0, NaN) は NaN を素通しし、ComputeProfile の maxIref 集約を
            //汚染して族全体を無言で消す。非有限値はここで 0 に正規化する
            I = double.IsFinite(intensity) ? Math.Max(0, intensity) : 0,
            IRef = double.IsFinite(iRef) ? Math.Max(0, iRef) : 0,
        };
    }

    #endregion

    #region systematic row (2N+1 波)

    /// <summary>
    /// 260807Cl 追加 (設計 Phase 2.5): systematic row {n·g : n = −N..+N} の x 非依存な前計算。
    /// U(m·g) と原子位置位相 e^{−2πi m g·r_a} はどちらも出射方向に依らないので、バンドごとに 1 回でよい
    /// (2 ビーム経路が点ごとに getU を呼んでいたのは snapshot 側の辞書キャッシュ頼み)。
    /// </summary>
    private sealed class RowPotential
    {
        public readonly int Order;      // N
        public readonly int BeamCount;  // 2N+1
        public readonly double GLength2; // |g|²
        private readonly Complex[] _u;       // _u[m + 2N] = U(m·g) (Real + i·Imag を畳んだ形)
        private readonly Complex[] _phase;   // _phase[(m + 2N) * _nSite + a] = e^{−2πi m g·r_a}
        private readonly int _nSite;

        public RowPotential(KikuchiPotentialSnapshot snap, in BandGeometry geo, int order)
        {
            Order = Math.Clamp(order, 1, Options.MaxRowOrder);
            BeamCount = 2 * Order + 1;
            GLength2 = geo.GLab.Length2;
            var (h, k, l) = geo.Family.Index;
            int span = 4 * Order + 1; // m = −2N..+2N
            _u = new Complex[span];
            _nSite = snap.Sites.Length;
            int nSite = _nSite;
            _phase = new Complex[span * nSite];
            for (int m = -2 * Order; m <= 2 * Order; m++)
            {
                int mi = m + 2 * Order;
                if (m == 0)
                    _u[mi] = Complex.Zero; // 対角は U'(0)+Q で別に入る (U(0) は k0 側に畳み込み済み)
                else
                {
                    var u = snap.GetU((m * h, m * k, m * l), m * geo.GLab);
                    _u[mi] = u.Real + Complex.ImaginaryOne * u.Imag;
                }
                for (int a = 0; a < nSite; a++)
                {
                    var s = snap.Sites[a];
                    //2 ビーム経路の q01 と同じ符号規約 (CreatePhaseFactors 準拠): e^{−2πi (h_j−h_i)·r_a}
                    var (sin, cos) = Math.SinCos(2 * Math.PI * m * (h * s.X + k * s.Y + l * s.Z));
                    _phase[mi * nSite + a] = new Complex(cos, -sin);
                }
            }
        }

        /// <summary>U(m·g)</summary>
        public Complex U(int m) => _u[m + 2 * Order];

        /// <summary>e^{−2πi m g·r_a}</summary>
        public Complex Phase(int m, int site) => _phase[(m + 2 * Order) * _nSite + site];
    }

    /// <summary>
    /// systematic row の 1 点計算 (2N+1 波)。規約はすべて 2 ビーム版 ComputePoint と同一:
    /// 行列 A[row + col·nb] = U(h_col − h_row)/P_col (非対角) / (i·U'(0) + Q_col)/P_col (対角)、
    /// 反転幾何 (n̂ = −d̂, k0 = −√(k_vac²+u0)·d̂)、源の運動量移行 q_i = q_0 − h_i、
    /// 厚み積分 F_{jj'}、I_ref = Bragg 結合を切った直進チャネル。
    /// ⚠ ±g が同じ基底に入るので**族としての ± 和は取らない** (直進チャネルも 1 本。ComputeProfile の refCount)。
    /// </summary>
    private static (double I, double IRef) ComputeRowPoint(KikuchiPotentialSnapshot snap, in BandGeometry geo,
        RowPotential row, double x, double thickness, Options opt)
    {
        var kVac = snap.KVac;
        var dHat = geo.Direction(x, opt.SampleAngle);
        var bHat = new Vector3DBase(0, 0, -1);
        int N = row.Order, nb = row.BeamCount;
        var gLab = geo.GLab;

        var nHat = -dHat;
        var k0 = -Math.Sqrt(kVac * kVac + snap.U0) * dHat;
        double k0Len2 = k0.Length2;

        // --- P_i = 2 n̂·(k0 + h_i), Q_i = |k0|² − |k0 + h_i|² ---
        var p = new double[nb];
        var bigQ = new double[nb];
        for (int i = 0; i < nb; i++)
        {
            var kh = k0 + (i - N) * gLab;
            p[i] = 2 * (nHat * kh);
            bigQ[i] = k0Len2 - kh.Length2;
        }
        //2 ビーム版と同じ扱い: ほぼ掠め角のビームは 2 波近似の外なので寄与 0 で返す
        if (p[N] <= 0)
            return (0, 0);
        for (int i = 0; i < nb; i++)
            if (p[i] <= 1e-3 * p[N])
                return (0, 0);

        // --- 固有値問題行列 (column-major。getEigenMatrix と同一の構成) ---
        var a = new Complex[nb * nb];
        var diagU = Complex.ImaginaryOne * snap.UPrime0;
        for (int col = 0; col < nb; col++)
        {
            var invP = 1.0 / p[col];
            for (int r = 0; r < nb; r++)
                a[r + col * nb] = r == col
                    ? (diagU + bigQ[col]) * invP
                    : opt.DisableBraggCoupling ? Complex.Zero : row.U(col - r) * invP;
        }

        // --- EVD と C⁻¹ (BetheMethod と同じ経路。native が無ければ MathNet) ---
        Complex[] eigVal, eigVec, eigInv;
        if (NativeWrapper.Enabled)
        {
            (eigVal, eigVec) = NativeWrapper.EigenSolver(nb, a);
            eigInv = NativeWrapper.Inverse(nb, eigVec);
        }
        else
        {
            var evd = new MathNet.Numerics.LinearAlgebra.Complex.DenseMatrix(nb, nb, a)
                .Evd(MathNet.Numerics.LinearAlgebra.Symmetricity.Asymmetric);
            eigVal = ((MathNet.Numerics.LinearAlgebra.Complex.DenseVector)evd.EigenValues).Values;
            eigVec = ((MathNet.Numerics.LinearAlgebra.Complex.DenseMatrix)evd.EigenVectors).Values;
            eigInv = ((MathNet.Numerics.LinearAlgebra.Complex.DenseMatrix)evd.EigenVectors.Inverse()).Values;
        }
        // α_j = [C⁻¹]_{j, N} (ψ0 = 中央ビーム = e_N)
        var alpha = new Complex[nb];
        for (int j = 0; j < nb; j++)
            alpha[j] = eigInv[j + N * nb];

        // --- 源密度行列 Q_ij = Σ_a Occ_a · coh_a(s_i², s_j², s_ij²) · e^{−2πi (h_j−h_i)·r_a} ---
        //     s_ij² = |q_i − q_j|²/4 = |h_j − h_i|²/4 = ((j−i)|g|)²/4
        var q0 = kVac * (dHat - bHat);
        var s2 = new double[nb];
        for (int i = 0; i < nb; i++)
            s2[i] = (q0 - (i - N) * gLab).Length2 / 4;
        if (opt.SwapSourceWeights)
            Array.Reverse(s2); // 診断: 源勾配の反転 (row では n → −n の鏡映。中央ビームは自分自身へ写る)

        int nSpec = snap.SpeciesCount;
        var coh = new double[nSpec];
        var qMat = new Complex[nb * nb];
        for (int i = 0; i < nb; i++)
            for (int j = 0; j < nb; j++)
            {
                if (opt.DiagonalSourceOnly && i != j)
                    continue; // 対角近似 (診断・単体テスト用)
                int m = j - i;
                var s2ij = m * m * row.GLength2 / 4;
                for (int sp = 0; sp < nSpec; sp++)
                    coh[sp] = snap.Kernel.SourceCoherence(sp, s2[i], s2[j], s2ij);
                Complex acc = Complex.Zero;
                for (int site = 0; site < snap.Sites.Length; site++)
                {
                    var st = snap.Sites[site];
                    acc += st.Occ * coh[st.AtomsIndex] * row.Phase(m, site);
                }
                qMat[i + j * nb] = acc;
            }

        // --- W_{jj'} = Σ_{ii'} C_{i,j} Q_{i,i'} C*_{i',j'} を 2 段の行列積で (O(nb³)) ---
        var t1 = new Complex[nb * nb]; // t1[i, j'] = Σ_{i'} Q_{i,i'} C*_{i',j'}
        for (int jp = 0; jp < nb; jp++)
            for (int i = 0; i < nb; i++)
            {
                Complex s = Complex.Zero;
                for (int ip = 0; ip < nb; ip++)
                    s += qMat[i + ip * nb] * Complex.Conjugate(eigVec[ip + jp * nb]);
                t1[i + jp * nb] = s;
            }

        // --- I(t) = Re Σ_{jj'} α_j α_j'* W_{jj'} F_{jj'}(t) ---
        Complex intensity = Complex.Zero;
        for (int j = 0; j < nb; j++)
            for (int jp = 0; jp < nb; jp++)
            {
                Complex w = Complex.Zero;
                for (int i = 0; i < nb; i++)
                    w += eigVec[i + j * nb] * t1[i + jp * nb];
                intensity += alpha[j] * Complex.Conjugate(alpha[jp]) * w * Fjj(eigVal[j], eigVal[jp], thickness);
            }

        // --- I_ref: Bragg 結合を切ると直進チャネルだけが残る (2 ビーム版と同式) ---
        var gRef = diagU / p[N];
        var iRef = (Fjj(gRef, gRef, thickness) * qMat[N + N * nb].Real).Real;

        //260806Cl /simplify2 (F-1) と同じ規律: 非有限値は描画層へ流さない
        return (double.IsFinite(intensity.Real) ? Math.Max(0, intensity.Real) : 0,
                double.IsFinite(iRef) ? Math.Max(0, iRef) : 0);
    }

    /// <summary>F_{jj'}(t) = [e^{λt}−1]/λ, λ = 2πi(γ_j − γ_j'*)。λ≈0 は t (EBSDSolverManaged と同一)</summary>
    private static Complex Fjj(in Complex gammaJ, in Complex gammaJp, double t)
    {
        var lam = 2 * Math.PI * Complex.ImaginaryOne * (gammaJ - Complex.Conjugate(gammaJp));
        return lam.MagnitudeSquared() < 1e-30 ? t : (Complex.Exp(lam * t) - Complex.One) / lam;
    }

    #endregion
}

//260806Cl /simplify 削除: KikuchiProfileCache (Phase 0 骨格)。参照ゼロのまま Form 側の freshness 判定と
//「stale の定義」が二重化し (cosε 量子化 vs 回転タプル)、キーが Options の診断フラグを含まないまま
//実装が乖離し始めていたため、実プロファイリングでキャッシュが要ると分かった時点で設計し直す。
//旧実装は git log (37b8a83 以前) 参照。
