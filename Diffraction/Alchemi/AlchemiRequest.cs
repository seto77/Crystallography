// 260807Cl 新規作成: ALCHEMI の run 単位の型 (A1′、設計 §5.2)。
// WinForms 非依存。GUI (FormALCHEMI, A4′) はこれらを組み立てて BetheMethod.RunAlchemi に渡すだけになる。
//
// 設計 §5.2 の AlchemiRequest / AlchemiResult / OrientationSample に対応する。v1 (1D forward) が
// 実際に使う項目だけを持たせ、2D・fit・自己吸収の項目は実装時に足す (空の器を先に置かない)。

using System;
using System.Collections.Generic;
using System.Linq;
using System.Security.Cryptography;
using System.Text;

namespace Crystallography;

/// <summary>260807Cl 追加: 走査する 1 方位 (設計 §5.2 の OrientationSample)。
/// ⚠方位は **BaseRotation ではなくビーム方向**で表す — Find_gVectors の BaseRotation を傾斜ごとに
/// 変えると gCache が毎回作り直しになって破滅的に遅くなる (設計 §8)。</summary>
/// <param name="Index">走査内の連番 (結果テンソルの第 1 添字)</param>
/// <param name="BeamDirection">入射方向の単位ベクトル (ReciPro 座標系。既定の正入射は (0,0,−1))</param>
/// <param name="TiltRad">走査中心からの符号付き傾斜角 [rad] (曲線の横軸。計算には使わない)</param>
/// <param name="Weight">ICP の規格化などで使う重み (設計 §3.6: per-maximum でなく参照方位集合の加重平均)</param>
public sealed record AlchemiOrientation(int Index, Vector3DBase BeamDirection, double TiltRad, double Weight = 1.0);

/// <summary>260807Cl 追加: 走査の組み立て。</summary>
public static class AlchemiScan
{
    /// <summary>1D ロッキング走査 (設計 §5.2 の Topology = OneDimensional)。
    /// <paramref name="tiltAxis"/> のまわりに <paramref name="center"/> を回した方向列を返す。</summary>
    /// <param name="center">走査中心のビーム方向 (null で (0,0,−1))</param>
    /// <param name="tiltAxis">傾斜軸 (center と直交していなくてよい — 直交成分だけを使う)</param>
    /// <param name="startRad">開始角 [rad]</param>
    /// <param name="endRad">終了角 [rad]</param>
    /// <param name="count">点数 (≥2)</param>
    public static AlchemiOrientation[] Rocking1D(Vector3DBase center, Vector3DBase tiltAxis, double startRad, double endRad, int count)
    {
        if (count < 2) throw new ArgumentOutOfRangeException(nameof(count), count, "1D rocking scan needs at least 2 points");
        var d0 = Vector3DBase.Normarize(center ?? new Vector3DBase(0, 0, -1));
        var axis = Vector3DBase.Normarize(tiltAxis ?? new Vector3DBase(1, 0, 0));
        //軸の d0 に平行な成分を落とす (残りが真の回転軸)
        axis = Vector3DBase.Normarize(axis - (axis * d0) * d0);
        if (double.IsNaN(axis.Length) || axis.Length2 < 1e-20)
            throw new ArgumentException("tiltAxis is parallel to the beam direction", nameof(tiltAxis));
        var perp = Vector3DBase.VectorProduct(axis, d0);//d0 ⊥ 面内のもう 1 本
        var result = new AlchemiOrientation[count];
        for (int i = 0; i < count; i++)
        {
            var theta = startRad + (endRad - startRad) * i / (count - 1);
            var (sin, cos) = Math.SinCos(theta);
            result[i] = new AlchemiOrientation(i, Vector3DBase.Normarize(cos * d0 + sin * perp), theta);
        }
        return result;
    }
}

/// <summary>260807Cl 追加: ALCHEMI run の要求 (設計 §5.2)。run 開始時に防御的コピーを取ってから使う。</summary>
/// <param name="IncidentEnergyKeV">加速電圧 [kV]</param>
/// <param name="BaseRotation">結晶方位 (走査中は固定。傾斜はビーム方向で表す)</param>
/// <param name="Orientations">走査する方位列</param>
/// <param name="ThicknessesNm">厚み [nm]</param>
/// <param name="Sites">サイト仮説の幾何部分</param>
/// <param name="Channels">イオン化チャネル (元素 × 殻)</param>
public sealed record AlchemiRequest(
    double IncidentEnergyKeV,
    Matrix3D BaseRotation,
    AlchemiOrientation[] Orientations,
    double[] ThicknessesNm,
    AlchemiSiteBasis[] Sites,
    IonizationChannelSpec[] Channels)
{
    /// <summary>試料表面 (から内部への) 法線。ReciPro の既定は (0,0,−1)。</summary>
    public Vector3DBase Surface { get; init; } = new Vector3DBase(0, 0, -1);

    /// <summary>物理モデル階層 (設計 §3.1)。既定 = v1 の公開モデル。</summary>
    public AlchemiModelTier ModelTier { get; init; } = AlchemiModelTier.LocalFormFactor;

    /// <summary>1 方位あたりの Bloch 波数の上限 (union はこれより増える)。</summary>
    public int MaxNumOfBloch { get; init; } = 300;

    /// <summary>dechannelling 項 (設計 §3.4)。v1 から必須なので既定 true。false は診断用。</summary>
    public bool IncludeDechannelledComponent { get; init; } = true;

    /// <summary>expanded-basis 診断の倍率 (設計 §5.4)。0 以下で診断を省略。</summary>
    public double ExpandedBasisFactor { get; init; } = 1.25;

    /// <summary>expanded-basis 診断の合否閾値 (Total の最大相対差)。設計 §5.4 の「最大相対差 ≤3e-3」。</summary>
    public double ExpandedBasisTolerance { get; init; } = 3e-3;

    /// <summary>方位ループの並列度 (-1 = 既定)。</summary>
    public int MaxDegreeOfParallelism { get; init; } = -1;

    /// <summary>native Eigen を使う (false で managed MathNet に固定。backend 一致検証用)。</summary>
    public bool UseNativeSolver { get; init; } = true;

    /// <summary>260813Cl 追加: **解決済みチャネルデータの差し替え (診断専用)**。null なら
    /// 通常どおり <see cref="IonizationDataProvider.Resolve(IonizationChannelSpec, double, IonizationFsTable)"/>
    /// が引く。<see cref="Channels"/> と同じ長さ・同じ順序でなければならない。
    ///
    /// ⚠ **用途は感度試験** — F(s) に既知の摂動 δF を載せた形状を渡し、同じ幾何・同じ
    /// Bloch 解の下で観測量 Y の変化を測る。F はμ を通してしか入らず、Bloch 係数は F に
    /// 依らないので、**2 回の run の差は δF の寄与だけを厳密に取り出す**。
    /// ⚠ 既定 null では 1 ビットも挙動が変わらない (`UseNativeSolver` と同じ位置づけ)。</summary>
    public IonizationData[] ChannelDataOverride { get; init; } = null;
}

/// <summary>260807Cl 追加: 結果テンソルの形 (flat 配列 [sample, thickness, site, channel] の添字計算)。</summary>
public sealed record AlchemiTensorShape(int OrientationCount, int ThicknessCount, int SiteCount, int ChannelCount)
{
    public int Length => OrientationCount * ThicknessCount * SiteCount * ChannelCount;

    /// <summary>flat 添字。channel が最内 = 同一方位・同一厚みのチャネル列が連続する。</summary>
    public int Index(int orientation, int thickness, int site, int channel)
        => ((orientation * ThicknessCount + thickness) * SiteCount + site) * ChannelCount + channel;
}

/// <summary>260809Cl 追加: 定量 fit に使ってよいかの**保証表示**。生の診断値 (bool) と分けるための三値。
/// v1 が常に <see cref="NotEvaluated"/> を返す理由は <see cref="AlchemiBasisDiagnostic.Eligibility"/>。</summary>
public enum AlchemiFitEligibility
{
    /// <summary>判定を出さない (v1)。「適格でない」ではなく「評価していない」。</summary>
    NotEvaluated,
    /// <summary>基底収束の診断が閾値内。</summary>
    Eligible,
    /// <summary>基底収束の診断が閾値を超えた。</summary>
    NotEligible,
}

/// <summary>260807Cl 追加: 基底 (FixedUnion) の診断 (設計 §5.4)。fit の可否判定に使うので結果に必ず保存する。</summary>
/// <param name="BeamCount">union の本数 (実際に解いた次元)</param>
/// <param name="CenterOnlyBeamCount">走査中心だけで Find_gVectors したときの本数</param>
/// <param name="MaxTiltRad">走査の全角幅 [rad]</param>
/// <param name="MaxMinExcitationErrorPerNm">union に入った g の「走査内 min|s|」の最大値 [nm⁻¹]。
/// 小さいほど「どの g も走査中どこかで励起される」= union が無駄に太っていない</param>
/// <param name="MaxShapeArgumentAngstromInv">この基底が F(s) に要求する s = |G|/2 の最大値 [Å⁻¹]
/// (G = g_h − g_g なので max|g| に等しい)。260807Cl 追加: 小さい単位胞を系統反射列条件で解くと
/// Find_gVectors が**長い 1 次元の列**を選ぶため、ビーム数が穏やかでも |g| が線形に伸びて
/// F(s) テーブルの収録上限 (s ≤ 8 Å⁻¹) を超え得る。**run 前に警告できる唯一の量**</param>
/// <param name="BasisHash">sorted hkl の SHA-256 先頭 16 桁 (基底が同一かの照合用)</param>
/// <param name="ExpandedBasisMaxRelDiff">1.25 倍基底との Total 最大相対差 (NaN = 診断未実施)</param>
/// <param name="AcceptedForFit">expanded-basis 診断が閾値内に収まったか (**生の診断結果**。
/// 公開表示に使ってよいかは <see cref="Eligibility"/> を見ること)</param>
/// <param name="Warnings">傾斜幅・基底膨張・F(s) 収録範囲などの警告</param>
public sealed record AlchemiBasisDiagnostic(
    int BeamCount, int CenterOnlyBeamCount, double MaxTiltRad,
    double MaxMinExcitationErrorPerNm, double MaxShapeArgumentAngstromInv, string BasisHash,
    double ExpandedBasisMaxRelDiff, bool AcceptedForFit, string[] Warnings)
{
    /// <summary>union が中心のみの基底より何本増えたか。</summary>
    public int AddedByUnion => BeamCount - CenterOnlyBeamCount;

    /// <summary>260809Cl 追加: 「fit 適格」と**保証表示してよいか**。GUI と CSV はこちらを使う。
    /// <para>
    /// ⚠ 現在は常に <see cref="AlchemiFitEligibility.NotEvaluated"/> を返す (作者決定 260809Cl)。
    /// 指示書 §2-8 の 3 点が未修正なため:
    /// ⑤ 分母が (方位×厚み×サイト×チャネル) テンソル全体の最大値なので σ の小さいチャネルが薄まる /
    /// ⑥ 分子が絶対収率なので ICP では落ちる共通スケール変化も数える /
    /// ⑦ Find_gVectors がエワルド球スクリーニング律速で「1.25 倍にしても基底が増えない」領域があり、
    ///    そこでは**基底を自分自身と比べて**しまうので <see cref="AcceptedForFit"/> = true が
    ///    「収束した」を意味しない (= 偽陽性)。
    /// </para>
    /// <para>
    /// 偽陽性がある状態で「適格」と保証表示するのは誤りの方向が悪いので、v1 は判定を出さない。
    /// **⑦ → ⑤ → ⑥ の順で直したら、ここを <c>AcceptedForFit ? Eligible : NotEligible</c> に戻すこと。**
    /// <see cref="AcceptedForFit"/> と <see cref="ExpandedBasisMaxRelDiff"/> は生の診断値として残してあるので、
    /// 回帰テスト (AlchemiCheck) はそのまま使える。
    /// </para></summary>
    public AlchemiFitEligibility Eligibility => AlchemiFitEligibility.NotEvaluated;

    /// <summary>sorted hkl から basis hash を作る。</summary>
    internal static string Hash(IEnumerable<(int H, int K, int L)> indices)
    {
        var sb = new StringBuilder();
        foreach (var (h, k, l) in indices.OrderBy(v => v.H).ThenBy(v => v.K).ThenBy(v => v.L))
            sb.Append(h).Append(',').Append(k).Append(',').Append(l).Append(';');
        return Convert.ToHexString(SHA256.HashData(Encoding.ASCII.GetBytes(sb.ToString())))[..16];
    }
}

/// <summary>260807Cl 追加: 進捗ステージ。表示文言は GUI 側で組む (backend は表示都合の文字列を持たない)。
/// 設計 §5.5 の 7 段階のうち、角度広がり畳み込みと曲線合成は A4′ (GUI) 側の後処理なのでここには無い。</summary>
public enum AlchemiStage
{
    ResolvingIonizationData = 1,
    BuildingUnionBasis = 2,
    BuildingMuMatrices = 3,
    SolvingOrientations = 4,
    ExpandedBasisCheck = 5,
}

/// <summary>260807Cl 追加: 進捗通知。</summary>
/// <param name="Stage">ステージ</param>
/// <param name="Fraction">当該ステージ内の進捗 0–1</param>
public sealed record AlchemiProgress(AlchemiStage Stage, double Fraction);

/// <summary>260807Cl 追加: ALCHEMI run の結果 (設計 §5.2)。構築後は不変。
/// yield の単位 = 入射電子 1 個あたりの発生イオン化数 (無次元)、Tracer 基底 (§3.2)。</summary>
public sealed class AlchemiResult
{
    public AlchemiTensorShape Shape { get; init; }
    /// <summary>動力学項 flat [sample, thickness, site, channel]</summary>
    public double[] Dynamic { get; init; }
    /// <summary>非チャネリング項 (同上)</summary>
    public double[] Dechannelled { get; init; }
    /// <summary>Dynamic + Dechannelled (同上)</summary>
    public double[] Total { get; init; }
    /// <summary>L_coh(t) flat [sample, thickness] — 診断・回帰用</summary>
    public double[] CoherentPathLengthNm { get; init; }
    public AlchemiBasisDiagnostic Basis { get; init; }
    /// <summary>解決済みチャネルデータ (σ・F(s) の provenance を含む)</summary>
    public IonizationData[] ChannelData { get; init; }
    public AlchemiOrientation[] Orientations { get; init; }
    public double[] ThicknessesNm { get; init; }
    public AlchemiSiteBasis[] Sites { get; init; }
    public AlchemiModelTier ModelTier { get; init; }
    public double IncidentEnergyKeV { get; init; }
    /// <summary>μ の対角 [nm²] [site, channel] (dechannelling と規格化の参照値)</summary>
    public double[] Mu00Nm2 { get; init; }
    /// <summary>単位胞体積 [nm³]</summary>
    public double UnitCellVolumeNm3 { get; init; }
    /// <summary>物理量と規格化 (設計 §3.6: 表示正規化とは別管理)</summary>
    public SignalQuantity Quantity { get; init; } = SignalQuantity.IonizationVacanciesGenerated;
    public SignalNormalization Normalization { get; init; } = SignalNormalization.PerIncidentElectron;
    /// <summary>サイト応答の線形合成が許されるか (Tracer 近似のみ true。設計 §3.5)</summary>
    public bool LinearCombinationValid { get; init; } = true;

    /// <summary>1 (サイト, チャネル) の曲線を厚み固定で切り出す (GUI・CSV export の入口)。</summary>
    public double[] Curve(double[] tensor, int thicknessIndex, int siteIndex, int channelIndex)
    {
        var curve = new double[Shape.OrientationCount];
        for (int o = 0; o < curve.Length; o++)
            curve[o] = tensor[Shape.Index(o, thicknessIndex, siteIndex, channelIndex)];
        return curve;
    }
}
