using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using V3 = OpenTK.Mathematics.Vector3d;
using System.Threading;

namespace Crystallography;

//Electron beam-specimen interactions and simulation methods in microscopy 2018
//Eqs (2.38), (2.41), (2.42) などを参考
public class MonteCarlo
{
    #region モデル
    public enum StoppingPowerModels
    {
        /// <summary>
        /// Joy &amp; Luo (1989) の経験的阻止能モデル。
        /// Bethe 式を低エネルギー側に拡張したもので、物質定数 k と平均イオン化ポテンシャル J を用いる。
        /// 計算が軽く実装も単純だが、低エネルギー領域 (&lt;5 keV) での精度はやや劣る。
        /// </summary>
        JoyLuo1989, // 260331Cl コメント追加

        /// <summary>
        /// Jablonski (2008) の修正阻止能モデル。
        /// TPP-2M で求めた非弾性平均自由行程 (IMFP) を基に阻止能を導出する。
        /// 密度 ρ やプラズマエネルギーなどの材料パラメータを反映するため、
        /// 特に低エネルギー領域や EBSD で扱う BSE エネルギー範囲での精度が JoyLuo1989 より高い。
        /// </summary>
        JablonskiModified2008, // 260331Cl コメント追加
    }

    public enum ElasticScatteringModels
    {
        /// <summary>
        /// 遮蔽 Rutherford 散乱モデル。
        /// 散乱角分布を解析的に計算でき非常に高速。軽元素〜中程度の原子番号では十分な精度がある。
        /// 重元素 (Au, W など) では大角散乱の確率を過小評価する傾向があり、
        /// 定量的な BSE 深さ分布の精度が必要な場合は MottNistSampler2023 が望ましい。
        /// </summary>
        ScreenedRutherford, // 260331Cl コメント追加

        /// <summary>
        /// NIST Electron Elastic-Scattering Cross-Section Database (SRD 64) に基づく Mott 散乱サンプラー。
        /// 部分波展開による正確な微分断面積テーブルから散乱角を逆関数法でサンプリングする。
        /// 重元素で Screened Rutherford との差が大きく、BSE 収率や深さ分布の精度が向上する。
        /// データファイル (E_ZZ.TXT) の読み込みが必要で、初回にやや時間がかかる。
        /// 混合物系では元素ごとの巨視的断面積で重み付けして散乱元素を選択する。
        /// </summary>
        MottNistSampler2023, // 260331Cl コメント追加
    }

    #region お蔵入り // (260401Ch) generated / external の source flag 比較経路は配布版では使わない
    // public enum ElasticSamplerDataSources
    // {
    //     /// <summary>260401Ch generated PCHIP を優先し、存在しなければ外部 TXT にフォールバックする既定動作</summary>
    //     Auto,
    //     /// <summary>260401Ch generated PCHIP のみを使う。未生成元素は Mott sampler が使えず Screened Rutherford に落ちる</summary>
    //     GeneratedOnly,
    //     /// <summary>260401Ch 外部 TXT テーブルのみを使う。圧縮前との比較ベンチ用</summary>
    //     ExternalTextOnly,
    // }
    #endregion

    public enum InelasticScatteringModels
    {
        /// <summary>
        /// 連続減速近似 (CSDA)。非弾性散乱を離散イベントとして扱わず、
        /// ステップごとに阻止能 × 飛行距離だけエネルギーを連続的に減少させる。
        /// 最も高速だが「最後の非弾性散乱」を定義できないため、
        /// EBSD の非弾性散乱深さ・エネルギー分布の解析には使用不可 (HasLastInelasticEvent = false)。
        /// </summary>
        ContinuousSlowingDownApproximation, // 260331Cl コメント追加

        /// <summary>
        /// 離散非弾性散乱モデル (平均損失)。非弾性散乱イベントごとに平均エネルギー損失を一定値として失う。
        /// イベント発生位置 (深さ) の統計は正しく再現されるが、エネルギー損失の確率的ばらつきがない。
        /// 計算速度は離散モデルの中で最速。深さ分布の概形を素早く確認するのに適する。
        /// </summary>
        DiscreteMeanLoss, // 260331Cl コメント追加

        /// <summary>
        /// 離散非弾性散乱モデル (簡易バルク DIIMFP 近似)。
        /// プラズモンピーク・低損失テール・高損失テールの 3 成分を混合した損失分布からサンプリングする。
        /// エネルギー損失の確率的ばらつきを物理的に最も妥当な形で再現するため、
        /// 「最後の非弾性散乱後のエネルギー」の分布解析に最適。計算コストは DiscreteMeanLoss よりやや大きい。
        /// </summary>
        DiscreteBulkDiimfpApproximation, // 260331Cl コメント追加

        /// <summary>
        /// 260922Cl 追加: 離散非弾性散乱モデル (価電子の拡張 Drude + 内殻)。事象の発生率 1/λ_in (TPP-2M) と平均損失 S·λ_in は
        /// <see cref="DiscreteBulkDiimfpApproximation"/> と同じで、<b>1 事象の損失分布の形</b>だけを物理的にしたもの。
        /// <para>価電子: 単極の拡張 Drude 損失関数 Im[−1/ε(q,ω)] = ω_p²γω/((ω²−ω_q²)²+γ²ω²)、ω_q = ω_p + q²/2 (原子単位、ω_p = TPP の E_p、γ = 4 eV) を
        /// q で積分した DIIMFP (Bethe ridge = 価電子の二体衝突の尾まで含む)。Si 20 keV で λ_vb 30.5 nm、平均 30 eV、プラズモン付近 (E_p ± 5 eV) が 67 %。</para>
        /// <para>内殻: 吸収端が <see cref="InelasticLocalizedLossEv"/> 以上の副殻 (Bote–Salvat)。副殻は Bote–Salvat の電離断面積の比で選び、
        /// 損失は吸収端 B〜E/2 で ∝ ω⁻²。内殻を選ぶ割合は 1 事象の平均損失が S·λ_in になるように決める (Si 20 keV で 11.7 %、
        /// Vos &amp; Winkelmann 2019 表 1 の 12.0 % と独立に一致)。</para>
        /// <para>旧モデルとの違い (Si 20 keV): 旧はプラズモン成分の重みが 0 に飽和し 30 eV 以上が 85 %、平均が目標の 0.89、負の損失 0.65 %。
        /// 検討の経緯は .project-guidance/ReciPro/ReciPro_EBSD_密度行列定式化_案d.md §2.4。</para>
        /// </summary>
        DiscreteDrudeValenceInnerShell,
    }
    #endregion

    #region const, static 定数

    public const double Th = 0.0000001; // 260401Cl 方向ベクトル回転時の特異点判定閾値。vZ ≈ -1 (ほぼ真下向き) のとき回転行列が退化するのを避ける

    private static readonly double log50 = Math.Log(50.0); // (260331Ch) NIST テーブルの最小エネルギー 50 eV の自然対数
    private static readonly double log20000 = Math.Log(20000.0); // (260331Ch) sampler 基準点 20 keV の自然対数。260603Cl: テーブル上限ではなく対数刻みの基準 (50eV-20keV を100分割)
    private static readonly double LogNistElasticEnergyStep = (log20000 - log50) / 100.0; // (260331Ch) 対数エネルギー軸の刻み。260603Cl: 20keV超も同刻みで等間隔延長 (blockIndex 0-110)


    // 260401Cl Jablonski (2008) 修正阻止能モデルのフィッティング定数 D1〜D5。
    // 阻止能 S(E) = D1·E^D2·ln(D3·E)·(ρ+D4)^D5 / λ_in  [eV/Å] の形で使用。
    private const double Jablonski2008D1 = 7.89271; // (260331Ch)
    private const double Jablonski2008D2 = 0.0117088; // (260331Ch)
    private const double Jablonski2008D3 = 0.0545126; // (260331Ch)
    private const double Jablonski2008D4 = -0.0254488; // (260331Ch)
    private const double Jablonski2008D5 = 0.326907; // (260331Ch)
    // private const int NistElasticEnergyCount = 101; // 260603Cl 変更前 (50eV-20keV, 101点)
    private const int NistElasticEnergyCount = 111; // 260603Cl 50eV-36.4keV へ拡張 (20keV超を NIST DCS から継ぎ足し)。log spaced, step は不変
    private static readonly double NistElasticMaxEnergyEv = Math.Exp(log50 + LogNistElasticEnergyStep * (NistElasticEnergyCount - 1)); // 260603Cl 追加: テーブル上限 = blockIndex 110 = 36411 eV
    #region お蔵入り // (260401Ch) オリジナル TXT の 2001 点 CDF は配布版ランタイムでは使わない
    // private const int NistElasticPhiCount = 2001; // (260331Ch) X = cos(theta) from +1 to -1 with step 0.001
    #endregion

    #endregion



    internal readonly Random Rnd = Random.Shared; //260922Cl private → internal (tools/EbsdSourceStudy の源定義の研究用。本番の挙動は不変)
    /// <summary>平均原子番号 (混合物の場合は重み付き平均) </summary>
    public readonly double Z;
    /// <summary>平均原子量 (g/mol) </summary>
    public readonly double A;
    /// <summary>密度 (g/cm³) </summary>
    public readonly double ρ;  // 260401Cl 密度 (g/cm³)
    /// <summary>260401Cl 入射電子のエネルギー (keV)</summary>
    public readonly double InitialKev;
    /// <summary>260401Cl 試料表面の傾斜角 (rad, X軸回り)</summary>
    public readonly double Tilt;
    /// <summary>260401Cl Screened Rutherford 遮蔽パラメータ α の係数: 0.0034 * Z^(2/3)。α = coeff0 / E で遮蔽効果を表す</summary>
    public readonly double coeff0;
    /// <summary>260401Cl 弾性散乱断面積 σ_E の係数: Ze²/(8πε₀)。Rutherford 散乱の微分断面積の前因子</summary>
    public readonly double coeff1;
    /// <summary>260401Cl 弾性平均自由行程 λ_el の係数: A/(N_A·ρ) [nm³]。λ_el = coeff2 / σ_E で求まる</summary>
    public readonly double coeff2;
    /// <summary>260401Cl Joy-Luo 阻止能 dE/ds の係数 (keV/nm 単位)。Bethe 式の前因子 -Z·N_A·ρ·e⁴/(4πε₀²·A) を含む</summary>
    public readonly double coeff3;
    /// <summary>260401Cl Joy-Luo 阻止能の低エネルギー補正係数。Z 依存の経験的パラメータ (k = 0.0299·ln(Z) + 0.7307)</summary>
    public readonly double k;
    /// <summary>260401Cl 平均イオン化ポテンシャル (eV)。Bethe 阻止能式で物質のエネルギー損失特性を決める定数</summary>
    public readonly double J;
    /// <summary>260401Cl tan(tilt): 試料表面境界の判定に使用 (Y·tan ≥ Z で試料内)</summary>
    public readonly double tan;
    /// <summary>260401Cl cos(tilt): 深さ d の更新時の射影成分 (d += s·(sin·vY - cos·vZ))</summary>
    public readonly double cos;
    /// <summary>260401Cl sin(tilt): 深さ d の更新時の射影成分</summary>
    public readonly double sin;
    /// <summary>(260331Ch) TPP-2M に使う平均価電子数 Nv</summary>
    public readonly double ValenceElectronCount;
    /// <summary>260921Cl 追加: 調査用の「干渉性区間」の累積 (<see cref="BackscatteredElectronDetail.CoherentInelasticCount"/> ほか) を集めるか。
    /// ⚠ 既定 false。true にすると非弾性散乱 1 回ごとに <see cref="CoherencePreservingAngularVariance"/> (Math.Log を含む) を
    /// MC の最内ループで呼ぶ。260920Cl に無条件で入れてしまい、誰も読まない量のために電子 250 万本ぶんの計算を払っていた。
    /// 調査を再開するときだけ true にすること。
    /// ⚠ false でも <see cref="BackscatteredElectronDetail.CoherentPathLengthNm"/> だけは集め続ける
    /// (加算 1 回なので止める価値が無い)。0 のままになるのはイベント数と角度分散の 2 つ。</summary>
    public static bool CollectCoherenceDiagnostics = false;
    /// <summary>(260331Ch) バンドギャップ Eg (eV)</summary>
    public readonly double BandGapEv;
    /// <summary>(260331Ch) 阻止能モデルの切替</summary>
    public readonly StoppingPowerModels StoppingPowerModel;
    /// <summary>(260331Ch) 弾性散乱モデルの切替</summary>
    public readonly ElasticScatteringModels ElasticScatteringModel;
    /// <summary>(260331Ch) 非弾性散乱モデルの切替</summary>
    public readonly InelasticScatteringModels InelasticScatteringModel;
    /// <summary>260401Cl シミュレーション打ち切りエネルギー (keV)。電子エネルギーがこれ以下になると追跡を終了する</summary>
    public readonly double ThresholdKev;
    /// <summary>260919Cl 追加: 弾性散乱のコヒーレンス破壊確率 1−exp(−2B s²) に使う組成平均の等方 B [nm²] (下限 1e-3 nm²)。
    /// Bloch 波が扱う干渉性 Bragg 散乱の割合 exp(−2B s²) は「イベント無し」と等価に扱い、残りの熱散漫成分だけが源の深さをリセットする。</summary>
    public readonly double MeanDebyeWallerBNm2;
    /// <summary>260919Cl 追加: 非弾性散乱が局在 (コヒーレンス破壊) と判定される運動量移行のしきい値 q_c = ω_p / v_F (プラズモンカットオフ) を k=1/λ 系 [nm⁻¹] (物理慣例の値 ÷ 2π) で保持。
    /// q &lt; q_c は集団励起 (プラズモン) / 非局在の電子正孔対で Bloch 状態を保つ。q &gt; q_c は単一粒子励起 (Bethe ridge) で反跳電子が事象を 1/q に局在させる。</summary>
    public readonly double InelasticLocalizationQNm;
    /// <summary>260919Cl 追加: これ以上の損失 [eV] は内殻・単一粒子励起として運動量移行によらず局在とみなす (損失スペクトルの高損失テール開始点 max(1.8 E_p, 30 eV) と同じ)。</summary>
    public readonly double InelasticLocalizedLossEv;
    /// <summary>(260331Ch) Mott/NIST sampler 用の元素組成と数密度</summary>
    private readonly ElasticSpecies[] ElasticComponents = [];
    /// <summary>(260331Ch) 混合系の全元素数密度合計 [1/nm³]。巨視的断面積 Σ = Σ_i(n_i·σ_i) から有効微視的断面積 σ_eff = Σ/n_total を逆算する際に使う</summary>
    private readonly double TotalElasticNumberDensityPerNm3;
    /// <summary>260401Cl TPP-2M (Tanuma-Powell-Penn, 2M パラメータ版) の材料定数。非弾性平均自由行程 λ_in(E) = E / {Ep²·[β·ln(γE) - C/E + D/E²]} で使用</summary>
    /// <summary>260401Cl 自由電子プラズマエネルギー Ep = 28.8·√(Nv·ρ/A) [eV]。価電子のプラズモン励起エネルギーに対応</summary>
    private readonly double TppPlasmaEnergyEv;
    /// <summary>260606Cl 追加: 自由電子プラズマエネルギー Ep [eV] の読み取り公開 (Beam Interaction の電子輸送スカラ表示用)。</summary>
    public double PlasmaEnergyEv => TppPlasmaEnergyEv;
    /// <summary>260401Cl TPP-2M の β パラメータ。ln(γE) の係数で、IMFP のエネルギー依存の主項を制御</summary>
    private readonly double TppBeta;
    /// <summary>260401Cl TPP-2M の γ パラメータ [1/eV]。γE が対数項の引数で、密度依存 (0.191·ρ^(-0.5))</summary>
    private readonly double TppGamma;
    /// <summary>260401Cl TPP-2M の C パラメータ。1/E 項の係数で低エネルギー側の IMFP 補正に寄与</summary>
    private readonly double TppC;
    /// <summary>260401Cl TPP-2M の D パラメータ。1/E² 項の係数でさらに低エネルギー側の補正</summary>
    private readonly double TppD;
    /// <summary>(260331Ch) 1 eV 刻みの輸送パラメータ cache</summary>
    private readonly TransportParameters[] TransportParameterCache = [];
    /// <summary>(260331Ch) NIST sampler の混合断面積/CDF cache</summary>
    private readonly MottElasticMixtureEntry[] MottElasticMixtureCache = [];
    /// <summary>(260331Ch) 軽量な bulk DIIMFP 近似 sampler cache</summary>
    private readonly BulkLossSamplerEntry[] BulkLossSamplerCache = [];
    /// <summary>260922Cl 追加: <see cref="InelasticScatteringModels.DiscreteDrudeValenceInnerShell"/> の損失サンプラー (対数等間隔のエネルギー点ごと)</summary>
    private readonly DrudeLossSamplerEntry[] DrudeLossSamplerCache = [];
    /// <summary>260922Cl 追加: <see cref="DrudeLossSamplerCache"/> のエネルギー点 [keV] (昇順、対数等間隔)</summary>
    private readonly double[] DrudeLossSamplerEnergiesKev = [];
    #region お蔵入り // (260401Ch) generated / external の source flag 比較経路は配布版では使わない
    // /// <summary>260401Ch generated / external / auto の sampler source flag。MC ベンチで圧縮版と元テーブルを比較するために使う</summary>
    // public readonly ElasticSamplerDataSources ElasticSamplerDataSource;
    //
    // /// <summary>(260331Ch) 配布物に同梱した sampler を既定で使う</summary>
    // public static string NistElasticSamplerDirectory { get; set; } = Path.Combine(AppContext.BaseDirectory, "NistElasticSampler");
    // /// <summary>260401Ch 元の TXT テーブル置き場。ベンチ時は source tree の NistElasticSampler_Original を指す</summary>
    // public static string NistElasticSamplerTextDirectory { get; set; } = Path.Combine(AppContext.BaseDirectory, "NistElasticSampler_Original");
    // /// <summary>260401Ch constructor で source 未指定時に使う既定値</summary>
    // public static ElasticSamplerDataSources DefaultElasticSamplerDataSource { get; set; } = ElasticSamplerDataSources.Auto;
    #endregion

    /// <summary>(260331Ch) NIST テーブル cache のスレッド同期用</summary>
    private static readonly Lock lockObj = new();
    /// <summary>260401Cl 原子番号→NIST弾性散乱テーブルの static cache。全 MonteCarlo インスタンスで共有</summary>
    #region お蔵入り // (260401Ch) source flag ごとの cache は配布版では使わない
    // private static readonly Dictionary<(int AtomicNumber, ElasticSamplerDataSources Source), NistElasticScatteringTable> NistElasticSamplerCache = [];
    #endregion
    private static readonly Dictionary<int, NistElasticScatteringTable> NistElasticSamplerCache = []; // (260401Ch) 配布版は generated data のみを原子番号ごとに cache する
    /// <summary>(260331Ch) Bohr 半径の二乗 a₀² から nm² への換算係数: 1 a₀² = 2.8002852×10⁻³ nm²</summary>
    private const double NistElasticCrossSectionUnitNm2 = 2.8002852E-3;

    #region class, record

    /// <summary>混合物中の各元素の原子番号・数密度・NIST散乱テーブルを保持する。Mott 散乱サンプリング時に元素選択と角度サンプリングに使用。</summary>
    private readonly record struct ElasticSpecies(int AtomicNumber, double NumberDensityPerNm3, NistElasticScatteringTable NistElasticTable); // (260331Ch) Mott/NIST sampler の毎イベント lock/dictionary lookup を避ける

    /// <summary>NIST SRD 64 の弾性散乱データを保持するテーブル。111 エネルギー点 (50 eV〜36.4 keV, 対数等間隔。50eV-20keVは sampler 由来、20keV超は DCS 由来。260603Cl 拡張) ごとに全断面積と累積角度分布を格納。</summary>
    private sealed class NistElasticScatteringTable // (260331Ch)
    {
        public int AtomicNumber;
        /// <summary>260401Cl PCHIP 補間で生成済みの高速ルックアップデータ (存在すれば Phi[] より優先)</summary>
        public NistElasticPchipRuntimeElement GeneratedPchipRuntimeElement;
        /// <summary>260401Cl 各エネルギーでの弾性散乱全断面積 σ [a₀² 単位]</summary>
        public readonly double[] SigmaA0Squared = new double[NistElasticEnergyCount];
        #region お蔵入り // (260401Ch) オリジナル TXT の 2001 点 CDF は配布版ランタイムでは使わない
        // /// <summary>260401Cl 累積散乱角分布 Φ(cosθ)。cosθ = 1〜-1 を 0.001 刻み 2001 点で格納。逆関数法で散乱角をサンプリング</summary>
        // public readonly double[][] Phi = new double[NistElasticEnergyCount][];
        #endregion
    }

    /// <summary>混合物系の弾性散乱における各エネルギーでの巨視的断面積と元素選択 CDF。散乱イベント発生時にどの元素で散乱するかを確率的に決定する。</summary>
    private sealed class MottElasticMixtureEntry(double totalMacroscopicCrossSectionPerNm, double[] cumulativeProbabilities) // (260331Ch)
    {
        /// <summary>260401Cl 全元素の巨視的弾性散乱断面積の合計 Σ_total = Σ_i(n_i·σ_i) [1/nm]。逆数が弾性平均自由行程</summary>
        public readonly double TotalMacroscopicCrossSectionPerNm = totalMacroscopicCrossSectionPerNm;
        /// <summary>260401Cl 元素選択用の累積確率。i 番目の元素が選ばれる確率 = n_i·σ_i / Σ_total</summary>
        public readonly double[] CumulativeProbabilities = cumulativeProbabilities;
    }

    /// <summary>DiscreteBulkDiimfpApproximation 用のエネルギー損失サンプラー。低損失・プラズモン・高損失テールの 3 成分混合 PDF から CDF を構築し、逆関数法でサンプリングする。</summary>
    // private sealed class BulkLossSamplerEntry(double minLossKev, double lossStepKev, double[] cumulativeProbabilities) // 260401Cl 旧シグネチャ
    private sealed class BulkLossSamplerEntry(double minLossKev, double lossStepKev, double[] cumulativeProbabilities, byte[] guideTable) // 260401Cl GuideTable 追加
    {
        /// <summary>260401Cl 損失スペクトルの最小エネルギー損失 (keV)。バンドギャップ以上の値</summary>
        public readonly double MinLossKev = minLossKev;
        /// <summary>260401Cl CDF の 1 ビンあたりのエネルギー幅 (keV)</summary>
        public readonly double LossStepKev = lossStepKev;
        /// <summary>260401Cl エネルギー損失分布の累積確率 (256 ビン)。逆関数法でサンプリング</summary>
        public readonly double[] CumulativeProbabilities = cumulativeProbabilities;
        /// <summary>260401Cl CDF バイナリサーチを O(1) に高速化するガイドテーブル (64 エントリ)</summary>
        public readonly byte[] GuideTable = guideTable;
    }

    /// <summary>後方散乱電子の詳細情報。EBSD パターン形成に寄与する電子の最後の非弾性散乱の深さ・エネルギー・方向を記録する。</summary>
    /// <param name="Depth">260401Cl 電子が試料表面を脱出した (または停止した) 時点での表面からの深さ (nm)</param>
    /// <param name="Direction">260401Cl 脱出時の進行方向の単位ベクトル</param>
    /// <param name="Energy">260401Cl 脱出時のエネルギー (keV)</param>
    /// <param name="TotalEnergyLoss">260401Cl 入射エネルギーからの総エネルギー損失 (keV)</param>
    /// <param name="HasLastInelasticEvent">260401Cl 離散非弾性散乱イベントが 1 回以上発生したか (CSDA モードでは常に false)</param>
    /// <param name="LastInelasticDepth">260401Cl 最後の非弾性散乱が起きた深さ (nm)。EBSD の情報深さに直結</param>
    /// <param name="LastInelasticEnergyBeforeLoss">260401Cl 最後の非弾性散乱直前のエネルギー (keV)</param>
    /// <param name="LastInelasticEnergyAfterLoss">260401Cl 最後の非弾性散乱直後のエネルギー (keV)。この電子が回折に寄与する</param>
    /// <param name="LastInelasticDirection">260401Cl 最後の非弾性散乱時点での進行方向。回折条件の評価に使用</param>
    /// <param name="HasLastDecoherenceEvent">260919Cl 追加: コヒーレンス破壊イベント (非弾性、または DW 確率で熱散漫と判定された弾性) が 1 回以上発生したか</param>
    /// <param name="LastDecoherenceDepth">260919Cl 追加: 最後のコヒーレンス破壊イベントの深さ (nm)。ここから表面までが Bloch 波で扱う干渉性の経路</param>
    // public readonly record struct BackscatteredElectronDetail( // 260919Cl 変更前 (HasLastDecoherenceEvent / LastDecoherenceDepth 無し)
    //     double Depth, V3 Direction, double Energy, double TotalEnergyLoss, bool HasLastInelasticEvent,
    //     double LastInelasticDepth, double LastInelasticEnergyBeforeLoss, double LastInelasticEnergyAfterLoss, V3 LastInelasticDirection);
    /// <param name="CoherentInelasticCount">260920Cl 追加 (調査用): 最後のコヒーレンス破壊イベント以降に起きた、コヒーレンスを壊さない非弾性散乱 (q &lt; q_c) の回数</param>
    /// <param name="CoherentAngularVarianceRad2">260920Cl 追加 (調査用): 同じ区間で累積した角度偏位の二乗平均 Σ⟨θ²⟩ [rad²]。
    ///   菊池パターンは「最後の破壊以降は完全に干渉性」と扱われているが、q &lt; q_c の非弾性散乱は角度を僅かに振る。
    ///   その累積が帯の角幅 (Si {220} で 2.56°) に対してどれだけかを測るための量。弾性散乱の角度は含めない
    ///   (干渉性の弾性散乱は Bloch 波計算そのものなので、ぼかしとして数えると二重計上になる)</param>
    /// <param name="CoherentPathLengthNm">260920Cl 追加 (調査用): 最後のコヒーレンス破壊イベント以降の実経路長 [nm]</param>
    //旧: …, bool HasLastDecoherenceEvent, double LastDecoherenceDepth); // 260920Cl 変更前
    public readonly record struct BackscatteredElectronDetail( // (260331Ch) EBSD 寄与電子の最後の非弾性散乱情報を後段で解析できるようにする
        double Depth, V3 Direction, double Energy, double TotalEnergyLoss, bool HasLastInelasticEvent,
        double LastInelasticDepth, double LastInelasticEnergyBeforeLoss, double LastInelasticEnergyAfterLoss, V3 LastInelasticDirection,
        bool HasLastDecoherenceEvent, double LastDecoherenceDepth, // 260919Cl 追加
        int CoherentInelasticCount, double CoherentAngularVarianceRad2, double CoherentPathLengthNm); // 260920Cl 追加 (調査用)

    /// <summary>あるエネルギーにおける電子輸送パラメータの一式。弾性・非弾性散乱のステップ長と方向・エネルギー損失の計算に使う。</summary>
    /// <param name="ScreeningParameter">260401Cl Screened Rutherford の遮蔽パラメータ α = coeff0/E。原子核電荷の遮蔽効果を表し、散乱角分布の前方集中度を制御</param>
    /// <param name="ElasticCrossSectionNm2">260401Cl 弾性散乱全断面積 σ_el [nm²]。単一散乱イベントの確率を決める</param>
    /// <param name="ElasticMeanFreePathNm">260401Cl 弾性平均自由行程 λ_el = 1/(n·σ_el) [nm]。連続する弾性散乱間の平均飛行距離</param>
    /// <param name="StoppingPowerKevPerNm">260401Cl 阻止能 dE/ds [keV/nm] (負値)。単位飛行距離あたりのエネルギー損失率</param>
    /// <param name="InelasticMeanFreePathNm">260401Cl 非弾性平均自由行程 λ_in [nm]。TPP-2M で計算。連続する非弾性散乱間の平均飛行距離</param>
    /// <param name="MeanInelasticLossKev">260401Cl 1 回の非弾性散乱あたりの平均エネルギー損失 &lt;ΔE&gt; = |dE/ds|·λ_in [keV]</param>
    /// <param name="TotalRate">260401Cl 弾性+非弾性の全散乱レート 1/λ_el + 1/λ_in [1/nm]。ホットループ内の除算を事前計算で除去</param>
    /// <param name="InverseTotalRate">260603Cl 全散乱レートの逆数 1/TotalRate [nm]。ステップ長 s = -ln(R)·InverseTotalRate でホットループ内の除算を乗算に置換</param>
    /// <param name="ElasticProbability">260401Cl 散乱イベントが弾性である確率 = ElasticRate / TotalRate。ホットループ内の乗算を除去</param>
    /// <param name="NearestNistElasticEnergyIndex">260401Cl NIST エネルギー点 (260603Cl 111 点に拡張) 上の最近傍インデックス。ホットループ内の Math.Log を事前計算で除去</param>
    internal readonly record struct TransportParameters( //260922Cl private → internal (tools/EbsdSourceStudy の源定義の研究用。本番の挙動は不変) // (260331Ch) 1 ステップで使う輸送パラメータをまとめて扱う
        double ScreeningParameter, double ElasticCrossSectionNm2, double ElasticMeanFreePathNm, double StoppingPowerKevPerNm, double InelasticMeanFreePathNm, double MeanInelasticLossKev,
        // double TotalRate, double ElasticProbability, int NearestNistElasticEnergyIndex); // 260401Cl 追加: ホットループの除算・Math.Log を排除 // 260603Cl 旧シグネチャ
        double TotalRate, double InverseTotalRate, double ElasticProbability, int NearestNistElasticEnergyIndex); // 260603Cl 追加: InverseTotalRate でステップ長の除算を排除
    #endregion

    /// <summary>コンストラクタ
    /// <para>⚠ <paramref name="z"/> と <paramref name="a"/> は内部で単独ではなく <c>ρ·z/a</c> という**比**
    /// (Bethe 阻止能の電子密度) として使われるので、<b>同じ重み = 原子 1 個あたりで揃えて渡すこと</b>。
    /// 結晶から作るなら <see cref="GetMeanAtomicParameters"/> の Z と A をそのまま渡せばよい。
    /// 単体元素では平均の取り方に依らないので、揃っていなくても露見しない (260921Cl の不具合はこれだった)。</para>
    /// <para>260921Cl: <paramref name="valenceElectronCount"/> と <paramref name="meanIonizationPotentialEv"/> は
    /// <b><paramref name="atoms"/> を渡していれば省略してよい</b> — 組成から正しい規約で導出する。
    /// 明示した場合はそちらが優先される (感度解析用の上書き)。</para></summary>
    /// <param name="z">原子番号 (単位無し)。化合物では**原子数平均** Σn_iZ_i/Σn_i</param>
    /// <param name="a">原子量 (g/mol)。化合物では**原子数平均** Σn_iA_i/Σn_i = 原子 1 個あたりの質量</param>
    /// <param name="_ρ">密度 (g/cm^3)</param>
    /// <param name="kev">入射電子エネルギー (kev)</param>
    /// <param name="tilt">試料表面の傾き (rad, X軸で回転)</param>
    /// <param name="thresholdKev">飛程計算を打ち切るエネルギー (kev)</param>
    /// <param name="stoppingPowerModel">阻止能モデル。null のときは DefaultStoppingPowerModel を使う。</param>
    /// <param name="elasticScatteringModel">弾性散乱モデル。null のときは DefaultElasticScatteringModel を使う。</param>
    /// <param name="inelasticScatteringModel">非弾性散乱モデル。null のときは DefaultInelasticScatteringModel を使う。</param>
    /// <param name="valenceElectronCount">TPP-2M に使う平均価電子数 Nv (原子数平均)。
    ///   null のとき: <paramref name="atoms"/> があれば組成から導出、無ければ <paramref name="z"/> からの簡易推定。</param>
    /// <param name="meanIonizationPotentialEv">260921Cl 追加: 平均イオン化ポテンシャル J [eV]。
    ///   null のとき: <paramref name="atoms"/> があれば **Bragg 則** (電子数重みで ln J を平均) で導出、
    ///   無ければ <paramref name="z"/> から <see cref="ElementIonizationPotentialEv"/>。
    ///   ⚠ J は Z ≤ 12 で式が切り替わる非線形関数なので、化合物で平均 Z から作ると大きく外れる。</param>
    /// <param name="bandGapEv">TPP-2M に使うバンドギャップ Eg (eV)。null のときは 0 eV。</param>
    /// <param name="atoms">Mott/NIST sampler 用の元素組成。null のときは平均 Z 近似のまま扱う。</param>
    public MonteCarlo(
        double z, double a, double _ρ, double kev, double tilt, double thresholdKev = 2,
        StoppingPowerModels stoppingPowerModel = StoppingPowerModels.JablonskiModified2008,
        ElasticScatteringModels elasticScatteringModel = ElasticScatteringModels.MottNistSampler2023,
        InelasticScatteringModels inelasticScatteringModel = InelasticScatteringModels.DiscreteBulkDiimfpApproximation,
        double? valenceElectronCount = null,
        double? bandGapEv = null,
        IEnumerable<Atoms> atoms = null,
        double? meanIonizationPotentialEv = null) //260921Cl 追加 (オーバーロードではなく既定引数)
    {
        Z = z; A = a; ρ = _ρ; 
        InitialKev = kev; Tilt = tilt; ThresholdKev = thresholdKev; 
        StoppingPowerModel = stoppingPowerModel;
        ElasticScatteringModel = elasticScatteringModel;
        InelasticScatteringModel = inelasticScatteringModel;
        #region お蔵入り // (260401Ch) generated / external の source flag 比較経路は配布版では使わない
        // ElasticSamplerDataSource = elasticSamplerDataSource ?? DefaultElasticSamplerDataSource;
        #endregion
        //260921Cl 変更: 組成 atoms があるなら Nv と J は**そこから導く**。旧実装は「呼び出し側が正しい重みで渡す」前提で、
        //  渡し忘れると平均 Z からの簡易推定 (= 化合物では誤り) に黙って落ちていた。これがまさに今回直した不具合の構造なので、
        //  既定を常に正しい側へ倒す。明示引数は上書きとして優先される (感度解析用。tools の --mc-nv / --mc-a)。
        var derived = atoms is null ? ((double, double, double, double)?)null : GetMeanAtomicParameters(atoms);
        // ValenceElectronCount = valenceElectronCount is > 0 ? valenceElectronCount.Value : EstimateValenceElectronCount(z); // (260331Ch) // 260921Cl 変更前
        ValenceElectronCount = valenceElectronCount is > 0 ? valenceElectronCount.Value
            : derived is { Item3: > 0 } d3 ? d3.Item3 : EstimateValenceElectronCount(z); // (260331Ch) / 260921Cl 組成から導く経路を追加
        BandGapEv = bandGapEv is >= 0 ? bandGapEv.Value : 0.0; // (260331Ch)
        ElasticComponents = atoms is null ? [] : BuildElasticSpecies(atoms, ρ); // (260401Ch) 配布版は generated data だけから Mott sampler を構築する
        // 260919Cl 追加: 組成平均の等方 B [nm²] (Atoms.Dsf は nm² 格納)。B 未設定の結晶でも零点振動相当の下限 1e-3 nm² を入れる
        {
            double sumB = 0, sumW = 0;
            if (atoms != null)
                foreach (var atomsItem in atoms)
                {
                    double b = atomsItem.Dsf?.BisoEffective ?? 0; // 選択規則は Dsf.BisoEffective に集約
                    double w = atomsItem.Occ * Math.Max(1, atomsItem.Atom?.Length ?? 1);
                    sumB += b * w; sumW += w;
                }
            MeanDebyeWallerBNm2 = Math.Max(DiffuseScatteringFactor.BisoFloorNm2, sumW > 0 ? sumB / sumW : 0);
        }
        for (int i = 0; i < ElasticComponents.Length; i++)
            TotalElasticNumberDensityPerNm3 += ElasticComponents[i].NumberDensityPerNm3;

        var nv = Math.Max(ValenceElectronCount, 0.1);
        var eg = Math.Max(BandGapEv, 0.0);
        var u = nv * ρ / A;
        TppPlasmaEnergyEv = 28.8 * Math.Sqrt(u); // (260331Ch)
        TppBeta = -0.10 + 0.944 / Math.Sqrt(TppPlasmaEnergyEv * TppPlasmaEnergyEv + eg * eg) + 0.069 * Math.Pow(ρ, 0.1); // (260331Ch)
        TppGamma = 0.191 * Math.Pow(ρ, -0.5); // (260331Ch)
        TppC = 1.97 - 0.91 * u; // (260331Ch)
        TppD = 53.4 - 20.8 * u; // (260331Ch)
        // 260919Cl 追加: プラズモンカットオフ q_c = ω_p / v_F。価電子密度 n_v [nm⁻³] = Nv·ρ·N_A/A、k_F = (3π² n_v)^{1/3}、
        //   q_c = E_p / ((ħ²/m_e)·k_F) (ħ²/m_e = 0.0761996 eV·nm²、物理慣例 rad/nm)。Si で ≈ 12 rad/nm (1.2 Å⁻¹)。InelasticDecoherenceProbability の q は k=1/λ 系なので 2π で割って格納 (Si ≈ 1.9 nm⁻¹)。
        {
            double atomsPerNm3 = A > 0 ? ρ / A * UniversalConstants.A * 1E-21 : 0; // ρ [g/cm³]·N_A/A → 1/nm³ (/simplify2: リテラルを UniversalConstants.A へ)
            double kF = Math.Cbrt(3 * Math.PI * Math.PI * nv * atomsPerNm3);
            InelasticLocalizationQNm = kF > 0 ? TppPlasmaEnergyEv / (0.0761996 * kF) / (2 * Math.PI) : 0; // (/simplify2) k_F が作れないときは q_c=0 = 常に局在 (従来どおり非弾性は必ずリセット)。+∞ だと逆に「絶対にリセットしない」になる // 2π で割り、k=1/λ (結晶学慣例) の q² = k²(θ²+θ_E²) と同じ単位にする
            double minLossEv = Math.Max(BandGapEv > 0 ? BandGapEv : 1.0, 0.5);
            InelasticLocalizedLossEv = Math.Max(1.8 * Math.Max(TppPlasmaEnergyEv, minLossEv + 0.5), 30.0); // CreateBulkLossSamplerEntry の tailOnsetEv と同じ定義
        }

        //散乱係数の計算中に出てくる定数
        coeff0 = 0.0034 * Math.Pow(Z, 2.0 / 3.0);
        //トータル散乱断面積の計算中に出てくる定数
        coeff1 = Z * UniversalConstants.e0 * UniversalConstants.e0 / (8.0 * Math.PI * UniversalConstants.ε0);
        //平均自由行程のところに出てくる定数
        coeff2 = A / UniversalConstants.A / ρ * 1E21;// / Math.PI;
        //阻止能の計算中に出てくる定数
        coeff3 = -Z * UniversalConstants.A * ρ * 1E3 / (A * 1E-3) * Math.Pow(UniversalConstants.e0, 4)
            / 4 / Math.PI / UniversalConstants.ε0 / UniversalConstants.ε0 / UniversalConstants.eV_joule * 1E-9 * 1E-3;

        //阻止能の計算中に出てくる物質依存の定数 k    Joy and Luo (1989)によれば 6C:0.77, 13Al: 0.815, 14Si: 0.822, 28Ni: 0.83, 29Cu: 0.83,  47Ag:0.852, 79Au: 0.851
        //取りあえず対数近似した値を使う
        k = 0.0299 * Math.Log(Z) + 0.7307;
        //阻止能の計算中に出てくる物質依存の定数 J (eV) Z<=12の時は J=11.5evにするらしい (Joy&Luo 1989)
        // J = Z <= 12 ? 11.5 * Z : (9.76 * Z + 58.5 / Math.Pow(Z, 0.19)); // 260921Cl 変更前 (化合物でも平均 Z から作っていた)
        //260921Cl 変更: 化合物は Bragg 則の J を使う。式は ElementIonizationPotentialEv に一本化。
        //  ⚠ 平均 Z からの ElementIonizationPotentialEv(Z) へ落ちるのは、組成も明示値も無いとき (= 単体元素) だけ。
        //  J は Z ≤ 12 で式が切り替わる非線形関数なので、化合物で平均 Z から作ると大きく外れる
        J = meanIonizationPotentialEv is > 0 ? meanIonizationPotentialEv.Value
            : derived is { Item4: > 0 } d4 ? d4.Item4 : ElementIonizationPotentialEv(Z);
        tan = Math.Tan(tilt);
        (sin, cos) = Math.SinCos(tilt);
        MottElasticMixtureCache = BuildMottElasticMixtureCache(); // (260331Ch)
        TransportParameterCache = BuildTransportParameterCache(); // (260331Ch)
        BulkLossSamplerCache = BuildBulkLossSamplerCache(); // (260331Ch) 簡易 bulk DIIMFP sampler
        if (InelasticScatteringModel == InelasticScatteringModels.DiscreteDrudeValenceInnerShell) //260922Cl 追加
            (DrudeLossSamplerEnergiesKev, DrudeLossSamplerCache) = BuildDrudeLossSamplerCache();
    }

    /// <summary>混合物の各元素の原子番号と重みから、TPP-2M に使う平均価電子数 Nv を重み付き平均で求める。
    /// ⚠ <paramref name="components"/> の Weight は**原子数** (占有率×多重度) であって質量比ではない。260921Cl 規約を明記。
    /// 【理由】Nv は <see cref="MonteCarlo"/> 内で必ず u = Nv·ρ/A という**比**として使われ (TPP-2M)、
    ///   この u は「単位質量あたりの価電子数 × 密度」= N_v/M·ρ でなければならない。A が原子数平均
    ///   (= 原子 1 個あたりの質量) なので、Nv も原子数平均 (= 原子 1 個あたりの価電子数) で揃える必要がある。
    ///   質量重みにすると Nv/A ≠ N_v/M になり、水素を含む化合物ほど大きく外れる
    ///   (Cu₂(OH)₃Cl で 1.63 倍、Mg(OH)₂ で 1.30 倍。単体元素では両者が一致するので露見しなかった)。
    /// 旧: 「質量比から…質量加重平均で推定する」(260331Ch〜260921Cl 手前)</summary>
    public static double EstimateAverageValenceElectronCount(IEnumerable<(int AtomicNumber, double Weight)> components)
    {
        double sum = 0, weightSum = 0;
        foreach (var (atomicNumber, weight) in components)
        {
            if (!(weight > 0))
                continue;
            sum += EstimateElementValenceElectronCount(atomicNumber) * weight; // 260921Cl 式は不変。変わったのは呼び出し側が渡す weight の規約 (質量 → 原子数)
            weightSum += weight;
        }
        return weightSum > 0 ? sum / weightSum : double.NaN;
    }

    /// <summary>260921Cl 追加: Joy &amp; Luo (1989) の平均イオン化ポテンシャル J [eV] を元素 1 つについて返す。
    /// Z ≤ 12 では J = 11.5·Z、それ以外は J = 9.76·Z + 58.5·Z^(−0.19)。
    /// ⚠ コンストラクタ (単体元素) と <see cref="GetMeanAtomicParameters"/> の Bragg 則 (化合物) の
    /// **両方から呼ぶ**。式が 2 箇所に分かれて食い違うのを防ぐために切り出してある (1 行だがインライン化しないこと)。</summary>
    public static double ElementIonizationPotentialEv(double z) => z <= 12 ? 11.5 * z : 9.76 * z + 58.5 / Math.Pow(z, 0.19);

    /// <summary>結晶の原子リスト(非対称単位)から、MC 計算に渡す物質パラメータを Multiplicity×Occ 加重で求める
    /// (Crystal.GetFormulaAndDensity と同じ規約)。260612Cl 追加 / 260921Cl 重み付けを修正・J を追加。
    /// (FormTrajectory/FormEBSD で重複していた sum1/sum2/sum3 集計を集約)
    ///
    /// <para>⚠⚠ <b>4 つの量はすべて「原子 1 個あたり」で揃える。</b> 260921Cl までは Z と Nv だけ質量加重で、
    /// A だけ原子数平均という不整合があった。<see cref="MonteCarlo"/> は受け取った 3 つを
    /// <c>ρ·Z/A</c> (Bethe 阻止能) と <c>Nv·ρ/A</c> (TPP-2M) という**比**として使うので、
    /// 重みが揃っていないと比そのものが狂う。単体元素では 3 つの平均が一致するため露見していなかった。
    /// 水素を含む化合物ほどひどく、Cu₂(OH)₃Cl では Z/A が 0.9226 (正しくは 0.4776 = 1.93 倍)、
    /// Nv·ρ/A が 1.63 倍だった。無水のケイ酸塩・酸化物は 1.06〜1.08 倍。</para>
    ///
    /// <para>⚠ J だけは平均 Z から作ってはいけない。<see cref="ElementIonizationPotentialEv"/> は Z に対して
    /// 非線形 (Z ≤ 12 で式が切り替わる) なので、元素ごとに J_i を出して**電子数重みで ln J を平均する**
    /// (Bragg 則)。Cu₂(OH)₃Cl では Bragg 則が 197.9 eV、平均 Z からだと 130 eV (原子数平均 11.33) か
    /// 246 eV (質量加重 21.89) になってどちらも外れる。</para></summary>
    /// <returns>Z = 原子数平均の原子番号、A = 原子数平均の原子量 [g/mol] (= 原子 1 個あたりの質量)、
    /// Nv = 原子数平均の価電子数、JEv = Bragg 則の平均イオン化ポテンシャル [eV]</returns>
    // 旧: public static (double Z, double A, double Nv) GetMeanAtomicParameters(IEnumerable<Atoms> atoms) // 260921Cl 手前 (Z/Nv が質量加重・J 無し)
    public static (double Z, double A, double Nv, double JEv) GetMeanAtomicParameters(IEnumerable<Atoms> atoms)
    {
        // 旧 (260612Cl〜260921Cl 手前): Z と Nv を質量重み w = A_i·n_i で、A だけ原子数 n_i で平均していた
        // double sumWZ = 0, sumW = 0, sumN = 0, sumNv = 0, sumNvW = 0;
        // foreach (var atom in atoms)
        // {
        //     var count = atom.Multiplicity * atom.Occ;
        //     var w = AtomStatic.AtomicWeight(atom.AtomicNumber) * count;
        //     sumWZ += w * atom.AtomicNumber;
        //     sumW += w;
        //     sumN += count;
        //     if (w > 0) { sumNv += EstimateElementValenceElectronCount(atom.AtomicNumber) * w; sumNvW += w; }
        // }
        // return (sumWZ / sumW, sumW / sumN, sumNvW > 0 ? sumNv / sumNvW : double.NaN);
        double sumN = 0, sumNZ = 0, sumNA = 0, sumNNv = 0, sumNZlnJ = 0;
        int onlyZ = -1; //単体元素の判定用 (-1 = まだ無し、-2 = 複数種)
        foreach (var atom in atoms)
        {
            var count = atom.Multiplicity * atom.Occ; //原子数 (占有率込み)
            if (!(count > 0)) continue;
            //260921Cl 追加: Z が 1..99 の外 (既定値のまま作られた原子・ダミーサイト・壊れた XML) を弾く。
            //  ⚠ 弾かないと ElementIonizationPotentialEv(0) = 0 → Math.Log(0) = -∞、count*z = 0 で
            //  0 × (-∞) = NaN となり、Bragg 則の J が NaN → 化合物なのに平均 Z 経路へ無警告で落ちる。
            //  Crystallography.Controls の AggregateElements も同じ範囲で弾いており、規約はこちらに揃える
            if (atom.AtomicNumber < 1 || atom.AtomicNumber > 99) continue;
            double z = atom.AtomicNumber;
            sumN += count;
            sumNZ += count * z;
            sumNA += count * AtomStatic.AtomicWeight(atom.AtomicNumber);
            sumNNv += count * EstimateElementValenceElectronCount(atom.AtomicNumber);
            sumNZlnJ += count * z * Math.Log(ElementIonizationPotentialEv(z)); //Bragg 則は電子数 n_i·Z_i 重み
            onlyZ = onlyZ == -1 || onlyZ == atom.AtomicNumber ? atom.AtomicNumber : -2;
        }
        if (!(sumN > 0) || !(sumNA > 0))
            return (double.NaN, double.NaN, double.NaN, double.NaN);
        //⚠ 元素 1 種類のときは exp(ln J) の往復を通さない。往復すると 1〜2 ulp ずれて
        //  「単体元素では値が一切変わらない」という保証が崩れる (Si/Fe/Au の既存結果を動かさないため)
        double jEv = onlyZ >= 0 ? ElementIonizationPotentialEv(onlyZ)
                   : sumNZ > 0 ? Math.Exp(sumNZlnJ / sumNZ) : double.NaN;
        return (sumNZ / sumN, sumNA / sumN, sumNNv / sumN, jEv);
    }

    /// <summary>
    /// 結晶の原子リストから Mott 散乱サンプリング用の元素別数密度 n_i [1/nm³] を構築する。
    /// n_i = ρ·N_A·(多重度×占有率)_i / (式量) で計算。
    /// </summary>
    private static ElasticSpecies[] BuildElasticSpecies(IEnumerable<Atoms> atoms, double density)
    {
        if (!(density > 0))
            return [];

        var counts = new Dictionary<int, double>();
        double totalWeightPerFormula = 0.0;
        foreach (var atom in atoms)
        {
            var count = atom.Multiplicity * atom.Occ;
            if (!(count > 0) || atom.AtomicNumber <= 0)
                continue;

            if (!counts.TryAdd(atom.AtomicNumber, count))
                counts[atom.AtomicNumber] += count;
            totalWeightPerFormula += AtomStatic.AtomicWeight(atom.AtomicNumber) * count;
        }

        if (!(totalWeightPerFormula > 0) || counts.Count == 0)
            return [];

        var species = new ElasticSpecies[counts.Count];
        int index = 0;
        foreach (var (atomicNumber, count) in counts)
        {
            var numberDensityPerNm3 = density * UniversalConstants.A * count / totalWeightPerFormula * 1E-21; // (260331Ch) ρ[g/cm3] から元素ごとの数密度 [1/nm3] を作る
            #region お蔵入り // (260401Ch) source flag 比較経路
            // species[index++] = new ElasticSpecies(atomicNumber, numberDensityPerNm3, GetNistElasticScatteringTable(atomicNumber, elasticSamplerDataSource));
            #endregion
            species[index++] = new ElasticSpecies(atomicNumber, numberDensityPerNm3, GetNistElasticScatteringTable(atomicNumber)); // (260401Ch) 配布版は generated data のみを使う
        }
        return species;
    }

    /// <summary>
    /// NIST テーブルの 111 エネルギー点 (260603Cl 拡張) ごとに、混合物系の巨視的弾性散乱断面積 Σ_total と元素選択 CDF を事前計算する。
    /// シミュレーション中のイベントごとの再計算を避けるためのキャッシュ。
    /// </summary>
    private MottElasticMixtureEntry[] BuildMottElasticMixtureCache()
    {
        if (ElasticScatteringModel != ElasticScatteringModels.MottNistSampler2023 ||
            ElasticComponents.Length == 0 || !(TotalElasticNumberDensityPerNm3 > 0))
            return [];

        // var entries = new MottElasticMixtureEntry[101]; // 260603Cl 変更前: テーブル拡張時にハードコード 101 が残ると GetLowerNistElasticEnergyIndex の返す index(最大 EnergyCount-2)で範囲外
        var entries = new MottElasticMixtureEntry[NistElasticEnergyCount]; // 260603Cl 111 点に追従 (SigmaA0Squared.Length と一致させる)
        var macroscopicCrossSections = new double[ElasticComponents.Length];
        for (int energyIndex = 0; energyIndex < entries.Length; energyIndex++)
        {
            double totalMacroscopicCrossSection = 0.0;
            for (int i = 0; i < ElasticComponents.Length; i++)
            {
                // var table = GetNistElasticScatteringTable(ElasticComponents[i].AtomicNumber); // (260331Ch) 旧実装
                var table = ElasticComponents[i].NistElasticTable;
                if (table is null)
                    return [];

                var sigmaNm2 = table.SigmaA0Squared[energyIndex] * NistElasticCrossSectionUnitNm2;
                var macroscopicCrossSection = ElasticComponents[i].NumberDensityPerNm3 * sigmaNm2;
                macroscopicCrossSections[i] = macroscopicCrossSection;
                totalMacroscopicCrossSection += macroscopicCrossSection;
            }

            if (!(totalMacroscopicCrossSection > 0))
                return [];

            var cumulativeProbabilities = new double[ElasticComponents.Length];
            double partial = 0.0;
            for (int i = 0; i < ElasticComponents.Length; i++)
            {
                partial += macroscopicCrossSections[i];
                cumulativeProbabilities[i] = partial / totalMacroscopicCrossSection;
            }
            cumulativeProbabilities[^1] = 1.0; // (260331Ch) 丸め誤差で 1 未満になるのを防ぐ
            entries[energyIndex] = new MottElasticMixtureEntry(totalMacroscopicCrossSection, cumulativeProbabilities);
        }
        return entries;
    }

    /// <summary>
    /// 1 eV 刻みで TransportParameters を事前計算し配列に格納する。
    /// シミュレーション中は配列インデックスアクセスだけで α, σ_el, λ_el, dE/ds, λ_in, &lt;ΔE&gt; を取得でき、毎ステップの再計算を回避する。
    /// </summary>
    private TransportParameters[] BuildTransportParameterCache()
    {
        int maxEnergyEv = Math.Max(1, (int)Math.Ceiling(InitialKev * 1000.0));
        var cache = new TransportParameters[maxEnergyEv + 1];
        cache[0] = ComputeTransportParameters(0.001); // (260331Ch) 0 eV は log の都合で扱いづらいので 1 eV 相当を入れる
        for (int energyEv = 1; energyEv < cache.Length; energyEv++)
            cache[energyEv] = ComputeTransportParameters(energyEv * 0.001); // (260331Ch) 1 eV 刻みなら近似誤差は十分小さい
        return cache;
    }

    /// <summary>
    /// DiscreteBulkDiimfpApproximation 用のエネルギー損失分布サンプラーを 10 eV 刻みで事前構築する。
    /// 各エネルギーでプラズモンピーク・低損失・高損失テールの 3 成分混合 PDF → CDF を作成。
    /// </summary>
    private BulkLossSamplerEntry[] BuildBulkLossSamplerCache()
    {
        const int energyStepEv = 10; // (260331Ch) 10 eV 刻みなら十分軽く、loss spectrum の変化も追いやすい
        int maxEnergyEv = Math.Max(energyStepEv, (int)Math.Ceiling(InitialKev * 1000.0));
        int entryCount = maxEnergyEv / energyStepEv + 1;
        var cache = new BulkLossSamplerEntry[entryCount];
        for (int entryIndex = 0; entryIndex < cache.Length; entryIndex++)
        {
            int energyEv = Math.Max(1, entryIndex * energyStepEv);
            var parameters = TransportParameterCache[Math.Min(energyEv, TransportParameterCache.Length - 1)];
            cache[entryIndex] = CreateBulkLossSamplerEntry(energyEv, parameters);
        }
        return cache;
    }

    /// <summary>
    /// 指定エネルギーでの非弾性エネルギー損失分布を構築する。
    /// 低損失成分 (∝ 1/(ω+0.35Ep)^1.6)、プラズモン成分 (ガウシアン @ Ep)、高損失テール (∝ 1/ω^2.2) の
    /// 3 成分を平均損失 &lt;ΔE&gt; を再現するように重み付け混合し、逆関数法用の CDF を返す。
    /// </summary>
    private BulkLossSamplerEntry CreateBulkLossSamplerEntry(int energyEv, TransportParameters parameters)
    {
        const int binCount = 256;
        double minLossEv = Math.Max(BandGapEv > 0 ? BandGapEv : 1.0, 0.5);
        double meanLossEv = parameters.MeanInelasticLossKev * 1000.0;
        if (!(meanLossEv > 0) || energyEv <= minLossEv + 1.0)
            return null;

        double ep = Math.Max(TppPlasmaEnergyEv, minLossEv + 0.5);
        double maxLossEv = Math.Min(energyEv - 0.5, Math.Max(Math.Max(12.0 * ep, 6.0 * meanLossEv), 250.0));
        if (!(maxLossEv > minLossEv))
            return null;

        double lossStepEv = (maxLossEv - minLossEv) / binCount;
        if (!(lossStepEv > 0))
            return null;

        var lowLossPdf = new double[binCount];
        var plasmonPdf = new double[binCount];
        var tailPdf = new double[binCount];
        double plasmonSigmaEv = Math.Max(1.5, 0.18 * ep);
        double tailOnsetEv = Math.Min(maxLossEv, Math.Max(1.8 * ep, 30.0));
        for (int i = 0; i < binCount; i++)
        {
            double lossEv = minLossEv + (i + 0.5) * lossStepEv;
            lowLossPdf[i] = 1.0 / Math.Pow(lossEv + 0.35 * ep, 1.6);
            double z = (lossEv - ep) / plasmonSigmaEv;
            plasmonPdf[i] = Math.Exp(-0.5 * z * z);
            tailPdf[i] = lossEv >= tailOnsetEv ? 1.0 / Math.Pow(lossEv, 2.2) : 0.0;
        }

        NormalizePdf(lowLossPdf);
        NormalizePdf(plasmonPdf);
        NormalizePdf(tailPdf);

        double muLow = MeanLossEv(lowLossPdf, minLossEv, lossStepEv);
        double muPlasmon = MeanLossEv(plasmonPdf, minLossEv, lossStepEv);
        double muTail = MeanLossEv(tailPdf, minLossEv, lossStepEv);
        var combinedPdf = new double[binCount];
        if (muTail <= muPlasmon + 1e-9 || meanLossEv <= muPlasmon)
        {
            double wLow = Clamp01((muPlasmon - meanLossEv) / Math.Max(muPlasmon - muLow, 1e-9));
            double wPlasmon = 1.0 - wLow;
            CombinePdfs(combinedPdf, lowLossPdf, wLow, plasmonPdf, wPlasmon, tailPdf, 0.0);
        }
        else
        {
            double wLow = Math.Clamp(0.12 + 0.18 * muPlasmon / Math.Max(meanLossEv, muPlasmon), 0.08, 0.30);
            double remainingMean = (meanLossEv - wLow * muLow) / Math.Max(1.0 - wLow, 1e-9);
            double t = Clamp01((remainingMean - muPlasmon) / Math.Max(muTail - muPlasmon, 1e-9));
            double wTail = (1.0 - wLow) * t;
            double wPlasmon = 1.0 - wLow - wTail;
            CombinePdfs(combinedPdf, lowLossPdf, wLow, plasmonPdf, wPlasmon, tailPdf, wTail);
        }

        NormalizePdf(combinedPdf);
        var cumulativeProbabilities = new double[binCount];
        double cumulative = 0.0;
        for (int i = 0; i < binCount; i++)
        {
            cumulative += combinedPdf[i];
            cumulativeProbabilities[i] = cumulative;
        }
        cumulativeProbabilities[^1] = 1.0;
        // 260401Cl 追加: CDF バイナリサーチを O(1) に高速化するガイドテーブル構築
        const int guideSize = 64;
        var guide = new byte[guideSize];
        int bin = 0;
        for (int k = 0; k < guideSize; k++)
        {
            double v = k / (double)(guideSize - 1);
            while (bin < cumulativeProbabilities.Length - 1 && cumulativeProbabilities[bin] < v)
                bin++;
            guide[k] = (byte)Math.Max(0, bin > 0 ? bin - 1 : 0);
        }
        // return new BulkLossSamplerEntry(minLossEv * 0.001, lossStepEv * 0.001, cumulativeProbabilities); // 260401Cl 旧
        return new BulkLossSamplerEntry(minLossEv * 0.001, lossStepEv * 0.001, cumulativeProbabilities, guide); // 260401Cl GuideTable 追加
    }

    /// <summary>3 つの確率密度関数を重み付き線形結合する。destination[i] = w1·pdf1[i] + w2·pdf2[i] + w3·pdf3[i]</summary>
    private static void CombinePdfs(double[] destination, double[] pdf1, double weight1, double[] pdf2, double weight2, double[] pdf3, double weight3)
    {
        for (int i = 0; i < destination.Length; i++)
            destination[i] = weight1 * pdf1[i] + weight2 * pdf2[i] + weight3 * pdf3[i];
    }

    /// <summary>確率密度関数の配列を合計 1 に正規化する。</summary>
    private static void NormalizePdf(double[] pdf)
    {
        double sum = 0.0;
        for (int i = 0; i < pdf.Length; i++)
            sum += pdf[i];
        if (!(sum > 0))
            return;
        for (int i = 0; i < pdf.Length; i++)
            pdf[i] /= sum;
    }

    /// <summary>離散化された PDF から期待値 Σ pdf[i]·(minLoss + (i+0.5)·step) を計算する [eV]。</summary>
    private static double MeanLossEv(double[] pdf, double minLossEv, double lossStepEv)
    {
        double sum = 0.0;
        for (int i = 0; i < pdf.Length; i++)
            sum += pdf[i] * (minLossEv + (i + 0.5) * lossStepEv);
        return sum;
    }

    private static double Clamp01(double value)
        => Math.Clamp(value, 0.0, 1.0);

    /// <summary>260922Cl 追加: <see cref="InelasticScatteringModels.DiscreteDrudeValenceInnerShell"/> の 1 エネルギー点ぶんのサンプラー。</summary>
    private sealed class DrudeLossSamplerEntry(double[] valenceLossEv, double[] valenceCdf, double[] shellEdgeEv, double[] shellCdf, double coreProbability, double omegaMaxEv)
    {
        /// <summary>価電子の損失の区間境界 [eV] (対数等間隔、長さ = ValenceCdf.Length + 1)</summary>
        public readonly double[] ValenceLossEv = valenceLossEv;
        /// <summary>価電子の損失の累積確率 (区間ごと、末尾 = 1)</summary>
        public readonly double[] ValenceCdf = valenceCdf;
        /// <summary>内殻の副殻の吸収端 [eV]</summary>
        public readonly double[] ShellEdgeEv = shellEdgeEv;
        /// <summary>副殻の累積確率 (数密度 × Bote–Salvat の電離断面積の比、末尾 = 1)。内殻が無ければ空</summary>
        public readonly double[] ShellCdf = shellCdf;
        /// <summary>1 事象が内殻である確率 (1 事象の平均損失 = S·λ_in になるように決めた値)</summary>
        public readonly double CoreProbability = coreProbability;
        /// <summary>損失の上限 E/2 [eV]</summary>
        public readonly double OmegaMaxEv = omegaMaxEv;
    }

    /// <summary>260922Cl 追加: DiscreteDrudeValenceInnerShell の価電子の減衰幅 γ [eV] (Si のプラズモンの幅 3〜5 eV の中央。3/5 eV でも λ_vb・30 eV 以上の割合の変化は ±5 % 程度)</summary>
    private const double DrudeDampingEv = 4.0;

    /// <summary>260922Cl 追加: DiscreteDrudeValenceInnerShell のサンプラーを、対数等間隔の 16 エネルギー点 (max(打ち切り, 0.5 keV)〜E0) で作る。
    /// 内殻の副殻は吸収端が <see cref="InelasticLocalizedLossEv"/> 以上のもの (それ未満は価電子として Drude 側が担う)。
    /// 組成が無いときは平均 Z を丸めた単体元素とみなす。</summary>
    private (double[] energiesKev, DrudeLossSamplerEntry[] entries) BuildDrudeLossSamplerCache()
    {
        const int energyCount = 16;
        double eMax = InitialKev, eMin = Math.Min(Math.Max(ThresholdKev, 0.5), eMax);
        var energies = new double[energyCount];
        for (int i = 0; i < energyCount; i++) energies[i] = eMin * Math.Pow(eMax / eMin, i / (energyCount - 1.0));
        energies[^1] = eMax;

        var shells = new List<(int z, int subshell, double edgeEv, double numberDensity)>();
        var species = new List<(int z, double n)>();
        if (ElasticComponents.Length > 0) foreach (var s in ElasticComponents) species.Add((s.AtomicNumber, s.NumberDensityPerNm3));
        else species.Add(((int)Math.Round(Z), 1.0));
        foreach (var (zi, ni) in species)
        {
            if (zi < 1 || zi > 99 || !(ni > 0)) continue;
            for (int ss = 1; ss <= BoteSalvat.SubshellCount(zi); ss++)
            {
                double edge = BoteSalvat.EdgeEv(zi, ss);
                if (edge >= InelasticLocalizedLossEv) shells.Add((zi, ss, edge, ni));
            }
        }

        var entries = new DrudeLossSamplerEntry[energyCount];
        System.Threading.Tasks.Parallel.For(0, energyCount, i =>
        {
            double eEv = energies[i] * 1000, omegaMax = eEv / 2;
            var (edges, cdfV, meanV) = DrudeValenceDiimfp(eEv, TppPlasmaEnergyEv, DrudeDampingEv);
            //副殻の重み = 数密度 × 電離断面積 (Bote–Salvat)。E/2 が吸収端以下の副殻は除く。損失は [B, E/2] で ∝ ω⁻² (平均 ln(ω_max/B)/(1/B − 1/ω_max))
            var shellEdges = new List<double>(); var shellW = new List<double>(); double meanCoreNum = 0, wSum = 0;
            foreach (var s in shells)
            {
                if (!(omegaMax > s.edgeEv * 1.01)) continue;
                double w = s.numberDensity * BoteSalvat.SigmaCm2(s.z, s.subshell, eEv);
                if (!(w > 0)) continue;
                shellEdges.Add(s.edgeEv); shellW.Add(w);
                meanCoreNum += w * Math.Log(omegaMax / s.edgeEv) / (1 / s.edgeEv - 1 / omegaMax); wSum += w;
            }
            double pCore = 0; var shellCdf = new double[wSum > 0 ? shellW.Count : 0];
            if (wSum > 0)
            {
                double meanCore = meanCoreNum / wSum, target = GetTransportParametersRef(energies[i]).MeanInelasticLossKev * 1000;
                pCore = meanCore > meanV && target > 0 ? Clamp01((target - meanV) / (meanCore - meanV)) : 0;
                double c = 0;
                for (int k = 0; k < shellCdf.Length; k++) { c += shellW[k] / wSum; shellCdf[k] = c; }
                shellCdf[^1] = 1;
            }
            entries[i] = new DrudeLossSamplerEntry(edges, cdfV, [.. shellEdges], shellCdf, pCore, omegaMax);
        });
        return (energies, entries);
    }

    /// <summary>260922Cl 追加: 価電子の拡張 Drude 模型の DIIMFP を、ω の対数等間隔の区間 (0.5 eV〜E/2、600 区間) ごとに積分して
    /// 区間境界・累積確率・平均損失を返す。dλ⁻¹/dω = 1/(π a0 T) ∫_{q−}^{q+} (dq/q) Im[−1/ε(q,ω)] (原子単位)、
    /// Im[−1/ε] = ω_p²γω/((ω²−ω_q²)²+γ²ω²)、ω_q = ω_p + q²/2、T = mv²/2、q± = k(E) ± k(E−ω) (相対論的運動量)。
    /// ω ≥ 10 E_p では Bethe ridge の Lorentzian が q 格子より細くなるので、ridge を解析的に積分した ω_p²/(4Tω(ω−ω_p)) (価電子の二体衝突) を使う。</summary>
    internal static (double[] edgesEv, double[] cdf, double meanEv) DrudeValenceDiimfp(double energyEv, double epEv, double gammaEv)
    {
        const double Ha = 27.211386, a0Nm = 0.0529177, mc2 = 510998.95, hbarcEvNm = 197.32698;
        const int nw = 600, nq = 2000;
        double wp = epEv / Ha, g = gammaEv / Ha;
        double gam = 1 + energyEv / mc2, T = 0.5 * mc2 * (1 - 1 / (gam * gam)) / Ha;
        double K(double eEv) => Math.Sqrt(eEv * (eEv + 2 * mc2)) / hbarcEvNm * a0Nm; //相対論的運動量 [1/a0]
        double k0 = K(energyEv), wMin = 0.5, wMax = energyEv / 2, ridgeFromEv = 10 * epEv;
        var edges = new double[nw + 1];
        for (int i = 0; i <= nw; i++) edges[i] = wMin * Math.Pow(wMax / wMin, (double)i / nw);
        var cdf = new double[nw]; double sum = 0, mom = 0;
        for (int i = 0; i < nw; i++)
        {
            double wEv = Math.Sqrt(edges[i] * edges[i + 1]), w = wEv / Ha, dens;
            if (wEv >= ridgeFromEv) dens = wp * wp / (4 * T * w * (w - wp));
            else
            {
                double k1 = K(energyEv - wEv), lqm = Math.Log(k0 - k1), lqp = Math.Log(k0 + k1), h = (lqp - lqm) / nq, s = 0;
                for (int j = 0; j <= nq; j++)
                {
                    double q = Math.Exp(lqm + j * h), wq = wp + 0.5 * q * q, d = w * w - wq * wq;
                    s += (j == 0 || j == nq ? 0.5 : 1) * wp * wp * g * w / (d * d + g * g * w * w);
                }
                dens = s * h / (Math.PI * T);
            }
            double r = dens * (edges[i + 1] - edges[i]) / Ha; //区間の相対確率
            sum += r; mom += r * wEv; cdf[i] = sum;
        }
        for (int i = 0; i < nw; i++) cdf[i] /= sum;
        cdf[^1] = 1;
        return (edges, cdf, mom / sum);
    }

    /// <summary>260922Cl 追加: DiscreteDrudeValenceInnerShell の 1 事象の損失 [keV] (対数で最も近いエネルギー点のサンプラーを使う)。</summary>
    private double SampleDrudeLossKev(double currentKev, double meanLossKev)
    {
        var energies = DrudeLossSamplerEnergiesKev;
        if (energies.Length == 0) return meanLossKev;
        int i = Array.BinarySearch(energies, currentKev);
        if (i < 0)
        {
            i = ~i;
            if (i >= energies.Length) i = energies.Length - 1;
            else if (i > 0 && currentKev * currentKev < energies[i - 1] * energies[i]) i--; //対数の中点で振り分け
        }
        var entry = DrudeLossSamplerCache[i];
        double lossEv;
        if (entry.ShellCdf.Length > 0 && Rnd.NextDouble() < entry.CoreProbability)
        {
            int s = LowerBound(entry.ShellCdf, Rnd.NextDouble());
            double b = entry.ShellEdgeEv[s], wMax = Math.Max(entry.OmegaMaxEv, b), u = Rnd.NextDouble();
            lossEv = 1 / (1 / b - u * (1 / b - 1 / wMax)); //[B, ω_max] で ∝ ω⁻² の逆関数
        }
        else
        {
            double u = Rnd.NextDouble();
            int j = LowerBound(entry.ValenceCdf, u);
            double c0 = j > 0 ? entry.ValenceCdf[j - 1] : 0, c1 = entry.ValenceCdf[j], f = c1 > c0 ? (u - c0) / (c1 - c0) : 0.5;
            lossEv = entry.ValenceLossEv[j] * Math.Pow(entry.ValenceLossEv[j + 1] / entry.ValenceLossEv[j], f); //区間内は対数で内挿 (必ず区間内 = 正)
        }
        return Math.Min(lossEv * 0.001, currentKev);
    }

    /// <summary>260922Cl 追加: 昇順の累積確率 cdf で cdf[k] ≥ u となる最小の k</summary>
    private static int LowerBound(double[] cdf, double u)
    {
        int k = Array.BinarySearch(cdf, u);
        if (k < 0) k = ~k;
        return Math.Min(k, cdf.Length - 1);
    }

    /// <summary>指定エネルギーでの弾性散乱パラメータ (遮蔽パラメータ α, 断面積 σ [nm²], 平均自由行程 λ [nm], 阻止能 dE/ds [keV/nm]) を返す。</summary>
    public (double ScreeningParameter, double CrossSection, double MeanFreePath, double StoppingPower) GetParameters(double kev)
        => GetParameters(GetTransportParameters(kev)); // (260331Ch)

    private static (double ScreeningParameter, double CrossSection, double MeanFreePath, double StoppingPower) GetParameters(TransportParameters parameters)
        => (parameters.ScreeningParameter, parameters.ElasticCrossSectionNm2, parameters.ElasticMeanFreePathNm, parameters.StoppingPowerKevPerNm); // (260331Ch)

    /// <summary>連続減速近似 (CSDA) の電子飛程 [µm]。阻止能を R = ∫_{cutoff}^{E} dE'/|dE/ds(E')| で積分した経路長。
    /// 後方散乱込みの侵入深さ近似である Kanaya-Okayama とは別物 (CSDA は経路長で常に KO より長め)。
    /// Bethe/Joy-Luo は数百 eV 以下で発散・非物理なので lowerCutoffKeV 以下は積分しない (UI で注記要)。
    /// 被積分 1/|dE/ds| が低 E で大きいので log-E 等間隔格子 + (非等間隔) 台形則。dE/ds が非有限/≧0 の区間は寄与に含めない。260606Cl 追加</summary>
    public double GetCsdaRangeMicron(double keV, double lowerCutoffKeV = 0.2, int samplesPerDecade = 64)
    {
        if (!(keV > lowerCutoffKeV) || !(lowerCutoffKeV > 0)) return double.NaN;
        double logLo = Math.Log10(lowerCutoffKeV), logHi = Math.Log10(keV);
        int n = Math.Max(2, (int)Math.Ceiling((logHi - logLo) * samplesPerDecade));
        double rangeNm = 0.0;
        double prevE = double.NaN, prevG = double.NaN; // g = 1/|dE/ds| [nm/keV]
        for (int i = 0; i <= n; i++)
        {
            double e = Math.Pow(10, logLo + (logHi - logLo) * i / n);
            double sp = GetParameters(e).StoppingPower; // keV/nm (損失なので負)
            double g = double.IsFinite(sp) && sp < 0 ? -1.0 / sp : double.NaN; // nm/keV
            if (i > 0 && double.IsFinite(g) && double.IsFinite(prevG))
                rangeNm += 0.5 * (g + prevG) * (e - prevE); // 台形 (e 格子は非等間隔)
            prevE = e; prevG = g;
        }
        return rangeNm > 0 ? rangeNm / 1000.0 : double.NaN; // nm → µm
    }

    /// <summary>CSDA の質量飛程 [g/cm²] = ρ·R。密度依存を外に出して可搬にする (R[µm]→cm は ×1e-4)。260606Cl 追加</summary>
    public double GetCsdaMassRangeGramPerCm2(double keV, double densityGramPerCm3, double lowerCutoffKeV = 0.2, int samplesPerDecade = 64)
    {
        double rMicron = GetCsdaRangeMicron(keV, lowerCutoffKeV, samplesPerDecade);
        return double.IsFinite(rMicron) ? densityGramPerCm3 * rMicron * 1e-4 : double.NaN; // µm → cm
    }

    /// <summary>元素 z 単体の電子 弾性散乱(全)断面積 σ_el [nm²] を加速電圧 keV で返す。260606Cl 追加
    /// 50eV≤E≤36411eV は NIST SRD64 Mott テーブル(per-element, log-E 線形補間)、範囲外およびデータ無しは screened
    /// Rutherford 近似 (GetParameters の混合物挙動と同じモデル切替。境界で不連続。TEM 域 100-300keV は常に Rutherford)。
    /// z&lt;1 / keV≤0 は NaN。混合物の有効値が要るときは GetParameters を使うこと (本メソッドは元素別表示用)。</summary>
    public static double ElasticCrossSectionNm2(int z, double kev)
    {
        if (z < 1 || !(kev > 0)) return double.NaN;
        double energyEv = kev * 1000.0;
        if (energyEv >= 50.0 && energyEv <= NistElasticMaxEnergyEv)
        {
            var sig = GetNistElasticScatteringTable(z)?.SigmaA0Squared; // NIST Mott (per-element, a₀²)
            if (sig != null && sig.Length == NistElasticEnergyCount)
            {
                int i = GetLowerNistElasticEnergyIndex(energyEv, out double frac);
                double a0sq = sig[i];
                if (frac > 0 && i < sig.Length - 1)
                    a0sq += frac * (sig[i + 1] - sig[i]); // 対数エネルギー軸で線形補間
                if (a0sq > 0)
                    return a0sq * NistElasticCrossSectionUnitNm2; // a₀² → nm²
            }
            // NIST 範囲内だがデータ無し → 下の Rutherford へフォールバック (GetParameters と同挙動)
        }
        return ScreenedRutherfordElasticCrossSectionNm2(z, kev);
    }

    /// <summary>screened Rutherford の単体元素 弾性散乱断面積 σ [nm²] (density 非依存、ComputeTransportParameters と同式)。260606Cl 追加</summary>
    private static double ScreenedRutherfordElasticCrossSectionNm2(int z, double kev)
    {
        double coeff0 = 0.0034 * Math.Pow(z, 2.0 / 3.0);
        double coeff1 = z * UniversalConstants.e0 * UniversalConstants.e0 / (8.0 * Math.PI * UniversalConstants.ε0);
        double mv2 = UniversalConstants.Convert.EnergyToElectronMass(kev) * UniversalConstants.Convert.EnergyToElectronVelositySquared(kev);
        double α = coeff0 / kev;
        double tmp = 2 * coeff1 / mv2;
        double σ = tmp * tmp * Math.PI / α / (α + 1) * 1E18; // nm²
        return double.IsFinite(σ) && σ > 0 ? σ : double.NaN;
    }

    /// <summary>キャッシュから輸送パラメータを取得する。キャッシュ範囲外のエネルギーでは都度計算にフォールバック。</summary>
    private TransportParameters GetTransportParameters(double kev)
    {
        int energyEv = (int)Math.Round(kev * 1000.0);
        if ((uint)energyEv < (uint)TransportParameterCache.Length)
            return TransportParameterCache[energyEv];
        return ComputeTransportParameters(kev);
    }

    /// <summary>
    /// 260603Cl 追加: ホットループ専用の輸送パラメータ参照取得。キャッシュ要素を <c>ref readonly</c> で返し、
    /// 毎イベントの構造体値コピー (8 double + int ≈ 68 B) を回避する。
    /// 飛程ループでは電子エネルギー e は InitialKev から単調減少し常にキャッシュ範囲内に収まるため、
    /// 範囲外は理論上発生しないが、保険として末尾要素にクランプする (都度計算フォールバックは ref で返せないため)。
    /// ビン丸めは値版 GetTransportParameters と同じ Math.Round を使い、旧挙動とビット単位で同一の輸送パラメータを返す。
    /// </summary>
    internal ref readonly TransportParameters GetTransportParametersRef(double kev) // 260603Cl 追加
    {
        int energyEv = (int)Math.Round(kev * 1000.0); // 260603Cl 値版と同一丸め (物理不変を厳密に保つ)
        if ((uint)energyEv >= (uint)TransportParameterCache.Length)
            energyEv = TransportParameterCache.Length - 1;
        return ref TransportParameterCache[energyEv];
    }

    /// <summary>
    /// 指定エネルギーでの全輸送パラメータを計算する。
    /// 弾性散乱: Screened Rutherford (α, σ_el, λ_el) または Mott/NIST (σ_el, λ_el)。
    /// 非弾性散乱: TPP-2M で λ_in を求め、阻止能との整合から平均損失 &lt;ΔE&gt; = |dE/ds|·λ_in を導出。
    /// </summary>
    private TransportParameters ComputeTransportParameters(double kev)
    {
        //電子の質量 (kg) × 電子の速度の2乗 (m^2/s^2)
        var mv2 = UniversalConstants.Convert.EnergyToElectronMass(kev) * UniversalConstants.Convert.EnergyToElectronVelositySquared(kev);
        //散乱係数 / トータル散乱断面積 / 平均自由行程
        double α;
        if (ElasticScatteringModel == ElasticScatteringModels.MottNistSampler2023 &&
            TryGetMottElasticTransport(kev, out double σ_E, out double λ_el))
        {
            α = double.NaN; // (260331Ch) Mott sampler では screening parameter を使わない
        }
        else
        {
            α = coeff0 / kev;
            var tmp = 2 * coeff1 / mv2;
            σ_E = tmp * tmp * Math.PI / α / (α + 1) * 1E18;
            //σ_E = 5.21E-21 * Z * Z / kev / kev * 12.56 / α / (α + 1) * Math.Pow((kev + 511) / (kev + 1022), 2);
            λ_el = coeff2 / σ_E; //λ_el = A / UniversalConstants.A / ρ / σ_E * 1E7;
        }
        //阻止能 (Joy and Luo 1989) (kev/nm単位)
        // var sp = coeff3 / mv2 * Math.Log(1.166 * k + 0.583 / UniversalConstants.eV_joule / J * mv2); // (260331Ch) 旧 Joy-Luo 直書き
        var sp = GetStoppingPower(kev, mv2); // (260331Ch)
        var λ_in = GetInelasticMeanFreePathAngstrom(kev * 1000.0) * 0.1; // (260331Ch) TPP-2M の IMFP [A] -> [nm]
        var meanLossKev = λ_in > 0 && sp < 0 ? -sp * λ_in : double.NaN; // (260331Ch) 平均的には <ΔE>/λ_in = stopping power を満たす
        // 260401Cl 追加: ホットループ内の除算・Math.Log を事前計算で排除
        var elasticRate = λ_el > 0 ? 1.0 / λ_el : 0.0;
        var inelasticRate = λ_in > 0 && meanLossKev > 0 ? 1.0 / λ_in : 0.0;
        var totalRate = elasticRate + inelasticRate;
        var inverseTotalRate = totalRate > 0 ? 1.0 / totalRate : 0.0; // 260603Cl 追加: ステップ長 s = -ln(R)·InverseTotalRate でホットループの除算を排除
        var elasticProbability = totalRate > 0 ? elasticRate / totalRate : 1.0;
        var energyEv = kev * 1000.0;
        var nearestNistIndex = energyEv >= 50.0 && energyEv <= NistElasticMaxEnergyEv ? GetNearestNistElasticEnergyIndex(energyEv) : 0; // 260603Cl 上限 20000→NistElasticMaxEnergyEv(36411eV)
        // return new TransportParameters(α, σ_E, λ_el, sp, λ_in, meanLossKev, totalRate, elasticProbability, nearestNistIndex); // 260603Cl 旧
        return new TransportParameters(α, σ_E, λ_el, sp, λ_in, meanLossKev, totalRate, inverseTotalRate, elasticProbability, nearestNistIndex); // 260603Cl InverseTotalRate 追加
    }

    /// <summary>
    /// Mott/NIST テーブルから弾性散乱全断面積と平均自由行程を取得する。
    /// 対数エネルギー軸上で隣接 2 点の線形補間を行う。適用範囲は 50 eV〜36.4 keV (260603Cl 20keV から拡張)。
    /// </summary>
    private bool TryGetMottElasticTransport(double kev, out double crossSectionNm2, out double meanFreePathNm)
    {
        crossSectionNm2 = meanFreePathNm = double.NaN;
        var energyEv = kev * 1000.0;
        if (MottElasticMixtureCache.Length == 0 || !(TotalElasticNumberDensityPerNm3 > 0) || energyEv < 50.0 || energyEv > NistElasticMaxEnergyEv) // 260603Cl 上限 20000→36411eV
            return false;

        int lowerIndex = GetLowerNistElasticEnergyIndex(energyEv, out var fraction);
        var totalMacroscopicCrossSection = MottElasticMixtureCache[lowerIndex].TotalMacroscopicCrossSectionPerNm;
        if (fraction > 0 && lowerIndex < MottElasticMixtureCache.Length - 1)
            totalMacroscopicCrossSection += fraction * (MottElasticMixtureCache[lowerIndex + 1].TotalMacroscopicCrossSectionPerNm - totalMacroscopicCrossSection); // (260331Ch)
        if (!(totalMacroscopicCrossSection > 0))
            return false;

        crossSectionNm2 = totalMacroscopicCrossSection / TotalElasticNumberDensityPerNm3;
        meanFreePathNm = 1.0 / totalMacroscopicCrossSection;
        return true;
    }

    /// <summary>
    /// 弾性散乱角 cosθ をサンプリングする。Mott/NIST テーブルが利用可能ならそちらを使い、
    /// なければ Screened Rutherford の解析式 cosθ = 1 - 2αR/(1+α-R) (R: 一様乱数) でサンプリング。
    /// </summary>
    // private double SampleElasticScatteringCosTheta(double kev, double α) // 260401Cl 旧シグネチャ
    internal double SampleElasticScatteringCosTheta(double kev, double α, int nistEnergyIndex) // 260401Cl nistEnergyIndex 追加
    {
        if (ElasticScatteringModel == ElasticScatteringModels.MottNistSampler2023 &&
            TrySampleMottElasticCosTheta(kev, nistEnergyIndex, out var cosTheta)) // 260401Cl
            return cosTheta;

        if (double.IsNaN(α))
            α = coeff0 / kev; // 260718Cl: Mott 経路では ScreeningParameter=NaN が渡ってくるため、PCHIP サンプリング失敗時の Rutherford フォールバックで cosθ=NaN → 方向ベクトルが無音で NaN 汚染されていた。非 Mott 分岐 (ComputeTransportParameters) と同式で α を再計算するガードを追加 (正常経路は bit 不変)

        var rnd = Rnd.NextDouble();
        return 1 - 2 * α * rnd / (1 + α - rnd);
    }

    /// <summary>
    /// NIST SRD 64 の累積角度分布 Φ(cosθ) テーブルから弾性散乱角 cosθ を逆関数法でサンプリングする。
    /// 混合物系では巨視的断面積の比で散乱元素を確率的に選択した後、その元素の Φ テーブルを使う。
    /// </summary>
    // private bool TrySampleMottElasticCosTheta(double kev, out double cosTheta) // 260401Cl 旧シグネチャ
    private bool TrySampleMottElasticCosTheta(double kev, int nistEnergyIndex, out double cosTheta) // 260401Cl nistEnergyIndex 追加
    {
        cosTheta = double.NaN;
        if (ElasticComponents.Length == 0 || MottElasticMixtureCache.Length == 0)
            return false;

        var energyEv = kev * 1000.0;
        if (energyEv < 50.0 || energyEv > NistElasticMaxEnergyEv) // 260603Cl 上限 20000→36411eV
            return false;

        // int energyIndex = GetNearestNistElasticEnergyIndex(energyEv); // 260401Cl 事前計算済み
        int energyIndex = nistEnergyIndex; // 260401Cl
        var mixture = MottElasticMixtureCache[energyIndex];
        var choice = Rnd.NextDouble();
        int speciesIndex;
        if (ElasticComponents.Length == 1)
            speciesIndex = 0; // (260331Ch) 単一元素では mixture CDF の binary search を省く
        else
        {
            speciesIndex = Array.BinarySearch(mixture.CumulativeProbabilities, choice);
            if (speciesIndex < 0)
                speciesIndex = ~speciesIndex;
            if (speciesIndex >= ElasticComponents.Length)
                speciesIndex = ElasticComponents.Length - 1;
        }

        // var table = GetNistElasticScatteringTable(ElasticComponents[speciesIndex].AtomicNumber); // (260331Ch) 旧実装
        var table = ElasticComponents[speciesIndex].NistElasticTable;
        if (table is null)
            return false;

        var target = Rnd.NextDouble();
        if (table.GeneratedPchipRuntimeElement is not null)
            return AtomStatic.TryEvaluateGeneratedNistElasticPchipCosTheta(table.GeneratedPchipRuntimeElement, energyIndex, target, out cosTheta); // (260401Ch) generated data の PCHIP 評価は AtomStatic 側へ集約
        #region お蔵入り // (260401Ch) オリジナル TXT の 2001 点 CDF 逆補間は配布版ランタイムでは使わない
        // var phi = table.Phi[energyIndex];
        // if (phi is null)
        //     return false;
        //
        // int upper = Array.BinarySearch(phi, target);
        // if (upper < 0)
        //     upper = ~upper;
        // upper = Math.Clamp(upper, 1, phi.Length - 1);
        // int lower = upper - 1;
        // var phi0 = phi[lower];
        // var phi1 = phi[upper];
        // var x0 = 1.0 - lower * 0.001;
        // var x1 = 1.0 - upper * 0.001;
        // cosTheta = Math.Clamp(phi1 > phi0 ? x0 + (target - phi0) * (x1 - x0) / (phi1 - phi0) : x0, -1.0, 1.0);
        // return true;
        #endregion
        return false; // (260401Ch) generated data が無ければ Mott sampler は失敗として上位で Screened Rutherford にフォールバック
    }

    /// <summary>NIST テーブルの 111 エネルギー点 (260603Cl 拡張) のうち、指定エネルギーに最も近いインデックスを返す (対数軸上で最近傍)。</summary>
    private static int GetNearestNistElasticEnergyIndex(double energyEv)
    {
        var logEnergy = Math.Log(energyEv);
        int lowerIndex = (int)Math.Floor((logEnergy - log50) / LogNistElasticEnergyStep);
        lowerIndex = Math.Clamp(lowerIndex, 0, NistElasticEnergyCount - 2); // 260603Cl 99→EnergyCount-2(=109) に連動
        var lowerEnergy = log50 + lowerIndex * LogNistElasticEnergyStep;
        var upperEnergy = lowerEnergy + LogNistElasticEnergyStep;
        return logEnergy - lowerEnergy < upperEnergy - logEnergy ? lowerIndex : lowerIndex + 1;
    }

    /// <summary>NIST テーブルの対数エネルギー軸上で下側インデックスと補間比率 fraction (0〜1) を返す。線形補間に使用。</summary>
    private static int GetLowerNistElasticEnergyIndex(double energyEv, out double fraction)
    {
        var logEnergy = Math.Log(energyEv);
        int lowerIndex = (int)Math.Floor((logEnergy - log50) / LogNistElasticEnergyStep);
        lowerIndex = Math.Clamp(lowerIndex, 0, NistElasticEnergyCount - 2); // 260603Cl 99→EnergyCount-2(=109) に連動
        var lowerEnergy = log50 + lowerIndex * LogNistElasticEnergyStep;
        var upperEnergy = lowerEnergy + LogNistElasticEnergyStep;
        fraction = (logEnergy - lowerEnergy) / (upperEnergy - lowerEnergy);
        return lowerIndex;
    }

    /// <summary>
    /// 指定原子番号の NIST 弾性散乱テーブルを取得する。
    /// 配布版では generated PCHIP データのみをロードし、static cache に保存する。
    /// </summary>
    private static NistElasticScatteringTable GetNistElasticScatteringTable(int atomicNumber)
    {
        lock (lockObj)
        {
            if (NistElasticSamplerCache.TryGetValue(atomicNumber, out var cached))
                return cached;
        }

        var table = TryLoadGeneratedNistElasticPchipTable(atomicNumber); // (260401Ch) standalone 配布では generated data のみを使う

        lock (lockObj)
            NistElasticSamplerCache[atomicNumber] = table;
        return table;
    }

    #region お蔵入り // (260401Ch) 配布版では外部 TXT 経路を使わない
    // // 260401Cl TryLoadNistElasticBinaryTable メソッドを除去 (BIN 形式読み取りを廃止)
    //
    // /// <summary>NIST SRD 64 オリジナルのテキスト形式 (*.TXT) から弾性散乱テーブルをロードする。PCHIP が利用不可な場合のフォールバック。</summary>
    // private static NistElasticScatteringTable TryLoadNistElasticTextTable(string path)
    // {
    //     if (!File.Exists(path))
    //         return null;
    //
    //     try
    //     {
    //         using var reader = new StreamReader(path);
    //         var table = new NistElasticScatteringTable();
    //         for (int i = 0; i < NistElasticEnergyCount; i++)
    //         {
    //             _ = reader.ReadLine(); // (260331Ch) NIST sampler 内のブロック番号。Fortran 版でも未使用
    //             table.SigmaA0Squared[i] = double.Parse(reader.ReadLine() ?? throw new InvalidDataException(path), CultureInfo.InvariantCulture);
    //             var phi = new double[NistElasticPhiCount];
    //             for (int j = 0; j < phi.Length; j++)
    //                 phi[j] = double.Parse(reader.ReadLine() ?? throw new InvalidDataException(path), CultureInfo.InvariantCulture);
    //             table.Phi[i] = phi;
    //         }
    //         return table;
    //     }
    //     catch
    //     {
    //         return null;
    //     }
    // }
    #endregion

    /// <summary>コード生成済みの PCHIP 補間データから弾性散乱テーブルを構築する。ファイル I/O 不要で最も高速。</summary>
    private static NistElasticScatteringTable TryLoadGeneratedNistElasticPchipTable(int atomicNumber)
    {
        if (!AtomStatic.TryGetGeneratedNistElasticPchipRuntimeElement(atomicNumber, out var runtimeElement) || runtimeElement is null)
            return null;

        if (runtimeElement.SigmaA0Squared is null || runtimeElement.SigmaA0Squared.Length != NistElasticEnergyCount)
            return null;

        var table = new NistElasticScatteringTable()
        {
            AtomicNumber = atomicNumber,
            GeneratedPchipRuntimeElement = runtimeElement,
        };
        Array.Copy(runtimeElement.SigmaA0Squared, table.SigmaA0Squared, NistElasticEnergyCount);
        return table;
    }

    /// <summary>
    /// 選択された阻止能モデルに基づき dE/ds [keV/nm] (負値) を返す。
    /// JoyLuo1989: 修正 Bethe 式。JablonskiModified2008: TPP-2M の IMFP を用いた経験式。
    /// </summary>
    private double GetStoppingPower(double kev, double mv2)
    {
        if (StoppingPowerModel == StoppingPowerModels.JoyLuo1989)
            return coeff3 / mv2 * Math.Log(1.166 * k + 0.583 / UniversalConstants.eV_joule / J * mv2);

        var energyEv = kev * 1000.0;
        var λ_in = GetInelasticMeanFreePathAngstrom(energyEv);
        if (!(λ_in > 0) || double.IsNaN(λ_in) || double.IsInfinity(λ_in))
            return coeff3 / mv2 * Math.Log(1.166 * k + 0.583 / UniversalConstants.eV_joule / J * mv2); // (260331Ch) TPP-2M が破綻したら Joy-Luo に戻す

        var spEvPerAngstrom = StoppingPowerModel switch
        {
            StoppingPowerModels.JablonskiModified2008
                => Jablonski2008D1 * Math.Pow(energyEv, Jablonski2008D2) * Math.Log(Jablonski2008D3 * energyEv) * Math.Pow(ρ + Jablonski2008D4, Jablonski2008D5) / λ_in,
            _ => throw new ArgumentOutOfRangeException(),
        };

        return -spEvPerAngstrom * 0.01; // (260331Ch) 1 eV/Å = 0.01 keV/nm
    }

    /// <summary>
    /// TPP-2M (Tanuma-Powell-Penn) 式で非弾性平均自由行程 λ_in [Å] を計算する。
    /// λ_in(E) = E / {Ep²·[β·ln(γE) - C/E + D/E²]}。Ep はプラズマエネルギー。
    /// </summary>
    // 260605Cl private→public 化 (IMFP universal curve のグラフ表示用)。旧: private double GetInelasticMeanFreePathAngstrom(double energyEv)
    public double GetInelasticMeanFreePathAngstrom(double energyEv)
    {
        var gammaE = TppGamma * energyEv;
        if (!(gammaE > 1.0))
            return double.NaN;

        var denominator = TppPlasmaEnergyEv * TppPlasmaEnergyEv * (TppBeta * Math.Log(gammaE) - TppC / energyEv + TppD / (energyEv * energyEv));
        return denominator > 0 ? energyEv / denominator : double.NaN;
    }

    /// <summary>
    /// 平均原子番号 (非整数可) から TPP-2M 用の価電子数 Nv を推定する。
    /// 既知元素表にあればその値を返し、なければ隣接 2 元素の線形補間。
    /// </summary>
    private static double EstimateValenceElectronCount(double atomicNumber)
    {
        if (AtomStatic.ElementInelasticParameters.TryGetValue((int)Math.Round(atomicNumber), out var exact))
            return exact.ValenceElectrons;

        int lower = (int)Math.Floor(atomicNumber), upper = (int)Math.Ceiling(atomicNumber);
        if (lower == upper)
            return EstimateElementValenceElectronCount(lower);

        var lowerValue = EstimateElementValenceElectronCount(lower);
        var upperValue = EstimateElementValenceElectronCount(upper);
        return lowerValue + (atomicNumber - lower) * (upperValue - lowerValue);
    }

    /// <summary>整数原子番号から価電子数を推定する。既知元素表→周期律に基づく簡易推定の順で決定。</summary>
    private static double EstimateElementValenceElectronCount(int atomicNumber)
    {
        if (AtomStatic.ElementInelasticParameters.TryGetValue(atomicNumber, out var exact))
            return exact.ValenceElectrons;

        return atomicNumber switch
        {
            <= 0 => 1.0,
            1 => 1.0,
            2 => 2.0,
            <= 10 => atomicNumber - 2.0,
            <= 18 => atomicNumber - 10.0,
            <= 36 => EstimatePeriodicValenceElectrons(atomicNumber, 18),
            <= 54 => EstimatePeriodicValenceElectrons(atomicNumber, 36),
            <= 86 => EstimatePeriodicValenceElectrons(atomicNumber, 54),
            _ => EstimatePeriodicValenceElectrons(atomicNumber, 86),
        };
    }

    /// <summary>周期表の族位置から価電子数を簡易推定する。previousNobleGas は直前の希ガスの原子番号。</summary>
    private static double EstimatePeriodicValenceElectrons(int atomicNumber, int previousNobleGas)
    {
        int groupLike = atomicNumber - previousNobleGas;
        return groupLike switch
        {
            <= 2 => groupLike,
            <= 12 => groupLike,
            _ => Math.Min(groupLike - 10, 8),
        };
    }

    /// <summary>
    /// 選択された非弾性散乱モデルに応じて 1 回の非弾性散乱でのエネルギー損失 ΔE [keV] をサンプリングする。
    /// DiscreteMeanLoss: 常に平均値。DiscreteBulkDiimfpApproximation: bulk DIIMFP 近似分布。
    /// </summary>
    /// <summary>
    /// 260919Cl 追加: 弾性散乱 (エネルギー e [keV]、散乱角 cosθ) が熱散漫 (非干渉) 成分か。確率は 1−exp(−2B s²)、s = k·sin(θ/2) → 2s² = k²(1−cosθ)。
    /// 干渉性 Bragg 散乱 (割合 exp(−2B s²)) は Bloch 波側が扱うので、ここでは「イベント無し」と同じ扱い。
    /// x &gt; 20 は確率 1 (乱数を省く)、x &lt; 0.03 は 1−e^{−x} ≈ x(1−x/2) (相対誤差 1.5e-4 未満) で Exp を省く。
    /// </summary>
    [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.AggressiveInlining)]
    internal bool IsThermalDiffuseElasticEvent(double e, double cosθ)
    {
        double x = MeanDebyeWallerBNm2 * UniversalConstants.Convert.EnergyToElectronWaveNumberSquared(e) * (1 - cosθ);
        if (x > 20) return true;
        double pDiffuse = x < 0.03 ? x * (1 - 0.5 * x) : 1 - Math.Exp(-x);
        return Rnd.NextDouble() < pDiffuse;
    }

    /// <summary>
    /// 260919Cl 追加: 損失 lossKev の非弾性散乱が電子のコヒーレンス (Bloch 状態) を壊す確率。
    /// (1) 損失が InelasticLocalizedLossEv 以上なら内殻・単一粒子励起として常に 1 (終状態が原子を特定する)。
    /// (2) それ未満の価電子励起は、運動量移行 q がプラズモンカットオフ q_c を超える割合を返す。角度分布は Bethe ridge θ_c=√(ΔE/E) で
    ///     打ち切った Lorentzian dP ∝ θdθ/(θ²+θ_E²) (θ_E=ΔE/2E)、q² = k²(θ²+θ_E²)。q &lt; q_c の集団励起・非局在の電子正孔対は
    ///     Bloch 状態を保つ (バンド内遷移、エネルギーフィルター実験でプラズモン損失電子にも菊池バンドが残る事実に対応)。
    /// </summary>
    public double InelasticDecoherenceProbability(double energyKev, double lossKev)
    {
        double dE = lossKev * 1000.0, E = energyKev * 1000.0;
        if (!(dE > 0) || !(E > dE)) return 1.0;
        if (dE >= InelasticLocalizedLossEv) return 1.0;
        double k2 = UniversalConstants.Convert.EnergyToElectronWaveNumberSquared(energyKev); // nm⁻² (k=1/λ 系)
        double thetaE2 = dE / (2 * E); thetaE2 *= thetaE2;
        double thetaC2 = dE / E;
        double thetaQ2 = InelasticLocalizationQNm * InelasticLocalizationQNm / k2 - thetaE2;
        if (thetaQ2 <= 0) return 1.0;
        if (thetaQ2 >= thetaC2) return 0.0;
        return 1.0 - Math.Log(1 + thetaQ2 / thetaE2) / Math.Log(1 + thetaC2 / thetaE2);
    }

    internal double SampleInelasticLossKev(double currentKev, double meanLossKev)
    {
        if (!(meanLossKev > 0))
            return 0.0;

        var lossKev = InelasticScatteringModel switch
        {
            InelasticScatteringModels.DiscreteMeanLoss => meanLossKev,
            // 260401Cl DiscreteExponentialLoss の分岐を除去
            InelasticScatteringModels.DiscreteBulkDiimfpApproximation => SampleBulkLossKev(currentKev, meanLossKev),
            InelasticScatteringModels.DiscreteDrudeValenceInnerShell => SampleDrudeLossKev(currentKev, meanLossKev), //260922Cl 追加
            _ => 0.0,
        };
        return Math.Min(lossKev, currentKev);
    }

    /// <summary>
    /// BulkLossSamplerCache の CDF テーブルから逆関数法でエネルギー損失をサンプリングする。
    /// 隣接ビン間の線形補間で連続的な損失値を生成。キャッシュ範囲外では平均損失にフォールバック。
    /// </summary>
    private double SampleBulkLossKev(double currentKev, double meanLossKev)
    {
        if (BulkLossSamplerCache.Length == 0)
            return meanLossKev;

        int energyIndex = (int)Math.Round(currentKev * 100.0); // (260331Ch) 10 eV 刻み cache
        if ((uint)energyIndex >= (uint)BulkLossSamplerCache.Length)
            energyIndex = BulkLossSamplerCache.Length - 1;
        var entry = BulkLossSamplerCache[energyIndex];
        if (entry == null)
            return meanLossKev;

        var target = Rnd.NextDouble();
        // int upper = Array.BinarySearch(entry.CumulativeProbabilities, target); // 260401Cl 旧: O(log 256) バイナリサーチ
        // 260401Cl ガイドテーブルで O(1) に高速化
        int upper;
        var guide = entry.GuideTable;
        if (guide is not null)
        {
            upper = guide[Math.Min((int)(target * (guide.Length - 1)), guide.Length - 1)];
            while (upper < entry.CumulativeProbabilities.Length - 1 && entry.CumulativeProbabilities[upper] < target)
                upper++;
        }
        else
        {
            upper = Array.BinarySearch(entry.CumulativeProbabilities, target);
            if (upper < 0) upper = ~upper;
        }
        upper = Math.Clamp(upper, 1, entry.CumulativeProbabilities.Length - 1);
        int lower = upper - 1;
        var p0 = entry.CumulativeProbabilities[lower];
        var p1 = entry.CumulativeProbabilities[upper];
        var l0 = entry.MinLossKev + (lower + 0.5) * entry.LossStepKev;
        var l1 = entry.MinLossKev + (upper + 0.5) * entry.LossStepKev;
        var lossKev = p1 > p0 ? l0 + (target - p0) * (l1 - l0) / (p1 - p0) : l0;
        return Math.Min(lossKev, currentKev);
    }

    /// <summary>
    /// 電子線の飛程を計算する. 電子は -Z軸に沿って入射し、試料と(0,0,0)の座標で衝突したあと、thresholdで指定したエネルギーまで減衰するか、
    /// 試料表面を脱出するまでの飛程を計算する。返り値は、座標 p (nm単位)と エネルギー e (kev単位) のタプル配列
    /// </summary>
    /// <returns>返り値は、座標 p (nm単位)と エネルギー e (kev単位) のタプル配列</returns>
    public List<(V3 p, double e)> GetTrajectories()
    {
        if (InelasticScatteringModel != InelasticScatteringModels.ContinuousSlowingDownApproximation)
            return GetTrajectoriesDiscreteInelastic(); // (260331Ch) 旧 CSDA 実装は残し、比較できるように分岐する

        var trajectory = new List<(V3 p, double e)>(100) { (new V3(0, 0, 0), InitialKev) };
        //var vec = new V3(0, 0, -1);
        double m11, m12, m13, m22, m23;
        double vX = 0, vY = 0, vZ = -1;
        int n = 0;

        //電子エネルギーがThresholdKev以下になるか、試料を脱出するまでループ
        while (trajectory[^1].e > ThresholdKev && trajectory[^1].p.Y * tan >= trajectory[^1].p.Z)
        {
            //パラメーター取得
            var (α, _, λ_el, sp) = GetParameters(trajectory[^1].e);
            //飛行距離 s
            var s = -λ_el * Math.Log(Rnd.NextDouble());
            if (n++ != 0)
            {
                double rnd3 = Rnd.NextDouble();
                double cosθ = SampleElasticScatteringCosTheta(trajectory[^1].e, α, GetNearestNistElasticEnergyIndex(trajectory[^1].e * 1000.0)), sinθ = Math.Sqrt(1 - cosθ * cosθ); // 260401Cl nistEnergyIndex 追加 (CSDA パスは性能非優先)
                double φ = 2 * Math.PI * rnd3;
                var (sinφ, cosφ) = Math.SinCos(φ);
                double sinθcosφ = sinθ * cosφ, sinθsinφ = sinθ * sinφ;

                var vZ1 = vZ + 1;
                if (vZ1 < Th)
                { vX = sinθcosφ; vY = sinθsinφ; vZ = -cosθ; }
                else
                {
                    m11 = 1 - vX * vX / vZ1;
                    m22 = 1 - vY * vY / vZ1;
                    m12 = -vX * vY / vZ1;
                    m13 = vX;
                    m23 = vY;

                    vX = m11 * sinθcosφ + m12 * sinθsinφ + m13 * cosθ;
                    vY = m12 * sinθcosφ + m22 * sinθsinφ + m23 * cosθ;
                    vZ = -m13 * sinθcosφ - m23 * sinθsinφ + vZ * cosθ;

                    if (n % 10 == 0)
                    {
                        var len = Math.Sqrt(vX * vX + vY * vY + vZ * vZ);
                        vX /= len; vY /= len; vZ /= len;
                    }
                }
            }
            trajectory.Add((trajectory[^1].p + s * new V3(vX, vY, vZ), trajectory[^1].e + s * sp));
        }
        return trajectory;
    }

    /// <summary>
    /// 離散非弾性散乱モデルでの飛程計算。弾性・非弾性を独立なポアソン過程として扱い、
    /// 次イベントまでの距離を s = -ln(R)/Σ_total でサンプリングし、弾性/非弾性を確率的に分岐する。
    /// </summary>
    private List<(V3 p, double e)> GetTrajectoriesDiscreteInelastic()
    {
        var trajectory = new List<(V3 p, double e)>(100) { (new V3(0, 0, 0), InitialKev) };
        double m11, m12, m13, m22, m23;
        double vX = 0, vY = 0, vZ = -1;
        int n = 0;

        while (trajectory[^1].e > ThresholdKev && trajectory[^1].p.Y * tan >= trajectory[^1].p.Z)
        {
            // var parameters = GetTransportParameters(trajectory[^1].e); // 260603Cl 旧: 値コピー (68B)
            ref readonly var parameters = ref GetTransportParametersRef(trajectory[^1].e); // 260603Cl ref readonly で構造体コピー回避
            // if (!(parameters.ElasticMeanFreePathNm > 0)) break; // 260401Cl TotalRate に統合
            // var elasticRate = 1.0 / parameters.ElasticMeanFreePathNm; // 260401Cl 事前計算済み
            // var inelasticRate = parameters.InelasticMeanFreePathNm > 0 && parameters.MeanInelasticLossKev > 0 // 260401Cl
            //     ? 1.0 / parameters.InelasticMeanFreePathNm : 0.0;
            // var totalRate = elasticRate + inelasticRate; // 260401Cl 事前計算済み
            if (!(parameters.TotalRate > 0)) // 260401Cl
                break;

            // var s = -Math.Log(Math.Max(Rnd.NextDouble(), double.Epsilon)) / parameters.TotalRate; // 260401Cl // 260603Cl 旧: 除算
            var s = -Math.Log(Math.Max(Rnd.NextDouble(), double.Epsilon)) * parameters.InverseTotalRate; // 260603Cl 除算→乗算
            var nextPoint = trajectory[^1].p + s * new V3(vX, vY, vZ);
            if (nextPoint.Y * tan < nextPoint.Z)
            {
                trajectory.Add((nextPoint, trajectory[^1].e)); // (260331Ch) 表面を抜けた時点では試料内の次イベントは未発生
                break;
            }

            bool isElastic = parameters.ElasticProbability >= 1.0 || Rnd.NextDouble() < parameters.ElasticProbability; // 260401Cl 事前計算済み
            double nextEnergy = trajectory[^1].e;
            if (isElastic)
            {
                double rnd3 = Rnd.NextDouble();
                double cosθ = SampleElasticScatteringCosTheta(trajectory[^1].e, parameters.ScreeningParameter, parameters.NearestNistElasticEnergyIndex), sinθ = Math.Sqrt(1 - cosθ * cosθ); // 260401Cl nistEnergyIndex 追加
                double φ = 2 * Math.PI * rnd3;
                var (sinφ, cosφ) = Math.SinCos(φ);
                double sinθcosφ = sinθ * cosφ, sinθsinφ = sinθ * sinφ;

                var vZ1 = vZ + 1;
                if (vZ1 < Th)
                { vX = sinθcosφ; vY = sinθsinφ; vZ = -cosθ; }
                else
                {
                    m11 = 1 - vX * vX / vZ1;
                    m22 = 1 - vY * vY / vZ1;
                    m12 = -vX * vY / vZ1;
                    m13 = vX;
                    m23 = vY;

                    vX = m11 * sinθcosφ + m12 * sinθsinφ + m13 * cosθ;
                    vY = m12 * sinθcosφ + m22 * sinθsinφ + m23 * cosθ;
                    vZ = -m13 * sinθcosφ - m23 * sinθsinφ + vZ * cosθ;

                    if (++n % 10 == 0)
                    {
                        var len = Math.Sqrt(vX * vX + vY * vY + vZ * vZ);
                        vX /= len; vY /= len; vZ /= len;
                    }
                }
            }
            else
            {
                nextEnergy = Math.Max(trajectory[^1].e - SampleInelasticLossKev(trajectory[^1].e, parameters.MeanInelasticLossKev), 0.0); // (260331Ch)
            }

            trajectory.Add((nextPoint, nextEnergy));
        }
        return trajectory;
    }

    /// <summary>
    /// 電子線の飛程を計算する. 電子は -Z軸に沿って入射し、試料と(0,0,0)の座標で衝突したあと、thresholdで指定したエネルギーまで減衰するか、
    /// 試料表面を脱出するまでの飛程を計算する。返り値は、座標 p (nm単位), 出射方向 v (単位ベクトル), エネルギー e (kev単位) のタプル
    /// </summary>
    /// <returns>返り値は、深さ (nm単位), 出射方向 (単位ベクトル), エネルギー e (kev単位) のタプル</returns>
    public (double d, V3 v, double e) GetBackscatteredElectrons()
    {
        var electron = GetBackscatteredElectronDetail();
        return (electron.Depth, electron.Direction, electron.Energy); // (260331Ch) 既存呼び出しは壊さず、詳細版へ寄せる
    }

    //260921Cl 修正: この summary は下の GetBackscatteredElectronDetail のもの。
    //  260920Cl に CoherencePreservingAngularVariance を間に挿入した際、doc を奪ってしまっていた
    /// <summary>260920Cl 追加 (調査用): コヒーレンスを壊さなかった非弾性散乱 1 回あたりの角度偏位の二乗平均 ⟨θ²⟩ [rad²]。
    /// 小角非弾性散乱の角度分布を Lorentzian dσ/dΩ ∝ 1/(θ² + θ_E²) とし、θ_E = ΔE/(2E) (特性角)、
    /// 上限を局在カットオフ q_c に対応する θ_c = q_c·λ とする。⟨θ²⟩ = (θ_c² − θ_E²·L)/L、L = ln(1 + θ_c²/θ_E²)。
    /// q &gt; q_c のイベントは別途コヒーレンス破壊として扱うので、ここでは θ ≤ θ_c の裾だけを積む</summary>
    public double CoherencePreservingAngularVariance(double energyKev, double lossKev)
    {
        if (!(lossKev > 0) || !(energyKev > 0) || !(InelasticLocalizationQNm > 0)) return 0;
        double k = UniversalConstants.Convert.EnergyToElectronWaveNumber(energyKev); //1/nm (k = 1/λ)
        if (!(k > 0)) return 0;
        double thetaE = lossKev / (2 * energyKev);
        //260921Cl 変更: θ_c を InelasticDecoherenceProbability (q² = k²(θ² + θ_E²)) と同じ定義に揃える。
        //  旧: thetaC = q_c / k は θ_E を無視しており、同じ q_c に対して 2 つの上限角が存在していた
        //  (Si 20 keV・ΔE 16.7 eV で約 6 % ずれ、低エネルギー側ほど拡大する)
        double tc2 = InelasticLocalizationQNm * InelasticLocalizationQNm / (k * k) - thetaE * thetaE;
        if (!(tc2 > 0) || !(thetaE > 0)) return 0;
        double thetaC = Math.Sqrt(tc2);
        if (thetaC <= thetaE) return 0;
        double L = Math.Log(1 + thetaC * thetaC / (thetaE * thetaE));
        return L > 0 ? (thetaC * thetaC - thetaE * thetaE * L) / L : 0;
    }

    /// <summary>
    /// 後方散乱電子の詳細情報を返す。軌跡座標は保持せず、脱出深さ・方向・エネルギーと
    /// 最後の非弾性散乱の情報のみを記録する軽量版。大量の電子統計を取る EBSD シミュレーション向け。
    /// </summary>
    public BackscatteredElectronDetail GetBackscatteredElectronDetail()
    {
        if (InelasticScatteringModel != InelasticScatteringModels.ContinuousSlowingDownApproximation)
            return GetBackscatteredElectronDetailDiscreteInelastic(); // (260331Ch)

        double e = InitialKev;
        double vX = 0, vY = 0, vZ = -1;
        double d = 0;// 260321Ch: 表面からの深さだけを直接追跡する
        bool hasLastDecoherenceEvent = false; double lastDecoherenceDepth = double.NaN; // 260919Cl 追加
        //260920Cl 追加 (調査用): 最後のコヒーレンス破壊以降の「干渉性区間」で溜まる経路長。破壊イベントのたびに 0 へ戻す
        //260921Cl: CSDA には離散非弾性イベントが無いので、イベント数と角度分散は常に 0。ローカルを置かない
        double coherentPathLength = 0;
        int n = 0;
        //電子エネルギーがThresholdKev以下になるか、試料を脱出するまでループ
        while (e > ThresholdKev)
        {
            //乱数発生
            double rnd1 = Rnd.NextDouble(), rnd3 = Rnd.NextDouble();
            //パラメーター取得
            var (α, _, λ_el, sp) = GetParameters(e);
            //飛行距離 s
            var s = -λ_el * Math.Log(rnd1);
            if (n++ != 0)
            {
                double cosθ = SampleElasticScatteringCosTheta(e, α, GetNearestNistElasticEnergyIndex(e * 1000.0)), sinθ = Math.Sqrt(1 - cosθ * cosθ); // 260401Cl nistEnergyIndex 追加 (CSDA パスは性能非優先)
                if (IsThermalDiffuseElasticEvent(e, cosθ)) { hasLastDecoherenceEvent = true; lastDecoherenceDepth = d; coherentPathLength = 0; } // 260919Cl 追加: 熱散漫成分なら源の深さをリセット / 260920Cl 干渉性区間の累積もリセット
                double φ = 2 * Math.PI * rnd3;
                var (sinφ, cosφ) = Math.SinCos(φ);
                double sinθcosφ = sinθ * cosφ, sinθsinφ = sinθ * sinφ;
                var vZ1 = vZ + 1;
                if (vZ1 < Th)
                { vX = sinθcosφ; vY = sinθsinφ; vZ = -cosθ; }
                else
                {
                    double m11 = 1 - vX * vX / vZ1, m22 = 1 - vY * vY / vZ1, m12 = -vX * vY / vZ1, m13 = vX, m23 = vY;
                    vX = m11 * sinθcosφ + m12 * sinθsinφ + m13 * cosθ;
                    vY = m12 * sinθcosφ + m22 * sinθsinφ + m23 * cosθ;
                    vZ = -m13 * sinθcosφ - m23 * sinθsinφ + vZ * cosθ;

                    if (n % 10 == 0)
                    {
                        var len = Math.Sqrt(vX * vX + vY * vY + vZ * vZ);
                        vX /= len; vY /= len; vZ /= len;
                    }
                }
            }
            var dtmp = d + s * (sin * vY - cos * vZ);// 260321Ch
            if (dtmp < 0)
                break;
            d = dtmp;
            coherentPathLength += s; // 260920Cl 追加 (調査用)

            e += s * sp;
        }
        return new BackscatteredElectronDetail(
            d,
            new V3(vX, vY, vZ),
            e,
            InitialKev - e,
            false,
            double.NaN,
            double.NaN,
            double.NaN,
            new V3(double.NaN, double.NaN, double.NaN), // (260331Ch) CSDA では最後の離散非弾性散乱は定義できない
            hasLastDecoherenceEvent, lastDecoherenceDepth, // 260919Cl 追加
            0, 0, coherentPathLength); // 260920Cl 追加 (調査用) / 260921Cl: CSDA は離散非弾性が無いので前 2 つは常に 0
    }

    /// <summary>
    /// 離散非弾性散乱モデルでの後方散乱電子詳細計算。弾性/非弾性をポアソン過程で分岐し、
    /// 非弾性イベント発生時に深さ・エネルギー・方向を記録して「最後の非弾性散乱」情報を更新する。
    /// </summary>
    private BackscatteredElectronDetail GetBackscatteredElectronDetailDiscreteInelastic()
    {
        double e = InitialKev;
        double vX = 0, vY = 0, vZ = -1;
        double d = 0;
        int n = 0;
        bool hasLastInelasticEvent = false;
        double lastInelasticDepth = double.NaN, lastInelasticEnergyBeforeLoss = double.NaN, lastInelasticEnergyAfterLoss = double.NaN;
        var lastInelasticDirection = new V3(double.NaN, double.NaN, double.NaN);
        bool hasLastDecoherenceEvent = false; double lastDecoherenceDepth = double.NaN; // 260919Cl 追加
        //260920Cl 追加 (調査用): 最後のコヒーレンス破壊以降の「干渉性区間」で溜まる量。破壊イベントのたびに 0 へ戻す
        int coherentInelasticCount = 0; double coherentAngularVariance = 0, coherentPathLength = 0;

        while (e > ThresholdKev)
        {
            // var parameters = GetTransportParameters(e); // 260603Cl 旧: 値コピー (68B)
            ref readonly var parameters = ref GetTransportParametersRef(e); // 260603Cl ref readonly で構造体コピー回避
            // if (!(parameters.ElasticMeanFreePathNm > 0)) break; // 260401Cl TotalRate に統合
            // var elasticRate = 1.0 / parameters.ElasticMeanFreePathNm; // 260401Cl 事前計算済み
            // var inelasticRate = parameters.InelasticMeanFreePathNm > 0 && parameters.MeanInelasticLossKev > 0 // 260401Cl
            //     ? 1.0 / parameters.InelasticMeanFreePathNm : 0.0;
            // var totalRate = elasticRate + inelasticRate; // 260401Cl 事前計算済み
            if (!(parameters.TotalRate > 0)) // 260401Cl
                break;

            // var s = -Math.Log(Math.Max(Rnd.NextDouble(), double.Epsilon)) / parameters.TotalRate; // 260401Cl // 260603Cl 旧: 除算
            var s = -Math.Log(Math.Max(Rnd.NextDouble(), double.Epsilon)) * parameters.InverseTotalRate; // 260603Cl 除算→乗算
            var dtmp = d + s * (sin * vY - cos * vZ); // (260331Ch) 現在の進行方向で次イベント位置まで進む
            if (dtmp < 0)
                break;
            d = dtmp;
            coherentPathLength += s; // 260920Cl 追加 (調査用)

            bool isElastic = parameters.ElasticProbability >= 1.0 || Rnd.NextDouble() < parameters.ElasticProbability; // 260401Cl 事前計算済み
            if (isElastic)
            {
                double rnd3 = Rnd.NextDouble();
                double cosθ = SampleElasticScatteringCosTheta(e, parameters.ScreeningParameter, parameters.NearestNistElasticEnergyIndex), sinθ = Math.Sqrt(1 - cosθ * cosθ); // 260401Cl nistEnergyIndex 追加
                if (IsThermalDiffuseElasticEvent(e, cosθ)) { hasLastDecoherenceEvent = true; lastDecoherenceDepth = d; coherentInelasticCount = 0; coherentAngularVariance = 0; coherentPathLength = 0; } // 260919Cl 追加: 熱散漫成分なら源の深さをリセット / 260920Cl 干渉性区間の累積もリセット
                double φ = 2 * Math.PI * rnd3;
                var (sinφ, cosφ) = Math.SinCos(φ);
                double sinθcosφ = sinθ * cosφ, sinθsinφ = sinθ * sinφ;
                var vZ1 = vZ + 1;
                if (vZ1 < Th)
                { vX = sinθcosφ; vY = sinθsinφ; vZ = -cosθ; }
                else
                {
                    double m11 = 1 - vX * vX / vZ1, m22 = 1 - vY * vY / vZ1, m12 = -vX * vY / vZ1, m13 = vX, m23 = vY;
                    vX = m11 * sinθcosφ + m12 * sinθsinφ + m13 * cosθ;
                    vY = m12 * sinθcosφ + m22 * sinθsinφ + m23 * cosθ;
                    vZ = -m13 * sinθcosφ - m23 * sinθsinφ + vZ * cosθ;

                    if (++n % 10 == 0)
                    {
                        var len = Math.Sqrt(vX * vX + vY * vY + vZ * vZ);
                        vX /= len; vY /= len; vZ /= len;
                    }
                }
            }
            else
            {
                hasLastInelasticEvent = true;
                lastInelasticDepth = d;
                // hasLastDecoherenceEvent = true; lastDecoherenceDepth = d; // 260919Cl 追加 (同日改訂前): 非弾性散乱は常にコヒーレンスを壊す扱いだった
                lastInelasticEnergyBeforeLoss = e;
                lastInelasticDirection = new V3(vX, vY, vZ);
                // e = Math.Max(e - SampleInelasticLossKev(e, parameters.MeanInelasticLossKev), 0.0); // (260331Ch) 260919Cl 変更前
                var lossKev = SampleInelasticLossKev(e, parameters.MeanInelasticLossKev); // 260919Cl 変更: 損失を先に決め、局在判定に使う
                // 260919Cl 追加: 非弾性散乱は運動量移行 (局在) で判定。内殻/高損失は常に、価電子励起は q > q_c の割合だけ源をリセットする
                var pDecoherence = InelasticDecoherenceProbability(e, lossKev);
                if (pDecoherence >= 1 || (pDecoherence > 0 && Rnd.NextDouble() < pDecoherence)) { hasLastDecoherenceEvent = true; lastDecoherenceDepth = d; coherentInelasticCount = 0; coherentAngularVariance = 0; coherentPathLength = 0; } // p=0/1 では乱数を引かない / 260920Cl 干渉性区間の累積もリセット
                //260921Cl 変更: 既定では計算しない (CollectCoherenceDiagnostics の doc 参照)。旧は無条件で Math.Log を呼んでいた
                else if (CollectCoherenceDiagnostics) { coherentInelasticCount++; coherentAngularVariance += CoherencePreservingAngularVariance(e, lossKev); } // 260920Cl 追加 (調査用)
                e = Math.Max(e - lossKev, 0.0);
                lastInelasticEnergyAfterLoss = e;
            }
        }
        return new BackscatteredElectronDetail(
            d,
            new V3(vX, vY, vZ),
            e,
            InitialKev - e,
            hasLastInelasticEvent,
            lastInelasticDepth,
            lastInelasticEnergyBeforeLoss,
            lastInelasticEnergyAfterLoss,
            lastInelasticDirection,
            hasLastDecoherenceEvent, lastDecoherenceDepth, // 260919Cl 追加
            coherentInelasticCount, coherentAngularVariance, coherentPathLength); // 260920Cl 追加 (調査用)
    }

}
