// 260807Cl 新規作成: ALCHEMI の **1D forward orientation engine** (A1′、設計 §5.3 / §5.4)。
//
// BetheMethod の partial として置いてあるのは、方位ループが AccVoltage / Crystal / BaseRotation /
// getEigenMatrix / getPotentialMatrix (いずれも private) を必要とするため。BetheMethod.cs 本体への
// 変更は `class` → `partial class` の 1 語だけで、ロジックは全部このファイルに閉じている
// (ロードマップの「solver 本体は新規ファイルに置き、BetheMethod.cs は薄い結線のみ」)。
//
// 既存 CBED / STEM / EBSD の worker・バッファ・総和順序には一切触れない (設計 §8 の別パス方式)。
// 骨格は cbed_DoWork (:289-576) から流用するが、**円形マスク・√N 正方格子の前提は持たない**
// (ALCHEMI は表面固定・入射傾斜・任意の 1D 配列)。
//
// 流れ:
//   1. チャネル解決 (IonizationDataProvider)               … run 中は immutable
//   2. FixedUnion 基底 (走査の全方位で Find_gVectors して hkl の**真の union**を取る)
//   3. μ を (サイト × チャネル) ぶん **1 回だけ**組む       … μ は方位非依存 (AlchemiMuBuilder)
//   4. 方位ループ: reset_gVectors → getEigenMatrix → EVD → α = C⁻¹e₀ → AlchemiReduction.Yield
//   5. expanded-basis 診断 (代表 3 方位を 1.25 倍基底で再計算して差を測る)
//
// 同期実行の core だけを持つ。BackgroundWorker + イベント (設計 §5.3 の bwALCHEMI) は GUI (A4′) を
// 作るときに薄く被せる — 作者決定の着手順が「物理コア → 検証ツール → GUI」なので、
// まず CLI (AlchemiCheck) から素直に叩ける形にしてある。

using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Threading;
using System.Threading.Tasks;
using static System.Buffers.ArrayPool<System.Numerics.Complex>;//BetheMethod.cs と同じ Shared (using 別名はファイル単位)
using DMat = MathNet.Numerics.LinearAlgebra.Complex.DenseMatrix;
using DVec = MathNet.Numerics.LinearAlgebra.Complex.DenseVector;

namespace Crystallography;

public partial class BetheMethod
{
    #region ALCHEMI (260807Cl 追加)

    /// <summary>260807Cl 追加: 1D forward ALCHEMI シミュレーション (同期実行)。
    /// 既存 worker とは完全に別経路で、CBED/STEM/EBSD のバッファ・数値には触れない。</summary>
    /// <param name="request">run 要求 (呼び出し後に書き換えても影響しないよう内部で必要分をコピーする)</param>
    /// <param name="progress">進捗通知 (null 可)。呼ばれるスレッドは不定</param>
    /// <param name="cancel">キャンセル</param>
    public AlchemiResult RunAlchemi(AlchemiRequest request, Action<AlchemiProgress> progress = null, CancellationToken cancel = default)
    {
        ArgumentNullException.ThrowIfNull(request);
        Validate(request);

        //--- run パラメータの確定 (private setter はここで代入する。設計 §2 の AccVoltage 注記) ---
        AccVoltage = request.IncidentEnergyKeV;
        BaseRotation = new Matrix3D(request.BaseRotation);
        MaxNumOfBloch = request.MaxNumOfBloch;
        Surface = new Vector3D(request.Surface.X, request.Surface.Y, request.Surface.Z);
        Thicknesses = [.. request.ThicknessesNm];
        var surface = new Vector3DBase(request.Surface.X, request.Surface.Y, request.Surface.Z);

        //--- 1. チャネル解決 (run 中は immutable。プロバイダ選択と範囲判定を実行中に持ち込まない) ---
        progress?.Invoke(new AlchemiProgress(AlchemiStage.ResolvingIonizationData, 0));
        var chData = new IonizationData[request.Channels.Length];
        for (int c = 0; c < chData.Length; c++)
            chData[c] = IonizationDataProvider.Resolve(request.Channels[c], request.IncidentEnergyKeV);
        progress?.Invoke(new AlchemiProgress(AlchemiStage.ResolvingIonizationData, 1));
        cancel.ThrowIfCancellationRequested();

        var kvac = UniversalConstants.Convert.EnergyToElectronWaveNumber(AccVoltage);
        var u0 = getU(AccVoltage).Real.Real;

        //--- 2. FixedUnion 基底 ---
        progress?.Invoke(new AlchemiProgress(AlchemiStage.BuildingUnionBasis, 0));
        var (union, diagnostic) = BuildFixedUnion(request, request.Orientations, surface, kvac, u0, request.MaxNumOfBloch);
        progress?.Invoke(new AlchemiProgress(AlchemiStage.BuildingUnionBasis, 1));
        cancel.ThrowIfCancellationRequested();

        //--- 3-4. 本計算 ---
        var shape = new AlchemiTensorShape(request.Orientations.Length, request.ThicknessesNm.Length,
            request.Sites.Length, request.Channels.Length);
        var (dyn, dech, lcoh, mu00) = SolveOrientations(request, chData, union, request.Orientations, shape,
            surface, kvac, u0, progress, cancel);

        //--- 5. expanded-basis 診断 (設計 §5.4: 中心・両端を 1.25 倍基底で再計算) ---
        var expandedMaxRel = double.NaN;
        var accepted = true;
        var warnings = new List<string>(diagnostic.Warnings);
        if (request.ExpandedBasisFactor > 1)
        {
            progress?.Invoke(new AlchemiProgress(AlchemiStage.ExpandedBasisCheck, 0));
            expandedMaxRel = CheckExpandedBasis(request, chData, surface, kvac, u0, dyn, dech, shape, cancel);
            accepted = expandedMaxRel <= request.ExpandedBasisTolerance;
            if (!accepted)
                warnings.Add($"expanded basis ({request.ExpandedBasisFactor:f2}x) changes the total yield by {expandedMaxRel:e2} "
                    + $"(> {request.ExpandedBasisTolerance:e1}) — not accepted for quantitative fitting (設計 §5.4)");
            progress?.Invoke(new AlchemiProgress(AlchemiStage.ExpandedBasisCheck, 1));
        }

        var total = new double[dyn.Length];
        for (int i = 0; i < total.Length; i++) total[i] = dyn[i] + dech[i];

        return new AlchemiResult
        {
            Shape = shape,
            Dynamic = dyn,
            Dechannelled = dech,
            Total = total,
            CoherentPathLengthNm = lcoh,
            Basis = diagnostic with { ExpandedBasisMaxRelDiff = expandedMaxRel, AcceptedForFit = accepted, Warnings = [.. warnings] },
            ChannelData = chData,
            Orientations = [.. request.Orientations],
            ThicknessesNm = [.. request.ThicknessesNm],
            Sites = [.. request.Sites],
            ModelTier = request.ModelTier,
            IncidentEnergyKeV = request.IncidentEnergyKeV,
            Mu00Nm2 = mu00,
            UnitCellVolumeNm3 = Crystal.Volume,
        };
    }

    private static void Validate(AlchemiRequest r)
    {
        if (r.Orientations is null || r.Orientations.Length == 0) throw new ArgumentException("ALCHEMI: no orientations");
        if (r.ThicknessesNm is null || r.ThicknessesNm.Length == 0) throw new ArgumentException("ALCHEMI: no thicknesses");
        if (r.Sites is null || r.Sites.Length == 0) throw new ArgumentException("ALCHEMI: no site hypotheses");
        if (r.Channels is null || r.Channels.Length == 0) throw new ArgumentException("ALCHEMI: no ionization channels");
        if (r.ThicknessesNm.Any(t => t < 0 || !double.IsFinite(t))) throw new ArgumentException("ALCHEMI: thickness must be finite and non-negative");
        //同一 (Z,Shell) の重複はチャネル表の作り間違いなので hard error (STEM-EDX と同じ規律)
        if (r.Channels.Distinct().Count() != r.Channels.Length) throw new ArgumentException("ALCHEMI: duplicate ionization channel");
        if (r.BaseRotation is null) throw new ArgumentException("ALCHEMI: BaseRotation is null");
    }

    /// <summary>
    /// 260807Cl 追加: FixedUnion 基底 (設計 §5.4)。
    /// **走査の全方位で Find_gVectors を呼んで hkl の真の union を取る** — 「基準方位で 1 回」は union ではなく、
    /// 落ちた g の誤差は方位相関を持つので回帰がサイト固有 ICP と誤認する (ランダムノイズより危険)。
    /// BaseRotation は固定なので gCache は共有され、2 方位目以降はエワルド球スクリーニングだけで済む。
    /// </summary>
    private (Beam[] Union, AlchemiBasisDiagnostic Diagnostic) BuildFixedUnion(
        AlchemiRequest request, AlchemiOrientation[] orientations, Vector3DBase surface,
        double kvac, double u0, int maxNumOfBloch)
    {
        var map = new Dictionary<(int H, int K, int L), Beam>(maxNumOfBloch * 2);
        int centerOnly = 0;
        var centerIndex = orientations.Length / 2;
        for (int i = 0; i < orientations.Length; i++)
        {
            var beams = Find_gVectors(BaseRotation, getVecK0(kvac, u0, orientations[i].BeamDirection, surface), surface, maxNumOfBloch);
            if (i == centerIndex) centerOnly = beams.Length;
            foreach (var b in beams)
                map.TryAdd(b.Index, b);//Vec / Ureal / Uimag は方位非依存なので最初の 1 本で足りる (P,Q は方位ごとに reset)
        }
        if (!map.ContainsKey((0, 0, 0)))
            throw new InvalidOperationException("ALCHEMI: the 000 beam is missing from the union basis — the incident boundary condition psi_g(0) = delta_g0 cannot be imposed");

        //決定的な順序: |g|² → h → k → l。000 は |g|=0 なので必ず先頭に来る (psi0 = e0 の前提)
        var union = map.Values
            .OrderBy(b => b.Vec.Length2).ThenBy(b => b.Index.H).ThenBy(b => b.Index.K).ThenBy(b => b.Index.L)
            .ToArray();

        //走査内 min|s| (union が無駄に太っていないかの指標)。P ≤ 0 はここで捕まえる
        var minS = new double[union.Length];
        Array.Fill(minS, double.PositiveInfinity);
        bool nonPositiveP = false;
        foreach (var o in orientations)
        {
            var vecK0 = getVecK0(kvac, u0, o.BeamDirection, surface);
            for (int g = 0; g < union.Length; g++)
            {
                var (q, p) = getQP(union[g].Vec, vecK0, surface);
                if (p <= 0) { nonPositiveP = true; continue; }
                var s = Math.Abs(Math.Sqrt(p * p / 4 + q) - p / 2);
                if (s < minS[g]) minS[g] = s;
            }
        }

        var tiltSpan = orientations.Max(o => o.TiltRad) - orientations.Min(o => o.TiltRad);
        var warnings = new List<string>();
        //設計 §5.4 の初期目安 (警告閾値であって最終判定ではない — 最終判定は expanded-basis 収束誤差)
        if (tiltSpan > 30e-3)
            warnings.Add($"tilt span {tiltSpan * 1e3:f1} mrad exceeds 30 mrad — outside the v1 FixedUnion guarantee (use TiledUnion, v1.1)");
        else if (tiltSpan > 10e-3)
            warnings.Add($"tilt span {tiltSpan * 1e3:f1} mrad is in the 10-30 mrad range — a true union plus the expanded-basis check is mandatory here");
        if (union.Length > 4 * Math.Max(centerOnly, 1))
            warnings.Add($"the union basis is {(double)union.Length / Math.Max(centerOnly, 1):f1}x the centre-only basis ({centerOnly} -> {union.Length} beams) — the scan range may be too wide for a fixed basis");
        if (nonPositiveP)
            warnings.Add("some union beams have P <= 0 (back-directed) at part of the scan and were skipped in the excitation-error diagnostic");

        var diagnostic = new AlchemiBasisDiagnostic(
            union.Length, centerOnly, tiltSpan,
            minS.Where(double.IsFinite).DefaultIfEmpty(double.NaN).Max(),
            AlchemiBasisDiagnostic.Hash(union.Select(b => b.Index)),
            double.NaN, true, [.. warnings]);
        return (union, diagnostic);
    }

    /// <summary>260807Cl 追加: 方位ループ本体。μ は方位非依存なのでここに入る前に 1 回だけ組む。</summary>
    private (double[] Dyn, double[] Dech, double[] Lcoh, double[] Mu00) SolveOrientations(
        AlchemiRequest request, IonizationData[] chData, Beam[] union, AlchemiOrientation[] orientations,
        AlchemiTensorShape shape, Vector3DBase surface, double kvac, double u0,
        Action<AlchemiProgress> progress, CancellationToken cancel)
    {
        int bLen = union.Length, tLen = request.ThicknessesNm.Length;
        int nSite = request.Sites.Length, nCh = request.Channels.Length;

        //--- 3. μ (サイト × チャネル) ---
        //μ は結晶固定量 = 方位に依らない。FixedUnion なら run 全体でここ 1 回きり (設計 §3.2 の batch 設計)
        progress?.Invoke(new AlchemiProgress(AlchemiStage.BuildingMuMatrices, 0));
        var muBuilder = new AlchemiMuBuilder(Crystal, [.. union.Select(b => b.Index)]);
        var mu = new Complex[nSite * nCh][];
        var mu00 = new double[nSite * nCh];
        for (int s = 0; s < nSite; s++)
            for (int c = 0; c < nCh; c++)
            {
                mu[s * nCh + c] = muBuilder.Build(chData[c], request.Sites[s], request.ModelTier);
                mu00[s * nCh + c] = muBuilder.Mu00(chData[c], request.Sites[s]);
                progress?.Invoke(new AlchemiProgress(AlchemiStage.BuildingMuMatrices, (double)(s * nCh + c + 1) / (nSite * nCh)));
            }
        cancel.ThrowIfCancellationRequested();

        //--- 4. 方位ループ ---
        uDictionary.Clear();
        var potentialMatrix = getPotentialMatrix(union);
        var psi0 = new Complex[bLen];
        psi0[0] = Complex.One;//境界条件 psi_g(0) = delta_g0 (union は 000 が先頭)

        var dyn = new double[shape.Length];
        var dech = new double[shape.Length];
        var lcoh = new double[shape.OrientationCount * tLen];
        int done = 0;
        var options = new ParallelOptions { CancellationToken = cancel };
        if (request.MaxDegreeOfParallelism > 0) options.MaxDegreeOfParallelism = request.MaxDegreeOfParallelism;

        Parallel.For(0, orientations.Length, options, oi =>
        {
            var o = orientations[oi];
            var vecK0 = getVecK0(kvac, u0, o.BeamDirection, surface);
            var eigenMatrix = Shared.Rent(bLen * bLen);
            var beams = System.Buffers.ArrayPool<Beam>.Shared.Rent(bLen);
            try
            {
                reset_gVectors(bLen, union, BaseRotation, vecK0, surface, ref beams);
                getEigenMatrix(bLen, beams, ref eigenMatrix, potentialMatrix);

                Complex[] eigenValues, eigenVectors, alpha;
                if (request.UseNativeSolver && EigenEnabled)
                {
                    (eigenValues, eigenVectors) = NativeWrapper.EigenSolver(bLen, eigenMatrix.AsSpan()[..(bLen * bLen)].ToArray());
                    alpha = NativeWrapper.PartialPivLuSolve(bLen, eigenVectors, psi0);
                }
                else
                {
                    var evd = new DMat(bLen, bLen, eigenMatrix.AsSpan()[..(bLen * bLen)].ToArray()).Evd(Symmetricity.Asymmetric);
                    eigenValues = ((DVec)evd.EigenValues).Values;
                    eigenVectors = ((DMat)evd.EigenVectors).Values;
                    alpha = ((DVec)((DMat)evd.EigenVectors).LU().Solve(new DVec((Complex[])psi0.Clone()))).Values;
                }

                var reduction = new AlchemiReduction(bLen, eigenValues, eigenVectors, alpha, request.ThicknessesNm);
                var l = reduction.CoherentPathLengthNm();//方位あたり 1 回 (以後キャッシュ)
                for (int t = 0; t < tLen; t++) lcoh[oi * tLen + t] = l[t];

                for (int s = 0; s < nSite; s++)
                    for (int c = 0; c < nCh; c++)
                    {
                        var y = reduction.Yield(mu[s * nCh + c], mu00[s * nCh + c], Crystal.Volume, request.IncludeDechannelledComponent);
                        for (int t = 0; t < tLen; t++)
                        {
                            var i = shape.Index(oi, t, s, c);
                            dyn[i] = y.Dynamic[t];
                            dech[i] = y.Dechannelled[t];
                        }
                    }
            }
            finally
            {
                Shared.Return(eigenMatrix);
                System.Buffers.ArrayPool<Beam>.Shared.Return(beams);
                progress?.Invoke(new AlchemiProgress(AlchemiStage.SolvingOrientations,
                    (double)Interlocked.Increment(ref done) / orientations.Length));
            }
        });
        return (dyn, dech, lcoh, mu00);
    }

    /// <summary>260807Cl 追加: expanded-basis 診断 (設計 §5.4)。
    /// 中心・両端の 3 方位を <see cref="AlchemiRequest.ExpandedBasisFactor"/> 倍の基底で解き直し、
    /// Total yield の最大相対差を返す。固定基底からの g 欠落は方位相関を持つので、
    /// この差が閾値を超えたら fit 不適格にする (曲線を混ぜて隠さない)。</summary>
    private double CheckExpandedBasis(AlchemiRequest request, IonizationData[] chData, Vector3DBase surface,
        double kvac, double u0, double[] dyn, double[] dech, AlchemiTensorShape shape, CancellationToken cancel)
    {
        var pick = new[] { 0, request.Orientations.Length / 2, request.Orientations.Length - 1 }.Distinct().ToArray();
        var subset = pick.Select((oi, i) => request.Orientations[oi] with { Index = i }).ToArray();
        var expandedMax = (int)Math.Round(request.MaxNumOfBloch * request.ExpandedBasisFactor);

        //⚠ union は「走査全体」で作る。部分集合だけで作ると基底が縮んで比較にならない
        var (expandedUnion, _) = BuildFixedUnion(request, request.Orientations, surface, kvac, u0, expandedMax);
        var subShape = new AlchemiTensorShape(subset.Length, shape.ThicknessCount, shape.SiteCount, shape.ChannelCount);
        var (exDyn, exDech, _, _) = SolveOrientations(request, chData, expandedUnion, subset, subShape,
            surface, kvac, u0, null, cancel);

        double worst = 0, scale = 0;
        for (int i = 0; i < pick.Length; i++)
            for (int t = 0; t < shape.ThicknessCount; t++)
                for (int s = 0; s < shape.SiteCount; s++)
                    for (int c = 0; c < shape.ChannelCount; c++)
                    {
                        var a = dyn[shape.Index(pick[i], t, s, c)] + dech[shape.Index(pick[i], t, s, c)];
                        var b = exDyn[subShape.Index(i, t, s, c)] + exDech[subShape.Index(i, t, s, c)];
                        worst = Math.Max(worst, Math.Abs(a - b));
                        scale = Math.Max(scale, Math.Abs(b));
                    }
        return scale > 0 ? worst / scale : worst;
    }

    #endregion
}
