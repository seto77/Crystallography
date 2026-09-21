#region using
using System;
using System.Threading.Tasks;
#endregion

namespace Crystallography;

/// <summary>
/// ラスター (width×height) のピクセル中心を表示パターン座標 (mm、検出器中心基準) へ写す係数と、検出器の物理サイズ。260726Cl 追加。
/// px_view = (2w+1-width)·ScaleW + OffX、py_view = (2h+1-height)·ScaleH + OffY (OffX/OffY = 表示のパン量)。
/// HalfWidth/HalfHeight は検出器の物理半幅・半高 (mm)。
/// ⚠ 260921Cl: MC のビニングを検出器から射出半球へ移したので、ビンの内挿にはもう使わない (旧: MC ビン補間の正規化と、検出器外の端ビン外挿)。
///   現在は A(E) の Ā を検出器面で (画素の立体角で重み付けして) 平均する範囲として使う (EbsdPatternComposer.MeanCoherentFraction)。
///   表示の視野 (パン・ズーム) ではなく検出器の寸法なので、Ā は視野に依存しない。
/// XMirror は左右反転トグル (±1)、DetX は検出器中心の X オフセット (mm)。
/// </summary>
public readonly record struct EbsdRasterView(
    double ScaleW, double ScaleH, double OffX, double OffY,
    double HalfWidth, double HalfHeight, double XMirror, double DetX);

/// <summary>
/// MasterPattern から表示用 EBSD パターン (現在の視野全体をカバーするラスター) を合成する。
/// 260726Cl 追加: FormEBSD.cs にあった「ピクセルごとの MasterPattern 参照テーブル構築 + 3 つの合成モデル」を
/// そのまま移設したもの (GUI 非依存の純計算)。UI から読む値 (検出器幾何・視野・トラックバー) は引数で受け取る。
///
/// 2 段構成:
///   ① <see cref="BuildLookupTable"/> — 視線方向 → Rosca-Lambert → 格子補間係数を 1 回だけ作る (エネルギー・深さに依存しない)
///   ② ApplySingleSliceModelN / ApplyWeightedModelN — 構築済みテーブルを使って強度を書き込む
/// モデルの意味は masterPatternCombinationModel と同じ: 0=current、1=globally normalized master、2=absolute MC × differential master。
///
/// <see cref="EbsdPatternProjector"/> との違い: あちらは実測との照合用 (検出器 native グリッド・単一スライス・方位を毎回変える)。
/// こちらは表示用 (視野全体のラスター・MC 重み付き全スライス合成・方位は結晶の現在値)。
/// </summary>
public sealed class EbsdPatternComposer
{
    //260727Cl (/simplify): Lambert 逆写像を SphereToRoscaLambertSquare へ委譲したため未参照になった定数を削除
    //旧: const double Inv_PI = 1 / Math.PI; / const double Half_PI = 0.5 * Math.PI;

    // 260325Cl: ピクセルごとの MasterPattern 参照テーブル (エネルギー・深さに依存しない)
    // 正方格子: idx[i] = 左上グリッドインデックス, wt[i*2] = fw, wt[i*2+1] = fh, posZ[i] = 半球
    // 六方格子: idx[i*3..i*3+2] = 3近傍インデックス, wt[i*3..i*3+2] = バリセントリック重み (260331Cl)
    int[] lookupIdx = [];
    float[] lookupWt = [];
    bool[] lookupPosZ = [];
    int lookupGridSize; // 260325Cl: Apply で idx+gridSize の復元に使用
    MasterPattern.Types lookupGridType; // 260331Cl: 六方格子モードかどうか

    // 260325Cl: DetTilt/SmpTilt 由来の回転係数キャッシュ (tilt 変更時のみ再計算)
    double yCoeffPy, zCoeffPy, yConst, zConst;

    /// <summary>DetTilt/SmpTilt/DetY/DetZ から回転係数を再計算する。260325Cl 追加 (260726Cl: FormEBSD.UpdateEbsdTiltCoeffs から移設)</summary>
    public void UpdateTiltCoefficients(double detTilt, double sampleTilt, double detY, double detZ)
    {
        var (sinDet, cosDet) = Math.SinCos(detTilt);
        var (sinSmp, cosSmp) = Math.SinCos(sampleTilt);
        yCoeffPy = cosSmp * cosDet + sinSmp * sinDet;
        zCoeffPy = -sinSmp * cosDet + cosSmp * sinDet;
        yConst = cosSmp * detY + sinSmp * detZ;
        zConst = -sinSmp * detY + cosSmp * detZ;
    }

    /// <summary>
    /// 検出器ジオメトリと結晶方位から、ピクセルごとの MasterPattern 参照テーブルを構築する。260325Cl 追加
    /// エネルギー・深さに依存しないため、畳み込み時は 1 回だけ呼べばよい。
    /// </summary>
    //260726Cl シグネチャ変更 (FormEBSD から移設): UI 直読だった MasterPattern・Crystal.RotationMatrix・検出器/視野の値を引数化。
    //旧: private unsafe void BuildEbsdLookupTable(int width, int height)
    public unsafe void BuildLookupTable(MasterPattern mp, Matrix3D rotation, int width, int height, in EbsdRasterView view)
    {
        ArgumentNullException.ThrowIfNull(mp);
        if (width <= 0) throw new ArgumentOutOfRangeException(nameof(width), "EBSD raster width must be positive."); //260725Ch: unsafe 配列長と除算の前提を入口で保証
        if (height <= 0) throw new ArgumentOutOfRangeException(nameof(height), "EBSD raster height must be positive."); //260725Ch
        var gridSize = mp.GridSize; //260725Ch: 1×1 以下では正方格子の gridMax-1 と bilinear の idx+1 が成立しない
        if (gridSize < 2) throw new InvalidOperationException("MasterPattern.GridSize must be at least 2."); //260725Ch
        //var totalPixels = width * height; //260725Ch 変更前
        var totalPixels = checked(width * height); //260725Ch: 将来ラスター上限が変わっても unsafe 配列長の整数オーバーフローを許さない
        lookupGridType = mp.GridType; // 260331Cl
        var isHexGrid = lookupGridType == MasterPattern.Types.Hexagonal; // 260331Cl

        // 260331Cl: 六方格子は 3 idx + 3 wt/pixel、正方格子は 1 idx + 2 wt/pixel
        //260727Cl: 上の checked(width*height) と対称にする (旧 unchecked。3 倍・2 倍でオーバーフローすると
        //  idxCount が小さな正数になり、再確保がスキップされて unsafe 書き込みが確保長を超え得た)
        var idxCount = isHexGrid ? checked(totalPixels * 3) : totalPixels;
        var wtCount = checked(totalPixels * (isHexGrid ? 3 : 2));
        // if (lookupIdx.Length != idxCount) // 260724Cl 変更前: idx 長のみで判定。六方⇔正方切替とラスター画素数変化の組合せで wt/posZ 長だけが不整合になり、unsafe ループが境界外へ書く恐れがあった
        if (lookupIdx.Length != idxCount || lookupWt.Length != wtCount || lookupPosZ.Length != totalPixels) // 260724Cl
        {
            //lookupIdx = new int[idxCount]; // 260402Cl 変更前
            //lookupWt = new float[wtCount]; // 260402Cl 変更前
            lookupIdx = GC.AllocateUninitializedArray<int>(idxCount); // 260402Cl 直後に全要素上書きされるため未初期化で確保
            lookupWt = GC.AllocateUninitializedArray<float>(wtCount); // 260402Cl
            lookupPosZ = GC.AllocateUninitializedArray<bool>(totalPixels); // 260402Cl
        }

        // 260325Cl: tilt 係数はキャッシュ済み (UpdateTiltCoefficients で更新)
        double yCoeffPyLocal = yCoeffPy, zCoeffPyLocal = zCoeffPy, yConstLocal = yConst, zConstLocal = zConst;

        var Ri = rotation.Inverse();
        double xm = view.XMirror; // 260718Cl: 左右反転トグル (既定 +1)。X 方向のみ符号を掛ける
        double ax = -xm * Ri.E11, ay = -xm * Ri.E21, az = -xm * Ri.E31, bx = Ri.E12 * yCoeffPyLocal + Ri.E13 * zCoeffPyLocal;
        double by = Ri.E22 * yCoeffPyLocal + Ri.E23 * zCoeffPyLocal, bz = Ri.E32 * yCoeffPyLocal + Ri.E33 * zCoeffPyLocal;
        // 260723Cl 変更: 検出器中心 X オフセット (DetX) を定数項に追加。視線ベクトルの X 成分は (ピクセル項) - DetX
        // (既存の Y/Z 定数項 +yConst/+zConst が検出器中心 -C 由来なのと同じ規約。lab X は試料/検出器傾斜 (X 軸回転) で不変なので Ri の第 1 列に掛かる)
        // double cx = Ri.E12 * yConst + Ri.E13 * zConst, cy = ..., cz = ...; // 260723Cl 変更前
        double cx = Ri.E12 * yConstLocal + Ri.E13 * zConstLocal - Ri.E11 * view.DetX;
        double cy = Ri.E22 * yConstLocal + Ri.E23 * zConstLocal - Ri.E21 * view.DetX;
        double cz = Ri.E32 * yConstLocal + Ri.E33 * zConstLocal - Ri.E31 * view.DetX;

        // double scaleW = DetR / width, scaleH = DetR / height; // 260723Cl 変更前: 画面幅=検出器直径 (2·DetR) 前提
        // double scaleW = DetHalfWidth / width, scaleH = DetHalfHeight / height; // 260723Cl 変更: ラスター全域=検出器の物理サイズ // 260724Cl 変更前
        // 260724Cl 変更: ラスター全域=現在の視野 (ClientSize×Resolution)。視野中心のずれ (viewPan) は定数項へ ax/bx 系数経由で加算
        double scaleW = view.ScaleW, scaleH = view.ScaleH;
        cx += ax * view.OffX + bx * view.OffY;
        cy += ay * view.OffX + by * view.OffY;
        cz += az * view.OffX + bz * view.OffY;
        double ax2 = ax * scaleW, ay2 = ay * scaleW, az2 = az * scaleW;
        double bx2 = bx * scaleH, by2 = by * scaleH, bz2 = bz * scaleH;

        lookupGridSize = gridSize; // 260325Cl: Apply 用にキャッシュ
        var startPxFactor = 1 - width; // (260325Ch) 列方向は等間隔なので、旧 pxFactor 再計算を増分更新へ置き換える
        double dxStep = 2.0 * ax2, dyStep = 2.0 * ay2, dzStep = 2.0 * az2; // (260325Ch)

        fixed (int* pIdx = lookupIdx)
        fixed (float* pWt = lookupWt)
        fixed (bool* pPosZ = lookupPosZ)
        {
            var pIdx0 = pIdx; var pWt0 = pWt; var pPosZ0 = pPosZ;

            if (isHexGrid) // 六方格子パス
            {
                var gs = gridSize; // Parallel.For のキャプチャ用ローカル
                Parallel.For(0, height, h =>
                {
                    double pyFactor = 2 * h + 1 - height;
                    double rowBx = bx2 * pyFactor + cx, rowBy = by2 * pyFactor + cy, rowBz = bz2 * pyFactor + cz;
                    int rowOffset = h * width;
                    double dx = ax2 * startPxFactor + rowBx, dy = ay2 * startPxFactor + rowBy, dz = az2 * startPxFactor + rowBz;

                    for (int w = 0; w < width; w++)
                    {
                        int i = rowOffset + w;
                        double invLen = 1.0 / Math.Sqrt(dx * dx + dy * dy + dz * dz);
                        double nx = dx * invLen, ny = dy * invLen, nz = dz * invLen;
                        pPosZ0[i] = nz >= 0;

                        // 球面→六方格子
                        var (hx, hy) = MasterPattern.SphereToRoscaLambertHex(nx, ny, Math.Abs(nz));
                        MasterPattern.GetHexBarycentricLookup(hx, hy, gs, out int idx0, out int idx1, out int idx2, out float bw0, out float bw1, out float bw2);

                        int i3 = i * 3;
                        pIdx0[i3] = idx0; pIdx0[i3 + 1] = idx1; pIdx0[i3 + 2] = idx2;
                        pWt0[i3] = bw0; pWt0[i3 + 1] = bw1; pWt0[i3 + 2] = bw2;

                        dx += dxStep; dy += dyStep; dz += dzStep;
                    }
                });
            }
            else // 正方格子パス
            {
                var sqLim = MasterPattern.SquareLimit;
                var invStep = gridSize / (2.0 * sqLim);
                var gridMax = gridSize - 1;

                Parallel.For(0, height, h =>
                {
                    double pyFactor = 2 * h + 1 - height;
                    double rowBx = bx2 * pyFactor + cx, rowBy = by2 * pyFactor + cy, rowBz = bz2 * pyFactor + cz;
                    int rowOffset = h * width;
                    double dx = ax2 * startPxFactor + rowBx, dy = ay2 * startPxFactor + rowBy, dz = az2 * startPxFactor + rowBz;

                    for (int w = 0; w < width; w++)
                    {
                        int i = rowOffset + w;
                        pPosZ0[i] = dz >= 0;

                        //260727Cl (/simplify): インライン展開していた Lambert 逆写像を、六方パス (SphereToRoscaLambertHex) と
                        //  同方針で共有ヘルパへ委譲した。EbsdPatternProjector.Project が 260724Cl に済ませた整理と同じもので、
                        //  数式は同一 (旧 edgeRadius=√(π/2·(1−|z|)) は SphereToRoscaLambertSquare の radialScale 経由の式と等価。
                        //  正規化もヘルパ内で行うので invLen は不要)。400,000 方向 × grid 256/512/1024 で照合し、
                        //  グリッド索引は完全一致・補間重みの差は最大 5.1e-07 (float 精度の約 9 倍で表示のみに効く桁)。
                        //旧: double absDx=|dx|, absDy=|dy|; invLen=1/√(dx²+dy²+dz²); absDzNorm=|dz|·invLen;
                        //    if (absDx<1e-15 && absDy<1e-15) (a,b)=(0,0);
                        //    else { edgeRadius=√(max(0, Half_PI·(1−absDzNorm)));
                        //           if (absDx>=absDy) { a=±edgeRadius(dx の符号); b=4a·Inv_PI·atan(dy/dx); }
                        //           else              { b=±edgeRadius(dy の符号); a=4b·Inv_PI·atan(dx/dy); } }
                        var (a, b) = MasterPattern.SphereToRoscaLambertSquare(dx, dy, Math.Abs(dz));

                        double gw = (a + sqLim) * invStep - 0.5, gh = (sqLim - b) * invStep - 0.5;
                        int w0 = (int)Math.Floor(gw), h0 = (int)Math.Floor(gh);
                        double fw = gw - w0, fh = gh - h0;
                        if (w0 < 0) { w0 = 0; fw = 0; } else if (w0 >= gridMax) { w0 = gridMax - 1; fw = 1; }
                        if (h0 < 0) { h0 = 0; fh = 0; } else if (h0 >= gridMax) { h0 = gridMax - 1; fh = 1; }

                        pIdx0[i] = h0 * gridSize + w0;
                        int i2 = i * 2;
                        pWt0[i2] = (float)fw;
                        pWt0[i2 + 1] = (float)fh;

                        dx += dxStep; dy += dyStep; dz += dzStep;
                    }
                });
            }
        }
    }

    #region masterPatternCombinationModel 0 — current
    /// <summary>構築済みルックアップテーブルを使い、指定 energy/depth の EBSD パターンを values に書き込む。260325Cl 追加</summary>
    //260727Cl (/simplify): 本文は ApplySingleSliceModel1 の planeScaleFactor=1.0 と 1 文字違い (scaleFactor 乗算の有無) の
    //  完全な写しだったので委譲に置き換えた。float 式の結果に 1.0 (double) を掛けても値は不変なのでビット一致。
    //  旧実装は約 45 行 (欠損スライスガード・hex 3 近傍バリセントリック・square バイリニアの 2 パス) — ApplySingleSliceModel1 を参照。
    public unsafe void ApplySingleSliceModel0(double[] values, int totalPixels, float[] posPlane, float[] negPlane) // 260325Cl: unsafe 化
        => ApplySingleSliceModel1(values, totalPixels, posPlane, negPlane);

    /// <summary>260727Cl 追加: weighted 合成の前提 — 分布の重み配列が mp と同じ (energy × depth) 格子であること — を検査する。
    /// 重み参照は wIdx = ei * dLen + di (dLen は mp 由来) を dist の配列へ素通しするので、格子がずれると
    /// dLen が減る方向では**例外も警告も出さずに別スライスの重み**を引き、増える方向では IndexOutOfRange になる。
    /// MC は MasterPattern の構築とは別のタイミングでも走る (Calc BSE・深さ格子や非晶質層を変えたときの再ビニング。
    /// 260921Cl: 射出半球のビニングになって検出器幾何の変更では再ビニングしなくなった) ので、
    /// 呼び出し側の規律だけに任せず入口で弾く。</summary>
    static void EnsureGridMatches(MasterPattern mp, EbsdMonteCarloDistribution dist)
    {
        ArgumentNullException.ThrowIfNull(mp);
        ArgumentNullException.ThrowIfNull(dist);
        if (!dist.MatchesGridOf(mp))
            throw new ArgumentException(
                $"The MC distribution grid ({dist.EnergyCount} energies x {dist.DepthCount} depths) does not match the MasterPattern ({mp.Energies.Length} x {mp.Depths.Length}).",
                nameof(dist));
    }

    /// <summary>260718Cl 追加: 全 energy×depth の plane 参照をピクセルループ外で一括取得する (Weighted 3 モデル共通)。
    /// 旧実装はピクセル×energy×depth 回 GetPlane を呼んでいた。</summary>
    static (float[][] pos, float[][] neg) GetAllPlanes(MasterPattern mp, int eLen, int dLen)
    {
        var pos = new float[eLen * dLen][];
        var neg = new float[eLen * dLen][];
        int requiredPlaneLength = checked(mp.GridSize * mp.GridSize); //260725Ch
        for (int ei = 0; ei < eLen; ei++)
            for (int di = 0; di < dLen; di++)
            {
                //pos[ei * dLen + di] = mp.GetPlane(MasterPattern.Hemisphere.PositiveZ, ei, di); //260725Ch 変更前
                //neg[ei * dLen + di] = mp.GetPlane(MasterPattern.Hemisphere.NegativeZ, ei, di);
                int index = ei * dLen + di;
                var plane = mp.GetPlane(MasterPattern.Hemisphere.PositiveZ, ei, di);
                pos[index] = plane != null && plane.Length >= requiredPlaneLength ? plane : null; //260725Ch: weighted unsafe 経路へ短い plane を渡さない
                plane = mp.GetPlane(MasterPattern.Hemisphere.NegativeZ, ei, di);
                neg[index] = plane != null && plane.Length >= requiredPlaneLength ? plane : null;
            }
        return (pos, neg);
    }

    /// <summary>260919Cl 追加: 各 (energy, depth) 平面の方向平均 (全セルの単純平均)。depthWidths を渡すと model 2 と同じ差分 max(0, M_d − M_{d−1})/Δt の平均。
    /// 表面非晶質層に源を持つ電子 (菊池変調なし・電子 1 本あたりの強度は結晶内の源と同じ) の寄与に使う。</summary>
    static (double[] pos, double[] neg) GetPlaneMeans(float[][] posPlanes, float[][] negPlanes, int dLen, double[] depthWidths)
    {
        double[] Means(float[][] planes)
        {
            var m = new double[planes.Length];
            Parallel.For(0, planes.Length, w =>
            {
                var p = planes[w];
                if (p == null || p.Length == 0) return;
                int di = w % dLen;
                var prev = depthWidths != null && di > 0 ? planes[w - 1] : null;
                double invWidth = depthWidths != null ? 1.0 / depthWidths[di] : 1.0;
                // (/simplify2) 六方格子の無効セルは plane に 0 が入るので、累積強度 p[i] > 0 のセルだけを平均に入れる (正方格子では全セル有効)
                double sum = 0; long n = 0;
                if (prev != null && prev.Length >= p.Length)
                    for (int i = 0; i < p.Length; i++) { if (p[i] > 0) { sum += Math.Max(0.0, p[i] - prev[i]); n++; } }
                else if (depthWidths != null)
                    for (int i = 0; i < p.Length; i++) { if (p[i] > 0) { sum += p[i]; n++; } }
                else
                    for (int i = 0; i < p.Length; i++) { if (p[i] > 0) { sum += p[i]; n++; } }
                m[w] = n > 0 ? sum / n * invWidth : 0;
            });
            return m;
        }
        return (Means(posPlanes), Means(negPlanes));
    }

    /// <summary>260919Cl 追加: 各 (energy, depth) の平面平均を、全ビン合算の深さ重み g で d 方向に平均した「エネルギーごとの基準値」(長さ eLen) にする。
    /// 非晶質層内の源の強度は「同じエネルギーの結晶内の源の平均強度」と定義し、ビンごとの深さフィットには依存させない (依存させると電子数の少ないビンで一様成分がムラになる)。</summary>
    static double[] CollapseToEnergyReference(double[] means, int dLen, double[] g)
    {
        int eLen = means.Length / dLen;
        var o = new double[eLen]; // (/simplify2) d 方向に複製せず、エネルギーごとの値だけ返す (消費側は [ei] で引く)
        for (int ei = 0; ei < eLen; ei++)
        {
            double num = 0, den = 0, plain = 0;
            for (int di = 0; di < dLen; di++)
            {
                int w = ei * dLen + di;
                double gw = g != null && w < g.Length ? g[w] : 0;
                num += means[w] * gw; den += gw; plain += means[w];
            }
            o[ei] = den > 0 ? num / den : plain / dLen;
        }
        return o;
    }

    /// <summary>260919Cl 追加: 同じ MasterPattern・同じ MC 分布・同じ差分フラグなら基準値を再計算しない (パン・ズームのたびに 40M 要素を舐めない)。</summary>
    //260921Cl 変更: Sphere (平面ごとの全球方向平均) を追加
    private (MasterPattern Mp, EbsdMonteCarloDistribution Dist, bool Differential, double[] Pos, double[] Neg, double[] Sphere) planeMeanCache;
    /// <summary>260920Cl 追加 (作者指示): ビームエネルギー [keV]。<see cref="CoherenceLossDecayKeV"/> の基準 E0。NaN なら master pattern の最大エネルギーを使う</summary>
    public double BeamEnergyKeV { get; set; } = double.NaN;

    /// <summary>260920Cl 追加 (作者指示): 損失依存のコントラスト係数 A(E) = exp(−(E0 − E)/E_c) の特性損失 E_c [keV]。
    /// NaN または 0 以下で無効 (= 全エネルギースライスが満額のコントラストを持つ従来動作)。
    /// 【なぜ要るか】master pattern のスライスはどれも自分のエネルギーで正しく計算されているが、合成の重みは MC の
    ///   後方散乱電子の**射出エネルギー分布そのもの**である。菊池バンドの角幅は λ(E) に比例するので、射出エネルギーの
    ///   中央値が 0.83·E0 (Si 20 kV 実測) であるぶんだけ合成パターンのバンドが広がる。実測 (Si004) の実効波長は
    ///   ほぼ λ(E0) に対応しており、エネルギーフィルタ EBSD の実験でも鋭い菊池コントラストは低損失電子に集中し、
    ///   大きく損失した電子はほとんど特徴のない背景しか作らないことが知られている。つまり「電子数」で重み付けするのは
    ///   正しいが、電子 1 個あたりの**変調の振幅**が損失とともに落ちる項が抜けていた。
    /// 【何をするか】コヒーレント成分を A(E) 倍し、失った (1 − A(E)) 分を同じスライスの方向平均 (= 平坦な台座) へ回す。
    ///   総量は保存するので明るさは変わらず、バンドのコントラストと幅だけが変わる。表面非晶質層 (BinAmorphousFraction) と
    ///   同じ配分機構で、両者は掛け合わさる。
    ///   260921Cl: 台座を 1 個のスカラーにすると画素ごとの Σ W·A のムラが明るさに残ったので、現行は
    ///   V = Σ W·M̄ + Ā·(Σ W)·(⟨I⟩_A − ⟨M̄⟩_A) (明るさは A(E) 無効時と同じ、コントラストは全面で Ā 倍)。式と Ā の定義は <see cref="MeanCoherentFraction"/> の doc。
    /// 【実測】`tools/EbsdProfileFit` のプロファイル一致法 (正本 §2.7)。実測 Si004 と合成パターンを同じ推定量で測り、
    ///   実効波長比 α = 実測の実効波長 / λ(E0) で比べた (grid 512、既定の 16 エネルギー × 40 深さグリッド):
    ///   実測 α = 1.017 (反射ごと中央値) / 1.019 (全画像 ZNCC)。合成は A(E) 無効で 1.130 / 1.166、
    ///   E_c = 1.0 keV で 1.030 / 1.036、E_c = 0.7 keV で 1.017 / 1.020。→ E_c = 0.71 / 0.68 keV が最良値。
    /// ⚠ 既定 0.7 keV は Si 20 kV の実測パターン 1 枚で校正した経験値であって、物質・加速電圧に依らない定数ではない。
    ///   E_c は「損失がどれだけ溜まると鋭いコントラストが失われるか」を 1 個のスカラーへ畳んだものなので、
    ///   阻止能 (物質の ρ·Z/A と E0) と劣化の経路長 (弾性・非弾性平均自由行程) の積のオーダーで決まるはずで、どちらも物質依存する。
    ///   無次元化の候補は E_c/ħω_p (Si 20 kV では ≈ 42 プラズモン) と E_c ≈ S·Λ。
    /// 【2 物質目 (260921Cl)】Botallackite Cu₂(OH)₃Cl 20 kV でも同じ手順で測った (正本 §2.7.1)。
    ///   合成は A(E) 無効で α̂ = 1.139、E_c = 0.7 keV で 1.018、0.25 keV で 1.001。実測 α = 1.001。
    ///   → E_c ≤ 0.25 keV。系統誤差 (格子定数の絶対スケールが α と縮退する) を最大に見積もっても 0.6 keV 以下で、
    ///   **Si の 0.70 keV より小さい**。プラズモン説 (1.12 keV 予想) と「阻止能 × 一定 Λ」説 (0.99 keV 予想) は
    ///   どちらも符号が逆で棄却方向。整合するのは「原子あたり阻止断面積 ÷ 弾性/TDS 断面積」(0.44〜0.51 keV 予想) だけ。
    ///   ⚠ 2 物質では指数も機構も決まらないので、**既定 0.7 keV は据え置き**、物質ごとに変える実装も入れない。
    ///   次の独立検証は Cu → ダイヤモンド → Au、および同一物質の温度依存 (正本 §6 P0-0 / 文献調査 md §11)。</summary>
    public double CoherenceLossDecayKeV { get; set; } = double.NaN;

    //260921Cl 変更前 (作者報告: A(E) を入れるとビン間隔を周期とするブロードな暗線が出る)。
    //  Smoothstep はビン中心で傾きを 0 にするので、直線的な勾配が「平坦→急変→平坦」の階段になり、
    //  背景平坦化の高域通過がそれをブロードな縞にしていた。BinSplineTaps へ置き換え。
    //static double SmoothBinFraction(double f)
    //{
    //    f = Math.Clamp(f, 0, 1);
    //    return f * f * (3 - 2 * f);
    //}

    //260921Cl (/simplify2): 下は旧 2 次版の doc。/// のままだと現行メソッドの /// と連結されて summary が 2 つになるので //// にした。
    //// <summary>260921Cl 追加 (作者報告: A(E) を入れるとビン間隔を周期とするブロードな暗線が複数出る):
    //// MC の 8×8 検出器ビンを内挿する重み。ビン中心を整数とする座標 b に対し、3 タップの **2 次 B スプライン**を返す。
    //// <para>【なぜ 2 タップではだめか】ビンの値は「滑らかな場を 8×8 で粗くサンプルしたもの」なので、内挿には
    //// (1) 継ぎ目で折れない (C1) ことと、(2) 直線的な勾配をそのまま再現することの**両方**が要る。2 タップでは両立しない:</para>
    //// <para>・素の双線形 … 直線は再現するが C0。継ぎ目で傾きが折れ、背景平坦化 (高域通過) が**細い暗線**にする
    ////   (260920Cl に作者が報告した症状。790x602 の実例で列 49/740・行 37/564 = 最外の継ぎ目に一致)。</para>
    //// <para>・Smoothstep 3f²−2f³ … C1 だがビン中心で傾きが 0 になるので、直線の勾配が「平坦→急変→平坦」の階段になる。
    ////   高域通過がこれを**ビン間隔 (視野幅/8) を周期とするブロードな縞**にする (260921Cl に作者が報告した症状)。</para>
    //// <para>2 次 B スプライン q(t) = [(0.5−t)²/2, 0.75−t², (0.5+t)²/2] (t = b − round(b) ∈ [−0.5, 0.5]) なら
    //// C1 かつ 1 次を厳密に再現する (Σq = 1、Σq·(中心位置) = b) ので、どちらの縞も原理的に出ない。</para>
    //// <para>⚠ これは内挿ではなく**近似**で、ビン中心の値は 0.75·v_i + 0.125·(v_{i−1} + v_{i+1}) になる。
    //// つまりビンごとの当てはめのばらつきを 1 ビンぶん均す。元の場 (取り出し角による射出分布の変化) は滑らかなはずで、
    //// 8×8 の当てはめノイズの方が偽物なので、これは副作用ではなく望ましい性質。</para>
    //// <para>⚠ なぜ縞が「A(E) を入れたときだけ」見えるか: A(E) は菊池コントラストを Ā ≈ 0.07 倍に潰す
    //// (Forsterite 20 kV・E_c = 0.7 keV の実測で、背景平坦化後の rms が 11.6 倍小さくなる)。
    //// 重み場由来の縞は A(E) の有無にかかわらず同じ振幅で存在しているが、表示が min/max で自動伸張されるため、
    //// コントラストが潰れたぶんだけ縞が相対的に前に出る。**縞は A(E) が作るのではなく、A(E) が露出させる。**</para>
    //// <para>検出器の外は端のビンへクランプ (従来と同じ外挿規約)。</para></summary>
    //260921Cl 変更前 (2 次 = C1)。C1 だと Laplacian が区画ごとに一定になり、背景平坦化がそれを 8x8 のブロックとして映す。
    //static void BinSplineTaps(double b, int binCount, out int i0, out int i1, out int i2, out double w0, out double w1, out double w2)
    //{
    //    int i = Math.Clamp((int)Math.Floor(b + 0.5), 0, binCount - 1);
    //    double t = Math.Clamp(b - i, -0.5, 0.5);
    //    double lo = 0.5 - t, hi = 0.5 + t;
    //    w0 = 0.5 * lo * lo; w1 = 0.75 - t * t; w2 = 0.5 * hi * hi;
    //    i0 = i > 0 ? i - 1 : 0; i1 = i; i2 = i < binCount - 1 ? i + 1 : binCount - 1;
    //}

    /// <summary>260921Cl: ビン中心を整数とする座標 b に対する 4 タップの **3 次 B スプライン** 重み。
    /// <para>【要件はなぜ C2 なのか】表示の背景平坦化は「原画像 − Gaussian ぼかし」で、ぼかし幅 σ が
    /// ビン間隔よりずっと小さいとき 原画像 − ぼかし ≈ (σ²/2)·∇²(原画像) になる。つまり画面に出るのは
    /// **重み場の 2 階微分**である。したがって内挿に必要なのは</para>
    /// <para>(1) 1 次を厳密に再現すること (さもないと直線勾配が階段になる)、
    /// (2) <b>2 階微分が連続 = C2</b> であること (さもないと 2 階微分の跳びがそのまま模様になる)。</para>
    /// <para>実際に観測された症状はこの順で説明がつく:</para>
    /// <para>・素の双線形 (C0) … ∇² が節点で δ 関数 → **継ぎ目の細い暗線** (260920Cl の症状)。</para>
    /// <para>・Smoothstep (C1 だが 1 次を再現しない) … 直線勾配が階段になる → **ビン間隔周期のブロードな縞** (260921Cl の症状)。</para>
    /// <para>・2 次 B スプライン (C1 + 1 次再現) … 縞も細線も消えるが、∇² が区画ごとに一定なので
    ///   **8×8 のブロック**が残る (背景成分だけを取り出すと明瞭に見える)。</para>
    /// <para>・3 次 B スプライン (C2 + 1 次再現) … ∇² が連続になるのでブロックの縁も消える。これが現行。</para>
    /// <para>基底は t = b − floor(b) ∈ [0,1) に対し
    /// [(1−t)³, 3t³−6t²+4, −3t³+3t²+3t+1, t³]/6、タップ位置は floor(b)−1 … floor(b)+2。
    /// Σw = 1、Σw·(タップ位置) = b (1 次を厳密再現)。</para>
    /// <para>⚠ 内挿ではなく**近似**で、ビン中心の値は (v_{i−1} + 4v_i + v_{i+1})/6 になる。ビンごとに独立な当てはめの
    /// 統計ノイズを 1 ビンぶん均す効果があり、元の場 (取り出し角による射出分布の変化) は滑らかなはずなので望ましい。</para>
    /// <para>⚠ なぜ縞が「A(E) を入れたときだけ」見えるか: A(E) は菊池コントラストを Ā ≈ 0.07 倍に潰す
    /// (Forsterite 20 kV・E_c = 0.7 keV の実測で、背景平坦化後の rms が 11.6 倍小さくなる)。
    /// 重み場由来の模様は A(E) の有無にかかわらず同じ振幅で存在するが、表示が min/max で自動伸張されるため、
    /// コントラストが潰れたぶんだけ前に出る。**A(E) は模様を作るのではなく、露出させる。**</para>
    /// <para>【端の扱い】添字だけをクランプする (= 制御点を端で複製する、B スプライン標準の境界条件)。重みは標準基底のままなので
    /// <b>必ず非負</b>、しかも複製した制御点の上でもやはり B スプラインなので <b>C2 のまま</b>。代償は端の 1 ビンで場が平らに寄ること。
    /// 射出半球をビニングする現行 (260921Cl) では、試料表面より上の方向は b ∈ [−0.5, binCount−0.5] に収まるので、
    /// この処理が効くのは地平線の近くだけ (地平線より下の画素は合成側で 0 にしている)。
    /// 260921Cl (Lambert 等積ディスク): 円板から外れたビンは分布側で内側から延長してある (EbsdMonteCarloDistribution の ctor の doc【縁のビン】) ので、
    /// 4×4 タップが円板の外へ届いても 0 を拾わない。
    /// b の ±3 丸めは極端な値 (視野を極端に広く取ったとき) の int 変換対策。⚠ NaN は Math.Clamp を素通りするので保険にはならない
    /// (260921Cl /simplify2 で訂正。合成側は !(oz &gt; 0) の行を先に落とすので NaN は届かない)。</para>
    /// <para>⚠ 検出器ビニングだった頃の教訓 (再び検出器を切るときのために残す): <b>b 自体をクランプしてはいけない</b>
    /// (画面座標での傾きがそこで折れて枠状の線になった)。<b>制御点の 1 次外挿もいけない</b> (端で重みが負になり、
    /// 合成側の weight ≤ 0 読み飛ばしと合わさって検出器の外の明るさが 42 % 暗くなった)。</para></summary>
    static void BinSplineTaps(double b, int binCount, out int i0, out int i1, out int i2, out int i3,
        out double w0, out double w1, out double w2, out double w3)
    {
        int last = binCount - 1;
        //⚠ b 自体はクランプしない。b を止めると画面座標で傾きが折れ、高域通過がそこを線にする
        //  (260921Cl に一度これで枠状の線を作った)。添字だけをクランプする = 制御点を端で複製する、
        //  という B スプライン標準の境界条件を使う。b ≤ −2 / b ≥ last+2 では 4 タップとも端のビンへ
        //  落ちて値が凍るが、その手前まで C2 で滑らかに漸近するので折れ目が出ない。
        //  評価範囲だけ有限へ丸めておく (視野を極端に広く取ったときの int 変換対策)。
        double bc = Math.Clamp(b, -3.0, last + 3.0);
        int i = (int)Math.Floor(bc);
        double t = bc - i, t2 = t * t, t3 = t2 * t, u = 1 - t;
        const double inv6 = 1.0 / 6.0;
        w0 = inv6 * u * u * u;
        w1 = inv6 * (3 * t3 - 6 * t2 + 4);
        w2 = inv6 * (-3 * t3 + 3 * t2 + 3 * t + 1);
        w3 = inv6 * t3;
        i0 = Math.Clamp(i - 1, 0, last); i1 = Math.Clamp(i, 0, last);
        i2 = Math.Clamp(i + 1, 0, last); i3 = Math.Clamp(i + 2, 0, last);
    }

    /// <summary>260921Cl 追加 (/simplify): 画素 → 試料系の<b>射出方向</b>。3 つの合成モデルで共有する (符号規約はここ 1 か所)。
    /// <para>⚠ <see cref="BuildLookupTable"/> が作る画素方向 (結晶回転を掛ける前) は射出方向の<b>逆向き</b>なので、ここで符号を反転している。
    /// 合成器の向きを +lab へ戻すと検出器交点が求まらず、−lab にすると旧・検出器ビニングの (px, py) が厳密に一致することを
    /// 両者の式を突き合わせて数値で確認した。tilt 係数 (<see cref="UpdateTiltCoefficients"/>) は試料傾斜も畳み込み済みなので、
    /// 結果はそのまま試料系 = <see cref="EbsdMonteCarloDistribution.DirectionToBinCoords"/> に渡せる座標系になる。</para></summary>
    readonly struct ExitRay(double xMirror, double detX, double yCoeff, double zCoeff, double yConst, double zConst)
    {
        public double X(double pxView) => xMirror * pxView + detX;
        public double Y(double pyView) => -(yCoeff * pyView + yConst);
        public double Z(double pyView) => -(zCoeff * pyView + zConst);

        /// <summary>260921Cl 追加: 照射点から検出器面までの垂直距離 D [mm] (= EbsdDetectorGeometry.CameraLength)。
        /// (X, Y, Z) は照射点から検出器上の点への物理的な変位 (の符号反転) で、tilt 係数 (yCoeff, zCoeff) は単位ベクトル、
        /// 定数 (yConst, zConst) は (detY, detZ) の回転なので長さは mm のまま。検出器面の法線 n = (0, zCoeff, −yCoeff) との内積は
        /// pyView に依らず zConst·yCoeff − yConst·zCoeff になる (試料傾斜 0 で |detZ·cosδ − detY·sinδ| = |n·C| と一致)。</summary>
        public double PlaneDistance => Math.Abs(zConst * yCoeff - yConst * zCoeff);
    }

    /// <summary>260921Cl 追加 (/simplify): 現在の tilt 係数と視野から <see cref="ExitRay"/> を作る。</summary>
    ExitRay CreateExitRay(in EbsdRasterView view) => new(view.XMirror, view.DetX, yCoeffPy, zCoeffPy, yConst, zConst);

    /// <summary>260921Cl 追加 (/simplify): 合成ループの作業領域。Parallel.For の localInit で worker ごとに 1 個だけ作る。</summary>
    //260921Cl 変更 (深さ写像 A2): 重み配列の内挿 (C, Bw) をやめ、パラメータの内挿用に N(E)・λ(E) を持たせる
    //旧: sealed class BinScratch(int nSlices)
    //旧: {
    //旧:     public readonly double[] Wv = new double[nSlices];
    //旧:     public readonly double[] C = new double[16];
    //旧:     public readonly double[][] Bw = new double[16][];
    //旧: }
    sealed class BinScratch(int nSlices, int eLen)
    {
        public readonly double[] Wv = new double[nSlices];
        /// <summary>内挿したエネルギー重み N(E) = Σ ω·G_b(E)</summary>
        public readonly double[] NE = new double[eLen];
        /// <summary>内挿した λ_d(E) (N で重み付けした平均)</summary>
        public readonly double[] LamE = new double[eLen];
    }

    /// <summary>260921Cl 追加 (深さ写像 A2): 画素の重みベクトル (エネルギー × 深さ) を、MC のビンパラメータを 4×4 タップで内挿して作る。
    /// 戻り値は同じ重み ω で平均した非晶質割合 fA (amorphousFraction が null、または電子のあるタップが無ければ 0)。
    /// <para>【何が変わったか】旧版はビンごとの<b>重み配列そのもの</b>を内挿していた。ところが MC の源深さはビンの中で
    /// <b>垂直深さ</b>として当てはめてあり、マスターパターンの深さ格子は<b>出射方向の経路長</b> (t = d/μ) なので、
    /// 画素の μ で換算しないと使えない。16×16 の射出半球ビン 1 個の中でも χ は ~11° 動くので (検出器下端では
    /// 1/μ が 2.4 → 4 と変わる)、ビン中心の μ で換算した配列を内挿するのではなく、パラメータを内挿してから
    /// <b>画素の μ で</b> <see cref="EbsdMonteCarloDistribution.FillPathLengthWeights"/> を呼ぶ。</para>
    /// <para>【内挿の重み】タップ t の重みは ω_t = c_t·s_b。model 2 (absolute) は s_b = F_b (電子の割合) なので、
    /// N(E) = Σ ω·G_b(E) は電子数の混合になり、総和は旧版の Σ c·F_b と同じ。model 0/1 は s_b = [F_b &gt; 0] で、
    /// 総和は旧版の「電子のあるビンの係数の和」と同じ (各ビンの重みは総和 1、空ビンは 0)。
    /// λ(E) は N で重み付けした平均 (指数の混合を平均で代表させる近似)。電子の無いビンは ω = 0 なので隣へ漏れない。</para>
    /// <para>⚠ 旧版とは μ = 1 でも一致しない (重みの内挿とパラメータの内挿は非線形の分だけ違う)。</para>
    /// <para>260921Cl (Lambert 等積ディスク): 内挿するのは分布の「場」(Flat* 配列)。F は被覆率で割ったビン全面あたりの値で、
    /// 射出半球の円板から外れたビンは内側から延長してある (EbsdMonteCarloDistribution の ctor の doc【縁のビン】)。</para></summary>
    //260921Cl シグネチャ変更: 非晶質割合はビンの配列 (double[,]) ではなく、分布の場 (FlatAmorphousFraction、縁で延長済み) を使う。
    //旧: static double EvaluatePathLengthWeights(EbsdMonteCarloDistribution dist, double bx, double by, double mu, bool absolute, bool sliceMass,
    //旧:     double[] depths, double[] depthWidths, BinScratch s, int eLen, int nSlices, double[,] amorphousFraction)
    static double EvaluatePathLengthWeights(EbsdMonteCarloDistribution dist, double bx, double by, double mu, bool absolute, bool sliceMass,
        double[] depths, double[] depthWidths, BinScratch s, int eLen, int nSlices, bool withAmorphous)
    {
        int binCount = dist.BinCount;
        BinSplineTaps(bx, binCount, out int x0, out int x1, out int x2, out int x3, out double qx0, out double qx1, out double qx2, out double qx3);
        BinSplineTaps(by, binCount, out int y0, out int y1, out int y2, out int y3, out double qy0, out double qy1, out double qy2, out double qy3);
        Span<int> xs = [x0, x1, x2, x3], ys = [y0, y1, y2, y3];
        Span<double> qx = [qx0, qx1, qx2, qx3], qy = [qy0, qy1, qy2, qy3];
        var n = s.NE.AsSpan(0, eLen); var lam = s.LamE.AsSpan(0, eLen);
        n.Clear(); lam.Clear();
        var flatG = dist.FlatEnergyDistribution; var flatL = dist.FlatLambdaNm; var flatF = dist.FlatFraction;
        var flatA = withAmorphous ? dist.FlatAmorphousFraction : null; //260921Cl 追加
        //260921Cl (Codex 指摘): 非晶質割合も強度と同じ重み ω で平均する (旧: fA = Σ c·fA_b)。
        //  旧式では電子の無い空ビン (fA_b = 0) が隣にあるだけで fA が下がり、全部が非晶質の領域でも結晶の変調が戻った。
        //  model 2 では電子数で重み付けした割合 (= その画素へ来る電子のうち層内に源を持つものの割合) になる
        double fANum = 0, omegaSum = 0;
        for (int ty = 0; ty < 4; ty++) //並びは旧 GatherBinTaps と同じ (x が速い)
            for (int tx = 0; tx < 4; tx++)
            {
                double c = qx[tx] * qy[ty];
                int ix = xs[tx], iy = ys[ty], b = ix * binCount + iy;
                double f = flatF[b];
                double omega = c * (absolute ? f : (f > 0 ? 1.0 : 0.0));
                if (!(omega > 0)) continue;
                //if (amorphousFraction != null) fANum += omega * amorphousFraction[ix, iy]; //260921Cl 変更前 (縁で延長した場を使う)
                if (flatA != null) fANum += omega * flatA[b];
                omegaSum += omega;
                int o = b * eLen;
                for (int e = 0; e < eLen; e++) { double g = omega * flatG[o + e]; n[e] += g; lam[e] += g * flatL[o + e]; }
            }
        for (int e = 0; e < eLen; e++) lam[e] = n[e] > 0 ? lam[e] / n[e] : EbsdMonteCarloDistribution.UniformDepthLambdaNm;
        EbsdMonteCarloDistribution.FillPathLengthWeights(s.Wv.AsSpan(0, nSlices), n, lam, mu, depths, depthWidths, sliceMass);
        return omegaSum > 0 ? fANum / omegaSum : 0;
    }

    /// <summary>260921Cl 追加 (深さ写像 A2、段階 4): <b>検出器の画素で平均した</b> model 2 の重み (エネルギー × 深さの区間質量、長さ eLen·dLen)。
    /// <see cref="EbsdMonteCarloDistribution.ComposeGlobalWeightedPattern"/> に渡すと、ZNCC 照合・E_c 較正用の 1 枚の合成が
    /// 「射出半球全体の平均」ではなく「この検出器画像の画素平均」になる。
    /// <para>【なぜ】A2 以降の重みは画素の μ で経路長へ換算するので方向に強く依存する。半球全体で平均すると、検出器が見ない
    /// 表面すれすれ (μ → 0、経路の長い) 方向まで混ざる。各画素の重みは表示合成 (<see cref="ApplyWeightedModel2"/>) と
    /// 同じ <see cref="EvaluatePathLengthWeights"/> で作るので、この平均は表示合成の重みの画素平均そのもの (F も 1 回だけ入る)。
    /// 試料表面より下を向く画素 (表示合成でも 0) は寄与 0 として分母に数える。</para>
    /// <para>⚠ それでも「重みを平均してから 1 枚の球へ縮約して投影」する近似は残る (重みもマスターの応答も画素方向に依存するため)。
    /// E_c の最終較正は画素ごとの合成で行うこと (Codex 指摘)。これは ZNCC の探索用の近似。</para>
    /// <para>260921Cl 追加 (作者判断): 各画素の重みに<b>画素の立体角</b> dΩ/dA ∝ cos³α (α = 検出器面の法線からの角度) を掛ける。
    /// 表示合成 (model 2) も同じ係数を掛けるので「表示合成の重みの画素平均」のまま、全体としては検出器の立体角で重み付けした
    /// 電子の平均 (= 検出器に当たる電子の平均。旧・検出器ビニングの全ビン和と同じ意味) になる。
    /// 射出半球のビンは等立体角なので、この係数を掛けないと検出器の端 (α が大きく、画素 1 個の立体角が小さい) を過大に数える。</para></summary>
    /// <param name="stride">画素を間引く間隔 (1 で全画素)。重みは画素間でなめらかなので 4 程度で十分。</param>
    public static double[] ComputeDetectorAverageSliceWeights(EbsdMonteCarloDistribution dist, EbsdDetectorGeometry detector, int stride = 4)
    {
        ArgumentNullException.ThrowIfNull(dist);
        ArgumentNullException.ThrowIfNull(detector);
        if (stride < 1) throw new ArgumentOutOfRangeException(nameof(stride));
        int eLen = dist.EnergyCount, dLen = dist.DepthCount, nSlices = eLen * dLen, binCount = dist.BinCount;
        var depths = dist.Depths; var widths = MasterPattern.ComputeDepthIntervals(depths);
        int rows = (detector.HeightPx + stride - 1) / stride;
        //行ごとに合計してから行の順に足す (並列の足し順で結果が揺れないように。E_c の測定を再現できるようにするため)
        var rowSums = new double[rows][]; var rowCounts = new long[rows];
        Parallel.For(0, rows, () => new BinScratch(nSlices, eLen), (ri, _, scratch) =>
        {
            int row = ri * stride;
            var sum = new double[nSlices]; long n = 0;
            for (int col = 0; col < detector.WidthPx; col += stride)
            {
                var v = detector.PixelToSampleDirection(col, row); //視線 = 射出の逆 (EbsdDetectorGeometry の doc)
                double sx = -v.X, sy = -v.Y, sz = -v.Z;
                n++; //260921Cl (Codex 指摘): 地平線より下の画素も数える (寄与 0)。表示合成のその画素の値も 0 なので、平均の定義を画素ごとの合成と揃える
                if (!(sz > 0)) continue; //試料表面より下へ向かう方向には電子が出てこない
                double mu = sz / Math.Sqrt(sx * sx + sy * sy + sz * sz);
                var (bx, by) = EbsdMonteCarloDistribution.DirectionToBinCoords(sx, sy, sz, binCount);
                //EvaluatePathLengthWeights(dist, bx, by, mu, absolute: true, sliceMass: true, depths, widths, scratch, eLen, nSlices, null); //260921Cl 変更前
                EvaluatePathLengthWeights(dist, bx, by, mu, absolute: true, sliceMass: true, depths, widths, scratch, eLen, nSlices, false);
                var wv = scratch.Wv;
                //for (int k = 0; k < nSlices; k++) sum[k] += wv[k]; //260921Cl 変更前
                //260921Cl 追加: 画素の立体角 cos³α (表示合成 model 2 と同じ係数)
                double cosA = detector.CameraLength / detector.PixelToLabPoint(col, row).Length, jac = cosA * cosA * cosA;
                for (int k = 0; k < nSlices; k++) sum[k] += wv[k] * jac;
            }
            rowSums[ri] = sum; rowCounts[ri] = n;
            return scratch;
        }, _ => { });
        var total = new double[nSlices];
        long count = 0;
        for (int ri = 0; ri < rows; ri++) { for (int k = 0; k < nSlices; k++) total[k] += rowSums[ri][k]; count += rowCounts[ri]; }
        if (count > 0) for (int k = 0; k < nSlices; k++) total[k] /= count;
        return total;
    }

    //260921Cl 削除 (深さ写像 A2): 重み配列の内挿は EvaluatePathLengthWeights (パラメータを内挿して画素の μ で換算) に置き換えたので未使用。
    ///// <summary>260921Cl 追加 (/simplify): ビン座標 (bx, by) の 4×4 タップを集め、係数 c と重み配列 bw を 16 個ずつ埋める。
    ///// 並びは x が速い (t = x + 4y) = 旧コードの c00, c10, c20, c30, c01, … と同じ順。
    ///// scalarField を渡すと同じ係数で内挿した値を返す (非晶質割合 fA 用。null なら 0)。</summary>
    //static double GatherBinTaps(double[,][] field, double[,] scalarField, double bx, double by, int binCount, double[][] bw, double[] c)
    //{
    //    BinSplineTaps(bx, binCount, out int x0, out int x1, out int x2, out int x3, out double qx0, out double qx1, out double qx2, out double qx3);
    //    BinSplineTaps(by, binCount, out int y0, out int y1, out int y2, out int y3, out double qy0, out double qy1, out double qy2, out double qy3);
    //    double s = 0;
    //    Put(0, x0, y0, qx0 * qy0); Put(1, x1, y0, qx1 * qy0); Put(2, x2, y0, qx2 * qy0); Put(3, x3, y0, qx3 * qy0);
    //    Put(4, x0, y1, qx0 * qy1); Put(5, x1, y1, qx1 * qy1); Put(6, x2, y1, qx2 * qy1); Put(7, x3, y1, qx3 * qy1);
    //    Put(8, x0, y2, qx0 * qy2); Put(9, x1, y2, qx1 * qy2); Put(10, x2, y2, qx2 * qy2); Put(11, x3, y2, qx3 * qy2);
    //    Put(12, x0, y3, qx0 * qy3); Put(13, x1, y3, qx1 * qy3); Put(14, x2, y3, qx2 * qy3); Put(15, x3, y3, qx3 * qy3);
    //    return s;
    //
    //    void Put(int t, int ix, int iy, double ct)
    //    {
    //        c[t] = ct; bw[t] = field[ix, iy];
    //        if (scalarField != null) s += ct * scalarField[ix, iy];
    //    }
    //}
    //
    ///// <summary>260921Cl 追加 (/simplify): 画素の重みベクトル wv[k] = Σ_t c[t]·bw[t][k] (t = 0..15 の順)。
    ///// <para>⚠ 足す順が旧コードの 16 項和と同じなので出力はビット一致する
    ///// (リファクタ前後で <c>tools/EbsdProfileFit --golden-compose</c> の 12 通りのハッシュが一致することを確認)。
    ///// ⚠ FMA や SIMD に置き換えると丸めが変わってビット一致しなくなる (それ自体は誤差 1 ulp 程度で害は無いが、検証の物差しを失う)。</para></summary>
    //static void EvaluateBinWeights(double[][] bw, double[] c, double[] wv, int n)
    //{
    //    var w = wv.AsSpan(0, n);
    //    w.Clear();
    //    for (int t = 0; t < 16; t++)
    //    {
    //        var b = bw[t].AsSpan(0, n);
    //        double ct = c[t];
    //        for (int k = 0; k < w.Length; k++) w[k] += ct * b[k];
    //    }
    //}

    /// <summary>260921Cl 追加 (/simplify): A(E) 有効時の 4 つの累算器と合成式を 1 か所に (旧: 3 モデルに複製)。
    /// W = その画素・スライスの実効重み、M̄ = 平面ごとの全球方向平均、A = A(E)。
    /// Dc = Σ W·M̄ (= A(E) 無効時の明るさ)、MbarA = Σ W·A·M̄、W = Σ W、SA = Σ W·A。
    /// 合成式 V = Σ W·M̄ + Ā·(Σ W)·(⟨I⟩_A − ⟨M̄⟩_A) の意味は <see cref="MeanCoherentFraction"/> の doc。</summary>
    struct CohAccum
    {
        public double Dc, MbarA, W, SA;
        /// <param name="w">実効重み W</param><param name="m">W·M̄ (model 2 では weight·M̄。M̄ が既に /Δt 済みのため)</param><param name="aE">A(E)</param>
        public void Add(double w, double m, double aE) { Dc += m; MbarA += m * aE; W += w; SA += w * aE; }
        /// <param name="sum">Σ W·A·I</param><param name="meanCohFraction">Ā</param>
        public readonly double Combine(double sum, double meanCohFraction) => SA > 1e-300 ? Dc + meanCohFraction * W / SA * (sum - MbarA) : Dc;
    }

    /// <summary>260921Cl 変更 (作者報告: Forsterite で A(E) を入れるとブロードな暗線が複数出る): 平均の「コヒーレント割合」
    /// Ā = Σ W·A(E) / Σ W (無次元、0〜1)。画素に依らない 1 個のスカラー。
    /// <para>260921Cl 変更 (作者判断): 平均は<b>検出器面で、画素の立体角で重み付けして</b>取る。W は各モデルの画素の重み
    /// (<see cref="EvaluatePathLengthWeights"/>、model 1 は平面の規格化係数も掛ける) × 画素の立体角 dΩ/dA ∝ cos³α。
    /// 旧 (全ビンの平均) は、検出器ビニングだった頃は「検出器に当たる電子の平均」で、E_c = 0.7 keV の較正 (正本 §2.7) もこの定義で行った。
    /// 射出半球のビニングにしてからは「射出半球全体の平均」になり、検出器が見ない方向 (表面すれすれ・後方) まで混ざっていた。
    /// 検出器面で積分すれば旧定義に戻り、表示の視野 (パン・ズーム) にも依存しない。
    /// 検出器面を 48 × 48 点で標本化する (重みは画素間でなめらか)。行ごとに足してから行の順に足すので並列の順序に依らない。</para>
    /// <para>【旧 IncoherentPedestal の何が足りなかったか】旧版は「A(E) で失った分」を強度の次元を持つ 1 個の台座 P として
    /// 足していた。ところが画素の明るさは Σ w(画素)·A(E)·I なので、A(E) を入れた瞬間に**コヒーレント項の総重み
    /// S(画素) = Σ w·A が画素ごとに変わり、それが明るさにそのまま掛かる**。BinWeights はビンごとに Σ = 1 へ正規化されて
    /// いるので A(E) 無効なら S ≡ 1 (= 明るさに 8×8 の痕跡はゼロ) だが、A(E) を入れると S は 8×8 の値を双線形補間した場になる。
    /// Forsterite 20 kV・BSE 45 万本の実測で S は max/min = 3.04、標準偏差/平均 = 33 % (E_c = 0.7 keV。E_c = 0.2 keV でも
    /// 2.67 / 28 %)。定数の台座ではこの面内ムラを打ち消せず、背景平坦化の高域通過を通すと
    /// **ビン間隔 (視野幅/8) を周期とするブロードな暗線**として出る。検証は tools/EbsdProfileFit --mcweights / --bin-artifact。</para>
    /// <para>【なぜ S が暴れるか】A(E) = exp(−(E0−E)/E_c) は E_c が小さいほど E ≈ E0 の最上段だけを見る。ところがビンごとの
    /// エネルギー分布は左右非対称ガウシアンの当てはめで、E0 は平均 (≈ 0.8 E0) から 2〜3σ も外側 = **裾の外挿**である。
    /// 裾の値は当てはめた平均・σ の小さな差に指数関数的に効くので、ビン間のばらつきが激増する。</para>
    /// <para>【新しい合成式】V = Σ W·M̄ + Ā·(Σ W)·(⟨I⟩_A − ⟨M̄⟩_A)。⟨·⟩_A は W·A(E) で重み付けした**正規化**平均、
    /// M̄ は (energy, depth) 平面ごとの全球方向平均、W はその画素の実効重み。第 1 項は A(E) 無効時の明るさそのもので、
    /// 第 2 項は正規化済みなので総重みのスケールに依らない。→ 明るさの 8×8 依存は完全に消え、コントラストの減衰は
    /// 全面で Ā に揃う。バンド幅を決める「エネルギースライス間の相対的な重み配分」は W·A のまま残るので、
    /// E_c の校正値 (正本 §2.7) はそのまま意味を持つ。</para>
    /// <para>⚠ コントラストの面内変化 (取り出し角による低損失電子の割合の違い) は意図的に捨てている。旧版が台座を
    /// スカラーにしたのと同じ判断: 実体は滑らかで小さいはずの効果なのに、8×8 のガウシアン裾の外挿では 3 倍ものムラとして
    /// 出てきてしまい、忠実な表現になっていないため。</para>
    /// <para>【方向平均】半球ごとの平均をそのまま使うと赤道で段差が出るので、両半球の平均 = 全球平均を使う</para></summary>
    //旧シグネチャ: static double IncoherentPedestal(EbsdMonteCarloDistribution dist, double[] cohA, double[] posMeans, double[] negMeans, int eLen, int dLen, bool differential, double[] depthWidths, double[] planeScaleFactors)
    //旧本体 (260920Cl、6e21ac1。/simplify2 でコメントとして復元):
    //{
    //    if (cohA == null || posMeans == null || negMeans == null) return 0; //A(E) 無効 = 台座なし (従来動作)
    //    var g = new double[eLen * dLen];
    //    int nb = 0;
    //    for (int bi = 0; bi < dist.BinCount; bi++)
    //        for (int bj = 0; bj < dist.BinCount; bj++)
    //        {
    //            var bw = differential ? dist.BinAbsoluteSliceWeights[bi, bj] : dist.BinWeights[bi, bj];
    //            if (bw == null) continue;
    //            nb++;
    //            for (int k = 0; k < g.Length && k < bw.Length; k++) g[k] += bw[k];
    //        }
    //    if (nb == 0) return 0;
    //    double p = 0;
    //    for (int ei = 0; ei < eLen; ei++)
    //    {
    //        double mean = 0.5 * (posMeans[ei] + negMeans[ei]); //全球の方向平均
    //        double wSum = 0;
    //        for (int di = 0; di < dLen; di++)
    //        {
    //            int k = ei * dLen + di;
    //            double w = g[k] / nb;
    //            if (differential && depthWidths != null) w /= depthWidths[di]; //model 2 は区間平均 ΔM/Δt に合わせる
    //            if (planeScaleFactors != null) w *= (uint)k < (uint)planeScaleFactors.Length ? planeScaleFactors[k] : 0.0; //model 1 の規格化係数
    //            wSum += w;
    //        }
    //        p += (1 - cohA[ei]) * mean * wSum;
    //    }
    //    return p;
    //}
    //260921Cl (/simplify) シグネチャ変更: differential は depthWidths != null から決まる (GetPlaneMeansCached と同じ規約) ので引数から外した
    //旧: static double MeanCoherentFraction(EbsdMonteCarloDistribution dist, double[] cohA, int eLen, int dLen, bool differential, double[] depthWidths, double[] planeScaleFactors)
    //260921Cl シグネチャ変更 (作者判断: 検出器の立体角で重み付け): 検出器面で積分するので視野 (検出器の寸法) と tilt 係数が要る → インスタンスメソッド。
    //  model の違いは absolute / sliceMass (EvaluatePathLengthWeights と同じ) と planeScaleFactors で渡す。
    //旧: static double MeanCoherentFraction(EbsdMonteCarloDistribution dist, double[] cohA,
    //旧:     int eLen, int dLen, double[] depthWidths, double[] planeScaleFactors)
    //旧: {
    //旧:     if (cohA == null) return 1; //A(E) 無効 (呼び出し側は使わない)
    //旧:     bool differential = depthWidths != null; //model 2 (区間平均 ΔM/Δt と絶対スライス重み) かどうか
    //旧:     var g = new double[eLen * dLen];
    //旧:     int nb = 0;
    //旧:     for (int bi = 0; bi < dist.BinCount; bi++)
    //旧:         for (int bj = 0; bj < dist.BinCount; bj++)
    //旧:         {
    //旧:             var bw = differential ? dist.BinAbsoluteSliceWeights[bi, bj] : dist.BinWeights[bi, bj];
    //旧:             if (bw == null) continue;
    //旧:             nb++;
    //旧:             for (int k = 0; k < g.Length && k < bw.Length; k++) g[k] += bw[k];
    //旧:         }
    //旧:     if (nb == 0) return 1;
    //旧:     double num = 0, den = 0;
    //旧:     for (int ei = 0; ei < eLen; ei++)
    //旧:         for (int di = 0; di < dLen; di++)
    //旧:         {
    //旧:             int k = ei * dLen + di;
    //旧:             double w = g[k] / nb;
    //旧:             //260921Cl 削除 (深さ写像 A2、Codex 指摘): Ā は「電子の割合」なので区間質量のまま集計する。Δt で割ると
    //旧:             //  不等間隔の深さ格子で、同じ区間を細分しただけで Ā が変わってしまう (等間隔なら定数倍で消えるので従来と同じ)。
    //旧:             //旧: if (differential) w /= depthWidths[di]; //model 2 は区間平均 ΔM/Δt に合わせる
    //旧:             if (planeScaleFactors != null) w *= (uint)k < (uint)planeScaleFactors.Length ? planeScaleFactors[k] : 0.0; //model 1 の規格化係数
    //旧:             num += w * cohA[ei];
    //旧:             den += w;
    //旧:         }
    //旧:     return den > 0 ? num / den : 1;
    //旧: }
    double MeanCoherentFraction(EbsdMonteCarloDistribution dist, double[] cohA, in EbsdRasterView view, bool absolute, bool sliceMass,
        double[] depths, double[] depthWidths, double[] planeScaleFactors, int eLen, int dLen)
    {
        if (cohA == null) return 1; //A(E) 無効 (呼び出し側は使わない)
        double hw = view.HalfWidth, hh = view.HalfHeight;
        if (!(hw > 0) || !(hh > 0)) return 1;
        var ray = CreateExitRay(view); //ラムダは in 引数を捕捉できないので値で持つ
        double planeD = ray.PlaneDistance;
        const int n = 48;
        int nSlices = eLen * dLen, binCount = dist.BinCount;
        //Ā は「電子の割合」なので区間質量のまま集計する (深さ写像 A2、Codex 指摘: Δt で割ると不等間隔格子で区間の細分に依存する)
        var rowNum = new double[n]; var rowDen = new double[n];
        Parallel.For(0, n, () => new BinScratch(nSlices, eLen), (iy, _, scratch) =>
        {
            double py = (2.0 * iy + 1 - n) / n * hh;
            double oy = ray.Y(py), oz = ray.Z(py);
            if (!(oz > 0)) return scratch; //試料表面より下へ向かう方向 (寄与 0)
            double num = 0, den = 0;
            for (int ix = 0; ix < n; ix++)
            {
                double ox = ray.X((2.0 * ix + 1 - n) / n * hw);
                double r = Math.Sqrt(ox * ox + oy * oy + oz * oz), cosA = planeD / r, jac = cosA * cosA * cosA; //画素の立体角 ∝ cos³α
                var (bx, by) = EbsdMonteCarloDistribution.DirectionToBinCoords(ox, oy, oz, binCount);
                EvaluatePathLengthWeights(dist, bx, by, oz / r, absolute, sliceMass, depths, depthWidths, scratch, eLen, nSlices, false);
                var wv = scratch.Wv;
                for (int ei = 0; ei < eLen; ei++)
                    for (int di = 0; di < dLen; di++)
                    {
                        int k = ei * dLen + di;
                        double w = wv[k] * jac;
                        if (planeScaleFactors != null) w *= (uint)k < (uint)planeScaleFactors.Length ? planeScaleFactors[k] : 0.0; //model 1 の規格化係数
                        num += w * cohA[ei];
                        den += w;
                    }
            }
            rowNum[iy] = num; rowDen[iy] = den;
            return scratch;
        }, _ => { });
        double sn = 0, sd = 0;
        for (int iy = 0; iy < n; iy++) { sn += rowNum[iy]; sd += rowDen[iy]; }
        return sd > 0 ? sn / sd : 1;
    }

    /// <summary>260921Cl 追加: 損失依存のコントラスト係数 A(E) = exp(−(E0 − E)/E_c) をエネルギースライスごとに返す。
    /// 全部が 1 とみなせる (= 実質無効) なら null を返し、呼び出し側がホットループの分岐ごと省けるようにする。
    /// <para>⚠ <b>A(E) の式と E0 のフォールバック規約は、ここ 1 箇所だけに置くこと。</b>
    /// 表示合成 (<see cref="ApplyWeightedModel"/> 系) と ZNCC の目的関数
    /// (<see cref="EbsdMonteCarloDistribution.ComposeGlobalWeightedPattern"/>) が**同じ重み**を使うことが機能の前提で、
    /// 260921Cl までは同じ式が 2 クラスに別々に書かれていた。片方だけ直すと、症状は
    /// 「探索結果が実測と微妙に合わない」だけなので発見が非常に遅れる。
    /// <see cref="MonteCarlo.ElementIonizationPotentialEv"/> を切り出したのと同じ理由。</para></summary>
    internal static double[] CoherenceFactors(double[] energies, double beamEnergyKeV, double decayKeV)
    {
        if (!(decayKeV > 0) || !double.IsFinite(decayKeV) || energies.Length == 0) return null;
        double e0 = beamEnergyKeV > 0 && double.IsFinite(beamEnergyKeV) ? beamEnergyKeV : energies.Max();
        var a = new double[energies.Length];
        bool any = false;
        for (int i = 0; i < a.Length; i++)
        {
            a[i] = Math.Exp(-Math.Max(0, e0 - energies[i]) / decayKeV);
            if (a[i] < 1 - 1E-12) any = true;
        }
        return any ? a : null;
    }

    /// <summary>260920Cl 追加 / 260921Cl 変更: インスタンスのプロパティ束縛を 1 箇所に閉じるだけの薄い層。
    /// 式そのものは <see cref="CoherenceFactors"/> にある</summary>
    private double[] BuildCoherenceFactors(MasterPattern mp) => CoherenceFactors(mp.Energies, BeamEnergyKeV, CoherenceLossDecayKeV);

    //260921Cl シグネチャ変更: 第 3 要素 sphere (平面ごと = 長さ eLen·dLen の全球方向平均) を追加。A(E) の合成 (CohAccum の M̄) が使う。
    //  ⚠ Pos/Neg は「エネルギーごと」へ畳んだ値 (非晶質層の一様成分用。深さのビン依存を持ち込まないための意図的な畳み込み)。
    //  A(E) の台座はコヒーレント項と (energy, depth) の 1 対 1 で打ち消し合う必要があるので、畳む前の平面ごとの値を使う。
    //旧: private (double[] pos, double[] neg) GetPlaneMeansCached(MasterPattern mp, EbsdMonteCarloDistribution dist, float[][] posPlanes, float[][] negPlanes, int dLen, double[] depthWidths)
    private (double[] pos, double[] neg, double[] sphere) GetPlaneMeansCached(MasterPattern mp, EbsdMonteCarloDistribution dist, float[][] posPlanes, float[][] negPlanes, int dLen, double[] depthWidths)
    {
        bool differential = depthWidths != null;
        if (ReferenceEquals(planeMeanCache.Mp, mp) && ReferenceEquals(planeMeanCache.Dist, dist) && planeMeanCache.Differential == differential && planeMeanCache.Pos != null)
            return (planeMeanCache.Pos, planeMeanCache.Neg, planeMeanCache.Sphere);
        var r = GetPlaneMeans(posPlanes, negPlanes, dLen, depthWidths);
        var g = differential ? dist.GlobalDepthSliceWeights : dist.GlobalDepthWeights;
        var pos = CollapseToEnergyReference(r.pos, dLen, g);
        var neg = CollapseToEnergyReference(r.neg, dLen, g);
        var sphere = new double[r.pos.Length]; //260921Cl 追加: 平面ごとの全球平均 (半球ごとだと赤道で段差になる)
        for (int k = 0; k < sphere.Length; k++) sphere[k] = 0.5 * (r.pos[k] + r.neg[k]);
        planeMeanCache = (mp, dist, differential, pos, neg, sphere);
        return (pos, neg, sphere);
    }

    /// <summary>
    /// 構築済みルックアップテーブルと MC フィッティング結果を使い、
    /// 全エネルギー・深さの加重平均 EBSD パターンを計算する。260325Cl 追加
    /// </summary>
    public unsafe void ApplyWeightedModel0(double[] values, int width, int height, MasterPattern mp, EbsdMonteCarloDistribution dist, in EbsdRasterView view)
    {
        EnsureGridMatches(mp, dist); //260727Cl
        double xm = view.XMirror; // 260718Cl: 左右反転。Parallel.For 前に UI スレッドで捕捉 (ワーカーから checkbox 直読は不可)
        int eLen = mp.Energies.Length, dLen = mp.Depths.Length;
        int binCount = dist.BinCount, gs = lookupGridSize;
        double scaleW = view.ScaleW, scaleH = view.ScaleH, viewOffX = view.OffX, viewOffY = view.OffY; // 260724Cl 追加: ラスター=視野全体化に伴い、検出器正規化±1 は物理位置から算出
        //260921Cl 変更 (旧: double halfW = view.HalfWidth, halfH = view.HalfHeight;)
        //  ビン座標は検出器の正規化座標ではなく試料系の射出方向から作る。画素 → 射出方向の規約は ExitRay 1 か所 (/simplify)
        var ray = CreateExitRay(view);
        //var binField = dist.BinWeights; //260921Cl (/simplify) //260921Cl 変更前 (深さ写像 A2)
        //260921Cl 変更 (深さ写像 A2): 重み配列ではなくビンのパラメータを内挿し、画素の μ で経路長へ換算する (EvaluatePathLengthWeights の doc)
        const bool absoluteWeights = false, sliceMassWeights = false; //model 0/1: 密度 × 区間幅、電子のあるビンは総和 1
        var depthGrid = mp.Depths; var depthGridWidths = mp.DepthIntervals;
        int nSlices = eLen * dLen; //260921Cl (/simplify)
        var (posPlanes, negPlanes) = GetAllPlanes(mp, eLen, dLen);//260718Cl
        //var amorphousFraction = dist.BinAmorphousFraction; // 260919Cl 追加: 表面非晶質層に源を持つ電子の割合 (ビンごと) //260921Cl 変更前 (EvaluatePathLengthWeights が分布の場を直接読む)
        bool hasAmorphous = dist.HasAmorphousLayer; // 260919Cl 追加 (/simplify: 層が無ければ fA の内挿も省く。260921Cl: 旧「双線形補間」→ 3 次 B スプライン)
        //260920Cl 追加: 損失依存のコントラスト係数 A(E)。⚠非晶質層が無いと dist.GlobalDepthWeights は null なので、
        //  台座の基準値 (方向平均) は深さの単純平均になる (CollapseToEnergyReference のフォールバック)。方向に依らない成分なので
        //  バンド幅にも ZNCC にも効かないが、絶対値を論じるときはここが MC 重みで積まれていないことに注意
        //  260921Cl: A(E) は台座をやめて平面ごとの方向平均 (planeMeans、エネルギーへ畳まない) を使うようになったので、この注意はもう当たらない
        //  (畳んだ値 posMeans/negMeans を使うのは非晶質層の変調なし成分だけで、そのときは GlobalDepthWeights がある)
        var cohA = BuildCoherenceFactors(mp);
        var (posMeans, negMeans, planeMeans) = hasAmorphous || cohA != null ? GetPlaneMeansCached(mp, dist, posPlanes, negPlanes, dLen, null) : (null, null, null); // 260919Cl 追加: 変調なし成分用の方向平均 / 260920Cl A(E) でも使う / 260921Cl planeMeans 追加
        //260920Cl 追加: 変調なし成分の基準は**全球**の方向平均。半球ごとの平均 (posMeans / negMeans) をそのまま使うと、
        //  パターンが赤道 (試料系 z = 0、ノモニック投影では直線) を跨ぐ所で段差になる。源の向きを失った電子に半球の区別は無い
        double[] sphereMeans = null;
        if (posMeans != null) { sphereMeans = new double[posMeans.Length]; for (int q = 0; q < sphereMeans.Length; q++) sphereMeans[q] = 0.5 * (posMeans[q] + negMeans[q]); }
        //double incoherentPedestal = IncoherentPedestal(dist, cohA, posMeans, negMeans, eLen, dLen, false, null, null); //260920Cl 追加 //260921Cl 変更前
        //double meanCohFraction = MeanCoherentFraction(dist, cohA, eLen, dLen, null, null); //260921Cl 変更: 定数の台座 → 全面共通のコントラスト減衰 Ā (理由は MeanCoherentFraction の doc) //260921Cl 変更前 (全ビン平均)
        double meanCohFraction = MeanCoherentFraction(dist, cohA, view, absoluteWeights, sliceMassWeights, depthGrid, depthGridWidths, null, eLen, dLen); //260921Cl 変更: 検出器の立体角で重み付けした平均
        bool hasA = cohA != null; //260921Cl 追加 (/simplify: 旧 planeMeansLocal は planeMeans の単なる別名だったので削除)

        //Array.Clear(values); //260725Ch: 下の Parallel.For が全画素を必ず代入するため、描画前の全配列ゼロクリアは不要

        // ピクセルごとの加重合計を並列で計算
        fixed (int* pIdx = lookupIdx)
        fixed (float* pWt = lookupWt)
        fixed (bool* pPosZ = lookupPosZ)
        fixed (double* pVal = values)
        {
            var pIdx0 = pIdx; var pWt0 = pWt; var pPosZ0 = pPosZ; var pVal0 = pVal;
            var isHexGrid = lookupGridType == MasterPattern.Types.Hexagonal; // 260331Cl

            Parallel.For(0, height, () => new BinScratch(nSlices, eLen), (h, _, scratch) => //260921Cl (/simplify): 作業領域は worker ごとに 1 個 (260921Cl 深さ写像 A2: eLen を追加)
            {
                //260921Cl 変更 (作者指示): 検出器面の正規化座標ではなく、**試料系の射出方向**を
                //  射出半球の等積格子へ写してビン座標にする (分布側とまったく同じ写像。260921Cl: Rosca-Lambert 正方 → Lambert 等積ディスク)。
                //  旧: double detNormY = ((2h+1-height)*scaleH + viewOffY)/halfH; double by = (1 - detNormY)*0.5*binCount - 0.5;
                double pyView = (2.0 * h + 1 - height) * scaleH + viewOffY;
                double oy = ray.Y(pyView), oz = ray.Z(pyView); //射出方向の Y, Z (試料系)。規約は ExitRay の doc
                //試料表面より下へ向かう方向には電子が出てこない。oz は h だけで決まるので行ごとに落とす
                if (!(oz > 0)) { for (int w = 0; w < width; w++) pVal0[h * width + w] = 0; return scratch; }

                for (int w = 0; w < width; w++)
                {
                    int i = h * width + w;

                    //260921Cl 変更: 旧 double detNormX = -xm * ((2w+1-width)*scaleW + viewOffX)/halfW; double bx = (detNormX + 1)*0.5*binCount - 0.5;
                    double pxView = (2.0 * w + 1 - width) * scaleW + viewOffX;
                    double ox = ray.X(pxView);
                    //⚠ 分布を作る側とまったく同じ写像を使う (EbsdMonteCarloDistribution.DirectionToBinCoords の doc)
                    var (bx, by) = EbsdMonteCarloDistribution.DirectionToBinCoords(ox, oy, oz, binCount);
                    //260921Cl (/simplify): 4×4 タップの収集・非晶質割合の内挿・重みベクトルの評価を共通化
                    //  (旧: 16 個の係数 c<x><y> と 16 本の bw<x><y> と 16 項和を、3 モデル × 六方/正方の 6 か所に複製していた)
                    //double fA = GatherBinTaps(binField, hasAmorphous ? amorphousFraction : null, bx, by, binCount, scratch.Bw, scratch.C); // 260919Cl 非晶質源の割合 //260921Cl 変更前 (深さ写像 A2)
                    //EvaluateBinWeights(scratch.Bw, scratch.C, scratch.Wv, nSlices); //260921Cl 変更前 (深さ写像 A2)
                    //260921Cl 変更 (深さ写像 A2): 画素の出射方向の μ = cos χ で、垂直深さの分布を経路長の分布へ換算して重みを作る
                    double mu = oz / Math.Sqrt(ox * ox + oy * oy + oz * oz);
                    //double fA = EvaluatePathLengthWeights(dist, bx, by, mu, absoluteWeights, sliceMassWeights, depthGrid, depthGridWidths, scratch, eLen, nSlices, hasAmorphous ? amorphousFraction : null); // 260919Cl 非晶質源の割合 //260921Cl 変更前
                    double fA = EvaluatePathLengthWeights(dist, bx, by, mu, absoluteWeights, sliceMassWeights, depthGrid, depthGridWidths, scratch, eLen, nSlices, hasAmorphous); // 260919Cl 非晶質源の割合 (260921Cl: 分布の場 = 縁で延長済み)
                    var wv = scratch.Wv;

                    // ルックアップテーブルからマスターパターン補間パラメータ取得
                    bool posZ = pPosZ0[i];

                    // 全エネルギー・深さで加重合計
                    double sum = 0;
                    double sumMean = 0; // 260919Cl 追加: 非晶質源 (変調なし) 用の方向平均強度
                    var acc = new CohAccum(); //260921Cl 追加: A(E) 有効時だけ使う累算器 (意味は CohAccum の doc)

                    if (isHexGrid) // 260331Cl: 六方格子
                    {
                        int i3 = i * 3;
                        int hIdx0 = pIdx0[i3], hIdx1 = pIdx0[i3 + 1], hIdx2 = pIdx0[i3 + 2];
                        float hw0 = pWt0[i3], hw1 = pWt0[i3 + 1], hw2 = pWt0[i3 + 2];
                        for (int ei = 0; ei < eLen; ei++)
                        {
                            //260921Cl 変更: aE は ei にしか依らないので di ループの外へ出す (値は完全に同一)。
                            //  旧は最内 (画素 × eLen × dLen) で毎回 null 判定していた
                            double aE = cohA == null ? 1.0 : cohA[ei];
                            //260921Cl: sphereMeans[ei] も ei にしか依らないので一緒に出す (aE を出したときの取りこぼし)
                            double sMean = sphereMeans == null ? 0 : sphereMeans[ei];
                            for (int di = 0; di < dLen; di++)
                            {
                                int wIdx = ei * dLen + di;
                                double weight = wv[wIdx]; //260921Cl (/simplify): 旧 16 項和は EvaluateBinWeights で画素ごとに 1 回だけ評価
                                if (weight < 1e-15) continue;
                                var plane = posZ ? posPlanes[wIdx] : negPlanes[wIdx];//260718Cl 事前展開した配列を参照
                                if (plane == null || plane.Length == 0) continue;
                                sum += weight * (hw0 * plane[hIdx0] + hw1 * plane[hIdx1] + hw2 * plane[hIdx2]) * aE; //260920Cl 変更: A(E) を末尾に掛ける (無効時は 1.0 なので丸めも含めて従来と同一) 
                                if (hasA) acc.Add(weight, weight * planeMeans[wIdx], aE); //260921Cl 追加
                                if (fA > 0) sumMean += weight * sMean; // 260919Cl 追加 / 260920Cl 変更: 半球平均 → 全球平均
                            }
                        }
                    }
                    else // 正方格子
                    {
                        int idx = pIdx0[i];
                        int i2 = i * 2;
                        float mpFw = pWt0[i2], mpFh = pWt0[i2 + 1];
                        float mpW0 = 1 - mpFw, mpW1 = mpFw;
                        float mpFh1 = 1 - mpFh;
                        for (int ei = 0; ei < eLen; ei++)
                        {
                            //260921Cl 変更: aE は ei にしか依らないので di ループの外へ出す (値は完全に同一)。
                            //  旧は最内 (画素 × eLen × dLen) で毎回 null 判定していた
                            double aE = cohA == null ? 1.0 : cohA[ei];
                            //260921Cl: sphereMeans[ei] も ei にしか依らないので一緒に出す (aE を出したときの取りこぼし)
                            double sMean = sphereMeans == null ? 0 : sphereMeans[ei];
                            for (int di = 0; di < dLen; di++)
                            {
                                int wIdx = ei * dLen + di;
                                double weight = wv[wIdx]; //260921Cl (/simplify): 旧 16 項和は EvaluateBinWeights で画素ごとに 1 回だけ評価
                                if (weight < 1e-15) continue;
                                var plane = posZ ? posPlanes[wIdx] : negPlanes[wIdx];//260718Cl 事前展開した配列を参照
                                if (plane == null || plane.Length == 0) continue;
                                double intensity = (mpW0 * plane[idx] + mpW1 * plane[idx + 1]) * mpFh1 + (mpW0 * plane[idx + gs] + mpW1 * plane[idx + gs + 1]) * mpFh;
                                sum += weight * intensity * aE; //260920Cl 変更: A(E) を末尾に掛ける (無効時は 1.0 なので丸めも含めて従来と同一) 
                                if (hasA) acc.Add(weight, weight * planeMeans[wIdx], aE); //260921Cl 追加
                                if (fA > 0) sumMean += weight * sMean; // 260919Cl 追加 / 260920Cl 変更: 半球平均 → 全球平均
                            }
                        }
                    }
                    // pVal0[i] = sum; // 260919Cl 変更前
                    //260920Cl 変更: A(E) で失った分を台座として足し戻す (総量保存)。 //260921Cl 変更前
                    //  ⚠台座は検出器位置に依らない 1 個のスカラーにしてある。理由は IncoherentPedestal の doc を参照 //260921Cl 変更前
                    // double coherentSum = sum + incoherentPedestal; //260920Cl //260921Cl 変更前
                    //260921Cl 変更: 定数の台座では Σ W·A の面内ムラが明るさに残り、8×8 のビン格子がブロードな暗線になった。
                    //  V = Σ W·M̄ + Ā·(Σ W)·(⟨I⟩_A − ⟨M̄⟩_A) へ。第 1 項は A(E) 無効時の明るさそのもの、第 2 項は正規化済み。
                    //  A(E) 無効時は従来どおり sum をそのまま使うので、丸めも含めて数値が一致する。詳細は MeanCoherentFraction の doc
                    double coherentSum = hasA ? acc.Combine(sum, meanCohFraction) : sum; //260921Cl (/simplify: 式は CohAccum.Combine に 1 か所)
                    pVal0[i] = fA > 0 ? (1 - fA) * coherentSum + fA * sumMean : coherentSum; // 260919Cl 変更: 非晶質層内の源は方向平均 (変調なし) で寄与
                }
                return scratch;
            }, _ => { });
        }
    }
    #endregion

    #region masterPatternCombinationModel 1 — globally normalized master
    double[] globalNormalizationFactors = []; // (260325Ch) model 1 用。各 energy/depth slice の全球積算強度を 1 にそろえる係数
    MasterPattern globalNormalizationSource = null; // (260325Ch) 現在の model 1 規格化係数が対応している MasterPattern

    /// <summary>model 1 用に、各 energy/depth slice の全球積算強度 ((+Z) + (-Z)) を 1 にそろえる係数を準備する。260325Ch 追加</summary>
    void EnsureGlobalNormalizationFactorsModel1(MasterPattern mp)
    {
        if (mp == null)
        {
            globalNormalizationFactors = [];
            globalNormalizationSource = null;
            return;
        }

        if (ReferenceEquals(globalNormalizationSource, mp)
            && globalNormalizationFactors.Length == mp.PlaneCount)
            return;

        var factors = new double[mp.PlaneCount];
        for (int planeIndex = 0; planeIndex < mp.PlaneCount; planeIndex++)
        {
            double globalSum = 0.0;

            var positivePlane = (uint)planeIndex < (uint)mp.PositivePlanes.Length ? mp.PositivePlanes[planeIndex] : null;
            if (positivePlane != null)
                for (int i = 0; i < positivePlane.Length; i++)
                    globalSum += positivePlane[i];

            var negativePlane = (uint)planeIndex < (uint)mp.NegativePlanes.Length ? mp.NegativePlanes[planeIndex] : null;
            if (negativePlane != null)
                for (int i = 0; i < negativePlane.Length; i++)
                    globalSum += negativePlane[i];

            factors[planeIndex] = globalSum > 1e-30 ? 1.0 / globalSum : 0.0; // (260325Ch) 全球積算強度が 0 の slice は 0 扱いにする
        }

        globalNormalizationFactors = factors;
        globalNormalizationSource = mp;
    }

    /// <summary>model 1 の単一スライス表示に使う規格化係数 (範囲外は 0)。260726Cl 追加:
    /// 呼び出し側 (FormEBSD) が係数配列そのものを持つ必要が無いよう、準備と参照をここへまとめた</summary>
    public double GetGlobalNormalizationFactorModel1(MasterPattern mp, int planeIndex)
    {
        EnsureGlobalNormalizationFactorsModel1(mp);
        return (uint)planeIndex < (uint)globalNormalizationFactors.Length ? globalNormalizationFactors[planeIndex] : 0.0; // (260325Ch)
    }

    /// <summary>model 1: 各 energy/depth slice の全球積算強度を 1 にそろえてから、単一スライスの EBSD パターンを描く。260325Ch 追加</summary>
    public unsafe void ApplySingleSliceModel1(double[] values, int totalPixels, float[] posPlane, float[] negPlane, double planeScaleFactor = 1.0)
    {
        var gs = lookupGridSize;
        int requiredPlaneLength = checked(gs * gs); //260725Ch
        //if (posPlane == null && negPlane == null) return; //260725Ch 変更前: 欠損スライスで古い values を保持
        if ((posPlane?.Length ?? 0) < requiredPlaneLength && (negPlane?.Length ?? 0) < requiredPlaneLength) { Array.Clear(values); return; } //260725Ch
        var scaleFactor = planeScaleFactor;
        var isHexGrid = lookupGridType == MasterPattern.Types.Hexagonal; // 260331Cl

        fixed (int* pIdx = lookupIdx)
        fixed (float* pWt = lookupWt)
        fixed (bool* pPosZ = lookupPosZ)
        fixed (double* pVal = values)
        fixed (float* pPos = posPlane ?? [])
        fixed (float* pNeg = negPlane ?? [])
        {
            var pIdx0 = pIdx; var pWt0 = pWt; var pPosZ0 = pPosZ;
            var pVal0 = pVal; var pPos0 = pPos; var pNeg0 = pNeg;
            //var hasPos = posPlane != null && posPlane.Length > 0; var hasNeg = negPlane != null && negPlane.Length > 0; //260725Ch 変更前
            var hasPos = posPlane != null && posPlane.Length >= requiredPlaneLength; //260725Ch
            var hasNeg = negPlane != null && negPlane.Length >= requiredPlaneLength;

            if (isHexGrid) // 260331Cl
            {
                Parallel.For(0, totalPixels, i =>
                {
                    float* plane = pPosZ0[i] ? pPos0 : pNeg0;
                    bool hasPlane = pPosZ0[i] ? hasPos : hasNeg;
                    if (!hasPlane) { pVal0[i] = 0; return; }
                    int i3 = i * 3;
                    pVal0[i] = scaleFactor * (pWt0[i3] * plane[pIdx0[i3]]
                             + pWt0[i3 + 1] * plane[pIdx0[i3 + 1]]
                             + pWt0[i3 + 2] * plane[pIdx0[i3 + 2]]);
                });
            }
            else
            {
                Parallel.For(0, totalPixels, i =>
                {
                    float* plane = pPosZ0[i] ? pPos0 : pNeg0;
                    bool hasPlane = pPosZ0[i] ? hasPos : hasNeg;
                    if (!hasPlane) { pVal0[i] = 0; return; }
                    int idx = pIdx0[i];
                    int i2 = i * 2;
                    float fw = pWt0[i2], fh = pWt0[i2 + 1];
                    float w0 = (1 - fw), w1 = fw;
                    pVal0[i] = scaleFactor * ((w0 * plane[idx] + w1 * plane[idx + 1]) * (1 - fh)
                             + (w0 * plane[idx + gs] + w1 * plane[idx + gs + 1]) * fh);
                });
            }
        }
    }

    /// <summary>model 1: 各 energy/depth slice の全球積算強度を 1 にそろえてから weighted 合成する。260325Ch 追加</summary>
    public unsafe void ApplyWeightedModel1(double[] values, int width, int height, MasterPattern mp, EbsdMonteCarloDistribution dist, in EbsdRasterView view)
    {
        EnsureGridMatches(mp, dist); //260727Cl
        double xm = view.XMirror; // 260718Cl: 左右反転 (UI スレッドで捕捉)
        int eLen = mp.Energies.Length, dLen = mp.Depths.Length;
        int binCount = dist.BinCount;
        var gs = lookupGridSize;
        double scaleW = view.ScaleW, scaleH = view.ScaleH, viewOffX = view.OffX, viewOffY = view.OffY; // 260724Cl 追加
        //260921Cl 変更 (旧: double halfW = view.HalfWidth, halfH = view.HalfHeight;)
        //  ビン座標は検出器の正規化座標ではなく試料系の射出方向から作る。画素 → 射出方向の規約は ExitRay 1 か所 (/simplify)
        var ray = CreateExitRay(view);
        //var binField = dist.BinWeights; //260921Cl (/simplify) //260921Cl 変更前 (深さ写像 A2)
        //260921Cl 変更 (深さ写像 A2): 重み配列ではなくビンのパラメータを内挿し、画素の μ で経路長へ換算する (EvaluatePathLengthWeights の doc)
        const bool absoluteWeights = false, sliceMassWeights = false; //model 0/1: 密度 × 区間幅、電子のあるビンは総和 1
        var depthGrid = mp.Depths; var depthGridWidths = mp.DepthIntervals;
        int nSlices = eLen * dLen; //260921Cl (/simplify)
        EnsureGlobalNormalizationFactorsModel1(mp); //260726Cl: 呼び出し側の Ensure 忘れを構造的に不可能にする (係数はキャッシュ済みなら再計算しない)
        var planeScaleFactors = globalNormalizationFactors;
        var (posPlanes, negPlanes) = GetAllPlanes(mp, eLen, dLen);//260718Cl
        //var amorphousFraction = dist.BinAmorphousFraction; // 260919Cl 追加: 表面非晶質層に源を持つ電子の割合 (ビンごと) //260921Cl 変更前 (EvaluatePathLengthWeights が分布の場を直接読む)
        bool hasAmorphous = dist.HasAmorphousLayer; // 260919Cl 追加 (/simplify: 層が無ければ fA の内挿も省く。260921Cl: 旧「双線形補間」→ 3 次 B スプライン)
        //260920Cl 追加: 損失依存のコントラスト係数 A(E)。⚠非晶質層が無いと dist.GlobalDepthWeights は null なので、
        //  台座の基準値 (方向平均) は深さの単純平均になる (CollapseToEnergyReference のフォールバック)。方向に依らない成分なので
        //  バンド幅にも ZNCC にも効かないが、絶対値を論じるときはここが MC 重みで積まれていないことに注意
        //  260921Cl: A(E) は台座をやめて平面ごとの方向平均 (planeMeans、エネルギーへ畳まない) を使うようになったので、この注意はもう当たらない
        //  (畳んだ値 posMeans/negMeans を使うのは非晶質層の変調なし成分だけで、そのときは GlobalDepthWeights がある)
        var cohA = BuildCoherenceFactors(mp);
        var (posMeans, negMeans, planeMeans) = hasAmorphous || cohA != null ? GetPlaneMeansCached(mp, dist, posPlanes, negPlanes, dLen, null) : (null, null, null); // 260919Cl 追加: 変調なし成分用の方向平均 / 260920Cl A(E) でも使う / 260921Cl planeMeans 追加
        //260920Cl 追加: 変調なし成分の基準は**全球**の方向平均。半球ごとの平均 (posMeans / negMeans) をそのまま使うと、
        //  パターンが赤道 (試料系 z = 0、ノモニック投影では直線) を跨ぐ所で段差になる。源の向きを失った電子に半球の区別は無い
        double[] sphereMeans = null;
        if (posMeans != null) { sphereMeans = new double[posMeans.Length]; for (int q = 0; q < sphereMeans.Length; q++) sphereMeans[q] = 0.5 * (posMeans[q] + negMeans[q]); }
        //double incoherentPedestal = IncoherentPedestal(dist, cohA, posMeans, negMeans, eLen, dLen, false, null, planeScaleFactors); //260920Cl 追加 //260921Cl 変更前
        //double meanCohFraction = MeanCoherentFraction(dist, cohA, eLen, dLen, null, planeScaleFactors); //260921Cl 変更: 定数の台座 → 全面共通のコントラスト減衰 Ā (理由は MeanCoherentFraction の doc) //260921Cl 変更前 (全ビン平均)
        double meanCohFraction = MeanCoherentFraction(dist, cohA, view, absoluteWeights, sliceMassWeights, depthGrid, depthGridWidths, planeScaleFactors, eLen, dLen); //260921Cl 変更: 検出器の立体角で重み付けした平均
        bool hasA = cohA != null; //260921Cl 追加 (/simplify: 旧 planeMeansLocal は planeMeans の単なる別名だったので削除)

        //Array.Clear(values); //260725Ch: 全画素上書きのため不要

        fixed (int* pIdx = lookupIdx)
        fixed (float* pWt = lookupWt)
        fixed (bool* pPosZ = lookupPosZ)
        fixed (double* pVal = values)
        {
            var pIdx0 = pIdx; var pWt0 = pWt; var pPosZ0 = pPosZ; var pVal0 = pVal;
            var isHexGrid = lookupGridType == MasterPattern.Types.Hexagonal; // 260331Cl

            Parallel.For(0, height, () => new BinScratch(nSlices, eLen), (h, _, scratch) => //260921Cl (/simplify): 作業領域は worker ごとに 1 個 (260921Cl 深さ写像 A2: eLen を追加)
            {
                //260921Cl 変更 (作者指示): 検出器面の正規化座標ではなく、**試料系の射出方向**を
                //  射出半球の等積格子へ写してビン座標にする (分布側とまったく同じ写像。260921Cl: Rosca-Lambert 正方 → Lambert 等積ディスク)。
                //  旧: double detNormY = ((2h+1-height)*scaleH + viewOffY)/halfH; double by = (1 - detNormY)*0.5*binCount - 0.5;
                double pyView = (2.0 * h + 1 - height) * scaleH + viewOffY;
                double oy = ray.Y(pyView), oz = ray.Z(pyView); //射出方向の Y, Z (試料系)。規約は ExitRay の doc
                //試料表面より下へ向かう方向には電子が出てこない。oz は h だけで決まるので行ごとに落とす
                if (!(oz > 0)) { for (int w = 0; w < width; w++) pVal0[h * width + w] = 0; return scratch; }

                for (int w = 0; w < width; w++)
                {
                    int i = h * width + w;
                    //260921Cl 変更: 旧 double detNormX = -xm * ((2w+1-width)*scaleW + viewOffX)/halfW; double bx = (detNormX + 1)*0.5*binCount - 0.5;
                    double pxView = (2.0 * w + 1 - width) * scaleW + viewOffX;
                    double ox = ray.X(pxView);
                    //⚠ 分布を作る側とまったく同じ写像を使う (EbsdMonteCarloDistribution.DirectionToBinCoords の doc)
                    var (bx, by) = EbsdMonteCarloDistribution.DirectionToBinCoords(ox, oy, oz, binCount);
                    //260921Cl (/simplify): 4×4 タップの収集・非晶質割合の内挿・重みベクトルの評価を共通化
                    //  (旧: 16 個の係数 c<x><y> と 16 本の bw<x><y> と 16 項和を、3 モデル × 六方/正方の 6 か所に複製していた)
                    //double fA = GatherBinTaps(binField, hasAmorphous ? amorphousFraction : null, bx, by, binCount, scratch.Bw, scratch.C); // 260919Cl 非晶質源の割合 //260921Cl 変更前 (深さ写像 A2)
                    //EvaluateBinWeights(scratch.Bw, scratch.C, scratch.Wv, nSlices); //260921Cl 変更前 (深さ写像 A2)
                    //260921Cl 変更 (深さ写像 A2): 画素の出射方向の μ = cos χ で、垂直深さの分布を経路長の分布へ換算して重みを作る
                    double mu = oz / Math.Sqrt(ox * ox + oy * oy + oz * oz);
                    //double fA = EvaluatePathLengthWeights(dist, bx, by, mu, absoluteWeights, sliceMassWeights, depthGrid, depthGridWidths, scratch, eLen, nSlices, hasAmorphous ? amorphousFraction : null); // 260919Cl 非晶質源の割合 //260921Cl 変更前
                    double fA = EvaluatePathLengthWeights(dist, bx, by, mu, absoluteWeights, sliceMassWeights, depthGrid, depthGridWidths, scratch, eLen, nSlices, hasAmorphous); // 260919Cl 非晶質源の割合 (260921Cl: 分布の場 = 縁で延長済み)
                    var wv = scratch.Wv;
                    bool posZ = pPosZ0[i];

                    double sum = 0;
                    double sumMean = 0; // 260919Cl 追加: 非晶質源 (変調なし) 用の方向平均強度
                    var acc = new CohAccum(); //260921Cl 追加: A(E) 有効時だけ使う累算器 (意味は CohAccum の doc)
                    if (isHexGrid) // 260331Cl
                    {
                        int i3 = i * 3;
                        int hIdx0 = pIdx0[i3], hIdx1 = pIdx0[i3 + 1], hIdx2 = pIdx0[i3 + 2];
                        float hw0 = pWt0[i3], hw1 = pWt0[i3 + 1], hw2 = pWt0[i3 + 2];
                        for (int ei = 0; ei < eLen; ei++)
                        {
                            //260921Cl 変更: aE は ei にしか依らないので di ループの外へ出す (値は完全に同一)。
                            //  旧は最内 (画素 × eLen × dLen) で毎回 null 判定していた
                            double aE = cohA == null ? 1.0 : cohA[ei];
                            //260921Cl: sphereMeans[ei] も ei にしか依らないので一緒に出す (aE を出したときの取りこぼし)
                            double sMean = sphereMeans == null ? 0 : sphereMeans[ei];
                            for (int di = 0; di < dLen; di++)
                            {
                                int wIdx = ei * dLen + di;
                                double weight = wv[wIdx]; //260921Cl (/simplify): 旧 16 項和は EvaluateBinWeights で画素ごとに 1 回だけ評価
                                if (weight < 1e-15) continue;
                                double planeScaleFactor = (uint)wIdx < (uint)planeScaleFactors.Length ? planeScaleFactors[wIdx] : 0.0;
                                if (planeScaleFactor < 1e-30) continue;
                                var plane = posZ ? posPlanes[wIdx] : negPlanes[wIdx];//260718Cl 事前展開した配列を参照
                                if (plane == null || plane.Length == 0) continue;
                                sum += weight * (hw0 * plane[hIdx0] + hw1 * plane[hIdx1] + hw2 * plane[hIdx2]) * planeScaleFactor * aE; //260920Cl 変更: A(E) を末尾に掛ける (無効時は 1.0 なので丸めも含めて従来と同一) 
                                if (hasA) { double W = weight * planeScaleFactor; acc.Add(W, W * planeMeans[wIdx], aE); } //260921Cl 追加
                                if (fA > 0) sumMean += weight * sMean * planeScaleFactor; // 260919Cl 追加 / 260920Cl 変更: 半球平均 → 全球平均
                            }
                        }
                    }
                    else
                    {
                        int idx = pIdx0[i];
                        int i2 = i * 2;
                        float mpFw = pWt0[i2], mpFh = pWt0[i2 + 1];
                        float mpW0 = 1 - mpFw, mpW1 = mpFw;
                        float mpFh1 = 1 - mpFh;
                        for (int ei = 0; ei < eLen; ei++)
                        {
                            //260921Cl 変更: aE は ei にしか依らないので di ループの外へ出す (値は完全に同一)。
                            //  旧は最内 (画素 × eLen × dLen) で毎回 null 判定していた
                            double aE = cohA == null ? 1.0 : cohA[ei];
                            //260921Cl: sphereMeans[ei] も ei にしか依らないので一緒に出す (aE を出したときの取りこぼし)
                            double sMean = sphereMeans == null ? 0 : sphereMeans[ei];
                            for (int di = 0; di < dLen; di++)
                            {
                                int wIdx = ei * dLen + di;
                                double weight = wv[wIdx]; //260921Cl (/simplify): 旧 16 項和は EvaluateBinWeights で画素ごとに 1 回だけ評価
                                if (weight < 1e-15) continue;
                                double planeScaleFactor = (uint)wIdx < (uint)planeScaleFactors.Length ? planeScaleFactors[wIdx] : 0.0;
                                if (planeScaleFactor < 1e-30) continue;
                                var plane = posZ ? posPlanes[wIdx] : negPlanes[wIdx];//260718Cl 事前展開した配列を参照
                                if (plane == null || plane.Length == 0) continue;
                                double intensity = (mpW0 * plane[idx] + mpW1 * plane[idx + 1]) * mpFh1
                                                 + (mpW0 * plane[idx + gs] + mpW1 * plane[idx + gs + 1]) * mpFh;
                                sum += weight * intensity * planeScaleFactor * aE; //260920Cl 変更: A(E) を末尾に掛ける (無効時は 1.0 なので丸めも含めて従来と同一) 
                                if (hasA) { double W = weight * planeScaleFactor; acc.Add(W, W * planeMeans[wIdx], aE); } //260921Cl 追加
                                if (fA > 0) sumMean += weight * sMean * planeScaleFactor; // 260919Cl 追加 / 260920Cl 変更: 半球平均 → 全球平均
                            }
                        }
                    }
                    // pVal0[i] = sum; // 260919Cl 変更前
                    //260920Cl 変更: A(E) で失った分を台座として足し戻す (総量保存)。 //260921Cl 変更前
                    //  ⚠台座は検出器位置に依らない 1 個のスカラーにしてある。理由は IncoherentPedestal の doc を参照 //260921Cl 変更前
                    // double coherentSum = sum + incoherentPedestal; //260920Cl //260921Cl 変更前
                    //260921Cl 変更: 定数の台座では Σ W·A の面内ムラが明るさに残り、8×8 のビン格子がブロードな暗線になった。
                    //  V = Σ W·M̄ + Ā·(Σ W)·(⟨I⟩_A − ⟨M̄⟩_A) へ。第 1 項は A(E) 無効時の明るさそのもの、第 2 項は正規化済み。
                    //  A(E) 無効時は従来どおり sum をそのまま使うので、丸めも含めて数値が一致する。詳細は MeanCoherentFraction の doc
                    double coherentSum = hasA ? acc.Combine(sum, meanCohFraction) : sum; //260921Cl (/simplify: 式は CohAccum.Combine に 1 か所)
                    pVal0[i] = fA > 0 ? (1 - fA) * coherentSum + fA * sumMean : coherentSum; // 260919Cl 変更: 非晶質層内の源は方向平均 (変調なし) で寄与
                }
                return scratch;
            }, _ => { });
        }
    }
    #endregion

    #region masterPatternCombinationModel 2 — absolute MC x differential master
    /// <summary>Model 2: depthIndex と depthIndex-1 の差分を取り、単一 depth slice の EBSD パターンとして描く。260325Ch 追加</summary>
    public unsafe void ApplySingleSliceModel2(double[] values, int totalPixels, float[] posPlane, float[] negPlane, float[] posPlanePrevious = null, float[] negPlanePrevious = null)
    {
        var gs = lookupGridSize;
        int requiredPlaneLength = checked(gs * gs); //260725Ch
        //if (posPlane == null && negPlane == null) return; //260725Ch 変更前: 欠損スライスで古い values を保持
        if ((posPlane?.Length ?? 0) < requiredPlaneLength && (negPlane?.Length ?? 0) < requiredPlaneLength) { Array.Clear(values); return; } //260725Ch
        var isHexGrid = lookupGridType == MasterPattern.Types.Hexagonal; // 260331Cl

        fixed (int* pIdx = lookupIdx)
        fixed (float* pWt = lookupWt)
        fixed (bool* pPosZ = lookupPosZ)
        fixed (double* pVal = values)
        fixed (float* pPos = posPlane ?? [])
        fixed (float* pNeg = negPlane ?? [])
        fixed (float* pPosPrev = posPlanePrevious ?? [])
        fixed (float* pNegPrev = negPlanePrevious ?? [])
        {
            var pIdx0 = pIdx; var pWt0 = pWt; var pPosZ0 = pPosZ;
            var pVal0 = pVal; var pPos0 = pPos; var pNeg0 = pNeg; var pPosPrev0 = pPosPrev; var pNegPrev0 = pNegPrev;
            //var hasPos = posPlane != null && posPlane.Length > 0; var hasNeg = negPlane != null && negPlane.Length > 0; //260725Ch 変更前
            //var hasPosPrev = posPlanePrevious != null && posPlanePrevious.Length > 0; var hasNegPrev = negPlanePrevious != null && negPlanePrevious.Length > 0;
            var hasPos = posPlane != null && posPlane.Length >= requiredPlaneLength; //260725Ch
            var hasNeg = negPlane != null && negPlane.Length >= requiredPlaneLength;
            var hasPosPrev = posPlanePrevious != null && posPlanePrevious.Length >= requiredPlaneLength;
            var hasNegPrev = negPlanePrevious != null && negPlanePrevious.Length >= requiredPlaneLength;

            if (isHexGrid) // 260331Cl
            {
                Parallel.For(0, totalPixels, i =>
                {
                    float* plane = pPosZ0[i] ? pPos0 : pNeg0;
                    float* planePrev = pPosZ0[i] ? pPosPrev0 : pNegPrev0;
                    bool hasPlane = pPosZ0[i] ? hasPos : hasNeg;
                    bool hasPlanePrev = pPosZ0[i] ? hasPosPrev : hasNegPrev;
                    if (!hasPlane) { pVal0[i] = 0; return; }
                    int i3 = i * 3;
                    double intensity = pWt0[i3] * plane[pIdx0[i3]]
                                     + pWt0[i3 + 1] * plane[pIdx0[i3 + 1]]
                                     + pWt0[i3 + 2] * plane[pIdx0[i3 + 2]];
                    if (hasPlanePrev)
                        intensity -= pWt0[i3] * planePrev[pIdx0[i3]]
                                   + pWt0[i3 + 1] * planePrev[pIdx0[i3 + 1]]
                                   + pWt0[i3 + 2] * planePrev[pIdx0[i3 + 2]];
                    pVal0[i] = Math.Max(0.0, intensity);
                });
            }
            else
            {
                Parallel.For(0, totalPixels, i =>
                {
                    float* plane = pPosZ0[i] ? pPos0 : pNeg0;
                    float* planePrev = pPosZ0[i] ? pPosPrev0 : pNegPrev0;
                    bool hasPlane = pPosZ0[i] ? hasPos : hasNeg;
                    bool hasPlanePrev = pPosZ0[i] ? hasPosPrev : hasNegPrev;
                    if (!hasPlane) { pVal0[i] = 0; return; }
                    int idx = pIdx0[i];
                    int i2 = i * 2;
                    float fw = pWt0[i2], fh = pWt0[i2 + 1];
                    float w0 = (1 - fw), w1 = fw;
                    double intensity = (w0 * plane[idx] + w1 * plane[idx + 1]) * (1 - fh)
                                     + (w0 * plane[idx + gs] + w1 * plane[idx + gs + 1]) * fh;
                    if (hasPlanePrev)
                        intensity -= (w0 * planePrev[idx] + w1 * planePrev[idx + 1]) * (1 - fh)
                                   + (w0 * planePrev[idx + gs] + w1 * planePrev[idx + gs + 1]) * fh;
                    pVal0[i] = Math.Max(0.0, intensity);
                });
            }
        }
    }

    /// <summary>model 2: absolute MC 重みと differential MasterPattern を掛け合わせて weighted 合成する。260325Ch 追加</summary>
    public unsafe void ApplyWeightedModel2(double[] values, int width, int height, MasterPattern mp, EbsdMonteCarloDistribution dist, in EbsdRasterView view)
    {
        EnsureGridMatches(mp, dist); //260727Cl
        double xm = view.XMirror; // 260718Cl: 左右反転 (UI スレッドで捕捉)
        int eLen = mp.Energies.Length, dLen = mp.Depths.Length;
        int binCount = dist.BinCount;
        var gs = lookupGridSize;
        double scaleW = view.ScaleW, scaleH = view.ScaleH, viewOffX = view.OffX, viewOffY = view.OffY; // 260724Cl 追加
        //260921Cl 変更 (旧: double halfW = view.HalfWidth, halfH = view.HalfHeight;)
        //  ビン座標は検出器の正規化座標ではなく試料系の射出方向から作る。画素 → 射出方向の規約は ExitRay 1 か所 (/simplify)
        var ray = CreateExitRay(view);
        //var binField = dist.BinAbsoluteSliceWeights; //260921Cl (/simplify) //260921Cl 変更前 (深さ写像 A2)
        //260921Cl 変更 (深さ写像 A2): 重み配列ではなくビンのパラメータを内挿し、画素の μ で経路長へ換算する (EvaluatePathLengthWeights の doc)
        const bool absoluteWeights = true, sliceMassWeights = true; //model 2: 区間質量、総和は電子の割合 F
        //260921Cl 追加 (作者判断): 画素の立体角 dΩ/dA = D/|r|³ = cos³α/D²。射出半球のビンは等立体角なので、ビンの重みは「立体角あたり」。
        //  検出器の画素 1 個が受ける電子はそれに画素の立体角を掛けたもの (旧・検出器面の等面積ビンでは暗黙に入っていた)。
        //  定数 1/D² は落として cos³α だけ掛ける (パターン中心で 1。表示は自動伸張、ZNCC は定数倍に不変)
        double planeD = ray.PlaneDistance;
        var depthGrid = mp.Depths; var depthGridWidths = mp.DepthIntervals;
        int nSlices = eLen * dLen; //260921Cl (/simplify)
        var (posPlanes, negPlanes) = GetAllPlanes(mp, eLen, dLen);//260718Cl
        //var amorphousFraction = dist.BinAmorphousFraction; // 260919Cl 追加: 表面非晶質層に源を持つ電子の割合 (ビンごと) //260921Cl 変更前 (EvaluatePathLengthWeights が分布の場を直接読む)
        bool hasAmorphous = dist.HasAmorphousLayer; // 260919Cl 追加 (/simplify: 層が無ければ fA の内挿も省く。260921Cl: 旧「双線形補間」→ 3 次 B スプライン)
        double[] posMeans = null, negMeans = null, planeMeans = null; // 260919Cl 追加: model 2 は差分 ΔM/Δt の平均 (depthWidths 確定後に取得) / 260921Cl planeMeans 追加
        //260726Cl 追加 (正本 §1.4): plane は累積 M(t) なので隣接差は区間積分。区間平均 R̄=ΔM/Δt にするため区間幅で割る
        //(MC 側の重みは区間質量なので割らない)。等間隔グリッドでは全体が定数倍だが、不等間隔では区間ごとの重み比が変わる
        var depthWidths = mp.DepthIntervals;
        //260920Cl 追加: 損失依存のコントラスト係数 A(E)。⚠非晶質層が無いと dist.GlobalDepthWeights は null なので、
        //  台座の基準値 (方向平均) は深さの単純平均になる (CollapseToEnergyReference のフォールバック)。方向に依らない成分なので
        //  バンド幅にも ZNCC にも効かないが、絶対値を論じるときはここが MC 重みで積まれていないことに注意
        //  260921Cl: A(E) は台座をやめて平面ごとの方向平均 (planeMeans、エネルギーへ畳まない) を使うようになったので、この注意はもう当たらない
        //  (畳んだ値 posMeans/negMeans を使うのは非晶質層の変調なし成分だけで、そのときは GlobalDepthWeights がある)
        var cohA = BuildCoherenceFactors(mp);
        if (hasAmorphous || cohA != null) (posMeans, negMeans, planeMeans) = GetPlaneMeansCached(mp, dist, posPlanes, negPlanes, dLen, depthWidths); // 260919Cl 追加: model 2 は差分 ΔM/Δt の平均 (/simplify: 以前は null 版を先に呼んでキャッシュを取りこぼしていた) / 260920Cl A(E) でも使う
        //260920Cl 追加: 変調なし成分の基準は**全球**の方向平均。半球ごとの平均 (posMeans / negMeans) をそのまま使うと、
        //  パターンが赤道 (試料系 z = 0、ノモニック投影では直線) を跨ぐ所で段差になる。源の向きを失った電子に半球の区別は無い
        double[] sphereMeans = null;
        if (posMeans != null) { sphereMeans = new double[posMeans.Length]; for (int q = 0; q < sphereMeans.Length; q++) sphereMeans[q] = 0.5 * (posMeans[q] + negMeans[q]); }
        //double incoherentPedestal = IncoherentPedestal(dist, cohA, posMeans, negMeans, eLen, dLen, true, depthWidths, null); //260920Cl 追加 //260921Cl 変更前
        //double meanCohFraction = MeanCoherentFraction(dist, cohA, eLen, dLen, depthWidths, null); //260921Cl 変更: 定数の台座 → 全面共通のコントラスト減衰 Ā (理由は MeanCoherentFraction の doc) //260921Cl 変更前 (全ビン平均)
        double meanCohFraction = MeanCoherentFraction(dist, cohA, view, absoluteWeights, sliceMassWeights, depthGrid, depthGridWidths, null, eLen, dLen); //260921Cl 変更: 検出器の立体角で重み付けした平均
        bool hasA = cohA != null; //260921Cl 追加 (/simplify: 旧 planeMeansLocal は planeMeans の単なる別名だったので削除)

        //Array.Clear(values); //260725Ch: 全画素上書きのため不要

        fixed (int* pIdx = lookupIdx)
        fixed (float* pWt = lookupWt)
        fixed (bool* pPosZ = lookupPosZ)
        fixed (double* pVal = values)
        {
            var pIdx0 = pIdx; var pWt0 = pWt; var pPosZ0 = pPosZ; var pVal0 = pVal;
            var isHexGrid = lookupGridType == MasterPattern.Types.Hexagonal; // 260331Cl

            Parallel.For(0, height, () => new BinScratch(nSlices, eLen), (h, _, scratch) => //260921Cl (/simplify): 作業領域は worker ごとに 1 個 (260921Cl 深さ写像 A2: eLen を追加)
            {
                //260921Cl 変更 (作者指示): 検出器面の正規化座標ではなく、**試料系の射出方向**を
                //  射出半球の等積格子へ写してビン座標にする (分布側とまったく同じ写像。260921Cl: Rosca-Lambert 正方 → Lambert 等積ディスク)。
                //  旧: double detNormY = ((2h+1-height)*scaleH + viewOffY)/halfH; double by = (1 - detNormY)*0.5*binCount - 0.5;
                double pyView = (2.0 * h + 1 - height) * scaleH + viewOffY;
                double oy = ray.Y(pyView), oz = ray.Z(pyView); //射出方向の Y, Z (試料系)。規約は ExitRay の doc
                //試料表面より下へ向かう方向には電子が出てこない。oz は h だけで決まるので行ごとに落とす
                if (!(oz > 0)) { for (int w = 0; w < width; w++) pVal0[h * width + w] = 0; return scratch; }

                for (int w = 0; w < width; w++)
                {
                    int i = h * width + w;
                    //260921Cl 変更: 旧 double detNormX = -xm * ((2w+1-width)*scaleW + viewOffX)/halfW; double bx = (detNormX + 1)*0.5*binCount - 0.5;
                    double pxView = (2.0 * w + 1 - width) * scaleW + viewOffX;
                    double ox = ray.X(pxView);
                    //⚠ 分布を作る側とまったく同じ写像を使う (EbsdMonteCarloDistribution.DirectionToBinCoords の doc)
                    var (bx, by) = EbsdMonteCarloDistribution.DirectionToBinCoords(ox, oy, oz, binCount);
                    //260921Cl (/simplify): 4×4 タップの収集・非晶質割合の内挿・重みベクトルの評価を共通化
                    //  (旧: 16 個の係数 c<x><y> と 16 本の bw<x><y> と 16 項和を、3 モデル × 六方/正方の 6 か所に複製していた)
                    //double fA = GatherBinTaps(binField, hasAmorphous ? amorphousFraction : null, bx, by, binCount, scratch.Bw, scratch.C); // 260919Cl 非晶質源の割合 //260921Cl 変更前 (深さ写像 A2)
                    //EvaluateBinWeights(scratch.Bw, scratch.C, scratch.Wv, nSlices); //260921Cl 変更前 (深さ写像 A2)
                    //260921Cl 変更 (深さ写像 A2): 画素の出射方向の μ = cos χ で、垂直深さの分布を経路長の分布へ換算して重みを作る
                    //double mu = oz / Math.Sqrt(ox * ox + oy * oy + oz * oz); //260921Cl 変更前 (|r| を画素の立体角にも使う)
                    double r = Math.Sqrt(ox * ox + oy * oy + oz * oz), mu = oz / r;
                    double cosA = planeD / r, pixelSolidAngle = cosA * cosA * cosA; //260921Cl 追加: 画素の立体角 ∝ cos³α (planeD の doc)
                    //double fA = EvaluatePathLengthWeights(dist, bx, by, mu, absoluteWeights, sliceMassWeights, depthGrid, depthGridWidths, scratch, eLen, nSlices, hasAmorphous ? amorphousFraction : null); // 260919Cl 非晶質源の割合 //260921Cl 変更前
                    double fA = EvaluatePathLengthWeights(dist, bx, by, mu, absoluteWeights, sliceMassWeights, depthGrid, depthGridWidths, scratch, eLen, nSlices, hasAmorphous); // 260919Cl 非晶質源の割合 (260921Cl: 分布の場 = 縁で延長済み)
                    var wv = scratch.Wv;
                    bool posZ = pPosZ0[i];

                    double sum = 0;
                    double sumMean = 0; // 260919Cl 追加: 非晶質源 (変調なし) 用の方向平均強度
                    var acc = new CohAccum(); //260921Cl 追加: A(E) 有効時だけ使う累算器 (意味は CohAccum の doc)
                    if (isHexGrid) // 260331Cl
                    {
                        int i3 = i * 3;
                        int hIdx0 = pIdx0[i3], hIdx1 = pIdx0[i3 + 1], hIdx2 = pIdx0[i3 + 2];
                        float hw0 = pWt0[i3], hw1 = pWt0[i3 + 1], hw2 = pWt0[i3 + 2];
                        for (int ei = 0; ei < eLen; ei++)
                        {
                            //260921Cl 変更: aE は ei にしか依らないので di ループの外へ出す (値は完全に同一)。
                            //  旧は最内 (画素 × eLen × dLen) で毎回 null 判定していた
                            double aE = cohA == null ? 1.0 : cohA[ei];
                            //260921Cl: sphereMeans[ei] も ei にしか依らないので一緒に出す (aE を出したときの取りこぼし)
                            double sMean = sphereMeans == null ? 0 : sphereMeans[ei];
                            for (int di = 0; di < dLen; di++)
                            {
                                int wIdx = ei * dLen + di;
                                double weight = wv[wIdx]; //260921Cl (/simplify): 旧 16 項和は EvaluateBinWeights で画素ごとに 1 回だけ評価
                                if (weight < 1e-15) continue;
                                var plane = posZ ? posPlanes[wIdx] : negPlanes[wIdx];//260718Cl 事前展開した配列を参照
                                if (plane == null || plane.Length == 0) continue;
                                var planePrevious = di > 0 ? posZ ? posPlanes[wIdx - 1] : negPlanes[wIdx - 1] : null;//260718Cl 事前展開した配列を参照 (di>0 なら wIdx-1 = 同 energy の di-1)
                                double intensity = hw0 * plane[hIdx0] + hw1 * plane[hIdx1] + hw2 * plane[hIdx2];
                                if (planePrevious != null && planePrevious.Length > 0)
                                    intensity -= hw0 * planePrevious[hIdx0] + hw1 * planePrevious[hIdx1] + hw2 * planePrevious[hIdx2];
                                sum += weight * Math.Max(0.0, intensity) / depthWidths[di] * aE; //260920Cl 変更: A(E) を末尾に掛ける (無効時は 1.0 なので丸めも含めて従来と同一) //260726Cl: 区間平均 ΔM/Δt
                                //if (hasA) acc.Add(weight / depthWidths[di], weight * planeMeans[wIdx], aE); //260921Cl 追加 (planeMeans は既に /Δt 済みなので m = weight·M̄) //260921Cl 変更前 (深さ写像 A2)
                                //260921Cl 変更 (深さ写像 A2、Codex 指摘): W は区間質量のまま (Δt で割ると不等間隔格子で Ā·W/SA が格子の切り方に依存する。等間隔なら従来と同じ)
                                if (hasA) acc.Add(weight, weight * planeMeans[wIdx], aE); //planeMeans は既に /Δt 済みなので m = weight·M̄
                                if (fA > 0) sumMean += weight * sMean; // 260919Cl 追加 / 260920Cl 変更: 半球平均 → 全球平均 (平均は既に /Δt 済み)
                            }
                        }
                    }
                    else
                    {
                        int idx = pIdx0[i];
                        int i2 = i * 2;
                        float mpFw = pWt0[i2], mpFh = pWt0[i2 + 1];
                        float mpW0 = 1 - mpFw, mpW1 = mpFw;
                        float mpFh1 = 1 - mpFh;
                        for (int ei = 0; ei < eLen; ei++)
                        {
                            //260921Cl 変更: aE は ei にしか依らないので di ループの外へ出す (値は完全に同一)。
                            //  旧は最内 (画素 × eLen × dLen) で毎回 null 判定していた
                            double aE = cohA == null ? 1.0 : cohA[ei];
                            //260921Cl: sphereMeans[ei] も ei にしか依らないので一緒に出す (aE を出したときの取りこぼし)
                            double sMean = sphereMeans == null ? 0 : sphereMeans[ei];
                            for (int di = 0; di < dLen; di++)
                            {
                                int wIdx = ei * dLen + di;
                                double weight = wv[wIdx]; //260921Cl (/simplify): 旧 16 項和は EvaluateBinWeights で画素ごとに 1 回だけ評価
                                if (weight < 1e-15) continue;
                                var plane = posZ ? posPlanes[wIdx] : negPlanes[wIdx];//260718Cl 事前展開した配列を参照
                                if (plane == null || plane.Length == 0) continue;
                                var planePrevious = di > 0 ? posZ ? posPlanes[wIdx - 1] : negPlanes[wIdx - 1] : null;//260718Cl 事前展開した配列を参照 (di>0 なら wIdx-1 = 同 energy の di-1)
                                double intensity = (mpW0 * plane[idx] + mpW1 * plane[idx + 1]) * mpFh1
                                                 + (mpW0 * plane[idx + gs] + mpW1 * plane[idx + gs + 1]) * mpFh;
                                if (planePrevious != null && planePrevious.Length > 0)
                                    intensity -= (mpW0 * planePrevious[idx] + mpW1 * planePrevious[idx + 1]) * mpFh1
                                             + (mpW0 * planePrevious[idx + gs] + mpW1 * planePrevious[idx + gs + 1]) * mpFh;
                                sum += weight * Math.Max(0.0, intensity) / depthWidths[di] * aE; //260920Cl 変更: A(E) を末尾に掛ける (無効時は 1.0 なので丸めも含めて従来と同一) //260726Cl: 区間平均 ΔM/Δt
                                //if (hasA) acc.Add(weight / depthWidths[di], weight * planeMeans[wIdx], aE); //260921Cl 追加 (planeMeans は既に /Δt 済みなので m = weight·M̄) //260921Cl 変更前 (深さ写像 A2)
                                //260921Cl 変更 (深さ写像 A2、Codex 指摘): W は区間質量のまま (Δt で割ると不等間隔格子で Ā·W/SA が格子の切り方に依存する。等間隔なら従来と同じ)
                                if (hasA) acc.Add(weight, weight * planeMeans[wIdx], aE); //planeMeans は既に /Δt 済みなので m = weight·M̄
                                if (fA > 0) sumMean += weight * sMean; // 260919Cl 追加 / 260920Cl 変更: 半球平均 → 全球平均 (平均は既に /Δt 済み)
                            }
                        }
                    }
                    // pVal0[i] = sum; // 260919Cl 変更前
                    //260920Cl 変更: A(E) で失った分を台座として足し戻す (総量保存)。 //260921Cl 変更前
                    //  ⚠台座は検出器位置に依らない 1 個のスカラーにしてある。理由は IncoherentPedestal の doc を参照 //260921Cl 変更前
                    // double coherentSum = sum + incoherentPedestal; //260920Cl //260921Cl 変更前
                    //260921Cl 変更: 定数の台座では Σ W·A の面内ムラが明るさに残り、8×8 のビン格子がブロードな暗線になった。
                    //  V = Σ W·M̄ + Ā·(Σ W)·(⟨I⟩_A − ⟨M̄⟩_A) へ。第 1 項は A(E) 無効時の明るさそのもの、第 2 項は正規化済み。
                    //  A(E) 無効時は従来どおり sum をそのまま使うので、丸めも含めて数値が一致する。詳細は MeanCoherentFraction の doc
                    double coherentSum = hasA ? acc.Combine(sum, meanCohFraction) : sum; //260921Cl (/simplify: 式は CohAccum.Combine に 1 か所)
                    //pVal0[i] = fA > 0 ? (1 - fA) * coherentSum + fA * sumMean : coherentSum; // 260919Cl 変更: 非晶質層内の源は方向平均 (変調なし) で寄与 //260921Cl 変更前
                    //260921Cl 変更: 画素の立体角を掛ける (値は重みについて 1 次同次なので、最後に 1 回掛ければ全項に掛けたのと同じ)
                    pVal0[i] = pixelSolidAngle * (fA > 0 ? (1 - fA) * coherentSum + fA * sumMean : coherentSum); // 260919Cl: 非晶質層内の源は方向平均 (変調なし) で寄与
                }
                return scratch;
            }, _ => { });
        }
    }
    #endregion
}
