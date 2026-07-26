#region using
using System;
using System.Threading.Tasks;
#endregion

namespace Crystallography;

/// <summary>
/// ラスター (width×height) のピクセル中心を表示パターン座標 (mm、検出器中心基準) へ写す係数と、検出器の物理サイズ。260726Cl 追加。
/// px_view = (2w+1-width)·ScaleW + OffX、py_view = (2h+1-height)·ScaleH + OffY (OffX/OffY = 表示のパン量)。
/// HalfWidth/HalfHeight は検出器の物理半幅・半高 (mm) で、MC ビン補間の正規化 (検出器外は端ビンへ外挿) に使う。
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

    /// <summary>
    /// 構築済みルックアップテーブルと MC フィッティング結果を使い、
    /// 全エネルギー・深さの加重平均 EBSD パターンを計算する。260325Cl 追加
    /// </summary>
    public unsafe void ApplyWeightedModel0(double[] values, int width, int height, MasterPattern mp, EbsdMonteCarloDistribution dist, in EbsdRasterView view)
    {
        double xm = view.XMirror; // 260718Cl: 左右反転。Parallel.For 前に UI スレッドで捕捉 (ワーカーから checkbox 直読は不可)
        int eLen = mp.Energies.Length, dLen = mp.Depths.Length;
        int binCount = dist.BinCount, gs = lookupGridSize;
        double scaleW = view.ScaleW, scaleH = view.ScaleH, viewOffX = view.OffX, viewOffY = view.OffY; // 260724Cl 追加: ラスター=視野全体化に伴い、検出器正規化±1 は物理位置から算出
        double halfW = view.HalfWidth, halfH = view.HalfHeight; // 260724Cl 追加
        var (posPlanes, negPlanes) = GetAllPlanes(mp, eLen, dLen);//260718Cl

        //Array.Clear(values); //260725Ch: 下の Parallel.For が全画素を必ず代入するため、描画前の全配列ゼロクリアは不要

        // ピクセルごとの加重合計を並列で計算
        fixed (int* pIdx = lookupIdx)
        fixed (float* pWt = lookupWt)
        fixed (bool* pPosZ = lookupPosZ)
        fixed (double* pVal = values)
        {
            var pIdx0 = pIdx; var pWt0 = pWt; var pPosZ0 = pPosZ; var pVal0 = pVal;
            var isHexGrid = lookupGridType == MasterPattern.Types.Hexagonal; // 260331Cl

            Parallel.For(0, height, h =>
            {
                // この行のピクセルの検出器 Y 座標 (ビン補間用)
                // 260325Cl: スクリーン h=0 → pyFactor≈-DetR (検出器底) → detNormY≈-1, 符号反転しない
                // double detNormY = (2.0 * h + 1 - height) / (double)height; // 260724Cl 変更前: ラスター=検出器全面が前提
                double detNormY = ((2.0 * h + 1 - height) * scaleH + viewOffY) / halfH; // 260724Cl: 物理位置/halfH (検出器外は端ビンへクランプ外挿)
                double by = (1 - detNormY) * 0.5 * binCount - 0.5;
                int bj0 = Math.Clamp((int)Math.Floor(by), 0, binCount - 2);
                double fy = Math.Clamp(by - bj0, 0, 1);

                for (int w = 0; w < width; w++)
                {
                    int i = h * width + w;

                    // 検出器 X 座標
                    // double detNormX = -xm * (2.0 * w + 1 - width) / (double)width; // 260325Cl: スクリーン X は検出器面 X と反転 (BuildLookupTable で -Ri.E11 を使用) / 260718Cl: 左右反転 xm を掛ける // 260724Cl 変更前
                    double detNormX = -xm * ((2.0 * w + 1 - width) * scaleW + viewOffX) / halfW; // 260724Cl: 物理位置/halfW
                    double bx = (detNormX + 1) * 0.5 * binCount - 0.5;
                    int bi0 = Math.Clamp((int)Math.Floor(bx), 0, binCount - 2);
                    double fx = Math.Clamp(bx - bi0, 0, 1);

                    // ビン重みのバイリニア補間係数
                    double c00 = (1 - fx) * (1 - fy), c10 = fx * (1 - fy), c01 = (1 - fx) * fy, c11 = fx * fy;

                    var bw00 = dist.BinWeights[bi0, bj0];
                    var bw10 = dist.BinWeights[bi0 + 1, bj0];
                    var bw01 = dist.BinWeights[bi0, bj0 + 1];
                    var bw11 = dist.BinWeights[bi0 + 1, bj0 + 1];

                    // ルックアップテーブルからマスターパターン補間パラメータ取得
                    bool posZ = pPosZ0[i];

                    // 全エネルギー・深さで加重合計
                    double sum = 0;

                    if (isHexGrid) // 260331Cl: 六方格子
                    {
                        int i3 = i * 3;
                        int hIdx0 = pIdx0[i3], hIdx1 = pIdx0[i3 + 1], hIdx2 = pIdx0[i3 + 2];
                        float hw0 = pWt0[i3], hw1 = pWt0[i3 + 1], hw2 = pWt0[i3 + 2];
                        for (int ei = 0; ei < eLen; ei++)
                            for (int di = 0; di < dLen; di++)
                            {
                                int wIdx = ei * dLen + di;
                                double weight = c00 * bw00[wIdx] + c10 * bw10[wIdx] + c01 * bw01[wIdx] + c11 * bw11[wIdx];
                                if (weight < 1e-15) continue;
                                var plane = posZ ? posPlanes[wIdx] : negPlanes[wIdx];//260718Cl 事前展開した配列を参照
                                if (plane == null || plane.Length == 0) continue;
                                sum += weight * (hw0 * plane[hIdx0] + hw1 * plane[hIdx1] + hw2 * plane[hIdx2]);
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
                            for (int di = 0; di < dLen; di++)
                            {
                                int wIdx = ei * dLen + di;
                                double weight = c00 * bw00[wIdx] + c10 * bw10[wIdx] + c01 * bw01[wIdx] + c11 * bw11[wIdx];
                                if (weight < 1e-15) continue;
                                var plane = posZ ? posPlanes[wIdx] : negPlanes[wIdx];//260718Cl 事前展開した配列を参照
                                if (plane == null || plane.Length == 0) continue;
                                double intensity = (mpW0 * plane[idx] + mpW1 * plane[idx + 1]) * mpFh1 + (mpW0 * plane[idx + gs] + mpW1 * plane[idx + gs + 1]) * mpFh;
                                sum += weight * intensity;
                            }
                    }
                    pVal0[i] = sum;
                }
            });
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
        double xm = view.XMirror; // 260718Cl: 左右反転 (UI スレッドで捕捉)
        int eLen = mp.Energies.Length, dLen = mp.Depths.Length;
        int binCount = dist.BinCount;
        var gs = lookupGridSize;
        double scaleW = view.ScaleW, scaleH = view.ScaleH, viewOffX = view.OffX, viewOffY = view.OffY; // 260724Cl 追加
        double halfW = view.HalfWidth, halfH = view.HalfHeight; // 260724Cl 追加
        EnsureGlobalNormalizationFactorsModel1(mp); //260726Cl: 呼び出し側の Ensure 忘れを構造的に不可能にする (係数はキャッシュ済みなら再計算しない)
        var planeScaleFactors = globalNormalizationFactors;
        var (posPlanes, negPlanes) = GetAllPlanes(mp, eLen, dLen);//260718Cl

        //Array.Clear(values); //260725Ch: 全画素上書きのため不要

        fixed (int* pIdx = lookupIdx)
        fixed (float* pWt = lookupWt)
        fixed (bool* pPosZ = lookupPosZ)
        fixed (double* pVal = values)
        {
            var pIdx0 = pIdx; var pWt0 = pWt; var pPosZ0 = pPosZ; var pVal0 = pVal;
            var isHexGrid = lookupGridType == MasterPattern.Types.Hexagonal; // 260331Cl

            Parallel.For(0, height, h =>
            {
                // double detNormY = (2.0 * h + 1 - height) / (double)height; // 260724Cl 変更前
                double detNormY = ((2.0 * h + 1 - height) * scaleH + viewOffY) / halfH; // 260724Cl: ラスター=視野全体化 (物理位置/halfH)
                double by = (1 - detNormY) * 0.5 * binCount - 0.5;
                int bj0 = Math.Clamp((int)Math.Floor(by), 0, binCount - 2);
                double fy = Math.Clamp(by - bj0, 0, 1);

                for (int w = 0; w < width; w++)
                {
                    int i = h * width + w;
                    // double detNormX = -xm * (2.0 * w + 1 - width) / (double)width; // 260718Cl: 左右反転 xm // 260724Cl 変更前
                    double detNormX = -xm * ((2.0 * w + 1 - width) * scaleW + viewOffX) / halfW; // 260724Cl
                    double bx = (detNormX + 1) * 0.5 * binCount - 0.5;
                    int bi0 = Math.Clamp((int)Math.Floor(bx), 0, binCount - 2);
                    double fx = Math.Clamp(bx - bi0, 0, 1);

                    double c00 = (1 - fx) * (1 - fy), c10 = fx * (1 - fy), c01 = (1 - fx) * fy, c11 = fx * fy;

                    double[] bw00 = dist.BinWeights[bi0, bj0], bw10 = dist.BinWeights[bi0 + 1, bj0], bw01 = dist.BinWeights[bi0, bj0 + 1], bw11 = dist.BinWeights[bi0 + 1, bj0 + 1];
                    bool posZ = pPosZ0[i];

                    double sum = 0;
                    if (isHexGrid) // 260331Cl
                    {
                        int i3 = i * 3;
                        int hIdx0 = pIdx0[i3], hIdx1 = pIdx0[i3 + 1], hIdx2 = pIdx0[i3 + 2];
                        float hw0 = pWt0[i3], hw1 = pWt0[i3 + 1], hw2 = pWt0[i3 + 2];
                        for (int ei = 0; ei < eLen; ei++)
                            for (int di = 0; di < dLen; di++)
                            {
                                int wIdx = ei * dLen + di;
                                double weight = c00 * bw00[wIdx] + c10 * bw10[wIdx] + c01 * bw01[wIdx] + c11 * bw11[wIdx];
                                if (weight < 1e-15) continue;
                                double planeScaleFactor = (uint)wIdx < (uint)planeScaleFactors.Length ? planeScaleFactors[wIdx] : 0.0;
                                if (planeScaleFactor < 1e-30) continue;
                                var plane = posZ ? posPlanes[wIdx] : negPlanes[wIdx];//260718Cl 事前展開した配列を参照
                                if (plane == null || plane.Length == 0) continue;
                                sum += weight * (hw0 * plane[hIdx0] + hw1 * plane[hIdx1] + hw2 * plane[hIdx2]) * planeScaleFactor;
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
                            for (int di = 0; di < dLen; di++)
                            {
                                int wIdx = ei * dLen + di;
                                double weight = c00 * bw00[wIdx] + c10 * bw10[wIdx]
                                              + c01 * bw01[wIdx] + c11 * bw11[wIdx];
                                if (weight < 1e-15) continue;
                                double planeScaleFactor = (uint)wIdx < (uint)planeScaleFactors.Length ? planeScaleFactors[wIdx] : 0.0;
                                if (planeScaleFactor < 1e-30) continue;
                                var plane = posZ ? posPlanes[wIdx] : negPlanes[wIdx];//260718Cl 事前展開した配列を参照
                                if (plane == null || plane.Length == 0) continue;
                                double intensity = (mpW0 * plane[idx] + mpW1 * plane[idx + 1]) * mpFh1
                                                 + (mpW0 * plane[idx + gs] + mpW1 * plane[idx + gs + 1]) * mpFh;
                                sum += weight * intensity * planeScaleFactor;
                            }
                    }
                    pVal0[i] = sum;
                }
            });
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
        double xm = view.XMirror; // 260718Cl: 左右反転 (UI スレッドで捕捉)
        int eLen = mp.Energies.Length, dLen = mp.Depths.Length;
        int binCount = dist.BinCount;
        var gs = lookupGridSize;
        double scaleW = view.ScaleW, scaleH = view.ScaleH, viewOffX = view.OffX, viewOffY = view.OffY; // 260724Cl 追加
        double halfW = view.HalfWidth, halfH = view.HalfHeight; // 260724Cl 追加
        var (posPlanes, negPlanes) = GetAllPlanes(mp, eLen, dLen);//260718Cl
        //260726Cl 追加 (正本 §1.4): plane は累積 M(t) なので隣接差は区間積分。区間平均 R̄=ΔM/Δt にするため区間幅で割る
        //(MC 側の重みは区間質量なので割らない)。等間隔グリッドでは全体が定数倍だが、不等間隔では区間ごとの重み比が変わる
        var depthWidths = mp.DepthIntervals;

        //Array.Clear(values); //260725Ch: 全画素上書きのため不要

        fixed (int* pIdx = lookupIdx)
        fixed (float* pWt = lookupWt)
        fixed (bool* pPosZ = lookupPosZ)
        fixed (double* pVal = values)
        {
            var pIdx0 = pIdx; var pWt0 = pWt; var pPosZ0 = pPosZ; var pVal0 = pVal;
            var isHexGrid = lookupGridType == MasterPattern.Types.Hexagonal; // 260331Cl

            Parallel.For(0, height, h =>
            {
                // double detNormY = (2.0 * h + 1 - height) / (double)height; // 260724Cl 変更前
                double detNormY = ((2.0 * h + 1 - height) * scaleH + viewOffY) / halfH; // 260724Cl: ラスター=視野全体化 (物理位置/halfH)
                double by = (1 - detNormY) * 0.5 * binCount - 0.5;
                int bj0 = Math.Clamp((int)Math.Floor(by), 0, binCount - 2);
                double fy = Math.Clamp(by - bj0, 0, 1);

                for (int w = 0; w < width; w++)
                {
                    int i = h * width + w;
                    // double detNormX = -xm * (2.0 * w + 1 - width) / (double)width; // 260718Cl: 左右反転 xm // 260724Cl 変更前
                    double detNormX = -xm * ((2.0 * w + 1 - width) * scaleW + viewOffX) / halfW; // 260724Cl
                    double bx = (detNormX + 1) * 0.5 * binCount - 0.5;
                    int bi0 = Math.Clamp((int)Math.Floor(bx), 0, binCount - 2);
                    double fx = Math.Clamp(bx - bi0, 0, 1);

                    double c00 = (1 - fx) * (1 - fy), c10 = fx * (1 - fy), c01 = (1 - fx) * fy, c11 = fx * fy;

                    double[] bw00 = dist.BinAbsoluteSliceWeights[bi0, bj0], bw10 = dist.BinAbsoluteSliceWeights[bi0 + 1, bj0], bw01 = dist.BinAbsoluteSliceWeights[bi0, bj0 + 1], bw11 = dist.BinAbsoluteSliceWeights[bi0 + 1, bj0 + 1];
                    bool posZ = pPosZ0[i];

                    double sum = 0;
                    if (isHexGrid) // 260331Cl
                    {
                        int i3 = i * 3;
                        int hIdx0 = pIdx0[i3], hIdx1 = pIdx0[i3 + 1], hIdx2 = pIdx0[i3 + 2];
                        float hw0 = pWt0[i3], hw1 = pWt0[i3 + 1], hw2 = pWt0[i3 + 2];
                        for (int ei = 0; ei < eLen; ei++)
                            for (int di = 0; di < dLen; di++)
                            {
                                int wIdx = ei * dLen + di;
                                double weight = c00 * bw00[wIdx] + c10 * bw10[wIdx] + c01 * bw01[wIdx] + c11 * bw11[wIdx];
                                if (weight < 1e-15) continue;
                                var plane = posZ ? posPlanes[wIdx] : negPlanes[wIdx];//260718Cl 事前展開した配列を参照
                                if (plane == null || plane.Length == 0) continue;
                                var planePrevious = di > 0 ? posZ ? posPlanes[wIdx - 1] : negPlanes[wIdx - 1] : null;//260718Cl 事前展開した配列を参照 (di>0 なら wIdx-1 = 同 energy の di-1)
                                double intensity = hw0 * plane[hIdx0] + hw1 * plane[hIdx1] + hw2 * plane[hIdx2];
                                if (planePrevious != null && planePrevious.Length > 0)
                                    intensity -= hw0 * planePrevious[hIdx0] + hw1 * planePrevious[hIdx1] + hw2 * planePrevious[hIdx2];
                                sum += weight * Math.Max(0.0, intensity) / depthWidths[di]; //260726Cl: 区間平均 ΔM/Δt
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
                            for (int di = 0; di < dLen; di++)
                            {
                                int wIdx = ei * dLen + di;
                                double weight = c00 * bw00[wIdx] + c10 * bw10[wIdx]
                                              + c01 * bw01[wIdx] + c11 * bw11[wIdx];
                                if (weight < 1e-15) continue;
                                var plane = posZ ? posPlanes[wIdx] : negPlanes[wIdx];//260718Cl 事前展開した配列を参照
                                if (plane == null || plane.Length == 0) continue;
                                var planePrevious = di > 0 ? posZ ? posPlanes[wIdx - 1] : negPlanes[wIdx - 1] : null;//260718Cl 事前展開した配列を参照 (di>0 なら wIdx-1 = 同 energy の di-1)
                                double intensity = (mpW0 * plane[idx] + mpW1 * plane[idx + 1]) * mpFh1
                                                 + (mpW0 * plane[idx + gs] + mpW1 * plane[idx + gs + 1]) * mpFh;
                                if (planePrevious != null && planePrevious.Length > 0)
                                    intensity -= (mpW0 * planePrevious[idx] + mpW1 * planePrevious[idx + 1]) * mpFh1
                                             + (mpW0 * planePrevious[idx + gs] + mpW1 * planePrevious[idx + gs + 1]) * mpFh;
                                sum += weight * Math.Max(0.0, intensity) / depthWidths[di]; //260726Cl: 区間平均 ΔM/Δt
                            }
                    }
                    pVal0[i] = sum;
                }
            });
        }
    }
    #endregion
}
