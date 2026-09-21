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

    /// <summary>260727Cl 追加: weighted 合成の前提 — 分布の重み配列が mp と同じ (energy × depth) 格子であること — を検査する。
    /// 重み参照は wIdx = ei * dLen + di (dLen は mp 由来) を dist の配列へ素通しするので、格子がずれると
    /// dLen が減る方向では**例外も警告も出さずに別スライスの重み**を引き、増える方向では IndexOutOfRange になる。
    /// MC は MasterPattern の構築とは別のタイミングでも走る (Calc BSE・検出器幾何変更時の再ビニング) ので、
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
    private (MasterPattern Mp, EbsdMonteCarloDistribution Dist, bool Differential, double[] Pos, double[] Neg) planeMeanCache;
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

    /// <summary>260920Cl 追加 (作者指示: 表面非晶質層を入れると格子状の暗線が出る件): 検出器 8x8 ビンの双線形補間に使う
    /// 小数部を Smoothstep でならす。素の Clamp(f, 0, 1) だと補間は C0 止まりで、ビンの継ぎ目ごとに傾きが折れる。
    /// とくに最外の継ぎ目では、外側で f が 0 または 1 に凍る (= 傾き 0) のに内側は傾きを持つので折れ方が最大になる。
    /// 【なぜ見えるか】コヒーレント項では重み場の折れ目が菊池模様に紛れて見えない。ところが「方向依存を持たない項」は
    ///   重み場をそのまま像にする。具体的には非晶質層の配分 fA と、その相手の方向平均 sumMean。背景平坦化の高域通過を
    ///   通すと、傾きの折れが細い暗線になる (790x602 の実例で列 49/740・行 37/564 = 最外の継ぎ目に一致)。
    /// 【なぜ Smoothstep か】3f²−2f³ はビン中心で傾きが 0 になるので、外側の凍った領域と滑らかに繋がり、
    ///   内側の継ぎ目も C1 になる。ビン中心での値は変わらないので、フィットの節点は動かない。
    /// ⚠これは fA だけでなく重み全体の補間に効くので、コヒーレント項の数値もビン中心以外でわずかに変わる。
    ///   元の場は 8x8 の統計フィットで、piecewise-bilinear である必然性は無いため、滑らかな内挿の方が素直と判断した</summary>
    static double SmoothBinFraction(double f)
    {
        f = Math.Clamp(f, 0, 1);
        return f * f * (3 - 2 * f);
    }

    /// <summary>260920Cl 追加 (作者指示): A(E) で失ったコントラスト分を戻す台座。検出器位置にも方向にも依らない 1 個のスカラー。
    /// 【なぜスカラーか】最初は画素ごとに、その場の MC ビン補間重みで作っていた。コヒーレント項では重み場の粗さが
    ///   菊池模様に紛れて見えないが、台座は方向依存が無いぶん重み場をそのまま像にしてしまう。MC 重みは 8x8 の検出器ビンを
    ///   双線形補間した piecewise-bilinear な場なので継ぎ目で傾きが折れ、特に最外の継ぎ目 (bx/by が 0 と BinCount−1。
    ///   そこから外側は fx/fy がクランプされて凍る) で折れ方が最大になる。背景平坦化の高域通過を通すと、そこが
    ///   格子状の暗線として現れた (790x602 の実例で列 49/740・行 37/564 = 最外の継ぎ目に一致。E_c を小さくするほど顕著)。
    ///   台座の検出器面内の変化は「取り出し角による射出エネルギー分布の違い」という小さく滑らかな効果でしかなく、
    ///   8x8 のビン格子ではそもそも忠実に表せない。全ビン平均の重みで 1 個のスカラーにすれば、物理を落とさずに折れ目が消える。
    /// 【方向平均】半球ごとの平均をそのまま使うと赤道で段差が出るので、両半球の平均 = 全球平均を使う</summary>
    static double IncoherentPedestal(EbsdMonteCarloDistribution dist, double[] cohA, double[] posMeans, double[] negMeans,
        int eLen, int dLen, bool differential, double[] depthWidths, double[] planeScaleFactors)
    {
        if (cohA == null || posMeans == null || negMeans == null) return 0; //A(E) 無効 = 台座なし (従来動作)
        var g = new double[eLen * dLen];
        int nb = 0;
        for (int bi = 0; bi < dist.BinCount; bi++)
            for (int bj = 0; bj < dist.BinCount; bj++)
            {
                var bw = differential ? dist.BinAbsoluteSliceWeights[bi, bj] : dist.BinWeights[bi, bj];
                if (bw == null) continue;
                nb++;
                for (int k = 0; k < g.Length && k < bw.Length; k++) g[k] += bw[k];
            }
        if (nb == 0) return 0;
        double p = 0;
        for (int ei = 0; ei < eLen; ei++)
        {
            double mean = 0.5 * (posMeans[ei] + negMeans[ei]); //全球の方向平均
            double wSum = 0;
            for (int di = 0; di < dLen; di++)
            {
                int k = ei * dLen + di;
                double w = g[k] / nb;
                if (differential && depthWidths != null) w /= depthWidths[di]; //model 2 は区間平均 ΔM/Δt に合わせる
                if (planeScaleFactors != null) w *= (uint)k < (uint)planeScaleFactors.Length ? planeScaleFactors[k] : 0.0; //model 1 の規格化係数
                wSum += w;
            }
            p += (1 - cohA[ei]) * mean * wSum;
        }
        return p;
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

    private (double[] pos, double[] neg) GetPlaneMeansCached(MasterPattern mp, EbsdMonteCarloDistribution dist, float[][] posPlanes, float[][] negPlanes, int dLen, double[] depthWidths)
    {
        bool differential = depthWidths != null;
        if (ReferenceEquals(planeMeanCache.Mp, mp) && ReferenceEquals(planeMeanCache.Dist, dist) && planeMeanCache.Differential == differential && planeMeanCache.Pos != null)
            return (planeMeanCache.Pos, planeMeanCache.Neg);
        var r = GetPlaneMeans(posPlanes, negPlanes, dLen, depthWidths);
        var g = differential ? dist.GlobalDepthSliceWeights : dist.GlobalDepthWeights;
        var pos = CollapseToEnergyReference(r.pos, dLen, g);
        var neg = CollapseToEnergyReference(r.neg, dLen, g);
        planeMeanCache = (mp, dist, differential, pos, neg);
        return (pos, neg);
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
        double halfW = view.HalfWidth, halfH = view.HalfHeight; // 260724Cl 追加
        var (posPlanes, negPlanes) = GetAllPlanes(mp, eLen, dLen);//260718Cl
        var amorphousFraction = dist.BinAmorphousFraction; // 260919Cl 追加: 表面非晶質層に源を持つ電子の割合 (ビンごと)
        bool hasAmorphous = dist.HasAmorphousLayer; // 260919Cl 追加 (/simplify: 層が無ければ fA の双線形補間も省く)
        //260920Cl 追加: 損失依存のコントラスト係数 A(E)。⚠非晶質層が無いと dist.GlobalDepthWeights は null なので、
        //  台座の基準値 (方向平均) は深さの単純平均になる (CollapseToEnergyReference のフォールバック)。方向に依らない成分なので
        //  バンド幅にも ZNCC にも効かないが、絶対値を論じるときはここが MC 重みで積まれていないことに注意
        var cohA = BuildCoherenceFactors(mp);
        var (posMeans, negMeans) = hasAmorphous || cohA != null ? GetPlaneMeansCached(mp, dist, posPlanes, negPlanes, dLen, null) : (null, null); // 260919Cl 追加: 変調なし成分用の方向平均 / 260920Cl A(E) でも使う
        //260920Cl 追加: 変調なし成分の基準は**全球**の方向平均。半球ごとの平均 (posMeans / negMeans) をそのまま使うと、
        //  パターンが赤道 (試料系 z = 0、ノモニック投影では直線) を跨ぐ所で段差になる。源の向きを失った電子に半球の区別は無い
        double[] sphereMeans = null;
        if (posMeans != null) { sphereMeans = new double[posMeans.Length]; for (int q = 0; q < sphereMeans.Length; q++) sphereMeans[q] = 0.5 * (posMeans[q] + negMeans[q]); }
        double incoherentPedestal = IncoherentPedestal(dist, cohA, posMeans, negMeans, eLen, dLen, false, null, null); //260920Cl 追加

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
                //260920Cl 変更: 素の双線形だと継ぎ目で傾きが折れる。Smoothstep でならす (理由は SmoothBinFraction の doc)
                //旧: double fy = Math.Clamp(by - bj0, 0, 1);
                double fy = SmoothBinFraction(by - bj0);

                for (int w = 0; w < width; w++)
                {
                    int i = h * width + w;

                    // 検出器 X 座標
                    // double detNormX = -xm * (2.0 * w + 1 - width) / (double)width; // 260325Cl: スクリーン X は検出器面 X と反転 (BuildLookupTable で -Ri.E11 を使用) / 260718Cl: 左右反転 xm を掛ける // 260724Cl 変更前
                    double detNormX = -xm * ((2.0 * w + 1 - width) * scaleW + viewOffX) / halfW; // 260724Cl: 物理位置/halfW
                    double bx = (detNormX + 1) * 0.5 * binCount - 0.5;
                    int bi0 = Math.Clamp((int)Math.Floor(bx), 0, binCount - 2);
                    //旧: double fx = Math.Clamp(bx - bi0, 0, 1); //260920Cl 変更: fy と同じく Smoothstep
                    double fx = SmoothBinFraction(bx - bi0);

                    // ビン重みのバイリニア補間係数
                    double c00 = (1 - fx) * (1 - fy), c10 = fx * (1 - fy), c01 = (1 - fx) * fy, c11 = fx * fy;
                    double fA = hasAmorphous ? c00 * amorphousFraction[bi0, bj0] + c10 * amorphousFraction[bi0 + 1, bj0] + c01 * amorphousFraction[bi0, bj0 + 1] + c11 * amorphousFraction[bi0 + 1, bj0 + 1] : 0; // 260919Cl 追加: 非晶質源の割合 (双線形)

                    var bw00 = dist.BinWeights[bi0, bj0];
                    var bw10 = dist.BinWeights[bi0 + 1, bj0];
                    var bw01 = dist.BinWeights[bi0, bj0 + 1];
                    var bw11 = dist.BinWeights[bi0 + 1, bj0 + 1];

                    // ルックアップテーブルからマスターパターン補間パラメータ取得
                    bool posZ = pPosZ0[i];

                    // 全エネルギー・深さで加重合計
                    double sum = 0;
                    double sumMean = 0; // 260919Cl 追加: 非晶質源 (変調なし) 用の方向平均強度

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
                                double weight = c00 * bw00[wIdx] + c10 * bw10[wIdx] + c01 * bw01[wIdx] + c11 * bw11[wIdx];
                                if (weight < 1e-15) continue;
                                var plane = posZ ? posPlanes[wIdx] : negPlanes[wIdx];//260718Cl 事前展開した配列を参照
                                if (plane == null || plane.Length == 0) continue;
                                sum += weight * (hw0 * plane[hIdx0] + hw1 * plane[hIdx1] + hw2 * plane[hIdx2]) * aE; //260920Cl 変更: A(E) を末尾に掛ける (無効時は 1.0 なので丸めも含めて従来と同一) 
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
                                double weight = c00 * bw00[wIdx] + c10 * bw10[wIdx] + c01 * bw01[wIdx] + c11 * bw11[wIdx];
                                if (weight < 1e-15) continue;
                                var plane = posZ ? posPlanes[wIdx] : negPlanes[wIdx];//260718Cl 事前展開した配列を参照
                                if (plane == null || plane.Length == 0) continue;
                                double intensity = (mpW0 * plane[idx] + mpW1 * plane[idx + 1]) * mpFh1 + (mpW0 * plane[idx + gs] + mpW1 * plane[idx + gs + 1]) * mpFh;
                                sum += weight * intensity * aE; //260920Cl 変更: A(E) を末尾に掛ける (無効時は 1.0 なので丸めも含めて従来と同一) 
                                if (fA > 0) sumMean += weight * sMean; // 260919Cl 追加 / 260920Cl 変更: 半球平均 → 全球平均
                            }
                        }
                    }
                    // pVal0[i] = sum; // 260919Cl 変更前
                    //260920Cl 変更: A(E) で失った分を台座として足し戻す (総量保存)。A(E) 無効時は 0 なので従来と数値が一致する。
                    //  ⚠台座は検出器位置に依らない 1 個のスカラーにしてある。理由は IncoherentPedestal の doc を参照
                    double coherentSum = sum + incoherentPedestal; //260920Cl
                    pVal0[i] = fA > 0 ? (1 - fA) * coherentSum + fA * sumMean : coherentSum; // 260919Cl 変更: 非晶質層内の源は方向平均 (変調なし) で寄与
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
        EnsureGridMatches(mp, dist); //260727Cl
        double xm = view.XMirror; // 260718Cl: 左右反転 (UI スレッドで捕捉)
        int eLen = mp.Energies.Length, dLen = mp.Depths.Length;
        int binCount = dist.BinCount;
        var gs = lookupGridSize;
        double scaleW = view.ScaleW, scaleH = view.ScaleH, viewOffX = view.OffX, viewOffY = view.OffY; // 260724Cl 追加
        double halfW = view.HalfWidth, halfH = view.HalfHeight; // 260724Cl 追加
        EnsureGlobalNormalizationFactorsModel1(mp); //260726Cl: 呼び出し側の Ensure 忘れを構造的に不可能にする (係数はキャッシュ済みなら再計算しない)
        var planeScaleFactors = globalNormalizationFactors;
        var (posPlanes, negPlanes) = GetAllPlanes(mp, eLen, dLen);//260718Cl
        var amorphousFraction = dist.BinAmorphousFraction; // 260919Cl 追加: 表面非晶質層に源を持つ電子の割合 (ビンごと)
        bool hasAmorphous = dist.HasAmorphousLayer; // 260919Cl 追加 (/simplify: 層が無ければ fA の双線形補間も省く)
        //260920Cl 追加: 損失依存のコントラスト係数 A(E)。⚠非晶質層が無いと dist.GlobalDepthWeights は null なので、
        //  台座の基準値 (方向平均) は深さの単純平均になる (CollapseToEnergyReference のフォールバック)。方向に依らない成分なので
        //  バンド幅にも ZNCC にも効かないが、絶対値を論じるときはここが MC 重みで積まれていないことに注意
        var cohA = BuildCoherenceFactors(mp);
        var (posMeans, negMeans) = hasAmorphous || cohA != null ? GetPlaneMeansCached(mp, dist, posPlanes, negPlanes, dLen, null) : (null, null); // 260919Cl 追加: 変調なし成分用の方向平均 / 260920Cl A(E) でも使う
        //260920Cl 追加: 変調なし成分の基準は**全球**の方向平均。半球ごとの平均 (posMeans / negMeans) をそのまま使うと、
        //  パターンが赤道 (試料系 z = 0、ノモニック投影では直線) を跨ぐ所で段差になる。源の向きを失った電子に半球の区別は無い
        double[] sphereMeans = null;
        if (posMeans != null) { sphereMeans = new double[posMeans.Length]; for (int q = 0; q < sphereMeans.Length; q++) sphereMeans[q] = 0.5 * (posMeans[q] + negMeans[q]); }
        double incoherentPedestal = IncoherentPedestal(dist, cohA, posMeans, negMeans, eLen, dLen, false, null, planeScaleFactors); //260920Cl 追加

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
                //260920Cl 変更: 素の双線形だと継ぎ目で傾きが折れる。Smoothstep でならす (理由は SmoothBinFraction の doc)
                //旧: double fy = Math.Clamp(by - bj0, 0, 1);
                double fy = SmoothBinFraction(by - bj0);

                for (int w = 0; w < width; w++)
                {
                    int i = h * width + w;
                    // double detNormX = -xm * (2.0 * w + 1 - width) / (double)width; // 260718Cl: 左右反転 xm // 260724Cl 変更前
                    double detNormX = -xm * ((2.0 * w + 1 - width) * scaleW + viewOffX) / halfW; // 260724Cl
                    double bx = (detNormX + 1) * 0.5 * binCount - 0.5;
                    int bi0 = Math.Clamp((int)Math.Floor(bx), 0, binCount - 2);
                    //旧: double fx = Math.Clamp(bx - bi0, 0, 1); //260920Cl 変更: fy と同じく Smoothstep
                    double fx = SmoothBinFraction(bx - bi0);

                    double c00 = (1 - fx) * (1 - fy), c10 = fx * (1 - fy), c01 = (1 - fx) * fy, c11 = fx * fy;
                    double fA = hasAmorphous ? c00 * amorphousFraction[bi0, bj0] + c10 * amorphousFraction[bi0 + 1, bj0] + c01 * amorphousFraction[bi0, bj0 + 1] + c11 * amorphousFraction[bi0 + 1, bj0 + 1] : 0; // 260919Cl 追加: 非晶質源の割合 (双線形)

                    double[] bw00 = dist.BinWeights[bi0, bj0], bw10 = dist.BinWeights[bi0 + 1, bj0], bw01 = dist.BinWeights[bi0, bj0 + 1], bw11 = dist.BinWeights[bi0 + 1, bj0 + 1];
                    bool posZ = pPosZ0[i];

                    double sum = 0;
                    double sumMean = 0; // 260919Cl 追加: 非晶質源 (変調なし) 用の方向平均強度
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
                                double weight = c00 * bw00[wIdx] + c10 * bw10[wIdx] + c01 * bw01[wIdx] + c11 * bw11[wIdx];
                                if (weight < 1e-15) continue;
                                double planeScaleFactor = (uint)wIdx < (uint)planeScaleFactors.Length ? planeScaleFactors[wIdx] : 0.0;
                                if (planeScaleFactor < 1e-30) continue;
                                var plane = posZ ? posPlanes[wIdx] : negPlanes[wIdx];//260718Cl 事前展開した配列を参照
                                if (plane == null || plane.Length == 0) continue;
                                sum += weight * (hw0 * plane[hIdx0] + hw1 * plane[hIdx1] + hw2 * plane[hIdx2]) * planeScaleFactor * aE; //260920Cl 変更: A(E) を末尾に掛ける (無効時は 1.0 なので丸めも含めて従来と同一) 
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
                                double weight = c00 * bw00[wIdx] + c10 * bw10[wIdx]
                                              + c01 * bw01[wIdx] + c11 * bw11[wIdx];
                                if (weight < 1e-15) continue;
                                double planeScaleFactor = (uint)wIdx < (uint)planeScaleFactors.Length ? planeScaleFactors[wIdx] : 0.0;
                                if (planeScaleFactor < 1e-30) continue;
                                var plane = posZ ? posPlanes[wIdx] : negPlanes[wIdx];//260718Cl 事前展開した配列を参照
                                if (plane == null || plane.Length == 0) continue;
                                double intensity = (mpW0 * plane[idx] + mpW1 * plane[idx + 1]) * mpFh1
                                                 + (mpW0 * plane[idx + gs] + mpW1 * plane[idx + gs + 1]) * mpFh;
                                sum += weight * intensity * planeScaleFactor * aE; //260920Cl 変更: A(E) を末尾に掛ける (無効時は 1.0 なので丸めも含めて従来と同一) 
                                if (fA > 0) sumMean += weight * sMean * planeScaleFactor; // 260919Cl 追加 / 260920Cl 変更: 半球平均 → 全球平均
                            }
                        }
                    }
                    // pVal0[i] = sum; // 260919Cl 変更前
                    //260920Cl 変更: A(E) で失った分を台座として足し戻す (総量保存)。A(E) 無効時は 0 なので従来と数値が一致する。
                    //  ⚠台座は検出器位置に依らない 1 個のスカラーにしてある。理由は IncoherentPedestal の doc を参照
                    double coherentSum = sum + incoherentPedestal; //260920Cl
                    pVal0[i] = fA > 0 ? (1 - fA) * coherentSum + fA * sumMean : coherentSum; // 260919Cl 変更: 非晶質層内の源は方向平均 (変調なし) で寄与
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
        EnsureGridMatches(mp, dist); //260727Cl
        double xm = view.XMirror; // 260718Cl: 左右反転 (UI スレッドで捕捉)
        int eLen = mp.Energies.Length, dLen = mp.Depths.Length;
        int binCount = dist.BinCount;
        var gs = lookupGridSize;
        double scaleW = view.ScaleW, scaleH = view.ScaleH, viewOffX = view.OffX, viewOffY = view.OffY; // 260724Cl 追加
        double halfW = view.HalfWidth, halfH = view.HalfHeight; // 260724Cl 追加
        var (posPlanes, negPlanes) = GetAllPlanes(mp, eLen, dLen);//260718Cl
        var amorphousFraction = dist.BinAmorphousFraction; // 260919Cl 追加: 表面非晶質層に源を持つ電子の割合 (ビンごと)
        bool hasAmorphous = dist.HasAmorphousLayer; // 260919Cl 追加 (/simplify: 層が無ければ fA の双線形補間も省く)
        double[] posMeans = null, negMeans = null; // 260919Cl 追加: model 2 は差分 ΔM/Δt の平均 (depthWidths 確定後に取得)
        //260726Cl 追加 (正本 §1.4): plane は累積 M(t) なので隣接差は区間積分。区間平均 R̄=ΔM/Δt にするため区間幅で割る
        //(MC 側の重みは区間質量なので割らない)。等間隔グリッドでは全体が定数倍だが、不等間隔では区間ごとの重み比が変わる
        var depthWidths = mp.DepthIntervals;
        //260920Cl 追加: 損失依存のコントラスト係数 A(E)。⚠非晶質層が無いと dist.GlobalDepthWeights は null なので、
        //  台座の基準値 (方向平均) は深さの単純平均になる (CollapseToEnergyReference のフォールバック)。方向に依らない成分なので
        //  バンド幅にも ZNCC にも効かないが、絶対値を論じるときはここが MC 重みで積まれていないことに注意
        var cohA = BuildCoherenceFactors(mp);
        if (hasAmorphous || cohA != null) (posMeans, negMeans) = GetPlaneMeansCached(mp, dist, posPlanes, negPlanes, dLen, depthWidths); // 260919Cl 追加: model 2 は差分 ΔM/Δt の平均 (/simplify: 以前は null 版を先に呼んでキャッシュを取りこぼしていた) / 260920Cl A(E) でも使う
        //260920Cl 追加: 変調なし成分の基準は**全球**の方向平均。半球ごとの平均 (posMeans / negMeans) をそのまま使うと、
        //  パターンが赤道 (試料系 z = 0、ノモニック投影では直線) を跨ぐ所で段差になる。源の向きを失った電子に半球の区別は無い
        double[] sphereMeans = null;
        if (posMeans != null) { sphereMeans = new double[posMeans.Length]; for (int q = 0; q < sphereMeans.Length; q++) sphereMeans[q] = 0.5 * (posMeans[q] + negMeans[q]); }
        double incoherentPedestal = IncoherentPedestal(dist, cohA, posMeans, negMeans, eLen, dLen, true, depthWidths, null); //260920Cl 追加

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
                //260920Cl 変更: 素の双線形だと継ぎ目で傾きが折れる。Smoothstep でならす (理由は SmoothBinFraction の doc)
                //旧: double fy = Math.Clamp(by - bj0, 0, 1);
                double fy = SmoothBinFraction(by - bj0);

                for (int w = 0; w < width; w++)
                {
                    int i = h * width + w;
                    // double detNormX = -xm * (2.0 * w + 1 - width) / (double)width; // 260718Cl: 左右反転 xm // 260724Cl 変更前
                    double detNormX = -xm * ((2.0 * w + 1 - width) * scaleW + viewOffX) / halfW; // 260724Cl
                    double bx = (detNormX + 1) * 0.5 * binCount - 0.5;
                    int bi0 = Math.Clamp((int)Math.Floor(bx), 0, binCount - 2);
                    //旧: double fx = Math.Clamp(bx - bi0, 0, 1); //260920Cl 変更: fy と同じく Smoothstep
                    double fx = SmoothBinFraction(bx - bi0);

                    double c00 = (1 - fx) * (1 - fy), c10 = fx * (1 - fy), c01 = (1 - fx) * fy, c11 = fx * fy;
                    double fA = hasAmorphous ? c00 * amorphousFraction[bi0, bj0] + c10 * amorphousFraction[bi0 + 1, bj0] + c01 * amorphousFraction[bi0, bj0 + 1] + c11 * amorphousFraction[bi0 + 1, bj0 + 1] : 0; // 260919Cl 追加: 非晶質源の割合 (双線形)

                    double[] bw00 = dist.BinAbsoluteSliceWeights[bi0, bj0], bw10 = dist.BinAbsoluteSliceWeights[bi0 + 1, bj0], bw01 = dist.BinAbsoluteSliceWeights[bi0, bj0 + 1], bw11 = dist.BinAbsoluteSliceWeights[bi0 + 1, bj0 + 1];
                    bool posZ = pPosZ0[i];

                    double sum = 0;
                    double sumMean = 0; // 260919Cl 追加: 非晶質源 (変調なし) 用の方向平均強度
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
                                double weight = c00 * bw00[wIdx] + c10 * bw10[wIdx] + c01 * bw01[wIdx] + c11 * bw11[wIdx];
                                if (weight < 1e-15) continue;
                                var plane = posZ ? posPlanes[wIdx] : negPlanes[wIdx];//260718Cl 事前展開した配列を参照
                                if (plane == null || plane.Length == 0) continue;
                                var planePrevious = di > 0 ? posZ ? posPlanes[wIdx - 1] : negPlanes[wIdx - 1] : null;//260718Cl 事前展開した配列を参照 (di>0 なら wIdx-1 = 同 energy の di-1)
                                double intensity = hw0 * plane[hIdx0] + hw1 * plane[hIdx1] + hw2 * plane[hIdx2];
                                if (planePrevious != null && planePrevious.Length > 0)
                                    intensity -= hw0 * planePrevious[hIdx0] + hw1 * planePrevious[hIdx1] + hw2 * planePrevious[hIdx2];
                                sum += weight * Math.Max(0.0, intensity) / depthWidths[di] * aE; //260920Cl 変更: A(E) を末尾に掛ける (無効時は 1.0 なので丸めも含めて従来と同一) //260726Cl: 区間平均 ΔM/Δt
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
                                sum += weight * Math.Max(0.0, intensity) / depthWidths[di] * aE; //260920Cl 変更: A(E) を末尾に掛ける (無効時は 1.0 なので丸めも含めて従来と同一) //260726Cl: 区間平均 ΔM/Δt
                                if (fA > 0) sumMean += weight * sMean; // 260919Cl 追加 / 260920Cl 変更: 半球平均 → 全球平均 (平均は既に /Δt 済み)
                            }
                        }
                    }
                    // pVal0[i] = sum; // 260919Cl 変更前
                    //260920Cl 変更: A(E) で失った分を台座として足し戻す (総量保存)。A(E) 無効時は 0 なので従来と数値が一致する。
                    //  ⚠台座は検出器位置に依らない 1 個のスカラーにしてある。理由は IncoherentPedestal の doc を参照
                    double coherentSum = sum + incoherentPedestal; //260920Cl
                    pVal0[i] = fA > 0 ? (1 - fA) * coherentSum + fA * sumMean : coherentSum; // 260919Cl 変更: 非晶質層内の源は方向平均 (変調なし) で寄与
                }
            });
        }
    }
    #endregion
}
