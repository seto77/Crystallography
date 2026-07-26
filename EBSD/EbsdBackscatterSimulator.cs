#region using
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using V3 = OpenTK.Mathematics.Vector3d;
using M3 = OpenTK.Mathematics.Matrix3d;
#endregion

namespace Crystallography;

/// <summary>
/// EBSD の統計解析に使う後方散乱電子 1 本ぶんの情報。260726Cl 追加:
/// <see cref="MonteCarlo.BackscatteredElectronDetail"/> に、ステレオ投影 (Schmidt) 済みの出射方向 <see cref="Position"/> を足したもの。
/// (FormEBSD が持っていた 10 要素タプルを名前付きにした。メンバー名はタプル要素名と同一)
/// </summary>
public readonly record struct EbsdBackscatteredElectron(
    double Depth, V3 Vec, PointD Position, double Energy, double TotalEnergyLoss,
    bool HasLastInelasticEvent, double LastInelasticDepth,
    double LastInelasticEnergyBeforeLoss, double LastInelasticEnergyAfterLoss, V3 LastInelasticDirection);

/// <summary>
/// モンテカルロによる後方散乱電子の飛程シミュレーション (大量の電子を並列に走らせ、脱出した電子を集める)。
/// 260726Cl 追加: FormEBSD.cs の RunBackscatterMonteCarlo をそのまま移設したもの (GUI 非依存)。
/// </summary>
public static class EbsdBackscatterSimulator
{
    /// <summary>モンテカルロによる飛程シミュレーション</summary>
    // private static ... RunBackscatterMonteCarlo(..., Action<int, int> reportProgress = null) // 260406Cl 旧シグネチャ: CancellationToken なし
    //260726Cl シグネチャ変更 (FormEBSD.RunBackscatterMonteCarlo から移設): 戻り値のタプル配列を EbsdBackscatteredElectron[] へ (メンバー名は同一)
    public static EbsdBackscatteredElectron[] Run(
        MonteCarlo monte, int loop, double energyThreshold, M3 sampleRotation,
        Action<int, int> reportProgress = null, CancellationToken cancellationToken = default) // 260406Cl CancellationToken 追加
    {
        ArgumentNullException.ThrowIfNull(monte);
        var bseLists = new ConcurrentBag<List<EbsdBackscatteredElectron>>(); // (260331Ch) 最後の非弾性散乱情報も保持する
        var reportStep = reportProgress == null ? int.MaxValue : Math.Max(1, loop / 100); // (260327Ch) UI へは 1% ごとにだけ流す
        int completed = 0;
        const int progressBatch = 1024; // 260603Cl 追加: 共有カウンタの Interlocked を worker ごとにまとめ、毎電子の cache-line 競合を解消する

        // 260603Cl 追加: バッチ加算後の進捗通知 (本体と localFinally の共通処理)。reportStep(~1%) 境界を跨いだ時だけ通知し、微妙な境界判定を一元化する
        void notifyProgress(int delta, int current)
        {
            if (reportProgress != null && (current >= loop || current / reportStep != (current - delta) / reportStep))
                reportProgress(Math.Min(current, loop), loop); // (260327Ch) MasterPattern 前処理時だけ進捗を通知する
        }

        Parallel.For(0, loop,
            new ParallelOptions { CancellationToken = cancellationToken }, // 260406Cl 追加: キャンセル対応
                                                                          // () => new List<(...)>(256), // 260603Cl 旧: thread-local は List のみ
            () => (list: new List<EbsdBackscatteredElectron>(256), pending: 0), // 260603Cl thread-local に進捗カウンタ pending を同梱
            (index, state, local) =>
            {
                var electron = monte.GetBackscatteredElectronDetail();
                if (electron.Energy > energyThreshold)
                    local.list.Add(new EbsdBackscatteredElectron(electron.Depth, electron.Direction, Stereonet.ConvertVectorToSchmidt(sampleRotation * electron.Direction), electron.Energy,
                        electron.TotalEnergyLoss, electron.HasLastInelasticEvent, electron.LastInelasticDepth, electron.LastInelasticEnergyBeforeLoss, electron.LastInelasticEnergyAfterLoss, electron.LastInelasticDirection)); // (260331Ch)

                // 260603Cl 旧: 毎電子で共有カウンタを Interlocked.Increment (全 worker が同一 cache-line を叩く)
                // var current = Interlocked.Increment(ref completed);
                // if (reportProgress != null && (current == loop || current % reportStep == 0))
                //     reportProgress(current, loop);
                if (++local.pending >= progressBatch) // 260603Cl progressBatch 電子貯めてから 1 回だけ共有カウンタへ加算
                {
                    int delta = local.pending;
                    local.pending = 0;
                    notifyProgress(delta, Interlocked.Add(ref completed, delta));
                }
                return local;
            },
            // localList => { if (localList.Count > 0) bseLists.Add(localList); }, // 260603Cl 旧
            local =>
            {
                if (local.pending > 0) // 260603Cl worker 終了時に端数電子を共有カウンタへ反映
                    notifyProgress(local.pending, Interlocked.Add(ref completed, local.pending));
                if (local.list.Count > 0)
                    bseLists.Add(local.list);
            });

        return [.. bseLists.SelectMany(localList => localList)];
    }
}
