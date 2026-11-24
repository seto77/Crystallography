#region using
using System;
using System.Collections.Generic; // 261031Cl 追加 (ReadElectronChunks の IEnumerable)
using System.Globalization;
using System.IO;
using System.IO.Compression;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Text.Json;
using System.Text.Json.Serialization; // 261028Cl 追加 (JsonNumberHandling)
#endregion

namespace Crystallography;

/// <summary>
/// 260921Cl 追加: 保存した MasterPattern の作成条件 (ファイルのヘッダに JSON で入る)。
/// 読み込み側はこれを今の結晶・条件と照合して、別の結晶や別の電圧で作ったパターンを黙って使わないようにする。
/// </summary>
public sealed class MasterPatternFileInfo
{
    /// <summary>ファイル形式の版。<see cref="MasterPatternFile.FormatVersion"/> と同じ値が入る</summary>
    public int FormatVersion { get; set; } = MasterPatternFile.FormatVersion;
    /// <summary>保存した日時 (UTC, ISO 8601)</summary>
    public string CreatedUtc { get; set; } = "";
    /// <summary>保存したアプリと版 (例 "ReciPro 4.947")</summary>
    public string Creator { get; set; } = "";
    /// <summary>結晶名 (表示用。照合には <see cref="CrystalFingerprint"/> を使う)</summary>
    public string CrystalName { get; set; } = "";
    /// <summary>結晶の同一性を表す文字列 (格子定数・空間群・原子の Z/座標/占有率/B)。<see cref="MasterPatternFile.CrystalFingerprint"/></summary>
    public string CrystalFingerprint { get; set; } = "";
    /// <summary>入射電子のエネルギー [keV]</summary>
    public double BeamEnergyKeV { get; set; }
    /// <summary>Bethe 計算のブロッホ波の最大数</summary>
    public int MaxNumOfBloch { get; set; }
    /// <summary>非局所吸収</summary>
    public bool UseNonLocalAbsorption { get; set; }
    /// <summary>非局所源 (TDS 背景)</summary>
    public bool IncludeTDSBackground { get; set; }
    /// <summary>吸収フラックスの再注入</summary>
    public bool IncludeAbsorbedFluxBackground { get; set; }
    /// <summary>自由記述 (作成条件の補足など)</summary>
    public string Note { get; set; } = "";
    /// <summary>保存した後方散乱電子の数 (0 なら電子の区画なし)</summary>
    public int BackscatteredElectronCount { get; set; }
    /// <summary>電子を作った MC の試料傾斜 [rad]。読み込み側はこれが今の試料傾斜と一致するときだけ電子を使い回す</summary>
    public double McSampleTiltRad { get; set; } = double.NaN;
    /// <summary>電子を作った MC の追跡打ち切りエネルギー [keV]</summary>
    public double McEnergyThresholdKeV { get; set; } = double.NaN;
    /// <summary>電子を作った MC の入射電子数</summary>
    public int McIncidentCount { get; set; }
    /// <summary>261110Cl 追加: ビーム選択の方式 (要求・実際・新方式の設定の版・従来へ戻した理由)。項目の無い古いファイルは従来 (Legacy)</summary>
    public string BeamSelectionRequested { get; set; } = "Legacy";
    public string BeamSelectionUsed { get; set; } = "Legacy";
    public string BeamSelectionSettingsVersion { get; set; } = "";
    public string BeamSelectionFallbackReason { get; set; } = "";
    public string BeamSelectionNative { get; set; } = ""; // 261110Cl 追加 (codex2 第 32 回): 新方式が使った native の版
    /// <summary>261118Cl 追加 (codex2 第 5 回): 後方 TDS 源の軸の方式。項目の無い古いファイルは "" (従来 = 結晶の z で作ったもの)</summary>
    public string BackSourceAxis { get; set; } = "";
    /// <summary>261118Cl 追加: 従来のビーム選択の候補の打ち切り方 (Legacy・Shell・ResearchHook・Diagnostic xM)。着衣の裾は "" (対象外)。
    /// 項目の無い古いファイルも "" = 未記録 (当時の既定は Legacy だが、Shell を選んで作っていても区別できない。codex2 第 15 回)。空のまま保持する</summary>
    public string BeamCandidateCut { get; set; } = "";

    /// <summary>261124Cl 追加 (形式 v2): マスターパターンを計算した物理の版 (<see cref="MasterPatternFile.PhysicsVersion"/>)。
    /// 項目の無い古いファイル (v1) は "" = 未記録。読み込み側は今の版と違えば知らせる (作者決定 11-24: 警告して読む)</summary>
    public string PhysicsVersion { get; set; } = "";
    /// <summary>261124Cl 追加 (形式 v2): 作った結晶 (<see cref="Crystal2"/> の MemoryPack + Brotli を base64)。
    /// 今の結晶と合わないとき、この結晶を一覧に足して読み込めるようにする (作者決定 11-24)。古いファイルは ""</summary>
    public string CrystalData { get; set; } = "";
    /// <summary>261124Cl 追加 (形式 v2): 格子 (本体を読まずに条件を見せるため見出しにも入れる。v1 のファイルでは 0・空)</summary>
    public int GridSize { get; set; }
    public int GridType { get; set; }
    public double[] Energies { get; set; } = [];
    public double[] Depths { get; set; } = [];

    /// <summary>261124Cl 追加 (codex2 査読 中 4): 浅い写し。保存は見出しを書き換えるので、画面が持つ記録を写しで保存し、失敗しても元を変えない</summary>
    public MasterPatternFileInfo Clone() => (MasterPatternFileInfo)MemberwiseClone();

    /// <summary>261124Cl 追加: 同梱の結晶 (<see cref="CrystalData"/>) を復元する。無い・読めないときは null</summary>
    public Crystal GetCrystal()
    {
        if (string.IsNullOrEmpty(CrystalData)) return null;
        try { return Crystal2.Deserialize(Convert.FromBase64String(CrystalData))?.ToCrystal(); }
        // catch (FormatException) { return null; } // 261124Cl 変更前
        catch (Exception) { return null; } // 261124Cl 変更 (codex2 査読 中 6): 不完全な Crystal2 は ToCrystal() でほかの例外も出す。どれも「戻せない」として扱う
    }

    /// <summary>構築リクエストから作成条件を作る</summary>
    // public static MasterPatternFileInfo FromBuildRequest(EBSD.MasterPatternBuildRequest request, double beamEnergyKeV, string creator = "") // 261110Cl 旧シグネチャ
    public static MasterPatternFileInfo FromBuildRequest(EBSD.MasterPatternBuildRequest request, double beamEnergyKeV, string creator = "",
        MasterPattern result = null) // 261110Cl result 追加: 構築の結果 (実際のビーム選択の方式を写す。null なら要求だけ)
    {
        ArgumentNullException.ThrowIfNull(request);
        return new MasterPatternFileInfo
        {
            CreatedUtc = DateTime.UtcNow.ToString("o", CultureInfo.InvariantCulture),
            Creator = creator ?? "",
            CrystalName = request.Crystal.Name ?? "",
            CrystalFingerprint = MasterPatternFile.CrystalFingerprint(request.Crystal),
            BeamEnergyKeV = beamEnergyKeV,
            MaxNumOfBloch = request.MaxNumOfBloch,
            UseNonLocalAbsorption = request.UseNonLocalAbsorption,
            IncludeTDSBackground = request.IncludeTDSBackground,
            IncludeAbsorbedFluxBackground = request.IncludeAbsorbedFluxBackground,
            BeamSelectionRequested = request.BeamSelection.ToString(), // 261110Cl 追加
            BeamSelectionUsed = result?.BeamSelectionUsed ?? "Legacy",
            BeamSelectionSettingsVersion = result?.BeamSelectionSettingsVersion ?? "",
            BeamSelectionFallbackReason = result?.BeamSelectionFallbackReason ?? "",
            BeamSelectionNative = result?.BeamSelectionNative ?? "",
            BackSourceAxis = result?.BackSourceAxis ?? "", // 261118Cl 追加 (codex2 第 5 回)
            BeamCandidateCut = result?.BeamCandidateCut ?? "", // 261118Cl 追加
            PhysicsVersion = MasterPatternFile.PhysicsVersion, // 261124Cl 追加: 今この版の物理で計算した
            CrystalData = Convert.ToBase64String(Crystal2.Serialize(request.Crystal.ToCrystal2())), // 261124Cl 追加: 作った結晶を同梱
        };
    }
}

/// <summary>260921Cl 追加: <see cref="MasterPatternFile.Load"/> の結果。BackscatteredElectrons は保存されていなければ空</summary>
public sealed record MasterPatternFileContent(MasterPattern MasterPattern, MasterPatternFileInfo Info, EbsdBackscatteredElectron[] BackscatteredElectrons);

/// <summary>
/// 260921Cl 追加: <see cref="MasterPattern"/> のファイル保存・読み込み。
/// <para>grid 512 の構築は結晶によって数分かかり、しかも MC の乱数で (energy × depth) 格子が毎回わずかに変わる。
/// 同じ条件で方位探索や表示を何度も試すときに、構築をやり直さずに済むようにする。</para>
/// <para>【形式】先頭 "RPMP" + 版 (int32) + ヘッダ長 (int32) + ヘッダ (UTF-8 JSON = <see cref="MasterPatternFileInfo"/>) +
/// 本体 (zlib 圧縮)。本体は grid・格子の種類・エネルギー列・深さ列と、平面ごとに「長さ (int32、0 は欠落) + float 列」を
/// +Z 半球 → −Z 半球の順に energy-major で並べる。数値はリトルエンディアン。</para>
/// <para>任意で、MC が出した後方散乱電子 (<see cref="EbsdBackscatteredElectron"/>) を平面のあとに float32 で入れられる。
/// MC 分布 (<see cref="EbsdMonteCarloDistribution"/>) そのものは保存しない: 源の深さモード・非晶質層・蛍光体重み・格子を
/// 変えると電子からの再ビニングが要るので、電子を持っておけば MC をやり直さずに済み、しかも MC の乱数による揺らぎも入らない
/// (同じファイルを読めば毎回まったく同じ分布になる。条件を比べる調整作業ではこれが大事)。</para>
/// <para>261124Cl 追加 (作者決定 11-24): 【形式 v2】"RPMP" + 版 (int32 = 2) + ヘッダ長 (int32) + ヘッダ (JSON。v2 では格子・エネルギー列・深さ列・
/// 物理の版・同梱の結晶も入る) + 区画の表 (区画の数 int32 + 区画ごとに 種類 int32・番号 int32・位置 int64・格納長 int64・元の長さ int64・方式 int32・CRC32 uint32)
/// + 区画。平面は 1 枚ずつ 1 区画 (種類 1 = +Z、2 = −Z、番号 = energy-major の平面番号、元の長さ 0 は欠落)、電子は <see cref="ElectronsPerBlock"/> 本ずつ
/// 列ごと (float32 × 15 列 + フラグ 1 列) に並べて 1 区画 (種類 3、番号 = 塊の番号)。方式 1 = 4 バイトの並べ替え (Blosc の shuffle) + Brotli、
/// 0 = そのまま。CRC32 は元のバイト列のもの。区画ごとに独立なので並列に縮め・戻せ、一部だけも読める。
/// 実測 (grid 512 の Si・Al₂O₃): 平面は v1 の zlib で 0.62〜0.66、v2 で 0.29〜0.43 (MasterPatternファイル仕様_提案_261124.md)。v1 のファイルもそのまま読める。</para>
/// </summary>
public static class MasterPatternFile
{
    /// <summary>ファイル形式の版。形式を変えたら上げ、<see cref="Load"/> で旧版を読めるようにすること</summary>
    /// <remarks>261031Cl 追加: 本体 (格子の見出し・平面の並び・電子の区画) は <see cref="ReadElectronChunks"/> でも読む。形式を変えたら両方を直すこと</remarks>
    // public const int FormatVersion = 1; // 261124Cl 変更前
    public const int FormatVersion = 2; // 261124Cl 変更: 区画の形式 (v2)。v1 も読む

    /// <summary>261124Cl 追加 (作者決定 11-24): マスターパターンを計算する物理の版。同じ入力でマスターパターンの強度が変わる直し
    /// (源・吸収・ビーム選択の既定・散乱因子の表など) を入れたら上げる。読み込み側は、ファイルの版が今と違えば警告して読む。
    /// 値は直しを入れた日 (yyMMdd)。最後の直しは 11-18 の後方源の軸 (既定で入射の軸)</summary>
    public const string PhysicsVersion = "261118";

    // 261124Cl 追加: 形式 v2 の区画の種類・方式・大きさ
    const int KindPlanePositive = 1, KindPlaneNegative = 2, KindElectrons = 3;
    const int CodecStored = 0, CodecShuffleBrotli = 1;
    /// <summary>261124Cl 追加: 電子の 1 区画の本数 (元の大きさ 8 MB)</summary>
    internal const int ElectronsPerBlock = 1 << 17;
    const int ElectronBytes = 61; // float32 × 15 + フラグ 1 バイト (v1 と同じ値の丸め)
    const long BatchBytes = 128L << 20; // 並列に縮め・戻すときに一度に持つ量の目安 (元の大きさ。縮める間は一時に約 3 倍を持つ)
    readonly record struct Block(int Kind, int Index, long Offset, long StoredLength, long RawLength, int Codec, uint Crc);
    const int BlockEntryBytes = 40; // 区画の表の 1 行 (int32 × 2 + int64 × 3 + int32 + uint32)
    const int MaxHeaderBytes = 1 << 24; // 261124Cl 追加 (codex2 査読 中 7): 見出しの上限。読む側 (ReadHeader) と書く側 (Save) で同じ値

    /// <summary>261110Cl 追加: 読み込んだ見出しのビーム選択の記録を MasterPattern に写す (項目の無い古いファイルは従来)</summary>
    static MasterPattern WithBeamSelection(MasterPattern mp, MasterPatternFileInfo info) => new(mp.GridSize, mp.Energies, mp.Depths, mp.PositivePlanes, mp.NegativePlanes, mp.GridType)
    {
        BeamSelectionRequested = string.IsNullOrEmpty(info.BeamSelectionRequested) ? "Legacy" : info.BeamSelectionRequested,
        BeamSelectionUsed = string.IsNullOrEmpty(info.BeamSelectionUsed) ? "Legacy" : info.BeamSelectionUsed,
        BeamSelectionSettingsVersion = info.BeamSelectionSettingsVersion ?? "", BeamSelectionFallbackReason = info.BeamSelectionFallbackReason ?? "",
        BeamSelectionNative = info.BeamSelectionNative ?? "", // 261110Cl 追加
        BackSourceAxis = info.BackSourceAxis ?? "", // 261118Cl 追加
        BeamCandidateCut = info.BeamCandidateCut ?? "", // 261118Cl 追加
    };

    /// <summary>既定の拡張子</summary>
    public const string Extension = ".rpmp";

    static ReadOnlySpan<byte> Magic => "RPMP"u8;

    /// <summary>261028Cl 追加: ヘッダ (JSON) の読み書きの設定。<see cref="MasterPatternFileInfo"/> の既定の McSampleTiltRad・McEnergyThresholdKeV は NaN で、
    /// 既定の設定の JsonSerializer は NaN を書けずに ArgumentException を出していた (電子なしの保存がすべて落ちた)。NaN・±Infinity は
    /// 文字列 ("NaN" など) で書いて読む。有限の数はこれまでどおり JSON の数で書く・読むので、既存のファイルはそのまま読め、ヘッダのバイト列も変わらない。
    /// ⚠ NaN の入ったヘッダは、この修正より前の版では読めない (その版では書けもしなかった)。</summary>
    static readonly JsonSerializerOptions HeaderJsonOptions = new() { NumberHandling = JsonNumberHandling.AllowNamedFloatingPointLiterals };

    /// <summary>MasterPattern を保存する (261124Cl 変更: 形式 v2。区画ごとに並列に縮める)。</summary>
    /// <param name="compression">261124Cl 変更: 区画の Brotli の強さ (NoCompression = 縮めない、Fastest = 品質 1 (既定。実測で zlib より小さく速い)、
    ///   Optimal = 4、SmallestSize = 9)。v1 では本体の zlib の圧縮レベルだった</param>
    /// <param name="backscatteredElectrons">一緒に保存する MC の後方散乱電子 (null なら保存しない)。
    ///   入れるときは info の McSampleTiltRad / McEnergyThresholdKeV / McIncidentCount も埋めること</param>
    public static void Save(MasterPattern masterPattern, string path, MasterPatternFileInfo info,
        EbsdBackscatteredElectron[] backscatteredElectrons = null, CompressionLevel compression = CompressionLevel.Fastest)
    {
        ArgumentNullException.ThrowIfNull(masterPattern);
        ArgumentNullException.ThrowIfNull(path);
        info = FillInfo(masterPattern, info, backscatteredElectrons, FormatVersion);
        info.GridSize = masterPattern.GridSize; info.GridType = (int)masterPattern.GridType; // v2 は格子を見出しに書く (本体を読まずに条件を見せる)
        info.Energies = masterPattern.Energies; info.Depths = masterPattern.Depths;
        int quality = compression switch { CompressionLevel.NoCompression => -1, CompressionLevel.Optimal => 4, CompressionLevel.SmallestSize => 9, _ => 1 };

        // 区画の一覧 (書く順 = +Z の平面 → −Z の平面 → 電子の塊)
        // 261124Cl 追加 (codex2 査読 2 往復目 中 1): PlaneCount は int の積。あふれる大きさ・自分で読めない件数は保存しない
        long planeCountL = (long)masterPattern.Energies.Length * masterPattern.Depths.Length;
        if (planeCountL > Array.MaxLength || 2 * planeCountL + ((long)(backscatteredElectrons?.Length ?? 0) + ElectronsPerBlock - 1) / ElectronsPerBlock > 1 << 24) // 261124Cl: 足す前に long に (3 往復目 中 3)
            throw new ArgumentException("Too many planes or electrons for one file.", nameof(masterPattern));
        if (masterPattern.GridSize <= 0 || (long)masterPattern.GridSize * masterPattern.GridSize > Array.MaxLength / sizeof(float) || !System.Enum.IsDefined(masterPattern.GridType)) // 261124Cl (3 往復目 中 2): 読む側 (ValidateBlocks) と同じ検査
            throw new ArgumentException("GridSize or GridType is not valid.", nameof(masterPattern));
        int planeCount = masterPattern.PlaneCount;
        var jobs = new List<(int Kind, int Index, float[] Plane, int ElectronStart, int ElectronCount)>();
        long hemi = (long)masterPattern.GridSize * masterPattern.GridSize; // 261124Cl 追加 (codex2 査読 中 1): 平面の長さは正方・六方とも GridSize²
        foreach (var (kind, planes) in new[] { (KindPlanePositive, masterPattern.PositivePlanes), (KindPlaneNegative, masterPattern.NegativePlanes) })
            for (int k = 0; k < planeCount; k++)
            {
                var plane = k < planes.Length ? planes[k] : null;
                if (plane is { Length: > 0 } && plane.Length != hemi) throw new ArgumentException($"Plane {k} has {plane.Length} values, not GridSize² = {hemi}.", nameof(masterPattern));
                jobs.Add((kind, k, plane, 0, 0));
            }
        int nBse = backscatteredElectrons?.Length ?? 0;
        // for (int s = 0, b = 0; s < nBse; s += ElectronsPerBlock, b++) // 261124Cl 変更前: s が int の上限の近くであふれる (3 往復目 中 3)
        //     jobs.Add((KindElectrons, b, null, s, Math.Min(ElectronsPerBlock, nBse - s)));
        for (long s = 0; s < nBse; s += ElectronsPerBlock)
            jobs.Add((KindElectrons, (int)(s / ElectronsPerBlock), null, (int)s, (int)Math.Min(ElectronsPerBlock, nBse - s)));

        var header = JsonSerializer.SerializeToUtf8Bytes(info, HeaderJsonOptions);
        if (header.Length > MaxHeaderBytes) throw new ArgumentException($"The header is {header.Length} bytes (limit {MaxHeaderBytes}); shorten the note.", nameof(info)); // 261124Cl 追加 (codex2 査読 中 7): 自分で読めないファイルを書かない

        //途中で失敗しても既存のファイルを壊さないよう、一時ファイルへ書いてから置き換える
        // var tmp = path + ".tmp"; // 261124Cl 変更前: 固定の名前は、同じパスへの保存が重なると取り合う
        var tmp = $"{path}.{Guid.NewGuid():N}.tmp"; // 261124Cl 変更 (codex2 査読 中 7): 保存ごとに違う名前。失敗したら消す
        try
        {
        using (var fs = new FileStream(tmp, FileMode.Create, FileAccess.ReadWrite, FileShare.None, 1 << 20))
        {
            // var header = JsonSerializer.SerializeToUtf8Bytes(info, HeaderJsonOptions); // 261124Cl 上へ移した (上限の検査)
            Span<byte> i32 = stackalloc byte[4];
            fs.Write(Magic);
            BitConverter.TryWriteBytes(i32, FormatVersion); fs.Write(i32);
            BitConverter.TryWriteBytes(i32, header.Length); fs.Write(i32);
            fs.Write(header);
            long tablePosition = fs.Position;
            fs.Write(new byte[4 + (long)BlockEntryBytes * jobs.Count]); // 区画の表は最後に書き戻す

            var blocks = new Block[jobs.Count];
            for (int start = 0; start < jobs.Count;)
            {
                // 元の大きさで BatchBytes まで束ね、並列に縮めてから順に書く
                int end = start; long bytes = 0;
                do { var j = jobs[end++]; bytes += j.Plane != null ? 4L * j.Plane.Length : (long)ElectronBytes * j.ElectronCount; } while (end < jobs.Count && bytes < BatchBytes);
                var encoded = new (byte[] Stored, long RawLength, int Codec, uint Crc)[end - start];
                try
                {
                    System.Threading.Tasks.Parallel.For(start, end, i =>
                    {
                        var j = jobs[i];
                        if (j.Kind == KindElectrons)
                        {
                            var raw = EncodeElectrons(backscatteredElectrons, j.ElectronStart, j.ElectronCount);
                            encoded[i - start] = EncodeBlock(raw, quality);
                        }
                        else encoded[i - start] = EncodeBlock(j.Plane is { Length: > 0 } ? MemoryMarshal.AsBytes(j.Plane.AsSpan()) : [], quality);
                    });
                }
                catch (AggregateException ex) when (ex.InnerExceptions.Count > 0) // 261124Cl 追加 (codex2 査読): 読み込みと同じく中の例外をそのまま呼び手へ
                { System.Runtime.ExceptionServices.ExceptionDispatchInfo.Capture(ex.InnerExceptions[0]).Throw(); }
                for (int i = start; i < end; i++)
                {
                    var e = encoded[i - start];
                    blocks[i] = new Block(jobs[i].Kind, jobs[i].Index, fs.Position, e.Stored.Length, e.RawLength, e.Codec, e.Crc);
                    fs.Write(e.Stored);
                }
                start = end;
            }

            fs.Position = tablePosition; // 区画の表
            using (var w = new BinaryWriter(fs, Encoding.UTF8, leaveOpen: true))
            {
                w.Write(blocks.Length);
                foreach (var b in blocks)
                { w.Write(b.Kind); w.Write(b.Index); w.Write(b.Offset); w.Write(b.StoredLength); w.Write(b.RawLength); w.Write(b.Codec); w.Write(b.Crc); }
            }
        }
        File.Move(tmp, path, overwrite: true);
        }
        catch
        {
            try { File.Delete(tmp); } catch (IOException) { } catch (UnauthorizedAccessException) { } // 261124Cl 追加 (codex2 査読 中 7): 失敗の一時ファイルを残さない
            throw;
        }
    }

    /// <summary>261124Cl 追加: 保存する MasterPattern から見出しの記録を埋める (v2 の <see cref="Save"/> と v1 の <see cref="SaveVersion1"/> で共用)</summary>
    static MasterPatternFileInfo FillInfo(MasterPattern masterPattern, MasterPatternFileInfo info, EbsdBackscatteredElectron[] backscatteredElectrons, int formatVersion)
    {
        info ??= new MasterPatternFileInfo();
        info.FormatVersion = formatVersion;
        // 261110Cl 追加 (codex2 第 32 回): ビーム選択の方式は保存する MasterPattern から書く (info だけでは新方式の像を従来として記録しうる)
        info.BeamSelectionRequested = masterPattern.BeamSelectionRequested; info.BeamSelectionUsed = masterPattern.BeamSelectionUsed;
        info.BeamSelectionSettingsVersion = masterPattern.BeamSelectionSettingsVersion; info.BeamSelectionFallbackReason = masterPattern.BeamSelectionFallbackReason;
        info.BeamSelectionNative = masterPattern.BeamSelectionNative;
        info.BackSourceAxis = masterPattern.BackSourceAxis; // 261118Cl 追加
        info.BeamCandidateCut = masterPattern.BeamCandidateCut; // 261118Cl 追加
        info.BackscatteredElectronCount = backscatteredElectrons?.Length ?? 0;
        if (string.IsNullOrEmpty(info.CreatedUtc)) info.CreatedUtc = DateTime.UtcNow.ToString("o", CultureInfo.InvariantCulture);
        return info;
    }

    /// <summary>261124Cl 変更: 旧形式 (v1) で保存する。261124Cl までの <see cref="Save"/> そのもの (版の数を 1 に固定しただけ)。
    /// v1 のファイルを読めることを試験で確かめるために残す。普段の保存は <see cref="Save"/> (v2)</summary>
    internal static void SaveVersion1(MasterPattern masterPattern, string path, MasterPatternFileInfo info,
        EbsdBackscatteredElectron[] backscatteredElectrons = null, CompressionLevel compression = CompressionLevel.Fastest)
    {
        ArgumentNullException.ThrowIfNull(masterPattern);
        ArgumentNullException.ThrowIfNull(path);
        info = FillInfo(masterPattern, info, backscatteredElectrons, 1); // 261124Cl: 見出しを埋める処理は FillInfo へ (v2 と共用)。旧: 以下の行をここに直書き
        /*
        info ??= new MasterPatternFileInfo();
        info.FormatVersion = FormatVersion;
        // 261110Cl 追加 (codex2 第 32 回): ビーム選択の方式は保存する MasterPattern から書く (info だけでは新方式の像を従来として記録しうる)
        info.BeamSelectionRequested = masterPattern.BeamSelectionRequested; info.BeamSelectionUsed = masterPattern.BeamSelectionUsed;
        info.BeamSelectionSettingsVersion = masterPattern.BeamSelectionSettingsVersion; info.BeamSelectionFallbackReason = masterPattern.BeamSelectionFallbackReason;
        info.BeamSelectionNative = masterPattern.BeamSelectionNative;
        info.BackSourceAxis = masterPattern.BackSourceAxis; // 261118Cl 追加
        info.BeamCandidateCut = masterPattern.BeamCandidateCut; // 261118Cl 追加
        info.BackscatteredElectronCount = backscatteredElectrons?.Length ?? 0;
        if (string.IsNullOrEmpty(info.CreatedUtc)) info.CreatedUtc = DateTime.UtcNow.ToString("o", CultureInfo.InvariantCulture);
        */

        //途中で失敗しても既存のファイルを壊さないよう、一時ファイルへ書いてから置き換える
        var tmp = path + ".tmp";
        using (var fs = new FileStream(tmp, FileMode.Create, FileAccess.Write, FileShare.None, 1 << 20))
        {
            // var header = JsonSerializer.SerializeToUtf8Bytes(info); // 261028Cl 変更前: NaN (電子なしのときの Mc* の既定) で例外
            var header = JsonSerializer.SerializeToUtf8Bytes(info, HeaderJsonOptions); // 261028Cl 変更
            fs.Write(Magic);
            Span<byte> i32 = stackalloc byte[4];
            // BitConverter.TryWriteBytes(i32, FormatVersion); fs.Write(i32); // 261124Cl 変更前
            BitConverter.TryWriteBytes(i32, 1); fs.Write(i32); // 261124Cl 変更: v1 に固定
            BitConverter.TryWriteBytes(i32, header.Length); fs.Write(i32);
            fs.Write(header);

            using var z = new ZLibStream(fs, compression, leaveOpen: true);
            using var w = new BinaryWriter(z, Encoding.UTF8, leaveOpen: true);
            w.Write(masterPattern.GridSize);
            w.Write((int)masterPattern.GridType);
            w.Write(masterPattern.Energies.Length);
            w.Write(masterPattern.Depths.Length);
            foreach (var e in masterPattern.Energies) w.Write(e);
            foreach (var d in masterPattern.Depths) w.Write(d);
            int planeCount = masterPattern.PlaneCount;
            foreach (var planes in new[] { masterPattern.PositivePlanes, masterPattern.NegativePlanes })
                for (int k = 0; k < planeCount; k++)
                {
                    var p = k < planes.Length ? planes[k] : null;
                    w.Write(p?.Length ?? 0);
                    if (p is { Length: > 0 }) w.Write(MemoryMarshal.AsBytes(p.AsSpan()));
                }
            //電子の区画 (件数はヘッダにもあるが、本体だけで読めるようにここにも書く)
            w.Write(backscatteredElectrons?.Length ?? 0);
            if (backscatteredElectrons != null)
                foreach (var e in backscatteredElectrons) WriteElectron(w, e);
        }
        File.Move(tmp, path, overwrite: true);
    }

    /// <summary>ヘッダ (作成条件) だけを読む。本体は読まないので速い</summary>
    public static MasterPatternFileInfo ReadInfo(string path)
    {
        using var fs = new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.Read);
        return ReadHeader(fs, out _);
    }

    /// <summary>MasterPattern (と、保存されていれば後方散乱電子) を読み込む。</summary>
    /// <param name="loadElectrons">261028Cl 追加: false なら後方散乱電子の区画を読まない (返す BackscatteredElectrons は空。保存されていた本数は
    ///   Info.BackscatteredElectronCount)。電子の配列 (MnAl 1e7 で約 545 万本 ≈ 0.7 GiB の連続の領域) を確保しないので、平面だけが要るときに
    ///   峰のメモリと LOH の断片化の心配が減る。既定 (true) は従来と同じ</param>
    // public static MasterPatternFileContent Load(string path) // 261028Cl 旧シグネチャ
    public static MasterPatternFileContent Load(string path, bool loadElectrons = true)
    {
        ArgumentNullException.ThrowIfNull(path);
        using var fs = new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.Read, 1 << 20);
        // var info = ReadHeader(fs); // 261124Cl 変更前
        var info = ReadHeader(fs, out int version);
        if (version >= 2) return LoadVersion2(fs, info, loadElectrons); // 261124Cl 追加: 形式 v2。以下は v1
        using var z = new ZLibStream(fs, CompressionMode.Decompress);
        using var r = new BinaryReader(z, Encoding.UTF8);
        int gridSize = r.ReadInt32();
        var gridType = (MasterPattern.Types)r.ReadInt32();
        int nE = r.ReadInt32(), nD = r.ReadInt32();
        if (gridSize <= 0 || nE < 0 || nD < 0 || !System.Enum.IsDefined(gridType)) throw new InvalidDataException("Corrupt MasterPattern file (bad grid header).");
        // 261124Cl 追加: v2 (ValidateBlocks) と同じ上限の検査。v1 では nE*nD の int の桁あふれと、格子と合わない平面の長さを通していた (codex2 第 1 回 C1-23)
        if ((long)gridSize * gridSize > Array.MaxLength / sizeof(float) || (long)nE * nD > Array.MaxLength) throw new InvalidDataException("Corrupt MasterPattern file (bad grid header).");
        var energies = new double[nE]; for (int i = 0; i < nE; i++) energies[i] = r.ReadDouble();
        var depths = new double[nD]; for (int i = 0; i < nD; i++) depths[i] = r.ReadDouble();
        int planeCount = nE * nD;
        var pos = new float[planeCount][];
        var neg = new float[planeCount][];
        foreach (var planes in new[] { pos, neg })
            for (int k = 0; k < planeCount; k++)
            {
                int len = r.ReadInt32();
                if (len < 0) throw new InvalidDataException("Corrupt MasterPattern file (negative plane length).");
                if (len == 0) continue;
                if (len != gridSize * gridSize) throw new InvalidDataException("Corrupt MasterPattern file (plane length does not match the grid)."); // 261124Cl 追加 (C1-23)
                var p = new float[len];
                r.BaseStream.ReadExactly(MemoryMarshal.AsBytes(p.AsSpan()));
                planes[k] = p;
            }
        // 261028Cl 追加: 電子を読まないときは、平面のあとの区画 (電子の本数と本体) を展開せずに返す
        if (!loadElectrons)
            // return new MasterPatternFileContent(new MasterPattern(gridSize, energies, depths, pos, neg, gridType), info, []); // 261110Cl 変更前
            return new MasterPatternFileContent(WithBeamSelection(new MasterPattern(gridSize, energies, depths, pos, neg, gridType), info), info, []); // 261110Cl: ビーム選択の記録を写す
        int nBse = r.ReadInt32();
        if (nBse < 0) throw new InvalidDataException("Corrupt MasterPattern file (negative electron count).");
        var bses = new EbsdBackscatteredElectron[nBse];
        for (int i = 0; i < nBse; i++) bses[i] = ReadElectron(r);
        // return new MasterPatternFileContent(new MasterPattern(gridSize, energies, depths, pos, neg, gridType), info, bses); // 261110Cl 変更前
        return new MasterPatternFileContent(WithBeamSelection(new MasterPattern(gridSize, energies, depths, pos, neg, gridType), info), info, bses); // 261110Cl
    }

    /// <summary>261031Cl 追加: 保存した後方散乱電子を、chunkSize 本ずつの配列 (塊) にしてファイルの並び順に返す。
    /// <see cref="Load"/> と違って電子の全部を 1 本の配列 (MnAl 1e7 で約 545 万本 ≈ 0.7 GiB の連続の領域) にしないので、GC のヒープの上限のもとや
    /// LOH が断片化した後でも読める。1 本ずつの読み方は <see cref="Load"/> と同じなので、塊をつなぐと Load の BackscatteredElectrons とビットで同じになる。
    /// <para>平面は配列を作らずに読み飛ばす。ただし本体は 1 本の zlib の流れなので、電子の区画へ行くには平面の分も展開する (平面の多いファイルでは
    /// その分の時間がかかる)。平面も要るときは <see cref="Load"/> (loadElectrons: false) を別に呼ぶ。作成条件 (MC の試料傾斜など) は <see cref="ReadInfo"/> で読む。</para>
    /// <para>ファイルは列挙を始めたとき (最初の MoveNext) に開き、最後まで読んだとき・foreach を途中で抜けたとき (列挙子の Dispose)・例外のときに閉じる。
    /// 列挙している間はファイルを (FileShare.Read で) 開いたままなので、その間は同じパスへの <see cref="Save"/> (一時ファイルからの置き換え) が
    /// 例外で失敗する (path.tmp が残る)。Load にも同じ窓はあるが、こちらは呼び手が塊を処理する時間だけ長くなる。
    /// 列挙するたびに初めから読み直す。返す塊は毎回新しい配列なので、持っておいてよい。電子が 0 本のファイルでは何も返さない。</para></summary>
    /// <param name="chunkSize">1 つの塊の本数 (1 以上。最後の塊は端数)。既定の 512 本 (69,632 B。今の構造体 136 B/本 のとき) は LOH の閾値 (85,000 B) より小さい。
    ///   全本数より大きくすると全部が 1 つの塊になる (Load と同じ大きさの配列になる)</param>
    /// <exception cref="ArgumentNullException">path が null (列挙の前、呼んだときに出す)</exception>
    /// <exception cref="ArgumentOutOfRangeException">chunkSize が 0 以下 (列挙の前、呼んだときに出す)</exception>
    public static IEnumerable<EbsdBackscatteredElectron[]> ReadElectronChunks(string path, int chunkSize = 512)
    {
        ArgumentNullException.ThrowIfNull(path);
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(chunkSize);
        return Chunks(); //引数の検査を列挙の前に済ませるため、反復子の本体は局所関数にする

        IEnumerable<EbsdBackscatteredElectron[]> Chunks()
        {
            using var fs = new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.Read, 1 << 20);
            // ReadHeader(fs); // 261124Cl 変更前
            var info = ReadHeader(fs, out int version);
            if (version >= 2) // 261124Cl 追加: 形式 v2 は電子の区画だけを順に戻す (平面は読まない)。以下は v1
            {
                foreach (var c in ElectronChunksVersion2(fs, info, chunkSize)) yield return c;
                yield break;
            }
            using var z = new ZLibStream(fs, CompressionMode.Decompress);
            using var r = new BinaryReader(z, Encoding.UTF8);
            //本体の書式は Load と同じに読む (形式を変えたら両方を直す。FormatVersion の remarks)。格子の見出しも Load と同じに検査する
            int gridSize = r.ReadInt32();
            var gridType = (MasterPattern.Types)r.ReadInt32();
            int nE = r.ReadInt32(), nD = r.ReadInt32();
            if (gridSize <= 0 || nE < 0 || nD < 0 || !System.Enum.IsDefined(gridType)) throw new InvalidDataException("Corrupt MasterPattern file (bad grid header).");
            // 261124Cl 追加: Load (v1) と同じ上限の検査 (codex2 第 3 回: こちらの入口は負の長さしか見ていなかった)
            if ((long)gridSize * gridSize > Array.MaxLength / sizeof(float) || (long)nE * nD > Array.MaxLength) throw new InvalidDataException("Corrupt MasterPattern file (bad grid header).");
            //for (int i = 0; i < nE + nD; i++) r.ReadDouble(); //エネルギー列・深さ列は読み捨てる // 261124Cl 変更前: nE + nD が int
            for (long i = 0; i < (long)nE + nD; i++) r.ReadDouble(); //エネルギー列・深さ列は読み捨てる // 261124Cl
            //平面は配列を作らず、長さ × 4 バイトを捨て緩衝 (64 KiB、LOH に入らない) へ読み捨てる (+Z → −Z の 2 × nE·nD 枚)
            var discard = new byte[1 << 16];
            for (long k = 0, planeCount = 2L * nE * nD; k < planeCount; k++)
            {
                int len = r.ReadInt32();
                if (len < 0) throw new InvalidDataException("Corrupt MasterPattern file (negative plane length).");
                if (len != 0 && len != gridSize * gridSize) throw new InvalidDataException("Corrupt MasterPattern file (plane length does not match the grid)."); // 261124Cl 追加 (Load と同じ)
                for (long left = 4L * len; left > 0; left -= discard.Length)
                    r.BaseStream.ReadExactly(discard, 0, (int)Math.Min(left, discard.Length));
            }
            int nBse = r.ReadInt32();
            if (nBse < 0) throw new InvalidDataException("Corrupt MasterPattern file (negative electron count).");
            for (int done = 0; done < nBse;)
            {
                var chunk = new EbsdBackscatteredElectron[Math.Min(chunkSize, nBse - done)];
                for (int i = 0; i < chunk.Length; i++) chunk[i] = ReadElectron(r);
                done += chunk.Length;
                yield return chunk;
            }
        }
    }

    //電子 1 本 = float32 × 15 + フラグ 1 バイト (61 バイト)。MC の統計誤差に比べて float の丸めは無視できる
    static void WriteElectron(BinaryWriter w, in EbsdBackscatteredElectron e)
    {
        w.Write((float)e.Depth); w.Write((float)e.Vec.X); w.Write((float)e.Vec.Y); w.Write((float)e.Vec.Z);
        w.Write((float)e.Position.X); w.Write((float)e.Position.Y); w.Write((float)e.Energy); w.Write((float)e.TotalEnergyLoss);
        w.Write((float)e.LastInelasticDepth); w.Write((float)e.LastInelasticEnergyBeforeLoss); w.Write((float)e.LastInelasticEnergyAfterLoss);
        w.Write((float)e.LastInelasticDirection.X); w.Write((float)e.LastInelasticDirection.Y); w.Write((float)e.LastInelasticDirection.Z);
        w.Write((float)e.LastDecoherenceDepth);
        w.Write((byte)((e.HasLastInelasticEvent ? 1 : 0) | (e.HasLastDecoherenceEvent ? 2 : 0)));
    }

    static EbsdBackscatteredElectron ReadElectron(BinaryReader r)
    {
        double F() => r.ReadSingle();
        double depth = F(); var vec = new OpenTK.Mathematics.Vector3d(F(), F(), F());
        var position = new PointD(F(), F()); double energy = F(), totalLoss = F();
        double lastInelDepth = F(), lastInelBefore = F(), lastInelAfter = F();
        var lastInelDir = new OpenTK.Mathematics.Vector3d(F(), F(), F());
        double lastDecoDepth = F();
        byte flags = r.ReadByte();
        return new EbsdBackscatteredElectron(depth, vec, position, energy, totalLoss, (flags & 1) != 0, lastInelDepth, lastInelBefore, lastInelAfter, lastInelDir,
            (flags & 2) != 0, lastDecoDepth);
    }

    // static MasterPatternFileInfo ReadHeader(Stream fs) // 261124Cl 旧シグネチャ
    static MasterPatternFileInfo ReadHeader(Stream fs, out int version) // 261124Cl 変更: 版も返す (v1 と v2 で本体の読み方が違う)
    {
        Span<byte> magic = stackalloc byte[4];
        fs.ReadExactly(magic);
        if (!magic.SequenceEqual(Magic)) throw new InvalidDataException("Not a ReciPro MasterPattern file.");
        Span<byte> i32 = stackalloc byte[4];
        // fs.ReadExactly(i32); int version = BitConverter.ToInt32(i32); // 261124Cl 変更前
        fs.ReadExactly(i32); version = BitConverter.ToInt32(i32);
        if (version < 1 || version > FormatVersion) throw new InvalidDataException($"Unsupported MasterPattern file version {version} (this build reads up to {FormatVersion}).");
        fs.ReadExactly(i32); int headerLength = BitConverter.ToInt32(i32);
        // if (headerLength < 0 || headerLength > 1 << 24) throw new InvalidDataException("Corrupt MasterPattern file (bad header length)."); // 261124Cl 変更前
        if (headerLength < 0 || headerLength > MaxHeaderBytes) throw new InvalidDataException("Corrupt MasterPattern file (bad header length)."); // 261124Cl 変更: 書く側と同じ上限
        var header = new byte[headerLength];
        fs.ReadExactly(header);
        // return JsonSerializer.Deserialize<MasterPatternFileInfo>(header) ?? new MasterPatternFileInfo(); // 261028Cl 変更前
        return JsonSerializer.Deserialize<MasterPatternFileInfo>(header, HeaderJsonOptions) ?? new MasterPatternFileInfo(); // 261028Cl 変更: "NaN" なども読む
    }

    #region 261124Cl 追加: 形式 v2 (区画・並べ替え + Brotli・CRC32)

    /// <summary>形式 v2 の本体を読む (見出しの直後から)。平面は区画ごとに並列に戻す</summary>
    static MasterPatternFileContent LoadVersion2(FileStream fs, MasterPatternFileInfo info, bool loadElectrons)
    {
        var gridType = (MasterPattern.Types)info.GridType;
        int nE = info.Energies?.Length ?? -1, nD = info.Depths?.Length ?? -1;
        if (info.GridSize <= 0 || nE < 0 || nD < 0 || !System.Enum.IsDefined(gridType)) throw new InvalidDataException("Corrupt MasterPattern file (bad grid header).");
        var blocks = ReadBlockTable(fs, info);
        ValidateBlocks(blocks, info, fs.Position, fs.Length); // 261124Cl 追加 (codex2 査読 高 2・3): 番号・本数・長さ・位置を確かめてから確保する
        int planeCount = nE * nD; // ValidateBlocks が Array.MaxLength 以下を確かめた
        var pos = new float[planeCount][];
        var neg = new float[planeCount][];
        ReadBlocks(fs, blocks.Where(b => b.Kind is KindPlanePositive or KindPlaneNegative), (b, raw) =>
        {
            if ((uint)b.Index >= (uint)planeCount || raw.Length % 4 != 0) throw new InvalidDataException("Corrupt MasterPattern file (bad plane block).");
            if (raw.Length == 0) return; // 欠落の平面
            var p = new float[raw.Length / 4];
            raw.AsSpan().CopyTo(MemoryMarshal.AsBytes(p.AsSpan()));
            (b.Kind == KindPlanePositive ? pos : neg)[b.Index] = p;
        });
        EbsdBackscatteredElectron[] bses = [];
        if (loadElectrons)
        {
            var eb = blocks.Where(b => b.Kind == KindElectrons).ToArray();
            long n = eb.Sum(b => b.RawLength / ElectronBytes);
            if (n != info.BackscatteredElectronCount || n > Array.MaxLength) throw new InvalidDataException("Corrupt MasterPattern file (electron count).");
            bses = new EbsdBackscatteredElectron[n];
            ReadBlocks(fs, eb, (b, raw) =>
            {
                int count = (int)(raw.Length / ElectronBytes);
                long start = (long)b.Index * ElectronsPerBlock;
                if (raw.Length % ElectronBytes != 0 || start + count > n) throw new InvalidDataException("Corrupt MasterPattern file (bad electron block).");
                for (int i = 0; i < count; i++) bses[start + i] = DecodeElectron(raw, count, i);
            });
        }
        return new MasterPatternFileContent(WithBeamSelection(new MasterPattern(info.GridSize, info.Energies, info.Depths, pos, neg, gridType), info), info, bses);
    }

    /// <summary>形式 v2 の電子の区画を順に戻し、chunkSize 本ずつの配列で返す (<see cref="ReadElectronChunks"/>)。平面の区画は読まない。
    /// 261124Cl 変更 (codex2 査読 中 2): 区画の境をまたいで塊を詰める (v1 と同じく、最後の塊だけが端数。旧は区画ごとに端数の塊を返していた)</summary>
    // static IEnumerable<EbsdBackscatteredElectron[]> ElectronChunksVersion2(FileStream fs, int chunkSize) // 261124Cl 旧シグネチャ
    static IEnumerable<EbsdBackscatteredElectron[]> ElectronChunksVersion2(FileStream fs, MasterPatternFileInfo info, int chunkSize)
    {
        var blocks = ReadBlockTable(fs, info);
        ValidateBlocks(blocks, info, fs.Position, fs.Length); // 261124Cl 追加 (codex2 査読 高 2): 連番・本数を確かめる (総数も見出しと照合)
        long remaining = info.BackscatteredElectronCount;
        EbsdBackscatteredElectron[] chunk = null;
        int filled = 0;
        foreach (var b in blocks.Where(b => b.Kind == KindElectrons).OrderBy(b => b.Index))
        {
            var raw = new byte[b.RawLength];
            DecodeBlock(b, ReadStored(fs, b), raw);
            int count = raw.Length / ElectronBytes;
            // for (int done = 0; done < count;) // 261124Cl 変更前 (区画ごとに端数の塊)
            // {
            //     var chunk = new EbsdBackscatteredElectron[Math.Min(chunkSize, count - done)];
            //     for (int i = 0; i < chunk.Length; i++) chunk[i] = DecodeElectron(raw, count, done + i);
            //     done += chunk.Length;
            //     yield return chunk;
            // }
            for (int i = 0; i < count; i++)
            {
                chunk ??= new EbsdBackscatteredElectron[(int)Math.Min(chunkSize, remaining)];
                chunk[filled++] = DecodeElectron(raw, count, i);
                remaining--;
                if (filled == chunk.Length) { yield return chunk; chunk = null; filled = 0; }
            }
        }
    }

    /// <summary>261124Cl 追加 (codex2 査読 高 2・3): 区画の表を確かめる。平面は ± それぞれ 0〜planeCount−1 の番号がちょうど 1 回ずつ、長さは 0 (欠落) か GridSize² × 4 バイト。
    /// 電子は番号が 0 からの連番で、最後の区画以外はちょうど <see cref="ElectronsPerBlock"/> 本、合計は見出しの本数。区画は表の後ろ・ファイルの中に収まり、互いに重ならない。
    /// 壊れた表から巨大な配列を確保したり、同じ平面を上書きして別の値を返したりしないため、読む前に呼ぶ</summary>
    static void ValidateBlocks(Block[] blocks, MasterPatternFileInfo info, long tableEnd, long fileLength)
    {
        static InvalidDataException Bad(string what) => new($"Corrupt MasterPattern file ({what}).");
        // 261124Cl 追加 (codex2 査読 2 往復目 高 3): 格子の検査を、欠落の判定より先に、Load と ReadElectronChunks の両方で (全平面が欠落でも巨大な格子を通さない)
        if (info.GridSize <= 0 || (long)info.GridSize * info.GridSize > Array.MaxLength / sizeof(float) || !System.Enum.IsDefined((MasterPattern.Types)info.GridType)
            || info.Energies == null || info.Depths == null)
            throw Bad("bad grid header");
        long planeCount = (long)info.Energies.Length * info.Depths.Length;
        if (planeCount > Array.MaxLength || 2 * planeCount > blocks.Length) throw Bad("too many planes for the block table"); // 平面は ± で 2·planeCount 区画あるはず
        long planeBytes = 4L * info.GridSize * info.GridSize;
        var seen = new[] { new bool[planeCount], new bool[planeCount] };
        var electrons = new List<Block>();
        foreach (var b in blocks)
        {
            if (b.Offset < tableEnd || b.StoredLength > fileLength - b.Offset) throw Bad("block outside the file");
            if (b.Codec is not (CodecStored or CodecShuffleBrotli)) throw Bad($"unsupported block codec {b.Codec}");
            if (b.Codec == CodecStored && b.StoredLength != b.RawLength) throw Bad("bad stored block length");
            if (b.Kind is KindPlanePositive or KindPlaneNegative)
            {
                var s = seen[b.Kind == KindPlanePositive ? 0 : 1];
                if ((ulong)b.Index >= (ulong)planeCount || s[b.Index]) throw Bad("duplicate or out-of-range plane block");
                s[b.Index] = true;
                if (b.RawLength != 0 && b.RawLength != planeBytes) throw Bad("plane length does not match the grid");
            }
            else if (b.Kind == KindElectrons) electrons.Add(b);
            else throw Bad($"unknown block kind {b.Kind}");
        }
        if (seen.Any(s => s.Contains(false))) throw Bad("missing plane block");
        electrons.Sort((a, b) => a.Index.CompareTo(b.Index));
        long total = 0;
        for (int i = 0; i < electrons.Count; i++)
        {
            var b = electrons[i];
            long count = b.RawLength / ElectronBytes;
            bool last = i == electrons.Count - 1;
            if (b.Index != i || b.RawLength % ElectronBytes != 0 || (last ? count is < 1 or > ElectronsPerBlock : count != ElectronsPerBlock))
                throw Bad("bad electron block sequence");
            total += count;
        }
        if (total != info.BackscatteredElectronCount || total > Array.MaxLength) throw Bad("electron count");
        var byOffset = blocks.Where(b => b.StoredLength > 0).OrderBy(b => b.Offset).ToArray();
        for (int i = 1; i < byOffset.Length; i++)
            if (byOffset[i].Offset < byOffset[i - 1].Offset + byOffset[i - 1].StoredLength) throw Bad("overlapping blocks");
    }

    // static Block[] ReadBlockTable(Stream fs) // 261124Cl 旧シグネチャ (同日)
    static Block[] ReadBlockTable(Stream fs, MasterPatternFileInfo info)
    {
        using var r = new BinaryReader(fs, Encoding.UTF8, leaveOpen: true);
        int n = r.ReadInt32();
        if (n < 0 || n > 1 << 24) throw new InvalidDataException("Corrupt MasterPattern file (bad block table).");
        // 261124Cl 追加 (codex2 査読 2 往復目 中 1): 件数は見出しだけで決まる (± の平面 + 電子の塊)。確保の前に照らす
        long expected = 2L * (info.Energies?.Length ?? 0) * (info.Depths?.Length ?? 0) + ((long)info.BackscatteredElectronCount + ElectronsPerBlock - 1) / ElectronsPerBlock; // 261124Cl: 足す前に long に (3 往復目 中 3)
        if (n != expected) throw new InvalidDataException("Corrupt MasterPattern file (block count does not match the header).");
        if ((long)n * BlockEntryBytes > fs.Length - fs.Position) throw new InvalidDataException("Corrupt MasterPattern file (block table beyond the end of the file)."); // 261124Cl 追加 (codex2 査読 高 3): 確保の前に
        var blocks = new Block[n];
        for (int i = 0; i < n; i++)
        {
            blocks[i] = new Block(r.ReadInt32(), r.ReadInt32(), r.ReadInt64(), r.ReadInt64(), r.ReadInt64(), r.ReadInt32(), r.ReadUInt32());
            if (blocks[i].Offset < 0 || blocks[i].StoredLength < 0 || blocks[i].RawLength < 0 || blocks[i].RawLength > Array.MaxLength || blocks[i].StoredLength > Array.MaxLength)
                throw new InvalidDataException("Corrupt MasterPattern file (bad block entry).");
        }
        return blocks;
    }

    static byte[] ReadStored(Stream fs, in Block b)
    {
        var stored = new byte[b.StoredLength];
        fs.Position = b.Offset;
        fs.ReadExactly(stored);
        return stored;
    }

    /// <summary>区画を BatchBytes ずつ順に読み、並列に戻して consume に渡す (consume は区画ごとに別のスレッドから呼ばれる)</summary>
    static void ReadBlocks(FileStream fs, IEnumerable<Block> blocks, Action<Block, byte[]> consume)
    {
        var list = blocks.OrderBy(b => b.Offset).ToArray();
        for (int start = 0; start < list.Length;)
        {
            int end = start; long bytes = 0;
            var stored = new List<byte[]>();
            // do { stored.Add(ReadStored(fs, list[end])); bytes += list[end].RawLength; end++; } while (end < list.Length && bytes < BatchBytes); // 261124Cl 変更前
            do { stored.Add(ReadStored(fs, list[end])); bytes += list[end].RawLength + list[end].StoredLength; end++; } while (end < list.Length && bytes < BatchBytes); // 261124Cl 変更 (codex2 査読 2 往復目 低): 格納長も数える
            try
            {
                System.Threading.Tasks.Parallel.For(start, end, i =>
                {
                    var raw = new byte[list[i].RawLength];
                    DecodeBlock(list[i], stored[i - start], raw);
                    consume(list[i], raw);
                });
            }
            catch (AggregateException ex) when (ex.InnerExceptions.Count > 0) // 「壊れている」の理由 (InvalidDataException) をそのまま呼び手へ
            { System.Runtime.ExceptionServices.ExceptionDispatchInfo.Capture(ex.InnerExceptions[0]).Throw(); }
            start = end;
        }
    }

    /// <summary>1 区画を縮める。並べ替え + Brotli で縮まなければそのまま格納する。CRC32 は元のバイト列のもの</summary>
    static (byte[] Stored, long RawLength, int Codec, uint Crc) EncodeBlock(ReadOnlySpan<byte> raw, int quality)
    {
        uint crc = Crc32(raw);
        if (raw.Length > 0 && quality >= 0)
        {
            var shuffled = new byte[raw.Length];
            Shuffle4(raw, shuffled);
            var dst = new byte[System.IO.Compression.BrotliEncoder.GetMaxCompressedLength(raw.Length)];
            if (System.IO.Compression.BrotliEncoder.TryCompress(shuffled, dst, out int n, quality, 22) && n < raw.Length)
                return (dst[..n], raw.Length, CodecShuffleBrotli, crc);
        }
        return (raw.ToArray(), raw.Length, CodecStored, crc);
    }

    /// <summary>1 区画を戻して raw (元の長さちょうど) に書き、CRC32 を照合する</summary>
    static void DecodeBlock(in Block b, ReadOnlySpan<byte> stored, Span<byte> raw)
    {
        if (b.Codec == CodecStored)
        {
            if (stored.Length != raw.Length) throw new InvalidDataException("Corrupt MasterPattern file (bad block length).");
            stored.CopyTo(raw);
        }
        else if (b.Codec == CodecShuffleBrotli)
        {
            var shuffled = new byte[raw.Length];
            if (!System.IO.Compression.BrotliDecoder.TryDecompress(stored, shuffled, out int n) || n != raw.Length)
                throw new InvalidDataException("Corrupt MasterPattern file (block cannot be decompressed).");
            Unshuffle4(shuffled, raw);
        }
        else throw new InvalidDataException($"Unsupported MasterPattern block codec {b.Codec}.");
        if (Crc32(raw) != b.Crc) throw new InvalidDataException($"Corrupt MasterPattern file (checksum mismatch in block {b.Kind}:{b.Index}).");
    }

    /// <summary>4 バイトの並べ替え (Blosc の shuffle)。i 番目の 4 バイトの j 位を j·m + i へ (m = 4 バイトの個数)。4 で割り切れない末尾はそのまま</summary>
    static void Shuffle4(ReadOnlySpan<byte> src, Span<byte> dst)
    {
        int m = src.Length / 4;
        for (int i = 0; i < m; i++)
            for (int j = 0; j < 4; j++) dst[j * m + i] = src[4 * i + j];
        src[(4 * m)..].CopyTo(dst[(4 * m)..]);
    }

    static void Unshuffle4(ReadOnlySpan<byte> src, Span<byte> dst)
    {
        int m = src.Length / 4;
        for (int i = 0; i < m; i++)
            for (int j = 0; j < 4; j++) dst[4 * i + j] = src[j * m + i];
        src[(4 * m)..].CopyTo(dst[(4 * m)..]);
    }

    static readonly uint[] crcTable = Enumerable.Range(0, 256).Select(n =>
    {
        uint c = (uint)n;
        for (int k = 0; k < 8; k++) c = (c & 1) != 0 ? 0xEDB88320u ^ (c >> 1) : c >> 1;
        return c;
    }).ToArray();

    /// <summary>CRC-32 (IEEE 802.3、zlib と同じ多項式)</summary>
    static uint Crc32(ReadOnlySpan<byte> data)
    {
        uint c = 0xFFFFFFFFu;
        foreach (var b in data) c = crcTable[(c ^ b) & 0xFF] ^ (c >> 8);
        return c ^ 0xFFFFFFFFu;
    }

    /// <summary>電子を列ごとに並べる (float32 × 15 列 + フラグ 1 列)。各列の値の丸めは v1 (<see cref="WriteElectron"/>) と同じ</summary>
    static byte[] EncodeElectrons(EbsdBackscatteredElectron[] electrons, int start, int count)
    {
        var raw = new byte[(long)ElectronBytes * count];
        var f = MemoryMarshal.Cast<byte, float>(raw.AsSpan(0, 60 * count));
        for (int i = 0; i < count; i++)
        {
            ref readonly var e = ref electrons[start + i];
            f[0 * count + i] = (float)e.Depth; f[1 * count + i] = (float)e.Vec.X; f[2 * count + i] = (float)e.Vec.Y; f[3 * count + i] = (float)e.Vec.Z;
            f[4 * count + i] = (float)e.Position.X; f[5 * count + i] = (float)e.Position.Y; f[6 * count + i] = (float)e.Energy; f[7 * count + i] = (float)e.TotalEnergyLoss;
            f[8 * count + i] = (float)e.LastInelasticDepth; f[9 * count + i] = (float)e.LastInelasticEnergyBeforeLoss; f[10 * count + i] = (float)e.LastInelasticEnergyAfterLoss;
            f[11 * count + i] = (float)e.LastInelasticDirection.X; f[12 * count + i] = (float)e.LastInelasticDirection.Y; f[13 * count + i] = (float)e.LastInelasticDirection.Z;
            f[14 * count + i] = (float)e.LastDecoherenceDepth;
            raw[60 * count + i] = (byte)((e.HasLastInelasticEvent ? 1 : 0) | (e.HasLastDecoherenceEvent ? 2 : 0));
        }
        return raw;
    }

    static EbsdBackscatteredElectron DecodeElectron(byte[] raw, int count, int i)
    {
        var f = MemoryMarshal.Cast<byte, float>(raw.AsSpan(0, 60 * count));
        byte flags = raw[60 * count + i];
        return new EbsdBackscatteredElectron(f[i], new OpenTK.Mathematics.Vector3d(f[count + i], f[2 * count + i], f[3 * count + i]),
            new PointD(f[4 * count + i], f[5 * count + i]), f[6 * count + i], f[7 * count + i], (flags & 1) != 0,
            f[8 * count + i], f[9 * count + i], f[10 * count + i], new OpenTK.Mathematics.Vector3d(f[11 * count + i], f[12 * count + i], f[13 * count + i]),
            (flags & 2) != 0, f[14 * count + i]);
    }

    #endregion

    /// <summary>結晶の同一性を表す文字列。格子定数 (nm, 6 桁)・角度 (度, 4 桁)・空間群の系列番号・
    /// 原子ごとの Z / 座標 (5 桁) / 占有率 (4 桁) / 等方 B (nm², 6 桁) を並べる。異方性の原子はその後ろに β の 6 成分 (有効数字 7 桁) を足す (261105Cl)。
    /// 名前は含めない (名前だけ変えた同じ結晶は同一とみなす)。</summary>
    public static string CrystalFingerprint(Crystal crystal)
    {
        ArgumentNullException.ThrowIfNull(crystal);
        var ci = CultureInfo.InvariantCulture;
        var sb = new StringBuilder();
        sb.Append(string.Format(ci, "cell={0:F6},{1:F6},{2:F6},{3:F4},{4:F4},{5:F4};sg={6};",
            crystal.A, crystal.B, crystal.C, crystal.Alpha * 180 / Math.PI, crystal.Beta * 180 / Math.PI, crystal.Gamma * 180 / Math.PI, crystal.SymmetrySeriesNumber));
        foreach (var a in crystal.Atoms)
        {
            //sb.Append(string.Format(ci, "{0}:{1:F5},{2:F5},{3:F5},{4:F4},{5:F6};", a.AtomicNumber, a.X, a.Y, a.Z, a.Occ, a.Dsf?.BisoEffective ?? 0)); // 261105Cl 変更前
            sb.Append(string.Format(ci, "{0}:{1:F5},{2:F5},{3:F5},{4:F4},{5:F6}", a.AtomicNumber, a.X, a.Y, a.Z, a.Occ, a.Dsf?.BisoEffective ?? 0));
            // 261105Cl 追加: 異方性の原子は β の 6 成分 (無次元, 11,22,33,12,23,31) も入れる。等価な B だけでは、非対角成分だけが違う結晶
            // (構造因子は違う) を区別できない (Codex2 r3)。また 261105Cl の温度因子の修正 (β の非対角の係数・B_eq の式) の前に作った異方性の結晶のファイルは、
            // ここが違うので読み込まれない (古い計算のパターンを受け入れない)。等方の原子は以前と同じ文字列なので、等方の結晶のファイルはそのまま読める。
            // 数値だけをカンマで足すので、指紋を読む側 (tools/EbsdOccupancy の ParseFingerprint) は先頭の 5 つだけを読んで動く
            if (a.Dsf is { UseIso: false } d)
                sb.Append(string.Format(ci, ",{0:E6},{1:E6},{2:E6},{3:E6},{4:E6},{5:E6}", d.B11, d.B22, d.B33, d.B12, d.B23, d.B31));
            sb.Append(';');
        }
        return sb.ToString();
    }

    /// <summary>作成条件が今の結晶・電圧と一致するか。一致しなければ理由を返す (一致なら null)</summary>
    /// <param name="beamSelectionUsed">261110Cl 追加 (codex2 第 34 回): null 以外なら、ビーム選択の方式 (Legacy・DressedTail) も一致を求める (較正や複数のマスターを組み合わせる道具用。表示は従来どおり方式を問わない)</param>
    /// <param name="beamSelectionSettingsVersion">261110Cl 追加: null 以外なら、新方式の設定の版 (DT1-261110、+custom-指紋) も一致を求める</param>
    // public static string CheckCompatibility(MasterPatternFileInfo info, Crystal crystal, double beamEnergyKeV) // 261110Cl 旧シグネチャ
    /// <param name="beamCandidateCut">261118Cl 追加 (codex2 第 15 回): null 以外なら、従来のビーム選択の候補の打ち切り方 (Legacy・Shell) の一致も求める。記録の無い ("") ファイルは一致としない</param>
    // public static string CheckCompatibility(MasterPatternFileInfo info, Crystal crystal, double beamEnergyKeV, string beamSelectionUsed = null, string beamSelectionSettingsVersion = null) // 261118Cl 旧シグネチャ
    public static string CheckCompatibility(MasterPatternFileInfo info, Crystal crystal, double beamEnergyKeV, string beamSelectionUsed = null, string beamSelectionSettingsVersion = null,
        string beamCandidateCut = null)
    {
        ArgumentNullException.ThrowIfNull(info);
        ArgumentNullException.ThrowIfNull(crystal);
        if (beamCandidateCut is "") throw new ArgumentException("Pass null to skip the beam-candidate-cut check; an empty requirement would match unrecorded files.", nameof(beamCandidateCut)); // 261118Cl 追加 (codex2 第 16 回)
        if (info.CrystalFingerprint != CrystalFingerprint(crystal))
            return $"The file was made for a different crystal (\"{info.CrystalName}\"), or the cell, atoms or B differ.";
        if (!double.IsFinite(info.BeamEnergyKeV) || info.BeamEnergyKeV <= 0) // 261124Cl 追加 (codex2 査読 中 5): NaN は下の比べで「一致」になってしまう
            return "The file does not record a valid beam energy.";
        if (!double.IsFinite(beamEnergyKeV) || beamEnergyKeV <= 0) // 261124Cl 追加 (codex2 査読 2 往復目 中 2): 比べ先が NaN でも「一致」になる
            return $"The current beam energy ({beamEnergyKeV:g6} keV) is not valid.";
        if (Math.Abs(info.BeamEnergyKeV - beamEnergyKeV) > 1E-6)
            return $"The file was made at {info.BeamEnergyKeV:g6} keV, not {beamEnergyKeV:g6} keV.";
        // 261110Cl 追加 (codex2 第 34 回): 方式と設定の版 (求められたときだけ。古いファイルは Legacy・版なしとして読まれている)
        if (beamSelectionUsed != null && (info.BeamSelectionUsed ?? "Legacy") != beamSelectionUsed)
            return $"The file was made with the {info.BeamSelectionUsed ?? "Legacy"} beam selection, not {beamSelectionUsed}.";
        if (beamSelectionSettingsVersion != null && (info.BeamSelectionSettingsVersion ?? "") != beamSelectionSettingsVersion)
            return $"The file was made with the beam-selection settings '{info.BeamSelectionSettingsVersion}', not '{beamSelectionSettingsVersion}'.";
        if (beamCandidateCut != null && (info.BeamCandidateCut ?? "") != beamCandidateCut) // 261118Cl 追加 (codex2 第 15 回)
            return string.IsNullOrEmpty(info.BeamCandidateCut)
                ? $"The file does not record the beam candidate cut (made before it was recorded, when the default was Legacy, or with the dressed-tail selection); {beamCandidateCut} was required."
                : $"The file was made with the {info.BeamCandidateCut} beam candidate cut, not {beamCandidateCut}.";
        return null;
    }
}
