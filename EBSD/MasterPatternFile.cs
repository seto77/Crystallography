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
/// </summary>
public static class MasterPatternFile
{
    /// <summary>ファイル形式の版。形式を変えたら上げ、<see cref="Load"/> で旧版を読めるようにすること</summary>
    /// <remarks>261031Cl 追加: 本体 (格子の見出し・平面の並び・電子の区画) は <see cref="ReadElectronChunks"/> でも読む。形式を変えたら両方を直すこと</remarks>
    public const int FormatVersion = 1;

    /// <summary>261110Cl 追加: 読み込んだ見出しのビーム選択の記録を MasterPattern に写す (項目の無い古いファイルは従来)</summary>
    static MasterPattern WithBeamSelection(MasterPattern mp, MasterPatternFileInfo info) => new(mp.GridSize, mp.Energies, mp.Depths, mp.PositivePlanes, mp.NegativePlanes, mp.GridType)
    {
        BeamSelectionRequested = string.IsNullOrEmpty(info.BeamSelectionRequested) ? "Legacy" : info.BeamSelectionRequested,
        BeamSelectionUsed = string.IsNullOrEmpty(info.BeamSelectionUsed) ? "Legacy" : info.BeamSelectionUsed,
        BeamSelectionSettingsVersion = info.BeamSelectionSettingsVersion ?? "", BeamSelectionFallbackReason = info.BeamSelectionFallbackReason ?? "",
    };

    /// <summary>既定の拡張子</summary>
    public const string Extension = ".rpmp";

    static ReadOnlySpan<byte> Magic => "RPMP"u8;

    /// <summary>261028Cl 追加: ヘッダ (JSON) の読み書きの設定。<see cref="MasterPatternFileInfo"/> の既定の McSampleTiltRad・McEnergyThresholdKeV は NaN で、
    /// 既定の設定の JsonSerializer は NaN を書けずに ArgumentException を出していた (電子なしの保存がすべて落ちた)。NaN・±Infinity は
    /// 文字列 ("NaN" など) で書いて読む。有限の数はこれまでどおり JSON の数で書く・読むので、既存のファイルはそのまま読め、ヘッダのバイト列も変わらない。
    /// ⚠ NaN の入ったヘッダは、この修正より前の版では読めない (その版では書けもしなかった)。</summary>
    static readonly JsonSerializerOptions HeaderJsonOptions = new() { NumberHandling = JsonNumberHandling.AllowNamedFloatingPointLiterals };

    /// <summary>MasterPattern を保存する。</summary>
    /// <param name="compression">本体の zlib 圧縮レベル。float の強度はあまり縮まないので既定は Fastest</param>
    /// <param name="backscatteredElectrons">一緒に保存する MC の後方散乱電子 (null なら保存しない)。
    ///   入れるときは info の McSampleTiltRad / McEnergyThresholdKeV / McIncidentCount も埋めること</param>
    public static void Save(MasterPattern masterPattern, string path, MasterPatternFileInfo info,
        EbsdBackscatteredElectron[] backscatteredElectrons = null, CompressionLevel compression = CompressionLevel.Fastest)
    {
        ArgumentNullException.ThrowIfNull(masterPattern);
        ArgumentNullException.ThrowIfNull(path);
        info ??= new MasterPatternFileInfo();
        info.FormatVersion = FormatVersion;
        info.BackscatteredElectronCount = backscatteredElectrons?.Length ?? 0;
        if (string.IsNullOrEmpty(info.CreatedUtc)) info.CreatedUtc = DateTime.UtcNow.ToString("o", CultureInfo.InvariantCulture);

        //途中で失敗しても既存のファイルを壊さないよう、一時ファイルへ書いてから置き換える
        var tmp = path + ".tmp";
        using (var fs = new FileStream(tmp, FileMode.Create, FileAccess.Write, FileShare.None, 1 << 20))
        {
            // var header = JsonSerializer.SerializeToUtf8Bytes(info); // 261028Cl 変更前: NaN (電子なしのときの Mc* の既定) で例外
            var header = JsonSerializer.SerializeToUtf8Bytes(info, HeaderJsonOptions); // 261028Cl 変更
            fs.Write(Magic);
            Span<byte> i32 = stackalloc byte[4];
            BitConverter.TryWriteBytes(i32, FormatVersion); fs.Write(i32);
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
        return ReadHeader(fs);
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
        var info = ReadHeader(fs);
        using var z = new ZLibStream(fs, CompressionMode.Decompress);
        using var r = new BinaryReader(z, Encoding.UTF8);
        int gridSize = r.ReadInt32();
        var gridType = (MasterPattern.Types)r.ReadInt32();
        int nE = r.ReadInt32(), nD = r.ReadInt32();
        if (gridSize <= 0 || nE < 0 || nD < 0 || !System.Enum.IsDefined(gridType)) throw new InvalidDataException("Corrupt MasterPattern file (bad grid header).");
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
            ReadHeader(fs);
            using var z = new ZLibStream(fs, CompressionMode.Decompress);
            using var r = new BinaryReader(z, Encoding.UTF8);
            //本体の書式は Load と同じに読む (形式を変えたら両方を直す。FormatVersion の remarks)。格子の見出しも Load と同じに検査する
            int gridSize = r.ReadInt32();
            var gridType = (MasterPattern.Types)r.ReadInt32();
            int nE = r.ReadInt32(), nD = r.ReadInt32();
            if (gridSize <= 0 || nE < 0 || nD < 0 || !System.Enum.IsDefined(gridType)) throw new InvalidDataException("Corrupt MasterPattern file (bad grid header).");
            for (int i = 0; i < nE + nD; i++) r.ReadDouble(); //エネルギー列・深さ列は読み捨てる
            //平面は配列を作らず、長さ × 4 バイトを捨て緩衝 (64 KiB、LOH に入らない) へ読み捨てる (+Z → −Z の 2 × nE·nD 枚)
            var discard = new byte[1 << 16];
            for (long k = 0, planeCount = 2L * nE * nD; k < planeCount; k++)
            {
                int len = r.ReadInt32();
                if (len < 0) throw new InvalidDataException("Corrupt MasterPattern file (negative plane length).");
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

    static MasterPatternFileInfo ReadHeader(Stream fs)
    {
        Span<byte> magic = stackalloc byte[4];
        fs.ReadExactly(magic);
        if (!magic.SequenceEqual(Magic)) throw new InvalidDataException("Not a ReciPro MasterPattern file.");
        Span<byte> i32 = stackalloc byte[4];
        fs.ReadExactly(i32); int version = BitConverter.ToInt32(i32);
        if (version < 1 || version > FormatVersion) throw new InvalidDataException($"Unsupported MasterPattern file version {version} (this build reads up to {FormatVersion}).");
        fs.ReadExactly(i32); int headerLength = BitConverter.ToInt32(i32);
        if (headerLength < 0 || headerLength > 1 << 24) throw new InvalidDataException("Corrupt MasterPattern file (bad header length).");
        var header = new byte[headerLength];
        fs.ReadExactly(header);
        // return JsonSerializer.Deserialize<MasterPatternFileInfo>(header) ?? new MasterPatternFileInfo(); // 261028Cl 変更前
        return JsonSerializer.Deserialize<MasterPatternFileInfo>(header, HeaderJsonOptions) ?? new MasterPatternFileInfo(); // 261028Cl 変更: "NaN" なども読む
    }

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
    public static string CheckCompatibility(MasterPatternFileInfo info, Crystal crystal, double beamEnergyKeV)
    {
        ArgumentNullException.ThrowIfNull(info);
        ArgumentNullException.ThrowIfNull(crystal);
        if (info.CrystalFingerprint != CrystalFingerprint(crystal))
            return $"The file was made for a different crystal (\"{info.CrystalName}\"), or the cell, atoms or B differ.";
        if (Math.Abs(info.BeamEnergyKeV - beamEnergyKeV) > 1E-6)
            return $"The file was made at {info.BeamEnergyKeV:g6} keV, not {beamEnergyKeV:g6} keV.";
        return null;
    }
}
