#region using
using System;
using System.Collections.Generic;
using System.IO;
using System.IO.Compression;
using System.Linq;
using System.Reflection;
using System.Security.Cryptography;
using System.Text;
#endregion

namespace Crystallography;

//260801Cl 追加: STEM-EDX (内殻イオン化) チャネルの型群・F(s,E0) テーブル・補間器。
//設計正本 = .project-guidance/ReciPro/ReciPro_STEM-EDX設計.md §5.1、データ契約 = tools/IonizationGen/prod/MANIFEST.md (codex 15-16巡)。
//テーブルは完全自前計算 (DHFS-KS23-semi-rel-fullrange-sym-v1)。OA2000/µSTEM のデータは一切含まれない。

#region 公開 enum / record (設計書 §5.1)

/// <summary>260801Cl 追加: イオン化殻。v1 のプロバイダが返すのは K / LTotal のみ (L1/L2/L3 は v2 で分離)。</summary>
public enum IonizationShell { K, LTotal, L1, L2, L3 }

/// <summary>260801Cl 追加: 元素×殻のチャネル指定。</summary>
public record IonizationChannelSpec(int Z, IonizationShell Shell);

/// <summary>260801Cl 追加: データ出所 (σ と形状で分離して保持する)。</summary>
public sealed record IonizationDataProvenance(string ModelId, string DatasetVersion, string Detail);

/// <summary>260801Cl 追加: 正規化イオン化形状 F(s)。F(0)=1。s の単位は nm⁻¹ (s=|G|/2)。batch 評価 (N² 内の virtual call 回避)。</summary>
public interface INormalizedIonizationShape
{
    void Evaluate(ReadOnlySpan<double> sPerNm, Span<double> values);
}

/// <summary>260801Cl 追加: run 開始時に immutable へ解決されたチャネルデータ (プロバイダ選択と範囲判定を実行中に持ち込まない)。</summary>
public sealed record IonizationData(
    IonizationChannelSpec Target,
    double EdgeEnergyKeV,                 // Bote/xion edge (LTotal は最小 subshell edge)
    double TotalCrossSectionNm2,          // Bote–Salvat (LTotal は開いている subshell の合算)
    INormalizedIonizationShape Shape,     // F(0)=1
    IonizationDataProvenance CrossSectionSource,
    IonizationDataProvenance ShapeSource);

/// <summary>260801Cl 追加: 物理シグナル量。深さ分解自己吸収 (v3) の伏線として v1 から区別する (設計書 §5.5)。</summary>
public enum SignalQuantity { IonizationVacanciesGenerated, XrayPhotonsGenerated, XrayPhotonsSelfAbsorbed, XrayPhotonsDetected }

/// <summary>260801Cl 追加: モデル上の規格化状態 (表示正規化と混同しない)。</summary>
public enum SignalNormalization { ModelAbsoluteNotAudited, PerIncidentElectron }

/// <summary>260801Cl 追加: 表示正規化 (GUI 専用)。</summary>
public enum DisplayNormalization { PerMaximum, Absolute }

/// <summary>260801Cl 追加: RunSTEM への EDX 要求。v0a は 1 チャネルのみ明示受理 (codex 16巡)。
/// HermitianTolerance は ±q 非 Hermitian 残差 (相対) の許容値。超過時は対称化せず hard fail (設計書 §3.4)。</summary>
public sealed record StemIonizationRequest(IonizationChannelSpec Channel, double HermitianTolerance = 0.01);

/// <summary>260801Cl 追加: STEM-EDX 結果 (v0a 内部形。v1 で StemSimulationResult へ統合予定)。
/// ResultSTEM と同一 run の worker 終端で同時公開される (RunId で対応検証可)。</summary>
public sealed class StemEdxResult
{
    public long RunId;
    public IonizationChannelSpec Channel;
    public IonizationData Data;
    public SignalQuantity Quantity = SignalQuantity.IonizationVacanciesGenerated;
    public SignalNormalization Normalization = SignalNormalization.ModelAbsoluteNotAudited;
    /// <summary>Image[thickness][defocus][pixel]</summary>
    public double[][][] Image;
    /// <summary>±q 対称化前の非 Hermitian 残差最大値 (相対)</summary>
    public double HermitianResidualMax;
    /// <summary>q=0 の虚部残差最大値 (相対)</summary>
    public double QZeroImagMax;
    /// <summary>clamp 前の最小画素値 (負値診断)</summary>
    public double MinPixelBeforeClamp;
    /// <summary>形状評価が s>4 Å⁻¹ の tail 外挿を使ったか (診断フラグ、silent extrapolation 禁止の契約)</summary>
    public bool UsedTailExtrapolation;
}

#endregion

#region scipy 互換 PCHIP (260801Cl 追加)

/// <summary>260801Cl 追加: scipy.interpolate.PchipInterpolator 互換の単調 3 次 Hermite 補間。
/// 導関数は Fritsch–Carlson 加重調和平均 + scipy 流エッジ処理 (_edge_case)。評価は PPoly と同じ
/// 局所 power 基底 Horner。Python golden vector (tools/IonizationGen/build/golden_v1.json) とロックする。</summary>
public static class Pchip
{
    /// <summary>ノード導関数 (scipy _find_derivatives 互換)。x は狭義単調増加、n≥2。</summary>
    public static double[] Derivatives(double[] x, double[] y)
    {
        int n = x.Length;
        var d = new double[n];
        if (n == 2)
        {
            d[0] = d[1] = (y[1] - y[0]) / (x[1] - x[0]);
            return d;
        }
        var h = new double[n - 1];
        var m = new double[n - 1];
        for (int k = 0; k < n - 1; k++)
        {
            h[k] = x[k + 1] - x[k];
            m[k] = (y[k + 1] - y[k]) / h[k];
        }
        for (int k = 1; k < n - 1; k++)
        {
            if (Math.Sign(m[k]) != Math.Sign(m[k - 1]) || m[k] == 0 || m[k - 1] == 0)
                d[k] = 0;
            else
            {
                double w1 = 2 * h[k] + h[k - 1], w2 = h[k] + 2 * h[k - 1];
                var whmean = (w1 / (w1 + w2)) / m[k - 1] + (w2 / (w1 + w2)) / m[k];
                d[k] = 1.0 / whmean;
            }
        }
        d[0] = EdgeCase(h[0], h[1], m[0], m[1]);
        d[n - 1] = EdgeCase(h[n - 2], h[n - 3], m[n - 2], m[n - 3]);
        return d;
    }

    private static double EdgeCase(double h0, double h1, double m0, double m1)
    {
        var d = ((2 * h0 + h1) * m0 - h0 * m1) / (h0 + h1);
        if (Math.Sign(d) != Math.Sign(m0)) return 0.0;
        if (Math.Sign(m0) != Math.Sign(m1) && Math.Abs(d) > 3.0 * Math.Abs(m0)) return 3.0 * m0;
        return d;
    }

    /// <summary>1 点評価。範囲外は端区間の 3 次式で外挿 (scipy extrapolate=True 相当)。</summary>
    public static double Evaluate(double[] x, double[] y, double[] d, double xq)
    {
        int n = x.Length;
        int i = Array.BinarySearch(x, xq);
        if (i < 0) i = ~i - 1;
        i = Math.Clamp(i, 0, n - 2);
        double hh = x[i + 1] - x[i], slope = (y[i + 1] - y[i]) / hh;
        double c0 = (d[i] + d[i + 1] - 2 * slope) / (hh * hh);
        double c1 = (3 * slope - 2 * d[i] - d[i + 1]) / hh;
        double s = xq - x[i];
        return ((c0 * s + c1) * s + d[i]) * s + y[i];
    }
}

#endregion

#region Bote–Salvat 断面積 (260801Cl 追加)

/// <summary>260801Cl 追加: Bote–Salvat 2008 電子衝撃イオン化総断面積 (K/L/M subshell, Z=1–99)。
/// 移植元 = usnistgov/BoteSalvatICX.jl (Unlicense) / xion.f (Bote, Salvat, Jablonski, Powell, ADNDT 95 (2009) 871)。
/// 係数は埋め込みリソース Crystallography.BoteSalvat.bin (tools/IonizationGen/pack_resource.py が bote_full.json から生成)。
/// Python 参照実装 = tools/IonizationGen/botesalvat.py (golden vector で照合)。</summary>
public static class BoteSalvat
{
    private const string ResourceName = "Crystallography.BoteSalvat.bin"; // csproj の LogicalName と一致させること
    private const int Magic = 0x45544F42; // "BOTE"
    private const double Rev = 5.10998918e5;   // 電子静止エネルギー [eV] (xion.f と同値)
    private const double A0Cm = 5.291772108e-9; // Bohr 半径 [cm]

    private sealed class Element
    {
        public double[] Be, Anlj, G, EdgeEv, A; // G は [nss*4]、A は [nss*5] row-major
    }

    private static volatile Element[] _elements; // [z-1]、非 null になった時点で全構築済み (NistElastic と同じ volatile+lock 公開)
    private static readonly object _sync = new();

    private static Element[] Load()
    {
        var el = _elements;
        if (el is not null) return el;
        lock (_sync)
        {
            if (_elements is not null) return _elements;
            using var stream = Assembly.GetExecutingAssembly().GetManifestResourceStream(ResourceName)
                ?? throw new InvalidOperationException($"Embedded resource not found: {ResourceName}");
            using var reader = new BinaryReader(stream);
            if (reader.ReadInt32() != Magic) throw new InvalidDataException("BoteSalvat.bin: bad magic");
            reader.ReadInt32(); // formatVersion
            var codec = reader.ReadInt32();
            if (codec != 1) throw new InvalidDataException($"BoteSalvat.bin: unknown codec {codec}");
            ReadString(reader); // source_ref
            ReadString(reader); // packer
            var sha = reader.ReadBytes(32);
            var compLen = reader.ReadInt32();
            var comp = reader.ReadBytes(compLen);
            using var ms = new MemoryStream(comp, writable: false);
            using var br = new BrotliStream(ms, CompressionMode.Decompress);
            using var payload = new MemoryStream();
            br.CopyTo(payload);
            var raw = payload.ToArray();
            if (!SHA256.HashData(raw).AsSpan().SequenceEqual(sha))
                throw new InvalidDataException("BoteSalvat.bin: payload SHA-256 mismatch");
            using var pr = new BinaryReader(new MemoryStream(raw, writable: false));
            var zCount = pr.ReadInt32();
            var arr = new Element[100];
            for (int i = 0; i < zCount; i++)
            {
                int z = pr.ReadInt32(), nss = pr.ReadInt32();
                var e = new Element
                {
                    Be = ReadDoubles(pr, nss),
                    Anlj = ReadDoubles(pr, nss),
                    G = ReadDoubles(pr, nss * 4),
                    EdgeEv = ReadDoubles(pr, nss),
                    A = ReadDoubles(pr, nss * 5),
                };
                arr[z] = e;
            }
            _elements = arr;
            return arr;
        }
    }

    internal static string ReadString(BinaryReader r)
    {
        var len = r.ReadInt32();
        return Encoding.UTF8.GetString(r.ReadBytes(len));
    }

    private static double[] ReadDoubles(BinaryReader r, int count)
    {
        var a = new double[count];
        for (int i = 0; i < count; i++) a[i] = r.ReadDouble();
        return a;
    }

    private static Element Get(int z)
        => (uint)z <= 99 && Load()[z] is Element e ? e
           : throw new ArgumentOutOfRangeException(nameof(z), z, "Bote–Salvat coefficients cover Z=1–99");

    /// <summary>Z の収録サブシェル数 (1〜9)。index 順 = K, L1, L2, L3, M1..M5。</summary>
    public static int SubshellCount(int z) => Get(z).EdgeEv.Length;

    /// <summary>吸収端エネルギー [eV]。subshell は 1 始まり (1=K, 2=L1, 3=L2, 4=L3, 5..9=M1..M5)。</summary>
    public static double EdgeEv(int z, int subshell) => Get(z).EdgeEv[CheckSubshell(z, subshell)];

    private static int CheckSubshell(int z, int subshell)
        => subshell >= 1 && subshell <= SubshellCount(z) ? subshell - 1
           : throw new ArgumentOutOfRangeException(nameof(subshell), subshell, $"Z={z} has {SubshellCount(z)} subshells");

    /// <summary>イオン化断面積 [cm²]。演算順は botesalvat.py sigma_cm2 と厳密一致 (golden vector 照合)。</summary>
    public static double SigmaCm2(int z, int subshell, double energyEv, double edgeEvOverride = double.NaN)
    {
        var el = Get(z);
        int ss = CheckSubshell(z, subshell);
        var edge = double.IsNaN(edgeEvOverride) ? el.EdgeEv[ss] : edgeEvOverride;
        var overv = energyEv / edge;
        if (overv <= 1.0) return 0.0;
        double xione;
        if (overv <= 16.0)
        {
            double a1 = el.A[ss * 5], a2 = el.A[ss * 5 + 1], a3 = el.A[ss * 5 + 2], a4 = el.A[ss * 5 + 3], a5 = el.A[ss * 5 + 4];
            var opu = 1.0 / (1.0 + overv);
            var ffitlo = a1 + a2 * overv + opu * (a3 + opu * opu * (a4 + opu * opu * a5));
            var r = ffitlo / overv;
            xione = (overv - 1.0) * (r * r);
        }
        else
        {
            var beta2 = (energyEv * (energyEv + 2.0 * Rev)) / ((energyEv + Rev) * (energyEv + Rev));
            var x = Math.Sqrt(energyEv * (energyEv + 2.0 * Rev)) / Rev;
            double g1 = el.G[ss * 4], g2 = el.G[ss * 4 + 1], g3 = el.G[ss * 4 + 2], g4 = el.G[ss * 4 + 3];
            var ffitup = (2.0 * Math.Log(x) - beta2) * (1.0 + g1 / x) + g2
                + g3 * Math.Sqrt(Rev / (energyEv + Rev)) + g4 / x;
            var factr = el.Anlj[ss] / beta2;
            xione = ((factr * overv) / (overv + el.Be[ss])) * ffitup;
        }
        return 4.0 * Math.PI * (A0Cm * A0Cm) * xione;
    }

    /// <summary>イオン化断面積 [nm²]。</summary>
    public static double SigmaNm2(int z, int subshell, double energyEv, double edgeEvOverride = double.NaN)
        => SigmaCm2(z, subshell, energyEv, edgeEvOverride) * 1e14;
}

#endregion

#region F(s,E0) テーブル (260801Cl 追加)

/// <summary>260801Cl 追加: 本番 F(s,E0) テーブル (dataset 1.0.0, 127ch) のリーダー。
/// フォーマット・契約 = tools/IonizationGen/pack_resource.py ヘッダコメント + prod/MANIFEST.md。
/// NistElasticPchipResource と同じ「blob 常駐 + チャネル単位 lazy Brotli decode + volatile 公開」。</summary>
public sealed class IonizationFsTable
{
    private const string ResourceName = "Crystallography.IonizationFsE0.bin"; // csproj の LogicalName と一致させること
    private const int Magic = 0x31534649; // "IFS1"
    public const int ShellCodeK = 0, ShellCodeL1 = 1, ShellCodeL23 = 2;
    public const double SMaxAngstromInv = 4.0;

    public int Method { get; }             // 1=float32 / 2=1e-6 量子化+delta+shuffle
    public int SCount { get; }             // 81
    public double SStep { get; }           // 0.05 Å⁻¹
    public string DatasetVersion { get; }
    public string ModelId { get; }
    public string BoteRef { get; }
    public string Packer { get; }
    public byte[] PayloadSha256 { get; }

    private readonly byte[] _blob;
    private readonly int _payloadStart;
    private readonly Dictionary<(int ShellCode, int Z), (int Offset, int Length)> _index = [];
    private readonly Dictionary<(int ShellCode, int Z), IonizationChannelTable> _cache = [];
    private readonly object _sync = new();
    internal readonly double[] SGrid;

    private static volatile IonizationFsTable _default;
    private static readonly object _defaultSync = new();

    /// <summary>埋め込みリソースから構築される既定インスタンス。</summary>
    public static IonizationFsTable Default
    {
        get
        {
            var t = _default;
            if (t is not null) return t;
            lock (_defaultSync)
            {
                if (_default is null)
                {
                    using var s = Assembly.GetExecutingAssembly().GetManifestResourceStream(ResourceName)
                        ?? throw new InvalidOperationException($"Embedded resource not found: {ResourceName}");
                    _default = new IonizationFsTable(s);
                }
                return _default;
            }
        }
    }

    /// <summary>任意ストリームから構築 (方式比較・破損テスト用)。ヘッダは厳格検査し、unknown codec/method は拒否する。</summary>
    public IonizationFsTable(Stream stream)
    {
        var blob = new byte[stream.Length];
        stream.ReadExactly(blob);
        _blob = blob;
        using var reader = new BinaryReader(new MemoryStream(blob, writable: false));
        if (reader.ReadInt32() != Magic) throw new InvalidDataException("IonizationFsE0.bin: bad magic");
        var formatVersion = reader.ReadInt32();
        if (formatVersion != 1) throw new InvalidDataException($"IonizationFsE0.bin: unknown format version {formatVersion}");
        var codec = reader.ReadInt32();
        if (codec != 1) throw new InvalidDataException($"IonizationFsE0.bin: unknown codec {codec}");
        Method = reader.ReadInt32();
        if (Method is not (1 or 2)) throw new InvalidDataException($"IonizationFsE0.bin: unknown method {Method}");
        SCount = reader.ReadInt32();
        SStep = reader.ReadDouble();
        if (SCount != 81 || SStep != 0.05) throw new InvalidDataException("IonizationFsE0.bin: unexpected s grid");
        reader.ReadInt32(); // schemaVersion
        DatasetVersion = BoteSalvat.ReadString(reader);
        ModelId = BoteSalvat.ReadString(reader);
        BoteRef = BoteSalvat.ReadString(reader);
        Packer = BoteSalvat.ReadString(reader);
        PayloadSha256 = reader.ReadBytes(32);
        reader.ReadBytes(32); // sourceSha256 (記録用)
        var channelCount = reader.ReadInt32();
        if (channelCount is <= 0 or > 10000) throw new InvalidDataException("IonizationFsE0.bin: bad channel count");
        long payloadLen = 0;
        for (int i = 0; i < channelCount; i++)
        {
            int shellCode = reader.ReadInt32(), z = reader.ReadInt32(), offset = reader.ReadInt32(), length = reader.ReadInt32();
            if (shellCode is < ShellCodeK or > ShellCodeL23 || length <= 0 || offset != payloadLen)
                throw new InvalidDataException("IonizationFsE0.bin: bad index entry"); // offset 連続 = 重複/オーバーラップ拒否
            if (!_index.TryAdd((shellCode, z), (offset, length)))
                throw new InvalidDataException($"IonizationFsE0.bin: duplicate channel ({shellCode},{z})");
            payloadLen += length;
        }
        _payloadStart = (int)reader.BaseStream.Position;
        if (_payloadStart + payloadLen != blob.Length) throw new InvalidDataException("IonizationFsE0.bin: payload length mismatch");
        SGrid = new double[SCount];
        for (int i = 0; i < SCount; i++) SGrid[i] = i * SStep;
    }

    public bool Contains(int shellCode, int z) => _index.ContainsKey((shellCode, z));

    public IonizationChannelTable GetChannel(int shellCode, int z)
    {
        lock (_sync)
        {
            if (_cache.TryGetValue((shellCode, z), out var hit)) return hit;
            if (!_index.TryGetValue((shellCode, z), out var entry))
                throw new NotSupportedException($"Ionization table has no channel shellCode={shellCode}, Z={z} (K: Z=6–50, L1/L23: Z=20–60)");
            var table = Decode(entry, shellCode, z);
            _cache.Add((shellCode, z), table); // lock 内構築 = ExecutionAndPublication (半初期化を見せない)
            return table;
        }
    }

    /// <summary>全チャネルを展開して canonical payload の SHA-256 を検証 (ハーネス用。runtime 起動時には呼ばない)。</summary>
    public bool VerifyPayloadHash()
    {
        using var sha = IncrementalHash.CreateHash(HashAlgorithmName.SHA256);
        foreach (var kv in _index.OrderBy(e => e.Value.Offset))
            sha.AppendData(DecompressBlob(kv.Value));
        return sha.GetHashAndReset().AsSpan().SequenceEqual(PayloadSha256);
    }

    private byte[] DecompressBlob((int Offset, int Length) entry)
    {
        using var ms = new MemoryStream(_blob, _payloadStart + entry.Offset, entry.Length, writable: false);
        using var br = new BrotliStream(ms, CompressionMode.Decompress);
        using var outMs = new MemoryStream();
        br.CopyTo(outMs);
        return outMs.ToArray();
    }

    private IonizationChannelTable Decode((int Offset, int Length) entry, int expectShell, int expectZ)
    {
        var raw = DecompressBlob(entry);
        using var r = new BinaryReader(new MemoryStream(raw, writable: false));
        int z = r.ReadInt32(), shellCode = r.ReadInt32();
        if (z != expectZ || shellCode != expectShell) throw new InvalidDataException("IonizationFsE0.bin: channel blob/index mismatch");
        var eth = r.ReadDouble();
        var rowCount = r.ReadInt32();
        if (rowCount is < 2 or > 1000) throw new InvalidDataException("IonizationFsE0.bin: bad row count");
        var e0 = new double[rowCount];
        for (int i = 0; i < rowCount; i++) e0[i] = r.ReadDouble();
        var u = new double[rowCount];
        for (int i = 0; i < rowCount; i++) u[i] = r.ReadDouble();
        var tailFlag = r.ReadBytes(rowCount);
        var tailA = new double[rowCount];
        for (int i = 0; i < rowCount; i++) tailA[i] = r.ReadDouble();
        var tailB = new double[rowCount];
        for (int i = 0; i < rowCount; i++) tailB[i] = r.ReadDouble();
        var f = new double[rowCount][];
        if (Method == 1)
        {
            for (int i = 0; i < rowCount; i++)
            {
                var row = new double[SCount];
                for (int j = 0; j < SCount; j++) row[j] = r.ReadSingle();
                f[i] = row;
            }
        }
        else // method 2: int32 量子化 + 行内 s 方向 delta + byte-plane shuffle (out[p*n+i] = raw[i*4+p])
        {
            int n = rowCount * SCount;
            var shuffled = r.ReadBytes(n * 4);
            if (shuffled.Length != n * 4) throw new InvalidDataException("IonizationFsE0.bin: truncated F block");
            var q = new int[n];
            for (int i = 0; i < n; i++)
                q[i] = shuffled[i] | (shuffled[n + i] << 8) | (shuffled[2 * n + i] << 16) | (shuffled[3 * n + i] << 24);
            for (int i = 0; i < rowCount; i++)
            {
                var row = new double[SCount];
                long acc = 0;
                for (int j = 0; j < SCount; j++)
                {
                    acc += q[i * SCount + j];
                    row[j] = acc * 1e-6;
                }
                f[i] = row;
            }
        }
        if (r.BaseStream.Position != raw.Length) throw new InvalidDataException("IonizationFsE0.bin: blob length mismatch");
        // 構造検証: F(0)=1 (method2 は量子化後も 1e6*1e-6)、tail ノード整合 a=F(4)·e^{4b} (量子化誤差ぶんの許容)
        var aTol = Method == 2 ? 1.1e-6 : 6e-8;
        for (int i = 0; i < rowCount; i++)
        {
            if (Math.Abs(f[i][0] - 1.0) > 1e-12) throw new InvalidDataException($"IonizationFsE0.bin: F(0)≠1 (Z={z}, row {i})");
            if (tailFlag[i] != 0 && Math.Abs(tailA[i] * Math.Exp(-4.0 * tailB[i]) - f[i][SCount - 1]) > aTol)
                throw new InvalidDataException($"IonizationFsE0.bin: tail/F(4) inconsistent (Z={z}, row {i})");
        }
        return new IonizationChannelTable(this, z, shellCode, eth, e0, u, tailFlag, tailA, tailB, f);
    }
}

/// <summary>260801Cl 追加: 1 チャネル分の展開済みテーブルと E0 補間 (契約 = prod/MANIFEST.md)。</summary>
public sealed class IonizationChannelTable
{
    private readonly IonizationFsTable _owner;
    public readonly int Z;
    public readonly int ShellCode;
    public readonly double EthKeV;
    public readonly double[] E0KeV;   // 厳密昇順
    public readonly double[] U;       // serialized row.u (4桁丸め) = 補間ノット契約値
    private readonly double[] _x;     // ln(u-1)
    private readonly byte[] _tailFlag;
    private readonly double[] _tailA, _tailB;
    private readonly double[][] _f;   // [row][81]
    private readonly List<(int Lo, int Hi)> _tailRuns = [];

    internal IonizationChannelTable(IonizationFsTable owner, int z, int shellCode, double eth,
        double[] e0, double[] u, byte[] tailFlag, double[] tailA, double[] tailB, double[][] f)
    {
        _owner = owner; Z = z; ShellCode = shellCode; EthKeV = eth;
        E0KeV = e0; U = u; _tailFlag = tailFlag; _tailA = tailA; _tailB = tailB; _f = f;
        _x = new double[u.Length];
        for (int i = 0; i < u.Length; i++) _x[i] = Math.Log(u[i] - 1.0);
        for (int i = 0; i < tailFlag.Length;)
        {
            if (tailFlag[i] != 0)
            {
                int k = i;
                while (k + 1 < tailFlag.Length && tailFlag[k + 1] != 0) k++;
                _tailRuns.Add((i, k));
                i = k + 1;
            }
            else i++;
        }
    }

    public int RowCount => E0KeV.Length;

    /// <summary>全 (row, s節点) の F 総和 (golden との構造照合用の診断値)。</summary>
    public double SumF()
    {
        var sum = 0.0;
        foreach (var row in _f)
            foreach (var v in row) sum += v;
        return sum;
    }

    /// <summary>E0 [keV] における 81 節点グリッドを契約どおり補間 (各 s 節点で x=ln(u−1) PCHIP、
    /// 全行正なら lnF・非正含みは符号付き F 直接)。E0 は 30–400 keV 限定 (外挿・clamp 禁止)。</summary>
    public double[] GridAt(double e0KeV)
    {
        if (!(e0KeV >= 30.0 && e0KeV <= 400.0))
            throw new ArgumentOutOfRangeException(nameof(e0KeV), e0KeV, "F(s,E0) table covers E0 = 30–400 keV only (no extrapolation)");
        var xq = Math.Log(e0KeV / EthKeV - 1.0);
        int rows = RowCount, sCount = _owner.SCount;
        var grid = new double[sCount];
        var col = new double[rows];
        for (int j = 0; j < sCount; j++)
        {
            var allPositive = true;
            for (int i = 0; i < rows; i++)
            {
                col[i] = _f[i][j];
                if (col[i] <= 0) allPositive = false;
            }
            if (allPositive)
            {
                for (int i = 0; i < rows; i++) col[i] = Math.Log(col[i]);
                grid[j] = Math.Exp(Pchip.Evaluate(_x, col, Pchip.Derivatives(_x, col), xq));
            }
            else
                grid[j] = Pchip.Evaluate(_x, col, Pchip.Derivatives(_x, col), xq);
        }
        grid[0] = 1.0; // s=0 は厳密 1 (契約)
        return grid;
    }

    /// <summary>s>4 tail の減衰係数 b̂(E0) (連続性アンカー方式、codex 16巡)。取得不能なら false。
    /// tail≠null の連続 E0 区間内のみで PCHIP。E0 を挟む行に null が絡む場合は不可。exact node はその行の b。</summary>
    public bool TryGetTailB(double e0KeV, out double bHat)
    {
        var hit = Array.BinarySearch(E0KeV, e0KeV);
        if (hit >= 0 && _tailFlag[hit] != 0) { bHat = _tailB[hit]; return true; }
        int i = hit >= 0 ? hit : ~hit - 1;
        i = Math.Clamp(i, 0, RowCount - 2);
        foreach (var (lo, hi) in _tailRuns)
            if (lo <= i && i + 1 <= hi && hi > lo)
            {
                var xs = _x[lo..(hi + 1)];
                var bs = _tailB[lo..(hi + 1)];
                bHat = Pchip.Evaluate(xs, bs, Pchip.Derivatives(xs, bs), Math.Log(e0KeV / EthKeV - 1.0));
                return true;
            }
        bHat = double.NaN;
        return false;
    }

    /// <summary>E0 を固定して解決した形状 (run-scoped)。</summary>
    public IonizationTableShape BuildShape(double e0KeV) => new(this, _owner.SGrid, GridAt(e0KeV), e0KeV);
}

/// <summary>260801Cl 追加: E0 解決済みの単一殻形状。s 方向は符号付き F 直接 PCHIP、
/// s>4 Å⁻¹ は連続性アンカー tail F(4)e^{−b̂(s−4)} (b̂ 不能なら hard fail)。入力 s の単位は nm⁻¹。</summary>
public sealed class IonizationTableShape : INormalizedIonizationShape
{
    private readonly double[] _sGrid, _grid, _deriv;
    private readonly double _f4, _bHat;
    private readonly bool _tailAvailable;
    private readonly int _z, _shellCode;
    private bool _usedTail;

    internal IonizationTableShape(IonizationChannelTable table, double[] sGrid, double[] grid, double e0KeV)
    {
        _sGrid = sGrid; _grid = grid;
        _deriv = Pchip.Derivatives(sGrid, grid);
        _f4 = grid[^1];
        _tailAvailable = table.TryGetTailB(e0KeV, out _bHat);
        _z = table.Z; _shellCode = table.ShellCode;
    }

    /// <summary>s>4 Å⁻¹ の tail 外挿を使ったか (診断)。</summary>
    public bool UsedTailExtrapolation => _usedTail;

    public void Evaluate(ReadOnlySpan<double> sPerNm, Span<double> values)
    {
        for (int k = 0; k < sPerNm.Length; k++)
        {
            var sA = sPerNm[k] * 0.1; // nm⁻¹ → Å⁻¹
            if (sA == 0.0)
                values[k] = 1.0;
            else if (sA <= IonizationFsTable.SMaxAngstromInv)
                values[k] = Pchip.Evaluate(_sGrid, _grid, _deriv, sA);
            else if (_tailAvailable)
            {
                values[k] = _f4 * Math.Exp(-_bHat * (sA - 4.0));
                _usedTail = true;
            }
            else
                throw new InvalidOperationException(
                    $"F(s) tail unavailable for s={sA:f3} Å⁻¹ > 4 (Z={_z}, shellCode={_shellCode}: null-tail rows bracket this E0). Reduce gMax or refuse the channel.");
        }
    }
}

/// <summary>260801Cl 追加: LTotal 合成形状 F_L = [σ_L1·F_L1 + (σ_L2+σ_L3)·F_L23]/Σσ (実行時 Bote 重み、MANIFEST 契約)。</summary>
public sealed class IonizationLTotalShape : INormalizedIonizationShape
{
    private readonly IonizationTableShape _l1, _l23;
    private readonly double _w1, _w23; // σ 重み (正規化済み)

    internal IonizationLTotalShape(IonizationTableShape l1, double sigmaL1, IonizationTableShape l23, double sigmaL23)
    {
        var total = sigmaL1 + sigmaL23;
        _l1 = l1; _l23 = l23;
        _w1 = sigmaL1 / total;
        _w23 = sigmaL23 / total;
    }

    public bool UsedTailExtrapolation => (_l1?.UsedTailExtrapolation ?? false) || (_l23?.UsedTailExtrapolation ?? false);

    public void Evaluate(ReadOnlySpan<double> sPerNm, Span<double> values)
    {
        Span<double> tmp = sPerNm.Length <= 256 ? stackalloc double[sPerNm.Length] : new double[sPerNm.Length];
        if (_w1 > 0) _l1.Evaluate(sPerNm, values);
        else values.Clear();
        for (int k = 0; k < values.Length; k++) values[k] *= _w1;
        if (_w23 > 0)
        {
            _l23.Evaluate(sPerNm, tmp);
            for (int k = 0; k < values.Length; k++) values[k] += _w23 * tmp[k];
        }
    }
}

#endregion

#region プロバイダ (260801Cl 追加)

/// <summary>260801Cl 追加: チャネル指定 → 解決済み IonizationData。run 開始時に 1 回だけ呼び、実行中は immutable を使う (設計書 §5.1)。</summary>
public static class IonizationDataProvider
{
    /// <summary>解決。E0 範囲外 (30–400 keV 以外) は ArgumentOutOfRangeException、
    /// 未収録 Z/殻は NotSupportedException。v1 で解決可能な殻は K / LTotal のみ。</summary>
    public static IonizationData Resolve(IonizationChannelSpec spec, double e0KeV, IonizationFsTable table = null)
    {
        ArgumentNullException.ThrowIfNull(spec);
        if (!(e0KeV >= 30.0 && e0KeV <= 400.0))
            throw new ArgumentOutOfRangeException(nameof(e0KeV), e0KeV, "STEM-EDX supports E0 = 30–400 keV only (F table range, no extrapolation)");
        table ??= IonizationFsTable.Default;
        var shapeProvenance = new IonizationDataProvenance(table.ModelId, table.DatasetVersion, "self-generated DHFS tables (tools/IonizationGen prod)");
        var sigmaProvenance = new IonizationDataProvenance("Bote-Salvat-2008", "xion.f/ADNDT95", table.BoteRef);
        var eV = e0KeV * 1e3;
        switch (spec.Shell)
        {
            case IonizationShell.K:
                {
                    var ch = table.GetChannel(IonizationFsTable.ShellCodeK, spec.Z);
                    var sigma = BoteSalvat.SigmaNm2(spec.Z, 1, eV);
                    return new IonizationData(spec, ch.EthKeV, sigma, ch.BuildShape(e0KeV), sigmaProvenance, shapeProvenance);
                }
            case IonizationShell.LTotal:
                {
                    var chL1 = table.GetChannel(IonizationFsTable.ShellCodeL1, spec.Z);
                    var chL23 = table.GetChannel(IonizationFsTable.ShellCodeL23, spec.Z);
                    // σ は各サブシェル自身の edge で計算 (MANIFEST 契約)。閉じている subshell は 0 で自然に落ちる
                    double s1 = BoteSalvat.SigmaNm2(spec.Z, 2, eV), s2 = BoteSalvat.SigmaNm2(spec.Z, 3, eV), s3 = BoteSalvat.SigmaNm2(spec.Z, 4, eV);
                    var total = s1 + s2 + s3;
                    if (total <= 0)
                        throw new NotSupportedException($"Z={spec.Z} LTotal: all L subshells closed at E0={e0KeV} keV");
                    var shape = new IonizationLTotalShape(chL1.BuildShape(e0KeV), s1, chL23.BuildShape(e0KeV), s2 + s3);
                    var edge = Math.Min(chL1.EthKeV, chL23.EthKeV);
                    return new IonizationData(spec, edge, total, shape, sigmaProvenance, shapeProvenance);
                }
            default:
                throw new NotSupportedException($"IonizationShell.{spec.Shell} is not available in v1 (use K or LTotal)");
        }
    }
}

#endregion
