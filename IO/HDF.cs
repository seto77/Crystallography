using PureHDF;
using System;
using System.Collections.Generic;
using System.Data;
using System.Linq;

namespace Crystallography;
//public class HDF// (260715Ch) 旧: OpenRead した NativeFile の所有権を保持・解放できなかった
public class HDF : IDisposable// (260715Ch)
{
    private readonly IDisposable file;// (260715Ch) H5DatasetAdv の遅延 Read が終わるまで NativeFile を保持する

    /// <summary>データ型が不明なときに使うダミークラス</summary>
    private class Unknown() { }

    /// <summary>IH5Datasetを使いやすくするため拡張したクラス (階層構造を把握するためにPathを追加した)</summary>
    public class H5DatasetAdv
    {
        private readonly IH5Dataset dataset;
        public string Path;
        public Type DataType;

        /// <param name="dataset"></param>
        /// <param name="path"></param>
        public H5DatasetAdv(IH5Dataset dataset, string path)
        {
            this.dataset = dataset;
            Path = path;
            DataType = Type.Class switch
            {
                H5DataTypeClass.String => typeof(string),
                H5DataTypeClass.VariableLength => typeof(string),
                H5DataTypeClass.FixedPoint => Type.Size switch
                {
                    1 => Type.FixedPoint.IsSigned ? typeof(sbyte) : typeof(byte),
                    2 => Type.FixedPoint.IsSigned ? typeof(short) : typeof(ushort),
                    4 => Type.FixedPoint.IsSigned ? typeof(int) : typeof(uint),
                    8 => Type.FixedPoint.IsSigned ? typeof(long) : typeof(ulong),
                    _ => typeof(Unknown),
                },
                H5DataTypeClass.FloatingPoint => Type.Size switch
                {
                    4 => typeof(float),
                    8 => typeof(double),
                    _ => typeof(Unknown),
                },
                _ => typeof(Unknown),
            };
        }

        /// <summary>Datasetの名前</summary>
        public string Name => dataset.Name;

        /// <summary>Datasetの次元など</summary>
        public IH5Dataspace Space => dataset.Space;

        /// <summary>データセットのデータ型など</summary>
        public IH5DataType Type => dataset.Type;

        /// <summary>データ読み込み</summary>
        /// <typeparam name="T"></typeparam>
        /// <returns></returns>
        public T Read<T>() => dataset.Read<T>();

        /// <summary>データ読み込み (多次元配列の一部を読み込むときに使用)</summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="rank"></param>
        /// <param name="starts"></param>
        /// <param name="blocks"></param>
        /// <returns></returns>
        public T Read<T>(ulong[] starts, ulong[] blocks) => dataset.Read<T>(new PureHDF.Selections.HyperslabSelection(starts.Length, starts, blocks), null, null);

        /// <summary>
        /// データ読み込み (多次元配列の一部を読み込むときに使用)
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="rank"></param>
        /// <param name="starts"></param>
        /// <param name="blocks"></param>
        /// <returns></returns>
        public T Read<T>(int[] starts, int[] blocks) => Read<T>([.. starts.Select(s => (ulong)s)], [.. blocks.Select(b => (ulong)b)]);


        /// <summary>数値配列(1次元)データを読み込み、double配列に変換して返す</summary>
        /// <returns></returns>
        public double[] ReadAsDoubleArray()
        {
            if (DataType == typeof(sbyte))
                return [.. Read<sbyte[]>().Select(e => (double)e)];
            else if (DataType == typeof(byte))
                return [.. Read<byte[]>().Select(e => (double)e)];
            else if (DataType == typeof(short))
                return [.. Read<short[]>().Select(e => (double)e)];
            else if (DataType == typeof(ushort))
                return [.. Read<ushort[]>().Select(e => (double)e)];
            else if (DataType == typeof(int))
                return [.. Read<int[]>().Select(e => (double)e)];
            else if (DataType == typeof(uint))
                return [.. Read<uint[]>().Select(e => (double)e)];
            else if (DataType == typeof(long))
                return [.. Read<long[]>().Select(e => (double)e)];
            else if (DataType == typeof(ulong))
                return [.. Read<ulong[]>().Select(e => (double)e)];
            else if (DataType == typeof(float))
                return [.. Read<float[]>().Select(e => (double)e)];
            else if (DataType == typeof(double))
                return [.. Read<double[]>()];
            else
                return null;
        }

        /// <summary>数値配列(2次元以上)データを読み込み、double配列(1次元)に変換して返す</summary>
        /// <returns></returns>
        public double[] ReadAsDoubleArray(ulong[] starts, ulong[] blocks)
        {
            if (DataType == typeof(sbyte))
                return [.. Read<sbyte[]>(starts,blocks).Select(e => (double)e)];
            else if (DataType == typeof(byte))
                return [.. Read<byte[]>(starts, blocks).Select(e => (double)e)];
            else if (DataType == typeof(short))
                return [.. Read<short[]>(starts, blocks).Select(e => (double)e)];
            else if (DataType == typeof(ushort))
                return [.. Read<ushort[]>(starts, blocks).Select(e => (double)e)];
            else if (DataType == typeof(int))
                return [.. Read<int[]>(starts, blocks).Select(e => (double)e)];
            else if (DataType == typeof(uint))
                return [.. Read<uint[]>(starts, blocks).Select(e => (double)e)];
            else if (DataType == typeof(long))
                return [.. Read<long[]>(starts, blocks).Select(e => (double)e)];
            else if (DataType == typeof(ulong))
                return [.. Read<ulong[]>(starts, blocks).Select(e => (double)e)];
            else if (DataType == typeof(float))
                return [.. Read<float[]>(starts, blocks).Select(e => (double)e)];
            else if (DataType == typeof(double))
                return [.. Read<double[]>(starts, blocks)];
            else
                return null;
        }

        /// <summary>数値配列(2次元以上)データを読み込み、double配列(1次元)に変換して返す</summary>
        /// <returns></returns>
        public double[] ReadAsDoubleArray(int[] starts, int[] blocks)
            =>ReadAsDoubleArray([.. starts.Select(s => (ulong)s)], [.. blocks.Select(b => (ulong)b)]);

        /// <summary>データ読み込み (文字列専用)</summary>
        /// <returns></returns>
        public string ReadStr() => dataset.Read<string>();

        /// <summary>データ(形式は問わないが、配列はNG)を読み込み文字列として返す</summary>
        /// <returns></returns>
        public string ReadAsStr()
        {
            if (DataType == typeof(string))
                return Read<string>();
            else if (DataType == typeof(sbyte))
                return Read<sbyte>().ToString();
            else if (DataType == typeof(byte))
                return Read<byte>().ToString();
            else if (DataType == typeof(short))
                return Read<short>().ToString();
            else if (DataType == typeof(ushort))
                return Read<ushort>().ToString();
            else if (DataType == typeof(int))
                return Read<int>().ToString();
            else if (DataType == typeof(uint))
                return Read<uint>().ToString();
            else if (DataType == typeof(long))
                return Read<long>().ToString();
            else if (DataType == typeof(ulong))
                return Read<ulong>().ToString();
            else if (DataType == typeof(float))
                return Read<float>().ToString();
            else if (DataType == typeof(double))
                return Read<double>().ToString();
            else
                return "";
        }
    }

    /// <summary>HDFファイルに含まれるDatasetのリスト</summary>
    public List<H5DatasetAdv> Datasets = [];

    /// <summary>pathで指定したデータセットを取得する</summary>
    /// <param name="path"></param>
    /// <returns></returns>
    public H5DatasetAdv GetDataset(string path)
    {
        var dataset = Datasets.FirstOrDefault(d => d.Path == path);
        return (dataset!=null && path == dataset.Path) ? dataset : null;
    }

    /// <summary>HDFクラスのコンストラクタ</summary>
    /// <param name="filename"></param>
    public HDF(string filename)
    {
        // 全てのデータセットを再帰的に取得するためのローカル関数
        void addDatasetRecursively (IH5Group group, string path)
        {
            List<IH5Dataset> datasets = [];
            List<IH5Group> groups = [];
            foreach (var obj in group.Children())
            {
                if (obj is IH5Group childGroup)
                    groups.Add(childGroup);
                else if (obj is IH5Dataset childDataset)
                    datasets.Add(childDataset);
            }
            Datasets.AddRange(datasets.Select(d => new H5DatasetAdv(d, $"{path}/{d.Name}")));
            
            foreach (var childGroup in groups)
                addDatasetRecursively(childGroup, $"{path}/{childGroup.Name}");
        }

        //var file = H5File.OpenRead(filename);// (260715Ch) 旧: ローカル変数のまま破棄されず、ファイルハンドルをリークした
        var openedFile = H5File.OpenRead(filename);// (260715Ch)
        file = openedFile;
        try
        {
            addDatasetRecursively(openedFile.Group("/"), "");
        }
        catch
        {
            openedFile.Dispose();// (260715Ch) 列挙中に例外が出たコンストラクタ経路でも所有ハンドルを解放する
            throw;
        }
    }

    /// <summary>HDF5 ファイルと、そこから取得した遅延読み込み Dataset の所有リソースを解放する。</summary>
    public void Dispose() => file.Dispose();// (260715Ch)

}
