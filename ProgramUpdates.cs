using System;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Net.Http;
using System.Reflection;
using System.Windows.Forms;

namespace Crystallography;

public static class ProgramUpdates
{
    private static readonly string UserAppDataPath = new DirectoryInfo(Application.UserAppDataPath).Parent.FullName + @"\";

    //260317Cl HttpClientはstaticで再利用するのがベストプラクティス
    private static readonly HttpClient httpClient = new();

    //260613Cl installerAssetName 追加 (デフォルト引数、既存呼び出しは無変更): アーキ別インストーラ資産名。
    //          空なら従来どおり {software}Setup.msi。ReciPro の arm64 ビルドは "ReciProSetup-arm64.msi" を渡す。
    //260805Cl portable 追加 (デフォルト引数): Portable ZIP 版は MSI 自動更新が成り立たないため、新版があれば
    //          最新リリースページの URL を返す (Path 空 = ダウンロード対象なしの印。呼び出し側はブラウザで URL を開く)。
    //旧シグネチャ: public static (string Title, string Message, bool NeedUpdate, string URL, string Path) Check(string software, string version)
    //旧シグネチャ: public static (string Title, string Message, bool NeedUpdate, string URL, string Path) Check(string software, string version, string installerAssetName = "")
    public static (string Title, string Message, bool NeedUpdate, string URL, string Path) Check(string software, string version, string installerAssetName = "", bool portable = false)
    {
        try
        {
            var ver = httpClient.GetByteArrayAsync($"https://raw.githubusercontent.com/seto77/{software}/master/{software}/Version.cs").Result;

            //V上手くダウンロードできなかった場合
            if (ver == null || ver.Length == 0)
                return ("Error!", $"An error occurred while trying to locate the update to {software}.\r\n " +
                    "This could be caused if you do not have an active internet connection, or host server may be down. ", false, "", "");

            var temp = System.Text.Encoding.UTF8.GetString(ver).Split(new[] { '\r', '\n' });
            var newVersion = temp.First(s => s.Contains(" ver", StringComparison.Ordinal));
            newVersion = newVersion.Substring(newVersion.IndexOf("ver") + 3, 5);

            if (double.Parse(newVersion, System.Globalization.CultureInfo.InvariantCulture) <=
                double.Parse(version.Substring(3, 5), System.Globalization.CultureInfo.InvariantCulture)) //260715Ch 小数点記号がカンマのカルチャでは解析に失敗するため、バージョン表記は常にドット区切り (InvariantCulture) として比較する
                return ("Update checked!", $"You are running the latest version of {software}. Thank you!", false, "", "");
            else
            {
                //260805Cl 追加: Portable ZIP 版はインストーラを扱わず、最新リリースページをブラウザで開く導線のみ提供する
                if (portable)
                    return ("Update checked!", $"New version {newVersion} is available.\r\n" +
                        $"If you press 'Yes', the release page of the latest {software} will open in your browser.", true,
                        $"https://github.com/seto77/{software}/releases/latest", "");

                //260613Cl アーキ別インストーラ対応 + 資産存在チェック。
                //          arm64 MSI は release 作成後に実機 smoke 合格を待って後付け添付されるため、
                //          「新版は出ているが当該アーキのインストーラが未添付」の時間窓がある。
                //          DL を 404 で失敗させる前に GET (ヘッダのみ) で存在を確認して案内する。
                var asset = string.IsNullOrEmpty(installerAssetName) ? software + "Setup.msi" : installerAssetName;
                var url = $"https://github.com/seto77/{software}/releases/download/v.{newVersion}/{asset}"; //260715Ch 更新インストーラを最初のリクエストから TLS で取得する
                using var request = new HttpRequestMessage(HttpMethod.Get, url); //260715Ch Send 後も request を確実に破棄する
                using var res = httpClient.Send(request, HttpCompletionOption.ResponseHeadersRead); //260715Ch ヘッダ確認後に response を破棄する
                if (!res.IsSuccessStatusCode)
                    return ("Update checked!", $"New version {newVersion} is available, but the installer for this architecture ({asset}) " +
                        $"has not been published yet. Please try again later, or download the latest package from\r\nhttps://github.com/seto77/{software}/releases/latest", false, "", "");
                return ($"Update checked!", $"Now, new version {newVersion} is available.\r\n" +
                     $"If you press 'Yes', the current {software} will be closed immediately and the installer of new {software} launched.", true,
                     url, UserAppDataPath + asset);
            }
        }
        catch
        {
            return ("Error!", "An error occurred while trying to locate the update to " + software + ".\r\n" +
                " This could be caused if you do not have an active internet connection, administrative" +
                " right to access to internet, or host server may be down. Sorry.", false, "", "");
        }
    }

    [System.Runtime.Versioning.SupportedOSPlatform("windows")]
    public static bool Execute(string path)
    {
        try
        {
            if (File.Exists(path))
            {
                Process.Start(new ProcessStartInfo(path) { UseShellExecute = true });
                return true;
            }
            return false;
        }
        catch
        {
            return false;
        }

    }

    //260317Cl DownloadProgressChangedEventArgs→引数に変更 (WebClient依存を除去)
    //public static ... ProgressMessage(DownloadProgressChangedEventArgs e, Stopwatch stopwatch)
    public static (long Current, long Total, long ElapsedMilliseconds, string Message)
        ProgressMessage(long bytesReceived, long totalBytesToReceive, Stopwatch stopwatch)
    {
        var receivedMb = bytesReceived / 1E6;
        var totalMb = totalBytesToReceive / 1E6;
        var message = $"Downloading setup file.  Received: {receivedMb:f1} MB / {totalMb:f1} MB.  ";
        return (bytesReceived, totalBytesToReceive, stopwatch.ElapsedMilliseconds, message);
    }

    /// <summary>260317Cl 追加 HttpClientでファイルをダウンロードし進捗を報告する</summary>
    public static async System.Threading.Tasks.Task DownloadFileWithProgressAsync(
        string url, string path, IProgress<(long Current, long Total, long ElapsedMilliseconds, string Message)> progress, Stopwatch stopwatch)
    {
        using var response = await httpClient.GetAsync(new Uri(url), HttpCompletionOption.ResponseHeadersRead);
        response.EnsureSuccessStatusCode();
        var totalBytes = response.Content.Headers.ContentLength ?? -1;
        using var contentStream = await response.Content.ReadAsStreamAsync();
        using var fileStream = new FileStream(path, FileMode.Create, FileAccess.Write, FileShare.None, 8192, true);
        var buffer = new byte[8192];
        long bytesRead = 0;
        int read;
        long counter = 0;
        while ((read = await contentStream.ReadAsync(buffer)) > 0)
        {
            await fileStream.WriteAsync(buffer.AsMemory(0, read));
            bytesRead += read;
            if (counter++ % 10 == 0)
                progress?.Report(ProgressMessage(bytesRead, totalBytes, stopwatch));
        }
    }

}
