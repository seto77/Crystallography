using System;
using System.Drawing;
using System.Drawing.Imaging;
using System.IO;
using System.Runtime.InteropServices;

namespace Crystallography
{
    public class ClipboardMetafileHelper
    {
        [DllImport("user32.dll")]
        private static extern bool OpenClipboard(IntPtr hWndNewOwner);

        [DllImport("user32.dll")]
        private static extern bool EmptyClipboard();

        [DllImport("user32.dll")]
        private static extern IntPtr SetClipboardData(uint uFormat, IntPtr hMem);

        [DllImport("user32.dll")]
        private static extern bool CloseClipboard();

        [DllImport("gdi32.dll")]
        private static extern IntPtr CopyEnhMetaFile(IntPtr hemfSrc, IntPtr hNULL);

        [DllImport("gdi32.dll")]
        private static extern bool DeleteEnhMetaFile(IntPtr hemf);

        // Metafile mf is set to an invalid state inside this function
        static public bool PutEnhMetafileOnClipboard(IntPtr hWnd, Metafile mf)
        {
            //bool bResult = false; // (260715Ch) 旧: 失敗経路で複製 EMF と OpenClipboard を後始末できなかった
            IntPtr hEMF = IntPtr.Zero, hEMF2 = IntPtr.Zero; // (260715Ch)
            bool clipboardOpened = false; // (260715Ch)
            try
            {
                hEMF = mf.GetHenhmetafile(); // invalidates mf
                if (hEMF == IntPtr.Zero)
                    return false;

                hEMF2 = CopyEnhMetaFile(hEMF, new IntPtr(0));
                if (hEMF2 == IntPtr.Zero || !OpenClipboard(hWnd))
                    return false;

                clipboardOpened = true; // (260715Ch)
                if (!EmptyClipboard())
                    return false;

                //IntPtr hRes = SetClipboardData(14 /*CF_ENHMETAFILE*/, hEMF2); // 旧: 失敗時も hEMF2 を解放せず、成功時との所有権を区別しなかった
                if (SetClipboardData(14 /*CF_ENHMETAFILE*/, hEMF2) != hEMF2)
                    return false;

                hEMF2 = IntPtr.Zero; // (260715Ch) 成功後の所有権は Windows clipboard に移る
                return true;
            }
            finally
            {
                if (clipboardOpened)
                    CloseClipboard(); // (260715Ch) EmptyClipboard / SetClipboardData の失敗時も必ず閉じる
                if (hEMF2 != IntPtr.Zero)
                    DeleteEnhMetaFile(hEMF2); // (260715Ch) 所有権移譲前の複製ハンドルだけを解放
                if (hEMF != IntPtr.Zero)
                    DeleteEnhMetaFile(hEMF); // (260715Ch)
            }
        }

        // 260504Cl 追加: 任意の描画アクションを EMF+ 化してクリップボードへ書き込む。
        // 既存 5 箇所 (ScalablePictureBox / FormStereonet / FormDiffractionSimulator / FormImageSimulator
        // / FormSymmetryInformation) で同じ HDC→Metafile→PutEnh… の手順が複製されていたので集約。
        // draw 引数では SmoothingMode / Clear など Graphics の初期状態を呼び出し側で設定する。
        public static bool PutDrawingOnClipboardAsEnhMetafile(IntPtr hWnd, Action<Graphics> draw)
            => SaveOrCopyDrawingAsEnhMetafile(hWnd, draw); // 260716Cl 旧本体は SaveOrCopyDrawingAsEnhMetafile へ移動 (EMF ファイル保存 sink を追加して一般化)

        // 260716Cl 追加: 描画アクションを EMF+ 録画し、filename が空ならクリップボードへ、指定があればファイルへ書き出す。
        // FormDiffractionSimulator / FormImageSimulator がファイル保存経路のために自前で複製していた
        // HDC→Metafile 生存管理 (コンストラクタ失敗時の ReleaseHdc を含む) をここへ集約。
        public static bool SaveOrCopyDrawingAsEnhMetafile(IntPtr hWnd, Action<Graphics> draw, string filename = "")
        {
            ArgumentNullException.ThrowIfNull(draw);
            using var refG = Graphics.FromHwnd(hWnd);
            IntPtr hdc = refG.GetHdc();
            using var ms = new MemoryStream();
            //using var mf = new Metafile(ms, hdc, EmfType.EmfPlusDual); // 旧: コンストラクタ例外時に HDC を解放できなかった
            Metafile mf;
            try
            {
                mf = new Metafile(ms, hdc, EmfType.EmfPlusDual); // (260715Ch)
            }
            finally
            {
                refG.ReleaseHdc(hdc); // (260715Ch)
            }
            using (mf) // (260715Ch)
            {
                using (var g = Graphics.FromImage(mf))
                    draw(g);
                if (string.IsNullOrEmpty(filename))
                    return PutEnhMetafileOnClipboard(hWnd, mf);
                using var fsm = new FileStream(filename, FileMode.Create, FileAccess.Write); // 260716Cl EMF 内容確定後にファイルへ書き出す
                fsm.Write(ms.GetBuffer(), 0, (int)ms.Length);
                return true;
            }
        }
    }
}
