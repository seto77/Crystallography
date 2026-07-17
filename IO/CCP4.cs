using System.IO;

namespace Crystallography;

public class CCP4
{
    private static readonly uint[] setbits =
                     [0x00000000, 0x00000001, 0x00000003, 0x00000007,
                          0x0000000F, 0x0000001F, 0x0000003F, 0x0000007F,
                          0x000000FF, 0x000001FF, 0x000003FF, 0x000007FF,
                          0x00000FFF, 0x00001FFF, 0x00003FFF, 0x00007FFF,
                          0x0000FFFF, 0x0001FFFF, 0x0003FFFF, 0x0007FFFF,
                          0x000FFFFF, 0x001FFFFF, 0x003FFFFF, 0x007FFFFF,
                          0x00FFFFFF, 0x01FFFFFF, 0x03FFFFFF, 0x07FFFFFF,
                          0x0FFFFFFF, 0x1FFFFFFF, 0x3FFFFFFF, 0x7FFFFFFF,
                          0xFFFFFFFF];

    private static uint shift_left(uint x, int n)
    { return (((x) & (uint)setbits[32 - (n)]) << (n)); }

    private static uint shift_right(uint x, int n)
    { return (((x) >> (n)) & (uint)setbits[32 - (n)]); }

    public static unsafe uint[] unpack(BinaryReader br, int x, int y)
    {
        int valids = 0, spillbits = 0, usedbits;
        int total = x * y;
        uint window = 0;
        uint spill = 0, pixel = 0, nextint;
        int bitnum;
        int pixnum;
        int[] bitdecode = [0, 4, 5, 6, 7, 8, 16, 32];
        uint[] img = new uint[total];

        while (pixel < total)
        {
            if (valids < 6)
            {
                if (spillbits > 0)
                {
                    window |= shift_left(spill, valids);
                    valids += spillbits;
                    spillbits = 0;
                }
                else
                {
                    spill = br.ReadByte();
                    spillbits = 8;
                }
            }
            else
            {
                pixnum = 1 << (int)(window & setbits[3]);
                window = shift_right(window, 3);
                bitnum = bitdecode[window & setbits[3]];
                window = shift_right(window, 3);
                valids -= 6;
                while ((pixnum > 0) && (pixel < total))
                {
                    if (valids < bitnum)
                    {
                        if (spillbits > 0)
                        {
                            window |= shift_left(spill, valids);
                            if ((32 - valids) > spillbits)
                            {
                                valids += spillbits;
                                spillbits = 0;
                            }
                            else
                            {
                                usedbits = 32 - valids;
                                spill = shift_right(spill, usedbits);
                                spillbits -= usedbits;
                                valids = 32;
                            }
                        }
                        else
                        {
                            spill = br.ReadByte();
                            spillbits = 8;
                        }
                    }
                    else
                    {
                        --pixnum;
                        if (bitnum == 0)
                            nextint = 0;
                        else
                        {
                            nextint = window & setbits[bitnum];
                            valids -= bitnum;
                            window = shift_right(window, bitnum);
                            if ((nextint & (1 << (bitnum - 1))) != 0)
                                nextint |= (~(setbits[bitnum]));
                        }
                        if (pixel > x)
                        {
                            img[pixel] = (nextint +
                                              (img[pixel - 1] + img[pixel - x + 1] +
                                               img[pixel - x] + img[pixel - x - 1] + 2) / 4);
                            ++pixel;
                        }
                        else if (pixel != 0)
                        {
                            img[pixel] = (img[pixel - 1] + nextint);
                            ++pixel;
                        }
                        else
                            img[pixel++] = nextint;
                    }
                }
            }
        }
        return img;
    }

    public static unsafe uint[] unpack_long(BinaryReader br, int x, int y)
    // void unpack_long(FILE *packfile, int x, int y, LONG *img)
    {
        int valids = 0, spillbits = 0, usedbits;
        int total = x * y;
        uint window = 0;
        uint spill = 0, pixel = 0, nextint;
        int bitnum;
        int pixnum;
        int[] bitdecode = [0, 4, 5, 6, 7, 8, 16, 32];
        uint[] img = new uint[total];

        //int valids = 0, spillbits = 0, usedbits, total = x * y;
        //long window = 0L, spill = 0, pixel = 0, nextint, bitnum, pixnum;
        //int bitdecode = new []{0, 4, 5, 6, 7, 8, 16, 32};

        while (pixel < total)
        {
            if (valids < 6)
            {
                if (spillbits > 0)
                {
                    window |= shift_left(spill, valids);
                    valids += spillbits;
                    spillbits = 0;
                }
                else
                {
                    spill = br.ReadByte();
                    spillbits = 8;
                }
            }
            else
            {
                pixnum = 1 << (int)(window & setbits[3]);
                window = shift_right(window, 3);
                bitnum = bitdecode[window & setbits[3]];
                window = shift_right(window, 3);
                valids -= 6;
                while ((pixnum > 0) && (pixel < total))
                {
                    if (valids < bitnum)
                    {
                        if (spillbits > 0)
                        {
                            window |= shift_left(spill, valids);
                            if ((32 - valids) > spillbits)
                            {
                                valids += spillbits;
                                spillbits = 0;
                            }
                            else
                            {
                                usedbits = 32 - valids;
                                spill = shift_right(spill, usedbits);
                                spillbits -= usedbits;
                                valids = 32;
                            }
                        }
                        else
                        {
                            spill = br.ReadByte();
                            spillbits = 8;
                        }
                    }
                    else
                    {
                        --pixnum;
                        if (bitnum == 0)
                            nextint = 0;
                        else
                        {
                            nextint = window & setbits[bitnum];
                            valids -= bitnum;
                            window = shift_right(window, bitnum);
                            if ((nextint & (1 << (bitnum - 1))) != 0)
                                nextint |= ~setbits[bitnum];
                        }
                        if (pixel > x)
                        {
                            img[pixel] = (nextint +
                                            (img[pixel - 1] + img[pixel - x + 1] +
                                           img[pixel - x] + img[pixel - x - 1] + 2) / 4);
                            ++pixel;
                        }
                        else if (pixel != 0)
                        {
                            img[pixel] = (img[pixel - 1] + nextint);
                            ++pixel;
                        }
                        else
                            img[pixel++] = nextint;
                    }
                }
            }
        }
        return img;
    }

    public static unsafe uint[] v2unpack(BinaryReader br, int x, int y)
    {
        int valids = 0, spillbits = 0, usedbits;
        int total = x * y;
        uint window = 0;
        uint spill = 0, pixel = 0, nextint;
        int bitnum;
        int pixnum;
        int[] bitdecode = [0, 4, 5, 6, 7, 8, 16, 32];
        uint[] img = new uint[total];

        while (pixel < total)
        {
            if (valids < 7)
            {
                if (spillbits > 0)
                {
                    window |= shift_left(spill, valids);
                    valids += spillbits;
                    spillbits = 0;
                }
                else
                {
                    spill = br.ReadByte();
                    spillbits = 8;
                }
            }
            else
            {
                pixnum = 1 << (int)(window & (int)setbits[3]);
                window = shift_right(window, 3);
                bitnum = bitdecode[window & setbits[4]];
                window = shift_right(window, 4);
                valids -= 7;
                while ((pixnum > 0) && (pixel < total))
                {
                    if (valids < bitnum)
                    {
                        if (spillbits > 0)
                        {
                            window |= shift_left(spill, valids);
                            if ((32 - valids) > spillbits)
                            {
                                valids += spillbits;
                                spillbits = 0;
                            }
                            else
                            {
                                usedbits = 32 - valids;
                                spill = shift_right(spill, usedbits);
                                spillbits -= usedbits;
                                valids = 32;
                            }
                        }
                        else
                        {
                            spill = br.ReadByte();
                            spillbits = 8;
                        }
                    }
                    else
                    {
                        --pixnum;
                        if (bitnum == 0)
                            nextint = 0;
                        else
                        {
                            nextint = window & setbits[bitnum];
                            valids -= bitnum;
                            window = shift_right(window, bitnum);
                            if ((nextint & (1 << (bitnum - 1))) != 0)
                                nextint |= ~setbits[bitnum];
                        }
                        if (pixel > x)
                        {
                            img[pixel] = (uint)(nextint +
                                          (img[pixel - 1] + img[pixel - x + 1] +
                                                           img[pixel - x] + img[pixel - x - 1] + 2) / 4);
                            ++pixel;
                        }
                        else if (pixel != 0)
                        {
                            img[pixel] = (uint)(img[pixel - 1] + nextint);
                            ++pixel;
                        }
                        else
                            img[pixel++] = (uint)nextint;
                    }
                }
            }
        }
        return img;
    }

    // unpack / unpack_long は marControl の pck.c (Dr. Claudio Klein, X-ray Research GmbH, v1.1, 30/10/1995) を C# へ移植した実装。原典 C ソースは冗長なため削除 (git 履歴で参照可)。
}
