using System;
using System.Buffers;
using System.Collections.Generic;
using System.Runtime.CompilerServices;

//
//https://github.com/nenoNaninu/NthItemUtils/blob/master/NthItemUtils/NthItemUtils.cs
//

namespace Crystallography;


[Obsolete("261111Cl: 4 アプリで未使用。将来削除予定")] // 261111Cl 追加
public readonly struct ItemWithIndex<T>(T item, int index)
{
    public readonly T Item = item;
    public readonly int Index = index;
}

public static class NthItemExtensions
{
    /// <param name="n">0 ~ source.Count - 1</param>
    [Obsolete("261111Cl: 4 アプリで未使用。将来削除予定")] // 261111Cl 追加
    public static ItemWithIndex<T> NthLargest<T>(this IReadOnlyList<T> source, int n) where T : IComparable<T>
    {
        return source.NthSmallest(source.Count - 1 - n);
    }

    /// <param name="n">0 ~ source.Count - 1</param>
    [Obsolete("261111Cl: 4 アプリで未使用。将来削除予定")] // 261111Cl 追加
    public static ItemWithIndex<T> NthSmallest<T>(this IReadOnlyList<T> source, int n) where T : IComparable<T>
    {
        var pool = ArrayPool<int>.Shared.Rent(source.Count);

        try
        {
            var indices = pool.AsSpan(0, source.Count);
            QuickSelect.Iota(indices);
            QuickSelect.Execute(source, indices, n);

            return new ItemWithIndex<T>(source[indices[n]], indices[n]);
        }
        finally
        {
            ArrayPool<int>.Shared.Return(pool);
        }
    }

    /// <param name="n">0 ~ source.Count - 1</param>
    [Obsolete("261111Cl: 4 アプリで未使用。将来削除予定")] // 261111Cl 追加
    public static ItemWithIndex<T> NthLargest<T>(this Span<T> source, int n) where T : IComparable<T>
    {
        return source.NthSmallest(source.Length - 1 - n);
    }

    /// <param name="n">0 ~ source.Length - 1</param>
    [Obsolete("261111Cl: 4 アプリで未使用。将来削除予定")] // 261111Cl 追加
    public static ItemWithIndex<T> NthSmallest<T>(this Span<T> source, int n) where T : IComparable<T>
    {
        var pool = ArrayPool<int>.Shared.Rent(source.Length);

        try
        {
            var indices = pool.AsSpan(0, source.Length);
            QuickSelect.Iota(indices);
            QuickSelect.Execute<T>(source, indices, n);

            return new ItemWithIndex<T>(source[indices[n]], indices[n]);
        }
        finally
        {
            ArrayPool<int>.Shared.Return(pool);
        }
    }

    /// <param name="n">0 ~ source.Count - 1</param>
    [Obsolete("261111Cl: 4 アプリで未使用。将来削除予定")] // 261111Cl 追加
    public static ItemWithIndex<T> NthLargest<T>(this ReadOnlySpan<T> source, int n) where T : IComparable<T>
    {
        return source.NthSmallest(source.Length - 1 - n);
    }

    /// <param name="n">0 ~ source.Length - 1</param>
    [Obsolete("261111Cl: 4 アプリで未使用。将来削除予定")] // 261111Cl 追加
    public static ItemWithIndex<T> NthSmallest<T>(this ReadOnlySpan<T> source, int n) where T : IComparable<T>
    {
        var pool = ArrayPool<int>.Shared.Rent(source.Length);

        try
        {
            var indices = pool.AsSpan(0, source.Length);
            QuickSelect.Iota(indices);
            QuickSelect.Execute(source, indices, n);

            return new ItemWithIndex<T>(source[indices[n]], indices[n]);
        }
        finally
        {
            ArrayPool<int>.Shared.Return(pool);
        }
    }

    /// <param name="n">0 ~ source.Count - 1</param>
    [Obsolete("261111Cl: 4 アプリで未使用。将来削除予定")] // 261111Cl 追加
    public static ItemWithIndex<T> NthLargest<T>(this IReadOnlyList<T> source, int n, IComparer<T> comparer)
    {
        return source.NthSmallest(source.Count - 1 - n, comparer);
    }

    /// <param name="n">0 ~ source.Count - 1</param>
    [Obsolete("261111Cl: 4 アプリで未使用。将来削除予定")] // 261111Cl 追加
    public static ItemWithIndex<T> NthSmallest<T>(this IReadOnlyList<T> source, int n, IComparer<T> comparer)
    {
        var pool = ArrayPool<int>.Shared.Rent(source.Count);

        try
        {
            var indices = pool.AsSpan(0, source.Count);
            QuickSelect.Iota(indices);
            QuickSelect.Execute(source, indices, n, comparer);

            return new ItemWithIndex<T>(source[indices[n]], indices[n]);
        }
        finally
        {
            ArrayPool<int>.Shared.Return(pool);
        }
    }

    /// <param name="n">0 ~ source.Count - 1</param>
    [Obsolete("261111Cl: 4 アプリで未使用。将来削除予定")] // 261111Cl 追加
    public static ItemWithIndex<T> NthLargest<T>(this Span<T> source, int n, IComparer<T> comparer)
    {
        return source.NthSmallest(source.Length - 1 - n, comparer);
    }

    /// <param name="n">0 ~ source.Length - 1</param>
    [Obsolete("261111Cl: 4 アプリで未使用。将来削除予定")] // 261111Cl 追加
    public static ItemWithIndex<T> NthSmallest<T>(this Span<T> source, int n, IComparer<T> comparer)
    {
        var pool = ArrayPool<int>.Shared.Rent(source.Length);

        try
        {
            var indices = pool.AsSpan(0, source.Length);
            QuickSelect.Iota(indices);
            QuickSelect.Execute<T>(source, indices, n, comparer);

            return new ItemWithIndex<T>(source[indices[n]], indices[n]);
        }
        finally
        {
            ArrayPool<int>.Shared.Return(pool);
        }
    }

    /// <param name="n">0 ~ source.Count - 1</param>
    [Obsolete("261111Cl: 4 アプリで未使用。将来削除予定")] // 261111Cl 追加
    public static ItemWithIndex<T> NthLargest<T>(this ReadOnlySpan<T> source, int n, IComparer<T> comparer)
    {
        return source.NthSmallest(source.Length - 1 - n, comparer);
    }

    /// <param name="n">0 ~ source.Length - 1</param>
    [Obsolete("261111Cl: 4 アプリで未使用。将来削除予定")] // 261111Cl 追加
    public static ItemWithIndex<T> NthSmallest<T>(this ReadOnlySpan<T> source, int n, IComparer<T> comparer)
    {
        var pool = ArrayPool<int>.Shared.Rent(source.Length);

        try
        {
            var indices = pool.AsSpan(0, source.Length);
            QuickSelect.Iota(indices);
            QuickSelect.Execute(source, indices, n, comparer);

            return new ItemWithIndex<T>(source[indices[n]], indices[n]);
        }
        finally
        {
            ArrayPool<int>.Shared.Return(pool);
        }
    }

    #region Max

    [Obsolete("261111Cl: 4 アプリで未使用。将来削除予定")] // 261111Cl 追加
    public static ItemWithIndex<T> MaxWithIndex<T>(this IReadOnlyList<T> source) where T : IComparable<T>
    {
        T maxValue = source[0];
        int index = 0;

        for (int i = 1; i < source.Count; i++)
        {
            if (Comparer<T>.Default.Compare(source[i], maxValue) > 0)
            {
                maxValue = source[i];
                index = i;
            }
        }

        return new ItemWithIndex<T>(maxValue, index);
    }

    [Obsolete("261111Cl: 4 アプリで未使用。将来削除予定")] // 261111Cl 追加
    public static ItemWithIndex<T> MaxWithIndex<T>(this Span<T> source) where T : IComparable<T>
    {
        T maxValue = source[0];
        int index = 0;

        for (int i = 1; i < source.Length; i++)
        {
            if (Comparer<T>.Default.Compare(source[i], maxValue) > 0)
            {
                maxValue = source[i];
                index = i;
            }
        }

        return new ItemWithIndex<T>(maxValue, index);
    }

    [Obsolete("261111Cl: 4 アプリで未使用。将来削除予定")] // 261111Cl 追加
    public static ItemWithIndex<T> MaxWithIndex<T>(this ReadOnlySpan<T> source) where T : IComparable<T>
    {
        T maxValue = source[0];
        int index = 0;

        for (int i = 1; i < source.Length; i++)
        {
            if (Comparer<T>.Default.Compare(source[i], maxValue) > 0)
            {
                maxValue = source[i];
                index = i;
            }
        }

        return new ItemWithIndex<T>(maxValue, index);
    }

    [Obsolete("261111Cl: 4 アプリで未使用。将来削除予定")] // 261111Cl 追加
    public static ItemWithIndex<T> MaxWithIndex<T>(this IReadOnlyList<T> source, IComparer<T> comparer)
    {
        T maxValue = source[0];
        int index = 0;

        for (int i = 1; i < source.Count; i++)
        {
            if (comparer.Compare(source[i], maxValue) > 0)
            {
                maxValue = source[i];
                index = i;
            }
        }

        return new ItemWithIndex<T>(maxValue, index);
    }

    [Obsolete("261111Cl: 4 アプリで未使用。将来削除予定")] // 261111Cl 追加
    public static ItemWithIndex<T> MaxWithIndex<T>(this Span<T> source, IComparer<T> comparer)
    {
        T maxValue = source[0];
        int index = 0;

        for (int i = 1; i < source.Length; i++)
        {
            if (comparer.Compare(source[i], maxValue) > 0)
            {
                maxValue = source[i];
                index = i;
            }
        }

        return new ItemWithIndex<T>(maxValue, index);
    }

    [Obsolete("261111Cl: 4 アプリで未使用。将来削除予定")] // 261111Cl 追加
    public static ItemWithIndex<T> MaxWithIndex<T>(this ReadOnlySpan<T> source, IComparer<T> comparer)
    {
        T maxValue = source[0];
        int index = 0;

        for (int i = 1; i < source.Length; i++)
        {
            if (comparer.Compare(source[i], maxValue) > 0)
            {
                maxValue = source[i];
                index = i;
            }
        }

        return new ItemWithIndex<T>(maxValue, index);
    }

    #endregion //Max

    #region Min
    [Obsolete("261111Cl: 4 アプリで未使用。将来削除予定")] // 261111Cl 追加
    public static ItemWithIndex<T> MinWithIndex<T>(this IReadOnlyList<T> source) where T : IComparable<T>
    {
        T minValue = source[0];
        int index = 0;

        for (int i = 1; i < source.Count; i++)
        {
            if (Comparer<T>.Default.Compare(source[i], minValue) < 0)
            {
                minValue = source[i];
                index = i;
            }
        }

        return new ItemWithIndex<T>(minValue, index);
    }

    [Obsolete("261111Cl: 4 アプリで未使用。将来削除予定")] // 261111Cl 追加
    public static ItemWithIndex<T> MinWithIndex<T>(this Span<T> source) where T : IComparable<T>
    {
        T minValue = source[0];
        int index = 0;

        for (int i = 1; i < source.Length; i++)
        {
            if (Comparer<T>.Default.Compare(source[i], minValue) < 0)
            {
                minValue = source[i];
                index = i;
            }
        }

        return new ItemWithIndex<T>(minValue, index);
    }

    [Obsolete("261111Cl: 4 アプリで未使用。将来削除予定")] // 261111Cl 追加
    public static ItemWithIndex<T> MinWithIndex<T>(this ReadOnlySpan<T> source) where T : IComparable<T>
    {
        T minValue = source[0];
        int index = 0;

        for (int i = 1; i < source.Length; i++)
        {
            if (Comparer<T>.Default.Compare(source[i], minValue) < 0)
            {
                minValue = source[i];
                index = i;
            }
        }

        return new ItemWithIndex<T>(minValue, index);
    }

    [Obsolete("261111Cl: 4 アプリで未使用。将来削除予定")] // 261111Cl 追加
    public static ItemWithIndex<T> MinWithIndex<T>(this IReadOnlyList<T> source, IComparer<T> comparer)
    {
        T minValue = source[0];
        int index = 0;

        for (int i = 1; i < source.Count; i++)
        {
            if (comparer.Compare(source[i], minValue) < 0)
            {
                minValue = source[i];
                index = i;
            }
        }

        return new ItemWithIndex<T>(minValue, index);
    }

    [Obsolete("261111Cl: 4 アプリで未使用。将来削除予定")] // 261111Cl 追加
    public static ItemWithIndex<T> MinWithIndex<T>(this Span<T> source, IComparer<T> comparer)
    {
        T minValue = source[0];
        int index = 0;

        for (int i = 1; i < source.Length; i++)
        {
            if (comparer.Compare(source[i], minValue) < 0)
            {
                minValue = source[i];
                index = i;
            }
        }

        return new ItemWithIndex<T>(minValue, index);
    }

    [Obsolete("261111Cl: 4 アプリで未使用。将来削除予定")] // 261111Cl 追加
    public static ItemWithIndex<T> MinWithIndex<T>(this ReadOnlySpan<T> source, IComparer<T> comparer)
    {
        T minValue = source[0];
        int index = 0;

        for (int i = 1; i < source.Length; i++)
        {
            if (comparer.Compare(source[i], minValue) < 0)
            {
                minValue = source[i];
                index = i;
            }
        }

        return new ItemWithIndex<T>(minValue, index);
    }

    #endregion // Min
}

public static class QuickSelect
{
    public static void Iota(Span<int> indices)
    {
        for (int i = 0; i < indices.Length; i++)
            indices[i] = i;
    }

    #region Swap
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void Swap(Span<int> indices, int i, int j)
    {
        if (i == j) return;

        (indices[j], indices[i]) = (indices[i], indices[j]);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void Swap<T>(Span<T> indices, int i, int j)
    {
        if (i == j) return;

        (indices[j], indices[i]) = (indices[i], indices[j]);
    }
    #endregion

    #region 261111Cl 変更前: 末尾固定ピボットの分割 (同値が多い・整列済みの入力で O(n²))
    //#region Partition
    //private static int Partition<T>(ReadOnlySpan<T> source, Span<int> indices, int begin, int count, IComparer<T> comparer)
    //{
    //    if (count <= 1)
    //        return begin;
//
    //    int l = begin;
    //    int r = begin + count - 2;
    //    int endIndex = begin + count - 1;
//
    //    T pivot = source[indices[endIndex]];
//
    //    while (l <= r)
    //    {
    //        while (l < endIndex && comparer.Compare(source[indices[l]], pivot) < 0)
    //            l++;
//
    //        while (begin <= r && comparer.Compare(source[indices[r]], pivot) >= 0)
    //            r--;
//
    //        if (r < l)
    //            break;
//
    //        Swap(indices, l, r);
//
    //        l++;
    //        r--;
    //    }
//
    //    Swap(indices, l, endIndex);
//
    //    return l;
    //}
//
    //private static int Partition<T>(ReadOnlySpan<T> source, Span<int> indices, int begin, int count, Comparison<T> comparison)
    //{
    //    if (count <= 1)
    //        return begin;
//
    //    int l = begin;
    //    int r = begin + count - 2;
    //    int endIndex = begin + count - 1;
//
    //    T pivot = source[indices[endIndex]];
//
    //    while (l <= r)
    //    {
    //        while (l < endIndex && comparison(source[indices[l]], pivot) < 0)
    //            l++;
//
    //        while (begin <= r && comparison(source[indices[r]], pivot) >= 0)
    //            r--;
//
    //        if (r < l)
    //            break;
//
    //        Swap(indices, l, r);
//
    //        l++;
    //        r--;
    //    }
//
    //    Swap(indices, l, endIndex);
    //    return l;
    //}
//
    //private static int Partition<T>(Span<T> source, int begin, int count, Comparison<T> comparison)
    //{
    //    if (count <= 1)
    //        return begin;
//
    //    int l = begin, r = begin + count - 2, endIndex = begin + count - 1;
//
    //    T pivot = source[endIndex];
//
    //    while (l <= r)
    //    {
    //        while (l < endIndex && comparison(source[l], pivot) < 0)
    //            l++;
//
    //        while (begin <= r && comparison(source[r], pivot) >= 0)
    //            r--;
//
    //        if (r < l)
    //            break;
//
    //        (source[l], source[r]) = (source[r], source[l]);
//
    //        l++;
    //        r--;
    //    }
//
    //    (source[l], source[endIndex]) = (source[endIndex], source[l]);
    //    
    //    return l;
    //}
//
//
    //private static int Partition<T>(IReadOnlyList<T> source, Span<int> indices, int begin, int count, IComparer<T> comparer)
    //{
    //    if (count <= 1)
    //        return begin;
//
    //    int l = begin;
    //    int r = begin + count - 2;
    //    int endIndex = begin + count - 1;
//
    //    T pivot = source[indices[endIndex]];
//
    //    while (l <= r)
    //    {
    //        while (l < endIndex && comparer.Compare(source[indices[l]], pivot) < 0)
    //            l++;
//
    //        while (begin <= r && comparer.Compare(source[indices[r]], pivot) >= 0)
    //            r--;
//
    //        if (r < l)
    //            break;
//
    //        Swap(indices, l, r);
//
    //        l++;
    //        r--;
    //    }
//
    //    Swap(indices, l, endIndex);
//
    //    return l;
    //}
//
    //private static int Partition<T>(IReadOnlyList<T> source, Span<int> indices, int begin, int count, Comparison<T> comparison)
    //{
    //    if (count <= 1)
    //        return begin;
//
    //    int l = begin;
    //    int r = begin + count - 2;
    //    int endIndex = begin + count - 1;
//
    //    T pivot = source[indices[endIndex]];
//
    //    while (l <= r)
    //    {
    //        while (l < endIndex && comparison(source[indices[l]], pivot) < 0)
    //            l++;
//
    //        while (begin <= r && comparison(source[indices[r]], pivot) >= 0)
    //            r--;
//
    //        if (r < l)
    //            break;
//
    //        Swap(indices, l, r);
//
    //        l++;
    //        r--;
    //    }
//
    //    Swap(indices, l, endIndex);
//
    //    return l;
    //}
    //#endregion
    #endregion

    #region Validation
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void Validation<T>(IReadOnlyList<T> source, Span<int> indices, int n)
    {
        if (source.Count <= n)
            throw new ArgumentException("n is bigger than source.Count");

        if (indices.Length <= n)
            throw new ArgumentException("n is bigger than indices.Length");
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void Validation<T>(ReadOnlySpan<T> source, Span<int> indices, int n)
    {
        if (source.Length <= n)
            throw new ArgumentException("n is bigger than source.Length");

        if (indices.Length <= n)
            throw new ArgumentException("n is bigger than indices.Length");
    }
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void Validation<T>(Span<T> source, Span<int> indices, int n)
    {
        if (source.Length <= n)
            throw new ArgumentException("n is bigger than source.Length");

        if (indices.Length <= n)
            throw new ArgumentException("n is bigger than indices.Length");
    }

    private static void Validation<T>(Span<T> source, int n)
    {
        if (source.Length <= n)
            throw new ArgumentException("n is bigger than source.Length");
    }
    #endregion

    #region Execute
    /// <param name="n">0 ~ source.Count - 1</param>
    [Obsolete("261111Cl: 4 アプリで未使用。将来削除予定")] // 261111Cl 追加
    public static void Execute<T>(IReadOnlyList<T> source, Span<int> indices, int n) where T : IComparable<T>
    {
        Validation(source, indices, n);
        //ExecuteCore(source, indices, n, 0, indices.Length, Comparer<T>.Default); // 261111Cl 変更前 (末尾固定ピボット)
        SelectCore(indices, n, new ListComparerOrder<T>(source, Comparer<T>.Default)); // 261111Cl 変更: 全順序 (値 → 添字) で選ぶ
    }

    /// <param name="n">0 ~ source.Count - 1</param>
    [Obsolete("261111Cl: 4 アプリで未使用。将来削除予定")] // 261111Cl 追加
    public static void Execute<T>(IReadOnlyList<T> source, Span<int> indices, int n, IComparer<T> comparer)
    {
        Validation(source, indices, n);
        //ExecuteCore(source, indices, n, 0, indices.Length, comparer); // 261111Cl 変更前 (末尾固定ピボット)
        SelectCore(indices, n, new ListComparerOrder<T>(source, comparer)); // 261111Cl 変更: 全順序 (値 → 添字) で選ぶ
    }

    /// <param name="n">0 ~ source.Count - 1</param>
    [Obsolete("261111Cl: 4 アプリで未使用。将来削除予定")] // 261111Cl 追加
    public static void Execute<T>(IReadOnlyList<T> source, Span<int> indices, int n, Comparison<T> comparison)
    {
        Validation(source, indices, n);
        //ExecuteCore(source, indices, n, 0, indices.Length, comparison); // 261111Cl 変更前 (末尾固定ピボット)
        SelectCore(indices, n, new ListComparisonOrder<T>(source, comparison)); // 261111Cl 変更: 全順序 (値 → 添字) で選ぶ
    }

    /// <param name="n">0 ~ source.Length - 1</param>
    [Obsolete("261111Cl: 4 アプリで未使用。将来削除予定")] // 261111Cl 追加
    public static void Execute<T>(ReadOnlySpan<T> source, Span<int> indices, int n) where T : IComparable<T>
    {
        Validation(source, indices, n);
        //ExecuteCore(source, indices, n, 0, indices.Length, Comparer<T>.Default); // 261111Cl 変更前 (末尾固定ピボット)
        SelectCore(indices, n, new SpanComparerOrder<T>(source, Comparer<T>.Default)); // 261111Cl 変更: 全順序 (値 → 添字) で選ぶ
    }


    /// <param name="n">0 ~ source.Length - 1</param>
    [Obsolete("261111Cl: 4 アプリで未使用。将来削除予定")] // 261111Cl 追加
    public static void Execute<T>(ReadOnlySpan<T> source, Span<int> indices, int n, IComparer<T> comparer)
    {
        Validation(source, indices, n);
        //ExecuteCore(source, indices, n, 0, indices.Length, comparer); // 261111Cl 変更前 (末尾固定ピボット)
        SelectCore(indices, n, new SpanComparerOrder<T>(source, comparer)); // 261111Cl 変更: 全順序 (値 → 添字) で選ぶ
    }

    /// <param name="n">0 ~ source.Length - 1</param>
    [Obsolete("261111Cl: 4 アプリで未使用。将来削除予定")] // 261111Cl 追加
    public static void Execute<T>(ReadOnlySpan<T> source, Span<int> indices, int n, Comparison<T> comparison)
    {
        Validation(source, indices, n);
        //ExecuteCore(source, indices, n, 0, indices.Length, comparison); // 261111Cl 変更前 (末尾固定ピボット)
        SelectCore(indices, n, new SpanComparisonOrder<T>(source, comparison)); // 261111Cl 変更: 全順序 (値 → 添字) で選ぶ
    }

    /// <param name="n">0 ~ source.Length - 1</param>
    [Obsolete("261111Cl: 4 アプリで未使用。将来削除予定")] // 261111Cl 追加
    public static void Execute<T>(Span<T> source, Span<int> indices, int n, Comparison<T> comparison)
    {
        Validation(source, indices, n);
        //ExecuteCore(source, indices, n, 0, indices.Length, comparison); // 261111Cl 変更前 (末尾固定ピボット)
        SelectCore(indices, n, new SpanComparisonOrder<T>(source, comparison)); // 261111Cl 変更: 全順序 (値 → 添字) で選ぶ
    }

    /// <param name="n">0 ~ source.Length - 1</param>
    public static void Execute<T>(Span<T> source, int n, Comparison<T> comparison)
    {
        Validation(source, n);
        //ExecuteCore(source, n, 0, comparison); // 261111Cl 変更前 (末尾固定ピボット)
        SelectInPlace(source, n, comparison); // 261111Cl 変更: 添字で選んでから並べ直す
    }
    #endregion

    #region 261111Cl 変更前: 旧 ExecuteCore (旧 Partition を繰り返す)
    //#region ExecuteCore
    //[MethodImpl(MethodImplOptions.AggressiveInlining)]
    //private static void ExecuteCore<T>(IReadOnlyList<T> source, Span<int> indices, int n, int begin, int count, IComparer<T> comparer)
    //{
    //    int partition = Partition(source, indices, begin, count, comparer);
//
    //    while (partition != n)
    //    {
    //        if (partition < n)
    //        {
    //            count = begin + count - partition - 1;
    //            begin = partition + 1;
    //        }
    //        else if (n < partition)
    //            count = partition - begin;
//
    //        partition = Partition(source, indices, begin, count, comparer);
    //    }
    //}
//
    //[MethodImpl(MethodImplOptions.AggressiveInlining)]
    //private static void ExecuteCore<T>(IReadOnlyList<T> source, Span<int> indices, int n, int begin, int count, Comparison<T> comparison)
    //{
    //    int partition = Partition(source, indices, begin, count, comparison);
//
    //    while (partition != n)
    //    {
    //        if (partition < n)
    //        {
    //            count = begin + count - partition - 1;
    //            begin = partition + 1;
    //        }
    //        else if (n < partition)
    //            count = partition - begin;
//
    //        partition = Partition(source, indices, begin, count, comparison);
    //    }
    //}
//
    //[MethodImpl(MethodImplOptions.AggressiveInlining)]
    //private static void ExecuteCore<T>(ReadOnlySpan<T> source, Span<int> indices, int n, int begin, int count, IComparer<T> comparer)
    //{
    //    int partition = Partition(source, indices, begin, count, comparer);
//
    //    while (partition != n)
    //    {
    //        if (partition < n)
    //        {
    //            count = begin + count - partition - 1;
    //            begin = partition + 1;
    //        }
    //        else if (n < partition)
    //            count = partition - begin;
//
    //        partition = Partition(source, indices, begin, count, comparer);
    //    }
    //}
//
    //[MethodImpl(MethodImplOptions.AggressiveInlining)]
    //private static void ExecuteCore<T>(ReadOnlySpan<T> source, Span<int> indices, int n, int begin, int count, Comparison<T> comparison)
    //{
    //    int partition = Partition(source, indices, begin, count, comparison);
//
    //    while (partition != n)
    //    {
    //        if (partition < n)
    //        {
    //            count = begin + count - partition - 1;
    //            begin = partition + 1;
    //        }
    //        else if (n < partition)
    //            count = partition - begin;
//
    //        partition = Partition(source, indices, begin, count, comparison);
    //    }
    //}
//
    //[MethodImpl(MethodImplOptions.AggressiveInlining)]
    //private static void ExecuteCore<T>(Span<T> source, Span<int> indices, int n, int begin, int count, Comparison<T> comparison)
    //{
    //    int partition = Partition(source, indices, begin, count, comparison);
//
    //    while (partition != n)
    //    {
    //        if (partition < n)
    //        {
    //            count = begin + count - partition - 1;
    //            begin = partition + 1;
    //        }
    //        else if (n < partition)
    //            count = partition - begin;
//
    //        partition = Partition(source, indices, begin, count, comparison);
    //    }
    //}
//
    //[MethodImpl(MethodImplOptions.AggressiveInlining)]
    //private static void ExecuteCore<T>(Span<T> source, int n, int begin, Comparison<T> comparison)
    //{
    //    int count = source.Length;
    //    int partition = Partition(source, begin, count, comparison);
//
    //    while (partition != n)
    //    {
    //        if (partition < n)
    //        {
    //            count = begin + count - partition - 1;
    //            begin = partition + 1;
    //        }
    //        else if (n < partition)
    //            count = partition - begin;
//
    //        partition = Partition(source, begin, count, comparison);
    //    }
    //}
    //#endregion
    #endregion

    #region 261111Cl 追加: 全順序による選択 (同値が多い入力・整列済みの入力でも速い)
    // 261111Cl 追加: 旧版は末尾固定のピボットで分割していたため、同値が多い入力 (例: 0 が過半の画像の |応答|) や
    //   整列済み (に近い) 入力で O(n²) になった (18 万要素・0 が 55% で 110 秒超)。ここでは比較を
    //   「値 → 同じ値なら添字 (元の位置) の小さい順」の全順序にして選ぶ。
    //   - 全順序なので同値の塊ができず、3 点の中央値のピボットで期待 O(n)。分割の回数が上限を超えたら、残りの区間をヒープソートする (最悪 O(n log n))。
    //   - 約束は旧版と同じ: indices[n] (その場版は source[n]) が小さい方から n 番目、前はそれ以下、後ろはそれ以上。
    //   - 同値が n の境界をまたぐとき、前に入るのは元の位置の小さいもの (= 安定な並べ替えの先頭 n 個と同じ集合)。
    //     旧版はここが分割の手順次第だった (ビーム選択で境界の同値の組がどれになるかが、旧版とは変わり得る)。

    private interface IIndexOrder
    {
        int Compare(int a, int b);
    }

    private readonly ref struct SpanComparisonOrder<T> : IIndexOrder
    {
        private readonly ReadOnlySpan<T> source;
        private readonly Comparison<T> comparison;
        public SpanComparisonOrder(ReadOnlySpan<T> source, Comparison<T> comparison) { this.source = source; this.comparison = comparison; }
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public int Compare(int a, int b) { int c = comparison(source[a], source[b]); return c != 0 ? c : a.CompareTo(b); }
    }

    private readonly ref struct SpanComparerOrder<T> : IIndexOrder
    {
        private readonly ReadOnlySpan<T> source;
        private readonly IComparer<T> comparer;
        public SpanComparerOrder(ReadOnlySpan<T> source, IComparer<T> comparer) { this.source = source; this.comparer = comparer; }
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public int Compare(int a, int b) { int c = comparer.Compare(source[a], source[b]); return c != 0 ? c : a.CompareTo(b); }
    }

    private readonly struct ListComparisonOrder<T>(IReadOnlyList<T> source, Comparison<T> comparison) : IIndexOrder
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public int Compare(int a, int b) { int c = comparison(source[a], source[b]); return c != 0 ? c : a.CompareTo(b); }
    }

    private readonly struct ListComparerOrder<T>(IReadOnlyList<T> source, IComparer<T> comparer) : IIndexOrder
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public int Compare(int a, int b) { int c = comparer.Compare(source[a], source[b]); return c != 0 ? c : a.CompareTo(b); }
    }

    /// <summary>261111Cl 追加: idx の中で order による n 番目を idx[n] に置く (前は小さく、後ろは大きい)。</summary>
    private static void SelectCore<TOrder>(Span<int> idx, int n, TOrder order) where TOrder : IIndexOrder, allows ref struct
    {
        int lo = 0, hi = idx.Length - 1;
        int budget = 2 * (System.Numerics.BitOperations.Log2((uint)Math.Max(idx.Length, 1)) + 2);
        while (hi - lo >= 16)
        {
            if (budget-- <= 0)
            {
                HeapSort(idx.Slice(lo, hi - lo + 1), order);
                return;
            }
            // 3 点の中央値をピボットにし、idx[lo] <= pivot <= idx[hi] を番兵にする
            int mid = lo + ((hi - lo) >> 1);
            if (order.Compare(idx[mid], idx[lo]) < 0) Swap(idx, mid, lo);
            if (order.Compare(idx[hi], idx[lo]) < 0) Swap(idx, hi, lo);
            if (order.Compare(idx[hi], idx[mid]) < 0) Swap(idx, hi, mid);
            Swap(idx, mid, hi - 1);
            int pivot = idx[hi - 1];
            int i = lo, j = hi - 1;
            while (true)
            {
                while (order.Compare(idx[++i], pivot) < 0) { }
                while (order.Compare(pivot, idx[--j]) < 0) { }
                if (i >= j) break;
                Swap(idx, i, j);
            }
            Swap(idx, i, hi - 1);
            if (i == n) return;
            if (n < i) hi = i - 1; else lo = i + 1;
        }
        // 小さい区間は挿入ソート
        var s = idx.Slice(lo, hi - lo + 1);
        for (int i = 1; i < s.Length; i++)
        {
            int v = s[i], j = i - 1;
            while (j >= 0 && order.Compare(v, s[j]) < 0) { s[j + 1] = s[j]; j--; }
            s[j + 1] = v;
        }
    }

    /// <summary>261111Cl 追加: SelectCore の深さの上限を超えたときの退避 (最悪 O(n log n))</summary>
    private static void HeapSort<TOrder>(Span<int> s, TOrder order) where TOrder : IIndexOrder, allows ref struct
    {
        int n = s.Length;
        for (int i = n / 2 - 1; i >= 0; i--) SiftDown(s, i, n, order);
        for (int end = n - 1; end > 0; end--)
        {
            (s[0], s[end]) = (s[end], s[0]);
            SiftDown(s, 0, end, order);
        }
    }

    /// <summary>261111Cl 追加</summary>
    private static void SiftDown<TOrder>(Span<int> s, int i, int n, TOrder order) where TOrder : IIndexOrder, allows ref struct
    {
        int v = s[i];
        while (true)
        {
            int c = 2 * i + 1;
            if (c >= n) break;
            if (c + 1 < n && order.Compare(s[c], s[c + 1]) < 0) c++;
            if (order.Compare(v, s[c]) >= 0) break;
            s[i] = s[c];
            i = c;
        }
        s[i] = v;
    }

    /// <summary>261111Cl 追加: Span&lt;T&gt; をその場で並べる版。添字で選んでから、選んだ順に要素を並べ直す</summary>
    private static void SelectInPlace<T>(Span<T> source, int n, Comparison<T> comparison)
    {
        int len = source.Length;
        var idxPool = ArrayPool<int>.Shared.Rent(len);
        var tmpPool = ArrayPool<T>.Shared.Rent(len);
        try
        {
            var idx = idxPool.AsSpan(0, len);
            Iota(idx);
            SelectCore(idx, n, new SpanComparisonOrder<T>(source, comparison));
            var tmp = tmpPool.AsSpan(0, len);
            for (int i = 0; i < len; i++)
                tmp[i] = source[idx[i]];
            tmp.CopyTo(source);
        }
        finally
        {
            ArrayPool<int>.Shared.Return(idxPool);
            ArrayPool<T>.Shared.Return(tmpPool, RuntimeHelpers.IsReferenceOrContainsReferences<T>());
        }
    }
    #endregion

}
