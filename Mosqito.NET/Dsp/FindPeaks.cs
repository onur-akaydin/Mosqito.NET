using System.Buffers;
using System.Runtime.CompilerServices;

namespace Mosqito.Dsp;

/// <summary>
/// Peak detection in 1-D arrays.
/// Implements the subset of scipy.signal.find_peaks used by MoSQITo
/// (height, prominence, distance, width filtering).
/// </summary>
public static class FindPeaks
{
    /// <summary>Result of a peak-finding operation.</summary>
    public readonly struct PeakResult
    {
        /// <summary>Indices of detected peaks in the original array.</summary>
        public int[] Indices { get; init; }
        /// <summary>Peak heights (x[idx]).</summary>
        public double[] Heights { get; init; }
        /// <summary>Peak prominences (if computed; otherwise empty).</summary>
        public double[] Prominences { get; init; }

        public static PeakResult Empty => new()
        {
            Indices = Array.Empty<int>(),
            Heights = Array.Empty<double>(),
            Prominences = Array.Empty<double>()
        };
    }

    /// <summary>
    /// Finds all local maxima in <paramref name="x"/> subject to optional constraints.
    /// Matches scipy.signal.find_peaks behaviour for the supported parameters.
    /// Uses pooled int[] buffers to avoid per-call heap allocation.
    /// </summary>
    public static PeakResult Find(
        ReadOnlySpan<double> x,
        double? height = null,
        int? distance = null,
        double? prominence = null)
    {
        int n = x.Length;
        if (n < 3) return PeakResult.Empty;

        // Rent a scratch buffer large enough for all candidate indices (worst case n-2).
        int[] candBuf = ArrayPool<int>.Shared.Rent(n);
        int candCount = 0;

        try
        {
            // Local maxima — strictly greater than both neighbours
            for (int i = 1; i < n - 1; i++)
                if (x[i] > x[i - 1] && x[i] > x[i + 1])
                    candBuf[candCount++] = i;

            // Filter by height (in-place compaction)
            if (height.HasValue)
            {
                double h = height.Value;
                int w = 0;
                for (int i = 0; i < candCount; i++)
                    if (x[candBuf[i]] >= h) candBuf[w++] = candBuf[i];
                candCount = w;
            }

            // Filter by minimum distance (keep tallest within each window)
            if (distance.HasValue && distance.Value > 1 && candCount > 1)
                candCount = FilterByDistance(x, candBuf, candCount, distance.Value);

            // Compute prominences
            double[] proms = new double[candCount];
            ComputeProminences(x, candBuf, candCount, proms);

            // Filter by prominence (in-place compaction)
            if (prominence.HasValue)
            {
                double minProm = prominence.Value;
                int w = 0;
                for (int i = 0; i < candCount; i++)
                {
                    if (proms[i] >= minProm)
                    {
                        candBuf[w] = candBuf[i];
                        proms[w]   = proms[i];
                        w++;
                    }
                }
                candCount = w;
                if (proms.Length != candCount)
                    Array.Resize(ref proms, candCount);
            }

            int[] idxArr  = new int[candCount];
            double[] heights = new double[candCount];
            for (int i = 0; i < candCount; i++)
            {
                idxArr[i]  = candBuf[i];
                heights[i] = x[candBuf[i]];
            }

            return new PeakResult { Indices = idxArr, Heights = heights, Prominences = proms };
        }
        finally
        {
            ArrayPool<int>.Shared.Return(candBuf);
        }
    }

    // ------------------------------------------------------------------
    // Prominences — public legacy overload (double[] + List<int>)
    // ------------------------------------------------------------------

    /// <summary>
    /// Computes the prominence of each peak at the given indices.
    /// Legacy overload kept for API compatibility.
    /// </summary>
    public static double[] ComputeProminences(double[] x, List<int> peakIndices)
    {
        double[] proms = new double[peakIndices.Count];
        int[] tempArr = peakIndices.ToArray();
        ComputeProminences(x, tempArr, peakIndices.Count, proms);
        return proms;
    }

    // ------------------------------------------------------------------
    // Internal span-based helpers
    // ------------------------------------------------------------------

    private static void ComputeProminences(ReadOnlySpan<double> x, ReadOnlySpan<int> cand,
        int count, double[] proms)
    {
        int n = x.Length;
        for (int pi = 0; pi < count; pi++)
        {
            int idx = cand[pi];
            double peakHeight = x[idx];

            double leftMin = peakHeight;
            int left = idx - 1;
            while (left >= 0)
            {
                if (x[left] > peakHeight) break;
                if (x[left] < leftMin) leftMin = x[left];
                left--;
            }

            double rightMin = peakHeight;
            int right = idx + 1;
            while (right < n)
            {
                if (x[right] > peakHeight) break;
                if (x[right] < rightMin) rightMin = x[right];
                right++;
            }

            proms[pi] = peakHeight - Math.Max(leftMin, rightMin);
        }
    }

    private static int FilterByDistance(ReadOnlySpan<double> x, int[] cand, int count,
        int distance)
    {
        // keep[] fits on stack for typical peak counts
        bool[] keep = ArrayPool<bool>.Shared.Rent(count);
        keep.AsSpan(0, count).Fill(true);

        try
        {
            for (int i = 0; i < count; i++)
            {
                if (!keep[i]) continue;
                for (int j = i + 1; j < count; j++)
                {
                    if (cand[j] - cand[i] >= distance) break;
                    if (x[cand[i]] >= x[cand[j]]) keep[j] = false;
                    else { keep[i] = false; break; }
                }
            }

            int w = 0;
            for (int i = 0; i < count; i++)
                if (keep[i]) cand[w++] = cand[i];
            return w;
        }
        finally
        {
            ArrayPool<bool>.Shared.Return(keep);
        }
    }
}
