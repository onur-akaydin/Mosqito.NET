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
    /// </summary>
    /// <param name="x">Input 1-D array.</param>
    /// <param name="height">Minimum peak height (inclusive). Null = no constraint.</param>
    /// <param name="distance">Minimum horizontal distance between peaks. Null = no constraint.</param>
    /// <param name="prominence">Minimum prominence. Null = no constraint.</param>
    public static PeakResult Find(
        ReadOnlySpan<double> x,
        double? height = null,
        int? distance = null,
        double? prominence = null)
    {
        // Copy to array so it can be captured in lambdas and passed to helper methods
        double[] xArr = x.ToArray();
        int n = xArr.Length;
        var candidates = new List<int>();

        // Find all local maxima (strictly greater than both neighbours)
        for (int i = 1; i < n - 1; i++)
        {
            if (xArr[i] > xArr[i - 1] && xArr[i] > xArr[i + 1])
                candidates.Add(i);
        }

        // Filter by height
        if (height.HasValue)
        {
            double h = height.Value;
            candidates.RemoveAll(i => xArr[i] < h);
        }

        // Filter by minimum distance (keep the tallest peak within each distance window)
        if (distance.HasValue && distance.Value > 1 && candidates.Count > 1)
        {
            candidates = FilterByDistance(xArr, candidates, distance.Value);
        }

        // Compute prominences and optionally filter
        double[] proms = ComputeProminences(xArr, candidates);

        if (prominence.HasValue)
        {
            double minProm = prominence.Value;
            var filtered = new List<int>();
            var filteredProms = new List<double>();
            for (int i = 0; i < candidates.Count; i++)
            {
                if (proms[i] >= minProm)
                {
                    filtered.Add(candidates[i]);
                    filteredProms.Add(proms[i]);
                }
            }
            candidates = filtered;
            proms = filteredProms.ToArray();
        }

        int[] idxArr = candidates.ToArray();
        double[] heights = new double[idxArr.Length];
        for (int i = 0; i < idxArr.Length; i++) heights[i] = xArr[idxArr[i]];

        return new PeakResult { Indices = idxArr, Heights = heights, Prominences = proms };
    }

    // ------------------------------------------------------------------
    // Prominences  (matches scipy.signal.peak_prominences)
    // ------------------------------------------------------------------

    /// <summary>
    /// Computes the prominence of each peak at the given indices.
    /// Prominence = peak height - the highest minimum within the surrounding contour bases.
    /// </summary>
    public static double[] ComputeProminences(double[] x, List<int> peakIndices)
    {
        int n = x.Length;
        double[] proms = new double[peakIndices.Count];

        for (int pi = 0; pi < peakIndices.Count; pi++)
        {
            int idx = peakIndices[pi];
            double peakHeight = x[idx];

            // Find left base: walk left until signal ≥ peak or array boundary
            double leftMin = peakHeight;
            int left = idx - 1;
            while (left >= 0)
            {
                if (x[left] > peakHeight) break;
                if (x[left] < leftMin) leftMin = x[left];
                left--;
            }

            // Find right base: walk right
            double rightMin = peakHeight;
            int right = idx + 1;
            while (right < n)
            {
                if (x[right] > peakHeight) break;
                if (x[right] < rightMin) rightMin = x[right];
                right++;
            }

            double baseHeight = Math.Max(leftMin, rightMin);
            proms[pi] = peakHeight - baseHeight;
        }

        return proms;
    }

    // ------------------------------------------------------------------
    // Distance filtering (keep tallest in each window)
    // ------------------------------------------------------------------
    private static List<int> FilterByDistance(
        double[] x, List<int> peaks, int distance)
    {
        bool[] keep = new bool[peaks.Count];
        keep.AsSpan().Fill(true);

        for (int i = 0; i < peaks.Count; i++)
        {
            if (!keep[i]) continue;
            for (int j = i + 1; j < peaks.Count; j++)
            {
                if (peaks[j] - peaks[i] >= distance) break;
                if (x[peaks[i]] >= x[peaks[j]]) keep[j] = false;
                else { keep[i] = false; break; }
            }
        }

        var result = new List<int>();
        for (int i = 0; i < peaks.Count; i++)
            if (keep[i]) result.Add(peaks[i]);
        return result;
    }
}
