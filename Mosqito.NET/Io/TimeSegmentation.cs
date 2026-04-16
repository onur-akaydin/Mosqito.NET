using System.Threading.Tasks;

namespace Mosqito.Io;

/// <summary>
/// Segments a time signal into overlapping blocks.
/// Matches MoSQITo <c>utils/time_segmentation.py</c>.
/// </summary>
public static class TimeSegmentation
{
    /// <summary>
    /// Segments a 1-D signal into overlapping blocks.
    /// </summary>
    public static (double[,] blockArray, double[] timeAxis) Segment(
        ReadOnlySpan<double> sig,
        double fs,
        int nPerSeg = 2048,
        int? noOverlap = null,
        bool isEcma = false)
    {
        int noverlap = noOverlap ?? (nPerSeg / 2);
        if (noverlap == 0) noverlap = nPerSeg;

        // ECMA path: prepend nPerSeg zeros into a new array.
        // Non-ECMA: operate directly on the original span (no copy).
        double[]? ecmaBuf = null;
        ReadOnlySpan<double> src;
        if (isEcma)
        {
            ecmaBuf = new double[sig.Length + nPerSeg];
            sig.CopyTo(ecmaBuf.AsSpan(nPerSeg));
            src = ecmaBuf;
        }
        else
        {
            src = sig;
        }

        int n = src.Length;

        // Count blocks
        int nBlocks = 0;
        { int l = 0; while (l * noverlap <= n - nPerSeg) { nBlocks++; l++; } }

        if (nBlocks == 0)
            return (new double[nPerSeg, 0], Array.Empty<double>());

        double[,] blockArray = new double[nPerSeg, nBlocks];
        double[] timeArray   = new double[nBlocks];

        // Compute per-block mean time via closed-form: mean of [start .. start+nPerSeg-1] / fs
        // = (start + (nPerSeg-1)/2.0) / fs
        double halfSeg = (nPerSeg - 1) * 0.5;

        // Parallelize block extraction — each block writes to a disjoint column.
        // Pass src as array for closure (span cannot be captured).
        if (isEcma)
        {
            double[] srcArr = ecmaBuf!;
            Parallel.For(0, nBlocks, block =>
            {
                int start = block * noverlap;
                timeArray[block] = (start + halfSeg) / fs;
                for (int s = 0; s < nPerSeg; s++)
                    blockArray[s, block] = srcArr[start + s];
            });
        }
        else
        {
            // Pin the span via ToArray only if parallelism is worthwhile (≥4 blocks),
            // otherwise do a single serial pass to avoid the allocation.
            if (nBlocks >= 4)
            {
                double[] srcArr = sig.ToArray();
                Parallel.For(0, nBlocks, block =>
                {
                    int start = block * noverlap;
                    timeArray[block] = (start + halfSeg) / fs;
                    for (int s = 0; s < nPerSeg; s++)
                        blockArray[s, block] = srcArr[start + s];
                });
            }
            else
            {
                for (int block = 0; block < nBlocks; block++)
                {
                    int start = block * noverlap;
                    timeArray[block] = (start + halfSeg) / fs;
                    for (int s = 0; s < nPerSeg; s++)
                        blockArray[s, block] = src[start + s];
                }
            }
        }

        return (blockArray, timeArray);
    }

    /// <summary>
    /// Returns each segment as a separate <c>double[]</c> rather than a 2-D block.
    /// Convenience overload used by loudness_zwst_perseg and sharpness_din_perseg.
    /// </summary>
    public static (IReadOnlyList<double[]> segments, double[] timeAxis) SegmentList(
        ReadOnlySpan<double> sig,
        double fs,
        int nPerSeg = 2048,
        int? noOverlap = null,
        bool isEcma = false)
    {
        var (blockArray, timeAxis) = Segment(sig, fs, nPerSeg, noOverlap, isEcma);
        int nBlocks  = blockArray.GetLength(1);
        int nSamples = blockArray.GetLength(0);

        var list = new List<double[]>(nBlocks);
        for (int b = 0; b < nBlocks; b++)
        {
            double[] seg = new double[nSamples];
            for (int s = 0; s < nSamples; s++) seg[s] = blockArray[s, b];
            list.Add(seg);
        }

        return (list, timeAxis);
    }
}
