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
    /// <param name="sig">Input signal (1-D).</param>
    /// <param name="fs">Sampling frequency [Hz].</param>
    /// <param name="nPerSeg">Block length in samples. Default 2048.</param>
    /// <param name="noOverlap">
    /// Number of overlapping samples between adjacent blocks.
    /// If null, defaults to <c>nPerSeg / 2</c>. Set to 0 for no overlap (step = nPerSeg).
    /// </param>
    /// <param name="isEcma">
    /// If true, zero-pads <c>nPerSeg</c> samples at the start (ECMA-418-2, Formula 16).
    /// </param>
    /// <returns>
    /// (<c>blockArray</c> of shape [nPerSeg, nSegments], <c>timeAxis</c> of length nSegments).
    /// </returns>
    public static (double[,] blockArray, double[] timeAxis) Segment(
        ReadOnlySpan<double> sig,
        double fs,
        int nPerSeg = 2048,
        int? noOverlap = null,
        bool isEcma = false)
    {
        int noverlap = noOverlap ?? (nPerSeg / 2);
        if (noverlap == 0) noverlap = nPerSeg; // step = nPerSeg (Python behaviour for noverlap=0)

        // ECMA path: prepend nPerSeg zeros
        double[] sigArr;
        if (isEcma)
        {
            sigArr = new double[sig.Length + nPerSeg];
            sig.CopyTo(sigArr.AsSpan(nPerSeg));
        }
        else
        {
            sigArr = sig.ToArray();
        }
        int n = sigArr.Length;

        // Build time axis for the (possibly padded) signal
        double[] time = new double[n];
        for (int i = 0; i < n; i++) time[i] = i / fs;

        // Count blocks
        int nBlocks = 0;
        int l = 0;
        while (l * noverlap <= n - nPerSeg) { nBlocks++; l++; }

        if (nBlocks == 0)
            return (new double[nPerSeg, 0], Array.Empty<double>());

        double[,] blockArray = new double[nPerSeg, nBlocks];
        double[] timeArray   = new double[nBlocks];

        for (int block = 0; block < nBlocks; block++)
        {
            int start = block * noverlap;
            for (int s = 0; s < nPerSeg; s++)
                blockArray[s, block] = sigArr[start + s];

            // Mean of the time values for this block
            double tSum = 0.0;
            for (int s = 0; s < nPerSeg; s++) tSum += time[start + s];
            timeArray[block] = tSum / nPerSeg;
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
