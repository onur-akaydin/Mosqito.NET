using System.Collections.Concurrent;
using System.Runtime.CompilerServices;

namespace Mosqito.Dsp;

/// <summary>
/// Window functions used throughout Mosqito.NET.
/// All Fill* methods write into a caller-supplied span to avoid allocation.
/// Allocating wrappers (Hann, Hamming, Blackman, VonHannEcma) are cached — callers
/// must treat the returned arrays as <b>read-only</b>; never mutate them.
/// </summary>
public static class Windows
{
    // Cached by (length, kind-tag). Returned arrays are immutable shared instances.
    private static readonly ConcurrentDictionary<(int n, int kind), double[]> _cache = new();

    private static double[] Cached(int n, int kind, Func<int, double[]> factory) =>
        _cache.GetOrAdd((n, kind), static (key, f) => f(key.n), factory);
    // ------------------------------------------------------------------
    // Hann (Hanning)  — matches numpy.hanning(N)
    // w[n] = 0.5 * (1 - cos(2π*n/(N-1)))
    // ------------------------------------------------------------------

    /// <summary>Fills <paramref name="output"/> with a Hann window of the same length.</summary>
    public static void FillHann(Span<double> output)
    {
        int n = output.Length;
        if (n == 1) { output[0] = 1.0; return; }
        double scale = 2.0 * Math.PI / (n - 1);
        for (int i = 0; i < n; i++)
            output[i] = 0.5 * (1.0 - Math.Cos(scale * i));
    }

    /// <summary>Returns a cached Hann window of length <paramref name="n"/> (read-only).</summary>
    public static double[] Hann(int n) => Cached(n, 0, static len => { var w = new double[len]; FillHann(w); return w; });

    // ------------------------------------------------------------------
    // Hamming  — matches numpy.hamming(N)
    // w[n] = 0.54 - 0.46 * cos(2π*n/(N-1))
    // ------------------------------------------------------------------

    /// <summary>Fills <paramref name="output"/> with a Hamming window.</summary>
    public static void FillHamming(Span<double> output)
    {
        int n = output.Length;
        if (n == 1) { output[0] = 1.0; return; }
        double scale = 2.0 * Math.PI / (n - 1);
        for (int i = 0; i < n; i++)
            output[i] = 0.54 - 0.46 * Math.Cos(scale * i);
    }

    /// <summary>Returns a cached Hamming window of length <paramref name="n"/> (read-only).</summary>
    public static double[] Hamming(int n) => Cached(n, 1, static len => { var w = new double[len]; FillHamming(w); return w; });

    // ------------------------------------------------------------------
    // Blackman  — matches numpy.blackman(N)
    // w[n] = 0.42 - 0.5*cos(2π*n/(N-1)) + 0.08*cos(4π*n/(N-1))
    // ------------------------------------------------------------------

    /// <summary>Fills <paramref name="output"/> with a Blackman window.</summary>
    public static void FillBlackman(Span<double> output)
    {
        int n = output.Length;
        if (n == 1) { output[0] = 1.0; return; }
        double s1 = 2.0 * Math.PI / (n - 1);
        double s2 = 4.0 * Math.PI / (n - 1);
        for (int i = 0; i < n; i++)
            output[i] = 0.42 - 0.5 * Math.Cos(s1 * i) + 0.08 * Math.Cos(s2 * i);
    }

    /// <summary>Returns a cached Blackman window of length <paramref name="n"/> (read-only).</summary>
    public static double[] Blackman(int n) => Cached(n, 2, static len => { var w = new double[len]; FillBlackman(w); return w; });

    // ------------------------------------------------------------------
    // Rectangular (boxcar)
    // ------------------------------------------------------------------

    /// <summary>Fills <paramref name="output"/> with ones (rectangular window).</summary>
    public static void FillRectangular(Span<double> output) => output.Fill(1.0);

    /// <summary>Allocates and returns a rectangular window of length <paramref name="n"/>.</summary>
    public static double[] Rectangular(int n)
    {
        double[] w = new double[n];
        w.AsSpan().Fill(1.0);
        return w;
    }

    // ------------------------------------------------------------------
    // Tukey (cosine-tapered)
    // w[n] = 0.5*(1 + cos(π*(2n/α/(N-1) - 1))) for n ≤ α*(N-1)/2
    //      = 1.0                                  for α*(N-1)/2 < n < (N-1)*(1-α/2)
    //      = 0.5*(1 + cos(π*(2n/α/(N-1) - 2/α + 1))) otherwise
    // ------------------------------------------------------------------

    /// <summary>Fills <paramref name="output"/> with a Tukey (cosine-tapered) window.</summary>
    /// <param name="alpha">Fraction of the window inside the taper (0 = rectangular, 1 = Hann).</param>
    public static void FillTukey(Span<double> output, double alpha = 0.5)
    {
        int n = output.Length;
        if (n == 1) { output[0] = 1.0; return; }
        if (alpha <= 0.0) { FillRectangular(output); return; }
        if (alpha >= 1.0) { FillHann(output); return; }

        double nm1 = n - 1;
        int width = (int)(alpha * nm1 / 2.0);
        for (int i = 0; i < n; i++)
        {
            if (i <= width)
                output[i] = 0.5 * (1.0 + Math.Cos(Math.PI * (2.0 * i / (alpha * nm1) - 1.0)));
            else if (i >= n - 1 - width)
                output[i] = 0.5 * (1.0 + Math.Cos(Math.PI * (2.0 * i / (alpha * nm1) - 2.0 / alpha + 1.0)));
            else
                output[i] = 1.0;
        }
    }

    // ------------------------------------------------------------------
    // Von-Hann variant used in ECMA-418-2 roughness (_von_hann_window.py)
    // w[n] = 0.5 * (1 - cos(2π*(n+1)/(N+1)))   (one-indexed, ECMA style)
    // ------------------------------------------------------------------

    /// <summary>
    /// ECMA-418-2 Von Hann window variant (one-indexed formula).
    /// w[n] = 0.5 * (1 - cos(2π*(n+1)/(N+1))).
    /// </summary>
    public static void FillVonHannEcma(Span<double> output)
    {
        int n = output.Length;
        double denom = n + 1.0;
        for (int i = 0; i < n; i++)
            output[i] = 0.5 * (1.0 - Math.Cos(2.0 * Math.PI * (i + 1) / denom));
    }

    /// <summary>Returns a cached ECMA-418-2 Von Hann window of length <paramref name="n"/> (read-only).</summary>
    public static double[] VonHannEcma(int n) => Cached(n, 3, static len => { var w = new double[len]; FillVonHannEcma(w); return w; });

    // ------------------------------------------------------------------
    // Helper: apply a window vector elementwise to a signal span (in-place).
    // ------------------------------------------------------------------

    /// <summary>
    /// Multiplies <paramref name="signal"/> elementwise by <paramref name="window"/> in-place.
    /// Both spans must have the same length.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static void Apply(Span<double> signal, ReadOnlySpan<double> window)
    {
        if (signal.Length != window.Length)
            throw new ArgumentException("Signal and window lengths must match.");
        VectorMath.Multiply(signal, window, signal);
    }

    /// <summary>
    /// Computes and returns the sum of the window values (for amplitude correction).
    /// </summary>
    public static double Sum(ReadOnlySpan<double> window)
    {
        double s = 0.0;
        for (int i = 0; i < window.Length; i++) s += window[i];
        return s;
    }

    /// <summary>
    /// Normalises <paramref name="window"/> by dividing every element by the window sum.
    /// This replicates MoSQITo's <c>window = window / sum(window)</c> amplitude correction.
    /// </summary>
    public static void NormaliseBySumInPlace(Span<double> window)
    {
        double s = Sum(window);
        if (s == 0.0) return;
        double inv = 1.0 / s;
        for (int i = 0; i < window.Length; i++) window[i] *= inv;
    }
}
