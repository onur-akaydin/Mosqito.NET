using System.Numerics;
using System.Runtime.CompilerServices;

namespace Mosqito.Dsp;

/// <summary>
/// SIMD-accelerated element-wise math operations on double arrays.
/// Uses System.Numerics.Vector&lt;double&gt; for vectorised inner loops.
/// </summary>
public static class VectorMath
{
    private static readonly int VWidth = Vector<double>.Count;

    // ------------------------------------------------------------------
    // Trapezoid integration  (numpy.trapz)
    // ------------------------------------------------------------------

    /// <summary>
    /// Numerical integration using the trapezoidal rule.
    /// Matches numpy.trapz(y, x=None, dx=1.0).
    /// </summary>
    public static double Trapz(ReadOnlySpan<double> y, double dx = 1.0)
    {
        int n = y.Length;
        if (n < 2) return 0.0;
        double sum = 0.0;
        for (int i = 0; i < n - 1; i++)
            sum += (y[i] + y[i + 1]) * 0.5 * dx;
        return sum;
    }

    /// <summary>
    /// Trapezoidal integration with explicit x axis.
    /// Matches numpy.trapz(y, x).
    /// </summary>
    public static double Trapz(ReadOnlySpan<double> y, ReadOnlySpan<double> x)
    {
        if (y.Length != x.Length) throw new ArgumentException("y and x must have the same length.");
        int n = y.Length;
        if (n < 2) return 0.0;
        double sum = 0.0;
        for (int i = 0; i < n - 1; i++)
            sum += (y[i] + y[i + 1]) * 0.5 * (x[i + 1] - x[i]);
        return sum;
    }

    // ------------------------------------------------------------------
    // Cumulative sum
    // ------------------------------------------------------------------

    /// <summary>Computes the cumulative sum of <paramref name="x"/> (numpy.cumsum).</summary>
    public static double[] Cumsum(ReadOnlySpan<double> x)
    {
        double[] result = new double[x.Length];
        double acc = 0.0;
        for (int i = 0; i < x.Length; i++) { acc += x[i]; result[i] = acc; }
        return result;
    }

    // ------------------------------------------------------------------
    // Element-wise squared SIMD
    // ------------------------------------------------------------------

    /// <summary>Computes element-wise x^2, writing into <paramref name="output"/>.</summary>
    public static void Square(ReadOnlySpan<double> x, Span<double> output)
    {
        if (output.Length < x.Length) throw new ArgumentException("Output too short.");
        int i = 0;
        int vLen = x.Length - x.Length % VWidth;
        for (; i < vLen; i += VWidth)
        {
            var v = new Vector<double>(x.Slice(i, VWidth));
            (v * v).CopyTo(output.Slice(i, VWidth));
        }
        for (; i < x.Length; i++) output[i] = x[i] * x[i];
    }

    /// <summary>Returns element-wise x^2 as a new array.</summary>
    public static double[] Square(ReadOnlySpan<double> x)
    {
        double[] result = new double[x.Length];
        Square(x, result);
        return result;
    }

    // ------------------------------------------------------------------
    // Element-wise sqrt
    // ------------------------------------------------------------------

    /// <summary>Computes element-wise sqrt(x), writing into <paramref name="output"/>.</summary>
    public static void Sqrt(ReadOnlySpan<double> x, Span<double> output)
    {
        if (output.Length < x.Length) throw new ArgumentException("Output too short.");
        // Vector<double> sqrt not available universally; scalar loop
        for (int i = 0; i < x.Length; i++) output[i] = Math.Sqrt(x[i]);
    }

    /// <summary>Returns element-wise sqrt as a new array.</summary>
    public static double[] Sqrt(ReadOnlySpan<double> x)
    {
        double[] result = new double[x.Length];
        Sqrt(x, result);
        return result;
    }

    // ------------------------------------------------------------------
    // Element-wise multiply (SIMD)
    // ------------------------------------------------------------------

    /// <summary>Element-wise a * b → output (SIMD).</summary>
    public static void Multiply(ReadOnlySpan<double> a, ReadOnlySpan<double> b, Span<double> output)
    {
        if (a.Length != b.Length) throw new ArgumentException("Lengths must match.");
        if (output.Length < a.Length) throw new ArgumentException("Output too short.");
        int i = 0, vLen = a.Length - a.Length % VWidth;
        for (; i < vLen; i += VWidth)
            (new Vector<double>(a.Slice(i, VWidth)) * new Vector<double>(b.Slice(i, VWidth)))
                .CopyTo(output.Slice(i, VWidth));
        for (; i < a.Length; i++) output[i] = a[i] * b[i];
    }

    /// <summary>Element-wise a + b → output (SIMD).</summary>
    public static void Add(ReadOnlySpan<double> a, ReadOnlySpan<double> b, Span<double> output)
    {
        if (a.Length != b.Length) throw new ArgumentException("Lengths must match.");
        if (output.Length < a.Length) throw new ArgumentException("Output too short.");
        int i = 0, vLen = a.Length - a.Length % VWidth;
        for (; i < vLen; i += VWidth)
            (new Vector<double>(a.Slice(i, VWidth)) + new Vector<double>(b.Slice(i, VWidth)))
                .CopyTo(output.Slice(i, VWidth));
        for (; i < a.Length; i++) output[i] = a[i] + b[i];
    }

    /// <summary>Scalar multiply: a * scalar → output (SIMD).</summary>
    public static void Scale(ReadOnlySpan<double> a, double scalar, Span<double> output)
    {
        if (output.Length < a.Length) throw new ArgumentException("Output too short.");
        var sv = new Vector<double>(scalar);
        int i = 0, vLen = a.Length - a.Length % VWidth;
        for (; i < vLen; i += VWidth)
            (new Vector<double>(a.Slice(i, VWidth)) * sv).CopyTo(output.Slice(i, VWidth));
        for (; i < a.Length; i++) output[i] = a[i] * scalar;
    }

    // ------------------------------------------------------------------
    // RMS
    // ------------------------------------------------------------------

    /// <summary>Root-mean-square of <paramref name="x"/>.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double Rms(ReadOnlySpan<double> x)
    {
        double sumSq = 0.0;
        for (int i = 0; i < x.Length; i++) sumSq += x[i] * x[i];
        return Math.Sqrt(sumSq / x.Length);
    }

    // ------------------------------------------------------------------
    // Max / ArgMax / ArgMin
    // ------------------------------------------------------------------

    /// <summary>Returns the index of the maximum value (matches numpy.argmax).</summary>
    public static int ArgMax(ReadOnlySpan<double> x)
    {
        int best = 0;
        for (int i = 1; i < x.Length; i++)
            if (x[i] > x[best]) best = i;
        return best;
    }

    /// <summary>Returns the index of the minimum value (matches numpy.argmin).</summary>
    public static int ArgMin(ReadOnlySpan<double> x)
    {
        int best = 0;
        for (int i = 1; i < x.Length; i++)
            if (x[i] < x[best]) best = i;
        return best;
    }

    // ------------------------------------------------------------------
    // Magnitude of complex array
    // ------------------------------------------------------------------

    /// <summary>Computes |x| for each element of a complex array.</summary>
    public static double[] Abs(ReadOnlySpan<System.Numerics.Complex> x)
    {
        double[] result = new double[x.Length];
        for (int i = 0; i < x.Length; i++) result[i] = x[i].Magnitude;
        return result;
    }

    // ------------------------------------------------------------------
    // Log10 (element-wise)
    // ------------------------------------------------------------------

    /// <summary>Computes log10 for each element.</summary>
    public static double[] Log10(ReadOnlySpan<double> x)
    {
        double[] result = new double[x.Length];
        for (int i = 0; i < x.Length; i++) result[i] = Math.Log10(x[i]);
        return result;
    }

    // ------------------------------------------------------------------
    // Clip / clamp
    // ------------------------------------------------------------------

    /// <summary>Clamps all values in <paramref name="x"/> to [min, max] in-place.</summary>
    public static void ClipInPlace(Span<double> x, double min, double max)
    {
        for (int i = 0; i < x.Length; i++)
            x[i] = Math.Max(min, Math.Min(max, x[i]));
    }

    // ------------------------------------------------------------------
    // Convolve (numpy.convolve, full mode, no FFT — for short kernels)
    // ------------------------------------------------------------------

    /// <summary>
    /// Full linear convolution. For short kernels (len ≤ 256); beyond that use FFT.
    /// </summary>
    public static double[] Convolve(ReadOnlySpan<double> a, ReadOnlySpan<double> b)
    {
        int na = a.Length, nb = b.Length;
        int nOut = na + nb - 1;
        double[] result = new double[nOut];
        for (int i = 0; i < na; i++)
            for (int j = 0; j < nb; j++)
                result[i + j] += a[i] * b[j];
        return result;
    }

    // ------------------------------------------------------------------
    // Min / Max of arrays
    // ------------------------------------------------------------------

    /// <summary>Element-wise minimum of two arrays.</summary>
    public static double[] Minimum(ReadOnlySpan<double> a, ReadOnlySpan<double> b)
    {
        if (a.Length != b.Length) throw new ArgumentException("Lengths must match.");
        double[] result = new double[a.Length];
        for (int i = 0; i < a.Length; i++) result[i] = Math.Min(a[i], b[i]);
        return result;
    }

    /// <summary>Element-wise maximum of two arrays.</summary>
    public static double[] Maximum(ReadOnlySpan<double> a, ReadOnlySpan<double> b)
    {
        if (a.Length != b.Length) throw new ArgumentException("Lengths must match.");
        double[] result = new double[a.Length];
        for (int i = 0; i < a.Length; i++) result[i] = Math.Max(a[i], b[i]);
        return result;
    }
}
