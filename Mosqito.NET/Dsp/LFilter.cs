using System.Buffers;
using System.Runtime.CompilerServices;

namespace Mosqito.Dsp;

/// <summary>
/// Direct-form IIR filter (scipy.signal.lfilter equivalent).
/// Coefficients are given as separate numerator <c>b</c> and denominator <c>a</c> arrays.
/// </summary>
public static class LFilter
{
    /// <summary>
    /// Applies a causal IIR filter defined by <c>b/a</c> polynomial coefficients.
    /// Matches scipy.signal.lfilter(b, a, x).
    /// </summary>
    /// <param name="b">Numerator polynomial coefficients (length M+1).</param>
    /// <param name="a">Denominator polynomial coefficients (length N+1, a[0] should be 1).</param>
    /// <param name="x">Input signal.</param>
    /// <returns>Filtered output, same length as x.</returns>
    public static double[] Apply(ReadOnlySpan<double> b, ReadOnlySpan<double> a,
        ReadOnlySpan<double> x)
    {
        double[] y = new double[x.Length];
        Apply(b, a, x, y);
        return y;
    }

    /// <summary>
    /// Applies a causal IIR filter, writing into <paramref name="output"/>.
    /// </summary>
    public static void Apply(ReadOnlySpan<double> b, ReadOnlySpan<double> a,
        ReadOnlySpan<double> x, Span<double> output)
    {
        int nb = b.Length, na = a.Length;
        int n = x.Length;
        if (output.Length < n) throw new ArgumentException("Output too short.", nameof(output));

        // Normalise by a[0]
        double a0 = a[0];
        if (a0 == 0.0) throw new ArgumentException("a[0] must not be zero.", nameof(a));

        int stateLen = Math.Max(nb, na) - 1;
        Span<double> zi = stateLen <= 64
            ? stackalloc double[stateLen]
            : new double[stateLen]; // rare; avoids alloc for large filters via ArrayPool in hot paths
        zi.Clear();

        ApplyCore(b, a, a0, x, output, zi, nb, na, n);
    }

    /// <summary>
    /// Applies IIR filter with an explicit initial state vector <paramref name="zi"/>
    /// (length max(nb,na)-1). Updates zi in-place with the final state.
    /// </summary>
    public static void Apply(ReadOnlySpan<double> b, ReadOnlySpan<double> a,
        ReadOnlySpan<double> x, Span<double> output, Span<double> zi)
    {
        int nb = b.Length, na = a.Length;
        int n = x.Length;
        double a0 = a[0];
        ApplyCore(b, a, a0, x, output, zi, nb, na, n);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void ApplyCore(ReadOnlySpan<double> b, ReadOnlySpan<double> a, double a0,
        ReadOnlySpan<double> x, Span<double> y, Span<double> zi, int nb, int na, int n)
    {
        double invA0 = 1.0 / a0;
        int stateLen = zi.Length;

        for (int i = 0; i < n; i++)
        {
            double xi = x[i];
            double yi = b[0] * xi + (stateLen > 0 ? zi[0] : 0.0);
            yi *= invA0;
            y[i] = yi;

            // Update state
            for (int j = 0; j < stateLen - 1; j++)
            {
                double bj = j + 1 < nb ? b[j + 1] : 0.0;
                double aj = j + 1 < na ? a[j + 1] : 0.0;
                zi[j] = bj * xi - aj * yi + zi[j + 1];
            }
            if (stateLen > 0)
            {
                double bLast = stateLen < nb ? b[stateLen] : 0.0;
                double aLast = stateLen < na ? a[stateLen] : 0.0;
                zi[stateLen - 1] = bLast * xi - aLast * yi;
            }
        }
    }

    /// <summary>
    /// Zero-phase IIR filter (forward + backward pass, matches scipy.signal.filtfilt).
    /// </summary>
    public static double[] FiltFilt(ReadOnlySpan<double> b, ReadOnlySpan<double> a,
        ReadOnlySpan<double> x)
    {
        int n = x.Length;
        int pad = Math.Min(3 * Math.Max(b.Length, a.Length), n - 1);
        if (pad < 1) pad = 1;
        int nPad = n + 2 * pad;

        double[] padded  = ArrayPool<double>.Shared.Rent(nPad);
        double[] forward = ArrayPool<double>.Shared.Rent(nPad);
        try
        {
            for (int i = 0; i < pad; i++)  padded[i] = 2.0 * x[0] - x[pad - i];
            for (int i = 0; i < n; i++)    padded[pad + i] = x[i];
            for (int i = 0; i < pad; i++)  padded[pad + n + i] = 2.0 * x[n - 1] - x[n - 2 - i];

            Apply(b, a, padded.AsSpan(0, nPad), forward.AsSpan(0, nPad));

            // Reverse and filter again
            for (int i = 0; i < nPad / 2; i++)
                (forward[i], forward[nPad - 1 - i]) = (forward[nPad - 1 - i], forward[i]);

            double[] backward = new double[nPad];
            Apply(b, a, forward.AsSpan(0, nPad), backward.AsSpan());

            double[] result = new double[n];
            for (int i = 0; i < n; i++)
                result[i] = backward[nPad - 1 - (pad + n - 1 - i)];
            return result;
        }
        finally
        {
            ArrayPool<double>.Shared.Return(padded);
            ArrayPool<double>.Shared.Return(forward);
        }
    }
}
