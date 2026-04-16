using System.Buffers;
using System.Runtime.CompilerServices;
using System.Threading.Tasks;

namespace Mosqito.Dsp;

/// <summary>
/// Second-order-section (biquad cascade) IIR filter processing.
/// SOS coefficients are [b0, b1, b2, a0, a1, a2] per row (a0 normalised to 1 by convention).
///
/// Hot path: Direct-Form-II-Transposed biquad cascade with zero-allocation processing.
/// Ported from the scipy sosfilt / sosfiltfilt behaviour used throughout MoSQITo.
/// </summary>
public static class SosFilter
{
    // ------------------------------------------------------------------
    // Forward filtering  (sosfilt)
    // ------------------------------------------------------------------

    /// <summary>
    /// Applies an SOS filter to <paramref name="input"/>, writing the result to
    /// <paramref name="output"/>. Filter state is initialised to zero.
    /// </summary>
    /// <param name="sos">SOS coefficients, shape (nSections, 6): [b0,b1,b2, a0,a1,a2].</param>
    /// <param name="input">Input signal.</param>
    /// <param name="output">Output signal (same length as input).</param>
    public static void Process(double[,] sos, ReadOnlySpan<double> input, Span<double> output)
    {
        int nSec = sos.GetLength(0);
        int n = input.Length;
        if (output.Length < n) throw new ArgumentException("Output too short.", nameof(output));

        // Allocate filter state on the stack for small section counts, pool for large.
        if (nSec <= 32)
        {
            Span<double> zi = stackalloc double[nSec * 2];
            zi.Clear();
            ProcessCore(sos, input, output, zi, nSec);
        }
        else
        {
            double[] ziArr = ArrayPool<double>.Shared.Rent(nSec * 2);
            ziArr.AsSpan(0, nSec * 2).Clear();
            ProcessCore(sos, input, output, ziArr, nSec);
            ArrayPool<double>.Shared.Return(ziArr);
        }
    }

    /// <summary>
    /// Applies an SOS filter to <paramref name="input"/>, returning the result.
    /// </summary>
    public static double[] Process(double[,] sos, ReadOnlySpan<double> input)
    {
        double[] result = new double[input.Length];
        Process(sos, input, result);
        return result;
    }

    /// <summary>
    /// Applies an SOS filter with an explicit initial state <paramref name="zi"/> and
    /// writes the final state back. zi has shape (nSections * 2).
    /// </summary>
    public static void ProcessZi(double[,] sos, ReadOnlySpan<double> input,
        Span<double> output, Span<double> zi)
    {
        int nSec = sos.GetLength(0);
        if (zi.Length < nSec * 2) throw new ArgumentException("zi too short.", nameof(zi));
        if (output.Length < input.Length) throw new ArgumentException("Output too short.", nameof(output));
        ProcessCore(sos, input, output, zi, nSec);
    }

    // ------------------------------------------------------------------
    // Zero-phase forward-backward filtering  (sosfiltfilt)
    // ------------------------------------------------------------------

    /// <summary>
    /// Zero-phase IIR filter (forward + backward pass).
    /// Matches scipy.signal.sosfiltfilt(sos, x).
    /// Adds a padded transient region (3 * max(len(sos)) or 3*(nSec*2+1) samples).
    /// </summary>
    public static double[] FiltFilt(double[,] sos, ReadOnlySpan<double> input)
    {
        int nSec = sos.GetLength(0);
        int n = input.Length;
        // Determine the required edge padding (approx. 3 * filter order)
        int pad = Math.Min(3 * nSec * 2, n - 1);
        if (pad < 1) pad = 1;

        // Build padded signal: [reflect left | signal | reflect right]
        int nPad = n + 2 * pad;
        double[] padded = ArrayPool<double>.Shared.Rent(nPad);
        double[] forward = ArrayPool<double>.Shared.Rent(nPad);
        double[] backward = ArrayPool<double>.Shared.Rent(nPad);
        double[] ziArr = ArrayPool<double>.Shared.Rent(nSec * 2);

        try
        {
            // Reflect edges
            for (int i = 0; i < pad; i++)
                padded[i] = 2.0 * input[0] - input[pad - i];
            for (int i = 0; i < n; i++)
                padded[pad + i] = input[i];
            for (int i = 0; i < pad; i++)
                padded[pad + n + i] = 2.0 * input[n - 1] - input[n - 2 - i];

            // Forward pass
            ziArr.AsSpan(0, nSec * 2).Clear();
            ProcessCore(sos, padded.AsSpan(0, nPad), forward.AsSpan(0, nPad), ziArr, nSec);

            // Reverse
            forward.AsSpan(0, nPad).Reverse();

            // Backward pass
            ziArr.AsSpan(0, nSec * 2).Clear();
            ProcessCore(sos, forward.AsSpan(0, nPad), backward.AsSpan(0, nPad), ziArr, nSec);

            // Reverse the backward output and extract centre region.
            // After the backward pass the output is still in reversed order,
            // so element i of the original signal sits at index (n+pad-1-i).
            double[] result = new double[n];
            for (int i = 0; i < n; i++)
                result[i] = backward[n + pad - 1 - i];

            return result;
        }
        finally
        {
            ArrayPool<double>.Shared.Return(padded);
            ArrayPool<double>.Shared.Return(forward);
            ArrayPool<double>.Shared.Return(backward);
            ArrayPool<double>.Shared.Return(ziArr);
        }
    }

    // ------------------------------------------------------------------
    // Internal biquad cascade core (Direct-Form-II-Transposed)
    // ------------------------------------------------------------------

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void ProcessCore(double[,] sos, ReadOnlySpan<double> x, Span<double> y,
        Span<double> zi, int nSec)
    {
        int n = x.Length;
        for (int i = 0; i < n; i++)
        {
            double sample = x[i];
            for (int s = 0; s < nSec; s++)
            {
                int z0 = s * 2;
                double b0 = sos[s, 0], b1 = sos[s, 1], b2 = sos[s, 2];
                // a0 is always 1 (normalised); a1, a2 are the denominator coefficients
                double a1 = sos[s, 4], a2 = sos[s, 5];

                double yn = b0 * sample + zi[z0];
                zi[z0]     = b1 * sample - a1 * yn + zi[z0 + 1];
                zi[z0 + 1] = b2 * sample - a2 * yn;
                sample = yn;
            }
            y[i] = sample;
        }
    }

    // ------------------------------------------------------------------
    // 2-D signal support (axis=0): filter each column independently
    // ------------------------------------------------------------------

    /// <summary>
    /// Applies an SOS filter along axis 0 of a 2-D array (rows = samples, cols = segments).
    /// Columns are processed in parallel — each column is independent.
    /// Returns a new 2-D array of the same shape.
    /// </summary>
    public static double[,] Process2D(double[,] sos, double[,] input)
    {
        int nSamples = input.GetLength(0);
        int nCols    = input.GetLength(1);
        double[,] output = new double[nSamples, nCols];

        Parallel.For(0, nCols,
            () => (col: new double[nSamples], colOut: new double[nSamples]),
            (c, _, state) =>
            {
                var (col, colOut) = state;
                for (int r = 0; r < nSamples; r++) col[r] = input[r, c];
                Process(sos, col, colOut);
                for (int r = 0; r < nSamples; r++) output[r, c] = colOut[r];
                return state;
            },
            _ => { });

        return output;
    }
}
