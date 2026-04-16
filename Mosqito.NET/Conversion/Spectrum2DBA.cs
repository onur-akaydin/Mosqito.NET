using System.Runtime.CompilerServices;

namespace Mosqito.Conversion;

/// <summary>
/// A-weighting of a linear-frequency spectrum according to IEC 61672:2014.
/// Ported from MoSQITo <c>spectrum2dBA</c>.
/// </summary>
public static class Spectrum2DBA
{
    // IEC 61672-1:2014 Table 1 — A-weighting corrections (dB) at 1/3-octave centre frequencies.
    private static readonly double[] AStandard =
    {
        -70.4, -63.4, -56.7, -50.5, -44.7, -39.4, -34.6, -30.2,
        -26.2, -22.5, -19.1, -16.1, -13.4, -10.9,  -8.6,  -6.6,
         -4.8,  -3.2,  -1.9,  -0.8,   0.0,   0.6,   1.0,   1.2,
          1.3,   1.2,   1.0,   0.5,  -0.1,  -1.1,  -2.5,  -4.3,
         -6.6,  -9.3,
    };

    private static readonly double[] FreqStandard =
    {
           10,  12.5,   16,   20,   25,  31.5,   40,   50,
           63,    80,  100,  125,  160,  200,  250,  315,
          400,   500,  630,  800, 1000, 1250, 1600, 2000,
         2500,  3150, 4000, 5000, 6300, 8000,10000,12500,
        16000, 20000,
    };

    /// <summary>
    /// Applies A-weighting to a spectrum in dB.
    /// </summary>
    /// <param name="spectrum">
    /// Input spectrum in dB, length N. Frequency axis is assumed to span linearly
    /// from 0 Hz to fs/2 Hz across <c>N</c> bins (matching numpy's
    /// <c>linspace(0, fs/2, N)</c>).
    /// </param>
    /// <param name="fs">Sampling frequency in Hz.</param>
    /// <returns>A-weighted spectrum in dB, same length as <paramref name="spectrum"/>.</returns>
    public static double[] Apply(ReadOnlySpan<double> spectrum, int fs)
    {
        double[] result = new double[spectrum.Length];
        Apply(spectrum, fs, result);
        return result;
    }

    /// <summary>
    /// Applies A-weighting to a spectrum in dB, writing into <paramref name="output"/>.
    /// </summary>
    public static void Apply(ReadOnlySpan<double> spectrum, int fs, Span<double> output)
    {
        if (output.Length < spectrum.Length)
            throw new ArgumentException("Output span is too short.", nameof(output));

        int n = spectrum.Length;
        double fMax = fs * 0.5;

        for (int i = 0; i < n; i++)
        {
            double freq = fMax * i / (n - 1 == 0 ? 1 : n - 1);
            double aPond = InterpolateA(freq);
            output[i] = spectrum[i] + aPond;
        }
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double InterpolateA(double freq)
    {
        // numpy.interp equivalent: clamp at boundaries, piecewise linear.
        if (freq <= FreqStandard[0]) return AStandard[0];
        if (freq >= FreqStandard[FreqStandard.Length - 1]) return AStandard[AStandard.Length - 1];

        int lo = 0, hi = FreqStandard.Length - 2;
        while (lo < hi)
        {
            int mid = (lo + hi) >> 1;
            if (FreqStandard[mid + 1] <= freq) lo = mid + 1;
            else hi = mid;
        }
        double t = (freq - FreqStandard[lo]) / (FreqStandard[lo + 1] - FreqStandard[lo]);
        return AStandard[lo] + t * (AStandard[lo + 1] - AStandard[lo]);
    }
}
