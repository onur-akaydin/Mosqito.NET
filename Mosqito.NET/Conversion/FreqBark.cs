using System.Runtime.CompilerServices;

namespace Mosqito.Conversion;

/// <summary>
/// Conversions between frequency in Hz and the Bark psychoacoustic scale.
///
/// Reference: E. Zwicker, H. Fastl — Psychoacoustics. Springer, 1990, Table 6.1.
/// Coefficients are linearly interpolated from the table values (identical to MoSQITo
/// <c>freq2bark</c> / <c>bark2freq</c>).
/// </summary>
public static class FreqBark
{
    // Zwicker Table 6.1 — 49 control points mapping Bark (0, 0.5, 1.0, …, 24.0) → Hz.
    private static readonly double[] BarkAxis = BuildBarkAxis();
    private static readonly double[] HzAxis =
    {
            0, 50, 100, 150, 200, 250, 300, 350, 400, 450,
            510, 570, 630, 700, 770, 840, 920, 1000, 1080, 1170,
            1270, 1370, 1480, 1600, 1720, 1850, 2000, 2150, 2320, 2500,
            2700, 2900, 3150, 3400, 3700, 4000, 4400, 4800, 5300, 5800,
            6400, 7000, 7700, 8500, 9500, 10500, 12000, 13500, 15500, 20000,
    };

    private static double[] BuildBarkAxis()
    {
        // 0.0, 0.5, 1.0, …, 24.0  (49 values)
        double[] b = new double[50];
        for (int i = 0; i < 50; i++) b[i] = i * 0.5;
        return b;
    }

    /// <summary>Converts a single frequency in Hz to Bark.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double Freq2Bark(double hz) => Interp(hz, HzAxis, BarkAxis);

    /// <summary>Converts a single Bark value to frequency in Hz.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double Bark2Freq(double bark) => Interp(bark, BarkAxis, HzAxis);

    /// <summary>
    /// Converts an array of Hz values to Bark, writing into <paramref name="output"/>.
    /// </summary>
    public static void Freq2Bark(ReadOnlySpan<double> hz, Span<double> output)
    {
        if (output.Length < hz.Length) throw new ArgumentException("Output span is too short.", nameof(output));
        for (int i = 0; i < hz.Length; i++)
            output[i] = Freq2Bark(hz[i]);
    }

    /// <summary>Converts an array of Hz values to Bark, returning a new array.</summary>
    public static double[] Freq2Bark(ReadOnlySpan<double> hz)
    {
        double[] result = new double[hz.Length];
        Freq2Bark(hz, result);
        return result;
    }

    /// <summary>
    /// Converts an array of Bark values to Hz, writing into <paramref name="output"/>.
    /// </summary>
    public static void Bark2Freq(ReadOnlySpan<double> bark, Span<double> output)
    {
        if (output.Length < bark.Length) throw new ArgumentException("Output span is too short.", nameof(output));
        for (int i = 0; i < bark.Length; i++)
            output[i] = Bark2Freq(bark[i]);
    }

    /// <summary>Converts an array of Bark values to Hz, returning a new array.</summary>
    public static double[] Bark2Freq(ReadOnlySpan<double> bark)
    {
        double[] result = new double[bark.Length];
        Bark2Freq(bark, result);
        return result;
    }

    // Piecewise linear interpolation matching numpy.interp behaviour (clamps at boundaries).
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double Interp(double x, double[] xp, double[] yp)
    {
        if (x <= xp[0]) return yp[0];
        if (x >= xp[xp.Length - 1]) return yp[yp.Length - 1];

        // Binary search for the interval.
        int lo = 0, hi = xp.Length - 2;
        while (lo < hi)
        {
            int mid = (lo + hi) >> 1;
            if (xp[mid + 1] <= x) lo = mid + 1;
            else hi = mid;
        }
        double t = (x - xp[lo]) / (xp[lo + 1] - xp[lo]);
        return yp[lo] + t * (yp[lo + 1] - yp[lo]);
    }
}
