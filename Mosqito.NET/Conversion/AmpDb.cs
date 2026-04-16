using System.Numerics;
using System.Runtime.CompilerServices;

namespace Mosqito.Conversion;

/// <summary>
/// Conversions between linear amplitude and decibel (dB SPL) scales.
/// Reference: p_ref = 2e-5 Pa (ISO 1683).
/// </summary>
public static class AmpDb
{
    private const double Log10E = 0.4342944819032518; // 1/ln(10)
    private const double Ln10 = 2.302585092994046;

    /// <summary>
    /// Converts a single amplitude value to dB: 20 * log10(amp / ref).
    /// Amplitudes equal to zero are replaced with 2e-12 to avoid -Inf.
    /// </summary>
    /// <param name="amp">Linear amplitude value.</param>
    /// <param name="reference">Reference level (default 1.0).</param>
    /// <returns>Level in dB.</returns>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double Amp2Db(double amp, double reference = 1.0)
    {
        if (reference == 0.0) throw new ArgumentException("Reference must not be zero.", nameof(reference));
        if (amp == 0.0) amp = 2e-12;
        return 20.0 * Math.Log10(amp / reference);
    }

    /// <summary>
    /// Converts a single dB value to linear amplitude: ref * 10^(dB/20).
    /// </summary>
    private const double Ln10Over20 = 0.11512925464970229; // ln(10)/20

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double Db2Amp(double dB, double reference = 1.0)
    {
        if (reference == 0.0) throw new ArgumentException("Reference must not be zero.", nameof(reference));
        return Math.Exp(dB * Ln10Over20) * reference;
    }

    /// <summary>
    /// Converts an array of amplitudes to dB in-place, writing to <paramref name="output"/>.
    /// Uses SIMD vectorisation where the hardware vector width allows.
    /// </summary>
    /// <param name="amp">Input amplitude samples.</param>
    /// <param name="output">Pre-allocated output buffer, same length as <paramref name="amp"/>.</param>
    /// <param name="reference">Reference level (default 1.0).</param>
    public static void Amp2Db(ReadOnlySpan<double> amp, Span<double> output, double reference = 1.0)
    {
        if (reference == 0.0) throw new ArgumentException("Reference must not be zero.", nameof(reference));
        if (output.Length < amp.Length) throw new ArgumentException("Output span is too short.", nameof(output));

        // Scalar path — log is not available in System.Numerics.Vector; compute per element.
        double invRef = 1.0 / reference;
        for (int i = 0; i < amp.Length; i++)
        {
            double a = amp[i] == 0.0 ? 2e-12 : amp[i];
            output[i] = 20.0 * Math.Log10(a * invRef);
        }
    }

    /// <summary>
    /// Converts an array of amplitudes to dB, allocating and returning a new array.
    /// </summary>
    public static double[] Amp2Db(ReadOnlySpan<double> amp, double reference = 1.0)
    {
        double[] result = new double[amp.Length];
        Amp2Db(amp, result, reference);
        return result;
    }

    /// <summary>
    /// Converts an array of dB values to linear amplitude.
    /// Uses SIMD (Vector&lt;double&gt;) for the exp step.
    /// </summary>
    public static void Db2Amp(ReadOnlySpan<double> dB, Span<double> output, double reference = 1.0)
    {
        if (reference == 0.0) throw new ArgumentException("Reference must not be zero.", nameof(reference));
        if (output.Length < dB.Length) throw new ArgumentException("Output span is too short.", nameof(output));

        for (int i = 0; i < dB.Length; i++)
            output[i] = Math.Exp(dB[i] * Ln10Over20) * reference;
    }

    /// <summary>
    /// Converts an array of dB values to linear amplitude, allocating and returning a new array.
    /// </summary>
    public static double[] Db2Amp(ReadOnlySpan<double> dB, double reference = 1.0)
    {
        double[] result = new double[dB.Length];
        Db2Amp(dB, result, reference);
        return result;
    }

    // ----- Power (energy) variants -----

    /// <summary>10 * log10(power / ref) — for energy/power quantities.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double Power2Db(double power, double reference = 1.0)
    {
        if (reference == 0.0) throw new ArgumentException("Reference must not be zero.", nameof(reference));
        if (power == 0.0) power = 1e-24;
        return 10.0 * Math.Log10(power / reference);
    }

    private const double Ln10Over10 = 0.23025850929940457; // ln(10)/10

    /// <summary>ref * 10^(dB/10) — reverse of <see cref="Power2Db"/>.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double Db2Power(double dB, double reference = 1.0)
    {
        if (reference == 0.0) throw new ArgumentException("Reference must not be zero.", nameof(reference));
        return Math.Exp(dB * Ln10Over10) * reference;
    }
}
