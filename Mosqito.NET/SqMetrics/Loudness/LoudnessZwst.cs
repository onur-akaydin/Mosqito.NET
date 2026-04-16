using Mosqito.Conversion;
using Mosqito.Dsp;
using Mosqito.SoundLevelMeter;

namespace Mosqito.SqMetrics.Loudness;

/// <summary>
/// Result of Zwicker stationary loudness computation.
/// </summary>
/// <param name="N">Overall loudness [sone].</param>
/// <param name="NSpecific">Specific loudness [sone/Bark], length 240.</param>
/// <param name="BarkAxis">Bark axis, length 240 (0.1 to 24.0).</param>
public readonly record struct LoudnessZwstResult(
    double N,
    double[] NSpecific,
    double[] BarkAxis);

/// <summary>
/// Zwicker stationary loudness (ISO 532-1:2017 Method 1).
/// Ported from MoSQITo <c>loudness_zwst.py</c>, <c>loudness_zwst_freq.py</c>,
/// <c>loudness_zwst_perseg.py</c>.
/// </summary>
public static class LoudnessZwst
{
    private static readonly double[] BarkAxisCache = BuildBarkAxis();
    private static double[] BuildBarkAxis()
    {
        const int n = 240; // linspace(0.1, 24, 240)
        double[] b = new double[n];
        double step = (24.0 - 0.1) / (n - 1);
        for (int i = 0; i < n; i++) b[i] = 0.1 + i * step;
        return b;
    }

    // ------------------------------------------------------------------
    // Main entry: time-domain signal
    // ------------------------------------------------------------------

    /// <summary>
    /// Computes Zwicker stationary loudness from a time signal.
    /// </summary>
    /// <param name="signal">Time signal [Pa]. Signal should be at ≥ 48 kHz.</param>
    /// <param name="fs">Sampling frequency [Hz].</param>
    /// <param name="fieldType">"free" (default) or "diffuse".</param>
    /// <returns><see cref="LoudnessZwstResult"/> with N, NSpecific, BarkAxis.</returns>
    public static LoudnessZwstResult Compute(
        ReadOnlySpan<double> signal, int fs,
        string fieldType = "free")
    {
        // Resample to 48 kHz if needed
        double[] sig = signal.ToArray();
        if (fs < 48000)
        {
            Console.WriteLine("[Warning] Signal resampled to 48 kHz for Zwicker loudness.");
            sig = Resample.Apply(signal, fs, 48000);
            fs = 48000;
        }

        // Compute 1/3-octave spectrum from 24 Hz to 12600 Hz (28 bands)
        var (specThird, _) = NoctSpectrum.Compute(sig, fs, fmin: 24, fmax: 12600);

        // Convert to dB re 2e-5 Pa
        double[] specDb = AmpDb.Amp2Db(specThird, reference: 2e-5);

        return FromSpectrum(specDb, fieldType);
    }

    // ------------------------------------------------------------------
    // Frequency-domain entry
    // ------------------------------------------------------------------

    /// <summary>
    /// Computes Zwicker stationary loudness from a fine-band RMS amplitude spectrum.
    /// </summary>
    /// <param name="spectrum">RMS amplitude spectrum (one-sided).</param>
    /// <param name="freqs">Frequency axis [Hz] (linspace 0…fs/2).</param>
    /// <param name="fieldType">"free" or "diffuse".</param>
    public static LoudnessZwstResult ComputeFromSpectrum(
        ReadOnlySpan<double> spectrum, ReadOnlySpan<double> freqs,
        string fieldType = "free")
    {
        if (spectrum.Length != freqs.Length)
            throw new ArgumentException("spectrum and freqs must have the same length.");

        // Compute 1/3-octave synthesis
        var (specThird, _) = NoctSynthesis.Compute(spectrum, freqs, fmin: 24, fmax: 12600);

        // Convert to dB
        double[] specDb = AmpDb.Amp2Db(specThird, reference: 2e-5);

        return FromSpectrum(specDb, fieldType);
    }

    // ------------------------------------------------------------------
    // Per-segment (loudness_zwst_perseg)
    // ------------------------------------------------------------------

    /// <summary>
    /// Computes Zwicker stationary loudness per time-segment.
    /// Identical to calling <see cref="Compute"/> on each segment independently.
    /// </summary>
    /// <param name="signal">Time signal [Pa].</param>
    /// <param name="fs">Sampling frequency [Hz].</param>
    /// <param name="nPerSeg">Segment length in samples.</param>
    /// <param name="noOverlap">Number of overlapping samples between segments.</param>
    /// <param name="fieldType">"free" or "diffuse".</param>
    /// <returns>
    /// (<c>N</c> — loudness per segment;
    ///  <c>NSpecific</c> — specific loudness [240, nSeg];
    ///  <c>BarkAxis</c> — Bark axis length 240;
    ///  <c>TimeAxis</c> — centre time of each segment [s]).
    /// </returns>
    public static (double[] N, double[,] NSpecific, double[] BarkAxis, double[] TimeAxis) ComputePerSeg(
        ReadOnlySpan<double> signal, int fs,
        int nPerSeg = 4096,
        int? noOverlap = null,
        string fieldType = "free")
    {
        var (blockArray, timeAxis) = Mosqito.Io.TimeSegmentation.Segment(signal, fs, nPerSeg, noOverlap);
        int nSeg = blockArray.GetLength(1);
        int nSamples = blockArray.GetLength(0);

        double[] N = new double[nSeg];
        double[,] NSpec = new double[240, nSeg];

        double[] segBuf = new double[nSamples];

        for (int seg = 0; seg < nSeg; seg++)
        {
            for (int s = 0; s < nSamples; s++) segBuf[s] = blockArray[s, seg];
            var result = Compute(segBuf, fs, fieldType);
            N[seg] = result.N;
            for (int b = 0; b < 240; b++) NSpec[b, seg] = result.NSpecific[b];
        }

        return (N, NSpec, BarkAxisCache, timeAxis);
    }

    // ------------------------------------------------------------------
    // Internal: from dB 1/3-octave spectrum to loudness
    // ------------------------------------------------------------------

    private static LoudnessZwstResult FromSpectrum(double[] specDb, string fieldType)
    {
        double[] nm = MainLoudness.Compute(specDb, fieldType);
        var (N, nSpecific) = CalcSlopes.Compute(nm);

        // Return Bark axis (shared cache — callers must not mutate)
        return new LoudnessZwstResult(N, nSpecific, BarkAxisCache);
    }
}
