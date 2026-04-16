using Mosqito.Dsp;

namespace Mosqito.SoundLevelMeter;

/// <summary>
/// Nth-octave band spectrum from a time-domain signal.
/// Ported from MoSQITo <c>noct_spectrum.py</c> and <c>_n_oct_time_filter.py</c>.
/// </summary>
public static class NoctSpectrum
{
    /// <summary>
    /// Computes the RMS level in each nth-octave band for a time-domain signal.
    /// </summary>
    /// <param name="sig">Time signal [Pa]. 1-D array.</param>
    /// <param name="fs">Sampling frequency [Hz].</param>
    /// <param name="fmin">Minimum centre frequency [Hz].</param>
    /// <param name="fmax">Maximum centre frequency [Hz].</param>
    /// <param name="n">Bands per octave (default 3).</param>
    /// <param name="G">Base system: 2 or 10 (default 10).</param>
    /// <param name="fr">Reference frequency [Hz] (default 1000).</param>
    /// <returns>
    /// (<c>spec</c> — RMS amplitude per band, length nBands;
    ///  <c>fPref</c> — preferred centre frequencies, same length).
    /// </returns>
    public static (double[] spec, double[] fPref) Compute(
        ReadOnlySpan<double> sig, double fs,
        double fmin, double fmax,
        int n = 3, int G = 10, double fr = 1000.0)
    {
        var (fcVec, fPref) = CenterFreq.Compute(fmin, fmax, n, G, fr);
        var (alphaVec, _, _) = FilterBandwidth.Compute(fcVec, n);

        double[] spec = new double[fcVec.Length];
        for (int i = 0; i < fcVec.Length; i++)
            spec[i] = TimeFilter(sig, fs, fcVec[i], alphaVec[i]);

        return (spec, fPref);
    }

    // ------------------------------------------------------------------
    // Per-band bandpass filter and RMS  (_n_oct_time_filter)
    // ------------------------------------------------------------------

    /// <summary>
    /// Applies a Butterworth bandpass filter and returns the RMS level.
    /// </summary>
    private static double TimeFilter(ReadOnlySpan<double> sig, double fs, double fc, double alpha,
        int filterOrder = 3)
    {
        double nyq = fs / 2.0;

        if (fc > 0.88 * nyq)
            throw new ArgumentException(
                $"Filter centre frequency {fc} Hz exceeds 0.88 * Nyquist ({0.88 * nyq:F1} Hz).");

        // Downsample when fc is very low relative to fs to avoid numerical issues
        double[] sigArr = sig.ToArray();
        double effectiveFs = fs;

        if (fc < fs / 200.0)
        {
            int q = 2;
            while (fc < effectiveFs / q / 200.0) q++;
            sigArr = Decimate.Apply(sigArr, q);
            effectiveFs = fs / q;
        }

        // Normalised cutoff frequencies
        double w1 = fc / (effectiveFs / 2.0) / alpha;
        double w2 = fc / (effectiveFs / 2.0) * alpha;

        // Clamp to valid range
        w1 = Math.Max(1e-6, Math.Min(w1, 0.9999));
        w2 = Math.Max(w1 + 1e-6, Math.Min(w2, 0.9999));

        double[,] sos = Butter.DesignBandpass(filterOrder, w1, w2);
        double[] filtered = SosFilter.Process(sos, sigArr);

        // RMS
        double sumSq = 0.0;
        for (int i = 0; i < filtered.Length; i++) sumSq += filtered[i] * filtered[i];
        return Math.Sqrt(sumSq / filtered.Length);
    }
}
