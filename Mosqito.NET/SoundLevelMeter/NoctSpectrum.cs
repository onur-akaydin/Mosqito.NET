using System.Buffers;
using System.Collections.Concurrent;
using Mosqito.Dsp;

namespace Mosqito.SoundLevelMeter;

/// <summary>
/// Nth-octave band spectrum from a time-domain signal.
/// Ported from MoSQITo <c>noct_spectrum.py</c> and <c>_n_oct_time_filter.py</c>.
/// </summary>
public static class NoctSpectrum
{
    // -----------------------------------------------------------------------
    // Plan cache — SOS matrices and decimation factors keyed on call parameters.
    // For a fixed (fs, fmin, fmax, n, G, fr) configuration the bandpass SOS
    // and decimation Cheby SOS are constant and computed once.
    // -----------------------------------------------------------------------

    private sealed record NoctKey(double Fs, double Fmin, double Fmax, int N, int G, double Fr);

    private sealed class BandPlan
    {
        public readonly int      DecimQ;     // 0 = no decimation
        public readonly double[,] DecimSos;  // Cheby LP SOS for decimation (null if DecimQ==0)
        public readonly double[,] BandSos;   // Butterworth bandpass SOS
        public BandPlan(int decimQ, double[,] decimSos, double[,] bandSos)
        { DecimQ = decimQ; DecimSos = decimSos; BandSos = bandSos; }
    }

    private sealed class NoctPlan
    {
        public readonly BandPlan[] Bands;
        public readonly double[]   FPref;
        public NoctPlan(BandPlan[] bands, double[] fPref) { Bands = bands; FPref = fPref; }
    }

    private static readonly ConcurrentDictionary<NoctKey, NoctPlan> PlanCache = new();

    // -----------------------------------------------------------------------
    // Public entry point
    // -----------------------------------------------------------------------

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
        var key  = new NoctKey(fs, fmin, fmax, n, G, fr);
        var plan = PlanCache.GetOrAdd(key, k => BuildPlan(k.Fs, k.Fmin, k.Fmax, k.N, k.G, k.Fr));

        double[] spec = new double[plan.Bands.Length];

        // Rent one scratch buffer for the bandpass-filter output of non-decimated bands.
        // Reused across all bands — eliminates one full-signal alloc per band.
        double[] scratch = ArrayPool<double>.Shared.Rent(sig.Length);
        try
        {
            for (int i = 0; i < plan.Bands.Length; i++)
            {
                BandPlan band = plan.Bands[i];
                if (band.DecimQ > 0)
                {
                    // Decimation: FiltFilt still allocates its output, but the SOS is cached.
                    double[] decimated = ApplyDecimate(sig, band.DecimQ, band.DecimSos);
                    int dLen = decimated.Length;
                    double[] decimOut = ArrayPool<double>.Shared.Rent(dLen);
                    try
                    {
                        SosFilter.Process(band.BandSos, decimated.AsSpan(0, dLen), decimOut.AsSpan(0, dLen));
                        spec[i] = Rms(decimOut.AsSpan(0, dLen));
                    }
                    finally { ArrayPool<double>.Shared.Return(decimOut); }
                }
                else
                {
                    // Non-decimated bands: zero-copy path using the shared scratch buffer.
                    SosFilter.Process(band.BandSos, sig, scratch.AsSpan(0, sig.Length));
                    spec[i] = Rms(scratch.AsSpan(0, sig.Length));
                }
            }
        }
        finally { ArrayPool<double>.Shared.Return(scratch); }

        return (spec, plan.FPref);
    }

    // -----------------------------------------------------------------------
    // Plan builder — runs once per unique (fs, fmin, fmax, n, G, fr) tuple
    // -----------------------------------------------------------------------

    private static NoctPlan BuildPlan(double fs, double fmin, double fmax, int n, int G, double fr)
    {
        var (fcVec, fPref) = CenterFreq.Compute(fmin, fmax, n, G, fr);
        var (alphaVec, _, _) = FilterBandwidth.Compute(fcVec, n);

        var bands = new BandPlan[fcVec.Length];
        for (int i = 0; i < fcVec.Length; i++)
        {
            double fc    = fcVec[i];
            double alpha = alphaVec[i];

            if (fc > 0.88 * (fs / 2.0))
                throw new ArgumentException(
                    $"Filter centre frequency {fc} Hz exceeds 0.88 * Nyquist ({0.88 * fs / 2.0:F1} Hz).");

            int    decimQ    = 0;
            double effectiveFs = fs;
            double[,] decimSos = null!;

            if (fc < fs / 200.0)
            {
                int q = 2;
                while (fc < fs / q / 200.0) q++;
                decimQ     = q;
                effectiveFs = fs / q;
                decimSos   = Decimate.DesignCheby1Lp(8, 0.05, 0.8 / q);
            }

            double nyq = effectiveFs / 2.0;
            double w1  = Math.Max(1e-6, Math.Min(fc / nyq / alpha, 0.9999));
            double w2  = Math.Max(w1 + 1e-6, Math.Min(fc / nyq * alpha, 0.9999));
            double[,] bandSos = Butter.DesignBandpass(3, w1, w2);

            bands[i] = new BandPlan(decimQ, decimSos, bandSos);
        }

        return new NoctPlan(bands, fPref);
    }

    // -----------------------------------------------------------------------
    // Helpers
    // -----------------------------------------------------------------------

    // Applies cached Cheby LP decimation without re-designing the filter each call.
    private static double[] ApplyDecimate(ReadOnlySpan<double> input, int q, double[,] decimSos)
    {
        double[] filtered = SosFilter.FiltFilt(decimSos, input);
        int nOut = (input.Length + q - 1) / q;
        double[] output = new double[nOut];
        for (int i = 0; i < nOut; i++) output[i] = filtered[i * q];
        return output;
    }

    private static double Rms(ReadOnlySpan<double> signal)
    {
        double sumSq = 0.0;
        for (int i = 0; i < signal.Length; i++) sumSq += signal[i] * signal[i];
        return Math.Sqrt(sumSq / signal.Length);
    }
}
