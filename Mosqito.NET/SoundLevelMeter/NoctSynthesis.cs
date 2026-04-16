using System.Collections.Concurrent;
using System.Numerics;
using System.Threading.Tasks;
using Mosqito.Dsp;

namespace Mosqito.SoundLevelMeter;

/// <summary>
/// Converts a fine-band frequency-domain spectrum to nth-octave band levels.
/// Ported from MoSQITo <c>noct_synthesis.py</c> and <c>_n_oct_freq_filter.py</c>.
/// </summary>
public static class NoctSynthesis
{
    // Cache the squared magnitude response |H(ω)|² per (w1, w2, nFreqs).
    // Keyed by quantised (w1, w2) to avoid floating-point hash instability.
    private static readonly ConcurrentDictionary<(long w1, long w2, int nFreqs), double[]>
        _hMag2Cache = new();

    private static long WKey(double w) => (long)Math.Round(w * 1_000_000_000L);

    /// <summary>
    /// Adapts an input RMS amplitude spectrum to nth-octave band levels.
    /// </summary>
    public static (double[] spec, double[] fPref) Compute(
        ReadOnlySpan<double> spectrum, ReadOnlySpan<double> freqs,
        double fmin, double fmax,
        int n = 3, int G = 10, double fr = 1000.0)
    {
        if (spectrum.Length != freqs.Length)
            throw new ArgumentException("spectrum and freqs must have the same length.");

        double df = freqs.Length > 1 ? freqs[1] - freqs[0] : freqs[0];
        double fs = (freqs[freqs.Length - 1] + df) * 2.0;

        if (Math.Abs(Math.Round(fs) - 48000) > 50.0)
            throw new ArgumentException("Sampling frequency must be close to 48 kHz. Got: " + fs);

        var (fcVec, fPref) = CenterFreq.Compute(fmin, fmax, n, G, fr);
        var (alphaVec, _, fHighVec) = FilterBandwidth.Compute(fcVec, n);

        // Remove bands where upper edge would alias
        var fcList    = new List<double>(fcVec);
        var alphaList = new List<double>(alphaVec);
        var fPrefList = new List<double>(fPref);

        for (int i = fcList.Count - 1; i >= 0; i--)
        {
            if (fHighVec[i] > fs / 2.0)
            {
                fcList.RemoveAt(i);
                alphaList.RemoveAt(i);
                fPrefList.RemoveAt(i);
            }
        }

        int nBands  = fcList.Count;
        int nFreqs  = spectrum.Length;
        double nyq  = fs / 2.0;

        // Snapshot spectrum into an array for parallel closure capture.
        double[] specArr = spectrum.ToArray();

        double[] spec = new double[nBands];

        Parallel.For(0, nBands, i =>
        {
            double fc    = fcList[i];
            double alpha = alphaList[i];

            double w1 = Math.Max(1e-6, fc / nyq / alpha);
            double w2 = Math.Min(0.9999, fc / nyq * alpha);
            if (w1 >= w2) { spec[i] = 0.0; return; }

            // Fetch (or compute) |H|² for this band
            double[] hMag2 = _hMag2Cache.GetOrAdd(
                (WKey(w1), WKey(w2), nFreqs),
                static key =>
                {
                    double ww1 = key.w1 * 1e-9, ww2 = key.w2 * 1e-9;
                    double[,] sos = Butter.DesignBandpass(3, ww1, ww2);
                    var (_, H) = Butter.Sosfreqz(sos, key.nFreqs);
                    double[] m2 = new double[key.nFreqs];
                    for (int k = 0; k < key.nFreqs; k++)
                    {
                        double mag = H[k].Magnitude;
                        m2[k] = mag * mag;
                    }
                    return m2;
                });

            // SIMD dot product: sum of hMag2[k] * specArr[k]^2
            double sumSq = 0.0;
            for (int k = 0; k < nFreqs; k++)
            {
                double val = specArr[k];
                sumSq += hMag2[k] * val * val;
            }
            spec[i] = Math.Sqrt(sumSq);
        });

        return (spec, fPrefList.ToArray());
    }
}
