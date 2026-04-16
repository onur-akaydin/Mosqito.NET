using System.Numerics;
using Mosqito.Dsp;

namespace Mosqito.SoundLevelMeter;

/// <summary>
/// Converts a fine-band frequency-domain spectrum to nth-octave band levels.
/// Ported from MoSQITo <c>noct_synthesis.py</c> and <c>_n_oct_freq_filter.py</c>.
/// </summary>
public static class NoctSynthesis
{
    /// <summary>
    /// Adapts an input RMS amplitude spectrum to nth-octave band levels.
    /// </summary>
    /// <param name="spectrum">One-sided RMS amplitude spectrum, length nFreqs.</param>
    /// <param name="freqs">Frequency axis [Hz], length nFreqs (linspace 0…fs/2).</param>
    /// <param name="fmin">Minimum band centre frequency [Hz].</param>
    /// <param name="fmax">Maximum band centre frequency [Hz].</param>
    /// <param name="n">Bands per octave (default 3).</param>
    /// <param name="G">Base: 2 or 10 (default 10).</param>
    /// <param name="fr">Reference frequency [Hz] (default 1000).</param>
    /// <returns>
    /// (<c>spec</c> — RMS amplitude per band, length nBands;
    ///  <c>fPref</c> — preferred centre frequencies).
    /// </returns>
    public static (double[] spec, double[] fPref) Compute(
        ReadOnlySpan<double> spectrum, ReadOnlySpan<double> freqs,
        double fmin, double fmax,
        int n = 3, int G = 10, double fr = 1000.0)
    {
        if (spectrum.Length != freqs.Length)
            throw new ArgumentException("spectrum and freqs must have the same length.");

        // FftFreq gives bins 0…(N/2-1)*df; the Nyquist is N/2*df = last + df.
        double df = freqs.Length > 1 ? freqs[1] - freqs[0] : freqs[0];
        double fs = (freqs[freqs.Length - 1] + df) * 2.0;

        if (Math.Abs(Math.Round(fs) - 48000) > 50.0)
            throw new ArgumentException("Sampling frequency must be close to 48 kHz. Got: " + fs);

        var (fcVec, fPref) = CenterFreq.Compute(fmin, fmax, n, G, fr);
        var (alphaVec, _, fHighVec) = FilterBandwidth.Compute(fcVec, n);

        // Remove bands where upper edge would alias
        var fcList     = new List<double>(fcVec);
        var alphaList  = new List<double>(alphaVec);
        var fPrefList  = new List<double>(fPref);

        for (int i = fcList.Count - 1; i >= 0; i--)
        {
            if (fHighVec[i] > fs / 2.0)
            {
                fcList.RemoveAt(i);
                alphaList.RemoveAt(i);
                fPrefList.RemoveAt(i);
            }
        }

        double[] spec = new double[fcList.Count];
        for (int i = 0; i < fcList.Count; i++)
            spec[i] = FreqFilter(spectrum, freqs, fs, fcList[i], alphaList[i]);

        return (spec, fPrefList.ToArray());
    }

    // ------------------------------------------------------------------
    // Frequency-domain bandpass weighting  (_n_oct_freq_filter)
    // ------------------------------------------------------------------

    private static double FreqFilter(
        ReadOnlySpan<double> spectrum, ReadOnlySpan<double> freqs,
        double fs, double fc, double alpha, int filterOrder = 3)
    {
        int nFreqs = spectrum.Length;
        double nyq = fs / 2.0;

        // Normalised cutoffs
        double w1 = Math.Max(1e-6, fc / nyq / alpha);
        double w2 = Math.Min(0.9999, fc / nyq * alpha);
        if (w1 >= w2) return 0.0;

        double[,] sos = Butter.DesignBandpass(filterOrder, w1, w2);

        // Frequency response at nFreqs points (matching sosfreqz(sos, worN=nFreqs))
        var (_, H) = Butter.Sosfreqz(sos, nFreqs);

        // Apply filter frequency response
        double sumSq = 0.0;
        for (int i = 0; i < nFreqs; i++)
        {
            double hMag = H[i].Magnitude;
            double val = hMag * spectrum[i];
            sumSq += val * val;
        }

        return Math.Sqrt(sumSq);
    }
}
