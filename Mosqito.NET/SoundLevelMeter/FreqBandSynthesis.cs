namespace Mosqito.SoundLevelMeter;

/// <summary>
/// Converts a fine-band dB spectrum to per-band levels between given frequency limits.
/// Ported from MoSQITo <c>sound_level_meter/freq_band_synthesis.py</c>.
/// </summary>
public static class FreqBandSynthesis
{
    /// <summary>
    /// Sums spectrum power within each band defined by <paramref name="fmin"/>/<paramref name="fmax"/> pairs.
    /// </summary>
    /// <param name="spectrum">dB spectrum (one-sided), length N.</param>
    /// <param name="freqs">Corresponding frequency axis [Hz], length N.</param>
    /// <param name="fmin">Lower band edges [Hz], length nBands.</param>
    /// <param name="fmax">Upper band edges [Hz], length nBands.</param>
    /// <returns>
    /// (<c>bandSpectrum</c> — power-sum dB level per band;
    ///  <c>centerFreqs</c> — (fmin+fmax)/2 per band).
    /// </returns>
    public static (double[] bandSpectrum, double[] centerFreqs) Compute(
        ReadOnlySpan<double> spectrum, ReadOnlySpan<double> freqs,
        ReadOnlySpan<double> fmin, ReadOnlySpan<double> fmax)
    {
        int n = freqs.Length;
        int nBands = fmin.Length;

        double[] bandSpectrum = new double[nBands];
        double[] centerFreqs  = new double[nBands];

        for (int b = 0; b < nBands; b++)
        {
            double lo = fmin[b];
            double hi = fmax[b];
            centerFreqs[b] = (lo + hi) / 2.0;

            double sumPow = 0.0;
            for (int i = 0; i < n; i++)
            {
                if (freqs[i] >= lo && freqs[i] <= hi)
                    sumPow += Math.Pow(10.0, spectrum[i] / 10.0);
            }

            bandSpectrum[b] = sumPow > 0.0 ? 10.0 * Math.Log10(sumPow) : -999.0;
        }

        return (bandSpectrum, centerFreqs);
    }
}
