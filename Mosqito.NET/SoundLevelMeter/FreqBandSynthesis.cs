using System.Buffers;

namespace Mosqito.SoundLevelMeter;

/// <summary>
/// Converts a fine-band dB spectrum to per-band levels between given frequency limits.
/// Ported from MoSQITo <c>sound_level_meter/freq_band_synthesis.py</c>.
/// </summary>
public static class FreqBandSynthesis
{
    private const double Ln10Over10 = 0.23025850929940457; // ln(10)/10

    /// <summary>
    /// Sums spectrum power within each band defined by <paramref name="fmin"/>/<paramref name="fmax"/> pairs.
    /// Assumes <paramref name="freqs"/> is monotonically non-decreasing (FFT frequency axis).
    /// </summary>
    public static (double[] bandSpectrum, double[] centerFreqs) Compute(
        ReadOnlySpan<double> spectrum, ReadOnlySpan<double> freqs,
        ReadOnlySpan<double> fmin, ReadOnlySpan<double> fmax)
    {
        int n      = freqs.Length;
        int nBands = fmin.Length;

        double[] bandSpectrum = new double[nBands];
        double[] centerFreqs  = new double[nBands];

        // Precompute linear power for every bin once — avoids O(nBands × n) transcendentals.
        double[] pow = ArrayPool<double>.Shared.Rent(n);
        try
        {
            for (int i = 0; i < n; i++)
                pow[i] = Math.Exp(spectrum[i] * Ln10Over10); // 10^(spec/10) = exp(spec*ln10/10)

            // Two-pointer sweep — O(n + nBands) instead of O(n × nBands).
            // Assumes fmin[] and fmax[] are sorted ascending and freqs[] is sorted ascending.
            // If either is not sorted, falls back to correct (but slower) linear scan.
            bool sorted = true;
            for (int b = 1; b < nBands && sorted; b++)
                if (fmin[b] < fmin[b - 1]) sorted = false;

            if (sorted)
            {
                int lo = 0;
                for (int b = 0; b < nBands; b++)
                {
                    double bandLo = fmin[b];
                    double bandHi = fmax[b];
                    centerFreqs[b] = (bandLo + bandHi) * 0.5;

                    // Advance lo to first freq >= bandLo
                    while (lo < n && freqs[lo] < bandLo) lo++;

                    double sumPow = 0.0;
                    for (int i = lo; i < n && freqs[i] <= bandHi; i++)
                        sumPow += pow[i];

                    bandSpectrum[b] = sumPow > 0.0 ? 10.0 * Math.Log10(sumPow) : -999.0;
                }
            }
            else
            {
                for (int b = 0; b < nBands; b++)
                {
                    double lo2 = fmin[b], hi2 = fmax[b];
                    centerFreqs[b] = (lo2 + hi2) * 0.5;
                    double sumPow = 0.0;
                    for (int i = 0; i < n; i++)
                        if (freqs[i] >= lo2 && freqs[i] <= hi2) sumPow += pow[i];
                    bandSpectrum[b] = sumPow > 0.0 ? 10.0 * Math.Log10(sumPow) : -999.0;
                }
            }
        }
        finally
        {
            ArrayPool<double>.Shared.Return(pow);
        }

        return (bandSpectrum, centerFreqs);
    }
}
