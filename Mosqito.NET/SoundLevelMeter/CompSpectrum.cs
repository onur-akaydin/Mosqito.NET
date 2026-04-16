using System.Numerics;
using Mosqito.Conversion;
using Mosqito.Dsp;

namespace Mosqito.SoundLevelMeter;

/// <summary>
/// Computes a one-sided FFT spectrum from a time signal.
/// Ported from MoSQITo <c>comp_spectrum.py</c>.
/// </summary>
public static class CompSpectrum
{
    /// <summary>
    /// Window type for the spectral analysis.
    /// </summary>
    public enum WindowType { Hanning, Blackman, Rectangular }

    /// <summary>
    /// Computes the one-sided spectrum of a time signal in Pa.
    /// </summary>
    /// <param name="signal">Time signal [Pa].</param>
    /// <param name="fs">Sampling frequency [Hz].</param>
    /// <param name="nfft">FFT length. -1 = use signal length (default).</param>
    /// <param name="window">Window type (default Hanning).</param>
    /// <param name="oneSided">Return one-sided spectrum (default true).</param>
    /// <param name="db">Return spectrum in dB (default true).</param>
    /// <returns>
    /// (<c>spectrum</c> — dB or amplitude spectrum;
    ///  <c>freqAxis</c> — frequency axis in Hz).
    /// </returns>
    public static (double[] spectrum, double[] freqAxis) Compute(
        ReadOnlySpan<double> signal, int fs,
        int nfft = -1,
        WindowType window = WindowType.Hanning,
        bool oneSided = true,
        bool db = true)
    {
        if (nfft < 0) nfft = signal.Length;

        // Build window
        double[] win = new double[nfft];
        switch (window)
        {
            case WindowType.Hanning:    Windows.FillHann(win); break;
            case WindowType.Blackman:   Windows.FillBlackman(win); break;
            case WindowType.Rectangular: Windows.FillRectangular(win); break;
        }

        // Amplitude correction: window / sum(window)
        Windows.NormaliseBySumInPlace(win);

        // Apply window and FFT
        double[] windowed = new double[nfft];
        int copyLen = Math.Min(nfft, signal.Length);
        for (int i = 0; i < copyLen; i++) windowed[i] = signal[i] * win[i];

        Complex[] full = Fft.Fft2(windowed);

        double[] spectrum;
        double[] freqAxis;

        if (oneSided)
        {
            int half = nfft / 2;
            spectrum  = new double[half];
            freqAxis  = new double[half];
            double fScale = (double)fs / nfft;
            for (int i = 0; i < half; i++)
            {
                spectrum[i] = full[i].Magnitude * 1.42;  // match Python factor
                freqAxis[i] = (i + 1) * fScale;
            }
        }
        else
        {
            int half = nfft / 2;
            spectrum = new double[nfft];
            freqAxis = new double[nfft];
            double fScale = (double)fs / nfft;
            for (int i = 0; i < nfft; i++) spectrum[i] = full[i].Magnitude * 1.42;
            for (int i = 0; i < half; i++) freqAxis[i] = (i + 1) * fScale;
            for (int i = half; i < nfft; i++) freqAxis[i] = (nfft - i + 1) * fScale;
        }

        if (db)
        {
            const double ref2e5 = 2e-5;
            for (int i = 0; i < spectrum.Length; i++)
                spectrum[i] = AmpDb.Amp2Db(spectrum[i], ref2e5);
        }

        return (spectrum, freqAxis);
    }
}
