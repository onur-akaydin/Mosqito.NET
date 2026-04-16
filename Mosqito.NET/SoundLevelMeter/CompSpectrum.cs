using System.Buffers;
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
    /// <summary>Window type for the spectral analysis.</summary>
    public enum WindowType { Hanning, Blackman, Rectangular }

    private const double Ref2e5 = 2e-5;
    private const double Ln10Over20 = 0.11512925464970229;

    /// <summary>
    /// Computes the one-sided spectrum of a time signal in Pa.
    /// </summary>
    public static (double[] spectrum, double[] freqAxis) Compute(
        ReadOnlySpan<double> signal, int fs,
        int nfft = -1,
        WindowType window = WindowType.Hanning,
        bool oneSided = true,
        bool db = true)
    {
        if (nfft < 0) nfft = signal.Length;

        // Get cached (read-only) raw window then compute normalisation scalar inline
        // so we never mutate the shared cached array.
        double[] win = window switch
        {
            WindowType.Hanning     => Windows.Hann(nfft),
            WindowType.Blackman    => Windows.Blackman(nfft),
            WindowType.Rectangular => Windows.Rectangular(nfft),
            _                      => Windows.Hann(nfft)
        };
        double winSum = Windows.Sum(win);
        double normFactor = winSum > 0.0 ? 1.0 / winSum : 1.0;

        // Pool the windowed buffer — saves one large allocation per call.
        double[] windowed = ArrayPool<double>.Shared.Rent(nfft);
        try
        {
            int copyLen = Math.Min(nfft, signal.Length);
            for (int i = 0; i < copyLen; i++) windowed[i] = signal[i] * win[i] * normFactor;
            for (int i = copyLen; i < nfft; i++) windowed[i] = 0.0;

            // Use Fft2(ReadOnlySpan<double>, Complex[]) with thread-local scratch
            Complex[] full = Fft.GetThreadBuf(nfft);
            Fft.Fft2(windowed.AsSpan(0, nfft), full);

            double[] spectrum;
            double[] freqAxis;

            if (oneSided)
            {
                int half = nfft / 2;
                spectrum = new double[half];
                freqAxis = new double[half];
                double fScale = (double)fs / nfft;
                const double factor = 1.42;

                if (db)
                {
                    double invRef = 1.0 / Ref2e5;
                    for (int i = 0; i < half; i++)
                    {
                        double amp = full[i].Magnitude * factor;
                        if (amp == 0.0) amp = 2e-12;
                        // 20 * log10(amp * invRef) = 20 * log10 / ln10 * ln(amp * invRef)
                        spectrum[i] = 20.0 * Math.Log10(amp * invRef);
                        freqAxis[i] = (i + 1) * fScale;
                    }
                }
                else
                {
                    for (int i = 0; i < half; i++)
                    {
                        spectrum[i] = full[i].Magnitude * factor;
                        freqAxis[i] = (i + 1) * fScale;
                    }
                }
            }
            else
            {
                int half = nfft / 2;
                spectrum = new double[nfft];
                freqAxis = new double[nfft];
                double fScale = (double)fs / nfft;
                const double factor = 1.42;

                if (db)
                {
                    double invRef = 1.0 / Ref2e5;
                    for (int i = 0; i < nfft; i++)
                    {
                        double amp = full[i].Magnitude * factor;
                        if (amp == 0.0) amp = 2e-12;
                        spectrum[i] = 20.0 * Math.Log10(amp * invRef);
                    }
                }
                else
                {
                    for (int i = 0; i < nfft; i++)
                        spectrum[i] = full[i].Magnitude * factor;
                }

                for (int i = 0; i < half; i++) freqAxis[i] = (i + 1) * fScale;
                for (int i = half; i < nfft; i++) freqAxis[i] = (nfft - i + 1) * fScale;
            }

            return (spectrum, freqAxis);
        }
        finally
        {
            ArrayPool<double>.Shared.Return(windowed);
        }
    }
}
