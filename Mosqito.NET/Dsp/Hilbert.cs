using System.Numerics;
using MathNet.Numerics.IntegralTransforms;

namespace Mosqito.Dsp;

/// <summary>
/// Computes the analytic signal via the Hilbert transform.
/// Matches scipy.signal.hilbert(x).
/// </summary>
public static class Hilbert
{
    /// <summary>
    /// Returns the analytic signal: <c>x + j*H(x)</c>.
    /// The real part is the original signal; the imaginary part is the
    /// Hilbert transform.
    /// </summary>
    /// <param name="x">Real input signal of length N.</param>
    /// <returns>Complex analytic signal of length N.</returns>
    public static Complex[] Analytic(ReadOnlySpan<double> x)
    {
        int n = x.Length;
        Complex[] X = new Complex[n];
        for (int i = 0; i < n; i++) X[i] = new Complex(x[i], 0.0);

        // Forward FFT
        Fourier.Forward(X, FourierOptions.AsymmetricScaling);

        // Build h = [1, 2, 2, …, 2, 1, 0, …, 0] (symmetric)
        // Positive-frequency bins doubled, DC and Nyquist left as-is, neg-freq zeroed.
        int half = n / 2;
        bool isEven = (n % 2 == 0);

        X[0] *= 1.0;          // DC: unchanged
        for (int i = 1; i < half; i++) X[i] *= 2.0;
        if (isEven)
        {
            X[half] *= 1.0;   // Nyquist: unchanged
            for (int i = half + 1; i < n; i++) X[i] = Complex.Zero;
        }
        else
        {
            X[half] *= 2.0;   // not Nyquist in odd case
            for (int i = half + 1; i < n; i++) X[i] = Complex.Zero;
        }

        // Inverse FFT
        Fourier.Inverse(X, FourierOptions.AsymmetricScaling);
        return X;
    }

    /// <summary>
    /// Returns just the Hilbert transform (imaginary part of analytic signal).
    /// </summary>
    public static double[] Transform(ReadOnlySpan<double> x)
    {
        var analytic = Analytic(x);
        double[] h = new double[x.Length];
        for (int i = 0; i < h.Length; i++) h[i] = analytic[i].Imaginary;
        return h;
    }

    /// <summary>
    /// Computes the instantaneous amplitude (envelope) of the signal.
    /// </summary>
    public static double[] Envelope(ReadOnlySpan<double> x)
    {
        var analytic = Analytic(x);
        double[] env = new double[x.Length];
        for (int i = 0; i < env.Length; i++) env[i] = analytic[i].Magnitude;
        return env;
    }
}
