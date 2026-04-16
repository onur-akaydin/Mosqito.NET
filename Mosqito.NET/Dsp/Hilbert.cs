using System.Numerics;
using MathNet.Numerics.IntegralTransforms;

namespace Mosqito.Dsp;

/// <summary>
/// Computes the analytic signal via the Hilbert transform.
/// Matches scipy.signal.hilbert(x).
/// Uses a thread-local <see cref="Complex"/> scratch buffer to avoid per-call allocation.
/// </summary>
public static class Hilbert
{
    /// <summary>
    /// Returns the analytic signal: <c>x + j*H(x)</c>.
    /// The real part is the original signal; the imaginary part is the Hilbert transform.
    /// </summary>
    public static Complex[] Analytic(ReadOnlySpan<double> x)
    {
        int n = x.Length;
        Complex[] X = new Complex[n];
        AnalyticCore(x, X);
        return X;
    }

    /// <summary>
    /// Computes the instantaneous amplitude (envelope) of the signal.
    /// </summary>
    public static double[] Envelope(ReadOnlySpan<double> x)
    {
        int n = x.Length;
        Complex[] scratch = Fft.GetThreadBuf(n);
        AnalyticCore(x, scratch);
        double[] env = new double[n];
        for (int i = 0; i < n; i++) env[i] = scratch[i].Magnitude;
        return env;
    }

    /// <summary>
    /// Computes the envelope and writes it into a caller-supplied span — zero heap allocation
    /// after the first call on each thread for a given signal length.
    /// </summary>
    public static void Envelope(ReadOnlySpan<double> x, Span<double> envOut)
    {
        int n = x.Length;
        if (envOut.Length < n) throw new ArgumentException("envOut too short.", nameof(envOut));
        Complex[] scratch = Fft.GetThreadBuf(n);
        AnalyticCore(x, scratch);
        for (int i = 0; i < n; i++) envOut[i] = scratch[i].Magnitude;
    }

    /// <summary>
    /// Returns just the Hilbert transform (imaginary part of analytic signal).
    /// </summary>
    public static double[] Transform(ReadOnlySpan<double> x)
    {
        int n = x.Length;
        Complex[] scratch = Fft.GetThreadBuf(n);
        AnalyticCore(x, scratch);
        double[] h = new double[n];
        for (int i = 0; i < n; i++) h[i] = scratch[i].Imaginary;
        return h;
    }

    // ------------------------------------------------------------------
    // Core analytic-signal computation in a caller-supplied buffer.
    // ------------------------------------------------------------------
    private static void AnalyticCore(ReadOnlySpan<double> x, Complex[] X)
    {
        int n = x.Length;
        for (int i = 0; i < n; i++) X[i] = new Complex(x[i], 0.0);

        Fourier.Forward(X, FourierOptions.AsymmetricScaling);

        int half = n / 2;
        bool isEven = (n % 2 == 0);
        // DC unchanged
        for (int i = 1; i < half; i++) X[i] *= 2.0;
        if (isEven)
        {
            // Nyquist unchanged
            for (int i = half + 1; i < n; i++) X[i] = Complex.Zero;
        }
        else
        {
            X[half] *= 2.0;
            for (int i = half + 1; i < n; i++) X[i] = Complex.Zero;
        }

        Fourier.Inverse(X, FourierOptions.AsymmetricScaling);
    }
}
