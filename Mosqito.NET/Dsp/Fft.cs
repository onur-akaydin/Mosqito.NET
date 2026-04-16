using System.Buffers;
using System.Numerics;
using MathNet.Numerics;
using MathNet.Numerics.IntegralTransforms;

namespace Mosqito.Dsp;

/// <summary>
/// FFT wrappers mirroring numpy/scipy conventions used throughout MoSQITo.
/// All transforms operate in-place or via caller-supplied output spans to minimise allocation.
/// </summary>
public static class Fft
{
    // ------------------------------------------------------------------
    // Forward real FFT  (matches numpy.fft.rfft)
    // Returns N/2+1 complex bins.
    // ------------------------------------------------------------------

    /// <summary>
    /// Computes the one-sided (real) FFT of <paramref name="input"/>.
    /// Output length is <c>input.Length / 2 + 1</c>.
    /// </summary>
    public static Complex[] Rfft(ReadOnlySpan<double> input)
    {
        int n = input.Length;
        // Pad / copy into a Complex array and forward-transform.
        Complex[] buf = new Complex[n];
        for (int i = 0; i < n; i++) buf[i] = new Complex(input[i], 0.0);
        Fourier.Forward(buf, FourierOptions.AsymmetricScaling);
        // Return one-sided half (0..N/2 inclusive).
        int half = n / 2 + 1;
        Complex[] result = new Complex[half];
        buf.AsSpan(0, half).CopyTo(result);
        return result;
    }

    /// <summary>
    /// Writes one-sided FFT output into <paramref name="output"/> (length must be ≥ N/2+1).
    /// </summary>
    public static void Rfft(ReadOnlySpan<double> input, Span<Complex> output)
    {
        int n = input.Length;
        int half = n / 2 + 1;
        if (output.Length < half) throw new ArgumentException("Output too short.", nameof(output));

        Complex[] buf = ArrayPool<Complex>.Shared.Rent(n);
        try
        {
            for (int i = 0; i < n; i++) buf[i] = new Complex(input[i], 0.0);
            // Fourier.Forward requires a contiguous array slice, not Span
            Complex[] slice = new Complex[n];
            Array.Copy(buf, slice, n);
            Fourier.Forward(slice, FourierOptions.AsymmetricScaling);
            for (int i = 0; i < half; i++) output[i] = slice[i];
        }
        finally
        {
            ArrayPool<Complex>.Shared.Return(buf);
        }
    }

    // ------------------------------------------------------------------
    // Full (two-sided) FFT  (matches numpy.fft.fft)
    // ------------------------------------------------------------------

    /// <summary>
    /// Computes the full complex FFT. Output length equals input length.
    /// </summary>
    public static Complex[] Fft2(ReadOnlySpan<double> input)
    {
        int n = input.Length;
        Complex[] buf = new Complex[n];
        for (int i = 0; i < n; i++) buf[i] = new Complex(input[i], 0.0);
        Fourier.Forward(buf, FourierOptions.AsymmetricScaling);
        return buf;
    }

    /// <summary>Writes full FFT result into <paramref name="output"/> (same length as input).</summary>
    public static void Fft2(ReadOnlySpan<double> input, Complex[] output)
    {
        int n = input.Length;
        if (output.Length < n) throw new ArgumentException("Output too short.", nameof(output));
        for (int i = 0; i < n; i++) output[i] = new Complex(input[i], 0.0);
        Fourier.Forward(output, FourierOptions.AsymmetricScaling);
    }

    // ------------------------------------------------------------------
    // Inverse real FFT  (matches numpy.fft.irfft)
    // Input is the one-sided spectrum of length N/2+1;
    // output is real of length N.
    // ------------------------------------------------------------------

    /// <summary>
    /// Inverse real FFT. <paramref name="spectrum"/> is the one-sided half (N/2+1 bins).
    /// <paramref name="nOutput"/> specifies the desired output length (default 2*(N/2+1-1)=N).
    /// </summary>
    public static double[] Irfft(ReadOnlySpan<Complex> spectrum, int nOutput = -1)
    {
        int half = spectrum.Length;
        if (nOutput < 0) nOutput = 2 * (half - 1);

        Complex[] buf = new Complex[nOutput];
        // Fill positive-frequency half.
        for (int i = 0; i < half && i < nOutput; i++) buf[i] = spectrum[i];
        // Mirror conjugate for negative frequencies.
        for (int i = half; i < nOutput; i++) buf[i] = Complex.Conjugate(buf[nOutput - i]);
        // Inverse FFT.
        Fourier.Inverse(buf, FourierOptions.AsymmetricScaling);
        double[] result = new double[nOutput];
        for (int i = 0; i < nOutput; i++) result[i] = buf[i].Real;
        return result;
    }

    // ------------------------------------------------------------------
    // Frequency axis helpers
    // ------------------------------------------------------------------

    /// <summary>
    /// Returns the one-sided frequency axis for an FFT of length <paramref name="n"/>
    /// at sample rate <paramref name="fs"/>.
    /// Length = n/2+1, values span [0 … fs/2].
    /// Matches numpy.fft.rfftfreq(n, 1/fs).
    /// </summary>
    public static double[] RfftFreq(int n, double fs)
    {
        int half = n / 2 + 1;
        double[] freq = new double[half];
        double scale = fs / n;
        for (int i = 0; i < half; i++) freq[i] = i * scale;
        return freq;
    }

    /// <summary>
    /// Returns the full (two-sided) frequency axis for an FFT of length <paramref name="n"/>
    /// at sample rate <paramref name="fs"/>. Matches numpy.fft.fftfreq(n, 1/fs).
    /// </summary>
    public static double[] FftFreq(int n, double fs)
    {
        double[] freq = new double[n];
        double scale = fs / n;
        int half = (n + 1) / 2;
        for (int i = 0; i < half; i++) freq[i] = i * scale;
        for (int i = half; i < n; i++) freq[i] = (i - n) * scale;
        return freq;
    }

    // ------------------------------------------------------------------
    // Helpers
    // ------------------------------------------------------------------

    /// <summary>Returns the smallest power of 2 that is ≥ <paramref name="n"/>.</summary>
    public static int NextPow2(int n)
    {
        if (n <= 1) return 1;
        return 1 << (int)Math.Ceiling(Math.Log2(n));
    }

    /// <summary>Computes the magnitude of each complex bin in-place.</summary>
    public static double[] Magnitude(ReadOnlySpan<Complex> spectrum)
    {
        double[] mag = new double[spectrum.Length];
        for (int i = 0; i < spectrum.Length; i++) mag[i] = spectrum[i].Magnitude;
        return mag;
    }
}
