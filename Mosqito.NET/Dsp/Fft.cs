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
    // Per-thread scratch Complex[] buffers keyed by length — allocated once per (thread, size).
    [ThreadStatic]
    private static Dictionary<int, Complex[]>? _threadBufs;

    /// <summary>
    /// Returns a thread-local <see cref="Complex"/> scratch buffer of exactly <paramref name="n"/>
    /// elements. The same array is reused on every call from the same thread with the same length.
    /// Callers must not hold a reference across an await or Parallel.For boundary.
    /// </summary>
    internal static Complex[] GetThreadBuf(int n)
    {
        _threadBufs ??= new Dictionary<int, Complex[]>();
        if (!_threadBufs.TryGetValue(n, out var buf))
        {
            buf = new Complex[n];
            _threadBufs[n] = buf;
        }
        return buf;
    }
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
    /// Uses a thread-local scratch buffer — do not call across await or Parallel.For boundaries
    /// without ensuring each worker thread uses its own call.
    /// </summary>
    public static void Rfft(ReadOnlySpan<double> input, Span<Complex> output)
    {
        int n = input.Length;
        int half = n / 2 + 1;
        if (output.Length < half) throw new ArgumentException("Output too short.", nameof(output));
        Complex[] buf = GetThreadBuf(n);
        for (int i = 0; i < n; i++) buf[i] = new Complex(input[i], 0.0);
        Fourier.Forward(buf, FourierOptions.AsymmetricScaling);
        for (int i = 0; i < half; i++) output[i] = buf[i];
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

        Complex[] buf = GetThreadBuf(nOutput);
        for (int i = 0; i < half && i < nOutput; i++) buf[i] = spectrum[i];
        for (int i = half; i < nOutput; i++) buf[i] = Complex.Conjugate(buf[nOutput - i]);
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
