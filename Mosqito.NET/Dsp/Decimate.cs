using System.Buffers;

namespace Mosqito.Dsp;

/// <summary>
/// Downsampling with anti-aliasing filter.
/// Matches scipy.signal.decimate(x, q, n=8, ftype='iir') — Chebyshev type I IIR
/// anti-aliasing at 0.8 * Nyquist / q (the scipy default).
/// </summary>
public static class Decimate
{
    /// <summary>
    /// Decimates <paramref name="input"/> by factor <paramref name="q"/>.
    /// Uses an 8th-order IIR Chebyshev-I anti-aliasing filter (scipy default)
    /// applied with zero-phase (forward-backward) filtering, matching
    /// scipy.signal.decimate behaviour.
    /// </summary>
    /// <param name="input">Input signal.</param>
    /// <param name="q">Decimation factor (integer ≥ 2).</param>
    /// <returns>Decimated output of length ceil(n / q).</returns>
    public static double[] Apply(ReadOnlySpan<double> input, int q)
    {
        if (q < 2) throw new ArgumentOutOfRangeException(nameof(q), "Decimation factor must be ≥ 2.");

        // Design Chebyshev-I LP anti-aliasing filter
        double wCutoff = 0.8 / q;   // normalised cutoff (scipy: Wn = 0.8/q relative to Nyquist)
        double[,] sos = ChebType1Sos(order: 8, wCutoff: wCutoff);

        // Forward-backward filter (zero phase)
        double[] filtered = SosFilter.FiltFilt(sos, input);

        // Downsample
        int nOut = (input.Length + q - 1) / q;
        double[] output = new double[nOut];
        for (int i = 0; i < nOut; i++)
            output[i] = filtered[i * q];
        return output;
    }

    /// <summary>Applies decimation along axis 0 of a 2-D array (each column independently).</summary>
    public static double[,] Apply2D(double[,] input, int q)
    {
        int nSamples = input.GetLength(0);
        int nCols    = input.GetLength(1);
        int nOut = (nSamples + q - 1) / q;

        double[,] output = new double[nOut, nCols];
        double[] col    = new double[nSamples];
        double[] colOut;

        for (int c = 0; c < nCols; c++)
        {
            for (int r = 0; r < nSamples; r++) col[r] = input[r, c];
            colOut = Apply(col, q);
            for (int r = 0; r < nOut; r++) output[r, c] = colOut[r];
        }
        return output;
    }

    // ------------------------------------------------------------------
    // Chebyshev Type-I LP filter design (8th order)
    // This is a simplified design matching the scipy defaults.
    // We implement via the standard analog prototype → bilinear approach.
    // ------------------------------------------------------------------
    private static double[,] ChebType1Sos(int order, double wCutoff)
    {
        // For scipy compat: ripple = 0.05 dB
        // Scipy actually uses 0.05 dB ripple in decimate().
        // We use Butterworth as an approximation for simplicity, unless wCutoff is within
        // valid range. For production accuracy against scipy we use the exact Chebyshev.
        // Here we implement a simple recursive Chebyshev-I LP prototype.

        double rp = 0.05; // passband ripple in dB (scipy default for decimate)
        return DesignCheby1Lp(order, rp, wCutoff);
    }

    /// <summary>
    /// Designs a Chebyshev Type-I low-pass digital filter as SOS sections.
    /// Matches scipy.signal.cheby1(N, rp, Wn, output='sos').
    /// </summary>
    public static double[,] DesignCheby1Lp(int N, double rpDb, double wn)
    {
        // Analog prototype poles for Chebyshev-I
        double eps = Math.Sqrt(Math.Pow(10.0, rpDb / 10.0) - 1.0);
        double mu = Math.Asinh(1.0 / eps) / N;
        double sinhMu = Math.Sinh(mu);
        double coshMu = Math.Cosh(mu);

        System.Numerics.Complex[] poles = new System.Numerics.Complex[N];
        for (int k = 1; k <= N; k++)
        {
            double theta = Math.PI * (2.0 * k - 1) / (2.0 * N);
            double re = -sinhMu * Math.Sin(theta);
            double im =  coshMu * Math.Cos(theta);
            poles[k - 1] = new System.Numerics.Complex(re, im);
        }

        // Gain of analog prototype
        double gainAnal;
        {
            System.Numerics.Complex prod = System.Numerics.Complex.One;
            foreach (var p in poles) prod *= -p;
            gainAnal = (N % 2 == 0 ? 1.0 / Math.Sqrt(1 + eps * eps) : 1.0) * prod.Real;
        }

        // Pre-warp cutoff
        double warped = 2.0 * Math.Tan(Math.PI * wn / 2.0);

        // Scale LP prototype to cutoff (LP→LP)
        System.Numerics.Complex[] polesW = new System.Numerics.Complex[N];
        for (int i = 0; i < N; i++) polesW[i] = poles[i] * warped;
        double gainW = gainAnal * Math.Pow(warped, N);

        // Bilinear transform
        System.Numerics.Complex[] zerosD = new System.Numerics.Complex[N];
        System.Numerics.Complex[] polesD = new System.Numerics.Complex[N];
        for (int i = 0; i < N; i++) zerosD[i] = new System.Numerics.Complex(-1, 0);
        double kD = gainW;
        System.Numerics.Complex kNum = System.Numerics.Complex.One,
                                kDen = System.Numerics.Complex.One;
        for (int i = 0; i < N; i++)
        {
            kNum *= 2.0 - polesW[i];
            polesD[i] = (2.0 + polesW[i]) / (2.0 - polesW[i]);
        }
        kD = gainW / kNum.Real;

        // Use Butter helper ZPK→SOS via the Butter.Sosfreqz pathway
        // Actually, build SOS directly from conjugate-pair grouping.
        return BuildSosFromZpk(zerosD, polesD, kD);
    }

    private static double[,] BuildSosFromZpk(
        System.Numerics.Complex[] z, System.Numerics.Complex[] p, double k)
    {
        // Group into sections using the same logic as Butter.Zpk2Sos
        // We call the internal logic via Butter public methods by hijacking DesignLowpass
        // at a near-0 frequency to get the SOS structure, then replace coefficients.
        // Alternatively: directly construct SOS from conjugate pairs.

        int n = p.Length;
        int nSec = (n + 1) / 2;
        double[,] sos = new double[nSec, 6];

        // Sort poles by imaginary part magnitude (complex first)
        var pList = new List<System.Numerics.Complex>(p);
        var zList = new List<System.Numerics.Complex>(z);
        pList.Sort((a, b) => Math.Abs(b.Imaginary).CompareTo(Math.Abs(a.Imaginary)));
        zList.Sort((a, b) => Math.Abs(b.Imaginary).CompareTo(Math.Abs(a.Imaginary)));
        while (zList.Count < nSec * 2) zList.Add(System.Numerics.Complex.Zero);

        for (int s = 0; s < nSec; s++)
        {
            System.Numerics.Complex p1 = pList[s * 2 < pList.Count ? s * 2 : 0];
            System.Numerics.Complex p2 = s * 2 + 1 < pList.Count ? pList[s * 2 + 1] :
                                         System.Numerics.Complex.Zero;
            System.Numerics.Complex z1 = zList[s * 2];
            System.Numerics.Complex z2 = s * 2 + 1 < zList.Count ? zList[s * 2 + 1] :
                                         System.Numerics.Complex.Zero;

            double gFactor = s == 0 ? k : 1.0;
            sos[s, 0] =  gFactor;
            sos[s, 1] = -gFactor * (z1.Real + z2.Real);
            sos[s, 2] =  gFactor * (z1.Real * z2.Real - z1.Imaginary * z2.Imaginary);
            sos[s, 3] = 1.0;
            sos[s, 4] = -(p1.Real + p2.Real);
            sos[s, 5] = p1.Real * p2.Real - p1.Imaginary * p2.Imaginary;
        }

        return sos;
    }
}
