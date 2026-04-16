using System.Buffers;
using System.Numerics;
using Mosqito.Dsp;

namespace Mosqito.SqMetrics.Loudness;

/// <summary>
/// Result of ECMA-418-2 loudness computation.
/// </summary>
public readonly record struct LoudnessEcmaResult(
    double N,               // Overall single-value loudness [sone_HMS]
    double[] NTime,         // Time-varying overall loudness [sone_HMS], length nTime
    double[][] NSpecific,   // Specific loudness [sone_HMS/bark], 53 bands × nTime each
    double[] BarkAxis,      // Bark axis (0.5 to 26.5, 53 values)
    double[][] TimeArray);  // Time axis per band (53 bands, nTime each)

/// <summary>
/// ECMA-418-2 (2nd Ed. 2022) specific and total loudness.
/// Ported from MoSQITo <c>loudness_ecma.py</c> and helper modules.
///
/// Pipeline (Section 5):
///   5.1.2 Windowing + zero-padding.
///   5.1.3 Outer and middle ear filter (8-section SOS).
///   5.1.4 Gammatone filter bank (53 complex 5th-order IIR filters).
///   5.1.5 Segmentation into blocks (sb/sh overlap-add).
///   5.1.6 Half-wave rectification (clip to 0).
///   5.1.7–5.1.9 RMS → compressive nonlinearity → subtract threshold.
///   Eq. 116  Time-dependent total loudness N(t) = sum_z N'(z,t) × Δz.
///   Eq. 117  Single-value N = (mean(N(t)^e))^(1/e), e = 1/log10(2).
/// </summary>
public static class LoudnessEcma
{
    // ------------------------------------------------------------------
    // Constants
    // ------------------------------------------------------------------
    private const int NBands    = 53;
    private const int Fs        = 48000;
    private const double DeltaZ = 0.5;  // Bark step
    private const int FilterOrder = 5;  // Gammatone filter order (k)

    // ------------------------------------------------------------------
    // ECMA-418-2 Table 1 — ear filter SOS (8 sections × 6 coefficients)
    // Format: [b0, b1, b2, a0, a1, a2]
    // ------------------------------------------------------------------
    internal static readonly double[,] EarSos =
    {
        { 1.015896, -1.925299,  0.922118,  1.0, -1.925299,  0.938014 },
        { 0.958943, -1.806088,  0.876439,  1.0, -1.806088,  0.835382 },
        { 0.961372, -1.763632,  0.821788,  1.0, -1.763632,  0.783160 },
        { 2.225804, -1.434650, -0.498204,  1.0, -1.434650,  0.727599 },
        { 0.471735, -0.366092,  0.244145,  1.0, -0.366092, -0.284120 },
        { 0.115267,  0.000000, -0.115267,  1.0, -1.796003,  0.805838 },
        { 0.988029, -1.912434,  0.926132,  1.0, -1.912434,  0.914161 },
        { 1.952238,  0.162320, -0.667994,  1.0,  0.162320,  0.284244 }
    };

    // ------------------------------------------------------------------
    // ECMA-418-2 ltq_z — threshold in quiet per band (53 values)
    // ------------------------------------------------------------------
    internal static readonly double[] LtqZ =
    {
        0.3310, 0.1625, 0.1051, 0.0757, 0.0576, 0.0453, 0.0365, 0.0298,
        0.0247, 0.0207, 0.0176, 0.0151, 0.0131, 0.0115, 0.0103, 0.0093,
        0.0086, 0.0081, 0.0077, 0.0074, 0.0073, 0.0072, 0.0071, 0.0072,
        0.0073, 0.0074, 0.0076, 0.0079, 0.0082, 0.0086, 0.0092, 0.0100,
        0.0109, 0.0122, 0.0138, 0.0157, 0.0172, 0.0180, 0.0180, 0.0177,
        0.0176, 0.0177, 0.0182, 0.0190, 0.0202, 0.0217, 0.0237, 0.0263,
        0.0296, 0.0339, 0.0398, 0.0485, 0.0622
    };

    // ------------------------------------------------------------------
    // Shared Bark axis
    // ------------------------------------------------------------------
    private static readonly double[] BarkAxisCache = Interp.Linspace(0.5, 26.5, NBands);

    // ------------------------------------------------------------------
    // Pre-computed gammatone filter coefficients (53 bands)
    // ------------------------------------------------------------------
    internal static readonly (Complex[] b, Complex[] a)[] GammatoneFilters = BuildAllGammatone();

    private static (Complex[] b, Complex[] a)[] BuildAllGammatone()
    {
        double[] centreFreq = BuildCentreFreq();
        var filters = new (Complex[] b, Complex[] a)[NBands];
        for (int z = 0; z < NBands; z++)
            filters[z] = GammatoneCoeffs(centreFreq[z], FilterOrder, Fs);
        return filters;
    }

    private static double[] BuildCentreFreq()
    {
        // ECMA-418-2, eq. 9: f(z) = (af_f0/c) * sinh(c * z)
        const double af_f0 = 81.9289;
        const double c = 0.1618;
        const double zStep = 0.5;

        double[] freq = new double[NBands];
        for (int z = 0; z < NBands; z++)
        {
            double zval = (z + 1) * zStep;
            freq[z] = (af_f0 / c) * Math.Sinh(c * zval);
        }
        return freq;
    }

    // ------------------------------------------------------------------
    // Gammatone filter design (ECMA-418-2, eq. 8–17)
    // ------------------------------------------------------------------
    private static (Complex[] b, Complex[] a) GammatoneCoeffs(double freq, int k, double fs)
    {
        const double af_f0 = 81.9289;
        const double c = 0.1618;

        // Bandwidth (eq. 10)
        double deltaF = Math.Sqrt(af_f0 * af_f0 + (c * freq) * (c * freq));

        // Time constant (eq. 8)
        // comb(2k-2, k-1) = C(2*(k-1), k-1)
        double binom = Binomial(2 * k - 2, k - 1);
        double tau = Math.Pow(2.0, -(2 * k - 1)) * binom / deltaF;

        // d coefficient
        double d = Math.Exp(-1.0 / (fs * tau));

        // am coefficients (eq. 14): a[0]=1, a[m] = (-d)^m * C(k, m) for m=1..k
        Complex[] amBase = new Complex[k + 1];
        amBase[0] = Complex.One;
        for (int m = 1; m <= k; m++)
            amBase[m] = Math.Pow(-d, m) * Binomial(k, m);

        // bm coefficients (eq. 15): b[m] = ((1-d)^k / denom) * d^m * em[m]
        int[] em = { 0, 1, 11, 11, 1 }; // length k for k=5 → indices 0..4
        double denom = 0.0;
        for (int i = 1; i < k; i++) denom += em[i] * Math.Pow(d, i);
        double scale = Math.Pow(1.0 - d, k) / denom;

        Complex[] bmBase = new Complex[k];
        for (int m = 0; m < k; m++)
            bmBase[m] = scale * Math.Pow(d, m) * em[m];

        // Modulate by positive exponential (eq. 16–17)
        // am_prim[m] = amBase[m] * exp(+j * 2π * freq * m / fs)
        // bm_prim[m] = bmBase[m] * exp(+j * 2π * freq * m / fs)
        Complex[] a = new Complex[k + 1];
        for (int m = 0; m <= k; m++)
        {
            double phase = 2.0 * Math.PI * freq * m / fs;
            a[m] = amBase[m] * new Complex(Math.Cos(phase), Math.Sin(phase));
        }

        Complex[] b = new Complex[k];
        for (int m = 0; m < k; m++)
        {
            double phase = 2.0 * Math.PI * freq * m / fs;
            b[m] = bmBase[m] * new Complex(Math.Cos(phase), Math.Sin(phase));
        }

        return (b, a);
    }

    // ------------------------------------------------------------------
    // Complex IIR filter (matches scipy.signal.lfilter with complex b, a)
    // Applied to real input; output is 2 × real(y).
    // ------------------------------------------------------------------
    internal static void ApplyGammatoneFilter(double[] signal, Complex[] b, Complex[] a, double[] output)
    {
        int n = signal.Length;
        int nb = b.Length; // k = 5
        int na = a.Length; // k+1 = 6

        Complex[] y = ArrayPool<Complex>.Shared.Rent(n);
        try
        {
            for (int i = 0; i < n; i++)
            {
                Complex acc = Complex.Zero;
                for (int k = 0; k < nb && i - k >= 0; k++)
                    acc += b[k] * signal[i - k];
                for (int m = 1; m < na && i - m >= 0; m++)
                    acc -= a[m] * y[i - m];
                y[i] = acc;
                output[i] = 2.0 * acc.Real;
            }
        }
        finally
        {
            ArrayPool<Complex>.Shared.Return(y);
        }
    }

    // ------------------------------------------------------------------
    // Preprocessing: fade-in window + zero padding (Section 5.1.2)
    // ------------------------------------------------------------------
    internal static (double[] padded, int nNew) Preprocessing(double[] signal, int sb, int sh)
    {
        int nSamples = signal.Length;

        // Zero padding
        int nZeroStart = sb;
        int nNew = (int)(sh * (Math.Ceiling((double)(nSamples + sh + sb) / sh) - 1));
        int nZeroEnd = nNew - nSamples;

        // Write directly into padded, applying the fade-in window during copy.
        // Avoids a full Clone() + Array.Copy().
        double[] padded = new double[nZeroStart + nSamples + nZeroEnd];
        const int nFadeIn = 240;
        for (int i = 0; i < nSamples; i++)
        {
            double v = signal[i];
            if (i < nFadeIn)
                v *= 0.5 - 0.5 * Math.Cos(Math.PI * i / nFadeIn);
            padded[nZeroStart + i] = v;
        }
        // trailing zeros already zero-initialized

        return (padded, nNew);
    }

    // ------------------------------------------------------------------
    // Segmentation into blocks (Section 5.1.5, eq. 18–20)
    // Returns blocks[nBlocks][sb] and times[nBlocks]
    // ------------------------------------------------------------------
    internal static (double[][] blocks, double[] times) Segment(
        double[] signal, int signalLength, int sb, int sh, int nNew, int iStart)
    {
        int nAll = signalLength;
        int lLast = (int)(Math.Ceiling((double)(nNew + sh) / sh) - 1);

        double[][] blocks = new double[lLast][];
        double[]   times  = new double[lLast];
        double tScale = 1.0 / Fs;
        double sbMinus1 = sb - 1.0;

        for (int l = 0; l < lLast; l++)
        {
            double linStart = l * sh + iStart;
            double linStop  = linStart + sb;

            double[] block = new double[sb];
            double tSum = 0.0;
            for (int k = 0; k < sb; k++)
            {
                int si = (int)(linStart + k * (linStop - linStart) / sbMinus1);
                block[k] = (si >= 0 && si < nAll) ? signal[si] : 0.0;
                tSum += si * tScale;
            }
            blocks[l] = block;
            times[l]  = tSum / sb;
        }

        return (blocks, times);
    }

    // ------------------------------------------------------------------
    // Half-wave rectification + RMS (eq. 22)
    // rms = sqrt(2 * mean(max(x,0)^2))
    // ------------------------------------------------------------------
    internal static double[] RectifiedRms(double[][] blocks, int sb)
    {
        int nBlocks = blocks.Length;
        double[] rms = new double[nBlocks];
        for (int l = 0; l < nBlocks; l++)
        {
            double sumSq = 0.0;
            double[] bl = blocks[l];
            for (int k = 0; k < sb; k++)
            {
                double v = bl[k] > 0 ? bl[k] : 0.0;
                sumSq += v * v;
            }
            rms[l] = Math.Sqrt(2.0 * sumSq / sb);
        }
        return rms;
    }

    // ------------------------------------------------------------------
    // Compressive nonlinearity (eq. 23)
    // ------------------------------------------------------------------
    internal static double[] Nonlinearity(double[] p)
    {
        const double p0    = 2e-5;
        const double cN    = 0.0211668;
        const double alpha = 1.50;

        double[] vi    = { 1.0, 0.6602, 0.0864, 0.6384, 0.0328, 0.4068, 0.2082, 0.3994, 0.6434 };
        double[] thresh = { 0.0, 15.0, 25.0, 35.0, 45.0, 55.0, 65.0, 75.0, 85.0 };

        double[] pti = new double[thresh.Length];
        for (int i = 0; i < pti.Length; i++)
            pti[i] = p0 * Math.Pow(10.0, thresh[i] / 20.0);

        double[] aPrime = new double[p.Length];
        for (int n = 0; n < p.Length; n++) aPrime[n] = 1.0;

        for (int i = 1; i < 9; i++)
        {
            double exp = (vi[i] - vi[i - 1]) / alpha;
            for (int n = 0; n < p.Length; n++)
                aPrime[n] *= Math.Pow(1.0 + Math.Pow(p[n] / pti[i], alpha), exp);
        }

        for (int n = 0; n < p.Length; n++)
            aPrime[n] *= cN * p[n] / p0;

        return aPrime;
    }

    // ------------------------------------------------------------------
    // Public entry point
    // ------------------------------------------------------------------

    /// <summary>
    /// Computes ECMA-418-2 specific and total loudness.
    /// </summary>
    /// <param name="signal">Time signal [Pa] at 48 kHz.</param>
    /// <param name="fs">Sampling frequency [Hz] — must be 48000 or will be resampled.</param>
    /// <param name="sb">Block size [samples] (default 2048).</param>
    /// <param name="sh">Hop size [samples] (default 1024).</param>
    /// <returns><see cref="LoudnessEcmaResult"/>.</returns>
    public static LoudnessEcmaResult Compute(
        ReadOnlySpan<double> signal, int fs, int sb = 2048, int sh = 1024)
    {
        double[] sig;
        if (fs != 48000)
        {
            MosqitoLog.Warn("[Warning] Signal resampled to 48 kHz for ECMA-418-2 loudness.");
            sig = Resample.Apply(signal, fs, 48000);
        }
        else
        {
            sig = signal.ToArray();
        }

        // 5.1.2 Preprocessing (windowing + zero padding)
        var (padded, nNew) = Preprocessing(sig, sb, sh);

        // 5.1.3 Outer and middle ear filter
        double[] earFiltered = SosFilter.Process(EarSos, padded);

        // i_start = sb_max - sb = 0 (since sb is scalar here)
        int iStart = 0;

        // 5.1.4–5.1.9 Per-band processing
        double[][] NSpec    = new double[NBands][];
        double[][] timeArrays = new double[NBands][];

        int earLen = earFiltered.Length;
        Parallel.For(0, NBands, z =>
        {
            // 5.1.4 Gammatone filter (complex IIR) — bandPass rented to avoid per-band heap alloc
            var (bGamma, aGamma) = GammatoneFilters[z];
            double[] bandPass = ArrayPool<double>.Shared.Rent(earLen);
            try
            {
                ApplyGammatoneFilter(earFiltered, bGamma, aGamma, bandPass);

                // 5.1.5 Segmentation — pass earLen so the over-sized rented buffer is bounded correctly
                var (blocks, times) = Segment(bandPass, earLen, sb, sh, nNew, iStart);

                // 5.1.6 Rectification + RMS (eq. 22)
                double[] rms = RectifiedRms(blocks, sb);

                // 5.1.7–5.1.8 Nonlinearity (eq. 23)
                double[] aPrime = Nonlinearity(rms);

                // 5.1.9 Apply threshold
                double ltq = LtqZ[z];
                for (int l = 0; l < aPrime.Length; l++)
                    if (aPrime[l] < ltq) aPrime[l] = ltq;

                // N'(z, t) = aPrime - ltq
                double[] nPrime = new double[aPrime.Length];
                for (int l = 0; l < aPrime.Length; l++)
                    nPrime[l] = aPrime[l] - ltq;

                NSpec[z]      = nPrime;
                timeArrays[z] = times;
            }
            finally
            {
                ArrayPool<double>.Shared.Return(bandPass);
            }
        });

        // Eq. 116: N_time = sum_z N'_z * delta_z
        // All bands should have the same length (same sb/sh)
        int nTime = NSpec[0].Length;
        double[] NTime = new double[nTime];
        for (int z = 0; z < NBands; z++)
        {
            double[] ns = NSpec[z];
            int nt = Math.Min(ns.Length, nTime);
            for (int t = 0; t < nt; t++)
                NTime[t] += ns[t] * DeltaZ;
        }

        // Eq. 117: N = (mean(N_time^e))^(1/e) where e = 1/log10(2)
        double e = 1.0 / Math.Log10(2.0);
        double sumPow = 0.0;
        for (int t = 0; t < nTime; t++)
            sumPow += Math.Pow(NTime[t], e);
        double N = Math.Pow(sumPow / nTime, 1.0 / e);

        return new LoudnessEcmaResult(N, NTime, NSpec, BarkAxisCache, timeArrays);
    }

    // ------------------------------------------------------------------
    // Helper: binomial coefficient C(n, k)
    // ------------------------------------------------------------------
    private static double Binomial(int n, int k)
    {
        if (k < 0 || k > n) return 0;
        if (k == 0 || k == n) return 1;
        k = Math.Min(k, n - k);
        double result = 1.0;
        for (int i = 0; i < k; i++)
        {
            result *= (n - i);
            result /= (i + 1);
        }
        return result;
    }
}
