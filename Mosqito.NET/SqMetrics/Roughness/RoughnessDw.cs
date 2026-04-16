using System.Buffers;
using System.Numerics;
using Mosqito.Conversion;
using Mosqito.Dsp;
using Mosqito.Io;
using Mosqito.SoundLevelMeter;
using Mosqito.SqMetrics.Loudness;

namespace Mosqito.SqMetrics.Roughness;

/// <summary>
/// Result of Daniel–Weber roughness computation.
/// </summary>
/// <param name="R">Roughness per segment [asper], length nSeg.</param>
/// <param name="RSpecific">Specific roughness [asper/Bark], shape (47, nSeg).</param>
/// <param name="BarkAxis">Bark axis, length 47 (0.5 … 23.5 Bark).</param>
/// <param name="TimeAxis">Centre time of each segment [s], length nSeg.</param>
public readonly record struct RoughnessDwResult(
    double[] R,
    double[,] RSpecific,
    double[] BarkAxis,
    double[] TimeAxis);

/// <summary>
/// Daniel and Weber psychoacoustical roughness (1997).
/// Ported from MoSQITo <c>roughness_dw.py</c> and <c>roughness_dw_freq.py</c>.
/// </summary>
public static class RoughnessDw
{
    // -----------------------------------------------------------------------
    // Bark channel centres: zi = (1..47)/2 → [0.5, 1.0, …, 23.5]
    // -----------------------------------------------------------------------
    private static readonly double[] BarkAxisCache = BuildBarkAxis();

    private static double[] BuildBarkAxis()
    {
        double[] b = new double[47];
        for (int i = 0; i < 47; i++) b[i] = (i + 1) / 2.0;
        return b;
    }

    // -----------------------------------------------------------------------
    // Time-domain entry point
    // -----------------------------------------------------------------------

    /// <summary>
    /// Computes Daniel–Weber roughness from a time signal.
    /// </summary>
    /// <param name="signal">Time signal [Pa] at <paramref name="fs"/> Hz.</param>
    /// <param name="fs">Sampling frequency [Hz].</param>
    /// <param name="overlap">Overlap fraction for 200 ms windows (default 0.5).</param>
    /// <returns><see cref="RoughnessDwResult"/>.</returns>
    public static RoughnessDwResult Compute(
        ReadOnlySpan<double> signal, int fs, double overlap = 0.5)
    {
        int nperseg = (int)(0.2 * fs);
        int noverlap = (int)(overlap * nperseg);

        var (blocks, timeAxis) = TimeSegmentation.Segment(signal, fs, nperseg, noverlap);
        int nSeg = blocks.GetLength(1);
        int nSamples = blocks.GetLength(0);

        // Frequency axis (one-sided, starting at fs/nperseg)
        int halfN = nperseg / 2;
        double[] freqAxis = new double[halfN];
        for (int i = 0; i < halfN; i++) freqAxis[i] = (i + 1) * (double)fs / nperseg;

        // Pre-compute weighting tables
        double[] gzi     = GziWeighting();
        double[,] hWeight = HWeighting(nperseg, fs);

        double[] R    = new double[nSeg];
        double[,] RSpec = new double[47, nSeg];
        double[] segBuf = new double[nSamples];

        for (int seg = 0; seg < nSeg; seg++)
        {
            for (int s = 0; s < nSamples; s++) segBuf[s] = blocks[s, seg];

            // Compute Blackman-windowed COMPLEX spectrum (one-sided).
            // Python passes the complex FFT to _roughness_dw_main_calc to preserve
            // phase information — this is essential for accurate excitation patterns.
            Complex[] specComplex = ComputeComplexSpectrum(segBuf, nperseg);

            double[] r = MainCalc(specComplex, freqAxis, fs, gzi, hWeight);
            R[seg] = r[47]; // last element is overall R; first 47 are R_spec
            for (int b = 0; b < 47; b++) RSpec[b, seg] = r[b];
        }

        return new RoughnessDwResult(R, RSpec, BarkAxisCache, timeAxis);
    }

    // -----------------------------------------------------------------------
    // Frequency-domain entry point
    // -----------------------------------------------------------------------

    /// <summary>
    /// Computes Daniel–Weber roughness from a fine-band amplitude spectrum.
    /// </summary>
    /// <param name="spectrum">One-sided amplitude spectrum (real magnitudes).</param>
    /// <param name="freqs">Frequency axis [Hz].</param>
    public static (double R, double[] RSpecific, double[] BarkAxis) ComputeFromSpectrum(
        ReadOnlySpan<double> spectrum, ReadOnlySpan<double> freqs)
    {
        if (spectrum.Length != freqs.Length)
            throw new ArgumentException("spectrum and freqs must have the same length.");

        int nperseg = spectrum.Length;
        // Derive fs from frequency resolution
        double df = 0.0;
        for (int i = 1; i < freqs.Length; i++) df += freqs[i] - freqs[i - 1];
        df /= (freqs.Length - 1);
        int fs = (int)(2 * nperseg * df + 0.5);

        double[] gzi      = GziWeighting();
        double[,] hWeight = HWeighting(2 * nperseg, fs);

        // For freq-domain path, input is amplitude spectrum (real, phase=0)
        Complex[] specComplex = new Complex[spectrum.Length];
        for (int i = 0; i < spectrum.Length; i++)
            specComplex[i] = new Complex(spectrum[i], 0.0);

        double[] freqArr = freqs.ToArray();

        double[] r = MainCalc(specComplex, freqArr, fs, gzi, hWeight);
        double[] rSpec = new double[47];
        for (int b = 0; b < 47; b++) rSpec[b] = r[b];

        return (r[47], rSpec, BarkAxisCache);
    }

    // -----------------------------------------------------------------------
    // Core algorithm
    // -----------------------------------------------------------------------

    /// <summary>
    /// Returns array of length 48: elements [0..46] = R_spec per channel; [47] = total R.
    /// </summary>
    private static double[] MainCalc(
        Complex[] spec, double[] freqAxis, int fs,
        double[] gzi, double[,] hWeight)
    {
        int m = spec.Length;  // one-sided length
        int n = 2 * m;        // two-sided

        // Rent all large temporary buffers from the pool.
        // Fourier.Forward/Inverse requires exact-length arrays, so we track sizes carefully.
        Complex[] spec2      = ArrayPool<Complex>.Shared.Rent(n);
        Complex[] specW      = ArrayPool<Complex>.Shared.Rent(n);
        double[]  module     = ArrayPool<double>.Shared.Rent(m);
        double[]  specDb     = ArrayPool<double>.Shared.Rent(m);
        Complex[] exc        = ArrayPool<Complex>.Shared.Rent(n);
        // ifftBuf and envBuf must be exactly n for in-place FFT — allocate normally.
        Complex[] ifftBuf    = new Complex[n];
        Complex[] envBuf     = new Complex[n];

        try
        {
        // Build two-sided complex spectrum: [spec, reversed(spec)]
        // Matches Python: concatenate((spec, spec[len(spec)::-1]))
        for (int i = 0; i < m; i++) spec2[i] = spec[i];
        for (int i = 0; i < m; i++) spec2[m + i] = spec[m - 1 - i];

        // Bark axis for the one-sided part
        double[] barkAxis = FreqBark.Freq2Bark(freqAxis);

        // Zwicker outer/inner ear transfer (a0 factor)
        double[] a0dB  = EarFilterCoeff(barkAxis);
        double[] a0Real = ArrayPool<double>.Shared.Rent(m);
        try
        {
            for (int i = 0; i < m; i++) a0Real[i] = AmpDb.Db2Amp(a0dB[i], reference: 1.0);

            // Apply a0 to first half; second half stays Complex.Zero
            for (int i = 0; i < m; i++) specW[i] = a0Real[i] * spec2[i];
            for (int i = m; i < n; i++) specW[i] = Complex.Zero;
        }
        finally { ArrayPool<double>.Shared.Return(a0Real); }

        for (int i = 0; i < m; i++) module[i] = specW[i].Magnitude;

        for (int i = 0; i < m; i++)
            specDb[i] = module[i] > 0 ? AmpDb.Amp2Db(module[i], reference: 2e-5) : -999.0;

        double[] threshold = LTQ.Compute(barkAxis, "roughness");

        // Find audible components
        var audibleIdx = new System.Collections.Generic.List<int>();
        for (int i = 0; i < m; i++)
            if (specDb[i] > threshold[i]) audibleIdx.Add(i);
        int nAud = audibleIdx.Count;

        // Terhardt slope parameters
        double s1 = -27.0;
        double[] s2 = new double[nAud];
        for (int k = 0; k < nAud; k++)
        {
            int idx = audibleIdx[k];
            double v = -24.0 - (230.0 / freqAxis[idx]) + (0.2 * specDb[idx]);
            s2[k] = Math.Min(v, 0.0);
        }

        // Channel centres in normalised bin space
        int nChannel = 47;
        double[] zi    = BarkAxisCache; // [0.5..23.5]
        double[] zbFreq = FreqBark.Bark2Freq(zi);
        double[] zb    = new double[nChannel];
        for (int i = 0; i < nChannel; i++) zb[i] = zbFreq[i] * n / fs;

        // Minimum excitation level per channel
        double[] minExcitDb = new double[nChannel];
        double[] nZ = new double[m];
        for (int i = 0; i < m; i++) nZ[i] = i + 1;
        for (int i = 0; i < nChannel; i++)
            minExcitDb[i] = Interp.Linear(zb[i], nZ, threshold);

        // Channel membership for each audible component
        double[] chLow  = new double[nAud];
        double[] chHigh = new double[nAud];
        for (int i = 0; i < nAud; i++)
        {
            double b = barkAxis[audibleIdx[i]];
            chLow[i]  = Math.Floor(2.0 * b) - 1.0;
            chHigh[i] = Math.Ceiling(2.0 * b) - 1.0;
        }

        // Excitation slopes [nAud x nChannel]
        double[,] slopes = new double[nAud, nChannel];
        for (int k = 0; k < nAud; k++)
        {
            double levDb = specDb[audibleIdx[k]];
            double b     = barkAxis[audibleIdx[k]];

            for (int j = 0; j <= (int)chLow[k]; j++)
            {
                double sl = s1 * (b - (j + 1) * 0.5) + levDb;
                if (sl > minExcitDb[j])
                    slopes[k, j] = AmpDb.Db2Amp(sl, reference: 2e-5);
            }
            for (int j = (int)chHigh[k]; j < nChannel; j++)
            {
                double sl = s2[k] * ((j + 1) * 0.5 - b) + levDb;
                if (sl > minExcitDb[j])
                    slopes[k, j] = AmpDb.Db2Amp(sl, reference: 2e-5);
            }
        }

        // Per-channel processing
        double[,] hBP      = new double[nChannel, n];
        double[]  modDepth = new double[nChannel];

        for (int ch = 0; ch < nChannel; ch++)
        {
            // Build excitation spectrum for this channel
            for (int i = 0; i < n; i++) exc[i] = Complex.Zero;
            for (int j = 0; j < nAud; j++)
            {
                int ind = audibleIdx[j];
                double ampl;
                if (chLow[j] == ch || chHigh[j] == ch)
                    ampl = 1.0;
                else if (chHigh[j] > ch)
                    ampl = module[ind] > 0 ? slopes[j, ch + 1] / module[ind] : 0.0;
                else
                    ampl = module[ind] > 0 ? slopes[j, ch - 1] / module[ind] : 0.0;

                exc[ind] = ampl * specW[ind];
            }

            // IFFT → temporal specific excitation (reuse ifftBuf, exact length n)
            exc.AsSpan(0, n).CopyTo(ifftBuf);
            MathNet.Numerics.IntegralTransforms.Fourier.Inverse(
                ifftBuf,
                MathNet.Numerics.IntegralTransforms.FourierOptions.AsymmetricScaling);

            double h0 = 0.0;
            for (int i = 0; i < n; i++)
            {
                double v = n * ifftBuf[i].Real;
                hBP[ch, i] = v;
                h0 += Math.Abs(v);
            }
            h0 /= n;

            // Envelope spectrum: FFT of (|temporal_excitation| - h0)
            for (int i = 0; i < n; i++) envBuf[i] = new Complex(Math.Abs(hBP[ch, i]) - h0, 0.0);
            MathNet.Numerics.IntegralTransforms.Fourier.Forward(
                envBuf,
                MathNet.Numerics.IntegralTransforms.FourierOptions.AsymmetricScaling);

            for (int i = 0; i < n; i++) envBuf[i] *= hWeight[ch, i];

            MathNet.Numerics.IntegralTransforms.Fourier.Inverse(
                envBuf,
                MathNet.Numerics.IntegralTransforms.FourierOptions.AsymmetricScaling);
            for (int i = 0; i < n; i++) hBP[ch, i] = 2.0 * envBuf[i].Real;

            double rms2 = 0.0;
            for (int i = 0; i < n; i++) { double v = hBP[ch, i]; rms2 += v * v; }
            double hBpRms = Math.Sqrt(rms2 / n);

            if (h0 > 0.0)
                modDepth[ch] = Math.Min(hBpRms / h0, 1.0);
        }

        // Cross-correlation between channels i and i+2
        double[] ki = new double[47];
        for (int i = 0; i < 45; i++)
        {
            bool allZeroA = true, allZeroB = true;
            for (int t = 0; t < n; t++)
            {
                if (hBP[i, t] != 0.0) allZeroA = false;
                if (hBP[i + 2, t] != 0.0) allZeroB = false;
            }
            if (!allZeroA && !allZeroB)
                ki[i] = Corrcoef(hBP, i, i + 2, n);
        }

        // Specific roughness — use x*x instead of Math.Pow(x, 2)
        double[] rSpec = new double[47];
        { double v = modDepth[0] * ki[0]; rSpec[0] = gzi[0] * v * v; }
        { double v = modDepth[1] * ki[1]; rSpec[1] = gzi[1] * v * v; }
        for (int i = 2; i < 45; i++)
        { double v = modDepth[i] * ki[i] * ki[i - 2]; rSpec[i] = gzi[i] * v * v; }
        { double v = modDepth[45] * ki[43]; rSpec[45] = gzi[45] * v * v; }
        { double v = modDepth[46] * ki[44]; rSpec[46] = gzi[46] * v * v; }

        double totalR = 0.0;
        for (int i = 0; i < 47; i++) totalR += rSpec[i];
        totalR *= 0.25;

        double[] result = new double[48];
        for (int i = 0; i < 47; i++) result[i] = rSpec[i];
        result[47] = totalR;
        return result;

        } // end try
        finally
        {
            ArrayPool<Complex>.Shared.Return(spec2);
            ArrayPool<Complex>.Shared.Return(specW);
            ArrayPool<double>.Shared.Return(module);
            ArrayPool<double>.Shared.Return(specDb);
            ArrayPool<Complex>.Shared.Return(exc);
        }
    }

    // -----------------------------------------------------------------------
    // Complex spectrum helper
    // -----------------------------------------------------------------------

    /// <summary>
    /// Computes the one-sided complex FFT spectrum using a Blackman window,
    /// normalised by sum(window) and scaled by 1.42 — matching Python's
    /// <c>comp_spectrum(sig, fs, window='blackman', db=False)</c>.
    /// Returns m = nfft/2 complex bins (indices 0..m-1, freqs 1*df..m*df).
    /// </summary>
    private static Complex[] ComputeComplexSpectrum(double[] seg, int nfft)
    {
        // Build and normalise window
        double[] win = new double[nfft];
        Windows.FillBlackman(win);
        Windows.NormaliseBySumInPlace(win);

        // Apply window
        Complex[] buf = new Complex[nfft];
        int copyLen = Math.Min(nfft, seg.Length);
        for (int i = 0; i < copyLen; i++) buf[i] = new Complex(seg[i] * win[i], 0.0);

        // Forward FFT (AsymmetricScaling = no scaling on forward, 1/N on inverse)
        MathNet.Numerics.IntegralTransforms.Fourier.Forward(
            buf,
            MathNet.Numerics.IntegralTransforms.FourierOptions.AsymmetricScaling);

        // Extract one-sided part (indices 0..m-1, matching Python's [0:nfft//2])
        // and apply 1.42 amplitude correction
        int m = nfft / 2;
        Complex[] result = new Complex[m];
        for (int i = 0; i < m; i++) result[i] = buf[i] * 1.42;
        return result;
    }

    // -----------------------------------------------------------------------
    // Weighting tables
    // -----------------------------------------------------------------------

    /// <summary>Gzi weighting (Aures) for 47 channels.</summary>
    private static double[] GziWeighting()
    {
        double[] grX = new double[25];
        for (int i = 0; i < 25; i++) grX[i] = i;

        double[] grY =
        {
            0.15, 0.26, 0.38, 0.47, 0.54, 0.65, 0.76, 0.83, 0.90, 0.98,
            0.98, 0.90, 0.80, 0.70, 0.62, 0.54, 0.49, 0.43, 0.39, 0.35,
            0.30, 0.30, 0.30, 0.30, 0.30,
        };

        double[] gzi = new double[47];
        for (int i = 0; i < 47; i++)
            gzi[i] = Interp.Linear((i + 1) / 2.0, grX, grY);
        return gzi;
    }

    /// <summary>
    /// H weighting functions for each 1-Bark-wide channel (shape: 47 × n).
    /// Ported from <c>_H_weighting.py</c>.
    /// </summary>
    private static double[,] HWeighting(int n, int fs)
    {
        double[,] H = new double[47, n];
        int cut = 2;

        // Helper: fill row from x/y table up to cutoff frequency.
        // Python: freq = j * fs/n (NOT (j-cut)*fs/n — that was a bug).
        // H21 and H42 share the same j range as H16 (502 Hz) in Python due to
        // variable reuse — replicated here intentionally to match Python output.
        void FillRow(int row, double[] hx, double[] hy, double maxFreq)
        {
            int last = (int)Math.Floor((maxFreq / fs) * n);
            for (int j = cut; j < last && j < n; j++)
            {
                double freq = j * (double)fs / n;   // match Python: freq = j * fs/n
                H[row, j] = Interp.Linear(freq, hx, hy);
            }
        }

        double[] H2_x  = { 0, 17, 23, 25, 32, 37, 48, 67, 90, 114, 171, 206, 247, 294, 358 };
        double[] H2_y  = { 0, 0.8, 0.95, 0.975, 1, 0.975, 0.9, 0.8, 0.7, 0.6, 0.4, 0.3, 0.2, 0.1, 0 };
        FillRow(1, H2_x, H2_y, 358);

        double[] H5_x  = { 0, 32, 43, 56, 69, 92, 120, 142, 165, 231, 277, 331, 397, 502 };
        double[] H5_y  = { 0, 0.8, 0.95, 1, 0.975, 0.9, 0.8, 0.7, 0.6, 0.4, 0.3, 0.2, 0.1, 0 };
        FillRow(4, H5_x, H5_y, 502);

        double[] H16_x = { 0, 23.5, 34, 47, 56, 63, 79, 100, 115, 135, 159, 172, 194, 215, 244, 290, 348, 415, 500, 645 };
        double[] H16_y = { 0, 0.4, 0.6, 0.8, 0.9, 0.95, 1, 0.975, 0.95, 0.9, 0.85, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0 };
        FillRow(15, H16_x, H16_y, 502);

        double[] H21_x = { 0, 19, 44, 52.5, 58, 75, 101.5, 114.5, 132.5, 143.5, 165.5, 197.5, 241, 290, 348, 415, 500, 645 };
        double[] H21_y = { 0, 0.4, 0.8, 0.9, 0.95, 1, 0.95, 0.9, 0.85, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0 };
        // Python reuses the j-range from H16 (502 Hz) for H21 — matching that here.
        FillRow(20, H21_x, H21_y, 502);

        double[] H42_x = { 0, 15, 41, 49, 53, 64, 71, 88, 94, 106, 115, 137, 180, 238, 290, 348, 415, 500, 645 };
        double[] H42_y = { 0, 0.4, 0.8, 0.9, 0.965, 0.99, 1, 0.95, 0.9, 0.85, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0 };
        // Python reuses the j-range from H16 (502 Hz) for H42 — matching that here.
        FillRow(41, H42_x, H42_y, 502);

        // Copy rows according to the article rules
        // H[0]=H[2]=H[3]=H[1]
        CopyRow(H, 1, 0, n);
        CopyRow(H, 1, 2, n);
        CopyRow(H, 1, 3, n);
        // i=5..15 → H[4]
        for (int i = 5; i < 15; i++) CopyRow(H, 4, i, n);
        // i=16..20 → H[15]
        for (int i = 16; i < 20; i++) CopyRow(H, 15, i, n);
        // i=21..41 → H[20]
        for (int i = 21; i < 41; i++) CopyRow(H, 20, i, n);
        // i=42..47 → H[41]
        for (int i = 42; i < 47; i++) CopyRow(H, 41, i, n);

        return H;
    }

    private static void CopyRow(double[,] H, int srcRow, int dstRow, int n)
    {
        for (int j = 0; j < n; j++) H[dstRow, j] = H[srcRow, j];
    }

    /// <summary>
    /// Ear filter coefficients (Zwicker a0) as dB values along a Bark axis.
    /// Linearly interpolated from figure 8.18 of Zwicker &amp; Fastl 1990.
    /// </summary>
    private static double[] EarFilterCoeff(ReadOnlySpan<double> barkAxis)
    {
        double[] xp =
        {
            0, 10, 12, 13, 14, 15, 16, 16.5, 17, 18, 18.5,
            19, 20, 21, 21.5, 22, 22.5, 23, 23.5, 24, 25, 26,
        };
        double[] yp =
        {
            0, 0, 1.15, 2.31, 3.85, 5.62, 6.92, 7.38, 6.92, 4.23, 2.31,
            0, -1.43, -2.59, -3.57, -5.19, -7.41, -11.3, -20, -40, -130, -999,
        };

        double[] result = new double[barkAxis.Length];
        for (int i = 0; i < barkAxis.Length; i++)
            result[i] = Interp.Linear(barkAxis[i], xp, yp);
        return result;
    }

    // -----------------------------------------------------------------------
    // Statistics helper
    // -----------------------------------------------------------------------

    /// <summary>Pearson correlation coefficient between rows <paramref name="r1"/> and <paramref name="r2"/>.</summary>
    private static double Corrcoef(double[,] mat, int r1, int r2, int len)
    {
        double mean1 = 0, mean2 = 0;
        for (int i = 0; i < len; i++) { mean1 += mat[r1, i]; mean2 += mat[r2, i]; }
        mean1 /= len; mean2 /= len;

        double cov = 0, var1 = 0, var2 = 0;
        for (int i = 0; i < len; i++)
        {
            double d1 = mat[r1, i] - mean1;
            double d2 = mat[r2, i] - mean2;
            cov  += d1 * d2;
            var1 += d1 * d1;
            var2 += d2 * d2;
        }
        double denom = Math.Sqrt(var1 * var2);
        return denom > 0 ? cov / denom : 0.0;
    }
}
