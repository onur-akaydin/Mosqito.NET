using System.Buffers;
using System.Numerics;
using System.Threading.Tasks;
using Mosqito.Dsp;
using Mosqito.SqMetrics.Loudness;

namespace Mosqito.SqMetrics.Roughness;

/// <summary>
/// ECMA-418-2 (2nd Ed., 2022) psychoacoustical roughness.
/// Ported from MoSQITo <c>roughness_ecma.py</c>.
/// </summary>
public static class RoughnessEcma
{
    private const int CBF = 53;
    private const int Fs  = 48000;

    // Roughness uses larger block/hop sizes than loudness
    private const int Sb = 16384;
    private const int Sh = 4096;

    // Envelope downsampling: 48000 / 32 = 1500 Hz
    private const int DownsamplingFactor = 32;
    // sbb = 512 (block length at 1500 Hz)
    private const int Sbb = 512;

    // Von Hann normalised by sqrt(0.375), cached as a static readonly
    private static readonly double[] _vonHann = BuildVonHann(Sbb);

    private static double[] BuildVonHann(int n)
    {
        double[] w = new double[n];
        double norm = 1.0 / Math.Sqrt(0.375);
        for (int k = 0; k < n; k++)
            w[k] = (0.5 - 0.5 * Math.Cos(2.0 * Math.PI * k / n)) * norm;
        return w;
    }

    // ------------------------------------------------------------------
    // Public entry point
    // ------------------------------------------------------------------

    /// <summary>
    /// Computes ECMA-418-2 roughness.
    /// </summary>
    public static (double R, double[] RTime, double[,] RSpecific, double[] BarkAxis, double[] TimeAxis)
        Compute(ReadOnlySpan<double> signal, int fs)
    {
        double[] sig;
        if (fs != Fs)
        {
            MosqitoLog.Warn("[Warning] Signal resampled to 48 kHz for ECMA-418-2 roughness.");
            sig = Resample.Apply(signal, fs, Fs);
        }
        else
        {
            sig = signal.ToArray();
        }

        double duration = sig.Length / (double)Fs;

        // STEP 1 — Preprocessing (fade-in + zero-padding)
        var (padded, nNew) = LoudnessEcma.Preprocessing(sig, Sb, Sh);

        // STEP 2 — Outer/middle ear filter
        double[] earFiltered = SosFilter.Process(LoudnessEcma.EarSos, padded);
        int earLen = earFiltered.Length;

        // STEP 3 — Per-band gammatone filtering + segmentation + loudness (parallel over z)
        int iStart = 0;
        double[][] nSpec      = new double[CBF][];
        double[][][] bandBlocks = new double[CBF][][];
        double[][] timeArrays  = new double[CBF][];

        Parallel.For(0, CBF,
            () => ArrayPool<double>.Shared.Rent(earLen),
            (z, _, bp) =>
            {
                var (bGamma, aGamma) = LoudnessEcma.GammatoneFilters[z];
                LoudnessEcma.ApplyGammatoneFilter(earFiltered, bGamma, aGamma, bp);

                var (blocks, times) = LoudnessEcma.Segment(bp, earLen, Sb, Sh, nNew, iStart);
                bandBlocks[z] = blocks;
                timeArrays[z] = times;

                double[] rms    = LoudnessEcma.RectifiedRms(blocks, Sb);
                double[] aPrime = LoudnessEcma.Nonlinearity(rms);
                double   ltq    = LoudnessEcma.LtqZ[z];
                double[] nPrime = new double[aPrime.Length];
                for (int l = 0; l < aPrime.Length; l++)
                {
                    if (aPrime[l] < ltq) aPrime[l] = ltq;
                    nPrime[l] = aPrime[l] - ltq;
                }
                nSpec[z] = nPrime;
                return bp;
            },
            bp => ArrayPool<double>.Shared.Return(bp));

        int L = nSpec[0].Length; // number of time blocks

        // N_specific max across bands for each time block
        double[] nSpecMax = new double[L];
        for (int l = 0; l < L; l++)
        {
            double maxV = 0.0;
            for (int z = 0; z < CBF; z++)
                if (nSpec[z][l] > maxV) maxV = nSpec[z][l];
            nSpecMax[l] = maxV;
        }

        // STEP 4 — Hilbert envelope + decimate to 1500 Hz (parallel over z)
        double[] envsDown = new double[L * CBF * Sbb];

        Parallel.For(0, CBF, z =>
        {
            for (int l = 0; l < L; l++)
            {
                double[] block = bandBlocks[z][l];
                // Hilbert.Envelope uses [ThreadStatic] scratch — thread-safe
                double[] env = Hilbert.Envelope(block);
                double[] d8  = DecimateSimple(env, 8);
                double[] d32 = DecimateSimple(d8, 4);
                int copyLen = Math.Min(d32.Length, Sbb);
                int baseIdx = (l * CBF + z) * Sbb;
                for (int k = 0; k < copyLen; k++) envsDown[baseIdx + k] = d32[k];
            }
        });

        // STEP 5 — Von Hann window and scaled power spectrum (parallel over l)
        double[] vhann     = _vonHann;
        double[,,] phiE0   = new double[L, CBF, 1];
        double[,,] phiERaw = new double[L, CBF, Sbb / 2];

        Parallel.For(0, L,
            () => (ewBuf: new double[Sbb], specBuf: new Complex[Sbb]),
            (l, _, state) =>
            {
                double[] ewBuf   = state.ewBuf;
                Complex[] specBuf = state.specBuf;
                for (int z = 0; z < CBF; z++)
                {
                    int baseIdx = (l * CBF + z) * Sbb;
                    double sumSq = 0.0;
                    for (int k = 0; k < Sbb; k++)
                    {
                        double v = envsDown[baseIdx + k] * vhann[k];
                        ewBuf[k] = v;
                        sumSq += v * v;
                    }
                    phiE0[l, z, 0] = sumSq;

                    for (int k = 0; k < Sbb; k++) specBuf[k] = new Complex(ewBuf[k], 0.0);
                    MathNet.Numerics.IntegralTransforms.Fourier.Forward(
                        specBuf,
                        MathNet.Numerics.IntegralTransforms.FourierOptions.AsymmetricScaling);

                    for (int k = 0; k < Sbb / 2; k++)
                    {
                        double mag = specBuf[k].Magnitude / 2.0 * Math.Sqrt(2.0);
                        phiERaw[l, z, k] = mag * mag;
                    }
                }
                return state;
            },
            _ => { });

        // Apply scaling: phi_E = (N'^2 / (N_max * phi_E0)) * dft — parallel over l
        double[,,] phiE = new double[L, CBF, Sbb / 2];
        Parallel.For(0, L, l =>
        {
            for (int z = 0; z < CBF; z++)
            {
                double den   = nSpecMax[l] * phiE0[l, z, 0];
                double scale = den != 0 ? nSpec[z][l] * nSpec[z][l] / den : 0.0;
                for (int k = 0; k < Sbb / 2; k++)
                    phiE[l, z, k] = scale * phiERaw[l, z, k];
            }
        });

        // STEP 6 — Noise reduction
        double[,,] PhiE = NoiseReduction(phiE, L, CBF, Sbb / 2);

        // STEP 7 — Centre frequencies for weighting
        double[] centreFreq = BuildCentreFreq();
        double[] fmax  = ComputeFmax(centreFreq);
        double[] rmax  = ComputeRmax(centreFreq);
        double[] q2Hi  = ComputeQ2High(centreFreq);
        double[] q2Lo  = ComputeQ2Low(centreFreq);

        // STEP 8 — Spectral weighting per time/band block
        double[,] amplitude = new double[L, CBF];
        double[] sliceBuf = new double[Sbb / 2];
        for (int l = 0; l < L; l++)
        for (int z = 0; z < CBF; z++)
        {
            for (int k = 0; k < Sbb / 2; k++) sliceBuf[k] = PhiE[l, z, k];

            var (fp, Ai) = PeakPicking(sliceBuf);
            int nPeak = fp.Length;

            if (nPeak == 0)
            {
                amplitude[l, z] = 0.0;
            }
            else
            {
                double[] aiTilde = new double[nPeak];
                for (int i = 0; i < nPeak; i++)
                    aiTilde[i] = HighModRateWeighting(fp[i], Ai[i], fmax[z], rmax[z], q2Hi[z]);

                var (modRate, aHat) = EstimateFundModRate(fp, aiTilde);
                amplitude[l, z] = LowModRateWeighting(modRate, aHat, fmax[z], q2Lo[z]);
            }
        }

        // Zero below threshold
        for (int l = 0; l < L; l++)
        for (int z = 0; z < CBF; z++)
            if (amplitude[l, z] < 0.074376) amplitude[l, z] = 0.0;

        // STEP 9 — Interpolation to 50 Hz (parallel over z)
        double[] timeAxis0 = timeArrays[0];
        var (amp50, t50) = Interpolation50(amplitude, timeAxis0, duration, L, CBF);

        // Clip negative
        for (int n = 0; n < amp50.GetLength(0); n++)
        for (int z = 0; z < CBF; z++)
            if (amp50[n, z] < 0.0) amp50[n, z] = 0.0;

        // STEP 10 — Non-linear transform
        double[,] rHat = NonLinearTransform(amp50);

        // STEP 11 — Low-pass filter
        double[,] rTimeSpec = LowPassFilter(rHat);

        // STEP 12 — Representative values
        int n50 = rTimeSpec.GetLength(0);
        double[] rSpec = new double[CBF];
        int startL = Math.Min(10, n50);
        for (int z = 0; z < CBF; z++)
        {
            double sum = 0.0;
            int cnt = 0;
            for (int n = startL; n < n50; n++) { sum += rTimeSpec[n, z]; cnt++; }
            rSpec[z] = cnt > 0 ? sum / cnt : 0.0;
        }

        double[] rTime = new double[n50];
        for (int n = 0; n < n50; n++)
        {
            double s = 0.0;
            for (int z = 0; z < CBF; z++) s += rTimeSpec[n, z];
            rTime[n] = 0.5 * s;
        }

        double R = Percentile90(rTime);
        double[] barkAxis = Interp.Linspace(0.5, 26.5, CBF);

        return (R, rTime, rTimeSpec, barkAxis, t50);
    }

    // -----------------------------------------------------------------------
    // Simple IIR decimation (matches scipy.signal.decimate step)
    // -----------------------------------------------------------------------

    private static double[] DecimateSimple(double[] x, int q) => Decimate.Apply(x, q);

    // -----------------------------------------------------------------------
    // Noise reduction (section 7.1.4)
    // -----------------------------------------------------------------------

    private static double[,,] NoiseReduction(double[,,] spectrum, int L, int cbf, int K)
    {
        double[,,] avg = new double[L, cbf, K];
        for (int l = 0; l < L; l++)
        {
            for (int k = 0; k < K; k++)
                avg[l, 0, k] = (spectrum[l, 0, k] + spectrum[l, 1, k]) / 2.0;
            for (int z = 1; z < cbf - 1; z++)
            for (int k = 0; k < K; k++)
                avg[l, z, k] = (spectrum[l, z - 1, k] + spectrum[l, z, k] + spectrum[l, z + 1, k]) / 3.0;
            for (int k = 0; k < K; k++)
                avg[l, cbf - 1, k] = (spectrum[l, cbf - 1, k] + spectrum[l, cbf - 2, k]) / 2.0;
        }

        double[,] s = new double[L, K];
        for (int l = 0; l < L; l++)
        for (int k = 0; k < K; k++)
        {
            double sum = 0.0;
            for (int z = 0; z < cbf; z++) sum += avg[l, z, k];
            s[l, k] = sum;
        }

        double[] sTilde = new double[L];
        for (int l = 0; l < L; l++)
        {
            int mLen = K - 2;
            double[] vals = new double[mLen];
            for (int k = 2; k < K; k++) vals[k - 2] = s[l, k];
            Array.Sort(vals);
            sTilde[l] = mLen % 2 == 1 ? vals[mLen / 2] : (vals[mLen / 2 - 1] + vals[mLen / 2]) / 2.0;
        }

        double[,] wTilde = new double[L, K];
        for (int l = 0; l < L; l++)
        {
            double st = sTilde[l] + 1e-9;
            for (int k = 0; k < K; k++)
            {
                double clip = Math.Min(0.1891 * Math.Exp(0.0120 * k), 1.0);
                wTilde[l, k] = 0.0856 * (s[l, k] / st) * clip;
            }
        }

        double[] wMax = new double[L];
        for (int l = 0; l < L; l++)
        {
            double m = 0.0;
            for (int k = 2; k < K; k++)
                if (wTilde[l, k] > m) m = wTilde[l, k];
            wMax[l] = m;
        }

        double[,] nsw = new double[L, K];
        for (int l = 0; l < L; l++)
        for (int k = 0; k < K; k++)
        {
            if (wTilde[l, k] >= 0.05 * wMax[l])
                nsw[l, k] = Math.Min(Math.Max(wTilde[l, k] - 0.1407, 0.0), 1.0);
        }

        double[,,] result = new double[L, cbf, K];
        for (int l = 0; l < L; l++)
        for (int z = 0; z < cbf; z++)
        for (int k = 0; k < K; k++)
            result[l, z, k] = avg[l, z, k] * nsw[l, k];

        return result;
    }

    // -----------------------------------------------------------------------
    // Peak picking (section 7.1.5.1)
    // -----------------------------------------------------------------------

    private static (double[] fp, double[] Ai) PeakPicking(double[] phiElz)
    {
        int n = phiElz.Length;

        double[] sub = new double[n - 2];
        for (int k = 2; k < n; k++) sub[k - 2] = phiElz[k];

        var peaks = FindPeaks.Find(sub, prominence: 0.0);
        if (peaks.Indices.Length == 0) return (Array.Empty<double>(), Array.Empty<double>());

        int[] idx    = peaks.Indices;
        double[] prom = peaks.Prominences;
        for (int i = 0; i < idx.Length; i++) idx[i] += 2; // shift back by 2

        // Amplitude condition: keep only where phiElz[idx] > 0.05 * max
        double maxAmp = 0.0;
        for (int i = 0; i < idx.Length; i++)
            if (phiElz[idx[i]] > maxAmp) maxAmp = phiElz[idx[i]];

        // In-place filter using write cursor — eliminates List<int>/List<double> allocations
        int keptCount = 0;
        for (int i = 0; i < idx.Length; i++)
        {
            if (phiElz[idx[i]] > 0.05 * maxAmp)
            {
                idx[keptCount]  = idx[i];
                prom[keptCount] = prom[i];
                keptCount++;
            }
        }

        if (keptCount == 0) return (Array.Empty<double>(), Array.Empty<double>());

        // Keep top 10 by prominence (stable sort preserving original relative order for ties)
        if (keptCount > 10)
        {
            int[] sortOrder = new int[keptCount];
            for (int i = 0; i < keptCount; i++) sortOrder[i] = i;
            Array.Sort(sortOrder, (a, b) =>
            {
                int cmp = prom[b].CompareTo(prom[a]); // descending
                return cmp != 0 ? cmp : a.CompareTo(b); // stable by original index
            });
            int[]    newIdx  = new int[10];
            double[] newProm = new double[10];
            for (int i = 0; i < 10; i++) { newIdx[i] = idx[sortOrder[i]]; newProm[i] = prom[sortOrder[i]]; }
            idx       = newIdx;
            prom      = newProm;
            keptCount = 10;
        }

        int nPeaks = keptCount;
        double deltaF = 1500.0 / 512.0;

        double[] fp = new double[nPeaks];
        double[] Ai = new double[nPeaks];
        for (int i = 0; i < nPeaks; i++)
            (fp[i], Ai[i]) = Refinement(idx[i], phiElz, deltaF);

        return (fp, Ai);
    }

    // -----------------------------------------------------------------------
    // Refinement (section 7.1.5.1)
    // -----------------------------------------------------------------------

    private static (double modRate, double amp) Refinement(int kpi, double[] phiElz, double deltaF)
    {
        int n = phiElz.Length;
        double amp;
        if (kpi == 0)
            amp = phiElz[0] + (kpi + 1 < n ? phiElz[1] : 0.0);
        else if (kpi == n - 1)
            amp = phiElz[kpi - 1] + phiElz[kpi];
        else
            amp = phiElz[kpi - 1] + phiElz[kpi] + phiElz[kpi + 1];

        double f;
        if (kpi == 0 || kpi == n - 1)
            f = kpi * deltaF;
        else
        {
            double denom = 2.0 * phiElz[kpi - 1] + 2.0 * phiElz[kpi + 1] - 4.0 * phiElz[kpi];
            double frac = denom != 0 ?
                (phiElz[kpi + 1] - phiElz[kpi - 1]) / denom : 0.0;
            f = (kpi - frac) * deltaF;
        }

        double modRate = f + Rho(f, deltaF);
        return (modRate, amp);
    }

    private static double Rho(double f, double deltaF)
    {
        double[] E =
        {
            0, 0.0457, 0.0907, 0.1346, 0.1765, 0.2157, 0.2515, 0.2828,
            0.3084, 0.3269, 0.3364, 0.3348, 0.3188, 0.2844, 0.2259, 0.1351,
            0.0000, -0.1351, -0.2259, -0.2844, -0.3188, -0.3348, -0.3364,
            -0.3269, -0.3084, -0.2828, -0.2515, -0.2157, -0.1765, -0.1346,
            -0.0907, -0.0457, 0.000, 0.000,
        };
        int n = E.Length;

        double[] B = new double[n];
        for (int theta = 0; theta < n; theta++)
            B[theta] = (Math.Floor(f / deltaF) + theta / 32.0) * deltaF - (f + E[theta]);

        int thetaMin = 0;
        double minAbs = Math.Abs(B[0]);
        for (int theta = 1; theta < n; theta++)
        {
            if (Math.Abs(B[theta]) < minAbs)
            { minAbs = Math.Abs(B[theta]); thetaMin = theta; }
        }

        int thetaCorr;
        if (thetaMin > 0 && B[thetaMin] * B[thetaMin - 1] < 0)
            thetaCorr = thetaMin;
        else
            thetaCorr = thetaMin + 1;
        if (thetaCorr >= n) thetaCorr = n - 1;

        return E[thetaCorr] -
               ((E[thetaCorr] - E[thetaCorr - 1]) *
                B[thetaCorr] / (B[thetaCorr] - B[thetaCorr - 1] + 1e-30));
    }

    // -----------------------------------------------------------------------
    // Spectral weighting helpers
    // -----------------------------------------------------------------------

    private static double[] BuildCentreFreq()
    {
        const double af = 81.9289, c = 0.1618, step = 0.5;
        double[] f = new double[CBF];
        for (int z = 0; z < CBF; z++)
            f[z] = (af / c) * Math.Sinh(c * (z + 1) * step);
        return f;
    }

    private static double[] ComputeFmax(double[] cf)
    {
        double[] r = new double[CBF];
        for (int z = 0; z < CBF; z++)
            r[z] = 72.6937 * (1.0 - 1.1739 * Math.Exp(-5.4583 * cf[z] / 1000.0));
        return r;
    }

    private static double[] ComputeRmax(double[] cf)
    {
        double[] r = new double[CBF];
        for (int z = 0; z < CBF; z++)
        {
            double r1 = cf[z] < 1000.0 ? 0.3560 : 0.8024;
            double r2 = cf[z] < 1000.0 ? 0.8049 : 0.9333;
            r[z] = 1.0 / (1.0 + r1 * Math.Pow(Math.Abs(Math.Log2(cf[z] / 1000.0)), r2));
        }
        return r;
    }

    private static double[] ComputeQ2High(double[] cf)
    {
        double[] q = new double[CBF];
        for (int z = 0; z < CBF; z++)
        {
            if (cf[z] / 1000.0 < Math.Pow(2, -3.4253))
                q[z] = 0.2471;
            else
                q[z] = 0.2471 + 0.0129 * Math.Pow(Math.Log2(cf[z] / 1000.0) + 3.4253, 2);
        }
        return q;
    }

    private static double[] ComputeQ2Low(double[] cf)
    {
        double[] q = new double[CBF];
        for (int z = 0; z < CBF; z++)
            q[z] = 1.0967 - 0.0640 * Math.Log2(cf[z] / 1000.0);
        return q;
    }

    private static double HighModRateWeighting(double modRate, double amp,
        double fmax, double rmax, double q2High)
    {
        if (modRate < fmax) return amp * rmax;
        double ratio = modRate / fmax - fmax / modRate;
        double G = 1.0 / Math.Pow(1.0 + Math.Pow(ratio * 1.2822, 2), q2High);
        return G * amp * rmax;
    }

    // -----------------------------------------------------------------------
    // Fundamental modulation rate estimation (section 7.1.5.3)
    // Preserves Python quirks:
    //   (1) candidateIdx sorted R order (singles first, then duplicates)
    //   (2) iPeakLocal local-as-global index into fp
    // -----------------------------------------------------------------------

    private static (double modRate, double[] aHat) EstimateFundModRate(
        double[] fp, double[] aiTilde)
    {
        int nPeak = fp.Length;
        double[] E = new double[nPeak];
        int[][] I = new int[nPeak][];

        // Scratch: reusable per-i0 arrays for R values (nPeak ≤ 10)
        double[] R = new double[nPeak];

        for (int i0 = 0; i0 < nPeak; i0++)
        {
            for (int j = 0; j < nPeak; j++)
                R[j] = Math.Round(fp[j] / (fp[i0] + 1e-30));

            // Gather unique R values sorted ascending (matches np.unique)
            // nPeak is tiny (≤10), so O(n²) dedup is fine
            double[] uniqueR = new double[nPeak];
            int nUnique = 0;
            for (int j = 0; j < nPeak; j++)
            {
                bool found = false;
                for (int u = 0; u < nUnique; u++) if (uniqueR[u] == R[j]) { found = true; break; }
                if (!found) uniqueR[nUnique++] = R[j];
            }
            Array.Sort(uniqueR, 0, nUnique);

            // Precompute counts and match lists for each unique R
            int[] counts = new int[nUnique];
            int[][] matchArrays = new int[nUnique][];
            for (int u = 0; u < nUnique; u++)
            {
                int cnt = 0;
                for (int j = 0; j < nPeak; j++) if (R[j] == uniqueR[u]) cnt++;
                counts[u] = cnt;
                int[] matches = new int[cnt];
                int mi = 0;
                for (int j = 0; j < nPeak; j++) if (R[j] == uniqueR[u]) matches[mi++] = j;
                matchArrays[u] = matches;
            }

            // candidateIdx: singles first, then best duplicate representative (Python ordering)
            var candidateIdx = new System.Collections.Generic.List<int>(nPeak);
            for (int u = 0; u < nUnique; u++)
                if (counts[u] == 1) candidateIdx.Add(matchArrays[u][0]);
            for (int u = 0; u < nUnique; u++)
            {
                if (counts[u] > 1)
                {
                    int[] matches = matchArrays[u];
                    int best = matches[0];
                    double bestErr = Math.Abs(fp[matches[0]] / (R[matches[0]] * fp[i0] + 1e-30) - 1);
                    for (int mi = 1; mi < matches.Length; mi++)
                    {
                        double err = Math.Abs(fp[matches[mi]] / (R[matches[mi]] * fp[i0] + 1e-30) - 1);
                        if (err < bestErr) { bestErr = err; best = matches[mi]; }
                    }
                    candidateIdx.Add(best);
                }
            }

            // Harmonic complex: within 4% of exact harmonic
            var hComplex = new System.Collections.Generic.List<int>(candidateIdx.Count);
            foreach (int ci in candidateIdx)
            {
                double relErr = Math.Abs(fp[ci] / (R[ci] * fp[i0] + 1e-30) - 1.0);
                if (relErr < 0.04) hComplex.Add(ci);
            }

            I[i0] = hComplex.ToArray();
            double sumE = 0.0;
            foreach (int ci in I[i0]) sumE += aiTilde[ci];
            E[i0] = sumE;
        }

        int iMax = 0;
        for (int i = 1; i < nPeak; i++) if (E[i] > E[iMax]) iMax = i;

        int[] iSet = I[iMax];
        if (iSet.Length == 0) iSet = new[] { iMax };

        double modRate = fp[iMax];

        // w_peak: centre-of-gravity distance weighting
        // NOTE: replicates Python quirk — i_peak is the LOCAL argmax index within I_max,
        // then f_p[i_peak] uses that local index directly as a GLOBAL index into f_p.
        double sumA = 0.0, cogFNum = 0.0;
        for (int k = 0; k < iSet.Length; k++)
        {
            double a = aiTilde[iSet[k]];
            sumA    += a;
            cogFNum += fp[iSet[k]] * a;
        }
        double cogF = sumA > 0 ? cogFNum / sumA : fp[iMax];

        // Find LOCAL argmax index within iSet (Python quirk: used as GLOBAL index into fp)
        int iPeakLocal = 0;
        for (int k = 1; k < iSet.Length; k++)
            if (aiTilde[iSet[k]] > aiTilde[iSet[iPeakLocal]]) iPeakLocal = k;
        // iPeakLocal as GLOBAL index into fp (matches Python: f_p[i_peak])
        double wPeak = 1.0 + 0.1 * Math.Pow(Math.Abs(cogF - fp[iPeakLocal]), 0.749);

        double[] aHat = new double[iSet.Length];
        for (int k = 0; k < iSet.Length; k++)
            aHat[k] = aiTilde[iSet[k]] * wPeak;

        return (modRate, aHat);
    }

    private static double LowModRateWeighting(double modRate, double[] aHat,
        double fmax, double q2Low)
    {
        double sum = 0.0;
        for (int k = 0; k < aHat.Length; k++) sum += aHat[k];

        if (modRate < fmax)
        {
            double ratio = modRate / fmax - fmax / modRate;
            double G = 1.0 / Math.Pow(1.0 + Math.Pow(ratio * 0.7066, 2), q2Low);
            return G * sum;
        }
        return sum;
    }

    // -----------------------------------------------------------------------
    // Interpolation to 50 Hz (section 7.1.7) — parallel over z
    // -----------------------------------------------------------------------

    private static (double[,] amp50, double[] t50) Interpolation50(
        double[,] amplitude, double[] timeAxis, double duration, int L, int cbf)
    {
        const int Rs50 = 50;
        int n50 = (int)(duration * Rs50);
        double[] t50 = new double[n50];
        for (int n = 0; n < n50; n++) t50[n] = n / (double)Rs50;

        double[,] amp50 = new double[n50, cbf];

        Parallel.For(0, cbf, z =>
        {
            double[] col = new double[L];
            for (int l = 0; l < L; l++) col[l] = amplitude[l, z];
            double[] interped = Interp.Pchip(t50, timeAxis, col);
            for (int n = 0; n < n50; n++) amp50[n, z] = interped[n];
        });

        return (amp50, t50);
    }

    // -----------------------------------------------------------------------
    // Non-linear transform (section 7.1.7)
    // -----------------------------------------------------------------------

    private static double[,] NonLinearTransform(double[,] rEst)
    {
        int n50 = rEst.GetLength(0);
        int cbf  = rEst.GetLength(1);
        const double cR = 0.045;

        double[] rSqMean  = new double[n50];
        double[] rLinMean = new double[n50];
        for (int n = 0; n < n50; n++)
        {
            double sumSq = 0.0, sumLin = 0.0;
            for (int z = 0; z < cbf; z++) { sumSq += rEst[n, z] * rEst[n, z]; sumLin += rEst[n, z]; }
            rSqMean[n]  = Math.Sqrt(sumSq / cbf);
            rLinMean[n] = sumLin / cbf;
        }

        double[] B = new double[n50];
        for (int n = 0; n < n50; n++)
            B[n] = rLinMean[n] != 0 ? rSqMean[n] / rLinMean[n] : 0.0;

        double[] E = new double[n50];
        for (int n = 0; n < n50; n++)
            E[n] = 0.25 * Math.Tanh(1.75 * (B[n] - 2.5)) + 0.7;

        double[,] rHat = new double[n50, cbf];
        for (int n = 0; n < n50; n++)
        for (int z = 0; z < cbf; z++)
            rHat[n, z] = cR * Math.Pow(rEst[n, z], E[n]);

        return rHat;
    }

    // -----------------------------------------------------------------------
    // Low-pass filter (section 7.1.7, Eq. 109–110)
    // Outer n loop stays serial — Python causal-recursion quirk (rHat[1,:] constant input).
    // -----------------------------------------------------------------------

    private static double[,] LowPassFilter(double[,] rHat)
    {
        // Matches _lowpass_filter.py exactly, including Python implementation quirks:
        //   tau[1, :] computed from second time frame is used constant for ALL n >= 1,
        //   and rHat[1, :] is the constant "new input" for all n >= 1.
        //
        // result[0, z] = rHat[0, z]
        // result[n, z] = rHat[1, z] * (1 - alpha_z) + rHat[n-1, z] * alpha_z,  n >= 1

        int n50 = rHat.GetLength(0);
        int cbf = rHat.GetLength(1);
        const double fs50 = 50.0;

        double[,] result = new double[n50, cbf];
        if (n50 == 0) return result;

        for (int z = 0; z < cbf; z++) result[0, z] = rHat[0, z];
        if (n50 < 2) return result;

        // Compute alpha from bark-axis comparison at time frame 1
        double[] alphaZ = new double[cbf];
        alphaZ[0] = Math.Exp(-1.0 / (fs50 * 0.5000));
        for (int z = 1; z < cbf; z++)
        {
            double tau = rHat[1, z] >= rHat[1, z - 1] ? 0.0625 : 0.5000;
            alphaZ[z] = Math.Exp(-1.0 / (fs50 * tau));
        }

        // Precompute rHat[1,z] * (1 - alphaZ[z]) to avoid repeated multiply in the n-loop
        double[] rHat1Scaled = new double[cbf];
        for (int z = 0; z < cbf; z++) rHat1Scaled[z] = rHat[1, z] * (1.0 - alphaZ[z]);

        // Serial outer n (causal Python quirk), vectorised inner z
        for (int n = 1; n < n50; n++)
        for (int z = 0; z < cbf; z++)
            result[n, z] = rHat1Scaled[z] + rHat[n - 1, z] * alphaZ[z];

        return result;
    }

    // -----------------------------------------------------------------------
    // 90th percentile
    // -----------------------------------------------------------------------

    private static double Percentile90(double[] values)
    {
        if (values.Length == 0) return 0.0;
        double[] sorted = (double[])values.Clone();
        Array.Sort(sorted);
        double idx = 0.9 * (sorted.Length - 1);
        int lo = (int)idx;
        int hi = lo + 1;
        if (hi >= sorted.Length) return sorted[lo];
        return sorted[lo] + (idx - lo) * (sorted[hi] - sorted[lo]);
    }
}
