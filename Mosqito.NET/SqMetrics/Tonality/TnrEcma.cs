using System.Threading.Tasks;
using Mosqito.Dsp;
using Mosqito.Io;
using Mosqito.SoundLevelMeter;

namespace Mosqito.SqMetrics.Tonality;

/// <summary>
/// Tone-to-Noise Ratio (TNR) per ECMA-74 Annex D and ECMA TR/108.
/// Ported from MoSQITo <c>tnr_ecma_st.py</c> and <c>_tnr_main_calc.py</c>.
/// </summary>
public static class TnrEcma
{
    // ------------------------------------------------------------------
    // Core calculation (operates on a dB spectrum)
    // ------------------------------------------------------------------

    // Variant used by ComputePerSeg to inject the joint scrambled smooth
    // instead of computing a per-segment smooth internally.
    private static TonalityResult MainCalcWithSmooth(
        ReadOnlySpan<double> specDb, ReadOnlySpan<double> freqAxis,
        ReadOnlySpan<double> filteredSmooth)
    {
        int n = specDb.Length;

        var filtFreqs = new List<double>();
        var filtSpec  = new List<double>();
        var filtSmooth = new List<double>();
        for (int i = 0; i < n; i++)
        {
            if (freqAxis[i] > 89.1 && freqAxis[i] < 11200.0)
            {
                filtFreqs.Add(freqAxis[i]);
                filtSpec.Add(specDb[i]);
                filtSmooth.Add(filteredSmooth[i]);
            }
        }
        if (filtFreqs.Count < 3) return new TonalityResult { TTotal = 0 };

        double[] fr     = filtFreqs.ToArray();
        double[] spec   = filtSpec.ToArray();
        double[] smooth = filtSmooth.ToArray();
        double deltaFreq = fr.Length > 1 ? fr[1] - fr[0] : 1.0;

        List<int> peaks = TonalityInternals.ScreeningForTonesWithSmooth(fr, spec, smooth, 90, 11200);
        int nbTones = peaks.Count;

        var tnrVals  = new List<double>();
        var promList = new List<bool>();
        var freqList = new List<double>();

        while (nbTones > 0)
        {
            int ind = peaks[0];
            int? indS = null;

            if (peaks.Count > 1)
            {
                int nbT = peaks.Count;
                (ind, indS, peaks, nbT) = TonalityInternals.FindHighestTone(fr, spec, peaks, nbT, ind);
                nbTones = nbT;
            }

            double fp = fr[ind];
            double ft = fp;
            double Lt, Ltot;
            double deltaFt;

            if (indS.HasValue)
            {
                int sIdx = indS.Value;
                double fs = fr[sIdx];
                double deltaF = 21.0 * Math.Pow(10.0, 1.2 * Math.Pow(Math.Abs(Math.Log10(fp / 212.0)), 1.8));
                if (Math.Abs(fs - fp) < deltaF)
                {
                    double Lp = TonalityInternals.PeakLevel(fr, spec, ind);
                    double Ls = TonalityInternals.PeakLevel(fr, spec, sIdx);
                    Lt = 10.0 * Math.Log10(Math.Pow(10.0, Lp / 10.0) + Math.Pow(10.0, Ls / 10.0));
                    var (f1, f2) = TonalityInternals.CriticalBand(fp);
                    int lo = TonalityInternals.ArgMinAbs(fr, f1);
                    int hi = TonalityInternals.ArgMinAbs(fr, f2);
                    Ltot = 10.0 * Math.Log10(TonalityInternals.Sum10Pow(spec.AsSpan(lo, hi - lo)));
                    peaks.Remove(sIdx);
                    nbTones--;
                    deltaFt = 2.0 * deltaFreq;
                    var (fc1, fc2) = TonalityInternals.CriticalBand(fp);
                    int loFull  = TonalityInternals.ArgMinAbs(freqAxis, fc1);
                    int hiFull  = TonalityInternals.ArgMinAbs(freqAxis, fc2);
                    double deltaFc   = fc2 - fc1;
                    double deltaFtot = freqAxis[hiFull] - freqAxis[loFull];
                    double Ln = 10.0 * Math.Log10(Math.Pow(10.0, Ltot / 10.0) - Math.Pow(10.0, Lt / 10.0))
                              + 10.0 * Math.Log10(deltaFc / (deltaFtot - deltaFt));
                    double delta = Lt - Ln;
                    AppendTone(delta, ft, fr, tnrVals, promList, freqList);
                }
                else
                {
                    Lt = spec[ind];
                    var (f1, f2) = TonalityInternals.CriticalBand(fp);
                    int lo = TonalityInternals.ArgMinAbs(fr, f1);
                    int hi = TonalityInternals.ArgMinAbs(fr, f2);
                    Ltot = 10.0 * Math.Log10(TonalityInternals.Sum10Pow(spec.AsSpan(lo, hi - lo)));
                    deltaFt = deltaFreq;
                    double deltaFc   = f2 - f1;
                    int loFull  = TonalityInternals.ArgMinAbs(freqAxis, f1);
                    int hiFull  = TonalityInternals.ArgMinAbs(freqAxis, f2);
                    double deltaFtot = freqAxis[hiFull] - freqAxis[loFull];
                    double Ln = 10.0 * Math.Log10(Math.Pow(10.0, Ltot / 10.0) - Math.Pow(10.0, Lt / 10.0))
                              + 10.0 * Math.Log10(deltaFc / (deltaFtot - deltaFt));
                    double delta = Lt - Ln;
                    AppendTone(delta, ft, fr, tnrVals, promList, freqList);
                }
            }
            else
            {
                Lt = TonalityInternals.PeakLevel(fr, spec, ind);
                var (f1, f2) = TonalityInternals.CriticalBand(fp);
                int lo = TonalityInternals.ArgMinAbs(fr, f1);
                int hi = TonalityInternals.ArgMinAbs(fr, f2);
                Ltot = 10.0 * Math.Log10(TonalityInternals.Sum10Pow(spec.AsSpan(lo, Math.Max(hi - lo, 1))));
                deltaFt = deltaFreq;
                double deltaFc   = f2 - f1;
                int loFull  = TonalityInternals.ArgMinAbs(freqAxis, f1);
                int hiFull  = TonalityInternals.ArgMinAbs(freqAxis, f2);
                double deltaFtot = freqAxis[hiFull] - freqAxis[loFull];
                double noise = Math.Pow(10.0, Ltot / 10.0) - Math.Pow(10.0, Lt / 10.0);
                double Ln;
                if (noise > 0 && (deltaFtot - deltaFt) > 0)
                    Ln = 10.0 * Math.Log10(noise) + 10.0 * Math.Log10(deltaFc / (deltaFtot - deltaFt));
                else
                    Ln = -200.0;
                double delta = Lt - Ln;
                AppendTone(delta, ft, fr, tnrVals, promList, freqList);
            }

            peaks.Remove(ind);
            nbTones = peaks.Count;
        }

        double tTotal = 0.0;
        double sum10 = 0.0;
        for (int i = 0; i < tnrVals.Count; i++)
            if (promList[i]) sum10 += Math.Pow(10.0, tnrVals[i] / 10.0);
        if (sum10 > 0) tTotal = 10.0 * Math.Log10(sum10);

        return new TonalityResult
        {
            TTotal = tTotal,
            Values = tnrVals.ToArray(),
            Prominence = promList.ToArray(),
            ToneFrequencies = freqList.ToArray()
        };
    }

    private static TonalityResult MainCalc(ReadOnlySpan<double> specDb, ReadOnlySpan<double> freqAxis)
    {
        int n = specDb.Length;

        // Limit to 89.1–11200 Hz
        var filtFreqs = new List<double>();
        var filtSpec  = new List<double>();
        for (int i = 0; i < n; i++)
        {
            if (freqAxis[i] > 89.1 && freqAxis[i] < 11200.0)
            {
                filtFreqs.Add(freqAxis[i]);
                filtSpec.Add(specDb[i]);
            }
        }
        if (filtFreqs.Count < 3) return new TonalityResult { TTotal = 0 };

        double[] fr  = filtFreqs.ToArray();
        double[] spec = filtSpec.ToArray();
        double deltaFreq = fr.Length > 1 ? fr[1] - fr[0] : 1.0;

        // Screening
        List<int> peaks = TonalityInternals.ScreeningForTones(fr, spec, 90, 11200);
        int nbTones = peaks.Count;

        var tnrVals  = new List<double>();
        var promList = new List<bool>();
        var freqList = new List<double>();

        while (nbTones > 0)
        {
            int ind = peaks[0];
            int? indS = null;

            if (peaks.Count > 1)
            {
                int nbT = peaks.Count;
                (ind, indS, peaks, nbT) = TonalityInternals.FindHighestTone(fr, spec, peaks, nbT, ind);
                nbTones = nbT;
            }

            double fp = fr[ind];
            double ft = fp;
            double Lt, Ltot;
            double deltaFt;

            if (indS.HasValue)
            {
                int sIdx = indS.Value;
                double fs = fr[sIdx];
                // proximity criterion
                double deltaF = 21.0 * Math.Pow(10.0, 1.2 * Math.Pow(Math.Abs(Math.Log10(fp / 212.0)), 1.8));
                if (Math.Abs(fs - fp) < deltaF)
                {
                    // Two tones are close — treat as one
                    double Lp = TonalityInternals.PeakLevel(fr, spec, ind);
                    double Ls = TonalityInternals.PeakLevel(fr, spec, sIdx);
                    Lt = 10.0 * Math.Log10(Math.Pow(10.0, Lp / 10.0) + Math.Pow(10.0, Ls / 10.0));

                    var (f1, f2) = TonalityInternals.CriticalBand(fp);
                    int lo = TonalityInternals.ArgMinAbs(fr, f1);
                    int hi = TonalityInternals.ArgMinAbs(fr, f2);
                    Ltot = 10.0 * Math.Log10(TonalityInternals.Sum10Pow(spec.AsSpan(lo, hi - lo)));

                    // Remove second tone from peaks
                    peaks.Remove(sIdx);
                    nbTones--;
                    deltaFt = 2.0 * deltaFreq;

                    // Compute Ln and TNR
                    var (fc1, fc2) = TonalityInternals.CriticalBand(fp);
                    int loFull  = TonalityInternals.ArgMinAbs(freqAxis, fc1);
                    int hiFull  = TonalityInternals.ArgMinAbs(freqAxis, fc2);
                    double deltaFc   = fc2 - fc1;
                    double deltaFtot = freqAxis[hiFull] - freqAxis[loFull];
                    double Ln = 10.0 * Math.Log10(Math.Pow(10.0, Ltot / 10.0) - Math.Pow(10.0, Lt / 10.0))
                              + 10.0 * Math.Log10(deltaFc / (deltaFtot - deltaFt));
                    double delta = Lt - Ln;
                    AppendTone(delta, ft, fr, tnrVals, promList, freqList);
                }
                else
                {
                    // Two tones are far — treat primary alone
                    Lt = spec[ind];
                    var (f1, f2) = TonalityInternals.CriticalBand(fp);
                    int lo = TonalityInternals.ArgMinAbs(fr, f1);
                    int hi = TonalityInternals.ArgMinAbs(fr, f2);
                    Ltot = 10.0 * Math.Log10(TonalityInternals.Sum10Pow(spec.AsSpan(lo, hi - lo)));
                    deltaFt = deltaFreq;

                    double deltaFc   = f2 - f1;
                    int loFull  = TonalityInternals.ArgMinAbs(freqAxis, f1);
                    int hiFull  = TonalityInternals.ArgMinAbs(freqAxis, f2);
                    double deltaFtot = freqAxis[hiFull] - freqAxis[loFull];
                    double Ln = 10.0 * Math.Log10(Math.Pow(10.0, Ltot / 10.0) - Math.Pow(10.0, Lt / 10.0))
                              + 10.0 * Math.Log10(deltaFc / (deltaFtot - deltaFt));
                    double delta = Lt - Ln;
                    AppendTone(delta, ft, fr, tnrVals, promList, freqList);
                }
            }
            else
            {
                // Single tone in critical band
                Lt = TonalityInternals.PeakLevel(fr, spec, ind);
                var (f1, f2) = TonalityInternals.CriticalBand(fp);
                int lo = TonalityInternals.ArgMinAbs(fr, f1);
                int hi = TonalityInternals.ArgMinAbs(fr, f2);
                Ltot = 10.0 * Math.Log10(TonalityInternals.Sum10Pow(spec.AsSpan(lo, Math.Max(hi - lo, 1))));
                deltaFt = deltaFreq;

                double deltaFc   = f2 - f1;
                int loFull  = TonalityInternals.ArgMinAbs(freqAxis, f1);
                int hiFull  = TonalityInternals.ArgMinAbs(freqAxis, f2);
                double deltaFtot = freqAxis[hiFull] - freqAxis[loFull];
                double noise = Math.Pow(10.0, Ltot / 10.0) - Math.Pow(10.0, Lt / 10.0);
                double Ln;
                if (noise > 0 && (deltaFtot - deltaFt) > 0)
                    Ln = 10.0 * Math.Log10(noise) + 10.0 * Math.Log10(deltaFc / (deltaFtot - deltaFt));
                else
                    Ln = -200.0;
                double delta = Lt - Ln;
                AppendTone(delta, ft, fr, tnrVals, promList, freqList);
            }

            peaks.Remove(ind);
            nbTones = peaks.Count;
        }

        // T-TNR = 10 * log10(sum of 10^(tnr_prominent / 10))
        double tTotal = 0.0;
        double sum10 = 0.0;
        for (int i = 0; i < tnrVals.Count; i++)
            if (promList[i]) sum10 += Math.Pow(10.0, tnrVals[i] / 10.0);
        if (sum10 > 0) tTotal = 10.0 * Math.Log10(sum10);

        return new TonalityResult
        {
            TTotal = tTotal,
            Values = tnrVals.ToArray(),
            Prominence = promList.ToArray(),
            ToneFrequencies = freqList.ToArray()
        };
    }

    private static void AppendTone(double delta, double f,
        ReadOnlySpan<double> fr, List<double> tnr, List<bool> prom, List<double> freqs)
    {
        if (delta <= 0.0) return;
        freqs.Add(f);
        tnr.Add(delta);
        if (f >= 89.1 && f < 1000.0)
            prom.Add(delta >= 8.0 + 8.33 * Math.Log10(1000.0 / f));
        else if (f >= 1000.0 && f <= 11200.0)
            prom.Add(delta >= 8.0);
        else
            prom.Add(false);
    }

    // ------------------------------------------------------------------
    // Public entry points
    // ------------------------------------------------------------------

    /// <summary>
    /// Computes TNR for a stationary time signal.
    /// </summary>
    /// <param name="signal">Time signal [Pa].</param>
    /// <param name="fs">Sampling frequency [Hz].</param>
    /// <param name="prominentOnly">If true, filter to prominent tones only (default true).</param>
    /// <returns><see cref="TonalityResult"/> for the signal.</returns>
    public static TonalityResult ComputeSt(ReadOnlySpan<double> signal, int fs,
        bool prominentOnly = true)
    {
        var (specDb, freqAxis) = CompSpectrum.Compute(signal, fs, db: true);
        var result = MainCalc(specDb, freqAxis);
        if (!prominentOnly) return result;

        // Filter to prominent tones
        var filtVals  = new List<double>();
        var filtProm  = new List<bool>();
        var filtFreqs = new List<double>();
        for (int i = 0; i < result.Prominence.Length; i++)
        {
            if (result.Prominence[i])
            {
                filtVals.Add(result.Values[i]);
                filtProm.Add(true);
                filtFreqs.Add(result.ToneFrequencies[i]);
            }
        }
        return new TonalityResult
        {
            TTotal = result.TTotal,
            Values = filtVals.ToArray(),
            Prominence = filtProm.ToArray(),
            ToneFrequencies = filtFreqs.ToArray()
        };
    }

    /// <summary>
    /// Computes TNR from a pre-computed dB spectrum.
    /// </summary>
    public static TonalityResult ComputeFreq(ReadOnlySpan<double> specDb,
        ReadOnlySpan<double> freqAxis, bool prominentOnly = true)
    {
        var result = MainCalc(specDb, freqAxis);
        if (!prominentOnly) return result;

        var filtVals  = new List<double>();
        var filtProm  = new List<bool>();
        var filtFreqs = new List<double>();
        for (int i = 0; i < result.Prominence.Length; i++)
            if (result.Prominence[i])
            {
                filtVals.Add(result.Values[i]);
                filtProm.Add(true);
                filtFreqs.Add(result.ToneFrequencies[i]);
            }
        return new TonalityResult
        {
            TTotal = result.TTotal,
            Values = filtVals.ToArray(),
            Prominence = filtProm.ToArray(),
            ToneFrequencies = filtFreqs.ToArray()
        };
    }

    /// <summary>
    /// Computes TNR per time segment using joint smooth spectrum (matches Python's
    /// <c>tnr_ecma_perseg</c> which processes all segments together via a 2D spectrum,
    /// accidentally using scrambled frequency references in the smooth computation).
    /// </summary>
    public static (TonalityResult[] Results, double[] TimeAxis)
        ComputePerSeg(ReadOnlySpan<double> signal, int fs,
            int nPerSeg = 4096, int? noOverlap = null, bool prominentOnly = true)
    {
        var (blockArray, timeAxis) = Mosqito.Io.TimeSegmentation.Segment(signal, fs, nPerSeg, noOverlap);
        int nSeg = blockArray.GetLength(1);
        int nSamp = blockArray.GetLength(0);

        if (nSeg == 0) return (Array.Empty<TonalityResult>(), timeAxis);

        // freqAxis depends only on fs and nSamp, not signal content — precompute analytically
        int nHalf = nSamp / 2;
        double fScale = (double)fs / nSamp;
        double[] freqAxis = new double[nHalf];
        for (int i = 0; i < nHalf; i++) freqAxis[i] = (i + 1) * fScale;

        // Compute spectra for all segments (parallel — CompSpectrum is thread-safe)
        var allSpecDb = new double[nSeg][];
        Parallel.For(0, nSeg,
            () => new double[nSamp],
            (s, pls, segBuf) =>
            {
                for (int i = 0; i < nSamp; i++) segBuf[i] = blockArray[i, s];
                var (sp, __) = CompSpectrum.Compute(segBuf, fs, db: true);
                allSpecDb[s] = sp;
                return segBuf;
            },
            _ => { });

        // Build filtered (89.1–11200 Hz) freq axis + spectra [nSeg, nfreqs]
        int nFull = freqAxis!.Length;
        var filtIdxList = new List<int>();
        for (int i = 0; i < nFull; i++)
            if (freqAxis[i] > 89.1 && freqAxis[i] < 11200.0)
                filtIdxList.Add(i);
        int nfreqs = filtIdxList.Count;
        int[] filtIdx = filtIdxList.ToArray();

        double[] filtFreqs = new double[nfreqs];
        for (int k = 0; k < nfreqs; k++)
            filtFreqs[k] = freqAxis[filtIdx[k]];

        double[,] filtSpecAll = new double[nSeg, nfreqs];
        for (int s = 0; s < nSeg; s++)
            for (int k = 0; k < nfreqs; k++)
                filtSpecAll[s, k] = allSpecDb[s][filtIdx[k]];

        // Compute joint scrambled smooth [nSeg, nfreqs]
        double[,] smooth = TonalityInternals.ComputeScrambledSmooth(filtSpecAll, filtFreqs);

        // Per-segment MainCalcWithSmooth — parallel; each thread owns its fullSmooth buffer
        var results = new TonalityResult[nSeg];
        Parallel.For(0, nSeg,
            () => new double[nFull],
            (s, _, fullSmooth) =>
            {
                Array.Fill(fullSmooth, -200.0);
                for (int k = 0; k < nfreqs; k++)
                    fullSmooth[filtIdx[k]] = smooth[s, k];

                var result = MainCalcWithSmooth(allSpecDb[s], freqAxis, fullSmooth);
                if (!prominentOnly) { results[s] = result; return fullSmooth; }

                var filtVals     = new List<double>();
                var filtProm     = new List<bool>();
                var filtFreqsOut = new List<double>();
                for (int i = 0; i < result.Prominence.Length; i++)
                    if (result.Prominence[i])
                    {
                        filtVals.Add(result.Values[i]);
                        filtProm.Add(true);
                        filtFreqsOut.Add(result.ToneFrequencies[i]);
                    }
                results[s] = new TonalityResult
                {
                    TTotal = result.TTotal,
                    Values = filtVals.ToArray(),
                    Prominence = filtProm.ToArray(),
                    ToneFrequencies = filtFreqsOut.ToArray()
                };
                return fullSmooth;
            },
            _ => { });

        return (results, timeAxis);
    }
}
