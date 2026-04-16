using System.Threading.Tasks;
using Mosqito.SoundLevelMeter;

namespace Mosqito.SqMetrics.Tonality;

/// <summary>
/// Prominence Ratio (PR) per ECMA-74 Annex D.
/// Ported from MoSQITo <c>pr_ecma_st.py</c> and <c>_pr_main_calc.py</c>.
/// </summary>
public static class PrEcma
{
    // ------------------------------------------------------------------
    // Core calculation (operates on a dB spectrum)
    // ------------------------------------------------------------------

    // Variant that accepts a pre-computed smooth (for joint per-segment processing)
    private static TonalityResult MainCalcWithSmooth(
        ReadOnlySpan<double> specDb, ReadOnlySpan<double> freqAxis,
        ReadOnlySpan<double> filteredSmooth)
    {
        int n = specDb.Length;

        var filtFreqsList  = new List<double>();
        var filtSpecList   = new List<double>();
        var filtSmoothList = new List<double>();
        for (int i = 0; i < n; i++)
        {
            if (freqAxis[i] > 89.1 && freqAxis[i] < 11200.0)
            {
                filtFreqsList.Add(freqAxis[i]);
                filtSpecList.Add(specDb[i]);
                filtSmoothList.Add(filteredSmooth[i]);
            }
        }
        if (filtFreqsList.Count < 3) return new TonalityResult { TTotal = 0 };

        double[] fr     = filtFreqsList.ToArray();
        double[] spec   = filtSpecList.ToArray();
        double[] smooth = filtSmoothList.ToArray();

        List<int> peaks = TonalityInternals.ScreeningForTonesWithSmooth(fr, spec, smooth, 90, 11200);
        int nbTones = peaks.Count;

        var prVals   = new List<double>();
        var promList = new List<bool>();
        var freqList = new List<double>();

        while (nbTones > 0)
        {
            int ind = peaks[0];

            if (peaks.Count > 1)
            {
                int nbT = peaks.Count;
                (ind, _, peaks, nbT) = TonalityInternals.FindHighestTone(fr, spec, peaks, nbT, ind);
                nbTones = nbT;
            }

            double ft = fr[ind];

            var (mf1, mf2) = TonalityInternals.CriticalBand(ft);
            int mLo = TonalityInternals.ArgMinAbs(fr, mf1);
            int mHi = TonalityInternals.ArgMinAbs(fr, mf2);
            double specSumM = TonalityInternals.Sum10Pow(spec.AsSpan(mLo, Math.Max(mHi - mLo, 1)));
            double Lm = specSumM > 0 ? 10.0 * Math.Log10(specSumM) : 0.0;

            var (lf1, lf2) = TonalityInternals.LowerCriticalBand(ft);
            int lLo = TonalityInternals.ArgMinAbs(fr, lf1);
            int lHi = TonalityInternals.ArgMinAbs(fr, lf2);
            double deltaF = lf2 - lf1;
            double specSumL = TonalityInternals.Sum10Pow(spec.AsSpan(lLo, Math.Max(lHi - lLo, 1)));
            double Ll = specSumL > 0 ? 10.0 * Math.Log10(specSumL) : 0.0;

            var (uf1, uf2) = TonalityInternals.UpperCriticalBand(ft);
            int uLo = TonalityInternals.ArgMinAbs(fr, uf1);
            int uHi = TonalityInternals.ArgMinAbs(fr, uf2);
            double specSumU = TonalityInternals.Sum10Pow(spec.AsSpan(uLo, Math.Max(uHi - uLo, 1)));
            double Lu = specSumU > 0 ? 10.0 * Math.Log10(specSumU) : 0.0;

            double delta;
            if (ft <= 171.4)
                delta = 10.0 * Math.Log10(Math.Pow(10.0, 0.1 * Lm))
                      - 10.0 * Math.Log10(((100.0 / deltaF) * Math.Pow(10.0, 0.1 * Ll)
                                           + Math.Pow(10.0, 0.1 * Lu)) * 0.5);
            else
                delta = 10.0 * Math.Log10(Math.Pow(10.0, 0.1 * Lm))
                      - 10.0 * Math.Log10((Math.Pow(10.0, 0.1 * Ll) + Math.Pow(10.0, 0.1 * Lu)) * 0.5);

            if (delta > 0.0)
            {
                freqList.Add(ft);
                prVals.Add(delta);
                if (ft >= 89.1 && ft <= 1000.0)
                    promList.Add(delta >= 9.0 + 10.0 * Math.Log10(1000.0 / ft));
                else if (ft > 1000.0)
                    promList.Add(delta >= 9.0);
                else
                    promList.Add(false);
            }

            peaks.RemoveAll(p => p >= mLo && p <= mHi);
            nbTones = peaks.Count;
        }

        double sum10 = 0.0;
        for (int i = 0; i < prVals.Count; i++)
            if (promList[i]) sum10 += Math.Pow(10.0, prVals[i] / 10.0);
        double tTotal = sum10 > 0 ? 10.0 * Math.Log10(sum10) : 0.0;

        return new TonalityResult
        {
            TTotal = tTotal,
            Values = prVals.ToArray(),
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

        double[] fr   = filtFreqs.ToArray();
        double[] spec = filtSpec.ToArray();

        // Screening
        List<int> peaks = TonalityInternals.ScreeningForTones(fr, spec, 90, 11200);
        int nbTones = peaks.Count;

        var prVals   = new List<double>();
        var promList = new List<bool>();
        var freqList = new List<double>();

        while (nbTones > 0)
        {
            int ind = peaks[0];

            if (peaks.Count > 1)
            {
                int nbT = peaks.Count;
                (ind, _, peaks, nbT) = TonalityInternals.FindHighestTone(fr, spec, peaks, nbT, ind);
                nbTones = nbT;
            }

            double ft = fr[ind];

            // Middle critical band
            var (mf1, mf2) = TonalityInternals.CriticalBand(ft);
            int mLo = TonalityInternals.ArgMinAbs(fr, mf1);
            int mHi = TonalityInternals.ArgMinAbs(fr, mf2);
            double specSumM = TonalityInternals.Sum10Pow(spec.AsSpan(mLo, Math.Max(mHi - mLo, 1)));
            double Lm = specSumM > 0 ? 10.0 * Math.Log10(specSumM) : 0.0;

            // Lower critical band
            var (lf1, lf2) = TonalityInternals.LowerCriticalBand(ft);
            int lLo = TonalityInternals.ArgMinAbs(fr, lf1);
            int lHi = TonalityInternals.ArgMinAbs(fr, lf2);
            double deltaF = lf2 - lf1;
            double specSumL = TonalityInternals.Sum10Pow(spec.AsSpan(lLo, Math.Max(lHi - lLo, 1)));
            double Ll = specSumL > 0 ? 10.0 * Math.Log10(specSumL) : 0.0;

            // Upper critical band
            var (uf1, uf2) = TonalityInternals.UpperCriticalBand(ft);
            int uLo = TonalityInternals.ArgMinAbs(fr, uf1);
            int uHi = TonalityInternals.ArgMinAbs(fr, uf2);
            double specSumU = TonalityInternals.Sum10Pow(spec.AsSpan(uLo, Math.Max(uHi - uLo, 1)));
            double Lu = specSumU > 0 ? 10.0 * Math.Log10(specSumU) : 0.0;

            double delta;
            if (ft <= 171.4)
            {
                delta = 10.0 * Math.Log10(Math.Pow(10.0, 0.1 * Lm))
                      - 10.0 * Math.Log10(((100.0 / deltaF) * Math.Pow(10.0, 0.1 * Ll)
                                           + Math.Pow(10.0, 0.1 * Lu)) * 0.5);
            }
            else
            {
                delta = 10.0 * Math.Log10(Math.Pow(10.0, 0.1 * Lm))
                      - 10.0 * Math.Log10((Math.Pow(10.0, 0.1 * Ll) + Math.Pow(10.0, 0.1 * Lu)) * 0.5);
            }

            if (delta > 0.0)
            {
                freqList.Add(ft);
                prVals.Add(delta);
                if (ft >= 89.1 && ft <= 1000.0)
                    promList.Add(delta >= 9.0 + 10.0 * Math.Log10(1000.0 / ft));
                else if (ft > 1000.0)
                    promList.Add(delta >= 9.0);
                else
                    promList.Add(false);
            }

            // Remove all peaks in middle critical band
            peaks.RemoveAll(p => p >= mLo && p <= mHi);
            nbTones = peaks.Count;
        }

        // T-PR
        double sum10 = 0.0;
        for (int i = 0; i < prVals.Count; i++)
            if (promList[i]) sum10 += Math.Pow(10.0, prVals[i] / 10.0);
        double tTotal = sum10 > 0 ? 10.0 * Math.Log10(sum10) : 0.0;

        return new TonalityResult
        {
            TTotal = tTotal,
            Values = prVals.ToArray(),
            Prominence = promList.ToArray(),
            ToneFrequencies = freqList.ToArray()
        };
    }

    // ------------------------------------------------------------------
    // Public entry points
    // ------------------------------------------------------------------

    /// <summary>
    /// Computes PR for a stationary time signal.
    /// </summary>
    public static TonalityResult ComputeSt(ReadOnlySpan<double> signal, int fs,
        bool prominentOnly = true)
    {
        var (specDb, freqAxis) = CompSpectrum.Compute(signal, fs, db: true);
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
    /// Computes PR from a pre-computed dB spectrum.
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
    /// Computes PR per time segment using joint smooth spectrum (matches Python's
    /// <c>pr_ecma_perseg</c> which processes all segments together via a 2D spectrum,
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

        int nHalf = nSamp / 2;
        double fScale = (double)fs / nSamp;
        double[] freqAxis = new double[nHalf];
        for (int i = 0; i < nHalf; i++) freqAxis[i] = (i + 1) * fScale;

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
