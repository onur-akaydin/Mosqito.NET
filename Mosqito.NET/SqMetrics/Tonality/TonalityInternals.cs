namespace Mosqito.SqMetrics.Tonality;

/// <summary>
/// Internal helper methods shared by TNR and PR.
/// Ported from MoSQITo tonality helper modules.
/// </summary>
internal static class TonalityInternals
{
    // ------------------------------------------------------------------
    // Critical band (ECMA-74 Annex D.8)
    // ------------------------------------------------------------------

    internal static (double f1, double f2) CriticalBand(double f0)
    {
        double deltaFc = 25.0 + 75.0 * Math.Pow(1.0 + 1.4 * Math.Pow(f0 / 1000.0, 2.0), 0.69);
        double f1, f2;
        if (f0 < 500.0)
        {
            f1 = f0 - deltaFc / 2.0;
            f2 = f0 + deltaFc / 2.0;
        }
        else
        {
            f1 = -deltaFc / 2.0 + Math.Sqrt(deltaFc * deltaFc + 4.0 * f0 * f0) / 2.0;
            f2 = f1 + deltaFc;
        }
        return (f1, f2);
    }

    internal static (double f1, double f2) LowerCriticalBand(double f0)
    {
        var (f2, _) = CriticalBand(f0);
        double c0, c1, c2;
        if (f0 < 171.4) { c0 = 20; c1 = 0; c2 = 0; }
        else if (f0 <= 1600) { c0 = -149.5; c1 = 1.001; c2 = -6.9e-5; }
        else { c0 = 6.8; c1 = 0.806; c2 = -8.2e-6; }
        double f1 = c0 + c1 * f0 + c2 * f0 * f0;
        return (f1, f2);
    }

    internal static (double f1, double f2) UpperCriticalBand(double f0)
    {
        var (_, f1) = CriticalBand(f0);
        double c0, c1, c2;
        if (f0 <= 1600) { c0 = 149.5; c1 = 1.035; c2 = 7.7e-5; }
        else { c0 = 3.3; c1 = 1.215; c2 = 2.16e-5; }
        double f2 = c0 + c1 * f0 + c2 * f0 * f0;
        return (f1, f2);
    }

    // ------------------------------------------------------------------
    // Lower threshold of hearing  (ECMA-74 Annex D.7.1)
    // ------------------------------------------------------------------

    internal static double[] LTH(ReadOnlySpan<double> freqs)
    {
        double[] lth = new double[freqs.Length];
        for (int i = 0; i < freqs.Length; i++)
        {
            double f = freqs[i];
            double fmean, fstd, a1, a2, a3, a4, a5;
            if (f >= 20 && f < 305)
            {
                fmean = 167.5; fstd = 87.3212;
                a1 = 1.415532; a2 = -2.451068; a3 = 1.498869; a4 = -6.983224; a5 = 8.621226;
            }
            else if (f >= 305 && f < 2230)
            {
                fmean = 1157.5; fstd = 488.582;
                a1 = 0.397994; a2 = -0.891839; a3 = -0.815138; a4 = -1.221319; a5 = -7.600754;
            }
            else if (f >= 2230 && f < 14000)
            {
                fmean = 7250; fstd = 3033.25;
                a1 = 1.584978; a2 = -2.766599; a3 = -6.906191; a4 = 10.138553; a5 = -3.149339;
            }
            else if (f >= 14000)
            {
                fmean = 16990; fstd = 4049;
                a1 = -5.775593; a2 = -9.200034; a3 = 26.591150; a4 = 52.167120; a5 = 15.615520;
            }
            else
            {
                lth[i] = 0; continue;
            }
            double ff = (f - fmean) / fstd;
            lth[i] = a1 * ff * ff * ff * ff + a2 * ff * ff * ff + a3 * ff * ff + a4 * ff + a5;
        }
        return lth;
    }

    // ------------------------------------------------------------------
    // Spectrum smoothing (1/N-octave energy-averaged)
    // ------------------------------------------------------------------

    internal static double[] SpectrumSmoothing(
        ReadOnlySpan<double> freqs, ReadOnlySpan<double> spec, int noct,
        double lowFreq, double highFreq)
    {
        int n = freqs.Length;
        double[] smooth = new double[n];

        // Build 1/noct-octave bands from lowFreq to highFreq.
        // Each band is [fc / r, fc * r] where r = 2^(1/(2*noct)).
        double r = Math.Pow(2.0, 1.0 / (2.0 * noct));

        // List of band limits [f_low, f_high]
        var bands = new List<(double fLow, double fHigh)>();
        double fc = lowFreq;
        while (fc <= highFreq)
        {
            bands.Add((fc / r, fc * r));
            fc *= Math.Pow(2.0, 1.0 / noct);
        }

        // For each band, energy-average the spectrum bins in that band
        double[] bandLevels = new double[bands.Count];
        int[] bandLowIdx  = new int[bands.Count];
        int[] bandHighIdx = new int[bands.Count];

        for (int b = 0; b < bands.Count; b++)
        {
            double sum = 0.0;
            int cnt = 0;
            int loIdx = n - 1, hiIdx = 0;
            for (int i = 0; i < n; i++)
            {
                if (freqs[i] >= bands[b].fLow && freqs[i] <= bands[b].fHigh)
                {
                    sum += Math.Pow(10.0, spec[i] / 10.0);
                    cnt++;
                    if (i < loIdx) loIdx = i;
                    if (i > hiIdx) hiIdx = i;
                }
            }
            bandLevels[b] = cnt > 0 ? 10.0 * Math.Log10(sum / cnt) : -1e10;
            bandLowIdx[b]  = loIdx < n ? loIdx : 0;
            bandHighIdx[b] = hiIdx >= 0 ? hiIdx : 0;
        }

        // Project band levels onto frequency axis
        // fill each freq bin with the level of the band that contains it
        double[] bStartFreq = new double[bands.Count];
        double[] bEndFreq   = new double[bands.Count];
        for (int b = 0; b < bands.Count; b++)
        {
            bStartFreq[b] = bands[b].fLow;
            bEndFreq[b]   = bands[b].fHigh;
        }

        for (int i = 0; i < n; i++)
        {
            // Find the band that contains freqs[i]
            double fi = freqs[i];
            smooth[i] = -200.0; // default: very low
            for (int b = 0; b < bands.Count; b++)
            {
                if (fi >= bStartFreq[b] && fi <= bEndFreq[b])
                {
                    smooth[i] = bandLevels[b];
                    break;
                }
            }
        }

        return smooth;
    }

    // ------------------------------------------------------------------
    // Scrambled smooth spectrum — replicates Python's multi-segment
    // _spectrum_smoothing indexing bug for PerSeg joint processing.
    //
    // Python stacks all nseg segments into spec_db[nseg, nfreqs] then
    // calls _spectrum_smoothing(freqs, spec_db.T, ...) where spec_db.T
    // has shape (nfreqs, nseg). After raveling both arrays, the indexing
    // of freqs_in (segment-major) and spec (freq-major) is inconsistent:
    //   spec.ravel()[j*nfreqs + k]  =  spec_db[(j*nfreqs+k)%nseg, (j*nfreqs+k)//nseg]
    // This makes the smooth for the 200/2000 Hz band-edge tones reference
    // out-of-band frequencies (≈102 Hz), giving near-zero smooth and
    // allowing those tones to clear the 6 dB threshold.
    // ------------------------------------------------------------------

    /// <summary>
    /// Computes the scrambled smooth spectrum matching Python's joint multi-segment
    /// <c>_spectrum_smoothing</c> behavior in <c>_screening_for_tones</c>.
    /// </summary>
    /// <param name="specDb">Filtered spectra [nseg, nfreqs].</param>
    /// <param name="freqs">Filtered frequency axis [nfreqs].</param>
    /// <returns>Smooth spectrum [nseg, nfreqs].</returns>
    internal static double[,] ComputeScrambledSmooth(
        double[,] specDb, ReadOnlySpan<double> freqs,
        int noct = 24, double lowFreq = 90.0, double highFreq = 11200.0)
    {
        int nseg   = specDb.GetLength(0);
        int nfreqs = specDb.GetLength(1);
        double[,] smooth = new double[nseg, nfreqs];
        for (int j = 0; j < nseg; j++)
            for (int k = 0; k < nfreqs; k++)
                smooth[j, k] = -200.0;

        double r  = Math.Pow(2.0, 1.0 / (2.0 * noct));
        double fc = lowFreq;
        while (fc <= highFreq)
        {
            double fLow  = fc / r;
            double fHigh = fc * r;

            // Find filtered bins in this band
            var bandBins = new List<int>();
            for (int k = 0; k < nfreqs; k++)
                if (freqs[k] >= fLow && freqs[k] <= fHigh)
                    bandBins.Add(k);

            if (bandBins.Count > 0)
            {
                for (int j = 0; j < nseg; j++)
                {
                    double sumPow = 0.0;
                    foreach (int k in bandBins)
                    {
                        int idx     = j * nfreqs + k;
                        int scrSeg  = idx % nseg;        // Python: spec_db.T.ravel() seg index
                        int scrFreq = idx / nseg;        // Python: spec_db.T.ravel() freq index
                        // scrFreq is guaranteed < nfreqs (proven: max idx = (nseg-1)*nfreqs + nfreqs-1,
                        // max scrFreq = (nseg*nfreqs-1)/nseg = nfreqs - 1/nseg < nfreqs)
                        sumPow += Math.Pow(10.0, specDb[scrSeg, scrFreq] / 10.0);
                    }
                    double bSmooth = 10.0 * Math.Log10(sumPow / bandBins.Count);
                    foreach (int k in bandBins)
                        smooth[j, k] = bSmooth;
                }
            }

            fc *= Math.Pow(2.0, 1.0 / noct);
        }
        return smooth;
    }

    /// <summary>
    /// Screening for tones using a pre-computed smooth spectrum (skips internal smooth
    /// computation). All other criteria (local max, 6 dB threshold, LTH, tonal width)
    /// are applied as in <see cref="ScreeningForTones"/>.
    /// </summary>
    internal static List<int> ScreeningForTonesWithSmooth(
        ReadOnlySpan<double> freqs, ReadOnlySpan<double> spec,
        ReadOnlySpan<double> smooth,
        double lowFreq = 90.0, double highFreq = 11200.0)
    {
        int n   = freqs.Length;
        double[] lth = LTH(freqs);

        var maxima = new List<int>();
        for (int i = 1; i < n - 1; i++)
            if (spec[i] > spec[i - 1] && spec[i] > spec[i + 1])
                maxima.Add(i);

        var cand = new List<int>();
        foreach (int m in maxima)
            if (spec[m] > smooth[m] + 6.0)
                cand.Add(m);

        var audible = new List<int>();
        foreach (int c in cand)
            if (spec[c] > lth[c] + 10.0)
                audible.Add(c);

        var tones = new List<int>();
        int j = 0;
        while (j < audible.Count)
        {
            int peakIdx = audible[j];
            int lowLim  = peakIdx;
            int highLim = peakIdx;

            int t = peakIdx + 1;
            while (t < n && spec[t] > smooth[t] + 6.0)
            {
                if (spec[t] > spec[peakIdx]) peakIdx = t;
                highLim = t;
                t++;
            }
            t = peakIdx - 1;
            while (t >= 0 && spec[t] > smooth[t] + 6.0)
            {
                if (spec[t] > spec[peakIdx]) peakIdx = t;
                lowLim = t;
                t--;
            }

            var (f1, f2) = CriticalBand(freqs[peakIdx]);
            double cbWidth = f2 - f1;
            double tWidth  = freqs[highLim] - freqs[lowLim];

            if (tWidth < cbWidth) tones.Add(peakIdx);

            while (j < audible.Count && audible[j] <= highLim) j++;
        }
        return tones;
    }

    // ------------------------------------------------------------------
    // Screening for tones (Sottek / Bray method — 'smoothed')
    // ------------------------------------------------------------------

    internal static List<int> ScreeningForTones(
        ReadOnlySpan<double> freqs, ReadOnlySpan<double> spec,
        double lowFreq = 90.0, double highFreq = 11200.0)
    {
        int n = freqs.Length;

        // Smoothed spectrum (1/24 octave)
        double[] smooth = SpectrumSmoothing(freqs, spec, 24, lowFreq, highFreq);
        double[] lth    = LTH(freqs);

        // Criterion 1: local maxima
        var maxima = new List<int>();
        for (int i = 1; i < n - 1; i++)
            if (spec[i] > spec[i - 1] && spec[i] > spec[i + 1])
                maxima.Add(i);

        // Criterion 2: at least 6 dB above smoothed
        var cand = new List<int>();
        foreach (int m in maxima)
            if (spec[m] > smooth[m] + 6.0)
                cand.Add(m);

        // Criterion 3: above threshold of hearing + 10 dB
        var audible = new List<int>();
        foreach (int c in cand)
            if (spec[c] > lth[c] + 10.0)
                audible.Add(c);

        // Tonal width criterion: width < half critical band
        var tones = new List<int>();
        int j = 0;
        while (j < audible.Count)
        {
            int peakIdx = audible[j];
            int lowLim  = peakIdx;
            int highLim = peakIdx;

            // Walk right
            int t = peakIdx + 1;
            while (t < n && spec[t] > smooth[t] + 6.0)
            {
                if (spec[t] > spec[peakIdx]) peakIdx = t;
                highLim = t;
                t++;
            }
            // Walk left
            t = peakIdx - 1;
            while (t >= 0 && spec[t] > smooth[t] + 6.0)
            {
                if (spec[t] > spec[peakIdx]) peakIdx = t;
                lowLim = t;
                t--;
            }

            var (f1, f2) = CriticalBand(freqs[peakIdx]);
            double cbWidth = f2 - f1;
            double tWidth  = freqs[highLim] - freqs[lowLim];

            if (tWidth < cbWidth) tones.Add(peakIdx);

            // Skip all candidates up to highLim
            while (j < audible.Count && audible[j] <= highLim) j++;
        }

        return tones;
    }

    // ------------------------------------------------------------------
    // Find highest tone in critical band
    // ------------------------------------------------------------------

    internal static (int indP, int? indS, List<int> peaks, int nbTones)
        FindHighestTone(ReadOnlySpan<double> freqs, ReadOnlySpan<double> spec,
            List<int> peaks, int nbTones, int ind)
    {
        double f = freqs[ind];
        var (f1, f2) = CriticalBand(f);

        int loIdx = ArgMinAbs(freqs, f1);
        int hiIdx = ArgMinAbs(freqs, f2);

        var inBand = new List<int>();
        foreach (int p in peaks)
            if (p > loIdx && p < hiIdx)
                inBand.Add(p);

        if (inBand.Count > 1)
        {
            // Sort by descending level — capture array to avoid span-in-lambda restriction
            double[] specArr = spec.ToArray();
            inBand.Sort((a, b) => specArr[b].CompareTo(specArr[a]));
            int indP = inBand[0];
            int indS = inBand[1];

            // Remove lower ones (3rd and beyond)
            for (int k = 2; k < inBand.Count; k++)
            {
                peaks.Remove(inBand[k]);
                nbTones--;
            }

            if (indP != ind)
            {
                // Re-screen with new primary
                return FindHighestTone(freqs, spec, peaks, nbTones, indP);
            }

            return (indP, indS, peaks, nbTones);
        }

        return (ind, null, peaks, nbTones);
    }

    // ------------------------------------------------------------------
    // Peak level (correct for spectrum resolution)
    // ------------------------------------------------------------------

    internal static double PeakLevel(ReadOnlySpan<double> freqs,
        ReadOnlySpan<double> spec, int peakIdx)
    {
        int n = spec.Length;
        double li = spec[peakIdx];
        double L  = spec[peakIdx];

        // Walk right
        int t = peakIdx + 1;
        if (t < n)
        {
            double lTemp = li;
            while (lTemp - Math.Abs(spec[t]) > 0)
            {
                if (li - spec[t] < 10.0)
                {
                    lTemp = spec[t];
                    L = 10.0 * Math.Log10(Math.Pow(10.0, L / 10.0) + Math.Pow(10.0, spec[t] / 10.0));
                    t++;
                    if (t >= n) { t = n - 1; lTemp = -1; }
                }
                else lTemp = -1;
            }
        }

        // Walk left
        t = peakIdx - 1;
        if (t >= 0)
        {
            double lTemp = li;
            while (lTemp - Math.Abs(spec[t]) > 0)
            {
                if (li - spec[t] < 10.0)
                {
                    lTemp = spec[t];
                    L = 10.0 * Math.Log10(Math.Pow(10.0, L / 10.0) + Math.Pow(10.0, spec[t] / 10.0));
                    t--;
                    if (t < 0) { t = 0; lTemp = -1; }
                }
                else lTemp = -1;
            }
        }

        return L;
    }

    // ------------------------------------------------------------------
    // ArgMin helpers
    // ------------------------------------------------------------------

    internal static int ArgMinAbs(ReadOnlySpan<double> arr, double target)
    {
        int best = 0;
        double bestDist = Math.Abs(arr[0] - target);
        for (int i = 1; i < arr.Length; i++)
        {
            double d = Math.Abs(arr[i] - target);
            if (d < bestDist) { bestDist = d; best = i; }
        }
        return best;
    }

    internal static double Sum10Pow(ReadOnlySpan<double> dBArr)
    {
        double s = 0;
        for (int i = 0; i < dBArr.Length; i++)
            s += Math.Pow(10.0, dBArr[i] / 10.0);
        return s;
    }
}
