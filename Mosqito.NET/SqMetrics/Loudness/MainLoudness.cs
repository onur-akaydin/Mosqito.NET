namespace Mosqito.SqMetrics.Loudness;

/// <summary>
/// Computes the core (main) loudness from a 1/3-octave band spectrum.
/// Exact C# port of MoSQITo <c>_main_loudness.py</c> which implements the
/// ISO 532-1:2017 / DIN 45631 algorithm.
///
/// Reference:
///   E. Zwicker and H. Fastl, "Program for calculating loudness according to
///   DIN 45631 (ISO 532-1:2017)", J.A.S.J (E) 12, 1 (1991).
/// </summary>
internal static class MainLoudness
{
    private const int NRap  = 8;   // number of RAP ranges
    private const int NLow  = 11;  // number of low-frequency bands (DLL applies)
    private const int NBands = 20; // number of critical band groups

    /// <summary>
    /// Computes core loudness from a 1/3-octave band spectrum.
    /// </summary>
    /// <param name="specThird">
    /// 1/3-octave band levels [dB ref. 2e-5 Pa].
    /// Length must be ≥ 28 (bands 25 Hz to 12500 Hz).
    /// </param>
    /// <param name="fieldType">"free" (default) or "diffuse".</param>
    /// <returns>
    /// Core loudness array <c>nm</c>, length 21 (20 critical bands + 1 zero padding).
    /// </returns>
    internal static double[] Compute(ReadOnlySpan<double> specThird, string fieldType = "free")
    {
        // ISO 532-1 A.3: levels in first 11 bands must not exceed 120 dB.
        for (int i = 0; i < 11 && i < specThird.Length; i++)
            if (specThird[i] > 120.0)
                throw new InvalidOperationException(
                    "1/3 octave band level exceeds 120 dB for which the Zwicker method is no longer valid.");

        // --- Step 1: Equal-loudness correction for bands 0-10 (below 315 Hz) ---
        // For each of the 11 low-frequency bands, find the dll correction to apply.
        double[] dllResult = new double[NLow];

        for (int band = 0; band < NLow; band++)
        {
            double spLevel = band < specThird.Length ? specThird[band] : 0.0;
            // Find which RAP range this band falls in (first RAP where spec > RAP - dll)
            double dll = ZwickerData.Dll[0, band]; // default to lowest range
            for (int r = 0; r < NRap - 1; r++)
            {
                double threshold = ZwickerData.Rap[r] - ZwickerData.Dll[r, band];
                if (spLevel > threshold)
                    dll = ZwickerData.Dll[r + 1, band];
                else
                    break;
            }
            // If level exceeds highest threshold, dll = 0
            {
                double threshold = ZwickerData.Rap[NRap - 1] - ZwickerData.Dll[NRap - 1, band];
                if (spLevel > threshold) dll = 0.0;
            }
            dllResult[band] = dll;
        }

        // xp = dll_result + spec_third_aux (eq-loudness corrected levels for first 11 bands)
        // Then compute intensities for first 11 bands
        double[] xp = new double[NLow];
        for (int b = 0; b < NLow; b++)
            xp[b] = dllResult[b] + (b < specThird.Length ? specThird[b] : 0.0);

        // ti = 10^(xp/10)
        double[] ti = new double[NLow];
        for (int b = 0; b < NLow; b++)
            ti[b] = Math.Pow(10.0, xp[b] / 10.0);

        // --- Step 2: Group first 11 bands into 3 critical bands (GI → LCB) ---
        double[] gi  = new double[3];
        gi[0] = ti[0] + ti[1] + ti[2] + ti[3] + ti[4] + ti[5];
        gi[1] = ti[6] + ti[7] + ti[8];
        gi[2] = ti[9] + ti[10];

        double[] lcb = new double[3];
        for (int k = 0; k < 3; k++)
            if (gi[k] > 0.0) lcb[k] = 10.0 * Math.Log10(gi[k]);

        // --- Step 3: Compute main loudness for 20 critical bands ---
        // LE = spec_third[8:] (20 values starting from band index 8)
        double[] le = new double[NBands];
        for (int i = 0; i < NBands; i++)
            le[i] = i + 8 < specThird.Length ? specThird[i + 8] : 0.0;

        // Replace first 3 with LCB
        le[0] = lcb[0];
        le[1] = lcb[1];
        le[2] = lcb[2];

        // Apply A0 (ear transmission)
        for (int i = 0; i < NBands; i++)
            le[i] -= ZwickerData.A0[i];

        // Apply DDF for diffuse field
        if (fieldType == "diffuse")
            for (int i = 0; i < NBands; i++)
                le[i] += ZwickerData.Ddf[i];

        // Compute core loudness nm for each critical band
        double[] nm = new double[NBands + 1]; // +1 zero padding (appended in Python)
        const double s = 0.25;

        for (int i = 0; i < NBands; i++)
        {
            double ltq = ZwickerData.Ltq[i];
            if (le[i] > ltq)
            {
                double leCorr = le[i] - ZwickerData.Dcb[i];
                double mp1 = 0.0635 * Math.Pow(10.0, 0.025 * ltq);
                double mp2 = Math.Pow(1.0 - s + s * Math.Pow(10.0, 0.1 * (leCorr - ltq)), 0.25) - 1.0;
                double n = mp1 * mp2;
                nm[i] = Math.Max(0.0, n);
            }
            // else nm[i] remains 0
        }

        // Correction of specific loudness in the lowest critical band
        // korry = 0.4 + 0.32 * nm[0]^0.2
        double korry = 0.4 + 0.32 * Math.Pow(nm[0], 0.2);
        if (korry <= 1.0) nm[0] *= korry;

        // Last element is 0 (padding, matches numpy.append)
        nm[NBands] = 0.0;

        return nm;
    }

    // ------------------------------------------------------------------
    // Batch overload for time-varying path (28×nTime → 21×nTime)
    // ------------------------------------------------------------------

    /// <summary>
    /// Batch variant: computes core loudness for every time frame.
    /// </summary>
    /// <param name="specThird">
    /// Level matrix (28, nTime) in dB re 2e-5 Pa.
    /// Row index = band index (0=25 Hz … 27=12500 Hz).
    /// </param>
    /// <param name="fieldType">"free" or "diffuse".</param>
    /// <returns>Core loudness matrix (21, nTime).</returns>
    internal static double[,] ComputeBatch(double[,] specThird, string fieldType = "free")
    {
        int nBands = specThird.GetLength(0); // 28
        int nTime  = specThird.GetLength(1);
        const int nmLen = NBands + 1;        // 21

        double[,] result = new double[nmLen, nTime];
        double[]  col    = new double[nBands];

        for (int t = 0; t < nTime; t++)
        {
            for (int b = 0; b < nBands; b++) col[b] = specThird[b, t];
            double[] nm = Compute(col, fieldType);
            for (int b = 0; b < nmLen; b++) result[b, t] = nm[b];
        }

        return result;
    }
}
