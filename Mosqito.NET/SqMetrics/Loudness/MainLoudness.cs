using System.Runtime.CompilerServices;

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
    private const int NRap   = 8;   // number of RAP ranges
    private const int NLow   = 11;  // number of low-frequency bands (DLL applies)
    private const int NBands = 20;  // number of critical band groups

    // ln(10)/10 — converts dB ratios to natural exponent: 10^(x/10) = exp(x * Ln10Over10)
    private const double Ln10Over10 = 0.23025850929940457;

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
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal static double[] Compute(ReadOnlySpan<double> specThird, string fieldType = "free")
    {
        // ISO 532-1 A.3: levels in first 11 bands must not exceed 120 dB.
        for (int i = 0; i < 11 && i < specThird.Length; i++)
            if (specThird[i] > 120.0)
                throw new InvalidOperationException(
                    "1/3 octave band level exceeds 120 dB for which the Zwicker method is no longer valid.");

        // --- Step 1: Equal-loudness correction for bands 0-10 (below 315 Hz) ---
        Span<double> dllResult = stackalloc double[NLow];

        for (int band = 0; band < NLow; band++)
        {
            double spLevel = band < specThird.Length ? specThird[band] : 0.0;
            double dll = ZwickerData.Dll[0, band];
            for (int r = 0; r < NRap - 1; r++)
            {
                double threshold = ZwickerData.Rap[r] - ZwickerData.Dll[r, band];
                if (spLevel > threshold)
                    dll = ZwickerData.Dll[r + 1, band];
                else
                    break;
            }
            {
                double threshold = ZwickerData.Rap[NRap - 1] - ZwickerData.Dll[NRap - 1, band];
                if (spLevel > threshold) dll = 0.0;
            }
            dllResult[band] = dll;
        }

        // xp = dll_result + spec_third, then ti = exp(xp * ln10/10)
        Span<double> xp = stackalloc double[NLow];
        for (int b = 0; b < NLow; b++)
            xp[b] = dllResult[b] + (b < specThird.Length ? specThird[b] : 0.0);

        Span<double> ti = stackalloc double[NLow];
        for (int b = 0; b < NLow; b++)
            ti[b] = Math.Exp(xp[b] * Ln10Over10);

        // --- Step 2: Group first 11 bands into 3 critical bands ---
        Span<double> gi = stackalloc double[3];
        gi[0] = ti[0] + ti[1] + ti[2] + ti[3] + ti[4] + ti[5];
        gi[1] = ti[6] + ti[7] + ti[8];
        gi[2] = ti[9] + ti[10];

        Span<double> lcb = stackalloc double[3];
        lcb.Clear();
        for (int k = 0; k < 3; k++)
            if (gi[k] > 0.0) lcb[k] = 10.0 * Math.Log10(gi[k]);

        // --- Step 3: Compute main loudness for 20 critical bands ---
        Span<double> le = stackalloc double[NBands];
        for (int i = 0; i < NBands; i++)
            le[i] = i + 8 < specThird.Length ? specThird[i + 8] : 0.0;

        le[0] = lcb[0];
        le[1] = lcb[1];
        le[2] = lcb[2];

        double[] a0  = ZwickerData.A0;
        for (int i = 0; i < NBands; i++) le[i] -= a0[i];

        if (fieldType == "diffuse")
        {
            double[] ddf = ZwickerData.Ddf;
            for (int i = 0; i < NBands; i++) le[i] += ddf[i];
        }

        double[]  ltqArr    = ZwickerData.Ltq10Pow; // precomputed: 0.0635*10^(0.025*Ltq[i])
        int[]     ltqInt    = ZwickerData.Ltq;
        double[]  dcb       = ZwickerData.Dcb;
        double[]  nm        = new double[NBands + 1]; // escapes; last element stays 0

        const double s = 0.25;

        for (int i = 0; i < NBands; i++)
        {
            int ltq = ltqInt[i];
            if (le[i] > ltq)
            {
                double leCorr = le[i] - dcb[i];
                double mp1    = ltqArr[i];  // 0.0635 * 10^(0.025*ltq), precomputed
                double inner  = 1.0 - s + s * Math.Exp((leCorr - ltq) * Ln10Over10);
                double mp2    = Math.Sqrt(Math.Sqrt(inner)) - 1.0; // x^0.25 = sqrt(sqrt(x))
                double n      = mp1 * mp2;
                nm[i] = Math.Max(0.0, n);
            }
        }

        // Correction of specific loudness in the lowest critical band
        double korry = 0.4 + 0.32 * Math.Pow(nm[0], 0.2);
        if (korry <= 1.0) nm[0] *= korry;

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
