using System.Numerics;
using System.Runtime.CompilerServices;

namespace Mosqito.SqMetrics.Loudness;

/// <summary>
/// Computes the specific loudness pattern from core loudness and integrates it to overall loudness.
/// Exact scalar C# port of MoSQITo <c>_calc_slopes.py</c> (Ernesto Avedillo 2022 version).
///
/// Produces 240 specific loudness values at 0.1 Bark resolution (0.1 to 24 Bark)
/// by attaching upper slopes between critical bands, then integrates.
/// </summary>
internal static class CalcSlopes
{
    private const int SpecLength = 240; // 0.1 to 24.0 Bark in 0.1 steps
    private const int NmWide     = 21;  // nm has 20 core bands + 1 zero padding
    private const int DecComp    = 8;   // decimal places for rounded comparisons

    /// <summary>
    /// Computes specific loudness and overall loudness from core loudness <paramref name="nm"/>.
    /// </summary>
    /// <param name="nm">Core loudness, length 21 (from <see cref="MainLoudness"/>).</param>
    /// <returns>
    /// (<c>N</c> overall loudness [sone];
    ///  <c>NSpecific</c> specific loudness [sone/Bark], length 240).
    /// </returns>
    internal static (double N, double[] NSpecific) Compute(ReadOnlySpan<double> nm)
    {
        double[] rns  = ZwickerData.Rns;
        double[] zup  = ZwickerData.Zup;
        double[,] uslR = ZwickerData.UslReshaped;

        // zupEa[i] = (int)(zup[i] * 10)  (int32 truncation, not round)
        // zupEa[21] = 0  (Python's appended 0, used as zupEa[-1] for i=0)
        Span<int> zupEa = stackalloc int[22]; // zero-initialised by runtime; zupEa[21]=0
        for (int i = 0; i < 21; i++) zupEa[i] = (int)(zup[i] * 10.0);

        double[] nSpec = new double[SpecLength];
        double   N     = 0.0;

        // Initial state: at Bark 0, level 0 (n2_array[-1] = nm[20] = 0)
        double n1 = 0.0;
        double z1 = 0.0;

        // Main loop: iterate over all 21 bands (including nm[20]=0 to terminate the last slope).
        // Python: for i in arange(nm_wide)  →  i = 0..20
        for (int i = 0; i < NmWide; i++)
        {
            // Lower boundary in 240-pt array:
            //   Python: j = zup_ea[i-1]   (i=0 → zup_ea[-1] = 0)
            int jLo = (i == 0) ? zupEa[21] : zupEa[i - 1];  // = 0 for i=0
            int jHi = zupEa[i];

            double nmI  = nm[i];
            double zupI = zup[i];

            // USL table column: Python uses i-1 (wraps to 20 when i=0)
            int col = (i == 0) ? 20 : i - 1;

            // Determine rns and usl based on current n1
            // Python: rns_values_specific[j] = rns[_get_rns_index(n2_array[j-1])]
            //         usl_array_specific[j]  = uslR[rns_idx, i-1]
            int    rnsIdx0 = GetRnsIndex(n1, rns);
            double rns0    = rns[rnsIdx0];
            double usl0    = uslR[rnsIdx0, col];

            bool n1BiggerNm = Math.Round(n1, DecComp) > Math.Round(nmI, DecComp);

            if (!n1BiggerNm)
            {
                // -------------------------------------------------------
                // FLAT: level non-decreasing → rectangular contribution
                // -------------------------------------------------------
                N += nmI * (zupI - z1);
                for (int j = jLo; j < jHi; j++) nSpec[j] = nmI;
                n1 = nmI;
                z1 = zupI;
            }
            else
            {
                // -------------------------------------------------------
                // SLOPE: level decreasing from n1 down toward nmI
                // -------------------------------------------------------
                // Python uses rns_values_specific[j-1] (value from previous position).
                // At the start of the slope, that equals rns[GetRnsIndex(n1)].
                double maxRnsNm = Math.Max(rns0, nmI);
                double z2seg    = Math.Min((n1 - maxRnsNm) / usl0 + z1, zupI);
                double n2seg    = n1 - (z2seg - z1) * usl0;
                N += (z2seg - z1) * (n1 + n2seg) / 2.0;

                // Locals for the inner walk
                double n1L   = n1;
                double z1L   = z1;
                double z2L   = z2seg;
                double n2L   = n2seg;
                double uslL  = usl0;

                // Walk through each 0.1-Bark position in [jLo, jHi)
                // Python: z_array = zup[i-1] + 0.1  (for i=0: zup[-1]+0.1=24.1, never reached since n1=0)
                // Position j corresponds to Bark (j+1)*0.1, so z_arr = (jLo+1)*0.1 initially.
                double zArr = (i == 0 ? zup[20] : zup[i - 1]) + 0.1;

                bool descending = true;

                for (int j = jLo; j < jHi && descending; j++)
                {
                    // Python copies arrays[j] from arrays[j-1] for all j != jLo — handled by locals.

                    bool zGteZ2 = Math.Round(z2L, DecComp) <= Math.Round(zArr, DecComp);

                    if (zGteZ2)
                    {
                        // ---- Hit end of current slope segment ----
                        n1L = n2L;
                        z1L = z2L;

                        int    ri2   = GetRnsIndexEq(n1L, rns);
                        double rnsL  = rns[ri2];
                        double uslL2 = uslR[ri2, col];

                        if (Math.Round(n1L, DecComp) <= Math.Round(nmI, DecComp))
                        {
                            // Reached target — fill rest of band flat at nmI
                            N += nmI * (zupI - z1L);
                            for (int jj = j; jj < jHi; jj++) nSpec[jj] = nmI;
                            n1 = nmI;
                            z1 = zupI;
                            descending = false;
                            break;
                        }
                        else
                        {
                            // Another slope segment with updated usl
                            double maxRns2 = Math.Max(rnsL, nmI);
                            z2L  = Math.Min((n1L - maxRns2) / uslL2 + z1L, zupI);
                            n2L  = n1L - (z2L - z1L) * uslL2;
                            N   += (z2L - z1L) * (n1L + n2L) / 2.0;
                            uslL = uslL2;

                            // NSpecific at this z position (on the new slope)
                            nSpec[j] = Math.Max(0.0, n1L - (zArr - z1L) * uslL);
                        }
                    }
                    else
                    {
                        // z < z2: still on current slope
                        nSpec[j] = Math.Max(0.0, n1L - (zArr - z1L) * uslL);
                    }

                    zArr += 0.1;
                }

                // After inner loop: update state from end of band
                if (descending)
                {
                    n1 = n2L;
                    z1 = z2L;
                }
            }
        }

        // Clamp and round overall loudness
        if (N < 0.0) N = 0.0;
        N = N <= 16.0
            ? Math.Floor(N * 1000.0 + 0.5) / 1000.0
            : Math.Floor(N * 100.0  + 0.5) / 100.0;

        // Clamp specific loudness — SIMD path (AVX2: 4 doubles/iter → 60 iters for 240 pts)
        {
            var zero = new Vector<double>(0.0);
            int vw = Vector<double>.Count;
            int i  = 0;
            for (; i <= SpecLength - vw; i += vw)
                Vector.Max(new Vector<double>(nSpec, i), zero).CopyTo(nSpec, i);
            for (; i < SpecLength; i++)
                if (nSpec[i] < 0.0) nSpec[i] = 0.0;
        }

        return (N, nSpec);
    }

    // ------------------------------------------------------------------
    // Batch overload for time-varying path  (21×nTime → N[nTime], NSpec[240,nTime])
    // ------------------------------------------------------------------

    /// <summary>
    /// Batch variant: computes specific loudness and overall loudness for each time frame.
    /// </summary>
    internal static (double[] N, double[,] NSpecific) ComputeBatch(double[,] nm)
    {
        int nTime = nm.GetLength(1);
        int nmLen = nm.GetLength(0); // 21

        double[]   Nout  = new double[nTime];
        double[,]  NSpec = new double[SpecLength, nTime];
        double[]   col   = new double[nmLen];

        for (int t = 0; t < nTime; t++)
        {
            for (int b = 0; b < nmLen; b++) col[b] = nm[b, t];
            var (n, spec) = Compute(col);
            Nout[t] = n;
            for (int b = 0; b < SpecLength; b++) NSpec[b, t] = spec[b];
        }

        return (Nout, NSpec);
    }

    // ------------------------------------------------------------------
    // RNS index search (matches Python _get_rns_index)
    // ------------------------------------------------------------------

    /// <summary>Returns count of RNS values strictly greater than <paramref name="nm"/>.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal static int GetRnsIndex(double nm, double[] rns)
    {
        // rns[] values are compile-time constants with ≤8 decimal places, so
        // Math.Round(rns[i], 8) == rns[i] exactly — no need to round the table side.
        double nmR = Math.Round(nm, DecComp);
        int cnt = 0;
        for (int i = 0; i < rns.Length; i++)
            if (nmR < rns[i]) cnt++;
        return Math.Min(cnt, 17);
    }

    /// <summary>Returns count of RNS values >= <paramref name="nm"/> (equal_too=True variant).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal static int GetRnsIndexEq(double nm, double[] rns)
    {
        double nmR = Math.Round(nm, DecComp);
        int cnt = 0;
        for (int i = 0; i < rns.Length; i++)
            if (nmR <= rns[i]) cnt++;
        return Math.Min(cnt, 17);
    }
}
