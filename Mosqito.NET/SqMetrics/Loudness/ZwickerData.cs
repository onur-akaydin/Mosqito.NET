namespace Mosqito.SqMetrics.Loudness;

/// <summary>
/// Reference data tables for the Zwicker stationary loudness model (ISO 532-1).
/// All tables are verbatim from MoSQITo <c>_main_loudness.py</c> and <c>_calc_slopes.py</c>.
/// </summary>
internal static class ZwickerData
{
    // ------------------------------------------------------------------
    // _main_loudness data tables
    // ------------------------------------------------------------------

    /// <summary>Ranges of 1/3 octave band levels for equal-loudness correction (RAP).</summary>
    internal static readonly int[] Rap = { 45, 55, 65, 71, 80, 90, 100, 120 };

    /// <summary>
    /// Reduction of 1/3 octave band levels at low frequencies (DLL).
    /// Shape [8 ranges, 11 bands].
    /// </summary>
    internal static readonly int[,] Dll =
    {
        { -32, -24, -16, -10, -5,  0, -7, -3, 0, -2, 0 },
        { -29, -22, -15, -10, -4,  0, -7, -2, 0, -2, 0 },
        { -27, -19, -14,  -9, -4,  0, -6, -2, 0, -2, 0 },
        { -25, -17, -12,  -9, -3,  0, -5, -2, 0, -2, 0 },
        { -23, -16, -11,  -7, -3,  0, -4, -1, 0, -1, 0 },
        { -20, -14, -10,  -6, -3,  0, -4, -1, 0, -1, 0 },
        { -18, -12,  -9,  -6, -2,  0, -3, -1, 0, -1, 0 },
        { -15, -10,  -8,  -4, -2,  0, -3, -1, 0, -1, 0 },
    };

    /// <summary>Critical band absolute threshold levels (LTQ) — 20 bands.</summary>
    internal static readonly int[] Ltq = { 30, 18, 12, 8, 7, 6, 5, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 };

    /// <summary>Transmission characteristic correction of the ear (A0) — 20 bands.</summary>
    internal static readonly double[] A0 =
    {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        -0.5, -1.6, -3.2, -5.4, -5.6, -4.0, -1.5, 2.0, 5.0, 12.0
    };

    /// <summary>Free-to-diffuse field level difference (DDF) — 20 bands.</summary>
    internal static readonly double[] Ddf =
    {
        0, 0, 0.5, 0.9, 1.2, 1.6, 2.3, 2.8, 3.0, 2.0,
        0, -1.4, -2.0, -1.9, -1.0, 0.5, 3.0, 4.0, 4.3, 4.0
    };

    /// <summary>Adaptation 1/3 oct → critical band (DCB) — 20 bands.</summary>
    internal static readonly double[] Dcb =
    {
        -0.25, -0.6, -0.8, -0.8, -0.5, 0,
         0.5,   1.1,  1.5,  1.7,  1.8, 1.8,
         1.7,   1.6,  1.4,  1.2,  0.8, 0.5,
         0.0,  -0.5
    };

    // ------------------------------------------------------------------
    // _calc_slopes data tables
    // ------------------------------------------------------------------

    /// <summary>Upper limits of approximated critical bands [Bark], 21 values.</summary>
    internal static readonly double[] Zup =
    {
        0.9, 1.8, 2.8, 3.5, 4.4, 5.4, 6.6, 7.9, 9.2, 10.6,
        12.3, 13.8, 15.2, 16.7, 18.1, 19.3, 20.6, 21.8, 22.7, 23.6, 24.0
    };

    /// <summary>Specific loudness range thresholds for upper slope steepness (RNS), 18 values.</summary>
    internal static readonly double[] Rns =
    {
        21.5, 18.0, 15.1, 11.5, 9.0, 6.1, 4.4, 3.1, 2.13,
        1.36,  0.82, 0.42, 0.30, 0.22, 0.15, 0.10, 0.035, 0.0
    };

    /// <summary>
    /// Steepness of upper slopes (USL) — shape [18, 8].
    /// Extended to [18, 21] by repeating column 7 for columns 8-20 (Ernesto Avedillo 2022).
    /// </summary>
    internal static readonly double[,] Usl =
    {
        { 13.0, 8.2,  6.3,  5.5, 5.5, 5.5, 5.5, 5.5 },
        {  9.0, 7.5,  6.0,  5.1, 4.5, 4.5, 4.5, 4.5 },
        {  7.8, 6.7,  5.6,  4.9, 4.4, 3.9, 3.9, 3.9 },
        {  6.2, 5.4,  4.6,  4.0, 3.5, 3.2, 3.2, 3.2 },
        {  4.5, 3.8,  3.6,  3.2, 2.9, 2.7, 2.7, 2.7 },
        {  3.7, 3.0,  2.8,  2.35,2.2, 2.2, 2.2, 2.2 },
        {  2.9, 2.3,  2.1,  1.9, 1.8, 1.7, 1.7, 1.7 },
        {  2.4, 1.7,  1.5,  1.35,1.3, 1.3, 1.3, 1.3 },
        {  1.95,1.45, 1.3,  1.15,1.1, 1.1, 1.1, 1.1 },
        {  1.5, 1.2,  0.94, 0.86,0.82,0.82,0.82,0.82},
        {  0.72,0.67, 0.64, 0.63,0.62,0.62,0.62,0.62},
        {  0.59,0.53, 0.51, 0.50,0.42,0.42,0.42,0.42},
        {  0.40,0.33, 0.26, 0.24,0.24,0.22,0.22,0.22},
        {  0.27,0.21, 0.20, 0.18,0.17,0.17,0.17,0.17},
        {  0.16,0.15, 0.14, 0.12,0.11,0.11,0.11,0.11},
        {  0.12,0.11, 0.10, 0.08,0.08,0.08,0.08,0.08},
        {  0.09,0.08, 0.07, 0.06,0.06,0.06,0.06,0.05},
        {  0.06,0.05, 0.03, 0.02,0.02,0.02,0.02,0.02},
    };

    // Precomputed USL extended to 21 columns (columns 8-20 = column 7)
    internal static readonly double[,] UslReshaped = BuildUslReshaped();

    private static double[,] BuildUslReshaped()
    {
        int rows = Usl.GetLength(0);   // 18
        int cols = 21;
        double[,] r = new double[rows, cols];
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < 8; j++) r[i, j] = Usl[i, j];
            for (int j = 8; j < cols; j++) r[i, j] = Usl[i, 7];
        }
        return r;
    }
}
