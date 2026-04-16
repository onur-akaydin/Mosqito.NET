using Mosqito.Dsp;

namespace Mosqito.SqMetrics.Loudness;

/// <summary>
/// Threshold of quiet (level threshold in quiet) as a function of Bark frequency.
/// Ported from MoSQITo <c>utils/LTQ.py</c>.
/// </summary>
public static class LTQ
{
    // Zwicker/Fastl Psychoacoustics fig. 2.1 — free-field threshold re 2e-5 Pa.
    private static readonly double[] ZwickerX =
    {
        2.40445500e-04, 2.97265560e-04, 3.83444580e-04, 4.87665300e-04,
        6.38024300e-04, 8.26911900e-04, 1.06166650e-03, 1.54817010e-03,
        2.08369200e-03, 2.92606650e-03, 4.05138200e-03, 6.02053940e-03,
        9.07457100e-03, 1.38729700e-02, 2.19208700e-02, 3.27319300e-02,
        5.14789440e-02, 7.68705100e-02, 1.13174880e-01, 1.81393830e-01,
        3.22530400e-01, 4.97813570e-01, 7.79302600e-01, 1.26089480e+00,
        1.83900010e+00, 2.49687100e+00, 3.30065520e+00, 4.01741500e+00,
        4.82492625e+00, 5.38932800e+00, 6.06309500e+00, 6.92502964e+00,
        8.12003875e+00, 9.61927600e+00, 1.15144812e+01, 1.26382940e+01,
        1.41298833e+01, 1.53680458e+01, 1.68414347e+01, 1.86183590e+01,
        1.99594350e+01, 2.10583273e+01, 2.17400030e+01, 2.22243315e+01,
        2.25462820e+01, 2.27627940e+01, 2.29925427e+01, 2.31538743e+01,
        2.32710993e+01, 2.33580350e+01, 2.34824357e+01,
    };

    private static readonly double[] ZwickerY =
    {
        73.28456, 69.49444, 65.17124, 61.262524, 57.471996, 53.918324,
        50.60151, 46.27743, 42.78269, 39.050854, 36.08868, 32.060455,
        28.624102, 25.48363, 22.461315, 20.090572, 18.015444, 16.47346,
        15.286756, 14.158645, 12.970594, 12.612313, 12.07634, 11.066555,
        10.412693, 9.108175, 7.8235803, 6.4902444, 5.551555, 5.00799,
        4.513699, 5.2524257, 6.6815033, 8.505031, 11.117457, 12.546872,
        13.827873, 14.665177, 15.699439, 17.867777, 20.875132, 24.474821,
        28.814922, 33.59928, 38.82787, 45.191124, 51.653038, 58.213753,
        65.02121, 75.87385, 89.63695,
    };

    // Threshold used in roughness calculation.
    private static readonly double[] RoughnessX =
    {
        0, 0.01, 0.17, 0.8, 1, 1.5, 2, 3.3, 4, 5, 6, 8, 10, 12,
        13.3, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 24.5, 25,
    };

    private static readonly double[] RoughnessY =
    {
        130, 70, 60, 30, 25, 20, 15, 10, 8.1, 6.3, 5, 3.5, 2.5, 1.7,
        0, -2.5, -4, -3.7, -1.5, 1.4, 3.8, 5, 7.5, 15, 48, 60, 130,
    };

    /// <summary>
    /// Computes the threshold in quiet (dB re 2e-5 Pa) at each Bark frequency.
    /// </summary>
    /// <param name="barkAxis">Frequency axis in Bark.</param>
    /// <param name="reference">"zwicker" (default) or "roughness".</param>
    public static double[] Compute(ReadOnlySpan<double> barkAxis, string reference = "zwicker")
    {
        double[] xp, yp;
        if (reference == "roughness")
        {
            xp = RoughnessX;
            yp = RoughnessY;
        }
        else
        {
            xp = ZwickerX;
            yp = ZwickerY;
        }

        double[] result = new double[barkAxis.Length];
        for (int i = 0; i < barkAxis.Length; i++)
            result[i] = Interp.Linear(barkAxis[i], xp, yp);
        return result;
    }
}
