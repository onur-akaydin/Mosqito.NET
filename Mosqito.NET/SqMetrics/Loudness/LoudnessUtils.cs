using Mosqito.Dsp;

namespace Mosqito.SqMetrics.Loudness;

/// <summary>
/// Equal-loudness contours (ISO 226:1987), sone↔phon conversions, and LTQ.
/// Ported from MoSQITo <c>loudness/utils/</c> helper modules.
/// </summary>
public static class LoudnessUtils
{
    // ------------------------------------------------------------------
    // Equal loudness contours — ISO 226:1987 (29 frequencies)
    // Ported from MoSQITo equal_loudness_contours.py (Jeff Tackett MATLAB original)
    // ------------------------------------------------------------------

    private static readonly double[] Freq29 =
    {
        20, 25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200, 250, 315, 400,
        500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300,
        8000, 10000, 12500
    };

    private static readonly double[] Af =
    {
        0.532, 0.506, 0.480, 0.455, 0.432, 0.409, 0.387, 0.367, 0.349,
        0.330, 0.315, 0.301, 0.288, 0.276, 0.267, 0.259, 0.253, 0.250,
        0.246, 0.244, 0.243, 0.243, 0.243, 0.242, 0.242, 0.245, 0.254,
        0.271, 0.301
    };

    private static readonly double[] Lu =
    {
        -31.6, -27.2, -23.0, -19.1, -15.9, -13.0, -10.3, -8.1, -6.2,
        -4.5, -3.1, -2.0, -1.1, -0.4, 0.0, 0.3, 0.5, 0.0, -2.7, -4.1,
        -1.0, 1.7, 2.5, 1.2, -2.1, -7.1, -11.2, -10.7, -3.1
    };

    private static readonly double[] Tf =
    {
        78.5, 68.7, 59.5, 51.1, 44.0, 37.5, 31.5, 26.5, 22.1, 17.9,
        14.4, 11.4, 8.6, 6.2, 4.4, 3.0, 2.2, 2.4, 3.5, 1.7, -1.3,
        -4.2, -6.0, -5.4, -1.5, 6.0, 12.6, 13.9, 12.3
    };

    /// <summary>
    /// Returns a 29-point equal-loudness contour at the given loudness level.
    /// Implements ISO 226:1987 Eq. (1) (sect. 4.1).
    /// </summary>
    /// <param name="phones">Loudness level [phons].</param>
    /// <returns>
    /// (<c>spl</c> — sound pressure levels [dB] at 29 frequencies;
    ///  <c>freqs</c> — the 29 frequency values in Hz).
    /// </returns>
    public static (double[] Spl, double[] Freqs) EqualLoudnessContours(double phones)
    {
        const int n = 29;
        double[] spl = new double[n];

        for (int i = 0; i < n; i++)
        {
            double tfLu = Tf[i] + Lu[i];
            double Af_val = 4.47e-3 * (Math.Pow(10.0, 0.025 * phones) - 1.15)
                          + Math.Pow(0.4 * Math.Pow(10.0, tfLu / 10.0 - 9.0), Af[i]);
            spl[i] = (10.0 / Af[i]) * Math.Log10(Af_val) - Lu[i] + 94.0;
        }

        return (spl, (double[])Freq29.Clone());
    }

    // ------------------------------------------------------------------
    // Sone → Phon (ISO 532-1 / Zwicker & Fastl)
    // ------------------------------------------------------------------

    /// <summary>
    /// Converts loudness in sones to loudness level in phons.
    /// Matches MoSQITo <c>sone_to_phon.py</c>.
    /// </summary>
    public static double SoneToPhon(double sone)
    {
        if (sone < 1.0)
        {
            double phon = 40.0 * Math.Pow(sone, 0.35);
            return Math.Max(phon, 3.0);
        }
        return 10.0 * Math.Log2(sone) + 40.0;
    }

    /// <summary>
    /// Converts loudness in sones to loudness level in phons.
    /// Matches MoSQITo <c>sone2phone.py</c> (Fastl &amp; Zwicker version).
    /// </summary>
    public static double Sone2Phone(double sone)
    {
        if (sone >= 1.0)
            return 40.0 + 10.0 * Math.Log10(sone) / Math.Log10(2.0);
        else
            return 40.0 * Math.Pow(sone + 0.0005, 0.35);
    }

    // ------------------------------------------------------------------
    // Phone → SPL (inverse of equal-loudness contour)
    // ------------------------------------------------------------------

    /// <summary>
    /// Converts phons to approximate dB SPL at a given frequency using ISO 226.
    /// </summary>
    public static double Phone2Spl(double phones, double freq)
    {
        var (spl, freqs) = EqualLoudnessContours(phones);
        return Interp.Linear(freq, freqs, spl);
    }
}
