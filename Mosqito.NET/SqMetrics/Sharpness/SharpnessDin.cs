using Mosqito.Dsp;
using Mosqito.SqMetrics.Loudness;
using Mosqito.SoundLevelMeter;
using Mosqito.Conversion;

namespace Mosqito.SqMetrics.Sharpness;

/// <summary>
/// Sharpness weighting function to use.
/// </summary>
public enum SharpnessWeighting
{
    /// <summary>DIN 45692 standard weighting (default).</summary>
    Din,
    /// <summary>Aures weighting (time- and loudness-dependent).</summary>
    Aures,
    /// <summary>Von Bismarck weighting.</summary>
    Bismarck,
    /// <summary>Zwicker &amp; Fastl weighting (table interpolation).</summary>
    Fastl
}

/// <summary>
/// DIN 45692 acoustic sharpness.
/// Ported from MoSQITo <c>sharpness_din*.py</c>.
///
/// Formula: S = 0.11 × (∑ N'(z) × g(z) × z × Δz) / N   [acum]
/// where z is the Bark axis [0.1, 24] in 0.1 Bark steps.
/// </summary>
public static class SharpnessDin
{
    // Bark axis (240 points, 0.1 to 24 Bark)
    private static readonly double[] BarkAxis = Interp.Linspace(0.1, 24.0, 240);

    // Fastl weighting table (x = Bark position, y = weight)
    private static readonly double[] FastlX =
    {
        0.0985854, 14.826764, 15.364039, 15.863559, 16.297388, 16.533346,
        16.844551, 17.052244, 17.335314, 17.599474, 17.816587, 18.014925,
        18.231972, 18.383130, 18.543644, 18.704224, 18.883648, 18.978205,
        19.167118, 19.346476, 19.441233, 19.554567, 19.668102, 19.762926,
        19.923306, 20.036907, 20.188398, 20.320976, 20.462912, 20.605180,
        20.785070, 20.908228, 21.078962, 21.268276, 21.410543, 21.533966,
        21.704834, 21.847036, 21.989570, 22.141394, 22.312195, 22.464752,
        22.626330, 22.759441, 22.893020, 23.054932, 23.236088, 23.388710,
        23.541600, 23.684800, 23.932978
    };

    private static readonly double[] FastlY =
    {
        0.9783246, 0.9970180, 1.0092967, 1.0169811, 1.0438102, 1.0713409,
        1.0891489, 1.1167798, 1.1441435, 1.1668463, 1.1944437, 1.2268357,
        1.2497056, 1.2775370, 1.3006072, 1.3284053, 1.3561363, 1.3794405,
        1.4118660, 1.4348694, 1.4723566, 1.4908663, 1.5235590, 1.5657737,
        1.5793887, 1.6168091, 1.6682787, 1.7150875, 1.7571354, 1.8228214,
        1.8836461, 1.9304883, 2.0102572, 2.0710485, 2.1367345, 2.2024875,
        2.2917116, 2.3526700, 2.4372666, 2.5123746, 2.5968711, 2.7239833,
        2.8226962, 2.9073262, 3.0250500, 3.1474010, 3.2980514, 3.4298910,
        3.5806417, 3.7125149, 3.9385738
    };

    // ------------------------------------------------------------------
    // From loudness (primary kernel)
    // ------------------------------------------------------------------

    /// <summary>
    /// Computes sharpness from Zwicker specific loudness (stationary or time-varying).
    /// Ported from MoSQITo <c>sharpness_din_from_loudness.py</c>.
    /// </summary>
    /// <param name="N">Overall loudness [sone]. Must be > 0.1 for valid result.</param>
    /// <param name="NSpecific">Specific loudness [sone/Bark], length 240 (one frame).</param>
    /// <param name="weighting">Weighting function.</param>
    /// <returns>Sharpness [acum].</returns>
    public static double FromLoudness(double N, ReadOnlySpan<double> NSpecific,
        SharpnessWeighting weighting = SharpnessWeighting.Din)
    {
        if (N < 0.1) return 0.0;
        const int nBark = 240;
        const double dz = 0.1;

        double sum = 0.0;
        for (int i = 0; i < nBark; i++)
        {
            double z = BarkAxis[i];
            double g = WeightingG(i, z, N, weighting);
            sum += NSpecific[i] * g * z * dz;
        }

        return 0.11 * sum / N;
    }

    /// <summary>
    /// Batch variant: computes sharpness for each time frame.
    /// </summary>
    /// <param name="N">Overall loudness per frame [sone], length nTime.</param>
    /// <param name="NSpecific">Specific loudness [sone/Bark] (240 × nTime).</param>
    /// <param name="weighting">Weighting function.</param>
    /// <returns>Sharpness per frame [acum], length nTime.</returns>
    public static double[] FromLoudness(ReadOnlySpan<double> N, double[,] NSpecific,
        SharpnessWeighting weighting = SharpnessWeighting.Din)
    {
        int nTime = N.Length;
        int nBark = NSpecific.GetLength(0);
        const double dz = 0.1;

        double[] S = new double[nTime];
        for (int t = 0; t < nTime; t++)
        {
            double nt = N[t];
            if (nt < 0.1) { S[t] = 0.0; continue; }

            double sum = 0.0;
            for (int i = 0; i < nBark; i++)
            {
                double z = BarkAxis[i];
                double g = WeightingG(i, z, nt, weighting);
                sum += NSpecific[i, t] * g * z * dz;
            }
            S[t] = 0.11 * sum / nt;
        }
        return S;
    }

    // ------------------------------------------------------------------
    // Stationary (time-domain entry)
    // ------------------------------------------------------------------

    /// <summary>
    /// Computes DIN 45692 sharpness from a time-domain signal (stationary).
    /// Ported from MoSQITo <c>sharpness_din_st.py</c>.
    /// </summary>
    public static double ComputeSt(ReadOnlySpan<double> signal, int fs,
        SharpnessWeighting weighting = SharpnessWeighting.Din,
        string fieldType = "free")
    {
        var r = LoudnessZwst.Compute(signal, fs, fieldType);
        return FromLoudness(r.N, r.NSpecific, weighting);
    }

    // ------------------------------------------------------------------
    // Time-varying entry
    // ------------------------------------------------------------------

    /// <summary>
    /// Computes DIN 45692 sharpness from a time-domain signal (time-varying).
    /// Ported from MoSQITo <c>sharpness_din_tv.py</c>.
    /// </summary>
    public static (double[] S, double[] TimeAxis) ComputeTv(
        ReadOnlySpan<double> signal, int fs,
        SharpnessWeighting weighting = SharpnessWeighting.Din,
        string fieldType = "free")
    {
        var (N, NSpec, _, timeAxis) = LoudnessZwtv.Compute(signal, fs, fieldType);
        double[] S = FromLoudness(N, NSpec, weighting);
        return (S, timeAxis);
    }

    // ------------------------------------------------------------------
    // Per-segment entry
    // ------------------------------------------------------------------

    /// <summary>
    /// Computes DIN 45692 sharpness per time segment.
    /// Ported from MoSQITo <c>sharpness_din_perseg.py</c>.
    /// </summary>
    public static (double[] S, double[] TimeAxis) ComputePerSeg(
        ReadOnlySpan<double> signal, int fs,
        int nPerSeg = 4096, int? noOverlap = null,
        SharpnessWeighting weighting = SharpnessWeighting.Din,
        string fieldType = "free")
    {
        var (N, NSpec, _, timeAxis) = LoudnessZwst.ComputePerSeg(signal, fs, nPerSeg, noOverlap, fieldType);
        double[] S = FromLoudness(N, NSpec, weighting);
        return (S, timeAxis);
    }

    // ------------------------------------------------------------------
    // Frequency-domain entry
    // ------------------------------------------------------------------

    /// <summary>
    /// Computes DIN 45692 sharpness from a one-sided amplitude spectrum.
    /// Ported from MoSQITo <c>sharpness_din_freq.py</c>.
    /// </summary>
    public static double ComputeFreq(
        ReadOnlySpan<double> spectrum, ReadOnlySpan<double> freqs,
        SharpnessWeighting weighting = SharpnessWeighting.Din,
        string fieldType = "free")
    {
        var r = LoudnessZwst.ComputeFromSpectrum(spectrum, freqs, fieldType);
        return FromLoudness(r.N, r.NSpecific, weighting);
    }

    // ------------------------------------------------------------------
    // Private: weighting function g(z)
    // ------------------------------------------------------------------

    private static double WeightingG(int barkIndex, double z, double N,
        SharpnessWeighting weighting)
    {
        switch (weighting)
        {
            case SharpnessWeighting.Din:
                if (z > 15.8)
                    return 0.15 * Math.Exp(0.42 * (z - 15.8)) + 0.85;
                return 1.0;

            case SharpnessWeighting.Aures:
            {
                double ln = Math.Log(N * 0.05 + 1.0);
                return 0.078 * (Math.Exp(0.171 * z) / z) * (N / ln);
            }

            case SharpnessWeighting.Bismarck:
                if (z > 15.0)
                    return 0.2 * Math.Exp(0.308 * (z - 15.0)) + 0.8;
                return 1.0;

            case SharpnessWeighting.Fastl:
                return Interp.Linear(z, FastlX, FastlY);

            default:
                throw new ArgumentOutOfRangeException(nameof(weighting));
        }
    }
}
