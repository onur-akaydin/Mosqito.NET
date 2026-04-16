using System.Numerics;
using System.Threading.Tasks;
using Mosqito.Dsp;
using Mosqito.SqMetrics.Loudness;

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
    private const int NBark = 240;
    private const double Dz   = 0.1;
    private static readonly int VW = Vector<double>.Count;

    // Bark axis (240 points, 0.1 to 24 Bark)
    private static readonly double[] BarkAxis = Interp.Linspace(0.1, 24.0, NBark);

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

    // Precomputed tables: BarkAxis[i] * Dz * G[i] — dot with NSpecific[*,t] → sum
    // Multiplied by 0.11 at the call site to avoid cluttering the table.
    private static readonly double[] BarkDzGDin      = BuildBarkDzG(SharpnessWeighting.Din);
    private static readonly double[] BarkDzGBismarck = BuildBarkDzG(SharpnessWeighting.Bismarck);
    private static readonly double[] BarkDzGFastl     = BuildBarkDzG(SharpnessWeighting.Fastl);
    // Aures base: BarkAxis[i] * Dz * 0.078 * (exp(0.171*z)/z) — multiply by N/ln(0.05N+1) at runtime
    private static readonly double[] BarkDzAuresBase  = BuildAuresBase();

    private static double[] BuildBarkDzG(SharpnessWeighting w)
    {
        var table = new double[NBark];
        for (int i = 0; i < NBark; i++)
        {
            double z = BarkAxis[i];
            double g = w switch
            {
                SharpnessWeighting.Din =>
                    z > 15.8 ? 0.15 * Math.Exp(0.42 * (z - 15.8)) + 0.85 : 1.0,
                SharpnessWeighting.Bismarck =>
                    z > 15.0 ? 0.2 * Math.Exp(0.308 * (z - 15.0)) + 0.8 : 1.0,
                SharpnessWeighting.Fastl =>
                    Interp.Linear(z, FastlX, FastlY),
                _ => 1.0
            };
            table[i] = z * Dz * g;
        }
        return table;
    }

    private static double[] BuildAuresBase()
    {
        var table = new double[NBark];
        for (int i = 0; i < NBark; i++)
        {
            double z = BarkAxis[i];
            table[i] = z * Dz * 0.078 * (Math.Exp(0.171 * z) / z); // = Dz * 0.078 * exp(0.171*z)
        }
        return table;
    }

    // ------------------------------------------------------------------
    // From loudness (primary kernel)
    // ------------------------------------------------------------------

    /// <summary>
    /// Computes sharpness from Zwicker specific loudness (stationary or time-varying).
    /// </summary>
    public static double FromLoudness(double N, ReadOnlySpan<double> NSpecific,
        SharpnessWeighting weighting = SharpnessWeighting.Din)
    {
        if (N < 0.1) return 0.0;

        double sum;
        if (weighting == SharpnessWeighting.Aures)
        {
            double nFactor = N / Math.Log(N * 0.05 + 1.0);
            sum = nFactor * DotSimd(NSpecific, BarkDzAuresBase);
        }
        else
        {
            double[] table = weighting switch
            {
                SharpnessWeighting.Din      => BarkDzGDin,
                SharpnessWeighting.Bismarck => BarkDzGBismarck,
                SharpnessWeighting.Fastl    => BarkDzGFastl,
                _                           => BarkDzGDin
            };
            sum = DotSimd(NSpecific, table);
        }

        return 0.11 * sum / N;
    }

    /// <summary>
    /// Batch variant: computes sharpness for each time frame — parallel over frames.
    /// </summary>
    public static double[] FromLoudness(ReadOnlySpan<double> N, double[,] NSpecific,
        SharpnessWeighting weighting = SharpnessWeighting.Din)
    {
        int nTime = N.Length;
        int nBark = NSpecific.GetLength(0);
        // Snapshot span to array so it can be captured in the Parallel.For lambda
        double[] NArr = N.ToArray();

        double[] S = new double[nTime];

        // Aures depends on N per frame — can't precompute a single table.
        // All other weightings use a static table — extract column then SIMD dot.
        bool isAures = weighting == SharpnessWeighting.Aures;
        double[] table = isAures ? BarkDzAuresBase : weighting switch
        {
            SharpnessWeighting.Din      => BarkDzGDin,
            SharpnessWeighting.Bismarck => BarkDzGBismarck,
            SharpnessWeighting.Fastl    => BarkDzGFastl,
            _                           => BarkDzGDin
        };

        // Parallel over time frames — each frame is independent.
        // Thread-local column buffer avoids NSpecific[i,t] strided access in SIMD.
        Parallel.For(0, nTime,
            () => new double[nBark],
            (t, _, col) =>
            {
                double nt = NArr[t];
                if (nt < 0.1) { S[t] = 0.0; return col; }

                // Extract column t into contiguous scratch for SIMD
                for (int i = 0; i < nBark; i++) col[i] = NSpecific[i, t];

                double sum;
                if (isAures)
                {
                    double nFactor = nt / Math.Log(nt * 0.05 + 1.0);
                    sum = nFactor * DotSimd(col, table);
                }
                else
                {
                    sum = DotSimd(col, table);
                }

                S[t] = 0.11 * sum / nt;
                return col;
            },
            _ => { });

        return S;
    }

    // ------------------------------------------------------------------
    // Stationary (time-domain entry)
    // ------------------------------------------------------------------

    /// <summary>Computes DIN 45692 sharpness from a time-domain signal (stationary).</summary>
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

    /// <summary>Computes DIN 45692 sharpness from a time-domain signal (time-varying).</summary>
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

    /// <summary>Computes DIN 45692 sharpness per time segment.</summary>
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

    /// <summary>Computes DIN 45692 sharpness from a one-sided amplitude spectrum.</summary>
    public static double ComputeFreq(
        ReadOnlySpan<double> spectrum, ReadOnlySpan<double> freqs,
        SharpnessWeighting weighting = SharpnessWeighting.Din,
        string fieldType = "free")
    {
        var r = LoudnessZwst.ComputeFromSpectrum(spectrum, freqs, fieldType);
        return FromLoudness(r.N, r.NSpecific, weighting);
    }

    // ------------------------------------------------------------------
    // SIMD dot product
    // ------------------------------------------------------------------

    private static double DotSimd(ReadOnlySpan<double> a, ReadOnlySpan<double> b)
    {
        int n    = a.Length;
        int vLen = n - n % VW;
        var acc  = Vector<double>.Zero;
        for (int i = 0; i < vLen; i += VW)
            acc += new Vector<double>(a.Slice(i, VW)) * new Vector<double>(b.Slice(i, VW));
        double sum = Vector.Dot(acc, Vector<double>.One);
        for (int i = vLen; i < n; i++) sum += a[i] * b[i];
        return sum;
    }
}
