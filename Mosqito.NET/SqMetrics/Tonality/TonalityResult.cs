namespace Mosqito.SqMetrics.Tonality;

/// <summary>
/// Result of a single-frame TNR or PR computation.
/// </summary>
public sealed class TonalityResult
{
    /// <summary>Overall (total) TNR or PR value [dB] (T-TNR / T-PR).</summary>
    public double TTotal { get; init; }

    /// <summary>Per-tone TNR or PR values [dB].</summary>
    public double[] Values { get; init; } = Array.Empty<double>();

    /// <summary>Prominence flag for each tone.</summary>
    public bool[] Prominence { get; init; } = Array.Empty<bool>();

    /// <summary>Frequencies of each detected tone [Hz].</summary>
    public double[] ToneFrequencies { get; init; } = Array.Empty<double>();
}
