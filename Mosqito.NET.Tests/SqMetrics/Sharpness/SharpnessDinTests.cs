using Mosqito.Dsp;
using Mosqito.Io;
using Mosqito.SoundLevelMeter;
using Mosqito.SqMetrics.Sharpness;
using Mosqito.Tests.Support;
using Xunit;
using Isoclose = Mosqito.Tests.Support.Isoclose;

namespace Mosqito.Tests.SqMetrics.Sharpness;

/// <summary>
/// Ported from MoSQITo tests/sq_metrics/sharpness/test_sharpness_din.py.
/// Reference: DIN 45692:2009-E, Chapter 6 compliance tolerance.
/// </summary>
public class SharpnessDinTests
{
    private const string Input = "tests/input";
    private const string BroadbandWav = Input + "/broadband_570.wav";
    private const double S_din = 2.85; // DIN 45692:2009-E reference value [acum]

    // DIN 45692 compliance tolerance: max(ref*0.95, ref-0.05) ≤ S ≤ min(ref*1.05, ref+0.05)
    private static readonly double Rtol = 0.05;
    private static readonly double Atol = 0.05;

    // ------------------------------------------------------------------
    // Test 1: Stationary (time-domain) with DIN weighting
    // Mirrors test_sharpness_din_st
    // ------------------------------------------------------------------
    [Fact]
    [Trait("Category", "Sharpness")]
    public void DinWeighting_BroadbandWav_MatchesReference_2p85()
    {
        var (sig, fs) = WavLoader.Read(BroadbandWav, wavCalib: 1.0);
        double S = SharpnessDin.ComputeSt(sig, fs, SharpnessWeighting.Din);
        Isoclose.Assert(S, S_din, rtol: Rtol, atol: Atol, context: "Sharpness DIN stationary");
    }

    // ------------------------------------------------------------------
    // Test 2: Per-segment with DIN weighting (all averaged segments)
    // Mirrors test_sharpness_din_perseg
    // ------------------------------------------------------------------
    [Fact]
    [Trait("Category", "Sharpness")]
    public void DinWeighting_PerSeg_AllSegmentsNearReference()
    {
        var (sig, fs) = WavLoader.Read(BroadbandWav, wavCalib: 1.0);
        var (S, _) = SharpnessDin.ComputePerSeg(sig, fs, nPerSeg: 1 << 14, weighting: SharpnessWeighting.Din);
        foreach (double s in S)
            Isoclose.Assert(s, S_din, rtol: Rtol, atol: Atol, context: "Sharpness DIN per-seg");
    }

    // ------------------------------------------------------------------
    // Test 3: Frequency-domain path
    // Mirrors test_sharpness_din_freq
    // ------------------------------------------------------------------
    [Fact]
    [Trait("Category", "Sharpness")]
    public void DinWeighting_FrequencyPath_MatchesReference_2p85()
    {
        var (sig, fs) = WavLoader.Read(BroadbandWav, wavCalib: 1.0);
        var (spectrum, freqs) = CompSpectrum.Compute(sig, fs, window: CompSpectrum.WindowType.Blackman, db: false);
        double S = SharpnessDin.ComputeFreq(spectrum, freqs, SharpnessWeighting.Din);
        Isoclose.Assert(S, S_din, rtol: Rtol, atol: Atol, context: "Sharpness DIN frequency path");
    }

    // ------------------------------------------------------------------
    // Test 4: Other weightings (structural / no-crash tests)
    // Mirrors test_sharpness_din_aures, _bismarck, _fastl
    // ------------------------------------------------------------------
    [Theory]
    [InlineData(SharpnessWeighting.Aures)]
    [InlineData(SharpnessWeighting.Bismarck)]
    [InlineData(SharpnessWeighting.Fastl)]
    [Trait("Category", "Sharpness")]
    public void AllWeightings_ProducePositiveResult(SharpnessWeighting weighting)
    {
        var (sig, fs) = WavLoader.Read(BroadbandWav, wavCalib: 1.0);
        double S = SharpnessDin.ComputeSt(sig, fs, weighting);
        Assert.True(S > 0.0, $"Expected S > 0 for {weighting}, got {S:G4}");
    }
}
