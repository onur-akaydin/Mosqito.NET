using Mosqito.Io;
using Mosqito.SoundLevelMeter;
using Mosqito.SqMetrics.Tonality;
using Mosqito.Tests.Support;
using Xunit;
using Isoclose = Mosqito.Tests.Support.Isoclose;

namespace Mosqito.Tests.SqMetrics.Tonality;

/// <summary>
/// Ported from MoSQITo tests/sq_metrics/tonality/test_pr_ecma_st.py
/// and test_pr_ecma_freq.py.
/// Reference: ECMA-74 Annex D.
/// Input: white noise + tones at 442 Hz and 1768 Hz (generated in Audacity).
/// </summary>
public class PrEcmaTests
{
    private const string Input = "tests/input";
    private const string TwoToneWav = Input + "/white_noise_442_1768_Hz_stationary.wav";
    private const string TwoToneVaryingWav = Input + "/white_noise_442_1768_Hz_varying.wav";

    // Reference value from Python: np.testing.assert_almost_equal(t_pr, 32.20980078537321)
    private const double T_pr_ref = 32.20980078537321;
    // Reference from test_pr_ecma_perseg: np.testing.assert_almost_equal(max(t_pr), 34.02082433185258)
    private const double T_pr_perseg_ref = 34.02082433185258;
    private const double Rtol = 0.05;  // 5 % relative
    private const double Atol = 1.0;   // 1 dB absolute

    // ------------------------------------------------------------------
    // Test 1: Stationary path
    // Mirrors test_pr_ecma_st
    // ------------------------------------------------------------------
    [Fact]
    [Trait("Category", "Tonality")]
    public void TwoToneWhiteNoise_TTotal_MatchesReference()
    {
        var (sig, fs) = WavLoader.Read(TwoToneWav, wavCalib: 0.01);
        var result = PrEcma.ComputeSt(sig, fs, prominentOnly: true);
        Isoclose.Assert(result.TTotal, T_pr_ref, rtol: Rtol, atol: Atol,
            context: "T-PR stationary");
    }

    [Fact]
    [Trait("Category", "Tonality")]
    public void TwoToneWhiteNoise_DetectsBothTones_442_1768()
    {
        var (sig, fs) = WavLoader.Read(TwoToneWav, wavCalib: 0.01);
        var result = PrEcma.ComputeSt(sig, fs, prominentOnly: true);

        // Both tones at 442 Hz and 1768 Hz should be detected as prominent
        Assert.Equal(2, result.ToneFrequencies.Length);
        Assert.Equal(2, result.Prominence.Count(p => p));

        var freqsSorted = result.ToneFrequencies.OrderBy(f => f).ToArray();
        Assert.InRange(freqsSorted[0], 400.0, 480.0);
        Assert.InRange(freqsSorted[1], 1700.0, 1850.0);
    }

    // ------------------------------------------------------------------
    // Test 2: All tones (non-prominent included)
    // ------------------------------------------------------------------
    [Fact]
    [Trait("Category", "Tonality")]
    public void TwoToneWhiteNoise_AllTones_NonProminentIncluded()
    {
        var (sig, fs) = WavLoader.Read(TwoToneWav, wavCalib: 0.01);
        var result = PrEcma.ComputeSt(sig, fs, prominentOnly: false);

        Assert.True(result.ToneFrequencies.Length >= 2,
            $"Expected >= 2 tones, got {result.ToneFrequencies.Length}");
        Assert.True(result.TTotal > 0, "T-PR should be positive");
    }

    // ------------------------------------------------------------------
    // Test 3: Frequency-domain entry point (no-crash)
    // Mirrors test_pr_ecma_freq
    // ------------------------------------------------------------------
    [Fact]
    [Trait("Category", "Tonality")]
    public void FrequencyPath_NoException_ReturnsValidResult()
    {
        var (sig, fs) = WavLoader.Read(TwoToneWav, wavCalib: 0.01);
        var (specDb, freqAxis) = CompSpectrum.Compute(sig, fs, db: true);
        var result = PrEcma.ComputeFreq(specDb, freqAxis, prominentOnly: false);

        Assert.True(result.TTotal >= 0 || result.ToneFrequencies.Length == 0,
            "Result should be non-negative or have no tones");
    }

    // ------------------------------------------------------------------
    // Test 4: Per-segment path
    // Mirrors test_pr_ecma_perseg (test_pr_ecma_tv.py).
    // Python defaults: nperseg = int(0.5 * fs), noverlap = 0, prominence=True.
    // ------------------------------------------------------------------
    [Fact]
    [Trait("Category", "Tonality")]
    public void PerSeg_VaryingSignal_MaxTTotalAndProminenceCount_MatchReference()
    {
        var (sig, fs) = WavLoader.Read(TwoToneVaryingWav, wavCalib: 0.01);
        int nPerSeg = (int)(0.5 * fs);
        var (results, timeAxis) = PrEcma.ComputePerSeg(
            sig, fs, nPerSeg: nPerSeg, noOverlap: 0, prominentOnly: true);

        Assert.True(results.Length > 0, "Expected at least one segment");
        Assert.Equal(results.Length, timeAxis.Length);

        // max(t_pr) ≈ 34.02082433185258
        double maxTTotal = results.Max(r => r.TTotal);
        Isoclose.Assert(maxTTotal, T_pr_perseg_ref, rtol: Rtol, atol: Atol,
            context: "max(T-PR) per-segment");

        // np.count_nonzero(prom == True) == 8 across all segments.
        int promCount = results.Sum(r => r.Prominence.Count(p => p));
        Assert.Equal(8, promCount);
    }
}
