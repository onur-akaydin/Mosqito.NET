using Mosqito.Io;
using Mosqito.SqMetrics.Loudness;
using Mosqito.Tests.Support;
using Xunit;
using Isoclose = Mosqito.Tests.Support.Isoclose;

namespace Mosqito.Tests.SqMetrics.Loudness;

/// <summary>
/// Ported from MoSQITo tests/sq_metrics/loudness/test_loudness_zwtv.py.
/// Reference: ISO 532-1:2017 Annex B4/B5.
/// </summary>
public class LoudnessZwtvTests
{
    private const string Input = "tests/input";
    private const string TonePulseWav = Input + "/Test signal 10 (tone pulse 1 kHz 10 ms 70 dB).wav";
    private const double WavCalib = 2.0 * 1.4142135623730951; // 2 * sqrt(2)

    // ------------------------------------------------------------------
    // Test 1: Structural correctness (dimensions) for tone pulse signal
    // Mirrors test_loudness_zwtv dimension assertions
    // ------------------------------------------------------------------
    [Fact]
    [Trait("Category", "Loudness")]
    public void TonePulse_OutputDimensions_AreConsistent()
    {
        var (sig, fs) = WavLoader.Read(TonePulseWav, WavCalib);
        var (N, NSpec, barkAxis, timeAxis) = LoudnessZwtv.Compute(sig, fs);

        Assert.Equal(N.Length, timeAxis.Length);
        Assert.Equal(NSpec.GetLength(1), timeAxis.Length);
        Assert.Equal(NSpec.GetLength(0), barkAxis.Length);
        Assert.Equal(240, barkAxis.Length);
    }

    // ------------------------------------------------------------------
    // Test 2: Tone pulse produces a positive loudness peak
    // ------------------------------------------------------------------
    [Fact]
    [Trait("Category", "Loudness")]
    public void TonePulse_PeakLoudness_IsPositive()
    {
        var (sig, fs) = WavLoader.Read(TonePulseWav, WavCalib);
        var (N, _, _, _) = LoudnessZwtv.Compute(sig, fs);

        double maxN = N.Max();
        Assert.True(maxN > 0.1, $"Expected peak N > 0.1 sone, got {maxN:G4}");
    }

    // ------------------------------------------------------------------
    // Test 3: Synthetic steady-state 1 kHz 60 dB sine — N should be near
    // the stationary value (LoudnessZwst gives 4.019 for this signal).
    // We check the time-averaged N over the second half of the signal.
    // ------------------------------------------------------------------
    [Fact]
    [Trait("Category", "Loudness")]
    public void SteadySine_1kHz_60dB_TimeAverageApproachesStationary()
    {
        // Generate 2 s of 1 kHz at 60 dB SPL re 20 µPa
        const int fs = 48000;
        const double d = 2.0;
        const double freq = 1000.0;
        const double spl = 60.0;
        int nSamples = (int)(fs * d);

        // RMS amplitude for 60 dB SPL
        double ampl = 2e-5 * Math.Pow(10.0, spl / 20.0);
        // Peak amplitude = ampl * sqrt(2) for sine
        double peak = ampl * Math.Sqrt(2.0);

        double[] sig = new double[nSamples];
        for (int i = 0; i < nSamples; i++)
            sig[i] = peak * Math.Sin(2.0 * Math.PI * freq * i / fs);

        var (N, _, _, _) = LoudnessZwtv.Compute(sig, fs);

        // Average over second half (steady state)
        int half = N.Length / 2;
        double avgN = 0;
        for (int i = half; i < N.Length; i++) avgN += N[i];
        avgN /= (N.Length - half);

        // Stationary target from ZwstTests is 4.019; time-varying may differ slightly
        // Use generous tolerance (10%, 0.5 sone) for the time-varying path
        Isoclose.Assert(avgN, 4.019, rtol: 0.10, atol: 0.5,
            context: "Zwtv steady-state 1 kHz 60 dB time-average");
    }
}
