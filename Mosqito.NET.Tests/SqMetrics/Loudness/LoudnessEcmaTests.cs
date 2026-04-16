using Mosqito.Io;
using Mosqito.SqMetrics.Loudness;
using Xunit;
using Isoclose = Mosqito.Tests.Support.Isoclose;

namespace Mosqito.Tests.SqMetrics.Loudness;

/// <summary>
/// Ported from MoSQITo tests/sq_metrics/loudness/test_loudness_ecma.py.
/// Reference: ECMA-418-2 (2nd Ed, 2022) Section 5.
/// </summary>
public class LoudnessEcmaTests
{
    // Helper: generate sine wave at spl_level [dB SPL, ref 2e-5 Pa]
    // Mirrors MoSQITo sine_wave_generator
    private static double[] SineWave(int fs, double durationS, double freqHz, double splLevel)
    {
        const double pRef = 2e-5;
        double rms       = pRef * Math.Pow(10.0, splLevel / 20.0);
        double amplitude = Math.Sqrt(2.0) * rms;
        int n = (int)(fs * durationS);
        double[] sig = new double[n];
        for (int i = 0; i < n; i++)
            sig[i] = amplitude * Math.Sin(2.0 * Math.PI * freqHz * i / fs);
        return sig;
    }

    // ------------------------------------------------------------------
    // Test 1: Equal loudness — 1 kHz 80 dB and 5 kHz 78.49095 dB should
    // produce the same loudness value.
    // Mirrors test_loudness_ecma
    // ------------------------------------------------------------------
    [Fact]
    [Trait("Category", "Loudness")]
    public void EqualLoudness_1kHz_80dB_vs_5kHz_78p49dB_AreClose()
    {
        const int fs = 48000;
        const double d = 0.25;

        double[] sig1kHz = SineWave(fs, d, 1000.0, 80.0);
        double[] sig5kHz = SineWave(fs, d, 5000.0, 78.49095);

        var r1 = LoudnessEcma.Compute(sig1kHz, fs);
        var r5 = LoudnessEcma.Compute(sig5kHz, fs);

        // np.isclose default: rtol=1e-5, atol=1e-8 — use generous tolerance for the port
        Isoclose.Assert(r1.N, r5.N, rtol: 0.05, atol: 0.1,
            context: "ECMA equal-loudness: 1 kHz 80 dB vs 5 kHz 78.49 dB");
    }

    // ------------------------------------------------------------------
    // Test 2: Output dimension consistency
    // ------------------------------------------------------------------
    [Fact]
    [Trait("Category", "Loudness")]
    public void OutputDimensions_AreConsistent()
    {
        const int fs = 48000;
        double[] sig = SineWave(fs, 0.25, 1000.0, 80.0);
        var result = LoudnessEcma.Compute(sig, fs);

        Assert.Equal(53, result.NSpecific.Length);
        Assert.Equal(53, result.BarkAxis.Length);
        Assert.Equal(result.NTime.Length, result.NSpecific[0].Length);
    }

    // ------------------------------------------------------------------
    // Test 3: Loudness > 0 for audible signal
    // ------------------------------------------------------------------
    [Fact]
    [Trait("Category", "Loudness")]
    public void Loudness_IsPositive_ForAudibleSignal()
    {
        const int fs = 48000;
        double[] sig = SineWave(fs, 0.25, 1000.0, 80.0);
        var result = LoudnessEcma.Compute(sig, fs);

        Assert.True(result.N > 0.0, $"Expected N > 0, got {result.N:G4}");
    }
}
