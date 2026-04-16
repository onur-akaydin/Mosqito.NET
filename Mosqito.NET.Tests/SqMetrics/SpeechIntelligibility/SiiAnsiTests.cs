using Mosqito.SqMetrics.SpeechIntelligibility;
using Xunit;

namespace Mosqito.Tests.SqMetrics.SpeechIntelligibility;

/// <summary>
/// Ported from MoSQITo tests/sq_metrics/speech_intelligibility/test_sii_ansi.py.
/// Reference: ANSI S3.5 — SII = 0.504 for given octave-band spectrum.
/// </summary>
public class SiiAnsiTests
{
    // Reference from test_main_sii: octave method, speech=[50,40,40,30,20,0],
    // noise=[70,65,45,25,1,-15] → SII = 0.504 (1% tolerance)
    private const double SII_ref = 0.504;

    private static readonly double[] NoiseBands  = { 70, 65, 45, 25, 1, -15 };
    private static readonly double[] SpeechBands = { 50, 40, 40, 30, 20, 0 };
    private static readonly double[] OctaveCenters = { 250, 500, 1000, 2000, 4000, 8000 };

    // ANSI S3.5 conformance: max(ref*0.999, ref-0.01) ≤ SII ≤ min(ref*1.01, ref+0.01)
    private static bool CheckCompliance(double sii) =>
        sii >= Math.Max(SII_ref * 0.999, SII_ref - 0.01) &&
        sii <= Math.Min(SII_ref * 1.01,  SII_ref + 0.01);

    [Fact]
    [Trait("Category", "SII")]
    public void MainSii_OctaveMethod_MatchesReference_0p504()
    {
        // Direct port of Python test_main_sii: _main_sii("octave", speech, noise, threshold=None)
        // must satisfy ANSI S3.5 compliance against SII_ref = 0.504.
        var (sii, _, _) = SiiAnsi.MainSiiForTest(
            "octave", SpeechBands, NoiseBands, threshold: null);

        Assert.True(
            CheckCompliance(sii),
            $"SII = {sii:G10} fails ANSI S3.5 compliance vs reference {SII_ref} " +
            $"(allowed: [{Math.Max(SII_ref * 0.999, SII_ref - 0.01):G6}, " +
            $"{Math.Min(SII_ref * 1.01, SII_ref + 0.01):G6}])");
    }

    [Fact]
    [Trait("Category", "SII")]
    public void AllMethods_DoNotThrow()
    {
        // Smoke test: all method+level combinations run without exception
        string[] methods = { "critical", "equally_critical", "third_octave", "octave" };
        string[] levels  = { "normal", "raised", "loud", "shout" };

        foreach (var m in methods)
        foreach (var l in levels)
        {
            var (sii, spec, freqs) = SiiAnsi.ComputeFromLevel(60, m, l);
            Assert.InRange(sii, 0.0, 1.0);
            Assert.Equal(freqs.Length, spec.Length);
        }
    }

    [Fact]
    [Trait("Category", "SII")]
    public void LevelEntry_ProducesPositiveSII_For60dBNormal()
    {
        var (sii, _, _) = SiiAnsi.ComputeFromLevel(60.0, "critical", "normal");
        Assert.True(sii > 0.0, $"SII = {sii} not positive");
        Assert.True(sii <= 1.0, $"SII = {sii} exceeds 1");
    }

    [Fact]
    [Trait("Category", "SII")]
    public void FreqEntry_MatchesLevelEntry_ForUniformSpectrum()
    {
        string method = "third_octave";
        string level  = "raised";
        double noiseLevel = 55.0;

        var (siiLevel, _, freqAxis) = SiiAnsi.ComputeFromLevel(noiseLevel, method, level);

        // Build matching uniform spectrum
        int nBands = freqAxis.Length;
        double bandLevel = 10.0 * Math.Log10(Math.Pow(10.0, noiseLevel / 10.0) / nBands);
        double[] noiseBands = new double[nBands];
        for (int i = 0; i < nBands; i++) noiseBands[i] = bandLevel;

        var (siiFreq, _, _) = SiiAnsi.ComputeFromSpectrum(noiseBands, freqAxis, method, level);

        Assert.Equal(siiLevel, siiFreq, 6);
    }
}
