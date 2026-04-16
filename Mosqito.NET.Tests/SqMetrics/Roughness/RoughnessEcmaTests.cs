using Mosqito.Io;
using Mosqito.SqMetrics.Roughness;
using Xunit;

namespace Mosqito.Tests.SqMetrics.Roughness;

/// <summary>
/// Ported from MoSQITo tests/sq_metrics/roughness/test_roughness_ecma.py.
/// Reference: ECMA-418-2 roughness — 1 kHz AM tone at 70 Hz mod, 60 dB → 1 asper ± 0.1.
/// </summary>
public class RoughnessEcmaTests
{
    private const int Fs = 48000;

    private static double[] MakeAmTone()
    {
        int n = Fs; // 1 second at 48 kHz
        double fmod = 70.0;
        double[] xmod = new double[n];
        for (int i = 0; i < n; i++)
            xmod[i] = Math.Sin(2.0 * Math.PI * fmod * i / Fs);
        var (yAm, _) = SignalGenerators.AmSine(xmod, fs: Fs, fc: 1000.0, splLevel: 60.0);
        return yAm;
    }

    [Fact]
    [Trait("Category", "Roughness")]
    public void AmTone_1kHz_70Hz_60dB_MatchesReference_1Asper()
    {
        double[] sig = MakeAmTone();
        var (R, _, _, _, _) = RoughnessEcma.Compute(sig, Fs);

        // ECMA-418-2 tolerance: ±0.1 asper around 1 asper reference
        Assert.True(R >= 0.9, $"R = {R:F3} below lower bound 0.9 asper");
        Assert.True(R <= 1.1, $"R = {R:F3} above upper bound 1.1 asper");
    }

    [Fact]
    [Trait("Category", "Roughness")]
    public void OutputDimensions_AreConsistent()
    {
        double[] sig = MakeAmTone();
        var (R, rTime, rSpecific, barkAxis, timeAxis) = RoughnessEcma.Compute(sig, Fs);

        Assert.True(R >= 0.0);
        Assert.Equal(53, barkAxis.Length);
        Assert.Equal(rTime.Length, timeAxis.Length);
        Assert.Equal(rTime.Length, rSpecific.GetLength(0));
        Assert.Equal(53, rSpecific.GetLength(1));
    }

    [Fact]
    [Trait("Category", "Roughness")]
    public void RTime_IsPositive_ForAudibleSignal()
    {
        double[] sig = MakeAmTone();
        var (_, rTime, _, _, _) = RoughnessEcma.Compute(sig, Fs);

        bool anyPositive = rTime.Any(v => v > 0.0);
        Assert.True(anyPositive, "All R_time values are zero — signal not detected as rough.");
    }
}
