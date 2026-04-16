using Mosqito.Io;
using Mosqito.SqMetrics.Roughness;
using Xunit;

namespace Mosqito.Tests.SqMetrics.Roughness;

/// <summary>
/// Ported from MoSQITo tests/sq_metrics/roughness/test_roughness_dw.py.
/// Reference: Daniel and Weber 1997 — 1 kHz AM tone at 70 Hz mod, 60 dB → 1 asper ± 17%.
/// </summary>
public class RoughnessDwTests
{
    // Daniel and Weber reference: 1 kHz 60 dB tone AM at 70 Hz → 1 asper
    private const double R_ref = 1.0;
    private const double Rtol  = 0.17; // 17% tolerance

    private static (double[] signal, int fs) MakeAmTone(int fs = 44100)
    {
        int n = fs; // 1 second
        double fmod = 70.0;
        double[] time   = new double[n];
        double[] xmod   = new double[n];
        for (int i = 0; i < n; i++)
        {
            time[i] = i / (double)fs;
            xmod[i] = Math.Sin(2.0 * Math.PI * fmod * time[i]);
        }
        var (yAm, _) = SignalGenerators.AmSine(xmod, fs: fs, fc: 1000.0, splLevel: 60.0);
        return (yAm, fs);
    }

    [Fact]
    [Trait("Category", "Roughness")]
    public void AmTone_TimeDomain_IsWithin17Pct_of_1Asper()
    {
        var (signal, fs) = MakeAmTone(44100);
        var result = RoughnessDw.Compute(signal, fs, overlap: 0);

        foreach (double r in result.R)
        {
            Assert.True(r >= R_ref * (1.0 - Rtol),
                $"R = {r:F3} below lower bound {R_ref * (1.0 - Rtol):F3}");
            Assert.True(r <= R_ref * (1.0 + Rtol),
                $"R = {r:F3} above upper bound {R_ref * (1.0 + Rtol):F3}");
        }
    }

    [Fact]
    [Trait("Category", "Roughness")]
    public void AmTone_FreqDomain_IsWithin17Pct_of_1Asper()
    {
        int fs = 48000;
        var (signal, _) = MakeAmTone(fs);
        int nperseg = (int)(0.2 * fs);

        // Build Blackman-windowed spectrum matching Python comp_spectrum(db=False)
        var (spec, freqs) = Mosqito.SoundLevelMeter.CompSpectrum.Compute(
            signal.AsSpan(0, nperseg), fs,
            nfft: nperseg,
            window: Mosqito.SoundLevelMeter.CompSpectrum.WindowType.Blackman,
            oneSided: true, db: false);

        var (R, _, _) = RoughnessDw.ComputeFromSpectrum(spec, freqs);

        Assert.True(R >= R_ref * (1.0 - Rtol),
            $"R = {R:F3} below lower bound {R_ref * (1.0 - Rtol):F3}");
        Assert.True(R <= R_ref * (1.0 + Rtol),
            $"R = {R:F3} above upper bound {R_ref * (1.0 + Rtol):F3}");
    }

    [Fact]
    [Trait("Category", "Roughness")]
    public void OutputDimensions_AreConsistent()
    {
        var (signal, fs) = MakeAmTone(44100);
        var result = RoughnessDw.Compute(signal, fs, overlap: 0.5);

        int nSeg = result.R.Length;
        Assert.True(nSeg > 1);
        Assert.Equal(47, result.RSpecific.GetLength(0));
        Assert.Equal(nSeg, result.RSpecific.GetLength(1));
        Assert.Equal(47, result.BarkAxis.Length);
        Assert.Equal(nSeg, result.TimeAxis.Length);
    }
}
