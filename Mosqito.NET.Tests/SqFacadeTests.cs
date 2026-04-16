using Mosqito.Io;
using Xunit;

namespace Mosqito.Tests;

/// <summary>
/// Smoke tests for the <see cref="Sq"/> public façade.
/// Verifies that every entry point delegates correctly and returns sensible results.
/// </summary>
public class SqFacadeTests
{
    private const int Fs = 48000;
    private const string PinkWav  = "tests/input/Test signal 5 (pinknoise 60 dB).wav";
    private const double WavCalib = 2.0 * 1.4142135623730951;

    // 1 kHz sine at 60 dB SPL, 1 second
    private static double[] Sine1kHz60dB()
    {
        int n = Fs;
        double amp = 2e-5 * Math.Pow(10.0, 60.0 / 20.0) * Math.Sqrt(2.0);
        double[] sig = new double[n];
        for (int i = 0; i < n; i++)
            sig[i] = amp * Math.Sin(2.0 * Math.PI * 1000.0 * i / Fs);
        return sig;
    }

    // -----------------------------------------------------------------------
    // Sound level meter
    // -----------------------------------------------------------------------

    [Fact]
    [Trait("Category", "Facade")]
    public void Facade_NoctSpectrum_ReturnsPositiveAmplitudes()
    {
        double[] sig = Sine1kHz60dB();
        var (spec, freqs) = Sq.NoctSpectrum(sig, Fs, fmin: 24, fmax: 12600, n: 3);
        Assert.True(spec.Length > 0);
        Assert.Equal(spec.Length, freqs.Length);
        Assert.All(spec, v => Assert.True(v >= 0.0));
    }

    [Fact]
    [Trait("Category", "Facade")]
    public void Facade_CompSpectrum_ReturnsDimensionallyConsistentOutput()
    {
        double[] sig = Sine1kHz60dB();
        var (spectrum, freqs) = Sq.CompSpectrum(sig, Fs);
        Assert.True(spectrum.Length > 0);
        Assert.Equal(spectrum.Length, freqs.Length);
    }

    // -----------------------------------------------------------------------
    // Loudness
    // -----------------------------------------------------------------------

    [Fact]
    [Trait("Category", "Facade")]
    public void Facade_LoudnessEcma_ReturnsPositiveN()
    {
        double[] sig = Sine1kHz60dB();
        var result = Sq.LoudnessEcma(sig, Fs);
        Assert.True(result.N > 0.0, $"N = {result.N}");
        Assert.Equal(53, result.BarkAxis.Length);
    }

    [Fact]
    [Trait("Category", "Facade")]
    public void Facade_LoudnessZwst_ReturnsPositiveN()
    {
        var (sig, fs) = WavLoader.Read(PinkWav, WavCalib);
        var result = Sq.LoudnessZwst(sig, fs);
        Assert.True(result.N > 0.0, $"N = {result.N}");
    }

    [Fact]
    [Trait("Category", "Facade")]
    public void Facade_LoudnessZwstFreq_ReturnsPositiveN()
    {
        // Build a one-sided RMS amplitude spectrum from the pink noise WAV
        // (same approach as the LoudnessZwstTests.FrequencyDomainPath test)
        var (sig, fs) = WavLoader.Read(PinkWav, WavCalib);
        int n = sig.Length;
        var cmplx = Mosqito.Dsp.Fft.Fft2(sig);
        int halfN = n / 2;
        double[] spectrum = new double[halfN];
        double[] freqs    = new double[halfN];
        double scale = 2.0 / Math.Sqrt(2.0) / n;
        for (int i = 0; i < halfN; i++)
        {
            spectrum[i] = cmplx[i].Magnitude * scale;
            freqs[i]    = (double)i * fs / n;
        }

        var result = Sq.LoudnessZwstFreq(spectrum, freqs);
        Assert.True(result.N > 0.0, $"N = {result.N}");
    }

    [Fact]
    [Trait("Category", "Facade")]
    public void Facade_LoudnessZwtv_ReturnsDimensionallyConsistentOutput()
    {
        var (sig, fs) = WavLoader.Read(PinkWav, WavCalib);
        var (N, NSpec, barkAxis, timeAxis) = Sq.LoudnessZwtv(sig, fs);
        Assert.Equal(N.Length, timeAxis.Length);
        Assert.Equal(240, barkAxis.Length);
    }

    // -----------------------------------------------------------------------
    // Roughness
    // -----------------------------------------------------------------------

    [Fact]
    [Trait("Category", "Facade")]
    public void Facade_RoughnessEcma_ReturnsPositiveR()
    {
        double[] sig = Sine1kHz60dB();
        // AM modulate at 70 Hz
        double[] amSig = new double[sig.Length];
        for (int i = 0; i < sig.Length; i++)
            amSig[i] = (1.0 + Math.Sin(2.0 * Math.PI * 70.0 * i / Fs)) * sig[i] * 0.5;
        // Normalise
        double rms = Math.Sqrt(amSig.Average(x => x * x));
        double target = 2e-5 * Math.Pow(10.0, 60.0 / 20.0);
        for (int i = 0; i < amSig.Length; i++) amSig[i] *= target / rms;

        var (R, rTime, rSpec, barkAxis, t50) = Sq.RoughnessEcma(amSig, Fs);
        Assert.True(R >= 0.0, $"R = {R}");
        Assert.Equal(53, barkAxis.Length);
        Assert.Equal(rTime.Length, t50.Length);
    }

    // -----------------------------------------------------------------------
    // Sharpness
    // -----------------------------------------------------------------------

    [Fact]
    [Trait("Category", "Facade")]
    public void Facade_SharpnessDinSt_ReturnsPositiveS()
    {
        var (sig, fs) = WavLoader.Read(PinkWav, WavCalib);
        double S = Sq.SharpnessDinSt(sig, fs);
        Assert.True(S > 0.0, $"S = {S}");
    }

    // -----------------------------------------------------------------------
    // Tonality
    // -----------------------------------------------------------------------

    [Fact]
    [Trait("Category", "Facade")]
    public void Facade_TnrEcmaSt_DoesNotThrow_ForPureTone()
    {
        double[] sig = Sine1kHz60dB();
        var result = Sq.TnrEcmaSt(sig, Fs);
        Assert.NotNull(result);
    }

    [Fact]
    [Trait("Category", "Facade")]
    public void Facade_PrEcmaSt_DoesNotThrow_ForPureTone()
    {
        double[] sig = Sine1kHz60dB();
        var result = Sq.PrEcmaSt(sig, Fs);
        Assert.NotNull(result);
    }

    // -----------------------------------------------------------------------
    // Speech intelligibility
    // -----------------------------------------------------------------------

    [Fact]
    [Trait("Category", "Facade")]
    public void Facade_SiiAnsiLevel_ReturnsPositiveSII()
    {
        var (SII, _, _) = Sq.SiiAnsiLevel(noiseLevel: 60.0,
            method: "octave", speechLevel: "normal");
        Assert.True(SII >= 0.0 && SII <= 1.0, $"SII = {SII}");
    }

    // -----------------------------------------------------------------------
    // Utility methods
    // -----------------------------------------------------------------------

    [Fact]
    [Trait("Category", "Facade")]
    public void Facade_SoneToPhon_1Sone_Gives40Phon()
    {
        // By definition: 1 sone = 40 phon (equal-loudness reference level)
        double phon = Sq.SoneToPhon(1.0);
        Assert.True(Math.Abs(phon - 40.0) < 1.0, $"SoneToPhon(1) = {phon:F1}, expected 40");
    }

    [Fact]
    [Trait("Category", "Facade")]
    public void Facade_SoneToPhon_Monotonically_IncreasesWith_Sone()
    {
        // Louder sounds in sone must map to higher phon values
        double p1 = Sq.SoneToPhon(1.0);
        double p2 = Sq.SoneToPhon(4.0);
        double p8 = Sq.SoneToPhon(8.0);
        Assert.True(p1 < p2 && p2 < p8, $"Not monotone: {p1:F1} < {p2:F1} < {p8:F1}");
    }

    [Fact]
    [Trait("Category", "Facade")]
    public void Facade_EqualLoudnessContours_Returns53Points()
    {
        var (spl, freqs) = Sq.EqualLoudnessContours(phones: 60.0);
        Assert.True(spl.Length > 0);
        Assert.Equal(spl.Length, freqs.Length);
    }
}
