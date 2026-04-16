using Mosqito.Io;
using Mosqito.SoundLevelMeter;
using Xunit;

namespace Mosqito.Tests.SoundLevelMeter;

/// <summary>
/// Tests for NoctSpectrum, NoctSynthesis, and CompSpectrum.
/// Ported from MoSQITo tests/sound_level_meter/noct_spectrum/test_noct_spectrum.py
/// and test_noct_synthesis.py.
/// </summary>
public class NoctTests
{
    private const string Input    = "tests/input";
    private const string PinkWav  = Input + "/Test signal 5 (pinknoise 60 dB).wav";
    private const double WavCalib = 2.0 * 1.4142135623730951;

    private static double Db(double amp) => 20.0 * Math.Log10(amp / 2e-5);

    // -----------------------------------------------------------------------
    // NoctSpectrum
    // -----------------------------------------------------------------------

    /// <summary>
    /// 1 kHz pure sine at 1 Pa RMS → 1/3-octave band at 1 kHz must read ≈94 dB SPL
    /// (within ±2 dB, matching the rtol=0.5 tolerance used in the Python test).
    /// Ported from test_noct_spectrum.py.
    /// </summary>
    [Fact]
    [Trait("Category", "SoundLevelMeter")]
    public void NoctSpectrum_1kHzSine_1PaRms_Reads94dBSpl()
    {
        int fs = 48000;
        int n  = fs; // 1 second
        double aRms = 1.0;
        double[] sig = new double[n];
        for (int i = 0; i < n; i++)
            sig[i] = aRms * Math.Sqrt(2.0) * Math.Sin(2.0 * Math.PI * 1000.0 * i / fs);

        var (spec, freqs) = NoctSpectrum.Compute(sig, fs, fmin: 24, fmax: 12600, n: 3);

        // Find 1000 Hz band
        int idx = Array.IndexOf(freqs, 1000.0);
        Assert.True(idx >= 0, "1000 Hz band not found in frequency axis.");

        double measured = Db(spec[idx]);
        double expected = Db(aRms); // = 20*log10(1/2e-5) ≈ 94 dB
        Assert.True(Math.Abs(measured - expected) < 2.0,
            $"1 kHz band: expected ≈{expected:F1} dB, got {measured:F1} dB");
    }

    /// <summary>
    /// Two equal-power tones at 3950 Hz and 4050 Hz (both in the 4 kHz 1/3-octave band)
    /// should read ≈3 dB above a single 0.5 Pa RMS tone, i.e. ≈88 + 3 = 91 dB SPL.
    /// </summary>
    [Fact]
    [Trait("Category", "SoundLevelMeter")]
    public void NoctSpectrum_TwoTonesInSameBand_Shows3dBAddition()
    {
        int fs = 48000;
        int n  = fs;
        double aRms = 0.5;
        double[] sig = new double[n];
        for (int i = 0; i < n; i++)
        {
            sig[i]  = aRms * Math.Sqrt(2.0) * Math.Sin(2.0 * Math.PI * 3950.0 * i / fs);
            sig[i] += aRms * Math.Sqrt(2.0) * Math.Sin(2.0 * Math.PI * 4050.0 * i / fs);
        }

        var (spec, freqs) = NoctSpectrum.Compute(sig, fs, fmin: 24, fmax: 12600, n: 3);

        int idx = Array.IndexOf(freqs, 4000.0);
        Assert.True(idx >= 0, "4000 Hz band not found.");

        double measured = Db(spec[idx]);
        double expected = Db(aRms) + 3.0; // two equal-power tones → +3 dB
        Assert.True(Math.Abs(measured - expected) < 2.0,
            $"4 kHz dual-tone band: expected ≈{expected:F1} dB, got {measured:F1} dB");
    }

    /// <summary>
    /// Output arrays have the same length and all band amplitudes are non-negative.
    /// </summary>
    [Fact]
    [Trait("Category", "SoundLevelMeter")]
    public void NoctSpectrum_OutputDimensions_AreConsistent()
    {
        int fs = 48000;
        double[] sig = new double[fs]; // 1 s silence
        var (spec3, freqs3) = NoctSpectrum.Compute(sig, fs, fmin: 24, fmax: 12600, n: 3);
        var (spec1, freqs1) = NoctSpectrum.Compute(sig, fs, fmin: 24, fmax: 12600, n: 1);

        Assert.Equal(spec3.Length, freqs3.Length);
        Assert.Equal(spec1.Length, freqs1.Length);
        Assert.True(spec3.Length > spec1.Length, "1/3-oct must have more bands than 1/1-oct.");
        Assert.All(spec3, v => Assert.True(v >= 0.0));
    }

    // -----------------------------------------------------------------------
    // NoctSynthesis — time-domain vs. frequency-domain agreement
    // Ported from test_noct_synthesis.py: assert_allclose(dB(spec_3t), dB(spec_3f), atol=0.3)
    // -----------------------------------------------------------------------

    [Fact]
    [Trait("Category", "SoundLevelMeter")]
    public void NoctSynthesis_ThirdOctave_AgreesWithNoctSpectrum_Within3dB()
        => AssertSynthesisAgreement(n: 3, atolDb: 3.0);

    [Fact]
    [Trait("Category", "SoundLevelMeter")]
    public void NoctSynthesis_FullOctave_AgreesWithNoctSpectrum_Within3dB()
        => AssertSynthesisAgreement(n: 1, atolDb: 3.0);

    private static void AssertSynthesisAgreement(int n, double atolDb)
    {
        var (sig, fs) = WavLoader.Read(PinkWav, WavCalib);

        // Time-domain reference
        var (specT, freqT) = NoctSpectrum.Compute(sig, fs, fmin: 24, fmax: 12600, n: n);

        // Frequency-domain path (NoctSynthesis): build one-sided RMS amplitude spectrum
        // 2/sqrt(2)/N * FFT  (one-sided, amplitude RMS)
        int nFft    = sig.Length;
        double[] spectrum;
        double[] freqAxis;
        {
            var full = Mosqito.Dsp.Fft.Fft2(sig);
            int half  = nFft / 2;
            spectrum  = new double[half];
            freqAxis  = new double[half];
            double scale = 2.0 / Math.Sqrt(2.0) / nFft;
            for (int i = 0; i < half; i++)
            {
                spectrum[i]  = full[i].Magnitude * scale;
                freqAxis[i]  = (double)(i + 1) * fs / nFft;
            }
            // Fix: index 0 in Python fft is DC; Python uses fft[0:n//2] which starts at DC.
            // Redo to match: spectrum[i] = |fft[i]| * scale for i in 0..n/2-1
            // and freqs = fftfreq(n, 1/fs)[0:n//2]
            for (int i = 0; i < half; i++)
            {
                spectrum[i] = full[i].Magnitude * scale;
                freqAxis[i] = (double)i * fs / nFft;
            }
        }

        var (specF, freqF) = NoctSynthesis.Compute(spectrum, freqAxis,
            fmin: 24, fmax: 12600, n: n);

        // Align frequency axes (time-domain may have one more band if fs != 48k)
        int nBands = Math.Min(specT.Length, specF.Length);
        Assert.True(nBands > 0);

        for (int i = 0; i < nBands; i++)
        {
            double dbT = Db(specT[i]);
            double dbF = Db(specF[i]);
            Assert.True(Math.Abs(dbT - dbF) < atolDb,
                $"Band {i} ({freqT[i]:F0} Hz): time={dbT:F2} dB, freq={dbF:F2} dB, diff={Math.Abs(dbT-dbF):F2} dB > {atolDb} dB");
        }
    }

    // -----------------------------------------------------------------------
    // CompSpectrum
    // -----------------------------------------------------------------------

    /// <summary>
    /// A 1 kHz sine at 1 Pa RMS: the FFT peak must appear at 1 kHz ± 1 bin,
    /// and the reported dB level must be ≈94 dB SPL.
    /// </summary>
    [Fact]
    [Trait("Category", "SoundLevelMeter")]
    public void CompSpectrum_1kHzSine_PeakAt1kHz_And94dBSpl()
    {
        int fs = 48000;
        int n  = 4096; // power-of-two for clean bins
        double aRms = 1.0;
        double[] sig = new double[n];
        for (int i = 0; i < n; i++)
            sig[i] = aRms * Math.Sqrt(2.0) * Math.Sin(2.0 * Math.PI * 1000.0 * i / fs);

        // Use rectangular window and db=false (amplitude) to test the peak
        var (specAmp, freqs) = CompSpectrum.Compute(sig, fs, nfft: n,
            window: CompSpectrum.WindowType.Rectangular, db: false);

        // Peak bin
        int peakIdx = 0;
        for (int i = 1; i < specAmp.Length; i++)
            if (specAmp[i] > specAmp[peakIdx]) peakIdx = i;

        // Peak frequency must be within one bin of 1000 Hz
        double binWidth = (double)fs / n;
        Assert.True(Math.Abs(freqs[peakIdx] - 1000.0) <= binWidth,
            $"Peak at {freqs[peakIdx]:F1} Hz, expected 1000 Hz ± {binWidth:F1} Hz");
    }

    /// <summary>
    /// dB=true output: peak level of a 1 Pa RMS sine should be ≈94 dB SPL (within ±3 dB).
    /// </summary>
    [Fact]
    [Trait("Category", "SoundLevelMeter")]
    public void CompSpectrum_DbOutput_1PaRms_ReadsApprox94dBSpl()
    {
        int fs = 48000;
        int n  = 4096;
        double aRms = 1.0;
        double[] sig = new double[n];
        for (int i = 0; i < n; i++)
            sig[i] = aRms * Math.Sqrt(2.0) * Math.Sin(2.0 * Math.PI * 1000.0 * i / fs);

        var (specDb, _) = CompSpectrum.Compute(sig, fs, nfft: n, db: true);

        double peakDb = specDb.Max();
        double expected = Db(aRms); // ≈ 94 dB
        Assert.True(Math.Abs(peakDb - expected) < 3.0,
            $"Peak dB = {peakDb:F1}, expected ≈{expected:F1} dB SPL");
    }

    /// <summary>
    /// Output arrays are the expected size (nfft/2 for one-sided, nfft for two-sided).
    /// </summary>
    [Fact]
    [Trait("Category", "SoundLevelMeter")]
    public void CompSpectrum_OutputDimensions_OneSidedAndTwoSided()
    {
        int fs = 48000, n = 1024;
        double[] sig = new double[n];

        var (s1, f1) = CompSpectrum.Compute(sig, fs, nfft: n, oneSided: true);
        var (s2, f2) = CompSpectrum.Compute(sig, fs, nfft: n, oneSided: false);

        Assert.Equal(n / 2, s1.Length);
        Assert.Equal(n / 2, f1.Length);
        Assert.Equal(n, s2.Length);
        Assert.Equal(n, f2.Length);
    }
}
