using Mosqito.Dsp;
using Mosqito.Io;
using Mosqito.SqMetrics.Loudness;
using Mosqito.Tests.Support;
using Xunit;
using Isoclose = Mosqito.Tests.Support.Isoclose;

namespace Mosqito.Tests.SqMetrics.Loudness;

/// <summary>
/// Ported from MoSQITo tests/sq_metrics/loudness/test_loudness_zwst.py.
/// Reference: ISO 532-1:2017 Annex B — compliance tolerance rtol=5%, atol=0.1.
/// </summary>
public class LoudnessZwstTests
{
    // ------------------------------------------------------------------
    // ISO 532-1 Annex B2 — third-octave levels for Test signal 1
    // (28 values, 25 Hz to 12.5 kHz third-octave bands)
    // ------------------------------------------------------------------
    private static readonly double[] TestSignal1Spectrum =
    {
        -60, -60, 78, 79, 89, 72, 80, 89, 75, 87, 85, 79, 86,
        80, 71, 70, 72, 71, 72, 74, 69, 65, 67, 77, 68, 58, 45, 30.0
    };

    // Reference paths (linked into test output via .csproj)
    private const string Input = "tests/input";
    private const string PinkNoiseWav = Input + "/Test signal 5 (pinknoise 60 dB).wav";
    private const string Tone44100Wav = Input + "/Test signal 3 (1 kHz 60 dB)_44100Hz.wav";
    private const string Signal1Csv   = Input + "/test_signal_1.csv";
    private const string Signal3Csv   = Input + "/test_signal_3.csv";
    private const string Signal5Csv   = Input + "/test_signal_5.csv";
    private const double WavCalib = 2.0 * 1.4142135623730951; // 2 * sqrt(2) as in Python

    // ------------------------------------------------------------------
    // Test 1: MainLoudness + CalcSlopes with direct third-octave input (ISO B2)
    // Mirrors test_loudness_zwst_3oct
    // ------------------------------------------------------------------
    [Fact]
    [Trait("Category", "Loudness")]
    public void ThirdOctaveSpectrum_MatchesIsoTarget_83_296()
    {
        const double N_iso = 83.296;
        double[] N_specif_iso = CsvLoader.LoadColumn(Signal1Csv, skipRows: 1);

        // Compute via MainLoudness → CalcSlopes (same path as Python _main_loudness + _calc_slopes)
        double[] nm = Mosqito.SqMetrics.Loudness.MainLoudness.Compute(TestSignal1Spectrum, "free");
        var (N, NSpecific) = CalcSlopes.Compute(nm);

        Isoclose.Assert(N, N_iso, rtol: 0.05, atol: 0.1, context: "N overall (Test signal 1)");
        Isoclose.Assert(NSpecific, N_specif_iso, rtol: 0.05, atol: 0.1, context: "N_specific (Test signal 1)");
    }

    // ------------------------------------------------------------------
    // Test 2: Time-domain WAV path (pink noise 60 dB, ISO B3)
    // Mirrors test_loudness_zwst_wav
    // ------------------------------------------------------------------
    [Fact]
    [Trait("Category", "Loudness")]
    public void PinkNoiseWav_MatchesIsoTarget_10_498()
    {
        const double N_iso = 10.498;
        double[] N_specif_iso = CsvLoader.LoadColumn(Signal5Csv, skipRows: 1);

        var (sig, fs) = WavLoader.Read(PinkNoiseWav, WavCalib);
        var result = LoudnessZwst.Compute(sig, fs);

        Isoclose.Assert(result.N, N_iso, rtol: 0.05, atol: 0.1, context: "N overall (Test signal 5, pink noise)");
        Isoclose.Assert(result.NSpecific, N_specif_iso, rtol: 0.05, atol: 0.1,
            context: "N_specific (Test signal 5, pink noise)");
    }

    // ------------------------------------------------------------------
    // Test 3: 44.1 kHz input with auto-resample (1 kHz tone 60 dB)
    // Mirrors test_loudness_zwst_44100Hz
    // ------------------------------------------------------------------
    [Fact]
    [Trait("Category", "Loudness")]
    public void Tone44100Hz_MatchesIsoTarget_4_019()
    {
        const double N_iso = 4.019;
        double[] N_specif_iso = CsvLoader.LoadColumn(Signal3Csv, skipRows: 1);

        var (sig, fs) = WavLoader.Read(Tone44100Wav, WavCalib);
        var result = LoudnessZwst.Compute(sig, fs);

        Isoclose.Assert(result.N, N_iso, rtol: 0.05, atol: 0.1, context: "N overall (Test signal 3, 1 kHz 60 dB 44.1 kHz)");
        Isoclose.Assert(result.NSpecific, N_specif_iso, rtol: 0.05, atol: 0.1,
            context: "N_specific (Test signal 3, 1 kHz 60 dB 44.1 kHz)");
    }

    // ------------------------------------------------------------------
    // Test 4: Per-segment — overall average must be within 5% of N_iso
    // Mirrors test_loudness_zwst_perseg
    // ------------------------------------------------------------------
    [Fact]
    [Trait("Category", "Loudness")]
    public void PinkNoise_PerSeg_AverageN_MatchesIsoTarget_10_498()
    {
        const double N_iso = 10.498;

        var (sig, fs) = WavLoader.Read(PinkNoiseWav, WavCalib);
        var (N, _, _, _) = LoudnessZwst.ComputePerSeg(sig, fs, nPerSeg: 8192 * 2, noOverlap: 4096);

        // Assert each segment within tolerance (matches np.testing.assert_allclose rtol=0.05)
        foreach (double n in N)
            Isoclose.Assert(n, N_iso, rtol: 0.05, atol: 0.1, context: "N per segment (pink noise)");
    }

    // ------------------------------------------------------------------
    // Test 5: Frequency-domain spectrum path (loudness_zwst_freq)
    // Mirrors test_loudness_zwst_freq
    // ------------------------------------------------------------------
    [Fact]
    [Trait("Category", "Loudness")]
    public void FrequencyDomainPath_PinkNoise_MatchesIsoTarget_10_498()
    {
        const double N_iso = 10.498;
        double[] N_specif_iso = CsvLoader.LoadColumn(Signal5Csv, skipRows: 1);

        var (sig, fs) = WavLoader.Read(PinkNoiseWav, WavCalib);

        // Build one-sided amplitude spectrum: |2 / sqrt(2) / n * FFT(sig)[0:n/2]|
        // Matches Python: spec = np.abs(2 / np.sqrt(2) / n * fft(sig)[0:n//2])
        int n = sig.Length;
        System.Numerics.Complex[] cmplx = Fft.Fft2(sig);
        int halfN = n / 2;
        double[] spectrum = new double[halfN];
        double scale = 2.0 / Math.Sqrt(2.0) / n;
        for (int i = 0; i < halfN; i++) spectrum[i] = cmplx[i].Magnitude * scale;

        double[] freqs = Fft.FftFreq(n, fs)[..halfN];

        var result = LoudnessZwst.ComputeFromSpectrum(spectrum, freqs);

        Isoclose.Assert(result.N, N_iso, rtol: 0.05, atol: 0.1, context: "N overall (freq path, pink noise)");
        Isoclose.Assert(result.NSpecific, N_specif_iso, rtol: 0.05, atol: 0.1,
            context: "N_specific (freq path, pink noise)");
    }
}
