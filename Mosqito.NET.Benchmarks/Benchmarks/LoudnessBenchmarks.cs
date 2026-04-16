using BenchmarkDotNet.Attributes;
using Mosqito.Dsp;
using Mosqito.Io;
using Mosqito.SqMetrics.Loudness;

namespace Mosqito.Benchmarks;

/// <summary>
/// Benchmarks for all Loudness algorithms.
/// WAV: "Test signal 5 (pinknoise 60 dB)" — the ISO 532-1 Annex B pink-noise reference.
/// Synthetic: 1 kHz sine, 60 dB SPL, 1 s @ 48 kHz.
/// </summary>
[MemoryDiagnoser]
public class LoudnessBenchmarks
{
    private const string WavPath = "tests/input/Test signal 5 (pinknoise 60 dB).wav";
    private const double WavCalib = 2.8284271247461903; // 2 * sqrt(2)

    private double[] _synthetic = null!;
    private double[] _wavSignal = null!;
    private int _wavFs;

    private double[] _syntheticSpectrum = null!;
    private double[] _syntheticFreqs = null!;
    private double[] _wavSpectrum = null!;
    private double[] _wavFreqs = null!;

    [GlobalSetup]
    public void Setup()
    {
        // Synthetic: 1 kHz sine, 60 dB, 1 s @ 48 kHz
        (_synthetic, _) = SignalGenerators.SineWave(fs: 48000, d: 1.0, freq: 1000, splLevel: 60);

        // WAV file
        (_wavSignal, _wavFs) = Load.FromWav(WavPath, WavCalib);

        // Pre-build frequency-domain inputs for the freq-path benchmarks
        (_syntheticSpectrum, _syntheticFreqs) = BuildSpectrum(_synthetic, 48000);
        (_wavSpectrum, _wavFreqs) = BuildSpectrum(_wavSignal, _wavFs);
    }

    private static (double[] spectrum, double[] freqs) BuildSpectrum(double[] signal, int fs)
    {
        int n = signal.Length;
        var cmplx = Fft.Fft2(signal);
        int halfN = n / 2;
        double[] spectrum = new double[halfN];
        double scale = 2.0 / Math.Sqrt(2.0) / n;
        for (int i = 0; i < halfN; i++)
            spectrum[i] = cmplx[i].Magnitude * scale;
        double[] freqs = Fft.FftFreq(n, fs)[..halfN];
        return (spectrum, freqs);
    }

    // ── LoudnessZwst ──────────────────────────────────────────────────────

    [Benchmark]
    public LoudnessZwstResult LoudnessZwst_Synthetic()
        => LoudnessZwst.Compute(_synthetic, 48000);

    [Benchmark]
    public LoudnessZwstResult LoudnessZwst_WavFile()
        => LoudnessZwst.Compute(_wavSignal, _wavFs);

    // ── LoudnessZwst (frequency domain) ───────────────────────────────────

    [Benchmark]
    public LoudnessZwstResult LoudnessZwstFreq_Synthetic()
        => LoudnessZwst.ComputeFromSpectrum(_syntheticSpectrum, _syntheticFreqs);

    [Benchmark]
    public LoudnessZwstResult LoudnessZwstFreq_WavFile()
        => LoudnessZwst.ComputeFromSpectrum(_wavSpectrum, _wavFreqs);

    // ── LoudnessZwtv ──────────────────────────────────────────────────────

    [Benchmark]
    public (double[] N, double[,] NSpecific, double[] BarkAxis, double[] TimeAxis) LoudnessZwtv_Synthetic()
        => LoudnessZwtv.Compute(_synthetic, 48000);

    [Benchmark]
    public (double[] N, double[,] NSpecific, double[] BarkAxis, double[] TimeAxis) LoudnessZwtv_WavFile()
        => LoudnessZwtv.Compute(_wavSignal, _wavFs);

    // ── LoudnessEcma ──────────────────────────────────────────────────────

    [Benchmark]
    public LoudnessEcmaResult LoudnessEcma_Synthetic()
        => LoudnessEcma.Compute(_synthetic, 48000);

    [Benchmark]
    public LoudnessEcmaResult LoudnessEcma_WavFile()
        => LoudnessEcma.Compute(_wavSignal, _wavFs);
}
