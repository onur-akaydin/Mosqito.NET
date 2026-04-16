using BenchmarkDotNet.Attributes;
using Mosqito.Io;
using Mosqito.SqMetrics.SpeechIntelligibility;

namespace Mosqito.Benchmarks;

/// <summary>
/// Benchmarks for ANSI S3.5 Speech Intelligibility Index (SII).
/// WAV: "white_noise_200_2000_Hz_stationary" — noise signal is the natural input for SII.
/// Synthetic: 500 Hz sine, 60 dB SPL, 1 s @ 48 kHz (broadband proxy).
/// </summary>
[MemoryDiagnoser]
public class SpeechIntelligibilityBenchmarks
{
    private const string WavPath = "tests/input/white_noise_200_2000_Hz_stationary.wav";
    private const double WavCalib = 0.01;

    private double[] _synthetic = null!;
    private double[] _wavSignal = null!;
    private int _wavFs;

    [GlobalSetup]
    public void Setup()
    {
        (_synthetic, _) = SignalGenerators.SineWave(fs: 48000, d: 1.0, freq: 500, splLevel: 60);
        (_wavSignal, _wavFs) = Load.FromWav(WavPath, WavCalib);
    }

    // ── SiiAnsi ───────────────────────────────────────────────────────────

    [Benchmark]
    public (double SII, double[] SIISpecific, double[] FreqAxis) SiiAnsi_Synthetic()
        => SiiAnsi.Compute(_synthetic, 48000, method: "octave", speechLevel: "normal");

    [Benchmark]
    public (double SII, double[] SIISpecific, double[] FreqAxis) SiiAnsi_WavFile()
        => SiiAnsi.Compute(_wavSignal, _wavFs, method: "octave", speechLevel: "normal");
}
