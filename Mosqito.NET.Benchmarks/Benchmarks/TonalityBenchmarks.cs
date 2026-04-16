using BenchmarkDotNet.Attributes;
using Mosqito.Io;
using Mosqito.SqMetrics.Tonality;

namespace Mosqito.Benchmarks;

/// <summary>
/// Benchmarks for ECMA-74 Tonality (TNR and PR).
/// WAV: "Test signal 10 (tone pulse 1 kHz 10 ms 70 dB)" — tonal content required for TNR/PR peaks.
/// Synthetic: 1 kHz sine, 70 dB SPL, 1 s @ 48 kHz.
/// </summary>
[MemoryDiagnoser]
public class TonalityBenchmarks
{
    private const string WavPath = "tests/input/Test signal 10 (tone pulse 1 kHz 10 ms 70 dB).wav";
    private const double WavCalib = 2.8284271247461903; // 2 * sqrt(2)

    private double[] _synthetic = null!;
    private double[] _wavSignal = null!;
    private int _wavFs;

    [GlobalSetup]
    public void Setup()
    {
        (_synthetic, _) = SignalGenerators.SineWave(fs: 48000, d: 1.0, freq: 1000, splLevel: 70);
        (_wavSignal, _wavFs) = Load.FromWav(WavPath, WavCalib);
    }

    // ── TnrEcma Stationary ────────────────────────────────────────────────

    [Benchmark]
    public TonalityResult TnrEcmaSt_Synthetic()
        => TnrEcma.ComputeSt(_synthetic, 48000, prominentOnly: false);

    [Benchmark]
    public TonalityResult TnrEcmaSt_WavFile()
        => TnrEcma.ComputeSt(_wavSignal, _wavFs, prominentOnly: false);

    // ── PrEcma Stationary ─────────────────────────────────────────────────

    [Benchmark]
    public TonalityResult PrEcmaSt_Synthetic()
        => PrEcma.ComputeSt(_synthetic, 48000, prominentOnly: false);

    [Benchmark]
    public TonalityResult PrEcmaSt_WavFile()
        => PrEcma.ComputeSt(_wavSignal, _wavFs, prominentOnly: false);
}
