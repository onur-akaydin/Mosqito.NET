using BenchmarkDotNet.Attributes;
using Mosqito.Io;
using Mosqito.SqMetrics.Sharpness;

namespace Mosqito.Benchmarks;

/// <summary>
/// Benchmarks for DIN 45692 Sharpness.
/// WAV: "broadband_570" — broadband signal that exercises the full sharpness algorithm.
/// Synthetic: 1 kHz sine, 60 dB SPL, 1 s @ 48 kHz.
/// </summary>
[MemoryDiagnoser]
public class SharpnessBenchmarks
{
    private const string WavPath = "tests/input/broadband_570.wav";
    private const double WavCalib = 1.0;

    private double[] _synthetic = null!;
    private double[] _wavSignal = null!;
    private int _wavFs;

    [GlobalSetup]
    public void Setup()
    {
        (_synthetic, _) = SignalGenerators.SineWave(fs: 48000, d: 1.0, freq: 1000, splLevel: 60);
        (_wavSignal, _wavFs) = Load.FromWav(WavPath, WavCalib);
    }

    // ── SharpnessDin Stationary ───────────────────────────────────────────

    [Benchmark]
    public double SharpnessDinSt_Synthetic()
        => SharpnessDin.ComputeSt(_synthetic, 48000);

    [Benchmark]
    public double SharpnessDinSt_WavFile()
        => SharpnessDin.ComputeSt(_wavSignal, _wavFs);

    // ── SharpnessDin Time-Varying ─────────────────────────────────────────

    [Benchmark]
    public (double[] S, double[] TimeAxis) SharpnessDinTv_Synthetic()
        => SharpnessDin.ComputeTv(_synthetic, 48000);

    [Benchmark]
    public (double[] S, double[] TimeAxis) SharpnessDinTv_WavFile()
        => SharpnessDin.ComputeTv(_wavSignal, _wavFs);
}
