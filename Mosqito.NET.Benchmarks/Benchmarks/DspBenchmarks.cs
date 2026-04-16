using System.Numerics;
using BenchmarkDotNet.Attributes;
using Mosqito.Dsp;
using Mosqito.Io;

namespace Mosqito.Benchmarks;

/// <summary>
/// Benchmarks for core DSP primitives: FFT and IIR filtering.
/// WAV: "Test signal 5 (pinknoise 60 dB)" — broadband signal for real-world filtering cost.
/// Synthetic: 1 kHz sine, 60 dB SPL, 1 s @ 48 kHz.
/// </summary>
[MemoryDiagnoser]
public class DspBenchmarks
{
    private const string WavPath = "tests/input/Test signal 5 (pinknoise 60 dB).wav";
    private const double WavCalib = 2.8284271247461903; // 2 * sqrt(2)

    [Params(1024, 8192, 65536)]
    public int N;

    private double[] _syntheticFull = null!;
    private double[] _wavSignalFull = null!;

    // Slices sized to N for FFT benchmarks (set in IterationSetup)
    private double[] _syntheticSlice = null!;
    private double[] _wavSlice = null!;

    private double[,] _sos = null!;

    [GlobalSetup]
    public void GlobalSetup()
    {
        // Synthetic: 1 kHz sine at 48 kHz, long enough for any N
        int len = Math.Max(65536 + 1, 48000);
        (_syntheticFull, _) = SignalGenerators.SineWave(fs: 48000, d: (double)len / 48000 + 0.1, freq: 1000, splLevel: 60);

        // WAV file
        var (wav, _) = Load.FromWav(WavPath, WavCalib);
        // Pad or trim to at least 65536 samples
        if (wav.Length >= 65536)
        {
            _wavSignalFull = wav;
        }
        else
        {
            _wavSignalFull = new double[65536];
            wav.CopyTo(_wavSignalFull, 0);
        }

        // 8th-order Butterworth lowpass at Wn = 0.5 (fs/4 cutoff) — used by SosFilter benchmarks
        _sos = Butter.DesignLowpass(order: 8, wCutoff: 0.5);
    }

    [IterationSetup]
    public void IterationSetup()
    {
        _syntheticSlice = _syntheticFull[..N];
        _wavSlice = _wavSignalFull[..N];
    }

    // ── FFT ───────────────────────────────────────────────────────────────

    [Benchmark]
    public Complex[] Fft_Synthetic() => Fft.Fft2(_syntheticSlice);

    [Benchmark]
    public Complex[] Fft_WavFile() => Fft.Fft2(_wavSlice);

    // ── SOS IIR Filter (full 1-second signal, N not applied) ──────────────

    [Benchmark]
    public double[] SosFilter_Synthetic() => SosFilter.Process(_sos, _syntheticFull[..48000]);

    [Benchmark]
    public double[] SosFilter_WavFile() => SosFilter.Process(_sos, _wavSignalFull[..Math.Min(48000, _wavSignalFull.Length)]);
}
