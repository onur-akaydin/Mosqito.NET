using BenchmarkDotNet.Attributes;
using Mosqito.Io;
using Mosqito.SoundLevelMeter;
using Mosqito.SqMetrics.Roughness;

namespace Mosqito.Benchmarks;

/// <summary>
/// Benchmarks for all Roughness algorithms.
/// WAV: "white_noise_442_1768_Hz_stationary" — band-limited noise representative of roughness input.
/// Synthetic: AM tone (1 kHz carrier, 70 Hz modulation, 60 dB) — the Daniel-Weber reference signal.
/// </summary>
[MemoryDiagnoser]
public class RoughnessBenchmarks
{
    private const string WavPath = "tests/input/white_noise_442_1768_Hz_stationary.wav";
    private const double WavCalib = 0.01;

    private double[] _syntheticAm = null!;
    private int _syntheticFs;
    private double[] _wavSignal = null!;
    private int _wavFs;

    private double[] _syntheticSpectrum = null!;
    private double[] _syntheticFreqs = null!;
    private double[] _wavSpectrum = null!;
    private double[] _wavFreqs = null!;

    [GlobalSetup]
    public void Setup()
    {
        // Synthetic: AM tone — Daniel-Weber reference (1 kHz @ 70 Hz mod, 60 dB, 1 s @ 44100 Hz)
        _syntheticFs = 44100;
        int n = _syntheticFs;
        double[] xmod = new double[n];
        for (int i = 0; i < n; i++)
            xmod[i] = Math.Sin(2.0 * Math.PI * 70.0 * i / _syntheticFs);
        (_syntheticAm, _) = SignalGenerators.AmSine(xmod, fs: _syntheticFs, fc: 1000.0, splLevel: 60.0);

        // WAV file
        (_wavSignal, _wavFs) = Load.FromWav(WavPath, WavCalib);

        // Pre-build frequency-domain inputs
        int nperseg = (int)(0.2 * _syntheticFs);
        (_syntheticSpectrum, _syntheticFreqs) = BuildBlackmanSpectrum(_syntheticAm, _syntheticFs, nperseg);

        int wavNperseg = (int)(0.2 * _wavFs);
        (_wavSpectrum, _wavFreqs) = BuildBlackmanSpectrum(_wavSignal, _wavFs, wavNperseg);
    }

    private static (double[] spectrum, double[] freqs) BuildBlackmanSpectrum(
        double[] signal, int fs, int nperseg)
    {
        var (spec, freqs) = CompSpectrum.Compute(
            signal.AsSpan(0, nperseg), fs,
            nfft: nperseg,
            window: CompSpectrum.WindowType.Blackman,
            oneSided: true, db: false);
        return (spec, freqs);
    }

    // ── RoughnessDw (time domain) ─────────────────────────────────────────

    [Benchmark]
    public RoughnessDwResult RoughnessDw_Synthetic()
        => RoughnessDw.Compute(_syntheticAm, _syntheticFs, overlap: 0);

    [Benchmark]
    public RoughnessDwResult RoughnessDw_WavFile()
        => RoughnessDw.Compute(_wavSignal, _wavFs, overlap: 0);

    // ── RoughnessDw (frequency domain) ────────────────────────────────────

    [Benchmark]
    public (double R, double[] RSpecific, double[] BarkAxis) RoughnessDwFreq_Synthetic()
        => RoughnessDw.ComputeFromSpectrum(_syntheticSpectrum, _syntheticFreqs);

    [Benchmark]
    public (double R, double[] RSpecific, double[] BarkAxis) RoughnessDwFreq_WavFile()
        => RoughnessDw.ComputeFromSpectrum(_wavSpectrum, _wavFreqs);

    // ── RoughnessEcma ─────────────────────────────────────────────────────

    [Benchmark]
    public (double R, double[] RTime, double[,] RSpecific, double[] BarkAxis, double[] TimeAxis) RoughnessEcma_Synthetic()
        => RoughnessEcma.Compute(_syntheticAm, _syntheticFs);

    [Benchmark]
    public (double R, double[] RTime, double[,] RSpecific, double[] BarkAxis, double[] TimeAxis) RoughnessEcma_WavFile()
        => RoughnessEcma.Compute(_wavSignal, _wavFs);
}
