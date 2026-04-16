using Mosqito.Conversion;
using Mosqito.SoundLevelMeter;
using Mosqito.SqMetrics.Loudness;
using Mosqito.SqMetrics.Roughness;
using Mosqito.SqMetrics.Sharpness;
using Mosqito.SqMetrics.SpeechIntelligibility;
using Mosqito.SqMetrics.Tonality;

namespace Mosqito;

/// <summary>
/// Single-entry-point façade for Mosqito.NET — mirrors the flat public namespace of
/// the MoSQITo Python library (<c>from mosqito import ...</c>).
///
/// All methods delegate directly to the strongly-typed classes in the
/// <c>Mosqito.SqMetrics.*</c>, <c>Mosqito.SoundLevelMeter</c>, and
/// <c>Mosqito.Conversion</c> namespaces.  Use those classes directly for
/// full IntelliSense; use this façade for exploratory / scripting work.
/// </summary>
public static class Sq
{
    // =========================================================================
    // Sound level meter
    // =========================================================================

    /// <inheritdoc cref="Mosqito.SoundLevelMeter.NoctSpectrum.Compute"/>
    public static (double[] Spec, double[] FreqAxis) NoctSpectrum(
        ReadOnlySpan<double> signal, double fs,
        double fmin, double fmax,
        int n = 3, int g = 10, double fr = 1000.0)
        => Mosqito.SoundLevelMeter.NoctSpectrum.Compute(signal, fs, fmin, fmax, n, g, fr);

    /// <inheritdoc cref="Mosqito.SoundLevelMeter.CompSpectrum.Compute"/>
    public static (double[] Spectrum, double[] FreqAxis) CompSpectrum(
        ReadOnlySpan<double> signal, int fs,
        int nfft = -1,
        Mosqito.SoundLevelMeter.CompSpectrum.WindowType window =
            Mosqito.SoundLevelMeter.CompSpectrum.WindowType.Hanning,
        bool oneSided = true,
        bool db = true)
        => Mosqito.SoundLevelMeter.CompSpectrum.Compute(signal, fs, nfft, window, oneSided, db);

    // =========================================================================
    // Loudness
    // =========================================================================

    /// <inheritdoc cref="LoudnessEcma.Compute"/>
    public static LoudnessEcmaResult LoudnessEcma(
        ReadOnlySpan<double> signal, int fs,
        int sb = 2048, int sh = 1024)
        => Mosqito.SqMetrics.Loudness.LoudnessEcma.Compute(signal, fs, sb, sh);

    /// <inheritdoc cref="LoudnessZwst.Compute"/>
    public static LoudnessZwstResult LoudnessZwst(
        ReadOnlySpan<double> signal, int fs,
        string fieldType = "free")
        => Mosqito.SqMetrics.Loudness.LoudnessZwst.Compute(signal, fs, fieldType);

    /// <inheritdoc cref="LoudnessZwst.ComputeFromSpectrum"/>
    public static LoudnessZwstResult LoudnessZwstFreq(
        ReadOnlySpan<double> spectrum, ReadOnlySpan<double> freqs,
        string fieldType = "free")
        => Mosqito.SqMetrics.Loudness.LoudnessZwst.ComputeFromSpectrum(spectrum, freqs, fieldType);

    /// <inheritdoc cref="LoudnessZwst.ComputePerSeg"/>
    public static (double[] N, double[,] NSpecific, double[] BarkAxis, double[] TimeAxis)
        LoudnessZwstPerSeg(
            ReadOnlySpan<double> signal, int fs,
            int nPerSeg = 4096, int? noOverlap = null,
            string fieldType = "free")
        => Mosqito.SqMetrics.Loudness.LoudnessZwst.ComputePerSeg(signal, fs, nPerSeg, noOverlap, fieldType);

    /// <inheritdoc cref="LoudnessZwtv.Compute"/>
    public static (double[] N, double[,] NSpecific, double[] BarkAxis, double[] TimeAxis)
        LoudnessZwtv(ReadOnlySpan<double> signal, int fs, string fieldType = "free")
        => Mosqito.SqMetrics.Loudness.LoudnessZwtv.Compute(signal, fs, fieldType);

    /// <inheritdoc cref="LoudnessUtils.EqualLoudnessContours"/>
    public static (double[] Spl, double[] Freqs) EqualLoudnessContours(double phones)
        => LoudnessUtils.EqualLoudnessContours(phones);

    /// <inheritdoc cref="LoudnessUtils.SoneToPhon"/>
    public static double SoneToPhon(double sone)
        => LoudnessUtils.SoneToPhon(sone);

    // =========================================================================
    // Roughness
    // =========================================================================

    /// <inheritdoc cref="RoughnessEcma.Compute"/>
    public static (double R, double[] RTime, double[,] RSpecific, double[] BarkAxis, double[] TimeAxis)
        RoughnessEcma(ReadOnlySpan<double> signal, int fs)
        => Mosqito.SqMetrics.Roughness.RoughnessEcma.Compute(signal, fs);

    /// <inheritdoc cref="RoughnessDw.Compute"/>
    public static RoughnessDwResult RoughnessDw(
        ReadOnlySpan<double> signal, int fs, double overlap = 0.5)
        => Mosqito.SqMetrics.Roughness.RoughnessDw.Compute(signal, fs, overlap);

    /// <inheritdoc cref="RoughnessDw.ComputeFromSpectrum"/>
    public static (double R, double[] RSpecific, double[] BarkAxis) RoughnessDwFreq(
        ReadOnlySpan<double> spectrum, ReadOnlySpan<double> freqs)
        => Mosqito.SqMetrics.Roughness.RoughnessDw.ComputeFromSpectrum(spectrum, freqs);

    // =========================================================================
    // Sharpness
    // =========================================================================

    /// <inheritdoc cref="SharpnessDin.FromLoudness(double,ReadOnlySpan{double},SharpnessWeighting)"/>
    public static double SharpnessDinFromLoudness(
        double N, ReadOnlySpan<double> nSpecific,
        SharpnessWeighting weighting = SharpnessWeighting.Din)
        => SharpnessDin.FromLoudness(N, nSpecific, weighting);

    /// <inheritdoc cref="SharpnessDin.ComputeSt"/>
    public static double SharpnessDinSt(
        ReadOnlySpan<double> signal, int fs,
        SharpnessWeighting weighting = SharpnessWeighting.Din,
        string fieldType = "free")
        => SharpnessDin.ComputeSt(signal, fs, weighting, fieldType);

    /// <inheritdoc cref="SharpnessDin.ComputeTv"/>
    public static (double[] S, double[] TimeAxis) SharpnessDinTv(
        ReadOnlySpan<double> signal, int fs,
        SharpnessWeighting weighting = SharpnessWeighting.Din,
        string fieldType = "free")
        => SharpnessDin.ComputeTv(signal, fs, weighting, fieldType);

    /// <inheritdoc cref="SharpnessDin.ComputePerSeg"/>
    public static (double[] S, double[] TimeAxis) SharpnessDinPerSeg(
        ReadOnlySpan<double> signal, int fs,
        int nPerSeg = 4096, int? noOverlap = null,
        SharpnessWeighting weighting = SharpnessWeighting.Din,
        string fieldType = "free")
        => SharpnessDin.ComputePerSeg(signal, fs, nPerSeg, noOverlap, weighting, fieldType);

    /// <inheritdoc cref="SharpnessDin.ComputeFreq"/>
    public static double SharpnessDinFreq(
        ReadOnlySpan<double> spectrum, ReadOnlySpan<double> freqs,
        SharpnessWeighting weighting = SharpnessWeighting.Din,
        string fieldType = "free")
        => SharpnessDin.ComputeFreq(spectrum, freqs, weighting, fieldType);

    // =========================================================================
    // Tonality (TNR + PR, ECMA-74)
    // =========================================================================

    /// <inheritdoc cref="TnrEcma.ComputeSt"/>
    public static TonalityResult TnrEcmaSt(
        ReadOnlySpan<double> signal, int fs, bool prominentOnly = true)
        => TnrEcma.ComputeSt(signal, fs, prominentOnly);

    /// <inheritdoc cref="TnrEcma.ComputeFreq"/>
    public static TonalityResult TnrEcmaFreq(
        ReadOnlySpan<double> specDb, ReadOnlySpan<double> freqAxis,
        bool prominentOnly = true)
        => TnrEcma.ComputeFreq(specDb, freqAxis, prominentOnly);

    /// <inheritdoc cref="TnrEcma.ComputePerSeg"/>
    public static (TonalityResult[] Results, double[] TimeAxis) TnrEcmaPerSeg(
        ReadOnlySpan<double> signal, int fs,
        int nPerSeg = 4096, int? noOverlap = null, bool prominentOnly = true)
        => TnrEcma.ComputePerSeg(signal, fs, nPerSeg, noOverlap, prominentOnly);

    /// <inheritdoc cref="PrEcma.ComputeSt"/>
    public static TonalityResult PrEcmaSt(
        ReadOnlySpan<double> signal, int fs, bool prominentOnly = true)
        => PrEcma.ComputeSt(signal, fs, prominentOnly);

    /// <inheritdoc cref="PrEcma.ComputeFreq"/>
    public static TonalityResult PrEcmaFreq(
        ReadOnlySpan<double> specDb, ReadOnlySpan<double> freqAxis,
        bool prominentOnly = true)
        => PrEcma.ComputeFreq(specDb, freqAxis, prominentOnly);

    /// <inheritdoc cref="PrEcma.ComputePerSeg"/>
    public static (TonalityResult[] Results, double[] TimeAxis) PrEcmaPerSeg(
        ReadOnlySpan<double> signal, int fs,
        int nPerSeg = 4096, int? noOverlap = null, bool prominentOnly = true)
        => PrEcma.ComputePerSeg(signal, fs, nPerSeg, noOverlap, prominentOnly);

    // =========================================================================
    // Speech intelligibility (SII, ANSI S3.5)
    // =========================================================================

    /// <inheritdoc cref="SiiAnsi.Compute"/>
    public static (double SII, double[] SIISpecific, double[] FreqAxis) SiiAnsi(
        ReadOnlySpan<double> noise, int fs,
        string method, string speechLevel,
        double[]? threshold = null)
        => Mosqito.SqMetrics.SpeechIntelligibility.SiiAnsi.Compute(
            noise, fs, method, speechLevel, threshold);

    /// <inheritdoc cref="SiiAnsi.ComputeFromSpectrum"/>
    public static (double SII, double[] SIISpecific, double[] FreqAxis) SiiAnsiFreq(
        ReadOnlySpan<double> spectrum, ReadOnlySpan<double> freqs,
        string method, string speechLevel,
        double[]? threshold = null)
        => Mosqito.SqMetrics.SpeechIntelligibility.SiiAnsi.ComputeFromSpectrum(
            spectrum, freqs, method, speechLevel, threshold);

    /// <inheritdoc cref="SiiAnsi.ComputeFromLevel"/>
    public static (double SII, double[] SIISpecific, double[] FreqAxis) SiiAnsiLevel(
        double noiseLevel,
        string method, string speechLevel,
        double[]? threshold = null)
        => Mosqito.SqMetrics.SpeechIntelligibility.SiiAnsi.ComputeFromLevel(
            noiseLevel, method, speechLevel, threshold);
}
