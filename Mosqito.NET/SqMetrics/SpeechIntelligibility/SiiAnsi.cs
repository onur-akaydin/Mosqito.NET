using Mosqito.Conversion;
using Mosqito.SoundLevelMeter;
using Mosqito.SqMetrics.Loudness;

namespace Mosqito.SqMetrics.SpeechIntelligibility;

/// <summary>
/// Speech Intelligibility Index (ANSI S3.5).
/// Ported from MoSQITo <c>sii_ansi.py</c>, <c>sii_ansi_freq.py</c>,
/// <c>sii_ansi_level.py</c>, <c>_main_sii.py</c>.
/// </summary>
public static class SiiAnsi
{
    // -----------------------------------------------------------------------
    // Public entry points
    // -----------------------------------------------------------------------

    /// <summary>Computes SII from a noise time signal.</summary>
    /// <param name="noise">Noise time signal [Pa].</param>
    /// <param name="fs">Sampling frequency [Hz].</param>
    /// <param name="method">Band method: "critical", "equally_critical", "third_octave", "octave".</param>
    /// <param name="speechLevel">Speech level: "normal", "raised", "loud", "shout".</param>
    /// <param name="threshold">Hearing threshold per band [dB], or null for zeros, or "zwicker".</param>
    public static (double SII, double[] SIISpecific, double[] FreqAxis) Compute(
        ReadOnlySpan<double> noise, int fs,
        string method, string speechLevel,
        double[]? threshold = null)
    {
        ValidateMethod(method);
        ValidateSpeechLevel(speechLevel);

        var (speechSpectrum, _) = GetSpeechData(method, speechLevel);
        var bandData = GetBandData(method);

        // FFT + band synthesis
        var (spec, freqs) = CompSpectrum.Compute(noise, fs,
            window: CompSpectrum.WindowType.Blackman, db: true);
        var (noiseSpectrum, _) = FreqBandSynthesis.Compute(
            spec, freqs, bandData.LowerFreqs, bandData.UpperFreqs);

        return MainSii(method, speechSpectrum, noiseSpectrum, threshold, bandData);
    }

    /// <summary>Computes SII from a noise dB spectrum.</summary>
    public static (double SII, double[] SIISpecific, double[] FreqAxis) ComputeFromSpectrum(
        ReadOnlySpan<double> spectrum, ReadOnlySpan<double> freqs,
        string method, string speechLevel,
        double[]? threshold = null)
    {
        ValidateMethod(method);
        ValidateSpeechLevel(speechLevel);

        var (speechSpectrum, _) = GetSpeechData(method, speechLevel);
        var bandData = GetBandData(method);

        // Re-synthesise into band levels if needed
        double[] noiseSpectrum;
        if (spectrum.Length == bandData.CenterFreqs.Length)
        {
            noiseSpectrum = spectrum.ToArray();
        }
        else
        {
            (noiseSpectrum, _) = FreqBandSynthesis.Compute(
                spectrum, freqs, bandData.LowerFreqs, bandData.UpperFreqs);
        }

        return MainSii(method, speechSpectrum, noiseSpectrum, threshold, bandData);
    }

    /// <summary>Computes SII from an overall noise level in dB (uniform spectrum assumed).</summary>
    public static (double SII, double[] SIISpecific, double[] FreqAxis) ComputeFromLevel(
        double noiseLevel,
        string method, string speechLevel,
        double[]? threshold = null)
    {
        ValidateMethod(method);
        ValidateSpeechLevel(speechLevel);

        var (speechSpectrum, _) = GetSpeechData(method, speechLevel);
        var bandData = GetBandData(method);
        int nBands = bandData.CenterFreqs.Length;

        double bandLevel = 10.0 * Math.Log10(Math.Pow(10.0, noiseLevel / 10.0) / nBands);
        double[] noiseSpectrum = new double[nBands];
        for (int i = 0; i < nBands; i++) noiseSpectrum[i] = bandLevel;

        return MainSii(method, speechSpectrum, noiseSpectrum, threshold, bandData);
    }

    // -----------------------------------------------------------------------
    // Core algorithm  (_main_sii)
    // -----------------------------------------------------------------------

    /// <summary>
    /// Direct port of <c>_main_sii(method, speech_spectrum, noise_spectrum, threshold)</c>.
    /// Exposes the core algorithm with both speech and noise band spectra supplied
    /// explicitly (bypassing the level-derived speech tables).
    /// </summary>
    internal static (double SII, double[] SIISpecific, double[] FreqAxis) MainSiiForTest(
        string method, double[] speechSpectrum, double[] noiseSpectrum,
        double[]? threshold = null)
    {
        ValidateMethod(method);
        return MainSii(method, speechSpectrum, noiseSpectrum, threshold, GetBandData(method));
    }

    private static (double SII, double[] SIISpecific, double[] FreqAxis) MainSii(
        string method, double[] speechSpectrum, double[] noiseSpectrum,
        double[]? threshold, BandData bd)
    {
        int nBands = bd.CenterFreqs.Length;

        // Threshold of hearing T
        double[] T;
        if (threshold == null)
        {
            T = new double[nBands];
        }
        else
        {
            T = threshold;
        }

        // STEP 3 — effective noise Z (upward spread of masking)
        double[] Z = new double[nBands];
        if (method == "octave")
        {
            // For octave, no spreading needed
            for (int i = 0; i < nBands; i++) Z[i] = noiseSpectrum[i];
        }
        else
        {
            double[] V = new double[nBands];
            double[] B = new double[nBands];
            for (int i = 0; i < nBands; i++)
            {
                V[i] = speechSpectrum[i] - 24.0;
                B[i] = Math.Max(noiseSpectrum[i], V[i]);
            }

            double[] C = new double[nBands];
            if (method == "third_octave")
            {
                for (int i = 0; i < nBands; i++)
                    C[i] = -80.0 + 0.6 * (B[i] + 10.0 * Math.Log10(bd.CenterFreqs[i]) - 6.353);

                for (int i = 0; i < nBands; i++)
                {
                    double s = 0.0;
                    for (int k = 0; k < i; k++)
                        s += Math.Pow(10.0, 0.1 * (B[k] + 3.32 * C[k] *
                            Math.Log10(0.89 * bd.CenterFreqs[i] / bd.CenterFreqs[k])));
                    Z[i] = 10.0 * Math.Log10(Math.Pow(10.0, 0.1 * noiseSpectrum[i]) + s);
                }
            }
            else
            {
                // critical / equally_critical
                for (int i = 0; i < nBands; i++)
                    C[i] = -80.0 + 0.6 * (B[i] + 10.0 * Math.Log10(bd.UpperFreqs[i] - bd.LowerFreqs[i]));

                for (int i = 0; i < nBands; i++)
                {
                    double s = 0.0;
                    for (int k = 0; k < i - 1; k++)
                        s += Math.Pow(10.0, 0.1 * (B[k] + 3.32 * C[k] *
                            Math.Log10(bd.CenterFreqs[i] / bd.CenterFreqs[k])));
                    Z[i] = 10.0 * Math.Log10(Math.Pow(10.0, 0.1 * noiseSpectrum[i]) + s);
                }
            }
            Z[0] = B[0]; // 4.3.2.4
        }

        // STEP 4
        double[] X = new double[nBands];
        for (int i = 0; i < nBands; i++) X[i] = bd.RefInternalNoise[i] + T[i];

        // STEP 5
        double[] D = new double[nBands];
        for (int i = 0; i < nBands; i++) D[i] = Math.Max(Z[i], X[i]);

        // STEP 6
        double[] L = new double[nBands];
        for (int i = 0; i < nBands; i++)
        {
            double v = 1.0 - (speechSpectrum[i] - bd.StdSpeechSpectrumNormal[i] - 10.0) / 160.0;
            L[i] = Math.Min(v, 1.0);
        }

        // STEP 7
        double[] K = new double[nBands];
        for (int i = 0; i < nBands; i++)
        {
            double v = (speechSpectrum[i] - D[i] + 15.0) / 30.0;
            K[i] = Math.Max(0.0, Math.Min(1.0, v));
        }

        double[] A = new double[nBands];
        for (int i = 0; i < nBands; i++) A[i] = L[i] * K[i];

        // STEP 8
        double sii = 0.0;
        double[] siiSpec = new double[nBands];
        for (int i = 0; i < nBands; i++)
        {
            siiSpec[i] = bd.Importance[i] * A[i];
            sii += siiSpec[i];
        }

        return (sii, siiSpec, bd.CenterFreqs);
    }

    // -----------------------------------------------------------------------
    // Band and speech data
    // -----------------------------------------------------------------------

    private readonly struct BandData
    {
        public readonly double[] CenterFreqs;
        public readonly double[] LowerFreqs;
        public readonly double[] UpperFreqs;
        public readonly double[] Importance;
        public readonly double[] RefInternalNoise;
        public readonly double[] StdSpeechSpectrumNormal;

        public BandData(double[] cf, double[] lf, double[] uf,
                         double[] imp, double[] rin, double[] ssn)
        {
            CenterFreqs = cf; LowerFreqs = lf; UpperFreqs = uf;
            Importance = imp; RefInternalNoise = rin; StdSpeechSpectrumNormal = ssn;
        }
    }

    private static BandData GetBandData(string method) => method switch
    {
        "critical"          => CriticalBandData(),
        "equally_critical"  => EqualCriticalBandData(),
        "third_octave"      => ThirdOctaveBandData(),
        "octave"            => OctaveBandData(),
        _                   => throw new ArgumentException($"Unknown method: {method}"),
    };

    private static (double[] speechSpectrum, double speechLevel) GetSpeechData(
        string method, string speechLevel) => method switch
    {
        "critical"          => CriticalBandSpeechData(speechLevel),
        "equally_critical"  => EqualCriticalBandSpeechData(speechLevel),
        "third_octave"      => ThirdOctaveBandSpeechData(speechLevel),
        "octave"            => OctaveBandSpeechData(speechLevel),
        _                   => throw new ArgumentException($"Unknown method: {method}"),
    };

    private static void ValidateMethod(string m)
    {
        if (m != "critical" && m != "equally_critical" && m != "third_octave" && m != "octave")
            throw new ArgumentException($"Method must be one of critical, equally_critical, third_octave, octave. Got: {m}");
    }

    private static void ValidateSpeechLevel(string s)
    {
        if (s != "normal" && s != "raised" && s != "loud" && s != "shout")
            throw new ArgumentException($"Speech level must be one of normal, raised, loud, shout. Got: {s}");
    }

    // -----------------------------------------------------------------------
    // Band procedure data  (_band_procedure_data.py)
    // -----------------------------------------------------------------------

    private static BandData CriticalBandData() => new(
        cf:  new[] { 150d, 250, 350, 450, 570, 700, 840, 1000, 1170, 1370, 1600, 1850, 2150, 2500, 2900, 3400, 4000, 4800, 5800, 7000, 8500 },
        lf:  new[] { 100d, 200, 300, 400, 510, 630, 770, 920, 1080, 1270, 1480, 1720, 2000, 2320, 2700, 3150, 3700, 4400, 5300, 6400, 7700 },
        uf:  new[] { 200d, 300, 400, 510, 630, 770, 920, 1080, 1270, 1480, 1720, 2000, 2320, 2700, 3150, 3700, 4400, 5300, 6400, 7700, 9500 },
        imp: new[] { 0.0103, 0.0261, 0.0419, 0.0577, 0.0577, 0.0577, 0.0577, 0.0577, 0.0577, 0.0577, 0.0577, 0.0577, 0.0577, 0.0577, 0.0577, 0.0577, 0.0577, 0.0460, 0.0343, 0.0226, 0.0110 },
        rin: new[] { 1.50, -3.90, -7.20, -8.90, -10.30, -11.40, -12.00, -12.50, -13.20, -14.00, -15.40, -16.90, -18.80, -21.20, -23.20, -24.90, -25.90, -24.20, -19.00, -11.70, -6.00 },
        ssn: new[] { 31.44, 34.75, 34.14, 34.58, 33.17, 30.64, 27.59, 25.01, 23.52, 22.28, 20.15, 18.29, 16.37, 13.80, 12.21, 11.09, 9.33, 5.84, 3.47, 1.78, -0.14 });

    private static BandData EqualCriticalBandData() => new(
        cf:  new[] { 350d, 450, 570, 700, 840, 1000, 1170, 1370, 1600, 1850, 2150, 2500, 2900, 3400, 4000, 4800, 5800 },
        lf:  new[] { 300d, 400, 510, 630, 770, 920, 1080, 1270, 1480, 1720, 2000, 2320, 2700, 3150, 3700, 4400, 5300 },
        uf:  new[] { 400d, 510, 630, 770, 920, 1080, 1270, 1480, 1720, 2000, 2320, 2700, 3150, 3700, 4400, 5300, 6400 },
        imp: new double[17].Select(_ => 0.0588).ToArray(),
        rin: new[] { -7.20, -8.90, -10.30, -11.40, -12.00, -12.50, -13.20, -14.00, -15.40, -16.90, -18.80, -21.20, -23.20, -24.90, -25.90, -24.20, -19.00 },
        ssn: new[] { 34.14, 34.58, 33.17, 30.64, 27.59, 25.01, 23.52, 22.28, 20.15, 18.29, 16.37, 13.80, 12.21, 11.09, 9.33, 5.84, 3.47 });

    private static BandData ThirdOctaveBandData() => new(
        cf:  new[] { 160d, 200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000 },
        lf:  new[] { 141d, 178, 224, 282, 355, 447, 562, 708, 891, 1122, 1413, 1778, 2239, 2818, 3548, 4467, 5623, 7079 },
        uf:  new[] { 178d, 224, 282, 355, 447, 562, 708, 891, 1122, 1413, 1778, 2239, 2818, 3548, 4467, 5623, 7079, 8913 },
        imp: new[] { 0.0083, 0.0095, 0.0150, 0.0289, 0.0440, 0.0578, 0.0653, 0.0711, 0.0818, 0.0844, 0.0882, 0.0898, 0.0868, 0.0844, 0.0771, 0.0527, 0.0364, 0.0185 },
        rin: new[] { 0.60, -1.70, -3.90, -6.10, -8.20, -9.70, -10.80, -11.90, -12.50, -13.50, -15.40, -17.70, -21.20, -24.20, -25.90, -23.60, -15.80, -7.10 },
        ssn: new[] { 32.41, 34.48, 34.75, 33.98, 34.59, 34.27, 32.06, 28.30, 25.01, 23.00, 20.15, 17.32, 13.18, 11.55, 9.33, 5.31, 2.59, 1.13 });

    private static BandData OctaveBandData() => new(
        cf:  new[] { 250d, 500, 1000, 2000, 4000, 8000 },
        lf:  new[] { 177d, 355, 710, 1420, 2840, 5680 },
        uf:  new[] { 355d, 710, 1420, 2840, 5680, 11360 },
        imp: new[] { 0.0617, 0.1671, 0.2373, 0.2648, 0.2142, 0.0549 },
        rin: new[] { -3.90, -9.70, -12.50, -17.70, -25.90, -7.10 },
        ssn: new[] { 34.75, 34.27, 25.01, 17.32, 9.33, 1.13 });

    // -----------------------------------------------------------------------
    // Speech data  (_speech_data.py)
    // -----------------------------------------------------------------------

    private static (double[] spectrum, double level) CriticalBandSpeechData(string speechLevel) => speechLevel switch
    {
        "normal" => (new[] { 31.44, 34.75, 34.14, 34.58, 33.17, 30.64, 27.59, 25.01, 23.52, 22.28, 20.15, 18.29, 16.37, 13.80, 12.21, 11.09, 9.33, 5.84, 3.47, 1.78, -0.14 }, 62.35),
        "raised" => (new[] { 34.06, 38.98, 38.62, 39.84, 39.44, 37.99, 35.85, 33.86, 32.56, 30.91, 28.58, 26.37, 24.34, 22.35, 21.04, 19.56, 16.78, 12.14, 9.04, 6.36, 3.44 }, 68.34),
        "loud"   => (new[] { 34.21, 41.55, 43.68, 44.08, 45.34, 45.22, 43.60, 42.16, 41.07, 39.68, 37.70, 35.62, 33.17, 30.98, 29.01, 27.71, 25.41, 19.20, 15.37, 12.61, 9.62 }, 74.85),
        "shout"  => (new[] { 28.69, 42.50, 47.14, 48.46, 50.17, 51.68, 51.43, 51.31, 49.40, 49.03, 47.65, 45.47, 43.13, 40.80, 39.15, 37.30, 34.41, 29.01, 25.17, 22.08, 18.76 }, 82.30),
        _        => throw new ArgumentException($"Unknown speech level: {speechLevel}"),
    };

    private static (double[] spectrum, double level) EqualCriticalBandSpeechData(string speechLevel) => speechLevel switch
    {
        "normal" => (new[] { 34.14, 34.58, 33.17, 30.64, 27.59, 25.01, 23.52, 22.28, 20.15, 18.29, 16.37, 13.80, 12.21, 11.09, 9.33, 5.84, 3.47 }, 62.35),
        "raised" => (new[] { 38.62, 39.84, 39.44, 37.99, 35.85, 33.86, 32.56, 30.91, 28.58, 26.37, 24.34, 22.35, 21.04, 19.56, 16.78, 12.14, 9.04 }, 68.34),
        "loud"   => (new[] { 43.68, 44.08, 45.34, 45.22, 43.60, 42.16, 41.07, 39.68, 37.70, 35.62, 33.17, 30.98, 29.01, 27.71, 25.41, 19.20, 15.37 }, 74.85),
        "shout"  => (new[] { 47.14, 48.46, 50.17, 51.68, 51.43, 51.31, 49.40, 49.03, 47.65, 45.47, 43.13, 40.80, 39.15, 37.30, 34.41, 29.01, 25.17 }, 82.30),
        _        => throw new ArgumentException($"Unknown speech level: {speechLevel}"),
    };

    private static (double[] spectrum, double level) ThirdOctaveBandSpeechData(string speechLevel) => speechLevel switch
    {
        "normal" => (new[] { 32.41, 34.48, 34.75, 33.98, 34.59, 34.27, 32.06, 28.30, 25.01, 23.00, 20.15, 17.32, 13.18, 11.55, 9.33, 5.31, 2.59, 1.13 }, 62.35),
        "raised" => (new[] { 33.81, 33.92, 38.98, 38.57, 39.11, 40.15, 38.78, 36.37, 33.86, 31.89, 28.58, 25.32, 22.35, 20.15, 16.78, 11.47, 7.67, 5.07 }, 68.34),
        "loud"   => (new[] { 35.29, 37.76, 41.55, 43.78, 43.40, 44.85, 45.55, 44.05, 42.16, 40.53, 37.70, 34.39, 30.98, 28.21, 25.41, 18.35, 13.87, 11.39 }, 74.85),
        "shout"  => (new[] { 30.77, 36.65, 42.50, 46.51, 47.40, 49.24, 51.21, 51.44, 51.31, 49.63, 47.65, 44.32, 40.80, 38.13, 34.41, 28.24, 23.45, 20.72 }, 82.30),
        _        => throw new ArgumentException($"Unknown speech level: {speechLevel}"),
    };

    private static (double[] spectrum, double level) OctaveBandSpeechData(string speechLevel) => speechLevel switch
    {
        "normal" => (new[] { 34.75, 34.27, 25.01, 17.32, 9.33, 1.13 }, 62.35),
        "raised" => (new[] { 38.98, 40.15, 33.86, 25.32, 16.78, 5.07 }, 68.34),
        "loud"   => (new[] { 41.55, 44.85, 42.16, 34.39, 25.41, 11.39 }, 74.85),
        "shout"  => (new[] { 42.50, 49.24, 51.31, 44.32, 34.41, 20.72 }, 82.30),
        _        => throw new ArgumentException($"Unknown speech level: {speechLevel}"),
    };
}
