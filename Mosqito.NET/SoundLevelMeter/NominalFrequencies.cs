namespace Mosqito.SoundLevelMeter;

/// <summary>
/// IEC 61260 / ANSI S1.1 nominal octave and 1/3-octave centre frequencies.
/// Ported from MoSQITo <c>_nominal_frequency.py</c>.
/// </summary>
internal static class NominalFrequencies
{
    internal static readonly double[] Octave = {
        31.5, 63.0, 125.0, 250.0, 500.0, 1000.0, 2000.0, 4000.0, 8000.0, 16000.0
    };

    internal static readonly double[] ThirdOctave = {
        25.0, 31.5, 40.0, 50.0, 63.0, 80.0, 100.0, 125.0, 160.0, 200.0,
        250.0, 315.0, 400.0, 500.0, 630.0, 800.0, 1000.0, 1250.0, 1600.0, 2000.0,
        2500.0, 3150.0, 4000.0, 5000.0, 6300.0, 8000.0, 10000.0, 12500.0, 16000.0, 20000.0
    };
}
