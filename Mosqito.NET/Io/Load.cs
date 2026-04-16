using Mosqito.Dsp;

namespace Mosqito.Io;

/// <summary>
/// Signal loader — reads .wav (and optionally .mat) files, auto-resamples to 48 kHz.
/// Matches MoSQITo <c>utils/load.py</c> behaviour exactly.
/// </summary>
public static class Load
{
    /// <summary>
    /// Loads a signal from a .wav file, optionally applies a Pa/FS calibration factor,
    /// and resamples to 48 kHz if necessary.
    /// </summary>
    /// <param name="file">Path to the .wav (or .mat) file.</param>
    /// <param name="wavCalib">
    /// Calibration factor [Pa/FS]. Null → 1.0 (prints a notice, matching Python behaviour).
    /// </param>
    /// <param name="matSignal">Variable name for the signal in .mat files (ignored for .wav).</param>
    /// <param name="matFs">Variable name for the sample rate in .mat files (ignored for .wav).</param>
    /// <returns>(<c>signal</c> in Pa, <c>fs</c> which is always 48000 after resampling).</returns>
    public static (double[] signal, int fs) FromFile(
        string file,
        double? wavCalib = null,
        string matSignal = "",
        string matFs = "")
    {
        string ext = Path.GetExtension(file).ToLowerInvariant();

        double[] signal;
        int fs;

        switch (ext)
        {
            case ".wav":
                (signal, fs) = WavLoader.Read(file, wavCalib);
                break;

            case ".mat":
                throw new NotSupportedException(
                    ".mat loading requires the optional MatFileHandler package. " +
                    "Reference Mosqito.Io.MatLoader and define MOSQITO_MAT.");

            case ".uff":
            case ".unv":
                throw new NotSupportedException(
                    ".uff/.unv loading is not supported in the runtime library. " +
                    "Use the Mosqito.NET.Tests UFF reader for test/validation use.");

            default:
                throw new ArgumentException(
                    $"ERROR: only .wav .mat or .uff files are supported. Got: {ext}", nameof(file));
        }

        // Resample to 48 kHz if needed (mandatory for ECMA-418-2 and Zwicker paths)
        if (fs != 48000)
        {
            signal = Resample.Apply(signal, fs, 48000);
            fs = 48000;
            Console.WriteLine("[Info] Signal resampled to 48 kHz to allow calculation.");
        }

        return (signal, fs);
    }

    /// <summary>
    /// Convenience overload: loads from a .wav file path only.
    /// </summary>
    public static (double[] signal, int fs) FromWav(string path, double? wavCalib = null)
        => FromFile(path, wavCalib);
}
