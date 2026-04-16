using NAudio.Wave;

namespace Mosqito.Io;

/// <summary>
/// WAV file reader backed by NAudio.Core.
/// Returns floating-point samples in Pa when a calibration factor is supplied.
/// Matches the scipy.io.wavfile + calibration logic of MoSQITo's load.py.
/// </summary>
internal static class WavLoader
{
    /// <summary>
    /// Reads a WAV file and returns the first channel as a <c>double[]</c>.
    /// </summary>
    /// <param name="path">Path to the .wav file.</param>
    /// <param name="wavCalib">
    /// Calibration factor in Pa/FS.  If null, 1.0 is used (no calibration).
    /// For 16-bit PCM: signal = raw_int16 / 32767.0 * wavCalib.
    /// For 32-bit PCM: signal = raw_int32 / 2147483647.0 * wavCalib.
    /// For float WAV: signal = raw_float * wavCalib.
    /// </param>
    /// <returns>(<c>signal</c>, <c>sampleRate</c>)</returns>
    public static (double[] signal, int sampleRate) Read(string path, double? wavCalib = null)
    {
        double calib = wavCalib ?? 1.0;

        using var reader = new AudioFileReader(path);
        int fs = reader.WaveFormat.SampleRate;
        int channels = reader.WaveFormat.Channels;
        long sampleFrames = reader.Length / (reader.WaveFormat.BitsPerSample / 8 * channels);

        float[] buffer = new float[sampleFrames * channels];
        int read = reader.Read(buffer, 0, buffer.Length);
        int nSamples = read / channels;

        double[] signal = new double[nSamples];

        // AudioFileReader normalises everything to float in [-1, 1]
        for (int i = 0; i < nSamples; i++)
            signal[i] = buffer[i * channels] * calib; // first channel only

        if (channels > 1)
            Console.WriteLine("[Info] Multichannel signal loaded. Keeping only first channel.");

        return (signal, fs);
    }
}
