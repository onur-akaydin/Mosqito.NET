using System.Numerics;
using MathNet.Numerics.IntegralTransforms;

namespace Mosqito.Dsp;

/// <summary>
/// FFT-based rational resampling.
/// Matches scipy.signal.resample(x, num) — truncated-spectrum method.
/// </summary>
public static class Resample
{
    /// <summary>
    /// Resamples <paramref name="input"/> from its current length to <paramref name="nOut"/> samples.
    /// Uses the FFT truncation/zero-padding method matching scipy.signal.resample.
    /// </summary>
    public static double[] Apply(ReadOnlySpan<double> input, int nOut)
    {
        int nIn = input.Length;
        if (nOut == nIn)
        {
            double[] copy = new double[nIn];
            input.CopyTo(copy);
            return copy;
        }

        // Forward FFT
        Complex[] spec = Fft.Fft2(input);

        // Resample in frequency domain
        Complex[] specR = new Complex[nOut];
        if (nOut > nIn)
        {
            // Zero-pad: copy nIn FFT bins, pad rest with zeros, scale
            int half = nIn / 2;
            // Positive freqs
            for (int i = 0; i <= half; i++) specR[i] = spec[i];
            // Negative freqs → at the end
            int negStart = nOut - (nIn - half - 1);
            for (int i = half + 1; i < nIn; i++)
                specR[negStart + (i - half - 1)] = spec[i];
            // Scale
            double scale = (double)nOut / nIn;
            for (int i = 0; i < nOut; i++) specR[i] *= scale;
        }
        else
        {
            // Truncate: keep nOut FFT bins, scale
            int half = nOut / 2;
            for (int i = 0; i <= half; i++) specR[i] = spec[i];
            int negStart = nOut - (nIn / 2 - half - 1 + 1);
            for (int i = half + 1; i < nOut; i++)
            {
                int srcIdx = nIn - (nOut - i);
                if (srcIdx >= 0 && srcIdx < nIn) specR[i] = spec[srcIdx];
            }
            double scale = (double)nOut / nIn;
            for (int i = 0; i < nOut; i++) specR[i] *= scale;
        }

        // Inverse FFT
        Fourier.Inverse(specR, FourierOptions.AsymmetricScaling);
        double[] result = new double[nOut];
        for (int i = 0; i < nOut; i++) result[i] = specR[i].Real;
        return result;
    }

    /// <summary>
    /// Resamples <paramref name="input"/> from <paramref name="fsIn"/> Hz to
    /// <paramref name="fsOut"/> Hz.
    /// </summary>
    public static double[] Apply(ReadOnlySpan<double> input, int fsIn, int fsOut)
    {
        if (fsIn == fsOut)
        {
            double[] copy = new double[input.Length];
            input.CopyTo(copy);
            return copy;
        }
        int nOut = (int)Math.Round((double)input.Length * fsOut / fsIn);
        return Apply(input, nOut);
    }
}
