namespace Mosqito.SoundLevelMeter;

/// <summary>
/// ANSI S1.1-1986 nth-octave filter bandwidth computation.
/// Ported from MoSQITo <c>_filter_bandwidth.py</c>.
/// </summary>
internal static class FilterBandwidth
{
    /// <summary>
    /// Computes the bandwidth parameters for nth-octave Butterworth filters.
    /// </summary>
    /// <param name="fc">Centre frequencies [Hz].</param>
    /// <param name="n">Bands per octave.</param>
    /// <param name="filterOrder">Butterworth filter order (ANSI N).</param>
    /// <returns>(alpha, fLow, fHigh) arrays of same length as fc.</returns>
    internal static (double[] alpha, double[] fLow, double[] fHigh) Compute(
        double[] fc, int n = 3, int filterOrder = 3)
    {
        int count = fc.Length;
        double b = 1.0 / n;
        double[] alpha = new double[count];
        double[] fLow  = new double[count];
        double[] fHigh = new double[count];

        double pow2b = Math.Pow(2.0, b / 2.0);

        for (int i = 0; i < count; i++)
        {
            double f1 = fc[i] / pow2b;         // ANSI eq5
            double f2 = fc[i] * pow2b;         // ANSI eq6
            double qr = fc[i] / (f2 - f1);     // ANSI eq7&8
            double qd = (Math.PI / 2.0 / filterOrder) /
                        Math.Sin(Math.PI / 2.0 / filterOrder) * qr; // ANSI eq9
            double a  = (1.0 + Math.Sqrt(1.0 + 4.0 * qd * qd)) / (2.0 * qd);

            alpha[i] = a;
            fLow[i]  = f1;
            fHigh[i] = f2;
        }

        return (alpha, fLow, fHigh);
    }
}
