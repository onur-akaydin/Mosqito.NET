namespace Mosqito.SoundLevelMeter;

/// <summary>
/// Computes exact and preferred nth-octave centre frequencies.
/// Ported from MoSQITo <c>_center_freq.py</c> (ANSI S1.1-1986).
/// </summary>
internal static class CenterFreq
{
    /// <summary>
    /// Computes exact and preferred nth-octave centre frequencies between fmin and fmax.
    /// </summary>
    /// <param name="fmin">Minimum frequency [Hz].</param>
    /// <param name="fmax">Maximum frequency [Hz].</param>
    /// <param name="n">Bands per octave (1 or 3).</param>
    /// <param name="G">Base system: 2 (binary) or 10 (decade).</param>
    /// <param name="fr">Reference frequency [Hz], default 1000 Hz.</param>
    /// <returns>(exactCentres, preferredCentres)</returns>
    internal static (double[] fExact, double[] fNom) Compute(
        double fmin, double fmax, int n = 3, int G = 10, double fr = 1000.0)
    {
        double b = 1.0 / n;
        double U = G switch {
            2  => Math.Pow(2.0, b),
            10 => Math.Pow(10.0, 3.0 * b / 10.0),
            _  => throw new ArgumentException("G must be 2 or 10.", nameof(G))
        };

        int kmin = (int)Math.Round(Math.Log(fmin / fr) / Math.Log(U), MidpointRounding.AwayFromZero);
        int kmax = (int)Math.Round(Math.Log(fmax / fr) / Math.Log(U), MidpointRounding.AwayFromZero);

        int count = kmax - kmin + 1;
        double[] fExact = new double[count];
        for (int i = 0; i < count; i++)
            fExact[i] = fr * Math.Pow(U, kmin + i);

        // Preferred (nominal) frequencies for n=1 or n=3
        double[] nominalTable = n == 1 ? NominalFrequencies.Octave : NominalFrequencies.ThirdOctave;
        double[] fNom = new double[count];

        if (n == 1 || n == 3)
        {
            // Find reference index in nominal table
            int iRef = Array.IndexOf(nominalTable, fr);
            if (iRef < 0) iRef = 0; // fallback

            for (int i = 0; i < count; i++)
            {
                int k = kmin + i;
                int nomIdx = iRef + k;
                if (nomIdx >= 0 && nomIdx < nominalTable.Length)
                    fNom[i] = nominalTable[nomIdx];
                else
                    fNom[i] = fExact[i];
            }
        }
        else
        {
            // For other n values, use exact
            Array.Copy(fExact, fNom, count);
        }

        return (fExact, fNom);
    }
}
