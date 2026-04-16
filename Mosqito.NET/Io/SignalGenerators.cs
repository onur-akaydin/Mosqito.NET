namespace Mosqito.Io;

/// <summary>
/// Synthetic signal generators.
/// Matches MoSQITo's utils/sine_wave_generator.py, am_sine_generator.py,
/// am_noise_generator.py, fm_sine_generator.py.
/// </summary>
public static class SignalGenerators
{
    private const double P_Ref = 2e-5; // 20 µPa — ISO 1683 reference pressure

    // ------------------------------------------------------------------
    // Sine wave  (sine_wave_generator.py)
    // ------------------------------------------------------------------

    /// <summary>
    /// Generates a sinusoidal pressure signal at a given SPL.
    /// </summary>
    /// <param name="fs">Sampling frequency [Hz].</param>
    /// <param name="d">Signal duration [s].</param>
    /// <param name="freq">Sine frequency [Hz].</param>
    /// <param name="splLevel">Sound Pressure Level [dB ref 20 µPa].</param>
    /// <returns>(<c>signal</c> [Pa], <c>time</c> [s])</returns>
    public static (double[] signal, double[] time) SineWave(int fs, double d, double freq, double splLevel)
    {
        double pressureRms = P_Ref * Math.Pow(10.0, splLevel / 20.0);
        double amplitude   = Math.Sqrt(2.0) * pressureRms;

        int n = (int)(d * fs);
        double[] signal = new double[n];
        double[] time   = new double[n];
        double dt = 1.0 / fs;
        for (int i = 0; i < n; i++)
        {
            time[i] = i * dt;
            signal[i] = amplitude * Math.Sin(2.0 * Math.PI * freq * time[i]);
        }
        return (signal, time);
    }

    // ------------------------------------------------------------------
    // AM sine carrier  (am_sine_generator.py)
    // ------------------------------------------------------------------

    /// <summary>
    /// Amplitude-modulated signal with a sinusoidal carrier.
    /// </summary>
    /// <param name="xmod">Modulating signal (length N).</param>
    /// <param name="fs">Sampling frequency [Hz].</param>
    /// <param name="fc">Carrier frequency [Hz]. Must be &lt; fs/2.</param>
    /// <param name="splLevel">Target SPL [dB ref 20 µPa RMS].</param>
    /// <param name="printM">Print modulation index to console.</param>
    /// <returns>(<c>yAm</c> [Pa], modulation index m)</returns>
    public static (double[] yAm, double m) AmSine(
        ReadOnlySpan<double> xmod, double fs, double fc, double splLevel, bool printM = false)
    {
        if (fc >= fs / 2) throw new ArgumentException("Carrier frequency must be < fs/2.", nameof(fc));

        int nt = xmod.Length;
        double T = nt / fs, dt = 1.0 / fs;

        // Carrier
        double[] t = new double[nt];
        double[] xc = new double[nt];
        for (int i = 0; i < nt; i++) { t[i] = i * dt; xc[i] = Math.Sin(2.0 * Math.PI * fc * t[i]); }

        // AM signal
        double[] yAm = new double[nt];
        for (int i = 0; i < nt; i++) yAm[i] = (1.0 + xmod[i]) * xc[i];

        // Modulation index
        double m = 0.0;
        for (int i = 0; i < nt; i++) m = Math.Max(m, Math.Abs(xmod[i]));

        if (printM) Console.WriteLine($"AM Modulation index = {m}");
        if (m > 1.0) Console.WriteLine("Warning ['am_sine_generator']: modulation index m > 1\n\tSignal is overmodulated!");

        // Normalise to target SPL
        double aRms = 20e-6 * Math.Pow(10.0, splLevel / 20.0);
        double stdY = Std(yAm);
        if (stdY > 0)
            for (int i = 0; i < nt; i++) yAm[i] *= aRms / stdY;

        return (yAm, m);
    }

    // ------------------------------------------------------------------
    // AM broadband noise carrier  (am_noise_generator.py)
    // ------------------------------------------------------------------

    /// <summary>
    /// Amplitude-modulated signal with a Gaussian noise carrier.
    /// Uses a thread-local, optionally seeded RNG to allow reproducibility.
    /// </summary>
    /// <param name="xmod">Modulating signal (length N).</param>
    /// <param name="splLevel">Target SPL [dB ref 20 µPa RMS].</param>
    /// <param name="seed">Optional RNG seed for reproducibility.</param>
    /// <param name="printM">Print modulation index.</param>
    /// <returns>(<c>yAm</c> [Pa], modulation index m)</returns>
    public static (double[] yAm, double m) AmNoise(
        ReadOnlySpan<double> xmod, double splLevel, int? seed = null, bool printM = false)
    {
        int nt = xmod.Length;
        var rng = seed.HasValue ? new Random(seed.Value) : new Random();

        // Gaussian noise carrier (Box-Muller)
        double[] xc = GaussianNoise(nt, rng);

        double[] yAm = new double[nt];
        for (int i = 0; i < nt; i++) yAm[i] = (1.0 + xmod[i]) * xc[i];

        double m = 0.0;
        for (int i = 0; i < nt; i++) m = Math.Max(m, Math.Abs(xmod[i]));

        if (printM) Console.WriteLine($"AM Modulation index = {m}");
        if (m > 1.0) Console.WriteLine("Warning ['am_noise_generator']: modulation index m > 1\n\tSignal is overmodulated!");

        double aRms = 20e-6 * Math.Pow(10.0, splLevel / 20.0);
        double stdY = Std(yAm);
        if (stdY > 0)
            for (int i = 0; i < nt; i++) yAm[i] *= aRms / stdY;

        return (yAm, m);
    }

    // ------------------------------------------------------------------
    // FM sine  (fm_sine_generator.py)
    // ------------------------------------------------------------------

    /// <summary>
    /// Frequency-modulated signal with sinusoidal carrier.
    /// </summary>
    /// <param name="xmod">Modulating signal (length N).</param>
    /// <param name="fs">Sampling frequency [Hz].</param>
    /// <param name="fc">Carrier frequency [Hz]. Must be &lt; fs/2.</param>
    /// <param name="k">Frequency sensitivity [Hz per unit amplitude].</param>
    /// <param name="splLevel">Target SPL [dB ref 20 µPa RMS].</param>
    /// <param name="printInfo">Print frequency deviation and modulation index.</param>
    /// <returns>(<c>yFm</c> [Pa], instantaneous freq, max freq deviation, modulation index)</returns>
    public static (double[] yFm, double[] instFreq, double fDelta, double modIdx)
        FmSine(ReadOnlySpan<double> xmod, double fs, double fc, double k, double splLevel,
               bool printInfo = false)
    {
        if (fc >= fs / 2) throw new ArgumentException("Carrier frequency must be < fs/2.", nameof(fc));

        int nt = xmod.Length;
        double dt = 1.0 / fs;

        // Instantaneous frequency
        double[] instFreq = new double[nt];
        for (int i = 0; i < nt; i++) instFreq[i] = fc + k * xmod[i];

        // FM signal via cumulative phase
        double[] cumPhase = new double[nt];
        double acc = 0.0;
        for (int i = 0; i < nt; i++) { acc += instFreq[i]; cumPhase[i] = acc; }

        double[] yFm = new double[nt];
        for (int i = 0; i < nt; i++)
            yFm[i] = Math.Sin(2.0 * Math.PI * cumPhase[i] * dt);

        // Max frequency deviation and modulation index
        double maxAbs = 0.0;
        for (int i = 0; i < nt; i++) maxAbs = Math.Max(maxAbs, Math.Abs(xmod[i]));
        double fDelta = k * maxAbs;

        double cumMax = 0.0;
        acc = 0.0;
        for (int i = 0; i < nt; i++)
        {
            acc += xmod[i];
            cumMax = Math.Max(cumMax, Math.Abs(2.0 * Math.PI * k * acc * dt));
        }
        double modIdx = cumMax;

        // Normalise to target SPL
        double aRms = 20e-6 * Math.Pow(10.0, splLevel / 20.0);
        double stdY = Std(yFm);
        if (stdY > 0)
            for (int i = 0; i < nt; i++) yFm[i] *= aRms / stdY;

        if (printInfo)
        {
            Console.WriteLine($"\tMax freq deviation: {fDelta} Hz");
            Console.WriteLine($"\tFM modulation index: {modIdx:F2}");
        }

        return (yFm, instFreq, fDelta, modIdx);
    }

    // ------------------------------------------------------------------
    // Helpers
    // ------------------------------------------------------------------

    private static double Std(double[] x)
    {
        double mean = 0.0;
        for (int i = 0; i < x.Length; i++) mean += x[i];
        mean /= x.Length;
        double sumSq = 0.0;
        for (int i = 0; i < x.Length; i++) sumSq += (x[i] - mean) * (x[i] - mean);
        return Math.Sqrt(sumSq / x.Length);
    }

    private static double[] GaussianNoise(int n, Random rng)
    {
        double[] result = new double[n];
        for (int i = 0; i < n - 1; i += 2)
        {
            // Box-Muller
            double u1 = rng.NextDouble();
            double u2 = rng.NextDouble();
            if (u1 < 1e-300) u1 = 1e-300;
            double mag = Math.Sqrt(-2.0 * Math.Log(u1));
            result[i]     = mag * Math.Cos(2.0 * Math.PI * u2);
            result[i + 1] = mag * Math.Sin(2.0 * Math.PI * u2);
        }
        if (n % 2 != 0)
        {
            double u1 = rng.NextDouble();
            double u2 = rng.NextDouble();
            if (u1 < 1e-300) u1 = 1e-300;
            result[n - 1] = Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
        }
        return result;
    }
}
