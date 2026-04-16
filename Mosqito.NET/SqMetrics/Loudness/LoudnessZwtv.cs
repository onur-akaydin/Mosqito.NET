using Mosqito.Dsp;

namespace Mosqito.SqMetrics.Loudness;

/// <summary>
/// Zwicker time-varying loudness (ISO 532-1:2017 Method 1, time-varying path).
/// Ported from MoSQITo <c>loudness_zwtv.py</c> and its helper modules.
///
/// Pipeline (at 48 kHz input):
///   1. ISO 532-1 Table A.1/A.2 third-octave filter bank → level vs. time (2 kHz).
///   2. MainLoudness per time frame → core loudness (21, nTime).
///   3. Nonlinear temporal decay (_nl_loudness).
///   4. CalcSlopes per time frame → N(t), N'(z, t).
///   5. Temporal weighting → filtered N(t).
///   6. Decimate by 4 → 500 Hz output rate (2 ms resolution).
/// </summary>
public static class LoudnessZwtv
{
    // ------------------------------------------------------------------
    // ISO 532-1 Table A.1 — reference SOS sections (3 sections × 6 coefficients)
    // Format: [b0, b1, b2, a0, a1, a2]
    // ------------------------------------------------------------------
    private static readonly double[,] RefSos =
    {
        { 1,  2,  1,  1, -2, 1 },
        { 1,  0, -1,  1, -2, 1 },
        { 1, -2,  1,  1, -2, 1 }
    };

    // ------------------------------------------------------------------
    // ISO 532-1 Table A.2 — per-band a1/a2 corrections (28 bands × 3 sections × 2 coeffs)
    // Stored as flat (28, 3, 2): [band][section][coeff_index], coeff_index 0=a1 1=a2
    // ------------------------------------------------------------------
    // ReSharper disable once InconsistentNaming
    private static readonly double[,,] DiffA = BuildDiffA();

    private static double[,,] BuildDiffA()
    {
        // diff values from _third_octave_levels.py: third_octave_filter[band][section][4] and [5]
        // positions [4]=a1 delta, [5]=a2 delta
        double[,,] d = new double[28, 3, 2];
        // Band 0
        d[0, 0, 0] = -6.70260e-004; d[0, 0, 1] = 6.59453e-004;
        d[0, 1, 0] = -3.75071e-004; d[0, 1, 1] = 3.61926e-004;
        d[0, 2, 0] = -3.06523e-004; d[0, 2, 1] = 2.97634e-004;
        // Band 1
        d[1, 0, 0] = -8.47258e-004; d[1, 0, 1] = 8.30131e-004;
        d[1, 1, 0] = -4.76448e-004; d[1, 1, 1] = 4.55616e-004;
        d[1, 2, 0] = -3.88773e-004; d[1, 2, 1] = 3.74685e-004;
        // Band 2
        d[2, 0, 0] = -1.07210e-003; d[2, 0, 1] = 1.04496e-003;
        d[2, 1, 0] = -6.06567e-004; d[2, 1, 1] = 5.73553e-004;
        d[2, 2, 0] = -4.94004e-004; d[2, 2, 1] = 4.71677e-004;
        // Band 3
        d[3, 0, 0] = -1.35836e-003; d[3, 0, 1] = 1.31535e-003;
        d[3, 1, 0] = -7.74327e-004; d[3, 1, 1] = 7.22007e-004;
        d[3, 2, 0] = -6.29154e-004; d[3, 2, 1] = 5.93771e-004;
        // Band 4
        d[4, 0, 0] = -1.72380e-003; d[4, 0, 1] = 1.65564e-003;
        d[4, 1, 0] = -9.91780e-004; d[4, 1, 1] = 9.08866e-004;
        d[4, 2, 0] = -8.03529e-004; d[4, 2, 1] = 7.47455e-004;
        // Band 5
        d[5, 0, 0] = -2.19188e-003; d[5, 0, 1] = 2.08388e-003;
        d[5, 1, 0] = -1.27545e-003; d[5, 1, 1] = 1.14406e-003;
        d[5, 2, 0] = -1.02976e-003; d[5, 2, 1] = 9.40900e-004;
        // Band 6
        d[6, 0, 0] = -2.79386e-003; d[6, 0, 1] = 2.62274e-003;
        d[6, 1, 0] = -1.64828e-003; d[6, 1, 1] = 1.44006e-003;
        d[6, 2, 0] = -1.32520e-003; d[6, 2, 1] = 1.18438e-003;
        // Band 7
        d[7, 0, 0] = -3.57182e-003; d[7, 0, 1] = 3.30071e-003;
        d[7, 1, 0] = -2.14252e-003; d[7, 1, 1] = 1.81258e-003;
        d[7, 2, 0] = -1.71397e-003; d[7, 2, 1] = 1.49082e-003;
        // Band 8
        d[8, 0, 0] = -4.58305e-003; d[8, 0, 1] = 4.15355e-003;
        d[8, 1, 0] = -2.80413e-003; d[8, 1, 1] = 2.28135e-003;
        d[8, 2, 0] = -2.23006e-003; d[8, 2, 1] = 1.87646e-003;
        // Band 9
        d[9, 0, 0] = -5.90655e-003; d[9, 0, 1] = 5.22622e-003;
        d[9, 1, 0] = -3.69947e-003; d[9, 1, 1] = 2.87118e-003;
        d[9, 2, 0] = -2.92205e-003; d[9, 2, 1] = 2.36178e-003;
        // Band 10
        d[10, 0, 0] = -7.65243e-003; d[10, 0, 1] = 6.57493e-003;
        d[10, 1, 0] = -4.92540e-003; d[10, 1, 1] = 3.61318e-003;
        d[10, 2, 0] = -3.86007e-003; d[10, 2, 1] = 2.97240e-003;
        // Band 11
        d[11, 0, 0] = -1.00023e-002; d[11, 0, 1] = 8.29610e-003;
        d[11, 1, 0] = -6.63788e-003; d[11, 1, 1] = 4.55999e-003;
        d[11, 2, 0] = -5.15982e-003; d[11, 2, 1] = 3.75306e-003;
        // Band 12
        d[12, 0, 0] = -1.31230e-002; d[12, 0, 1] = 1.04220e-002;
        d[12, 1, 0] = -9.02274e-003; d[12, 1, 1] = 5.73132e-003;
        d[12, 2, 0] = -6.94543e-003; d[12, 2, 1] = 4.71734e-003;
        // Band 13
        d[13, 0, 0] = -1.73693e-002; d[13, 0, 1] = 1.30947e-002;
        d[13, 1, 0] = -1.24176e-002; d[13, 1, 1] = 7.20526e-003;
        d[13, 2, 0] = -9.46002e-003; d[13, 2, 1] = 5.93145e-003;
        // Band 14
        d[14, 0, 0] = -2.31934e-002; d[14, 0, 1] = 1.64308e-002;
        d[14, 1, 0] = -1.73009e-002; d[14, 1, 1] = 9.04761e-003;
        d[14, 2, 0] = -1.30358e-002; d[14, 2, 1] = 7.44926e-003;
        // Band 15
        d[15, 0, 0] = -3.13292e-002; d[15, 0, 1] = 2.06370e-002;
        d[15, 1, 0] = -2.44342e-002; d[15, 1, 1] = 1.13731e-002;
        d[15, 2, 0] = -1.82108e-002; d[15, 2, 1] = 9.36778e-003;
        // Band 16
        d[16, 0, 0] = -4.28261e-002; d[16, 0, 1] = 2.59325e-002;
        d[16, 1, 0] = -3.49619e-002; d[16, 1, 1] = 1.43046e-002;
        d[16, 2, 0] = -2.57855e-002; d[16, 2, 1] = 1.17912e-002;
        // Band 17
        d[17, 0, 0] = -5.91733e-002; d[17, 0, 1] = 3.25054e-002;
        d[17, 1, 0] = -5.06072e-002; d[17, 1, 1] = 1.79513e-002;
        d[17, 2, 0] = -3.69401e-002; d[17, 2, 1] = 1.48094e-002;
        // Band 18
        d[18, 0, 0] = -8.26348e-002; d[18, 0, 1] = 4.05894e-002;
        d[18, 1, 0] = -7.40348e-002; d[18, 1, 1] = 2.24476e-002;
        d[18, 2, 0] = -5.34977e-002; d[18, 2, 1] = 1.85371e-002;
        // Band 19
        d[19, 0, 0] = -1.17018e-001; d[19, 0, 1] = 5.08116e-002;
        d[19, 1, 0] = -1.09516e-001; d[19, 1, 1] = 2.81387e-002;
        d[19, 2, 0] = -7.85097e-002; d[19, 2, 1] = 2.32872e-002;
        // Band 20
        d[20, 0, 0] = -1.67714e-001; d[20, 0, 1] = 6.37872e-002;
        d[20, 1, 0] = -1.63378e-001; d[20, 1, 1] = 3.53729e-002;
        d[20, 2, 0] = -1.16419e-001; d[20, 2, 1] = 2.93723e-002;
        // Band 21
        d[21, 0, 0] = -2.42528e-001; d[21, 0, 1] = 7.98576e-002;
        d[21, 1, 0] = -2.45161e-001; d[21, 1, 1] = 4.43370e-002;
        d[21, 2, 0] = -1.73972e-001; d[21, 2, 1] = 3.70015e-002;
        // Band 22
        d[22, 0, 0] = -3.53142e-001; d[22, 0, 1] = 9.96330e-002;
        d[22, 1, 0] = -3.69163e-001; d[22, 1, 1] = 5.53535e-002;
        d[22, 2, 0] = -2.61399e-001; d[22, 2, 1] = 4.65428e-002;
        // Band 23
        d[23, 0, 0] = -5.16316e-001; d[23, 0, 1] = 1.24177e-001;
        d[23, 1, 0] = -5.55473e-001; d[23, 1, 1] = 6.89403e-002;
        d[23, 2, 0] = -3.93998e-001; d[23, 2, 1] = 5.86715e-002;
        // Band 24
        d[24, 0, 0] = -7.56635e-001; d[24, 0, 1] = 1.55023e-001;
        d[24, 1, 0] = -8.34281e-001; d[24, 1, 1] = 8.58123e-002;
        d[24, 2, 0] = -5.94547e-001; d[24, 2, 1] = 7.43960e-002;
        // Band 25
        d[25, 0, 0] = -1.10165e+000; d[25, 0, 1] = 1.91713e-001;
        d[25, 1, 0] = -1.23939e+000; d[25, 1, 1] = 1.05243e-001;
        d[25, 2, 0] = -8.91666e-001; d[25, 2, 1] = 9.40354e-002;
        // Band 26
        d[26, 0, 0] = -1.58477e+000; d[26, 0, 1] = 2.39049e-001;
        d[26, 1, 0] = -1.80505e+000; d[26, 1, 1] = 1.28794e-001;
        d[26, 2, 0] = -1.32500e+000; d[26, 2, 1] = 1.21333e-001;
        // Band 27
        d[27, 0, 0] = -2.50630e+000; d[27, 0, 1] = 1.42308e-001;
        d[27, 1, 0] = -2.19464e+000; d[27, 1, 1] = 2.76470e-001;
        d[27, 2, 0] = -1.90231e+000; d[27, 2, 1] = 1.47304e-001;
        return d;
    }

    // ISO 532-1 Table A.2 — filter gain values (28 bands)
    private static readonly double[] FilterGain =
    {
        4.30764e-011, 8.59340e-011, 1.71424e-010, 3.41944e-010, 6.82035e-010,
        1.36026e-009, 2.71261e-009, 5.40870e-009, 1.07826e-008, 2.14910e-008,
        4.28228e-008, 8.54316e-008, 1.70009e-007, 3.38215e-007, 6.71990e-007,
        1.33531e-006, 2.65172e-006, 5.25477e-006, 1.03780e-005, 2.04870e-005,
        4.05198e-005, 7.97914e-005, 1.56511e-004, 3.04954e-004, 5.99157e-004,
        1.16544e-003, 2.27488e-003, 3.91006e-003
    };

    // Pre-built SOS matrices for all 28 bands (static readonly, computed once)
    private static readonly double[][,] BandSos = BuildBandSos();

    private static double[][,] BuildBandSos()
    {
        double[][,] all = new double[28][,];
        for (int b = 0; b < 28; b++)
        {
            double[,] s = new double[3, 6];
            for (int sec = 0; sec < 3; sec++)
            {
                // b coefficients unchanged from reference
                s[sec, 0] = RefSos[sec, 0];
                s[sec, 1] = RefSos[sec, 1];
                s[sec, 2] = RefSos[sec, 2];
                // a0 always 1
                s[sec, 3] = 1.0;
                // a1 = ref_a1 - diff_a1
                s[sec, 4] = RefSos[sec, 4] - DiffA[b, sec, 0];
                // a2 = ref_a2 - diff_a2
                s[sec, 5] = RefSos[sec, 5] - DiffA[b, sec, 1];
            }
            all[b] = s;
        }
        return all;
    }

    // ------------------------------------------------------------------
    // Shared Bark axis (for public API)
    // ------------------------------------------------------------------
    private static readonly double[] BarkAxisCache = Interp.Linspace(0.1, 24.0, 240);

    // ------------------------------------------------------------------
    // Public entry point
    // ------------------------------------------------------------------

    /// <summary>
    /// Computes Zwicker time-varying loudness from a time signal.
    /// </summary>
    /// <param name="signal">Time signal [Pa]. Must be at ≥ 48 kHz (resampled if lower).</param>
    /// <param name="fs">Sampling frequency [Hz].</param>
    /// <param name="fieldType">"free" (default) or "diffuse".</param>
    /// <returns>
    /// (<c>N</c> — overall loudness vs. time [sone];
    ///  <c>NSpecific</c> — specific loudness [sone/Bark] shape (240, nTime);
    ///  <c>BarkAxis</c> — Bark axis length 240;
    ///  <c>TimeAxis</c> — time axis [s]).
    /// </returns>
    public static (double[] N, double[,] NSpecific, double[] BarkAxis, double[] TimeAxis)
        Compute(ReadOnlySpan<double> signal, int fs, string fieldType = "free")
    {
        // Resample to 48 kHz if needed
        double[] sig48;
        if (fs != 48000)
        {
            Console.WriteLine("[Warning] Signal resampled to 48 kHz for Zwicker time-varying loudness.");
            sig48 = Resample.Apply(signal, fs, 48000);
        }
        else
        {
            sig48 = signal.ToArray();
        }

        // Step 1: ISO 532-1 third-octave filter bank at 2 kHz output rate
        var (specThird, timeAxis) = ThirdOctaveLevels(sig48);

        // Step 2: Core loudness per time frame (21, nTime)
        double[,] coreLoudness = MainLoudness.ComputeBatch(specThird, fieldType);

        // Step 3: Nonlinear temporal decay
        coreLoudness = NlLoudness(coreLoudness);

        // Step 4: Specific loudness and overall loudness per time frame
        var (N, NSpec) = CalcSlopes.ComputeBatch(coreLoudness);

        // Step 5: Temporal weighting on total loudness
        double[] filtN = TemporalWeighting(N);

        // Step 6: Decimate by 4 (2 kHz → 500 Hz = 2 ms resolution)
        const int decFactor = 4;
        double[] Nout     = Decimate1D(filtN, decFactor);
        double[,] NSpecOut = Decimate2D(NSpec, decFactor);
        double[] tAxisOut  = Decimate1D(timeAxis, decFactor);

        return (Nout, NSpecOut, BarkAxisCache, tAxisOut);
    }

    // ------------------------------------------------------------------
    // Step 1: ISO 532-1 Third-octave filter bank
    // ------------------------------------------------------------------

    private static (double[,] specThird, double[] timeAxis) ThirdOctaveLevels(double[] sig)
    {
        const int nBands     = 28;
        const int fs         = 48000;
        const int decFactor  = fs / 2000; // = 24
        const double tinyVal = 1e-12;
        const double iRef    = 4e-10;

        int nTime = sig.Length / decFactor;
        double[,] levels = new double[nBands, nTime];

        double[] filtered = new double[sig.Length];
        double[] smoothed = new double[sig.Length];

        for (int b = 0; b < nBands; b++)
        {
            // Apply ISO SOS filter
            SosFilter.Process(BandSos[b], sig, filtered);

            // Scale by filter gain
            double gain = FilterGain[b];
            for (int i = 0; i < filtered.Length; i++) filtered[i] *= gain;

            // Square and smooth (3× first-order LP)
            double centerFreq = Math.Pow(10.0, (b - 16.0) / 10.0) * 1000.0;
            SquareAndSmooth(filtered, centerFreq, fs, smoothed);

            // Decimate and convert to SPL dB
            for (int t = 0; t < nTime; t++)
            {
                double val = smoothed[t * decFactor];
                levels[b, t] = 10.0 * Math.Log10((val + tinyVal) / iRef);
            }
        }

        // Time axis
        double[] tAxis = new double[nTime];
        double dt = (double)sig.Length / fs / nTime;
        for (int t = 0; t < nTime; t++) tAxis[t] = t * dt;

        return (levels, tAxis);
    }

    /// <summary>
    /// Squares the signal then applies 3 cascaded first-order LP filters in-place.
    /// Matches MoSQITo <c>_square_and_smooth</c>.
    /// </summary>
    private static void SquareAndSmooth(
        ReadOnlySpan<double> input, double centerFreq, int fs, Span<double> output)
    {
        // Frequency-dependent time constant
        double tau = centerFreq <= 1000.0
            ? 2.0 / (3.0 * centerFreq)
            : 2.0 / 3000.0;

        double a1 = Math.Exp(-1.0 / (fs * tau));
        double b0 = 1.0 - a1;

        // Square
        for (int i = 0; i < input.Length; i++) output[i] = input[i] * input[i];

        // 3 cascaded first-order LP: y[n] = b0*x[n] + a1*y[n-1]
        for (int pass = 0; pass < 3; pass++)
        {
            double y = 0.0;
            for (int i = 0; i < output.Length; i++)
            {
                y = b0 * output[i] + a1 * y;
                output[i] = y;
            }
        }
    }

    // ------------------------------------------------------------------
    // Step 3: Nonlinear temporal decay (_nl_loudness)
    // ------------------------------------------------------------------

    /// <summary>
    /// Simulates the nonlinear temporal decay of the hearing system.
    /// Ported from MoSQITo <c>_nonlinear_decay.py</c>.
    /// </summary>
    private static double[,] NlLoudness(double[,] coreLoudness)
    {
        int nBark = coreLoudness.GetLength(0); // 21
        int nTime = coreLoudness.GetLength(1);

        const int    sampleRate = 2000;
        const int    nlIter     = 24;
        const double tShort     = 0.005;
        const double tLong      = 0.015;
        const double tVar       = 0.075;

        double deltaT  = 1.0 / (sampleRate * nlIter);
        double P       = (tVar + tLong) / (tVar * tShort);
        double Q       = 1.0 / (tShort * tVar);
        double disc    = P * P / 4.0 - Q;
        double sqrtD   = Math.Sqrt(disc);
        double lambda1 = -P / 2.0 + sqrtD;
        double lambda2 = -P / 2.0 - sqrtD;
        double den     = tVar * (lambda1 - lambda2);
        double e1      = Math.Exp(lambda1 * deltaT);
        double e2      = Math.Exp(lambda2 * deltaT);

        double B0 = (e1 - e2) / den;
        double B1 = ((tVar * lambda2 + 1.0) * e1 - (tVar * lambda1 + 1.0) * e2) / den;
        double B2 = ((tVar * lambda1 + 1.0) * e1 - (tVar * lambda2 + 1.0) * e2) / den;
        double B3 = (tVar * lambda1 + 1.0) * (tVar * lambda2 + 1.0) * (e1 - e2) / den;
        double B4 = Math.Exp(-deltaT / tLong);
        double B5 = Math.Exp(-deltaT / tVar);

        // Output — we process per-Bark row independently
        double[,] result = new double[nBark, nTime];

        int nExpanded = nTime * nlIter;

        // Temporary buffers (per-Bark row)
        double[] ui = new double[nExpanded];
        double[] uo = new double[nExpanded];
        double[] u2 = new double[nExpanded];

        for (int b = 0; b < nBark; b++)
        {
            // Build expanded signal with linear interpolation between frames
            for (int t = 0; t < nTime; t++)
            {
                double cur  = coreLoudness[b, t];
                double next = (t < nTime - 1) ? coreLoudness[b, t + 1] : cur;
                double delta = (next - cur) / nlIter;
                for (int k = 0; k < nlIter; k++)
                    ui[t * nlIter + k] = cur + k * delta;
            }

            // Initialise uo and u2
            Array.Copy(ui, uo, nExpanded);
            Array.Clear(u2, 0, nExpanded);
            if (coreLoudness[b, 0] >= 1e-5)
                u2[0] = coreLoudness[b, 0] * (1.0 - B5);

            // Nonlinear IIR loop
            for (int col = 1; col < nExpanded; col++)
            {
                double uiCur  = ui[col];
                double uoLast = uo[col - 1];
                double u2Last = u2[col - 1];

                double uo2a = uoLast * B2 - u2Last * B3;
                if (uoLast > u2Last && uo2a >= uiCur)
                {
                    uo[col] = uo2a;
                }
                else
                {
                    double uo2b = uoLast * B4;
                    if (uoLast <= u2Last && uo2b >= uiCur)
                        uo[col] = uo2b;
                    // else uo[col] = uiCur (already set by Array.Copy)
                }

                u2[col] = uo[col]; // default: u2 = uo

                double u22 = uoLast * B0 - u2Last * B1;
                if (uiCur < uoLast && uoLast > u2Last && u22 <= uo[col])
                    u2[col] = u22;

                double u2_2 = (u2Last - uiCur) * B5 + uiCur;
                bool notSpecialCase = !(Math.Abs(uiCur - uoLast) < 1e-5 && uo[col] <= u2Last);
                if (uiCur >= uoLast && notSpecialCase)
                    u2[col] = u2_2;
            }

            // Extract first inner sample of each frame
            for (int t = 0; t < nTime; t++)
                result[b, t] = uo[t * nlIter];
        }

        return result;
    }

    // ------------------------------------------------------------------
    // Step 5: Temporal weighting (_temporal_weighting + _lowpass_intp)
    // ------------------------------------------------------------------

    /// <summary>
    /// Temporal weighting: two LP filters (3.5 ms and 70 ms), weighted 0.47/0.53.
    /// Ported from MoSQITo <c>_temporal_weighting.py</c>.
    /// </summary>
    private static double[] TemporalWeighting(double[] loudness)
    {
        const int sampleRate = 2000;
        double[] lp1 = LowpassIntp(loudness, 3.5e-3, sampleRate);
        double[] lp2 = LowpassIntp(loudness, 70.0e-3, sampleRate);

        double[] result = new double[loudness.Length];
        for (int i = 0; i < result.Length; i++)
            result[i] = 0.47 * lp1[i] + 0.53 * lp2[i];
        return result;
    }

    /// <summary>
    /// 1st-order low-pass with linear interpolation for increased precision.
    /// Ported from MoSQITo <c>_lowpass_intp.py</c>.
    /// </summary>
    private static double[] LowpassIntp(double[] loudness, double tau, int sampleRate)
    {
        int n       = loudness.Length;
        const int lpIter = 24;
        double a1  = Math.Exp(-1.0 / (sampleRate * lpIter * tau));
        double b0  = 1.0 - a1;

        // Expand by lpIter with linear interpolation between samples
        double[] expanded = new double[n * lpIter];
        for (int t = 0; t < n; t++)
        {
            double cur  = loudness[t];
            double next = (t < n - 1) ? loudness[t + 1] : cur;
            double delta = (next - cur) / lpIter;
            for (int k = 0; k < lpIter; k++)
                expanded[t * lpIter + k] = cur + k * delta;
        }

        // Apply 1st-order LP: y[n] = b0*x[n] + a1*y[n-1]
        double y = 0.0;
        for (int i = 0; i < expanded.Length; i++)
        {
            y = b0 * expanded[i] + a1 * y;
            expanded[i] = y;
        }

        // Return first sample of each expanded block
        double[] result = new double[n];
        for (int t = 0; t < n; t++) result[t] = expanded[t * lpIter];
        return result;
    }

    // ------------------------------------------------------------------
    // Decimation helpers
    // ------------------------------------------------------------------

    private static double[] Decimate1D(double[] x, int factor)
    {
        int n = x.Length / factor;
        double[] result = new double[n];
        for (int i = 0; i < n; i++) result[i] = x[i * factor];
        return result;
    }

    private static double[,] Decimate2D(double[,] x, int factor)
    {
        int rows = x.GetLength(0);
        int cols = x.GetLength(1) / factor;
        double[,] result = new double[rows, cols];
        for (int r = 0; r < rows; r++)
            for (int c = 0; c < cols; c++)
                result[r, c] = x[r, c * factor];
        return result;
    }
}
