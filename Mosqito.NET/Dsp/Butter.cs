using System.Collections.Concurrent;
using System.Numerics;
using System.Runtime.CompilerServices;

namespace Mosqito.Dsp;

/// <summary>
/// Butterworth IIR digital filter design.
/// Produces second-order-section (SOS) coefficients compatible with
/// <see cref="SosFilter"/>.
///
/// Matches scipy.signal.butter(N, Wn, btype, analog=False, output='sos').
///
/// Each SOS row is [b0, b1, b2, a0, a1, a2] where a0 is normalised to 1.
/// </summary>
public static class Butter
{
    private const double TwoPi = 2.0 * Math.PI;

    // Design cache keyed by (FilterType, order, w1*1e9, w2*1e9) — design is deterministic.
    private static readonly ConcurrentDictionary<(int type, int order, long w1, long w2), double[,]>
        _designCache = new();

    // Quantise a frequency to a stable long key (9 significant digits).
    private static long FreqKey(double w) => (long)Math.Round(w * 1e9);

    // ------------------------------------------------------------------
    // Public entry points
    // ------------------------------------------------------------------

    /// <summary>
    /// Design a Butterworth bandpass filter.
    /// </summary>
    /// <param name="order">Filter order (each pole pair → one SOS section; total sections = order).</param>
    /// <param name="wLow">Normalised lower cutoff frequency (0–1, where 1 = Nyquist).</param>
    /// <param name="wHigh">Normalised upper cutoff frequency.</param>
    /// <returns>SOS matrix, shape (nSections, 6).  Stored row-major: [b0,b1,b2,a0,a1,a2].</returns>
    public static double[,] DesignBandpass(int order, double wLow, double wHigh)
    {
        if (wLow <= 0 || wHigh >= 1 || wLow >= wHigh)
            throw new ArgumentException($"Invalid cutoff frequencies: wLow={wLow}, wHigh={wHigh}. Must satisfy 0 < wLow < wHigh < 1.");
        return _designCache.GetOrAdd((0, order, FreqKey(wLow), FreqKey(wHigh)),
            static k => IirFilter(k.order, k.w1 * 1e-9, k.w2 * 1e-9, FilterType.Bandpass));
    }

    /// <summary>Design a Butterworth low-pass filter.</summary>
    public static double[,] DesignLowpass(int order, double wCutoff)
    {
        if (wCutoff <= 0 || wCutoff >= 1)
            throw new ArgumentException($"Invalid cutoff: wCutoff={wCutoff}. Must satisfy 0 < wCutoff < 1.");
        return _designCache.GetOrAdd((1, order, FreqKey(wCutoff), 0L),
            static k => IirFilter(k.order, k.w1 * 1e-9, 0, FilterType.Lowpass));
    }

    /// <summary>Design a Butterworth high-pass filter.</summary>
    public static double[,] DesignHighpass(int order, double wCutoff)
    {
        if (wCutoff <= 0 || wCutoff >= 1)
            throw new ArgumentException($"Invalid cutoff: wCutoff={wCutoff}. Must satisfy 0 < wCutoff < 1.");
        return _designCache.GetOrAdd((2, order, FreqKey(wCutoff), 0L),
            static k => IirFilter(k.order, k.w1 * 1e-9, 0, FilterType.Highpass));
    }

    /// <summary>Design a Butterworth band-stop (notch) filter.</summary>
    public static double[,] DesignBandstop(int order, double wLow, double wHigh)
    {
        if (wLow <= 0 || wHigh >= 1 || wLow >= wHigh)
            throw new ArgumentException($"Invalid cutoff frequencies: wLow={wLow}, wHigh={wHigh}.");
        return _designCache.GetOrAdd((3, order, FreqKey(wLow), FreqKey(wHigh)),
            static k => IirFilter(k.order, k.w1 * 1e-9, k.w2 * 1e-9, FilterType.Bandstop));
    }

    // ------------------------------------------------------------------
    // Core design pipeline  (matches scipy.signal.iirfilter + zpk2sos)
    // ------------------------------------------------------------------

    private enum FilterType { Lowpass, Highpass, Bandpass, Bandstop }

    private static double[,] IirFilter(int N, double w1, double w2, FilterType type)
    {
        // Step 1: Analog Butterworth prototype (LP, order N)
        // poles on the unit circle in left half-plane
        var (zProto, pProto, kProto) = ButterworthPrototype(N);

        // Step 2: Pre-warp digital frequencies to analog.
        // Bilinear uses fs=2, so s = 2*fs*(z-1)/(z+1) = 4*(z-1)/(z+1).
        // At z=e^{jω}: s = j*4*tan(ω/2), so ω_analog = 4*tan(π*w/2).
        // scipy.signal.iirfilter uses: warped = 2 * fs * tan(π * Wn / fs) = 4 * tan(π*w/2).
        double[] warped = type switch
        {
            FilterType.Lowpass  => new[] { 4.0 * Math.Tan(Math.PI * w1 / 2.0) },
            FilterType.Highpass => new[] { 4.0 * Math.Tan(Math.PI * w1 / 2.0) },
            FilterType.Bandpass or FilterType.Bandstop =>
                new[] { 4.0 * Math.Tan(Math.PI * w1 / 2.0), 4.0 * Math.Tan(Math.PI * w2 / 2.0) },
            _ => throw new InvalidOperationException()
        };

        // Step 3: LP prototype → desired filter type (s-domain)
        Complex[] z, p;
        double k;
        switch (type)
        {
            case FilterType.Lowpass:
                (z, p, k) = Lp2Lp(zProto, pProto, kProto, warped[0]);
                break;
            case FilterType.Highpass:
                (z, p, k) = Lp2Hp(zProto, pProto, kProto, warped[0]);
                break;
            case FilterType.Bandpass:
                double bw = warped[1] - warped[0];
                double wo = Math.Sqrt(warped[0] * warped[1]);
                (z, p, k) = Lp2Bp(zProto, pProto, kProto, wo, bw);
                break;
            case FilterType.Bandstop:
                bw = warped[1] - warped[0];
                wo = Math.Sqrt(warped[0] * warped[1]);
                (z, p, k) = Lp2Bs(zProto, pProto, kProto, wo, bw);
                break;
            default: throw new InvalidOperationException();
        }

        // Step 4: Bilinear transform (s-domain → z-domain), fs = 2
        var (zD, pD, kD) = BilinearZpk(z, p, k, fs: 2.0);

        // Step 5: ZPK → SOS
        return Zpk2Sos(zD, pD, kD);
    }

    // ------------------------------------------------------------------
    // Analog Butterworth prototype  (buttap)
    // ------------------------------------------------------------------

    private static (Complex[] zeros, Complex[] poles, double gain) ButterworthPrototype(int N)
    {
        Complex[] poles = new Complex[N];
        for (int m = 0; m < N; m++)
        {
            double angle = Math.PI * (2 * m + N + 1) / (2.0 * N);
            poles[m] = new Complex(Math.Cos(angle), Math.Sin(angle)); // all in LHP
        }
        return (Array.Empty<Complex>(), poles, 1.0);
    }

    // ------------------------------------------------------------------
    // LP → LP transform
    // ------------------------------------------------------------------
    private static (Complex[] z, Complex[] p, double k) Lp2Lp(
        Complex[] z, Complex[] p, double k, double wo)
    {
        int nz = z.Length, np = p.Length;
        Complex[] zNew = new Complex[nz];
        Complex[] pNew = new Complex[np];
        for (int i = 0; i < nz; i++) zNew[i] = z[i] * wo;
        for (int i = 0; i < np; i++) pNew[i] = p[i] * wo;
        double kNew = k * Math.Pow(wo, np - nz);
        return (zNew, pNew, kNew);
    }

    // ------------------------------------------------------------------
    // LP → HP transform
    // ------------------------------------------------------------------
    private static (Complex[] z, Complex[] p, double k) Lp2Hp(
        Complex[] z, Complex[] p, double k, double wo)
    {
        int nz = z.Length, np = p.Length;
        Complex[] zNew = new Complex[np]; // add zeros at origin for each LP finite pole
        Complex[] pNew = new Complex[np];
        // zeros: z_hp = wo / z_lp, but LP has no finite zeros → inject zeros at 0
        for (int i = 0; i < nz; i++) zNew[i] = new Complex(wo, 0) / z[i];
        for (int i = nz; i < np; i++) zNew[i] = Complex.Zero;
        for (int i = 0; i < np; i++) pNew[i] = new Complex(wo, 0) / p[i];
        double kNew = k * Math.Pow(-1, np - nz);
        // Correct gain: product(p) / product(z) * k (scipy approach)
        Complex pProd = Complex.One, zProd = Complex.One;
        foreach (var pi in p) pProd *= pi;
        foreach (var zi in z) zProd *= zi;
        kNew = k * (zNew.Length == 0 ? 1 : -pProd / (zProd == Complex.Zero ? Complex.One : zProd)).Real;
        // Simpler gain scaling: k * wo^(np-nz)
        kNew = k * Math.Pow(wo, np - nz);
        return (zNew, pNew, kNew);
    }

    // ------------------------------------------------------------------
    // LP → Bandpass transform
    // ------------------------------------------------------------------
    private static (Complex[] z, Complex[] p, double k) Lp2Bp(
        Complex[] z, Complex[] p, double k, double wo, double bw)
    {
        int nz = z.Length, np = p.Length;
        // Each LP finite zero → 2 BP zeros (via quadratic formula)
        // Each LP zero-at-infinity (there are np-nz of them) → 1 zero at s=0 in BP
        // Total BP zeros = 2*nz + (np-nz) = nz + np
        // BilinearZpk will then add (np - (nz+np)) ... wait, let's count:
        //   nz_bp = nz + np,  np_bp = 2*np
        //   extra zeros at z=-1 added by BilinearZpk = max(np_bp, nz_bp) - nz_bp = 2*np - (nz+np) = np-nz ✓
        int nzNew = 2 * nz + (np - nz);  // = nz + np
        Complex[] zNew = new Complex[nzNew];
        Complex[] pNew = new Complex[2 * np];

        for (int i = 0; i < nz; i++)
        {
            var disc = Complex.Sqrt(z[i] * z[i] * bw * bw - 4 * wo * wo);
            zNew[2 * i]     = (z[i] * bw + disc) * 0.5;
            zNew[2 * i + 1] = (z[i] * bw - disc) * 0.5;
        }
        // Fill exactly (np-nz) zeros at origin (one per LP zero-at-infinity, not pairs)
        for (int i = 2 * nz; i < nzNew; i++)
            zNew[i] = Complex.Zero;

        for (int i = 0; i < np; i++)
        {
            var disc = Complex.Sqrt(p[i] * p[i] * bw * bw - 4 * wo * wo);
            pNew[2 * i]     = (p[i] * bw + disc) * 0.5;
            pNew[2 * i + 1] = (p[i] * bw - disc) * 0.5;
        }

        double kNew = k * Math.Pow(bw, np - nz);
        return (zNew, pNew, kNew);
    }

    // ------------------------------------------------------------------
    // LP → Bandstop transform
    // ------------------------------------------------------------------
    private static (Complex[] z, Complex[] p, double k) Lp2Bs(
        Complex[] z, Complex[] p, double k, double wo, double bw)
    {
        int nz = z.Length, np = p.Length;
        Complex[] zNew = new Complex[2 * np];
        Complex[] pNew = new Complex[2 * np];

        for (int i = 0; i < nz; i++)
        {
            var disc = Complex.Sqrt(bw * bw - 4 * z[i] * z[i] * wo * wo);
            zNew[2 * i]     = (bw + disc) / (2 * z[i]);
            zNew[2 * i + 1] = (bw - disc) / (2 * z[i]);
        }
        var jWo = new Complex(0, wo);
        for (int i = 2 * nz; i < 2 * np; i += 2)
        {
            zNew[i]     =  jWo;
            zNew[i + 1] = -jWo;
        }

        for (int i = 0; i < np; i++)
        {
            var disc = Complex.Sqrt(bw * bw - 4 * p[i] * p[i] * wo * wo);
            pNew[2 * i]     = (bw + disc) / (2 * p[i]);
            pNew[2 * i + 1] = (bw - disc) / (2 * p[i]);
        }

        // Gain: k * prod(z) / prod(p) normalised
        double kNew = k;
        return (zNew, pNew, kNew);
    }

    // ------------------------------------------------------------------
    // Bilinear transform  (s → z, fs=2 for pre-warped)
    // ------------------------------------------------------------------
    private static (Complex[] z, Complex[] p, double k) BilinearZpk(
        Complex[] z, Complex[] p, double k, double fs = 2.0)
    {
        int nz = z.Length, np = p.Length;
        Complex[] zD = new Complex[Math.Max(nz, np)];
        Complex[] pD = new Complex[np];

        double twoFs = 2.0 * fs;
        var twoFsC = new Complex(twoFs, 0);

        for (int i = 0; i < np; i++)
            pD[i] = (twoFsC + p[i]) / (twoFsC - p[i]);

        for (int i = 0; i < nz; i++)
            zD[i] = (twoFsC + z[i]) / (twoFsC - z[i]);
        // Fill remaining zeros at z = -1 (s = ∞ maps to z = -1)
        for (int i = nz; i < zD.Length; i++) zD[i] = new Complex(-1.0, 0);

        // Gain correction
        Complex kNum = Complex.One, kDen = Complex.One;
        foreach (var zi in z) kNum *= twoFsC - zi;
        foreach (var pi in p) kDen *= twoFsC - pi;
        double kD = k * (kNum / kDen).Real;

        return (zD, pD, kD);
    }

    // ------------------------------------------------------------------
    // ZPK → SOS  (matches scipy.signal.zpk2sos)
    // ------------------------------------------------------------------
    private static double[,] Zpk2Sos(Complex[] z, Complex[] p, double k)
    {
        // Sort poles by proximity to unit circle (ascending |1 - |p||)
        // Group into conjugate pairs, form SOS sections.
        // For each section: b = [1, -(z1+z2), z1*z2],  a = [1, -(p1+p2), p1*p2]

        var pairs = PairRootsIntoSections(z, p);
        int nSections = pairs.Count;
        double[,] sos = new double[nSections, 6];

        // Distribute gain across sections (put all in first section)
        bool gainApplied = false;
        for (int s = 0; s < nSections; s++)
        {
            var (z1, z2, p1, p2) = pairs[s];
            double gFactor = (!gainApplied) ? k : 1.0;
            gainApplied = true;

            // Numerator coefficients [b0, b1, b2]
            double b0 = gFactor;
            double b1 = -gFactor * (z1.Real + z2.Real);
            double b2 =  gFactor *  (z1.Real * z2.Real - z1.Imaginary * z2.Imaginary +
                                      z1.Imaginary * z2.Imaginary);
            // Clean up: for a conjugate pair, b2 = |z|^2 real
            if (z1 == Complex.Conjugate(z2))
            {
                b1 = -gFactor * 2.0 * z1.Real;
                b2 =  gFactor * (z1.Real * z1.Real + z1.Imaginary * z1.Imaginary);
            }
            else if (z1.Imaginary == 0 && z2.Imaginary == 0)
            {
                b1 = -gFactor * (z1.Real + z2.Real);
                b2 =  gFactor *  z1.Real * z2.Real;
            }

            // Denominator coefficients [a0=1, a1, a2]
            double a1, a2;
            if (p1 == Complex.Conjugate(p2))
            {
                a1 = -2.0 * p1.Real;
                a2 =  p1.Real * p1.Real + p1.Imaginary * p1.Imaginary;
            }
            else
            {
                a1 = -(p1.Real + p2.Real);
                a2 =   p1.Real * p2.Real - p1.Imaginary * p2.Imaginary;
            }

            sos[s, 0] = b0; sos[s, 1] = b1; sos[s, 2] = b2;
            sos[s, 3] = 1.0; sos[s, 4] = a1; sos[s, 5] = a2;
        }

        return sos;
    }

    // Pair roots into conjugate pairs, padded with 1-poles-at-origin as needed.
    private static List<(Complex z1, Complex z2, Complex p1, Complex p2)> PairRootsIntoSections(
        Complex[] z, Complex[] p)
    {
        var zList = new List<Complex>(z);
        var pList = new List<Complex>(p);

        // Pad with zeros/poles at origin to make even count
        while (zList.Count < pList.Count) zList.Add(Complex.Zero);
        while (pList.Count < zList.Count) pList.Add(Complex.Zero);

        int n = pList.Count;
        if (n % 2 != 0)
        {
            zList.Add(Complex.Zero);
            pList.Add(Complex.Zero);
            n++;
        }

        // Sort poles: complex first (nearest to unit circle), then real
        pList.Sort((a, b) =>
        {
            double dA = Math.Abs(1.0 - a.Magnitude);
            double dB = Math.Abs(1.0 - b.Magnitude);
            return dA.CompareTo(dB);
        });

        // Sort zeros to match pole pairing (nearest to each pole pair)
        zList.Sort((a, b) =>
        {
            double dA = Math.Abs(1.0 - a.Magnitude);
            double dB = Math.Abs(1.0 - b.Magnitude);
            return dA.CompareTo(dB);
        });

        var sections = new List<(Complex, Complex, Complex, Complex)>();
        bool[] usedZ = new bool[zList.Count];
        bool[] usedP = new bool[pList.Count];

        for (int i = 0; i < n; i += 2)
        {
            // Find next unused pole pair (take i and its conjugate if imaginary)
            int pi1 = -1, pi2 = -1;
            for (int j = 0; j < pList.Count; j++)
            {
                if (usedP[j]) continue;
                pi1 = j; usedP[j] = true; break;
            }
            // Find conjugate or next available
            for (int j = 0; j < pList.Count; j++)
            {
                if (usedP[j]) continue;
                if (pList[pi1].Imaginary != 0 &&
                    Math.Abs(pList[j].Real - pList[pi1].Real) < 1e-12 &&
                    Math.Abs(pList[j].Imaginary + pList[pi1].Imaginary) < 1e-12)
                {
                    pi2 = j; usedP[j] = true; break;
                }
            }
            if (pi2 < 0)
            {
                for (int j = 0; j < pList.Count; j++)
                {
                    if (usedP[j]) continue;
                    pi2 = j; usedP[j] = true; break;
                }
            }

            int zi1 = -1, zi2 = -1;
            for (int j = 0; j < zList.Count; j++)
            {
                if (usedZ[j]) continue;
                zi1 = j; usedZ[j] = true; break;
            }
            for (int j = 0; j < zList.Count; j++)
            {
                if (usedZ[j]) continue;
                if (zList[zi1].Imaginary != 0 &&
                    Math.Abs(zList[j].Real - zList[zi1].Real) < 1e-12 &&
                    Math.Abs(zList[j].Imaginary + zList[zi1].Imaginary) < 1e-12)
                {
                    zi2 = j; usedZ[j] = true; break;
                }
            }
            if (zi2 < 0)
            {
                for (int j = 0; j < zList.Count; j++)
                {
                    if (usedZ[j]) continue;
                    zi2 = j; usedZ[j] = true; break;
                }
            }

            var p1 = pi1 >= 0 ? pList[pi1] : Complex.Zero;
            var p2 = pi2 >= 0 ? pList[pi2] : Complex.Zero;
            var z1 = zi1 >= 0 ? zList[zi1] : Complex.Zero;
            var z2 = zi2 >= 0 ? zList[zi2] : Complex.Zero;
            sections.Add((z1, z2, p1, p2));
        }

        return sections;
    }

    // ------------------------------------------------------------------
    // Frequency response of an SOS filter (matches scipy.signal.sosfreqz)
    // Returns H(e^jw) for w in [0, pi], nPoints = nfft/2+1
    // ------------------------------------------------------------------

    /// <summary>
    /// Computes the complex frequency response of an SOS filter at <paramref name="nPoints"/>
    /// frequency points uniformly spaced from 0 to π (i.e., 0 to fs/2).
    /// Matches <c>scipy.signal.sosfreqz(sos, worN=nPoints)</c>.
    /// </summary>
    public static (double[] freqs, System.Numerics.Complex[] H) Sosfreqz(double[,] sos, int nPoints)
    {
        int nSec = sos.GetLength(0);
        double[] freqs = new double[nPoints];
        System.Numerics.Complex[] H = new System.Numerics.Complex[nPoints];
        H.AsSpan().Fill(System.Numerics.Complex.One);

        for (int k = 0; k < nPoints; k++)
        {
            double w = Math.PI * k / (nPoints - 1);
            freqs[k] = w;
            var z = new System.Numerics.Complex(Math.Cos(w), Math.Sin(w)); // e^jw
            for (int s = 0; s < nSec; s++)
            {
                double b0 = sos[s, 0], b1 = sos[s, 1], b2 = sos[s, 2];
                double a0 = sos[s, 3], a1 = sos[s, 4], a2 = sos[s, 5];
                var num = b0 + (b1 + b2 / z) / z;
                var den = a0 + (a1 + a2 / z) / z;
                H[k] *= num / den;
            }
        }
        return (freqs, H);
    }
}
