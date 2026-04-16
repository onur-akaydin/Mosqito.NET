using MathNet.Numerics.Interpolation;

namespace Mosqito.Dsp;

/// <summary>
/// 1-D interpolation utilities matching scipy.interpolate.interp1d and
/// pchip_interpolate behaviour used throughout MoSQITo.
/// </summary>
public static class Interp
{
    // ------------------------------------------------------------------
    // Linear interpolation  (numpy.interp / interp1d(kind='linear'))
    // ------------------------------------------------------------------

    /// <summary>
    /// Piecewise-linear interpolation with boundary clamping.
    /// Matches numpy.interp(x, xp, fp) exactly (clamps outside range).
    /// </summary>
    /// <param name="x">Query point.</param>
    /// <param name="xp">Monotonically increasing x-coordinates of data points.</param>
    /// <param name="yp">y-coordinates of data points (same length as xp).</param>
    public static double Linear(double x, ReadOnlySpan<double> xp, ReadOnlySpan<double> yp)
    {
        if (x <= xp[0]) return yp[0];
        if (x >= xp[xp.Length - 1]) return yp[yp.Length - 1];
        int lo = BinarySearch(xp, x);
        double t = (x - xp[lo]) / (xp[lo + 1] - xp[lo]);
        return yp[lo] + t * (yp[lo + 1] - yp[lo]);
    }

    /// <summary>
    /// Piecewise-linear interpolation for a batch of query points.
    /// Output written into <paramref name="output"/> (must have same length as <paramref name="x"/>).
    /// </summary>
    public static void Linear(ReadOnlySpan<double> x, ReadOnlySpan<double> xp,
        ReadOnlySpan<double> yp, Span<double> output)
    {
        if (output.Length < x.Length) throw new ArgumentException("Output too short.", nameof(output));
        for (int i = 0; i < x.Length; i++)
            output[i] = Linear(x[i], xp, yp);
    }

    /// <summary>Returns a new array with linearly interpolated values.</summary>
    public static double[] Linear(ReadOnlySpan<double> x, ReadOnlySpan<double> xp,
        ReadOnlySpan<double> yp)
    {
        double[] result = new double[x.Length];
        Linear(x, xp, yp, result);
        return result;
    }

    // ------------------------------------------------------------------
    // Cubic spline  (interp1d(kind='cubic'))
    // ------------------------------------------------------------------

    /// <summary>
    /// Cubic spline interpolation at a single query point.
    /// Backed by MathNet.Numerics NaturalSpline (not-a-knot end conditions).
    /// </summary>
    public static double Cubic(double x, double[] xp, double[] yp)
    {
        var spline = CubicSpline.InterpolateNatural(xp, yp);
        return spline.Interpolate(x);
    }

    /// <summary>
    /// Cubic spline interpolation for a batch of query points.
    /// </summary>
    public static double[] Cubic(ReadOnlySpan<double> x, double[] xp, double[] yp)
    {
        var spline = CubicSpline.InterpolateNatural(xp, yp);
        double[] result = new double[x.Length];
        for (int i = 0; i < x.Length; i++) result[i] = spline.Interpolate(x[i]);
        return result;
    }

    // ------------------------------------------------------------------
    // PCHIP  (scipy.interpolate.pchip_interpolate)
    // Implements the Fritsch-Carlson monotone cubic Hermite interpolant,
    // matching scipy.interpolate.PchipInterpolator exactly including
    // endpoint derivative conditions and extrapolation behaviour.
    // ------------------------------------------------------------------

    /// <summary>
    /// Monotone cubic (PCHIP) interpolation at a single query point.
    /// Matches scipy.interpolate.pchip_interpolate(xi, yi, x).
    /// </summary>
    public static double Pchip(double x, double[] xi, double[] yi)
    {
        double[] d = PchipDerivatives(xi, yi);
        return PchipEvalOne(x, xi, yi, d);
    }

    /// <summary>
    /// Monotone cubic (PCHIP) interpolation for a batch of query points.
    /// </summary>
    public static double[] Pchip(ReadOnlySpan<double> x, double[] xi, double[] yi)
    {
        double[] d = PchipDerivatives(xi, yi);
        double[] result = new double[x.Length];
        for (int i = 0; i < x.Length; i++)
            result[i] = PchipEvalOne(x[i], xi, yi, d);
        return result;
    }

    /// <summary>
    /// Computes the PCHIP derivative (slope) at each data point using the
    /// Fritsch-Carlson formula. Matches scipy's PchipInterpolator._find_derivatives.
    /// </summary>
    private static double[] PchipDerivatives(double[] xi, double[] yi)
    {
        int n = xi.Length;
        double[] h = new double[n - 1];
        double[] delta = new double[n - 1];
        for (int k = 0; k < n - 1; k++)
        {
            h[k] = xi[k + 1] - xi[k];
            delta[k] = (yi[k + 1] - yi[k]) / h[k];
        }

        double[] d = new double[n];

        // Interior points: harmonic mean formula (Fritsch-Carlson)
        for (int k = 1; k < n - 1; k++)
        {
            if (delta[k - 1] * delta[k] <= 0.0)
            {
                d[k] = 0.0;
            }
            else
            {
                double w1 = 2.0 * h[k] + h[k - 1];
                double w2 = h[k] + 2.0 * h[k - 1];
                d[k] = (w1 + w2) / (w1 / delta[k - 1] + w2 / delta[k]);
            }
        }

        // Endpoint derivatives (non-centred, one-sided, matches scipy _edge_case)
        d[0]     = PchipEdge(h[0], h[1], delta[0], delta[1]);
        d[n - 1] = PchipEdge(h[n - 2], h[n - 3], delta[n - 2], delta[n - 3]);

        return d;
    }

    private static double PchipEdge(double h0, double h1, double m0, double m1)
    {
        double d = ((2.0 * h0 + h1) * m0 - h0 * m1) / (h0 + h1);
        if (Math.Sign(d) != Math.Sign(m0))
            return 0.0;
        if (Math.Sign(m0) != Math.Sign(m1) && Math.Abs(d) > 3.0 * Math.Abs(m0))
            return 3.0 * m0;
        return d;
    }

    private static double PchipEvalOne(double x, double[] xi, double[] yi, double[] d)
    {
        int n = xi.Length;

        // Extrapolation: use the polynomial at the nearest endpoint interval
        int lo;
        if (x <= xi[0])
            lo = 0;
        else if (x >= xi[n - 1])
            lo = n - 2;
        else
        {
            lo = 0;
            int hi = n - 2;
            while (lo < hi)
            {
                int mid = (lo + hi) >> 1;
                if (xi[mid + 1] <= x) lo = mid + 1;
                else hi = mid;
            }
        }

        double h = xi[lo + 1] - xi[lo];
        double t = (x - xi[lo]) / h;

        // Cubic Hermite basis functions
        double t2 = t * t;
        double t3 = t2 * t;
        double h00 = 2 * t3 - 3 * t2 + 1;
        double h10 = t3 - 2 * t2 + t;
        double h01 = -2 * t3 + 3 * t2;
        double h11 = t3 - t2;

        return h00 * yi[lo] + h10 * h * d[lo] + h01 * yi[lo + 1] + h11 * h * d[lo + 1];
    }

    // ------------------------------------------------------------------
    // Linspace helper (numpy.linspace)
    // ------------------------------------------------------------------

    /// <summary>
    /// Returns <paramref name="num"/> evenly-spaced values in [<paramref name="start"/>,
    /// <paramref name="stop"/>]. Matches numpy.linspace(start, stop, num).
    /// </summary>
    public static double[] Linspace(double start, double stop, int num)
    {
        if (num <= 0) return Array.Empty<double>();
        if (num == 1) return new[] { start };
        double[] result = new double[num];
        double step = (stop - start) / (num - 1);
        for (int i = 0; i < num; i++) result[i] = start + i * step;
        result[num - 1] = stop; // ensure exact endpoint
        return result;
    }

    /// <summary>Returns N evenly-spaced integers in [start, stop).</summary>
    public static double[] Arange(double start, double stop, double step = 1.0)
    {
        int n = (int)Math.Ceiling((stop - start) / step);
        if (n < 0) n = 0;
        double[] result = new double[n];
        for (int i = 0; i < n; i++) result[i] = start + i * step;
        return result;
    }

    // ------------------------------------------------------------------
    // Binary search helper
    // ------------------------------------------------------------------
    private static int BinarySearch(ReadOnlySpan<double> xp, double x)
    {
        int lo = 0, hi = xp.Length - 2;
        while (lo < hi)
        {
            int mid = (lo + hi) >> 1;
            if (xp[mid + 1] <= x) lo = mid + 1;
            else hi = mid;
        }
        return lo;
    }
}
