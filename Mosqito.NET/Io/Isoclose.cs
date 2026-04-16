using System.Runtime.CompilerServices;

namespace Mosqito.Io;

/// <summary>
/// ISO 532-1 section 5.1–style tolerance comparison for validation.
/// Ported from MoSQITo <c>isoclose</c>.
///
/// The test compares <c>actual</c> against the range:
/// <code>
///   lower = min(desired * (1 - |rtol|), desired - |atol|)
///   upper = max(desired * (1 + |rtol|), desired + |atol|)
/// </code>
/// and passes when <c>lower &lt;= actual &lt;= upper</c> for every element.
/// </summary>
public static class Isoclose
{
    /// <summary>
    /// Checks whether <paramref name="actual"/> is within tolerance of <paramref name="desired"/>.
    /// </summary>
    /// <param name="actual">Computed value.</param>
    /// <param name="desired">Reference value.</param>
    /// <param name="rtol">Relative tolerance (default 1e-7).</param>
    /// <param name="atol">Absolute tolerance (default 0).</param>
    /// <returns><see langword="true"/> if within tolerance.</returns>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static bool Check(double actual, double desired, double rtol = 1e-7, double atol = 0.0)
    {
        double lower = Math.Min(desired * (1.0 - Math.Abs(rtol)), desired - Math.Abs(atol));
        double upper = Math.Max(desired * (1.0 + Math.Abs(rtol)), desired + Math.Abs(atol));
        return actual >= lower && actual <= upper;
    }

    /// <summary>
    /// Checks whether every element of <paramref name="actual"/> is within tolerance of the
    /// corresponding element of <paramref name="desired"/>.
    /// </summary>
    public static bool Check(ReadOnlySpan<double> actual, ReadOnlySpan<double> desired,
        double rtol = 1e-7, double atol = 0.0)
    {
        if (actual.Length != desired.Length)
            throw new ArgumentException("Arrays must have the same length.");

        double absRtol = Math.Abs(rtol);
        double absAtol = Math.Abs(atol);

        for (int i = 0; i < actual.Length; i++)
        {
            double d = desired[i];
            double lower = Math.Min(d * (1.0 - absRtol), d - absAtol);
            double upper = Math.Max(d * (1.0 + absRtol), d + absAtol);
            if (actual[i] < lower || actual[i] > upper) return false;
        }
        return true;
    }

    /// <summary>
    /// Returns (isClose, maxAbsError, maxRelError) for diagnostic use.
    /// </summary>
    public static (bool IsClose, double MaxAbsError, double MaxRelError) CheckDetailed(
        ReadOnlySpan<double> actual, ReadOnlySpan<double> desired,
        double rtol = 1e-7, double atol = 0.0)
    {
        if (actual.Length != desired.Length)
            throw new ArgumentException("Arrays must have the same length.");

        double absRtol = Math.Abs(rtol);
        double absAtol = Math.Abs(atol);
        bool allClose = true;
        double maxAbs = 0.0, maxRel = 0.0;

        for (int i = 0; i < actual.Length; i++)
        {
            double d = desired[i];
            double lower = Math.Min(d * (1.0 - absRtol), d - absAtol);
            double upper = Math.Max(d * (1.0 + absRtol), d + absAtol);
            if (actual[i] < lower || actual[i] > upper) allClose = false;

            double absErr = Math.Abs(actual[i] - d);
            double relErr = d != 0.0 ? absErr / Math.Abs(d) : absErr;
            if (absErr > maxAbs) maxAbs = absErr;
            if (relErr > maxRel) maxRel = relErr;
        }
        return (allClose, maxAbs, maxRel);
    }

    /// <summary>
    /// Asserts that <paramref name="actual"/> is within tolerance of <paramref name="desired"/>.
    /// Throws <see cref="IsoCloseException"/> on failure with a diagnostic message.
    /// </summary>
    public static void Assert(double actual, double desired, double rtol = 1e-7, double atol = 0.0,
        string? context = null)
    {
        if (!Check(actual, desired, rtol, atol))
        {
            double absErr = Math.Abs(actual - desired);
            double relErr = desired != 0.0 ? absErr / Math.Abs(desired) : absErr;
            string msg = context is null
                ? $"Isoclose assertion failed: actual={actual:G10} desired={desired:G10} |err|={absErr:G4} rel={relErr:G4} (rtol={rtol}, atol={atol})"
                : $"[{context}] Isoclose assertion failed: actual={actual:G10} desired={desired:G10} |err|={absErr:G4} rel={relErr:G4} (rtol={rtol}, atol={atol})";
            throw new IsoCloseException(msg);
        }
    }

    /// <summary>
    /// Asserts that every element of <paramref name="actual"/> is within tolerance.
    /// </summary>
    public static void Assert(ReadOnlySpan<double> actual, ReadOnlySpan<double> desired,
        double rtol = 1e-7, double atol = 0.0, string? context = null)
    {
        var (isClose, maxAbs, maxRel) = CheckDetailed(actual, desired, rtol, atol);
        if (!isClose)
        {
            string msg = context is null
                ? $"Isoclose assertion failed: maxAbsErr={maxAbs:G4} maxRelErr={maxRel:G4} (rtol={rtol}, atol={atol})"
                : $"[{context}] Isoclose assertion failed: maxAbsErr={maxAbs:G4} maxRelErr={maxRel:G4} (rtol={rtol}, atol={atol})";
            throw new IsoCloseException(msg);
        }
    }
}

/// <summary>Exception thrown by <see cref="Isoclose.Assert"/> on tolerance failure.</summary>
public sealed class IsoCloseException : Exception
{
    public IsoCloseException(string message) : base(message) { }
}
