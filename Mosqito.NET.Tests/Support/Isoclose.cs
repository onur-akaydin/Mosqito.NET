using Xunit;

namespace Mosqito.Tests.Support;

/// <summary>
/// ISO 532-1 section 5.1 tolerance assertion for xUnit tests.
/// Mirrors MoSQITo's <c>isoclose(actual, desired, rtol, atol)</c>.
/// </summary>
public static class Isoclose
{
    /// <summary>
    /// Checks and asserts that <paramref name="actual"/> is within tolerance of
    /// <paramref name="desired"/>. Throws xUnit assertion failure on mismatch.
    /// </summary>
    public static void Assert(double actual, double desired, double rtol = 0.05, double atol = 0.1,
        string? context = null)
    {
        double absRtol = Math.Abs(rtol);
        double absAtol = Math.Abs(atol);
        double lower = Math.Min(desired * (1.0 - absRtol), desired - absAtol);
        double upper = Math.Max(desired * (1.0 + absRtol), desired + absAtol);
        
        string msg = context is null
            ? $"Isoclose: actual={actual:G10} not in [{lower:G6}, {upper:G6}] (desired={desired:G10}, rtol={rtol}, atol={atol})"
            : $"[{context}] Isoclose: actual={actual:G10} not in [{lower:G6}, {upper:G6}] (desired={desired:G10}, rtol={rtol}, atol={atol})";

        Xunit.Assert.True(actual >= lower && actual <= upper, msg);
    }

    /// <summary>
    /// Checks that every element of <paramref name="actual"/> is within tolerance of
    /// the corresponding element of <paramref name="desired"/>.
    /// </summary>
    public static void Assert(double[] actual, double[] desired, double rtol = 0.05, double atol = 0.1,
        string? context = null)
    {
        Xunit.Assert.Equal(desired.Length, actual.Length);
        double absRtol = Math.Abs(rtol);
        double absAtol = Math.Abs(atol);

        for (int i = 0; i < actual.Length; i++)
        {
            double lower = Math.Min(desired[i] * (1.0 - absRtol), desired[i] - absAtol);
            double upper = Math.Max(desired[i] * (1.0 + absRtol), desired[i] + absAtol);
            string msg = context is null
                ? $"Isoclose[{i}]: actual={actual[i]:G10} not in [{lower:G6}, {upper:G6}] (desired={desired[i]:G10})"
                : $"[{context}][{i}] Isoclose: actual={actual[i]:G10} not in [{lower:G6}, {upper:G6}] (desired={desired[i]:G10})";
            Xunit.Assert.True(actual[i] >= lower && actual[i] <= upper, msg);
        }
    }

    /// <summary>
    /// Returns <see langword="true"/> if <paramref name="actual"/> is within tolerance of
    /// <paramref name="desired"/> (no assertion thrown).
    /// </summary>
    public static bool Check(double actual, double desired, double rtol = 0.05, double atol = 0.1)
    {
        double lower = Math.Min(desired * (1.0 - Math.Abs(rtol)), desired - Math.Abs(atol));
        double upper = Math.Max(desired * (1.0 + Math.Abs(rtol)), desired + Math.Abs(atol));
        return actual >= lower && actual <= upper;
    }
}
