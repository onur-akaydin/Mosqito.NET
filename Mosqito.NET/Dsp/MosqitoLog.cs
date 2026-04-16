namespace Mosqito.Dsp;

/// <summary>
/// Library-wide logging callback. Default is silent; wire up to emit warnings:
/// <code>MosqitoLog.OnWarning = Console.WriteLine;</code>
/// </summary>
public static class MosqitoLog
{
    public static Action<string>? OnWarning;
    internal static void Warn(string msg) => OnWarning?.Invoke(msg);
}
