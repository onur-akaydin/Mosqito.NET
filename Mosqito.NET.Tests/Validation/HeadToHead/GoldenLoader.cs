using System.Text.Json;

namespace Mosqito.Tests.Validation.HeadToHead;

/// <summary>
/// Loads a single head-to-head golden JSON file and provides typed accessors.
/// Golden files live in ReferenceData/Golden/&lt;wav-stem&gt;/&lt;function&gt;.json.
/// NaN is stored as JSON null; error records contain an "error" key.
/// </summary>
public sealed class GoldenFile
{
    private readonly JsonDocument _doc;

    private GoldenFile(JsonDocument doc) => _doc = doc;

    public static GoldenFile Load(string wavStem, string functionName)
    {
        string path = Path.Combine("ReferenceData", "Golden", wavStem, functionName + ".json");
        if (!File.Exists(path))
            throw new FileNotFoundException($"Golden file not found: {path}");
        string json = File.ReadAllText(path);
        return new GoldenFile(JsonDocument.Parse(json));
    }

    /// <summary>Returns true if the Python generator encountered an error for this (wav, func) pair.</summary>
    public bool IsError => _doc.RootElement.TryGetProperty("error", out _);

    /// <summary>Error message when IsError is true, or null.</summary>
    public string? ErrorMessage =>
        IsError ? _doc.RootElement.GetProperty("error").GetString() : null;

    /// <summary>Gets a named scalar (double) from the "scalars" object.</summary>
    public double Scalar(string key)
    {
        var scalars = _doc.RootElement.GetProperty("scalars");
        var el = scalars.GetProperty(key);
        return el.ValueKind == JsonValueKind.Null ? double.NaN : el.GetDouble();
    }

    /// <summary>Gets a named 1-D array from the "arrays_1d" object. JSON null elements become NaN.</summary>
    public double[] Array1D(string key)
    {
        var arrays = _doc.RootElement.GetProperty("arrays_1d");
        if (arrays.ValueKind == JsonValueKind.Null)
            throw new KeyNotFoundException($"arrays_1d is null; key '{key}' not available");
        var arr = arrays.GetProperty(key);
        if (arr.ValueKind == JsonValueKind.Null) return Array.Empty<double>();
        return arr.EnumerateArray()
            .Select(e => e.ValueKind == JsonValueKind.Null ? double.NaN : e.GetDouble())
            .ToArray();
    }

    /// <summary>Gets a named metadata value (int) from the "metadata" object.</summary>
    public int MetaInt(string key)
    {
        var meta = _doc.RootElement.GetProperty("metadata");
        return meta.GetProperty(key).GetInt32();
    }

    /// <summary>Gets a named metadata value (string) from the "metadata" object.</summary>
    public string MetaString(string key)
    {
        var meta = _doc.RootElement.GetProperty("metadata");
        return meta.GetProperty(key).GetString() ?? "";
    }

    /// <summary>Gets a named bool[] from arrays_1d (encoded as a JSON array of booleans).</summary>
    public bool[] BoolArray(string key)
    {
        var arrays = _doc.RootElement.GetProperty("arrays_1d");
        if (arrays.ValueKind == JsonValueKind.Null)
            return Array.Empty<bool>();
        var arr = arrays.GetProperty(key);
        if (arr.ValueKind == JsonValueKind.Null) return Array.Empty<bool>();
        return arr.EnumerateArray()
            .Select(e => e.GetBoolean())
            .ToArray();
    }
}
