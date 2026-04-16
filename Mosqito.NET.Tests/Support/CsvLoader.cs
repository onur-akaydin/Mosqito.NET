using System.Globalization;

namespace Mosqito.Tests.Support;

/// <summary>
/// Lightweight CSV loader for reference data files (no external package).
/// </summary>
public static class CsvLoader
{
    /// <summary>
    /// Reads a single-column CSV with optional header row.
    /// Matches numpy.genfromtxt(file, skip_header=1) behaviour.
    /// </summary>
    public static double[] LoadColumn(string path, int skipRows = 1, int column = 0, char delimiter = ',')
    {
        var lines = File.ReadAllLines(path)
            .Skip(skipRows)
            .Where(l => !string.IsNullOrWhiteSpace(l))
            .ToArray();

        double[] data = new double[lines.Length];
        for (int i = 0; i < lines.Length; i++)
        {
            string[] parts = lines[i].Split(delimiter);
            string cell = parts.Length > column ? parts[column].Trim() : "0";
            data[i] = double.Parse(cell, CultureInfo.InvariantCulture);
        }
        return data;
    }

    /// <summary>
    /// Reads an N-column CSV, returning a 2-D jagged array (rows × columns).
    /// </summary>
    public static double[][] LoadMatrix(string path, int skipRows = 0, char delimiter = ',')
    {
        return File.ReadAllLines(path)
            .Skip(skipRows)
            .Where(l => !string.IsNullOrWhiteSpace(l))
            .Select(l => l.Split(delimiter)
                          .Select(c => double.Parse(c.Trim(), CultureInfo.InvariantCulture))
                          .ToArray())
            .ToArray();
    }
}
