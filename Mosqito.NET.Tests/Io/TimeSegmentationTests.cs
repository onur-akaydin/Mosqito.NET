using Mosqito.Io;
using Xunit;

namespace Mosqito.Tests.Io;

/// <summary>
/// Tests for TimeSegmentation.
/// Ported from MoSQITo tests/utils/test_time_segmentation.py.
/// </summary>
public class TimeSegmentationTests
{
    // Reference values from test_time_segmentation.py:
    //   signal: fs=48000, d=1s (48000 samples)
    //   time_segmentation(signal, fs, 8192, 2048, is_ecma=True)  → 24 blocks
    //   time_segmentation(signal, fs, 8192, 2048, is_ecma=False) → 20 blocks

    private static double[] MakeSine(int fs = 48000, double freq = 40.0, double splLevel = 60.0)
    {
        int n = fs; // 1 second
        double pRef = 2e-5;
        double amp = pRef * Math.Pow(10.0, splLevel / 20.0) * Math.Sqrt(2.0);
        double[] sig = new double[n];
        for (int i = 0; i < n; i++)
            sig[i] = amp * Math.Sin(2.0 * Math.PI * freq * i / fs);
        return sig;
    }

    [Fact]
    [Trait("Category", "Utils")]
    public void EcmaMode_Produces24Blocks_For1sAt48kHz()
    {
        double[] sig = MakeSine();
        var (blocks, timeAxis) = TimeSegmentation.Segment(sig, fs: 48000,
            nPerSeg: 8192, noOverlap: 2048, isEcma: true);

        Assert.Equal(24, blocks.GetLength(1));
        Assert.Equal(24, timeAxis.Length);
    }

    [Fact]
    [Trait("Category", "Utils")]
    public void StandardMode_Produces20Blocks_For1sAt48kHz()
    {
        double[] sig = MakeSine();
        var (blocks, timeAxis) = TimeSegmentation.Segment(sig, fs: 48000,
            nPerSeg: 8192, noOverlap: 2048, isEcma: false);

        Assert.Equal(20, blocks.GetLength(1));
        Assert.Equal(20, timeAxis.Length);
    }

    [Fact]
    [Trait("Category", "Utils")]
    public void BlockShape_IsNPerSeg_By_NBlocks()
    {
        double[] sig = new double[48000];
        var (blocks, _) = TimeSegmentation.Segment(sig, fs: 48000,
            nPerSeg: 4096, noOverlap: 2048, isEcma: false);

        Assert.Equal(4096, blocks.GetLength(0));
        Assert.True(blocks.GetLength(1) > 0);
    }

    [Fact]
    [Trait("Category", "Utils")]
    public void TimeAxis_IsMonotonicallyIncreasing()
    {
        double[] sig = MakeSine();
        var (_, timeAxis) = TimeSegmentation.Segment(sig, fs: 48000,
            nPerSeg: 4096, noOverlap: 2048, isEcma: false);

        for (int i = 1; i < timeAxis.Length; i++)
            Assert.True(timeAxis[i] > timeAxis[i - 1],
                $"Time axis not increasing at index {i}: {timeAxis[i - 1]:G4} → {timeAxis[i]:G4}");
    }

    [Fact]
    [Trait("Category", "Utils")]
    public void DefaultOverlap_IsHalfNPerSeg()
    {
        double[] sig = new double[48000];
        int nPerSeg = 4096;
        var (b1, _) = TimeSegmentation.Segment(sig, fs: 48000, nPerSeg: nPerSeg,
            noOverlap: nPerSeg / 2, isEcma: false);
        var (b2, _) = TimeSegmentation.Segment(sig, fs: 48000, nPerSeg: nPerSeg,
            noOverlap: null, isEcma: false);

        Assert.Equal(b1.GetLength(1), b2.GetLength(1));
    }

    [Fact]
    [Trait("Category", "Utils")]
    public void SegmentList_ReturnsCorrectNumberOfSegments()
    {
        double[] sig = MakeSine();
        var (segments, timeAxis) = TimeSegmentation.SegmentList(sig, fs: 48000,
            nPerSeg: 8192, noOverlap: 2048, isEcma: true);

        Assert.Equal(24, segments.Count);
        Assert.Equal(24, timeAxis.Length);
        Assert.All(segments, seg => Assert.Equal(8192, seg.Length));
    }
}
