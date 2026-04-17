# Mosqito.NET

A complete, native C# (.NET 8) port of the [MoSQITo](https://github.com/Eomys/MoSQITo) Python psychoacoustics library.

Mosqito.NET implements the same sound quality metrics as the original Python library — loudness, roughness, sharpness, tonality, and speech intelligibility — with a flat `Sq.*` façade that mirrors MoSQITo's public API.

## Acknowledgements

Heartfelt thanks to the authors, maintainers, and contributors of the original MoSQITo Python library — in particular the team at [Eomys](https://eomys.com) — for their open, well-documented reference implementations of ISO 532-1, ECMA-418-2, DIN 45692, ANSI S3.5, ECMA-74, and the Daniel & Weber roughness model. Mosqito.NET would not exist without their work, and all scientific credit for the underlying algorithms belongs to them.

## Disclaimer — AI-assisted port

Mosqito.NET was produced largely with AI assistance (Claude), working from the original Python MoSQITo sources as the reference implementation. The test suite (65 ported unit tests plus 168 head-to-head comparisons against Python golden outputs, all passing) is the primary guarantor of correctness, but users should still validate results against their own reference data before relying on Mosqito.NET in production or safety-critical contexts. Bug reports and corrections are very welcome.

## Attribution, trademark, and endorsement

- **Mosqito.NET is an independent port and is not affiliated with or endorsed by the MoSQITo project or Eomys.**
- "MoSQITo" is a name used by the original project; "Mosqito.NET" refers only to this C# port and is not a claim over the original name.
- The Apache-2.0 license (§6) does not grant trademark rights; this port is licensed under Apache-2.0 for the code only.

## Standards implemented

| Standard | Metric |
|---|---|
| ISO 532-1:2017 | Zwicker loudness (stationary + time-varying) |
| ECMA-418-2:2022 | Loudness (hearing model), Roughness |
| DIN 45692:2009 | Sharpness (DIN, Bismarck, Aures, Fastl weightings) |
| Daniel & Weber (1997) | Roughness (time-domain and frequency-domain) |
| ECMA-74 Annex D | Tone-to-Noise Ratio (TNR), Prominence Ratio (PR) |
| ANSI S3.5-1997 | Speech Intelligibility Index (SII) |
| IEC 61260 | N-th octave spectral analysis |

## Quick start

```csharp
using Mosqito;
using Mosqito.Io;

// Load a WAV file
var (signal, fs) = WavLoader.Read("my_signal.wav");

// Zwicker stationary loudness (ISO 532-1)
var loudness = Sq.LoudnessZwst(signal, fs);
Console.WriteLine($"N = {loudness.N:F3} sone");

// DIN sharpness
double S = Sq.SharpnessDinSt(signal, fs);
Console.WriteLine($"S = {S:F3} acum");

// ECMA-418-2 roughness
var (R, RTime, RSpec, bark, time) = Sq.RoughnessEcma(signal, fs);
Console.WriteLine($"R = {R:F3} asper");

// Tone-to-noise ratio (ECMA-74)
var tnr = Sq.TnrEcmaSt(signal, fs);
Console.WriteLine($"T_total = {tnr.TTotal:F3}");
```

## API overview

All functions are available through the `Sq` static façade class. Each function also has a strongly-typed class in the corresponding `Mosqito.SqMetrics.*` namespace.

### Sound level meter

| Method | Description |
|---|---|
| `Sq.CompSpectrum` | Complex (amplitude or dB) FFT spectrum |
| `Sq.NoctSpectrum` | N-th octave spectrum in dB |
| `Sq.NoctSynthesis` | N-th octave synthesis (amplitude bands) |
| `Sq.FreqBandSynthesis` | Arbitrary frequency-band integration |

### Loudness

| Method | Description |
|---|---|
| `Sq.LoudnessZwst` | ISO 532-1 stationary loudness from time-domain signal |
| `Sq.LoudnessZwstFreq` | ISO 532-1 stationary loudness from spectrum |
| `Sq.LoudnessZwstPerSeg` | ISO 532-1 stationary loudness, segmented |
| `Sq.LoudnessZwtv` | ISO 532-1 time-varying loudness |
| `Sq.LoudnessEcma` | ECMA-418-2 loudness |
| `Sq.EqualLoudnessContours` | ISO 226 equal-loudness contours |
| `Sq.SoneToPhon` | Sone ↔ phon conversion |

### Roughness

| Method | Description |
|---|---|
| `Sq.RoughnessDw` | Daniel & Weber roughness from time-domain signal |
| `Sq.RoughnessDwFreq` | Daniel & Weber roughness from spectrum |
| `Sq.RoughnessEcma` | ECMA-418-2 roughness |

### Sharpness

| Method | Description |
|---|---|
| `Sq.SharpnessDinSt` | DIN 45692 stationary sharpness |
| `Sq.SharpnessDinTv` | DIN 45692 time-varying sharpness |
| `Sq.SharpnessDinPerSeg` | DIN 45692 sharpness, segmented |
| `Sq.SharpnessDinFreq` | DIN 45692 sharpness from spectrum |

### Tonality

| Method | Description |
|---|---|
| `Sq.TnrEcmaSt` | ECMA-74 Tone-to-Noise Ratio, stationary |
| `Sq.TnrEcmaFreq` | ECMA-74 TNR from spectrum |
| `Sq.TnrEcmaPerSeg` | ECMA-74 TNR, segmented |
| `Sq.PrEcmaSt` | ECMA-74 Prominence Ratio, stationary |
| `Sq.PrEcmaFreq` | ECMA-74 PR from spectrum |
| `Sq.PrEcmaPerSeg` | ECMA-74 PR, segmented |

### Speech intelligibility

| Method | Description |
|---|---|
| `Sq.SiiAnsi` | ANSI S3.5 Speech Intelligibility Index from signal |
| `Sq.SiiAnsiFreq` | ANSI S3.5 SII from spectrum |

## Performance

All public APIs were accelerated via `Parallel.For`, SIMD (`System.Numerics.Vector<double>`), `ArrayPool<T>`, plan caches, and allocation removal across two optimisation passes. Numbers below are from BenchmarkDotNet on an AMD Ryzen 7 4800H (8 cores / 16 threads), .NET 8 AVX2, comparing the unoptimised initial implementation to the current build.

### Loudness

| Benchmark | Before | After | Speedup | Alloc Δ |
|---|---:|---:|:---:|:---:|
| LoudnessZwst (synthetic) | 29.94 ms | 5.70 ms | **5.3×** | −65% |
| LoudnessZwst (WAV file) | 312.66 ms | 61.45 ms | **5.1×** | −69% |
| LoudnessZwstFreq (WAV file) | 572.06 ms | 3.57 ms | **160×** | −99% |
| LoudnessZwtv (synthetic) | 39.65 ms | 15.08 ms | **2.6×** | −7% |
| LoudnessZwtv (WAV file) | 475.30 ms | 201.90 ms | **2.4×** | — |
| LoudnessEcma (synthetic) | 84.96 ms | 34.29 ms | **2.5×** | −44% |
| LoudnessEcma (WAV file) | 1,047.49 ms | 358.89 ms | **2.9×** | −18% |

`LoudnessZwstFreq` gains are dominated by plan-caching the critical-band filter-bank design (previously recomputed from scratch on every call).

### Roughness

| Benchmark | Before | After | Speedup | Alloc Δ |
|---|---:|---:|:---:|:---:|
| RoughnessDw (synthetic) | 3,732.9 ms | 687.4 ms | **5.4×** | — |
| RoughnessDw (WAV file) | 11,338.7 ms | 2,216.8 ms | **5.1×** | — |
| RoughnessDwFreq (synthetic) | 749.3 ms | 150.7 ms | **5.0×** | — |
| RoughnessDwFreq (WAV file) | 772.8 ms | 150.8 ms | **5.1×** | — |
| RoughnessEcma (synthetic) | 1,977.0 ms | 464.4 ms | **4.3×** | −24% |
| RoughnessEcma (WAV file) | 4,676.5 ms | 979.3 ms | **4.8×** | −28% |

### Sharpness

| Benchmark | Before | After | Speedup | Alloc Δ |
|---|---:|---:|:---:|:---:|
| SharpnessDinSt (synthetic) | 30.73 ms | 6.14 ms | **5.0×** | −63% |
| SharpnessDinSt (WAV file) | 154.82 ms | 30.31 ms | **5.1×** | −61% |
| SharpnessDinTv (synthetic) | 40.95 ms | 14.79 ms | **2.8×** | −6% |
| SharpnessDinTv (WAV file) | 217.34 ms | 84.07 ms | **2.6×** | −9% |

### Tonality (TNR / PR) and Speech Intelligibility (SII)

Single-call benchmarks are FFT-dominated; the parallel gains apply to the `PerSeg` variants on longer signals. Allocation reductions from `stackalloc` locals and removed LINQ are present in all call paths.

| Benchmark | Before | After | Alloc Δ |
|---|---:|---:|:---:|
| TnrEcmaSt (synthetic) | 34.17 ms | 34.20 ms | −18% |
| TnrEcmaSt (WAV file) | 33.92 ms | 35.45 ms | −19% |
| PrEcmaSt (synthetic) | 34.15 ms | 34.44 ms | −18% |
| PrEcmaSt (WAV file) | 34.54 ms | 34.36 ms | −19% |
| SiiAnsi (synthetic) | 31.30 ms | 30.19 ms | −20% |
| SiiAnsi (WAV file) | 240.46 ms | 239.31 ms | −15% |

## Dependencies

- [MathNet.Numerics](https://numerics.mathdotnet.com/) 5.0 — numerical routines (FFT, interpolation, filtering)
- [NAudio](https://github.com/naudio/NAudio) 2.2 — WAV file I/O

## Test results

The library ships a comprehensive test suite with **233 tests, all passing** (as of 2026-04-16, run time ~7 minutes).

```
Test Run Successful.
Total tests:  233
Passed:       233
Failed:         0
```

### Test tiers

#### Tier 1 — Ported MoSQITo unit tests (65 cases)

Direct ports of the original Python MoSQITo test suite, validating each module against ISO/ECMA/DIN reference values.

| Test class | Cases | What it validates |
|---|---|---|
| `LoudnessZwstTests` | 5 | ISO 532-1 stationary loudness — time-domain, freq-domain, per-segment, 44.1 kHz resample, 1/3-octave input |
| `LoudnessZwtvTests` | 3 | ISO 532-1 time-varying loudness — dimensions, peak detection, time-average convergence |
| `LoudnessEcmaTests` | 3 | ECMA-418-2 loudness — dimensions, equal-loudness contour, 1 kHz vs 5 kHz comparison |
| `SharpnessDinTests` | 6 | DIN 45692 sharpness — all 4 weightings, freq-domain path, per-segment |
| `RoughnessDwTests` | 3 | Daniel-Weber roughness — AM tone ~1 asper (time & freq domain), dimensions |
| `RoughnessEcmaTests` | 3 | ECMA-418-2 roughness — AM tone 1 kHz/70 Hz/60 dB ~1 asper, R_time positive, dimensions |
| `TnrEcmaTests` | 5 | ECMA TNR — 2-tone detection, t_total value, per-segment, freq path, non-prominent tones |
| `PrEcmaTests` | 5 | ECMA PR — same 5 patterns as TNR |
| `SiiAnsiTests` | 4 | ANSI S3.5 SII — octave-method reference 0.504, level vs freq entry equivalence, all methods |
| `NoctTests` | 8 | N-th octave spectrum & synthesis — 1/3 & full octave, 94 dB SPL @ 1 kHz, 3 dB band addition |
| `TimeSegmentationTests` | 6 | I/O segmentation — block count, shape, overlap, monotonic time axis, ECMA 24-block mode |
| `SqFacadeTests` | 14 | Facade smoke tests for every public `Sq.*` method |

#### Tier 2 — Head-to-head Python golden tests (168 cases)

Each of the **24 public functions** is run against **7 WAV files** and compared to outputs pre-generated by the original Python MoSQITo library. Tolerance is **rtol = 1%** for psychoacoustic metrics, **atol = 0.5 dB** for dB-domain arrays.

**Functions covered:** `loudness_zwst`, `loudness_zwst_freq`, `loudness_zwst_perseg`, `loudness_zwtv`, `loudness_ecma`, `sharpness_din_st`, `sharpness_din_tv`, `sharpness_din_perseg`, `sharpness_din_freq`, `roughness_dw`, `roughness_dw_freq`, `roughness_ecma`, `tnr_ecma_st`, `pr_ecma_st`, `tnr_ecma_freq`, `pr_ecma_freq`, `tnr_ecma_perseg`, `pr_ecma_perseg`, `sii_ansi`, `sii_ansi_freq`, `noct_spectrum`, `noct_synthesis`, `comp_spectrum`, `freq_band_synthesis`

**WAV signals used:**

| Signal | Sample rate | Calibration |
|---|---|---|
| `white_noise_200_2000_Hz_stationary` | 48 kHz | 0.01 |
| `white_noise_442_1768_Hz_stationary` | 48 kHz | 0.01 |
| `white_noise_442_1768_Hz_varying` | 48 kHz | 0.01 |
| `broadband_570` | 48 kHz | 1.0 |
| `Test signal 3 (1 kHz 60 dB)` | 44.1 kHz | 2√2 |
| `Test signal 5 (pinknoise 60 dB)` | 48 kHz | 2√2 |
| `Test signal 10 (tone pulse 1 kHz 10 ms 70 dB)` | 48 kHz | 2√2 |

### Running the tests

```bash
dotnet test
```

> Note: the Daniel-Weber roughness head-to-head tests are computationally intensive (~2–3 minutes for the long WAV files). The full suite completes in approximately 7 minutes.

## License

Mosqito.NET is derived from [MoSQITo](https://github.com/Eomys/MoSQITo) (Apache License 2.0) and is distributed under the same license. See [NOTICE.md](NOTICE.md) for full attribution.
