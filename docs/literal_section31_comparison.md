# Literal Section 3.1 comparison: zsim=5.76, zQSO=6.14

The density and physical-length transformations printed in
[Ďurovčíková et al. (2024), Section 3.1](https://arxiv.org/html/2401.10328v1#S3.SS1)
were applied before a fresh thermochemistry and Lyα calculation. On four matched
sightlines, this does **not** reproduce the downloaded array. Both treatments
of the unspecified background-rate rescaling give the same conclusion.

This experiment tests the printed redshift transformations within the current,
previously benchmarked solver. It retains that solver's other physics and unit
conventions; it is not an independently reconstructed implementation of every
detail of the publication or its inference pipeline.

## Configuration

- Reference: `/Users/ssatyavolu/Downloads/qso_z5.76_m-26.99_tq1e+00_los.npy`
- SHA256: `3cdf6d41a196f54bc4486c5a33404cfa4908d0fdcb01b196a2cf0a3ae02f26c7`
- Reference shape: `(945, 891, 7)`; rerun rows: **0, 314, 629, 944**.
- Corresponding original global input sightlines: **0, 335, 672, 999**.
- Requested redshifts: **zsim=5.76, zQSO=6.14**, with no substitution of the
  exact observed redshift or a neighbouring snapshot redshift.
- Magnitude: **M1450=-26.99**, represented as float32 as in the validated legacy
  run; lifetime and on-time **1 Myr**, duty cycle **1**.
- Initial temperature, xHI, and peculiar velocity of each cell are unchanged.
- Both initial nH and nHI are multiplied together before RT. Their ratio is
  preserved, but the evolved xHI is allowed to change.
- Runtime expansion factor, thermochemistry cosmological terms, and spectral
  Hubble conversion all use **6.14** in the rescaled cases.

Let A=(1+6.14)/(1+5.76). The literal printed transformation is:

| Quantity | Multiplier |
| --- | ---: |
| Physical density relative to zsim=5.76 | A³ = **1.1782963910525566** |
| Physical lengths relative to zsim=5.76 | A = **1.0562130177514792** |
| Comoving lengths supplied to this solver at zQSO=6.14 | A² = **1.1155859388676865** |

The last row follows because the C solver converts comoving lengths to proper
lengths by dividing by `(1+runtime_z)`. Thus
`(A² / (1+6.14)) / (1 / (1+5.76)) = A`. The cell spacing, inner radius, and
spectral integration box length are all multiplied by A². The original spectral
boxsize `102441.268` becomes `114282.03814057629` in the solver's existing units.

The original text input densities were actually converted at **5.755954**,
as checked against the binary snapshot. They are first normalized to the
explicitly requested **5.76**, a factor **1.0017977136384093**, and then multiplied
by A³. The total multiplier on the existing nH and nHI files is therefore
**1.1804146305448402**. This small normalization is recorded explicitly rather
than silently treating the rounded filename as the original conversion redshift.

The `rho` input is **comoving**, so it remains unchanged: its physical density
rescaling is supplied by the runtime scale factor in the optical-depth formula.
Multiplying this input by A³ as well would apply the density change twice.

## Background variants and retained conventions

Section 3.1 does not specify how the ionizing/heating background changes under
the rescaling. Two variants were run on all four sightlines:

1. **Original background:** retain the source snapshot's background files.
2. **Recalculated background:** recompute GammaHI, GammaHeI, GammaHeII, and
   background heating from the rescaled nH/nHI and unchanged initial temperature,
   using the existing `backgroundgamma.py` formulas ported in the input generator.

Both use the established legacy choices (`legacy_density=1` for already prepared
physical nH, `legacy_alphahe2_bug=1`, X=0.76, 891 cells, 256 redward spectral bins,
existing SED, and spectral cosmology h=0.678, OmegaM=0.308). These choices are
stated so the experiment is reproducible; they are not inferred uniquely from
the paper's prose. Changing those choices is outside this rescaling comparison.

## Lyα results

Flux is `F=exp(-tau)`. All differences below are **absolute normalized-flux
differences**, not relative percentages. Values describe only the four rerun
sightlines, not the full 945-row ensemble.

| Comparison | Original background | Recalculated background |
| --- | ---: | ---: |
| Mean absolute difference, corresponding saved pixels (all 891) | 0.05160 | 0.04956 |
| Maximum absolute difference, corresponding saved pixels | 0.52291 | 0.52327 |
| Mean absolute difference, individual profiles in 600 km/s bins, 0–4800 km/s blueward | 0.07289 | 0.07259 |
| Maximum absolute difference, individual profiles in those bins | 0.36838 | 0.36744 |
| Mean absolute difference of the four-sightline mean, same 600 km/s bins | **0.05000** | **0.04952** |
| Maximum absolute difference of the four-sightline mean, same bins | **0.12942** | **0.12842** |

The mean-profile difference is thus about **five percentage points** on the
continuum-normalized flux scale after this binning. The rescaled mean is lower
by 0.04319 (original background) or 0.04268 (recalculated background) averaged
over 0–4800 km/s.

Raw absorption features also shift because velocity spacing changes from
**26.01815814** to **28.14046507 km/s** per bin. On a common velocity grid over
0–5000 km/s blueward, mean absolute differences are 0.18595 and 0.18574 before
600 km/s binning. These larger values include displacement of narrow features;
they should not be confused with a 19% change in the mean transmitted flux.
The spectral grid uses the exact Hubble spacing in `src/lya_spectra.c`, removing
the 256-bin redward offset. The reference grid uses the previously verified
spectral runtime redshift 6.604205; the new grid uses 6.14 and the rescaled box.
Common-grid comparisons interpolate only within overlapping support. The
600 km/s values integrate piecewise-linear flux on each spectrum's native grid.

Changing the background treatment alone changes flux by only **0.000571** on
average (maximum **0.003614**) over the common 0–5000 km/s interval. It does not
remove the discrepancy.

![Mean Lyα profiles in 600 km/s bins](../outputs/section31_zsim5p76_zqso6p14_m26p99/mean_lya_600kms.png)

![Individual Lyα profiles and residuals](../outputs/section31_zsim5p76_zqso6p14_m26p99/lya_comparison.png)

## Verification and artifacts

- The optional C `comoving_length_scale` argument defaults to 1. A fresh control
  run with the established original settings reproduced reference row 0 to a
  maximum absolute flux difference of **6.23e-10** after this code change.
- The rescaled spectra were independently reconstructed using the new final
  gas fields, original comoving rho/vpec, runtime z=6.14, and the scaled box.
  The largest absolute flux discrepancy against C was **1.11e-9** across the
  eight rescaled runs.
- Initial xHI invariance and unchanged T/vpec/rho are checked before execution.
- `make optimized` succeeded with the two existing warnings for the deliberately
  selectable legacy helium comparison. Invalid/nonfinite length scale is rejected.
- No original legacy input or downloaded array was overwritten.

Repeat from the repository root:

```bash
make optimized
/Users/ssatyavolu/miniconda3/bin/python tools/compare_literal_section31.py
```

To regenerate just diagnostics and figures from completed RT outputs:

```bash
/Users/ssatyavolu/miniconda3/bin/python tools/compare_literal_section31.py --analyze-only
```

Artifacts are in `outputs/section31_zsim5p76_zqso6p14_m26p99/`:

- `configuration.json`, `manifest.json`, `completed_jobs.json`: factors, exact
  commands, provenance, and timings.
- `comparison.json`: all gas-field, flux, regression, and independent-spectrum metrics.
- `fixed_background_sample_los.npy`, `recomputed_background_sample_los.npy`:
  fresh outputs, each shape `(4, 891, 7)`, in reference-row order 0, 314, 629, 944.
- `comparison_arrays.npz`: matched reference and new arrays, row mapping,
  native/common velocity grids, and rebinned flux arrays.
- `lya_comparison.png/.pdf`, `mean_lya_600kms.png/.pdf`: comparison figures.
- Per-variant/per-row subdirectories contain the prepared inputs, output text,
  and complete solver logs.
