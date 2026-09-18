# Downloaded z5.76, M1450=-26.99, 1 Myr benchmark

The current solver reproduces the tested sightlines in the reference supplied
on 2026-09-07. This is a successful reproduction of this candidate analysis
output; the user's identification of it as a final paper input remains tentative.

## Reference and recovered configuration

```text
/Users/ssatyavolu/Downloads/qso_z5.76_m-26.99_tq1e+00_los.npy
shape: (945, 891, 7), float64
SHA256: 3cdf6d41a196f54bc4486c5a33404cfa4908d0fdcb01b196a2cf0a3ae02f26c7
```

| Setting | Reproducing value |
| --- | --- |
| Input snapshot | `dom_zloc_M3_spec_z5.756.dat`, z=5.755954 |
| Input grids | `f1kz5.76M3LOS_*` and `l1kz5.76M3LOS_*` |
| Runtime redshift for RT and spectral synthesis | **6.604205** |
| Density rescaling relative to original snapshot inputs | **None: factor 1** |
| Initial background ionization/heating | Original snapshot backgrounds |
| Magnitude | -26.99, rounded to float32 as in the legacy executable |
| Lifetime / on-time / duty cycle | 1 Myr / 1 Myr / 1 |
| Evolved cells / redward spectral pixels | 891 / 256 |
| Spectral boxsize override | 102441.268 |
| Existing legacy density and helium settings | Retained (`legacy_density=1`, `legacy_alphahe2_bug=1`) |

The filename's redshift labels the input snapshot. The reproducing runtime
redshift is different, and the chemistry densities have not been rescaled
to it. These are recovered numerical settings, not an inference about why
they were chosen historically.

## Verification

Four complete RT runs reproduce the evolved fractions, temperature and optical
depth. They sample both input batches and widely separated retained rows:

| Reference row | Global input LOS | Maximum absolute flux error | Maximum absolute temperature error (K) |
| --- | --- | --- | --- |
| 0 | 0 | 6.23e-10 | 2.00e-7 |
| 314 | 335 | 3.19e-9 | 8.40e-7 |
| 629 | 672 | 1.51e-10 | 4.60e-7 |
| 944 | 999 | 9.30e-10 | 9.50e-7 |

Across those runs, the maximum absolute xHI error is **8.31e-16** and maximum
relative temperature error is **1.81e-11**. The comparison is numerical, not
bitwise equality. Relative differences in extremely small species fractions
can be larger; all seven columns' detailed errors are in `benchmark.json`.

Separately, all **945 rows and all 891 spectral bins** were reconstructed
using the reference's own final xHI/T and the original rho/velocity inputs:
841,995 saved bins in total. The maximum relative tau error is **1.21e-7**;
maximum absolute flux error is **1.75e-8**. This independently verifies the
spectral configuration and inferred retained-row mapping throughout the file.
The retained mapping was inferred by monotonic alignment of far-field initial
T/xHI and then validated spectrally. Duplicate input field content does not
uniquely identify which historical duplicate LOS ID was retained.

The original `.dat` fields were also used to regenerate all nine input quantities
for the four sampled LOS. T and peculiar velocities agree exactly. Number-density
and rho differences are below 5e-16 relative; background differences are below
1.44e-12. In particular, the physical nH conversion independently gives
**z=5.755954** in every sampled input. Density rescaling was therefore not hidden
in the sampled text inputs.

The full 945-row RT array has **not** been regenerated. The four complete RT runs
plus all-row spectral reconstruction establish the tested scope of this benchmark.
Only the supplied 1 Myr reference is validated here.

## Density-rescaling controls

### Direct recheck of the paper

[Table 1 in the 2024 arXiv PDF, page 5](https://arxiv.org/pdf/2401.10328v1#page=5)
identifies M1450=-26.99 with J1319+0950 at **z=6.1347 +/- 0.0005**, measured
from [C II]. Section 3.1 (PDF page 9) explicitly distinguishes the quasar
redshift z from the input snapshot redshift z_sim. Its stated density factor
for this downloaded file's snapshot is
`((1+6.1347)/(1+5.755954))^3 = 1.1777879269`; using the nearest available
snapshot 6.148904 as quasar redshift gives 1.1848362795.

There is a separate ambiguity in the literal length prescription: Section
3.1 prints the physical-length multiplier `(1+z)/(1+z_sim)`. If this is
intended as direct conversion of a fixed comoving skewer's proper lengths
from snapshot to quasar redshift, `L_proper=L_comoving/(1+z)` would instead
give `(1+z_sim)/(1+z)`. Both the v1 PDF text and v2 HTML contain the printed
factor. This may be a wording/typographical issue or require clarification
of the actual coordinate transformation; the historical implementation has
not been recovered. It must not be silently treated as a verified length
convention or as proof of an error in the published calculation. It does
not change the independently checked quasar redshift or stated density
factor. A subsequent [literal Section 3.1 experiment](literal_section31_comparison.md)
tests both printed factors with zsim=5.76 and zQSO=6.14 on four complete RT
sightlines and compares both background treatments directly with the download.

Section 3.3 further states that snapshot redshift acts as a proxy for the
volume-averaged neutral fraction. Thus a filename snapshot label differing
from the quasar redshift is expected in that method; the unexplained
difference here is the recovered numerical configuration, not that naming
convention itself.

### Numerical controls

Two additional LOS 0 runs multiply both initial nH and nHI, preserving cellwise
xHI, and otherwise retain the same legacy settings and original backgrounds:

| Runtime redshift | nH/nHI factor | Mean absolute flux error versus download | Maximum absolute flux error |
| --- | --- | --- | --- |
| 6.604205 | 1 (reproduction) | 5.79e-12 | 6.23e-10 |
| 6.604205 | 1.4259404783 | 0.03853 | 0.20215 |
| 6.148904 | 1.1848362795 | 0.02151 | 0.14064 |

These controls show the sensitivity to those specific alternative configurations.
They are **not** presented as a validated reconstruction of the paper's method:
backgrounds remain fixed and the existing legacy physics conventions are retained.
The strong evidence for the download's actual configuration is the match of its
full evolved fields to the unrescaled run, together with the `.dat` input check.

## Repeat the four-sightline RT benchmark

From the repository root, using the NumPy-enabled Python environment:

```bash
make optimized
/Users/ssatyavolu/miniconda3/bin/python tools/run_durovcikova_grid.py \
  --rte ./rte1d_opt \
  --input-dir /Users/ssatyavolu/projects/fireqso_dominika \
  --input-prefix f1kz5.76M3LOS \
  --second-input-prefix l1kz5.76M3LOS \
  --output-dir outputs/downloaded_m26p99_z5p76/reproduced_from_original_grids \
  --z-sim 6.604205 --z-qso 5.76 \
  --m1450 -26.99 --legacy-float-magnitude \
  --tq-grid 1 --t-on 1 --fduty 1 \
  --n-cells 891 --n-los 500 \
  --legacy-density 1 --num-pixels-redward 256 \
  --legacy-alphahe2-bug 1 --boxsize-override 102441.268 \
  --los 0 335 672 999
```

In this legacy mode, the current driver's `--z-sim` sets the **runtime** redshift,
while `--z-qso` supplies the output label. The input filenames and their physical
density conversion remain at 5.755954. These argument names must not be read as a
physically consistent snapshot/quasar rescaling in this reproduction command.

The executed benchmark extracted individual input columns into separate files;
those inputs and the exact executed argument lists are preserved in the artifacts.

To reconstruct all saved spectra using the reference's own final fields:

```bash
/Users/ssatyavolu/miniconda3/bin/python tools/audit_aton_redshifts.py \
  --archive-dir /Users/ssatyavolu/Downloads \
  --input-dir /Users/ssatyavolu/projects/fireqso_dominika \
  --indices-file outputs/downloaded_m26p99_z5p76/inferred_kept_indices.txt \
  --label 5.76 --m1450 -26.99 --lifetimes 1 \
  --trial-redshifts 6.604205 5.755954 6.1347 6.148904 \
  --output outputs/downloaded_m26p99_z5p76/all_rows_spectra.json
```

## Artifacts

All artifacts are under `outputs/downloaded_m26p99_z5p76/`:

- `benchmark.json`: reference identity, configuration, full RT sample errors and spectral audit.
- `reproduced_sample.npy`: four reproduced rows, in reference-row order 0, 314, 629, 944.
- `matched_samples.npz`: reference/reproduced arrays and explicit row/input mappings.
- `reproduction_comparison.png` and `.pdf`: overlaid flux profiles and residuals.
- `input_dat_verification.json`: original binary input comparison.
- `row0_rt_manifest.json`, `sampled_rt_manifest.json`: exact commands and timing.
- `row0_rt_comparisons.json`: comparisons of the two rescaling controls and reproduction.
- `inferred_kept_indices.txt`, `inferred_deleted_indices.txt`, `alignment.json`: row mapping.
- Individual case folders: extracted inputs, RT output text and logs.
