# Reproducing qso_z5.11_m-26.76_tq1e-02_los.npy

## Result

The attached array

```text
/Users/ssatyavolu/Downloads/qso_z5.11_m-26.76_tq1e-02_los.npy
```

has SHA-256
`da24e26ed48e93467f42762ca8d11ba63a24314dd64d06ad4ee7a2dcf043ef6b`
and shape `(960, 891, 7)`.

It is reproducible, but it was not generated with the literal Section 3.1
redshift-rescaling procedure. Its recovered configuration is:

| Setting | Recovered value |
| --- | --- |
| Quasar represented by the magnitude | J1306+0356, observed z=6.0330 |
| M1450 | -26.76, rounded to legacy float32 precision |
| Lifetime and on-time | 0.01 Myr |
| Input fields | z5.11 `nH`, `nHI`, T, rho, peculiar velocity, and backgrounds |
| Exact input snapshot redshift | 5.109794 |
| Runtime RT and Ly-alpha redshift | **6.604205** |
| Output redshift label | 5.11 |
| Redshift density rescaling | none |
| Length rescaling | none |
| Saved cells and Ly-alpha offset | 891 and 256 |
| Retained sightlines | 960 of 1000 |

The filename's `z5.11` is the input snapshot/neutral-fraction label. It is
not the quasar redshift or the runtime cosmology used by this file.

This also establishes that the file uses z5.11 density and peculiar-velocity
fields. It does not combine z5.11 xHI/T with z5.95 density/velocity fields.

## Evidence

Five complete RT reruns (input LOS 0, 37, 39, 500, and 999) reproduce rows
from both input halves and rows on both sides of a discarded sightline. Across
those runs, the largest absolute xHI discrepancy is `2.25e-15`, the largest
temperature discrepancy is `1.00e-7 K`, and the largest absolute flux
discrepancy is `2.23e-9`.

An independent spectral reconstruction also checks 128 Ly-alpha bins in each
of all 960 rows. At runtime redshift 6.604205, its maximum relative optical-
depth discrepancy is `9.21e-8`; no row exceeds `1e-6`. Fitting redshift and
optical-depth amplitude to row 0 gives

```text
fitted spectral redshift = 6.604204999974313
fitted amplitude         = 0.999999999998365
RMS log-tau residual     = 6.02e-10
```

The current paper-style command at observed z=6.033 does not match the
attachment. For row 0, its mean absolute flux discrepancy is `0.07168` and
maximum discrepancy is `0.85777`. A trial at the closest z5.95 snapshot with
no rescaling also fails, with mean absolute flux discrepancy `0.02226`.

Machine-readable comparisons are under
`outputs/attached_m2676_validation/`.

## Exact archive-reproduction instructions

Run all commands from this repository. First build the solver:

```bash
make optimized
```

Infer the 960 retained input sightlines from the attached array:

```bash
/Users/ssatyavolu/miniconda3/bin/python tools/infer_retained_los.py \
  /Users/ssatyavolu/Downloads/qso_z5.11_m-26.76_tq1e-02_los.npy \
  --input-dir /Users/ssatyavolu/projects/fireqso_dominika \
  --snapshot-label 5.11 \
  --first-cell 591 --last-cell 891 \
  --output outputs/attached_m2676_validation/inferred_kept_indices.txt
```

Run the retained sightlines with the recovered legacy configuration:

```bash
/Users/ssatyavolu/miniconda3/bin/python tools/run_durovcikova_grid.py \
  --rte ./rte1d_opt \
  --input-dir /Users/ssatyavolu/projects/fireqso_dominika \
  --input-prefix f1kz5.11M3LOS \
  --second-input-prefix l1kz5.11M3LOS \
  --output-dir outputs/attached_m2676_validation/archive_candidate \
  --z-sim 6.604205 --z-qso 5.11 \
  --m1450 -26.76 --legacy-float-magnitude \
  --tq-grid 0.01 --t-on 0.01 --fduty 1 \
  --n-cells 891 --n-los 500 \
  --legacy-density 1 --density-multiplier 1 \
  --comoving-length-scale 1 \
  --num-pixels-redward 256 --legacy-alphahe2-bug 1 \
  --boxsize-override 102441.268 \
  --los-file outputs/attached_m2676_validation/inferred_kept_indices.txt \
  --workers 6 --resume
```

Assemble the text files in the recovered row order:

```bash
/Users/ssatyavolu/miniconda3/bin/python tools/assemble_outputs_to_npy.py \
  --input-dir outputs/attached_m2676_validation/archive_candidate \
  --indices-file outputs/attached_m2676_validation/inferred_kept_indices.txt \
  --output outputs/attached_m2676_validation/qso_z5.11_m-26.76_tq1e-02_los.npy \
  --z-qso 5.11 --m1450 -26.76 --tq 0.01 --t-on 0.01 --fduty 1 \
  --n-cells 891
```

The resulting array should have shape `(960, 891, 7)`.

## Instructions for new Section 3.1-style calculations

Archive reproduction and a new paper-style calculation are different tasks.
For a new model of J1306+0356, keep the runtime cosmology at its observed
redshift, z=6.0330. For each selected ATON snapshot, rescale its physical
density and length to z=6.0330. The high-level driver performs those steps:

```bash
/Users/ssatyavolu/miniconda3/bin/python tools/run_aton_section31_grid.py run \
  --rte ./rte1d_opt \
  --input-dir /Users/ssatyavolu/projects/fireqso_dominika \
  --output-root outputs/j1306_z6p033_m2676/section31_grid \
  --z-qso 6.033 --m1450 -26.76 \
  --snapshots 5.11 --tq-grid 0.01 \
  --los 0:999 --workers 6 --resume
```

That output is expected to differ from the attached file because it corrects
the attached file's recovered runtime-redshift and rescaling conventions.

## Limits relevant to future changes

1. The archived `.npy` file contains no run metadata. Its snapshot label,
   runtime redshift, density convention, and retained LOS mapping must be
   recovered from its numerical contents and surviving inputs.
2. Historical arrays cannot validate the literal paper rescaling when their
   recovered configuration did not use it. They validate the legacy solver
   and archive-reproduction path only.
3. The high-level Section 3.1 driver intentionally fixes runtime cosmology to
   the observed quasar redshift and applies snapshot-dependent density and
   length factors. It therefore will not numerically reproduce legacy arrays
   generated at fixed runtime z=6.604205 without those factors.
4. The legacy archives use float32 magnitude precision, the historical helium
   condition, a 256-pixel redward Ly-alpha offset, and filtered sightline sets.
   Omitting any of these can prevent a close numerical reproduction.
5. The local `z5.11` input family supplies all spatial fields used by this
   attachment, not only xHI and temperature. Future controlled experiments
   that deliberately mix xHI/T from one snapshot with density/velocity from
   another require a new explicit mixed-input generator and should use a
   separate label; the current drivers do not silently mix snapshots.
