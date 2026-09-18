# Ly-alpha impact of the archived redshift mismatch

The Downloads file `qso_z6.15_m-27.09_tq1e+00_los.npy` is byte-for-byte
identical to the earlier ATON reference, with SHA-256
`29c12c30582eba0cdedffcd839861e11ac2aad166eb7c463ee1de29830e8da3e`.
It therefore retains the z6.15-input / z6.604205-runtime mismatch.

## Matched calculation

24 distinct initial temperature profiles were sampled uniformly from the
464 unique profiles represented among the 924 retained archive rows, using
NumPy's default random generator with seed 20260906. Each selected skewer
was evolved to 1 Myr with M1450 = float32(-27.09), changing only the runtime
redshift from 6.604205 to 6.148904. Density, backgrounds, grid, helium rate
condition, and the other legacy settings were retained. Thus "corrected"
below means **corrected for this redshift mismatch**, not a validation of
all the legacy physics and unit conventions.

Comparisons use F = exp(-tau), without noise or instrumental broadening.
Both spectra use the same nominal z6.148904 analysis coordinate:

```text
r_pMpc = (distance_ckpc - distance_ckpc[256]) / (1000 * (1 + 6.148904))
```

This follows the user's `proximity_zone_p1d/p1d/simulation.py` convention.
The main statistics include only 0 <= r < 10 proper Mpc; redward pixels
are excluded. These are comparisons of matched native pixel indices on
the shared analysis coordinate, not separately remapped velocity grids.
The reported sample means weight the 24 distinct profiles equally; they
are not measurements from rerunning all 924 rows.

## Transmission profiles

| Distance (proper Mpc) | Archive mean F | Corrected mean F | Archive excess relative to corrected |
| --- | --- | --- | --- |
| 0–2 | 0.844171 | 0.837029 | 0.85% |
| 2–5 | 0.356758 | 0.344326 | 3.61% |
| 5–10 | 0.072652 | 0.068399 | 6.22% |
| 0–10 | 0.312729 | 0.305441 | 2.39% |

Across the selected blueward pixels:

- Mean absolute flux difference: **0.01417**, or 1.417 percentage points
  on the continuum-normalized transmission scale.
- RMS flux difference: 0.03802.
- 95th percentile of absolute pixel differences: 0.06901.
- 7.85% of pixels differ by more than 0.05; 2.71% differ by more than 0.1.
- Median largest pixel difference per skewer: 0.22530. The largest difference
  in the whole sample is 0.65521, on input LOS 266 (archive row 241).
- The per-skewer mean absolute difference ranges from 0.00452 to 0.02635;
  its 16th–84th percentile interval is 0.00699–0.02224.
- After a 20-Angstrom observed-frame boxcar, the mean absolute difference
  is 0.00787, excluding positions without full boxcar coverage.

The mean profile shifts by a few percent in transmission, while individual
absorption features can change much more. A single multiplicative flux
correction would not recover the changed feature shapes.

## Proximity-zone edge

Rp is measured at the first blueward pixel where the flux, smoothed over
20 observed Angstroms, is <= 0.1. Wavelengths are computed from the shared
distance coordinate with Planck18, as in the user's analysis code. The
smoothing implementation and all 48 Rp measurements were checked against
`proximity_zone_p1d/p1d/proximity.py` and matched exactly.

- Median paired Rp excess in the archive: **0.01607 proper Mpc** (one pixel).
- 16th–84th percentile paired excess: 0–0.15942 proper Mpc.
- Mean paired excess: 0.22098 proper Mpc, driven by two large outliers.
- Input LOS 334 / archive row 307: Rp changes from **6.139 to 4.211 pMpc**.
- Input LOS 587 / archive row 539: Rp changes from **9.610 to 7.248 pMpc**.

In the outliers, an earlier absorption trough crosses the 10% smoothed
flux threshold after the redshift correction. This changes which trough
sets Rp, so a modest flux change can cause a large change in the measured
edge. All 24 pairs have measurable boundaries within the saved spectra.

## Saved results and reproduction

All output paths below are relative to `outputs/lya_redshift_impact/`:

- `lya_profile_impact.png` and `.pdf`: mean profiles, difference, and a
  representative individual skewer.
- `lya_rp_outlier.png`: raw and smoothed profiles for the largest Rp change.
- `profile_impact.json`: full summary statistics and attachment identity.
- `per_sightline.csv`: individual flux metrics and Rp values.
- `corrected_redshift_subset.npy`: **24**, not 924, corrected skewers, shape
  `(24, 891, 7)`; ordering follows `sample_rows.txt` and `sample_input_los.txt`.
- `matched_flux_profiles.npz`: distance, native and smoothed fluxes for both
  calculations, with archive-row and global-input-LOS mappings.
- `run_manifest.json`: selection and simulation settings.
- `run_commands.json`: all executed commands.
- `correct_z6p15/`: individual text outputs and logs.

To repeat the RT runs using the selected IDs:

```bash
python tools/run_durovcikova_grid.py \
  --rte ./rte1d_opt \
  --input-dir /Users/ssatyavolu/projects/fireqso_dominika \
  --input-prefix f1kz6.15M3LOS --second-input-prefix l1kz6.15M3LOS \
  --output-dir outputs/lya_redshift_impact/correct_z6p15 \
  --z-sim 6.148904 --z-qso 6.15 \
  --m1450 -27.09 --legacy-float-magnitude \
  --tq-grid 1 --t-on 1 --fduty 1 \
  --n-cells 891 --n-los 500 --legacy-density 1 \
  --num-pixels-redward 256 --legacy-alphahe2-bug 1 \
  --boxsize-override 102441.268 \
  --los-file outputs/lya_redshift_impact/sample_input_los.txt
```

To recompute statistics and the main figure:

```bash
python tools/compare_lya_profiles.py \
  --reference /Users/ssatyavolu/projects/proximity-p1d/proximity_zone_p1d/data/simulations/ATON/qso_z6.15_m-27.09_tq1e+00_los.npy \
  --attachment /Users/ssatyavolu/Downloads/qso_z6.15_m-27.09_tq1e+00_los.npy \
  --output-dir outputs/lya_redshift_impact
```

The analysis helper requires NumPy, Astropy, and Matplotlib. They are
available in `/Users/ssatyavolu/miniconda3/bin/python` on this machine.
