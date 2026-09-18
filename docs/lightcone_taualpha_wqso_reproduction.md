# Lightcone taualpha_wqso reproduction

The requested path is a directory of 250 sightlines:

```text
/Users/ssatyavolu/projects/meanfreepath/mean_free_path/lightcone/taualpha_wqso/
```

The concrete reference selected for the audit is
`z6.14LOS0tq1e+00mag-27.00_tep_1e+00Myr1.0fd.txt`, shape `(8071, 7)`.
Its SHA256 is
`75123042b5473d28eae1060147e0da7eba30e2e7a94339433835c6c231472f69`.

The labelled parameters are **z=6.14, M1450=-27.00, tq=1 Myr**.
The requested alternative **z=6.15, M1450=-27.09** is tested separately, as
is **z=6.148904** (the exact snapshot behind the FIRE filename label 6.15).
All three runs use the same verified lightcone initial fields. They do not
generate a new lightcone starting at a different quasar redshift.

## Complete-run result

The labelled **6.14, -27.00** parameters reproduce LOS 0. The requested
**6.15, -27.09** parameters do not reproduce the saved array. All three
configurations were evolved through **all 8071 cells**, followed by full
spectral synthesis, with the same 1 Myr lifetime.

Differences below are absolute continuum-normalized flux differences,
`F=exp(-tau)`, at corresponding saved pixel indices:

| Runtime z / M1450 | Mean absolute flux difference, all pixels | Maximum absolute flux difference, all pixels | Mean absolute flux difference, first 635 pixels |
| --- | ---: | ---: | ---: |
| **6.14 / -27.00** | **4.87e-10** | **2.76e-8** | **4.54e-10** |
| 6.15 / -27.09 | 0.001235 | 0.03733 | 0.009698 |
| 6.148904 / -27.09 | 0.001197 | 0.03744 | 0.009684 |

For the matching control, the maximum absolute xHI difference is **9.00e-16**,
and maximum absolute T difference is **5.20e-7 K**. Thus the match includes
the evolved gas fields, not only the flux. Independently recomputing the
optical depth from each candidate's stored final fields reproduces its C
flux to at most **1.85e-9** across the three runs and every output bin.

Changing runtime redshift also changes the velocity assigned to a pixel.
After interpolation onto a common 0–5000 km/s grid, the requested 6.15/-27.09
case has mean absolute difference **0.01402**, maximum **0.09334**. The
6.148904/-27.09 case gives **0.01383** and **0.08722**. These metrics include
small displacements of sharp absorption features; they use a different
comparison coordinate from the native-pixel table above.

The whole-lightcone mean is diluted by the long foreground segment where
the quasar has little effect, so the first-635-pixel and common-velocity
metrics are included alongside the whole-array numbers.

![Complete reproduction and requested-parameter comparison](../outputs/lightcone_taualpha_wqso/comparison_los0.png)

The completed configurations each took approximately 8 minutes in the
three-process comparison. Their measured runtime redshift is **6.14** for
the matching control; no 6.604205 setting is needed for this reference.

## Recovering the inputs

These are long lightcones, with initial gas conditions interpolated across
snapshots along the line of sight. The local redshift grid in
`zs8072_z6.14.txt` has 8071 entries, decreasing from 6.14 to 4.3351674572.
The relevant stitched snapshots are in `lightcone/stitchlos2/`.

For each cell, the reconstruction interpolates:

- Physical nH and peculiar velocity linearly in snapshot redshift.
- xHI, temperature, and comoving rho logarithmically in snapshot redshift.
- nHI as interpolated nH multiplied by interpolated xHI.

Against the supplied `f1kz6.14M3LOS_*` inputs, maximum relative differences are
3.33e-16 for nH and 4.22e-15 for nHI and temperature. The reconstructed initial
xHI and T also match the companion `taualpha_noqso/z6.14LOS0.txt` fields to
their text precision.

The root-directory `f1kz6.14M3LOS_rho.txt` and `..._vpec.txt` each have **8192**
rows and retain an older stitching layout. Using them directly fails farther
along the sightline. Their matching values are regenerated from `stitchlos2`
using the same interpolation as the verified nH/nHI/T. The original root
files are preserved.

The surviving lightcone input generator uses nH and nHI conversion without
the extra X=0.76 factor that appeared in the old FIRE inputs. No additional
division by 0.76 is applied in the reproducing calculation. The original
background-ionization and heating files are retained.

These input fields differ from the FIRE `f1kz6.15M3LOS_*` fields. In particular,
the lightcone inputs are not a single fixed-z=6.15 snapshot copied along the
entire skewer.

## RT and spectral conventions

| Setting | Value |
| --- | --- |
| Labelled runtime redshift | 6.14 |
| Labelled magnitude | -27.00 |
| Lifetime / on-time / duty cycle | 1 Myr / 1 Myr / 1 |
| Evolved cells | 8071 |
| Cell spacing / inner radius | 114.88970588235293, existing comoving units |
| Spectral boxsize override | 927389.7 |
| Extra redward spectral pixels | **0**, unlike the FIRE archive's 256 |
| Additional RT density multiplier | 1 |
| Background files | Original lightcone inputs |
| Legacy helium-coefficient condition | Retained |

The boxsize and cell count must be supplied separately. A fit to the no-quasar
reference gives boxsize 927389.699999999. Using `8071 * cellsize` changes the
velocity spacing and does not reproduce the stored spectra.

The measured spectral convention uses a fixed runtime z=6.14, despite the
changing local redshift of the input fields. Holding the reference's evolved
xHI/T fixed and varying this spectral redshift rejects 6.148904, 6.15, and
6.604205 as exact reproductions. The initial sampled spectral trials are
saved in `spectral_trials_los0.json`.

The first 891-cell gas-only probe already reproduced the archived xHI to
5.31e-16 and T to 5.20e-7 K at the labelled parameters. The probe's optical
depth is deliberately excluded from that comparison because the rest of the
lightcone was truncated. The complete runs remove that truncation.

## Code and repeatability

The solver previously had fixed 2048-cell input and background arrays. Their
capacity is now controlled by `RTE_MAX_CELLS`, defaulting to **2048**. The audit
compiles a separate executable with `-DRTE_MAX_CELLS=8192`, preserving the
ordinary build's capacity and physics. Both builds reject out-of-range cell
counts. The compiler retains the two existing warnings for the deliberately
selectable legacy helium comparison.

Run from the 1D-radiative-transfer repository:

```bash
/Users/ssatyavolu/miniconda3/bin/python tools/reproduce_lightcone_wqso.py --los 0
/Users/ssatyavolu/miniconda3/bin/python tools/summarize_lightcone_wqso.py --los 0
```

The first command reconstructs the inputs, builds the long-sightline executable,
and runs the labelled parameters and both requested-redshift interpretations.
The second writes the comparison report and figures. Optical-depth synthesis
is independently checked using a Python implementation on every output bin.

All artifacts are under `outputs/lightcone_taualpha_wqso/`:

- `reference_los0.json`, `input_reconstruction_los0.json`: provenance and input checks.
- `full_manifest_los0.json`, `full_results_los0.json`: exact commands and results.
- `full_inputs/los0/`: the nine input columns used in the complete runs.
- `full_runs/*/los0/`: logs, text output, `(8071,7)` NumPy output, and comparisons.
- `benchmark_los0.json`, `comparison_arrays_los0.npz`: consolidated results.
- `comparison_los0.png/.pdf`: native-profile overlays and common-velocity residuals.

The audit concerns complete RT runs of **LOS 0**. It is not a rerun of all 250
sightlines in the directory, nor a regeneration of the underlying cosmological
snapshots.
