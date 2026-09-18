# FIRE z6.15 with matching runtime cosmology, with and without a quasar

This calculation changes the archived FIRE runtime redshift from 6.604205 to
**6.148904**, the exact redshift of the snapshot labelled **6.15**. It follows
the same conventions used to reproduce the z5.95 archive. The paired arrays
use the original archive's 924 retained sightlines in the same order.

## Outputs

The output directory is `outputs/corrected_z6p15_m2709/`:

- `qso_z6.15_m-27.09_tq1e+00_los.npy`: quasar at M1450=-27.09, evolved for 1 Myr.
- `noqso_z6.15_los.npy`: initial snapshot gas before quasar turn-on.
- `input_los_ids.txt`: original global input LOS ID for each array row.
- `representative_input_los_ids.txt`: the computed representative of each row.
- `manifest.json`: exact configuration, input paths and SHA256 hashes.
- `commands.json`: all executable commands and individual output locations.
- `results.json`, `qso_results.json`, `noqso_results.json`: completion, file
  hashes, and numerical validation. A combined `results.json` is written
  only after both complete arrays have been generated.

Each completed NumPy array has shape **(924, 891, 7)**, with columns:

```text
distance_ckpc, xHI, temperature_K, xHeI, xHeII, xHeIII, tau_Lya
```

The no-quasar calculation synthesizes the spectrum directly from the initial
hydrogen fraction and temperature. It does not evolve a background-only
model for another 1 Myr and has no quasar lifetime or magnitude dependence.
The snapshot's existing IGM ionization state is retained. Its helium columns
are the solver's initialization (HeII=HeIII=1e-20), not separately extracted
simulation helium fractions.

## Configuration and scope

| Setting | Value |
| --- | --- |
| Initial snapshot and runtime redshift | **6.148904** |
| Filename redshift label | 6.15 |
| Quasar magnitude | -27.09, rounded to float32 as in the legacy runner |
| Quasar lifetime / on-time / duty cycle | 1 Myr / 1 Myr / 1 |
| Density treatment | Original nH/nHI, no additional multiplier |
| Background ionization and heating | Original FIRE inputs |
| Legacy helium coefficient condition | Enabled |
| Cells | 891 |
| Extra redward spectral pixels | 256 |
| Spectral boxsize | 102441.268 |
| Comoving length multiplier | 1 |

The existing FIRE density conversion includes an extra X=0.76 factor in
both nH and nHI. **This factor is retained**, as in the successful z5.95
reproduction. Here "corrected" refers specifically to the runtime redshift;
this is not a simultaneous change of the density normalization or helium
physics. The spectral comoving rho is also retained independently of nH.

All nine original input fields are compared when finding duplicate sightlines.
There are **469 exactly distinct input sets** among the 924 rows. Each set
is computed once, and exact duplicates are inserted back into their original
positions. No approximate matching or averaging is used.

The selected LOS IDs come from the previously verified archive mapping in
`outputs/recover_z6p15/inferred_kept_indices_z6.15_m27.09_from_target.txt`.
This retains the archive's selection, including its duplicate rows; it does
not claim that the rows are 924 independent gas realizations.

## Coordinates

The seven-column layout preserves the archive's coordinate convention:
the gas columns use the distance in column 0, while tau has a **256-pixel
redward offset**. Gas fields and tau at the same array index therefore do
not have the same spatial origin. Both paired arrays use the same grid.

For the spectral velocity coordinate used by the C code:

```python
import numpy as np

z = 6.148904
a = 1 / (1 + z)
n = 891
box = 102441.268
wmax = (a * 3.2407789e-18 * 3.085678e21 / 1e5
        * np.sqrt(0.308 / a**3 + 0.692) * box / 0.678)
velocity_kms = (np.arange(n) - 256) * wmax / n
```

The usual analysis distance convention is
`(data[0,:,0] - data[0,256,0]) / (1000 * (1 + z))` in proper Mpc.

## Repeat

From the repository root, using Python with NumPy and SciPy installed:

```bash
make optimized
/Users/ssatyavolu/miniconda3/bin/python tools/generate_fire_snapshot_pair.py --workers 6
```

The batch validates each job and resumes completed ones only when the
configuration, executable, command and output hashes still match. Use a new
`--output-dir` for a different configuration. `--tq` changes the quasar
lifetime and on-time together; it does not change the no-quasar definition.

The C runner now accepts a final optional `initial_only` argument, default 0.
With `initial_only=1`, it skips time evolution and writes
`noqso_z6.15LOS<ID>.txt`. `tq=0` is allowed in this mode; a positive on-time
placeholder is still required. The existing grid runner and assembler also
accept `--initial-only`.

## Validation

Every generated sightline is checked for finite values, physical ionization
fractions, positive temperature, and nonnegative optical depth. An independent
Python spectral calculation checks all 891 saved bins of every distinct
sightline in each mode. The no-quasar fields are additionally checked against
the input nHI/nH and temperature.

The no-quasar array preserves initial xHI to **5.00e-14 absolute** and
temperature to **5.01e-13 relative**. Its maximum independent flux discrepancy
is **2.63e-9**. A fresh corrected-redshift quasar run for input LOS 5 is
**identical in all seven saved columns** to the pre-existing validated run,
confirming that the optional initial-only path leaves that calculation unchanged.
A fresh z5.95, 0.1 Myr control also reproduces its previous output exactly.

An additional helium normalization diagnostic found a pre-existing limitation:
the evolved helium columns do not sum exactly to one. For LOS 5 the maximum
excess is 0.026986, identical in the earlier and current z6.15 calculations;
the z5.95 control also has a smaller excess (0.005562). The legacy outputs
are preserved rather than renormalized, and this diagnostic is recorded in
each quasar completion record. These arrays should not be described as a
validation or correction of the helium solver. The numerical comparisons
are in `checks/legacy_helium_closure.json`.

Full quasar-batch completion and errors are recorded in `qso_results.json`
and `results.json`; these records are the authority for the completed run.
