# ULAS J0148+0600 at z=5.98 and M1450=-27.39

This recipe applies the ATON line-of-sight procedure described in
[Section 3.1 of Ďurovčíková et al. (2024)](https://arxiv.org/html/2401.10328v1#S3.SS1)
to ULAS J0148+0600. The quasar is not one of the
18 objects in that paper's Table 1 (the paper selects `z >= 6.0`), so this is
an application of its simulation method to J0148 rather than a reproduction
of one of that paper's saved J0148 model files.

## Snapshot and rescaling

Use the snapshot labelled `z5.95`, whose exact redshift is

```text
z_snapshot = 5.945518
```

and use the actual quasar redshift for the solver's cosmology and spectrum:

```text
z_runtime = z_qso = 5.98
M1450 = -27.39
```

Section 3.1 gives the physical-length and physical-density factors

```text
A = (1 + z_qso) / (1 + z_snapshot) = 1.0049646405062949
A^3 = 1.0149679868517030
```

Thus both physical `nH` and `nHI` are multiplied by `A^3`. Their ratio, the
initial `xHI`, is unchanged. Initial `xHI`, temperature, and peculiar velocity
come from the z=5.945518 snapshot.

The solver stores its cell spacing in comoving units and converts it to a
physical length using the runtime redshift. To obtain the paper's physical
factor `A`, pass

```text
comoving_length_scale = A^2 = 1.0099539286679464
```

because `(A^2 / (1+z_qso)) / (1/(1+z_snapshot)) = A`.

The original `fireqso_dominika` nH/nHI text files contain the known extra
factor `X=0.76`. With these specific legacy files, `--legacy-density 0`
removes that factor and `--density-multiplier A^3` performs the redshift
rescaling. The total nH/nHI multiplier relative to those files is
`A^3/0.76 = 1.3354841932259249`.

If using newly generated inputs that already omit the erroneous 0.76 factor,
use `--legacy-density 1` with the same `--density-multiplier A^3` instead.

## 1. Build the optimized solver

Run from this repository:

```bash
cd /Users/ssatyavolu/projects/proximity/1D-radiative-transfer
make optimized
```

The input directory must contain both 500-sightline sets and all nine fields:

```text
/Users/ssatyavolu/projects/fireqso_dominika/f1kz5.95M3LOS_{nh,nh1,T,vpec,rho,Gammabg,Gammabghe1,Gammabghe2,Heatbg}.txt
/Users/ssatyavolu/projects/fireqso_dominika/l1kz5.95M3LOS_{nh,nh1,T,vpec,rho,Gammabg,Gammabghe1,Gammabghe2,Heatbg}.txt
```

## 2. Inspect one command without running it

```bash
/Users/ssatyavolu/miniconda3/bin/python tools/run_durovcikova_grid.py \
  --rte ./rte1d_opt \
  --input-dir /Users/ssatyavolu/projects/fireqso_dominika \
  --input-prefix f1kz5.95M3LOS \
  --second-input-prefix l1kz5.95M3LOS \
  --output-dir outputs/j0148_z5p98_m2739/qso \
  --z-sim 5.98 --z-qso 5.98 --m1450 -27.39 \
  --tq-grid 1 --n-cells 891 --n-los 500 \
  --legacy-density 0 --density-multiplier 1.014967986851703 \
  --comoving-length-scale 1.0099539286679464 \
  --num-pixels-redward 256 --legacy-alphahe2-bug 1 \
  --boxsize-override 102441.268 --los 0 --dry-run
```

Here `--z-sim` is retained as the driver's historical option name; it is the
solver's runtime cosmology redshift. The source snapshot redshift is represented
by the selected input files and the explicit density/length multipliers.

## 3. Run every sightline and lifetime

The following runs all 1000 input sightlines for the five paper-grid lifetimes,
0.01, 0.1, 1, 10, and 100 Myr. With `t_on` omitted, each run uses
`t_on=tQ`; the duty cycle defaults to one. Six independent solver processes
run concurrently. Reduce `--workers` if memory or CPU pressure is too high.

```bash
/Users/ssatyavolu/miniconda3/bin/python tools/run_durovcikova_grid.py \
  --rte ./rte1d_opt \
  --input-dir /Users/ssatyavolu/projects/fireqso_dominika \
  --input-prefix f1kz5.95M3LOS \
  --second-input-prefix l1kz5.95M3LOS \
  --output-dir outputs/j0148_z5p98_m2739/qso \
  --z-sim 5.98 --z-qso 5.98 --m1450 -27.39 \
  --tq-grid 0.01 0.1 1 10 100 \
  --n-cells 891 --n-los 500 \
  --legacy-density 0 --density-multiplier 1.014967986851703 \
  --comoving-length-scale 1.0099539286679464 \
  --num-pixels-redward 256 --legacy-alphahe2-bug 1 \
  --boxsize-override 102441.268 \
  --los 0:999 --workers 6 --resume
```

Each run writes one text file. Global LOS 0–499 are read from `f1k`; LOS
500–999 are read from `l1k`, with 500 subtracted from the input column index.
The output filenames retain the global LOS IDs. `--resume` validates and skips
complete outputs, making the command safe to repeat after an interruption.

`--legacy-alphahe2-bug 1` retains the legacy helium condition used by the
archived calculations. Set it to zero only if deliberately creating a new,
corrected helium model that will not be compared numerically with that grid.

## 4. Generate the matched no-quasar sightlines

This writes the initial snapshot state and its Ly-alpha spectrum without
evolving the gas. It retains the same density and length transformations.

```bash
/Users/ssatyavolu/miniconda3/bin/python tools/run_durovcikova_grid.py \
  --rte ./rte1d_opt \
  --input-dir /Users/ssatyavolu/projects/fireqso_dominika \
  --input-prefix f1kz5.95M3LOS \
  --second-input-prefix l1kz5.95M3LOS \
  --output-dir outputs/j0148_z5p98_m2739/noqso \
  --z-sim 5.98 --z-qso 5.98 --m1450 -27.39 \
  --tq-grid 1 --n-cells 891 --n-los 500 \
  --legacy-density 0 --density-multiplier 1.014967986851703 \
  --comoving-length-scale 1.0099539286679464 \
  --num-pixels-redward 256 --legacy-alphahe2-bug 1 \
  --boxsize-override 102441.268 \
  --los 0:999 --workers 6 --resume --initial-only
```

No-quasar here means immediately before quasar turn-on. It retains the
snapshot's ATON ionization state; it is not a second one-Myr evolution with
only the ionizing background.

## 5. Assemble the text files into NumPy arrays

For example, assemble the 1 Myr quasar result:

```bash
/Users/ssatyavolu/miniconda3/bin/python tools/assemble_outputs_to_npy.py \
  --input-dir outputs/j0148_z5p98_m2739/qso \
  --los 0:999 \
  --output outputs/j0148_z5p98_m2739/qso_z5.98_m-27.39_tq1e+00_los.npy \
  --z-qso 5.98 --m1450 -27.39 --tq 1 --n-cells 891
```

Repeat with `--tq 0.01`, `0.1`, `10`, and `100`, changing the output filename
to `tq1e-02`, `tq1e-01`, `tq1e+01`, and `tq1e+02`, respectively.

Assemble the no-quasar result:

```bash
/Users/ssatyavolu/miniconda3/bin/python tools/assemble_outputs_to_npy.py \
  --input-dir outputs/j0148_z5p98_m2739/noqso \
  --los 0:999 \
  --output outputs/j0148_z5p98_m2739/noqso_z5.98_los.npy \
  --z-qso 5.98 --m1450 -27.39 --tq 0 --n-cells 891 --initial-only
```

Successful complete arrays have shape `(1000, 891, 7)`. The assembler stops
if any text output is absent or empty, rather than silently changing the
sightline selection. If a run failed, inspect its file under the `logs/`
subdirectory and repeat Step 3 or 4 with `--resume`.

The columns are

```text
distance_ckpc, xHI, temperature_K, xHeI, xHeII, xHeIII, tau_Lya
```

As in the archived FIRE files, the optical-depth coordinate includes a
256-pixel redward offset. For an analysis distance anchored at Ly-alpha,
use `(distance_ckpc - distance_ckpc[256]) / (1000 * (1 + 5.98))`.

## Using another neutral-fraction snapshot

For a snapshot at exact redshift `z_snapshot`, select its `f1k`/`l1k` input
prefixes and recompute

```text
A = (1 + 5.98) / (1 + z_snapshot)
density_multiplier = A^3
comoving_length_scale = A^2
```

The snapshot supplies initial xHI and temperature. Multiplying nH and nHI
together preserves xHI while moving the physical density normalization to
z=5.98. Keep `--z-sim 5.98 --z-qso 5.98`; do not replace the runtime redshift
with the snapshot redshift.
