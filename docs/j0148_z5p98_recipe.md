# ULAS J0148+0600: Section 3.1 grid

This recipe generates the quasar-on Ly-alpha sightline grid for ULAS
J0148+0600 with

```text
z_qso = 5.98
M1450 = -27.39
```

It applies the ATON line-of-sight procedure described in
[Section 3.1 of Durovcikova et al. (2024)](https://arxiv.org/html/2401.10328v1#S3.SS1).
J0148 is below that paper's `z >= 6.0` sample cut, so this applies the paper's
method to J0148 rather than reproducing a saved J0148 grid from the paper.

The complete grid contains eight input neutral fractions and five quasar
lifetimes, giving 40 `.npy` files. This recipe does not generate no-quasar
files.

## What the driver calculates automatically

The eight ATON inputs are:

| Snapshot label | Exact snapshot redshift | Paper xHI |
| --- | ---: | ---: |
| `5.11` | 5.109794 | 0.000037 |
| `5.76` | 5.755954 | 0.07 |
| `5.95` | 5.945518 | 0.13 |
| `6.15` | 6.148904 | 0.21 |
| `6.60` | 6.604205 | 0.37 |
| `7.14` | 7.139393 | 0.54 |
| `8.15` | 8.150086 | 0.75 |
| `9.55` | 9.545964 | 0.89 |

For each selected snapshot, `tools/run_aton_section31_grid.py` reads its exact
redshift and calculates

```text
A = (1 + z_qso) / (1 + z_snapshot)
density multiplier = A^3
comoving-length multiplier = A^2
```

The runtime cosmology and quasar spectrum always use `z_qso=5.98`. The initial
`xHI`, temperature, peculiar velocity, density structure, and background rates
come from the selected ATON snapshot. Scaling `nH` and `nHI` together preserves
that snapshot's `xHI`.

The driver also corrects the known `0.76` normalization error in the existing
`fireqso_dominika` `nH` and `nHI` files. All of these transformations happen
inside the driver; do not supply density or length multipliers on the command
line and do not edit the input files.

The five lifetimes are `0.01`, `0.1`, `1`, `10`, and `100` Myr. For one
snapshot, the driver finishes all sightlines at one lifetime before starting
the next lifetime. The `--workers` option parallelizes sightlines only within
the current lifetime.

## 1. Build the faster solver

Run:

```bash
cd /Users/ssatyavolu/projects/proximity/1D-radiative-transfer
make -B faster
```

This builds `rte1d_faster`. It caches the three photoionization cross-section
functions and reuses the fixed frequency grid and spectral weights. It also
keeps the helium fractions normalized and retries a chemistry update with a
smaller local timestep if the proposed temperature or abundance is
nonphysical. This is required by the strongly rescaled `z=9.55` model. Solver
errors now return a failure status, and the Python runner stops if a solver
does not create a complete 891-row output.

Use `-B` here because an older `rte1d_faster` executable may already exist.
Do not rebuild this executable while another grid command is still launching
new sightlines from it; wait for that command to finish first.

The input directory must contain the `f1k` and `l1k` files for every snapshot
listed above.

## 2. Check available disk space

The full grid produces 40,000 sightline text files. At the measured size of
about 113 KB per sightline, these occupy approximately 4.2 GiB. The 40 final
NumPy arrays add approximately 1.9 GiB, for a total of about 6.1 GiB. Logs add
only a small amount.

Keep at least 10 GiB free before running the grid so macOS and the assembly
step have working space. Check with:

```bash
df -h /Users/ssatyavolu/projects/proximity/1D-radiative-transfer
```

Do not delete the text outputs until their NumPy arrays have been assembled
and verified.

## 3. Run two snapshots at a time on the M3 Mac

Use one Screen session for one snapshot. Start a session, for example:

```bash
screen -S j0148_511
```

Inside that session, run the following. Change only the `SNAPSHOT` value for
each session:

```bash
cd /Users/ssatyavolu/projects/proximity/1D-radiative-transfer

SNAPSHOT=5.11

/Users/ssatyavolu/miniconda3/bin/python tools/run_aton_section31_grid.py run \
  --input-dir /Users/ssatyavolu/projects/fireqso_dominika \
  --output-root outputs/j0148_z5p98_m2739/section31_grid \
  --z-qso 5.98 --m1450 -27.39 \
  --rte ./rte1d_faster \
  --snapshots "$SNAPSHOT" \
  --workers 3 --resume
```

Detach from the running session with `Ctrl-A`, then `D`.

For the first batch, the complete z5.11 command is:

```bash
cd /Users/ssatyavolu/projects/proximity/1D-radiative-transfer

/Users/ssatyavolu/miniconda3/bin/python tools/run_aton_section31_grid.py run \
  --rte ./rte1d_faster \
  --input-dir /Users/ssatyavolu/projects/fireqso_dominika \
  --output-root outputs/j0148_z5p98_m2739/section31_grid \
  --z-qso 5.98 --m1450 -27.39 \
  --snapshots 5.11 \
  --workers 3 --resume
```

In the second Screen session, the complete z5.76 command is:

```bash
cd /Users/ssatyavolu/projects/proximity/1D-radiative-transfer

/Users/ssatyavolu/miniconda3/bin/python tools/run_aton_section31_grid.py run \
  --rte ./rte1d_faster \
  --input-dir /Users/ssatyavolu/projects/fireqso_dominika \
  --output-root outputs/j0148_z5p98_m2739/section31_grid \
  --z-qso 5.98 --m1450 -27.39 \
  --snapshots 5.76 \
  --workers 3 --resume
```

There are eight snapshot jobs in total. On the 8-core M3 MacBook Air with
16 GB RAM, run two Screen sessions with `--workers 3` in each session. This
gives six concurrent solver processes. Run the snapshots in four batches:

| Batch | First Screen | Second Screen |
| ---: | ---: | ---: |
| 1 | `5.11` | `5.76` |
| 2 | `5.95` | `6.15` |
| 3 | `6.60` | `7.14` |
| 4 | `8.15` | `9.55` |

At the measured faster-solver rate, the complete grid should take approximately
24--32 hours if the Mac remains awake and the runs continue without long
interruptions. Sustained thermal throttling can increase this estimate.

Each snapshot command runs its five lifetimes sequentially in this order:

```text
0.01 -> 0.1 -> 1 -> 10 -> 100 Myr
```

`--resume` validates existing text outputs and skips complete sightlines. If a
session is interrupted, rerun the same command.

### Existing `z=9.55` results made with the older executable

The older solver completed all 1,000 sightlines only for `tQ=0.01` Myr. At
`tQ=0.1`, `1`, `10`, and `100` Myr it exited during the helium update and made
no sightline files, although its zero exit status made the batch look
successful. Do not assemble that directory.

After all currently running jobs have stopped, rebuild the solver, archive the
old `z=9.55` directory, and rerun all five lifetimes:

```bash
cd /Users/ssatyavolu/projects/proximity/1D-radiative-transfer
make -B faster

mv outputs/j0148_z5p98_m2739/section31_grid/qso/xhi_z9.55 \
   outputs/j0148_z5p98_m2739/section31_grid/qso/xhi_z9.55_before_chemistry_fix

/Users/ssatyavolu/miniconda3/bin/python tools/run_aton_section31_grid.py run \
  --rte ./rte1d_faster \
  --input-dir /Users/ssatyavolu/projects/fireqso_dominika \
  --output-root outputs/j0148_z5p98_m2739/section31_grid \
  --z-qso 5.98 --m1450 -27.39 \
  --snapshots 9.55 \
  --workers 6 --resume
```

For a single `z=9.55` job on the 16-core/32-thread ThinkStation, use
`--workers 16` if at least 20 GiB of memory is free.

That command is sufficient to replace the missing `z=9.55` portion of the
grid. It leaves the seven completed snapshot directories as legacy-chemistry
results. The correction changes Ly-alpha optical depth very little in most
pixels, but representative single-sightline tests found local differences up
to about 1% at `z=5.95` and 3% at `z=7.14`. For one numerically uniform,
scientifically corrected grid, archive the complete `qso/` directory and
rerun all eight snapshots with the rebuilt executable. The Section 3.1
density, length, quasar-redshift, and magnitude settings are unchanged; this
choice concerns only the chemistry integrator.

To inspect the Screen sessions, use:

```bash
screen -ls
screen -r j0148_511
```

## 4. Assemble the 40 NumPy arrays

After all eight snapshot jobs finish, run one assembly command outside Screen:

```bash
cd /Users/ssatyavolu/projects/proximity/1D-radiative-transfer

/Users/ssatyavolu/miniconda3/bin/python tools/run_aton_section31_grid.py assemble \
  --output-root outputs/j0148_z5p98_m2739/section31_grid \
  --z-qso 5.98 --m1450 -27.39
```

The arrays are written under `section31_grid/npy/`. Their names include both
the neutral-fraction snapshot and the actual quasar redshift, for example:

```text
qso_z6.15_zqso5.98_m-27.39_tq1e+00_los.npy
```

Each array has shape `(N, 891, 7)` with columns

```text
distance_ckpc, xHI, temperature_K, xHeI, xHeII, xHeIII, tau_Lya
```

For this grid, `N` should be 1,000. The assembler records any absent global LOS
ids in a matching `deleted_indices_*.txt` file, but a nonempty list now means
the run is incomplete and should be resumed before using the array.
