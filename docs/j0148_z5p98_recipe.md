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
make faster
```

This builds `rte1d_faster`. It caches the three photoionization cross-section
functions and reuses the fixed frequency grid and spectral weights. The cache
keys and cached values retain the exact floating-point inputs and results.
Representative J0148 sightlines ran about 1.5 times faster than `rte1d_fast`.
Tests at `tQ=0.01`, `1`, and `100` Myr, using both 500-sightline input blocks
and multiple ATON snapshots, produced byte-for-byte identical text outputs.
A known numerically failing sightline also failed identically in both builds.

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
  --snapshots 6.15 \
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

Here `N <= 1000`. As in the archived FIRE analysis, the solver can reject a
small number of numerically failed sightlines. The assembler excludes those
rows and writes their global LOS ids to a matching `deleted_indices_*.txt`
file beside the array. It still stops if it finds no valid sightlines.
