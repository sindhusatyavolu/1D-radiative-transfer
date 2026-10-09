# ULAS J0148+0600: corrected ATON grid

This recipe generates the quasar-on grid for ULAS J0148+0600 using

```text
z_qso = 5.98
M1450 = -27.39
tQ = 0.01, 0.1, 1, 10, 100 Myr
```

The eight ATON snapshots supply the neutral-fraction and temperature fields:

| Label | Exact snapshot redshift | Paper xHI |
| --- | ---: | ---: |
| `5.11` | 5.109794 | 0.000037 |
| `5.76` | 5.755954 | 0.07 |
| `5.95` | 5.945518 | 0.13 |
| `6.15` | 6.148904 | 0.21 |
| `6.60` | 6.604205 | 0.37 |
| `7.14` | 7.139393 | 0.54 |
| `8.15` | 8.150086 | 0.75 |
| `9.55` | 9.545964 | 0.89 |

The complete run produces 40 NumPy files. It does not include no-quasar files.

## Scientific convention

The production defaults use one consistent coordinate system:

- The comoving ATON grid is unchanged for every snapshot.
- RT geometry, the Hubble velocity conversion, and Ly-alpha synthesis use
  `z_qso=5.98`.
- With `A=(1+z_qso)/(1+z_snapshot)`, the legacy text-file `nH` and `nHI` are
  multiplied by `A^3/0.76`. The division by `0.76` removes the documented
  preprocessing error in these particular input files.
- `xHI`, temperature, peculiar velocity, and the comoving mass-density field
  retain their snapshot values.
- The `rho` files store total gas mass density. Ly-alpha synthesis multiplies
  them by `X=0.76` when converting to hydrogen number density.
- The Ly-alpha calculation treats the stored spacing as ckpc and `rho` as the
  Gadget density unit proportional to `h^2`. It uses `H0=100h` for the Hubble
  velocity. This removes two inherited `h` factors that assumed ckpc/h even
  though the preprocessing had already divided the grid spacing by `h`.
- The cell width is calculated from the ATON resolution as
  `78.125/0.678 = 115.228613569 ckpc`, and the seconds-to-Myr conversion is
  used consistently for `H0` in the thermochemistry clock.
- Background photoionization rates are multiplied by the same total density
  factor, and background heating by its square. This keeps the supplied
  pre-quasar ionization and thermal state consistent after density rescaling.
- The Ly-alpha spectral box is exactly `891 * cellsize`; the corrected helium
  condition is used, and the integrator reaches the full requested lifetime.

The driver computes these factors. Do not edit the input files or add manual
density or length arguments. The detailed derivation and numerical checks are
in [section31_scientific_audit.md](section31_scientific_audit.md).

## Build

Use the same checked-out commit for the entire grid:

```bash
cd /path/to/1D-radiative-transfer
git rev-parse HEAD
make -B faster
```

Do not rebuild `rte1d_faster` while a grid is running. Every snapshot output
directory receives a `.run_configuration.json` containing the executable hash
and all physics settings. `--resume` refuses to combine incompatible outputs.

## Run on the ThinkStation

Use a new output root. The old `section31_grid` directories contain the
superseded length/background convention and must not be resumed into this run.

```bash
cd /path/to/1D-radiative-transfer

SNAPSHOT=5.11

python tools/run_aton_section31_grid.py run \
  --rte ./rte1d_faster \
  --input-dir /path/to/aton_sightlines \
  --output-root outputs/j0148_z5p98_m2739/section31_physical_v1 \
  --z-qso 5.98 --m1450 -27.39 \
  --snapshots "$SNAPSHOT" \
  --workers 4 --resume
```

Run one snapshot per Screen session. On the 16-core/32-thread ThinkStation,
four sessions with four workers each use the 16 physical cores:

```text
5.11  5.76  5.95  6.15
6.60  7.14  8.15  9.55
```

Start the second group as slots become free, or run four at once if memory and
cooling remain comfortable. Each command runs the five lifetimes sequentially.
If a session is interrupted, rerun exactly the same command with `--resume`.

For a small preflight before the full run, add `--los 0` and use a separate
output root such as `/tmp/j0148_preflight`.

## Assemble and verify

After all eight snapshot commands finish:

```bash
cd /path/to/1D-radiative-transfer

python tools/run_aton_section31_grid.py assemble \
  --output-root outputs/j0148_z5p98_m2739/section31_physical_v1 \
  --z-qso 5.98 --m1450 -27.39
```

Assembly now fails if any sightline is absent or invalid. Each successful file
has shape `(1000, 891, 7)` and a JSON metadata sidecar. A filename is, for
example,

```text
qso_z6.15_zqso5.98_m-27.39_tq1e+00_los.npy
```

The historical seven-column format is preserved:

```text
stored_grid_coordinate_comoving_kpc, xHI, temperature_K,
xHeI, xHeII, xHeIII, tau_Lya
```

Column 0 has the same 115.228613569 ckpc spacing in all 40 arrays. Columns 1--5
are gas-cell fields at that stored coordinate. Column 6 is a spectrum with 256
pixels redward of the quasar; its zero-Hubble-velocity origin is therefore
array index 256. Plot the Ly-alpha optical depth against proper Mpc using

```python
r_tau_pMpc = (data[:, 0] - data[256, 0]) / (1000 * (1 + z_qso))
```

The C solver uses `z_qso` for the same conversion to Hubble velocity. This
offset is part of the old file format and does not remove pixels. It is also
recorded in each JSON metadata sidecar. For this quasar the coordinate and
solver spacing is about 11.49 km/s per pixel; the exact value is included in
the sidecar.

Keep at least 10 GiB free: the text files and final arrays together require
about 6 GiB. Retain the text files until the NumPy files and metadata have been
checked and copied safely.
