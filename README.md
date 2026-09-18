# 1D Radiative Transfer

Minimal modular C version of the 1D quasar proximity-zone radiative-transfer code.

## Build

Requires a C compiler and `make`.

```bash
make
```

This creates the executable:

```bash
./rte1d
```

To remove the executable:

```bash
make clean
```

## Inputs

Run the code from the directory containing the input sightline files. The current file naming convention is inherited from the original code and expects files like:

```text
f1kz5.95M3LOS_nh.txt
f1kz5.95M3LOS_nh1.txt
f1kz5.95M3LOS_T.txt
f1kz5.95M3LOS_vpec.txt
f1kz5.95M3LOS_rho.txt
f1kz5.95M3LOS_Gammabg.txt
f1kz5.95M3LOS_Gammabghe1.txt
f1kz5.95M3LOS_Gammabghe2.txt
f1kz5.95M3LOS_Heatbg.txt
```

By default the executable reads from `example_input_files` using the prefix
derived from `z` in `src/params.c`, for example `f1kz5.95M3LOS`.

## Run

```bash
./rte1d tq_Myr t_on_Myr fduty LOS_index
```

Example:

```bash
./rte1d 10 1 1 0
```

Optional arguments allow one-off runs with a different quasar/snapshot:

```bash
./rte1d tq_Myr t_on_Myr fduty LOS_index M1450 z_sim z_qso input_dir input_prefix output_dir
```

Example:

```bash
./rte1d 10 10 1 0 -27.0 5.945518 6.10 example_input_files f1kz5.95M3LOS outputs
```

Arguments:

- `tq_Myr`: quasar lifetime in Myr
- `t_on_Myr`: quasar on-time in Myr
- `fduty`: duty cycle, must be in `(0, 1]`
- `LOS_index`: sightline index, currently `0` to `499`

## Output

The code writes one text file in the run directory. Each row contains:

```text
distance_ckpc xHI temperature_K xHeI xHeII xHeIII tau_Lya
```

The output filename follows the original convention, for example:

```text
z5.95LOS0tq10mag-27.00_tep_1Myr1fd.txt
```

## Regenerate Inputs

The helper below regenerates the direct skewer files and the equilibrium
background files from a `dom_zloc_*.dat` file:

```bash
python tools/generate_inputs_from_dom.py \
  example_input_files/dom_zloc_M3_spec_z5.946.dat \
  --output-dir example_input_files/regenerated_direct \
  --prefix f1kz5.95M3LOS \
  --with-backgrounds \
  --compare \
  --write-halo-catalog
```

The direct files are `nh`, `nh1`, `T`, `vpec`, and `rho`. The generated
background files are `Gammabg`, `Gammabghe1`, `Gammabghe2`, and `Heatbg`.

## Batch Runs

Use the Python driver to loop over lifetimes and sightlines while keeping the C
code as a single-run executable:

```bash
python tools/run_durovcikova_grid.py \
  --input-dir example_input_files/regenerated_direct \
  --input-prefix f1kz5.95M3LOS \
  --output-dir outputs/z6p10_xHI0p21 \
  --z-sim 5.945518 \
  --z-qso 6.10 \
  --m1450 -27.0 \
  --tq-grid 0.01 0.1 1 10 \
  --los 0:499
```

Use `--dry-run` to print the commands without running them.

## Reproducing the archived ATON z6.15 outputs

The [ULAS J0148+0600 recipe](docs/j0148_z5p98_recipe.md) gives complete
commands for applying the Section 3.1 snapshot rescaling to a quasar at
z=5.98 and M1450=-27.39, including all five lifetimes, 1000 sightlines,
matched pre-quasar spectra, and NumPy assembly.

The [corrected-redshift FIRE pair](docs/corrected_fire_z6p15_pair.md) generates
the z6.15, M1450=-27.09, 1 Myr array with runtime cosmology at the actual
snapshot redshift 6.148904, plus a matched pre-quasar array. It preserves
the original 924-row selection and the legacy density conventions used
in the z5.95 reproduction. Run `tools/generate_fire_snapshot_pair.py` after
building the optimized solver; both outputs include configuration and
validation records.

The [newly downloaded -26.99, z5.76, 1 Myr benchmark](docs/downloaded_m26p99_benchmark.md)
also reproduces with the current solver: four complete RT sightlines and
spectral reconstruction of all 945 rows. It uses original z5.755954 densities
and runtime z6.604205 without density rescaling. The notes give exact commands,
errors, artifacts, and the remaining uncertainty about final-paper provenance.

The [literal Section 3.1 comparison](docs/literal_section31_comparison.md)
reruns four matched sightlines with zsim=5.76, zQSO=6.14, M1450=-26.99,
and tq=1 Myr using the printed density and physical-length multipliers.
It includes both background treatments, spectra on a common velocity grid,
600 km/s flux averages, and reproducible commands. The optional
`--comoving-length-scale` driver argument scales cell spacing, inner radius,
and spectral box length together; it defaults to 1 and does not scale density.

The [m26p58_rescaledrho reproduction](docs/rescaledrho_m26p58_reproduction.md)
checks a different z5.76, 1 Myr file. Four full RT sightlines match when RT
densities are rescaled to the observed quasar redshift 7.0848 and the legacy
input X factor is removed, while runtime geometry and spectral synthesis
remain at the snapshot redshift 5.755954.

The [taualpha_wqso lightcone audit](docs/lightcone_taualpha_wqso_reproduction.md)
checks the separate mean-free-path lightcones, labelled z6.14, M1450=-27.00,
and tq=1 Myr, against the requested z6.15, M1450=-27.09 parameters. It reconstructs
the stitched inputs and uses a separate build with `RTE_MAX_CELLS=8192` for the
8071-cell sightlines; the ordinary build retains its 2048-cell capacity.

The [paper mapping audit](docs/aton_paper_mapping.md) checks whether the
filename redshift instead labels a neutral-fraction snapshot for a quasar
at another redshift. That interpretation does not resolve the -27.09
archive discrepancy: the paper's corresponding quasar is at z6.1225
(nearest snapshot 6.148904), while most local files synthesize spectra at
6.604205. The z5.95 reproduction is a baseline check, not validation of
the paper's cross-redshift rescaling.

The archived `qso_z6.15_m-27.09_tq1e+00_los.npy` combines inputs from
`dom_zloc_M3_spec_z6.149.dat` (density conversion at **z = 6.148904**) with
an RT/spectrum calculation at **z = 6.604205**. Using 6.148904 for the RT run
therefore does not reproduce that archive. The z5.95 reference uses 5.945518
for both steps.

See [the reproduction notes](docs/aton_z6p15_reproduction.md) for the measured
evidence, complete run/assembly commands, and the legacy settings required.
The batch driver's `--legacy-float-magnitude` option also reproduces the old
runner's float32 magnitude rounding without changing the normal C defaults.

The [Ly-alpha impact comparison](docs/lya_redshift_impact.md) measures the
effect of correcting only this redshift in 24 matched skewers, including
mean transmission profiles, individual flux differences, and proximity-zone
boundaries.

The [legacy configuration and lifetime audit](docs/aton_lifetime_audit.md)
explains how independent hard-coded settings could create the mismatch and
checks all saved spectral bins for both redshifts at all five lifetimes.
