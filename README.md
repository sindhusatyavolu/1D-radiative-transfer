# 1D Radiative Transfer

Used Codex to make the main C version of the 1D quasar proximity-zone radiative-transfer code modular and optimized. 

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


