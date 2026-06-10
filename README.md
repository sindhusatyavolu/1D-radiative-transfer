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

At present, redshift and other defaults are set in `src/params.c`.

## Run

```bash
./rte1d tq_Myr t_on_Myr fduty LOS_index
```

Example:

```bash
./rte1d 10 1 1 0
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
