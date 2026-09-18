# Archived ATON z6.15 reproduction

## Finding

`qso_z6.15_m-27.09_tq1e+00_los.npy` uses the **z6.148904 input snapshot**
but the **z6.604205 RT and spectrum redshift**. Its shape is `(924, 891, 7)`.
The redshift in its filename does not describe the redshift used everywhere
in the calculation. It can identify the input snapshot/neutral-fraction
model rather than the quasar redshift. The follow-up
[paper mapping audit](aton_paper_mapping.md) nevertheless finds that
6.604205 does not match the paper's quasar identified by M1450 = -27.09,
even under a nearest-snapshot approximation. A stale compiled redshift or
an interim common-reference grid are possible explanations; the exact
source and historical intent have not been identified.

This explains why a run at the expected redshift does not match, while the
z5.95 control reproduces much more closely. Reproducing this historical
configuration and running a consistent z6.15 model are different objectives.
The normal solver defaults have not been changed to encode this archive's
redshift mismatch.

The archive is in:

```text
/Users/ssatyavolu/projects/proximity-p1d/proximity_zone_p1d/data/simulations/ATON/
```

## Evidence

1. All nine regenerated `example_input_files/z6p15_f1/f1kz6.15M3LOS_*.txt`
   grids were compared with the files in `fireqso_dominika`. The temperature
   and velocity arrays are identical; density differences are below 6e-16
   relative, and background differences below 4.4e-12 relative. The input
   regeneration does not explain the output mismatch.
2. Recomputing optical depths from the reference's *own final* temperature
   and xHI, plus the original input rho and peculiar velocities, fits
   z = 6.604205. This isolates the spectral calculation from the chemistry
   solver. For the z5.95 control the same procedure fits z = 5.945518.
3. At z = 6.604205, all 924 retained input mappings pass an optical-depth
   check at 29 bins per sightline, for each lifetime 0.01, 0.1, 1, 10, and
   100 Myr. The largest relative discrepancy is 5.83e-8. This validates the
   redshift and the existing inferred mapping on the sampled bins; it is
   not a full RT rerun of all 4,620 rows. See
   `outputs/diagnose_z6p15/all_rows_spectral_check.json`.
4. Full 1 Myr RT runs at z = 6.604205 reproduce representative sightlines
   from both f1k and l1k inputs. For input LOS 0, the mean absolute flux
   error falls from 0.015061 to 1.58e-8 when only the RT redshift is changed.
5. The archived monolithic runner reads `float Mag = atof(argv[5])`, whereas
   the current runner stores M1450 as a double. Float32(-27.09) is
   -27.09000015258789. Matching this precision reduces LOS 0's mean absolute
   flux error to 1.02e-11 (maximum 1.88e-9); maximum absolute xHI error is
   1.01e-13 and temperature error is 6.01e-7 K. The result is numerically
   very close but not bit-for-bit identical.

With both the redshift and magnitude precision matched, the four verified
full RT runs give:

| Input LOS | Reference row | Mean absolute flux error | Maximum absolute flux error |
| --- | --- | --- | --- |
| 0 | 0 | 1.02e-11 | 1.88e-9 |
| 5 | 5 | 1.34e-11 | 3.31e-9 |
| 7 | 6 | 2.79e-12 | 1.53e-10 |
| 999 | 923 | 2.42e-12 | 9.55e-11 |

Across these runs, maximum absolute xHI error is 1.01e-13 and maximum
relative temperature error is 3.03e-9 (an absolute difference of 0.000217 K
in a roughly 71,578 K cell). All seven columns' comparison metrics are
saved in `outputs/diagnose_z6p15/verified_rt_samples.json`.

The detailed LOS 0 comparison is in
`outputs/diagnose_z6p15/los0_report.json`; the figure is
`outputs/diagnose_z6p15/los0_comparison.png`. The z5.95 control is in
`outputs/diagnose_z6p15/z5p95_control.json`.

## Settings

| Setting | Archive reproduction value |
| --- | --- |
| Input file | `dom_zloc_M3_spec_z6.149.dat` |
| Redshift for input density conversion | **6.148904** |
| RT redshift (`--z-sim` in the current batch driver) | **6.604205** |
| Output label (`--z-qso`, with legacy density enabled) | 6.15 |
| Input prefixes | `f1kz6.15M3LOS`, `l1kz6.15M3LOS` |
| M1450 | -27.09 with `--legacy-float-magnitude` |
| Lifetime, on-time, duty cycle | 1 Myr, 1 Myr, 1 |
| Number of evolved cells | 891 |
| Columns per input file | 500 |
| Legacy density handling | 1; no division by X or redshift rescaling |
| Redward spectral pixels | 256 |
| Legacy helium recombination condition | 1 |
| Spectral boxsize override | 102441.268 |

Do **not** regenerate the input densities at 6.604205: the mismatch belongs
to the historical RT calculation. If using `generate_inputs_from_dom.py`,
pass `--z-sim 6.148904` explicitly; its default is 5.945518, and changing
the filename/prefix alone does not change that default.

The output columns containing physical fields and the optical-depth column
also have different coordinate origins in the archive: the optical-depth
calculation includes a 256-pixel redward offset. The saved array retains
only its first 891 optical-depth bins. The comparison plot accounts for
this offset on its flux panel.

## Reproduce selected sightlines

Run from the repository root with a Python environment containing NumPy.
On this machine, `/Users/ssatyavolu/miniconda3/bin/python` has the required
packages; the system `python3` does not have NumPy.

```bash
make optimized

python tools/run_durovcikova_grid.py \
  --rte ./rte1d_opt \
  --input-dir /Users/ssatyavolu/projects/fireqso_dominika \
  --input-prefix f1kz6.15M3LOS \
  --second-input-prefix l1kz6.15M3LOS \
  --output-dir outputs/reproduced_z6p15 \
  --z-sim 6.604205 --z-qso 6.15 \
  --m1450 -27.09 --legacy-float-magnitude \
  --tq-grid 1 --t-on 1 --fduty 1 \
  --n-cells 891 --n-los 500 \
  --legacy-density 1 --num-pixels-redward 256 \
  --legacy-alphahe2-bug 1 --boxsize-override 102441.268 \
  --los 0 5 7 999
```

These input LOS IDs correspond to reference rows 0, 5, 6, and 923.
Reference row numbers stop corresponding directly to input LOS IDs after
discarded sightlines; f1k/l1k input IDs are global 0–999.

The verification runs are saved under
`outputs/diagnose_z6p15/redshift6604_floatmag/`.

## Assemble all 924 rows

To run the full set, replace `--los 0 5 7 999` in the command above with:

```bash
  --los-file outputs/recover_z6p15/inferred_kept_indices_z6.15_m27.09_from_target.txt
```

After the runs finish, assemble them in reference order:

```bash
python tools/assemble_outputs_to_npy.py \
  --input-dir outputs/reproduced_z6p15 \
  --indices-file outputs/recover_z6p15/inferred_kept_indices_z6.15_m27.09_from_target.txt \
  --output outputs/reproduced_z6p15/qso_z6.15_m-27.09_tq1e+00_los.npy \
  --z-qso 6.15 --m1450 -27.09 --tq 1 --t-on 1 --fduty 1 \
  --n-cells 891 --output-los-mode input
```

Only representative full RT runs were performed during this diagnosis;
the complete 924-row `.npy` has not been regenerated. The retained list
is inferred from the archive, not recovered from original run metadata.
Some original input sightlines are duplicates, so identical field data
cannot uniquely establish which duplicate supplied a reference row.

## Repeat the independent redshift diagnostic

This helper requires NumPy and SciPy. It ports the existing C spectral
formula and holds all target fields fixed while varying only redshift.

```bash
python tools/diagnose_aton_reproduction.py \
  /Users/ssatyavolu/projects/proximity-p1d/proximity_zone_p1d/data/simulations/ATON/qso_z6.15_m-27.09_tq1e+00_los.npy \
  --input-dir /Users/ssatyavolu/projects/fireqso_dominika \
  --input-prefix f1kz6.15M3LOS --input-los 0 --target-row 0 \
  --redshifts 6.148904 6.604205 \
  --compare reproduced outputs/diagnose_z6p15/redshift6604_floatmag/z6.15LOS0tq1mag-27.09_tep_1Myr1fd.txt \
  --report outputs/diagnose_z6p15/recheck.json
```

Add `--plot PATH.png` to save a comparison figure (requires Matplotlib).
