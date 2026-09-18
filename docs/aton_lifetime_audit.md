# Legacy configuration and lifetime audit

The follow-up [paper mapping audit](aton_paper_mapping.md) distinguishes
the input snapshot label from the quasar redshift. Their difference alone
is not a mistake under the paper's method. However, the inferred 6.604205
still does not match the paper's -27.09 quasar (z6.1225; nearest snapshot
6.148904). The historical scenario below is one possibility, not proof of
an accidental configuration. The z5.95 control validates its original
same-redshift configuration, not redshift matching to that quasar.

## How the mismatch could occur

The legacy code did not derive its runtime redshift from the input data or
the output filename. Three independent pieces of configuration had to be
kept consistent:

1. **Compiled redshift.** In the surviving
   `/Users/ssatyavolu/projects/fireqso_dominika/dominika/magcode.c`, line 34
   reads `#define z 6.604205 //5.945518 //end redshift`. This constant controls
   the scale factor used by the RT calculation and the spectrum calculation.
   The `ascale` routine returns `1/(1+z)` independently of its time argument.
2. **Input filenames.** The same file has nine separately hard-coded
   `fopen` paths starting at line 1345. There is no check that the redshift
   encoded in these filenames agrees with `z`, and the text grids do not
   supply a redshift header that the reader verifies.
3. **Output labels.** The C text-output filename uses LOS, lifetime,
   magnitude, on-time and duty cycle. The surviving `save_spectra.py`,
   line 84, independently hard-codes the redshift label in `np.save(...)`.
   Renaming the NumPy output does not change the calculation that generated it.

The batch command in `dominika/job.sh`, line 18, supplies lifetime, LOS,
on-time, duty cycle and magnitude, but **no redshift argument**. The script
runs an already compiled `./mrte` without a compilation step.

A plausible sequence is therefore to adapt a z6.6 runner to read z6.15
inputs, leaving the compiled redshift at 6.604205, and later save the arrays
under a z6.15 label. Correcting the source's redshift but running an older
executable compiled with those z6.15 input paths would have the same effect.
All lifetime jobs using that executable would inherit the same mismatch.

This is a reconstruction of **how it could happen**, not a recovered record
of the specific historical edit. The currently surviving
`z6.14/magcode.c` has the correct `#define z 6.148904`, and the surviving
`dominika/magcode.c` uses z6.6 input paths. Neither has been identified as
the exact source used for the affected June 2023 archive. The configuration
actually present in the archive is established numerically by the checks
below.

## Validation method and scope

The audit concerns the ten ATON arrays named
`qso_z{6.15,5.95}_m-27.09_tq{1e-02,1e-01,1e+00,1e+01,1e+02}_los.npy`
in `proximity-p1d/proximity_zone_p1d/data/simulations/ATON`.

For each file, the diagnostic holds its final xHI and temperature fixed,
uses the original input rho and peculiar velocities, and independently
fits the spectral redshift on row 0. It then checks the full 891-bin
optical-depth profile of every retained row at the matching redshift.
Thus the spectral test covers all retained rows, not a sample of bins.

The z6.15 retained-index mapping is the previously recovered list. The
z5.95 mapping was inferred from the 0.01 Myr archive using the far-field
temperature and xHI profiles, with the existing monotonic alignment helper.
Some input skewers are duplicates, so the numerical mapping validates the
input field content without uniquely identifying the original duplicate ID.

Full RT runs were also performed for z6.15 input LOS 0 at 0.01, 0.1, 10,
and 100 Myr, using runtime z6.604205 and archived float32 magnitude
precision. The previously verified 1 Myr run supplies the fifth comparison.
These check the evolved fields as well as the spectral redshift.

## Results

| Lifetime (Myr) | z6.15 file: fitted spectral redshift | z5.95 file: fitted spectral redshift |
| --- | --- | --- |
| 0.01 | 6.604205 | 5.945518 |
| 0.1 | 6.604205 | 5.945518 |
| 1 | 6.604205 | 5.945518 |
| 10 | 6.604205 | 5.945518 |
| 100 | 6.604205 | 5.945518 |

Every z6.15 file has 924 retained rows and every z5.95 file has 730.
All 891 saved spectral bins were checked in every row: **7,368,570 bins**
across the ten files. All five z6.15 archives contain the mismatch; none
of the five z5.95 archives shows this particular redshift mismatch.

At z6.604205 the largest relative optical-depth residual over all the
z6.15 files is 1.86e-7. At z5.945518 the largest residual over the z5.95
files is 1.06e-6 (one row of the 100 Myr file); these checks are numerical
reconstructions from saved fields, not bitwise equality tests.
That largest z5.95 relative residual occurs at a very small optical depth,
tau = 3.34e-6; the maximum absolute flux residual on that entire row is
only 4.08e-11. An independent fit of that row gives z = 5.945518025,
also consistent with the intended redshift. Its details are saved in
`z5p95_worst_residual.json` alongside the other audit results.

The independent complete RT reruns for z6.15 LOS 0 give:

| Lifetime (Myr) | Mean absolute flux error against archive | Maximum absolute flux error |
| --- | --- | --- |
| 0.01 | 3.24e-12 | 8.29e-10 |
| 0.1 | 3.28e-12 | 5.19e-10 |
| 1 | 1.02e-11 | 1.87e-9 |
| 10 | 3.62e-12 | 3.24e-10 |
| 100 | 3.29e-11 | 3.76e-9 |

The maximum relative temperature discrepancy across these RT comparisons
is 1.48e-10. These additional runs confirm that the inferred z6.604205
configuration reproduces the evolved gas fields as well as the spectra
at each lifetime. They are one LOS per lifetime, alongside the all-row
spectral checks; a full RT rerun of every sightline was not performed.

The *presence* of the mismatch is common to all five lifetimes. The flux
impact percentages measured previously for 1 Myr should not be assumed
to hold at the other lifetimes without matched corrected-redshift runs.

## Saved results

Machine-readable results are saved under `outputs/lifetime_redshift_audit/`:

- `all_lifetimes_full_spectra.json`: fitted redshifts and all-row spectral errors.
- `z6p15_rt_all_lifetimes.json`: full RT versus archive, all seven columns.
- `inferred_kept_indices_z5p95.txt`: the inferred z5.95 row mapping.
- `rt_commands.json` and `z6p15_rt6604/`: full RT commands, outputs and logs.

To repeat the spectral audit:

```bash
python tools/audit_aton_lifetimes.py \
  --archive-dir /Users/ssatyavolu/projects/proximity-p1d/proximity_zone_p1d/data/simulations/ATON \
  --input-dir /Users/ssatyavolu/projects/fireqso_dominika \
  --kept-z615 outputs/recover_z6p15/inferred_kept_indices_z6.15_m27.09_from_target.txt \
  --kept-z595 outputs/lifetime_redshift_audit/inferred_kept_indices_z5p95.txt \
  --output outputs/lifetime_redshift_audit/all_lifetimes_full_spectra.json
```

This helper requires NumPy and SciPy, available in
`/Users/ssatyavolu/miniconda3/bin/python` on this machine.
