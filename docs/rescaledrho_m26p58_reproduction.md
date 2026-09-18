# Reproduction of one m26p58_rescaledrho file

The selected file **does contain RT density rescaling to the observed quasar
redshift**, together with removal of the legacy input factor X=0.76. Four complete
sightlines reproduce numerically with these densities while retaining the
**snapshot redshift for runtime geometry and Lyα synthesis**.

## File and quasar

```text
/Users/ssatyavolu/projects/fireqso_dominika/updated_aug23/m26p58_rescaledrho/qso_z5.76_m-26.58_tq1e+00_los.npy
shape: (994, 2048, 7)
SHA256: 4956349d8965c9e4cf210cbb921d42827e5f53bdd7ad5e8ce411e26f2982d22e
```

The M1450=-26.58 entry is **J1120+0641, z=7.0848**, in
[Table 1 of Ďurovčíková et al. (2024)](https://arxiv.org/html/2401.10328v1#S2.T1).
The input snapshot is `dom_zloc_M3_spec_z5.756.dat`, whose physical input-density
conversion uses **zsim=5.755954**. The original f1k/l1k text inputs were checked
against the binary file for all four selected sightlines, including background
rates, before the reruns. The archive's own deleted-index list supplies the
retained-row mapping.

## Recovered configuration

Both nH and nHI passed to RT are multiplied by

```text
(1 / 0.76) * [(1 + 7.0848) / (1 + 5.755954)]^3
= 1.3157894736842106 * 1.7137512755321869
= 2.2549358888581406
```

The factor **1.7137513** is the redshift density rescaling. The separate **1/0.76**
undoes the factor introduced in the original input generator, as documented in
the opening comment of `fireqso_dominika/get_los.py`. Their product is the factor
on the existing text files; it should not be interpreted as a 2.255-fold
cosmological rescaling alone.

| Component | Value that reproduces this file |
| --- | --- |
| Input snapshot | 5.755954 |
| Density target redshift | **7.0848**, the observed quasar redshift |
| Initial nH, nHI | Original text inputs × **2.2549358888581406** |
| Initial xHI, T, peculiar velocity | Original snapshot values |
| Background ionization/heating | Original snapshot files |
| Runtime redshift | **5.755954** |
| Comoving cell size / inner radius | Original 114.88970588235293 |
| Spectral boxsize | Original 235294.118 |
| Comoving rho supplied to spectra | Original input, multiplier **1** |
| Spectral scale factor and Hubble conversion | Evaluated at **5.755954** |
| M1450 | -26.58 represented as legacy float32 |
| Lifetime, on-time, duty cycle | 1 Myr, 1 Myr, 1 |
| Cells / redward spectral offset | 2048 / 256 |
| Legacy helium coefficient option | Enabled, as in the other reproductions |

The matched runtime redshift is **not 6.604205**, which reproduced the earlier
downloaded -26.99 file. These are distinct numerical configurations.

This behavior is consistent with the current code's separation of variables:
`main.c` uses `z_qso` in the density multiplier, whereas `ascale()` in
`src/utils.c` and spectral synthesis in `src/lya_spectra.c` use `z`. Updating
the density target alone does not update those geometry/spectral conversions.
For the explicit audit, inputs were pre-multiplied and the additional C density
scaling disabled to prevent applying the factors twice.

## Evidence

| Reference row | Global input LOS | Maximum absolute flux difference | Maximum absolute T difference (K) |
| --- | --- | ---: | ---: |
| 0 | 0 | 1.26e-9 | 6.33e-6 |
| 331 | 334 | 1.08e-10 | 3.88e-6 |
| 662 | 668 | 7.26e-10 | 1.99e-6 |
| 993 | 999 | 1.49e-9 | 2.25e-6 |

The maximum absolute xHI difference across the four complete RT runs is
**1.00e-13**. Independently reconstructing all 2048 spectral bins from each
rerun's final fields agrees with the C spectra to **1.50e-9** in flux.

Row 0 controls distinguish the relevant density choices:

| Alternative | Maximum absolute flux difference from reference |
| --- | ---: |
| Original densities without either factor | 0.39978 |
| Redshift factor to 7.0848, retaining the input 0.76 | 0.17378 |
| Redshift factor to nearest snapshot 7.139393, removing 0.76 | 0.01888 |
| Redshift factor to observed 7.0848, removing 0.76 | **1.26e-9** |

Thus the observed redshift, rather than the nearest simulated redshift,
reproduces the density treatment of the tested sightlines.

![Matched profiles and residuals](../outputs/rescaledrho_m26p58_z5p76/reproduction_comparison.png)

## Comparison with Section 3.1

The redshift-dependent **RT density multiplier** agrees with the one printed in
[Section 3.1](https://arxiv.org/html/2401.10328v1#S3.SS1). The extra X correction
addresses an input-unit factor and is separate from that prescription.

The geometry and spectral conversion in the matched run remain at the input
snapshot redshift. A further row 0 control retains the matched RT densities,
changes runtime redshift to 7.0848, and applies the paper's literal physical
length multiplier `(1+7.0848)/(1+5.755954)`. To implement this in the C code's
comoving variables, cell size, inner radius, and boxsize are multiplied by the
square of that ratio, **1.432073123013158**. This control does not reproduce the
stored final fields (maximum relative T difference **0.3291**). The accompanying
figure compares the spectra on their respective velocity axes; it does not
assume corresponding pixel indices represent the same velocity.

![Literal length/runtime change on row 0](../outputs/rescaledrho_m26p58_z5p76/literal_section31_comparison.png)

This is a check of one file using four complete sightlines, plus row 0 controls.
It does not establish the configuration of every other file in the folder or
which files were used in the final published inference.

## Reproduce and inspect

```bash
make optimized
/Users/ssatyavolu/miniconda3/bin/python tools/audit_rescaledrho_reproduction.py \
  --cases density_observed_Xcorrected --rows 0 331 662 993
```

For the controls and figures:

```bash
/Users/ssatyavolu/miniconda3/bin/python tools/audit_rescaledrho_reproduction.py \
  --cases unscaled density_nearest density_nearest_Xcorrected density_observed literal_observed_Xcorrected --rows 0
/Users/ssatyavolu/miniconda3/bin/python tools/summarize_rescaledrho_reproduction.py
```

All artifacts are under `outputs/rescaledrho_m26p58_z5p76/`:

- `benchmark.json`: exact recovered configuration, four complete comparisons,
  controls, and independent spectral validation.
- `reproduced_sample_los.npy`: four reproduced rows, shape `(4, 2048, 7)`.
- `matched_samples.npz`: reference and reproduced arrays, explicit row mapping,
  literal-control output, and native velocity coordinates.
- `reference.json`, `kept_indices.txt`, `input_verification_los*.json`: provenance
  and binary/text input verification.
- `manifest_*.json`, `results_*.json`, and per-case directories: exact commands,
  timings, prepared inputs, full logs, output text, and individual comparisons.
- `reproduction_comparison.png/.pdf`, `literal_section31_comparison.png/.pdf`:
  inspected figures.
