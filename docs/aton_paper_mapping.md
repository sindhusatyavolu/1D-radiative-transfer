# ATON archive versus the paper's quasar redshift mapping

## Conclusion

**New reference, 2026-09-07:** The user supplied a separate downloaded
`qso_z5.76_m-26.99_tq1e+00_los.npy`, tentatively identified as an original-analysis
output. Its [dedicated benchmark](downloaded_m26p99_benchmark.md) reproduces four
complete RT sightlines and spectrally reconstructs all 945 rows. It uses original
z5.755954 densities with runtime z6.604205, without redshift density rescaling.
This extends the evidence beyond the earlier one-row -26.99 z5.11 spectral check;
it does not retrospectively validate that other file's gas evolution or prove
which file versions entered the published analysis.

**Scope clarification:** This is primarily an audit of the 41 files in
`proximity-p1d/proximity_zone_p1d/data/simulations/ATON`, with nine supplementary
spectral checks under `fireqso_dominika/updated_aug23`. It is not an audit of
all `fireqso_dominika` outputs or of the final published analysis. Mapping a
filename magnitude to a Table 1 object is a hypothesis about the file's
intended use, not recovered run metadata. The conclusions about disagreement
with the paper are conditional on these being the final quasar-matched
outputs, with no further relevant processing.

Interpreting the filename redshift as the input snapshot, and hence as an
index into the neutral-fraction grid, is consistent with the construction
in [Durovcikova et al. (2024), Section 3.1](https://arxiv.org/html/2401.10328v1#S3.SS1).
Different input and quasar redshifts are therefore not evidence of an error
by themselves. However, that interpretation does **not** explain the
specific `qso_z6.15_m-27.09_...` mismatch in the local proximity-p1d archive.

The paper's Table 1 associates M1450 = -27.09 with **J1509-1749 at z =
6.1225**. The closest locally available simulation redshift is **6.148904**,
not 6.604205. The same table associates -26.99 with J1319+0950 at 6.1347,
also closest to 6.148904. These associations are conditional on the filename
magnitudes identifying the paper's objects; the arrays contain no object ID.
The table entries agree in arXiv versions 1 and 2.

For J1509, the paper's construction would hold the quasar redshift at
6.1225 (or 6.148904 under the suggested nearest-snapshot approximation)
while varying the input snapshot. The checked local arrays do not do this.
This establishes a discrepancy in these copies, not which copies or later
processing were used to obtain the published results.

## Spectral audit of all 41 files in the local archive

Directory:

```text
/Users/ssatyavolu/projects/proximity-p1d/proximity_zone_p1d/data/simulations/ATON
```

All eight snapshot labels have five M1450 = -27.09 lifetimes: 0.01, 0.1,
1, 10, and 100 Myr. The fitted spectral redshift is the same at all five
lifetimes within each family:

| Filename snapshot label | Fitted spectral redshift | Files checked |
| --- | --- | --- |
| 5.11 | 6.604205 | 5 |
| 5.76 | 6.604205 | 5 |
| 5.95 | 5.945518 | 5 |
| 6.15 | 6.604205 | 5 |
| 6.60 | 6.604205 | 5 |
| 7.14 | 6.604205 | 5 |
| 8.15 | 6.604205 | 5 |
| 9.55 | 6.604205 | 5 |

The additional file `qso_z5.11_m-26.99_tq1e+02_los.npy` also fits 6.604205.
Only its first retained spectrum was checked; its RT density rescaling and
evolved-field history have not been reconstructed. The result therefore
does not establish that it shares the full chemistry configuration inferred
for the -27.09 z6.15 archive. No -26.99 files were found under
`fireqso_dominika`, including a search of hidden and ignored paths.

This new audit checks **the first retained sightline of each file**, fitting
redshift and an independent multiplicative optical-depth amplitude. It uses
the archive's own final xHI and temperature, plus the original input rho
and peculiar velocities. Every match is to input LOS 0. The fit uses every
seventh bin, then validates all 891 saved bins at the fitted parameters.
Multiple starting redshifts avoid local minima in the resolved spectra.

The largest relative optical-depth residual across these full row profiles
is **2.11e-8**. Every fitted amplitude is unity within **3.33e-10**. Testing
the actual quasar redshift or nearest snapshot, while allowing an arbitrary
amplitude, does not reproduce the profiles. For the 1 Myr z6.15 row, the RMS
log-optical-depth residual at z6.148904 remains 0.915 after optimizing that
amplitude.

The spectral calculation uses the independently established legacy boxsize
102441.268 and 256 redward pixels. These are fits under the recovered
legacy coordinate and spectral conventions; no claim of model-independent
redshift determination is made.

The earlier [lifetime audit](aton_lifetime_audit.md) separately checks every
row and every saved bin in the ten z6.15/z5.95 files. The new 41-file audit
does not constitute a full RT rerun or an all-row check of the other families.

## Density rescaling and the original z6.15 mismatch

The .npy columns are distance, xHI, temperature, three helium fractions,
and optical depth. They do not directly store total hydrogen density.
Consequently, unit optical-depth amplitude alone must not be interpreted
as proof of the density used during the chemistry calculation.

That evidence instead comes from the previously completed full RT reruns:
the z6.15 archive's evolved fields and spectra reproduce with the original
physical nH/nHI inputs at 6.148904, runtime redshift 6.604205, and **no
redshift density rescaling**. Across LOS 0 at all five lifetimes, the maximum
absolute flux discrepancy is 3.76e-9 and maximum relative temperature
discrepancy is 1.48e-10. See `outputs/lifetime_redshift_audit/z6p15_rt_all_lifetimes.json`
and the commands alongside it.

Relevant density factors, holding the existing legacy X convention fixed:

| Conversion of z6.148904 input | Required nH and nHI factor |
| --- | --- |
| To the observed quasar z6.1225 | 0.9889605741 |
| To the nearest-snapshot quasar z6.148904 | 1 |
| To the inferred runtime z6.604205 | 1.2034915735 |

Thus the detected 6.604205 configuration is neither the quasar's nearest
snapshot nor a consistently density-rescaled version at 6.604205. A
consistent nearest-snapshot run for the z6.15 input is exactly the same-z
case previously compared against the archive.

The z5.95 reproduction should be described narrowly: it validates an
unrescaled z5.945518 baseline. If intended for J1509, that input would need
a density factor of 1.0904467737 to reach the nearest-snapshot quasar
redshift 6.148904, together with corresponding target-redshift geometry
and spectral synthesis. Its numerical reproduction does not validate that
additional mapping.

## Later local files and provenance limits

There are later 2048-cell arrays in `fireqso_dominika/updated_aug23`, including
folders named `dense_rescaled` and `m26p58_rescaledrho`. Those names alone
do not establish that the entire paper procedure was applied.

As a supplementary check, the first row of all eight 1 Myr files in
`updated_aug23/m26p58_rescaledrho` reconstructs spectrally at its **input
snapshot redshift**, with maximum relative tau error below 4.56e-9 when
using the full-box legacy convention 235294.1176470588 and 256 redward pixels.
The fitted redshifts are 5.109794, 5.755954, 5.945518, 6.148904, 6.604205,
7.139393, 8.150086, and 9.545964. The small fitted amplitude offset from
unity (about 5e-9) is negligible numerical disagreement. M1450 = -26.58
corresponds to J1120+0641 at z7.0848, nearest local snapshot 7.139393.
That eight-file check concerns spectra only. A subsequent
[full RT check of the z5.76 file](rescaledrho_m26p58_reproduction.md) reproduces
four complete sightlines with RT densities rescaled to the **observed 7.0848**
and the legacy input X=0.76 factor removed, while runtime geometry and spectral
synthesis remain at **5.755954**. No final-paper provenance is established.

The later `updated_aug23/dense_rescaled/qso_z5.95_m-27.09_tq1e+00_los.npy`
also has 2048 cells and different retained rows from the 891-cell archive.
Its first row fits spectral redshift 5.945518, with maximum relative tau
error 5.48e-9 under the full-box convention. It therefore does not provide
evidence of spectral synthesis at J1509's redshift merely because its
folder name contains `rescaled`.

The historical intent of the fixed 6.604205 grid is therefore unresolved.
Possibilities include an interim common-reference grid or a configuration
mistake. The earlier stale-executable explanation remains a hypothesis,
not a recovered history. Identifying the actual final paper inputs and any
subsequent transformations requires their run or assembly provenance.

## Reproduce the audit

```bash
/Users/ssatyavolu/miniconda3/bin/python tools/check_aton_paper_mapping.py \
  --archive-dir /Users/ssatyavolu/projects/proximity-p1d/proximity_zone_p1d/data/simulations/ATON \
  --input-dir /Users/ssatyavolu/projects/fireqso_dominika \
  --output outputs/paper_mapping_audit/archive_mapping.json
```

For the supplementary later-folder check, use that folder as `--archive-dir`,
add `--boxsize 235294.1176470588`, and choose a separate output path.
Machine-readable results are in `outputs/paper_mapping_audit/`.
