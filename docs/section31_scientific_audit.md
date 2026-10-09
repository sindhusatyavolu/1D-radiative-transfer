# Scientific audit of the ATON Section 3.1 workflow

## Decision

New production grids use a fixed comoving skewer and evaluate the physics at
the quasar redshift. For an ATON snapshot at `z_sim` and a quasar at `z_qso`,

```text
A = (1 + z_qso) / (1 + z_sim)
comoving length multiplier = 1
physical nH and nHI multiplier relative to the legacy text files = A^3 / 0.76
background Gamma multiplier = A^3 / 0.76
background heating multiplier = (A^3 / 0.76)^2
Ly-alpha density = X times the supplied total-gas rho
Ly-alpha spatial unit = stored ckpc (not ckpc/h)
Ly-alpha Hubble constant = 100h km/s/Mpc
runtime redshift = z_qso
```

The comoving mass-density field, peculiar velocities, temperature, and initial
ionization fractions remain those of the selected snapshot.

## Output-coordinate convention

The assembled arrays retain the historical `(1000, 891, 7)` layout. Columns
1--5 are gas-cell fields at the comoving coordinate in column 0. The Ly-alpha
optical depth in column 6 is a spectral-bin field with 256 pixels redward of
the quasar, so its zero-Hubble-velocity origin is index 256. Its proper-distance
coordinate is

```text
r_tau = (column0 - column0[256]) / (1 + z_qso).
```

Every assembled array has a JSON sidecar that records this origin and formula.
The offset is inherited from the legacy format; it is distinct from the
redshift-rescaling issue and does not change the 891-pixel array length.

## Why the length scale stays fixed

For a fixed comoving cell width `dx_com`, its proper width at the quasar is

```text
dx_proper(z_qso) = dx_com / (1 + z_qso).
```

Defining `A=(1+z_qso)/(1+z_sim)=a_sim/a_qso` is algebraically correct. The
direction of the mapping is the important part. For the same comoving cell,

```text
dx_proper(z_qso) / dx_proper(z_sim) = a_qso / a_sim = 1/A.
```

Multiplication by `A` describes the reverse proper-length mapping. Section 3.2
says that CROC uses the same expansion factor, but does not add information
that resolves this direction issue.

The solver already performs this conversion through its runtime scale factor.
It also uses `H(z_qso) dx_com/(1+z_qso)` for the Hubble velocity grid. Therefore
an additional comoving multiplier would count the redshift conversion again.

[Durovcikova et al. (2024), Section 3.1](https://arxiv.org/html/2401.10328v1#S3.SS1)
prints a physical-length factor `(1+z_qso)/(1+z_sim)`. Implementing that
sentence literally in this comoving solver requires a comoving factor `A^2`.
That was the convention used by the previous wrapper. It is opposite to the
standard relation between proper and comoving separation for a fixed skewer;
see [Hogg (1999), Distance measures in cosmology](https://arxiv.org/abs/astro-ph/9905116).
It also does not reproduce the recovered historical ATON arrays, whose geometry
was left at the snapshot redshift. The production choice here follows the
coordinate relation and is explicitly recorded as a scientific correction,
rather than being presented as exact reproduction of the archived paper files.

## Density and background normalization

The binary ATON `cube_D` field is already comoving hydrogen number density.
The opening note in `fireqso_dominika/get_los.py` records that the legacy text
generator multiplied both `nH` and `nHI` by an extra hydrogen mass fraction
`X=0.76`. Removing it preserves `xHI=nHI/nH`.

The desired physical density is proportional to `(1+z_qso)^3`, giving the
combined multiplier `A^3/0.76` relative to those text files. The `rho` input is
the correctly converted comoving total-gas mass density and remains unchanged.
When converting it to hydrogen number density, Ly-alpha synthesis must include
the hydrogen mass fraction `X=0.76`. The legacy code omitted this factor. The
solver applies the scale-factor dependence at `z_qso` internally.

The supplied background rates were derived cell by cell from the legacy
`nH`, `nHI`, and temperature under photoionization equilibrium. At fixed
fractions and temperature, all three photoionization rates scale linearly with
density and the saved background heating term scales as density squared.
Scaling these fields with the corrected density prevents the input neutral
fraction from drifting solely because of the remapping.

## Ly-alpha `h` units

The ATON grid spacing is `78.125 ckpc/h`. With the simulation value `h=0.678`,
this is `115.228613569 ckpc`. The legacy constant used `78.125/0.68 =
114.889705882 ckpc`; production now uses the precise simulation value. The
`rho` field uses the Gadget mass-density unit proportional to `h^2`. The
legacy Ly-alpha code nevertheless divided the density particle count by
another `h`, used `H0=100 km/s/Mpc`, and divided the Hubble width by `h`.
Those operations belong to a grid still expressed in ckpc/h and are
inconsistent with the stored ckpc coordinate.

The production calculation now converts `rho` with its `h^2` unit directly and
uses `H0=100h km/s/Mpc`. At `z_qso=5.98`, one 115.2286 ckpc cell therefore
corresponds to about 11.49 km/s of Hubble velocity. Holding the older cell
width fixed, the legacy expression gave 24.928 km/s instead of 11.459 km/s,
larger by `1/h^2=2.1754`. This correction makes the internal
spectral grid agree with the proper-distance and velocity axes used when the
arrays are compared to data.

## Numerical checks

The audit used snapshot `z_sim=5.755954`, quasar `z_qso=7.0848`, `M1450=-26.58`,
`tQ=1 Myr`, and LOS 0.

- Changing only the previous `A^2` comoving multiplier to a fixed comoving grid
  changed the transmitted flux by mean `|dF|=0.07390`, maximum `0.81570`; 38.2%
  of saved pixels changed by more than 0.01.
- Scaling the backgrounds consistently changed the same representative quasar
  profile by mean `|dF|=0.01165`, maximum `0.14719`.
- Correcting the total-gas `rho` to hydrogen normalization multiplies every
  Ly-alpha optical depth by exactly `0.76` when the legacy spatial-unit
  convention is held fixed; for the representative profile this changed mean
  `|dF|=0.02706` and maximum `0.10064`.
- Correcting the remaining Ly-alpha `h` factors changed a `z_qso=5.98`,
  `z_sim=5.95`, `tQ=0.01 Myr` control by mean `|dF|=0.07366` and maximum
  `0.77636`; 66.9% of pixels changed by more than 0.01. Re-enabling the legacy
  branch reproduced the pre-correction output exactly.
- The legacy final-time condition could stop one full global step early. A
  `tQ=100 Myr` high-neutral control stopped at 99 Myr. The corrected condition
  reaches 100 Myr; in that sightline the change had mean `|dF|=3.99e-4` and
  maximum `0.00749`.
- Replacing the legacy spectral box `102441.268 ckpc` with exactly
  `891*cellsize` changed mean flux by `9.43e-5` and maximum flux by `0.00387`
  when the legacy cell width was held fixed.
- Replacing the rounded `78.125/0.68` cell width with `78.125/0.678` changed
  mean `|dF|=6.27e-4` and maximum `0.0132` in a `z_qso=5.98`, `tQ=0.01 Myr`
  control. Correcting the Myr conversion for the cosmological clock had a
  negligible flux effect in that short-lifetime control (`<1.6e-7` maximum).
- The optimized and non-cached executables compiled from the same source were
  bitwise identical in all seven output columns.
- The final high-level command, including the corrected Ly-alpha units, agrees
  with an independently pre-scaled input control to maximum flux difference
  `1.7e-12`.
- With a negligible quasar, the median `xHI` after 1 Myr changed by 0.23% for
  the `z=5.11` snapshot and less than `1e-8` fractionally for the `z=9.55`
  snapshot when both were mapped to `z_qso=5.98`.
- A final production-default smoke test for the most demanding corner,
  `z_sim=9.55`, `z_qso=5.98`, and `tQ=100 Myr`, reached the exact requested
  final time and produced a finite 891-row output.

Audit artifacts are stored under `outputs/scientific_convention_audit/`.

## Remaining model approximation

Using a snapshot from another redshift as a neutral-fraction model is a
phenomenological redshift-matching procedure: it preserves that snapshot's
ionization topology, temperature, density contrast, and peculiar velocity
while changing the mean physical density and quasar/cosmological conversion.
It is not a self-consistent hydrodynamic simulation evolved to `z_qso`.

The 1D post-processing also treats the cosmological redshift as fixed during a
quasar episode, including a 100 Myr episode. This is inherited from the model
and should be stated when the grid is used. These approximations are distinct
from the corrected coordinate and density bookkeeping above.
