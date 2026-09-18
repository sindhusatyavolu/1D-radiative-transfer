#!/usr/bin/env python3
"""Summarize completed, matched RT runs for the selected rescaledrho file."""
import json
import os

import numpy as np

from audit_rescaledrho_reproduction import OUT, TARGET, ZS, ZOBS, DOBS, BOX


def velocity_grid(z, length_scale=1.):
    a = 1/(1+z)
    wmax = a * 3.2407789e-18 * 3.085678e21/1e5 * np.sqrt(.308/a**3+.692) * BOX*length_scale/.678
    return (np.arange(2048)-256) * wmax/2048


def main():
    rows = np.array([0, 331, 662, 993])
    kept = np.loadtxt(OUT/'kept_indices.txt', dtype=int)
    reference = np.load(TARGET, mmap_mode='r')[rows]
    comparisons, reproduced = [], []
    for row in rows:
        record = json.loads((OUT/'density_observed_Xcorrected'/f'row{row}'/'comparison.json').read_text())
        result = np.loadtxt(record['output'])
        if not np.all(np.isfinite(result)):
            raise ValueError('Nonfinite RT output')
        if record['comparison']['max_absolute_flux_error'] > 1e-7:
            raise ValueError(f'Reproduction failed for row {row}')
        if record['independent_spectral_max_flux_difference'] > 1e-7:
            raise ValueError(f'Spectral check failed for row {row}')
        comparisons.append(record)
        reproduced.append(result)
    reproduced = np.array(reproduced)
    alternative = json.loads((OUT/'literal_observed_Xcorrected'/'row0'/'comparison.json').read_text())
    literal = np.loadtxt(alternative['output'])
    control_cases = ['unscaled', 'density_nearest', 'density_nearest_Xcorrected',
                     'density_observed', 'literal_observed_Xcorrected']
    controls = [json.loads((OUT/key/'row0'/'comparison.json').read_text()) for key in control_cases]
    vref = velocity_grid(ZS)
    vliteral = velocity_grid(ZOBS, ((1+ZOBS)/(1+ZS))**2)
    np.save(OUT/'reproduced_sample_los.npy', reproduced)
    np.savez_compressed(OUT/'matched_samples.npz', target_rows=rows, global_input_los=kept[rows],
                        reference=reference, reproduced=reproduced,
                        literal_section31_row0=literal, reference_velocity_kms=vref,
                        literal_velocity_kms=vliteral)
    report = dict(reference=json.loads((OUT/'reference.json').read_text()),
                  configuration=dict(input_snapshot_z=ZS, runtime_z=ZS, density_target_z=ZOBS,
                                     density_redshift_factor=DOBS, legacy_X_correction=1/.76,
                                     total_nh_nh1_multiplier=DOBS/.76,
                                     comoving_rho_multiplier=1., comoving_length_multiplier=1.,
                                     spectral_boxsize=BOX, redward_pixels=256,
                                     initial_T_xHI_vpec_unchanged=True,
                                     backgrounds='Original source-snapshot backgrounds',
                                     legacy_alphahe2_bug=1, magnitude_as_float32=float(np.float32(-26.58))),
                  full_RT_samples=comparisons, row0_controls=controls,
                  maximum_absolute_flux_difference=max(r['comparison']['max_absolute_flux_error'] for r in comparisons),
                  maximum_absolute_xHI_difference=max(r['comparison']['max_absolute_error'][1] for r in comparisons),
                  maximum_absolute_T_difference_K=max(r['comparison']['max_absolute_error'][2] for r in comparisons),
                  maximum_independent_spectral_flux_difference=max(r['independent_spectral_max_flux_difference'] for r in comparisons),
                  scope='Four complete 2048-cell RT sightlines from one 994-row file. Other rescaledrho files and full ensemble RT not checked here.')
    (OUT/'benchmark.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    os.environ.setdefault('MPLCONFIGDIR', str(OUT/'.matplotlib'))
    os.environ.setdefault('XDG_CACHE_HOME', str(OUT/'.cache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(4, 2, figsize=(12, 9), sharex=True, constrained_layout=True)
    for k, row in enumerate(rows):
        fref, frep = np.exp(-reference[k,:,6]), np.exp(-reproduced[k,:,6])
        axes[k,0].plot(vref, fref, color='black', lw=1.5, label='rescaledrho reference')
        axes[k,0].plot(vref, frep, color='#0072b2', lw=1., ls='--', label='Reproduction')
        axes[k,1].plot(vref, frep-fref, color='#0072b2', lw=.9)
        axes[k,0].set_ylabel(f'Row {row}\nTransmitted flux')
        axes[k,1].set_ylabel('Reproduction − reference')
        axes[k,0].set_ylim(-.03, 1.03)
        for ax in axes[k]:
            ax.set_xlim(-500, 5000)
            ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=9)
    for ax in axes[-1]:
        ax.set_xlabel('Velocity blueward of Lyα (km/s)')
    fig.suptitle('rescaledrho: zsim = 5.76, M1450 = −26.58, tQ = 1 Myr')
    for ext in ('png','pdf'):
        fig.savefig(OUT/f'reproduction_comparison.{ext}', dpi=160)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(10, 4.5), constrained_layout=True)
    ax.plot(vref, np.exp(-reference[0,:,6]), color='black', lw=1.5, label='rescaledrho reference')
    ax.plot(vref, np.exp(-reproduced[0,:,6]), color='#0072b2', lw=1., ls='--', label='Matched reproduction')
    ax.plot(vliteral, np.exp(-literal[:,6]), color='#d55e00', lw=1., label='Literal Section 3.1; zQSO = 7.0848')
    ax.set(xlim=(-500,5000), ylim=(-.03,1.03), xlabel='Velocity blueward of Lyα (km/s)',
           ylabel='Transmitted flux', title='Same initial skewer, row 0; density corrected in both reruns')
    ax.legend(fontsize=9)
    ax.grid(alpha=.2)
    for ext in ('png','pdf'):
        fig.savefig(OUT/f'literal_section31_comparison.{ext}', dpi=160)
    plt.close(fig)
    print(json.dumps({k:v for k,v in report.items() if k.startswith('maximum')}, indent=2))


if __name__ == '__main__':
    main()
