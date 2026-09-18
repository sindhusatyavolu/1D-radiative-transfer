#!/usr/bin/env python3
"""Make a reviewable comparison of the complete lightcone reproduction."""
import argparse
import json
import os
from pathlib import Path

import numpy as np

from audit_lightcone_wqso import SOURCE, OUT, write_json
from reproduce_lightcone_wqso import BOX


def velocity_grid(z,n):
    a=1/(1+z)
    span=a*3.2407789e-18*3.085678e21/1e5*np.sqrt(.308/a**3+.692)*BOX/.678
    return np.arange(n)*span/n


def metrics(d):
    return dict(mean_absolute=float(np.mean(abs(d))),max_absolute=float(np.max(abs(d))),
                rms=float(np.sqrt(np.mean(d*d))),mean_signed=float(np.mean(d)))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--los',type=int,default=0)
    args=parser.parse_args()
    los=args.los
    meta=json.loads((OUT/f'reference_los{los}.json').read_text())
    results=json.loads((OUT/f'full_results_los{los}.json').read_text())
    ref=np.loadtxt(meta['path'])
    n=len(ref)
    vref=velocity_grid(6.14,n)
    common=np.arange(0.,5000.1,5.)
    fref=np.exp(-ref[:,6])
    arrays=dict(reference=ref,reference_velocity_kms=vref,common_velocity_kms=common)
    report=dict(reference=meta,
                input_reconstruction=json.loads((OUT/f'input_reconstruction_los{los}.json').read_text()),
                convention='Same 8071-cell lightcone initial fields, 1 Myr, no extra density rescale or X factor; zero redward offset',
                cases=results)
    for result in results:
        a=np.loadtxt(result['output'])
        v=velocity_grid(result['runtime_redshift'],n)
        arrays[result['case']]=a
        arrays[result['case']+'_velocity_kms']=v
        if not np.isfinite(a).all():
            raise ValueError('Nonfinite result')
        if result['case']=='labelled_z614_m27':
            if result['comparison']['max_absolute_flux_error']>1e-6:
                raise ValueError('Labelled configuration failed to reproduce')
        if result['independent_spectral_max_flux_error']>1e-6:
            raise ValueError('Independent spectral validation failed')
        f=np.exp(-a[:,6])
        result['first635_same_pixel_flux']=metrics(f[:635]-fref[:635])
        result['common_velocity_0_5000_kms_flux']=metrics(np.interp(common,v,f)-np.interp(common,vref,fref))
    write_json(OUT/f'benchmark_los{los}.json',report)
    np.savez_compressed(OUT/f'comparison_arrays_los{los}.npz',**arrays)
    os.environ.setdefault('MPLCONFIGDIR',str(OUT/'.matplotlib'))
    os.environ.setdefault('XDG_CACHE_HOME',str(OUT/'.cache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(3,1,figsize=(10,9),constrained_layout=True)
    axes[0].plot(vref,fref,color='black',lw=1.8,label='Archived lightcone')
    colors={'labelled_z614_m27':'#0072b2','requested_z615_m2709':'#d55e00',
            'snapshot_z6148904_m2709':'#009e73'}
    labels={'labelled_z614_m27':'z=6.14, M=−27.00: reproduction',
            'requested_z615_m2709':'z=6.15, M=−27.09: requested',
            'snapshot_z6148904_m2709':'z=6.148904, M=−27.09: exact snapshot z'}
    for result in results:
        key=result['case']
        a=arrays[key]
        v=arrays[key+'_velocity_kms']
        f=np.exp(-a[:,6])
        axes[0].plot(v,f,color=colors[key],lw=.9,ls='--' if key.startswith('labelled') else '-',label=labels[key])
        if key.startswith('labelled'):
            axes[2].plot(vref,f-fref,color=colors[key],lw=1.)
        else:
            axes[1].plot(common,np.interp(common,v,f)-np.interp(common,vref,fref),color=colors[key],lw=1.,label=labels[key])
    for ax in axes:
        ax.set_xlim(0,5000)
        ax.grid(alpha=.2)
        ax.set_xlabel('Velocity blueward of the quasar (km/s)')
    axes[0].set_ylim(-.03,1.03)
    axes[0].set_ylabel('Transmitted Lyα flux')
    axes[0].legend(fontsize=8,loc='upper right')
    axes[1].set_ylabel('Requested − archive\n(common velocity grid)')
    axes[2].set_ylabel('Reproduction − archive\n(native pixels)')
    axes[1].axhline(0,color='gray',lw=.7)
    axes[2].axhline(0,color='gray',lw=.7)
    fig.suptitle(f'taualpha_wqso: complete 8071-cell RT runs; LOS {los}, tQ = 1 Myr')
    for ext in ('png','pdf'):
        fig.savefig(OUT/f'comparison_los{los}.{ext}',dpi=170)
    plt.close(fig)
    for result in results:
        print(result['case'],json.dumps(dict(all_pixels=result['comparison']['max_absolute_flux_error'],
              first635=result['first635_same_pixel_flux'],common=result['common_velocity_0_5000_kms_flux'])))


if __name__=='__main__':
    main()
