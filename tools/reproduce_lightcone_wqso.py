#!/usr/bin/env python3
"""Regenerate one complete archived lightcone and test the requested FIRE parameters."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import subprocess
import time

import numpy as np

from audit_lightcone_wqso import ROOT, SOURCE, OUT, FIELDS, write_json, gas_metrics
from diagnose_aton_reproduction import comparison_metrics, optical_depth

ZS = np.array([4.278114,4.379683,4.598184,4.971176,5.109794,5.256697,
               5.412706,5.578765,5.755954,6.148904])
BOX = 927389.7


def reconstruct_inputs(los, n):
    z = np.loadtxt(SOURCE/'zs8072_z6.14.txt')[:n]
    k = np.clip(np.searchsorted(ZS,z)-1,0,len(ZS)-2)
    w = (z-ZS[k])/(ZS[k+1]-ZS[k])
    i = np.arange(n)
    fields = {}
    samples = {}
    for key in ('nh','nh1','T','rho','vpec'):
        samples[key] = np.stack([np.loadtxt(SOURCE/'stitchlos2'/f'z{zz:.2f}M3LOS_{key}long.txt',
                                           usecols=los,max_rows=n) for zz in ZS])
    for key in ('nh','xHI','T','rho','vpec'):
        log = key in ('xHI','T','rho')
        a = samples['nh1']/samples['nh'] if key=='xHI' else samples[key]
        if log:
            a = np.log10(a)
        a = (1-w)*a[k,i]+w*a[k+1,i]
        fields[key] = 10**a if log else a
    fields['nh1'] = fields.pop('xHI')*fields['nh']
    grids = {key:np.loadtxt(SOURCE/f'f1kz6.14M3LOS_{key}.txt',usecols=los,max_rows=n) for key in FIELDS}
    verified = {}
    for key in ('nh','nh1','T'):
        error = float(np.max(abs(fields[key]/grids[key]-1)))
        if error > 1e-9:
            raise ValueError(f'Reconstructed {key} differs from supplied lightcone input: {error}')
        verified[key] = error
    # Existing rho/vpec files have 8192 rows and use the older untrimmed stitching.
    # Regenerate the fields for the same 8072-cell stitching used by nH/nHI/T.
    for key in ('rho','vpec'):
        grids[key] = fields[key]
    write_json(OUT/f'input_reconstruction_los{los}.json',dict(
        local_redshift_file=str(SOURCE/'zs8072_z6.14.txt'),
        snapshot_redshifts=ZS.tolist(),initial_field_max_relative_difference=verified,
        method='Linear interpolation of physical nH and vpec; log interpolation of xHI, T, comoving rho, using stitchlos2 snapshots'))
    return grids


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--los',type=int,default=0)
    args=parser.parse_args()
    los=args.los
    OUT.mkdir(parents=True,exist_ok=True)
    target_path=SOURCE/'taualpha_wqso'/f'z6.14LOS{los}tq1e+00mag-27.00_tep_1e+00Myr1.0fd.txt'
    target=np.loadtxt(target_path)
    n=len(target)
    grids=reconstruct_inputs(los,n)
    inp=OUT/'full_inputs'/f'los{los}'
    inp.mkdir(parents=True,exist_ok=True)
    for key,value in grids.items():
        np.savetxt(inp/f'los{los}_{key}.txt',value,fmt='%.18e')
    build=['cc','-O3','-march=native','-std=c11','-DRTE_MAX_CELLS=8192','main.c']
    build += [f'src/{name}.c' for name in ['evolve_cell','lya_spectra','params','photoion','proximity_zone',
                                         'rates','save_outputs','thermochemistry','utils']]
    build += ['-lm','-o',str(OUT/'rte1d_lightcone')]
    subprocess.run(build,cwd=ROOT,check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    write_json(OUT/f'reference_los{los}.json',dict(path=str(target_path),shape=list(target.shape),
        sha256=hashlib.sha256(target_path.read_bytes()).hexdigest(),build_command=build,
        n_cells=n,spectral_boxsize=BOX,redward_pixels=0,quasar_lifetime_Myr=1,
        scope='Full RT and spectrum for this one sightline; all candidates use the same available lightcone initial fields.'))
    cases=[('labelled_z614_m27',6.14,-27.),
           ('requested_z615_m2709',6.15,float(np.float32(-27.09))),
           ('snapshot_z6148904_m2709',6.148904,float(np.float32(-27.09)))]
    jobs=[]
    for name,z,mag in cases:
        work=OUT/'full_runs'/name/f'los{los}'
        work.mkdir(parents=True,exist_ok=True)
        cmd=[str(OUT/'rte1d_lightcone'),'1','1','1','0',str(mag),str(z),'6.14',str(inp),
             f'los{los}',str(work),str(n),'1',str(los),'1','0','1',str(BOX),'1']
        jobs.append(dict(case=name,input_los=los,runtime_redshift=z,magnitude=mag,
                         command=cmd,log=str(work/'run.log'),
                         output=str(work/f'z6.14LOS{los}tq1mag{mag:.2f}_tep_1Myr1fd.txt')))
    write_json(OUT/f'full_manifest_los{los}.json',jobs)

    def run(job):
        start=time.monotonic()
        with open(job['log'],'w') as log:
            subprocess.run(job['command'],stdout=log,stderr=subprocess.STDOUT,check=True)
        candidate=np.loadtxt(job['output'])
        np.save(Path(job['log']).parent/'los.npy',candidate)
        result=dict(job,elapsed_seconds=time.monotonic()-start,
                    comparison=comparison_metrics(target,candidate),gas_comparison=gas_metrics(target,candidate))
        spectral_error=0.
        for start_bin in range(0,n,128):
            bins=np.arange(start_bin,min(start_bin+128,n))
            tau=optical_depth(candidate,grids['rho'],grids['vpec'],job['runtime_redshift'],bins,BOX,0)
            spectral_error=max(spectral_error,float(np.max(abs(np.exp(-tau)-np.exp(-candidate[bins,6])))))
        result['independent_spectral_max_flux_error']=spectral_error
        write_json(Path(job['log']).parent/'comparison.json',result)
        print(job['case'],'elapsed',result['elapsed_seconds'],'max flux difference',
              result['comparison']['max_absolute_flux_error'],'median xHI ratio',
              result['gas_comparison']['median_ratio'][0],flush=True)
        return result
    results=[]
    with ThreadPoolExecutor(max_workers=3) as pool:
        for future in as_completed([pool.submit(run,job) for job in jobs]):
            results.append(future.result())
            write_json(OUT/f'full_results_los{los}.json',results)


if __name__=='__main__':
    main()
