#!/usr/bin/env python3
"""Compare near-quasar RT fields with the archived taualpha_wqso lightcones.

These probes evolve the first 891 cells. Their truncated optical depths are
not compared with the long lightcone; only the causally independent RT fields
are compared at this stage.
"""
import json
from pathlib import Path
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCE = Path('/Users/ssatyavolu/projects/meanfreepath/mean_free_path/lightcone')
OUT = ROOT/'outputs/lightcone_taualpha_wqso'
FIELDS = ('nh','nh1','T','vpec','rho','Gammabg','Gammabghe1','Gammabghe2','Heatbg')
CASES = {
    'label_z614_m27_legacy_density': (6.14, -27., 1.),
    'requested_z615_m2709_legacy_density': (6.148904, float(np.float32(-27.09)), 1.),
    'label_z614_m27_Xcorrected': (6.14, -27., 1/.76),
}


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def gas_metrics(reference, reproduced):
    delta = abs(reproduced[:,1:6]-reference[:,1:6])
    rel = delta/np.maximum(abs(reference[:,1:6]), 1e-40)
    return dict(columns=['xHI','T','xHeI','xHeII','xHeIII'],
                maximum_absolute_difference=delta.max(axis=0).tolist(),
                maximum_relative_difference=rel.max(axis=0).tolist(),
                median_relative_difference=np.median(rel,axis=0).tolist(),
                median_ratio=np.median(reproduced[:,1:6]/np.maximum(reference[:,1:6],1e-40),axis=0).tolist())


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    n, los = 891, 0
    target_path = SOURCE/'taualpha_wqso'/f'z6.14LOS{los}tq1e+00mag-27.00_tep_1e+00Myr1.0fd.txt'
    target = np.loadtxt(target_path)
    original = {key: np.loadtxt(SOURCE/f'f1kz6.14M3LOS_{key}.txt',usecols=los,max_rows=n) for key in FIELDS}
    jobs = []
    for name,(z,mag,dscale) in CASES.items():
        work = OUT/'probes'/name
        inp,rt = work/'inputs',work/'rt'
        inp.mkdir(parents=True, exist_ok=True)
        rt.mkdir(exist_ok=True)
        for key,value in original.items():
            np.savetxt(inp/f'los{los}_{key}.txt',value*dscale if key in ('nh','nh1') else value,fmt='%.18e')
        cmd=[str(ROOT/'rte1d_opt'),'1','1','1','0',str(mag),str(z),'6.14',str(inp),f'los{los}',str(rt),str(n),'1',str(los),'1','0','1']
        jobs.append(dict(case=name,runtime_redshift=z,magnitude=mag,density_multiplier=dscale,
                         input_los=los,target=str(target_path),target_shape=list(target.shape),n_rt_cells=n,
                         command=cmd,log=str(work/'run.log'),
                         output=str(rt/f'z6.14LOS{los}tq1mag{mag:.2f}_tep_1Myr1fd.txt')))
    write_json(OUT/'probe_manifest.json',jobs)

    def run(job):
        start=time.monotonic()
        with open(job['log'],'w') as log:
            subprocess.run(job['command'],stdout=log,stderr=subprocess.STDOUT,check=True)
        a=np.loadtxt(job['output'])
        result=dict(job,elapsed_seconds=time.monotonic()-start,gas_comparison=gas_metrics(target[:n],a))
        write_json(Path(job['log']).parent/'comparison.json',result)
        print(job['case'],json.dumps(result['gas_comparison']),flush=True)
        return result
    done=[]
    with ThreadPoolExecutor(max_workers=3) as pool:
        for f in as_completed([pool.submit(run,j) for j in jobs]):
            done.append(f.result())
            write_json(OUT/'probe_results.json',done)


if __name__=='__main__':
    main()
