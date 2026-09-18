#!/usr/bin/env python3
"""Full RT candidate checks for the z5.76, -26.58, 1 Myr rescaledrho array."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import subprocess
import time

import numpy as np

from diagnose_aton_reproduction import comparison_metrics, optical_depth
from generate_inputs_from_dom import (add_backgrounds, field_block, halo_pixel_offset,
                                      header_values, shift_to_halo, hydrogen_number_density_scale,
                                      RHO_SCALE)

ROOT = Path(__file__).resolve().parents[1]
SOURCE = Path('/Users/ssatyavolu/projects/fireqso_dominika')
ARCHIVE = SOURCE/'updated_aug23/m26p58_rescaledrho'
TARGET = ARCHIVE/'qso_z5.76_m-26.58_tq1e+00_los.npy'
OUT = ROOT/'outputs/rescaledrho_m26p58_z5p76'
ZS, ZOBS, ZNEAR = 5.755954, 7.0848, 7.139393
BOX = 235294.118
FIELDS = ('nh', 'nh1', 'T', 'vpec', 'rho', 'Gammabg', 'Gammabghe1', 'Gammabghe2', 'Heatbg')
DOBS, DNEAR = ((1+ZOBS)/(1+ZS))**3, ((1+ZNEAR)/(1+ZS))**3
CASES = {
    'unscaled': (1., ZS, 1., False, 1),
    'density_nearest': (DNEAR, ZS, 1., False, 1),
    'density_nearest_Xcorrected': (DNEAR/.76, ZS, 1., False, 1),
    'density_observed': (DOBS, ZS, 1., False, 1),
    'density_observed_Xcorrected': (DOBS/.76, ZS, 1., False, 1),
    'unscaled_Xcorrected': (1/.76, ZS, 1., False, 1),
    'density_nearest_helium_fixed': (DNEAR, ZS, 1., False, 0),
    'density_nearest_bg_recomputed': (DNEAR, ZS, 1., True, 1),
    'literal_observed': (DOBS, ZOBS, ((1+ZOBS)/(1+ZS))**2, False, 1),
    'literal_observed_Xcorrected': (DOBS/.76, ZOBS, ((1+ZOBS)/(1+ZS))**2, False, 1),
    'literal_nearest': (DNEAR, ZNEAR, ((1+ZNEAR)/(1+ZS))**2, False, 1),
}


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def source_grids(los):
    prefix = f"{'f' if los<500 else 'l'}1kz5.76M3LOS"
    grids = {key: np.loadtxt(SOURCE/f'{prefix}_{key}.txt', usecols=los%500) for key in FIELDS}
    raw = np.memmap(SOURCE/'dom_zloc_M3_spec_z5.756.dat', dtype='<f8', mode='r')
    offset = halo_pixel_offset(header_values(raw, los)[0])
    direct = {key: shift_to_halo(field_block(raw, los, idx), offset)
              for key, idx in [('T',1), ('vpec',2), ('nh',3), ('nh1',4)]}
    direct['rho'] = direct['nh']*RHO_SCALE
    for key in ['nh','nh1']:
        direct[key] *= hydrogen_number_density_scale(ZS)
    add_backgrounds(direct)
    errors = {}
    for key in FIELDS:
        np.testing.assert_allclose(grids[key], direct[key], rtol=1e-10, atol=0)
        errors[key] = float(np.max(np.abs(grids[key]/np.where(direct[key]!=0,direct[key],1)-1))) if key!='vpec' else float(np.max(abs(grids[key]-direct[key])))
    return grids, errors


def prepare(case, row, los, original):
    density_factor, zrun, length_factor, recalc, helium_bug = CASES[case]
    grids = {key: value.copy() for key,value in original.items()}
    for key in ['nh','nh1']:
        grids[key] *= density_factor
    np.testing.assert_allclose(grids['nh1']/grids['nh'], original['nh1']/original['nh'], rtol=5e-16)
    if recalc:
        add_backgrounds(grids)
    work = OUT/case/f'row{row}'
    inp, rt = work/'inputs', work/'rt'
    inp.mkdir(parents=True, exist_ok=True)
    rt.mkdir(exist_ok=True)
    for key,value in grids.items():
        np.savetxt(inp/f'los{los}_{key}.txt', value, fmt='%.18e')
    cmd = [str(ROOT/'rte1d_opt'), '1', '1', '1', '0', f'{float(np.float32(-26.58)):.17g}',
           str(zrun), '5.76', str(inp), f'los{los}', str(rt), '2048', '1', str(los),
           '1', '256', str(helium_bug), str(BOX), f'{length_factor:.17g}']
    return dict(case=case, target_row=row, global_input_los=los, runtime_redshift=zrun,
                density_factor=density_factor, comoving_length_factor=length_factor,
                backgrounds_recomputed=recalc, legacy_alphahe2_bug=helium_bug,
                command=cmd, log=str(work/'run.log'),
                output=str(rt/f'z5.76LOS{los}tq1mag-26.58_tep_1Myr1fd.txt'))


def run(job):
    start = time.monotonic()
    with open(job['log'],'w') as log:
        subprocess.run(job['command'], stdout=log, stderr=subprocess.STDOUT, check=True)
    result = dict(job, elapsed_seconds=time.monotonic()-start)
    candidate = np.loadtxt(job['output'])
    target = np.load(TARGET, mmap_mode='r')[job['target_row']]
    result['comparison'] = comparison_metrics(target, candidate)
    result['median_final_xHI_ratio'] = float(np.median(candidate[:,1]/target[:,1]))
    result['median_final_temperature_ratio'] = float(np.median(candidate[:,2]/target[:,2]))
    inp = Path(job['command'][8])
    rho, velocity = [np.loadtxt(inp/f"los{job['global_input_los']}_{key}.txt") for key in ('rho','vpec')]
    tau = optical_depth(candidate,rho,velocity,job['runtime_redshift'],np.arange(2048),BOX*job['comoving_length_factor'],256)
    result['independent_spectral_max_flux_difference'] = float(np.max(abs(np.exp(-tau)-np.exp(-candidate[:,6]))))
    write_json(Path(job['log']).parent/'comparison.json', result)
    print(f"{job['case']} row {job['target_row']}, {result['elapsed_seconds']:.1f}s: "
          f"max flux diff {result['comparison']['max_absolute_flux_error']:.6g}, "
          f"max relative T diff {result['comparison']['max_relative_error'][2]:.6g}, "
          f"median xHI ratio {result['median_final_xHI_ratio']:.6g}", flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cases', nargs='+', choices=CASES, default=['unscaled','density_nearest','density_nearest_Xcorrected'])
    parser.add_argument('--rows', nargs='+', type=int, default=[0])
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    target = np.load(TARGET, mmap_mode='r')
    deleted = np.loadtxt(ARCHIVE/'deleted_indices_los_z5.76_mfp_t1e+00Myr_mag-26.58.txt', ndmin=1).astype(int)
    kept = np.delete(np.arange(1000), deleted)
    if target.shape != (len(kept),2048,7):
        raise ValueError(f'Unexpected shape {target.shape}')
    np.savetxt(OUT/'kept_indices.txt', kept, fmt='%d')
    write_json(OUT/'reference.json', dict(path=str(TARGET), shape=list(target.shape),
               sha256=hashlib.sha256(TARGET.read_bytes()).hexdigest(), z_sim=ZS,
               observed_quasar_z=ZOBS, nearest_quasar_z=ZNEAR, magnitude=-26.58, tq_Myr=1))
    jobs = []
    for row in args.rows:
        los = int(kept[row])
        grids, errors = source_grids(los)
        write_json(OUT/f'input_verification_los{los}.json',errors)
        for case in args.cases:
            jobs.append(prepare(case,row,los,grids))
    name = '_'.join(args.cases)+'_rows'+'_'.join(map(str,args.rows))
    write_json(OUT/f'manifest_{name}.json',jobs)
    completed = []
    with ThreadPoolExecutor(max_workers=3) as pool:
        for future in as_completed([pool.submit(run,job) for job in jobs]):
            completed.append(future.result())
            write_json(OUT/f'results_{name}.json',completed)


if __name__ == '__main__':
    main()
