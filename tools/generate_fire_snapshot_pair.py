#!/usr/bin/env python3
"""Generate matched FIRE arrays with a quasar and before quasar turn-on.

Defaults correct only the archived z6.15 runtime redshift, retaining the
conventions validated against z5.95. Exact duplicate inputs are run once,
then restored to the archive's complete row order. Completed jobs are
resumable only when their input, executable, and command hashes match.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import subprocess
import time

import numpy as np

from diagnose_aton_reproduction import optical_depth

ROOT = Path(__file__).resolve().parents[1]
SOURCE = Path('/Users/ssatyavolu/projects/fireqso_dominika')
FIELDS = ('nh', 'nh1', 'T', 'vpec', 'rho', 'Gammabg', 'Gammabghe1', 'Gammabghe2', 'Heatbg')
DEFAULT_INDICES = ROOT/'outputs/recover_z6p15/inferred_kept_indices_z6.15_m27.09_from_target.txt'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path, value):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def validate_array(a, n):
    if a.shape != (n, 7) or not np.isfinite(a).all():
        raise ValueError(f'Invalid output shape or values: {a.shape}')
    if (a[:, 2] <= 0).any() or (a[:, 6] < 0).any():
        raise ValueError('Nonpositive temperature or negative optical depth')
    if (a[:, [1, 3, 4, 5]] < -1e-12).any() or (a[:, [1, 3, 4, 5]] > 1+1e-12).any():
        raise ValueError('Ionization fraction outside [0, 1]')
    # The legacy evolved helium columns do not close exactly to unity.
    # Record their closure error separately; do not silently renormalize them.


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, default=SOURCE)
    parser.add_argument('--indices-file', type=Path, default=DEFAULT_INDICES)
    parser.add_argument('--output-dir', type=Path, default=ROOT/'outputs/corrected_z6p15_m2709')
    parser.add_argument('--rte', type=Path, default=ROOT/'rte1d_opt')
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--tq', type=float, default=1.)
    parser.add_argument('--prepare-only', action='store_true')
    args = parser.parse_args()
    if args.workers < 1 or not np.isfinite(args.tq) or args.tq <= 0:
        parser.error('workers and tq must be positive')
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    rte = args.rte.resolve()
    ids = np.atleast_1d(np.loadtxt(args.indices_file, dtype=int))
    if len(ids) == 0 or (ids < 0).any() or (ids >= 1000).any():
        raise ValueError('Expected global input LOS IDs from 0 to 999')
    n, z, label, box, extra = 891, 6.148904, 6.15, 102441.268, 256
    mag = float(np.float32(-27.09))
    grids, inputs = {}, []
    for key in FIELDS:
        arrays = []
        for prefix in ('f1kz6.15M3LOS', 'l1kz6.15M3LOS'):
            path = args.input_dir.resolve()/f'{prefix}_{key}.txt'
            a = np.loadtxt(path, max_rows=n)
            if a.shape != (n, 500):
                raise ValueError(f'Unexpected shape {a.shape}: {path}')
            arrays.append(a)
            inputs.append(dict(path=str(path), sha256=digest(path)))
        grids[key] = np.concatenate(arrays, axis=1)
        if not np.isfinite(grids[key][:, ids]).all():
            raise ValueError(f'Nonfinite selected input {key}')
    if (grids['nh'][:, ids] <= 0).any() or (grids['T'][:, ids] <= 0).any():
        raise ValueError('Selected density/temperature must be positive')
    xhi = grids['nh1'][:, ids]/grids['nh'][:, ids]
    if (xhi < 0).any() or (xhi > 1).any():
        raise ValueError('Selected initial xHI outside [0, 1]')

    config = dict(
        snapshot_redshift=z, runtime_redshift=z, filename_redshift=label,
        m1450=-27.09, runtime_magnitude=mag, tq_Myr=args.tq, duty_cycle=1,
        n_cells=n, n_rows=len(ids), legacy_density=1, density_multiplier=1.,
        legacy_alphahe2_bug=1, num_pixels_redward=extra,
        spectral_boxsize=box, comoving_length_scale=1,
        columns=['distance_ckpc', 'xHI', 'temperature_K', 'xHeI', 'xHeII', 'xHeIII', 'tau_Lya'],
        no_quasar_definition='Initial snapshot gas, no time evolution; includes the existing IGM ionization state.',
        initial_helium='Solver initialization: xHeII=xHeIII=1e-20; xHeI=1-2e-20.',
        spectral_coordinate='tau has a 256-pixel redward offset relative to the gas columns; retained for archive compatibility.',
        density_note='Original FIRE nH/nHI, including the extra 0.76 factor; this run corrects runtime redshift only.',
        input_files=inputs, indices_file=str(args.indices_file.resolve()),
        indices_sha256=digest(args.indices_file), executable_sha256=digest(rte),
    )
    config_id = hashlib.sha256(json.dumps(config, sort_keys=True).encode()).hexdigest()
    manifest_path = out/'manifest.json'
    if manifest_path.exists() and json.loads(manifest_path.read_text())['configuration_sha256'] != config_id:
        raise ValueError('Output directory contains a different configuration; use another --output-dir')

    groups = {}
    for los in ids:
        h = hashlib.sha256(np.stack([grids[key][:, los] for key in FIELDS]).tobytes()).hexdigest()
        groups.setdefault(h, []).append(int(los))
    reps = [group[0] for group in groups.values()]
    representative = {los: group[0] for group in groups.values() for los in group}
    row_reps = np.array([representative[int(los)] for los in ids])
    np.savetxt(out/'input_los_ids.txt', ids, fmt='%d')
    np.savetxt(out/'representative_input_los_ids.txt', row_reps, fmt='%d')
    # Preserve the physical and spectral inputs alongside the output metadata.
    np.savez(out/'source_inputs.npz', **grids)
    manifest = dict(configuration_sha256=config_id, config=config,
                    unique_inputs=len(reps), input_los_ids=ids.tolist(),
                    representative_by_row=row_reps.tolist(),
                    deduplication='Exact SHA256 of all nine 891-cell input arrays; no approximate matching.')
    write_json(manifest_path, manifest)
    jobs = []
    for h, group in groups.items():
        los = group[0]
        inp = out/'inputs'/f'los{los}'
        work = out/'runs'/f'los{los}'
        inp.mkdir(parents=True, exist_ok=True)
        work.mkdir(parents=True, exist_ok=True)
        for key in FIELDS:
            np.savetxt(inp/f'los{los}_{key}.txt', grids[key][:, los], fmt='%.18e')
        for initial in (True, False):
            mode = 'noqso' if initial else 'qso'
            command = [str(rte), '0' if initial else f'{args.tq:g}',
                       '1' if initial else f'{args.tq:g}', '1', '0', f'{mag:.17g}',
                       str(z), str(label), str(inp), f'los{los}', str(work),
                       str(n), '1', str(los), '1', str(extra), '1', str(box), '1', str(int(initial))]
            filename = (f'noqso_z{label:.2f}LOS{los}.txt' if initial else
                        f'z{label:.2f}LOS{los}tq{args.tq:g}mag-27.09_tep_{args.tq:g}Myr1fd.txt')
            jobs.append(dict(input_los=los, identical_input_ids=group, mode=mode,
                             input_sha256=h, configuration_sha256=config_id,
                             command=command, output=str(work/filename),
                             log=str(work/f'{mode}.log'), completion=str(work/f'{mode}.json')))
    write_json(out/'commands.json', jobs)
    print(f'Prepared {len(ids)} rows, {len(reps)} distinct inputs, {args.workers} workers', flush=True)
    if args.prepare_only:
        return

    def run(job):
        completion, output = Path(job['completion']), Path(job['output'])
        if completion.exists() and output.exists():
            saved = json.loads(completion.read_text())
            if (saved['configuration_sha256'] == config_id and saved['command'] == job['command']
                    and saved['output_sha256'] == digest(output)):
                validate_array(np.loadtxt(output), n)
                return saved
        start = time.monotonic()
        with open(job['log'], 'w') as log:
            subprocess.run(job['command'], check=True, stdout=log, stderr=subprocess.STDOUT)
        a = np.loadtxt(output)
        validate_array(a, n)
        los = job['input_los']
        check = {'helium_fraction_sum_max_absolute_error': float(np.max(abs(a[:, 3:6].sum(axis=1)-1)))}
        if job['mode'] == 'noqso':
            if check['helium_fraction_sum_max_absolute_error'] > 1e-10:
                raise ValueError(f'Initial helium fractions do not sum to one for LOS {los}')
            check['max_absolute_initial_xHI_error'] = float(np.max(abs(a[:, 1]-grids['nh1'][:, los]/grids['nh'][:, los])))
            check['max_relative_initial_temperature_error'] = float(np.max(abs(a[:, 2]/grids['T'][:, los]-1)))
            if check['max_absolute_initial_xHI_error'] > 1e-12 or check['max_relative_initial_temperature_error'] > 1e-12:
                raise ValueError(f'Initial fields changed for LOS {los}: {check}')
        # Independent synthesis checks every retained spectral bin.
        spectral_error = 0.
        for start_bin in range(0, n, 128):
            bins = np.arange(start_bin, min(n, start_bin+128))
            tau = optical_depth(a, grids['rho'][:, los], grids['vpec'][:, los], z, bins, box, extra)
            spectral_error = max(spectral_error, float(np.max(abs(np.exp(-tau)-np.exp(-a[bins, 6])))))
        check['independent_spectral_max_flux_error'] = spectral_error
        if spectral_error > 1e-6:
            raise ValueError(f'Spectral mismatch for LOS {los}: {spectral_error}')
        result = dict(job, elapsed_seconds=time.monotonic()-start, validation=check,
                      output_sha256=digest(output))
        write_json(completion, result)
        return result

    results = {}
    for mode in ('noqso', 'qso'):
        stage = [job for job in jobs if job['mode'] == mode]
        started = time.monotonic()
        completed = []
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = [pool.submit(run, job) for job in stage]
            try:
                for future in as_completed(futures):
                    completed.append(future.result())
                    write_json(out/f'{mode}_progress.json', dict(
                        completed=len(completed), total=len(stage), elapsed_seconds=time.monotonic()-started))
                    if len(completed) % 10 == 0 or len(completed) == len(stage):
                        print(f'{mode}: {len(completed)}/{len(stage)}, {time.monotonic()-started:.1f} s', flush=True)
            except BaseException:
                for future in futures:
                    future.cancel()
                raise
        by_los = {result['input_los']: np.loadtxt(result['output']) for result in completed}
        a = np.stack([by_los[int(los)] for los in row_reps])
        name = (f'noqso_z{label:.2f}_los.npy' if mode == 'noqso' else
                f'qso_z{label:.2f}_m-27.09_tq{args.tq:.0e}_los.npy')
        temporary = out/(name + '.partial')
        with temporary.open('wb') as stream:
            np.save(stream, a)
        temporary.replace(out/name)
        validation = dict(
            shape=list(a.shape), all_finite=bool(np.isfinite(a).all()),
            unique_inputs=len(stage), represented_rows=len(ids),
            independent_spectral_max_flux_error=max(r['validation']['independent_spectral_max_flux_error'] for r in completed),
            helium_fraction_sum_max_absolute_error=float(np.max(abs(a[:, :, 3:6].sum(axis=2)-1))),
        )
        if mode == 'noqso':
            for key in ('max_absolute_initial_xHI_error', 'max_relative_initial_temperature_error'):
                validation[key] = max(r['validation'][key] for r in completed)
        results[mode] = dict(path=str(out/name), sha256=digest(out/name), validation=validation,
                             elapsed_seconds=time.monotonic()-started)
        write_json(out/f'{mode}_results.json', dict(summary=results[mode], runs=completed))
        print(f'Wrote {name}: {a.shape}', flush=True)
    write_json(out/'results.json', results)


if __name__ == '__main__':
    main()
