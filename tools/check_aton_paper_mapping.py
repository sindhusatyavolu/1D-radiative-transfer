#!/usr/bin/env python3
"""Test the paper's fixed-quasar-redshift interpretation of the local ATON grid.

Fits spectral redshift and an independent optical-depth amplitude using each
archive's saved final xHI/T and the original input rho/vpec. This is a spectral
configuration check, not a full thermochemistry or provenance validation.
"""

import argparse
import json
import re
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from diagnose_aton_reproduction import optical_depth


SNAPSHOTS = {"5.11": 5.110, "5.76": 5.755954, "5.95": 5.945518,
             "6.15": 6.148904, "6.60": 6.604205, "7.14": 7.139393,
             "8.15": 8.150, "9.55": 9.546}
# Table 1, https://arxiv.org/html/2401.10328v1#S2.T1
QUASARS = {"-27.09": ("J1509-1749", 6.1225),
           "-26.99": ("J1319+0950", 6.1347),
           "-26.58": ("J1120+0641", 7.0848)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive-dir", type=Path, required=True)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--boxsize", type=float, default=102441.268)
    parser.add_argument("--extra", type=int, default=256)
    args = parser.parse_args()
    report = {"paper": "https://arxiv.org/html/2401.10328v1",
              "scope": "First retained row in each file; all saved bins at fitted parameters",
              "results": []}
    cache = {}
    pattern = re.compile(r"qso_z([\d.]+)_m(-[\d.]+)_tq([\de+.-]+)_los.npy")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    for path in sorted(args.archive_dir.glob("qso_*_los.npy")):
        match = pattern.fullmatch(path.name)
        if not match or match[2] not in QUASARS:
            continue
        label, magnitude, lifetime = match.groups()
        target = np.load(path, mmap_mode="r")
        row = target[0]
        n = len(row)
        key = (label, n)
        if key not in cache:
            cache[key] = [np.hstack([
                np.loadtxt(args.input_dir / f"{part}1kz{label}M3LOS_{suffix}.txt", max_rows=n)
                for part in ("f", "l")]) for suffix in ("T", "rho", "vpec")]
        temp, rho, velocity = cache[key]
        # Initial temperature is a candidate selector only. Agreement with the
        # entire reconstructed spectrum must validate the chosen field content.
        score = np.mean((np.log(temp[-150:]) - np.log(row[-150:, 2, None]))**2, axis=0)
        candidates = list(dict.fromkeys([0, *np.argsort(score)[:4].tolist()]))
        bins = np.arange(0, n, 7)
        quasar, qso_z = QUASARS[magnitude]
        nearest = min(SNAPSHOTS.values(), key=lambda z: abs(z - qso_z))
        fits = []
        for los in candidates:
            def residual(p):
                tau = optical_depth(row, rho[:, los], velocity[:, los], p[0],
                                    bins, args.boxsize, args.extra)
                return np.log(tau / row[bins, 6]) + p[1]

            for start in dict.fromkeys([SNAPSHOTS[label], 6.604205, nearest, qso_z]):
                fit = least_squares(residual, [start, 0.0],
                                    bounds=([3, -5], [12, 5]),
                                    xtol=1e-11, ftol=1e-11, gtol=1e-11)
                rms = float(np.sqrt(np.mean(fit.fun**2)))
                fits.append((rms, los, fit))
                if rms < 1e-6:
                    break
            if rms < 1e-6:
                break
        rms, los, fit = min(fits, key=lambda item: item[0])
        full_bins = np.arange(n)
        tau = optical_depth(row, rho[:, los], velocity[:, los], fit.x[0],
                            full_bins, args.boxsize, args.extra) * np.exp(fit.x[1])
        # Let amplitude float even when testing the quasar redshift, to check
        # whether a density normalization alone could explain the discrepancy.
        trials = []
        for z in (qso_z, nearest, SNAPSHOTS[label], 6.604205):
            trial = optical_depth(row, rho[:, los], velocity[:, los], z,
                                  full_bins, args.boxsize, args.extra)
            log_ratio = np.log(trial / row[:, 6])
            trials.append({"redshift": z,
                           "rms_log_tau_unscaled": float(np.sqrt(np.mean(log_ratio**2))),
                           "best_amplitude": float(np.exp(-np.mean(log_ratio))),
                           "rms_log_tau_free_amplitude": float(np.std(log_ratio))})
        item = {"file": str(path.resolve()), "shape": list(target.shape),
                "snapshot_label": label, "snapshot_z": SNAPSHOTS[label],
                "m1450": float(magnitude), "tq_Myr": float(lifetime),
                "quasar": quasar, "observed_quasar_z": qso_z,
                "nearest_available_quasar_z": nearest,
                "required_density_factor_to_observed_z": ((1+qso_z)/(1+SNAPSHOTS[label]))**3,
                "required_density_factor_to_nearest_z": ((1+nearest)/(1+SNAPSHOTS[label]))**3,
                "input_los_candidate": los, "fit_success": bool(fit.success),
                "fitted_spectral_z": float(fit.x[0]),
                "fitted_tau_amplitude": float(np.exp(fit.x[1])),
                "fit_rms_log_tau": rms,
                "all_bins_max_relative_tau_error": float(np.max(abs(tau/row[:, 6]-1))),
                "all_bins_max_absolute_flux_error": float(np.max(abs(np.exp(-tau)-np.exp(-row[:, 6])))),
                "fixed_redshift_trials": trials}
        report["results"].append(item)
        args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
        print(f"{path.name}: LOS={los}, z={fit.x[0]:.9f}, amplitude={np.exp(fit.x[1]):.9f}, "
              f"max rel tau={item['all_bins_max_relative_tau_error']:.3g}", flush=True)


if __name__ == "__main__":
    main()
