#!/usr/bin/env python3
"""Check the spectral redshift of every retained skewer across ATON lifetimes.

Fit the redshift of the first row from its stored final fields, then use
that redshift to reconstruct tau for every row and selected spectral bin.
The indices file maps archive rows to global input IDs (f1k then l1k).
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from diagnose_aton_reproduction import optical_depth


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive-dir", type=Path, required=True)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--indices-file", type=Path, required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--trial-redshifts", type=float, nargs="+", required=True)
    parser.add_argument("--lifetimes", type=float, nargs="+", default=[0.01, 0.1, 1, 10, 100])
    parser.add_argument("--m1450", type=float, default=-27.09)
    parser.add_argument("--bin-stride", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.bin_stride < 1:
        parser.error("--bin-stride must be positive")
    kept = np.loadtxt(args.indices_file, dtype=int, ndmin=1)
    rho, vel = [np.hstack([
        np.loadtxt(args.input_dir / f"{part}1kz{args.label}M3LOS_{suffix}.txt", max_rows=891)
        for part in ["f", "l"]
    ]) for suffix in ["rho", "vpec"]]
    if np.any(kept < 0) or np.any(kept >= rho.shape[1]):
        raise ValueError("indices must be global input IDs")
    bins = np.arange(0, 891, args.bin_stride)
    reports = []
    args.output.parent.mkdir(parents=True, exist_ok=True)
    for tq in args.lifetimes:
        path = args.archive_dir / f"qso_z{args.label}_m{args.m1450:.2f}_tq{tq:.0e}_los.npy"
        target = np.load(path, mmap_mode="r")
        if target.shape != (len(kept), 891, 7):
            raise ValueError(f"unexpected target shape: {path}: {target.shape}")

        def residual(z):
            predicted = optical_depth(target[0], rho[:, kept[0]], vel[:, kept[0]],
                                      z, bins, 102441.268, 256)
            return np.log(predicted / target[0, bins, 6])

        trials = [{"redshift": z, "row0_rms_log_tau_error": float(np.sqrt(np.mean(residual(z)**2)))}
                  for z in args.trial_redshifts]
        best = min(trials, key=lambda trial: trial["row0_rms_log_tau_error"])
        fit = least_squares(lambda p: residual(p[0]), [best["redshift"]],
                            bounds=([0.0], [20.0]), xtol=1e-13, ftol=1e-13, gtol=1e-13)
        if not fit.success:
            raise RuntimeError(f"redshift fit failed: {path}")
        # Use the listed snapshot redshift for the complete reconstruction;
        # retain the continuous fit separately as an independent check.
        z = best["redshift"]
        max_relative = []
        max_flux = []
        for row, los in enumerate(kept):
            predicted = optical_depth(target[row], rho[:, los], vel[:, los],
                                      z, bins, 102441.268, 256)
            observed = target[row, bins, 6]
            if np.any(observed <= 0) or not np.all(np.isfinite(predicted)):
                raise ValueError(f"invalid optical depths in {path}, row {row}")
            max_relative.append(float(np.max(abs(predicted / observed - 1))))
            max_flux.append(float(np.max(abs(np.exp(-predicted) - np.exp(-observed)))))
        worst = int(np.argmax(max_relative))
        result = {
            "label": args.label, "tq_Myr": tq, "path": str(path),
            "shape": list(target.shape), "input_mapping": str(args.indices_file),
            "rows_checked": len(kept), "bins_per_row": len(bins),
            "redshift_trials": trials, "fitted_redshift_row0": float(fit.x[0]),
            "reconstruction_redshift": z,
            "max_relative_tau_error": max(max_relative),
            "median_of_row_max_relative_tau_error": float(np.median(max_relative)),
            "max_absolute_flux_error": max(max_flux),
            "worst_tau_row": worst, "worst_tau_input_los": int(kept[worst]),
        }
        reports.append(result)
        args.output.write_text(json.dumps(reports, indent=2, allow_nan=False) + "\n")
        print(json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
