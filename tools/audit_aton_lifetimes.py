#!/usr/bin/env python3
"""Check every saved Ly-alpha bin of the z6.15 and z5.95 lifetime grids.

Uses target final temperature/xHI and original rho/velocity to isolate the
redshift in the spectrum calculation. Kept-index lists map archive rows to
global f1k/l1k input columns and must be supplied explicitly.
"""

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
    parser.add_argument("--kept-z615", type=Path, required=True)
    parser.add_argument("--kept-z595", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = {"method": "Every archive row and all 891 saved bins; target final fields held fixed",
              "m1450": -27.09, "snapshots": {}}
    for label, expected, matched, kept_path in [
        ("6.15", 6.148904, 6.604205, args.kept_z615),
        ("5.95", 5.945518, 5.945518, args.kept_z595),
    ]:
        kept = np.loadtxt(kept_path, dtype=int)
        grids = [np.hstack([
            np.loadtxt(args.input_dir / f"{part}1kz{label}M3LOS_{suffix}.txt", max_rows=891)
            for part in ("f", "l")
        ]) for suffix in ("rho", "vpec")]
        bins = np.arange(891)
        results = []
        for tq in (0.01, 0.1, 1, 10, 100):
            path = args.archive_dir / f"qso_z{label}_m-27.09_tq{tq:.0e}_los.npy"
            data = np.load(path, mmap_mode="r")
            if data.shape != (len(kept), 891, 7):
                raise ValueError(f"unexpected shape {data.shape} for {path}")
            if not np.all(np.isfinite(data)) or np.any(data[:, :, 6] <= 0):
                raise ValueError(f"nonfinite fields or nonpositive optical depths in {path}")

            def residual(z):
                tau = optical_depth(data[0], grids[0][:, kept[0]], grids[1][:, kept[0]],
                                    z, bins, 102441.268, 256)
                return np.log(tau / data[0, :, 6])

            fit = least_squares(lambda p: residual(p[0]), [matched],
                                bounds=([0], [20]), xtol=1e-13, ftol=1e-13, gtol=1e-13)
            errors = []
            for row, los in enumerate(kept):
                tau = optical_depth(data[row], grids[0][:, los], grids[1][:, los],
                                    matched, bins, 102441.268, 256)
                errors.append(float(np.max(np.abs(tau / data[row, :, 6] - 1))))
            worst_row = int(np.argmax(errors))
            item = {
                "tq_Myr": tq, "file": str(path.resolve()), "rows": len(data),
                "bins_per_row": len(bins), "input_redshift": expected,
                "checked_spectrum_redshift": matched,
                "row0_fitted_redshift": float(fit.x[0]), "fit_success": bool(fit.success),
                "row0_rms_log_tau_at_input_redshift": float(np.sqrt(np.mean(residual(expected)**2))),
                "row0_rms_log_tau_at_6p604205": float(np.sqrt(np.mean(residual(6.604205)**2))),
                "all_rows_max_relative_tau_error": max(errors),
                "rows_above_1e_minus_6": int(np.sum(np.asarray(errors) > 1e-6)),
                "worst_row": worst_row, "worst_input_los": int(kept[worst_row]),
            }
            results.append(item)
            print(f"z{label} tq={tq:g} Myr: fitted z={fit.x[0]:.12f}, "
                  f"rows={len(data)}, bins={len(bins)}, max relative tau error={max(errors):.3e}, "
                  f"rows above 1e-6={item['rows_above_1e_minus_6']}", flush=True)
            report["snapshots"][label] = results
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
