#!/usr/bin/env python3
"""Infer the RT redshift of an archived LOS from its stored fields and optical depth.

Recompute src/lya_spectra.c using the target's final xHI and temperature, plus
the original input rho and vpec. This isolates spectral redshift from the
chemistry solver. Requires NumPy and SciPy; plots also require Matplotlib.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares


def optical_depth(target, rho, velocity, redshift, bins, boxsize, extra):
    """Port the existing C spectrum formula, including its unit conventions."""
    n_cells = len(target)
    scale = 1.0 / (1.0 + redshift)
    hubble_h = 0.678
    unit_length = 3.085678e21
    b = np.sqrt(2.0 * 1.3806e-16 * target[:, 2] / 1.6726e-24)
    damping = 6.265e8 * 1215.6e-8 / (4.0 * np.pi) / b
    nhi = rho * target[:, 1] * 1.989e43 / 1.6726e-24 / hubble_h
    wmax = (
        scale * 3.2407789e-18 * unit_length / 1e5
        * np.sqrt(0.308 / scale**3 + 0.692) * boxsize / hubble_h
    )
    dv = (
        (np.arange(n_cells)[None, :] + extra - bins[:, None]) * wmax / n_cells
        + velocity
    ) * 1e5
    t0 = np.maximum(dv**2 / b**2, 1e-20)
    t1 = np.exp(-t0)
    t2 = 1.5 / t0
    profile = np.where(
        t0 < 1e-6,
        t1,
        t1 - damping / np.sqrt(np.pi) / t0
        * (t1**2 * (4.0 * t0**2 + 7.0 * t0 + 4.0 + t2) - t2 - 1.0),
    )
    factor = (
        hubble_h**2 / unit_length**2 * 0.41615 * np.pi * 1215.6e-8
        * np.sqrt(3.0 * 6.65245e-25 / (8.0 * np.pi))
        * 2.9979e10 / scale**2 / np.sqrt(np.pi)
    )
    return np.sum(nhi * profile / b * boxsize / n_cells, axis=1) * factor


def comparison_metrics(target, candidate):
    if candidate.shape != target.shape or not np.all(np.isfinite(candidate)):
        raise ValueError("comparison must have the target shape and finite values")
    delta = np.abs(candidate - target)
    relative = delta / np.maximum(np.abs(target), 1e-40)
    flux_delta = np.abs(np.exp(-candidate[:, 6]) - np.exp(-target[:, 6]))
    return {
        "columns": ["distance_ckpc", "xHI", "T", "xHeI", "xHeII", "xHeIII", "tau"],
        "median_relative_error": np.median(relative, axis=0).tolist(),
        "max_relative_error": np.max(relative, axis=0).tolist(),
        "max_absolute_error": np.max(delta, axis=0).tolist(),
        "mean_absolute_flux_error": float(np.mean(flux_delta)),
        "max_absolute_flux_error": float(np.max(flux_delta)),
        "exactly_equal": bool(np.array_equal(candidate, target)),
    }


def plot_comparisons(target, comparisons, extra, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 1, figsize=(10, 9), constrained_layout=True)
    distance = target[:, 0] / 1000.0
    spectral_distance = distance - extra * (distance[1] - distance[0])
    series = [("Reference", target, {"color": "black", "lw": 1.8})]
    series.extend((label, data, {"lw": 1.0, "ls": "--"}) for label, data in comparisons)
    for label, data, style in series:
        axes[0].semilogy(distance, data[:, 1], label=label, **style)
        axes[1].plot(distance, data[:, 2] / 1e4, **style)
        axes[2].plot(spectral_distance, np.exp(-data[:, 6]), **style)
    axes[0].set_ylabel("Neutral hydrogen fraction")
    axes[0].legend(fontsize=9)
    axes[1].set_ylabel("Temperature (10⁴ K)")
    axes[2].set_ylabel("Transmitted flux")
    for ax in axes[:2]:
        ax.set_xlabel("Gas distance (comoving Mpc)")
    axes[2].set_xlabel("Spectral coordinate (comoving Mpc; 256-pixel offset removed)"
                       if extra == 256 else "Spectral coordinate (comoving Mpc)")
    for ax in axes:
        ax.grid(alpha=0.2)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", type=Path)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--input-prefix", required=True)
    parser.add_argument("--input-los", type=int, required=True,
                        help="Column within the chosen f1k/l1k input file")
    parser.add_argument("--target-row", type=int, default=0)
    parser.add_argument("--redshifts", type=float, nargs="+", required=True)
    parser.add_argument("--boxsize", type=float, default=102441.268)
    parser.add_argument("--num-pixels-redward", type=int, default=256)
    parser.add_argument("--bin-stride", type=int, default=3)
    parser.add_argument("--compare", nargs=2, action="append", default=[],
                        metavar=("LABEL", "TEXT_OUTPUT"))
    parser.add_argument("--report", type=Path)
    parser.add_argument("--plot", type=Path)
    args = parser.parse_args()
    if args.bin_stride < 1:
        parser.error("--bin-stride must be positive")

    target = np.load(args.target, mmap_mode="r")[args.target_row]
    if target.ndim != 2 or target.shape[1] != 7 or not np.all(np.isfinite(target)):
        raise ValueError("expected a finite target row with seven columns")
    grids = [
        np.loadtxt(args.input_dir / f"{args.input_prefix}_{suffix}.txt",
                   usecols=args.input_los, max_rows=len(target))
        for suffix in ("rho", "vpec")
    ]
    bins = np.arange(0, len(target), args.bin_stride)
    observed = target[bins, 6]
    if np.any(observed <= 0):
        raise ValueError("redshift fitting requires positive optical depths")

    def residual(redshift):
        predicted = optical_depth(target, *grids, redshift, bins,
                                  args.boxsize, args.num_pixels_redward)
        return np.log(predicted / observed)

    trials = [{"redshift": z, "rms_log_tau_error": float(np.sqrt(np.mean(residual(z)**2)))}
              for z in args.redshifts]
    best = min(trials, key=lambda item: item["rms_log_tau_error"])
    fit = least_squares(lambda p: residual(p[0]), [best["redshift"]],
                        bounds=([0.0], [20.0]), xtol=1e-13, ftol=1e-13, gtol=1e-13)
    report = {
        "target": str(args.target.resolve()), "target_row": args.target_row,
        "input_prefix": args.input_prefix, "input_los": args.input_los,
        "candidate_redshifts": trials, "fitted_redshift": float(fit.x[0]),
        "fit_success": bool(fit.success),
        "rms_log_tau_error": float(np.sqrt(np.mean(fit.fun**2))),
        "comparisons": {},
    }
    comparisons = []
    for label, path in args.compare:
        data = np.loadtxt(path)
        report["comparisons"][label] = comparison_metrics(target, data)
        comparisons.append((label, data))
    text = json.dumps(report, indent=2, allow_nan=False) + "\n"
    print(text, end="")
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(text)
    if args.plot:
        plot_comparisons(target, comparisons, args.num_pixels_redward, args.plot)


if __name__ == "__main__":
    main()
