#!/usr/bin/env python3
"""Bin the frozen archived-versus-literal Section 3.1 comparison."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--bin-width", type=float, default=600.0)
    parser.add_argument("--velocity-min", type=float, default=0.0)
    parser.add_argument("--velocity-max", type=float, default=4800.0)
    return parser.parse_args()


def average_in_bins(
    velocity: np.ndarray, flux: np.ndarray, edges: np.ndarray
) -> np.ndarray:
    """Average a spectrum using piecewise-linear integration on its native grid."""
    if edges[0] < velocity[0] or edges[-1] > velocity[-1]:
        raise ValueError("Requested bins extend beyond the native velocity grid")
    result = []
    for left, right in zip(edges[:-1], edges[1:]):
        inside = (velocity > left) & (velocity < right)
        nodes = np.r_[left, velocity[inside], right]
        values = np.interp(nodes, velocity, flux)
        integral = np.sum(np.diff(nodes) * (values[:-1] + values[1:]) / 2.0)
        result.append(integral / (right - left))
    return np.asarray(result)


def main() -> None:
    args = parse_args()
    if args.bin_width <= 0.0 or args.velocity_max <= args.velocity_min:
        raise SystemExit("Invalid velocity interval or bin width")
    span = args.velocity_max - args.velocity_min
    number_of_bins = round(span / args.bin_width)
    if not np.isclose(number_of_bins * args.bin_width, span):
        raise SystemExit("Velocity interval must contain an integer number of bins")

    with np.load(args.input) as saved:
        reference = saved["reference"][0]
        reproduced = saved["reproduced"][0]
        literal = saved["literal_section31_row0"]
        reference_velocity = saved["reference_velocity_kms"]
        literal_velocity = saved["literal_velocity_kms"]

    reference_flux = np.exp(-reference[:, 6])
    reproduced_flux = np.exp(-reproduced[:, 6])
    literal_flux = np.exp(-literal[:, 6])
    reproduction_native_error = float(
        np.max(np.abs(reproduced_flux - reference_flux))
    )
    if reproduction_native_error > 1e-7:
        raise SystemExit(
            "Frozen matched reproduction is inconsistent with the archived profile: "
            f"maximum flux difference {reproduction_native_error}"
        )

    edges = np.linspace(
        args.velocity_min, args.velocity_max, number_of_bins + 1
    )
    centers = (edges[:-1] + edges[1:]) / 2.0
    archived_binned = average_in_bins(reference_velocity, reference_flux, edges)
    reproduced_binned = average_in_bins(reference_velocity, reproduced_flux, edges)
    literal_binned = average_in_bins(literal_velocity, literal_flux, edges)
    reproduction_delta = reproduced_binned - archived_binned
    literal_delta = literal_binned - archived_binned

    report = {
        "source_npz": str(args.input.resolve()),
        "source_npz_sha256": hashlib.sha256(args.input.read_bytes()).hexdigest(),
        "source_row": 0,
        "scope": "One archived sightline; each curve integrated on its native velocity grid.",
        "bin_width_km_s": args.bin_width,
        "bin_edges_km_s": edges.tolist(),
        "bin_centers_km_s": centers.tolist(),
        "archived_flux": archived_binned.tolist(),
        "matched_reproduction_flux": reproduced_binned.tolist(),
        "literal_section31_flux": literal_binned.tolist(),
        "literal_minus_archived": literal_delta.tolist(),
        "matched_reproduction_max_native_flux_difference": reproduction_native_error,
        "matched_reproduction_max_binned_flux_difference": float(
            np.max(np.abs(reproduction_delta))
        ),
        "literal_mean_absolute_binned_flux_difference": float(
            np.mean(np.abs(literal_delta))
        ),
        "literal_max_absolute_binned_flux_difference": float(
            np.max(np.abs(literal_delta))
        ),
        "literal_mean_signed_binned_flux_difference": float(np.mean(literal_delta)),
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "literal_section31_comparison_600kms.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )

    os.environ.setdefault("MPLCONFIGDIR", str(args.output_dir / ".matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(args.output_dir / ".cache"))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(
        2,
        1,
        figsize=(9, 7),
        sharex=True,
        gridspec_kw={"height_ratios": [2, 1]},
        constrained_layout=True,
    )
    axes[0].plot(
        centers, archived_binned, "o-", color="black", lw=1.5,
        label="Archived rescaledrho file",
    )
    axes[0].plot(
        centers, reproduced_binned, "s--", color="#0072b2", lw=1.2,
        label="Matched reproduction",
    )
    axes[0].plot(
        centers, literal_binned, "o-", color="#d55e00", lw=1.5,
        label="Literal Section 3.1",
    )
    axes[0].set(ylabel="Mean transmitted flux", ylim=(-0.03, 1.03))
    axes[0].legend(fontsize=9)

    axes[1].plot(
        centers, reproduction_delta, "s--", color="#0072b2", lw=1.2,
        label="Reproduction − archived",
    )
    axes[1].plot(
        centers, literal_delta, "o-", color="#d55e00", lw=1.5,
        label="Literal − archived",
    )
    axes[1].axhline(0.0, color="black", lw=0.8)
    axes[1].set(
        xlabel="Velocity blueward of Lyα (km s$^{-1}$)",
        ylabel="Flux difference",
        xlim=(args.velocity_min, args.velocity_max),
    )
    axes[1].legend(fontsize=9)
    for axis in axes:
        axis.grid(alpha=0.2)
    fig.suptitle(
        "Same sightline after 600 km s$^{-1}$ binning\n"
        "$z_{sim}=5.755954$, $z_{QSO}=7.0848$, "
        "$M_{1450}=-26.58$, $t_Q=1$ Myr"
    )
    for suffix in ("png", "pdf"):
        fig.savefig(
            args.output_dir / f"literal_section31_comparison_600kms.{suffix}",
            dpi=180,
        )
    plt.close(fig)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
