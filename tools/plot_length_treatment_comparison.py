#!/usr/bin/env python3
"""Compare fixed-comoving and literal Section 3.1 length treatments."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np


OMEGA_M = 0.308
H = 0.678


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixed", type=Path, required=True)
    parser.add_argument("--literal", type=Path, required=True)
    parser.add_argument("--fixed-config", type=Path, required=True)
    parser.add_argument("--literal-config", type=Path, required=True)
    parser.add_argument("--z-sim", type=float, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def max_abs(values: np.ndarray) -> float:
    return float(np.max(np.abs(values)))


def main() -> None:
    args = parse_args()
    fixed_config = json.loads(args.fixed_config.read_text())
    literal_config = json.loads(args.literal_config.read_text())

    fixed_without_length = dict(fixed_config)
    literal_without_length = dict(literal_config)
    fixed_scale = fixed_without_length.pop("comoving_length_scale")
    literal_scale = literal_without_length.pop("comoving_length_scale")
    if fixed_without_length != literal_without_length:
        raise SystemExit("run configurations differ in more than comoving_length_scale")

    fixed = np.loadtxt(args.fixed)
    literal = np.loadtxt(args.literal)
    if fixed.shape != literal.shape or fixed.ndim != 2 or fixed.shape[1] != 7:
        raise SystemExit(f"unexpected output shapes: {fixed.shape}, {literal.shape}")

    z_qso = float(fixed_config["quasar_redshift"])
    origin = int(fixed_config["num_pixels_redward"])
    expansion = (1.0 + z_qso) / (1.0 + args.z_sim)
    expected_literal_scale = expansion**2
    if not np.isclose(fixed_scale, 1.0, rtol=0.0, atol=1e-14):
        raise SystemExit(f"fixed-comoving scale is {fixed_scale}, expected 1")
    if not np.isclose(literal_scale, expected_literal_scale, rtol=1e-13):
        raise SystemExit(
            f"literal scale is {literal_scale}, expected A^2={expected_literal_scale}"
        )

    fixed_dx_comoving_kpc = float(np.median(np.diff(fixed[:, 0])))
    literal_dx_comoving_kpc = float(np.median(np.diff(literal[:, 0])))
    fixed_dx_proper_kpc = fixed_dx_comoving_kpc / (1.0 + z_qso)
    literal_dx_proper_kpc = literal_dx_comoving_kpc / (1.0 + z_qso)
    snapshot_dx_proper_kpc = fixed_dx_comoving_kpc / (1.0 + args.z_sim)

    hubble = 100.0 * H * np.sqrt(OMEGA_M * (1.0 + z_qso) ** 3 + 1.0 - OMEGA_M)
    fixed_r_pmpc = (fixed[:, 0] - fixed[origin, 0]) / (1000.0 * (1.0 + z_qso))
    literal_r_pmpc = (literal[:, 0] - literal[origin, 0]) / (1000.0 * (1.0 + z_qso))
    fixed_velocity = hubble * fixed_r_pmpc
    literal_velocity = hubble * literal_r_pmpc
    fixed_flux = np.exp(-fixed[:, 6])
    literal_flux = np.exp(-literal[:, 6])

    common_velocity = np.arange(
        max(fixed_velocity[0], literal_velocity[0]),
        min(fixed_velocity[-1], literal_velocity[-1]),
        2.0,
    )
    fixed_common = np.interp(common_velocity, fixed_velocity, fixed_flux)
    literal_common = np.interp(common_velocity, literal_velocity, literal_flux)
    common_delta = literal_common - fixed_common
    plotted_velocity = (common_velocity >= -500.0) & (common_velocity <= 5000.0)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    report = {
        "controlled_difference": "comoving_length_scale only",
        "z_sim": args.z_sim,
        "z_qso": z_qso,
        "A": expansion,
        "density_multiplier_in_both_runs": fixed_config["density_multiplier"],
        "fixed_comoving_scale": fixed_scale,
        "paper_literal_comoving_scale_A_squared": literal_scale,
        "pixel_count_in_both_runs": int(fixed.shape[0]),
        "fixed_comoving_pixel_kpc": fixed_dx_comoving_kpc,
        "paper_literal_comoving_pixel_kpc": literal_dx_comoving_kpc,
        "snapshot_proper_pixel_kpc": snapshot_dx_proper_kpc,
        "fixed_proper_pixel_at_zqso_kpc": fixed_dx_proper_kpc,
        "paper_literal_proper_pixel_at_zqso_kpc": literal_dx_proper_kpc,
        "fixed_proper_ratio_to_snapshot": fixed_dx_proper_kpc / snapshot_dx_proper_kpc,
        "paper_literal_proper_ratio_to_snapshot": literal_dx_proper_kpc
        / snapshot_dx_proper_kpc,
        "hubble_zqso_km_s_Mpc": hubble,
        "fixed_velocity_pixel_km_s": float(np.median(np.diff(fixed_velocity))),
        "paper_literal_velocity_pixel_km_s": float(
            np.median(np.diff(literal_velocity))
        ),
        "same_pixel_flux_mean_absolute_difference": float(
            np.mean(np.abs(literal_flux - fixed_flux))
        ),
        "same_pixel_flux_max_absolute_difference": max_abs(literal_flux - fixed_flux),
        "common_velocity_flux_mean_absolute_difference": float(
            np.mean(np.abs(common_delta))
        ),
        "common_velocity_flux_max_absolute_difference": max_abs(common_delta),
        "minus500_to_5000_km_s_flux_mean_absolute_difference": float(
            np.mean(np.abs(common_delta[plotted_velocity]))
        ),
        "minus500_to_5000_km_s_flux_max_absolute_difference": max_abs(
            common_delta[plotted_velocity]
        ),
        "final_xHI_mean_absolute_difference": float(
            np.mean(np.abs(literal[:, 1] - fixed[:, 1]))
        ),
        "final_xHI_max_absolute_difference": max_abs(literal[:, 1] - fixed[:, 1]),
        "final_temperature_mean_absolute_difference_K": float(
            np.mean(np.abs(literal[:, 2] - fixed[:, 2]))
        ),
        "final_temperature_max_absolute_difference_K": max_abs(
            literal[:, 2] - fixed[:, 2]
        ),
        "fixed_output": str(args.fixed.resolve()),
        "paper_literal_output": str(args.literal.resolve()),
    }
    (args.output_dir / "length_treatment_comparison.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )

    os.environ.setdefault("MPLCONFIGDIR", str(args.output_dir / ".matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(args.output_dir / ".cache"))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    pixel = np.arange(fixed.shape[0]) - origin
    colors = {"fixed": "#1f77b4", "literal": "#d95f02"}
    fig, axes = plt.subplots(2, 2, figsize=(13, 8.5), constrained_layout=True)

    axes[0, 0].plot(pixel, fixed_flux, color=colors["fixed"], lw=1.0, label="Fixed comoving grid")
    axes[0, 0].plot(pixel, literal_flux, color=colors["literal"], lw=1.0, label="Literal paper length factor")
    axes[0, 0].set(xlabel="Pixel offset from Lyα origin", ylabel="Transmitted flux", xlim=(-40, 430), ylim=(-0.03, 1.03))
    axes[0, 0].legend(fontsize=9)

    axes[0, 1].plot(fixed_velocity, fixed_flux, color=colors["fixed"], lw=1.0, label=f"Fixed: {report['fixed_velocity_pixel_km_s']:.2f} km s$^{{-1}}$/pixel")
    axes[0, 1].plot(literal_velocity, literal_flux, color=colors["literal"], lw=1.0, label=f"Literal: {report['paper_literal_velocity_pixel_km_s']:.2f} km s$^{{-1}}$/pixel")
    axes[0, 1].set(xlabel="Hubble velocity from Lyα origin (km s$^{-1}$)", ylabel="Transmitted flux", xlim=(-500, 5000), ylim=(-0.03, 1.03))
    axes[0, 1].legend(fontsize=9)

    axes[1, 0].plot(common_velocity[plotted_velocity], common_delta[plotted_velocity], color="#6a3d9a", lw=1.0)
    axes[1, 0].axhline(0.0, color="black", lw=0.7)
    axes[1, 0].set(xlabel="Common Hubble-velocity grid (km s$^{-1}$)", ylabel="Flux: literal − fixed", xlim=(-500, 5000))

    offsets = np.arange(-40, 431)
    axes[1, 1].plot(offsets, offsets * fixed_dx_proper_kpc / 1000.0, color=colors["fixed"], lw=2.0, label=f"Fixed: {fixed_dx_proper_kpc:.2f} proper kpc/pixel")
    axes[1, 1].plot(offsets, offsets * literal_dx_proper_kpc / 1000.0, color=colors["literal"], lw=2.0, label=f"Literal: {literal_dx_proper_kpc:.2f} proper kpc/pixel")
    axes[1, 1].set(xlabel="Pixel offset from Lyα origin", ylabel="Proper distance at $z_{QSO}$ (pMpc)", xlim=(-40, 430))
    axes[1, 1].legend(fontsize=9)

    for ax in axes.flat:
        ax.grid(alpha=0.18)
    fig.suptitle(
        f"Controlled length-treatment test: $z_{{sim}}={args.z_sim:.6f}$, "
        f"$z_{{QSO}}={z_qso:.4f}$, $M_{{1450}}={fixed_config['m1450']:.2f}$, "
        "$t_Q=1$ Myr\n"
        "Same input skewer, density scaling, background, cosmology, and Lyα units"
    )
    for suffix in ("png", "pdf"):
        fig.savefig(args.output_dir / f"length_treatment_comparison.{suffix}", dpi=180)
    plt.close(fig)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
