#!/usr/bin/env python3
"""Audit the native z=7.14, M1450=-26.58 density and archived profile."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np


H = 0.678
OMEGA_B = 0.0478
OMEGA_M = 0.308
HYDROGEN_MASS_FRACTION = 0.76
G_SI = 6.67430e-11
PROTON_MASS_KG = 1.67262192369e-27
Z_SIM = 7.139393
Z_QSO = 7.0848
N_CELLS = 2048
N_LOS = 1000
EXTRA = 256
RECORD_BYTES = 114_760
HEADER_BYTES = 64
BLOCK_BYTES = N_CELLS * 8


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dom-file", type=Path, required=True)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--qso-redshift-run", type=Path, required=True)
    parser.add_argument("--archive-compatible-run", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def cosmic_comoving_hydrogen_density_m3() -> float:
    hubble_si = 100.0 * H * 1000.0 / 3.085677581491367e22
    critical_density = 3.0 * hubble_si**2 / (8.0 * np.pi * G_SI)
    return (
        HYDROGEN_MASS_FRACTION
        * OMEGA_B
        * critical_density
        / PROTON_MASS_KG
    )


def halo_centered_mean_density(dom_file: Path) -> np.ndarray:
    data = np.memmap(dom_file, dtype="<f8", mode="r")
    expected_values = N_LOS * RECORD_BYTES // 8
    if data.size != expected_values:
        raise ValueError(f"Unexpected binary size: {data.size}, expected {expected_values}")
    mean_density = np.zeros(N_CELLS)
    for sightline in range(N_LOS):
        record_start = sightline * RECORD_BYTES // 8
        halo_coordinate = float(data[record_start + 1])
        offset = int(halo_coordinate * (160.0 - 0.078125) / 0.078125)
        density_start = (
            sightline * RECORD_BYTES + HEADER_BYTES + 3 * BLOCK_BYTES
        ) // 8
        density = np.asarray(data[density_start : density_start + N_CELLS])
        mean_density += np.r_[density[offset:], density[:offset]]
    return mean_density / N_LOS


def legacy_velocity_grid(redshift: float) -> np.ndarray:
    """Velocity grid used by the archived calculation."""
    scale_factor = 1.0 / (1.0 + redshift)
    expansion = np.sqrt(OMEGA_M / scale_factor**3 + 1.0 - OMEGA_M)
    box_comoving_kpc = 235294.118
    width_km_s = (
        scale_factor * (100.0 / H) * expansion * box_comoving_kpc / 1000.0
    )
    return (np.arange(N_CELLS) - EXTRA) * width_km_s / N_CELLS


def bin_flux(velocity: np.ndarray, flux: np.ndarray, edges: np.ndarray) -> np.ndarray:
    result = []
    for left, right in zip(edges[:-1], edges[1:]):
        inside = (velocity > left) & (velocity < right)
        nodes = np.r_[left, velocity[inside], right]
        values = np.interp(nodes, velocity, flux)
        integral = np.sum(np.diff(nodes) * (values[:-1] + values[1:]) / 2.0)
        result.append(integral / (right - left))
    return np.asarray(result)


def difference_metrics(candidate: np.ndarray, reference: np.ndarray) -> dict[str, float]:
    difference = candidate - reference
    return {
        "mean_absolute": float(np.mean(np.abs(difference))),
        "maximum_absolute": float(np.max(np.abs(difference))),
        "mean_signed": float(np.mean(difference)),
    }


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    density = halo_centered_mean_density(args.dom_file)
    cosmic_hydrogen = cosmic_comoving_hydrogen_density_m3()
    corrected_ratio = density / cosmic_hydrogen
    legacy_ratio = HYDROGEN_MASS_FRACTION * corrected_ratio
    proper_pixel_mpc = (78.125 / H) / (1.0 + Z_SIM) / 1000.0
    distance_pmpc = np.arange(N_CELLS) * proper_pixel_mpc

    density_regions = {}
    for left, right in ((0, 1), (1, 5), (5, 10), (10, 20), (20, 28)):
        selected = (distance_pmpc >= left) & (distance_pmpc < right)
        density_regions[f"{left}_{right}_proper_Mpc"] = {
            "cube_D_as_nH_over_cosmic_nH": float(np.mean(corrected_ratio[selected])),
            "legacy_times_0p76_over_cosmic_nH": float(np.mean(legacy_ratio[selected])),
        }

    archived = np.load(args.archive, mmap_mode="r")[0]
    qso_redshift_run = np.loadtxt(args.qso_redshift_run)
    archive_compatible = np.loadtxt(args.archive_compatible_run)
    for name, values in {
        "archive": archived,
        "qso_redshift_run": qso_redshift_run,
        "archive_compatible_run": archive_compatible,
    }.items():
        if values.shape != (N_CELLS, 7):
            raise ValueError(f"Unexpected {name} shape: {values.shape}")

    archive_velocity = legacy_velocity_grid(Z_SIM)
    qso_velocity = legacy_velocity_grid(Z_QSO)
    edges = np.arange(0.0, 4800.0 + 600.0, 600.0)
    centers = (edges[:-1] + edges[1:]) / 2.0
    archive_binned = bin_flux(archive_velocity, np.exp(-archived[:, 6]), edges)
    qso_binned = bin_flux(qso_velocity, np.exp(-qso_redshift_run[:, 6]), edges)
    compatible_binned = bin_flux(
        archive_velocity, np.exp(-archive_compatible[:, 6]), edges
    )

    report = {
        "quasar": {
            "name": "J1120+0641",
            "M1450": -26.58,
            "observed_redshift": Z_QSO,
            "native_snapshot_redshift": Z_SIM,
        },
        "density_normalization": {
            "omega_b": OMEGA_B,
            "h": H,
            "hydrogen_mass_fraction": HYDROGEN_MASS_FRACTION,
            "cosmic_comoving_nH_m_minus_3": cosmic_hydrogen,
            "regions": density_regions,
            "interpretation": (
                "Outside the halo environment, cube_D itself is near the cosmic "
                "hydrogen mean; applying another factor 0.76 makes it about 0.76 "
                "times the cosmic hydrogen mean."
            ),
        },
        "profiles": {
            "bin_edges_km_s": edges.tolist(),
            "archive_flux": archive_binned.tolist(),
            "native_at_observed_qso_redshift_flux": qso_binned.tolist(),
            "archive_compatible_snapshot_runtime_flux": compatible_binned.tolist(),
            "native_at_observed_qso_redshift_vs_archive": difference_metrics(
                qso_binned, archive_binned
            ),
            "archive_compatible_vs_archive": difference_metrics(
                compatible_binned, archive_binned
            ),
        },
        "settings": {
            "common": {
                "snapshot": "z7.14 input fields",
                "density_multiplier_relative_to_legacy_text": 1.0 / 0.76,
                "redshift_density_multiplier": 1.0,
                "comoving_length_multiplier": 1.0,
                "M1450_as_float32": float(np.float32(-26.58)),
                "tQ_Myr": 1.0,
                "background": "unchanged snapshot files",
                "legacy_helium_branch": True,
                "legacy_spectral_units": True,
            },
            "native_at_observed_qso_redshift": {
                "runtime_redshift": Z_QSO,
                "purpose": "No explicit A or A^3 remapping; run the quasar at its observed redshift.",
            },
            "archive_compatible": {
                "runtime_redshift": Z_SIM,
                "purpose": "Recovered archived behavior; runtime geometry stays at the snapshot redshift.",
            },
        },
        "files": {
            "dom_file": str(args.dom_file.resolve()),
            "dom_file_sha256": sha256(args.dom_file),
            "archive": str(args.archive.resolve()),
            "archive_sha256": sha256(args.archive),
            "qso_redshift_run": str(args.qso_redshift_run.resolve()),
            "archive_compatible_run": str(args.archive_compatible_run.resolve()),
        },
    }
    (args.output_dir / "native_baseline_audit.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )

    os.environ.setdefault("MPLCONFIGDIR", str(args.output_dir / ".matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(args.output_dir / ".cache"))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 1, figsize=(9, 7), constrained_layout=True)
    axes[0].plot(distance_pmpc, corrected_ratio, color="#0072b2", lw=1.2,
                 label=r"$\mathrm{cube\_D}$ as $n_H$ (current input correction)")
    axes[0].plot(distance_pmpc, legacy_ratio, color="#d55e00", lw=1.2,
                 label=r"$\mathrm{cube\_D}\times0.76$ (legacy text input)")
    axes[0].axhline(1.0, color="black", ls="--", lw=1.0, label="Cosmic mean hydrogen density")
    axes[0].set(xlabel="Proper distance from host halo (pMpc)",
                ylabel=r"Mean $n_H/\bar n_H$", yscale="log", ylim=(0.4, 150))
    axes[0].legend(fontsize=9)
    far = distance_pmpc >= 1.0
    axes[1].plot(distance_pmpc[far], corrected_ratio[far], color="#0072b2", lw=1.2)
    axes[1].plot(distance_pmpc[far], legacy_ratio[far], color="#d55e00", lw=1.2)
    axes[1].axhline(1.0, color="black", ls="--", lw=1.0)
    axes[1].set(xlabel="Proper distance from host halo (pMpc)",
                ylabel=r"Mean $n_H/\bar n_H$", ylim=(0.5, 1.5))
    for axis in axes:
        axis.grid(alpha=0.2)
    fig.suptitle("Hydrogen-density normalization of the native z=7.14 skewers")
    for suffix in ("png", "pdf"):
        fig.savefig(args.output_dir / f"density_normalization.{suffix}", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(
        2, 1, figsize=(9, 7), sharex=True,
        gridspec_kw={"height_ratios": [2, 1]}, constrained_layout=True
    )
    axes[0].plot(centers, archive_binned, "o-", color="black", lw=1.6,
                 label="Archived native file")
    axes[0].plot(centers, qso_binned, "s-", color="#0072b2", lw=1.4,
                 label=r"Native run at observed $z_{QSO}=7.0848$")
    axes[0].plot(centers, compatible_binned, "^-", color="#009e73", lw=1.3,
                 label=r"Archive-compatible run at $z_{sim}=7.139393$")
    axes[0].set(ylabel="Mean transmitted flux", ylim=(-0.03, 1.03))
    axes[0].legend(fontsize=9)
    axes[1].plot(centers, qso_binned - archive_binned, "s-", color="#0072b2", lw=1.4,
                 label=r"Observed-$z_{QSO}$ run − archive")
    axes[1].plot(centers, compatible_binned - archive_binned, "^-", color="#009e73", lw=1.3,
                 label="Archive-compatible run − archive")
    axes[1].axhline(0.0, color="black", lw=0.8)
    axes[1].set(xlabel=r"Velocity blueward of Ly$\alpha$ (km s$^{-1}$)",
                ylabel="Flux difference", xlim=(0, 4800))
    axes[1].legend(fontsize=9)
    for axis in axes:
        axis.grid(alpha=0.2)
    fig.suptitle(
        "Native baseline, 600 km s$^{-1}$ bins: J1120+0641, "
        "$M_{1450}=-26.58$, $t_Q=1$ Myr"
    )
    for suffix in ("png", "pdf"):
        fig.savefig(args.output_dir / f"native_profile_comparison_600kms.{suffix}", dpi=180)
    plt.close(fig)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
