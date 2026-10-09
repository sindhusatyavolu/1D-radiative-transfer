#!/usr/bin/env python3
"""Run or assemble the eight-snapshot ATON grid from Durovcikova et al. III.1."""

from __future__ import annotations

import argparse
import json
import signal
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class Snapshot:
    label: str
    redshift: float
    paper_xhi: float


# Section III.1 lists eight ATON neutral-fraction models.  The exact redshifts
# below are the redshifts represented by the corresponding local input files.
PAPER_SNAPSHOTS = (
    Snapshot("5.11", 5.109794, 3.7e-5),
    Snapshot("5.76", 5.755954, 0.07),
    Snapshot("5.95", 5.945518, 0.13),
    Snapshot("6.15", 6.148904, 0.21),
    Snapshot("6.60", 6.604205, 0.37),
    Snapshot("7.14", 7.139393, 0.54),
    Snapshot("8.15", 8.150086, 0.75),
    Snapshot("9.55", 9.545964, 0.89),
)

ATON_CELL_COMOVING_KPC = 78.125 / 0.678

FIELDS = (
    "nh",
    "nh1",
    "T",
    "vpec",
    "rho",
    "Gammabg",
    "Gammabghe1",
    "Gammabghe2",
    "Heatbg",
)


def run_child(command: list[str]) -> None:
    """Run a child driver and allow it to clean up its solvers on Ctrl-C."""
    process = subprocess.Popen(command)
    try:
        returncode = process.wait()
    except KeyboardInterrupt:
        if process.poll() is None:
            try:
                process.send_signal(signal.SIGINT)
            except ProcessLookupError:
                pass
        try:
            process.wait(timeout=10)
        except subprocess.TimeoutExpired:
            process.terminate()
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait()
        raise
    if returncode != 0:
        raise SystemExit(f"child command failed with exit status {returncode}")


def add_grid_selection(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--z-qso", type=float, required=True)
    parser.add_argument("--m1450", type=float, required=True)
    parser.add_argument("--tq-grid", type=float, nargs="+", default=[0.01, 0.1, 1, 10, 100])
    parser.add_argument(
        "--snapshots",
        nargs="+",
        choices=[snapshot.label for snapshot in PAPER_SNAPSHOTS],
        default=[snapshot.label for snapshot in PAPER_SNAPSHOTS],
        help="Snapshot labels to process; the default is the eight-point paper grid",
    )
    parser.add_argument("--los", nargs="+", default=["0:999"])
    parser.add_argument("--n-cells", type=int, default=891)
    parser.add_argument("--initial-only", action="store_true")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="action", required=True)

    run = subparsers.add_parser("run", help="Run the text-output sightline grid")
    add_grid_selection(run)
    run.add_argument("--rte", type=Path, default=Path("./rte1d_opt"))
    run.add_argument("--input-dir", type=Path, required=True)
    run.add_argument("--workers", type=int, default=1)
    run.add_argument("--resume", action="store_true")
    run.add_argument("--dry-run", action="store_true")
    run.add_argument("--n-los", type=int, default=500)
    run.add_argument("--split-los-at", type=int, default=500)
    run.add_argument("--legacy-density", type=int, choices=(0, 1), default=0)
    run.add_argument("--num-pixels-redward", type=int, default=256)
    run.add_argument("--legacy-alphahe2-bug", type=int, choices=(0, 1), default=0)
    run.add_argument(
        "--length-treatment", choices=("fixed-comoving", "paper-literal"),
        default="fixed-comoving",
        help=("fixed-comoving uses the simulation's comoving grid at z_qso; paper-literal "
              "implements the physical-length factor printed in Section 3.1"),
    )
    run.add_argument(
        "--background-treatment", choices=("equilibrium", "snapshot"),
        default="equilibrium",
        help="equilibrium rescales the supplied background rates with the corrected density",
    )
    run.add_argument(
        "--spectral-density-treatment", choices=("hydrogen", "legacy-total-gas"),
        default="hydrogen",
        help="hydrogen applies X=0.76 to total-gas rho during Ly-alpha synthesis",
    )
    run.add_argument(
        "--spectral-unit-treatment", choices=("consistent-ckpc", "legacy"),
        default="consistent-ckpc",
        help=("consistent-ckpc uses H0=100h with the stored ckpc grid and Gadget "
              "h^2 density units; legacy retains the old h-factor mismatch"),
    )
    run.add_argument(
        "--boxsize-override", type=float, default=0.0,
        help="0 uses n_cells times cellsize; use 102441.268 only for legacy archive tests",
    )
    run.add_argument(
        "--cellsize-comoving-kpc", type=float, default=ATON_CELL_COMOVING_KPC,
        help="ATON comoving cell size after converting 78.125 ckpc/h with h=0.678",
    )
    run.add_argument(
        "--cosmology-time-treatment", choices=("exact", "legacy"), default="exact",
        help="exact derives H0 in Myr^-1 from the code's unit constants",
    )

    assemble = subparsers.add_parser("assemble", help="Assemble completed text files into .npy")
    add_grid_selection(assemble)
    assemble.add_argument(
        "--allow-missing", action="store_true",
        help="assemble incomplete grids and record skipped LOS ids (disabled by default)",
    )

    return parser.parse_args()


def selected_snapshots(labels: list[str]) -> list[Snapshot]:
    wanted = set(labels)
    return [snapshot for snapshot in PAPER_SNAPSHOTS if snapshot.label in wanted]


def preflight_inputs(input_dir: Path, snapshots: list[Snapshot]) -> None:
    missing = []
    for snapshot in snapshots:
        for half in ("f1k", "l1k"):
            for field in FIELDS:
                path = input_dir / f"{half}z{snapshot.label}M3LOS_{field}.txt"
                if not path.is_file():
                    missing.append(path)
    if missing:
        lines = "\n".join(str(path) for path in missing[:20])
        suffix = "" if len(missing) <= 20 else f"\n... and {len(missing) - 20} more"
        raise SystemExit(f"missing {len(missing)} required input files:\n{lines}{suffix}")


def run_grid(args: argparse.Namespace, snapshots: list[Snapshot]) -> None:
    script = Path(__file__).with_name("run_durovcikova_grid.py")
    if args.cellsize_comoving_kpc <= 0:
        raise SystemExit("--cellsize-comoving-kpc must be positive")
    preflight_inputs(args.input_dir, snapshots)
    family = "noqso" if args.initial_only else "qso"

    for snapshot in snapshots:
        expansion = (1.0 + args.z_qso) / (1.0 + snapshot.redshift)
        comoving_length_scale = (
            1.0 if args.length_treatment == "fixed-comoving" else expansion ** 2
        )
        rescale_background = int(args.background_treatment == "equilibrium")
        apply_spectral_hydrogen_fraction = int(
            args.spectral_density_treatment == "hydrogen"
        )
        correct_spectral_units = int(args.spectral_unit_treatment == "consistent-ckpc")
        output_dir = args.output_root / family / f"xhi_z{snapshot.label}"
        print(
            f"snapshot z{snapshot.label}: z_snapshot={snapshot.redshift:.6f}, "
            f"paper_xHI={snapshot.paper_xhi:g}, density={expansion ** 3:.15g}, "
            f"comoving_length={comoving_length_scale:.15g}, "
            f"background={args.background_treatment}, "
            f"spectral_density={args.spectral_density_treatment}, "
            f"spectral_units={args.spectral_unit_treatment}, "
            f"helium_bug={args.legacy_alphahe2_bug}",
            flush=True,
        )
        # Complete one lifetime before starting the next. This lets users run
        # one snapshot per batch/screen while --workers parallelizes only the
        # sightlines belonging to the current lifetime.
        tqs = [args.tq_grid[0]] if args.initial_only else args.tq_grid
        for tq in tqs:
            if not args.initial_only:
                print(f"snapshot z{snapshot.label}: starting tQ={tq:g} Myr", flush=True)
            command = [
                sys.executable,
                str(script),
                "--rte", str(args.rte),
                "--input-dir", str(args.input_dir),
                "--input-prefix", f"f1kz{snapshot.label}M3LOS",
                "--second-input-prefix", f"l1kz{snapshot.label}M3LOS",
                "--output-dir", str(output_dir),
                # Runtime cosmology always stays at the observed quasar redshift.
                "--z-sim", f"{args.z_qso:.12g}",
                "--z-qso", f"{args.z_qso:.12g}",
                "--m1450", f"{args.m1450:.12g}",
                "--tq-grid", f"{tq:g}",
                "--n-cells", str(args.n_cells),
                "--n-los", str(args.n_los),
                "--split-los-at", str(args.split_los_at),
                "--legacy-density", str(args.legacy_density),
                "--density-multiplier", f"{expansion ** 3:.17g}",
                "--comoving-length-scale", f"{comoving_length_scale:.17g}",
                "--rescale-background", str(rescale_background),
                "--apply-spectral-hydrogen-fraction", str(apply_spectral_hydrogen_fraction),
                "--correct-spectral-units", str(correct_spectral_units),
                "--cellsize-override", f"{args.cellsize_comoving_kpc:.17g}",
                "--correct-cosmology-time-units",
                str(int(args.cosmology_time_treatment == "exact")),
                "--num-pixels-redward", str(args.num_pixels_redward),
                "--legacy-alphahe2-bug", str(args.legacy_alphahe2_bug),
                "--boxsize-override", f"{args.boxsize_override:.12g}",
                "--los", *args.los,
                "--workers", str(args.workers),
            ]
            if args.resume:
                command.append("--resume")
            if args.dry_run:
                command.append("--dry-run")
            if args.initial_only:
                command.append("--initial-only")

            run_child(command)


def assemble_grid(args: argparse.Namespace, snapshots: list[Snapshot]) -> None:
    script = Path(__file__).with_name("assemble_outputs_to_npy.py")
    family = "noqso" if args.initial_only else "qso"
    tqs = [0.0] if args.initial_only else args.tq_grid

    for snapshot in snapshots:
        input_dir = args.output_root / family / f"xhi_z{snapshot.label}"
        run_configuration_path = input_dir / ".run_configuration.json"
        if not run_configuration_path.is_file():
            raise SystemExit(
                f"missing verified run configuration: {run_configuration_path}; "
                "refusing to assemble an untracked production grid"
            )
        run_configuration = json.loads(run_configuration_path.read_text())
        required_configuration = {
            "schema_version",
            "rte_sha256",
            "input_signatures",
            "runtime_redshift",
            "quasar_redshift",
            "m1450",
            "n_cells",
            "initial_only",
            "input_prefix",
            "second_input_prefix",
            "num_pixels_redward",
            "comoving_length_scale",
            "density_multiplier",
            "rescale_background",
            "apply_spectral_hydrogen_fraction",
            "correct_spectral_units",
            "cellsize_override_comoving_kpc",
            "correct_cosmology_time_units",
            "legacy_alphahe2_bug",
            "boxsize_override",
        }
        missing_configuration = sorted(required_configuration - run_configuration.keys())
        if missing_configuration:
            raise SystemExit(
                f"incomplete run configuration {run_configuration_path}: missing "
                + ", ".join(missing_configuration)
            )
        expected_configuration = {
            "runtime_redshift": args.z_qso,
            "quasar_redshift": args.z_qso,
            "m1450": args.m1450,
            "n_cells": args.n_cells,
            "initial_only": args.initial_only,
            "input_prefix": f"f1kz{snapshot.label}M3LOS",
            "second_input_prefix": f"l1kz{snapshot.label}M3LOS",
        }
        mismatches = []
        for key, expected in expected_configuration.items():
            actual = run_configuration.get(key)
            if isinstance(expected, float):
                matches = bool(np.isclose(float(actual), expected, rtol=0.0, atol=1e-12))
            else:
                matches = actual == expected
            if not matches:
                mismatches.append(f"  {key}: run={actual!r}, assembly={expected!r}")
        if mismatches:
            raise SystemExit(
                f"assembly arguments do not match {run_configuration_path}:\n"
                + "\n".join(mismatches)
            )
        for tq in tqs:
            if args.initial_only:
                filename = f"noqso_z{snapshot.label}_zqso{args.z_qso:.2f}_los.npy"
                deleted_filename = (
                    f"deleted_indices_noqso_z{snapshot.label}_zqso{args.z_qso:.2f}.txt"
                )
            else:
                filename = (
                    f"qso_z{snapshot.label}_zqso{args.z_qso:.2f}_m{args.m1450:.2f}_"
                    f"tq{tq:.0e}_los.npy"
                )
                deleted_filename = (
                    f"deleted_indices_qso_z{snapshot.label}_zqso{args.z_qso:.2f}_"
                    f"m{args.m1450:.2f}_tq{tq:.0e}.txt"
                )
            output_path = args.output_root / "npy" / filename
            deleted_path = args.output_root / "npy" / deleted_filename
            command = [
                sys.executable,
                str(script),
                "--input-dir", str(input_dir),
                "--los", *args.los,
                "--output", str(output_path),
                "--z-qso", f"{args.z_qso:.12g}",
                "--m1450", f"{args.m1450:.12g}",
                "--tq", f"{tq:g}",
                "--n-cells", str(args.n_cells),
                "--deleted-indices-output",
                str(deleted_path),
            ]
            if args.allow_missing:
                command.append("--allow-missing")
            if args.initial_only:
                command.append("--initial-only")
            run_child(command)

            assembled = np.load(output_path, mmap_mode="r")
            distance_spacing = (
                float(assembled[0, 1, 0] - assembled[0, 0, 0])
                if assembled.shape[1] > 1 else None
            )
            hubble_h = 0.678
            omega_m = 0.308
            expansion_rate_factor = np.sqrt(
                omega_m * (1.0 + args.z_qso) ** 3 + (1.0 - omega_m)
            )
            hubble_z_km_s_Mpc = 100.0 * hubble_h * expansion_rate_factor
            correct_spectral_units = bool(run_configuration["correct_spectral_units"])
            solver_hubble_z_km_s_Mpc = (
                hubble_z_km_s_Mpc if correct_spectral_units
                else 100.0 / hubble_h * expansion_rate_factor
            )
            coordinate_velocity_spacing = (
                hubble_z_km_s_Mpc * distance_spacing
                / (1000.0 * (1.0 + args.z_qso))
                if distance_spacing is not None else None
            )
            solver_velocity_spacing = (
                solver_hubble_z_km_s_Mpc * distance_spacing
                / (1000.0 * (1.0 + args.z_qso))
                if distance_spacing is not None else None
            )
            lya_origin_index = int(run_configuration["num_pixels_redward"])
            if not 0 <= lya_origin_index < assembled.shape[1]:
                raise SystemExit(
                    f"invalid Ly-alpha origin index {lya_origin_index} for "
                    f"assembled shape {assembled.shape}"
                )

            metadata = {
                "schema_version": 2,
                "npy_file": output_path.name,
                "shape": list(assembled.shape),
                "dtype": str(assembled.dtype),
                "distance_spacing_comoving_kpc": distance_spacing,
                "hubble_parameter_km_s_Mpc": float(hubble_z_km_s_Mpc),
                "coordinate_hubble_velocity_spacing_km_s": (
                    None if coordinate_velocity_spacing is None
                    else float(coordinate_velocity_spacing)
                ),
                "lya_hubble_velocity_spacing_km_s": (
                    None if solver_velocity_spacing is None
                    else float(solver_velocity_spacing)
                ),
                "spectral_units_consistent_with_column_0": correct_spectral_units,
                "spectral_cosmology": {
                    "h": hubble_h,
                    "omega_m": omega_m,
                    "omega_lambda": 1.0 - omega_m,
                },
                "lya_origin_index": lya_origin_index,
                "lya_distance_comoving_kpc": (
                    "column_0 - column_0[lya_origin_index]"
                ),
                "lya_distance_proper_Mpc": (
                    "(column_0 - column_0[lya_origin_index]) / "
                    "(1000 * (1 + quasar_redshift))"
                ),
                "format_note": (
                    "Columns 1-5 are gas-cell fields at the stored column-0 coordinate; "
                    "column 6 is a spectral-bin field whose zero-Hubble-velocity origin "
                    "is lya_origin_index. This preserves the historical seven-column format."
                ),
                "columns": [
                    {
                        "index": 0,
                        "name": "stored_grid_coordinate",
                        "unit": "comoving kpc",
                    },
                    {"index": 1, "name": "xHI", "unit": "dimensionless"},
                    {"index": 2, "name": "temperature", "unit": "K"},
                    {"index": 3, "name": "xHeI", "unit": "dimensionless"},
                    {"index": 4, "name": "xHeII", "unit": "dimensionless"},
                    {"index": 5, "name": "xHeIII", "unit": "dimensionless"},
                    {"index": 6, "name": "lya_optical_depth", "unit": "dimensionless"},
                ],
                "snapshot_label": snapshot.label,
                "snapshot_redshift": snapshot.redshift,
                "paper_volume_averaged_xHI": snapshot.paper_xhi,
                "quasar_redshift": args.z_qso,
                "m1450": args.m1450,
                "quasar_lifetime_Myr": tq,
                "initial_only": args.initial_only,
                "allow_missing": args.allow_missing,
                "deleted_indices_file": deleted_path.name,
                "run_configuration": run_configuration,
            }
            metadata_path = output_path.with_suffix(".json")
            metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
            print(f"wrote {metadata_path}")


def main() -> None:
    args = parse_args()
    snapshots = selected_snapshots(args.snapshots)
    if args.action == "run":
        run_grid(args, snapshots)
    else:
        assemble_grid(args, snapshots)


if __name__ == "__main__":
    main()
