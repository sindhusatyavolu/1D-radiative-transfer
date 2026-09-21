#!/usr/bin/env python3
"""Run or assemble the eight-snapshot ATON grid from Durovcikova et al. III.1."""

from __future__ import annotations

import argparse
import signal
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


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
        raise subprocess.CalledProcessError(returncode, command)


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
    run.add_argument("--legacy-alphahe2-bug", type=int, choices=(0, 1), default=1)
    run.add_argument("--boxsize-override", type=float, default=102441.268)

    assemble = subparsers.add_parser("assemble", help="Assemble completed text files into .npy")
    add_grid_selection(assemble)

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
    preflight_inputs(args.input_dir, snapshots)
    family = "noqso" if args.initial_only else "qso"

    for snapshot in snapshots:
        expansion = (1.0 + args.z_qso) / (1.0 + snapshot.redshift)
        output_dir = args.output_root / family / f"xhi_z{snapshot.label}"
        print(
            f"snapshot z{snapshot.label}: z_snapshot={snapshot.redshift:.6f}, "
            f"paper_xHI={snapshot.paper_xhi:g}, density={expansion ** 3:.15g}, "
            f"comoving_length={expansion ** 2:.15g}",
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
                "--comoving-length-scale", f"{expansion ** 2:.17g}",
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
            command = [
                sys.executable,
                str(script),
                "--input-dir", str(input_dir),
                "--los", *args.los,
                "--output", str(args.output_root / "npy" / filename),
                "--z-qso", f"{args.z_qso:.12g}",
                "--m1450", f"{args.m1450:.12g}",
                "--tq", f"{tq:g}",
                "--n-cells", str(args.n_cells),
                "--allow-missing",
                "--deleted-indices-output",
                str(args.output_root / "npy" / deleted_filename),
            ]
            if args.initial_only:
                command.append("--initial-only")
            run_child(command)


def main() -> None:
    args = parse_args()
    snapshots = selected_snapshots(args.snapshots)
    if args.action == "run":
        run_grid(args, snapshots)
    else:
        assemble_grid(args, snapshots)


if __name__ == "__main__":
    main()
