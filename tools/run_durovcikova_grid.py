#!/usr/bin/env python3
"""Run a Durovcikova-style post-processing grid with the 1D RTE executable.

The C executable remains a single-sightline runner. This script handles the
outer loops over snapshots, neutral-fraction labels, quasar properties, tq, and
LOS index.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import math
import struct
import subprocess
from pathlib import Path


def parse_los(values: list[str]) -> list[int]:
    los: list[int] = []
    for value in values:
        if ":" in value:
            start_s, stop_s = value.split(":", 1)
            start = int(start_s)
            stop = int(stop_s)
            los.extend(range(start, stop + 1))
        else:
            los.append(int(value))
    return sorted(set(los))


def load_los_file(path: Path) -> list[int]:
    los = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            los.append(int(line))
    return los


def run_command(command: list[str], dry_run: bool, log: Path | None = None) -> None:
    print(" ".join(command))
    if not dry_run:
        if log is None:
            subprocess.run(command, check=True)
        else:
            log.parent.mkdir(parents=True, exist_ok=True)
            with log.open("w") as stream:
                subprocess.run(command, check=True, stdout=stream, stderr=subprocess.STDOUT)


def output_name(output_los: int, z_qso: float, tq: float, magnitude: float,
                t_on: float, fduty: float, initial_only: bool) -> str:
    if initial_only:
        return f"noqso_z{z_qso:.2f}LOS{output_los}.txt"
    return (f"z{z_qso:.2f}LOS{output_los}tq{tq:g}mag{magnitude:.2f}"
            f"_tep_{t_on:g}Myr{fduty:g}fd.txt")


def complete_output(path: Path, n_cells: int) -> bool:
    if not path.exists() or path.stat().st_size == 0:
        return False
    rows = 0
    try:
        with path.open() as stream:
            for line in stream:
                values = line.split()
                if len(values) != 7 or not all(math.isfinite(float(value)) for value in values):
                    return False
                rows += 1
    except (OSError, ValueError):
        return False
    return rows == n_cells


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rte", type=Path, default=Path("./rte1d"))
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--input-prefix", required=True)
    parser.add_argument("--second-input-prefix", help="Prefix for global LOS >= --split-los-at")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--z-sim", type=float, required=True,
                        help="Runtime cosmology/spectral redshift (historical option name)")
    parser.add_argument("--z-qso", type=float, required=True,
                        help="Quasar redshift and output filename label")
    parser.add_argument("--m1450", type=float, required=True)
    parser.add_argument(
        "--legacy-float-magnitude",
        action="store_true",
        help="Round M1450 to float32 before computing luminosity, as in the archived C runner",
    )
    parser.add_argument("--tq-grid", type=float, nargs="+", required=True)
    parser.add_argument("--los", nargs="+", help="LOS ids or inclusive ranges like 0:499")
    parser.add_argument("--los-file", type=Path, help="File containing one input LOS id per line")
    parser.add_argument("--t-on", type=float, default=None)
    parser.add_argument("--fduty", type=float, default=1.0)
    parser.add_argument("--n-cells", type=int, default=891)
    parser.add_argument("--n-los", type=int, default=500)
    parser.add_argument("--split-los-at", type=int, default=500)
    parser.add_argument(
        "--legacy-density", type=int, choices=(0, 1), default=1,
        help="1 keeps supplied nH/nHI; 0 removes the legacy X=0.76 factor and applies z-qso/z-sim cubed",
    )
    parser.add_argument("--num-pixels-redward", type=int, default=256)
    parser.add_argument("--legacy-alphahe2-bug", type=int, choices=(0, 1), default=1)
    parser.add_argument("--boxsize-override", type=float, default=102441.268)
    parser.add_argument(
        "--comoving-length-scale", type=float, default=1.0,
        help="Multiply RT cell size, inner radius, and spectral box length; density is separate",
    )
    parser.add_argument(
        "--density-multiplier", type=float, default=1.0,
        help="Additional multiplier applied to physical nH and nHI after legacy density handling",
    )
    parser.add_argument(
        "--output-los-mode",
        choices=("input", "row"),
        default="input",
        help="Use input LOS id or row number in output filenames",
    )
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--workers", type=int, default=1,
                        help="Number of independent sightlines to run concurrently")
    parser.add_argument("--log-dir", type=Path,
                        help="Write one log per run; defaults to OUTPUT_DIR/logs when workers > 1")
    parser.add_argument("--resume", action="store_true",
                        help="Skip existing finite seven-column outputs with exactly --n-cells rows")
    parser.add_argument(
        "--initial-only", action="store_true",
        help="Synthesize the pre-quasar snapshot without evolving it; ignores the supplied tq grid",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.workers < 1:
        raise SystemExit("--workers must be positive")
    if args.density_multiplier <= 0:
        raise SystemExit("--density-multiplier must be positive")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if args.los_file is not None:
        los_values = load_los_file(args.los_file)
    elif args.los is not None:
        los_values = parse_los(args.los)
    else:
        raise SystemExit("provide --los or --los-file")
    rte = args.rte.resolve() if not args.rte.is_absolute() else args.rte
    magnitude = args.m1450
    if args.legacy_float_magnitude:
        magnitude = struct.unpack("f", struct.pack("f", magnitude))[0]

    jobs: list[tuple[list[str], Path | None]] = []
    for tq in ([0.0] if args.initial_only else args.tq_grid):
        t_on = args.t_on if args.t_on is not None else (1.0 if args.initial_only else tq)
        for row, los in enumerate(los_values):
            input_los = los
            input_prefix = args.input_prefix
            if args.second_input_prefix is not None and los >= args.split_los_at:
                input_los = los - args.split_los_at
                input_prefix = args.second_input_prefix
            output_los = los if args.output_los_mode == "input" else row
            command = [
                str(rte),
                f"{tq:g}",
                f"{t_on:g}",
                f"{args.fduty:g}",
                str(input_los),
                f"{magnitude:.17g}",
                f"{args.z_sim:.12g}",
                f"{args.z_qso:.12g}",
                str(args.input_dir),
                input_prefix,
                str(args.output_dir),
                str(args.n_cells),
                str(args.n_los),
                str(output_los),
                str(args.legacy_density),
                str(args.num_pixels_redward),
                str(args.legacy_alphahe2_bug),
                f"{args.boxsize_override:.12g}",
                f"{args.comoving_length_scale:.17g}",
                str(int(args.initial_only)),
                f"{args.density_multiplier:.17g}",
            ]
            log_dir = args.log_dir
            if log_dir is None and args.workers > 1:
                log_dir = args.output_dir / "logs"
            log = None if log_dir is None else log_dir / f"LOS{output_los}_tq{tq:g}.log"
            output = args.output_dir / output_name(
                output_los, args.z_qso, tq, magnitude, t_on, args.fduty, args.initial_only)
            if args.resume and complete_output(output, args.n_cells):
                print(f"skip complete {output}")
                continue
            jobs.append((command, log))

    if args.dry_run or args.workers == 1:
        for command, log in jobs:
            run_command(command, args.dry_run, log)
        return
    if not jobs:
        print("all requested outputs are already complete")
        return
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(run_command, command, False, log): (command, log)
                   for command, log in jobs}
        completed = 0
        for future in as_completed(futures):
            command, log = futures[future]
            try:
                future.result()
            except Exception as error:
                raise RuntimeError(f"Run failed; see {log}: {' '.join(command)}") from error
            completed += 1
            print(f"completed {completed}/{len(jobs)}", flush=True)


if __name__ == "__main__":
    main()
