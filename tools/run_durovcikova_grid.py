#!/usr/bin/env python3
"""Run a Durovcikova-style post-processing grid with the 1D RTE executable.

The C executable remains a single-sightline runner. This script handles the
outer loops over snapshots, neutral-fraction labels, quasar properties, tq, and
LOS index.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
import math
import struct
import subprocess
import threading
import time
from pathlib import Path


_active_processes: set[subprocess.Popen] = set()
_active_processes_lock = threading.Lock()

INPUT_FIELDS = (
    "nh", "nh1", "T", "vpec", "rho", "Gammabg", "Gammabghe1",
    "Gammabghe2", "Heatbg",
)
CONFIGURATION_FILENAME = ".run_configuration.json"


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


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def input_signatures(input_dir: Path, prefixes: list[str]) -> list[dict[str, object]]:
    signatures = []
    for prefix in prefixes:
        for field in INPUT_FIELDS:
            path = (input_dir / f"{prefix}_{field}.txt").resolve()
            stat = path.stat()
            signatures.append({
                "path": str(path),
                "size": stat.st_size,
                "mtime_ns": stat.st_mtime_ns,
            })
    return signatures


def output_configuration(args: argparse.Namespace, rte: Path) -> dict[str, object]:
    prefixes = [args.input_prefix]
    if args.second_input_prefix is not None:
        prefixes.append(args.second_input_prefix)
    return {
        "schema_version": 2,
        "rte_path": str(rte),
        "rte_sha256": sha256_file(rte),
        "input_signatures": input_signatures(args.input_dir, prefixes),
        "input_prefix": args.input_prefix,
        "second_input_prefix": args.second_input_prefix,
        "runtime_redshift": args.z_sim,
        "quasar_redshift": args.z_qso,
        "m1450": args.m1450,
        "legacy_float_magnitude": args.legacy_float_magnitude,
        "fduty": args.fduty,
        "t_on_mode": "match_lifetime" if args.t_on is None else "fixed",
        "fixed_t_on_Myr": args.t_on,
        "n_cells": args.n_cells,
        "n_los_per_input_file": args.n_los,
        "split_los_at": args.split_los_at,
        "legacy_density": args.legacy_density,
        "density_multiplier": args.density_multiplier,
        "rescale_background": args.rescale_background,
        "apply_spectral_hydrogen_fraction": args.apply_spectral_hydrogen_fraction,
        "correct_spectral_units": args.correct_spectral_units,
        "cellsize_override_comoving_kpc": args.cellsize_override,
        "correct_cosmology_time_units": args.correct_cosmology_time_units,
        "num_pixels_redward": args.num_pixels_redward,
        "legacy_alphahe2_bug": args.legacy_alphahe2_bug,
        "boxsize_override": args.boxsize_override,
        "comoving_length_scale": args.comoving_length_scale,
        "initial_only": args.initial_only,
        "output_los_mode": args.output_los_mode,
    }


def ensure_output_configuration(output_dir: Path, configuration: dict[str, object]) -> None:
    path = output_dir / CONFIGURATION_FILENAME
    if path.exists():
        existing = json.loads(path.read_text())
        if existing != configuration:
            keys = sorted(set(existing) | set(configuration))
            differences = [
                f"  {key}: existing={existing.get(key)!r} requested={configuration.get(key)!r}"
                for key in keys if existing.get(key) != configuration.get(key)
            ]
            raise SystemExit(
                f"refusing to mix incompatible runs in {output_dir}:\n" + "\n".join(differences)
            )
        return

    existing_outputs = list(output_dir.glob("*.txt"))
    if existing_outputs:
        raise SystemExit(
            f"refusing to resume untracked outputs in {output_dir}; choose a new output directory "
            "so the scientific configuration can be verified"
        )
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(configuration, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def run_command(
    command: list[str],
    dry_run: bool,
    log: Path | None = None,
    expected_output: Path | None = None,
    n_cells: int | None = None,
) -> None:
    print(" ".join(command))
    if dry_run:
        return

    stream = None
    process = None
    try:
        if log is not None:
            log.parent.mkdir(parents=True, exist_ok=True)
            stream = log.open("w")
        process = subprocess.Popen(
            command,
            stdout=stream,
            stderr=subprocess.STDOUT if stream is not None else None,
        )
        with _active_processes_lock:
            _active_processes.add(process)
        returncode = process.wait()
        if returncode != 0:
            raise subprocess.CalledProcessError(returncode, command)
    finally:
        if process is not None:
            with _active_processes_lock:
                _active_processes.discard(process)
        if stream is not None:
            stream.close()

    if expected_output is not None:
        if n_cells is None:
            raise ValueError("n_cells is required when validating solver output")
        if not complete_output(expected_output, n_cells):
            raise RuntimeError(
                f"solver returned success but did not create a complete output: {expected_output}"
            )


def terminate_active_processes() -> None:
    """Stop solver children when the driver is interrupted or one job fails."""
    with _active_processes_lock:
        processes = list(_active_processes)
    for process in processes:
        if process.poll() is None:
            process.terminate()
    deadline = time.monotonic() + 5
    for process in processes:
        if process.poll() is not None:
            continue
        try:
            process.wait(timeout=max(0, deadline - time.monotonic()))
        except subprocess.TimeoutExpired:
            process.kill()
            process.wait()


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
        "--rescale-background", type=int, choices=(0, 1), default=0,
        help=("1 scales photoionization rates with the total applied density factor and "
              "background heating with its square, preserving the input equilibrium"),
    )
    parser.add_argument(
        "--apply-spectral-hydrogen-fraction", type=int, choices=(0, 1), default=0,
        help=("1 multiplies total-gas rho by X=0.76 when synthesizing HI Ly-alpha; "
              "0 retains the legacy normalization"),
    )
    parser.add_argument(
        "--correct-spectral-units", type=int, choices=(0, 1), default=0,
        help=("1 treats cellsize as ckpc and rho as Gadget 1e10 Msun h^2/kpc^3 "
              "when computing Ly-alpha; 0 retains the legacy h-factor mismatch"),
    )
    parser.add_argument(
        "--cellsize-override", type=float, default=0.0,
        help="Comoving cell size in kpc; 0 retains the compiled legacy value",
    )
    parser.add_argument(
        "--correct-cosmology-time-units", type=int, choices=(0, 1), default=0,
        help="1 derives H0 in Myr^-1 from the code's s-to-Myr conversion",
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
    if args.cellsize_override < 0:
        raise SystemExit("--cellsize-override must be non-negative")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if args.los_file is not None:
        los_values = load_los_file(args.los_file)
    elif args.los is not None:
        los_values = parse_los(args.los)
    else:
        raise SystemExit("provide --los or --los-file")
    rte = args.rte.resolve() if not args.rte.is_absolute() else args.rte
    if not rte.is_file():
        raise SystemExit(f"RTE executable does not exist: {rte}")
    if not args.dry_run:
        ensure_output_configuration(args.output_dir, output_configuration(args, rte))
    magnitude = args.m1450
    if args.legacy_float_magnitude:
        magnitude = struct.unpack("f", struct.pack("f", magnitude))[0]

    jobs: list[tuple[list[str], Path | None, Path]] = []
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
                str(args.rescale_background),
                str(args.apply_spectral_hydrogen_fraction),
                str(args.correct_spectral_units),
                f"{args.cellsize_override:.17g}",
                str(args.correct_cosmology_time_units),
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
            jobs.append((command, log, output))

    if args.dry_run or args.workers == 1:
        for command, log, output in jobs:
            run_command(command, args.dry_run, log, output, args.n_cells)
        return
    if not jobs:
        print("all requested outputs are already complete")
        return
    pool = ThreadPoolExecutor(max_workers=args.workers)
    futures = {
        pool.submit(run_command, command, False, log, output, args.n_cells):
        (command, log, output)
        for command, log, output in jobs
    }
    try:
        completed = 0
        for future in as_completed(futures):
            command, log, _ = futures[future]
            try:
                future.result()
            except Exception as error:
                raise SystemExit(
                    f"run failed; see {log}: {' '.join(command)}\n{error}"
                ) from None
            completed += 1
            print(f"completed {completed}/{len(jobs)}", flush=True)
    except BaseException:
        for future in futures:
            future.cancel()
        terminate_active_processes()
        pool.shutdown(wait=True, cancel_futures=True)
        raise
    else:
        pool.shutdown(wait=True)


if __name__ == "__main__":
    main()
