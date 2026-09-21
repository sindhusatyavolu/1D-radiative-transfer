#!/usr/bin/env python3
"""Stack RTE text outputs into a qso_..._los.npy array."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def load_indices(path: Path) -> list[int]:
    return [int(line) for line in path.read_text().splitlines() if line.strip()]


def parse_los(values: list[str]) -> list[int]:
    result = []
    for value in values:
        if ":" in value:
            start, stop = map(int, value.split(":", 1))
            result.extend(range(start, stop + 1))
        else:
            result.append(int(value))
    return sorted(set(result))


def output_name(output_los: int, z_qso: float, tq: float, m1450: float, t_on: float, fduty: float) -> str:
    return f"z{z_qso:.2f}LOS{output_los}tq{tq:g}mag{m1450:.2f}_tep_{t_on:g}Myr{fduty:g}fd.txt"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--indices-file", type=Path)
    source.add_argument("--los", nargs="+", help="LOS ids or inclusive ranges such as 0:999")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--z-qso", type=float, required=True)
    parser.add_argument("--m1450", type=float, required=True)
    parser.add_argument("--tq", type=float, required=True)
    parser.add_argument("--t-on", type=float, default=None)
    parser.add_argument("--fduty", type=float, default=1.0)
    parser.add_argument("--n-cells", type=int, default=891)
    parser.add_argument("--output-los-mode", choices=("input", "row"), default="input")
    parser.add_argument("--initial-only", action="store_true",
                        help="Read noqso_z{z}LOS{id}.txt pre-quasar outputs")
    parser.add_argument(
        "--allow-missing",
        action="store_true",
        help="Skip missing/invalid sightlines and record their LOS ids instead of aborting",
    )
    parser.add_argument(
        "--deleted-indices-output",
        type=Path,
        help="Write skipped input LOS ids here; intended for legacy-style grids",
    )
    args = parser.parse_args()

    t_on = args.tq if args.t_on is None else args.t_on
    input_indices = load_indices(args.indices_file) if args.indices_file else parse_los(args.los)
    rows = []
    missing = []
    invalid = []
    deleted = []
    for row, los in enumerate(input_indices):
        output_los = los if args.output_los_mode == "input" else row
        name = (f"noqso_z{args.z_qso:.2f}LOS{output_los}.txt" if args.initial_only else
                output_name(output_los, args.z_qso, args.tq, args.m1450, t_on, args.fduty))
        path = args.input_dir / name
        if not path.exists() or path.stat().st_size == 0:
            missing.append(path)
            deleted.append(los)
            continue
        try:
            values = np.loadtxt(path)
        except (OSError, ValueError):
            invalid.append(path)
            deleted.append(los)
            continue
        if values.shape != (args.n_cells, 7) or not np.all(np.isfinite(values)):
            invalid.append(path)
            deleted.append(los)
            continue
        rows.append(values)

    if missing or invalid:
        print(f"missing={len(missing)} invalid={len(invalid)}")
        for path in (missing + invalid)[:20]:
            print(path)
        if not args.allow_missing:
            raise SystemExit(1)
    if not rows:
        raise SystemExit("no valid sightline outputs found")
    if args.deleted_indices_output is not None:
        args.deleted_indices_output.parent.mkdir(parents=True, exist_ok=True)
        if deleted:
            np.savetxt(args.deleted_indices_output, np.asarray(deleted, dtype=int), fmt="%d")
        else:
            args.deleted_indices_output.write_text("")
        print(f"wrote {args.deleted_indices_output} deleted={len(deleted)}")

    stack = np.asarray(rows, dtype=np.float64)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.save(args.output, stack)
    print(f"wrote {args.output} shape={stack.shape}")


if __name__ == "__main__":
    main()
