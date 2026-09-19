#!/usr/bin/env python3
"""Infer the monotonically retained input LOS indices of an archived array."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from infer_deleted_indices_from_target import alignment_cost, monotonic_alignment


def load_grid(input_dir: Path, label: str, field: str, n_cells: int) -> np.ndarray:
    halves = [
        np.loadtxt(
            input_dir / f"{prefix}z{label}M3LOS_{field}.txt",
            max_rows=n_cells,
        )
        for prefix in ("f1k", "l1k")
    ]
    return np.hstack(halves)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", type=Path)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--snapshot-label", required=True, help="Input label such as 5.11")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--first-cell", type=int, default=591)
    parser.add_argument("--last-cell", type=int, default=891)
    args = parser.parse_args()

    target = np.load(args.target, mmap_mode="r")
    if target.ndim != 3 or target.shape[2] != 7:
        raise SystemExit(f"expected target shape (rows, cells, 7), got {target.shape}")
    n_cells = target.shape[1]
    if not 0 <= args.first_cell < args.last_cell <= n_cells:
        raise SystemExit("require 0 <= --first-cell < --last-cell <= target cells")

    temperature = load_grid(args.input_dir, args.snapshot_label, "T", n_cells)
    nh = load_grid(args.input_dir, args.snapshot_label, "nh", n_cells)
    nhi = load_grid(args.input_dir, args.snapshot_label, "nh1", n_cells)
    xhi = np.divide(nhi, nh, out=np.zeros_like(nhi), where=nh != 0)

    cost = alignment_cost(
        target,
        temperature,
        xhi,
        args.first_cell,
        args.last_cell,
    )
    kept = np.asarray(monotonic_alignment(cost), dtype=int)
    chosen = cost[np.arange(len(kept)), kept]
    kept_set = set(kept.tolist())
    deleted = [los for los in range(temperature.shape[1]) if los not in kept_set]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(args.output, kept, fmt="%d")
    report = {
        "target": str(args.target.resolve()),
        "target_shape": list(target.shape),
        "snapshot_label": args.snapshot_label,
        "matching_cells": [args.first_cell, args.last_cell],
        "retained": len(kept),
        "deleted": len(deleted),
        "deleted_indices": deleted,
        "median_chosen_cost": float(np.median(chosen)),
        "maximum_chosen_cost": float(np.max(chosen)),
    }
    report_path = args.report or args.output.with_suffix(".json")
    report_path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps(report, indent=2, allow_nan=False))
    print(f"wrote {args.output}")
    print(f"wrote {report_path}")


if __name__ == "__main__":
    main()
