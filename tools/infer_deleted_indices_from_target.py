#!/usr/bin/env python3
"""Infer retained LOS indices by aligning a target npy file to input grids."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def load_full_grid(input_dir: Path, suffix: str) -> np.ndarray:
    first = np.loadtxt(input_dir / f"f1kz6.15M3LOS_{suffix}.txt")
    second = np.loadtxt(input_dir / f"l1kz6.15M3LOS_{suffix}.txt")
    return np.hstack([first, second])


def alignment_cost(
    target: np.ndarray,
    input_t: np.ndarray,
    input_xhi: np.ndarray,
    first_cell: int,
    last_cell: int,
) -> np.ndarray:
    cells = np.arange(first_cell, last_cell)
    input_log_t = np.log10(np.clip(input_t[:last_cell, :].T[:, cells], 1.0, 1e9))
    input_log_xhi = np.log10(np.clip(input_xhi[:last_cell, :].T[:, cells], 1e-12, 1.0))
    target_log_t = np.log10(np.clip(target[:, cells, 2], 1.0, 1e9))
    target_log_xhi = np.log10(np.clip(target[:, cells, 1], 1e-12, 1.0))

    cost = np.empty((target.shape[0], input_t.shape[1]), dtype=np.float64)
    for row in range(target.shape[0]):
        d_t = np.mean(np.abs(input_log_t - target_log_t[row]), axis=1)
        d_xhi = np.mean(np.abs(input_log_xhi - target_log_xhi[row]), axis=1)
        cost[row] = d_t + 0.25 * d_xhi
    return cost


def monotonic_alignment(cost: np.ndarray) -> list[int]:
    n_rows, n_los = cost.shape
    inf = 1e100
    dp = np.full((n_rows + 1, n_los + 1), inf)
    prev = np.zeros((n_rows + 1, n_los + 1), dtype=np.uint8)
    dp[0, 0] = 0.0

    for los in range(n_los):
        skip = dp[:, los] < dp[:, los + 1]
        dp[skip, los + 1] = dp[skip, los]
        prev[skip, los + 1] = 2

        values = dp[:n_rows, los] + cost[:, los]
        keep = values < dp[1:, los + 1]
        dp[1:, los + 1][keep] = values[keep]
        prev[1:, los + 1][keep] = 1

    row = n_rows
    los = n_los
    kept = []
    while los > 0:
        if prev[row, los] == 1:
            kept.append(los - 1)
            row -= 1
            los -= 1
        else:
            los -= 1

    kept.reverse()
    if len(kept) != n_rows:
        raise RuntimeError(f"expected {n_rows} matched rows, got {len(kept)}")
    return kept


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", type=Path, help="Target qso_..._los.npy file")
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("/Users/ssatyavolu/projects/fireqso_dominika"),
        help="Directory containing f1/l1 z6.15 input text files",
    )
    parser.add_argument("--output-dir", type=Path, default=Path("outputs/recover_z6p15"))
    parser.add_argument("--first-cell", type=int, default=591)
    parser.add_argument("--last-cell", type=int, default=891)
    args = parser.parse_args()

    target = np.load(args.target, mmap_mode="r")
    input_t = load_full_grid(args.input_dir, "T")
    nh = load_full_grid(args.input_dir, "nh")
    nh1 = load_full_grid(args.input_dir, "nh1")
    input_xhi = nh1 / nh

    cost = alignment_cost(target, input_t, input_xhi, args.first_cell, args.last_cell)
    kept = monotonic_alignment(cost)
    deleted = [los for los in range(input_t.shape[1]) if los not in set(kept)]

    args.output_dir.mkdir(parents=True, exist_ok=True)
    kept_path = args.output_dir / "inferred_kept_indices_z6.15_m27.09_from_target.txt"
    deleted_path = args.output_dir / "inferred_deleted_indices_z6.15_m27.09_from_target.txt"
    np.savetxt(kept_path, np.asarray(kept, dtype=int), fmt="%d")
    np.savetxt(deleted_path, np.asarray(deleted, dtype=int), fmt="%d")

    print(f"target_shape={target.shape}")
    print(f"kept={len(kept)} deleted={len(deleted)}")
    print(f"first_deleted={deleted[:80]}")
    print(f"wrote {kept_path}")
    print(f"wrote {deleted_path}")


if __name__ == "__main__":
    main()
