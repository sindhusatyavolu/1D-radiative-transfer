#!/usr/bin/env python3
"""Regenerate direct RTE input grids from a dom_zloc binary skewer file.

This currently handles the fields that are direct reads or constant unit
conversions from the skewer file: nh, nh1, T, vpec, and rho.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


NCELLS = 2048
NLOS = 500
RECORD_BYTES = 114_760
HEADER_BYTES = 64
BLOCK_BYTES = NCELLS * 8
NBOX_MPC_H = 160.0
DR_MPC_H = 78.125e-3

FIELD_T = 1
FIELD_VPEC = 2
FIELD_DENSITY = 3
FIELD_NEUTRAL = 4
TIMEFAC = 3.15e13
KB = 1.38e-16
HELIUM_MASS_FRACTION = 0.24
XHE2_INITIAL = 1e-20
XHE3_INITIAL = 1e-20

# Inferred from get_los.py and the existing example files:
#   nh  = density_field * (frac_H / 1e6) * (1 + z_sim)^3
#   nh1 = neutral_field * (frac_H / 1e6) * (1 + z_sim)^3
#   rho = density_field * RHO_SCALE
HYDROGEN_MASS_FRACTION = 0.76
RHO_SCALE = 7.024355665137852e-9


def field_block(data: np.memmap, los: int, field_index: int) -> np.ndarray:
    start_byte = los * RECORD_BYTES + HEADER_BYTES + field_index * BLOCK_BYTES
    start = start_byte // 8
    return np.asarray(data[start : start + NCELLS])


def header_values(data: np.memmap, los: int) -> tuple[float, float]:
    start = los * RECORD_BYTES // 8
    zpos = float(data[start + 1])
    mhalo = float(data[start + 2])
    return zpos, mhalo


def halo_pixel_offset(zpos: float) -> int:
    return int(zpos * (NBOX_MPC_H - DR_MPC_H) / DR_MPC_H)


def shift_to_halo(raw: np.ndarray, offset: int) -> np.ndarray:
    if offset == 0:
        return raw.copy()
    return np.concatenate((raw[offset:], raw[:offset]))


def infer_shift_from_reference(raw_t: np.ndarray, reference_t: np.ndarray) -> int:
    matches = np.where(reference_t == raw_t[0])[0]
    for shift in matches:
        if np.array_equal(raw_t, np.concatenate((reference_t[shift:], reference_t[:shift]))):
            return int(shift)
    raise ValueError("could not infer circular shift from reference T column")


def get_offsets_and_masses(data: np.memmap) -> tuple[list[int], np.ndarray]:
    offsets = []
    masses = np.empty(NLOS)
    for los in range(NLOS):
        zpos, mhalo = header_values(data, los)
        offsets.append(halo_pixel_offset(zpos))
        masses[los] = mhalo
    return offsets, masses


def hydrogen_number_density_scale(z_sim: float) -> float:
    return HYDROGEN_MASS_FRACTION / 1e6 * (1.0 + z_sim) ** 3


def regenerate(data: np.memmap, offsets: list[int], z_sim: float) -> dict[str, np.ndarray]:
    grids = {
        "nh": np.empty((NCELLS, NLOS)),
        "nh1": np.empty((NCELLS, NLOS)),
        "T": np.empty((NCELLS, NLOS)),
        "vpec": np.empty((NCELLS, NLOS)),
        "rho": np.empty((NCELLS, NLOS)),
    }

    for los, offset in enumerate(offsets):
        density = shift_to_halo(field_block(data, los, FIELD_DENSITY), offset)
        neutral = shift_to_halo(field_block(data, los, FIELD_NEUTRAL), offset)

        grids["nh"][:, los] = density * hydrogen_number_density_scale(z_sim)
        grids["nh1"][:, los] = neutral * hydrogen_number_density_scale(z_sim)
        grids["T"][:, los] = shift_to_halo(field_block(data, los, FIELD_T), offset)
        grids["vpec"][:, los] = shift_to_halo(field_block(data, los, FIELD_VPEC), offset)
        grids["rho"][:, los] = density * RHO_SCALE

    return grids


def alpha_h2(T: np.ndarray) -> np.ndarray:
    lam = 2.0 * 157807.0 / T
    return TIMEFAC * 1.269e-13 * lam**1.503 / (1.0 + (lam / 0.522) ** 0.47) ** 1.923


def collion_h1(T: np.ndarray) -> np.ndarray:
    lam = 2.0 * 157807.0 / T
    return (
        TIMEFAC
        * 21.11
        * T**-1.5
        * np.exp(-lam / 2.0)
        * lam**-1.089
        / (1.0 + (lam / 0.354) ** 0.874) ** 1.101
    )


def collion_he1(T: np.ndarray) -> np.ndarray:
    lam = 2.0 * 285335.0 / T
    return (
        TIMEFAC
        * 32.38
        * T**-1.5
        * np.exp(-lam / 2.0)
        * lam**-1.146
        / (1.0 + (lam / 0.416) ** 0.987) ** 1.056
    )


def collion_he2(T: np.ndarray) -> np.ndarray:
    lam = 2.0 * 631515.0 / T
    return (
        TIMEFAC
        * 19.95
        * T**-1.5
        * np.exp(-lam / 2.0)
        * lam**-1.089
        / (1.0 + (lam / 0.553) ** 0.735) ** 1.275
    )


def alpha_he2(T: np.ndarray) -> np.ndarray:
    lam = 570670.0 / T
    return TIMEFAC * 3.0e-14 * lam**0.654


def alpha_he3(T: np.ndarray) -> np.ndarray:
    lam = 1263030.0 / T
    return TIMEFAC * 2.538e-13 * lam**1.503 / (1.0 + (lam / 0.522) ** 0.47) ** 1.923


def add_backgrounds(grids: dict[str, np.ndarray]) -> None:
    """Add background ionization and heating files recovered from backgroundgamma.py."""
    nh = grids["nh"]
    nh1 = grids["nh1"]
    T = grids["T"]

    nhe = nh * HELIUM_MASS_FRACTION / (4.0 * (1.0 - HELIUM_MASS_FRACTION))
    xhe2 = XHE2_INITIAL
    xhe3 = XHE3_INITIAL
    xhe1 = 1.0 - xhe2 - xhe3

    with np.errstate(divide="ignore", invalid="ignore"):
        xh1 = nh1 / nh
        xh2 = 1.0 - xh1
        ne = xh2 * nh + xhe2 * nhe + 2.0 * xhe3 * nhe
        nh2 = xh2 * nh

        grids["Gammabg"] = np.maximum(
            1e-40,
            (alpha_h2(T) * ne * xh2 / xh1 - ne * collion_h1(T)) / TIMEFAC,
        )
        grids["Gammabghe1"] = np.maximum(
            1e-40,
            (alpha_he2(T) * ne * xhe2 / xhe1 - ne * collion_he1(T)) / TIMEFAC,
        )
        grids["Gammabghe2"] = np.maximum(
            1e-40,
            (alpha_he3(T) * ne * xhe3 / xhe2 - ne * collion_he2(T)) / TIMEFAC,
        )

    T5 = T / 100000.0
    lamb = 2.0 * 157807.0 / T
    f1 = 1.778e-29 * TIMEFAC
    f2 = 7.5e-19 * TIMEFAC
    f3 = 1.42e-27 * TIMEFAC

    rec_h = f1 * T * lamb**1.965 / (1.0 + (lamb / 0.541) ** 0.502) ** 2.697 * ne * nh2
    coll_exc_h = f2 * (1.0 + T5**0.5) ** -1 * np.exp(-118355.0 / T) * ne * nh1
    free_free = 1.5 * f3 * T**0.5 * nh2 * ne
    grids["Heatbg"] = np.maximum(1e-40, rec_h + coll_exc_h + free_free)


def write_grids(grids: dict[str, np.ndarray], output_dir: Path, prefix: str) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    for name, grid in grids.items():
        np.savetxt(output_dir / f"{prefix}_{name}.txt", grid, fmt="%.18e")


def compare_grids(grids: dict[str, np.ndarray], reference_dir: Path, prefix: str) -> None:
    for name, grid in grids.items():
        reference = np.loadtxt(reference_dir / f"{prefix}_{name}.txt")
        finite = np.isfinite(grid) & np.isfinite(reference)
        nan_masks_match = np.array_equal(np.isnan(grid), np.isnan(reference))
        diff = np.abs(grid[finite] - reference[finite])
        rel = diff / (np.abs(reference[finite]) + 1e-300)
        allclose = nan_masks_match and np.allclose(grid[finite], reference[finite], rtol=0, atol=1e-14)
        print(
            f"{name:10s} finite={finite.sum()} nan_masks_match={nan_masks_match} "
            f"max_abs={diff.max():.6e} max_rel={rel.max():.6e} allclose={allclose}"
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("dom_file", type=Path)
    parser.add_argument("--reference-dir", type=Path, default=Path("example_input_files"))
    parser.add_argument("--output-dir", type=Path, default=Path("example_input_files/regenerated_direct"))
    parser.add_argument("--prefix", default="f1kz5.95M3LOS")
    parser.add_argument(
        "--z-sim",
        type=float,
        default=5.945518,
        help="snapshot redshift used for the physical hydrogen density conversion",
    )
    parser.add_argument("--compare", action="store_true")
    parser.add_argument("--write-halo-catalog", action="store_true")
    parser.add_argument(
        "--with-backgrounds",
        action="store_true",
        help="also generate Gammabg, Gammabghe1, Gammabghe2, and Heatbg",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    data = np.memmap(args.dom_file, dtype="<f8", mode="r")

    offsets, masses = get_offsets_and_masses(data)
    grids = regenerate(data, offsets, args.z_sim)
    if args.with_backgrounds:
        add_backgrounds(grids)
    write_grids(grids, args.output_dir, args.prefix)

    if args.write_halo_catalog:
        catalog = np.column_stack((np.arange(NLOS), offsets, masses))
        np.savetxt(
            args.output_dir / f"{args.prefix}_halo_catalog.txt",
            catalog,
            fmt=["%d", "%d", "%.18e"],
            header="LOS_index halo_pixel_offset Mhalo",
        )

    print(f"wrote regenerated direct inputs to {args.output_dir}")
    if args.compare:
        compare_grids(grids, args.reference_dir, args.prefix)


if __name__ == "__main__":
    main()
