#!/usr/bin/env python3
"""Measure the Ly-alpha effect of changing only the archived RT redshift.

Read matched RT outputs and compare F=exp(-tau) on the shared z6.15 analysis
coordinate. Use the same coordinate and observed-frame 20-Angstrom boxcar
conventions as proximity_zone_p1d/p1d/{simulation,proximity}.py.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from astropy.cosmology import Planck18


def smooth_flux(wave, flux, width=20.0):
    """Average a centred observed-wavelength boxcar; exclude incomplete edges."""
    w = wave[::-1]
    f = flux[..., ::-1]
    lo = np.searchsorted(w, w - width / 2, side="left")
    hi = np.searchsorted(w, w + width / 2, side="right")
    cumulative = np.concatenate((np.zeros((*f.shape[:-1], 1)), np.cumsum(f, axis=-1)), axis=-1)
    result = (cumulative[..., hi] - cumulative[..., lo]) / (hi - lo)
    result[..., (w - width / 2 < w[0]) | (w + width / 2 > w[-1])] = np.nan
    return result[..., ::-1]


def describe(values):
    values = np.asarray(values)
    values = values[np.isfinite(values)]
    if not len(values):
        return {"count": 0}
    return {"count": len(values), "mean": float(np.mean(values)),
            "p16": float(np.percentile(values, 16)), "median": float(np.median(values)),
            "p84": float(np.percentile(values, 84)), "min": float(np.min(values)),
            "max": float(np.max(values))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--attachment", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    out = args.output_dir
    manifest = json.loads((out / "run_manifest.json").read_text())
    rows = np.asarray(manifest["rows"])
    ids = np.asarray(manifest["global_input_los"])
    reference = np.load(args.reference, mmap_mode="r")[rows]
    corrected = np.stack([
        np.loadtxt(out / "correct_z6p15" / f"z6.15LOS{los}tq1mag-27.09_tep_1Myr1fd.txt")
        for los in ids
    ])
    if corrected.shape != reference.shape or not np.all(np.isfinite(corrected)):
        raise ValueError("corrected runs must be finite and have the reference shape")
    if not np.array_equal(corrected[:, :, 0], reference[:, :, 0]):
        raise ValueError("reference and corrected coordinates differ")
    if np.any(corrected[:, :, 6] < 0):
        raise ValueError("negative corrected optical depths")
    z = manifest["rt_redshift"]
    origin = manifest["num_pixels_redward"]
    r = (reference[0, :, 0] - reference[0, origin, 0]) / (1000 * (1 + z))
    velocity = Planck18.H(z).value * r
    wave = 1215.67 * (1 + z) * np.exp(-velocity / 299792.458)
    archive_flux = np.exp(-reference[:, :, 6])
    correct_flux = np.exp(-corrected[:, :, 6])
    delta = archive_flux - correct_flux
    sm_old = smooth_flux(wave, archive_flux)
    sm_new = smooth_flux(wave, correct_flux)
    blue = (r >= 0) & (r < 10)
    per_los_mae = np.mean(np.abs(delta[:, blue]), axis=1)
    report = {
        "manifest": manifest,
        "definition": "delta_F = archive F - corrected-redshift F; F=exp(-tau)",
        "coordinate": "Both arrays use the same nominal z=6.148904 distance, origin at pixel 256",
        "scope": f"Redshift-only correction, retaining all other legacy settings; {len(rows)} distinct input profiles",
        "blueward_0_10_pMpc": {
            "mean_absolute_flux_difference": float(np.mean(np.abs(delta[:, blue]))),
            "rms_flux_difference": float(np.sqrt(np.mean(delta[:, blue]**2))),
            "signed_mean_flux_difference": float(np.mean(delta[:, blue])),
            "absolute_pixel_difference_p95": float(np.percentile(np.abs(delta[:, blue]), 95)),
            "max_absolute_pixel_difference": float(np.max(np.abs(delta[:, blue]))),
            "pixel_fraction_abs_difference_gt_0p05": float(np.mean(np.abs(delta[:, blue]) > 0.05)),
            "pixel_fraction_abs_difference_gt_0p1": float(np.mean(np.abs(delta[:, blue]) > 0.1)),
            "per_sightline_mean_absolute_difference": describe(per_los_mae),
            "per_sightline_max_absolute_difference": describe(np.max(np.abs(delta[:, blue]), axis=1)),
            "smoothed_20A_mean_absolute_difference": float(np.nanmean(np.abs((sm_old-sm_new)[:, blue]))),
        },
        "radial_bins": [],
    }
    for lo, hi in [(0, 2), (2, 5), (5, 10), (0, 10)]:
        mask = (r >= lo) & (r < hi)
        old_mean = float(np.mean(archive_flux[:, mask]))
        new_mean = float(np.mean(correct_flux[:, mask]))
        report["radial_bins"].append({
            "range_pMpc": [lo, hi], "archive_mean_flux": old_mean,
            "corrected_mean_flux": new_mean, "archive_minus_corrected": old_mean-new_mean,
            "relative_bias_percent": 100 * (old_mean-new_mean) / new_mean,
            "mean_absolute_pixel_difference": float(np.mean(np.abs(delta[:, mask]))),
        })

    def rp(smoothed):
        indices = np.flatnonzero((r >= 0) & np.isfinite(smoothed) & (smoothed <= 0.1))
        return float(r[indices[0]]) if len(indices) else np.nan

    rp_old = np.array([rp(s) for s in sm_old])
    rp_new = np.array([rp(s) for s in sm_new])
    report["Rp_20A_threshold_0p1"] = {
        "archive_pMpc": describe(rp_old), "corrected_pMpc": describe(rp_new),
        "archive_minus_corrected_pMpc": describe(rp_old-rp_new),
        "absolute_difference_pMpc": describe(abs(rp_old-rp_new)),
        "unmeasurable_archive": int(np.sum(~np.isfinite(rp_old))),
        "unmeasurable_corrected": int(np.sum(~np.isfinite(rp_new))),
    }
    if args.attachment:
        with args.reference.open("rb") as file:
            archive_hash = hashlib.file_digest(file, "sha256").hexdigest()
        with args.attachment.open("rb") as file:
            attachment_hash = hashlib.file_digest(file, "sha256").hexdigest()
        report["attachment"] = {"path": str(args.attachment), "sha256": attachment_hash,
                                "archive_sha256": archive_hash, "identical_to_archive": archive_hash == attachment_hash}
    (out / "profile_impact.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    np.save(out / "corrected_redshift_subset.npy", corrected)
    np.savez(out / "matched_flux_profiles.npz", target_rows=rows, global_input_los=ids,
             distance_pMpc=r, archive_flux=archive_flux, corrected_flux=correct_flux,
             archive_smoothed_20A=sm_old, corrected_smoothed_20A=sm_new)
    with (out / "per_sightline.csv").open("w") as file:
        writer = csv.writer(file)
        writer.writerow(["target_row", "input_los", "mean_abs_delta_F_0_10pMpc", "max_abs_delta_F_0_10pMpc", "Rp_archive_pMpc", "Rp_corrected_pMpc"])
        writer.writerows(zip(rows, ids, per_los_mae, np.max(abs(delta[:, blue]), axis=1), rp_old, rp_new))
    plot(out, r, archive_flux, correct_flux, sm_old, sm_new, rows, ids, per_los_mae)
    print(json.dumps(report, indent=2))


def plot(out, r, old, new, sm_old, sm_new, rows, ids, mae):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    old_color, new_color = "#b54837", "#2167a5"
    fig, axes = plt.subplots(3, 1, figsize=(10, 9), constrained_layout=True)
    axes[0].plot(r, np.mean(old, axis=0), color=old_color, alpha=0.45, lw=0.7)
    axes[0].plot(r, np.mean(new, axis=0), color=new_color, alpha=0.45, lw=0.7)
    axes[0].plot(r, np.mean(sm_old, axis=0), color=old_color, label="Archive: RT z=6.604205")
    axes[0].plot(r, np.mean(sm_new, axis=0), color=new_color, label="Redshift corrected: RT z=6.148904")
    axes[0].set_ylabel("Mean transmitted flux")
    axes[0].set_title(f"{len(rows)} distinct input skewers · M1450 = −27.09 · tQ = 1 Myr")
    axes[0].legend(fontsize=9)
    axes[0].text(0.98, 0.95, "Solid: 20 Å observed-frame smoothing\nFaint: native pixels", transform=axes[0].transAxes, ha="right", va="top", fontsize=8)
    difference = old-new
    axes[1].plot(r, 100*np.mean(difference, axis=0), color="#4a376c", lw=0.7, alpha=0.5)
    axes[1].plot(r, 100*np.mean(sm_old-sm_new, axis=0), color="#4a376c", label="Smoothed difference of means")
    axes[1].axhline(0, color="gray", lw=0.8)
    axes[1].set_ylabel("Archive − corrected\n(flux percentage points)")
    axes[1].legend(fontsize=9)
    selected = np.argsort(mae)[len(mae)//2]
    axes[2].plot(r, old[selected], color=old_color, lw=1, label="Archive")
    axes[2].plot(r, new[selected], color=new_color, lw=1, label="Redshift corrected")
    axes[2].set_title(f"Example nearest the median mean absolute difference: input LOS {ids[selected]} / archive row {rows[selected]}", fontsize=10)
    axes[2].set_ylabel("Transmitted flux")
    axes[2].legend(fontsize=9)
    for ax in axes:
        ax.set_xlim(0, 10)
        ax.set_xlabel("Distance from quasar (proper Mpc, shared z=6.148904 analysis coordinate)")
        ax.grid(alpha=0.2)
    fig.savefig(out / "lya_profile_impact.png", dpi=180)
    fig.savefig(out / "lya_profile_impact.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
