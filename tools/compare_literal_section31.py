#!/usr/bin/env python3
"""Run the printed Section 3.1 redshift transformations on four matched skewers.

Requires the previously validated downloaded_m26p99_z5p76 benchmark inputs.
The two background variants expose a choice not specified in Section 3.1.
Other physics and unit conventions remain those of the validated legacy run.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

import numpy as np

from diagnose_aton_reproduction import comparison_metrics, optical_depth
from generate_inputs_from_dom import add_backgrounds

ROOT = Path(__file__).resolve().parents[1]
BENCH = ROOT / "outputs/downloaded_m26p99_z5p76"
OUT = ROOT / "outputs/section31_zsim5p76_zqso6p14_m26p99"
TARGET = Path("/Users/ssatyavolu/Downloads/qso_z5.76_m-26.99_tq1e+00_los.npy")
ZS, ZQ, ZINPUT, ZOLD = 5.76, 6.14, 5.755954, 6.604205
A = (1 + ZQ) / (1 + ZS)
D, C = A**3, A**2
BOX, NCELLS, EXTRA = 102441.268, 891, 256
SAMPLES = [(0, 0, "original_density_runtime6604"), (314, 335, "sample_row314"),
           (629, 672, "sample_row629"), (944, 999, "sample_row944")]
FIELDS = ("nh", "nh1", "T", "vpec", "rho", "Gammabg", "Gammabghe1", "Gammabghe2", "Heatbg")


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def prepare():
    OUT.mkdir(parents=True, exist_ok=True)
    jobs = []
    normalization = ((1 + ZS) / (1 + ZINPUT))**3
    for variant in ("fixed_background", "recomputed_background", "baseline_regression"):
        for row, los, folder in SAMPLES:
            if variant == "baseline_regression" and row != 0:
                continue
            source = BENCH / folder / "inputs"
            original = {key: np.loadtxt(source / f"los{los}_{key}.txt") for key in FIELDS}
            grids = {key: value.copy() for key, value in original.items()}
            scale, zrun = (1., ZOLD) if variant == "baseline_regression" else (C, ZQ)
            if variant != "baseline_regression":
                # The label is 5.76 but existing physical input densities use 5.755954.
                # Normalize to the explicitly requested 5.76, then apply the printed D.
                for key in ("nh", "nh1"):
                    grids[key] *= normalization * D
                # rho is comoving. Its physical A^3 factor is supplied by runtime a.
                if variant == "recomputed_background":
                    add_backgrounds(grids)
            np.testing.assert_allclose(grids["nh1"] / grids["nh"], original["nh1"] / original["nh"], rtol=5e-16)
            for key in ("T", "vpec", "rho"):
                np.testing.assert_array_equal(grids[key], original[key])
            work = OUT / variant / f"row{row}"
            inp, rt = work / "inputs", work / "rt"
            inp.mkdir(parents=True, exist_ok=True)
            rt.mkdir(exist_ok=True)
            for key, value in grids.items():
                if not np.all(np.isfinite(value)):
                    raise ValueError(f"nonfinite {variant} {row} {key}")
                np.savetxt(inp / f"los{los}_{key}.txt", value, fmt="%.18e")
            cmd = [str(ROOT / "rte1d_opt"), "1", "1", "1", "0", f"{float(np.float32(-26.99)):.17g}",
                   str(zrun), str(ZQ), str(inp), f"los{los}", str(rt), str(NCELLS), "1", str(los),
                   "1", str(EXTRA), "1", str(BOX), f"{scale:.17g}"]
            jobs.append(dict(variant=variant, target_row=row, global_input_los=los,
                             command=cmd, log=str(work / "run.log"),
                             output=str(rt / f"z6.14LOS{los}tq1mag-26.99_tep_1Myr1fd.txt")))
    config = dict(target=str(TARGET), target_sha256=hashlib.sha256(TARGET.read_bytes()).hexdigest(),
                  z_sim=ZS, z_qso=ZQ, input_density_conversion_redshift=ZINPUT, m1450=-26.99,
                  magnitude_as_float32=float(np.float32(-26.99)), tq_Myr=1,
                  density_factor_from_requested_zsim=D, physical_length_factor=A,
                  comoving_length_factor=C, input_density_normalization_to_zsim=normalization,
                  total_density_factor_applied_to_existing_inputs=normalization * D,
                  baseline_runtime_redshift=ZOLD, baseline_boxsize=BOX, scaled_boxsize=BOX*C,
                  initial_xHI_T_vpec_unchanged=True, comoving_rho_unchanged=True,
                  backgrounds="Both retained source rates and rates recalculated at scaled density are tested",
                  retained_legacy_settings=dict(legacy_density=1, legacy_alphahe2_bug=1, redward_pixels=EXTRA,
                                               X=.76, spectral_h=.678, spectral_omega_m=.308),
                  scope="Four complete RT sightlines, 891 cells each; not a full 945-row regeneration",
                  method_source="https://arxiv.org/html/2401.10328v1#S3.SS1")
    # Check physical geometry against the source convention, including both RT and spectra.
    np.testing.assert_allclose((C / (1+ZQ)) / (1/(1+ZS)), A, rtol=5e-16)
    write_json(OUT / "configuration.json", config)
    write_json(OUT / "manifest.json", jobs)
    return jobs


def run_job(job):
    start = time.monotonic()
    with open(job["log"], "w") as log:
        subprocess.run(job["command"], stdout=log, stderr=subprocess.STDOUT, check=True)
    result = dict(job, elapsed_seconds=time.monotonic()-start)
    print(f"Finished {job['variant']} row {job['target_row']}: {result['elapsed_seconds']:.1f}s", flush=True)
    return result


def velocity_grid(z, box):
    a = 1/(1+z)
    wmax = a * 3.2407789e-18 * 3.085678e21 / 1e5 * np.sqrt(.308/a**3 + .692) * box/.678
    # Positive toward the forest; bin EXTRA is the zero-velocity anchor of the C convolution.
    return (np.arange(NCELLS)-EXTRA) * wmax/NCELLS


def flux_metrics(delta):
    return dict(mean_absolute=float(np.abs(delta).mean()), max_absolute=float(np.abs(delta).max()),
                rms=float(np.sqrt(np.mean(delta**2))), mean_signed=float(delta.mean()))


def rebin_flux(velocity, spectra, edges):
    """Average piecewise-linear flux over each bin without extrapolating."""
    if edges[0] < velocity[0] or edges[-1] > velocity[-1]:
        raise ValueError("Rebinning edges exceed the available velocity range")
    result = []
    for spectrum in spectra:
        bins = []
        for left, right in zip(edges[:-1], edges[1:]):
            nodes = np.r_[left, velocity[(velocity>left)&(velocity<right)], right]
            values = np.interp(nodes, velocity, spectrum)
            bins.append(float(np.sum(np.diff(nodes)*(values[:-1]+values[1:])/2)/(right-left)))
        result.append(bins)
    return np.array(result)


def analyze():
    jobs = json.loads((OUT / "manifest.json").read_text())
    rows = np.array([item[0] for item in SAMPLES])
    reference = np.load(TARGET, mmap_mode="r")[rows]
    vref, vnew = velocity_grid(ZOLD, BOX), velocity_grid(ZQ, BOX*C)
    common = np.arange(max(vref[0], vnew[0]), min(vref[-1], vnew[-1]), 10.)
    report = dict(config=json.loads((OUT / "configuration.json").read_text()),
                  velocity_bin_kms_reference=float(vref[1]-vref[0]),
                  velocity_bin_kms_rescaled=float(vnew[1]-vnew[0]),
                  common_velocity_definition="Positive blueward of Ly-alpha; exact C Hubble bin spacing, 256-bin offset removed",
                  common_grid_step_kms=10., variants={}, baseline_regression={})
    arrays = dict(target_rows=rows, global_input_los=np.array([item[1] for item in SAMPLES]),
                  reference=reference, reference_velocity_kms=vref, rescaled_velocity_kms=vnew,
                  common_velocity_kms=common)
    ref_interp = np.array([np.interp(common, vref, np.exp(-r[:, 6])) for r in reference])
    arrays["reference_flux_common"] = ref_interp
    edges = np.arange(0., 4801., 600.)
    arrays["blueward_bin_edges_kms"] = edges
    arrays["reference_flux_600kms"] = rebin_flux(vref, np.exp(-reference[:,:,6]), edges)
    for variant in ("fixed_background", "recomputed_background", "baseline_regression"):
        selected = [job for job in jobs if job["variant"] == variant]
        results, candidate = [], []
        for job in selected:
            data = np.loadtxt(job["output"])
            row_index = list(rows).index(job["target_row"])
            metrics = comparison_metrics(reference[row_index], data)
            inp = Path(job["command"][8])
            rho, vpec = [np.loadtxt(inp/f"los{job['global_input_los']}_{key}.txt") for key in ("rho", "vpec")]
            z, box = (ZOLD, BOX) if variant == "baseline_regression" else (ZQ, BOX*C)
            tau = optical_depth(data, rho, vpec, z, np.arange(NCELLS), box, EXTRA)
            spectral_error = float(np.max(np.abs(np.exp(-tau)-np.exp(-data[:, 6]))))
            if spectral_error > 1e-6:
                raise ValueError(f"Independent spectral check failed: {spectral_error}")
            if variant == "baseline_regression" and metrics["max_absolute_flux_error"] > 1e-7:
                raise ValueError("Default length scale changed baseline")
            results.append(dict(target_row=job["target_row"], global_input_los=job["global_input_los"],
                                same_cell_metrics=metrics, independent_spectral_max_flux_error=spectral_error))
            candidate.append(data)
        if variant == "baseline_regression":
            report[variant] = results[0]
            continue
        candidate = np.array(candidate)
        arrays[variant] = candidate
        np.save(OUT / f"{variant}_sample_los.npy", candidate)
        flux = np.array([np.interp(common, vnew, np.exp(-r[:, 6])) for r in candidate])
        arrays[f"{variant}_flux_common"] = flux
        arrays[f"{variant}_flux_600kms"] = rebin_flux(vnew, np.exp(-candidate[:,:,6]), edges)
        delta = flux-ref_interp
        regions = {"whole_overlap": np.ones(len(common), dtype=bool),
                   "blueward_0_5000_kms": (common>=0)&(common<=5000),
                   "redward_minus5000_0_kms": (common>=-5000)&(common<0)}
        report["variants"][variant] = dict(per_sightline=results,
            aggregate_same_cell_flux=flux_metrics(np.exp(-candidate[:, :, 6])-np.exp(-reference[:, :, 6])),
            common_velocity_flux={key: flux_metrics(delta[:,mask]) for key,mask in regions.items()},
            common_velocity_mean_profile={key: flux_metrics(delta.mean(axis=0)[mask]) for key,mask in regions.items()},
            blueward_0_4800_flux_600kms=flux_metrics(arrays[f"{variant}_flux_600kms"]-arrays["reference_flux_600kms"]),
            blueward_0_4800_mean_profile_600kms=flux_metrics((arrays[f"{variant}_flux_600kms"]-arrays["reference_flux_600kms"]).mean(axis=0)),
            final_xHI_median_ratio=float(np.median(candidate[:, :, 1]/reference[:, :, 1])),
            final_T_median_ratio=float(np.median(candidate[:, :, 2]/reference[:, :, 2])))
        for k, result in enumerate(results):
            result["blueward_common_velocity_flux"] = flux_metrics(delta[k, regions["blueward_0_5000_kms"]])
    report["background_sensitivity_blueward"] = flux_metrics(
        (arrays["recomputed_background_flux_common"]-arrays["fixed_background_flux_common"])[:,(common>=0)&(common<=5000)])
    np.savez_compressed(OUT/"comparison_arrays.npz", **arrays)
    write_json(OUT/"comparison.json", report)
    plot(arrays)
    print(json.dumps({"factors": [D, A, C], "variants": {
        k: v["common_velocity_flux"] for k,v in report["variants"].items()},
        "background_sensitivity_blueward": report["background_sensitivity_blueward"]}, indent=2), flush=True)


def plot(arrays):
    os.environ.setdefault("MPLCONFIGDIR", str(OUT / ".matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(OUT / ".cache"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(4, 2, figsize=(13, 10), sharex=True, sharey="col", constrained_layout=True)
    variants = [("fixed_background", "Section 3.1; original background", "#d95f02"),
                ("recomputed_background", "Section 3.1; recalculated background", "#1b9e77")]
    for k, row in enumerate(arrays["target_rows"]):
        axes[k, 0].plot(arrays["reference_velocity_kms"], np.exp(-arrays["reference"][k,:,6]),
                        color="black", lw=1.1, label="Downloaded file")
        for key, label, color in variants:
            axes[k, 0].plot(arrays["rescaled_velocity_kms"], np.exp(-arrays[key][k,:,6]),
                            color=color, lw=.9, alpha=.9, label=label)
            axes[k, 1].plot(arrays["common_velocity_kms"], arrays[key+"_flux_common"][k]-arrays["reference_flux_common"][k],
                            color=color, lw=.9)
        axes[k, 0].set_ylabel(f"Row {row}\nTransmitted flux")
        axes[k, 1].set_ylabel("Rescaled − downloaded")
        for ax in axes[k]:
            ax.grid(alpha=.15)
            ax.set_xlim(-500, 5000)
            ax.axvline(0, color="gray", lw=.7, ls=":")
        axes[k, 1].axhline(0, color="black", lw=.7)
    axes[0, 0].legend(fontsize=8, loc="upper right")
    axes[0, 0].set_ylim(-.03, 1.03)
    for ax in axes[-1]:
        ax.set_xlabel("Velocity blueward of Lyα (km/s)")
    fig.suptitle("Literal Section 3.1 rescaling: zsim = 5.76 → zQSO = 6.14, M1450 = −26.99, tQ = 1 Myr")
    for suffix in ("png", "pdf"):
        fig.savefig(OUT/f"lya_comparison.{suffix}", dpi=160)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(8, 4.5), constrained_layout=True)
    centers = (arrays["blueward_bin_edges_kms"][:-1]+arrays["blueward_bin_edges_kms"][1:])/2
    ax.plot(centers, arrays["reference_flux_600kms"].mean(axis=0), "o-", color="black", label="Downloaded file")
    for key, label, color in variants:
        ax.plot(centers, arrays[key+"_flux_600kms"].mean(axis=0), "o-", color=color, label=label)
    ax.set(xlabel="Velocity blueward of Lyα (km/s)", ylabel="Mean transmitted flux",
           title="Four matched sightlines; flux averaged in 600 km/s bins", ylim=(-.02, 1.02))
    ax.grid(alpha=.2)
    ax.legend(fontsize=9)
    for suffix in ("png", "pdf"):
        fig.savefig(OUT/f"mean_lya_600kms.{suffix}", dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analyze-only", action="store_true")
    args = parser.parse_args()
    if not args.analyze_only:
        jobs = prepare()
        completed = []
        with ThreadPoolExecutor(max_workers=3) as pool:
            for future in as_completed([pool.submit(run_job, job) for job in jobs]):
                completed.append(future.result())
                write_json(OUT/"completed_jobs.json", completed)
    analyze()


if __name__ == "__main__":
    main()
