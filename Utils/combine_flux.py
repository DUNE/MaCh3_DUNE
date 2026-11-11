#!/usr/bin/env python3
"""
Combine 18 separate per-position MaCh3 samples (each with a 1D
TrueNeutrinoEnergy binning) into a single 2D (energy x off-axis position)
+1sigma/nominal ratio plot per flux systematic, for thesis use.

Assumed file structure (per MaCh3 sample, one file per sample OR all
samples' outputs living in one file under a directory per sample -- see
--file-per-sample / --infile below):

    <SystName>/<SampleDir>/Variation_N     (TH1D, 1D in energy)

where <SampleDir> is one of the 18 sample names below, and Variation_N
follows the same -3,-1,0,+1,+3 sigma convention as before (nominal=2,
+1sigma=3 by default).

Two supported layouts:

1) ONE file, containing every sample's directory tree:
    python combine_position_ratios.py \
        --infile all_samples_template.root \
        --outdir thesis_plots

2) ONE file PER sample (e.g. produced by 18 separate MaCh3 runs):
    python combine_position_ratios.py \
        --file-pattern "eventrate_{sample}_template.root" \
        --outdir thesis_plots

The --positions list below gives the y-axis ordering/value for each
sample; edit if your sample names or positions differ.
"""

import argparse
import os
import re

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import uproot


# Ordered (sample_name, off_axis_position_in_metres) pairs.
# Edit this list if your sample naming or position values differ.
SAMPLE_POSITIONS = [
    ("OffAxis0m_numuCC_numode", 0.0),
    ("OffAxis1.75m_numuCC_numode", 1.75),
    ("OffAxis2m_numuCC_numode", 2.0),
    ("OffAxis4m_numuCC_numode", 4.0),
    ("OffAxis5.75m_numuCC_numode", 5.75),
    ("OffAxis8m_numuCC_numode", 8.0),
    ("OffAxis9.75m_numuCC_numode", 9.75),
    ("OffAxis12m_numuCC_numode", 12.0),
    ("OffAxis13.75m_numuCC_numode", 13.75),
    ("OffAxis16m_numuCC_numode", 16.0),
    ("OffAxis17.75m_numuCC_numode", 17.75),
    ("OffAxis20m_numuCC_numode", 20.0),
    ("OffAxis21.75m_numuCC_numode", 21.75),
    ("OffAxis24m_numuCC_numode", 24.0),
    ("OffAxis25.75m_numuCC_numode", 25.75),
    ("OffAxis26.25m_numuCC_numode", 26.25),
    ("OffAxis28.25m_numuCC_numode", 28.25),
    ("OffAxis28.5m_numuCC_numode", 28.5),
]


def strip_cycle(key):
    return re.sub(r";\d+$", "", key)


def discover_systematics(f, sample_names):
    """Find top-level directories that contain at least one known sample dir."""
    systematics = []
    for key in f.keys(recursive=False):
        name = strip_cycle(key)
        obj = f[key]
        if not isinstance(obj, uproot.reading.ReadOnlyDirectory):
            continue
        for sample_name in sample_names:
            try:
                f[f"{name}/{sample_name}"]
                systematics.append(name)
                break
            except Exception:
                continue
    return sorted(systematics)


def th1_to_arrays(th1):
    values, edges = th1.to_numpy()
    return values, edges


def stack_ratio_curves(per_sample):
    """Given a list of (position, edges, centers, ratio_curve), interpolate
    each onto the common (finest) energy grid and stack into a 2D array,
    returning (ratio_grid, x_edges, y_edges)."""
    if not per_sample:
        return None, None, None

    all_edges = [edges for _, edges, _, _ in per_sample]
    x_edges = max(all_edges, key=len)
    x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])

    ratio_rows = []
    used_positions = []
    for position, edges, centers, ratio_curve in per_sample:
        if len(edges) == len(x_edges) and np.allclose(edges, x_edges):
            resampled = ratio_curve
        else:
            resampled = np.interp(x_centers, centers, ratio_curve,
                                   left=np.nan, right=np.nan)
        ratio_rows.append(resampled)
        used_positions.append(position)

    ratio_grid = np.array(ratio_rows)

    positions = np.array(used_positions, dtype=float)
    if len(positions) == 1:
        y_edges = np.array([positions[0] - 0.5, positions[0] + 0.5])
    else:
        mids = (positions[:-1] + positions[1:]) / 2.0
        first_edge = positions[0] - (mids[0] - positions[0])
        last_edge = positions[-1] + (positions[-1] - mids[-1])
        y_edges = np.concatenate([[first_edge], mids, [last_edge]])

    return ratio_grid, x_edges, y_edges


def build_2d_grid(f, syst_name, nominal_key, sigma_key):
    """Stack each sample's 1D nominal/+1sigma histograms into a 2D ratio grid.

    Samples are allowed to have different energy binning (common in this
    dataset -- some positions use a coarser array). Rather than requiring
    identical bin edges, each sample's own ratio curve (sigma/nominal vs.
    its own bin centers) is linearly interpolated onto a single common
    energy grid (the finest -- i.e. most-binned -- of the samples found)
    before stacking into rows.
    """
    per_sample = []  # list of (position, edges, centers, ratio_curve)

    for sample_name, position in SAMPLE_POSITIONS:
        base = f"{syst_name}/{sample_name}"
        try:
            nom_h = f[f"{base}/{nominal_key}"]
            sig_h = f[f"{base}/{sigma_key}"]
        except Exception:
            print(f"  ...missing {sample_name} for {syst_name}, skipping this position")
            continue

        nom_vals, edges = th1_to_arrays(nom_h)
        sig_vals, _ = th1_to_arrays(sig_h)

        with np.errstate(divide="ignore", invalid="ignore"):
            ratio_curve = np.where(nom_vals != 0, sig_vals / nom_vals, np.nan)

        centers = 0.5 * (edges[:-1] + edges[1:])
        per_sample.append((position, edges, centers, ratio_curve))

    return stack_ratio_curves(per_sample)


def plot_ratio(name, ratio_grid, x_edges, y_edges, outdir,
               x_label, y_label, vmin, vmax, cmap, dpi):
    fig, ax = plt.subplots(figsize=(7, 5.5))

    norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=1.0, vmax=vmax)
    # ratio_grid has shape (n_positions, n_energy_bins) -> matches (y, x) for pcolormesh
    mesh = ax.pcolormesh(x_edges, y_edges, ratio_grid, cmap=cmap, norm=norm, shading="auto")

    ax.set_title(f"{name} (+1$\\sigma$ / Nominal)")
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    fig.colorbar(mesh, ax=ax)
    fig.tight_layout()

    safe_name = re.sub(r"[^\w\-.]", "_", name)
    outpath_pdf = os.path.join(outdir, f"{safe_name}_p1sigma_ratio.pdf")
    outpath_png = os.path.join(outdir, f"{safe_name}_p1sigma_ratio.png")
    fig.savefig(outpath_pdf)
    fig.savefig(outpath_png, dpi=dpi)
    plt.close(fig)

    print(f"Wrote {outpath_pdf} and {outpath_png}")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--infile", help="Single ROOT file containing all samples' directories")
    parser.add_argument("--file-pattern",
                         help="Pattern with a {sample} placeholder for one-file-per-sample "
                              "layout, e.g. 'eventrate_{sample}_template.root'")
    parser.add_argument("--outdir", default="thesis_plots")
    parser.add_argument("--nominal-variation", type=int, default=2)
    parser.add_argument("--sigma-variation", type=int, default=3)
    parser.add_argument("--x-label", default="Neutrino Energy [GeV]")
    parser.add_argument("--y-label", default="Off-axis [m]")
    parser.add_argument("--vmin", type=float, default=0.98)
    parser.add_argument("--vmax", type=float, default=1.02)
    parser.add_argument("--cmap", default="viridis")
    parser.add_argument("--dpi", type=int, default=200)
    parser.add_argument("--only", nargs="*", default=None,
                         help="Optional list of systematic names to restrict to")
    args = parser.parse_args()

    if not args.infile and not args.file_pattern:
        parser.error("Must supply either --infile or --file-pattern")

    os.makedirs(args.outdir, exist_ok=True)

    nominal_key = f"Variation_{args.nominal_variation}"
    sigma_key = f"Variation_{args.sigma_variation}"
    sample_names = [s for s, _ in SAMPLE_POSITIONS]

    if args.infile:
        with uproot.open(args.infile) as f:
            systematics = discover_systematics(f, sample_names)
            if args.only:
                systematics = [s for s in systematics if s in args.only]

            if not systematics:
                print("No systematics found matching any known sample directories. "
                      "Top-level keys in file:")
                for k in f.keys(recursive=False):
                    print(f"  {strip_cycle(k)}")
                return

            print(f"Found {len(systematics)} flux systematics: {systematics}")

            for name in systematics:
                print(f"Processing {name}...")
                ratio_grid, x_edges, y_edges = build_2d_grid(
                    f, name, nominal_key, sigma_key
                )
                if ratio_grid is None:
                    print(f"  ...no usable positions found for {name}, skipping")
                    continue
                plot_ratio(
                    name, ratio_grid, x_edges, y_edges, args.outdir,
                    args.x_label, args.y_label, args.vmin, args.vmax, args.cmap, args.dpi,
                )
    else:
        # one-file-per-sample layout: open each sample's file, look for
        # every systematic directory inside it, aggregate by systematic name
        print("NOTE: one-file-per-sample layout requires opening all 18 files "
              "to discover systematics — this may take a moment.")
        per_sample_files = {}
        for sample_name, _ in SAMPLE_POSITIONS:
            path = args.file_pattern.format(sample=sample_name)
            if os.path.exists(path):
                per_sample_files[sample_name] = uproot.open(path)
            else:
                print(f"  ...file not found for {sample_name}: {path}")

        if not per_sample_files:
            print("No sample files found matching --file-pattern.")
            return

        # collect the union of systematic names across all opened files
        all_systs = set()
        for f in per_sample_files.values():
            for key in f.keys(recursive=False):
                obj = f[key]
                if isinstance(obj, uproot.reading.ReadOnlyDirectory):
                    all_systs.add(strip_cycle(key))
        systematics = sorted(all_systs)

        if args.only:
            systematics = [s for s in systematics if s in args.only]

        print(f"Found {len(systematics)} flux systematics: {systematics}")

        for name in systematics:
            print(f"Processing {name}...")
            per_sample = []
            for sample_name, position in SAMPLE_POSITIONS:
                f = per_sample_files.get(sample_name)
                if f is None:
                    continue
                try:
                    nom_h = f[f"{name}/{sample_name}/{nominal_key}"]
                    sig_h = f[f"{name}/{sample_name}/{sigma_key}"]
                except Exception:
                    continue
                nom_vals, edges = th1_to_arrays(nom_h)
                sig_vals, _ = th1_to_arrays(sig_h)
                with np.errstate(divide="ignore", invalid="ignore"):
                    ratio_curve = np.where(nom_vals != 0, sig_vals / nom_vals, np.nan)
                centers = 0.5 * (edges[:-1] + edges[1:])
                per_sample.append((position, edges, centers, ratio_curve))

            ratio_grid, x_edges, y_edges = stack_ratio_curves(per_sample)
            if ratio_grid is None:
                print(f"  ...no usable positions found for {name}, skipping")
                continue

            plot_ratio(
                name, ratio_grid, x_edges, y_edges, args.outdir,
                args.x_label, args.y_label, args.vmin, args.vmax, args.cmap, args.dpi,
            )

        for f in per_sample_files.values():
            f.close()


if __name__ == "__main__":
    main()