#!/usr/bin/env python3
"""
Compute and plot, for each sample:
  - Posterior predictive mean       (from combining toy histograms)
  - Posterior standard deviation    (from combining toy histograms)
  - Asimov prediction                (read directly from file)
  - Statistical uncertainty          (bin errors on the Asimov histogram)

Statistical uncertainty and posterior standard deviation share a colour
scale and colour map so they can be directly compared.

Usage:
    python plot_toy_mean_stddev_2D.py Posteriorpredictive_out_2D.root
    python plot_toy_mean_stddev_2D.py Posteriorpredictive_out_2D.root --sample FD_BeamFD_FHC_numu
"""

import argparse
import os
import re
from pathlib import Path
from collections import defaultdict

import numpy as np
import uproot
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# ---------------------------------------------------------------------------
# Thesis-quality plot styling
# ---------------------------------------------------------------------------
mpl.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
    "font.size": 12,
    "axes.labelsize": 13,
    "axes.titlesize": 13,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "legend.fontsize": 11,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})

CMAP_MEAN   = "viridis"
CMAP_ASIMOV = "viridis"
CMAP_UNC    = "magma"   # shared between posterior stddev and stat error

TOY_PATTERN = re.compile(
    r"^(?P<ND_or_FD>ND|FD)_(?P<sample>.+)_TrueNeutrinoEnergy_vs_RecoNeutrinoEnergy_posterior_toy_(?P<toy>\d+)$"
)


def find_toy_histograms(f):
    """Group toy histogram paths by sample key: dict[sample_key] -> dict[toy_idx] -> path"""
    toys_by_sample = defaultdict(dict)
    for full_key in f.keys(recursive=True, cycle=False):
        name = full_key.split("/")[-1]
        m = TOY_PATTERN.match(name)
        if not m:
            continue
        sample_key = f"{m.group('ND_or_FD')}_{m.group('sample')}"
        toys_by_sample[sample_key][int(m.group("toy"))] = full_key
    return toys_by_sample


def stack_toys(f, toy_paths):
    """Read all toys for a sample and stack into a 3D array (n_toys, nx, ny)."""
    toy_indices = sorted(toy_paths.keys())
    arrays = []
    xedges = yedges = None

    for toy_idx in toy_indices:
        h = f[toy_paths[toy_idx]]
        values, xe, ye = h.to_numpy()
        if xedges is None:
            xedges, yedges = xe, ye
        elif not (np.allclose(xe, xedges) and np.allclose(ye, yedges)):
            raise ValueError(
                f"Toy {toy_idx} has different binning to earlier toys — "
                "cannot safely combine into mean/stddev."
            )
        arrays.append(values)

    stack = np.stack(arrays, axis=0)
    return stack, xedges, yedges


def get_asimov(f, ND_or_FD, cleanTitle):
    """
    Return (values, xedges, yedges, stat_errors) for the Asimov histogram,
    or None if not found.
    """
    path = (
        f"Asimov/{ND_or_FD}/2D/{cleanTitle}/"
        f"{ND_or_FD}_{cleanTitle}_TrueNeutrinoEnergy_vs_RecoNeutrinoEnergy_Asimov"
    )
    try:
        h = f[path]
    except KeyError:
        print(f"[WARN] Could not find Asimov histogram at {path}")
        return None

    values, xedges, yedges = h.to_numpy()
    stat_err = h.errors()  # uses Sumw2 if tracked, else sqrt(N)
    return values, xedges, yedges, stat_err


def plot_map(values, xedges, yedges, title, cmap, cbar_label, outname,
             log=False, vmin=None, vmax=None, mask_zero=None):
    fig, ax = plt.subplots(figsize=(6, 5))

    if mask_zero is None:
        mask_zero = values == 0
    plot_vals = np.where(mask_zero, np.nan, values)

    cmap_obj = mpl.colormaps[cmap].copy()
    cmap_obj.set_bad(color="white")

    norm = None
    if log:
        lo = vmin if vmin is not None else (
            np.nanmin(plot_vals[plot_vals > 0]) if np.any(plot_vals > 0) else 1e-6
        )
        hi = vmax if vmax is not None else np.nanmax(plot_vals)
        norm = LogNorm(vmin=lo, vmax=hi)

    mesh = ax.pcolormesh(
        xedges, yedges, plot_vals.T, cmap=cmap_obj, norm=norm,
        vmin=None if log else vmin, vmax=None if log else vmax,
        shading="auto",
    )
    cbar = fig.colorbar(mesh, ax=ax)
    cbar.set_label(cbar_label)

    ax.set_xlabel(r"True Neutrino Energy $E_\nu^{\mathrm{true}}$ [GeV]")
    ax.set_ylabel(r"Reconstructed Neutrino Energy $E_\nu^{\mathrm{reco}}$ [GeV]")
    ax.set_title(title)

    fig.tight_layout()
    fig.savefig(outname)
    plt.close(fig)
    print(f"[INFO] Saved {outname}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("infile", help="Posterior predictive ROOT file")
    parser.add_argument("--outdir", default=None,
                         help="Output directory (default: same directory as infile)")
    parser.add_argument("--sample", default=None,
                         help="Restrict to one sample key (default: all samples found)")
    parser.add_argument("--log", action="store_true", help="Use log colour scale")
    args = parser.parse_args()

    if args.outdir is None:
        args.outdir = str(Path(args.infile).resolve().parent)
    os.makedirs(args.outdir, exist_ok=True)

    f = uproot.open(args.infile)
    toys_by_sample = find_toy_histograms(f)

    if not toys_by_sample:
        print("[ERROR] No posterior toy histograms found in file.")
        return

    samples = [args.sample] if args.sample else sorted(toys_by_sample.keys())

    for sample_key in samples:
        if sample_key not in toys_by_sample:
            print(f"[WARN] Sample '{sample_key}' not found, skipping. "
                  f"Available: {sorted(toys_by_sample.keys())}")
            continue

        m = TOY_PATTERN.match(  # reuse pattern's groups by matching any one toy path's basename
            next(iter(toys_by_sample[sample_key].values())).split("/")[-1]
        )
        ND_or_FD  = m.group("ND_or_FD")
        cleanTitle = m.group("sample")

        toy_paths = toys_by_sample[sample_key]
        n_toys = len(toy_paths)
        print(f"[INFO] {sample_key}: combining {n_toys} toys...")

        stack, xedges, yedges = stack_toys(f, toy_paths)
        mean_vals   = np.mean(stack, axis=0)
        stddev_vals = np.std(stack, axis=0, ddof=0)

        sample_label = sample_key.replace("_", r"\_")
        zero_mask = mean_vals == 0

        # -----------------------------------------------------------------
        # Posterior mean
        # -----------------------------------------------------------------
        plot_map(mean_vals, xedges, yedges,
                  title=rf"{sample_label}: Posterior Predictive Mean ({n_toys} toys)",
                  cmap=CMAP_MEAN, cbar_label=r"Mean events / bin",
                  outname=os.path.join(args.outdir, f"{sample_key}_toymean.pdf"),
                  log=args.log, mask_zero=zero_mask)

        # -----------------------------------------------------------------
        # Asimov + statistical uncertainty
        # -----------------------------------------------------------------
        asimov = get_asimov(f, ND_or_FD, cleanTitle)

        stat_vals = None
        if asimov is not None:
            asimov_vals, axedges, ayedges, stat_err = asimov

            plot_map(asimov_vals, axedges, ayedges,
                      title=rf"{sample_label}: Asimov Prediction",
                      cmap=CMAP_ASIMOV, cbar_label=r"Events / bin",
                      outname=os.path.join(args.outdir, f"{sample_key}_asimov.pdf"),
                      log=args.log, mask_zero=(asimov_vals == 0))

            stat_vals = stat_err
        else:
            print(f"[WARN] Skipping Asimov/stat-error plots for {sample_key} "
                  "(no Asimov histogram found).")

        # -----------------------------------------------------------------
        # Posterior stddev + statistical error, shared colour scale
        # -----------------------------------------------------------------
        if stat_vals is not None:
            combined = np.concatenate([
                stddev_vals[~zero_mask].ravel(),
                stat_vals[~zero_mask].ravel(),
            ])
            shared_vmin = float(np.min(combined)) if combined.size else 0.0
            shared_vmax = float(np.max(combined)) if combined.size else 1.0

            plot_map(stddev_vals, xedges, yedges,
                      title=rf"{sample_label}: Posterior Standard Deviation ({n_toys} toys)",
                      cmap=CMAP_UNC, cbar_label=r"Uncertainty [events]",
                      outname=os.path.join(args.outdir, f"{sample_key}_toystddev.pdf"),
                      log=args.log, vmin=shared_vmin, vmax=shared_vmax, mask_zero=zero_mask)

            plot_map(stat_vals, xedges, yedges,
                      title=rf"{sample_label}: Statistical Uncertainty",
                      cmap=CMAP_UNC, cbar_label=r"Uncertainty [events]",
                      outname=os.path.join(args.outdir, f"{sample_key}_staterr.pdf"),
                      log=args.log, vmin=shared_vmin, vmax=shared_vmax, mask_zero=zero_mask)
        else:
            # Fall back to unscaled stddev plot if no Asimov/stat error available
            plot_map(stddev_vals, xedges, yedges,
                      title=rf"{sample_label}: Posterior Standard Deviation ({n_toys} toys)",
                      cmap=CMAP_UNC, cbar_label=r"$\sigma_{\mathrm{posterior}}$ [events]",
                      outname=os.path.join(args.outdir, f"{sample_key}_toystddev.pdf"),
                      log=args.log, mask_zero=zero_mask)


if __name__ == "__main__":
    main()