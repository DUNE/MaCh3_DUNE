#!/usr/bin/env python3
"""
Compare two run plans -- e.g. On-Axis vs On+Off-Axis ND running -- computed
directly from each run plan's raw Posteriorpredictive_out_2D.root file.
No MakePosteriorPredictiveSummary2D.C step needed for either file: this
reads the individual toy histograms and Asimov reference straight out of
the raw job output, same as plot_toy_mean_stddev_2D.py, then ratios the
results between the two files.

For every sample present in BOTH files, computes and plots:

    mean ratio      = mean_A / mean_B
    rel_post ratio  = (stddev_A/mean_A) / (stddev_B/mean_B)
    rel_stat ratio  = (staterr_A/mean_A) / (staterr_B/mean_B)

For all three: ratio < 1 means run plan A is smaller/tighter in that bin,
ratio > 1 means run plan B is. All three use a diverging, log-scaled
colour map centred on 1 so an improvement and a degradation of the same
factor get equal visual weight.

Usage:
    python plot_run_plan_relunc_ratio_2D.py \
        onaxis_out_2D.root      "On-Axis" \
        onoffaxis_out_2D.root   "On+Off-Axis" \
        --outdir thesis_plots
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
from matplotlib.colors import TwoSlopeNorm

# ---------------------------------------------------------------------------
# Thesis-quality plot styling (matches plot_toy_mean_stddev_2D.py)
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

CMAP_RATIO = "RdBu_r"

TOY_PATTERN = re.compile(
    r"^(?P<ND_or_FD>ND|FD)_(?P<sample>.+)_TrueNeutrinoEnergy_vs_RecoNeutrinoEnergy_posterior_toy_(?P<toy>\d+)$"
)


# ------------------------------------------------------------
# Reading helpers (same logic as plot_toy_mean_stddev_2D.py)
# ------------------------------------------------------------

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
    """Return (values, xedges, yedges, stat_errors) for the Asimov histogram, or None."""
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
    stat_err = h.errors()
    return values, xedges, yedges, stat_err


def compute_sample_quantities(f, sample_key):
    """
    Return dict with mean, stddev, stat_err, ND_or_FD, cleanTitle, x/y edges
    for one sample key, or None if the sample or its Asimov reference is
    missing.
    """
    toys_by_sample = find_toy_histograms(f)
    if sample_key not in toys_by_sample:
        return None

    toy_paths = toys_by_sample[sample_key]
    m = TOY_PATTERN.match(next(iter(toy_paths.values())).split("/")[-1])
    ND_or_FD = m.group("ND_or_FD")
    cleanTitle = m.group("sample")

    stack, xedges, yedges = stack_toys(f, toy_paths)
    mean_vals = np.mean(stack, axis=0)
    stddev_vals = np.std(stack, axis=0, ddof=0)

    asimov = get_asimov(f, ND_or_FD, cleanTitle)
    if asimov is None:
        print(f"[WARN] No Asimov reference for {sample_key}, skipping "
              "(needed for statistical uncertainty).")
        return None
    _, _, _, stat_err = asimov

    if stat_err.shape != mean_vals.shape:
        print(f"[WARN] Asimov/toy binning mismatch for {sample_key}, skipping.")
        return None

    return dict(mean=mean_vals, stddev=stddev_vals, stat_err=stat_err,
                xedges=xedges, yedges=yedges, n_toys=len(toy_paths))


# ------------------------------------------------------------
# Plotting
# ------------------------------------------------------------

def diverging_log_norm(values, mask, center=1.0, min_span=0.2):
    finite = values[mask & np.isfinite(values) & (values > 0)]
    if finite.size == 0:
        return None
    lo, hi = np.nanpercentile(finite, [2, 98])
    span = max(abs(np.log(lo / center)), abs(np.log(hi / center)), min_span)
    return TwoSlopeNorm(vmin=center * np.exp(-span), vcenter=center, vmax=center * np.exp(span))


def plot_ratio_map(ratio, xedges, yedges, title, subtitle, cbar_label, outname,
                    mask_zero, min_span=0.2):
    finite = ratio[~mask_zero & np.isfinite(ratio) & (ratio > 0)]
    if finite.size == 0:
        print(f"[WARN] All-NaN/zero ratio, skipping {outname}")
        return

    norm = diverging_log_norm(ratio, ~mask_zero, min_span=min_span)
    plot_vals = np.where(mask_zero, np.nan, ratio)

    cmap_obj = mpl.colormaps[CMAP_RATIO].copy()
    cmap_obj.set_bad(color="white")

    fig, ax = plt.subplots(figsize=(6.5, 5.3))
    mesh = ax.pcolormesh(xedges, yedges, plot_vals.T, cmap=cmap_obj, norm=norm, shading="auto")

    ax.set_xlabel(r"True Neutrino Energy $E_\nu^{\mathrm{true}}$ [GeV]")
    ax.set_ylabel(r"Reconstructed Neutrino Energy $E_\nu^{\mathrm{reco}}$ [GeV]")
    ax.set_title(title, pad=26)
    ax.text(0.5, 1.04, subtitle, transform=ax.transAxes, ha="center", va="bottom",
             fontsize=10, style="italic")

    cbar = fig.colorbar(mesh, ax=ax)
    cbar.set_label(cbar_label)

    fig.tight_layout()
    fig.savefig(outname)
    plt.close(fig)
    print(f"[INFO] Saved {outname}")


def plot_two_panel_ratio(nd_entry, fd_entry, suptitle, cbar_label, outname, min_span=0.2):
    """
    nd_entry / fd_entry: dict with 'ratio', 'xedges', 'yedges', 'mask_zero',
    or None. Builds one figure, ND panel on top, FD panel underneath,
    sharing a single diverging log colour scale and colorbar computed from
    the combined data of both panels.
    """
    panels = [(det, entry) for det, entry in (("ND", nd_entry), ("FD", fd_entry)) if entry is not None]
    if len(panels) < 2:
        print(f"[WARN] Need both ND and FD to build combined panel for '{outname}', skipping")
        return

    all_finite = []
    for _, entry in panels:
        vals = entry["ratio"]
        mask = ~entry["mask_zero"]
        all_finite.append(vals[mask & np.isfinite(vals) & (vals > 0)])
    all_finite = np.concatenate(all_finite) if all_finite else np.array([])
    if all_finite.size == 0:
        print(f"[WARN] No finite values for combined plot '{outname}', skipping")
        return

    lo, hi = np.nanpercentile(all_finite, [2, 98])
    span = max(abs(np.log(lo)), abs(np.log(hi)), min_span)
    norm = TwoSlopeNorm(vmin=np.exp(-span), vcenter=1.0, vmax=np.exp(span))

    cmap_obj = mpl.colormaps[CMAP_RATIO].copy()
    cmap_obj.set_bad(color="white")

    fig, axes = plt.subplots(2, 1, figsize=(6.5, 9.5), sharex=True)

    mesh = None
    for ax, (det, entry) in zip(axes, panels):
        plot_vals = np.where(entry["mask_zero"], np.nan, entry["ratio"])
        mesh = ax.pcolormesh(entry["xedges"], entry["yedges"], plot_vals.T,
                              cmap=cmap_obj, norm=norm, shading="auto")
        ax.set_ylabel(r"Reconstructed Neutrino Energy $E_\nu^{\mathrm{reco}}$ [GeV]")
        ax.set_title(det, fontsize=12, loc="left")

    axes[-1].set_xlabel(r"True Neutrino Energy $E_\nu^{\mathrm{true}}$ [GeV]")
    fig.suptitle(suptitle, fontsize=13)

    fig.subplots_adjust(right=0.86, hspace=0.15, top=0.93)
    cax = fig.add_axes([0.88, 0.11, 0.03, 0.78])
    cbar = fig.colorbar(mesh, cax=cax)
    cbar.set_label(cbar_label)

    fig.savefig(outname)
    plt.close(fig)
    print(f"[INFO] Saved {outname}")


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("file_a", help="raw Posteriorpredictive_out_2D.root for run plan A")
    parser.add_argument("label_a", help="short label for run plan A, e.g. 'On-Axis'")
    parser.add_argument("file_b", help="raw Posteriorpredictive_out_2D.root for run plan B")
    parser.add_argument("label_b", help="short label for run plan B, e.g. 'On+Off-Axis'")
    parser.add_argument("--outdir", default="thesis_plots", help="output directory for PDFs")
    parser.add_argument("--sample", default=None, help="Restrict to one sample key")
    args = parser.parse_args()

    for path in (args.file_a, args.file_b):
        if not Path(path).exists():
            print(f"[ERROR] File not found: {path}")
            return

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    f_a = uproot.open(args.file_a)
    f_b = uproot.open(args.file_b)

    samples_a = set(find_toy_histograms(f_a).keys())
    samples_b = set(find_toy_histograms(f_b).keys())
    common = sorted(samples_a & samples_b) if args.sample is None else [args.sample]

    for s in samples_b - samples_a:
        print(f"[WARN] '{s}' present in '{args.label_b}' but not '{args.label_a}', skipping")
    for s in samples_a - samples_b:
        print(f"[WARN] '{s}' present in '{args.label_a}' but not '{args.label_b}', skipping")

    n_plots = 0
    relpost_by_det = {"ND": {}, "FD": {}}  # cleanTitle -> entry dict, per detector

    for sample_key in common:
        print(f"[INFO] {sample_key}: computing ratio {args.label_a} / {args.label_b}...")

        qa = compute_sample_quantities(f_a, sample_key)
        qb = compute_sample_quantities(f_b, sample_key)
        if qa is None or qb is None:
            print(f"[WARN] Incomplete data for {sample_key} in one or both files, skipping")
            continue

        if qa["mean"].shape != qb["mean"].shape:
            print(f"[WARN] Binning mismatch for {sample_key} "
                  f"({qa['mean'].shape} vs {qb['mean'].shape}), skipping")
            continue

        xedges, yedges = qa["xedges"], qa["yedges"]
        zero_mask = (qa["mean"] == 0) | (qb["mean"] == 0)
        sample_label = sample_key.replace("_", r"\_")

        with np.errstate(divide="ignore", invalid="ignore"):
            rel_post_a = np.where(qa["mean"] > 0, qa["stddev"] / qa["mean"], np.nan)
            rel_post_b = np.where(qb["mean"] > 0, qb["stddev"] / qb["mean"], np.nan)
            rel_stat_a = np.where(qa["mean"] > 0, qa["stat_err"] / qa["mean"], np.nan)
            rel_stat_b = np.where(qb["mean"] > 0, qb["stat_err"] / qb["mean"], np.nan)

            mean_ratio     = np.where(qb["mean"] > 0, qa["mean"] / qb["mean"], np.nan)
            rel_post_ratio = np.where(rel_post_b > 0, rel_post_a / rel_post_b, np.nan)
            rel_stat_ratio = np.where(rel_stat_b > 0, rel_stat_a / rel_stat_b, np.nan)

        det, _, cleanTitle = sample_key.partition("_")
        if det in relpost_by_det:
            relpost_by_det[det][cleanTitle] = dict(
                ratio=rel_post_ratio, xedges=xedges, yedges=yedges, mask_zero=zero_mask)

        common_sub = rf"{args.label_a} / {args.label_b}"

        plot_ratio_map(
            mean_ratio, xedges, yedges,
            title=rf"{sample_label}: Predicted Event Rate Ratio",
            subtitle=rf"Posterior mean, {common_sub}",
            cbar_label=r"Mean$_A$ / Mean$_B$",
            outname=str(outdir / f"{sample_key}_mean_ratio.pdf"),
            mask_zero=zero_mask, min_span=0.1)
        n_plots += 1

        plot_ratio_map(
            rel_post_ratio, xedges, yedges,
            title=rf"{sample_label}: Relative Posterior Uncertainty Ratio",
            subtitle=rf"$(\sigma_{{\mathrm{{post}}}}/\mu)$, {common_sub}",
            cbar_label=r"$(\sigma_{\mathrm{post}}/\mu)_A \,/\, (\sigma_{\mathrm{post}}/\mu)_B$",
            outname=str(outdir / f"{sample_key}_relpost_ratio.pdf"),
            mask_zero=zero_mask)
        n_plots += 1

        plot_ratio_map(
            rel_stat_ratio, xedges, yedges,
            title=rf"{sample_label}: Relative Statistical Uncertainty Ratio",
            subtitle=rf"$(\sigma_{{\mathrm{{stat}}}}/\mu)$, {common_sub}",
            cbar_label=r"$(\sigma_{\mathrm{stat}}/\mu)_A \,/\, (\sigma_{\mathrm{stat}}/\mu)_B$",
            outname=str(outdir / f"{sample_key}_relstat_ratio.pdf"),
            mask_zero=zero_mask)
        n_plots += 1

    # -----------------------------------------------------------------
    # Combined ND-on-top / FD-underneath figure(s) for the relative
    # posterior uncertainty ratio. Pairs samples with the same cleanTitle
    # across detectors; if there's exactly one ND and one FD sample
    # overall (the common single-sample-per-detector case), pairs those
    # regardless of title, since that's almost always the intended
    # near/far comparison.
    # -----------------------------------------------------------------
    nd_titles = set(relpost_by_det["ND"].keys())
    fd_titles = set(relpost_by_det["FD"].keys())
    matched_titles = sorted(nd_titles & fd_titles)

    common_sub = rf"{args.label_a} / {args.label_b}"
    cbar_label = r"$(\sigma_{\mathrm{post}}/\mu)_A \,/\, (\sigma_{\mathrm{post}}/\mu)_B$"

    if matched_titles:
        for cleanTitle in matched_titles:
            title_label = cleanTitle.replace("_", r"\_")
            plot_two_panel_ratio(
                relpost_by_det["ND"][cleanTitle], relpost_by_det["FD"][cleanTitle],
                suptitle=rf"{title_label}: Relative Posterior Uncertainty Ratio ({common_sub})",
                cbar_label=cbar_label,
                outname=str(outdir / f"{cleanTitle}_ND_FD_relpost_ratio.pdf"))
    elif len(nd_titles) == 1 and len(fd_titles) == 1:
        nd_title = next(iter(nd_titles))
        fd_title = next(iter(fd_titles))
        plot_two_panel_ratio(
            relpost_by_det["ND"][nd_title], relpost_by_det["FD"][fd_title],
            suptitle=rf"Relative Posterior Uncertainty Ratio ({common_sub})",
            cbar_label=cbar_label,
            outname=str(outdir / "ND_FD_relpost_ratio_combined.pdf"))
    elif nd_titles and fd_titles:
        print("[WARN] Multiple ND and FD samples with no matching titles -- "
              "not auto-combining into a single ND/FD figure. Use --sample to "
              "narrow to one ND and one FD sample if you want the combined plot.")

    if n_plots == 0:
        print("[WARN] No comparison plots were made -- check that both files contain "
              "posterior toy histograms and matching Asimov references for the same samples.")
    else:
        print(f"[INFO] Saved {n_plots} plot(s) to {outdir}/")


if __name__ == "__main__":
    main()