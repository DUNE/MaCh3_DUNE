#!/usr/bin/env python3
"""
compare_template_params.py
============================

Compares template (systematic) parameters across run plans, relative to
a reference run plan (default: "On Axis only"), reading
mcmc_diagnostics_*.root files produced by the MaCh3 diagnostics generator.

Real structure (per the C++ generator):
    xsec_eventrate_histo   TH2D, binned in (xsec_var1, xsec_var2)
                            -- here: (TrueNeutrinoEnergy, Enubias)
    xsec_param_mean        TH2D, SAME (q0, q3) binning as above. Each
                            template parameter's posterior mean is written
                            into whichever 2D bin its kinematic cut maps
                            to; all other bins stay at 0.
    xsec_param_stddev      TH2D, same idea, posterior std dev.

Produces, per histogram (mean, stddev, eventrate):
    - a grid of 2D ratio heatmaps (one per run plan) vs the reference
And, for mean/stddev specifically, also:
    - a 1D "per active parameter" ratio line plot, recovered by flattening
      just the nonzero bins (in a consistent order defined by the
      reference file) rather than assuming any particular 1D layout.

Usage
-----
1. Edit `run_plan_files` and `REFERENCE_LABEL` below.
2. python compare_template_params.py
"""

import logging
import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.backends.backend_pdf import PdfPages

log = logging.getLogger(__name__)

# ============================================================================
# EDIT THIS
# ============================================================================

run_plan_files = {
    "On Axis only":
        "/home/abipeake/scratch/July2026/OnAxis_fulllsysts_adaptive_post/mcmc_diagnostics_2026827_162830.root",
    "50% On Axis + 50% Off Axis":
        "/home/abipeake/scratch/July2026/Asimov_allsystematics_onplusoffaxis_postadaptive/mcmc_diagnostics_2026825_43951.root",
    "75% On Axis + 25% Off Axis":
        "/home/abipeake/scratch/July2026/Runplan_75onaxis_25offaxis_post/mcmc_diagnostics_2026824_51125.root",
    "25% On Axis + 75% Off Axis":
        "/home/abipeake/scratch/July2026/Runplan_25onaxis_75offaxis_post/mcmc_diagnostics_2026824_64824.root",
    "Equal time at each position":
        "/home/abipeake/scratch/July2026/Runplan_equal_allposition_post/mcmc_diagnostics_2026824_52541.root",
    "50% On Axis + 50% Off Axis (12m)":
        "/home/abipeake/scratch/July2026/Runplan_50onaxis_50upto12mpost/mcmc_diagnostics_2026824_65745.root",
}

REFERENCE_LABEL = "On Axis only"

XLABEL = r"True Neutrino Energy (GeV)"   # xsec_var1
YLABEL = r"$E_\nu$ Bias (GeV)"           # xsec_var2

PARAMS_PER_PAGE = 30   # split the 1D ratio plot into pages of this many params
HIST_2D_SPECS = [
    ("xsec_param_mean",     "Posterior Mean"),
    ("xsec_param_stddev",   "Posterior Std Dev"),
    ("xsec_eventrate_histo","Event Rate"),
]

# Okabe-Ito colourblind-safe palette (verified via CVD simulation)
PALETTE = [
    "#0072B2", "#E69F00", "#009E73", "#D55E00",
    "#CC79A7", "#56B4E9", "#F0E442",
]

THESIS_RC = {
    "font.family":         "serif",
    "font.serif":          ["STIXTwo Text", "DejaVu Serif", "Times New Roman"],
    "mathtext.fontset":    "stix",
    "font.size":           11,
    "axes.titlesize":      13,
    "axes.labelsize":      12,
    "xtick.labelsize":     9,
    "ytick.labelsize":     10,
    "legend.fontsize":     10,
    "pdf.fonttype":        42,
    "ps.fonttype":         42,
    "figure.dpi":          150,
    "xtick.direction":     "in",
    "ytick.direction":     "in",
    "ytick.right":         True,
    "axes.grid":           True,
    "grid.alpha":          0.3,
    "grid.linewidth":      0.6,
    "axes.linewidth":      0.8,
}


# ============================================================================
# Extraction
# ============================================================================

def extract_2d_hist(f, hist_name, verbose=True):
    """Return (values, xedges, yedges) for a 2D histogram."""
    h = f[hist_name]
    vals = h.values()
    xedges = h.axes[0].edges()
    yedges = h.axes[1].edges()
    if verbose:
        nonzero = np.count_nonzero(vals)
        print(f"    {hist_name}: shape={vals.shape}, {nonzero} nonzero bins, "
              f"x range=[{xedges[0]:.3g}, {xedges[-1]:.3g}], "
              f"y range=[{yedges[0]:.3g}, {yedges[-1]:.3g}]")
    return vals, xedges, yedges


def load_run_plan_hists(path):
    """Load all HIST_2D_SPECS histograms for one file. Returns
    {hist_name: (values, xedges, yedges)}, skipping any that fail."""
    f = uproot.open(path)
    out = {}
    for hist_name, _ in HIST_2D_SPECS:
        try:
            out[hist_name] = extract_2d_hist(f, hist_name)
        except Exception as exc:
            print(f"    could not load {hist_name}: {exc}")
    return out


# ============================================================================
# 2D ratio heatmaps
# ============================================================================

def _sanitize_filename(name):
    """Make a run-plan label safe to use in a filename."""
    import re
    return re.sub(r"[^A-Za-z0-9_.-]", "_", name)


def plot_2d_ratio(ax, ratio, xedges, yedges, title, norm, fontsize=13, cmap="coolwarm"):
    """One run-plan/reference ratio heatmap on `ax`, using a pre-computed
    shared `norm` so every page in the output is on the same colour scale
    and genuinely comparable to the others. Diverging colormap: coolwarm
    (cool blue = below reference, warm red = above), CVD-checked -- worst-
    case adjacent-step distance under simulated deuteranomaly/protanomaly
    is 14.7, comfortably above the ~10 threshold where colours start being
    genuinely risky to confuse."""
    mesh = ax.pcolormesh(xedges, yedges, ratio.T, shading="auto", cmap=cmap, norm=norm)
    ax.set_title(title, fontsize=fontsize, pad=8)
    return mesh


def _compute_shared_norm(ratios, vmax_pct=98):
    """Symmetric-around-1 TwoSlopeNorm computed from the combined finite
    values of every ratio array, so all pages/panels share one colour
    scale instead of each being scaled to its own range independently."""
    all_finite = np.concatenate([r[np.isfinite(r)] for r in ratios if np.isfinite(r).any()])
    if all_finite.size == 0:
        return mcolors.TwoSlopeNorm(vmin=0.5, vcenter=1.0, vmax=1.5)
    spread = max(np.percentile(np.abs(all_finite - 1.0), vmax_pct), 1e-6)
    return mcolors.TwoSlopeNorm(vmin=1 - spread, vcenter=1.0, vmax=1 + spread)


def plot_all_2d_ratios(results, hist_name, hist_label, reference_label,
                        xlabel=XLABEL, ylabel=YLABEL, save_prefix=None):
    """
    results: {label: {hist_name: (vals, xe, ye), ...}, ...}

    Writes one standalone PDF per non-reference run plan:
        "{save_prefix}_{run_plan}.pdf"
    All files share one colour scale (computed across every run plan
    together first) so they're genuinely comparable to each other even
    though each lives in its own file.
    """
    if save_prefix is None:
        save_prefix = f"template_{hist_name}_ratio"
    if reference_label not in results or hist_name not in results[reference_label]:
        print(f"ERROR: reference '{reference_label}' missing {hist_name}.")
        return
    ref_vals, ref_xe, ref_ye = results[reference_label][hist_name]

    others = [(label, r[hist_name]) for label, r in results.items()
              if label != reference_label and hist_name in r]
    if not others:
        print(f"No other run plans have {hist_name} loaded.")
        return

    # Compute every ratio up front so the shared colour scale reflects all
    # of them together, not just whichever file happened to render last.
    ratios = {}
    for label, (vals, xe, ye) in others:
        if vals.shape != ref_vals.shape or not (np.allclose(xe, ref_xe) and np.allclose(ye, ref_ye)):
            log.warning("'%s' skipped for %s: binning mismatch vs reference", label, hist_name)
            continue
        with np.errstate(divide="ignore", invalid="ignore"):
            ratios[label] = np.where(ref_vals != 0, vals / ref_vals, np.nan)

    if not ratios:
        print(f"No comparable (matching binning) run plans for {hist_name}.")
        return

    norm = _compute_shared_norm(list(ratios.values()))

    saved = []
    for label, ratio in ratios.items():
        fig, ax = plt.subplots(figsize=(8, 6.5))
        mesh = plot_2d_ratio(ax, ratio, ref_xe, ref_ye,
                              f"{hist_label}: {label} / {reference_label}", norm)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        fig.colorbar(mesh, ax=ax, label=f"{hist_label} ratio")
        out_path = f"{save_prefix}_{_sanitize_filename(label)}.pdf"
        fig.savefig(out_path, bbox_inches="tight")
        plt.close(fig)
        saved.append(out_path)

    print(f"Saved {len(saved)} file(s) for {hist_name}:")
    for p in saved:
        print(f"    {p}")


# ============================================================================
# 1D "per active parameter" ratio lines (recovered from the sparse 2D grid)
# ============================================================================

def get_active_bin_order(ref_vals, atol=0.0):
    """
    Returns the flat (row-major) indices of nonzero bins in the reference
    array, in a fixed order -- this defines "parameter 0, 1, 2, ..." for
    every file, so ratios line up bin-for-bin even though the underlying
    histogram is a sparse 2D grid, not a flat 1D list.
    """
    flat = ref_vals.flatten(order="C")
    active_idx = np.flatnonzero(np.abs(flat) > atol)
    return active_idx


def format_bin_label(ix, iy, xedges, yedges, xname=r"$E_\nu$", yname="Bias"):
    """Human-readable kinematic region for one 2D bin, e.g.
    'Eν [0.5, 1)\\nBias [-0.1, 0)'."""
    return (f"{xname} [{xedges[ix]:.2g}, {xedges[ix+1]:.2g})\n"
            f"{yname} [{yedges[iy]:.2g}, {yedges[iy+1]:.2g})")


def plot_param_ratios_1d(results, hist_name_mean, hist_name_std, reference_label,
                          save_path="template_param_ratios.pdf",
                          params_per_page=PARAMS_PER_PAGE):
    """
    Multi-page PDF: the 1D per-active-parameter ratio lines, split into
    pages of `params_per_page` parameters each (so a file with hundreds
    of active parameters stays readable instead of one very long plot).
    """
    if reference_label not in results:
        print(f"ERROR: reference '{reference_label}' failed to load.")
        return
    ref = results[reference_label]
    if hist_name_mean not in ref or hist_name_std not in ref:
        print(f"ERROR: reference is missing {hist_name_mean} or {hist_name_std}.")
        return

    ref_mean_vals, xe, ye = ref[hist_name_mean]
    ref_std_vals, _, _ = ref[hist_name_std]
    active_idx = get_active_bin_order(ref_mean_vals)
    if active_idx.size == 0:
        print("Reference has no nonzero parameter bins -- nothing to plot.")
        return
    ref_mean_1d = ref_mean_vals.flatten(order="C")[active_idx]
    ref_std_1d = ref_std_vals.flatten(order="C")[active_idx]
    nparams = len(active_idx)

    bin_labels = []
    for flat_i in active_idx:
        ix, iy = np.unravel_index(flat_i, ref_mean_vals.shape, order="C")
        bin_labels.append(format_bin_label(ix, iy, xe, ye))

    # Precompute each non-reference run plan's full ratio arrays once
    per_plan_ratios = {}
    color_of = {}
    color_i = 0
    for label, r in results.items():
        if label == reference_label:
            continue
        if hist_name_mean not in r or hist_name_std not in r:
            continue
        mean_vals, mxe, mye = r[hist_name_mean]
        std_vals, sxe, sye = r[hist_name_std]
        if mean_vals.shape != ref_mean_vals.shape or not (np.allclose(mxe, xe) and np.allclose(mye, ye)):
            log.warning("Skipping '%s' in 1D plot: binning mismatch vs reference", label)
            continue
        mean_1d = mean_vals.flatten(order="C")[active_idx]
        std_1d = std_vals.flatten(order="C")[active_idx]
        with np.errstate(divide="ignore", invalid="ignore"):
            mean_ratio = np.where(ref_mean_1d != 0, mean_1d / ref_mean_1d, np.nan)
            std_ratio = np.where(ref_std_1d != 0, std_1d / ref_std_1d, np.nan)
        per_plan_ratios[label] = (mean_ratio, std_ratio)
        color_of[label] = PALETTE[color_i % len(PALETTE)]
        color_i += 1

    if not per_plan_ratios:
        print("No comparable run plans for the 1D parameter plot.")
        return

    npages = int(np.ceil(nparams / params_per_page))
    with PdfPages(save_path) as pdf:
        for page in range(npages):
            lo, hi = page * params_per_page, min((page + 1) * params_per_page, nparams)
            x = np.arange(hi - lo)
            page_labels = bin_labels[lo:hi]

            fig, (ax_mean, ax_std) = plt.subplots(
                2, 1, figsize=(max(8, (hi - lo) * 0.5), 7.5), sharex=True,
                gridspec_kw={"hspace": 0.12},
            )
            for label, (mean_ratio, std_ratio) in per_plan_ratios.items():
                ax_mean.plot(x, mean_ratio[lo:hi], "o-", ms=4.5, lw=1.3,
                             color=color_of[label], label=label, alpha=0.9)
                ax_std.plot(x, std_ratio[lo:hi], "o-", ms=4.5, lw=1.3,
                            color=color_of[label], alpha=0.9)

            for ax in (ax_mean, ax_std):
                ax.axhline(1.0, color="0.3", lw=1.0, ls="--", zorder=0)

            ax_mean.set_ylabel(f"Mean / {reference_label}")
            ax_std.set_ylabel(f"Std Dev / {reference_label}")
            ax_std.set_xlabel("Template Parameter")
            title = f"Template Parameter Ratios (relative to {reference_label})"
            if npages > 1:
                title += f"  \u2014  page {page + 1}/{npages} (params {lo}\u2013{hi - 1})"
            ax_mean.set_title(title, pad=10, fontsize=12)
            ax_mean.legend(loc="upper right", ncol=2, fontsize=8, framealpha=0.9)
            ax_std.set_xticks(x)
            ax_std.set_xticklabels(page_labels, rotation=90, fontsize=6.5, ha="center")
            ax_mean.set_xlim(-0.5, len(x) - 0.5)

            fig.align_ylabels([ax_mean, ax_std])
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

    print(f"Saved -> {save_path}  ({npages} page(s), {nparams} parameters)")


# ============================================================================
# Main
# ============================================================================

def main():
    logging.basicConfig(level=logging.INFO, format="%(levelname)s  %(message)s")
    plt.rcParams.update(THESIS_RC)

    results = {}
    for label, path in run_plan_files.items():
        print(f"[{label}]  {path}")
        results[label] = load_run_plan_hists(path)

    if REFERENCE_LABEL not in results or not results[REFERENCE_LABEL]:
        print(f"\nERROR: reference '{REFERENCE_LABEL}' failed to load anything.")
        return

    # 2D ratio heatmaps for every histogram
    for hist_name, hist_label in HIST_2D_SPECS:
        plot_all_2d_ratios(results, hist_name, hist_label, REFERENCE_LABEL)

    # 1D "per active parameter" ratio lines for mean/stddev
    plot_param_ratios_1d(results, "xsec_param_mean", "xsec_param_stddev", REFERENCE_LABEL)


if __name__ == "__main__":
    main()
