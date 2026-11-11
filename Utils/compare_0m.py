#!/usr/bin/env python3
"""
compare_onaxis_standalone.py
=============================

Self-contained script: compares the on-axis (0m) posterior across several
run plans, for TrueNeutrinoEnergy and RecoNeutrinoEnergy, each saved as
its own PDF. No other local files needed -- just uproot/numpy/matplotlib.

Handles ROOT files with inconsistent internal layout (toys/Asimov nested
at different depths, sample names embedded in flat histogram names or in
subdirectory names) by recursively scanning for matching histograms
rather than assuming one fixed structure.

Usage
-----
1. Edit `run_plan_files` below with your current file paths.
2. python compare_onaxis_standalone.py
"""

import re
import logging
from dataclasses import dataclass
from typing import Optional, Sequence, List

import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

log = logging.getLogger(__name__)

# ============================================================================
# EDIT THIS: your run plans, one ROOT path each
# ============================================================================

run_plan_files = {
    "On Axis only":
        "/home/abipeake/scratch/July2026/OnAxis_templatecovmat_new/Posteriorpredictive_out_posteriorpredictive_withtoys_ND.root",
    "50% On Axis + 50% Off Axis":
        "/home/abipeake/scratch/July2026/Asimov_allsystematics_onplusoffaxis_postadaptive/Posteriorpredictive_out_posteriorpredictive_withtoys.root",
    "75% On Axis + 25% Off Axis":
        "/home/abipeake/scratch/July2026/Runplan_75onaxis_25offaxis_post/Posteriorpredictive_out_posteriorpredictive_withtoys.root",
    "25% On Axis + 75% Off Axis":
        "/home/abipeake/scratch/July2026/Runplan_25onaxis_75offaxis_post/Posteriorpredictive_out_posteriorpredictive_withtoys.root",
    "Equal time at each position":
        "/home/abipeake/scratch/July2026/Runplan_equal_allposition_post/Posteriorpredictive_out_posteriorpredictive_withtoys.root",
    "50% On Axis + 50% Off Axis (12m)":
        "/home/abipeake/scratch/July2026/Runplan_50onaxis_50upto12mpost/Posteriorpredictive_out_posteriorpredictive_withtoys_NDnews.root",
}

VARIABLES = ("TrueNeutrinoEnergy", "RecoNeutrinoEnergy")

VAR_META = {
    "TrueNeutrinoEnergy": {"xlabel": "True Neutrino Energy (GeV)"},
    "RecoNeutrinoEnergy": {"xlabel": "Reconstructed Neutrino Energy (GeV)"},
    "Enubias":            {"xlabel": "E_{\u03bd} bias (GeV)"},
}

# Okabe-Ito colourblind-safe palette, verified via CVD simulation to keep
# every pair of colours distinguishable even under full colour-vision
# deficiency (deuteranopia/protanopia/tritanopia).
PALETTE: List[str] = [
    "#0072B2",  # blue
    "#E69F00",  # orange
    "#009E73",  # bluish green
    "#D55E00",  # vermillion
    "#CC79A7",  # reddish purple
    "#56B4E9",  # sky blue
    "#F0E442",  # yellow
]

THESIS_RC = {
    "font.family":         "serif",
    "font.serif":          ["STIXTwo Text", "DejaVu Serif", "Times New Roman"],
    "mathtext.fontset":    "stix",
    "font.size":           11,
    "axes.titlesize":      15,
    "axes.labelsize":      13,
    "xtick.labelsize":     10,
    "ytick.labelsize":     10,
    "legend.fontsize":     10,
    "pdf.fonttype":        42,
    "ps.fonttype":         42,
    "figure.dpi":          150,
    "lines.linewidth":     1.8,
    "xtick.direction":     "in",
    "ytick.direction":     "in",
    "xtick.top":           True,
    "ytick.right":         True,
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "axes.grid":           True,
    "grid.alpha":          0.3,
    "grid.linewidth":      0.6,
    "axes.linewidth":      0.8,
}

# Patterns that count as "on axis" (0m), checked in order, case-insensitive.
ON_AXIS_PATTERNS = [
    re.compile(r"^onaxis", re.IGNORECASE),
    re.compile(r"offaxis0(\.0+)?m", re.IGNORECASE),
]


# ============================================================================
# ROOT parsing (recursive, layout-agnostic)
# ============================================================================

def safe_get(f, path):
    try:
        return f[path]
    except Exception:
        return None


def is_2d(h):
    return len(h.axes) == 2


def get_1d(h):
    return h.values(), h.axes[0].edges()


def strip_cycle(key):
    return key.split(";")[0]


def is_dir(obj):
    return hasattr(obj, "keys") and not hasattr(obj, "axes")


def collect_all_leaves(node, max_depth=10, _depth=0):
    """Recursively collect every histogram under `node`, however deep it's
    nested, as {leaf_name: histogram_object}. Handles files that mix
    nesting depths inconsistently (e.g. the first toy of a variable
    sitting directly under {det}/, the rest one level deeper)."""
    leaves = {}
    if _depth > max_depth:
        return leaves
    try:
        keys = node.keys(recursive=False)
    except Exception:
        return leaves
    for k in keys:
        ks = strip_cycle(k)
        try:
            obj = node[k]
        except Exception:
            continue
        if is_dir(obj):
            leaves.update(collect_all_leaves(obj, max_depth, _depth + 1))
        else:
            leaves[ks] = obj
    return leaves


def _toy_regex(det, var):
    return re.compile(
        rf"^{re.escape(det)}_(?P<sample>.+)_{re.escape(var)}_posterior_toy"
        rf"(?:_(?P<variant>[A-Za-z]+))?_(?P<idx>\d+)$"
    )


def _asimov_regex(det, var, suffix):
    return re.compile(
        rf"^{re.escape(det)}_(?P<sample>.+)_{re.escape(var)}_{re.escape(suffix)}$"
    )


def find_detectors(f):
    top = [strip_cycle(k) for k in f.keys(recursive=False)]
    return [k for k in top if k not in ("Asimov", "AsimovFakeData")
            and safe_get(f, k) is not None]


def discover_subsamples(f, det, var):
    leaves = collect_all_leaves(safe_get(f, det))
    toy_re = _toy_regex(det, var)
    samples = set()
    for name in leaves:
        m = toy_re.match(name)
        if m:
            samples.add(m.group("sample"))
    return sorted(samples)


def find_on_axis_sample(samples):
    """Return the on-axis sample from a list, or None if not found."""
    if len(samples) == 1:
        return samples[0]  # flat single-sample file
    for pat in ON_AXIS_PATTERNS:
        matches = [s for s in samples if pat.search(s)]
        if matches:
            if len(matches) > 1:
                log.warning("Multiple on-axis-like matches %s -- using first", matches)
            return matches[0]
    return None


def find_asimov_hist(f, top, det, var, sample):
    leaves = collect_all_leaves(safe_get(f, top))
    rx = _asimov_regex(det, var, top)
    matches = [(name, obj) for name, obj in leaves.items()
               if (m := rx.match(name)) and m.group("sample") == sample]
    if not matches:
        return None
    return matches[0][1]


def collect_toys(f, det, var, sample, edges, preferred_variant="combined"):
    leaves = collect_all_leaves(safe_get(f, det))
    toy_re = _toy_regex(det, var)
    tagged, untagged = [], []
    for name, obj in leaves.items():
        m = toy_re.match(name)
        if not m or m.group("sample") != sample:
            continue
        if m.group("variant") == preferred_variant:
            tagged.append((obj, name))
        elif m.group("variant") is None:
            untagged.append((obj, name))
    use = tagged if tagged else untagged

    toys = []
    for h, key in use:
        try:
            if is_2d(h):
                continue
            v, e = get_1d(h)
            if not np.allclose(e, edges):
                continue
            toys.append(v / np.diff(e))
        except Exception:
            pass
    return toys


def resolve_on_axis(path, var):
    """Open `path`, find (det, sample) for the on-axis position. Prints
    what it finds along the way. Returns (det, sample) or (None, None)."""
    try:
        f = uproot.open(path)
    except Exception as exc:
        print(f"    [could not open file: {exc}]")
        return None, None, None

    dets = find_detectors(f)
    print(f"    detectors found: {dets}")
    search_order = (["ND"] if "ND" in dets else []) + [d for d in dets if d != "ND"]

    for det in search_order:
        samples = discover_subsamples(f, det, var)
        print(f"    {det}/{var} samples: {samples}")
        if not samples:
            continue
        on_axis = find_on_axis_sample(samples)
        if on_axis:
            print(f"    -> using det={det!r}, sample={on_axis!r}")
            return f, det, on_axis

    print(f"    -> could not resolve an on-axis sample for {var}")
    return f, None, None


def extract_fit_result(f, det, var, sample, percentiles=(5, 95)):
    h_asimov = find_asimov_hist(f, "Asimov", det, var, sample)
    if h_asimov is None:
        raise ValueError(f"No Asimov histogram for {det}/{sample}/{var}")
    asimov_vals, edges = get_1d(h_asimov)
    cv = asimov_vals / np.diff(edges)

    toys = collect_toys(f, det, var, sample, edges)
    if not toys:
        raise ValueError(f"No toys found for {det}/{sample}/{var}")
    toys_arr = np.vstack(toys)

    mean = np.mean(toys_arr, axis=0)
    std = np.std(toys_arr, axis=0)
    p_lo = np.percentile(toys_arr, percentiles[0], axis=0)
    p_hi = np.percentile(toys_arr, percentiles[1], axis=0)
    return dict(edges=edges, cv=cv, fit=mean, fit_sigma=std, fit_p05=p_lo, fit_p95=p_hi)


# ============================================================================
# Comparison plot
# ============================================================================

@dataclass
class FitResult:
    label: str
    bin_edges: np.ndarray
    cv: np.ndarray
    fit: np.ndarray
    fit_sigma: Optional[np.ndarray] = None
    fit_p05: Optional[np.ndarray] = None
    fit_p95: Optional[np.ndarray] = None
    color: Optional[str] = None

    @property
    def bin_centers(self):
        return 0.5 * (self.bin_edges[:-1] + self.bin_edges[1:])


def _step_xy(edges, values):
    x = np.repeat(edges, 2)[1:-1]
    y = np.repeat(values, 2)
    return x, y


def _band_step(edges, lo, hi):
    x = np.repeat(edges, 2)[1:-1]
    return x, np.repeat(lo, 2), np.repeat(hi, 2)


def plot_comparison(fits: Sequence[FitResult], xlabel, ylabel, title,
                     reference=None, show_cv=True, save_path=None):
    for i, fit in enumerate(fits):
        if fit.color is None:
            fit.color = PALETTE[i % len(PALETTE)]
    if reference is None:
        reference = fits[0]

    fig, (ax_main, ax_prec, ax_ratio) = plt.subplots(
        3, 1, figsize=(8, 9), sharex=True,
        gridspec_kw={"height_ratios": [3, 1.1, 1.1], "hspace": 0.08},
    )

    for fit in fits:
        x, y_fit = _step_xy(fit.bin_edges, fit.fit)
        ax_main.plot(x, y_fit, color=fit.color, lw=1.8, label=fit.label, zorder=3)
        if show_cv:
            x_cv, y_cv = _step_xy(fit.bin_edges, fit.cv)
            ax_main.plot(x_cv, y_cv, color=fit.color, lw=1.1, ls=(0, (4, 2)), alpha=0.75, zorder=2)
        if fit.fit_p05 is not None and fit.fit_p95 is not None:
            xb, lo, hi = _band_step(fit.bin_edges, fit.fit_p05, fit.fit_p95)
            ax_main.fill_between(xb, lo, hi, color=fit.color, alpha=0.12, linewidth=0, zorder=1)
        if fit.fit_sigma is not None:
            xb, lo, hi = _band_step(fit.bin_edges, fit.fit - fit.fit_sigma, fit.fit + fit.fit_sigma)
            ax_main.fill_between(xb, lo, hi, color=fit.color, alpha=0.25, linewidth=0, zorder=1.5)

    ax_main.set_ylabel(ylabel)
    ax_main.set_title(title, pad=10)
    ax_main.set_ylim(bottom=0)
    ax_main.legend(loc="upper right")

    for fit in fits:
        if fit.fit_sigma is None:
            continue
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(fit.fit != 0, fit.fit_sigma / fit.fit, np.nan)
        x, y = _step_xy(fit.bin_edges, ratio)
        ax_prec.plot(x, y, color=fit.color, lw=1.6)
    ax_prec.axhline(0, color="0.6", lw=0.8, ls=":")
    ax_prec.set_ylabel(r"$\sigma_{\mathrm{Fit}}/\mu$")

    ref_y = reference.fit
    for fit in fits:
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(ref_y != 0, fit.fit / ref_y, np.nan)
        x, y = _step_xy(fit.bin_edges, ratio)
        ax_ratio.plot(x, y, color=fit.color, lw=1.6)
        if fit.fit_sigma is not None:
            with np.errstate(divide="ignore", invalid="ignore"):
                lo = np.where(ref_y != 0, (fit.fit - fit.fit_sigma) / ref_y, np.nan)
                hi = np.where(ref_y != 0, (fit.fit + fit.fit_sigma) / ref_y, np.nan)
            xb, lo_s, hi_s = _band_step(fit.bin_edges, lo, hi)
            ax_ratio.fill_between(xb, lo_s, hi_s, color=fit.color, alpha=0.18, linewidth=0)
    ax_ratio.axhline(1.0, color="0.3", lw=0.9, ls="--")
    ax_ratio.set_ylabel(f"Fit / {reference.label}")
    ax_ratio.set_xlabel(xlabel)

    all_edges = np.concatenate([fit.bin_edges for fit in fits])
    ax_main.set_xlim(all_edges.min(), all_edges.max())
    for ax in (ax_main, ax_prec, ax_ratio):
        ax.tick_params(which="both", direction="in", top=True, right=True)
    fig.align_ylabels([ax_main, ax_prec, ax_ratio])

    if save_path:
        fig.savefig(save_path, bbox_inches="tight")
    plt.close(fig)


# ============================================================================
# Main
# ============================================================================

def main():
    logging.basicConfig(level=logging.INFO, format="%(levelname)s  %(message)s")
    plt.rcParams.update(THESIS_RC)

    for var in VARIABLES:
        print(f"\n=== {var} ===")
        fits = []
        for label, path in run_plan_files.items():
            print(f"  [{label}]  {path}")
            f, det, sample = resolve_on_axis(path, var)
            if det is None:
                print(f"    SKIPPING '{label}'")
                continue
            try:
                r = extract_fit_result(f, det, var, sample)
                fits.append(FitResult(
                    label=label, bin_edges=r["edges"], cv=r["cv"], fit=r["fit"],
                    fit_sigma=r["fit_sigma"], fit_p05=r["fit_p05"], fit_p95=r["fit_p95"],
                ))
            except Exception as exc:
                print(f"    SKIPPING '{label}' -- extraction failed: {exc}")

        if not fits:
            print(f"  No fits loaded for {var} -- nothing to plot.")
            continue

        save_path = f"onaxis_comparison_{var}.pdf"
        xlabel = VAR_META.get(var, {}).get("xlabel", var)
        plot_comparison(
            fits, xlabel=xlabel, ylabel="Events / GeV",
            title=f"On-Axis {var} \u2014 Run Plan Comparison",
            reference=fits[0], save_path=save_path,
        )
        print(f"  Saved -> {save_path}  ({len(fits)}/{len(run_plan_files)} run plans)")


if __name__ == "__main__":
    main()