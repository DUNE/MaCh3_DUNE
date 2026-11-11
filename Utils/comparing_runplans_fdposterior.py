#!/usr/bin/env python3
"""
plot_posterior_all.py
======================

One-stop script for posterior predictive plots from the neutrino fit's
ROOT output. Combines three things that used to be separate files:

  1. ROOT parsing layer   - ND (18-position, per-sample subdirectory) and
                             FD (flat, single-sample) layouts, auto-detected.
  2. Single-file plotter  - the standard 3-panel PDF (per detector/sample/
                             variable), same as before.
  3. Multi-fit comparison - overlay several run plans' posteriors on one
                             figure, for thesis-quality comparison plots.

Directory layouts handled (auto-detected per detector/variable):

  (A) Nested, per-sample subdirectory (e.g. ND with 18 OffAxis positions):

        Asimov/{det}/{var}/{sample}/{det}_{sample}_{var}_Asimov
        AsimovFakeData/{det}/{var}/{sample}/..._AsimovFakeData
        {det}/{var}/{sample}/..._posterior_toy[_{variant}]_{NNN}

  (B) Flat, single-sample-per-detector (e.g. FD):

        Asimov/{det}/{var}/{det}_{sample}_{var}_Asimov
        AsimovFakeData/{det}/{var}/{det}_{sample}_{var}_AsimovFakeData
        {det}/{var}/{det}_{sample}_{var}_posterior_toy[_{variant}]_{NNN}

  {variant}, when present, is one of "xsecOnly", "fluxOnly", "combined".

Usage
-----
Generate the standard multi-page PDF for one file (unchanged behaviour):

    python plot_posterior_all.py file.root [data.root] [output.pdf]

Build a multi-run-plan comparison figure (edit `run_plan_files` below, or
import the functions into your own script/notebook):

    python plot_posterior_all.py --compare
"""

from __future__ import annotations

import os
import re
import sys
import logging
from dataclasses import dataclass
from typing import List, Optional, Sequence

import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.backends.backend_pdf import PdfPages

log = logging.getLogger(__name__)


# ============================================================================
# 1. STYLE
# ============================================================================

# Preferred toy variant when a file tags toys by variant (layout A).
# Files with untagged toys (layout B) use all toys found, regardless of
# this setting.
TOY_VARIANT = "combined"

# Okabe-Ito colourblind-safe palette (used for the single-file PDF plots)
_C = {
    "blue":   "#0072B2",
    "orange": "#E69F00",
    "red":    "#D55E00",
    "green":  "#009E73",
    "sky":    "#56B4E9",
    "purple": "#CC79A7",
    "yellow": "#F0E442",
    "black":  "#000000",
}

# A wider, muted palette for overlaying many run plans at once (comparison
# plot). Extendable; cycles if you have more fits than colours.
# PALETTE: List[str] = [
#     "#2E5C8A",  # deep steel blue
#     "#D98E04",  # warm amber
#     "#8C1C13",  # muted brick red
#     "#3F7D5C",  # forest green
#     "#6A4C93",  # muted violet
#     "#B0763A",  # ochre
#     "#4C6EF5",  # indigo
#     "#A63A50",  # dusty rose
# ]

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
    "font.family":          "serif",
    "font.serif":           ["STIXTwo Text", "DejaVu Serif", "Times New Roman"],
    "mathtext.fontset":     "stix",
    "font.size":            11,
    "axes.titlesize":       11,
    "axes.labelsize":       11,
    "xtick.labelsize":      10,
    "ytick.labelsize":      10,
    "legend.fontsize":      9,
    "pdf.fonttype":         42,
    "ps.fonttype":          42,
    "figure.dpi":           150,
    "lines.linewidth":      1.8,
    "patch.linewidth":      1.0,
    "xtick.direction":      "in",
    "ytick.direction":      "in",
    "xtick.top":            True,
    "ytick.right":          True,
    "xtick.minor.visible":  True,
    "ytick.minor.visible":  True,
    "axes.grid":            True,
    "grid.alpha":           0.3,
    "grid.linewidth":       0.6,
    "axes.linewidth":       0.8,
}


def set_thesis_style() -> None:
    plt.rcParams.update(THESIS_RC)


VAR_META = {
    "TrueNeutrinoEnergy":  {"xlabel": "True Neutrino Energy (GeV)"},
    "RecoNeutrinoEnergy":  {"xlabel": "Reconstructed Neutrino Energy (GeV)"},
    "Enubias":             {"xlabel": "E_{\u03bd} bias (GeV)"},
}


# ============================================================================
# 2. ROOT PARSING LAYER (layout-agnostic: handles ND-nested and FD-flat)
# ============================================================================

def safe_get(f, path):
    """Return f[path] or None, logging any exception."""
    try:
        return f[path]
    except Exception as exc:
        log.debug("safe_get(%s): %s", path, exc)
        return None


def is_2d(h):
    return len(h.axes) == 2


def get_1d(h):
    return h.values(), h.axes[0].edges()


def get_2d(h):
    return h.values(), h.axes[0].edges(), h.axes[1].edges()


def extend(v):
    """Repeat last element so step plots close correctly."""
    return np.append(v, v[-1])


def mask_empty(v):
    v = v.astype(float)
    v[v == 0] = np.nan
    return v


def strip_cycle(key):
    """Strip ROOT cycle suffix (';N') from a key string."""
    return key.split(";")[0]


def is_dir(obj):
    """True if a keys()/values() ROOT object is a directory, not a histogram."""
    return hasattr(obj, "keys") and not hasattr(obj, "axes")


def collect_all_leaves(node, max_depth=10, _depth=0):
    """
    Recursively collect every histogram under `node`, regardless of how
    deep it's nested, as {leaf_name: histogram_object}.

    This is deliberately layout-agnostic: real files have been observed to
    mix depths inconsistently (e.g. the first toy of a variable sitting
    directly under {det}/, with the rest nested one level deeper under
    {det}/{var}/; or Asimov histograms sitting directly under Asimov/
    rather than under Asimov/{det}/{var}/ as the directory names imply).
    Recursing and matching on the leaf name (which already encodes
    det/sample/var/toy-index) sidesteps needing to know the exact depth.
    """
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


# Cache of (id(uproot_file), path) -> {leaf_name: hist_obj}, so repeated
# discover_subsamples/find_asimov_hist/collect_toys calls against the same
# file don't re-walk a large subtree (e.g. an 18-position ND file with
# thousands of toys) over and over.
_LEAF_CACHE = {}


def _leaves_under(f, path):
    key = (id(f), path)
    if key not in _LEAF_CACHE:
        node = safe_get(f, path)
        _LEAF_CACHE[key] = collect_all_leaves(node) if node is not None else {}
    return _LEAF_CACHE[key]


def _toy_regex(det, var):
    # Matches: {det}_{sample}_{var}_posterior_toy[_{variant}]_{NNN}
    return re.compile(
        rf"^{re.escape(det)}_(?P<sample>.+)_{re.escape(var)}_posterior_toy"
        rf"(?:_(?P<variant>[A-Za-z]+))?_(?P<idx>\d+)$"
    )


def _asimov_regex(det, var, suffix):
    # suffix is "Asimov" or "AsimovFakeData"
    return re.compile(
        rf"^{re.escape(det)}_(?P<sample>.+)_{re.escape(var)}_{re.escape(suffix)}$"
    )


def discover_subsamples(f, det, variables):
    """
    Collect the set of subsample names available for a detector by
    recursively scanning every leaf under {det}/ and parsing the sample
    out of any histogram name matching the known toy pattern, regardless
    of nesting depth.
    """
    leaves = _leaves_under(f, det)
    samples = set()
    for var in variables:
        toy_re = _toy_regex(det, var)
        for name in leaves:
            m = toy_re.match(name)
            if m:
                samples.add(m.group("sample"))
    return sorted(samples)


def find_asimov_hist(f, top, det, var, sample):
    """
    Find the Asimov / AsimovFakeData histogram for (det, var, sample) by
    recursively scanning every leaf under `top` ("Asimov" or
    "AsimovFakeData"), regardless of nesting depth.
    """
    leaves = _leaves_under(f, top)
    rx = _asimov_regex(det, var, top)
    matches = [(name, obj) for name, obj in leaves.items()
               if (m := rx.match(name)) and m.group("sample") == sample]
    if not matches:
        return None
    if len(matches) > 1:
        log.warning("Multiple %s matches for %s/%s/%s -- using first",
                    top, det, sample, var)
    return matches[0][1]


def collect_toys(f, det, var, sample, edges, preferred_variant=TOY_VARIANT):
    """
    Collect all posterior-toy histograms for (det, var, sample) by
    recursively scanning every leaf under {det}/, regardless of nesting
    depth. Returns a list of per-bin value arrays (divided by bin width).
    """
    leaves = _leaves_under(f, det)
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
                log.debug("Edge mismatch for toy %s -- skipped.", key)
                continue
            toys.append(v / np.diff(e))
        except Exception as exc:
            log.debug("Could not read toy %s: %s", key, exc)
    return toys


def load_data_histograms(data_root_path):
    """Open a data ROOT file, return { name: (values/bin_width, edges) }."""
    if data_root_path is None or not os.path.exists(data_root_path):
        return {}
    try:
        df = uproot.open(data_root_path)
    except Exception as exc:
        log.error("Could not open data file %s: %s", data_root_path, exc)
        return {}

    sample_folder = safe_get(df, "SampleFolder")
    if sample_folder is not None:
        data_hists = {}
        for k in sample_folder.keys(recursive=False):
            ks = strip_cycle(k)
            if not ks.startswith("data_"):
                continue
            try:
                h = sample_folder[k]
                if is_2d(h):
                    continue
                vals, edges = get_1d(h)
                data_hists[ks] = (vals / np.diff(edges), edges)
            except Exception as exc:
                log.warning("Could not read data histogram %s: %s", ks, exc)
        return data_hists

    data_hists = {}
    for k in df.keys(recursive=False):
        ks = strip_cycle(k)
        try:
            h = df[k]
            if is_2d(h):
                continue
            vals, edges = get_1d(h)
            data_hists[ks] = (vals / np.diff(edges), edges)
        except Exception as exc:
            log.warning("Could not read data histogram %s: %s", ks, exc)
    return data_hists


def _friendly_label(name):
    for prefix in ("data_", "hOffAxis", "h"):
        if name.startswith(prefix):
            name = name[len(prefix):]
            break
    label = name.replace("_", " ").strip()
    return f"Fake Data: {label}" if label else "Fake Data"


# ----------------------------------------------------------------------
# High-level extraction: mean/std/percentiles of the toys, in one call
# ----------------------------------------------------------------------

def extract_fit_result(root_path, det, var, sample=None, percentiles=(5, 95)):
    """
    Open `root_path` and return a dict with everything needed to build a
    FitResult for the comparison plot:
        edges, cv, fit (toy mean), fit_sigma (toy std), fit_p05, fit_p95

    If `sample` is None and there's exactly one sample for this det/var,
    it's used automatically; otherwise pass `sample=...` explicitly
    (e.g. "OffAxis0m_numuCC_numode" for the on-axis ND position).
    """
    f = uproot.open(root_path)

    if sample is None:
        samples = discover_subsamples(f, det, [var])
        if len(samples) != 1:
            raise ValueError(
                f"Expected exactly one sample for {det}/{var} in {root_path}, "
                f"found {samples}. Pass `sample=...` explicitly."
            )
        sample = samples[0]

    h_asimov = find_asimov_hist(f, "Asimov", det, var, sample)
    if h_asimov is None:
        raise ValueError(f"No Asimov histogram for {det}/{sample}/{var} in {root_path}")
    asimov_vals, edges = get_1d(h_asimov)
    bin_widths = np.diff(edges)
    cv = asimov_vals / bin_widths

    toys = collect_toys(f, det, var, sample, edges)
    if not toys:
        raise ValueError(f"No toys found for {det}/{sample}/{var} in {root_path}")
    toys_arr = np.vstack(toys)

    mean = np.mean(toys_arr, axis=0)
    std = np.std(toys_arr, axis=0)
    p_lo = np.percentile(toys_arr, percentiles[0], axis=0)
    p_hi = np.percentile(toys_arr, percentiles[1], axis=0)

    return dict(edges=edges, cv=cv, fit=mean, fit_sigma=std,
                fit_p05=p_lo, fit_p95=p_hi, sample=sample)


# ============================================================================
# 3. SINGLE-FILE 3-PANEL PDF PLOTTER (unchanged behaviour from before)
# ============================================================================

def _toy_envelope(ax, toys, edges, percentiles=(5, 95)):
    lo = np.percentile(toys, percentiles[0], axis=0)
    hi = np.percentile(toys, percentiles[1], axis=0)
    ax.fill_between(edges, extend(lo), extend(hi),
                     step="post", alpha=0.12, color=_C["black"],
                     label=f"{percentiles[0]}\u2013{percentiles[1]}th pct")


def plot_1d(pdf, asimov, mean, std, edges, title, xlabel,
            toys=None, data_hists=None, is_reco=False,
            fakedata_asimov=None,
            no_mc_prediction=False,
            no_fakedata_asimov_line=False,
            sigma_blue=False,
            sigma_toy_envelope=False,
            xlim_min=None):

    fig, (ax1, ax2, ax3) = plt.subplots(
        3, 1, figsize=(6.5, 7.5),
        gridspec_kw={"height_ratios": [3, 1, 1]},
        sharex=True,
    )
    fig.subplots_adjust(hspace=0.05)

    e = edges
    m = extend(mean)
    u = extend(mean + std)
    l = extend(mean - std)
    a = extend(asimov)

    if toys is not None and len(toys):
        _toy_envelope(ax1, toys, e)

    ax1.fill_between(e, l, u, step="post", alpha=0.25,
                      color=_C["blue"], label=r"$\pm1\sigma$ (posterior)")
    ax1.step(e, m, where="post", label="Fit Result",
              linewidth=1.8, color=_C["blue"])

    if not no_mc_prediction:
        ax1.step(e, a, where="post", label="CV Prediction",
                  linewidth=1.8, color=_C["orange"])

    if fakedata_asimov is not None and not no_fakedata_asimov_line:
        ax1.step(e, extend(fakedata_asimov), where="post",
                  label="Fake Data Asimov", linewidth=1.8,
                  color=_C["red"], linestyle=(0, (5, 2)))

    if is_reco and data_hists:
        data_colors = [_C["green"], _C["purple"], _C["sky"], _C["yellow"]]
        markers = ["o", "s", "^", "D"]
        for idx, (name, (dvals, dedges)) in enumerate(data_hists.items()):
            centres = 0.5 * (dedges[:-1] + dedges[1:])
            ax1.errorbar(
                centres, dvals,
                yerr=np.sqrt(np.maximum(dvals * np.diff(dedges), 0)) / np.diff(dedges),
                fmt=markers[idx % len(markers)],
                markersize=4, linewidth=1.2, capsize=3,
                label=_friendly_label(name),
                color=data_colors[idx % len(data_colors)],
                zorder=5,
            )

    ax1.set_ylabel("Events / GeV")
    ax1.set_title(title)
    ax1.legend(loc="best", framealpha=0.9, edgecolor="0.7")
    plt.setp(ax1.get_xticklabels(), visible=False)

    rel_unc = np.divide(std, mean, out=np.zeros_like(std), where=mean != 0)
    sigma_color = _C["blue"] if sigma_blue else _C["red"]

    if sigma_toy_envelope and toys is not None and len(toys):
        toy_rel_each = np.divide(
            np.std(toys, axis=0, keepdims=True) * np.ones_like(toys),
            np.maximum(toys, 1e-30),
        )
        lo = np.percentile(toy_rel_each, 5, axis=0)
        hi = np.percentile(toy_rel_each, 95, axis=0)
        ax2.fill_between(e, extend(lo), extend(hi),
                          step="post", alpha=0.15, color="gray",
                          label="5\u201395th pct")

    ax2.step(e, extend(rel_unc), where="post", color=sigma_color,
              label=r"$\sigma_\mathrm{Fit}\,/\,\mu$", linewidth=1.8)
    ax2.axhline(0, color="black", linewidth=0.6, linestyle=":", alpha=0.5)
    ax2.set_ylabel(r"$\sigma_\mathrm{Fit}\,/\,\mu$")
    ax2.legend(loc="upper right", framealpha=0.9, edgecolor="0.7")
    plt.setp(ax2.get_xticklabels(), visible=False)

    bias    = np.divide(asimov, mean,       out=np.ones_like(asimov), where=mean != 0)
    bias_up = np.divide(asimov, mean + std, out=np.ones_like(asimov), where=(mean + std) != 0)
    bias_dn = np.divide(asimov, mean - std, out=np.ones_like(asimov), where=(mean - std) != 0)

    if fakedata_asimov is not None and not no_fakedata_asimov_line:
        fd_bias = np.divide(fakedata_asimov, mean,
                             out=np.ones_like(fakedata_asimov), where=mean != 0)
        fit_mean_val = np.nanmean(mean)
        fd_mean_val = np.nanmean(fakedata_asimov)
        ratio_means = (fit_mean_val - fd_mean_val) / fd_mean_val if fd_mean_val != 0 else np.nan
        bias_panel_title = (
            f"(Fit Result \u2212 Fake Data Asimov) / Fake Data Asimov  "
            f"[mean: {ratio_means:+.3f}]"
        )
        ax3.step(e, extend(fd_bias), where="post", color=_C["red"],
                  linewidth=1.8, linestyle=(0, (5, 2)), label="Fake Data Asimov / Fit Result")
        ax3.set_title(bias_panel_title, fontsize=8, loc="left")

    ax3.fill_between(e, extend(bias_dn), extend(bias_up),
                      step="post", alpha=0.20, color=_C["orange"], label=r"$\pm1\sigma$ band")
    ax3.step(e, extend(bias), where="post", color=_C["orange"],
              linewidth=1.8, label="MC/Fit")
    ax3.axhline(1, color="black", linewidth=0.8, linestyle="--", alpha=0.7)
    ax3.set_ylabel("MC/Fit")
    ax3.set_xlabel(xlabel)
    ax3.legend(loc="upper right", framealpha=0.9, edgecolor="0.7")

    if xlim_min is not None:
        ax3.set_xlim(left=xlim_min)

    fig.align_ylabels([ax1, ax2, ax3])
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def plot_2d(ax, v, xe, ye, title, cmap="viridis", norm=None):
    mesh = ax.pcolormesh(xe, ye, v.T, shading="auto", cmap=cmap, norm=norm)
    ax.set_title(title)
    ax.set_xlabel("Reconstructed Neutrino Energy (GeV)")
    ax.set_ylabel("Lepton Energy (GeV)")
    return mesh


def save_2d(pdf, v, xe, ye, title, cbar_label, cmap="viridis", norm=None):
    fig, ax = plt.subplots(figsize=(6.5, 5.0))
    mesh = plot_2d(ax, v, xe, ye, title, cmap, norm)
    cb = fig.colorbar(mesh, ax=ax)
    cb.set_label(cbar_label)
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def plot_2d_summaries(f, detectors, pdf):
    log.info("Building 2-D summary pages\u2026")
    for det in detectors:
        node = safe_get(f, f"{det}/2D/summary")
        if node is None:
            continue
        keys = [strip_cycle(k) for k in node.keys()]
        samples = {k.replace("_2D_posteriorMean", "")
                   for k in keys if "_2D_posteriorMean" in k}
        for s in sorted(samples):
            try:
                h_mean = node[f"{s}_2D_posteriorMean"]
                h_post = node[f"{s}_2D_posteriorErr"]
                h_stat = node[f"{s}_2D_statErr"]
            except Exception as exc:
                log.warning("Skipping %s/%s: %s", det, s, exc)
                continue
            mean, xe, ye = get_2d(h_mean)
            post = h_post.values()
            stat = h_stat.values()
            mean = mask_empty(mean)
            post = mask_empty(post)
            stat = mask_empty(stat)
            ratio = np.where(stat > 0, post / stat, np.nan)
            if not np.isfinite(ratio).any():
                log.warning("All-NaN ratio for %s/%s \u2013 skipping.", det, s)
                continue
            rel_post = np.where(mean > 0, post / mean, np.nan)
            rel_stat = np.where(mean > 0, stat / mean, np.nan)
            vmax = np.nanpercentile(ratio, 95)
            norm = mcolors.TwoSlopeNorm(vmin=0, vcenter=1, vmax=max(2.0, float(vmax)))
            save_2d(pdf, mean,     xe, ye, f"{s} Fit Result",         "Events")
            save_2d(pdf, rel_post, xe, ye, f"{s} Posterior \u03c3 / N", "\u03c3 / N")
            save_2d(pdf, rel_stat, xe, ye, f"{s} Stat \u03c3 / N",      "\u03c3 / N")
            save_2d(pdf, ratio,    xe, ye, f"{s} Posterior / Stat \u03c3", "\u03c3_post / \u03c3_stat",
                    cmap="RdBu_r", norm=norm)


def run_single_file_pdf(root_path, data_root_path=None, out_path=None):
    """Generate the standard multi-page PDF for one ROOT file."""
    root_dir = os.path.dirname(os.path.abspath(root_path))
    if data_root_path is None:
        data_root_path = os.path.join(root_dir, "fake_data.root")
    if out_path is None:
        out_path = os.path.join(root_dir, "plots.pdf")

    set_thesis_style()
    f = uproot.open(root_path)

    top_keys = [strip_cycle(k) for k in f.keys(recursive=False)]
    EXCLUDE_TOP = {"Asimov", "AsimovFakeData"}
    detectors = [k for k in top_keys
                 if k not in EXCLUDE_TOP and safe_get(f, k) is not None]
    log.info("Found detectors: %s", detectors)

    asimov_dir = safe_get(f, "Asimov")
    data_hists = load_data_histograms(data_root_path)
    variables = list(VAR_META.keys())

    with PdfPages(out_path) as pdf:
        for det in detectors:
            if asimov_dir is None:
                log.warning("No Asimov directory \u2013 skipping 1-D plots.")
                break

            samples = discover_subsamples(f, det, variables)
            if not samples:
                log.warning("No subsamples found under detector %s \u2013 skipping.", det)
                continue
            log.info("Detector %s: %d subsample(s): %s", det, len(samples), samples)

            for sample in samples:
                sample_label = sample.replace("_", " ").strip()
                for var in variables:
                    h_asimov = find_asimov_hist(f, "Asimov", det, var, sample)
                    if h_asimov is None:
                        log.warning("No Asimov histogram for %s/%s/%s \u2013 skipping.",
                                    det, sample, var)
                        continue
                    if is_2d(h_asimov):
                        continue

                    asimov_vals, edges = get_1d(h_asimov)
                    bin_widths = np.diff(edges)
                    asimov_vals = asimov_vals / bin_widths

                    h_fd = find_asimov_hist(f, "AsimovFakeData", det, var, sample)
                    fakedata_asimov_vals = None
                    if h_fd is not None and not is_2d(h_fd):
                        fd_vals, _ = get_1d(h_fd)
                        fakedata_asimov_vals = fd_vals / bin_widths

                    toys = collect_toys(f, det, var, sample, edges)
                    if not toys:
                        log.warning("No toys found for %s/%s/%s.", det, sample, var)
                        continue

                    toys_arr = np.vstack(toys)
                    mean = np.mean(toys_arr, axis=0)
                    std = np.std(toys_arr, axis=0)

                    xlabel = VAR_META[var]["xlabel"]
                    is_reco = (var == "RecoNeutrinoEnergy")
                    is_true = (var == "TrueNeutrinoEnergy")
                    is_bias = (var == "Enubias")
                    title = f"{var}  \u2014  {det} / {sample_label}"

                    common_kwargs = dict(
                        asimov=asimov_vals, mean=mean, std=std, edges=edges,
                        title=title, xlabel=xlabel, toys=toys_arr,
                        data_hists=data_hists if is_reco else None,
                        is_reco=is_reco, fakedata_asimov=fakedata_asimov_vals,
                    )

                    plot_1d(pdf, **common_kwargs)

                    if is_reco:
                        plot_1d(pdf, **common_kwargs, no_mc_prediction=True, sigma_blue=True)
                        clean_kwargs = dict(common_kwargs)
                        clean_kwargs.update(data_hists=None, fakedata_asimov=None, toys=None,
                                             title=f"{title} (Fit Result vs CV Prediction)")
                        plot_1d(pdf, **clean_kwargs)

                    if is_true:
                        plot_1d(pdf, **common_kwargs, no_fakedata_asimov_line=True,
                                sigma_blue=True, sigma_toy_envelope=True)

                    if is_bias:
                        plot_1d(pdf, **common_kwargs, xlim_min=-0.2)

                    log.info("1-D done: %s / %s / %s  (%d toys)", det, sample, var, len(toys))

        plot_2d_summaries(f, detectors, pdf)

    log.info("Saved \u2192 %s", out_path)


# ============================================================================
# 4. MULTI-FIT COMPARISON OVERLAY (multiple run plans on one figure)
# ============================================================================

@dataclass
class FitResult:
    """One posterior-predictive fit to overlay on the comparison plot."""
    label: str
    bin_edges: np.ndarray
    cv: np.ndarray
    fit: np.ndarray
    fit_sigma: Optional[np.ndarray] = None
    fit_p05: Optional[np.ndarray] = None
    fit_p95: Optional[np.ndarray] = None
    color: Optional[str] = None
    show_band: bool = True

    @property
    def bin_centers(self) -> np.ndarray:
        return 0.5 * (self.bin_edges[:-1] + self.bin_edges[1:])

    @classmethod
    def from_root(cls, root_path, det, var, sample=None, label=None, **kwargs):
        """Build a FitResult directly from a ROOT file via extract_fit_result."""
        r = extract_fit_result(root_path, det, var, sample=sample)
        return cls(
            label=label or r["sample"], bin_edges=r["edges"], cv=r["cv"],
            fit=r["fit"], fit_sigma=r["fit_sigma"],
            fit_p05=r["fit_p05"], fit_p95=r["fit_p95"], **kwargs,
        )


def _step_xy(edges, values):
    x = np.repeat(edges, 2)[1:-1]
    y = np.repeat(values, 2)
    return x, y


def _band_step(edges, lo, hi):
    x = np.repeat(edges, 2)[1:-1]
    return x, np.repeat(lo, 2), np.repeat(hi, 2)


def plot_comparison(
    fits: Sequence[FitResult],
    xlabel: str = "",
    ylabel: str = "Events",
    title: str = "",
    reference: Optional[FitResult] = None,
    show_cv: bool = True,
    colors: Optional[Sequence[str]] = None,
    figsize=(8, 9),
    save_path: Optional[str] = None,
):
    """Build the 3-panel run-plan comparison figure."""
    if colors is None:
        colors = PALETTE
    for i, fit in enumerate(fits):
        if fit.color is None:
            fit.color = colors[i % len(colors)]
    if reference is None:
        reference = fits[0]

    fig, axes = plt.subplots(
        3, 1, figsize=figsize, sharex=True,
        gridspec_kw={"height_ratios": [3, 1.1, 1.1], "hspace": 0.08},
    )
    ax_main, ax_prec, ax_ratio = axes

    for fit in fits:
        x, y_fit = _step_xy(fit.bin_edges, fit.fit)
        ax_main.plot(x, y_fit, color=fit.color, lw=1.8, label=fit.label, zorder=3)

        if show_cv:
            x_cv, y_cv = _step_xy(fit.bin_edges, fit.cv)
            ax_main.plot(x_cv, y_cv, color=fit.color, lw=1.1, ls=(0, (4, 2)),
                         alpha=0.75, zorder=2)

        if fit.show_band:
            if fit.fit_p05 is not None and fit.fit_p95 is not None:
                xb, lo, hi = _band_step(fit.bin_edges, fit.fit_p05, fit.fit_p95)
                ax_main.fill_between(xb, lo, hi, color=fit.color, alpha=0.12,
                                      linewidth=0, zorder=1)
            if fit.fit_sigma is not None:
                xb, lo, hi = _band_step(fit.bin_edges, fit.fit - fit.fit_sigma,
                                         fit.fit + fit.fit_sigma)
                ax_main.fill_between(xb, lo, hi, color=fit.color, alpha=0.25,
                                      linewidth=0, zorder=1.5)

    ax_main.set_ylabel(ylabel)
    if title:
        ax_main.set_title(title, fontsize=15, pad=10)
    ax_main.set_ylim(bottom=0)
    ax_main.legend(loc="upper right", ncol=1)

    any_sigma = False
    for fit in fits:
        if fit.fit_sigma is None:
            continue
        any_sigma = True
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(fit.fit != 0, fit.fit_sigma / fit.fit, np.nan)
        x, y = _step_xy(fit.bin_edges, ratio)
        ax_prec.plot(x, y, color=fit.color, lw=1.6)
    ax_prec.axhline(0, color="0.6", lw=0.8, ls=":")
    ax_prec.set_ylabel(r"$\sigma_{\mathrm{Fit}}/\mu$")
    if not any_sigma:
        ax_prec.text(0.5, 0.5, "no fit_sigma provided", transform=ax_prec.transAxes,
                     ha="center", va="center", color="0.5", fontsize=10)

    ref_x, ref_y = reference.bin_centers, reference.fit
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
    for ax in axes:
        ax.tick_params(which="both", direction="in", top=True, right=True)
    fig.align_ylabels(axes)

    if save_path:
        fig.savefig(save_path)
    return fig, axes


def run_comparison(run_plan_files, det, var, sample=None, reference_label=None,
                    save_path="run_plan_comparison.png"):
    """
    Convenience wrapper: run_plan_files is {label: root_path, ...}.
    Builds one FitResult per file and plots the comparison.

    A file's value can also be a dict to override det/sample for just that
    entry, e.g. when some run plans only have an FD directory (no ND):

        run_plan_files = {
            "On Axis only": {"path": "...withtoysFD1.root", "det": "FD"},
            "75% On Axis + 25% Off Axis": "...withtoys.root",  # uses default det
            ...
        }

    When overriding with a dict and `sample` isn't given, it falls back
    to auto-detection (works when the file has exactly one sample for
    that det/var, as with a flat single-sample FD file).
    """
    set_thesis_style()
    fits = []
    for label, entry in run_plan_files.items():
        if isinstance(entry, dict):
            path = entry["path"]
            this_det = entry.get("det", det)
            this_sample = entry.get("sample", None)  # None -> auto-detect for this file
        else:
            path = entry
            this_det = det
            this_sample = sample
        try:
            fits.append(FitResult.from_root(path, det=this_det, var=var,
                                             sample=this_sample, label=label))
        except Exception as exc:
            log.warning("Skipping '%s': %s", label, exc)
    if not fits:
        raise RuntimeError("No fits could be loaded -- check run_plan_files paths/sample.")

    reference = fits[0]
    if reference_label:
        for fit in fits:
            if fit.label == reference_label:
                reference = fit
                break

    xlabel = VAR_META.get(var, {}).get("xlabel", var)
    fig, axes = plot_comparison(
        fits, xlabel=xlabel, ylabel="Events / GeV",
        title=f"{det} {var} \u2014 Run Plan Comparison",
        reference=reference, save_path=save_path,
    )
    log.info("Saved \u2192 %s", save_path)
    return fig, axes


# ============================================================================
# 5. CLI
# ============================================================================

def main():
    logging.basicConfig(level=logging.INFO, format="%(levelname)s  %(message)s")

    if len(sys.argv) < 2:
        print("Usage:")
        print("  python plot_posterior_all.py file.root [data.root] [output.pdf]")
        print("  python plot_posterior_all.py --compare   # edit run_plan_files below first")
        sys.exit(1)

    if sys.argv[1] == "--compare":
        # ---- EDIT THIS: your six run plans, one FD ROOT path each ----
        base = "/home/abipeake/scratch/July2026"
        run_plan_files = {
            "On Axis only":
                #f"/home/abipeake/scratch/March2026/July8th_onaxis_fullsysts_somemc/Posteriorpredictive_out_posteriorpredictive_withtoys.root",
                #f"{base}/Template_0m_asimov_post/Posteriorpredictive_out_posteriorpredictive_withtoysFD1.root",
                f"/home/abipeake/scratch/July2026/OnAxis_templatecovmat_new/Posteriorpredictive_out_posteriorpredictive_withtoys.root",
            "50% On Axis + 50% Off Axis":
                f"{base}/Asimov_allsystematics_onplusoffaxis_postadaptive/Posteriorpredictive_out_posteriorpredictive_withtoysFD.root",
            "75% On Axis + 25% Off Axis":
                f"{base}/Runplan_75onaxis_25offaxis_post/Posteriorpredictive_out_posteriorpredictive_withtoys.root",
            "25% On Axis + 75% Off Axis":
                f"{base}/Runplan_25onaxis_75offaxis_post/Posteriorpredictive_out_posteriorpredictive_withtoys.root",
            "Equal time at each position":
                f"{base}/Runplan_equal_allposition_post/Posteriorpredictive_out_posteriorpredictive_withtoysFD.root",
            "50% On Axis + 50% Off Axis (12m)":
                f"{base}/Runplan_50onaxis_50upto12mpost/Posteriorpredictive_out_posteriorpredictive_withtoys.root",
        }

        # One comparison figure per variable, each saved as its own PDF.
        for var in ("TrueNeutrinoEnergy", "RecoNeutrinoEnergy", "Enubias"):
            save_path = f"run_plan_comparison_{var}.pdf"
            log.info("Building comparison for %s ...", var)
            try:
                run_comparison(
                    run_plan_files, det="FD", var=var,
                    sample=None,   # auto-detected: each FD file has exactly one sample
                    save_path=save_path,
                )
            except Exception as exc:
                log.error("Failed to build comparison for %s: %s", var, exc)
        return

    root_path = sys.argv[1]
    root_dir = os.path.dirname(os.path.abspath(root_path))
    data_root_path = os.path.join(root_dir, "fake_data.root")
    out_path = os.path.join(root_dir, "plots.pdf")
    for arg in sys.argv[2:]:
        if arg.endswith(".root"):
            data_root_path = arg
        elif arg.endswith(".pdf"):
            out_path = arg

    run_single_file_pdf(root_path, data_root_path, out_path)


if __name__ == "__main__":
    main()
