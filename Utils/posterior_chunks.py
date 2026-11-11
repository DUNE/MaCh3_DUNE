#!/usr/bin/env python3
"""
Posterior predictive plots, split into N consecutive chunks of toy index,
from a ROOT file produced by the neutrino fit.

Toy index (posterior_toy_NNN) is monotonically increasing with MCMC chain
step (see MakeSpectaVariations / the sort of chosenEntries in the C++ code),
so "1st chunk" = earliest post-burn-in steps, "last chunk" = latest steps.
This is a chain-position split, NOT a split by RNG draw order (that
information is discarded before the toys are written to file).

Usage
-----
    python plot_posterior_chunks.py file.root [--data data.root] \
        [--out plots.pdf] [--chunk-size 10000] [--n-chunks 5]
"""

import os
import re
import sys
import argparse
import logging
import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.backends.backend_pdf import PdfPages

log = logging.getLogger(__name__)

# ── thesis-quality style ────────────────────────────────────────────────────
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
# CHUNK_COLORS = [_C["blue"], _C["orange"], _C["green"], _C["purple"], _C["red"],
#                 _C["sky"], _C["yellow"]]
CHUNK_COLORS = [plt.cm.tab10(i) for i in range(20)]
THESIS_RC = {
    "font.family":        "serif",
    "font.serif":         ["STIXTwo Text", "DejaVu Serif", "Times New Roman"],
    "mathtext.fontset":   "stix",
    "font.size":          11,
    "axes.titlesize":     11,
    "axes.labelsize":     11,
    "xtick.labelsize":    10,
    "ytick.labelsize":    10,
    "legend.fontsize":    8,
    "pdf.fonttype":       42,
    "ps.fonttype":        42,
    "figure.dpi":         150,
    "lines.linewidth":    1.6,
    "xtick.direction":    "in",
    "ytick.direction":    "in",
    "xtick.top":          True,
    "ytick.right":        True,
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "axes.grid":          True,
    "grid.alpha":         0.3,
    "grid.linewidth":     0.6,
    "axes.linewidth":     0.8,
}
plt.rcParams.update(THESIS_RC)

VAR_META = {
    "TrueNeutrinoEnergy":  {"xlabel": "True Neutrino Energy (GeV)"},
    "RecoNeutrinoEnergy":  {"xlabel": "Reconstructed Neutrino Energy (GeV)"},
    "Enubias":             {"xlabel": "E_{ν} bias (GeV)"},
}

TOY_IDX_RE = re.compile(r"_posterior_toy_(\d+)$")

# ── helpers ──────────────────────────────────────────────────────────────────

def safe_get(f, path):
    try:
        return f[path]
    except Exception as exc:
        log.debug("safe_get(%s): %s", path, exc)
        return None


def is_2d(h):
    return len(h.axes) == 2


def get_1d(h):
    return h.values(), h.axes[0].edges()


def extend(v):
    return np.append(v, v[-1])


def strip_cycle(key):
    return key.split(";")[0]


def toy_index(key):
    """Extract the integer toy index from a histogram key, or None."""
    m = TOY_IDX_RE.search(key)
    return int(m.group(1)) if m else None


# ── toy collection ──────────────────────────────────────────────────────────

def collect_toys_by_index(f, det, var, edges, toy_marker):
    """
    Return dict {toy_index: values_per_bin_width} for every toy histogram
    found for this detector/variable, read from wherever it lives in the file.
    """
    out = {}

    def _try_add(h, key):
        idx = toy_index(key)
        if idx is None:
            log.debug("Could not parse toy index from %s - skipped.", key)
            return
        try:
            if is_2d(h):
                return
            v, e = get_1d(h)
            if not np.allclose(e, edges):
                log.debug("Edge mismatch for toy %s - skipped.", key)
                return
            out[idx] = v / np.diff(e)
        except Exception as exc:
            log.debug("Could not read toy %s: %s", key, exc)

    var_subdir = safe_get(f, f"{det}/{var}")
    if var_subdir is not None:
        for k in var_subdir.keys(recursive=False):
            ks = strip_cycle(k)
            if toy_marker in ks:
                _try_add(var_subdir[k], ks)

    det_dir = safe_get(f, det)
    if det_dir is not None:
        for k in det_dir.keys(recursive=False):
            ks = strip_cycle(k)
            if toy_marker in ks:
                _try_add(det_dir[k], ks)

    return out


def build_chunks(toys_by_idx, chunk_size, n_chunks):
    """
    Split toys into n_chunks consecutive ranges of chunk_size indices each:
    [0, chunk_size), [chunk_size, 2*chunk_size), ... last chunk.
    Returns list of (label, np.ndarray[n_toys_in_chunk, n_bins]) - chunks
    with zero toys found are skipped with a warning.
    """
    if not toys_by_idx:
        return []

    all_idx = sorted(toys_by_idx.keys())
    max_idx = all_idx[-1]

    chunks = []
    for c in range(n_chunks):
        lo = c * chunk_size
        hi = (c + 1) * chunk_size - 1
        if c == n_chunks - 1:
            # last chunk: extend to whatever the highest available index is,
            # so "last 10k" still works if the total isn't an exact multiple
            hi = max(hi, max_idx)
        sel = [toys_by_idx[i] for i in all_idx if lo <= i <= hi]
        if not sel:
            log.warning("No toys found in range [%d, %d] - skipping chunk.", lo, hi)
            continue
        label = f"Throws {lo}\u2013{hi}"
        chunks.append((label, np.vstack(sel)))

    return chunks


def _percentile_band(arr, percentiles):
    lo = np.percentile(arr, percentiles[0], axis=0)
    hi = np.percentile(arr, percentiles[1], axis=0)
    return lo, hi


# ── plotting ─────────────────────────────────────────────────────────────────

def plot_1d_chunks(pdf, asimov, chunks, edges, title, xlabel,
                    fakedata_asimov=None, xlim_min=None,
                    percentiles=(5, 95), show_percentiles=True):
    """
    Top panel: mean +/- 1 sigma envelope per chunk, overlaid, with an
    additional (lighter) percentile envelope per chunk if requested.
    Bottom panel: each chunk's mean relative to the full-sample (all chunks
    combined) mean, to visualise drift/convergence across the chain.
    """
    e = edges
    a = extend(asimov)

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(6.5, 6.5),
        gridspec_kw={"height_ratios": [3, 1.2]},
        sharex=True,
    )
    fig.subplots_adjust(hspace=0.06)

    # overall mean across all chunks combined, used as the reference line
    combined = np.vstack([arr for _, arr in chunks])
    overall_mean = np.mean(combined, axis=0)

    for (label, arr), color in zip(chunks, CHUNK_COLORS):
        mean = np.mean(arr, axis=0)
        std = np.std(arr, axis=0)
        m = extend(mean)
        u = extend(mean + std)
        l = extend(mean - std)

        if show_percentiles:
            plo, phi = _percentile_band(arr, percentiles)
            ax1.fill_between(e, extend(plo), extend(phi), step="post",
                              alpha=0.08, color=color)

        ax1.fill_between(e, l, u, step="post", alpha=0.18, color=color)
        ax1.step(e, m, where="post", linewidth=1.6, color=color,
                  label=f"{label} (n={arr.shape[0]})")

        rel = np.divide(mean, overall_mean,
                         out=np.ones_like(mean), where=overall_mean != 0)
        ax2.step(e, extend(rel), where="post", linewidth=1.6, color=color)

    if not np.allclose(asimov, 0):
        ax1.step(e, a, where="post", label="CV Prediction",
                  linewidth=1.6, color=_C["black"], linestyle=":")

    if fakedata_asimov is not None:
        ax1.step(e, extend(fakedata_asimov), where="post",
                  label="Fake Data ", linewidth=1.6,
                  color=_C["black"], linestyle=(0, (5, 2)))

    ax1.set_ylabel("Events / GeV")
    ax1.set_title(title)

    handles, labels = ax1.get_legend_handles_labels()
    if show_percentiles:
        proxy = Patch(facecolor=_C["black"], alpha=0.08,
                       label=f"{percentiles[0]}\u2013{percentiles[1]}th pct (per chunk)")
        handles.append(proxy)
        labels.append(proxy.get_label())
    proxy_sigma = Patch(facecolor=_C["black"], alpha=0.18, label=r"$\pm1\sigma$ (per chunk)")
    handles.append(proxy_sigma)
    labels.append(proxy_sigma.get_label())
    ax1.legend(handles, labels, loc="best", framealpha=0.9, edgecolor="0.7", ncol=1)
    plt.setp(ax1.get_xticklabels(), visible=False)

    ax2.axhline(1, color="black", linewidth=0.7, linestyle="--", alpha=0.6)
    ax2.set_ylabel("Chunk mean /\nfull-sample mean")
    ax2.set_xlabel(xlabel)

    if xlim_min is not None:
        ax2.set_xlim(left=xlim_min)

    fig.align_ylabels([ax1, ax2])
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def plot_1d_single_chunk(pdf, asimov, label, arr, edges, title, xlabel,
                          fakedata_asimov=None, xlim_min=None, color=_C["blue"],
                          percentiles=(5, 95), show_percentiles=True):
    """
    One page per chunk: mean +/- 1 sigma band (and optional percentile band)
    for that chunk only, against the CV prediction and (if available) the
    fake-data Asimov, plus a bottom MC/Fit-style ratio panel.
    """
    e = edges
    mean = np.mean(arr, axis=0)
    std = np.std(arr, axis=0)
    m = extend(mean)
    u = extend(mean + std)
    l = extend(mean - std)
    a = extend(asimov)

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(6.5, 6.5),
        gridspec_kw={"height_ratios": [3, 1.2]},
        sharex=True,
    )
    fig.subplots_adjust(hspace=0.06)

    if show_percentiles:
        plo, phi = _percentile_band(arr, percentiles)
        ax1.fill_between(e, extend(plo), extend(phi), step="post",
                          alpha=0.12, color=color,
                          label=f"{percentiles[0]}\u2013{percentiles[1]}th pct")

    ax1.fill_between(e, l, u, step="post", alpha=0.25, color=color,
                      label=r"$\pm1\sigma$ (posterior)")
    ax1.step(e, m, where="post", linewidth=1.8, color=color,
              label=f"Fit Result ({label}, n={arr.shape[0]})")

    if not np.allclose(asimov, 0):
        ax1.step(e, a, where="post", label="CV Prediction",
                  linewidth=1.6, color=_C["orange"])

    if fakedata_asimov is not None:
        ax1.step(e, extend(fakedata_asimov), where="post",
                  label="Fake Data ", linewidth=1.6,
                  color=_C["black"], linestyle=(0, (5, 2)))

    ax1.set_ylabel("Events / GeV")
    ax1.set_title(f"{title}  [{label}]")
    ax1.legend(loc="best", framealpha=0.9, edgecolor="0.7")
    plt.setp(ax1.get_xticklabels(), visible=False)

    bias = np.divide(asimov, mean, out=np.ones_like(asimov), where=mean != 0)
    bias_up = np.divide(asimov, mean + std, out=np.ones_like(asimov), where=(mean + std) != 0)
    bias_dn = np.divide(asimov, mean - std, out=np.ones_like(asimov), where=(mean - std) != 0)

    ax2.fill_between(e, extend(bias_dn), extend(bias_up),
                      step="post", alpha=0.20, color=_C["orange"])
    ax2.step(e, extend(bias), where="post", color=_C["orange"], linewidth=1.6, label="CV/Fit")
    ax2.axhline(1, color="black", linewidth=0.7, linestyle="--", alpha=0.6)
    ax2.set_ylabel("CV/Fit")
    ax2.set_xlabel(xlabel)
    ax2.legend(loc="upper right", framealpha=0.9, edgecolor="0.7")

    if xlim_min is not None:
        ax2.set_xlim(left=xlim_min)

    fig.align_ylabels([ax1, ax2])
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root_path", help="Input ROOT file with posterior toys")
    parser.add_argument("--data", dest="data_root_path", default=None,
                         help="Optional fake-data ROOT file (unused here, kept for parity)")
    parser.add_argument("--out", dest="out_path", default=None,
                         help="Output PDF path (default: <input dir>/plots_chunks.pdf)")
    parser.add_argument("--chunk-size", type=int, default=10000,
                         help="Number of toys per chunk (default: 10000)")
    parser.add_argument("--n-chunks", type=int, default=5,
                         help="Number of chunks (default: 5)")
    parser.add_argument("--mode", choices=["overlay", "separate", "both"],
                         default="both",
                         help="'overlay' = one page with all chunks overlaid; "
                              "'separate' = one page per chunk; "
                              "'both' = overlay page then separate pages (default)")
    parser.add_argument("--percentiles", type=float, nargs=2, default=[5, 95],
                         metavar=("LOW", "HIGH"),
                         help="Percentile band to shade in addition to +/-1 sigma "
                              "(default: 5 95)")
    parser.add_argument("--no-percentiles", action="store_true",
                         help="Disable the percentile band entirely")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(levelname)s  %(message)s")

    root_dir = os.path.dirname(os.path.abspath(args.root_path))
    out_path = args.out_path or os.path.join(root_dir, "plots_chunks.pdf")

    f = uproot.open(args.root_path)

    top_keys = [strip_cycle(k) for k in f.keys(recursive=False)]
    detectors = [k for k in top_keys if k not in ("Asimov", "AsimovFakeData")
                 and safe_get(f, k) is not None]
    log.info("Found detectors: %s", detectors)

    asimov_dir = safe_get(f, "Asimov")
    fakedata_asimov_dir = safe_get(f, "AsimovFakeData")

    variables = list(VAR_META.keys())

    with PdfPages(out_path) as pdf:
        for det in detectors:
            for var in variables:

                if asimov_dir is None:
                    log.warning("No Asimov directory in file - skipping 1-D plots.")
                    break

                asimov_keys = [
                    strip_cycle(k) for k in asimov_dir.keys(recursive=False)
                    if strip_cycle(k).startswith(det + "_")
                    and f"_{var}_Asimov" in strip_cycle(k)
                ]
                if not asimov_keys:
                    log.warning("No Asimov histogram for %s/%s - skipping.", det, var)
                    continue

                h_asimov = asimov_dir[asimov_keys[0]]
                if is_2d(h_asimov):
                    log.warning("Asimov hist for %s/%s is 2-D - skipping.", det, var)
                    continue

                asimov_vals, edges = get_1d(h_asimov)
                bin_widths = np.diff(edges)
                asimov_vals = asimov_vals / bin_widths

                fakedata_asimov_vals = None
                if fakedata_asimov_dir is not None:
                    fd_keys = [
                        strip_cycle(k) for k in fakedata_asimov_dir.keys(recursive=False)
                        if strip_cycle(k).startswith(det + "_")
                        and f"_{var}_AsimovFakeData" in strip_cycle(k)
                    ]
                    if fd_keys:
                        h_fd = fakedata_asimov_dir[fd_keys[0]]
                        if not is_2d(h_fd):
                            fd_vals, _ = get_1d(h_fd)
                            fakedata_asimov_vals = fd_vals / bin_widths

                toy_marker = f"_{var}_posterior_toy_"
                toys_by_idx = collect_toys_by_index(f, det, var, edges, toy_marker)

                if not toys_by_idx:
                    log.warning("No toys found for %s/%s.", det, var)
                    continue

                chunks = build_chunks(toys_by_idx, args.chunk_size, args.n_chunks)
                if not chunks:
                    log.warning("No usable chunks for %s/%s.", det, var)
                    continue

                xlabel = VAR_META[var]["xlabel"]
                is_bias = (var == "Enubias")
                xlim_min = -0.2 if is_bias else None

                if args.mode in ("overlay", "both"):
                    plot_1d_chunks(
                        pdf, asimov=asimov_vals, chunks=chunks, edges=edges,
                        title=f"{det} / {var}", xlabel=xlabel,
                        fakedata_asimov=fakedata_asimov_vals,
                        xlim_min=xlim_min,
                        percentiles=tuple(args.percentiles),
                        show_percentiles=not args.no_percentiles,
                    )

                if args.mode in ("separate", "both"):
                    for (label, arr), color in zip(chunks, CHUNK_COLORS):
                        plot_1d_single_chunk(
                            pdf, asimov=asimov_vals, label=label, arr=arr,
                            edges=edges, title=f"{det} / {var}", xlabel=xlabel,
                            fakedata_asimov=fakedata_asimov_vals,
                            xlim_min=xlim_min, color=color,
                            percentiles=tuple(args.percentiles),
                            show_percentiles=not args.no_percentiles,
                        )

                log.info("Done: %s / %s  (%d chunks, %d total toys)",
                          det, var, len(chunks), sum(a.shape[0] for _, a in chunks))

    log.info("Saved -> %s", out_path)


if __name__ == "__main__":
    main()
