#!/usr/bin/env python3
"""
Thesis-subfigure export: for each detector position (subsample), write out
exactly two plots —

    1. RecoNeutrinoEnergy  (Fit Result, CV Prediction, Fake Data Asimov,
                             ±1σ posterior band, 5–95th pct toy envelope,
                             data points if available)
    2. TrueNeutrinoEnergy  (Fit Result, CV Prediction, ±1σ posterior band,
                             5–95th pct toy envelope, per-bin error bars;
                             no Fake Data Asimov line)

named sensibly by position, so they can be dropped straight into a thesis
as subfigures. No Enubias plot, no extra variants, no combined multi-page
PDF, no 2-D summary pages — just these two per position.

Output layout
-------------
    <output_dir>/<detector>/<position>__Reco.pdf
    <output_dir>/<detector>/<position>__True.pdf

Usage
-----
    python thesis_figures.py file.root [data.root] [output_dir]

    - file.root    : required, the fit-output ROOT file
    - data.root    : optional, defaults to "<file.root's dir>/fake data.root"
    - output_dir   : optional, defaults to "<file.root's dir>/thesis_figures/"
"""

import os
import re
import sys
import logging
import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt

log = logging.getLogger(__name__)

# ── configuration ───────────────────────────────────────────────────────────

TOY_VARIANT = "combined"

_C = {
    "blue":   "#0072B2",
    "orange": "#E69F00",
    "red":    "#D55E00",
    "green":  "#009E73",
    "purple": "#CC79A7",
    "sky":    "#56B4E9",
    "yellow": "#F0E442",
    "black":  "#000000",
}

THESIS_RC = {
    "font.family":        "serif",
    "font.serif":         ["STIXTwo Text", "DejaVu Serif", "Times New Roman"],
    "mathtext.fontset":   "stix",
    "font.size":          11,
    "axes.titlesize":     11,
    "axes.labelsize":     11,
    "xtick.labelsize":    10,
    "ytick.labelsize":    10,
    "legend.fontsize":    9,
    "pdf.fonttype":       42,
    "ps.fonttype":        42,
    "figure.dpi":         150,
    "lines.linewidth":    1.8,
    "patch.linewidth":    1.0,
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
    "TrueNeutrinoEnergy": {"xlabel": "True Neutrino Energy (GeV)"},
    "RecoNeutrinoEnergy": {"xlabel": "Reconstructed Neutrino Energy (GeV)"},
}

# ── helpers ──────────────────────────────────────────────────────────────────

def slugify(text):
    text = re.sub(r"[^\w\-]+", "_", text)
    return re.sub(r"_+", "_", text).strip("_")


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


def find_hist_in_dir(node, prefix):
    if node is None:
        return None
    matches = [k for k in node.keys(recursive=False)
               if strip_cycle(k).startswith(prefix)]
    if not matches:
        return None
    if len(matches) > 1:
        log.warning("Multiple keys match prefix '%s': %s -- using first",
                     prefix, matches)
    return node[matches[0]]


def discover_subsamples(f, det, variables):
    samples = set()
    for var in variables:
        node = safe_get(f, f"{det}/{var}")
        if node is None:
            continue
        for k in node.keys(recursive=False):
            samples.add(strip_cycle(k))
    return sorted(samples)


def load_data_histograms(data_root_path):
    if data_root_path is None:
        return {}
    try:
        df = uproot.open(data_root_path)
    except Exception as exc:
        log.error("Could not open data file %s: %s", data_root_path, exc)
        return {}

    sample_folder = safe_get(df, "SampleFolder")
    node = sample_folder if sample_folder is not None else df
    prefix_required = sample_folder is not None

    data_hists = {}
    for k in node.keys(recursive=False):
        ks = strip_cycle(k)
        if prefix_required and not ks.startswith("data_"):
            continue
        try:
            h = node[k]
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


def _toy_envelope(ax, toys, edges, percentiles=(5, 95)):
    lo = np.percentile(toys, percentiles[0], axis=0)
    hi = np.percentile(toys, percentiles[1], axis=0)
    ax.fill_between(edges, extend(lo), extend(hi),
                    step="post", alpha=0.12, color=_C["black"],
                    label=f"{percentiles[0]}–{percentiles[1]}th pct")


# ── the two plot styles ─────────────────────────────────────────────────────

def plot_reco(output_dir, filename, asimov, mean, std, edges, title, xlabel,
              toys, fakedata_asimov, data_hists):
    """Reco-style figure: Fit Result, CV Prediction, Fake Data Asimov,
    ±1σ band, 5–95th pct toy envelope, data points if available."""

    fig, (ax1, ax2, ax3) = plt.subplots(
        3, 1, figsize=(6.5, 7.5),
        gridspec_kw={"height_ratios": [3, 1, 1]}, sharex=True,
    )
    fig.subplots_adjust(hspace=0.05)

    e = edges
    m, u, l, a = extend(mean), extend(mean + std), extend(mean - std), extend(asimov)

    _toy_envelope(ax1, toys, e)
    ax1.fill_between(e, l, u, step="post", alpha=0.25, color=_C["blue"],
                     label=r"$\pm1\sigma$ (posterior)")
    ax1.step(e, m, where="post", label="Fit Result", linewidth=1.8, color=_C["blue"])
    ax1.step(e, a, where="post", label="CV Prediction", linewidth=1.8, color=_C["orange"])

    if fakedata_asimov is not None:
        ax1.step(e, extend(fakedata_asimov), where="post", label="Fake Data Asimov",
                 linewidth=1.8, color=_C["red"], linestyle=(0, (5, 2)))

    if data_hists:
        colors  = [_C["green"], _C["purple"], _C["sky"], _C["yellow"]]
        markers = ["o", "s", "^", "D"]
        for idx, (name, (dvals, dedges)) in enumerate(data_hists.items()):
            centres = 0.5 * (dedges[:-1] + dedges[1:])
            ax1.errorbar(
                centres, dvals,
                yerr=np.sqrt(np.maximum(dvals * np.diff(dedges), 0)) / np.diff(dedges),
                fmt=markers[idx % len(markers)], markersize=4, linewidth=1.2,
                capsize=3, label=_friendly_label(name),
                color=colors[idx % len(colors)], zorder=5,
            )

    ax1.set_ylabel("Events / GeV")
    ax1.set_title(title)
    ax1.legend(loc="best", framealpha=0.9, edgecolor="0.7")
    plt.setp(ax1.get_xticklabels(), visible=False)

    rel_unc = np.divide(std, mean, out=np.zeros_like(std), where=mean != 0)
    ax2.step(e, extend(rel_unc), where="post", color=_C["red"],
             label=r"$\sigma_\mathrm{Fit}\,/\,\mu$", linewidth=1.8)
    ax2.axhline(0, color="black", linewidth=0.6, linestyle=":", alpha=0.5)
    ax2.set_ylabel(r"$\sigma_\mathrm{Fit}\,/\,\mu$")
    ax2.legend(loc="upper right", framealpha=0.9, edgecolor="0.7")
    plt.setp(ax2.get_xticklabels(), visible=False)

    bias    = np.divide(asimov, mean,       out=np.ones_like(asimov), where=mean != 0)
    bias_up = np.divide(asimov, mean + std, out=np.ones_like(asimov), where=(mean + std) != 0)
    bias_dn = np.divide(asimov, mean - std, out=np.ones_like(asimov), where=(mean - std) != 0)

    if fakedata_asimov is not None:
        fd_bias = np.divide(fakedata_asimov, mean, out=np.ones_like(fakedata_asimov), where=mean != 0)
        fit_mean_val = np.nanmean(mean)
        fd_mean_val  = np.nanmean(fakedata_asimov)
        ratio_means  = (fit_mean_val - fd_mean_val) / fd_mean_val if fd_mean_val != 0 else np.nan
        ax3.step(e, extend(fd_bias), where="post", color=_C["red"], linewidth=1.8,
                 linestyle=(0, (5, 2)), label="Fake Data Asimov / Fit Result")
        ax3.set_title(f"(Fit Result \u2212 Fake Data Asimov) / Fake Data Asimov  "
                      f"[mean: {ratio_means:+.3f}]", fontsize=8, loc="left")

    ax3.fill_between(e, extend(bias_dn), extend(bias_up), step="post", alpha=0.20,
                     color=_C["orange"], label=r"$\pm1\sigma$ band")
    ax3.step(e, extend(bias), where="post", color=_C["orange"], linewidth=1.8, label="MC/Fit")
    ax3.axhline(1, color="black", linewidth=0.8, linestyle="--", alpha=0.7)
    ax3.set_ylabel("MC/Fit")
    ax3.set_xlabel(xlabel)
    ax3.legend(loc="upper right", framealpha=0.9, edgecolor="0.7")

    fig.align_ylabels([ax1, ax2, ax3])
    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, filename + ".pdf")
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    log.info("Wrote %s", out_path)


def plot_true(output_dir, filename, asimov, mean, std, edges, title, xlabel, toys):
    """True-style figure: Fit Result, CV Prediction, ±1σ band, 5–95th pct
    toy envelope, per-bin error bars. No Fake Data Asimov line."""

    fig, (ax1, ax2, ax3) = plt.subplots(
        3, 1, figsize=(6.5, 7.5),
        gridspec_kw={"height_ratios": [3, 1, 1]}, sharex=True,
    )
    fig.subplots_adjust(hspace=0.05)

    e = edges
    m, u, l, a = extend(mean), extend(mean + std), extend(mean - std), extend(asimov)
    centres = 0.5 * (e[:-1] + e[1:])

    _toy_envelope(ax1, toys, e)
    ax1.fill_between(e, l, u, step="post", alpha=0.25, color=_C["blue"],
                     label=r"$\pm1\sigma$ (posterior)")
    ax1.step(e, m, where="post", label="Fit Result", linewidth=1.8, color=_C["blue"])
    ax1.step(e, a, where="post", label="CV Prediction", linewidth=1.8, color=_C["orange"])
    ax1.errorbar(centres, mean, yerr=std, fmt="none", ecolor=_C["blue"],
                 elinewidth=1.2, capsize=3, alpha=0.9, zorder=6, label="_nolegend_")

    ax1.set_ylabel("Events / GeV")
    ax1.set_title(title)
    ax1.legend(loc="best", framealpha=0.9, edgecolor="0.7")
    plt.setp(ax1.get_xticklabels(), visible=False)

    rel_unc = np.divide(std, mean, out=np.zeros_like(std), where=mean != 0)
    toy_rel_each = np.divide(np.std(toys, axis=0, keepdims=True) * np.ones_like(toys),
                              np.maximum(toys, 1e-30))
    lo = np.percentile(toy_rel_each, 5, axis=0)
    hi = np.percentile(toy_rel_each, 95, axis=0)
    ax2.fill_between(e, extend(lo), extend(hi), step="post", alpha=0.15,
                     color="gray", label="5–95th pct")
    ax2.step(e, extend(rel_unc), where="post", color=_C["blue"],
             label=r"$\sigma_\mathrm{Fit}\,/\,\mu$", linewidth=1.8)
    ax2.axhline(0, color="black", linewidth=0.6, linestyle=":", alpha=0.5)
    ax2.set_ylabel(r"$\sigma_\mathrm{Fit}\,/\,\mu$")
    ax2.legend(loc="upper right", framealpha=0.9, edgecolor="0.7")
    plt.setp(ax2.get_xticklabels(), visible=False)

    bias    = np.divide(asimov, mean,       out=np.ones_like(asimov), where=mean != 0)
    bias_up = np.divide(asimov, mean + std, out=np.ones_like(asimov), where=(mean + std) != 0)
    bias_dn = np.divide(asimov, mean - std, out=np.ones_like(asimov), where=(mean - std) != 0)

    ax3.fill_between(e, extend(bias_dn), extend(bias_up), step="post", alpha=0.20,
                     color=_C["orange"], label=r"$\pm1\sigma$ band")
    ax3.step(e, extend(bias), where="post", color=_C["orange"], linewidth=1.8, label="MC/Fit")
    ax3.axhline(1, color="black", linewidth=0.8, linestyle="--", alpha=0.7)
    ax3.set_ylabel("MC/Fit")
    ax3.set_xlabel(xlabel)
    ax3.legend(loc="upper right", framealpha=0.9, edgecolor="0.7")

    fig.align_ylabels([ax1, ax2, ax3])
    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, filename + ".pdf")
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    log.info("Wrote %s", out_path)


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    logging.basicConfig(level=logging.INFO, format="%(levelname)s  %(message)s")

    if len(sys.argv) < 2:
        print("Usage: thesis_figures.py file.root [data.root] [output_dir]")
        sys.exit(1)

    root_path = sys.argv[1]
    root_dir  = os.path.dirname(os.path.abspath(root_path))

    data_root_path = os.path.join(root_dir, "fake data.root")
    output_dir = os.path.join(root_dir, "thesis_figures")
    for arg in sys.argv[2:]:
        if arg.endswith(".root"):
            data_root_path = arg
        else:
            output_dir = arg

    f = uproot.open(root_path)

    top_keys = [strip_cycle(k) for k in f.keys(recursive=False)]
    EXCLUDE_TOP = {"Asimov", "AsimovFakeData"}
    detectors = [k for k in top_keys
                 if k not in EXCLUDE_TOP and safe_get(f, k) is not None]
    log.info("Found detectors: %s", detectors)

    asimov_dir = safe_get(f, "Asimov")
    if asimov_dir is None:
        log.error("No Asimov directory in file – nothing to plot.")
        return

    data_hists = load_data_histograms(data_root_path)
    variables = list(VAR_META.keys())

    for det in detectors:
        samples = discover_subsamples(f, det, variables)
        if not samples:
            log.warning("No subsamples found under detector %s – skipping.", det)
            continue
        log.info("Detector %s: %d position(s): %s", det, len(samples), samples)

        det_dir = os.path.join(output_dir, slugify(det))

        for sample in samples:
            sample_label = sample.replace("_", " ").strip()
            position_slug = slugify(sample_label)

            per_var = {}
            for var in variables:
                asimov_node = safe_get(f, f"Asimov/{det}/{var}/{sample}")
                h_asimov = find_hist_in_dir(asimov_node, f"{det}_{sample}_{var}_Asimov")
                if h_asimov is None or is_2d(h_asimov):
                    log.warning("No usable Asimov histogram for %s/%s/%s – skipping.",
                                det, sample, var)
                    continue

                asimov_vals, edges = get_1d(h_asimov)
                bin_widths = np.diff(edges)
                asimov_vals = asimov_vals / bin_widths

                fakedata_node = safe_get(f, f"AsimovFakeData/{det}/{var}/{sample}")
                h_fd = find_hist_in_dir(fakedata_node, f"{det}_{sample}_{var}_AsimovFakeData")
                fakedata_asimov_vals = None
                if h_fd is not None and not is_2d(h_fd):
                    fd_vals, _ = get_1d(h_fd)
                    fakedata_asimov_vals = fd_vals / bin_widths

                toy_node = safe_get(f, f"{det}/{var}/{sample}")
                toy_prefix = f"{det}_{sample}_{var}_posterior_toy_{TOY_VARIANT}_"
                toys = []
                if toy_node is not None:
                    for k in toy_node.keys(recursive=False):
                        ks = strip_cycle(k)
                        if not ks.startswith(toy_prefix):
                            continue
                        try:
                            h = toy_node[k]
                            if is_2d(h):
                                continue
                            v, e = get_1d(h)
                            if not np.allclose(e, edges):
                                continue
                            toys.append(v / np.diff(e))
                        except Exception as exc:
                            log.debug("Could not read toy %s: %s", ks, exc)

                if not toys:
                    log.warning("No '%s' toys found for %s/%s/%s – skipping.",
                                TOY_VARIANT, det, sample, var)
                    continue

                toys_arr = np.vstack(toys)
                per_var[var] = dict(
                    asimov=asimov_vals,
                    mean=np.mean(toys_arr, axis=0),
                    std=np.std(toys_arr, axis=0),
                    edges=edges,
                    toys=toys_arr,
                    fakedata_asimov=fakedata_asimov_vals,
                )

            title_suffix = f"{det} / {sample_label}"

            if "RecoNeutrinoEnergy" in per_var:
                d = per_var["RecoNeutrinoEnergy"]
                plot_reco(
                    det_dir, f"{position_slug}__Reco",
                    asimov=d["asimov"], mean=d["mean"], std=d["std"], edges=d["edges"],
                    title=f"RecoNeutrinoEnergy  —  {title_suffix}",
                    xlabel=VAR_META["RecoNeutrinoEnergy"]["xlabel"],
                    toys=d["toys"], fakedata_asimov=d["fakedata_asimov"],
                    data_hists=data_hists,
                )

            if "TrueNeutrinoEnergy" in per_var:
                d = per_var["TrueNeutrinoEnergy"]
                plot_true(
                    det_dir, f"{position_slug}__True",
                    asimov=d["asimov"], mean=d["mean"], std=d["std"], edges=d["edges"],
                    title=f"TrueNeutrinoEnergy  —  {title_suffix}",
                    xlabel=VAR_META["TrueNeutrinoEnergy"]["xlabel"],
                    toys=d["toys"],
                )

            log.info("Done: %s / %s", det, sample)

    log.info("All thesis figures written under → %s", output_dir)


if __name__ == "__main__":
    main()
