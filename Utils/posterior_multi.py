#!/usr/bin/env python3
"""
Posterior predictive plots from TWO ROOT files (e.g. "On Axis" and
"On+Off Axis") produced by the neutrino fit. Both files must share the
same internal structure as the original single-file version of this
script (Asimov/, AsimovFakeData/, <Detector>/<Var>/..._posterior_toy_,
<Detector>/2D/summary/...).

Usage
-----
    python plot_posterior_multi.py onaxis.root onoffaxis.root [data.root] [output.pdf]

Optional labels can be set via the FILE_LABELS list below if you want
something other than "On Axis" / "On+Off Axis" in plot titles/logs.
"""

import os
import sys
import logging
import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.backends.backend_pdf import PdfPages

log = logging.getLogger(__name__)

# ── labels for the two input files, in the order given on the command line ──
FILE_LABELS = ["On Axis", "On+Off Axis"]

# ── thesis-quality style ──────────────────────────────────────────────────────
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
    "TrueNeutrinoEnergy":  {"xlabel": "True Neutrino Energy (GeV)"},
    "RecoNeutrinoEnergy":  {"xlabel": "Reconstructed Neutrino Energy (GeV)"},
    "Enubias":             {"xlabel": "E_{ν} bias (GeV)"},
}

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


def get_2d(h):
    return h.values(), h.axes[0].edges(), h.axes[1].edges()


def extend(v):
    return np.append(v, v[-1])


def mask_empty(v):
    v = v.astype(float)
    v[v == 0] = np.nan
    return v


def strip_cycle(key):
    return key.split(";")[0]


# ── data histogram loader ────────────────────────────────────────────────────

def load_data_histograms(data_root_path):
    if data_root_path is None:
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
                log.info("Loaded data histogram: %s", ks)
            except Exception as exc:
                log.warning("Could not read data histogram %s: %s", ks, exc)
        return data_hists

    log.info("No SampleFolder found; loading all top-level histograms from %s",
             data_root_path)
    data_hists = {}
    for k in df.keys(recursive=False):
        ks = strip_cycle(k)
        try:
            h = df[k]
            if is_2d(h):
                continue
            vals, edges = get_1d(h)
            data_hists[ks] = (vals / np.diff(edges), edges)
            log.info("Loaded data histogram: %s", ks)
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


# ── 1-D posterior predictive ─────────────────────────────────────────────────

def _toy_envelope(ax, toys, edges, percentiles=(5, 95)):
    lo = np.percentile(toys, percentiles[0], axis=0)
    hi = np.percentile(toys, percentiles[1], axis=0)
    ax.fill_between(edges, extend(lo), extend(hi),
                    step="post", alpha=0.12, color=_C["black"],
                    label=f"{percentiles[0]}–{percentiles[1]}th pct")


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
        data_colors  = [_C["green"], _C["purple"], _C["sky"], _C["yellow"]]
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
                         label="5–95th pct")

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
        fd_mean_val  = np.nanmean(fakedata_asimov)
        ratio_means  = (fit_mean_val - fd_mean_val) / fd_mean_val if fd_mean_val != 0 else np.nan
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
             linewidth=1.8, label="CV/Fit")
    ax3.axhline(1, color="black", linewidth=0.8, linestyle="--", alpha=0.7)
    ax3.set_ylabel("MC/Fit")
    ax3.set_xlabel(xlabel)
    ax3.legend(loc="upper right", framealpha=0.9, edgecolor="0.7")

    if xlim_min is not None:
        ax3.set_xlim(left=xlim_min)

    fig.align_ylabels([ax1, ax2, ax3])
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


# ── 2-D plots ────────────────────────────────────────────────────────────────

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


def plot_2d_summaries(f, detectors, pdf, file_label=""):
    log.info("[%s] Building 2-D summary pages…", file_label)

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
                log.warning("All-NaN ratio for %s/%s – skipping 2-D pages.", det, s)
                continue

            rel_post = np.where(mean > 0, post / mean, np.nan)
            rel_stat = np.where(mean > 0, stat / mean, np.nan)

            vmax = np.nanpercentile(ratio, 95)
            norm = mcolors.TwoSlopeNorm(vmin=0, vcenter=1, vmax=max(2.0, float(vmax)))

            prefix = f"[{file_label}] " if file_label else ""
            save_2d(pdf, mean,     xe, ye, f"{prefix}{s} Fit Result",         "Events")
            save_2d(pdf, rel_post, xe, ye, f"{prefix}{s} Posterior σ / N",    "σ / N")
            save_2d(pdf, rel_stat, xe, ye, f"{prefix}{s} Stat σ / N",         "σ / N")
            save_2d(pdf, ratio,    xe, ye, f"{prefix}{s} Posterior / Stat σ", "σ_post / σ_stat",
                    cmap="RdBu_r", norm=norm)


# ── per-file 1-D + 2-D loop ─────────────────────────────────────────────────

def process_file(f, data_hists, pdf, file_label=""):
    """Run the full 1-D + 2-D plotting loop against an already-open uproot
    file, tagging every plot title with file_label so pages from different
    inputs (e.g. On Axis vs On+Off Axis) stay distinguishable in one PDF."""

    top_keys = [strip_cycle(k) for k in f.keys(recursive=False)]
    detectors = [k for k in top_keys if k != "Asimov" and safe_get(f, k) is not None]
    log.info("[%s] Found detectors: %s", file_label, detectors)

    asimov_dir = safe_get(f, "Asimov")
    variables = list(VAR_META.keys())
    prefix = f"[{file_label}] " if file_label else ""

    for det in detectors:
        for var in variables:

            if asimov_dir is None:
                log.warning("[%s] No Asimov directory in file – skipping 1-D plots.", file_label)
                break

            asimov_keys = [
                strip_cycle(k) for k in asimov_dir.keys(recursive=False)
                if strip_cycle(k).startswith(det + "_")
                and f"_{var}_Asimov" in strip_cycle(k)
            ]

            if not asimov_keys:
                log.warning("[%s] No Asimov histogram for %s/%s – skipping.", file_label, det, var)
                continue

            if len(asimov_keys) > 1:
                log.warning(
                    "[%s] Multiple Asimov keys for %s/%s – using first: %s",
                    file_label, det, var, asimov_keys,
                )

            h_asimov = asimov_dir[asimov_keys[0]]
            if is_2d(h_asimov):
                log.warning("[%s] Asimov hist for %s/%s is 2-D – skipping.", file_label, det, var)
                continue

            asimov_vals, edges = get_1d(h_asimov)
            bin_widths = np.diff(edges)
            asimov_vals = asimov_vals / bin_widths

            fakedata_asimov_dir = safe_get(f, "AsimovFakeData")
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
            toys = []

            def _try_add_toy(h, key):
                try:
                    if is_2d(h):
                        return
                    v, e = get_1d(h)
                    if not np.allclose(e, edges):
                        log.debug("[%s] Edge mismatch for toy %s – skipped.", file_label, key)
                        return
                    toys.append(v / np.diff(e))
                except Exception as exc:
                    log.debug("[%s] Could not read toy %s: %s", file_label, key, exc)

            var_subdir = safe_get(f, f"{det}/{var}")
            if var_subdir is not None:
                for k in var_subdir.keys(recursive=False):
                    ks = strip_cycle(k)
                    if toy_marker in ks:
                        _try_add_toy(var_subdir[k], ks)

            det_dir = safe_get(f, det)
            if det_dir is not None:
                for k in det_dir.keys(recursive=False):
                    ks = strip_cycle(k)
                    if toy_marker in ks:
                        _try_add_toy(det_dir[k], ks)

            if not toys:
                log.warning("[%s] No toys found for %s/%s.", file_label, det, var)
                continue

            toys_arr = np.vstack(toys)
            mean = np.mean(toys_arr, axis=0)
            std  = np.std(toys_arr,  axis=0)

            xlabel  = VAR_META[var]["xlabel"]
            is_reco = (var == "RecoNeutrinoEnergy")
            is_true = (var == "TrueNeutrinoEnergy")
            is_bias = (var == "Enubias")

            common_kwargs = dict(
                asimov=asimov_vals, mean=mean, std=std, edges=edges,
                title=f"{prefix}{var}", xlabel=xlabel,
                toys=toys_arr,
                data_hists=data_hists if is_reco else None,
                is_reco=is_reco,
                fakedata_asimov=fakedata_asimov_vals,
            )

            plot_1d(pdf, **common_kwargs)
            plot_1d(pdf, **common_kwargs)

            if is_reco:
                plot_1d(pdf, **common_kwargs,
                        no_mc_prediction=True,
                        sigma_blue=True)

                clean_kwargs = dict(common_kwargs)
                clean_kwargs.update(
                    data_hists=None,
                    fakedata_asimov=None,
                    toys=None,
                    title=f"{prefix}{var} (Fit Result vs CV Prediction)",
                )
                plot_1d(pdf, **clean_kwargs)

            if is_true:
                plot_1d(pdf, **common_kwargs,
                        no_fakedata_asimov_line=True,
                        sigma_blue=True,
                        sigma_toy_envelope=True)

            if is_bias:
                plot_1d(pdf, **common_kwargs,
                        xlim_min=-0.2)

            log.info("[%s] 1-D done: %s / %s  (%d toys)", file_label, det, var, len(toys))

    plot_2d_summaries(f, detectors, pdf, file_label=file_label)


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    logging.basicConfig(level=logging.INFO, format="%(levelname)s  %(message)s")

    if len(sys.argv) < 3:
        print("Usage: plot_posterior_multi.py file1.root file2.root [data.root] [output.pdf]")
        sys.exit(1)

    root_path_1 = sys.argv[1]
    root_path_2 = sys.argv[2]
    root_dir = os.path.dirname(os.path.abspath(root_path_1))

    data_root_path = os.path.join(root_dir, "fake data.root")
    out_path = os.path.join(root_dir, "plots.pdf")
    for arg in sys.argv[3:]:
        if arg.endswith(".root"):
            data_root_path = arg
        elif arg.endswith(".pdf"):
            out_path = arg

    root_paths = [root_path_1, root_path_2]
    labels = FILE_LABELS if len(FILE_LABELS) == len(root_paths) else \
        [os.path.basename(p) for p in root_paths]

    # Same fake-data/data file overlaid on both inputs' Reco panels.
    # If each file needs its own data overlay instead, load per-file inside
    # the loop below using a matching data path convention.
    data_hists = load_data_histograms(data_root_path)
    if data_hists:
        log.info("Data histograms loaded: %s", list(data_hists.keys()))
    else:
        log.info("No data histograms loaded (no data file provided or none found).")

    with PdfPages(out_path) as pdf:
        for label, path in zip(labels, root_paths):
            log.info("Opening %s (%s)", path, label)
            f = uproot.open(path)
            process_file(f, data_hists, pdf, file_label=label)

    log.info("Saved → %s", out_path)


if __name__ == "__main__":
    main()