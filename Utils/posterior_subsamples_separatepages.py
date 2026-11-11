#!/usr/bin/env python3
"""
Posterior predictive plots from a ROOT file produced by the neutrino fit.

This is a variant of plot_posterior.py that writes each plot out as its own
PDF file (organised into a folder tree by detector/subsample) instead of
appending every plot as a page of one big combined PDF. Useful if you want
to grab, review, or share individual plots without paging through the whole
document.

Handles files where each detector (ND/FD) contains multiple subsamples
(e.g. the 18 OffAxisND positions), written by MakeSpectaVariations as:

    Asimov/{ND_or_FD}/{var}/{cleanTitle}/{ND_or_FD}_{cleanTitle}_{var}_Asimov
    AsimovFakeData/{ND_or_FD}/{var}/{cleanTitle}/..._AsimovFakeData
    {ND_or_FD}/{var}/{cleanTitle}/..._posterior_toy_{variant}_{NNN}

where {variant} is one of "xsecOnly", "fluxOnly", "combined".

Output layout
-------------
    <output_dir>/<detector>/<subsample>/<variable>__<variant>.pdf
    <output_dir>/<detector>/2Dsummary/<subsample>__<plot_type>.pdf

Usage
-----
    python plot_posterior_pages.py file.root [data.root] [output_dir]

    - file.root    : required, the fit-output ROOT file
    - data.root    : optional, defaults to "<file.root's dir>/fake data.root"
    - output_dir   : optional, defaults to "<file.root's dir>/plots/"
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
import matplotlib.colors as mcolors

log = logging.getLogger(__name__)

# ── configuration ───────────────────────────────────────────────────────────

# Which toy parameter-variation variant to use for the posterior mean/std and
# envelope. Must match one of the tags written by RunToyVariant in the C++
# generator: "xsecOnly", "fluxOnly", "combined".
TOY_VARIANT = "combined"

# ── thesis-quality style ──────────────────────────────────────────────────────
# Okabe-Ito colourblind-safe palette
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
    # fonts
    "font.family":        "serif",
    "font.serif":         ["STIXTwo Text", "DejaVu Serif", "Times New Roman"],
    "mathtext.fontset":   "stix",
    "font.size":          11,
    "axes.titlesize":     11,
    "axes.labelsize":     11,
    "xtick.labelsize":    10,
    "ytick.labelsize":    10,
    "legend.fontsize":    9,
    # embedded fonts for thesis PDF submission
    "pdf.fonttype":       42,
    "ps.fonttype":        42,
    # figure
    "figure.dpi":         150,
    # lines
    "lines.linewidth":    1.8,
    "patch.linewidth":    1.0,
    # ticks
    "xtick.direction":    "in",
    "ytick.direction":    "in",
    "xtick.top":          True,
    "ytick.right":        True,
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    # grid
    "axes.grid":          True,
    "grid.alpha":         0.3,
    "grid.linewidth":     0.6,
    # layout
    "axes.linewidth":     0.8,
}

plt.rcParams.update(THESIS_RC)

# ── variable metadata ────────────────────────────────────────────────────────

VAR_META = {
    "TrueNeutrinoEnergy":  {"xlabel": "True Neutrino Energy (GeV)"},
    "RecoNeutrinoEnergy":  {"xlabel": "Reconstructed Neutrino Energy (GeV)"},
    "Enubias":             {"xlabel": "E_{ν} bias (GeV)"},
}

# ── helpers ──────────────────────────────────────────────────────────────────

def slugify(text):
    """Make a string safe to use as a filename/directory component."""
    text = re.sub(r"[^\w\-]+", "_", text)
    return re.sub(r"_+", "_", text).strip("_")


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


def find_hist_in_dir(node, prefix):
    """
    Return the histogram object inside `node` whose (cycle-stripped) key
    starts with `prefix`, or None if not found. Warns on ambiguous matches.
    """
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
    """
    Collect the set of subsample ('cleanTitle') names available for a
    detector by scanning the toy-histogram directories {det}/{var}/*.
    """
    samples = set()
    for var in variables:
        node = safe_get(f, f"{det}/{var}")
        if node is None:
            continue
        for k in node.keys(recursive=False):
            samples.add(strip_cycle(k))
    return sorted(samples)


# ── data histogram loader ────────────────────────────────────────────────────

def load_data_histograms(data_root_path):
    """
    Open a data ROOT file and return a dict of
    { histogram_name: (values/bin_width, edges) }.

    Tries SampleFolder/data_* first; if no SampleFolder is found,
    falls back to loading all 1-D histograms from the top level.
    """
    if data_root_path is None:
        return {}

    try:
        df = uproot.open(data_root_path)
    except Exception as exc:
        log.error("Could not open data file %s: %s", data_root_path, exc)
        return {}

    # ── try SampleFolder first ────────────────────────────────────────────
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

    # ── fall back to top-level histograms ─────────────────────────────────
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
    """Turn a raw histogram key into a readable legend label."""
    for prefix in ("data_", "hOffAxis", "h"):
        if name.startswith(prefix):
            name = name[len(prefix):]
            break
    label = name.replace("_", " ").strip()
    return f"Fake Data: {label}" if label else "Fake Data"


# ── 1-D posterior predictive ─────────────────────────────────────────────────

def _toy_envelope(ax, toys, edges, percentiles=(5, 95)):
    """Draw a shaded percentile envelope (individual throws are not shown)."""
    lo = np.percentile(toys, percentiles[0], axis=0)
    hi = np.percentile(toys, percentiles[1], axis=0)
    ax.fill_between(edges, extend(lo), extend(hi),
                    step="post", alpha=0.12, color=_C["black"],
                    label=f"{percentiles[0]}–{percentiles[1]}th pct")


def plot_1d(output_dir, filename, asimov, mean, std, edges, title, xlabel,
            toys=None, data_hists=None, is_reco=False,
            fakedata_asimov=None,
            # variant flags
            no_mc_prediction=False,
            no_fakedata_asimov_line=False,
            sigma_blue=False,
            sigma_toy_envelope=False,
            xlim_min=None):

    # thesis print width: 6.5" wide, proportional height
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

    # ── top panel ────────────────────────────────────────────────────────────
    if toys is not None and len(toys):
        _toy_envelope(ax1, toys, e)

    ax1.fill_between(e, l, u, step="post", alpha=0.25,
                     color=_C["blue"], label=r"$\pm1\sigma$ (posterior)")
    ax1.step(e, m, where="post", label="Best Fit",
             linewidth=1.8, color=_C["blue"])

    if not no_mc_prediction:
        ax1.step(e, a, where="post", label="Nominal",
                 linewidth=1.8, color=_C["orange"])

    if fakedata_asimov is not None and not no_fakedata_asimov_line:
        ax1.step(e, extend(fakedata_asimov), where="post",
                 label="Fake Data ", linewidth=1.8,
                 color=_C["red"], linestyle=(0, (5, 2)))

    # Overlay data histograms only on the RecoNeutrinoEnergy panel
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
    # remove x tick labels on top two panels (shared axis handles it)
    plt.setp(ax1.get_xticklabels(), visible=False)

    # ── middle panel: posterior relative uncertainty ──────────────────────
    rel_unc = np.divide(std, mean, out=np.zeros_like(std), where=mean != 0)
    sigma_color = _C["blue"] if sigma_blue else _C["red"]

    # Optional toy envelope on the sigma panel
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

    # ── bottom panel: MC / Fit ────────────────────────────────────────────
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
            f"(Best Fit \u2212 Fake Data ) / Fake Data   "
            f"[mean: {ratio_means:+.3f}]"
        )
        ax3.step(e, extend(fd_bias), where="post", color=_C["red"],
                 linewidth=1.8, linestyle=(0, (5, 2)), label="Fake Data  / Best Fit")
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

    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, filename + ".pdf")
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    log.info("Wrote %s", out_path)


# ── 2-D plots ────────────────────────────────────────────────────────────────

def plot_2d(ax, v, xe, ye, title, cmap="viridis", norm=None):
    mesh = ax.pcolormesh(xe, ye, v.T, shading="auto", cmap=cmap, norm=norm)
    ax.set_title(title)
    ax.set_xlabel("Reconstructed Neutrino Energy (GeV)")
    ax.set_ylabel("Lepton Energy (GeV)")
    return mesh


def save_2d(output_dir, filename, v, xe, ye, title, cbar_label, cmap="viridis", norm=None):
    fig, ax = plt.subplots(figsize=(6.5, 5.0))
    mesh = plot_2d(ax, v, xe, ye, title, cmap, norm)
    cb = fig.colorbar(mesh, ax=ax)
    cb.set_label(cbar_label)

    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, filename + ".pdf")
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    log.info("Wrote %s", out_path)


# ── 2-D summary pages ────────────────────────────────────────────────────────
# NOTE: the C++ generator shown alongside this script only writes the 1-D
# TrueNeutrinoEnergy / RecoNeutrinoEnergy / Enubias histograms, so this
# function is left untouched apart from the output-path handling. If your
# 2D/summary histograms also gained a per-subsample directory level, this
# will need the same treatment as discover_subsamples()/find_hist_in_dir()
# above.

def plot_2d_summaries(f, detectors, output_dir):
    log.info("Building 2-D summary pages…")

    for det in detectors:
        node = safe_get(f, f"{det}/2D/summary")
        if node is None:
            continue

        det_dir = os.path.join(output_dir, slugify(det), "2Dsummary")

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

            s_slug = slugify(s)
            save_2d(det_dir, f"{s_slug}__fitResult",   mean,     xe, ye, f"{s} Best Fit",         "Events")
            save_2d(det_dir, f"{s_slug}__postSigmaN",  rel_post, xe, ye, f"{s} Posterior σ / N",    "σ / N")
            save_2d(det_dir, f"{s_slug}__statSigmaN",  rel_stat, xe, ye, f"{s} Stat σ / N",         "σ / N")
            save_2d(det_dir, f"{s_slug}__postOverStat", ratio,   xe, ye, f"{s} Posterior / Stat σ", "σ_post / σ_stat",
                    cmap="RdBu_r", norm=norm)


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    logging.basicConfig(level=logging.INFO, format="%(levelname)s  %(message)s")

    if len(sys.argv) < 2:
        print("Usage: plot_posterior_pages.py file.root [data.root] [output_dir]")
        sys.exit(1)

    root_path = sys.argv[1]
    root_dir  = os.path.dirname(os.path.abspath(root_path))

    # Parse optional arguments: data file (.root) and output directory
    data_root_path = os.path.join(root_dir, "fake data.root")
    output_dir = os.path.join(root_dir, "plots")
    for arg in sys.argv[2:]:
        if arg.endswith(".root"):
            data_root_path = arg
        else:
            output_dir = arg

    f = uproot.open(root_path)

    top_keys = [strip_cycle(k) for k in f.keys(recursive=False)]
    # "Asimov" and "AsimovFakeData" are auxiliary top-level dirs, not detectors
    EXCLUDE_TOP = {"Asimov", "AsimovFakeData"}
    detectors = [k for k in top_keys
                 if k not in EXCLUDE_TOP and safe_get(f, k) is not None]
    log.info("Found detectors: %s", detectors)

    asimov_dir = safe_get(f, "Asimov")

    # Load data histograms from the (optional) data ROOT file
    data_hists = load_data_histograms(data_root_path)
    if data_hists:
        log.info("Data histograms loaded: %s", list(data_hists.keys()))
    else:
        log.info("No data histograms loaded (no data file provided or none found).")

    variables = list(VAR_META.keys())

    # ── 1-D loop: detector → subsample → variable ───────────────────────
    for det in detectors:

        if asimov_dir is None:
            log.warning("No Asimov directory in file – skipping 1-D plots.")
            break

        samples = discover_subsamples(f, det, variables)
        if not samples:
            log.warning("No subsamples found under detector %s – skipping.", det)
            continue
        log.info("Detector %s: %d subsample(s): %s", det, len(samples), samples)

        for sample in samples:
            sample_label = sample.replace("_", " ").strip()
            sample_dir = os.path.join(output_dir, slugify(det), slugify(sample))

            for var in variables:

                # ── Asimov (nominal MC prediction) ─────────────────
                asimov_node = safe_get(f, f"Asimov/{det}/{var}/{sample}")
                h_asimov = find_hist_in_dir(
                    asimov_node, f"{det}_{sample}_{var}_Asimov")
                if h_asimov is None:
                    log.warning("No Asimov histogram for %s/%s/%s – skipping.",
                                det, sample, var)
                    continue
                if is_2d(h_asimov):
                    log.warning("Asimov hist for %s/%s/%s is 2-D – skipping.",
                                det, sample, var)
                    continue

                asimov_vals, edges = get_1d(h_asimov)
                bin_widths = np.diff(edges)
                asimov_vals = asimov_vals / bin_widths

                # ── FakeData Asimov ─────────────────────────────────
                fakedata_node = safe_get(f, f"AsimovFakeData/{det}/{var}/{sample}")
                h_fd = find_hist_in_dir(
                    fakedata_node, f"{det}_{sample}_{var}_AsimovFakeData")
                fakedata_asimov_vals = None
                if h_fd is not None and not is_2d(h_fd):
                    fd_vals, _ = get_1d(h_fd)
                    fakedata_asimov_vals = fd_vals / bin_widths

                # ── collect toys (filtered to TOY_VARIANT) ──────────
                toy_node = safe_get(f, f"{det}/{var}/{sample}")
                toy_prefix = f"{det}_{sample}_{var}_posterior_toy_{TOY_VARIANT}_"
                toys = []

                def _try_add_toy(h, key):
                    try:
                        if is_2d(h):
                            return
                        v, e = get_1d(h)
                        if not np.allclose(e, edges):
                            log.debug("Edge mismatch for toy %s – skipped.", key)
                            return
                        toys.append(v / np.diff(e))
                    except Exception as exc:
                        log.debug("Could not read toy %s: %s", key, exc)

                if toy_node is not None:
                    for k in toy_node.keys(recursive=False):
                        ks = strip_cycle(k)
                        if ks.startswith(toy_prefix):
                            _try_add_toy(toy_node[k], ks)

                if not toys:
                    log.warning("No '%s' toys found for %s/%s/%s.",
                                TOY_VARIANT, det, sample, var)
                    continue

                toys_arr = np.vstack(toys)
                mean = np.mean(toys_arr, axis=0)
                std  = np.std(toys_arr,  axis=0)

                xlabel  = VAR_META[var]["xlabel"]
                is_reco = (var == "RecoNeutrinoEnergy")
                is_true = (var == "TrueNeutrinoEnergy")
                is_bias = (var == "Enubias")

                title = f"{var}  —  {det} / {sample_label}"

                common_kwargs = dict(
                    asimov=asimov_vals, mean=mean, std=std, edges=edges,
                    title=title, xlabel=xlabel,
                    toys=toys_arr,
                    data_hists=data_hists if is_reco else None,
                    is_reco=is_reco,
                    fakedata_asimov=fakedata_asimov_vals,
                )

                # ── standard plot ─────────────────────────────────
                plot_1d(sample_dir, f"{var}__standard", **common_kwargs)

                # ── RecoNeutrinoEnergy variants ────────────────────
                if is_reco:
                    plot_1d(sample_dir, f"{var}__fitOnly", **common_kwargs,
                            no_mc_prediction=True,
                            sigma_blue=True)

                    # "clean comparison": just Best Fit vs Nominal,
                    # no data points, no Fake Data  line, no toy
                    # envelope shading.
                    clean_kwargs = dict(common_kwargs)
                    clean_kwargs.update(
                        data_hists=None,
                        fakedata_asimov=None,
                        toys=None,
                        title=f"{title} (Best Fit vs Nominal)",
                    )
                    plot_1d(sample_dir, f"{var}__fitVsCV", **clean_kwargs)

                # ── TrueNeutrinoEnergy variant ─────────────────────
                if is_true:
                    plot_1d(sample_dir, f"{var}__toyEnvelope", **common_kwargs,
                            no_fakedata_asimov_line=True,
                            sigma_blue=True,
                            sigma_toy_envelope=True)

                # ── Enubias variant ─────────────────────────────────
                if is_bias:
                    plot_1d(sample_dir, f"{var}__trimmed", **common_kwargs,
                            xlim_min=-0.2)

                log.info("1-D done: %s / %s / %s  (%d '%s' toys)",
                         det, sample, var, len(toys), TOY_VARIANT)

    # ── 2-D summaries ─────────────────────────────────────────────────────
    plot_2d_summaries(f, detectors, output_dir)

    log.info("All plots written under → %s", output_dir)


if __name__ == "__main__":
    main()