# #!/usr/bin/env python3
# """
# Posterior predictive plots from a ROOT file produced by the neutrino fit.

# Usage
# -----
#     python plot_posterior.py file.root [data.root] [output.pdf]
# """

# import sys
# import logging
# import numpy as np
# import uproot
# import matplotlib
# matplotlib.use("Agg")

# import matplotlib.pyplot as plt
# import matplotlib.colors as mcolors
# from matplotlib.backends.backend_pdf import PdfPages

# log = logging.getLogger(__name__)

# # ── variable metadata ────────────────────────────────────────────────────────

# VAR_META = {
#     "TrueNeutrinoEnergy":  {"xlabel": "True Neutrino Energy (GeV)"},
#     "RecoNeutrinoEnergy":  {"xlabel": "Reconstructed Neutrino Energy (GeV)"},
#     "Enubias":             {"xlabel": "E_{ν} bias (GeV)"},
# }

# # ── helpers ──────────────────────────────────────────────────────────────────

# def safe_get(f, path):
#     """Return f[path] or None, logging any exception."""
#     try:
#         return f[path]
#     except Exception as exc:
#         log.debug("safe_get(%s): %s", path, exc)
#         return None


# def is_2d(h):
#     return len(h.axes) == 2


# def get_1d(h):
#     return h.values(), h.axes[0].edges()


# def get_2d(h):
#     return h.values(), h.axes[0].edges(), h.axes[1].edges()


# def extend(v):
#     """Repeat last element so step plots close correctly."""
#     return np.append(v, v[-1])


# def mask_empty(v):
#     v = v.astype(float)
#     v[v == 0] = np.nan
#     return v


# def strip_cycle(key):
#     """Strip ROOT cycle suffix (';N') from a key string."""
#     return key.split(";")[0]


# # ── data histogram loader ────────────────────────────────────────────────────

# def load_data_histograms(data_root_path):
#     """
#     Open a data ROOT file and return a dict of
#     { histogram_name: (values/bin_width, edges) }.

#     Tries SampleFolder/data_* first; if no SampleFolder is found,
#     falls back to loading all 1-D histograms from the top level.
#     """
#     if data_root_path is None:
#         return {}

#     try:
#         df = uproot.open(data_root_path)
#     except Exception as exc:
#         log.error("Could not open data file %s: %s", data_root_path, exc)
#         return {}

#     # ── try SampleFolder first ────────────────────────────────────────────
#     sample_folder = safe_get(df, "SampleFolder")
#     if sample_folder is not None:
#         data_hists = {}
#         for k in sample_folder.keys(recursive=False):
#             ks = strip_cycle(k)
#             if not ks.startswith("data_"):
#                 continue
#             try:
#                 h = sample_folder[k]
#                 if is_2d(h):
#                     continue
#                 vals, edges = get_1d(h)
#                 data_hists[ks] = (vals / np.diff(edges), edges)
#                 log.info("Loaded data histogram: %s", ks)
#             except Exception as exc:
#                 log.warning("Could not read data histogram %s: %s", ks, exc)
#         return data_hists

#     # ── fall back to top-level histograms ─────────────────────────────────
#     log.info("No SampleFolder found; loading all top-level histograms from %s",
#              data_root_path)
#     data_hists = {}
#     for k in df.keys(recursive=False):
#         ks = strip_cycle(k)
#         try:
#             h = df[k]
#             if is_2d(h):
#                 continue
#             vals, edges = get_1d(h)
#             data_hists[ks] = (vals / np.diff(edges), edges)
#             log.info("Loaded data histogram: %s", ks)
#         except Exception as exc:
#             log.warning("Could not read data histogram %s: %s", ks, exc)

#     return data_hists


# def _friendly_label(name):
#     """Turn a raw histogram key into a readable legend label."""
#     for prefix in ("data_", "hOffAxis", "h"):
#         if name.startswith(prefix):
#             name = name[len(prefix):]
#             break
#     label = name.replace("_", " ").strip()
#     return f"Fake Data: {label}" if label else "Fake Data"


# # ── 1-D posterior predictive ─────────────────────────────────────────────────

# def _toy_envelope(ax, toys, edges, percentiles=(5, 95)):
#     """Draw a shaded percentile envelope (individual throws are not shown)."""
#     lo = np.percentile(toys, percentiles[0], axis=0)
#     hi = np.percentile(toys, percentiles[1], axis=0)
#     ax.fill_between(edges, extend(lo), extend(hi),
#                     step="post", alpha=0.15, color="gray",
#                     label=f"{percentiles[0]}–{percentiles[1]}th pct")


# def plot_1d(pdf, asimov, mean, std, edges, title, xlabel,
#             toys=None, data_hists=None, is_reco=False,
#             fakedata_asimov=None):

#     fig, (ax1, ax2, ax3) = plt.subplots(
#         3, 1, figsize=(9, 8),
#         gridspec_kw={"height_ratios": [3, 1, 1]},
#         sharex=True,
#     )
#     fig.subplots_adjust(hspace=0.08)

#     e = edges
#     m = extend(mean)
#     u = extend(mean + std)
#     l = extend(mean - std)
#     a = extend(asimov)

#     # ── top panel ────────────────────────────────────────────────────────────
#     if toys is not None and len(toys):
#         _toy_envelope(ax1, toys, e)

#     ax1.fill_between(e, l, u, step="post", alpha=0.35, label="±1σ (posterior)")
#     ax1.step(e, m, where="post", label="Fit Result",    linewidth=2, color="steelblue")
#     ax1.step(e, a, where="post", label="MC Prediction", linewidth=2, color="darkorange")

#     if fakedata_asimov is not None:
#         ax1.step(e, extend(fakedata_asimov), where="post",
#                  label="Fake Data Asimov", linewidth=2,
#                  color="firebrick", linestyle="--")

#     # Overlay data histograms only on the RecoNeutrinoEnergy panel
#     if is_reco and data_hists:
#         markers = ["o", "s", "^", "D"]
#         for idx, (name, (dvals, dedges)) in enumerate(data_hists.items()):
#             centres = 0.5 * (dedges[:-1] + dedges[1:])
#             ax1.errorbar(
#                 centres, dvals,
#                 yerr=np.sqrt(np.maximum(dvals * np.diff(dedges), 0)) / np.diff(dedges),
#                 fmt=markers[idx % len(markers)],
#                 markersize=4, linewidth=1,
#                 label=_friendly_label(name),
#                 color=f"C{idx + 3}",
#                 zorder=5,
#             )

#     ax1.set_ylabel("Events / GeV")
#     ax1.set_title(title)
#     ax1.legend(fontsize=7)
#     ax1.grid(alpha=0.3)

#     # ── middle panel: posterior relative uncertainty ──────────────────────
#     rel_unc = np.divide(std, mean, out=np.zeros_like(std), where=mean != 0)

#     ax2.step(e, extend(rel_unc), where="post", color="firebrick",
#              label=r"$\sigma_\mathrm{Fit}\,/\,\mu$", linewidth=1.5)
#     ax2.axhline(0, color="black", linewidth=0.5, linestyle=":", alpha=0.4)
#     ax2.set_ylabel(r"$\sigma_\mathrm{Fit}\,/\,\mu$")
#     ax2.legend(fontsize=7, loc="upper right")
#     ax2.grid(alpha=0.3)

#     # ── bottom panel: bias (MC Prediction / Fit Result) ───────────────────
#     bias    = np.divide(asimov, mean,       out=np.ones_like(asimov), where=mean != 0)
#     bias_up = np.divide(asimov, mean + std, out=np.ones_like(asimov), where=(mean + std) != 0)
#     bias_dn = np.divide(asimov, mean - std, out=np.ones_like(asimov), where=(mean - std) != 0)

#     if fakedata_asimov is not None:
#         fd_bias = np.divide(fakedata_asimov, mean,
#                             out=np.ones_like(fakedata_asimov), where=mean != 0)
#         ax3.step(e, extend(fd_bias), where="post", color="firebrick",
#                 linewidth=1.5, linestyle="--", label="Fake Data Asimov / Fit Result")

#     ax3.fill_between(e, extend(bias_dn), extend(bias_up),
#                      step="post", alpha=0.25, color="darkorange", label="±1σ band")
#     ax3.step(e, extend(bias), where="post", color="darkorange",
#              linewidth=1.5, label="MC Pred. / Fit Result")
#     ax3.axhline(1, color="black", linewidth=0.8, linestyle="--", alpha=0.6)
#     ax3.set_ylabel("Bias")
#     ax3.set_xlabel(xlabel)
#     ax3.legend(fontsize=7, loc="upper right")
#     ax3.grid(alpha=0.3)

#     pdf.savefig(fig)
#     plt.close(fig)


# # ── 2-D plots ────────────────────────────────────────────────────────────────

# def plot_2d(ax, v, xe, ye, title, cmap="viridis", norm=None):
#     mesh = ax.pcolormesh(xe, ye, v.T, shading="auto", cmap=cmap, norm=norm)
#     ax.set_title(title)
#     ax.set_xlabel("Reconstructed Neutrino Energy (GeV)")
#     ax.set_ylabel("Lepton Energy (GeV)")
#     return mesh


# def save_2d(pdf, v, xe, ye, title, cbar_label, cmap="viridis", norm=None):
#     fig, ax = plt.subplots(figsize=(9, 6))
#     mesh = plot_2d(ax, v, xe, ye, title, cmap, norm)
#     cb = fig.colorbar(mesh, ax=ax)
#     cb.set_label(cbar_label)
#     pdf.savefig(fig)
#     plt.close(fig)


# # ── 2-D summary pages ────────────────────────────────────────────────────────

# def plot_2d_summaries(f, detectors, pdf):
#     log.info("Building 2-D summary pages…")

#     for det in detectors:
#         node = safe_get(f, f"{det}/2D/summary")
#         if node is None:
#             continue

#         keys = [strip_cycle(k) for k in node.keys()]
#         samples = {k.replace("_2D_posteriorMean", "")
#                    for k in keys if "_2D_posteriorMean" in k}

#         for s in sorted(samples):
#             try:
#                 h_mean = node[f"{s}_2D_posteriorMean"]
#                 h_post = node[f"{s}_2D_posteriorErr"]
#                 h_stat = node[f"{s}_2D_statErr"]
#             except Exception as exc:
#                 log.warning("Skipping %s/%s: %s", det, s, exc)
#                 continue

#             mean, xe, ye = get_2d(h_mean)
#             post = h_post.values()
#             stat = h_stat.values()

#             mean = mask_empty(mean)
#             post = mask_empty(post)
#             stat = mask_empty(stat)

#             ratio = np.where(stat > 0, post / stat, np.nan)

#             if not np.isfinite(ratio).any():
#                 log.warning("All-NaN ratio for %s/%s – skipping 2-D pages.", det, s)
#                 continue

#             rel_post = np.where(mean > 0, post / mean, np.nan)
#             rel_stat = np.where(mean > 0, stat / mean, np.nan)

#             vmax = np.nanpercentile(ratio, 95)
#             norm = mcolors.TwoSlopeNorm(vmin=0, vcenter=1, vmax=max(2.0, float(vmax)))

#             save_2d(pdf, mean,     xe, ye, f"{s} Fit Result",         "Events")
#             save_2d(pdf, rel_post, xe, ye, f"{s} Posterior σ / N",    "σ / N")
#             save_2d(pdf, rel_stat, xe, ye, f"{s} Stat σ / N",         "σ / N")
#             save_2d(pdf, ratio,    xe, ye, f"{s} Posterior / Stat σ", "σ_post / σ_stat",
#                     cmap="RdBu_r", norm=norm)


# # ── main ─────────────────────────────────────────────────────────────────────

# def main():
#     logging.basicConfig(level=logging.INFO, format="%(levelname)s  %(message)s")

#     if len(sys.argv) < 2:
#         print("Usage: plot_posterior.py file.root [data.root] [output.pdf]")
#         sys.exit(1)

#     root_path = sys.argv[1]

#     # Parse optional arguments: data file (.root) and output (.pdf)
#     data_root_path = "fake data.root"
#     out_path = "plots.pdf"
#     for arg in sys.argv[2:]:
#         if arg.endswith(".root"):
#             data_root_path = arg
#         elif arg.endswith(".pdf"):
#             out_path = arg

#     f = uproot.open(root_path)

#     top_keys = [strip_cycle(k) for k in f.keys(recursive=False)]
#     detectors = [k for k in top_keys if k != "Asimov" and safe_get(f, k) is not None]
#     log.info("Found detectors: %s", detectors)

#     asimov_dir = safe_get(f, "Asimov")

#     # Load data histograms from the (optional) data ROOT file
#     data_hists = load_data_histograms(data_root_path)
#     if data_hists:
#         log.info("Data histograms loaded: %s", list(data_hists.keys()))
#     else:
#         log.info("No data histograms loaded (no data file provided or none found).")

#     variables = list(VAR_META.keys())

#     with PdfPages(out_path) as pdf:

#         # ── 1-D loop ──────────────────────────────────────────────────────
#         for det in detectors:
#             for var in variables:

#                 if asimov_dir is None:
#                     log.warning("No Asimov directory in file – skipping 1-D plots.")
#                     break

#                 asimov_keys = [
#                     strip_cycle(k) for k in asimov_dir.keys(recursive=False)
#                     if strip_cycle(k).startswith(det + "_")
#                     and f"_{var}_Asimov" in strip_cycle(k)
#                 ]

#                 if not asimov_keys:
#                     log.warning("No Asimov histogram for %s/%s – skipping.", det, var)
#                     continue

#                 if len(asimov_keys) > 1:
#                     log.warning(
#                         "Multiple Asimov keys for %s/%s – using first: %s",
#                         det, var, asimov_keys,
#                     )

#                 h_asimov = asimov_dir[asimov_keys[0]]
#                 if is_2d(h_asimov):
#                     log.warning("Asimov hist for %s/%s is 2-D – skipping.", det, var)
#                     continue

#                 asimov_vals, edges = get_1d(h_asimov)
#                 bin_widths = np.diff(edges)
#                 asimov_vals = asimov_vals / bin_widths

#                 # load FakeData Asimov
#                 fakedata_asimov_dir = safe_get(f, "AsimovFakeData")
#                 fakedata_asimov_vals = None
#                 if fakedata_asimov_dir is not None:
#                     fd_keys = [
#                         strip_cycle(k) for k in fakedata_asimov_dir.keys(recursive=False)
#                         if strip_cycle(k).startswith(det + "_")
#                         and f"_{var}_AsimovFakeData" in strip_cycle(k)
#                     ]
#                     if fd_keys:
#                         h_fd = fakedata_asimov_dir[fd_keys[0]]
#                         if not is_2d(h_fd):
#                             fd_vals, _ = get_1d(h_fd)
#                             fakedata_asimov_vals = fd_vals / bin_widths

#                 # ── collect toys ──────────────────────────────────────────
#                 toy_marker = f"_{var}_posterior_toy_"
#                 toys = []

#                 def _try_add_toy(h, key):
#                     try:
#                         if is_2d(h):
#                             return
#                         v, e = get_1d(h)
#                         if not np.allclose(e, edges):
#                             log.debug("Edge mismatch for toy %s – skipped.", key)
#                             return
#                         toys.append(v / np.diff(e))
#                     except Exception as exc:
#                         log.debug("Could not read toy %s: %s", key, exc)

#                 var_subdir = safe_get(f, f"{det}/{var}")
#                 if var_subdir is not None:
#                     for k in var_subdir.keys(recursive=False):
#                         ks = strip_cycle(k)
#                         if toy_marker in ks:
#                             _try_add_toy(var_subdir[k], ks)

#                 det_dir = safe_get(f, det)
#                 if det_dir is not None:
#                     for k in det_dir.keys(recursive=False):
#                         ks = strip_cycle(k)
#                         if toy_marker in ks:
#                             _try_add_toy(det_dir[k], ks)

#                 if not toys:
#                     log.warning("No toys found for %s/%s.", det, var)
#                     continue

#                 toys_arr = np.vstack(toys)
#                 mean = np.mean(toys_arr, axis=0)
#                 std  = np.std(toys_arr,  axis=0)

#                 xlabel  = VAR_META[var]["xlabel"]
#                 is_reco = (var == "RecoNeutrinoEnergy")

#                 plot_1d(pdf, asimov_vals, mean, std, edges,
#                         title=f"{var}",
#                         xlabel=xlabel,
#                         toys=toys_arr,
#                         data_hists=data_hists if is_reco else None,
#                         is_reco=is_reco)
#                 plot_1d(pdf, asimov_vals, mean, std, edges,
#                         title=f"{var}",
#                         xlabel=xlabel,
#                         toys=toys_arr,
#                         data_hists=data_hists if is_reco else None,
#                         is_reco=is_reco,
#                         fakedata_asimov=fakedata_asimov_vals)

#                 log.info("1-D done: %s / %s  (%d toys)", det, var, len(toys))

#         # ── 2-D summaries ─────────────────────────────────────────────────
#         plot_2d_summaries(f, detectors, pdf)

#     log.info("Saved → %s", out_path)


# if __name__ == "__main__":
#     main()

#!/usr/bin/env python3
"""
Posterior predictive plots from a ROOT file produced by the neutrino fit.

Usage
-----
    python plot_posterior.py file.root [data.root] [output.pdf]
"""
#!/usr/bin/env python3
"""
Posterior predictive plots from a ROOT file produced by the neutrino fit.

Usage
-----
    python plot_posterior.py file.root [data.root] [output.pdf]
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


def plot_1d(pdf, asimov, mean, std, edges, title, xlabel,
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
    ax1.step(e, m, where="post", label="Fit Result",
             linewidth=1.8, color=_C["blue"])

    if not no_mc_prediction:
        ax1.step(e, a, where="post", label="CV Prediction",
                 linewidth=1.8, color=_C["orange"])

    if fakedata_asimov is not None and not no_fakedata_asimov_line:
        ax1.step(e, extend(fakedata_asimov), where="post",
                 label="Fake Data", linewidth=1.8,
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
            f"(Fit Result \u2212 Fake Data ) / Fake Data   "
            f"[mean: {ratio_means:+.3f}]"
        )
        ax3.step(e, extend(fd_bias), where="post", color=_C["red"],
                 linewidth=1.8, linestyle=(0, (5, 2)), label="Fake Data  / Fit Result")
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


# ── 2-D summary pages ────────────────────────────────────────────────────────

def plot_2d_summaries(f, detectors, pdf):
    log.info("Building 2-D summary pages…")

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

            save_2d(pdf, mean,     xe, ye, f"{s} Fit Result",         "Events")
            save_2d(pdf, rel_post, xe, ye, f"{s} Posterior σ / N",    "σ / N")
            save_2d(pdf, rel_stat, xe, ye, f"{s} Stat σ / N",         "σ / N")
            save_2d(pdf, ratio,    xe, ye, f"{s} Posterior / Stat σ", "σ_post / σ_stat",
                    cmap="RdBu_r", norm=norm)


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    logging.basicConfig(level=logging.INFO, format="%(levelname)s  %(message)s")

    if len(sys.argv) < 2:
        print("Usage: plot_posterior.py file.root [data.root] [output.pdf]")
        sys.exit(1)

    root_path = sys.argv[1]
    root_dir  = os.path.dirname(os.path.abspath(root_path))

    # Parse optional arguments: data file (.root) and output (.pdf)
    data_root_path = os.path.join(root_dir, "fake data.root")
    out_path = os.path.join(root_dir, "plots.pdf")
    for arg in sys.argv[2:]:
        if arg.endswith(".root"):
            data_root_path = arg
        elif arg.endswith(".pdf"):
            out_path = arg

    f = uproot.open(root_path)

    top_keys = [strip_cycle(k) for k in f.keys(recursive=False)]
    detectors = [k for k in top_keys if k != "Asimov" and safe_get(f, k) is not None]
    log.info("Found detectors: %s", detectors)

    asimov_dir = safe_get(f, "Asimov")

    # Load data histograms from the (optional) data ROOT file
    data_hists = load_data_histograms(data_root_path)
    if data_hists:
        log.info("Data histograms loaded: %s", list(data_hists.keys()))
    else:
        log.info("No data histograms loaded (no data file provided or none found).")

    variables = list(VAR_META.keys())

    with PdfPages(out_path) as pdf:

        # ── 1-D loop ──────────────────────────────────────────────────────
        for det in detectors:
            for var in variables:

                if asimov_dir is None:
                    log.warning("No Asimov directory in file – skipping 1-D plots.")
                    break

                asimov_keys = [
                    strip_cycle(k) for k in asimov_dir.keys(recursive=False)
                    if strip_cycle(k).startswith(det + "_")
                    and f"_{var}_Asimov" in strip_cycle(k)
                ]

                if not asimov_keys:
                    log.warning("No Asimov histogram for %s/%s – skipping.", det, var)
                    continue

                if len(asimov_keys) > 1:
                    log.warning(
                        "Multiple Asimov keys for %s/%s – using first: %s",
                        det, var, asimov_keys,
                    )

                h_asimov = asimov_dir[asimov_keys[0]]
                if is_2d(h_asimov):
                    log.warning("Asimov hist for %s/%s is 2-D – skipping.", det, var)
                    continue

                asimov_vals, edges = get_1d(h_asimov)
                bin_widths = np.diff(edges)
                asimov_vals = asimov_vals / bin_widths

                # load FakeData Asimov
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

                # ── collect toys ──────────────────────────────────────────
                toy_marker = f"_{var}_posterior_toy_"
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
                    log.warning("No toys found for %s/%s.", det, var)
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
                    title=f"{var}", xlabel=xlabel,
                    toys=toys_arr,
                    data_hists=data_hists if is_reco else None,
                    is_reco=is_reco,
                    fakedata_asimov=fakedata_asimov_vals,
                )

                # ── standard plot ─────────────────────────────────────────
                plot_1d(pdf, **common_kwargs)

                # ── standard plot + fakedata asimov line ──────────────────
                plot_1d(pdf, **common_kwargs)

                # ── RecoNeutrinoEnergy variant ────────────────────────────
                # No MC Prediction, blue sigma line, bias title with mean ratio
                if is_reco:
                    plot_1d(pdf, **common_kwargs,
                            no_mc_prediction=True,
                            sigma_blue=True)

                # ── TrueNeutrinoEnergy variant ────────────────────────────
                # No Fake Data Asimov line, blue sigma line, percentile band on sigma panel
                if is_true:
                    plot_1d(pdf, **common_kwargs,
                            no_fakedata_asimov_line=True,
                            sigma_blue=True,
                            sigma_toy_envelope=True)

                # ── Enubias variant ───────────────────────────────────────
                # x-axis lower limit = -0.2
                if is_bias:
                    plot_1d(pdf, **common_kwargs,
                            xlim_min=-0.2)

                log.info("1-D done: %s / %s  (%d toys)", det, var, len(toys))

        # ── 2-D summaries ─────────────────────────────────────────────────
        plot_2d_summaries(f, detectors, pdf)

    log.info("Saved → %s", out_path)


if __name__ == "__main__":
    main()