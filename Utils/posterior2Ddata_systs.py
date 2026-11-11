#!/usr/bin/env python3
"""
Posterior predictive plots from a ROOT file produced by the neutrino fit.

Usage
-----
    python posterior2Ddata_systs.py file.root [data.root] [output.pdf] [--fakedata=True|False]

Notes
-----
Histograms in the fit output are nested as:
    <Asimov|AsimovFakeData|det>/<var>/<sample>/<det>_<sample>_<var>_<Asimov|AsimovFakeData|posterior_toy_...>
e.g.
    Asimov/ND/TrueNeutrinoEnergy/OffAxis0m_numuCC_numode/ND_OffAxis0m_numuCC_numode_TrueNeutrinoEnergy_Asimov
    ND/TrueNeutrinoEnergy/OffAxis0m_numuCC_numode/ND_OffAxis0m_numuCC_numode_TrueNeutrinoEnergy_posterior_toy_combined_996

All histogram lookups below use a recursive key search rather than a fixed
directory depth, so this remains robust if another sample-level or
detector-level directory is added later.

Toy pairing
-----------
Toy histograms are keyed by an explicit toy *index* parsed from the trailing
integer in the histogram name (e.g. "..._combined_996" -> 996), NOT by
whatever order uproot's recursive key() listing happens to return. This
matters anywhere code assumes toys_xsec[p], toys_flux[p], toys_combined[p]
are the same underlying posterior draw p — ROOT TKey ordering is not
guaranteed to be numeric or consistent across variants/directories, so
relying on positional (list) order can silently mispair draws without
raising any error. See `_collect_toys_indexed` / `paired_toy_arrays`.
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

# Toy variants written by the C++ posterior-predictive app.
# Must match the `variant` tags passed to MakeSpectaVariations().
TOY_VARIANTS = ["xsecOnly", "fluxOnly", "combined"]

# Trailing integer in a toy histogram key, e.g. "..._combined_996" -> 996.
TOY_INDEX_RE = re.compile(r'_(\d+)$')

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


def check_covariance_decomposition(toys_xsec, toys_flux, toys_combined, det, var,
                                    n_bins_to_check=5):
    """
    Sanity check: does Var(combined) match Var(xsec) + Var(flux) + 2*Cov(xsec,flux),
    computed directly from paired toy indices? Logs per-bin comparison.

    Callers MUST pass arrays that are already row-aligned by toy index (see
    `paired_toy_arrays`) — this function no longer assumes list order implies
    pairing, but it also can't fix mispaired input, so alignment has to
    happen upstream.
    """
    if toys_xsec is None or toys_flux is None or toys_combined is None:
        log.info("Skipping covariance check for %s/%s — missing a variant.", det, var)
        return

    if not (len(toys_xsec) == len(toys_flux) == len(toys_combined)):
        log.error(
            "Toy count mismatch for %s/%s after pairing: xsec=%d flux=%d "
            "combined=%d — this should not happen if paired_toy_arrays() was "
            "used; skipping check.",
            det, var, len(toys_xsec), len(toys_flux), len(toys_combined),
        )
        return

    n = len(toys_xsec)
    xsec, flux, comb = toys_xsec, toys_flux, toys_combined

    nbins = xsec.shape[1]
    step = max(1, nbins // n_bins_to_check)

    log.info("Covariance decomposition check for %s/%s (n_paired_toys=%d):", det, var, n)
    for k in range(0, nbins, step):
        xsec_dev = xsec[:, k] - xsec[:, k].mean()
        flux_dev = flux[:, k] - flux[:, k].mean()

        var_xsec = xsec_dev.var()
        var_flux = flux_dev.var()
        cross_cov = np.mean(xsec_dev * flux_dev)

        predicted = var_xsec + var_flux + 2 * cross_cov
        measured  = comb[:, k].var()

        corr = cross_cov / np.sqrt(var_xsec * var_flux) if var_xsec > 0 and var_flux > 0 else np.nan

        log.info(
            "  bin %3d: Var(xsec)=%.3g Var(flux)=%.3g Cov=%.3g "
            "corr=%.3f | predicted=%.3g measured=%.3g ratio=%.3f",
            k, var_xsec, var_flux, cross_cov, corr,
            predicted, measured,
            predicted / measured if measured != 0 else np.nan,
        )


def find_hists_recursive(d, det, marker):
    """
    Recursively search directory `d` for histogram keys whose leaf name
    starts with "{det}_" and contains `marker`, regardless of how many
    subdirectory levels deep they live (handles det/var/sample/hist
    nesting without hardcoding depth).

    Returns a list of full key paths (usable directly as d[path]). Order is
    whatever uproot/ROOT returns — NOT assumed to be numeric or stable
    across variants. Anything that needs cross-variant pairing must parse
    and sort by toy index explicitly (see `_collect_toys_indexed`).
    """
    if d is None:
        return []
    try:
        all_keys = d.keys(recursive=True)
    except Exception as exc:
        log.debug("Could not recursively list keys under %r: %s", d, exc)
        return []

    matches = []
    for k in all_keys:
        ks = strip_cycle(k)
        leaf = ks.split("/")[-1]
        if leaf.startswith(det + "_") and marker in leaf:
            matches.append(ks)
    return matches


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
    base = name.split("/")[-1]
    for prefix in ("data_", "hOffAxis", "h"):
        if base.startswith(prefix):
            base = base[len(prefix):]
            break
    label = base.replace("_", " ").strip()
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
        ax1.step(e, a, where="post", label="MC Prediction",
                  linewidth=1.8, color=_C["orange"])

    if fakedata_asimov is not None and not no_fakedata_asimov_line:
        ax1.step(e, extend(fakedata_asimov), where="post",
                  label="Fake Data Asimov", linewidth=1.8,
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


def plot_group_bands(pdf, mean, edges, title, xlabel,
                      toys_xsec=None, toys_flux=None, toys_combined=None):
    """
    Overlay +/-1 sigma error bands from Xsec-only, Flux-only, and combined
    toy variations, each group's own mean line, a relative-uncertainty
    ratio panel, and a fraction-of-combined-variance panel (which group
    dominates the total uncertainty in each bin).

    Note: this function only needs each group internally consistent
    (mean/std computed within the group's own rows) — it does NOT require
    cross-variant row alignment, unlike the covariance decomposition check.
    """
    fig, (ax1, ax2, ax3) = plt.subplots(
        3, 1, figsize=(6.5, 8.0),
        gridspec_kw={"height_ratios": [3, 1.2, 1.2]},
        sharex=True,
    )
    fig.subplots_adjust(hspace=0.08)

    e = edges
    m = extend(mean)

    bands = {
        "Xsec only": (toys_xsec,     _C["orange"]),
        "Flux only": (toys_flux,     _C["green"]),
        "Xsec+Flux": (toys_combined, _C["blue"]),
    }

    missing = [label for label, (toys, _) in bands.items()
               if toys is None or len(toys) == 0]

    var_combined = None
    if toys_combined is not None and len(toys_combined):
        var_combined = np.var(toys_combined, axis=0)

    rel_uncs = {}
    for label, (toys, color) in bands.items():
        if toys is None or len(toys) == 0:
            continue
        ntoys = len(toys)
        std = np.std(toys, axis=0)
        grp_mean = np.mean(toys, axis=0)

        u = extend(grp_mean + std)
        l = extend(grp_mean - std)
        ax1.fill_between(e, l, u, step="post", alpha=0.22, color=color,
                          label=f"{label} " + r"$\pm1\sigma$" + f" (N={ntoys})")
        # group's own mean, so shape differences are visible, not just width
        ax1.step(e, extend(grp_mean), where="post", color=color,
                  linewidth=1.2, linestyle="--", alpha=0.9)

        rel_uncs[label] = (std, color)

        # fraction of combined variance contributed by this group
        if var_combined is not None and label != "Xsec+Flux":
            frac = np.divide(np.var(toys, axis=0), var_combined,
                              out=np.zeros_like(var_combined),
                              where=var_combined != 0)
            ax3.step(e, extend(frac), where="post", color=color,
                      linewidth=1.6, label=f"{label} / Combined")

    ax1.step(e, m, where="post", label="Fit Result", linewidth=1.8, color="black")
    ax1.set_ylabel("Events / GeV")
    ax1.set_title(title)
    if missing:
        note = "Missing: " + ", ".join(missing)
        ax1.text(0.98, 0.02, note, transform=ax1.transAxes,
                  ha="right", va="bottom", fontsize=7, color="gray", style="italic")
    ax1.legend(loc="best", framealpha=0.9, edgecolor="0.7", fontsize=8)
    plt.setp(ax1.get_xticklabels(), visible=False)

    # ── middle panel: fractional uncertainty per group ─────────────────────
    for label, (std, color) in rel_uncs.items():
        rel = np.divide(std, mean, out=np.zeros_like(std), where=mean != 0)
        ax2.step(e, extend(rel), where="post", color=color,
                  linewidth=1.6, label=label)
    ax2.axhline(0, color="black", linewidth=0.6, linestyle=":", alpha=0.5)
    ax2.set_ylabel(r"$\sigma\,/\,\mu$")
    ax2.legend(loc="upper right", framealpha=0.9, edgecolor="0.7", fontsize=7)
    plt.setp(ax2.get_xticklabels(), visible=False)

    # ── bottom panel: variance fraction ─────────────────────────────────────
    ax3.axhline(1, color="black", linewidth=0.6, linestyle=":", alpha=0.5)
    ax3.set_ylabel("Var(group) /\nVar(combined)")
    ax3.set_xlabel(xlabel)
    if var_combined is None:
        ax3.text(0.5, 0.5, "No combined toys \u2014 fraction unavailable",
                  transform=ax3.transAxes, ha="center", va="center",
                  fontsize=8, color="gray")
    else:
        ax3.legend(loc="upper right", framealpha=0.9, edgecolor="0.7", fontsize=7)

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


# ── toy collection ───────────────────────────────────────────────────────────

def _try_add_toy(h, key, edges, bucket):
    """
    Read histogram `h` (from key `key`) and, if valid, insert it into
    `bucket` — a {toy_index: value_array} dict — keyed by the toy index
    parsed from the trailing integer in `key`.

    This is the load-bearing fix: toys are identified by their numeric
    index, not by the order they were appended while iterating keys().
    ROOT/uproot key ordering is not guaranteed to be numeric or consistent
    across variant subdirectories, so building arrays via list-append +
    np.vstack (as before) can silently mispair xsecOnly/fluxOnly/combined
    draws that are supposed to correspond to the same posterior sample.
    """
    m = TOY_INDEX_RE.search(key)
    if not m:
        log.debug("Could not parse toy index from %s – skipped.", key)
        return
    idx = int(m.group(1))

    try:
        if is_2d(h):
            return
        v, e = get_1d(h)
        if not np.allclose(e, edges):
            log.debug("Edge mismatch for toy %s – skipped.", key)
            return
        if idx in bucket:
            log.warning("Duplicate toy index %d seen for %s – overwriting "
                        "previous entry.", idx, key)
        bucket[idx] = v / np.diff(e)
    except Exception as exc:
        log.debug("Could not read toy %s: %s", key, exc)


def _collect_toys_indexed(f, det, var, edges):
    """
    Return {variant: {toy_index: value_array}} for each entry in
    TOY_VARIANTS, scanning det/var recursively for histogram keys tagged
    "_{var}_posterior_toy_{variant}_". This is the canonical, index-keyed
    collection that both `collect_toys_by_variant` (unpaired, per-group
    plots) and `paired_toy_arrays` (paired, cross-variant checks) build on.
    """
    toys_by_variant = {v: {} for v in TOY_VARIANTS}

    var_subdir = safe_get(f, f"{det}/{var}")
    if var_subdir is None:
        return toys_by_variant

    for variant in TOY_VARIANTS:
        toy_marker = f"_{var}_posterior_toy_{variant}_"
        for full_path in find_hists_recursive(var_subdir, det, toy_marker):
            _try_add_toy(var_subdir[full_path], full_path, edges,
                         toys_by_variant[variant])

    return toys_by_variant


def collect_toys_by_variant(f, det, var, edges):
    """
    Return {variant: np.vstack(toys) or None} for each entry in
    TOY_VARIANTS. Rows are sorted by toy index within each variant (so
    output is deterministic and reproducible run-to-run), but rows are
    NOT guaranteed aligned across variants — a given row index means
    "the toy with the k-th smallest index available in this variant",
    which may differ between xsecOnly/fluxOnly/combined if any variant is
    missing toys the others have.

    Safe to use for: mean/std bands, percentile envelopes, per-group
    plots (plot_1d, plot_group_bands) — anything that only needs each
    variant's array to be internally self-consistent.

    NOT safe to use for: anything that pairs rows across variants by
    position (e.g. covariance/correlation between xsec and flux toys) —
    use `paired_toy_arrays` for that instead.
    """
    indexed = _collect_toys_indexed(f, det, var, edges)
    result = {}
    for v, idx_map in indexed.items():
        result[v] = None if not idx_map else np.vstack(
            [idx_map[i] for i in sorted(idx_map)]
        )
    return result


def paired_toy_arrays(f, det, var, edges, variants=TOY_VARIANTS):
    """
    Like `collect_toys_by_variant`, but restricted to the intersection of
    toy indices present in ALL requested `variants`, with matching rows
    guaranteed to correspond to the same posterior draw across variants.

    Use this wherever pairing matters, e.g. feeding
    `check_covariance_decomposition`, which assumes
    toys_xsec[p] <-> toys_flux[p] <-> toys_combined[p] are the same draw.

    Returns (arrays, sorted_common_indices) where `arrays` is
    {variant: np.vstack(...) or None}. If fewer than 2 variants have any
    toys, or the intersection of indices is empty, all entries are None
    and sorted_common_indices is [].
    """
    indexed = _collect_toys_indexed(f, det, var, edges)
    present = {v: indexed[v] for v in variants if indexed.get(v)}

    if len(present) < 2:
        return {v: None for v in variants}, []

    common_idx = set.intersection(*(set(m) for m in present.values()))
    if not common_idx:
        log.warning("%s/%s: no toy index common to variants %s — cannot pair.",
                    det, var, list(present))
        return {v: None for v in variants}, []

    for v, idx_map in present.items():
        dropped = set(idx_map) - common_idx
        if dropped:
            log.warning(
                "%s/%s [%s]: dropping %d toy(s) not present in all paired "
                "variants: %s%s",
                det, var, v, len(dropped), sorted(dropped)[:10],
                " ..." if len(dropped) > 10 else "",
            )

    sorted_idx = sorted(common_idx)
    arrays = {
        v: (np.vstack([indexed[v][i] for i in sorted_idx]) if v in present else None)
        for v in variants
    }
    return arrays, sorted_idx


# ── CLI ──────────────────────────────────────────────────────────────────────

def _str2bool(s):
    return s.lower() not in ("false", "0", "no", "off")


def parse_args():
    p = argparse.ArgumentParser(
        description="Posterior predictive plots from a ROOT file produced by the neutrino fit."
    )
    p.add_argument("root_file", help="Path to the posterior-predictive ROOT file")
    p.add_argument("data_file", nargs="?", default=None,
                    help="Optional data ROOT file (ignored if --fakedata=False)")
    p.add_argument("output_pdf", nargs="?", default=None,
                    help="Output PDF path (default: <root_dir>/plots.pdf)")
    p.add_argument("--fakedata", type=_str2bool, default=True,
                    help="Overlay fake data points and Fake Data Asimov line "
                         "(default: True). Use --fakedata=False for an "
                         "Asimov-only fit with no fake data comparison.")
    return p.parse_args()


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    logging.basicConfig(level=logging.INFO, format="%(levelname)s  %(message)s")

    args = parse_args()

    root_path = args.root_file
    root_dir  = os.path.dirname(os.path.abspath(root_path))

    use_fakedata = args.fakedata
    log.info("Fake data comparison: %s", "ON" if use_fakedata else "OFF")

    if use_fakedata:
        data_root_path = args.data_file or os.path.join(root_dir, "fake_data.root")
    else:
        data_root_path = None

    out_path = args.output_pdf or os.path.join(root_dir, "plots.pdf")

    f = uproot.open(root_path)

    top_keys = [strip_cycle(k) for k in f.keys(recursive=False)]
    detectors = [
        k for k in top_keys
        if k not in ("Asimov", "AsimovFakeData") and safe_get(f, k) is not None
    ]
    log.info("Found detectors: %s", detectors)

    asimov_dir = safe_get(f, "Asimov")

    # Load data histograms from the (optional) data ROOT file
    data_hists = load_data_histograms(data_root_path) if use_fakedata else {}
    if data_hists:
        log.info("Data histograms loaded: %s", list(data_hists.keys()))
    elif use_fakedata:
        log.info("No data histograms loaded (no data file provided or none found).")

    variables = list(VAR_META.keys())

    with PdfPages(out_path) as pdf:

        # ── 1-D loop ──────────────────────────────────────────────────────
        for det in detectors:
            for var in variables:

                if asimov_dir is None:
                    log.warning("No Asimov directory in file – skipping 1-D plots.")
                    break

                asimov_keys = find_hists_recursive(asimov_dir, det, f"_{var}_Asimov")

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

                # load FakeData Asimov (skipped entirely for Asimov-only fits)
                fakedata_asimov_vals = None
                if use_fakedata:
                    fakedata_asimov_dir = safe_get(f, "AsimovFakeData")
                    if fakedata_asimov_dir is not None:
                        fd_keys = find_hists_recursive(fakedata_asimov_dir, det, f"_{var}_AsimovFakeData")
                        if fd_keys:
                            h_fd = fakedata_asimov_dir[fd_keys[0]]
                            if not is_2d(h_fd):
                                fd_vals, _ = get_1d(h_fd)
                                fakedata_asimov_vals = fd_vals / bin_widths

                # ── collect toys, split by variant (index-keyed, not order-keyed) ──
                toys_by_variant = collect_toys_by_variant(f, det, var, edges)
                toys_arr      = toys_by_variant["combined"]
                toys_xsec_arr = toys_by_variant["xsecOnly"]
                toys_flux_arr = toys_by_variant["fluxOnly"]

                if toys_arr is None:
                    log.warning("No combined toys found for %s/%s.", det, var)
                    continue

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

                # ── standard combined/total plot ──────────────────────────
                plot_1d(pdf, **common_kwargs)

                # ── Xsec-only systematics, separate page ───────────────────
                if toys_xsec_arr is not None:
                    mean_xsec = np.mean(toys_xsec_arr, axis=0)
                    std_xsec  = np.std(toys_xsec_arr, axis=0)
                    plot_1d(pdf,
                            asimov=asimov_vals, mean=mean_xsec, std=std_xsec,
                            edges=edges, title=f"{var} — Xsec-only systematics",
                            xlabel=xlabel, toys=toys_xsec_arr,
                            data_hists=None, is_reco=False,
                            fakedata_asimov=None)
                else:
                    log.warning("No Xsec-only toys for %s/%s — skipping Xsec page.", det, var)

                # ── Flux-only systematics, separate page ────────────────────
                if toys_flux_arr is not None:
                    mean_flux = np.mean(toys_flux_arr, axis=0)
                    std_flux  = np.std(toys_flux_arr, axis=0)
                    plot_1d(pdf,
                            asimov=asimov_vals, mean=mean_flux, std=std_flux,
                            edges=edges, title=f"{var} — Flux-only systematics",
                            xlabel=xlabel, toys=toys_flux_arr,
                            data_hists=None, is_reco=False,
                            fakedata_asimov=None)
                else:
                    log.warning("No Flux-only toys for %s/%s — skipping Flux page.", det, var)

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

                # ── group-decomposed error bands (Xsec / Flux / combined) ──
                if toys_xsec_arr is not None or toys_flux_arr is not None:
                    plot_group_bands(pdf, mean, edges,
                                      title=f"{var} — error band by group",
                                      xlabel=xlabel,
                                      toys_xsec=toys_xsec_arr,
                                      toys_flux=toys_flux_arr,
                                      toys_combined=toys_arr)
                else:
                    log.warning(
                        "No Xsec-only/Flux-only toys for %s/%s — rerun the "
                        "posterior-predictive app to generate them.", det, var)

                log.info("1-D done: %s / %s  (%d combined toys)", det, var, len(toys_arr))

                # ── covariance decomposition check — REQUIRES paired rows ──
                # Uses paired_toy_arrays() (index-intersection across
                # variants), not the toys_xsec_arr/toys_flux_arr/toys_arr
                # from collect_toys_by_variant() above, since those are only
                # guaranteed self-consistent within a variant, not aligned
                # across variants.
                if toys_xsec_arr is not None and toys_flux_arr is not None:
                    paired, common_idx = paired_toy_arrays(
                        f, det, var, edges,
                        variants=["xsecOnly", "fluxOnly", "combined"],
                    )
                    if all(paired[v] is not None for v in
                           ("xsecOnly", "fluxOnly", "combined")):
                        log.info("%s/%s: %d toy indices common across all "
                                 "three variants for paired covariance check.",
                                 det, var, len(common_idx))
                        check_covariance_decomposition(
                            paired["xsecOnly"], paired["fluxOnly"],
                            paired["combined"], det, var,
                        )
                    else:
                        log.warning(
                            "%s/%s: could not build a fully-paired toy set "
                            "across xsecOnly/fluxOnly/combined — skipping "
                            "covariance decomposition check.", det, var,
                        )

        # ── 2-D summaries ─────────────────────────────────────────────────
        plot_2d_summaries(f, detectors, pdf)

    log.info("Saved → %s", out_path)


if __name__ == "__main__":
    main()