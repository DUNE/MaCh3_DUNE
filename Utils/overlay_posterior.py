#!/usr/bin/env python3
"""
Overlay posterior-predictive plots from several fit ROOT files on one set of axes.

Same three-panel layout as plot_posterior.py:
    top     Events / GeV  - Fit Result + ±1σ (posterior) [+ 5–95th pct], CV, Fake Data Asimov
    middle  σ_Fit / μ     - one line per fit
    bottom  ratio         - CV/Fit (default), FakeData/Fit, or Fit/Fit_ref, with ±1σ band

Each input is   file.root[:Label[:Detector]]
    Label     legend name (default: file name)
    Detector  only use this top-level detector dir from that file. Give it when the
              files name their detectors differently and you want them on one plot.

Usage
-----
    # same detector/variable in both files -> one overlay page per (detector, variable)
    python overlay_posterior.py fitA.root:"Fake-data fit" fitB.root:"Nominal fit" -o overlay.pdf

    # force a pairing of differently named detector dirs
    python overlay_posterior.py fitA.root:"FD fit":ND fitB.root:"OffAxis0m":ND_OffAxis0m

    # shape-only comparison (area-normalise each fit), compare to the first fit
    python overlay_posterior.py a.root b.root --normalise --ratio first

    # restrict variables, drop the percentile band
    python overlay_posterior.py a.root b.root --vars TrueNeutrinoEnergy --no-pct
"""

import argparse
import logging
import os
import re

import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

log = logging.getLogger(__name__)

# ── style (identical to plot_posterior.py) ───────────────────────────────────
THESIS_RC = {
    "font.family": "serif",
    "font.serif": ["STIXTwo Text", "DejaVu Serif", "Times New Roman"],
    "mathtext.fontset": "stix",
    "font.size": 11, "axes.titlesize": 11, "axes.labelsize": 11,
    "xtick.labelsize": 10, "ytick.labelsize": 10, "legend.fontsize": 8,
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "figure.dpi": 150,
    "lines.linewidth": 1.8, "patch.linewidth": 1.0,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
    "xtick.minor.visible": True, "ytick.minor.visible": True,
    "axes.grid": True, "grid.alpha": 0.3, "grid.linewidth": 0.6,
    "axes.linewidth": 0.8,
}
plt.rcParams.update(THESIS_RC)

# One Okabe-Ito colour per fit. Orange is left out: it is kept for a shared CV line.
FIT_COLOURS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#56B4E9", "#000000"]
CV_COLOUR = "#E69F00"
LINESTYLES = ["-", (0, (6, 2)), (0, (2, 1.5)), (0, (6, 2, 1.5, 2))]

VAR_META = {
    "TrueNeutrinoEnergy": {"xlabel": "True Neutrino Energy (GeV)"},
    "RecoNeutrinoEnergy": {"xlabel": "Reconstructed Neutrino Energy (GeV)"},
    "Enubias":            {"xlabel": r"$E_\nu$ bias (GeV)", "xlim_min": -0.2},
}


# ── ROOT helpers (same conventions as plot_posterior.py) ─────────────────────
def strip_cycle(k):
    return k.split(";")[0]


def safe_get(f, path):
    try:
        return f[path]
    except Exception as exc:
        log.debug("safe_get(%s): %s", path, exc)
        return None


def is_2d(h):
    return len(h.axes) == 2


def extend(v):
    return np.append(v, v[-1])


TOY_RE = re.compile(r"_posterior_toy_(\d+)")


def _th1s(directory):
    """All 1-D histograms under `directory` (recursive) -> {path_without_cycle: key}."""
    if directory is None or not hasattr(directory, "classnames"):
        return {}
    out = {}
    for k, cls in directory.classnames(recursive=True).items():
        if cls.startswith("TH1"):
            out.setdefault(strip_cycle(k), k)
    return out


def _has_sample(path, sample):
    """True if `sample` appears in `path` as a whole token (OffAxis2m != OffAxis20m)."""
    return re.search(rf"(^|[_/]){re.escape(sample)}($|[_/])", path) is not None


def list_samples(f, det, var):
    """Sample sub-folders of <det>/<var> that hold toys, in file order."""
    d = safe_get(f, f"{det}/{var}")
    if d is None:
        return []
    samples = []
    for path in _th1s(d):
        if "/" in path and TOY_RE.search(path):
            s = path.split("/")[0]
            if s not in samples:
                samples.append(s)
    return samples


def _cv_like(f, top, det, var, suffix, samples, exclude=None):
    """
    Find the CV-type histogram(s) in directory `top` for this det/var, summed over
    `samples` (or the single total histogram when samples is empty).
    Returns (values, edges) or None.
    """
    hists = {p: k for p, k in _th1s(safe_get(f, top)).items()
             if f"_{var}{suffix}" in p.split("/")[-1]
             and not (exclude and exclude in p.split("/")[-1])}
    if not hists:
        return None
    with_det = {p: k for p, k in hists.items() if p.split("/")[-1].startswith(det + "_") or
                p.startswith(det + "/")}
    hists = with_det or hists
    d = f[top]

    def read(k):
        h = d[k]
        return h.values().astype(float), h.axes[0].edges()

    if not samples:
        # prefer a histogram that does not belong to any sample folder
        p = sorted(hists, key=len)[0]
        return read(hists[p])

    total, edges = None, None
    for s in samples:
        match = [p for p in hists if _has_sample(p, s)]
        if not match:
            # no per-sample histogram: fall back to one total histogram if that is all there is
            if len(hists) == 1:
                log.debug("%s: no per-sample %s for %s; using the only one", top, var, s)
                return read(next(iter(hists.values())))
            log.warning("%s: no histogram for sample %s / %s", top, s, var)
            return None
        v, e = read(hists[sorted(match, key=len)[0]])
        if edges is not None and (len(e) != len(edges) or not np.allclose(e, edges)):
            log.warning("%s: binning differs between samples for %s – can't sum.", top, var)
            return None
        total = v if total is None else total + v
        edges = e
    return total, edges


def load_fit(f, det, var, sample=None):
    """
    Everything needed for one (detector, variable[, sample]) from an open ROOT file.

    sample=None  -> toys stored directly in <det>/<var> (old layout); if the toys live in
                    per-sample folders instead, all samples are summed toy-by-toy.
    sample=name  -> only <det>/<var>/<name>.
    Returns a dict of per-GeV arrays, or None if nothing usable was found.
    """
    # ── toys, grouped by sample and toy index ────────────────────────────────
    by_sample = {}          # sample (None = no folder) -> {toy_idx: (values, edges)}
    for top in (f"{det}/{var}", det):
        d = safe_get(f, top)
        if d is None:
            continue
        paths = _th1s(d)
        if top == det:      # old layout also kept toys straight in <det>; don't recurse there
            paths = {p: k for p, k in paths.items() if "/" not in p}
        for p, k in paths.items():
            name = p.split("/")[-1]
            m = TOY_RE.search(name)
            if not m or f"_{var}_posterior_toy_" not in name:
                continue
            s = p.split("/")[0] if "/" in p else None
            if sample is not None and s != sample:
                continue
            try:
                h = d[k]
                by_sample.setdefault(s, {})[int(m.group(1))] = (
                    h.values().astype(float), h.axes[0].edges())
            except Exception as exc:
                log.debug("toy %s unreadable: %s", p, exc)
    if not by_sample:
        return None

    samples = [s for s in by_sample if s is not None]
    if sample is None and None in by_sample:
        samples, by_sample = [], {None: by_sample[None]}    # total already stored

    # sum toy i across samples (same throw), keeping indices present in every sample
    idx = sorted(set.intersection(*(set(t) for t in by_sample.values())))
    if not idx:
        log.warning("%s/%s: no toy index common to all samples.", det, var)
        return None
    edges = next(iter(by_sample.values()))[idx[0]][1]
    toys = []
    for i in idx:
        parts = [by_sample[s][i] for s in by_sample]
        if any(len(e) != len(edges) or not np.allclose(e, edges) for _, e in parts):
            log.warning("%s/%s: binning differs between samples – can't sum.", det, var)
            return None
        toys.append(sum(v for v, _ in parts))
    w = np.diff(edges)
    toys = np.vstack(toys) / w

    # ── CV and fake-data Asimov, matched to the same samples ─────────────────
    cv = _cv_like(f, "Asimov", det, var, "_Asimov", samples, exclude="AsimovFakeData")
    if cv is None or len(cv[1]) != len(edges) or not np.allclose(cv[1], edges):
        log.warning("%s/%s%s: no CV histogram with matching binning.",
                    det, var, f"/{sample}" if sample else "")
        return None
    fd = _cv_like(f, "AsimovFakeData", det, var, "_AsimovFakeData", samples)
    fd_vals = None
    if fd is not None and len(fd[1]) == len(edges) and np.allclose(fd[1], edges):
        fd_vals = fd[0] / w

    return dict(edges=edges, cv=cv[0] / w, fakedata=fd_vals, toys=toys,
                mean=toys.mean(0), std=toys.std(0),
                p5=np.percentile(toys, 5, 0), p95=np.percentile(toys, 95, 0),
                samples=samples)


def print_structure(f, path, max_show=4):
    """Short tree of a file, collapsing long runs of toys, for --list."""
    print(f"\n== {path}")
    groups = {}
    for p in _th1s(f):
        key = TOY_RE.sub("_posterior_toy_<N>", p)
        groups[key] = groups.get(key, 0) + 1
    shown = {}
    for key, n in sorted(groups.items()):
        folder = key.rsplit("/", 1)[0] if "/" in key else ""
        shown[folder] = shown.get(folder, 0) + 1
        if shown[folder] <= max_show:
            print(f"  {key}" + (f"   (x{n})" if n > 1 else ""))
        elif shown[folder] == max_show + 1:
            print(f"  {folder}/ ...")


def list_detectors(f):
    return [strip_cycle(k) for k in f.keys(recursive=False)
            if strip_cycle(k) not in ("Asimov", "AsimovFakeData")
            and hasattr(safe_get(f, strip_cycle(k)), "keys")]


def normalise(fit):
    """Area-normalise every histogram in a fit to the posterior-mean integral."""
    w = np.diff(fit["edges"])
    n = np.sum(fit["mean"] * w)
    if n <= 0:
        return fit
    out = dict(fit)
    for k in ("cv", "fakedata", "toys", "mean", "std", "p5", "p95"):
        if out[k] is not None:
            out[k] = out[k] / n
    return out


def _div(a, b):
    return np.divide(a, b, out=np.full_like(a, np.nan, dtype=float), where=b != 0)


# ── the overlay plot ─────────────────────────────────────────────────────────
def plot_overlay(pdf, fits, title, xlabel, *, show_pct=True, show_cv=True,
                 show_fakedata=True, ratio="cv", xlim_min=None, norm=False):
    """
    fits: list of (label, fit_dict). Every fit must share the binning of the first;
    fits with different edges are skipped with a warning.
    """
    ref_edges = fits[0][1]["edges"]
    ok = []
    for label, fit in fits:
        if len(fit["edges"]) != len(ref_edges) or not np.allclose(fit["edges"], ref_edges):
            log.warning("'%s' has different binning for %s – left off this page.", label, title)
            continue
        ok.append((label, fit))
    fits = ok
    if len(fits) < 1:
        return
    e = ref_edges

    fig, (ax1, ax2, ax3) = plt.subplots(
        3, 1, figsize=(6.5, 7.5), sharex=True,
        gridspec_kw={"height_ratios": [3, 1, 1]})
    fig.subplots_adjust(hspace=0.05)

    # If every CV is the same, draw it once in orange; otherwise one per fit.
    cvs = [f["cv"] for _, f in fits]
    shared_cv = all(np.allclose(c, cvs[0], rtol=1e-6, atol=0) for c in cvs[1:])

    top_handles = []
    ref_mean, ref_std = fits[0][1]["mean"], fits[0][1]["std"]

    for i, (label, fit) in enumerate(fits):
        c = FIT_COLOURS[i % len(FIT_COLOURS)]
        m, s = fit["mean"], fit["std"]

        # top panel
        if show_pct:
            ax1.fill_between(e, extend(fit["p5"]), extend(fit["p95"]), step="post",
                             color=c, alpha=0.10, lw=0)
        ax1.fill_between(e, extend(m - s), extend(m + s), step="post",
                         color=c, alpha=0.28, lw=0)
        ax1.step(e, extend(m), where="post", color=c, lw=1.8)
        top_handles.append((Line2D([], [], color=c, lw=1.8), f"{label}: Fit Result"))

        if show_cv and not shared_cv:
            ax1.step(e, extend(fit["cv"]), where="post", color=c, lw=1.3, ls=(0, (1, 1.2)))
            top_handles.append((Line2D([], [], color=c, lw=1.3, ls=(0, (1, 1.2))),
                                f"{label}: CV Prediction"))
        if show_fakedata and fit["fakedata"] is not None:
            ax1.step(e, extend(fit["fakedata"]), where="post", color=c, lw=1.5, ls=(0, (5, 2)))
            top_handles.append((Line2D([], [], color=c, lw=1.5, ls=(0, (5, 2))),
                                f"{label}: Fake Data Asimov"))

        # middle panel: fractional posterior width
        ax2.step(e, extend(_div(s, m)), where="post", color=c, lw=1.6, label=label)

        # bottom panel
        if ratio == "first":
            if i == 0:
                continue
            num = m
            den, den_s = ref_mean, ref_std
            r_label = f"{label} / {fits[0][0]}"
        elif ratio == "fakedata":
            if fit["fakedata"] is None:
                continue
            num, den, den_s = fit["fakedata"], m, s
            r_label = f"{label}: Fake Data / Fit"
        else:  # cv
            num, den, den_s = fit["cv"], m, s
            r_label = f"{label}: CV / Fit"

        r = _div(num, den)
        r_up = _div(num, den - den_s)   # num/(μ-σ) is the upper edge
        r_dn = _div(num, den + den_s)
        ax3.fill_between(e, extend(r_dn), extend(r_up), step="post",
                         color=c, alpha=0.22, lw=0)
        ax3.step(e, extend(r), where="post", color=c, lw=1.6, label=r_label)

    if show_cv and shared_cv:
        ax1.step(e, extend(cvs[0]), where="post", color=CV_COLOUR, lw=1.6, zorder=4)
        top_handles.insert(0, (Line2D([], [], color=CV_COLOUR, lw=1.6), "CV Prediction"))

    # band key, drawn in grey so it reads for every colour
    top_handles.append((Patch(color="0.5", alpha=0.35, lw=0), r"$\pm1\sigma$ (posterior)"))
    if show_pct:
        top_handles.append((Patch(color="0.5", alpha=0.15, lw=0), "5–95th pct"))

    # cosmetics
    ax1.set_title(title)
    ax1.set_ylabel("Normalised events / GeV" if norm else "Events / GeV")
    ax1.ticklabel_format(axis="y", style="sci", scilimits=(-3, 4))
    h, lab = zip(*top_handles)
    ax1.legend(h, lab, loc="upper right", framealpha=0.9, edgecolor="0.7",
               ncol=2 if len(h) > 6 else 1)
    # headroom so the legend sits above the peak
    ymax = max(np.nanmax(f["p95"] if show_pct else f["mean"] + f["std"]) for _, f in fits)
    ax1.set_ylim(0, ymax * (1.15 + 0.06 * len(h) / (2 if len(h) > 6 else 1)))

    ax2.set_ylabel(r"$\sigma_\mathrm{Fit}\,/\,\mu$")
    ax2.set_ylim(bottom=0)
    ax2.legend(loc="upper right", framealpha=0.9, edgecolor="0.7")

    ax3.axhline(1, color="black", lw=0.8, ls="--", alpha=0.7)
    ax3.set_ylabel({"cv": "MC/Fit", "fakedata": "FD/Fit", "first": "Fit/Ref"}[ratio])
    ax3.set_xlabel(xlabel)
    if ax3.get_legend_handles_labels()[0]:
        ax3.legend(loc="upper right", framealpha=0.9, edgecolor="0.7")

    ax3.set_xlim(e[0] if xlim_min is None else xlim_min, e[-1])
    fig.align_ylabels([ax1, ax2, ax3])
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


# ── main ─────────────────────────────────────────────────────────────────────
def parse_input(spec):
    """file.root[:Label[:Detector[:Sample]]]"""
    parts = spec.split(":")
    path = parts[0]
    get = lambda i: parts[i] if len(parts) > i and parts[i] else None
    label = get(1) or os.path.splitext(os.path.basename(path))[0]
    return path, label, get(2), get(3)


def pretty_sample(s):
    return s.replace("_", " ") if s else ""


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+", help="file.root[:Label[:Detector[:Sample]]]")
    ap.add_argument("-o", "--output", default="overlay.pdf")
    ap.add_argument("--sample", default=None,
                    help="sample folder to use in files that have per-sample folders, "
                         "e.g. OffAxis0m_numuCC_numode (on-axis). 'all' sums every sample "
                         "(the default). Files without sample folders use their total.")
    ap.add_argument("--list", action="store_true",
                    help="print each file's histogram layout and exit")
    ap.add_argument("--vars", nargs="+", default=list(VAR_META),
                    help="variables to plot (default: %(default)s)")
    ap.add_argument("--ratio", choices=["cv", "fakedata", "first"], default="cv",
                    help="bottom panel: CV/Fit, FakeData/Fit, or each Fit / first Fit")
    ap.add_argument("--normalise", action="store_true",
                    help="area-normalise each fit (shape comparison when rates differ)")
    ap.add_argument("--no-pct", action="store_true", help="hide the 5–95th pct bands")
    ap.add_argument("--no-cv", action="store_true", help="hide CV prediction lines")
    ap.add_argument("--no-fakedata", action="store_true", help="hide Fake Data Asimov lines")
    ap.add_argument("-v", "--verbose", action="store_true")
    a = ap.parse_args()

    logging.basicConfig(level=logging.DEBUG if a.verbose else logging.INFO,
                        format="%(levelname)s  %(message)s")
    logging.getLogger("fontTools").setLevel(logging.WARNING)   # PDF font-subsetting chatter

    inputs = []
    for spec in a.inputs:
        path, label, det, sample = parse_input(spec)
        f = uproot.open(path)
        if a.list:
            print_structure(f, path)
            continue
        dets = [det] if det else list_detectors(f)
        log.info("%s  (%s): detectors %s", label, path, dets)
        inputs.append((label, f, dets, sample or a.sample))
    if a.list:
        return

    # Pages: if every input names a detector explicitly, they are paired into one
    # page per variable; otherwise one page per detector found in all files.
    forced = all(parse_input(s)[2] for s in a.inputs)
    if forced:
        groups = [[(lab, f, d[0], s) for lab, f, d, s in inputs]]
    else:
        common = [d for d in inputs[0][2] if all(d in ds for _, _, ds, _ in inputs[1:])]
        if not common:
            log.error("No detector directory is common to all files. "
                      "Pair them explicitly with file.root:Label:Detector.")
            return
        groups = [[(lab, f, d, s) for lab, f, _, s in inputs] for d in common]

    n_pages = 0
    with PdfPages(a.output) as pdf:
        for group in groups:
            dets = sorted({d for _, _, d, _ in group})
            det_title = dets[0] if len(dets) == 1 else " vs ".join(dets)
            for var in a.vars:
                fits, used = [], set()
                for lab, f, d, s in group:
                    have = list_samples(f, d, var)
                    if s and s != "all" and have:
                        if s not in have:
                            log.warning("'%s' has no sample %s for %s/%s (has: %s) – skipped.",
                                        lab, s, d, var, ", ".join(have))
                            continue
                        fit = load_fit(f, d, var, s)
                        used.add(s)
                    else:
                        fit = load_fit(f, d, var)          # total, or sum of all samples
                        if have:
                            used.add("all samples")
                    if fit is None:
                        log.warning("No CV/toys for %s / %s in '%s' – skipped.", d, var, lab)
                        continue
                    fits.append((lab, normalise(fit) if a.normalise else fit))
                if not fits:
                    continue
                meta = VAR_META.get(var, {})
                title = f"{var} — {det_title}"
                if used:
                    title += " / " + ", ".join(pretty_sample(u) for u in sorted(used))
                plot_overlay(pdf, fits, title,
                             meta.get("xlabel", var),
                             show_pct=not a.no_pct, show_cv=not a.no_cv,
                             show_fakedata=not a.no_fakedata, ratio=a.ratio,
                             xlim_min=meta.get("xlim_min"), norm=a.normalise)
                n_pages += 1
                log.info("page: %s / %s  (%s)", det_title, var,
                         ", ".join(f"{l}: {len(x['toys'])} toys" for l, x in fits))

    log.info("Saved %d page(s) → %s", n_pages, a.output)


if __name__ == "__main__":
    main()