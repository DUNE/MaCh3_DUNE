"""
plot_sigma_variations.py
========================
Thesis-quality sigma variation plots from JUNE25th_ssigmavar.root

Structure assumed
-----------------
  <Systematic>/
      <Sample>/
          Variation_0   → -3σ
          Variation_1   → -1σ
          Variation_2   → nominal
          Variation_3   → +1σ
          Variation_4   → +3σ

Output
------
  Single multipage PDF with one page per systematic × sample.
  Each page has two panels:
    TOP    – absolute histograms (nominal + ±1σ + ±3σ lines)
             with smooth stat error band on the nominal
    BOTTOM – ratio to nominal with shaded ±1σ / ±3σ bands
             and smooth stat uncertainty on the ratio
             (y-axis auto-scales to data each plot)

Usage
-----
  python plot_sigma_variations.py                          # all systematics, all samples
  python plot_sigma_variations.py --systematics EMEnergyScale TrackedMuonEnergyScale
  python plot_sigma_variations.py --samples OffAxis0m_numuCC_numode
  python plot_sigma_variations.py --out my_results.pdf
  python plot_sigma_variations.py --xlabel "Muon Momentum (GeV/c)"
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.backends.backend_pdf import PdfPages

# ── Try uproot (pure-Python); fall back to ROOT if available ──────────────────
try:
    import uproot
    _BACKEND = "uproot"
except ImportError:
    try:
        import ROOT  # noqa: F401
        _BACKEND = "ROOT"
    except ImportError:
        sys.exit(
            "ERROR: neither 'uproot' nor 'ROOT' is importable.\n"
            "  Install uproot:  pip install uproot\n"
            "  or activate your ROOT environment."
        )

# ─────────────────────────────────────────────────────────────────────────────
# Variation index → label / style
# ─────────────────────────────────────────────────────────────────────────────
VARIATIONS = {
    0: dict(label=r"$-3\sigma$", color="#d62728", ls="--",  lw=1.4, zorder=2),
    1: dict(label=r"$-1\sigma$", color="#ff7f0e", ls="-.",  lw=1.4, zorder=3),
    2: dict(label="Nominal",     color="#1f77b4", ls="-",   lw=2.0, zorder=5),
    3: dict(label=r"$+1\sigma$", color="#2ca02c", ls="-.",  lw=1.4, zorder=3),
    4: dict(label=r"$+3\sigma$", color="#9467bd", ls="--",  lw=1.4, zorder=2),
}

BAND_1S_COLOR = "#2ca02c"
BAND_3S_COLOR = "#9467bd"
BAND_1S_ALPHA = 0.18
BAND_3S_ALPHA = 0.10


# ─────────────────────────────────────────────────────────────────────────────
# Publication style
# ─────────────────────────────────────────────────────────────────────────────
def apply_thesis_style():
    plt.rcParams.update({
        "font.family":         "serif",
        "font.serif":          ["Computer Modern Roman", "DejaVu Serif"],
        "font.size":           12,
        "axes.titlesize":      13,
        "axes.labelsize":      13,
        "xtick.labelsize":     11,
        "ytick.labelsize":     11,
        "legend.fontsize":     10,
        "figure.dpi":          150,
        "savefig.dpi":         300,
        "savefig.bbox":        "tight",
        "axes.linewidth":      1.2,
        "xtick.major.width":   1.1,
        "ytick.major.width":   1.1,
        "xtick.minor.width":   0.8,
        "ytick.minor.width":   0.8,
        "xtick.direction":     "in",
        "ytick.direction":     "in",
        "xtick.top":           True,
        "ytick.right":         True,
        "xtick.minor.visible": True,
        "ytick.minor.visible": True,
        "axes.grid":           True,
        "grid.alpha":          0.3,
        "grid.linewidth":      0.6,
    })


# ─────────────────────────────────────────────────────────────────────────────
# Smooth error band  (fill_between at bin centres — no pointy box edges)
# ─────────────────────────────────────────────────────────────────────────────
def draw_smooth_errors(ax, edges, counts, errors, color, alpha=0.35, zorder=6):
    """
    Draw a smooth filled band [counts-errors, counts+errors] using fill_between
    evaluated at bin centres.  Much smoother than per-bin rectangles.
    """
    centres = 0.5 * (edges[:-1] + edges[1:])
    valid   = np.isfinite(counts) & np.isfinite(errors) & (errors > 0)
    y_lo    = np.where(valid, counts - errors, np.nan)
    y_hi    = np.where(valid, counts + errors, np.nan)
    ax.fill_between(centres, y_lo, y_hi,
                    alpha=alpha, color=color, linewidth=0, zorder=zorder)


# ─────────────────────────────────────────────────────────────────────────────
# Bin-width division helper
# ─────────────────────────────────────────────────────────────────────────────
def per_bw(counts, errors, edges):
    """Divide counts and errors by bin widths."""
    widths = np.diff(edges)
    return counts / widths, errors / widths


# ─────────────────────────────────────────────────────────────────────────────
# ROOT-file helpers  (uproot backend)
# ─────────────────────────────────────────────────────────────────────────────
def load_file_uproot(path):
    return uproot.open(path)

def list_systematics_uproot(f):
    return [k.rstrip(";1") for k in f.keys() if "/" not in k]

def list_samples_uproot(f, systematic):
    d = f[systematic]
    return [k.rstrip(";1") for k in d.keys() if "/" not in k]

def load_histograms_uproot(f, systematic, sample):
    base = f[f"{systematic}/{sample}"]
    hists = {}
    for idx in range(5):
        key = f"Variation_{idx}"
        if key not in [k.rstrip(";1") for k in base.keys()]:
            continue
        h = base[key]
        counts, edges = h.to_numpy()
        errors = (np.sqrt(h.variances())
                  if h.variances() is not None
                  else np.sqrt(np.abs(counts)))
        hists[idx] = (edges, counts, errors)
    return hists


# ── ROOT-backend equivalents ─────────────────────────────────────────────────
def load_file_ROOT(path):
    import ROOT
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        sys.exit(f"Cannot open {path}")
    return f

def list_systematics_ROOT(f):
    return [k.GetName() for k in f.GetListOfKeys()
            if k.GetClassName() == "TDirectoryFile"]

def list_samples_ROOT(f, systematic):
    d = f.Get(systematic)
    return [k.GetName() for k in d.GetListOfKeys()
            if k.GetClassName() == "TDirectoryFile"]

def load_histograms_ROOT(f, systematic, sample):
    import ROOT
    base = f.Get(f"{systematic}/{sample}")
    hists = {}
    for idx in range(5):
        h = base.Get(f"Variation_{idx}")
        if not h:
            continue
        n = h.GetNbinsX()
        edges  = np.array([h.GetBinLowEdge(i+1) for i in range(n)]
                          + [h.GetBinLowEdge(n) + h.GetBinWidth(n)])
        counts = np.array([h.GetBinContent(i+1) for i in range(n)])
        errors = np.array([h.GetBinError(i+1)   for i in range(n)])
        hists[idx] = (edges, counts, errors)
    return hists


# ─────────────────────────────────────────────────────────────────────────────
# Unified dispatch
# ─────────────────────────────────────────────────────────────────────────────
def load_file(path):
    return load_file_uproot(path) if _BACKEND == "uproot" else load_file_ROOT(path)

def list_systematics(f):
    return list_systematics_uproot(f) if _BACKEND == "uproot" else list_systematics_ROOT(f)

def list_samples(f, systematic):
    return list_samples_uproot(f, systematic) if _BACKEND == "uproot" else list_samples_ROOT(f, systematic)

def load_histograms(f, systematic, sample):
    return (load_histograms_uproot(f, systematic, sample)
            if _BACKEND == "uproot"
            else load_histograms_ROOT(f, systematic, sample))


# ─────────────────────────────────────────────────────────────────────────────
# Pretty label helpers
# ─────────────────────────────────────────────────────────────────────────────
_SYS_LABELS = {
    "ContainedMuonSqrtEnergyScale":    r"Contained $\mu$ $\sqrt{E}$ Scale",
    "ContainedMuonInvSqrtEnergyScale": r"Contained $\mu$ $1/\sqrt{E}$ Scale",
    "TrackedMuonEnergyScale":          r"Tracked $\mu$ Energy Scale",
    "TrackedMuonSqrtEnergyScale":      r"Tracked $\mu$ $\sqrt{E}$ Scale",
    "TrackedMuonInvSqrtEnergyScale":   r"Tracked $\mu$ $1/\sqrt{E}$ Scale",
    "EMEnergyScale":                   r"EM Energy Scale",
    "EMSqrtEnergyScale":               r"EM $\sqrt{E}$ Scale",
    "EMInvSqrtEnergyScale":            r"EM $1/\sqrt{E}$ Scale",
    "ChgHadEnergyScale":               r"Charged Hadron Energy Scale",
    "ChgHadSqrtEnergyScale":           r"Charged Hadron $\sqrt{E}$ Scale",
    "ChgHadInvSqrtEnergyScale":        r"Charged Hadron $1/\sqrt{E}$ Scale",
    "NeutronEnergyScale":              r"Neutron Energy Scale",
    "NeutronSqrtEnergyScale":          r"Neutron $\sqrt{E}$ Scale",
    "NeutronInvSqrtEnergyScale":       r"Neutron $1/\sqrt{E}$ Scale",
    "TotalEnergyScale":                r"Total Energy Scale",
    "TotalSqrtEnergyScale":            r"Total $\sqrt{E}$ Scale",
    "TotalInvSqrtEnergyScale":         r"Total $1/\sqrt{E}$ Scale",
    "MuonEnergyResolution":            r"Muon Energy Resolution",
    "EMEnergyResolution":              r"EM Energy Resolution",
    "ChgHadEnergyResolution":          r"Charged Hadron Energy Resolution",
    "NeutronEnergyResolution":         r"Neutron Energy Resolution",
}

def sys_label(name):
    return _SYS_LABELS.get(name, name.replace("_", " "))

def sample_label(name):
    return name.replace("_", " ")


# ─────────────────────────────────────────────────────────────────────────────
# Core plotting function  (writes one figure to an open PdfPages)
# ─────────────────────────────────────────────────────────────────────────────
def make_plot(pdf, hists, systematic, sample, x_label="Reconstructed Energy (GeV)"):
    """
    Draw one two-panel figure and save it as the next page in `pdf`.

    Changes vs original:
      • All counts/errors are divided by bin width before plotting.
      • Ratio y-axis auto-scales to the actual data range (+ 5 % padding).
      • Error bands use smooth fill_between at bin centres (no pointy boxes).
    """
    if 2 not in hists:
        print(f"  [WARN] No nominal (Variation_2) for {systematic}/{sample} — skipping.")
        return

    nom_edges, _nom_counts, _nom_errors = hists[2]

    # ── Divide every variation by bin width ───────────────────────────────────
    bw_hists = {}
    for idx, (edges, counts, errors) in hists.items():
        c, e = per_bw(counts, errors, edges)
        bw_hists[idx] = (edges, c, e)

    nom_edges, nom_counts, nom_errors = bw_hists[2]
    bin_centres = 0.5 * (nom_edges[:-1] + nom_edges[1:])
    safe_nom    = np.where(nom_counts > 0, nom_counts, np.nan)

    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1,
        figsize=(8, 7),
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05},
        sharex=True,
    )

    # ── TOP panel ─────────────────────────────────────────────────────────────

    # ±3σ shaded envelope
    if 0 in bw_hists and 4 in bw_hists:
        ax_top.fill_between(bin_centres,
                            bw_hists[0][1], bw_hists[4][1],
                            alpha=BAND_3S_ALPHA, color=BAND_3S_COLOR,
                            label=r"$\pm3\sigma$ envelope", zorder=1)

    # ±1σ shaded envelope
    if 1 in bw_hists and 3 in bw_hists:
        ax_top.fill_between(bin_centres,
                            bw_hists[1][1], bw_hists[3][1],
                            alpha=BAND_1S_ALPHA, color=BAND_1S_COLOR,
                            label=r"$\pm1\sigma$ envelope", zorder=1)

    # Variation lines
    for idx, (edges, counts, errors) in sorted(bw_hists.items()):
        style = VARIATIONS[idx]
        ax_top.step(edges, np.append(counts, counts[-1]),
                    where="post",
                    color=style["color"], ls=style["ls"],
                    lw=style["lw"], label=style["label"],
                    zorder=style["zorder"])

    # Smooth stat error band on nominal
    draw_smooth_errors(ax_top, nom_edges, nom_counts, nom_errors,
                       color=VARIATIONS[2]["color"], alpha=0.30, zorder=6)

    ax_top.set_ylabel("Event Rate / Bin Width (a.u.)", labelpad=8)
    ax_top.set_xlim(nom_edges[0], nom_edges[-1])
    ymax = max(h[1].max() for h in bw_hists.values()) * 1.05
    ax_top.set_ylim(bottom=0, top=ymax)

    handles, labels = ax_top.get_legend_handles_labels()
    ax_top.legend(handles, labels, framealpha=0.85, edgecolor="0.7",
                  ncol=2, loc="upper right", fontsize=9)

    ax_top.set_title(
        f"{sys_label(systematic)}\n"
        r"$\bf{" + sample_label(sample).replace(" ", r"\ ") + r"}$",
        pad=6, fontsize=11,
    )

    # ── BOTTOM panel (ratio) ──────────────────────────────────────────────────

    # Collect all finite ratio values to auto-scale the y-axis
    all_ratios = []

    # Shaded ratio bands
    if 0 in bw_hists and 4 in bw_hists:
        r_lo = bw_hists[0][1] / safe_nom
        r_hi = bw_hists[4][1] / safe_nom
        ax_bot.fill_between(bin_centres, r_lo, r_hi,
                            alpha=BAND_3S_ALPHA, color=BAND_3S_COLOR, zorder=1)
        all_ratios.extend([r_lo, r_hi])

    if 1 in bw_hists and 3 in bw_hists:
        r_lo = bw_hists[1][1] / safe_nom
        r_hi = bw_hists[3][1] / safe_nom
        ax_bot.fill_between(bin_centres, r_lo, r_hi,
                            alpha=BAND_1S_ALPHA, color=BAND_1S_COLOR, zorder=1)
        all_ratios.extend([r_lo, r_hi])

    # Ratio lines
    for idx, (edges, counts, _) in sorted(bw_hists.items()):
        if idx == 2:
            continue
        style = VARIATIONS[idx]
        ratio = counts / safe_nom
        ax_bot.step(edges, np.append(ratio, ratio[-1]),
                    where="post",
                    color=style["color"], ls=style["ls"],
                    lw=style["lw"], zorder=style["zorder"])
        all_ratios.append(ratio)

    # Unity line
    ax_bot.axhline(1.0, color=VARIATIONS[2]["color"], lw=2.0, ls="-", zorder=5)

    # Smooth stat uncertainty on the ratio
    nom_rel_err  = nom_errors / safe_nom
    ratio_ones   = np.ones_like(nom_counts)
    draw_smooth_errors(ax_bot, nom_edges, ratio_ones, nom_rel_err,
                       color=VARIATIONS[2]["color"], alpha=0.22, zorder=4)
    all_ratios.append(ratio_ones + nom_rel_err)
    all_ratios.append(ratio_ones - nom_rel_err)

    # ── Auto-scale ratio y-axis ───────────────────────────────────────────────
    finite_vals = np.concatenate([r[np.isfinite(r)] for r in all_ratios if len(r)])
    if finite_vals.size:
        r_min, r_max = finite_vals.min(), finite_vals.max()
        pad = max(0.05 * (r_max - r_min), 0.02)   # at least 2 % padding
        ax_bot.set_ylim(r_min - pad, r_max + pad)
    else:
        ax_bot.set_ylim(0.7, 1.3)                  # sensible fallback

    ax_bot.yaxis.set_major_locator(ticker.MaxNLocator(nbins=5, symmetric=True))
    ax_bot.yaxis.set_minor_locator(ticker.AutoMinorLocator())

    ax_bot.set_ylabel("Ratio to\nNominal", labelpad=8)
    ax_bot.set_xlabel(x_label, labelpad=8)

    # Experiment watermark — edit or remove as needed
    fig.text(0.13, 0.91, "DUNE ND", fontsize=9,
             fontstyle="italic", color="0.45", transform=fig.transFigure)

    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(
        description="Plot sigma variations — all pages in one multipage PDF."
    )
    p.add_argument("rootfile",
                   nargs="?",
                   default="/scratch/abipeake/MaCh3DUNE_LukesVersion/MaCh3_DUNE/ehadressigmavar.root",
                   help="Path to the .root file (default: JUNE25th_ssigmavar.root)")
    p.add_argument("--systematics", nargs="+", default=None, metavar="SYS",
                   help="Only plot these systematics (default: all)")
    p.add_argument("--samples", nargs="+", default=None, metavar="SAMPLE",
                   help="Only plot these samples (default: all)")
    p.add_argument("--out", default="sigma_variationsehadres.pdf",
                   help="Output PDF filename (default: sigma_variations.pdf)")
    p.add_argument("--xlabel", default="Reconstructed Neutrino Energy (GeV)",
                   help="x-axis label")
    return p.parse_args()


def main():
    args = parse_args()
    apply_thesis_style()

    rootfile = Path(args.rootfile)
    if not rootfile.exists():
        sys.exit(f"ERROR: file not found: {rootfile}")

    print(f"Backend : {_BACKEND}")
    print(f"File    : {rootfile}")
    print(f"Output  : {args.out}")

    f = load_file(str(rootfile))
    systematics = list_systematics(f)

    if args.systematics:
        missing = set(args.systematics) - set(systematics)
        if missing:
            print(f"[WARN] Systematics not found in file: {missing}")
        systematics = [s for s in args.systematics if s in systematics]

    print(f"Systematics ({len(systematics)}): {systematics}\n")

    n_pages = 0
    with PdfPages(args.out) as pdf:
        for systematic in systematics:
            samples = list_samples(f, systematic)
            if args.samples:
                samples = [s for s in args.samples if s in samples]

            for sample in samples:
                print(f"  {systematic} / {sample}")
                hists = load_histograms(f, systematic, sample)
                if not hists:
                    print("    [WARN] No histograms — skipping.")
                    continue
                make_plot(pdf, hists, systematic, sample, x_label=args.xlabel)
                n_pages += 1

    print(f"\nDone. {n_pages} pages → {args.out}")


if __name__ == "__main__":
    main()


