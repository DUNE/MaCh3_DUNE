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

Output modes
------------
  Default  – single multipage PDF (one page per systematic × sample)
  --individual – one PDF/PNG per plot, named <Systematic>_<Sample>.<ext>
                 suitable for \includegraphics in LaTeX subfigures

Each page/file has two panels:
  TOP    – absolute histograms (nominal + ±1σ + ±3σ lines)
           with smooth stat error band on the nominal
  BOTTOM – ratio to nominal with shaded ±1σ / ±3σ bands,
           smooth stat uncertainty on the ratio,
           and full ratio-line stat propagation
           (y-axis auto-scales; override with --ratio-range)

Usage
-----
  python plot_sigma_variations.py                             # all systematics, all samples
  python plot_sigma_variations.py --systematics EMEnergyScale TrackedMuonEnergyScale
  python plot_sigma_variations.py --samples OffAxis0m_numuCC_numode
  python plot_sigma_variations.py --out my_results.pdf
  python plot_sigma_variations.py --individual --out-dir plots/
  python plot_sigma_variations.py --figwidth 3.5             # half-page column (inches)
  python plot_sigma_variations.py --xlabel "Muon Momentum (GeV/c)"
  python plot_sigma_variations.py --ratio-range 0.8 1.2
  python plot_sigma_variations.py --watermark "DUNE ND Preliminary"
  python plot_sigma_variations.py --no-title
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
    2: dict(label="Nominal",     color="#1f77b4", ls="-",   lw=2.2, zorder=6),  # drawn last, thicker
    3: dict(label=r"$+1\sigma$", color="#2ca02c", ls="-.",  lw=1.4, zorder=3),
    4: dict(label=r"$+3\sigma$", color="#9467bd", ls="--",  lw=1.4, zorder=2),
}

BAND_1S_COLOR = "#2ca02c"
BAND_3S_COLOR = "#9467bd"
BAND_1S_ALPHA = 0.18
BAND_3S_ALPHA = 0.10

# Clip ratio lines at this many sigma to prevent outlier bins blowing the scale
RATIO_CLIP = 5.0


# ─────────────────────────────────────────────────────────────────────────────
# Publication style — sizes tuned for thesis subfigures
# ─────────────────────────────────────────────────────────────────────────────
def apply_thesis_style(base_fontsize=9):
    plt.rcParams.update({
        "font.family":         "serif",
        "font.serif":          ["Computer Modern Roman", "DejaVu Serif"],
        "font.size":           base_fontsize,
        "axes.titlesize":      base_fontsize + 1,
        "axes.labelsize":      base_fontsize + 1,
        "xtick.labelsize":     base_fontsize - 1,
        "ytick.labelsize":     base_fontsize - 1,
        "legend.fontsize":     base_fontsize - 1,
        "figure.dpi":          150,
        "savefig.dpi":         300,
        "savefig.bbox":        "tight",
        "axes.linewidth":      1.0,
        "xtick.major.width":   1.0,
        "ytick.major.width":   1.0,
        "xtick.minor.width":   0.7,
        "ytick.minor.width":   0.7,
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

def safe_filename(systematic, sample):
    """Sanitise names for use as a filename stem."""
    s = f"{systematic}__{sample}"
    return s.replace("/", "_").replace(" ", "_")


# ─────────────────────────────────────────────────────────────────────────────
# Core plotting function
# ─────────────────────────────────────────────────────────────────────────────
def make_figure(hists, systematic, sample,
                x_label="Reconstructed Neutrino Energy (GeV)",
                figwidth=3.5,
                ratio_range=None,
                show_title=True,
                watermark=None):
    """
    Build and return a matplotlib Figure for one systematic × sample.

    Parameters
    ----------
    hists       : dict {variation_index: (edges, counts, errors)}
    figwidth    : float  — physical width in inches (match your LaTeX textwidth)
    ratio_range : (float, float) or None  — manual y-limits for ratio panel
    show_title  : bool   — set False when LaTeX captions replace the title
    watermark   : str or None — text drawn in bottom-left of the figure
    """
    if 2 not in hists:
        print(f"  [WARN] No nominal (Variation_2) for {systematic}/{sample} — skipping.")
        return None

    # ── Divide every variation by bin width ───────────────────────────────────
    bw_hists = {}
    for idx, (edges, counts, errors) in hists.items():
        c, e = per_bw(counts, errors, edges)
        bw_hists[idx] = (edges, c, e)

    nom_edges, nom_counts, nom_errors = bw_hists[2]
    bin_centres = 0.5 * (nom_edges[:-1] + nom_edges[1:])
    safe_nom    = np.where(nom_counts > 0, nom_counts, np.nan)

    # Aspect ratio: taller than wide suits a two-panel layout in a column
    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1,
        figsize=(figwidth, figwidth * 1.35),
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

    # Non-nominal variation lines (drawn first so nominal sits on top)
    for idx in [0, 1, 3, 4]:
        if idx not in bw_hists:
            continue
        edges, counts, _ = bw_hists[idx]
        style = VARIATIONS[idx]
        ax_top.step(edges, np.append(counts, counts[-1]),
                    where="post",
                    color=style["color"], ls=style["ls"],
                    lw=style["lw"], label=style["label"],
                    zorder=style["zorder"])

    # Smooth stat error band on nominal (drawn before the nominal line)
    draw_smooth_errors(ax_top, nom_edges, nom_counts, nom_errors,
                       color=VARIATIONS[2]["color"], alpha=0.30, zorder=5)

    # Nominal line — drawn last so it sits on top of everything
    ax_top.step(nom_edges, np.append(nom_counts, nom_counts[-1]),
                where="post",
                color=VARIATIONS[2]["color"],
                ls=VARIATIONS[2]["ls"],
                lw=VARIATIONS[2]["lw"],
                label=VARIATIONS[2]["label"],
                zorder=VARIATIONS[2]["zorder"])

    ax_top.set_ylabel("Event Rate / Bin Width (a.u.)", labelpad=6)
    ax_top.set_xlim(nom_edges[0], nom_edges[-1])
    ymax = max(h[1][np.isfinite(h[1])].max() for h in bw_hists.values()) * 1.05
    ax_top.set_ylim(bottom=0, top=ymax)

    # Legend: single column for narrow figures to avoid overcrowding
    handles, labels = ax_top.get_legend_handles_labels()
    ax_top.legend(handles, labels, framealpha=0.85, edgecolor="0.7",
                  ncol=1, loc="upper right")

    if show_title:
        ax_top.set_title(
            f"{sys_label(systematic)} — {sample_label(sample)}",
            pad=4,
        )

    # ── BOTTOM panel (ratio) ──────────────────────────────────────────────────
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

    # Non-nominal ratio lines with propagated stat uncertainty bands
    for idx in [0, 1, 3, 4]:
        if idx not in bw_hists:
            continue
        edges, counts, errors = bw_hists[idx]
        style = VARIATIONS[idx]
        ratio = counts / safe_nom

        # Propagate both nominal and variation statistical errors into ratio
        # σ(ratio) = ratio * sqrt((σ_var/var)² + (σ_nom/nom)²)
        with np.errstate(invalid="ignore", divide="ignore"):
            rel_err_var = np.where(counts > 0, errors / counts, np.nan)
            rel_err_nom = np.where(nom_counts > 0, nom_errors / nom_counts, np.nan)
            ratio_err   = np.abs(ratio) * np.sqrt(rel_err_var**2 + rel_err_nom**2)

        # Clip extreme ratio values from outlier bins
        ratio_clipped = np.where(np.abs(ratio - 1.0) < RATIO_CLIP, ratio, np.nan)

        ax_bot.step(edges, np.append(ratio_clipped, ratio_clipped[-1]),
                    where="post",
                    color=style["color"], ls=style["ls"],
                    lw=style["lw"], zorder=style["zorder"])

        draw_smooth_errors(ax_bot, edges, ratio_clipped, ratio_err,
                           color=style["color"], alpha=0.15, zorder=2)

        all_ratios.append(ratio_clipped)

    # Unity line (nominal reference)
    ax_bot.axhline(1.0, color=VARIATIONS[2]["color"], lw=2.0, ls="-", zorder=5)

    # Smooth stat uncertainty band on nominal ratio (pure nominal stat error)
    nom_rel_err = nom_errors / safe_nom
    ratio_ones  = np.ones_like(nom_counts)
    draw_smooth_errors(ax_bot, nom_edges, ratio_ones, nom_rel_err,
                       color=VARIATIONS[2]["color"], alpha=0.22, zorder=4)
    all_ratios.append(ratio_ones + nom_rel_err)
    all_ratios.append(ratio_ones - nom_rel_err)

    # ── Ratio y-axis ──────────────────────────────────────────────────────────
    if ratio_range is not None:
        ax_bot.set_ylim(*ratio_range)
    else:
        finite_vals = np.concatenate([r[np.isfinite(r)] for r in all_ratios if len(r)])
        if finite_vals.size:
            r_min, r_max = finite_vals.min(), finite_vals.max()
            pad = max(0.05 * (r_max - r_min), 0.02)
            ax_bot.set_ylim(r_min - pad, r_max + pad)
        else:
            ax_bot.set_ylim(0.7, 1.3)

    ax_bot.yaxis.set_major_locator(ticker.MaxNLocator(nbins=4, symmetric=True))
    ax_bot.yaxis.set_minor_locator(ticker.AutoMinorLocator())

    ax_bot.set_ylabel("Ratio to\nNominal", labelpad=6)
    ax_bot.set_xlabel(x_label, labelpad=6)

    # ── Watermark (axis-relative so it survives any figsize) ─────────────────
    if watermark:
        ax_top.text(0.02, 0.97, watermark,
                    transform=ax_top.transAxes,
                    fontsize=plt.rcParams["font.size"] - 1,
                    fontstyle="italic", color="0.45",
                    va="top", ha="left")

    return fig


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(
        description="Plot sigma variations — multipage PDF or individual files."
    )
    p.add_argument("rootfile",
                   nargs="?",
                   default="JULY2nd_sigmavar.root",
                   help="Path to the .root file (default: JUNE25th_ssigmavar.root)")
    p.add_argument("--systematics", nargs="+", default=None, metavar="SYS",
                   help="Only plot these systematics (default: all)")
    p.add_argument("--samples", nargs="+", default=None, metavar="SAMPLE",
                   help="Only plot these samples (default: all)")

    # Output mode
    mode = p.add_mutually_exclusive_group()
    mode.add_argument("--out", default="sigma_variations.pdf",
                      help="Output multipage PDF filename (default: sigma_variations.pdf)")
    mode.add_argument("--individual", action="store_true",
                      help="Save one file per plot instead of a multipage PDF")
    p.add_argument("--out-dir", default="plots",
                   help="Output directory for individual files (default: plots/)")
    p.add_argument("--fmt", default="pdf", choices=["pdf", "png", "svg"],
                   help="File format for individual mode (default: pdf)")

    # Figure geometry
    p.add_argument("--figwidth", type=float, default=3.5,
                   help="Figure width in inches — match your LaTeX column width "
                        "(default: 3.5 for a half-page subfigure, use 6.5 for full-width)")

    # Axis / style options
    p.add_argument("--xlabel", default="Reconstructed Neutrino Energy (GeV)",
                   help="x-axis label")
    p.add_argument("--ratio-range", nargs=2, type=float, default=None,
                   metavar=("LO", "HI"),
                   help="Fix ratio y-axis range, e.g. --ratio-range 0.8 1.2")
    p.add_argument("--no-title", action="store_true",
                   help="Omit the figure title (use LaTeX captions instead)")
    p.add_argument("--watermark", default="DUNE ND",
                   help="Watermark text drawn on the top panel (default: 'DUNE ND'); "
                        "pass empty string to disable")

    return p.parse_args()


def main():
    args = parse_args()
    apply_thesis_style(base_fontsize=9 if args.figwidth <= 4.0 else 11)

    rootfile = Path(args.rootfile)
    if not rootfile.exists():
        sys.exit(f"ERROR: file not found: {rootfile}")

    watermark   = args.watermark if args.watermark else None
    show_title  = not args.no_title
    ratio_range = tuple(args.ratio_range) if args.ratio_range else None

    print(f"Backend  : {_BACKEND}")
    print(f"File     : {rootfile}")
    print(f"Figwidth : {args.figwidth} in")

    f = load_file(str(rootfile))
    systematics = list_systematics(f)

    if args.systematics:
        missing = set(args.systematics) - set(systematics)
        if missing:
            print(f"[WARN] Systematics not found in file: {missing}")
        systematics = [s for s in args.systematics if s in systematics]

    print(f"Systematics ({len(systematics)}): {systematics}\n")

    # ── Individual-file mode ──────────────────────────────────────────────────
    if args.individual:
        out_dir = Path(args.out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        print(f"Output   : individual {args.fmt.upper()} files → {out_dir}/\n")
        n_files = 0
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
                fig = make_figure(
                    hists, systematic, sample,
                    x_label=args.xlabel,
                    figwidth=args.figwidth,
                    ratio_range=ratio_range,
                    show_title=show_title,
                    watermark=watermark,
                )
                if fig is None:
                    continue
                stem = safe_filename(systematic, sample)
                fpath = out_dir / f"{stem}.{args.fmt}"
                fig.savefig(fpath, bbox_inches="tight")
                plt.close(fig)
                n_files += 1
        print(f"\nDone. {n_files} files → {out_dir}/")

    # ── Multipage PDF mode (default) ─────────────────────────────────────────
    else:
        print(f"Output   : {args.out}\n")
        n_pages = 0
        metadata = {
            "Title":   "Sigma Variations",
            "Author":  "",          # fill in your name
            "Subject": "Systematic uncertainties",
            "Creator": "plot_sigma_variations.py",
        }
        with PdfPages(args.out, metadata=metadata) as pdf:
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
                    fig = make_figure(
                        hists, systematic, sample,
                        x_label=args.xlabel,
                        figwidth=args.figwidth,
                        ratio_range=ratio_range,
                        show_title=show_title,
                        watermark=watermark,
                    )
                    if fig is None:
                        continue
                    pdf.savefig(fig, bbox_inches="tight")
                    plt.close(fig)
                    n_pages += 1
        print(f"\nDone. {n_pages} pages → {args.out}")


if __name__ == "__main__":
    main()