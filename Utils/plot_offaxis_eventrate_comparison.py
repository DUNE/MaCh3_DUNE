#!/usr/bin/env python3
"""
plot_offaxis_eventrate_comparison.py
===========================

Produces thesis-quality event-rate plots from a pair of ROOT files:
one containing reconstructed-energy histograms, one containing
true-energy histograms, for the same set of DUNE-PRISM off-axis
samples.

For every off-axis distance found (matched across both files by the
histogram name `hOffAxis<X>m_numuCC_numode`), this script produces:

  1. An individual, two-curve plot (reconstructed energy vs. true
     energy, overlaid) for that sample, sized and styled to drop
     straight into a LaTeX `subfigure`/`subcaptionbox` grid.
  2. One combined plot with every off-axis sample's reco and true
     curves overlaid on the same axes (never stacked/summed).
     Colour encodes off-axis distance (colourbar); linestyle encodes
     reconstructed (solid) vs. true (dashed) energy.

Each histogram is bin-width-divided (so the y-axis is a rate density,
e.g. events/GeV, independent of binning) and then normalised so that
samples sit on a comparable scale and their *shapes* can be compared
directly -- by default each histogram is normalised to unit area
(--normalize area) independently of its reco/true partner, but you can
also normalise to peak height (--normalize max) or turn normalisation
off entirely (--normalize none) to see absolute rates.

Dependencies
------------
    pip install uproot numpy matplotlib

Usage
-----
    python plot_offaxis_eventrate_comparison \\
        eventrate_offaxis_reco.root eventrate_offaxis_true.root \\
        --outdir thesis_plots \\
        --xlabel "Neutrino Energy (GeV)" \\
        --normalize area

Notes
-----
- Output is written as both .pdf (vector, for LaTeX \\includegraphics) and
  .png (for quick previews/slides) at 300 dpi.
- If a system LaTeX installation is available, set --usetex to render
  labels with true LaTeX (matches thesis fonts exactly). This is off by
  default since compute-cluster nodes often lack a LaTeX install.
- Adjust --xlabel to match the actual axis quantity/units of your
  histograms. --ylabel is set automatically based on --normalize unless
  you override it.
- If a given off-axis distance is present in only one of the two files,
  it is skipped (with a warning) since there is nothing to overlay it
  against.
"""

import argparse
import os
import re
import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

try:
    import uproot
except ImportError:
    sys.exit(
        "This script requires uproot. Install it with:\n"
        "    pip install uproot\n"
    )


# --------------------------------------------------------------------------
# Styling
# --------------------------------------------------------------------------
def set_publication_style(usetex: bool = False):
    """Configure matplotlib rcParams for a clean, thesis-ready look."""
    mpl.rcParams.update({
        "text.usetex": usetex,
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman", "Times New Roman", "DejaVu Serif"],
        "mathtext.fontset": "cm",
        "font.size": 11,
        "axes.labelsize": 12,
        "axes.titlesize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 9,
        "legend.frameon": False,
        "axes.linewidth": 1.0,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.major.size": 5,
        "xtick.minor.size": 3,
        "ytick.major.size": 5,
        "ytick.minor.size": 3,
        "xtick.minor.visible": True,
        "ytick.minor.visible": True,
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.05,
        "axes.grid": False,
        "lines.linewidth": 1.4,
    })


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
NAME_RE = re.compile(r"hOffAxis([\d\.]+)m_numuCC_numode")


def parse_offaxis_distance(key_name: str):
    """Extract the off-axis distance (float, metres) from a histogram key name."""
    m = NAME_RE.search(key_name)
    if not m:
        return None
    return float(m.group(1))


def sanitize(distance_str: str) -> str:
    """Turn e.g. '1.75' into '1p75' for safe filenames."""
    return distance_str.replace(".", "p")


def get_hist_arrays(hist):
    """Return (edges, values, errors) as numpy arrays from an uproot TH1."""
    values, edges = hist.to_numpy()
    errors = hist.errors() if hasattr(hist, "errors") else np.sqrt(np.abs(values))
    return edges, values, errors


def step_plot(ax, edges, values, errors=None, fill=False, fill_alpha=0.35, **kwargs):
    """Draw a ROOT-style step histogram (with optional error bars and fill)."""
    centers = 0.5 * (edges[:-1] + edges[1:])
    widths = np.diff(edges)

    line_kwargs = dict(kwargs)
    color = line_kwargs.pop("color", None)

    ax.stairs(values, edges, color=color, fill=False, **line_kwargs)
    if fill:
        ax.stairs(values, edges, color=color, fill=True, alpha=fill_alpha, baseline=0)
    if errors is not None:
        ax.errorbar(
            centers, values, yerr=errors, xerr=widths / 2.0,
            fmt="none", ecolor=color, elinewidth=0.8, capsize=0, alpha=0.9,
        )
    return centers, widths


def bin_width_divide_and_normalize(values, edges, errors, mode="area"):
    """
    Divide out the bin width (values -> a rate density, e.g. events/GeV)
    and then place the histogram on a common scale.

    mode:
      "area" -- normalise so the histogram integrates to 1 (shape comparison)
      "max"  -- normalise so the peak bin is 1
      "none" -- only divide by bin width, no further normalisation
    """
    widths = np.diff(edges)
    density = values / widths
    density_err = errors / widths

    if mode == "area":
        integral = np.sum(values)  # raw (pre-density) event count
        if integral > 0:
            density = density / integral
            density_err = density_err / integral
    elif mode == "max":
        peak = density.max() if density.size else 0
        if peak > 0:
            density = density / peak
            density_err = density_err / peak
    elif mode == "none":
        pass
    else:
        raise ValueError(f"Unknown normalisation mode: {mode}")

    return density, density_err


def default_ylabel(mode):
    return {
        "area": "Normalised events / GeV (area = 1)",
        "max": "Events / GeV (normalised to peak)",
        "none": "Events / GeV",
    }[mode]


def load_samples(infile, pattern):
    """Load every histogram matching `pattern` from `infile` into a dict
    keyed by distance_str, each entry holding edges/values/errors."""
    f = uproot.open(infile)
    samples = {}
    for key in f.keys(cycle=False):
        if not pattern.fullmatch(key):
            continue
        dist = parse_offaxis_distance(key)
        if dist is None:
            continue
        hist = f[key]
        edges, values, errors = get_hist_arrays(hist)
        distance_str = NAME_RE.search(key).group(1)
        samples[distance_str] = {
            "distance": dist,
            "distance_str": distance_str,
            "key": key,
            "edges": edges,
            "values": values,
            "errors": errors,
        }
    return samples


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("reco_file", help="Path to the ROOT file with reconstructed-energy histograms")
    parser.add_argument("true_file", help="Path to the ROOT file with true-energy histograms")
    parser.add_argument("--outdir", default="thesis_plots",
                         help="Directory to write plots to (default: thesis_plots)")
    parser.add_argument("--xlabel", default="Neutrino Energy (GeV)",
                         help="X-axis label (edit to match your histograms' actual quantity)")
    parser.add_argument("--ylabel", default=None,
                         help="Y-axis label (default: chosen automatically from --normalize)")
    parser.add_argument("--normalize", choices=["area", "max", "none"], default="area",
                         help="How to place each histogram on a common scale after dividing "
                              "out bin width: 'area' (unit integral, default), 'max' "
                              "(unit peak height), or 'none' (rate density only). Reco and "
                              "true histograms are normalised independently of each other.")
    parser.add_argument("--reco-label", default="Reconstructed energy",
                         help="Legend label for the reconstructed-energy curve")
    parser.add_argument("--true-label", default="True energy",
                         help="Legend label for the true-energy curve")
    parser.add_argument("--reco-color", default="#1f77b4",
                         help="Fixed colour for the reconstructed-energy curve in individual plots")
    parser.add_argument("--true-color", default="#d62728",
                         help="Fixed colour for the true-energy curve in individual plots")
    parser.add_argument("--reco-linestyle", default="-",
                         help="Line style for reconstructed-energy curves (used in both "
                              "individual and combined plots)")
    parser.add_argument("--true-linestyle", default="--",
                         help="Line style for true-energy curves (used in both individual "
                              "and combined plots)")
    parser.add_argument("--fill", dest="fill", action="store_true", default=False,
                         help="Fill under the curves in individual plots (default: off, "
                              "since two overlapping fills get busy)")
    parser.add_argument("--no-fill", dest="fill", action="store_false",
                         help="Disable the fill under the curves in individual plots")
    parser.add_argument("--pattern", default=r"hOffAxis[\d\.]+m_numuCC_numode",
                         help="Regex used to select histograms from the ROOT files")
    parser.add_argument("--figwidth", type=float, default=4.2,
                         help="Width (inches) of each individual subfigure plot")
    parser.add_argument("--figheight", type=float, default=3.2,
                         help="Height (inches) of each individual subfigure plot")
    parser.add_argument("--combined-width", type=float, default=7.0,
                         help="Width (inches) of the combined superposed plot")
    parser.add_argument("--combined-height", type=float, default=5.2,
                         help="Height (inches) of the combined superposed plot")
    parser.add_argument("--logy", action="store_true",
                         help="Use a log scale on the y-axis")
    parser.add_argument("--usetex", action="store_true",
                         help="Use a real LaTeX installation to render text (requires latex on PATH)")
    parser.add_argument("--cmap", default="cividis",
                         help="Colormap used to colour-code samples by off-axis distance in "
                              "the combined plot. 'cividis' (default) is designed specifically "
                              "to be perceptually uniform for all common forms of colour "
                              "vision deficiency; 'viridis' is a good alternative if you want "
                              "more contrast between colours")
    parser.add_argument("--combined-linewidth", type=float, default=1.5,
                         help="Line width for curves in the combined plot")
    parser.add_argument("--combined-alpha", type=float, default=0.9,
                         help="Line transparency in the combined plot (helps show overlap)")
    parser.add_argument("--inset", dest="inset", action="store_true", default=True,
                         help="Add a zoomed inset over the crowded low-energy peak in the "
                              "combined plot (default: on)")
    parser.add_argument("--no-inset", dest="inset", action="store_false",
                         help="Disable the zoomed inset")
    parser.add_argument("--inset-xmin", type=float, default=0.0,
                         help="Lower x-axis bound (GeV) of the inset zoom region")
    parser.add_argument("--inset-xmax", type=float, default=2.0,
                         help="Upper x-axis bound (GeV) of the inset zoom region")
    parser.add_argument("--inset-loc", default="upper right",
                         help="Location of the inset within the combined plot "
                              "(matplotlib legend-style location string)")
    parser.add_argument("--inset-size", default="48%",
                         help="Width/height of the inset as a percentage of the main axes")
    args = parser.parse_args()
    if args.ylabel is None:
        args.ylabel = default_ylabel(args.normalize)

    set_publication_style(usetex=args.usetex)
    os.makedirs(args.outdir, exist_ok=True)
    indiv_dir = os.path.join(args.outdir, "individual")
    os.makedirs(indiv_dir, exist_ok=True)

    # ---- Load histograms from both files -----------------------------------
    pattern = re.compile(args.pattern)
    reco_samples = load_samples(args.reco_file, pattern)
    true_samples = load_samples(args.true_file, pattern)

    reco_keys = set(reco_samples.keys())
    true_keys = set(true_samples.keys())
    common_keys = reco_keys & true_keys
    only_reco = reco_keys - true_keys
    only_true = true_keys - reco_keys

    if only_reco:
        print(f"Warning: distances only found in reco file, skipping: {sorted(only_reco)}")
    if only_true:
        print(f"Warning: distances only found in true file, skipping: {sorted(only_true)}")
    if not common_keys:
        sys.exit("No off-axis distances matched between the two files -- nothing to plot.")

    samples = []
    for distance_str in common_keys:
        reco = reco_samples[distance_str]
        true = true_samples[distance_str]
        samples.append({
            "distance": reco["distance"],
            "distance_str": distance_str,
            "reco": reco,
            "true": true,
        })
    samples.sort(key=lambda s: s["distance"])

    print(f"Found {len(samples)} matched off-axis samples (reco + true):")
    for s in samples:
        print(f"  {s['distance']:>6.2f} m   (reco: {s['reco']['key']}, true: {s['true']['key']})")

    # ---- Bin-width divide + normalise every histogram (reco and true, ------
    # ---- independently of each other) onto a common scale ------------------
    for s in samples:
        for kind in ("reco", "true"):
            h = s[kind]
            density, density_err = bin_width_divide_and_normalize(
                h["values"], h["edges"], h["errors"], mode=args.normalize
            )
            h["density"] = density
            h["density_err"] = density_err

    # ---- Individual plots: reco vs. true overlaid for each distance --------
    for s in samples:
        fig, ax = plt.subplots(figsize=(args.figwidth, args.figheight))

        reco, true = s["reco"], s["true"]
        step_plot(ax, reco["edges"], reco["density"], reco["density_err"],
                  color=args.reco_color, linestyle=args.reco_linestyle,
                  fill=args.fill, fill_alpha=0.35, label=args.reco_label)
        step_plot(ax, true["edges"], true["density"], true["density_err"],
                  color=args.true_color, linestyle=args.true_linestyle,
                  fill=args.fill, fill_alpha=0.35, label=args.true_label)

        ax.set_xlabel(args.xlabel)
        ax.set_ylabel(args.ylabel)
        ax.set_title(f"Off-axis position: {s['distance']:g} m", pad=8,
                      color="#333333", fontweight="bold")
        xlo = min(reco["edges"][0], true["edges"][0])
        xhi = max(reco["edges"][-1], true["edges"][-1])
        ax.set_xlim(xlo, xhi)
        ax.set_ylim(bottom=0 if not args.logy else None)
        if args.logy:
            ax.set_yscale("log")
        ax.grid(True, alpha=0.2, linewidth=0.6)
        ax.set_axisbelow(True)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
        ax.legend(loc="best")
        fig.tight_layout()

        stem = f"offaxis_{sanitize(s['distance_str'])}m"
        fig.savefig(os.path.join(indiv_dir, f"{stem}.pdf"))
        fig.savefig(os.path.join(indiv_dir, f"{stem}.png"))
        plt.close(fig)

    print(f"\nSaved {len(samples)} individual subfigure-ready plots to '{indiv_dir}/'")

    # ---- Combined plot: every sample's reco + true curves overlaid, --------
    # ---- never summed. Colour = distance, linestyle = reco/true. -----------
    fig, ax = plt.subplots(figsize=(args.combined_width, args.combined_height))

    dmin = min(s["distance"] for s in samples)
    dmax = max(s["distance"] for s in samples)
    norm = mpl.colors.Normalize(vmin=dmin, vmax=dmax)
    cmap = mpl.cm.get_cmap(args.cmap)

    for s in samples:
        color = cmap(norm(s["distance"]))
        reco, true = s["reco"], s["true"]
        step_plot(ax, reco["edges"], reco["density"], errors=None, color=color,
                  linestyle=args.reco_linestyle, linewidth=args.combined_linewidth,
                  alpha=args.combined_alpha)
        step_plot(ax, true["edges"], true["density"], errors=None, color=color,
                  linestyle=args.true_linestyle, linewidth=args.combined_linewidth,
                  alpha=args.combined_alpha)

    ax.set_xlabel(args.xlabel)
    ax.set_ylabel(args.ylabel)
    xlo = min(min(s["reco"]["edges"][0], s["true"]["edges"][0]) for s in samples)
    xhi = max(max(s["reco"]["edges"][-1], s["true"]["edges"][-1]) for s in samples)
    ax.set_xlim(xlo, xhi)
    ax.set_ylim(bottom=0 if not args.logy else None)
    if args.logy:
        ax.set_yscale("log")

    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label("Off-axis position (m)")

    # Linestyle legend (colour is taken by the colourbar, so explain
    # reco vs. true separately via proxy artists in grey).
    style_handles = [
        Line2D([0], [0], color="0.25", linestyle=args.reco_linestyle,
               linewidth=args.combined_linewidth, label=args.reco_label),
        Line2D([0], [0], color="0.25", linestyle=args.true_linestyle,
               linewidth=args.combined_linewidth, label=args.true_label),
    ]
    ax.legend(handles=style_handles, loc="upper left")

    # ---- Zoomed inset over the crowded low-energy peak ---------------------
    if args.inset:
        axins = inset_axes(ax, width=args.inset_size, height=args.inset_size,
                            loc=args.inset_loc, borderpad=1.6)
        for s in samples:
            color = cmap(norm(s["distance"]))
            reco, true = s["reco"], s["true"]
            step_plot(axins, reco["edges"], reco["density"], errors=None, color=color,
                      linestyle=args.reco_linestyle, linewidth=args.combined_linewidth,
                      alpha=args.combined_alpha)
            step_plot(axins, true["edges"], true["density"], errors=None, color=color,
                      linestyle=args.true_linestyle, linewidth=args.combined_linewidth,
                      alpha=args.combined_alpha)

        axins.set_xlim(args.inset_xmin, args.inset_xmax)
        ymax = 0.0
        for s in samples:
            for kind in ("reco", "true"):
                h = s[kind]
                mask = (h["edges"][:-1] >= args.inset_xmin) & (h["edges"][:-1] <= args.inset_xmax)
                candidate = h["density"][mask].max() if mask.any() else h["density"].max()
                ymax = max(ymax, candidate)
        axins.set_ylim(0, ymax * 1.08)
        axins.tick_params(labelsize=7)
        axins.set_facecolor("white")
        for spine in axins.spines.values():
            spine.set_linewidth(0.8)
        mark_inset(ax, axins, loc1=2, loc2=3, fc="none", ec="0.4", lw=0.7, alpha=0.8)

    fig.tight_layout()
    combined_stem = "offaxis_all_superposed"
    fig.savefig(os.path.join(args.outdir, f"{combined_stem}.pdf"))
    fig.savefig(os.path.join(args.outdir, f"{combined_stem}.png"))
    plt.close(fig)

    print(f"Saved combined superposed plot to '{args.outdir}/{combined_stem}.pdf'")
    print("\nDone. LaTeX snippet for the individual plots as a subfigure grid:\n")
    print(r"""\begin{figure}[htbp]
    \centering""")
    for i, s in enumerate(samples):
        stem = f"offaxis_{sanitize(s['distance_str'])}m"
        print(r"    \begin{subfigure}[b]{0.3\textwidth}")
        print(r"        \centering")
        print(rf"        \includegraphics[width=\textwidth]{{individual/{stem}.pdf}}")
        print(rf"        \caption{{{s['distance']:g}\,m off-axis}}")
        print(r"    \end{subfigure}")
        if (i + 1) % 3 == 0 and i != len(samples) - 1:
            print(r"    \\[1em]")
    print(r"""    \caption{Predicted $\nu_\mu$ CC event rate for each off-axis position sampled by DUNE-PRISM, in neutrino-beam mode, comparing reconstructed and true neutrino energy.}
    \label{fig:offaxis-eventrates}
\end{figure}""")


if __name__ == "__main__":
    main()