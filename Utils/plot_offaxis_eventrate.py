#!/usr/bin/env python3
"""
plot_offaxis_eventrates.py
===========================

Produces thesis-quality event-rate plots from `eventrate_offaxis.root`.

For every histogram matching `hOffAxis<X>m_numuCC_numode` in the input
ROOT file, this script produces:

  1. An individual, colourful plot of that sample, sized and styled to
     drop straight into a LaTeX `subfigure`/`subcaptionbox` grid.
  2. One combined plot with every off-axis sample overlaid on the same
     axes (never stacked/summed), colour-coded by off-axis distance with
     a colourbar, for use as a single full-width thesis figure.

Each histogram is bin-width-divided (so the y-axis is a rate density,
e.g. events/GeV, independent of binning) and then normalised so that all
samples sit on the same scale and their *shapes* can be compared directly
-- by default each is normalised to unit area (--normalize area), but
you can also normalise to peak height (--normalize max) or turn
normalisation off entirely (--normalize none) to see absolute rates.

Dependencies
------------
    pip install uproot numpy matplotlib

Usage
-----
    python plot_offaxis_eventrates.py eventrate_offaxis.root \
        --outdir thesis_plots \
        --xlabel "Reconstructed Neutrino Energy (GeV)" \
        --normalize area

Notes
-----
- Output is written as both .pdf (vector, for LaTeX \\includegraphics) and
  .png (for quick previews/slides) at 300 dpi.
- If a system LaTeX installation is available, set --usetex to render
  labels with true LaTeX (matches thesis fonts exactly). This is off by
  default since compute-cluster nodes often lack a LaTeX install.
- Adjust --xlabel to match the actual axis quantity/units of your
  histograms (inferred here as reconstructed neutrino energy, but you
  know your file best). --ylabel is set automatically based on
  --normalize unless you override it.
"""

import argparse
import os
import re
import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

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
    and then place all samples on a common scale.

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


def colourful_palette(n, cmap_name="turbo"):
    """n evenly-spaced, vivid colours from a perceptually broad colormap."""
    cmap = mpl.cm.get_cmap(cmap_name)
    if n == 1:
        return [cmap(0.5)]
    # avoid the very dark/near-black ends of turbo/jet-like maps
    return [cmap(0.08 + 0.87 * i / (n - 1)) for i in range(n)]


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("infile", help="Path to eventrate_offaxis.root")
    parser.add_argument("--outdir", default="thesis_plots",
                         help="Directory to write plots to (default: thesis_plots)")
    parser.add_argument("--xlabel", default="Reconstructed Neutrino Energy (GeV)",
                         help="X-axis label (edit to match your histogram's actual quantity)")
    parser.add_argument("--ylabel", default=None,
                         help="Y-axis label (default: chosen automatically from --normalize)")
    parser.add_argument("--normalize", choices=["area", "max", "none"], default="area",
                         help="How to place samples on a common scale after dividing out "
                              "bin width: 'area' (unit integral, default), 'max' "
                              "(unit peak height), or 'none' (rate density only)")
    parser.add_argument("--indiv-cmap", default="turbo",
                         help="Colormap used to give each individual plot a distinct, "
                              "vivid colour")
    parser.add_argument("--fill", dest="fill", action="store_true", default=True,
                         help="Fill under the curve in individual plots (default: on)")
    parser.add_argument("--no-fill", dest="fill", action="store_false",
                         help="Disable the fill under the curve in individual plots")
    parser.add_argument("--pattern", default=r"hOffAxis[\d\.]+m_numuCC_numode",
                         help="Regex used to select histograms from the ROOT file")
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
    parser.add_argument("--cmap", default="viridis",
                         help="Colormap used to colour-code samples by off-axis distance")
    args = parser.parse_args()
    if args.ylabel is None:
        args.ylabel = default_ylabel(args.normalize)

    set_publication_style(usetex=args.usetex)
    os.makedirs(args.outdir, exist_ok=True)
    indiv_dir = os.path.join(args.outdir, "individual")
    os.makedirs(indiv_dir, exist_ok=True)

    # ---- Load histograms -------------------------------------------------
    f = uproot.open(args.infile)
    pattern = re.compile(args.pattern)

    samples = []  # list of dicts: {distance, label, edges, values, errors}
    for key in f.keys(cycle=False):
        if not pattern.fullmatch(key):
            continue
        dist = parse_offaxis_distance(key)
        if dist is None:
            continue
        hist = f[key]
        edges, values, errors = get_hist_arrays(hist)
        samples.append({
            "distance": dist,
            "distance_str": NAME_RE.search(key).group(1),
            "key": key,
            "edges": edges,
            "values": values,
            "errors": errors,
        })

    if not samples:
        sys.exit(f"No histograms matching pattern '{args.pattern}' found in {args.infile}")

    samples.sort(key=lambda s: s["distance"])
    print(f"Found {len(samples)} off-axis samples:")
    for s in samples:
        print(f"  {s['distance']:>6.2f} m   ({s['key']})")

    # ---- Bin-width divide + normalise every sample onto a common scale -----
    for s in samples:
        density, density_err = bin_width_divide_and_normalize(
            s["values"], s["edges"], s["errors"], mode=args.normalize
        )
        s["density"] = density
        s["density_err"] = density_err

    # ---- Individual plots (colourful) --------------------------------------
    colours = colourful_palette(len(samples), args.indiv_cmap)
    for s, colour in zip(samples, colours):
        fig, ax = plt.subplots(figsize=(args.figwidth, args.figheight))
        step_plot(ax, s["edges"], s["density"], s["density_err"],
                  color=colour, fill=args.fill, fill_alpha=0.4)
        ax.set_xlabel(args.xlabel)
        ax.set_ylabel(args.ylabel)
        ax.set_title(f"Off-axis position: {s['distance']:g} m", pad=8,
                      color="#333333", fontweight="bold")
        ax.set_xlim(s["edges"][0], s["edges"][-1])
        ax.set_ylim(bottom=0 if not args.logy else None)
        if args.logy:
            ax.set_yscale("log")
        ax.grid(True, alpha=0.2, linewidth=0.6)
        ax.set_axisbelow(True)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
        fig.tight_layout()

        stem = f"offaxis_{sanitize(s['distance_str'])}m"
        fig.savefig(os.path.join(indiv_dir, f"{stem}.pdf"))
        fig.savefig(os.path.join(indiv_dir, f"{stem}.png"))
        plt.close(fig)

    print(f"\nSaved {len(samples)} individual subfigure-ready plots to '{indiv_dir}/'")

    # ---- Combined plot: every sample overlaid, never summed ----------------
    fig, ax = plt.subplots(figsize=(args.combined_width, args.combined_height))

    dmin = min(s["distance"] for s in samples)
    dmax = max(s["distance"] for s in samples)
    norm = mpl.colors.Normalize(vmin=dmin, vmax=dmax)
    cmap = mpl.cm.get_cmap(args.cmap)

    for s in samples:
        color = cmap(norm(s["distance"]))
        step_plot(ax, s["edges"], s["density"], errors=None, color=color,
                  label=f"{s['distance']:g} m")

    ax.set_xlabel(args.xlabel)
    ax.set_ylabel(args.ylabel)
    ax.set_xlim(0,5)
    #ax.set_xlim(samples[0]["edges"][0], samples[0]["edges"][-1])
    ax.set_ylim(bottom=0 if not args.logy else None)
    if args.logy:
        ax.set_yscale("log")

    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label("Off-axis position (m)")

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
    print(r"""    \caption{Predicted $\nu_\mu$ CC event rate for each off-axis position sampled by DUNE-PRISM, in neutrino-beam mode.}
    \label{fig:offaxis-eventrates}
\end{figure}""")


if __name__ == "__main__":
    main()