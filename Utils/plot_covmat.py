"""
plot_covariance.py
------------------
Thesis-quality plots of a covariance (or correlation) matrix from a ROOT or
numpy file.

Produces:
  1. Covariance matrix heatmap
  2. Correlation matrix heatmap  (derived automatically)

Usage:
  # From numpy binary (output of combine_covariance.py):
  python plot_covariance.py --input combined_covmat.npy

  # From a ROOT file (TMatrixT<double> or TH2):
  python plot_covariance.py --input combined_covmat.root --key NDCovMatrix

  # Custom output name / format:
  python plot_covariance.py --input combined_covmat.npy --output cov_plot.pdf

  # Supply your own axis labels:
  python plot_covariance.py --input combined_covmat.npy --labels labels.txt

  # Plot only the correlation matrix:
  python plot_covariance.py --input combined_covmat.npy --mode correlation
"""

import argparse
import os
import sys

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable

import yaml

def labels_from_binning_yaml(path, sample_name):
    with open(path) as f:
        cfg = yaml.safe_load(f)

    edges = cfg[sample_name]["Binning"]["VarBins"][0]

    labels = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        labels.append(f"{lo:.2g}–{hi:.2g}")

    return labels


# ---------------------------------------------------------------------------
# Colour palette — bespoke blue/white/red diverging map
# ---------------------------------------------------------------------------

def make_bwr_cmap() -> LinearSegmentedColormap:
    """
    A clean, perceptually-balanced blue–white–red diverging colormap.
    Deeper blue = strong negative correlation, deep red = strong positive,
    white = zero. Avoids the harshness of the default 'bwr'.
    """
    colors = [
        (0.017, 0.169, 0.490),   # deep navy blue
        (0.196, 0.416, 0.722),   # mid blue
        (0.573, 0.745, 0.894),   # light blue
        (1.000, 1.000, 1.000),   # white (zero)
        (0.957, 0.647, 0.510),   # light red/salmon
        (0.796, 0.196, 0.125),   # mid red
        (0.490, 0.012, 0.012),   # deep crimson
    ]
    cmap = LinearSegmentedColormap.from_list("thesis_bwr", colors, N=512)
    return cmap


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------

def read_matrix(path: str, key: str | None) -> np.ndarray:
    ext = os.path.splitext(path)[-1].lower()

    if ext == ".npy":
        mat = np.load(path)
        print(f"Loaded numpy matrix: shape {mat.shape}")
        return mat

    if ext == ".root":
        try:
            import uproot
        except ImportError:
            sys.exit("uproot is required for ROOT files: pip install uproot")
        with uproot.open(path) as f:
            if key is None:
                sys.exit("--key required when reading from a ROOT file.")
            bare = key.split(";")[0]
            obj = None
            for k in [key, bare, bare + ";1"]:
                if k in f:
                    obj = f[k]
                    break
            if obj is None:
                sys.exit(f"Key '{key}' not found. Available: {list(f.keys())}")
            class_name = obj.classname if hasattr(obj, "classname") else ""
            if "TH2" in class_name:
                mat = obj.values()
            elif "TMatrix" in class_name:
                mat = np.array(obj.member("fElements")).reshape(
                    obj.member("fNrows"), obj.member("fNcols")
                )
            else:
                sys.exit(f"Unsupported ROOT class: {class_name}")
        print(f"Loaded ROOT matrix ({class_name}): shape {mat.shape}")
        return mat

    sys.exit(f"Unsupported file extension '{ext}'. Use .npy or .root.")


def covariance_to_correlation(cov: np.ndarray) -> np.ndarray:
    std = np.sqrt(np.diag(cov))
    with np.errstate(invalid="ignore", divide="ignore"):
        corr = cov / np.outer(std, std)
    corr = np.clip(corr, -1.0, 1.0)
    np.fill_diagonal(corr, 1.0)
    return corr


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

THESIS_STYLE = {
    # Font — use a serif for body, matching common thesis styles
    "font.family":        "serif",
    "font.serif":         ["Times New Roman", "DejaVu Serif", "serif"],
    "mathtext.fontset":   "stix",
    # Sizes
    "font.size":          11,
    "axes.titlesize":     13,
    "axes.labelsize":     12,
    "xtick.labelsize":    9,
    "ytick.labelsize":    9,
    "figure.dpi":         150,
    # Lines / ticks
    "axes.linewidth":     0.8,
    "xtick.major.width":  0.8,
    "ytick.major.width":  0.8,
    "xtick.direction":    "out",
    "ytick.direction":    "out",
    # Grid off by default
    "axes.grid":          False,
    # Layout
    "figure.constrained_layout.use": True,
}


def plot_matrix(
    mat: np.ndarray,
    title: str,
    cbar_label: str,
    labels: list[str] | None,
    output_path: str,
    vmin: float | None = None,
    vmax: float | None = None,
    symmetric_norm: bool = True,
) -> None:
    """Render one matrix as a heatmap and save to file."""

    n = mat.shape[0]
    cmap = make_bwr_cmap()

    # Normalisation: symmetric around zero so white = 0
    abs_max = np.nanmax(np.abs(mat))
    v = vmax if vmax is not None else abs_max
    if symmetric_norm:
        norm = TwoSlopeNorm(vmin=-v, vcenter=0.0, vmax=v)
    else:
        norm = mpl.colors.Normalize(vmin=vmin or mat.min(), vmax=vmax or mat.max())

    # Figure sizing: scale with matrix dimension, but cap for readability
    cell = max(0.45, min(0.75, 12.0 / n))
    fig_size = min(n * cell + 2.5, 16.0)

    with mpl.rc_context(THESIS_STYLE):
        fig, ax = plt.subplots(figsize=(fig_size, fig_size * 0.88))

        # origin="lower" places row 0 at the bottom, so the diagonal runs
        # from bottom-left to top-right (standard mathematical convention).
        im = ax.imshow(mat, cmap=cmap, norm=norm, aspect="equal",
                       interpolation="nearest", origin="lower")

        # --- Axis ticks & labels ---
        if labels:
            ax.set_xticks(range(n))
            ax.set_yticks(range(n))
            ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
            ax.set_yticklabels(labels, fontsize=8)
        else:
            # Show every 5th index
            step = max(1, n // 10)
            ticks = list(range(0, n, step))
            ax.set_xticks(ticks)
            ax.set_yticks(ticks)
            ax.set_xticklabels([str(t) for t in ticks])
            ax.set_yticklabels([str(t) for t in ticks])

        ax.set_xlabel(r"Reconstructed Neutrino Energy [GeV]")
        ax.set_ylabel(r"Reconstructed Neutrino Energy [GeV]")
        ax.set_title(title, pad=10, fontweight="bold")

        # Thin grid lines between cells for small matrices
        if n <= 40:
            ax.set_xticks(np.arange(-0.5, n, 1), minor=True)
            ax.set_yticks(np.arange(-0.5, n, 1), minor=True)
            ax.grid(which="minor", color="white", linewidth=0.4, alpha=0.6)
            ax.tick_params(which="minor", bottom=False, left=False)

        # Diagonal highlight
        for i in range(n):
            ax.add_patch(mpl.patches.Rectangle(
                (i - 0.5, i - 0.5), 1, 1,
                fill=False, edgecolor="black", linewidth=0.5, alpha=0.5
            ))

        # --- Colourbar ---
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="4%", pad=0.12)
        cb = fig.colorbar(im, cax=cax)
        cb.set_label(cbar_label, labelpad=8)
        cb.outline.set_linewidth(0.8)
        cb.ax.tick_params(labelsize=9, width=0.8)

        # Neat formatting for colourbar ticks
        cb.locator = ticker.MaxNLocator(nbins=8, symmetric=True)
        cb.update_ticks()

        fig.savefig(output_path, dpi=300, bbox_inches="tight",
                    facecolor="white", edgecolor="none")
        plt.close(fig)

    print(f"  Saved: {output_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Plot a covariance/correlation matrix to thesis-quality figures."
    )
    parser.add_argument("--input", required=True,
                        help="Input file: .npy or .root")
    parser.add_argument("--key", default=None,
                        help="ROOT key name (required for .root input).")
    parser.add_argument("--output", default=None,
                        help="Output file base (without extension). "
                             "Default: same name as input.")
    parser.add_argument("--format", default="pdf", choices=["pdf", "png", "svg"],
                        help="Output format. Default: pdf (vector, ideal for thesis).")
    parser.add_argument("--labels", default=None,
                        help="Path to a plain-text file with one label per line "
                             "(one per matrix row/column).")
    parser.add_argument("--mode", default="both",
                        choices=["covariance", "correlation", "both"],
                        help="Which plot(s) to produce. Default: both.")
    parser.add_argument("--title-prefix", default="",
                        help="Optional prefix for plot titles.")
    args = parser.parse_args()

    # Load matrix
    cov = read_matrix(args.input, args.key)
    n = cov.shape[0]

    edges = [
    0, 0.5, 0.6556, 0.8111, 0.9667, 1.1222,
    1.2778, 1.4333, 1.5889, 1.7444, 1.9,
    2.0556, 2.2111, 2.3667, 2.5222, 2.6778,
    2.8333, 2.9889, 3.2, 3.4, 3.6, 3.8,
    4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0
    ]

    default_labels = [
            f"{lo:.2g}–{hi:.2g}"
            for lo, hi in zip(edges[:-1], edges[1:])
        ]

    if len(default_labels) == n:
        labels = default_labels
    else:
        print(
            f"WARNING: {len(default_labels)} default energy-bin labels for a "
            f"{n}×{n} matrix — labels don't correspond 1:1 to matrix bins "
            "(likely multiple samples/selections concatenated). "
            "Falling back to index ticks."
        )
        labels = None

    # Override with user-supplied labels if requested
    if args.labels:
        with open(args.labels) as fh:
            user_labels = [ln.strip() for ln in fh if ln.strip()]

        if len(user_labels) == n:
            labels = user_labels
        else:
            print(
                f"WARNING: {len(user_labels)} labels for {n}×{n} matrix — "
                "using default labels instead."
            )

    # # Labels
    # labels = None
    # if args.labels:
    #     with open(args.labels) as fh:
    #         labels = [ln.strip() for ln in fh if ln.strip()]
    #     if len(labels) != n:
    #         print(f"  WARNING: {len(labels)} labels for {n}×{n} matrix — ignoring labels.")
    #         labels = None

    # Output base
    base = args.output or os.path.splitext(args.input)[0]
    fmt = args.format
    prefix = (args.title_prefix + " ") if args.title_prefix else ""

    # --- Covariance plot ---
    if args.mode in ("covariance", "both"):
        plot_matrix(
            mat=cov,
            title=f"{prefix}Covariance Matrix",
            cbar_label="Covariance",
            labels=labels,
            output_path=f"{base}_covariance.{fmt}",
        )

    # --- Correlation plot ---
    if args.mode in ("correlation", "both"):
        corr = covariance_to_correlation(cov)
        plot_matrix(
            mat=corr,
            title=f"{prefix}Correlation Matrix",
            cbar_label="Pearson correlation coefficient",
            labels=labels,
            output_path=f"{base}_correlation.{fmt}",
            vmin=-1.0,
            vmax=1.0,
        )

    print("\nDone.")


if __name__ == "__main__":
    main()