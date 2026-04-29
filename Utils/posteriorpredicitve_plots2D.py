#!/usr/bin/env python3
"""
posteriorpredictive_plots.py
----------------------------------------
Plots Asimov and posterior-toy spectra from a ROOT output file.
Also plots 2D summary histograms: posterior uncertainty and statistical uncertainty.

Usage:
    python utils/posteriorpredictive_plots.py Posteriorpredictive_out.root [output.pdf]

Requirements:
    pip install uproot numpy matplotlib
"""

import sys
import re
import numpy as np
import uproot
import pathlib
import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.backends.backend_pdf import PdfPages


# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

def get_hist_values(obj):
    values = obj.values()
    axes = obj.axes

    if len(axes) == 1:
        return values, axes[0].edges()
    elif len(axes) == 2:
        return values, (axes[0].edges(), axes[1].edges())
    else:
        raise ValueError(f"Unsupported histogram dimension: {len(axes)}")


def safe_get_dir(f, path):
    try:
        return f[path]
    except Exception:
        return None


def extend_for_step(values, edges):
    return np.append(values, values[-1]), edges


def plot_2d_hist(ax, values, x_edges, y_edges, title, cbar_label,
                 cmap="viridis", norm=None,
                 xlabel="Reco Energy (GeV)",
                 ylabel="Lepton Energy (GeV)"):

    mesh = ax.pcolormesh(
        x_edges, y_edges, values.T,
        cmap=cmap, norm=norm, shading="flat"
    )
    ax.set_title(title, fontsize=11)
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.tick_params(labelsize=8)
    return mesh


def save_single_2d(pdf_pages, values, x_edges, y_edges,
                   page_title, plot_title, cbar_label,
                   cmap="viridis", norm=None):

    fig, ax = plt.subplots(figsize=(9, 6))
    fig.subplots_adjust(left=0.10, right=0.82, top=0.88, bottom=0.12)
    fig.suptitle(page_title, fontsize=12, fontweight="bold")

    mesh = plot_2d_hist(ax, values, x_edges, y_edges,
                        plot_title, cbar_label,
                        cmap=cmap, norm=norm)

    cax = fig.add_axes([0.84, 0.12, 0.03, 0.76])
    cb = fig.colorbar(mesh, cax=cax)
    cb.set_label(cbar_label, fontsize=9)
    cb.ax.tick_params(labelsize=8)

    pdf_pages.savefig(fig)
    plt.close(fig)


# ------------------------------------------------------------
# 2D summaries
# ------------------------------------------------------------

def plot_2d_summaries(f, detectors, pdf_pages):

    print("\n[INFO] Generating 2D summary plots...")

    for det in detectors:
        summary_path = f"{det}/2D/summary"
        summary_node = safe_get_dir(f, summary_path)
        if summary_node is None:
            print(f"[WARN] Missing {summary_path}")
            continue

        keys = [k.split(";")[0] for k in summary_node.keys()]

        suffix = re.compile(r"_2D_(posteriorMean|posteriorErr|statErr)$")
        samples = set()

        for k in keys:
            m = suffix.search(k)
            if m:
                samples.add(k[:m.start()])

        for sample in sorted(samples):

            mean_key = sample + "_2D_posteriorMean"
            post_key = sample + "_2D_posteriorErr"
            stat_key = sample + "_2D_statErr"

            if any(k not in keys for k in [mean_key, post_key, stat_key]):
                continue

            h_mean = summary_node[mean_key]
            h_post = summary_node[post_key]
            h_stat = summary_node[stat_key]

            mean_vals = h_mean.values()
            post_vals = h_post.values()
            stat_vals = h_stat.values()

            x_edges = h_mean.axes[0].edges()
            y_edges = h_mean.axes[1].edges()

            asimov_vals = None
            asimov_node = safe_get_dir(f, f"Asimov/{det}/2D")

            if asimov_node is not None:
                akey = f"{sample}_2D_Asimov"
                akeys = [k.split(";")[0] for k in asimov_node.keys()]
                if akey in akeys:
                    asimov_vals = asimov_node[akey].values()

            with np.errstate(divide="ignore", invalid="ignore"):
                rel_post = np.where(mean_vals > 0, post_vals / mean_vals, np.nan)
                rel_stat = np.where(mean_vals > 0, stat_vals / mean_vals, np.nan)
                ratio = np.where(stat_vals > 0, post_vals / stat_vals, np.nan)

            vmax = max(2.0, np.nanpercentile(ratio, 95))
            norm_ratio = mcolors.TwoSlopeNorm(vmin=0, vcenter=1, vmax=vmax)

            save_single_2d(pdf_pages, mean_vals, x_edges, y_edges,
                           f"{sample} — Posterior Mean",
                           "Posterior Mean",
                           "Events")

            save_single_2d(pdf_pages, rel_post, x_edges, y_edges,
                           f"{sample} — Posterior Uncertainty",
                           "σ_post / N",
                           "Rel unc")

            save_single_2d(pdf_pages, rel_stat, x_edges, y_edges,
                           f"{sample} — Statistical Uncertainty",
                           "1/sqrt(N)",
                           "Rel stat")

            save_single_2d(pdf_pages, ratio, x_edges, y_edges,
                           f"{sample} — σ_post / σ_stat",
                           "Posterior / Stat",
                           "Ratio",
                           cmap="RdBu_r",
                           norm=norm_ratio)

            if asimov_vals is not None:
                with np.errstate(divide="ignore", invalid="ignore"):
                    ma = np.where(asimov_vals > 0,
                                  mean_vals / asimov_vals,
                                  np.nan)

                vmax_ma = max(1.5, np.nanpercentile(ma, 97))
                norm_ma = mcolors.TwoSlopeNorm(vmin=0, vcenter=1, vmax=vmax_ma)

                save_single_2d(pdf_pages, ma, x_edges, y_edges,
                               f"{sample} — Mean / Asimov",
                               "Mean / Asimov",
                               "Ratio",
                               cmap="RdBu_r",
                               norm=norm_ma)


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():

    if len(sys.argv) < 2:
        print("Usage: python posteriorpredictive_plots.py input.root [out.pdf]")
        sys.exit(1)

    root_file = sys.argv[1]
    output_pdf = sys.argv[2] if len(sys.argv) > 2 else "SpectraPlots.pdf"

    if not pathlib.Path(root_file).exists():
        print("[ERROR] File not found")
        sys.exit(1)

    f = uproot.open(root_file)

    detectors = [d for d in ["ND", "FD", "Other"] if safe_get_dir(f, d)]

    variables = ["RecoNeutrinoEnergy", "TrueNeutrinoEnergy"]

    pdf = PdfPages(output_pdf)

    # --------------------------------------------------------
    # 1D plots
    # --------------------------------------------------------
    for det in detectors:
        for var in variables:

            asimov_node = safe_get_dir(f, "Asimov")
            if asimov_node is None:
                continue

            asimov_keys = [k for k in asimov_node.keys()
                           if det in k and var in k and "Asimov" in k]

            if not asimov_keys:
                continue

            asimov_vals, edges = get_hist_values(asimov_node[asimov_keys[0]])

            # force 1D projection if histogram is 2D
            if isinstance(edges, tuple):
                # collapse to 1D before plotting
                asimov_vals = np.sum(asimov_vals, axis=1)
                edges = edges[0]

            if isinstance(edges, tuple):
                x_edges, y_edges = edges
                xw = np.diff(x_edges)
                yw = np.diff(y_edges)
                bin_area = xw[:, None] * yw[None, :]
                asimov_vals = asimov_vals / bin_area
            else:
                bw = np.diff(edges)
                asimov_vals = asimov_vals / bw

            toy_node = safe_get_dir(f, f"{det}/{var}")
            if toy_node is None:
                continue

            toy_keys = [k for k in toy_node.keys()
                        if "posterior_toy_" in k]

            toys = []

            for k in toy_keys:
                h = toy_node[k]
                v = h.values()
                bw = np.diff(h.axes[0].edges())
                if len(v) == len(asimov_vals):
                    toys.append(v / bw)

            if not toys:
                continue

            toys = np.vstack(toys)
            mean = np.mean(toys, axis=0)
            std = np.std(toys, axis=0)

            edges_plot = edges if not isinstance(edges, tuple) else edges[0]

            mean_p, _ = extend_for_step(mean, edges_plot)
            up, _ = extend_for_step(mean + std, edges_plot)
            lo, _ = extend_for_step(mean - std, edges_plot)

            asimov_p, _ = extend_for_step(asimov_vals, edges_plot)

            fig, ax = plt.subplots()
            ax.step(edges_plot, mean_p, where="post", label="Posterior")
            ax.step(edges_plot, asimov_p, where="post", label="Asimov")
            ax.fill_between(edges_plot, lo, up, step="post", alpha=0.3)
            ax.legend()

            pdf.savefig(fig)
            plt.close(fig)

    # --------------------------------------------------------
    # 2D plots
    # --------------------------------------------------------
    plot_2d_summaries(f, detectors, pdf)

    pdf.close()
    print(f"[INFO] Saved {output_pdf}")


if __name__ == "__main__":
    main()