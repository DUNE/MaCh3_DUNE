#!/usr/bin/env python3
"""
plot_posterior_predictive_summary.py
-------------------------------------
Reads the *_summary.root file produced by MakePosteriorPredictiveSummary2D.C
and makes styled 2D plots, per detector/sample, of:

  1. Posterior mean prediction                        (events / bin)
  2. Relative posterior (systematic) uncertainty       sigma_post / mean
  3. Relative statistical uncertainty                  sigma_stat / mean  (~ 1/sqrt(N))
  4. Ratio of posterior to statistical uncertainty      sigma_post / sigma_stat
     (>1 => nuisance parameters dominate over counting stats in that bin)
  5. Posterior mean / Asimov                            (bias check, if an
     Asimov reference is present for that sample)

Usage:
    python plot_posterior_predictive_summary.py Posteriorpredictive_out_2D_summary.root [out.pdf]

Requirements:
    pip install uproot numpy matplotlib
"""

import sys
import pathlib
import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.backends.backend_pdf import PdfPages

DETECTORS = ["ND", "FD"]


# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

def safe_get(f, path):
    try:
        return f[path]
    except Exception:
        return None


def list_samples(summary_node):
    """Given a {det}/2D/summary directory, find distinct sample names by
    stripping the _2D_posteriorMean suffix off matching keys."""
    keys = [k.split(";")[0] for k in summary_node.keys()]
    suffix = "_2D_posteriorMean"
    return sorted(k[: -len(suffix)] for k in keys if k.endswith(suffix))


def plot_2d(ax, values, x_edges, y_edges, title, cmap="viridis", norm=None):
    mesh = ax.pcolormesh(x_edges, y_edges, values.T, cmap=cmap, norm=norm, shading="flat")
    ax.set_title(title, fontsize=11)
    ax.set_xlabel("True Neutrino Energy (GeV)", fontsize=9)
    ax.set_ylabel("Reco Neutrino Energy (GeV)", fontsize=9)
    ax.tick_params(labelsize=8)
    return mesh


def save_page(pdf, values, x_edges, y_edges, page_title, panel_title, cbar_label,
              cmap="viridis", norm=None):
    fig, ax = plt.subplots(figsize=(9, 6))
    fig.subplots_adjust(left=0.10, right=0.82, top=0.88, bottom=0.12)
    fig.suptitle(page_title, fontsize=12, fontweight="bold")

    mesh = plot_2d(ax, values, x_edges, y_edges, panel_title, cmap=cmap, norm=norm)

    cax = fig.add_axes([0.84, 0.12, 0.03, 0.76])
    cb = fig.colorbar(mesh, cax=cax)
    cb.set_label(cbar_label, fontsize=9)
    cb.ax.tick_params(labelsize=8)

    pdf.savefig(fig)
    plt.close(fig)


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():
    if len(sys.argv) < 2:
        print("Usage: python plot_posterior_predictive_summary.py summary.root [out.pdf]")
        sys.exit(1)

    in_file = sys.argv[1]
    out_pdf = sys.argv[2] if len(sys.argv) > 2 else "PosteriorPredictiveSummary.pdf"

    if not pathlib.Path(in_file).exists():
        print(f"[ERROR] File not found: {in_file}")
        sys.exit(1)

    f = uproot.open(in_file)
    pdf = PdfPages(out_pdf)
    n_plots = 0

    for det in DETECTORS:
        summary_node = safe_get(f, f"{det}/2D/summary")
        if summary_node is None:
            continue

        samples = list_samples(summary_node)
        if not samples:
            continue

        asimov_node = safe_get(f, f"Asimov/{det}/2D")
        asimov_keys = [k.split(";")[0] for k in asimov_node.keys()] if asimov_node is not None else []

        summary_keys = [k.split(";")[0] for k in summary_node.keys()]

        for sample in samples:
            mean_key = sample + "_2D_posteriorMean"
            post_key = sample + "_2D_posteriorErr"
            stat_key = sample + "_2D_statErr"

            if any(k not in summary_keys for k in (mean_key, post_key, stat_key)):
                print(f"[WARN] Incomplete set for {det}/{sample}, skipping")
                continue

            h_mean = summary_node[mean_key]
            h_post = summary_node[post_key]
            h_stat = summary_node[stat_key]

            mean_vals = h_mean.values()
            post_vals = h_post.values()
            stat_vals = h_stat.values()
            x_edges = h_mean.axes[0].edges()
            y_edges = h_mean.axes[1].edges()

            with np.errstate(divide="ignore", invalid="ignore"):
                rel_post = np.where(mean_vals > 0, post_vals / mean_vals, np.nan)
                rel_stat = np.where(mean_vals > 0, stat_vals / mean_vals, np.nan)
                ratio    = np.where(stat_vals > 0, post_vals / stat_vals, np.nan)

            page_title = f"{det} -- {sample}"

            # 1. posterior mean prediction (log scale -- spectra fall off steeply)
            positive_mean = mean_vals[mean_vals > 0]
            mean_norm = (mcolors.LogNorm(vmin=positive_mean.min(), vmax=mean_vals.max())
                         if positive_mean.size else None)
            save_page(pdf, mean_vals, x_edges, y_edges, page_title,
                      "Posterior Mean Prediction", "Events / bin",
                      cmap="viridis", norm=mean_norm)
            n_plots += 1

            # 2. relative posterior (systematic) uncertainty
            save_page(pdf, rel_post, x_edges, y_edges, page_title,
                      r"Relative Posterior (Systematic) Uncertainty:  $\sigma_{\mathrm{post}}/\mu$",
                      r"$\sigma_{\mathrm{post}} / \mu$", cmap="viridis")
            n_plots += 1

            # 3. relative statistical uncertainty (~ 1/sqrt(N))
            save_page(pdf, rel_stat, x_edges, y_edges, page_title,
                      r"Relative Statistical Uncertainty:  $\sigma_{\mathrm{stat}}/\mu \approx 1/\sqrt{N}$",
                      r"$\sigma_{\mathrm{stat}} / \mu$", cmap="viridis")
            n_plots += 1

            # 4. ratio: posterior / statistical uncertainty
            finite_ratio = ratio[np.isfinite(ratio)]
            vmax_ratio = max(2.0, np.percentile(finite_ratio, 95)) if finite_ratio.size else 2.0
            norm_ratio = mcolors.TwoSlopeNorm(vmin=0, vcenter=1, vmax=vmax_ratio)
            save_page(pdf, ratio, x_edges, y_edges, page_title,
                      r"$\sigma_{\mathrm{post}} / \sigma_{\mathrm{stat}}$   (> 1: systematics dominate)",
                      "Ratio", cmap="RdBu_r", norm=norm_ratio)
            n_plots += 1

            # 5. mean / Asimov, if a reference is present
            akey = sample + "_2D_Asimov"
            if akey in asimov_keys:
                asimov_vals = asimov_node[akey].values()
                with np.errstate(divide="ignore", invalid="ignore"):
                    ma = np.where(asimov_vals > 0, mean_vals / asimov_vals, np.nan)
                finite_ma = ma[np.isfinite(ma)]
                vmax_ma = max(1.5, np.percentile(finite_ma, 97)) if finite_ma.size else 1.5
                norm_ma = mcolors.TwoSlopeNorm(vmin=0, vcenter=1, vmax=vmax_ma)
                save_page(pdf, ma, x_edges, y_edges, page_title,
                          "Posterior Mean / Asimov", "Ratio",
                          cmap="RdBu_r", norm=norm_ma)
                n_plots += 1
            else:
                print(f"[WARN] No Asimov reference for {det}/{sample}, skipping Mean/Asimov plot")

    pdf.close()

    if n_plots == 0:
        print("[WARN] No plots were made -- check that the input file has the expected "
              "{det}/2D/summary layout (run MakePosteriorPredictiveSummary2D.C first).")
    else:
        print(f"[INFO] Saved {n_plots} plots to {out_pdf}")


if __name__ == "__main__":
    main()