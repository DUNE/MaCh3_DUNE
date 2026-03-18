#!/usr/bin/env python3
"""
posteriorpredictive_comparison.py
----------------------------------------
Compare posterior predictive spectra from two ROOT output files.

Usage:
    python posteriorpredictive_comparison.py on+off.root on.root [output.pdf]

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
    """Return (values, edges) from an uproot TH1 object."""
    return obj.values(), obj.axis().edges()

def safe_get_dir(f, path):
    """Safely get a directory from uproot file."""
    try:
        return f[path]
    except Exception:
        return None

def extend_for_step(values, edges):
    """
    Extend a values array so that step(..., where='post') drawn against
    the full edges array will render every bin including the last one.

    Parameters
    ----------
    values : 1D array of length N
    edges  : 1D array of length N+1

    Returns
    -------
    values_ext : 1D array of length N+1 (last value repeated)
    edges      : unchanged, returned for convenience
    """
    return np.append(values, values[-1]), edges

def load_toys(f, det, var):
    """Return mean, std, and edges of posterior toys for a given detector and variable."""
    toy_dir = f"{det}/{var}"
    toy_node = safe_get_dir(f, toy_dir)
    if not toy_node:
        return None, None, None

    all_keys = list(toy_node.keys())
    toy_keys = [k for k in all_keys if f"{var}_posterior_toy_" in k and k.startswith(f"{det}_")]
    if not toy_keys:
        return None, None, None

    toy_keys = sorted(toy_keys,
                      key=lambda x: int(re.search(r"(\d+)$", x).group(1))
                                    if re.search(r"(\d+)$", x) else 0)

    toys  = []
    edges = None
    for key in toy_keys:
        try:
            h = toy_node[key]
            vals, e = get_hist_values(h)
            bin_widths = np.diff(e)
            if edges is None:
                edges = e
            if len(vals) == len(bin_widths):
                toys.append(vals / bin_widths)
        except Exception:
            continue

    if not toys:
        return None, None, None

    toys = np.vstack(toys)
    return np.mean(toys, axis=0), np.std(toys, axis=0), edges


def load_2d_summary(f, det):
    """
    Return a dict of 2D summary arrays for each sample found under
    <det>/2D/summary.  Keys are sample names; values are dicts with
    keys 'mean', 'post_err', 'stat_err', 'x_edges', 'y_edges'.
    """
    summary_path = f"{det}/2D/summary"
    summary_node = safe_get_dir(f, summary_path)
    if summary_node is None:
        return {}

    all_keys_clean = [k.split(";")[0] for k in summary_node.keys()]
    suffix_pattern = re.compile(r"_2D_(posteriorMean|posteriorErr|statErr)$")

    sample_names = set()
    for k in all_keys_clean:
        m = suffix_pattern.search(k)
        if m:
            sample_names.add(k[: m.start()])

    results = {}
    for sample in sorted(sample_names):
        mean_key = sample + "_2D_posteriorMean"
        post_key = sample + "_2D_posteriorErr"
        stat_key = sample + "_2D_statErr"

        if any(k not in all_keys_clean for k in (mean_key, post_key, stat_key)):
            continue

        try:
            h_mean = summary_node[mean_key]
            h_post = summary_node[post_key]
            h_stat = summary_node[stat_key]
        except Exception as e:
            print(f"[WARN] Could not read 2D summaries for {sample}: {e}")
            continue

        results[sample] = {
            "mean":     h_mean.values(),
            "post_err": h_post.values(),
            "stat_err": h_stat.values(),
            "x_edges":  h_mean.axes[0].edges(),
            "y_edges":  h_mean.axes[1].edges(),
        }

    return results


def save_single_2d(pdf_pages, values, x_edges, y_edges,
                   page_title, plot_title, cbar_label,
                   cmap="viridis", norm=None,
                   xlabel="Reconstructed Neutrino Energy (GeV)",
                   ylabel="Lepton Energy (GeV)"):
    """Save one 2D histogram as a single full-page figure."""
    fig, ax = plt.subplots(figsize=(9, 6))
    fig.subplots_adjust(left=0.10, right=0.82, top=0.88, bottom=0.12)
    fig.suptitle(page_title, fontsize=12, fontweight="bold")

    mesh = ax.pcolormesh(x_edges, y_edges, values.T,
                         cmap=cmap, norm=norm, shading="flat")
    ax.set_title(plot_title, fontsize=11)
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.tick_params(labelsize=8)

    cax = fig.add_axes([0.84, 0.12, 0.03, 0.76])
    cb  = fig.colorbar(mesh, cax=cax)
    cb.set_label(cbar_label, fontsize=9)
    cb.ax.tick_params(labelsize=8)

    pdf_pages.savefig(fig)
    plt.close(fig)


# ------------------------------------------------------------
# 2D comparison plots
# ------------------------------------------------------------

def plot_2d_comparisons(fA, fB, detectors, pdf_pages,
                        label_a="On+Off-axis", label_b="On-axis"):
    """
    For each detector and each sample found in both files, produce:
        Page 1 : Posterior Mean  — file A
        Page 2 : Posterior Mean  — file B
        Page 3 : Mean A / Mean B  (ratio, diverging around 1)
        Page 4 : Posterior σ     — file A
        Page 5 : Posterior σ     — file B
        Page 6 : Posterior σ A / Posterior σ B
        Page 7 : Relative posterior uncertainty (σ/mean) A
        Page 8 : Relative posterior uncertainty (σ/mean) B
        Page 9 : (σ/mean) A / (σ/mean) B
    """
    print("\n[INFO] Generating 2D comparison plots...")

    for det in detectors:
        summaries_a = load_2d_summary(fA, det)
        summaries_b = load_2d_summary(fB, det)

        common_samples = sorted(set(summaries_a) & set(summaries_b))
        if not common_samples:
            print(f"[WARN] No common 2D summary samples for detector {det}")
            continue

        for sample in common_samples:
            sa = summaries_a[sample]
            sb = summaries_b[sample]

            x_edges = sa["x_edges"]
            y_edges = sa["y_edges"]

            mean_a     = sa["mean"]
            post_err_a = sa["post_err"]
            mean_b     = sb["mean"]
            post_err_b = sb["post_err"]

            with np.errstate(invalid="ignore", divide="ignore"):
                rel_a = np.where(mean_a > 0, post_err_a / mean_a, np.nan)
                rel_b = np.where(mean_b > 0, post_err_b / mean_b, np.nan)

            # Ratios — default to NaN where denominator is zero
            with np.errstate(invalid="ignore", divide="ignore"):
                ratio_mean     = np.where(mean_b     > 0, mean_a     / mean_b,     np.nan)
                ratio_post_err = np.where(post_err_b > 0, post_err_a / post_err_b, np.nan)
                ratio_rel      = np.where(rel_b      > 0, rel_a      / rel_b,      np.nan)

            def diverging_norm(arr, vcenter=1.0, percentile=97):
                finite = arr[np.isfinite(arr)]
                if finite.size == 0:
                    return mcolors.TwoSlopeNorm(vmin=0.5, vcenter=vcenter, vmax=1.5)
                vmin = max(0.0, float(np.nanpercentile(arr, 100 - percentile)))
                vmax = float(np.nanpercentile(arr, percentile))
                vmax = max(vmax, vcenter + 1e-6)
                vmin = min(vmin, vcenter - 1e-6)
                return mcolors.TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)

            # ---- Posterior Mean ----
            vmax_mean = max(float(np.nanmax(mean_a)), float(np.nanmax(mean_b)), 1e-6)
            norm_mean = mcolors.Normalize(vmin=0, vmax=vmax_mean)

            save_single_2d(pdf_pages, mean_a, x_edges, y_edges,
                           page_title=f"{sample} — Posterior Mean ({label_a})",
                           plot_title=f"Posterior Mean — {label_a}",
                           cbar_label="Events / bin", cmap="Blues", norm=norm_mean)

            save_single_2d(pdf_pages, mean_b, x_edges, y_edges,
                           page_title=f"{sample} — Posterior Mean ({label_b})",
                           plot_title=f"Posterior Mean — {label_b}",
                           cbar_label="Events / bin", cmap="Blues", norm=norm_mean)

            save_single_2d(pdf_pages, ratio_mean, x_edges, y_edges,
                           page_title=f"{sample} — Posterior Mean Ratio",
                           plot_title=f"Posterior Mean: {label_a} / {label_b}  (1 = identical)",
                           cbar_label=f"{label_a} / {label_b}",
                           cmap="RdBu_r", norm=diverging_norm(ratio_mean))

            # ---- Posterior σ ----
            vmax_err = max(float(np.nanmax(post_err_a)), float(np.nanmax(post_err_b)), 1e-6)
            norm_err = mcolors.Normalize(vmin=0, vmax=vmax_err)

            save_single_2d(pdf_pages, post_err_a, x_edges, y_edges,
                           page_title=f"{sample} — Posterior σ ({label_a})",
                           plot_title=f"Posterior σ — {label_a}",
                           cbar_label="σ (events / bin)", cmap="Reds", norm=norm_err)

            save_single_2d(pdf_pages, post_err_b, x_edges, y_edges,
                           page_title=f"{sample} — Posterior σ ({label_b})",
                           plot_title=f"Posterior σ — {label_b}",
                           cbar_label="σ (events / bin)", cmap="Reds", norm=norm_err)

            save_single_2d(pdf_pages, ratio_post_err, x_edges, y_edges,
                           page_title=f"{sample} — Posterior σ Ratio",
                           plot_title=f"Posterior σ: {label_a} / {label_b}  (1 = identical)",
                           cbar_label=f"{label_a} / {label_b}",
                           cmap="RdBu_r", norm=diverging_norm(ratio_post_err))

            # ---- Relative uncertainty σ/mean ----
            vmax_rel = max(float(np.nanmax(rel_a[np.isfinite(rel_a)])) if np.any(np.isfinite(rel_a)) else 1.0,
                           float(np.nanmax(rel_b[np.isfinite(rel_b)])) if np.any(np.isfinite(rel_b)) else 1.0,
                           1e-6)
            norm_rel = mcolors.Normalize(vmin=0, vmax=vmax_rel)

            save_single_2d(pdf_pages, rel_a, x_edges, y_edges,
                           page_title=f"{sample} — Relative Posterior Uncertainty ({label_a})",
                           plot_title=f"σ / Mean — {label_a}",
                           cbar_label="σ / Mean", cmap="Purples", norm=norm_rel)

            save_single_2d(pdf_pages, rel_b, x_edges, y_edges,
                           page_title=f"{sample} — Relative Posterior Uncertainty ({label_b})",
                           plot_title=f"σ / Mean — {label_b}",
                           cbar_label="σ / Mean", cmap="Purples", norm=norm_rel)

            save_single_2d(pdf_pages, ratio_rel, x_edges, y_edges,
                           page_title=f"{sample} — Relative Uncertainty Ratio",
                           plot_title=f"(σ/Mean): {label_a} / {label_b}  (1 = identical)",
                           cbar_label=f"{label_a} / {label_b}",
                           cmap="RdBu_r", norm=diverging_norm(ratio_rel))

            print(f"[INFO] Added 9 2D comparison pages for {sample}")

    print("[INFO] 2D comparison plots done.")


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():
    if len(sys.argv) < 3:
        print("Usage: python posteriorpredictive_comparison.py on+off.root on.root [output.pdf]")
        sys.exit(1)

    fileA_path = sys.argv[1]
    fileB_path = sys.argv[2]
    output_pdf = sys.argv[3] if len(sys.argv) > 3 else "PosteriorComparison.pdf"

    for path in [fileA_path, fileB_path]:
        if not pathlib.Path(path).exists():
            print(f"[ERROR] File not found: {path}")
            sys.exit(1)

    try:
        fA = uproot.open(fileA_path)
        print(f"[INFO] Opened {fileA_path}")
    except Exception as e:
        print(f"[ERROR] Failed to open {fileA_path}: {e}")
        sys.exit(1)

    try:
        fB = uproot.open(fileB_path)
        print(f"[INFO] Opened {fileB_path}")
    except Exception as e:
        print(f"[ERROR] Failed to open {fileB_path}: {e}")
        sys.exit(1)

    label_a = "On+Off Axis"
    label_b = "On Axis"

    detectors = ["ND", "FD", "Other"]
    variables  = ["RecoNeutrinoEnergy", "TrueNeutrinoEnergy"]

    pdf = PdfPages(output_pdf)
    plots_created = 0

    # ----------------------------------------------------------
    # 1D spectra comparison
    # ----------------------------------------------------------
    for det in detectors:
        for var in variables:
            mean_a, std_a, edges_a = load_toys(fA, det, var)
            mean_b, std_b, edges_b = load_toys(fB, det, var)

            if mean_a is None or mean_b is None:
                print(f"[WARN] Skipping {det}/{var} — missing data in one or both files")
                continue

            if not np.allclose(edges_a, edges_b):
                print(f"[WARN] Bin edges differ for {det}/{var}, using file A binning")
            edges = edges_a

            # Extend all arrays so the final bin is rendered
            mean_a_plot,       _ = extend_for_step(mean_a,            edges)
            mean_b_plot,       _ = extend_for_step(mean_b,            edges)
            upper_a_plot,      _ = extend_for_step(mean_a + std_a,    edges)
            lower_a_plot,      _ = extend_for_step(mean_a - std_a,    edges)
            upper_b_plot,      _ = extend_for_step(mean_b + std_b,    edges)
            lower_b_plot,      _ = extend_for_step(mean_b - std_b,    edges)

            # ---- Event rate plot ----
            fig, (ax_top, ax_ratio) = plt.subplots(
                2, 1, figsize=(7, 5),
                gridspec_kw={"height_ratios": [3, 1]},
                sharex=True
            )
            ax_top.set_title(f"{det} — {var}", fontsize=12)

            ax_top.fill_between(edges, lower_a_plot, upper_a_plot,
                                step="post", alpha=0.3, color="red",
                                label=f"{label_a} ±1σ")
            ax_top.step(edges, mean_a_plot, where="post",
                        color="red", lw=1.5, label=f"{label_a} mean")

            ax_top.fill_between(edges, lower_b_plot, upper_b_plot,
                                step="post", alpha=0.3, color="blue",
                                label=f"{label_b} ±1σ")
            ax_top.step(edges, mean_b_plot, where="post",
                        color="blue", lw=1.5, label=f"{label_b} mean")

            ax_top.set_ylabel("Events / GeV", fontsize=11)
            ax_top.legend(frameon=False, loc="best")
            ax_top.grid(alpha=0.3)
            ax_top.set_xlim(edges[0], edges[-1])

            # Ratio: A / B
            ratio = np.divide(mean_a, mean_b,
                              out=np.ones_like(mean_a), where=mean_b != 0)
            ratio_plot, _ = extend_for_step(ratio, edges)

            ax_ratio.step(edges, ratio_plot, where="post", color="black", lw=1.2)
            ax_ratio.axhline(1.0, color="gray", lw=0.8, ls="--")
            ax_ratio.set_ylabel(f"{label_a} / {label_b}", fontsize=10)
            ax_ratio.set_xlabel(f"{var} (GeV)", fontsize=11)
            ax_ratio.grid(alpha=0.3)
            ax_ratio.set_xlim(edges[0], edges[-1])

            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
            plots_created += 1

            # ---- Relative uncertainty plot ----
            rel_a = np.divide(std_a, mean_a, out=np.zeros_like(std_a), where=mean_a != 0)
            rel_b = np.divide(std_b, mean_b, out=np.zeros_like(std_b), where=mean_b != 0)

            rel_a_plot,  _ = extend_for_step(rel_a, edges)
            rel_b_plot,  _ = extend_for_step(rel_b, edges)

            fig, (ax_top, ax_ratio) = plt.subplots(
                2, 1, figsize=(7, 5),
                gridspec_kw={"height_ratios": [3, 1]},
                sharex=True
            )
            ax_top.set_title(f"{det} — {var} Relative Uncertainty", fontsize=12)

            ax_top.step(edges, rel_a_plot, where="post",
                        color="red",  lw=1.5, label=label_a)
            ax_top.step(edges, rel_b_plot, where="post",
                        color="blue", lw=1.5, label=label_b)

            ax_top.set_ylabel("σ / Mean", fontsize=11)
            ax_top.legend(frameon=False, loc="best")
            ax_top.grid(alpha=0.3)
            ax_top.set_xlim(edges[0], edges[-1])

            rel_ratio = np.divide(rel_a, rel_b,
                                  out=np.ones_like(rel_a), where=rel_b != 0)
            rel_ratio_plot, _ = extend_for_step(rel_ratio, edges)

            ax_ratio.step(edges, rel_ratio_plot, where="post", color="black", lw=1.2)
            ax_ratio.axhline(1.0, color="gray", lw=0.8, ls="--")
            ax_ratio.set_ylabel(f"{label_a} / {label_b}", fontsize=10)
            ax_ratio.set_xlabel(f"{var} (GeV)", fontsize=11)
            ax_ratio.grid(alpha=0.3)
            ax_ratio.set_xlim(edges[0], edges[-1])

            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
            plots_created += 1

            print(f"[INFO] 1D plots added for {det}/{var}")

    # ----------------------------------------------------------
    # 2D comparison plots
    # ----------------------------------------------------------
    plot_2d_comparisons(fA, fB, detectors, pdf, label_a=label_a, label_b=label_b)

    pdf.close()

    if plots_created > 0:
        print(f"[INFO] Saved plots to {output_pdf}")
    else:
        print("[WARN] No plots created — check that your ROOT files contain the expected data")


# ------------------------------------------------------------
# Entry point
# ------------------------------------------------------------
if __name__ == "__main__":
    main()