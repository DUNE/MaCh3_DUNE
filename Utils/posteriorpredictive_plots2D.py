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
    """Return (values, edges) from an uproot TH1 object."""
    return obj.values(), obj.axis().edges()


def safe_get_dir(f, path):
    """Safely get a directory from an uproot file; returns None on failure."""
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
    values : 1D array of length N  (one entry per bin)
    edges  : 1D array of length N+1

    Returns
    -------
    values_ext : 1D array of length N+1  (last value repeated)
    edges      : unchanged, returned for convenience
    """
    return np.append(values, values[-1]), edges


def plot_2d_hist(ax, values, x_edges, y_edges, title, cbar_label,
                 cmap="viridis", norm=None,
                 xlabel="Reconstructed Neutrino Energy (GeV)",
                 ylabel="Lepton Energy (GeV)"):
    """
    Draw a 2D histogram on ax using pcolormesh.
    values : 2D numpy array shape (nX, nY)
    Returns the mappable for colorbar attachment.
    """
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
    """Save one 2D histogram as a single full-page figure."""
    fig, ax = plt.subplots(figsize=(9, 6))
    fig.subplots_adjust(left=0.10, right=0.82, top=0.88, bottom=0.12)
    fig.suptitle(page_title, fontsize=12, fontweight="bold")
    mesh = plot_2d_hist(ax, values, x_edges, y_edges, plot_title, cbar_label,
                        cmap=cmap, norm=norm)
    cax = fig.add_axes([0.84, 0.12, 0.03, 0.76])
    cb = fig.colorbar(mesh, cax=cax)
    cb.set_label(cbar_label, fontsize=9)
    cb.ax.tick_params(labelsize=8)
    pdf_pages.savefig(fig)
    plt.close(fig)


# ------------------------------------------------------------
# 2D summary plots — one plot per page
# ------------------------------------------------------------

def plot_2d_summaries(f, detectors, pdf_pages):
    """
    For each detector read the three summary 2D histograms written by the C++ code:
        <det>/2D/summary/<det>_<title>_2D_posteriorMean
        <det>/2D/summary/<det>_<title>_2D_posteriorErr
        <det>/2D/summary/<det>_<title>_2D_statErr

    Produces five separate pages per sample:
        Page 1: Posterior Mean
        Page 2: Relative posterior uncertainty sigma / N
        Page 3: Relative statistical uncertainty 1/sqrt(N)
        Page 4: Posterior sigma / Statistical sigma  (systematic inflation)
        Page 5: Posterior Mean / Asimov  (bin-by-bin ratio to nominal)
    """
    print("\n[INFO] Generating 2D summary plots...")

    for det in detectors:
        summary_path = f"{det}/2D/summary"
        summary_node = safe_get_dir(f, summary_path)
        if summary_node is None:
            print(f"[WARN] No 2D summary directory found at '{summary_path}' — skipping")
            continue

        all_keys = list(summary_node.keys())
        all_keys_clean = [k.split(";")[0] for k in all_keys]
        print(f"[DEBUG] Keys in {summary_path}: {all_keys_clean}")

        suffix_pattern = re.compile(r"_2D_(posteriorMean|posteriorErr|statErr)$")
        sample_names = set()
        for k in all_keys_clean:
            m = suffix_pattern.search(k)
            if m:
                sample_names.add(k[: m.start()])

        if not sample_names:
            print(f"[WARN] No recognised summary histograms in {summary_path}")
            continue

        for sample in sorted(sample_names):
            mean_key = sample + "_2D_posteriorMean"
            post_key = sample + "_2D_posteriorErr"
            stat_key = sample + "_2D_statErr"

            missing = [k for k in (mean_key, post_key, stat_key)
                       if k not in all_keys_clean]
            if missing:
                print(f"[WARN] Missing histograms for sample '{sample}': {missing}")
                continue

            try:
                h_mean = summary_node[mean_key]
                h_post = summary_node[post_key]
                h_stat = summary_node[stat_key]
            except Exception as e:
                print(f"[WARN] Could not read histograms for {sample}: {e}")
                continue

            mean_vals = h_mean.values()
            post_vals = h_post.values()
            stat_vals = h_stat.values()
            x_edges   = h_mean.axes[0].edges()
            y_edges   = h_mean.axes[1].edges()

            # Load matching Asimov 2D histogram
            asimov_2d_node = safe_get_dir(f, f"Asimov/{det}/2D")
            asimov_vals_2d = None
            if asimov_2d_node is not None:
                asimov_2d_key = f"{sample}_2D_Asimov"
                asimov_keys_clean = [k.split(";")[0] for k in asimov_2d_node.keys()]
                if asimov_2d_key in asimov_keys_clean:
                    try:
                        asimov_vals_2d = asimov_2d_node[asimov_2d_key].values()
                    except Exception as e:
                        print(f"[WARN] Could not load Asimov 2D for {sample}: {e}")
                else:
                    print(f"[WARN] Asimov 2D key '{asimov_2d_key}' not found — "
                          f"available: {asimov_keys_clean}")
            else:
                print(f"[WARN] Asimov/{det}/2D directory not found")

            print(f"[DEBUG] {sample}: shape={mean_vals.shape}, "
                  f"posteriorErr max={np.nanmax(post_vals):.2f}, "
                  f"statErr max={np.nanmax(stat_vals):.2f}")

            with np.errstate(invalid="ignore", divide="ignore"):
                rel_post = np.where(mean_vals > 0, post_vals / mean_vals, np.nan)
                rel_stat = np.where(mean_vals > 0, stat_vals / mean_vals, np.nan)

            with np.errstate(invalid="ignore", divide="ignore"):
                ratio_vals = np.where(stat_vals > 0, post_vals / stat_vals, np.nan)

            ratio_finite = ratio_vals[np.isfinite(ratio_vals)]
            vmax = max(2.0, float(np.nanpercentile(ratio_vals, 95))) \
                   if ratio_finite.size > 0 else 2.0
            norm_ratio = mcolors.TwoSlopeNorm(vmin=0.0, vcenter=1.0, vmax=vmax)

            save_single_2d(
                pdf_pages, mean_vals, x_edges, y_edges,
                page_title=f"{sample} — Posterior Mean",
                plot_title="Posterior Mean",
                cbar_label="Events / bin",
                cmap="Blues"
            )

            save_single_2d(
                pdf_pages, rel_post, x_edges, y_edges,
                page_title=f"{sample} — Relative Posterior Uncertainty",
                plot_title="Posterior σ / N  (relative uncertainty)",
                cbar_label="σₚₒₛₜ / N",
                cmap="Reds"
            )

            save_single_2d(
                pdf_pages, rel_stat, x_edges, y_edges,
                page_title=f"{sample} — Relative Statistical Uncertainty",
                plot_title="Statistical σ / N = 1/√N  (relative stat. uncertainty)",
                cbar_label="σₛₜₐₜ / N",
                cmap="Greens"
            )

            save_single_2d(
                pdf_pages, ratio_vals, x_edges, y_edges,
                page_title=f"{sample} — Posterior σ / Statistical σ",
                plot_title="Posterior σ / Statistical σ   (>1 ⟹ syst. dominated)",
                cbar_label="σₚₒₛₜ / σₛₜₐₜ",
                cmap="RdBu_r",
                norm=norm_ratio
            )

            if asimov_vals_2d is not None:
                with np.errstate(invalid="ignore", divide="ignore"):
                    mean_over_asimov = np.where(
                        asimov_vals_2d > 0,
                        mean_vals / asimov_vals_2d,
                        np.nan
                    )
                ratio_finite_ma = mean_over_asimov[np.isfinite(mean_over_asimov)]
                vmax_ma = max(1.5, float(np.nanpercentile(mean_over_asimov, 97))) \
                          if ratio_finite_ma.size > 0 else 1.5
                norm_ma = mcolors.TwoSlopeNorm(vmin=0.0, vcenter=1.0, vmax=vmax_ma)
                save_single_2d(
                    pdf_pages, mean_over_asimov, x_edges, y_edges,
                    page_title=f"{sample} — Posterior Mean / Asimov",
                    plot_title="Posterior Mean / Asimov  (1 = nominal)",
                    cbar_label="Posterior Mean / Asimov",
                    cmap="RdBu_r",
                    norm=norm_ma
                )
            else:
                print(f"[WARN] Skipping Posterior Mean / Asimov page for {sample} "
                      "(Asimov 2D not found)")

            print(f"[INFO] Added pages for {sample}")

    print("[INFO] 2D summary plots done.")


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():
    if len(sys.argv) < 2:
        print("Usage: python utils/posteriorpredictive_plots.py "
              "Posteriorpredictive_out.root [output.pdf]")
        sys.exit(1)

    root_file  = sys.argv[1]
    output_pdf = sys.argv[2] if len(sys.argv) > 2 else "SpectraPlots.pdf"

    if not pathlib.Path(root_file).exists():
        print(f"[ERROR] Input file not found: {root_file}")
        sys.exit(1)

    print(f"[INFO] Opening ROOT file: {root_file}")
    f = uproot.open(root_file)

    # Detect which detectors are actually present
    detectors = []
    for candidate in ["ND", "FD", "Other"]:
        if safe_get_dir(f, candidate) is not None:
            detectors.append(candidate)
    print(f"[INFO] Found detectors: {detectors}")

    variables = ["RecoNeutrinoEnergy", "TrueNeutrinoEnergy"]

    pdf_pages = PdfPages(output_pdf)

    # ----------------------------------------------------------
    # First pass: determine consistent y-axis limits per detector
    # ----------------------------------------------------------
    y_limits     = {}
    ratio_limits = {}

    for det in detectors:
        max_y     = 0.0
        max_ratio = 0.0

        for var in variables:
            asimov_node = safe_get_dir(f, "Asimov")
            if asimov_node is None:
                continue

            asimov_keys = [k for k in asimov_node.keys()
                           if f"{det}_" in k and var in k and "_Asimov" in k]
            if not asimov_keys:
                continue

            asimov_vals, edges = get_hist_values(asimov_node[asimov_keys[0]])
            bin_widths  = np.diff(edges)
            asimov_vals = asimov_vals / bin_widths

            toy_node = safe_get_dir(f, f"{det}/{var}")
            if toy_node is None:
                continue

            toy_keys = [k for k in toy_node.keys()
                        if k.startswith(f"{det}_") and f"{var}_posterior_toy_" in k]
            if not toy_keys:
                continue

            toys = []
            for key in toy_keys:
                try:
                    h  = toy_node[key]
                    bw = np.diff(h.axis().edges())
                    v  = h.values()
                    if len(v) == len(asimov_vals):
                        toys.append(v / bw)
                except Exception:
                    pass

            if toys:
                toys      = np.vstack(toys)
                mean_toys = np.mean(toys, axis=0)
                std_toys  = np.std(toys,  axis=0)
                max_y     = max(max_y,
                                np.max(asimov_vals),
                                np.max(mean_toys + std_toys))
                rm = np.divide(std_toys, mean_toys,
                               out=np.zeros_like(std_toys), where=mean_toys != 0)
                ra = np.divide(std_toys, asimov_vals,
                               out=np.zeros_like(std_toys), where=asimov_vals != 0)
                max_ratio = max(max_ratio, np.max(rm), np.max(ra))

        y_limits[det]     = max_y     * 1.1
        ratio_limits[det] = max_ratio * 1.4

    # ----------------------------------------------------------
    # Second pass: 1D spectra plots
    # ----------------------------------------------------------
    for det in detectors:
        for var in variables:
            asimov_node = safe_get_dir(f, "Asimov")
            if asimov_node is None:
                print("[WARN] Missing Asimov directory")
                continue

            asimov_keys = [k for k in asimov_node.keys()
                           if f"{det}_" in k and var in k and "_Asimov" in k]
            if not asimov_keys:
                print(f"[WARN] No Asimov histograms found for {det}/{var}")
                print(f"[DEBUG] Available keys: {list(asimov_node.keys())[:5]}...")
                continue

            print(f"[DEBUG] Using Asimov histogram: {asimov_keys[0]}")
            asimov_vals, edges = get_hist_values(asimov_node[asimov_keys[0]])
            bin_widths      = np.diff(edges)
            asimov_vals     = asimov_vals / bin_widths
            asimov_stat_err = np.sqrt(asimov_vals * bin_widths) / bin_widths

            print(f"[DEBUG] Asimov: min={np.min(asimov_vals):.2f}, "
                  f"max={np.max(asimov_vals):.2f}, "
                  f"sum={np.sum(asimov_vals * bin_widths):.2f}")

            toy_dir  = f"{det}/{var}"
            toy_node = safe_get_dir(f, toy_dir)
            if toy_node is None:
                print(f"[WARN] Missing {toy_dir}")
                continue

            all_keys = list(toy_node.keys())
            print(f"[DEBUG] Found {len(all_keys)} total keys in {toy_dir}")
            print(f"[DEBUG] First 5 keys: {all_keys[:5]}")

            toy_keys = [k for k in all_keys
                        if k.startswith(f"{det}_") and f"{var}_posterior_toy_" in k]

            seen = set()
            toy_keys = [k for k in toy_keys if not (k in seen or seen.add(k))]
            print(f"[DEBUG] {len(toy_keys)} unique posterior_toy histograms")

            if not toy_keys:
                print(f"[WARN] No posterior toy histograms for {toy_dir}")
                continue

            toy_keys = sorted(
                toy_keys,
                key=lambda x: int(re.search(r"(\d+)$", x).group(1))
                              if re.search(r"(\d+)$", x) else 0
            )

            toys = []
            for key in toy_keys:
                try:
                    h  = toy_node[key]
                    bw = np.diff(h.axis().edges())
                    v  = h.values()
                    if len(v) == len(asimov_vals):
                        if not np.allclose(h.axis().edges(), edges):
                            print(f"[WARN] Bin edges don't match for {key}")
                        toys.append(v / bw)
                except Exception as e:
                    print(f"[WARN] Problem reading {toy_dir}/{key}: {e}")

            if not toys:
                print(f"[WARN] No valid toy histograms for {toy_dir}")
                continue

            toys      = np.vstack(toys)
            mean_toys = np.mean(toys, axis=0)
            std_toys  = np.std(toys,  axis=0)

            print(f"[DEBUG] Posterior mean: min={np.min(mean_toys):.2f}, "
                  f"max={np.max(mean_toys):.2f}, "
                  f"sum={np.sum(mean_toys * bin_widths):.2f}")
            print(f"[DEBUG] Number of toys used: {len(toys)}")

            # ----------------------------------------------------------------
            # Extend all arrays so the final bin is rendered by step(..., where='post')
            # Each array gets the last value repeated; edges already has N+1 entries.
            # ----------------------------------------------------------------
            asimov_plot,    _ = extend_for_step(asimov_vals,              edges)
            mean_plot,      _ = extend_for_step(mean_toys,                edges)
            std_upper_plot, _ = extend_for_step(mean_toys + std_toys,     edges)
            std_lower_plot, _ = extend_for_step(mean_toys - std_toys,     edges)
            asimov_stat_plot, _ = extend_for_step(asimov_stat_err,        edges)

            bin_centers = 0.5 * (edges[:-1] + edges[1:])

            # ---- Plot ----
            fig, (ax_top, ax_bot) = plt.subplots(
                2, 1, figsize=(6.5, 5.5),
                gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05},
                sharex=True
            )

            ax_top.set_title(f"{det} — {var}", fontsize=12)
            ax_top.fill_between(edges, std_lower_plot, std_upper_plot,
                                step="post", color="red", alpha=0.3, label="±1σ")
            ax_top.step(edges, mean_plot,    where="post",
                        color="red",   lw=1.0, label="Posterior Mean")
            ax_top.step(edges, asimov_plot,  where="post",
                        color="black", lw=1.0, label="Asimov")
            ax_top.errorbar(bin_centers, asimov_vals, yerr=asimov_stat_err,
                            fmt="none", ecolor="black", elinewidth=0.8,
                            capsize=0, alpha=0.8, label="Asimov stat. error")
            ax_top.set_ylabel("Events / GeV")
            ax_top.legend(frameon=False)
            ax_top.grid(alpha=0.3)
            ax_top.set_xlim(edges[0], edges[-1])
            if det in y_limits:
                ax_top.set_ylim(0, y_limits[det])

            # Ratio panel — extend likewise
            ratio_mean   = np.divide(std_toys, mean_toys,
                                     out=np.zeros_like(std_toys), where=mean_toys != 0)
            ratio_asimov = np.divide(std_toys, asimov_vals,
                                     out=np.zeros_like(std_toys), where=asimov_vals != 0)
            ratio_stat   = np.divide(asimov_stat_err, asimov_vals,
                                     out=np.zeros_like(asimov_stat_err), where=asimov_vals != 0)

            ratio_mean_plot,   _ = extend_for_step(ratio_mean,   edges)
            ratio_asimov_plot, _ = extend_for_step(ratio_asimov, edges)
            ratio_stat_plot,   _ = extend_for_step(ratio_stat,   edges)

            ax_bot.step(edges, ratio_mean_plot,   where="post",
                        color="red",   lw=1.3, label="σ / Posterior Mean")
            ax_bot.step(edges, ratio_asimov_plot, where="post",
                        color="blue",  lw=1.3, label="σ / Asimov")
            ax_bot.step(edges, ratio_stat_plot,   where="post",
                        color="green", lw=1.3, label="Stat / Asimov")
            ax_bot.set_ylabel("Relative σ", fontsize=10)
            ax_bot.set_xlabel(f"{var} (GeV)")
            ax_bot.legend(frameon=False, fontsize=8)
            ax_bot.grid(alpha=0.3)
            if det in ratio_limits:
                ax_bot.set_ylim(0, ratio_limits[det])

            pdf_pages.savefig(fig)
            plt.close(fig)
            print(f"[INFO] Added plot for {det}/{var} ({len(toy_keys)} toys)")

    # ----------------------------------------------------------
    # 2D summary plots
    # ----------------------------------------------------------
    plot_2d_summaries(f, detectors, pdf_pages)

    pdf_pages.close()
    print(f"[INFO] Saved plots to {output_pdf}")


# ------------------------------------------------------------
# Entry point
# ------------------------------------------------------------
if __name__ == "__main__":
    main()