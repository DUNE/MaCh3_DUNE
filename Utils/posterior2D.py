#!/usr/bin/env python3

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


# ============================================================
# HELPERS
# ============================================================

def safe_get(f, path):
    try:
        return f[path]
    except Exception:
        return None


def is_2d(h):
    return len(h.axes) == 2


def get_1d(h):
    return h.values(), h.axes[0].edges()


def get_2d(h):
    return h.values(), h.axes[0].edges(), h.axes[1].edges()


def extend(v):
    return np.append(v, v[-1])


def mask_empty(v):
    v = v.astype(float)
    v[v == 0] = np.nan
    return v


# ============================================================
# 1D POSTERIOR PREDICTIVE
# ============================================================

def plot_1d(pdf, asimov, mean, std, edges, title, xlabel, toys=None):

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(7, 5),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05}
    )

    e = edges

    m = extend(mean)
    u = extend(mean + std)
    l = extend(mean - std)
    a = extend(asimov)

    # top
    # Plot a sample of individual toys to see the spread (plot first so they're behind)
    if toys is not None:
        n_toys_to_plot = min(5000, len(toys))  # Plot up to 50 toys
        for toy in toys[::max(1, len(toys)//n_toys_to_plot)]:
            ax1.step(e, extend(toy), where="post", alpha=0.15, color="gray", linewidth=0.8)
    
    ax1.fill_between(e, l, u, step="post", alpha=0.3, label="±1σ")
    ax1.step(e, m, where="post", label="Posterior Mean", linewidth=2, color="blue")
    ax1.step(e, a, where="post", label="Asimov", linewidth=2, color="orange")

    ax1.set_ylabel("Events / GeV")
    ax1.set_title(title)
    ax1.legend()
    ax1.grid(alpha=0.3)

    # ratio
    ratio = np.divide(std, mean, out=np.zeros_like(std), where=mean != 0)
    r = extend(ratio)

    ax2.step(e, r, where="post", color="red")
    ax2.set_ylabel("σ / μ")
    ax2.set_xlabel(xlabel)
    ax2.grid(alpha=0.3)

    pdf.savefig(fig)
    plt.close(fig)


# ============================================================
# 2D PLOTS
# ============================================================

def plot_2d(ax, v, xe, ye, title, cmap="viridis", norm=None):

    mesh = ax.pcolormesh(
        xe, ye, v.T,
        shading="auto",
        cmap=cmap,
        norm=norm
    )

    ax.set_title(title)
    ax.set_xlabel("Reconstructed Neutrino Energy (GeV)")
    ax.set_ylabel("Lepton Energy (GeV)")

    return mesh


def save_2d(pdf, v, xe, ye, title, cbar_label,
            cmap="viridis", norm=None):

    fig, ax = plt.subplots(figsize=(9, 6))
    mesh = plot_2d(ax, v, xe, ye, title, cmap, norm)

    cb = fig.colorbar(mesh, ax=ax)
    cb.set_label(cbar_label)

    pdf.savefig(fig)
    plt.close(fig)


# ============================================================
# 2D SUMMARY
# ============================================================

def plot_2d_summaries(f, detectors, pdf):

    print("[INFO] 2D summaries...")

    for det in detectors:

        node = safe_get(f, f"{det}/2D/summary")
        if node is None:
            continue

        keys = [k.split(";")[0] for k in node.keys()]

        samples = set()
        for k in keys:
            if "_2D_posteriorMean" in k:
                samples.add(k.replace("_2D_posteriorMean", ""))

        for s in sorted(samples):

            try:
                h_mean = node[f"{s}_2D_posteriorMean"]
                h_post = node[f"{s}_2D_posteriorErr"]
                h_stat = node[f"{s}_2D_statErr"]
            except Exception:
                continue

            mean, xe, ye = get_2d(h_mean)
            post = h_post.values()
            stat = h_stat.values()

            # mask empty bins
            mean = mask_empty(mean)
            post = mask_empty(post)
            stat = mask_empty(stat)

            rel_post = np.where(mean > 0, post / mean, np.nan)
            rel_stat = np.where(mean > 0, stat / mean, np.nan)
            ratio = np.where(stat > 0, post / stat, np.nan)

            vmax = np.nanpercentile(ratio, 95) if np.isfinite(ratio).any() else 2
            norm = mcolors.TwoSlopeNorm(vmin=0, vcenter=1, vmax=max(2, vmax))

            save_2d(pdf, mean, xe, ye,
                    f"{s} Posterior Mean", "Events")

            save_2d(pdf, rel_post, xe, ye,
                    f"{s} Posterior σ / N", "σ/N")

            save_2d(pdf, rel_stat, xe, ye,
                    f"{s} Stat σ / N", "σ/N")

            save_2d(pdf, ratio, xe, ye,
                    f"{s} Posterior / Stat", "σ_post / σ_stat",
                    cmap="RdBu_r", norm=norm)


# ============================================================
# MAIN
# ============================================================

def main():

    if len(sys.argv) < 2:
        print("Usage: script file.root [output.pdf]")
        sys.exit(1)

    f = uproot.open(sys.argv[1])
    out = sys.argv[2] if len(sys.argv) > 2 else "plots.pdf"

    detectors = [d for d in ["ND", "FD"] if safe_get(f, d)]
    print("[INFO] Found detectors:", detectors)

    pdf = PdfPages(out)

    variables = ["TrueNeutrinoEnergy","RecoNeutrinoEnergy", "Enubias"]

    # ========================================================
    # 1D LOOP
    # ========================================================
    for det in detectors:
        for var in variables:

            asimov_dir = safe_get(f, "Asimov")
            if asimov_dir is None:
                continue

            # find correct 1D hist only
            asimov_keys = [
                k for k in asimov_dir.keys()
                if f"{det}_" in k and var in k and "_Asimov" in k
            ]

            if not asimov_keys:
                continue

            h_asimov = asimov_dir[asimov_keys[0]]

            if is_2d(h_asimov):
                continue

            asimov_vals, edges = get_1d(h_asimov)
            bin_widths = np.diff(edges)
            asimov_vals = asimov_vals / bin_widths
            toy_dir = safe_get(f, det)
            if toy_dir is None:
                continue

            toy_keys = [
                k for k in toy_dir.keys()
                if f"{var}_posterior_toy_" in k
            ]

            toys = []

            for k in toy_keys:
                try:
                    h = toy_dir[k]

                    if is_2d(h):
                        continue

                    v, e = get_1d(h)

                    if not np.allclose(e, edges):
                        continue

                    toys.append(v / np.diff(e))

                except Exception:
                    continue

            if not toys:
                continue

            toys = np.vstack(toys)
            mean = np.mean(toys, axis=0)
            std = np.std(toys, axis=0)

            plot_1d(pdf,
                    asimov_vals,
                    mean,
                    std,
                    edges,
                    f"{det} — {var}",
                    f"{var} (GeV)",
                    toys=toys)

            print(f"[INFO] 1D done: {det}/{var}")

    # ========================================================
    # 2D SUMMARY
    # ========================================================
    plot_2d_summaries(f, detectors, pdf)

    pdf.close()
    print("[INFO] Saved to", out)


if __name__ == "__main__":
    main()