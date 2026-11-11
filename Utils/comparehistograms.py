#!/usr/bin/env python3
"""
Overlay two histograms from two ROOT files, divided by bin width.

Usage
-----
    python compare_histograms.py file1.root file2.root [output.pdf]
"""

import sys
import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HIST_NAME = "hOffAxis0m_numuCC_numode"


def load_hist(path):
    f = uproot.open(path)
    h = f[HIST_NAME]
    vals, edges = h.values(), h.axes[0].edges()
    widths = np.diff(edges)
    centres = 0.5 * (edges[:-1] + edges[1:])
    return vals / widths, edges, centres


def main():
    if len(sys.argv) < 3:
        print("Usage: compare_histograms.py file1.root file2.root [output.pdf]")
        sys.exit(1)

    path1 = sys.argv[1]
    path2 = sys.argv[2]
    out_path = sys.argv[3] if len(sys.argv) > 3 else "comparison.pdf"

    vals1, edges1, centres1 = load_hist(path1)
    vals2, edges2, centres2 = load_hist(path2)

    fig, ax = plt.subplots(figsize=(8, 5))

    ax.step(edges1, np.append(vals1, vals1[-1]), where="post",
            label=path1, linewidth=2, color="steelblue")
    ax.step(edges2, np.append(vals2, vals2[-1]), where="post",
            label=path2, linewidth=2, color="darkorange")

    ax.set_xlabel("Reconstructed Neutrino Energy (GeV)")
    ax.set_ylabel("Events / GeV")
    ax.set_title(HIST_NAME)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(out_path)
    print(f"Saved → {out_path}")


if __name__ == "__main__":
    main()