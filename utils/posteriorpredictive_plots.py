#!/usr/bin/env python3
"""
posteriorpredictive_plots.py
----------------------------------------
Plots Asimov and posterior-toy spectra from a ROOT output file.

Usage:
    python utils/posteriorpredictive_plots.py Posteriorpredictive_out.root [output.pdf]

Requirements:
    pip install uproot numpy matplotlib
"""

import sys
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import uproot
import pathlib

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

# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():
    if len(sys.argv) < 2:
        print("Usage: python utils/posteriorpredictive_plots.py Posteriorpredictive_out.root [output.pdf]")
        sys.exit(1)

    root_file = sys.argv[1]
    output_pdf = sys.argv[2] if len(sys.argv) > 2 else "SpectraPlots.pdf"

    print(f"[INFO] Opening ROOT file: {root_file}")
    f = uproot.open(root_file)

    detectors = ["ND", "Other"]
    variables = ["RecoNeutrinoEnergy", "TrueNeutrinoEnergy"]

    pdf_pages = PdfPages(output_pdf)

    for det in detectors:
        for var in variables:
            # --------------------------------------------------------
            # Asimov histograms
            # --------------------------------------------------------
            asimov_dir = "Asimov"
            asimov_node = safe_get_dir(f, asimov_dir)
            if not asimov_node:
                print(f"[WARN] Missing {asimov_dir}")
                continue

            # Match Asimov histogram by detector & variable
            asimov_keys = [
                k for k in asimov_node.keys()
                if f"{det}_" in k and var in k and "_Asimov" in k
            ]
            if not asimov_keys:
                print(f"[WARN] No Asimov histograms found for {det}/{var}")
                continue
            h_asimov = asimov_node[asimov_keys[0]]
            asimov_vals, edges = get_hist_values(h_asimov)


            bin_widths = np.diff(edges)
            asimov_vals = asimov_vals / bin_widths


            # --------------------------------------------------------
            # Posterior toy histograms
            # --------------------------------------------------------
            if var == "RecoNeutrinoEnergy":
                # Only use the subdirectory for Reco
                toy_dir = f"{det}/{var}"
                toy_node = safe_get_dir(f, toy_dir)
                if not toy_node:
                    print(f"[WARN] Missing {toy_dir}")
                    continue
                toy_keys = [k for k in toy_node.keys() if "_posterior_toy_" in k]
            else:
                # TrueNeutrinoEnergy: keep existing logic
                toy_dir = f"{det}/{var}"
                toy_node = safe_get_dir(f, toy_dir)
                if not toy_node:
                    print(f"[WARN] Missing {toy_dir}")
                    continue
                toy_keys = [k for k in toy_node.keys() if "_posterior_toy_" in k or var in k]

            if not toy_keys:
                print(f"[WARN] No posterior toy histograms for {toy_dir}")
                continue

            # Sort toy histograms by index
            toy_keys = sorted(
                toy_keys,
                key=lambda x: int(re.search(r"(\d+)$", x).group(1)) if re.search(r"(\d+)$", x) else 0
            )

            toys = []
            for key in toy_keys:
                try:
                    vals = toy_node[key].values()
                    if len(vals) == len(asimov_vals):
                        toys.append(vals / bin_widths)
                except Exception:
                    print(f"[WARN] Problem reading {toy_dir}/{key}")

            if not toys:
                print(f"[WARN] No valid toy histograms for {toy_dir}")
                continue

            toys = np.vstack(toys)
            mean_toys = np.mean(toys, axis=0)
            std_toys = np.std(toys, axis=0)

            # --------------------------------------------------------
            # Plot: step-style histograms + ratio panel
            # --------------------------------------------------------
            fig, (ax_top, ax_bottom) = plt.subplots(
                2, 1,
                figsize=(6.5, 5.5),
                gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05},
                sharex=True
            )

            # --- Top pad ---
            ax_top.set_title(f"{det} — {var}", fontsize=12)
            ax_top.fill_between(
                edges[:-1],
                mean_toys - std_toys,
                mean_toys + std_toys,
                step="post",
                color="red",
                alpha=0.3,
                label="±1σ"
            )
            ax_top.step(edges[:-1], mean_toys, where="post", color="red", lw=0.5, label="Posterior Mean")
            ax_top.step(edges[:-1], asimov_vals, where="post", color="black", lw=0.5, label="Asimov")
            ax_top.set_ylabel("Events / bin")
            ax_top.legend(frameon=False)
            ax_top.grid(alpha=0.3)

            # --- Bottom pad: relative uncertainty ---
            ratio = np.divide(std_toys, mean_toys, out=np.zeros_like(std_toys), where=mean_toys != 0)
            ratio2 = np.divide(std_toys, asimov_vals, out=np.zeros_like(std_toys), where=asimov_vals != 0)
            ax_bottom.step(edges[:-1], ratio, where="post", color="red", lw=1.3)
            ax_bottom.step(edges[:-1], ratio2, where="post", color="blue", lw=1.3)
            ax_bottom.set_ylabel("σ / Mean", fontsize=10)
            ax_bottom.set_xlabel(var)
            ax_bottom.grid(alpha=0.3)
            ax_bottom.set_ylim(0, max(ratio) * 1.4)

            plt.tight_layout()
            pdf_pages.savefig()
            plt.close()

            print(f"[INFO] Added plot for {det}/{var} ({len(toy_keys)} toys)")

    pdf_pages.close()
    print(f"[INFO] Saved plots to {output_pdf}")

# ------------------------------------------------------------
# Entry point
# ------------------------------------------------------------
if __name__ == "__main__":
    main()
