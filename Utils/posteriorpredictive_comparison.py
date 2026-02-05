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

def load_toys(f, det, var):
    """Return mean and std of posterior toys for a given detector and variable."""
    toy_dir = f"{det}/{var}"
    toy_node = safe_get_dir(f, toy_dir)
    if not toy_node:
        return None, None, None
    
    all_keys = list(toy_node.keys())
    toy_keys = [k for k in all_keys if f"{var}_posterior_toy_" in k and k.startswith(f"{det}_")]
    if not toy_keys:
        return None, None, None

    # Sort by trailing number
    toy_keys = sorted(toy_keys, key=lambda x: int(re.search(r"(\d+)$", x).group(1)) if re.search(r"(\d+)$", x) else 0)
    
    toys = []
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
    
    if len(toys) == 0:
        return None, None, None
    
    toys = np.vstack(toys)
    mean = np.mean(toys, axis=0)
    std = np.std(toys, axis=0)
    return mean, std, edges

def load_asimov(f, det, var):
    """Return Asimov histogram for a given detector and variable."""
    asimov_node = safe_get_dir(f, "Asimov")
    if not asimov_node:
        return None, None
    keys = [k for k in asimov_node.keys() if f"{det}_" in k and var in k and "_Asimov" in k]
    if not keys:
        return None, None
    
    try:
        h = asimov_node[keys[0]]
        vals, edges = h.to_numpy(flow=True)
        bin_widths = np.diff(edges)
        vals = vals / bin_widths
        return vals, edges
    except Exception as e:
        print(f"[WARN] Failed to load Asimov for {det}/{var}: {e}")
        return None, None

# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():
    if len(sys.argv) < 3:
        print("Usage: python posteriorpredictive_comparison.py on+off.root on.root [output.pdf]")
        sys.exit(1)

    fileA_path = sys.argv[1]  # On+Off axis
    fileB_path = sys.argv[2]  # On-axis only
    output_pdf = sys.argv[3] if len(sys.argv) > 3 else "PosteriorComparison.pdf"

    # Validate file existence
    for path in [fileA_path, fileB_path]:
        if not pathlib.Path(path).exists():
            print(f"[ERROR] File not found: {path}")
            sys.exit(1)

    # Open ROOT files with error handling
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

    detectors = ["ND", "Other"]
    variables = ["RecoNeutrinoEnergy", "TrueNeutrinoEnergy"]

    pdf = PdfPages(output_pdf)
    plots_created = 0

    for det in detectors:
        for var in variables:
            # Load toys
            mean_onoff, std_onoff, edges_onoff = load_toys(fA, det, var)  # On+Off-axis
            mean_on, std_on, edges_on = load_toys(fB, det, var)  # On-axis only

            if mean_onoff is None or mean_on is None or edges_onoff is None or edges_on is None:
                print(f"[WARN] Skipping {det}/{var}, missing data")
                continue

            if not np.allclose(edges_onoff, edges_on):
                print(f"[WARN] Bin edges differ for {det}/{var}, using On+Off-axis binning")
            edges = edges_onoff

            # -----------------------------
            # Event rate plot
            # -----------------------------
            fig, (ax_top, ax_ratio) = plt.subplots(2, 1, figsize=(7,5), gridspec_kw={"height_ratios":[3,1]}, sharex=True)
            ax_top.set_title(f"{det} — {var}", fontsize=12)

            # On+Off-axis in red
            ax_top.fill_between(edges[:-1], mean_onoff-std_onoff, mean_onoff+std_onoff, 
                               step="post", alpha=0.3, color="red", label="On+Off-axis ±1σ")
            ax_top.plot(edges[:-1], mean_onoff, drawstyle="steps-post", 
                       color="red", lw=1.5, label="On+Off-axis mean")

            # On-axis in blue
            ax_top.fill_between(edges[:-1], mean_on-std_on, mean_on+std_on, 
                               step="post", alpha=0.3, color="blue", label="On-axis ±1σ")
            ax_top.plot(edges[:-1], mean_on, drawstyle="steps-post", 
                       color="blue", lw=1.5, label="On-axis mean")

            ax_top.set_ylabel("Events / bin", fontsize=11)
            ax_top.legend(frameon=False, loc='best')
            ax_top.grid(alpha=0.3)

            # Ratio plot: (On+Off-axis) / (On-axis)
            ratio = np.divide(mean_onoff, mean_on, out=np.ones_like(mean_onoff), where=mean_on!=0)
            ax_ratio.plot(edges[:-1], ratio, drawstyle="steps-post", color="black", lw=1.2)
            ax_ratio.axhline(1.0, color="gray", lw=0.8, ls="--")
            ax_ratio.set_ylabel("(On+Off) / On", fontsize=10)
            ax_ratio.set_xlabel(var, fontsize=11)
            ax_ratio.grid(alpha=0.3)

            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)
            plots_created += 1

            # -----------------------------
            # Relative uncertainty plot
            # -----------------------------
            fig, (ax_top, ax_ratio) = plt.subplots(2, 1, figsize=(7,5), gridspec_kw={"height_ratios":[3,1]}, sharex=True)
            ax_top.set_title(f"{det} — {var} Relative Uncertainty", fontsize=12)

            # Calculate relative uncertainties
            rel_onoff = np.divide(std_onoff, mean_onoff, out=np.zeros_like(std_onoff), where=mean_onoff!=0)
            rel_on = np.divide(std_on, mean_on, out=np.zeros_like(std_on), where=mean_on!=0)

            # On+Off-axis in red
            ax_top.plot(edges[:-1], rel_onoff, drawstyle="steps-post", 
                       color="red", lw=1.5, label="On+Off-axis")
            # On-axis in blue
            ax_top.plot(edges[:-1], rel_on, drawstyle="steps-post", 
                       color="blue", lw=1.5, label="On-axis")
            
            ax_top.set_ylabel("σ / Mean", fontsize=11)
            ax_top.legend(frameon=False, loc='best')
            ax_top.grid(alpha=0.3)

            # Ratio of relative uncertainties: (On+Off-axis) / (On-axis)
            rel_ratio = np.divide(rel_onoff, rel_on, out=np.ones_like(rel_onoff), where=rel_on!=0)
            ax_ratio.plot(edges[:-1], rel_ratio, drawstyle="steps-post", color="black", lw=1.2)
            ax_ratio.axhline(1.0, color="gray", lw=0.8, ls="--")
            ax_ratio.set_ylabel("(On+Off) / On", fontsize=10)
            ax_ratio.set_xlabel(var, fontsize=11)
            ax_ratio.grid(alpha=0.3)

            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)
            plots_created += 1

            print(f"[INFO] Plots added for {det}/{var}")

    pdf.close()
    
    if plots_created > 0:
        print(f"[INFO] Saved {plots_created} plots to {output_pdf}")
    else:
        print(f"[WARN] No plots created - check that your ROOT files contain the expected data")

# ------------------------------------------------------------
# Entry point
# ------------------------------------------------------------
if __name__ == "__main__":
    main()