# #!/usr/bin/env python3
# """
# posteriorpredictive_plots.py
# ----------------------------------------
# Plots Asimov and posterior-toy spectra from a ROOT output file.

# Usage:
#     python utils/posteriorpredictive_plots.py Posteriorpredictive_out.root [output.pdf]

# Requirements:
#     pip install uproot numpy matplotlib
# """

# import sys
# import re
# import numpy as np

# import uproot
# import pathlib
# import matplotlib
# matplotlib.use("Agg")

# import matplotlib.pyplot as plt
# from matplotlib.backends.backend_pdf import PdfPages

# import ROOT
# ROOT.gROOT.SetBatch(True)

# import numpy as np

# # ------------------------------------------------------------
# # Helpers
# # ------------------------------------------------------------

# def get_hist_values(obj):
#     """Return (values, edges) from an uproot TH1 object."""
#     return obj.values(), obj.axis().edges()

# def safe_get_dir(f, path):
#     """Safely get a directory from uproot file."""
#     try:
#         return f[path]
#     except Exception:
#         return None

# # ------------------------------------------------------------
# # Main
# # ------------------------------------------------------------

# def main():
#     if len(sys.argv) < 2:
#         print("Usage: python utils/posteriorpredictive_plots.py Posteriorpredictive_out.root [output.pdf]")
#         sys.exit(1)

#     root_file = sys.argv[1]
#     output_pdf = sys.argv[2] if len(sys.argv) > 2 else "SpectraPlots.pdf"

#     print(f"[INFO] Opening ROOT file: {root_file}")
#     f = uproot.open(root_file)

#     detectors = ["ND", "Other"]
#     variables = ["RecoNeutrinoEnergy", "TrueNeutrinoEnergy"]

#     pdf_pages = PdfPages(output_pdf)

#     for det in detectors:
#         for var in variables:
#             # --------------------------------------------------------
#             # Asimov histograms
#             # --------------------------------------------------------
#             asimov_dir = "Asimov"
#             asimov_node = safe_get_dir(f, asimov_dir)
#             if not asimov_node:
#                 print(f"[WARN] Missing {asimov_dir}")
#                 continue

#             # Match Asimov histogram by detector & variable
#             asimov_keys = [
#                 k for k in asimov_node.keys()
#                 if f"{det}_" in k and var in k and "_Asimov" in k
#             ]
#             if not asimov_keys:
#                 print(f"[WARN] No Asimov histograms found for {det}/{var}")
#                 continue
#             h_asimov = asimov_node[asimov_keys[0]]
#             asimov_vals, edges = get_hist_values(h_asimov)


#             bin_widths = np.diff(edges)
#             asimov_vals = asimov_vals / bin_widths


#             # --------------------------------------------------------
#             # Posterior toy histograms
#             # --------------------------------------------------------
#             if var == "RecoNeutrinoEnergy":
#                 # Only use the subdirectory for Reco
#                 toy_dir = f"{det}/{var}"
#                 print("In directoy {}".format(toy_dir))
#                 toy_node = safe_get_dir(f, toy_dir)
#                 if not toy_node:
#                     print(f"[WARN] Missing {toy_dir}")
#                     continue
#                 toy_keys = [k for k in toy_node.keys() if "_OffAxisND_posterior_toy_" in k]
#             if var == "TrueNeutrinoEnergy":
#                 # Only use the subdirectory for Reco
#                 toy_dir = f"{det}/{var}"
#                 print("In directoy {}".format(toy_dir))
#                 toy_node = safe_get_dir(f, toy_dir)
#                 if not toy_node:
#                     print(f"[WARN] Missing {toy_dir}")
#                     continue
#                 toy_keys = [k for k in toy_node.keys() if "_OffAxisND_posterior_toy_" in k]
#             else:
#                 print("found another direcotry...")
#                 # # TrueNeutrinoEnergy: keep existing logic
#                 # toy_dir = f"{det}/{var}"
#                 # toy_node = safe_get_dir(f, toy_dir)
#                 # if not toy_node:
#                 #     print(f"[WARN] Missing {toy_dir}")
#                 #     continue
#                 # toy_keys = [k for k in toy_node.keys() if "_OffAxisND_posterior_toy_" in k or var in k]

#             if not toy_keys:
#                 print(f"[WARN] No posterior toy histograms for {toy_dir}")
#                 continue

#             # Sort toy histograms by index
#             toy_keys = sorted(
#                 toy_keys,
#                 key=lambda x: int(re.search(r"(\d+)$", x).group(1)) if re.search(r"(\d+)$", x) else 0
#             )

#             toys = []
#             for key in toy_keys:
#                 try:
#                     vals = toy_node[key].values()
#                     if len(vals) == len(asimov_vals):
#                         toys.append(vals / bin_widths)
#                 except Exception:
#                     print(f"[WARN] Problem reading {toy_dir}/{key}")

#             if not toys:
#                 print(f"[WARN] No valid toy histograms for {toy_dir}")
#                 continue

#             toys = np.vstack(toys)
#             mean_toys = np.mean(toys, axis=0)
#             std_toys = np.std(toys, axis=0)

#             # --------------------------------------------------------
#             # Plot: step-style histograms + ratio panel
#             # --------------------------------------------------------
#             fig, (ax_top, ax_bottom) = plt.subplots(
#                 2, 1,
#                 figsize=(6.5, 5.5),
#                 gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05},
#                 sharex=True
#             )
           

#             # --- Top pad ---
#             ax_top.set_title(f"{det} — {var}", fontsize=12)
#             ax_top.fill_between(
#                 edges[:-1],
#                 mean_toys - std_toys,
#                 mean_toys + std_toys,
#                 step="post",
#                 color="red",
#                 alpha=0.3,
#                 label="±1σ"
#             )
#             ax_top.step(edges[:-1], mean_toys, where="post", color="red", lw=0.5, label="Posterior Mean")
#             ax_top.step(edges[:-1], asimov_vals, where="post", color="black", lw=0.5, label="Asimov")
#             ax_top.set_ylabel("Events / bin")
#             ax_top.legend(frameon=False)
#             ax_top.grid(alpha=0.3)

#             # --- Bottom pad: relative uncertainty ---
#             ratio = np.divide(std_toys, mean_toys, out=np.zeros_like(std_toys), where=mean_toys != 0)
#             ratio2 = np.divide(std_toys, asimov_vals, out=np.zeros_like(std_toys), where=asimov_vals != 0)
#             ax_bottom.step(edges[:-1], ratio, where="post", color="red", lw=1.3)
#             ax_bottom.step(edges[:-1], ratio2, where="post", color="blue", lw=1.3)
#             ax_bottom.set_ylabel("σ / Mean", fontsize=10)
#             ax_bottom.set_xlabel(var)
#             ax_bottom.grid(alpha=0.3)
#             ax_bottom.set_ylim(0, max(ratio) * 1.4)

            
#             pdf_pages.savefig(fig)
#             plt.close(fig)


#             print(f"[INFO] Added plot for {det}/{var} ({len(toy_keys)} toys)")

#     pdf_pages.close()
#     print(f"[INFO] Saved plots to {output_pdf}")

# # ------------------------------------------------------------
# # Entry point
# # ------------------------------------------------------------
# if __name__ == "__main__":
#     main()




# #!/usr/bin/env python3
# """
# posteriorpredictive_plots.py
# ----------------------------------------
# Plots Asimov and posterior-toy spectra from a ROOT output file.

# Usage:
#     python utils/posteriorpredictive_plots.py Posteriorpredictive_out.root [output.pdf]

# Requirements:
#     pip install uproot numpy matplotlib
# """

# import sys
# import re
# import numpy as np
# import uproot
# import pathlib
# import matplotlib
# matplotlib.use("Agg")

# import matplotlib.pyplot as plt
# from matplotlib.backends.backend_pdf import PdfPages

# import ROOT
# ROOT.gROOT.SetBatch(True)

# # ------------------------------------------------------------
# # Helpers
# # ------------------------------------------------------------

# def get_hist_values(obj):
#     """Return (values, edges) from an uproot TH1 object."""
#     return obj.values(), obj.axis().edges()

# def safe_get_dir(f, path):
#     """Safely get a directory from uproot file."""
#     try:
#         return f[path]
#     except Exception:
#         return None

# # ------------------------------------------------------------
# # Main
# # ------------------------------------------------------------

# def main():
#     if len(sys.argv) < 2:
#         print("Usage: python utils/posteriorpredictive_plots.py Posteriorpredictive_out.root [output.pdf]")
#         sys.exit(1)

#     root_file = sys.argv[1]
#     output_pdf = sys.argv[2] if len(sys.argv) > 2 else "SpectraPlots.pdf"

#     # Check if input file exists
#     if not pathlib.Path(root_file).exists():
#         print(f"[ERROR] Input file not found: {root_file}")
#         sys.exit(1)

#     print(f"[INFO] Opening ROOT file: {root_file}")
#     f = uproot.open(root_file)

#     detectors = ["ND", "Other"]
#     variables = ["RecoNeutrinoEnergy", "TrueNeutrinoEnergy"]

#     pdf_pages = PdfPages(output_pdf)

#     # First pass: determine y-axis limits for each detector
#     y_limits = {}
#     ratio_limits = {}
#     for det in detectors:
#         max_y = 0
#         max_ratio = 0
#         for var in variables:
#             # Get Asimov
#             asimov_dir = "Asimov"
#             asimov_node = safe_get_dir(f, asimov_dir)
#             if not asimov_node:
#                 continue
            
#             asimov_keys = [
#                 k for k in asimov_node.keys()
#                 if f"{det}_" in k and var in k and "_Asimov" in k
#             ]
#             if not asimov_keys:
#                 continue
            
#             h_asimov = asimov_node[asimov_keys[0]]
#             asimov_vals, edges = get_hist_values(h_asimov)
#             bin_widths = np.diff(edges)
#             asimov_vals = asimov_vals / bin_widths
            
#             # Get toys
#             toy_dir = f"{det}/{var}"
#             toy_node = safe_get_dir(f, toy_dir)
#             if not toy_node:
#                 continue
            
#             all_keys = list(toy_node.keys())
#             toy_keys = [k for k in all_keys if f"{det}_OffAxisND_{var}_posterior_toy_" in k]
            
#             if not toy_keys:
#                 continue
            
#             toys = []
#             for key in toy_keys:
#                 try:
#                     toy_hist = toy_node[key]
#                     vals = toy_hist.values()
#                     toy_edges = toy_hist.axis().edges()
#                     toy_bin_widths = np.diff(toy_edges)
#                     if len(vals) == len(asimov_vals):
#                         toys.append(vals / toy_bin_widths)
#                 except Exception:
#                     pass
            
#             if toys:
#                 toys = np.vstack(toys)
#                 mean_toys = np.mean(toys, axis=0)
#                 std_toys = np.std(toys, axis=0)
                
#                 # Track maximum y value
#                 max_y = max(max_y, np.max(asimov_vals), np.max(mean_toys + std_toys))
                
#                 # Track maximum ratio value
#                 ratio_mean = np.divide(std_toys, mean_toys, out=np.zeros_like(std_toys), where=mean_toys != 0)
#                 ratio_asimov = np.divide(std_toys, asimov_vals, out=np.zeros_like(std_toys), where=asimov_vals != 0)
#                 max_ratio = max(max_ratio, np.max(ratio_mean), np.max(ratio_asimov))
        
#         y_limits[det] = max_y * 1.1  # Add 10% padding
#         ratio_limits[det] = max_ratio * 1.4  # Add 40% padding

#     # Second pass: create plots with consistent y-axis
#     for det in detectors:
#         for var in variables:
#             # --------------------------------------------------------
#             # Asimov histograms
#             # --------------------------------------------------------
#             asimov_dir = "Asimov"
#             asimov_node = safe_get_dir(f, asimov_dir)
#             if not asimov_node:
#                 print(f"[WARN] Missing {asimov_dir}")
#                 continue

#             # Match Asimov histogram by detector & variable
#             asimov_keys = [
#                 k for k in asimov_node.keys()
#                 if f"{det}_" in k and var in k and "_Asimov" in k
#             ]
#             if not asimov_keys:
#                 print(f"[WARN] No Asimov histograms found for {det}/{var}")
#                 print(f"[DEBUG] Available keys in {asimov_dir}: {list(asimov_node.keys())[:5]}...")
#                 continue
#             print(f"[DEBUG] Using Asimov histogram: {asimov_keys[0]}")
#             h_asimov = asimov_node[asimov_keys[0]]
#             asimov_vals, edges = get_hist_values(h_asimov)

#             bin_widths = np.diff(edges)
#             asimov_vals = asimov_vals / bin_widths
            
#             print(f"[DEBUG] Asimov: min={np.min(asimov_vals):.2f}, max={np.max(asimov_vals):.2f}, sum={np.sum(asimov_vals * bin_widths):.2f}")

#             # --------------------------------------------------------
#             # Posterior toy histograms
#             # --------------------------------------------------------
#             toy_dir = f"{det}/{var}"
#             print(f"[INFO] In directory {toy_dir}")
#             toy_node = safe_get_dir(f, toy_dir)
#             if not toy_node:
#                 print(f"[WARN] Missing {toy_dir}")
#                 continue
            
#             # Find posterior toy histograms in this directory
#             all_keys = list(toy_node.keys())
#             print(f"[DEBUG] Found {len(all_keys)} total keys in {toy_dir}")
#             print(f"[DEBUG] First 5 keys: {all_keys[:5]}")
            
#             # Search for histograms matching pattern: ND_OffAxisND_{var}_posterior_toy_XXX
#             toy_keys = [k for k in all_keys if f"{det}_OffAxisND_{var}_posterior_toy_" in k]
#             print(f"[DEBUG] Found {len(toy_keys)} posterior_toy histograms matching pattern")
            
#             # Remove duplicates by name (keep first occurrence)
#             seen = set()
#             unique_toy_keys = []
#             for k in toy_keys:
#                 if k not in seen:
#                     seen.add(k)
#                     unique_toy_keys.append(k)
#             toy_keys = unique_toy_keys
#             print(f"[DEBUG] After removing duplicates: {len(toy_keys)} unique histograms")

#             if not toy_keys:
#                 print(f"[WARN] No posterior toy histograms for {toy_dir}")
#                 continue

#             # Sort toy histograms by index
#             toy_keys = sorted(
#                 toy_keys,
#                 key=lambda x: int(re.search(r"(\d+)$", x).group(1)) if re.search(r"(\d+)$", x) else 0
#             )

#             toys = []
#             for key in toy_keys:
#                 try:
#                     toy_hist = toy_node[key]
#                     vals = toy_hist.values()
#                     toy_edges = toy_hist.axis().edges()
#                     toy_bin_widths = np.diff(toy_edges)
                    
#                     # Check if binning matches
#                     if len(vals) == len(asimov_vals):
#                         if not np.allclose(toy_edges, edges):
#                             print(f"[WARN] Bin edges don't match for {key}")
#                         toys.append(vals / toy_bin_widths)
#                 except Exception as e:
#                     print(f"[WARN] Problem reading {toy_dir}/{key}: {e}")

#             if not toys:
#                 print(f"[WARN] No valid toy histograms for {toy_dir}")
#                 continue

#             toys = np.vstack(toys)
#             mean_toys = np.mean(toys, axis=0)
#             std_toys = np.std(toys, axis=0)
            
#             print(f"[DEBUG] Posterior mean: min={np.min(mean_toys):.2f}, max={np.max(mean_toys):.2f}, sum={np.sum(mean_toys * bin_widths):.2f}")
#             print(f"[DEBUG] Number of toys used: {len(toys)}")

#             # --------------------------------------------------------
#             # Plot: step-style histograms + ratio panel
#             # --------------------------------------------------------
#             fig, (ax_top, ax_bottom) = plt.subplots(
#                 2, 1,
#                 figsize=(6.5, 5.5),
#                 gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05},
#                 sharex=True
#             )

#             # --- Top pad ---
#             ax_top.set_title(f"{det} — {var}", fontsize=12)
#             ax_top.fill_between(
#                 edges[:-1],
#                 mean_toys - std_toys,
#                 mean_toys + std_toys,
#                 step="post",
#                 color="red",
#                 alpha=0.3,
#                 label="±1σ"
#             )
#             ax_top.step(edges[:-1], mean_toys, where="post", color="red", lw=1.0, label="Posterior Mean")
#             ax_top.step(edges[:-1], asimov_vals, where="post", color="black", lw=1.0, label="Asimov")
#             ax_top.set_ylabel("Events / bin")
#             ax_top.legend(frameon=False)
#             ax_top.grid(alpha=0.3)
            
#             # Set consistent y-axis limit for this detector
#             if det in y_limits:
#                 ax_top.set_ylim(0, y_limits[det])

#             # --- Bottom pad: relative uncertainty ---
#             ratio_mean = np.divide(std_toys, mean_toys, out=np.zeros_like(std_toys), where=mean_toys != 0)
#             ratio_asimov = np.divide(std_toys, asimov_vals, out=np.zeros_like(std_toys), where=asimov_vals != 0)
            
#             ax_bottom.step(edges[:-1], ratio_mean, where="post", color="red", lw=1.3, label="σ / Posterior Mean")
#             ax_bottom.step(edges[:-1], ratio_asimov, where="post", color="blue", lw=1.3, label="σ / Asimov")
#             ax_bottom.set_ylabel("Relative σ", fontsize=10)
#             ax_bottom.set_xlabel(var)
#             ax_bottom.legend(frameon=False, fontsize=8)
#             ax_bottom.grid(alpha=0.3)
            
#             # Set consistent ratio y-axis limit for this detector
#             if det in ratio_limits:
#                 ax_bottom.set_ylim(0, ratio_limits[det])

#             pdf_pages.savefig(fig)
#             plt.close(fig)

#             print(f"[INFO] Added plot for {det}/{var} ({len(toy_keys)} toys)")

#     pdf_pages.close()
#     print(f"[INFO] Saved plots to {output_pdf}")

# # ------------------------------------------------------------
# # Entry point
# # ------------------------------------------------------------
# if __name__ == "__main__":
#     main()

# # # ------------------------------------------------------------
# # # Helpers
# # # ------------------------------------------------------------

# # def get_hist_values(obj):
# #     """Return (values, edges) from an uproot TH1 object."""
# #     return obj.values(), obj.axis().edges()

# # def safe_get_dir(f, path):
# #     """Safely get a directory from uproot file."""
# #     try:
# #         return f[path]
# #     except Exception:
# #         return None

# # # ------------------------------------------------------------
# # # Main
# # # ------------------------------------------------------------

# # def main():
# #     if len(sys.argv) < 2:
# #         print("Usage: python utils/posteriorpredictive_plots.py Posteriorpredictive_out.root [output.pdf]")
# #         sys.exit(1)

# #     root_file = sys.argv[1]
# #     output_pdf = sys.argv[2] if len(sys.argv) > 2 else "SpectraPlots.pdf"

# #     # Check if input file exists
# #     if not pathlib.Path(root_file).exists():
# #         print(f"[ERROR] Input file not found: {root_file}")
# #         sys.exit(1)

# #     print(f"[INFO] Opening ROOT file: {root_file}")
# #     f = uproot.open(root_file)

# #     detectors = ["ND", "Other"]
# #     variables = ["RecoNeutrinoEnergy", "TrueNeutrinoEnergy"]

# #     pdf_pages = PdfPages(output_pdf)

# #     # First pass: determine y-axis limits for each detector
# #     y_limits = {}
# #     for det in detectors:
# #         max_y = 0
# #         for var in variables:
# #             # Get Asimov
# #             asimov_dir = "Asimov"
# #             asimov_node = safe_get_dir(f, asimov_dir)
# #             if not asimov_node:
# #                 continue
            
# #             asimov_keys = [
# #                 k for k in asimov_node.keys()
# #                 if f"{det}_" in k and var in k and "_Asimov" in k
# #             ]
# #             if not asimov_keys:
# #                 continue
            
# #             h_asimov = asimov_node[asimov_keys[0]]
# #             asimov_vals, edges = get_hist_values(h_asimov)
# #             bin_widths = np.diff(edges)
# #             asimov_vals = asimov_vals / bin_widths
            
# #             # Get toys
# #             toy_dir = f"{det}/{var}"
# #             toy_node = safe_get_dir(f, toy_dir)
# #             if not toy_node:
# #                 continue
            
# #             all_keys = list(toy_node.keys())
# #             toy_keys = [k for k in all_keys if f"{det}_OffAxisND_{var}_posterior_toy_" in k]
            
# #             if not toy_keys:
# #                 continue
            
# #             toys = []
# #             for key in toy_keys:
# #                 try:
# #                     toy_hist = toy_node[key]
# #                     vals = toy_hist.values()
# #                     toy_edges = toy_hist.axis().edges()
# #                     toy_bin_widths = np.diff(toy_edges)
# #                     if len(vals) == len(asimov_vals):
# #                         toys.append(vals / toy_bin_widths)
# #                 except Exception:
# #                     pass
            
# #             if toys:
# #                 toys = np.vstack(toys)
# #                 mean_toys = np.mean(toys, axis=0)
# #                 std_toys = np.std(toys, axis=0)
                
# #                 # Track maximum y value
# #                 max_y = max(max_y, np.max(asimov_vals), np.max(mean_toys + std_toys))
        
# #         y_limits[det] = max_y * 1.1  # Add 10% padding

# #     # Second pass: create plots with consistent y-axis
# #     for det in detectors:
# #         for var in variables:
# #             # --------------------------------------------------------
# #             # Asimov histograms
# #             # --------------------------------------------------------
# #             asimov_dir = "Asimov"
# #             asimov_node = safe_get_dir(f, asimov_dir)
# #             if not asimov_node:
# #                 print(f"[WARN] Missing {asimov_dir}")
# #                 continue

# #             # Match Asimov histogram by detector & variable
# #             asimov_keys = [
# #                 k for k in asimov_node.keys()
# #                 if f"{det}_" in k and var in k and "_Asimov" in k
# #             ]
# #             if not asimov_keys:
# #                 print(f"[WARN] No Asimov histograms found for {det}/{var}")
# #                 print(f"[DEBUG] Available keys in {asimov_dir}: {list(asimov_node.keys())[:5]}...")
# #                 continue
# #             print(f"[DEBUG] Using Asimov histogram: {asimov_keys[0]}")
# #             h_asimov = asimov_node[asimov_keys[0]]
# #             asimov_vals, edges = get_hist_values(h_asimov)

# #             bin_widths = np.diff(edges)
# #             asimov_vals = asimov_vals / bin_widths
            
# #             print(f"[DEBUG] Asimov: min={np.min(asimov_vals):.2f}, max={np.max(asimov_vals):.2f}, sum={np.sum(asimov_vals * bin_widths):.2f}")

# #             # --------------------------------------------------------
# #             # Posterior toy histograms
# #             # --------------------------------------------------------
# #             toy_dir = f"{det}/{var}"
# #             print(f"[INFO] In directory {toy_dir}")
# #             toy_node = safe_get_dir(f, toy_dir)
# #             if not toy_node:
# #                 print(f"[WARN] Missing {toy_dir}")
# #                 continue
            
# #             # Find posterior toy histograms in this directory
# #             all_keys = list(toy_node.keys())
# #             print(f"[DEBUG] Found {len(all_keys)} total keys in {toy_dir}")
# #             print(f"[DEBUG] First 5 keys: {all_keys[:5]}")
            
# #             # Search for histograms matching pattern: ND_OffAxisND_{var}_posterior_toy_XXX
# #             toy_keys = [k for k in all_keys if f"{det}_OffAxisND_{var}_posterior_toy_" in k]
# #             print(f"[DEBUG] Found {len(toy_keys)} posterior_toy histograms matching pattern")
            
# #             # Remove duplicates by name (keep first occurrence)
# #             seen = set()
# #             unique_toy_keys = []
# #             for k in toy_keys:
# #                 if k not in seen:
# #                     seen.add(k)
# #                     unique_toy_keys.append(k)
# #             toy_keys = unique_toy_keys
# #             print(f"[DEBUG] After removing duplicates: {len(toy_keys)} unique histograms")

# #             if not toy_keys:
# #                 print(f"[WARN] No posterior toy histograms for {toy_dir}")
# #                 continue

# #             # Sort toy histograms by index
# #             toy_keys = sorted(
# #                 toy_keys,
# #                 key=lambda x: int(re.search(r"(\d+)$", x).group(1)) if re.search(r"(\d+)$", x) else 0
# #             )

# #             toys = []
# #             for key in toy_keys:
# #                 try:
# #                     toy_hist = toy_node[key]
# #                     vals = toy_hist.values()
# #                     toy_edges = toy_hist.axis().edges()
# #                     toy_bin_widths = np.diff(toy_edges)
                    
# #                     # Check if binning matches
# #                     if len(vals) == len(asimov_vals):
# #                         if not np.allclose(toy_edges, edges):
# #                             print(f"[WARN] Bin edges don't match for {key}")
# #                         toys.append(vals / toy_bin_widths)
# #                 except Exception as e:
# #                     print(f"[WARN] Problem reading {toy_dir}/{key}: {e}")

# #             if not toys:
# #                 print(f"[WARN] No valid toy histograms for {toy_dir}")
# #                 continue

# #             toys = np.vstack(toys)
# #             mean_toys = np.mean(toys, axis=0)
# #             std_toys = np.std(toys, axis=0)
            
# #             print(f"[DEBUG] Posterior mean: min={np.min(mean_toys):.2f}, max={np.max(mean_toys):.2f}, sum={np.sum(mean_toys * bin_widths):.2f}")
# #             print(f"[DEBUG] Number of toys used: {len(toys)}")

# #             # --------------------------------------------------------
# #             # Plot: step-style histograms + ratio panel
# #             # --------------------------------------------------------
# #             fig, (ax_top, ax_bottom) = plt.subplots(
# #                 2, 1,
# #                 figsize=(6.5, 5.5),
# #                 gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05},
# #                 sharex=True
# #             )

# #             # --- Top pad ---
# #             ax_top.set_title(f"{det} — {var}", fontsize=12)
# #             ax_top.fill_between(
# #                 edges[:-1],
# #                 mean_toys - std_toys,
# #                 mean_toys + std_toys,
# #                 step="post",
# #                 color="red",
# #                 alpha=0.3,
# #                 label="±1σ"
# #             )
# #             ax_top.step(edges[:-1], mean_toys, where="post", color="red", lw=1.0, label="Posterior Mean")
# #             ax_top.step(edges[:-1], asimov_vals, where="post", color="black", lw=1.0, label="Asimov")
# #             ax_top.set_ylabel("Events / bin")
# #             ax_top.legend(frameon=False)
# #             ax_top.grid(alpha=0.3)
            
# #             # Set consistent y-axis limit for this detector
# #             if det in y_limits:
# #                 ax_top.set_ylim(0, y_limits[det])

# #             # --- Bottom pad: relative uncertainty ---
# #             ratio_mean = np.divide(std_toys, mean_toys, out=np.zeros_like(std_toys), where=mean_toys != 0)
#             ratio_asimov = np.divide(std_toys, asimov_vals, out=np.zeros_like(std_toys), where=asimov_vals != 0)
            
#             ax_bottom.step(edges[:-1], ratio_mean, where="post", color="red", lw=1.3, label="σ / Posterior Mean")
#             ax_bottom.step(edges[:-1], ratio_asimov, where="post", color="blue", lw=1.3, label="σ / Asimov")
#             ax_bottom.set_ylabel("Relative σ", fontsize=10)
#             ax_bottom.set_xlabel(var)
#             ax_bottom.legend(frameon=False, fontsize=8)
#             ax_bottom.grid(alpha=0.3)
#             max_ratio = max(np.max(ratio_mean), np.max(ratio_asimov))
#             ax_bottom.set_ylim(0, max_ratio * 1.4)

#             pdf_pages.savefig(fig)
#             plt.close(fig)

#             print(f"[INFO] Added plot for {det}/{var} ({len(toy_keys)} toys)")

#     pdf_pages.close()
#     print(f"[INFO] Saved plots to {output_pdf}")

# # ------------------------------------------------------------
# # Entry point
# # ------------------------------------------------------------
# if __name__ == "__main__":
#     main()

# #     for det in detectors:
# #         for var in variables:
# #             # --------------------------------------------------------
# #             # Asimov histograms
# #             # --------------------------------------------------------
# #             asimov_dir = "Asimov"
# #             asimov_node = safe_get_dir(f, asimov_dir)
# #             if not asimov_node:
# #                 print(f"[WARN] Missing {asimov_dir}")
# #                 continue

# #             # Match Asimov histogram by detector & variable
# #             asimov_keys = [
# #                 k for k in asimov_node.keys()
# #                 if f"{det}_" in k and var in k and "_Asimov" in k
# #             ]
# #             if not asimov_keys:
# #                 print(f"[WARN] No Asimov histograms found for {det}/{var}")
# #                 print(f"[DEBUG] Available keys in {asimov_dir}: {list(asimov_node.keys())[:5]}...")
# #                 continue
# #             print(f"[DEBUG] Using Asimov histogram: {asimov_keys[0]}")
# #             h_asimov = asimov_node[asimov_keys[0]]
# #             asimov_vals, edges = get_hist_values(h_asimov)

# #             bin_widths = np.diff(edges)
# #             asimov_vals = asimov_vals / bin_widths
            
# #             print(f"[DEBUG] Asimov: min={np.min(asimov_vals):.2f}, max={np.max(asimov_vals):.2f}, sum={np.sum(asimov_vals * bin_widths):.2f}")

# #             # --------------------------------------------------------
# #             # Posterior toy histograms
# #             # --------------------------------------------------------
# #             toy_dir = f"{det}/{var}"
# #             print(f"[INFO] In directory {toy_dir}")
# #             toy_node = safe_get_dir(f, toy_dir)
# #             if not toy_node:
# #                 print(f"[WARN] Missing {toy_dir}")
# #                 continue
            
# #             # Find posterior toy histograms in this directory
# #             all_keys = list(toy_node.keys())
# #             print(f"[DEBUG] Found {len(all_keys)} total keys in {toy_dir}")
# #             print(f"[DEBUG] First 5 keys: {all_keys[:5]}")
            
# #             # Search for histograms matching pattern: ND_OffAxisND_{var}_posterior_toy_XXX
# #             toy_keys = [k for k in all_keys if f"{det}_OffAxisND_{var}_posterior_toy_" in k]
# #             print(f"[DEBUG] Found {len(toy_keys)} posterior_toy histograms matching pattern")
            
# #             # Remove duplicates by name (keep first occurrence)
# #             seen = set()
# #             unique_toy_keys = []
# #             for k in toy_keys:
# #                 if k not in seen:
# #                     seen.add(k)
# #                     unique_toy_keys.append(k)
# #             toy_keys = unique_toy_keys
# #             print(f"[DEBUG] After removing duplicates: {len(toy_keys)} unique histograms")

# #             if not toy_keys:
# #                 print(f"[WARN] No posterior toy histograms for {toy_dir}")
# #                 continue

# #             # Sort toy histograms by index
# #             toy_keys = sorted(
# #                 toy_keys,
# #                 key=lambda x: int(re.search(r"(\d+)$", x).group(1)) if re.search(r"(\d+)$", x) else 0
# #             )

# #             toys = []
# #             for key in toy_keys:
# #                 try:
# #                     toy_hist = toy_node[key]
# #                     vals = toy_hist.values()
# #                     toy_edges = toy_hist.axis().edges()
# #                     toy_bin_widths = np.diff(toy_edges)
                    
# #                     # Check if binning matches
# #                     if len(vals) == len(asimov_vals):
# #                         if not np.allclose(toy_edges, edges):
# #                             print(f"[WARN] Bin edges don't match for {key}")
# #                         toys.append(vals / toy_bin_widths)
# #                 except Exception as e:
# #                     print(f"[WARN] Problem reading {toy_dir}/{key}: {e}")

# #             if not toys:
# #                 print(f"[WARN] No valid toy histograms for {toy_dir}")
# #                 continue

# #             toys = np.vstack(toys)
# #             mean_toys = np.mean(toys, axis=0)
# #             std_toys = np.std(toys, axis=0)
            
# #             print(f"[DEBUG] Posterior mean: min={np.min(mean_toys):.2f}, max={np.max(mean_toys):.2f}, sum={np.sum(mean_toys * bin_widths):.2f}")
# #             print(f"[DEBUG] Number of toys used: {len(toys)}")

# #             # --------------------------------------------------------
# #             # Plot: step-style histograms + ratio panel
# #             # --------------------------------------------------------
# #             fig, (ax_top, ax_bottom) = plt.subplots(
# #                 2, 1,
# #                 figsize=(6.5, 5.5),
# #                 gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05},
# #                 sharex=True
# #             )

# #             # --- Top pad ---
# #             ax_top.set_title(f"{det} — {var}", fontsize=12)
# #             ax_top.fill_between(
# #                 edges[:-1],
# #                 mean_toys - std_toys,
# #                 mean_toys + std_toys,
# #                 step="post",
# #                 color="red",
# #                 alpha=0.3,
# #                 label="±1σ"
# #             )
# #             ax_top.step(edges[:-1], mean_toys, where="post", color="red", lw=1.0, label="Posterior Mean")
# #             ax_top.step(edges[:-1], asimov_vals, where="post", color="black", lw=1.0, label="Asimov")
# #             ax_top.set_ylabel("Events / bin")
# #             ax_top.legend(frameon=False)
# #             ax_top.grid(alpha=0.3)

# #             # --- Bottom pad: relative uncertainty ---
# #             ratio_mean = np.divide(std_toys, mean_toys, out=np.zeros_like(std_toys), where=mean_toys != 0)
# #             ratio_asimov = np.divide(std_toys, asimov_vals, out=np.zeros_like(std_toys), where=asimov_vals != 0)
            
# #             ax_bottom.step(edges[:-1], ratio_mean, where="post", color="red", lw=1.3, label="σ / Posterior Mean")
# #             ax_bottom.step(edges[:-1], ratio_asimov, where="post", color="blue", lw=1.3, label="σ / Asimov")
# #             ax_bottom.set_ylabel("Relative σ", fontsize=10)
# #             ax_bottom.set_xlabel(var)
# #             ax_bottom.legend(frameon=False, fontsize=8)
# #             ax_bottom.grid(alpha=0.3)
# #             max_ratio = max(np.max(ratio_mean), np.max(ratio_asimov))
# #             ax_bottom.set_ylim(0, max_ratio * 1.4)

# #             pdf_pages.savefig(fig)
# #             plt.close(fig)

# #             print(f"[INFO] Added plot for {det}/{var} ({len(toy_keys)} toys)")

# #     pdf_pages.close()
# #     print(f"[INFO] Saved plots to {output_pdf}")

# # # ------------------------------------------------------------
# # # Entry point
# # # ------------------------------------------------------------
# # if __name__ == "__main__":
# #     main()




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
import uproot
import pathlib
import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import ROOT
ROOT.gROOT.SetBatch(True)

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

    # Check if input file exists
    if not pathlib.Path(root_file).exists():
        print(f"[ERROR] Input file not found: {root_file}")
        sys.exit(1)

    print(f"[INFO] Opening ROOT file: {root_file}")
    f = uproot.open(root_file)

    detectors = ["ND", "Other"]
    variables = ["RecoNeutrinoEnergy", "TrueNeutrinoEnergy"]

    pdf_pages = PdfPages(output_pdf)

    # First pass: determine y-axis limits for each detector
    y_limits = {}
    ratio_limits = {}
    for det in detectors:
        max_y = 0
        max_ratio = 0
        for var in variables:
            # Get Asimov
            asimov_dir = "Asimov"
            asimov_node = safe_get_dir(f, asimov_dir)
            if not asimov_node:
                continue
            
            asimov_keys = [
                k for k in asimov_node.keys()
                if f"{det}_" in k and var in k and "_Asimov" in k
            ]
            if not asimov_keys:
                continue
            
            h_asimov = asimov_node[asimov_keys[0]]
            asimov_vals, edges = get_hist_values(h_asimov)
            bin_widths = np.diff(edges)
            asimov_vals = asimov_vals / bin_widths
            
            # Get toys
            toy_dir = f"{det}/{var}"
            toy_node = safe_get_dir(f, toy_dir)
            if not toy_node:
                continue
            
            all_keys = list(toy_node.keys())
            toy_keys = [k for k in all_keys if f"{var}_posterior_toy_" in k and k.startswith(f"{det}_")]
            
            if not toy_keys:
                continue
            
            toys = []
            for key in toy_keys:
                try:
                    toy_hist = toy_node[key]
                    vals = toy_hist.values()
                    toy_edges = toy_hist.axis().edges()
                    toy_bin_widths = np.diff(toy_edges)
                    if len(vals) == len(asimov_vals):
                        toys.append(vals / toy_bin_widths)
                except Exception:
                    pass
            
            if toys:
                toys = np.vstack(toys)
                mean_toys = np.mean(toys, axis=0)
                std_toys = np.std(toys, axis=0)
                
                # Track maximum y value
                max_y = max(max_y, np.max(asimov_vals), np.max(mean_toys + std_toys))
                
                # Track maximum ratio value
                ratio_mean = np.divide(std_toys, mean_toys, out=np.zeros_like(std_toys), where=mean_toys != 0)
                ratio_asimov = np.divide(std_toys, asimov_vals, out=np.zeros_like(std_toys), where=asimov_vals != 0)
                max_ratio = max(max_ratio, np.max(ratio_mean), np.max(ratio_asimov))
        
        y_limits[det] = max_y * 1.1  # Add 10% padding
        ratio_limits[det] = max_ratio * 1.4  # Add 40% padding

    # Second pass: create plots with consistent y-axis
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
                print(f"[DEBUG] Available keys in {asimov_dir}: {list(asimov_node.keys())[:5]}...")
                continue
            print(f"[DEBUG] Using Asimov histogram: {asimov_keys[0]}")
            h_asimov = asimov_node[asimov_keys[0]]
            asimov_vals, edges = get_hist_values(h_asimov)

            bin_widths = np.diff(edges)
            asimov_vals = asimov_vals / bin_widths
            
            print(f"[DEBUG] Asimov: min={np.min(asimov_vals):.2f}, max={np.max(asimov_vals):.2f}, sum={np.sum(asimov_vals * bin_widths):.2f}")

            # --------------------------------------------------------
            # Posterior toy histograms
            # --------------------------------------------------------
            toy_dir = f"{det}/{var}"
            print(f"[INFO] In directory {toy_dir}")
            toy_node = safe_get_dir(f, toy_dir)
            if not toy_node:
                print(f"[WARN] Missing {toy_dir}")
                continue
            
            # Find posterior toy histograms in this directory
            all_keys = list(toy_node.keys())
            print(f"[DEBUG] Found {len(all_keys)} total keys in {toy_dir}")
            print(f"[DEBUG] First 5 keys: {all_keys[:5]}")
            
            # Search for histograms matching pattern: {det}_*_{var}_posterior_toy_XXX
            # More flexible pattern to handle both ND_OffAxisND and Other_FHC_numu naming
            toy_keys = [k for k in all_keys if f"{var}_posterior_toy_" in k and k.startswith(f"{det}_")]
            print(f"[DEBUG] Found {len(toy_keys)} posterior_toy histograms matching pattern")
            
            # Remove duplicates by name (keep first occurrence)
            seen = set()
            unique_toy_keys = []
            for k in toy_keys:
                if k not in seen:
                    seen.add(k)
                    unique_toy_keys.append(k)
            toy_keys = unique_toy_keys
            print(f"[DEBUG] After removing duplicates: {len(toy_keys)} unique histograms")

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
                    toy_hist = toy_node[key]
                    vals = toy_hist.values()
                    toy_edges = toy_hist.axis().edges()
                    toy_bin_widths = np.diff(toy_edges)
                    
                    # Check if binning matches
                    if len(vals) == len(asimov_vals):
                        if not np.allclose(toy_edges, edges):
                            print(f"[WARN] Bin edges don't match for {key}")
                        toys.append(vals / toy_bin_widths)
                except Exception as e:
                    print(f"[WARN] Problem reading {toy_dir}/{key}: {e}")

            if not toys:
                print(f"[WARN] No valid toy histograms for {toy_dir}")
                continue

            toys = np.vstack(toys)
            mean_toys = np.mean(toys, axis=0)
            std_toys = np.std(toys, axis=0)
            
            print(f"[DEBUG] Posterior mean: min={np.min(mean_toys):.2f}, max={np.max(mean_toys):.2f}, sum={np.sum(mean_toys * bin_widths):.2f}")
            print(f"[DEBUG] Number of toys used: {len(toys)}")

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
            ax_top.step(edges[:-1], mean_toys, where="post", color="red", lw=1.0, label="Posterior Mean")
            ax_top.step(edges[:-1], asimov_vals, where="post", color="black", lw=1.0, label="Asimov")
            ax_top.set_ylabel("Events / bin")
            ax_top.legend(frameon=False)
            ax_top.grid(alpha=0.3)
            
            # Set consistent y-axis limit for this detector
            if det in y_limits:
                ax_top.set_ylim(0, y_limits[det])

            # --- Bottom pad: relative uncertainty ---
            ratio_mean = np.divide(std_toys, mean_toys, out=np.zeros_like(std_toys), where=mean_toys != 0)
            ratio_asimov = np.divide(std_toys, asimov_vals, out=np.zeros_like(std_toys), where=asimov_vals != 0)
            
            ax_bottom.step(edges[:-1], ratio_mean, where="post", color="red", lw=1.3, label="σ / Posterior Mean")
            ax_bottom.step(edges[:-1], ratio_asimov, where="post", color="blue", lw=1.3, label="σ / Asimov")
            ax_bottom.set_ylabel("Relative σ", fontsize=10)
            ax_bottom.set_xlabel(var)
            ax_bottom.legend(frameon=False, fontsize=8)
            ax_bottom.grid(alpha=0.3)
            
            # Set consistent ratio y-axis limit for this detector
            if det in ratio_limits:
                ax_bottom.set_ylim(0, ratio_limits[det])

            pdf_pages.savefig(fig)
            plt.close(fig)

            print(f"[INFO] Added plot for {det}/{var} ({len(toy_keys)} toys)")

    pdf_pages.close()
    print(f"[INFO] Saved plots to {output_pdf}")

# ------------------------------------------------------------
# Entry point
# ------------------------------------------------------------
if __name__ == "__main__":
    main()