# # import ROOT

# # # Open ROOT file and get histogram
# # file = ROOT.TFile("/scratch/abipeake/MaCh3_DUNE_merged/MaCh3_DUNE/enubiasbinning.root", "READ")
# # hist = file.Get("ND_FHC_CCnumu_generic2D")

# # # Total event rate
# # total_events = hist.Integral()
# # print(f"Total event rate: {total_events}\n")

# # fraction_of_event_rate = 0.001 * total_events

# # # Axis titles
# # x_title = hist.GetXaxis().GetTitle()
# # y_title = hist.GetYaxis().GetTitle()
# # print(f"X-axis: {x_title}")
# # print(f"Y-axis: {y_title}\n")

# # # Bin counts
# # nx = hist.GetNbinsX()
# # ny = hist.GetNbinsY()

# # for ix in range(1, nx+1):
# #     x_low = hist.GetXaxis().GetBinLowEdge(ix)
# #     x_high = hist.GetXaxis().GetBinUpEdge(ix)
# #     print(f"{x_title} bin {ix}: [{x_low:.2f}, {x_high:.2f}]")
# #     print(f"{y_title:>6} | {'Bin Edges':>30} | {'Events':>10}")
# #     print("-" * 50)

# #     # Variables to track merging
# #     merged_y_edges = []  # store individual bin edges
# #     merged_events = 0

# #     for iy in range(1, ny+1):
# #         y_low = hist.GetYaxis().GetBinLowEdge(iy)
# #         y_high = hist.GetYaxis().GetBinUpEdge(iy)
# #         events_in_bin = hist.GetBinContent(ix, iy)

# #         if events_in_bin >= fraction_of_event_rate:
# #             # Flush any accumulated merged bin first
# #             if merged_events > 0:
# #                 print(f"  Merged {y_title} bin {merged_y_edges}: {merged_events:.2f} events")
# #                 merged_events = 0
# #                 merged_y_edges = []

# #             # Print this large bin as is
# #             print(f"  {y_title} bin [{y_low:.2f}, {y_high:.2f}]: {events_in_bin:.2f} events")

# #         else:
# #             # Bin is below threshold: merge it
# #             merged_events += events_in_bin
# #             merged_y_edges.append(f"{y_low:.2f}")

# #             # Always include the upper edge of the last bin for clarity
# #             y_high_last = y_high

# #             # Flush merged bin if threshold reached
# #             if merged_events >= fraction_of_event_rate:
# #                 # Add upper edge of last bin to list
# #                 merged_y_edges.append(f"{y_high_last:.2f}")
# #                 print(f"  Merged {y_title} bin [{', '.join(merged_y_edges)}]: {merged_events:.2f} events")
# #                 merged_events = 0
# #                 merged_y_edges = []

# #     # Flush any remaining merged bin at the end of the column
# #     if merged_events > 0:
# #         merged_y_edges.append(f"{y_high_last:.2f}")
# #         print(f"  Merged {y_title} bin [{', '.join(merged_y_edges)}]: {merged_events:.2f} events")

# #     print("\n")  # blank line between x bins
# import ROOT

# # Open ROOT file and get histogram
# file = ROOT.TFile("/scratch/abipeake/MaCh3_DUNE_merged/MaCh3_DUNE/enubiasbinning.root", "READ")
# hist = file.Get("ND_FHC_CCnumu_generic2D")

# # Total event rate
# total_events = hist.Integral()
# print(f"Total event rate: {total_events:.2f}\n")

# # Fraction threshold for merging
# fraction_of_event_rate = 0.001 * total_events
# print(f"Merging threshold per column: {fraction_of_event_rate:.2f} events\n")

# # Axis titles
# x_title = hist.GetXaxis().GetTitle()
# y_title = hist.GetYaxis().GetTitle()

# nx = hist.GetNbinsX()
# ny = hist.GetNbinsY()

# # === Store all results here ===
# all_new_y_edges = []  # vector of vectors: [[edges for x1], [edges for x2], ...]

# for ix in range(1, nx + 1):
#     x_low = hist.GetXaxis().GetBinLowEdge(ix)
#     x_high = hist.GetXaxis().GetBinUpEdge(ix)

#     merged_events = 0
#     new_y_edges = [hist.GetYaxis().GetBinLowEdge(1)]  # start at lowest y edge

#     for iy in range(1, ny + 1):
#         y_low = hist.GetYaxis().GetBinLowEdge(iy)
#         y_high = hist.GetYaxis().GetBinUpEdge(iy)
#         events_in_bin = hist.GetBinContent(ix, iy)

#         if events_in_bin >= fraction_of_event_rate:
#             # Flush any active merged bin
#             if merged_events > 0:
#                 new_y_edges.append(y_high_prev)
#                 merged_events = 0
#             # Standalone bin
#             new_y_edges.append(y_high)
#         else:
#             # Add small bin to merged total
#             merged_events += events_in_bin
#             y_high_prev = y_high
#             if merged_events >= fraction_of_event_rate:
#                 new_y_edges.append(y_high)
#                 merged_events = 0

#     # Flush any remaining merged bin
#     if merged_events > 0 and new_y_edges[-1] != hist.GetYaxis().GetBinUpEdge(ny):
#         new_y_edges.append(hist.GetYaxis().GetBinUpEdge(ny))

#     # Save this column’s y-edges
#     all_new_y_edges.append(new_y_edges)

#     # Print summary table
#     print("=" * 90)
#     print(f"{x_title} bin {ix}: [{x_low:.2f}, {x_high:.2f}]")
#     print(f"New {y_title} bin edges ({len(new_y_edges) - 1} bins):")
#     print(", ".join([f"{edge:.2f}" for edge in new_y_edges]))
#     print("\n")

# # === Print the vector of vectors at the end ===
# print("\n==================== Summary: Vector of Vectors ====================")
# print("[")
# for i, edges in enumerate(all_new_y_edges, start=1):
#     formatted = ", ".join([f"{edge:.4f}" for edge in edges])
#     print(f"  [{formatted}],  # x-bin {i}")
# print("]")


import ROOT

# -------------------------------
# Configuration
# -------------------------------
ROOT_FILE = "/scratch/abipeake/MaCh3_DUNE_merged/MaCh3_DUNE/enubiasbinning.root"
HIST_NAME = "ND_FHC_CCnumu_generic2D"
THRESHOLD_FRACTION = 0.01  # fraction of total events
OUTPUT_YAML = "TrueNeutrinoEnergy_Enubias_merged_global.yaml"
PARAM_LIST_FILE = "TrueNeutrinoEnergy_Enubias_merged_global_params.txt"

# -------------------------------
# Open ROOT file
# -------------------------------
file = ROOT.TFile(ROOT_FILE, "READ")
if not file.IsOpen():
    raise RuntimeError(f"ROOT file not found: {ROOT_FILE}")

hist = file.Get(HIST_NAME)
if not hist or not hist.InheritsFrom("TH2"):
    raise RuntimeError(f"Histogram '{HIST_NAME}' not found or not TH2.")

nx = hist.GetNbinsX()
ny = hist.GetNbinsY()


# Compute threshold for removing events
# -------------------------------
total_events = sum(hist.GetBinContent(ix, iy) for ix in range(1, nx+1) for iy in range(1, ny+1))
global_threshold = THRESHOLD_FRACTION * total_events
print(f"Total events: {total_events:.1f}, Global threshold: {global_threshold:.1f}")

# -------------------------------
# Prepare outputs
# -------------------------------
all_new_y_edges = []
yaml_lines = ["Systematics:"]
param_lines = ["ParameterName, ene_min, ene_max, enubias_min, enubias_max"]
count = 1

# -------------------------------
# Merge bins along Y-axis for each X-bin
# -------------------------------
for ix in range(1, nx+1):
    x_low = hist.GetXaxis().GetBinLowEdge(ix)
    x_high = hist.GetXaxis().GetBinUpEdge(ix)

    merged_events = 0
    new_y_edges = [hist.GetYaxis().GetBinLowEdge(1)]
    y_low_merge = new_y_edges[0]

    for iy in range(1, ny+1):
        y_high = hist.GetYaxis().GetBinUpEdge(iy)
        events_in_bin = hist.GetBinContent(ix, iy)

        merged_events += events_in_bin

        # Store merged bin if  event threshold passed
        if merged_events >= global_threshold:
            new_y_edges.append(y_high)
            merged_events = 0
            y_low_merge = y_high

    # Add last bin edge if needed
    if new_y_edges[-1] != hist.GetYaxis().GetBinUpEdge(ny):
        new_y_edges.append(hist.GetYaxis().GetBinUpEdge(ny))

    all_new_y_edges.append(new_y_edges)

    # -------------------------------
    # Generate YAML & param list
    # -------------------------------
    for j in range(len(new_y_edges)-1):
        y_low_bin = new_y_edges[j]
        y_high_bin = new_y_edges[j+1]
        param_name = f"enubias_enu_{count}"
        count += 1

        yaml_block = f"""- Systematic:
    SampleNames: ["ND*", "FD*"]
    Error: 0.5
    FlatPrior: false
    KinematicCuts:
    - TrueNeutrinoEnergy:
      - {x_low}
      - {x_high}
    - Enubias:
      - {y_low_bin}
      - {y_high_bin}
    Names:
      FancyName: {param_name}
      ParameterName: {param_name}
    ParameterBounds:
    - 0.0
    - 4.0
    ParameterValues:
      Generated: 1.0
      PreFitValue: 1.0
    StepScale:
      MCMC: 1.0
    Type: Norm
    ParameterGroup: Xsec"""
        yaml_lines.append(yaml_block)
        param_lines.append(f"{param_name}, {x_low}, {x_high}, {y_low_bin}, {y_high_bin}")

# -------------------------------
# Write output files
# -------------------------------
with open(OUTPUT_YAML, "w") as f:
    f.write("\n".join(yaml_lines))

with open(PARAM_LIST_FILE, "w") as f:
    f.write("\n".join(param_lines))

# -------------------------------
# Print summary
# -------------------------------
print(f"YAML config generated: {OUTPUT_YAML}")
print(f"Parameter list generated: {PARAM_LIST_FILE}")
print("\nVector of vectors of new y-edges per x-bin:")
print("[")
for i, edges in enumerate(all_new_y_edges, start=1):
    formatted = ", ".join([f"{e:.4f}" for e in edges])
    print(f"  [{formatted}],  # x-bin {i}")
print("]")
