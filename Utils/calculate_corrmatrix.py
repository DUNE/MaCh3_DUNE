import ROOT
import numpy as np

# -----------------------------
# Configuration
# -----------------------------
input_file = "/scratch/abipeake/Off_axis_chains_new/ORIGINAL_q0q3_afteradaptive_moresteps/ORIGINAL_q0q3_afteradaptive_moresteps_chain_4_job_0.root"
tree_name = "posteriors"
output_file = "nonfrozen_corrmatrix.root"

# List of all parameters
all_params = [f"xsec_{i}" for i in range(1, 3222)]


# -----------------------------
# Open input ROOT file
# -----------------------------
f = ROOT.TFile(input_file, "READ")
tree = f.Get(tree_name)

n_entries_to_check = min(1000, tree.GetEntries())

# -----------------------------
# Check frozen parameters efficiently
# -----------------------------
non_frozen_params = set()
for i, event in enumerate(tree):
    if i >= n_entries_to_check:
        break
    for param in all_params:
        if getattr(event, param) != 0:
            non_frozen_params.add(param)
            
non_frozen_params = list(non_frozen_params)
print(f"Non-frozen parameters ({len(non_frozen_params)}): {non_frozen_params}")

# -----------------------------
# Extract data for non-frozen parameters
# -----------------------------
data = np.zeros((tree.GetEntries(), len(non_frozen_params)))

for idx, param in enumerate(non_frozen_params):
    for entry_idx, event in enumerate(tree):
        data[entry_idx, idx] = getattr(event, param)

# -----------------------------
# Compute covariance and correlation matrices
# -----------------------------
cov_matrix = np.cov(data, rowvar=False)
corr_matrix = np.corrcoef(data, rowvar=False)

# -----------------------------
# Write matrices to ROOT file
# -----------------------------
root_file = ROOT.TFile(output_file, "RECREATE")
n = len(non_frozen_params)

# Covariance
cov_hist = ROOT.TH2D("Covariance", "Covariance Matrix", n, 0, n, n, 0, n)
for i in range(n):
    for j in range(n):
        cov_hist.SetBinContent(i+1, j+1, cov_matrix[i, j])
root_file.WriteTObject(cov_hist)

# Correlation
corr_hist = ROOT.TH2D("Correlation", "Correlation Matrix", n, 0, n, n, 0, n)
for i in range(n):
    for j in range(n):
        corr_hist.SetBinContent(i+1, j+1, corr_matrix[i, j])
root_file.WriteTObject(corr_hist)

root_file.Close()
print(f"Covariance and correlation matrices saved to {output_file}")
