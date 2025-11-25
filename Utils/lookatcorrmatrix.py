import ROOT
import numpy as np

# Open the ROOT file
f = ROOT.TFile("cov_corr_tmatrix.root")

# Load the correlation matrix
corr_sym = f.Get("correlation_matrix")  # TMatrixDSym

n = corr_sym.GetNrows()

# Convert to NumPy array
corr_matrix = np.array([[corr_sym[i][j] for j in range(n)] for i in range(n)])

print("Shape:", corr_matrix.shape)
print("First 10x10 block:")
print(corr_matrix[:100, :100])  # Print a small block

# Assuming corr_matrix is your NumPy array from TMatrixDSym
corr_matrix = np.nan_to_num(corr_matrix, nan=0.0)

# Optional: check result
print("Number of NaNs after replacement:", np.isnan(corr_matrix).sum())
print(corr_matrix[:100, :100]) 

