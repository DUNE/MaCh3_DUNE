import argparse
from pathlib import Path
from typing import List
import fnmatch
import uproot
import numpy as np
import ROOT

def save_cov_corr_as_tmatrix_sym(file_path: Path,
                                 branch_names: List[str],
                                 output_file: Path,
                                 step: int = 1,
                                 tree_name="posteriors"):
    print("[DEBUG] Script starting...")
    if not file_path.exists():
        raise FileNotFoundError(f"File {file_path} does not exist.")
    print(f"[DEBUG] Input ROOT file exists: {file_path}")

    tree: uproot.TTree = uproot.open(f"{file_path}:{tree_name}")
    
    print(f"[DEBUG] Opened tree '{tree_name}' with {tree.num_entries} entries")

    # Match branches (supports wildcards)
    to_read = []
    for b in branch_names:
        matches = [key for key in tree.keys() if fnmatch.fnmatch(key, b)]
        if not matches:
            print(f"[Warning] No matches for branch pattern: {b}")
        else:
            to_read.extend(matches)

    if not to_read:
        raise ValueError("No branches matched the provided names.")
    print(f"[INFO] Matched {len(to_read)} branches. First 10:")
    for b in to_read[:10]:
        print(f"  {b}")
    if len(to_read) > 10:
        print(f"  ... ({len(to_read)-10} more)")

    # Read in chunks and subsample every 'step' entries
    chunks = []
    for batch in tree.iterate(to_read, library="np", step_size=10000):
        batch_subsample = {k: v[::step] for k, v in batch.items()}
        chunks.append(batch_subsample)
    print(f"[DEBUG] Loaded {len(chunks)} chunks")

    # Concatenate chunks into single array
    data = {k: np.concatenate([chunk[k] for chunk in chunks]) for k in to_read}
    matrix_data = np.vstack([data[b] for b in to_read]).T
    matrix_data = np.nan_to_num(matrix_data, nan=0.0, posinf=0.0, neginf=0.0)
    print(f"[DEBUG] matrix_data shape: {matrix_data.shape}")

    n = len(to_read)
    if matrix_data.shape[0] == 0:
        raise ValueError("No data loaded after subsampling — cannot build matrices.")

    # --- Covariance & correlation ---
    cov_matrix = np.cov(matrix_data, rowvar=False)
    corr_matrix = np.zeros_like(cov_matrix)
    diag = np.sqrt(np.diag(cov_matrix))
    for i in range(n):
        for j in range(n):
            denom = diag[i] * diag[j]
            corr_matrix[i, j] = cov_matrix[i, j] / denom if denom > 0 else 0.0
    print(f"[DEBUG] Covariance and correlation matrices computed")
    print(f"[DEBUG] Covariance matrix shape: {cov_matrix.shape}")
    print(f"[DEBUG] Correlation matrix shape: {corr_matrix.shape}")

    # --- Convert to TMatrixDSym ---
    cov_sym = ROOT.TMatrixDSym(n)
    corr_sym = ROOT.TMatrixDSym(n)
    for i in range(n):
        for j in range(n):
            cov_sym[i][j] = cov_matrix[i, j]
            corr_sym[i][j] = corr_matrix[i, j]

    # --- Save to ROOT file ---
    f = ROOT.TFile(str(output_file), "RECREATE")
    cov_sym.Write("covariance_matrix")
    corr_sym.Write("correlation_matrix")

    # Save branch names as single TObjArray
    branch_array = ROOT.TObjArray(len(to_read))
    branch_array.SetName("branch_names")
    for b in to_read:
        branch_array.Add(ROOT.TObjString(b))
    branch_array.Write()
    f.Close()

    print(f"[OK] Saved covariance and correlation matrices to {output_file}")
    print(f"[OK] Number of branches: {n}")
    print("[INFO] First 10 branch names with σ:")
    for i, b in enumerate(to_read[:10]):
        print(f"  {i:3d}: {b} (σ = {diag[i]:.3e})")
    if len(to_read) > 10:
        print(f"  ... ({len(to_read)-10} more)")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute covariance/correlation matrices and save as TMatrixDSym in a ROOT file.")
    parser.add_argument("file", type=Path, help="Path to the input ROOT file.")
    parser.add_argument("branches", nargs='+', help="List of branch names (supports wildcards).")
    parser.add_argument("-o", "--output", type=Path, default=Path("cov_corr_tmatrix.root"),
                        help="Output ROOT file")
    parser.add_argument("--step", type=int, default=1,
                        help="Take every N-th entry (default=1, i.e., all entries)")
    args = parser.parse_args()

    save_cov_corr_as_tmatrix_sym(args.file, args.branches, args.output, step=args.step)
