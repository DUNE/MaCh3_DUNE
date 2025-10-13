import argparse
from pathlib import Path
from typing import List
import fnmatch
import uproot
import numpy as np
import ROOT

import os


def save_cov_corr_as_tmatrix_sym(file_path: Path,
                                 branch_patterns: List[str],
                                 output_file: Path,
                                 step: int = 1,
                                 tree_name="posteriors"):
    print("[DEBUG] Starting covariance/correlation generation...")

    if output_file.exists():
        os.remove(output_file)


    if not file_path.exists():
        raise FileNotFoundError(f"[ERROR] Input ROOT file does not exist: {file_path}")
    print(f"[DEBUG] Input ROOT file exists: {file_path}")

    tree = uproot.open(f"{file_path}:{tree_name}")
    print(f"[DEBUG] Opened tree '{tree_name}' with {tree.num_entries} entries")



    # Match branches using wildcards
    to_read = []
    for pattern in branch_patterns:
        matches = [b for b in tree.keys() if fnmatch.fnmatch(b, pattern)]
        if not matches:
            print(f"[WARNING] No matches for branch pattern: {pattern}")
        else:
            to_read.extend(matches)

    if not to_read:
        raise ValueError("[ERROR] No branches matched the provided patterns")

    print(f"[INFO] Matched {len(to_read)} branches. First 10:")
    for b in to_read[:10]:
        print(f"  {b}")
    if len(to_read) > 10:
        print(f"  ... ({len(to_read)-10} more)")

    # Load data in chunks and subsample
    chunks = []
    step_size = 10000
    print(f"[DEBUG] Reading data in chunks of {step_size} entries and taking every {step} step")
    for batch in tree.iterate(to_read, library="np", step_size=step_size):
        batch_sub = {k: v[::step] for k, v in batch.items()}
        chunks.append(batch_sub)

    # Concatenate all chunks
    data = {k: np.concatenate([chunk[k] for chunk in chunks]) for k in to_read}
    matrix_data = np.vstack([data[b] for b in to_read]).T
    print(f"[DEBUG] Raw matrix_data shape: {matrix_data.shape}")

    # Replace NaNs/Infs with zeros
    matrix_data = np.nan_to_num(matrix_data, nan=0.0, posinf=0.0, neginf=0.0)

    n = len(to_read)
    print(f"[DEBUG] Number of branches used: {n}")

    # Covariance and correlation
    cov_matrix = np.cov(matrix_data, rowvar=False)
    print(f"[DEBUG] Covariance matrix shape: {cov_matrix.shape}")

    corr_matrix = np.zeros_like(cov_matrix)
    diag = np.sqrt(np.diag(cov_matrix))
    for i in range(n):
        for j in range(n):
            denom = diag[i] * diag[j]
            corr_matrix[i, j] = cov_matrix[i, j] / denom if denom > 0 else 0.0
    print(f"[DEBUG] Correlation matrix shape: {corr_matrix.shape}")

    # Convert to TMatrixDSym
    cov_sym = ROOT.TMatrixDSym(n)
    corr_sym = ROOT.TMatrixDSym(n)
    for i in range(n):
        for j in range(n):
            cov_sym[i][j] = cov_matrix[i, j]
            corr_sym[i][j] = corr_matrix[i, j]

    # Write to ROOT
    f = ROOT.TFile(str(output_file), "RECREATE")
    cov_sym.Write("covariance_matrix")
    corr_sym.Write("correlation_matrix")

    # Correct way to create a single TObjArray
    branch_array = ROOT.TObjArray()  # no length argument
    for b in to_read:
        branch_array.Add(ROOT.TObjString(b))
    branch_array.Write("paramNames")  # write once

    f.Close()
    print(f"[OK] Saved covariance and correlation matrices to {output_file}")
    print("[INFO] First 10 branch names:")
    for i, b in enumerate(to_read[:10]):
        print(f"  {i:3d}: {b} (σ = {diag[i]:.3e})")
    if len(to_read) > 10:
        print(f"  ... ({len(to_read)-10} more)")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute covariance/correlation matrices and save as TMatrixDSym in ROOT.")
    parser.add_argument("file", type=Path, help="Input ROOT file with tree 'posteriors'")
    parser.add_argument("branches", nargs='+', help="Branch names (supports wildcards)")
    parser.add_argument("-o", "--output", type=Path, default=Path("cov_corr_tmatrix.root"), help="Output ROOT file")
    parser.add_argument("--step", type=int, default=1, help="Subsample every Nth entry")
    args = parser.parse_args()

    save_cov_corr_as_tmatrix_sym(args.file, args.branches, args.output, step=args.step)
