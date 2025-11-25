import argparse
from pathlib import Path
from typing import List
import fnmatch
import ROOT

import uproot
import numpy as np
import argparse
from pathlib import Path
from typing import List
import fnmatch
import uproot
import numpy as np


import argparse
from pathlib import Path
from typing import List
import fnmatch
import uproot
import numpy as np

import argparse
from pathlib import Path
from typing import List
import fnmatch
import uproot
import numpy as np
import ROOT

def save_cov_corr_as_tmatrix_sym(file_path: Path, branch_names: List[str], output_file: Path, step: int = 1, tree_name="posteriors"):
    if not file_path.exists():
        raise FileNotFoundError(f"File {file_path} does not exist.")

    tree: uproot.TTree = uproot.open(f"{file_path}:{tree_name}")

    # Match branches (supports wildcards)
    to_read = []
    for b in branch_names:
        matches = fnmatch.filter(tree.keys(), b)
        if not matches:
            continue
        to_read.extend(matches)

    if not to_read:
        raise ValueError("No branches matched the provided names.")

    # Read in chunks and subsample every 'step' entries
    chunks = []
    for batch in tree.iterate(to_read, library="np", step_size=10000):
        batch_subsample = {k: v[::step] for k, v in batch.items()}
        chunks.append(batch_subsample)

    # Concatenate chunks
    data = {k: np.concatenate([chunk[k] for chunk in chunks]) for k in to_read}

    # Build matrix
    matrix_data = np.vstack([data[b] for b in to_read]).T

    # Remove constant or invalid columns
    #valid_cols = (np.std(matrix_data, axis=0) > 0) & (~np.all(np.isnan(matrix_data), axis=0))
    #matrix_data = matrix_data[:, valid_cols]
    #valid_branches = [b for i, b in enumerate(to_read) if valid_cols[i]]

    # Build matrix
    matrix_data = np.vstack([data[b] for b in to_read]).T

    # Replace NaNs/Infs with 0
    matrix_data = np.nan_to_num(matrix_data, nan=0.0, posinf=0.0, neginf=0.0)

    # Keep all branches
    valid_branches = to_read
    n = len(valid_branches)


    # Replace NaNs/Infs with 0
    matrix_data = np.nan_to_num(matrix_data, nan=0.0, posinf=0.0, neginf=0.0)

    n = len(valid_branches)
    if n == 0:
        raise ValueError("No valid branches left after filtering; cannot build matrices.")

    # Covariance and correlation
    cov_matrix = np.cov(matrix_data, rowvar=False)
    corr_matrix = np.corrcoef(matrix_data, rowvar=False)

    # Convert to TMatrixDSym
    cov_sym = ROOT.TMatrixDSym(n)
    corr_sym = ROOT.TMatrixDSym(n)

    for i in range(n):
        for j in range(n):
            cov_sym[i][j] = cov_matrix[i, j]
            corr_sym[i][j] = corr_matrix[i, j]

    # Save to ROOT file
    f = ROOT.TFile(str(output_file), "RECREATE")
    cov_sym.Write("covariance_matrix")
    corr_sym.Write("correlation_matrix")

    # Save branch names
    branch_array = ROOT.TObjArray()
    for b in valid_branches:
        branch_array.Add(ROOT.TObjString(b))
    branch_array.Write("branch_names")
    f.Close()

    print(f"Saved covariance and correlation matrices as TMatrixDSym to {output_file}")
    print(f"Number of branches: {n}")

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
