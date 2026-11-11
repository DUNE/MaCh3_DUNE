#!/usr/bin/env python3
"""
Merge MaCh3 covariance matrix ROOT files using the Chan et al. parallel
covariance formula. This is statistically exact — it is equivalent to
having run all throws in a single job.

Usage:
    python merge_covmat.py --input-dir covmat_0m_2D/ --output merged_covmat.root
    python merge_covmat.py --input-dir covmat_0m_2D/ --output merged_covmat.root --pattern "covmat_*.root"
"""

import argparse
import glob
import os
import sys

import numpy as np

try:
    import ROOT
    ROOT.gROOT.SetBatch(True)
except ImportError:
    sys.exit("ERROR: PyROOT not available. Make sure your ROOT environment is sourced.")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def tmatrix_to_numpy(tmat):
    """Convert a ROOT TMatrixT<double> to a numpy array."""
    nrows = tmat.GetNrows()
    ncols = tmat.GetNcols()
    arr = np.zeros((nrows, ncols))
    for i in range(nrows):
        for j in range(ncols):
            arr[i, j] = tmat[i][j]
    return arr


def numpy_to_tmatrix(arr):
    """Convert a 2D numpy array to a ROOT TMatrixD."""
    nrows, ncols = arr.shape
    tmat = ROOT.TMatrixD(nrows, ncols)
    for i in range(nrows):
        for j in range(ncols):
            tmat[i][j] = arr[i, j]
    return tmat


def numpy_to_tmatrix_1d(arr):
    """Convert a 1D numpy array to a ROOT TMatrixD column vector (n x 1)."""
    n = arr.shape[0]
    tmat = ROOT.TMatrixD(n, 1)
    for i in range(n):
        tmat[i][0] = arr[i]
    return tmat


def read_file(path):
    """
    Read one covmat ROOT file and return a dict with:
        n_throws  : int   — number of throws (inferred from bin_0 histogram)
        mean      : np.array (n,)
        cov       : np.array (n, n)
        cv        : np.array (n,) or None
        bin_hists : list of TH1D clones (one per bin)
    """
    f = ROOT.TFile.Open(path, "READ")
    if not f or f.IsZombie():
        raise RuntimeError(f"Cannot open {path}")

    cov  = tmatrix_to_numpy(f.Get("NDCovMatrix"))
    mean = tmatrix_to_numpy(f.Get("NDThrowsMean")).flatten()
    cv_obj = f.Get("NDCV")
    cv = tmatrix_to_numpy(cv_obj).flatten() if cv_obj else None

    # Infer n_throws from bin_0 histogram entries
    bin0 = f.Get("bin_0")
    n_throws = int(bin0.GetEntries()) if bin0 else None

    # Clone all bin histograms so they outlive the TFile
    bin_hists = []
    i = 0
    while True:
        h = f.Get(f"bin_{i}")
        if not h:
            break
        clone = h.Clone(f"bin_{i}")
        clone.SetDirectory(0)
        bin_hists.append(clone)
        i += 1

    f.Close()
    return dict(n_throws=n_throws, mean=mean, cov=cov, cv=cv, bin_hists=bin_hists)


# ---------------------------------------------------------------------------
# Chan et al. parallel covariance combination
# ---------------------------------------------------------------------------

def combine_two(nA, meanA, covA, nB, meanB, covB):
    """
    Combine two batches A and B into a single (mean, cov) using the
    Chan / parallel formula. This is algebraically equivalent to computing
    the sample covariance over all nA + nB throws directly.

        delta  = meanB - meanA
        nAB    = nA + nB
        meanAB = meanA + delta * nB / nAB
        covAB  = ( (nA-1)*covA + (nB-1)*covB
                   + delta[:,None]*delta[None,:] * nA*nB/nAB ) / (nAB - 1)
    """
    nAB   = nA + nB
    delta = meanB - meanA
    meanAB = meanA + delta * (nB / nAB)
    outer  = np.outer(delta, delta) * (nA * nB / nAB)
    covAB  = ((nA - 1) * covA + (nB - 1) * covB + outer) / (nAB - 1)
    return nAB, meanAB, covAB


def merge_all(records):
    """Reduce a list of (n, mean, cov) records into one combined (n, mean, cov)."""
    n_acc, mean_acc, cov_acc = records[0]
    for n, mean, cov in records[1:]:
        n_acc, mean_acc, cov_acc = combine_two(n_acc, mean_acc, cov_acc,
                                               n, mean, cov)
    return n_acc, mean_acc, cov_acc


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Merge MaCh3 covmat ROOT files.")
    parser.add_argument("--input-dir",  required=True,  help="Directory containing covmat_*.root files")
    parser.add_argument("--output",     required=True,  help="Output ROOT file path")
    parser.add_argument("--pattern",    default="covmat_*.root",
                        help="Glob pattern for input files (default: covmat_*.root)")
    parser.add_argument("--n-throws",   type=int, default=50000,
                        help="Throws per job if histogram entry count is unavailable (default: 50000)")
    args = parser.parse_args()

    # ---- find input files --------------------------------------------------
    pattern = os.path.join(args.input_dir, args.pattern)
    paths   = sorted(glob.glob(pattern))
    if not paths:
        sys.exit(f"ERROR: No files matched '{pattern}'")
    print(f"Found {len(paths)} input files.")

    # ---- read all files ----------------------------------------------------
    records   = []   # (n_throws, mean, cov)
    cv_master = None
    bin_hists_master = None

    for p in paths:
        print(f"  Reading {p} ...", end=" ", flush=True)
        d = read_file(p)

        n = d["n_throws"] if d["n_throws"] is not None else args.n_throws
        records.append((n, d["mean"], d["cov"]))
        print(f"n_throws={n}  matrix shape={d['cov'].shape}")

        # NDCV should be identical across files — just keep the first one
        if cv_master is None and d["cv"] is not None:
            cv_master = d["cv"]

        # Accumulate bin histograms by adding
        if bin_hists_master is None:
            bin_hists_master = d["bin_hists"]
        else:
            for i, h in enumerate(d["bin_hists"]):
                bin_hists_master[i].Add(h)

    # ---- combine -----------------------------------------------------------
    print("\nCombining with Chan et al. parallel formula ...")
    n_total, mean_combined, cov_combined = merge_all(records)
    print(f"Total throws combined: {n_total}")

    # Sanity checks
    if np.any(np.diag(cov_combined) < 0):
        print("WARNING: Negative diagonal entries found — check inputs!")
    if not np.allclose(cov_combined, cov_combined.T, atol=1e-10):
        print("WARNING: Combined matrix is not symmetric — forcing symmetry.")
        cov_combined = 0.5 * (cov_combined + cov_combined.T)

    # ---- write output ------------------------------------------------------
    print(f"\nWriting output to {args.output} ...")
    out = ROOT.TFile.Open(args.output, "RECREATE")

    numpy_to_tmatrix(cov_combined).Write("NDCovMatrix")
    numpy_to_tmatrix_1d(mean_combined).Write("NDThrowsMean")

    if cv_master is not None:
        numpy_to_tmatrix_1d(cv_master).Write("NDCV")

    if bin_hists_master:
        for h in bin_hists_master:
            h.Write()

    out.Close()
    print("Done.")

    # ---- quick report ------------------------------------------------------
    print(f"\n{'='*50}")
    print(f"  Files merged      : {len(paths)}")
    print(f"  Total throws      : {n_total}")
    print(f"  Matrix dimension  : {cov_combined.shape[0]} x {cov_combined.shape[1]}")
    diag = np.diag(cov_combined)
    print(f"  Diag range        : [{diag.min():.4e}, {diag.max():.4e}]")
    print(f"  Output            : {args.output}")
    print(f"{'='*50}\n")


if __name__ == "__main__":
    main()