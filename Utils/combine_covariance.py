"""
combine_covariance.py
----------------------
Combines covariance matrices (FRACTIONAL or ABSOLUTE) from multiple ROOT
files into a single pooled covariance matrix, equivalent to having run all
throws as one job.

Statistical basis: for independent throws all drawn with the same N per
job, a simple (unweighted) average of the per-job covariance matrices is
the exact pooled estimator. This holds whether the matrix is expressed in
absolute bin-content units (Cov[x_i, x_j]) or fractional/relative units
(Cov[(x_i/nom_i - 1), (x_j/nom_j - 1)]) -- the only thing that changes is
which sanity check makes sense on the diagonal:

  - fractional: sqrt(diag) is a fractional uncertainty per bin, so it
    should typically be << 1 (i.e. << 100%). A huge value usually means
    an absolute matrix was fed in by mistake.
  - absolute: sqrt(diag) is in the same units as the bin content itself,
    so there's no universal "too big" threshold -- instead this script
    just reports the scale so you can eyeball consistency rather than
    apply a fixed cutoff.

Assumes:
  - All ROOT files contain a TMatrixT<double> (or TH2) holding the
    covariance matrix under the same key
  - All throws had the same number of throws (N) per job -- if not, pass
    --weights to do a proper N-weighted average instead of a plain mean

Throw-count consistency check:
  Each input file is expected to contain a per-bin throw-distribution
  histogram (default key "bin_0", e.g. MaCh3's "Bin 0 throw distribution").
  This script reads that histogram's entry count (TH1's fEntries) from
  every file as a proxy for "how many throws went into this file" and
  warns if they don't all match -- catching cases where a stray file with
  a different throw count silently entered the average. Disable with
  --skip-throw-check, or point it at a different histogram with
  --throw-count-key.

Usage:
  # Single matrix
  python combine_covariance.py --pattern "covmat_*.root" --key NDCovMatrix_relative \
      --type fractional --output combined_fraccov.root

  # Both fractional and absolute in one pass -> one output file, two keys
  python combine_covariance.py --pattern "covmat_*.root" --both \
      --output combined_covmat.root
"""

import argparse
import glob
import sys

import numpy as np

try:
    import uproot
except ImportError:
    sys.exit("uproot is required: pip install uproot")


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------

def read_matrix_from_root(filepath: str, key: str) -> np.ndarray:
    """
    Read a covariance matrix from a ROOT file.
    Supports TMatrixT<double> / TMatrixD / TMatrixDSym and TH2 histograms.
    """
    with uproot.open(filepath) as f:
        bare_key = key.split(";")[0]
        obj = None
        for k in [key, bare_key, bare_key + ";1"]:
            if k in f:
                obj = f[k]
                break
        if obj is None:
            available = list(f.keys())
            raise KeyError(
                f"Key '{key}' not found in {filepath}.\n"
                f"Available keys: {available}"
            )

        class_name = obj.classname if hasattr(obj, "classname") else type(obj).__name__

        if "TH2" in class_name:
            return obj.values()

        if "TMatrix" in class_name:
            elements = obj.member("fElements")
            nrow = obj.member("fNrows")
            ncol = obj.member("fNcols")
            return np.array(elements).reshape(nrow, ncol)

        try:
            return np.array(obj.to_numpy())
        except Exception:
            pass

        raise TypeError(
            f"Don't know how to read class '{class_name}' from {filepath}.\n"
            f"Supported types: TMatrixT<double>, TMatrixD, TMatrixDSym, TH2*."
        )


def read_throw_count(filepath: str, hist_key: str) -> "int | None":
    """
    Read a proxy for "number of throws in this file" from a TH1 histogram
    that gets one fill per throw (default: "bin_0", MaCh3's per-bin throw
    distribution). Returns None if the histogram or its entry count isn't
    available, so the caller can decide whether to warn or skip.
    """
    try:
        with uproot.open(filepath) as f:
            bare_key = hist_key.split(";")[0]
            obj = None
            for k in [hist_key, bare_key, bare_key + ";1"]:
                if k in f:
                    obj = f[k]
                    break
            if obj is None:
                return None
            entries = obj.member("fEntries")
            return int(round(entries))
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Combining
# ---------------------------------------------------------------------------

def combine_equal_n(matrices: list) -> np.ndarray:
    """Simple average -- exact pooled estimator when all throws have the same N."""
    return np.stack(matrices, axis=0).mean(axis=0)


def combine_weighted(matrices: list, weights: list) -> np.ndarray:
    """N-weighted average -- use when jobs had different throw counts."""
    w = np.array(weights, dtype=float)
    w = w / w.sum()
    stacked = np.stack(matrices, axis=0)
    return np.tensordot(w, stacked, axes=(0, 0))


def symmetrise(mat: np.ndarray) -> np.ndarray:
    return 0.5 * (mat + mat.T)


def check_positive_semidefinite(mat: np.ndarray) -> None:
    eigvals = np.linalg.eigvalsh(mat)
    min_eig = eigvals.min()
    if min_eig < -1e-10:
        print(
            f"  WARNING: negative eigenvalue detected ({min_eig:.3e}). "
            "This may indicate numerical issues in the input matrices.",
            file=sys.stderr,
        )
    else:
        print(f"  Smallest eigenvalue: {min_eig:.3e}  (positive semi-definite \u2713)")


def check_fractional_scale(mat: np.ndarray, warn_threshold: float = 4.0) -> None:
    """
    Sanity check for FRACTIONAL covariance matrices: sqrt(diag) is a
    fractional (not percent) uncertainty per bin, so a value of 1.0 already
    means 100%. sqrt(diag) >> 1 almost certainly means an absolute-scale
    matrix was loaded by mistake.
    """
    diag = np.diag(mat)
    frac_unc = np.sqrt(np.clip(diag, 0, None))

    print(f"  Fractional uncertainty per bin (sqrt of diagonal): "
          f"min={frac_unc.min():.3g}, median={np.median(frac_unc):.3g}, "
          f"max={frac_unc.max():.3g}")

    n_over = int((frac_unc > warn_threshold).sum())
    if n_over:
        print(
            f"  WARNING: {n_over} bin(s) have fractional uncertainty > "
            f"{warn_threshold:g} (i.e. > {warn_threshold*100:.0f}%). This is "
            "unusually large for a fractional covariance matrix -- double check "
            "the input files actually hold FRACTIONAL (not absolute) "
            "covariances before trusting this output.",
            file=sys.stderr,
        )


def check_absolute_scale(mat: np.ndarray) -> None:
    """
    Sanity check for ABSOLUTE covariance matrices: sqrt(diag) is in the same
    units as the bin content, so there's no fixed pass/fail threshold --
    just report it so it can be compared by eye against the CV spectrum's
    bin contents.
    """
    diag = np.diag(mat)
    abs_unc = np.sqrt(np.clip(diag, 0, None))

    print(f"  Absolute uncertainty per bin (sqrt of diagonal, same units as "
          f"bin content): min={abs_unc.min():.3g}, median={np.median(abs_unc):.3g}, "
          f"max={abs_unc.max():.3g}")

    if np.allclose(abs_unc, 0):
        print(
            "  WARNING: all diagonal entries are ~0. This looks like an empty "
            "or degenerate matrix -- check the input files.",
            file=sys.stderr,
        )
    print(
        "  (No fixed threshold applies to absolute-scale uncertainties -- "
        "sanity check these against your CV spectrum's bin contents by eye.)"
    )


def check_throw_counts(filepaths: list, hist_key: str) -> None:
    """
    Reads a throw-count proxy from each file (entry count of `hist_key`,
    a per-bin throw-distribution histogram) and warns if they differ,
    since combine_equal_n's plain average is only exact when every file
    used the same N.
    """
    counts = {}
    missing = []
    for fp in filepaths:
        n = read_throw_count(fp, hist_key)
        if n is None:
            missing.append(fp)
        else:
            counts[fp] = n

    if missing:
        print(
            f"  Could not read throw count (key '{hist_key}') from {len(missing)} "
            f"file(s); skipping the check for those. Missing from: "
            f"{missing[:3]}{' ...' if len(missing) > 3 else ''}",
            file=sys.stderr,
        )

    if not counts:
        print("  No throw counts could be read -- skipping consistency check.",
              file=sys.stderr)
        return

    unique_counts = set(counts.values())
    if len(unique_counts) == 1:
        n = next(iter(unique_counts))
        print(f"  Throw count consistent across {len(counts)} file(s): {n} throws each \u2713")
    else:
        print(
            f"  WARNING: inconsistent throw counts across input files "
            f"(from key '{hist_key}'): {sorted(unique_counts)}. "
            "combine_equal_n's plain average is only the correct pooled "
            "estimator when every file has the SAME N. Either remove the "
            "mismatched file(s) or use --weights with each file's actual "
            "throw count for a proper N-weighted average.",
            file=sys.stderr,
        )
        for fp, n in counts.items():
            print(f"    {n:>8d}  {fp}", file=sys.stderr)


# ---------------------------------------------------------------------------
# Writing
# ---------------------------------------------------------------------------

def save_multi_to_root(mats: dict, output_path: str) -> None:
    """
    Write multiple named matrices into a single ROOT file, one object per
    key. Tries TMatrixT-compatible output first via PyROOT (all objects in
    one TFile), falls back to TH2D via boost-histogram/uproot (also all in
    one file), then numpy per-matrix as a last resort.
    """
    try:
        import ROOT
        f_out = ROOT.TFile(output_path, "RECREATE")
        for key, mat in mats.items():
            n = mat.shape[0]
            tmat = ROOT.TMatrixD(n, n)
            for i in range(n):
                for j in range(n):
                    tmat[i][j] = mat[i, j]
            tmat.Write(key)
        f_out.Close()
        print(f"  Saved {len(mats)} matrice(s) as TMatrixT<double> in {output_path}  (via PyROOT \u2713)")
        for key in mats:
            print(f"    key: '{key}'")
        return
    except ImportError:
        pass
    except Exception as e:
        print(f"  PyROOT write failed ({e}), trying boost-histogram fallback...")

    try:
        import boost_histogram as bh
        with uproot.recreate(output_path) as f:
            for key, mat in mats.items():
                n = mat.shape[0]
                h = bh.Histogram(
                    bh.axis.Regular(n, 0, n),
                    bh.axis.Regular(n, 0, n),
                    storage=bh.storage.Double(),
                )
                h.view()[:] = mat
                f[key] = h
        print(f"  Saved {len(mats)} matrice(s) as TH2D in {output_path}")
        for key in mats:
            print(f"    key: '{key}'")
        print("  (TH2D stores the same NxN values; readable by ROOT and uproot)")
        return
    except Exception as e:
        print(f"  TH2D write failed ({e}), falling back to numpy...")

    for key, mat in mats.items():
        npy_path = output_path.replace(".root", f"_{key}.npy")
        np.save(npy_path, mat)
        print(f"  Saved '{key}' as numpy binary to {npy_path}")


def save_to_root(mat: np.ndarray, output_path: str, key: str) -> None:
    """Single-matrix convenience wrapper around save_multi_to_root."""
    save_multi_to_root({key: mat}, output_path)


def save_to_npy(mat: np.ndarray, output_path: str) -> None:
    np.save(output_path, mat)
    print(f"  Saved to {output_path}")


def save_to_csv(mat: np.ndarray, output_path: str) -> None:
    np.savetxt(output_path, mat, delimiter=",")
    print(f"  Saved to {output_path}")


# ---------------------------------------------------------------------------
# Core combination routine (shared by single-key and --both modes)
# ---------------------------------------------------------------------------

def combine_one_key(filepaths: list, key: str, matrix_type: str,
                     weights, scale_warn_threshold: float) -> np.ndarray:
    matrices = []
    for fp in filepaths:
        print(f"  Reading {fp} [{key}] ...", end=" ", flush=True)
        mat = read_matrix_from_root(fp, key)
        print(f"shape {mat.shape}")
        matrices.append(mat)

    shapes = [m.shape for m in matrices]
    if len(set(shapes)) != 1:
        sys.exit(f"Inconsistent matrix shapes for key '{key}': {set(shapes)}")
    if matrices[0].ndim != 2 or matrices[0].shape[0] != matrices[0].shape[1]:
        sys.exit(f"Expected square matrices for key '{key}', got {matrices[0].shape}")

    print(f"\nCombining {len(matrices)} {matrix_type} covariance matrices "
          f"(key '{key}') of shape {matrices[0].shape} ...")

    if weights is not None:
        print(f"  Using N-weighted average with weights: {weights}")
        combined = symmetrise(combine_weighted(matrices, weights))
    else:
        print("  Using plain (equal-N) average.")
        combined = symmetrise(combine_equal_n(matrices))

    print("Checking result:")
    check_positive_semidefinite(combined)
    if matrix_type == "fractional":
        check_fractional_scale(combined, warn_threshold=scale_warn_threshold)
    else:
        check_absolute_scale(combined)
    print()

    return combined


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Combine covariance matrices (fractional and/or absolute) from multiple ROOT files."
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--files", nargs="+", metavar="FILE",
                       help="Explicit list of ROOT files.")
    group.add_argument("--pattern", metavar="GLOB",
                       help='Glob pattern, e.g. "covmat_*.root".')

    parser.add_argument("--both", action="store_true",
                        help="Combine BOTH the fractional and absolute matrices in one pass "
                             "and write them into a single output ROOT file under two keys. "
                             "Use --frac-key/--abs-key to set the key names (defaults match "
                             "MaCh3 ND covariance output: NDCovMatrix_relative / NDCovMatrix).")
    parser.add_argument("--frac-key", default="NDCovMatrix_relative",
                        help="[--both mode] Key of the fractional matrix. Default: NDCovMatrix_relative")
    parser.add_argument("--abs-key", default="NDCovMatrix",
                        help="[--both mode] Key of the absolute matrix. Default: NDCovMatrix")

    parser.add_argument("--key",
                        help="[single-matrix mode] ROOT key of the covariance matrix object, "
                             "e.g. NDCovMatrix_relative (fractional) or NDCovMatrix (absolute). "
                             "Required unless --both is set.")
    parser.add_argument("--type", choices=["fractional", "absolute"], default="fractional",
                        help="[single-matrix mode] Which sanity check to apply. Default: fractional.")

    parser.add_argument("--output", default="combined_covmat.root",
                        help="Output file (.root / .npy / .csv). Default: combined_covmat.root. "
                             "Note: .npy/.csv only support single-matrix mode.")
    parser.add_argument("--scale-warn-threshold", type=float, default=4.0,
                        help="[fractional only] Warn if any bin's fractional uncertainty "
                             "(sqrt of diagonal) exceeds this value. Default: 4.0 (400%%).")
    parser.add_argument("--weights", nargs="+", type=float, default=None,
                        help="Per-file throw counts (or any relative weights) if files do "
                             "NOT all have equal N. Must match the number of files, in the "
                             "same order. If omitted, a plain (equal-N) average is used.")

    parser.add_argument("--throw-count-key", default="bin_0",
                        help="Histogram key used to cross-check that every input file has "
                             "the same number of throws (via its entry count). Default: bin_0.")
    parser.add_argument("--skip-throw-check", action="store_true",
                        help="Skip the throw-count consistency check across input files.")

    args = parser.parse_args()

    if not args.both and not args.key:
        parser.error("--key is required unless --both is set.")

    if args.files:
        filepaths = args.files
    else:
        filepaths = sorted(glob.glob(args.pattern))
        if not filepaths:
            sys.exit(f"No files matched: {args.pattern}")

    print(f"Found {len(filepaths)} ROOT file(s).")
    for fp in filepaths:
        print(f"  {fp}")

    if args.weights is not None and len(args.weights) != len(filepaths):
        sys.exit(
            f"--weights has {len(args.weights)} entries but {len(filepaths)} files were found. "
            "These must match 1:1, in order."
        )

    if not args.skip_throw_check:
        print(f"\nChecking throw-count consistency (key '{args.throw_count_key}') ...")
        check_throw_counts(filepaths, args.throw_count_key)

    print()

    if args.both:
        frac_combined = combine_one_key(
            filepaths, args.frac_key, "fractional", args.weights, args.scale_warn_threshold
        )
        abs_combined = combine_one_key(
            filepaths, args.abs_key, "absolute", args.weights, args.scale_warn_threshold
        )

        ext = args.output.lower().rsplit(".", 1)[-1]
        if ext != "root":
            print(f"--both requires a .root output (got '.{ext}'); writing to "
                  f"'{args.output.rsplit('.', 1)[0]}.root' instead.")
            args.output = args.output.rsplit(".", 1)[0] + ".root"

        save_multi_to_root(
            {args.frac_key: frac_combined, args.abs_key: abs_combined},
            args.output,
        )
    else:
        combined = combine_one_key(
            filepaths, args.key, args.type, args.weights, args.scale_warn_threshold
        )

        ext = args.output.lower().rsplit(".", 1)[-1]
        if ext == "root":
            save_to_root(combined, args.output, key=args.key)
        elif ext == "csv":
            save_to_csv(combined, args.output)
        elif ext == "npy":
            save_to_npy(combined, args.output)
        else:
            print(f"Unknown extension '.{ext}', saving as numpy.")
            save_to_npy(combined, args.output + ".npy")

    print("\nDone.")


if __name__ == "__main__":
    main()