#!/usr/bin/env python3
"""
Find (RecoE, TrueE) bins with no events across all off-axis samples.

A cross-section parameter tied to a bin only matters if at least one sample
puts events in that bin, so the bins to drop are those empty in EVERY sample.

Usage
-----
  python3 find_empty_bins.py EventRates.root
  python3 find_empty_bins.py EventRates.root --match "numuCC" --threshold 1e-3 --plot

Input: a ROOT file holding one 2D histogram per sample (RecoE on x, TrueE on y,
as in OffAxisND.yaml). Every TH2 found (searched recursively) is used, unless
you narrow it with --match (a regex on the histogram path). If your histograms
are TrueE on x instead, pass --true-on-x.

Output (prefix set with --out, default "binmask"):
  binmask_summary.txt   counts, plus which samples populate each kept bin
  binmask_bins.csv      one row per bin: indices, edges, total events, kept flag
  binmask_kept.txt      flat indices of bins to keep (for building the
                        reduced covariance / parameter list)
  binmask.png           (with --plot) map of kept / dropped bins
"""
import argparse
import re
import sys

import numpy as np
import uproot


def find_th2(f, pattern):
    out = {}
    for path, cls in f.classnames().items():
        if cls.startswith("TH2") and (pattern is None or re.search(pattern, path)):
            out[path.split(";")[0]] = f[path]
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("rootfile")
    ap.add_argument("--match", default=None,
                    help="regex: only use histograms whose path matches")
    ap.add_argument("--threshold", type=float, default=0.0,
                    help="bin counts as empty if its summed events are <= this "
                         "(default 0, i.e. strictly empty)")
    ap.add_argument("--true-on-x", action="store_true",
                    help="histograms have TrueE on x and RecoE on y")
    ap.add_argument("--out", default="binmask", help="output file prefix")
    ap.add_argument("--plot", action="store_true")
    a = ap.parse_args()

    f = uproot.open(a.rootfile)
    hists = find_th2(f, a.match)
    if not hists:
        sys.exit("No TH2 histograms found (check --match). "
                 f"File contains: {sorted(set(f.classnames().values()))}")

    names = list(hists)
    first = hists[names[0]]
    xe, ye = first.axis(0).edges(), first.axis(1).edges()
    for n in names[1:]:
        h = hists[n]
        if not (np.allclose(h.axis(0).edges(), xe) and np.allclose(h.axis(1).edges(), ye)):
            sys.exit(f"Binning of {n} differs from {names[0]}; can't combine.")

    # stack: (sample, x, y); flow bins excluded
    rates = np.stack([hists[n].values(flow=False) for n in names])
    if a.true_on_x:
        rates = rates.transpose(0, 2, 1)
        xe, ye = ye, xe
    reco_e, true_e = xe, ye
    total = rates.sum(axis=0)                     # (reco, true)
    nsamp_filled = (rates > a.threshold).sum(axis=0)
    keep = total > a.threshold

    nreco, ntrue = total.shape
    nbins = total.size
    nkeep = int(keep.sum())

    # Flat index: reco is the fast index (ireco + nreco*itrue). Change here if
    # your parameter ordering differs.
    flat = lambda ir, it: ir + nreco * it

    with open(f"{a.out}_bins.csv", "w") as fo:
        fo.write("flat_index,reco_bin,true_bin,reco_lo,reco_hi,true_lo,true_hi,"
                 "total_events,n_samples_with_events,keep\n")
        for it in range(ntrue):
            for ir in range(nreco):
                fo.write(f"{flat(ir,it)},{ir},{it},{reco_e[ir]:g},{reco_e[ir+1]:g},"
                         f"{true_e[it]:g},{true_e[it+1]:g},{total[ir,it]:.6g},"
                         f"{nsamp_filled[ir,it]},{int(keep[ir,it])}\n")

    kept_idx = sorted(flat(ir, it) for ir, it in zip(*np.nonzero(keep)))
    with open(f"{a.out}_kept.txt", "w") as fo:
        fo.write("\n".join(map(str, kept_idx)) + "\n")

    # Bins kept only because of a tiny contribution are worth a look
    nonzero = total[keep]
    frac = nonzero / total.sum()
    lines = [
        f"File: {a.rootfile}",
        f"Samples used ({len(names)}):", *[f"  {n}" for n in names],
        "",
        f"Binning: {nreco} reco x {ntrue} true = {nbins} bins",
        f"Empty threshold: <= {a.threshold:g} events (summed over samples)",
        f"Bins with events (keep): {nkeep}",
        f"Empty bins (drop):       {nbins - nkeep}  ({100*(nbins-nkeep)/nbins:.1f}%)",
        "",
        "Kept bins holding < 0.01% / 0.1% / 1% of all events: "
        f"{(frac < 1e-4).sum()} / {(frac < 1e-3).sum()} / {(frac < 1e-2).sum()}",
        "",
        "Per-sample populated bins:",
        *[f"  {n}: {(rates[i] > a.threshold).sum()}" for i, n in enumerate(names)],
        "",
        "Kept bins per true-energy bin:",
        *[f"  TrueE [{true_e[it]:g}, {true_e[it+1]:g}): {keep[:, it].sum()}/{nreco}"
          for it in range(ntrue)],
    ]
    with open(f"{a.out}_summary.txt", "w") as fo:
        fo.write("\n".join(lines) + "\n")
    print("\n".join(lines[: lines.index("Per-sample populated bins:")]))
    print(f"\nWrote {a.out}_summary.txt, {a.out}_bins.csv, {a.out}_kept.txt")

    if a.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm
        fig, axs = plt.subplots(1, 2, figsize=(13, 5.5), constrained_layout=True)
        masked = np.ma.masked_where(~keep, total)
        m = axs[0].pcolormesh(reco_e, true_e, masked.T,
                              norm=LogNorm(vmin=max(masked.min(), 1e-6), vmax=masked.max()))
        fig.colorbar(m, ax=axs[0], label="events (all samples)")
        axs[0].set_title("Summed event rate (white = empty)")
        m2 = axs[1].pcolormesh(reco_e, true_e, nsamp_filled.T, cmap="viridis")
        fig.colorbar(m2, ax=axs[1], label="# samples with events")
        axs[1].set_title(f"Kept {nkeep} / {nbins} bins")
        for ax in axs:
            ax.set_xlabel("Reco $E_\\nu$ [GeV]")
            ax.set_ylabel("True $E_\\nu$ [GeV]")
        fig.savefig(f"{a.out}.png", dpi=130)
        print(f"Wrote {a.out}.png")


if __name__ == "__main__":
    main()