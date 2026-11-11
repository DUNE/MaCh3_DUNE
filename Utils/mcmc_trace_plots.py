#!/usr/bin/env python3
"""
mcmc_trace_plots.py
====================

Make thesis-quality trace plots from an MCMC chain stored in a ROOT file.

WHAT IT DOES
------------
- Opens a ROOT file with `uproot` (no PyROOT/ROOT install needed).
- Reads a TTree where each entry is one MCMC sample.
- Optionally splits samples into multiple chains/walkers using an
  integer-valued branch (e.g. "chain", "walker", "fChain", ...).
- Draws one trace-plot panel per parameter (value vs. iteration),
  stacked in a single tidy figure, styled for print/PDF inclusion
  in a LaTeX thesis.
- Saves both a vector PDF (for the actual thesis) and a PNG preview.

QUICK START
-----------
    # See what's inside your file/tree first:
    python mcmc_trace_plots.py --file chain.root --tree Chain --list-branches

    # Make trace plots for specific parameters:
    python mcmc_trace_plots.py --file chain.root --tree Chain \
        --params "m_h" "sigma" "alpha_s" \
        --labels "$m_H$ [GeV]" "$\\sigma$" "$\\alpha_s$" \
        --burnin 2000 --outdir plots

    # Auto-detect (numeric, non-meta) branches and plot all of them:
    python mcmc_trace_plots.py --file chain.root --tree Chain --outdir plots

    # Only plot steps 0 to 1,000,000, in black:
    python mcmc_trace_plots.py --file chain.root --tree Chain \
        --step-range 0 1000000 --black --outdir plots

    # Try it right now with no ROOT file, using synthetic demo data:
    python mcmc_trace_plots.py --demo --outdir plots

REQUIREMENTS
------------
    pip install uproot awkward matplotlib numpy
"""

import argparse
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker


# --------------------------------------------------------------------------
# Style: tuned to look right when embedded at ~0.9\textwidth in a LaTeX thesis
# --------------------------------------------------------------------------
def set_thesis_style(fontsize=9):
    plt.rcParams.update({
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "font.family": "serif",
        "font.serif": ["DejaVu Serif", "Computer Modern Roman"],
        "mathtext.fontset": "cm",
        "axes.formatter.use_mathtext": True,
        "axes.labelsize": fontsize + 1,
        "axes.titlesize": fontsize + 1,
        "xtick.labelsize": fontsize,
        "ytick.labelsize": fontsize,
        "legend.fontsize": fontsize,
        "axes.linewidth": 0.7,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.minor.visible": True,
        "ytick.minor.visible": True,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "grid.linewidth": 0.5,
        "legend.frameon": False,
        "pdf.fonttype": 42,   # embed fonts properly (no Type-3 bitmap fonts)
        "ps.fonttype": 42,
        "axes.unicode_minus": False,
    })


# Colorblind-friendly palette (Okabe-Ito), used per chain
PALETTE = [
    "#0072B2", "#D55E00", "#009E73", "#CC79A7",
    "#E69F00", "#56B4E9", "#F0E442", "#000000",
]

# Branches that are almost never physics parameters -> excluded from
# auto-detection by default (extend this list to match your ntuple).
META_BRANCH_PATTERN = re.compile(
    r"^(step|iter|iteration|index|idx|entry|walker|chain|weight|weights|"
    r"logl|loglike|loglikelihood|log_likelihood|logpost|logposterior|"
    r"posterior|prior|neff|accept|acceptance|fChain|fWalker)$",
    re.IGNORECASE,
)


def list_branches(tree):
    print(f"Branches in tree '{tree.name}' ({tree.num_entries} entries):")
    for name, branch in tree.items():
        print(f"  {name:30s} {branch.typename}")


def load_chain(file_path, tree_name, params, chain_branch, param_prefix=None):
    """Read requested branches (plus optional chain-id branch) from ROOT file."""
    import uproot

    f = uproot.open(file_path)
    tree = f[tree_name]

    if params is None:
        candidates = []
        for name, branch in tree.items():
            if param_prefix is not None:
                # explicit prefix given -> just match on that, nothing else
                if not name.startswith(param_prefix):
                    continue
            else:
                if META_BRANCH_PATTERN.match(name):
                    continue
            if chain_branch is not None and name == chain_branch:
                continue
            # keep only numeric, flat (non-jagged) branches
            if "[" in branch.typename or "vector" in branch.typename.lower():
                continue
            candidates.append(name)
        params = candidates
        if not params:
            raise RuntimeError(
                "Could not auto-detect any parameter branches. "
                "Pass --list-branches to inspect the tree and use --params explicitly."
            )
        print(f"[auto-detect] Using branches: {params}")

    to_read = list(params) + ([chain_branch] if chain_branch else [])
    arrays = tree.arrays(to_read, library="np")

    data = {p: np.asarray(arrays[p], dtype=float) for p in params}
    chain_id = np.asarray(arrays[chain_branch]).astype(int) if chain_branch else None
    return data, chain_id, tree, params


def make_demo_data(n_chains=4, n_steps=5000, seed=0):
    """Synthetic correlated MCMC-like chains, for a no-file dry run."""
    rng = np.random.default_rng(seed)
    params = ["mu", "sigma", "rho"]
    data = {p: np.empty(n_chains * n_steps) for p in params}
    chain_id = np.empty(n_chains * n_steps, dtype=int)

    idx = 0
    for c in range(n_chains):
        mu = np.zeros(n_steps)
        sigma = np.zeros(n_steps)
        rho = np.zeros(n_steps)
        mu[0], sigma[0], rho[0] = rng.normal(0, 2, 3)
        for t in range(1, n_steps):
            mu[t] = 0.995 * mu[t - 1] + rng.normal(0, 0.15) + (2.0 if t < 300 else 0.0) * 0.0
            sigma[t] = abs(0.98 * sigma[t - 1] + rng.normal(0, 0.1)) + 0.5
            rho[t] = np.tanh(0.97 * np.arctanh(np.clip(rho[t - 1], -0.99, 0.99)) + rng.normal(0, 0.05))
        # burn-in drift toward the true stationary region
        burn = min(400, n_steps)
        mu[:burn] += np.linspace(4.0, 0.0, burn) * (1 - c * 0.1)
        sl = slice(idx, idx + n_steps)
        data["mu"][sl], data["sigma"][sl], data["rho"][sl] = mu, sigma, rho
        chain_id[sl] = c
        idx += n_steps

    return data, chain_id, params


def trace_plot(data, chain_id, params, labels=None, burnin=0, outdir="plots",
                basename="mcmc_trace", max_points=None, alpha=0.85, linewidth=0.6,
                step_range=None, black=False):
    os.makedirs(outdir, exist_ok=True)
    labels = labels or params
    if len(labels) != len(params):
        raise ValueError("--labels must have the same length as --params")

    n_params = len(params)
    fig_height = 1.55 * n_params + 0.4
    fig, axes = plt.subplots(
        n_params, 1, figsize=(6.3, fig_height), sharex=True,
        gridspec_kw={"hspace": 0.12},
    )
    if n_params == 1:
        axes = [axes]

    chain_ids = np.unique(chain_id) if chain_id is not None else [0]

    for ax, p, lab in zip(axes, params, labels):
        y_full = data[p]
        n_total = len(y_full)

        if chain_id is None:
            groups = [(0, np.arange(n_total))]
        else:
            groups = [(cid, np.where(chain_id == cid)[0]) for cid in chain_ids]

        for i, (cid, sel) in enumerate(groups):
            y = y_full[sel]
            x = np.arange(len(y))

            # Restrict to the requested step range (in steps *within each chain*),
            # applied before any thinning so --max-points still behaves sensibly.
            if step_range is not None:
                start, end = step_range
                end = min(end, len(x))
                if start >= len(x):
                    continue
                x, y = x[start:end], y[start:end]

            if max_points is not None and len(x) > max_points:
                thin = len(x) // max_points
                x, y = x[::thin], y[::thin]

            color = "black" if black else PALETTE[i % len(PALETTE)]
            ax.plot(x, y, color=color, lw=linewidth, alpha=alpha,
                    label=f"chain {cid}" if chain_id is not None else None)

        if burnin > 0:
            ax.axvspan(0, burnin, color="0.85", zorder=0, lw=0)

        ax.set_ylabel("", rotation=0, ha="right", va="center", labelpad=14)
        ax.yaxis.set_major_locator(mticker.MaxNLocator(nbins=4))

    axes[-1].set_xlabel("Step")
    if chain_id is not None and len(chain_ids) > 1 and not black:
        axes[0].legend(loc="upper right", ncol=min(len(chain_ids), 4),
                        bbox_to_anchor=(1.0, 1.35), handlelength=1.2, columnspacing=1.0)

    fig.align_ylabels(axes)

    pdf_path = os.path.join(outdir, f"{basename}.pdf")
    png_path = os.path.join(outdir, f"{basename}.png")
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, bbox_inches="tight")
    plt.close(fig)
    return pdf_path, png_path


def main():
    ap = argparse.ArgumentParser(
        description="Make thesis-quality MCMC trace plots from a ROOT file.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument("--file", help="Path to the ROOT file.")
    ap.add_argument("--tree", help="Name of the TTree inside the ROOT file.")
    ap.add_argument("--params", nargs="+", default=None,
                     help="Branch names to plot. If omitted, auto-detect numeric branches.")
    ap.add_argument("--param-prefix", default=None,
                     help="When --params is omitted, only auto-detect branches starting with "
                          "this prefix (e.g. 'param_'). Overrides the generic meta-branch filter.")
    ap.add_argument("--labels", nargs="+", default=None,
                     help="Axis labels for each parameter (LaTeX ok), same order/length as --params. "
                          "If omitted, the raw branch names are used (prefix included).")
    ap.add_argument("--chain-branch", default=None,
                     help="Integer branch identifying which chain/walker each row belongs to "
                          "(omit if your tree has a single chain).")
    ap.add_argument("--burnin", type=int, default=0,
                     help="Number of initial iterations to shade as burn-in.")
    ap.add_argument("--step-range", type=int, nargs=2, default=None,
                     metavar=("START", "END"),
                     help="Only plot steps in [START, END) within each chain, "
                          "e.g. --step-range 0 1000000.")
    ap.add_argument("--black", action="store_true",
                     help="Draw all traces in plain black instead of the colored palette "
                          "(legend is suppressed since chains become indistinguishable).")
    ap.add_argument("--max-points", type=int, default=20000,
                     help="Thin each trace to at most this many points (keeps PDFs small/fast). "
                          "Use 0 to disable thinning.")
    ap.add_argument("--outdir", default="plots", help="Output directory for figures.")
    ap.add_argument("--basename", default="mcmc_trace", help="Output file base name.")
    ap.add_argument("--list-branches", action="store_true",
                     help="Just print the branches in --tree and exit.")
    ap.add_argument("--demo", action="store_true",
                     help="Ignore --file/--tree and run on synthetic demo data.")
    args = ap.parse_args()

    set_thesis_style()

    if args.demo:
        data, chain_id, params = make_demo_data()
        labels = [r"$\mu$", r"$\sigma$", r"$\rho$"]
        pdf_path, png_path = trace_plot(
            data, chain_id, params, labels=labels, burnin=400,
            outdir=args.outdir, basename=args.basename,
            max_points=(args.max_points or None),
            step_range=args.step_range, black=args.black,
        )
        print(f"Wrote {pdf_path}\nWrote {png_path}")
        return

    if not args.file or not args.tree:
        ap.error("--file and --tree are required (or pass --demo to try synthetic data).")

    import uproot
    f = uproot.open(args.file)
    tree = f[args.tree]

    if args.list_branches:
        list_branches(tree)
        return

    data, chain_id, tree, params = load_chain(
        args.file, args.tree, args.params, args.chain_branch, args.param_prefix
    )
    labels = args.labels
    if labels is None and args.param_prefix:
        # Nicer default labels: strip the shared prefix, e.g. "param_mh" -> "mh"
        labels = [p[len(args.param_prefix):] for p in params]
    pdf_path, png_path = trace_plot(
        data, chain_id, params, labels=labels, burnin=args.burnin,
        outdir=args.outdir, basename=args.basename,
        max_points=(args.max_points or None),
        step_range=args.step_range, black=args.black,
    )
    print(f"Wrote {pdf_path}\nWrote {png_path}")


if __name__ == "__main__":
    sys.exit(main())