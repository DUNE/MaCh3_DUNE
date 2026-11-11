#!/usr/bin/env python3
"""
Batch-generate +1sigma/nominal ratio plots for flux systematics, in the
style of Fig. 1.3 (DecayPipeGeoBField +1sigma/Nominal), for thesis use.

File structure assumed (as produced by e.g. MaCh3 eventrate/template dumps):

    <SystName>/                (TDirectoryFile, one per systematic)
        <SampleDir>/            (TDirectoryFile, one per sample, e.g. "OffAxis_numuCC_numode")
            Variation_0         (TH2D, e.g. -3 sigma)
            Variation_1         (TH2D, e.g. -1 sigma)
            Variation_2         (TH2D, nominal)
            Variation_3         (TH2D, +1 sigma)
            Variation_4         (TH2D, e.g. +3 sigma)

By default Variation_2 = nominal and Variation_3 = +1sigma, matching the
5-point (-3,-1,0,+1,+3) sigma dial convention. Override with
--nominal-variation / --sigma-variation if your convention differs.

Usage:
    python plot_flux_syst_ratios.py \
        --infile eventrate_onaxis_template.root \
        --sample-dir OffAxis_numuCC_numode \
        --outdir thesis_plots
"""

import argparse
import os
import re

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import uproot


def strip_cycle(key):
    """Strip uproot's trailing ';cycle' from a key name."""
    return re.sub(r";\d+$", "", key)


def discover_systematics(f, sample_dir):
    """Find top-level directories that contain <sample_dir>/Variation_*."""
    systematics = []
    skipped = []

    for key in f.keys(recursive=False):
        name = strip_cycle(key)
        obj = f[key]
        if not isinstance(obj, uproot.reading.ReadOnlyDirectory):
            continue

        sample_path = f"{name}/{sample_dir}"
        try:
            f[sample_path]
        except Exception:
            skipped.append(name)
            continue

        systematics.append(name)

    if skipped:
        print(f"Skipping {len(skipped)} top-level entries with no "
              f"'{sample_dir}' subdirectory: {skipped}")

    return sorted(systematics)


def th2_to_arrays(th2):
    values, xedges, yedges = th2.to_numpy()
    return values, xedges, yedges


def plot_ratio(name, nominal_h, sigma_h, outdir, x_label, y_label, z_label,
               vmin, vmax, cmap, dpi):
    nom_vals, xedges, yedges = th2_to_arrays(nominal_h)
    sig_vals, xedges2, yedges2 = th2_to_arrays(sigma_h)

    if not (np.allclose(xedges, xedges2) and np.allclose(yedges, yedges2)):
        print(f"SKIPPING {name}: nominal and +1sigma histograms have different binning")
        return

    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(nom_vals != 0, sig_vals / nom_vals, np.nan)
    # ratio has shape (n_x_bins, n_y_bins)

    # Trim away rows/columns that are entirely NaN (i.e. no real data --
    # outside the detector's actual coverage) so the plot auto-scales to
    # just the populated region instead of showing blank space.
    valid_x = ~np.all(np.isnan(ratio), axis=1)
    valid_y = ~np.all(np.isnan(ratio), axis=0)

    if not valid_x.any() or not valid_y.any():
        print(f"SKIPPING {name}: no non-NaN bins found (nominal is all zero)")
        return

    x_lo, x_hi = np.where(valid_x)[0][[0, -1]]
    y_lo, y_hi = np.where(valid_y)[0][[0, -1]]

    ratio = ratio[x_lo:x_hi + 1, y_lo:y_hi + 1]
    xedges = xedges[x_lo:x_hi + 2]
    yedges = yedges[y_lo:y_hi + 2]

    fig, ax = plt.subplots(figsize=(7, 5.5))

    norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=1.0, vmax=vmax)
    mesh = ax.pcolormesh(
        xedges, yedges, ratio.T, cmap=cmap, norm=norm, shading="auto"
    )

    ax.set_title(f"{name} (+1$\\sigma$ / Nominal)")
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    cbar = fig.colorbar(mesh, ax=ax)
    cbar.set_label(z_label)
    fig.tight_layout()

    safe_name = re.sub(r"[^\w\-.]", "_", name)
    outpath_pdf = os.path.join(outdir, f"{safe_name}_p1sigma_ratio.pdf")
    outpath_png = os.path.join(outdir, f"{safe_name}_p1sigma_ratio.png")
    fig.savefig(outpath_pdf)
    fig.savefig(outpath_png, dpi=dpi)
    plt.close(fig)

    print(f"Wrote {outpath_pdf} and {outpath_png}")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--infile", required=True, help="ROOT file with the systematic dirs")
    parser.add_argument("--sample-dir", required=True,
                         help="Name of the sample subdirectory inside each systematic "
                              "dir, e.g. 'OffAxis_numuCC_numode'")
    parser.add_argument("--outdir", default="thesis_plots", help="Output directory")
    parser.add_argument("--nominal-variation", type=int, default=2,
                         help="Variation_N index that is nominal (default: 2)")
    parser.add_argument("--sigma-variation", type=int, default=3,
                         help="Variation_N index that is +1 sigma (default: 3)")
    parser.add_argument("--x-label", default="Neutrino Energy [GeV]")
    parser.add_argument("--y-label", default="Off-axis [m]")
    parser.add_argument("--z-label", default="+1$\\sigma$ / Nominal",
                         help="Label for the colorbar (z-axis)")
    parser.add_argument("--vmin", type=float, default=0.98)
    parser.add_argument("--vmax", type=float, default=1.02)
    parser.add_argument("--cmap", default="viridis",
                         help="Matplotlib colormap (viridis approximates ROOT's kBird "
                              "blue-green-yellow scheme; try RdYlBu_r for red-blue diverging)")
    parser.add_argument("--dpi", type=int, default=200)
    parser.add_argument("--only", nargs="*", default=None,
                         help="Optional list of systematic names to restrict to")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    nominal_key = f"Variation_{args.nominal_variation}"
    sigma_key = f"Variation_{args.sigma_variation}"

    with uproot.open(args.infile) as f:
        systematics = discover_systematics(f, args.sample_dir)

        if args.only:
            systematics = [s for s in systematics if s in args.only]

        if not systematics:
            print(f"No systematics found with a '{args.sample_dir}' subdirectory. "
                  "Top-level keys in file:")
            for k in f.keys(recursive=False):
                print(f"  {strip_cycle(k)}")
            return

        print(f"Found {len(systematics)} flux systematics: {systematics}")

        for name in systematics:
            base = f"{name}/{args.sample_dir}"
            try:
                nominal_h = f[f"{base}/{nominal_key}"]
                sigma_h = f[f"{base}/{sigma_key}"]
            except Exception as exc:
                print(f"SKIPPING {name}: could not find {nominal_key}/{sigma_key} "
                      f"under {base} ({exc})")
                continue

            plot_ratio(
                name, nominal_h, sigma_h, args.outdir,
                args.x_label, args.y_label, args.z_label,
                args.vmin, args.vmax, args.cmap, args.dpi,
            )


if __name__ == "__main__":
    main()