#!/usr/bin/env python3
"""
plot_sigma_variations.py

Make thesis-quality 1D sigma-variation plots from a ROOT file structured as:

    tmplw_e29_eb10/OffAxis0m_numuCC_numode/Variation_0   (-3 sigma)
    tmplw_e29_eb10/OffAxis0m_numuCC_numode/Variation_1   (-1 sigma)
    tmplw_e29_eb10/OffAxis0m_numuCC_numode/Variation_2   (nominal)
    tmplw_e29_eb10/OffAxis0m_numuCC_numode/Variation_3   (+1 sigma)
    tmplw_e29_eb10/OffAxis0m_numuCC_numode/Variation_4   (+3 sigma)

Each top-level TDirectoryFile is one systematic parameter (an energy bin of a
binned normalisation systematic). The physical meaning of each directory name
(true-neutrino-energy range + Enubias range) is read from the corresponding
MaCh3/generic systematics YAML config, matched via the `FancyName` field.

Requires: uproot, matplotlib, numpy, pyyaml
    pip install uproot matplotlib numpy pyyaml

Usage:
    python plot_sigma_variations.py \
        --root-file templates.root \
        --yaml-file systematics.yaml \
        --subdir OffAxis0m_numuCC_numode \
        --outdir plots \
        [--prefix tmplw_e29_eb] \
        [--ratio] \
        [--bin-width-norm] \
        [--formats pdf png]

By default --prefix is 'tmplw_', so ALL top-level directories are plotted
(every true-energy bin, every bias bin). Narrow --prefix to e.g.
'tmplw_e29_eb' to restrict to a single true-energy bin's sweep.

If --yaml-file is omitted, raw directory names are used as labels instead.
"""

import argparse
import os
import re
import sys
import warnings

warnings.filterwarnings("ignore", message=".*not compatible with tight_layout.*")

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import uproot

try:
    import yaml
    HAVE_YAML = True
except ImportError:
    HAVE_YAML = False


# ---------------------------------------------------------------------------
# Sigma-level bookkeeping. Variation_0..4 -> sigma value, label, color, style.
# Edit here if your ordering/convention differs.
# ---------------------------------------------------------------------------
VARIATION_SIGMA = {
    "Variation_0": -3,
    "Variation_1": -1,
    "Variation_2": 0,   # nominal
    "Variation_3": +1,
    "Variation_4": +3,
}

# Consistent, thesis-friendly color scheme: symmetric sigma levels share hue,
# darker = larger |sigma|. Nominal is black.
SIGMA_STYLE = {
    -3: dict(color="#08306b", ls="--", lw=1.6, label=r"$-3\sigma$"),
    -1: dict(color="#4292c6", ls="--", lw=1.6, label=r"$-1\sigma$"),
     0: dict(color="black",   ls="-",  lw=2.2, label="Nominal"),
    +1: dict(color="#fb6a4a", ls="--", lw=1.6, label=r"$+1\sigma$"),
    +3: dict(color="#a50f15", ls="--", lw=1.6, label=r"$+3\sigma$"),
}


def load_yaml_mapping(yaml_path):
    """Parse the systematics YAML into {FancyName: {cut_name: [lo, hi], ...}}."""
    if not HAVE_YAML:
        print("WARNING: pyyaml not installed, cannot read --yaml-file. "
              "Falling back to raw directory names.", file=sys.stderr)
        return {}

    with open(yaml_path) as f:
        data = yaml.safe_load(f)

    mapping = {}
    for entry in data.get("Systematics", []):
        s = entry.get("Systematic", {})
        name = s.get("Names", {}).get("FancyName")
        if name is None:
            continue
        cuts = {}
        for cut in s.get("KinematicCuts", []) or []:
            cuts.update(cut)
        mapping[name] = cuts
    return mapping


def fmt_num(x):
    """Compact, thesis-friendly number formatting."""
    if x == int(x):
        return f"{int(x)}"
    return f"{x:g}"


def make_label(dirname, mapping):
    """
    Build a human-readable label for a directory name using the YAML mapping.
    Falls back to the raw directory name if no mapping entry is found.
    """
    cuts = mapping.get(dirname)
    if not cuts:
        return dirname

    parts = []
    if "TrueNeutrinoEnergy" in cuts:
        lo, hi = cuts["TrueNeutrinoEnergy"]
        parts.append(rf"$E_\nu^{{true}}\in[{fmt_num(lo)},{fmt_num(hi)}]$ GeV")
    if "Enubias" in cuts:
        lo, hi = cuts["Enubias"]
        parts.append(rf"bias$\in[{fmt_num(lo)},{fmt_num(hi)}]$")
    # Any other kinematic cuts present in the YAML get appended generically.
    for key, (lo, hi) in cuts.items():
        if key in ("TrueNeutrinoEnergy", "Enubias"):
            continue
        parts.append(rf"{key}$\in[{fmt_num(lo)},{fmt_num(hi)}]$")

    return "\n".join(parts) if parts else dirname


def th1_to_arrays(hist):
    """uproot TH1 -> (bin_edges, values), stripping under/overflow."""
    values, edges = hist.to_numpy()
    return edges, values


def discover_directories(root_file, prefix):
    """
    Top-level directory keys matching the given prefix, naturally sorted.

    Names of the form '..._e<E>_eb<B>' are sorted numerically by (E, B) so
    the full sweep across true-energy bins and bias bins comes out in a
    sensible order rather than alphabetically (e.g. e9 before e10).
    Falls back to a generic trailing-number sort for any other naming.
    """
    names = set()
    for key in root_file.keys():
        top = key.split(";")[0].split("/")[0]
        if top.startswith(prefix):
            names.add(top)

    e_eb_pattern = re.compile(r"e(\d+)_eb(\d+)$")

    def natural_key(s):
        m = e_eb_pattern.search(s)
        if m:
            return (0, int(m.group(1)), int(m.group(2)), s)
        m = re.search(r"(\d+)$", s)
        if m:
            return (1, 0, int(m.group(1)), s)
        return (2, 0, 0, s)

    return sorted(names, key=natural_key)


def plot_one(dirname, hists, label, outdir, formats, ratio, subdir, bin_width_norm=False):
    """Make a single sigma-variation figure for one directory."""
    fig_h = 6.0 if ratio else 5.0
    fig, axes = plt.subplots(
        2 if ratio else 1, 1, sharex=True, figsize=(7, fig_h),
        gridspec_kw=({"height_ratios": [3, 1], "hspace": 0.06} if ratio else {}),
        squeeze=False,
    )
    ax = axes[0, 0]
    ax_ratio = axes[1, 0] if ratio else None

    nominal_edges, nominal_vals = None, None
    # For the ratio panel we always want the raw (un-normalised) values, since
    # bin width cancels in a variation/nominal ratio and normalising both
    # numerator and denominator is equivalent to not normalising at all.
    nominal_vals_raw = None
    for var_name in ["Variation_0", "Variation_1", "Variation_2", "Variation_3", "Variation_4"]:
        if var_name not in hists:
            continue
        sigma = VARIATION_SIGMA[var_name]
        edges, vals = th1_to_arrays(hists[var_name])
        plot_vals = vals / np.diff(edges) if bin_width_norm else vals
        style = SIGMA_STYLE[sigma]
        ax.step(edges, np.append(plot_vals, plot_vals[-1]), where="post", **style)
        if sigma == 0:
            nominal_edges, nominal_vals = edges, plot_vals
            nominal_vals_raw = vals

    ax.set_ylabel("Events / GeV" if bin_width_norm else "Events / bin")
    ax.legend(frameon=False, fontsize=9, ncol=1, loc="best")
    ax.set_title(label, fontsize=11)
    ax.tick_params(direction="in", top=True, right=True, which="both")
    ax.margins(x=0)

    if ratio and nominal_vals_raw is not None:
        for var_name in ["Variation_0", "Variation_1", "Variation_3", "Variation_4"]:
            if var_name not in hists:
                continue
            sigma = VARIATION_SIGMA[var_name]
            edges, vals = th1_to_arrays(hists[var_name])
            with np.errstate(divide="ignore", invalid="ignore"):
                r = np.where(nominal_vals_raw != 0, vals / nominal_vals_raw, np.nan)
            style = dict(SIGMA_STYLE[sigma])
            style.pop("label", None)
            ax_ratio.step(edges, np.append(r, r[-1]), where="post", **style)
        ax_ratio.axhline(1.0, color="black", lw=1.0)
        ax_ratio.set_ylabel("Var / Nom")
        ax_ratio.set_xlabel("Reconstructed energy [GeV]")  # <-- adjust to your axis quantity
        ax_ratio.tick_params(direction="in", top=True, right=True, which="both")
        ax_ratio.margins(x=0)
    else:
        ax.set_xlabel("Reconstructed energy [GeV]")  # <-- adjust to your axis quantity

    fig.tight_layout()
    if ratio:
        # Re-apply spacing tight_layout tends to disturb on the shared-x panel.
        fig.subplots_adjust(hspace=0.06)

    os.makedirs(outdir, exist_ok=True)
    safe_name = dirname
    for fmt in formats:
        fig.savefig(os.path.join(outdir, f"{safe_name}.{fmt}"), dpi=300)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--root-file", required=True, help="Path to the ROOT file with the templates.")
    p.add_argument("--yaml-file", default=None,
                   help="Path to the systematics YAML config (for energy-range labels).")
    p.add_argument("--subdir", default="OffAxis0m_numuCC_numode",
                   help="Name of the nested TDirectoryFile holding Variation_* histograms.")
    p.add_argument("--prefix", default="tmplw_",
                   help="Top-level directory name prefix to select. Default 'tmplw_' picks up "
                        "ALL energy bins (e0..e29 etc). Narrow it, e.g. --prefix tmplw_e29_eb, "
                        "to plot only one true-energy bin's sweep.")
    p.add_argument("--outdir", default="plots", help="Output directory for figures.")
    p.add_argument("--formats", nargs="+", default=["pdf", "png"],
                   help="Output file formats, e.g. --formats pdf png.")
    p.add_argument("--ratio", action="store_true",
                   help="Add a variation/nominal ratio sub-panel.")
    p.add_argument("--bin-width-norm", action="store_true",
                   help="Divide bin contents by bin width (e.g. plot Events/GeV instead of "
                        "Events/bin). The ratio panel is unaffected since bin width cancels.")
    p.add_argument("--style", default=None,
                   help="Optional matplotlib style file/sheet to apply (e.g. a thesis .mplstyle).")
    args = p.parse_args()

    if args.style:
        plt.style.use(args.style)
    else:
        mpl.rcParams.update({
            "font.size": 11,
            "axes.linewidth": 1.0,
            "xtick.minor.visible": True,
            "ytick.minor.visible": True,
            "figure.dpi": 120,
        })

    mapping = load_yaml_mapping(args.yaml_file) if args.yaml_file else {}

    with uproot.open(args.root_file) as rf:
        dirnames = discover_directories(rf, args.prefix)
        if not dirnames:
            print(f"No directories found with prefix '{args.prefix}' in {args.root_file}",
                  file=sys.stderr)
            sys.exit(1)

        print(f"Found {len(dirnames)} directories to plot.")
        for dirname in dirnames:
            subpath = f"{dirname}/{args.subdir}"
            try:
                subdir_obj = rf[subpath]
            except KeyError:
                print(f"  [skip] {subpath} not found", file=sys.stderr)
                continue

            hists = {}
            for var_name in VARIATION_SIGMA:
                try:
                    hists[var_name] = subdir_obj[var_name]
                except KeyError:
                    pass

            if not hists:
                print(f"  [skip] no Variation_* histograms in {subpath}", file=sys.stderr)
                continue

            label = make_label(dirname, mapping)
            plot_one(dirname, hists, label, args.outdir, args.formats, args.ratio, args.subdir,
                     bin_width_norm=args.bin_width_norm)
            print(f"  [ok]  {dirname}  ->  {args.outdir}/{dirname}.{{{','.join(args.formats)}}}")

    print("Done.")


if __name__ == "__main__":
    main()