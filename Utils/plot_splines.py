#!/usr/bin/env python3
"""
plot_splines.py

Read a MaCh3-style systematics YAML config plus its corresponding ROOT
spline file, and produce one publication-quality figure per cross-section
parameter, overlaying every associated TSpline3 (i.e. every mode/sample/bin
belonging to that parameter) on a single set of axes.

Usage
-----
    python plot_splines.py --yaml config.yaml --root splines.root --outdir plots/

    # Only plot specific parameters (by FancyName or SplineName):
    python plot_splines.py --yaml config.yaml --root splines.root \
        --params FrInel_N MaCCQE

    # Restrict to a ParameterGroup:
    python plot_splines.py --yaml config.yaml --root splines.root --group Xsec

Dependencies
------------
    pip install uproot pyyaml matplotlib numpy --break-system-packages

Notes on matching splines to parameters
----------------------------------------
Each systematic in the YAML has a short `SplineInformation.SplineName`
(e.g. "frineln"). In the ROOT file, spline keys are typically named like:

    dev_<SplineName>_<something>_sp_<paramindex>_<binnumber>

e.g.  dev_frineln_ccqe_sp_2_12

This script finds every TSpline3 key in the ROOT file whose name contains
"_<SplineName>_" (or starts with "<SplineName>_" / "dev_<SplineName>")
and groups them together under that parameter's FancyName. If your file
uses a different naming convention, adjust `SPLINE_NAME_PATTERN` below or
pass --key-pattern to override it.
"""

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize

try:
    import uproot
except ImportError:
    sys.exit(
        "uproot is required. Install with:\n"
        "    pip install uproot --break-system-packages"
    )


# ----------------------------------------------------------------------
# Thesis-friendly matplotlib defaults
# ----------------------------------------------------------------------
plt.rcParams.update({
    "font.family": "serif",
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "legend.fontsize": 8,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "axes.axisbelow": True,
})


def load_yaml_params(yaml_path, group_filter=None):
    """Parse the MaCh3 YAML config and return a list of dicts describing
    each systematic: {fancy_name, param_name, spline_name, group}."""
    with open(yaml_path, "r") as f:
        cfg = yaml.safe_load(f)

    params = []
    for entry in cfg.get("Systematics", []):
        sysdict = entry.get("Systematic", entry)  # tolerate either shape
        spline_info = sysdict.get("SplineInformation", {})
        spline_name = spline_info.get("SplineName")
        if not spline_name:
            continue  # not a spline-type systematic

        names = sysdict.get("Names", {})
        fancy = names.get("FancyName", spline_name)
        pname = names.get("ParameterName", fancy)
        group = sysdict.get("ParameterGroup", "Unknown")

        if group_filter and group != group_filter:
            continue

        params.append({
            "fancy_name": fancy,
            "param_name": pname,
            "spline_name": spline_name,
            "group": group,
            "modes": spline_info.get("Mode", []),
        })
    return params


def find_matching_keys(root_keys, spline_name):
    """Return all ROOT key names that belong to a given short SplineName.

    Matches keys like:
        dev_<spline_name>_..._sp_<i>_<j>
        <spline_name>_...
    using a boundary-safe regex so e.g. 'maqe' doesn't match 'maqeshape'.
    """
    pattern = re.compile(rf"(^|_){re.escape(spline_name)}(_|$)", re.IGNORECASE)
    return [k for k in root_keys if pattern.search(k)]


def _get_member(obj, name):
    """Try several ways of pulling a named member off an uproot dynamic
    object, since the API differs slightly between record-array items
    and top-level model objects."""
    try:
        return obj.member(name)
    except Exception:
        pass
    try:
        return obj[name]
    except Exception:
        pass
    try:
        return getattr(obj, name)
    except Exception:
        pass
    return None


def get_spline_xy(spline_obj, n_points=400):
    """Evaluate a TSpline3 (via uproot) over its knot range.

    Newer uproot TSpline3 models (Model_TSpline3_v2) do NOT expose flat
    fXvals/fYvals/fB/fC/fD arrays. Instead the knots and per-interval
    cubic coefficients live inside `fPoly`, an array of TSplinePoly3
    records, each with members fX, fY, fB, fC, fD such that for
    dx = x - fX_i (with fX_i <= x < fX_{i+1}):

        y(x) = fY_i + fB_i*dx + fC_i*dx^2 + fD_i*dx^3

    fNp gives the number of knots; fXmin/fXmax give the overall range.
    """
    poly = _get_member(spline_obj, "fPoly")
    if poly is None or len(poly) == 0:
        raise RuntimeError("Spline has no 'fPoly' knot array; cannot evaluate.")

    def col(name):
        # poly may be an awkward array of records, a list of objects,
        # or a numpy structured array - try each access style.
        try:
            return np.asarray(poly[name])
        except Exception:
            pass
        return np.asarray([_get_member(p, name) for p in poly])

    x_knots = col("fX")
    y_knots = col("fY")
    b = col("fB")
    c = col("fC")
    d = col("fD")

    if x_knots is None or len(x_knots) == 0:
        raise RuntimeError("Could not extract 'fX' knot positions from fPoly.")

    order = np.argsort(x_knots)
    x_knots, y_knots, b, c, d = (arr[order] for arr in (x_knots, y_knots, b, c, d))

    x_lo = float(_get_member(spline_obj, "fXmin") or x_knots[0])
    x_hi = float(_get_member(spline_obj, "fXmax") or x_knots[-1])
    x_dense = np.linspace(x_lo, x_hi, n_points)

    y_dense = np.empty_like(x_dense)
    n_knots = len(x_knots)
    for i, xv in enumerate(x_dense):
        idx = np.searchsorted(x_knots, xv, side="right") - 1
        idx = int(np.clip(idx, 0, n_knots - 2 if n_knots > 1 else 0))
        dx = xv - x_knots[idx]
        y_dense[i] = y_knots[idx] + b[idx] * dx + c[idx] * dx**2 + d[idx] * dx**3

    return x_dense, y_dense


def sanitize_filename(name):
    return re.sub(r"[^\w\-.]", "_", name)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--yaml", required=True, help="Path to MaCh3 systematics YAML config")
    ap.add_argument("--root", required=True, help="Path to ROOT file containing TSpline3 objects")
    ap.add_argument("--outdir", default="spline_plots", help="Output directory for figures")
    ap.add_argument("--group", default=None, help="Only plot parameters with this ParameterGroup (e.g. Xsec)")
    ap.add_argument("--params", nargs="*", default=None,
                     help="Only plot these parameters (match FancyName, ParameterName, or SplineName)")
    ap.add_argument("--cmap", default="turbo", help="Matplotlib colormap for the many curves (default: turbo)")
    ap.add_argument("--format", default="pdf", choices=["pdf", "png", "both"],
                     help="Output file format (default: pdf, thesis-ready vector)")
    ap.add_argument("--max-curves", type=int, default=None,
                     help="Optional cap on number of splines drawn per parameter (useful if there are hundreds)")
    ap.add_argument("--debug", action="store_true",
                     help="On the first spline-evaluation failure, dump the raw fPoly structure to help diagnose the ROOT/uproot format")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Reading YAML config: {args.yaml}")
    params = load_yaml_params(args.yaml, group_filter=args.group)
    if args.params:
        wanted = set(p.lower() for p in args.params)
        params = [
            p for p in params
            if p["fancy_name"].lower() in wanted
            or p["param_name"].lower() in wanted
            or p["spline_name"].lower() in wanted
        ]
    print(f"Found {len(params)} spline-type parameter(s) to process.")

    print(f"Opening ROOT file: {args.root}")
    f = uproot.open(args.root)
    all_keys = [k.split(";")[0] for k in f.keys()]  # strip ROOT cycle numbers e.g. ';1'
    all_keys = list(dict.fromkeys(all_keys))  # dedupe, keep order

    n_ok, n_empty, n_fail = 0, 0, 0
    debug_dumped = False

    for p in params:
        fancy = p["fancy_name"]
        spline_name = p["spline_name"]
        matches = find_matching_keys(all_keys, spline_name)

        if not matches:
            print(f"  [!] No ROOT keys found for '{fancy}' (SplineName='{spline_name}') - skipping")
            n_empty += 1
            continue

        if args.max_curves and len(matches) > args.max_curves:
            matches = matches[: args.max_curves]

        fig, ax = plt.subplots(figsize=(7, 5))
        cmap = cm.get_cmap(args.cmap)
        norm = Normalize(vmin=0, vmax=max(len(matches) - 1, 1))

        n_drawn = 0
        for i, key in enumerate(matches):
            try:
                spline_obj = f[key]
                x, y = get_spline_xy(spline_obj)
            except Exception as e:
                print(f"      could not evaluate {key}: {e}")
                if args.debug and not debug_dumped:
                    debug_dumped = True
                    print(f"      --- DEBUG dump for {key} ---")
                    try:
                        poly = _get_member(spline_obj, "fPoly")
                        print(f"      type(fPoly) = {type(poly)}")
                        print(f"      len(fPoly)  = {len(poly)}")
                        first = poly[0]
                        print(f"      type(fPoly[0]) = {type(first)}")
                        print(f"      dir(fPoly[0])  = {[a for a in dir(first) if not a.startswith('_')]}")
                        try:
                            print(f"      fPoly[0] fields = {first.fields}")
                        except Exception:
                            pass
                        try:
                            print(f"      fPoly[0].tolist() = {first.tolist()}")
                        except Exception:
                            pass
                    except Exception as debug_e:
                        print(f"      (debug dump itself failed: {debug_e})")
                    print(f"      --- end debug dump ---")
                continue
            ax.plot(x, y, color=cmap(norm(i)), linewidth=1.3, alpha=0.85)
            n_drawn += 1

        if n_drawn == 0:
            plt.close(fig)
            print(f"  [!] Found keys for '{fancy}' but none could be evaluated - skipping")
            n_fail += 1
            continue

        ax.set_title(f"{fancy}  ({n_drawn} splines)")
        ax.set_xlabel("Parameter variation (units of $1\\sigma$)")
        ax.set_ylabel("Weight / normalisation")
        ax.axhline(1.0, color="grey", linewidth=0.8, linestyle="--", alpha=0.6)

        sm = cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, pad=0.02)
        cbar.set_label("Spline index (bin / mode / sample)")

        fig.tight_layout()

        stem = outdir / f"spline_{sanitize_filename(fancy)}"
        if args.format in ("pdf", "both"):
            fig.savefig(f"{stem}.pdf")
        if args.format in ("png", "both"):
            fig.savefig(f"{stem}.png")
        plt.close(fig)

        print(f"  [OK] {fancy}: drew {n_drawn}/{len(matches)} splines -> {stem}.{args.format}")
        n_ok += 1

    print("\nDone.")
    print(f"  Parameters plotted:        {n_ok}")
    print(f"  Parameters with no keys:   {n_empty}")
    print(f"  Parameters that failed:    {n_fail}")
    print(f"  Output directory:          {outdir.resolve()}")


if __name__ == "__main__":
    main()