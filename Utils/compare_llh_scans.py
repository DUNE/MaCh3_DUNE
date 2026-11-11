"""
compare_llh_scans.py

Compare 1D likelihood scans (TH1D) between two MaCh3 ROOT files.

Outputs:
- Individual plots
- One multipage PDF per LLH directory
- CSV summary
"""

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


# ------------------------------------------------------------------
# Style
# ------------------------------------------------------------------

C_COV = "#1f77b4"
C_NOCOV = "#d62728"
LW = 2.0

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 12,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linestyle": "--",
})

DIRS_DEFAULT = [
    "xsec_LLH",
    "Sample_LLH",
    "Total_LLH",
]


# ------------------------------------------------------------------
# ROOT helpers
# ------------------------------------------------------------------

def th1_to_xy(h):
    vals = h.values()
    edges = h.axis().edges()
    centres = 0.5 * (edges[:-1] + edges[1:])
    return centres, vals


def load_dir(filepath, dirname):

    import uproot

    result = {}

    with uproot.open(filepath) as f:

        if dirname not in f:
            print(f"[warn] {dirname} missing")
            return result

        d = f[dirname]

        for key in d.keys():

            name = key.split(";")[0]

            try:
                obj = d[key]

                if "TH1" not in obj.classname:
                    continue

                param = name.removesuffix("_full")

                result[param] = th1_to_xy(obj)

            except Exception as e:
                print(f"skip {name}: {e}")

    return result


# ------------------------------------------------------------------
# Plotting
# ------------------------------------------------------------------
def make_figure(
    param,
    x_cov,
    y_cov,
    x_nocov,
    y_nocov,
    dirname,
):
    fig, (ax_cov, ax_nocov, ax_overlay, ax_diff) = plt.subplots(
        1,
        4,
        figsize=(20, 4.5),
        gridspec_kw={"width_ratios": [2, 2, 3, 2]},
    )

    # --- cov only ---
    ax_cov.plot(x_cov, y_cov, color=C_COV, lw=LW)
    ax_cov.axhline(1.0, color="grey", ls=":")
    ax_cov.axhline(4.0, color="grey", ls="-.")
    ax_cov.set_xlabel("Parameter value")
    ax_cov.set_ylabel(r"$\Delta(-2\ln\mathcal{L})$")
    ax_cov.set_title(f"{param}\n(with covariance)")

    # --- nocov only ---
    ax_nocov.plot(x_nocov, y_nocov, color=C_NOCOV, lw=LW, linestyle="--")
    ax_nocov.axhline(1.0, color="grey", ls=":")
    ax_nocov.axhline(4.0, color="grey", ls="-.")
    ax_nocov.set_xlabel("Parameter value")
    ax_nocov.set_ylabel(r"$\Delta(-2\ln\mathcal{L})$")
    ax_nocov.set_title(f"{param}\n(without  covariance)")

    # --- overlay ---
    ax_overlay.plot(x_cov, y_cov, color=C_COV, lw=LW, label="with covariance")
    ax_overlay.plot(x_nocov, y_nocov, color=C_NOCOV, lw=LW, linestyle="--", label="without covariance")
    ax_overlay.axhline(1.0, color="grey", ls=":")
    ax_overlay.axhline(4.0, color="grey", ls="-.")
    ax_overlay.set_xlabel("Parameter value")
    ax_overlay.set_ylabel(r"$\Delta(-2\ln\mathcal{L})$")
    ax_overlay.set_title(f"{param}\n(overlay)")
    ax_overlay.legend(frameon=False)
    ax_overlay.text(
        0.98, 0.97, dirname,
        transform=ax_overlay.transAxes,
        ha="right", va="top", fontsize=9, color="grey",
    )

    # --- difference ---
    interp = np.interp(x_cov, x_nocov, y_nocov, left=np.nan, right=np.nan)
    diff = y_cov - interp

    ax_diff.plot(x_cov, diff, color="black")
    ax_diff.axhline(0, color="black", ls="--")
    ax_diff.fill_between(
        x_cov, diff, 0,
        where=(diff > 0), color=C_COV, alpha=0.25,
        label="Covariance larger ΔLLH",
    )
    ax_diff.fill_between(
        x_cov, diff, 0,
        where=(diff < 0), color=C_NOCOV, alpha=0.25,
        label="No covariance larger ΔLLH",
    )
    ax_diff.set_xlabel("Parameter value")
    ax_diff.set_ylabel(r"$\Delta_{\rm cov}-\Delta_{\rm no\,cov}$")
    ax_diff.set_title(f"{param}\n(difference)")
    ax_diff.legend(frameon=False, fontsize=9)

    fig.tight_layout()

    return fig


# ------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------

def best_fit_and_1sigma(x, y):
    idx = np.argmin(y)
    best = x[idx]
    # find crossings of y=1
    above = np.where(y <= 1.0)[0]
    if len(above) > 0:
        lo = x[above[0]]
        hi = x[above[-1]]
    else:
        lo, hi = np.nan, np.nan
    return best, lo, hi

# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------

def main():

    ap = argparse.ArgumentParser()

    ap.add_argument("--with_cov", required=True)
    ap.add_argument("--without_cov", required=True)

    ap.add_argument(
        "--outdir",
        default="./llh_plots",
    )

    ap.add_argument(
        "--dirs",
        nargs="+",
        default=DIRS_DEFAULT,
    )

    ap.add_argument(
        "--fmt",
        default="pdf",
    )

    ap.add_argument(
        "--multipage",
        action="store_true",
    )

    args = ap.parse_args()

    try:
        import uproot
    except ImportError:
        sys.exit("Install uproot")

    outdir = Path(args.outdir)

    outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    summary = []

    total = 0

    for dirname in args.dirs:

        print(f"\nLoading {dirname}")

        pdf = None

        if args.multipage:

            pdf = PdfPages(
                outdir /
                f"{dirname}_all.pdf"
            )

        cov = load_dir(
            args.with_cov,
            dirname,
        )

        nocov = load_dir(
            args.without_cov,
            dirname,
        )

        common = sorted(
            set(cov)
            &
            set(nocov)
        )

        ddir = outdir / dirname

        ddir.mkdir(
            exist_ok=True
        )

        for param in common:

            fig = make_figure(
                param,
                *cov[param],
                *nocov[param],
                dirname,
            )

            path = (
                ddir
                /
                f"{param}.{args.fmt}"
            )

            fig.savefig(
                path,
                bbox_inches="tight",
            )

            if pdf:
                pdf.savefig(fig)

            plt.close(fig)

            summary.append({
                "directory": dirname,
                "parameter": param,
                "bestfit_cov":
                    best_fit_and_1sigma(
                        *cov[param]
                    ),
                "bestfit_nocov":
                    best_fit_and_1sigma(
                        *nocov[param]
                    ),
            })

            total += 1

        if pdf:
            pdf.close()

    csv_path = (
        outdir
        /
        "llh_comparison_summary.csv"
    )

    with open(
        csv_path,
        "w",
        newline="",
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=summary[0].keys(),
        )

        writer.writeheader()
        writer.writerows(summary)

    print(f"\nDone — {total} plots")


if __name__ == "__main__":
    main()