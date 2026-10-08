#!/usr/bin/env python3
"""Plot PDSP fit-performance scan degradation from performance_summary.csv."""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = Path(__file__).resolve().parent
PLOT_MACRO = SCRIPT_DIR / "plot_pdsp_fit_performance_scan.C"


def root_quote(value: Path | str) -> str:
    return str(value).replace("\\", "\\\\").replace('"', '\\"')


def infer_scan_target(manifest: dict[str, object], row: dict[str, str]) -> tuple[str, float]:
    process_values = manifest.get("process_values", {})
    if isinstance(process_values, dict) and process_values:
        target, value = next(iter(process_values.items()))
        return str(target), float(value)

    parameter_values = manifest.get("parameter_values", {})
    if isinstance(parameter_values, dict) and parameter_values:
        target, value = next(iter(parameter_values.items()))
        return str(target), float(value)

    if "scan_target" in manifest and "injected_value" in manifest:
        return str(manifest["scan_target"]), float(manifest["injected_value"])

    return row.get("process") or row.get("parameter") or "unknown", float(row["injected"])


def build_summary_from_recovery_files(summary: Path) -> bool:
    scan_dir = summary.parent
    recovery_files = sorted(scan_dir.glob("*/*_recovery.csv"))
    if not recovery_files:
        return False

    rows: list[dict[str, object]] = []
    for recovery_file in recovery_files:
        manifest_file = recovery_file.parent / "manifest.json"
        manifest: dict[str, object] = {}
        if manifest_file.exists():
            with manifest_file.open("r", encoding="utf-8") as handle:
                manifest = json.load(handle)

        with recovery_file.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            for recovery_row in reader:
                target, injected_value = infer_scan_target(manifest, recovery_row)
                row: dict[str, object] = {
                    "scan_target": target,
                    "injected_value": injected_value,
                    "status": manifest.get("status", "unknown"),
                }
                row.update(recovery_row)
                rows.append(row)

    if not rows:
        return False

    fieldnames = [
        "scan_target",
        "injected_value",
        "label",
        "status",
        "parameter",
        "process",
        "injected",
        "posterior_mean",
        "posterior_rms",
        "abs_error",
        "rel_error",
        "pull",
        "entries_used",
        "tolerance",
        "pass",
    ]
    with summary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Built missing summary from {len(recovery_files)} recovery file(s): {summary}")
    return True


def summary_targets(summary: Path) -> set[str]:
    if not summary.exists():
        return set()
    with summary.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return {row.get("scan_target", "") for row in reader if row.get("scan_target", "")}


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Make discrepancy/degradation plots from a PDSP fit performance scan."
    )
    parser.add_argument("--summary", default=Path("PDSPFitPerformanceScan/performance_summary.csv"), type=Path,
                        help="Input CSV from run_pdsp_fit_performance_scan.py")
    parser.add_argument("--output-dir", default=Path("PDSPFitPerformanceScan/plots"), type=Path)
    parser.add_argument("--tolerance", type=float, default=0.10,
                        help="Horizontal reference line for acceptable relative discrepancy")
    parser.add_argument("--target", default="",
                        help="Optional scan target to plot, e.g. Abs, CEx, Pion, or Abs_TrueEBin_0")
    parser.add_argument("--rebuild-summary", action="store_true",
                        help="Rebuild performance_summary.csv from per-point *_recovery.csv files before plotting")
    parser.add_argument("--root", default=shutil.which("root") or "root")
    args = parser.parse_args()

    summary = args.summary if args.summary.is_absolute() else REPO_ROOT / args.summary
    output_dir = args.output_dir if args.output_dir.is_absolute() else REPO_ROOT / args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    if args.rebuild_summary:
        build_summary_from_recovery_files(summary)

    if not summary.exists() and not build_summary_from_recovery_files(summary):
        raise FileNotFoundError(f"Missing performance summary CSV: {summary}")

    if args.target and args.target not in summary_targets(summary):
        if build_summary_from_recovery_files(summary):
            if args.target not in summary_targets(summary):
                raise ValueError(f"No rows found for target '{args.target}' in {summary}")

    macro_call = (
        f'{PLOT_MACRO}("'
        f'{root_quote(summary)}","{root_quote(output_dir)}",'
        f'{args.tolerance},"{root_quote(args.target)}")'
    )
    command = [args.root, "-l", "-b", "-q", macro_call]
    completed = subprocess.run(command, cwd=REPO_ROOT, check=False)
    return completed.returncode


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(1)
