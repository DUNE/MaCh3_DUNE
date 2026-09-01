#!/usr/bin/env python3
"""Scan PDSP injected Generator normalisations and summarise fit recovery."""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import json
import shutil
import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[0]
sys.path.insert(0, str(SCRIPT_DIR))

from run_pdsp_generator_consistency import (  # noqa: E402
    DEFAULT_DIAG_CONFIG,
    DEFAULT_FIT_CONFIG,
    DEFAULT_XSEC_CONFIG,
    available_table,
    check_apps,
    read_lines,
    resolve_apps,
    run_case,
)


SUMMARY_MACRO = SCRIPT_DIR / "summarise_pdsp_fit_recovery.C"


def parse_csv(raw: str) -> list[str]:
    return [item.strip() for item in raw.split(",") if item.strip()]


def scan_values(start: float, stop: float, step: float) -> list[float]:
    if step <= 0:
        raise ValueError("--step must be positive")
    values = []
    current = start
    epsilon = step / 1000.0
    while current <= stop + epsilon:
        values.append(round(current, 10))
        current += step
    return values


def value_label(value: float) -> str:
    return f"{value:g}".replace("-", "m").replace(".", "p")


def run_recovery_summary(
    *,
    root_exe: str,
    fit_output: Path,
    xsec_config: Path,
    output_csv: Path,
    label: str,
    targets: list[str],
    burn_in: int,
    tolerance: float,
    log_file: Path,
    dry_run: bool,
) -> dict[str, object]:
    targets_csv = ",".join(targets)
    macro_call = (
        f'{SUMMARY_MACRO}("'
        f'{fit_output}","{xsec_config}","{output_csv}","{label}",'
        f'"{targets_csv}",{burn_in},{tolerance})'
    )
    command = [root_exe, "-l", "-b", "-q", macro_call]
    if dry_run:
        return {
            "command": command,
            "log_file": str(log_file),
            "status": "dry-run",
            "returncode": None,
        }

    with log_file.open("w", encoding="utf-8") as log:
        completed = subprocess.run(
            command,
            cwd=REPO_ROOT,
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
    return {
        "command": command,
        "log_file": str(log_file),
        "status": "ok" if completed.returncode == 0 else "failed",
        "returncode": completed.returncode,
    }


def read_recovery_rows(csv_path: Path) -> list[dict[str, str]]:
    with csv_path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_combined_summary(output_dir: Path, rows: list[dict[str, object]], tolerance: float) -> None:
    summary_path = output_dir / "performance_summary.csv"
    if not rows:
        return

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
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    safe_rows = []
    by_target: dict[str, list[dict[str, object]]] = {}
    for row in rows:
        by_target.setdefault(str(row["scan_target"]), []).append(row)

    for target, target_rows in by_target.items():
        values = sorted({float(row["injected_value"]) for row in target_rows})
        safe_limit = None
        first_failed = None
        for value in values:
            value_rows = [row for row in target_rows if float(row["injected_value"]) == value]
            passed = all(str(row["pass"]).lower() == "true" for row in value_rows)
            if passed:
                safe_limit = value
            elif first_failed is None:
                first_failed = value
        safe_rows.append({
            "scan_target": target,
            "tolerance": tolerance,
            "safe_limit": "" if safe_limit is None else f"{safe_limit:g}",
            "first_failed": "" if first_failed is None else f"{first_failed:g}",
        })

    with (output_dir / "safe_region_summary.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["scan_target", "tolerance", "safe_limit", "first_failed"])
        writer.writeheader()
        writer.writerows(safe_rows)


def run_scan_case(
    *,
    target: str,
    target_type: str,
    value: float,
    fit_template: list[str],
    xsec_template: list[str],
    output_dir: Path,
    apps: dict[str, Path | str],
    diag_config: Path,
    workflow: str,
    extra_args: list[str],
    no_summarise: bool,
    root_exe: str,
    burn_in: int,
    tolerance: float,
    dry_run: bool,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    label = f"{target}_scan_{value_label(value)}"
    process_values = {target: value} if target_type == "process" else {}
    parameter_values = {target: value} if target_type == "parameter" else {}
    manifest = run_case(
        label=label,
        fit_template=fit_template,
        xsec_template=xsec_template,
        output_dir=output_dir,
        apps=apps,
        diag_config=diag_config,
        workflow=workflow,
        process_values=process_values,
        parameter_values=parameter_values,
        extra_args=extra_args,
        dry_run=dry_run,
    )
    manifest["scan_target"] = target
    manifest["scan_target_type"] = target_type
    manifest["injected_value"] = value

    recovery_rows: list[dict[str, object]] = []
    if no_summarise or manifest["status"] not in ("ok", "dry-run"):
        return manifest, recovery_rows

    case_dir = output_dir / label
    recovery_csv = case_dir / f"{label}_recovery.csv"
    summary_result = run_recovery_summary(
        root_exe=root_exe,
        fit_output=Path(str(manifest["fit_output"])),
        xsec_config=Path(str(manifest["xsec_config"])),
        output_csv=recovery_csv,
        label=label,
        targets=[target],
        burn_in=burn_in,
        tolerance=tolerance,
        log_file=case_dir / f"{label}_recovery.log",
        dry_run=dry_run,
    )
    manifest["recovery_summary"] = str(recovery_csv)
    manifest["recovery_stage"] = summary_result
    if summary_result["status"] == "ok":
        for row in read_recovery_rows(recovery_csv):
            row["scan_target"] = target
            row["injected_value"] = value
            row["status"] = manifest["status"]
            recovery_rows.append(row)

    return manifest, recovery_rows


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Scan injected PDSP Generator normalisations up to 2x and summarise fit recovery."
    )
    parser.add_argument("--fit-config", default=DEFAULT_FIT_CONFIG, type=Path)
    parser.add_argument("--xsec-config", default=DEFAULT_XSEC_CONFIG, type=Path)
    parser.add_argument("--diag-config", default=DEFAULT_DIAG_CONFIG, type=Path)
    parser.add_argument("--output-dir", default=Path("PDSPFitPerformanceScan"), type=Path)
    parser.add_argument("--process", action="append",
                        help="Process to scan, e.g. Abs, CEx, Pion. Can be repeated.")
    parser.add_argument("--parameter", action="append",
                        help="Exact PDSPFitModel parameter to scan. Can be repeated.")
    parser.add_argument("--start", type=float, default=1.1)
    parser.add_argument("--stop", type=float, default=2.0)
    parser.add_argument("--step", type=float, default=0.1)
    parser.add_argument("--values", help="Comma-separated explicit injected values, overrides start/stop/step.")
    parser.add_argument("--include-nominal", action="store_true",
                        help="Also include injected value 1.0 in the scan.")
    parser.add_argument("--workflow", choices=["fit", "full"], default="fit",
                        help="Use 'fit' for recovery scans, or 'full' to also make predictive plots.")
    parser.add_argument("--jobs", type=int, default=1,
                        help="Number of local scan points to run concurrently.")
    parser.add_argument("--burn-in", type=int, default=10000)
    parser.add_argument("--tolerance", type=float, default=0.10,
                        help="Allowed absolute relative bias before a point is marked failed.")
    parser.add_argument("--no-summarise", action="store_true",
                        help="Run fits only; skip ROOT posterior recovery summaries.")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--extra-arg", action="append", default=[],
                        help="Extra argument passed through to Fit.")
    args = parser.parse_args()

    fit_config = args.fit_config if args.fit_config.is_absolute() else REPO_ROOT / args.fit_config
    xsec_config = args.xsec_config if args.xsec_config.is_absolute() else REPO_ROOT / args.xsec_config
    diag_config = args.diag_config if args.diag_config.is_absolute() else REPO_ROOT / args.diag_config
    output_dir = args.output_dir if args.output_dir.is_absolute() else REPO_ROOT / args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    fit_template = read_lines(fit_config)
    xsec_template = read_lines(xsec_config)

    known_processes = sorted({process for process, _ in available_table(xsec_template)})
    targets = args.process or []
    parameters = args.parameter or []
    if not targets and not parameters:
        targets = [process for process in known_processes if process != "Impure"]

    unknown_processes = sorted(set(targets) - set(known_processes))
    if unknown_processes:
        raise ValueError(f"Unknown process(es): {', '.join(unknown_processes)}")

    values = [float(value) for value in parse_csv(args.values)] if args.values else scan_values(args.start, args.stop, args.step)
    if not args.include_nominal:
        values = [value for value in values if abs(value - 1.0) > 1e-12]
    elif all(abs(value - 1.0) > 1e-12 for value in values):
        values = [1.0, *values]

    apps = resolve_apps()
    check_apps(apps)
    root_exe = shutil.which("root") or "root"

    manifests = []
    recovery_rows: list[dict[str, object]] = []
    scan_targets = [(target, "process") for target in targets] + [(parameter, "parameter") for parameter in parameters]

    tasks = [(target, target_type, value) for target, target_type in scan_targets for value in values]
    if args.jobs < 1:
        raise ValueError("--jobs must be at least 1")

    if args.jobs == 1:
        for target, target_type, value in tasks:
            manifest, rows = run_scan_case(
                target=target,
                target_type=target_type,
                value=value,
                fit_template=fit_template,
                xsec_template=xsec_template,
                output_dir=output_dir,
                apps=apps,
                diag_config=diag_config,
                workflow=args.workflow,
                extra_args=args.extra_arg,
                no_summarise=args.no_summarise,
                root_exe=root_exe,
                burn_in=args.burn_in,
                tolerance=args.tolerance,
                dry_run=args.dry_run,
            )
            manifests.append(manifest)
            recovery_rows.extend(rows)
            print(f"{manifest['label']}: {manifest['status']} -> {manifest['fit_output']}")
    else:
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as executor:
            futures = [
                executor.submit(
                    run_scan_case,
                    target=target,
                    target_type=target_type,
                    value=value,
                    fit_template=fit_template,
                    xsec_template=xsec_template,
                    output_dir=output_dir,
                    apps=apps,
                    diag_config=diag_config,
                    workflow=args.workflow,
                    extra_args=args.extra_arg,
                    no_summarise=args.no_summarise,
                    root_exe=root_exe,
                    burn_in=args.burn_in,
                    tolerance=args.tolerance,
                    dry_run=args.dry_run,
                )
                for target, target_type, value in tasks
            ]
            for future in concurrent.futures.as_completed(futures):
                manifest, rows = future.result()
                manifests.append(manifest)
                recovery_rows.extend(rows)
                print(f"{manifest['label']}: {manifest['status']} -> {manifest['fit_output']}")

    with (output_dir / "manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifests, handle, indent=2, sort_keys=True)
        handle.write("\n")

    if recovery_rows:
        write_combined_summary(output_dir, recovery_rows, args.tolerance)

    failed = [manifest for manifest in manifests if manifest["status"] == "failed"]
    return 1 if failed else 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(1)
