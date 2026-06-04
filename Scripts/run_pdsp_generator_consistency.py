#!/usr/bin/env python3
"""Run PDSP generator-value consistency checks from copied configs.

This script intentionally uses only the Python standard library so it works in
MaCh3 runtime environments that do not provide PyYAML.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_FIT_CONFIG = Path("Configs/FitterConfig_PDSP.yaml")
DEFAULT_XSEC_CONFIG = Path("Configs/CovObjs/PDSPFitModel.yaml")
DEFAULT_DIAG_CONFIG = Path("Configs/PDSPDiagConfig.yaml")

PARAMETER_RE = re.compile(r"^(\s*)ParameterName:\s*(\S+)\s*$")
GENERATOR_RE = re.compile(r"^(\s*)Generator:\s*([-+0-9.eE]+)\s*$")


def parse_assignment(raw: str) -> tuple[str, float]:
    if "=" not in raw:
        raise argparse.ArgumentTypeError(f"Expected NAME=VALUE, got '{raw}'")
    name, value = raw.split("=", 1)
    name = name.strip()
    if not name:
        raise argparse.ArgumentTypeError(f"Missing name in '{raw}'")
    try:
        return name, float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"Generator value must be numeric in '{raw}'") from exc


def process_name(parameter_name: str) -> str:
    return parameter_name.split("_", 1)[0]


def read_lines(path: Path) -> list[str]:
    return path.read_text(encoding="utf-8").splitlines(keepends=True)


def write_lines(path: Path, lines: list[str]) -> None:
    path.write_text("".join(lines), encoding="utf-8")


def format_float(value: float) -> str:
    return f"{value:g}"


def extract_parameters(xsec_lines: list[str]) -> list[str]:
    parameters = []
    for line in xsec_lines:
        match = PARAMETER_RE.match(line.rstrip("\n"))
        if match:
            parameters.append(match.group(2))
    return parameters


def set_generators(
    xsec_lines: list[str],
    *,
    default_value: float,
    process_values: dict[str, float],
    parameter_values: dict[str, float],
) -> tuple[list[str], dict[str, float]]:
    parameters = extract_parameters(xsec_lines)
    known_parameters = set(parameters)
    known_processes = {process_name(name) for name in parameters}

    unknown_processes = sorted(set(process_values) - known_processes)
    unknown_parameters = sorted(set(parameter_values) - known_parameters)
    if unknown_processes:
        raise ValueError(f"Unknown process(es): {', '.join(unknown_processes)}")
    if unknown_parameters:
        raise ValueError(f"Unknown parameter(s): {', '.join(unknown_parameters)}")

    changed: dict[str, float] = {}
    current_parameter: str | None = None
    output = []

    for line in xsec_lines:
        stripped = line.rstrip("\n")
        parameter_match = PARAMETER_RE.match(stripped)
        if parameter_match:
            current_parameter = parameter_match.group(2)
            output.append(line)
            continue

        generator_match = GENERATOR_RE.match(stripped)
        if generator_match and current_parameter is not None:
            value = default_value
            proc = process_name(current_parameter)
            if proc in process_values:
                value = process_values[proc]
            if current_parameter in parameter_values:
                value = parameter_values[current_parameter]
            changed[current_parameter] = float(value)
            newline = "\n" if line.endswith("\n") else ""
            output.append(f"{generator_match.group(1)}Generator: {format_float(value)}{newline}")
            continue

        output.append(line)

    missing = sorted(set(parameters) - set(changed))
    if missing:
        raise ValueError(f"No Generator line found for parameter(s): {', '.join(missing)}")

    return output, changed


def replace_simple_key(lines: list[str], key: str, value: str, indent: str = "  ") -> list[str]:
    pattern = re.compile(rf"^{re.escape(indent)}{re.escape(key)}:\s*.*$")
    replaced = False
    output = []
    for line in lines:
        if not replaced and pattern.match(line.rstrip("\n")):
            newline = "\n" if line.endswith("\n") else ""
            output.append(f"{indent}{key}: {value}{newline}")
            replaced = True
        else:
            output.append(line)
    if not replaced:
        raise ValueError(f"Could not find '{key}' in fit config")
    return output


def replace_xsec_cov_file(lines: list[str], xsec_config: Path) -> list[str]:
    output = []
    index = 0
    replaced = False
    while index < len(lines):
        line = lines[index]
        if not replaced and re.match(r"^    XsecCovFile:\s*", line):
            newline = "\n" if line.endswith("\n") else ""
            output.append(f'    XsecCovFile: ["{xsec_config}"]{newline}')
            index += 1
            while index < len(lines):
                stripped = lines[index].strip()
                if stripped == "" or stripped.startswith("#"):
                    output.append(lines[index])
                    index += 1
                    continue
                if re.match(r"^    [A-Za-z0-9_]+:", lines[index]):
                    break
                if lines[index].startswith("      ") or stripped in ("]", "]"):
                    index += 1
                    continue
                break
            replaced = True
            continue
        output.append(line)
        index += 1
    if not replaced:
        raise ValueError("Could not find 'XsecCovFile' in fit config")
    return output


def make_run_config(fit_lines: list[str], xsec_config: Path, output_file: Path) -> list[str]:
    lines = replace_simple_key(fit_lines, "OutputFile", f'"{output_file}"')
    lines = replace_simple_key(lines, "Data", "false")
    lines = replace_simple_key(lines, "StatOnly", "true")
    lines = replace_xsec_cov_file(lines, xsec_config)
    lines = replace_simple_key(lines, "XsecAsimovTune", '"Generator"', indent="    ")
    return lines


def app_path(raw: str) -> Path | str:
    candidate = Path(raw)
    if candidate.is_absolute() or "/" in raw:
        return candidate

    for relative in (Path("build/bin") / raw, Path("build/Apps") / raw):
        full = REPO_ROOT / relative
        if full.exists():
            return full
    return raw


def run_command(command: list[str], log_file: Path, dry_run: bool, cwd: Path = REPO_ROOT) -> dict[str, object]:
    if dry_run:
        return {
            "command": command,
            "cwd": str(cwd),
            "log_file": str(log_file),
            "returncode": None,
            "status": "dry-run",
        }

    with log_file.open("w", encoding="utf-8") as log:
        completed = subprocess.run(
            command,
            cwd=cwd,
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )

    return {
        "command": command,
        "cwd": str(cwd),
        "log_file": str(log_file),
        "returncode": completed.returncode,
        "status": "ok" if completed.returncode == 0 else "failed",
    }


def resolve_apps() -> dict[str, Path | str]:
    return {
        "Fit": app_path("Fit"),
        "ProcessMCMC": app_path("ProcessMCMC"),
        "PredictivePDSP": app_path("PredictivePDSP"),
        "PredictivePlotting": app_path("PredictivePlotting"),
    }


def check_apps(apps: dict[str, Path | str]) -> None:
    missing = []
    for name, path in apps.items():
        if isinstance(path, Path) and not path.exists():
            missing.append(f"{name}: {path}")
    if missing:
        raise FileNotFoundError("Missing executable(s): " + ", ".join(missing))


def run_case(
    *,
    label: str,
    fit_template: list[str],
    xsec_template: list[str],
    output_dir: Path,
    apps: dict[str, Path | str],
    diag_config: Path,
    workflow: str,
    process_values: dict[str, float],
    parameter_values: dict[str, float],
    extra_args: list[str],
    dry_run: bool,
) -> dict[str, object]:
    case_dir = output_dir / label
    case_dir.mkdir(parents=True, exist_ok=True)

    xsec_path = case_dir / "PDSPFitModel.yaml"
    fit_path = case_dir / "FitterConfig_PDSP.yaml"
    fit_output = case_dir / f"{label}.root"
    posterior_output = case_dir / f"{label}_posterior_predictive.root"
    prior_output = case_dir / f"{label}_prior_predictive.root"
    overlay_outputs = [
        case_dir / "Overlay_Predictive.pdf",
        case_dir / "Overlay_Predictive_norm.pdf",
        case_dir / "Overlay_Predictive_pi.pdf",
    ]

    xsec_lines, generators = set_generators(
        list(xsec_template),
        default_value=1.0,
        process_values=process_values,
        parameter_values=parameter_values,
    )
    fit_lines = make_run_config(list(fit_template), xsec_path, fit_output)

    write_lines(xsec_path, xsec_lines)
    write_lines(fit_path, fit_lines)

    stages: list[tuple[str, list[str], Path, Path]] = [
        ("fit", [str(apps["Fit"]), str(fit_path), *extra_args], case_dir / f"{label}_fit.log", REPO_ROOT),
    ]
    if workflow == "full":
        stages.extend([
            (
                "process_mcmc",
                [str(apps["ProcessMCMC"]), str(diag_config), str(fit_output)],
                case_dir / f"{label}_process_mcmc.log",
                REPO_ROOT,
            ),
            (
                "posterior_predictive",
                [
                    str(apps["PredictivePDSP"]),
                    str(fit_path),
                    f"General:OutputFile:{posterior_output}",
                    f"Predictive:PosteriorFile:{fit_output}",
                    "Predictive:PriorPredictive:False",
                ],
                case_dir / f"{label}_posterior_predictive.log",
                REPO_ROOT,
            ),
            (
                "prior_predictive",
                [
                    str(apps["PredictivePDSP"]),
                    str(fit_path),
                    f"General:OutputFile:{prior_output}",
                    f"Predictive:PosteriorFile:{fit_output}",
                    "Predictive:PriorPredictive:True",
                ],
                case_dir / f"{label}_prior_predictive.log",
                REPO_ROOT,
            ),
            (
                "predictive_plotting",
                [str(apps["PredictivePlotting"]), str(diag_config), str(posterior_output), str(prior_output)],
                case_dir / f"{label}_predictive_plotting.log",
                case_dir,
            ),
        ])

    stage_results = []
    for stage_name, command, log_file, cwd in stages:
        result = run_command(command, log_file, dry_run, cwd)
        result["name"] = stage_name
        stage_results.append(result)
        if result["status"] == "failed":
            break

    failed = [stage for stage in stage_results if stage["status"] == "failed"]
    status = "failed" if failed else stage_results[-1]["status"]
    returncode = failed[0]["returncode"] if failed else stage_results[-1]["returncode"]

    manifest = {
        "label": label,
        "status": status,
        "returncode": returncode,
        "fit_config": str(fit_path),
        "xsec_config": str(xsec_path),
        "fit_output": str(fit_output),
        "posterior_predictive_output": str(posterior_output),
        "prior_predictive_output": str(prior_output),
        "overlay_outputs": [str(path) for path in overlay_outputs],
        "workflow": workflow,
        "stages": stage_results,
        "process_values": process_values,
        "parameter_values": parameter_values,
        "generators": generators,
    }
    with (case_dir / "manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest


def available_table(xsec_lines: list[str]) -> list[tuple[str, str]]:
    return [(process_name(name), name) for name in extract_parameters(xsec_lines)]


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Run PDSP Generator-varied fake-data checks using copied "
            "FitterConfig/PDSPFitModel YAML files."
        )
    )
    parser.add_argument("--fit-config", default=DEFAULT_FIT_CONFIG, type=Path)
    parser.add_argument("--xsec-config", default=DEFAULT_XSEC_CONFIG, type=Path)
    parser.add_argument("--diag-config", default=DEFAULT_DIAG_CONFIG, type=Path)
    parser.add_argument("--output-dir", default=Path("PDSPGeneratorConsistency"), type=Path)
    parser.add_argument("--workflow", choices=["fit", "full"], default="full",
                        help="Use 'fit' for only Fit, or 'full' for Fit plus ProcessMCMC/PredictivePDSP/PredictivePlotting")
    parser.add_argument("--set", dest="process_sets", action="append", default=[], type=parse_assignment,
                        metavar="PROCESS=VALUE", help="Set Generator for every systematic in PROCESS")
    parser.add_argument("--parameter", action="append", default=[], type=parse_assignment,
                        metavar="PARAMETER=VALUE", help="Set Generator for one exact systematic parameter")
    parser.add_argument("--all-systematics", type=float,
                        help="Run one varied case where every systematic Generator is set to this value")
    parser.add_argument("--label", default=None, help="Label for a single custom varied case")
    parser.add_argument("--nominal", action="store_true",
                        help="Also run the nominal all-Generator=1 case")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--list", action="store_true", help="List available PDSP processes and parameters")
    parser.add_argument("--extra-arg", action="append", default=[],
                        help="Extra argument passed through to the MaCh3 executable")
    args = parser.parse_args()

    fit_config = args.fit_config if args.fit_config.is_absolute() else REPO_ROOT / args.fit_config
    xsec_config = args.xsec_config if args.xsec_config.is_absolute() else REPO_ROOT / args.xsec_config
    diag_config = args.diag_config if args.diag_config.is_absolute() else REPO_ROOT / args.diag_config
    output_dir = args.output_dir if args.output_dir.is_absolute() else REPO_ROOT / args.output_dir

    fit_template = read_lines(fit_config)
    xsec_template = read_lines(xsec_config)

    if args.list:
        print("Available PDSP systematics:")
        for proc, param in available_table(xsec_template):
            print(f"  {proc:8s} {param}")
        return 0

    if output_dir.exists() and not output_dir.is_dir():
        raise ValueError(f"Output path exists and is not a directory: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    apps = resolve_apps()
    check_apps(apps)

    process_values = dict(args.process_sets)
    parameter_values = dict(args.parameter)

    cases: list[tuple[str, dict[str, float], dict[str, float]]] = []
    if args.nominal:
        cases.append(("nominal", {}, {}))

    if args.all_systematics is not None:
        cases.append((args.label or f"all_generator_{args.all_systematics:g}", {}, {
            name: args.all_systematics for _, name in available_table(xsec_template)
        }))

    if process_values or parameter_values:
        if args.label:
            cases.append((args.label, process_values, parameter_values))
        else:
            for proc, value in process_values.items():
                cases.append((f"{proc}_generator_{value:g}", {proc: value}, {}))
            for param, value in parameter_values.items():
                cases.append((f"{param}_generator_{value:g}", {}, {param: value}))

    if not cases:
        parser.error("Nothing to run. Use --set PROCESS=VALUE, --parameter NAME=VALUE, --all-systematics VALUE, or --nominal.")

    summary_rows = []
    manifests = []
    for label, proc_values, param_values in cases:
        manifest = run_case(
            label=label,
            fit_template=fit_template,
            xsec_template=xsec_template,
            output_dir=output_dir,
            apps=apps,
            diag_config=diag_config,
            workflow=args.workflow,
            process_values=proc_values,
            parameter_values=param_values,
            extra_args=args.extra_arg,
            dry_run=args.dry_run,
        )
        manifests.append(manifest)
        summary_rows.append({
            "label": manifest["label"],
            "status": manifest["status"],
            "returncode": manifest["returncode"],
            "fit_output": manifest["fit_output"],
            "posterior_predictive_output": manifest["posterior_predictive_output"],
            "prior_predictive_output": manifest["prior_predictive_output"],
            "process_values": json.dumps(manifest["process_values"], sort_keys=True),
            "parameter_values": json.dumps(manifest["parameter_values"], sort_keys=True),
        })
        print(f"{manifest['label']}: {manifest['status']} -> {manifest['fit_output']}")

    summary_csv = output_dir / "summary.csv"
    with summary_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0]))
        writer.writeheader()
        writer.writerows(summary_rows)

    with (output_dir / "manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifests, handle, indent=2, sort_keys=True)
        handle.write("\n")

    failed = [m for m in manifests if m["status"] == "failed"]
    return 1 if failed else 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(1)
