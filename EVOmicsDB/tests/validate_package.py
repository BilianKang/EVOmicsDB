#!/usr/bin/env python3
"""Validate this distribution's structure and bounded regression checks.

Run from any directory: python3 tests/validate_package.py [--rscript /path/to/Rscript]
This does not install software, contact external services, rerun study datasets,
or establish complete manuscript or production-platform reproducibility.
"""

from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import subprocess
import sys
from typing import Callable


ROOT = Path(__file__).resolve().parents[1]
EXPECTED_IDS = set(range(1, 111))
EXPECTED_COUNTS = {
    "mRNA": (20614, 23, 20614, 0),
    "lncRNA": (13760, 23, 13760, 0),
    "miRNA": (4216, 25, 7076, 2860),
    "Proteome": (10258, 24, 10258, 0),
    "Phosphoproteome": (3341, 3, 3341, 0),
    "Metabolites": (2764, 3, 2789, 25),
}
COUNT_COLUMNS = (
    "distinct_identifiers",
    "matrix_entries",
    "case_sensitive_identifiers",
    "case_only_duplicates_removed",
)


def require(condition: bool, message: str) -> None:
    """Assertions must remain enabled even when Python runs with -O."""
    if not condition:
        raise ValueError(message)


def source_files(suffixes: set[str]) -> list[Path]:
    paths = sorted(
        path for path in ROOT.rglob("*")
        if path.is_file()
        and path.suffix in suffixes
        and not any(part in {".git", "__pycache__"} for part in path.relative_to(ROOT).parts)
    )
    require(bool(paths), f"No source files found for {sorted(suffixes)}")
    return paths


def python_syntax() -> dict:
    paths = source_files({".py"})
    for path in paths:
        ast.parse(path.read_bytes(), filename=str(path.relative_to(ROOT)))
    return {"files": len(paths), "method": "ast.parse; no source execution"}


def command(argv: list[str], timeout: int = 300) -> dict:
    env = os.environ.copy()
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    result = subprocess.run(
        argv, cwd=ROOT, env=env, capture_output=True, text=True,
        errors="replace", timeout=timeout, check=False,
    )
    output = (result.stdout + result.stderr).strip()
    require(result.returncode == 0, f"Exit {result.returncode}: {output[-8000:]}")
    return {"exit_code": result.returncode, "output": output[-8000:]}


def r_syntax(rscript: str) -> dict:
    paths = source_files({".r", ".R"})
    result = command([
        rscript, "--vanilla", "-e",
        "for (p in commandArgs(TRUE)) parse(file=p); cat('PASS: R syntax\\n')",
        *[str(path) for path in paths],
    ])
    return {"files": len(paths), **result}


def dataset_configuration() -> dict:
    version = (ROOT / "VERSION").read_text(encoding="utf-8").strip()
    require(version == "1.1.0", f"Expected release 1.1.0, found {version!r}")
    registry = json.loads((ROOT / "config/differential_analysis.json").read_text())
    rows = registry["datasets"]
    require(isinstance(rows, list), "Configuration datasets must be a list")
    ids = [int(row["dataset_id"]) for row in rows]
    require(len(ids) == 110, f"Expected 110 configurations, found {len(ids)}")
    require(len(set(ids)) == 110, "Configuration contains duplicate dataset IDs")
    require(set(ids) == EXPECTED_IDS, "Configuration IDs must be exactly 1 through 110")
    with (ROOT / "data/entries.tsv").open(newline="", encoding="utf-8-sig") as handle:
        entries = list(csv.DictReader(handle, delimiter="\t"))
    entry_ids = [int(row["dataset_id"]) for row in entries]
    require(len(entry_ids) == 110, f"Expected 110 entry records, found {len(entry_ids)}")
    require(len(set(entry_ids)) == 110, "Entry table contains duplicate dataset IDs")
    require(set(entry_ids) == set(ids), "Entry table and configuration dataset IDs differ")
    by_id = {int(row["dataset_id"]): row for row in entries}
    sample_columns = sum(int(row["Sample Size"]) for row in entries)
    require(sample_columns == 7120, "Registry sample-column total must equal 7,120")
    assets = json.loads((ROOT / "data/input_assets.json").read_text())["files"]
    require(isinstance(assets, list) and bool(assets), "External-input identities must be a nonempty list")
    assets_by_path = {}
    registered_inputs = {}
    for asset in assets:
        name = asset["path"]
        relative = PurePosixPath(name)
        require(not relative.is_absolute() and ".." not in relative.parts and str(relative) == name,
                f"Unsafe external-input path: {name}")
        require(re.fullmatch(r"[0-9a-fA-F]{64}", asset["sha256"]) is not None,
                f"Invalid external-input SHA256: {name}")
        require(type(asset["size"]) is int and asset["size"] >= 0,
                f"Invalid external-input size: {name}")
        if name in assets_by_path:
            previous = assets_by_path[name]
            require((asset["sha256"].lower(), asset["size"]) ==
                    (previous["sha256"].lower(), previous["size"]),
                    f"Conflicting identities for shared external input: {name}")
        assets_by_path[name] = asset
        if asset.get("category") == "registered":
            key = (int(asset["dataset_id"]), asset["field"])
            require(key not in registered_inputs, f"Duplicate registered input identity: {key}")
            registered_inputs[key] = name

    module_path = ROOT / "scripts/differential_config.py"
    spec = importlib.util.spec_from_file_location("evomics_validation_config", module_path)
    require(spec is not None and spec.loader is not None, "Cannot import configuration module")
    module = importlib.util.module_from_spec(spec)
    # Dataclasses resolve annotations using the module registration.
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    require(str(registry.get("schema_version")) == module.SCHEMA_VERSION,
            "Configuration schema differs from the runtime configuration parser")
    matching_fields = {
        "data_file": "Data File",
        "molecule_type": "Molecule Type",
        "omics_layer": "Omics Layer",
        "detection_method": "Detection Method",
        "quantification_method": "Quantification Method",
        "analysis_profile": "analysis_profile",
        "input_scale": "input_scale",
        "normalization_todo": "normalization_todo",
        "transform_todo": "transform_todo",
        "group_source": "group_source",
    }
    for row in rows:
        parsed = module.DifferentialDatasetConfig.from_dict(row)
        r_args = parsed.r_arguments()
        require(isinstance(r_args, dict) and bool(r_args), f"ID {parsed.dataset_id}: no R arguments")
        require(r_args["--input_scale"] == parsed.input_scale, "R argument input scale mismatch")
        require(r_args["--transform_todo"] == parsed.transform_todo, "R argument transform mismatch")
        entry = by_id[parsed.dataset_id]
        for config_field, entry_field in matching_fields.items():
            require(str(row.get(config_field) or "") == str(entry.get(entry_field) or ""),
                    f"ID {parsed.dataset_id}: registry/configuration mismatch for {config_field}")
        for field in ("data_file", "group_file", "probe_map_file"):
            name = row.get(field)
            if name:
                require(name in assets_by_path,
                        f"ID {parsed.dataset_id}: missing external-input identity for {field}")
        if row.get("group_file"):
            require(row["group_file"] == entry["group_file"],
                    f"ID {parsed.dataset_id}: group-file path differs")
        for field, entry_field in (("data_file", "Data File"), ("sample_metadata_file", "Metadata File")):
            require(registered_inputs.get((parsed.dataset_id, field)) == entry[entry_field],
                    f"ID {parsed.dataset_id}: registered input identity differs for {field}")
        if row["input_scale"] == "log2":
            require(row["transform_todo"] == "none",
                    f"ID {parsed.dataset_id}: log2 inputs must not be transformed again")
    id12 = next(row for row in rows if int(row["dataset_id"]) == 12)
    require(id12["input_scale"] == "log2", "ID 12 must have input_scale=log2")
    require(id12["transform_todo"] == "none", "ID 12 must have transform_todo=none")
    return {"configurations": len(rows), "entry_records": len(entries),
            "ids": "1–110, unique and matching", "id12": "log2 input; no additional transform",
            "sample_columns": sample_columns, "external_input_records": len(assets),
            "unique_external_input_paths": len(assets_by_path),
            "runtime_input_files_checked": False}


def figure1d_counts() -> dict:
    path = ROOT / "data/figures/Figure1d/counts.tsv"
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    require(len(rows) == 6, f"Expected six Figure 1d categories, found {len(rows)}")
    found = {}
    for row in rows:
        category = row["molecular_category"]
        require(category not in found, f"Duplicate Figure 1d category: {category}")
        found[category] = tuple(int(row[column]) for column in COUNT_COLUMNS)
    require(found == EXPECTED_COUNTS, f"Figure 1d counts differ from retained source: {found}")
    return {"categories": len(rows), "numeric_cells_checked": len(rows) * len(COUNT_COLUMNS),
            "scope": "Retained aggregate values; not recomputed from external matrices"}


def provenance() -> dict:
    document = json.loads((ROOT / "docs/provenance.json").read_text())
    rows = document["retained_files"]
    require(isinstance(rows, list) and bool(rows), "Provenance retained_files must be a nonempty list")
    required = {"source_archive", "source_path", "path", "original_sha256", "current_sha256", "change"}
    targets = set()
    for index, row in enumerate(rows, start=1):
        require(isinstance(row, dict), f"Provenance row {index} must be an object")
        require(required <= set(row), f"Provenance row {index} missing {sorted(required - set(row))}")
        for key in required:
            require(isinstance(row[key], str) and bool(row[key].strip()), f"Provenance row {index}: empty {key}")
        for key in ("original_sha256", "current_sha256"):
            require(re.fullmatch(r"[0-9a-fA-F]{64}", row[key]) is not None, f"Invalid {key} for {row['path']}")
        relative = PurePosixPath(row["path"])
        require(not relative.is_absolute() and ".." not in relative.parts,
                f"Unsafe mapped path: {row['path']}")
        require(str(relative) == row["path"], f"Noncanonical mapped path: {row['path']}")
        require(row["path"] not in targets, f"Duplicate provenance target: {row['path']}")
        targets.add(row["path"])
        path = (ROOT / relative).resolve()
        require(path.is_relative_to(ROOT), f"Mapped path leaves package: {row['path']}")
        require(path.is_file(), f"Missing retained source file: {row['path']}")
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        require(digest == row["current_sha256"].lower(), f"SHA256 mismatch: {row['path']}")
    return {"mapped_files": len(rows), "current_sha256": "verified for every mapping",
            "original_sha256": "format checked; original archives not required or reverified"}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rscript", default="Rscript", help="Rscript executable or absolute path")
    args = parser.parse_args()
    sys.dont_write_bytecode = True
    rscript = shutil.which(args.rscript)
    results = []

    def run_check(name: str, function: Callable[[], dict]) -> None:
        try:
            detail = function()
            results.append({"check": name, "status": "PASS", **detail})
        except Exception as exc:
            results.append({"check": name, "status": "FAIL", "error": f"{type(exc).__name__}: {exc}"})

    run_check("python_syntax", python_syntax)
    if rscript is None:
        results.append({"check": "r_executable", "status": "FAIL", "error": f"Rscript not found: {args.rscript}"})
    else:
        run_check("r_syntax", lambda: r_syntax(rscript))
        regressions = (
            ("identifier_namespace", "test_identifier_namespace.R", "gsea.r"),
            ("regulator_resolution", "test_regulator_resolution.R", "mirna_resolution_utils.R"),
            ("go_classification", "test_go_classification.R", "go_network.r"),
        )
        for name, test_file, source_file in regressions:
            run_check(name, lambda test_file=test_file, source_file=source_file: command([
                rscript, "--vanilla", str(ROOT / "tests" / test_file), str(ROOT / "R" / source_file),
            ]))
    run_check("external_kegg_restore", lambda: command([sys.executable, "-B", str(ROOT / "tests/test_external_kegg.py")]))
    run_check("dataset_configuration", dataset_configuration)
    run_check("figure1d_retained_counts", figure1d_counts)
    run_check("provenance_current_hashes", provenance)
    passed = all(result["status"] == "PASS" for result in results)
    report = {
        "status": "PASS" if passed else "FAIL",
        "scope": "Package syntax, three targeted R regressions, synthetic external-input restoration, configuration/entry consistency, retained Figure 1d aggregate values and mapped file integrity.",
        "not_assessed": ["Complete real-data analysis reruns", "All manuscript panels", "External input availability", "Production web-platform deployment", "Original archive hashes"],
        "checks": results,
    }
    print(json.dumps(report, indent=2, ensure_ascii=False))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
