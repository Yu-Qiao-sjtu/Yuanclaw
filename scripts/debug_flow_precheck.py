#!/usr/bin/env python3
"""Static precheck for YuanSeq business flow wiring.

Checks:
1) App sources required modules.
2) UI defines required input IDs for DE flow.
3) DE module references those IDs.
4) Optional dataset sanity summary for CSV counts matrix.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

REQUIRED_MODULES = [
    "modules/data_input.R",
    "modules/differential_analysis.R",
    "modules/kegg_enrichment.R",
    "modules/go_analysis.R",
    "modules/gsea_analysis.R",
    "modules/tf_activity.R",
    "modules/pathway_activity.R",
    "modules/chip_analysis.R",
    "modules/ai_interpretation.R",
]

REQUIRED_UI_INPUTS = [
    "species_select",
    "pval_cutoff",
    "log2fc_cutoff",
    "deg_pval_cutoff",
    "deg_log2fc_cutoff",
    "chip_pval_cutoff",
    "chip_log2fc_cutoff",
    "analyze",
]

REQUIRED_DE_INPUTS = [
    "control_group",
    "treat_group",
    "species_select",
    "pval_cutoff",
    "log2fc_cutoff",
]


def read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def check_app_module_wiring(root: Path) -> list[str]:
    issues = []
    app = read(root / "inst/shiny/app.R")
    for m in REQUIRED_MODULES:
        snippet = f'source("{m}")'
        if snippet not in app:
            issues.append(f"app.R missing module source: {snippet}")
    return issues


def check_ui_inputs(root: Path) -> list[str]:
    issues = []
    ui = read(root / "modules/ui_theme.R")
    for input_id in REQUIRED_UI_INPUTS:
        if f'"{input_id}"' not in ui:
            issues.append(f"ui_theme.R missing input id: {input_id}")
    return issues


def check_de_inputs(root: Path) -> list[str]:
    issues = []
    de = read(root / "modules/differential_analysis.R")
    for input_id in REQUIRED_DE_INPUTS:
        if f"input${input_id}" not in de:
            issues.append(f"differential_analysis.R missing input usage: input${input_id}")
    return issues


def dataset_summary(csv_path: Path) -> list[str]:
    out = []
    with csv_path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        rows = 0
        zero_all = 0
        for row in reader:
            rows += 1
            vals = [float(x) for x in row[1:]]
            if all(v == 0 for v in vals):
                zero_all += 1
    out.append(f"dataset: {csv_path}")
    out.append(f"columns: {len(header)} ({', '.join(header[:5])}{'...' if len(header) > 5 else ''})")
    out.append(f"rows: {rows}")
    out.append(f"all-zero rows: {zero_all}")
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description="YuanSeq business-flow static precheck")
    parser.add_argument("--repo", default=".", help="repo root path")
    parser.add_argument("--csv", default=None, help="optional count matrix csv for quick summary")
    args = parser.parse_args()

    root = Path(args.repo).resolve()

    issues = []
    issues.extend(check_app_module_wiring(root))
    issues.extend(check_ui_inputs(root))
    issues.extend(check_de_inputs(root))

    print("=== Business Flow Precheck ===")
    if issues:
        print("❌ Issues found:")
        for i in issues:
            print(f" - {i}")
    else:
        print("✅ Static wiring checks passed")

    csv_failed = False
    if args.csv:
        try:
            for line in dataset_summary(Path(args.csv).resolve()):
                print(line)
        except Exception as exc:
            csv_failed = True
            print(f"⚠️ dataset summary failed: {exc}")

    if issues:
        return 2
    if csv_failed:
        return 4
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
