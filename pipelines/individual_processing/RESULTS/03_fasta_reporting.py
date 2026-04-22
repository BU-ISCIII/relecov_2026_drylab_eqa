#!/usr/bin/env python3

import csv
import json
from pathlib import Path


BASE_DIR = Path.cwd()
REPORTS_DIR = BASE_DIR / "FASTA_ANALYSIS_REPORTS"
OUTPUT_JSON = BASE_DIR / "calculated_values.json"

TABLE_MAP = {
    "wrong_nt": ["Wrong nucleotide"],
    "ambiguity2nt": ["Nucleotide instead of ambiguity"],
    "nt2ambigity": ["Ambiguity instead of nucleotide"],
    "ns2nt": ["Nucleotide stretch instead of stretch of Ns"],
    "nt2ns": ["Stretch of Ns instead of nucleotide"],
    "insertions": ["Insertion relative to gold standard"],
    "deletions": ["Deletion relative to gold standard"],
}


def load_existing_output(path: Path):
    if not path.exists():
        return {}
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def sample_code_from_row(row):
    sample_id = row["Sample_ID"]
    return sample_id.split("_")[0]


def read_report_rows():
    rows = []
    csv_files = sorted(
        path for path in REPORTS_DIR.iterdir()
        if path.is_file() and path.name.endswith("_fasta_analysis.csv")
    )

    for csv_file in csv_files:
        with open(csv_file, "r", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            rows.extend(list(reader))
    return rows


def build_discrepancy_breakdown(rows):
    grouped = {}

    for row in rows:
        sample_code = sample_code_from_row(row)
        grouped.setdefault(sample_code, []).append(row)

    discrepancy_breakdown = {}
    for sample_code, sample_rows in grouped.items():
        discrepancy_breakdown[sample_code] = {"discrepancy_breakdown": {}}
        for json_key, valid_labels in TABLE_MAP.items():
            count = sum(1 for row in sample_rows if row["Resultado"] in valid_labels)
            discrepancy_breakdown[sample_code]["discrepancy_breakdown"][json_key] = count

    return discrepancy_breakdown


def main():
    rows = read_report_rows()
    discrepancy_breakdown = build_discrepancy_breakdown(rows)
    output_payload = load_existing_output(OUTPUT_JSON)
    output_payload.update(discrepancy_breakdown)

    with open(OUTPUT_JSON, "w", encoding="utf-8") as handle:
        json.dump(output_payload, handle, indent=4)


if __name__ == "__main__":
    main()
