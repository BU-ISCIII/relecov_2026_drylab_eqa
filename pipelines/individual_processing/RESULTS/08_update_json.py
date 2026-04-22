#!/usr/bin/env python3

import argparse
import json
from pathlib import Path


VARIANTS_DISCREPANCY_KEYS = (
    "wrong_nt",
    "insertions",
    "deletions",
    "missing",
    "denovo",
)

CONSENSUS_PREFIXES = ("SARS", "FLU")
VARIANTS_PREFIXES = ("SARS",)


def load_json(path: Path):
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def write_json(path: Path, payload):
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)


def normalize_consensus_breakdown(raw_breakdown: dict) -> dict:
    raw_breakdown = raw_breakdown or {}
    return {
        "wrong_nt": int(raw_breakdown.get("wrong_nt", 0)),
        "ambiguity2nt": int(raw_breakdown.get("ambiguity2nt", 0)),
        "nt2ambiguity": int(raw_breakdown.get("nt2ambigity", raw_breakdown.get("nt2ambiguity", 0))),
        "ns2nt": int(raw_breakdown.get("ns2nt", 0)),
        "nt2ns": int(raw_breakdown.get("nt2ns", 0)),
        "insertions": int(raw_breakdown.get("insertions", 0)),
        "deletions": int(raw_breakdown.get("deletions", 0)),
    }


def build_consensus_payload(consolidated_entry: dict, original_consensus: dict) -> dict:
    breakdown = normalize_consensus_breakdown(consolidated_entry.get("discrepancy_breakdown", {}))
    updated = dict(original_consensus)

    if "genome_identity_pct" in updated:
        updated["genome_identity_pct"] = consolidated_entry.get("genome_identity_pct")

    if "total_discrepancies" in updated:
        original_breakdown = original_consensus.get("discrepancy_breakdown", {})
        total_discrepancies = sum(
            int(breakdown.get(key, 0))
            for key in original_breakdown.keys()
        )
        updated["total_discrepancies"] = total_discrepancies

    if "discrepancy_breakdown" in updated and isinstance(updated["discrepancy_breakdown"], dict):
        updated_breakdown = dict(updated["discrepancy_breakdown"])
        for key in list(updated_breakdown.keys()):
            if key in breakdown:
                updated_breakdown[key] = breakdown[key]
        updated["discrepancy_breakdown"] = updated_breakdown

    return updated


def build_variants_payload(consolidated_entry: dict, original_variants: dict) -> dict:
    consolidated_entry = consolidated_entry or {}
    consolidated_values = {
        "successful_hits": consolidated_entry.get("successful_hits"),
        "high_and_low_freq": consolidated_entry.get("high_and_low_freq"),
        "high_freq_only": consolidated_entry.get("high_freq_only"),
        "low_freq_only": consolidated_entry.get("low_freq_only"),
        "number_of_variants_in_consensus_vcf": consolidated_entry.get("number_of_variants_in_consensus_vcf"),
        "number_of_variants_with_effect_vcf": consolidated_entry.get("number_of_variants_with_effect_vcf"),
        "wrong_nt": int(consolidated_entry.get("wrong_nt", 0)) if consolidated_entry.get("wrong_nt") is not None else None,
        "insertions": int(consolidated_entry.get("insertions", 0)) if consolidated_entry.get("insertions") is not None else None,
        "deletions": int(consolidated_entry.get("deletions", 0)) if consolidated_entry.get("deletions") is not None else None,
        "missing": int(consolidated_entry.get("missing", 0)) if consolidated_entry.get("missing") is not None else None,
        "denovo": int(consolidated_entry.get("denovo", 0)) if consolidated_entry.get("denovo") is not None else None,
    }

    discrepancy_values = [
        consolidated_values[key] for key in VARIANTS_DISCREPANCY_KEYS if isinstance(consolidated_values.get(key), int)
    ]
    consolidated_values["total_discrepancies"] = sum(discrepancy_values) if discrepancy_values else None

    updated = dict(original_variants)
    for key in list(updated.keys()):
        if key in consolidated_values:
            updated[key] = consolidated_values[key]

    return updated


def should_update_consensus(sample_code: str) -> bool:
    return sample_code.startswith(CONSENSUS_PREFIXES)


def should_update_variants(sample_code: str) -> bool:
    return sample_code.startswith(VARIANTS_PREFIXES)


def update_lab_json(lab_payload: dict, consolidated_payload: dict) -> tuple[int, int]:
    components = lab_payload.get("components", {})
    consensus_source = consolidated_payload.get("consensus", {})
    variants_source = consolidated_payload.get("variants", {})

    updated_consensus = 0
    updated_variants = 0

    for component in components.values():
        samples = component.get("samples", {})

        for sample_code, sample_entry in samples.items():
            if should_update_consensus(sample_code) and sample_code in consensus_source:
                consensus = sample_entry.get("consensus", {})
                sample_entry["consensus"] = build_consensus_payload(consensus_source[sample_code], consensus)
                updated_consensus += 1

            if should_update_variants(sample_code) and sample_code in variants_source:
                variants = sample_entry.get("variants", {})
                sample_entry["variants"] = build_variants_payload(variants_source[sample_code], variants)
                updated_variants += 1

    return updated_consensus, updated_variants


def parse_args():
    parser = argparse.ArgumentParser(
        description="Update lab JSON consensus/variants fields using consolidated_json_reports.json"
    )
    parser.add_argument("lab_json", help="Path to the original lab JSON to update in place")
    parser.add_argument(
        "consolidated_json",
        help="Path to consolidated_json_reports.json",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    lab_json_path = Path(args.lab_json)
    consolidated_json_path = Path(args.consolidated_json)

    lab_payload = load_json(lab_json_path)
    consolidated_payload = load_json(consolidated_json_path)

    updated_consensus, updated_variants = update_lab_json(lab_payload, consolidated_payload)
    write_json(lab_json_path, lab_payload)

    print(
        f"Updated {updated_consensus} consensus entries and {updated_variants} variants entries in {lab_json_path}"
    )


if __name__ == "__main__":
    main()
