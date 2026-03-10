#!/usr/bin/env python3
import argparse
import json
import logging
from pathlib import Path
from statistics import median
from typing import Any, Dict, Iterable, List, Optional, Tuple

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Field groups (kept compatible with create_merged_json.py)
# ---------------------------------------------------------------------------
FIELD_GROUPS: Dict[str, List[str]] = {
    "Pre-filled fields": [
        "organism",
        "collecting_lab_sample_id",
        "enrichment_protocol",
        "enrichment_panel_version",
        "library_layout",
        "sequence_file_R1",
        "sequence_file_R2",
        "Sequence file R2",
        "sequencing_instrument_platform",
        "sequencing_instrument_model",
    ],
    "Laboratory fields": ["submitting_institution_id"],
    "De-hosting fields": [
        "dehosting_method_software_name",
        "dehosting_method_software_version",
        "per_reads_host",
    ],
    "Pre-processing fields": [
        "read_length",
        "preprocessing_software_name",
        "preprocessing_software_version",
        "preprocessing_params",
        "number_of_reads_sequenced",
        "pass_reads",
    ],
    "Mapping fields": [
        "reference_genome_accession",
        "mapping_software_name",
        "mapping_software_version",
        "mapping_params",
        "depth_of_coverage_threshold",
        "per_reads_virus",
    ],
    "Bioinformatics protocol fields": [
        "bioinformatics_protocol_software_name",
        "bioinformatics_protocol_software_version",
        "commercial_open_source_both",
        "bioinformatics_analysis_date",
    ],
    "Assembly fields": ["assembly", "assembly_version", "assembly_params"],
    "Variant calling fields": [
        "vcf_filename",
        "variant_calling_software_name",
        "variant_calling_software_version",
        "variant_calling_params",
        "number_of_variants_in_consensus",
        "number_of_variants_with_effect",
    ],
    "Consensus analysis fields": [
        "consensus_sequence_name",
        "consensus_sequence_filename",
        "consensus_sequence_md5",
        "consensus_sequence_software_name",
        "consensus_sequence_software_version",
        "consensus_params",
        "consensus_genome_length",
    ],
    "SARS-CoV-2 QC metrics": [
        "number_of_sgene_frameshifts",
        "per_ldmutations",
        "per_sgene_ambiguous",
        "per_sgene_coverage",
    ],
    "QC metrics fields": [
        "depth_of_coverage_value",
        "per_genome_greater_10x",
        "per_Ns",
        "number_of_Ns",
        "ns_per_100_kbp",
        "qc_test",
        "number_of_unambiguous_bases",
    ],
    "Lineage assignment fields": [
        "variant_name",
        "variant_designation",
        "lineage_assignment",
        "lineage_assignment_software_name",
        "lineage_assignment_software_version",
        "lineage_algorithm_software_version",
        "lineage_assignment_scorpio_version",
        "lineage_assignment_constellation_version",
        "lineage_assignment_date",
        "lineage_assignment_database_version",
    ],
    "Clade assignment fields": [
        "clade_assignment",
        "clade_assignment_software_name",
        "clade_assignment_software_version",
        "clade_assignment_software_database_version",
        "clade_assignment_date",
    ],
    "Type assignment fields": [
        "type_assignment",
        "type_assignment_software_name",
        "subtype_assignment_software_version",
        "type_assignment_software_database_version",
    ],
    "Subtype assignment fields": [
        "subtype_assignment",
        "subtype_assignment_software_name",
        "subtype_assignment_software_version",
        "subtype_assignment_software_database_version",
    ],
}

METADATA_METRICS_FIELDS = [
    "read_length",
    "reference_genome_accession",
    "number_of_reads_sequenced",
    "pass_reads",
    "per_reads_host",
    "per_reads_virus",
    "per_unmapped",
    "depth_of_coverage_value",
    "per_genome_greater_10x",
    "per_Ns",
    "number_of_Ns",
    "ns_per_100_kbp",
    "number_of_sgene_frameshifts",
    "number_of_unambiguous_bases",
    "per_ldmutations",
    "per_sgene_ambiguous",
    "per_sgene_coverage",
]

SOFTWARE_BENCHMARKING_FIELDS = [
    "bioinformatics_protocol_software_name",
    "bioinformatics_protocol_software_version",
    "commercial_open_source_both",
    "dehosting_method_software_name",
    "dehosting_method_software_version",
    "preprocessing_software_name",
    "preprocessing_software_version",
    "preprocessing_params",
    "mapping_software_name",
    "mapping_software_version",
    "mapping_params",
    "assembly",
    "assembly_version",
    "assembly_params",
    "variant_calling_software_name",
    "variant_calling_software_version",
    "variant_calling_params",
    "consensus_sequence_software_name",
    "consensus_sequence_software_version",
    "consensus_params",
    "depth_of_coverage_threshold",
    "clade_assignment_software_name",
    "clade_assignment_software_version",
    "clade_assignment_software_database_version",
    "lineage_assignment_software_name",
    "lineage_assignment_software_version",
    "lineage_algorithm_software_version",
    "lineage_assignment_scorpio_version",
    "lineage_assignment_constellation_version",
    "lineage_assignment_database_version",
    "type_assignment_software_name",
    "type_assignment_software_database_version",
    "subtype_assignment_software_name",
    "subtype_assignment_software_version",
    "subtype_assignment_software_database_version",
]

# ---------------------------------------------------------------------------
# Generic helpers
# ---------------------------------------------------------------------------
def load_json(path: str | Path) -> Any:
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def dump_json(data: Any, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2, ensure_ascii=False)
        handle.write("\n")


def is_meaningful(value: Any) -> bool:
    if value is None:
        return False
    if isinstance(value, str):
        stripped = value.strip()
        if stripped in {"", "000", "0000", "Not Provided [SNOMED:434941000124101]"}:
            return False
        return True
    return True


def normalize_label(value: Any) -> Optional[str]:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    low = text.lower()
    if low == "pass":
        return "Pass"
    if low == "fail":
        return "Fail"
    return text

def normalize(value: Any) -> Optional[str]:
    if value is None:
        return None
    text = str(value).strip()
    if text.lower() in {"", "na", "n/a", "null", "none"}:
        return None
    return text

def numeric_median(values: Iterable[Any]) -> Optional[float]:
    clean: List[float] = []
    for value in values:
        if value is None:
            continue
        try:
            clean.append(float(value))
        except (TypeError, ValueError):
            continue
    return median(clean) if clean else None


def count_existing(rows: Iterable[Dict[str, Any]], field_name: str) -> int:
    return sum(1 for row in rows if is_meaningful(row.get(field_name)))


def derive_type_subtype(row: Dict[str, Any]) -> Optional[str]:
    virus_type = row.get("type_assignment")
    subtype = row.get("subtype_assignment")
    if is_meaningful(virus_type) and is_meaningful(subtype):
        return f"{virus_type}/{subtype}"
    if is_meaningful(virus_type):
        return str(virus_type)
    if is_meaningful(subtype):
        return str(subtype)
    return None


def classification_outcome(lineage_match: Optional[bool], clade_match: Optional[bool]) -> Dict[str, Any]:
    n_matches = 0
    n_discrepancies = 0
    if lineage_match is not None:
        if lineage_match:
            n_matches += 1
        else:
            n_discrepancies += 1
    if clade_match is not None:
        if clade_match:
            n_matches += 1
        else:
            n_discrepancies += 1
    return {
        "lineage_match": lineage_match,
        "clade_match": clade_match,
        "number_matches": n_matches,
        "number_discrepancies": n_discrepancies,
    }

def applicable_group_fields(non_evaluable_groups: List[str]) -> Dict[str, List[str]]:
    groups: Dict[str, List[str]] = {}
    for group_name, group_fields in FIELD_GROUPS.items():
        # Pre-filled metadata are never evaluated for completeness
        if group_name == "Pre-filled fields":
            continue
        if group_name in non_evaluable_groups:
            continue
        groups[group_name] = group_fields
    return groups


def completeness_from_row(
    row: Optional[Dict[str, Any]],
    non_evaluable_groups: List[str],
) -> Tuple[int, int, List[str], List[str]]:
    groups = applicable_group_fields(non_evaluable_groups)

    if row is None:
        # No metadata row available: no group has been started, therefore nothing counts
        # towards expected fields or incompleteness.
        return 0, 0, [], []

    total_expected = 0
    filled_total = 0
    missing: List[str] = []
    partial_incomplete_groups: List[str] = []

    for group_name, group_fields in groups.items():
        group_filled_fields: List[str] = []
        group_missing_fields: List[str] = []

        for field in group_fields:
            if is_meaningful(row.get(field)):
                group_filled_fields.append(field)
            else:
                group_missing_fields.append(field)

        # If the whole block is empty, ignore it completely:
        # - it does not count toward total_expected_fields
        # - its empty fields are not counted as incomplete_fields
        # - it is not considered a primary incompleteness driver
        if len(group_filled_fields) == 0:
            continue

        # The block has been started, therefore all of its fields become expected.
        total_expected += len(group_fields)
        filled_total += len(group_filled_fields)

        # Only partially filled groups contribute missing fields and incompleteness drivers.
        if group_missing_fields:
            missing.extend(group_missing_fields)
            partial_incomplete_groups.append(group_name)

    return total_expected, filled_total, sorted(set(missing)), sorted(set(partial_incomplete_groups))


def extract_subset(row: Optional[Dict[str, Any]], fields: List[str]) -> Dict[str, Any]:
    if row is None:
        return {}
    return {field: row.get(field) for field in fields if is_meaningful(row.get(field))}


def unify_component_name(component_info: Dict[str, Any]) -> str:
    virus = component_info.get("virus", "")
    platform = component_info.get("sequencing_instrument_platform", "")
    if virus and platform:
        return f"{virus}, {platform}"
    return component_info.get("component_code", "Unknown")


def safe_int(value: Any) -> Optional[int]:
    if value is None:
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None

# ---------------------------------------------------------------------------
# Flexible comparison loaders
# ---------------------------------------------------------------------------
def index_comparison_payload(payload: Any) -> Dict[str, Dict[str, Any]]:
    if payload is None:
        return {}

    if isinstance(payload, list):
        out: Dict[str, Dict[str, Any]] = {}
        for entry in payload:
            sample_id = entry.get("sample_id") or entry.get("collecting_lab_sample_id")
            if sample_id:
                out[str(sample_id)] = entry
        return out

    if isinstance(payload, dict):
        if "components" in payload and isinstance(payload["components"], dict):
            out: Dict[str, Dict[str, Any]] = {}
            for comp in payload["components"].values():
                samples = comp.get("samples", {})
                if isinstance(samples, dict):
                    for sid, entry in samples.items():
                        out[str(sid)] = entry
                elif isinstance(samples, list):
                    for entry in samples:
                        sid = entry.get("sample_id") or entry.get("collecting_lab_sample_id")
                        if sid:
                            out[str(sid)] = entry
            return out
        if all(isinstance(v, dict) for v in payload.values()):
            return {str(k): v for k, v in payload.items()}
    return {}


def pick_sample_metrics(index: Dict[str, Dict[str, Any]], sample_id: str, comp_code: Optional[str] = None) -> Dict[str, Any]:
    if sample_id in index:
        return index[sample_id]
    if comp_code and comp_code in index and isinstance(index[comp_code], dict):
        comp_payload = index[comp_code]
        if sample_id in comp_payload:
            return comp_payload[sample_id]
        if "samples" in comp_payload and isinstance(comp_payload["samples"], dict):
            return comp_payload["samples"].get(sample_id, {})
    return {}


# ---------------------------------------------------------------------------
# Builders
# ---------------------------------------------------------------------------
def build_lab_json(
    expected_data: Dict[str, Any],
    metadata_rows: List[Dict[str, Any]],
    consensus_comparison: Optional[Any] = None,
    variant_comparison: Optional[Any] = None,
    lab_name_override: Optional[str] = None,
    figures_root: Optional[str] = None,
) -> Dict[str, Any]:
    metadata_by_sample = {row["collecting_lab_sample_id"]: row for row in metadata_rows}
    consensus_index = index_comparison_payload(consensus_comparison)
    variant_index = index_comparison_payload(variant_comparison)

    institution_ids = sorted({row.get("submitting_institution_id") for row in metadata_rows if row.get("submitting_institution_id")})
    if len(institution_ids) != 1:
        raise ValueError(f"Expected exactly one submitting_institution_id in metadata, found: {institution_ids}")
    lab_id = institution_ids[0]

    output: Dict[str, Any] = {
        "lab": {
            "submitting_institution_id": lab_id,
            "laboratory_name": lab_name_override or lab_id,
        },
        "metadata": {},
        "components": {},
    }

    all_expected_total = 0
    all_filled_total = 0
    global_incomplete_counter: Dict[str, int] = {}

    for comp_code, comp_expected in expected_data["components"].items():
        comp_rows: List[Dict[str, Any]] = []
        component_out: Dict[str, Any] = {
            "display_name": unify_component_name(comp_expected),
            "component_code": comp_code,
            "virus": comp_expected.get("virus"),
            "sequencing_instrument_platform": comp_expected.get("sequencing_instrument_platform"),
            "library_layout": comp_expected.get("library_layout"),
            "enrichment_protocol": comp_expected.get("enrichment_protocol"),
            "enrichment_panel": comp_expected.get("enrichment_panel"),
            "enrichment_panel_version": comp_expected.get("enrichment_panel_version"),
            "metadata": {},
            "samples": {},
        }

        sample_consensus_identities: List[float] = []
        sample_consensus_discrepancies: List[float] = []
        total_classification_matches = 0

        for sample_id, sample_expected in comp_expected["samples"].items():
            row = metadata_by_sample.get(sample_id)
            if row is None:
                continue

            comp_rows.append(row)

            non_evaluable_groups = sample_expected.get("non_evaluable_metadata", [])
            total_expected_fields, filled_fields, incomplete_fields, incomplete_groups = completeness_from_row(row, non_evaluable_groups)
            all_expected_total += total_expected_fields
            all_filled_total += filled_fields
            for group_name in incomplete_groups:
                global_incomplete_counter[group_name] = global_incomplete_counter.get(group_name, 0) + 1

            expected_lineage = normalize(sample_expected.get("expected_lineage"))
            expected_clade = normalize(sample_expected.get("expected_clade"))
            expected_qc = normalize_label(sample_expected.get("expected_qc"))

            if comp_expected.get("virus") == "SARS-CoV-2":
                reported_lineage = normalize(row.get("lineage_assignment") if row else None)
            else:
                reported_lineage = derive_type_subtype(row) if row else None
            reported_clade = normalize(row.get("clade_assignment") if row else None)
            lineage_match = None if expected_lineage is None else (reported_lineage == expected_lineage)
            clade_match = None if expected_clade is None else (reported_clade == expected_clade)
            classification_summary = classification_outcome(lineage_match, clade_match)
            total_classification_matches += classification_summary["number_matches"]

            reported_qc = normalize_label(row.get("qc_test"))
            qc_match = None if expected_qc is None else (reported_qc == expected_qc)

            consensus_metrics = pick_sample_metrics(consensus_index, sample_id, comp_code)
            variant_metrics = pick_sample_metrics(variant_index, sample_id, comp_code)

            consensus_breakdown = {
                "wrong_nt": consensus_metrics.get("wrong_nt") or consensus_metrics.get("substitutions") or 0,
                "ambiguity2nt": consensus_metrics.get("ambiguity2nt") or 0,
                "nt2ambiguity": consensus_metrics.get("nt2ambiguity") or consensus_metrics.get("excess_ambiguous") or 0,
                "ns2nt": consensus_metrics.get("ns2nt") or consensus_metrics.get("missing_Ns") or 0,
                "nt2ns": consensus_metrics.get("nt2ns") or consensus_metrics.get("excess_Ns") or 0,
                "insertions": consensus_metrics.get("insertions") or 0,
                "deletions": consensus_metrics.get("deletions") or 0,
            }

            consensus_total_discrepancies = sum(consensus_breakdown.values())

            consensus_block: Dict[str, Any] = {
                "genome_identity_pct": consensus_metrics.get("genome_identity_pct"),
                "total_discrepancies": consensus_total_discrepancies,
                "consensus_genome_length": row.get("consensus_genome_length"),
                "discrepancy_breakdown": consensus_breakdown,
            }
            if consensus_block.get("genome_identity_pct") is not None:
                sample_consensus_identities.append(float(consensus_block["genome_identity_pct"]))
            if consensus_block.get("total_discrepancies") is not None:
                sample_consensus_discrepancies.append(float(consensus_block["total_discrepancies"]))

            variant_wrong_nt = variant_metrics.get("wrong_nt") or 0
            variant_insertions = variant_metrics.get("insertions") or 0
            variant_deletions = variant_metrics.get("deletions") or 0

            variant_total_discrepancies = (
                variant_wrong_nt +
                variant_insertions +
                variant_deletions
            )

            number_of_variants_in_consensus = safe_int(row.get("number_of_variants_in_consensus"))
            number_of_variants_in_consensus_vcf = safe_int(variant_metrics.get("number_of_variants_in_consensus_vcf"))

            number_of_variants_with_effect = safe_int(row.get("number_of_variants_with_effect"))
            number_of_variants_with_effect_vcf = safe_int(variant_metrics.get("number_of_variants_with_effect_vcf"))

            if number_of_variants_in_consensus is not None and number_of_variants_in_consensus_vcf is not None:
                discrepancies_variants = number_of_variants_in_consensus_vcf - number_of_variants_in_consensus
            else:
                discrepancies_variants = None

            if number_of_variants_with_effect is not None and number_of_variants_with_effect_vcf is not None:
                discrepancies_variants_effect = number_of_variants_with_effect_vcf - number_of_variants_with_effect
            else:
                discrepancies_variants_effect = None

            variants_block: Dict[str, Any] = {
                "number_of_variants_in_consensus": number_of_variants_in_consensus,
                "number_of_variants_in_consensus_vcf": number_of_variants_in_consensus_vcf,
                "number_of_variants_with_effect": number_of_variants_with_effect,
                "discrepancies_in_reported_variants": discrepancies_variants,
            }

            if comp_expected.get("virus") == "SARS-CoV-2":
                variants_block.update({
                    "number_of_variants_with_effect_vcf": number_of_variants_with_effect_vcf,
                    "discrepancies_in_reported_variants_effect": discrepancies_variants_effect,
                })

            if comp_expected.get("virus") != "Influenza virus":
                variant_wrong_nt = variant_metrics.get("wrong_nt")
                variant_insertions = variant_metrics.get("insertions")
                variant_deletions = variant_metrics.get("deletions")

                variants_block.update({
                    "total_discrepancies": variant_total_discrepancies,
                    "wrong_nt": variant_wrong_nt,
                    "insertions": variant_insertions,
                    "deletions": variant_deletions,
                })

            classification_block = {
                "expected_lineage": expected_lineage,
                "lineage_assignment": reported_lineage,
                "expected_clade": expected_clade,
                "clade_assignment": reported_clade,
                **classification_summary,
                "variant_designation": row.get("variant_designation"),
                "variant_name": row.get("variant_name"),
            }

            sample_out: Dict[str, Any] = {
                "collecting_lab_sample_id": sample_id,
                "source": sample_expected.get("source"),
                "ref_sample": sample_expected.get("ref_sample"),
                "sequencing_instrument_platform": sample_expected.get("sequencing_instrument_platform"),
                "enrichment_protocol": sample_expected.get("enrichment_protocol"),
                "enrichment_panel": sample_expected.get("enrichment_panel"),
                "enrichment_panel_version": sample_expected.get("enrichment_panel_version"),
                "key_feature": sample_expected.get("key_feature"),
                "expected_qc": expected_qc,
                "total_expected_fields": total_expected_fields,
                "filled_fields": filled_fields,
                "incomplete_fields": incomplete_fields if incomplete_fields else [],
                "incomplete_groups": incomplete_groups if incomplete_groups else [],
                "consensus": consensus_block,
                "variants": variants_block,
                "classification": classification_block,
                "qc_test": reported_qc,
                "qc_match": qc_match,
                "metadata_metrics": extract_subset(row, METADATA_METRICS_FIELDS),
                "software_benchmarking": extract_subset(row, SOFTWARE_BENCHMARKING_FIELDS),
            }
            component_out["samples"][sample_id] = sample_out
        
        if not component_out["samples"]:
            continue

        component_total_expected = sum(v["total_expected_fields"] for v in component_out["samples"].values())
        component_filled = sum(v["filled_fields"] for v in component_out["samples"].values())
        component_completeness = (100.0 * component_filled / component_total_expected) if component_total_expected else None

        component_incomplete_counter: Dict[str, int] = {}
        for sample in component_out["samples"].values():
            for group_name in sample.get("incomplete_groups", []):
                component_incomplete_counter[group_name] = component_incomplete_counter.get(group_name, 0) + 1

        component_out["metadata"] = {
            "fasta_submitted": count_existing(comp_rows, "consensus_sequence_filename"),
            "fasta_expected": comp_expected.get("fasta_expected", len(comp_expected["samples"])),
            "vcf_submitted": count_existing(comp_rows, "vcf_filename"),
            "vcf_expected": comp_expected.get("vcf_expected", len(comp_expected["samples"])),
            "completeness_pct": component_completeness,
            "primary_incompleteness_drivers": [
                name for name, _ in sorted(component_incomplete_counter.items(), key=lambda kv: (-kv[1], kv[0]))[:10]
            ],
        }

        component_out["total_number_discrepancies"] = sum(v["consensus"].get("total_discrepancies") or 0 for v in component_out["samples"].values())
        component_out["median_genome_identity_pct"] = numeric_median(sample_consensus_identities)
        component_out["total_classification_matches"] = total_classification_matches

        output["components"][comp_code] = component_out

    global_completeness = (100.0 * all_filled_total / all_expected_total) if all_expected_total else None
    primary_incompleteness = sorted(global_incomplete_counter.items(), key=lambda kv: (-kv[1], kv[0]))
    output["metadata"] = {
        "completeness_pct": global_completeness,
        "primary_incompleteness_drivers": [name for name, _ in primary_incompleteness[:10]],
    }
    return output


# ---------------------------------------------------------------------------
# Input discovery helpers (legacy/iterative functionality preserved)
# ---------------------------------------------------------------------------
def find_lab_dirs(root_folder: Path) -> List[Path]:
    return [p for p in root_folder.iterdir() if p.is_dir() and (p / "RESULTS").is_dir()]


def find_validated_metadata(results_dir: Path) -> Optional[Path]:
    candidates = list(results_dir.glob("validated_metadata*.json")) + list(results_dir.glob("*validated_metadata*.json"))
    return candidates[0] if candidates else None


def find_optional_json(path_or_none: Optional[str]) -> Optional[Any]:
    if not path_or_none:
        return None
    p = Path(path_or_none)
    return load_json(p) if p.exists() else None


def infer_lab_name(metadata_rows: List[Dict[str, Any]], fallback: Optional[str] = None) -> str:
    return fallback or metadata_rows[0].get("submitting_institution_id") or "Unknown laboratory"


def validate_load_inputs(root_folder: Optional[str], expected_data_path: str, output_path: str, heading_file: Optional[str]) -> Tuple[Dict[str, Any], Optional[List[str]]]:
    if root_folder and not Path(root_folder).is_dir():
        raise ValueError(f"Given root_folder {root_folder} is not a valid folder")
    out = Path(output_path)
    if not out.exists():
        out.mkdir(parents=True, exist_ok=True)
    expected_data = load_json(expected_data_path)
    requested_fields = None
    if heading_file:
        requested_fields = list(load_json(heading_file).keys())
    return expected_data, requested_fields


# ---------------------------------------------------------------------------
# CLI workflows
# ---------------------------------------------------------------------------
def process_single_metadata_file(
    metadata_path: Path,
    expected_data: Dict[str, Any],
    output_path: Path,
    consensus_path: Optional[str] = None,
    variant_path: Optional[str] = None,
    lab_name: Optional[str] = None,
    figures_root: Optional[str] = None,
) -> Path:
    metadata_rows = load_json(metadata_path)
    consensus_data = find_optional_json(consensus_path)
    variant_data = find_optional_json(variant_path)
    out = build_lab_json(
        expected_data=expected_data,
        metadata_rows=metadata_rows,
        consensus_comparison=consensus_data,
        variant_comparison=variant_data,
        lab_name_override=lab_name,
        figures_root=figures_root,
    )
    dump_json(out, output_path)
    return output_path


def iterative_process(
    root_folder: str,
    expected_data_path: str,
    output_path: str,
    heading_file: Optional[str] = None,
    lab_cod: Optional[str] = None,
    lab_name: Optional[str] = None,
    consensus_comparison: Optional[str] = None,
    variant_comparison: Optional[str] = None,
) -> List[Path]:
    expected_data, _ = validate_load_inputs(root_folder, expected_data_path, output_path, heading_file)
    root = Path(root_folder)
    output_dir = Path(output_path)
    generated: List[Path] = []

    for folder in find_lab_dirs(root):
        results_dir = folder / "RESULTS"
        metadata_file = find_validated_metadata(results_dir)
        if not metadata_file:
            log.warning("Skipped invalid folder %s missing validated metadata JSON", folder.name)
            continue
        try:
            metadata_rows = load_json(metadata_file)
            inferred_lab_id = metadata_rows[0].get("submitting_institution_id") or folder.name
            if lab_cod and inferred_lab_id != lab_cod:
                continue
            final_lab_name = lab_name or inferred_lab_id
            figures_root = f"figures/labs/{inferred_lab_id}"
            output_file = output_dir / f"lab_{inferred_lab_id}.json"
            out = build_lab_json(
                expected_data=expected_data,
                metadata_rows=metadata_rows,
                consensus_comparison=find_optional_json(consensus_comparison),
                variant_comparison=find_optional_json(variant_comparison),
                lab_name_override=final_lab_name,
                figures_root=figures_root,
            )
            dump_json(out, output_file)
            generated.append(output_file)
        except Exception as exc:
            log.exception("Skipped %s due to processing error: %s", folder, exc)
    return generated


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Merged RELECOV per-lab JSON builder")

    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--metadata", help="Direct path to validated metadata JSON (single-lab mode)")
    mode.add_argument("-f", "--root_folder", help="Folder containing one subfolder per lab with RESULTS/validated_metadata.json")

    parser.add_argument("-e", "--expected_data", required=True, help="Path to expected_data.json")
    parser.add_argument("-o", "--output", required=True, help="Output file path (single mode) or output directory (folder mode)")
    parser.add_argument("--consensus-comparison", help="Optional path to consensus comparison JSON")
    parser.add_argument("--variant-comparison", help="Optional path to variant comparison JSON")
    parser.add_argument("--lab_cod", default=None, help="Restrict folder mode to a single lab code")
    parser.add_argument("--lab_name", default=None, help="Override full laboratory name")
    parser.add_argument("--heading_file", default=None, help="Optional heading/template JSON (kept for backward compatibility)")
    return parser.parse_args()


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args = parse_args()

    if args.metadata:
        expected_data = load_json(args.expected_data)
        output_path = Path(args.output)
        process_single_metadata_file(
            metadata_path=Path(args.metadata),
            expected_data=expected_data,
            output_path=output_path,
            consensus_path=args.consensus_comparison,
            variant_path=args.variant_comparison,
            lab_name=args.lab_name,
        )
        print(f"Generated {output_path}")
    else:
        generated = iterative_process(
            root_folder=args.root_folder,
            expected_data_path=args.expected_data,
            output_path=args.output,
            heading_file=args.heading_file,
            lab_cod=args.lab_cod,
            lab_name=args.lab_name,
            consensus_comparison=args.consensus_comparison,
            variant_comparison=args.variant_comparison,
        )
        print(f"Generated {len(generated)} file(s)")
        for path in generated:
            print(path)


if __name__ == "__main__":
    main()