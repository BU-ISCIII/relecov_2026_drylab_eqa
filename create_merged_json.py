import copy
import json
import logging
from collections import Counter

from pathlib import Path

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Schema: lab-specific fields overlaid on each expected sample entry.
# Values are None until filled from validated_metadata where available.
# ---------------------------------------------------------------------------

SAMPLES_BY_COMPONENT = {
    "SARS1": ["SARS1", "SARS2", "SARS3", "SARS4", "SARS5"],
    "SARS2": ["SARS6", "SARS7", "SARS8", "SARS9", "SARS10"],
    "FLU1": ["FLU1", "FLU2", "FLU3", "FLU4", "FLU5"],
    "FLU2": ["FLU6", "FLU7", "FLU8", "FLU9", "FLU10"]
}

FIELDS_CLASSIFICATION = {
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
        "sequencing_instrument_model"
    ],
    "Laboratory fields":  [
       "submitting_institution_id"
    ],
    "De-hosting fields": [
        "dehosting_method_software_name",
        "dehosting_method_software_version",
        "per_reads_host"
    ],
    "Pre-processing fields": [
        "read_length",
        "preprocessing_software_name",
        "preprocessing_software_version",
        "preprocessing_params",
        "number_of_reads_sequenced",
        "pass_reads"
    ],
    "Mapping fields": [
        "reference_genome_accession",
        "mapping_software_name",
        "mapping_software_version",
        "mapping_params",
        "depth_of_coverage_threshold",
        "per_reads_virus"
    ],
    "Bioinformatics protocol fields": [
        "bioinformatics_protocol_software_name",
        "bioinformatics_protocol_software_version",
        "commercial_open_source_both",
        "bioinformatics_analysis_date"
     ],
    "Assembly fields": [
        "assembly",
        "assembly_version",
        "assembly_params"
    ],
    "Variant calling fields": [
        "vcf_filename",
        "variant_calling_software_name",
        "variant_calling_software_version",
        "variant_calling_params",
        "number_of_variants_in_consensus",
        "number_of_variants_with_effect"
    ],
    "Consensus analysis fields": [
        "consensus_sequence_name",
        "consensus_sequence_filename",
        "consensus_sequence_md5",
        "consensus_sequence_software_name",
        "consensus_sequence_software_version",
        "consensus_params",
        "consensus_genome_length"
    ],
    "SARS-CoV-2 QC metrics": [
        "number_of_sgene_frameshifts",
        "per_ldmutations",
        "per_sgene_ambiguous",
        "per_sgene_coverage"
    ],
     "QC metrics fields": [
        "depth_of_coverage_value",
        "per_genome_greater_10x",
        "per_Ns",
        "number_of_Ns",
        "ns_per_100_kbp",
        "qc_test",
        "number_of_unambiguous_bases"
    ],
    "Lineage assignment fields" :[
        "variant_name",
        "variant_designation",
        "lineage_assignment",
        "lineage_assignment_software_name",
        "lineage_assignment_software_version",
        "lineage_algorithm_software_version",
        "lineage_assignment_scorpio_version",
        "lineage_assignment_constellation_version",
        "lineage_assignment_date",
        "lineage_assignment_file",
        "lineage_assignment_database_version"
    ],
    "Clade assignment fields": [
        "clade_assignment",
        "clade_assignment_software_name",
        "clade_assignment_software_version",
        "clade_assignment_software_database_version",
        "clade_assignment_date"
    ],
    "Type assignment fields": [
        "type_assignment",
        "type_assignment_software_name",
        "subtype_assignment_software_version",
        "type_assignment_software_database_version"
    ],
    "Subtype assignment fields": [
        "subtype_assignment",
        "subtype_assignment_software_name",
        "subtype_assignment_software_version",
        "subtype_assignment_software_database_version"
    ]
}


LAB_SAMPLE_FIELDS = {
    "total_expected_fields": None,
    "filled_fields": None,
    "incomplete_fields": None,
    "metadata_complete_pct": None,
    "reported_qc": None,
    "classification": {
        "expected_lineage": None,
        "reported_lineage": None,
        "expected_clade": None,
        "reported_clade": None,
        "concordance": None,
    },
    "consensus": None,
    "variants": None,
    "metadata_metrics": {
        "pct_genome_covered_ge_10x": None,
        "mean_depth": None,
        "pct_genome_masked_Ns": None,
        "pct_virus_reads": None,
        "pct_host_reads": None,
        "total_variants_reported": None,
        "variants_with_predicted_effect": None,
    },
    "metadata_metrics_network": None,
}

CONSENSUS_TEMPLATE = {
    "genome_identity_pct": None,
        "total_discrepancies": None,
        "indel_events": None,
        "ambiguous_bases_pct": None,
        "genome_completeness_pct": None,
        "discrepancy_breakdown": {
          "substitutions": None,
          "excess_Ns": None,
          "missing_Ns": None,
          "insertions": None,
          "deletions": None,
        }
}

VARIANTS_TEMPLATE = {
    "variants": {
        "sensitivity_pct": None,
        "precision_pct": None,
        "tp": None,
        "fp": None,
        "fn": None
    },
}

CLASSIFICATION_TEMPLATE = {
    "expected_lineage": None,
    "reported_lineage": None,
    "expected_clade": None,
    "reported_clade": None,
    "concordance": None
}

TEMPLATE_MAP = {
    "consensus": CONSENSUS_TEMPLATE,
    "variants": VARIANTS_TEMPLATE,
    "classification": CLASSIFICATION_TEMPLATE,
    #"metadata_metrics": METADATA_METRICS_TEMPLATE,

}

# Mapping: validated_metadata field → destination path inside sample entry
METADATA_SAMPLE_MAP = {
    "qc_test":                         ("reported_qc",),
    "per_genome_greater_10x":          ("metadata_metrics", "pct_genome_covered_ge_10x"),
    "depth_of_coverage_value":         ("metadata_metrics", "mean_depth"),
    "per_Ns":                          ("metadata_metrics", "pct_genome_masked_Ns"),
    "number_of_variants_in_consensus": ("metadata_metrics", "total_variants_reported"),
    "lineage_assignment":              ("classification", "reported_lineage"),
    "clade_assignment":                ("classification", "reported_clade"),
}

# Administrative fields in validated_metadata that are not bioinformatics metadata
_ADMIN_FIELDS = {"batch_id"}

# ---------------------------------------------------------------------------
# Iterative processing helpers (legacy / future use)
# ---------------------------------------------------------------------------

def process_sample_metadata(sample, requested_fields):
    """Fill the dictionary for each sample if possible"""
    sample_dict = copy.deepcopy(LAB_SAMPLE_FIELDS)
    num_total = len(requested_fields)
    match_fields = [x for x in sample.keys() if x in requested_fields]
    pct_comp = len(match_fields) / num_total
    # Fill the dictionary
    sample_dict["expected_fields"] = num_total
    sample_dict["filled_fields"] = match_fields
    sample_dict["metadata_complete_pct"] = pct_comp
    return sample_dict

def clean_meta_value(value):
    if isinstance(value, str):
        return value.split(" [")[0]
    else:
        return value

"""def autofill_matching_fields(data_to_fill, validated_metadata):
    #Fill all the fields with the same name in both dictionarieswith the values from validated_metadata
    def _rec_fill(data_to_fill, validated_metadata, parent=None):
        for k, v in data_to_fill.items():
            new_parent = ".".join(parent, k)
            if k in validated_metadata:
                if isinstance(v, dict):
                    _rec_fill(v, new_parent)
                    print("a")
        return 2
    filled_data = {}
    _rec_fill(data_to_fill, validated_metadata)
    return filled_data"""

def safe_match(sub, value):
    """Check if given substring is inside value, but only if value is string. Return False otherwise"""
    if isinstance(value, str):
        if sub in value:
            return True
        else:
            return False
    else:
        return False
    
def evaluate_group_completeness(sample_metadata, groups_to_evaluate):
    valid_groupdict = {
        k:v for k,v in FIELDS_CLASSIFICATION.items() if k in groups_to_evaluate
    }
    merged_dict = {
        group: {f: sample_metadata.get(f) for f in fields}
        for group, fields in valid_groupdict.items()
    }
    groups_completeness= {}
    fieldmatch_map = {
        "pct_params_missing": "_params",
        "pct_sw_version_missing": "software_version",
        "pct_db_version_missing": ("database_version", "constellation_version", "scorpio_version"),
    }
    for group, fieldict in merged_dict.items():
        fvalues = fieldict.values()
        missing_fields = [k for k,v in fieldict.items() if v is None]
        pct_fields_missing = (len(missing_fields) / len(fieldict)) * 100
        
        completness_groupdict = {"missing_fields": missing_fields, "pct_fields_missing": pct_fields_missing}
        for fieldname, regex in fieldmatch_map.items():
            if isinstance(regex, str):
                match_fields = [v for v in fvalues if safe_match(regex, v)]
            elif isinstance(regex, tuple):
                match_fields = [v for v in fvalues if any(safe_match(db, v) for db in regex)]
            
            if match_fields:
                completness_groupdict[fieldname] = sum([bool(v is None) for v in match_fields]) / len(match_fields)
            else:
                completness_groupdict[fieldname] = None
        groups_completeness[group] = completness_groupdict
    
    return groups_completeness

def get_sample_component(sample_id):
    comp_match = [k for k,v in SAMPLES_BY_COMPONENT.items() if sample_id in v]
    if not comp_match:
        raise ValueError(f"Sample name {sample_id} does not match with any defined sample")
    return comp_match[0]

def get_groups_by_component(sample_id, sample_component, expected_sample_data):
    try:
        invalid_groups = expected_sample_data["non_evaluable_metadata"]
    except KeyError:
        raise KeyError(
            f"Error getting `non_evaluable_metadata``in expected_data.json for sample {sample_id} with component {sample_component}"
        )
    return [k for k in FIELDS_CLASSIFICATION.keys() if k not in invalid_groups]

def get_sample_completeness(sample, groups_to_evaluate):
    sample_completeness = evaluate_group_completeness(sample, groups_to_evaluate)
    total_expected_fields = sum(
        [len(FIELDS_CLASSIFICATION[group]) for group in sample_completeness.keys()]
    )
    filled_fields = sum([len(v) for v in sample_completeness.values()])
    sample_completeness["total_expected_fields"] = total_expected_fields
    sample_completeness["filled_fields"] = filled_fields
    return sample_completeness

def get_lab_completeness(samples_completeness: dict) -> dict:
    """Aggregate per-sample completeness into lab-level percentages."""
    aggregated = {}

    for _, groups in samples_completeness.items():
        for group, metrics in groups.items():
            if group not in aggregated:
                aggregated[group] = {}
            for metric, value in metrics.items():
                if value is None:
                    continue
                aggregated[group].setdefault(metric, []).append(value)

    lab_completeness = {}
    for group, metrics in aggregated.items():
        lab_completeness[group] = {
            metric: sum(values) / len(values)
            for metric, values in metrics.items()
        }
        # Insert None if metric was always None or not found
        all_metrics = {
            metric
            for groups in samples_completeness.values()
            for metric in groups.get(group, {})
        }
        for metric in all_metrics - lab_completeness[group].keys():
            lab_completeness[group][metric] = None

    return lab_completeness

def lineage_concordance(expected_lineage, reported_lineage, expected_clade, reported_clade):
    """Compute lineage/clade concordance label."""
    if expected_lineage is None or reported_lineage is None:
        return None
    lineage_match = str(expected_lineage).strip().lower() == str(reported_lineage).strip().lower()
    clade_match = (
        expected_clade is None
        or (
            reported_clade is not None
            and str(expected_clade).strip().lower() == str(reported_clade).strip().lower()
        )
    )
    if lineage_match and clade_match:
        return "Exact"
    if lineage_match:
        return "Partial"
    return "Discordant"

def extract_sample_bioresults(data_dir, result_group):
    return


def evaluate_classification(sample, expected_sample_data, sample_component):
    if "SARS" in sample_component:
        groups_to_check = ("Lineage assignment fields", "Clade assignment fields")
    elif "FLU" in sample_component:
        groups_to_check = ("Type assignment fields", "Subtype assignment fields")
    else:
        raise ValueError(f"invalid sample_component {sample_component} for sample {sample.get('collecting_lab_sample_id')}")

    fields_to_extract = [field for group in groups_to_check for field in FIELDS_CLASSIFICATION[group]]
    reported_meta = {}
    for field in fields_to_extract:
        reported_meta[field] = sample.get(field)
    if "SARS" in sample_component:
        expected_lineage = expected_sample_data.get("expected_lineage")
        reported_lineage = sample.get("lineage_assignment")
        expected_clade = expected_sample_data.get("expected_clade")
        reported_clade = sample.get("clade_assignment")
    elif "FLU" in sample_component:
        expected_lineage = expected_sample_data.get("expected_lineage")
        reported_lineage = sample.get("type_assignment")
        expected_clade = expected_sample_data.get("expected_clade")
        reported_clade = sample.get("clade_assignment")
    classification_dict = {
        "expected_lineage": expected_lineage,
        "reported_lineage": reported_lineage,
        "lineage_match": bool(expected_lineage == reported_lineage),
        "expected_clade": expected_clade,
        "reported_clade": reported_clade,
        "clade_match": bool(expected_lineage == reported_lineage),
    }
    match_num = sum([classification_dict["lineage_match"], classification_dict["clade_match"]])
    classification_dict["number_matches"] = match_num
    classification_dict["number_discrepances"] = 2 - match_num
    return classification_dict

def process_lab_metadata(
        folder: Path,
        requested_fields: list,
        expected_data: dict,
    ):
    files_to_process = {}
    valid_meta_files = [x.as_posix() for x in folder.glob("**/RESULTS/validated_metadata.json")]
    if not valid_meta_files:
        log.error(f"No metadata files found in RESULTS for folder {folder.name}")
        raise ValueError(f"No metadata files found in RESULTS for folder {folder.name}. Skipped")
    metadata_file = valid_meta_files[0]
    with open(metadata_file, "r") as f:
        validated_meta = json.load(f) # Nota: This json is a list of dicts
    processed_metadata = {}
    """pre_filled_data = autofill_matching_fields(
        data_to_fill=expected_data,
        validated_metadata=validated_meta
    )"""
    sample_completeness = {}
    samples_results = {}
    total_fields = len(requested_fields)
    for sample in validated_meta:
        sample_id = sample.get("collecting_lab_sample_id")
        sample_component = get_sample_component(sample_id)
        expected_sample_data = expected_data["components"][sample_component][sample_id]
        
        groups_to_evaluate = get_groups_by_component(sample_id, sample_component, expected_sample_data)
        sample_completeness[sample_id] = get_sample_completeness(sample, groups_to_evaluate)
        if "classification" in groups_to_evaluate:
            samples_results["classification"] = evaluate_classification(
                sample, expected_sample_data, sample_component
            )

    lab_completeness = get_lab_completeness(sample_completeness)
    files_to_process["metadata"] = valid_meta_files[0]

    return lab_completeness

def iterative_process(
        root_folder: str,
        expected_data_path: str,
        output_path: str,
        heading_file: str,
        lab_cod: str,
        lab_name: str,
    ):
    expected_data, requested_fields = validate_load_inputs(
        root_folder=root_folder,
        expected_data_path=expected_data_path,
        output_path=output_path,
        heading_file=heading_file
    )

    lab_folder_list = [p for p in Path(root_folder).iterdir() if p.is_dir()]
    for folder in lab_folder_list:
        if (folder / "RESULTS").is_dir():
            try:
                lab_completeness = process_lab_metadata(
                    folder=folder,
                    requested_fields=requested_fields,
                    expected_data=expected_data
                )
            except Exception as e:
                print(f"[red]Skipped due to error processing {folder}: {e}")
                continue
        else:
            log.error(f"Skipped invalid folder {folder.name} missing RESULTS subfolder")

def validate_load_inputs(
    root_folder: str,
    expected_data_path: str,
    output_path: str,
    heading_file: str,
):
    if not Path(root_folder).is_dir():
        raise ValueError(f"Given root_folder {root_folder} is not a valid folder")
    if not Path(output_path).exists():
        raise ValueError(f"Given output_path {output_path} does not exist")
    with open(expected_data_path) as f:
        expected_data = json.load(f)
    if not expected_data:
        raise ValueError(f"Expected data file {expected_data_path} is empty")
    with open(heading_file) as f:
        requested_fields = json.load(f).keys()
    if not requested_fields:
        raise ValueError(f"Heading file {heading_file} is empty")

    

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Merge expected_data + validated_metadata into a per-lab JSON report"
    )
    parser.add_argument("-f", "--root_folder", required=True, help="Path to folder containing all laboratories")
    parser.add_argument("-e", "--expected_data", required=True, help="Path to expected_data.json")
    parser.add_argument("-o", "--output", required=True, help="Output JSON path")
    parser.add_argument("-h", "--heading_file", required=True, help="Json file containing the map of fields-columns asked in the excel template. e.g. `template_heading_eqa.json`")
    parser.add_argument("--lab_cod", default=None, help="Lab code (derived from metadata if omitted)")
    parser.add_argument("--lab_name", default=None, help="Full laboratory name")
    args = parser.parse_args()

    iterative_process(
        root_folder=args.root_folder,
        expected_data_path=args.expected_data,
        output_path=args.output,
        heading_file=args.heading_file,
        lab_cod=args.lab_cod,
        lab_name=args.lab_name,
    )

