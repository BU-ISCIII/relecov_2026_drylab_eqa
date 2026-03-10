#!/usr/bin/env python3
import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median
from typing import Any, Dict, Iterable, List, Optional

FIGURE_PATHS = {
    "consensus_summary": "figures/network/consensus_summary.png",
    "variant_summary": "figures/network/variant_summary.png",
    "sars_variant_reporting_summary": "figures/network/sars_variant_reporting_summary.png",
    "influenza_variant_reporting_summary": "figures/network/influenza_variant_reporting_summary.png",
    "influenza_reference_summary": "figures/network/influenza_reference_summary.png",
    "classification_summary": "figures/network/classification_summary.png",
    "metadata_completeness": "figures/network/metadata_completeness.png",
    "qc_match_rate_by_component": "figures/network/qc_match_rate_by_component.png",
}


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
        if stripped.lower() in {"", "na", "n/a", "null", "none"}:
            return False
        if stripped in {"000", "0000", "Not Provided [SNOMED:434941000124101]"}:
            return False
        return True
    return True


def safe_number(value: Any) -> Optional[float]:
    if not is_meaningful(value):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def safe_int(value: Any) -> Optional[int]:
    num = safe_number(value)
    if num is None:
        return None
    try:
        return int(num)
    except (TypeError, ValueError):
        return None


def pct(numerator: int, denominator: int) -> Optional[float]:
    if denominator == 0:
        return None
    return round(100.0 * numerator / denominator, 2)


def median_or_none(values: Iterable[Any]) -> Optional[float]:
    nums = [safe_number(v) for v in values]
    nums = [v for v in nums if v is not None]
    if not nums:
        return None
    return round(float(median(nums)), 4)


def min_or_none(values: Iterable[Any]) -> Optional[float]:
    nums = [safe_number(v) for v in values]
    nums = [v for v in nums if v is not None]
    if not nums:
        return None
    return round(float(min(nums)), 4)


def max_or_none(values: Iterable[Any]) -> Optional[float]:
    nums = [safe_number(v) for v in values]
    nums = [v for v in nums if v is not None]
    if not nums:
        return None
    return round(float(max(nums)), 4)


def get_workflow_signature(sample: Dict[str, Any]) -> Optional[str]:
    sb = sample.get("software_benchmarking", {})

    name = sb.get("bioinformatics_protocol_software_name")
    version = sb.get("bioinformatics_protocol_software_version")

    return software_signature(name, version)


def sample_is_evaluable(expected_sample: Dict[str, Any], group_name: str) -> bool:
    return group_name not in expected_sample.get("non_evaluable_metadata", [])


def load_lab_jsons(input_dir: Path) -> List[Dict[str, Any]]:
    files = sorted(input_dir.glob("*.json"))
    out = []
    for path in files:
        if path.name == "expected_data.json":
            continue
        out.append(load_json(path))
    return out


def software_key(name, version=None):
    if not is_meaningful(name):
        return None
    return (str(name).strip(), str(version).strip() if is_meaningful(version) else None)


def software_signature(name: Any, version: Any = None) -> Optional[tuple]:
    if not is_meaningful(name):
        return None
    clean_name = str(name).strip()
    clean_version = str(version).strip() if is_meaningful(version) else None
    return (clean_name, clean_version)


def most_common_or_none(values: Iterable[Any]) -> Any:
    clean = [v for v in values if is_meaningful(v)]
    if not clean:
        return None
    return Counter(clean).most_common(1)[0][0]


def collect_software_groups(
    participating_labs: List[Dict[str, Any]],
    comp_code: str,
    name_field: str,
    version_field: Optional[str] = None,
) -> Dict[tuple, List[Dict[str, Any]]]:
    groups = defaultdict(list)

    for lab in participating_labs:
        lab_id = lab["lab"]["submitting_institution_id"]
        comp = lab["components"][comp_code]

        for sample_id, sample in comp.get("samples", {}).items():
            sb = sample.get("software_benchmarking", {})
            key = software_signature(
                sb.get(name_field),
                sb.get(version_field) if version_field else None,
            )
            if key is None:
                continue

            groups[key].append({
                "lab_id": lab_id,
                "sample_id": sample_id,
                "sample": sample,
            })

    return groups


def build_software_entries(
    groups: Dict[tuple, List[Dict[str, Any]]],
    metrics_builder,
) -> List[Dict[str, Any]]:
    entries = []

    for (name, version), records in groups.items():
        entry = {
            "name": name,
            "n_labs": len({r["lab_id"] for r in records}),
        }
        if version is not None:
            entry["version"] = version

        entry.update(metrics_builder(records))
        entries.append(entry)

    return sorted(entries, key=lambda x: (x.get("name") or "", x.get("version") or ""))


def build_general(expected_data: Dict[str, Any], labs: List[Dict[str, Any]]) -> Dict[str, Any]:
    components_expected = expected_data["components"]
    total_invited = expected_data.get("general", {}).get("total_invited", expected_data.get("total_invited"))
    total_participants = len(labs)

    participation_per_component: Dict[str, int] = {}
    components_analysed_per_lab: List[int] = []
    total_fasta_submitted = 0
    total_fasta_expected = 0
    total_vcf_submitted = 0
    total_vcf_expected = 0

    lab_metadata_completeness: List[float] = []
    global_driver_counter: Counter = Counter()

    consensus_illumina_identity: List[float] = []
    consensus_nanopore_identity: List[float] = []

    sars_variant_discrepancies_illumina: List[float] = []
    sars_variant_discrepancies_nanopore: List[float] = []
    sars_variant_reporting_modes: Counter = Counter()

    influenza_variant_reporting_modes: Counter = Counter()
    influenza_reference_by_segment: Dict[str, set] = {seg: set() for seg in ["PB1", "PB2", "PA", "HA", "NP", "NA", "M", "NS"]}

    sars_lineage_matches = 0
    sars_lineage_total = 0
    sars_clade_matches = 0
    sars_clade_total = 0
    flu_type_matches = 0
    flu_type_total = 0
    flu_clade_matches = 0
    flu_clade_total = 0

    full_match_pct_per_component: List[float] = []

    qc_matches = 0
    qc_discrepancies = 0
    qc_total = 0

    software_names_count = 0
    software_names_total = 0
    software_version_count = 0
    software_version_total = 0
    coverage_threshold_count = 0
    coverage_threshold_total = 0
    frequency_threshold_count = 0
    frequency_threshold_total = 0
    reference_genome_count = 0
    reference_genome_total = 0

    metadata_incomplete_samples = 0
    metadata_evaluable_samples = 0

    all_workflows = set()
    consensus_softwares = set()
    variant_softwares = set()
    lineage_softwares = set()

    components_out: Dict[str, Any] = {}

    for lab in labs:
        lab_components = lab.get("components", {})
        components_analysed_per_lab.append(len(lab_components))

        meta_pct = safe_number(lab.get("metadata", {}).get("completeness_pct"))
        if meta_pct is not None:
            lab_metadata_completeness.append(meta_pct)

        for drv in lab.get("metadata", {}).get("primary_incompleteness_drivers", []):
            global_driver_counter[drv] += 1

        for comp_code, comp_expected in components_expected.items():
            if comp_code not in lab_components:
                continue

            comp = lab_components[comp_code]
            participation_per_component[comp_code] = participation_per_component.get(comp_code, 0) + 1

            total_fasta_submitted += min(
                safe_int(comp.get("metadata", {}).get("fasta_submitted")) or 0,
                safe_int(comp_expected.get("fasta_expected")) or 0,
            )
            total_fasta_expected += safe_int(comp_expected.get("fasta_expected")) or 0

            total_vcf_submitted += min(
                safe_int(comp.get("metadata", {}).get("vcf_submitted")) or 0,
                safe_int(comp_expected.get("vcf_expected")) or 0,
            )
            total_vcf_expected += safe_int(comp_expected.get("vcf_expected")) or 0

            for sample_id, sample in comp.get("samples", {}).items():
                expected_sample = comp_expected["samples"][sample_id]

                software_names_total += 1
                if is_meaningful(sample.get("software_benchmarking", {}).get("bioinformatics_protocol_software_name")):
                    software_names_count += 1

                software_version_total += 1
                if is_meaningful(sample.get("software_benchmarking", {}).get("bioinformatics_protocol_software_version")):
                    software_version_count += 1

                if sample_is_evaluable(expected_sample, "Mapping fields"):
                    coverage_threshold_total += 1
                    if is_meaningful(sample.get("software_benchmarking", {}).get("depth_of_coverage_threshold")):
                        coverage_threshold_count += 1

                    reference_genome_total += 1
                    if is_meaningful(sample.get("metadata_metrics", {}).get("reference_genome_accession")):
                        reference_genome_count += 1

                if sample_is_evaluable(expected_sample, "Variant calling fields"):
                    frequency_threshold_total += 1
                    if is_meaningful(sample.get("software_benchmarking", {}).get("variant_calling_params")):
                        frequency_threshold_count += 1

                total_expected_fields = safe_number(sample.get("total_expected_fields"))
                filled_fields = safe_number(sample.get("filled_fields"))
                if total_expected_fields is not None and total_expected_fields > 0:
                    metadata_evaluable_samples += 1
                    if filled_fields is None or filled_fields < total_expected_fields:
                        metadata_incomplete_samples += 1

                wf = get_workflow_signature(sample)
                if wf:
                    all_workflows.add(wf)

                sig = software_signature(
                    sample.get("software_benchmarking", {}).get("consensus_sequence_software_name"),
                    sample.get("software_benchmarking", {}).get("consensus_sequence_software_version"),
                )

                if sig:
                    consensus_softwares.add(sig)

                sig = software_signature(
                    sample.get("software_benchmarking", {}).get("variant_calling_software_name"),
                    sample.get("software_benchmarking", {}).get("variant_calling_software_version"),
                )

                if sig:
                    variant_softwares.add(sig)

                if comp_expected.get("virus") == "SARS-CoV-2":
                    sig = software_signature(
                        sample.get("software_benchmarking", {}).get("lineage_assignment_software_name"),
                        sample.get("software_benchmarking", {}).get("lineage_assignment_software_version"),
                    )

                    if sig:
                        lineage_softwares.add(sig)
                else:
                    sig = software_signature(
                        sample.get("software_benchmarking", {}).get("type_assignment_software_name"),
                        None,
                    )
                    if sig:
                        lineage_softwares.add(sig)

                    sig = software_signature(
                        sample.get("software_benchmarking", {}).get("subtype_assignment_software_name"),
                        sample.get("software_benchmarking", {}).get("subtype_assignment_software_version"),
                    )
                    if sig:
                        lineage_softwares.add(sig)

                qc_val = sample.get("qc_match")
                if qc_val is True:
                    qc_matches += 1
                    qc_total += 1
                elif qc_val is False:
                    qc_discrepancies += 1
                    qc_total += 1

                if sample_is_evaluable(expected_sample, "Consensus analysis fields"):
                    gi = safe_number(sample.get("consensus", {}).get("genome_identity_pct"))
                    if gi is not None:
                        if comp_expected.get("sequencing_instrument_platform") == "Illumina":
                            consensus_illumina_identity.append(gi)
                        else:
                            consensus_nanopore_identity.append(gi)

                cls = sample.get("classification", {})
                lineage_match = cls.get("lineage_match")
                clade_match = cls.get("clade_match")
                if comp_expected.get("virus") == "SARS-CoV-2":
                    if lineage_match is not None:
                        sars_lineage_total += 1
                        if lineage_match:
                            sars_lineage_matches += 1
                    if clade_match is not None:
                        sars_clade_total += 1
                        if clade_match:
                            sars_clade_matches += 1
                else:
                    if lineage_match is not None:
                        flu_type_total += 1
                        if lineage_match:
                            flu_type_matches += 1
                    if clade_match is not None:
                        flu_clade_total += 1
                        if clade_match:
                            flu_clade_matches += 1

                var = sample.get("variants", {})
                if comp_expected.get("virus") == "SARS-CoV-2" and sample_is_evaluable(expected_sample, "Variant calling fields"):
                    td = safe_number(var.get("total_discrepancies"))
                    if td is not None:
                        if comp_expected.get("sequencing_instrument_platform") == "Illumina":
                            sars_variant_discrepancies_illumina.append(td)
                        else:
                            sars_variant_discrepancies_nanopore.append(td)

                    if var.get("high_and_low_freq"):
                        sars_variant_reporting_modes["high_and_low_freq"] += 1
                    elif var.get("high_freq_only"):
                        sars_variant_reporting_modes["high_freq_only"] += 1
                    elif var.get("low_freq_only"):
                        sars_variant_reporting_modes["low_freq_only"] += 1

                if comp_expected.get("virus") == "Influenza virus":
                    if var.get("high_and_low_freq"):
                        influenza_variant_reporting_modes["high_and_low_freq"] += 1
                    elif var.get("high_freq_only"):
                        influenza_variant_reporting_modes["high_freq_only"] += 1
                    elif var.get("low_freq_only"):
                        influenza_variant_reporting_modes["low_freq_only"] += 1

                    ref = sample.get("metadata_metrics", {}).get("reference_genome_accession")
                    if is_meaningful(ref):
                        ref_str = str(ref)
                        for seg in influenza_reference_by_segment:
                            if seg in ref_str.upper():
                                influenza_reference_by_segment[seg].add(ref_str)

    for comp_code, comp_expected in components_expected.items():
        participating_labs = [lab for lab in labs if comp_code in lab.get("components", {})]
        if not participating_labs:
            continue

        comp_obj = {
            "name": f'{comp_expected.get("virus")}, {comp_expected.get("sequencing_instrument_platform")}',
            "total_labs": len(participating_labs),
            "total_fasta": sum(min(safe_int(lab["components"][comp_code].get("metadata", {}).get("fasta_submitted")) or 0,
                                   safe_int(comp_expected.get("fasta_expected")) or 0) for lab in participating_labs),
            "total_vcf": sum(min(safe_int(lab["components"][comp_code].get("metadata", {}).get("vcf_submitted")) or 0,
                                 safe_int(comp_expected.get("vcf_expected")) or 0) for lab in participating_labs),
            "metadata_completeness_median": median_or_none(
                [lab["components"][comp_code].get("metadata", {}).get("completeness_pct") for lab in participating_labs]
            ),
            "metadata_completeness_min_pct": min_or_none(
                [lab["components"][comp_code].get("metadata", {}).get("completeness_pct") for lab in participating_labs]
            ),
            "metadata_completeness_max_pct": max_or_none(
                [lab["components"][comp_code].get("metadata", {}).get("completeness_pct") for lab in participating_labs]
            ),
        }

        consensus_vals = []
        discrepancy_vals = []
        sample_entries = []
        disc_breakdown_all = defaultdict(list)
        breakdown_per_sample = defaultdict(lambda: defaultdict(list))

        for lab in participating_labs:
            comp = lab["components"][comp_code]
            for sample_id, sample in comp.get("samples", {}).items():
                expected_sample = comp_expected["samples"][sample_id]
                if not sample_is_evaluable(expected_sample, "Consensus analysis fields"):
                    continue
                cons = sample.get("consensus", {})
                gi = safe_number(cons.get("genome_identity_pct"))
                td = safe_number(cons.get("total_discrepancies"))
                if gi is not None:
                    consensus_vals.append(gi)
                if td is not None:
                    discrepancy_vals.append(td)
                bd = cons.get("discrepancy_breakdown", {})
                for key in ["wrong_nt", "ambiguity2nt", "nt2ambiguity", "ns2nt", "nt2ns", "insertions", "deletions"]:
                    val = safe_number(bd.get(key))
                    if val is not None:
                        disc_breakdown_all[key].append(val)
                        breakdown_per_sample[sample_id][key].append(val)

        for sample_id, expected_sample in comp_expected["samples"].items():
            if not sample_is_evaluable(expected_sample, "Consensus analysis fields"):
                continue
            gis = []
            tds = []
            sample_bd = breakdown_per_sample[sample_id]
            for lab in participating_labs:
                sample = lab["components"][comp_code]["samples"].get(sample_id)
                if not sample:
                    continue
                cons = sample.get("consensus", {})
                gi = safe_number(cons.get("genome_identity_pct"))
                td = safe_number(cons.get("total_discrepancies"))
                if gi is not None:
                    gis.append(gi)
                if td is not None:
                    tds.append(td)
            sample_entries.append({
                "collecting_lab_sample_id": sample_id,
                "median_identity_pct": median_or_none(gis),
                "median_discrepancies": median_or_none(tds),
                "min": min_or_none(tds),
                "max": max_or_none(tds),
                "wrong_nt": median_or_none(sample_bd.get("wrong_nt", [])),
                "ambiguity2nt": median_or_none(sample_bd.get("ambiguity2nt", [])),
                "nt2ambigity": median_or_none(sample_bd.get("nt2ambiguity", [])),
                "ns2nt": median_or_none(sample_bd.get("ns2nt", [])),
                "nt2ns": median_or_none(sample_bd.get("nt2ns", [])),
                "insertions": median_or_none(sample_bd.get("insertions", [])),
                "deletions": median_or_none(sample_bd.get("deletions", [])),
            })

        comp_obj["consensus"] = {
            "median_identity": median_or_none(consensus_vals),
            "median_discrepancies": median_or_none(discrepancy_vals),
            "total_median_discrepancies": round(sum(v for v in [median_or_none(discrepancy_vals)] if v is not None), 4) if discrepancy_vals else None,
            "min_discrepancies": min_or_none(discrepancy_vals),
            "max_discrepancies": max_or_none(discrepancy_vals),
            "median_identity_pct": median_or_none(consensus_vals),
            "identity_pct_min": min_or_none(consensus_vals),
            "identity_pct_max": max_or_none(consensus_vals),
            "samples": sample_entries,
            "discrepancy_breakdown": {
                key: {
                    "median": median_or_none(vals),
                    "min": min_or_none(vals),
                    "max": max_or_none(vals),
                } for key, vals in disc_breakdown_all.items()
            },
            "most_frequent_discrepancy_pattern": max(disc_breakdown_all.items(), key=lambda kv: len(kv[1]))[0] if disc_breakdown_all else None,
            "fig_discrepancies_boxplot_by_sample": f"figures/{comp_code}/consensus_discrepancies_boxplot_by_sample.png",
            "fig_discrepancies_stacked_by_sample": f"figures/{comp_code}/consensus_discrepancies_stacked_by_sample.png",
            "fig_discrepancy_type_boxplot": f"figures/{comp_code}/consensus_discrepancy_type_boxplot.png",
        }

        variant_discs = []
        variant_breakdown_all = defaultdict(list)
        variant_samples = []
        for sample_id, expected_sample in comp_expected["samples"].items():
            if not sample_is_evaluable(expected_sample, "Variant calling fields"):
                continue
            tds = []
            bd = defaultdict(list)
            for lab in participating_labs:
                sample = lab["components"][comp_code]["samples"].get(sample_id)
                if not sample:
                    continue
                var = sample.get("variants", {})
                td = safe_number(var.get("total_discrepancies"))
                if td is not None:
                    variant_discs.append(td)
                    tds.append(td)
                for key in ["wrong_nt", "insertions", "deletions"]:
                    v = safe_number(var.get(key))
                    if v is not None:
                        variant_breakdown_all[key].append(v)
                        bd[key].append(v)
            variant_samples.append({
                "collecting_lab_sample_id": sample_id,
                "median_discrepancies": median_or_none(tds),
                "min": min_or_none(tds),
                "max": max_or_none(tds),
                "wrong_nt": median_or_none(bd.get("wrong_nt", [])),
                "insertions": median_or_none(bd.get("insertions", [])),
                "deletions": median_or_none(bd.get("deletions", [])),
            })

        comp_obj["variant"] = {
            "median_discrepancies": median_or_none(variant_discs),
            "total_median_discrepancies": round(sum(v for v in [median_or_none(variant_discs)] if v is not None), 4) if variant_discs else None,
            "min_discrepancies": min_or_none(variant_discs),
            "max_discrepancies": max_or_none(variant_discs),
            "samples": variant_samples,
            "discrepancy_breakdown": {
                key: {
                    "median": median_or_none(vals),
                    "min": min_or_none(vals),
                    "max": max_or_none(vals),
                } for key, vals in variant_breakdown_all.items()
            },
            "fig_discrepancies_boxplot_by_sample": f"figures/{comp_code}/variant_discrepancies_boxplot_by_sample.png",
            "fig_discrepancy_type_boxplot": f"figures/{comp_code}/variant_discrepancy_type_boxplot.png",
        }

        classification_matches_per_lab = []
        per_sample_cls = []

        lineage_hits = 0
        lineage_total = 0
        clade_hits = 0
        clade_total = 0
        discordant_evaluations = 0
        classification_evaluable_total = 0

        # 1. Total matches per lab across all samples of the component
        for lab in participating_labs:
            lab_total_matches = 0

            for sample_id in comp_expected["samples"].keys():
                sample = lab["components"][comp_code]["samples"].get(sample_id)
                if not sample:
                    continue

                cls = sample.get("classification", {})
                nm = safe_int(cls.get("number_matches"))
                nd = safe_int(cls.get("number_discrepancies"))
                lm = cls.get("lineage_match")
                cm = cls.get("clade_match")

                if nm is not None:
                    lab_total_matches += nm

                if lm is not None:
                    lineage_total += 1
                    if lm is True:
                        lineage_hits += 1

                if cm is not None:
                    clade_total += 1
                    if cm is True:
                        clade_hits += 1

                if nd is not None and (lm is not None or cm is not None):
                    classification_evaluable_total += 1
                    if nd == 2:
                        discordant_evaluations += 1

            classification_matches_per_lab.append(lab_total_matches)

        # 2. Per-sample lineage/clade hit percentage
        for sample_id in comp_expected["samples"].keys():
            sample_lineage_hits = 0
            sample_lineage_total = 0
            sample_clade_hits = 0
            sample_clade_total = 0

            for lab in participating_labs:
                sample = lab["components"][comp_code]["samples"].get(sample_id)
                if not sample:
                    continue

                cls = sample.get("classification", {})
                lm = cls.get("lineage_match")
                cm = cls.get("clade_match")

                if lm is not None:
                    sample_lineage_total += 1
                    if lm is True:
                        sample_lineage_hits += 1

                if cm is not None:
                    sample_clade_total += 1
                    if cm is True:
                        sample_clade_hits += 1

            per_sample_cls.append({
                "collecting_lab_sample_id": sample_id,
                "lineage_hit_pct": pct(sample_lineage_hits, sample_lineage_total) if sample_lineage_total else None,
                "clade_hit_pct": pct(sample_clade_hits, sample_clade_total) if sample_clade_total else None,
            })

        comp_obj["typing"] = {
            "total_classification_matches_median": median_or_none(classification_matches_per_lab),
            "total_classification_matches_min": min_or_none(classification_matches_per_lab),
            "total_classification_matches_max": max_or_none(classification_matches_per_lab),
            "lineage_hit_pct": pct(lineage_hits, lineage_total),
            "clade_hit_pct": pct(clade_hits, clade_total),
            "discordance_pct": pct(discordant_evaluations, classification_evaluable_total),
            "samples": per_sample_cls,
            "fig_stacked_bar_by_sample": f"figures/{comp_code}/typing_outcome_stackedbar_by_sample.png",
        }


        qc_matches_comp = qc_disc_comp = qc_total_comp = 0
        qc_samples = []
        for sample_id, expected_sample in comp_expected["samples"].items():
            sample_qc_matches = sample_qc_disc = sample_qc_total = 0
            expected_qc = expected_sample.get("expected_qc")
            for lab in participating_labs:
                sample = lab["components"][comp_code]["samples"].get(sample_id)
                if not sample:
                    continue
                q = sample.get("qc_match")
                if q is True:
                    sample_qc_matches += 1
                    qc_matches_comp += 1
                    sample_qc_total += 1
                    qc_total_comp += 1
                elif q is False:
                    sample_qc_disc += 1
                    qc_disc_comp += 1
                    sample_qc_total += 1
                    qc_total_comp += 1
            qc_samples.append({
                "sample_id": sample_id,
                "gold_standard_qc": expected_qc,
                "match_rate_pct": pct(sample_qc_matches, sample_qc_total),
                "matches": sample_qc_matches,
                "discrepancies": sample_qc_disc,
                "total_evaluations": sample_qc_total,
            })

        comp_obj["qc"] = {
            "match_rate_pct": pct(qc_matches_comp, qc_total_comp),
            "matches": qc_matches_comp,
            "discrepancies": qc_disc_comp,
            "total_evaluations": qc_total_comp,
            "samples": qc_samples,
            "fig_qc_match_by_sample": f"figures/{comp_code}/qc_match_by_sample.png",
        }

        metric_map = {
            "per_genome_greater_10x": "per_genome_greater_10x",
            "depth_of_coverage_value": "depth_of_coverage_value",
            "per_Ns": "per_Ns",
            "per_reads_virus": "per_reads_virus",
            "per_reads_host": "per_reads_host",
            "number_of_variants_in_consensus": "number_of_variants_in_consensus",
        }
        if comp_expected.get("virus") == "SARS-CoV-2":
            metric_map["number_of_variants_with_effect"] = "number_of_variants_with_effect"

        meta_metric_samples = []
        for sample_id, expected_sample in comp_expected["samples"].items():
            entry = {"sample_id": sample_id}
            has_any = False
            for out_key, in_key in metric_map.items():
                vals = []
                for lab in participating_labs:
                    sample = lab["components"][comp_code]["samples"].get(sample_id)
                    if not sample:
                        continue
                    val = sample.get("metadata_metrics", {}).get(in_key)
                    if is_meaningful(val):
                        vals.append(val)
                if vals:
                    has_any = True
                    entry[out_key] = {
                        "median": median_or_none(vals),
                        "min": min_or_none(vals),
                        "max": max_or_none(vals),
                    }
            if has_any:
                meta_metric_samples.append(entry)
        comp_obj["metadata_metrics"] = {"samples": meta_metric_samples}

        benchmarking = {}

        # -------------------------------------------------
        # bioinformatics_protocol
        # -------------------------------------------------
        groups = collect_software_groups(
            participating_labs,
            comp_code,
            "bioinformatics_protocol_software_name",
            "bioinformatics_protocol_software_version",
        )

        def build_bioinfo_metrics(records):
            identities = []
            discrepancies = []
            metadata_completeness = []
            clade_hits = 0
            clade_total = 0
            lineage_hits = 0
            lineage_total = 0

            for r in records:
                sample = r["sample"]
                cons = sample.get("consensus", {})
                cls = sample.get("classification", {})

                gi = safe_number(cons.get("genome_identity_pct"))
                td = safe_number(cons.get("total_discrepancies"))
                ff = safe_number(sample.get("filled_fields"))
                te = safe_number(sample.get("total_expected_fields"))

                if gi is not None:
                    identities.append(gi)
                if td is not None:
                    discrepancies.append(td)
                if ff is not None and te is not None and te > 0:
                    metadata_completeness.append(100.0 * ff / te)

                cm = cls.get("clade_match")
                if cm is not None:
                    clade_total += 1
                    if cm is True:
                        clade_hits += 1

                lm = cls.get("lineage_match")
                if lm is not None:
                    lineage_total += 1
                    if lm is True:
                        lineage_hits += 1

            out = {
                "median_identity_pct": median_or_none(identities),
                "median_discrepancies": median_or_none(discrepancies),
                "median_metadata_completeness_pct": median_or_none(metadata_completeness),
                "clade_hit_pct": pct(clade_hits, clade_total),
                "lineage_hit_pct": pct(lineage_hits, lineage_total),
            }

            return out

        entries = build_software_entries(groups, build_bioinfo_metrics)
        if entries:
            benchmarking["bioinformatics_protocol"] = {
                "total_number": len(entries),
                "softwares": entries,
                "fig_metric_boxplots": f"figures/{comp_code}/bioinformatics_protocol_metric_boxplots_by_pipeline.png",
            }

        # -------------------------------------------------
        # dehosting
        # -------------------------------------------------
        groups = collect_software_groups(
            participating_labs,
            comp_code,
            "dehosting_method_software_name",
            "dehosting_method_software_version",
        )

        def build_dehosting_metrics(records):
            per_reads_host = []
            for r in records:
                val = safe_number(r["sample"].get("metadata_metrics", {}).get("per_reads_host"))
                if val is not None:
                    per_reads_host.append(val)

            return {
                "per_reads_host": median_or_none(per_reads_host),
            }

        entries = build_software_entries(groups, build_dehosting_metrics)
        if entries:
            benchmarking["dehosting"] = {
                "total_number": len(entries),
                "softwares": entries,
                "fig_metric_boxplots": f"figures/{comp_code}/dehosting_metric_boxplots_by_pipeline.png",
            }

        # -------------------------------------------------
        # preprocessing
        # -------------------------------------------------
        groups = collect_software_groups(
            participating_labs,
            comp_code,
            "preprocessing_software_name",
            "preprocessing_software_version",
        )

        def build_preprocessing_metrics(records):
            params = []
            reads_seq = []
            pass_reads = []

            for r in records:
                sample = r["sample"]
                params_val = sample.get("software_benchmarking", {}).get("preprocessing_params")
                if is_meaningful(params_val):
                    params.append(params_val)

                rs = safe_number(sample.get("metadata_metrics", {}).get("number_of_reads_sequenced"))
                pr = safe_number(sample.get("metadata_metrics", {}).get("pass_reads"))

                if rs is not None:
                    reads_seq.append(rs)
                if pr is not None:
                    pass_reads.append(pr)

            return {
                "params": most_common_or_none(params),
                "number_of_reads_sequenced": median_or_none(reads_seq),
                "pass_reads": median_or_none(pass_reads),
            }

        entries = build_software_entries(groups, build_preprocessing_metrics)
        if entries:
            benchmarking["preprocessing"] = {
                "total_number": len(entries),
                "softwares": entries,
                "fig_metric_boxplots": f"figures/{comp_code}/preprocessing_metric_boxplots_by_pipeline.png",
            }

        # -------------------------------------------------
        # mapping
        # -------------------------------------------------
        groups = collect_software_groups(
            participating_labs,
            comp_code,
            "mapping_software_name",
            "mapping_software_version",
        )

        def build_mapping_metrics(records):
            params = []
            coverage_thresholds = []
            virus_reads = []

            for r in records:
                sample = r["sample"]
                sb = sample.get("software_benchmarking", {})
                mm = sample.get("metadata_metrics", {})

                mp = sb.get("mapping_params")
                dct = sb.get("depth_of_coverage_threshold")

                if is_meaningful(mp):
                    params.append(mp)
                if is_meaningful(dct):
                    coverage_thresholds.append(dct)

                prv = safe_number(mm.get("per_reads_virus"))
                if prv is not None:
                    virus_reads.append(prv)

            return {
                "params": most_common_or_none(params),
                "depth_of_coverage_threshold": most_common_or_none(coverage_thresholds),
                "per_reads_virus": median_or_none(virus_reads),
            }

        entries = build_software_entries(groups, build_mapping_metrics)
        if entries:
            benchmarking["mapping"] = {
                "total_number": len(entries),
                "softwares": entries,
                "fig_metric_boxplots": f"figures/{comp_code}/mapping_metric_boxplots_by_pipeline.png",
            }

        # -------------------------------------------------
        # assembly
        # -------------------------------------------------
        groups = collect_software_groups(
            participating_labs,
            comp_code,
            "assembly",
            "assembly_version",
        )

        def build_assembly_metrics(records):
            params = []
            genome_lengths = []
            identities = []
            discrepancies = []

            for r in records:
                sample = r["sample"]
                sb = sample.get("software_benchmarking", {})
                cons = sample.get("consensus", {})

                ap = sb.get("assembly_params")
                if is_meaningful(ap):
                    params.append(ap)

                cgl = safe_number(cons.get("consensus_genome_length"))
                gi = safe_number(cons.get("genome_identity_pct"))
                td = safe_number(cons.get("total_discrepancies"))

                if cgl is not None:
                    genome_lengths.append(cgl)
                if gi is not None:
                    identities.append(gi)
                if td is not None:
                    discrepancies.append(td)

            return {
                "params": most_common_or_none(params),
                "consensus_genome_length": median_or_none(genome_lengths),
                "median_identity_pct": median_or_none(identities),
                "median_discrepancies": median_or_none(discrepancies),
            }

        entries = build_software_entries(groups, build_assembly_metrics)
        if entries:
            benchmarking["assembly"] = {
                "total_number": len(entries),
                "softwares": entries,
                "fig_metric_boxplots": f"figures/{comp_code}/assembly_metric_boxplots_by_pipeline.png",
            }

        # -------------------------------------------------
        # consensus_software
        # -------------------------------------------------
        groups = collect_software_groups(
            participating_labs,
            comp_code,
            "consensus_sequence_software_name",
            "consensus_sequence_software_version",
        )

        def build_consensus_metrics(records):
            params = []
            genome_lengths = []
            identities = []
            discrepancies = []

            for r in records:
                sample = r["sample"]
                sb = sample.get("software_benchmarking", {})
                cons = sample.get("consensus", {})

                cp = sb.get("consensus_params")
                if is_meaningful(cp):
                    params.append(cp)

                cgl = safe_number(cons.get("consensus_genome_length"))
                gi = safe_number(cons.get("genome_identity_pct"))
                td = safe_number(cons.get("total_discrepancies"))

                if cgl is not None:
                    genome_lengths.append(cgl)
                if gi is not None:
                    identities.append(gi)
                if td is not None:
                    discrepancies.append(td)

            return {
                "params": most_common_or_none(params),
                "consensus_genome_length": median_or_none(genome_lengths),
                "median_identity_pct": median_or_none(identities),
                "median_discrepancies": median_or_none(discrepancies),
            }

        entries = build_software_entries(groups, build_consensus_metrics)
        if entries:
            benchmarking["consensus_software"] = {
                "total_number": len(entries),
                "softwares": entries,
                "fig_metric_boxplots": f"figures/{comp_code}/consensus_metric_boxplots_by_pipeline.png",
            }

        # -------------------------------------------------
        # variant_calling
        # -------------------------------------------------
        groups = collect_software_groups(
            participating_labs,
            comp_code,
            "variant_calling_software_name",
            "variant_calling_software_version",
        )

        def build_variant_metrics(records):
            params = []
            n_variants = []
            n_effect = []
            discrepancies = []

            for r in records:
                sample = r["sample"]
                sb = sample.get("software_benchmarking", {})
                var = sample.get("variants", {})

                vp = sb.get("variant_calling_params")
                if is_meaningful(vp):
                    params.append(vp)

                nvic = safe_number(var.get("number_of_variants_in_consensus"))
                td = safe_number(var.get("total_discrepancies"))

                if nvic is not None:
                    n_variants.append(nvic)
                if td is not None:
                    discrepancies.append(td)

                if comp_expected.get("virus") == "SARS-CoV-2":
                    nvw = safe_number(var.get("number_of_variants_with_effect"))
                    if nvw is not None:
                        n_effect.append(nvw)

            out = {
                "params": most_common_or_none(params),
                "number_of_variants_in_consensus": median_or_none(n_variants),
                "median_discrepancies": median_or_none(discrepancies),
            }

            if comp_expected.get("virus") == "SARS-CoV-2":
                out["number_of_variants_with_effect"] = median_or_none(n_effect)

            return out

        entries = build_software_entries(groups, build_variant_metrics)
        if entries:
            benchmarking["variant_calling"] = {
                "total_number": len(entries),
                "softwares": entries,
                "fig_metric_boxplots": f"figures/{comp_code}/variant_calling_metric_boxplots_by_pipeline.png",
            }

        # -------------------------------------------------
        # clade_assignment
        # -------------------------------------------------
        groups = collect_software_groups(
            participating_labs,
            comp_code,
            "clade_assignment_software_name",
            "clade_assignment_software_version",
        )

        def build_clade_metrics(records):
            db_versions = []
            hits = 0
            total = 0

            for r in records:
                sample = r["sample"]
                sb = sample.get("software_benchmarking", {})
                cls = sample.get("classification", {})

                dbv = sb.get("clade_assignment_software_database_version")
                if is_meaningful(dbv):
                    db_versions.append(dbv)

                cm = cls.get("clade_match")
                if cm is not None:
                    total += 1
                    if cm is True:
                        hits += 1

            return {
                "database_version": most_common_or_none(db_versions),
                "clade_discordance_pct": pct(total - hits, total) if total else None,
                "clade_hit_pct": pct(hits, total),
            }

        entries = build_software_entries(groups, build_clade_metrics)
        if entries:
            benchmarking["clade_assignment"] = {
                "total_number": len(entries),
                "softwares": entries,
                "fig_metric_boxplots": f"figures/{comp_code}/clade_assignment_metric_boxplots_by_pipeline.png",
            }

        # -------------------------------------------------
        # lineage_assignment (only SARS)
        # -------------------------------------------------
        if comp_expected.get("virus") == "SARS-CoV-2":
            groups = collect_software_groups(
                participating_labs,
                comp_code,
                "lineage_assignment_software_name",
                "lineage_assignment_software_version",
            )

            def build_lineage_metrics(records):
                db_versions = []
                hits = 0
                total = 0

                for r in records:
                    sample = r["sample"]
                    sb = sample.get("software_benchmarking", {})
                    cls = sample.get("classification", {})

                    dbv = sb.get("lineage_assignment_database_version")
                    if is_meaningful(dbv):
                        db_versions.append(dbv)

                    lm = cls.get("lineage_match")
                    if lm is not None:
                        total += 1
                        if lm is True:
                            hits += 1

                return {
                    "database_version": most_common_or_none(db_versions),
                    "lineage_discordance_pct": pct(total - hits, total) if total else None,
                    "lineage_hit_pct": pct(hits, total),
                }

            entries = build_software_entries(groups, build_lineage_metrics)
            if entries:
                benchmarking["lineage_assignment"] = {
                    "total_number": len(entries),
                    "softwares": entries,
                    "fig_metric_boxplots": f"figures/{comp_code}/lineage_assignment_metric_boxplots_by_pipeline.png",
                }

        # -------------------------------------------------
        # type_assignment (only FLU)
        # -------------------------------------------------
        if comp_expected.get("virus") != "SARS-CoV-2":
            groups = collect_software_groups(
                participating_labs,
                comp_code,
                "type_assignment_software_name",
                None,
            )

            def build_type_metrics(records):
                db_versions = []
                hits = 0
                total = 0

                for r in records:
                    sample = r["sample"]
                    sb = sample.get("software_benchmarking", {})
                    cls = sample.get("classification", {})

                    dbv = sb.get("type_assignment_software_database_version")
                    if is_meaningful(dbv):
                        db_versions.append(dbv)

                    lm = cls.get("lineage_match")
                    if lm is not None:
                        total += 1
                        if lm is True:
                            hits += 1

                return {
                    "database_version": most_common_or_none(db_versions),
                    "type_discordance_pct": pct(total - hits, total) if total else None,
                    "type_hit_pct": pct(hits, total),
                }

            entries = build_software_entries(groups, build_type_metrics)
            if entries:
                benchmarking["type_assignment"] = {
                    "total_number": len(entries),
                    "softwares": entries,
                    "fig_metric_boxplots": f"figures/{comp_code}/type_assignment_metric_boxplots_by_pipeline.png",
                }

        # -------------------------------------------------
        # subtype_assignment (only FLU)
        # -------------------------------------------------
        if comp_expected.get("virus") != "SARS-CoV-2":
            groups = collect_software_groups(
                participating_labs,
                comp_code,
                "subtype_assignment_software_name",
                "subtype_assignment_software_version",
            )

            def build_subtype_metrics(records):
                db_versions = []
                hits = 0
                total = 0

                for r in records:
                    sample = r["sample"]
                    sb = sample.get("software_benchmarking", {})
                    cls = sample.get("classification", {})

                    dbv = sb.get("subtype_assignment_software_database_version")
                    if is_meaningful(dbv):
                        db_versions.append(dbv)

                    lm = cls.get("lineage_match")
                    if lm is not None:
                        total += 1
                        if lm is True:
                            hits += 1

                return {
                    "database_version": most_common_or_none(db_versions),
                    "subtype_discordance_pct": pct(total - hits, total) if total else None,
                    "subtype_hit_pct": pct(hits, total),
                }

            entries = build_software_entries(groups, build_subtype_metrics)
            if entries:
                benchmarking["subtype_assignment"] = {
                    "total_number": len(entries),
                    "softwares": entries,
                    "fig_metric_boxplots": f"figures/{comp_code}/subtype_assignment_metric_boxplots_by_pipeline.png",
                }

        comp_obj["benchmarking"] = benchmarking
        components_out[comp_code] = comp_obj

    out = {
        "eqa_name": expected_data.get("general", {}).get("eqa_name", "RELECOV 2026 Dry Lab EQA"),
        "total_invited": total_invited,
        "total_participants": total_participants,
        "total_participants_pct": pct(total_participants, total_invited) if total_invited else None,
        "participation_per_component": participation_per_component,
        "median_components_analysed_per_lab": median_or_none(components_analysed_per_lab),
        "total_fasta_submitted": total_fasta_submitted,
        "total_fasta_expected": total_fasta_expected,
        "total_vcf_submitted": total_vcf_submitted,
        "total_vcf_expected": total_vcf_expected,
        "submission_rates_pct": {
            "fasta": pct(total_fasta_submitted, total_fasta_expected),
            "vcf": pct(total_vcf_submitted, total_vcf_expected),
        },
        "figures": FIGURE_PATHS,
        "general_results": {
            "consensus": {
                "median_identity_illumina_pct": median_or_none(consensus_illumina_identity),
                "median_identity_nanopore_pct": median_or_none(consensus_nanopore_identity),
            },
            "sars_variants": {
                "median_discrepancy_illumina": median_or_none(sars_variant_discrepancies_illumina),
                "median_discrepancy_nanopore": median_or_none(sars_variant_discrepancies_nanopore),
                "high_and_low_freq_pct": pct(sars_variant_reporting_modes["high_and_low_freq"], sum(sars_variant_reporting_modes.values())),
                "low_freq_only_pct": pct(sars_variant_reporting_modes["low_freq_only"], sum(sars_variant_reporting_modes.values())),
                "high_freq_only_pct": pct(sars_variant_reporting_modes["high_freq_only"], sum(sars_variant_reporting_modes.values())),
            },
            "influenza_variants": {
                "high_and_low_freq_pct": pct(influenza_variant_reporting_modes["high_and_low_freq"], sum(influenza_variant_reporting_modes.values())),
                "low_freq_only_pct": pct(influenza_variant_reporting_modes["low_freq_only"], sum(influenza_variant_reporting_modes.values())),
                "high_freq_only_pct": pct(influenza_variant_reporting_modes["high_freq_only"], sum(influenza_variant_reporting_modes.values())),
                "total_distinct_references": len(set().union(*influenza_reference_by_segment.values())) if influenza_reference_by_segment else 0,
                "total_distinct_references_PB1": len(influenza_reference_by_segment["PB1"]),
                "total_distinct_references_PB2": len(influenza_reference_by_segment["PB2"]),
                "total_distinct_references_PA": len(influenza_reference_by_segment["PA"]),
                "total_distinct_references_HA": len(influenza_reference_by_segment["HA"]),
                "total_distinct_references_NP": len(influenza_reference_by_segment["NP"]),
                "total_distinct_references_NA": len(influenza_reference_by_segment["NA"]),
                "total_distinct_references_M": len(influenza_reference_by_segment["M"]),
                "total_distinct_references_NS": len(influenza_reference_by_segment["NS"]),
            },
            "classification": {
                "sars_cov_2_concordance_pct": pct(sars_lineage_matches, sars_lineage_total),
                "influenza_type_concordance_pct": pct(flu_type_matches, flu_type_total),
                "sars_clade_concordance_pct": pct(sars_clade_matches, sars_clade_total),
                "flu_clade_concordance_pct": pct(flu_clade_matches, flu_clade_total),
                "median_full_match_pct": median_or_none(full_match_pct_per_component),
                "min_full_match_pct": min_or_none(full_match_pct_per_component),
            },
        },
        "metadata_completeness": {
            "median_pct": median_or_none(lab_metadata_completeness),
            "min_pct": min_or_none(lab_metadata_completeness),
            "max_pct": max_or_none(lab_metadata_completeness),
            "software_names_pct": pct(software_names_count, software_names_total),
            "software_version_pct": pct(software_version_count, software_version_total),
            "coverage_threshold_pct": pct(coverage_threshold_count, coverage_threshold_total),
            "frequency_threshold_pct": pct(frequency_threshold_count, frequency_threshold_total),
            "reference_genome_pct": pct(reference_genome_count, reference_genome_total),
            "incomplete_parameters_pct": pct(metadata_incomplete_samples, metadata_evaluable_samples),
            "free_text_predefine_pct": None,
            "inconsistent_tool_version": None,
            "total_workflows": len(all_workflows),
            "total_consensus_softwares": len(consensus_softwares),
            "total_variant_softwares": len(variant_softwares),
            "total_lineage_softwares": len(lineage_softwares),
            "primary_incompleteness_drivers": [name for name, _ in global_driver_counter.most_common(10)],
        },
        "qc": {
            "match_rate_pct": pct(qc_matches, qc_total),
            "matches": qc_matches,
            "discrepancies": qc_discrepancies,
            "total_evaluations": qc_total,
        },
        "components": components_out,
    }
    return out


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build general.json from lab JSON files")
    parser.add_argument("--expected-data", required=True, help="Path to expected_data.json")
    parser.add_argument("--labs-dir", required=True, help="Directory containing per-lab JSON files")
    parser.add_argument("--output", required=True, help="Output general.json path")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    expected_data = load_json(args.expected_data)
    labs = load_lab_jsons(Path(args.labs_dir))
    general = build_general(expected_data, labs)
    dump_json(general, args.output)
    print(f"Generated {args.output}")


if __name__ == "__main__":
    main()
