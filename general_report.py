#!/usr/bin/env python3
import argparse
from html import parser
import json
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median
from typing import Any, Dict, Iterable, List, Optional

import matplotlib.pyplot as plt
from matplotlib.cbook import boxplot_stats

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

CBF_COLORS = {
    "match": "#0072B2",
    "discrepancy": "#E69F00",
    "high_freq_only": "#D55E00",
    "low_freq_only": "#56B4E9",
    "high_and_low_freq": "#009E73",
    "box_sars1": "#56B4E9",
    "box_sars2": "#0072B2",
    "box_flu1": "#E69F00",
    "box_flu2": "#CC79A7",
    "box_default": "#999999",
    "median": "#000000",
    "outlier": "#CC79A7",
}

COMPONENT_BOX_COLORS = {
    "SARS1": CBF_COLORS["box_sars1"],
    "SARS2": CBF_COLORS["box_sars2"],
    "FLU1": CBF_COLORS["box_flu1"],
    "FLU2": CBF_COLORS["box_flu2"],
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


def summarize_numeric_values(values: Iterable[Any], total_count: Optional[int] = None) -> Dict[str, Optional[float]]:
    raw_values = list(values)
    nums = [safe_number(v) for v in raw_values]
    observed_nums = [v for v in nums if v is not None]
    observed_count = len(observed_nums)
    denominator = total_count if total_count is not None else len(raw_values)
    zero_count = sum(1 for v in observed_nums if v == 0)

    return {
        "median": median_or_none(observed_nums),
        "min": min_or_none(observed_nums),
        "max": max_or_none(observed_nums),
        "observed_count": observed_count,
        "zero_count": zero_count,
        "missing_or_null_count": max(denominator - observed_count, 0),
    }


def dominant_metric_key(metric_summaries: Dict[str, Dict[str, Any]]) -> Optional[str]:
    best_key = None
    best_rank = None

    for key, summary in metric_summaries.items():
        median = safe_number(summary.get("median"))
        max_value = safe_number(summary.get("max"))
        if median is None:
            continue

        rank = (median, max_value if max_value is not None else float("-inf"))
        if best_rank is None or rank > best_rank:
            best_key = key
            best_rank = rank

    return best_key


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


def ensure_network_figures_dir(figures_dir: str | Path) -> Path:
    output_dir = Path(figures_dir) / "network"
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


def style_boxplot(bp: Dict[str, Any], labels: List[str]) -> None:
    for patch, label in zip(bp["boxes"], labels):
        patch.set_facecolor(COMPONENT_BOX_COLORS.get(label, CBF_COLORS["box_default"]))
        patch.set_edgecolor("#333333")
        patch.set_alpha(0.65)

    for median_line in bp["medians"]:
        median_line.set_color(CBF_COLORS["median"])
        median_line.set_linewidth(2)

    for whisker in bp["whiskers"]:
        whisker.set_color("#444444")
    for cap in bp["caps"]:
        cap.set_color("#444444")
    for flier in bp["fliers"]:
        flier.set_marker("o")
        flier.set_markerfacecolor("white")
        flier.set_markeredgecolor("#444444")
        flier.set_markersize(5)


def classification_hits_discrepancies_from_component(
    comp_data: Dict[str, Any],
    mode: str,
) -> tuple[int, int]:
    """
    mode:
      - 'lineage_type'
      - 'clade'
    """
    typing = comp_data.get("typing", {})
    samples = typing.get("samples", [])
    total_labs = comp_data.get("total_labs", 0)

    total_hits = 0
    total_discrepancies = 0

    for sample in samples:
        if mode == "lineage_type":
            hit_pct = sample.get("lineage_hit_pct")
        elif mode == "clade":
            hit_pct = sample.get("clade_hit_pct")
        else:
            raise ValueError(f"Unknown classification mode: {mode}")

        if hit_pct is None or total_labs == 0:
            continue

        sample_hits = round(total_labs * hit_pct / 100.0)
        sample_discrepancies = total_labs - sample_hits

        total_hits += sample_hits
        total_discrepancies += sample_discrepancies

    return total_hits, total_discrepancies


def make_stacked_classification_plot(
    general_data: Dict[str, Any],
    figures_dir: str | Path,
    mode: str,
    output_filename: str,
    title: str,
) -> str:
    components = general_data.get("components", {})

    component_names = []
    hit_counts = []
    discrepancy_counts = []

    for comp_code, comp_data in components.items():
        hits, discrepancies = classification_hits_discrepancies_from_component(comp_data, mode)

        component_names.append(comp_code)
        hit_counts.append(hits)
        discrepancy_counts.append(discrepancies)

    output_dir = ensure_network_figures_dir(figures_dir)
    output_path = output_dir / output_filename

    plt.figure(figsize=(10, 6))
    x_positions = list(range(len(component_names)))

    plt.bar(x_positions, hit_counts, label="Match", color=CBF_COLORS["match"])
    plt.bar(x_positions, discrepancy_counts, bottom=hit_counts, label="Discrepancy", color=CBF_COLORS["discrepancy"])

    plt.xticks(x_positions, component_names)
    plt.xlabel("Component")
    plt.ylabel("Total number of assignments")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()

    return str(output_path)


def qc_hits_discrepancies_from_component(comp_data: Dict[str, Any]) -> tuple[int, int]:
    """
    Calculate total QC matches and discrepancies across all samples of one component.
    """
    qc = comp_data.get("qc", {})
    samples = qc.get("samples", [])

    total_matches = 0
    total_discrepancies = 0

    for sample in samples:
        total_matches += safe_int(sample.get("matches")) or 0
        total_discrepancies += safe_int(sample.get("discrepancies")) or 0

    return total_matches, total_discrepancies


def collect_consensus_discrepancies_by_component(labs: List[Dict[str, Any]]) -> Dict[str, List[float]]:
    """
    Collect sample-level consensus discrepancy counts grouped by component.
    Only observations with a non-null genome identity are considered evaluable.
    """
    discrepancies_by_component: Dict[str, List[float]] = defaultdict(list)

    for lab in labs:
        for comp_code, comp in lab.get("components", {}).items():
            for sample in comp.get("samples", {}).values():
                consensus = sample.get("consensus", {})
                genome_identity_pct = safe_number(consensus.get("genome_identity_pct"))
                total_discrepancies = safe_number(consensus.get("total_discrepancies"))

                if genome_identity_pct is None or total_discrepancies is None:
                    continue

                discrepancies_by_component[comp_code].append(total_discrepancies)

    return discrepancies_by_component


def make_consensus_summary_plot(
    labs: List[Dict[str, Any]],
    figures_dir: str | Path,
    output_filename: str = "consensus_summary.png",
    title: str = "Consensus genome discrepancies relative to the gold standard",
) -> str:
    output_dir = ensure_network_figures_dir(figures_dir)
    output_path = output_dir / output_filename

    discrepancies_by_component = collect_consensus_discrepancies_by_component(labs)
    component_order = ["SARS1", "SARS2", "FLU1", "FLU2"]
    component_names = [comp for comp in component_order if discrepancies_by_component.get(comp)]
    data = [list(discrepancies_by_component[comp]) for comp in component_names]

    if not data:
        return str(output_path)

    fixed_y_upper = 500.0
    plotted_data = [list(vals) for vals in data]
    outlier_annotations = []

    for idx, values in enumerate(data):
        outliers_above_limit = sorted([value for value in values if value > fixed_y_upper], reverse=True)
        if not outliers_above_limit:
            continue

        plotted_data[idx] = [value for value in values if value <= fixed_y_upper]
        outlier_annotations.append((idx + 1, outliers_above_limit))

    plt.figure(figsize=(10, 6))
    bp = plt.boxplot(plotted_data, labels=component_names, showfliers=True, patch_artist=True)
    style_boxplot(bp, component_names)
    plt.xlabel("Component")
    plt.ylabel("Consensus discrepancies")
    plt.title(title)
    plt.ylim(0, fixed_y_upper)
    plt.xlim(0.5, len(component_names) + 0.5)

    for x_pos, outliers in outlier_annotations:
        display_value = outliers[0]
        y_marker = fixed_y_upper * 0.95
        y_text = fixed_y_upper * 0.91
        component_label = component_names[x_pos - 1]
        outlier_color = COMPONENT_BOX_COLORS.get(component_label, CBF_COLORS["outlier"])
        plt.text(
            x_pos,
            y_marker,
            "*",
            ha="center",
            va="center",
            fontsize=18,
            color=outlier_color,
            fontweight="bold",
        )
        plt.text(
            x_pos,
            y_text,
            f"Outlier: {display_value:g}",
            ha="center",
            va="top",
            fontsize=9,
            color=outlier_color,
            fontweight="bold",
        )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()

    return str(output_path)


def collect_sars_variant_discrepancies_by_component(labs: List[Dict[str, Any]]) -> Dict[str, List[float]]:
    """
    Collect sample-level variant discrepancy counts for SARS-CoV-2 components.
    """
    discrepancies_by_component: Dict[str, List[float]] = defaultdict(list)

    for lab in labs:
        for comp_code, comp in lab.get("components", {}).items():
            if comp_code not in {"SARS1", "SARS2"}:
                continue

            for sample in comp.get("samples", {}).values():
                variants = sample.get("variants", {})
                total_discrepancies = safe_number(variants.get("total_discrepancies"))
                if total_discrepancies is None:
                    continue

                discrepancies_by_component[comp_code].append(total_discrepancies)

    return discrepancies_by_component


def make_variant_summary_plot(
    labs: List[Dict[str, Any]],
    figures_dir: str | Path,
    output_filename: str = "variant_summary.png",
    title: str = "SARS-CoV-2 variant discrepancies relative to the reference set",
) -> str:
    output_dir = ensure_network_figures_dir(figures_dir)
    output_path = output_dir / output_filename

    discrepancies_by_component = collect_sars_variant_discrepancies_by_component(labs)
    component_order = ["SARS1", "SARS2"]
    component_names = [comp for comp in component_order if discrepancies_by_component.get(comp)]
    data = [list(discrepancies_by_component[comp]) for comp in component_names]

    if not data:
        return str(output_path)

    fixed_y_upper = 210.0
    plotted_data = [list(vals) for vals in data]
    outlier_annotations = []

    for idx, values in enumerate(data):
        outliers_above_limit = sorted([value for value in values if value > fixed_y_upper], reverse=True)
        if not outliers_above_limit:
            continue

        plotted_data[idx] = [value for value in values if value <= fixed_y_upper]
        outlier_annotations.append((idx + 1, outliers_above_limit))

    plt.figure(figsize=(8, 6))
    bp = plt.boxplot(plotted_data, labels=component_names, showfliers=True, patch_artist=True)
    style_boxplot(bp, component_names)
    plt.xlabel("Component")
    plt.ylabel("Variant discrepancies")
    plt.title(title)
    plt.ylim(0, fixed_y_upper)
    plt.xlim(0.5, len(component_names) + 0.5)

    for x_pos, outliers in outlier_annotations:
        display_value = outliers[0]
        y_marker = fixed_y_upper * 0.95
        y_text = fixed_y_upper * 0.91
        component_label = component_names[x_pos - 1]
        outlier_color = COMPONENT_BOX_COLORS.get(component_label, CBF_COLORS["outlier"])
        plt.text(
            x_pos,
            y_marker,
            "*",
            ha="center",
            va="center",
            fontsize=18,
            color=outlier_color,
            fontweight="bold",
        )
        plt.text(
            x_pos,
            y_text,
            f"Outlier: {display_value:g}",
            ha="center",
            va="top",
            fontsize=9,
            color=outlier_color,
            fontweight="bold",
        )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()

    return str(output_path)


def make_sars_variant_reporting_summary_plot(
    general_data: Dict[str, Any],
    figures_dir: str | Path,
    output_filename: str = "sars_variant_reporting_summary.png",
    title: str = "SARS-CoV-2 variant reporting practices across the network",
) -> str:
    output_dir = ensure_network_figures_dir(figures_dir)
    output_path = output_dir / output_filename

    sars_variants = general_data.get("general_results", {}).get("sars_variants", {})
    categories = [
        ("High + low freq", safe_number(sars_variants.get("high_and_low_freq_pct"))),
        ("Low freq only", safe_number(sars_variants.get("low_freq_only_pct"))),
        ("High freq only", safe_number(sars_variants.get("high_freq_only_pct"))),
    ]

    labels = [label for label, value in categories if value is not None]
    values = [value for _, value in categories if value is not None]

    if not values:
        return str(output_path)

    plt.figure(figsize=(8, 6))
    bars = plt.bar(
        labels,
        values,
        color=[
            CBF_COLORS["high_and_low_freq"],
            CBF_COLORS["low_freq_only"],
            CBF_COLORS["high_freq_only"],
        ][: len(values)],
    )
    plt.xlabel("Reporting mode")
    plt.ylabel("Laboratories (%)")
    plt.title(title)
    plt.ylim(0, 100)

    for bar, value in zip(bars, values):
        plt.text(
            bar.get_x() + bar.get_width() / 2,
            value + 1.5,
            f"{value:.2f}%",
            ha="center",
            va="bottom",
            fontsize=9,
        )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()

    return str(output_path)


def make_influenza_variant_reporting_summary_plot(
    general_data: Dict[str, Any],
    figures_dir: str | Path,
    output_filename: str = "influenza_variant_reporting_summary.png",
    title: str = "Influenza variant reporting practices across the network",
) -> str:
    output_dir = ensure_network_figures_dir(figures_dir)
    output_path = output_dir / output_filename

    influenza_variants = general_data.get("general_results", {}).get("influenza_variants", {})
    categories = [
        ("High + low freq", safe_number(influenza_variants.get("high_and_low_freq_pct"))),
        ("Low freq only", safe_number(influenza_variants.get("low_freq_only_pct"))),
        ("High freq only", safe_number(influenza_variants.get("high_freq_only_pct"))),
    ]

    labels = [label for label, value in categories if value is not None]
    values = [value for _, value in categories if value is not None]

    if not values:
        return str(output_path)

    plt.figure(figsize=(8, 6))
    bars = plt.bar(
        labels,
        values,
        color=[
            CBF_COLORS["high_and_low_freq"],
            CBF_COLORS["low_freq_only"],
            CBF_COLORS["high_freq_only"],
        ][: len(values)],
    )
    plt.xlabel("Reporting mode")
    plt.ylabel("Laboratories (%)")
    plt.title(title)
    plt.ylim(0, 100)

    for bar, value in zip(bars, values):
        plt.text(
            bar.get_x() + bar.get_width() / 2,
            value + 1.5,
            f"{value:.2f}%",
            ha="center",
            va="bottom",
            fontsize=9,
        )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()

    return str(output_path)


def make_qc_match_rate_by_component_plot(
    general_data: Dict[str, Any],
    figures_dir: str | Path,
    output_filename: str = "qc_match_rate_by_component.png",
    title: str = "QC concordance by component relative to the gold standard",
) -> str:
    components = general_data.get("components", {})

    component_names = []
    match_counts = []
    discrepancy_counts = []

    for comp_code, comp_data in components.items():
        matches, discrepancies = qc_hits_discrepancies_from_component(comp_data)

        component_names.append(comp_code)
        match_counts.append(matches)
        discrepancy_counts.append(discrepancies)

    output_dir = ensure_network_figures_dir(figures_dir)
    output_path = output_dir / output_filename

    plt.figure(figsize=(10, 6))
    x_positions = list(range(len(component_names)))

    plt.bar(x_positions, match_counts, label="Match", color=CBF_COLORS["match"])
    plt.bar(x_positions, discrepancy_counts, bottom=match_counts, label="Discrepancy", color=CBF_COLORS["discrepancy"])

    plt.xticks(x_positions, component_names)
    plt.xlabel("Component")
    plt.ylabel("Total number of QC evaluations")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()

    return str(output_path)


def generate_network_figures(
    general_data: Dict[str, Any],
    labs: List[Dict[str, Any]],
    figures_dir: str | Path = "figures",
) -> Dict[str, str]:
    outputs = {}

    outputs["consensus_summary"] = make_consensus_summary_plot(
        labs=labs,
        figures_dir=figures_dir,
    )

    outputs["variant_summary"] = make_variant_summary_plot(
        labs=labs,
        figures_dir=figures_dir,
    )

    outputs["sars_variant_reporting_summary"] = make_sars_variant_reporting_summary_plot(
        general_data=general_data,
        figures_dir=figures_dir,
    )

    outputs["influenza_variant_reporting_summary"] = make_influenza_variant_reporting_summary_plot(
        general_data=general_data,
        figures_dir=figures_dir,
    )

    outputs["classification_summary_lineage_type"] = make_stacked_classification_plot(
        general_data=general_data,
        figures_dir=figures_dir,
        mode="lineage_type",
        output_filename="classification_summary_lineage_type.png",
        title="Network-level lineage/type assignment performance summary",
    )

    outputs["classification_summary_clade"] = make_stacked_classification_plot(
        general_data=general_data,
        figures_dir=figures_dir,
        mode="clade",
        output_filename="classification_summary_clade.png",
        title="Network-level clade assignment performance summary",
    )

    outputs["metadata_completeness_distribution"] = make_metadata_completeness_distribution_plot(
        labs=labs,
        figures_dir=figures_dir,
    )

    outputs["qc_match_rate_by_component"] = make_qc_match_rate_by_component_plot(
        general_data=general_data,
        figures_dir=figures_dir,
    )

    return outputs


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


def collect_metadata_completeness_by_component(labs: List[Dict[str, Any]]) -> Dict[str, List[float]]:
    """
    Collect sample-level metadata completeness percentages grouped by component.
    Each value corresponds to one sample from one laboratory.
    """
    completeness_by_component: Dict[str, List[float]] = defaultdict(list)

    for lab in labs:
        for comp_code, comp in lab.get("components", {}).items():
            for sample in comp.get("samples", {}).values():
                total_expected = safe_number(sample.get("total_expected_fields"))
                filled = safe_number(sample.get("filled_fields"))

                if total_expected is None or total_expected == 0:
                    continue
                if filled is None:
                    continue

                completeness_pct = 100.0 * filled / total_expected
                completeness_by_component[comp_code].append(completeness_pct)

    return completeness_by_component


def make_metadata_completeness_distribution_plot(
    labs: List[Dict[str, Any]],
    figures_dir: str | Path,
) -> str:
    """
    Create a boxplot of sample-level metadata completeness percentages by component.
    """
    output_dir = ensure_network_figures_dir(figures_dir)
    output_path = output_dir / "metadata_completeness_distribution.png"

    completeness_by_component = collect_metadata_completeness_by_component(labs)

    component_order = ["SARS1", "SARS2", "FLU1", "FLU2"]
    component_names = [comp for comp in component_order if comp in completeness_by_component]
    data = [completeness_by_component[comp] for comp in component_names]

    if not data:
        return str(output_path)

    plt.figure(figsize=(10, 6))
    bp = plt.boxplot(data, labels=component_names, patch_artist=True)
    style_boxplot(bp, component_names)

    plt.xlabel("Component")
    plt.ylabel("Metadata completeness (%)")
    plt.title("Distribution of metadata completeness across participating laboratories")
    plt.ylim(0, 100)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()

    return str(output_path)


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
    influenza_variants_in_consensus: List[float] = []
    influenza_variants_in_consensus_vcf: List[float] = []
    influenza_discrepancies_in_reported_variants: List[float] = []
    influenza_variants_in_vcf: List[float] = []
    distinct_sars_references: set = set()
    distinct_influenza_references: set = set()

    sars_lineage_matches = 0
    sars_lineage_total = 0
    sars_clade_matches = 0
    sars_clade_total = 0
    flu_type_matches = 0
    flu_type_total = 0
    flu_clade_matches = 0
    flu_clade_total = 0

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

                coverage_threshold_total += 1
                if is_meaningful(sample.get("software_benchmarking", {}).get("depth_of_coverage_threshold")):
                    coverage_threshold_count += 1

                reference_genome_total += 1
                if is_meaningful(sample.get("metadata_metrics", {}).get("reference_genome_accession")):
                    reference_genome_count += 1

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

                sig = software_signature(
                    sample.get("software_benchmarking", {}).get("clade_assignment_software_name"),
                    sample.get("software_benchmarking", {}).get("clade_assignment_software_version"),
                )

                if sig:
                    variant_softwares.add(sig)

                expected_qc = expected_sample.get("expected_qc")
                reported_qc = sample.get("qc_test")
                qc_val = sample.get("qc_match")

                if expected_qc is not None and reported_qc is not None and qc_val is not None:
                    if qc_val is True:
                        qc_matches += 1
                        qc_total += 1
                    elif qc_val is False:
                        qc_discrepancies += 1
                        qc_total += 1

                gi = safe_number(sample.get("consensus", {}).get("genome_identity_pct"))
                if gi is not None:
                    if comp_expected.get("sequencing_instrument_platform") == "Illumina":
                        consensus_illumina_identity.append(gi)
                    else:
                        consensus_nanopore_identity.append(gi)

                cls = sample.get("classification", {})

                expected_lineage = cls.get("expected_lineage")
                expected_clade = cls.get("expected_clade")
                lineage_match = cls.get("lineage_match")
                clade_match = cls.get("clade_match")

                if comp_expected.get("virus") == "SARS-CoV-2":
                    if expected_lineage is not None and lineage_match is not None:
                        sars_lineage_total += 1
                        if lineage_match:
                            sars_lineage_matches += 1

                    if expected_clade is not None and clade_match is not None:
                        sars_clade_total += 1
                        if clade_match:
                            sars_clade_matches += 1
                else:
                    if expected_lineage is not None and lineage_match is not None:
                        flu_type_total += 1
                        if lineage_match:
                            flu_type_matches += 1

                    if expected_clade is not None and clade_match is not None:
                        flu_clade_total += 1
                        if clade_match:
                            flu_clade_matches += 1

                var = sample.get("variants", {})
                if comp_expected.get("virus") == "SARS-CoV-2":
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
                    nivc = safe_number(var.get("number_of_variants_in_consensus"))
                    nivcv = safe_number(var.get("number_of_variants_in_consensus_vcf"))
                    dirv = safe_number(var.get("discrepancies_in_reported_variants"))
                    niv = safe_number(var.get("number_of_variants_in_vcf"))

                    if nivc is not None:
                        influenza_variants_in_consensus.append(nivc)
                    if nivcv is not None:
                        influenza_variants_in_consensus_vcf.append(nivcv)
                    if dirv is not None:
                        influenza_discrepancies_in_reported_variants.append(dirv)
                    if niv is not None:
                        influenza_variants_in_vcf.append(niv)

                    if var.get("high_and_low_freq"):
                        influenza_variant_reporting_modes["high_and_low_freq"] += 1
                    elif var.get("high_freq_only"):
                        influenza_variant_reporting_modes["high_freq_only"] += 1
                    elif var.get("low_freq_only"):
                        influenza_variant_reporting_modes["low_freq_only"] += 1

                ref = sample.get("metadata_metrics", {}).get("reference_genome_accession")

                if is_meaningful(ref):
                    ref_list = [r.strip().upper() for r in str(ref).split(",") if r.strip()]

                    if comp_expected.get("virus") == "SARS-CoV-2":
                        for accession in ref_list:
                            distinct_sars_references.add(accession)
                    elif comp_expected.get("virus") == "Influenza virus":
                        for accession in ref_list:
                            distinct_influenza_references.add(accession)

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
                cons = sample.get("consensus", {})
                gi = safe_number(cons.get("genome_identity_pct"))
                if gi is None:
                    continue
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
            gis = []
            tds = []
            sample_bd = defaultdict(list)
            for lab in participating_labs:
                sample = lab["components"][comp_code]["samples"].get(sample_id)
                if not sample:
                    gis.append(None)
                    tds.append(None)
                    for key in ["wrong_nt", "ambiguity2nt", "nt2ambiguity", "ns2nt", "nt2ns", "insertions", "deletions"]:
                        sample_bd[key].append(None)
                    continue
                cons = sample.get("consensus", {})
                gi = safe_number(cons.get("genome_identity_pct"))
                if gi is None:
                    gis.append(None)
                    tds.append(None)
                    for key in ["wrong_nt", "ambiguity2nt", "nt2ambiguity", "ns2nt", "nt2ns", "insertions", "deletions"]:
                        sample_bd[key].append(None)
                    continue
                td = safe_number(cons.get("total_discrepancies"))
                gis.append(gi)
                tds.append(td)
                bd = cons.get("discrepancy_breakdown", {})
                for key in ["wrong_nt", "ambiguity2nt", "nt2ambiguity", "ns2nt", "nt2ns", "insertions", "deletions"]:
                    sample_bd[key].append(bd.get(key))

            identity_summary = summarize_numeric_values(gis, total_count=len(participating_labs))
            discrepancy_summary = summarize_numeric_values(tds, total_count=len(participating_labs))
            breakdown_summary = {
                key: summarize_numeric_values(vals, total_count=len(participating_labs))
                for key, vals in sample_bd.items()
            }
            sample_entries.append({
                "collecting_lab_sample_id": sample_id,
                "median_identity_pct": identity_summary["median"],
                "identity_pct_support": identity_summary,
                "median_discrepancies": discrepancy_summary["median"],
                "min": discrepancy_summary["min"],
                "max": discrepancy_summary["max"],
                "total_discrepancies_support": discrepancy_summary,
                "wrong_nt": breakdown_summary["wrong_nt"]["median"],
                "ambiguity2nt": breakdown_summary["ambiguity2nt"]["median"],
                "nt2ambiguity": breakdown_summary["nt2ambiguity"]["median"],
                "ns2nt": breakdown_summary["ns2nt"]["median"],
                "nt2ns": breakdown_summary["nt2ns"]["median"],
                "insertions": breakdown_summary["insertions"]["median"],
                "deletions": breakdown_summary["deletions"]["median"],
                "discrepancy_breakdown_support": breakdown_summary,
            })

        total_consensus_slots = len(participating_labs) * len(comp_expected["samples"])
        consensus_identity_summary = summarize_numeric_values(consensus_vals, total_count=total_consensus_slots)
        consensus_discrepancy_summary = summarize_numeric_values(discrepancy_vals, total_count=total_consensus_slots)
        consensus_breakdown_summary = {
            key: summarize_numeric_values(vals, total_count=total_consensus_slots)
            for key, vals in disc_breakdown_all.items()
        }
        comp_obj["consensus"] = {
            "median_identity": consensus_identity_summary["median"],
            "median_discrepancies": consensus_discrepancy_summary["median"],
            "total_median_discrepancies": round(sum(v for v in [consensus_discrepancy_summary["median"]] if v is not None), 4) if consensus_discrepancy_summary["median"] is not None else None,
            "min_discrepancies": consensus_discrepancy_summary["min"],
            "max_discrepancies": consensus_discrepancy_summary["max"],
            "median_identity_pct": consensus_identity_summary["median"],
            "identity_pct_min": consensus_identity_summary["min"],
            "identity_pct_max": consensus_identity_summary["max"],
            "identity_pct_summary": consensus_identity_summary,
            "total_discrepancies_summary": consensus_discrepancy_summary,
            "samples": sample_entries,
            "discrepancy_breakdown": consensus_breakdown_summary,
            "dominant_discrepancy_pattern": dominant_metric_key(consensus_breakdown_summary),
            "most_frequent_discrepancy_pattern": dominant_metric_key(consensus_breakdown_summary),
            "fig_discrepancies_boxplot_by_sample": f"figures/{comp_code}/consensus_discrepancies_boxplot_by_sample.png",
            "fig_discrepancies_stacked_by_sample": f"figures/{comp_code}/consensus_discrepancies_stacked_by_sample.png",
            "fig_discrepancy_type_boxplot": f"figures/{comp_code}/consensus_discrepancy_type_boxplot.png",
        }

        if comp_expected.get("virus") == "SARS-CoV-2":
            variant_discs = []
            variant_successful_hits = []
            variant_breakdown_all = defaultdict(list)
            variant_samples = []

            for sample_id, expected_sample in comp_expected["samples"].items():
                tds = []
                successful_hits_vals = []
                bd = defaultdict(list)

                for lab in participating_labs:
                    sample = lab["components"][comp_code]["samples"].get(sample_id)
                    if not sample:
                        tds.append(None)
                        successful_hits_vals.append(None)
                        for key in ["wrong_nt", "insertions", "deletions", "missing", "denovo"]:
                            bd[key].append(None)
                        continue

                    var = sample.get("variants", {})

                    td = safe_number(var.get("total_discrepancies"))
                    sh = safe_number(var.get("successful_hits"))

                    tds.append(td)
                    if td is not None:
                        variant_discs.append(td)

                    successful_hits_vals.append(sh)
                    if sh is not None:
                        variant_successful_hits.append(sh)

                    for key in ["wrong_nt", "insertions", "deletions", "missing", "denovo"]:
                        v = safe_number(var.get(key))
                        bd[key].append(v)
                        if v is not None:
                            variant_breakdown_all[key].append(v)

                discrepancy_summary = summarize_numeric_values(tds, total_count=len(participating_labs))
                successful_hits_summary = summarize_numeric_values(successful_hits_vals, total_count=len(participating_labs))
                breakdown_summary = {
                    key: summarize_numeric_values(vals, total_count=len(participating_labs))
                    for key, vals in bd.items()
                }

                variant_samples.append({
                    "collecting_lab_sample_id": sample_id,
                    "median_discrepancies": discrepancy_summary["median"],
                    "min": discrepancy_summary["min"],
                    "max": discrepancy_summary["max"],
                    "total_discrepancies_support": discrepancy_summary,
                    "median_successful_hits": successful_hits_summary["median"],
                    "min_successful_hits": successful_hits_summary["min"],
                    "max_successful_hits": successful_hits_summary["max"],
                    "successful_hits_support": successful_hits_summary,
                    "wrong_nt": breakdown_summary["wrong_nt"]["median"],
                    "insertions": breakdown_summary["insertions"]["median"],
                    "deletions": breakdown_summary["deletions"]["median"],
                    "missing": breakdown_summary["missing"]["median"],
                    "denovo": breakdown_summary["denovo"]["median"],
                    "discrepancy_breakdown_support": breakdown_summary,
                })

            variant_discrepancy_summary = summarize_numeric_values(variant_discs)
            successful_hits_component_summary = summarize_numeric_values(variant_successful_hits)
            variant_breakdown_summary = {
                key: summarize_numeric_values(vals)
                for key, vals in variant_breakdown_all.items()
            }
            comp_obj["variant"] = {
                "median_discrepancies": variant_discrepancy_summary["median"],
                "total_median_discrepancies": round(sum(v for v in [variant_discrepancy_summary["median"]] if v is not None), 4) if variant_discrepancy_summary["median"] is not None else None,
                "min_discrepancies": variant_discrepancy_summary["min"],
                "max_discrepancies": variant_discrepancy_summary["max"],
                "total_discrepancies_summary": variant_discrepancy_summary,
                "median_successful_hits": successful_hits_component_summary["median"],
                "min_successful_hits": successful_hits_component_summary["min"],
                "max_successful_hits": successful_hits_component_summary["max"],
                "successful_hits_summary": successful_hits_component_summary,
                "samples": variant_samples,
                "discrepancy_breakdown": variant_breakdown_summary,
                "dominant_discrepancy_pattern": dominant_metric_key(variant_breakdown_summary),
                "most_frequent_discrepancy_pattern": dominant_metric_key(variant_breakdown_summary),
                "fig_discrepancies_boxplot_by_sample": f"figures/{comp_code}/variant_discrepancies_boxplot_by_sample.png",
                "fig_discrepancy_type_boxplot": f"figures/{comp_code}/variant_discrepancy_type_boxplot.png",
            }
        else:
            variant_in_consensus_all = []
            variant_in_consensus_vcf_all = []
            discrepancy_reported_all = []
            variant_in_vcf_all = []
            variant_samples = []

            for sample_id, expected_sample in comp_expected["samples"].items():
                sample_variants_in_consensus = []
                sample_variants_in_consensus_vcf = []
                sample_discrepancies_in_reported_variants = []
                sample_variants_in_vcf = []

                for lab in participating_labs:
                    sample = lab["components"][comp_code]["samples"].get(sample_id)
                    if not sample:
                        sample_variants_in_consensus.append(None)
                        sample_variants_in_consensus_vcf.append(None)
                        sample_discrepancies_in_reported_variants.append(None)
                        sample_variants_in_vcf.append(None)
                        continue

                    var = sample.get("variants", {})

                    nivc = safe_number(var.get("number_of_variants_in_consensus"))
                    nivcv = safe_number(var.get("number_of_variants_in_consensus_vcf"))
                    dirv = safe_number(var.get("discrepancies_in_reported_variants"))
                    niv = safe_number(var.get("number_of_variants_in_vcf"))

                    sample_variants_in_consensus.append(nivc)
                    if nivc is not None:
                        variant_in_consensus_all.append(nivc)
                    sample_variants_in_consensus_vcf.append(nivcv)
                    if nivcv is not None:
                        variant_in_consensus_vcf_all.append(nivcv)
                    sample_discrepancies_in_reported_variants.append(dirv)
                    if dirv is not None:
                        discrepancy_reported_all.append(dirv)
                    sample_variants_in_vcf.append(niv)
                    if niv is not None:
                        variant_in_vcf_all.append(niv)

                variants_in_consensus_summary = summarize_numeric_values(sample_variants_in_consensus, total_count=len(participating_labs))
                variants_in_consensus_vcf_summary = summarize_numeric_values(sample_variants_in_consensus_vcf, total_count=len(participating_labs))
                reported_variant_discrepancies_summary = summarize_numeric_values(sample_discrepancies_in_reported_variants, total_count=len(participating_labs))
                variants_in_vcf_summary = summarize_numeric_values(sample_variants_in_vcf, total_count=len(participating_labs))
                variant_samples.append({
                    "collecting_lab_sample_id": sample_id,
                    "variants_in_consensus": variants_in_consensus_summary,
                    "variants_in_consensus_vcf": variants_in_consensus_vcf_summary,
                    "reported_variant_discrepancies": reported_variant_discrepancies_summary,
                    "discrepancies_in_reported_variants": reported_variant_discrepancies_summary,
                    "variants_in_vcf": variants_in_vcf_summary,
                })

            variants_in_consensus_summary = summarize_numeric_values(variant_in_consensus_all)
            variants_in_consensus_vcf_summary = summarize_numeric_values(variant_in_consensus_vcf_all)
            reported_variant_discrepancies_summary = summarize_numeric_values(discrepancy_reported_all)
            variants_in_vcf_summary = summarize_numeric_values(variant_in_vcf_all)
            comp_obj["variant"] = {
                "median_variants_in_consensus": variants_in_consensus_summary["median"],
                "min_variants_in_consensus": variants_in_consensus_summary["min"],
                "max_variants_in_consensus": variants_in_consensus_summary["max"],
                "variants_in_consensus_summary": variants_in_consensus_summary,
                "median_variants_in_consensus_vcf": variants_in_consensus_vcf_summary["median"],
                "min_variants_in_consensus_vcf": variants_in_consensus_vcf_summary["min"],
                "max_variants_in_consensus_vcf": variants_in_consensus_vcf_summary["max"],
                "variants_in_consensus_vcf_summary": variants_in_consensus_vcf_summary,
                "median_discrepancies_in_reported_variants": reported_variant_discrepancies_summary["median"],
                "min_discrepancies_in_reported_variants": reported_variant_discrepancies_summary["min"],
                "max_discrepancies_in_reported_variants": reported_variant_discrepancies_summary["max"],
                "reported_variant_discrepancies_summary": reported_variant_discrepancies_summary,
                "median_variants_in_vcf": variants_in_vcf_summary["median"],
                "min_variants_in_vcf": variants_in_vcf_summary["min"],
                "max_variants_in_vcf": variants_in_vcf_summary["max"],
                "variants_in_vcf_summary": variants_in_vcf_summary,
                "samples": variant_samples,
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
                "total_distinct_references": len(distinct_sars_references),
            },
            "influenza_variants": {
                "high_and_low_freq_pct": pct(influenza_variant_reporting_modes["high_and_low_freq"], sum(influenza_variant_reporting_modes.values())),
                "low_freq_only_pct": pct(influenza_variant_reporting_modes["low_freq_only"], sum(influenza_variant_reporting_modes.values())),
                "high_freq_only_pct": pct(influenza_variant_reporting_modes["high_freq_only"], sum(influenza_variant_reporting_modes.values())),
                "total_distinct_fragments": len(distinct_influenza_references),
                "total_distinct_references": round(len(distinct_influenza_references) / 8.0, 4) if distinct_influenza_references else 0,
                "median_variants_in_consensus": median_or_none(influenza_variants_in_consensus),
                "min_variants_in_consensus": min_or_none(influenza_variants_in_consensus),
                "max_variants_in_consensus": max_or_none(influenza_variants_in_consensus),
                "median_variants_in_consensus_vcf": median_or_none(influenza_variants_in_consensus_vcf),
                "min_variants_in_consensus_vcf": min_or_none(influenza_variants_in_consensus_vcf),
                "max_variants_in_consensus_vcf": max_or_none(influenza_variants_in_consensus_vcf),
                "median_discrepancies_in_reported_variants": median_or_none(influenza_discrepancies_in_reported_variants),
                "min_discrepancies_in_reported_variants": min_or_none(influenza_discrepancies_in_reported_variants),
                "max_discrepancies_in_reported_variants": max_or_none(influenza_discrepancies_in_reported_variants),
                "median_variants_in_vcf": median_or_none(influenza_variants_in_vcf),
                "min_variants_in_vcf": min_or_none(influenza_variants_in_vcf),
                "max_variants_in_vcf": max_or_none(influenza_variants_in_vcf),
            },
            "classification": {
                "sars_cov_2_concordance_pct": pct(sars_lineage_matches, sars_lineage_total),
                "influenza_type_concordance_pct": pct(flu_type_matches, flu_type_total),
                "sars_clade_concordance_pct": pct(sars_clade_matches, sars_clade_total),
                "flu_clade_concordance_pct": pct(flu_clade_matches, flu_clade_total),
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
    parser.add_argument("--figures-dir", default="figures", help="Base directory where figures will be written")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    expected_data = load_json(args.expected_data)
    labs = load_lab_jsons(Path(args.labs_dir))
    general = build_general(expected_data, labs)

    generated_figures = generate_network_figures(
        general_data=general,
        labs=labs,
        figures_dir=args.figures_dir,
    )

    # Update figure paths in general.json
    general.setdefault("figures", {})
    general["figures"].update(generated_figures)

    dump_json(general, args.output)
    print(f"Generated {args.output}")

if __name__ == "__main__":
    main()
