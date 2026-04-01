#!/usr/bin/env python3
import argparse
from html import parser
import json
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median
from typing import Any, Dict, Iterable, List, Optional

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import transforms as mtransforms
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
    "null": "#999999",
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

LAB_HIGHLIGHT_FACE_COLOR = "#000000"
LAB_HIGHLIGHT_EDGE_COLOR = (1.0, 1.0, 1.0, 0.8)

COMPONENT_CONSENSUS_SAMPLE_Y_LIMITS = {
    "SARS1": None,
    "SARS2": None,
    "FLU1": None,
    "FLU2": 410.0,
}

COMPONENT_CONSENSUS_TYPE_BOXPLOT_Y_LIMITS = {
    "FLU2": 400.0,
}

CONSENSUS_DISCREPANCY_TYPE_ORDER = [
    "wrong_nt",
    "ambiguity2nt",
    "nt2ambiguity",
    "ns2nt",
    "nt2ns",
    "insertions",
    "deletions",
]

CONSENSUS_DISCREPANCY_TYPE_LABELS = {
    "wrong_nt": "Wrong nucleotide",
    "ambiguity2nt": "Ambiguity instead\nof nucleotide",
    "nt2ambiguity": "Nucleotide instead\nof ambiguity",
    "ns2nt": "Stretch of Ns\ninstead of nucleotide",
    "nt2ns": "Nucleotide stretch\ninstead of Ns",
    "insertions": "Insertion relative\nto gold standard",
    "deletions": "Deletion relative\nto gold standard",
}

CONSENSUS_DISCREPANCY_TYPE_COLORS = {
    "wrong_nt": "#0072B2",
    "ambiguity2nt": "#56B4E9",
    "nt2ambiguity": "#009E73",
    "ns2nt": "#F0E442",
    "nt2ns": "#E69F00",
    "insertions": "#D55E00",
    "deletions": "#CC79A7",
}

VARIANT_DISCREPANCY_TYPE_ORDER = [
    "wrong_nt",
    "insertions",
    "deletions",
    "missing",
    "denovo",
]

VARIANT_DISCREPANCY_TYPE_LABELS = {
    "wrong_nt": "Wrong nucleotide",
    "insertions": "Insertions",
    "deletions": "Deletions",
    "missing": "Missing expected\nvariants",
    "denovo": "De novo\nvariants",
}

VARIANT_DISCREPANCY_TYPE_COLORS = {
    "wrong_nt": "#0072B2",
    "insertions": "#E69F00",
    "deletions": "#CC79A7",
    "missing": "#56B4E9",
    "denovo": "#D55E00",
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


def software_signature(name: Any, version: Any = None, database_version: Any = None) -> Optional[tuple]:
    if not is_meaningful(name):
        return None
    clean_name = str(name).strip()
    if "custom" in clean_name.lower():
        clean_version = None
        clean_database_version = None
    else:
        clean_version = str(version).strip() if is_meaningful(version) else None
        clean_database_version = str(database_version).strip() if is_meaningful(database_version) else None
    return (clean_name, clean_version, clean_database_version)


def most_common_or_none(values: Iterable[Any]) -> Any:
    clean = [v for v in values if is_meaningful(v)]
    if not clean:
        return None
    return Counter(clean).most_common(1)[0][0]


def ensure_network_figures_dir(figures_dir: str | Path) -> Path:
    output_dir = Path(figures_dir) / "network"
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


def ensure_component_figures_dir(figures_dir: str | Path, comp_code: str) -> Path:
    output_dir = Path(figures_dir) / comp_code
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


def ensure_lab_component_figures_dir(figures_dir: str | Path, lab_code: str, comp_code: str) -> Path:
    output_dir = Path(figures_dir) / "labs" / lab_code / comp_code
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


def get_lab_identifier(lab: Dict[str, Any]) -> str:
    lab_meta = lab.get("lab", {})
    return (
        str(lab_meta.get("lab_cod")).strip()
        if is_meaningful(lab_meta.get("lab_cod"))
        else str(lab_meta.get("submitting_institution_id", "unknown_lab")).strip()
    )


def style_boxplot_axes(ax: Any) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(True)
    ax.spines["bottom"].set_visible(True)


def style_percent_boxplot_axis(ax: Any) -> None:
    ax.set_ylim(0, 110)
    ax.set_yticks(np.arange(0, 101, 20))


def annotate_group_n_labs(ax: Any, positions: List[float], counts: List[int], y_offset: float = -0.22) -> None:
    for x_pos, count in zip(positions, counts):
        ax.text(
            x_pos,
            y_offset,
            f"n={count}",
            transform=ax.get_xaxis_transform(),
            ha="center",
            va="top",
            fontsize=8,
            color="#444444",
            clip_on=False,
        )


def style_boxplot(bp: Dict[str, Any], labels: List[str], ax: Optional[Any] = None) -> None:
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
    style_boxplot_axes(ax if ax is not None else plt.gca())


def add_component_boxplot_points(
    ax: Any,
    bp: Dict[str, Any],
    data: List[List[float]],
    positions: List[float],
    labels: List[str],
) -> None:
    jittered_positions = build_jittered_positions(data, positions)
    for xs, values, label in zip(jittered_positions, data, labels):
        if not values:
            continue
        ax.scatter(
            xs,
            values,
            s=16,
            color=COMPONENT_BOX_COLORS.get(label, CBF_COLORS["box_default"]),
            alpha=0.35,
            edgecolors="none",
            zorder=3,
        )
    align_boxplot_fliers(bp, data, jittered_positions)


def add_colored_boxplot_points(
    ax: Any,
    bp: Dict[str, Any],
    data: List[List[float]],
    positions: List[float],
    colors: List[str],
) -> None:
    jittered_positions = build_jittered_positions(data, positions)
    for xs, values, color in zip(jittered_positions, data, colors):
        if not values:
            continue
        ax.scatter(
            xs,
            values,
            s=16,
            color=color,
            alpha=0.35,
            edgecolors="none",
            zorder=3,
        )
    align_boxplot_fliers(bp, data, jittered_positions)


def build_jittered_positions(
    data: List[List[float]],
    positions: List[float],
    jitter_width: float = 0.04,
) -> List[List[float]]:
    rng = np.random.default_rng(42)
    jittered_positions: List[List[float]] = []
    for pos, values in zip(positions, data):
        if not values:
            jittered_positions.append([])
            continue
        jitter = rng.uniform(-jitter_width, jitter_width, size=len(values))
        jittered_positions.append(list(np.full(len(values), pos) + jitter))
    return jittered_positions


def align_boxplot_fliers(
    bp: Dict[str, Any],
    data: List[List[float]],
    jittered_positions: List[List[float]],
) -> None:
    for flier, values, xs in zip(bp["fliers"], data, jittered_positions):
        flier_y = list(flier.get_ydata())
        if not flier_y or not values:
            continue

        centered_x = list(flier.get_xdata())
        used_indices = set()
        aligned_x = []

        for y_value in flier_y:
            matched_index = None

            for idx, value in enumerate(values):
                if idx in used_indices:
                    continue
                if value == y_value or np.isclose(value, y_value, rtol=0.0, atol=1e-9):
                    matched_index = idx
                    break

            if matched_index is None:
                aligned_x.append(centered_x[0] if centered_x else 0.0)
                continue

            used_indices.add(matched_index)
            aligned_x.append(xs[matched_index])

        flier.set_xdata(aligned_x)


def add_lab_result_diamond(
    ax: Any,
    positions: List[float],
    values: List[Optional[float]],
    y_upper: Optional[float] = None,
) -> None:
    xs = []
    ys = []
    clipped_annotations: List[tuple[float, float]] = []
    for pos, value in zip(positions, values):
        numeric_value = safe_number(value)
        if numeric_value is None:
            continue
        plotted_value = numeric_value
        if y_upper is not None and numeric_value > y_upper:
            plotted_value = y_upper * 0.97
            clipped_annotations.append((pos, numeric_value))
        xs.append(pos)
        ys.append(plotted_value)

    if not xs:
        return []

    ax.scatter(
        xs,
        ys,
        s=64,
        marker="D",
        facecolor=LAB_HIGHLIGHT_FACE_COLOR,
        edgecolor=LAB_HIGHLIGHT_EDGE_COLOR,
        linewidth=1.4,
        zorder=6,
    )
    return clipped_annotations


def collect_classification_sample_outcomes(
    general_data: Dict[str, Any],
    labs: List[Dict[str, Any]],
    comp_code: str,
    mode: str,
) -> Dict[str, Any]:
    comp_data = general_data.get("components", {}).get(comp_code, {})
    sample_entries = comp_data.get("typing", {}).get("samples", [])
    participating_labs = [lab for lab in labs if comp_code in lab.get("components", {})]
    total_labs = len(participating_labs)

    if mode == "lineage_type":
        expected_key = "expected_lineage"
        assignment_key = "lineage_assignment"
        match_key = "lineage_match"
    elif mode == "clade":
        expected_key = "expected_clade"
        assignment_key = "clade_assignment"
        match_key = "clade_match"
    else:
        raise ValueError(f"Unknown classification mode: {mode}")

    sample_outcomes = []
    total_hits = 0
    total_discrepancies = 0
    total_nulls = 0

    for sample_entry in sample_entries:
        sample_id = sample_entry.get("collecting_lab_sample_id")
        if not sample_id:
            continue

        expected_value = None
        for lab in participating_labs:
            sample = lab.get("components", {}).get(comp_code, {}).get("samples", {}).get(sample_id)
            if not sample:
                continue
            candidate = sample.get("classification", {}).get(expected_key)
            if is_meaningful(candidate):
                expected_value = candidate
                break

        if not is_meaningful(expected_value):
            continue

        hits = 0
        discrepancies = 0
        nulls = 0

        for lab in participating_labs:
            sample = lab.get("components", {}).get(comp_code, {}).get("samples", {}).get(sample_id)
            cls = sample.get("classification", {}) if sample else {}
            assignment_value = cls.get(assignment_key)
            match_value = cls.get(match_key)

            if not is_meaningful(assignment_value):
                nulls += 1
            elif match_value is True:
                hits += 1
            else:
                discrepancies += 1

        total_hits += hits
        total_discrepancies += discrepancies
        total_nulls += nulls
        sample_outcomes.append({
            "collecting_lab_sample_id": sample_id,
            "hits": hits,
            "discrepancies": discrepancies,
            "nulls": nulls,
            "total": total_labs,
            "hit_pct": (100.0 * hits / total_labs) if total_labs else 0.0,
            "discrepancy_pct": (100.0 * discrepancies / total_labs) if total_labs else 0.0,
            "null_pct": (100.0 * nulls / total_labs) if total_labs else 0.0,
        })

    return {
        "samples": sample_outcomes,
        "hits": total_hits,
        "discrepancies": total_discrepancies,
        "nulls": total_nulls,
        "total_possible": total_labs * len(sample_outcomes),
        "total_labs": total_labs,
    }


def make_stacked_classification_plot(
    general_data: Dict[str, Any],
    labs: List[Dict[str, Any]],
    figures_dir: str | Path,
    mode: str,
    output_filename: str,
    title: str,
) -> str:
    components = general_data.get("components", {})

    component_names = []
    hit_pcts = []
    discrepancy_pcts = []
    null_pcts = []

    for comp_code in components.keys():
        outcomes = collect_classification_sample_outcomes(general_data, labs, comp_code, mode)
        hits = outcomes["hits"]
        discrepancies = outcomes["discrepancies"]
        nulls = outcomes["nulls"]
        total_possible = outcomes["total_possible"]

        component_names.append(comp_code)
        hit_pcts.append(100.0 * hits / total_possible if total_possible else 0.0)
        discrepancy_pcts.append(100.0 * discrepancies / total_possible if total_possible else 0.0)
        null_pcts.append(100.0 * nulls / total_possible if total_possible else 0.0)

    output_dir = ensure_network_figures_dir(figures_dir)
    output_path = output_dir / output_filename

    plt.figure(figsize=(10, 6))
    x_positions = list(range(len(component_names)))

    plt.bar(x_positions, hit_pcts, label="Match", color=CBF_COLORS["match"])
    plt.bar(x_positions, discrepancy_pcts, bottom=hit_pcts, label="Discrepancy", color=CBF_COLORS["discrepancy"])
    plt.bar(
        x_positions,
        null_pcts,
        bottom=np.array(hit_pcts) + np.array(discrepancy_pcts),
        label="Not provided",
        color=CBF_COLORS["null"],
    )

    plt.xticks(x_positions, component_names)
    plt.xlabel("Component")
    plt.ylabel("Assignments (%)")
    plt.ylim(0, 100)
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()

    return str(output_path)


def make_combined_classification_summary_plot(
    general_data: Dict[str, Any],
    labs: List[Dict[str, Any]],
    figures_dir: str | Path,
    output_filename: str = "classification_summary.png",
) -> str:
    components = general_data.get("components", {})
    component_names = list(components.keys())

    lineage_hit_pcts = []
    lineage_discrepancy_pcts = []
    lineage_null_pcts = []
    clade_hit_pcts = []
    clade_discrepancy_pcts = []
    clade_null_pcts = []

    for comp_code in components.keys():
        lineage_outcomes = collect_classification_sample_outcomes(general_data, labs, comp_code, "lineage_type")
        clade_outcomes = collect_classification_sample_outcomes(general_data, labs, comp_code, "clade")
        lineage_hits = lineage_outcomes["hits"]
        lineage_discrepancies = lineage_outcomes["discrepancies"]
        lineage_nulls = lineage_outcomes["nulls"]
        lineage_total_possible = lineage_outcomes["total_possible"]
        clade_hits = clade_outcomes["hits"]
        clade_discrepancies = clade_outcomes["discrepancies"]
        clade_nulls = clade_outcomes["nulls"]
        clade_total_possible = clade_outcomes["total_possible"]

        component_names.append(comp_code)
        lineage_hit_pcts.append(100.0 * lineage_hits / lineage_total_possible if lineage_total_possible else 0.0)
        lineage_discrepancy_pcts.append(100.0 * lineage_discrepancies / lineage_total_possible if lineage_total_possible else 0.0)
        lineage_null_pcts.append(100.0 * lineage_nulls / lineage_total_possible if lineage_total_possible else 0.0)
        clade_hit_pcts.append(100.0 * clade_hits / clade_total_possible if clade_total_possible else 0.0)
        clade_discrepancy_pcts.append(100.0 * clade_discrepancies / clade_total_possible if clade_total_possible else 0.0)
        clade_null_pcts.append(100.0 * clade_nulls / clade_total_possible if clade_total_possible else 0.0)

    # Remove duplicated names introduced while collecting counts.
    component_names = list(components.keys())
    x_positions = list(range(len(component_names)))

    output_dir = ensure_network_figures_dir(figures_dir)
    output_path = output_dir / output_filename

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

    lineage_match_bars = axes[0].bar(x_positions, lineage_hit_pcts, color=CBF_COLORS["match"], label="Match")
    lineage_discrepancy_bars = axes[0].bar(x_positions, lineage_discrepancy_pcts, bottom=lineage_hit_pcts, color=CBF_COLORS["discrepancy"], label="Discrepancy")
    lineage_null_bars = axes[0].bar(
        x_positions,
        lineage_null_pcts,
        bottom=np.array(lineage_hit_pcts) + np.array(lineage_discrepancy_pcts),
        color=CBF_COLORS["null"],
        label="Not provided",
    )
    axes[0].set_xticks(x_positions)
    axes[0].set_xticklabels(component_names)
    axes[0].set_xlabel("Component")
    axes[0].set_ylabel("Assignments (%)")
    axes[0].set_ylim(0, 100)
    axes[0].set_title("A. Lineage/type assignments")

    clade_match_bars = axes[1].bar(x_positions, clade_hit_pcts, color=CBF_COLORS["match"], label="Match")
    clade_discrepancy_bars = axes[1].bar(x_positions, clade_discrepancy_pcts, bottom=clade_hit_pcts, color=CBF_COLORS["discrepancy"], label="Discrepancy")
    clade_null_bars = axes[1].bar(
        x_positions,
        clade_null_pcts,
        bottom=np.array(clade_hit_pcts) + np.array(clade_discrepancy_pcts),
        color=CBF_COLORS["null"],
        label="Not provided",
    )
    axes[1].set_xticks(x_positions)
    axes[1].set_xticklabels(component_names)
    axes[1].set_xlabel("Component")
    axes[1].set_ylim(0, 100)
    axes[1].set_title("B. Clade assignments")

    for ax, match_bars, discrepancy_bars, null_bars, match_pcts, discrepancy_pcts, null_pcts in [
        (axes[0], lineage_match_bars, lineage_discrepancy_bars, lineage_null_bars, lineage_hit_pcts, lineage_discrepancy_pcts, lineage_null_pcts),
        (axes[1], clade_match_bars, clade_discrepancy_bars, clade_null_bars, clade_hit_pcts, clade_discrepancy_pcts, clade_null_pcts),
    ]:
        for bar, value in zip(match_bars, match_pcts):
            if value <= 0:
                continue
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                value / 2,
                f"{value:.1f}%",
                ha="center",
                va="center",
                fontsize=8,
                color="white",
                fontweight="bold",
            )

        for bar, match_value, discrepancy_value in zip(discrepancy_bars, match_pcts, discrepancy_pcts):
            if discrepancy_value <= 0:
                continue
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                match_value + discrepancy_value / 2,
                f"{discrepancy_value:.1f}%",
                ha="center",
                va="center",
                fontsize=8,
                color="white",
                fontweight="bold",
            )

        for bar, match_value, discrepancy_value, null_value in zip(null_bars, match_pcts, discrepancy_pcts, null_pcts):
            if null_value <= 0:
                continue
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                match_value + discrepancy_value + null_value / 2,
                f"{null_value:.1f}%",
                ha="center",
                va="center",
                fontsize=8,
                color="white",
                fontweight="bold",
            )

    fig.suptitle("Distribution of classification outcomes across participating laboratories")
    fig.legend(
        handles=[
            plt.Rectangle((0, 0), 1, 1, color=CBF_COLORS["match"]),
            plt.Rectangle((0, 0), 1, 1, color=CBF_COLORS["discrepancy"]),
            plt.Rectangle((0, 0), 1, 1, color=CBF_COLORS["null"]),
        ],
        labels=["Match", "Discrepancy", "Not provided"],
        loc="lower center",
        bbox_to_anchor=(0.5, 0.02),
        ncol=3,
        frameon=False,
    )
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    return str(output_path)


def make_component_typing_outcome_stacked_bar_by_sample(
    general_data: Dict[str, Any],
    labs: List[Dict[str, Any]],
    comp_code: str,
    figures_dir: str | Path,
    output_filename: str = "typing_outcome_stackedbar_by_sample.png",
) -> str:
    output_dir = ensure_component_figures_dir(figures_dir, comp_code)
    output_path = output_dir / output_filename

    lineage_outcomes = collect_classification_sample_outcomes(general_data, labs, comp_code, "lineage_type")
    clade_outcomes = collect_classification_sample_outcomes(general_data, labs, comp_code, "clade")
    if not lineage_outcomes["samples"] and not clade_outcomes["samples"]:
        return str(output_path)
    width = 0.62

    max_samples = max(len(lineage_outcomes["samples"]), len(clade_outcomes["samples"]), 1)
    fig, axes = plt.subplots(1, 2, figsize=(max(12, max_samples * 1.5), 6.6), sharey=True)

    panel_specs = [
        ("A. Lineage/Subtype assignments", lineage_outcomes),
        ("B. Clade assignments", clade_outcomes),
    ]

    for ax, (title, outcomes) in zip(axes, panel_specs):
        if not outcomes["samples"]:
            ax.set_visible(False)
            continue
        sample_names = [sample["collecting_lab_sample_id"] for sample in outcomes["samples"]]
        match_rates = [sample["hit_pct"] for sample in outcomes["samples"]]
        discrepancy_rates = [sample["discrepancy_pct"] for sample in outcomes["samples"]]
        null_rates = [sample["null_pct"] for sample in outcomes["samples"]]
        x_positions = np.arange(len(sample_names))

        match_bars = ax.bar(
            x_positions,
            match_rates,
            width=width,
            color=CBF_COLORS["match"],
            label="Match",
        )
        discrepancy_bars = ax.bar(
            x_positions,
            discrepancy_rates,
            width=width,
            bottom=match_rates,
            color=CBF_COLORS["discrepancy"],
            label="Discrepancy",
        )
        null_bars = ax.bar(
            x_positions,
            null_rates,
            width=width,
            bottom=np.array(match_rates) + np.array(discrepancy_rates),
            color=CBF_COLORS["null"],
            label="Not provided",
        )
        ax.set_xticks(x_positions)
        ax.set_xticklabels(sample_names, rotation=0, ha="center")
        ax.set_xlabel("Sample")
        ax.set_ylim(0, 100)
        ax.set_title(title)

        for bar, value in zip(match_bars, match_rates):
            if value is None or np.isnan(value) or value <= 0:
                continue
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                value / 2,
                f"{value:.1f}%",
                ha="center",
                va="center",
                fontsize=8,
                color="white",
                fontweight="bold",
            )

        for bar, match_value, discrepancy_value in zip(discrepancy_bars, match_rates, discrepancy_rates):
            if (
                match_value is None
                or discrepancy_value is None
                or np.isnan(match_value)
                or np.isnan(discrepancy_value)
                or discrepancy_value <= 0
            ):
                continue
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                match_value + discrepancy_value / 2,
                f"{discrepancy_value:.1f}%",
                ha="center",
                va="center",
                fontsize=8,
                color="white",
                fontweight="bold",
            )

        for bar, match_value, discrepancy_value, null_value in zip(null_bars, match_rates, discrepancy_rates, null_rates):
            if null_value is None or np.isnan(null_value) or null_value <= 0:
                continue
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                match_value + discrepancy_value + null_value / 2,
                f"{null_value:.1f}%",
                ha="center",
                va="center",
                fontsize=8,
                color="white",
                fontweight="bold",
            )

    for ax in axes:
        if ax.get_visible():
            ax.set_ylabel("Assignments (%)")
            break

    fig.legend(
        handles=[
            plt.Rectangle((0, 0), 1, 1, color=CBF_COLORS["match"]),
            plt.Rectangle((0, 0), 1, 1, color=CBF_COLORS["discrepancy"]),
            plt.Rectangle((0, 0), 1, 1, color=CBF_COLORS["null"]),
        ],
        labels=["Match", "Discrepancy", "Not provided"],
        loc="lower center",
        bbox_to_anchor=(0.5, 0.02),
        ncol=3,
        frameon=False,
    )
    fig.suptitle(f"{comp_code} classification outcome distribution by sample")
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    return str(output_path)


def make_component_qc_match_by_sample_plot(
    general_data: Dict[str, Any],
    comp_code: str,
    figures_dir: str | Path,
    output_filename: str = "qc_match_by_sample.png",
) -> str:
    output_dir = ensure_component_figures_dir(figures_dir, comp_code)
    output_path = output_dir / output_filename

    comp_data = general_data.get("components", {}).get(comp_code, {})
    samples = comp_data.get("qc", {}).get("samples", [])
    if not samples:
        return str(output_path)

    evaluable_samples = [
        sample for sample in samples
        if (sample.get("collecting_lab_sample_id") or sample.get("sample_id"))
        and safe_number(sample.get("match_rate_pct")) is not None
    ]
    if not evaluable_samples:
        return str(output_path)

    sample_names = [
        sample.get("collecting_lab_sample_id") or sample.get("sample_id")
        for sample in evaluable_samples
    ]

    match_rates = []
    discrepancy_rates = []

    for sample in evaluable_samples:
        match_rate = safe_number(sample.get("match_rate_pct"))
        match_rates.append(match_rate if match_rate is not None else np.nan)
        discrepancy_rates.append(100.0 - match_rate if match_rate is not None else np.nan)

    x_positions = np.arange(len(sample_names))
    width = 0.62

    plt.figure(figsize=(max(9, len(sample_names) * 1.35), 6))
    match_bars = plt.bar(
        x_positions,
        match_rates,
        width=width,
        color=CBF_COLORS["match"],
        label="Match",
    )
    discrepancy_bars = plt.bar(
        x_positions,
        discrepancy_rates,
        width=width,
        bottom=match_rates,
        color=CBF_COLORS["discrepancy"],
        label="Discrepancy",
    )

    plt.xticks(x_positions, sample_names, rotation=0, ha="center")
    plt.xlabel("Sample")
    plt.ylabel("QC evaluations (%)")
    plt.ylim(0, 100)
    plt.title(f"{comp_code} QC concordance by sample")
    plt.legend(frameon=False, loc="lower center", bbox_to_anchor=(0.5, -0.28), ncol=2)

    for bar, value in zip(match_bars, match_rates):
        if value is None or np.isnan(value) or value <= 0:
            continue
        plt.text(
            bar.get_x() + bar.get_width() / 2,
            value / 2,
            f"{value:.1f}%",
            ha="center",
            va="center",
            fontsize=8,
            color="white",
            fontweight="bold",
        )

    for bar, match_value, discrepancy_value in zip(discrepancy_bars, match_rates, discrepancy_rates):
        if (
            match_value is None
            or discrepancy_value is None
            or np.isnan(match_value)
            or np.isnan(discrepancy_value)
            or discrepancy_value <= 0
        ):
            continue
        plt.text(
            bar.get_x() + bar.get_width() / 2,
            match_value + discrepancy_value / 2,
            f"{discrepancy_value:.1f}%",
            ha="center",
            va="center",
            fontsize=8,
            color="white",
            fontweight="bold",
        )

    plt.tight_layout(rect=(0, 0.08, 1, 1))
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()

    return str(output_path)


def format_software_label(name: Any, version: Any, database_version: Any = None) -> str:
    clean_name = str(name).strip() if is_meaningful(name) else "Unknown"
    clean_version = str(version).strip() if is_meaningful(version) else "Version N/A"
    clean_database_version = str(database_version).strip() if is_meaningful(database_version) else None
    name_lower = clean_name.lower()
    version_lower = clean_version.lower()

    if name_lower != "bwa mem":
        clean_name = "\n".join(clean_name.split())

    if clean_version in {
        "DRAGEN Microbial Amplicon ad hoc version",
        "DRAGEN Microbial Amplicon ad-hoc version",
        "DRAGEN Microbial Amplicon",
    }:
        clean_version = "DRAGEN version"
        version_lower = clean_version.lower()

    if "custom" in name_lower:
        label = clean_name
    elif "dragen" in name_lower and "dragen" in version_lower:
        label = clean_name
    else:
        label = f"{clean_name}\n{clean_version}"

    if "custom" not in name_lower and clean_database_version is not None:
        formatted_database_version = clean_database_version.replace(",", ",\n")
        label = f"{label}\nDB version: {formatted_database_version}"
    return label


def format_benchmark_group_label(
    name: Any,
    version: Any,
    database_version: Any = None,
    n_labs: Any = None,
) -> str:
    label = format_software_label(name, version, database_version)
    n_labs_int = safe_int(n_labs)
    if n_labs_int is not None:
        label = f"{label}\n(n={n_labs_int})"
    return label


def extract_variant_reporting_mode_pct(sample: Dict[str, Any], mode_key: str) -> Optional[float]:
    var = sample.get("variants", {})
    mode_keys = ["high_and_low_freq", "high_freq_only", "low_freq_only"]
    mode_values = [var.get(key) for key in mode_keys]
    if all(value is None for value in mode_values):
        return None
    return 100.0 if var.get(mode_key) is True else 0.0


def make_component_bioinformatics_protocol_metric_boxplots(
    labs: List[Dict[str, Any]],
    comp_code: str,
    figures_dir: str | Path,
    output_filename: str = "bioinformatics_protocol_metric_boxplots_by_pipeline.png",
) -> str:
    output_dir = ensure_component_figures_dir(figures_dir, comp_code)
    output_path = output_dir / output_filename

    groups = collect_software_groups(
        labs,
        comp_code,
        "bioinformatics_protocol_software_name",
        "bioinformatics_protocol_software_version",
    )
    if not groups:
        return str(output_path)

    ordered_groups = sorted(
        groups.items(),
        key=lambda item: (
            str(item[0][0]).lower() if item[0][0] is not None else "",
            str(item[0][1]).lower() if item[0][1] is not None else "",
        ),
    )

    group_data = []

    for (name, version, database_version), records in ordered_groups:
        identities = []
        discrepancies = []
        metadata_completeness = []
        exact_classification = []
        clade_hits = 0
        clade_total = 0
        lineage_hits = 0
        lineage_total = 0

        for record in records:
            sample = record["sample"]
            cons = sample.get("consensus", {})
            cls = sample.get("classification", {})

            genome_identity = safe_number(cons.get("genome_identity_pct"))
            total_discrepancies = safe_number(cons.get("total_discrepancies"))
            filled_fields = safe_number(sample.get("filled_fields"))
            total_expected_fields = safe_number(sample.get("total_expected_fields"))
            number_matches = safe_number(cls.get("number_matches"))
            number_discrepancies = safe_number(cls.get("number_discrepancies"))

            if genome_identity is not None:
                identities.append(genome_identity)
            if total_discrepancies is not None:
                discrepancies.append(total_discrepancies)
            if (
                filled_fields is not None
                and total_expected_fields is not None
                and total_expected_fields > 0
            ):
                metadata_completeness.append(100.0 * filled_fields / total_expected_fields)
            if (
                number_matches is not None
                and number_discrepancies is not None
                and (number_matches + number_discrepancies) > 0
            ):
                exact_classification.append(
                    100.0 * number_matches / (number_matches + number_discrepancies)
                )

            clade_match = cls.get("clade_match")
            if clade_match is not None:
                clade_total += 1
                if clade_match is True:
                    clade_hits += 1

            lineage_match = cls.get("lineage_match")
            if lineage_match is not None:
                lineage_total += 1
                if lineage_match is True:
                    lineage_hits += 1

        if not any([identities, discrepancies, metadata_completeness, exact_classification]):
            continue

        group_n_labs = len({record["lab_id"] for record in records})

        group_data.append({
            "label": format_benchmark_group_label(name, version, database_version, group_n_labs),
            "panels": [
                identities,
                discrepancies,
                metadata_completeness,
                exact_classification,
            ],
            "lineage_hit_pct": pct(lineage_hits, lineage_total),
            "clade_hit_pct": pct(clade_hits, clade_total),
            "n_labs": group_n_labs,
        })

    discrepancy_y_limit = 500.0 if comp_code == "FLU2" else None
    metric_panels = [
        ("A. Genome identity", "Genome identity (%)", 0, None),
        ("B. Discrepancies", "Consensus discrepancies", 1, discrepancy_y_limit),
        ("C. Metadata completeness", "Metadata completeness (%)", 2, None),
        ("D. Exact classification concordance", "Exact classification concordance (%)", 3, None),
    ]

    if not group_data:
        return str(output_path)

    fig, axes = plt.subplots(2, 2, figsize=(max(12, len(group_data) * 1.6), 10))
    axes = axes.flatten()

    for ax, (title, ylabel, panel_idx, y_limit) in zip(axes, metric_panels):
        panel_groups = [
            group for group in group_data
            if group["panels"][panel_idx]
        ]
        if not panel_groups:
            ax.set_visible(False)
            continue

        panel_labels = [group["label"] for group in panel_groups]
        panel_data = [list(group["panels"][panel_idx]) for group in panel_groups]
        panel_display_labels = list(panel_labels)
        plotted_data = [list(values) for values in panel_data]
        outlier_annotations = []

        if y_limit is not None:
            for idx, values in enumerate(panel_data):
                outliers_above_limit = sorted([value for value in values if value > y_limit], reverse=True)
                if not outliers_above_limit:
                    continue

                plotted_values = [value for value in values if value <= y_limit]
                if not plotted_values:
                    continue

                plotted_data[idx] = plotted_values
                outlier_annotations.append((idx + 1, outliers_above_limit[0]))

        use_broken_identity_axis = panel_idx == 0 and comp_code in {"SARS1", "FLU1"}

        if use_broken_identity_axis:
            original_spec = ax.get_subplotspec()
            ax.remove()
            inner_gs = original_spec.subgridspec(2, 1, height_ratios=[5.2, 0.8], hspace=0.09)
            ax_upper = fig.add_subplot(inner_gs[0])
            ax_lower = fig.add_subplot(inner_gs[1], sharex=ax_upper)

            for subax in (ax_upper, ax_lower):
                bp = subax.boxplot(
                    plotted_data,
                    labels=panel_display_labels,
                    showfliers=True,
                    patch_artist=True,
                )
                style_boxplot(bp, [comp_code] * len(panel_labels), ax=subax)
                add_component_boxplot_points(
                    subax,
                    bp,
                    plotted_data,
                    list(range(1, len(panel_labels) + 1)),
                    [comp_code] * len(panel_labels),
                )
                subax.set_xlim(0.5, len(panel_labels) + 0.5)

            ax_lower.set_ylim(0, 5)
            ax_upper.set_ylim(90, 102)
            ax_lower.set_yticks([0, 5])
            ax_upper.set_yticks([90, 95, 100])
            ax_upper.spines["bottom"].set_visible(False)
            ax_lower.spines["top"].set_visible(False)
            ax_upper.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
            ax_lower.tick_params(axis="x", rotation=0, labelsize=8)
            ax_upper.set_title(title)
            ax_upper.set_ylabel("")
            ax_lower.set_ylabel("")

            fig.tight_layout()

            upper_pos = ax_upper.get_position()
            lower_pos = ax_lower.get_position()
            y_label_center = (upper_pos.y1 + lower_pos.y0) / 2.0 + 0.04
            fig.text(
                upper_pos.x0 - 0.03,
                y_label_center,
                ylabel,
                rotation=90,
                va="center",
                ha="center",
            )

            slash_x = (-0.02, 0.02)
            slash_y_offset = 0.1
            slash_gap = 0.8
            upper_transform = mtransforms.blended_transform_factory(ax_upper.transAxes, ax_upper.transData)

            ax_upper.plot(
                slash_x,
                [89.5 - slash_y_offset + slash_gap / 2.0, 90 + slash_y_offset + slash_gap / 2.0],
                transform=upper_transform,
                color="#444444",
                linewidth=1.5,
                solid_capstyle="butt",
                clip_on=False,
            )
            ax_upper.plot(
                slash_x,
                [89.5 - slash_y_offset - slash_gap / 2.0, 90 + slash_y_offset - slash_gap / 2.0],
                transform=upper_transform,
                color="#444444",
                linewidth=1.5,
                solid_capstyle="butt",
                clip_on=False,
            )
            continue

        bp = ax.boxplot(
            plotted_data,
            labels=panel_display_labels,
            showfliers=True,
            patch_artist=True,
        )
        style_boxplot(bp, [comp_code] * len(panel_labels), ax=ax)
        add_component_boxplot_points(
            ax,
            bp,
            plotted_data,
            list(range(1, len(panel_labels) + 1)),
            [comp_code] * len(panel_labels),
        )
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.tick_params(axis="x", rotation=0, labelsize=8)
        if "Reads" in ylabel:
            ax.ticklabel_format(axis="y", style="plain", useOffset=False)

        if ylabel.endswith("(%)"):
            style_percent_boxplot_axis(ax)
        elif y_limit is not None:
            ax.set_ylim(0, y_limit)
            annotate_outlier_caps(
                ax,
                outlier_annotations,
                y_limit,
                COMPONENT_BOX_COLORS.get(comp_code, CBF_COLORS["outlier"]),
            )

        if panel_idx == 1:
            secondary_ax = ax.twinx()
            lineage_color = CBF_COLORS["match"]
            clade_color = CBF_COLORS["box_flu2"]
            x_positions = list(range(1, len(panel_groups) + 1))
            lineage_points = [safe_number(group.get("lineage_hit_pct")) for group in panel_groups]
            clade_points = [safe_number(group.get("clade_hit_pct")) for group in panel_groups]

            secondary_ax.plot(
                x_positions,
                lineage_points,
                linestyle="--",
                linewidth=1.4,
                marker="o",
                markersize=5,
                color=lineage_color,
                label="Lineage/Type accuracy",
                zorder=4,
            )
            secondary_ax.plot(
                x_positions,
                clade_points,
                linestyle="--",
                linewidth=1.4,
                marker="s",
                markersize=5,
                color=clade_color,
                label="Clade accuracy",
                zorder=4,
            )
            style_percent_boxplot_axis(secondary_ax)
            secondary_ax.set_ylabel("Classification accuracy (%)")
            secondary_ax.spines["top"].set_visible(False)
            secondary_ax.spines["left"].set_visible(False)
            secondary_ax.spines["right"].set_visible(True)
            secondary_ax.tick_params(axis="y", colors="#444444")
            secondary_ax.legend(
                loc="upper center",
                bbox_to_anchor=(0.5, -0.34),
                borderaxespad=0.0,
                frameon=False,
                fontsize=8,
                ncol=2,
            )

    fig.suptitle(f"{comp_code} performance metrics by bioinformatics protocol")
    fig.tight_layout(rect=(0, 0.14, 1, 1))
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    return str(output_path)


def collect_benchmark_group_records(
    labs: List[Dict[str, Any]],
    comp_code: str,
    name_field: str,
    version_field: Optional[str] = None,
    db_version_field: Optional[str] = None,
) -> List[tuple[str, List[Dict[str, Any]]]]:
    groups = collect_software_groups(labs, comp_code, name_field, version_field, db_version_field)
    ordered_groups = sorted(
        groups.items(),
        key=lambda item: (
            str(item[0][0]).lower() if item[0][0] is not None else "",
            str(item[0][1]).lower() if item[0][1] is not None else "",
            str(item[0][2]).lower() if len(item[0]) > 2 and item[0][2] is not None else "",
        ),
    )
    return [
        (
            format_benchmark_group_label(
                name,
                version,
                db_version,
                len({record["lab_id"] for record in records}),
            ),
            records,
        )
        for (name, version, db_version), records in ordered_groups
    ]


def make_component_benchmark_metric_boxplots(
    labs: List[Dict[str, Any]],
    comp_code: str,
    figures_dir: str | Path,
    benchmark_key: str,
) -> Optional[str]:
    variant_calling_panels = (
        [
            ("A. Allele frequency patterns", "Number of different allele frequency patterns", lambda s: extract_variant_reporting_mode_pct(s, "high_and_low_freq")),
            ("B. Reported variants with AF >=75%", "Number of variants", lambda s: safe_number(s.get("variants", {}).get("number_of_variants_in_consensus"))),
            ("C. Variants with AF >=75% in VCF", "Number of variants", lambda s: safe_number(s.get("variants", {}).get("number_of_variants_in_consensus_vcf"))),
            ("D. Variants with effect", "Number of variants", lambda s: safe_number(s.get("variants", {}).get("number_of_variants_with_effect"))),
            ("E. Metadata-VCF discrepancies", "Number of discrepancies", lambda s: safe_number(s.get("variants", {}).get("discrepancies_in_reported_variants"))),
            ("F. Total variants in VCF", "Number of variants in VCF", lambda s: safe_number(s.get("variants", {}).get("number_of_variants_in_vcf"))),
        ]
        if comp_code.startswith("FLU")
        else [
            ("A. Allele frequency patterns", "Number of different allele frequency patterns", lambda s: extract_variant_reporting_mode_pct(s, "high_and_low_freq")),
            ("B. Discrepancies in reported variants with AF>=75% in VCF", "Number of discrepancies", lambda s: safe_number(s.get("variants", {}).get("discrepancies_in_reported_variants"))),
            ("C. Discrepancy in reported variants with effect", "Number of discrepancies", lambda s: safe_number(s.get("variants", {}).get("discrepancies_in_reported_variants_effect"))),
            ("D. Successful hits", "Number of successful hits", lambda s: safe_number(s.get("variants", {}).get("successful_hits"))),
            ("E. Total discrepancies", "Number of discrepancies", lambda s: safe_number(s.get("variants", {}).get("total_discrepancies"))),
        ]
    )

    benchmark_configs = {
        "dehosting": {
            "name_field": "dehosting_method_software_name",
            "version_field": "dehosting_method_software_version",
            "output_filename": "dehosting_metric_boxplots_by_pipeline.png",
            "title": "performance metrics by de-hosting software",
            "panels": [
                ("A. Host reads", "Host reads (%)", lambda s: safe_number(s.get("metadata_metrics", {}).get("per_reads_host"))),
            ],
        },
        "preprocessing": {
            "name_field": "preprocessing_software_name",
            "version_field": "preprocessing_software_version",
            "output_filename": "preprocessing_metric_boxplots_by_pipeline.png",
            "title": "performance metrics by pre-processing software",
            "panels": [
                ("A. Reads sequenced", "Reads sequenced", lambda s: safe_number(s.get("metadata_metrics", {}).get("number_of_reads_sequenced"))),
                ("B. Reads passing filters", "Reads passing filters", lambda s: safe_number(s.get("metadata_metrics", {}).get("pass_reads"))),
            ],
        },
        "mapping": {
            "name_field": "mapping_software_name",
            "version_field": "mapping_software_version",
            "output_filename": "mapping_metric_boxplots_by_pipeline.png",
            "title": "performance metrics by mapping software",
            "panels": [
                ("A. Viral reads", "Reads virus (%)", lambda s: safe_number(s.get("metadata_metrics", {}).get("per_reads_virus"))),
            ],
        },
        "assembly": {
            "name_field": "assembly",
            "version_field": "assembly_version",
            "output_filename": "assembly_metric_boxplots_by_pipeline.png",
            "title": "performance metrics by assembly software",
            "panels": [
                ("A. Genome length", "Consensus genome length", lambda s: safe_number(s.get("consensus", {}).get("consensus_genome_length"))),
                ("B. Genome identity", "Genome identity (%)", lambda s: safe_number(s.get("consensus", {}).get("genome_identity_pct"))),
                ("C. Discrepancies", "Consensus discrepancies", lambda s: safe_number(s.get("consensus", {}).get("total_discrepancies"))),
            ],
        },
        "consensus_software": {
            "name_field": "consensus_sequence_software_name",
            "version_field": "consensus_sequence_software_version",
            "output_filename": "consensus_metric_boxplots_by_pipeline.png",
            "title": "performance metrics by consensus software",
            "panels": [
                ("A. Genome length", "Consensus genome length", lambda s: safe_number(s.get("consensus", {}).get("consensus_genome_length"))),
                ("B. Genome identity", "Genome identity (%)", lambda s: safe_number(s.get("consensus", {}).get("genome_identity_pct"))),
                ("C. Discrepancies", "Consensus discrepancies", lambda s: safe_number(s.get("consensus", {}).get("total_discrepancies"))),
            ],
        },
        "variant_calling": {
            "name_field": "variant_calling_software_name",
            "version_field": "variant_calling_software_version",
            "output_filename": "variant_calling_metric_boxplots_by_pipeline.png",
            "title": "performance metrics by variant calling software",
            "panels": variant_calling_panels,
        },
        "clade_assignment": {
            "name_field": "clade_assignment_software_name",
            "version_field": "clade_assignment_software_version",
            "db_version_field": "clade_assignment_software_database_version",
            "output_filename": "clade_assignment_metric_boxplots_by_pipeline.png",
            "title": "performance metrics by clade assignment software",
            "panels": [
                ("A. Clade concordance", "Clade concordance (%)", lambda s: 100.0 if s.get("classification", {}).get("clade_match") is True else (0.0 if s.get("classification", {}).get("clade_match") is False else None)),
                ("B. Clade discrepancy", "Clade discrepancy (%)", lambda s: 0.0 if s.get("classification", {}).get("clade_match") is True else (100.0 if s.get("classification", {}).get("clade_match") is False else None)),
            ],
        },
        "lineage_assignment": {
            "name_field": "lineage_assignment_software_name",
            "version_field": "lineage_assignment_software_version",
            "db_version_field": "lineage_assignment_database_version",
            "output_filename": "lineage_assignment_metric_boxplots_by_pipeline.png",
            "title": "performance metrics by lineage assignment software",
            "panels": [
                ("A. Lineage concordance", "Lineage concordance (%)", lambda s: 100.0 if s.get("classification", {}).get("lineage_match") is True else (0.0 if s.get("classification", {}).get("lineage_match") is False else None)),
                ("B. Lineage discrepancy", "Lineage discrepancy (%)", lambda s: 0.0 if s.get("classification", {}).get("lineage_match") is True else (100.0 if s.get("classification", {}).get("lineage_match") is False else None)),
            ],
        },
        "type_assignment": {
            "name_field": "type_assignment_software_name",
            "version_field": None,
            "db_version_field": "type_assignment_software_database_version",
            "output_filename": "type_assignment_metric_boxplots_by_pipeline.png",
            "title": "performance metrics by type assignment software",
            "panels": [
                ("A. Type concordance", "Type concordance (%)", lambda s: 100.0 if s.get("classification", {}).get("lineage_match") is True else (0.0 if s.get("classification", {}).get("lineage_match") is False else None)),
                ("B. Type discrepancy", "Type discrepancy (%)", lambda s: 0.0 if s.get("classification", {}).get("lineage_match") is True else (100.0 if s.get("classification", {}).get("lineage_match") is False else None)),
            ],
        },
        "subtype_assignment": {
            "name_field": "subtype_assignment_software_name",
            "version_field": "subtype_assignment_software_version",
            "db_version_field": "subtype_assignment_software_database_version",
            "output_filename": "subtype_assignment_metric_boxplots_by_pipeline.png",
            "title": "performance metrics by subtype assignment software",
            "panels": [
                ("A. Subtype concordance", "Subtype concordance (%)", lambda s: 100.0 if s.get("classification", {}).get("lineage_match") is True else (0.0 if s.get("classification", {}).get("lineage_match") is False else None)),
                ("B. Subtype discrepancy", "Subtype discrepancy (%)", lambda s: 0.0 if s.get("classification", {}).get("lineage_match") is True else (100.0 if s.get("classification", {}).get("lineage_match") is False else None)),
            ],
        },
    }

    config = benchmark_configs.get(benchmark_key)
    if not config:
        return None

    output_dir = ensure_component_figures_dir(figures_dir, comp_code)
    output_path = output_dir / config["output_filename"]
    groups = collect_benchmark_group_records(
        labs,
        comp_code,
        config["name_field"],
        config["version_field"],
        config.get("db_version_field"),
    )
    if not groups:
        return str(output_path)

    group_metric_values = []
    for label, records in groups:
        per_panel_values = []
        has_any = False
        for _, _, extractor in config["panels"]:
            values = []
            for record in records:
                value = extractor(record["sample"])
                if value is not None:
                    values.append(value)
            if values:
                has_any = True
            per_panel_values.append(values)
        if has_any:
            group_metric_values.append({
                "label": label,
                "panels": per_panel_values,
                "records": records,
            })

    metric_panels = []
    for panel_idx, (panel_title, ylabel, _) in enumerate(config["panels"]):
        panel_groups = [
            group for group in group_metric_values
            if group["panels"][panel_idx]
        ]
        if panel_groups:
            metric_panels.append((panel_title, ylabel, panel_idx, panel_groups))

    if not group_metric_values or not metric_panels:
        return str(output_path)

    n_panels = len(metric_panels)
    stacked_vertical_benchmark = benchmark_key in {"lineage_assignment", "type_assignment", "clade_assignment"} and n_panels == 2
    ncols = 1 if n_panels == 1 or stacked_vertical_benchmark else 2
    nrows = int(np.ceil(n_panels / ncols))
    max_groups_in_panel = max(len(panel_groups) for _, _, _, panel_groups in metric_panels)
    if stacked_vertical_benchmark:
        fig_width = max(18, max_groups_in_panel * 2.3)
        fig_height = 5.2 * nrows
    else:
        fig_width = max(14, max_groups_in_panel * 2.1)
        fig_height = 4.8 * nrows
    fig, axes = plt.subplots(nrows, ncols, figsize=(fig_width, fig_height))
    axes = np.atleast_1d(axes).flatten()

    for ax, (title, ylabel, panel_idx, panel_groups) in zip(axes, metric_panels):
        panel_labels = [group["label"] for group in panel_groups]
        panel_data = [group["panels"][panel_idx] for group in panel_groups]

        if benchmark_key == "variant_calling" and panel_idx == 0:
            mode_keys = ["high_and_low_freq", "low_freq_only", "high_freq_only"]
            mode_labels = ["High and low frequency", "Low frequency only", "High frequency only"]
            mode_colors = [
                CBF_COLORS["high_and_low_freq"],
                CBF_COLORS["low_freq_only"],
                CBF_COLORS["high_freq_only"],
            ]
            x_positions = np.arange(len(panel_labels))
            bottoms = np.zeros(len(panel_labels))

            for mode_key, mode_label, mode_color in zip(mode_keys, mode_labels, mode_colors):
                counts = [
                    sum(
                        1
                        for record in group["records"]
                        if record["sample"].get("variants", {}).get(mode_key) is True
                    )
                    for group in panel_groups
                ]
                ax.bar(
                    x_positions,
                    counts,
                    bottom=bottoms,
                    color=mode_color,
                    edgecolor="#4A4A4A",
                    linewidth=0.8,
                    label=mode_label,
                )
                bottoms = bottoms + np.array(counts)

            ax.set_title(title)
            ax.set_ylabel("Samples (n)")
            ax.set_xticks(x_positions)
            ax.set_xticklabels(panel_labels, fontsize=8)
            ax.legend(
                loc="upper left",
                bbox_to_anchor=(1.02, 1.0),
                borderaxespad=0.0,
                frameon=False,
                fontsize=8,
            )
            continue

        bp = ax.boxplot(
            panel_data,
            labels=panel_labels,
            showfliers=True,
            patch_artist=True,
        )
        style_boxplot(bp, [comp_code] * len(panel_labels), ax=ax)
        add_component_boxplot_points(
            ax,
            bp,
            panel_data,
            list(range(1, len(panel_labels) + 1)),
            [comp_code] * len(panel_labels),
        )
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.tick_params(
            axis="x",
            rotation=0,
            labelsize=8,
            pad=10 if stacked_vertical_benchmark else 4,
        )
        if ylabel.endswith("(%)"):
            style_percent_boxplot_axis(ax)
        if "Reads" in ylabel:
            ax.ticklabel_format(axis="y", style="plain", useOffset=False)

    for ax in axes[len(metric_panels):]:
        ax.set_visible(False)

    fig.suptitle(f"{comp_code} {config['title']}")
    fig.tight_layout()
    if benchmark_key == "variant_calling" and ncols > 1:
        fig.subplots_adjust(wspace=0.38, hspace=0.35, top=0.92)
    if stacked_vertical_benchmark:
        fig.subplots_adjust(hspace=0.45, top=0.92, bottom=0.1)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

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


def collect_consensus_discrepancies_by_sample(
    labs: List[Dict[str, Any]],
    comp_code: str,
) -> Dict[str, List[float]]:
    discrepancies_by_sample: Dict[str, List[float]] = defaultdict(list)

    for lab in labs:
        comp = lab.get("components", {}).get(comp_code)
        if not comp:
            continue

        for sample_id, sample in comp.get("samples", {}).items():
            consensus = sample.get("consensus", {})
            genome_identity_pct = safe_number(consensus.get("genome_identity_pct"))
            total_discrepancies = safe_number(consensus.get("total_discrepancies"))

            if genome_identity_pct is None or total_discrepancies is None:
                continue

            discrepancies_by_sample[sample_id].append(total_discrepancies)

    return discrepancies_by_sample


def make_component_consensus_discrepancies_boxplot_by_sample(
    general_data: Dict[str, Any],
    labs: List[Dict[str, Any]],
    comp_code: str,
    figures_dir: str | Path,
    output_filename: str = "consensus_discrepancies_boxplot_by_sample.png",
) -> str:
    output_dir = ensure_component_figures_dir(figures_dir, comp_code)
    output_path = output_dir / output_filename

    comp_data = general_data.get("components", {}).get(comp_code, {})
    sample_order = [
        sample.get("collecting_lab_sample_id")
        for sample in comp_data.get("consensus", {}).get("samples", [])
        if sample.get("collecting_lab_sample_id")
    ]
    discrepancies_by_sample = collect_consensus_discrepancies_by_sample(labs, comp_code)

    sample_names = [sample_id for sample_id in sample_order if discrepancies_by_sample.get(sample_id)]
    data = [discrepancies_by_sample[sample_id] for sample_id in sample_names]

    if not data:
        return str(output_path)

    y_limit = COMPONENT_CONSENSUS_SAMPLE_Y_LIMITS.get(comp_code)
    plotted_data = [list(values) for values in data]
    outlier_annotations = []

    if y_limit is not None:
        for idx, values in enumerate(data):
            outliers_above_limit = sorted([value for value in values if value > y_limit], reverse=True)
            if not outliers_above_limit:
                continue

            excluded_outlier = outliers_above_limit[0]
            removed = False
            filtered_values = []
            for value in values:
                if not removed and value == excluded_outlier:
                    removed = True
                    continue
                filtered_values.append(value)

            if removed and filtered_values:
                plotted_data[idx] = filtered_values
                outlier_annotations.append((idx + 1, excluded_outlier))

    plotted_max = max(max(values) for values in plotted_data if values)
    y_upper = y_limit if y_limit is not None else (plotted_max * 1.15 if plotted_max > 0 else 1.0)

    plt.figure(figsize=(max(8, len(sample_names) * 1.15), 6))
    bp = plt.boxplot(
        plotted_data,
        labels=sample_names,
        showfliers=True,
        patch_artist=True,
    )
    style_boxplot(bp, [comp_code] * len(sample_names), ax=plt.gca())
    add_component_boxplot_points(
        plt.gca(),
        bp,
        plotted_data,
        list(range(1, len(sample_names) + 1)),
        [comp_code] * len(sample_names),
    )
    plt.xlabel("Sample")
    plt.ylabel("Consensus discrepancies")
    plt.title(f"{comp_code} consensus discrepancies by sample")
    plt.xticks(rotation=0, ha="center")

    plt.ylim(0, y_upper)

    if outlier_annotations:
        outlier_color = COMPONENT_BOX_COLORS.get(comp_code, CBF_COLORS["outlier"])
        for x_pos, display_value in outlier_annotations:
            y_marker = y_upper * 0.95
            y_text = y_upper * 0.91
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


def make_component_consensus_discrepancies_stacked_by_sample(
    general_data: Dict[str, Any],
    labs: List[Dict[str, Any]],
    comp_code: str,
    figures_dir: str | Path,
    output_filename: str = "consensus_discrepancies_stacked_by_sample.png",
) -> str:
    output_dir = ensure_component_figures_dir(figures_dir, comp_code)
    output_path = output_dir / output_filename

    comp_data = general_data.get("components", {}).get(comp_code, {})
    sample_names = [
        sample.get("collecting_lab_sample_id")
        for sample in comp_data.get("consensus", {}).get("samples", [])
        if sample.get("collecting_lab_sample_id")
    ]
    if not sample_names:
        return str(output_path)

    stacked_values = {key: [] for key in CONSENSUS_DISCREPANCY_TYPE_ORDER}
    for sample_id in sample_names:
        totals = {key: 0.0 for key in CONSENSUS_DISCREPANCY_TYPE_ORDER}
        for lab in labs:
            comp = lab.get("components", {}).get(comp_code)
            if not comp:
                continue

            sample = comp.get("samples", {}).get(sample_id)
            if not sample:
                continue

            consensus = sample.get("consensus", {})
            gi = safe_number(consensus.get("genome_identity_pct"))
            if gi is None:
                continue

            discrepancy_breakdown = consensus.get("discrepancy_breakdown", {})
            for key in CONSENSUS_DISCREPANCY_TYPE_ORDER:
                value = safe_number(discrepancy_breakdown.get(key))
                if value is not None:
                    totals[key] += value

        for key in CONSENSUS_DISCREPANCY_TYPE_ORDER:
            stacked_values[key].append(totals[key])

    x_positions = list(range(len(sample_names)))
    bottoms = [0.0] * len(sample_names)

    plt.figure(figsize=(max(8, len(sample_names) * 1.15), 6))

    for key in CONSENSUS_DISCREPANCY_TYPE_ORDER:
        values = stacked_values[key]
        plt.bar(
            x_positions,
            values,
            bottom=bottoms,
            label=CONSENSUS_DISCREPANCY_TYPE_LABELS.get(key, key),
            color=CONSENSUS_DISCREPANCY_TYPE_COLORS[key],
        )
        bottoms = [bottom + value for bottom, value in zip(bottoms, values)]

    plt.xticks(x_positions, sample_names, rotation=0, ha="center")
    plt.xlabel("Sample")
    plt.ylabel("Total consensus discrepancies")
    plt.title(f"{comp_code} consensus discrepancy types by sample")
    plt.legend(
        title="Discrepancy type",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=False,
    )
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()

    return str(output_path)


def make_component_consensus_discrepancy_type_boxplot(
    general_data: Dict[str, Any],
    labs: List[Dict[str, Any]],
    comp_code: str,
    figures_dir: str | Path,
    output_filename: str = "consensus_discrepancy_type_boxplot.png",
) -> str:
    output_dir = ensure_component_figures_dir(figures_dir, comp_code)
    output_path = output_dir / output_filename

    labels = []
    data = []
    used_keys = []
    outlier_annotations = []
    y_limit = COMPONENT_CONSENSUS_TYPE_BOXPLOT_Y_LIMITS.get(comp_code)
    for key in CONSENSUS_DISCREPANCY_TYPE_ORDER:
        values = []
        for lab in labs:
            comp = lab.get("components", {}).get(comp_code)
            if not comp:
                continue

            for sample in comp.get("samples", {}).values():
                consensus = sample.get("consensus", {})
                gi = safe_number(consensus.get("genome_identity_pct"))
                if gi is None:
                    continue

                discrepancy_breakdown = consensus.get("discrepancy_breakdown", {})
                value = safe_number(discrepancy_breakdown.get(key))
                if value is not None:
                    values.append(value)

        if values:
            plotted_values = list(values)
            if y_limit is not None:
                outliers_above_limit = sorted([value for value in plotted_values if value > y_limit], reverse=True)
                if outliers_above_limit:
                    plotted_values = [value for value in plotted_values if value <= y_limit]
                    if plotted_values:
                        outlier_annotations.append((len(labels) + 1, outliers_above_limit))

            if not plotted_values:
                continue

            used_keys.append(key)
            labels.append(CONSENSUS_DISCREPANCY_TYPE_LABELS.get(key, key))
            data.append(plotted_values)

    if not data:
        return str(output_path)

    plt.figure(figsize=(max(9, len(labels) * 1.35), 6))
    bp = plt.boxplot(
        data,
        labels=labels,
        showfliers=True,
        patch_artist=True,
    )

    for patch, key in zip(bp["boxes"], used_keys):
        patch.set_facecolor(CONSENSUS_DISCREPANCY_TYPE_COLORS.get(key, CBF_COLORS["box_default"]))
        patch.set_edgecolor("#333333")
        patch.set_alpha(0.75)

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
    style_boxplot_axes(plt.gca())
    add_colored_boxplot_points(
        plt.gca(),
        bp,
        data,
        list(range(1, len(labels) + 1)),
        [CONSENSUS_DISCREPANCY_TYPE_COLORS.get(key, CBF_COLORS["box_default"]) for key in used_keys],
    )

    plt.xticks(rotation=20, ha="center")
    plt.xlabel("Discrepancy type")
    plt.ylabel("Number of discrepancies per sample")
    plt.title(f"{comp_code} consensus discrepancy types")
    if y_limit is not None:
        plt.ylim(0, y_limit)
        for x_pos, outliers in outlier_annotations:
            outlier_color = CONSENSUS_DISCREPANCY_TYPE_COLORS.get(used_keys[x_pos - 1], CBF_COLORS["outlier"])
            plt.text(
                x_pos,
                y_limit * 0.95,
                "*",
                ha="center",
                va="center",
                fontsize=18,
                color=outlier_color,
                fontweight="bold",
            )
            plt.text(
                x_pos,
                y_limit * 0.91,
                f"Outlier: {outliers[0]:g}",
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


def make_component_variant_discrepancies_stacked_by_sample(
    general_data: Dict[str, Any],
    labs: List[Dict[str, Any]],
    comp_code: str,
    figures_dir: str | Path,
    output_filename: str = "variant_discrepancies_stacked_by_sample.png",
) -> str:
    output_dir = ensure_component_figures_dir(figures_dir, comp_code)
    output_path = output_dir / output_filename

    comp_data = general_data.get("components", {}).get(comp_code, {})
    sample_names = [
        sample.get("collecting_lab_sample_id")
        for sample in comp_data.get("variant", {}).get("samples", [])
        if sample.get("collecting_lab_sample_id")
    ]
    if not sample_names:
        return str(output_path)

    stacked_values = {key: [] for key in VARIANT_DISCREPANCY_TYPE_ORDER}
    for sample_id in sample_names:
        totals = {key: 0.0 for key in VARIANT_DISCREPANCY_TYPE_ORDER}
        for lab in labs:
            comp = lab.get("components", {}).get(comp_code)
            if not comp:
                continue

            sample = comp.get("samples", {}).get(sample_id)
            if not sample:
                continue

            variants = sample.get("variants", {})
            total_discrepancies = safe_number(variants.get("total_discrepancies"))
            if total_discrepancies is None:
                continue

            for key in VARIANT_DISCREPANCY_TYPE_ORDER:
                value = safe_number(variants.get(key))
                if value is not None:
                    totals[key] += value

        for key in VARIANT_DISCREPANCY_TYPE_ORDER:
            stacked_values[key].append(totals[key])

    x_positions = list(range(len(sample_names)))
    bottoms = [0.0] * len(sample_names)

    plt.figure(figsize=(max(8, len(sample_names) * 1.15), 6))
    for key in VARIANT_DISCREPANCY_TYPE_ORDER:
        values = stacked_values[key]
        plt.bar(
            x_positions,
            values,
            bottom=bottoms,
            label=VARIANT_DISCREPANCY_TYPE_LABELS.get(key, key),
            color=VARIANT_DISCREPANCY_TYPE_COLORS[key],
        )
        bottoms = [bottom + value for bottom, value in zip(bottoms, values)]

    plt.legend(
        title="Discrepancy type",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=False,
    )

    plt.xticks(x_positions, sample_names, rotation=0, ha="center")
    plt.xlabel("Sample")
    plt.ylabel("Total variant discrepancies")
    plt.title(f"{comp_code} variant discrepancy types by sample")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()

    return str(output_path)


def make_component_variant_discrepancy_type_boxplot(
    labs: List[Dict[str, Any]],
    comp_code: str,
    figures_dir: str | Path,
    output_filename: str = "variant_discrepancy_type_boxplot.png",
) -> str:
    output_dir = ensure_component_figures_dir(figures_dir, comp_code)
    output_path = output_dir / output_filename

    labels = []
    data = []
    used_keys = []

    for key in VARIANT_DISCREPANCY_TYPE_ORDER:
        values = []
        for lab in labs:
            comp = lab.get("components", {}).get(comp_code)
            if not comp:
                continue

            for sample in comp.get("samples", {}).values():
                variants = sample.get("variants", {})
                total_discrepancies = safe_number(variants.get("total_discrepancies"))
                if total_discrepancies is None:
                    continue

                value = safe_number(variants.get(key))
                if value is not None:
                    values.append(value)

        if values:
            used_keys.append(key)
            labels.append(VARIANT_DISCREPANCY_TYPE_LABELS.get(key, key))
            data.append(values)

    if not data:
        return str(output_path)

    broken_axis_ranges = {
        "SARS1": ((0, 80), (450, None)),
        "SARS2": ((0, 40), (75, None)),
    }
    broken_axis_spec = broken_axis_ranges.get(comp_code)

    def style_variant_bp(ax, bp_obj):
        for patch, key in zip(bp_obj["boxes"], used_keys):
            patch.set_facecolor(VARIANT_DISCREPANCY_TYPE_COLORS.get(key, CBF_COLORS["box_default"]))
            patch.set_edgecolor("#333333")
            patch.set_alpha(0.75)

        for median_line in bp_obj["medians"]:
            median_line.set_color(CBF_COLORS["median"])
            median_line.set_linewidth(2)

        for whisker in bp_obj["whiskers"]:
            whisker.set_color("#444444")
        for cap in bp_obj["caps"]:
            cap.set_color("#444444")
        for flier in bp_obj["fliers"]:
            flier.set_marker("o")
            flier.set_markerfacecolor("white")
            flier.set_markeredgecolor("#444444")
            flier.set_markersize(5)
        style_boxplot_axes(ax)
        add_colored_boxplot_points(
            ax,
            bp_obj,
            data,
            list(range(1, len(labels) + 1)),
            [VARIANT_DISCREPANCY_TYPE_COLORS.get(key, CBF_COLORS["box_default"]) for key in used_keys],
        )

    if broken_axis_spec is not None:
        (lower_min, lower_max), (upper_min, upper_max) = broken_axis_spec
        max_value = max(max(values) for values in data if values)
        if upper_max is None:
            upper_max = max_value * 1.08 if max_value > 0 else upper_min + 1

        fig, (ax_upper, ax_lower) = plt.subplots(
            2,
            1,
            sharex=True,
            figsize=(max(9, len(labels) * 1.35), 6.8),
            gridspec_kw={"height_ratios": [1.2, 3.5], "hspace": 0.05},
        )

        for ax in (ax_upper, ax_lower):
            bp = ax.boxplot(
                data,
                labels=labels,
                showfliers=True,
                patch_artist=True,
            )
            style_variant_bp(ax, bp)
            ax.set_xlim(0.5, len(labels) + 0.5)

        ax_lower.set_ylim(lower_min, lower_max)
        ax_upper.set_ylim(upper_min, upper_max)
        ax_upper.spines["bottom"].set_visible(False)
        ax_lower.spines["top"].set_visible(False)
        ax_upper.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
        ax_lower.tick_params(axis="x", rotation=0, labelsize=10)
        for tick in ax_lower.get_xticklabels():
            tick.set_ha("center")
        ax_lower.set_xlabel("Discrepancy type")
        ax_lower.set_ylabel("Number of discrepancies")
        ax_upper.set_title(f"{comp_code} variant discrepancy types")

        upper_pos = ax_upper.get_position()
        lower_pos = ax_lower.get_position()
        x_center = upper_pos.x0
        y_center = (upper_pos.y0 + lower_pos.y1) / 2.0
        dx = 0.012
        dy = 0.01
        gap = 0.012

        fig.add_artist(plt.Line2D(
            [x_center - dx, x_center + dx],
            [y_center - dy + gap / 2.0, y_center + dy + gap / 2.0],
            transform=fig.transFigure,
            color="#444444",
            linewidth=1.2,
            solid_capstyle="butt",
            clip_on=False,
        ))
        fig.add_artist(plt.Line2D(
            [x_center - dx, x_center + dx],
            [y_center - dy - gap / 2.0, y_center + dy - gap / 2.0],
            transform=fig.transFigure,
            color="#444444",
            linewidth=1.2,
            solid_capstyle="butt",
            clip_on=False,
        ))

        fig.tight_layout()
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.figure(figsize=(max(9, len(labels) * 1.35), 6))
        bp = plt.boxplot(
            data,
            labels=labels,
            showfliers=True,
            patch_artist=True,
        )

        style_variant_bp(plt.gca(), bp)

        plt.xticks(rotation=0, ha="center")
        plt.xlabel("Discrepancy type")
        plt.ylabel("Number of discrepancies")
        plt.title(f"{comp_code} variant discrepancy types")
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close()

    return str(output_path)


INFLUENZA_REPORTING_METRICS = [
    ("number_of_variants_in_consensus", ">75% AF in metadata", "#0072B2"),
    ("number_of_variants_in_consensus_vcf", ">75% AF in VCF files", "#009E73"),
    ("discrepancies_in_reported_variants", "Metadata-VCF discrepancies", "#E69F00"),
]

INFLUENZA_TOTAL_VCF_COLOR = "#CC79A7"


def collect_influenza_variant_reporting_by_sample(
    labs: List[Dict[str, Any]],
    comp_code: str,
) -> Dict[str, Dict[str, List[float]]]:
    metrics_by_sample: Dict[str, Dict[str, List[float]]] = defaultdict(
        lambda: defaultdict(list)
    )

    for lab in labs:
        comp = lab.get("components", {}).get(comp_code)
        if not comp:
            continue

        for sample_id, sample in comp.get("samples", {}).items():
            variants = sample.get("variants", {})
            for metric_key, _, _ in INFLUENZA_REPORTING_METRICS:
                value = safe_number(variants.get(metric_key))
                if value is not None:
                    metrics_by_sample[sample_id][metric_key].append(value)

            total_vcf = safe_number(variants.get("number_of_variants_in_vcf"))
            if total_vcf is not None:
                metrics_by_sample[sample_id]["number_of_variants_in_vcf"].append(total_vcf)

    return metrics_by_sample


def style_boxplot_with_color(bp: Dict[str, Any], color: str, ax: Optional[Any] = None) -> None:
    for patch in bp["boxes"]:
        patch.set_facecolor(color)
        patch.set_edgecolor("#333333")
        patch.set_alpha(0.75)

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
    style_boxplot_axes(ax if ax is not None else plt.gca())


def trim_boxplot_extreme_outliers(
    data: List[List[float]],
) -> tuple[List[List[float]], List[tuple[int, float]]]:
    trimmed_data: List[List[float]] = []
    outlier_annotations: List[tuple[int, float]] = []

    for idx, values in enumerate(data, start=1):
        plotted_values = list(values)
        if len(values) >= 4:
            stats = boxplot_stats(values)[0]
            whisker_high = safe_number(stats.get("whishi"))
            max_value = max(values)
            if (
                whisker_high is not None
                and whisker_high > 0
                and max_value > whisker_high * 2
            ):
                candidate_values = [value for value in values if value <= whisker_high]
                if candidate_values:
                    plotted_values = candidate_values
                    outlier_annotations.append((idx, max_value))

        trimmed_data.append(plotted_values)

    return trimmed_data, outlier_annotations


def add_boxplot_points(
    ax: Any,
    bp: Dict[str, Any],
    data: List[List[float]],
    positions: List[float],
    color: str,
) -> None:
    jittered_positions = build_jittered_positions(data, positions)
    for xs, values in zip(jittered_positions, data):
        if not values:
            continue
        ax.scatter(
            xs,
            values,
            s=16,
            color=color,
            alpha=0.45,
            edgecolors="none",
            zorder=3,
        )
    align_boxplot_fliers(bp, data, jittered_positions)


def annotate_outlier_caps(
    ax: Any,
    annotations: List[tuple[float, float]],
    y_upper: float,
    color: str,
) -> None:
    for x_pos, display_value in annotations:
        ax.text(
            x_pos,
            y_upper * 0.95,
            "*",
            ha="center",
            va="center",
            fontsize=18,
            color=color,
            fontweight="bold",
        )
        ax.text(
            x_pos,
            y_upper * 0.91,
            f"Outlier:\n{display_value:g}",
            ha="center",
            va="top",
            fontsize=9,
            color=color,
            fontweight="bold",
        )


def make_component_influenza_variant_reporting_summary(
    general_data: Dict[str, Any],
    labs: List[Dict[str, Any]],
    comp_code: str,
    figures_dir: str | Path,
    output_filename: str = "influenza_variant_reporting_summary_by_sample.png",
) -> str:
    output_dir = ensure_component_figures_dir(figures_dir, comp_code)
    output_path = output_dir / output_filename

    comp_data = general_data.get("components", {}).get(comp_code, {})
    samples = comp_data.get("variant", {}).get("samples", [])
    if not samples:
        return str(output_path)

    metrics_by_sample = collect_influenza_variant_reporting_by_sample(labs, comp_code)
    sample_names = [
        sample.get("collecting_lab_sample_id")
        for sample in samples
        if sample.get("collecting_lab_sample_id")
        and any(metrics_by_sample.get(sample.get("collecting_lab_sample_id"), {}).get(metric_key, []) for metric_key in [
            "number_of_variants_in_consensus",
            "number_of_variants_in_consensus_vcf",
            "discrepancies_in_reported_variants",
            "number_of_variants_in_vcf",
        ])
    ]
    if not sample_names:
        return str(output_path)

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(max(12, len(sample_names) * 1.8), 6.5),
        gridspec_kw={"width_ratios": [2.6, 1.2]},
    )

    base_positions = np.arange(1, len(sample_names) + 1)
    offsets = [-0.27, 0.0, 0.27]
    box_width = 0.22

    panel_a_max = 0.0
    panel_a_annotations = []

    for offset, (metric_key, label, color) in zip(offsets, INFLUENZA_REPORTING_METRICS):
        metric_data = [
            list(metrics_by_sample.get(sample_id, {}).get(metric_key, []))
            for sample_id in sample_names
        ]
        positions = [pos + offset for pos in base_positions]

        trimmed_data, outlier_annotations = trim_boxplot_extreme_outliers(metric_data)
        valid_data = [values for values in trimmed_data if values]
        if not valid_data:
            continue

        bp = axes[0].boxplot(
            trimmed_data,
            positions=positions,
            widths=box_width,
            showfliers=True,
            patch_artist=True,
        )
        style_boxplot_with_color(bp, color, ax=axes[0])
        add_boxplot_points(axes[0], bp, trimmed_data, positions, color)

        panel_a_max = max(panel_a_max, max(max(values) for values in valid_data))
        panel_a_annotations.extend((positions[idx - 1], value, color) for idx, value in outlier_annotations)

    panel_a_upper = panel_a_max * 1.18 if panel_a_max > 0 else 1.0
    axes[0].set_ylim(0, panel_a_upper)
    for x_pos, display_value, color in panel_a_annotations:
        annotate_outlier_caps(axes[0], [(x_pos, display_value)], panel_a_upper, color)

    axes[0].set_xticks(base_positions)
    axes[0].set_xticklabels(sample_names, rotation=45, ha="right")
    axes[0].set_xlim(0.5, len(sample_names) + 0.5)
    axes[0].set_ylabel("Number of variants")
    axes[0].set_title("A. High-frequency variant reporting")
    axes[0].legend(
        handles=[
            plt.Line2D([0], [0], color=color, lw=8, alpha=0.75)
            for _, _, color in INFLUENZA_REPORTING_METRICS
        ],
        labels=[label for _, label, _ in INFLUENZA_REPORTING_METRICS],
        frameon=False,
        loc="upper right",
    )

    total_vcf_data = [
        list(metrics_by_sample.get(sample_id, {}).get("number_of_variants_in_vcf", []))
        for sample_id in sample_names
    ]
    trimmed_total_vcf_data, total_vcf_outliers = trim_boxplot_extreme_outliers(total_vcf_data)
    valid_total_vcf_data = [values for values in trimmed_total_vcf_data if values]
    if valid_total_vcf_data:
        bp_total = axes[1].boxplot(
            trimmed_total_vcf_data,
            positions=base_positions,
            widths=0.5,
            showfliers=True,
            patch_artist=True,
        )
        style_boxplot_with_color(bp_total, INFLUENZA_TOTAL_VCF_COLOR, ax=axes[1])
        add_boxplot_points(
            axes[1],
            bp_total,
            trimmed_total_vcf_data,
            list(base_positions),
            INFLUENZA_TOTAL_VCF_COLOR,
        )

        panel_b_max = max(max(values) for values in valid_total_vcf_data)
        panel_b_upper = panel_b_max * 1.18 if panel_b_max > 0 else 1.0
        axes[1].set_ylim(0, panel_b_upper)
        annotate_outlier_caps(
            axes[1],
            [(base_positions[idx - 1], value) for idx, value in total_vcf_outliers],
            panel_b_upper,
            INFLUENZA_TOTAL_VCF_COLOR,
        )

    axes[1].set_xticks(base_positions)
    axes[1].set_xticklabels(sample_names, rotation=45, ha="right")
    axes[1].set_xlim(0.5, len(sample_names) + 0.5)
    axes[1].set_ylabel("Total variants in VCF")
    axes[1].set_title("B. Total variants in VCF")

    fig.suptitle(f"{comp_code} influenza variant reporting summary by sample")
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    return str(output_path)


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
    style_boxplot(bp, component_names, ax=plt.gca())
    add_component_boxplot_points(
        plt.gca(),
        bp,
        plotted_data,
        list(range(1, len(component_names) + 1)),
        component_names,
    )
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
    style_boxplot(bp, component_names, ax=plt.gca())
    add_component_boxplot_points(
        plt.gca(),
        bp,
        plotted_data,
        list(range(1, len(component_names) + 1)),
        component_names,
    )
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
        ("High and low frequency", safe_number(sars_variants.get("high_and_low_freq_pct"))),
        ("Low frequency only", safe_number(sars_variants.get("low_freq_only_pct"))),
        ("High frequency only", safe_number(sars_variants.get("high_freq_only_pct"))),
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
        ("High and low frequency", safe_number(influenza_variants.get("high_and_low_freq_pct"))),
        ("Low frequency only", safe_number(influenza_variants.get("low_freq_only_pct"))),
        ("High frequency only", safe_number(influenza_variants.get("high_freq_only_pct"))),
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
    match_rates = []
    discrepancy_rates = []

    for comp_code, comp_data in components.items():
        matches, discrepancies = qc_hits_discrepancies_from_component(comp_data)
        total = matches + discrepancies
        match_rate = 100.0 * matches / total if total else 0.0
        discrepancy_rate = 100.0 * discrepancies / total if total else 0.0

        component_names.append(comp_code)
        match_rates.append(match_rate)
        discrepancy_rates.append(discrepancy_rate)

    output_dir = ensure_network_figures_dir(figures_dir)
    output_path = output_dir / output_filename

    plt.figure(figsize=(10, 6))
    x_positions = list(range(len(component_names)))

    match_bars = plt.bar(x_positions, match_rates, label="Match", color=CBF_COLORS["match"])
    discrepancy_bars = plt.bar(x_positions, discrepancy_rates, bottom=match_rates, label="Discrepancy", color=CBF_COLORS["discrepancy"])

    plt.xticks(x_positions, component_names)
    plt.xlabel("Component")
    plt.ylabel("QC evaluations (%)")
    plt.title(title)
    plt.ylim(0, 100)
    plt.legend(frameon=False, loc="lower center", bbox_to_anchor=(0.5, -0.22), ncol=2)

    for bar, value in zip(match_bars, match_rates):
        if value <= 0:
            continue
        plt.text(
            bar.get_x() + bar.get_width() / 2,
            value / 2,
            f"{value:.1f}%",
            ha="center",
            va="center",
            fontsize=8,
            color="white",
            fontweight="bold",
        )

    for bar, match_value, discrepancy_value in zip(discrepancy_bars, match_rates, discrepancy_rates):
        if discrepancy_value <= 0:
            continue
        plt.text(
            bar.get_x() + bar.get_width() / 2,
            match_value + discrepancy_value / 2,
            f"{discrepancy_value:.1f}%",
            ha="center",
            va="center",
            fontsize=8,
            color="white",
            fontweight="bold",
        )

    plt.tight_layout(rect=(0, 0.08, 1, 1))
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

    outputs["classification_summary"] = make_combined_classification_summary_plot(
        general_data=general_data,
        labs=labs,
        figures_dir=figures_dir,
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


def generate_component_figures(
    general_data: Dict[str, Any],
    labs: List[Dict[str, Any]],
    figures_dir: str | Path = "figures",
) -> None:
    for comp_code in general_data.get("components", {}).keys():
        make_component_consensus_discrepancies_boxplot_by_sample(
            general_data=general_data,
            labs=labs,
            comp_code=comp_code,
            figures_dir=figures_dir,
        )
        make_component_consensus_discrepancies_stacked_by_sample(
            general_data=general_data,
            labs=labs,
            comp_code=comp_code,
            figures_dir=figures_dir,
        )
        make_component_consensus_discrepancy_type_boxplot(
            general_data=general_data,
            labs=labs,
            comp_code=comp_code,
            figures_dir=figures_dir,
        )
        if comp_code.startswith("SARS"):
            make_component_variant_discrepancies_stacked_by_sample(
                general_data=general_data,
                labs=labs,
                comp_code=comp_code,
                figures_dir=figures_dir,
            )
            make_component_variant_discrepancy_type_boxplot(
                labs=labs,
                comp_code=comp_code,
                figures_dir=figures_dir,
            )
        if comp_code.startswith("FLU"):
            make_component_influenza_variant_reporting_summary(
                general_data=general_data,
                labs=labs,
                comp_code=comp_code,
                figures_dir=figures_dir,
            )
        make_component_typing_outcome_stacked_bar_by_sample(
            general_data=general_data,
            labs=labs,
            comp_code=comp_code,
            figures_dir=figures_dir,
        )
        make_component_qc_match_by_sample_plot(
            general_data=general_data,
            comp_code=comp_code,
            figures_dir=figures_dir,
        )
        make_component_bioinformatics_protocol_metric_boxplots(
            labs=labs,
            comp_code=comp_code,
            figures_dir=figures_dir,
        )
        for benchmark_key in [
            "dehosting",
            "preprocessing",
            "mapping",
            "assembly",
            "consensus_software",
            "variant_calling",
            "clade_assignment",
            "lineage_assignment",
            "type_assignment",
            "subtype_assignment",
        ]:
            make_component_benchmark_metric_boxplots(
                labs=labs,
                comp_code=comp_code,
                figures_dir=figures_dir,
                benchmark_key=benchmark_key,
            )


def collect_lab_consensus_metric_distribution_data(
    general_data: Dict[str, Any],
    labs: List[Dict[str, Any]],
    lab: Dict[str, Any],
    comp_code: str,
    metric_key: str,
    y_limit: Optional[float] = None,
) -> Dict[str, Any]:
    sample_order = list(lab.get("components", {}).get(comp_code, {}).get("samples", {}).keys())
    sample_names = []
    network_data = []
    lab_values = []
    outlier_annotations: List[tuple[int, float]] = []
    lab_outlier_annotations: List[tuple[int, float]] = []

    for sample_id in sample_order:
        sample_values = []
        for network_lab in labs:
            sample = network_lab.get("components", {}).get(comp_code, {}).get("samples", {}).get(sample_id)
            if not sample:
                continue
            value = safe_number(sample.get("consensus", {}).get(metric_key))
            if value is not None:
                sample_values.append(value)

        lab_sample = lab.get("components", {}).get(comp_code, {}).get("samples", {}).get(sample_id)
        lab_value = safe_number(lab_sample.get("consensus", {}).get(metric_key)) if lab_sample else None

        if not sample_values and lab_value is None:
            continue

        if y_limit is not None and sample_values:
            outliers_above_limit = sorted([value for value in sample_values if value > y_limit], reverse=True)
            plotted_values = [value for value in sample_values if value <= y_limit]
            if outliers_above_limit and plotted_values:
                sample_values = plotted_values
                outlier_annotations.append((len(sample_names) + 1, outliers_above_limit[0]))
        if y_limit is not None and lab_value is not None and lab_value > y_limit:
            lab_outlier_annotations.append((len(sample_names) + 1, lab_value))

        sample_names.append(sample_id)
        network_data.append(sample_values)
        lab_values.append(lab_value)

    return {
        "sample_names": sample_names,
        "network_data": network_data,
        "lab_values": lab_values,
        "outlier_annotations": outlier_annotations,
        "lab_outlier_annotations": lab_outlier_annotations,
        "has_lab_values": any(value is not None for value in lab_values),
    }

def make_lab_consensus_distribution_panel_plot(
    general_data: Dict[str, Any],
    labs: List[Dict[str, Any]],
    lab: Dict[str, Any],
    comp_code: str,
    figures_dir: str | Path,
    output_filename: str = "consensus_distribution_panel.png",
) -> str:
    lab_code = get_lab_identifier(lab)
    output_dir = ensure_lab_component_figures_dir(figures_dir, lab_code, comp_code)
    output_path = output_dir / output_filename

    discrepancy_data = collect_lab_consensus_metric_distribution_data(
        general_data=general_data,
        labs=labs,
        lab=lab,
        comp_code=comp_code,
        metric_key="total_discrepancies",
        y_limit=COMPONENT_CONSENSUS_SAMPLE_Y_LIMITS.get(comp_code),
    )
    identity_data = collect_lab_consensus_metric_distribution_data(
        general_data=general_data,
        labs=labs,
        lab=lab,
        comp_code=comp_code,
        metric_key="genome_identity_pct",
    )

    if (
        (not discrepancy_data["sample_names"] or not discrepancy_data["has_lab_values"])
        and (not identity_data["sample_names"] or not identity_data["has_lab_values"])
    ):
        return str(output_path)

    max_samples = max(
        len(discrepancy_data["sample_names"]),
        len(identity_data["sample_names"]),
        1,
    )
    fig, axes = plt.subplots(1, 2, figsize=(max(12, max_samples * 2.0), 6))
    panel_specs = [
        (
            axes[0],
            "A. Consensus discrepancies",
            "Consensus discrepancies",
            discrepancy_data,
            COMPONENT_CONSENSUS_SAMPLE_Y_LIMITS.get(comp_code),
            False,
        ),
        (
            axes[1],
            "B. Genome identity",
            "Genome identity (%)",
            identity_data,
            None,
            True,
        ),
    ]

    for ax, title, ylabel, panel_data, y_limit, percent_axis in panel_specs:
        sample_names = panel_data["sample_names"]
        if not sample_names or not panel_data["has_lab_values"]:
            ax.set_visible(False)
            continue

        network_data = panel_data["network_data"]
        lab_values = panel_data["lab_values"]
        bp = ax.boxplot(
            network_data,
            labels=sample_names,
            showfliers=True,
            patch_artist=True,
        )
        style_boxplot(bp, [comp_code] * len(sample_names), ax=ax)
        add_component_boxplot_points(
            ax,
            bp,
            network_data,
            list(range(1, len(sample_names) + 1)),
            [comp_code] * len(sample_names),
        )
        add_lab_result_diamond(
            ax,
            list(range(1, len(sample_names) + 1)),
            lab_values,
            y_upper=y_limit,
        )

        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.tick_params(axis="x", rotation=0, labelsize=9)

        if percent_axis:
            style_percent_boxplot_axis(ax)
        elif y_limit is not None:
            ax.set_ylim(0, y_limit)
            combined_annotations = list(panel_data["outlier_annotations"])
            for annotation in panel_data.get("lab_outlier_annotations", []):
                if annotation not in combined_annotations:
                    combined_annotations.append(annotation)
            annotate_outlier_caps(
                ax,
                combined_annotations,
                y_limit,
                COMPONENT_BOX_COLORS.get(comp_code, CBF_COLORS["outlier"]),
            )

    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    return str(output_path)


def make_lab_consensus_discrepancy_breakdown_plot(
    lab: Dict[str, Any],
    comp_code: str,
    figures_dir: str | Path,
    output_filename: str = "consensus_discrepancy_breakdown_by_sample.png",
) -> str:
    lab_code = get_lab_identifier(lab)
    output_dir = ensure_lab_component_figures_dir(figures_dir, lab_code, comp_code)
    output_path = output_dir / output_filename

    comp = lab.get("components", {}).get(comp_code, {})
    samples = comp.get("samples", {})

    sample_names = []
    stacked_values = {key: [] for key in CONSENSUS_DISCREPANCY_TYPE_ORDER}

    for sample_id, sample in samples.items():
        breakdown = sample.get("consensus", {}).get("discrepancy_breakdown", {})
        raw_values = [
            safe_number(breakdown.get(key))
            for key in CONSENSUS_DISCREPANCY_TYPE_ORDER
        ]
        if all(value is None for value in raw_values):
            continue

        sample_counts = {}
        for key in CONSENSUS_DISCREPANCY_TYPE_ORDER:
            value = safe_number(breakdown.get(key))
            sample_counts[key] = 0.0 if value is None else value

        sample_names.append(sample_id)
        for key in CONSENSUS_DISCREPANCY_TYPE_ORDER:
            stacked_values[key].append(sample_counts[key])

    if not sample_names:
        return str(output_path)

    fig, ax = plt.subplots(figsize=(max(10, len(sample_names) * 1.8), 6))
    x_positions = np.arange(len(sample_names))
    bottoms = np.zeros(len(sample_names))

    for key in CONSENSUS_DISCREPANCY_TYPE_ORDER:
        values = stacked_values[key]
        ax.bar(
            x_positions,
            values,
            bottom=bottoms,
            color=CONSENSUS_DISCREPANCY_TYPE_COLORS[key],
            edgecolor="#4A4A4A",
            linewidth=0.8,
            label=CONSENSUS_DISCREPANCY_TYPE_LABELS[key].replace("\n", " "),
        )
        bottoms = bottoms + np.array(values)

    ax.set_title(f"{comp_code} discrepancy type breakdown for {lab_code}")
    ax.set_ylabel("Discrepancies (n)")
    ax.set_xticks(x_positions)
    ax.set_xticklabels(sample_names, fontsize=9)
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.18),
        ncol=3,
        frameon=False,
        fontsize=8,
    )

    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    return str(output_path)


def get_variant_reporting_mode_label(sample: Dict[str, Any]) -> Optional[str]:
    variants = sample.get("variants", {})
    if variants.get("high_and_low_freq") is True:
        return "High and low frequency"
    if variants.get("high_freq_only") is True:
        return "High frequency only"
    if variants.get("low_freq_only") is True:
        return "Low frequency only"
    return None


def collect_lab_classification_distribution_data(
    general_data: Dict[str, Any],
    labs: List[Dict[str, Any]],
    lab: Dict[str, Any],
    comp_code: str,
    mode: str,
) -> Dict[str, Any]:
    network_samples = collect_classification_sample_outcomes(general_data, labs, comp_code, mode).get("samples", [])
    lab_samples = lab.get("components", {}).get(comp_code, {}).get("samples", {})
    sample_names: List[str] = []
    match_rates: List[float] = []
    discrepancy_rates: List[float] = []
    null_rates: List[float] = []
    lab_values: List[float] = []

    for sample in network_samples:
        sample_id = sample.get("collecting_lab_sample_id")
        if not sample_id:
            continue

        lab_sample = lab_samples.get(sample_id, {})
        cls = lab_sample.get("classification", {})
        if mode == "lineage_type":
            lab_assignment = cls.get("lineage_assignment")
            lab_match = cls.get("lineage_match")
        else:
            lab_assignment = cls.get("clade_assignment")
            lab_match = cls.get("clade_match")

        sample_names.append(sample_id)
        match_pct = safe_number(sample.get("hit_pct")) or 0.0
        discrepancy_pct = safe_number(sample.get("discrepancy_pct")) or 0.0
        null_pct = safe_number(sample.get("null_pct")) or 0.0
        match_rates.append(match_pct)
        discrepancy_rates.append(discrepancy_pct)
        null_rates.append(null_pct)

        if not is_meaningful(lab_assignment):
            lab_values.append(match_pct + discrepancy_pct + null_pct / 2.0)
        elif lab_match is True:
            lab_values.append(match_pct / 2.0)
        else:
            lab_values.append(match_pct + discrepancy_pct / 2.0)

    return {
        "sample_names": sample_names,
        "match_rates": match_rates,
        "discrepancy_rates": discrepancy_rates,
        "null_rates": null_rates,
        "lab_values": lab_values,
        "lab_positions": [
            (idx - 0.22) if value is not None else idx
            for idx, value in enumerate(lab_values)
        ],
        "has_lab_values": any(value is not None for value in lab_values),
    }


def make_lab_classification_dimension_concordance_plot(
    general_data: Dict[str, Any],
    labs: List[Dict[str, Any]],
    lab: Dict[str, Any],
    comp_code: str,
    figures_dir: str | Path,
    output_filename: str = "classification_dimension_concordance.png",
) -> str:
    lab_code = get_lab_identifier(lab)
    output_dir = ensure_lab_component_figures_dir(figures_dir, lab_code, comp_code)
    output_path = output_dir / output_filename

    panel_specs = [
        (
            "A. Lineage/Subtype assignments",
            collect_lab_classification_distribution_data(general_data, labs, lab, comp_code, "lineage_type"),
        ),
        (
            "B. Clade assignments",
            collect_lab_classification_distribution_data(general_data, labs, lab, comp_code, "clade"),
        ),
    ]
    visible_panels = [spec for spec in panel_specs if spec[1]["sample_names"] and spec[1]["has_lab_values"]]
    if not visible_panels:
        return str(output_path)

    max_samples = max(len(panel_data["sample_names"]) for _, panel_data in visible_panels)
    fig_width = max(12, max_samples * 1.5)
    fig, axes = plt.subplots(1, 2, figsize=(fig_width, 6.6), sharey=True)
    if not isinstance(axes, np.ndarray):
        axes = np.array([axes])

    legend_handles = None
    width = 0.62

    for ax, (title, panel_data) in zip(axes, panel_specs):
        sample_names = panel_data["sample_names"]
        if not sample_names or not panel_data["has_lab_values"]:
            ax.set_visible(False)
            continue

        x_positions = np.arange(len(sample_names))
        match_bars = ax.bar(
            x_positions,
            panel_data["match_rates"],
            width=width,
            color=CBF_COLORS["match"],
            label="Match",
        )
        discrepancy_bars = ax.bar(
            x_positions,
            panel_data["discrepancy_rates"],
            width=width,
            bottom=panel_data["match_rates"],
            color=CBF_COLORS["discrepancy"],
            label="Discrepancy",
        )
        null_bars = ax.bar(
            x_positions,
            panel_data["null_rates"],
            width=width,
            bottom=np.array(panel_data["match_rates"]) + np.array(panel_data["discrepancy_rates"]),
            color=CBF_COLORS["null"],
            label="Not provided",
        )
        legend_handles = (match_bars[0], discrepancy_bars[0], null_bars[0])

        add_lab_result_diamond(
            ax=ax,
            positions=panel_data["lab_positions"],
            values=panel_data["lab_values"],
            y_upper=None,
        )

        ax.set_xticks(x_positions)
        ax.set_xticklabels(sample_names, rotation=45, ha="right")
        ax.set_xlabel("Sample")
        ax.set_ylim(0, 100)
        ax.set_title(title)

        for bar, value in zip(match_bars, panel_data["match_rates"]):
            if value <= 0:
                continue
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                value / 2,
                f"{value:.1f}%",
                ha="center",
                va="center",
                fontsize=8,
                color="white",
                fontweight="bold",
            )

        for bar, match_value, discrepancy_value in zip(
            discrepancy_bars,
            panel_data["match_rates"],
            panel_data["discrepancy_rates"],
        ):
            if discrepancy_value <= 0:
                continue
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                match_value + discrepancy_value / 2,
                f"{discrepancy_value:.1f}%",
                ha="center",
                va="center",
                fontsize=8,
                color="white",
                fontweight="bold",
            )

        for bar, match_value, discrepancy_value, null_value in zip(
            null_bars,
            panel_data["match_rates"],
            panel_data["discrepancy_rates"],
            panel_data["null_rates"],
        ):
            if null_value <= 0:
                continue
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                match_value + discrepancy_value + null_value / 2,
                f"{null_value:.1f}%",
                ha="center",
                va="center",
                fontsize=8,
                color="white",
                fontweight="bold",
            )

    axes[0].set_ylabel("Assignments (%)")
    if legend_handles is not None:
        fig.legend(
            handles=list(legend_handles),
            labels=["Match", "Discrepancy", "Not provided"],
            loc="lower center",
            bbox_to_anchor=(0.5, 0.02),
            ncol=3,
            frameon=False,
        )
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    return str(output_path)


WORKFLOW_POSITIONING_METRICS = [
    (
        "A. Total consensus discrepancies",
        "Total discrepancies",
        lambda comp: comp.get("total_number_discrepancies_consensus"),
        CBF_COLORS["discrepancy"],
        None,
    ),
    (
        "B. Median genome identity",
        "Genome identity (%)",
        lambda comp: comp.get("median_genome_identity_pct"),
        CBF_COLORS["match"],
        (0, 100),
    ),
    (
        "C. Total classification matches",
        "Classification matches",
        lambda comp: comp.get("total_classification_matches"),
        CBF_COLORS["high_and_low_freq"],
        None,
    ),
    (
        "D. Metadata completeness",
        "Metadata completeness (%)",
        lambda comp: comp.get("metadata", {}).get("completeness_pct"),
        CBF_COLORS["box_flu2"],
        (0, 100),
    ),
]


def collect_lab_workflow_positioning_data(
    labs: List[Dict[str, Any]],
    lab: Dict[str, Any],
    comp_code: str,
    extractor,
) -> Dict[str, Any]:
    network_values = []
    for network_lab in labs:
        comp = network_lab.get("components", {}).get(comp_code)
        if not comp:
            continue
        value = safe_number(extractor(comp))
        if value is not None:
            network_values.append(value)

    lab_comp = lab.get("components", {}).get(comp_code, {})
    lab_value = safe_number(extractor(lab_comp))

    return {
        "network_values": network_values,
        "lab_value": lab_value,
    }


def make_lab_workflow_positioning_boxplot(
    labs: List[Dict[str, Any]],
    lab: Dict[str, Any],
    comp_code: str,
    figures_dir: str | Path,
    output_filename: str = "workflow_positioning_boxplots.png",
) -> str:
    lab_code = get_lab_identifier(lab)
    output_dir = ensure_lab_component_figures_dir(figures_dir, lab_code, comp_code)
    output_path = output_dir / output_filename
    component_color = COMPONENT_BOX_COLORS.get(comp_code, CBF_COLORS["box_default"])

    panel_data_specs = []
    for title, ylabel, extractor, color, y_limits in WORKFLOW_POSITIONING_METRICS:
        metric_data = collect_lab_workflow_positioning_data(
            labs=labs,
            lab=lab,
            comp_code=comp_code,
            extractor=extractor,
        )
        if not metric_data["network_values"] or metric_data["lab_value"] is None:
            continue
        panel_data_specs.append((title, ylabel, metric_data, color, y_limits))

    if not panel_data_specs:
        return str(output_path)

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    axes = axes.flatten()

    for ax, (title, ylabel, metric_data, color, y_limits) in zip(axes, panel_data_specs):
        data = [metric_data["network_values"]]
        bp = ax.boxplot(
            data,
            patch_artist=True,
            widths=0.5,
            showfliers=True,
        )
        for patch in bp["boxes"]:
            patch.set_facecolor(component_color)
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
        style_boxplot_axes(ax)

        add_colored_boxplot_points(
            ax,
            bp,
            data,
            [1],
            [component_color],
        )
        add_lab_result_diamond(
            ax=ax,
            positions=[0.84],
            values=[metric_data["lab_value"]],
            y_upper=y_limits[1] if y_limits is not None else None,
        )

        ax.set_xticks([1])
        ax.set_xticklabels(["Network"])
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        if y_limits is not None:
            if ylabel.endswith("(%)"):
                style_percent_boxplot_axis(ax)
            else:
                ax.set_ylim(*y_limits)

    for ax in axes[len(panel_data_specs):]:
        ax.set_visible(False)

    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    return str(output_path)


def collect_lab_qc_match_distribution_data(
    general_data: Dict[str, Any],
    lab: Dict[str, Any],
    comp_code: str,
) -> Dict[str, Any]:
    comp_data = general_data.get("components", {}).get(comp_code, {})
    network_samples = comp_data.get("qc", {}).get("samples", [])
    lab_samples = lab.get("components", {}).get(comp_code, {}).get("samples", {})

    sample_names: List[str] = []
    match_rates: List[float] = []
    discrepancy_rates: List[float] = []
    lab_values: List[float] = []
    lab_positions: List[float] = []

    for sample in network_samples:
        sample_id = sample.get("collecting_lab_sample_id") or sample.get("sample_id")
        if not sample_id:
            continue

        network_match_pct = safe_number(sample.get("match_rate_pct"))
        lab_sample = lab_samples.get(sample_id, {})
        lab_match = lab_sample.get("qc_match")

        if network_match_pct is None or lab_match is None:
            continue

        sample_names.append(sample_id)
        match_rates.append(network_match_pct)
        discrepancy_rates.append(100.0 - network_match_pct)

        if lab_match is True:
            lab_values.append(network_match_pct / 2.0)
        else:
            lab_values.append(network_match_pct + (100.0 - network_match_pct) / 2.0)
        lab_positions.append(len(sample_names) - 1 - 0.22)

    return {
        "sample_names": sample_names,
        "match_rates": match_rates,
        "discrepancy_rates": discrepancy_rates,
        "lab_values": lab_values,
        "lab_positions": lab_positions,
        "has_lab_values": any(value is not None for value in lab_values),
    }


def make_lab_qc_match_rate_plot(
    general_data: Dict[str, Any],
    lab: Dict[str, Any],
    comp_code: str,
    figures_dir: str | Path,
    output_filename: str = "qc_match_rate.png",
) -> str:
    lab_code = get_lab_identifier(lab)
    output_dir = ensure_lab_component_figures_dir(figures_dir, lab_code, comp_code)
    output_path = output_dir / output_filename

    panel_data = collect_lab_qc_match_distribution_data(
        general_data=general_data,
        lab=lab,
        comp_code=comp_code,
    )
    if not panel_data["sample_names"] or not panel_data["has_lab_values"]:
        return str(output_path)

    sample_names = panel_data["sample_names"]
    x_positions = np.arange(len(sample_names))
    width = 0.62

    fig, ax = plt.subplots(figsize=(max(10, len(sample_names) * 1.5), 6.4))

    match_bars = ax.bar(
        x_positions,
        panel_data["match_rates"],
        width=width,
        color=CBF_COLORS["match"],
        label="Match",
    )
    discrepancy_bars = ax.bar(
        x_positions,
        panel_data["discrepancy_rates"],
        width=width,
        bottom=panel_data["match_rates"],
        color=CBF_COLORS["discrepancy"],
        label="Discrepancy",
    )
    add_lab_result_diamond(
        ax=ax,
        positions=panel_data["lab_positions"],
        values=panel_data["lab_values"],
        y_upper=None,
    )

    ax.set_xticks(x_positions)
    ax.set_xticklabels(sample_names, rotation=45, ha="right")
    ax.set_xlabel("Sample")
    ax.set_ylabel("QC outcomes (%)")
    ax.set_ylim(0, 100)
    ax.set_title(f"{comp_code} QC concordance across participating laboratories")

    for bar, value in zip(match_bars, panel_data["match_rates"]):
        if value <= 0:
            continue
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            value / 2,
            f"{value:.1f}%",
            ha="center",
            va="center",
            fontsize=8,
            color="white",
            fontweight="bold",
        )

    for bar, match_value, discrepancy_value in zip(
        discrepancy_bars,
        panel_data["match_rates"],
        panel_data["discrepancy_rates"],
    ):
        if discrepancy_value <= 0:
            continue
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            match_value + discrepancy_value / 2,
            f"{discrepancy_value:.1f}%",
            ha="center",
            va="center",
            fontsize=8,
            color="white",
            fontweight="bold",
        )

    ax.legend(
        loc="lower center",
        bbox_to_anchor=(0.5, -0.3),
        ncol=2,
        frameon=False,
    )
    fig.tight_layout(rect=(0, 0.14, 1, 1))
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    return str(output_path)


def make_lab_metadata_metrics_panel_plot(
    labs: List[Dict[str, Any]],
    lab: Dict[str, Any],
    comp_code: str,
    figures_dir: str | Path,
    output_filename: str = "metadata_metrics_panel.png",
) -> str:
    panel_specs = [
        (
            "A. Genome >10x",
            "Genome >10x (%)",
            lambda s: s.get("metadata_metrics", {}).get("per_genome_greater_10x"),
            100.0,
            True,
        ),
        (
            "B. Depth of coverage",
            "Depth of coverage",
            lambda s: s.get("metadata_metrics", {}).get("depth_of_coverage_value"),
            None,
            False,
        ),
        (
            "C. Ns",
            "Ns (%)",
            lambda s: s.get("metadata_metrics", {}).get("per_Ns"),
            100.0,
            True,
        ),
        (
            "D. Viral reads",
            "Reads virus (%)",
            lambda s: s.get("metadata_metrics", {}).get("per_reads_virus"),
            100.0,
            True,
        ),
        (
            "E. Host reads",
            "Reads host (%)",
            lambda s: s.get("metadata_metrics", {}).get("per_reads_host"),
            100.0,
            True,
        ),
    ]

    return make_lab_variant_boxplot_panel_figure(
        labs=labs,
        lab=lab,
        comp_code=comp_code,
        figures_dir=figures_dir,
        output_filename=output_filename,
        panel_specs=panel_specs,
    )


def collect_lab_variant_metric_distribution_data(
    labs: List[Dict[str, Any]],
    lab: Dict[str, Any],
    comp_code: str,
    extractor,
    y_limit: Optional[float] = None,
) -> Dict[str, Any]:
    sample_order = list(lab.get("components", {}).get(comp_code, {}).get("samples", {}).keys())
    sample_names = []
    network_data = []
    lab_values = []
    outlier_annotations: List[tuple[int, float]] = []
    lab_outlier_annotations: List[tuple[int, float]] = []

    for sample_id in sample_order:
        sample_values = []
        for network_lab in labs:
            sample = network_lab.get("components", {}).get(comp_code, {}).get("samples", {}).get(sample_id)
            if not sample:
                continue
            value = safe_number(extractor(sample))
            if value is not None:
                sample_values.append(value)

        lab_sample = lab.get("components", {}).get(comp_code, {}).get("samples", {}).get(sample_id)
        lab_value = safe_number(extractor(lab_sample)) if lab_sample else None

        if not sample_values and lab_value is None:
            continue

        if y_limit is not None and sample_values:
            outliers_above_limit = sorted([value for value in sample_values if value > y_limit], reverse=True)
            plotted_values = [value for value in sample_values if value <= y_limit]
            if outliers_above_limit and plotted_values:
                sample_values = plotted_values
                outlier_annotations.append((len(sample_names) + 1, outliers_above_limit[0]))
        if y_limit is not None and lab_value is not None and lab_value > y_limit:
            lab_outlier_annotations.append((len(sample_names) + 1, lab_value))

        sample_names.append(sample_id)
        network_data.append(sample_values)
        lab_values.append(lab_value)

    return {
        "sample_names": sample_names,
        "network_data": network_data,
        "lab_values": lab_values,
        "outlier_annotations": outlier_annotations,
        "lab_outlier_annotations": lab_outlier_annotations,
        "has_lab_values": any(value is not None for value in lab_values),
    }


def make_lab_variant_boxplot_panel_figure(
    labs: List[Dict[str, Any]],
    lab: Dict[str, Any],
    comp_code: str,
    figures_dir: str | Path,
    output_filename: str,
    panel_specs: List[tuple[str, str, Any, Optional[float], bool]],
) -> str:
    lab_code = get_lab_identifier(lab)
    output_dir = ensure_lab_component_figures_dir(figures_dir, lab_code, comp_code)
    output_path = output_dir / output_filename

    panel_data_specs = []
    max_samples = 1
    for title, ylabel, extractor, y_limit, percent_axis in panel_specs:
        panel_data = collect_lab_variant_metric_distribution_data(
            labs=labs,
            lab=lab,
            comp_code=comp_code,
            extractor=extractor,
            y_limit=y_limit,
        )
        if panel_data["sample_names"]:
            if not panel_data["has_lab_values"]:
                continue
            max_samples = max(max_samples, len(panel_data["sample_names"]))
            panel_data_specs.append((title, ylabel, panel_data, y_limit, percent_axis))

    if not panel_data_specs:
        return str(output_path)

    n_panels = len(panel_data_specs)
    ncols = 1 if n_panels == 1 else 2
    nrows = int(np.ceil(n_panels / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(max(12, max_samples * 2.0), 4.8 * nrows))
    axes = np.atleast_1d(axes).flatten()

    for ax, (title, ylabel, panel_data, y_limit, percent_axis) in zip(axes, panel_data_specs):
        sample_names = panel_data["sample_names"]
        network_data = panel_data["network_data"]
        lab_values = panel_data["lab_values"]
        bp = ax.boxplot(
            network_data,
            labels=sample_names,
            showfliers=True,
            patch_artist=True,
        )
        style_boxplot(bp, [comp_code] * len(sample_names), ax=ax)
        add_component_boxplot_points(
            ax,
            bp,
            network_data,
            list(range(1, len(sample_names) + 1)),
            [comp_code] * len(sample_names),
        )
        add_lab_result_diamond(
            ax,
            list(range(1, len(sample_names) + 1)),
            lab_values,
            y_upper=y_limit,
        )

        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.tick_params(axis="x", rotation=0, labelsize=9)

        if percent_axis:
            style_percent_boxplot_axis(ax)
        elif y_limit is not None:
            ax.set_ylim(0, y_limit)
            combined_annotations = list(panel_data["outlier_annotations"])
            for annotation in panel_data.get("lab_outlier_annotations", []):
                if annotation not in combined_annotations:
                    combined_annotations.append(annotation)
            annotate_outlier_caps(
                ax,
                combined_annotations,
                y_limit,
                COMPONENT_BOX_COLORS.get(comp_code, CBF_COLORS["outlier"]),
            )

    for ax in axes[len(panel_data_specs):]:
        ax.set_visible(False)

    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    return str(output_path)


def make_lab_variant_metrics_distribution_plot(
    labs: List[Dict[str, Any]],
    lab: Dict[str, Any],
    comp_code: str,
    figures_dir: str | Path,
    output_filename: str = "variant_metrics_distribution.png",
) -> str:
    if comp_code.startswith("SARS"):
        panel_specs = [
            (
                "A. Total discrepancies",
                "Total discrepancies",
                lambda s: s.get("variants", {}).get("total_discrepancies"),
                None,
                False,
            ),
            (
                "B. Successful hits",
                "Successful hits",
                lambda s: s.get("variants", {}).get("successful_hits"),
                None,
                False,
            ),
        ]
    else:
        panel_specs = [
            (
                "A. Reported variants (AF >=75%)",
                "Reported variants (AF >=75%)",
                lambda s: s.get("variants", {}).get("number_of_variants_in_consensus"),
                None,
                False,
            ),
            (
                "B. Variants in VCF (AF >=75%)",
                "Variants in VCF (AF >=75%)",
                lambda s: s.get("variants", {}).get("number_of_variants_in_consensus_vcf"),
                None,
                False,
            ),
            (
                "C. Metadata-VCF discrepancies",
                "Metadata-VCF discrepancies",
                lambda s: s.get("variants", {}).get("discrepancies_in_reported_variants"),
                None,
                False,
            ),
            (
                "D. Total variants in VCF",
                "Total variants in VCF",
                lambda s: s.get("variants", {}).get("number_of_variants_in_vcf"),
                None,
                False,
            ),
        ]

    return make_lab_variant_boxplot_panel_figure(
        labs=labs,
        lab=lab,
        comp_code=comp_code,
        figures_dir=figures_dir,
        output_filename=output_filename,
        panel_specs=panel_specs,
    )


def make_lab_variant_metadata_vs_vcf_distribution_plot(
    labs: List[Dict[str, Any]],
    lab: Dict[str, Any],
    comp_code: str,
    figures_dir: str | Path,
    output_filename: str = "variant_metadata_vs_vcf_distribution.png",
) -> str:
    if comp_code.startswith("SARS"):
        panel_specs = [
            (
                "A. Reported variants (AF >=75%)",
                "Number of reported variants",
                lambda s: s.get("variants", {}).get("number_of_variants_in_consensus"),
                None,
                False,
            ),
            (
                "B. Reported variants with effect",
                "Number of reported variants",
                lambda s: s.get("variants", {}).get("number_of_variants_with_effect"),
                None,
                False,
            ),
            (
                "C. Variants in VCF (AF >=75%)",
                "Number of variants",
                lambda s: s.get("variants", {}).get("number_of_variants_in_consensus_vcf"),
                None,
                False,
            ),
            (
                "D. Variants with effect VCF",
                "Number of variants",
                lambda s: s.get("variants", {}).get("number_of_variants_with_effect_vcf"),
                None,
                False,
            ),
            (
                "E. Metadata-VCF discrepancies (AF >=75%)",
                "Number of metadata-VCF discrepancies",
                lambda s: s.get("variants", {}).get("discrepancies_in_reported_variants"),
                None,
                False,
            ),
            (
                "F. Metadata-VCF discrepancies variants with effect",
                "Number of metadata-VCF discrepancies",
                lambda s: s.get("variants", {}).get("discrepancies_in_reported_variants_effect"),
                None,
                False,
            ),
        ]
    else:
        panel_specs = [
            (
                "A. Reported variants (AF >=75%)",
                "Reported variants (AF >=75%)",
                lambda s: s.get("variants", {}).get("number_of_variants_in_consensus"),
                None,
                False,
            ),
            (
                "B. Variants in VCF (AF >=75%)",
                "Variants in VCF (AF >=75%)",
                lambda s: s.get("variants", {}).get("number_of_variants_in_consensus_vcf"),
                None,
                False,
            ),
            (
                "C. Variants with effect",
                "Variants with effect",
                lambda s: s.get("variants", {}).get("number_of_variants_with_effect"),
                None,
                False,
            ),
            (
                "D. Metadata-VCF discrepancies",
                "Metadata-VCF discrepancies",
                lambda s: s.get("variants", {}).get("discrepancies_in_reported_variants"),
                None,
                False,
            ),
            (
                "E. Total variants in VCF",
                "Total variants in VCF",
                lambda s: s.get("variants", {}).get("number_of_variants_in_vcf"),
                None,
                False,
            ),
        ]

    return make_lab_variant_boxplot_panel_figure(
        labs=labs,
        lab=lab,
        comp_code=comp_code,
        figures_dir=figures_dir,
        output_filename=output_filename,
        panel_specs=panel_specs,
    )


def generate_individual_lab_figures(
    general_data: Dict[str, Any],
    labs: List[Dict[str, Any]],
    figures_dir: str | Path = "figures",
) -> None:
    for lab in labs:
        lab_code = get_lab_identifier(lab)
        if not lab_code:
            continue

        for comp_code in lab.get("components", {}).keys():
            comp = lab.get("components", {}).get(comp_code, {})
            make_lab_consensus_distribution_panel_plot(
                general_data=general_data,
                labs=labs,
                lab=lab,
                comp_code=comp_code,
                figures_dir=figures_dir,
            )
            make_lab_consensus_discrepancy_breakdown_plot(
                lab=lab,
                comp_code=comp_code,
                figures_dir=figures_dir,
            )
            make_lab_classification_dimension_concordance_plot(
                general_data=general_data,
                labs=labs,
                lab=lab,
                comp_code=comp_code,
                figures_dir=figures_dir,
            )
            make_lab_workflow_positioning_boxplot(
                labs=labs,
                lab=lab,
                comp_code=comp_code,
                figures_dir=figures_dir,
            )
            make_lab_qc_match_rate_plot(
                general_data=general_data,
                lab=lab,
                comp_code=comp_code,
                figures_dir=figures_dir,
            )
            make_lab_metadata_metrics_panel_plot(
                labs=labs,
                lab=lab,
                comp_code=comp_code,
                figures_dir=figures_dir,
            )
            if (safe_int(comp.get("metadata", {}).get("vcf_submitted")) or 0) >= 1:
                if comp_code.startswith("SARS"):
                    make_lab_variant_metrics_distribution_plot(
                        labs=labs,
                        lab=lab,
                        comp_code=comp_code,
                        figures_dir=figures_dir,
                        output_filename="variant_metadata_vs_vcf_distribution.png",
                    )
                    make_lab_variant_metadata_vs_vcf_distribution_plot(
                        labs=labs,
                        lab=lab,
                        comp_code=comp_code,
                        figures_dir=figures_dir,
                        output_filename="variant_metrics_distribution.png",
                    )
                else:
                    make_lab_variant_metadata_vs_vcf_distribution_plot(
                        labs=labs,
                        lab=lab,
                        comp_code=comp_code,
                        figures_dir=figures_dir,
                        output_filename="variant_metrics_distribution.png",
                    )


def collect_software_groups(
    participating_labs: List[Dict[str, Any]],
    comp_code: str,
    name_field: str,
    version_field: Optional[str] = None,
    db_version_field: Optional[str] = None,
) -> Dict[tuple, List[Dict[str, Any]]]:
    groups = defaultdict(list)

    for lab in participating_labs:
        lab_id = lab["lab"]["submitting_institution_id"]
        comp = lab.get("components", {}).get(comp_code)
        if not comp:
            continue

        for sample_id, sample in comp.get("samples", {}).items():
            sb = sample.get("software_benchmarking", {})
            key = software_signature(
                sb.get(name_field),
                sb.get(version_field) if version_field else None,
                sb.get(db_version_field) if db_version_field else None,
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

    for (name, version, database_version), records in groups.items():
        entry = {
            "name": name,
            "n_labs": len({r["lab_id"] for r in records}),
        }
        if version is not None:
            entry["version"] = version
        if database_version is not None:
            entry["database_version"] = database_version

        entry.update(metrics_builder(records))
        entries.append(entry)

    return sorted(entries, key=lambda x: (x.get("name") or "", x.get("version") or "", x.get("database_version") or ""))


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
    style_boxplot(bp, component_names, ax=plt.gca())
    add_component_boxplot_points(
        plt.gca(),
        bp,
        data,
        list(range(1, len(component_names) + 1)),
        component_names,
    )

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
    params_slots_total = 0
    params_slots_missing = 0

    all_workflows = set()
    consensus_softwares = set()
    variant_softwares = set()
    lineage_assignment_softwares = set()
    clade_assignment_softwares = set()
    type_assignment_softwares = set()
    subtype_assignment_softwares = set()

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

                sb = sample.get("software_benchmarking", {})

                software_name_slots = [
                    is_meaningful(sb.get("bioinformatics_protocol_software_name")),
                    is_meaningful(sb.get("dehosting_method_software_name")),
                    is_meaningful(sb.get("preprocessing_software_name")),
                    is_meaningful(sb.get("mapping_software_name")) or is_meaningful(sb.get("assembly")),
                    is_meaningful(sb.get("variant_calling_software_name")),
                    is_meaningful(sb.get("consensus_sequence_software_name")),
                    is_meaningful(sb.get("clade_assignment_software_name")),
                ]

                software_version_slots = [
                    is_meaningful(sb.get("bioinformatics_protocol_software_version")),
                    is_meaningful(sb.get("dehosting_method_software_version")),
                    is_meaningful(sb.get("preprocessing_software_version")),
                    is_meaningful(sb.get("mapping_software_version")) or is_meaningful(sb.get("assembly_version")),
                    is_meaningful(sb.get("variant_calling_software_version")),
                    is_meaningful(sb.get("consensus_sequence_software_version")),
                    is_meaningful(sb.get("clade_assignment_software_version")),
                ]

                if comp_expected.get("virus") == "SARS-CoV-2":
                    software_name_slots.append(is_meaningful(sb.get("lineage_assignment_software_name")))
                    software_version_slots.append(is_meaningful(sb.get("lineage_assignment_software_version")))
                else:
                    software_name_slots.append(is_meaningful(sb.get("type_assignment_software_name")))
                    software_name_slots.append(is_meaningful(sb.get("subtype_assignment_software_name")))
                    software_version_slots.append(is_meaningful(sb.get("type_assignment_software_version")))
                    software_version_slots.append(is_meaningful(sb.get("subtype_assignment_software_version")))

                software_names_total += len(software_name_slots)
                software_names_count += sum(1 for filled in software_name_slots if filled)

                software_version_total += len(software_version_slots)
                software_version_count += sum(1 for filled in software_version_slots if filled)

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

                params_slots_total += 4

                if not is_meaningful(sb.get("preprocessing_params")):
                    params_slots_missing += 1

                if not (
                    is_meaningful(sb.get("mapping_params"))
                    or is_meaningful(sb.get("assembly_params"))
                ):
                    params_slots_missing += 1

                if not is_meaningful(sb.get("variant_calling_params")):
                    params_slots_missing += 1

                if not is_meaningful(sb.get("consensus_params")):
                    params_slots_missing += 1

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
                        lineage_assignment_softwares.add(sig)
                else:
                    sig = software_signature(
                        sample.get("software_benchmarking", {}).get("type_assignment_software_name"),
                        None,
                    )
                    if sig:
                        type_assignment_softwares.add(sig)

                    sig = software_signature(
                        sample.get("software_benchmarking", {}).get("subtype_assignment_software_name"),
                        sample.get("software_benchmarking", {}).get("subtype_assignment_software_version"),
                    )
                    if sig:
                        subtype_assignment_softwares.add(sig)

                sig = software_signature(
                    sample.get("software_benchmarking", {}).get("clade_assignment_software_name"),
                    sample.get("software_benchmarking", {}).get("clade_assignment_software_version"),
                )

                if sig:
                    clade_assignment_softwares.add(sig)

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

        workflow_total_discrepancies_per_lab = [
            safe_number(lab["components"][comp_code].get("total_number_discrepancies_consensus"))
            for lab in participating_labs
            if safe_number(lab["components"][comp_code].get("total_number_discrepancies_consensus")) is not None
        ]
        workflow_median_identity_per_lab = [
            safe_number(lab["components"][comp_code].get("median_genome_identity_pct"))
            for lab in participating_labs
            if safe_number(lab["components"][comp_code].get("median_genome_identity_pct")) is not None
        ]

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
            "workflow_total_discrepancies_median": median_or_none(workflow_total_discrepancies_per_lab),
            "workflow_total_discrepancies_min": min_or_none(workflow_total_discrepancies_per_lab),
            "workflow_total_discrepancies_max": max_or_none(workflow_total_discrepancies_per_lab),
            "workflow_median_identity_pct_median": median_or_none(workflow_median_identity_per_lab),
            "workflow_median_identity_pct_min": min_or_none(workflow_median_identity_per_lab),
            "workflow_median_identity_pct_max": max_or_none(workflow_median_identity_per_lab),
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
            variant_in_consensus_all = []
            variant_in_consensus_vcf_all = []
            variant_with_effect_all = []
            discrepancy_reported_all = []
            variant_with_effect_vcf_all = []
            discrepancy_reported_effect_all = []
            variant_breakdown_all = defaultdict(list)
            variant_samples = []

            for sample_id, expected_sample in comp_expected["samples"].items():
                tds = []
                successful_hits_vals = []
                variants_in_consensus_vals = []
                variants_in_consensus_vcf_vals = []
                variants_with_effect_vals = []
                discrepancies_in_reported_variants_vals = []
                variants_with_effect_vcf_vals = []
                discrepancies_in_reported_variants_effect_vals = []
                bd = defaultdict(list)

                for lab in participating_labs:
                    sample = lab["components"][comp_code]["samples"].get(sample_id)
                    if not sample:
                        tds.append(None)
                        successful_hits_vals.append(None)
                        variants_in_consensus_vals.append(None)
                        variants_in_consensus_vcf_vals.append(None)
                        variants_with_effect_vals.append(None)
                        discrepancies_in_reported_variants_vals.append(None)
                        variants_with_effect_vcf_vals.append(None)
                        discrepancies_in_reported_variants_effect_vals.append(None)
                        for key in ["wrong_nt", "insertions", "deletions", "missing", "denovo"]:
                            bd[key].append(None)
                        continue

                    var = sample.get("variants", {})

                    td = safe_number(var.get("total_discrepancies"))
                    sh = safe_number(var.get("successful_hits"))
                    nivc = safe_number(var.get("number_of_variants_in_consensus"))
                    nivcv = safe_number(var.get("number_of_variants_in_consensus_vcf"))
                    nwee = safe_number(var.get("number_of_variants_with_effect"))
                    dirv = safe_number(var.get("discrepancies_in_reported_variants"))
                    nweev = safe_number(var.get("number_of_variants_with_effect_vcf"))
                    dirve = safe_number(var.get("discrepancies_in_reported_variants_effect"))

                    tds.append(td)
                    if td is not None:
                        variant_discs.append(td)

                    successful_hits_vals.append(sh)
                    if sh is not None:
                        variant_successful_hits.append(sh)

                    variants_in_consensus_vals.append(nivc)
                    if nivc is not None:
                        variant_in_consensus_all.append(nivc)

                    variants_in_consensus_vcf_vals.append(nivcv)
                    if nivcv is not None:
                        variant_in_consensus_vcf_all.append(nivcv)

                    variants_with_effect_vals.append(nwee)
                    if nwee is not None:
                        variant_with_effect_all.append(nwee)

                    discrepancies_in_reported_variants_vals.append(dirv)
                    if dirv is not None:
                        discrepancy_reported_all.append(dirv)

                    variants_with_effect_vcf_vals.append(nweev)
                    if nweev is not None:
                        variant_with_effect_vcf_all.append(nweev)

                    discrepancies_in_reported_variants_effect_vals.append(dirve)
                    if dirve is not None:
                        discrepancy_reported_effect_all.append(dirve)

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
                variants_in_consensus_summary = summarize_numeric_values(variants_in_consensus_vals, total_count=len(participating_labs))
                variants_in_consensus_vcf_summary = summarize_numeric_values(variants_in_consensus_vcf_vals, total_count=len(participating_labs))
                variants_with_effect_summary = summarize_numeric_values(variants_with_effect_vals, total_count=len(participating_labs))
                discrepancies_in_reported_variants_summary = summarize_numeric_values(discrepancies_in_reported_variants_vals, total_count=len(participating_labs))
                variants_with_effect_vcf_summary = summarize_numeric_values(variants_with_effect_vcf_vals, total_count=len(participating_labs))
                discrepancies_in_reported_variants_effect_summary = summarize_numeric_values(discrepancies_in_reported_variants_effect_vals, total_count=len(participating_labs))

                variant_samples.append({
                    "collecting_lab_sample_id": sample_id,
                    "variants_in_consensus": variants_in_consensus_summary,
                    "variants_in_consensus_vcf": variants_in_consensus_vcf_summary,
                    "variants_with_effect": variants_with_effect_summary,
                    "discrepancies_in_reported_variants": discrepancies_in_reported_variants_summary,
                    "variants_with_effect_vcf": variants_with_effect_vcf_summary,
                    "discrepancies_in_reported_variants_effect": discrepancies_in_reported_variants_effect_summary,
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
            variants_in_consensus_component_summary = summarize_numeric_values(variant_in_consensus_all)
            variants_in_consensus_vcf_component_summary = summarize_numeric_values(variant_in_consensus_vcf_all)
            variants_with_effect_component_summary = summarize_numeric_values(variant_with_effect_all)
            discrepancies_in_reported_variants_component_summary = summarize_numeric_values(discrepancy_reported_all)
            variants_with_effect_vcf_component_summary = summarize_numeric_values(variant_with_effect_vcf_all)
            discrepancies_in_reported_variants_effect_component_summary = summarize_numeric_values(discrepancy_reported_effect_all)
            variant_breakdown_summary = {
                key: summarize_numeric_values(vals)
                for key, vals in variant_breakdown_all.items()
            }
            comp_obj["variant"] = {
                "median_variants_in_consensus": variants_in_consensus_component_summary["median"],
                "min_variants_in_consensus": variants_in_consensus_component_summary["min"],
                "max_variants_in_consensus": variants_in_consensus_component_summary["max"],
                "variants_in_consensus_summary": variants_in_consensus_component_summary,
                "median_variants_in_consensus_vcf": variants_in_consensus_vcf_component_summary["median"],
                "min_variants_in_consensus_vcf": variants_in_consensus_vcf_component_summary["min"],
                "max_variants_in_consensus_vcf": variants_in_consensus_vcf_component_summary["max"],
                "variants_in_consensus_vcf_summary": variants_in_consensus_vcf_component_summary,
                "median_variants_with_effect": variants_with_effect_component_summary["median"],
                "min_variants_with_effect": variants_with_effect_component_summary["min"],
                "max_variants_with_effect": variants_with_effect_component_summary["max"],
                "variants_with_effect_summary": variants_with_effect_component_summary,
                "median_discrepancies_in_reported_variants": discrepancies_in_reported_variants_component_summary["median"],
                "min_discrepancies_in_reported_variants": discrepancies_in_reported_variants_component_summary["min"],
                "max_discrepancies_in_reported_variants": discrepancies_in_reported_variants_component_summary["max"],
                "discrepancies_in_reported_variants_summary": discrepancies_in_reported_variants_component_summary,
                "median_variants_with_effect_vcf": variants_with_effect_vcf_component_summary["median"],
                "min_variants_with_effect_vcf": variants_with_effect_vcf_component_summary["min"],
                "max_variants_with_effect_vcf": variants_with_effect_vcf_component_summary["max"],
                "variants_with_effect_vcf_summary": variants_with_effect_vcf_component_summary,
                "median_discrepancies_in_reported_variants_effect": discrepancies_in_reported_variants_effect_component_summary["median"],
                "min_discrepancies_in_reported_variants_effect": discrepancies_in_reported_variants_effect_component_summary["min"],
                "max_discrepancies_in_reported_variants_effect": discrepancies_in_reported_variants_effect_component_summary["max"],
                "discrepancies_in_reported_variants_effect_summary": discrepancies_in_reported_variants_effect_component_summary,
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
                "fig_discrepancies_stacked_by_sample": f"figures/{comp_code}/variant_discrepancies_stacked_by_sample.png",
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
                "fig_reporting_summary_by_sample": f"figures/{comp_code}/influenza_variant_reporting_summary_by_sample.png",
            }

        classification_matches_per_lab = []
        per_sample_cls = []

        lineage_hits = 0
        lineage_discrepancies = 0
        lineage_null = 0
        clade_hits = 0
        clade_discrepancies = 0
        clade_null = 0

        for lab in participating_labs:
            lab_total_matches = 0

            for sample_id in comp_expected["samples"].keys():
                sample = lab["components"][comp_code]["samples"].get(sample_id)
                if not sample:
                    continue

                cls = sample.get("classification", {})
                nm = safe_int(cls.get("number_matches"))
                lm = cls.get("lineage_match")
                cm = cls.get("clade_match")
                expected_lineage = cls.get("expected_lineage")
                lineage_assignment = cls.get("lineage_assignment")
                expected_clade = cls.get("expected_clade")
                clade_assignment = cls.get("clade_assignment")

                if nm is not None:
                    lab_total_matches += nm

                if is_meaningful(expected_lineage):
                    if not is_meaningful(lineage_assignment):
                        lineage_null += 1
                    elif lm is True:
                        lineage_hits += 1
                    else:
                        lineage_discrepancies += 1

                if is_meaningful(expected_clade):
                    if not is_meaningful(clade_assignment):
                        clade_null += 1
                    elif cm is True:
                        clade_hits += 1
                    else:
                        clade_discrepancies += 1

            classification_matches_per_lab.append(lab_total_matches)

        for sample_id in comp_expected["samples"].keys():
            sample_lineage_hits = 0
            sample_lineage_discrepancies = 0
            sample_lineage_null = 0
            sample_clade_hits = 0
            sample_clade_discrepancies = 0
            sample_clade_null = 0

            for lab in participating_labs:
                sample = lab["components"][comp_code]["samples"].get(sample_id)
                if not sample:
                    continue

                cls = sample.get("classification", {})
                lm = cls.get("lineage_match")
                cm = cls.get("clade_match")
                expected_lineage = cls.get("expected_lineage")
                lineage_assignment = cls.get("lineage_assignment")
                expected_clade = cls.get("expected_clade")
                clade_assignment = cls.get("clade_assignment")

                if is_meaningful(expected_lineage):
                    if not is_meaningful(lineage_assignment):
                        sample_lineage_null += 1
                    elif lm is True:
                        sample_lineage_hits += 1
                    else:
                        sample_lineage_discrepancies += 1

                if is_meaningful(expected_clade):
                    if not is_meaningful(clade_assignment):
                        sample_clade_null += 1
                    elif cm is True:
                        sample_clade_hits += 1
                    else:
                        sample_clade_discrepancies += 1

            sample_lineage_total = sample_lineage_hits + sample_lineage_discrepancies + sample_lineage_null
            sample_clade_total = sample_clade_hits + sample_clade_discrepancies + sample_clade_null

            per_sample_cls.append({
                "collecting_lab_sample_id": sample_id,
                "lineage_hits": sample_lineage_hits,
                "lineage_total": sample_lineage_total,
                "lineage_discrepancies": sample_lineage_discrepancies,
                "lineage_null": sample_lineage_null,
                "lineage_hit_pct": pct(sample_lineage_hits, sample_lineage_total) if sample_lineage_total else None,
                "clade_hits": sample_clade_hits,
                "clade_total": sample_clade_total,
                "clade_discrepancies": sample_clade_discrepancies,
                "clade_null": sample_clade_null,
                "clade_hit_pct": pct(sample_clade_hits, sample_clade_total) if sample_clade_total else None,
            })

        lineage_total = lineage_hits + lineage_discrepancies + lineage_null
        clade_total = clade_hits + clade_discrepancies + clade_null
        comp_obj["typing"] = {
            "total_classification_matches_median": median_or_none(classification_matches_per_lab),
            "total_classification_matches_min": min_or_none(classification_matches_per_lab),
            "total_classification_matches_max": max_or_none(classification_matches_per_lab),
            "lineage_hits": lineage_hits,
            "lineage_discrepancies": lineage_discrepancies,
            "lineage_null": lineage_null,
            "lineage_total": lineage_total,
            "lineage_hit_pct": pct(lineage_hits, lineage_total),
            "lineage_discrepancies_pct": pct(lineage_discrepancies + lineage_null, lineage_total),
            "clade_hits": clade_hits,
            "clade_discrepancies": clade_discrepancies,
            "clade_null": clade_null,
            "clade_total": clade_total,
            "clade_hit_pct": pct(clade_hits, clade_total),
            "clade_discrepancies_pct": pct(clade_discrepancies + clade_null, clade_total),
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
                "collecting_lab_sample_id": sample_id,
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
            high_and_low_freq = []
            high_freq_only = []
            low_freq_only = []
            n_variants = []
            n_variants_vcf = []
            n_effect = []
            discrepancies_reported = []

            for r in records:
                sample = r["sample"]
                sb = sample.get("software_benchmarking", {})
                var = sample.get("variants", {})

                vp = sb.get("variant_calling_params")
                if is_meaningful(vp):
                    params.append(vp)

                hal = extract_variant_reporting_mode_pct(sample, "high_and_low_freq")
                hfo = extract_variant_reporting_mode_pct(sample, "high_freq_only")
                lfo = extract_variant_reporting_mode_pct(sample, "low_freq_only")
                nvic = safe_number(var.get("number_of_variants_in_consensus"))
                nivcv = safe_number(var.get("number_of_variants_in_consensus_vcf"))
                dirv = safe_number(var.get("discrepancies_in_reported_variants"))

                if hal is not None:
                    high_and_low_freq.append(hal)
                if hfo is not None:
                    high_freq_only.append(hfo)
                if lfo is not None:
                    low_freq_only.append(lfo)
                if nvic is not None:
                    n_variants.append(nvic)
                if nivcv is not None:
                    n_variants_vcf.append(nivcv)
                if dirv is not None:
                    discrepancies_reported.append(dirv)

                nvw = safe_number(var.get("number_of_variants_with_effect"))
                if nvw is not None:
                    n_effect.append(nvw)

            out = {
                "params": most_common_or_none(params),
                "high_and_low_freq_pct": median_or_none(high_and_low_freq),
                "high_freq_only_pct": median_or_none(high_freq_only),
                "low_freq_only_pct": median_or_none(low_freq_only),
                "number_of_variants_in_consensus": median_or_none(n_variants),
                "number_of_variants_in_consensus_vcf": median_or_none(n_variants_vcf),
                "number_of_variants_with_effect": median_or_none(n_effect),
                "discrepancies_in_reported_variants": median_or_none(discrepancies_reported),
            }

            if comp_expected.get("virus") == "SARS-CoV-2":
                n_effect_vcf = []
                discrepancies_reported_effect = []
                successful_hits = []
                total_discrepancies = []

                for r in records:
                    var = r["sample"].get("variants", {})
                    nweev = safe_number(var.get("number_of_variants_with_effect_vcf"))
                    dirve = safe_number(var.get("discrepancies_in_reported_variants_effect"))
                    sh = safe_number(var.get("successful_hits"))
                    td = safe_number(var.get("total_discrepancies"))

                    if nweev is not None:
                        n_effect_vcf.append(nweev)
                    if dirve is not None:
                        discrepancies_reported_effect.append(dirve)
                    if sh is not None:
                        successful_hits.append(sh)
                    if td is not None:
                        total_discrepancies.append(td)

                out["number_of_variants_with_effect_vcf"] = median_or_none(n_effect_vcf)
                out["discrepancies_in_reported_variants_effect"] = median_or_none(discrepancies_reported_effect)
                out["successful_hits"] = median_or_none(successful_hits)
                out["total_discrepancies"] = median_or_none(total_discrepancies)
            else:
                n_variants_vcf_total = []

                for r in records:
                    niv = safe_number(r["sample"].get("variants", {}).get("number_of_variants_in_vcf"))
                    if niv is not None:
                        n_variants_vcf_total.append(niv)

                out["number_of_variants_in_vcf"] = median_or_none(n_variants_vcf_total)

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
            "clade_assignment_software_database_version",
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
                "lineage_assignment_database_version",
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
                "type_assignment_software_database_version",
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
                "subtype_assignment_software_database_version",
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
                "distinct_references": sorted(distinct_sars_references),
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
                "sars_lineage_concordance_pct": pct(sars_lineage_matches, sars_lineage_total),
                "influenza_type_concordance_pct": pct(flu_type_matches, flu_type_total),
                "sars_clade_concordance_pct": pct(sars_clade_matches, sars_clade_total),
                "influenza_clade_concordance_pct": pct(flu_clade_matches, flu_clade_total),
                "sars_cov_2_concordance_pct": pct(sars_lineage_matches, sars_lineage_total),
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
            "variant_calling_params_pct": pct(frequency_threshold_count, frequency_threshold_total),
            "reference_genome_pct": pct(reference_genome_count, reference_genome_total),
            "incomplete_parameters_pct": pct(params_slots_missing, params_slots_total),
            "filled_parameters_pct": pct(params_slots_total - params_slots_missing, params_slots_total),
            "total_workflows": len(all_workflows),
            "total_consensus_softwares": len(consensus_softwares),
            "total_variant_softwares": len(variant_softwares),
            "total_lineage_assignment_softwares": len(lineage_assignment_softwares),
            "total_clade_assignment_softwares": len(clade_assignment_softwares),
            "total_type_assignment_softwares": len(type_assignment_softwares),
            "total_subtype_assignment_softwares": len(subtype_assignment_softwares),
            "total_classification_softwares": len(
                lineage_assignment_softwares
                | clade_assignment_softwares
                | type_assignment_softwares
                | subtype_assignment_softwares
            ),
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
    parser.add_argument("--figures-dir", default="./figures", help="Base directory where all figures will be written")
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

    generate_component_figures(
        general_data=general,
        labs=labs,
        figures_dir=args.figures_dir,
    )

    generate_individual_lab_figures(
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
