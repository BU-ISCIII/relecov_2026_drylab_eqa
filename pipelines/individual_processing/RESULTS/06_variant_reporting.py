#!/usr/bin/env python3

import argparse as _argparse
import json as _json
import sys as _sys
from pathlib import Path as _Path

import pandas as _pd


def _detect_lab_code(base_dir):
    match = None
    for value in (base_dir.name, base_dir.parent.name, str(base_dir)):
        match = __import__("re").search(r"COD-\d+", value)
        if match:
            return match.group(0)
    return None


def _load_expected_variants(expected_data_path):
    with open(expected_data_path, "r", encoding="utf-8") as handle:
        expected_data = _json.load(handle)

    expected_variants = {}
    for component in expected_data.get("components", {}).values():
        for sample_code, sample_data in component.get("samples", {}).items():
            if "expected_variats" in sample_data:
                expected_variants[sample_code] = sample_data["expected_variats"]

    return expected_variants


def _run_csv_only_variant_reporting():
    resultados_map = {
        "Wrong nucleotide": "wrong_nt",
        "Insertion relative to gold standard": "insertions",
        "Deletion relative to gold standard": "deletions",
        "Missing variant": "missing",
        "De novo reported variant": "denovo",
    }

    parser = _argparse.ArgumentParser(
        description="Regenerar variants_report.json contando solo los CSV de variant_analysis_result."
    )
    parser.add_argument("-b", "--base-dir", type=str, default=".", help="Carpeta RESULTS del laboratorio")
    parser.add_argument(
        "--expected-data",
        type=str,
        default="/data/ucct/bi/pipelines/relecov_2026_drylab_eqa/expected_data.json",
        help="JSON con expected_variats por muestra",
    )
    args = parser.parse_args()

    base_dir = _Path(args.base_dir).resolve()
    variant_analysis_dir = base_dir / "variant_analysis_result"
    output_json = base_dir / "variants_report.json"
    lab_code = _detect_lab_code(base_dir)
    expected_variants = _load_expected_variants(_Path(args.expected_data))

    report_payload = {}
    if output_json.exists():
        with open(output_json, "r", encoding="utf-8") as handle:
            try:
                report_payload = _json.load(handle)
            except _json.JSONDecodeError:
                report_payload = {}

    csv_files = [
        path for path in sorted(variant_analysis_dir.iterdir())
        if (
            path.is_file()
            and path.name.startswith("COD")
            and path.name.endswith("_COMBINADO_v2.csv")
        )
    ]

    if not csv_files:
        raise FileNotFoundError(
            f"No COD*_COMBINADO_v2.csv files found in {variant_analysis_dir}"
        )

    print("Reading variant CSV files:")
    for csv_file in csv_files:
        print(f"  - {csv_file.name}")
        df = _pd.read_csv(csv_file)

        required_columns = {"EQA", "RESULTADOS_enviados"}
        missing_columns = required_columns - set(df.columns)
        if missing_columns:
            raise ValueError(f"{csv_file} is missing required columns: {sorted(missing_columns)}")

        if lab_code and "COD_LAB" in df.columns:
            df = df[df["COD_LAB"].astype(str) == lab_code].copy()

        for sample, sample_df in df.groupby("EQA"):
            sample = str(sample)
            report_payload.setdefault(sample, {})

            counts = sample_df["RESULTADOS_enviados"].value_counts(dropna=False).to_dict()
            for result_label, json_key in resultados_map.items():
                report_payload[sample][json_key] = int(counts.get(result_label, 0))

            expected = expected_variants.get(sample)
            wrong_nt = report_payload[sample].get("wrong_nt", 0) or 0
            missing = report_payload[sample].get("missing", 0) or 0
            if expected is not None:
                report_payload[sample]["successful_hits"] = int(expected) - int(wrong_nt) - int(missing)

    with open(output_json, "w", encoding="utf-8") as handle:
        _json.dump(report_payload, handle, indent=4)

    print(f"Updated {output_json}")


if __name__ == "__main__":
    _run_csv_only_variant_reporting()
    _sys.exit(0)

import argparse
import json
import os.path
import re
from pathlib import Path

import numpy as np
import pandas as pd


parser = argparse.ArgumentParser(
    description="Regenerar variants_report.json a partir de los mismos inputs del script de comparación, sin tocar los CSV."
)
parser.add_argument("-b", "--base-dir", type=str, default=".", help="Carpeta base donde buscar")
args = parser.parse_args()


def add_eqa_column(df, samples_list, sample_col="SAMPLE"):
    df = df.copy()
    df["EQA"] = pd.Series(pd.NA, index=df.index, dtype="object")

    for eqaid in samples_list:
        df.loc[df[sample_col].astype(str).str.contains(eqaid, na=False), "EQA"] = eqaid

    if "POS" in df.columns:
        df = df.sort_values(by=["EQA", "POS"], ascending=[True, True])

    return df


def comparar_alt(row, col_ref="REF_gold", col_gold="ALT_gold", col_comp="ALT_comp"):
    ref = str(row[col_ref]) if pd.notna(row[col_ref]) else ""
    alt_gold = str(row[col_gold]) if pd.notna(row[col_gold]) else ""
    alt_comp = str(row[col_comp]) if pd.notna(row[col_comp]) else ""

    has_gold = bool(alt_gold)
    has_comp = bool(alt_comp)

    if has_gold and not has_comp:
        return "Missing variant"

    if has_comp and not has_gold:
        if len(alt_comp) > len(ref):
            return "Insertion relative to gold standard"
        if len(ref) > len(alt_comp):
            return "Deletion relative to gold standard"
        if len(ref) == 1 and len(alt_comp) == 1:
            return "De novo reported variant"
        return "De novo reported variant"

    if has_gold and has_comp:
        if alt_gold == alt_comp:
            return "Match"
        if len(ref) == 1 and len(alt_gold) == 1 and len(alt_comp) == 1:
            return "Wrong nucleotide"
        if len(alt_comp) > len(alt_gold):
            return "Insertion relative to gold standard"
        if len(alt_comp) < len(alt_gold):
            return "Deletion relative to gold standard"
        return "Wrong nucleotide"

    return "Unknown"


def comparar_gold_vs(df_comp, col_suffix, var_gold, extra_cols_gold=None, extra_cols_comp=None):
    if extra_cols_gold is None:
        extra_cols_gold = []
    if extra_cols_comp is None:
        extra_cols_comp = []

    if df_comp.empty:
        merged = var_gold.copy()
        merged["SAMPLE_ID"] = pd.NA
        for col in [f"ALT_{col_suffix}", f"CHROM_{col_suffix}"] + extra_cols_comp:
            merged[col] = pd.NA
    else:
        eqa_validos = df_comp["EQA"].dropna().unique()

        gold_filtered = var_gold[var_gold["EQA"].isin(eqa_validos)].copy()
        df_comp = df_comp.copy()

        gold_filtered = gold_filtered.rename(
            columns={
                "SAMPLE": "SAMPLE_gold",
                "CHROM": "CHROM_gold",
                "REF": "REF_gold",
                "ALT": "ALT_gold",
            }
        )

        df_comp = df_comp.rename(
            columns={
                "SAMPLE": f"SAMPLE_{col_suffix}",
                "CHROM": f"CHROM_{col_suffix}",
                "REF": f"REF_{col_suffix}",
                "ALT": f"ALT_{col_suffix}",
            }
        )

        merged = gold_filtered.merge(
            df_comp,
            left_on=["EQA", "POS", "REF_gold"],
            right_on=["EQA", "POS", f"REF_{col_suffix}"],
            how="outer",
            indicator=True,
        )

        merged["REF_gold"] = merged["REF_gold"].combine_first(merged[f"REF_{col_suffix}"])
        merged["SAMPLE_ID"] = merged.get(f"SAMPLE_{col_suffix}").combine_first(merged.get("SAMPLE_gold"))

    diff = merged[
        (merged["_merge"] != "both") | (merged["ALT_gold"] != merged.get(f"ALT_{col_suffix}", pd.NA))
    ].copy()

    for col in [f"ALT_{col_suffix}", f"CHROM_{col_suffix}"]:
        if col not in diff.columns:
            diff[col] = pd.NA

    for col in extra_cols_gold:
        col_gold = f"{col}_gold"
        if col_gold not in diff.columns:
            diff[col_gold] = diff[col] if col in diff.columns else pd.NA

    for col in extra_cols_comp:
        col_comp = f"{col}_{col_suffix}"
        if col_comp not in diff.columns:
            diff[col_comp] = diff[col] if col in diff.columns else pd.NA

    diff["POS"] = pd.to_numeric(diff["POS"], errors="coerce").astype("Int64")

    diff_final = diff[
        ["SAMPLE_ID", "EQA", "POS", "REF_gold", "ALT_gold", f"ALT_{col_suffix}", f"CHROM_{col_suffix}"]
        + [f"{col}_gold" for col in extra_cols_gold]
        + [f"{col}_{col_suffix}" for col in extra_cols_comp]
    ].copy()

    diff_final[f"RESULTADOS_{col_suffix}"] = diff_final.apply(
        comparar_alt,
        axis=1,
        col_ref="REF_gold",
        col_gold="ALT_gold",
        col_comp=f"ALT_{col_suffix}",
    )

    return diff_final.sort_values(by=["EQA", "POS"])


def calculate_values_eqa(merged_all, vlt_lab, vlt_gold, expected_samples, output_json):
    variants_dict = {}
    resultados_map = {
        "Wrong nucleotide": "wrong_nt",
        "Insertion relative to gold standard": "insertions",
        "Deletion relative to gold standard": "deletions",
        "Missing variant": "missing",
        "De novo reported variant": "denovo",
    }

    df = merged_all.copy()
    vlt_df = vlt_lab.copy()
    vlt_gold_df = vlt_gold.copy()
    af_threshold = 0.75

    lab_samples = set(vlt_df["SAMPLE"].dropna().astype(str)) if not vlt_df.empty else set()
    gold_samples = set(vlt_gold_df["SAMPLE"].dropna().astype(str)) if not vlt_gold_df.empty else set()

    for sample in map(str, expected_samples):
        variants_dict[sample] = {}

        has_vcf = sample in lab_samples
        has_gold = sample in gold_samples
        is_calculable = has_vcf and has_gold

        variants_dict[sample]["successful_hits"] = None
        variants_dict[sample]["high_and_low_freq"] = None
        variants_dict[sample]["high_freq_only"] = None
        variants_dict[sample]["low_freq_only"] = None
        variants_dict[sample]["number_of_variants_in_consensus_vcf"] = None
        variants_dict[sample]["number_of_variants_with_effect_vcf"] = None

        for json_equivalent in resultados_map.values():
            variants_dict[sample][json_equivalent] = 0 if is_calculable else None

    for sample in map(str, expected_samples):
        has_vcf = sample in lab_samples
        has_gold = sample in gold_samples
        if not (has_vcf and has_gold):
            continue

        vlt_df_sample = vlt_df[vlt_df["SAMPLE"].astype(str) == sample]
        gold_sample_df = vlt_gold_df[vlt_gold_df["SAMPLE"].astype(str) == sample]

        set1 = set(map(tuple, vlt_df_sample[["POS", "REF", "ALT"]].values))
        set2 = set(map(tuple, gold_sample_df[["POS", "REF", "ALT"]].values))
        variants_dict[sample]["successful_hits"] = len(set1 & set2)

    for sample in map(str, expected_samples):
        if sample not in lab_samples:
            continue

        group = vlt_df[vlt_df["SAMPLE"].astype(str) == sample].copy()
        af_values = pd.to_numeric(group["AF"], errors="coerce")

        has_low_freq = (af_values < af_threshold).any()
        has_high_freq = (af_values > af_threshold).any()

        variants_dict[sample]["high_and_low_freq"] = bool(has_high_freq and has_low_freq)
        variants_dict[sample]["high_freq_only"] = bool(has_high_freq and not has_low_freq)
        variants_dict[sample]["low_freq_only"] = bool(has_low_freq and not has_high_freq)

        group_high_af = group[af_values > af_threshold].copy()
        variants_dict[sample]["number_of_variants_in_consensus_vcf"] = len(group_high_af)

        df_variants_effect = group_high_af[group_high_af["EFFECT"] == "missense_variant"]
        variants_dict[sample]["number_of_variants_with_effect_vcf"] = len(df_variants_effect)

    if not df.empty:
        for sample, reported_df in df.groupby("EQA"):
            sample = str(sample)
            if sample not in variants_dict:
                continue

            has_vcf = sample in lab_samples
            has_gold = sample in gold_samples
            if not (has_vcf and has_gold):
                continue

            resultados_counts = reported_df["RESULTADOS_enviados"].value_counts(dropna=False).to_dict()
            for key, json_equivalent in resultados_map.items():
                variants_dict[sample][json_equivalent] = int(resultados_counts.get(key, 0))

    if os.path.isfile(output_json):
        try:
            with open(output_json, "r", encoding="utf-8") as handle:
                existing = json.load(handle)
            if existing:
                existing.update(variants_dict)
                variants_dict = existing
        except json.JSONDecodeError:
            pass

    with open(output_json, "w", encoding="utf-8") as handle:
        json.dump(variants_dict, handle, indent=4)


def get_pattern_vlt(gold_standard_filename):
    organism = "Influenza" if "FLU" in gold_standard_filename else "coronavirus_2"
    tech = "Illumina" if "1" in gold_standard_filename else "Nanopore"
    return f"*{organism}*{tech}*csv"


def main():
    base_dir = Path(args.base_dir).resolve()
    cod_dir = Path(f"{args.base_dir}/../").resolve().name

    vcf_dir = base_dir / "variants_long_table"
    gold_dir = base_dir / "gold_standard/variants_long_table"
    output_json = base_dir / "variants_report.json"

    all_variants = {}

    for gold_standard in gold_dir.rglob("*.csv"):
        if "FLU" in gold_standard.name:
            continue

        var_gold = pd.read_csv(gold_standard)
        if "DP" in var_gold.columns:
            var_gold.rename(columns={"DP": "DP_gold"}, inplace=True)
        if "AF" in var_gold.columns:
            var_gold.rename(columns={"AF": "AF_gold"}, inplace=True)

        var_enviados = pd.DataFrame()
        samples = list(set(var_gold.SAMPLE.values))
        variant_long_table_pattern = get_pattern_vlt(gold_standard.name)

        for file in vcf_dir.rglob(variant_long_table_pattern):
            if "full" in file.name:
                continue
            var_enviados = pd.read_csv(file)

        if var_enviados.empty:
            continue

        var_gold = add_eqa_column(var_gold, samples)
        var_enviados = add_eqa_column(var_enviados, samples)

        for df in (var_gold, var_enviados):
            if not df.empty:
                for col in ["EQA", "POS", "ALT", "REF", "CHROM"]:
                    if col in df.columns:
                        df[col] = df[col].astype(str)

        diff_enviados = comparar_gold_vs(
            var_enviados,
            "enviados",
            var_gold,
            extra_cols_gold=["DP", "AF"],
            extra_cols_comp=["DP", "AF"],
        )

        merged_all = diff_enviados
        match = re.search(r"(COD-\d+)", cod_dir)
        sample_tag = match.group(1) if match else base_dir.name
        merged_all["COD_LAB"] = sample_tag

        final_cols = [
            "COD_LAB",
            "EQA",
            "POS",
            "REF_gold",
            "ALT_gold",
            "ALT_enviados",
            "DP_gold",
            "DP_enviados",
            "AF_gold",
            "AF_enviados",
            "RESULTADOS_enviados",
        ]
        for col in final_cols:
            if col not in merged_all.columns:
                merged_all[col] = pd.NA

        merged_all = merged_all[final_cols]
        merged_all.rename(columns={"REF_gold": "REF (NC_045512.2)", "ALT_gold": "Gold_Standard"}, inplace=True)
        merged_all = merged_all.sort_values(by=["EQA", "POS"])

        tmp_output = base_dir / ".variants_report.partial.json"
        calculate_values_eqa(merged_all, var_enviados, var_gold, samples, tmp_output)

        if tmp_output.exists():
            with open(tmp_output, "r", encoding="utf-8") as handle:
                partial = json.load(handle)
            all_variants.update(partial)
            tmp_output.unlink()

    with open(output_json, "w", encoding="utf-8") as handle:
        json.dump(all_variants, handle, indent=4)


if __name__ == "__main__":
    main()
