#!/usr/bin/env python3

import argparse
import re
from pathlib import Path

import pandas as pd


COMPONENT_SAMPLES = {
    "SARS1": [f"SARS{i}" for i in range(1, 6)],
    "SARS2": [f"SARS{i}" for i in range(6, 11)],
}

FINAL_COLUMNS = [
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


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Compare submitted variant calls against variants_long_table gold standards. "
            "Supports both pre-generated variants_long_table CSVs and per-sample VCF/iVar files."
        )
    )
    parser.add_argument(
        "--mode",
        choices=("auto", "long-table", "input-pair", "per-sample-vcf"),
        default="auto",
        help=(
            "Comparison mode. auto uses input-pair with two input files, per-sample-vcf when "
            "--component is provided, otherwise long-table."
        ),
    )
    parser.add_argument(
        "input_files",
        nargs="*",
        help="Optional explicit inputs: gold variants_long_table followed by submitted long-table, VCF, or iVar TSV",
    )
    parser.add_argument("-b", "--base-dir", default=".", help="RESULTS directory for long-table mode")
    parser.add_argument("-o", "--output-dir", default=None, help="Output directory for COMBINADO_v2 CSVs")
    parser.add_argument(
        "-c",
        "--component",
        choices=sorted(COMPONENT_SAMPLES),
        help="SARS component for per-sample VCF mode: SARS1 for SARS1-SARS5, SARS2 for SARS6-SARS10",
    )
    parser.add_argument(
        "-r",
        "--raw-dir",
        default="../RAW",
        help="Directory containing per-sample VCF/iVar files for per-sample VCF mode",
    )
    parser.add_argument(
        "-g",
        "--gold-dir",
        default="gold_standard/variants_long_table",
        help="Directory containing variants_long_table gold CSVs",
    )
    parser.add_argument("--gold-file", help="Explicit gold variants_long_table CSV for per-sample VCF mode")
    parser.add_argument(
        "--sample",
        help="Optional fallback only for a single VCF/iVar file whose filename does not contain SARSx/FLUx",
    )
    parser.add_argument(
        "--vcf-template",
        default="{sample}.medaka.vcf",
        help="Per-sample VCF filename template relative to --raw-dir",
    )
    return parser.parse_args()


def cod_tag_from_path(path: Path) -> str:
    match = re.search(r"(COD-\d+)", str(path.resolve()))
    return match.group(1) if match else path.resolve().name


def info_value(info: str, key: str):
    if not info or info == ".":
        return pd.NA
    for item in str(info).split(";"):
        if item == key:
            return True
        if item.startswith(f"{key}="):
            return item.split("=", 1)[1]
    return pd.NA


def format_value(format_keys: str, sample_values: str, key: str):
    if not format_keys or not sample_values or format_keys == "." or sample_values == ".":
        return pd.NA
    keys = str(format_keys).split(":")
    values = str(sample_values).split(":")
    return dict(zip(keys, values)).get(key, pd.NA)


def first_available(*values):
    for value in values:
        if pd.notna(value) and value != "":
            return value
    return pd.NA


def normalize_ivar_ref_alt(ref, alt):
    ref = str(ref)
    alt = str(alt)
    if alt.startswith("+"):
        return ref, f"{ref}{alt[1:]}"
    if alt.startswith("-"):
        return f"{ref}{alt[1:]}", ref
    return ref, alt


def parse_ivar_tsv(tsv_path: Path, sample: str) -> pd.DataFrame:
    df = pd.read_csv(tsv_path, sep="\t")
    required_columns = {"REGION", "POS", "REF", "ALT"}
    missing_columns = required_columns - set(df.columns)
    if missing_columns:
        raise ValueError(f"{tsv_path} is missing required iVar TSV columns: {sorted(missing_columns)}")

    rows = []
    for _, row in df.iterrows():
        ref, alt = normalize_ivar_ref_alt(row["REF"], row["ALT"])
        rows.append(
            {
                "SAMPLE": sample,
                "EQA": sample,
                "CHROM": row.get("REGION", pd.NA),
                "POS": row.get("POS", pd.NA),
                "REF": ref,
                "ALT": alt,
                "DP": first_available(row.get("TOTAL_DP", pd.NA), row.get("ALT_DP", pd.NA)),
                "AF": row.get("ALT_FREQ", pd.NA),
            }
        )
    return pd.DataFrame(rows)


def is_ivar_tsv(path: Path) -> bool:
    with open(path, "r", encoding="utf-8") as handle:
        first_line = handle.readline().strip()
    return first_line.startswith("REGION\tPOS\tREF\tALT")


def first_nonempty_line(path: Path) -> str:
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if line:
                return line
    return ""


def is_vcf(path: Path) -> bool:
    line = first_nonempty_line(path)
    return line.startswith("##fileformat=VCF") or line.startswith("#CHROM")


def read_table_auto(path: Path) -> pd.DataFrame:
    if path.suffix.lower() in {".tsv", ".tab"}:
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path, sep=None, engine="python")


def looks_like_long_table(path: Path) -> bool:
    try:
        df = read_table_auto(path)
    except Exception:
        return False
    return {"SAMPLE", "POS", "REF", "ALT"}.issubset(df.columns)


def detect_sample_from_path(path: Path) -> str | None:
    match = re.search(r"(SARS\d+|FLU\d+)", path.name, re.IGNORECASE)
    return match.group(1).upper() if match else None


def parse_vcf(vcf_path: Path, sample: str) -> pd.DataFrame:
    if is_ivar_tsv(vcf_path):
        return parse_ivar_tsv(vcf_path, sample)

    rows = []
    header = None
    with open(vcf_path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.lstrip("#").split("\t")
                continue
            if line.startswith("#"):
                continue
            if header is None:
                raise ValueError(f"VCF header not found before variants in {vcf_path}")

            record = dict(zip(header, line.split("\t")))
            info = record.get("INFO", "")
            fmt = record.get("FORMAT", "")
            sample_col = header[-1] if len(header) > 8 else None
            sample_values = record.get(sample_col, "") if sample_col else ""
            dp = first_available(
                info_value(info, "DP"),
                format_value(fmt, sample_values, "DP"),
                format_value(fmt, sample_values, "ALT_DP"),
            )
            af = first_available(
                info_value(info, "AF"),
                format_value(fmt, sample_values, "AF"),
                format_value(fmt, sample_values, "ALT_FREQ"),
            )

            for alt in str(record.get("ALT", "")).split(","):
                rows.append(
                    {
                        "SAMPLE": sample,
                        "EQA": sample,
                        "CHROM": record.get("CHROM", pd.NA),
                        "POS": record.get("POS", pd.NA),
                        "REF": record.get("REF", pd.NA),
                        "ALT": alt,
                        "DP": dp,
                        "AF": af,
                    }
                )
    return pd.DataFrame(rows)


def add_eqa_column(df: pd.DataFrame, samples: list[str], sample_col="SAMPLE") -> pd.DataFrame:
    df = df.copy()
    df["EQA"] = pd.Series(pd.NA, index=df.index, dtype="object")
    sample_values = df[sample_col].astype(str)
    for sample in samples:
        # Avoid matching SARS1 inside SARS10 while keeping historical substring behaviour.
        df.loc[sample_values.str.contains(rf"{re.escape(sample)}(?!\d)", regex=True, na=False), "EQA"] = sample
    if "POS" in df.columns:
        df = df.sort_values(by=["EQA", "POS"], ascending=[True, True])
    return df


def compare_alt(row, col_ref="REF_gold", col_gold="ALT_gold", col_comp="ALT_enviados"):
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


def compare_gold_vs_calls(gold: pd.DataFrame, calls: pd.DataFrame) -> pd.DataFrame:
    gold = gold.copy()
    calls = calls.copy()

    for df in (gold, calls):
        for col in ["EQA", "POS", "REF", "ALT", "CHROM"]:
            if col in df.columns:
                df[col] = df[col].astype(str)

    gold = gold.rename(
        columns={
            "SAMPLE": "SAMPLE_gold",
            "CHROM": "CHROM_gold",
            "REF": "REF_gold",
            "ALT": "ALT_gold",
            "DP": "DP_gold",
            "AF": "AF_gold",
        }
    )
    calls = calls.rename(
        columns={
            "SAMPLE": "SAMPLE_enviados",
            "CHROM": "CHROM_enviados",
            "REF": "REF_enviados",
            "ALT": "ALT_enviados",
            "DP": "DP_enviados",
            "AF": "AF_enviados",
        }
    )

    merged = gold.merge(
        calls,
        left_on=["EQA", "POS", "REF_gold"],
        right_on=["EQA", "POS", "REF_enviados"],
        how="outer",
        indicator=True,
    )
    merged["REF_gold"] = merged["REF_gold"].combine_first(merged["REF_enviados"])

    diff = merged[(merged["_merge"] != "both") | (merged["ALT_gold"] != merged["ALT_enviados"])].copy()
    for col in ["ALT_gold", "ALT_enviados", "DP_gold", "DP_enviados", "AF_gold", "AF_enviados"]:
        if col not in diff.columns:
            diff[col] = pd.NA

    diff["POS"] = pd.to_numeric(diff["POS"], errors="coerce").astype("Int64")
    diff["RESULTADOS_enviados"] = diff.apply(compare_alt, axis=1)
    return diff.sort_values(by=["EQA", "POS"])


def finalize_and_write(comparison: pd.DataFrame, output_file: Path, cod_lab: str):
    comparison = comparison.copy()
    comparison["COD_LAB"] = cod_lab
    for col in FINAL_COLUMNS:
        if col not in comparison.columns:
            comparison[col] = pd.NA
    comparison = comparison[FINAL_COLUMNS].rename(
        columns={"REF_gold": "REF (NC_045512.2)", "ALT_gold": "Gold_Standard"}
    )
    output_file.parent.mkdir(parents=True, exist_ok=True)
    comparison.to_csv(output_file, index=False)
    print(f"Output written: {output_file}")


def find_gold_file(component: str, gold_dir: Path, explicit_gold_file: str | None) -> Path:
    if explicit_gold_file:
        gold_file = Path(explicit_gold_file)
        if not gold_file.exists():
            raise FileNotFoundError(f"Gold file not found: {gold_file}")
        return gold_file

    candidates = sorted(gold_dir.glob(f"*{component}*variant_long_table*.csv"))
    if not candidates:
        raise FileNotFoundError(f"No gold variants_long_table CSV found for {component} in {gold_dir}")
    if len(candidates) > 1:
        print(f"WARNING: multiple gold files found for {component}; using {candidates[0]}")
    return candidates[0]


def load_component_vcfs(raw_dir: Path, samples: list[str], vcf_template: str) -> pd.DataFrame:
    dfs = []
    for sample in samples:
        vcf_path = raw_dir / vcf_template.format(sample=sample)
        if not vcf_path.exists():
            fallback = None
            for pattern in [f"{sample}_variants.vcf", f"*{sample}.vcf", f"*{sample}*.vcf"]:
                matches = sorted(raw_dir.glob(pattern))
                if matches:
                    fallback = matches[0]
                    break
            if fallback is None:
                print(f"WARNING: VCF not found for {sample}: {vcf_path}")
                continue
            print(f"WARNING: {vcf_path.name} not found; using fallback {fallback.name}")
            vcf_path = fallback
        dfs.append(parse_vcf(vcf_path, sample))

    if not dfs:
        return pd.DataFrame(columns=["SAMPLE", "EQA", "CHROM", "POS", "REF", "ALT", "DP", "AF"])
    return pd.concat(dfs, ignore_index=True)


def load_gold_for_samples(gold_file: Path, samples: list[str]) -> pd.DataFrame:
    gold = pd.read_csv(gold_file)
    if "SAMPLE" not in gold.columns:
        raise ValueError(f"Gold file is missing SAMPLE column: {gold_file}")
    gold = gold.copy()
    gold["SAMPLE"] = gold["SAMPLE"].astype(str)
    gold = gold[gold["SAMPLE"].isin(samples)].copy()
    gold["EQA"] = gold["SAMPLE"]
    return gold


def resolve_path(path_value: str, base_dir: Path) -> Path:
    path = Path(path_value)
    return path.resolve() if path.is_absolute() else (base_dir / path).resolve()


def load_gold_long_table(gold_file: Path) -> pd.DataFrame:
    gold = read_table_auto(gold_file)
    if "SAMPLE" not in gold.columns:
        raise ValueError(f"Gold file is missing SAMPLE column: {gold_file}")
    return gold


def prepare_gold_for_comparison(gold: pd.DataFrame, samples: list[str]) -> pd.DataFrame:
    gold = gold.copy()
    gold["SAMPLE"] = gold["SAMPLE"].astype(str)
    if samples:
        gold = gold[gold["SAMPLE"].isin(samples)].copy()
    return add_eqa_column(gold, sorted(set(gold["SAMPLE"].dropna().astype(str))))


def load_submitted_input(
    submitted_path: Path,
    gold_samples: list[str],
    sample_arg: str | None,
    vcf_template: str,
) -> tuple[pd.DataFrame, list[str]]:
    if submitted_path.is_dir():
        calls = load_component_vcfs(submitted_path, gold_samples, vcf_template)
        samples = sorted(set(calls["EQA"].dropna().astype(str)))
        return calls, samples

    if looks_like_long_table(submitted_path):
        calls = read_table_auto(submitted_path)
        calls = add_eqa_column(calls, gold_samples)
        samples = sorted(set(calls["EQA"].dropna().astype(str)))
        return calls, samples

    if is_ivar_tsv(submitted_path) or is_vcf(submitted_path):
        sample = sample_arg or detect_sample_from_path(submitted_path)
        if not sample:
            raise ValueError(
                "Could not infer sample from this single VCF/iVar filename. "
                "Use a filename containing SARSx/FLUx, pass a directory to process all samples, "
                "or use --sample only for this exceptional single-file case."
            )
        calls = parse_vcf(submitted_path, sample)
        return calls, [sample]

    raise ValueError(
        f"Could not detect submitted input format for {submitted_path}. "
        "Expected a variants_long_table CSV/TSV, VCF, iVar TSV, or a directory of VCF/iVar files."
    )


def run_input_pair_mode(args):
    if len(args.input_files) != 2:
        raise ValueError("input-pair mode requires exactly two files: gold_long_table submitted_calls")

    base_dir = Path(args.base_dir).resolve()
    output_dir = resolve_path(args.output_dir, base_dir) if args.output_dir else base_dir / "variant_analysis_result"
    gold_file = resolve_path(args.input_files[0], base_dir)
    submitted_file = resolve_path(args.input_files[1], base_dir)

    gold_raw = load_gold_long_table(gold_file)
    gold_samples = sorted(set(gold_raw["SAMPLE"].dropna().astype(str)))
    if args.component:
        gold_samples = [sample for sample in COMPONENT_SAMPLES[args.component] if sample in gold_samples]

    calls, submitted_samples = load_submitted_input(
        submitted_file,
        gold_samples,
        args.sample,
        args.vcf_template,
    )
    gold = prepare_gold_for_comparison(gold_raw, submitted_samples)

    comparison = compare_gold_vs_calls(gold, calls)
    cod_lab = cod_tag_from_path(output_dir)
    output_file = output_dir / f"{cod_lab}_diferencias_{gold_file.name}_COMBINADO_v2.csv"
    finalize_and_write(comparison, output_file, cod_lab)


def run_per_sample_vcf_mode(args):
    if not args.component:
        raise ValueError("--component is required for per-sample-vcf mode")

    base_dir = Path(args.base_dir).resolve()
    raw_dir = resolve_path(args.raw_dir, base_dir)
    gold_dir = resolve_path(args.gold_dir, base_dir)
    output_dir = resolve_path(args.output_dir, base_dir) if args.output_dir else base_dir / "variant_analysis_result"
    samples = COMPONENT_SAMPLES[args.component]

    gold_file = find_gold_file(args.component, gold_dir, args.gold_file)
    gold = load_gold_for_samples(gold_file, samples)
    calls = load_component_vcfs(raw_dir, samples, args.vcf_template)
    comparison = compare_gold_vs_calls(gold, calls)
    cod_lab = cod_tag_from_path(output_dir)
    output_file = output_dir / f"{cod_lab}_diferencias_{gold_file.name}_COMBINADO_v2.csv"
    finalize_and_write(comparison, output_file, cod_lab)


def get_pattern_vlt(gold_standard_filename: str) -> str:
    organism = "Influenza" if "FLU" in gold_standard_filename else "coronavirus_2"
    tech = "Illumina" if "1" in gold_standard_filename else "Nanopore"
    return f"*{organism}*{tech}*csv"


def run_long_table_mode(args):
    base_dir = Path(args.base_dir).resolve()
    output_dir = resolve_path(args.output_dir, base_dir) if args.output_dir else base_dir / "variant_analysis_result"
    vcf_dir = base_dir / "variants_long_table"
    gold_dir = resolve_path(args.gold_dir, base_dir) if args.gold_dir else base_dir / "gold_standard/variants_long_table"
    cod_lab = cod_tag_from_path(base_dir.parent)

    output_dir.mkdir(parents=True, exist_ok=True)
    written = 0

    for gold_file in sorted(gold_dir.rglob("*.csv")):
        if "FLU" in gold_file.name:
            continue

        gold = pd.read_csv(gold_file)
        if "SAMPLE" not in gold.columns:
            print(f"WARNING: skipping {gold_file}, missing SAMPLE column")
            continue

        samples = sorted(set(gold["SAMPLE"].dropna().astype(str)))
        pattern = get_pattern_vlt(gold_file.name)
        submitted_files = [path for path in sorted(vcf_dir.rglob(pattern)) if "full" not in path.name]
        if not submitted_files:
            print(f"WARNING: no submitted variants_long_table found for {gold_file.name} with pattern {pattern}")
            continue
        if len(submitted_files) > 1:
            print(f"WARNING: multiple submitted long-table files for {gold_file.name}; using {submitted_files[0]}")

        calls = pd.read_csv(submitted_files[0])
        gold = add_eqa_column(gold, samples)
        calls = add_eqa_column(calls, samples)

        comparison = compare_gold_vs_calls(gold, calls)
        output_file = output_dir / f"{cod_lab}_diferencias_{gold_file.name}_COMBINADO_v2.csv"
        finalize_and_write(comparison, output_file, cod_lab)
        written += 1

    if written == 0:
        raise FileNotFoundError("No variant comparison CSVs were generated")


def main():
    args = parse_args()
    mode = args.mode
    if mode == "auto":
        if len(args.input_files) == 2:
            mode = "input-pair"
        elif args.component:
            mode = "per-sample-vcf"
        else:
            mode = "long-table"

    if mode == "input-pair":
        run_input_pair_mode(args)
    elif mode == "per-sample-vcf":
        run_per_sample_vcf_mode(args)
    else:
        run_long_table_mode(args)


if __name__ == "__main__":
    main()
