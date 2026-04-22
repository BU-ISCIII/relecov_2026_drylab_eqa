#!/usr/bin/env python3

from Bio import SeqIO
import csv
import re
from pathlib import Path
from collections import defaultdict


BASE_DIR = Path.cwd()

SARS_DIR = BASE_DIR / "ALIGNMENTS_SARS"
FLU_DIR = BASE_DIR / "ALIGNMENTS_FLU"
FASTA_ANALYSIS = BASE_DIR / "FASTA_ANALYSIS_REPORTS"
INDIVIDUAL_REPORTS = FASTA_ANALYSIS

FASTA_ANALYSIS.mkdir(exist_ok=True)
INDIVIDUAL_REPORTS.mkdir(exist_ok=True)


AMBIGUITY_MAP = {
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"G", "C"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"},
    "B": {"C", "G", "T"},
    "D": {"A", "G", "T"},
    "H": {"A", "C", "T"},
    "V": {"A", "C", "G"},
}


def base_category(base):
    base = base.upper()
    if base == "N":
        return "N"
    if base == "-":
        return "gap"
    return "nt"


def flu_bases_are_equivalent(gold_base, sample_base):
    gold = gold_base.upper()
    sample = sample_base.upper()

    if gold == sample:
        return True

    if gold in AMBIGUITY_MAP:
        if sample in AMBIGUITY_MAP:
            return True
        if sample in AMBIGUITY_MAP[gold]:
            return True

    return False


def clasificar(ref_bases, sample_bases):
    r = ref_bases.upper()
    s = sample_bases.upper()
    nts = set("ATGC")
    ambigs = set("RYSWKMBDHV")

    if r == s:
        return ""

    if len(r) == 1:
        if s == "N":
            return "Stretch of Ns instead of nucleotide"
        if s == "-":
            return "Deletion relative to gold standard"
        if r == "-":
            return "Insertion relative to gold standard"
        if r == "N":
            return "Nucleotide stretch instead of stretch of Ns"
        if s in nts and r in nts:
            return "Wrong nucleotide"
        if s in ambigs and r in nts:
            return "Ambiguity instead of nucleotide"
        if r in ambigs and s in nts:
            return "Nucleotide instead of ambiguity"
        return "Other: Review"

    if set(r) <= {"-"} and set(s) <= {"N"}:
        return ""
    if set(r) <= {"N"} and set(s) <= {"-"}:
        return ""
    if "-" in s and any(base in "ATGCRYSWKMBDHV" for base in r) and "N" not in r:
        return "Deletion relative to gold standard"
    if ("-" not in s and "N" not in s) and "-" in r:
        return "Insertion relative to gold standard"
    if "N" in s and not any(x in r for x in ["N", "-"]):
        return "Stretch of Ns instead of nucleotide"
    if "N" in r and not any(x in s for x in ["N", "-"]):
        return "Nucleotide stretch instead of stretch of Ns"
    if "N" in r and "-" in s:
        return "Deletion relative to gold standard"
    if all(base in nts for base in r) and all(base in nts for base in s):
        return "Wrong nucleotide"
    if all(base in nts for base in r) and all(base in ambigs for base in s):
        return "Ambiguity instead of nucleotide"
    if all(base in ambigs for base in r) and all(base in nts for base in s):
        return "Nucleotide instead of ambiguity"
    return "Other: Review"


def append_grouped_result(output_rows, cod, eqa, sample_id, positions, start, end, ref_str, sample_str, result):
    if not result:
        return

    pos_str = positions[start - 1] if start == end else f"{positions[start - 1]}-{positions[end - 1]}"

    if (
        output_rows
        and output_rows[-1][0] == cod
        and output_rows[-1][1] == eqa
        and output_rows[-1][2] == sample_id
        and output_rows[-1][6] == result
        and start > 1
    ):
        prev_start, prev_end = output_rows[-1][3].split("-") if "-" in output_rows[-1][3] else (output_rows[-1][3], output_rows[-1][3])
        previous_position = positions[start - 2]
        if prev_end == previous_position:
            merged_pos = f"{prev_start}-{positions[end - 1]}" if prev_start != positions[end - 1] else prev_start
            output_rows[-1][3] = merged_pos
            output_rows[-1][4] += ref_str
            output_rows[-1][5] += sample_str
            return

    output_rows.append([cod, eqa, sample_id, pos_str, ref_str, sample_str, result])


def sort_key(row):
    cod_num = int(re.search(r"(\d+)", row[0]).group(1)) if re.search(r"(\d+)", row[0]) else 99999
    eqa_num = re.search(r"(\d+)", row[1])
    eqa_key = int(eqa_num.group(1)) if eqa_num else 99999
    return (cod_num, eqa_key)


def write_rows_by_lab(rows_by_lab, filename_template, headers):
    for cod_lab, rows in rows_by_lab.items():
        filename = INDIVIDUAL_REPORTS / filename_template.format(cod_lab=cod_lab)
        with open(filename, "w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(headers)
            writer.writerows(rows)


def process_sars():
    output_rows = []

    for fasta_file in SARS_DIR.glob("*_aligned.fasta"):
        filename = fasta_file.name
        eqa_match = re.search(r"(EQA[_-]SARS\d+)", filename, re.IGNORECASE)
        eqa = eqa_match.group(1) if eqa_match else "NA"

        records = list(SeqIO.parse(fasta_file, "fasta"))
        if len(records) < 4:
            continue

        ref_seq = str(records[0].seq)
        first_nt = next((i for i, b in enumerate(ref_seq) if b.upper() in "ATGCRYKMSWBDHV"), None)
        last_nt = next((i for i, b in reversed(list(enumerate(ref_seq))) if b.upper() in "ATGCRYKMSWBDHV"), None)

        ref_seq_positions = []
        ref_counter = 0
        for i, base in enumerate(ref_seq):
            if first_nt <= i <= last_nt:
                if base.upper() in "ATGCRYKMSWBDHV":
                    ref_counter += 1
                    ref_seq_positions.append(str(ref_counter))
                else:
                    ref_seq_positions.append(str(ref_counter))
            else:
                ref_seq_positions.append("NA")

        gold_seq = str(records[2].seq)

        for record in records[3:]:
            sample_id = record.id
            seq = str(record.seq)
            cod_match = re.search(r"COD[_-]?(\d+)", sample_id, re.IGNORECASE)
            cod = "COD-" + cod_match.group(1) if cod_match else "NA"

            diffs = [
                (i + 1, r.upper(), s.upper())
                for i, (r, s) in enumerate(zip(gold_seq, seq))
                if first_nt <= i <= last_nt and r.upper() != s.upper()
            ]
            if not diffs:
                continue

            start = diffs[0][0]
            ref_group = [diffs[0][1]]
            sample_group = [diffs[0][2]]

            for j in range(1, len(diffs)):
                pos, r_base, s_base = diffs[j]
                prev_pos, prev_r, prev_s = diffs[j - 1]

                if (
                    pos == prev_pos + 1
                    and base_category(r_base) == base_category(prev_r)
                    and base_category(s_base) == base_category(prev_s)
                ):
                    ref_group.append(r_base)
                    sample_group.append(s_base)
                else:
                    end = prev_pos
                    ref_str = "".join(ref_group)
                    sample_str = "".join(sample_group)
                    result = clasificar(ref_str, sample_str)
                    append_grouped_result(output_rows, cod, eqa, sample_id, ref_seq_positions, start, end, ref_str, sample_str, result)
                    start = pos
                    ref_group = [r_base]
                    sample_group = [s_base]

            end = diffs[-1][0]
            ref_str = "".join(ref_group)
            sample_str = "".join(sample_group)
            result = clasificar(ref_str, sample_str)
            append_grouped_result(output_rows, cod, eqa, sample_id, ref_seq_positions, start, end, ref_str, sample_str, result)

    output_rows.sort(key=sort_key)
    rows_by_lab = defaultdict(list)
    for row in output_rows:
        rows_by_lab[row[0]].append(row)

    write_rows_by_lab(
        rows_by_lab,
        "{cod_lab}_informe_SARS-CoV-2_2026_fasta_analysis.csv",
        ["COD_LAB", "EQA", "Sample_ID", "POS (NC_045512.2)", "Gold_LOW", "Sample", "Resultado"],
    )


def process_flu():
    output_rows = []

    for fasta_file in FLU_DIR.glob("*_aligned.fasta"):
        filename = fasta_file.name
        eqa_match = re.search(r"(EQA[_-]FLU\d+)", filename, re.IGNORECASE)
        eqa = eqa_match.group(1) if eqa_match else "NA"

        records = list(SeqIO.parse(fasta_file, "fasta"))
        if len(records) < 2:
            continue

        gold_record = next((r for r in records if "Gold_Standard" in r.id), None)
        if gold_record is None:
            continue

        gold_seq = str(gold_record.seq)
        first_nt = next((i for i, b in enumerate(gold_seq) if b.upper() != "-"), None)
        last_nt = next((i for i, b in reversed(list(enumerate(gold_seq))) if b.upper() != "-"), None)

        gold_positions = []
        gold_counter = 0
        for i, base in enumerate(gold_seq):
            if first_nt <= i <= last_nt:
                if base.upper() != "-":
                    gold_counter += 1
                    gold_positions.append(str(gold_counter))
                else:
                    gold_positions.append(str(gold_counter))
            else:
                gold_positions.append("NA")

        for record in records:
            sample_id = record.id
            if "Gold_Standard" in sample_id:
                continue

            seq = str(record.seq)
            cod_match = re.search(r"COD[_-]?(\d+)", sample_id, re.IGNORECASE)
            cod = "COD-" + cod_match.group(1) if cod_match else "NA"

            diffs = [
                (i + 1, r.upper(), s.upper())
                for i, (r, s) in enumerate(zip(gold_seq, seq))
                if first_nt <= i <= last_nt and not flu_bases_are_equivalent(r, s)
            ]
            if not diffs:
                continue

            start = diffs[0][0]
            ref_group = [diffs[0][1]]
            sample_group = [diffs[0][2]]

            for j in range(1, len(diffs)):
                pos, r_base, s_base = diffs[j]
                prev_pos, prev_r, prev_s = diffs[j - 1]

                if (
                    pos == prev_pos + 1
                    and base_category(r_base) == base_category(prev_r)
                    and base_category(s_base) == base_category(prev_s)
                ):
                    ref_group.append(r_base)
                    sample_group.append(s_base)
                else:
                    end = prev_pos
                    ref_str = "".join(ref_group)
                    sample_str = "".join(sample_group)
                    result = clasificar(ref_str, sample_str)
                    append_grouped_result(output_rows, cod, eqa, sample_id, gold_positions, start, end, ref_str, sample_str, result)
                    start = pos
                    ref_group = [r_base]
                    sample_group = [s_base]

            end = diffs[-1][0]
            ref_str = "".join(ref_group)
            sample_str = "".join(sample_group)
            result = clasificar(ref_str, sample_str)
            append_grouped_result(output_rows, cod, eqa, sample_id, gold_positions, start, end, ref_str, sample_str, result)

    output_rows.sort(key=sort_key)
    rows_by_lab = defaultdict(list)
    for row in output_rows:
        rows_by_lab[row[0]].append(row)

    write_rows_by_lab(
        rows_by_lab,
        "{cod_lab}_informe_INFLUENZA_2026_fasta_analysis.csv",
        ["COD_LAB", "EQA", "Sample_ID", "POS", "Gold", "Sample", "Resultado"],
    )


def main():
    process_sars()
    process_flu()


if __name__ == "__main__":
    main()
