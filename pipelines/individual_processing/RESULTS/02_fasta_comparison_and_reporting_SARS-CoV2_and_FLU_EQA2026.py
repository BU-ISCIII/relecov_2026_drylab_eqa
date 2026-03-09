from Bio import SeqIO
import csv
import re
from pathlib import Path
from collections import defaultdict
import pandas as pd
import json


# ==========================================================
# PATHS
# ==========================================================

BASE_DIR = Path.cwd()

SARS_DIR = BASE_DIR / "ALIGNMENTS_SARS"
FLU_DIR = BASE_DIR / "ALIGNMENTS_FLU"  
FASTA_ANALYSIS = BASE_DIR / "FASTA_ANALYSIS_REPORTS"
INDIVIDUAL_REPORTS = FASTA_ANALYSIS 

FASTA_ANALYSIS.mkdir(exist_ok=True)
INDIVIDUAL_REPORTS.mkdir(exist_ok=True)

# ==========================================================
# FUNCIONES
# ==========================================================

def calculate_and_write_discrepancies(df: pd.DataFrame, output_path: str):
    discrepancy_breakdown = {}
    table_map = {
              "wrong_nt": "Wrong nucleotide",
              "ambiguity2nt": "Nucleotide instread of ambiguity",
              "nt2ambigity": "Ambiguity instead of nucleotide",
              "ns2nt": "Nucleotide stretch instead of stretch of Ns",
              "nt2ns": "Stretch of Ns instead of nucleotide" ,
              "insertions": "Insertion relative to gold standard",
              "deletions": "Deletion relative to gold standard"
            }
    sample_names = list(df["Sample_ID"])
    for sample in sample_names:
        if sample not in discrepancy_breakdown:
            discrepancy_breakdown[sample] = {"discrepancy_breakdown": {}}
        df_by_sample = df.groupby(sample)
        for json_key, df_key in table_map.items():
            discrepancy_breakdown[sample]["discrepancy_breakdown"][json_key] = df_by_sample.count(df_by_sample["Resultado"] == df_key) or 0
    
    if not Path(output_path).is_file:
        dict_discrepancies = {}
    else:
        with open(output_path, 'r') as f:
            dict_discrepancies = json.load(f)
    
    dict_discrepancies.update(discrepancy_breakdown)
    
    with open(output_path, 'w') as f:
        json.dump(dict_discrepancies, f)

def base_category(base):
    base = base.upper()
    if base == "N":
        return "N"
    elif base == "-":
        return "gap"
    else:
        return "nt"


def clasificar(ref_bases, sample_bases):
    r = ref_bases.upper()
    s = sample_bases.upper()

    # Si son iguales → no hay cambio
    if r == s:
        return ""

    if len(r) == 1:
        if s == "N":
            return "Stretch of Ns instead of nucleotide"
        elif s == "-":
            return "Deletion relative to gold standard"
        elif r == "-":
            return "Insertion relative to gold standard"
        elif r == "N":
            return "Nucleotide stretch instead of stretch of Ns"
        elif s in "ATGC" and r in "ATGC":
            return "Wrong nucleotide"
        elif s in "RYSWKMBDHV" and r in "ATGC":
            return "Ambiguity instead of nucleotide"
        elif r in "RYSWKMBDHV" and s in "ATGC":
            return "Nucleotide instead of ambiguity"
        else:
            return "Other: Review"

    else:
        if set(r) <= {"-"} and set(s) <= {"N"}:
            return "Nucleotide instead of gap"
        elif set(r) <= {"N"} and set(s) <= {"-"}:
            return "Deletion relative to gold standard"
        elif "-" in s and any(base in "ATGCRYSWKMBDHV" for base in r) and "N" not in r:
            return "Deletion relative to gold standard"
        elif ("-" not in s and "N" not in s) and "-" in r:
            return "Insertion relative to gold standard"
        elif "N" in s and not any(x in r for x in ["N", "-"]):
            return "Stretch of Ns instead of nucleotide"
        elif "N" in r and not any(x in s for x in ["N", "-"]):
            return "Nucleotide stretch instead of stretch of Ns"
        elif "N" in r and "-" in s:
            return "Deletion relative to gold standard"
        elif all(base in "ATGCRYSWKMBDHV" for base in r) and all(base in "ATGCRYSWKMBDHV" for base in s):
            return "Wrong nucleotide"
        else:
            return "Other: Review"


def sort_key(row):
    cod_num = int(re.search(r"(\d+)", row[0]).group(1)) if re.search(r"(\d+)", row[0]) else 99999
    eqa_num = re.search(r"(\d+)", row[1])
    eqa_key = int(eqa_num.group(1)) if eqa_num else 99999
    return (cod_num, eqa_key)

# ==========================================================
# PROCESAMIENTO SARS
# ==========================================================

output_rows = []

for fasta_file in SARS_DIR.glob("*_aligned.fasta"):

    filename = fasta_file.name
    eqa_match = re.search(r"(EQA[_-]SARS\d+)", filename, re.IGNORECASE)
    eqa = eqa_match.group(1) if eqa_match else "NA"

    records = list(SeqIO.parse(fasta_file, "fasta"))

    # estructura mínima esperada, al menos ref numbering +  2 gold + sample
    if len(records) < 4:
        continue

    # ======================================================
    #  NUMERACIÓN SEGÚN WUHAN
    # ======================================================

    ref_seq = str(records[0].seq)

    # Encontrar índices del primer y último nucleótido real
    first_nt = next((i for i, b in enumerate(ref_seq) if b.upper() in "ATGCRYKMSWBDHV"), None)
    last_nt  = next((i for i, b in reversed(list(enumerate(ref_seq))) if b.upper() in "ATGCRYKMSWBDHV"), None)

    ref_seq_positions = []
    ref_counter = 0
    for i, base in enumerate(ref_seq):
        if first_nt <= i <= last_nt:
            if base.upper() in "ATGCRYKMSWBDHV":
                ref_counter += 1
                ref_seq_positions.append(str(ref_counter)) # gaps dentro de la región real
            else:
                ref_seq_positions.append(str(ref_counter))
        else:
            ref_seq_positions.append("NA")  # antes del primer y después del último nucleótido

    
    # ======================================================
    # 2GOLD LOW (LC)
    # ======================================================
    # SARS-CoV-2, Comparación SOLO contra Gold Standard LOW (LC), 
    # el uso de HC resulta muy muy sucio para los informes      
    gold_seq = str(records[2].seq)


    # ======================================================
    # MUESTRAS
    # ======================================================

    for record in records[3:]:
        
        sample_id = record.id
        seq = str(record.seq)

        cod_match = re.search(r"COD[_-]?(\d+)", sample_id, re.IGNORECASE)
        cod = "COD-" + cod_match.group(1) if cod_match else "NA"

        # Solo comparar dentro de la región numerada de la referencia
        diffs = [
            (i+1, r.upper(), s.upper()) # i+1 porque enumerate empieza en 0
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
            prev_pos, prev_r, prev_s = diffs[j-1]

            if (pos == prev_pos + 1 and
                base_category(r_base) == base_category(prev_r) and
                base_category(s_base) == base_category(prev_s)):

                ref_group.append(r_base)
                sample_group.append(s_base)
            else:
                end = prev_pos
                ref_str = "".join(ref_group)
                sample_str = "".join(sample_group)
                result = clasificar(ref_str, sample_str)
                pos_str = ref_seq_positions[start-1] if start == end else f"{ref_seq_positions[start-1]}-{ref_seq_positions[end-1]}"
                output_rows.append([cod, eqa, sample_id, pos_str, ref_str, sample_str, result])
                start = pos
                ref_group = [r_base]
                sample_group = [s_base]

        # último bloque
        end = diffs[-1][0]
        ref_str = "".join(ref_group)
        sample_str = "".join(sample_group)
        result = clasificar(ref_str, sample_str)
        pos_str = ref_seq_positions[start-1] if start == end else f"{ref_seq_positions[start-1]}-{ref_seq_positions[end-1]}"
        output_rows.append([cod, eqa, sample_id, pos_str, ref_str, sample_str, result])

# ==========================================================
# ORDENAR GLOBAL SARS
# ==========================================================
output_rows.sort(key=sort_key)

# ==========================================================
# ESCRIBIR CSV SARS
# ==========================================================

# Guardar informes individuales SARS
rows_by_lab = defaultdict(list)
for row in output_rows:
    rows_by_lab[row[0]].append(row)
for cod_lab, rows in rows_by_lab.items():
    df = pd.DataFrame(rows, columns=[
            "COD_LAB",
            "EQA",
            "Sample_ID",
            "POS (NC_045512.2)",
            "Gold_LOW",
            "Sample",
            "Resultado"
            ])
    calculate_and_write_discrepancies(df, "calculated_values.json")
    filename = INDIVIDUAL_REPORTS / f"{cod_lab}_informe_SARS-CoV-2_2026_fasta_analysis.csv"
    with open(filename, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "COD_LAB",
            "EQA",
            "Sample_ID",
            "POS (NC_045512.2)",
            "Gold_LOW",
            "Sample",
            "Resultado"
            ])
        writer.writerows(rows)

# ==========================================================
# PROCESAMIENTO FLU
# ==========================================================

output_rows_flu = []

for fasta_file in FLU_DIR.glob("*_aligned.fasta"):

    filename = fasta_file.name
    eqa_match = re.search(r"(EQA[_-]FLU\d+)", filename, re.IGNORECASE)
    eqa = eqa_match.group(1) if eqa_match else "NA"

    records = list(SeqIO.parse(fasta_file, "fasta"))
    if len(records) < 3:
        continue # referencia + gold + muestra


    # ======================================================
    # NUMERACIÓN SEGÚN REFERENCIA TOP
    # ======================================================
    ref_seq = str(records[0].seq)
    first_nt = next((i for i, b in enumerate(ref_seq) if b.upper() in "ATGCRYKMSWBDHV"), None)
    last_nt  = next((i for i, b in reversed(list(enumerate(ref_seq))) if b.upper() in "ATGCRYKMSWBDHV"), None)

    ref_positions = []
    ref_counter = 0
    for i, base in enumerate(ref_seq):
        if first_nt <= i <= last_nt:
            if base.upper() in "ATGCRYKMSWBDHV":
                ref_counter += 1
                ref_positions.append(str(ref_counter))
            else:
                ref_positions.append(str(ref_counter))
        else:
            ref_positions.append("NA")

    # ======================================================
    # GOLD STANDARD (único)
    # ======================================================
    gold_seq = str(records[1].seq)

    # ======================================================
    # MUESTRAS
    # ======================================================

    for record in records[2:]:
        sample_id = record.id
        seq = str(record.seq)
        cod_match = re.search(r"COD[_-]?(\d+)", sample_id, re.IGNORECASE)
        cod = "COD-" + cod_match.group(1) if cod_match else "NA"

        diffs = [
            (i+1, r.upper(), s.upper())
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
            prev_pos, prev_r, prev_s = diffs[j-1]

            if (pos == prev_pos + 1 and
                base_category(r_base) == base_category(prev_r) and
                base_category(s_base) == base_category(prev_s)):

                ref_group.append(r_base)
                sample_group.append(s_base)
            else:
                end = prev_pos
                ref_str = "".join(ref_group)
                sample_str = "".join(sample_group)
                result = clasificar(ref_str, sample_str)
                pos_str = ref_positions[start-1] if start == end else f"{ref_positions[start-1]}-{ref_positions[end-1]}"
                output_rows_flu.append([cod, eqa, sample_id, pos_str, ref_str, sample_str, result])
                start = pos
                ref_group = [r_base]
                sample_group = [s_base]

        # último bloque
        end = diffs[-1][0]
        ref_str = "".join(ref_group)
        sample_str = "".join(sample_group)
        result = clasificar(ref_str, sample_str)
        pos_str = ref_positions[start-1] if start == end else f"{ref_positions[start-1]}-{ref_positions[end-1]}"
        output_rows_flu.append([cod, eqa, sample_id, pos_str, ref_str, sample_str, result])

# ==========================================================
# ORDENAR FLU
# ==========================================================
output_rows_flu.sort(key=sort_key)

# ==========================================================
# ESCRIBIR CSV FLU
# ==========================================================

# Guardar informes individuales FLU
rows_by_lab_flu = defaultdict(list)
for row in output_rows_flu:
    rows_by_lab_flu[row[0]].append(row)
for cod_lab, rows in rows_by_lab_flu.items():
    df = pd.DataFrame(rows, columns=[
            "COD_LAB",
            "EQA",
            "Sample_ID",
            "POS (NC_045512.2)",
            "Gold_LOW",
            "Sample",
            "Resultado"
            ])
    calculate_and_write_discrepancies(df, "calculated_values.json")
    filename = INDIVIDUAL_REPORTS / f"{cod_lab}_informe_INFLUENZA_2026_fasta_analysis.csv"
    with open(filename, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "COD_LAB",
            "EQA",
            "Sample_ID",
            "POS",
            "Gold",
            "Sample",
            "Resultado"
            ])
        writer.writerows(rows)
