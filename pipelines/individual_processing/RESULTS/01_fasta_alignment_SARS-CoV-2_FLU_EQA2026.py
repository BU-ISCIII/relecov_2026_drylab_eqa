#!/usr/bin/env python3

import re
import json
import subprocess
from pathlib import Path
from Bio import SeqIO


# ==========================================================
# BASE PATHS
# ==========================================================

BASE_DIR = Path.cwd()
GOLD_DIR = BASE_DIR / "gold_standard"

ALIGN_DIR_SARS = BASE_DIR / "ALIGNMENTS_SARS"
ALIGN_DIR_FLU = BASE_DIR / "ALIGNMENTS_FLU"
FAST_LOG_DIR = BASE_DIR / "fasta_files_submitted_logs"
GLOBAL_ERROR_LOG = FAST_LOG_DIR / "global_fasta_issues.log"


WUHAN_REF = Path("/data/ucct/bi/references/virus/coronaviridae/NC_045512.2.fasta")

NUMBERING_DIR = BASE_DIR / "influenza_reference_numbering"
AH1_FILE = NUMBERING_DIR / "AH1_N1_numbering.fasta"
AH3_FILE = NUMBERING_DIR / "AH3_N2_numbering.fasta"
AH5_FILE = NUMBERING_DIR / "AH5_N1_numbering_2025.fasta"

AH5_samples = {"FLU1"}
AH1_samples = {"FLU2", "FLU7", "FLU9"}
AH3_samples = {"FLU4", "FLU5", "FLU8"}

ALIGN_DIR_SARS.mkdir(exist_ok=True)
ALIGN_DIR_FLU.mkdir(exist_ok=True)
FAST_LOG_DIR.mkdir(exist_ok=True)


# Limpiar log global al inicio
if GLOBAL_ERROR_LOG.exists():
    GLOBAL_ERROR_LOG.unlink()

# ==========================================================
# SEGMENT MAP
# ==========================================================

# Detecta nombres de segmento desde headers
SEGMENT_NAME_MAP = {
    "SEG1": "PB2", "SEGMENT 1": "PB2",
    "SEG2": "PB1", "SEGMENT 2": "PB1",
    "SEG3": "PA",  "SEGMENT 3": "PA",
    "SEG4": "HA",  "SEGMENT 4": "HA",
    "SEG5": "NP",  "SEGMENT 5": "NP",
    "SEG6": "NA",  "SEGMENT 6": "NA",
    "SEG7": "MP",  "SEGMENT 7": "MP",
    "SEG8": "NS",  "SEGMENT 8": "NS",
    "PB2": "PB2",
    "PB1": "PB1",
    "PA": "PA",
    "HA": "HA",
    "NP": "NP",
    "NA": "NA",
    "MP": "MP",
    "NS": "NS",
}

# Conversión nombre → número para numbering
SEGMENT_TO_NUMBER = {
    "PB2": 1,
    "PB1": 2,
    "PA": 3,
    "HA": 4,
    "NP": 5,
    "NA": 6,
    "MP": 7,   # MP corresponde a segment 7
    "NS": 8,
}


# ==========================================================
# UTILIDADES
# ==========================================================

def log_and_print(msg, log):
    print(msg)
    log.write(msg + "\n")


def build_sequence_id(sample_id, submitter_id, original_id):
    def clean(text):
        return re.sub(r"[^\w\-.]", "_", str(text))
    return f"{clean(sample_id)}_{clean(submitter_id)}_{clean(original_id)}"


def find_fasta_recursively(base_dir: Path, filename: str, log):

    matches = list(base_dir.rglob(filename))
    if len(matches) == 1:
        return matches[0]

    if len(matches) > 1:
        log_and_print(f"⚠ Multiple exact matches for {filename}, using first: {matches[0]}", log)
        return matches[0]

    partial = list(base_dir.rglob(f"*{filename}*"))
    if len(partial) == 1:
        log_and_print(f"⚠ Partial match used: {partial[0]}", log)
        return partial[0]

    if len(partial) > 1:
        log_and_print(f"⚠ Multiple partial matches, using first: {partial[0]}", log)
        return partial[0]

    return None


def record_matches_segment(record, segment_number):
    """
    Devuelve True si el record corresponde al número de segmento indicado.
    Detecta:
        - 'segment 4'
        - 'Segment 4'
        - 'Seg 4'
        - 'Seg 4/'
    """

    desc = record.description.lower()

    patterns = [
        rf"segment\s+{segment_number}\b",
        rf"seg\s+{segment_number}\b",
    ]

    for pattern in patterns:
        if re.search(pattern, desc):
            return True

    return False

def log_global_error(cod_name, sample_id, segment, message):
    with open(GLOBAL_ERROR_LOG, "a") as gf:
        line = f"{cod_name}\t{sample_id}\t{segment or '-'}\t{message}\n"
        gf.write(line)

# ==========================================================
# SARS INITIALIZATION
# ==========================================================
# Crea los archivos iniciales donde se van a ir añadiendo las secuencias para alinear
# Teniendo wuham en el top y los 2 gold standards después

def initialize_sars_temporals():

    print("Inicializando temporales SARS...")

    if not WUHAN_REF.exists():
        print("⚠ Wuhan reference not found")
        return

    wuhan_records = list(SeqIO.parse(WUHAN_REF, "fasta"))

    for i in range(1, 11):

        tmp_file = ALIGN_DIR_SARS / f"EQA_SARS{i}_tmp.fasta"

        with open(tmp_file, "w") as out_f:

            for record in wuhan_records:
                record.id = "NC_045512.2_Wuhan"
                record.description = ""
                SeqIO.write(record, out_f, "fasta")

            comp_dir = (
                GOLD_DIR / "consensus" 
            )

            for fasta in comp_dir.glob(f"*SARS{i}_*.fa*"):
                for rec in SeqIO.parse(fasta, "fasta"):
                    SeqIO.write(rec, out_f, "fasta")


# ==========================================================
# FLU INITIALIZATION
# ==========================================================
# Crea los archivos iniciales donde se van a ir añadiendo las secuencias para alinear
# Teniendo las refrencias de gripe que correspondan seguido del gold standards

def initialize_flu_tmp_if_needed(tmp_file, flu_sample, segment):

    if tmp_file.exists():
        return

    print(f"Inicializando {tmp_file.name}")

    if flu_sample in AH5_samples:
        numbering_file = AH5_FILE
    elif flu_sample in AH1_samples:
        numbering_file = AH1_FILE
    elif flu_sample in AH3_samples:
        numbering_file = AH3_FILE
    else:
        numbering_file = None

    with open(tmp_file, "w") as out_f:

        ## Numbering reference
        #if numbering_file and numbering_file.exists():
        #    for rec in SeqIO.parse(numbering_file, "fasta"):
        #        if rec.id.upper().endswith(f"_{segment}"):
        #            SeqIO.write(rec, out_f, "fasta")
        #            break
        # Numbering reference (robust segment detection)
        # Numbering reference
        if numbering_file and numbering_file.exists():
        
            if segment not in SEGMENT_TO_NUMBER:
                raise ValueError(f"Unknown influenza segment name: {segment}")

            segment_number = SEGMENT_TO_NUMBER[segment]

            found_numbering = False

            for rec in SeqIO.parse(numbering_file, "fasta"):
            
                if record_matches_segment(rec, segment_number):
                
                    #rec.id = f"Numbering_reference_{segment}"
                    #rec.description = ""
                    SeqIO.write(rec, out_f, "fasta")
                    found_numbering = True
                    break
                
            if not found_numbering:
                print(f"WARNING: No numbering reference found for segment {segment}")    

        # Gold standard
        comp = "INFL1" if int(flu_sample.replace("FLU", "")) <= 5 else "INFL2"
        comp_dir = GOLD_DIR / "consensus"

        for fasta in comp_dir.glob(f"*{flu_sample}*{segment}*.fa*"):
            for rec in SeqIO.parse(fasta, "fasta"):
                rec.id = f"{flu_sample}_{segment}_Gold_Standard"
                rec.description = ""
                SeqIO.write(rec, out_f, "fasta")


# ==========================================================
# JSON PROCESSING
# ==========================================================

def process_json_and_append(cod_dir: Path):

    cod_dir = Path(cod_dir)
    json_files = list(cod_dir.glob("*.json"))
    if not json_files:
        print(f"No JSON in {cod_dir.name}")
        return

    log_file = FAST_LOG_DIR / f"{cod_dir.name}_processing.log"

    with open(log_file, "w") as log:

        log_and_print(f"\n===== {cod_dir.name} =====", log)

        with open(json_files[0]) as jf:
            data = json.load(jf)

        for entry in data:

            submitter = entry.get("submitting_institution_id")
            organism = entry.get("organism", "")
            sample_id = entry.get("collecting_lab_sample_id")
            consensus_file = entry.get("consensus_sequence_filename")
            consensus_name = entry.get("consensus_sequence_name", "")

            if not consensus_file:
                continue

            fasta_path = find_fasta_recursively(cod_dir / "consensus_files", f"{sample_id}.fa", log)

            if not fasta_path:
                log_and_print(f"❌ FASTA not found: {consensus_file}", log)
                log_global_error(cod_dir.name, sample_id, None, f"FASTA not found: {consensus_file}")
                continue

            log_and_print(f"✓ FASTA found: {fasta_path}", log)

            # ==================================================
            # SARS-CoV-2
            # ==================================================
            if "Severe acute respiratory syndrome coronavirus 2" in organism:

                match = re.search(r"SARS(\d+)", sample_id)
                if not match:
                    continue

                number = int(match.group(1))
                tmp_file = ALIGN_DIR_SARS / f"EQA_SARS{number}_tmp.fasta"

                with open(tmp_file, "a") as out_f:
                    for record in SeqIO.parse(fasta_path, "fasta"):
                        record.id = build_sequence_id(sample_id, submitter, record.id)
                        record.description = ""
                        SeqIO.write(record, out_f, "fasta")

                log_and_print(f"✓ SARS{number} appended to EQA_SARS{number}_tmp.fasta", log)

            # ==================================================
            # Influenza
            # ==================================================
            elif "Influenza virus" in organism:

                match = re.search(r"FLU(\d+)", sample_id)
                if not match:
                    continue

                number = int(match.group(1))
                flu_sample = f"FLU{number}"

                expected_headers = [
                    h.strip() for h in consensus_name.split(",") if h.strip()
                ]

                fasta_records = list(SeqIO.parse(fasta_path, "fasta"))
                record_dict = {rec.id: rec for rec in fasta_records}

                for header in expected_headers:

                    if header not in record_dict:
                        log_and_print(f"❌ Missing header {header}", log)
                        log_global_error(cod_dir.name, sample_id, None, f"Missing header: {header}")
                        continue

                    rec = record_dict[header]
                    header_upper = header.upper()

                    segment_detected = None
                    #for key, val in SEGMENT_MAP.items():
                    for key, val in SEGMENT_NAME_MAP.items():
                        if key in header_upper:
                            segment_detected = val
                            break

                    if not segment_detected:
                        log_and_print(f"❌ Cannot detect segment in {header}", log)
                        log_global_error(cod_dir.name, sample_id, header, "Cannot detect segment")
                        continue

                    tmp_file = ALIGN_DIR_FLU / f"EQA_{flu_sample}_{segment_detected}_tmp.fasta"

                    initialize_flu_tmp_if_needed(tmp_file, flu_sample, segment_detected)

                    rec.id = build_sequence_id(sample_id, submitter, rec.id)
                    rec.description = ""

                    with open(tmp_file, "a") as out_f:
                        SeqIO.write(rec, out_f, "fasta")

                    log_and_print(f"✓ FLU {flu_sample} {segment_detected} appended to EQA_{flu_sample}_{segment_detected}_tmp.fasta", log)


# ==========================================================
# MAFFT
# ==========================================================

def run_mafft_alignments():

    print("\nRunning MAFFT for SARS...")
    for tmp_file in ALIGN_DIR_SARS.glob("*_tmp.fasta"):
        aligned = tmp_file.with_name(tmp_file.name.replace("_tmp", "_aligned"))
        with open(aligned, "w") as out_f:
            subprocess.run(["mafft", "--auto", str(tmp_file)],
                           stdout=out_f,
                           stderr=subprocess.DEVNULL)
        tmp_file.unlink()

    print("\nRunning MAFFT for Influenza...")
    for tmp_file in ALIGN_DIR_FLU.glob("*_tmp.fasta"):
        aligned = tmp_file.with_name(tmp_file.name.replace("_tmp", "_aligned"))
        with open(aligned, "w") as out_f:
            subprocess.run(["mafft", "--auto", str(tmp_file)],
                           stdout=out_f,
                           stderr=subprocess.DEVNULL)
        tmp_file.unlink()


# ==========================================================
# MAIN
# ==========================================================

def main():

    initialize_sars_temporals()

    process_json_and_append("/data/ucct/bi/research/20260217_EQA_SERVIFICATION_TEST/RESULTS/")

    run_mafft_alignments()

    print("\nProceso finalizado correctamente.")


if __name__ == "__main__":
    main()