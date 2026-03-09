from Bio import AlignIO
import sys
import os.path
import json
from pathlib import Path
from statistics import median

BASE_DIR = Path.cwd()

WANTED_COMPARISONS = ["SARS1_gold_standard_2026_LC", "_Gold_Standard"]

def find_reference_id(fasta_alignment_file):
    """
    Find last reference in the alinment
    """
    references_ids = []
    with open(fasta_alignment_file, 'r') as f:
        for line in f:
            if line.startswith(">"):
                references_ids.append(line.lstrip(">").rstrip())
    return references_ids[-1]

def write_identity(output_file, sample_name, record_id, identity_value):
    input_dict = {}
    if os.path.exists(output_file):
        with open(output_file, "r") as f:
            input_dict = json.load(f)
    if sample_name not in input_dict:
        input_dict[sample_name] = {
            "genome_identity_pct": []
        }
    if "genome_identity_pct" not in input_dict[sample_name]:
        input_dict[sample_name]["genome_identity_pct"] = []
    input_dict[sample_name]["genome_identity_pct"].append(identity_value)
    
    with open(output_file, "w") as f:
        json.dump(input_dict, f)

def identity_vs_reference(ref_seq, query_seq):
        matches = 0
        compared = 0

        for r, q in zip(ref_seq, query_seq):
            # Ignore positions where reference has gap
            if r == "-":
                continue

            # Ignore positions where query also has gap
            if q == "-":
                continue

            compared += 1
            if r == q:
                matches += 1

        return matches / compared if compared > 0 else 0

def write_median_identity(output_file):
    with open(output_file, "r") as f:
        full_file = json.read(f)
    for sample, values in full_file.items():
        all_values_identity = values["genome_identity_pct"]
        median_value = median(all_values_identity)
        full_file[sample]["genome_identity_pct"] = median_value

    with open(output_file, "w") as f:
        json.dump(full_file, f)


def main(alignment_file_path, reference_id, output_file):

    # Read alignment
    alignment = AlignIO.read(alignment_file_path, "fasta")

    # Find reference sequence
    reference = None
    for record in alignment:
        if record.id == reference_id:
            reference = record
            break

    if reference is None:
        print(f"Reference '{reference_id}' not found in alignment.")
        sys.exit(1)

    ref_seq = str(reference.seq)

    print(f"\nReference: {reference_id}\n")

    for record in alignment:
        if record.id == reference_id:
            continue
        if not any([record.id in wanted for wanted in WANTED_COMPARISONS]):
            continue

        query_seq = str(record.seq)
        record_id = record.id
        print(f"Comparing {reference_id} against {record_id}")
        identity = identity_vs_reference(ref_seq, query_seq)

        print(f"{record.id}\t{identity:.4f}\t({identity*100:.2f}%)")

        sample_name = alignment_file_path.name.split("_")[1]
        write_identity(output_file, sample_name, record_id, identity)
    write_median_identity(output_file)


if __name__ == '__main__':
    alignment_files = BASE_DIR.rglob("ALIGNMENT*/*.fasta")
    for alignment_file in alignment_files:
        ref_id = find_reference_id(alignment_file)
        main(alignment_file, ref_id, "identity_values.json")
