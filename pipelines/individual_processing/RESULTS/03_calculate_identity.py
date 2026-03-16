from Bio import AlignIO
import sys
import os.path
import json
from pathlib import Path
from statistics import median

BASE_DIR = Path.cwd()

WANTED_COMPARISONS = ["_LC", "_Gold_Standard"]

def find_reference_id(fasta_alignment_file):
    """
    Find the reference sequence based on WANTED_COMPARISONS
    """
    with open(fasta_alignment_file) as f:
        for line in f:
            if line.startswith(">"):
                rid = line.lstrip(">").strip()

                if any(tag in rid for tag in WANTED_COMPARISONS):
                    print(f"Found reference: {rid} in {fasta_alignment_file}")
                    return rid

    raise ValueError(f"No reference found in {fasta_alignment_file}")

def reference_is_fully_masked(ref_seq):
    """
    Returns True if the gold standard contains only Ns or gaps.
    """
    for base in ref_seq.upper():
        if base not in {"N", "-"}:
            return False
    return True

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
    elif isinstance(input_dict[sample_name]["genome_identity_pct"], float):
        input_dict[sample_name]["genome_identity_pct"] = []
    input_dict[sample_name]["genome_identity_pct"].append(identity_value)
    
    with open(output_file, "w") as f:
        json.dump(input_dict, f)
        
def identity_vs_reference(ref_seq, query_seq):
        matches = 0
        compared = 0

        for r, q in zip(ref_seq, query_seq):
            # Ignore positions where reference has gap
            if r == "-" and q == "-":
                continue
            if r == "-" and q == "N":
                continue
            if r == "N" and q == "-":
                continue

            compared += 1
            if r == q:
                matches += 1

        return matches / compared if compared > 0 else 0
        
def write_median_identity(output_file):
    with open(output_file, "r") as f:
        full_file = json.load(f)
    for sample, values in full_file.items():
        all_values_identity = values["genome_identity_pct"]
    
        if not isinstance(all_values_identity, list):
            all_values_identity = [all_values_identity]

        clean_values = [v for v in all_values_identity if isinstance(v, (int, float))]
        if not isinstance(all_values_identity, list):
            all_values_identity = [all_values_identity]

        clean_values = [v for v in all_values_identity if isinstance(v, (int, float))]

        if len(clean_values) == 0:
            full_file[sample]["genome_identity_pct"] = None
        else:
            median_value = median(clean_values)
            full_file[sample]["genome_identity_pct"] = median_value * 100

        if "discrepancy_breakdown" not in full_file[sample]:
            # Zero values for the discrepancy breakdown. assuming if there is an identity AND they're not here
            full_file[sample]["discrepancy_breakdown"] = {
                "wrong_nt": 0,
                "ambiguity2nt": 0,
                "nt2ambigity": 0,
                "ns2nt": 0,
                "nt2ns": 0,
                "insertions": 0,
                "deletions": 0
            }
            print(f"0 discrepancies assumed for sample {sample}")

    with open(output_file, "w") as f:
        json.dump(full_file, f, indent=4)


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

    # Check if gold standard is fully masked
    reference_masked = reference_is_fully_masked(ref_seq)

    if reference_masked:
        print(f"{reference_id} is fully masked (only Ns/gaps). Identity will be NA.")

    for record in alignment:
        if record.id == reference_id:
            continue

        # skip alternative reference-like sequences
        if any(tag in record.id for tag in WANTED_COMPARISONS):
            continue
        
        if "COD" not in record.id and ("SARS" not in record.id or "FLU" not in record.id):
            print(f"Skipping {record.id} as it does not match expected sample naming conventions.")
            continue

        query_seq = str(record.seq)
        record_id = record.id

        if reference_masked:
            identity = None
            print(f"{record.id}\tNA\t(NA)")
        else:
            print(f"Comparing {reference_id} against {record_id}")
            identity = identity_vs_reference(ref_seq, query_seq)
            print(f"{record.id}\t{identity:.4f}\t({identity*100:.2f}%)")

        sample_name = alignment_file_path.name.split("_")[1]
        write_identity(output_file, sample_name, record_id, identity)


if __name__ == '__main__':
    alignment_files = BASE_DIR.rglob("ALIGNMENT*/*.fasta")
    for alignment_file in alignment_files:
        ref_id = find_reference_id(alignment_file)
        main(alignment_file, ref_id, "calculated_values.json")
    write_median_identity("calculated_values.json")