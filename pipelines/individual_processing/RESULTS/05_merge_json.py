import json

def read_json(path):
    with open(path, "r") as f:
        file = json.load(f)
    return file

def write_json(path, dictionary):
    with open(path, "w") as f:
        json.dump(dictionary, f)

def main(consensus_json="calculated_values.json", variants_json="variants_report.json"):
    consensus_json_value = read_json(consensus_json)
    variants_json_value = read_json(variants_json)

    merged_dict = {
        "consensus": consensus_json_value,
        "variants": variants_json_value
    }

    write_json("consolidated_json_reports.json", merged_dict)

if __name__ == "__main__":
    main()