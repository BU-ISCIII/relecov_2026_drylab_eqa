import pandas as pd
import re
import sys

from collections import OrderedDict

## THIS FILE COULD BE REPURPOSED OUTSIDE OF EQA


class VCFFile:
    TYPE_COMPARISONS={
        "String": str,
        "Integer": lambda x: int(float(x)),
        "Float": float
    }
    
    def __init__(self, path:str):
        self.read_vcf_file(path)
        self.replace_NA_values("0")
        if self.headers['FORMAT']:
            self.validate_format_data()
        
        # Non-necessary, just convenient, values
        self.format_keys = [key["ID"] for key in self.headers["FORMAT"]]
        self.info_keys = [key["ID"] for key in self.headers["INFO"]]

    def read_vcf_file(self, path):
        self.headers = OrderedDict()   
        tabular_data = []
        with open(path, "r") as f:
            full_file = f.read().splitlines()
        for line in full_file:
            if line.startswith("##"):
                line = line.replace("##", "")
                line = line.split("=")
                header_name = line[0]
                if header_name not in self.headers:
                    self.headers[header_name] = []
                rest_of_values = line[1] if not line[1].startswith("<") else self.read_dict_from_line(line[1:])
                self.headers[header_name].append(rest_of_values)
            else:
                tabular_data.append(line.split("\t"))
        self.headers["FORMAT"] = self.headers.get("FORMAT", [])
        self.data = pd.DataFrame(tabular_data[1:])
        self.data.columns = tabular_data[0]
        self.sample_values_column = tabular_data[0][-1] #FIXME I DON'T THINK THIS IS CORRECT

    @staticmethod
    def read_dict_from_line(line: list) -> dict:
        """
        Method to read a dictionary from the header line. Assumes comma-separated keys and key=value format.

        Args:
            line (list): _description_
        """
        line = "=".join(line)
        line = re.split(',(?=(?:[^"]*"[^"]*")*[^"]*$)', line)
        return {key.lstrip("<"): value.rstrip(">") for key, value in [element.split("=", maxsplit=1) for element in line]}

    def validate_format_data(self):
        v_format = self.headers["FORMAT"]
        for index, line in self.data.iterrows():
            order = line["FORMAT"].split(":")
            values_sample = line[self.sample_values_column].split(":")
            for i, element in enumerate(order):
                eval_format = next(filter(lambda x: x["ID"] == element, v_format))
                values_comma_separated = values_sample[i].split(",")
                values_comma_separated = [self.TYPE_COMPARISONS[eval_format["Type"]](value) for value in values_comma_separated]
                values_sample[i] = ",".join([str(value) for value in values_comma_separated])
            self.data.at[index, self.sample_values_column] = ":".join(
                str(value) for value in values_sample
            )
        print("Done validating")
            
    def write(self, output_path):
        with open(output_path, "w") as f:
            for key, value in self.headers.items():
                for header_element in value:
                    if isinstance(header_element, dict):
                        f.write(f"##{key}=<{','.join(f'{k}={v}' for k, v in header_element.items())}>\n")
                    else:
                        f.write(f"##{key}={header_element}\n")
            f.write(self.data.to_csv(sep="\t", index=False))

    def add_format_value(self, format_dict: dict, values=None, default=""):
        assert values or default, f"Could not add format {format_dict.get('ID', '')} due to lacking values list or default"
        assert not isinstance(default, list), "Please use valid default values; can only provide single numbers or str"
        new_format = {}
        for mandatory_key in ["ID", "Type", "Number", "Description"]:
            assert format_dict.get(mandatory_key, None), f"Could not find value for mandatory key {mandatory_key}"
            new_format[mandatory_key] = format_dict[mandatory_key]
        assert "'" in new_format["Description"] or '"' in new_format["Description"], "Description field must be encased in quotes"
        assert all([format_dict["ID"] != format_v["ID"] for format_v in self.headers["FORMAT"]]), "Format header is already present"
        self.headers["FORMAT"].append(new_format)
        if "FORMAT" not in self.data:
            self.data["FORMAT"] = pd.Series([format_dict["ID"]] * len(self.data))
            self.sample_values_column = "SAMPLE"
            self.data[self.sample_values_column] = pd.Series(values)
        else:
            self.data.FORMAT = self.data.FORMAT.apply(lambda x: f"{x}:{format_dict['ID']}")
            self.data[self.sample_values_column] = self.data[self.sample_values_column] + ":" + pd.Series(values)
        self.validate_format_data()

    # This function could call to a config, but https://www.youtube.com/watch?v=waEC-8GFTP4&list=RDwaEC-8GFTP4&start_radio=1&pp=ygUeYWluJ3Qgbm9ib2R5IGdvdCB0aW1lIGZvciB0aGF0oAcB
    # Also updating this function should be pretty easy - new required header? add case. New VCF format? add "ifs" within header keys
    def find_equivalent(self, header_name):
        values = []
        format_dict = {}
        if header_name in self.format_keys:
            return {}, []
        match header_name:
            case "DP":
                if "DP" in self.info_keys:
                    for index, row in self.data.iterrows():
                        row_info = row['INFO'].split(";")
                        dp = next(filter(lambda x: x.startswith("DP="), row_info))
                        values.append(dp.split("=")[1])
                    format_dict = {
                        "ID": "DP",
                        "Number": "1",
                        "Type": "Integer",
                        "Description": '"Total Depth"'
                    }
                elif "DP" not in self.format_keys:
                    format_dict = {
                        "ID": "DP",
                        "Number": "1",
                        "Type": "Integer",
                        "Description": '"Total Depth"'
                    }
                    print("WARNING: DP value was not found - Setting it to '0'")
                    values.append(["0"] * len(self.data))

            case "REF_DP":
                format_dict = {
                        "ID": "REF_DP",
                        "Number": "1",
                        "Type": "Integer",
                        "Description": '"Depth of reference base"'
                    }
                if "AD" in self.format_keys:
                    # Nanopore, Allelic Depths (Ref base)
                    for index, row in self.data.iterrows():
                        index_AD = row["FORMAT"].split(":").index("AD")
                        value = row[self.sample_values_column].split(":")[index_AD].split(",")[0]
                        if value == ".":
                            # e.g. Octopus defaults missing values in haploid samples to missing, or '.'
                            value == "0"
                        values.append(value)
                elif "DP" in self.info_keys and "ALT_DP" in self.format_keys:
                    # ORDER MATTERS - ASSUMING DP is already inserted
                    for index, row in self.data.iterrows():
                        index_alt_dp = row["FORMAT"].split(":").index("ALT_DP")
                        index_dp = row["FORMAT"].split(":").index("DP")
                        value = str(int(row[self.sample_values_column].split(":")[index_dp].split(",")[0]) - int(row[self.sample_values_column].split(":")[index_alt_dp].split(",")[0]))
                        values.append(value)
                elif "DP4" in self.info_keys:
                    for index, row in self.data.iterrows():
                        row_info = row['INFO'].split(";")
                        dp4 = next(filter(lambda x: x.startswith("DP4="), row_info))
                        dp4_values = [int(value) for value in dp4.split("=")[1].split(",")[0:2]]
                        values.append(str(sum(dp4_values)))
                elif "RO" in self.info_keys:
                    for index, row in self.data.iterrows():
                        row_info = row['INFO'].split(";")
                        dep = next(filter(lambda x: x.startswith("RO="), row_info))
                        dep_values = dep.split("=")[1]
                        values.append(dep_values)
                else:
                    print("WARNING: REF_DP value was not found - Setting it to '0'")
                    values.append(["0"] * len(self.data))


            case "ALT_DP":
                if "AD" in self.format_keys:
                    # Nanopore, Allelic Depths (Alternate base)
                    for index, row in self.data.iterrows():
                        index_AD = row["FORMAT"].split(":").index("AD")
                        value = row[self.sample_values_column].split(":")[index_AD].split(",")[1]
                        if value == ".":
                            # e.g. Octopus defaults missing values in haploid samples to missing, or '.'
                            value == "0"
                        values.append(value)
                elif "DP4" in self.info_keys:
                    for index, row in self.data.iterrows():
                        row_info = row['INFO'].split(";")
                        dp4 = next(filter(lambda x: x.startswith("DP4="), row_info))
                        dp4_values = [int(value) for value in dp4.split("=")[1].split(",")[2:4]]
                        values.append(str(sum(dp4_values)))
                elif "AO" in self.info_keys:
                    for index, row in self.data.iterrows():
                        row_info = row['INFO'].split(";")
                        dep = next(filter(lambda x: x.startswith("AO="), row_info))
                        dep_values = dep.split("=")[1]
                        values.append(dep_values)
                else:
                    print("WARNING: ALT_DP value was not found - Setting it to '0'")
                    values.append(["0"] * len(self.data))
                format_dict = {
                        "ID": "ALT_DP",
                        "Number": "1",
                        "Type": "Integer",
                        "Description": '"Depth of alternate base"'
                    }
                
        return format_dict, values
    
    def replace_NA_values(self, replacement):
        self.data[self.sample_values_column] = self.data[self.sample_values_column].apply(lambda x: x.replace("NA", replacement))
        self.data[self.sample_values_column] = self.data[self.sample_values_column].apply(lambda x: x.replace(".", "0"))
        self.data["INFO"] = self.data["INFO"].apply(lambda x: x.replace("=NA", f"={replacement}"))

    def rename_chrom(self, replacement):
        if "contig" in self.headers:
            self.headers['contig'][0]['ID'] = replacement
        else:
            self.headers['contig'] = [{"ID": replacement}]
        replacement = pd.Series([replacement] * len(self.data["#CHROM"]))
        self.data["#CHROM"] = replacement

if __name__ == '__main__':
    vcf = VCFFile(sys.argv[1])
    reference = sys.argv[3]
    format_headers_needed = ["DP","REF_DP", "ALT_DP"]
    for header in format_headers_needed:
        format_dict, format_values = vcf.find_equivalent(header)
        try:
            vcf.add_format_value(format_dict, format_values)
        except AssertionError as e:
            print(e)

    if reference:
        vcf.rename_chrom(reference)

    vcf.write(sys.argv[2])