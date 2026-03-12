# relecov_2026_drylab_eqa
Scripts to process the results from the 2026 RELECOV DryLab EQA

# How to run create_merged_json.py.

```
python3 relecov_2026_drylab_eqa/create_merged_json.py \
    --root_folder /path/to/root_folder \
    --expected_data relecov_2026_drylab_eqa/expected_data.json \
    ---output_folder /path/to/root_folder/merged_json_results

# --root_folder: Folder containing all the processed labs: COD-XXXX, COD-YYYY...
#   Expected structure for each COD-XXXX folder:
#       * * RESULTS/
#           - validated_metadata.json
#           - consolidated_json_reports.json
#           - **any other files
```