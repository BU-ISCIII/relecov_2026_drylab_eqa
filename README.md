# relecov_2026_drylab_eqa

Scripts to process, summarise, and render the results from the 2026 RELECOV Dry-Lab EQA.

## Main scripts

This repository currently has three top-level Python entry points:

- `create_merged_json.py`: builds one merged JSON per laboratory from the processed lab folders.
- `general_report.py`: aggregates all laboratory JSONs into `general.json` and generates the network and individual figures.
- `render_reports.py`: fills the markdown Jinja template and optionally exports the resulting reports to PDF.

## 1. `create_merged_json.py`

Use this script to generate one consolidated JSON file per laboratory from the processed results folders.

Example:

```bash
python3 create_merged_json.py \
  --root_folder /path/to/root_folder \
  --expected_data expected_data.json \
  --output_folder /path/to/root_folder/merged_json_results
```

Important arguments:

- `--root_folder`: folder containing all processed lab folders such as `COD-XXXX`, `COD-YYYY`, etc.
- `--expected_data`: path to `expected_data.json`.
- `--output_folder`: destination directory for the generated `lab_*.json` files.
- `--lab_cod`: optional, restrict processing to a single laboratory.
- `--lab_name`: optional, override the full laboratory name.

Expected lab folder structure:

```text
COD-XXXX/
  RESULTS/
    validated_metadata.json
    consolidated_json_reports.json
    ...
```

Output:

- one `lab_<LAB>.json` file per processed laboratory

## 2. `general_report.py`

Use this script to generate the network-level summary JSON and all figures used by the report template.

Example:

```bash
python3 general_report.py \
  --expected-data expected_data.json \
  --labs-dir merged_results_per_lab/ \
  --output general.json
```

Important arguments:

- `--expected-data`: path to `expected_data.json`
- `--labs-dir`: directory containing the per-lab merged JSON files
- `--output`: output path for the generated `general.json`
- `--figures-dir`: optional, root folder for generated figures. Default: `./figures`

Output:

- `general.json`
- network figures under `figures/<COMPONENT>/` and `figures/network/`
- individual figures under `figures/labs/<LAB>/<COMPONENT>/`

## 3. `render_reports.py`

Use this script to render the markdown Jinja template into a general report and, optionally, individual laboratory reports. It can also export the rendered reports to PDF using the custom stylesheet in `report_pdf.css`.

Example: general report only

```bash
python3 render_reports.py \
  --template report_template.md \
  --general-json general.json
```

Example: general report plus individual lab reports

```bash
python3 render_reports.py \
  --template report_template.md \
  --general-json general.json \
  --labs-dir merged_json_results/
```

Example: markdown only

```bash
python3 render_reports.py \
  --template report_template.md \
  --general-json general.json \
  --labs-dir merged_json_results/ \
  --markdown-only
```

Example: PDF only

```bash
python3 render_reports.py \
  --template report_template.md \
  --general-json general.json \
  --labs-dir merged_json_results/ \
  --pdf-only
```

Important arguments:

- `--template`: markdown Jinja template to render
- `--general-json`: path to `general.json`
- `--labs-dir`: optional directory containing compatible `lab_*.json` files
- `--output-dir`: optional output directory. Default: `rendered_reports`
- `--css`: optional stylesheet path. Default: `report_pdf.css`
- `--markdown-only`: render only markdown files
- `--pdf-only`: render only PDF files
- `--keep-html`: keep the intermediate HTML files used to create the PDFs

Output:

- markdown files under `rendered_reports/markdown/`
- PDF files under `rendered_reports/pdf/`
- optional intermediate HTML files under `rendered_reports/html/`

## Suggested execution order

If you want to generate the full reporting output from processed lab folders, the usual order is:

1. Run `create_merged_json.py` to build the individual `lab_*.json` files.
2. Run `general_report.py` to generate `general.json` and all figures.
3. Run `render_reports.py` to render the markdown and PDF reports.

## Notes

- `render_reports.py` expects individual report JSONs compatible with the template, that is, files with top-level `lab`, `metadata`, and `components` sections.
- PDF generation currently uses a Chrome/Chromium executable available in `PATH`.
- The report template is `report_template.md` and the PDF stylesheet is `report_pdf.css`.
