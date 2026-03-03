## Individual EQA preprocessing/processing

### Steps

1. Download this folder where you have the data - The VCF, fasta and validated JSON files. It doesn't matter the actual organisation of the folder - The scripts will search recursively and create symlinks.
2. Run `bash lablog`. This will:
    - Move the contents of the folder to `RAW`
    - Create the base SNPEff database + empty config file
3. Move to ANALYSIS: `cd ANALYSIS`
4. Run the lablog (`bash lablog`). This will read the validated JSON and:
    - Create the necessary scripts (`_[0-9]...`) to run the preprocessing
    - Generate the supporting metadata tsv files
5. Run each step sequentially:
    - Activate micromamba environment with pandas installed and load singularity: `micromamba activate biopython_v1.84 && module load singularity`
    - Run `bash _00_create_vcf_files.sh`: This will run a python script that will take the VCF files in `RAW`, append columns with a specific name (DP, REF_DP and ALT_DP), and ensure format of the fields is actually correct (Integers where specified, no "0.0"). The resulting VCF files are dumped in 00-vcf.
    - Run `bash _01_bcftools_sort.sh`: Uses BCFTools to sort the VCF files and tabix to generate an index.
    - Run `bash _02_snpeff_build.sh`: Only needs to be run once; it generates a SNPEff database based off the references indicated in the metadata. Downloads references in GBK format and creates the database with `snpEff build`
    - Run `bash _03_snpeff_annotate.sh`: Using the downloaded databases, annotate the VCF files. **PLEASE NOTE**: vcf to reference in the metadata happens in array order; if this is not correct, modify the file `vcf_to_ref.tsv`, run `_99_rewrite_vcf_refs.sh` and then run step _03
    - Run `bash _04_bcftools_query_ivar_header.sh`: Uses BCFtools query tool to extract specific headers. Since we added extra columns in step 0, this should always work
    - Run `bash _05_running_snpsif.sh`: This will extract the fields with SnpSift
    - Run `bash _06_variant_long_table.sh`: This will launch the python script `make_variants_long_table.py` to create the variants_long_table, but with a catch: After the python script execution, the file is divided based off organism + technology (e.g. SARS + Illumina, or SARS + Nanopore)
6. This should have generated the necessary files under `variant_long_table` (VCF file output) and `00-consensus` (Fasta file "output")
7. Go to `RESULTS`: `cd ../RESULTS` and run the lablog: `bash lablog`
8. Activate the MAFFT module: `module load MAFFT/7.475-gompi-2020b-with-extensions` (And the previous modules/micromamba environments)
9. Run the steps:
    - Run `bash _01_generate_fasta_analysis_reports.sh`: This generates the analysis reports based off the consensus fasta files
    - Run `bash _02_vcf_comparison.sh`: This generates the analysis reports (CSV files) based off the VCF results, compared to the gold standard

From here, the analysis needs to be performed - (Probably in `global_processing`, running the plot/report tools for all the results combined?)

### Possible choke points:
- If `_00` fails immediately, you may need to activate an environment with pandas 
- If `_01` fails, it may be an issue with the column value re-writing happening in step _00.
- if `_02` fails, sometimes it's a race condition on trying to write to the DOC/snpeff.config file or an issue with NCBI; re-run. Don't worry, the databases already downloaded will not be overwritten.
- If `_04` fails, it probably is an issue with the column appending happening in step 01

From this point forward, if everything went correctly it's very difficult for anything to fail


### Notes

Please note:
- These scripts are specific to a SLURM cluster. They could be re-used, but the `lablog` files need to be modified to use the infrastructure that you need.
- The gold standard files could not be published, so the lablog within analysis includes full path within the filesystem used for the analysis, not actual re-usable paths.
- The `_02` script only needs to be run once per analysis - Although repeated runs will not hurt since it will detect that the database is already created
- All the scripts except _05 generates log files under the folder `logs` - _05 needs the stdout output to generate files, so a log could not be created.