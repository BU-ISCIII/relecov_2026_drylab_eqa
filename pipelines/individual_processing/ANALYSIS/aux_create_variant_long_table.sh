export BASE_RESULTS_FOLDER=../RAW/
export BASE_RESULTS_FOlDER_RP=$( realpath ../ )
validated_metadata_lab_file=$( find $BASE_RESULTS_FOLDER -name "*.json" )

rm -r variant_long_table
mkdir variant_long_table

echo "CALLING PYTHON SCRIPT TO CREATE VARIANTS LONG TABLE"

python3 make_variants_long_table.py -bd ./queried_vcf -sd ./snpsift -of variants_long_table.csv -vc ivar

echo "FINISHED CREATING MAIN VARIANTS LONG TABLE FILE"

# Correct variant calling (I DON'T CORRECT SAMPLE NAME HERE BECAUSE E.G. FLU10 WILL BE PICKED UP WITH FLU1)
while IFS=$'\t' read -r sample_id consensus_fasta_file vcf_files reference vcaller; do
    sed -i "/${sample_id}_/ s/ivar/${vcaller}/g" variants_long_table.csv
done < samples_files_map.tsv

headers=$( head -n1 variants_long_table.csv )

while IFS=$'\t' read -r sample_id org sequencing; do
    if [[ "$sample_id" == *"FLU"* ]]; then continue; fi
    if ! [ -f "variant_long_table/${org}_${sequencing}_variants_long_table.csv" ]; then echo $headers >> variant_long_table/"${org}_${sequencing}_variants_long_table.csv"; fi
    escaped=$(printf '%s\n' "$sample_id" | sed 's/[][\/.^$*+?(){}|]/\\&/g')
    grep "${sample_id}_" variants_long_table.csv | sed -E "s/${escaped}_[0-9]+/${sample_id}/g" >> "variant_long_table/${org}_${sequencing}_variants_long_table.csv"
done < sample_per_comp.tsv

mv variants_long_table.csv variant_long_table/full_variants_long_table.csv

echo "VARIANTS LONG TABLE GENERATED CORRECTLY"