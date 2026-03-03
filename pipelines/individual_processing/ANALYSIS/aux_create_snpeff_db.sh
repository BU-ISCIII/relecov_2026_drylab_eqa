#!/bin/bash

ref_id=$1

# Add genome file to config
if grep -q $ref_id ../DOC/snpeff.config;then 
    echo "Reference is already downloaded"
    exit 0
fi

mkdir ../DOC/snpeff_db/$ref_id

URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${ref_id}&rettype=gbwithparts&retmode=text"
OUTPUT="../DOC/snpeff_db/${ref_id}/genes.gbk"
MAX_RETRIES=20
DEFAULT_WAIT=15

for ((i=1; i<=MAX_RETRIES; i++)); do
    echo "Attempt $i..."

    # Capture headers and status
    response=$( wget -S -O "${OUTPUT}.tmp" "${URL}" 2>&1 )
    status=$( echo "$response" |  awk '/^  HTTP/{print $2}' )

    if [ "$status" = "200" ]; then
        echo "Download successful"
        mv "$OUTPUT.tmp" "$OUTPUT"
        break
    fi

    if [ $status -eq 429 ]; then

        echo "Rate limited. Waiting $retry_after seconds..."
        sleep $DEFAULT_WAIT
        continue
    else
        echo "Failed with HTTP $status"
        exit 1
    fi
    echo "Max retries reached"
    exit 1
done



echo $ref_id".genome: " $ref_id >> ../DOC/snpeff.config

snpEff build -config ../DOC/snpeff.config -dataDir ../DOC/snpeff_db -v $ref_id