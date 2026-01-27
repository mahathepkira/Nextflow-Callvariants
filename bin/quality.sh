#!/bin/bash


echo "===== running ====="

if [ -z "$1" ]; then
    echo "Usage: bash quality.sh <input_vcf.gz>"
    exit 1
fi


vcf="$1"
prefix=$(basename "$vcf" .vcf.gz)

vcftools --gzvcf "$vcf" --freq --out "${prefix}"
vcftools --gzvcf "$vcf" --missing-site --out "${prefix}"
vcftools --gzvcf "$vcf" --TsTv-summary --out "${prefix}"

zcat "$vcf" | vcf-annotate --fill-type | grep -oP "TYPE=\\w+" | sort | uniq -c > "${prefix}.summary"

echo "===== finished ====="


