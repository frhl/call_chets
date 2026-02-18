#!/usr/bin/env bash
# Prepare phenotype file for REGENIE
# REGENIE expects FID and IID columns (family ID and individual ID)

set -euo pipefail

input_dir="input"
pheno_file="${input_dir}/simulated.phenos.with_covariates.tsv"
output_file="${input_dir}/phenotypes.regenie.txt"

echo "Preparing REGENIE phenotype file..."

# Add FID column to match PLINK .fam format (FID=0 for all samples)
# Keep columns 1-17 (IID through PC10), exclude "Unnamed: 11"
awk 'BEGIN {FS=OFS="\t"}
NR==1 {
    # Header: add FID before IID, then columns 1-17
    printf "FID"
    for(i=1; i<=17; i++) printf "\t%s", $i
    printf "\n"
    next
}
{
    # Data rows: FID=0 to match .fam format, then columns 1-17
    printf "0"
    for(i=1; i<=17; i++) printf "\t%s", $i
    printf "\n"
}' "${pheno_file}" > "${output_file}"

echo "Created ${output_file}"
echo "Columns: FID IID phenotype age age2 sex age_sex age2_sex PC1-PC10"
