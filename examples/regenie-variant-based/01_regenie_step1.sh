#!/usr/bin/env bash
# REGENIE Step 1: Fit the null model using whole-genome regression
# This step fits a polygenic model using common variants to get predictions
# for use in Step 2 variant testing

set -euo pipefail

# Input/Output paths
input_dir="input"
output_dir="output"
mkdir -p "${output_dir}"

# Input files
bed_prefix="${input_dir}/simulated"
pheno_file="${input_dir}/phenotypes.regenie.txt"
sample_list="${input_dir}/sample_list.txt"

# Output prefix
out_prefix="${output_dir}/regenie_step1"

# Phenotype and covariates
pheno_col="phenotype"
covariates="age,age2,sex,age_sex,age2_sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"
cat_covariates="sex"

# REGENIE parameters
n_threads=4
bsize=100  # Block size for ridge regression (smaller to avoid low-variance SNPs)
maf=0.01    # Minimum MAF for Step 1

# Docker image
docker_image="ghcr.io/rgcgithub/regenie/regenie:v4.1.gz"

echo "======================================"
echo "REGENIE Step 1: Null Model Fitting"
echo "======================================"
echo "Input genotypes: ${bed_prefix}.bed/.bim/.fam"
echo "Phenotype file: ${pheno_file}"
echo "Output prefix: ${out_prefix}"
echo "Phenotype: ${pheno_col}"
echo "Covariates: ${covariates}"
echo ""

# Pull docker image if not present
echo "Checking for Docker image..."
if ! docker image inspect "${docker_image}" >/dev/null 2>&1; then
    echo "Pulling Docker image ${docker_image}..."
    docker pull "${docker_image}"
fi

# Create a keep file (FID IID format) from .fam file
echo "Creating keep file..."
# .fam file has FID=0, IID=sample_X, so we need to match that format
awk '{print 0, $1}' "${sample_list}" > "${output_dir}/keep_file.txt"

# Copy and filter genotype files using PLINK2 docker to remove low-MAF variants
echo "Filtering genotypes (MAF > ${maf}) using PLINK2..."

# First copy files (resolving symlinks) to output directory
cp -L "${bed_prefix}.bed" "${output_dir}/simulated_unfiltered.bed"
cp -L "${bed_prefix}.bim" "${output_dir}/simulated_unfiltered.bim"
cp -L "${bed_prefix}.fam" "${output_dir}/simulated_unfiltered.fam"

# Then filter using PLINK2
docker run --rm \
  -v "$(pwd)/${output_dir}":/data \
  biocontainer/plink2:alpha2.3_jan2020 \
  plink2 --bfile /data/simulated_unfiltered \
    --maf ${maf} \
    --make-bed \
    --out /data/simulated \
    --allow-extra-chr

# Clean up unfiltered files
rm -f "${output_dir}"/simulated_unfiltered.*

# Run REGENIE Step 1
echo ""
echo "Running REGENIE Step 1..."
docker run --rm \
  -v "$(pwd)/${input_dir}":/data/input:ro \
  -v "$(pwd)/${output_dir}":/data/output \
  "${docker_image}" \
  regenie \
    --step 1 \
    --bed /data/output/simulated \
    --phenoFile /data/input/phenotypes.regenie.txt \
    --phenoCol "${pheno_col}" \
    --covarFile /data/input/phenotypes.regenie.txt \
    --covarColList "${covariates}" \
    --catCovarList "${cat_covariates}" \
    --keep /data/output/keep_file.txt \
    --bsize ${bsize} \
    --threads ${n_threads} \
    --qt \
    --apply-rint \
    --out /data/output/regenie_step1

echo ""
echo "======================================"
echo "REGENIE Step 1 Complete!"
echo "======================================"
echo "Output files:"
echo "  - ${out_prefix}.log"
echo "  - ${out_prefix}.loco (leave-one-chromosome-out predictions)"
echo "  - ${out_prefix}_pred.list (list of prediction files)"
echo ""
echo "Next step: Run 02_regenie_step2.sh for variant testing"
