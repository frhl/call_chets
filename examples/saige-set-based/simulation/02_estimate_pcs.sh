#!/usr/bin/env bash
# Estimate PCs from genotype data using PLINK2 (via Docker)
#
# This script:
# 1. Converts VCF to PLINK format (if needed)
# 2. Prunes variants for LD
# 3. Estimates principal components

set -euo pipefail

# Configuration
PLINK2_IMAGE="biocontainer/plink2:alpha2.3_jan2020"
INPUT_DIR="$(cd "$(dirname "$0")/.." && pwd)/input"
VCF_FILE="${INPUT_DIR}/simulated.vcf.gz"
PLINK_PREFIX="${INPUT_DIR}/simulated"
OUTPUT_DIR="${INPUT_DIR}/pca"
N_PCS=10  # Number of PCs to compute

# Create output directory
mkdir -p ${OUTPUT_DIR}

echo "=========================================="
echo "Estimating Principal Components"
echo "=========================================="
echo "Using Docker image: ${PLINK2_IMAGE}"
echo ""

# Step 1: Convert VCF to PLINK (if not already done)
if [[ ! -f "${PLINK_PREFIX}.bed" ]]; then
    echo "Step 1: Converting VCF to PLINK format..."
    docker run --rm \
        -v ${INPUT_DIR}:/data \
        -w /data \
        ${PLINK2_IMAGE} \
        plink2 --vcf simulated.vcf.gz \
               --make-bed \
               --out simulated \
               --allow-extra-chr
    echo "  Done!"
    echo ""
else
    echo "Step 1: PLINK files already exist, skipping conversion"
    echo ""
fi

# Step 2: Prune variants for LD
echo "Step 2: LD pruning..."
docker run --rm \
    -v ${INPUT_DIR}:/data \
    -w /data \
    ${PLINK2_IMAGE} \
    plink2 --bfile simulated \
           --indep-pairwise 50 5 0.2 \
           --out pca/ld_pruned \
           --allow-extra-chr

echo "  Variants passing LD pruning: $(wc -l < ${OUTPUT_DIR}/ld_pruned.prune.in)"
echo ""

# Step 3: Calculate PCs
echo "Step 3: Calculating ${N_PCS} principal components..."
docker run --rm \
    -v ${INPUT_DIR}:/data \
    -w /data \
    ${PLINK2_IMAGE} \
    plink2 --bfile simulated \
           --extract pca/ld_pruned.prune.in \
           --pca ${N_PCS} \
           --out pca/pca_results \
           --allow-extra-chr

echo "  Done!"
echo ""

# Step 4: Format PCs for SAIGE
echo "Step 4: Formatting PCs for SAIGE..."

# PLINK2 PCA output: FID IID PC1 PC2 ... PC10
# Need to create a file with: IID PC1 PC2 ... PC10

awk 'NR > 1 {printf "%s", $2; for(i=3; i<=NF; i++) printf "\t%s", $i; printf "\n"}' \
    ${OUTPUT_DIR}/pca_results.eigenvec > ${OUTPUT_DIR}/pcs.txt

# Add header
(echo -e "IID\t$(seq -s '\t' -f 'PC%.0f' 1 ${N_PCS})" && cat ${OUTPUT_DIR}/pcs.txt) \
    > ${OUTPUT_DIR}/pcs.with_header.txt

mv ${OUTPUT_DIR}/pcs.with_header.txt ${OUTPUT_DIR}/pcs.txt

echo "  Done!"
echo ""

echo "=========================================="
echo "PC Estimation Complete!"
echo "=========================================="
echo ""
echo "Generated files in ${OUTPUT_DIR}:"
echo "  - pcs.txt (PCs for SAIGE, with header)"
echo "  - pca_results.eigenvec (full PLINK PCA output)"
echo "  - pca_results.eigenval (eigenvalues)"
echo "  - ld_pruned.prune.in (LD-pruned variant list)"
echo ""
echo "PC summary:"
head -5 ${OUTPUT_DIR}/pcs.txt | column -t
echo "..."
echo ""
echo "To use these PCs in SAIGE:"
echo "  1. Merge with phenotype file"
echo "  2. Add 'PC1,PC2,...,PC10' to covariates list in step1"
echo ""
