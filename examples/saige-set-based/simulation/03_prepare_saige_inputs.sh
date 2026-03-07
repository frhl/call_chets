#!/usr/bin/env bash
# Prepare SAIGE input files from simulated data
#
# This script:
# 1. Creates SAIGE group file (variant-gene mapping)
# 2. Prepares phenotype file with covariates and PCs
# 3. Creates gene set files for SAIGE

set -euo pipefail

# Configuration
INPUT_DIR="$(cd "$(dirname "$0")/.." && pwd)/input"
GENE_BOUNDARIES="${INPUT_DIR}/gene_boundaries.tsv"
VARIANT_ANNOTATIONS="${INPUT_DIR}/variant_annotations.tsv"
PHENOTYPE_FILE="${INPUT_DIR}/simulated.phenos.tsv"
VCF_FILE="${INPUT_DIR}/simulated.vcf.gz"
OUTPUT_DIR="${INPUT_DIR}"
PCS_FILE="${INPUT_DIR}/pca/pcs.txt"
ANNOTATION_FILTER="all"  # Options: pLoF, synonymous, all
SEED=42

echo "=========================================="
echo "Preparing SAIGE Input Files"
echo "=========================================="
echo ""

# Check if input files exist
if [[ ! -f "${GENE_BOUNDARIES}" ]]; then
    echo "Error: Gene boundaries file not found: ${GENE_BOUNDARIES}"
    exit 1
fi

if [[ ! -f "${VARIANT_ANNOTATIONS}" ]]; then
    echo "Error: Variant annotations file not found: ${VARIANT_ANNOTATIONS}"
    exit 1
fi

if [[ ! -f "${PHENOTYPE_FILE}" ]]; then
    echo "Error: Phenotype file not found: ${PHENOTYPE_FILE}"
    exit 1
fi

# Run the Python script (using conda env)
conda run -n pop_sim python "$(dirname "$0")/03_prepare_saige_inputs.py" \
    --gene_boundaries "${GENE_BOUNDARIES}" \
    --variant_annotations "${VARIANT_ANNOTATIONS}" \
    --phenotype_file "${PHENOTYPE_FILE}" \
    --output_dir "${OUTPUT_DIR}" \
    --annotation_filter "${ANNOTATION_FILTER}" \
    --add_covariates \
    --pcs_file "${PCS_FILE}" \
    --vcf_file "${VCF_FILE}" \
    --seed ${SEED}

echo ""
echo "=========================================="
echo "Done!"
echo "=========================================="
echo ""
echo "Next steps:"
echo "  1. Review generated files in ${OUTPUT_DIR}/"
echo "  2. Run SAIGE step 1: ../03_saige_step1.sh"
echo "  3. Run SAIGE step 2: ../04_saige_step2.sh"
echo ""
