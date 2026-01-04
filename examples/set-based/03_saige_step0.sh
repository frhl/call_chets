#!/usr/bin/env bash
# SAIGE Step 0: Generate GRM (Genetic Relationship Matrix) and sparse GRM
#
# This step is OPTIONAL for simulated unrelated samples.
# For real data with relatedness, this creates the GRM needed for mixed models.
#
# Note: For simplicity with simulated data, you can skip this step
# and not provide GRM to Step 1/2 (SAIGE will assume unrelated samples)

set -euo pipefail

# Configuration
DOCKER_IMAGE="wzhou88/saige:0.5.1"
INPUT_DIR="$(pwd)/input"
OUTPUT_DIR="$(pwd)/output"
PLINK_PREFIX="simulated"  # Assumes PLINK files exist

# Create output directory
mkdir -p ${OUTPUT_DIR}

echo "=========================================="
echo "SAIGE Step 0: Generate GRM"
echo "=========================================="
echo "NOTE: This step is OPTIONAL for simulated unrelated data"
echo ""

# Check if PLINK files exist (we need to convert VCF to PLINK first)
if [[ ! -f "${INPUT_DIR}/${PLINK_PREFIX}.bed" ]]; then
    echo "ERROR: PLINK files not found!"
    echo "You need to convert VCF to PLINK format first."
    echo ""
    echo "Using plink2:"
    echo "  plink2 --vcf ${INPUT_DIR}/simulated.vcf \\"
    echo "         --make-bed \\"
    echo "         --out ${INPUT_DIR}/${PLINK_PREFIX}"
    echo ""
    echo "Or skip Step 0 entirely for unrelated samples (recommended for simulated data)"
    exit 1
fi

echo "Running SAIGE Step 0..."
echo ""

# Run SAIGE createSparseGRM.R
docker run --rm \
    -v ${INPUT_DIR}:/input \
    -v ${OUTPUT_DIR}:/output \
    ${DOCKER_IMAGE} \
    createSparseGRM.R \
    --plinkFile=/input/${PLINK_PREFIX} \
    --nThreads=4 \
    --outputPrefix=/output/sparseGRM \
    --numRandomMarkerforSparseKin=1000 \
    --relatednessCutoff=0.125

echo ""
echo "=========================================="
echo "Step 0 Complete!"
echo "=========================================="
echo ""
echo "Generated files in ${OUTPUT_DIR}:"
echo "  - sparseGRM.mtx (sparse GRM matrix)"
echo "  - sparseGRM.mtx.sampleIDs.txt (sample IDs)"
echo ""
echo "These files will be used in Step 1 and Step 2"
echo ""
