#!/usr/bin/env bash
# SAIGE Step 1: Fit null model
#
# This step fits the null model (without genetic variants) to:
# - Estimate variance components
# - Calculate residuals
# - Prepare for Step 2 association testing

set -euo pipefail

# Configuration
DOCKER_IMAGE="wzhou88/saige:0.5.1"
INPUT_DIR="$(pwd)/input"
OUTPUT_DIR="$(pwd)/output"
PLINK_PREFIX="simulated"

# SAIGE parameters
PHENOTYPE_FILE="${INPUT_DIR}/simulated.phenos.with_covariates.tsv"
PHENOTYPE_COL="phenotype"
TRAIT_TYPE="quantitative"  # or "binary"
COVARIATES="age,age2,sex,age_sex,age2_sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"
USE_GRM=true  # Set to true if you ran step0

# Create output directory
mkdir -p ${OUTPUT_DIR}

echo "=========================================="
echo "SAIGE Step 1: Fit Null Model"
echo "=========================================="
echo ""

# Check if phenotype file exists
if [[ ! -f "${PHENOTYPE_FILE}" ]]; then
    echo "ERROR: Phenotype file not found: ${PHENOTYPE_FILE}"
    echo ""
    echo "Run: python 01_prepare_saige_inputs.py --add_covariates"
    exit 1
fi

# Check if PLINK files exist
if [[ ! -f "${INPUT_DIR}/${PLINK_PREFIX}.bed" ]]; then
    echo "ERROR: PLINK files not found: ${INPUT_DIR}/${PLINK_PREFIX}.bed"
    echo ""
    echo "Convert VCF to PLINK first using plink2"
    exit 1
fi

echo "Configuration:"
echo "  Phenotype file: ${PHENOTYPE_FILE}"
echo "  Trait type: ${TRAIT_TYPE}"
echo "  Covariates: ${COVARIATES:-none}"
echo "  Use GRM: ${USE_GRM}"
echo ""

# Check if variance ratio PLINK file exists
VR_PLINK="${PLINK_PREFIX}_vr"
if [[ ! -f "${INPUT_DIR}/${VR_PLINK}.bed" ]]; then
    echo "WARNING: Variance ratio PLINK file not found: ${INPUT_DIR}/${VR_PLINK}.bed"
    echo "Run: bash 02b_prepare_variance_ratio_markers.sh"
    echo "Using main PLINK file instead (may not work well for set-based tests)"
    VR_PLINK="${PLINK_PREFIX}"
else
    echo "✓ Variance ratio PLINK file detected: ${VR_PLINK}.bed"
    echo "  Will estimate categorical variance ratios for MAC categories"
fi
echo ""

# Build SAIGE command
SAIGE_CMD="step1_fitNULLGLMM.R \
    --plinkFile=/input/${VR_PLINK} \
    --phenoFile=/input/simulated.phenos.with_covariates.tsv \
    --phenoCol=${PHENOTYPE_COL} \
    --traitType=${TRAIT_TYPE} \
    --outputPrefix=/output/null_model \
    --nThreads=4 \
    --LOCO=FALSE \
    --isCateVarianceRatio=TRUE \
    --IsOverwriteVarianceRatioFile=TRUE"

# Add covariates if specified
if [[ -n "${COVARIATES}" ]]; then
    SAIGE_CMD="${SAIGE_CMD} --covarColList=${COVARIATES}"
fi

# Add GRM if using
# Add GRM if using
if [[ "${USE_GRM}" == "true" ]]; then
    # Find the sparse GRM file (name varies based on parameters)
    SPARSE_GRM_FILE=$(find "${OUTPUT_DIR}" -name "sparseGRM*.mtx" | head -n 1)

    if [[ -z "${SPARSE_GRM_FILE}" ]]; then
        echo "ERROR: GRM requested but not found in ${OUTPUT_DIR}. Run step0 first or set USE_GRM=false"
        exit 1
    else
        SPARSE_GRM_BASENAME=$(basename "${SPARSE_GRM_FILE}")

        echo "Using GRM: ${SPARSE_GRM_BASENAME}"

        SAIGE_CMD="${SAIGE_CMD} \
            --sparseGRMFile=/output/${SPARSE_GRM_BASENAME} \
            --sparseGRMSampleIDFile=/output/${SPARSE_GRM_BASENAME}.sampleIDs.txt \
            --useSparseGRMtoFitNULL=TRUE \
            --useSparseGRMforVarRatio=TRUE \
            --isCateVarianceRatio=TRUE"
    fi
else
    echo "Not using GRM (USE_GRM=false)"
fi

echo "Running SAIGE Step 1..."
echo ""

# Run SAIGE Step 1
docker run --rm \
    -v ${INPUT_DIR}:/input \
    -v ${OUTPUT_DIR}:/output \
    ${DOCKER_IMAGE} \
    ${SAIGE_CMD}

echo ""
echo "=========================================="
echo "Step 1 Complete!"
echo "=========================================="
echo ""
echo "Generated files in ${OUTPUT_DIR}:"
echo "  - null_model.rda (fitted null model)"
echo "  - null_model.varianceRatio.txt (variance ratio for Step 2)"
echo ""
echo "These files will be used in Step 2 for association testing"
echo ""
