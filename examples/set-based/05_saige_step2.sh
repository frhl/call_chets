#!/usr/bin/env bash
# SAIGE Step 2: Set-based association testing
#
# This step performs set-based (gene-level) burden tests using:
# - Null model from Step 1
# - VCF file with variants
# - Group file defining variant sets (genes)

set -euo pipefail

# Configuration
DOCKER_IMAGE="wzhou88/saige:0.5.1"
INPUT_DIR="$(pwd)/input"
STEP1_DIR="$(pwd)/output"
OUTPUT_DIR="$(pwd)/output"

# Input files
VCF_FILE="${INPUT_DIR}/simulated.vcf.gz"
GROUP_FILE="${INPUT_DIR}/genesets_all.txt"
NULL_MODEL="${STEP1_DIR}/null_model.rda"
VAR_RATIO="${STEP1_DIR}/null_model.varianceRatio.txt"
SAMPLE_FILE="/input/sample_list.txt"

# SAIGE parameters
USE_GRM=true        # Use GRM from step0
ANNOTATION="pLoF,synonymous"  # Annotation column to use
MAX_MAF=0.05       # Maximum MAF for variants to include
LOCO=FALSE         # Leave-one-chromosome-out (set FALSE for single chromosome)

# Create output directory
mkdir -p ${OUTPUT_DIR}

echo "=========================================="
echo "SAIGE Step 2: Set-Based Testing"
echo "=========================================="
echo ""

# Check required files
if [[ ! -f "${NULL_MODEL}" ]]; then
    echo "ERROR: Null model not found: ${NULL_MODEL}"
    echo "Run Step 1 first!"
    exit 1
fi



if [[ ! -f "${GROUP_FILE}" ]]; then
    echo "ERROR: Group file not found: ${GROUP_FILE}"
    echo "Run: python 01_prepare_saige_inputs.py"
    exit 1
fi

echo "Configuration:"
echo "  VCF file: ${VCF_FILE}"
echo "  Group file: ${GROUP_FILE}"
echo "  Null model: ${NULL_MODEL}"
echo "  Max MAF: ${MAX_MAF}"
echo "  Use GRM: ${USE_GRM}"
echo ""

# Create sample list if it doesn't exist
SAMPLE_LIST_HOST="${INPUT_DIR}/sample_list.txt"
if [[ ! -f "${SAMPLE_LIST_HOST}" ]]; then
    echo "Creating sample list from phenotypes..."
    awk 'NR>1 {print $1}' "${INPUT_DIR}/simulated.phenos.with_covariates.tsv" > "${SAMPLE_LIST_HOST}"
    echo "Created: ${SAMPLE_LIST_HOST}"
fi

# Check/Prepare VCF (Compress and Index)
VCF_RAW="${INPUT_DIR}/simulated.vcf" # Keep for reference or fallback logic if needed, though we use GZ
VCF_GZ="${INPUT_DIR}/simulated.vcf.gz"

if [[ -f "${VCF_GZ}" ]]; then
    echo "Found compressed VCF: ${VCF_GZ}"
elif [[ -f "${VCF_RAW}" ]]; then
    echo "Compressing VCF (required for SAIGE)..."
    if command -v bcftools &> /dev/null; then
        bcftools view -Oz -o "${VCF_GZ}" "${VCF_RAW}"
        echo "Created: ${VCF_GZ}"
    else
        echo "ERROR: bcftools not found. Cannot compress VCF."
        echo "Install bcftools or provide a bgzipped VCF."
        exit 1
    fi
else
    echo "ERROR: VCF file not found (checked ${VCF_GZ} and ${VCF_RAW})"
    exit 1
fi

if [[ ! -f "${VCF_GZ}.csi" ]]; then
    echo "Indexing VCF..."
    if command -v bcftools &> /dev/null; then
        bcftools index "${VCF_GZ}"
        echo "Index created: ${VCF_GZ}.csi"
    else
        echo "ERROR: bcftools not found. Cannot index VCF."
        exit 1
    fi
fi


# Build SAIGE command
SAIGE_CMD="step2_SPAtests.R \
    --vcfFile=/input/simulated.vcf.gz \
    --vcfField=GT \
    --chrom=1 \
    --minMAF=0 \
    --minMAC=0.5 \
    --sampleFile=${SAMPLE_FILE} \
    --GMMATmodelFile=/output/null_model.rda \
    --varianceRatioFile=/output/null_model.varianceRatio.txt \
    --groupFile=/input/genesets_all.txt \
    --annotation_in_groupTest=${ANNOTATION} \
    --maxMAF_in_groupTest=${MAX_MAF} \
    --SAIGEOutputFile=/output/results.txt \
    --LOCO=${LOCO} \
    --is_output_moreDetails=TRUE \
    --is_fastTest=TRUE"

# Add GRM if using
# Add GRM if using
if [[ "${USE_GRM}" == "true" ]]; then
    # Find the sparse GRM file (name varies based on parameters)
    SPARSE_GRM_FILE=$(find "${OUTPUT_DIR}" -name "sparseGRM*.mtx" | head -n 1)

    if [[ -z "${SPARSE_GRM_FILE}" ]]; then
        echo "WARNING: GRM requested but not found in ${OUTPUT_DIR}. Proceeding without GRM."
    else
        SPARSE_GRM_BASENAME=$(basename "${SPARSE_GRM_FILE}")
        SPARSE_GRM_SAMPLE_FILE="${SPARSE_GRM_FILE}.sampleIDs.txt"

        echo "Using GRM: ${SPARSE_GRM_BASENAME}"

        SAIGE_CMD="${SAIGE_CMD} \
            --sparseGRMFile=/output/${SPARSE_GRM_BASENAME} \
            --sparseGRMSampleIDFile=/output/${SPARSE_GRM_BASENAME}.sampleIDs.txt"
    fi
else
    echo "Not using GRM (USE_GRM=false)"
fi

echo "Running SAIGE Step 2..."
echo ""

# Run SAIGE Step 2
docker run --rm \
    -v ${INPUT_DIR}:/input \
    -v ${OUTPUT_DIR}:/output \
    ${DOCKER_IMAGE} \
    ${SAIGE_CMD}

echo ""
echo "=========================================="
echo "Step 2 Complete!"
echo "=========================================="
echo ""

# Clean up temporary files created by SAIGE
echo "Cleaning up temporary files..."
rm -f ${OUTPUT_DIR}/results.txt_P*Mat_Chunk_*.bin
rm -f ${OUTPUT_DIR}/results.txt.index
rm -f ${OUTPUT_DIR}/results.txt.singleAssoc.txt_temp
echo "  Removed temporary binary and index files"
echo ""

echo "Generated files in ${OUTPUT_DIR}:"
echo "  - results.txt (association test results)"
echo "  - results.txt.singleAssoc.txt (single variant results)"
echo ""
echo "Results contain p-values for each gene/set"
echo ""
