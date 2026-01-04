#!/usr/bin/env bash
# Encode simulated VCF with dominance mode using Docker
# This transforms the genotype dosages to orthogonal dominance encoding

set -euo pipefail

# Configuration
DOCKER_IMAGE="fhlassen/call_chets:1.0.4"
INPUT_VCF="input/simulated.vcf.gz"
OUTPUT_VCF="input/simulated.nonadditive.vcf.gz"
MODE="nonadditive"
MIN_HOM_COUNT=1  # Minimum number of homozygous alternate alleles required
SCALE_MODE="--scale-per-variant"  # Options: --scale-per-variant, --scale-globally, --scale-by-group <file>

echo "=========================================="
echo "Encoding VCF with nonadditive Mode"
echo "=========================================="
echo "Input VCF: ${INPUT_VCF}"
echo "Output VCF: ${OUTPUT_VCF}"
echo "Mode: ${MODE}"
echo "Min homozygous count: ${MIN_HOM_COUNT}"
echo "Scaling: ${SCALE_MODE}"
echo ""

# Check if input file exists
if [ ! -f "${INPUT_VCF}" ]; then
    echo "Error: Input VCF file not found: ${INPUT_VCF}"
    echo "Please run simulation/01_simulate.sh first."
    exit 1
fi

# Check if Docker is available
if ! command -v docker &> /dev/null; then
    echo "Error: Docker is not installed or not in PATH"
    exit 1
fi

# Get absolute path for Docker volume mounting
WORK_DIR=$(cd "$(dirname "${INPUT_VCF}")/.." && pwd)
INPUT_FILE=$(basename "${INPUT_VCF}")
OUTPUT_FILE=$(basename "${OUTPUT_VCF}")

echo "Running Docker container..."
echo ""

# Run encoding using Docker
# Mount the current directory to /data in the container
docker run --rm \
    -v "${WORK_DIR}:/data" \
    "${DOCKER_IMAGE}" \
    ./bin/recode \
        --input "/data/input/${INPUT_FILE}" \
        --mode "${MODE}" \
        --min-hom-count ${MIN_HOM_COUNT} \
        ${SCALE_MODE} \
        --set-variant-id \
        --all-info \
    | bgzip > "${OUTPUT_VCF}"

# Index the output VCF
if command -v bcftools &> /dev/null; then
    echo ""
    echo "Indexing output VCF..."
    bcftools index -f "${OUTPUT_VCF}"
    echo "Created: ${OUTPUT_VCF}.csi"
elif command -v tabix &> /dev/null; then
    echo ""
    echo "Indexing output VCF..."
    tabix -f -p vcf "${OUTPUT_VCF}"
    echo "Created: ${OUTPUT_VCF}.tbi"
else
    echo ""
    echo "WARNING: bcftools/tabix not found. VCF not indexed."
fi

echo ""
echo "=========================================="
echo "Encoding complete!"
echo "=========================================="
echo ""
echo "Output file: ${OUTPUT_VCF}"
echo ""
echo "Next steps:"
echo "  1. Run: ./02_prepare_vr.sh"
echo "  2. Run: ./03_saige_step0.sh"
echo ""
