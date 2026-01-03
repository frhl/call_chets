#!/usr/bin/env bash
# Prepare markers for categorical variance ratio estimation (via Docker)
#
# For set-based tests, SAIGE needs multiple variance ratios for different MAC categories:
# - Category 1: 10 <= MAC < 20
# - Category 2: MAC >= 20
#
# This script extracts random markers from each category for variance ratio estimation

set -euo pipefail

# Configuration
PLINK2_IMAGE="biocontainer/plink2:alpha2.3_jan2020"
INPUT_DIR="$(pwd)/input"
PLINK_PREFIX="simulated"
N_MARKERS_PER_CATEGORY=1000  # Number of random markers per MAC category

echo "=========================================="
echo "Preparing Variance Ratio Markers"
echo "=========================================="
echo "Using Docker image: ${PLINK2_IMAGE}"
echo ""

# Check if PLINK files exist
if [[ ! -f "${INPUT_DIR}/${PLINK_PREFIX}.bed" ]]; then
    echo "ERROR: PLINK files not found: ${INPUT_DIR}/${PLINK_PREFIX}.bed"
    echo "Run step 02_estimate_pcs.sh first to create PLINK files"
    exit 1
fi

# Step 1: Calculate allele counts
echo "Step 1: Calculating allele frequencies..."
docker run --rm \
    -v ${INPUT_DIR}:/data \
    -w /data \
    ${PLINK2_IMAGE} \
    plink2 --bfile ${PLINK_PREFIX} \
           --freq counts \
           --out ${PLINK_PREFIX} \
           --allow-extra-chr
echo "  Done!"
echo ""

# Step 2: Extract marker IDs for each MAC category
echo "Step 2: Extracting markers for MAC categories..."

# Check if we have the .acount file
ACOUNT_FILE="${INPUT_DIR}/${PLINK_PREFIX}.acount"
if [[ ! -f "${ACOUNT_FILE}" ]]; then
    echo "ERROR: Allele count file not found: ${ACOUNT_FILE}"
    exit 1
fi

# Category 1: 10 <= MAC < 20
echo "  Category 1 (10 <= MAC < 20)..."
N_CAT1=$(awk 'NR>1 {mac1=$5; mac2=2*$6-$5; if ((mac1 >= 10 && mac1 < 20) || (mac2 >= 10 && mac2 < 20)) print $2}' ${ACOUNT_FILE} | wc -l)
echo "    Found ${N_CAT1} markers in category 1"

if [ ${N_CAT1} -lt ${N_MARKERS_PER_CATEGORY} ]; then
    echo "    WARNING: Only ${N_CAT1} markers available, using all"
    N_TO_EXTRACT_CAT1=${N_CAT1}
else
    N_TO_EXTRACT_CAT1=${N_MARKERS_PER_CATEGORY}
fi

# Category 2: MAC >= 20
echo "  Category 2 (MAC >= 20)..."
N_CAT2=$(awk 'NR>1 {mac1=$5; mac2=2*$6-$5; if (mac1 >= 20 && mac2 >= 20) print $2}' ${ACOUNT_FILE} | wc -l)
echo "    Found ${N_CAT2} markers in category 2"

if [ ${N_CAT2} -lt ${N_MARKERS_PER_CATEGORY} ]; then
    echo "    WARNING: Only ${N_CAT2} markers available, using all"
    N_TO_EXTRACT_CAT2=${N_CAT2}
else
    N_TO_EXTRACT_CAT2=${N_MARKERS_PER_CATEGORY}
fi

# Extract random markers from each category
MARKER_LIST="${INPUT_DIR}/${PLINK_PREFIX}.vr_markers.list"

# Use sort -R (random sort) instead of shuf for macOS compatibility
cat <(awk 'NR>1 {mac1=$5; mac2=2*$6-$5; if ((mac1 >= 10 && mac1 < 20) || (mac2 >= 10 && mac2 < 20)) print $2}' ${ACOUNT_FILE} | sort -R | head -n ${N_TO_EXTRACT_CAT1}) \
    <(awk 'NR>1 {mac1=$5; mac2=2*$6-$5; if (mac1 >= 20 && mac2 >= 20) print $2}' ${ACOUNT_FILE} | sort -R | head -n ${N_TO_EXTRACT_CAT2}) \
    > ${MARKER_LIST}

N_TOTAL=$(wc -l < ${MARKER_LIST})
echo "  Total markers extracted: ${N_TOTAL}"
echo ""

# Step 3: Create PLINK file with variance ratio markers
echo "Step 3: Creating PLINK file for variance ratio estimation..."
docker run --rm \
    -v ${INPUT_DIR}:/data \
    -w /data \
    ${PLINK2_IMAGE} \
    plink2 --bfile ${PLINK_PREFIX} \
           --extract ${PLINK_PREFIX}.vr_markers.list \
           --make-bed \
           --out ${PLINK_PREFIX}_vr \
           --allow-extra-chr
echo "  Done!"
echo ""

echo "=========================================="
echo "Variance Ratio Markers Prepared!"
echo "=========================================="
echo ""
echo "Created files:"
echo "  - ${INPUT_DIR}/${PLINK_PREFIX}_vr.bed/bim/fam (PLINK file for variance ratio)"
echo "  - ${MARKER_LIST} (list of markers)"
echo ""
echo "Total markers: ${N_TOTAL}"
echo "  - Category 1 (10 <= MAC < 20): ${N_TO_EXTRACT_CAT1} markers"
echo "  - Category 2 (MAC >= 20): ${N_TO_EXTRACT_CAT2} markers"
echo ""
