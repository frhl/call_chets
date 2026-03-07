#!/usr/bin/env bash
# REGENIE Step 2: Test variants for association with the phenotype
# This step uses the predictions from Step 1 to test individual variants
# Runs both additive and dominance encodings

set -euo pipefail

# Input/Output paths
input_dir="input"
output_dir="output"
mkdir -p "${output_dir}"

# Input files
pheno_file="${input_dir}/phenotypes.regenie.txt"
pred_list="${output_dir}/regenie_step1_pred.list"
loco_file="${output_dir}/regenie_step1_1.loco"

# Phenotype and covariates
pheno_col="phenotype"
covariates="age,age2,sex,age_sex,age2_sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"
cat_covariates="sex"

# REGENIE parameters
bsize=400  # Block size for testing
minMAC=0.5  # Minimum minor allele count

# Docker image
docker_image="ghcr.io/rgcgithub/regenie/regenie:v4.1.gz"

# Encodings to test
encodings=("additive" "dominance")

echo "======================================"
echo "REGENIE Step 2: Variant Testing"
echo "======================================"
echo "Predictions: ${pred_list}"
echo "Phenotype: ${pheno_col}"
echo "Encodings: ${encodings[@]}"
echo ""

# Check that Step 1 outputs exist
if [[ ! -f "${pred_list}" ]]; then
    echo "ERROR: Step 1 prediction list not found: ${pred_list}"
    echo "Please run 01_regenie_step1.sh first"
    exit 1
fi

if [[ ! -f "${loco_file}" ]]; then
    echo "ERROR: Step 1 LOCO file not found: ${loco_file}"
    echo "Please run 01_regenie_step1.sh first"
    exit 1
fi

# Pull docker image if not present
echo "Checking for Docker image..."
if ! docker image inspect "${docker_image}" >/dev/null 2>&1; then
    echo "Pulling Docker image ${docker_image}..."
    docker pull "${docker_image}"
fi

# Loop through encodings
for encoding in "${encodings[@]}"; do

    echo ""
    echo "======================================"
    echo "Running REGENIE Step 2: ${encoding} encoding"
    echo "======================================"

    # Set input file based on encoding
    bgen_file="${output_dir}/simulated.${encoding}.bgen"
    sample_file="${output_dir}/simulated.${encoding}.sample"
    out_prefix="${output_dir}/regenie_step2.${encoding}"

    # Check if BGEN file exists
    if [[ ! -f "${bgen_file}" ]]; then
        echo "ERROR: BGEN file not found: ${bgen_file}"
        echo "Please run 01b_convert_vcf_to_bgen.sh first"
        exit 1
    fi

    echo "Input BGEN: ${bgen_file}"
    echo "Output: ${out_prefix}_${pheno_col}.regenie"

    # Run REGENIE Step 2
    docker run --rm \
      -v "$(pwd)/${input_dir}":/data/input:ro \
      -v "$(pwd)/${output_dir}":/data/output \
      "${docker_image}" \
      regenie \
        --step 2 \
        --bgen /data/output/simulated.${encoding}.bgen \
        --sample /data/output/simulated.${encoding}.sample \
        --ref-first \
        --phenoFile /data/input/phenotypes.regenie.txt \
        --phenoCol "${pheno_col}" \
        --covarFile /data/input/phenotypes.regenie.txt \
        --covarColList "${covariates}" \
        --catCovarList "${cat_covariates}" \
        --pred /data/output/regenie_step1_pred.list \
        --bsize ${bsize} \
        --minMAC ${minMAC} \
        --qt \
        --apply-rint \
        --out /data/output/regenie_step2.${encoding}

    echo "  Created: ${out_prefix}_${pheno_col}.regenie"

done

echo ""
echo "======================================"
echo "REGENIE Step 2 Complete!"
echo "======================================"
echo "Output files:"
for encoding in "${encodings[@]}"; do
    echo "  - output/regenie_step2.${encoding}_${pheno_col}.regenie"
done
echo ""
echo "Results contain variant-level association statistics for both encodings"
echo "Next step: Visualize results or compare with SAIGE output"
