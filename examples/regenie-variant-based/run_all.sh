#!/usr/bin/env bash
# Run the complete REGENIE variant-level analysis pipeline

set -euo pipefail

echo "======================================"
echo "REGENIE Variant-Level Analysis"
echo "======================================"
echo ""
echo "This script will run:"
echo "  1. Prepare phenotype file"
echo "  2. REGENIE Step 1 (null model fitting)"
echo "  3. Convert VCF to BGEN (additive & dominance)"
echo "  4. REGENIE Step 2 (variant testing for both encodings)"
echo ""

# Step 0: Prepare phenotypes
echo "Step 0: Preparing phenotype file..."
./00_prepare_phenotypes.sh
echo ""

# Step 1: Fit null model
echo "Step 1: Fitting null model..."
./01_regenie_step1.sh
echo ""

# Step 1b: Convert VCF to BGEN
echo "Step 1b: Converting VCF to BGEN..."
./01b_convert_vcf_to_bgen.sh
echo ""

# Step 2: Test variants
echo "Step 2: Testing variants..."
./02_regenie_step2.sh
echo ""

echo "======================================"
echo "Pipeline Complete!"
echo "======================================"
echo ""
echo "Results are in the output/ directory:"
echo "  - regenie_step2.additive_phenotype.regenie"
echo "  - regenie_step2.dominance_phenotype.regenie"
echo ""
echo "You can now visualize the results or compare with SAIGE output from ../saige-set-based/"
