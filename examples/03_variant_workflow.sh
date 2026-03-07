#!/bin/bash
set -e

# PATH setup
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ "$(basename "$SCRIPT_DIR")" == "examples" ]]; then
    ROOT_DIR="$(dirname "$SCRIPT_DIR")"
else
    ROOT_DIR="$SCRIPT_DIR"
fi

mkdir -p "$ROOT_DIR/examples/output/variant_level"
BIN_DIR="$ROOT_DIR/bin"
OUT_DIR="$ROOT_DIR/examples/output/variant_level"
IN_DIR="$ROOT_DIR/examples/input"

echo "Step 1: Prepare Additive VCF"
# Input is already additive. Copy it.
cp "$IN_DIR/phased.vcf.gz" "$OUT_DIR/encoded_additive.vcf.gz"

echo "Step 2: Orthogonalize to Dominance VCF"
# Apply orthogonalize directly to variants
# Using --scale-per-variant as we are testing variants individually
"$BIN_DIR/orthogonalize" \
    --input "$IN_DIR/phased.vcf.gz" \
    --mode dominance \
    --scale-per-variant \
    | bgzip > "$OUT_DIR/encoded_dominance.vcf.gz"

echo "Step 3: Run GWAS on Variants"
Rscript "$ROOT_DIR/examples/run_gwas.R" \
    --additive "$OUT_DIR/encoded_additive.vcf.gz" \
    --dominance "$OUT_DIR/encoded_dominance.vcf.gz" \
    --phenotype "$IN_DIR/phenotypes.txt" \
    --out "$OUT_DIR/gwas_results"

echo "Done. Results in $OUT_DIR"
