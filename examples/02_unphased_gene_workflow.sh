#!/bin/bash
set -e

# PATH setup
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ "$(basename "$SCRIPT_DIR")" == "examples" ]]; then
    ROOT_DIR="$(dirname "$SCRIPT_DIR")"
else
    ROOT_DIR="$SCRIPT_DIR"
fi

mkdir -p "$ROOT_DIR/examples/output/unphased_gene"
BIN_DIR="$ROOT_DIR/bin"
OUT_DIR="$ROOT_DIR/examples/output/unphased_gene"
IN_DIR="$ROOT_DIR/examples/input"

echo "Step 1: Extract Genotypes (Unphased)"
# We filter for GT="alt" to reduce file size.
if command -v bcftools &> /dev/null; then
    bcftools query -i'GT="alt"' -f'[%SAMPLE %CHROM:%POS:%REF:%ALT %GT\n]' "$IN_DIR/unphased.vcf.gz" | gzip > "$OUT_DIR/unphased_sites.txt.gz"
else
    echo "bcftools not found."
    exit 1
fi

echo "Step 2: Run interpret_phase (--unphased)"
# This treats all variants in a gene as potentially contributing, disregarding phase.
"$BIN_DIR/interpret_phase" \
    --geno "$OUT_DIR/unphased_sites.txt.gz" \
    --gene-map "$IN_DIR/gene_map.txt" \
    --unphased --show-variants \
    > "$OUT_DIR/unphased_results.txt"

echo "Step 3: Encode to VCF (Dominance)"
"$BIN_DIR/make_pseudo_vcf" \
    --input "$OUT_DIR/unphased_results.txt" \
    --samples "$IN_DIR/samples.txt" \
    --mode dominance \
    | bgzip > "$OUT_DIR/encoded_dominance.vcf.gz"

echo "Step 4: Create Additive VCF for comparison"
"$BIN_DIR/make_pseudo_vcf" \
    --input "$OUT_DIR/unphased_results.txt" \
    --samples "$IN_DIR/samples.txt" \
    --mode additive \
    | bgzip > "$OUT_DIR/encoded_additive.vcf.gz"

echo "Step 5: Run GWAS on Pseudo Variants"
Rscript "$ROOT_DIR/examples/run_gwas.R" \
    --additive "$OUT_DIR/encoded_additive.vcf.gz" \
    --dominance "$OUT_DIR/encoded_dominance.vcf.gz" \
    --phenotype "$IN_DIR/phenotypes.txt" \
    --out "$OUT_DIR/gwas_results"

echo "Done. Results in $OUT_DIR"
