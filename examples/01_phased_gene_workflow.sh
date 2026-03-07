#!/bin/bash
set -e

# PATH setup
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ "$(basename "$SCRIPT_DIR")" == "examples" ]]; then
    ROOT_DIR="$(dirname "$SCRIPT_DIR")"
else
    ROOT_DIR="$SCRIPT_DIR"
fi

mkdir -p "$ROOT_DIR/examples/output/phased_gene"
BIN_DIR="$ROOT_DIR/bin"
OUT_DIR="$ROOT_DIR/examples/output/phased_gene"
IN_DIR="$ROOT_DIR/examples/input"

echo "Step 1: Extract Phased Genotypes"
# We filter for GT="alt" to only Output variants with at least one alternate allele (heterozygous or homozygous alt)
if command -v bcftools &> /dev/null; then
    bcftools query -i'GT="alt"' -f'[%SAMPLE %CHROM:%POS:%REF:%ALT %GT\n]' "$IN_DIR/phased.vcf.gz" | gzip > "$OUT_DIR/phased_sites.txt.gz"
else
    echo "bcftools not found. Please install bcftools to run this step."
    exit 1
fi

echo "Step 2: Run interpret_phase (Compound Hets)"
# We need to map variants to genes.
"$BIN_DIR/interpret_phase" \
    --geno "$OUT_DIR/phased_sites.txt.gz" \
    --gene-map "$IN_DIR/gene_map.txt" \
    --show-variants \
    > "$OUT_DIR/chet_results.txt"

echo "Step 3: Encode to VCF (Dominance/Non-additive)"
# This creates a VCF with orthogonal encoding
"$BIN_DIR/make_pseudo_vcf" \
    --input "$OUT_DIR/chet_results.txt" \
    --samples "$IN_DIR/samples.txt" \
    --mode dominance \
    | bgzip > "$OUT_DIR/encoded_dominance.vcf.gz"

echo "Step 4: Create Additive VCF for comparison"
"$BIN_DIR/make_pseudo_vcf" \
    --input "$OUT_DIR/chet_results.txt" \
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
