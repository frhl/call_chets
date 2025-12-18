#!/bin/bash
set -e

# PATH setup
mkdir -p examples/output
BIN_DIR="bin"

echo "Step 1: Extract Genotypes (Unphased)"
if command -v bcftools &> /dev/null; then
    bcftools query -f'[%SAMPLE %CHROM:%POS:%REF:%ALT %GT\n]' examples/unphased.vcf.gz | gzip > examples/output/unphased_sites.txt.gz
else
    echo "bcftools not found."
    exit 1
fi

echo "Step 2: Run interpret_phase (--unphased)"
# This treats all variants in a gene as potentially contributing, disregarding phase.
$BIN_DIR/interpret_phase \
    --geno examples/output/unphased_sites.txt.gz \
    --gene-map examples/gene_map.txt \
    --unphased --show-variants \
    > examples/output/unphased_results.txt

echo "Step 3: Encode to VCF"
$BIN_DIR/make_pseudo_vcf \
    --input examples/output/unphased_results.txt \
    --samples examples/samples.txt \
    --mode additive \
    | bgzip > examples/output/encoded_unphased.vcf.gz

echo "Done."
