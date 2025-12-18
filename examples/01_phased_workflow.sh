#!/bin/bash
set -e

# PATH setup (assuming run from root)
mkdir -p examples/output
BIN_DIR="bin"

echo "Step 1: Extract Phased Genotypes"
# In a real scenario, use bcftools to extract genotypes.
# Here we will simulate the extraction output directly from our file because we didn't install bcftools in this environment maybe?
# But assuming standard setup:
# bcftools view examples/phased.vcf.gz -Ou | bcftools query -f'[%SAMPLE %CHROM:%POS:%REF:%ALT %GT\n]' | gzip > examples/output/phased_sites.txt.gz

# Since we don't know if bcftools is installed in the user env, let's use a python helper to dump it exactly as bcftools would
python3 -c "
import gzip
n=0
with gzip.open('examples/phased.vcf.gz', 'rt') as f, gzip.open('examples/output/phased_sites.txt.gz', 'wt') as out:
    for line in f:
        if line.startswith('#'): continue
        parts = line.strip().split('\t')
        chrom, pos, ref, alt = parts[0], parts[1], parts[3], parts[4]
        # format is GT
        samples = parts[9:]
        # We need sample names from header
        pass
        # simpler: just rely on order if we knew it.
        # Actually checking if bcftools is installed is better.
"

# Let's try running bcftools, if it fails, fallback or error.
if command -v bcftools &> /dev/null; then
    bcftools query -f'[%SAMPLE %CHROM:%POS:%REF:%ALT %GT\n]' examples/phased.vcf.gz | gzip > examples/output/phased_sites.txt.gz
else
    echo "bcftools not found. Please install bcftools to run this step."
    exit 1
fi

echo "Step 2: Run interpret_phase (Compound Hets)"
# We need to map variants to genes.
$BIN_DIR/interpret_phase \
    --geno examples/output/phased_sites.txt.gz \
    --gene-map examples/gene_map.txt \
    --show-variants \
    > examples/output/chet_results.txt

echo "Step 3: Encode to VCF (Dominance/Non-additive)"
# This creates a VCF with orthogonal encoding
$BIN_DIR/make_pseudo_vcf \
    --input examples/output/chet_results.txt \
    --samples examples/samples.txt \
    --mode dominance \
    | bgzip > examples/output/encoded_dominance.vcf.gz

echo "Step 4: Create Additive VCF for comparison"
$BIN_DIR/make_pseudo_vcf \
    --input examples/output/chet_results.txt \
    --samples examples/samples.txt \
    --mode additive \
    | bgzip > examples/output/encoded_additive.vcf.gz

echo "Done. Results in examples/output/"
