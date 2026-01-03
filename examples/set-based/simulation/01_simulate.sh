#!/usr/bin/env bash
# Simulate genotypes, gene boundaries, variant annotations, and phenotypes
# Using msprime/tskit for population genetic simulation

set -euo pipefail

# Configuration
CONDA_ENV="pop_sim"
N_SAMPLES=5000
TARGET_VARIANTS=10000
VARIANTS_PER_GENE=50
MAF_MIN=0.01
MAF_MAX=0.05
MAF_FRACTION=0.7
CAUSAL_GENE_FRACTION=0.1  # 10% of genes are causal
CAUSAL_VARIANTS_PER_GENE=5  # Average causal variants per causal gene
H2=0.10
ARCHITECTURE="additive"
SEED=42
OUTPUT_DIR="../input"
OUTPUT_PREFIX="simulated"
CHROMOSOME="1"

echo "=========================================="
echo "Running Genotype Simulation"
echo "=========================================="

# Run the simulation using conda environment
# Use -u for unbuffered Python output, --no-capture-output for conda
conda run --no-capture-output -n ${CONDA_ENV} python -u 01_simulate.py \
    --n_samples ${N_SAMPLES} \
    --target_variants ${TARGET_VARIANTS} \
    --variants_per_gene ${VARIANTS_PER_GENE} \
    --maf_min ${MAF_MIN} \
    --maf_max ${MAF_MAX} \
    --maf_fraction ${MAF_FRACTION} \
    --causal_gene_fraction ${CAUSAL_GENE_FRACTION} \
    --causal_variants_per_gene ${CAUSAL_VARIANTS_PER_GENE} \
    --h2 ${H2} \
    --architecture ${ARCHITECTURE} \
    --seed ${SEED} \
    --output_dir ${OUTPUT_DIR} \
    --output_prefix ${OUTPUT_PREFIX} \
    --chromosome ${CHROMOSOME}

# Compress VCF
if command -v bgzip &> /dev/null; then
    echo "Compressing VCF..."
    bgzip -f "${OUTPUT_DIR}/${OUTPUT_PREFIX}.vcf"
    echo "Created: ${OUTPUT_DIR}/${OUTPUT_PREFIX}.vcf.gz"
else
    echo "WARNING: bgzip not found. VCF remains uncompressed."
fi

# Index VCF
if command -v bcftools &> /dev/null; then
    echo "Indexing VCF..."
    bcftools index "${OUTPUT_DIR}/${OUTPUT_PREFIX}.vcf.gz"
    echo "Created: ${OUTPUT_DIR}/${OUTPUT_PREFIX}.vcf.gz.csi"
elif command -v tabix &> /dev/null; then
    echo "Indexing VCF (tabix)..."
    tabix -p vcf "${OUTPUT_DIR}/${OUTPUT_PREFIX}.vcf.gz"
    echo "Created: ${OUTPUT_DIR}/${OUTPUT_PREFIX}.vcf.gz.tbi"
else
    echo "WARNING: bcftools/tabix not found. VCF not indexed."
fi

echo ""
echo "=========================================="
echo "Simulation complete!"
echo "=========================================="
echo ""
echo "Generated files:"
echo "  - Genotypes (TSV): ${OUTPUT_DIR}/${OUTPUT_PREFIX}.genotypes.tsv.gz"
echo "  - Genotypes (VCF): ${OUTPUT_DIR}/${OUTPUT_PREFIX}.vcf.gz"
echo "  - Variant info: ${OUTPUT_DIR}/${OUTPUT_PREFIX}.variants.tsv"
echo "  - Phenotypes: ${OUTPUT_DIR}/${OUTPUT_PREFIX}.phenotypes.tsv"
echo "  - Gene boundaries: ${OUTPUT_DIR}/gene_boundaries.tsv"
echo "  - Variant annotations: ${OUTPUT_DIR}/variant_annotations.tsv"
echo ""


