#!/usr/bin/env bash
# Simulate genotypes, gene boundaries, variant annotations, and phenotypes
# Using msprime/tskit for population genetic simulation

set -euo pipefail

# Configuration
CONDA_ENV="pop_sim"
N_SAMPLES=5000
N_GENES=100                      # Number of independent genes
VARIANTS_PER_GENE=200            # Variants per gene
GENE_LENGTH_MIN=30000            # Min gene length in bp
GENE_LENGTH_MAX=40000            # Max gene length in bp
INTERGENIC_SPACING=1000000       # Space between genes (1 Mb)
RECOMBINATION_RATE=1e-7          # Recombination rate per bp
MUTATION_RATE=2e-7               # Mutation rate per bp
POPULATION_SIZE=10000            # Effective population size
CAUSAL_GENE_FRACTION=0.07        # 7% of genes are causal
CAUSAL_VARIANTS_PER_GENE=1       # Causal variants per causal gene
PLOF_MAF_MIN=0.02                # Minimum MAF for pLoF variant selection
PLOF_MAF_MAX=0.05                # Maximum MAF for pLoF variant selection
H2=0.01
ARCHITECTURE="recessive"
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
    --n_genes ${N_GENES} \
    --variants_per_gene ${VARIANTS_PER_GENE} \
    --gene_length_min ${GENE_LENGTH_MIN} \
    --gene_length_max ${GENE_LENGTH_MAX} \
    --intergenic_spacing ${INTERGENIC_SPACING} \
    --recombination_rate ${RECOMBINATION_RATE} \
    --mutation_rate ${MUTATION_RATE} \
    --population_size ${POPULATION_SIZE} \
    --causal_gene_fraction ${CAUSAL_GENE_FRACTION} \
    --causal_variants_per_gene ${CAUSAL_VARIANTS_PER_GENE} \
    --pLoF_maf_min ${PLOF_MAF_MIN} \
    --pLoF_maf_max ${PLOF_MAF_MAX} \
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


