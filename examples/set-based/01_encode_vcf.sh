#!/usr/bin/env bash
# Encode VCF to the dominance-deviation (orthogonal non-additive effects)

set -euo pipefail

# Config
MIN_HOM_COUNT=5  # Min minor homozygous alleles required per variant
MIN_HET_COUNT=10  # Min heterozygous alleles required per variant

# Paths
INPUT="input/simulated.vcf.gz"
OUTPUT_NONADD="input/simulated.nonadditive.vcf.gz"

# --min-het-count ${MIN_HET_COUNT} \
# Convert to nonadditive encoding using Docker 
docker run --rm -v "$(pwd):/data" fhlassen/call_chets:1.0.11 \
    ./bin/recode \
        --input "/data/${INPUT}" \
        --mode nonadditive \
        --min-hom-count ${MIN_HOM_COUNT} \
        --scale-per-variant \
        --set-variant-id \
        --all-info \
    | bgzip > "${OUTPUT_NONADD}"

# Index output
bcftools index -f "${OUTPUT_NONADD}" 2>/dev/null || tabix -f -p vcf "${OUTPUT_NONADD}" 2>/dev/null || true

echo "✓ Created: ${OUTPUT_NONADD}"

INPUT="input/simulated.vcf.gz"
OUTPUT_RECESSIVE="input/simulated.recessive.vcf.gz"

# Encode [0,1,2]->[0,0,2] using Docker 
docker run --rm -v "$(pwd):/data" fhlassen/call_chets:1.0.11 \
    ./bin/recode \
        --input "/data/${INPUT}" \
        --mode recessive \
        --min-hom-count ${MIN_HOM_COUNT} \
        --set-variant-id \
    | bgzip > "${OUTPUT_RECESSIVE}"

# Index output
bcftools index -f "${OUTPUT_RECESSIVE}" 2>/dev/null || tabix -f -p vcf "${OUTPUT_RECESSIVE}" 2>/dev/null || true

echo "✓ Created: ${OUTPUT_RECESSIVE}"



