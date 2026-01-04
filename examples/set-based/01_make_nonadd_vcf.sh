#!/usr/bin/env bash
# Encode VCF to dominance mode (orthogonal non-additive effects)

set -euo pipefail

# Config
MIN_HOM_COUNT=5  # Min homozygous alternate alleles required per variant

# Paths
INPUT="input/simulated.vcf.gz"
OUTPUT="input/simulated.nonadditive.vcf.gz"

# Encode using Docker
docker run --rm -v "$(pwd):/data" fhlassen/call_chets:latest \
    ./bin/recode \
        --input "/data/${INPUT}" \
        --mode nonadditive \
        --min-hom-count ${MIN_HOM_COUNT} \
        --scale-per-variant \
        --set-variant-id \
    | bgzip > "${OUTPUT}"

# Index output
bcftools index -f "${OUTPUT}" 2>/dev/null || tabix -f -p vcf "${OUTPUT}" 2>/dev/null || true

echo "✓ Created: ${OUTPUT}"
