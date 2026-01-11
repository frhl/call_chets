#!/usr/bin/env bash
# SAIGE Step 0: Generate sparse GRM (OPTIONAL for unrelated samples)

set -euo pipefail

# Config
PLINK="simulated"

# Create sparse GRM
mkdir -p output
docker run --rm -v "$(pwd)/input:/input" -v "$(pwd)/output:/output" \
    wzhou88/saige:0.5.1 \
    createSparseGRM.R \
        --plinkFile=/input/${PLINK} \
        --nThreads=4 \
        --outputPrefix=/output/sparseGRM \
        --numRandomMarkerforSparseKin=1000 \
        --relatednessCutoff=0.125

echo "✓ Created: output/sparseGRM.mtx"
