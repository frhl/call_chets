#!/usr/bin/env bash
# SAIGE Step 1: Fit null model with covariates

set -euo pipefail

# Config
PLINK="simulated_vr"  # Use variance ratio markers
PHENO="simulated.phenos.with_covariates.tsv"
COVARS="age,age2,sex,age_sex,age2_sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"
GRM=$(find output -name "sparseGRM*.mtx" | head -n 1 | xargs basename)

# Build command
CMD="step1_fitNULLGLMM.R \
    --plinkFile=/input/${PLINK} \
    --phenoFile=/input/${PHENO} \
    --phenoCol=phenotype \
    --traitType=quantitative \
    --covarColList=${COVARS} \
    --sparseGRMFile=/output/${GRM} \
    --sparseGRMSampleIDFile=/output/${GRM}.sampleIDs.txt \
    --useSparseGRMtoFitNULL=TRUE \
    --useSparseGRMforVarRatio=TRUE \
    --outputPrefix=/output/null_model \
    --nThreads=4 \
    --LOCO=FALSE \
    --isCateVarianceRatio=TRUE \
    --IsOverwriteVarianceRatioFile=TRUE"

# Fit null model
mkdir -p output
docker run --rm -v "$(pwd)/input:/input" -v "$(pwd)/output:/output" \
    wzhou88/saige:0.5.1 ${CMD}

echo "✓ Created: output/null_model.rda"
