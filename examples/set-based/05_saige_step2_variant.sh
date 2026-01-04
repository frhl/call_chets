#!/usr/bin/env bash
# SAIGE Step 2: Single variant association testing

set -euo pipefail

# Config
VCF="simulated.vcf.gz"
VCF_REC="simulated.recessive.vcf.gz"
VCF_NONADD="simulated.nonadditive.vcf.gz"
GRM=$(find output -name "sparseGRM*.mtx" | head -n 1 | xargs basename)

# Create sample list
awk 'NR>1 {print $1}' input/simulated.phenos.with_covariates.tsv > input/sample_list.txt

# Build command for standard (additive) SAIGE testing
CMD_ADD="step2_SPAtests.R \
    --vcfFile=/input/${VCF} \
    --vcfField=GT \
    --chrom=1 \
    --minMAF=0 \
    --minMAC=2 \
    --sampleFile=/input/sample_list.txt \
    --GMMATmodelFile=/output/null_model.rda \
    --varianceRatioFile=/output/null_model.varianceRatio.txt \
    --sparseGRMFile=/output/${GRM} \
    --sparseGRMSampleIDFile=/output/${GRM}.sampleIDs.txt \
    --SAIGEOutputFile=/output/saige.step2.additive.variant.txt \
    --LOCO=FALSE \
    --is_output_moreDetails=TRUE \
    --is_fastTest=TRUE"

# Run additive association test
docker run --rm -v "$(pwd)/input:/input" -v "$(pwd)/output:/output" \
    wzhou88/saige:0.5.1 ${CMD_ADD}

rm -f output/saige.step2.additive.variant.txt.index
echo "✓ Created: output/saige.step2.additive.variant.txt"

# Build command for standard (recessive [0,0,2]) SAIGE testing
CMD_REC="step2_SPAtests.R \
    --vcfFile=/input/${VCF_REC} \
    --vcfField=DS \
    --chrom=1 \
    --minMAF=0 \
    --minMAC=2 \
    --sampleFile=/input/sample_list.txt \
    --GMMATmodelFile=/output/null_model.rda \
    --varianceRatioFile=/output/null_model.varianceRatio.txt \
    --sparseGRMFile=/output/${GRM} \
    --sparseGRMSampleIDFile=/output/${GRM}.sampleIDs.txt \
    --SAIGEOutputFile=/output/saige.step2.recessive.variant.txt \
    --LOCO=FALSE \
    --is_output_moreDetails=TRUE \
    --is_fastTest=TRUE"

# Run additive association test
docker run --rm -v "$(pwd)/input:/input" -v "$(pwd)/output:/output" \
    wzhou88/saige:0.5.1 ${CMD_REC}

rm -f output/saige.step2.recessive.variant.txt.index
echo "✓ Created: output/saige.step2.recessive.variant.txt"


# Build command for non-additive SAIGE testing
CMD_NONADD="step2_SPAtests.R \
    --vcfFile=/input/${VCF_NONADD} \
    --vcfField=DS \
    --chrom=1 \
    --minMAF=0 \
    --minMAC=2 \
    --sampleFile=/input/sample_list.txt \
    --GMMATmodelFile=/output/null_model.rda \
    --varianceRatioFile=/output/null_model.varianceRatio.txt \
    --sparseGRMFile=/output/${GRM} \
    --sparseGRMSampleIDFile=/output/${GRM}.sampleIDs.txt \
    --SAIGEOutputFile=/output/saige.step2.nonadditive.variant.txt \
    --LOCO=FALSE \
    --is_output_moreDetails=TRUE \
    --dosage_zerod_MAC_cutoff=0 \
    --dosage_zerod_cutoff=0 \
    --SPAcutoff=0.5 \
    --pCutoffforFirth=0.10 \
    --is_Firth_beta=TRUE \
    --is_fastTest=FALSE"

# Run non-additive association test
docker run --rm -v "$(pwd)/input:/input" -v "$(pwd)/output:/output" \
    wzhou88/saige:0.5.1 ${CMD_NONADD}

rm -f output/saige.step2.nonadditive.variant.txt.index
echo "✓ Created: output/saige.step2.nonadditive.variant.txt"


