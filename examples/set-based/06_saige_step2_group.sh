#!/usr/bin/env bash
# SAIGE Step 2: Set-based (gene-level) association testing

set -euo pipefail

# Config
VCF="simulated.vcf.gz"
VCF_NONADD="simulated.nonadditive.vcf.gz"
GROUPS="genesets_all.txt"
ANNOT="pLoF,synonymous"
MAX_MAF="0.001,0.05,0.1"
GRM=$(find output -name "sparseGRM*.mtx" | head -n 1 | xargs basename)

# Create sample list (if needed)
[[ -f input/sample_list.txt ]] || \
    awk 'NR>1 {print $1}' input/simulated.phenos.with_covariates.tsv > input/sample_list.txt

# Build command for standard (additive) SAIGE testing
CMD_ADD="step2_SPAtests.R \
    --vcfFile=/input/${VCF} \
    --vcfField=GT \
    --chrom=1 \
    --minMAF=0 \
    --minMAC=0.5 \
    --AlleleOrder=ref-first \
    --sampleFile=/input/sample_list.txt \
    --GMMATmodelFile=/output/null_model.rda \
    --varianceRatioFile=/output/null_model.varianceRatio.txt \
    --sparseGRMFile=/output/${GRM} \
    --sparseGRMSampleIDFile=/output/${GRM}.sampleIDs.txt \
    --groupFile=/input/${GROUPS} \
    --annotation_in_groupTest=${ANNOT} \
    --maxMAF_in_groupTest=${MAX_MAF} \
    --SAIGEOutputFile=/output/saige.step2.additive.group.txt \
    --LOCO=FALSE \
    --is_output_markerList_in_groupTest=TRUE \
    --is_output_moreDetails=TRUE \
    --is_fastTest=TRUE"

# Run additive gene-based test
docker run --rm -v "$(pwd)/input:/input" -v "$(pwd)/output:/output" \
    wzhou88/saige:0.5.1 ${CMD_ADD}

# Clean up temp files
rm -f output/saige.step2.additive.group.txt.{index,singleAssoc.txt,singleAssoc.txt_temp}
echo "✓ Created: output/saige.step2.additive.group.txt"

# Build command for non-additive SAIGE testing
# Notes: we use DS instead of GT. Also, that 
# we use the WEIGHTS in the group file, that we 
# have calculated manually.
CMD_NONADD="step2_SPAtests.R \
    --vcfFile=/input/${VCF_NONADD} \
    --vcfField=DS \
    --chrom=1 \
    --minMAF=0 \
    --minMAC=0.5 \
    --AlleleOrder=ref-first \
    --sampleFile=/input/sample_list.txt \
    --GMMATmodelFile=/output/null_model.rda \
    --varianceRatioFile=/output/null_model.varianceRatio.txt \
    --sparseGRMFile=/output/${GRM} \
    --sparseGRMSampleIDFile=/output/${GRM}.sampleIDs.txt \
    --groupFile=/input/${GROUPS} \
    --annotation_in_groupTest=${ANNOT} \
    --maxMAF_in_groupTest=0.50 \
    --SAIGEOutputFile=/output/saige.step2.nonadditive.group.txt \
    --LOCO=FALSE \
    --is_output_markerList_in_groupTest=TRUE \
    --is_output_moreDetails=TRUE \
    --is_fastTest=FALSE"

# Run non-additive gene-based test
docker run --rm -v "$(pwd)/input:/input" -v "$(pwd)/output:/output" \
    wzhou88/saige:0.5.1 ${CMD_NONADD}

# Clean up temp files
rm -f output/saige.step2.nonadditive.group.txt.{index,singleAssoc.txt,singleAssoc.txt_temp}
echo "✓ Created: output/saige.step2.nonadditive.group.txt"
