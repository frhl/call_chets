#!/bin/bash

# Test script for min-hom-count logic with minor allele check
# We want to ensure that variants are filtered if EITHER hom-ref OR hom-alt count is below threshold.

# Create input VCF
cat > tests/test_minor.vcf << EOL
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1	Sample2	Sample3	Sample4
chr1	100	VarLowRef	A	T	.	PASS	.	GT	0/1	0/1	1/1	1/1
chr1	200	VarLowAlt	A	T	.	PASS	.	GT	0/0	0/0	0/1	0/1
chr1	300	VarGoodie	A	T	.	PASS	.	GT	0/0	0/1	0/1	1/1
EOL

# VarLowRef: 0 hom-ref (aa), 2 het (Aa), 2 hom-alt (AA). min(0,2) = 0. Should be SKIPPED if min-hom-count=1.
# VarLowAlt: 2 hom-ref (aa), 2 het (Aa), 0 hom-alt (AA). min(2,0) = 0. Should be SKIPPED.
# VarGoodie: 1 hom-ref (aa), 2 het (Aa), 1 hom-alt (AA). min(1,1) = 1. Should be KEPT.

echo "Running recode with --min-hom-count 1..."
./bin/recode --input tests/test_minor.vcf --mode dominance --min-hom-count 1 > tests/output_minor.vcf 2> tests/recode.log

echo "Checking VarLowRef (Pos 100) (should be skipped)..."
if grep -q "chr1:100:A:T" tests/output_minor.vcf; then
  echo "FAIL: VarLowRef (Pos 100) was kept but has 0 hom-ref."
else
  echo "PASS: VarLowRef (Pos 100) was skipped."
fi

echo "Checking VarLowAlt (Pos 200) (should be skipped)..."
if grep -q "chr1:200:A:T" tests/output_minor.vcf; then
  echo "FAIL: VarLowAlt (Pos 200) was kept but has 0 hom-alt."
else
  echo "PASS: VarLowAlt (Pos 200) was skipped."
fi

echo "Checking VarGoodie (Pos 300) (should be kept)..."
if grep -q "chr1:300:A:T" tests/output_minor.vcf; then
  echo "PASS: VarGoodie (Pos 300) was kept."
else
  echo "FAIL: VarGoodie (Pos 300) was skipped but has 1 hom-ref and 1 hom-alt."
fi
