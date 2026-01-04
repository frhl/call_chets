#!/bin/bash

# Test script for min-het-count logic
# We want to ensure that variants are filtered if het count is below threshold.

# Create input VCF
cat > tests/test_het_count.vcf << EOL
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1	Sample2	Sample3	Sample4	Sample5	Sample6
chr1	100	NoHets	A	T	.	PASS	.	GT	0/0	0/0	0/0	1/1	1/1	1/1
chr1	200	OneHet	A	T	.	PASS	.	GT	0/0	0/0	0/1	1/1	1/1	1/1
chr1	300	TwoHets	A	T	.	PASS	.	GT	0/0	0/1	0/1	1/1	1/1	1/1
chr1	400	ThreeHets	A	T	.	PASS	.	GT	0/0	0/1	0/1	0/1	1/1	1/1
EOL

# NoHets: 3 hom-ref (aa), 0 het (Aa), 3 hom-alt (AA). Should be SKIPPED if min-het-count=1.
# OneHet: 2 hom-ref (aa), 1 het (Aa), 3 hom-alt (AA). Should be KEPT if min-het-count=1, SKIPPED if min-het-count=2.
# TwoHets: 1 hom-ref (aa), 2 het (Aa), 3 hom-alt (AA). Should be KEPT if min-het-count=2.
# ThreeHets: 1 hom-ref (aa), 3 het (Aa), 2 hom-alt (AA). Should be KEPT if min-het-count=2.

echo "========================================="
echo "Test 1: --min-het-count 1 (default)"
echo "========================================="
./bin/recode --input tests/test_het_count.vcf --mode dominance --min-het-count 1 > tests/output_het1.vcf 2> tests/recode_het1.log

echo "Checking NoHets (Pos 100) (should be skipped)..."
if grep -q "chr1:100:A:T" tests/output_het1.vcf; then
  echo "FAIL: NoHets (Pos 100) was kept but has 0 heterozygotes."
else
  echo "PASS: NoHets (Pos 100) was skipped."
fi

echo "Checking OneHet (Pos 200) (should be kept)..."
if grep -q "chr1:200:A:T" tests/output_het1.vcf; then
  echo "PASS: OneHet (Pos 200) was kept."
else
  echo "FAIL: OneHet (Pos 200) was skipped but has 1 heterozygote."
fi

echo "Checking TwoHets (Pos 300) (should be kept)..."
if grep -q "chr1:300:A:T" tests/output_het1.vcf; then
  echo "PASS: TwoHets (Pos 300) was kept."
else
  echo "FAIL: TwoHets (Pos 300) was skipped but has 2 heterozygotes."
fi

echo "Checking ThreeHets (Pos 400) (should be kept)..."
if grep -q "chr1:400:A:T" tests/output_het1.vcf; then
  echo "PASS: ThreeHets (Pos 400) was kept."
else
  echo "FAIL: ThreeHets (Pos 400) was skipped but has 3 heterozygotes."
fi

echo ""
echo "========================================="
echo "Test 2: --min-het-count 2"
echo "========================================="
./bin/recode --input tests/test_het_count.vcf --mode dominance --min-het-count 2 > tests/output_het2.vcf 2> tests/recode_het2.log

echo "Checking NoHets (Pos 100) (should be skipped)..."
if grep -q "chr1:100:A:T" tests/output_het2.vcf; then
  echo "FAIL: NoHets (Pos 100) was kept but has 0 heterozygotes."
else
  echo "PASS: NoHets (Pos 100) was skipped."
fi

echo "Checking OneHet (Pos 200) (should be skipped)..."
if grep -q "chr1:200:A:T" tests/output_het2.vcf; then
  echo "FAIL: OneHet (Pos 200) was kept but has only 1 heterozygote."
else
  echo "PASS: OneHet (Pos 200) was skipped."
fi

echo "Checking TwoHets (Pos 300) (should be kept)..."
if grep -q "chr1:300:A:T" tests/output_het2.vcf; then
  echo "PASS: TwoHets (Pos 300) was kept."
else
  echo "FAIL: TwoHets (Pos 300) was skipped but has 2 heterozygotes."
fi

echo "Checking ThreeHets (Pos 400) (should be kept)..."
if grep -q "chr1:400:A:T" tests/output_het2.vcf; then
  echo "PASS: ThreeHets (Pos 400) was kept."
else
  echo "FAIL: ThreeHets (Pos 400) was skipped but has 3 heterozygotes."
fi

# Clean up
rm -f tests/test_het_count.vcf tests/output_het1.vcf tests/output_het2.vcf tests/recode_het1.log tests/recode_het2.log
