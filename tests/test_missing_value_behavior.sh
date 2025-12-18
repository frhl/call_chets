#!/bin/bash
# Test how missing values are handled in recode

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
RECODE="${SCRIPT_DIR}/../bin/recode"

# Create VCF with missing values
cat << EOF > "${SCRIPT_DIR}/test_missing.vcf"
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1	Sample2	Sample3
chr1	100	.	A	T	.	PASS	.	GT	0/0	./.	1/1
EOF

echo "Testing missing value behavior..."
$RECODE --input "${SCRIPT_DIR}/test_missing.vcf" --mode recessive

rm "${SCRIPT_DIR}/test_missing.vcf"
