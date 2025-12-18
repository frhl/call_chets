#!/bin/bash
# Test Float Dosage Rounding Behavior

set -e

# Setup paths
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT_DIR=$(dirname "$SCRIPT_DIR")
TEST_DIR=$(mktemp -d)
BIN_ORTH="$ROOT_DIR/bin/orthogonalize"
echo "Temp Dir: $TEST_DIR"

# Create VCF with float dosages around 0.5 and 1.5
# 6 samples:
# S1: 0.49 -> Should round to 0
# S2: 0.50 -> Should round to 1 (usually)
# S3: 0.51 -> Should round to 1
# S4: 1.49 -> Should round to 1
# S5: 1.50 -> Should round to 2
# S6: 1.51 -> Should round to 2

cat <<EOF > $TEST_DIR/input.vcf
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=DS,Number=1,Type=Float,Description="Dosage">
##contig=<ID=1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4	S5	S6
1	100	1:100:A:T	A	T	.	PASS	.	DS	0.49	0.50	0.51	1.49	1.50	1.51
EOF

# Compress
gzip $TEST_DIR/input.vcf

# Run Orthogonalize in Dominance Mode
# We expect rounding to happen internally.
# Output will be dominance encoded.
# We also want to see if we can infer what the internal integer genotype was.
# Dominance encoding:
# 0 -> -h*a
# 1 -> 2*a*r
# 2 -> -h*r
# Since we have specific counts, frequencies will be determined by the rounded values.
# If S1=0, S2=1, S3=1, S4=1, S5=2, S6=2. (Proposed)
# Counts: 0:1, 1:3, 2:2. Total 6.
# r = 1/6 ~ 0.166
# h = 3/6 = 0.5
# a = 2/6 ~ 0.333

echo "Running Orthogonalize..."
$BIN_ORTH --input $TEST_DIR/input.vcf.gz --mode dominance --scale-per-variant | gzip > $TEST_DIR/output.vcf.gz

echo "Analyzing Output..."
# Decompress and look at DS values
gunzip -c $TEST_DIR/output.vcf.gz | grep -v "^#"

echo ""
echo "Recessive Mode Check (Clearer to see rounding)"
# Recessive mode: 0->0, 1->0, 2->2 (scaled? no, 2.0).
# If unscaled:
# 0 -> 0
# 1 -> 0
# 2 -> 2
# So:
# 0.49 -> 0 -> 0
# 0.50 -> 1 -> 0
# 0.51 -> 1 -> 0
# 1.49 -> 1 -> 0
# 1.50 -> 2 -> 2
# 1.51 -> 2 -> 2
# If S2 becomes 1, it gets output 0. If S5 becomes 2, it gets output 2.

$BIN_ORTH --input $TEST_DIR/input.vcf.gz --mode recessive | gzip > $TEST_DIR/output_rec.vcf.gz
gunzip -c $TEST_DIR/output_rec.vcf.gz | grep -v "^#"

# Cleanup
rm -rf $TEST_DIR
