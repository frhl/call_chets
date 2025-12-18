#!/bin/bash
# Test Orthogonality of Additive and Dominance Encodings

set -e

# Setup paths
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT_DIR=$(dirname "$SCRIPT_DIR")
TEST_DIR=$(mktemp -d)
BIN_ORTH="$ROOT_DIR/bin/orthogonalize"
SCRIPT_CHECK="$SCRIPT_DIR/check_orthogonality.R"

echo "Running Orthogonality Test..."
echo "Temp Dir: $TEST_DIR"

# 1. Generate Toy Data
# We can use the examples/make_toy_data.py but tweak output path
cat <<EOF > $TEST_DIR/make_data.py
import gzip
import random
def write_vcf(filename):
    with gzip.open(filename, 'wt') as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("##contig=<ID=1,length=1000>\n")
        f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4\ts5\n")
        # 5 samples.
        # Var 1: 0, 1, 2, 0, 1
        # Genotypes: 0|0, 0|1, 1|1, 0|0, 0|1
        f.write("1\t100\t1:100:A:T\tA\tT\t.\tPASS\t.\tGT\t0|0\t0|1\t1|1\t0|0\t0|1\n")
        
        # Var 2 (Another freq): 0, 0, 0, 1, 2
        f.write("1\t200\t1:200:G:C\tG\tC\t.\tPASS\t.\tGT\t0|0\t0|0\t0|0\t1|0\t1|1\n")

if __name__ == "__main__":
    write_vcf("$TEST_DIR/input.vcf.gz")
EOF

python3 $TEST_DIR/make_data.py

# 2. Additive VCF is just the input
ADD_VCF="$TEST_DIR/input.vcf.gz"

# 3. Create Dominance VCF
DOM_VCF="$TEST_DIR/dominance.vcf.gz"
$BIN_ORTH --input $ADD_VCF --mode dominance --scale-per-variant | gzip > $DOM_VCF

# 4. Run Check
Rscript $SCRIPT_CHECK --additive $ADD_VCF --dominance $DOM_VCF

# Cleanup
rm -rf $TEST_DIR
echo "Test Passed."
