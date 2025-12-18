#!/bin/bash
set -e

# Determine the project root directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Check if we are in the examples directory or root
if [[ "$(basename "$SCRIPT_DIR")" == "examples" ]]; then
    ROOT_DIR="$(dirname "$SCRIPT_DIR")"
else
    ROOT_DIR="$SCRIPT_DIR"
fi

# Run the python script to generate data
# We cd to ROOT_DIR so that the relative paths inside the python script (examples/...) are valid
cd "$ROOT_DIR"
python3 examples/make_toy_data.py

echo "Generated: $ROOT_DIR/examples/input/phased.vcf.gz"
echo "Generated: $ROOT_DIR/examples/input/unphased.vcf.gz"
echo "Generated: $ROOT_DIR/examples/input/gene_map.txt"
echo "Generated: $ROOT_DIR/examples/input/phenotypes.txt"

# Create a sample list (just ID column) from the input VCF
gunzip -c "$ROOT_DIR/examples/input/phased.vcf.gz" | grep "#CHROM" | cut -f 10- | tr '\t' '\n' > "$ROOT_DIR/examples/input/samples.txt"
echo "Generated: $ROOT_DIR/examples/input/samples.txt"
