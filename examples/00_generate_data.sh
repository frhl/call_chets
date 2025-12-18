#!/bin/bash
set -e

# Run the python script to generate data
python3 examples/make_toy_data.py

echo "Generated: examples/phased.vcf.gz"
echo "Generated: examples/unphased.vcf.gz"
echo "Generated: examples/gene_map.txt"
echo "Generated: examples/phenotypes.txt"

# Create a sample list (just ID column)
gunzip -c examples/phased.vcf.gz | grep "#CHROM" | cut -f 10- | tr '\t' '\n' > examples/samples.txt
echo "Generated: examples/samples.txt"
