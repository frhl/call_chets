# call_chets

[![Tests](https://github.com/frhl/call_chets/actions/workflows/tests.yml/badge.svg)](https://github.com/frhl/call_chets/actions/workflows/tests.yml)
[![Docker](https://github.com/frhl/call_chets/actions/workflows/docker-image.yml/badge.svg)](https://github.com/frhl/call_chets/actions/workflows/docker-image.yml)

C++ tools for encoding biallelic variants in VCFs and performing allelic recoding for non-additive (dominance deviation) analysis.

## Tools

*   `interpret_phase`: Identify compound heterozygous (and homozygous) variants from phased genotypes.
*   `make_pseudo_vcf`: Convert phased results into a pseudo-variant (biallelic) VCF.
*   `recode`: Orthogonalize or recode VCFs for non-additive (`[0,1,2]->[-ha, -ar, -hr]`) or recessive (`[0,1,2]->[0,0,1]`) genotypes.

## Quick Start (Docker)

The easiest way to run the tools is via Docker.

```bash
# Pull the image
docker pull fhlassen/call_chets:latest

# Run example (mount current directory to /data)
docker run -v $PWD:/data fhlassen/call_chets interpret_phase --help
```

## Installation (Source)

Requirements: `htslib`, `zlib`, `g++`, `make`.

```bash
make
sudo make install
```

## Usage Workflow

Note: If running via Docker, prefix commands with `docker run -v $PWD:/data fhlassen/call_chets`.

### 1. Prepare Data
Extract phased genotypes from a VCF:
```bash
bcftools view tests/trio.vcf --max-af 0.01 -Ou | \
  bcftools query -i'GT="alt"' -f'[%SAMPLE %CHROM:%POS:%REF:%ALT %GT\n]' | \
  gzip > trio.phased_sites.txt.gz
```

### 2. Call Compound Hets
Run `interpret_phase` with your genotypes and a gene map:
```bash
interpret_phase \
  --geno genotypes.txt.gz \
  --gene-map tests/gene_map.txt \
  --verbose > results.txt
```

### 3. Convert to VCF
Create a VCF with non-additive dosages:
```bash
make_pseudo_vcf \
  --input results.txt \
  --samples tests/samples.txt \
  --mode dominance \
  --min-ac 1 | bgzip > output.vcf.gz
```

### 4. Recode Existing VCFs
Orthogonalize or recode any VCF:
```bash
recode \
  --input variants.vcf.gz \
  --mode dominance | bgzip > recoded.vcf.gz
```
