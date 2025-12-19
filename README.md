# Create a VCF with nonadditive encoding

[![Tests](https://github.com/frhl/call_chets/actions/workflows/tests.yml/badge.svg)](https://github.com/frhl/call_chets/actions/workflows/tests.yml)
[![Docker](https://github.com/frhl/call_chets/actions/workflows/docker-image.yml/badge.svg)](https://github.com/frhl/call_chets/actions/workflows/docker-image.yml)
[![Docker Hub](https://img.shields.io/docker/pulls/frhl/call_chets.svg)](https://hub.docker.com/r/fhlassen/call_chets)

Efficient C++ tools for detecting biallelic (compound heterozygous or homozygous) variants and/or recoding them to either recessive or nonadditive genotypes.

## Tools
All binaries are compiled to `bin/`. Legacy names are preserved as symlinks but deprecated.
- `interpret_phase` (formerly `call_chets`): Core tool to identify compound heterozygous variants from phased data.
- `make_pseudo_vcf` (formerly `encode_vcf`): Converts phased results into a VCF.
- `recode` (formerly `orthogonalize`): Transforms VCFs with orthogonal dominance or recessive encodings.
- `filter_pp` (formerly `filter_vcf_by_pp`): Filters VCFs based on posterior probabilities.

## Installation

### From Source
Requirements: `htslib` (bundled or system), `zlib`, `g++`, `make`.

```bash
make
sudo make install # Optional: installs to /usr/local/bin
```

### Docker
You can pull the pre-built image from Docker Hub (replace `<USERNAME>` with your Docker Hub username):
```bash
docker pull <USERNAME>/call_chets:latest
# or for a specific version
docker pull <USERNAME>/call_chets:1.0.1
```

Or build locally:
```bash
docker build -t call_chets .
```

## Quick Start
*Note: Examples assume you ran `make install`. If not, prefix commands with `bin/` (e.g., `bin/interpret_phase`).*

### 1/3. Prepare Data
Extract phased genotypes from your VCF:
```bash
bcftools view trio.vcf --max-af 0.01 -Ou | \
  bcftools query -i'GT="alt"' -f'[%SAMPLE %CHROM:%POS:%REF:%ALT %GT\n]' | \
  gzip > trio.phased_sites.txt.gz
```

### 2/3. Call Variants
Run `interpret_phase` using a gene mapping file:
```bash
interpret_phase \
  --geno genotypes.txt.gz \
  --gene-map gene_map.txt \
  --verbose > results.txt
```

For unphased data, add the `--unphased` flag:
```bash
interpret_phase \
  --geno genotypes.txt.gz \
  --gene-map gene_map.txt \
  --unphased > results.txt
```

### 3/3. Convert to VCF
Create a VCF with nonadditive (domiannce) dosages:
```bash
make_pseudo_vcf \
  --input results.txt \
  --samples samples.txt \
  --mode dominance \
  --min-ac 1 | bgzip > output.vcf.gz
```

### 4. Recode
You can also recode any VCF directly (either with or without phased information) into orthogonal or recessive encodings:
```bash
recode \
  --input any_variant_file.vcf.gz \
  --mode dominance | bgzip > recoded.vcf.gz
```
