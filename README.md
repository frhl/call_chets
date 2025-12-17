# Call Compound Heterozygous Variants

Efficient C++ tools for detecting compound heterozygous variants and counting cis/trans configurations from phased genetic data.

## Tools
All binaries are compiled to `bin/`. Legacy names are preserved as symlinks but deprecated.
- `interpret_phase` (formerly `call_chets`): Core tool to identify compound heterozygous variants.
- `make_pseudo_vcf` (formerly `encode_vcf`): Converts phase results into VCF format with calculated dosages.
- `orthogonalize` (formerly `transform`): Transforms VCFs with orthogonal dominance encodings.
- `filter_pp` (formerly `filter_vcf_by_pp`): Filters VCFs based on posterior probabilities.

## Installation

### From Source
Requirements: `htslib` (bundled or system), `zlib`, `g++`, `make`.

```bash
make
sudo make install # Optional: installs to /usr/local/bin
```

### Docker
```bash
# Build
./build_docker.sh

# Run
docker run -it call_chets:latest
```

## Quick Start
*Note: Examples assume you ran `make install`. If not, prefix commands with `bin/` (e.g., `bin/interpret_phase`).*

### 1. Prepare Data
Extract phased genotypes from your VCF:
```bash
bcftools view input.vcf.gz -Ou | \
  bcftools query -f'[%SAMPLE %CHROM:%POS:%REF:%ALT %GT\n]' | \
  gzip > genotypes.txt.gz
```

### 2. Call Variants
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

### 3. Convert to VCF (Optional)
Create a VCF with dominance dosages:
```bash
make_pseudo_vcf \
  --input results.txt \
  --samples samples.txt \
  --mode dominance \
  --min-ac 1 | bgzip > output.vcf.gz
```

### 4. Orthogonalize (Optional)
You can also orthogonalize any VCF (legacy `transform` tool) directly:
```bash
orthogonalize \
  --input any_variant_file.vcf.gz \
  --mode dominance | bgzip > orthogonalized.vcf.gz
```