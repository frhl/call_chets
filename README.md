# Call Compound Heterozygous Variants

A collection of C++ tools for detecting compound heterozygous variants from phased genetic data.

## Installation

### From Source
1. Install required C++ libraries (see `Dockerfile` for dependencies)
2. Install [BCFtools](https://samtools.github.io/bcftools/howtos/install.html)
3. Compile:
   ```
   make
   ```

### Docker
```bash
# Build locally
./build_docker.sh

# Or pull from DockerHub
docker pull fhlassen/call_chets:0.3.3

# Run container
docker run -it call_chets:0.3.3
```

## Usage

### 1. Create Phased Sites File
Filter variants and extract phased genotypes:
```bash
bcftools view trio.vcf --max-af 0.01 -Ou | bcftools query -i'GT="alt"' -f'[%SAMPLE %CHROM:%POS:%REF:%ALT %GT\n]' | gzip > trio.phased_sites.txt.gz
```

### 2. Call Variants
Identify compound heterozygous variants per gene:
```bash
./call_chets \
  --geno trio.phased_sites.txt.gz \
  --gene-map gene_map.txt \
  --verbose > trio.result.txt
```

Optional: Add variant information with `--info-map`:
```bash
./call_chets \
  --geno trio.phased_sites.txt.gz \
  --gene-map gene_map.txt \
  --info-map info_map.txt > trio.result.txt
```

### 3. Generate VCF
Convert results to VCF format:
```bash
./encode_vcf \
  --input trio.result.txt \
  --samples samples.txt \
  --mode additive \
  --min-ac 1 | bgzip > trio.result.vcf.gz
```