# call_chets

[![Tests](https://github.com/frhl/call_chets/actions/workflows/tests.yml/badge.svg)](https://github.com/frhl/call_chets/actions/workflows/tests.yml)
[![Docker](https://github.com/frhl/call_chets/actions/workflows/docker-image.yml/badge.svg)](https://github.com/frhl/call_chets/actions/workflows/docker-image.yml)

C++ tools for encoding biallelic variants in VCFs and performing allelic recoding for non-additive (dominance deviation) analysis.

## Overview

Standard genetic association studies assume additive effects, where each copy of an allele contributes equally to a phenotype. However, many genetic effects are non-additive -- for example, recessive variants only manifest when two copies are present. `call_chets` provides a suite of tools to:

1. **Identify compound heterozygotes**: Detect individuals carrying two different rare variants on separate haplotypes within the same gene (compound hets) versus the same haplotype (cis pairs), using phased or unphased genotype data.
2. **Generate pseudo-variant VCFs**: Collapse gene-level compound het and homozygous calls into biallelic VCFs suitable for downstream GWAS tools.
3. **Recode genotypes for non-additive testing**: Orthogonalize standard VCFs into dominance deviation or recessive encodings for variant-level analysis.

The output VCFs integrate directly with standard GWAS software such as [REGENIE](https://github.com/rgcgithub/regenie) and [SAIGE](https://github.com/saigegit/SAIGE).

## Table of Contents

- [Tools](#tools)
- [Quick Start (Docker)](#quick-start-docker)
- [Installation (Source)](#installation-source)
- [Usage Workflow](#usage-workflow)
  - [Gene-Level Pipeline](#gene-level-pipeline-interpret_phase--make_pseudo_vcf)
  - [Variant-Level Pipeline](#variant-level-pipeline-recode)
- [Practical Examples](#practical-examples)
- [Detailed Tool Reference](#detailed-tool-reference)
- [Related Resources](#related-resources)

## Tools

| Tool | Description |
|------|-------------|
| `interpret_phase` | Identify compound heterozygous and homozygous variants from phased (or unphased) genotypes within gene regions. |
| `make_pseudo_vcf` | Convert `interpret_phase` output into pseudo-variant biallelic VCFs with additive, dominance, or recessive dosage encodings. |
| `recode` | Orthogonalize or recode existing VCFs for non-additive (`dominance`) or `recessive` genotype encodings. |

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

There are two main pipelines depending on whether you are working at the gene level or variant level.

> **Note**: If running via Docker, prefix commands with `docker run -v $PWD:/data fhlassen/call_chets`.

### Gene-Level Pipeline (`interpret_phase` + `make_pseudo_vcf`)

Use this pipeline to identify compound heterozygotes and homozygous carriers within genes, then generate gene-level pseudo-variant VCFs for set-based or burden-style GWAS.

#### 1. Prepare genotype input
Extract phased genotypes from a VCF:
```bash
bcftools view input.vcf.gz --max-af 0.01 -Ou | \
  bcftools query -i'GT="alt"' -f'[%SAMPLE %CHROM:%POS:%REF:%ALT %GT\n]' | \
  gzip > genotypes.txt.gz
```

#### 2. Call compound heterozygotes
Run `interpret_phase` with a variant-to-gene mapping file:
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

#### 3. Convert to pseudo-variant VCF
Create a VCF encoding gene-level dosages:
```bash
# Additive encoding (standard 0/1/2 dosages)
make_pseudo_vcf \
  --input results.txt \
  --samples samples.txt \
  --mode additive \
  --min-ac 1 | bgzip > additive.vcf.gz

# Dominance encoding (orthogonalized heterozygote deviation)
make_pseudo_vcf \
  --input results.txt \
  --samples samples.txt \
  --mode dominance \
  --min-ac 1 | bgzip > dominance.vcf.gz
```

#### 4. Run GWAS
Test both additive and dominance VCFs jointly (e.g. `Y ~ Additive + Dominance`) using your preferred GWAS tool.

### Variant-Level Pipeline (`recode`)

Use this pipeline to recode individual variant genotypes for non-additive testing, without gene-level collapsing.

```bash
# Dominance encoding (per-variant scaling)
recode \
  --input variants.vcf.gz \
  --mode dominance \
  --scale-per-variant | bgzip > dominance.vcf.gz

# Recessive encoding (het=0, hom_alt=2)
recode \
  --input variants.vcf.gz \
  --mode recessive | bgzip > recessive.vcf.gz
```

For gene-based scaling (comparable betas within genes):
```bash
recode \
  --input variants.vcf.gz \
  --mode dominance \
  --scale-by-group gene_map.txt | bgzip > dominance.vcf.gz
```

## Practical Examples

The [`examples/`](examples/) directory contains complete, runnable workflows with synthetic data. These demonstrate end-to-end pipelines from data preparation through GWAS:

| Example | Description |
|---------|-------------|
| [`00_generate_data.sh`](examples/00_generate_data.sh) | Generate synthetic phased/unphased VCFs and phenotypes |
| [`01_phased_gene_workflow.sh`](examples/01_phased_gene_workflow.sh) | Full gene-level pipeline with phased data (extract genotypes, call compound hets, make pseudo-VCFs, run GWAS) |
| [`02_unphased_gene_workflow.sh`](examples/02_unphased_gene_workflow.sh) | Gene-level pipeline for unphased data (burden/collapsing) |
| [`03_variant_workflow.sh`](examples/03_variant_workflow.sh) | Variant-level orthogonalization and GWAS |

Integration examples with popular GWAS tools:

| Example | Description |
|---------|-------------|
| [`examples/regenie-variant-based/`](examples/regenie-variant-based/) | Variant-level association testing with REGENIE (additive + dominance) |
| [`examples/saige-set-based/`](examples/saige-set-based/) | Gene-based burden testing with SAIGE (simulation, null model, group tests) |

Each integration example includes its own README with detailed instructions, parameter documentation, and expected outputs. To get started:

```bash
cd examples
./00_generate_data.sh           # generate synthetic data
./01_phased_gene_workflow.sh    # run the phased gene-level pipeline
```

See [`examples/README.md`](examples/README.md) for full details.

## Detailed Tool Reference

### `interpret_phase`

Identifies compound heterozygous (chet) and homozygous (hom) variant configurations within gene regions.

```
Usage: interpret_phase --geno <file> --gene-map <file> [options]

Required:
  --geno/-g <file>           Phased/unphased genotype file (gzipped)
  --gene-map/-m <file>       Variant-to-gene mapping file

Optional:
  --unphased                 Run in unphased mode (het/hom burden only)
  --info-map/-i <file>       Variant info file (AF, AC, etc.)
  --score-map/-p <file>      Variant score file for weighted collapsing
  --show-variants/-sv        Include detailed variant info in output
  --verbose/-v               Enable verbose logging

Score/Collapse Options (requires --score-map):
  --haplotype-collapse/-hc   Haplotype score rule [product|min|max|additive]
  --gene-collapse/-gc        Gene score rule [product|min|max|additive]
```

### `make_pseudo_vcf`

Converts `interpret_phase` output into a biallelic VCF with dosage encodings.

```
Usage: make_pseudo_vcf --input <file> --samples <file> --mode <mode>

Required:
  --input/-i <file>          Output from interpret_phase
  --samples/-s <file>        Sample list (one per line, no header)
  --mode/-m <mode>           Encoding mode:
                               additive / 012   - standard 0, 1, 2 dosages
                               dominance        - orthogonal heterozygote deviation
                               recessive / 001  - recessive (0 and 2 only)

Optional:
  --min-ac <n>               Minimum allele count filter (sum of DS >= n)
  --max-ac <n>               Maximum allele count filter (sum of DS < n)
  --all-info                 Include detailed INFO fields
```

### `recode`

Orthogonalizes or recodes an existing VCF for non-additive analysis.

```
Usage: recode --input <file.vcf.gz> [options]

Required:
  --input/-i <file>          Input VCF/BCF (.vcf, .vcf.gz, .bcf)

Mode:
  --mode/-m <mode>           dominance (default) or recessive

Scaling (mutually exclusive):
  --scale-per-variant        Scale each variant independently to [0,2]
  --scale-globally           Scale using global min/max across all variants
  --scale-by-group <file>    Scale within groups (tab-separated: variant, gene)

Filters:
  --min-hom-count <n>        Min minor homozygous count (default: 1)
  --min-het-count <n>        Min heterozygous count (default: 1)

Other:
  --set-variant-id           Set IDs to chr:pos:ref:alt format
  --all-info                 Include frequency/scaling info in output
```

## Related Resources

- **Reproducibility repository**: [github.com/frhl/nonadditivity](https://github.com/frhl/nonadditivity) -- scripts to reproduce the analyses presented in the accompanying manuscript.
