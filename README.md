# arcade

[![Tests](https://github.com/frhl/call_chets/actions/workflows/tests.yml/badge.svg)](https://github.com/frhl/call_chets/actions/workflows/tests.yml)
[![Docker](https://github.com/frhl/call_chets/actions/workflows/docker-image.yml/badge.svg)](https://github.com/frhl/call_chets/actions/workflows/docker-image.yml)

**ARCADE** — Allelic Recoding via gram-sChmidt for Dominance Effects

C++ tools for encoding biallelic variants in VCFs and performing allelic recoding for non-additive (dominance deviation) analysis.

## Tools

| Tool | Description |
|------|-------------|
| `interpret_phase` | Identify compound heterozygous and homozygous variants from phased or unphased genotypes within gene regions |
| `make_pseudo_vcf` | Convert `interpret_phase` output into pseudo-variant biallelic VCFs with additive, dominance, or recessive dosage encodings |
| `recode` | Orthogonalize or recode existing VCFs for non-additive (dominance) or recessive genotype encodings |

## Quick start

```bash
# Docker
docker pull fhlassen/arcade:latest
docker run -v $PWD:/data fhlassen/arcade interpret_phase --help

# From source
make
sudo make install
```

See [Install](https://frhl.github.io/call_chets/install/) for full installation instructions.

## Usage

There are two main pipelines:

**Gene-level** (`interpret_phase` + `make_pseudo_vcf`) — identify compound heterozygotes within genes and generate gene-level pseudo-variant VCFs for burden-style GWAS:

```bash
interpret_phase --geno genotypes.txt.gz --gene-map gene_map.txt > results.txt
make_pseudo_vcf --input results.txt --samples samples.txt --mode dominance | bgzip > dominance.vcf.gz
```

**Variant-level** (`recode`) — recode individual variant genotypes for non-additive testing:

```bash
recode --input variants.vcf.gz --mode dominance --scale-per-variant | bgzip > dominance.vcf.gz
```

Output VCFs integrate directly with [REGENIE](https://github.com/rgcgithub/regenie) and [SAIGE](https://github.com/saigegit/SAIGE).

## Documentation

Full documentation, CLI reference, and worked examples are available at **[frhl.github.io/call_chets](https://frhl.github.io/call_chets/)**.

## Examples

The [`examples/`](examples/) directory contains complete, runnable workflows with synthetic data:

- [`01_phased_gene_workflow.sh`](examples/01_phased_gene_workflow.sh) — Gene-level pipeline with phased data
- [`02_unphased_gene_workflow.sh`](examples/02_unphased_gene_workflow.sh) — Gene-level pipeline for unphased data
- [`03_variant_workflow.sh`](examples/03_variant_workflow.sh) — Variant-level orthogonalization

Integration examples with [REGENIE](examples/regenie-variant-based/) and [SAIGE](examples/saige-set-based/) are also provided.

## Citation

> Lassen, F.H. et al. (2025). *Genome-wide analysis of non-additive genetic effects*. (in preparation)

See also the [reproducibility repository](https://github.com/frhl/nonadditivity).

## License

[MIT](LICENSE)
