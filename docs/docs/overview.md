## Overview

Standard genetic association studies assume **additive effects**, where each copy of an allele contributes equally to a phenotype. However, many genetic effects are **non-additive** -- for example, recessive variants only manifest when two copies are present.

**arcade** provides a suite of tools to detect and test for these non-additive effects.

### What are compound heterozygotes?

A **compound heterozygote** is an individual who carries two different rare variants on **separate haplotypes** within the same gene. If phasing information is available, compound heterozygotes (variants on different haplotypes) can be distinguished from **cis pairs** (variants on the same haplotype).

Compound heterozygotes can mimic homozygous loss-of-function by disrupting both copies of a gene, even when neither variant alone is common enough to produce homozygotes.

### Genetic encodings

**arcade** supports three genotype encodings:

| Encoding | Hom Ref (0/0) | Het (0/1) | Hom Alt (1/1) | Tests for |
|----------|:---:|:---:|:---:|-----------|
| **Additive** | 0 | 1 | 2 | Linear dose-response |
| **Dominance** | 0 | d | 0 | Heterozygote deviation (orthogonalized) |
| **Recessive** | 0 | 0 | 2 | Recessive effects |

The **dominance** encoding captures the deviation of heterozygotes from the additive expectation. It is orthogonalized so that additive and dominance effects can be estimated independently in a joint model (`Y ~ Additive + Dominance`).

### Two pipelines

#### Gene-level pipeline (`interpret_phase` + `make_pseudo_vcf`)

For identifying compound heterozygotes and homozygous carriers within genes, then generating gene-level pseudo-variant VCFs for set-based or burden-style GWAS.

```
VCF → bcftools query → interpret_phase → make_pseudo_vcf → GWAS
```

1. Extract genotypes from a VCF using `bcftools query`
2. Run `interpret_phase` to identify compound hets and homozygous carriers per gene
3. Convert results to pseudo-variant VCFs using `make_pseudo_vcf` (additive + dominance)
4. Run GWAS jointly: `Y ~ Additive + Dominance`

#### Variant-level pipeline (`recode`)

For recoding individual variant genotypes for non-additive testing, without gene-level collapsing.

```
VCF → recode → GWAS
```

1. Take an existing VCF with standard genotypes
2. Run `recode` to produce dominance or recessive encodings
3. Run GWAS jointly: `Y ~ Additive + Dominance`

### Tools

| Tool | Description |
|------|-------------|
| `interpret_phase` | Identify compound heterozygous and homozygous variants from phased (or unphased) genotypes within gene regions |
| `make_pseudo_vcf` | Convert `interpret_phase` output into pseudo-variant biallelic VCFs with additive, dominance, or recessive dosage encodings |
| `recode` | Orthogonalize or recode existing VCFs for non-additive (dominance) or recessive genotype encodings |
| `filter_pp` | Filter VCF genotypes by posterior probability threshold |
| `count_by_gene` | Count genotypes per gene from a VCF and mapping file |

### Downstream integration

The output VCFs from **arcade** integrate directly with standard GWAS software:

- [**REGENIE**](https://github.com/rgcgithub/regenie) -- Variant-level association testing (see [REGENIE example](regenie.md))
- [**SAIGE**](https://github.com/saigegit/SAIGE) -- Set-based burden testing (see [SAIGE example](saige.md))

### Related resources

- [Reproducibility repository](https://github.com/frhl/nonadditivity) -- Scripts to reproduce the analyses in the accompanying manuscript
