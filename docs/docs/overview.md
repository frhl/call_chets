## Overview

Standard genetic association studies assume **additive effects**, where each copy of an allele contributes equally to a phenotype. However, many genetic effects are **non-additive** -- for example, recessive variants only manifest when two copies are present.

**arcade** provides a suite of tools to detect and test for these non-additive effects.

### The model

We model the phenotype as a combination of additive and non-additive effects [^1] [^2] [^3] [^4] [^5]. Non-additive effects (also called dominance deviation) capture any deviation from a purely additive genetic model:

$$
y = X^A \beta^A + X^D \beta^D + \varepsilon
$$

where $y$ is the phenotype vector for $n$ individuals, $X^A$ and $X^D$ are the additive and non-additive genotype encoding matrices ($n \times m$ for $m$ loci), $\beta^A$ and $\beta^D$ are the corresponding effect size vectors, and $\varepsilon$ is the residual noise with variance $1 - h^2_A - h^2_D$.

The key requirement is that $X^A$ and $X^D$ are **orthogonal** ($X^A \cdot X^D = 0$), so that additive and non-additive effects can be estimated independently in a joint test.

### Genetic encodings

#### Under Hardy-Weinberg Equilibrium

Assuming HWE, the non-additive encoding maps genotypes as:

$$
[0, 1, 2] \rightarrow \left[-\frac{p}{1-p},\; 1,\; -\frac{1-p}{p}\right]
$$

where $p$ is the minor allele frequency [^1] [^2] [^3] [^4] [^5]. This encoding is uncorrelated with the additive encoding by construction.

#### General case (without assuming HWE)

To extend this approach beyond common variants [^1], **arcade** uses observed genotype proportions directly. Let $r$, $h$, and $a$ be the proportions of homozygous reference, heterozygous, and homozygous alternate genotypes respectively, with $r + h + a = 1$. At the gene level, these correspond to wildtype, monoallelic, and biallelic carrier proportions.

We start with the raw additive and non-additive vectors:

$$
{X^A}' = \begin{bmatrix}0\\1\\2\end{bmatrix}, \quad
{X^D}' = \begin{bmatrix}0\\1\\0\end{bmatrix}
$$

and define an inner product weighted by genotype proportions:

$$
\langle X, Y \rangle = r\, X_0 Y_0 + h\, X_1 Y_1 + a\, X_2 Y_2
$$

Applying the **Gram-Schmidt process** to orthogonalize these vectors yields the non-additive encoding:

$$
X^D = \frac{1}{\sqrt{har(a + r - (a - r)^2)}}
\begin{bmatrix}
-ha\\
2ar\\
-hr
\end{bmatrix}
$$

Under HWE (substituting $r = p^2$, $h = 2pq$, $a = q^2$), this recovers the standard HWE encoding above. In both cases, orthogonality holds: $X^A \cdot X^D = 0$. The resulting values are linearly rescaled to the $[0, 2]$ dosage range for compatibility with standard VCF tools, which preserves orthogonality.

### Encoding summary

| Encoding | Hom Ref (0/0) | Het (0/1) | Hom Alt (1/1) | Tests for |
|----------|:---:|:---:|:---:|-----------|
| **Additive** | 0 | 1 | 2 | Linear dose-response |
| **Non-additive** | 0 | d | 0 | Heterozygote deviation (orthogonalized) |
| **Recessive** | 0 | 0 | 2 | Recessive effects |

The **non-additive** encoding captures the deviation of heterozygotes from the additive expectation. The value $d$ depends on the allele frequency (or genotype proportions) and is computed by **arcade** automatically.

### Two pipelines

#### Gene-level pipeline (`interpret_phase` + `make_pseudo_vcf`)

For identifying compound heterozygotes and homozygous carriers within genes, then generating gene-level pseudo-variant VCFs for set-based or burden-style GWAS.

A **compound heterozygote** is an individual who carries two different rare variants on **separate haplotypes** within the same gene. These can mimic homozygous loss-of-function by disrupting both copies of a gene, even when neither variant alone is common enough to produce homozygotes. With phasing information, compound heterozygotes (variants on different haplotypes) can be distinguished from **in cis pairs** (variants on the same haplotype).

```
VCF → bcftools query → interpret_phase → make_pseudo_vcf → GWAS
```

1. Extract genotypes from a VCF using `bcftools query`
2. Run `interpret_phase` to identify compound hets and homozygous carriers per gene
3. Convert results to pseudo-variant VCFs using `make_pseudo_vcf` (additive + non-additive)
4. Run GWAS jointly: `Y ~ Additive + Non-additive`

#### Variant-level pipeline (`recode`)

For recoding individual variant genotypes for non-additive testing, without gene-level collapsing.

```
VCF → recode → GWAS
```

1. Take an existing VCF with standard genotypes
2. Run `recode` to produce non-additive or recessive encodings
3. Run GWAS jointly: `Y ~ Additive + Non-additive`

### Tools

| Tool | Description |
|------|-------------|
| `interpret_phase` | Identify compound heterozygous and homozygous variants from phased (or unphased) genotypes within gene regions |
| `make_pseudo_vcf` | Convert `interpret_phase` output into pseudo-variant biallelic VCFs with additive, non-additive, or recessive dosage encodings |
| `recode` | Orthogonalize or recode existing VCFs for non-additive or recessive genotype encodings |
| `filter_pp` | Filter VCF genotypes by posterior probability threshold |
| `count_by_gene` | Count genotypes per gene from a VCF and mapping file |

### Downstream integration

The output VCFs from **arcade** integrate directly with standard GWAS software:

- [**REGENIE**](https://github.com/rgcgithub/regenie) -- Variant-level association testing (see [REGENIE example](regenie.md))
- [**SAIGE**](https://github.com/saigegit/SAIGE) -- Variant-level and set-based burden testing (see [SAIGE example](saige.md))

### Related resources

- [Reproducibility repository](https://github.com/frhl/nonadditivity) -- Scripts to reproduce the analyses in the accompanying manuscript

### References

[^1]: Palmer DS, Zhou W, Abbott L, et al. Analysis of genetic dominance in the UK Biobank. *Science* 379(6639):1341-1348 (2023). [doi:10.1126/science.abn8455](https://doi.org/10.1126/science.abn8455)
[^2]: Vitezica ZG, Varona L, Legarra A. On the Additive and Dominant Variance and Covariance of Individuals Within the Genomic Selection Scope. *Genetics* 195(4):1223-1230 (2013). [doi:10.1534/genetics.113.155176](https://doi.org/10.1534/genetics.113.155176)
[^3]: Zhu Z, Bakshi A, Vinkhuyzen AAE, et al. Dominance Genetic Variation Contributes Little to the Missing Heritability for Human Complex Traits. *Am J Hum Genet* 96(3):377-385 (2015). [doi:10.1016/j.ajhg.2015.01.001](https://doi.org/10.1016/j.ajhg.2015.01.001)
[^4]: Pazokitoroudi A, Chiu AM, Burch KS, Pasaniuc B, Sankararaman S. Quantifying the contribution of dominance deviation effects to complex trait variation in biobank-scale data. *Am J Hum Genet* 108(5):799-808 (2021). [doi:10.1016/j.ajhg.2021.03.018](https://doi.org/10.1016/j.ajhg.2021.03.018)
[^5]: Hivert V, Sidorenko J, Rohart F, et al. Estimation of non-additive genetic variance in human complex traits from a large sample of unrelated individuals. *Am J Hum Genet* 108(5):786-798 (2021). [doi:10.1016/j.ajhg.2021.02.014](https://doi.org/10.1016/j.ajhg.2021.02.014)
