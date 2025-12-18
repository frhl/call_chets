# Practical Examples

This directory contains workflows demonstrating the usage of `call_chets` and related tools.

## Folder Structure

*   `input/`: Contains generated synthetic data (VCFs, Gene Map, Phenotypes).
*   `output/`: Contains results from the workflows.
    *   `phased_gene/`: Results from the Phased Gene Workflow.
    *   `unphased_gene/`: Results from the Unphased Gene Workflow.
    *   `variant_level/`: Results from the Variant Level Workflow.

## 0. Generate Data
Generates synthetic VCFs (`phased.vcf.gz`, `unphased.vcf.gz`) and phenotypes (`Y_rec`, `Y_add`, `Y_null`) in `examples/input/`.
```bash
./00_generate_data.sh
```

## 1. Phased Gene Workflow
Standard pipeline for phased data. Calls compound heterozygotes and aggregates variants by gene.
1.  Extract genotypes.
2.  Run `interpret_phase` (call_chets).
3.  Create pseudo-variant VCFs (Additive & Dominance).
4.  Run GWAS on pseudo-variants (Genes).
```bash
./01_phased_gene_workflow.sh
```

## 2. Unphased Gene Workflow (Collapsing)
Pipeline for unphased data (burden/collapsing).
1.  Extract genotypes.
2.  Run `interpret_phase --unphased`.
3.  Create pseudo-variant VCFs.
4.  Run GWAS on pseudo-variants (Genes).
```bash
./02_unphased_gene_workflow.sh
```

## 3. Variant Level Workflow
Pipeline for analyzing individual variants directly using the `orthogonalize` tool.
1.  Prepare Additive VCF (standard).
2.  Create Dominance VCF using `orthogonalize --scale-per-variant`.
3.  Run GWAS on individual variants.
```bash
./03_variant_workflow.sh
```

## Analysis details
The GWAS step (`run_gwas.R`) performs a joint regression: `Y ~ Additive + Dominance`.
*   **Additive VCF**: Standard allele counts (0, 1, 2).
*   **Dominance VCF**: Encodes heterozygote deviation (orthogonalized).
Significant signals in both terms (especially Dominance) indicate non-additive/recessive effects.


### Dose-Response Estimation
For each variant, the script calculates the estimated genetic effect (G) for dosage levels 0, 1, and 2 by combining the additive and dominance contributions.

The effect for a genome dosage *d* is calculated as:
`Effect(d) = (Beta_Additive * d) + (Beta_Dominance * X_d)`

Where:
*   `Beta_Additive` and `Beta_Dominance` are the effect sizes estimated from the GWAS.
*   `X_d` is the specific dominance encoding value for that dosage level (empirically determined from the data).

Standard errors are calculated by combining the variances of both terms, assuming orthogonality. The output plots visualize these estimates with their 95% confidence intervals, helping to distinguish between Additive (linear) and Recessive (hockey-stick) architectures.
