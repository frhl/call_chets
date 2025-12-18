# Practical Examples

This directory contains scripts to generate synthetic data and run the `call_chets` pipeline in both phased and unphased modes. It also includes an R script to verify the orthogonality of the dominance encoding.

## 1. Generate Data
Run the generation script to create synthetic VCFs (`phased.vcf.gz`, `unphased.vcf.gz`), a gene map, and simulated phenotypes.
```bash
./00_generate_data.sh
```
*Requires `python3`.*

## 2. Phased Analysis Workflow
Demonstrates the standard pipeline:
1. Extract genotypes (simulated `bcftools query`).
2. Run `interpret_phase` to call compound heterozygotes.
3. specific `make_pseudo_vcf` to create `dominance` (orthogonal non-additive) and `additive` VCFs.
```bash
./01_phased_workflow.sh
```

## 3. Unphased Analysis Workflow
Demonstrates the unphased pipeline (burden/gene-level counts):
1. Extract genotypes.
2. Run `interpret_phase --unphased`.
3. Create VCF.
```bash
./02_unphased_workflow.sh
```

## 4. GWAS & Orthogonality Check (R)
Analyzes the output from the Phased workflow.
- Loads the simulated phenotypes (recessive trait).
- Loads the `additive` and `dominance` VCFs.
- Checks correlation (should be ~0 aka orthogonal).
- Runs a linear regression to show how the dominance term captures the recessive signal.
```bash
Rscript 03_run_gwas.R
```
*Requires `data.table` package in R.*

## Expected Results
Input data simulates a recessive effect (risk only when 2 alleles are present).
- **Correlation**: The additive and dominance encodings should have 0 correlation.
- **Regression**: Both Additive and Dominance terms should be significant. The Dominance term captures the deviation of the heterozygotes from the additive expectation, which is the hallmark of a recessive effect.
