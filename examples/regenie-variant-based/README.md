# REGENIE Variant-Level Analysis Example

This directory demonstrates a variant-level genetic association analysis using REGENIE with **both additive and dominance encodings**. The analysis uses the same simulated data from the `../saige-set-based/` example, allowing for direct comparison between SAIGE and REGENIE results.

## Overview

REGENIE is a method for rare and common variant association testing that uses a two-step approach:

1. **Step 1 (Null Model)**: Fits a whole-genome regression model using common variants to capture population structure and polygenic effects
2. **Step 2 (Variant Testing)**: Tests individual variants for association with the phenotype, using predictions from Step 1 to control for confounding

This example tests **two genetic encodings**:
- **Additive** (standard 0/1/2 genotype encoding)
- **Dominance** (tests for heterozygote effects using dosages)

## Quick Start

```bash
# Run the complete pipeline
./run_all.sh
```

Or run steps individually:

```bash
# Step 0: Prepare phenotype file in REGENIE format
./00_prepare_phenotypes.sh

# Step 1: Fit null model using whole-genome regression
./01_regenie_step1.sh

# Step 1b: Convert VCF to BGEN format (both encodings)
./01b_convert_vcf_to_bgen.sh

# Step 2: Test variants for association (both encodings)
./02_regenie_step2.sh
```

## Requirements

- **Docker**: REGENIE runs in a Docker container
- **Docker image**: `ghcr.io/rgcgithub/regenie/regenie:v4.1.gz`
- The scripts will automatically pull the image if not present

## Input Data

The analysis uses simulated data from `../saige-set-based/`:

- **Genotypes**: `simulated.bed/bim/fam` (PLINK format for Step 1)
- **Genotypes**: `simulated.vcf.gz` (VCF format for Step 2)
- **Phenotypes**: `simulated.phenos.with_covariates.tsv`
  - Quantitative trait: `phenotype`
  - Covariates: age, age², sex, age×sex, age²×sex, PC1-PC10
- **Sample list**: `sample_list.txt` (5000 samples)

## Directory Structure

```
regenie-variant-based/
├── README.md                    # This file
├── run_all.sh                   # Run complete pipeline
├── 00_prepare_phenotypes.sh     # Format phenotypes for REGENIE
├── 01_regenie_step1.sh          # REGENIE Step 1: null model
├── 01b_convert_vcf_to_bgen.sh   # Convert VCF to BGEN (both encodings)
├── 02_regenie_step2.sh          # REGENIE Step 2: variant testing (both encodings)
├── input/                       # Input data (symlinked from saige-set-based)
│   ├── simulated.bed/bim/fam   # Genotypes (PLINK format for Step 1)
│   ├── simulated.vcf.gz        # Genotypes additive (VCF format)
│   ├── simulated.nonadditive.vcf.gz  # Genotypes dominance (VCF format)
│   ├── simulated.phenos.with_covariates.tsv
│   ├── phenotypes.regenie.txt  # Created by 00_prepare_phenotypes.sh
│   └── sample_list.txt
└── output/                      # Results
    ├── regenie_step1.log       # Step 1 log
    ├── regenie_step1_pred.list # Prediction file list
    ├── regenie_step1_1.loco    # LOCO predictions
    ├── simulated.additive.bgen # BGEN for additive encoding
    ├── simulated.dominance.bgen # BGEN for dominance encoding
    ├── regenie_step2.additive_phenotype.regenie  # Additive results
    └── regenie_step2.dominance_phenotype.regenie # Dominance results
```

## Output Files

### Step 1 Output
- `regenie_step1.log`: Log file with model fitting details
- `regenie_step1_pred.list`: List of prediction files
- `regenie_step1_1.loco`: Leave-one-chromosome-out (LOCO) predictions
- `simulated.bed/bim/fam`: Filtered PLINK files (MAF > 0.01)

### Step 1b Output (BGEN Conversion)
- `simulated.additive.bgen`: BGEN file for additive encoding (20,000 variants)
- `simulated.additive.sample`: Sample file for additive BGEN
- `simulated.dominance.bgen`: BGEN file for dominance encoding (1,160 variants)
- `simulated.dominance.sample`: Sample file for dominance BGEN

### Step 2 Output
- `regenie_step2.additive_phenotype.regenie`: Additive encoding results (~20,000 variants)
- `regenie_step2.dominance_phenotype.regenie`: Dominance encoding results (~1,160 variants)

Association results contain:
```
CHROM  GENPOS  ID       ALLELE0  ALLELE1  A1FREQ   INFO  N      TEST   BETA      SE       CHISQ    LOG10P  EXTRA
1      1234    var_0    A        G        0.016    1     5000   ADD    0.123     0.045    7.51     3.45    ...
```

Columns:
- `CHROM`: Chromosome
- `GENPOS`: Genomic position
- `ID`: Variant ID
- `ALLELE0`: Reference allele
- `ALLELE1`: Alternate allele
- `A1FREQ`: Alternate allele frequency
- `INFO`: Imputation info score (1 for genotyped data)
- `N`: Sample size
- `TEST`: Test performed (ADD = additive)
- `BETA`: Effect size estimate
- `SE`: Standard error
- `CHISQ`: Chi-square statistic
- `LOG10P`: -log10(p-value)

## Genetic Encodings

This pipeline tests two different genetic encodings:

### Additive Encoding (Standard)
- **Genotype encoding**: 0/1/2 (number of alternate alleles)
- **VCF field**: `GT` (genotype field)
- **Use case**: Standard GWAS, tests for linear dose-response relationship
- **Interpretation**: Each copy of the alternate allele has the same effect

### Dominance Encoding
- **Genotype encoding**: 0/1/0 coded as dosages (DS field)
  - Homozygous reference (0/0) → 0
  - Heterozygous (0/1) → 1
  - Homozygous alternate (1/1) → 0
- **VCF field**: `DS` (dosage field)
- **Use case**: Tests for heterozygote-specific effects
- **Interpretation**: Heterozygotes have a distinct effect from homozygotes
- **Note**: Fewer variants tested due to lower MAC in dominance encoding

## REGENIE Parameters

### Step 1 (Null Model)
```bash
--step 1                # Step 1: fit null model
--bed simulated         # PLINK format genotypes (filtered, MAF > 0.01)
--bsize 100            # Block size for ridge regression
--qt                   # Quantitative trait
--apply-rint           # Apply rank-inverse normal transformation
--threads 4            # Number of threads
```

### Step 1b (BGEN Conversion)
```bash
# Additive encoding
plink2 --vcf simulated.vcf.gz \
  --export bgen-1.3 'bits=16' ref-first

# Dominance encoding
plink2 --vcf simulated.nonadditive.vcf.gz dosage=DS \
  --import-dosage-certainty 1 \
  --hard-call-threshold 0 \
  --export bgen-1.3 'bits=16' ref-first
```

### Step 2 (Variant Testing)
```bash
--step 2               # Step 2: test variants
--bgen simulated.{encoding}.bgen  # BGEN format genotypes
--sample {encoding}.sample         # Sample file
--ref-first            # Reference allele is first
--pred step1_pred.list # Predictions from Step 1
--bsize 400           # Block size for testing
--minMAC 0.5          # Minimum minor allele count
--qt                  # Quantitative trait
--apply-rint          # Apply rank-inverse normal transformation
```

## Comparison with SAIGE

This example uses the same simulated data as `../saige-set-based/`, but focuses on **variant-level** testing rather than **gene-based** testing:

| Feature | SAIGE (set-based) | REGENIE (variant-level) |
|---------|-------------------|-------------------------|
| **Test unit** | Gene/variant sets | Individual variants |
| **Step 1** | Fit null mixed model | Whole-genome regression |
| **Step 2** | Gene burden tests | Variant association tests |
| **Output** | Gene-level p-values | Variant-level p-values |
| **Use case** | Rare variant aggregation | Common/rare variants |

You can compare results from both methods:
- SAIGE variant results: `../saige-set-based/output/variant_results.txt`
- REGENIE variant results: `output/regenie_step2_phenotype.regenie`

## Advanced Usage

### Customizing Parameters

Edit the scripts to modify:
- **Block sizes**: `bsize` parameter (affects speed/memory)
- **Covariates**: Modify `covariates` and `cat_covariates` variables
- **Filters**: Adjust `minMAC` to change minor allele count threshold
- **Phenotypes**: Test different phenotypes by changing `pheno_col`

### Running on Real Data

To adapt this pipeline for real data:

1. Replace symbolic links in `input/` with your data:
   - PLINK files (`.bed/.bim/.fam`) for Step 1
   - VCF file (`.vcf.gz`) for Step 2
   - Phenotype file with FID, IID, phenotype, and covariates

2. Update covariate lists in scripts:
   - Modify `covariates` and `cat_covariates` variables
   - Ensure categorical covariates are properly coded

3. Adjust REGENIE parameters:
   - `--bsize`: Larger values use more memory but may be faster
   - `--minMAC`: Filter out very rare variants
   - `--bt`: Use for binary traits instead of `--qt`

## References

- **REGENIE**: Mbatchou et al. (2021) Nature Genetics
  - Paper: https://doi.org/10.1038/s41588-021-00870-7
  - GitHub: https://github.com/rgcgithub/regenie
  - Docker: https://github.com/rgcgithub/regenie/pkgs/container/regenie%2Fregenie

- **Simulated Data**: See `../saige-set-based/README.md` for simulation details

## Troubleshooting

### Docker Issues
- Ensure Docker is installed and running
- On Mac/Windows, allocate sufficient memory to Docker (≥8GB recommended)

### Step 2 Fails
- Make sure Step 1 completed successfully
- Check that `regenie_step1_pred.list` and `regenie_step1_1.loco` exist in `output/`

### Memory Issues
- Reduce `bsize` parameter
- Process fewer variants or samples

## Related Examples

- **SAIGE set-based**: `../saige-set-based/` - Gene-based burden testing
- **Comparison**: Both examples use the same simulated data for direct comparison
