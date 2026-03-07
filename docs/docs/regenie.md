## REGENIE integration

This page demonstrates variant-level genetic association analysis using [REGENIE](https://github.com/rgcgithub/regenie) with both additive and dominance encodings produced by **arcade**.

The full example is in [`examples/regenie-variant-based/`](https://github.com/frhl/call_chets/tree/main/examples/regenie-variant-based).

### Overview

REGENIE uses a two-step approach:

1. **Step 1 (Null Model)**: Fits a whole-genome regression model using common variants to capture population structure and polygenic effects
2. **Step 2 (Variant Testing)**: Tests individual variants for association, using predictions from Step 1 to control for confounding

This example tests two genetic encodings:

- **Additive**: Standard 0/1/2 genotype encoding
- **Dominance**: Heterozygote deviation dosages from `recode`

### Quick start

```bash
cd examples/regenie-variant-based

# Run the complete pipeline
./run_all.sh
```

Or step by step:

```bash
./00_prepare_phenotypes.sh     # Format phenotypes for REGENIE
./01_regenie_step1.sh          # Step 1: fit null model
./01b_convert_vcf_to_bgen.sh   # Convert VCFs to BGEN format
./02_regenie_step2.sh          # Step 2: test both encodings
```

### Requirements

- **Docker** (REGENIE runs in a container)
- Docker image: `ghcr.io/rgcgithub/regenie/regenie:v4.1.gz`

### Step 0: Prepare phenotypes

REGENIE expects `FID` and `IID` columns. The script reformats the phenotype file:

```bash
./00_prepare_phenotypes.sh
```

Output columns: `FID IID phenotype age age2 sex age_sex age2_sex PC1-PC10`

### Step 1: Fit null model

```bash
./01_regenie_step1.sh
```

Key REGENIE parameters:

```bash
regenie \
    --step 1 \
    --bed simulated \
    --bsize 100 \
    --qt \
    --apply-rint \
    --threads 4
```

The genotypes are first filtered to `MAF > 0.01` using PLINK2. This step produces leave-one-chromosome-out (LOCO) predictions for use in Step 2.

### Step 1b: Convert VCF to BGEN

REGENIE Step 2 requires BGEN format. This script converts both encodings:

```bash
./01b_convert_vcf_to_bgen.sh
```

**Additive encoding** -- standard conversion from VCF `GT` field:

```bash
plink2 --vcf simulated.vcf.gz \
    --export bgen-1.3 'bits=16' ref-first
```

**Dominance encoding** -- uses `DS` field from the `recode`-transformed VCF:

```bash
plink2 --vcf simulated.nonadditive.vcf.gz dosage=DS \
    --import-dosage-certainty 1 \
    --hard-call-threshold 0 \
    --export bgen-1.3 'bits=16' ref-first
```

The dominance VCF is created by `recode`:

```bash
recode \
    --input simulated.vcf.gz \
    --mode nonadditive \
    --scale-globally \
    --min-hom-count 5 \
    --set-variant-id \
    --all-info \
    | bgzip > simulated.nonadditive.vcf.gz
```

### Step 2: Test variants

```bash
./02_regenie_step2.sh
```

Both encodings are tested in a loop:

```bash
for encoding in additive dominance; do
    regenie \
        --step 2 \
        --bgen simulated.${encoding}.bgen \
        --sample simulated.${encoding}.sample \
        --ref-first \
        --pred regenie_step1_pred.list \
        --bsize 400 \
        --minMAC 0.5 \
        --qt \
        --apply-rint
done
```

### Output format

Association results contain:

```
CHROM  GENPOS  ID       ALLELE0  ALLELE1  A1FREQ  INFO  N     TEST  BETA   SE     CHISQ  LOG10P
1      1234    var_0    A        G        0.016   1     5000  ADD   0.123  0.045  7.51   3.45
```

| Column | Description |
|--------|-------------|
| `CHROM` | Chromosome |
| `GENPOS` | Genomic position |
| `ID` | Variant ID |
| `A1FREQ` | Alternate allele frequency |
| `N` | Sample size |
| `TEST` | Test performed (`ADD` = additive) |
| `BETA` | Effect size estimate |
| `SE` | Standard error |
| `LOG10P` | -log10(p-value) |

### Adapting for real data

1. Replace input files with your own PLINK (`.bed/.bim/.fam`) and VCF (`.vcf.gz`) files
2. Update covariate lists (`covariates` and `cat_covariates` variables in scripts)
3. Adjust REGENIE parameters:
    - `--bsize`: Larger values use more memory but may be faster
    - `--minMAC`: Filter out very rare variants
    - `--bt`: Use for binary traits instead of `--qt`

### References

- Mbatchou et al. (2021) *Nature Genetics*. [doi:10.1038/s41588-021-00870-7](https://doi.org/10.1038/s41588-021-00870-7)
- REGENIE GitHub: [github.com/rgcgithub/regenie](https://github.com/rgcgithub/regenie)
