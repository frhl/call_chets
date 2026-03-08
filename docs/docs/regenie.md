## REGENIE integration

This page demonstrates variant-level genetic association testing using [REGENIE](https://github.com/rgcgithub/regenie) with both additive and non-additive encodings produced by **arcade**.

The full example is in [`examples/regenie-variant-based/`](https://github.com/frhl/call_chets/tree/main/examples/regenie-variant-based).

> **Note**: REGENIE currently supports **variant-level testing only** with arcade's non-additive encodings. For set-based (gene-level) burden testing, see the [SAIGE integration](saige.md).

### Overview

REGENIE uses a two-step approach:

1. **Step 1 (Null Model)**: Fits a whole-genome regression model using common variants to capture population structure and polygenic effects
2. **Step 2 (Variant Testing)**: Tests individual variants for association, using predictions from Step 1 to control for confounding

This example tests two genetic encodings:

- **Additive**: Standard 0/1/2 genotype encoding
- **Non-additive**: Heterozygote deviation dosages from `recode`

### Quick start

```bash
cd examples/regenie-variant-based
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

### Step 1: Encode genotypes with `recode`

Create the non-additive VCF from a standard additive VCF:

```bash
recode \
    --input simulated.vcf.gz \
    --mode nonadditive \
    --scale-per-variant \
    --min-hom-count 5 \
    --set-variant-id \
    --all-info \
    | bgzip > simulated.nonadditive.vcf.gz
```

Note the use of `--scale-per-variant` for variant-level analysis, which scales each variant's non-additive dosage independently.

### Step 2: Convert to BGEN

REGENIE Step 2 requires BGEN format. Convert both encodings:

**Additive** — standard conversion from VCF `GT` field:

```bash
plink2 --vcf simulated.vcf.gz \
    --export bgen-1.3 'bits=16' ref-first
```

**Non-additive** — uses `DS` field from the `recode`-transformed VCF:

```bash
plink2 --vcf simulated.nonadditive.vcf.gz dosage=DS \
    --import-dosage-certainty 1 \
    --hard-call-threshold 0 \
    --export bgen-1.3 'bits=16' ref-first
```

### Step 3: Fit null model

```bash
regenie \
    --step 1 \
    --bed simulated \
    --bsize 100 \
    --qt \
    --apply-rint \
    --threads 4
```

### Step 4: Test variants

Both encodings are tested against the same null model:

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

### Adapting for real data

1. Replace input files with your own PLINK (`.bed/.bim/.fam`) and VCF (`.vcf.gz`) files
2. Update covariate lists in the scripts
3. Adjust REGENIE parameters:
    - `--bsize`: Larger values use more memory but may be faster
    - `--minMAC`: Filter out very rare variants

### References

- Mbatchou et al. (2021) *Nature Genetics*. [doi:10.1038/s41588-021-00870-7](https://doi.org/10.1038/s41588-021-00870-7)
- REGENIE GitHub: [github.com/rgcgithub/regenie](https://github.com/rgcgithub/regenie)
