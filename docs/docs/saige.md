## SAIGE integration

This page demonstrates genetic association testing using [SAIGE](https://github.com/saigegit/SAIGE) with genotype encodings produced by **arcade**. SAIGE supports both **variant-level** and **set-based (gene-level)** testing.

The full example is in [`examples/saige-set-based/`](https://github.com/frhl/call_chets/tree/main/examples/saige-set-based).

### Overview

The SAIGE pipeline has three main stages:

1. **Encode** — Use `recode` to produce non-additive and recessive VCFs from a standard additive VCF
2. **Null model** — Fit SAIGE null model with sparse GRM
3. **Association testing** — Run variant-level and/or gene-based (set-based) burden tests

Three genetic encodings are tested:

- **Additive**: Standard 0/1/2 genotype (`GT` field)
- **Non-additive**: Orthogonalized heterozygote deviation (`DS` field from `recode`)
- **Recessive**: Binary recessive encoding, 0/0/2 (`DS` field from `recode`)

### Quick start

```bash
cd examples/saige-set-based

# Encode VCFs
./01_encode_vcf.sh

# Prepare variance ratio markers
./02_prepare_vr.sh

# Run SAIGE pipeline
./03_saige_step0.sh              # Sparse GRM
./04_saige_step1.sh              # Null model
./05_saige_step2_variant.sh      # Variant-level tests
./06_saige_step2_group.sh        # Gene-based (set-based) tests
```

### Requirements

- **Docker** (SAIGE runs in a container)
- Docker image: `wzhou88/saige:0.5.1`

### Step 1: Encode VCFs with `recode`

Create non-additive and recessive VCFs from the standard additive input:

**Non-additive encoding:**

```bash
recode \
    --input simulated.vcf.gz \
    --mode nonadditive \
    --min-hom-count 5 \
    --max-maf 0.05 \
    --scale-globally \
    --set-variant-id \
    --all-info \
    | bgzip > simulated.nonadditive.vcf.gz
```

**Recessive encoding:**

```bash
recode \
    --input simulated.vcf.gz \
    --mode recessive \
    --min-hom-count 5 \
    --set-variant-id \
    | bgzip > simulated.recessive.vcf.gz
```

Note the use of `--scale-globally` for the non-additive encoding — this ensures comparable dosages across variants, which is required for set-based analysis. For variant-level only analysis, `--scale-per-variant` is faster.

### Step 2: Fit null model

First create a sparse GRM, then fit the null model:

```bash
# Sparse GRM
createSparseGRM.R \
    --plinkFile=simulated \
    --nThreads=4 \
    --outputPrefix=sparseGRM \
    --numRandomMarkerforSparseKin=1000 \
    --relatednessCutoff=0.125

# Null model
step1_fitNULLGLMM.R \
    --plinkFile=simulated_vr \
    --phenoFile=simulated.phenos.with_covariates.tsv \
    --phenoCol=phenotype \
    --traitType=quantitative \
    --covarColList=age,age2,sex,age_sex,age2_sex,PC1,...,PC10 \
    --sparseGRMFile=sparseGRM.mtx \
    --useSparseGRMtoFitNULL=TRUE \
    --useSparseGRMforVarRatio=TRUE \
    --LOCO=FALSE \
    --isCateVarianceRatio=TRUE
```

To test multiple phenotypes, change `--phenoCol` to the column name of the phenotype you wish to test. Each phenotype requires its own null model.

### Step 3: Variant-level testing

All three encodings are tested against the same null model:

**Additive** (standard `GT` field):

```bash
step2_SPAtests.R \
    --vcfFile=simulated.vcf.gz \
    --vcfField=GT \
    --minMAC=2 \
    --SAIGEOutputFile=saige.step2.additive.variant.txt
```

**Recessive** (`DS` field, 0/0/2 encoding):

```bash
step2_SPAtests.R \
    --vcfFile=simulated.recessive.vcf.gz \
    --vcfField=DS \
    --minMAC=2 \
    --SAIGEOutputFile=saige.step2.recessive.variant.txt
```

**Non-additive** (`DS` field, orthogonalized):

```bash
step2_SPAtests.R \
    --vcfFile=simulated.nonadditive.vcf.gz \
    --vcfField=DS \
    --minMAC=2 \
    --SAIGEOutputFile=saige.step2.nonadditive.variant.txt
```

### Step 4: Set-based (gene-level) testing

Gene-based burden tests use a group file that maps variants to genes with annotations. Both additive and non-additive encodings are tested:

**Additive burden** (`GT` field, all annotations):

```bash
step2_SPAtests.R \
    --vcfFile=simulated.vcf.gz \
    --vcfField=GT \
    --groupFile=genesets_all.txt \
    --annotation_in_groupTest=pLoF,synonymous \
    --maxMAF_in_groupTest=0.50 \
    --r.corr=1 \
    --SAIGEOutputFile=saige.step2.additive.group.txt
```

**Non-additive burden** (`DS` field, tested per annotation):

```bash
# pLoF variants
step2_SPAtests.R \
    --vcfFile=simulated.nonadditive.vcf.gz \
    --vcfField=DS \
    --groupFile=genesets_all.txt \
    --annotation_in_groupTest=pLoF \
    --maxMAF_in_groupTest=0.50 \
    --r.corr=1 \
    --SAIGEOutputFile=saige.step2.nonadditive.group.pLoF.txt

# Synonymous variants
step2_SPAtests.R \
    --vcfFile=simulated.nonadditive.vcf.gz \
    --vcfField=DS \
    --groupFile=genesets_all.txt \
    --annotation_in_groupTest=synonymous \
    --maxMAF_in_groupTest=0.50 \
    --r.corr=1 \
    --SAIGEOutputFile=saige.step2.nonadditive.group.synonymous.txt
```

> **Note**: For the non-additive encoding, annotations (e.g. pLoF, synonymous) should be tested separately due to a labelling issue in SAIGE when combining annotations with `DS` field input.

### Adapting for real data

1. Replace input VCFs and PLINK files with your own data
2. Update `--phenoCol` for each phenotype you wish to test — each phenotype needs a separate null model
3. Update `--covarColList` with your covariates
4. Prepare a group file mapping variants to genes with annotations for set-based testing
5. Adjust parameters:
    - `--minMAC`: Filter out very rare variants
    - `--maxMAF_in_groupTest`: MAF threshold for set-based tests

### References

- Zhou et al. (2018) *Nature Genetics*. [doi:10.1038/s41588-018-0184-y](https://doi.org/10.1038/s41588-018-0184-y)
- SAIGE GitHub: [github.com/saigegit/SAIGE](https://github.com/saigegit/SAIGE)
