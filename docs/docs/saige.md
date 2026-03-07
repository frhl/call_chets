## SAIGE integration

This page demonstrates set-based (gene-level) genetic association testing using [SAIGE](https://github.com/saigegit/SAIGE) with genotype encodings produced by **arcade**.

The full example is in [`examples/saige-set-based/`](https://github.com/frhl/call_chets/tree/main/examples/saige-set-based).

### Overview

The SAIGE pipeline performs:

1. **Simulation** -- Generate realistic genotype data with `msprime`, assign gene structures and variant annotations, simulate phenotypes
2. **Encoding** -- Use `recode` to produce dominance and recessive VCFs
3. **Null model** -- Fit SAIGE null model with sparse GRM
4. **Association testing** -- Run variant-level and gene-based (set-based) burden tests

This example tests **three genetic encodings**:

- **Additive**: Standard 0/1/2 genotype (GT field)
- **Dominance (non-additive)**: Orthogonalized heterozygote deviation (DS field from `recode`)
- **Recessive**: Binary recessive encoding, 0/0/2 (DS field from `recode`)

### Quick start

```bash
cd examples/saige-set-based

# Generate simulated data
cd simulation/
./01_simulate.sh
./02_estimate_pcs.sh
./03_prepare_saige_inputs.sh
cd ..

# Encode VCFs
./01_encode_vcf.sh

# Run SAIGE
./03_saige_step0.sh     # Sparse GRM
./04_saige_step1.sh     # Null model
./05_saige_step2_variant.sh   # Variant-level tests
./06_saige_step2_group.sh     # Gene-based tests
```

### Requirements

- **Docker** (SAIGE and PLINK2 run in containers)
- Docker images:
    - `wzhou88/saige:0.5.1`
    - `biocontainer/plink2:alpha2.3_jan2020`
- **Conda** with `msprime` (for simulation only)

### Step 1: Encode VCFs with `recode`

The encoding step uses `recode` to create dominance and recessive VCFs from the additive input:

```bash
# Dominance (non-additive) encoding
recode \
    --input simulated.vcf.gz \
    --mode nonadditive \
    --min-hom-count 5 \
    --max-maf 0.05 \
    --scale-globally \
    --set-variant-id \
    --all-info \
    | bgzip > simulated.nonadditive.vcf.gz

# Recessive encoding
recode \
    --input simulated.vcf.gz \
    --mode recessive \
    --min-hom-count 5 \
    --set-variant-id \
    | bgzip > simulated.recessive.vcf.gz
```

Note the use of `--scale-globally` for set-based analysis, which ensures comparable dosages across variants within a gene. For variant-level analysis, use `--scale-per-variant` instead.

### Step 2: Sparse GRM

```bash
./03_saige_step0.sh
```

Creates a sparse genetic relationship matrix using 1,000 random markers:

```bash
createSparseGRM.R \
    --plinkFile=simulated \
    --nThreads=4 \
    --outputPrefix=sparseGRM \
    --numRandomMarkerforSparseKin=1000 \
    --relatednessCutoff=0.125
```

### Step 3: Fit null model

```bash
./04_saige_step1.sh
```

```bash
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

The variance ratio markers are pre-selected from two MAC categories (10-20 and >=20) by `02_prepare_vr.sh`.

### Step 4: Variant-level association testing

```bash
./05_saige_step2_variant.sh
```

Tests three encodings against the same null model:

**Additive** (standard GT field):
```bash
step2_SPAtests.R \
    --vcfFile=simulated.vcf.gz \
    --vcfField=GT \
    --minMAC=2 \
    --SAIGEOutputFile=saige.step2.additive.variant.txt
```

**Recessive** (DS field, 0/0/2 encoding):
```bash
step2_SPAtests.R \
    --vcfFile=simulated.recessive.vcf.gz \
    --vcfField=DS \
    --minMAC=2 \
    --SAIGEOutputFile=saige.step2.recessive.variant.txt
```

**Non-additive** (DS field, orthogonalized dominance):
```bash
step2_SPAtests.R \
    --vcfFile=simulated.nonadditive.vcf.gz \
    --vcfField=DS \
    --minMAC=2 \
    --SAIGEOutputFile=saige.step2.nonadditive.variant.txt
```

### Step 5: Gene-based (set-based) testing

```bash
./06_saige_step2_group.sh
```

Tests gene-level burden using variant annotation groups (pLoF, synonymous):

**Additive burden** (GT field):
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

**Non-additive burden** (DS field, tested per annotation):
```bash
step2_SPAtests.R \
    --vcfFile=simulated.nonadditive.vcf.gz \
    --vcfField=DS \
    --groupFile=genesets_all.txt \
    --annotation_in_groupTest=pLoF \
    --r.corr=1 \
    --SAIGEOutputFile=saige.step2.nonadditive.group.pLoF.txt
```

### Simulation details

The simulation pipeline uses `msprime` to generate realistic population genetic data:

| Parameter | Default |
|-----------|---------|
| Samples | 5,000 diploid individuals |
| Variants | 1,000 |
| Variants per gene | 50 (= 20 genes) |
| MAF range | 1-5% (70% of variants) |
| Causal gene fraction | 10% |
| Causal variants per gene | 5 |
| Heritability (h^2) | 0.1 |
| Architecture | Additive |

Supported genetic architectures: `additive`, `recessive`, `dominant`.

To customize:
```bash
cd simulation/
conda run -n pop_sim python 01_simulate.py \
    --n_samples 10000 \
    --h2 0.2 \
    --architecture recessive \
    --seed 123
```

### Comparison with REGENIE

This example uses the same simulated data as the [REGENIE example](regenie.md), allowing direct comparison:

| Feature | SAIGE (set-based) | REGENIE (variant-level) |
|---------|-------------------|-------------------------|
| Test unit | Gene/variant sets | Individual variants |
| Step 1 | Fit null mixed model | Whole-genome regression |
| Step 2 | Gene burden tests | Variant association tests |
| Output | Gene-level p-values | Variant-level p-values |

### References

- Zhou et al. (2018) *Nature Genetics*. [doi:10.1038/s41588-018-0184-y](https://doi.org/10.1038/s41588-018-0184-y)
- SAIGE GitHub: [github.com/saigegit/SAIGE](https://github.com/saigegit/SAIGE)
