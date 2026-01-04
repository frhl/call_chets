# Simulation Pipeline

This directory contains scripts for generating simulated genetic data for SAIGE gene-based testing.

## Pipeline Overview

Run the scripts in order:

### 1. Generate Simulated Data
```bash
./01_simulate.sh
```

**What it does:**
- Simulates 5000 samples with 1000 variants using msprime
- Creates 20 genes (50 variants each)
- Implements gene-based causal architecture (10% causal genes by default)
- Generates VCF, phenotypes, and annotation files

**Key parameters** (edit in script):
- `CAUSAL_GENE_FRACTION=0.1` - Fraction of genes that are causal
- `CAUSAL_VARIANTS_PER_GENE=5` - Causal variants per causal gene
- `H2=0.1` - Heritability

**Output:**
- `../input/simulated.vcf` - Genotypes in VCF format
- `../input/simulated.phenos.tsv` - Phenotype data
- `../input/meta/gene_boundaries.tsv` - Gene definitions
- `../input/meta/variant_annotations.tsv` - Variant annotations (pLoF/synonymous)

---

### 2. Estimate Principal Components
```bash
./02_estimate_pcs.sh
```

**What it does:**
- Converts VCF to PLINK format
- LD prunes variants (r² < 0.2)
- Computes 10 PCs using PLINK2

**Output:**
- `../input/pca/pcs.txt` - PCs formatted for SAIGE
- `../input/pca/pca_results.eigenvec` - Full PLINK output
- `../input/pca/ld_pruned.prune.in` - Pruned variant list

---

### 3. Prepare SAIGE Inputs
```bash
./03_prepare_saige_inputs.sh
```

**What it does:**
- Creates SAIGE group file (variant-gene mapping)
- Merges phenotypes with covariates and PCs
- Creates gene set files for burden testing

**Key parameters** (edit in script):
- `ANNOTATION_FILTER="pLoF"` - Which variants to test (pLoF, synonymous, or all)

**Output:**
- `../input/meta/variant_gene_mapping.txt` - Variant to gene mapping
- `../input/meta/simulated.phenos.with_covariates.tsv` - Phenotypes with covariates and PCs
- `../input/meta/genesets_pLoF.txt` - Gene variant sets for testing

---

## Complete Workflow Example

```bash
# 1. Generate data
cd simulation/
./01_simulate.sh

# 2. Estimate PCs
./02_estimate_pcs.sh

# 3. Prepare SAIGE inputs
./03_prepare_saige_inputs.sh

# 4. Run SAIGE (from parent directory)
cd ..
./03_saige_step1.sh  # Fit null model
./04_saige_step2.sh  # Run gene-based tests
```

---

## Customizing the Simulation

### High Power Scenario (easier to detect)
Edit `01_simulate.sh`:
```bash
CAUSAL_GENE_FRACTION=0.2    # More causal genes
CAUSAL_VARIANTS_PER_GENE=10 # More causal variants
H2=0.2                       # Higher heritability
```

### Low Power Scenario (harder to detect)
Edit `01_simulate.sh`:
```bash
CAUSAL_GENE_FRACTION=0.05   # Fewer causal genes
CAUSAL_VARIANTS_PER_GENE=2  # Fewer causal variants
H2=0.05                      # Lower heritability
```

---

## Docker Requirements

The pipeline uses Docker containers:
- **plink2**: `biocontainer/plink2:alpha2.3_jan2020`
- **SAIGE**: `wzhou88/saige:0.5.1` (used in parent directory scripts)

---

## Notes

- All output goes to `../input/` directory
- Simulation uses conda environment `pop_sim` (requires msprime, numpy, pandas, scipy)
- The gene-based causal model is more realistic than random variant selection
- PCs help control for population structure even in simulated data
