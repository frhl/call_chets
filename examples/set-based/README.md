# Set-Based Genetic Simulation Pipeline

This directory contains a complete pipeline for simulating genotype data, gene structures, variant annotations, and phenotypes for power analysis and testing of set-based genetic association methods.

## Overview

The simulation pipeline performs the following steps:

1. **Genotype Simulation**: Uses `msprime` to simulate realistic population genetic data
2. **Gene Structure**: Organizes variants into contiguous "genes" (sets of 50 variants each)
3. **Variant Annotation**: Assigns variants as either "pLoF" (causal) or "synonymous" (non-causal)
4. **Phenotype Simulation**: Generates quantitative phenotypes based on causal variants

## Quick Start

```bash
# Run the complete simulation pipeline
cd simulation/
./01_simulate.sh
./02_estimate_pcs.sh
./03_prepare_saige_inputs.sh

# Then run SAIGE
cd ..
./03_saige_step1.sh
./04_saige_step2.sh
```

## Directory Structure

### `simulation/` - Data generation pipeline

All simulation scripts have been moved to the `simulation/` directory:

- `01_simulate.sh` - Generate simulated genotypes and phenotypes
- `01_simulate.py` - Python simulation engine (called by bash script)
- `02_estimate_pcs.sh` - Estimate principal components using PLINK2
- `03_prepare_saige_inputs.sh` - Prepare SAIGE input files
- `03_prepare_saige_inputs.py` - Python helper for SAIGE inputs

See `simulation/README.md` for detailed documentation.

### Main directory - SAIGE analysis scripts

- `03_saige_step1.sh` - Fit SAIGE null model
- `04_saige_step2.sh` - Run gene-based burden tests
- `05_plot_results.sh` - Visualize results

### Output Files

After running the simulation, you will find:

**In `input/` directory:**
- `simulated.genotypes.tsv.gz` - Genotype matrix in TSV format (samples × variants)
- `simulated.vcf` - Genotype matrix in VCF format (standard variant call format)
- `simulated.variants.tsv` - Variant information (position, MAF, allele counts)
- `simulated.phenotypes.tsv` - Simulated phenotypes for all samples

**In `input/helpers/` directory:**
- `gene_boundaries.tsv` - Gene definitions (20 genes, 50 variants each)
- `variant_annotations.tsv` - Variant annotations (pLoF vs synonymous)

## Simulation Parameters

Default parameters in `simulation/01_simulate.sh`:

- **N_SAMPLES**: 5000 diploid individuals
- **TARGET_VARIANTS**: 1000 variants
- **VARIANTS_PER_GENE**: 50 variants per gene (= 20 genes total)
- **MAF_RANGE**: 1-5% (70% of variants in this range)
- **CAUSAL_GENE_FRACTION**: 0.1 (10% of genes harbor causal variants)
- **CAUSAL_VARIANTS_PER_GENE**: 5 (average causal variants per causal gene)
- **H2**: 0.1 (heritability)
- **ARCHITECTURE**: additive genetic model
- **SEED**: 42 (for reproducibility)

## Customization

You can modify parameters by editing `simulation/01_simulate.sh` or by calling the Python script directly:

```bash
cd simulation/
conda run -n pop_sim python 01_simulate.py \
    --n_samples 10000 \
    --target_variants 2000 \
    --variants_per_gene 100 \
    --causal_gene_fraction 0.2 \
    --causal_variants_per_gene 10 \
    --h2 0.2 \
    --architecture recessive \
    --seed 123 \
    --output_dir ../input
```

### Available Parameters

```
--n_samples N                    Number of diploid samples (default: 5000)
--target_variants N              Number of variants to simulate (default: 1000)
--variants_per_gene N            Variants per gene (default: 50)
--maf_min FLOAT                  Minimum MAF for enrichment (default: 0.01)
--maf_max FLOAT                  Maximum MAF for enrichment (default: 0.05)
--maf_fraction FLOAT             Fraction in MAF range (default: 0.7)
--causal_gene_fraction FLOAT     Fraction of genes with causal variants (default: 0.1)
--causal_variants_per_gene INT   Causal variants per causal gene (default: 5)
--h2 FLOAT                       Heritability (default: 0.1)
--architecture STR               Genetic model: additive, recessive, dominant
--seed INT                       Random seed (default: 42)
--output_dir PATH                Output directory (default: ../input)
--output_prefix STR              Output file prefix (default: simulated)
--chromosome STR                 Chromosome name for VCF (default: 1)
```

## Requirements

- **Conda environment**: `pop_sim` (must have msprime, tskit, pandas, numpy)
- Python 3.8+

To create the conda environment:

```bash
conda create -n pop_sim python=3.10
conda activate pop_sim
conda install -c conda-forge msprime tskit pandas numpy
```

## Simulation Details

### Genotype Simulation

Uses `msprime` coalescent simulator with:
- Effective population size: 10,000
- Sequence length: 1 Mb
- Mutation rate: 2×10⁻⁸ per bp per generation
- Recombination rate: 1×10⁻⁸ per bp per generation

### Gene Structure

- Variants are organized into contiguous chunks
- Default: 1000 variants → 20 genes of 50 variants each
- Gene boundaries are based on physical position

### Variant Annotations (Gene-Based Causal Model)

The simulation uses a realistic gene-based causal architecture:

1. **Select causal genes**: A fraction of genes (default: 10%) are designated as causal
2. **Select causal variants**: Within each causal gene, a specific number of variants (default: 5) are designated as causal
3. **Assign annotations**:
   - **pLoF variants**: Causal variants within causal genes
   - **Synonymous variants**: All non-causal variants

This approach is more realistic than random variant selection because:
- Diseases typically result from disruption of specific genes
- Concentrates signal in specific genes for better power
- Makes interpretation easier (know which genes should be significant)

### Phenotype Simulation

Phenotypes are simulated as:

```
Y = G + E

where:
G = genetic component (from causal variants)
E = environmental component (random noise)

var(G) / [var(G) + var(E)] = h²
```

Supported genetic architectures:
- **Additive**: Effect proportional to allele count [0, 1, 2]
- **Recessive**: Only homozygous alt (2) has effect
- **Dominant**: Heterozygous (1) and homozygous alt (2) have same effect

## Output File Formats

### genotypes.tsv.gz
- Rows: samples (sample_0, sample_1, ...)
- Columns: variants (var_0, var_1, ...)
- Values: genotypes [0, 1, 2]

### simulated.vcf
Standard VCF v4.2 format with:
- INFO fields: AF (allele frequency), MAF (minor allele frequency)
- FORMAT: GT (genotype)
- Genotypes encoded as: 0/0 (hom ref), 0/1 (het), 1/1 (hom alt)
- Chromosome name configurable via `--chromosome` parameter

### variants.tsv
```
variant_id  position  n_hom_ref  n_het  n_hom_alt  af      maf
var_0       1234      4850       140    10         0.016   0.016
var_1       5678      4900       95     5          0.0105  0.0105
```

### phenotypes.tsv
```
sample_id    phenotype
sample_0     1.234
sample_1     -0.567
```

### gene_boundaries.tsv
```
gene_id  gene_index  start_variant_idx  end_variant_idx  n_variants  start_position  end_position  variants
gene_0   0           0                  49               50          1234            98765         var_0,var_1,...
gene_1   1           50                 99               50          100234          198765        var_50,var_51,...
```

### variant_annotations.tsv
```
variant_id  position  maf     gene_id  annotation   is_causal
var_0       1234      0.016   gene_0   pLoF         True
var_1       5678      0.0105  gene_0   synonymous   False
```

## Reference

This simulation pipeline was inspired by the approach used in:
`/Users/flassen/Projects/11_wes_ko_ukbb_nexus/wes_ko_ukbb_nexus/scripts/simulation/saige/power/01_make_power_phenos.R`

The phenotype simulation logic follows similar principles for controlling heritability and genetic architecture.
