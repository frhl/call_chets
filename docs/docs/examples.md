## Examples

The [`examples/`](https://github.com/frhl/call_chets/tree/main/examples) directory contains complete, runnable workflows with synthetic data. These demonstrate end-to-end pipelines from data generation through GWAS.

### Getting started

```bash
cd examples

# Generate synthetic data (phased/unphased VCFs, gene map, phenotypes)
./00_generate_data.sh

# Run one of the pipelines
./01_phased_gene_workflow.sh
```

### Generate synthetic data

`00_generate_data.sh` creates synthetic test data in `examples/input/`:

- `phased.vcf.gz` -- Phased genotypes (with `|` separator)
- `unphased.vcf.gz` -- Unphased genotypes (with `/` separator)
- `gene_map.txt` -- Variant-to-gene mapping
- `phenotypes.txt` -- Simulated phenotypes (`Y_rec`, `Y_add`, `Y_null`)
- `samples.txt` -- Sample ID list

---

### Workflow 1: Phased gene-level pipeline

`01_phased_gene_workflow.sh` -- Standard pipeline for phased data. Calls compound heterozygotes and aggregates variants by gene.

**Step 1: Extract phased genotypes**
```bash
bcftools query -i'GT="alt"' \
    -f'[%SAMPLE %CHROM:%POS:%REF:%ALT %GT\n]' \
    phased.vcf.gz | gzip > phased_sites.txt.gz
```

**Step 2: Call compound heterozygotes**
```bash
interpret_phase \
    --geno phased_sites.txt.gz \
    --gene-map gene_map.txt \
    --show-variants \
    > chet_results.txt
```

**Step 3: Create dominance pseudo-variant VCF**
```bash
make_pseudo_vcf \
    --input chet_results.txt \
    --samples samples.txt \
    --mode dominance \
    | bgzip > encoded_dominance.vcf.gz
```

**Step 4: Create additive pseudo-variant VCF**
```bash
make_pseudo_vcf \
    --input chet_results.txt \
    --samples samples.txt \
    --mode additive \
    | bgzip > encoded_additive.vcf.gz
```

**Step 5: Run joint GWAS**
```bash
Rscript run_gwas.R \
    --additive encoded_additive.vcf.gz \
    --dominance encoded_dominance.vcf.gz \
    --phenotype phenotypes.txt \
    --out gwas_results
```

The GWAS step fits: `Y ~ Additive + Dominance`

---

### Workflow 2: Unphased gene-level pipeline

`02_unphased_gene_workflow.sh` -- Pipeline for unphased data (burden/collapsing). Same structure as workflow 1, but uses the `--unphased` flag:

```bash
interpret_phase \
    --geno unphased_sites.txt.gz \
    --gene-map gene_map.txt \
    --unphased --show-variants \
    > unphased_results.txt
```

In unphased mode, `interpret_phase` performs a het/hom burden collapse without distinguishing compound heterozygotes from cis pairs.

---

### Workflow 3: Variant-level pipeline

`03_variant_workflow.sh` -- For analyzing individual variants directly using the `recode` tool, without gene-level collapsing.

**Step 1: Prepare additive VCF**

The input VCF already contains standard additive genotypes.

**Step 2: Create dominance VCF**
```bash
recode \
    --input phased.vcf.gz \
    --mode dominance \
    --scale-per-variant \
    | bgzip > encoded_dominance.vcf.gz
```

**Step 3: Run joint GWAS**
```bash
Rscript run_gwas.R \
    --additive phased.vcf.gz \
    --dominance encoded_dominance.vcf.gz \
    --phenotype phenotypes.txt \
    --out gwas_results
```

---

### GWAS analysis details

The included `run_gwas.R` script performs a joint regression: `Y ~ Additive + Dominance`.

- **Additive VCF**: Standard allele counts (0, 1, 2)
- **Dominance VCF**: Orthogonalized heterozygote deviation

Significant signals in both terms (especially Dominance) indicate non-additive/recessive effects.

#### Dose-response estimation

For each variant, the script calculates the estimated genetic effect for dosage levels 0, 1, and 2 by combining the additive and dominance contributions:

```
Effect(d) = Beta_Additive * d + Beta_Dominance * X_d
```

Where `X_d` is the dominance encoding value for dosage level `d` (empirically determined from the data).

Standard errors are calculated by combining the variances of both terms, assuming orthogonality. The output plots visualize these estimates with their 95% confidence intervals, helping to distinguish between additive (linear) and recessive (hockey-stick) architectures.

---

### Integration examples

For complete integration examples with popular GWAS tools, see:

- [**REGENIE** variant-level testing](regenie.md) -- Additive + dominance testing with REGENIE
- [**SAIGE** set-based testing](saige.md) -- Gene-based burden testing with SAIGE
