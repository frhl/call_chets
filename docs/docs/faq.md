## Frequently asked questions

### General

*When should I use the gene-level pipeline vs. the variant-level pipeline?*

Use the **gene-level pipeline** (`interpret_phase` + `make_pseudo_vcf`) when you want to:

- Identify compound heterozygotes within genes
- Aggregate rare variants by gene for burden-style tests
- Test gene-level non-additive effects

Use the **variant-level pipeline** (`recode`) when you want to:

- Test individual variants for non-additive effects
- Orthogonalize existing VCFs for dominance deviation analysis
- Keep variant-level resolution (no gene collapsing)

---

*What is the difference between `dominance` and `recessive` mode?*

- **Dominance** mode produces an orthogonalized encoding that captures the deviation of heterozygotes from the additive expectation. The encoding is continuous and orthogonal to the additive component, allowing joint testing (`Y ~ Additive + Dominance`).
- **Recessive** mode produces a binary encoding where only homozygous alternate genotypes have a non-zero value (het=0, hom_alt=2). This directly tests for recessive effects.

---

*Can I use arcade with unphased data?*

Yes. Use `interpret_phase --unphased` for the gene-level pipeline. In unphased mode, the tool performs a het/hom burden collapse without distinguishing compound heterozygotes from cis pairs. The variant-level pipeline (`recode`) works with both phased and unphased VCFs.

---

*What is the difference between `--scale-per-variant`, `--scale-globally`, and `--scale-by-group`?*

These options control how the dominance dosages are scaled in `recode`:

- `--scale-per-variant`: Each variant is scaled independently to [0, 2]. Use this for **variant-level** tests where each variant is tested separately.
- `--scale-globally`: All variants are scaled using the global min/max across the dataset. Use this for **set-based** tests (e.g., SAIGE gene-based) where dosages must be comparable across variants.
- `--scale-by-group <file>`: Variants are scaled within groups defined by a mapping file. Use this when you want comparable betas **within genes** but not necessarily across genes.

---

### Input / Output

*What VCF fields does arcade use?*

- `recode` reads the `GT` (genotype) field from the input VCF and writes `DS` (dosage) to the output VCF.
- `make_pseudo_vcf` writes both `GT` and `DS` fields.
- Downstream tools (REGENIE, SAIGE) should read the `DS` field for non-additive encodings.

---

*Why does the dominance VCF have fewer variants than the additive VCF?*

The `recode` tool applies filters (`--min-hom-count`, `--min-het-count`, `--max-maf`) that remove variants with insufficient genotype counts. Dominance encoding requires heterozygous individuals, so very rare variants with no heterozygotes are excluded.

---

*What format does `interpret_phase` expect for input genotypes?*

The genotype file should be produced by `bcftools query`:

```bash
bcftools query -i'GT="alt"' \
    -f'[%SAMPLE %CHROM:%POS:%REF:%ALT %GT\n]' \
    input.vcf.gz | gzip > genotypes.txt.gz
```

This produces space-separated lines: `SAMPLE_ID CHROM:POS:REF:ALT GENOTYPE`

---

### GWAS integration

*How do I run a joint additive + dominance test?*

Produce both additive and dominance VCFs, then include both as predictors in your association model. For example:

- With **REGENIE**: Run Step 2 separately on each encoding, then combine results
- With **SAIGE**: Run `step2_SPAtests.R` separately with `--vcfField=GT` (additive) and `--vcfField=DS` (dominance)
- With **R**: `lm(Y ~ Additive + Dominance)`

The orthogonalized dominance encoding ensures that the additive and dominance components are independent.

---

*Should I use `--scale-per-variant` or `--scale-globally` for REGENIE?*

For **variant-level** testing with REGENIE, use `--scale-per-variant`. This scales each variant's dominance dosage independently, which is appropriate when each variant is tested separately.

For **set-based** analysis, use `--scale-globally` so that dosages are comparable across variants within a set.

---

### Docker

*How do I run arcade via Docker?*

```bash
docker pull fhlassen/arcade:latest

# Run any tool
docker run -v $PWD:/data fhlassen/arcade \
    interpret_phase --help

# Process files in current directory
docker run -v $PWD:/data fhlassen/arcade \
    recode \
        --input /data/variants.vcf.gz \
        --mode dominance \
        --scale-per-variant \
    | bgzip > dominance.vcf.gz
```

---

*I get "bcftools not found" when running the examples.*

The gene-level pipeline requires [bcftools](https://samtools.github.io/bcftools/) for the genotype extraction step. Install it via your package manager:

```bash
# macOS
brew install bcftools

# Ubuntu/Debian
sudo apt-get install bcftools
```
