## Documentation

### `interpret_phase`

Identifies compound heterozygous (chet) and homozygous (hom) variant configurations within gene regions from phased or unphased genotype data.

#### Usage

```bash
interpret_phase --geno <file> --gene-map <file> [options]
```

#### Required options

| Option | Description |
|--------|-------------|
| `--geno` / `-g` `<file>` | Phased or unphased genotype file (gzipped). Produced by `bcftools query` |
| `--gene-map` / `-m` `<file>` | Tab-separated variant-to-gene mapping file |

#### Optional options

| Option | Description |
|--------|-------------|
| `--unphased` | Run in unphased mode (het/hom burden only, no compound het calling) |
| `--info-map` / `-i` `<file>` | Variant info file (AF, AC, etc.) to annotate output |
| `--score-map` / `-p` `<file>` | Variant score file for weighted collapsing |
| `--show-variants` / `-sv` | Include detailed variant info in output |
| `--verbose` / `-v` | Enable verbose logging to stderr |

#### Score/collapse options

These options require `--score-map`:

| Option | Description |
|--------|-------------|
| `--haplotype-collapse` / `-hc` | Haplotype score collapse rule: `product` (default), `min`, `max`, `additive` |
| `--gene-collapse` / `-gc` | Gene score collapse rule: `product` (default), `min`, `max`, `additive` |

#### Input file formats

**Genotype file** (`--geno`): Tab/space-separated, one line per sample-variant pair. Produced by:

```bash
bcftools query -i'GT="alt"' -f'[%SAMPLE %CHROM:%POS:%REF:%ALT %GT\n]' input.vcf.gz | gzip > genotypes.txt.gz
```

Format: `SAMPLE_ID CHROM:POS:REF:ALT GENOTYPE`

- Phased genotypes use `|` separator (e.g., `0|1`, `1|0`)
- Unphased genotypes use `/` separator (e.g., `0/1`)

**Gene map file** (`--gene-map`): Tab-separated, no header.

```
CHROM:POS:REF:ALT    GENE_NAME
1:12345:A:G          BRCA1
1:12456:C:T          BRCA1
```

**Score map file** (`--score-map`): Tab-separated, no header.

```
CHROM:POS:REF:ALT    SCORE
1:12345:A:G          0.95
1:12456:C:T          0.80
```

---

### `make_pseudo_vcf`

Converts `interpret_phase` output into a biallelic VCF with dosage encodings. Each gene becomes a pseudo-variant in the output VCF.

#### Usage

```bash
make_pseudo_vcf --input <file> --samples <file> --mode <mode> [options]
```

#### Required options

| Option | Description |
|--------|-------------|
| `--input` / `-i` `<file>` | Output file from `interpret_phase` |
| `--samples` / `-s` `<file>` | Sample list (one sample ID per line, no header) |
| `--mode` / `-m` `<mode>` | Encoding mode (see below) |

#### Encoding modes

| Mode | Aliases | Description |
|------|---------|-------------|
| `additive` | `012` | Standard 0, 1, 2 dosages |
| `dominance` | `nonadditive` | Non-additive: orthogonalized heterozygote deviation |
| `recessive` | `001` | Recessive encoding (0 and 2 only) |

#### Optional options

| Option | Description |
|--------|-------------|
| `--min-ac` `<n>` | Minimum allele count filter (sum of DS >= n) |
| `--max-ac` `<n>` | Maximum allele count filter (sum of DS < n) |
| `--all-info` | Include detailed INFO fields (variant lists, counts) |
| `--global-dom-dosage` | Use global min/max dominance dosage for scaling |
| `--no-dosage-scaling` | Disable scaling of dosages to [0, 2] (dominance mode only) |
| `--scaling-factor` `<f>` | Apply additional scaling factor (default: 1.0) |
| `--force-chr-out-name` `<name>` | Force output chromosome name |
| `--suffix` `<str>` | Append suffix to output gene names |

---

### `recode`

Orthogonalizes or recodes an existing VCF for non-additive analysis. Works directly on VCF/BCF files without needing the gene-level pipeline.

#### Usage

```bash
recode --input <file.vcf.gz> [options]
```

#### Required options

| Option | Description |
|--------|-------------|
| `--input` / `-i` `<file>` | Input VCF/BCF (`.vcf`, `.vcf.gz`, `.bcf`) |

#### Mode

| Option | Description |
|--------|-------------|
| `--mode` / `-m` `<mode>` | `nonadditive` (default, alias: `dominance`) or `recessive` |

#### Scaling options

Mutually exclusive -- choose one:

| Option | Description |
|--------|-------------|
| `--scale-per-variant` | Scale each variant independently to [0, 2]. Use for variant-level tests |
| `--scale-globally` | Scale using global min/max across all variants. Use for set-based tests |
| `--scale-by-group` `<file>` | Scale within groups defined by a tab-separated file (variant, gene). Produces comparable betas within genes |
| `--scale-factor` `<f>` | Apply additional scaling factor (default: 1.0). Can be combined with any scaling option above |

#### Filter options

| Option | Description | Default |
|--------|-------------|---------|
| `--min-hom-count` `<n>` | Minimum minor homozygous count to include variant | 1 |
| `--min-het-count` `<n>` | Minimum heterozygous count to include variant | 1 |
| `--min-aac` / `--max-aac` `<n>` | Min/max alternate allele count | -- |
| `--min-mac` / `--max-mac` `<n>` | Min/max minor allele count | -- |
| `--min-aaf` / `--max-aaf` `<f>` | Min/max alternate allele frequency | -- |
| `--min-maf` / `--max-maf` `<f>` | Min/max minor allele frequency | -- |

#### Other options

| Option | Description |
|--------|-------------|
| `--set-variant-id` | Set variant IDs to `chr:pos:ref:alt` format |
| `--all-info` | Include frequency/scaling info in INFO fields |

