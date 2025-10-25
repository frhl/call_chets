# call_chets --unphased Mode: Behavior Reference

## Quick Comparison

| Feature | Phased Mode (default) | Unphased Mode (`--unphased`) |
|---------|----------------------|------------------------------|
| **Input formats accepted** | `1\|0`, `0\|1`, `1\|1` | `0/1`, `1/0`, `0\|1`, `1\|0`, `1/1`, `1\|1` |
| **Phase information** | Used | **Ignored** |
| **Output call types** | `het`, `cis`, `chet`, `hom` | `het`, `hom` only |
| **Compound het detection** | Yes (different haplotypes) | No |
| **Cis variant detection** | Yes (same haplotype) | No |
| **Scoring/pathogenicity** | Supported (with `--score-map`) | **Not supported** (ignored with warning) |

## Genotype Interpretation Examples

### Input: Two heterozygous variants in same gene

#### Phased Mode
```
SAMPLE1  chr1:1000:A:T  1|0
SAMPLE1  chr1:2000:C:G  0|1
```
**Output:** `SAMPLE1  chr1  GENE1  chet  2`
- Variants on **different haplotypes** → compound heterozygote

```
SAMPLE1  chr1:1000:A:T  1|0
SAMPLE1  chr1:2000:C:G  1|0
```
**Output:** `SAMPLE1  chr1  GENE1  cis  1`
- Variants on **same haplotype** → cis configuration

#### Unphased Mode
```
SAMPLE1  chr1:1000:A:T  0/1  (or 0|1 or 1|0)
SAMPLE1  chr1:2000:C:G  0/1  (or 0|1 or 1|0)
```
**Output:** `SAMPLE1  chr1  GENE1  het  1`
- Phase **ignored** → simple heterozygote

### Input: Homozygous variant

#### Both Modes (identical behavior)
```
SAMPLE2  chr1:1000:A:T  1/1  (or 1|1)
```
**Output:** `SAMPLE2  chr1  GENE1  hom  2`

## Normalization Rules (Unphased Mode)

When `--unphased` is active:

| Input Genotype | Normalized To | Call Type | Dosage |
|----------------|---------------|-----------|--------|
| `0/1` | `0/1` | het | 1 |
| `1/0` | `0/1` | het | 1 |
| `0\|1` | `0/1` | het | 1 |
| `1\|0` | `0/1` | het | 1 |
| `1/1` | `1/1` | hom | 2 |
| `1\|1` | `1/1` | hom | 2 |

**Key point:** Phase separator (`|` vs `/`) and allele order are both ignored in unphased mode.

## Use Cases

### When to use Phased Mode (default)
- You have phased genotype data (from statistical phasing or trio data)
- You want to detect compound heterozygotes
- You're interested in cis vs trans configurations
- You need haplotype-specific scoring

### When to use Unphased Mode (`--unphased`)
- You have unphased genotype data (standard VCF without phasing)
- You only care about variant presence (het vs hom)
- You want simpler gene-level burden testing
- Your data mixes phased and unphased variants

## Command Examples

### Standard phased analysis
```bash
./call_chets --geno phased_genotypes.txt.gz \
             --gene-map gene_mapping.txt.gz \
             --score-map pathogenicity_scores.txt.gz
```

### Unphased analysis (ignoring phase)
```bash
./call_chets --geno genotypes.txt.gz \
             --gene-map gene_mapping.txt.gz \
             --unphased
```

### Unphased with variant display
```bash
./call_chets --geno genotypes.txt.gz \
             --gene-map gene_mapping.txt.gz \
             --info-map variant_annotations.txt.gz \
             --unphased \
             --show-variants
```

## Warnings and Restrictions

### Unphased Mode Restrictions
These arguments are **ignored** in `--unphased` mode (with warnings):
- `--score-map` / `-p`
- `--haplotype-collapse` / `-hc`
- `--gene-collapse` / `-gc`
- `--show-haplotype-scores` / `-shs`

### Phased Mode Restrictions
- Unphased genotypes (`0/1`) are **rejected** with warnings
- Mixed data will cause errors

## Real-World Example

### Phased Data
```
# Input (phased)
SAMPLE1  chr1:1000:A:T  1|0
SAMPLE1  chr1:2000:C:G  0|1

# Phased mode output
SAMPLE1  chr1  GENE1  chet  2

# Unphased mode output (--unphased)
SAMPLE1  chr1  GENE1  het  1
```

### Impact on Association Testing
- **Phased mode:** Tests for recessive effects (2 copies needed)
- **Unphased mode:** Tests for any variant presence (simpler burden test)

Choose based on your genetic hypothesis!
