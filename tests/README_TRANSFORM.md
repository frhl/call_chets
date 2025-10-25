# Transform Unit Tests

## Overview

This test suite validates the `transform` binary which converts VCF genotype/dosage data into dominance or recessive encoding formats.

## Test Files

### Input VCF Files

- **test_transform_simple.vcf**: Basic VCF with GT fields, 4 samples, 3 variants
  - Variant 1 (chr1:1000): Genotypes `0/0, 0/1, 1/1, 0/1`
  - Variant 2 (chr1:2000): Genotypes `0/1, 1/1, 0/0, 1/1`
  - Variant 3 (chr1:3000): Genotypes `1/1, 0/1, 0/1, 1/1`

- **test_transform_dosage.vcf**: VCF with DS (dosage) fields instead of GT
  - Same variants as simple.vcf but with dosage values (0.0, 1.0, 2.0)

- **test_transform_no_homalt.vcf**: Edge case VCF
  - Variant 1: No homozygous alternate alleles (should be skipped in dominance mode)
  - Variant 2: Has homozygous alternates (should be processed)

### Gene Mapping File

- **test_transform_gene_map.txt**: Maps variants to genes for gene-specific scaling
  ```
  variant           gene
  chr1:1000:A:T    GENE1
  chr1:2000:C:G    GENE1
  chr1:3000:G:A    GENE2
  ```

## Running Tests

```bash
cd /Users/flassen/Projects/06_call_chets/call_chets/tests
./run_transform_tests.sh
```

## Test Coverage

### Basic Dominance Mode Tests (7 tests)
1. **dominance_mode_basic_gt**: Basic dominance encoding with GT field
2. **dominance_mode_dosage**: Dominance encoding with DS field
3. **dominance_mode_scaled**: With dosage scaling to [0,2] range
4. **dominance_mode_variant_id**: With variant ID formatting
5. **dominance_mode_all_info**: With additional frequency information
6. **dominance_mode_gene_map**: With gene-specific scaling
7. **dominance_mode_full**: All options combined

### Recessive Mode Tests (3 tests)
8. **recessive_mode_basic**: Basic recessive encoding (het=0, hom=2)
9. **recessive_mode_dosage**: Recessive with DS field
10. **recessive_mode_variant_id**: Recessive with variant ID formatting

### Edge Case Tests (2 tests)
11. **no_homozygous_alternate**: Validates skipping of variants without homozygous alternates
12. **scaling_factor**: Tests custom scaling factor application

### Error Handling Tests (3 tests)
13. **invalid_mode**: Validates rejection of invalid modes
14. **missing_input**: Validates error handling for missing input file
15. **missing_gene_map**: Validates error handling for missing gene map

## Dominance Encoding Formula

For a variant with genotype frequencies:
- `r` = frequency of ref/ref (0/0)
- `h` = frequency of het (0/1)
- `a` = frequency of alt/alt (1/1)

Dosage values:
- **aa (0/0)**: `-h * a`
- **Aa (0/1)**: `2 * a * r`
- **AA (1/1)**: `-h * r`

### Example Calculation

For chr1:1000 with genotypes `0/0, 0/1, 1/1, 0/1`:
- Counts: aa=1, Aa=2, AA=1 (n=4)
- Frequencies: r=0.25, h=0.5, a=0.25
- Expected dosages:
  - aa: `-0.5 * 0.25 = -0.125`
  - Aa: `2 * 0.25 * 0.25 = 0.125`
  - AA: `-0.5 * 0.25 = -0.125`

## Key Validation Features

The updated transform.cpp now validates inputs **before** producing any output:

1. ✅ Input file must exist and be readable
2. ✅ VCF must have valid header
3. ✅ VCF must have at least one sample
4. ✅ VCF must have GT or DS format fields
5. ✅ VCF must contain at least one variant
6. ✅ Gene map file (if specified) must exist and contain data
7. ✅ For dominance mode: at least one variant must have homozygous alternate alleles

**This ensures no output is generated if validation fails.**

## Expected Behavior

- **Dominance mode**: Requires variants with AA (1/1) genotypes; skips others
- **Recessive mode**: Sets heterozygotes to 0, keeps homozygotes as 2
- **Scaling**: When enabled, transforms dosages to [0,2] range
- **Gene map**: Enables per-gene scaling instead of global scaling
- **Error handling**: Returns non-zero exit code and prints clear error messages
- **Missing values**: Missing genotypes (`./.`) or hemizygous genotypes are excluded from frequency calculations
  - Frequencies (r, h, a) are calculated using only non-missing samples as denominator
  - Missing samples output dosage of 0
  - A note is printed if any missing values are encountered
  - Example: For 4 samples with genotypes `0/0, 0/1, ./., 1/1`:
    - Only 3 samples used for frequency calculation (n=3)
    - r = 1/3, h = 1/3, a = 1/3 (not r = 1/4, h = 1/4, a = 1/4)
- **Low variance warning**: When scaling is enabled, warns if any two genotype classes have very similar scaled dosages
  - Checks pairwise differences: |aa - Aa|, |aa - AA|, |Aa - AA|
  - Warns if minimum difference < 0.0001 (threshold)
  - Example: variant with scaled dosages aa=1.0000, Aa=1.00005, AA=0
    - |aa - Aa| = 0.00005 < 0.0001 → Warning triggered
  - Such variants have limited power in association testing due to lack of contrast between genotypes

## Test Output Format

Tests show:
- ✓ Green checkmarks for passed tests
- ✗ Red X for failed tests
- Blue sample output showing first variant
- Summary statistics at the end
