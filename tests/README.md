# call_chets Unit Tests

This directory contains unit tests for the `call_chets` program, specifically testing the `--unphased` mode functionality.

## Test Files

### Input Files
- **gene_map.txt.gz**: Gene mapping file with 5 variants mapping to 3 genes
- **genotypes_unphased.txt.gz**: Test genotypes in unphased format (0/1, 1/1)
- **genotypes_phased.txt.gz**: Test genotypes in phased format (0|1, 1|0, 1|1)
- **genotypes_mixed.txt.gz**: Test genotypes with mixed phased and unphased formats

### Expected Output Files
- **expected_unphased_mode_unphased_data.txt**: Expected output for unphased mode with unphased input
- **expected_unphased_mode_phased_data.txt**: Expected output for unphased mode with phased input
- **expected_phased_mode_phased_data.txt**: Expected output for phased mode (default) with phased input

## Test Scenarios

### Test 1: Unphased mode with unphased data
Tests that `--unphased` correctly processes unphased genotypes (0/1, 1/1)

**Expected behavior:**
- Heterozygous variants (0/1, 1/0) → `het`, dosage=1
- Homozygous variants (1/1) → `hom`, dosage=2

### Test 2: Unphased mode with phased data
Tests that `--unphased` accepts phased genotypes but treats them as unphased

**Expected behavior:**
- Phased hets (0|1, 1|0) → `het`, dosage=1
- Homozygous (1|1) → `hom`, dosage=2
- Phase information is **ignored**

### Test 3: Unphased mode with mixed data
Tests that `--unphased` accepts both phased and unphased formats in the same file

**Expected behavior:**
- All formats accepted
- Phase information ignored
- Same output as Test 2

### Test 4: Phased mode with phased data (default behavior)
Tests normal phased mode processing

**Expected behavior:**
- Compound heterozygotes (variants on different haplotypes) → `chet`, dosage=2
- Homozygous variants → `hom`, dosage=2
- Cis variants (multiple on same haplotype) → `cis`, dosage=1
- Single heterozygous variants → `het`, dosage=1

### Test 5: Phased mode with unphased data (should fail/warn)
Tests that phased mode rejects unphased data

**Expected behavior:**
- Warnings about unphased genotypes
- Unphased variants skipped or error

### Test 6-9: Additional phased mode features
Tests --show-variants, --score-map, --info-map, and --show-haplotype-scores flags

### Test 10: Phased mode with mixed phased/unphased input
Tests that phased mode correctly handles files with BOTH phased and unphased genotypes

**Expected behavior:**
- Warns about unphased genotypes (1/0, 0/1, 1/1)
- Skips unphased variants
- Only processes phased genotypes (0|1, 1|0, 1|1)

### Test 11: Unphased mode with mixed phased/unphased input
Tests that unphased mode accepts BOTH phased and unphased formats in the same file

**Expected behavior:**
- Accepts both 0|1 and 0/1 formats
- Treats all as unphased (ignores phase information)
- Processes all valid heterozygous and homozygous variants

### Test 12: Hemizygous/haploid genotypes
Tests handling of hemizygous genotypes (e.g., "./1" or "1/.", common for X/Y chromosomes in males)

**Expected behavior:**
- Warns about unphased/invalid genotype format
- Skips hemizygous genotypes (containing missing data ".")
- Processes any valid diploid genotypes in the file

### Test 13: Missing genotypes
Tests handling of missing data (./., .|., .)

**Expected behavior:**
- Warns about missing/invalid genotypes
- Skips missing data
- Processes any valid genotypes in the file

## Error Handling Tests (14-19)

These tests verify that the program properly fails with clear error messages for invalid input:

### Test 14: Empty genotype file
**Expected:** Program fails with "Genotype file is empty" error

### Test 15: Empty gene map file
**Expected:** Program fails with "Mapping file is empty" error

### Test 16: Genotype file with only header (no data)
**Expected:** Program fails with "No valid genotype data found" error

### Test 17: Gene map with only header (no data)
**Expected:** Program fails with "No matching variants found" error

### Test 18: All genotype lines have wrong column count
**Expected:** Program fails with "No valid genotype data found" error after reporting column count errors

### Test 19: Gene map with wrong column count
**Expected:** Program fails with error about insufficient columns in mapping file

**Note:** If only SOME genotype lines have wrong column counts, the program will warn and skip those lines but continue processing valid lines (defensive behavior).

## Running the Tests

### Prerequisites
1. Compile `call_chets`:
   ```bash
   cd /Users/flassen/Projects/06_call_chets/call_chets
   make  # or your build command
   ```

### Run all tests
```bash
cd /Users/flassen/Projects/06_call_chets/call_chets/tests
./run_tests.sh
```

### Run a single test manually
```bash
cd /Users/flassen/Projects/06_call_chets/call_chets
./call_chets --geno tests/genotypes_unphased.txt.gz \
             --gene-map tests/gene_map.txt.gz \
             --unphased
```

## Test Data Details

### Sample Data Structure

**SAMPLE1** (GENE1):
- chr1:1000:A:T - het or chet depending on mode
- chr1:2000:C:G - het or chet depending on mode

**SAMPLE2** (GENE1):
- chr1:1000:A:T - hom (1/1 or 1|1)

**SAMPLE3** (GENE2):
- chr1:3000:G:A - het or chet
- chr2:1000:T:C - het or chet

**SAMPLE4** (GENE3):
- chr2:2000:A:G - hom

### Expected Call Types by Mode

| Sample | Gene  | Phased Mode | Unphased Mode |
|--------|-------|-------------|---------------|
| SAMPLE1| GENE1 | chet        | het           |
| SAMPLE2| GENE1 | hom         | hom           |
| SAMPLE3| GENE2 | chet        | het           |
| SAMPLE4| GENE3 | hom         | hom           |

## Interpreting Results

### Success
```
====================================
  Test Summary
====================================
Tests passed: 19
Tests failed: 0
====================================
All tests passed!
```

### Failure
If tests fail, the script will show:
- Which test failed
- Expected output
- Actual output
- Differences between them

## Adding New Tests

To add a new test:

1. Create input genotype file in `tests/` directory
2. Gzip the file: `gzip filename.txt`
3. Create expected output file
4. Add test case to `run_tests.sh` using `run_test` function
5. Run `./run_tests.sh` to verify

## Notes

- All genotype files must be gzipped (`.gz` extension)
- Gene map file must also be gzipped
- Output format: `SAMPLE\tCHROM\tGENE\tCALL\tDOSAGE`
- The `--unphased` mode **ignores** all scoring/pathogenicity arguments
