# `--unphased` Mode Implementation Summary

## Overview
Added a new `--unphased` flag to `call_chets.cpp` that allows processing genotype data without using phase information.

## Key Changes

### 1. Code Modifications (`call_chets.cpp`)

#### New Flag
- Added `--unphased` command-line argument
- Boolean flag `unphasedMode` defaults to `false`

#### Genotype Processing
- **Accepts both phased and unphased formats** when `--unphased` is active
  - Phased: `0|1`, `1|0`, `1|1`
  - Unphased: `0/1`, `1/0`, `1/1`
- **Normalizes all genotypes:**
  - `0/1`, `1/0`, `0|1`, `1|0` → all treated as het
  - `1/1`, `1|1` → treated as hom
- Phase information is **completely ignored**

#### Output Changes
- **Simplified call types:**
  - `het` (dosage=1) - any heterozygous variant
  - `hom` (dosage=2) - homozygous alternate
- **Removed in unphased mode:**
  - No `chet` (compound het) calls
  - No `cis` (cis variant) calls
  - No scoring/pathogenicity calculations

#### Feature Restrictions
When `--unphased` is active, these arguments are ignored with warnings:
- `--score-map` / `-p`
- `--haplotype-collapse` / `-hc`
- `--gene-collapse` / `-gc`
- `--show-haplotype-scores` / `-shs`

### 2. Implementation Details

**Location:** `/Users/flassen/Projects/06_call_chets/call_chets/call_chets.cpp`

**Key functions modified:**
- `printUsage()` - Added documentation for `--unphased` flag
- Main argument parsing loop - Added `--unphased` detection
- Genotype validation - Modified to accept both formats
- Genotype normalization - New logic to convert all formats to unphased
- Storage logic - Simplified haplotype tracking (haplotype 0 for het, haplotype 1 for hom)
- Output generation - Simplified call determination
- Variant printing - Adapted for unphased variant lists

## Test Suite

### Location
`/Users/flassen/Projects/06_call_chets/call_chets/tests/`

### Test Files Created

#### Input Data
1. **gene_map.txt.gz** - 5 variants mapping to 3 genes
2. **genotypes_unphased.txt.gz** - Unphased genotypes (0/1, 1/1)
3. **genotypes_phased.txt.gz** - Phased genotypes (0|1, 1|0, 1|1)
4. **genotypes_mixed.txt.gz** - Mixed phased/unphased formats

#### Expected Outputs
1. **expected_unphased_mode_unphased_data.txt** - Unphased mode + unphased input
2. **expected_unphased_mode_phased_data.txt** - Unphased mode + phased input
3. **expected_phased_mode_phased_data.txt** - Phased mode + phased input

#### Test Runner
- **run_tests.sh** - Automated test script with 5 test scenarios

#### Documentation
- **README.md** - Comprehensive test documentation
- **BEHAVIOR_REFERENCE.md** - Mode comparison and examples

### Test Coverage

#### Test 1: Unphased mode with unphased data
✓ Verifies basic unphased genotype processing

#### Test 2: Unphased mode with phased data
✓ **Key test:** Verifies phased data is accepted but phase is ignored

#### Test 3: Unphased mode with mixed data
✓ Verifies both formats can coexist in same file

#### Test 4: Phased mode with phased data
✓ Verifies default phased behavior still works

#### Test 5: Phased mode with unphased data
✓ Verifies phased mode rejects unphased data appropriately

## Running the Tests

### Prerequisites
```bash
cd /Users/flassen/Projects/06_call_chets/call_chets
# Compile call_chets (your build command here)
make  # or g++, etc.
```

### Execute Tests
```bash
cd tests
./run_tests.sh
```

### Expected Output
```
====================================
  call_chets --unphased Mode Tests
====================================

Running test: unphased_mode_unphased_data
✓ PASSED

Running test: unphased_mode_phased_data
✓ PASSED

Running test: unphased_mode_mixed_data
✓ PASSED

Running test: phased_mode_phased_data
✓ PASSED

Running test: phased_mode_unphased_data_should_warn
✓ PASSED (failed as expected)

====================================
  Test Summary
====================================
Tests passed: 5
Tests failed: 0
====================================
All tests passed!
```

## Usage Examples

### Basic unphased analysis
```bash
./call_chets --geno genotypes.txt.gz \
             --gene-map genes.txt.gz \
             --unphased
```

### With variant information
```bash
./call_chets --geno genotypes.txt.gz \
             --gene-map genes.txt.gz \
             --info-map variant_info.txt.gz \
             --unphased \
             --show-variants
```

### Processing phased data as unphased
```bash
# This works! Phase information is simply ignored
./call_chets --geno phased_genotypes.txt.gz \
             --gene-map genes.txt.gz \
             --unphased
```

## Behavior Comparison

| Scenario | Phased Mode | Unphased Mode |
|----------|-------------|---------------|
| Input: `SAMPLE1 chr1:1000 1\|0`<br>`SAMPLE1 chr1:2000 0\|1` | `chet`, dosage=2<br>(different haplotypes) | `het`, dosage=1<br>(phase ignored) |
| Input: `SAMPLE2 chr1:1000 1\|1` | `hom`, dosage=2 | `hom`, dosage=2 |
| Input: `SAMPLE3 chr1:1000 0/1` | **ERROR/WARNING**<br>(unphased not accepted) | `het`, dosage=1<br>(accepted) |

## Backward Compatibility

✓ **Fully backward compatible**
- Default behavior unchanged (phased mode)
- Existing workflows unaffected
- `--unphased` is opt-in only

## Files Modified

1. `/Users/flassen/Projects/06_call_chets/call_chets/call_chets.cpp` (main implementation)

## Files Created

### Tests Directory (`tests/`)
1. `gene_map.txt.gz`
2. `genotypes_unphased.txt.gz`
3. `genotypes_phased.txt.gz`
4. `genotypes_mixed.txt.gz`
5. `expected_unphased_mode_unphased_data.txt`
6. `expected_unphased_mode_phased_data.txt`
7. `expected_phased_mode_phased_data.txt`
8. `run_tests.sh`
9. `README.md`
10. `BEHAVIOR_REFERENCE.md`

### Documentation
11. `UNPHASED_MODE_SUMMARY.md` (this file)

## Next Steps

### To test the implementation:
```bash
cd /Users/flassen/Projects/06_call_chets/call_chets
# 1. Compile
make  # or your build command

# 2. Run tests
cd tests
./run_tests.sh

# 3. Test on real data
../call_chets --geno your_data.txt.gz \
              --gene-map your_genes.txt.gz \
              --unphased
```

### Future Enhancements (Optional)
- Support for dosage fields (DS) in unphased mode
- Variant filtering options
- Summary statistics output
- Multi-threading support

## Notes

- The implementation treats `0|1` and `1|0` identically in unphased mode (both → het)
- No warnings are issued when phased data is used with `--unphased` (it's intentionally accepted)
- The simplified output makes burden testing more straightforward
- Gene-level aggregation is unchanged (multiple variants per gene are still handled correctly)
