#!/bin/bash

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Set paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CALL_CHETS="../call_chets"

# Track test results
TESTS_PASSED=0
TESTS_FAILED=0

# Function to run a test
run_test() {
    local test_name=$1
    local genotype_file=$2
    local expected_output=$3
    local extra_args=$4

    echo -e "\n${YELLOW}Running test: ${test_name}${NC}"

    # Run call_chets
    output_file="${SCRIPT_DIR}/output_${test_name}.txt"

    if $CALL_CHETS --geno "${SCRIPT_DIR}/${genotype_file}" \
                   --gene-map "${SCRIPT_DIR}/gene_map.txt.gz" \
                   $extra_args > "$output_file" 2>/dev/null; then

        # Compare with expected output
        if diff -q "$output_file" "${SCRIPT_DIR}/${expected_output}" > /dev/null 2>&1; then
            echo -e "${GREEN}✓ PASSED${NC}"
            ((TESTS_PASSED++))
        else
            echo -e "${RED}✗ FAILED - Output differs from expected${NC}"
            echo "Expected:"
            cat "${SCRIPT_DIR}/${expected_output}"
            echo ""
            echo "Got:"
            cat "$output_file"
            ((TESTS_FAILED++))
        fi
    else
        echo -e "${RED}✗ FAILED - call_chets returned non-zero exit code${NC}"
        ((TESTS_FAILED++))
    fi

    # Clean up output file
    rm -f "$output_file"
}

# Function to run a test that should fail
run_failing_test() {
    local test_name=$1
    local genotype_file=$2
    local extra_args=$3

    echo -e "\n${YELLOW}Running test (should fail): ${test_name}${NC}"

    # Run call_chets - expecting it to fail
    if ! $CALL_CHETS --geno "${SCRIPT_DIR}/${genotype_file}" \
                     --gene-map "${SCRIPT_DIR}/gene_map.txt.gz" \
                     $extra_args > /dev/null 2>&1; then
        echo -e "${GREEN}✓ PASSED (failed as expected)${NC}"
        ((TESTS_PASSED++))
    else
        echo -e "${RED}✗ FAILED - Should have failed but succeeded${NC}"
        ((TESTS_FAILED++))
    fi
}

# Function to run a test with custom gene map that should fail
run_failing_test_custom_map() {
    local test_name=$1
    local genotype_file=$2
    local gene_map_file=$3
    local extra_args=$4

    echo -e "\n${YELLOW}Running test (should fail): ${test_name}${NC}"

    # Run call_chets - expecting it to fail
    if ! $CALL_CHETS --geno "${SCRIPT_DIR}/${genotype_file}" \
                     --gene-map "${SCRIPT_DIR}/${gene_map_file}" \
                     $extra_args > /dev/null 2>&1; then
        echo -e "${GREEN}✓ PASSED (failed as expected)${NC}"
        ((TESTS_PASSED++))
    else
        echo -e "${RED}✗ FAILED - Should have failed but succeeded${NC}"
        ((TESTS_FAILED++))
    fi
}

echo "======================================"
echo "  call_chets --unphased Mode Tests"
echo "======================================"

# Check if call_chets exists
if [ ! -f "$CALL_CHETS" ]; then
    echo -e "${RED}Error: call_chets binary not found at $CALL_CHETS${NC}"
    echo "Please compile call_chets first"
    exit 1
fi

# Test 1: Unphased mode with unphased data
run_test "unphased_mode_unphased_data" "genotypes_unphased.txt.gz" "expected_unphased_mode_unphased_data.txt" "--unphased"

# Test 2: Unphased mode with phased data (should work and treat as unphased)
run_test "unphased_mode_phased_data" "genotypes_phased.txt.gz" "expected_unphased_mode_phased_data.txt" "--unphased"

# Test 3: Unphased mode with mixed data (should work)
run_test "unphased_mode_mixed_data" "genotypes_mixed.txt.gz" "expected_unphased_mode_phased_data.txt" "--unphased"

# Test 4: Phased mode (default) with phased data
run_test "phased_mode_phased_data" "genotypes_phased.txt.gz" "expected_phased_mode_phased_data.txt" ""

# Test 5: Phased mode with unphased data (should warn and skip unphased variants)
run_failing_test "phased_mode_unphased_data_should_warn" "genotypes_unphased.txt.gz" ""

# Test 6: Phased mode with --show-variants
run_test "phased_mode_show_variants" "genotypes_phased.txt.gz" "expected_phased_mode_show_variants.txt" "--show-variants"

# Test 7: Phased mode with --score-map
run_test "phased_mode_score_map" "genotypes_phased.txt.gz" "expected_phased_mode_score_map.txt" "--score-map score_map.txt.gz"

# Test 8: Phased mode with --info-map
run_test "phased_mode_info_map" "genotypes_phased.txt.gz" "expected_phased_mode_info_map.txt" "--info-map info_map.txt.gz"

# Test 9: Phased mode with --score-map and --show-haplotype-scores
run_test "phased_mode_haplotype_scores" "genotypes_phased.txt.gz" "expected_phased_mode_haplotype_scores.txt" "--score-map score_map.txt.gz --show-haplotype-scores"

# Test 10: Phased mode with mixed phased/unphased input (should skip unphased)
run_test "phased_mode_mixed_input" "genotypes_mixed_phased_unphased.txt.gz" "expected_phased_mode_mixed_input.txt" ""

# Test 11: Unphased mode with mixed phased/unphased input (should accept both)
run_test "unphased_mode_mixed_input" "genotypes_mixed_phased_unphased.txt.gz" "expected_unphased_mode_mixed_input.txt" "--unphased"

# Test 12: Haploid genotypes (should skip with warnings and process diploid)
run_test "haploid_genotypes" "genotypes_haploid.txt.gz" "expected_haploid_genotypes.txt" ""

# Test 13: Missing genotypes (should skip missing and process valid)
run_test "missing_genotypes" "genotypes_missing.txt.gz" "expected_missing_genotypes.txt" ""

echo ""
echo "======================================"
echo "  Error Handling Tests"
echo "======================================"

# Test 14: Empty genotype file (should fail)
run_failing_test "error_empty_genotype_file" "genotypes_empty.txt.gz" ""

# Test 15: Empty gene map file (should fail)
run_failing_test_custom_map "error_empty_gene_map" "genotypes_phased.txt.gz" "gene_map_empty.txt.gz" ""

# Test 16: Genotype file with only header (should fail)
run_failing_test "error_genotype_no_data" "genotypes_no_data.txt.gz" ""

# Test 17: Gene map with only header (should fail)
run_failing_test_custom_map "error_gene_map_no_data" "genotypes_phased.txt.gz" "gene_map_no_data.txt.gz" ""

# Test 18: ALL genotype lines have wrong column count (should fail)
run_failing_test "error_all_wrong_columns" "genotypes_all_wrong_columns.txt.gz" ""

# Test 19: Gene map with wrong column count (should fail)
run_failing_test_custom_map "error_gene_map_wrong_columns" "genotypes_phased.txt.gz" "gene_map_wrong_columns.txt.gz" ""

echo ""

echo "======================================

"
echo "  Test Summary"
echo "======================================"
echo -e "${GREEN}Tests passed: ${TESTS_PASSED}${NC}"
echo -e "${RED}Tests failed: ${TESTS_FAILED}${NC}"
echo "======================================"

if [ $TESTS_FAILED -eq 0 ]; then
    echo -e "${GREEN}All tests passed!${NC}"
    exit 0
else
    echo -e "${RED}Some tests failed!${NC}"
    exit 1
fi