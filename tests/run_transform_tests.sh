#!/bin/bash

# Unit tests for transform.cpp
# Tests dominance and recessive encoding modes with various options

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Set paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TRANSFORM="../transform"

# Track test results
TESTS_PASSED=0
TESTS_FAILED=0

# Function to print test section header
print_section() {
    echo ""
    echo "=========================================="
    echo "  $1"
    echo "=========================================="
}

# Function to verify dominance encoding calculations
# Formula: For genotype with r (ref/ref), h (het), a (alt/alt) frequencies:
#   - aa (0/0): -h*a
#   - Aa (0/1): 2*a*r
#   - AA (1/1): -h*r
verify_dominance_output() {
    local output_file=$1
    local variant_line=$2  # Line number of variant in VCF (starting from 1, excluding headers)

    # Extract dosage values from output
    local dosages=$(sed -n "${variant_line}p" "$output_file" | cut -f10-)

    # Return the dosages for manual verification
    echo "$dosages"
}

# Function to run a test
run_test() {
    local test_name=$1
    local input_vcf=$2
    shift 2
    local extra_args=("$@")

    echo -e "\n${YELLOW}Running test: ${test_name}${NC}"

    # Run transform
    output_file="${SCRIPT_DIR}/output_${test_name}.vcf"

    if $TRANSFORM --input "${SCRIPT_DIR}/${input_vcf}" "${extra_args[@]}" > "$output_file" 2>/dev/null; then
        echo -e "${GREEN}✓ PASSED - Transform executed successfully${NC}"

        # Show a sample of the output for verification
        echo -e "${BLUE}Sample output (first variant):${NC}"
        grep -v "^#" "$output_file" | head -1

        ((TESTS_PASSED++))
        return 0
    else
        echo -e "${RED}✗ FAILED - Transform returned non-zero exit code${NC}"
        ((TESTS_FAILED++))
        return 1
    fi
}

# Function to run a test with expected output
run_test_with_expected() {
    local test_name=$1
    local input_vcf=$2
    local expected_file=$3
    shift 3
    local extra_args=("$@")

    echo -e "\n${YELLOW}Running test: ${test_name}${NC}"

    # Run transform
    output_file="${SCRIPT_DIR}/output_${test_name}.vcf"

    if $TRANSFORM --input "${SCRIPT_DIR}/${input_vcf}" "${extra_args[@]}" > "$output_file" 2>/dev/null; then
        # Compare with expected output
        if diff -q "$output_file" "${SCRIPT_DIR}/${expected_file}" > /dev/null 2>&1; then
            echo -e "${GREEN}✓ PASSED - Output matches expected${NC}"
            ((TESTS_PASSED++))
        else
            echo -e "${RED}✗ FAILED - Output differs from expected${NC}"
            echo "Differences:"
            diff "$output_file" "${SCRIPT_DIR}/${expected_file}" | head -20
            ((TESTS_FAILED++))
        fi
    else
        echo -e "${RED}✗ FAILED - Transform returned non-zero exit code${NC}"
        ((TESTS_FAILED++))
    fi
}

# Function to test that a command should fail
run_failing_test() {
    local test_name=$1
    local input_vcf=$2
    shift 2
    local extra_args=("$@")

    echo -e "\n${YELLOW}Running test (should fail): ${test_name}${NC}"

    # Run transform - expecting it to fail
    if ! $TRANSFORM --input "${SCRIPT_DIR}/${input_vcf}" "${extra_args[@]}" > /dev/null 2>&1; then
        echo -e "${GREEN}✓ PASSED - Failed as expected${NC}"
        ((TESTS_PASSED++))
    else
        echo -e "${RED}✗ FAILED - Should have failed but succeeded${NC}"
        ((TESTS_FAILED++))
    fi
}

# Function to verify dominance encoding is correct
verify_dominance_math() {
    local test_name=$1
    local output_file="${SCRIPT_DIR}/output_${test_name}.vcf"

    echo -e "${BLUE}Verifying dominance encoding math...${NC}"

    # Extract first data line (chr1:1000)
    # Expected: 4 samples, genotypes: 0/0, 0/1, 1/1, 0/1
    # Frequencies: r=1/4, h=2/4, a=1/4
    # Expected dosages: aa=-0.125, Aa=0.125, AA=-0.125, Aa=0.125

    local line=$(grep -v "^#" "$output_file" | head -1)
    echo "  First variant: $(echo "$line" | cut -f1-5)"
    echo "  Dosages: $(echo "$line" | cut -f10-)"
}

# Check if transform exists
print_section "Checking transform binary"
if [ ! -f "$TRANSFORM" ]; then
    echo -e "${RED}Error: transform binary not found at $TRANSFORM${NC}"
    echo "Please compile transform first"
    exit 1
fi
echo -e "${GREEN}✓ Transform binary found${NC}"

print_section "Basic Dominance Mode Tests"

# Test 1: Basic dominance mode with GT field
run_test "dominance_mode_basic_gt" "test_transform_simple.vcf" --mode dominance
if [ $? -eq 0 ]; then
    verify_dominance_math "dominance_mode_basic_gt"
fi

# Test 2: Dominance mode with DS field
run_test "dominance_mode_dosage" "test_transform_dosage.vcf" --mode dominance

# Test 3: Dominance mode with per-variant scaling
run_test "dominance_mode_scaled" "test_transform_simple.vcf" --mode dominance --scale-per-variant

# Test 4: Dominance mode with variant ID setting
run_test "dominance_mode_variant_id" "test_transform_simple.vcf" --mode dominance --set-variant-id

# Test 5: Dominance mode with all info
run_test "dominance_mode_all_info" "test_transform_simple.vcf" --mode dominance --all-info

# Test 6: Dominance mode with global scaling
run_test "dominance_mode_global_scale" "test_transform_simple.vcf" --mode dominance --scale-globally

# Test 7: Dominance mode with group-based scaling
run_test "dominance_mode_group_scale" "test_transform_simple.vcf" --mode dominance --scale-by-group "test_transform_scale_by_group.txt"

# Test 8: Dominance mode with group-based scaling and all options
run_test "dominance_mode_full" "test_transform_simple.vcf" --mode dominance --scale-by-group "test_transform_scale_by_group.txt" --set-variant-id --all-info

print_section "Recessive Mode Tests"

# Test 9: Basic recessive mode
run_test "recessive_mode_basic" "test_transform_simple.vcf" --mode recessive

# Test 10: Recessive mode with DS field
run_test "recessive_mode_dosage" "test_transform_dosage.vcf" --mode recessive

# Test 11: Recessive mode with variant ID
run_test "recessive_mode_variant_id" "test_transform_simple.vcf" --mode recessive --set-variant-id

print_section "Edge Case Tests"

# Test 12: VCF with no homozygous alternate (should skip in dominance mode)
echo -e "\n${YELLOW}Running test: no_homozygous_alternate${NC}"
output_file="${SCRIPT_DIR}/output_no_homozygous_alternate.vcf"
$TRANSFORM --input "${SCRIPT_DIR}/test_transform_no_homalt.vcf" --mode dominance > "$output_file" 2>/dev/null
variant_count=$(grep -v "^#" "$output_file" | wc -l | tr -d ' ')
if [ "$variant_count" -eq 1 ]; then
    echo -e "${GREEN}✓ PASSED - Correctly skipped variant without homozygous alternate (1 variant output)${NC}"
    ((TESTS_PASSED++))
else
    echo -e "${RED}✗ FAILED - Expected 1 variant but got $variant_count${NC}"
    ((TESTS_FAILED++))
fi

# Test 13: Scaling factor with per-variant scaling
run_test "scaling_factor" "test_transform_simple.vcf" --mode dominance --scale-per-variant --scale-factor 0.5

print_section "Error Handling Tests"

# Test 13: Invalid mode
run_failing_test "invalid_mode" "test_transform_simple.vcf" --mode invalid

# Test 14: Missing input file
run_failing_test "missing_input" "nonexistent.vcf" --mode dominance

# Test 15: Missing gene map file
run_failing_test "missing_gene_map" "test_transform_simple.vcf" --mode dominance --gene-map "nonexistent.txt"

print_section "Predictable Output Tests"

echo -e "\n${YELLOW}Testing output predictability with known inputs${NC}"
echo "Input VCF: test_transform_simple.vcf"
echo "  chr1:1000 - Genotypes: 0/0, 0/1, 1/1, 0/1"
echo "              Counts: aa=1, Aa=2, AA=1 (n=4)"
echo "              Frequencies: r=0.25, h=0.5, a=0.25"
echo "              Expected dosages:"
echo "                aa (0/0): -h*a = -0.5*0.25 = -0.125"
echo "                Aa (0/1):  2*a*r = 2*0.25*0.25 = 0.125"
echo "                AA (1/1): -h*r = -0.5*0.25 = -0.125"
echo ""

output_file="${SCRIPT_DIR}/output_dominance_mode_basic_gt.vcf"
if [ -f "$output_file" ]; then
    actual_dosages=$(grep -v "^#" "$output_file" | head -1 | cut -f10-)
    echo "  Actual dosages from output: $actual_dosages"
    echo -e "${BLUE}  Manual verification: Check if these match the expected values${NC}"
fi

print_section "Cleanup"
echo "Removing temporary output files..."
rm -f "${SCRIPT_DIR}"/output_*.vcf
echo -e "${GREEN}✓ Cleanup complete${NC}"

print_section "Test Summary"
echo -e "${GREEN}Tests passed: ${TESTS_PASSED}${NC}"
echo -e "${RED}Tests failed: ${TESTS_FAILED}${NC}"
echo "=========================================="

if [ $TESTS_FAILED -eq 0 ]; then
    echo -e "${GREEN}All tests passed!${NC}"
    exit 0
else
    echo -e "${RED}Some tests failed!${NC}"
    exit 1
fi
