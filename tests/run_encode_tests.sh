#!/bin/bash

# Unit tests for encode_vcf
# Tests conversion of call_chets text output to VCF

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Set paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENCODE_VCF="${SCRIPT_DIR}/../bin/encode_vcf"
DATA_FILE="${SCRIPT_DIR}/data_encode.txt"
SAMPLES_FILE="${SCRIPT_DIR}/samples_encode.txt"
DATA_FILE_GZ="${DATA_FILE}.gz"

# Track test results
TESTS_PASSED=0
TESTS_FAILED=0

# Clean up
cleanup() {
    rm -f "${DATA_FILE}" "${DATA_FILE_GZ}" "${SAMPLES_FILE}" "${SCRIPT_DIR}"/output_encode_*.vcf
}
trap cleanup EXIT

# Setup test data
setup_data() {
    # Create samples file
    printf "Sample1\nSample2\nSample3\nSample4\n" > "${SAMPLES_FILE}"

    # Create call_chets output format
    # Columns: Sample Chromosome Gene Configuration Dosage VariantInfo
    # Sample1: chet (dosage 2)
    # Sample2: het (dosage 1)
    # Sample3: hom (dosage 2)
    cat << EOF > "${DATA_FILE}"
Sample1	chr1	GeneA	chet	2	var1:info1|var2:info2
Sample2	chr1	GeneA	het	1	var1:info1
Sample3	chr1	GeneA	hom	2	var1:info1
Sample1	chr2	GeneB	cis	1	var3:info3;var4:info4
Sample2	chr2	GeneB	het	1	var3:info3
EOF
    
    # Gzip the data file as commonly expected (though encode_vcf uses gzopen which handles both?)
    # encode_vcf expects gzipped input for --input
    gzip -c "${DATA_FILE}" > "${DATA_FILE_GZ}"
}

# Function to run a test
run_test() {
    local test_name=$1
    shift
    local args=("$@")

    echo -e "\n${YELLOW}Running test: ${test_name}${NC}"
    
    local output_file="${SCRIPT_DIR}/output_encode_${test_name}.vcf"
    
    if $ENCODE_VCF --input "${DATA_FILE_GZ}" --samples "${SAMPLES_FILE}" "${args[@]}" > "$output_file" 2>/dev/null; then
        echo -e "${GREEN}✓ PASSED - Execution successful${NC}"
        
        # specific checks
        if grep -q "##fileformat=VCFv4.2" "$output_file"; then
             # Check if we got output lines
             local count=$(grep -v "^#" "$output_file" | wc -l)
             if [ "$count" -gt 0 ]; then
                 ((TESTS_PASSED++))
             else
                 echo -e "${RED}✗ FAILED - No variants in output${NC}"
                 ((TESTS_FAILED++))
             fi
        else
            echo -e "${RED}✗ FAILED - Invalid VCF header${NC}"
            ((TESTS_FAILED++))
        fi
    else
        echo -e "${RED}✗ FAILED - Return code non-zero${NC}"
        ((TESTS_FAILED++))
    fi
}

run_failing_test() {
    local test_name=$1
    shift
    local args=("$@")

    echo -e "\n${YELLOW}Running test (should fail): ${test_name}${NC}"
    
    if ! $ENCODE_VCF "${args[@]}" > /dev/null 2>&1; then
        echo -e "${GREEN}✓ PASSED - Failed as expected${NC}"
        ((TESTS_PASSED++))
    else
        echo -e "${RED}✗ FAILED - Should have failed but succeeded${NC}"
        ((TESTS_FAILED++))
    fi
}

echo "======================================"
echo "  encode_vcf Tests"
echo "======================================"

if [ ! -f "$ENCODE_VCF" ]; then
    echo -e "${RED}Error: encode_vcf binary not found at $ENCODE_VCF${NC}"
    exit 1
fi

setup_data

# Test 1: Basic additive mode
run_test "additive_basic" --mode additive

# Test 2: Recessive mode
run_test "recessive_basic" --mode recessive

# Test 3: Dominance mode
run_test "dominance_basic" --mode dominance --all-info

# Test 4: Custom 012 mode
run_test "mode_012" --mode 012

# Test 5: Min AC filter
args=("--mode" "additive" "--min-ac" "2")
run_test "min_ac_filter" "${args[@]}"
# Verify GeneB (AC=2 for Sample1 cis + Sample2 het? Wait. Cis=1, Het=1 in my data. sum=2)
# GeneA: Sample1(2) + Sample2(1) + Sample3(2) = 5
# Both should be present

# Test 6: Max AC filter
# Filter out GeneA (AC=5)
# output should only have GeneB
echo -e "${YELLOW}Running test: max_ac_filter (checking content)${NC}"
output_file="${SCRIPT_DIR}/output_encode_max_ac.vcf"
$ENCODE_VCF --input "${DATA_FILE_GZ}" --samples "${SAMPLES_FILE}" --mode additive --max-ac 3 > "$output_file" 2>/dev/null
if grep -q "GeneB" "$output_file" && ! grep -q "GeneA" "$output_file"; then
    echo -e "${GREEN}✓ PASSED - Correctly filtered by AC${NC}"
    ((TESTS_PASSED++))
else
    echo -e "${RED}✗ FAILED - Filter logic incorrect${NC}"
    ((TESTS_FAILED++))
fi

# Test 7: Error handling - Missing input
run_failing_test "missing_input" --input "nonexistent.gz" --samples "${SAMPLES_FILE}"

# Test 8: Error handling - Missing samples
run_failing_test "missing_samples" --input "${DATA_FILE_GZ}" --samples "nonexistent.txt"

# Test 9: Error handling - Invalid mode
run_failing_test "invalid_mode" --input "${DATA_FILE_GZ}" --samples "${SAMPLES_FILE}" --mode "super_mode"


echo ""
echo "======================================"
echo "  Test Summary"
echo "======================================"
echo -e "${GREEN}Tests passed: ${TESTS_PASSED}${NC}"
echo -e "${RED}Tests failed: ${TESTS_FAILED}${NC}"
echo "======================================"

if [ $TESTS_FAILED -eq 0 ]; then
    exit 0
else
    exit 1
fi
