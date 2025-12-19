#!/bin/bash

# Test for variant counting bug fix
# This test ensures that all processed variants are properly accounted for
# Either written to output or tracked in discard counters

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Set paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RECODE="${SCRIPT_DIR}/../bin/recode"

# Track test results
TESTS_PASSED=0
TESTS_FAILED=0

print_section() {
    echo ""
    echo "=========================================="
    echo "  $1"
    echo "=========================================="
}

# Create test VCF with multiple variants
create_test_vcf() {
    cat << 'EOF' > "${SCRIPT_DIR}/test_variant_counting.vcf"
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DS,Number=1,Type=Float,Description="Dosage">
##contig=<ID=chr1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1	Sample2	Sample3
chr1	1000	.	A	T	.	PASS	.	GT	0/0	0/1	1/1
chr1	2000	.	C	G	.	PASS	.	GT	0/0	0/1	1/1
chr1	3000	.	T	A	.	PASS	.	GT	0/1	1/1	1/1
chr1	4000	.	G	C	.	PASS	.	GT	0/0	0/0	1/1
chr1	5000	.	A	G	.	PASS	.	GT	0/0	1/1	1/1
EOF
}

# Create test VCF with no homozygous alternates (for dominance mode, should be discarded)
create_no_homalt_vcf() {
    cat << 'EOF' > "${SCRIPT_DIR}/test_variant_counting_no_homalt.vcf"
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1	Sample2	Sample3
chr1	1000	.	A	T	.	PASS	.	GT	0/0	0/1	1/1
chr1	2000	.	C	G	.	PASS	.	GT	0/0	0/1	0/1
chr1	3000	.	T	A	.	PASS	.	GT	0/0	0/0	0/1
chr1	4000	.	G	C	.	PASS	.	GT	0/1	1/1	1/1
EOF
}

# Test that variant counts add up correctly
test_variant_counting() {
    local test_name=$1
    local input_vcf=$2
    local mode=$3
    local expected_written=$4
    local expected_discarded=$5

    echo -e "\n${YELLOW}Running test: ${test_name}${NC}"

    output_file="${SCRIPT_DIR}/output_${test_name}.vcf"
    stderr_file="${SCRIPT_DIR}/stderr_${test_name}.txt"

    # Run recode and capture stderr
    $RECODE --input "${SCRIPT_DIR}/${input_vcf}" --mode "$mode" > "$output_file" 2> "$stderr_file"

    # Count variants in input (excluding header lines)
    input_count=$(grep -v "^#" "${SCRIPT_DIR}/${input_vcf}" | wc -l | tr -d ' ')

    # Count variants in output
    output_count=$(grep -v "^#" "$output_file" | wc -l | tr -d ' ')

    # Extract processing statistics from stderr
    processed=$(grep "Variants processed" "$stderr_file" | awk -F: '{print $2}' | tr -d ' ')
    kept=$(grep "Variants kept" "$stderr_file" | awk -F: '{print $2}' | tr -d ' ')

    # Extract discard count if present
    discarded=0
    if grep -q "Variants discarded" "$stderr_file"; then
        discarded=$(grep "Variants discarded" "$stderr_file" | head -1 | awk -F: '{print $2}' | tr -d ' ')
    fi

    echo "  Input variants:     $input_count"
    echo "  Processed:          $processed"
    echo "  Written (kept):     $kept"
    echo "  Discarded:          $discarded"
    echo "  Output line count:  $output_count"

    # Verify accounting: processed = kept + discarded
    sum=$((kept + discarded))

    local passed=true

    # Check 1: Processed should equal input
    if [ "$processed" -ne "$input_count" ]; then
        echo -e "${RED}  ✗ FAIL: Processed ($processed) != Input count ($input_count)${NC}"
        passed=false
    fi

    # Check 2: Kept + Discarded should equal Processed
    if [ "$sum" -ne "$processed" ]; then
        echo -e "${RED}  ✗ FAIL: Kept ($kept) + Discarded ($discarded) = $sum != Processed ($processed)${NC}"
        passed=false
    fi

    # Check 3: Output count should match kept count
    if [ "$output_count" -ne "$kept" ]; then
        echo -e "${RED}  ✗ FAIL: Output count ($output_count) != Kept ($kept)${NC}"
        passed=false
    fi

    # Check 4: If expected values provided, verify them
    if [ -n "$expected_written" ] && [ "$kept" -ne "$expected_written" ]; then
        echo -e "${RED}  ✗ FAIL: Expected $expected_written written, got $kept${NC}"
        passed=false
    fi

    if [ -n "$expected_discarded" ] && [ "$discarded" -ne "$expected_discarded" ]; then
        echo -e "${RED}  ✗ FAIL: Expected $expected_discarded discarded, got $discarded${NC}"
        passed=false
    fi

    if [ "$passed" = true ]; then
        echo -e "${GREEN}  ✓ PASS: All counts match correctly${NC}"
        ((TESTS_PASSED++))
    else
        echo -e "${RED}  ✗ FAIL: Count mismatch detected${NC}"
        echo "  Stderr output:"
        cat "$stderr_file"
        ((TESTS_FAILED++))
    fi
}

# Check if recode binary exists
print_section "Checking recode binary"
if [ ! -f "$RECODE" ]; then
    echo -e "${RED}Error: recode binary not found at $RECODE${NC}"
    echo "Please compile recode first with: make"
    exit 1
fi
echo -e "${GREEN}✓ Recode binary found${NC}"

# Create test data
print_section "Creating test data"
create_test_vcf
create_no_homalt_vcf
echo -e "${GREEN}✓ Test VCF files created${NC}"

print_section "Variant Counting Tests"

# Test 1: Recessive mode - all variants should be kept (no AA requirement)
test_variant_counting "recessive_all_kept" "test_variant_counting.vcf" "recessive" 5 0

# Test 2: Dominance mode - all variants should be kept (all have AA)
test_variant_counting "dominance_all_kept" "test_variant_counting.vcf" "dominance" 5 0

# Test 3: Dominance mode with some variants lacking homozygous alternates
test_variant_counting "dominance_some_discarded" "test_variant_counting_no_homalt.vcf" "dominance" 2 2

# Test 4: Recessive mode with no-homalt VCF - all should still be kept
test_variant_counting "recessive_no_homalt" "test_variant_counting_no_homalt.vcf" "recessive" 4 0

print_section "Cleanup"
echo "Removing temporary files..."
rm -f "${SCRIPT_DIR}"/output_*.vcf
rm -f "${SCRIPT_DIR}"/stderr_*.txt
rm -f "${SCRIPT_DIR}"/test_variant_counting*.vcf
echo -e "${GREEN}✓ Cleanup complete${NC}"

print_section "Test Summary"
echo -e "${GREEN}Tests passed: ${TESTS_PASSED}${NC}"
echo -e "${RED}Tests failed: ${TESTS_FAILED}${NC}"
echo "=========================================="

if [ $TESTS_FAILED -eq 0 ]; then
    echo -e "${GREEN}All variant counting tests passed!${NC}"
    exit 0
else
    echo -e "${RED}Some tests failed!${NC}"
    exit 1
fi
