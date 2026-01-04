#!/bin/bash

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RECODE="${SCRIPT_DIR}/../bin/recode"
INPUT_VCF="${SCRIPT_DIR}/test_min_hom_count.vcf"

# Setup: Create VCF with variants 0, 1, 2 HomAlt
# Var1: 0 HomAlt (S1=0/0, S2=0/0, S3=0/1, S4=0/1) -> AA=0
# Var2: 1 HomAlt (S1=0/0, S2=1/1, S3=0/1, S4=0/0) -> AA=1
# Var3: 2 HomAlt (S1=1/1, S2=1/1, S3=0/1, S4=0/0) -> AA=2
cat << EOF > "$INPUT_VCF"
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4
chr1	100	Var1_0AA	A	T	.	PASS	.	GT	0/0	0/0	0/1	0/1
chr1	200	Var2_1AA	A	T	.	PASS	.	GT	0/0	1/1	0/1	0/0
chr1	300	Var3_2AA	A	T	.	PASS	.	GT	1/1	1/1	0/1	0/0
EOF

TESTS_PASSED=0
TESTS_FAILED=0

run_test() {
    local name=$1
    local min_hom=$2
    local expected_count=$3
    local should_fail=$4

    echo -e "\n${YELLOW}Running Test: $name (Min Hom: $min_hom)${NC}"
    
    OUTPUT_FILE="${SCRIPT_DIR}/output_${name}.vcf"
    
    CMD="$RECODE --input $INPUT_VCF --mode dominance --min-hom-count $min_hom"
    
    if [ "$should_fail" = "true" ]; then
        if $CMD > "$OUTPUT_FILE" 2>/dev/null; then
            echo -e "${RED}FAILED: Command succeeded but should have failed${NC}"
            ((TESTS_FAILED++))
        else
            echo -e "${GREEN}PASSED: Command failed as expected${NC}"
            ((TESTS_PASSED++))
        fi
        return
    fi

    $CMD > "$OUTPUT_FILE" 2>/dev/null
    
    if [ $? -ne 0 ]; then
         echo -e "${RED}FAILED: Command returned error${NC}"
         ((TESTS_FAILED++))
         return
    fi
 
    count=$(grep -v "^#" "$OUTPUT_FILE" | wc -l | tr -d ' ')
    
    if [ "$count" -eq "$expected_count" ]; then
        echo -e "${GREEN}PASSED: Got $count variants (Expected: $expected_count)${NC}"
        ((TESTS_PASSED++))
    else
        echo -e "${RED}FAILED: Got $count variants (Expected: $expected_count)${NC}"
        ((TESTS_FAILED++))
    fi
}

# Test 1: Default (min 1)
run_test "Default" 1 2 "false"

# Test 2: Min 2
run_test "Min2" 2 1 "false"

# Test 3: Min 3
run_test "Min3" 3 0 "false"

# Test 4: Min 0 (Should fail)
run_test "Min0" 0 0 "true"

echo "--------------------------------"
if [ $TESTS_FAILED -eq 0 ]; then
    echo -e "${GREEN}ALL TESTS PASSED ($TESTS_PASSED)${NC}"
    exit 0
else
    echo -e "${RED}SOME TESTS FAILED ($TESTS_FAILED)${NC}"
    exit 1
fi
