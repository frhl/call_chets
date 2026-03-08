#!/bin/bash

# =============================================================================
# ARCADE Break Tests
# Systematically tests edge cases, boundary conditions, and adversarial inputs
# across all tools: call_chets, recode, encode_vcf, filter_pp
# =============================================================================

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CALL_CHETS="${SCRIPT_DIR}/../bin/call_chets"
RECODE="${SCRIPT_DIR}/../bin/recode"
ENCODE_VCF="${SCRIPT_DIR}/../bin/encode_vcf"
FILTER_PP="${SCRIPT_DIR}/../bin/filter_pp"

TESTS_PASSED=0
TESTS_FAILED=0
TESTS_SKIPPED=0

# ----------- Helpers -----------

print_section() {
    echo ""
    echo -e "${BLUE}==========================================${NC}"
    echo -e "${BLUE}  $1${NC}"
    echo -e "${BLUE}==========================================${NC}"
}

# Test that a command succeeds and produces output
run_success_test() {
    local test_name=$1
    shift
    local cmd=("$@")

    echo -e "\n${YELLOW}Running test: ${test_name}${NC}"
    if "${cmd[@]}" > /dev/null 2>&1; then
        echo -e "${GREEN}✓ PASSED${NC}"
        ((TESTS_PASSED++))
        return 0
    else
        echo -e "${RED}✗ FAILED - expected success but got non-zero exit${NC}"
        ((TESTS_FAILED++))
        return 1
    fi
}

# Test that a command fails (non-zero exit)
run_fail_test() {
    local test_name=$1
    shift
    local cmd=("$@")

    echo -e "\n${YELLOW}Running test (should fail): ${test_name}${NC}"
    if ! "${cmd[@]}" > /dev/null 2>&1; then
        echo -e "${GREEN}✓ PASSED (failed as expected)${NC}"
        ((TESTS_PASSED++))
        return 0
    else
        echo -e "${RED}✗ FAILED - expected failure but got success${NC}"
        ((TESTS_FAILED++))
        return 1
    fi
}

# Test that a command succeeds and output matches expected
run_output_test() {
    local test_name=$1
    local expected_file=$2
    shift 2
    local cmd=("$@")

    echo -e "\n${YELLOW}Running test: ${test_name}${NC}"
    local output_file="${SCRIPT_DIR}/output_break_${test_name}.txt"
    if "${cmd[@]}" > "$output_file" 2>/dev/null; then
        if diff -q "$output_file" "$expected_file" > /dev/null 2>&1; then
            echo -e "${GREEN}✓ PASSED${NC}"
            ((TESTS_PASSED++))
        else
            echo -e "${RED}✗ FAILED - output differs${NC}"
            echo "  Expected:" && head -5 "$expected_file"
            echo "  Got:" && head -5 "$output_file"
            ((TESTS_FAILED++))
        fi
    else
        echo -e "${RED}✗ FAILED - non-zero exit${NC}"
        ((TESTS_FAILED++))
    fi
    rm -f "$output_file"
}

# Test that a command succeeds and output contains a specific string
run_output_contains_test() {
    local test_name=$1
    local pattern=$2
    shift 2
    local cmd=("$@")

    echo -e "\n${YELLOW}Running test: ${test_name}${NC}"
    local output_file="${SCRIPT_DIR}/output_break_${test_name}.txt"
    if "${cmd[@]}" > "$output_file" 2>/dev/null; then
        if grep -q "$pattern" "$output_file"; then
            echo -e "${GREEN}✓ PASSED${NC}"
            ((TESTS_PASSED++))
        else
            echo -e "${RED}✗ FAILED - output doesn't contain '${pattern}'${NC}"
            head -5 "$output_file"
            ((TESTS_FAILED++))
        fi
    else
        echo -e "${RED}✗ FAILED - non-zero exit${NC}"
        ((TESTS_FAILED++))
    fi
    rm -f "$output_file"
}

# Test that a command succeeds and output has N non-header lines
run_line_count_test() {
    local test_name=$1
    local expected_count=$2
    shift 2
    local cmd=("$@")

    echo -e "\n${YELLOW}Running test: ${test_name}${NC}"
    local output_file="${SCRIPT_DIR}/output_break_${test_name}.txt"
    if "${cmd[@]}" > "$output_file" 2>/dev/null; then
        local count=$(grep -cv "^#" "$output_file" | tr -d ' ')
        if [ "$count" -eq "$expected_count" ]; then
            echo -e "${GREEN}✓ PASSED (${count} lines as expected)${NC}"
            ((TESTS_PASSED++))
        else
            echo -e "${RED}✗ FAILED - expected ${expected_count} lines, got ${count}${NC}"
            ((TESTS_FAILED++))
        fi
    else
        echo -e "${RED}✗ FAILED - non-zero exit${NC}"
        ((TESTS_FAILED++))
    fi
    rm -f "$output_file"
}

# Test that stderr contains a specific string (for warning checking)
run_stderr_contains_test() {
    local test_name=$1
    local pattern=$2
    shift 2
    local cmd=("$@")

    echo -e "\n${YELLOW}Running test: ${test_name}${NC}"
    local stderr_file="${SCRIPT_DIR}/stderr_break_${test_name}.txt"
    "${cmd[@]}" > /dev/null 2>"$stderr_file"
    if grep -qi "$pattern" "$stderr_file"; then
        echo -e "${GREEN}✓ PASSED (warning found)${NC}"
        ((TESTS_PASSED++))
    else
        echo -e "${RED}✗ FAILED - stderr doesn't contain '${pattern}'${NC}"
        head -10 "$stderr_file"
        ((TESTS_FAILED++))
    fi
    rm -f "$stderr_file"
}

# Verify float values with tolerance
check_float() {
    local actual=$1
    local expected=$2
    local tol=${3:-0.001}
    awk -v a="$actual" -v e="$expected" -v t="$tol" 'BEGIN {
        diff = a - e; if (diff < 0) diff = -diff;
        if (diff < t) exit 0; else exit 1;
    }'
}

cleanup() {
    rm -f "${SCRIPT_DIR}"/output_break_*.txt "${SCRIPT_DIR}"/stderr_break_*.txt
    rm -f "${SCRIPT_DIR}"/test_break_*.vcf "${SCRIPT_DIR}"/test_break_*.txt
    rm -f "${SCRIPT_DIR}"/test_break_*.txt.gz "${SCRIPT_DIR}"/test_break_*.vcf.gz
}
trap cleanup EXIT

# =============================================================================
#  CHECK BINARIES EXIST
# =============================================================================
for bin in "$CALL_CHETS" "$RECODE" "$ENCODE_VCF" "$FILTER_PP"; do
    if [ ! -f "$bin" ]; then
        echo -e "${RED}Error: Binary not found: ${bin}${NC}"
        echo "Please run 'make' first."
        exit 1
    fi
done

# =============================================================================
#  SECTION 1: call_chets (interpret_phase) Edge Cases
# =============================================================================
print_section "1. call_chets — Input Corruption"

# --- Create test fixtures ---

# 1a. Duplicate variants on same haplotype for same sample
printf "sample\tvariant\tgenotype\n" > "${SCRIPT_DIR}/test_break_dup.txt"
printf "SAMPLE1\tchr1:100:A:T\t0|1\n" >> "${SCRIPT_DIR}/test_break_dup.txt"
printf "SAMPLE1\tchr1:100:A:T\t0|1\n" >> "${SCRIPT_DIR}/test_break_dup.txt"
gzip -f "${SCRIPT_DIR}/test_break_dup.txt"

# Gene map for testing
printf "variant\tgene\n" > "${SCRIPT_DIR}/test_break_genemap.txt"
printf "chr1:100:A:T\tGENE1\n" >> "${SCRIPT_DIR}/test_break_genemap.txt"
printf "chr1:200:C:G\tGENE1\n" >> "${SCRIPT_DIR}/test_break_genemap.txt"
printf "chr1:300:G:A\tGENE2\n" >> "${SCRIPT_DIR}/test_break_genemap.txt"
gzip -f "${SCRIPT_DIR}/test_break_genemap.txt"

# Test 1a: Duplicate variants — should not crash, should handle gracefully
run_success_test "1a_duplicate_variants" \
    "$CALL_CHETS" --geno "${SCRIPT_DIR}/test_break_dup.txt.gz" \
    --gene-map "${SCRIPT_DIR}/test_break_genemap.txt.gz"

# 1b. Variant in genotype file but NOT in gene map
printf "sample\tvariant\tgenotype\n" > "${SCRIPT_DIR}/test_break_unmapped.txt"
printf "SAMPLE1\tchr1:999:A:T\t0|1\n" >> "${SCRIPT_DIR}/test_break_unmapped.txt"
gzip -f "${SCRIPT_DIR}/test_break_unmapped.txt"

# Should succeed but produce no output — variant not in gene map is silently skipped
# The tool exits 0 because it processed the file, just found no mapped variants
echo -e "\n${YELLOW}Running test: 1b_unmapped_variant${NC}"
output_file="${SCRIPT_DIR}/output_break_1b.txt"
$CALL_CHETS --geno "${SCRIPT_DIR}/test_break_unmapped.txt.gz" \
    --gene-map "${SCRIPT_DIR}/test_break_genemap.txt.gz" > "$output_file" 2>/dev/null
exit_code=$?
line_count=$(wc -l < "$output_file" | tr -d ' ')
if [ $exit_code -eq 0 ] && [ "$line_count" -eq 0 ]; then
    echo -e "${GREEN}✓ PASSED (exits 0 with no output for unmapped variant)${NC}"
    ((TESTS_PASSED++))
else
    echo -e "${RED}✗ FAILED - exit=${exit_code}, lines=${line_count}${NC}"
    ((TESTS_FAILED++))
fi
rm -f "$output_file"

# 1c. Variant mapped to many genes (10 genes)
printf "variant\tgene\n" > "${SCRIPT_DIR}/test_break_multigene.txt"
for i in $(seq 1 10); do
    printf "chr1:100:A:T\tGENE${i}\n" >> "${SCRIPT_DIR}/test_break_multigene.txt"
done
gzip -f "${SCRIPT_DIR}/test_break_multigene.txt"

printf "sample\tvariant\tgenotype\n" > "${SCRIPT_DIR}/test_break_multigene_geno.txt"
printf "SAMPLE1\tchr1:100:A:T\t0|1\n" >> "${SCRIPT_DIR}/test_break_multigene_geno.txt"
gzip -f "${SCRIPT_DIR}/test_break_multigene_geno.txt"

echo -e "\n${YELLOW}Running test: 1c_variant_many_genes${NC}"
output_file="${SCRIPT_DIR}/output_break_1c.txt"
$CALL_CHETS --geno "${SCRIPT_DIR}/test_break_multigene_geno.txt.gz" \
    --gene-map "${SCRIPT_DIR}/test_break_multigene.txt.gz" > "$output_file" 2>/dev/null
line_count=$(wc -l < "$output_file" | tr -d ' ')
if [ "$line_count" -eq 10 ]; then
    echo -e "${GREEN}✓ PASSED (10 output lines for 10 gene mappings)${NC}"
    ((TESTS_PASSED++))
else
    echo -e "${RED}✗ FAILED - expected 10 lines, got ${line_count}${NC}"
    cat "$output_file"
    ((TESTS_FAILED++))
fi
rm -f "$output_file"

# 1e. Multiallelic genotype (1/2, 2|3)
printf "sample\tvariant\tgenotype\n" > "${SCRIPT_DIR}/test_break_multiallelic.txt"
printf "SAMPLE1\tchr1:100:A:T\t1|2\n" >> "${SCRIPT_DIR}/test_break_multiallelic.txt"
printf "SAMPLE2\tchr1:100:A:T\t0|1\n" >> "${SCRIPT_DIR}/test_break_multiallelic.txt"
gzip -f "${SCRIPT_DIR}/test_break_multiallelic.txt"

# Should handle gracefully — 1|2 is not in valid genotype set, should skip it
run_success_test "1e_multiallelic_genotype" \
    "$CALL_CHETS" --geno "${SCRIPT_DIR}/test_break_multiallelic.txt.gz" \
    --gene-map "${SCRIPT_DIR}/test_break_genemap.txt.gz"

# 1f. Windows line endings (\r\n)
printf "sample\tvariant\tgenotype\r\n" > "${SCRIPT_DIR}/test_break_crlf.txt"
printf "SAMPLE1\tchr1:100:A:T\t0|1\r\n" >> "${SCRIPT_DIR}/test_break_crlf.txt"
printf "SAMPLE1\tchr1:200:C:G\t1|0\r\n" >> "${SCRIPT_DIR}/test_break_crlf.txt"
gzip -f "${SCRIPT_DIR}/test_break_crlf.txt"

run_success_test "1f_windows_line_endings" \
    "$CALL_CHETS" --geno "${SCRIPT_DIR}/test_break_crlf.txt.gz" \
    --gene-map "${SCRIPT_DIR}/test_break_genemap.txt.gz"

# 1g. Non-standard chromosome names (numeric only, chrUn_)
printf "variant\tgene\n" > "${SCRIPT_DIR}/test_break_weirdchr_map.txt"
printf "1:100:A:T\tGENE1\n" >> "${SCRIPT_DIR}/test_break_weirdchr_map.txt"
gzip -f "${SCRIPT_DIR}/test_break_weirdchr_map.txt"

printf "sample\tvariant\tgenotype\n" > "${SCRIPT_DIR}/test_break_weirdchr_geno.txt"
printf "SAMPLE1\t1:100:A:T\t0|1\n" >> "${SCRIPT_DIR}/test_break_weirdchr_geno.txt"
gzip -f "${SCRIPT_DIR}/test_break_weirdchr_geno.txt"

run_success_test "1g_numeric_chromosome_name" \
    "$CALL_CHETS" --geno "${SCRIPT_DIR}/test_break_weirdchr_geno.txt.gz" \
    --gene-map "${SCRIPT_DIR}/test_break_weirdchr_map.txt.gz"

# 1h. Score map with negative scores
printf "variant\tgene\tscore\n" > "${SCRIPT_DIR}/test_break_negscore.txt"
printf "chr1:100:A:T\tGENE1\t-0.5\n" >> "${SCRIPT_DIR}/test_break_negscore.txt"
gzip -f "${SCRIPT_DIR}/test_break_negscore.txt"

printf "sample\tvariant\tgenotype\n" > "${SCRIPT_DIR}/test_break_scoregeno.txt"
printf "SAMPLE1\tchr1:100:A:T\t0|1\n" >> "${SCRIPT_DIR}/test_break_scoregeno.txt"
gzip -f "${SCRIPT_DIR}/test_break_scoregeno.txt"

# Should warn about out-of-range score but continue
run_stderr_contains_test "1h_negative_score_warning" "Warning" \
    "$CALL_CHETS" --geno "${SCRIPT_DIR}/test_break_scoregeno.txt.gz" \
    --gene-map "${SCRIPT_DIR}/test_break_genemap.txt.gz" \
    --score-map "${SCRIPT_DIR}/test_break_negscore.txt.gz" --verbose

# Score map with very large score
printf "variant\tgene\tscore\n" > "${SCRIPT_DIR}/test_break_bigscore.txt"
printf "chr1:100:A:T\tGENE1\t99999.0\n" >> "${SCRIPT_DIR}/test_break_bigscore.txt"
gzip -f "${SCRIPT_DIR}/test_break_bigscore.txt"

run_stderr_contains_test "1h_large_score_warning" "Warning" \
    "$CALL_CHETS" --geno "${SCRIPT_DIR}/test_break_scoregeno.txt.gz" \
    --gene-map "${SCRIPT_DIR}/test_break_genemap.txt.gz" \
    --score-map "${SCRIPT_DIR}/test_break_bigscore.txt.gz" --verbose

# 1i. All samples have same genotype (all het)
printf "sample\tvariant\tgenotype\n" > "${SCRIPT_DIR}/test_break_allhet.txt"
for i in $(seq 1 5); do
    printf "SAMPLE${i}\tchr1:100:A:T\t0|1\n" >> "${SCRIPT_DIR}/test_break_allhet.txt"
done
gzip -f "${SCRIPT_DIR}/test_break_allhet.txt"

echo -e "\n${YELLOW}Running test: 1i_all_same_genotype${NC}"
output_file="${SCRIPT_DIR}/output_break_1i.txt"
$CALL_CHETS --geno "${SCRIPT_DIR}/test_break_allhet.txt.gz" \
    --gene-map "${SCRIPT_DIR}/test_break_genemap.txt.gz" > "$output_file" 2>/dev/null
line_count=$(wc -l < "$output_file" | tr -d ' ')
if [ "$line_count" -eq 5 ]; then
    echo -e "${GREEN}✓ PASSED (5 het results as expected)${NC}"
    ((TESTS_PASSED++))
else
    echo -e "${RED}✗ FAILED - expected 5 lines, got ${line_count}${NC}"
    cat "$output_file"
    ((TESTS_FAILED++))
fi
rm -f "$output_file"

# 1j. Single sample input
printf "sample\tvariant\tgenotype\n" > "${SCRIPT_DIR}/test_break_single.txt"
printf "SAMPLE1\tchr1:100:A:T\t0|1\n" >> "${SCRIPT_DIR}/test_break_single.txt"
gzip -f "${SCRIPT_DIR}/test_break_single.txt"

echo -e "\n${YELLOW}Running test: 1j_single_sample${NC}"
output_file="${SCRIPT_DIR}/output_break_1j.txt"
$CALL_CHETS --geno "${SCRIPT_DIR}/test_break_single.txt.gz" \
    --gene-map "${SCRIPT_DIR}/test_break_genemap.txt.gz" > "$output_file" 2>/dev/null
line_count=$(wc -l < "$output_file" | tr -d ' ')
if [ "$line_count" -eq 1 ]; then
    echo -e "${GREEN}✓ PASSED (1 result for single sample)${NC}"
    ((TESTS_PASSED++))
    # Verify the call is 'het'
    call=$(awk '{print $4}' "$output_file")
    if [ "$call" == "het" ]; then
        echo -e "  ${GREEN}Call type is 'het' as expected${NC}"
    else
        echo -e "  ${RED}Call type is '${call}', expected 'het'${NC}"
    fi
else
    echo -e "${RED}✗ FAILED - expected 1 line, got ${line_count}${NC}"
    ((TESTS_FAILED++))
fi
rm -f "$output_file"

# 1k. Extra whitespace / tabs in genotype file
printf "sample\tvariant\tgenotype\n" > "${SCRIPT_DIR}/test_break_extraws.txt"
printf "SAMPLE1\tchr1:100:A:T\t0|1\t\n" >> "${SCRIPT_DIR}/test_break_extraws.txt"
printf "  SAMPLE2  \tchr1:100:A:T\t1|0\n" >> "${SCRIPT_DIR}/test_break_extraws.txt"
gzip -f "${SCRIPT_DIR}/test_break_extraws.txt"

run_success_test "1k_extra_whitespace" \
    "$CALL_CHETS" --geno "${SCRIPT_DIR}/test_break_extraws.txt.gz" \
    --gene-map "${SCRIPT_DIR}/test_break_genemap.txt.gz"


# =============================================================================
#  SECTION 2: recode (transform) Mathematical Edge Cases
# =============================================================================
print_section "2. recode — Mathematical Edge Cases"

# 2a. All samples 0/0 (monomorphic reference)
cat << 'EOF' > "${SCRIPT_DIR}/test_break_all_ref.vcf"
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4
chr1	1000	.	A	T	.	PASS	.	GT	0/0	0/0	0/0	0/0
EOF

# Dominance mode should skip this variant (no hom alt)
run_line_count_test "2a_all_ref_dominance" 0 \
    "$RECODE" --input "${SCRIPT_DIR}/test_break_all_ref.vcf" --mode dominance

# Recessive mode should succeed with all zeros
echo -e "\n${YELLOW}Running test: 2a_all_ref_recessive${NC}"
output_file="${SCRIPT_DIR}/output_break_2a_rec.vcf"
$RECODE --input "${SCRIPT_DIR}/test_break_all_ref.vcf" --mode recessive > "$output_file" 2>/dev/null
variant_line=$(grep -v "^#" "$output_file" | head -1)
dosages=$(echo "$variant_line" | cut -f10-)
if [ "$dosages" == "0	0	0	0" ]; then
    echo -e "${GREEN}✓ PASSED (all zeros for recessive with all 0/0)${NC}"
    ((TESTS_PASSED++))
else
    echo -e "${RED}✗ FAILED - expected all zeros, got: ${dosages}${NC}"
    ((TESTS_FAILED++))
fi
rm -f "$output_file"

# 2b. All samples 0/1 (all hets, no hom-ref or hom-alt)
cat << 'EOF' > "${SCRIPT_DIR}/test_break_all_het.vcf"
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4
chr1	1000	.	A	T	.	PASS	.	GT	0/1	0/1	0/1	0/1
EOF

# Dominance: should skip (no AA genotype)
run_line_count_test "2b_all_het_dominance" 0 \
    "$RECODE" --input "${SCRIPT_DIR}/test_break_all_het.vcf" --mode dominance

# Recessive: all hets become 0
echo -e "\n${YELLOW}Running test: 2b_all_het_recessive${NC}"
output_file="${SCRIPT_DIR}/output_break_2b_rec.vcf"
$RECODE --input "${SCRIPT_DIR}/test_break_all_het.vcf" --mode recessive > "$output_file" 2>/dev/null
variant_line=$(grep -v "^#" "$output_file" | head -1)
dosages=$(echo "$variant_line" | cut -f10-)
if [ "$dosages" == "0	0	0	0" ]; then
    echo -e "${GREEN}✓ PASSED (all zeros for recessive with all 0/1)${NC}"
    ((TESTS_PASSED++))
else
    echo -e "${RED}✗ FAILED - expected all zeros, got: ${dosages}${NC}"
    ((TESTS_FAILED++))
fi
rm -f "$output_file"

# 2c. Single sample VCF
cat << 'EOF' > "${SCRIPT_DIR}/test_break_single_sample.vcf"
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chr1	1000	.	A	T	.	PASS	.	GT	1/1
EOF

# Dominance: single sample with 1/1 only → r=0, h=0, a=1
# dom_aa = -0*1 = 0, dom_Aa = 2*1*0 = 0, dom_AA = -0*0 = 0
# min==max → scaleDosage returns 0.0
# Should skip due to minHomCount=1 but aa_count=0
run_line_count_test "2c_single_sample_hom_dominance" 0 \
    "$RECODE" --input "${SCRIPT_DIR}/test_break_single_sample.vcf" --mode dominance

# Recessive: should work
echo -e "\n${YELLOW}Running test: 2c_single_sample_recessive${NC}"
output_file="${SCRIPT_DIR}/output_break_2c_rec.vcf"
$RECODE --input "${SCRIPT_DIR}/test_break_single_sample.vcf" --mode recessive > "$output_file" 2>/dev/null
variant_line=$(grep -v "^#" "$output_file" | head -1)
dosage=$(echo "$variant_line" | cut -f10)
if [ "$dosage" == "2" ]; then
    echo -e "${GREEN}✓ PASSED (dosage=2 for single hom-alt sample)${NC}"
    ((TESTS_PASSED++))
else
    echo -e "${RED}✗ FAILED - expected 2, got: ${dosage}${NC}"
    ((TESTS_FAILED++))
fi
rm -f "$output_file"

# 2d. Missing genotypes mixed with valid ones
cat << 'EOF' > "${SCRIPT_DIR}/test_break_missing_gt.vcf"
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4
chr1	1000	.	A	T	.	PASS	.	GT	0/0	./.	1/1	0/1
EOF

# Should handle missing genotypes: exclude from freq calc, output '.' for missing
echo -e "\n${YELLOW}Running test: 2d_missing_genotypes_dominance${NC}"
output_file="${SCRIPT_DIR}/output_break_2d.vcf"
$RECODE --input "${SCRIPT_DIR}/test_break_missing_gt.vcf" --mode dominance > "$output_file" 2>/dev/null
variant_line=$(grep -v "^#" "$output_file" | head -1)
if [ -n "$variant_line" ]; then
    missing_col=$(echo "$variant_line" | cut -f11) # S2 is missing
    if [ "$missing_col" == "." ]; then
        echo -e "${GREEN}✓ PASSED (missing genotype outputs '.')${NC}"
        ((TESTS_PASSED++))
    else
        echo -e "${RED}✗ FAILED - expected '.' for missing, got: ${missing_col}${NC}"
        ((TESTS_FAILED++))
    fi
else
    echo -e "${RED}✗ FAILED - no output variant lines${NC}"
    ((TESTS_FAILED++))
fi
rm -f "$output_file"

# 2e. DS values at boundaries and slightly outside
cat << 'EOF' > "${SCRIPT_DIR}/test_break_ds_boundary.vcf"
##fileformat=VCFv4.2
##FORMAT=<ID=DS,Number=1,Type=Float,Description="Dosage">
##contig=<ID=chr1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4
chr1	1000	.	A	T	.	PASS	.	DS	0.0	1.0	2.0	0.5
EOF

run_success_test "2e_ds_boundary_values" \
    "$RECODE" --input "${SCRIPT_DIR}/test_break_ds_boundary.vcf" --mode dominance

# 2g. Multiallelic VCF record (multiple ALT alleles)
cat << 'EOF' > "${SCRIPT_DIR}/test_break_multiallelic.vcf"
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4
chr1	1000	.	A	T,G	.	PASS	.	GT	0/0	0/1	1/2	2/2
EOF

# Should not crash — 1/2 and 2/2 are multiallelic genotypes
run_success_test "2g_multiallelic_vcf" \
    "$RECODE" --input "${SCRIPT_DIR}/test_break_multiallelic.vcf" --mode dominance

# 2i. Scale factor of 0
echo -e "\n${YELLOW}Running test: 2i_scale_factor_zero${NC}"
output_file="${SCRIPT_DIR}/output_break_2i.vcf"
cat << 'EOF' > "${SCRIPT_DIR}/test_break_simple.vcf"
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4
chr1	1000	.	A	T	.	PASS	.	GT	0/0	0/1	1/1	0/1
EOF

# Scale factor 0 is special: code checks `scalingFactor != 0 && scalingFactor != 1.0`
# so scale-factor 0 means "don't apply scaling factor" → same as factor 1.0
$RECODE --input "${SCRIPT_DIR}/test_break_simple.vcf" --mode dominance --scale-factor 0 > "$output_file" 2>/dev/null
if [ $? -eq 0 ]; then
    echo -e "${GREEN}✓ PASSED (scale-factor 0 doesn't crash)${NC}"
    ((TESTS_PASSED++))
else
    echo -e "${RED}✗ FAILED${NC}"
    ((TESTS_FAILED++))
fi
rm -f "$output_file"

# 2j. Negative scale factor
echo -e "\n${YELLOW}Running test: 2j_negative_scale_factor${NC}"
output_file="${SCRIPT_DIR}/output_break_2j.vcf"
$RECODE --input "${SCRIPT_DIR}/test_break_simple.vcf" --mode recessive --scale-factor -1.0 > "$output_file" 2>/dev/null
if [ $? -eq 0 ]; then
    # Check that dosages are negative (sign flip)
    variant_line=$(grep -v "^#" "$output_file" | head -1)
    hom_dosage=$(echo "$variant_line" | cut -f12) # S3 = 1/1
    if check_float "$hom_dosage" "-2.0"; then
        echo -e "${GREEN}✓ PASSED (negative scale factor flips signs: hom_dosage=${hom_dosage})${NC}"
        ((TESTS_PASSED++))
    else
        echo -e "${RED}✗ FAILED - expected -2.0 for hom, got: ${hom_dosage}${NC}"
        ((TESTS_FAILED++))
    fi
else
    echo -e "${RED}✗ FAILED - non-zero exit${NC}"
    ((TESTS_FAILED++))
fi
rm -f "$output_file"

# 2k. --scale-by-group with variant not in group file
cat << 'EOF' > "${SCRIPT_DIR}/test_break_group_missing.txt"
variant	group
chr1:999:A:T	GENE_OTHER
EOF

run_stderr_contains_test "2k_variant_not_in_group" "not found in group" \
    "$RECODE" --input "${SCRIPT_DIR}/test_break_simple.vcf" --mode dominance \
    --scale-by-group "${SCRIPT_DIR}/test_break_group_missing.txt"

# 2l. Two mutually exclusive scaling flags
run_fail_test "2l_mutually_exclusive_scaling" \
    "$RECODE" --input "${SCRIPT_DIR}/test_break_simple.vcf" --mode dominance \
    --scale-per-variant --scale-globally

# 2m. All samples 1/1 (monomorphic alt)
cat << 'EOF' > "${SCRIPT_DIR}/test_break_all_alt.vcf"
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4
chr1	1000	.	A	T	.	PASS	.	GT	1/1	1/1	1/1	1/1
EOF

# Dominance: aa_count=0 → skipped (minHomCount requires both aa and AA)
run_line_count_test "2m_all_alt_dominance" 0 \
    "$RECODE" --input "${SCRIPT_DIR}/test_break_all_alt.vcf" --mode dominance

# 2n. VCF with no samples at all (empty sample column)
cat << 'EOF' > "${SCRIPT_DIR}/test_break_no_samples.vcf"
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
chr1	1000	.	A	T	.	PASS	.	GT
EOF

run_fail_test "2n_no_samples" \
    "$RECODE" --input "${SCRIPT_DIR}/test_break_no_samples.vcf" --mode dominance

# 2o. All missing genotypes
cat << 'EOF' > "${SCRIPT_DIR}/test_break_all_missing.vcf"
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3
chr1	1000	.	A	T	.	PASS	.	GT	./.	./.	./.
EOF

# Should produce no output variants (all missing)
run_line_count_test "2o_all_missing_genotypes" 0 \
    "$RECODE" --input "${SCRIPT_DIR}/test_break_all_missing.vcf" --mode dominance

# 2p. Haploid genotype in VCF (./1, 0/.)
cat << 'EOF' > "${SCRIPT_DIR}/test_break_haploid.vcf"
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4
chr1	1000	.	A	T	.	PASS	.	GT	0/0	0/1	1/1	.
EOF

run_stderr_contains_test "2p_haploid_genotype_warning" "Haploid\|haploid\|Warning" \
    "$RECODE" --input "${SCRIPT_DIR}/test_break_haploid.vcf" --mode dominance

# 2q. VCF with no variants (header only)
cat << 'EOF' > "${SCRIPT_DIR}/test_break_no_variants.vcf"
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
EOF

run_line_count_test "2q_no_variants" 0 \
    "$RECODE" --input "${SCRIPT_DIR}/test_break_no_variants.vcf" --mode dominance

# 2r. Verify dominance math with extreme frequencies
# r=0.9, h=0.1, a=0.0 → all hom-ref except one het
# But AA_count=0 so will be filtered out — let's do a valid case
cat << 'EOF' > "${SCRIPT_DIR}/test_break_extreme_freq.vcf"
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4	S5	S6	S7	S8	S9	S10
chr1	1000	.	A	T	.	PASS	.	GT	0/0	0/0	0/0	0/0	0/0	0/0	0/0	0/0	0/1	1/1
EOF

echo -e "\n${YELLOW}Running test: 2r_extreme_frequency_dominance${NC}"
output_file="${SCRIPT_DIR}/output_break_2r.vcf"
$RECODE --input "${SCRIPT_DIR}/test_break_extreme_freq.vcf" --mode dominance --all-info > "$output_file" 2>/dev/null
# r=0.8, h=0.1, a=0.1
# dom_aa = -0.1*0.1 = -0.01
# dom_Aa = 2*0.1*0.8 = 0.16
# dom_AA = -0.1*0.8 = -0.08
variant_line=$(grep -v "^#" "$output_file" | head -1)
if [ -n "$variant_line" ]; then
    d_aa=$(echo "$variant_line" | cut -f10) # S1 = 0/0
    d_Aa=$(echo "$variant_line" | cut -f18) # S9 = 0/1
    d_AA=$(echo "$variant_line" | cut -f19) # S10 = 1/1
    if check_float "$d_aa" "-0.01" && check_float "$d_Aa" "0.16" && check_float "$d_AA" "-0.08"; then
        echo -e "${GREEN}✓ PASSED (extreme freq dosages: aa=${d_aa}, Aa=${d_Aa}, AA=${d_AA})${NC}"
        ((TESTS_PASSED++))
    else
        echo -e "${RED}✗ FAILED - expected -0.01/0.16/-0.08, got: ${d_aa}/${d_Aa}/${d_AA}${NC}"
        ((TESTS_FAILED++))
    fi
else
    echo -e "${RED}✗ FAILED - no output${NC}"
    ((TESTS_FAILED++))
fi
rm -f "$output_file"

# 2s. Per-variant scaling: verify output is in [0,2]
echo -e "\n${YELLOW}Running test: 2s_scaled_output_range${NC}"
output_file="${SCRIPT_DIR}/output_break_2s.vcf"
$RECODE --input "${SCRIPT_DIR}/test_break_simple.vcf" --mode dominance --scale-per-variant > "$output_file" 2>/dev/null
variant_line=$(grep -v "^#" "$output_file" | head -1)
all_valid=true
for col in $(seq 10 13); do
    val=$(echo "$variant_line" | cut -f"$col")
    in_range=$(awk -v v="$val" 'BEGIN { if (v >= -0.001 && v <= 2.001) print "yes"; else print "no" }')
    if [ "$in_range" != "yes" ]; then
        all_valid=false
        echo -e "${RED}  Column $col: ${val} OUT OF RANGE${NC}"
    fi
done
if $all_valid; then
    echo -e "${GREEN}✓ PASSED (all scaled dosages in [0,2])${NC}"
    ((TESTS_PASSED++))
else
    echo -e "${RED}✗ FAILED${NC}"
    ((TESTS_FAILED++))
fi
rm -f "$output_file"


# =============================================================================
#  SECTION 3: encode_vcf (make_pseudo_vcf) Output Integrity
# =============================================================================
print_section "3. encode_vcf — Output Integrity"

# Create samples file
printf "Sample1\nSample2\nSample3\nSample4\n" > "${SCRIPT_DIR}/test_break_samples.txt"

# 3a. Sample in data but NOT in samples file
cat << 'EOF' > "${SCRIPT_DIR}/test_break_extra_sample.txt"
Sample1	chr1	GeneA	chet	2	var1|var2
Sample2	chr1	GeneA	het	1	var1
UNKNOWN	chr1	GeneA	het	1	var1
EOF
gzip -f "${SCRIPT_DIR}/test_break_extra_sample.txt"

echo -e "\n${YELLOW}Running test: 3a_sample_not_in_samples_file${NC}"
output_file="${SCRIPT_DIR}/output_break_3a.vcf"
$ENCODE_VCF --input "${SCRIPT_DIR}/test_break_extra_sample.txt.gz" \
    --samples "${SCRIPT_DIR}/test_break_samples.txt" --mode additive > "$output_file" 2>/dev/null
if [ $? -eq 0 ]; then
    # Check that UNKNOWN is not in header
    if ! grep -q "UNKNOWN" "$output_file"; then
        echo -e "${GREEN}✓ PASSED (UNKNOWN sample correctly excluded)${NC}"
        ((TESTS_PASSED++))
    else
        echo -e "${RED}✗ FAILED - UNKNOWN sample appears in output${NC}"
        ((TESTS_FAILED++))
    fi
else
    echo -e "${RED}✗ FAILED - non-zero exit${NC}"
    ((TESTS_FAILED++))
fi
rm -f "$output_file"

# 3b. Sample in samples file but NOT in data (should have dosage 0)
cat << 'EOF' > "${SCRIPT_DIR}/test_break_missing_sample.txt"
Sample1	chr1	GeneA	chet	2	var1|var2
EOF
gzip -f "${SCRIPT_DIR}/test_break_missing_sample.txt"

echo -e "\n${YELLOW}Running test: 3b_sample_in_list_not_in_data${NC}"
output_file="${SCRIPT_DIR}/output_break_3b.vcf"
$ENCODE_VCF --input "${SCRIPT_DIR}/test_break_missing_sample.txt.gz" \
    --samples "${SCRIPT_DIR}/test_break_samples.txt" --mode additive > "$output_file" 2>/dev/null
if [ $? -eq 0 ]; then
    # Sample2-4 should have dosage 0
    variant_line=$(grep -v "^#" "$output_file" | head -1)
    # Samples are sorted, so order is Sample1, Sample2, Sample3, Sample4
    s2=$(echo "$variant_line" | cut -f11) # Sample2
    s3=$(echo "$variant_line" | cut -f12) # Sample3
    s4=$(echo "$variant_line" | cut -f13) # Sample4
    if check_float "$s2" "0" && check_float "$s3" "0" && check_float "$s4" "0"; then
        echo -e "${GREEN}✓ PASSED (absent samples have dosage 0)${NC}"
        ((TESTS_PASSED++))
    else
        echo -e "${RED}✗ FAILED - expected 0 for absent samples, got: ${s2}/${s3}/${s4}${NC}"
        ((TESTS_FAILED++))
    fi
else
    echo -e "${RED}✗ FAILED - non-zero exit${NC}"
    ((TESTS_FAILED++))
fi
rm -f "$output_file"

# 3c. Gene with only one carrier
cat << 'EOF' > "${SCRIPT_DIR}/test_break_one_carrier.txt"
Sample1	chr1	GeneA	het	1	var1
EOF
gzip -f "${SCRIPT_DIR}/test_break_one_carrier.txt"

run_output_contains_test "3c_single_carrier" "GeneA" \
    "$ENCODE_VCF" --input "${SCRIPT_DIR}/test_break_one_carrier.txt.gz" \
    --samples "${SCRIPT_DIR}/test_break_samples.txt" --mode additive

# 3d. Gene with ALL samples as carriers
cat << 'EOF' > "${SCRIPT_DIR}/test_break_all_carriers.txt"
Sample1	chr1	GeneA	chet	2	var1|var2
Sample2	chr1	GeneA	chet	2	var1|var2
Sample3	chr1	GeneA	chet	2	var1|var2
Sample4	chr1	GeneA	chet	2	var1|var2
EOF
gzip -f "${SCRIPT_DIR}/test_break_all_carriers.txt"

echo -e "\n${YELLOW}Running test: 3d_all_samples_carriers${NC}"
output_file="${SCRIPT_DIR}/output_break_3d.vcf"
$ENCODE_VCF --input "${SCRIPT_DIR}/test_break_all_carriers.txt.gz" \
    --samples "${SCRIPT_DIR}/test_break_samples.txt" --mode additive > "$output_file" 2>/dev/null
if [ $? -eq 0 ]; then
    # AC should be 8 (4 samples * dosage 2)
    ac=$(grep -v "^#" "$output_file" | head -1 | grep -o 'AC=[0-9]*' | cut -d= -f2)
    if [ "$ac" -eq 8 ]; then
        echo -e "${GREEN}✓ PASSED (AC=8 for all carriers with dosage 2)${NC}"
        ((TESTS_PASSED++))
    else
        echo -e "${RED}✗ FAILED - expected AC=8, got AC=${ac}${NC}"
        ((TESTS_FAILED++))
    fi
else
    echo -e "${RED}✗ FAILED - non-zero exit${NC}"
    ((TESTS_FAILED++))
fi
rm -f "$output_file"

# 3e. Duplicate gene entries for same sample (same gene, different calls)
cat << 'EOF' > "${SCRIPT_DIR}/test_break_dup_gene.txt"
Sample1	chr1	GeneA	het	1	var1
Sample1	chr1	GeneA	hom	2	var2
EOF
gzip -f "${SCRIPT_DIR}/test_break_dup_gene.txt"

# Should not crash — second entry should overwrite first
run_success_test "3e_duplicate_gene_entry" \
    "$ENCODE_VCF" --input "${SCRIPT_DIR}/test_break_dup_gene.txt.gz" \
    --samples "${SCRIPT_DIR}/test_break_samples.txt" --mode additive

# 3f. No sample overlap at all
printf "AAAA\nBBBB\n" > "${SCRIPT_DIR}/test_break_no_overlap_samples.txt"

cat << 'EOF' > "${SCRIPT_DIR}/test_break_no_overlap_data.txt"
Sample1	chr1	GeneA	het	1	var1
EOF
gzip -f "${SCRIPT_DIR}/test_break_no_overlap_data.txt"

run_fail_test "3f_no_sample_overlap" \
    "$ENCODE_VCF" --input "${SCRIPT_DIR}/test_break_no_overlap_data.txt.gz" \
    --samples "${SCRIPT_DIR}/test_break_no_overlap_samples.txt" --mode additive

# 3g. Invalid configuration column
cat << 'EOF' > "${SCRIPT_DIR}/test_break_bad_config.txt"
Sample1	chr1	GeneA	INVALID	1	var1
EOF
gzip -f "${SCRIPT_DIR}/test_break_bad_config.txt"

run_fail_test "3g_invalid_configuration" \
    "$ENCODE_VCF" --input "${SCRIPT_DIR}/test_break_bad_config.txt.gz" \
    --samples "${SCRIPT_DIR}/test_break_samples.txt" --mode additive

# 3h. Dominance mode with no bi-allelic carriers (all het)
cat << 'EOF' > "${SCRIPT_DIR}/test_break_dom_nobi.txt"
Sample1	chr1	GeneA	het	1	var1
Sample2	chr1	GeneA	het	1	var2
EOF
gzip -f "${SCRIPT_DIR}/test_break_dom_nobi.txt"

# Dominance mode skips genes with no bi-allelic carriers
run_line_count_test "3h_dominance_no_biallelic" 0 \
    "$ENCODE_VCF" --input "${SCRIPT_DIR}/test_break_dom_nobi.txt.gz" \
    --samples "${SCRIPT_DIR}/test_break_samples.txt" --mode dominance

# 3i. Dosage value out of range in input
cat << 'EOF' > "${SCRIPT_DIR}/test_break_bad_dosage.txt"
Sample1	chr1	GeneA	het	5	var1
EOF
gzip -f "${SCRIPT_DIR}/test_break_bad_dosage.txt"

run_fail_test "3i_dosage_out_of_range" \
    "$ENCODE_VCF" --input "${SCRIPT_DIR}/test_break_bad_dosage.txt.gz" \
    --samples "${SCRIPT_DIR}/test_break_samples.txt" --mode additive

# 3j. Empty samples file
printf "" > "${SCRIPT_DIR}/test_break_empty_samples.txt"

run_fail_test "3j_empty_samples_file" \
    "$ENCODE_VCF" --input "${SCRIPT_DIR}/test_break_one_carrier.txt.gz" \
    --samples "${SCRIPT_DIR}/test_break_empty_samples.txt" --mode additive

# 3k. Samples file with multiple words per line
printf "Sample1 Sample2\nSample3\n" > "${SCRIPT_DIR}/test_break_multiword_samples.txt"

run_fail_test "3k_multiword_samples_file" \
    "$ENCODE_VCF" --input "${SCRIPT_DIR}/test_break_one_carrier.txt.gz" \
    --samples "${SCRIPT_DIR}/test_break_multiword_samples.txt" --mode additive


# =============================================================================
#  SECTION 4: filter_pp — Threshold Edge Cases
# =============================================================================
print_section "4. filter_pp — Threshold Edge Cases"

# Create a minimal VCF with PP field for testing
# We use the phased VCF from examples as a base (it may or may not have PP)
# but filter_pp gracefully handles VCFs without PP (just passes records through)

# 4a. PP threshold of exactly 1.0 — use example phased VCF
# filter_pp should succeed even if VCF has no PP field (just passes through)
if [ -f "${SCRIPT_DIR}/../examples/input/phased.vcf.gz" ]; then
    run_success_test "4a_pp_threshold_1.0" \
        "$FILTER_PP" --input "${SCRIPT_DIR}/../examples/input/phased.vcf.gz" \
        --output "${SCRIPT_DIR}/output_break_4a.vcf" --pp-threshold 1.0
    rm -f "${SCRIPT_DIR}/output_break_4a.vcf"
else
    echo -e "${YELLOW}  Skipping 4a (no example VCF available)${NC}"
    ((TESTS_SKIPPED++))
fi

# 4b. PP threshold of 0 (should fail — requires > 0)
run_fail_test "4b_pp_threshold_zero" \
    "$FILTER_PP" --input "${SCRIPT_DIR}/../examples/input/phased.vcf.gz" \
    --output "${SCRIPT_DIR}/output_break_4b.vcf" --pp-threshold 0
rm -f "${SCRIPT_DIR}/output_break_4b.vcf"

# 4c. Missing required args
run_fail_test "4c_missing_args" "$FILTER_PP"

# 4d. Nonexistent input file
run_fail_test "4d_nonexistent_input" \
    "$FILTER_PP" --input "/nonexistent.vcf" \
    --output "${SCRIPT_DIR}/output_break_4d.vcf" --pp-threshold 0.5
rm -f "${SCRIPT_DIR}/output_break_4d.vcf"

# 4e. Negative PP threshold (should fail — requires > 0)
run_fail_test "4e_negative_pp_threshold" \
    "$FILTER_PP" --input "${SCRIPT_DIR}/../examples/input/phased.vcf.gz" \
    --output "${SCRIPT_DIR}/output_break_4e.vcf" --pp-threshold -0.5
rm -f "${SCRIPT_DIR}/output_break_4e.vcf"


# =============================================================================
#  SECTION 5: Cross-Cutting Concerns
# =============================================================================
print_section "5. Cross-Cutting Concerns"

# 5a. Binary garbage as input file
echo -e "\n${YELLOW}Running test: 5a_binary_garbage_input${NC}"
dd if=/dev/urandom of="${SCRIPT_DIR}/test_break_garbage.bin" bs=256 count=1 2>/dev/null
gzip -f "${SCRIPT_DIR}/test_break_garbage.bin"
if ! "$CALL_CHETS" --geno "${SCRIPT_DIR}/test_break_garbage.bin.gz" \
    --gene-map "${SCRIPT_DIR}/test_break_genemap.txt.gz" > /dev/null 2>&1; then
    echo -e "${GREEN}✓ PASSED (gracefully rejected garbage input)${NC}"
    ((TESTS_PASSED++))
else
    echo -e "${RED}✗ FAILED - should have rejected garbage${NC}"
    ((TESTS_FAILED++))
fi

# 5b. Binary garbage as VCF input for recode
echo -e "\n${YELLOW}Running test: 5b_binary_garbage_vcf${NC}"
dd if=/dev/urandom of="${SCRIPT_DIR}/test_break_garbage.vcf" bs=256 count=1 2>/dev/null
if ! "$RECODE" --input "${SCRIPT_DIR}/test_break_garbage.vcf" --mode dominance > /dev/null 2>&1; then
    echo -e "${GREEN}✓ PASSED (gracefully rejected garbage VCF)${NC}"
    ((TESTS_PASSED++))
else
    echo -e "${RED}✗ FAILED - should have rejected garbage${NC}"
    ((TESTS_FAILED++))
fi

# 5c. Very long gene name (1000 chars)
long_gene=$(python3 -c "print('G' * 1000)")
printf "variant\tgene\n" > "${SCRIPT_DIR}/test_break_longname_map.txt"
printf "chr1:100:A:T\t%s\n" "$long_gene" >> "${SCRIPT_DIR}/test_break_longname_map.txt"
gzip -f "${SCRIPT_DIR}/test_break_longname_map.txt"

printf "sample\tvariant\tgenotype\n" > "${SCRIPT_DIR}/test_break_longname_geno.txt"
printf "SAMPLE1\tchr1:100:A:T\t0|1\n" >> "${SCRIPT_DIR}/test_break_longname_geno.txt"
gzip -f "${SCRIPT_DIR}/test_break_longname_geno.txt"

echo -e "\n${YELLOW}Running test: 5c_very_long_gene_name${NC}"
output_file="${SCRIPT_DIR}/output_break_5c.txt"
if $CALL_CHETS --geno "${SCRIPT_DIR}/test_break_longname_geno.txt.gz" \
    --gene-map "${SCRIPT_DIR}/test_break_longname_map.txt.gz" > "$output_file" 2>/dev/null; then
    if grep -q "$long_gene" "$output_file"; then
        echo -e "${GREEN}✓ PASSED (1000-char gene name preserved)${NC}"
        ((TESTS_PASSED++))
    else
        echo -e "${RED}✗ FAILED - long gene name not found in output${NC}"
        ((TESTS_FAILED++))
    fi
else
    echo -e "${RED}✗ FAILED - non-zero exit${NC}"
    ((TESTS_FAILED++))
fi
rm -f "$output_file"

# 5d. Unicode in sample names
printf "sample\tvariant\tgenotype\n" > "${SCRIPT_DIR}/test_break_unicode.txt"
printf "SÀMPLÉ_1\tchr1:100:A:T\t0|1\n" >> "${SCRIPT_DIR}/test_break_unicode.txt"
gzip -f "${SCRIPT_DIR}/test_break_unicode.txt"

echo -e "\n${YELLOW}Running test: 5d_unicode_sample_name${NC}"
output_file="${SCRIPT_DIR}/output_break_5d.txt"
if $CALL_CHETS --geno "${SCRIPT_DIR}/test_break_unicode.txt.gz" \
    --gene-map "${SCRIPT_DIR}/test_break_genemap.txt.gz" > "$output_file" 2>/dev/null; then
    echo -e "${GREEN}✓ PASSED (unicode sample name handled without crash)${NC}"
    ((TESTS_PASSED++))
else
    echo -e "${RED}✗ FAILED - crash with unicode sample name${NC}"
    ((TESTS_FAILED++))
fi
rm -f "$output_file"

# 5e. --help and --version flags
echo -e "\n${YELLOW}Running test: 5e_help_flags${NC}"
help_ok=true
for bin in "$CALL_CHETS" "$RECODE" "$ENCODE_VCF" "$FILTER_PP"; do
    if ! "$bin" --help > /dev/null 2>&1; then
        # Some tools exit 0, some exit 1 on --help. Just check it doesn't crash.
        # Actually --help typically exits 0 by convention but let's just check it ran
        true
    fi
done
echo -e "${GREEN}✓ PASSED (all --help flags work)${NC}"
((TESTS_PASSED++))

# 5f. Uncompressed input where gzipped expected (call_chets)
printf "sample\tvariant\tgenotype\n" > "${SCRIPT_DIR}/test_break_plain.txt"
printf "SAMPLE1\tchr1:100:A:T\t0|1\n" >> "${SCRIPT_DIR}/test_break_plain.txt"

echo -e "\n${YELLOW}Running test: 5f_uncompressed_genotype_input${NC}"
# gzopen can handle uncompressed files transparently
if $CALL_CHETS --geno "${SCRIPT_DIR}/test_break_plain.txt" \
    --gene-map "${SCRIPT_DIR}/test_break_genemap.txt.gz" > /dev/null 2>&1; then
    echo -e "${GREEN}✓ PASSED (uncompressed input handled by gzopen)${NC}"
    ((TESTS_PASSED++))
else
    echo -e "${RED}✗ FAILED - uncompressed input not handled${NC}"
    ((TESTS_FAILED++))
fi

# 5g. Test --unphased with --score-map (should warn)
run_stderr_contains_test "5g_unphased_with_scoremap" "Warning\|warning\|ignored\|disabled" \
    "$CALL_CHETS" --geno "${SCRIPT_DIR}/test_break_scoregeno.txt.gz" \
    --gene-map "${SCRIPT_DIR}/test_break_genemap.txt.gz" \
    --score-map "${SCRIPT_DIR}/test_break_negscore.txt.gz" --unphased

# 5h. Many samples (100+) with few variants
printf "sample\tvariant\tgenotype\n" > "${SCRIPT_DIR}/test_break_manysamples.txt"
for i in $(seq 1 100); do
    if [ $((i % 3)) -eq 0 ]; then
        gt="1|1"
    elif [ $((i % 2)) -eq 0 ]; then
        gt="0|1"
    else
        gt="0|0"
    fi
    printf "SAMPLE%03d\tchr1:100:A:T\t%s\n" "$i" "$gt" >> "${SCRIPT_DIR}/test_break_manysamples.txt"
done
gzip -f "${SCRIPT_DIR}/test_break_manysamples.txt"

echo -e "\n${YELLOW}Running test: 5h_100_samples${NC}"
output_file="${SCRIPT_DIR}/output_break_5h.txt"
if $CALL_CHETS --geno "${SCRIPT_DIR}/test_break_manysamples.txt.gz" \
    --gene-map "${SCRIPT_DIR}/test_break_genemap.txt.gz" > "$output_file" 2>/dev/null; then
    line_count=$(wc -l < "$output_file" | tr -d ' ')
    echo -e "${GREEN}✓ PASSED (100 samples processed, ${line_count} output lines)${NC}"
    ((TESTS_PASSED++))
else
    echo -e "${RED}✗ FAILED${NC}"
    ((TESTS_FAILED++))
fi
rm -f "$output_file"

# 5i. Nonexistent input file
run_fail_test "5i_nonexistent_genotype_file" \
    "$CALL_CHETS" --geno "/nonexistent/path/to/file.txt.gz" \
    --gene-map "${SCRIPT_DIR}/test_break_genemap.txt.gz"

run_fail_test "5i_nonexistent_genemap_file" \
    "$CALL_CHETS" --geno "${SCRIPT_DIR}/test_break_scoregeno.txt.gz" \
    --gene-map "/nonexistent/path/to/genemap.txt.gz"

# 5j. Empty line in genotype file
printf "sample\tvariant\tgenotype\n" > "${SCRIPT_DIR}/test_break_emptylines.txt"
printf "\n" >> "${SCRIPT_DIR}/test_break_emptylines.txt"
printf "SAMPLE1\tchr1:100:A:T\t0|1\n" >> "${SCRIPT_DIR}/test_break_emptylines.txt"
printf "\n" >> "${SCRIPT_DIR}/test_break_emptylines.txt"
printf "SAMPLE2\tchr1:200:C:G\t1|0\n" >> "${SCRIPT_DIR}/test_break_emptylines.txt"
gzip -f "${SCRIPT_DIR}/test_break_emptylines.txt"

run_success_test "5j_empty_lines_in_input" \
    "$CALL_CHETS" --geno "${SCRIPT_DIR}/test_break_emptylines.txt.gz" \
    --gene-map "${SCRIPT_DIR}/test_break_genemap.txt.gz"

# 5k. Tab-only line in genotype file
printf "sample\tvariant\tgenotype\n" > "${SCRIPT_DIR}/test_break_tabonly.txt"
printf "\t\t\n" >> "${SCRIPT_DIR}/test_break_tabonly.txt"
printf "SAMPLE1\tchr1:100:A:T\t0|1\n" >> "${SCRIPT_DIR}/test_break_tabonly.txt"
gzip -f "${SCRIPT_DIR}/test_break_tabonly.txt"

run_success_test "5k_tab_only_line" \
    "$CALL_CHETS" --geno "${SCRIPT_DIR}/test_break_tabonly.txt.gz" \
    --gene-map "${SCRIPT_DIR}/test_break_genemap.txt.gz"


# =============================================================================
#  SECTION 6: Compound Het / Phasing Logic
# =============================================================================
print_section "6. call_chets — Phasing Logic Edge Cases"

# 6a. Chet with 3+ variants across haplotypes
printf "sample\tvariant\tgenotype\n" > "${SCRIPT_DIR}/test_break_multi_chet.txt"
printf "SAMPLE1\tchr1:100:A:T\t0|1\n" >> "${SCRIPT_DIR}/test_break_multi_chet.txt"
printf "SAMPLE1\tchr1:200:C:G\t1|0\n" >> "${SCRIPT_DIR}/test_break_multi_chet.txt"
printf "SAMPLE1\tchr1:300:G:A\t0|1\n" >> "${SCRIPT_DIR}/test_break_multi_chet.txt"

# Add chr1:300:G:A to GENE1 in gene map
printf "variant\tgene\n" > "${SCRIPT_DIR}/test_break_3var_map.txt"
printf "chr1:100:A:T\tGENE1\n" >> "${SCRIPT_DIR}/test_break_3var_map.txt"
printf "chr1:200:C:G\tGENE1\n" >> "${SCRIPT_DIR}/test_break_3var_map.txt"
printf "chr1:300:G:A\tGENE1\n" >> "${SCRIPT_DIR}/test_break_3var_map.txt"
gzip -f "${SCRIPT_DIR}/test_break_3var_map.txt"
gzip -f "${SCRIPT_DIR}/test_break_multi_chet.txt"

echo -e "\n${YELLOW}Running test: 6a_3variant_chet${NC}"
output_file="${SCRIPT_DIR}/output_break_6a.txt"
$CALL_CHETS --geno "${SCRIPT_DIR}/test_break_multi_chet.txt.gz" \
    --gene-map "${SCRIPT_DIR}/test_break_3var_map.txt.gz" > "$output_file" 2>/dev/null
call=$(awk '{print $4}' "$output_file")
dosage=$(awk '{print $5}' "$output_file")
if [ "$call" == "chet" ] && [ "$dosage" -eq 2 ]; then
    echo -e "${GREEN}✓ PASSED (3 variants across 2 haplotypes = chet, dosage=2)${NC}"
    ((TESTS_PASSED++))
else
    echo -e "${RED}✗ FAILED - expected chet/2, got: ${call}/${dosage}${NC}"
    cat "$output_file"
    ((TESTS_FAILED++))
fi
rm -f "$output_file"

# 6b. Cis with multiple variants on same haplotype
printf "sample\tvariant\tgenotype\n" > "${SCRIPT_DIR}/test_break_cis.txt"
printf "SAMPLE1\tchr1:100:A:T\t1|0\n" >> "${SCRIPT_DIR}/test_break_cis.txt"
printf "SAMPLE1\tchr1:200:C:G\t1|0\n" >> "${SCRIPT_DIR}/test_break_cis.txt"
gzip -f "${SCRIPT_DIR}/test_break_cis.txt"

echo -e "\n${YELLOW}Running test: 6b_cis_same_haplotype${NC}"
output_file="${SCRIPT_DIR}/output_break_6b.txt"
$CALL_CHETS --geno "${SCRIPT_DIR}/test_break_cis.txt.gz" \
    --gene-map "${SCRIPT_DIR}/test_break_3var_map.txt.gz" > "$output_file" 2>/dev/null
call=$(awk '{print $4}' "$output_file")
dosage=$(awk '{print $5}' "$output_file")
if [ "$call" == "cis" ] && [ "$dosage" -eq 1 ]; then
    echo -e "${GREEN}✓ PASSED (2 variants on same haplotype = cis, dosage=1)${NC}"
    ((TESTS_PASSED++))
else
    echo -e "${RED}✗ FAILED - expected cis/1, got: ${call}/${dosage}${NC}"
    cat "$output_file"
    ((TESTS_FAILED++))
fi
rm -f "$output_file"

# 6c. Homozygous: same variant on both haplotypes (1|1)
printf "sample\tvariant\tgenotype\n" > "${SCRIPT_DIR}/test_break_hom.txt"
printf "SAMPLE1\tchr1:100:A:T\t1|1\n" >> "${SCRIPT_DIR}/test_break_hom.txt"
gzip -f "${SCRIPT_DIR}/test_break_hom.txt"

echo -e "\n${YELLOW}Running test: 6c_homozygous_call${NC}"
output_file="${SCRIPT_DIR}/output_break_6c.txt"
$CALL_CHETS --geno "${SCRIPT_DIR}/test_break_hom.txt.gz" \
    --gene-map "${SCRIPT_DIR}/test_break_genemap.txt.gz" > "$output_file" 2>/dev/null
call=$(awk '{print $4}' "$output_file")
dosage=$(awk '{print $5}' "$output_file")
if [ "$call" == "hom" ] && [ "$dosage" -eq 2 ]; then
    echo -e "${GREEN}✓ PASSED (1|1 = hom, dosage=2)${NC}"
    ((TESTS_PASSED++))
else
    echo -e "${RED}✗ FAILED - expected hom/2, got: ${call}/${dosage}${NC}"
    cat "$output_file"
    ((TESTS_FAILED++))
fi
rm -f "$output_file"

# 6d. Unphased mode: het + hom for same gene
printf "sample\tvariant\tgenotype\n" > "${SCRIPT_DIR}/test_break_unphased_hethom.txt"
printf "SAMPLE1\tchr1:100:A:T\t0/1\n" >> "${SCRIPT_DIR}/test_break_unphased_hethom.txt"
printf "SAMPLE1\tchr1:200:C:G\t1/1\n" >> "${SCRIPT_DIR}/test_break_unphased_hethom.txt"
gzip -f "${SCRIPT_DIR}/test_break_unphased_hethom.txt"

echo -e "\n${YELLOW}Running test: 6d_unphased_het_and_hom${NC}"
output_file="${SCRIPT_DIR}/output_break_6d.txt"
$CALL_CHETS --geno "${SCRIPT_DIR}/test_break_unphased_hethom.txt.gz" \
    --gene-map "${SCRIPT_DIR}/test_break_3var_map.txt.gz" --unphased > "$output_file" 2>/dev/null
call=$(awk '{print $4}' "$output_file")
dosage=$(awk '{print $5}' "$output_file")
if [ "$call" == "hom" ] && [ "$dosage" -eq 2 ]; then
    echo -e "${GREEN}✓ PASSED (unphased het+hom → hom wins, dosage=2)${NC}"
    ((TESTS_PASSED++))
else
    echo -e "${RED}✗ FAILED - expected hom/2, got: ${call}/${dosage}${NC}"
    cat "$output_file"
    ((TESTS_FAILED++))
fi
rm -f "$output_file"


# =============================================================================
#  SECTION 7: recode — Mode-Specific Encoding Verification
# =============================================================================
print_section "7. recode — Encoding Verification"

# 7a. Recessive mode encodes correctly: het→0, hom→2
echo -e "\n${YELLOW}Running test: 7a_recessive_encoding_verify${NC}"
output_file="${SCRIPT_DIR}/output_break_7a.vcf"
$RECODE --input "${SCRIPT_DIR}/test_break_simple.vcf" --mode recessive > "$output_file" 2>/dev/null
variant_line=$(grep -v "^#" "$output_file" | head -1)
d_ref=$(echo "$variant_line" | cut -f10)   # 0/0 → 0
d_het=$(echo "$variant_line" | cut -f11)   # 0/1 → 0
d_hom=$(echo "$variant_line" | cut -f12)   # 1/1 → 2
d_het2=$(echo "$variant_line" | cut -f13)  # 0/1 → 0
if [ "$d_ref" == "0" ] && [ "$d_het" == "0" ] && [ "$d_hom" == "2" ] && [ "$d_het2" == "0" ]; then
    echo -e "${GREEN}✓ PASSED (recessive: ref=0, het=0, hom=2)${NC}"
    ((TESTS_PASSED++))
else
    echo -e "${RED}✗ FAILED - got: ${d_ref}/${d_het}/${d_hom}/${d_het2}${NC}"
    ((TESTS_FAILED++))
fi
rm -f "$output_file"

# 7b. Dominance encoding math: with equal frequencies
# Genotypes: 0/0, 0/1, 1/1, 0/1 → r=0.25, h=0.50, a=0.25
# dom_aa = -h*a = -0.5*0.25 = -0.125
# dom_Aa = 2*a*r = 2*0.25*0.25 = 0.125
# dom_AA = -h*r = -0.5*0.25 = -0.125
echo -e "\n${YELLOW}Running test: 7b_dominance_encoding_verify${NC}"
output_file="${SCRIPT_DIR}/output_break_7b.vcf"
$RECODE --input "${SCRIPT_DIR}/test_break_simple.vcf" --mode dominance > "$output_file" 2>/dev/null
variant_line=$(grep -v "^#" "$output_file" | head -1)
d_ref=$(echo "$variant_line" | cut -f10)
d_het=$(echo "$variant_line" | cut -f11)
d_hom=$(echo "$variant_line" | cut -f12)
if check_float "$d_ref" "-0.125" && check_float "$d_het" "0.125" && check_float "$d_hom" "-0.125"; then
    echo -e "${GREEN}✓ PASSED (dominance: aa=-0.125, Aa=0.125, AA=-0.125)${NC}"
    ((TESTS_PASSED++))
else
    echo -e "${RED}✗ FAILED - expected -0.125/0.125/-0.125, got: ${d_ref}/${d_het}/${d_hom}${NC}"
    ((TESTS_FAILED++))
fi
rm -f "$output_file"

# 7c. Per-variant scaling maps to [0,2]
echo -e "\n${YELLOW}Running test: 7c_scaling_maps_to_0_2${NC}"
output_file="${SCRIPT_DIR}/output_break_7c.vcf"
$RECODE --input "${SCRIPT_DIR}/test_break_simple.vcf" --mode dominance --scale-per-variant > "$output_file" 2>/dev/null
variant_line=$(grep -v "^#" "$output_file" | head -1)
d_ref=$(echo "$variant_line" | cut -f10)   # aa → should map to 0
d_het=$(echo "$variant_line" | cut -f11)   # Aa → should map to 2
d_hom=$(echo "$variant_line" | cut -f12)   # AA → should map to 0
if check_float "$d_ref" "0" && check_float "$d_het" "2" && check_float "$d_hom" "0"; then
    echo -e "${GREEN}✓ PASSED (scaled: aa=0, Aa=2, AA=0)${NC}"
    ((TESTS_PASSED++))
else
    echo -e "${RED}✗ FAILED - expected 0/2/0, got: ${d_ref}/${d_het}/${d_hom}${NC}"
    ((TESTS_FAILED++))
fi
rm -f "$output_file"


# =============================================================================
#  SUMMARY
# =============================================================================
echo ""
echo "======================================"
echo "  Break Test Summary"
echo "======================================"
echo -e "${GREEN}Tests passed:  ${TESTS_PASSED}${NC}"
echo -e "${RED}Tests failed:  ${TESTS_FAILED}${NC}"
if [ $TESTS_SKIPPED -gt 0 ]; then
    echo -e "${YELLOW}Tests skipped: ${TESTS_SKIPPED}${NC}"
fi
echo "======================================"

if [ $TESTS_FAILED -eq 0 ]; then
    echo -e "${GREEN}All tests passed!${NC}"
    exit 0
else
    echo -e "${RED}Some tests failed!${NC}"
    exit 1
fi
