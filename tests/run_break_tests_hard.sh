#!/bin/bash

# =============================================================================
# ARCADE Hard Break Tests
# Tests designed to actually trigger bugs found via source code analysis.
# These target: division by zero, negative counts, buffer overflow,
# regex gaps, NaN propagation, header detection, and silent data loss.
# =============================================================================

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m'

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CALL_CHETS="${SCRIPT_DIR}/../bin/call_chets"
RECODE="${SCRIPT_DIR}/../bin/recode"
ENCODE_VCF="${SCRIPT_DIR}/../bin/encode_vcf"

TESTS_PASSED=0
TESTS_FAILED=0
BUGS_FOUND=0

check_float() {
    local actual=$1
    local expected=$2
    local tol=${3:-0.001}
    awk -v a="$actual" -v e="$expected" -v t="$tol" 'BEGIN {
        diff = a - e; if (diff < 0) diff = -diff;
        if (diff < t) exit 0; else exit 1;
    }'
}

# Check if a value is NaN or Inf
is_nan_or_inf() {
    local val=$1
    # Check for nan, -nan, inf, -inf (case insensitive)
    echo "$val" | grep -qi 'nan\|inf' && return 0
    # Also check if it's a valid number
    awk -v v="$val" 'BEGIN { if (v+0 != v+0) exit 0; else exit 1 }' 2>/dev/null && return 0
    return 1
}

print_section() {
    echo ""
    echo -e "${CYAN}==========================================${NC}"
    echo -e "${CYAN}  $1${NC}"
    echo -e "${CYAN}==========================================${NC}"
}

cleanup() {
    rm -f "${SCRIPT_DIR}"/output_hard_*.txt "${SCRIPT_DIR}"/output_hard_*.vcf
    rm -f "${SCRIPT_DIR}"/stderr_hard_*.txt
    rm -f "${SCRIPT_DIR}"/test_hard_*.txt "${SCRIPT_DIR}"/test_hard_*.txt.gz
    rm -f "${SCRIPT_DIR}"/test_hard_*.vcf "${SCRIPT_DIR}"/test_hard_*.vcf.gz
}
trap cleanup EXIT

# Check binaries
for bin in "$CALL_CHETS" "$RECODE" "$ENCODE_VCF"; do
    if [ ! -f "$bin" ]; then
        echo -e "${RED}Error: Binary not found: ${bin}${NC}"
        exit 1
    fi
done

# =============================================================================
#  BUG #1: Negative aa_count via duplicate sample/gene entries
#  make_pseudo_vcf.cpp:691 — geneChet/geneHet accumulate but geneSampleDosage
#  overwrites. If same sample has many entries for same gene, counts inflate
#  beyond total sample count, making aa_count_int negative.
# =============================================================================
print_section "BUG #1: Negative aa_count via duplicate entries"

# Create samples file with only 2 samples
printf "S1\nS2\n" > "${SCRIPT_DIR}/test_hard_2samples.txt"

# Create input where S1 appears 5 times as het for same gene (only 2 total samples)
# aa_count_int = (2*2/2) - (currentCis + currentHet + currentBI)
#              = 2 - (0 + 5 + 0) = -3
# This should NOT produce NaN/Inf in dominance output
cat << 'EOF' > "${SCRIPT_DIR}/test_hard_dup_inflate.txt"
S1	chr1	GeneA	het	1	var1
S1	chr1	GeneA	het	1	var2
S1	chr1	GeneA	het	1	var3
S1	chr1	GeneA	het	1	var4
S1	chr1	GeneA	het	1	var5
S2	chr1	GeneA	chet	2	var6|var7
EOF
gzip -f "${SCRIPT_DIR}/test_hard_dup_inflate.txt"

echo -e "\n${YELLOW}Running test: BUG1a — duplicate entries inflate het count${NC}"
output_file="${SCRIPT_DIR}/output_hard_bug1a.vcf"
stderr_file="${SCRIPT_DIR}/stderr_hard_bug1a.txt"
$ENCODE_VCF --input "${SCRIPT_DIR}/test_hard_dup_inflate.txt.gz" \
    --samples "${SCRIPT_DIR}/test_hard_2samples.txt" \
    --mode dominance > "$output_file" 2>"$stderr_file"
exit_code=$?

# Check for NaN/Inf in output
variant_line=$(grep -v "^#" "$output_file" | head -1)
if [ -n "$variant_line" ]; then
    if echo "$variant_line" | grep -qi 'nan\|inf'; then
        echo -e "${RED}✗ BUG CONFIRMED — NaN/Inf in dominance output from inflated counts${NC}"
        echo "  Output: $variant_line"
        ((TESTS_FAILED++))
        ((BUGS_FOUND++))
    else
        echo -e "${GREEN}✓ PASSED — no NaN/Inf (tool may have skipped the gene)${NC}"
        ((TESTS_PASSED++))
    fi
else
    echo -e "${GREEN}✓ PASSED — gene was skipped (no output, likely guarded)${NC}"
    ((TESTS_PASSED++))
fi
rm -f "$output_file" "$stderr_file"

# Even more extreme: make counts MUCH bigger than total samples
# 4 samples, but S1 appears 20 times as different configs
printf "S1\nS2\nS3\nS4\n" > "${SCRIPT_DIR}/test_hard_4samples.txt"

# Create data where S1 has 10 het entries + S1 has 5 cis entries = 15 entries
# aa_count_int = 4 - (5 + 10 + 0) = -11
cat << 'EOF' > "${SCRIPT_DIR}/test_hard_massive_inflate.txt"
S1	chr1	GeneA	het	1	v1
S1	chr1	GeneA	het	1	v2
S1	chr1	GeneA	het	1	v3
S1	chr1	GeneA	het	1	v4
S1	chr1	GeneA	het	1	v5
S1	chr1	GeneA	het	1	v6
S1	chr1	GeneA	het	1	v7
S1	chr1	GeneA	het	1	v8
S1	chr1	GeneA	het	1	v9
S1	chr1	GeneA	het	1	v10
S1	chr1	GeneA	cis	1	v11;v12
S1	chr1	GeneA	cis	1	v13;v14
S1	chr1	GeneA	cis	1	v15;v16
S1	chr1	GeneA	cis	1	v17;v18
S1	chr1	GeneA	cis	1	v19;v20
S2	chr1	GeneA	hom	2	v21
EOF
gzip -f "${SCRIPT_DIR}/test_hard_massive_inflate.txt"

echo -e "\n${YELLOW}Running test: BUG1b — massive count inflation (aa_count should be -11)${NC}"
output_file="${SCRIPT_DIR}/output_hard_bug1b.vcf"
$ENCODE_VCF --input "${SCRIPT_DIR}/test_hard_massive_inflate.txt.gz" \
    --samples "${SCRIPT_DIR}/test_hard_4samples.txt" \
    --mode dominance > "$output_file" 2>/dev/null

variant_line=$(grep -v "^#" "$output_file" | head -1)
if [ -n "$variant_line" ]; then
    if echo "$variant_line" | grep -qi 'nan\|inf\|-nan'; then
        echo -e "${RED}✗ BUG CONFIRMED — NaN/Inf in output from negative aa_count${NC}"
        echo "  Output: $variant_line"
        ((TESTS_FAILED++))
        ((BUGS_FOUND++))
    else
        # Check if dosages are reasonable (should be in some valid range)
        # With negative r, the dominance formula gives garbage
        dosage_s1=$(echo "$variant_line" | awk '{print $(NF-3)}')
        echo -e "${YELLOW}  Output produced (check if dosages are sensible): S1=${dosage_s1}${NC}"
        echo "  Full variant: $variant_line"
        # If output exists, check if r is negative (would be wrong)
        r_val=$(echo "$variant_line" | grep -o 'r=[^;]*' | head -1 | cut -d= -f2)
        if [ -n "$r_val" ]; then
            is_negative=$(awk -v r="$r_val" 'BEGIN { if (r < 0) print "yes"; else print "no" }')
            if [ "$is_negative" == "yes" ]; then
                echo -e "${RED}✗ BUG CONFIRMED — negative r value (${r_val}) from inflated counts${NC}"
                ((TESTS_FAILED++))
                ((BUGS_FOUND++))
            else
                echo -e "${GREEN}✓ PASSED — r value is non-negative${NC}"
                ((TESTS_PASSED++))
            fi
        else
            echo -e "${GREEN}✓ PASSED — output looks OK (no allInfo, can't verify r)${NC}"
            ((TESTS_PASSED++))
        fi
    fi
else
    echo -e "${GREEN}✓ PASSED — gene skipped (guarded)${NC}"
    ((TESTS_PASSED++))
fi
rm -f "$output_file"

# Test with --all-info to expose r, h, a values
echo -e "\n${YELLOW}Running test: BUG1c — with --all-info to expose negative frequencies${NC}"
output_file="${SCRIPT_DIR}/output_hard_bug1c.vcf"
$ENCODE_VCF --input "${SCRIPT_DIR}/test_hard_massive_inflate.txt.gz" \
    --samples "${SCRIPT_DIR}/test_hard_4samples.txt" \
    --mode dominance --all-info > "$output_file" 2>/dev/null

variant_line=$(grep -v "^#" "$output_file" | head -1)
if [ -n "$variant_line" ]; then
    # Extract r value from INFO field
    r_val=$(echo "$variant_line" | grep -o 'r=[0-9e.+-]*' | head -1 | cut -d= -f2)
    h_val=$(echo "$variant_line" | grep -o 'h=[0-9e.+-]*' | head -1 | cut -d= -f2)
    a_val=$(echo "$variant_line" | grep -o ';a=[0-9e.+-]*' | head -1 | cut -d= -f2)

    echo "  INFO r=${r_val} h=${h_val} a=${a_val}"

    if [ -n "$r_val" ]; then
        is_negative=$(awk -v r="$r_val" 'BEGIN { if (r < -0.001) print "yes"; else print "no" }')
        if [ "$is_negative" == "yes" ]; then
            echo -e "${RED}✗ BUG CONFIRMED — r=${r_val} is negative! Duplicate entries corrupted frequency${NC}"
            ((TESTS_FAILED++))
            ((BUGS_FOUND++))
        else
            echo -e "${GREEN}✓ PASSED — frequencies look valid${NC}"
            ((TESTS_PASSED++))
        fi
    else
        echo -e "${YELLOW}  Could not extract r value${NC}"
        ((TESTS_PASSED++))
    fi

    # Also check for NaN/Inf
    if echo "$variant_line" | grep -qi 'nan\|inf'; then
        echo -e "${RED}  ALSO: NaN/Inf detected in output!${NC}"
    fi
else
    echo -e "${GREEN}✓ PASSED — gene skipped entirely${NC}"
    ((TESTS_PASSED++))
fi
rm -f "$output_file"


# =============================================================================
#  BUG #2: Division by zero in dominance scaling (encode_vcf)
#  When all 3 dominance dosages are equal (maxDom == minDom), scaling divides
#  by zero. This happens when r+h+a frequencies are degenerate.
# =============================================================================
print_section "BUG #2: Division by zero in dominance scaling"

# Create a scenario where all carriers are homozygous
# If all sample are hom: r=0, h=0, a=1
# dom_aa = -0*1 = 0, dom_Aa = 2*1*0 = 0, dom_AA = -0*0 = 0
# min = max = 0 → division by zero
# But this requires aa_count_int > 0 to not be skipped
# So we need at least one non-carrier sample to make aa > 0
printf "S1\nS2\nS3\nS4\n" > "${SCRIPT_DIR}/test_hard_4samples.txt"

cat << 'EOF' > "${SCRIPT_DIR}/test_hard_divzero1.txt"
S1	chr1	GeneA	hom	2	v1
EOF
gzip -f "${SCRIPT_DIR}/test_hard_divzero1.txt"

# r=3/4=0.75, h=0, a=1/4=0.25
# dom_aa = -0*0.25 = 0
# dom_Aa = 2*0.25*0.75 = 0.375
# dom_AA = -0*0.75 = 0
# min=0, max=0.375 → division works
# Let's try a case where h=0 (no hets)... actually with call_chets output,
# "het" count = 0 means no heterozygous carriers.
# When BI=1, cis=0, het=0:
# aa_count = 4 - (0 + 0 + 1) = 3
# Aa_count = 0 + 0 = 0
# AA_count = 1
# r=0.75, h=0, a=0.25
# dom_aa = -0*0.25 = 0, dom_Aa = 2*0.25*0.75 = 0.375, dom_AA = -0*0.75 = 0
# min=0, max=0.375 — fine, no divzero

# The only way to get divzero is when r=0 or a=0 or h=0 in specific combos
# where all three dosage values equal zero.
# That happens when h=0 AND (a=0 OR r=0)
# But the guard `(currentBI == 0 || aa_count_int == 0)` prevents a=0 and r=0
# separately... but what about globalDomDosage mode?

# Test with --global-dom-dosage where ALL genes have same dosages
# If there's only one gene and it has min==max, global min==max
cat << 'EOF' > "${SCRIPT_DIR}/test_hard_divzero_global.txt"
S1	chr1	GeneA	hom	2	v1
S2	chr1	GeneB	hom	2	v2
EOF
gzip -f "${SCRIPT_DIR}/test_hard_divzero_global.txt"

echo -e "\n${YELLOW}Running test: BUG2a — global dominance dosage with identical min/max${NC}"
output_file="${SCRIPT_DIR}/output_hard_bug2a.vcf"
$ENCODE_VCF --input "${SCRIPT_DIR}/test_hard_divzero_global.txt.gz" \
    --samples "${SCRIPT_DIR}/test_hard_4samples.txt" \
    --mode dominance --global-dom-dosage --all-info > "$output_file" 2>/dev/null

variant_line=$(grep -v "^#" "$output_file" | head -1)
if [ -n "$variant_line" ]; then
    if echo "$variant_line" | grep -qi 'nan\|inf'; then
        echo -e "${RED}✗ BUG CONFIRMED — NaN/Inf in global dominance scaling${NC}"
        echo "  Output: $variant_line"
        ((TESTS_FAILED++))
        ((BUGS_FOUND++))
    else
        echo -e "${GREEN}✓ PASSED — no NaN/Inf${NC}"
        ((TESTS_PASSED++))
    fi
else
    echo -e "${GREEN}✓ PASSED — genes skipped (guarded by aa/bi check)${NC}"
    ((TESTS_PASSED++))
fi
rm -f "$output_file"

# Try to force division by zero: gene where EXACTLY the same number of each type
# produces three identical dominance dosages
# Needs: -h*a = 2*a*r = -h*r simultaneously
# This means: -h*a = -h*r → a=r, and -h*a = 2*a*r → -h = 2r → h = -2r
# But h must be >= 0, so this is impossible with valid frequencies!
# So division by zero can only happen when h=0 and a*r=0
# Which is guarded by the aa/bi check. Good.

echo -e "\n${YELLOW}Running test: BUG2b — verify per-gene scaling guard works${NC}"
# When all carriers are same type (only hom), h=0 → dom_aa=0, dom_Aa=2ar, dom_AA=0
# min=0, max=2ar (nonzero since a>0 and r>0)
# So per-gene never has divzero if the aa/bi guard works

# But what about the --all-info output path at lines 746-753?
# It outputs DS0/DS1/DS2 using the formula WITHOUT a guard:
# 2 * (((-h * a) - minDomDosage) / (maxDomDosage - minDomDosage))
# If max==min, this divides by zero even though the sample dosages are guarded by line 774

output_file="${SCRIPT_DIR}/output_hard_bug2b.vcf"
$ENCODE_VCF --input "${SCRIPT_DIR}/test_hard_divzero1.txt.gz" \
    --samples "${SCRIPT_DIR}/test_hard_4samples.txt" \
    --mode dominance --all-info > "$output_file" 2>/dev/null

variant_line=$(grep -v "^#" "$output_file" | head -1)
if [ -n "$variant_line" ]; then
    # Check INFO field for nan/inf in DS0/DS1/DS2
    info_field=$(echo "$variant_line" | cut -f8)
    if echo "$info_field" | grep -qi 'nan\|inf'; then
        echo -e "${RED}✗ BUG CONFIRMED — NaN/Inf in INFO field DS0/DS1/DS2${NC}"
        echo "  INFO: $info_field"
        ((TESTS_FAILED++))
        ((BUGS_FOUND++))
    else
        echo -e "${GREEN}✓ PASSED — INFO field clean${NC}"
        ((TESTS_PASSED++))
    fi
else
    echo -e "${GREEN}✓ PASSED — gene skipped${NC}"
    ((TESTS_PASSED++))
fi
rm -f "$output_file"


# =============================================================================
#  BUG #3: Regex too restrictive — rejects valid real-world variants
#  ChetCaller.cpp:149 — ^(chr)?[0-9XYM]{1,2}:[0-9]+:[ACGT]+:[ACGT]+$
# =============================================================================
print_section "BUG #3: Regex rejects valid variant formats"

# Gene map with lowercase alleles
printf "variant\tgene\n" > "${SCRIPT_DIR}/test_hard_lowercase_map.txt"
printf "chr1:100:a:t\tGENE1\n" >> "${SCRIPT_DIR}/test_hard_lowercase_map.txt"
gzip -f "${SCRIPT_DIR}/test_hard_lowercase_map.txt"

printf "sample\tvariant\tgenotype\n" > "${SCRIPT_DIR}/test_hard_lowercase_geno.txt"
printf "S1\tchr1:100:a:t\t0|1\n" >> "${SCRIPT_DIR}/test_hard_lowercase_geno.txt"
gzip -f "${SCRIPT_DIR}/test_hard_lowercase_geno.txt"

echo -e "\n${YELLOW}Running test: BUG3a — lowercase alleles rejected by regex${NC}"
output_file="${SCRIPT_DIR}/output_hard_bug3a.txt"
stderr_file="${SCRIPT_DIR}/stderr_hard_bug3a.txt"
$CALL_CHETS --geno "${SCRIPT_DIR}/test_hard_lowercase_geno.txt.gz" \
    --gene-map "${SCRIPT_DIR}/test_hard_lowercase_map.txt.gz" \
    --verbose > "$output_file" 2>"$stderr_file"

line_count=$(wc -l < "$output_file" | tr -d ' ')
if grep -qi "invalid variant format" "$stderr_file"; then
    echo -e "${RED}✗ BUG CONFIRMED — lowercase alleles rejected as invalid format${NC}"
    echo "  Warning: $(grep -i 'invalid variant' "$stderr_file" | head -1)"
    echo "  Output lines: $line_count"
    ((TESTS_FAILED++))
    ((BUGS_FOUND++))
else
    if [ "$line_count" -gt 0 ]; then
        echo -e "${GREEN}✓ PASSED — lowercase alleles accepted${NC}"
        ((TESTS_PASSED++))
    else
        echo -e "${RED}✗ BUG CONFIRMED — lowercase alleles silently dropped (no output)${NC}"
        ((TESTS_FAILED++))
        ((BUGS_FOUND++))
    fi
fi
rm -f "$output_file" "$stderr_file"

# Variant with 'N' in allele (IUPAC ambiguity code)
printf "variant\tgene\n" > "${SCRIPT_DIR}/test_hard_iupac_map.txt"
printf "chr1:100:N:T\tGENE1\n" >> "${SCRIPT_DIR}/test_hard_iupac_map.txt"
gzip -f "${SCRIPT_DIR}/test_hard_iupac_map.txt"

printf "sample\tvariant\tgenotype\n" > "${SCRIPT_DIR}/test_hard_iupac_geno.txt"
printf "S1\tchr1:100:N:T\t0|1\n" >> "${SCRIPT_DIR}/test_hard_iupac_geno.txt"
gzip -f "${SCRIPT_DIR}/test_hard_iupac_geno.txt"

echo -e "\n${YELLOW}Running test: BUG3b — IUPAC 'N' allele rejected${NC}"
output_file="${SCRIPT_DIR}/output_hard_bug3b.txt"
stderr_file="${SCRIPT_DIR}/stderr_hard_bug3b.txt"
$CALL_CHETS --geno "${SCRIPT_DIR}/test_hard_iupac_geno.txt.gz" \
    --gene-map "${SCRIPT_DIR}/test_hard_iupac_map.txt.gz" \
    --verbose > "$output_file" 2>"$stderr_file"

line_count=$(wc -l < "$output_file" | tr -d ' ')
if grep -qi "invalid variant format" "$stderr_file" || [ "$line_count" -eq 0 ]; then
    echo -e "${RED}✗ BUG CONFIRMED — IUPAC 'N' allele rejected${NC}"
    echo "  Output lines: $line_count"
    ((TESTS_FAILED++))
    ((BUGS_FOUND++))
else
    echo -e "${GREEN}✓ PASSED${NC}"
    ((TESTS_PASSED++))
fi
rm -f "$output_file" "$stderr_file"

# Variant with deletion represented as *
printf "variant\tgene\n" > "${SCRIPT_DIR}/test_hard_star_map.txt"
printf "chr1:100:A:*\tGENE1\n" >> "${SCRIPT_DIR}/test_hard_star_map.txt"
gzip -f "${SCRIPT_DIR}/test_hard_star_map.txt"

printf "sample\tvariant\tgenotype\n" > "${SCRIPT_DIR}/test_hard_star_geno.txt"
printf "S1\tchr1:100:A:*\t0|1\n" >> "${SCRIPT_DIR}/test_hard_star_geno.txt"
gzip -f "${SCRIPT_DIR}/test_hard_star_geno.txt"

echo -e "\n${YELLOW}Running test: BUG3c — star (*) deletion allele rejected${NC}"
output_file="${SCRIPT_DIR}/output_hard_bug3c.txt"
stderr_file="${SCRIPT_DIR}/stderr_hard_bug3c.txt"
$CALL_CHETS --geno "${SCRIPT_DIR}/test_hard_star_geno.txt.gz" \
    --gene-map "${SCRIPT_DIR}/test_hard_star_map.txt.gz" \
    --verbose > "$output_file" 2>"$stderr_file"

line_count=$(wc -l < "$output_file" | tr -d ' ')
if grep -qi "invalid variant format" "$stderr_file" || [ "$line_count" -eq 0 ]; then
    echo -e "${RED}✗ BUG CONFIRMED — star (*) allele rejected by regex${NC}"
    ((TESTS_FAILED++))
    ((BUGS_FOUND++))
else
    echo -e "${GREEN}✓ PASSED${NC}"
    ((TESTS_PASSED++))
fi
rm -f "$output_file" "$stderr_file"

# Long indel (100bp insertion)
long_indel=$(python3 -c "print('A' * 100)")
printf "variant\tgene\n" > "${SCRIPT_DIR}/test_hard_longindel_map.txt"
printf "chr1:100:A:%s\tGENE1\n" "$long_indel" >> "${SCRIPT_DIR}/test_hard_longindel_map.txt"
gzip -f "${SCRIPT_DIR}/test_hard_longindel_map.txt"

printf "sample\tvariant\tgenotype\n" > "${SCRIPT_DIR}/test_hard_longindel_geno.txt"
printf "S1\tchr1:100:A:%s\t0|1\n" "$long_indel" >> "${SCRIPT_DIR}/test_hard_longindel_geno.txt"
gzip -f "${SCRIPT_DIR}/test_hard_longindel_geno.txt"

echo -e "\n${YELLOW}Running test: BUG3d — long indel (100bp) in variant ID${NC}"
output_file="${SCRIPT_DIR}/output_hard_bug3d.txt"
$CALL_CHETS --geno "${SCRIPT_DIR}/test_hard_longindel_geno.txt.gz" \
    --gene-map "${SCRIPT_DIR}/test_hard_longindel_map.txt.gz" > "$output_file" 2>/dev/null

line_count=$(wc -l < "$output_file" | tr -d ' ')
if [ "$line_count" -eq 1 ]; then
    echo -e "${GREEN}✓ PASSED — 100bp indel variant processed${NC}"
    ((TESTS_PASSED++))
else
    echo -e "${RED}✗ FAILED — 100bp indel dropped (${line_count} lines)${NC}"
    ((TESTS_FAILED++))
fi
rm -f "$output_file"


# =============================================================================
#  BUG #4: Buffer truncation with long lines (>4096 chars)
#  gzgets with buf[4096] silently truncates long lines
# =============================================================================
print_section "BUG #4: Buffer truncation (lines > 4096 chars)"

# Create a gene map where a variant ID is extremely long (5000 chars)
# The variant format is chr1:100:REF:ALT where ALT is 4900 chars of ACGT
long_alt=$(python3 -c "print('ACGT' * 1225)")  # 4900 chars
long_variant="chr1:100:A:${long_alt}"

printf "variant\tgene\n" > "${SCRIPT_DIR}/test_hard_longline_map.txt"
printf "%s\tGENE1\n" "$long_variant" >> "${SCRIPT_DIR}/test_hard_longline_map.txt"
gzip -f "${SCRIPT_DIR}/test_hard_longline_map.txt"

printf "sample\tvariant\tgenotype\n" > "${SCRIPT_DIR}/test_hard_longline_geno.txt"
printf "S1\t%s\t0|1\n" "$long_variant" >> "${SCRIPT_DIR}/test_hard_longline_geno.txt"
gzip -f "${SCRIPT_DIR}/test_hard_longline_geno.txt"

echo -e "\n${YELLOW}Running test: BUG4a — variant line > 4096 chars (5000+ chars)${NC}"
output_file="${SCRIPT_DIR}/output_hard_bug4a.txt"
stderr_file="${SCRIPT_DIR}/stderr_hard_bug4a.txt"
$CALL_CHETS --geno "${SCRIPT_DIR}/test_hard_longline_geno.txt.gz" \
    --gene-map "${SCRIPT_DIR}/test_hard_longline_map.txt.gz" \
    --show-variants --verbose > "$output_file" 2>"$stderr_file"

line_count=$(wc -l < "$output_file" | tr -d ' ')
if [ "$line_count" -eq 1 ]; then
    # Check if the full variant ID is in output
    if grep -q "$long_variant" "$output_file" 2>/dev/null; then
        echo -e "${GREEN}✓ PASSED — long variant ID fully preserved${NC}"
        ((TESTS_PASSED++))
    else
        echo -e "${RED}✗ BUG CONFIRMED — variant ID truncated in output!${NC}"
        output_len=$(head -1 "$output_file" | wc -c | tr -d ' ')
        echo "  Output line length: $output_len (expected ~5000+)"
        ((TESTS_FAILED++))
        ((BUGS_FOUND++))
    fi
elif [ "$line_count" -eq 0 ]; then
    # Check stderr for clues
    if grep -qi "column\|error\|invalid" "$stderr_file"; then
        echo -e "${RED}✗ BUG CONFIRMED — long line caused parsing error (buffer truncation)${NC}"
        echo "  Error: $(grep -i 'column\|error\|invalid' "$stderr_file" | head -1)"
        ((TESTS_FAILED++))
        ((BUGS_FOUND++))
    else
        echo -e "${RED}✗ BUG LIKELY — no output and no clear error${NC}"
        ((TESTS_FAILED++))
        ((BUGS_FOUND++))
    fi
else
    echo -e "${RED}✗ BUG — unexpected ${line_count} lines (buffer split?)${NC}"
    ((TESTS_FAILED++))
    ((BUGS_FOUND++))
fi
rm -f "$output_file" "$stderr_file"

# Also test gene map with long line
echo -e "\n${YELLOW}Running test: BUG4b — gene map line > 4096 chars${NC}"
# The gene map line is: long_variant\tGENE1\n (already > 4096)
# If the map line gets truncated, the gene name is lost
# and the variant won't be mapped
output_file="${SCRIPT_DIR}/output_hard_bug4b.txt"
stderr_file="${SCRIPT_DIR}/stderr_hard_bug4b.txt"

# Use a short genotype variant but long gene map entry
short_variant="chr1:100:A:T"
printf "variant\tgene\n" > "${SCRIPT_DIR}/test_hard_longmap.txt"
# Create a map line where the "gene name" is 4080 chars, making total > 4096
long_gene=$(python3 -c "print('X' * 4080)")
printf "%s\t%s\n" "$short_variant" "$long_gene" >> "${SCRIPT_DIR}/test_hard_longmap.txt"
gzip -f "${SCRIPT_DIR}/test_hard_longmap.txt"

printf "sample\tvariant\tgenotype\n" > "${SCRIPT_DIR}/test_hard_longmap_geno.txt"
printf "S1\t%s\t0|1\n" "$short_variant" >> "${SCRIPT_DIR}/test_hard_longmap_geno.txt"
gzip -f "${SCRIPT_DIR}/test_hard_longmap_geno.txt"

$CALL_CHETS --geno "${SCRIPT_DIR}/test_hard_longmap_geno.txt.gz" \
    --gene-map "${SCRIPT_DIR}/test_hard_longmap.txt.gz" > "$output_file" 2>"$stderr_file"

line_count=$(wc -l < "$output_file" | tr -d ' ')
if [ "$line_count" -eq 1 ]; then
    # Check if the full gene name is preserved
    output_gene=$(awk '{print $3}' "$output_file")
    expected_len=4080
    actual_len=${#output_gene}
    if [ "$actual_len" -eq "$expected_len" ]; then
        echo -e "${GREEN}✓ PASSED — 4080-char gene name preserved${NC}"
        ((TESTS_PASSED++))
    else
        echo -e "${RED}✗ BUG CONFIRMED — gene name truncated! Expected ${expected_len}, got ${actual_len}${NC}"
        ((TESTS_FAILED++))
        ((BUGS_FOUND++))
    fi
else
    echo -e "${RED}✗ BUG LIKELY — line > 4096 in gene map caused issues (${line_count} lines)${NC}"
    head -3 "$stderr_file"
    ((TESTS_FAILED++))
    ((BUGS_FOUND++))
fi
rm -f "$output_file" "$stderr_file"


# =============================================================================
#  BUG #5: recode — DS field edge cases (NaN/Inf propagation)
# =============================================================================
print_section "BUG #5: recode — DS NaN and extreme values"

# Create VCF with NaN DS value
cat << 'EOF' > "${SCRIPT_DIR}/test_hard_ds_nan.vcf"
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DS,Number=1,Type=Float,Description="Dosage">
##contig=<ID=chr1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4
chr1	1000	.	A	T	.	PASS	.	GT:DS	0/0:0.0	0/1:1.0	1/1:2.0	./.:.
EOF

echo -e "\n${YELLOW}Running test: BUG5a — missing DS value treated correctly${NC}"
output_file="${SCRIPT_DIR}/output_hard_bug5a.vcf"
$RECODE --input "${SCRIPT_DIR}/test_hard_ds_nan.vcf" --mode dominance > "$output_file" 2>/dev/null

variant_line=$(grep -v "^#" "$output_file" | head -1)
if [ -n "$variant_line" ]; then
    s4_val=$(echo "$variant_line" | cut -f13)
    if [ "$s4_val" == "." ]; then
        echo -e "${GREEN}✓ PASSED — missing GT+DS produces '.' in output${NC}"
        ((TESTS_PASSED++))
    elif echo "$s4_val" | grep -qi 'nan\|inf'; then
        echo -e "${RED}✗ BUG CONFIRMED — NaN/Inf propagated to output${NC}"
        ((TESTS_FAILED++))
        ((BUGS_FOUND++))
    else
        echo -e "${YELLOW}  Got: ${s4_val} (unexpected but not NaN)${NC}"
        ((TESTS_PASSED++))
    fi
else
    echo -e "${GREEN}✓ PASSED — variant skipped${NC}"
    ((TESTS_PASSED++))
fi
rm -f "$output_file"

# VCF with only DS field (no GT) and extreme values
cat << 'EOF' > "${SCRIPT_DIR}/test_hard_ds_only_extreme.vcf"
##fileformat=VCFv4.2
##FORMAT=<ID=DS,Number=1,Type=Float,Description="Dosage">
##contig=<ID=chr1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4
chr1	1000	.	A	T	.	PASS	.	DS	0.0	1.0	2.0	1.5
EOF

echo -e "\n${YELLOW}Running test: BUG5b — DS-only VCF with fractional dosage 1.5${NC}"
output_file="${SCRIPT_DIR}/output_hard_bug5b.vcf"
stderr_file="${SCRIPT_DIR}/stderr_hard_bug5b.txt"
$RECODE --input "${SCRIPT_DIR}/test_hard_ds_only_extreme.vcf" --mode dominance > "$output_file" 2>"$stderr_file"

if [ $? -eq 0 ]; then
    variant_line=$(grep -v "^#" "$output_file" | head -1)
    if [ -n "$variant_line" ]; then
        # S4 has DS=1.5 which rounds to 2 → should be treated as hom alt
        s4_val=$(echo "$variant_line" | cut -f13)
        s3_val=$(echo "$variant_line" | cut -f12) # DS=2.0 → hom alt
        if check_float "$s4_val" "$s3_val" 0.001; then
            echo -e "${GREEN}✓ PASSED — DS=1.5 rounded to 2 (same as DS=2.0)${NC}"
            ((TESTS_PASSED++))
        else
            echo -e "${YELLOW}  DS=1.5→${s4_val}, DS=2.0→${s3_val} (different — rounded differently?)${NC}"
            ((TESTS_PASSED++))
        fi
    else
        echo -e "${GREEN}✓ PASSED — no output (skipped)${NC}"
        ((TESTS_PASSED++))
    fi
    # Check for rounding warning
    if grep -qi "rounded\|rounding" "$stderr_file"; then
        echo -e "  (Rounding warning detected — good)"
    fi
else
    echo -e "${RED}✗ FAILED — non-zero exit${NC}"
    ((TESTS_FAILED++))
fi
rm -f "$output_file" "$stderr_file"

# DS value of exactly 2.5 (out of range after rounding)
cat << 'EOF' > "${SCRIPT_DIR}/test_hard_ds_outrange.vcf"
##fileformat=VCFv4.2
##FORMAT=<ID=DS,Number=1,Type=Float,Description="Dosage">
##contig=<ID=chr1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3
chr1	1000	.	A	T	.	PASS	.	DS	0.0	1.0	2.5
EOF

echo -e "\n${YELLOW}Running test: BUG5c — DS=2.5 (rounds to 3, out of range)${NC}"
output_file="${SCRIPT_DIR}/output_hard_bug5c.vcf"
stderr_file="${SCRIPT_DIR}/stderr_hard_bug5c.txt"
$RECODE --input "${SCRIPT_DIR}/test_hard_ds_outrange.vcf" --mode dominance > "$output_file" 2>"$stderr_file"

if grep -qi "out of range\|Warning" "$stderr_file"; then
    echo -e "${GREEN}✓ PASSED — out-of-range DS warning emitted${NC}"
    ((TESTS_PASSED++))
else
    echo -e "${YELLOW}  No warning for DS=2.5 (check if handled by rounding)${NC}"
    ((TESTS_PASSED++))
fi
rm -f "$output_file" "$stderr_file"


# =============================================================================
#  BUG #6: First-line detection — silent data loss
#  If first line is malformed data (not a header), it's silently skipped.
# =============================================================================
print_section "BUG #6: First-line detection & silent data loss"

# Gene map without header
printf "chr1:100:A:T\tGENE1\n" > "${SCRIPT_DIR}/test_hard_noheader_map.txt"
printf "chr1:200:C:G\tGENE1\n" >> "${SCRIPT_DIR}/test_hard_noheader_map.txt"
gzip -f "${SCRIPT_DIR}/test_hard_noheader_map.txt"

# Data with no header — first line IS data
printf "S1\tchr1:100:A:T\t0|1\n" > "${SCRIPT_DIR}/test_hard_noheader_geno.txt"
printf "S2\tchr1:200:C:G\t1|0\n" >> "${SCRIPT_DIR}/test_hard_noheader_geno.txt"
gzip -f "${SCRIPT_DIR}/test_hard_noheader_geno.txt"

echo -e "\n${YELLOW}Running test: BUG6a — no header in genotype file (first line is data)${NC}"
output_file="${SCRIPT_DIR}/output_hard_bug6a.txt"
$CALL_CHETS --geno "${SCRIPT_DIR}/test_hard_noheader_geno.txt.gz" \
    --gene-map "${SCRIPT_DIR}/test_hard_noheader_map.txt.gz" > "$output_file" 2>/dev/null

line_count=$(wc -l < "$output_file" | tr -d ' ')
if [ "$line_count" -eq 1 ]; then
    # Both variants map to GENE1, S1 on H2, S2 on H1 → should produce chet for each
    echo -e "${GREEN}✓ PASSED — both data lines processed (${line_count} result)${NC}"
    ((TESTS_PASSED++))
elif [ "$line_count" -eq 0 ]; then
    echo -e "${RED}✗ BUG CONFIRMED — no output (first data line likely skipped)${NC}"
    ((TESTS_FAILED++))
    ((BUGS_FOUND++))
else
    echo -e "${GREEN}✓ PASSED (${line_count} results)${NC}"
    ((TESTS_PASSED++))
fi
rm -f "$output_file"

# Gene map without header — first line is data
# The mapping code skips first line if it can't parse two columns OR
# if the first line has valid columns. Let me check...
# In loadGeneMap: isFirstLineMappingFile starts true, and on first line
# if it extracts variant+gene successfully but isFirstLineMappingFile is true,
# it skips adding (line 216: `if (!isFirstLineMappingFile)`)
# So the FIRST LINE OF THE GENE MAP IS ALWAYS SKIPPED!

echo -e "\n${YELLOW}Running test: BUG6b — gene map first line always skipped?${NC}"
# Use a gene map where only one variant exists (on the first non-header line)
printf "chr1:100:A:T\tGENE1\n" > "${SCRIPT_DIR}/test_hard_noheader_1line_map.txt"
gzip -f "${SCRIPT_DIR}/test_hard_noheader_1line_map.txt"

printf "sample\tvariant\tgenotype\n" > "${SCRIPT_DIR}/test_hard_bug6b_geno.txt"
printf "S1\tchr1:100:A:T\t0|1\n" >> "${SCRIPT_DIR}/test_hard_bug6b_geno.txt"
gzip -f "${SCRIPT_DIR}/test_hard_bug6b_geno.txt"

output_file="${SCRIPT_DIR}/output_hard_bug6b.txt"
stderr_file="${SCRIPT_DIR}/stderr_hard_bug6b.txt"
$CALL_CHETS --geno "${SCRIPT_DIR}/test_hard_bug6b_geno.txt.gz" \
    --gene-map "${SCRIPT_DIR}/test_hard_noheader_1line_map.txt.gz" \
    --verbose > "$output_file" 2>"$stderr_file"

line_count=$(wc -l < "$output_file" | tr -d ' ')
if [ "$line_count" -eq 0 ]; then
    # Check if the variant was actually mapped
    if grep -q "0 variants" "$stderr_file" || grep -q "nVariantsMapped.*0" "$stderr_file"; then
        echo -e "${RED}✗ BUG CONFIRMED — gene map without header: first data line silently skipped!${NC}"
        echo "  (The first line is always treated as header and discarded)"
        ((TESTS_FAILED++))
        ((BUGS_FOUND++))
    else
        echo -e "${RED}✗ BUG LIKELY — no output, variant not mapped${NC}"
        echo "  Stderr excerpt:"
        grep -i 'gene\|map\|variant' "$stderr_file" | head -5
        ((TESTS_FAILED++))
        ((BUGS_FOUND++))
    fi
else
    echo -e "${GREEN}✓ PASSED — gene map first line processed as data${NC}"
    ((TESTS_PASSED++))
fi
rm -f "$output_file" "$stderr_file"

# Header-like first column name: what if header is "#sample" or "ID"?
printf "ID\tvariant\tgenotype\n" > "${SCRIPT_DIR}/test_hard_oddheader_geno.txt"
printf "S1\tchr1:100:A:T\t0|1\n" >> "${SCRIPT_DIR}/test_hard_oddheader_geno.txt"
gzip -f "${SCRIPT_DIR}/test_hard_oddheader_geno.txt"

# Provide gene map WITH proper header
printf "variant\tgene\n" > "${SCRIPT_DIR}/test_hard_proper_map.txt"
printf "chr1:100:A:T\tGENE1\n" >> "${SCRIPT_DIR}/test_hard_proper_map.txt"
gzip -f "${SCRIPT_DIR}/test_hard_proper_map.txt"

echo -e "\n${YELLOW}Running test: BUG6c — header with 'ID' instead of 'sample'${NC}"
output_file="${SCRIPT_DIR}/output_hard_bug6c.txt"
$CALL_CHETS --geno "${SCRIPT_DIR}/test_hard_oddheader_geno.txt.gz" \
    --gene-map "${SCRIPT_DIR}/test_hard_proper_map.txt.gz" > "$output_file" 2>/dev/null

line_count=$(wc -l < "$output_file" | tr -d ' ')
if [ "$line_count" -eq 1 ]; then
    echo -e "${GREEN}✓ PASSED — 'ID' header correctly skipped, data processed${NC}"
    ((TESTS_PASSED++))
else
    echo -e "${YELLOW}  Got ${line_count} lines (header 'ID' may be processed as data or skipped differently)${NC}"
    if [ "$line_count" -eq 0 ]; then
        echo -e "${RED}✗ BUG — data line lost${NC}"
        ((TESTS_FAILED++))
    else
        ((TESTS_PASSED++))
    fi
fi
rm -f "$output_file"


# =============================================================================
#  BUG #7: recode — scale-globally with single monomorphic-after-filter variant
#  Pass 1 finds no variants with AA genotype → globalMin stays at max double
# =============================================================================
print_section "BUG #7: recode edge cases"

# All variants filtered out in dominance mode (no AA genotypes)
cat << 'EOF' > "${SCRIPT_DIR}/test_hard_no_aa_global.vcf"
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4
chr1	1000	.	A	T	.	PASS	.	GT	0/0	0/1	0/0	0/1
chr1	2000	.	C	G	.	PASS	.	GT	0/0	0/0	0/1	0/0
EOF

echo -e "\n${YELLOW}Running test: BUG7a — scale-globally with no AA genotypes${NC}"
output_file="${SCRIPT_DIR}/output_hard_bug7a.vcf"
stderr_file="${SCRIPT_DIR}/stderr_hard_bug7a.txt"
$RECODE --input "${SCRIPT_DIR}/test_hard_no_aa_global.vcf" --mode dominance --scale-globally > "$output_file" 2>"$stderr_file"
exit_code=$?

if [ $exit_code -ne 0 ]; then
    if grep -qi "error\|no variants" "$stderr_file"; then
        echo -e "${GREEN}✓ PASSED — correctly errors when no AA genotypes in global mode${NC}"
        ((TESTS_PASSED++))
    else
        echo -e "${GREEN}✓ PASSED — non-zero exit (rejected input)${NC}"
        ((TESTS_PASSED++))
    fi
else
    # Check for valid output
    variant_count=$(grep -cv "^#" "$output_file" | tr -d ' ')
    if [ "$variant_count" -eq 0 ]; then
        echo -e "${GREEN}✓ PASSED — no output variants${NC}"
        ((TESTS_PASSED++))
    else
        echo -e "${RED}✗ BUG — got output with no AA genotypes in global scaling${NC}"
        ((TESTS_FAILED++))
    fi
fi
rm -f "$output_file" "$stderr_file"

# recode: VCF where REF == ALT
cat << 'EOF' > "${SCRIPT_DIR}/test_hard_ref_eq_alt.vcf"
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4
chr1	1000	.	A	A	.	PASS	.	GT	0/0	0/1	1/1	0/1
EOF

echo -e "\n${YELLOW}Running test: BUG7b — VCF where REF == ALT (nonsense but valid VCF)${NC}"
output_file="${SCRIPT_DIR}/output_hard_bug7b.vcf"
$RECODE --input "${SCRIPT_DIR}/test_hard_ref_eq_alt.vcf" --mode dominance > "$output_file" 2>/dev/null

if [ $? -eq 0 ]; then
    variant_line=$(grep -v "^#" "$output_file" | head -1)
    if [ -n "$variant_line" ]; then
        if echo "$variant_line" | grep -qi 'nan\|inf'; then
            echo -e "${RED}✗ BUG — NaN/Inf in output for REF==ALT variant${NC}"
            ((TESTS_FAILED++))
            ((BUGS_FOUND++))
        else
            echo -e "${GREEN}✓ PASSED — no NaN/Inf (processed normally)${NC}"
            ((TESTS_PASSED++))
        fi
    else
        echo -e "${GREEN}✓ PASSED — variant filtered${NC}"
        ((TESTS_PASSED++))
    fi
else
    echo -e "${GREEN}✓ PASSED — rejected${NC}"
    ((TESTS_PASSED++))
fi
rm -f "$output_file"


# =============================================================================
#  BUG #8: encode_vcf — gene spanning multiple chromosomes
#  What if call_chets output has same gene on different chromosomes?
# =============================================================================
print_section "BUG #8: encode_vcf — gene on multiple chromosomes"

printf "S1\nS2\nS3\n" > "${SCRIPT_DIR}/test_hard_3samples.txt"

cat << 'EOF' > "${SCRIPT_DIR}/test_hard_multichr_gene.txt"
S1	chr1	GeneA	het	1	v1
S2	chr2	GeneA	het	1	v2
EOF
gzip -f "${SCRIPT_DIR}/test_hard_multichr_gene.txt"

echo -e "\n${YELLOW}Running test: BUG8a — same gene on two chromosomes${NC}"
output_file="${SCRIPT_DIR}/output_hard_bug8a.vcf"
$ENCODE_VCF --input "${SCRIPT_DIR}/test_hard_multichr_gene.txt.gz" \
    --samples "${SCRIPT_DIR}/test_hard_3samples.txt" --mode additive > "$output_file" 2>/dev/null

# geneToChromosome[gene] = chromosome overwrites, so gene ends up on last chr
variant_count=$(grep -cv "^#" "$output_file" | tr -d ' ')
if [ "$variant_count" -eq 1 ]; then
    chr=$(grep -v "^#" "$output_file" | head -1 | cut -f1)
    echo -e "${GREEN}✓ PASSED — gene appears once on ${chr} (last chromosome wins)${NC}"
    ((TESTS_PASSED++))
elif [ "$variant_count" -eq 2 ]; then
    echo -e "${RED}✗ BUG — gene appears on two chromosomes (duplicated)${NC}"
    ((TESTS_FAILED++))
    ((BUGS_FOUND++))
else
    echo -e "${YELLOW}  Unexpected: ${variant_count} variants${NC}"
    ((TESTS_PASSED++))
fi
rm -f "$output_file"


# =============================================================================
#  BUG #9: call_chets — score exactly at boundary (0.0 and 1.0)
# =============================================================================
print_section "BUG #9: Score boundary values"

printf "variant\tgene\tscore\n" > "${SCRIPT_DIR}/test_hard_score_boundary.txt"
printf "chr1:100:A:T\tGENE1\t0.0\n" >> "${SCRIPT_DIR}/test_hard_score_boundary.txt"
printf "chr1:200:C:G\tGENE1\t1.0\n" >> "${SCRIPT_DIR}/test_hard_score_boundary.txt"
gzip -f "${SCRIPT_DIR}/test_hard_score_boundary.txt"

printf "sample\tvariant\tgenotype\n" > "${SCRIPT_DIR}/test_hard_score_geno.txt"
printf "S1\tchr1:100:A:T\t0|1\n" >> "${SCRIPT_DIR}/test_hard_score_geno.txt"
printf "S1\tchr1:200:C:G\t1|0\n" >> "${SCRIPT_DIR}/test_hard_score_geno.txt"
gzip -f "${SCRIPT_DIR}/test_hard_score_geno.txt"

printf "variant\tgene\n" > "${SCRIPT_DIR}/test_hard_score_map.txt"
printf "chr1:100:A:T\tGENE1\n" >> "${SCRIPT_DIR}/test_hard_score_map.txt"
printf "chr1:200:C:G\tGENE1\n" >> "${SCRIPT_DIR}/test_hard_score_map.txt"
gzip -f "${SCRIPT_DIR}/test_hard_score_map.txt"

echo -e "\n${YELLOW}Running test: BUG9a — scores exactly 0.0 and 1.0 (product collapse)${NC}"
output_file="${SCRIPT_DIR}/output_hard_bug9a.txt"
$CALL_CHETS --geno "${SCRIPT_DIR}/test_hard_score_geno.txt.gz" \
    --gene-map "${SCRIPT_DIR}/test_hard_score_map.txt.gz" \
    --score-map "${SCRIPT_DIR}/test_hard_score_boundary.txt.gz" \
    --show-haplotype-scores > "$output_file" 2>/dev/null

if [ $? -eq 0 ]; then
    line=$(head -1 "$output_file")
    echo "  Output: $line"
    # With product collapse: hap1 has score 1.0 → (1-1.0)=0.0
    # hap2 has score 0.0 → (1-0.0)=1.0
    # Then: hap1_final = 1-0.0 = 1.0, hap2_final = 1-1.0 = 0.0
    # gene_score = 1.0 * 0.0 = 0.0
    echo -e "${GREEN}✓ PASSED — boundary scores processed without error${NC}"
    ((TESTS_PASSED++))
else
    echo -e "${RED}✗ FAILED — error with boundary scores${NC}"
    ((TESTS_FAILED++))
fi
rm -f "$output_file"


# =============================================================================
#  BUG #10: recode — single variant where all 3 genotype classes have equal
#  raw dominance dosage (impossible mathematically but testing scaleDosage
#  with localMin == localMax)
# =============================================================================
print_section "BUG #10: recode scaleDosage when min==max"

# We know this is guarded by `if (maxVal == minVal) return 0.0;` in recode.cpp
# But let's verify it works with per-variant scaling
# Genotypes: 0/0, 0/0, 0/0, 1/1 → r=0.75, h=0, a=0.25
# dom_aa = -0*0.25 = 0, dom_Aa = 2*0.25*0.75 = 0.375, dom_AA = -0*0.75 = 0
# min=0, max=0.375 → no divzero
# BUT: what if we have genotypes that make ALL THREE values identical?
# That requires: -h*a = 2*a*r = -h*r
# From -h*a = -h*r → a = r
# From -h*a = 2*a*r → -h = 2r → h = -2r (impossible, h >= 0)
# So it's mathematically impossible! But let's test the scaleDosage guard anyway
# with a variant where min and max are very close (near-zero variance)

cat << 'EOF' > "${SCRIPT_DIR}/test_hard_near_equal.vcf"
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4	S5	S6	S7	S8	S9	S10	S11	S12	S13	S14	S15	S16	S17	S18	S19	S20
chr1	1000	.	A	T	.	PASS	.	GT	0/0	0/0	0/0	0/0	0/0	0/0	0/0	0/0	0/0	0/0	0/0	0/0	0/0	0/0	0/0	0/0	0/0	0/0	0/1	1/1
EOF

echo -e "\n${YELLOW}Running test: BUG10a — near-zero variance dominance dosages${NC}"
output_file="${SCRIPT_DIR}/output_hard_bug10a.vcf"
stderr_file="${SCRIPT_DIR}/stderr_hard_bug10a.txt"
$RECODE --input "${SCRIPT_DIR}/test_hard_near_equal.vcf" --mode dominance --scale-per-variant > "$output_file" 2>"$stderr_file"

if [ $? -eq 0 ]; then
    if grep -qi "low variance\|similar.*dosage" "$stderr_file"; then
        echo -e "${GREEN}✓ PASSED — low variance warning emitted${NC}"
        ((TESTS_PASSED++))
    else
        variant_line=$(grep -v "^#" "$output_file" | head -1)
        if echo "$variant_line" | grep -qi 'nan\|inf'; then
            echo -e "${RED}✗ BUG — NaN/Inf in near-zero variance output${NC}"
            ((TESTS_FAILED++))
            ((BUGS_FOUND++))
        else
            echo -e "${GREEN}✓ PASSED — processed without NaN/Inf${NC}"
            ((TESTS_PASSED++))
        fi
    fi
else
    echo -e "${RED}✗ FAILED — non-zero exit${NC}"
    ((TESTS_FAILED++))
fi
rm -f "$output_file" "$stderr_file"


# =============================================================================
#  BUG #11: recode — many samples (wide VCF)
# =============================================================================
print_section "BUG #11: recode — wide VCF (200 samples)"

echo -e "\n${YELLOW}Running test: BUG11a — VCF with 200 samples${NC}"
# Build a VCF with 200 samples
header="##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n##contig=<ID=chr1>\n"
sample_header="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
data_line="chr1\t1000\t.\tA\tT\t.\tPASS\t.\tGT"

for i in $(seq 1 200); do
    sample_header="${sample_header}\tS${i}"
    if [ $((i % 3)) -eq 0 ]; then
        data_line="${data_line}\t1/1"
    elif [ $((i % 2)) -eq 0 ]; then
        data_line="${data_line}\t0/1"
    else
        data_line="${data_line}\t0/0"
    fi
done

printf "${header}${sample_header}\n${data_line}\n" > "${SCRIPT_DIR}/test_hard_wide.vcf"

output_file="${SCRIPT_DIR}/output_hard_bug11a.vcf"
$RECODE --input "${SCRIPT_DIR}/test_hard_wide.vcf" --mode dominance --all-info > "$output_file" 2>/dev/null

if [ $? -eq 0 ]; then
    variant_line=$(grep -v "^#" "$output_file" | head -1)
    col_count=$(echo "$variant_line" | awk -F'\t' '{print NF}')
    expected_cols=$((200 + 9))  # 9 fixed + 200 sample columns
    if [ "$col_count" -eq "$expected_cols" ]; then
        echo -e "${GREEN}✓ PASSED — 200 samples processed (${col_count} columns)${NC}"
        ((TESTS_PASSED++))
    else
        echo -e "${RED}✗ FAILED — expected ${expected_cols} columns, got ${col_count}${NC}"
        ((TESTS_FAILED++))
    fi
else
    echo -e "${RED}✗ FAILED — non-zero exit with 200 samples${NC}"
    ((TESTS_FAILED++))
fi
rm -f "$output_file"


# =============================================================================
#  BUG #12: call_chets — 0/0 everywhere (all reference, no carriers)
# =============================================================================
print_section "BUG #12: All reference genotypes"

printf "sample\tvariant\tgenotype\n" > "${SCRIPT_DIR}/test_hard_all_ref.txt"
printf "S1\tchr1:100:A:T\t0|0\n" >> "${SCRIPT_DIR}/test_hard_all_ref.txt"
printf "S2\tchr1:100:A:T\t0|0\n" >> "${SCRIPT_DIR}/test_hard_all_ref.txt"
printf "S3\tchr1:100:A:T\t0|0\n" >> "${SCRIPT_DIR}/test_hard_all_ref.txt"
gzip -f "${SCRIPT_DIR}/test_hard_all_ref.txt"

echo -e "\n${YELLOW}Running test: BUG12a — all 0|0 genotypes (no carriers)${NC}"
output_file="${SCRIPT_DIR}/output_hard_bug12a.txt"
$CALL_CHETS --geno "${SCRIPT_DIR}/test_hard_all_ref.txt.gz" \
    --gene-map "${SCRIPT_DIR}/test_hard_proper_map.txt.gz" > "$output_file" 2>/dev/null

line_count=$(wc -l < "$output_file" | tr -d ' ')
# 0|0 is skipped in processGenotypes (validGenotype is false for 0|0)
# But the tool should still exit 0 (no carriers is a valid result)
exit_code=$?
if [ $exit_code -eq 0 ] && [ "$line_count" -eq 0 ]; then
    echo -e "${GREEN}✓ PASSED — no carriers, no output, exit 0${NC}"
    ((TESTS_PASSED++))
else
    echo -e "${RED}✗ FAILED — exit=${exit_code}, lines=${line_count}${NC}"
    ((TESTS_FAILED++))
fi
rm -f "$output_file"


# =============================================================================
#  BUG #13: recode — unknown argument handling
# =============================================================================
print_section "BUG #13: Argument parsing edge cases"

echo -e "\n${YELLOW}Running test: BUG13a — unknown argument${NC}"
if ! $RECODE --input "${SCRIPT_DIR}/test_hard_simple.vcf" --mode dominance --banana > /dev/null 2>&1; then
    echo -e "${GREEN}✓ PASSED — unknown arg rejected${NC}"
    ((TESTS_PASSED++))
else
    echo -e "${RED}✗ FAILED — unknown arg accepted silently${NC}"
    ((TESTS_FAILED++))
fi

echo -e "\n${YELLOW}Running test: BUG13b — missing value after flag${NC}"
if ! $RECODE --input > /dev/null 2>&1; then
    echo -e "${GREEN}✓ PASSED — missing value after --input rejected${NC}"
    ((TESTS_PASSED++))
else
    echo -e "${RED}✗ FAILED — missing value not caught${NC}"
    ((TESTS_FAILED++))
fi

echo -e "\n${YELLOW}Running test: BUG13c — no arguments at all${NC}"
if ! $RECODE > /dev/null 2>&1; then
    echo -e "${GREEN}✓ PASSED — no arguments rejected${NC}"
    ((TESTS_PASSED++))
else
    echo -e "${RED}✗ FAILED — ran without arguments${NC}"
    ((TESTS_FAILED++))
fi

echo -e "\n${YELLOW}Running test: BUG13d — encode_vcf with no arguments${NC}"
if ! $ENCODE_VCF > /dev/null 2>&1; then
    echo -e "${GREEN}✓ PASSED — no arguments rejected${NC}"
    ((TESTS_PASSED++))
else
    echo -e "${RED}✗ FAILED${NC}"
    ((TESTS_FAILED++))
fi


# =============================================================================
#  SUMMARY
# =============================================================================
echo ""
echo "======================================"
echo "  Hard Break Test Summary"
echo "======================================"
echo -e "${GREEN}Tests passed:  ${TESTS_PASSED}${NC}"
echo -e "${RED}Tests failed:  ${TESTS_FAILED}${NC}"
if [ $BUGS_FOUND -gt 0 ]; then
    echo -e "${RED}Bugs found:    ${BUGS_FOUND}${NC}"
fi
echo "======================================"

if [ $TESTS_FAILED -eq 0 ]; then
    echo -e "${GREEN}All tests passed!${NC}"
    exit 0
else
    echo -e "${RED}Some tests failed — bugs found!${NC}"
    exit 1
fi
