#!/bin/bash

# Path to recode binary
RECODE=../bin/recode
INPUT=trio.vcf

# Ensure binary exists
if [ ! -f "$RECODE" ]; then
    echo "Error: $RECODE not found"
    exit 1
fi

echo "Running tests on $INPUT"

# Function to count variants in output
RECODE_CMD="$RECODE --input $INPUT --mode recessive"

# Function to count variants in output
count_variants() {
    $RECODE_CMD "$@" 2>/dev/null | grep -v "^#" | wc -l
}

# Baseline: 12 variants in total in trio.vcf? 
# Let's count them first.
TOTAL_VAR=$(grep -v "^#" $INPUT | wc -l)
echo "Total variants in input: $TOTAL_VAR"

# Test 1: min-aac
echo "Test 1: Filter by min-aac"
C1=$(count_variants --min-aac 100) # Should be 0
if [ "$C1" -eq 0 ]; then echo "PASS: min-aac 100 -> 0"; else echo "FAIL: min-aac 100 -> $C1"; fi

echo "Test 2: Filter by max-aac"
C2=$(count_variants --max-aac 0) # Should be 0 (assuming all have some alts)
if [ "$C2" -eq 0 ]; then echo "PASS: max-aac 0 -> 0"; else echo "FAIL: max-aac 0 -> $C2"; fi

# Let's look at specific variants
# Variant at 302: pos 302. AAC=8.
# Recode output includes INFO field, but recode doesn't output AAC/AF in INFO by default unless --all-info is on?
# Actually recode outputs AC and AN in INFO.
# Let's verify AC=8 for the first variant.
$RECODE_CMD --all-info 2>/dev/null | grep "302" | grep -q "AC=8"
if [ $? -eq 0 ]; then echo "PASS: Variant 302 has AC=8"; else echo "FAIL: Variant 302 AC check"; fi

# Filter min-aac 9 should exclude pos 302
$RECODE_CMD --min-aac 9 2>/dev/null | grep -q "302"
if [ $? -ne 0 ]; then echo "PASS: Variant 302 excluded with min-aac 9"; else echo "FAIL: Variant 302 included with min-aac 9"; fi

# Filter min-aac 8 should include pos 302
$RECODE_CMD --min-aac 8 2>/dev/null | grep -q "302"
if [ $? -eq 0 ]; then echo "PASS: Variant 302 included with min-aac 8"; else echo "FAIL: Variant 302 excluded with min-aac 8"; fi

# Test 3: MAC
# Total alleles = 12. AC=8. Ref=4. MAC=4.
# max-mac 3 should exclude
$RECODE_CMD --max-mac 3 2>/dev/null | grep -q "302"
if [ $? -ne 0 ]; then echo "PASS: Variant 302 excluded with max-mac 3"; else echo "FAIL: Variant 302 included with max-mac 3"; fi

# max-mac 4 should include
$RECODE_CMD --max-mac 4 2>/dev/null | grep -q "302"
if [ $? -eq 0 ]; then echo "PASS: Variant 302 included with max-mac 4"; else echo "FAIL: Variant 302 excluded with max-mac 4"; fi

# Test 4: MAF
# MAF = 4/12 = 0.3333...
# min-maf 0.34 should exclude (0.333 < 0.34)
$RECODE_CMD --min-maf 0.34 2>/dev/null | grep -q "302"
if [ $? -ne 0 ]; then echo "PASS: Variant 302 excluded with min-maf 0.34"; else echo "FAIL: Variant 302 included with min-maf 0.34"; fi

# min-maf 0.33 should include (0.333 > 0.33)
$RECODE_CMD --min-maf 0.33 2>/dev/null | grep -q "302"
if [ $? -eq 0 ]; then echo "PASS: Variant 302 included with min-maf 0.33"; else echo "FAIL: Variant 302 excluded with min-maf 0.33"; fi

# Test 5: AAF
# AAF = 8/12 = 0.6666...
# max-aaf 0.6 should exclude
$RECODE_CMD --max-aaf 0.6 2>/dev/null | grep -q "302"
if [ $? -ne 0 ]; then echo "PASS: Variant 302 excluded with max-aaf 0.6"; else echo "FAIL: Variant 302 included with max-aaf 0.6"; fi

echo "Done."
