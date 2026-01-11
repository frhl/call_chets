#!/usr/bin/env bash
# Extract markers for categorical variance ratio estimation
# Categories: 10<=MAC<20 and MAC>=20

set -euo pipefail

# Config
N_MARKERS=1000  # Markers per MAC category
PLINK="simulated"

# Calculate allele counts
docker run --rm -v "$(pwd)/input:/data" -w /data biocontainer/plink2:alpha2.3_jan2020 \
    plink2 --bfile ${PLINK} --freq counts --out ${PLINK} --allow-extra-chr

# Extract random markers from each MAC category
awk 'NR>1 {m1=$5; m2=2*$6-$5; if ((m1>=10 && m1<20) || (m2>=10 && m2<20)) print $2}' \
    input/${PLINK}.acount | sort -R | head -n ${N_MARKERS} > input/${PLINK}_cat1.tmp
awk 'NR>1 {m1=$5; m2=2*$6-$5; if (m1>=20 && m2>=20) print $2}' \
    input/${PLINK}.acount | sort -R | head -n ${N_MARKERS} > input/${PLINK}_cat2.tmp
cat input/${PLINK}_cat{1,2}.tmp > input/${PLINK}.vr_markers.list
rm input/${PLINK}_cat{1,2}.tmp

# Create PLINK file subset
docker run --rm -v "$(pwd)/input:/data" -w /data biocontainer/plink2:alpha2.3_jan2020 \
    plink2 --bfile ${PLINK} --extract ${PLINK}.vr_markers.list \
           --make-bed --out ${PLINK}_vr --allow-extra-chr

echo "✓ Created: input/${PLINK}_vr.{bed,bim,fam}"
