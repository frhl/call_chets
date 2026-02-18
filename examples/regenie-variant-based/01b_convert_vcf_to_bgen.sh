#!/usr/bin/env bash
# Convert VCF files to BGEN format for REGENIE Step 2
# Converts both additive and dominance encodings

set -euo pipefail

input_dir="input"
output_dir="output"
mkdir -p "${output_dir}"

# Docker image
plink2_image="biocontainer/plink2:alpha2.3_jan2020"

echo "======================================"
echo "Convert VCF to BGEN"
echo "======================================"
echo ""

# Pull docker image if not present
echo "Checking for PLINK2 Docker image..."
if ! docker image inspect "${plink2_image}" >/dev/null 2>&1; then
    echo "Pulling Docker image ${plink2_image}..."
    docker pull "${plink2_image}"
fi

# Convert additive encoding (standard GT field)
echo "Converting additive encoding..."
vcf_add_src="../saige-set-based/input/simulated.vcf.gz"
vcf_add_tmp="${output_dir}/simulated.vcf.gz"
out_prefix_add="${output_dir}/simulated.additive"

# Copy VCF to output directory (resolve symlinks) for Docker access
echo "  Copying VCF from saige-set-based to output directory..."
cp -L "${vcf_add_src}" "${vcf_add_tmp}"
cp -L "${vcf_add_src}.csi" "${vcf_add_tmp}.csi" 2>/dev/null || true

docker run --rm \
  -v "$(pwd)/${output_dir}":/data \
  "${plink2_image}" \
  plink2 \
    --vcf /data/simulated.vcf.gz \
    --export bgen-1.3 'bits=16' ref-first \
    --out /data/simulated.additive

# Clean up temporary VCF
rm -f "${vcf_add_tmp}" "${vcf_add_tmp}.csi"

echo "  Created: ${out_prefix_add}.bgen"

# Convert dominance encoding (create 0/1/0 dosages from GT field)
echo "Converting dominance encoding..."
vcf_dom_src="../saige-set-based/input/simulated.vcf.gz"
vcf_dom_tmp="${output_dir}/simulated.dominance.vcf"
out_prefix_dom="${output_dir}/simulated.dominance"

# Copy VCF to output directory and create dominance dosages from GT
# Dominance encoding: 0/0 -> 0, 0/1 -> 1, 1/1 -> 0
echo "  Creating dominance dosages from GT field (heterozygotes = 1, homozygotes = 0)..."
bcftools view -H "${vcf_dom_src}" | \
  awk 'BEGIN {FS=OFS="\t"} {
    # Print header columns
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tDS", $1, $2, $3, $4, $5, $6, $7, $8
    # Process genotypes: het (0/1 or 1/0) -> 1, hom (0/0 or 1/1) -> 0
    for (i=10; i<=NF; i++) {
      split($i, gt, "/")
      if (gt[1] != gt[2]) {
        printf "\t1"  # Heterozygote
      } else {
        printf "\t0"  # Homozygote
      }
    }
    printf "\n"
  }' > "${vcf_dom_tmp}.body"

# Create VCF header (exclude existing FORMAT lines and #CHROM column header)
bcftools view -h "${vcf_dom_src}" | \
  grep -v "^##FORMAT=<ID=GT" | \
  sed '/^##FORMAT/d' | \
  grep -v "^#CHROM" > "${vcf_dom_tmp}.header"

# Add DS FORMAT line
echo '##FORMAT=<ID=DS,Number=1,Type=Float,Description="Dominance dosage: 1 for heterozygotes, 0 for homozygotes">' >> "${vcf_dom_tmp}.header"

# Add column header line
bcftools view -h "${vcf_dom_src}" | grep "^#CHROM" | \
  awk 'BEGIN {FS=OFS="\t"} {
    for (i=1; i<=9; i++) printf "%s\t", $i
    for (i=10; i<NF; i++) printf "%s\t", $i
    printf "%s\n", $NF
  }' >> "${vcf_dom_tmp}.header"

# Combine header and body
cat "${vcf_dom_tmp}.header" "${vcf_dom_tmp}.body" | \
  bgzip -c > "${vcf_dom_tmp}.gz"

bcftools index "${vcf_dom_tmp}.gz"

# Convert to BGEN using PLINK2
docker run --rm \
  -v "$(pwd)/${output_dir}":/data \
  "${plink2_image}" \
  plink2 \
    --vcf /data/simulated.dominance.vcf.gz dosage=DS \
    --import-dosage-certainty 1 \
    --hard-call-threshold 0 \
    --export bgen-1.3 'bits=16' ref-first \
    --out /data/simulated.dominance

# Clean up temporary files
rm -f "${vcf_dom_tmp}.header" "${vcf_dom_tmp}.body" "${vcf_dom_tmp}" "${vcf_dom_tmp}.gz" "${vcf_dom_tmp}.gz.csi"

echo "  Created: ${out_prefix_dom}.bgen"

echo ""
echo "======================================"
echo "BGEN Conversion Complete!"
echo "======================================"
echo "Output files:"
ls -lh "${output_dir}"/*.bgen 2>/dev/null || echo "  (No BGEN files created yet)"
echo ""
echo "Next step: Run 02_regenie_step2.sh for variant testing with both encodings"
