# Call compound heterozygous and homozygous sites by genes using phased genotypes
This repository contains C++ scripts tailored for analyses with phased data, specifically involving variant call format (VCF) files.

### Installation
1. Ensure you have the necessary cpp libraries. Check the provided `Dockerfile` for the complete list.
2. Install [BCFtools](https://samtools.github.io/bcftools/howtos/install.html).
3. To compile the scripts, navigate to the repository's directory and run:
```
make
```

### Usage

**Step 1: Create a Phased Sites File**. 
Generate a file containing phased sites per gene for each sample. Filter to allele frequency less than 1% using BCFtools.

```
bcftools view trio.vcf --max-af 0.01 -Ou | bcftools query -i'GT="alt"' -f'[%SAMPLE %CHROM:%POS:%REF:%ALT %GT\n]' | gzip > trio.phased_sites.txt.gz
```


**Step 2: Call bi/mono allelic variants**. 
This step calls compound heterozygous, homozygous, heterozygous, and cis variants by a group (typically genes). You can repeat this step based on the `--map` argument.

```
./call_chets \
   --geno trio.phased_sites.txt.gz \
   --map gene_map.txt > test/trio.result.txt
```

**Step 3: Create VCF**
Convert the results into a VCF file with additive or recessive encoding. Genes with variants on either the paternal or maternal haplotype  will be encoded with a dosage of 1, while genes with variants on both haplotypes will be encoded as a dosage of 2. Use `--mode recessive` to only keep sites with both haplotypes affected.
```
./encode_vcf \
  --input trio.result.txt \
  --samples test/samples.txt \
  --mode additive \
  --min-ac 1 | bgzip > trio.result.vcf.gz
```


### Mapping file
To create a mapping file, run VEP on your variants of interest. After that, use a bash script to extract the variant ID and the corresponding gene ID. Here's an example format of the mapping file:
```
Variant Gene Info
20:1665:T:C geneA splice_variant
20:1869:A:T geneA intron_variant
20:2041:G:A geneA intron_variant
20:2220:G:A geneA frameshift_variant
```

### Updates
- 13-09-23: A single variant can now map to multiple genes

### `run_call_chets.sh`: A wrapper that combines all the above 
That script combines several tools to a pipeline for CompHet calling. In particular, for a given chromosome:
1. Calls BCFtools (alternative to `get_non_ref_sites`) to extract non-ref genotypes for a phased BCF. 
2. Then uses `prepare_genemap.py` to prepare gene-variant maps according to each consequence (e.g. pLoF + damaging_missense).
3. Based on that, it calls `call_chets` to detect CompHet events from the output of step-1.
4. Finally, it runs `encode_vcf` to create a VCF based on each consequence and type of effects, e.g. additive or recessive.

As a last step, we could combine any chrom-based files to one genome-wide by running
```
for consq in pLoF pLoF_damaging damaging_missense synonymous; do
    for mode in additive recessive; do
        echo FIXTHIS.chr{1..22}.$mode.$consq.vcf.gz| tr ' ' '\n' > files.txt 
        bcftools concat -f files.txt -Oz -o FIXTHIS.chrALL.$mode.$consq.vcf.gz
    done
done
```






