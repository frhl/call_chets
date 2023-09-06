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
Generate a file containing phased sites per gene for each sample. Filter to have allele frequency (AF) less than 1% using BCFtools. Not filtering these might result in large output files.

```
bcftools view trio.vcf --max-af 0.01 -Ou | bcftools query -i'GT="alt"' -f'[%SAMPLE %CHROM:%POS:%REF:%ALT %GT\n]' | gzip > trio.phased_sites.txt.gz
```


**Step 2: Call bi and mono-allelic variants**. 
This step calls compound heterozygous, homozygous, heterozygous, and cis variants. Variants are grouped by gene and list samples with at least one alternate allele. You can repeat this step based on the `--map` argument.

```
./call_chets --geno trio.phased_sites.txt.gz --map gene_map.txt > test/trio.result.txt
```

**Step 3: Create VCF**
Here, you'll convert the results into a VCF file with additive encoding. Genes with variants on one haplotype will be encoded as '1', while those on both haplotypes will be '2'. To only consider compound heterozygous and homozygous calls, use the `--mode recessive` argument.
```
./encode_vcf --input trio.result.txt --samples test/samples.txt --mode additive | bgzip > trio.result.vcf.gz && bcftools index trio.result.vcf.gz
```
Remember, if combining multiple chromosomes, ensure they're in the correct order. Otherwise, tools like plink/bcftools might return errors during downstream processes.


### Mapping file
To create a mapping file, run VEP on your variants of interest. After that, use a bash script to extract the variant ID and the corresponding gene ID. Here's an example format of the mapping file:
```
tbd
```





