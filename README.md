# chet tools
This repository contains C++ scripts tailored for analyses with phased data, specifically involving variant call format (VCF) files.

### Dependencies:
* htslib:
* g++
* zlib1g-dev
* liblzma-dev
* libbz2-dev
* libcurl4-openssl-dev
* libhts-dev

### Installation
Scripts can be compiled using the provided Makefile. Ensure you have all the necessary dependencies installed or the docker.

### Usage

**Step 1**. We need to create a file of phased sites per gene per sample. First, we filter to AF<1%. If these are not excluded, the resulting file may be very large. We will use BCFtools to create a file in the right format:
```
bcftools view trio.vcf --max-af 0.01 -Ou | bcftools query -i'GT="alt"' -f'[%SAMPLE %CHROM:%POS:%REF:%ALT %GT\n]' | gzip > trio.phased_sites.txt.gz
```


**Step 2**. Call compound heterozygous, homozygous, heterozygous, cis variants. To do this, we group by gene and list samples with at least one alternate allele. This step can be
rinsed and repeated depending on the `--map` input. See below to find examples of mapping format.
```
./call_chets --geno trio.phased_sites.txt.gz --map gene_map.txt > test/trio.result.txt
```

**Step 3**. Convert into a VCF file using an additive encoding with the `--mode` argument. Here, genes with variants on a single haplotype will be encoded as 1 and genes with variants on both haplotypes will be encoded as 2. Use `--mode recessive` to only consider compound heterozygous and homozygous calls.
```
./encode_vcf --input trio.result.txt --samples test/samples.txt --mode additive | bgzip > trio.result.vcf.gz
zcat trio.result.vcf.gz | cut -f1-10 | head
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chr21>
##FORMAT=<ID=DS,Number=1,Type=Float,Description="Dosage">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1
chr21	1	ENSG00000141956	A	B	.	.	.	DS	0
chr21	2	ENSG00000141959	A	B	.	.	.	DS	0
chr21	3	ENSG00000142149	A	B	.	.	.	DS	0
chr21	4	ENSG00000142156	A	B	.	.	.	DS	0
chr21	5	ENSG00000142166	A	B	.	.	.	DS	0

```

### Mapping file
The mapping file can be obtained from running VEP on the variants of interest and extracting (using bash) the variant ID and the gene ID. Here is an example of a mapping file:
```
tbd
```





