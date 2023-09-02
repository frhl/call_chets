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

**Step 0**. Filter to AAF<1%. Ensure that these common variants are excluded before proceeding.
```
cd test
bcftools view --max-af 0.01 unfiltered.trio.vcf -o trio.vcf
```

**Step 1**. List sites by sample with at least one alternate allele.
```
./get_non_ref_sites trio.vcf trio.phased_sites.txt.gz
```

Note, that step0+step1 can also be accomplished directly with the following command.
```
bcftools view trio.vcf --max-af 0.01 -Ou | bcftools query -i'GT="alt"' -f'[%SAMPLE %CHROM:%POS:%REF:%ALT %GT\n]' | gzip > trio.phased_sites.txt.gz
```




**Step 2**. Group by gene and list samples with at least one alternate allele. This step can be
rinsed and repeated depending on the gene_map input. 
```
./call_chets trio.phased_sites.txt.gz gene_map.txt > test/trio.result.txt
./call_chets trio.phased_sites.txt.gz gene_map.info.txt | grep chet
Sample1 chr21 ENSG00000274391 chet 2 chr21:10542449:G:A-splice_donor_variant;chr21:10605416:G:T-splice_acceptor_variant
Sample2 chr21 ENSG00000177398 chet 2 chr21:42076249:T:C-splice_donor_variant;chr21:42084101:C:T-stop_gained
```

**Step 3**. Convert into a VCF file.
```
./encode_vcf trio.result.txt test/samples.txt | bgzip > trio.result.vcf.gz
```

### Mapping file
The mapping file can be obtained from running VEP on the variants of interest and extracting (using bash) the variant ID and the gene ID.




