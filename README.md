# chettools

* requires htslib

### Usage

0. Filter to AAF<1%. Ensure that these common variants are excluded before proceeding.
```
bcftools view --max-af 0.01 unfiltered.trio.vcf -o trio.vcf
```


1. List sites by sample with at least one alternate allele.
```
./get_non_ref_sites test/trio.vcf test/trio.phased_sites.txt.gz
```

2. Group by gene and list samples with at least one alternate allele. 
```
./call_chets test/trio.phased_sites.txt.gz test/gene_map.txt > test/trio.result.txt
cat test/trio.result.txt | grep -Ew "(chet)|(cis)" | head -n2
# HG0010020geneAchet220:1665:T:C-0|1;20:1869:A:T-1|0;20:2041:G:A-1|0;20:2220:G:A-1|0;20:2564:A:G-1|0
# HG0010020geneAchet22020geneBcis120:302:T:TA-1|0;20:3936:A:G-1|0
```

3. Convert into a VCF file.
```
./encode_vcf test/trio.result.txt test/samples.txt | bgzip > test/trio.result.vcf.gz
bcftools view test/trio.result.vcf.gz

```

