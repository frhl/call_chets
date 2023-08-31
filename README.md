# chet tools
Tools to deal with phased data.

# How to run

0. Filter to AAF<1%. Ensure that these common variants are excluded before proceeding.
```
bcftools view --max-af 0.01 unfiltered.trio.vcf -o trio.vcf
```

1. List sites by sample with at least one alternate allele.
```
./get_non_ref_sites test/trio.vcf test/trio.phased_sites.txt.gz
```

2. Group by gene and list samples with at least one alternate allele. This step can be
rinsed and repeated depending on the gene_map input. 
```
./call_chets test/trio.phased_sites.txt.gz test/gene_map.txt > test/trio.result.txt
```

3. Convert into a VCF file.
```
./encode_vcf test/trio.result.txt test/samples.txt | bgzip > test/trio.result.vcf.gz
```

