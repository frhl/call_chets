# chettools

* requires htslib

### Usage

0. Filter to AAF<1%. Ensure that these common variants are excluded before proceeding.
```
bcftools view --max-af 0.01 unfiltered.trio.vcf -o trio.vcf
```


1. List sites by sample with at least one alternate allele.
```
./get_non_ref_sites test/trio.vcf phased_sites.txt.gz
```

2. Group by gene and list samples with at least one alternate allele. 
```
./call_chets phased_sites.txt.gz test/gene_map.txt | gzip -c > chets.txt.gz
# zcat chets.txt.gz | head
# HG00100	20	geneA	chet	2	20:1665:T:C-0|1;20:1869:A:T-1|0;20:2041:G:A-1|0;20:2220:G:A-1|0;20:2564:A:G-1|0
# HG00100	20	geneB	cis	1	20:302:T:TA-1|0;20:3587:G:A-1|0;20:3936:A:G-1|0
# HG00100	20	geneC	cis	1	20:828:T:C-1|0;20:834:G:A-1|0
# HG00101	20	geneA	hom+het	2	20:1665:T:C-1|0;20:2041:G:A-1|1;20:2220:G:A-0|1;20:2564:A:G-1|1
```

3. Convert into a VCF file.
```
./encode_vcf chets.txt test/samples.txt > my.vcf
bcftools view my.vcf
# ##fileformat=VCFv4.2
# ##FILTER=<ID=PASS,Description="All filters passed">
# ##contig=<ID=20>
# ##FORMAT=<ID=DS,Number=1,Type=Float,Description="Dosage">
# ##bcftools_viewVersion=1.17+htslib-1.17
# ##bcftools_viewCommand=view my.vcf; Date=Fri Aug 11 16:34:23 2023
# #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG00100	HG00101	HG00102	HG00200	HG00201	HG00202	HG00300
# 20	1	geneA	A	B	.	.	.	DS	2	2	2	1	1	2	0
# 20	2	geneB	A	B	.	.	.	DS	1	2	2	1	1	2	0
# 20	3	geneC	A	B	.	.	.	DS	1	2	2	2	1	2	0

```

