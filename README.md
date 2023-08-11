# call_chets
Discard homozygous ref genotypes and output a list of variants and haplotypes.

todo: duplicate lines with alt/ref flipped

### Setup
* requires htslib


### Usage

1. List sites by sample with at least one alternate allele.
```
./get_non_ref_sites test/out.vcf 0.50 test/trio.vcf > phased_sites.txt
```

2. Group by gene and list samples with at least one alternate allele. 
```
./call_chets phased_sites.txt test/gene_map.txt > chets.txt
# head chets.txt
# HG00100	geneA	chet	2	20:1665:T:C-0|1;20:1869:A:T-1|0;20:2041:G:A-1|0;20:2220:G:A-1|0;20:2564:A:G-1|0
# HG00100	geneB	het	1	20:302:T:TA-1|0;20:3587:G:A-1|0;20:3936:A:G-1|0
# HG00100	geneC	het	1	20:828:T:C-1|0;20:834:G:A-1|0
# HG00101	geneA	hom+het	2	20:1665:T:C-1|0;20:2041:G:A-1|1;20:2220:G:A-0|1;20:2564:A:G-1|1
```

3. Convert into a VCF file.
```
./encode_vcf chets.txt tests/samples.txt
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG00100	HG00101	HG00102	HG00200	HG00201	HG00202	HG00300
#.	.	geneA	.	.	.	.	.	.	2	2	2	1	1	2	0
#.	.	geneB	.	.	.	.	.	.	1	2	2	1	1	2	0
#.	.	geneC	.	.	.	.	.	.	1	2	2	2	1	2	0

```

