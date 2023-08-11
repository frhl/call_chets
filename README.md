# call_chets
Discard homozygous ref genotypes and output a list of variants and haplotypes.

todo: duplicate lines with alt/ref flipped

### Setup
* requires htslib


### Usage

1. Get a list of gene-sample sites.
```
./get_phased_sites test/out.vcf 0.50 test/trio.vcf > phased_sites.txt
```

2. Get dosage by gene and call compound heterozygotes
```
./call_chets phased_sites.txt test/gene_map.txt > chets.txt
```

3. Convert gene dosages into VCF
```

```

