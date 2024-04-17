# Call compound heterozygous variants
This repository contains C++ scripts tailored for analyses with phased data, specifically involving variant call format (VCF) files.

## Installation
1. Ensure you have the necessary cpp libraries. Check the provided `Dockerfile` for the complete list.
2. Install [BCFtools](https://samtools.github.io/bcftools/howtos/install.html).
3. To compile the scripts, navigate to the repository's directory and run:
```
make
```

## Usage

### Step 1: Create a Phased Sites File
Generate a file containing phased sites per gene for each sample. First, filter to allele frequency less than 1% using BCFtools. Since a line will be produced for each variant-sample pair a lower MAF threshold will result in a smaller file.

```
bcftools view trio.vcf --max-af 0.01 -Ou | bcftools query -i'GT="alt"' -f'[%SAMPLE %CHROM:%POS:%REF:%ALT %GT\n]' | gzip > trio.phased_sites.txt.gz
```


### Step 2: Call bi/mono allelic variants
Call compound heterozygous and other variant types per group (e.g., genes). Repeat this step as needed based on the `--gene-map` argument. Note, that the --gene-map refers to a file with two columns for the variant IDs (preferably CHR:POS:REF:ALT) and the gene id (ensembl gene id).

```
./call_chets \
   --geno trio.phased_sites.txt.gz \
   --verbose \
   --gene-map gene_map.txt > test/trio.result.txt
```

For additional output information, provide a mapping file with variant, gene, and info columns for the `--info-map` argument. This approach adds an extra column to the output, displaying the info in a phase-aware manner. Info map contains three columns, variant, gene id and info string.

```
./call_chets \
   --geno trio.phased_sites.txt.gz \
   --gene-map gene_map.txt \
   --info-map info_map.txt  > test/trio.result.txt
```


Incorporate prior scores of pathogenicity and model the probability of complete knockout using the `--score-map`. Adjust the `--haplotype-collapse argument` to specify handling of variants on each haplotype (options: `product` (default), `burden`, or `max`). The process returns a score equivalent to the probability of both haplotypes being affected:

```
./call_chets \
   --geno trio.phased_sites.txt.gz \
   --gene-map gene_map.txt \
   --score-map score_map.txt \
   --info-map info_map.txt \
   --haplotype-collapse 'product' \
   --verbose | head
```


This is an example on output when all these arguments are specified and `| grep chet`. Sample, gene, configuration, dosage, probability of knockout based on alpha missense score and variant configuration (phase seperated by '|'). Here we have the sample, chromosome, gene id, variant configuration, gene collapsing rule, gene collapse score (here, probability of knockout of both haplotypes), and variants involved (seperated by "|" to indicate phase):
```
s1	chr21	ENSG00000215455	chet	2	g=product	0.0516	chr21:44539434:T:C:synonymous_variant|chr21:44539742:T:A:missense_variant
s2	chr21	ENSG00000205726	chet	2	g=product	0.48188	chr21:33829705:G:A:missense_variant|chr21:33875474:C:T:missense_variant
s3	chr21	ENSG00000215455	chet	2	g=product	1	chr21:44540070:ACAG:A:inframe_deletion|chr21:44539434:T:C:synonymous_variant
```
We can also take a look at cis (`| grep cis`) variants,  add the `--show-haplotype-scores` argument and see how the knockout scores are calculated for each haplotype:
```
5830593	chr21	ENSG00000142185	cis	1	g=product	0	h=product	0	0.896241	chr21:44414064:G:T:missense_variant;chr21:44425816:G:A:missense_variant
```

Two variants are found here (See below). We define the probability of haplotype knockout is defined as probability of haplotype knockout is `P(hap_ko)=1-((1-0.6391)*(1-0.7125))=0.896241`. If more variants are specified, these will also be included in the calculation. Other calculations (burden/min/max) can be invoked with the `--haplotype-collapse/-hc` and `--gene-collapse/-gc` arguments.

```
variant gene p(deleterious)
chr21:44414064:G:T	ENSG00000142185	0.6391
chr21:44425816:G:A	ENSG00000142185	0.7125
```


### Step 3: Create VCF
Convert the results into a VCF file with additive or recessive encoding. Genes with variants on either the paternal or maternal haplotype  will be encoded with a dosage of 1, while genes with variants on both haplotypes will be encoded as a dosage of 2. Use `--mode recessive` to only keep sites with both haplotypes affected.
```
./encode_vcf \
  --input trio.result.txt \
  --samples test/samples.txt \
  --mode additive \
  --min-ac 1 | bgzip > trio.result.vcf.gz
```

### Notes on testing for dominance deviation
TODO

## Haplotype and gene collapsing schemes:
Encode a model with probability of both haplotypes being knocked out:
```
./encode_vcf \
  --input trio.result.txt \
  --samples test/samples.txt \
  --mode additive \
  --min-ac 1 | bgzip > trio.result.vcf.gz
```

## Input file structure

### Gene Mapping file
Used with the `--gene-map` argument. Generate a mapping file by running VEP on your variants, then extract variant and gene IDs using a script. Example:
```
Variant Gene
20:1665:T:C geneA
20:1869:A:T geneA
20:2041:G:A geneA
20:2220:G:A geneA
```

### Info mapping file
Relates to  the `--info-map` argument. This file is generated similar to above. The 3rd column 'info' is keyed by the variant and gene.
```
varid	gene_id	csqs
chr21:10538674:C:A	ENSG00000274391	splice_region_variant
chr21:10538675:C:G	ENSG00000274391	splice_region_variant
chr21:10538680:G:C	ENSG00000274391	splice_acceptor_variant
chr21:10538734:GGT:G	ENSG00000274391	splice_donor_variant
```

### score mapping file
For the `--score-map argument`. The third column is indexed by variant and gene, with aggregation logic defined by the `--haplotype-collapse` argument:
```
varid	gene_id	am_pathogenicity
chr21:10541120:C:T	ENSG00000274391	0.0901
chr21:10541137:G:C	ENSG00000274391	0.1261
chr21:10541140:A:G	ENSG00000274391	0.0865
chr21:10541149:C:T	ENSG00000274391	0.087
chr21:10541158:A:C	ENSG00000274391	0.0729
..
```



## `run_call_chets.sh`: A wrapper that combines all the above (Needs to be tested with new version!)
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

