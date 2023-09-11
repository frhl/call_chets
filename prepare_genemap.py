#!/usr/bin/env python
# coding: utf-8
"""
A script converts the gene annotation file to a variant-gene index to use in `call_hets`, based on a specific consequence.

Run as
`python prepare_genemap.py -a brava_GNH_44k_chr22.txt -c pLoF X Y -o gene_map.pLoF_X_Y.txt`

###################
GK - Sept 1st, 2023

TODO: what if one variant is assigned to several genes?
"""
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Convert the gene annotation file to a variant-gene index, based on a specific consequence.")
parser.add_argument("--anno", "-a", help="file with gene annotation (as in BRaVa/SAIGE)", required=True, type=str)
parser.add_argument("--consq", "-c", nargs='+', help="class of consequence to process (list if many)", required=True, type=str)
parser.add_argument("--out", "-o", help="filename for the output", required=True, type=str)
args = parser.parse_args()

print("Will process the set of", args.consq)
gene_data_all = {} # this dict will contain one DF per gene, for all genes in the analysis
r=0
with open( args.anno, 'r') as fin:
    for line in fin:
        tokens = line.split()
        if r%2==0 :
            gene_set_vars = tokens[2:]
        else:
            gene_set_anns = tokens[2:]
            gene_data_all[ tokens[0] ] = pd.DataFrame( {"SNPID":gene_set_vars, "Class": gene_set_anns } ) #.set_index('SNPID')

        r += 1

print( f"Variant annotation loaded for {len(gene_data_all)} genes." )

snps_processed = 0
gens_processed = 0
weird = 0
with open(args.out, 'w') as fout:
    for gene in gene_data_all:
        tmp_df = gene_data_all[ gene ].query( 'Class in @args.consq' )
        if len(tmp_df)>0:
            # print one line per variant
            for x in tmp_df.iterrows():
                fout.write( f'{x[1].SNPID} {gene}\n')
            snps_processed += len(tmp_df)
            gens_processed += 1
        else:
            weird += 1
            # print( f"No {args.consq} snps found in {gene}.")

if weird>0: print(f"Genes without any variants in {args.consq}:{weird}.")

print( f"All done, {snps_processed} variants were processed across {gens_processed} genes." )