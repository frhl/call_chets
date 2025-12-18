import gzip
import random
import sys
import os

def write_vcf(filename, samples, variants, phased=True):
    with gzip.open(filename, 'wt') as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("##contig=<ID=1,length=1000>\n")
        f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
        
        for chrom, pos, ref, alt, genos in variants:
            # genos is a list of (a, b) tuples
            gt_strs = []
            for a, b in genos:
                sep = "|" if phased else "/"
                gt_strs.append(f"{a}{sep}{b}")
            
            f.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t" + "\t".join(gt_strs) + "\n")

def simulate_data():
    random.seed(42)
    n_samples = 500
    samples = [f"sample_{i}" for i in range(n_samples)]
    
    # Gene A: 2 variants. 
    variants = []
    # Variant 1: 1:100:A:T (Common)
    genos1 = []
    # Variant 2: 1:200:G:C (Rare)
    genos2 = []
    
    # Store allele counts for phenotype sim
    sample_genotypes = [] 
    
    gene_map = [
        ("1:100:A:T", "GeneA"),
        ("1:200:G:C", "GeneA")
    ]
    
    for i in range(n_samples):
        # Var 1
        g1 = sorted([1 if random.random() < 0.3 else 0, 1 if random.random() < 0.3 else 0])
        genos1.append(tuple(g1))
        
        # Var 2
        g2 = sorted([1 if random.random() < 0.1 else 0, 1 if random.random() < 0.1 else 0])
        genos2.append(tuple(g2))
        
        # Total alt count (Gene Burden)
        total_alt = sum(g1) + sum(g2)
        sample_genotypes.append(total_alt)

    variants.append(("1", 100, "A", "T", genos1))
    variants.append(("1", 200, "G", "C", genos2))
    
    # Write Phased VCF
    write_vcf("examples/phased.vcf.gz", samples, variants, phased=True)
    
    # Write Unphased VCF
    write_vcf("examples/unphased.vcf.gz", samples, variants, phased=False)
    
    # Write Gene Map
    with open("examples/gene_map.txt", "w") as f:
        f.write("variant\tgene\n")
        for var, gene in gene_map:
            f.write(f"{var}\t{gene}\n")
            
    # Simulate Phenotype - Carrier Recessive Effect
    # Y = 1 if user has >= 2 alleles (homozygous or CH), else 0
    # Additive effect 0.
    with open("examples/phenotypes.txt", "w") as f:
        f.write("IID\tY\n")
        for i, s in enumerate(samples):
            # Recessive trait: specific effect only when >= 2 copies
            is_recessive_state = 1 if sample_genotypes[i] >= 2 else 0
            
            # Use a simpler continuous phenotype for regression
            # Mean shift for recessive state
            y = 2.0 * is_recessive_state + random.gauss(0, 1.0)
            f.write(f"{s}\t{y:.4f}\n")

if __name__ == "__main__":
    simulate_data()
