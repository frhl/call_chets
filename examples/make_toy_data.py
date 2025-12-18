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
            
            var_id = f"{chrom}:{pos}:{ref}:{alt}"
            f.write(f"{chrom}\t{pos}\t{var_id}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t" + "\t".join(gt_strs) + "\n")

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
        g1 = [1 if random.random() < 0.3 else 0, 1 if random.random() < 0.3 else 0]
        genos1.append(tuple(g1))
        
        # Var 2
        g2 = [1 if random.random() < 0.1 else 0, 1 if random.random() < 0.1 else 0]
        genos2.append(tuple(g2))
        
        # Total alt count (Gene Burden)
        total_alt = sum(g1) + sum(g2)
        sample_genotypes.append(total_alt)

    variants.append(("1", 100, "A", "T", genos1))
    variants.append(("1", 200, "G", "C", genos2))
    
    # Write Phased VCF
    write_vcf("examples/input/phased.vcf.gz", samples, variants, phased=True)
    
    # Write Unphased VCF
    write_vcf("examples/input/unphased.vcf.gz", samples, variants, phased=False)
    
    # Write Gene Map
    with open("examples/input/gene_map.txt", "w") as f:
        f.write("variant\tgene\n")
        for var, gene in gene_map:
            f.write(f"{var}\t{gene}\n")
            
    # Simulate Phenotypes
    # 1. Y_rec: Carrier Recessive Effect (Risk when >= 2 alleles) -> Expect Additive + Dominance signal
    # 2. Y_add: Pure Additive Effect -> Expect Additive signal, Null Dominance
    # 3. Y_null: No Effect -> Expect Null Additive, Null Dominance
    
    with open("examples/input/phenotypes.txt", "w") as f:
        f.write("IID\tY_rec\tY_add\tY_null\n")
        for i, s in enumerate(samples):
            # 1. Recessive
            is_recessive_state = 1 if sample_genotypes[i] >= 2 else 0
            y_rec = 2.0 * is_recessive_state + random.gauss(0, 1.0)
            
            # 2. Additive (linear effect of allele count)
            # genotype is 0, 1, 2, 3, 4 (since 2 variants)
            # Let's say effect size is 0.5 per allele
            y_add = 0.5 * sample_genotypes[i] + random.gauss(0, 1.0)
            
            # 3. Null
            y_null = random.gauss(0, 1.0)
            
            f.write(f"{s}\t{y_rec:.4f}\t{y_add:.4f}\t{y_null:.4f}\n")

if __name__ == "__main__":
    simulate_data()
