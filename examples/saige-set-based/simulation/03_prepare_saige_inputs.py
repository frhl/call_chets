#!/usr/bin/env python3
"""
Prepare input files for SAIGE from simulated data.

This script:
1. Converts gene boundaries to SAIGE group file format
2. Prepares phenotype file in SAIGE format
3. Optionally generates simple covariates
"""

import argparse
import pandas as pd
import numpy as np
from pathlib import Path


def create_saige_group_file(gene_boundaries_file, variant_annotations_file, output_file):
    """
    Create SAIGE group file from gene boundaries and variant annotations.

    SAIGE group file format:
    variant_id   gene_id   annotation
    var_0        gene_0    pLoF
    var_1        gene_0    synonymous
    """
    print("Creating SAIGE group file...")

    # Read variant annotations
    annotations = pd.read_csv(variant_annotations_file, sep='\t')

    # SAIGE group file: variant_id, gene_id, annotation
    group_df = annotations[['variant_id', 'gene_id', 'annotation']].copy()

    # Save
    group_df.to_csv(output_file, sep='\t', index=False, header=False)
    print(f"  Saved SAIGE group file to: {output_file}")
    print(f"  Total variants: {len(group_df)}")
    print(f"  Total genes: {group_df['gene_id'].nunique()}")

    return group_df


def create_saige_phenotype_file(phenotype_file, output_file, add_covariates=True,
                                 pcs_file=None, seed=42):
    """
    Create SAIGE-formatted phenotype file.

    SAIGE phenotype file format (tab-separated):
    IID   phenotype   covariate1   covariate2   ...

    For quantitative traits, phenotype should be pre-residualized or
    covariates should be included.
    """
    print("\nCreating SAIGE phenotype file...")

    # Read phenotypes
    pheno_df = pd.read_csv(phenotype_file, sep='\t')

    # Rename sample_id to IID (SAIGE expects IID column)
    pheno_df = pheno_df.rename(columns={'sample_id': 'IID'})

    if add_covariates:
        print("  Adding simple covariates...")
        np.random.seed(seed)
        n_samples = len(pheno_df)

        # Add simple covariates
        # Age: random between 40-70
        pheno_df['age'] = np.random.uniform(40, 70, n_samples)
        pheno_df['age2'] = pheno_df['age'] ** 2

        # Sex: random binary (0=female, 1=male)
        pheno_df['sex'] = np.random.binomial(1, 0.5, n_samples)

        # Age-sex interaction
        pheno_df['age_sex'] = pheno_df['age'] * pheno_df['sex']
        pheno_df['age2_sex'] = pheno_df['age2'] * pheno_df['sex']

    # Merge with PCs if provided
    if pcs_file and Path(pcs_file).exists():
        print(f"  Merging with PCs from: {pcs_file}")
        pcs_df = pd.read_csv(pcs_file, sep='\t')

        # Merge on IID
        pheno_df = pheno_df.merge(pcs_df, on='IID', how='left')
        print(f"  Added {len([c for c in pcs_df.columns if c.startswith('PC')])} PCs")

    # Save
    pheno_df.to_csv(output_file, sep='\t', index=False)
    print(f"  Saved SAIGE phenotype file to: {output_file}")
    print(f"  Samples: {len(pheno_df)}")
    print(f"  Columns: {', '.join(pheno_df.columns)}")

    return pheno_df


def beta_weight(maf, a=1, b=25):
    """
    Calculate Beta weights based on MAF.
    Weight = dbeta(MAF, a, b) using numpy (no scipy dependency).
    Beta PDF: f(x; a, b) = x^(a-1) * (1-x)^(b-1) / B(a, b)
    For a=1, b=25: f(x) = (1-x)^24 / B(1, 25) = 25 * (1-x)^24
    """
    import math
    maf = np.array(maf)
    # Beta function B(a, b) = Gamma(a) * Gamma(b) / Gamma(a+b)
    # For a=1: B(1, b) = 1/b
    beta_func = math.gamma(a) * math.gamma(b) / math.gamma(a + b)
    # Avoid division by zero for MAF=0 or MAF=1
    weights = np.where(
        (maf > 0) & (maf < 1),
        (maf ** (a - 1)) * ((1 - maf) ** (b - 1)) / beta_func,
        0.0
    )
    return weights


def create_saige_geneset_file(gene_boundaries_file, variant_annotations_file,
                               output_file, annotation_filter='pLoF', vcf_file=None):
    """
    Create SAIGE geneset file - one file per gene with variant IDs, annotations, and weights.

    Format (3 lines per gene):
    gene_id var var_1 var_2 ...
    gene_id anno anno_1 anno_2 ...
    gene_id weight weight_1 weight_2 ...
    """
    print(f"\nCreating SAIGE geneset file (filter: {annotation_filter})...")

    # Read annotations
    annotations = pd.read_csv(variant_annotations_file, sep='\t')

    # If VCF file is provided, map variant IDs to chr:pos:ref:alt format
    if vcf_file and Path(vcf_file).exists():
        print(f"  Reading VCF to create chr:pos:ref:alt IDs: {vcf_file}")
        import gzip

        # Read VCF and create mapping from old ID to chr:pos:ref:alt
        var_id_map = {}
        open_func = gzip.open if vcf_file.endswith('.gz') else open
        with open_func(vcf_file, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                chrom, pos, var_id, ref, alt = fields[0:5]
                # Create chr:pos:ref:alt format
                new_id = f"{chrom}:{pos}:{ref}:{alt}"
                var_id_map[var_id] = new_id

        # Replace variant_id in annotations
        annotations['variant_id'] = annotations['variant_id'].map(var_id_map)
        print(f"  Mapped {len(var_id_map)} variant IDs to chr:pos:ref:alt format")

    # Calculate weights from existing MAF
    if 'maf' in annotations.columns:
        print("  Using existing MAF column for weights")
        annotations['maf'] = annotations['maf'].fillna(0)
        annotations['weight'] = beta_weight(annotations['maf'])
    else:
        print("  WARNING: 'maf' column not found. Weights will be set to 1.0")
        annotations['weight'] = 1.0

    # Filter by annotation if specified
    if annotation_filter:
        annotations = annotations[annotations['annotation'] == annotation_filter].copy()
        print(f"  Filtered to {len(annotations)} {annotation_filter} variants")

    # Group by gene
    # We need to aggregate variant_id, annotation, and weight
    grouped = annotations.groupby('gene_id').agg({
        'variant_id': list,
        'annotation': list,
        'weight': list
    }).reset_index()

    print(f"  Writing weighted group file to: {output_file}")
    with open(output_file, 'w') as f:
        for _, row in grouped.iterrows():
            gene_id = row['gene_id']
            vars_list = row['variant_id']
            annos_list = row['annotation']
            weights_list = row['weight']
            
            # Format: 
            # gene var var1 var2 ...
            # gene anno ann1 ann2 ...
            # gene weight w1 w2 ...
            
            # Line 1: Variants
            f.write(f"{gene_id} var " + " ".join(vars_list) + "\n")
            
            # Line 2: Annotations
            f.write(f"{gene_id} anno " + " ".join(annos_list) + "\n")
            
            # Line 3: Weights (formatted to reasonable precision)
            weights_str = " ".join([f"{w:.6f}" for w in weights_list])
            f.write(f"{gene_id} weight " + weights_str + "\n")

    print(f"  Saved SAIGE geneset file to: {output_file}")
    print(f"  Total genes: {len(grouped)}")

    return grouped


def main():
    parser = argparse.ArgumentParser(
        description='Prepare SAIGE input files from simulated data'
    )
    parser.add_argument('--gene_boundaries', type=str,
                       default='input/gene_boundaries.tsv',
                       help='Gene boundaries file')
    parser.add_argument('--variant_annotations', type=str,
                       default='input/variant_annotations.tsv',
                       help='Variant annotations file')
    parser.add_argument('--phenotype_file', type=str,
                       default='input/simulated.phenos.tsv',
                       help='Phenotype file')
    parser.add_argument('--output_dir', type=str, default='input',
                       help='Output directory for SAIGE files')
    parser.add_argument('--annotation_filter', type=str, default='pLoF',
                       help='Annotation type to include in geneset (pLoF, synonymous, or all)')
    parser.add_argument('--add_covariates', action='store_true', default=False,
                       help='Add simple covariates to phenotype file')
    parser.add_argument('--pcs_file', type=str, default='input/pca/pcs.txt',
                       help='PCs file to merge with phenotypes (optional)')
    parser.add_argument('--vcf_file', type=str, default=None,
                       help='VCF file to extract chr:pos:ref:alt variant IDs (optional)')
    parser.add_argument('--seed', type=int, default=42,
                       help='Random seed for covariate generation')

    args = parser.parse_args()

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("Preparing SAIGE Input Files")
    print("=" * 60)

    # 1. Create SAIGE group file (variant to gene mapping)
    # Note: This file seems redundant if we are using the set-based test file 
    # but SAIGE step 2 might not need this one specifically if we pass the groupFile. 
    # However, sometimes it's used for other checks. We'll leave it as is for now.
    group_file = output_dir / 'variant_gene_mapping.txt'
    create_saige_group_file(
        args.variant_annotations,
        args.variant_annotations,
        group_file
    )

    # 2. Create SAIGE phenotype file
    pheno_file = output_dir / 'simulated.phenos.with_covariates.tsv'
    create_saige_phenotype_file(
        args.phenotype_file,
        pheno_file,
        add_covariates=args.add_covariates,
        pcs_file=args.pcs_file,
        seed=args.seed
    )

    # 3. Create SAIGE geneset file
    if args.annotation_filter == 'all':
        annotation_filter = None
    else:
        annotation_filter = args.annotation_filter

    geneset_file = output_dir / f'genesets_{args.annotation_filter}.txt'

    create_saige_geneset_file(
        args.gene_boundaries,
        args.variant_annotations,
        geneset_file,
        annotation_filter=annotation_filter,
        vcf_file=args.vcf_file
    )

    print("\n" + "=" * 60)
    print("SAIGE input files prepared successfully!")
    print("=" * 60)
    print(f"\nOutput files in {output_dir}:")
    print(f"  - variant_gene_mapping.txt (variant-gene mapping)")
    print(f"  - simulated.phenos.with_covariates.tsv (phenotypes and covariates)")
    print(f"  - genesets_{args.annotation_filter}.txt (gene variant sets with weights)")


if __name__ == '__main__':
    main()
