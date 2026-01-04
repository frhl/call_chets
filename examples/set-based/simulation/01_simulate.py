#!/usr/bin/env python3
"""
Simulate genotypes using msprime/tskit for power analysis.

This script simulates:
- 5000 samples from a single population
- ~1000 variants with a large fraction having MAF 1-5%
- Variants organized into 20 "genes" of 50 variants each
- Variant annotations (pLoF vs synonymous)
- Quantitative phenotypes based on causal variants
"""

import argparse
import numpy as np
import pandas as pd
import msprime
import tskit
from pathlib import Path


def simulate_gene(
    gene_idx,
    n_samples=5000,
    gene_length=35_000,
    target_variants=100,
    recombination_rate=1e-8,
    mutation_rate=2e-8,
    population_size=10_000,
    seed=42
):
    """
    Simulate a single gene independently using msprime.

    Parameters:
    -----------
    gene_idx : int
        Gene index for seed offset
    n_samples : int
        Number of samples to simulate
    gene_length : int
        Length of gene region (bp)
    target_variants : int
        Target number of variants for this gene
    recombination_rate : float
        Recombination rate per base pair per generation
    mutation_rate : float
        Mutation rate per base pair per generation
    population_size : int
        Effective population size
    seed : int
        Base random seed

    Returns:
    --------
    variant_stats_df : pd.DataFrame
        Variant statistics for this gene
    genotype_matrix : np.ndarray
        Diploid genotype matrix for this gene
    """
    # Use gene-specific seed to ensure independence
    gene_seed = seed + gene_idx * 1000

    # Simulate ancestry
    ts = msprime.sim_ancestry(
        samples=n_samples,
        population_size=population_size,
        recombination_rate=recombination_rate,
        sequence_length=gene_length,
        random_seed=gene_seed
    )

    # Add mutations
    ts = msprime.sim_mutations(
        ts,
        rate=mutation_rate,
        random_seed=gene_seed
    )

    # Convert to diploid genotypes
    haplotype_matrix = ts.genotype_matrix().T
    n_haplotypes, n_variants = haplotype_matrix.shape
    n_diploid_samples = n_haplotypes // 2

    genotype_matrix = np.zeros((n_diploid_samples, n_variants), dtype=int)
    for i in range(n_diploid_samples):
        hap1 = haplotype_matrix[2*i, :]
        hap2 = haplotype_matrix[2*i+1, :]
        genotype_matrix[i, :] = hap1 + hap2

    # Compute variant stats
    variant_stats = []
    for var_idx, variant in enumerate(ts.variants()):
        genotypes = genotype_matrix[:, var_idx]

        n_0 = np.sum(genotypes == 0)
        n_1 = np.sum(genotypes == 1)
        n_2 = np.sum(genotypes == 2)

        total_alleles = 2 * n_diploid_samples
        alt_allele_count = n_1 + 2 * n_2
        af = alt_allele_count / total_alleles
        maf = min(af, 1 - af)

        variant_stats.append({
            'position': int(variant.site.position),
            'n_hom_ref': n_0,
            'n_het': n_1,
            'n_hom_alt': n_2,
            'af': af,
            'maf': maf
        })

    variant_stats_df = pd.DataFrame(variant_stats)

    # Track original count
    n_variants_original = len(variant_stats_df)

    # Downsample to target number of variants if we have more than target
    if n_variants_original > target_variants:
        # Randomly sample target_variants
        np.random.seed(gene_seed + 1)  # Different seed for sampling
        selected_indices = np.random.choice(n_variants_original, size=target_variants, replace=False)
        selected_indices = np.sort(selected_indices)  # Keep genomic order

        variant_stats_df = variant_stats_df.iloc[selected_indices].reset_index(drop=True)
        genotype_matrix = genotype_matrix[:, selected_indices]

    n_variants_final = len(variant_stats_df)

    return variant_stats_df, genotype_matrix, ts.num_mutations, n_variants_original, n_variants_final


def simulate_independent_genes(
    n_genes=20,
    n_samples=5000,
    variants_per_gene=100,
    gene_length_min=30_000,
    gene_length_max=40_000,
    intergenic_spacing=1_000_000,
    recombination_rate=1e-8,
    mutation_rate=2e-8,
    population_size=10_000,
    seed=42,
    chromosome="1"
):
    """
    Simulate multiple genes independently (zero LD between genes).

    Each gene is simulated separately with msprime, then combined
    with large genomic spacing to ensure independence.

    Parameters:
    -----------
    n_genes : int
        Number of independent genes to simulate
    n_samples : int
        Number of samples
    variants_per_gene : int
        Target variants per gene
    gene_length_min : int
        Minimum gene length (bp)
    gene_length_max : int
        Maximum gene length (bp)
    intergenic_spacing : int
        Space between genes on chromosome (bp)
    recombination_rate : float
        Recombination rate per bp
    mutation_rate : float
        Mutation rate per bp
    population_size : int
        Effective population size
    seed : int
        Random seed

    Returns:
    --------
    variant_stats_df : pd.DataFrame
        Combined variant statistics
    genotype_matrix : np.ndarray
        Combined genotype matrix
    gene_boundaries : pd.DataFrame
        Gene boundary information
    """
    print(f"Simulating {n_genes} independent genes...", flush=True)
    print(f"  Samples: {n_samples}", flush=True)
    print(f"  Variants per gene: ~{variants_per_gene}", flush=True)
    print(f"  Gene length: {gene_length_min:,} - {gene_length_max:,} bp", flush=True)
    print(f"  Intergenic spacing: {intergenic_spacing:,} bp", flush=True)
    print(f"  Population size: {population_size:,}", flush=True)
    print(f"  Mutation rate: {mutation_rate}", flush=True)
    print(f"  Recombination rate: {recombination_rate}", flush=True)
    print(f"  Random seed: {seed}", flush=True)

    # Set seed for reproducible gene lengths
    np.random.seed(seed)

    all_variants = []
    all_genotypes = []
    gene_boundaries = []
    current_position = 0

    for gene_idx in range(n_genes):
        # Random gene length for this gene
        gene_length = np.random.randint(gene_length_min, gene_length_max + 1)

        print(f"\n  Simulating gene_{gene_idx} ({gene_length:,} bp)...", flush=True)

        # Simulate this gene
        variant_stats, genotypes, n_mutations, n_original, n_final = simulate_gene(
            gene_idx=gene_idx,
            n_samples=n_samples,
            gene_length=gene_length,
            target_variants=variants_per_gene,
            recombination_rate=recombination_rate,
            mutation_rate=mutation_rate,
            population_size=population_size,
            seed=seed
        )

        # Adjust positions to place gene on chromosome
        gene_start = current_position
        gene_end = current_position + gene_length
        variant_stats['position'] = variant_stats['position'] + gene_start
        variant_stats['gene_id'] = f'gene_{gene_idx}'

        # Print variant count information
        if n_original > n_final:
            print(f"    Generated {n_original} variants → downsampled to {n_final} (MAF range: {variant_stats['maf'].min():.4f} - {variant_stats['maf'].max():.4f})", flush=True)
        else:
            print(f"    Generated {n_final} variants (MAF range: {variant_stats['maf'].min():.4f} - {variant_stats['maf'].max():.4f})", flush=True)

        # Store gene boundaries
        gene_boundaries.append({
            'gene_id': f'gene_{gene_idx}',
            'gene_index': gene_idx,
            'start_position': gene_start,
            'end_position': gene_end,
            'length': gene_length,
            'n_variants': len(variant_stats)
        })

        # Accumulate
        all_variants.append(variant_stats)
        all_genotypes.append(genotypes)

        # Move to next gene position
        current_position = gene_end + intergenic_spacing

    # Combine all genes
    print(f"\n  Combining {n_genes} genes...", flush=True)
    variant_stats_df = pd.concat(all_variants, ignore_index=True)
    genotype_matrix = np.hstack(all_genotypes)

    # Assign variant IDs in chr:pos:ref:alt format
    variant_stats_df['variant_id'] = [
        f"{chromosome}:{pos}:A:T"
        for pos in variant_stats_df['position']
    ]

    gene_boundaries_df = pd.DataFrame(gene_boundaries)

    print(f"\n  Total variants: {len(variant_stats_df)}", flush=True)
    print(f"  Total genes: {len(gene_boundaries_df)}", flush=True)
    print(f"  Genotype matrix: {genotype_matrix.shape[0]} samples × {genotype_matrix.shape[1]} variants", flush=True)

    return variant_stats_df, genotype_matrix, gene_boundaries_df


def compute_variant_stats(ts):
    """
    Compute allele frequencies and other statistics for variants.
    Converts haploid samples to diploid genotypes.

    Returns:
    --------
    variant_stats : pd.DataFrame
        DataFrame with variant statistics including position, MAF, etc.
    genotype_matrix : np.ndarray
        Diploid genotype matrix (individuals x variants)
    """
    print("\nComputing variant statistics...")

    # Get haplotype matrix
    haplotype_matrix = ts.genotype_matrix().T  # Transpose to haplotypes x variants

    n_haplotypes, n_variants = haplotype_matrix.shape
    n_diploid_samples = n_haplotypes // 2

    print(f"  Haplotypes: {n_haplotypes}")
    print(f"  Diploid individuals: {n_diploid_samples}")
    print(f"  Variants: {n_variants}")

    # Convert haplotypes to diploid genotypes
    # Combine pairs of haplotypes (0,1), (2,3), (4,5), etc.
    genotype_matrix = np.zeros((n_diploid_samples, n_variants), dtype=int)

    for i in range(n_diploid_samples):
        hap1 = haplotype_matrix[2*i, :]     # First haplotype
        hap2 = haplotype_matrix[2*i+1, :]   # Second haplotype
        genotype_matrix[i, :] = hap1 + hap2  # Genotype is sum of haplotypes

    print(f"  Genotype matrix shape: {n_diploid_samples} samples x {n_variants} variants")

    # Compute allele frequencies
    variant_stats = []

    for var_idx, variant in enumerate(ts.variants()):
        genotypes = genotype_matrix[:, var_idx]

        # Count genotypes
        n_0 = np.sum(genotypes == 0)  # Homozygous ref
        n_1 = np.sum(genotypes == 1)  # Heterozygous
        n_2 = np.sum(genotypes == 2)  # Homozygous alt

        # Compute allele frequency
        total_alleles = 2 * n_diploid_samples
        alt_allele_count = n_1 + 2 * n_2
        af = alt_allele_count / total_alleles
        maf = min(af, 1 - af)

        variant_stats.append({
            'variant_id': f'var_{var_idx}',
            'position': int(variant.site.position),
            'n_hom_ref': n_0,
            'n_het': n_1,
            'n_hom_alt': n_2,
            'af': af,
            'maf': maf
        })

    variant_stats_df = pd.DataFrame(variant_stats)

    print(f"  Total variants: {len(variant_stats_df)}")
    print(f"  MAF range: [{variant_stats_df['maf'].min():.4f}, {variant_stats_df['maf'].max():.4f}]")
    print(f"  Mean MAF: {variant_stats_df['maf'].mean():.4f}")

    return variant_stats_df, genotype_matrix


def filter_variants_by_maf(variant_stats_df, genotype_matrix,
                           target_n=1000, maf_range=(0.01, 0.05),
                           target_fraction=0.7, seed=42):
    """
    Filter variants to get target number with specified MAF distribution.

    Parameters:
    -----------
    variant_stats_df : pd.DataFrame
        Variant statistics
    genotype_matrix : np.ndarray
        Full genotype matrix
    target_n : int
        Target number of variants
    maf_range : tuple
        MAF range for enrichment (min, max)
    target_fraction : float
        Target fraction of variants in MAF range
    seed : int
        Random seed

    Returns:
    --------
    filtered_stats : pd.DataFrame
        Filtered variant statistics
    filtered_genotypes : np.ndarray
        Filtered genotype matrix
    """
    print(f"\nFiltering variants to target {target_n} variants...")
    print(f"  Target fraction in MAF range [{maf_range[0]}, {maf_range[1]}]: {target_fraction}")

    np.random.seed(seed)

    # Separate variants by MAF
    in_range = (variant_stats_df['maf'] >= maf_range[0]) & (variant_stats_df['maf'] < maf_range[1])
    out_range = ~in_range

    variants_in_range = variant_stats_df[in_range]
    variants_out_range = variant_stats_df[out_range]

    print(f"  Variants in MAF range: {len(variants_in_range)}")
    print(f"  Variants out of range: {len(variants_out_range)}")

    # Determine how many variants to select from each group
    n_from_range = int(target_n * target_fraction)
    n_from_out = target_n - n_from_range

    # Sample from each group
    if len(variants_in_range) < n_from_range:
        print(f"  WARNING: Only {len(variants_in_range)} variants in range, adjusting...")
        n_from_range = len(variants_in_range)
        n_from_out = target_n - n_from_range

    if len(variants_out_range) < n_from_out:
        print(f"  WARNING: Only {len(variants_out_range)} variants out of range, adjusting...")
        n_from_out = len(variants_out_range)

    # Random sample from each group
    selected_in_range = variants_in_range.sample(n=n_from_range, random_state=seed)
    selected_out_range = variants_out_range.sample(n=n_from_out, random_state=seed)

    # Combine and sort by position
    filtered_stats = pd.concat([selected_in_range, selected_out_range])
    filtered_stats = filtered_stats.sort_values('position').reset_index(drop=True)

    # Update variant IDs to be sequential
    filtered_stats['variant_id'] = [f'var_{i}' for i in range(len(filtered_stats))]

    # Extract corresponding genotypes
    variant_indices = filtered_stats.index.values
    filtered_genotypes = genotype_matrix[:, variant_indices]

    # Final statistics
    n_in_range_final = np.sum((filtered_stats['maf'] >= maf_range[0]) &
                               (filtered_stats['maf'] < maf_range[1]))

    print(f"\nFiltered variant statistics:")
    print(f"  Total variants: {len(filtered_stats)}")
    print(f"  Variants in MAF range [{maf_range[0]}, {maf_range[1]}]: {n_in_range_final} ({n_in_range_final/len(filtered_stats)*100:.1f}%)")
    print(f"  MAF range: [{filtered_stats['maf'].min():.4f}, {filtered_stats['maf'].max():.4f}]")
    print(f"  Mean MAF: {filtered_stats['maf'].mean():.4f}")

    return filtered_stats, filtered_genotypes


def create_gene_boundaries(variant_stats_df, variants_per_gene=50,
                          avg_intergenic_distance=250_000,
                          gene_length=50_000,
                          output_file=None, seed=42):
    """
    Assign variants to genes with random spacing between genes.

    Parameters:
    -----------
    variant_stats_df : pd.DataFrame
        Variant statistics
    variants_per_gene : int
        Target number of variants per gene
    avg_intergenic_distance : int
        Average distance between genes in bp (default: 250kb)
    gene_length : int
        Approximate length of each gene in bp (default: 50kb)
    output_file : str or None
        Output file path
    seed : int
        Random seed for spacing

    Returns:
    --------
    gene_df : pd.DataFrame
        DataFrame with gene boundaries
    variant_stats_df : pd.DataFrame
        Updated variant stats with gene assignments
    """
    print(f"\nCreating gene boundaries with spacing...")
    print(f"  Target variants per gene: {variants_per_gene}")
    print(f"  Average intergenic distance: {avg_intergenic_distance:,} bp")
    print(f"  Gene length: {gene_length:,} bp")

    np.random.seed(seed)

    n_variants = len(variant_stats_df)
    n_genes = n_variants // variants_per_gene

    print(f"  Total variants: {n_variants}")
    print(f"  Number of genes: {n_genes}")

    # Define gene regions with random spacing
    genes = []
    current_pos = 0

    for gene_idx in range(n_genes):
        # Add random intergenic distance (exponential distribution centered on avg)
        if gene_idx > 0:
            intergenic_gap = int(np.random.exponential(avg_intergenic_distance))
            current_pos += intergenic_gap

        gene_start = current_pos
        gene_end = current_pos + gene_length

        # Find variants within this gene region
        variants_in_gene = variant_stats_df[
            (variant_stats_df['position'] >= gene_start) &
            (variant_stats_df['position'] < gene_end)
        ]

        if len(variants_in_gene) > 0:
            genes.append({
                'gene_id': f'gene_{gene_idx}',
                'gene_index': gene_idx,
                'start_position': gene_start,
                'end_position': gene_end,
                'n_variants': len(variants_in_gene),
                'variants': ','.join(variants_in_gene['variant_id'].values)
            })

        # Move to end of current gene
        current_pos = gene_end

    gene_df = pd.DataFrame(genes)

    # Assign gene IDs to variants
    variant_stats_df = variant_stats_df.copy()
    variant_stats_df['gene_id'] = 'intergenic'

    for _, gene in gene_df.iterrows():
        mask = (
            (variant_stats_df['position'] >= gene['start_position']) &
            (variant_stats_df['position'] < gene['end_position'])
        )
        variant_stats_df.loc[mask, 'gene_id'] = gene['gene_id']

    # Summary statistics
    variants_per_gene_actual = gene_df['n_variants'].values
    intergenic_variants = (variant_stats_df['gene_id'] == 'intergenic').sum()

    print(f"\n  Gene statistics:")
    print(f"    Genes created: {len(gene_df)}")
    print(f"    Variants per gene: {variants_per_gene_actual.mean():.1f} ± {variants_per_gene_actual.std():.1f}")
    print(f"    Range: [{variants_per_gene_actual.min()}, {variants_per_gene_actual.max()}]")
    print(f"    Intergenic variants: {intergenic_variants} ({intergenic_variants/len(variant_stats_df)*100:.1f}%)")

    # Calculate actual spacing
    if len(gene_df) > 1:
        actual_spacing = np.diff([gene['start_position'] for _, gene in gene_df.iterrows()])
        print(f"    Intergenic spacing: {actual_spacing.mean():,.0f} ± {actual_spacing.std():,.0f} bp")

    if output_file:
        gene_df.to_csv(output_file, sep='\t', index=False)
        print(f"  Saved to: {output_file}")

    return gene_df, variant_stats_df


def create_variant_annotations(variant_stats_df, gene_df,
                                causal_gene_fraction=0.1,
                                causal_variants_per_gene=5,
                                pLoF_maf_min=0.0,
                                pLoF_maf_max=1.0,
                                seed=42,
                                output_file=None):
    """
    Annotate variants as pLoF (causal) or synonymous (non-causal).

    Uses a gene-based approach:
    1. Select a fraction of genes to be causal
    2. Within each causal gene, select N variants as causal (pLoF) from MAF window

    Parameters:
    -----------
    variant_stats_df : pd.DataFrame
        Variant statistics
    gene_df : pd.DataFrame
        Gene boundaries
    causal_gene_fraction : float
        Fraction of genes that contain causal variants (default: 0.1 = 10%)
    causal_variants_per_gene : int
        Average number of causal variants per causal gene (default: 5)
    pLoF_maf_min : float
        Minimum MAF for selecting pLoF variants (default: 0.0)
    pLoF_maf_max : float
        Maximum MAF for selecting pLoF variants (default: 1.0)
    seed : int
        Random seed
    output_file : str or None
        Output file path

    Returns:
    --------
    annotation_df : pd.DataFrame
        Variant annotations
    """
    print(f"\nCreating variant annotations (gene-based causal model)...")
    print(f"  Causal gene fraction: {causal_gene_fraction} ({causal_gene_fraction*100:.1f}%)")
    print(f"  Causal variants per gene: ~{causal_variants_per_gene}")
    print(f"  pLoF MAF window: [{pLoF_maf_min}, {pLoF_maf_max}]")

    np.random.seed(seed)

    n_genes = len(gene_df)
    n_causal_genes = max(1, int(n_genes * causal_gene_fraction))

    # Step 1: Pre-filter genes to only those with variants in MAF window
    eligible_gene_indices = []
    for gene_idx in range(n_genes):
        gene_id = f'gene_{gene_idx}'
        gene_variants = variant_stats_df[variant_stats_df['gene_id'] == gene_id]

        # Check if this gene has at least one variant in MAF window
        maf_filtered = gene_variants[
            (gene_variants['maf'] >= pLoF_maf_min) &
            (gene_variants['maf'] <= pLoF_maf_max)
        ]

        if len(maf_filtered) > 0:
            eligible_gene_indices.append(gene_idx)

    print(f"  Total genes: {n_genes}")
    print(f"  Genes with variants in MAF window [{pLoF_maf_min}, {pLoF_maf_max}]: {len(eligible_gene_indices)}")

    # Check if we have enough eligible genes
    if len(eligible_gene_indices) < n_causal_genes:
        print(f"  WARNING: Only {len(eligible_gene_indices)} genes have variants in MAF window, but {n_causal_genes} causal genes requested.")
        print(f"  Using all {len(eligible_gene_indices)} eligible genes as causal.")
        n_causal_genes = len(eligible_gene_indices)

    # Step 2: Select causal genes from eligible genes only
    causal_gene_indices = np.random.choice(eligible_gene_indices, size=n_causal_genes, replace=False)
    causal_genes = set([f'gene_{idx}' for idx in causal_gene_indices])

    print(f"  Causal genes selected: {n_causal_genes}")
    print(f"  Causal gene IDs: {sorted(causal_genes)}")

    # Step 3: For each causal gene, select causal variants from MAF window
    causal_variant_ids = set()

    for gene_idx in causal_gene_indices:
        gene_id = f'gene_{gene_idx}'

        # Get variants in this gene (using gene_id column)
        gene_variants = variant_stats_df[variant_stats_df['gene_id'] == gene_id]

        # Filter variants by MAF window for pLoF selection
        maf_filtered_variants = gene_variants[
            (gene_variants['maf'] >= pLoF_maf_min) &
            (gene_variants['maf'] <= pLoF_maf_max)
        ]

        # Select causal variants within this gene from MAF window
        # (guaranteed to have at least 1 variant due to pre-filtering)
        n_to_select = min(causal_variants_per_gene, len(maf_filtered_variants))
        selected_variant_ids = np.random.choice(
            maf_filtered_variants['variant_id'].values,
            size=n_to_select,
            replace=False
        )
        causal_variant_ids.update(selected_variant_ids)

    # Step 3: Annotate all variants
    annotations = []
    for idx, row in variant_stats_df.iterrows():
        # Use gene_id from variant_stats_df (assigned by create_gene_boundaries)
        gene_id = row.get('gene_id', 'intergenic')

        # Determine annotation
        if gene_id == 'intergenic':
            annotation_type = 'intergenic'
            is_causal = False
        else:
            is_causal = row['variant_id'] in causal_variant_ids
            annotation_type = 'pLoF' if is_causal else 'synonymous'

        annotation = {
            'variant_id': row['variant_id'],
            'position': row['position'],
            'maf': row['maf'],
            'n_hom_ref': row['n_hom_ref'],
            'n_het': row['n_het'],
            'n_hom_alt': row['n_hom_alt'],
            'gene_id': gene_id,
            'annotation': annotation_type,
            'is_causal': is_causal
        }
        annotations.append(annotation)

    annotation_df = pd.DataFrame(annotations)

    # Summary statistics
    n_plof = (annotation_df['annotation'] == 'pLoF').sum()
    n_synonymous = (annotation_df['annotation'] == 'synonymous').sum()
    n_intergenic = (annotation_df['annotation'] == 'intergenic').sum()
    causal_genes_final = annotation_df[annotation_df['is_causal']]['gene_id'].nunique()

    print(f"\n  Summary:")
    print(f"    Total variants: {len(annotation_df)}")
    print(f"    pLoF (causal): {n_plof} ({n_plof/len(annotation_df)*100:.1f}%)")
    print(f"    Synonymous (non-causal): {n_synonymous} ({n_synonymous/len(annotation_df)*100:.1f}%)")
    print(f"    Intergenic: {n_intergenic} ({n_intergenic/len(annotation_df)*100:.1f}%)")
    print(f"    Genes with causal variants: {causal_genes_final}")

    # Show causal variants per gene
    causal_per_gene = annotation_df[annotation_df['is_causal']].groupby('gene_id').size()
    print(f"    Causal variants per causal gene:")
    for gene_id in sorted(causal_per_gene.index):
        print(f"      {gene_id}: {causal_per_gene[gene_id]} variants")

    if output_file:
        annotation_df.to_csv(output_file, sep='\t', index=False)
        print(f"  Saved to: {output_file}")

    return annotation_df


def simulate_phenotypes(genotype_matrix, annotation_df,
                        h2_total=0.1, architecture='additive',
                        seed=42):
    """
    Simulate quantitative phenotypes based on causal variants.

    Parameters:
    -----------
    genotype_matrix : np.ndarray
        Genotype matrix (samples x variants)
    annotation_df : pd.DataFrame
        Variant annotations with is_causal column
    h2_total : float
        Total heritability
    architecture : str
        Genetic architecture ('additive', 'recessive', 'dominant')
    seed : int
        Random seed

    Returns:
    --------
    phenotypes : pd.DataFrame
        Simulated phenotypes
    """
    print(f"\nSimulating phenotypes...")
    print(f"  Heritability (h²): {h2_total}")
    print(f"  Architecture: {architecture}")

    np.random.seed(seed)

    n_samples = genotype_matrix.shape[0]

    # Get causal variants
    causal_indices = annotation_df[annotation_df['is_causal']].index.values
    n_causal = len(causal_indices)

    print(f"  Causal variants: {n_causal}")
    print(f"  Causal variant indices: {causal_indices[:10]}..." if len(causal_indices) > 10 else f"  Causal variant indices: {causal_indices}")

    if n_causal == 0:
        print("  WARNING: No causal variants, returning random phenotypes")
        phenotypes = np.random.randn(n_samples)
    else:
        # Extract causal genotypes
        causal_genotypes = genotype_matrix[:, causal_indices]

        # Apply genetic architecture
        if architecture == 'recessive':
            # Only homozygous alt (2) contributes
            encoded_genotypes = (causal_genotypes == 2).astype(float)
        elif architecture == 'dominant':
            # Het (1) or hom alt (2) contributes equally
            encoded_genotypes = (causal_genotypes >= 1).astype(float)
        else:  # additive
            # Use raw genotypes [0, 1, 2]
            encoded_genotypes = causal_genotypes.astype(float)

        # Standardize genotypes (mean 0, variance 1)
        encoded_genotypes = (encoded_genotypes - encoded_genotypes.mean(axis=0)) / (encoded_genotypes.std(axis=0) + 1e-8)

        # Effect sizes (equal effects per variant)
        # Variance of genetic component should be h2_total
        effect_sizes = np.random.randn(n_causal)

        # Compute genetic component
        genetic_component = encoded_genotypes @ effect_sizes

        # Scale genetic component to have variance = h2_total
        genetic_var = np.var(genetic_component)
        print(f"  Genetic variance before scaling: {genetic_var:.6f}")

        if genetic_var > 0:
            genetic_component = genetic_component * np.sqrt(h2_total / genetic_var)
        else:
            print("  WARNING: Genetic variance is zero! Causal variants have no variation.")

        # Environmental component
        env_var = 1 - h2_total
        environmental_component = np.random.randn(n_samples) * np.sqrt(env_var)

        # Total phenotype
        phenotypes = genetic_component + environmental_component

        # Calculate realized heritability
        realized_h2 = np.var(genetic_component) / np.var(phenotypes)
        print(f"  Realized heritability: {realized_h2:.4f} (target: {h2_total})")

        # Standardize final phenotype
        phenotypes = (phenotypes - phenotypes.mean()) / phenotypes.std()

        print(f"  Phenotype mean: {phenotypes.mean():.4f}")
        print(f"  Phenotype std: {phenotypes.std():.4f}")

    # Create output dataframe
    phenotype_df = pd.DataFrame({
        'sample_id': [f'sample_{i}' for i in range(n_samples)],
        'phenotype': phenotypes
    })

    return phenotype_df


def save_genotypes(genotype_matrix, variant_stats_df, sample_ids, output_prefix):
    """
    Save genotypes in simple format (TSV).

    Parameters:
    -----------
    genotype_matrix : np.ndarray
        Genotype matrix (samples x variants)
    variant_stats_df : pd.DataFrame
        Variant statistics
    sample_ids : list
        Sample IDs
    output_prefix : str
        Output file prefix
    """
    print(f"\nSaving genotypes...")

    # Create genotype dataframe
    genotype_df = pd.DataFrame(
        genotype_matrix,
        index=sample_ids,
        columns=variant_stats_df['variant_id'].values
    )

    output_file = f"{output_prefix}.genotypes.tsv.gz"
    genotype_df.to_csv(output_file, sep='\t', compression='gzip')
    print(f"  Genotypes saved to: {output_file}")

    # Save variant info
    output_file = f"{output_prefix}.variants.tsv"
    variant_stats_df.to_csv(output_file, sep='\t', index=False)
    print(f"  Variant info saved to: {output_file}")


def save_vcf(genotype_matrix, variant_stats_df, sample_ids, output_file, chromosome="1"):
    """
    Save genotypes in VCF format.

    Parameters:
    -----------
    genotype_matrix : np.ndarray
        Genotype matrix (samples x variants)
    variant_stats_df : pd.DataFrame
        Variant statistics
    sample_ids : list
        Sample IDs
    output_file : str
        Output VCF file path
    chromosome : str
        Chromosome name (default: "1")
    """
    print(f"\nSaving VCF...")

    with open(output_file, 'w') as f:
        # Write VCF header
        f.write("##fileformat=VCFv4.2\n")
        f.write("##source=msprime_simulation\n")
        f.write(f"##contig=<ID={chromosome}>\n")
        f.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n')
        f.write('##INFO=<ID=MAF,Number=1,Type=Float,Description="Minor Allele Frequency">\n')
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')

        # Write column header
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
        for sample_id in sample_ids:
            f.write(f"\t{sample_id}")
        f.write("\n")

        # Write variant records
        for idx, row in variant_stats_df.iterrows():
            var_id = row['variant_id']
            pos = row['position']
            af = row['af']
            maf = row['maf']

            # Get genotypes for this variant
            genotypes = genotype_matrix[:, idx]

            # Write variant line
            info = f"AF={af:.4f};MAF={maf:.4f}"
            f.write(f"{chromosome}\t{pos}\t{var_id}\tA\tT\t.\tPASS\t{info}\tGT")

            # Write genotypes
            for gt in genotypes:
                if gt == 0:
                    f.write("\t0/0")
                elif gt == 1:
                    f.write("\t0/1")
                elif gt == 2:
                    f.write("\t1/1")
                else:
                    f.write("\t./.")  # Missing
            f.write("\n")

    print(f"  VCF saved to: {output_file}")


def main():
    import sys
    print("Python script started...", flush=True)
    sys.stdout.flush()

    parser = argparse.ArgumentParser(
        description='Simulate genotypes using msprime for power analysis'
    )
    parser.add_argument('--n_samples', type=int, default=5000,
                       help='Number of samples to simulate (default: 5000)')
    parser.add_argument('--n_genes', type=int, default=20,
                       help='Number of independent genes to simulate (default: 20)')
    parser.add_argument('--variants_per_gene', type=int, default=100,
                       help='Target variants per gene (default: 100)')
    parser.add_argument('--gene_length_min', type=int, default=30_000,
                       help='Minimum gene length in bp (default: 30000)')
    parser.add_argument('--gene_length_max', type=int, default=40_000,
                       help='Maximum gene length in bp (default: 40000)')
    parser.add_argument('--intergenic_spacing', type=int, default=1_000_000,
                       help='Space between genes in bp (default: 1000000)')
    parser.add_argument('--recombination_rate', type=float, default=1e-8,
                       help='Recombination rate per bp (default: 1e-8)')
    parser.add_argument('--mutation_rate', type=float, default=2e-8,
                       help='Mutation rate per bp (default: 2e-8)')
    parser.add_argument('--population_size', type=int, default=10_000,
                       help='Effective population size (default: 10000)')
    parser.add_argument('--causal_gene_fraction', type=float, default=0.1,
                       help='Fraction of genes containing causal variants (default: 0.1)')
    parser.add_argument('--causal_variants_per_gene', type=int, default=1,
                       help='Number of causal variants per causal gene (default: 1)')
    parser.add_argument('--pLoF_maf_min', type=float, default=0.0,
                       help='Minimum MAF for selecting pLoF (causal) variants (default: 0.0)')
    parser.add_argument('--pLoF_maf_max', type=float, default=1.0,
                       help='Maximum MAF for selecting pLoF (causal) variants (default: 1.0)')
    parser.add_argument('--h2', type=float, default=0.1,
                       help='Heritability for phenotype simulation (default: 0.1)')
    parser.add_argument('--architecture', type=str, default='additive',
                       choices=['additive', 'recessive', 'dominant'],
                       help='Genetic architecture (default: additive)')
    parser.add_argument('--seed', type=int, default=42,
                       help='Random seed (default: 42)')
    parser.add_argument('--output_dir', type=str, default='input',
                       help='Output directory (default: input)')
    parser.add_argument('--output_prefix', type=str, default='simulated',
                       help='Output file prefix (default: simulated)')
    parser.add_argument('--chromosome', type=str, default='1',
                       help='Chromosome name for VCF (default: 1)')
    parser.add_argument('--force', action='store_true',
                       help='Force re-simulation even if genotypes already exist')

    args = parser.parse_args()

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("Genotype Simulation for Power Analysis")
    print("=" * 60)

    # Check if genotypes already exist
    output_prefix = output_dir / args.output_prefix
    genotype_file = output_dir / f'{args.output_prefix}.genotypes.tsv.gz'
    variant_file = output_dir / f'{args.output_prefix}.variants.tsv'
    annotation_file = output_dir / 'variant_annotations.tsv'
    gene_file = output_dir / 'gene_boundaries.tsv'
    params_file = output_dir / 'simulation_params.json'

    # Current genotype simulation parameters (excludes h2 and architecture - those are phenotype-only)
    current_params = {
        'n_samples': args.n_samples,
        'n_genes': args.n_genes,
        'variants_per_gene': args.variants_per_gene,
        'gene_length_min': args.gene_length_min,
        'gene_length_max': args.gene_length_max,
        'intergenic_spacing': args.intergenic_spacing,
        'recombination_rate': args.recombination_rate,
        'mutation_rate': args.mutation_rate,
        'population_size': args.population_size,
        'causal_gene_fraction': args.causal_gene_fraction,
        'causal_variants_per_gene': args.causal_variants_per_gene,
        'pLoF_maf_min': args.pLoF_maf_min,
        'pLoF_maf_max': args.pLoF_maf_max,
        'seed': args.seed,
    }

    # Check if parameters match saved parameters
    params_match = False
    if params_file.exists():
        import json
        with open(params_file, 'r') as f:
            saved_params = json.load(f)
        params_match = (saved_params == current_params)

    skip_genotype_sim = (
        not args.force and
        genotype_file.exists() and
        variant_file.exists() and
        annotation_file.exists() and
        gene_file.exists() and
        params_match
    )

    if skip_genotype_sim:
        print("\n✓ Genotype files already exist with matching parameters.", flush=True)
        print("  Skipping genotype simulation.", flush=True)
        print("  (Use --force to re-simulate genotypes)", flush=True)
        print("\nLoading existing genotype data...", flush=True)

        # Load existing data
        import gzip
        filtered_stats = pd.read_csv(variant_file, sep='\t')
        annotation_df = pd.read_csv(annotation_file, sep='\t')
        gene_df = pd.read_csv(gene_file, sep='\t')

        # Load genotypes
        genotype_df = pd.read_csv(genotype_file, sep='\t', index_col=0)
        filtered_genotypes = genotype_df.values

        print(f"  Loaded {len(filtered_stats)} variants", flush=True)
        print(f"  Loaded {filtered_genotypes.shape[0]} samples", flush=True)
        print(f"  Loaded {len(gene_df)} genes", flush=True)

    else:
        if args.force:
            print("\n⚠ Force mode enabled - re-simulating all genotypes", flush=True)
        elif params_file.exists() and not params_match:
            import json
            print("\n⚠ Genotype simulation parameters changed - re-simulating genotypes", flush=True)
            print("  Changed parameters:", flush=True)
            with open(params_file, 'r') as f:
                saved_params = json.load(f)
            for key in current_params:
                if current_params[key] != saved_params.get(key):
                    print(f"    {key}: {saved_params.get(key)} → {current_params[key]}", flush=True)

        # Step 1: Simulate independent genes
        filtered_stats, filtered_genotypes, gene_df = simulate_independent_genes(
            n_genes=args.n_genes,
            n_samples=args.n_samples,
            variants_per_gene=args.variants_per_gene,
            gene_length_min=args.gene_length_min,
            gene_length_max=args.gene_length_max,
            intergenic_spacing=args.intergenic_spacing,
            recombination_rate=args.recombination_rate,
            mutation_rate=args.mutation_rate,
            population_size=args.population_size,
            seed=args.seed,
            chromosome=args.chromosome
        )

        # Save gene boundaries
        gene_df.to_csv(output_dir / 'gene_boundaries.tsv', sep='\t', index=False)
        print(f"\n  Gene boundaries saved", flush=True)

        # Step 2: Create variant annotations
        annotation_df = create_variant_annotations(
            filtered_stats,
            gene_df,
            causal_gene_fraction=args.causal_gene_fraction,
            causal_variants_per_gene=args.causal_variants_per_gene,
            pLoF_maf_min=args.pLoF_maf_min,
            pLoF_maf_max=args.pLoF_maf_max,
            seed=args.seed,
            output_file=output_dir / 'variant_annotations.tsv'
        )

        # Save genotypes (only if we just simulated them)
        sample_ids = [f'sample_{i}' for i in range(args.n_samples)]
        save_genotypes(filtered_genotypes, filtered_stats, sample_ids, output_prefix)

        # Save VCF
        vcf_file = output_dir / f'{args.output_prefix}.vcf'
        save_vcf(filtered_genotypes, filtered_stats, sample_ids, vcf_file, chromosome=args.chromosome)

        # Save simulation parameters for future reference
        import json
        with open(params_file, 'w') as f:
            json.dump(current_params, f, indent=2)
        print(f"\n  Simulation parameters saved to: {params_file}", flush=True)

    # Step 6: Simulate phenotypes (always regenerate - it's fast)
    print("\n" + "=" * 60, flush=True)
    print("Regenerating phenotypes with current parameters", flush=True)
    print("=" * 60, flush=True)

    # IMPORTANT: Ensure annotation_df index matches genotype matrix columns
    # Reset index to ensure it's 0, 1, 2, ... matching genotype column order
    annotation_df = annotation_df.reset_index(drop=True)

    # Debug: verify alignment
    print(f"\nDebug - verifying data alignment:", flush=True)
    print(f"  Genotype matrix shape: {filtered_genotypes.shape}", flush=True)
    print(f"  Annotation rows: {len(annotation_df)}", flush=True)
    print(f"  Causal variants: {annotation_df['is_causal'].sum()}", flush=True)
    causal_genes = annotation_df[annotation_df['is_causal']]['gene_id'].unique()
    print(f"  Causal genes: {sorted(causal_genes)}", flush=True)

    sample_ids = [f'sample_{i}' for i in range(args.n_samples)]
    phenotype_df = simulate_phenotypes(
        filtered_genotypes,
        annotation_df,
        h2_total=args.h2,
        architecture=args.architecture,
        seed=args.seed
    )

    # Save phenotypes
    phenotype_file = output_dir / f'{args.output_prefix}.phenos.tsv'
    phenotype_df.to_csv(phenotype_file, sep='\t', index=False)
    print(f"\n  Phenotypes saved to: {phenotype_file}", flush=True)

    print("\n" + "=" * 60)
    print("Simulation completed successfully!")
    print("=" * 60)
    print(f"\nOutput files:")
    print(f"  - Genotypes (TSV): {genotype_file}")
    if not skip_genotype_sim:
        print(f"  - Genotypes (VCF): {vcf_file}")
    print(f"  - Variant info: {variant_file}")
    print(f"  - Phenotypes: {phenotype_file} ← UPDATED")
    print(f"  - Gene boundaries: {gene_file}")
    print(f"  - Variant annotations: {annotation_file}")
    print(f"  - Parameters: {params_file}")


if __name__ == '__main__':
    main()
