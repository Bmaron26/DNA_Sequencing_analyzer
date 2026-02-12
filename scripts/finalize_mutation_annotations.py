#!/usr/bin/env python3
"""
Final mutation annotation update and filtering.

This script:
1. Loads mutation data from all result directories
2. Updates gene names using NCBI gene mapping
3. Filters out unwanted mutations (synonymous, non-coding, ancestral)
4. Creates comprehensive output tables

Usage:
    python finalize_mutation_annotations.py \
        --results <dir1> <dir2> <dir3> \
        --gene-mapping <gene_mapping_complete.csv> \
        --output <output_dir> \
        --exclude-positions 1994565 \
        --exclude-genes agrC
"""

import argparse
import re
from pathlib import Path
from typing import Dict, List, Set, Optional
from collections import defaultdict

import pandas as pd
import numpy as np


def find_sample_mutations(results_dir: Path) -> List[tuple]:
    """Find all FreeBayes mutation CSV files, including nested directories."""
    samples = []

    for item in results_dir.iterdir():
        if not item.is_dir():
            continue

        sample_name = item.name

        # Direct sample directory
        mut_file = item / f"{sample_name}_mutations_freebayes.csv"
        if mut_file.exists():
            samples.append((sample_name, mut_file))
        else:
            # Check for nested directories
            for subitem in item.iterdir():
                if subitem.is_dir():
                    sub_name = subitem.name
                    sub_mut_file = subitem / f"{sub_name}_mutations_freebayes.csv"
                    if sub_mut_file.exists():
                        samples.append((sub_name, sub_mut_file))

    return samples


def extract_treatment(sample_name: str) -> str:
    """Extract treatment from sample name."""
    if '_' in sample_name:
        parts = sample_name.rsplit('_', 1)
        if len(parts) == 2 and parts[1] and parts[1][-1].isdigit():
            return re.sub(r'\d+$', '', sample_name)
    return re.sub(r'\d+$', '', sample_name)


def is_synonymous(effect: str) -> bool:
    """Check if mutation effect is synonymous."""
    if pd.isna(effect):
        return False
    effect_lower = str(effect).lower()
    return any(term in effect_lower for term in ['synonymous', 'silent', 'syn_coding'])


def is_coding_region(location_type: str) -> bool:
    """Check if mutation is in coding region."""
    if pd.isna(location_type):
        return True  # If unknown, assume coding
    loc_lower = str(location_type).lower().strip()
    return loc_lower == 'coding'


def load_gene_mapping(mapping_file: Path) -> Dict[str, dict]:
    """Load gene mapping from CSV file."""
    df = pd.read_csv(mapping_file)

    # Create lookup by locus_tag
    locus_to_info = {}

    for _, row in df.iterrows():
        locus_tag = row.get('locus_tag', '')
        if not locus_tag or pd.isna(locus_tag):
            continue

        locus_to_info[locus_tag] = {
            'ncbi_gene_name': row.get('ncbi_gene_name', row.get('gene_best', '')),
            'refseq_id': row.get('refseq_primary', row.get('refseq_id', '')),
            'ncbi_description': row.get('ncbi_description', row.get('clean_description', '')),
            'gene_best': row.get('gene_best', row.get('ncbi_gene_name', ''))
        }

    return locus_to_info


def main():
    parser = argparse.ArgumentParser(
        description='Finalize mutation annotations and create filtered tables'
    )
    parser.add_argument('--results', nargs='+', required=True,
                        help='Result directories to process')
    parser.add_argument('--labels', nargs='+', default=None,
                        help='Labels for each directory (Control, Single, Combination)')
    parser.add_argument('--gene-mapping', '-g', required=True,
                        help='Gene mapping CSV from NCBI parsing')
    parser.add_argument('--output', '-o', required=True,
                        help='Output directory')
    parser.add_argument('--exclude-genes', nargs='+', default=['agrC'],
                        help='Gene names to exclude (default: agrC)')
    parser.add_argument('--exclude-positions', nargs='+', type=int, default=[1994565],
                        help='Positions to exclude (default: 1994565)')

    args = parser.parse_args()

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    exclude_genes = set(g.lower() for g in args.exclude_genes)
    exclude_positions = set(args.exclude_positions)

    print("="*70)
    print("FINALIZE MUTATION ANNOTATIONS")
    print("="*70)

    # Load gene mapping
    print(f"\n1. Loading gene mapping: {args.gene_mapping}")
    gene_mapping = load_gene_mapping(Path(args.gene_mapping))
    print(f"   Loaded mapping for {len(gene_mapping)} genes")

    # Set up labels
    labels = args.labels if args.labels else ['Unknown'] * len(args.results)
    if len(labels) < len(args.results):
        labels.extend(['Unknown'] * (len(args.results) - len(labels)))

    # Load all mutations
    print("\n2. Loading mutations from all directories...")
    all_mutations = []

    for results_dir, label in zip(args.results, labels):
        results_path = Path(results_dir)
        if not results_path.exists():
            print(f"   Warning: {results_path} not found")
            continue

        # Determine treatment type
        dir_name = results_path.name.lower()
        if 'npsa' in dir_name or 'control' in label.lower():
            treatment_type = 'Control'
        elif 'combination' in dir_name or 'comb' in dir_name or 'combination' in label.lower():
            treatment_type = 'Combination'
        else:
            treatment_type = 'Single'

        samples = find_sample_mutations(results_path)
        print(f"   {results_path.name}: {len(samples)} samples ({treatment_type})")

        for sample_name, mut_file in samples:
            df = pd.read_csv(mut_file)

            if df.empty:
                continue

            # Add sample info
            df['sample'] = sample_name
            df['treatment'] = extract_treatment(sample_name)
            df['treatment_type'] = treatment_type

            all_mutations.append(df)

    if not all_mutations:
        print("\nNo mutations found!")
        return

    # Combine all mutations
    print("\n3. Combining all mutations...")
    combined_df = pd.concat(all_mutations, ignore_index=True)
    print(f"   Total raw mutations: {len(combined_df)}")

    # Detect column names
    pos_col = next((c for c in ['POS', 'position', 'pos'] if c in combined_df.columns), None)
    effect_col = next((c for c in ['EFFECT', 'effect', 'Effect'] if c in combined_df.columns), None)
    location_col = next((c for c in ['LOCATION', 'location_type', 'location'] if c in combined_df.columns), None)
    locus_col = next((c for c in ['LOCUS_TAG', 'locus_tag', 'gene_id'] if c in combined_df.columns), None)
    gene_col = next((c for c in ['GENE', 'gene_name', 'gene'] if c in combined_df.columns), None)
    product_col = next((c for c in ['PRODUCT', 'product'] if c in combined_df.columns), None)

    # Update gene names from mapping
    print("\n4. Updating gene names from NCBI mapping...")

    ncbi_gene_names = []
    ncbi_descriptions = []
    refseq_ids = []
    gene_best_names = []

    for _, row in combined_df.iterrows():
        locus = row.get(locus_col, '') if locus_col else ''

        if locus and str(locus) in gene_mapping:
            info = gene_mapping[str(locus)]
            ncbi_gene_names.append(info.get('ncbi_gene_name', ''))
            ncbi_descriptions.append(info.get('ncbi_description', ''))
            refseq_ids.append(info.get('refseq_id', ''))
            gene_best_names.append(info.get('gene_best', ''))
        else:
            # Fall back to original annotation
            orig_gene = row.get(gene_col, '') if gene_col else ''
            ncbi_gene_names.append('')
            ncbi_descriptions.append('')
            refseq_ids.append('')
            # Use original if available
            if orig_gene and not pd.isna(orig_gene):
                gene_best_names.append(str(orig_gene))
            else:
                gene_best_names.append(str(locus) if locus else 'unknown')

    combined_df['ncbi_gene'] = ncbi_gene_names
    combined_df['ncbi_description'] = ncbi_descriptions
    combined_df['refseq_id'] = refseq_ids
    combined_df['gene_final'] = gene_best_names

    # Override gene_final with ncbi_gene if available
    combined_df['gene_final'] = combined_df.apply(
        lambda row: row['ncbi_gene'] if row['ncbi_gene'] and str(row['ncbi_gene']).strip()
                    else row['gene_final'],
        axis=1
    )

    n_updated = (combined_df['ncbi_gene'] != '').sum()
    print(f"   Updated {n_updated} / {len(combined_df)} mutations with NCBI gene names")

    # Save full annotated table (before filtering)
    full_file = output_dir / 'all_mutations_annotated.csv'
    combined_df.to_csv(full_file, index=False)
    print(f"\n   Saved full annotated table: {full_file}")

    # Apply filters
    print("\n5. Applying filters...")
    print(f"   - Excluding synonymous mutations")
    print(f"   - Excluding non-coding mutations (keeping only location_type='coding')")
    print(f"   - Excluding genes: {args.exclude_genes}")
    print(f"   - Excluding positions: {args.exclude_positions}")

    initial_count = len(combined_df)
    filter_stats = {'initial': initial_count}

    # Filter 1: Synonymous mutations
    if effect_col:
        mask_syn = combined_df[effect_col].apply(is_synonymous)
        n_syn = mask_syn.sum()
        combined_df = combined_df[~mask_syn]
        filter_stats['synonymous_removed'] = n_syn
        print(f"      Removed {n_syn} synonymous mutations")

    # Filter 2: Non-coding mutations
    if location_col:
        mask_noncoding = ~combined_df[location_col].apply(is_coding_region)
        n_noncoding = mask_noncoding.sum()
        combined_df = combined_df[~mask_noncoding]
        filter_stats['noncoding_removed'] = n_noncoding
        print(f"      Removed {n_noncoding} non-coding mutations")

    # Filter 3: Excluded genes
    def should_exclude_gene(row):
        gene_final = str(row.get('gene_final', '')).lower()
        ncbi_gene = str(row.get('ncbi_gene', '')).lower()
        orig_gene = str(row.get(gene_col, '')).lower() if gene_col else ''

        for excl in exclude_genes:
            if excl in gene_final or excl in ncbi_gene or excl in orig_gene:
                return True
        return False

    mask_excluded_genes = combined_df.apply(should_exclude_gene, axis=1)
    excluded_genes_df = combined_df[mask_excluded_genes]
    n_excluded_genes = mask_excluded_genes.sum()
    combined_df = combined_df[~mask_excluded_genes]
    filter_stats['excluded_genes_removed'] = n_excluded_genes
    print(f"      Removed {n_excluded_genes} mutations in excluded genes")

    # Filter 4: Excluded positions
    if pos_col and exclude_positions:
        mask_excluded_pos = combined_df[pos_col].isin(exclude_positions)
        n_excluded_pos = mask_excluded_pos.sum()
        combined_df = combined_df[~mask_excluded_pos]
        filter_stats['excluded_positions_removed'] = n_excluded_pos
        print(f"      Removed {n_excluded_pos} mutations at excluded positions")

    filter_stats['final'] = len(combined_df)
    total_removed = initial_count - len(combined_df)
    print(f"\n   Total: {initial_count} -> {len(combined_df)} ({total_removed} removed, {100*total_removed/initial_count:.1f}%)")

    # Save filtered table
    filtered_file = output_dir / 'mutations_filtered.csv'
    combined_df.to_csv(filtered_file, index=False)
    print(f"\n   Saved filtered table: {filtered_file}")

    # Save excluded mutations for reference
    if not excluded_genes_df.empty:
        excluded_file = output_dir / 'mutations_excluded_genes.csv'
        excluded_genes_df.to_csv(excluded_file, index=False)
        print(f"   Saved excluded gene mutations: {excluded_file}")

    # Save filter statistics
    stats_df = pd.DataFrame([filter_stats])
    stats_file = output_dir / 'filter_statistics.csv'
    stats_df.to_csv(stats_file, index=False)

    # Create summary tables
    print("\n6. Creating summary tables...")

    # Summary by treatment
    treatment_summary = combined_df.groupby('treatment').agg({
        pos_col: 'count' if pos_col else 'size',
        'sample': 'nunique'
    }).reset_index()
    treatment_summary.columns = ['treatment', 'n_mutations', 'n_samples']
    treatment_summary['mutations_per_sample'] = (treatment_summary['n_mutations'] /
                                                  treatment_summary['n_samples']).round(2)
    treatment_summary.to_csv(output_dir / 'summary_by_treatment.csv', index=False)
    print(f"   Saved: summary_by_treatment.csv")

    # Summary by gene
    gene_summary = combined_df.groupby('gene_final').agg({
        pos_col: 'count' if pos_col else 'size',
        'sample': 'nunique',
        'treatment': lambda x: ', '.join(sorted(set(x)))
    }).reset_index()
    gene_summary.columns = ['gene', 'n_mutations', 'n_samples', 'treatments']
    gene_summary = gene_summary.sort_values('n_mutations', ascending=False)
    gene_summary.to_csv(output_dir / 'summary_by_gene.csv', index=False)
    print(f"   Saved: summary_by_gene.csv")

    # Gene x Treatment matrix
    gene_treatment = combined_df.groupby(['gene_final', 'treatment']).size().unstack(fill_value=0)
    gene_treatment['total'] = gene_treatment.sum(axis=1)
    gene_treatment = gene_treatment.sort_values('total', ascending=False)
    gene_treatment.to_csv(output_dir / 'gene_treatment_matrix.csv')
    print(f"   Saved: gene_treatment_matrix.csv")

    # Print top mutated genes
    print("\n7. Top 20 mutated genes (filtered):")
    print("-"*60)
    for _, row in gene_summary.head(20).iterrows():
        print(f"   {row['gene']:20s} {row['n_mutations']:4d} mutations in {row['n_samples']:2d} samples")

    print("\n" + "="*70)
    print("COMPLETE")
    print("="*70)
    print(f"\nOutput files in: {output_dir}")
    print("\n  Full data:")
    print("    - all_mutations_annotated.csv    : All mutations with NCBI gene names")
    print("\n  Filtered data (recommended for analysis):")
    print("    - mutations_filtered.csv         : Filtered (no syn/non-coding/ancestral)")
    print("    - mutations_excluded_genes.csv   : Mutations that were excluded")
    print("\n  Summary tables:")
    print("    - summary_by_treatment.csv       : Mutation counts by treatment")
    print("    - summary_by_gene.csv            : Mutation counts by gene")
    print("    - gene_treatment_matrix.csv      : Gene × Treatment matrix")
    print("    - filter_statistics.csv          : Filtering statistics")


if __name__ == '__main__':
    main()
