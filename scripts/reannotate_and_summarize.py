#!/usr/bin/env python3
"""
Re-annotate mutations using a better-annotated reference genome and create summaries.

This script:
1. Parses a well-annotated GFF file (e.g., NCTC8325) to extract gene names
2. Re-annotates existing FreeBayes mutation CSV files
3. Combines mutations from multiple result directories
4. Generates summary statistics

Usage:
    python reannotate_and_summarize.py --gff <gff_file> --results <dir1> <dir2> ... --output <output_dir>
"""

import sys
import os
import re
import argparse
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
from collections import defaultdict
from dataclasses import dataclass

import pandas as pd
import numpy as np


@dataclass
class GeneInfo:
    """Store gene annotation information."""
    gene_name: str  # Short name like 'graS', 'vraG'
    locus_tag: str  # Like 'SAOUHSC_00006'
    product: str    # Full description
    start: int
    end: int
    strand: str
    feature_type: str  # CDS, gene, etc.


def parse_gff_attributes(attr_string: str) -> Dict[str, str]:
    """Parse GFF3 attribute string into dictionary."""
    attributes = {}
    if not attr_string or attr_string == '.':
        return attributes

    for item in attr_string.split(';'):
        item = item.strip()
        if '=' in item:
            key, value = item.split('=', 1)
            # URL decode common patterns
            value = value.replace('%2C', ',').replace('%3B', ';').replace('%25', '%')
            attributes[key] = value

    return attributes


def parse_gff_file(gff_path: Path) -> Tuple[Dict[int, GeneInfo], List[GeneInfo]]:
    """
    Parse GFF file and extract gene information.

    Returns:
        position_to_gene: Dict mapping each position to the gene it falls within
        all_genes: List of all gene records
    """
    genes = []

    print(f"   Parsing GFF file: {gff_path}")

    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            seqid, source, feature_type, start, end, score, strand, phase, attributes = parts

            # Focus on CDS and gene features
            if feature_type not in ['CDS', 'gene', 'mRNA', 'tRNA', 'rRNA', 'ncRNA']:
                continue

            attrs = parse_gff_attributes(attributes)

            # Extract gene name - try multiple attribute names
            gene_name = ''
            for key in ['gene', 'Name', 'gene_name', 'symbol']:
                if key in attrs:
                    gene_name = attrs[key]
                    break

            # Extract locus tag
            locus_tag = attrs.get('locus_tag', attrs.get('ID', ''))

            # Extract product
            product = attrs.get('product', attrs.get('note', ''))

            # If no gene name, try to extract from product or locus_tag
            if not gene_name and product:
                # Try to find gene name pattern in product (e.g., "GraS" in "sensor histidine kinase GraS")
                match = re.search(r'\b([A-Z][a-z]{2,3}[A-Z0-9]?)\b', product)
                if match:
                    gene_name = match.group(1).lower()

            if not gene_name:
                gene_name = locus_tag

            try:
                gene_info = GeneInfo(
                    gene_name=gene_name,
                    locus_tag=locus_tag,
                    product=product,
                    start=int(start),
                    end=int(end),
                    strand=strand,
                    feature_type=feature_type
                )
                genes.append(gene_info)
            except ValueError:
                continue

    print(f"   Found {len(genes)} gene features")

    # Build position lookup - prioritize CDS over gene
    # First, create a simple interval lookup
    cds_genes = [g for g in genes if g.feature_type == 'CDS']
    other_genes = [g for g in genes if g.feature_type != 'CDS']

    # For efficiency, we'll create ranges and do lookups on demand
    return genes, cds_genes + other_genes


def find_gene_at_position(position: int, genes: List[GeneInfo]) -> Optional[GeneInfo]:
    """Find the gene that contains a given position."""
    for gene in genes:
        if gene.start <= position <= gene.end:
            return gene
    return None


def extract_short_gene_name(gene_name: str, product: str, locus_tag: str) -> str:
    """
    Extract short gene name from available annotations.

    Priority:
    1. If gene_name looks like a proper gene name (3-4 letters), use it
    2. Try to extract from product description
    3. Fall back to locus_tag
    """
    # Clean up gene name
    if gene_name:
        gene_name = gene_name.strip()
        # Check if it's already a good short name (e.g., graS, vraG, mprF)
        if re.match(r'^[a-zA-Z]{2,4}[A-Z0-9]?$', gene_name):
            return gene_name
        # Handle names like "SAOUHSC_00006" - not useful
        if gene_name.startswith('SAOUHSC') or gene_name.startswith('SAO'):
            gene_name = ''

    # Try to extract from product
    if product:
        product = str(product)
        # Common patterns for gene names in product descriptions
        patterns = [
            r'\b([A-Z][a-z]{2,3}[A-Z0-9])\b',  # GraS, VraG, MprF
            r'\b([a-z]{3,4}[A-Z0-9]?)\b',       # graS, vraG
            r'(?:protein|subunit|factor)\s+([A-Z][a-z]*[A-Z0-9]*)',  # protein GraS
        ]
        for pattern in patterns:
            match = re.search(pattern, product)
            if match:
                name = match.group(1)
                # Validate it looks like a gene name
                if 2 <= len(name) <= 6 and not name.isupper():
                    return name

    # If we have a gene name that's not a locus tag, use it
    if gene_name and not gene_name.startswith('SAOUHSC'):
        return gene_name

    # Fall back to locus tag, but clean it up
    if locus_tag:
        # Remove common prefixes
        short = re.sub(r'^SAOUHSC_', '', locus_tag)
        short = re.sub(r'^SAO_', '', short)
        return short if short != locus_tag else locus_tag

    return gene_name or 'unknown'


def find_sample_mutations(results_dir: Path) -> List[Tuple[str, Path]]:
    """
    Find all FreeBayes mutation CSV files in a results directory.
    Handles nested directories (e.g., BmKn_comb/BmKn_Smp1).
    """
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


def reannotate_mutation_file(
    mut_file: Path,
    sample_name: str,
    genes: List[GeneInfo],
    treatment_type: str
) -> pd.DataFrame:
    """Re-annotate a mutation file with better gene names."""

    df = pd.read_csv(mut_file)

    if df.empty:
        return df

    # Detect position column
    pos_col = None
    for col in ['POS', 'position', 'pos', 'Position']:
        if col in df.columns:
            pos_col = col
            break

    if pos_col is None:
        print(f"   Warning: No position column found in {mut_file}")
        return df

    # Add new annotation columns
    new_gene_names = []
    new_short_names = []
    new_products = []
    new_locus_tags = []

    for _, row in df.iterrows():
        pos = row.get(pos_col)

        if pd.isna(pos):
            new_gene_names.append('')
            new_short_names.append('')
            new_products.append('')
            new_locus_tags.append('')
            continue

        gene = find_gene_at_position(int(pos), genes)

        if gene:
            new_gene_names.append(gene.gene_name)
            new_locus_tags.append(gene.locus_tag)
            new_products.append(gene.product)

            # Get short name
            short_name = extract_short_gene_name(
                gene.gene_name,
                gene.product,
                gene.locus_tag
            )
            new_short_names.append(short_name)
        else:
            # Use existing annotation if available
            existing_gene = row.get('GENE', row.get('gene_name', ''))
            existing_product = row.get('PRODUCT', row.get('product', ''))
            existing_locus = row.get('LOCUS_TAG', row.get('locus_tag', ''))

            new_gene_names.append(existing_gene if not pd.isna(existing_gene) else '')
            new_products.append(existing_product if not pd.isna(existing_product) else '')
            new_locus_tags.append(existing_locus if not pd.isna(existing_locus) else '')

            short_name = extract_short_gene_name(
                str(existing_gene) if not pd.isna(existing_gene) else '',
                str(existing_product) if not pd.isna(existing_product) else '',
                str(existing_locus) if not pd.isna(existing_locus) else ''
            )
            new_short_names.append(short_name)

    # Add new columns
    df['gene_short'] = new_short_names
    df['gene_name_new'] = new_gene_names
    df['product_new'] = new_products
    df['locus_tag_new'] = new_locus_tags

    # Add sample info
    df['sample'] = sample_name
    df['treatment_type'] = treatment_type

    # Extract treatment from sample name
    treatment = extract_treatment(sample_name)
    df['treatment'] = treatment

    return df


def extract_treatment(sample_name: str) -> str:
    """Extract treatment from sample name."""
    # Handle combinations like 'Mel_Pex1' -> 'Mel_Pex'
    # Handle singles like 'Mel1' -> 'Mel'
    # Handle NPSA like 'NPSA1' -> 'NPSA'

    if '_' in sample_name:
        parts = sample_name.rsplit('_', 1)
        if len(parts) == 2 and parts[1] and parts[1][-1].isdigit():
            return re.sub(r'\d+$', '', sample_name)
    return re.sub(r'\d+$', '', sample_name)


def create_summary_statistics(
    all_mutations: pd.DataFrame,
    output_dir: Path
):
    """Create summary statistics and visualizations."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    # Summary by treatment type
    print("\n   Creating summary statistics...")

    # Count mutations per sample
    sample_counts = all_mutations.groupby(['sample', 'treatment', 'treatment_type']).size().reset_index(name='n_mutations')
    sample_counts.to_csv(output_dir / 'mutations_per_sample.csv', index=False)

    # Summary by treatment
    treatment_summary = sample_counts.groupby('treatment').agg({
        'n_mutations': ['mean', 'std', 'min', 'max', 'count']
    }).round(2)
    treatment_summary.columns = ['mean', 'std', 'min', 'max', 'n_samples']
    treatment_summary.to_csv(output_dir / 'mutations_by_treatment.csv')

    # Summary by treatment type (Single vs Combination vs Control)
    type_summary = sample_counts.groupby('treatment_type').agg({
        'n_mutations': ['mean', 'std', 'min', 'max', 'sum', 'count']
    }).round(2)
    type_summary.columns = ['mean', 'std', 'min', 'max', 'total_mutations', 'n_samples']
    type_summary.to_csv(output_dir / 'mutations_by_type.csv')

    print(f"\n   Summary by treatment type:")
    print(type_summary.to_string())

    # Gene frequency analysis
    if 'gene_short' in all_mutations.columns:
        gene_freq = all_mutations.groupby('gene_short').agg({
            'sample': 'nunique',
            'treatment': lambda x: ', '.join(sorted(set(x)))
        }).reset_index()
        gene_freq.columns = ['gene', 'n_samples', 'treatments']
        gene_freq = gene_freq.sort_values('n_samples', ascending=False)
        gene_freq.to_csv(output_dir / 'gene_frequency.csv', index=False)

        print(f"\n   Top 20 most frequently mutated genes:")
        print(gene_freq.head(20).to_string(index=False))

    # Effect distribution
    effect_col = None
    for col in ['EFFECT', 'effect', 'Effect']:
        if col in all_mutations.columns:
            effect_col = col
            break

    if effect_col:
        effect_counts = all_mutations[effect_col].value_counts()
        effect_counts.to_csv(output_dir / 'effect_distribution.csv')

        # Plot effect distribution
        fig, ax = plt.subplots(figsize=(12, 6))
        effect_counts.head(15).plot(kind='bar', ax=ax, color='steelblue')
        ax.set_ylabel('Count')
        ax.set_xlabel('Effect')
        ax.set_title('Mutation Effect Distribution')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig(output_dir / 'effect_distribution.png', dpi=150)
        plt.close()

    # Create boxplot of mutations by treatment
    if not sample_counts.empty:
        fig, ax = plt.subplots(figsize=(14, 8))

        treatments = sorted(sample_counts['treatment'].unique())
        data = [sample_counts[sample_counts['treatment'] == t]['n_mutations'].values
                for t in treatments]

        bp = ax.boxplot(data, labels=treatments, patch_artist=True)

        # Color by treatment type
        colors = {'Single': '#729ECE', 'Combination': '#FF9E4A', 'Control': '#67BF5C'}
        for i, treatment in enumerate(treatments):
            treatment_type = sample_counts[sample_counts['treatment'] == treatment]['treatment_type'].iloc[0]
            bp['boxes'][i].set_facecolor(colors.get(treatment_type, '#888888'))

        ax.set_ylabel('Number of Mutations')
        ax.set_xlabel('Treatment')
        ax.set_title('Mutation Counts by Treatment')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig(output_dir / 'mutations_by_treatment_boxplot.png', dpi=150)
        plt.close()


def main():
    parser = argparse.ArgumentParser(
        description='Re-annotate mutations and create summary'
    )
    parser.add_argument('--gff', required=True, type=str,
                        help='Path to well-annotated GFF file')
    parser.add_argument('--results', nargs='+', required=True,
                        help='Result directories to process')
    parser.add_argument('--output', '-o', required=True, type=str,
                        help='Output directory')
    parser.add_argument('--labels', nargs='+', default=None,
                        help='Labels for each result directory (e.g., Control Single Combination)')
    parser.add_argument('--exclude-genes', nargs='+', default=[],
                        help='Gene names to exclude (e.g., agrC agrA). Case-insensitive.')
    parser.add_argument('--exclude-positions', nargs='+', type=int, default=[],
                        help='Specific positions to exclude')

    args = parser.parse_args()

    # Convert exclude genes to lowercase for case-insensitive matching
    exclude_genes = set(g.lower() for g in args.exclude_genes)

    gff_path = Path(args.gff)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("="*70)
    print("MUTATION RE-ANNOTATION AND SUMMARY")
    print("="*70)

    # Parse GFF file
    print("\n1. Parsing GFF annotation file...")
    all_genes, sorted_genes = parse_gff_file(gff_path)

    # Process each results directory
    print("\n2. Processing mutation files...")

    all_mutations = []

    labels = args.labels if args.labels else ['Unknown'] * len(args.results)
    if len(labels) < len(args.results):
        labels.extend(['Unknown'] * (len(args.results) - len(labels)))

    for results_dir, label in zip(args.results, labels):
        results_path = Path(results_dir)
        if not results_path.exists():
            print(f"   Warning: Directory not found: {results_path}")
            continue

        print(f"\n   Processing: {results_path.name} ({label})")

        # Determine treatment type from label or directory name
        dir_name = results_path.name.lower()
        if 'npsa' in dir_name or 'control' in label.lower():
            treatment_type = 'Control'
        elif 'combination' in dir_name or 'comb' in dir_name or 'combination' in label.lower():
            treatment_type = 'Combination'
        else:
            treatment_type = 'Single'

        samples = find_sample_mutations(results_path)
        print(f"   Found {len(samples)} samples")

        for sample_name, mut_file in samples:
            df = reannotate_mutation_file(mut_file, sample_name, sorted_genes, treatment_type)
            if not df.empty:
                all_mutations.append(df)

                # Save re-annotated file
                output_sample_dir = output_dir / 'reannotated' / treatment_type
                output_sample_dir.mkdir(parents=True, exist_ok=True)
                df.to_csv(output_sample_dir / f'{sample_name}_mutations_reannotated.csv', index=False)

    if not all_mutations:
        print("\nNo mutations found!")
        return

    # Combine all mutations
    print("\n3. Combining all mutations...")
    combined_df = pd.concat(all_mutations, ignore_index=True)

    # Filter out excluded genes
    initial_count = len(combined_df)
    if exclude_genes:
        print(f"\n   Filtering out excluded genes: {', '.join(args.exclude_genes)}")

        # Create mask for genes to exclude (case-insensitive)
        def should_exclude_gene(row):
            gene_short = str(row.get('gene_short', '')).lower()
            gene_name = str(row.get('gene_name_new', '')).lower()
            gene_orig = str(row.get('GENE', row.get('gene_name', ''))).lower()

            for excl in exclude_genes:
                if excl in gene_short or excl in gene_name or excl in gene_orig:
                    return True
            return False

        exclude_mask = combined_df.apply(should_exclude_gene, axis=1)
        excluded_df = combined_df[exclude_mask]
        combined_df = combined_df[~exclude_mask]

        # Save excluded mutations for reference
        if not excluded_df.empty:
            excluded_df.to_csv(output_dir / 'excluded_mutations.csv', index=False)
            print(f"   Excluded {len(excluded_df)} mutations (saved to excluded_mutations.csv)")

    # Filter out excluded positions
    if args.exclude_positions:
        print(f"\n   Filtering out excluded positions: {args.exclude_positions}")
        pos_col = None
        for col in ['POS', 'position', 'pos']:
            if col in combined_df.columns:
                pos_col = col
                break

        if pos_col:
            before = len(combined_df)
            combined_df = combined_df[~combined_df[pos_col].isin(args.exclude_positions)]
            print(f"   Excluded {before - len(combined_df)} mutations at specified positions")

    print(f"\n   Total mutations after filtering: {len(combined_df)} (removed {initial_count - len(combined_df)})")

    combined_df.to_csv(output_dir / 'all_mutations_combined.csv', index=False)
    print(f"   Total samples: {combined_df['sample'].nunique()}")

    # Create summary statistics
    print("\n4. Creating summary statistics...")
    create_summary_statistics(combined_df, output_dir)

    # Create gene-centric summary
    print("\n5. Creating gene-centric summary...")

    # Genes by treatment type
    gene_by_type = combined_df.groupby(['gene_short', 'treatment_type']).size().unstack(fill_value=0)
    gene_by_type['total'] = gene_by_type.sum(axis=1)
    gene_by_type = gene_by_type.sort_values('total', ascending=False)
    gene_by_type.to_csv(output_dir / 'genes_by_treatment_type.csv')

    # Genes by treatment
    gene_by_treatment = combined_df.groupby(['gene_short', 'treatment']).size().unstack(fill_value=0)
    gene_by_treatment['total'] = gene_by_treatment.sum(axis=1)
    gene_by_treatment = gene_by_treatment.sort_values('total', ascending=False)
    gene_by_treatment.to_csv(output_dir / 'genes_by_treatment.csv')

    print(f"\n" + "="*70)
    print("SUMMARY COMPLETE")
    print("="*70)
    print(f"\nOutput files saved to: {output_dir}")
    print("\nFiles created:")
    print("  - all_mutations_combined.csv      : All mutations with new annotations")
    print("  - mutations_per_sample.csv        : Mutation counts per sample")
    print("  - mutations_by_treatment.csv      : Summary statistics by treatment")
    print("  - mutations_by_type.csv           : Summary by treatment type")
    print("  - gene_frequency.csv              : Gene mutation frequency")
    print("  - genes_by_treatment_type.csv     : Gene × treatment type matrix")
    print("  - genes_by_treatment.csv          : Gene × treatment matrix")
    print("  - effect_distribution.csv/.png    : Effect type distribution")
    print("  - mutations_by_treatment_boxplot.png : Boxplot visualization")
    print("  - reannotated/                    : Individual re-annotated files")


if __name__ == '__main__':
    main()
