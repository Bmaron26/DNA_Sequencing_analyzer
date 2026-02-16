#!/usr/bin/env python3
"""
Annotate and filter mutation summary table.

This script:
1. Reads an existing mutation summary CSV (e.g., all_mutations_summary_freebayes.csv)
2. Updates gene names using a gene_lookup_table.csv (by locus_tag)
3. Creates filtered versions:
   - Full annotated summary (with gene names resolved)
   - Filtered summary (excluding specified genes like agrC, and excluding NPSA strains)
   - NPSA-only summary (separate table for control samples)

Usage:
    python annotate_and_filter_summary.py \
        --input /mnt/c/Users/hayouka-lab/Documents/Bar/3rd_evolution/WGS/results_HC_freebayes/all_mutations_summary_freebayes.csv \
        --gene-lookup /mnt/c/Users/hayouka-lab/Documents/Bar/3rd_evolution/WGS/mutation_analysis_visualization_new/gene_lookup_table.csv \
        --output-dir /mnt/c/Users/hayouka-lab/Documents/Bar/3rd_evolution/WGS/results_HC_freebayes \
        --exclude-genes agrC \
        --exclude-strains NPSA
"""

import argparse
import re
import sys
from pathlib import Path
from typing import Dict, List, Set

import pandas as pd


def load_gene_lookup(lookup_path: Path) -> Dict[str, dict]:
    """
    Load gene lookup table keyed by locus_tag.

    Expected columns: locus_tag, product, refseq_id, gene_final, gene_alternative, gene_name
    """
    df = pd.read_csv(lookup_path, sep=None, engine='python')  # auto-detect delimiter

    lookup = {}
    for _, row in df.iterrows():
        locus_tag = str(row.get('locus_tag', '')).strip()
        if not locus_tag or locus_tag == 'nan':
            continue

        lookup[locus_tag] = {
            'gene_final': str(row.get('gene_final', '')).strip(),
            'gene_name': str(row.get('gene_name', '')).strip(),
            'gene_alternative': str(row.get('gene_alternative', '')).strip(),
            'product_lookup': str(row.get('product', '')).strip(),
            'refseq_id_lookup': str(row.get('refseq_id', '')).strip(),
        }

    return lookup


def resolve_gene_name(row: pd.Series, locus_col: str, gene_col: str,
                      lookup: Dict[str, dict]) -> str:
    """Resolve the best gene name using the lookup table."""
    locus_tag = ''
    if locus_col and locus_col in row.index:
        locus_tag = str(row[locus_col]).strip()
        if locus_tag == 'nan':
            locus_tag = ''

    # Try lookup by locus_tag
    if locus_tag and locus_tag in lookup:
        info = lookup[locus_tag]
        if info['gene_final'] and info['gene_final'] != 'nan':
            return info['gene_final']
        if info['gene_name'] and info['gene_name'] != 'nan':
            return info['gene_name']

    # Fallback to existing gene column
    if gene_col and gene_col in row.index:
        gene = str(row[gene_col]).strip()
        if gene and gene != 'nan':
            return gene

    # Fallback to locus_tag
    if locus_tag:
        return locus_tag

    return ''


def should_exclude_gene(gene_name: str, exclude_genes: Set[str]) -> bool:
    """Check if a gene should be excluded (case-insensitive partial match)."""
    if not gene_name:
        return False
    gene_lower = gene_name.lower()
    for excl in exclude_genes:
        if excl in gene_lower:
            return True
    return False


def should_exclude_strain(sample_name: str, exclude_prefixes: List[str]) -> bool:
    """Check if a sample should be excluded based on strain prefix."""
    if not sample_name:
        return False
    name_upper = str(sample_name).upper()
    return any(name_upper.startswith(prefix.upper()) for prefix in exclude_prefixes)


def main():
    parser = argparse.ArgumentParser(
        description='Annotate and filter mutation summary table'
    )
    parser.add_argument('--input', '-i', required=True,
                        help='Input mutation summary CSV file')
    parser.add_argument('--gene-lookup', '-g', required=True,
                        help='Path to gene_lookup_table.csv')
    parser.add_argument('--output-dir', '-o', default=None,
                        help='Output directory (defaults to input file directory)')
    parser.add_argument('--exclude-genes', nargs='+', default=['agrC'],
                        help='Gene names to exclude (case-insensitive, default: agrC)')
    parser.add_argument('--exclude-strains', nargs='+', default=['NPSA'],
                        help='Strain prefixes to exclude (default: NPSA)')

    args = parser.parse_args()

    input_file = Path(args.input)
    output_dir = Path(args.output_dir) if args.output_dir else input_file.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    exclude_genes = set(g.lower() for g in args.exclude_genes)
    exclude_strains = args.exclude_strains

    print("=" * 70)
    print("ANNOTATE AND FILTER MUTATION SUMMARY")
    print("=" * 70)

    # 1. Load gene lookup table
    lookup_path = Path(args.gene_lookup)
    print(f"\n1. Loading gene lookup table: {lookup_path}")
    gene_lookup = load_gene_lookup(lookup_path)
    print(f"   Loaded {len(gene_lookup)} locus_tag entries")

    # 2. Load input summary
    print(f"\n2. Loading mutation summary: {input_file}")
    df = pd.read_csv(input_file)
    print(f"   Loaded {len(df)} mutations")
    print(f"   Columns: {', '.join(df.columns)}")

    # Check for sample column
    sample_col = next((c for c in ['sample', 'Sample', 'SAMPLE'] if c in df.columns), None)
    if sample_col:
        samples = df[sample_col].unique()
        print(f"   Samples found: {', '.join(str(s) for s in samples)}")
    else:
        print("   WARNING: No 'sample' column found!")

    # 3. Detect column names
    locus_col = next((c for c in ['locus_tag', 'LOCUS_TAG', 'gene_id'] if c in df.columns), None)
    gene_col = next((c for c in ['gene_name', 'gene', 'GENE'] if c in df.columns), None)

    print(f"\n3. Detected columns:")
    print(f"   Locus tag column: {locus_col}")
    print(f"   Gene name column: {gene_col}")
    print(f"   Sample column: {sample_col}")

    # 4. Resolve gene names using lookup table
    print(f"\n4. Resolving gene names from lookup table...")

    df['gene_name_resolved'] = df.apply(
        lambda row: resolve_gene_name(row, locus_col, gene_col, gene_lookup),
        axis=1
    )

    # Count how many were resolved via lookup
    n_total = len(df)
    n_resolved = (df['gene_name_resolved'] != '').sum()
    n_from_lookup = 0
    if locus_col:
        for _, row in df.iterrows():
            lt = str(row.get(locus_col, '')).strip()
            if lt != 'nan' and lt in gene_lookup:
                n_from_lookup += 1

    print(f"   Gene names resolved: {n_resolved} / {n_total}")
    print(f"   Matched via lookup table: {n_from_lookup}")

    # 5. Reorder columns - put sample and gene_name_resolved first
    cols = list(df.columns)
    priority_cols = []
    if sample_col and sample_col in cols:
        cols.remove(sample_col)
        priority_cols.append(sample_col)
    if 'gene_name_resolved' in cols:
        cols.remove('gene_name_resolved')
        priority_cols.append('gene_name_resolved')

    df = df[priority_cols + cols]

    # 6. Save full annotated summary
    full_output = output_dir / 'mutations_annotated_full.csv'
    df.to_csv(full_output, index=False)
    print(f"\n5. Saved full annotated summary: {full_output}")
    print(f"   ({len(df)} mutations)")

    # 7. Create filtered version (exclude specific genes and strains)
    print(f"\n6. Creating filtered summaries...")
    print(f"   Excluding genes: {', '.join(args.exclude_genes)}")
    print(f"   Excluding strains: {', '.join(exclude_strains)}")

    # Mark mutations for exclusion
    df['_exclude_gene'] = df['gene_name_resolved'].apply(
        lambda g: should_exclude_gene(str(g), exclude_genes)
    )

    if sample_col:
        df['_exclude_strain'] = df[sample_col].apply(
            lambda s: should_exclude_strain(str(s), exclude_strains)
        )
    else:
        df['_exclude_strain'] = False

    # Count exclusions
    n_excluded_gene = df['_exclude_gene'].sum()
    n_excluded_strain = df['_exclude_strain'].sum()

    print(f"\n   Mutations excluded by gene ({', '.join(args.exclude_genes)}): {n_excluded_gene}")
    print(f"   Mutations excluded by strain ({', '.join(exclude_strains)}): {n_excluded_strain}")

    # Create NPSA-only table (control samples)
    if sample_col:
        npsa_df = df[df['_exclude_strain']].copy()
        npsa_df = npsa_df.drop(columns=['_exclude_gene', '_exclude_strain'])

        npsa_output = output_dir / 'mutations_NPSA_only.csv'
        npsa_df.to_csv(npsa_output, index=False)
        print(f"\n   Saved NPSA-only summary: {npsa_output}")
        print(f"   ({len(npsa_df)} mutations from {npsa_df[sample_col].nunique()} NPSA samples)")

        # List NPSA samples
        if len(npsa_df) > 0:
            print(f"   NPSA samples:")
            for sample in sorted(npsa_df[sample_col].unique()):
                count = len(npsa_df[npsa_df[sample_col] == sample])
                print(f"      {sample}: {count} mutations")

    # Create filtered table (exclude both genes and strains)
    filtered_df = df[~df['_exclude_gene'] & ~df['_exclude_strain']].copy()
    filtered_df = filtered_df.drop(columns=['_exclude_gene', '_exclude_strain'])

    filtered_output = output_dir / 'mutations_filtered.csv'
    filtered_df.to_csv(filtered_output, index=False)
    print(f"\n   Saved filtered summary: {filtered_output}")
    print(f"   ({len(filtered_df)} mutations, excluding {args.exclude_genes} and {exclude_strains})")

    # Also create a version that only excludes genes but keeps all strains
    genes_only_filtered = df[~df['_exclude_gene']].copy()
    genes_only_filtered = genes_only_filtered.drop(columns=['_exclude_gene', '_exclude_strain'])

    genes_filtered_output = output_dir / 'mutations_filtered_genes_only.csv'
    genes_only_filtered.to_csv(genes_filtered_output, index=False)
    print(f"\n   Saved gene-filtered summary (all strains): {genes_filtered_output}")
    print(f"   ({len(genes_only_filtered)} mutations, excluding only {args.exclude_genes})")

    # 8. Print quick overview by sample
    print(f"\n7. Quick overview by sample:")
    print("-" * 70)
    print(f"{'Sample':<25} {'Total':>8} {'Filtered':>10} {'Excluded':>10}")
    print("-" * 70)

    if sample_col:
        for sample in sorted(df[sample_col].unique()):
            total = len(df[df[sample_col] == sample])
            is_npsa = should_exclude_strain(sample, exclude_strains)

            if is_npsa:
                filtered = 0
                excluded = total
                marker = " [NPSA - separate table]"
            else:
                sample_df = df[df[sample_col] == sample]
                excluded = sample_df['_exclude_gene'].sum()
                filtered = total - excluded
                marker = ""

            print(f"   {sample:<25} {total:>6} {filtered:>10} {excluded:>10}{marker}")

    print("-" * 70)

    # 9. Show excluded gene details
    excluded_genes_df = df[df['_exclude_gene'] & ~df['_exclude_strain']]
    if len(excluded_genes_df) > 0:
        print(f"\n8. Excluded gene mutations ({', '.join(args.exclude_genes)}) detail:")
        print("-" * 70)

        # Save excluded mutations to separate file
        excluded_output = output_dir / 'mutations_excluded_genes.csv'
        excluded_save = excluded_genes_df.drop(columns=['_exclude_gene', '_exclude_strain'])
        excluded_save.to_csv(excluded_output, index=False)
        print(f"   Saved to: {excluded_output}")

        if sample_col:
            for sample in sorted(excluded_genes_df[sample_col].unique()):
                sample_excluded = excluded_genes_df[excluded_genes_df[sample_col] == sample]
                genes = sample_excluded['gene_name_resolved'].unique()
                print(f"   {sample}: {len(sample_excluded)} mutations in {', '.join(genes)}")

    print(f"\n" + "=" * 70)
    print("DONE!")
    print("=" * 70)
    print(f"\nOutput files in: {output_dir}")
    print(f"\n  mutations_annotated_full.csv      : All mutations with gene names resolved")
    print(f"  mutations_filtered.csv            : Filtered (no {'/'.join(args.exclude_genes)}, no {'/'.join(exclude_strains)})")
    print(f"  mutations_filtered_genes_only.csv : Filtered genes only (all strains)")
    print(f"  mutations_NPSA_only.csv           : NPSA control samples only")
    print(f"  mutations_excluded_genes.csv      : Excluded gene mutations (for reference)")


if __name__ == '__main__':
    main()
