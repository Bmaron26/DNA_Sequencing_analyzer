#!/usr/bin/env python3
"""
Create mutation summary with sample identification and gene name lookup.

This script:
1. Reads all per-sample *_mutations.csv files from a results directory
2. Adds the sample name (strain) as a column derived from the filename
3. Resolves gene names using a gene_lookup_table.csv (by locus_tag)
4. Outputs:
   - all_mutations_summary_with_samples.csv  : Full summary with sample column
   - mutations_summary_filtered_no_NPSA.csv  : Filtered summary excluding NPSA strains

Usage:
    python create_hc_summary.py \
        --input-dir /mnt/c/Users/hayouka-lab/Documents/Bar/3rd_evolution/WGS/results_HC \
        --gene-lookup /mnt/c/Users/hayouka-lab/Documents/Bar/3rd_evolution/WGS/mutation_analysis_visualization_new/gene_lookup_table.csv \
        --output-dir /mnt/c/Users/hayouka-lab/Documents/Bar/3rd_evolution/WGS/results_HC
"""

import argparse
import re
import sys
from pathlib import Path
from typing import Dict, Optional

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


def extract_sample_name(filename: str) -> str:
    """
    Extract sample name from a mutations CSV filename.

    E.g. 'Mel_HC_3_mutations.csv' -> 'Mel_HC_3'
    """
    return re.sub(r'_mutations\.csv$', '', filename)


def resolve_gene_name(row: pd.Series, locus_col: Optional[str],
                      gene_col: Optional[str], lookup: Dict[str, dict]) -> str:
    """Resolve the best gene name using the lookup table, falling back to existing annotation."""
    locus_tag = ''
    if locus_col and locus_col in row.index:
        locus_tag = str(row[locus_col]).strip()
        if locus_tag == 'nan':
            locus_tag = ''

    # Try lookup
    if locus_tag and locus_tag in lookup:
        info = lookup[locus_tag]
        # Prefer gene_final from lookup table
        if info['gene_final'] and info['gene_final'] != 'nan':
            return info['gene_final']
        if info['gene_name'] and info['gene_name'] != 'nan':
            return info['gene_name']

    # Fallback to existing gene column
    if gene_col and gene_col in row.index:
        gene = str(row[gene_col]).strip()
        if gene and gene != 'nan':
            return gene

    # Fallback to locus_tag itself
    if locus_tag:
        return locus_tag

    return ''


def main():
    parser = argparse.ArgumentParser(
        description='Create mutation summary with sample column and gene name lookup'
    )
    parser.add_argument('--input-dir', '-i', required=True,
                        help='Directory containing per-sample *_mutations.csv files')
    parser.add_argument('--gene-lookup', '-g', required=True,
                        help='Path to gene_lookup_table.csv')
    parser.add_argument('--output-dir', '-o', default=None,
                        help='Output directory (defaults to input-dir)')
    parser.add_argument('--exclude-strains', nargs='+', default=['NPSA'],
                        help='Strain prefixes to exclude in filtered output (default: NPSA)')

    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir) if args.output_dir else input_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    exclude_prefixes = [p.upper() for p in args.exclude_strains]

    print("=" * 70)
    print("CREATE MUTATION SUMMARY WITH SAMPLE IDENTIFICATION")
    print("=" * 70)

    # 1. Load gene lookup table
    lookup_path = Path(args.gene_lookup)
    print(f"\n1. Loading gene lookup table: {lookup_path}")
    gene_lookup = load_gene_lookup(lookup_path)
    print(f"   Loaded {len(gene_lookup)} locus_tag entries")

    # 2. Find per-sample mutation CSV files
    print(f"\n2. Scanning for mutation CSV files in: {input_dir}")
    mutation_files = sorted(input_dir.glob("*_mutations.csv"))

    # Exclude any existing summary files
    mutation_files = [f for f in mutation_files
                      if 'summary' not in f.name.lower()
                      and 'filtered' not in f.name.lower()
                      and 'all_mutations' not in f.name.lower()]

    if not mutation_files:
        print("   No per-sample *_mutations.csv files found!")
        sys.exit(1)

    print(f"   Found {len(mutation_files)} sample files:")
    for f in mutation_files:
        print(f"      {f.name}")

    # 3. Read and combine all files
    print(f"\n3. Reading and combining mutation files...")
    all_dfs = []

    for mut_file in mutation_files:
        sample_name = extract_sample_name(mut_file.name)
        df = pd.read_csv(mut_file)

        if df.empty:
            print(f"   {sample_name}: 0 mutations (skipped)")
            continue

        # Ensure sample column has the correct sample name from filename
        df['sample'] = sample_name

        all_dfs.append(df)
        print(f"   {sample_name}: {len(df)} mutations")

    if not all_dfs:
        print("\nNo mutations found in any file!")
        sys.exit(1)

    combined = pd.concat(all_dfs, ignore_index=True)
    print(f"\n   Total mutations across all samples: {len(combined)}")
    print(f"   Unique samples: {combined['sample'].nunique()}")

    # 4. Resolve gene names using lookup table
    print(f"\n4. Resolving gene names from lookup table...")

    # Detect column names
    locus_col = next((c for c in ['locus_tag', 'LOCUS_TAG', 'gene_id'] if c in combined.columns), None)
    gene_col = next((c for c in ['gene', 'gene_name', 'GENE'] if c in combined.columns), None)

    print(f"   Locus tag column: {locus_col}")
    print(f"   Gene column: {gene_col}")

    combined['gene_name_resolved'] = combined.apply(
        lambda row: resolve_gene_name(row, locus_col, gene_col, gene_lookup), axis=1
    )

    n_resolved = (combined['gene_name_resolved'] != '').sum()
    n_from_lookup = 0
    if locus_col:
        for _, row in combined.iterrows():
            lt = str(row.get(locus_col, '')).strip()
            if lt != 'nan' and lt in gene_lookup:
                n_from_lookup += 1

    print(f"   Resolved gene names: {n_resolved} / {len(combined)}")
    print(f"   Matched via lookup table: {n_from_lookup}")

    # 5. Reorder columns so sample is first
    cols = list(combined.columns)
    # Move 'sample' to first position
    if 'sample' in cols:
        cols.remove('sample')
    cols = ['sample'] + cols
    # Put gene_name_resolved early
    if 'gene_name_resolved' in cols:
        cols.remove('gene_name_resolved')
        # Insert after gene-related columns
        insert_idx = 1
        for target in ['gene', 'gene_name', 'GENE', 'locus_tag', 'LOCUS_TAG']:
            if target in cols:
                insert_idx = cols.index(target) + 1
                break
        cols.insert(insert_idx, 'gene_name_resolved')

    combined = combined[cols]

    # 6. Write full summary
    full_output = output_dir / 'all_mutations_summary_with_samples.csv'
    combined.to_csv(full_output, index=False)
    print(f"\n5. Saved full summary: {full_output}")
    print(f"   ({len(combined)} mutations, {combined['sample'].nunique()} samples)")

    # 7. Create filtered summary (exclude NPSA strains)
    print(f"\n6. Creating filtered summary (excluding {', '.join(exclude_prefixes)} strains)...")

    def should_exclude(sample_name: str) -> bool:
        name_upper = sample_name.upper()
        return any(name_upper.startswith(prefix) for prefix in exclude_prefixes)

    mask_exclude = combined['sample'].apply(should_exclude)
    n_excluded = mask_exclude.sum()
    excluded_samples = combined.loc[mask_exclude, 'sample'].unique()
    filtered = combined[~mask_exclude].copy()

    print(f"   Excluded {n_excluded} mutations from samples: {', '.join(excluded_samples) if len(excluded_samples) > 0 else 'none'}")
    print(f"   Remaining: {len(filtered)} mutations from {filtered['sample'].nunique()} samples")

    filtered_output = output_dir / 'mutations_summary_filtered_no_NPSA.csv'
    filtered.to_csv(filtered_output, index=False)
    print(f"   Saved filtered summary: {filtered_output}")

    # 8. Print quick overview
    print(f"\n7. Quick overview:")
    print("-" * 60)
    print(f"{'Sample':<25} {'Mutations':>10}")
    print("-" * 60)
    for sample in sorted(combined['sample'].unique()):
        count = len(combined[combined['sample'] == sample])
        marker = " [EXCLUDED]" if should_exclude(sample) else ""
        print(f"   {sample:<25} {count:>6}{marker}")
    print("-" * 60)

    print(f"\n" + "=" * 70)
    print("DONE!")
    print("=" * 70)
    print(f"\nOutput files:")
    print(f"  - {full_output.name:<50} : All mutations with sample column")
    print(f"  - {filtered_output.name:<50} : Filtered (no {'/'.join(exclude_prefixes)} strains)")


if __name__ == '__main__':
    main()
