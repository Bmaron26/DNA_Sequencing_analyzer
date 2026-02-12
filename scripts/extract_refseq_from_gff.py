#!/usr/bin/env python3
"""
Extract RefSeq IDs and gene information from GFF3 annotation file.

This script:
1. Parses GFF3 file to extract locus_tag, gene name, product, RefSeq ID
2. Creates a mapping table for all genes
3. Optionally queries NCBI to get official gene names

Usage:
    python extract_refseq_from_gff.py --gff <gff3_file> --output <output_dir>
"""

import argparse
import re
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from collections import defaultdict

import pandas as pd


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
            value = value.replace('%2C', ',').replace('%3B', ';').replace('%25', '%').replace('%3D', '=')
            attributes[key] = value

    return attributes


def extract_refseq_from_note(note: str) -> List[str]:
    """Extract RefSeq IDs from Note field."""
    refseq_ids = []

    if not note:
        return refseq_ids

    # Pattern for RefSeq IDs: WP_XXXXXXXXX.X, NP_XXXXXXXXX.X, YP_XXXXXXXXX.X
    patterns = [
        r'RefSeq:([WNY]P_\d+\.\d+)',  # RefSeq:WP_000691584.1
        r'\b([WNY]P_\d+\.\d+)\b',      # Just the ID
    ]

    for pattern in patterns:
        matches = re.findall(pattern, note)
        refseq_ids.extend(matches)

    return list(set(refseq_ids))  # Remove duplicates


def extract_uniparc_from_note(note: str) -> str:
    """Extract UniParc ID from Note field."""
    if not note:
        return ''

    match = re.search(r'UniParc:(UPI[A-Z0-9]+)', note)
    return match.group(1) if match else ''


def extract_uniref_from_note(note: str) -> List[str]:
    """Extract UniRef IDs from Note field."""
    if not note:
        return []

    matches = re.findall(r'UniRef\d+_([A-Z0-9_]+)', note)
    return matches


def parse_gff3_file(gff_path: Path) -> List[dict]:
    """
    Parse GFF3 file and extract gene/CDS information with RefSeq IDs.
    """
    genes = []

    print(f"   Parsing GFF3 file: {gff_path}")

    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            seqid, source, feature_type, start, end, score, strand, phase, attributes = parts

            # Focus on CDS and gene features
            if feature_type not in ['CDS', 'gene']:
                continue

            attrs = parse_gff_attributes(attributes)

            # Extract basic info
            locus_tag = attrs.get('locus_tag', '')
            gene_name = attrs.get('gene', attrs.get('Name', ''))
            product = attrs.get('product', '')
            note = attrs.get('Note', attrs.get('note', ''))
            feature_id = attrs.get('ID', '')

            # Extract RefSeq and other IDs from Note
            refseq_ids = extract_refseq_from_note(note)
            uniparc = extract_uniparc_from_note(note)
            uniref_ids = extract_uniref_from_note(note)

            if locus_tag or gene_name:  # Only include if we have some identifier
                genes.append({
                    'locus_tag': locus_tag,
                    'gene_name': gene_name,
                    'product': product,
                    'feature_type': feature_type,
                    'feature_id': feature_id,
                    'start': int(start),
                    'end': int(end),
                    'strand': strand,
                    'refseq_ids': ','.join(refseq_ids) if refseq_ids else '',
                    'refseq_primary': refseq_ids[0] if refseq_ids else '',
                    'uniparc': uniparc,
                    'uniref_ids': ','.join(uniref_ids) if uniref_ids else '',
                    'note': note
                })

    return genes


def create_gene_refseq_table(genes: List[dict]) -> pd.DataFrame:
    """Create a clean table mapping genes to RefSeq IDs."""

    # Group by locus_tag, preferring CDS over gene features
    gene_map = {}

    for g in genes:
        locus = g['locus_tag']
        if not locus:
            continue

        if locus not in gene_map:
            gene_map[locus] = g
        elif g['feature_type'] == 'CDS' and gene_map[locus]['feature_type'] != 'CDS':
            # Prefer CDS features
            gene_map[locus] = g
        elif g['refseq_primary'] and not gene_map[locus]['refseq_primary']:
            # Prefer entries with RefSeq
            gene_map[locus] = g

    # Convert to DataFrame
    df = pd.DataFrame(list(gene_map.values()))

    # Reorder columns
    col_order = ['locus_tag', 'gene_name', 'refseq_primary', 'product',
                 'start', 'end', 'strand', 'refseq_ids', 'uniparc', 'uniref_ids']
    available_cols = [c for c in col_order if c in df.columns]
    df = df[available_cols]

    return df


def create_refseq_batch_file(df: pd.DataFrame, output_path: Path):
    """Create a file with RefSeq IDs for batch NCBI query."""

    refseq_ids = df[df['refseq_primary'] != '']['refseq_primary'].unique()

    with open(output_path, 'w') as f:
        for refseq in refseq_ids:
            f.write(f"{refseq}\n")

    return len(refseq_ids)


def query_ncbi_batch(refseq_ids: List[str], output_path: Path, batch_size: int = 100):
    """
    Query NCBI for gene information using RefSeq IDs.
    Requires biopython: pip install biopython
    """
    try:
        from Bio import Entrez
        Entrez.email = "your.email@example.com"  # Required by NCBI
    except ImportError:
        print("   Biopython not installed. Install with: pip install biopython")
        print("   Skipping NCBI query. Use the RefSeq IDs manually on NCBI website.")
        return None

    results = []
    total = len(refseq_ids)

    print(f"   Querying NCBI for {total} RefSeq IDs...")

    for i in range(0, total, batch_size):
        batch = refseq_ids[i:i+batch_size]
        query = ' OR '.join([f"{rid}[Accession]" for rid in batch])

        try:
            # Search for the proteins
            handle = Entrez.esearch(db="protein", term=query, retmax=batch_size)
            record = Entrez.read(handle)
            handle.close()

            if record['IdList']:
                # Fetch details
                handle = Entrez.efetch(db="protein", id=record['IdList'],
                                       rettype="gb", retmode="text")
                # Parse would go here...
                handle.close()

            time.sleep(0.5)  # Be nice to NCBI

        except Exception as e:
            print(f"   Warning: NCBI query failed: {e}")
            continue

        print(f"   Processed {min(i+batch_size, total)}/{total}")

    return results


def main():
    parser = argparse.ArgumentParser(
        description='Extract RefSeq IDs from GFF3 annotation file'
    )
    parser.add_argument('--gff', '-g', required=True,
                        help='GFF3 annotation file')
    parser.add_argument('--output', '-o', required=True,
                        help='Output directory')
    parser.add_argument('--query-ncbi', action='store_true',
                        help='Query NCBI for gene names (requires biopython)')

    args = parser.parse_args()

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("="*70)
    print("EXTRACT REFSEQ IDS FROM GFF3")
    print("="*70)

    # Parse GFF3
    print("\n1. Parsing GFF3 file...")
    genes = parse_gff3_file(Path(args.gff))
    print(f"   Found {len(genes)} gene/CDS features")

    # Create gene table
    print("\n2. Creating gene-RefSeq mapping table...")
    df = create_gene_refseq_table(genes)

    # Statistics
    n_with_refseq = (df['refseq_primary'] != '').sum()
    n_with_gene_name = (df['gene_name'] != '').sum()

    print(f"   Total unique genes: {len(df)}")
    print(f"   Genes with RefSeq ID: {n_with_refseq} ({100*n_with_refseq/len(df):.1f}%)")
    print(f"   Genes with gene name: {n_with_gene_name} ({100*n_with_gene_name/len(df):.1f}%)")

    # Save full table
    full_table = output_dir / 'gene_refseq_mapping.csv'
    df.to_csv(full_table, index=False)
    print(f"\n   Saved full table: {full_table}")

    # Save table of genes WITH RefSeq IDs
    df_with_refseq = df[df['refseq_primary'] != '']
    refseq_table = output_dir / 'genes_with_refseq.csv'
    df_with_refseq.to_csv(refseq_table, index=False)
    print(f"   Saved genes with RefSeq: {refseq_table}")

    # Save table of genes WITHOUT good names (need lookup)
    df_need_names = df[(df['gene_name'] == '') | (df['gene_name'].str.len() > 10)]
    need_names_table = output_dir / 'genes_needing_names.csv'
    df_need_names.to_csv(need_names_table, index=False)
    print(f"   Saved genes needing names: {need_names_table} ({len(df_need_names)} genes)")

    # Create batch file for NCBI
    print("\n3. Creating RefSeq batch file for NCBI lookup...")
    batch_file = output_dir / 'refseq_ids_for_ncbi.txt'
    n_refseq = create_refseq_batch_file(df, batch_file)
    print(f"   Saved {n_refseq} RefSeq IDs to: {batch_file}")

    # Print genes that need names
    print("\n4. Genes without proper names (first 20):")
    print("-"*70)
    for _, row in df_need_names.head(20).iterrows():
        print(f"   {row['locus_tag']}: RefSeq={row['refseq_primary'] or 'N/A'}, "
              f"product={row['product'][:50]}...")

    # Instructions
    print("\n" + "="*70)
    print("HOW TO GET GENE NAMES FROM NCBI")
    print("="*70)
    print("""
Option 1: NCBI Batch Entrez (Recommended for many genes)
---------------------------------------------------------
1. Go to: https://www.ncbi.nlm.nih.gov/sites/batchentrez
2. Select Database: "Protein"
3. Upload file: refseq_ids_for_ncbi.txt
4. Click "Retrieve"
5. Download results as "Summary" or "Tabular"

Option 2: NCBI Protein Search (For few genes)
---------------------------------------------
1. Go to: https://www.ncbi.nlm.nih.gov/protein/
2. Search for RefSeq ID (e.g., WP_000691584.1)
3. Look for "gene" field in the record

Option 3: UniProt ID Mapping (Best for standardized names)
----------------------------------------------------------
1. Go to: https://www.uniprot.org/id-mapping
2. Select "From": RefSeq Protein
3. Select "To": UniProtKB
4. Paste RefSeq IDs
5. Submit and download mapping

The downloaded table will have official gene names like 'rpoC', 'graS', etc.
""")

    # Query NCBI if requested
    if args.query_ncbi:
        print("\n5. Querying NCBI (this may take a while)...")
        refseq_list = df[df['refseq_primary'] != '']['refseq_primary'].tolist()
        ncbi_results = query_ncbi_batch(refseq_list, output_dir / 'ncbi_results.csv')


if __name__ == '__main__':
    main()
