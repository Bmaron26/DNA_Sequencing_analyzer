#!/usr/bin/env python3
"""
Update mutation annotations using BLAST results.

This script:
1. Parses BLAST results (outfmt 6)
2. Extracts gene names from BLAST hits
3. Updates mutation CSV files with verified gene names
4. Creates a gene name mapping table

Usage:
    python update_annotations_from_blast.py --blast <blast_results.tsv> --mutations <mutations.csv> --output <output_dir>
"""

import argparse
import re
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from collections import defaultdict

import pandas as pd


def parse_blast_results(blast_file: Path) -> Dict[str, List[dict]]:
    """
    Parse BLAST results in outfmt 6 format.

    Expected columns: qseqid sseqid pident length evalue stitle

    Returns dict mapping query_id -> list of hits
    """
    hits = defaultdict(list)

    with open(blast_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 6:
                continue

            qseqid = parts[0]
            sseqid = parts[1]
            pident = float(parts[2])
            length = int(parts[3])
            evalue = float(parts[4])
            stitle = parts[5] if len(parts) > 5 else ''

            # Extract position from query ID (format: pos123456_genename|...)
            pos_match = re.match(r'pos(\d+)', qseqid)
            position = int(pos_match.group(1)) if pos_match else None

            hits[qseqid].append({
                'position': position,
                'subject_id': sseqid,
                'percent_identity': pident,
                'alignment_length': length,
                'evalue': evalue,
                'subject_title': stitle,
                'query_id': qseqid
            })

    return dict(hits)


def extract_gene_name_from_title(stitle: str) -> Tuple[str, str, str]:
    """
    Extract gene name from BLAST subject title.

    UniProt format examples:
    - "sp|P0A0E0|AGRC_STAAN Accessory gene regulator protein C OS=..."
    - "tr|A0A0H3JQR1|A0A0H3JQR1_STAAU GraS protein OS=..."
    - "DNA-directed RNA polymerase subunit beta' OS=Staphylococcus aureus"

    Returns: (short_gene_name, uniprot_id, full_description)
    """
    stitle = str(stitle)

    # Try to extract UniProt ID and gene name
    # Pattern: sp|XXXXX|GENENAME_SPECIES description
    uniprot_match = re.match(r'(sp|tr)\|([^|]+)\|(\w+)_\w+\s+(.+?)(?:\s+OS=|$)', stitle)
    if uniprot_match:
        db_type = uniprot_match.group(1)
        uniprot_id = uniprot_match.group(2)
        gene_code = uniprot_match.group(3)
        description = uniprot_match.group(4)

        # Gene code might be like "AGRC" -> "agrC" or "RPOC" -> "rpoC"
        if gene_code.isupper() and len(gene_code) <= 6:
            # Convert to standard format (lowercase with capital at end if present)
            short_name = gene_code[:-1].lower() + gene_code[-1]
        else:
            short_name = gene_code

        return short_name, uniprot_id, description

    # Try pattern: just description with gene name at start or end
    # e.g., "DNA-directed RNA polymerase subunit beta'"
    # Look for common gene name patterns in the title
    gene_patterns = [
        r'\b([a-z]{3,4}[A-Z])\b',  # rpoC, graS
        r'\b([A-Z][a-z]{2,3}[A-Z])\b',  # RpoC, GraS
        r'^(\w{3,6})\s+protein',  # "GraS protein"
        r'protein\s+(\w{3,6})(?:\s|$)',  # "protein GraS"
    ]

    for pattern in gene_patterns:
        match = re.search(pattern, stitle)
        if match:
            gene = match.group(1)
            # Standardize: lowercase first letters, keep last capital
            if len(gene) <= 5:
                standardized = gene[:-1].lower() + gene[-1] if gene[-1].isupper() else gene.lower()
                return standardized, '', stitle.split(' OS=')[0] if ' OS=' in stitle else stitle

    # Fall back to first word or UniProt ID
    first_word = stitle.split()[0] if stitle else ''
    return first_word, '', stitle.split(' OS=')[0] if ' OS=' in stitle else stitle


def get_best_gene_name(hits: List[dict], min_identity: float = 50.0, min_length: int = 50) -> dict:
    """
    Get the best gene name from BLAST hits.

    Prioritizes:
    1. High identity matches (>90%)
    2. SwissProt (sp|) over TrEMBL (tr|)
    3. Longer alignments
    """
    if not hits:
        return None

    # Filter by quality
    good_hits = [h for h in hits if h['percent_identity'] >= min_identity
                 and h['alignment_length'] >= min_length]

    if not good_hits:
        good_hits = hits  # Use all if none pass filter

    # Sort by identity, then by whether it's SwissProt
    def sort_key(h):
        is_swissprot = 'sp|' in h['subject_id'] or h['subject_title'].startswith('sp|')
        return (-h['percent_identity'], -is_swissprot, -h['alignment_length'])

    good_hits.sort(key=sort_key)

    best = good_hits[0]
    gene_name, uniprot_id, description = extract_gene_name_from_title(best['subject_title'])

    return {
        'blast_gene_name': gene_name,
        'blast_uniprot_id': uniprot_id,
        'blast_description': description,
        'blast_identity': best['percent_identity'],
        'blast_evalue': best['evalue'],
        'blast_subject_id': best['subject_id'],
        'position': best['position']
    }


def main():
    parser = argparse.ArgumentParser(
        description='Update mutation annotations using BLAST results'
    )
    parser.add_argument('--blast', '-b', required=True,
                        help='BLAST results file (outfmt 6)')
    parser.add_argument('--mutations', '-m', required=True,
                        help='Mutations CSV file to update')
    parser.add_argument('--output', '-o', required=True,
                        help='Output directory')
    parser.add_argument('--min-identity', type=float, default=50.0,
                        help='Minimum percent identity (default: 50)')
    parser.add_argument('--min-length', type=int, default=50,
                        help='Minimum alignment length (default: 50)')

    args = parser.parse_args()

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("="*70)
    print("UPDATE ANNOTATIONS FROM BLAST RESULTS")
    print("="*70)

    # Parse BLAST results
    print(f"\n1. Parsing BLAST results: {args.blast}")
    blast_hits = parse_blast_results(Path(args.blast))
    print(f"   Found hits for {len(blast_hits)} queries")

    # Get best gene name for each position
    print("\n2. Extracting best gene names...")
    position_to_gene = {}
    gene_mapping = []

    for query_id, hits in blast_hits.items():
        best = get_best_gene_name(hits, args.min_identity, args.min_length)
        if best and best['position']:
            position_to_gene[best['position']] = best
            gene_mapping.append({
                'query_id': query_id,
                **best
            })

    print(f"   Mapped {len(position_to_gene)} positions to gene names")

    # Save gene mapping
    mapping_df = pd.DataFrame(gene_mapping)
    mapping_file = output_dir / 'blast_gene_mapping.csv'
    mapping_df.to_csv(mapping_file, index=False)
    print(f"   Saved mapping to: {mapping_file}")

    # Print summary of found genes
    print("\n   Top BLAST-verified genes:")
    if not mapping_df.empty:
        gene_counts = mapping_df['blast_gene_name'].value_counts().head(20)
        for gene, count in gene_counts.items():
            print(f"      {gene}: {count} mutations")

    # Load and update mutations file
    print(f"\n3. Updating mutations file: {args.mutations}")
    mutations_df = pd.read_csv(args.mutations)
    print(f"   Loaded {len(mutations_df)} mutations")

    # Detect position column
    pos_col = None
    for col in ['POS', 'position', 'pos', 'Position']:
        if col in mutations_df.columns:
            pos_col = col
            break

    if pos_col is None:
        print("   ERROR: No position column found!")
        return

    # Add BLAST annotation columns
    blast_gene_names = []
    blast_uniprot_ids = []
    blast_descriptions = []
    blast_identities = []
    blast_verified = []

    for _, row in mutations_df.iterrows():
        pos = row[pos_col]
        if pd.notna(pos) and int(pos) in position_to_gene:
            info = position_to_gene[int(pos)]
            blast_gene_names.append(info['blast_gene_name'])
            blast_uniprot_ids.append(info['blast_uniprot_id'])
            blast_descriptions.append(info['blast_description'])
            blast_identities.append(info['blast_identity'])
            blast_verified.append(True)
        else:
            blast_gene_names.append('')
            blast_uniprot_ids.append('')
            blast_descriptions.append('')
            blast_identities.append('')
            blast_verified.append(False)

    mutations_df['blast_gene'] = blast_gene_names
    mutations_df['blast_uniprot'] = blast_uniprot_ids
    mutations_df['blast_description'] = blast_descriptions
    mutations_df['blast_identity'] = blast_identities
    mutations_df['blast_verified'] = blast_verified

    # Create a "best" gene name column that prefers BLAST results
    def get_best_name(row):
        if row.get('blast_gene') and str(row['blast_gene']).strip():
            return row['blast_gene']
        for col in ['gene_short', 'gene_name', 'GENE', 'gene_name_new']:
            if col in row and pd.notna(row[col]) and str(row[col]).strip():
                return str(row[col]).strip()
        return 'unknown'

    mutations_df['gene_best'] = mutations_df.apply(get_best_name, axis=1)

    # Save updated mutations
    updated_file = output_dir / 'mutations_blast_annotated.csv'
    mutations_df.to_csv(updated_file, index=False)
    print(f"   Saved updated mutations to: {updated_file}")

    # Create summary by gene
    print("\n4. Creating gene summary...")
    verified = mutations_df[mutations_df['blast_verified'] == True]
    print(f"   BLAST-verified mutations: {len(verified)} / {len(mutations_df)}")

    gene_summary = mutations_df.groupby('gene_best').agg({
        pos_col: 'count',
        'blast_verified': 'sum',
        'sample': 'nunique' if 'sample' in mutations_df.columns else 'count'
    }).reset_index()
    gene_summary.columns = ['gene', 'total_mutations', 'blast_verified', 'n_samples']
    gene_summary = gene_summary.sort_values('total_mutations', ascending=False)

    summary_file = output_dir / 'gene_summary_blast.csv'
    gene_summary.to_csv(summary_file, index=False)
    print(f"   Saved gene summary to: {summary_file}")

    print("\n" + "="*70)
    print("ANNOTATION UPDATE COMPLETE")
    print("="*70)
    print(f"\nOutput files:")
    print(f"  - {mapping_file}: BLAST hit to gene name mapping")
    print(f"  - {updated_file}: Mutations with BLAST annotations")
    print(f"  - {summary_file}: Gene summary with counts")
    print(f"\nNew columns added to mutations file:")
    print(f"  - blast_gene: Gene name from BLAST")
    print(f"  - blast_uniprot: UniProt ID")
    print(f"  - blast_description: Full description")
    print(f"  - blast_identity: Percent identity")
    print(f"  - blast_verified: Whether BLAST annotation exists")
    print(f"  - gene_best: Best available gene name (BLAST preferred)")


if __name__ == '__main__':
    main()
