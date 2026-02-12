#!/usr/bin/env python3
"""
Parse NCBI protein summary results and create gene name mapping.

This script:
1. Parses NCBI protein summary format
2. Extracts RefSeq IDs and protein descriptions
3. Extracts gene names from descriptions
4. Merges with gene-refseq mapping to create final annotation table

Usage:
    python parse_ncbi_results.py --ncbi <protein_result.txt> --mapping <gene_refseq_mapping.csv> --output <output_dir>
"""

import argparse
import re
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import pandas as pd


def extract_gene_name_from_description(description: str) -> Tuple[str, str]:
    """
    Extract short gene name from NCBI protein description.

    Examples:
    - "pseudouridine-5'-phosphate glycosidase" -> ("psuG", "pseudouridine-5'-phosphate glycosidase")
    - "urease accessory protein UreG" -> ("ureG", "urease accessory protein")
    - "phenylalanine--tRNA ligase subunit alpha" -> ("pheS", "phenylalanine--tRNA ligase subunit alpha")
    - "DUF2273 domain-containing protein" -> ("DUF2273", "DUF2273 domain-containing protein")

    Returns: (short_gene_name, cleaned_description)
    """
    description = description.strip()

    # Remove "MULTISPECIES: " prefix if present
    description = re.sub(r'^MULTISPECIES:\s*', '', description)

    # Remove organism suffix like "[Staphylococcus]"
    description = re.sub(r'\s*\[.*?\]\s*$', '', description)

    short_name = ''

    # Pattern 1: Gene name at end (e.g., "urease accessory protein UreG")
    match = re.search(r'\b([A-Z][a-z]{2,3}[A-Z0-9]?)\s*$', description)
    if match:
        name = match.group(1)
        # Convert to standard format: UreG -> ureG
        short_name = name[0].lower() + name[1:]
        description = description[:match.start()].strip()

    # Pattern 2: Gene name with number at end (e.g., "ribosomal protein S17")
    if not short_name:
        match = re.search(r'\b([A-Z]+)(\d+)\s*$', description)
        if match:
            letters = match.group(1).lower()
            numbers = match.group(2)
            if len(letters) <= 3:
                short_name = letters + numbers

    # Pattern 3: Enzyme with standard name (e.g., "phenylalanine--tRNA ligase")
    if not short_name:
        # Common patterns for standard gene names
        enzyme_patterns = {
            r'phenylalanine.*tRNA ligase.*alpha': 'pheS',
            r'phenylalanine.*tRNA ligase.*beta': 'pheT',
            r'urease accessory protein': 'ure',
            r'transcriptional.*regulator\s+(\w+)': r'\1',
            r'repressor\s+(\w+)': r'\1',
            r'(\w+)\s+synthase': r'\1S',
            r'(\w+)\s+kinase': r'\1K',
            r'(\w+)\s+dehydrogenase': r'\1DH',
        }
        for pattern, replacement in enzyme_patterns.items():
            match = re.search(pattern, description, re.IGNORECASE)
            if match:
                if '\\1' in replacement:
                    short_name = re.sub(pattern, replacement, description, flags=re.IGNORECASE)
                    short_name = short_name.split()[0] if short_name else ''
                else:
                    short_name = replacement
                break

    # Pattern 4: DUF or COG family proteins
    if not short_name:
        match = re.search(r'(DUF\d+|COG\d+)', description, re.IGNORECASE)
        if match:
            short_name = match.group(1)

    # Pattern 5: HIT family, ABC transporter, etc.
    if not short_name:
        match = re.search(r'^(\w+)\s+family', description, re.IGNORECASE)
        if match:
            short_name = match.group(1)

    # Pattern 6: 30S/50S ribosomal protein
    if not short_name:
        match = re.search(r'(\d+)S ribosomal protein\s+([SL]\d+)', description, re.IGNORECASE)
        if match:
            short_name = 'rp' + match.group(2)

    # Fallback: use first word if it looks like a gene name
    if not short_name:
        first_word = description.split()[0] if description else ''
        if re.match(r'^[a-zA-Z]{2,6}$', first_word):
            short_name = first_word.lower()

    # If still no name, mark as hypothetical or use description
    if not short_name:
        if 'hypothetical' in description.lower():
            short_name = 'hyp'
        else:
            # Use first ~10 chars of description
            short_name = re.sub(r'[^a-zA-Z0-9]', '', description)[:10]

    return short_name, description


def parse_ncbi_summary(ncbi_file: Path) -> List[dict]:
    """
    Parse NCBI protein summary format.

    Format:
    1. MULTISPECIES: description [Organism]
    307 aa protein
    WP_000002068.1 GI:445924213

    (blank line)

    2. MULTISPECIES: description [Organism]
    ...
    """
    entries = []

    with open(ncbi_file, 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()

    # Split by entry number pattern
    entry_pattern = re.compile(r'(\d+)\.\s+(.+?)(?=\n\d+\.\s+|\Z)', re.DOTALL)

    for match in entry_pattern.finditer(content):
        entry_num = match.group(1)
        entry_text = match.group(2).strip()

        lines = entry_text.split('\n')
        if len(lines) < 2:
            continue

        # First line: description
        description = lines[0].strip()

        # Find RefSeq ID line
        refseq_id = ''
        protein_size = ''

        for line in lines[1:]:
            line = line.strip()

            # Check for RefSeq ID (WP_, NP_, YP_)
            refseq_match = re.search(r'([WNY]P_\d+\.\d+)', line)
            if refseq_match:
                refseq_id = refseq_match.group(1)

            # Check for protein size
            size_match = re.search(r'(\d+)\s*aa\s*protein', line)
            if size_match:
                protein_size = size_match.group(1)

        if refseq_id:
            short_name, clean_desc = extract_gene_name_from_description(description)

            entries.append({
                'refseq_id': refseq_id,
                'ncbi_description': description,
                'ncbi_gene_name': short_name,
                'clean_description': clean_desc,
                'protein_size_aa': protein_size
            })

    return entries


def main():
    parser = argparse.ArgumentParser(
        description='Parse NCBI protein results and create gene name mapping'
    )
    parser.add_argument('--ncbi', '-n', required=True,
                        help='NCBI protein summary results file')
    parser.add_argument('--mapping', '-m', required=False,
                        help='Gene-RefSeq mapping CSV (from extract_refseq_from_gff.py)')
    parser.add_argument('--mutations', required=False,
                        help='Mutations CSV file to update')
    parser.add_argument('--output', '-o', required=True,
                        help='Output directory')

    args = parser.parse_args()

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("="*70)
    print("PARSE NCBI PROTEIN RESULTS")
    print("="*70)

    # Parse NCBI results
    print(f"\n1. Parsing NCBI results: {args.ncbi}")
    ncbi_entries = parse_ncbi_summary(Path(args.ncbi))
    print(f"   Parsed {len(ncbi_entries)} protein entries")

    # Create NCBI results DataFrame
    ncbi_df = pd.DataFrame(ncbi_entries)

    # Save NCBI parsed results
    ncbi_parsed_file = output_dir / 'ncbi_parsed.csv'
    ncbi_df.to_csv(ncbi_parsed_file, index=False)
    print(f"   Saved parsed NCBI results: {ncbi_parsed_file}")

    # Show some examples
    print("\n   Sample gene names extracted:")
    for _, row in ncbi_df.head(15).iterrows():
        print(f"      {row['refseq_id']}: {row['ncbi_gene_name']} <- {row['clean_description'][:50]}...")

    # Merge with gene mapping if provided
    if args.mapping:
        print(f"\n2. Merging with gene mapping: {args.mapping}")
        mapping_df = pd.read_csv(args.mapping)

        # Merge on RefSeq ID
        merged_df = mapping_df.merge(
            ncbi_df[['refseq_id', 'ncbi_gene_name', 'ncbi_description', 'clean_description']],
            left_on='refseq_primary',
            right_on='refseq_id',
            how='left'
        )

        # Create best gene name column
        def get_best_name(row):
            # Priority: ncbi_gene_name > gene_name from GFF > locus_tag
            if pd.notna(row.get('ncbi_gene_name')) and str(row['ncbi_gene_name']).strip():
                return str(row['ncbi_gene_name']).strip()
            if pd.notna(row.get('gene_name')) and str(row['gene_name']).strip():
                name = str(row['gene_name']).strip()
                if len(name) <= 10:  # Reasonable gene name length
                    return name
            return row.get('locus_tag', 'unknown')

        merged_df['gene_best'] = merged_df.apply(get_best_name, axis=1)

        # Save merged mapping
        merged_file = output_dir / 'gene_mapping_complete.csv'
        merged_df.to_csv(merged_file, index=False)
        print(f"   Saved complete gene mapping: {merged_file}")

        # Statistics
        n_with_ncbi_name = merged_df['ncbi_gene_name'].notna().sum()
        print(f"   Genes with NCBI names: {n_with_ncbi_name} / {len(merged_df)}")

    # Update mutations if provided
    if args.mutations:
        print(f"\n3. Updating mutations file: {args.mutations}")
        mut_df = pd.read_csv(args.mutations)

        # Create RefSeq to gene name lookup
        refseq_to_gene = dict(zip(ncbi_df['refseq_id'], ncbi_df['ncbi_gene_name']))
        refseq_to_desc = dict(zip(ncbi_df['refseq_id'], ncbi_df['clean_description']))

        # Try to match mutations to NCBI names
        # First, check if there's a refseq column in mutations
        refseq_cols = [c for c in mut_df.columns if 'refseq' in c.lower()]

        if refseq_cols:
            for col in refseq_cols:
                mut_df['ncbi_gene'] = mut_df[col].map(refseq_to_gene)
                mut_df['ncbi_desc'] = mut_df[col].map(refseq_to_desc)

        # Save updated mutations
        mut_updated_file = output_dir / 'mutations_ncbi_annotated.csv'
        mut_df.to_csv(mut_updated_file, index=False)
        print(f"   Saved updated mutations: {mut_updated_file}")

    # Create a simple lookup table for manual use
    print("\n4. Creating simple lookup table...")
    lookup_df = ncbi_df[['refseq_id', 'ncbi_gene_name', 'clean_description']].copy()
    lookup_df.columns = ['RefSeq_ID', 'Gene_Name', 'Description']
    lookup_file = output_dir / 'refseq_gene_lookup.csv'
    lookup_df.to_csv(lookup_file, index=False)
    print(f"   Saved lookup table: {lookup_file}")

    print("\n" + "="*70)
    print("COMPLETE")
    print("="*70)
    print(f"\nOutput files:")
    print(f"  - ncbi_parsed.csv: Raw NCBI parsing results")
    print(f"  - refseq_gene_lookup.csv: Simple RefSeq -> Gene name lookup")
    if args.mapping:
        print(f"  - gene_mapping_complete.csv: Full mapping with locus_tag, RefSeq, gene name")
    if args.mutations:
        print(f"  - mutations_ncbi_annotated.csv: Mutations with NCBI gene names")

    print("\nYou can now use 'refseq_gene_lookup.csv' or 'gene_mapping_complete.csv'")
    print("to annotate your mutation data with proper gene names.")


if __name__ == '__main__':
    main()
