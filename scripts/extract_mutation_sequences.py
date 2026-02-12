#!/usr/bin/env python3
"""
Extract sequences around mutations for BLAST verification.

This script:
1. Reads mutation CSV files
2. Extracts flanking DNA sequences from the reference genome
3. Translates to protein sequences where applicable
4. Outputs sequences in FASTA format for BLAST verification

Usage:
    python extract_mutation_sequences.py --mutations <csv> --reference <fasta> --output <dir>
"""

import sys
import argparse
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import re

import pandas as pd


def parse_fasta(fasta_path: Path) -> Dict[str, str]:
    """Parse a FASTA file and return dict of sequences."""
    sequences = {}
    current_id = None
    current_seq = []

    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.upper())

        if current_id:
            sequences[current_id] = ''.join(current_seq)

    return sequences


# Genetic code for translation
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


def reverse_complement(seq: str) -> str:
    """Return reverse complement of DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))


def translate(dna_seq: str) -> str:
    """Translate DNA to protein."""
    protein = []
    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i+3]
        aa = CODON_TABLE.get(codon, 'X')
        if aa == '*':
            break
        protein.append(aa)
    return ''.join(protein)


def extract_sequence(
    genome: Dict[str, str],
    chrom: str,
    position: int,
    flank_size: int = 500
) -> Tuple[str, int, int]:
    """Extract sequence around a position with flanking regions."""

    # Find the right chromosome
    seq = None
    for name, s in genome.items():
        if chrom in name or name in chrom:
            seq = s
            break

    if seq is None:
        # Try first chromosome
        seq = list(genome.values())[0]

    start = max(0, position - flank_size - 1)
    end = min(len(seq), position + flank_size)

    return seq[start:end], start + 1, end


def main():
    parser = argparse.ArgumentParser(
        description='Extract sequences around mutations for BLAST verification'
    )
    parser.add_argument('--mutations', '-m', required=True,
                        help='Mutation CSV file (combined or individual)')
    parser.add_argument('--reference', '-r', required=True,
                        help='Reference genome FASTA file')
    parser.add_argument('--output', '-o', required=True,
                        help='Output directory')
    parser.add_argument('--flank', '-f', type=int, default=500,
                        help='Flanking region size in bp (default: 500)')
    parser.add_argument('--top', '-t', type=int, default=None,
                        help='Only process top N unique positions')
    parser.add_argument('--genes', '-g', nargs='+', default=None,
                        help='Only extract specific genes (by any name column)')

    args = parser.parse_args()

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("="*70)
    print("MUTATION SEQUENCE EXTRACTION FOR BLAST")
    print("="*70)

    # Load reference genome
    print(f"\n1. Loading reference genome: {args.reference}")
    genome = parse_fasta(Path(args.reference))
    total_len = sum(len(s) for s in genome.values())
    print(f"   Loaded {len(genome)} contigs, total length: {total_len:,} bp")

    # Load mutations
    print(f"\n2. Loading mutations: {args.mutations}")
    df = pd.read_csv(args.mutations)
    print(f"   Loaded {len(df)} mutations")

    # Detect column names
    pos_col = None
    for col in ['POS', 'position', 'pos', 'Position']:
        if col in df.columns:
            pos_col = col
            break

    chrom_col = None
    for col in ['CHROM', 'chromosome', 'chrom', 'Chromosome']:
        if col in df.columns:
            chrom_col = col
            break

    gene_cols = []
    for col in ['gene_name', 'GENE', 'gene_short', 'gene_name_new', 'product', 'PRODUCT']:
        if col in df.columns:
            gene_cols.append(col)

    ref_col = None
    for col in ['REF', 'reference', 'ref', 'Reference']:
        if col in df.columns:
            ref_col = col
            break

    alt_col = None
    for col in ['ALT', 'alternative', 'alt', 'Alternative']:
        if col in df.columns:
            alt_col = col
            break

    effect_col = None
    for col in ['EFFECT', 'effect', 'Effect']:
        if col in df.columns:
            effect_col = col
            break

    # Filter by genes if specified
    if args.genes:
        print(f"\n   Filtering for genes: {args.genes}")
        mask = pd.Series([False] * len(df))
        for col in gene_cols:
            for gene in args.genes:
                mask |= df[col].astype(str).str.lower().str.contains(gene.lower(), na=False)
        df = df[mask]
        print(f"   {len(df)} mutations match specified genes")

    # Get unique positions
    df_unique = df.drop_duplicates(subset=[pos_col])
    print(f"\n3. Processing {len(df_unique)} unique positions")

    if args.top:
        df_unique = df_unique.head(args.top)
        print(f"   (Limited to top {args.top})")

    # Extract sequences
    print("\n4. Extracting sequences...")

    dna_records = []
    protein_records = []
    summary_data = []

    for idx, row in df_unique.iterrows():
        pos = int(row[pos_col])
        chrom = row[chrom_col] if chrom_col and not pd.isna(row[chrom_col]) else 'contig_1'

        # Get gene names
        gene_names = []
        for col in gene_cols:
            val = row.get(col, '')
            if not pd.isna(val) and str(val).strip():
                gene_names.append(str(val).strip())
        gene_str = '_'.join(set(gene_names)) if gene_names else 'unknown'
        gene_str = re.sub(r'[^\w\-]', '_', gene_str)[:50]  # Clean for FASTA header

        ref = row.get(ref_col, 'N') if ref_col else 'N'
        alt = row.get(alt_col, 'N') if alt_col else 'N'
        effect = row.get(effect_col, 'unknown') if effect_col else 'unknown'

        # Extract DNA sequence
        seq, start, end = extract_sequence(genome, chrom, pos, args.flank)

        # Create FASTA header
        header = f"pos{pos}_{gene_str}|{chrom}:{start}-{end}|{ref}>{alt}|{effect}"

        dna_records.append(f">{header}\n{seq}")

        # Translate if it's a coding mutation
        if 'missense' in str(effect).lower() or 'frameshift' in str(effect).lower():
            # Try to translate the region
            protein = translate(seq)
            if len(protein) > 10:
                protein_records.append(f">{header}\n{protein}")

        # Summary
        summary_data.append({
            'position': pos,
            'chromosome': chrom,
            'gene_names': ', '.join(set(gene_names)),
            'ref': ref,
            'alt': alt,
            'effect': effect,
            'seq_start': start,
            'seq_end': end,
            'seq_length': len(seq)
        })

    # Write DNA FASTA
    dna_fasta = output_dir / 'mutation_sequences_dna.fasta'
    with open(dna_fasta, 'w') as f:
        f.write('\n'.join(dna_records))
    print(f"\n   DNA sequences saved: {dna_fasta}")

    # Write protein FASTA
    if protein_records:
        protein_fasta = output_dir / 'mutation_sequences_protein.fasta'
        with open(protein_fasta, 'w') as f:
            f.write('\n'.join(protein_records))
        print(f"   Protein sequences saved: {protein_fasta}")

    # Write summary
    summary_df = pd.DataFrame(summary_data)
    summary_file = output_dir / 'extraction_summary.csv'
    summary_df.to_csv(summary_file, index=False)
    print(f"   Summary saved: {summary_file}")

    # Instructions
    print("\n" + "="*70)
    print("BLAST INSTRUCTIONS")
    print("="*70)
    print("""
To verify gene annotations:

1. ONLINE BLAST (recommended for small numbers):
   - Go to https://blast.ncbi.nlm.nih.gov/Blast.cgi
   - Use 'Nucleotide BLAST' for DNA or 'Protein BLAST' for protein
   - Paste sequences from the FASTA files
   - Select database: 'nr' (non-redundant) or 'refseq_protein'
   - Organism: Staphylococcus aureus (taxid:1280)

2. COMMAND LINE BLAST (for many sequences):
   # Install BLAST+
   sudo apt-get install ncbi-blast+

   # Download S. aureus protein database
   update_blastdb.pl --decompress swissprot

   # Run BLAST
   blastp -query mutation_sequences_protein.fasta \\
          -db swissprot \\
          -outfmt "6 qseqid sseqid pident length evalue stitle" \\
          -max_target_seqs 5 \\
          -evalue 1e-10 \\
          -out blast_results.tsv

3. UNIPROT SEARCH:
   - Go to https://www.uniprot.org/blast
   - Paste protein sequence
   - Filter by "Staphylococcus aureus"
   - Get standardized gene names

Output files:
   - mutation_sequences_dna.fasta    : DNA sequences for nucleotide BLAST
   - mutation_sequences_protein.fasta: Protein sequences for protein BLAST
   - extraction_summary.csv          : Summary of extracted regions
""")


if __name__ == '__main__':
    main()
