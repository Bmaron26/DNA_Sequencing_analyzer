#!/usr/bin/env python3
"""
Parse breseq results and optionally compare with FreeBayes output.

Breseq outputs results in GenomeDiff (.gd) format, which contains:
- SNP: Single nucleotide polymorphisms
- SUB: Substitutions (multiple bases)
- DEL: Deletions
- INS: Insertions
- MOB: Mobile element insertions
- AMP: Amplifications (gene duplications)
- CON: Gene conversions
- INV: Inversions

This script parses the .gd files and creates a CSV summary compatible
with the FreeBayes mutation tables.

Usage:
    python parse_breseq_results.py \
        --breseq-dir /path/to/breseq/results \
        --output /path/to/output \
        --freebayes /path/to/freebayes/summary.csv  # optional comparison
"""

import argparse
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass

import pandas as pd


@dataclass
class BreseqMutation:
    """Represents a mutation from breseq output."""
    sample: str
    mutation_type: str  # SNP, DEL, INS, MOB, AMP, SUB, CON, INV
    seq_id: str
    position: int
    size: int  # For deletions/insertions
    ref: str
    alt: str
    gene_name: str
    gene_product: str
    annotation: str
    frequency: float
    evidence_id: str


def parse_gd_file(gd_path: Path, sample_name: str) -> List[BreseqMutation]:
    """
    Parse a breseq GenomeDiff (.gd) file.

    GenomeDiff format:
    - Lines starting with # are comments
    - Each mutation line has tab-separated fields
    - First field is mutation type (SNP, DEL, INS, MOB, etc.)
    """
    mutations = []

    with open(gd_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split('\t')
            if len(parts) < 3:
                continue

            mut_type = parts[0]

            # Parse based on mutation type
            if mut_type == 'SNP':
                # SNP evidence_id seq_id position new_seq [annotations...]
                if len(parts) >= 5:
                    mutation = parse_snp(parts, sample_name)
                    if mutation:
                        mutations.append(mutation)

            elif mut_type == 'DEL':
                # DEL evidence_id seq_id position size [annotations...]
                if len(parts) >= 5:
                    mutation = parse_del(parts, sample_name)
                    if mutation:
                        mutations.append(mutation)

            elif mut_type == 'INS':
                # INS evidence_id seq_id position new_seq [annotations...]
                if len(parts) >= 5:
                    mutation = parse_ins(parts, sample_name)
                    if mutation:
                        mutations.append(mutation)

            elif mut_type == 'SUB':
                # SUB evidence_id seq_id position size new_seq [annotations...]
                if len(parts) >= 6:
                    mutation = parse_sub(parts, sample_name)
                    if mutation:
                        mutations.append(mutation)

            elif mut_type == 'MOB':
                # MOB evidence_id seq_id position repeat_name strand duplication_size [annotations...]
                if len(parts) >= 6:
                    mutation = parse_mob(parts, sample_name)
                    if mutation:
                        mutations.append(mutation)

            elif mut_type == 'AMP':
                # AMP evidence_id seq_id position size new_copy_number [annotations...]
                if len(parts) >= 6:
                    mutation = parse_amp(parts, sample_name)
                    if mutation:
                        mutations.append(mutation)

    return mutations


def parse_annotations(parts: List[str], start_idx: int) -> Tuple[str, str, str, float]:
    """Parse annotation fields from GenomeDiff line."""
    gene_name = ''
    gene_product = ''
    annotation = ''
    frequency = 1.0

    for i in range(start_idx, len(parts)):
        field = parts[i]
        if '=' in field:
            key, value = field.split('=', 1)
            if key == 'gene_name':
                gene_name = value
            elif key == 'gene_product':
                gene_product = value
            elif key == 'annotation':
                annotation = value
            elif key == 'frequency':
                try:
                    frequency = float(value)
                except ValueError:
                    pass

    return gene_name, gene_product, annotation, frequency


def parse_snp(parts: List[str], sample_name: str) -> Optional[BreseqMutation]:
    """Parse SNP mutation."""
    try:
        evidence_id = parts[1]
        seq_id = parts[2]
        position = int(parts[3])
        new_seq = parts[4]

        gene_name, gene_product, annotation, frequency = parse_annotations(parts, 5)

        return BreseqMutation(
            sample=sample_name,
            mutation_type='SNP',
            seq_id=seq_id,
            position=position,
            size=1,
            ref='',  # breseq doesn't always provide ref
            alt=new_seq,
            gene_name=gene_name,
            gene_product=gene_product,
            annotation=annotation,
            frequency=frequency,
            evidence_id=evidence_id
        )
    except (IndexError, ValueError):
        return None


def parse_del(parts: List[str], sample_name: str) -> Optional[BreseqMutation]:
    """Parse deletion mutation."""
    try:
        evidence_id = parts[1]
        seq_id = parts[2]
        position = int(parts[3])
        size = int(parts[4])

        gene_name, gene_product, annotation, frequency = parse_annotations(parts, 5)

        return BreseqMutation(
            sample=sample_name,
            mutation_type='DEL',
            seq_id=seq_id,
            position=position,
            size=size,
            ref=f'{size}bp',
            alt='-',
            gene_name=gene_name,
            gene_product=gene_product,
            annotation=annotation,
            frequency=frequency,
            evidence_id=evidence_id
        )
    except (IndexError, ValueError):
        return None


def parse_ins(parts: List[str], sample_name: str) -> Optional[BreseqMutation]:
    """Parse insertion mutation."""
    try:
        evidence_id = parts[1]
        seq_id = parts[2]
        position = int(parts[3])
        new_seq = parts[4]

        gene_name, gene_product, annotation, frequency = parse_annotations(parts, 5)

        return BreseqMutation(
            sample=sample_name,
            mutation_type='INS',
            seq_id=seq_id,
            position=position,
            size=len(new_seq),
            ref='-',
            alt=new_seq,
            gene_name=gene_name,
            gene_product=gene_product,
            annotation=annotation,
            frequency=frequency,
            evidence_id=evidence_id
        )
    except (IndexError, ValueError):
        return None


def parse_sub(parts: List[str], sample_name: str) -> Optional[BreseqMutation]:
    """Parse substitution mutation (multiple bases)."""
    try:
        evidence_id = parts[1]
        seq_id = parts[2]
        position = int(parts[3])
        size = int(parts[4])
        new_seq = parts[5]

        gene_name, gene_product, annotation, frequency = parse_annotations(parts, 6)

        return BreseqMutation(
            sample=sample_name,
            mutation_type='SUB',
            seq_id=seq_id,
            position=position,
            size=size,
            ref=f'{size}bp',
            alt=new_seq,
            gene_name=gene_name,
            gene_product=gene_product,
            annotation=annotation,
            frequency=frequency,
            evidence_id=evidence_id
        )
    except (IndexError, ValueError):
        return None


def parse_mob(parts: List[str], sample_name: str) -> Optional[BreseqMutation]:
    """Parse mobile element insertion."""
    try:
        evidence_id = parts[1]
        seq_id = parts[2]
        position = int(parts[3])
        repeat_name = parts[4]
        strand = parts[5] if len(parts) > 5 else ''
        dup_size = int(parts[6]) if len(parts) > 6 else 0

        gene_name, gene_product, annotation, frequency = parse_annotations(parts, 7)

        return BreseqMutation(
            sample=sample_name,
            mutation_type='MOB',
            seq_id=seq_id,
            position=position,
            size=dup_size,
            ref='-',
            alt=f'{repeat_name}({strand})',
            gene_name=gene_name,
            gene_product=gene_product,
            annotation=f'IS insertion: {repeat_name}; {annotation}',
            frequency=frequency,
            evidence_id=evidence_id
        )
    except (IndexError, ValueError):
        return None


def parse_amp(parts: List[str], sample_name: str) -> Optional[BreseqMutation]:
    """Parse amplification (gene duplication)."""
    try:
        evidence_id = parts[1]
        seq_id = parts[2]
        position = int(parts[3])
        size = int(parts[4])
        copy_number = parts[5] if len(parts) > 5 else '2'

        gene_name, gene_product, annotation, frequency = parse_annotations(parts, 6)

        return BreseqMutation(
            sample=sample_name,
            mutation_type='AMP',
            seq_id=seq_id,
            position=position,
            size=size,
            ref='1x',
            alt=f'{copy_number}x',
            gene_name=gene_name,
            gene_product=gene_product,
            annotation=f'Amplification {size}bp to {copy_number} copies; {annotation}',
            frequency=frequency,
            evidence_id=evidence_id
        )
    except (IndexError, ValueError):
        return None


def find_gd_files(breseq_dir: Path) -> List[Tuple[str, Path]]:
    """Find all output.gd files in breseq results directory."""
    results = []

    for sample_dir in breseq_dir.iterdir():
        if not sample_dir.is_dir():
            continue

        # breseq output structure: sample/output/output.gd
        gd_file = sample_dir / 'output' / 'output.gd'
        if gd_file.exists():
            results.append((sample_dir.name, gd_file))

    return results


def mutations_to_dataframe(mutations: List[BreseqMutation]) -> pd.DataFrame:
    """Convert list of mutations to DataFrame."""
    if not mutations:
        return pd.DataFrame()

    data = []
    for m in mutations:
        data.append({
            'sample': m.sample,
            'type': m.mutation_type,
            'seq_id': m.seq_id,
            'position': m.position,
            'size': m.size,
            'ref': m.ref,
            'alt': m.alt,
            'gene_name': m.gene_name,
            'gene_product': m.gene_product,
            'annotation': m.annotation,
            'frequency': m.frequency,
            'evidence_id': m.evidence_id
        })

    return pd.DataFrame(data)


def compare_with_freebayes(breseq_df: pd.DataFrame,
                           freebayes_df: pd.DataFrame,
                           position_tolerance: int = 10) -> pd.DataFrame:
    """
    Compare breseq mutations with FreeBayes results.

    Returns a DataFrame with comparison results.
    """
    comparison = []

    for _, breseq_row in breseq_df.iterrows():
        sample = breseq_row['sample']
        pos = breseq_row['position']
        mut_type = breseq_row['type']

        # Look for matching mutation in FreeBayes results
        sample_fb = freebayes_df[freebayes_df['sample'] == sample]

        # Check for position match (with tolerance for indels)
        pos_col = next((c for c in ['POS', 'position', 'Position']
                       if c in freebayes_df.columns), None)

        if pos_col:
            matches = sample_fb[
                (sample_fb[pos_col] >= pos - position_tolerance) &
                (sample_fb[pos_col] <= pos + position_tolerance)
            ]

            found_in_freebayes = len(matches) > 0
            freebayes_pos = matches[pos_col].tolist() if found_in_freebayes else []
        else:
            found_in_freebayes = False
            freebayes_pos = []

        comparison.append({
            'sample': sample,
            'breseq_type': mut_type,
            'breseq_position': pos,
            'breseq_size': breseq_row['size'],
            'breseq_gene': breseq_row['gene_name'],
            'breseq_annotation': breseq_row['annotation'],
            'found_in_freebayes': found_in_freebayes,
            'freebayes_positions': str(freebayes_pos) if freebayes_pos else ''
        })

    return pd.DataFrame(comparison)


def main():
    parser = argparse.ArgumentParser(
        description='Parse breseq results and compare with FreeBayes'
    )
    parser.add_argument('--breseq-dir', '-b', required=True,
                        help='Directory containing breseq results')
    parser.add_argument('--output', '-o', required=True,
                        help='Output directory')
    parser.add_argument('--freebayes', '-f', default=None,
                        help='FreeBayes summary CSV for comparison (optional)')
    parser.add_argument('--position-tolerance', '-t', type=int, default=10,
                        help='Position tolerance for matching (default: 10bp)')

    args = parser.parse_args()

    breseq_dir = Path(args.breseq_dir)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("PARSE BRESEQ RESULTS")
    print("=" * 70)

    # Find all GD files
    print(f"\n1. Scanning for breseq results in: {breseq_dir}")
    gd_files = find_gd_files(breseq_dir)

    if not gd_files:
        print("   No breseq output.gd files found!")
        print("   Expected structure: breseq_dir/<sample>/output/output.gd")
        return

    print(f"   Found {len(gd_files)} samples:")
    for name, path in gd_files:
        print(f"      {name}")

    # Parse all GD files
    print(f"\n2. Parsing GenomeDiff files...")
    all_mutations = []

    for sample_name, gd_file in gd_files:
        mutations = parse_gd_file(gd_file, sample_name)
        print(f"   {sample_name}: {len(mutations)} mutations")
        all_mutations.extend(mutations)

    if not all_mutations:
        print("\n   No mutations found!")
        return

    # Convert to DataFrame
    df = mutations_to_dataframe(all_mutations)

    # Summary by type
    print(f"\n3. Mutation summary by type:")
    type_counts = df['type'].value_counts()
    for mut_type, count in type_counts.items():
        type_names = {
            'SNP': 'Single nucleotide polymorphism',
            'DEL': 'Deletion',
            'INS': 'Insertion',
            'SUB': 'Substitution (multi-base)',
            'MOB': 'Mobile element insertion',
            'AMP': 'Amplification/Duplication'
        }
        print(f"   {mut_type}: {count:4d} - {type_names.get(mut_type, '')}")

    # Save breseq summary
    breseq_output = output_dir / 'breseq_mutations_summary.csv'
    df.to_csv(breseq_output, index=False)
    print(f"\n4. Saved breseq summary: {breseq_output}")

    # Highlight structural variants (these are what FreeBayes often misses)
    structural = df[df['type'].isin(['DEL', 'INS', 'MOB', 'AMP', 'SUB'])]
    if len(structural) > 0:
        print(f"\n5. Structural variants (often missed by FreeBayes):")
        print("-" * 70)
        for _, row in structural.iterrows():
            print(f"   [{row['sample']}] {row['type']} at {row['position']}: "
                  f"{row['size']}bp in {row['gene_name'] or 'intergenic'}")
            if row['annotation']:
                print(f"      {row['annotation'][:60]}...")
        print("-" * 70)

        sv_output = output_dir / 'breseq_structural_variants.csv'
        structural.to_csv(sv_output, index=False)
        print(f"   Saved: {sv_output}")

    # Compare with FreeBayes if provided
    if args.freebayes:
        print(f"\n6. Comparing with FreeBayes results...")
        fb_path = Path(args.freebayes)

        if not fb_path.exists():
            print(f"   Warning: FreeBayes file not found: {fb_path}")
        else:
            fb_df = pd.read_csv(fb_path)
            print(f"   Loaded {len(fb_df)} FreeBayes mutations")

            comparison_df = compare_with_freebayes(
                df, fb_df, args.position_tolerance
            )

            # Summary
            found = comparison_df['found_in_freebayes'].sum()
            total = len(comparison_df)
            missed = total - found

            print(f"\n   Comparison summary:")
            print(f"   - Total breseq mutations: {total}")
            print(f"   - Found in FreeBayes: {found}")
            print(f"   - Missed by FreeBayes: {missed}")

            # Show missed mutations
            missed_df = comparison_df[~comparison_df['found_in_freebayes']]
            if len(missed_df) > 0:
                print(f"\n   Mutations found by breseq but missed by FreeBayes:")
                for _, row in missed_df.iterrows():
                    print(f"      [{row['sample']}] {row['breseq_type']} at {row['breseq_position']}: "
                          f"{row['breseq_gene'] or 'intergenic'}")

            # Save comparison
            comp_output = output_dir / 'breseq_freebayes_comparison.csv'
            comparison_df.to_csv(comp_output, index=False)
            print(f"\n   Saved comparison: {comp_output}")

    print(f"\n" + "=" * 70)
    print("DONE!")
    print("=" * 70)
    print(f"\nOutput files in: {output_dir}")
    print(f"  - breseq_mutations_summary.csv      : All breseq mutations")
    print(f"  - breseq_structural_variants.csv    : DEL/INS/MOB/AMP only")
    if args.freebayes:
        print(f"  - breseq_freebayes_comparison.csv   : Comparison with FreeBayes")


if __name__ == '__main__':
    main()
