#!/usr/bin/env python3
"""
Convert VCF files to annotated mutation CSV files.

This script parses VCF files from FreeBayes, annotates mutations using a GFF3 file,
and creates per-sample mutation CSV files compatible with downstream analysis scripts.
"""

import argparse
import os
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass
import csv


@dataclass
class GeneFeature:
    """Represents a gene/CDS feature from GFF."""
    seqid: str
    feature_type: str
    start: int
    end: int
    strand: str
    gene_name: str
    locus_tag: str
    product: str
    attributes: Dict[str, str]


@dataclass
class Mutation:
    """Represents a mutation from VCF."""
    chrom: str
    pos: int
    ref: str
    alt: str
    qual: float
    info: Dict[str, str]
    # Annotation fields
    gene: Optional[str] = None
    locus_tag: Optional[str] = None
    product: Optional[str] = None
    effect: Optional[str] = None
    aa_change: Optional[str] = None
    location_type: Optional[str] = None
    codon_pos: Optional[int] = None


def parse_gff(gff_path: str) -> Tuple[List[GeneFeature], Dict[str, List[GeneFeature]]]:
    """Parse GFF3 file and return list of gene features."""
    features = []
    features_by_seqid = {}

    print(f"   Parsing GFF file: {gff_path}")

    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            seqid, source, ftype, start, end, score, strand, phase, attrs_str = parts

            # Only process gene and CDS features
            if ftype not in ['gene', 'CDS']:
                continue

            # Parse attributes
            attrs = {}
            for attr in attrs_str.split(';'):
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    attrs[key] = value

            # Extract gene name and locus tag
            gene_name = attrs.get('gene', attrs.get('Name', ''))
            locus_tag = attrs.get('locus_tag', attrs.get('ID', ''))
            product = attrs.get('product', '')

            feature = GeneFeature(
                seqid=seqid,
                feature_type=ftype,
                start=int(start),
                end=int(end),
                strand=strand,
                gene_name=gene_name,
                locus_tag=locus_tag,
                product=product,
                attributes=attrs
            )

            features.append(feature)

            if seqid not in features_by_seqid:
                features_by_seqid[seqid] = []
            features_by_seqid[seqid].append(feature)

    # Sort features by position for each sequence
    for seqid in features_by_seqid:
        features_by_seqid[seqid].sort(key=lambda x: x.start)

    print(f"   Found {len(features)} gene/CDS features")
    return features, features_by_seqid


def parse_vcf(vcf_path: str) -> List[Mutation]:
    """Parse VCF file and return list of mutations."""
    mutations = []

    with open(vcf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue

            chrom = parts[0]
            pos = int(parts[1])
            ref = parts[3]
            alt = parts[4]

            try:
                qual = float(parts[5]) if parts[5] != '.' else 0.0
            except ValueError:
                qual = 0.0

            # Parse INFO field
            info = {}
            if len(parts) > 7 and parts[7] != '.':
                for item in parts[7].split(';'):
                    if '=' in item:
                        key, value = item.split('=', 1)
                        info[key] = value
                    else:
                        info[item] = 'True'

            # Handle multiple alt alleles (take first one)
            if ',' in alt:
                alt = alt.split(',')[0]

            mutation = Mutation(
                chrom=chrom,
                pos=pos,
                ref=ref,
                alt=alt,
                qual=qual,
                info=info
            )
            mutations.append(mutation)

    return mutations


def load_reference(ref_path: str) -> Dict[str, str]:
    """Load reference genome from FASTA file."""
    sequences = {}
    current_seqid = None
    current_seq = []

    with open(ref_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seqid:
                    sequences[current_seqid] = ''.join(current_seq)
                current_seqid = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)

        if current_seqid:
            sequences[current_seqid] = ''.join(current_seq)

    return sequences


def get_codon_and_aa(seq: str, pos: int, strand: str, ref: str, alt: str) -> Tuple[Optional[int], Optional[str]]:
    """Calculate codon position and amino acid change."""
    # This is a simplified version - full implementation would need
    # to handle frameshifts, indels, etc.

    codon_table = {
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
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }

    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}

    # Only handle SNPs for now
    if len(ref) != 1 or len(alt) != 1:
        return None, None

    codon_pos = pos % 3
    if codon_pos == 0:
        codon_pos = 3

    return codon_pos, None  # Simplified - would need gene coordinates for full AA calculation


def annotate_mutation(mutation: Mutation, features_by_seqid: Dict[str, List[GeneFeature]],
                      reference: Dict[str, str]) -> Mutation:
    """Annotate a mutation with gene information."""

    pos = mutation.pos
    chrom = mutation.chrom

    if chrom not in features_by_seqid:
        mutation.location_type = 'intergenic'
        return mutation

    features = features_by_seqid[chrom]

    # Find overlapping features
    overlapping_genes = []
    overlapping_cds = []

    for feature in features:
        if feature.start <= pos <= feature.end:
            if feature.feature_type == 'gene':
                overlapping_genes.append(feature)
            elif feature.feature_type == 'CDS':
                overlapping_cds.append(feature)

    if overlapping_cds:
        # Mutation is in a coding region
        cds = overlapping_cds[0]
        mutation.gene = cds.gene_name or cds.locus_tag
        mutation.locus_tag = cds.locus_tag
        mutation.product = cds.product
        mutation.location_type = 'coding'

        # Determine if synonymous or non-synonymous
        if len(mutation.ref) == 1 and len(mutation.alt) == 1:
            # SNP - try to determine effect
            codon_pos, aa_change = get_codon_and_aa(
                reference.get(chrom, ''),
                pos,
                cds.strand,
                mutation.ref,
                mutation.alt
            )
            mutation.codon_pos = codon_pos

            # For now, mark as missense (would need full codon analysis)
            mutation.effect = 'missense_variant'
        elif len(mutation.ref) != len(mutation.alt):
            # Indel
            if (len(mutation.alt) - len(mutation.ref)) % 3 == 0:
                mutation.effect = 'inframe_indel'
            else:
                mutation.effect = 'frameshift_variant'

    elif overlapping_genes:
        # In a gene but not in CDS (e.g., intron, UTR)
        gene = overlapping_genes[0]
        mutation.gene = gene.gene_name or gene.locus_tag
        mutation.locus_tag = gene.locus_tag
        mutation.product = gene.product
        mutation.location_type = 'genic'
        mutation.effect = 'non_coding_variant'
    else:
        # Intergenic
        mutation.location_type = 'intergenic'
        mutation.effect = 'intergenic_variant'

        # Find nearest genes
        nearest_upstream = None
        nearest_downstream = None

        for feature in features:
            if feature.feature_type != 'gene':
                continue
            if feature.end < pos:
                if nearest_upstream is None or feature.end > nearest_upstream.end:
                    nearest_upstream = feature
            elif feature.start > pos:
                if nearest_downstream is None or feature.start < nearest_downstream.start:
                    nearest_downstream = feature

        # Set gene name as "upstream_gene/downstream_gene"
        if nearest_upstream and nearest_downstream:
            up_name = nearest_upstream.gene_name or nearest_upstream.locus_tag
            down_name = nearest_downstream.gene_name or nearest_downstream.locus_tag
            mutation.gene = f"{up_name}/{down_name}"

    return mutation


def write_mutations_csv(mutations: List[Mutation], output_path: str, sample_name: str):
    """Write mutations to CSV file."""

    fieldnames = [
        'sample', 'chrom', 'pos', 'ref', 'alt', 'qual',
        'gene', 'locus_tag', 'product', 'effect', 'aa_change',
        'location_type', 'codon_pos'
    ]

    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for mut in mutations:
            row = {
                'sample': sample_name,
                'chrom': mut.chrom,
                'pos': mut.pos,
                'ref': mut.ref,
                'alt': mut.alt,
                'qual': mut.qual,
                'gene': mut.gene or '',
                'locus_tag': mut.locus_tag or '',
                'product': mut.product or '',
                'effect': mut.effect or '',
                'aa_change': mut.aa_change or '',
                'location_type': mut.location_type or '',
                'codon_pos': mut.codon_pos or ''
            }
            writer.writerow(row)


def process_vcf_directory(vcf_dir: str, gff_path: str, ref_path: str, output_dir: str,
                          exclude_genes: List[str] = None, exclude_positions: List[int] = None):
    """Process all VCF files in a directory."""

    print("=" * 70)
    print("VCF TO MUTATIONS CSV CONVERTER")
    print("=" * 70)

    # Parse GFF
    print("\n1. Parsing GFF annotation file...")
    features, features_by_seqid = parse_gff(gff_path)

    # Load reference
    print("\n2. Loading reference genome...")
    reference = load_reference(ref_path)
    print(f"   Loaded {len(reference)} sequence(s)")

    # Find VCF files
    print("\n3. Processing VCF files...")
    vcf_dir = Path(vcf_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    vcf_files = list(vcf_dir.glob("*.vcf"))
    print(f"   Found {len(vcf_files)} VCF files")

    if not vcf_files:
        print("   No VCF files found!")
        return

    all_mutations = []

    for vcf_path in sorted(vcf_files):
        sample_name = vcf_path.stem
        print(f"\n   Processing: {sample_name}")

        # Parse VCF
        mutations = parse_vcf(str(vcf_path))
        print(f"      Raw variants: {len(mutations)}")

        # Annotate mutations
        for mut in mutations:
            annotate_mutation(mut, features_by_seqid, reference)

        # Filter excluded genes and positions
        if exclude_genes or exclude_positions:
            original_count = len(mutations)
            filtered = []
            for mut in mutations:
                skip = False
                if exclude_genes and mut.gene:
                    if any(eg.lower() in mut.gene.lower() for eg in exclude_genes):
                        skip = True
                if exclude_positions and mut.pos in exclude_positions:
                    skip = True
                if not skip:
                    filtered.append(mut)
            mutations = filtered
            print(f"      After filtering: {len(mutations)} (removed {original_count - len(mutations)})")

        # Write to CSV
        output_path = output_dir / f"{sample_name}_mutations.csv"
        write_mutations_csv(mutations, str(output_path), sample_name)
        print(f"      Wrote: {output_path.name}")

        all_mutations.extend(mutations)

    # Write summary
    print("\n4. Writing summary files...")

    # All mutations summary
    summary_path = output_dir / "all_mutations_summary.csv"
    write_mutations_csv(all_mutations, str(summary_path), "all")
    print(f"   Wrote: {summary_path.name} ({len(all_mutations)} mutations)")

    # Mutation counts per sample
    counts_path = output_dir / "mutation_counts.csv"
    sample_counts = {}
    for mut in all_mutations:
        # Get sample from the mutation's source (need to track this)
        pass

    print("\n" + "=" * 70)
    print("DONE!")
    print("=" * 70)


def main():
    parser = argparse.ArgumentParser(
        description='Convert VCF files to annotated mutation CSV files'
    )
    parser.add_argument('--vcf-dir', required=True,
                        help='Directory containing VCF files')
    parser.add_argument('--gff', required=True,
                        help='GFF3 annotation file')
    parser.add_argument('--ref', required=True,
                        help='Reference FASTA file')
    parser.add_argument('--output', required=True,
                        help='Output directory for CSV files')
    parser.add_argument('--exclude-genes', nargs='+', default=[],
                        help='Genes to exclude from output')
    parser.add_argument('--exclude-positions', nargs='+', type=int, default=[],
                        help='Positions to exclude from output')

    args = parser.parse_args()

    process_vcf_directory(
        vcf_dir=args.vcf_dir,
        gff_path=args.gff,
        ref_path=args.ref,
        output_dir=args.output,
        exclude_genes=args.exclude_genes,
        exclude_positions=args.exclude_positions
    )


if __name__ == '__main__':
    main()
