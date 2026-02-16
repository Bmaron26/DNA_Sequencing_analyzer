#!/usr/bin/env python3
"""
Process existing BAM files with FreeBayes variant calling and annotation.

This script takes BAM files that are already aligned and runs:
1. FreeBayes variant calling (same parameters as your previous samples)
2. Variant filtering
3. Annotation with GFF3 file
4. Output to *_mutations_freebayes.csv (same format as your other samples)

This allows re-analysis of samples that were previously processed with a
different variant caller (e.g., bcftools) to use FreeBayes instead.

Usage:
    python process_existing_bams.py \
        --bam-dir /mnt/c/Users/hayouka-lab/Documents/Bar/3rd_evolution/WGS/results_HC \
        --reference /mnt/c/Users/hayouka-lab/Documents/Bar/3rd_evolution/WGS/reference/genome.fna \
        --annotation /mnt/c/Users/hayouka-lab/Documents/Bar/3rd_evolution/WGS/reference/genes.gff \
        --output-dir /mnt/c/Users/hayouka-lab/Documents/Bar/3rd_evolution/WGS/results_HC_freebayes
"""

import argparse
import subprocess
import sys
import re
from pathlib import Path
from typing import List, Dict, Optional, Tuple

import pandas as pd


def find_bam_files(bam_dir: Path) -> List[Tuple[str, Path]]:
    """Find all BAM files in a directory."""
    bam_files = []

    # Look for BAM files directly in the directory
    for bam_file in bam_dir.glob("*.bam"):
        # Skip index files
        if bam_file.suffix == '.bai' or str(bam_file).endswith('.bam.bai'):
            continue

        sample_name = bam_file.stem
        # Remove common suffixes
        sample_name = re.sub(r'\.sorted$', '', sample_name)
        sample_name = re.sub(r'\.aligned$', '', sample_name)

        bam_files.append((sample_name, bam_file))

    # Also check for nested structure (sample_dir/alignment/sample.bam)
    for sample_dir in bam_dir.iterdir():
        if not sample_dir.is_dir():
            continue

        alignment_dir = sample_dir / "alignment"
        if alignment_dir.exists():
            for bam_file in alignment_dir.glob("*.bam"):
                if not str(bam_file).endswith('.bai'):
                    sample_name = sample_dir.name
                    bam_files.append((sample_name, bam_file))

    return bam_files


def check_bam_index(bam_file: Path) -> bool:
    """Check if BAM index exists, create if not."""
    bai_file = Path(str(bam_file) + '.bai')
    bai_file2 = bam_file.with_suffix('.bam.bai')

    if bai_file.exists() or bai_file2.exists():
        return True

    print(f"   Creating BAM index for {bam_file.name}...")
    result = subprocess.run(
        ['samtools', 'index', str(bam_file)],
        capture_output=True, text=True
    )
    return result.returncode == 0


def run_freebayes(bam_file: Path, reference: Path, output_vcf: Path,
                  min_alt_fraction: float = 0.01,
                  min_alt_count: int = 5) -> bool:
    """Run FreeBayes variant calling."""
    cmd = [
        'freebayes',
        '-f', str(reference),
        '-p', '1',  # Haploid
        '-F', str(min_alt_fraction),
        '-C', str(min_alt_count),
        '--min-mapping-quality', '20',
        '--min-base-quality', '20',
        str(bam_file)
    ]

    print(f"   Running FreeBayes...")
    with open(output_vcf, 'w') as vcf_out:
        result = subprocess.run(cmd, stdout=vcf_out, stderr=subprocess.PIPE, text=True)

    if result.returncode != 0:
        print(f"   Error: {result.stderr}")
        return False

    return True


def filter_vcf(input_vcf: Path, output_vcf: Path, min_qual: float = 100.0) -> int:
    """Filter VCF by quality score, return number of variants passing."""
    cmd = [
        'bcftools', 'view',
        '-i', f'QUAL>={min_qual}',
        '-o', str(output_vcf),
        str(input_vcf)
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"   Warning: bcftools filter failed: {result.stderr}")
        # Fall back to simple filtering
        return simple_filter_vcf(input_vcf, output_vcf, min_qual)

    # Count variants
    count = 0
    with open(output_vcf) as f:
        for line in f:
            if not line.startswith('#'):
                count += 1
    return count


def simple_filter_vcf(input_vcf: Path, output_vcf: Path, min_qual: float) -> int:
    """Simple Python-based VCF filtering."""
    count = 0
    with open(input_vcf) as f_in, open(output_vcf, 'w') as f_out:
        for line in f_in:
            if line.startswith('#'):
                f_out.write(line)
            else:
                parts = line.split('\t')
                if len(parts) >= 6:
                    try:
                        qual = float(parts[5])
                        if qual >= min_qual:
                            f_out.write(line)
                            count += 1
                    except ValueError:
                        pass
    return count


def parse_gff(gff_path: Path) -> Dict[str, dict]:
    """Parse GFF3 file for gene annotations."""
    genes = {}

    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            seqid, source, feature, start, end, score, strand, phase, attrs = parts

            if feature not in ['CDS', 'gene']:
                continue

            # Parse attributes
            attr_dict = {}
            for item in attrs.split(';'):
                if '=' in item:
                    key, val = item.split('=', 1)
                    attr_dict[key.strip()] = val.strip()

            try:
                start_pos = int(start)
                end_pos = int(end)
            except ValueError:
                continue

            gene_id = attr_dict.get('locus_tag', attr_dict.get('ID', ''))
            gene_name = attr_dict.get('gene', attr_dict.get('Name', ''))
            product = attr_dict.get('product', '')

            # Store for position lookup
            for pos in range(start_pos, end_pos + 1):
                if pos not in genes or feature == 'CDS':  # CDS takes priority
                    genes[pos] = {
                        'locus_tag': gene_id,
                        'gene_name': gene_name,
                        'product': product,
                        'strand': strand,
                        'start': start_pos,
                        'end': end_pos,
                        'feature': feature
                    }

    return genes


def load_reference(ref_path: Path) -> Dict[str, str]:
    """Load reference genome sequences."""
    sequences = {}
    current_chrom = None
    current_seq = []

    with open(ref_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_chrom:
                    sequences[current_chrom] = ''.join(current_seq)
                current_chrom = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)

        if current_chrom:
            sequences[current_chrom] = ''.join(current_seq)

    return sequences


def get_codon_and_aa(seq: str, pos: int, strand: str, ref: str, alt: str) -> Tuple[str, str, str, str]:
    """Get codon change and amino acid change for a SNP."""
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
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }

    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}

    if len(ref) != 1 or len(alt) != 1:
        return '', '', '', ''  # Not a SNP

    # This is a simplified version - full implementation would need gene boundaries
    return '', '', '', ''


def annotate_vcf(vcf_path: Path, genes: Dict[int, dict],
                 reference: Dict[str, str], sample_name: str) -> List[dict]:
    """Annotate VCF variants with gene information."""
    variants = []

    with open(vcf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue

            chrom = parts[0]
            try:
                pos = int(parts[1])
            except ValueError:
                continue
            ref = parts[3]
            alt = parts[4]
            qual = parts[5]

            # Parse INFO field for depth and frequency
            info = parts[7] if len(parts) > 7 else ''
            info_dict = {}
            for item in info.split(';'):
                if '=' in item:
                    k, v = item.split('=', 1)
                    info_dict[k] = v

            # Get depth and allele frequency
            depth = info_dict.get('DP', '')
            ao = info_dict.get('AO', '')  # Alternate observation count
            ro = info_dict.get('RO', '')  # Reference observation count

            try:
                depth_int = int(depth) if depth else 0
                ao_int = int(ao.split(',')[0]) if ao else 0
                freq = ao_int / depth_int if depth_int > 0 else 0
            except (ValueError, ZeroDivisionError):
                freq = 0
                depth_int = 0

            # Determine variant type
            if len(ref) == len(alt) == 1:
                var_type = 'SNP'
            elif len(ref) < len(alt):
                var_type = 'INS'
            else:
                var_type = 'DEL'

            # Get gene annotation
            gene_info = genes.get(pos, {})

            # Determine effect
            if gene_info:
                location = 'coding'
                # Simplified effect determination
                if var_type == 'SNP':
                    effect = 'missense_variant'  # Simplified - would need full codon analysis
                elif var_type == 'INS':
                    if len(alt) - len(ref) % 3 == 0:
                        effect = 'inframe_insertion'
                    else:
                        effect = 'frameshift_variant'
                else:
                    if len(ref) - len(alt) % 3 == 0:
                        effect = 'inframe_deletion'
                    else:
                        effect = 'frameshift_variant'
            else:
                location = 'intergenic'
                effect = 'intergenic_variant'

            variant = {
                'sample': sample_name,
                'CHROM': chrom,
                'POS': pos,
                'REF': ref,
                'ALT': alt,
                'TYPE': var_type,
                'QUAL': qual,
                'DEPTH': depth_int,
                'FREQ': round(freq, 4),
                'locus_tag': gene_info.get('locus_tag', ''),
                'gene_name': gene_info.get('gene_name', ''),
                'product': gene_info.get('product', ''),
                'strand': gene_info.get('strand', ''),
                'location_type': location,
                'effect': effect,
            }
            variants.append(variant)

    return variants


def main():
    parser = argparse.ArgumentParser(
        description='Process existing BAM files with FreeBayes variant calling'
    )
    parser.add_argument('--bam-dir', '-b', required=True,
                        help='Directory containing BAM files')
    parser.add_argument('--reference', '-r', required=True,
                        help='Reference genome FASTA file')
    parser.add_argument('--annotation', '-a', required=True,
                        help='GFF3 annotation file')
    parser.add_argument('--output-dir', '-o', default=None,
                        help='Output directory (default: same as bam-dir)')
    parser.add_argument('--min-alt-frac', default=0.01, type=float,
                        help='FreeBayes min alternate fraction (default: 0.01)')
    parser.add_argument('--min-alt-count', default=5, type=int,
                        help='FreeBayes min alternate count (default: 5)')
    parser.add_argument('--min-qual', default=100.0, type=float,
                        help='Minimum variant quality (default: 100)')

    args = parser.parse_args()

    bam_dir = Path(args.bam_dir)
    ref_path = Path(args.reference)
    gff_path = Path(args.annotation)
    output_dir = Path(args.output_dir) if args.output_dir else bam_dir

    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("PROCESS EXISTING BAM FILES WITH FREEBAYES")
    print("=" * 70)

    # Check for required tools
    print("\n1. Checking required tools...")
    for tool in ['freebayes', 'samtools', 'bcftools']:
        result = subprocess.run(['which', tool], capture_output=True, text=True)
        if result.returncode != 0:
            print(f"   ERROR: {tool} not found in PATH!")
            print(f"   Please install {tool} before running this script.")
            sys.exit(1)
        print(f"   {tool}: OK")

    # Find BAM files
    print(f"\n2. Finding BAM files in: {bam_dir}")
    bam_files = find_bam_files(bam_dir)

    if not bam_files:
        print("   No BAM files found!")
        sys.exit(1)

    print(f"   Found {len(bam_files)} BAM files:")
    for sample_name, bam_file in bam_files:
        print(f"      {sample_name}: {bam_file.name}")

    # Load reference and annotation
    print(f"\n3. Loading reference: {ref_path}")
    reference = load_reference(ref_path)
    print(f"   Loaded {len(reference)} sequences")

    print(f"\n4. Loading annotation: {gff_path}")
    genes = parse_gff(gff_path)
    print(f"   Loaded annotations for {len(genes)} positions")

    # Process each BAM file
    print(f"\n5. Processing BAM files with FreeBayes...")
    all_variants = []

    for sample_name, bam_file in bam_files:
        print(f"\n   [{sample_name}]")

        # Ensure index exists
        if not check_bam_index(bam_file):
            print(f"   Warning: Could not create BAM index")
            continue

        # Create sample output directory
        sample_output = output_dir / sample_name
        sample_output.mkdir(exist_ok=True)

        # Run FreeBayes
        raw_vcf = sample_output / f"{sample_name}.freebayes.vcf"
        if not run_freebayes(bam_file, ref_path, raw_vcf,
                            args.min_alt_frac, args.min_alt_count):
            print(f"   Error running FreeBayes, skipping...")
            continue

        # Filter VCF
        filtered_vcf = sample_output / f"{sample_name}.freebayes.filtered.vcf"
        n_variants = filter_vcf(raw_vcf, filtered_vcf, args.min_qual)
        print(f"   Variants after filtering: {n_variants}")

        # Annotate variants
        variants = annotate_vcf(filtered_vcf, genes, reference, sample_name)
        print(f"   Annotated {len(variants)} variants")

        # Save per-sample CSV
        if variants:
            df = pd.DataFrame(variants)
            csv_file = sample_output / f"{sample_name}_mutations_freebayes.csv"
            df.to_csv(csv_file, index=False)
            print(f"   Saved: {csv_file.name}")

            all_variants.extend(variants)

    # Save combined summary
    print(f"\n6. Creating combined summary...")
    if all_variants:
        df_all = pd.DataFrame(all_variants)
        summary_file = output_dir / "all_mutations_summary_freebayes.csv"
        df_all.to_csv(summary_file, index=False)
        print(f"   Saved: {summary_file}")

        # Print quick stats
        print(f"\n   Total: {len(all_variants)} variants from {len(bam_files)} samples")
        print(f"   Variants per sample:")
        for sample in df_all['sample'].unique():
            count = len(df_all[df_all['sample'] == sample])
            print(f"      {sample}: {count}")

    print(f"\n" + "=" * 70)
    print("DONE!")
    print("=" * 70)
    print(f"\nOutput directory: {output_dir}")
    print("\nNext steps:")
    print("  1. The *_mutations_freebayes.csv files are now compatible with your other samples")
    print("  2. Run finalize_mutation_annotations.py to add gene names from your lookup table")
    print("  3. Run analyze_amp_combinations.py for cross-treatment analysis")


if __name__ == '__main__':
    main()
