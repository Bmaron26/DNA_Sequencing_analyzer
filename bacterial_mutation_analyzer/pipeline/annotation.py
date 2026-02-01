"""
Variant annotation module for functional annotation of mutations.

Provides:
- Gene-level annotation (which gene is affected)
- Effect prediction (synonymous, nonsynonymous, etc.)
- Protein change annotation (amino acid changes)
- Intergenic variant annotation
"""

import os
from pathlib import Path
from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional, Tuple
import logging
from collections import defaultdict

from ..config import AnnotationConfig
from ..utils.file_handlers import GffParser, GenbankParser, VcfParser, GffFeature

logger = logging.getLogger(__name__)


# Codon table for standard genetic code
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

# Amino acid properties for effect classification
AA_PROPERTIES = {
    'A': 'hydrophobic', 'V': 'hydrophobic', 'L': 'hydrophobic',
    'I': 'hydrophobic', 'M': 'hydrophobic', 'F': 'hydrophobic',
    'W': 'hydrophobic', 'P': 'hydrophobic',
    'G': 'small', 'S': 'polar', 'T': 'polar', 'C': 'polar',
    'Y': 'polar', 'N': 'polar', 'Q': 'polar',
    'K': 'positive', 'R': 'positive', 'H': 'positive',
    'D': 'negative', 'E': 'negative',
    '*': 'stop',
}


@dataclass
class AnnotatedVariant:
    """Annotated variant with functional information."""
    chromosome: str
    position: int
    reference: str
    alternative: str
    variant_type: str
    quality: float
    depth: int
    allele_frequency: float

    # Annotation fields
    gene_id: Optional[str] = None
    gene_name: Optional[str] = None
    locus_tag: Optional[str] = None
    product: Optional[str] = None
    feature_type: Optional[str] = None
    strand: Optional[str] = None

    # Effect fields
    effect: str = "unknown"
    effect_impact: str = "unknown"
    codon_change: Optional[str] = None
    amino_acid_change: Optional[str] = None
    amino_acid_position: Optional[int] = None
    protein_position: Optional[str] = None

    # Location fields
    distance_to_feature: int = 0
    location_type: str = "unknown"  # coding, intergenic, upstream, downstream

    # Additional info
    info: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for export."""
        return {
            'chromosome': self.chromosome,
            'position': self.position,
            'reference': self.reference,
            'alternative': self.alternative,
            'variant_type': self.variant_type,
            'quality': self.quality,
            'depth': self.depth,
            'allele_frequency': round(self.allele_frequency, 4),
            'gene_id': self.gene_id or '',
            'gene_name': self.gene_name or '',
            'locus_tag': self.locus_tag or '',
            'product': self.product or '',
            'feature_type': self.feature_type or '',
            'strand': self.strand or '',
            'effect': self.effect,
            'effect_impact': self.effect_impact,
            'codon_change': self.codon_change or '',
            'amino_acid_change': self.amino_acid_change or '',
            'amino_acid_position': self.amino_acid_position or '',
            'location_type': self.location_type,
            'distance_to_feature': self.distance_to_feature,
        }


class VariantAnnotator:
    """
    Annotates variants with gene and functional information.
    """

    def __init__(self, config: Optional[AnnotationConfig] = None,
                 output_dir: str = "annotation"):
        self.config = config or AnnotationConfig()
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.features: List[GffFeature] = []
        self.features_by_location: Dict[str, List[GffFeature]] = defaultdict(list)
        self.reference_sequences: Dict[str, str] = {}

    def load_annotations(self, annotation_file: str) -> int:
        """
        Load gene annotations from GFF or GenBank file.

        Args:
            annotation_file: Path to annotation file

        Returns:
            Number of features loaded
        """
        path = Path(annotation_file)
        suffix = ''.join(path.suffixes).lower()

        if any(ext in suffix for ext in ['.gff', '.gff3']):
            parser = GffParser(annotation_file)
            self.features = parser.parse()
        elif any(ext in suffix for ext in ['.gbk', '.gb', '.genbank']):
            parser = GenbankParser(annotation_file)
            self.features = parser.parse()
        else:
            raise ValueError(f"Unsupported annotation format: {suffix}")

        # Index features by location for fast lookup
        for feature in self.features:
            self.features_by_location[feature.seqid].append(feature)

        # Sort features by position for efficient searching
        for seqid in self.features_by_location:
            self.features_by_location[seqid].sort(key=lambda f: f.start)

        logger.info(f"Loaded {len(self.features)} features from {annotation_file}")
        return len(self.features)

    def load_reference(self, reference_file: str) -> None:
        """Load reference genome for codon analysis."""
        from ..utils.file_handlers import FastaReader

        reader = FastaReader(reference_file)
        for record in reader:
            self.reference_sequences[record.id] = record.sequence.upper()

        logger.info(f"Loaded {len(self.reference_sequences)} reference sequences")

    def annotate_variants(self, vcf_file: str,
                         annotation_file: Optional[str] = None,
                         reference_file: Optional[str] = None) -> List[AnnotatedVariant]:
        """
        Annotate variants from VCF file.

        Args:
            vcf_file: Path to VCF file
            annotation_file: Path to annotation file (optional if already loaded)
            reference_file: Path to reference FASTA (optional, for protein effects)

        Returns:
            List of annotated variants
        """
        if annotation_file:
            self.load_annotations(annotation_file)

        if reference_file:
            self.load_reference(reference_file)

        # Parse VCF
        parser = VcfParser(vcf_file)
        variants = parser.parse()

        annotated = []
        for variant in variants:
            ann_variant = self._annotate_single(variant)
            annotated.append(ann_variant)

        logger.info(f"Annotated {len(annotated)} variants")
        return annotated

    def _annotate_single(self, variant) -> AnnotatedVariant:
        """Annotate a single variant."""
        ann = AnnotatedVariant(
            chromosome=variant.chrom,
            position=variant.pos,
            reference=variant.ref,
            alternative=variant.alt,
            variant_type=variant.variant_type,
            quality=variant.qual,
            depth=variant.depth,
            allele_frequency=variant.allele_frequency,
        )

        # Find overlapping features
        overlapping = self._find_overlapping_features(variant.chrom, variant.pos)

        if overlapping:
            # Use the most specific feature (prefer CDS over gene)
            feature = self._select_best_feature(overlapping)
            ann = self._annotate_with_feature(ann, feature, variant)
        else:
            # Find nearest feature for intergenic variants
            nearest, distance = self._find_nearest_feature(variant.chrom, variant.pos)
            if nearest:
                ann.gene_id = nearest.gene_id
                ann.gene_name = nearest.gene_name
                ann.locus_tag = nearest.locus_tag
                ann.product = nearest.product
                ann.distance_to_feature = distance

                if variant.pos < nearest.start:
                    ann.location_type = "upstream"
                    ann.effect = "upstream_gene_variant"
                else:
                    ann.location_type = "downstream"
                    ann.effect = "downstream_gene_variant"

                ann.effect_impact = "MODIFIER"
            else:
                ann.location_type = "intergenic"
                ann.effect = "intergenic_region"
                ann.effect_impact = "MODIFIER"

        return ann

    def _find_overlapping_features(self, chrom: str, position: int) -> List[GffFeature]:
        """Find all features overlapping a position."""
        overlapping = []
        for feature in self.features_by_location.get(chrom, []):
            if feature.start <= position <= feature.end:
                overlapping.append(feature)
        return overlapping

    def _find_nearest_feature(self, chrom: str, position: int) -> Tuple[Optional[GffFeature], int]:
        """Find the nearest gene/CDS feature."""
        features = self.features_by_location.get(chrom, [])

        nearest = None
        min_distance = float('inf')

        for feature in features:
            if feature.feature_type not in ['gene', 'CDS', 'mRNA']:
                continue

            if feature.start <= position <= feature.end:
                return feature, 0

            if position < feature.start:
                distance = feature.start - position
            else:
                distance = position - feature.end

            if distance < min_distance:
                min_distance = distance
                nearest = feature

            # Early termination for sorted features
            if feature.start > position and distance > min_distance:
                break

        return nearest, int(min_distance) if nearest else -1

    def _select_best_feature(self, features: List[GffFeature]) -> GffFeature:
        """Select the most informative feature from overlapping ones."""
        # Priority: CDS > mRNA > gene > other
        priority = {'CDS': 0, 'mRNA': 1, 'gene': 2}

        features.sort(key=lambda f: priority.get(f.feature_type, 3))
        return features[0]

    def _annotate_with_feature(self, ann: AnnotatedVariant, feature: GffFeature,
                               variant) -> AnnotatedVariant:
        """Annotate variant with overlapping feature information."""
        ann.gene_id = feature.gene_id
        ann.gene_name = feature.gene_name
        ann.locus_tag = feature.locus_tag
        ann.product = feature.product
        ann.feature_type = feature.feature_type
        ann.strand = feature.strand
        ann.location_type = "coding" if feature.feature_type == 'CDS' else "genic"

        # Calculate protein effect for CDS features
        if feature.feature_type == 'CDS' and self.config.protein_effect:
            ann = self._calculate_protein_effect(ann, feature, variant)
        elif feature.feature_type == 'CDS':
            ann.effect = "coding_sequence_variant"
            ann.effect_impact = "MODERATE"
        else:
            ann.effect = "genic_variant"
            ann.effect_impact = "LOW"

        return ann

    def _calculate_protein_effect(self, ann: AnnotatedVariant, feature: GffFeature,
                                  variant) -> AnnotatedVariant:
        """Calculate the protein-level effect of a variant."""
        chrom = variant.chrom
        pos = variant.pos
        ref = variant.ref
        alt = variant.alt

        # Get reference sequence
        ref_seq = self.reference_sequences.get(chrom)
        if not ref_seq:
            ann.effect = "coding_sequence_variant"
            ann.effect_impact = "MODERATE"
            return ann

        # Handle SNPs
        if len(ref) == 1 and len(alt) == 1:
            ann = self._annotate_snp(ann, feature, pos, ref, alt, ref_seq)
        elif len(ref) < len(alt):
            # Insertion
            insert_len = len(alt) - len(ref)
            if insert_len % 3 == 0:
                ann.effect = "inframe_insertion"
                ann.effect_impact = "MODERATE"
            else:
                ann.effect = "frameshift_variant"
                ann.effect_impact = "HIGH"
        elif len(ref) > len(alt):
            # Deletion
            del_len = len(ref) - len(alt)
            if del_len % 3 == 0:
                ann.effect = "inframe_deletion"
                ann.effect_impact = "MODERATE"
            else:
                ann.effect = "frameshift_variant"
                ann.effect_impact = "HIGH"
        else:
            ann.effect = "complex_variant"
            ann.effect_impact = "MODERATE"

        return ann

    def _annotate_snp(self, ann: AnnotatedVariant, feature: GffFeature,
                      pos: int, ref: str, alt: str, ref_seq: str) -> AnnotatedVariant:
        """Annotate a single nucleotide polymorphism."""
        # Calculate position within CDS
        if feature.strand == '+':
            cds_pos = pos - feature.start
        else:
            cds_pos = feature.end - pos

        # Calculate codon position (0-based)
        codon_pos = cds_pos % 3
        codon_start = pos - codon_pos if feature.strand == '+' else pos + codon_pos - 2

        # Extract reference codon
        try:
            if feature.strand == '+':
                ref_codon = ref_seq[codon_start - 1:codon_start + 2]
                alt_codon = ref_codon[:codon_pos] + alt + ref_codon[codon_pos + 1:]
            else:
                ref_codon = self._reverse_complement(ref_seq[codon_start - 1:codon_start + 2])
                codon_pos_rev = 2 - codon_pos
                alt_codon = ref_codon[:codon_pos_rev] + self._complement(alt) + ref_codon[codon_pos_rev + 1:]

            # Translate codons
            ref_aa = CODON_TABLE.get(ref_codon.upper(), 'X')
            alt_aa = CODON_TABLE.get(alt_codon.upper(), 'X')

            ann.codon_change = f"{ref_codon}>{alt_codon}"
            ann.amino_acid_position = (cds_pos // 3) + 1
            ann.amino_acid_change = f"{ref_aa}{ann.amino_acid_position}{alt_aa}"
            ann.protein_position = f"p.{ref_aa}{ann.amino_acid_position}{alt_aa}"

            if ref_aa == alt_aa:
                ann.effect = "synonymous_variant"
                ann.effect_impact = "LOW"
            elif alt_aa == '*':
                ann.effect = "stop_gained"
                ann.effect_impact = "HIGH"
            elif ref_aa == '*':
                ann.effect = "stop_lost"
                ann.effect_impact = "HIGH"
            elif ref_aa == 'M' and ann.amino_acid_position == 1:
                ann.effect = "start_lost"
                ann.effect_impact = "HIGH"
            else:
                # Check if conservative or non-conservative
                if AA_PROPERTIES.get(ref_aa) == AA_PROPERTIES.get(alt_aa):
                    ann.effect = "missense_variant_conservative"
                    ann.effect_impact = "MODERATE"
                else:
                    ann.effect = "missense_variant"
                    ann.effect_impact = "MODERATE"

        except (IndexError, KeyError) as e:
            logger.debug(f"Could not determine protein effect: {e}")
            ann.effect = "coding_sequence_variant"
            ann.effect_impact = "MODERATE"

        return ann

    def _reverse_complement(self, seq: str) -> str:
        """Get reverse complement of sequence."""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        return ''.join(complement.get(b, 'N') for b in reversed(seq.upper()))

    def _complement(self, base: str) -> str:
        """Get complement of a single base."""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        return complement.get(base.upper(), 'N')

    def get_summary(self, annotated_variants: List[AnnotatedVariant]) -> Dict[str, Any]:
        """Generate summary statistics of annotated variants."""
        summary = {
            'total_variants': len(annotated_variants),
            'by_effect': defaultdict(int),
            'by_impact': defaultdict(int),
            'by_type': defaultdict(int),
            'by_location': defaultdict(int),
            'genes_affected': set(),
            'high_impact_variants': [],
        }

        for v in annotated_variants:
            summary['by_effect'][v.effect] += 1
            summary['by_impact'][v.effect_impact] += 1
            summary['by_type'][v.variant_type] += 1
            summary['by_location'][v.location_type] += 1

            if v.gene_id:
                summary['genes_affected'].add(v.gene_id)

            if v.effect_impact == 'HIGH':
                summary['high_impact_variants'].append(v.to_dict())

        summary['genes_affected'] = list(summary['genes_affected'])
        summary['by_effect'] = dict(summary['by_effect'])
        summary['by_impact'] = dict(summary['by_impact'])
        summary['by_type'] = dict(summary['by_type'])
        summary['by_location'] = dict(summary['by_location'])

        return summary
