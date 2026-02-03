"""
Gene mapper module for cross-referencing genes between different annotations.

Useful for mapping mutations from one reference to another (e.g., de novo assembly
to a well-annotated reference like SH1000).
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass
from collections import defaultdict

from ..utils.file_handlers import GffParser, GffFeature

logger = logging.getLogger(__name__)


@dataclass
class GeneMapping:
    """Represents a mapping between genes in two annotations."""
    source_locus_tag: str
    source_gene_name: str
    source_product: str
    target_locus_tag: str
    target_gene_name: str
    target_product: str
    match_type: str  # 'gene_name', 'product', 'locus_tag'
    confidence: float  # 0-1, how confident the match is


class GeneMapper:
    """
    Maps genes between two different genome annotations.

    Useful for cross-referencing mutations found with one reference
    to gene annotations in another (better annotated) reference.
    """

    def __init__(self):
        self.source_features: List[GffFeature] = []
        self.target_features: List[GffFeature] = []
        self.gene_mappings: Dict[str, GeneMapping] = {}  # source_locus_tag -> mapping
        self.product_index: Dict[str, List[GffFeature]] = defaultdict(list)
        self.gene_name_index: Dict[str, List[GffFeature]] = defaultdict(list)

    def load_source_annotation(self, gff_path: str) -> int:
        """Load the source annotation (your original reference)."""
        parser = GffParser(gff_path)
        self.source_features = [f for f in parser.parse()
                                if f.feature_type in ('gene', 'CDS')]
        logger.info(f"Loaded {len(self.source_features)} source features from {gff_path}")
        return len(self.source_features)

    def load_target_annotation(self, gff_path: str) -> int:
        """Load the target annotation (e.g., SH1000 for lookup)."""
        parser = GffParser(gff_path)
        self.target_features = [f for f in parser.parse()
                                if f.feature_type in ('gene', 'CDS')]

        # Build indexes for fast lookup
        self.product_index.clear()
        self.gene_name_index.clear()

        for feature in self.target_features:
            # Index by product (normalized)
            if feature.product:
                key = self._normalize_product(feature.product)
                self.product_index[key].append(feature)

            # Index by gene name
            if feature.gene_name:
                key = feature.gene_name.lower()
                self.gene_name_index[key].append(feature)

        logger.info(f"Loaded {len(self.target_features)} target features from {gff_path}")
        return len(self.target_features)

    def _normalize_product(self, product: str) -> str:
        """Normalize product string for matching."""
        # Remove common variations
        normalized = product.lower()
        # Remove "putative", "probable", "hypothetical" prefixes
        for prefix in ['putative ', 'probable ', 'hypothetical ']:
            normalized = normalized.replace(prefix, '')
        # Remove whitespace variations
        normalized = ' '.join(normalized.split())
        return normalized

    def build_mappings(self) -> int:
        """
        Build mappings between source and target annotations.

        Returns:
            Number of successful mappings
        """
        self.gene_mappings.clear()
        matched = 0

        for source_feature in self.source_features:
            mapping = self._find_best_match(source_feature)
            if mapping:
                source_key = source_feature.locus_tag or source_feature.gene_id
                if source_key:
                    self.gene_mappings[source_key] = mapping
                    matched += 1

        logger.info(f"Created {matched} gene mappings out of {len(self.source_features)} source features")
        return matched

    def _find_best_match(self, source_feature: GffFeature) -> Optional[GeneMapping]:
        """Find the best matching target feature for a source feature."""
        candidates = []

        # Try matching by gene name first (highest confidence)
        if source_feature.gene_name:
            key = source_feature.gene_name.lower()
            if key in self.gene_name_index:
                for target in self.gene_name_index[key]:
                    candidates.append((target, 'gene_name', 1.0))

        # Try matching by product
        if source_feature.product:
            key = self._normalize_product(source_feature.product)
            if key in self.product_index:
                for target in self.product_index[key]:
                    # Lower confidence for product match
                    candidates.append((target, 'product', 0.8))

            # Try partial product match
            for product_key, targets in self.product_index.items():
                if len(key) > 10 and (key in product_key or product_key in key):
                    for target in targets:
                        candidates.append((target, 'product_partial', 0.5))

        if not candidates:
            return None

        # Sort by confidence and take the best
        candidates.sort(key=lambda x: x[2], reverse=True)
        best_target, match_type, confidence = candidates[0]

        return GeneMapping(
            source_locus_tag=source_feature.locus_tag or '',
            source_gene_name=source_feature.gene_name or '',
            source_product=source_feature.product or '',
            target_locus_tag=best_target.locus_tag or '',
            target_gene_name=best_target.gene_name or '',
            target_product=best_target.product or '',
            match_type=match_type,
            confidence=confidence
        )

    def get_target_locus_tag(self, source_locus_tag: str) -> Optional[str]:
        """Get the target locus tag for a source locus tag."""
        mapping = self.gene_mappings.get(source_locus_tag)
        return mapping.target_locus_tag if mapping else None

    def get_target_gene_name(self, source_locus_tag: str) -> Optional[str]:
        """Get the target gene name for a source locus tag."""
        mapping = self.gene_mappings.get(source_locus_tag)
        return mapping.target_gene_name if mapping else None

    def get_mapping(self, source_locus_tag: str) -> Optional[GeneMapping]:
        """Get the full mapping for a source locus tag."""
        return self.gene_mappings.get(source_locus_tag)

    def add_target_annotations_to_variants(self, variants: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Add target reference annotations to a list of variants.

        Args:
            variants: List of variant dictionaries from pipeline results

        Returns:
            Variants with added target annotation columns
        """
        for variant in variants:
            source_locus = variant.get('locus_tag', '')

            mapping = self.gene_mappings.get(source_locus)

            if mapping:
                variant['target_locus_tag'] = mapping.target_locus_tag
                variant['target_gene_name'] = mapping.target_gene_name
                variant['target_product'] = mapping.target_product
                variant['mapping_confidence'] = mapping.confidence
                variant['mapping_type'] = mapping.match_type
            else:
                # Try matching by gene name from variant
                gene_name = variant.get('gene_name', '').lower()
                if gene_name and gene_name in self.gene_name_index:
                    target = self.gene_name_index[gene_name][0]
                    variant['target_locus_tag'] = target.locus_tag or ''
                    variant['target_gene_name'] = target.gene_name or ''
                    variant['target_product'] = target.product or ''
                    variant['mapping_confidence'] = 0.9
                    variant['mapping_type'] = 'gene_name_direct'
                else:
                    variant['target_locus_tag'] = ''
                    variant['target_gene_name'] = ''
                    variant['target_product'] = ''
                    variant['mapping_confidence'] = 0
                    variant['mapping_type'] = 'no_match'

        return variants

    def export_mapping_table(self, output_path: str) -> str:
        """Export the gene mapping table to CSV."""
        import csv

        with open(output_path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow([
                'source_locus_tag', 'source_gene', 'source_product',
                'target_locus_tag', 'target_gene', 'target_product',
                'match_type', 'confidence'
            ])

            for source_tag, mapping in self.gene_mappings.items():
                writer.writerow([
                    mapping.source_locus_tag,
                    mapping.source_gene_name,
                    mapping.source_product,
                    mapping.target_locus_tag,
                    mapping.target_gene_name,
                    mapping.target_product,
                    mapping.match_type,
                    mapping.confidence
                ])

        logger.info(f"Exported {len(self.gene_mappings)} mappings to {output_path}")
        return output_path


def create_gene_mapping(source_gff: str, target_gff: str) -> GeneMapper:
    """
    Create a gene mapper between two annotation files.

    Args:
        source_gff: Path to source GFF (your reference)
        target_gff: Path to target GFF (e.g., SH1000)

    Returns:
        Configured GeneMapper with mappings built
    """
    mapper = GeneMapper()
    mapper.load_source_annotation(source_gff)
    mapper.load_target_annotation(target_gff)
    mapper.build_mappings()
    return mapper
