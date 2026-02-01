"""
Mutation parser for extracting and formatting mutation data.
"""

import csv
import json
from pathlib import Path
from typing import List, Dict, Any, Optional
from dataclasses import dataclass, asdict
import logging

logger = logging.getLogger(__name__)


@dataclass
class Mutation:
    """Represents a single mutation."""
    chromosome: str
    position: int
    reference: str
    alternative: str
    variant_type: str
    quality: float
    depth: int
    allele_frequency: float
    gene_id: str = ""
    gene_name: str = ""
    locus_tag: str = ""
    product: str = ""
    effect: str = ""
    effect_impact: str = ""
    amino_acid_change: str = ""
    codon_change: str = ""
    location_type: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


class MutationParser:
    """
    Parser for mutation data from various sources.
    """

    def __init__(self):
        self.mutations: List[Mutation] = []

    def parse_from_result(self, result: Dict[str, Any]) -> List[Mutation]:
        """
        Parse mutations from pipeline result.

        Args:
            result: Pipeline result dictionary

        Returns:
            List of Mutation objects
        """
        self.mutations = []

        for v in result.get('variants', []):
            mutation = Mutation(
                chromosome=v.get('chromosome', ''),
                position=v.get('position', 0),
                reference=v.get('reference', ''),
                alternative=v.get('alternative', ''),
                variant_type=v.get('variant_type', ''),
                quality=v.get('quality', 0.0),
                depth=v.get('depth', 0),
                allele_frequency=v.get('allele_frequency', 0.0),
                gene_id=v.get('gene_id', ''),
                gene_name=v.get('gene_name', ''),
                locus_tag=v.get('locus_tag', ''),
                product=v.get('product', ''),
                effect=v.get('effect', ''),
                effect_impact=v.get('effect_impact', ''),
                amino_acid_change=v.get('amino_acid_change', ''),
                codon_change=v.get('codon_change', ''),
                location_type=v.get('location_type', ''),
            )
            self.mutations.append(mutation)

        return self.mutations

    def export_csv(self, output_path: str, mutations: Optional[List[Mutation]] = None) -> str:
        """
        Export mutations to CSV format.

        Args:
            output_path: Output file path
            mutations: List of mutations (uses self.mutations if None)

        Returns:
            Output file path
        """
        mutations = mutations or self.mutations

        if not mutations:
            logger.warning("No mutations to export")
            return output_path

        fieldnames = [
            'chromosome', 'position', 'reference', 'alternative',
            'variant_type', 'quality', 'depth', 'allele_frequency',
            'gene_name', 'locus_tag', 'product',
            'effect', 'effect_impact', 'amino_acid_change', 'codon_change',
            'location_type'
        ]

        with open(output_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
            writer.writeheader()
            for mutation in mutations:
                writer.writerow(mutation.to_dict())

        logger.info(f"Exported {len(mutations)} mutations to {output_path}")
        return output_path

    def export_tsv(self, output_path: str, mutations: Optional[List[Mutation]] = None) -> str:
        """
        Export mutations to TSV format.

        Args:
            output_path: Output file path
            mutations: List of mutations (uses self.mutations if None)

        Returns:
            Output file path
        """
        mutations = mutations or self.mutations

        if not mutations:
            logger.warning("No mutations to export")
            return output_path

        fieldnames = [
            'chromosome', 'position', 'reference', 'alternative',
            'variant_type', 'quality', 'depth', 'allele_frequency',
            'gene_name', 'locus_tag', 'product',
            'effect', 'effect_impact', 'amino_acid_change', 'codon_change',
            'location_type'
        ]

        with open(output_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames,
                                   delimiter='\t', extrasaction='ignore')
            writer.writeheader()
            for mutation in mutations:
                writer.writerow(mutation.to_dict())

        logger.info(f"Exported {len(mutations)} mutations to {output_path}")
        return output_path

    def export_json(self, output_path: str, mutations: Optional[List[Mutation]] = None) -> str:
        """
        Export mutations to JSON format.

        Args:
            output_path: Output file path
            mutations: List of mutations (uses self.mutations if None)

        Returns:
            Output file path
        """
        mutations = mutations or self.mutations

        data = [m.to_dict() for m in mutations]

        with open(output_path, 'w') as f:
            json.dump(data, f, indent=2)

        logger.info(f"Exported {len(mutations)} mutations to {output_path}")
        return output_path

    def filter_by_impact(self, impact: str,
                        mutations: Optional[List[Mutation]] = None) -> List[Mutation]:
        """Filter mutations by impact level."""
        mutations = mutations or self.mutations
        return [m for m in mutations if m.effect_impact == impact]

    def filter_by_type(self, variant_type: str,
                      mutations: Optional[List[Mutation]] = None) -> List[Mutation]:
        """Filter mutations by variant type."""
        mutations = mutations or self.mutations
        return [m for m in mutations if m.variant_type == variant_type]

    def filter_by_gene(self, gene_name: str,
                      mutations: Optional[List[Mutation]] = None) -> List[Mutation]:
        """Filter mutations by gene name or locus tag."""
        mutations = mutations or self.mutations
        gene_lower = gene_name.lower()
        return [m for m in mutations
                if gene_lower in (m.gene_name.lower(), m.locus_tag.lower(), m.gene_id.lower())]

    def filter_by_region(self, chromosome: str, start: int, end: int,
                        mutations: Optional[List[Mutation]] = None) -> List[Mutation]:
        """Filter mutations by genomic region."""
        mutations = mutations or self.mutations
        return [m for m in mutations
                if m.chromosome == chromosome and start <= m.position <= end]

    def filter_by_frequency(self, min_freq: float = 0.0, max_freq: float = 1.0,
                           mutations: Optional[List[Mutation]] = None) -> List[Mutation]:
        """Filter mutations by allele frequency."""
        mutations = mutations or self.mutations
        return [m for m in mutations
                if min_freq <= m.allele_frequency <= max_freq]

    def get_genes_with_mutations(self,
                                mutations: Optional[List[Mutation]] = None) -> Dict[str, int]:
        """Get count of mutations per gene."""
        mutations = mutations or self.mutations
        gene_counts = {}

        for m in mutations:
            gene = m.gene_name or m.locus_tag or m.gene_id
            if gene:
                gene_counts[gene] = gene_counts.get(gene, 0) + 1

        return dict(sorted(gene_counts.items(), key=lambda x: x[1], reverse=True))
