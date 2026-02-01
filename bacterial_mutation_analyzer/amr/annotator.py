"""
AMR annotator that combines database lookups with species-specific analysis.
"""

import json
from pathlib import Path
from typing import List, Dict, Any, Optional
from dataclasses import dataclass, field
import logging

from .database import AMRDatabase, AMRHit
from .species_config import SpeciesConfig, get_species_config

logger = logging.getLogger(__name__)


@dataclass
class AMRReport:
    """Complete AMR analysis report."""
    sample_id: str
    species: str
    total_variants: int = 0
    known_resistance_mutations: int = 0
    mutations_in_resistance_genes: int = 0
    amp_related_mutations: int = 0

    # Detailed hits
    known_hits: List[Dict[str, Any]] = field(default_factory=list)
    resistance_gene_hits: List[Dict[str, Any]] = field(default_factory=list)
    amp_related_hits: List[Dict[str, Any]] = field(default_factory=list)

    # Summaries
    by_drug_class: Dict[str, int] = field(default_factory=dict)
    by_mechanism: Dict[str, int] = field(default_factory=dict)
    affected_genes: List[str] = field(default_factory=list)

    # Interpretation
    resistance_profile: Dict[str, str] = field(default_factory=dict)
    warnings: List[str] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return {
            'sample_id': self.sample_id,
            'species': self.species,
            'summary': {
                'total_variants': self.total_variants,
                'known_resistance_mutations': self.known_resistance_mutations,
                'mutations_in_resistance_genes': self.mutations_in_resistance_genes,
                'amp_related_mutations': self.amp_related_mutations,
            },
            'known_hits': self.known_hits,
            'resistance_gene_hits': self.resistance_gene_hits,
            'amp_related_hits': self.amp_related_hits,
            'by_drug_class': self.by_drug_class,
            'by_mechanism': self.by_mechanism,
            'affected_genes': self.affected_genes,
            'resistance_profile': self.resistance_profile,
            'warnings': self.warnings,
            'notes': self.notes,
        }

    def save(self, filepath: str) -> None:
        """Save report to JSON file."""
        with open(filepath, 'w') as f:
            json.dump(self.to_dict(), f, indent=2)


class AMRAnnotator:
    """
    Annotates variants with AMR information.

    Combines:
    - Known mutation database (CARD, ResFinder)
    - Species-specific resistance gene information
    - AMP-specific analysis
    """

    def __init__(self, species: str = "", cache_dir: Optional[str] = None):
        self.database = AMRDatabase(cache_dir)
        self.species_config: Optional[SpeciesConfig] = None

        if species:
            self.set_species(species)

        # Load built-in data
        self.database.load_builtin_data()

    def set_species(self, species: str) -> bool:
        """
        Set species for analysis.

        Args:
            species: Species name

        Returns:
            True if species config found, False otherwise
        """
        self.species_config = get_species_config(species)

        if self.species_config:
            logger.info(f"Loaded species config for: {self.species_config.species}")
            return True
        else:
            logger.warning(f"No species config found for: {species}")
            return False

    def load_custom_species(self, config_path: str) -> None:
        """Load custom species configuration from file."""
        self.species_config = SpeciesConfig.load(config_path)
        logger.info(f"Loaded custom species config: {self.species_config.species}")

    def annotate_variants(self, variants: List[Dict[str, Any]],
                         sample_id: str = "sample") -> AMRReport:
        """
        Annotate variants with AMR information.

        Args:
            variants: List of variant dictionaries
            sample_id: Sample identifier

        Returns:
            AMRReport with complete analysis
        """
        species_name = self.species_config.species if self.species_config else "Unknown"

        report = AMRReport(
            sample_id=sample_id,
            species=species_name,
            total_variants=len(variants),
        )

        affected_genes = set()
        from collections import defaultdict
        drug_class_counts = defaultdict(int)
        mechanism_counts = defaultdict(int)

        for variant in variants:
            hit = self.database.annotate_variant(variant, self.species_config)

            gene = hit.gene
            if gene:
                affected_genes.add(gene)

            # Categorize hits
            if hit.is_known_resistance:
                report.known_resistance_mutations += 1
                report.known_hits.append(hit.to_dict())

                if hit.resistance_class:
                    drug_class_counts[hit.resistance_class] += 1
                    # Add to resistance profile
                    report.resistance_profile[hit.resistance_class] = "resistant"

            if hit.is_resistance_gene:
                report.mutations_in_resistance_genes += 1
                if not hit.is_known_resistance:
                    report.resistance_gene_hits.append(hit.to_dict())

            if hit.mechanism:
                mechanism_counts[hit.mechanism] += 1

            # Check AMP-specific
            if self._is_amp_related(hit):
                report.amp_related_mutations += 1
                if hit.to_dict() not in report.amp_related_hits:
                    report.amp_related_hits.append(hit.to_dict())

        report.affected_genes = sorted(affected_genes)
        report.by_drug_class = dict(drug_class_counts)
        report.by_mechanism = dict(mechanism_counts)

        # Generate interpretation
        self._generate_interpretation(report)

        return report

    def _is_amp_related(self, hit: AMRHit) -> bool:
        """Check if mutation is AMP-related."""
        # Check resistance class
        if hit.resistance_class and 'amp' in hit.resistance_class.lower():
            return True

        # Check gene against species config
        if self.species_config:
            gene = hit.gene.lower() if hit.gene else ''

            amp_genes = [g.lower() for g in self.species_config.amp_resistance_genes]
            if gene in amp_genes:
                return True

            membrane_genes = [g.lower() for g in self.species_config.membrane_genes]
            if gene in membrane_genes:
                return True

        return False

    def _generate_interpretation(self, report: AMRReport) -> None:
        """Generate interpretation notes and warnings."""
        # Known resistance mutations
        if report.known_resistance_mutations > 0:
            report.notes.append(
                f"Found {report.known_resistance_mutations} known resistance mutation(s)"
            )

            for drug_class, count in report.by_drug_class.items():
                if drug_class != 'AMP':
                    report.notes.append(
                        f"  - {drug_class}: {count} mutation(s)"
                    )

        # AMP resistance
        if report.amp_related_mutations > 0:
            report.notes.append(
                f"Found {report.amp_related_mutations} mutation(s) potentially affecting AMP susceptibility"
            )

            # Highlight specific genes
            amp_genes = set()
            for hit in report.amp_related_hits:
                if hit.get('gene'):
                    amp_genes.add(hit['gene'])

            if amp_genes:
                report.notes.append(
                    f"  - Affected genes: {', '.join(sorted(amp_genes))}"
                )

        # Novel mutations in resistance genes
        if report.mutations_in_resistance_genes > report.known_resistance_mutations:
            novel_count = report.mutations_in_resistance_genes - report.known_resistance_mutations
            report.warnings.append(
                f"Found {novel_count} novel mutation(s) in known resistance genes - "
                f"may represent emerging resistance"
            )

        # Membrane modifications (important for AMP)
        membrane_count = report.by_mechanism.get('membrane modification', 0)
        if membrane_count > 0:
            report.notes.append(
                f"Found {membrane_count} mutation(s) affecting membrane - "
                f"relevant for cationic AMP resistance"
            )

        # LPS modifications (for Gram-negatives)
        lps_count = report.by_mechanism.get('LPS modification', 0)
        if lps_count > 0:
            report.notes.append(
                f"Found {lps_count} mutation(s) affecting LPS modification - "
                f"relevant for polymyxin/colistin and AMP resistance"
            )

        # Check for convergent evolution hints
        if len(report.known_hits) > 3:
            report.warnings.append(
                "Multiple resistance mutations detected - verify if sample is clonal or mixed"
            )

    def annotate_comparison_results(self, comparison_results: Dict[str, Any],
                                   sample_results: Dict[str, Dict]) -> Dict[str, Any]:
        """
        Annotate multi-sample comparison results with AMR information.

        Args:
            comparison_results: Output from MultiSampleComparison
            sample_results: Dictionary of sample_id -> pipeline results

        Returns:
            Enriched comparison results with AMR annotations
        """
        enriched = comparison_results.copy()

        # Annotate convergent mutations
        if 'convergent_mutations' in enriched:
            for mut in enriched['convergent_mutations']:
                gene = mut.get('gene_name', '')
                aa_change = mut.get('amino_acid_change', '')

                # Check if known resistance
                if gene and aa_change:
                    normalized = self.database._normalize_aa_change(aa_change)
                    known = self.database.query_mutation(gene, normalized)

                    if known:
                        mut['amr_annotation'] = {
                            'is_known_resistance': True,
                            'drug_class': known.drug_class,
                            'antibiotics': known.antibiotics,
                            'evidence': known.evidence_level,
                        }
                    elif self.database.is_resistance_gene(gene):
                        mut['amr_annotation'] = {
                            'is_known_resistance': False,
                            'is_resistance_gene': True,
                            'drug_class': self.database.get_gene_drug_class(gene),
                        }

                # Check AMP relevance
                if self.species_config:
                    if gene.lower() in [g.lower() for g in self.species_config.amp_resistance_genes]:
                        mut['amp_relevant'] = True
                        mut['amp_notes'] = "Gene associated with AMP resistance"

        # Add AMR summary for each group
        if 'group_summaries' in enriched:
            for group_name, group_data in enriched['group_summaries'].items():
                amr_summary = {
                    'resistance_genes_mutated': [],
                    'known_resistance_mutations': [],
                    'amp_related_genes': [],
                }

                for gene in group_data.get('top_mutated_genes', {}).keys():
                    if self.database.is_resistance_gene(gene):
                        amr_summary['resistance_genes_mutated'].append(gene)

                    if self.species_config:
                        if gene.lower() in [g.lower() for g in self.species_config.amp_resistance_genes]:
                            amr_summary['amp_related_genes'].append(gene)

                group_data['amr_summary'] = amr_summary

        return enriched

    def generate_amr_report(self, variants: List[Dict[str, Any]],
                           sample_id: str,
                           output_path: str) -> str:
        """
        Generate and save AMR report.

        Args:
            variants: Variant list
            sample_id: Sample ID
            output_path: Output file path

        Returns:
            Path to generated report
        """
        report = self.annotate_variants(variants, sample_id)
        report.save(output_path)
        logger.info(f"AMR report saved to {output_path}")
        return output_path


def annotate_with_amr(result_path: str, species: str,
                     output_path: Optional[str] = None) -> AMRReport:
    """
    Convenience function to annotate pipeline result with AMR info.

    Args:
        result_path: Path to pipeline result JSON
        species: Species name
        output_path: Optional output path for report

    Returns:
        AMRReport
    """
    with open(result_path, 'r') as f:
        result = json.load(f)

    annotator = AMRAnnotator(species)
    report = annotator.annotate_variants(
        result.get('variants', []),
        result.get('sample_name', 'sample')
    )

    if output_path:
        report.save(output_path)

    return report
