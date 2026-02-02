"""
Multi-sample comparison module for analyzing mutations across samples and groups.

Supports:
- Within-group comparison (parallel evolution)
- Between-group comparison (treatment effects)
- Convergent evolution detection
- Population vs single colony interpretation
"""

import json
from pathlib import Path
from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional, Set, Tuple
from collections import defaultdict
import logging

from ..experiment import Experiment, SampleMetadata, SampleOrigin

logger = logging.getLogger(__name__)


@dataclass
class MutationOccurrence:
    """Tracks occurrence of a specific mutation across samples."""
    chromosome: str
    position: int
    reference: str
    alternative: str
    gene_name: str = ""
    locus_tag: str = ""
    effect: str = ""
    amino_acid_change: str = ""

    # Occurrence tracking
    samples: List[str] = field(default_factory=list)
    groups: List[str] = field(default_factory=list)
    frequencies: Dict[str, float] = field(default_factory=dict)  # sample_id -> freq

    @property
    def mutation_id(self) -> str:
        """Unique identifier for this mutation."""
        return f"{self.chromosome}:{self.position}:{self.reference}>{self.alternative}"

    @property
    def sample_count(self) -> int:
        return len(self.samples)

    @property
    def group_count(self) -> int:
        return len(set(self.groups))

    @property
    def is_convergent(self) -> bool:
        """Check if mutation appears in multiple independent samples."""
        return self.sample_count > 1

    @property
    def mean_frequency(self) -> float:
        """Mean allele frequency across samples."""
        if not self.frequencies:
            return 0.0
        return sum(self.frequencies.values()) / len(self.frequencies)

    def to_dict(self) -> Dict[str, Any]:
        return {
            'mutation_id': self.mutation_id,
            'chromosome': self.chromosome,
            'position': self.position,
            'reference': self.reference,
            'alternative': self.alternative,
            'gene_name': self.gene_name,
            'locus_tag': self.locus_tag,
            'effect': self.effect,
            'amino_acid_change': self.amino_acid_change,
            'sample_count': self.sample_count,
            'samples': self.samples,
            'groups': list(set(self.groups)),
            'group_count': self.group_count,
            'frequencies': self.frequencies,
            'mean_frequency': round(self.mean_frequency, 4),
            'is_convergent': self.is_convergent,
        }


@dataclass
class GroupComparison:
    """Comparison results between two groups."""
    group1: str
    group2: str
    shared_mutations: List[MutationOccurrence] = field(default_factory=list)
    unique_to_group1: List[MutationOccurrence] = field(default_factory=list)
    unique_to_group2: List[MutationOccurrence] = field(default_factory=list)
    shared_genes: Set[str] = field(default_factory=set)
    unique_genes_group1: Set[str] = field(default_factory=set)
    unique_genes_group2: Set[str] = field(default_factory=set)

    def to_dict(self) -> Dict[str, Any]:
        return {
            'group1': self.group1,
            'group2': self.group2,
            'shared_mutation_count': len(self.shared_mutations),
            'unique_to_group1_count': len(self.unique_to_group1),
            'unique_to_group2_count': len(self.unique_to_group2),
            'shared_mutations': [m.to_dict() for m in self.shared_mutations],
            'unique_to_group1': [m.to_dict() for m in self.unique_to_group1],
            'unique_to_group2': [m.to_dict() for m in self.unique_to_group2],
            'shared_genes': list(self.shared_genes),
            'unique_genes_group1': list(self.unique_genes_group1),
            'unique_genes_group2': list(self.unique_genes_group2),
        }


@dataclass
class ConvergentMutation:
    """A mutation showing parallel/convergent evolution."""
    mutation: MutationOccurrence
    convergence_score: float  # Higher = more convergent
    groups_affected: List[str] = field(default_factory=list)
    replicates_affected: int = 0
    interpretation: str = ""

    def to_dict(self) -> Dict[str, Any]:
        result = self.mutation.to_dict()
        result.update({
            'convergence_score': round(self.convergence_score, 3),
            'groups_affected': self.groups_affected,
            'replicates_affected': self.replicates_affected,
            'interpretation': self.interpretation,
        })
        return result


class MultiSampleComparison:
    """
    Compares mutations across multiple samples and treatment groups.
    """

    def __init__(self, experiment: Optional[Experiment] = None):
        self.experiment = experiment
        self.sample_results: Dict[str, Dict[str, Any]] = {}
        self.mutation_index: Dict[str, MutationOccurrence] = {}
        self.gene_mutations: Dict[str, List[MutationOccurrence]] = defaultdict(list)
        self.ancestral_mutations: Set[str] = set()  # Mutation IDs present in ancestral/WT
        self.ancestral_sample_id: Optional[str] = None

    def load_sample_result(self, sample_id: str, result_path: str,
                          group: Optional[str] = None) -> None:
        """
        Load a sample's analysis result.

        Args:
            sample_id: Sample identifier
            result_path: Path to pipeline result JSON
            group: Treatment group (optional if in experiment)
        """
        with open(result_path, 'r') as f:
            result = json.load(f)

        self.sample_results[sample_id] = result

        # Get group from experiment or parameter
        if self.experiment and sample_id in self.experiment.samples:
            group = self.experiment.samples[sample_id].group
        elif group is None:
            group = "unknown"

        # Index mutations
        for variant in result.get('variants', []):
            self._index_mutation(sample_id, group, variant)

        logger.info(f"Loaded {len(result.get('variants', []))} variants from {sample_id}")

    def load_all_samples(self, results_dir: str) -> int:
        """
        Load all sample results from experiment.

        Args:
            results_dir: Base directory containing sample results

        Returns:
            Number of samples loaded
        """
        if not self.experiment:
            raise ValueError("No experiment defined")

        loaded = 0
        for sample_id, sample in self.experiment.samples.items():
            result_path = Path(results_dir) / sample_id / f"{sample_id}_pipeline_result.json"

            if result_path.exists():
                self.load_sample_result(sample_id, str(result_path))
                loaded += 1
            else:
                logger.warning(f"Result not found for {sample_id}: {result_path}")

        return loaded

    def _index_mutation(self, sample_id: str, group: str, variant: Dict) -> None:
        """Index a mutation for cross-sample lookup."""
        mut_id = f"{variant.get('chromosome', '')}:{variant.get('position', '')}:" \
                 f"{variant.get('reference', '')}>{variant.get('alternative', '')}"

        if mut_id in self.mutation_index:
            # Add to existing mutation
            occ = self.mutation_index[mut_id]
            if sample_id not in occ.samples:
                occ.samples.append(sample_id)
                occ.groups.append(group)
            occ.frequencies[sample_id] = variant.get('allele_frequency', 0)
        else:
            # Create new mutation occurrence
            occ = MutationOccurrence(
                chromosome=variant.get('chromosome', ''),
                position=variant.get('position', 0),
                reference=variant.get('reference', ''),
                alternative=variant.get('alternative', ''),
                gene_name=variant.get('gene_name', ''),
                locus_tag=variant.get('locus_tag', ''),
                effect=variant.get('effect', ''),
                amino_acid_change=variant.get('amino_acid_change', ''),
                samples=[sample_id],
                groups=[group],
                frequencies={sample_id: variant.get('allele_frequency', 0)},
            )
            self.mutation_index[mut_id] = occ

            # Index by gene
            gene = variant.get('gene_name') or variant.get('locus_tag') or ''
            if gene:
                self.gene_mutations[gene].append(occ)

    def get_convergent_mutations(self, min_samples: int = 2,
                                within_group: bool = True) -> List[ConvergentMutation]:
        """
        Find mutations that appear in multiple samples (convergent evolution).

        Args:
            min_samples: Minimum number of samples for convergence
            within_group: If True, only count convergence within same group

        Returns:
            List of convergent mutations sorted by score
        """
        convergent = []

        for mut_id, occ in self.mutation_index.items():
            if within_group:
                # Count samples per group
                group_counts = defaultdict(int)
                for g in occ.groups:
                    group_counts[g] += 1

                # Check if convergent within any group
                max_within_group = max(group_counts.values()) if group_counts else 0

                if max_within_group >= min_samples:
                    # Calculate convergence score
                    score = self._calculate_convergence_score(occ, within_group=True)

                    convergent.append(ConvergentMutation(
                        mutation=occ,
                        convergence_score=score,
                        groups_affected=list(set(occ.groups)),
                        replicates_affected=max_within_group,
                        interpretation=self._interpret_convergence(occ, score),
                    ))
            else:
                if occ.sample_count >= min_samples:
                    score = self._calculate_convergence_score(occ, within_group=False)

                    convergent.append(ConvergentMutation(
                        mutation=occ,
                        convergence_score=score,
                        groups_affected=list(set(occ.groups)),
                        replicates_affected=occ.sample_count,
                        interpretation=self._interpret_convergence(occ, score),
                    ))

        # Sort by convergence score
        convergent.sort(key=lambda x: x.convergence_score, reverse=True)

        return convergent

    def _calculate_convergence_score(self, occ: MutationOccurrence,
                                     within_group: bool = True) -> float:
        """
        Calculate a convergence score for a mutation.

        Score considers:
        - Number of independent samples
        - Mean allele frequency
        - Effect impact
        - Whether in same or different groups
        """
        score = 0.0

        # Base score from sample count
        score += occ.sample_count * 2

        # Bonus for high frequency (likely under selection)
        if occ.mean_frequency > 0.5:
            score += 2
        elif occ.mean_frequency > 0.25:
            score += 1

        # Bonus for functional effects
        effect_lower = occ.effect.lower()
        if any(e in effect_lower for e in ['missense', 'nonsense', 'frameshift', 'stop']):
            score += 3
        elif 'synonymous' not in effect_lower and occ.effect:
            score += 1

        # Cross-group convergence is more significant
        if not within_group and occ.group_count > 1:
            score += occ.group_count * 2

        return score

    def _interpret_convergence(self, occ: MutationOccurrence, score: float) -> str:
        """Generate interpretation of convergent mutation."""
        parts = []

        if occ.sample_count >= 4:
            parts.append("Strong parallel evolution signal")
        elif occ.sample_count >= 2:
            parts.append("Parallel evolution detected")

        if occ.mean_frequency > 0.8:
            parts.append("near-fixation suggests strong selection")
        elif occ.mean_frequency > 0.5:
            parts.append("high frequency suggests positive selection")

        if occ.group_count > 1:
            parts.append(f"present in {occ.group_count} treatment groups")

        effect = occ.effect.lower()
        if 'missense' in effect:
            parts.append("amino acid change may affect protein function")
        elif 'frameshift' in effect or 'stop' in effect:
            parts.append("likely loss-of-function mutation")

        return "; ".join(parts) if parts else "Recurrent mutation"

    def compare_groups(self, group1: str, group2: str) -> GroupComparison:
        """
        Compare mutations between two treatment groups.

        Args:
            group1: First group name
            group2: Second group name

        Returns:
            GroupComparison with shared and unique mutations
        """
        group1_mutations = set()
        group2_mutations = set()
        group1_genes = set()
        group2_genes = set()

        for mut_id, occ in self.mutation_index.items():
            unique_groups = set(occ.groups)

            gene = occ.gene_name or occ.locus_tag

            if group1 in unique_groups:
                group1_mutations.add(mut_id)
                if gene:
                    group1_genes.add(gene)

            if group2 in unique_groups:
                group2_mutations.add(mut_id)
                if gene:
                    group2_genes.add(gene)

        shared_ids = group1_mutations & group2_mutations
        unique1_ids = group1_mutations - group2_mutations
        unique2_ids = group2_mutations - group1_mutations

        return GroupComparison(
            group1=group1,
            group2=group2,
            shared_mutations=[self.mutation_index[m] for m in shared_ids],
            unique_to_group1=[self.mutation_index[m] for m in unique1_ids],
            unique_to_group2=[self.mutation_index[m] for m in unique2_ids],
            shared_genes=group1_genes & group2_genes,
            unique_genes_group1=group1_genes - group2_genes,
            unique_genes_group2=group2_genes - group1_genes,
        )

    def get_group_summary(self, group_name: str) -> Dict[str, Any]:
        """
        Get mutation summary for a treatment group.

        Considers population samples may have mixed variants.
        """
        samples = []
        if self.experiment:
            samples = [s.sample_id for s in self.experiment.get_samples_by_group(group_name)]
        else:
            # Infer from loaded results
            for mut_id, occ in self.mutation_index.items():
                for i, g in enumerate(occ.groups):
                    if g == group_name and occ.samples[i] not in samples:
                        samples.append(occ.samples[i])

        if not samples:
            return {'group': group_name, 'error': 'No samples found'}

        # Collect mutations for this group
        group_mutations = []
        for mut_id, occ in self.mutation_index.items():
            if group_name in occ.groups:
                group_mutations.append(occ)

        # Count by category
        total_unique = len(group_mutations)
        shared_across_replicates = sum(
            1 for m in group_mutations
            if sum(1 for g in m.groups if g == group_name) > 1
        )

        # Gene-level analysis
        genes_mutated = set()
        gene_mutation_counts = defaultdict(int)
        for m in group_mutations:
            gene = m.gene_name or m.locus_tag
            if gene:
                genes_mutated.add(gene)
                gene_mutation_counts[gene] += 1

        # Effect distribution
        effect_counts = defaultdict(int)
        for m in group_mutations:
            effect_counts[m.effect] += 1

        # Population consideration
        population_samples = []
        if self.experiment:
            population_samples = [
                s.sample_id for s in self.experiment.get_samples_by_group(group_name)
                if s.is_population
            ]

        return {
            'group': group_name,
            'sample_count': len(samples),
            'samples': samples,
            'total_unique_mutations': total_unique,
            'shared_across_replicates': shared_across_replicates,
            'convergence_rate': round(shared_across_replicates / max(total_unique, 1), 3),
            'genes_mutated': len(genes_mutated),
            'top_mutated_genes': dict(sorted(
                gene_mutation_counts.items(),
                key=lambda x: x[1],
                reverse=True
            )[:10]),
            'effect_distribution': dict(effect_counts),
            'population_samples': population_samples,
            'population_note': (
                f"{len(population_samples)} samples are from pooled colonies; "
                "low-frequency variants may represent subpopulation mutations"
                if population_samples else None
            ),
        }

    def get_all_group_summaries(self) -> Dict[str, Dict[str, Any]]:
        """Get summaries for all groups."""
        groups = set()

        if self.experiment:
            groups = set(self.experiment.groups.keys())
        else:
            for occ in self.mutation_index.values():
                groups.update(occ.groups)

        return {g: self.get_group_summary(g) for g in groups}

    def find_treatment_specific_mutations(self, treatment_group: str,
                                         control_group: str = "Control") -> List[MutationOccurrence]:
        """
        Find mutations specific to a treatment group (not in control).

        Args:
            treatment_group: Treatment group name
            control_group: Control group name

        Returns:
            List of treatment-specific mutations
        """
        comparison = self.compare_groups(treatment_group, control_group)
        return comparison.unique_to_group1

    def get_gene_mutation_summary(self, gene_name: str) -> Dict[str, Any]:
        """Get summary of all mutations in a specific gene."""
        mutations = self.gene_mutations.get(gene_name, [])

        if not mutations:
            return {'gene': gene_name, 'mutation_count': 0}

        samples_affected = set()
        groups_affected = set()
        effects = defaultdict(int)
        aa_changes = []

        for m in mutations:
            samples_affected.update(m.samples)
            groups_affected.update(m.groups)
            effects[m.effect] += 1
            if m.amino_acid_change:
                aa_changes.append(m.amino_acid_change)

        return {
            'gene': gene_name,
            'mutation_count': len(mutations),
            'samples_affected': len(samples_affected),
            'sample_list': list(samples_affected),
            'groups_affected': list(groups_affected),
            'effects': dict(effects),
            'amino_acid_changes': aa_changes,
            'mutations': [m.to_dict() for m in mutations],
        }

    def generate_comparison_report(self, output_path: str) -> str:
        """
        Generate a comprehensive comparison report.

        Args:
            output_path: Output file path (JSON)

        Returns:
            Path to generated report
        """
        report = {
            'experiment': self.experiment.name if self.experiment else 'Unknown',
            'samples_analyzed': len(self.sample_results),
            'total_unique_mutations': len(self.mutation_index),
            'group_summaries': self.get_all_group_summaries(),
            'convergent_mutations': [
                cm.to_dict() for cm in self.get_convergent_mutations(min_samples=2)
            ],
        }

        # Add pairwise comparisons if multiple groups
        groups = list(report['group_summaries'].keys())
        if len(groups) >= 2:
            comparisons = []
            for i, g1 in enumerate(groups):
                for g2 in groups[i+1:]:
                    comp = self.compare_groups(g1, g2)
                    comparisons.append({
                        'groups': [g1, g2],
                        'shared_count': len(comp.shared_mutations),
                        'unique_to_first': len(comp.unique_to_group1),
                        'unique_to_second': len(comp.unique_to_group2),
                        'shared_genes': list(comp.shared_genes)[:20],  # Limit for readability
                    })
            report['pairwise_comparisons'] = comparisons

        with open(output_path, 'w') as f:
            json.dump(report, f, indent=2)

        logger.info(f"Comparison report saved to {output_path}")
        return output_path

    def export_mutation_matrix(self, output_path: str) -> str:
        """
        Export a sample x mutation presence matrix.

        Useful for downstream analysis (clustering, heatmaps).

        Args:
            output_path: Output CSV path

        Returns:
            Path to generated file
        """
        import csv

        samples = list(self.sample_results.keys())
        mutations = list(self.mutation_index.keys())

        with open(output_path, 'w', newline='') as f:
            writer = csv.writer(f)

            # Header
            header = ['mutation_id', 'chromosome', 'position', 'ref', 'alt',
                     'gene', 'effect'] + samples
            writer.writerow(header)

            # Data rows
            for mut_id in mutations:
                occ = self.mutation_index[mut_id]
                row = [
                    mut_id,
                    occ.chromosome,
                    occ.position,
                    occ.reference,
                    occ.alternative,
                    occ.gene_name or occ.locus_tag,
                    occ.effect,
                ]

                # Add presence (frequency) for each sample
                for sample in samples:
                    freq = occ.frequencies.get(sample, 0)
                    row.append(round(freq, 4) if freq > 0 else 0)

                writer.writerow(row)

        logger.info(f"Mutation matrix exported to {output_path}")
        return output_path

    def set_ancestral_sample(self, sample_id: str, result_path: Optional[str] = None) -> int:
        """
        Set a sample as the ancestral/wildtype reference for filtering.

        Mutations present in the ancestral sample will be excluded from
        evolved sample analysis (considered as pre-existing variants).

        Args:
            sample_id: Sample ID of the ancestral strain
            result_path: Path to result JSON (optional if already loaded)

        Returns:
            Number of ancestral mutations indexed
        """
        # Load if not already loaded
        if sample_id not in self.sample_results:
            if result_path:
                self.load_sample_result(sample_id, result_path, group="Ancestor")
            else:
                raise ValueError(f"Ancestral sample {sample_id} not loaded and no path provided")

        self.ancestral_sample_id = sample_id
        self.ancestral_mutations.clear()

        # Index all mutations from ancestral sample
        result = self.sample_results[sample_id]
        for variant in result.get('variants', []):
            mut_id = f"{variant.get('chromosome', '')}:{variant.get('position', '')}:" \
                     f"{variant.get('reference', '')}>{variant.get('alternative', '')}"
            self.ancestral_mutations.add(mut_id)

        logger.info(f"Set {sample_id} as ancestral sample with {len(self.ancestral_mutations)} variants")
        return len(self.ancestral_mutations)

    def get_novel_mutations(self, sample_id: str) -> List[Dict[str, Any]]:
        """
        Get mutations in a sample that are NOT present in the ancestral sample.

        These represent de novo mutations that arose during evolution.

        Args:
            sample_id: Sample to analyze

        Returns:
            List of novel variant dictionaries
        """
        if not self.ancestral_mutations:
            logger.warning("No ancestral sample set - returning all mutations")
            return self.sample_results.get(sample_id, {}).get('variants', [])

        if sample_id not in self.sample_results:
            raise ValueError(f"Sample {sample_id} not loaded")

        novel = []
        for variant in self.sample_results[sample_id].get('variants', []):
            mut_id = f"{variant.get('chromosome', '')}:{variant.get('position', '')}:" \
                     f"{variant.get('reference', '')}>{variant.get('alternative', '')}"

            if mut_id not in self.ancestral_mutations:
                novel.append(variant)

        return novel

    def get_filtered_mutation_index(self) -> Dict[str, MutationOccurrence]:
        """
        Get mutation index excluding ancestral mutations.

        Returns:
            Dictionary of mutation_id -> MutationOccurrence for non-ancestral mutations
        """
        if not self.ancestral_mutations:
            return self.mutation_index

        return {
            mut_id: occ for mut_id, occ in self.mutation_index.items()
            if mut_id not in self.ancestral_mutations
        }

    def export_novel_mutations(self, output_path: str,
                               exclude_ancestral_sample: bool = True) -> str:
        """
        Export only novel (non-ancestral) mutations to CSV.

        Args:
            output_path: Output CSV path
            exclude_ancestral_sample: Whether to exclude ancestral sample from output

        Returns:
            Path to generated file
        """
        import csv

        samples = [s for s in self.sample_results.keys()
                   if not exclude_ancestral_sample or s != self.ancestral_sample_id]

        filtered_mutations = self.get_filtered_mutation_index()

        with open(output_path, 'w', newline='') as f:
            writer = csv.writer(f)

            # Header
            header = ['mutation_id', 'chromosome', 'position', 'ref', 'alt',
                     'gene', 'effect', 'amino_acid_change', 'is_novel'] + samples
            writer.writerow(header)

            # Data rows - only non-ancestral mutations
            for mut_id, occ in filtered_mutations.items():
                row = [
                    mut_id,
                    occ.chromosome,
                    occ.position,
                    occ.reference,
                    occ.alternative,
                    occ.gene_name or occ.locus_tag,
                    occ.effect,
                    occ.amino_acid_change,
                    'yes',  # All are novel since we filtered
                ]

                # Add presence (frequency) for each sample
                for sample in samples:
                    freq = occ.frequencies.get(sample, 0)
                    row.append(round(freq, 4) if freq > 0 else 0)

                writer.writerow(row)

        logger.info(f"Novel mutations exported to {output_path} "
                    f"({len(filtered_mutations)} mutations, excluding {len(self.ancestral_mutations)} ancestral)")
        return output_path

    def generate_novel_mutations_report(self, output_path: str) -> str:
        """
        Generate a report focusing on novel (non-ancestral) mutations.

        Useful for identifying mutations that arose during evolution experiments.

        Args:
            output_path: Output JSON path

        Returns:
            Path to generated report
        """
        filtered_mutations = self.get_filtered_mutation_index()

        # Group by sample
        sample_novel_counts = {}
        for sample_id in self.sample_results.keys():
            if sample_id == self.ancestral_sample_id:
                continue
            novel = self.get_novel_mutations(sample_id)
            sample_novel_counts[sample_id] = {
                'total_mutations': len(self.sample_results[sample_id].get('variants', [])),
                'novel_mutations': len(novel),
                'ancestral_mutations': len(self.sample_results[sample_id].get('variants', [])) - len(novel),
            }

        # Find convergent novel mutations
        novel_convergent = []
        for mut_id, occ in filtered_mutations.items():
            # Exclude ancestral sample from count
            non_anc_samples = [s for s in occ.samples if s != self.ancestral_sample_id]
            if len(non_anc_samples) >= 2:
                score = self._calculate_convergence_score(occ, within_group=True)
                novel_convergent.append({
                    'mutation_id': mut_id,
                    'gene': occ.gene_name or occ.locus_tag,
                    'amino_acid_change': occ.amino_acid_change,
                    'effect': occ.effect,
                    'samples': non_anc_samples,
                    'sample_count': len(non_anc_samples),
                    'convergence_score': round(score, 2),
                    'interpretation': self._interpret_convergence(occ, score),
                })

        # Sort by convergence score
        novel_convergent.sort(key=lambda x: x['convergence_score'], reverse=True)

        report = {
            'ancestral_sample': self.ancestral_sample_id,
            'ancestral_variant_count': len(self.ancestral_mutations),
            'total_samples_analyzed': len(self.sample_results) - (1 if self.ancestral_sample_id else 0),
            'total_novel_mutations': len(filtered_mutations),
            'sample_summaries': sample_novel_counts,
            'convergent_novel_mutations': novel_convergent,
            'note': 'Novel mutations are those NOT present in the ancestral/WT sample. '
                    'These represent de novo mutations that arose during evolution.',
        }

        with open(output_path, 'w') as f:
            json.dump(report, f, indent=2)

        logger.info(f"Novel mutations report saved to {output_path}")
        return output_path
