"""
Statistical analysis of mutations.
"""

from typing import List, Dict, Any, Optional, Tuple
from collections import defaultdict
from dataclasses import dataclass
import math
import logging

logger = logging.getLogger(__name__)


@dataclass
class MutationSpectrum:
    """6-class mutation spectrum."""
    c_to_a: int = 0
    c_to_g: int = 0
    c_to_t: int = 0
    t_to_a: int = 0
    t_to_c: int = 0
    t_to_g: int = 0

    def total(self) -> int:
        return sum([self.c_to_a, self.c_to_g, self.c_to_t,
                   self.t_to_a, self.t_to_c, self.t_to_g])

    def frequencies(self) -> Dict[str, float]:
        total = self.total()
        if total == 0:
            return {k: 0.0 for k in ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']}

        return {
            'C>A': self.c_to_a / total,
            'C>G': self.c_to_g / total,
            'C>T': self.c_to_t / total,
            'T>A': self.t_to_a / total,
            'T>C': self.t_to_c / total,
            'T>G': self.t_to_g / total,
        }

    def to_dict(self) -> Dict[str, int]:
        return {
            'C>A': self.c_to_a,
            'C>G': self.c_to_g,
            'C>T': self.c_to_t,
            'T>A': self.t_to_a,
            'T>C': self.t_to_c,
            'T>G': self.t_to_g,
        }


class MutationStatistics:
    """
    Calculate various statistics from mutation data.
    """

    # Complement bases
    COMPLEMENT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}

    def __init__(self, variants: List[Dict[str, Any]]):
        self.variants = variants
        self._snps = None
        self._indels = None

    @property
    def snps(self) -> List[Dict]:
        """Get SNP variants."""
        if self._snps is None:
            self._snps = [v for v in self.variants if v.get('variant_type') == 'SNP']
        return self._snps

    @property
    def indels(self) -> List[Dict]:
        """Get insertion/deletion variants."""
        if self._indels is None:
            self._indels = [v for v in self.variants
                          if v.get('variant_type') in ['insertion', 'deletion']]
        return self._indels

    def calculate_spectrum(self) -> MutationSpectrum:
        """
        Calculate 6-class mutation spectrum.

        Normalizes to pyrimidine reference (C or T).
        """
        spectrum = MutationSpectrum()

        for v in self.snps:
            ref = v.get('reference', '').upper()
            alt = v.get('alternative', '').upper()

            if len(ref) != 1 or len(alt) != 1:
                continue

            # Normalize to pyrimidine reference
            if ref in 'AG':
                ref = self.COMPLEMENT[ref]
                alt = self.COMPLEMENT[alt]

            if ref == 'C':
                if alt == 'A':
                    spectrum.c_to_a += 1
                elif alt == 'G':
                    spectrum.c_to_g += 1
                elif alt == 'T':
                    spectrum.c_to_t += 1
            elif ref == 'T':
                if alt == 'A':
                    spectrum.t_to_a += 1
                elif alt == 'C':
                    spectrum.t_to_c += 1
                elif alt == 'G':
                    spectrum.t_to_g += 1

        return spectrum

    def calculate_ti_tv_ratio(self) -> float:
        """
        Calculate transition/transversion ratio.

        Transitions: A<->G, C<->T
        Transversions: A<->C, A<->T, G<->C, G<->T
        """
        transitions = 0
        transversions = 0

        ti_pairs = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}

        for v in self.snps:
            ref = v.get('reference', '').upper()
            alt = v.get('alternative', '').upper()

            if len(ref) != 1 or len(alt) != 1:
                continue

            if (ref, alt) in ti_pairs:
                transitions += 1
            else:
                transversions += 1

        if transversions == 0:
            return float('inf') if transitions > 0 else 0.0

        return transitions / transversions

    def calculate_variant_density(self, genome_size: int) -> float:
        """
        Calculate variants per kilobase.

        Args:
            genome_size: Size of genome in base pairs

        Returns:
            Variants per kb
        """
        if genome_size == 0:
            return 0.0
        return len(self.variants) / (genome_size / 1000)

    def calculate_dnds_ratio(self) -> Optional[float]:
        """
        Calculate dN/dS ratio (nonsynonymous/synonymous).

        Requires annotated variants with effect information.
        """
        nonsynonymous = 0
        synonymous = 0

        for v in self.variants:
            effect = v.get('effect', '').lower()

            if 'synonymous' in effect and 'nonsynonymous' not in effect:
                synonymous += 1
            elif any(e in effect for e in ['missense', 'nonsense', 'frameshift',
                                           'stop_gained', 'stop_lost', 'start_lost']):
                nonsynonymous += 1

        if synonymous == 0:
            return None if nonsynonymous == 0 else float('inf')

        return nonsynonymous / synonymous

    def count_by_effect(self) -> Dict[str, int]:
        """Count variants by effect type."""
        counts = defaultdict(int)
        for v in self.variants:
            effect = v.get('effect', 'unknown')
            counts[effect] += 1
        return dict(sorted(counts.items(), key=lambda x: x[1], reverse=True))

    def count_by_impact(self) -> Dict[str, int]:
        """Count variants by impact level."""
        counts = defaultdict(int)
        for v in self.variants:
            impact = v.get('effect_impact', 'unknown')
            counts[impact] += 1
        return dict(counts)

    def count_by_gene(self) -> Dict[str, int]:
        """Count variants by gene."""
        counts = defaultdict(int)
        for v in self.variants:
            gene = v.get('gene_name') or v.get('locus_tag') or v.get('gene_id')
            if gene:
                counts[gene] += 1
        return dict(sorted(counts.items(), key=lambda x: x[1], reverse=True))

    def calculate_indel_size_distribution(self) -> Dict[str, Any]:
        """Calculate insertion and deletion size distributions."""
        insertions = []
        deletions = []

        for v in self.indels:
            ref = v.get('reference', '')
            alt = v.get('alternative', '')
            vtype = v.get('variant_type', '')

            if vtype == 'insertion':
                size = len(alt) - len(ref)
                insertions.append(size)
            elif vtype == 'deletion':
                size = len(ref) - len(alt)
                deletions.append(size)

        result = {
            'insertions': {
                'count': len(insertions),
                'sizes': insertions[:100],  # Limit for output
            },
            'deletions': {
                'count': len(deletions),
                'sizes': deletions[:100],
            }
        }

        if insertions:
            result['insertions']['mean_size'] = sum(insertions) / len(insertions)
            result['insertions']['max_size'] = max(insertions)

        if deletions:
            result['deletions']['mean_size'] = sum(deletions) / len(deletions)
            result['deletions']['max_size'] = max(deletions)

        return result

    def calculate_frequency_distribution(self) -> Dict[str, Any]:
        """Analyze allele frequency distribution."""
        frequencies = [v.get('allele_frequency', 0) for v in self.variants
                      if v.get('allele_frequency') is not None]

        if not frequencies:
            return {}

        # Bin frequencies
        bins = {
            '0-10%': 0,
            '10-25%': 0,
            '25-50%': 0,
            '50-75%': 0,
            '75-90%': 0,
            '90-100%': 0,
        }

        for f in frequencies:
            if f <= 0.1:
                bins['0-10%'] += 1
            elif f <= 0.25:
                bins['10-25%'] += 1
            elif f <= 0.5:
                bins['25-50%'] += 1
            elif f <= 0.75:
                bins['50-75%'] += 1
            elif f <= 0.9:
                bins['75-90%'] += 1
            else:
                bins['90-100%'] += 1

        return {
            'distribution': bins,
            'mean': sum(frequencies) / len(frequencies),
            'median': sorted(frequencies)[len(frequencies) // 2],
            'min': min(frequencies),
            'max': max(frequencies),
        }

    def get_summary(self) -> Dict[str, Any]:
        """Generate comprehensive summary statistics."""
        spectrum = self.calculate_spectrum()

        summary = {
            'total_variants': len(self.variants),
            'snps': len(self.snps),
            'indels': len(self.indels),
            'ti_tv_ratio': round(self.calculate_ti_tv_ratio(), 3),
            'mutation_spectrum': spectrum.to_dict(),
            'mutation_spectrum_frequencies': {k: round(v, 4)
                                              for k, v in spectrum.frequencies().items()},
            'by_effect': self.count_by_effect(),
            'by_impact': self.count_by_impact(),
            'genes_affected': len(self.count_by_gene()),
            'frequency_distribution': self.calculate_frequency_distribution(),
            'indel_sizes': self.calculate_indel_size_distribution(),
        }

        dnds = self.calculate_dnds_ratio()
        if dnds is not None:
            summary['dnds_ratio'] = round(dnds, 3) if dnds != float('inf') else 'inf'

        return summary

    def compare_samples(self, other: 'MutationStatistics') -> Dict[str, Any]:
        """
        Compare mutation statistics between two samples.

        Args:
            other: Another MutationStatistics instance

        Returns:
            Dictionary with comparison results
        """
        comparison = {
            'sample1': {
                'total_variants': len(self.variants),
                'snps': len(self.snps),
                'ti_tv': self.calculate_ti_tv_ratio(),
            },
            'sample2': {
                'total_variants': len(other.variants),
                'snps': len(other.snps),
                'ti_tv': other.calculate_ti_tv_ratio(),
            }
        }

        # Compare spectra
        spectrum1 = self.calculate_spectrum()
        spectrum2 = other.calculate_spectrum()

        freq1 = spectrum1.frequencies()
        freq2 = spectrum2.frequencies()

        comparison['spectrum_difference'] = {
            k: round(freq1[k] - freq2[k], 4) for k in freq1
        }

        # Find unique genes
        genes1 = set(self.count_by_gene().keys())
        genes2 = set(other.count_by_gene().keys())

        comparison['shared_genes'] = list(genes1 & genes2)
        comparison['unique_to_sample1'] = list(genes1 - genes2)
        comparison['unique_to_sample2'] = list(genes2 - genes1)

        return comparison
