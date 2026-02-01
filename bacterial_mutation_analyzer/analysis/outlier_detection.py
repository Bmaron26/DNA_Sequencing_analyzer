"""
Outlier detection module for identifying anomalies in mutation data.

Detects:
- Unusual mutation clustering (hypermutation)
- Coverage anomalies
- Quality metric outliers
- Unexpected variant patterns
"""

from typing import List, Dict, Any, Optional, Tuple
from dataclasses import dataclass, field
from collections import defaultdict
import logging
import math

logger = logging.getLogger(__name__)


@dataclass
class Outlier:
    """Represents a detected outlier or anomaly."""
    outlier_type: str
    severity: str  # 'high', 'medium', 'low', 'info'
    message: str
    value: Any
    expected_range: Optional[Tuple[float, float]] = None
    location: Optional[str] = None
    details: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            'type': self.outlier_type,
            'severity': self.severity,
            'message': self.message,
            'value': self.value,
            'expected_range': self.expected_range,
            'location': self.location,
            'details': self.details,
        }


class OutlierDetector:
    """
    Detects outliers and anomalies in WGS mutation analysis results.
    """

    def __init__(self,
                 coverage_zscore_threshold: float = 3.0,
                 cluster_window_size: int = 1000,
                 cluster_min_variants: int = 5,
                 hypermutation_threshold: int = 10,
                 quality_thresholds: Optional[Dict[str, float]] = None):
        """
        Initialize outlier detector.

        Args:
            coverage_zscore_threshold: Z-score threshold for coverage outliers
            cluster_window_size: Window size (bp) for mutation clustering
            cluster_min_variants: Minimum variants to consider a cluster
            hypermutation_threshold: Variants per kb to flag as hypermutation
            quality_thresholds: Custom quality metric thresholds
        """
        self.coverage_zscore_threshold = coverage_zscore_threshold
        self.cluster_window_size = cluster_window_size
        self.cluster_min_variants = cluster_min_variants
        self.hypermutation_threshold = hypermutation_threshold

        # Default quality thresholds
        self.quality_thresholds = quality_thresholds or {
            'min_q30_rate': 0.7,
            'min_mapping_rate': 0.8,
            'min_coverage_breadth': 0.9,
            'min_read_retention': 0.5,
            'min_ti_tv_ratio': 1.0,
            'max_duplication_rate': 0.5,
        }

    def detect_all(self, result: Dict[str, Any]) -> List[Outlier]:
        """
        Run all outlier detection methods.

        Args:
            result: Pipeline result dictionary

        Returns:
            List of detected outliers
        """
        outliers = []

        # QC outliers
        if result.get('qc_stats'):
            outliers.extend(self.detect_qc_outliers(result['qc_stats']))

        # Alignment outliers
        if result.get('alignment_stats'):
            outliers.extend(self.detect_alignment_outliers(result['alignment_stats']))

        # Variant outliers
        if result.get('variant_stats'):
            outliers.extend(self.detect_variant_outliers(result['variant_stats']))

        # Mutation clustering
        if result.get('variants'):
            outliers.extend(self.detect_mutation_clusters(result['variants']))

        # Annotation anomalies
        if result.get('annotation_summary'):
            outliers.extend(self.detect_annotation_anomalies(
                result['annotation_summary'], result.get('variants', [])
            ))

        # Sort by severity
        severity_order = {'high': 0, 'medium': 1, 'low': 2, 'info': 3}
        outliers.sort(key=lambda x: severity_order.get(x.severity, 4))

        return outliers

    def detect_qc_outliers(self, qc_stats: Dict[str, Any]) -> List[Outlier]:
        """Detect quality control outliers."""
        outliers = []

        # Check Q30 rate
        q30_rate = qc_stats.get('q30_rate', 1.0)
        if q30_rate < self.quality_thresholds['min_q30_rate']:
            severity = 'high' if q30_rate < 0.5 else 'medium'
            outliers.append(Outlier(
                outlier_type='qc_quality',
                severity=severity,
                message=f"Low Q30 rate: {q30_rate:.1%} (expected ≥{self.quality_thresholds['min_q30_rate']:.0%})",
                value=q30_rate,
                expected_range=(self.quality_thresholds['min_q30_rate'], 1.0),
            ))

        # Check read retention
        retention = qc_stats.get('read_retention_rate', 1.0)
        if retention < self.quality_thresholds['min_read_retention']:
            severity = 'high' if retention < 0.3 else 'medium'
            outliers.append(Outlier(
                outlier_type='qc_retention',
                severity=severity,
                message=f"High read loss: {(1-retention):.1%} of reads filtered out",
                value=retention,
                expected_range=(self.quality_thresholds['min_read_retention'], 1.0),
            ))

        # Check duplication rate
        dup_rate = qc_stats.get('duplication_rate', 0.0)
        if dup_rate > self.quality_thresholds['max_duplication_rate']:
            outliers.append(Outlier(
                outlier_type='qc_duplication',
                severity='medium',
                message=f"High duplication rate: {dup_rate:.1%}",
                value=dup_rate,
                expected_range=(0.0, self.quality_thresholds['max_duplication_rate']),
            ))

        # Check GC content (bacteria typically 25-75%)
        gc_content = qc_stats.get('gc_content', 0.5)
        if gc_content < 0.2 or gc_content > 0.8:
            outliers.append(Outlier(
                outlier_type='qc_gc_content',
                severity='low',
                message=f"Unusual GC content: {gc_content:.1%} (expected 25-75%)",
                value=gc_content,
                expected_range=(0.25, 0.75),
            ))

        return outliers

    def detect_alignment_outliers(self, align_stats: Dict[str, Any]) -> List[Outlier]:
        """Detect alignment outliers."""
        outliers = []

        # Check mapping rate
        mapping_rate = align_stats.get('mapping_rate', 1.0)
        if mapping_rate < self.quality_thresholds['min_mapping_rate']:
            severity = 'high' if mapping_rate < 0.5 else 'medium'
            outliers.append(Outlier(
                outlier_type='alignment_rate',
                severity=severity,
                message=f"Low mapping rate: {mapping_rate:.1%} - reads may not match reference",
                value=mapping_rate,
                expected_range=(self.quality_thresholds['min_mapping_rate'], 1.0),
                details={'possible_causes': [
                    'Wrong reference genome',
                    'Contamination',
                    'Highly divergent sample',
                    'Poor quality reads'
                ]}
            ))

        # Check coverage breadth
        coverage_breadth = align_stats.get('coverage_breadth', 1.0)
        if coverage_breadth < self.quality_thresholds['min_coverage_breadth']:
            severity = 'high' if coverage_breadth < 0.7 else 'medium'
            outliers.append(Outlier(
                outlier_type='coverage_breadth',
                severity=severity,
                message=f"Incomplete genome coverage: {coverage_breadth:.1%}",
                value=coverage_breadth,
                expected_range=(self.quality_thresholds['min_coverage_breadth'], 1.0),
            ))

        # Check coverage uniformity
        uniformity = align_stats.get('coverage_uniformity', 1.0)
        if uniformity < 0.5:
            outliers.append(Outlier(
                outlier_type='coverage_uniformity',
                severity='medium',
                message=f"Uneven coverage distribution (uniformity: {uniformity:.1%})",
                value=uniformity,
                expected_range=(0.7, 1.0),
            ))

        # Check mean coverage
        mean_coverage = align_stats.get('mean_coverage', 0)
        if mean_coverage < 10:
            outliers.append(Outlier(
                outlier_type='low_coverage',
                severity='high',
                message=f"Very low coverage: {mean_coverage:.1f}x (recommend ≥30x)",
                value=mean_coverage,
                expected_range=(30, None),
            ))
        elif mean_coverage < 30:
            outliers.append(Outlier(
                outlier_type='low_coverage',
                severity='medium',
                message=f"Low coverage: {mean_coverage:.1f}x (recommend ≥30x)",
                value=mean_coverage,
                expected_range=(30, None),
            ))

        return outliers

    def detect_variant_outliers(self, var_stats: Dict[str, Any]) -> List[Outlier]:
        """Detect variant calling outliers."""
        outliers = []

        # Check Ti/Tv ratio
        ti_tv = var_stats.get('ti_tv_ratio', 2.0)
        total_snps = var_stats.get('snps', 0)

        if total_snps >= 50:  # Only meaningful with enough SNPs
            if ti_tv < self.quality_thresholds['min_ti_tv_ratio']:
                outliers.append(Outlier(
                    outlier_type='ti_tv_ratio',
                    severity='medium',
                    message=f"Unusual Ti/Tv ratio: {ti_tv:.2f} (expected >1.5 for real variants)",
                    value=ti_tv,
                    expected_range=(1.5, 3.0),
                    details={
                        'interpretation': 'Low Ti/Tv may indicate sequencing artifacts or contamination'
                    }
                ))
            elif ti_tv > 4.0:
                outliers.append(Outlier(
                    outlier_type='ti_tv_ratio',
                    severity='low',
                    message=f"High Ti/Tv ratio: {ti_tv:.2f}",
                    value=ti_tv,
                    expected_range=(1.5, 3.0),
                    details={
                        'interpretation': 'High Ti/Tv may indicate selection bias or few mutations'
                    }
                ))

        # Check variant quality
        mean_qual = var_stats.get('mean_quality', 100)
        if mean_qual < 30:
            outliers.append(Outlier(
                outlier_type='variant_quality',
                severity='medium',
                message=f"Low mean variant quality: {mean_qual:.1f}",
                value=mean_qual,
                expected_range=(50, None),
            ))

        return outliers

    def detect_mutation_clusters(self, variants: List[Dict]) -> List[Outlier]:
        """
        Detect clusters of mutations that may indicate hypermutation
        or recombination events.
        """
        outliers = []

        # Group variants by chromosome
        by_chrom = defaultdict(list)
        for v in variants:
            chrom = v.get('chromosome', '')
            pos = v.get('position', 0)
            by_chrom[chrom].append({'pos': pos, 'variant': v})

        # Find clusters
        for chrom, chrom_variants in by_chrom.items():
            # Sort by position
            chrom_variants.sort(key=lambda x: x['pos'])
            positions = [v['pos'] for v in chrom_variants]

            if len(positions) < self.cluster_min_variants:
                continue

            # Sliding window approach
            i = 0
            while i < len(positions):
                start_pos = positions[i]
                j = i

                # Find all variants within window
                while j < len(positions) and positions[j] - start_pos <= self.cluster_window_size:
                    j += 1

                count = j - i
                if count >= self.cluster_min_variants:
                    end_pos = positions[j - 1]
                    region_size = max(end_pos - start_pos, 1)
                    density = count / (region_size / 1000)  # per kb

                    severity = 'medium'
                    if density >= self.hypermutation_threshold:
                        severity = 'high'

                    # Get affected genes
                    affected_genes = set()
                    for idx in range(i, j):
                        gene = (chrom_variants[idx]['variant'].get('gene_name') or
                               chrom_variants[idx]['variant'].get('locus_tag'))
                        if gene:
                            affected_genes.add(gene)

                    outliers.append(Outlier(
                        outlier_type='mutation_cluster',
                        severity=severity,
                        message=f"Mutation cluster: {count} variants in {region_size:,} bp region",
                        value=count,
                        location=f"{chrom}:{start_pos}-{end_pos}",
                        details={
                            'density_per_kb': round(density, 2),
                            'affected_genes': list(affected_genes),
                            'region_size': region_size,
                        }
                    ))

                    i = j  # Skip past this cluster
                else:
                    i += 1

        return outliers

    def detect_annotation_anomalies(self, annotation_summary: Dict,
                                   variants: List[Dict]) -> List[Outlier]:
        """Detect anomalies in annotation patterns."""
        outliers = []

        # Check for genes with many mutations
        genes_affected = annotation_summary.get('genes_affected', [])
        high_impact_count = annotation_summary.get('by_impact', {}).get('HIGH', 0)

        if high_impact_count > 20:
            outliers.append(Outlier(
                outlier_type='high_impact_count',
                severity='info',
                message=f"High number of high-impact variants: {high_impact_count}",
                value=high_impact_count,
            ))

        # Count mutations per gene
        gene_mutation_counts = defaultdict(int)
        for v in variants:
            gene = v.get('gene_id') or v.get('locus_tag') or v.get('gene_name')
            if gene:
                gene_mutation_counts[gene] += 1

        # Find genes with unusually high mutation counts
        if gene_mutation_counts:
            counts = list(gene_mutation_counts.values())
            mean_count = sum(counts) / len(counts)
            std_count = math.sqrt(sum((c - mean_count) ** 2 for c in counts) / len(counts))

            for gene, count in gene_mutation_counts.items():
                if std_count > 0:
                    zscore = (count - mean_count) / std_count
                    if zscore > 3:
                        outliers.append(Outlier(
                            outlier_type='gene_mutation_hotspot',
                            severity='info',
                            message=f"Gene {gene} has {count} mutations (z-score: {zscore:.1f})",
                            value=count,
                            location=gene,
                            details={'zscore': round(zscore, 2), 'mean': round(mean_count, 2)}
                        ))

        return outliers

    def detect_coverage_outliers(self, coverage_data: Dict[str, List[float]],
                                window_size: int = 1000) -> List[Outlier]:
        """
        Detect regions with unusual coverage.

        Args:
            coverage_data: Dictionary with chromosome names and coverage lists
            window_size: Size of coverage windows
        """
        outliers = []

        for chrom, coverages in coverage_data.items():
            if len(coverages) < 10:
                continue

            mean_cov = sum(coverages) / len(coverages)
            if mean_cov == 0:
                continue

            std_cov = math.sqrt(sum((c - mean_cov) ** 2 for c in coverages) / len(coverages))

            if std_cov == 0:
                continue

            # Find regions with unusual coverage
            for i, cov in enumerate(coverages):
                zscore = (cov - mean_cov) / std_cov

                if abs(zscore) > self.coverage_zscore_threshold:
                    start_pos = i * window_size
                    end_pos = start_pos + window_size

                    if zscore > 0:
                        outliers.append(Outlier(
                            outlier_type='high_coverage_region',
                            severity='low',
                            message=f"High coverage region: {cov:.0f}x (mean: {mean_cov:.0f}x)",
                            value=cov,
                            location=f"{chrom}:{start_pos}-{end_pos}",
                            details={
                                'zscore': round(zscore, 2),
                                'possible_causes': ['Duplication', 'Repetitive region', 'Amplification']
                            }
                        ))
                    else:
                        outliers.append(Outlier(
                            outlier_type='low_coverage_region',
                            severity='medium' if cov < 5 else 'low',
                            message=f"Low coverage region: {cov:.0f}x (mean: {mean_cov:.0f}x)",
                            value=cov,
                            location=f"{chrom}:{start_pos}-{end_pos}",
                            details={
                                'zscore': round(zscore, 2),
                                'possible_causes': ['Deletion', 'GC bias', 'Unmappable region']
                            }
                        ))

        return outliers
