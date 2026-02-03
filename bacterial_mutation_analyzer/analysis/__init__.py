"""
Analysis modules for mutation statistics, outlier detection, and multi-sample comparison.
"""

from .mutation_parser import MutationParser
from .outlier_detection import OutlierDetector
from .statistics import MutationStatistics
from .comparison import MultiSampleComparison, MutationOccurrence, GroupComparison
from .gene_mapper import GeneMapper, create_gene_mapping

__all__ = [
    "MutationParser",
    "OutlierDetector",
    "MutationStatistics",
    "MultiSampleComparison",
    "MutationOccurrence",
    "GroupComparison",
    "GeneMapper",
    "create_gene_mapping",
]
