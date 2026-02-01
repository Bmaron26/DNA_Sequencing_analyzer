"""
Bacterial Mutation Analyzer

A comprehensive pipeline for whole genome sequencing (WGS) analysis
of bacterial samples to identify mutations, particularly useful for
studying antimicrobial resistance evolution.

Main Features:
- Quality control of FASTQ reads
- Read alignment to reference genome
- Variant calling and filtering
- Mutation annotation using reference annotations
- Outlier and anomaly detection
- CSV/TSV output with mutation details
- Interactive visualizations
- Multi-sample comparison with treatment groups
- AMR database integration (CARD, ResFinder)
- Species-specific resistance gene analysis
"""

__version__ = "1.1.0"
__author__ = "DNA Sequencing Lab"

from .config import Config, PipelineConfig
from .pipeline.runner import PipelineRunner
from .experiment import Experiment, SampleMetadata, TreatmentGroup, SampleOrigin

__all__ = [
    "Config",
    "PipelineConfig",
    "PipelineRunner",
    "Experiment",
    "SampleMetadata",
    "TreatmentGroup",
    "SampleOrigin",
    "__version__",
]
