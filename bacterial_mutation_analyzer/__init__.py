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
"""

__version__ = "1.0.0"
__author__ = "DNA Sequencing Lab"

from .config import Config, PipelineConfig
from .pipeline.runner import PipelineRunner

__all__ = ["Config", "PipelineConfig", "PipelineRunner", "__version__"]
