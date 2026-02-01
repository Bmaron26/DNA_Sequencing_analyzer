"""
Pipeline modules for bacterial mutation analysis.

This package contains the core analysis steps:
- QC: Quality control and read preprocessing
- Alignment: Read alignment to reference genome
- Variant Calling: Identification of SNPs and indels
- Annotation: Functional annotation of variants
- Assembly: De novo assembly (optional)
"""

from .qc import QualityControl
from .alignment import Aligner
from .variant_calling import VariantCaller
from .annotation import VariantAnnotator
from .runner import PipelineRunner

__all__ = [
    "QualityControl",
    "Aligner",
    "VariantCaller",
    "VariantAnnotator",
    "PipelineRunner",
]
