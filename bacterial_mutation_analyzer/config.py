"""
Configuration management for the bacterial mutation analysis pipeline.
"""

import os
import yaml
from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any
import shutil


@dataclass
class QCConfig:
    """Quality control configuration."""
    min_quality: int = 20
    min_length: int = 50
    adapter_trimming: bool = True
    quality_encoding: str = "phred33"
    trim_front: int = 0
    trim_tail: int = 0
    max_n_ratio: float = 0.1
    complexity_filter: bool = True


@dataclass
class AlignmentConfig:
    """Alignment configuration."""
    aligner: str = "bwa"  # Options: bwa, minimap2, bowtie2
    threads: int = 4
    min_mapping_quality: int = 20
    mark_duplicates: bool = True
    remove_duplicates: bool = False
    max_insert_size: int = 1000


@dataclass
class VariantCallingConfig:
    """Variant calling configuration."""
    caller: str = "bcftools"  # Options: bcftools, freebayes, gatk
    min_base_quality: int = 20
    min_mapping_quality: int = 20
    min_depth: int = 10
    min_variant_frequency: float = 0.1
    min_strand_bias: float = 0.01
    ploidy: int = 1  # Bacteria are haploid
    call_indels: bool = True
    max_indel_length: int = 50


@dataclass
class AnnotationConfig:
    """Annotation configuration."""
    annotation_format: str = "auto"  # Options: auto, gff, gbk, gtf
    upstream_distance: int = 100
    downstream_distance: int = 0
    include_intergenic: bool = True
    protein_effect: bool = True


@dataclass
class FilterConfig:
    """Variant filtering configuration."""
    min_qual: float = 30.0
    min_depth: int = 10
    max_depth_factor: float = 3.0  # Max depth = mean_depth * factor
    min_variant_frequency: float = 0.1
    min_strand_depth: int = 2
    exclude_regions: List[str] = field(default_factory=list)
    include_only_regions: List[str] = field(default_factory=list)


@dataclass
class OutputConfig:
    """Output configuration."""
    output_dir: str = "results"
    prefix: str = "sample"
    formats: List[str] = field(default_factory=lambda: ["csv", "tsv", "vcf"])
    generate_html_report: bool = True
    keep_intermediate: bool = False
    compress_output: bool = True


@dataclass
class VisualizationConfig:
    """Visualization configuration."""
    enabled: bool = True
    plot_formats: List[str] = field(default_factory=lambda: ["png", "html"])
    dpi: int = 150
    coverage_plot: bool = True
    mutation_spectrum: bool = True
    circos_plot: bool = True
    quality_plots: bool = True
    interactive: bool = True


@dataclass
class OutlierConfig:
    """Outlier detection configuration."""
    enabled: bool = True
    coverage_zscore_threshold: float = 3.0
    cluster_distance_threshold: int = 100
    hypermutation_threshold: int = 10
    quality_drop_threshold: float = 0.5


@dataclass
class PipelineConfig:
    """Complete pipeline configuration."""
    qc: QCConfig = field(default_factory=QCConfig)
    alignment: AlignmentConfig = field(default_factory=AlignmentConfig)
    variant_calling: VariantCallingConfig = field(default_factory=VariantCallingConfig)
    annotation: AnnotationConfig = field(default_factory=AnnotationConfig)
    filtering: FilterConfig = field(default_factory=FilterConfig)
    output: OutputConfig = field(default_factory=OutputConfig)
    visualization: VisualizationConfig = field(default_factory=VisualizationConfig)
    outlier: OutlierConfig = field(default_factory=OutlierConfig)

    @classmethod
    def from_yaml(cls, yaml_path: str) -> "PipelineConfig":
        """Load configuration from a YAML file."""
        with open(yaml_path, 'r') as f:
            data = yaml.safe_load(f)
        return cls.from_dict(data)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "PipelineConfig":
        """Create configuration from dictionary."""
        config = cls()

        if 'qc' in data:
            config.qc = QCConfig(**data['qc'])
        if 'alignment' in data:
            config.alignment = AlignmentConfig(**data['alignment'])
        if 'variant_calling' in data:
            config.variant_calling = VariantCallingConfig(**data['variant_calling'])
        if 'annotation' in data:
            config.annotation = AnnotationConfig(**data['annotation'])
        if 'filtering' in data:
            config.filtering = FilterConfig(**data['filtering'])
        if 'output' in data:
            config.output = OutputConfig(**data['output'])
        if 'visualization' in data:
            config.visualization = VisualizationConfig(**data['visualization'])
        if 'outlier' in data:
            config.outlier = OutlierConfig(**data['outlier'])

        return config

    def to_yaml(self, yaml_path: str) -> None:
        """Save configuration to a YAML file."""
        from dataclasses import asdict
        data = asdict(self)
        with open(yaml_path, 'w') as f:
            yaml.dump(data, f, default_flow_style=False, sort_keys=False)

    def to_dict(self) -> Dict[str, Any]:
        """Convert configuration to dictionary."""
        from dataclasses import asdict
        return asdict(self)


class Config:
    """Global configuration and tool detection."""

    REQUIRED_TOOLS = {
        'bwa': 'BWA aligner for read mapping',
        'samtools': 'SAMtools for BAM file manipulation',
        'bcftools': 'BCFtools for variant calling',
    }

    OPTIONAL_TOOLS = {
        'fastp': 'Fast all-in-one FASTQ preprocessor',
        'fastqc': 'Quality control tool for sequencing data',
        'minimap2': 'Fast aligner for long reads',
        'freebayes': 'Bayesian variant caller',
        'snpeff': 'Variant annotation tool',
        'multiqc': 'Aggregate QC reports',
        'picard': 'BAM manipulation tools',
    }

    def __init__(self):
        self.available_tools: Dict[str, str] = {}
        self.missing_tools: Dict[str, str] = {}
        self._detect_tools()

    def _detect_tools(self) -> None:
        """Detect available bioinformatics tools."""
        all_tools = {**self.REQUIRED_TOOLS, **self.OPTIONAL_TOOLS}

        for tool, description in all_tools.items():
            path = shutil.which(tool)
            if path:
                self.available_tools[tool] = path
            else:
                self.missing_tools[tool] = description

    def check_requirements(self) -> tuple[bool, List[str]]:
        """Check if all required tools are available."""
        missing = []
        for tool in self.REQUIRED_TOOLS:
            if tool not in self.available_tools:
                missing.append(f"{tool}: {self.REQUIRED_TOOLS[tool]}")
        return len(missing) == 0, missing

    def get_tool_path(self, tool: str) -> Optional[str]:
        """Get the path to a specific tool."""
        return self.available_tools.get(tool)

    def is_available(self, tool: str) -> bool:
        """Check if a tool is available."""
        return tool in self.available_tools

    def get_status_report(self) -> str:
        """Generate a status report of available tools."""
        lines = ["Tool Availability Report", "=" * 40]

        lines.append("\nRequired Tools:")
        for tool, desc in self.REQUIRED_TOOLS.items():
            status = "✓" if tool in self.available_tools else "✗"
            lines.append(f"  {status} {tool}: {desc}")

        lines.append("\nOptional Tools:")
        for tool, desc in self.OPTIONAL_TOOLS.items():
            status = "✓" if tool in self.available_tools else "✗"
            lines.append(f"  {status} {tool}: {desc}")

        return "\n".join(lines)


def create_default_config(output_path: str = "config.yaml") -> None:
    """Create a default configuration file."""
    config = PipelineConfig()
    config.to_yaml(output_path)
