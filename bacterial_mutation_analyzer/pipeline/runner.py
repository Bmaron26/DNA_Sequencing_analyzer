"""
Pipeline runner that orchestrates the complete WGS analysis workflow.
"""

import os
import json
import logging
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Any, Optional, Tuple
from dataclasses import dataclass, field, asdict

from ..config import PipelineConfig
from .qc import QualityControl, QCStats
from .alignment import Aligner, AlignmentStats
from .variant_calling import VariantCaller, VariantStats
from .annotation import VariantAnnotator, AnnotatedVariant

logger = logging.getLogger(__name__)


@dataclass
class PipelineResult:
    """Complete results from pipeline execution."""
    sample_name: str
    start_time: str
    end_time: str
    status: str = "completed"

    # Input files
    input_fastq: List[str] = field(default_factory=list)
    reference: str = ""
    annotation: str = ""

    # Output files
    filtered_fastq: List[str] = field(default_factory=list)
    bam_file: str = ""
    vcf_file: str = ""
    output_dir: str = ""

    # Statistics
    qc_stats: Optional[Dict[str, Any]] = None
    alignment_stats: Optional[Dict[str, Any]] = None
    variant_stats: Optional[Dict[str, Any]] = None
    annotation_summary: Optional[Dict[str, Any]] = None

    # Annotated variants
    variants: List[Dict[str, Any]] = field(default_factory=list)

    # Warnings and errors
    warnings: List[str] = field(default_factory=list)
    errors: List[str] = field(default_factory=list)

    # Outliers detected
    outliers: List[Dict[str, Any]] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return asdict(self)

    def save(self, output_path: str) -> None:
        """Save results to JSON file."""
        with open(output_path, 'w') as f:
            json.dump(self.to_dict(), f, indent=2, default=str)


class PipelineRunner:
    """
    Orchestrates the complete bacterial WGS mutation analysis pipeline.

    Steps:
    1. Input validation
    2. Quality control and preprocessing
    3. Read alignment to reference
    4. Variant calling
    5. Variant annotation
    6. Outlier detection
    7. Report generation
    """

    def __init__(self, config: Optional[PipelineConfig] = None):
        self.config = config or PipelineConfig()
        self.logger = logging.getLogger(__name__)

    def run(self, fastq_files: List[str], reference: str,
           annotation: Optional[str] = None,
           sample_name: str = "sample",
           output_dir: Optional[str] = None) -> PipelineResult:
        """
        Execute the complete analysis pipeline.

        Args:
            fastq_files: List of FASTQ files (1 for single-end, 2 for paired-end)
            reference: Path to reference genome FASTA
            annotation: Path to annotation file (GFF/GBK)
            sample_name: Sample identifier for output naming
            output_dir: Output directory (default: results/{sample_name})

        Returns:
            PipelineResult with all statistics and outputs
        """
        start_time = datetime.now()

        # Setup output directory
        if output_dir is None:
            output_dir = os.path.join(self.config.output.output_dir, sample_name)
        os.makedirs(output_dir, exist_ok=True)

        # Initialize result
        result = PipelineResult(
            sample_name=sample_name,
            start_time=start_time.isoformat(),
            end_time="",
            input_fastq=fastq_files,
            reference=reference,
            annotation=annotation or "",
            output_dir=output_dir,
        )

        try:
            # Step 1: Input validation
            self.logger.info("Step 1: Validating inputs")
            self._validate_inputs(fastq_files, reference, annotation, result)

            # Step 2: Quality control
            self.logger.info("Step 2: Quality control")
            filtered_fastq, qc_stats = self._run_qc(fastq_files, sample_name, output_dir)
            result.filtered_fastq = filtered_fastq
            result.qc_stats = qc_stats.to_dict()
            result.warnings.extend(qc_stats.warnings)

            # Step 3: Alignment
            self.logger.info("Step 3: Aligning reads")
            bam_file, align_stats = self._run_alignment(
                filtered_fastq, reference, sample_name, output_dir
            )
            result.bam_file = bam_file
            result.alignment_stats = align_stats.to_dict()
            result.warnings.extend(align_stats.warnings)

            # Step 4: Variant calling
            self.logger.info("Step 4: Calling variants")
            vcf_file, var_stats = self._run_variant_calling(
                bam_file, reference, sample_name, output_dir
            )
            result.vcf_file = vcf_file
            result.variant_stats = var_stats.to_dict()
            result.warnings.extend(var_stats.warnings)

            # Step 5: Annotation
            if annotation:
                self.logger.info("Step 5: Annotating variants")
                annotated_variants, ann_summary = self._run_annotation(
                    vcf_file, annotation, reference, output_dir
                )
                result.variants = [v.to_dict() for v in annotated_variants]
                result.annotation_summary = ann_summary
            else:
                self.logger.info("Step 5: Skipping annotation (no annotation file provided)")
                # Still parse variants without annotation
                caller = VariantCaller(output_dir=output_dir)
                result.variants = caller.parse_variants(vcf_file)

            # Step 6: Outlier detection
            if self.config.outlier.enabled:
                self.logger.info("Step 6: Detecting outliers")
                outliers = self._detect_outliers(result)
                result.outliers = outliers

            result.status = "completed"

        except Exception as e:
            self.logger.error(f"Pipeline failed: {str(e)}")
            result.status = "failed"
            result.errors.append(str(e))
            raise

        finally:
            result.end_time = datetime.now().isoformat()

            # Save result summary
            result_path = os.path.join(output_dir, f"{sample_name}_pipeline_result.json")
            result.save(result_path)
            self.logger.info(f"Pipeline results saved to {result_path}")

        return result

    def _validate_inputs(self, fastq_files: List[str], reference: str,
                        annotation: Optional[str], result: PipelineResult) -> None:
        """Validate input files."""
        from ..utils.validators import InputValidator

        validator = InputValidator()
        validation = validator.validate_all(fastq_files, reference, annotation)

        if not validation.is_valid:
            raise ValueError(f"Input validation failed: {validation.message}")

        result.warnings.extend(validation.warnings)

    def _run_qc(self, fastq_files: List[str], sample_name: str,
                output_dir: str) -> Tuple[List[str], QCStats]:
        """Run quality control."""
        qc_dir = os.path.join(output_dir, "qc")
        qc = QualityControl(config=self.config.qc, output_dir=qc_dir)

        filtered_fastq, stats = qc.process(fastq_files, prefix=sample_name)
        qc.generate_report(stats, sample_name)

        return filtered_fastq, stats

    def _run_alignment(self, fastq_files: List[str], reference: str,
                       sample_name: str, output_dir: str) -> Tuple[str, AlignmentStats]:
        """Run read alignment."""
        align_dir = os.path.join(output_dir, "alignment")
        aligner = Aligner(config=self.config.alignment, output_dir=align_dir)

        bam_file, stats = aligner.align(fastq_files, reference, prefix=sample_name)

        return bam_file, stats

    def _run_variant_calling(self, bam_file: str, reference: str,
                            sample_name: str, output_dir: str) -> Tuple[str, VariantStats]:
        """Run variant calling."""
        var_dir = os.path.join(output_dir, "variants")
        caller = VariantCaller(
            config=self.config.variant_calling,
            filter_config=self.config.filtering,
            output_dir=var_dir
        )

        vcf_file, stats = caller.call_variants(bam_file, reference, prefix=sample_name)

        return vcf_file, stats

    def _run_annotation(self, vcf_file: str, annotation_file: str,
                        reference: str, output_dir: str) -> Tuple[List[AnnotatedVariant], Dict]:
        """Run variant annotation."""
        ann_dir = os.path.join(output_dir, "annotation")
        annotator = VariantAnnotator(config=self.config.annotation, output_dir=ann_dir)

        annotated = annotator.annotate_variants(vcf_file, annotation_file, reference)
        summary = annotator.get_summary(annotated)

        return annotated, summary

    def _detect_outliers(self, result: PipelineResult) -> List[Dict[str, Any]]:
        """Detect potential outliers and anomalies."""
        outliers = []

        # Check QC outliers
        if result.qc_stats:
            qc = result.qc_stats
            if qc.get('q30_rate', 1) < 0.7:
                outliers.append({
                    'type': 'qc_quality',
                    'severity': 'high',
                    'message': f"Very low Q30 rate: {qc['q30_rate']:.1%}",
                    'value': qc['q30_rate']
                })

            retention = qc.get('read_retention_rate', 1)
            if retention < 0.5:
                outliers.append({
                    'type': 'qc_retention',
                    'severity': 'high',
                    'message': f"Very high read loss: {(1-retention):.1%} reads filtered",
                    'value': retention
                })

        # Check alignment outliers
        if result.alignment_stats:
            align = result.alignment_stats
            if align.get('mapping_rate', 1) < 0.8:
                outliers.append({
                    'type': 'alignment_rate',
                    'severity': 'high',
                    'message': f"Low mapping rate: {align['mapping_rate']:.1%}",
                    'value': align['mapping_rate']
                })

            if align.get('coverage_breadth', 1) < 0.9:
                outliers.append({
                    'type': 'coverage_breadth',
                    'severity': 'medium',
                    'message': f"Incomplete genome coverage: {align['coverage_breadth']:.1%}",
                    'value': align['coverage_breadth']
                })

        # Check variant outliers
        if result.variant_stats:
            var = result.variant_stats
            if var.get('ti_tv_ratio', 2) < 1.0:
                outliers.append({
                    'type': 'ti_tv_ratio',
                    'severity': 'medium',
                    'message': f"Unusual Ti/Tv ratio: {var['ti_tv_ratio']:.2f} (may indicate artifacts)",
                    'value': var['ti_tv_ratio']
                })

        # Check for mutation clusters (potential hypermutation)
        if result.variants:
            clusters = self._find_mutation_clusters(result.variants)
            for cluster in clusters:
                outliers.append({
                    'type': 'mutation_cluster',
                    'severity': 'medium',
                    'message': f"Mutation cluster: {cluster['count']} variants in {cluster['region']}",
                    'value': cluster
                })

        # Check for high-impact variants
        high_impact = [v for v in result.variants if v.get('effect_impact') == 'HIGH']
        if len(high_impact) > 10:
            outliers.append({
                'type': 'high_impact_count',
                'severity': 'info',
                'message': f"High number of high-impact variants: {len(high_impact)}",
                'value': len(high_impact)
            })

        return outliers

    def _find_mutation_clusters(self, variants: List[Dict],
                                window_size: int = 1000,
                                min_variants: int = 5) -> List[Dict]:
        """Find clusters of mutations that may indicate hypermutation."""
        clusters = []

        # Group variants by chromosome
        by_chrom = {}
        for v in variants:
            chrom = v.get('chromosome', '')
            if chrom not in by_chrom:
                by_chrom[chrom] = []
            by_chrom[chrom].append(v)

        for chrom, chrom_variants in by_chrom.items():
            # Sort by position
            chrom_variants.sort(key=lambda x: x.get('position', 0))

            # Sliding window to find clusters
            i = 0
            while i < len(chrom_variants):
                start_pos = chrom_variants[i].get('position', 0)
                j = i

                # Find all variants within window
                while j < len(chrom_variants):
                    if chrom_variants[j].get('position', 0) - start_pos > window_size:
                        break
                    j += 1

                # Check if this is a cluster
                count = j - i
                if count >= min_variants:
                    end_pos = chrom_variants[j-1].get('position', 0)
                    clusters.append({
                        'chromosome': chrom,
                        'start': start_pos,
                        'end': end_pos,
                        'count': count,
                        'region': f"{chrom}:{start_pos}-{end_pos}",
                        'density': count / max(end_pos - start_pos, 1) * 1000
                    })
                    i = j
                else:
                    i += 1

        return clusters


def run_pipeline(fastq_files: List[str], reference: str,
                annotation: Optional[str] = None,
                sample_name: str = "sample",
                output_dir: Optional[str] = None,
                config: Optional[PipelineConfig] = None) -> PipelineResult:
    """
    Convenience function to run the pipeline.

    Args:
        fastq_files: Input FASTQ file(s)
        reference: Reference genome FASTA
        annotation: Annotation file (GFF/GBK)
        sample_name: Sample name for outputs
        output_dir: Output directory
        config: Pipeline configuration

    Returns:
        PipelineResult with all outputs and statistics
    """
    runner = PipelineRunner(config)
    return runner.run(fastq_files, reference, annotation, sample_name, output_dir)
