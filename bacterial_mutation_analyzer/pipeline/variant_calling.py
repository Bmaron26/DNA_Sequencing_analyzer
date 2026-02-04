"""
Variant calling module for identifying mutations.

Supports multiple variant callers:
- bcftools (default): Fast and memory-efficient
- FreeBayes: Bayesian variant caller, good for low-frequency variants
- GATK: Industry standard, requires more setup
"""

import os
import subprocess
from pathlib import Path
from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional, Tuple
import logging

from ..config import VariantCallingConfig, FilterConfig
from ..utils.file_handlers import VcfParser, VcfVariant

logger = logging.getLogger(__name__)


@dataclass
class VariantStats:
    """Variant calling statistics."""
    total_variants: int = 0
    snps: int = 0
    insertions: int = 0
    deletions: int = 0
    complex_variants: int = 0
    mnps: int = 0
    filtered_variants: int = 0
    passed_variants: int = 0
    transitions: int = 0
    transversions: int = 0
    ti_tv_ratio: float = 0.0
    het_hom_ratio: float = 0.0
    mean_depth: float = 0.0
    mean_quality: float = 0.0
    warnings: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            'total_variants': self.total_variants,
            'snps': self.snps,
            'insertions': self.insertions,
            'deletions': self.deletions,
            'complex_variants': self.complex_variants,
            'mnps': self.mnps,
            'filtered_variants': self.filtered_variants,
            'passed_variants': self.passed_variants,
            'transitions': self.transitions,
            'transversions': self.transversions,
            'ti_tv_ratio': round(self.ti_tv_ratio, 3),
            'mean_depth': round(self.mean_depth, 2),
            'mean_quality': round(self.mean_quality, 2),
            'warnings': self.warnings,
        }


class VariantCaller:
    """
    Variant caller for identifying mutations from aligned reads.
    """

    # Transition pairs (purine-purine or pyrimidine-pyrimidine)
    TRANSITIONS = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}

    def __init__(self, config: Optional[VariantCallingConfig] = None,
                 filter_config: Optional[FilterConfig] = None,
                 output_dir: str = "variants"):
        self.config = config or VariantCallingConfig()
        self.filter_config = filter_config or FilterConfig()
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self._check_tools()

    def _check_tools(self) -> Dict[str, bool]:
        """Check availability of variant calling tools."""
        tools = {}
        for tool in ['bcftools', 'freebayes']:
            try:
                result = subprocess.run(
                    [tool, '--version'],
                    capture_output=True,
                    text=True
                )
                tools[tool] = True
            except FileNotFoundError:
                tools[tool] = False
        return tools

    def call_variants(self, bam_file: str, reference: str,
                     prefix: str = "sample",
                     min_alternate_fraction: float = 0.01,
                     min_alternate_count: int = 5,
                     min_qual_filter: Optional[float] = None) -> Tuple[str, VariantStats]:
        """
        Call variants from aligned BAM file.

        Args:
            bam_file: Path to sorted, indexed BAM file
            reference: Path to reference FASTA file
            prefix: Output file prefix
            min_alternate_fraction: For FreeBayes - min fraction of reads with alt (default: 0.01)
            min_alternate_count: For FreeBayes - min number of reads with alt (default: 5)
            min_qual_filter: Override minimum quality filter (default: use config)

        Returns:
            Tuple of (VCF file path, variant statistics)
        """
        # Add caller name to output for distinction
        caller_name = self.config.caller
        raw_vcf = self.output_dir / f"{prefix}.{caller_name}.raw.vcf"
        filtered_vcf = self.output_dir / f"{prefix}.{caller_name}.filtered.vcf"

        logger.info(f"Calling variants using {caller_name}")

        if self.config.caller == 'bcftools':
            self._call_bcftools(bam_file, reference, raw_vcf)
        elif self.config.caller == 'freebayes':
            self._call_freebayes(bam_file, reference, raw_vcf,
                               min_alternate_fraction=min_alternate_fraction,
                               min_alternate_count=min_alternate_count)
        else:
            self._call_bcftools(bam_file, reference, raw_vcf)

        # Filter variants
        logger.info("Filtering variants")
        filter_qual = min_qual_filter if min_qual_filter is not None else self.filter_config.min_qual
        self._filter_variants(raw_vcf, filtered_vcf, min_qual=filter_qual)

        # Calculate statistics
        stats = self._calculate_stats(filtered_vcf)

        # Compress and index
        self._compress_and_index(filtered_vcf)

        return str(filtered_vcf), stats

    def _call_bcftools(self, bam_file: str, reference: str, output_vcf: Path) -> None:
        """Call variants using bcftools."""
        # Generate pileup
        mpileup_cmd = [
            'bcftools', 'mpileup',
            '-f', reference,
            '-q', str(self.config.min_mapping_quality),
            '-Q', str(self.config.min_base_quality),
            '-d', '10000',  # Max depth
            '-a', 'FORMAT/AD,FORMAT/DP,FORMAT/SP,INFO/AD',  # Annotations
            bam_file
        ]

        # Call variants
        call_cmd = [
            'bcftools', 'call',
            '-m',  # Multiallelic caller
            '-v',  # Output variant sites only
            '--ploidy', str(self.config.ploidy),
            '-Ov',  # Output VCF format
            '-o', str(output_vcf)
        ]

        # Pipe mpileup to call
        logger.debug(f"Running: {' '.join(mpileup_cmd)} | {' '.join(call_cmd)}")

        mpileup_proc = subprocess.Popen(
            mpileup_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )

        call_proc = subprocess.Popen(
            call_cmd,
            stdin=mpileup_proc.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )

        mpileup_proc.stdout.close()
        stdout, stderr = call_proc.communicate()

        if call_proc.returncode != 0:
            raise RuntimeError(f"bcftools call failed: {stderr.decode()}")

    def _call_freebayes(self, bam_file: str, reference: str, output_vcf: Path,
                        min_alternate_fraction: float = 0.01,
                        min_alternate_count: int = 5,
                        pooled_discrete: bool = False,
                        pooled_continuous: bool = True) -> None:
        """
        Call variants using FreeBayes.

        Optimized for pooled/population samples (e.g., 10 colonies pooled).

        Args:
            bam_file: Path to BAM file
            reference: Path to reference FASTA
            output_vcf: Output VCF path
            min_alternate_fraction: Minimum fraction of reads supporting alt allele (default: 0.01 = 1%)
            min_alternate_count: Minimum number of reads supporting alt allele (default: 5)
            pooled_discrete: Use pooled-discrete mode (known number of samples)
            pooled_continuous: Use pooled-continuous mode (unknown mixture)
        """
        cmd = [
            'freebayes',
            '-f', reference,
            '-p', str(self.config.ploidy),
            '--min-base-quality', str(self.config.min_base_quality),
            '--min-mapping-quality', str(self.config.min_mapping_quality),
            '--min-coverage', str(self.config.min_depth),
            '--min-alternate-fraction', str(min_alternate_fraction),
            '--min-alternate-count', str(min_alternate_count),
            bam_file
        ]

        # Add pooled mode options for population samples
        if pooled_continuous:
            cmd.insert(-1, '--pooled-continuous')
        elif pooled_discrete:
            cmd.insert(-1, '--pooled-discrete')

        logger.info(f"FreeBayes command: {' '.join(cmd)}")

        with open(output_vcf, 'w') as f:
            result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE)

        if result.returncode != 0:
            stderr_msg = result.stderr.decode('utf-8', errors='replace') if result.stderr else ''
            raise RuntimeError(f"FreeBayes failed: {stderr_msg}")

    def _filter_variants(self, input_vcf: Path, output_vcf: Path,
                        min_qual: Optional[float] = None) -> None:
        """Filter variants based on quality criteria."""
        # Build filter expression
        filters = []

        qual_threshold = min_qual if min_qual is not None else self.filter_config.min_qual
        if qual_threshold > 0:
            filters.append(f'QUAL>={qual_threshold}')

        if self.filter_config.min_depth > 0:
            filters.append(f'INFO/DP>={self.filter_config.min_depth}')

        filter_expr = ' && '.join(filters) if filters else 'QUAL>=0'

        cmd = [
            'bcftools', 'filter',
            '-i', filter_expr,
            '-o', str(output_vcf),
            str(input_vcf)
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            # If bcftools filter fails, just copy the file
            logger.warning(f"bcftools filter failed, using unfiltered variants: {result.stderr}")
            import shutil
            shutil.copy(input_vcf, output_vcf)

    def _compress_and_index(self, vcf_file: Path) -> None:
        """Compress VCF with bgzip and index with tabix."""
        try:
            # Compress
            subprocess.run(['bgzip', '-f', str(vcf_file)], check=True, capture_output=True)
            # Index
            subprocess.run(['tabix', '-p', 'vcf', f'{vcf_file}.gz'], check=True, capture_output=True)
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            logger.warning(f"Could not compress/index VCF: {e}")

    def _calculate_stats(self, vcf_file: Path) -> VariantStats:
        """Calculate variant statistics from VCF file."""
        stats = VariantStats()

        parser = VcfParser(str(vcf_file))
        variants = parser.parse()

        depths = []
        qualities = []

        for variant in variants:
            stats.total_variants += 1

            # Count by type
            if len(variant.ref) == 1 and len(variant.alt) == 1:
                stats.snps += 1
                # Count transitions/transversions
                if (variant.ref.upper(), variant.alt.upper()) in self.TRANSITIONS:
                    stats.transitions += 1
                else:
                    stats.transversions += 1
            elif len(variant.ref) < len(variant.alt):
                stats.insertions += 1
            elif len(variant.ref) > len(variant.alt):
                stats.deletions += 1
            else:
                if len(variant.ref) > 1:
                    stats.mnps += 1
                else:
                    stats.complex_variants += 1

            # Track quality metrics
            if variant.qual > 0:
                qualities.append(variant.qual)
            if variant.depth > 0:
                depths.append(variant.depth)

            # Count filtered
            if variant.filter == 'PASS' or variant.filter == '.':
                stats.passed_variants += 1
            else:
                stats.filtered_variants += 1

        # Calculate summary statistics
        if depths:
            stats.mean_depth = sum(depths) / len(depths)
        if qualities:
            stats.mean_quality = sum(qualities) / len(qualities)
        if stats.transversions > 0:
            stats.ti_tv_ratio = stats.transitions / stats.transversions

        # Add warnings
        if stats.ti_tv_ratio < 1.5 and stats.snps > 100:
            stats.warnings.append(
                f"Low Ti/Tv ratio ({stats.ti_tv_ratio:.2f}), may indicate sequencing artifacts"
            )
        if stats.mean_quality < 50:
            stats.warnings.append(
                f"Low mean variant quality ({stats.mean_quality:.1f})"
            )

        return stats

    def parse_variants(self, vcf_file: str) -> List[Dict[str, Any]]:
        """
        Parse VCF file and return structured variant data.

        Args:
            vcf_file: Path to VCF file

        Returns:
            List of variant dictionaries
        """
        # Handle compressed VCF
        if vcf_file.endswith('.gz'):
            vcf_path = vcf_file
        elif os.path.exists(f"{vcf_file}.gz"):
            vcf_path = f"{vcf_file}.gz"
        else:
            vcf_path = vcf_file

        parser = VcfParser(vcf_path)
        variants = parser.parse()

        result = []
        for v in variants:
            variant_dict = {
                'chromosome': v.chrom,
                'position': v.pos,
                'reference': v.ref,
                'alternative': v.alt,
                'quality': v.qual,
                'filter': v.filter,
                'type': v.variant_type,
                'depth': v.depth,
                'allele_frequency': v.allele_frequency,
            }

            # Add INFO fields
            for key, value in v.info.items():
                if key not in variant_dict:
                    variant_dict[f'info_{key}'] = value

            result.append(variant_dict)

        return result
