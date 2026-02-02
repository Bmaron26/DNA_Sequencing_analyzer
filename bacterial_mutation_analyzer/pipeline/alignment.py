"""
Read alignment module for mapping reads to reference genome.

Supports multiple aligners:
- BWA-MEM (default): Best for Illumina short reads
- Minimap2: Fast aligner, good for long reads
- Bowtie2: Alternative short read aligner
"""

import os
import subprocess
import tempfile
from pathlib import Path
from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional, Tuple
import logging
import re

from ..config import AlignmentConfig

logger = logging.getLogger(__name__)


@dataclass
class AlignmentStats:
    """Alignment statistics."""
    total_reads: int = 0
    mapped_reads: int = 0
    properly_paired: int = 0
    singletons: int = 0
    duplicates: int = 0
    mapping_rate: float = 0.0
    mean_coverage: float = 0.0
    median_coverage: float = 0.0
    coverage_uniformity: float = 0.0
    mean_insert_size: float = 0.0
    insert_size_sd: float = 0.0
    mean_mapping_quality: float = 0.0
    reference_length: int = 0
    covered_bases: int = 0
    coverage_breadth: float = 0.0
    warnings: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            'total_reads': self.total_reads,
            'mapped_reads': self.mapped_reads,
            'properly_paired': self.properly_paired,
            'singletons': self.singletons,
            'duplicates': self.duplicates,
            'mapping_rate': round(self.mapping_rate, 4),
            'mean_coverage': round(self.mean_coverage, 2),
            'median_coverage': round(self.median_coverage, 2),
            'coverage_uniformity': round(self.coverage_uniformity, 4),
            'mean_insert_size': round(self.mean_insert_size, 2),
            'insert_size_sd': round(self.insert_size_sd, 2),
            'mean_mapping_quality': round(self.mean_mapping_quality, 2),
            'reference_length': self.reference_length,
            'covered_bases': self.covered_bases,
            'coverage_breadth': round(self.coverage_breadth, 4),
            'warnings': self.warnings,
        }


class Aligner:
    """
    Read aligner for mapping FASTQ reads to reference genome.
    """

    def __init__(self, config: Optional[AlignmentConfig] = None, output_dir: str = "alignment"):
        self.config = config or AlignmentConfig()
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self._check_tools()

    def _check_tools(self) -> Dict[str, bool]:
        """Check availability of alignment tools."""
        tools = {}
        for tool in ['bwa', 'minimap2', 'bowtie2', 'samtools']:
            try:
                result = subprocess.run(
                    [tool, '--version'] if tool != 'bwa' else [tool],
                    capture_output=True
                )
                tools[tool] = True
            except FileNotFoundError:
                tools[tool] = False
                if tool == self.config.aligner:
                    logger.warning(f"Configured aligner '{tool}' not found")
        return tools

    def index_reference(self, reference: str) -> str:
        """
        Index reference genome for alignment.

        Args:
            reference: Path to reference FASTA file

        Returns:
            Path to indexed reference
        """
        ref_path = Path(reference)

        # Check if index already exists
        if self.config.aligner == 'bwa':
            index_suffix = '.bwt'
        elif self.config.aligner == 'minimap2':
            index_suffix = '.mmi'
        elif self.config.aligner == 'bowtie2':
            index_suffix = '.1.bt2'
        else:
            index_suffix = '.bwt'

        index_path = ref_path.with_suffix(ref_path.suffix + index_suffix)

        if index_path.exists():
            logger.info(f"Using existing index: {index_path}")
            return str(ref_path)

        logger.info(f"Indexing reference genome: {reference}")

        if self.config.aligner == 'bwa':
            cmd = ['bwa', 'index', reference]
        elif self.config.aligner == 'minimap2':
            cmd = ['minimap2', '-d', str(index_path), reference]
        elif self.config.aligner == 'bowtie2':
            cmd = ['bowtie2-build', reference, str(ref_path)]
        else:
            cmd = ['bwa', 'index', reference]

        result = subprocess.run(cmd, capture_output=True)

        if result.returncode != 0:
            stderr_msg = result.stderr.decode('utf-8', errors='replace') if result.stderr else ''
            raise RuntimeError(f"Failed to index reference: {stderr_msg}")

        # Also create samtools faidx
        subprocess.run(['samtools', 'faidx', reference], capture_output=True)

        return str(ref_path)

    def align(self, fastq_files: List[str], reference: str,
             prefix: str = "sample") -> Tuple[str, AlignmentStats]:
        """
        Align reads to reference genome.

        Args:
            fastq_files: List of FASTQ files (1 or 2 for paired-end)
            reference: Path to reference FASTA file
            prefix: Output file prefix

        Returns:
            Tuple of (BAM file path, alignment statistics)
        """
        # Index reference if needed
        indexed_ref = self.index_reference(reference)

        is_paired = len(fastq_files) == 2
        sam_file = self.output_dir / f"{prefix}.sam"
        bam_file = self.output_dir / f"{prefix}.bam"
        sorted_bam = self.output_dir / f"{prefix}.sorted.bam"

        # Run alignment
        logger.info(f"Aligning reads using {self.config.aligner}")

        if self.config.aligner == 'bwa':
            self._align_bwa(fastq_files, indexed_ref, sam_file, is_paired)
        elif self.config.aligner == 'minimap2':
            self._align_minimap2(fastq_files, indexed_ref, sam_file, is_paired)
        elif self.config.aligner == 'bowtie2':
            self._align_bowtie2(fastq_files, indexed_ref, sam_file, is_paired)
        else:
            self._align_bwa(fastq_files, indexed_ref, sam_file, is_paired)

        # Convert SAM to BAM
        logger.info("Converting SAM to BAM")
        self._sam_to_bam(sam_file, bam_file)

        # Sort BAM
        logger.info("Sorting BAM file")
        self._sort_bam(bam_file, sorted_bam)

        # Index BAM
        logger.info("Indexing BAM file")
        subprocess.run(['samtools', 'index', str(sorted_bam)], check=True)

        # Mark duplicates if configured
        if self.config.mark_duplicates:
            logger.info("Marking duplicates")
            dedup_bam = self.output_dir / f"{prefix}.dedup.bam"
            self._mark_duplicates(sorted_bam, dedup_bam)
            sorted_bam = dedup_bam

        # Calculate statistics
        stats = self._calculate_stats(sorted_bam, reference)

        # Cleanup intermediate files
        for f in [sam_file, bam_file]:
            if f.exists():
                f.unlink()

        return str(sorted_bam), stats

    def _align_bwa(self, fastq_files: List[str], reference: str,
                   output: Path, is_paired: bool) -> None:
        """Align using BWA-MEM."""
        cmd = [
            'bwa', 'mem',
            '-t', str(self.config.threads),
            '-M',  # Mark shorter split hits as secondary
            reference,
        ] + fastq_files

        with open(output, 'w') as f:
            result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE)

        if result.returncode != 0:
            raise RuntimeError(f"BWA alignment failed: {result.stderr.decode('utf-8', errors='replace')}")

    def _align_minimap2(self, fastq_files: List[str], reference: str,
                        output: Path, is_paired: bool) -> None:
        """Align using minimap2."""
        cmd = [
            'minimap2',
            '-ax', 'sr',  # Short read preset
            '-t', str(self.config.threads),
            reference,
        ] + fastq_files

        with open(output, 'w') as f:
            result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE)

        if result.returncode != 0:
            raise RuntimeError(f"Minimap2 alignment failed: {result.stderr.decode('utf-8', errors='replace')}")

    def _align_bowtie2(self, fastq_files: List[str], reference: str,
                       output: Path, is_paired: bool) -> None:
        """Align using Bowtie2."""
        ref_base = str(Path(reference))

        if is_paired:
            cmd = [
                'bowtie2',
                '-p', str(self.config.threads),
                '-x', ref_base,
                '-1', fastq_files[0],
                '-2', fastq_files[1],
                '-S', str(output),
            ]
        else:
            cmd = [
                'bowtie2',
                '-p', str(self.config.threads),
                '-x', ref_base,
                '-U', fastq_files[0],
                '-S', str(output),
            ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            raise RuntimeError(f"Bowtie2 alignment failed: {result.stderr}")

    def _sam_to_bam(self, sam_file: Path, bam_file: Path) -> None:
        """Convert SAM to BAM."""
        cmd = ['samtools', 'view', '-bS', '-@', str(self.config.threads), str(sam_file)]

        with open(bam_file, 'wb') as f:
            result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE)

        if result.returncode != 0:
            raise RuntimeError(f"SAM to BAM conversion failed: {result.stderr.decode()}")

    def _sort_bam(self, input_bam: Path, output_bam: Path) -> None:
        """Sort BAM file."""
        cmd = [
            'samtools', 'sort',
            '-@', str(self.config.threads),
            '-o', str(output_bam),
            str(input_bam)
        ]

        result = subprocess.run(cmd, capture_output=True)

        if result.returncode != 0:
            raise RuntimeError(f"BAM sorting failed: {result.stderr.decode()}")

    def _mark_duplicates(self, input_bam: Path, output_bam: Path) -> None:
        """Mark duplicate reads using samtools."""
        # Use samtools markdup
        fixmate_bam = input_bam.with_suffix('.fixmate.bam')

        # First fix mate information
        cmd1 = [
            'samtools', 'fixmate',
            '-m', '-@', str(self.config.threads),
            str(input_bam), str(fixmate_bam)
        ]
        subprocess.run(cmd1, check=True, capture_output=True)

        # Sort by position
        sorted_fixmate = fixmate_bam.with_suffix('.sorted.bam')
        cmd2 = [
            'samtools', 'sort',
            '-@', str(self.config.threads),
            '-o', str(sorted_fixmate),
            str(fixmate_bam)
        ]
        subprocess.run(cmd2, check=True, capture_output=True)

        # Mark duplicates
        cmd3 = [
            'samtools', 'markdup',
            '-@', str(self.config.threads),
        ]

        if self.config.remove_duplicates:
            cmd3.append('-r')

        cmd3.extend([str(sorted_fixmate), str(output_bam)])

        subprocess.run(cmd3, check=True, capture_output=True)

        # Index the output
        subprocess.run(['samtools', 'index', str(output_bam)], check=True)

        # Cleanup
        fixmate_bam.unlink()
        sorted_fixmate.unlink()

    def _calculate_stats(self, bam_file: Path, reference: str) -> AlignmentStats:
        """Calculate alignment statistics from BAM file."""
        stats = AlignmentStats()

        # Get flagstat
        result = subprocess.run(
            ['samtools', 'flagstat', str(bam_file)],
            capture_output=True, text=True
        )

        if result.returncode == 0:
            stats = self._parse_flagstat(result.stdout, stats)

        # Get coverage statistics
        result = subprocess.run(
            ['samtools', 'depth', '-a', str(bam_file)],
            capture_output=True, text=True
        )

        if result.returncode == 0:
            stats = self._parse_depth(result.stdout, stats)

        # Get insert size for paired-end
        result = subprocess.run(
            ['samtools', 'stats', str(bam_file)],
            capture_output=True, text=True
        )

        if result.returncode == 0:
            stats = self._parse_samtools_stats(result.stdout, stats)

        # Add warnings
        if stats.mapping_rate < 0.9:
            stats.warnings.append(f"Low mapping rate: {stats.mapping_rate:.1%}")
        if stats.mean_coverage < 30:
            stats.warnings.append(f"Low coverage: {stats.mean_coverage:.1f}x")
        if stats.coverage_breadth < 0.95:
            stats.warnings.append(f"Incomplete coverage: {stats.coverage_breadth:.1%} of genome covered")

        return stats

    def _parse_flagstat(self, flagstat_output: str, stats: AlignmentStats) -> AlignmentStats:
        """Parse samtools flagstat output."""
        for line in flagstat_output.split('\n'):
            if 'in total' in line:
                match = re.search(r'(\d+)', line)
                if match:
                    stats.total_reads = int(match.group(1))
            elif 'mapped (' in line and 'primary' not in line:
                match = re.search(r'(\d+)', line)
                if match:
                    stats.mapped_reads = int(match.group(1))
            elif 'properly paired' in line:
                match = re.search(r'(\d+)', line)
                if match:
                    stats.properly_paired = int(match.group(1))
            elif 'singletons' in line:
                match = re.search(r'(\d+)', line)
                if match:
                    stats.singletons = int(match.group(1))
            elif 'duplicates' in line:
                match = re.search(r'(\d+)', line)
                if match:
                    stats.duplicates = int(match.group(1))

        if stats.total_reads > 0:
            stats.mapping_rate = stats.mapped_reads / stats.total_reads

        return stats

    def _parse_depth(self, depth_output: str, stats: AlignmentStats) -> AlignmentStats:
        """Parse samtools depth output."""
        depths = []
        covered = 0

        for line in depth_output.split('\n'):
            if not line.strip():
                continue
            parts = line.split('\t')
            if len(parts) >= 3:
                depth = int(parts[2])
                depths.append(depth)
                if depth > 0:
                    covered += 1

        if depths:
            stats.mean_coverage = sum(depths) / len(depths)
            sorted_depths = sorted(depths)
            stats.median_coverage = sorted_depths[len(sorted_depths) // 2]
            stats.reference_length = len(depths)
            stats.covered_bases = covered
            stats.coverage_breadth = covered / len(depths)

            # Calculate uniformity (1 - coefficient of variation)
            if stats.mean_coverage > 0:
                import math
                variance = sum((d - stats.mean_coverage) ** 2 for d in depths) / len(depths)
                std_dev = math.sqrt(variance)
                cv = std_dev / stats.mean_coverage
                stats.coverage_uniformity = max(0, 1 - cv)

        return stats

    def _parse_samtools_stats(self, stats_output: str, stats: AlignmentStats) -> AlignmentStats:
        """Parse samtools stats output."""
        for line in stats_output.split('\n'):
            if line.startswith('SN\tinsert size average:'):
                match = re.search(r'([\d.]+)', line.split(':')[1])
                if match:
                    stats.mean_insert_size = float(match.group(1))
            elif line.startswith('SN\tinsert size standard deviation:'):
                match = re.search(r'([\d.]+)', line.split(':')[1])
                if match:
                    stats.insert_size_sd = float(match.group(1))
            elif line.startswith('SN\taverage quality:'):
                match = re.search(r'([\d.]+)', line.split(':')[1])
                if match:
                    stats.mean_mapping_quality = float(match.group(1))

        return stats

    def get_coverage_data(self, bam_file: str, window_size: int = 1000) -> Dict[str, List[float]]:
        """
        Get coverage data in windows for plotting.

        Args:
            bam_file: Path to BAM file
            window_size: Window size for coverage calculation

        Returns:
            Dictionary with chromosome names as keys and coverage lists as values
        """
        coverage_data = {}

        result = subprocess.run(
            ['samtools', 'depth', '-a', bam_file],
            capture_output=True, text=True
        )

        if result.returncode != 0:
            return coverage_data

        current_chrom = None
        current_depths = []
        window_depths = []

        for line in result.stdout.split('\n'):
            if not line.strip():
                continue

            parts = line.split('\t')
            if len(parts) < 3:
                continue

            chrom = parts[0]
            depth = int(parts[2])

            if chrom != current_chrom:
                if current_chrom is not None and window_depths:
                    current_depths.append(sum(window_depths) / len(window_depths))
                    coverage_data[current_chrom] = current_depths

                current_chrom = chrom
                current_depths = []
                window_depths = []

            window_depths.append(depth)

            if len(window_depths) >= window_size:
                current_depths.append(sum(window_depths) / len(window_depths))
                window_depths = []

        # Don't forget the last chromosome
        if current_chrom is not None and (window_depths or current_depths):
            if window_depths:
                current_depths.append(sum(window_depths) / len(window_depths))
            coverage_data[current_chrom] = current_depths

        return coverage_data
