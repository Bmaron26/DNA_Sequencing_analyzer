"""
Quality control module for FASTQ preprocessing.

Supports multiple QC tools:
- fastp (preferred): Fast all-in-one FASTQ preprocessor
- Built-in: Python-based QC when external tools unavailable
"""

import os
import json
import subprocess
import tempfile
from pathlib import Path
from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional, Tuple
import logging

from ..config import QCConfig
from ..utils.file_handlers import FastqReader

logger = logging.getLogger(__name__)


@dataclass
class QCStats:
    """Quality control statistics."""
    total_reads: int = 0
    total_bases: int = 0
    reads_after_filter: int = 0
    bases_after_filter: int = 0
    q20_rate: float = 0.0
    q30_rate: float = 0.0
    gc_content: float = 0.0
    mean_length: float = 0.0
    mean_quality: float = 0.0
    adapter_trimmed_reads: int = 0
    low_quality_reads: int = 0
    n_content: float = 0.0
    duplication_rate: float = 0.0
    warnings: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            'total_reads': self.total_reads,
            'total_bases': self.total_bases,
            'reads_after_filter': self.reads_after_filter,
            'bases_after_filter': self.bases_after_filter,
            'q20_rate': round(self.q20_rate, 4),
            'q30_rate': round(self.q30_rate, 4),
            'gc_content': round(self.gc_content, 4),
            'mean_length': round(self.mean_length, 2),
            'mean_quality': round(self.mean_quality, 2),
            'adapter_trimmed_reads': self.adapter_trimmed_reads,
            'low_quality_reads': self.low_quality_reads,
            'n_content': round(self.n_content, 4),
            'duplication_rate': round(self.duplication_rate, 4),
            'read_retention_rate': round(self.reads_after_filter / max(self.total_reads, 1), 4),
            'warnings': self.warnings,
        }


class QualityControl:
    """
    Quality control processor for FASTQ files.

    Performs adapter trimming, quality filtering, and generates QC reports.
    """

    def __init__(self, config: Optional[QCConfig] = None, output_dir: str = "qc_output"):
        self.config = config or QCConfig()
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self._fastp_available = self._check_fastp()

    def _check_fastp(self) -> bool:
        """Check if fastp is available."""
        try:
            result = subprocess.run(
                ['fastp', '--version'],
                capture_output=True,
                text=True
            )
            return result.returncode == 0
        except FileNotFoundError:
            return False

    def process(self, fastq_files: List[str], prefix: str = "sample") -> Tuple[List[str], QCStats]:
        """
        Process FASTQ files for quality control.

        Args:
            fastq_files: List of input FASTQ file paths (1 for single-end, 2 for paired-end)
            prefix: Output file prefix

        Returns:
            Tuple of (output FASTQ files, QC statistics)
        """
        is_paired = len(fastq_files) == 2

        if self._fastp_available:
            logger.info("Using fastp for quality control")
            return self._process_with_fastp(fastq_files, prefix, is_paired)
        else:
            logger.info("fastp not available, using built-in QC")
            return self._process_builtin(fastq_files, prefix, is_paired)

    def _process_with_fastp(self, fastq_files: List[str], prefix: str,
                           is_paired: bool) -> Tuple[List[str], QCStats]:
        """Process using fastp."""
        output_files = []
        json_report = self.output_dir / f"{prefix}_fastp.json"
        html_report = self.output_dir / f"{prefix}_fastp.html"

        cmd = [
            'fastp',
            '--json', str(json_report),
            '--html', str(html_report),
            '--qualified_quality_phred', str(self.config.min_quality),
            '--length_required', str(self.config.min_length),
            '--thread', '4',
        ]

        if self.config.adapter_trimming:
            cmd.append('--detect_adapter_for_pe' if is_paired else '--detect_adapter')

        if self.config.trim_front > 0:
            cmd.extend(['--trim_front1', str(self.config.trim_front)])
        if self.config.trim_tail > 0:
            cmd.extend(['--trim_tail1', str(self.config.trim_tail)])

        if self.config.complexity_filter:
            cmd.append('--low_complexity_filter')

        if is_paired:
            out1 = self.output_dir / f"{prefix}_R1_filtered.fastq.gz"
            out2 = self.output_dir / f"{prefix}_R2_filtered.fastq.gz"
            cmd.extend([
                '-i', fastq_files[0],
                '-I', fastq_files[1],
                '-o', str(out1),
                '-O', str(out2),
            ])
            output_files = [str(out1), str(out2)]
        else:
            out1 = self.output_dir / f"{prefix}_filtered.fastq.gz"
            cmd.extend([
                '-i', fastq_files[0],
                '-o', str(out1),
            ])
            output_files = [str(out1)]

        logger.info(f"Running fastp: {' '.join(cmd)}")

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            logger.error(f"fastp failed: {result.stderr}")
            raise RuntimeError(f"fastp failed: {result.stderr}")

        # Parse JSON report
        stats = self._parse_fastp_json(json_report)

        return output_files, stats

    def _parse_fastp_json(self, json_path: Path) -> QCStats:
        """Parse fastp JSON report."""
        with open(json_path, 'r') as f:
            data = json.load(f)

        stats = QCStats()

        # Before filtering stats
        before = data.get('summary', {}).get('before_filtering', {})
        stats.total_reads = before.get('total_reads', 0)
        stats.total_bases = before.get('total_bases', 0)

        # After filtering stats
        after = data.get('summary', {}).get('after_filtering', {})
        stats.reads_after_filter = after.get('total_reads', 0)
        stats.bases_after_filter = after.get('total_bases', 0)
        stats.q20_rate = after.get('q20_rate', 0)
        stats.q30_rate = after.get('q30_rate', 0)
        stats.gc_content = after.get('gc_content', 0)

        # Read length
        read1_len = after.get('read1_mean_length', 0)
        read2_len = after.get('read2_mean_length', 0)
        if read2_len > 0:
            stats.mean_length = (read1_len + read2_len) / 2
        else:
            stats.mean_length = read1_len

        # Filtering stats
        filtering = data.get('filtering_result', {})
        stats.low_quality_reads = filtering.get('low_quality_reads', 0)
        stats.adapter_trimmed_reads = filtering.get('adapter_trimmed_reads', 0)

        # Duplication
        dup = data.get('duplication', {})
        stats.duplication_rate = dup.get('rate', 0)

        # Generate warnings
        if stats.q30_rate < 0.8:
            stats.warnings.append(f"Low Q30 rate: {stats.q30_rate:.1%}")
        if stats.reads_after_filter / max(stats.total_reads, 1) < 0.8:
            stats.warnings.append(
                f"High read loss during filtering: "
                f"{(1 - stats.reads_after_filter / max(stats.total_reads, 1)):.1%}"
            )
        if stats.duplication_rate > 0.3:
            stats.warnings.append(f"High duplication rate: {stats.duplication_rate:.1%}")

        return stats

    def _process_builtin(self, fastq_files: List[str], prefix: str,
                        is_paired: bool) -> Tuple[List[str], QCStats]:
        """
        Built-in Python QC when external tools unavailable.
        Note: This is slower but doesn't require external dependencies.
        """
        import gzip
        from collections import defaultdict

        stats = QCStats()
        output_files = []

        for i, input_file in enumerate(fastq_files):
            suffix = f"_R{i+1}" if is_paired else ""
            output_file = self.output_dir / f"{prefix}{suffix}_filtered.fastq.gz"
            output_files.append(str(output_file))

            reader = FastqReader(input_file)

            with gzip.open(output_file, 'wt') as out_f:
                for record in reader:
                    stats.total_reads += 1
                    stats.total_bases += record.length

                    # Quality filter
                    if record.mean_quality < self.config.min_quality:
                        stats.low_quality_reads += 1
                        continue

                    # Length filter
                    if record.length < self.config.min_length:
                        continue

                    # N content filter
                    n_ratio = record.sequence.upper().count('N') / record.length
                    if n_ratio > self.config.max_n_ratio:
                        continue

                    # Quality trimming from ends
                    seq, qual = self._trim_quality(
                        record.sequence, record.quality,
                        self.config.min_quality
                    )

                    if len(seq) < self.config.min_length:
                        continue

                    # Count Q20/Q30 bases
                    for q in qual:
                        score = ord(q) - 33
                        if score >= 20:
                            stats.q20_rate += 1
                        if score >= 30:
                            stats.q30_rate += 1

                    stats.reads_after_filter += 1
                    stats.bases_after_filter += len(seq)
                    stats.gc_content += sum(1 for b in seq.upper() if b in 'GC') / len(seq)

                    # Write filtered read
                    out_f.write(f"@{record.id}\n{seq}\n+\n{qual}\n")

        # Calculate final statistics
        if stats.reads_after_filter > 0:
            stats.mean_length = stats.bases_after_filter / stats.reads_after_filter
            stats.gc_content /= stats.reads_after_filter
            stats.q20_rate /= stats.bases_after_filter
            stats.q30_rate /= stats.bases_after_filter

        # Generate warnings
        if stats.q30_rate < 0.8:
            stats.warnings.append(f"Low Q30 rate: {stats.q30_rate:.1%}")
        retention = stats.reads_after_filter / max(stats.total_reads, 1)
        if retention < 0.8:
            stats.warnings.append(f"High read loss during filtering: {(1 - retention):.1%}")

        return output_files, stats

    def _trim_quality(self, sequence: str, quality: str,
                     min_qual: int) -> Tuple[str, str]:
        """Trim low-quality bases from read ends."""
        # Trim from 3' end
        end = len(sequence)
        while end > 0 and ord(quality[end - 1]) - 33 < min_qual:
            end -= 1

        # Trim from 5' end
        start = 0
        while start < end and ord(quality[start]) - 33 < min_qual:
            start += 1

        return sequence[start:end], quality[start:end]

    def generate_report(self, stats: QCStats, prefix: str) -> str:
        """Generate a text QC report."""
        report_path = self.output_dir / f"{prefix}_qc_report.txt"

        lines = [
            "=" * 60,
            "Quality Control Report",
            "=" * 60,
            "",
            f"Total reads:              {stats.total_reads:,}",
            f"Total bases:              {stats.total_bases:,}",
            f"Reads after filtering:    {stats.reads_after_filter:,}",
            f"Bases after filtering:    {stats.bases_after_filter:,}",
            f"Read retention rate:      {stats.reads_after_filter / max(stats.total_reads, 1):.1%}",
            "",
            f"Mean read length:         {stats.mean_length:.1f} bp",
            f"Mean quality score:       {stats.mean_quality:.1f}",
            f"Q20 rate:                 {stats.q20_rate:.1%}",
            f"Q30 rate:                 {stats.q30_rate:.1%}",
            f"GC content:               {stats.gc_content:.1%}",
            "",
            f"Adapter trimmed reads:    {stats.adapter_trimmed_reads:,}",
            f"Low quality reads:        {stats.low_quality_reads:,}",
            f"Duplication rate:         {stats.duplication_rate:.1%}",
            "",
        ]

        if stats.warnings:
            lines.append("Warnings:")
            for warning in stats.warnings:
                lines.append(f"  - {warning}")
            lines.append("")

        lines.append("=" * 60)

        report_content = "\n".join(lines)

        with open(report_path, 'w') as f:
            f.write(report_content)

        return str(report_path)
