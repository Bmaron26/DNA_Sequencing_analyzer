"""
Input validation utilities for the bacterial mutation analysis pipeline.
"""

import os
import gzip
from pathlib import Path
from typing import List, Tuple, Optional, Dict, Any
from dataclasses import dataclass, field


@dataclass
class ValidationResult:
    """Result of a validation check."""
    is_valid: bool
    message: str
    details: Dict[str, Any] = field(default_factory=dict)
    warnings: List[str] = field(default_factory=list)


class InputValidator:
    """Validates input files for the pipeline."""

    FASTQ_EXTENSIONS = ['.fastq', '.fq', '.fastq.gz', '.fq.gz']
    FASTA_EXTENSIONS = ['.fasta', '.fa', '.fna', '.fasta.gz', '.fa.gz', '.fna.gz']
    ANNOTATION_EXTENSIONS = ['.gff', '.gff3', '.gff.gz', '.gff3.gz', '.gbk', '.gb', '.genbank']

    def __init__(self):
        self.errors: List[str] = []
        self.warnings: List[str] = []

    def validate_all(self, fastq_files: List[str], reference: str,
                    annotation: Optional[str] = None) -> ValidationResult:
        """Validate all input files."""
        self.errors = []
        self.warnings = []

        # Validate FASTQ files
        fastq_result = self.validate_fastq_files(fastq_files)
        if not fastq_result.is_valid:
            self.errors.append(fastq_result.message)

        # Validate reference
        ref_result = self.validate_reference(reference)
        if not ref_result.is_valid:
            self.errors.append(ref_result.message)

        # Validate annotation if provided
        if annotation:
            anno_result = self.validate_annotation(annotation)
            if not anno_result.is_valid:
                self.errors.append(anno_result.message)
            self.warnings.extend(anno_result.warnings)

        # Combine results
        is_valid = len(self.errors) == 0
        if is_valid:
            message = "All input files validated successfully"
        else:
            message = f"Validation failed: {'; '.join(self.errors)}"

        return ValidationResult(
            is_valid=is_valid,
            message=message,
            details={
                'fastq_files': fastq_files,
                'reference': reference,
                'annotation': annotation,
            },
            warnings=self.warnings
        )

    def validate_fastq_files(self, fastq_files: List[str]) -> ValidationResult:
        """Validate FASTQ input files."""
        if not fastq_files:
            return ValidationResult(
                is_valid=False,
                message="No FASTQ files provided"
            )

        # Check if files exist and have correct extension
        for fq_file in fastq_files:
            path = Path(fq_file)
            if not path.exists():
                return ValidationResult(
                    is_valid=False,
                    message=f"FASTQ file not found: {fq_file}"
                )

            # Check extension
            suffix = ''.join(path.suffixes).lower()
            if not any(suffix.endswith(ext) for ext in self.FASTQ_EXTENSIONS):
                return ValidationResult(
                    is_valid=False,
                    message=f"Invalid FASTQ extension: {fq_file}. Expected: {self.FASTQ_EXTENSIONS}"
                )

            # Validate FASTQ format
            result = validate_fastq(fq_file)
            if not result.is_valid:
                return result

        # Check paired-end consistency
        if len(fastq_files) == 2:
            warnings = self._check_paired_end_naming(fastq_files[0], fastq_files[1])
            return ValidationResult(
                is_valid=True,
                message="FASTQ files validated successfully",
                details={'paired': True, 'files': fastq_files},
                warnings=warnings
            )

        return ValidationResult(
            is_valid=True,
            message="FASTQ file(s) validated successfully",
            details={'paired': len(fastq_files) == 2, 'files': fastq_files}
        )

    def _check_paired_end_naming(self, file1: str, file2: str) -> List[str]:
        """Check if paired-end files follow naming conventions."""
        warnings = []
        name1 = Path(file1).stem.replace('.fastq', '').replace('.fq', '')
        name2 = Path(file2).stem.replace('.fastq', '').replace('.fq', '')

        # Common paired-end naming patterns
        patterns = [
            ('_R1', '_R2'),
            ('_1', '_2'),
            ('.R1', '.R2'),
            ('.1', '.2'),
            ('_R1_001', '_R2_001'),
        ]

        matched = False
        for p1, p2 in patterns:
            if p1 in name1 and p2 in name2:
                matched = True
                break
            if p2 in name1 and p1 in name2:
                warnings.append(f"Read files may be in reverse order: {file1} and {file2}")
                matched = True
                break

        if not matched:
            warnings.append(
                f"Paired-end files don't follow standard naming convention. "
                f"Make sure {file1} is R1 and {file2} is R2"
            )

        return warnings

    def validate_reference(self, reference_path: str) -> ValidationResult:
        """Validate reference genome file."""
        path = Path(reference_path)

        if not path.exists():
            return ValidationResult(
                is_valid=False,
                message=f"Reference file not found: {reference_path}"
            )

        suffix = ''.join(path.suffixes).lower()
        if not any(suffix.endswith(ext) for ext in self.FASTA_EXTENSIONS):
            return ValidationResult(
                is_valid=False,
                message=f"Invalid reference extension: {reference_path}. Expected FASTA format"
            )

        return validate_reference(reference_path)

    def validate_annotation(self, annotation_path: str) -> ValidationResult:
        """Validate annotation file."""
        path = Path(annotation_path)

        if not path.exists():
            return ValidationResult(
                is_valid=False,
                message=f"Annotation file not found: {annotation_path}"
            )

        suffix = ''.join(path.suffixes).lower()
        valid_ext = any(suffix.endswith(ext) for ext in self.ANNOTATION_EXTENSIONS)

        if not valid_ext:
            return ValidationResult(
                is_valid=False,
                message=f"Invalid annotation extension: {annotation_path}. "
                       f"Expected: {self.ANNOTATION_EXTENSIONS}"
            )

        warnings = []
        details = {}

        # Check format
        if any(ext in suffix for ext in ['.gff', '.gff3']):
            details['format'] = 'gff'
            result = self._validate_gff(annotation_path)
            if not result.is_valid:
                return result
            details.update(result.details)
        elif any(ext in suffix for ext in ['.gbk', '.gb', '.genbank']):
            details['format'] = 'genbank'
            result = self._validate_genbank(annotation_path)
            if not result.is_valid:
                return result
            details.update(result.details)

        return ValidationResult(
            is_valid=True,
            message="Annotation file validated successfully",
            details=details,
            warnings=warnings
        )

    def _validate_gff(self, gff_path: str) -> ValidationResult:
        """Validate GFF file format."""
        opener = gzip.open if gff_path.endswith('.gz') else open
        mode = 'rt' if gff_path.endswith('.gz') else 'r'

        feature_count = 0
        gene_count = 0
        cds_count = 0
        seqids = set()

        try:
            with opener(gff_path, mode) as f:
                for line_num, line in enumerate(f, 1):
                    if line_num > 1000:  # Sample first 1000 lines
                        break

                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue

                    parts = line.split('\t')
                    if len(parts) < 9:
                        return ValidationResult(
                            is_valid=False,
                            message=f"Invalid GFF format at line {line_num}: expected 9 columns, got {len(parts)}"
                        )

                    feature_count += 1
                    seqids.add(parts[0])

                    if parts[2].lower() == 'gene':
                        gene_count += 1
                    elif parts[2].lower() == 'cds':
                        cds_count += 1

        except Exception as e:
            return ValidationResult(
                is_valid=False,
                message=f"Error reading GFF file: {str(e)}"
            )

        return ValidationResult(
            is_valid=True,
            message="GFF file is valid",
            details={
                'feature_count': feature_count,
                'gene_count': gene_count,
                'cds_count': cds_count,
                'chromosomes': list(seqids)
            }
        )

    def _validate_genbank(self, gbk_path: str) -> ValidationResult:
        """Validate GenBank file format."""
        try:
            from Bio import SeqIO
            records = list(SeqIO.parse(gbk_path, "genbank"))

            if not records:
                return ValidationResult(
                    is_valid=False,
                    message="No records found in GenBank file"
                )

            total_features = sum(len(r.features) for r in records)

            return ValidationResult(
                is_valid=True,
                message="GenBank file is valid",
                details={
                    'record_count': len(records),
                    'total_features': total_features,
                    'record_ids': [r.id for r in records]
                }
            )

        except ImportError:
            return ValidationResult(
                is_valid=False,
                message="BioPython is required to parse GenBank files"
            )
        except Exception as e:
            return ValidationResult(
                is_valid=False,
                message=f"Error reading GenBank file: {str(e)}"
            )


def validate_fastq(fastq_path: str, sample_reads: int = 100) -> ValidationResult:
    """
    Validate FASTQ file format by sampling reads.

    Args:
        fastq_path: Path to FASTQ file
        sample_reads: Number of reads to sample for validation

    Returns:
        ValidationResult with validation status and details
    """
    path = Path(fastq_path)
    opener = gzip.open if path.suffix == '.gz' else open
    mode = 'rt' if path.suffix == '.gz' else 'r'

    read_count = 0
    total_bases = 0
    quality_scores = []
    warnings = []

    try:
        with opener(fastq_path, mode) as f:
            while read_count < sample_reads:
                # Read 4 lines (one FASTQ record)
                header = f.readline().strip()
                if not header:
                    break

                sequence = f.readline().strip()
                plus_line = f.readline().strip()
                quality = f.readline().strip()

                # Validate header
                if not header.startswith('@'):
                    return ValidationResult(
                        is_valid=False,
                        message=f"Invalid FASTQ: header doesn't start with '@' at read {read_count + 1}"
                    )

                # Validate plus line
                if not plus_line.startswith('+'):
                    return ValidationResult(
                        is_valid=False,
                        message=f"Invalid FASTQ: separator doesn't start with '+' at read {read_count + 1}"
                    )

                # Validate sequence and quality lengths match
                if len(sequence) != len(quality):
                    return ValidationResult(
                        is_valid=False,
                        message=f"Invalid FASTQ: sequence length ({len(sequence)}) != quality length ({len(quality)}) at read {read_count + 1}"
                    )

                # Validate sequence characters
                if not all(c in 'ACGTNacgtn' for c in sequence):
                    invalid_chars = set(c for c in sequence if c not in 'ACGTNacgtn')
                    return ValidationResult(
                        is_valid=False,
                        message=f"Invalid FASTQ: unexpected characters in sequence: {invalid_chars}"
                    )

                read_count += 1
                total_bases += len(sequence)

                # Sample quality scores
                for q in quality:
                    quality_scores.append(ord(q) - 33)

    except Exception as e:
        return ValidationResult(
            is_valid=False,
            message=f"Error reading FASTQ file: {str(e)}"
        )

    if read_count == 0:
        return ValidationResult(
            is_valid=False,
            message="FASTQ file appears to be empty"
        )

    # Calculate quality statistics
    mean_quality = sum(quality_scores) / len(quality_scores) if quality_scores else 0
    mean_length = total_bases / read_count

    # Check for potential issues
    if mean_quality < 20:
        warnings.append(f"Low mean quality score: {mean_quality:.1f}")
    if mean_length < 50:
        warnings.append(f"Short mean read length: {mean_length:.1f}")

    return ValidationResult(
        is_valid=True,
        message="FASTQ file validated successfully",
        details={
            'sampled_reads': read_count,
            'mean_quality': round(mean_quality, 2),
            'mean_length': round(mean_length, 2),
        },
        warnings=warnings
    )


def validate_reference(reference_path: str) -> ValidationResult:
    """
    Validate reference genome FASTA file.

    Args:
        reference_path: Path to reference FASTA file

    Returns:
        ValidationResult with validation status and details
    """
    path = Path(reference_path)
    opener = gzip.open if path.suffix == '.gz' else open
    mode = 'rt' if path.suffix == '.gz' else 'r'

    sequences = {}
    current_id = None
    current_seq = []
    warnings = []

    try:
        with opener(reference_path, mode) as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_id is not None:
                        sequences[current_id] = len(''.join(current_seq))
                    current_id = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line)

            # Don't forget the last sequence
            if current_id is not None:
                sequences[current_id] = len(''.join(current_seq))

    except Exception as e:
        return ValidationResult(
            is_valid=False,
            message=f"Error reading reference file: {str(e)}"
        )

    if not sequences:
        return ValidationResult(
            is_valid=False,
            message="No sequences found in reference file"
        )

    total_length = sum(sequences.values())

    # Check for typical bacterial genome size
    if total_length < 100000:
        warnings.append(f"Reference genome is very small ({total_length:,} bp). Is this complete?")
    elif total_length > 15000000:
        warnings.append(f"Reference genome is very large ({total_length:,} bp). Is this bacterial?")

    return ValidationResult(
        is_valid=True,
        message="Reference file validated successfully",
        details={
            'sequence_count': len(sequences),
            'total_length': total_length,
            'sequences': sequences,
        },
        warnings=warnings
    )
