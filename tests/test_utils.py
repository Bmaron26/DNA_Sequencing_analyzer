"""Tests for utility modules."""

import pytest
import tempfile
import os
from pathlib import Path

from bacterial_mutation_analyzer.utils.file_handlers import (
    FastqReader, FastaReader, FastqRecord, FastaRecord
)
from bacterial_mutation_analyzer.utils.validators import (
    InputValidator, validate_fastq, validate_reference
)


class TestFastqReader:
    """Tests for FastqReader."""

    def test_read_fastq(self, tmp_path):
        """Test reading a simple FASTQ file."""
        fastq_content = """@read1
ACGTACGTACGT
+
IIIIIIIIIIII
@read2
GCTAGCTAGCTA
+
HHHHHHHHHHHH
"""
        fastq_file = tmp_path / "test.fastq"
        fastq_file.write_text(fastq_content)

        reader = FastqReader(str(fastq_file))
        records = list(reader)

        assert len(records) == 2
        assert records[0].id == "read1"
        assert records[0].sequence == "ACGTACGTACGT"
        assert records[1].id == "read2"

    def test_fastq_record_properties(self):
        """Test FastqRecord properties."""
        record = FastqRecord(
            id="test",
            sequence="ACGTACGT",
            quality="IIIIIIII"
        )

        assert record.length == 8
        assert record.gc_content == 0.5
        assert record.mean_quality == 40.0  # 'I' = 40


class TestFastaReader:
    """Tests for FastaReader."""

    def test_read_fasta(self, tmp_path):
        """Test reading a simple FASTA file."""
        fasta_content = """>seq1 description
ACGTACGTACGT
GCTAGCTAGCTA
>seq2
AAAATTTTGGGG
"""
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(fasta_content)

        reader = FastaReader(str(fasta_file))
        records = list(reader)

        assert len(records) == 2
        assert records[0].id == "seq1"
        assert records[0].description == "description"
        assert len(records[0].sequence) == 24
        assert records[1].id == "seq2"


class TestValidators:
    """Tests for input validators."""

    def test_validate_fastq_valid(self, tmp_path):
        """Test validation of valid FASTQ file."""
        fastq_content = """@read1
ACGTACGTACGT
+
IIIIIIIIIIII
"""
        fastq_file = tmp_path / "test.fastq"
        fastq_file.write_text(fastq_content)

        result = validate_fastq(str(fastq_file))
        assert result.is_valid

    def test_validate_fastq_invalid(self, tmp_path):
        """Test validation of invalid FASTQ file."""
        # Invalid: header doesn't start with @
        fastq_content = """read1
ACGTACGTACGT
+
IIIIIIIIIIII
"""
        fastq_file = tmp_path / "test.fastq"
        fastq_file.write_text(fastq_content)

        result = validate_fastq(str(fastq_file))
        assert not result.is_valid

    def test_validate_reference_valid(self, tmp_path):
        """Test validation of valid reference FASTA."""
        fasta_content = """>chromosome1
ACGTACGTACGTACGTACGT
GCTAGCTAGCTAGCTAGCTA
"""
        fasta_file = tmp_path / "reference.fasta"
        fasta_file.write_text(fasta_content)

        result = validate_reference(str(fasta_file))
        assert result.is_valid
        assert result.details['sequence_count'] == 1


class TestInputValidator:
    """Tests for InputValidator."""

    def test_validate_all(self, tmp_path):
        """Test complete input validation."""
        # Create test files
        fastq_content = "@read1\nACGT\n+\nIIII\n"
        fasta_content = ">seq1\nACGTACGT\n"

        fastq_file = tmp_path / "reads.fastq"
        fasta_file = tmp_path / "ref.fasta"

        fastq_file.write_text(fastq_content)
        fasta_file.write_text(fasta_content)

        validator = InputValidator()
        result = validator.validate_all(
            fastq_files=[str(fastq_file)],
            reference=str(fasta_file)
        )

        assert result.is_valid


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
