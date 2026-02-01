"""
File handlers for various bioinformatics file formats.
"""

import gzip
import re
from pathlib import Path
from typing import Iterator, Dict, List, Optional, Tuple, Any
from dataclasses import dataclass, field
from collections import defaultdict

try:
    from Bio import SeqIO
    from Bio.SeqFeature import SeqFeature
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False


@dataclass
class FastqRecord:
    """Represents a single FASTQ record."""
    id: str
    sequence: str
    quality: str
    description: str = ""

    @property
    def mean_quality(self) -> float:
        """Calculate mean quality score (Phred33)."""
        if not self.quality:
            return 0.0
        scores = [ord(c) - 33 for c in self.quality]
        return sum(scores) / len(scores)

    @property
    def length(self) -> int:
        """Get sequence length."""
        return len(self.sequence)

    @property
    def gc_content(self) -> float:
        """Calculate GC content."""
        if not self.sequence:
            return 0.0
        gc = sum(1 for base in self.sequence.upper() if base in 'GC')
        return gc / len(self.sequence)


class FastqReader:
    """Reader for FASTQ files (gzipped or plain)."""

    def __init__(self, filepath: str):
        self.filepath = Path(filepath)
        self.is_gzipped = self.filepath.suffix == '.gz'

    def __iter__(self) -> Iterator[FastqRecord]:
        """Iterate over FASTQ records."""
        opener = gzip.open if self.is_gzipped else open
        mode = 'rt' if self.is_gzipped else 'r'

        with opener(self.filepath, mode) as f:
            while True:
                header = f.readline().strip()
                if not header:
                    break
                if not header.startswith('@'):
                    raise ValueError(f"Invalid FASTQ format: expected '@', got '{header[0]}'")

                parts = header[1:].split(None, 1)
                read_id = parts[0]
                description = parts[1] if len(parts) > 1 else ""

                sequence = f.readline().strip()
                plus_line = f.readline().strip()
                quality = f.readline().strip()

                if not plus_line.startswith('+'):
                    raise ValueError(f"Invalid FASTQ format: expected '+', got '{plus_line}'")

                yield FastqRecord(
                    id=read_id,
                    sequence=sequence,
                    quality=quality,
                    description=description
                )

    def get_stats(self, sample_size: int = 10000) -> Dict[str, Any]:
        """Calculate statistics from a sample of reads."""
        stats = {
            'total_reads': 0,
            'total_bases': 0,
            'mean_length': 0.0,
            'mean_quality': 0.0,
            'gc_content': 0.0,
            'q20_bases': 0,
            'q30_bases': 0,
        }

        lengths = []
        qualities = []
        gc_contents = []

        for i, record in enumerate(self):
            if sample_size and i >= sample_size:
                break

            stats['total_reads'] += 1
            stats['total_bases'] += record.length
            lengths.append(record.length)
            qualities.append(record.mean_quality)
            gc_contents.append(record.gc_content)

            for q in record.quality:
                score = ord(q) - 33
                if score >= 20:
                    stats['q20_bases'] += 1
                if score >= 30:
                    stats['q30_bases'] += 1

        if lengths:
            stats['mean_length'] = sum(lengths) / len(lengths)
            stats['mean_quality'] = sum(qualities) / len(qualities)
            stats['gc_content'] = sum(gc_contents) / len(gc_contents)

        return stats


@dataclass
class FastaRecord:
    """Represents a single FASTA record."""
    id: str
    sequence: str
    description: str = ""

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def gc_content(self) -> float:
        if not self.sequence:
            return 0.0
        gc = sum(1 for base in self.sequence.upper() if base in 'GC')
        return gc / len(self.sequence)


class FastaReader:
    """Reader for FASTA files (gzipped or plain)."""

    def __init__(self, filepath: str):
        self.filepath = Path(filepath)
        self.is_gzipped = self.filepath.suffix == '.gz'
        self._sequences: Dict[str, FastaRecord] = {}

    def __iter__(self) -> Iterator[FastaRecord]:
        """Iterate over FASTA records."""
        opener = gzip.open if self.is_gzipped else open
        mode = 'rt' if self.is_gzipped else 'r'

        with opener(self.filepath, mode) as f:
            header = None
            sequence_parts = []

            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if header is not None:
                        parts = header[1:].split(None, 1)
                        yield FastaRecord(
                            id=parts[0],
                            sequence=''.join(sequence_parts),
                            description=parts[1] if len(parts) > 1 else ""
                        )
                    header = line
                    sequence_parts = []
                else:
                    sequence_parts.append(line)

            if header is not None:
                parts = header[1:].split(None, 1)
                yield FastaRecord(
                    id=parts[0],
                    sequence=''.join(sequence_parts),
                    description=parts[1] if len(parts) > 1 else ""
                )

    def load_all(self) -> Dict[str, FastaRecord]:
        """Load all sequences into memory."""
        if not self._sequences:
            for record in self:
                self._sequences[record.id] = record
        return self._sequences

    def get_sequence(self, seq_id: str) -> Optional[FastaRecord]:
        """Get a specific sequence by ID."""
        if not self._sequences:
            self.load_all()
        return self._sequences.get(seq_id)


@dataclass
class GffFeature:
    """Represents a GFF feature."""
    seqid: str
    source: str
    feature_type: str
    start: int
    end: int
    score: Optional[float]
    strand: str
    phase: Optional[int]
    attributes: Dict[str, str]

    @property
    def gene_id(self) -> Optional[str]:
        return self.attributes.get('ID') or self.attributes.get('gene_id')

    @property
    def gene_name(self) -> Optional[str]:
        return self.attributes.get('Name') or self.attributes.get('gene_name') or self.attributes.get('gene')

    @property
    def product(self) -> Optional[str]:
        return self.attributes.get('product') or self.attributes.get('Product')

    @property
    def locus_tag(self) -> Optional[str]:
        return self.attributes.get('locus_tag') or self.attributes.get('ID')


class GffParser:
    """Parser for GFF3 annotation files."""

    def __init__(self, filepath: str):
        self.filepath = Path(filepath)
        self.features: List[GffFeature] = []
        self.features_by_type: Dict[str, List[GffFeature]] = defaultdict(list)
        self.features_by_location: Dict[str, List[GffFeature]] = defaultdict(list)
        self._parsed = False

    def parse(self) -> List[GffFeature]:
        """Parse the GFF file."""
        if self._parsed:
            return self.features

        opener = gzip.open if self.filepath.suffix == '.gz' else open
        mode = 'rt' if self.filepath.suffix == '.gz' else 'r'

        with opener(self.filepath, mode) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                parts = line.split('\t')
                if len(parts) < 9:
                    continue

                # Parse attributes
                attributes = {}
                for attr in parts[8].split(';'):
                    attr = attr.strip()
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        attributes[key] = value.replace('%20', ' ').replace('%2C', ',')

                feature = GffFeature(
                    seqid=parts[0],
                    source=parts[1],
                    feature_type=parts[2],
                    start=int(parts[3]),
                    end=int(parts[4]),
                    score=float(parts[5]) if parts[5] != '.' else None,
                    strand=parts[6],
                    phase=int(parts[7]) if parts[7] != '.' else None,
                    attributes=attributes
                )

                self.features.append(feature)
                self.features_by_type[feature.feature_type].append(feature)
                self.features_by_location[feature.seqid].append(feature)

        self._parsed = True
        return self.features

    def get_features_at_position(self, seqid: str, position: int) -> List[GffFeature]:
        """Get all features overlapping a specific position."""
        if not self._parsed:
            self.parse()

        overlapping = []
        for feature in self.features_by_location.get(seqid, []):
            if feature.start <= position <= feature.end:
                overlapping.append(feature)
        return overlapping

    def get_genes(self) -> List[GffFeature]:
        """Get all gene features."""
        if not self._parsed:
            self.parse()
        return self.features_by_type.get('gene', []) + self.features_by_type.get('CDS', [])

    def get_nearest_feature(self, seqid: str, position: int,
                           feature_types: Optional[List[str]] = None) -> Tuple[Optional[GffFeature], int]:
        """Find the nearest feature to a position."""
        if not self._parsed:
            self.parse()

        nearest = None
        min_distance = float('inf')

        features = self.features_by_location.get(seqid, [])
        for feature in features:
            if feature_types and feature.feature_type not in feature_types:
                continue

            # Check if position is within feature
            if feature.start <= position <= feature.end:
                return feature, 0

            # Calculate distance
            if position < feature.start:
                distance = feature.start - position
            else:
                distance = position - feature.end

            if distance < min_distance:
                min_distance = distance
                nearest = feature

        return nearest, int(min_distance) if nearest else -1


class GenbankParser:
    """Parser for GenBank annotation files using BioPython."""

    def __init__(self, filepath: str):
        self.filepath = Path(filepath)
        self.records = []
        self.features: List[GffFeature] = []
        self._parsed = False

    def parse(self) -> List[GffFeature]:
        """Parse the GenBank file and convert to GFF-like features."""
        if self._parsed:
            return self.features

        if not BIOPYTHON_AVAILABLE:
            raise ImportError("BioPython is required for GenBank parsing")

        for record in SeqIO.parse(str(self.filepath), "genbank"):
            self.records.append(record)

            for feature in record.features:
                if feature.type in ['source', 'misc_feature']:
                    continue

                # Extract attributes
                attributes = {}
                for key, values in feature.qualifiers.items():
                    attributes[key] = values[0] if len(values) == 1 else ';'.join(values)

                # Get ID
                if 'locus_tag' in attributes:
                    attributes['ID'] = attributes['locus_tag']
                elif 'gene' in attributes:
                    attributes['ID'] = attributes['gene']

                gff_feature = GffFeature(
                    seqid=record.id,
                    source='genbank',
                    feature_type=feature.type,
                    start=int(feature.location.start) + 1,  # Convert to 1-based
                    end=int(feature.location.end),
                    score=None,
                    strand='+' if feature.location.strand == 1 else '-',
                    phase=None,
                    attributes=attributes
                )
                self.features.append(gff_feature)

        self._parsed = True
        return self.features


@dataclass
class VcfVariant:
    """Represents a VCF variant."""
    chrom: str
    pos: int
    id: str
    ref: str
    alt: str
    qual: float
    filter: str
    info: Dict[str, Any]
    format_fields: List[str] = field(default_factory=list)
    samples: Dict[str, Dict[str, Any]] = field(default_factory=dict)

    @property
    def variant_type(self) -> str:
        """Determine variant type (SNP, insertion, deletion, complex)."""
        if len(self.ref) == len(self.alt) == 1:
            return 'SNP'
        elif len(self.ref) < len(self.alt):
            return 'insertion'
        elif len(self.ref) > len(self.alt):
            return 'deletion'
        else:
            return 'complex'

    @property
    def depth(self) -> int:
        """Get total depth from INFO field."""
        return int(self.info.get('DP', 0))

    @property
    def allele_frequency(self) -> float:
        """Get allele frequency."""
        af = self.info.get('AF', self.info.get('VAF', None))
        if af is not None:
            return float(af)
        # Calculate from AD if available
        if 'AD' in self.info:
            ad = self.info['AD']
            if isinstance(ad, str):
                ad = [int(x) for x in ad.split(',')]
            if sum(ad) > 0:
                return ad[1] / sum(ad)
        return 0.0


class VcfParser:
    """Parser for VCF files."""

    def __init__(self, filepath: str):
        self.filepath = Path(filepath)
        self.header_lines: List[str] = []
        self.sample_names: List[str] = []
        self.variants: List[VcfVariant] = []
        self._parsed = False

    def parse(self) -> List[VcfVariant]:
        """Parse the VCF file."""
        if self._parsed:
            return self.variants

        opener = gzip.open if self.filepath.suffix == '.gz' else open
        mode = 'rt' if self.filepath.suffix == '.gz' else 'r'

        with opener(self.filepath, mode) as f:
            for line in f:
                line = line.strip()
                if line.startswith('##'):
                    self.header_lines.append(line)
                    continue
                elif line.startswith('#CHROM'):
                    parts = line.split('\t')
                    if len(parts) > 9:
                        self.sample_names = parts[9:]
                    continue

                parts = line.split('\t')
                if len(parts) < 8:
                    continue

                # Parse INFO field
                info = {}
                for item in parts[7].split(';'):
                    if '=' in item:
                        key, value = item.split('=', 1)
                        # Try to parse numeric values
                        try:
                            if ',' in value:
                                value = [float(v) if '.' in v else int(v) for v in value.split(',')]
                            elif '.' in value:
                                value = float(value)
                            else:
                                value = int(value)
                        except ValueError:
                            pass
                        info[key] = value
                    else:
                        info[item] = True

                # Parse FORMAT and samples
                format_fields = []
                samples = {}
                if len(parts) > 8:
                    format_fields = parts[8].split(':')
                    for i, sample_name in enumerate(self.sample_names):
                        if len(parts) > 9 + i:
                            sample_values = parts[9 + i].split(':')
                            sample_data = {}
                            for j, field_name in enumerate(format_fields):
                                if j < len(sample_values):
                                    sample_data[field_name] = sample_values[j]
                            samples[sample_name] = sample_data

                variant = VcfVariant(
                    chrom=parts[0],
                    pos=int(parts[1]),
                    id=parts[2] if parts[2] != '.' else '',
                    ref=parts[3],
                    alt=parts[4],
                    qual=float(parts[5]) if parts[5] != '.' else 0.0,
                    filter=parts[6],
                    info=info,
                    format_fields=format_fields,
                    samples=samples
                )
                self.variants.append(variant)

        self._parsed = True
        return self.variants

    def __iter__(self) -> Iterator[VcfVariant]:
        """Iterate over variants."""
        if not self._parsed:
            self.parse()
        return iter(self.variants)

    def filter_variants(self, min_qual: float = 0, min_depth: int = 0,
                       variant_types: Optional[List[str]] = None) -> List[VcfVariant]:
        """Filter variants based on criteria."""
        if not self._parsed:
            self.parse()

        filtered = []
        for variant in self.variants:
            if variant.qual < min_qual:
                continue
            if variant.depth < min_depth:
                continue
            if variant_types and variant.variant_type not in variant_types:
                continue
            filtered.append(variant)

        return filtered
