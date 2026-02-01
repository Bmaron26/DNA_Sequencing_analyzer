"""
Sample metadata module for tracking experimental design.

Supports:
- Treatment groups (e.g., different AMPs)
- Sample origin (single colony vs population)
- Replicate tracking
- Time points
"""

import json
import csv
from pathlib import Path
from dataclasses import dataclass, field, asdict
from typing import List, Dict, Any, Optional
from enum import Enum
from datetime import datetime
import logging

logger = logging.getLogger(__name__)


class SampleOrigin(Enum):
    """Sample origin type."""
    SINGLE_COLONY = "single_colony"
    POPULATION = "population"  # Multiple colonies (~10)
    MIXED = "mixed"
    UNKNOWN = "unknown"

    @classmethod
    def from_string(cls, value: str) -> "SampleOrigin":
        """Create from string value."""
        value_lower = value.lower().replace(" ", "_").replace("-", "_")
        for member in cls:
            if member.value == value_lower:
                return member
        return cls.UNKNOWN


@dataclass
class TreatmentGroup:
    """Represents a treatment/experimental group."""
    name: str
    description: str = ""
    treatment_type: str = ""  # e.g., "AMP", "antibiotic", "control"
    treatment_agent: str = ""  # e.g., "Melittin", "Cecropin"
    concentration: Optional[str] = None
    duration: Optional[str] = None
    metadata: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            'name': self.name,
            'description': self.description,
            'treatment_type': self.treatment_type,
            'treatment_agent': self.treatment_agent,
            'concentration': self.concentration,
            'duration': self.duration,
            'metadata': self.metadata,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "TreatmentGroup":
        return cls(
            name=data.get('name', ''),
            description=data.get('description', ''),
            treatment_type=data.get('treatment_type', ''),
            treatment_agent=data.get('treatment_agent', ''),
            concentration=data.get('concentration'),
            duration=data.get('duration'),
            metadata=data.get('metadata', {}),
        )


@dataclass
class SampleMetadata:
    """Metadata for a single sample."""
    sample_id: str
    group: str  # Treatment group name
    replicate: int = 1
    origin: SampleOrigin = SampleOrigin.UNKNOWN
    colony_count: int = 1  # Number of colonies if population
    passage_number: Optional[int] = None
    timepoint: Optional[str] = None
    parent_sample: Optional[str] = None  # For tracking lineage
    species: str = ""
    strain: str = ""
    fastq_r1: str = ""
    fastq_r2: str = ""
    notes: str = ""
    metadata: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            'sample_id': self.sample_id,
            'group': self.group,
            'replicate': self.replicate,
            'origin': self.origin.value,
            'colony_count': self.colony_count,
            'passage_number': self.passage_number,
            'timepoint': self.timepoint,
            'parent_sample': self.parent_sample,
            'species': self.species,
            'strain': self.strain,
            'fastq_r1': self.fastq_r1,
            'fastq_r2': self.fastq_r2,
            'notes': self.notes,
            'metadata': self.metadata,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "SampleMetadata":
        origin = data.get('origin', 'unknown')
        if isinstance(origin, str):
            origin = SampleOrigin.from_string(origin)

        return cls(
            sample_id=data.get('sample_id', ''),
            group=data.get('group', ''),
            replicate=data.get('replicate', 1),
            origin=origin,
            colony_count=data.get('colony_count', 1),
            passage_number=data.get('passage_number'),
            timepoint=data.get('timepoint'),
            parent_sample=data.get('parent_sample'),
            species=data.get('species', ''),
            strain=data.get('strain', ''),
            fastq_r1=data.get('fastq_r1', ''),
            fastq_r2=data.get('fastq_r2', ''),
            notes=data.get('notes', ''),
            metadata=data.get('metadata', {}),
        )

    @property
    def is_population(self) -> bool:
        """Check if sample is from population (multiple colonies)."""
        return self.origin == SampleOrigin.POPULATION or self.colony_count > 1


@dataclass
class Experiment:
    """
    Represents a complete experiment with multiple samples and groups.
    """
    name: str
    description: str = ""
    species: str = ""
    strain: str = ""
    reference_genome: str = ""
    annotation_file: str = ""
    created_date: str = field(default_factory=lambda: datetime.now().isoformat())
    groups: Dict[str, TreatmentGroup] = field(default_factory=dict)
    samples: Dict[str, SampleMetadata] = field(default_factory=dict)
    metadata: Dict[str, Any] = field(default_factory=dict)

    def add_group(self, group: TreatmentGroup) -> None:
        """Add a treatment group."""
        self.groups[group.name] = group
        logger.info(f"Added group: {group.name}")

    def add_sample(self, sample: SampleMetadata) -> None:
        """Add a sample to the experiment."""
        if sample.group and sample.group not in self.groups:
            logger.warning(f"Sample {sample.sample_id} references unknown group: {sample.group}")

        # Inherit species/strain from experiment if not set
        if not sample.species and self.species:
            sample.species = self.species
        if not sample.strain and self.strain:
            sample.strain = self.strain

        self.samples[sample.sample_id] = sample
        logger.info(f"Added sample: {sample.sample_id} to group: {sample.group}")

    def get_samples_by_group(self, group_name: str) -> List[SampleMetadata]:
        """Get all samples in a specific group."""
        return [s for s in self.samples.values() if s.group == group_name]

    def get_group_names(self) -> List[str]:
        """Get list of all group names."""
        return list(self.groups.keys())

    def get_replicates(self, group_name: str) -> Dict[int, SampleMetadata]:
        """Get samples by replicate number for a group."""
        samples = self.get_samples_by_group(group_name)
        return {s.replicate: s for s in samples}

    def to_dict(self) -> Dict[str, Any]:
        return {
            'name': self.name,
            'description': self.description,
            'species': self.species,
            'strain': self.strain,
            'reference_genome': self.reference_genome,
            'annotation_file': self.annotation_file,
            'created_date': self.created_date,
            'groups': {k: v.to_dict() for k, v in self.groups.items()},
            'samples': {k: v.to_dict() for k, v in self.samples.items()},
            'metadata': self.metadata,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "Experiment":
        exp = cls(
            name=data.get('name', ''),
            description=data.get('description', ''),
            species=data.get('species', ''),
            strain=data.get('strain', ''),
            reference_genome=data.get('reference_genome', ''),
            annotation_file=data.get('annotation_file', ''),
            created_date=data.get('created_date', datetime.now().isoformat()),
            metadata=data.get('metadata', {}),
        )

        for name, group_data in data.get('groups', {}).items():
            exp.groups[name] = TreatmentGroup.from_dict(group_data)

        for sample_id, sample_data in data.get('samples', {}).items():
            exp.samples[sample_id] = SampleMetadata.from_dict(sample_data)

        return exp

    def save(self, filepath: str) -> None:
        """Save experiment to JSON file."""
        with open(filepath, 'w') as f:
            json.dump(self.to_dict(), f, indent=2)
        logger.info(f"Saved experiment to {filepath}")

    @classmethod
    def load(cls, filepath: str) -> "Experiment":
        """Load experiment from JSON file."""
        with open(filepath, 'r') as f:
            data = json.load(f)
        return cls.from_dict(data)

    @classmethod
    def from_csv(cls, filepath: str, name: str = "experiment") -> "Experiment":
        """
        Load experiment from CSV sample sheet.

        Expected columns:
        - sample_id (required)
        - group (required)
        - replicate
        - origin (single_colony/population)
        - colony_count
        - fastq_r1
        - fastq_r2
        - species
        - strain
        - passage_number
        - timepoint
        - notes
        """
        exp = cls(name=name)

        with open(filepath, 'r') as f:
            reader = csv.DictReader(f)

            for row in reader:
                # Create group if not exists
                group_name = row.get('group', 'default')
                if group_name and group_name not in exp.groups:
                    exp.groups[group_name] = TreatmentGroup(
                        name=group_name,
                        treatment_agent=row.get('treatment_agent', ''),
                    )

                # Create sample
                origin = SampleOrigin.from_string(row.get('origin', 'unknown'))

                sample = SampleMetadata(
                    sample_id=row.get('sample_id', ''),
                    group=group_name,
                    replicate=int(row.get('replicate', 1)),
                    origin=origin,
                    colony_count=int(row.get('colony_count', 1) or 1),
                    passage_number=int(row['passage_number']) if row.get('passage_number') else None,
                    timepoint=row.get('timepoint'),
                    parent_sample=row.get('parent_sample'),
                    species=row.get('species', ''),
                    strain=row.get('strain', ''),
                    fastq_r1=row.get('fastq_r1', ''),
                    fastq_r2=row.get('fastq_r2', ''),
                    notes=row.get('notes', ''),
                )
                exp.add_sample(sample)

        return exp

    def to_csv(self, filepath: str) -> None:
        """Export experiment to CSV sample sheet."""
        fieldnames = [
            'sample_id', 'group', 'replicate', 'origin', 'colony_count',
            'species', 'strain', 'passage_number', 'timepoint',
            'parent_sample', 'fastq_r1', 'fastq_r2', 'notes'
        ]

        with open(filepath, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()

            for sample in self.samples.values():
                row = {
                    'sample_id': sample.sample_id,
                    'group': sample.group,
                    'replicate': sample.replicate,
                    'origin': sample.origin.value,
                    'colony_count': sample.colony_count,
                    'species': sample.species,
                    'strain': sample.strain,
                    'passage_number': sample.passage_number or '',
                    'timepoint': sample.timepoint or '',
                    'parent_sample': sample.parent_sample or '',
                    'fastq_r1': sample.fastq_r1,
                    'fastq_r2': sample.fastq_r2,
                    'notes': sample.notes,
                }
                writer.writerow(row)


def create_example_experiment() -> Experiment:
    """Create an example experiment for testing."""
    exp = Experiment(
        name="AMP_Evolution_Study",
        description="Evolution of S. aureus under AMP pressure",
        species="Staphylococcus aureus",
        strain="ATCC 29213",
        reference_genome="reference/s_aureus.fasta",
        annotation_file="reference/s_aureus.gff",
    )

    # Add treatment groups
    exp.add_group(TreatmentGroup(
        name="Control",
        description="No treatment control",
        treatment_type="control",
    ))

    exp.add_group(TreatmentGroup(
        name="Melittin",
        description="Melittin-evolved strains",
        treatment_type="AMP",
        treatment_agent="Melittin",
        concentration="0.5x MIC",
        duration="30 passages",
    ))

    exp.add_group(TreatmentGroup(
        name="Cecropin",
        description="Cecropin A-evolved strains",
        treatment_type="AMP",
        treatment_agent="Cecropin A",
        concentration="0.5x MIC",
        duration="30 passages",
    ))

    # Add samples
    for i in range(1, 7):
        exp.add_sample(SampleMetadata(
            sample_id=f"Mel{i}",
            group="Melittin",
            replicate=i,
            origin=SampleOrigin.POPULATION,
            colony_count=10,
            passage_number=30,
            fastq_r1=f"fastq/Mel{i}_R1.fastq.gz",
            fastq_r2=f"fastq/Mel{i}_R2.fastq.gz",
        ))

    for i in range(1, 7):
        exp.add_sample(SampleMetadata(
            sample_id=f"Cec{i}",
            group="Cecropin",
            replicate=i,
            origin=SampleOrigin.POPULATION,
            colony_count=10,
            passage_number=30,
            fastq_r1=f"fastq/Cec{i}_R1.fastq.gz",
            fastq_r2=f"fastq/Cec{i}_R2.fastq.gz",
        ))

    for i in range(1, 4):
        exp.add_sample(SampleMetadata(
            sample_id=f"Ctrl{i}",
            group="Control",
            replicate=i,
            origin=SampleOrigin.POPULATION,
            colony_count=10,
            passage_number=30,
            fastq_r1=f"fastq/Ctrl{i}_R1.fastq.gz",
            fastq_r2=f"fastq/Ctrl{i}_R2.fastq.gz",
        ))

    return exp
