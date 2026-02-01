"""
Antimicrobial Resistance (AMR) database integration module.

Supports:
- CARD (Comprehensive Antibiotic Resistance Database)
- ResFinder
- NCBI AMRFinderPlus
- Custom species-specific databases
"""

from .database import AMRDatabase, AMRHit, AMRMutation
from .annotator import AMRAnnotator
from .species_config import SpeciesConfig, get_species_config

__all__ = [
    "AMRDatabase",
    "AMRHit",
    "AMRMutation",
    "AMRAnnotator",
    "SpeciesConfig",
    "get_species_config",
]
