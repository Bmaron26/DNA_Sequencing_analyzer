"""
Species-specific configuration for AMR analysis.

Provides species-specific:
- Known resistance genes
- Important regulatory regions
- Typical mutation hotspots
- Resistance mechanisms
"""

from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional
import json
from pathlib import Path


@dataclass
class ResistanceGene:
    """Known resistance gene information."""
    gene_name: str
    locus_tag: str = ""
    product: str = ""
    resistance_class: str = ""  # e.g., "beta-lactam", "aminoglycoside", "AMP"
    mechanism: str = ""  # e.g., "efflux", "modification", "target alteration"
    known_mutations: List[str] = field(default_factory=list)  # e.g., ["S84L", "E88K"]
    notes: str = ""


@dataclass
class SpeciesConfig:
    """Species-specific configuration for AMR analysis."""
    species: str
    strain: str = ""
    ncbi_taxid: str = ""

    # Reference genome info
    genome_size: int = 0
    gc_content: float = 0.0
    chromosome_count: int = 1

    # Known resistance genes
    resistance_genes: Dict[str, ResistanceGene] = field(default_factory=dict)

    # AMP resistance specific
    amp_resistance_genes: List[str] = field(default_factory=list)
    membrane_genes: List[str] = field(default_factory=list)
    regulatory_genes: List[str] = field(default_factory=list)

    # Important regions
    hypermutable_regions: List[str] = field(default_factory=list)
    essential_genes: List[str] = field(default_factory=list)

    # AMR databases to query
    amr_databases: List[str] = field(default_factory=lambda: ["card", "resfinder"])

    # Additional metadata
    metadata: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            'species': self.species,
            'strain': self.strain,
            'ncbi_taxid': self.ncbi_taxid,
            'genome_size': self.genome_size,
            'gc_content': self.gc_content,
            'chromosome_count': self.chromosome_count,
            'resistance_genes': {
                k: {
                    'gene_name': v.gene_name,
                    'locus_tag': v.locus_tag,
                    'product': v.product,
                    'resistance_class': v.resistance_class,
                    'mechanism': v.mechanism,
                    'known_mutations': v.known_mutations,
                    'notes': v.notes,
                }
                for k, v in self.resistance_genes.items()
            },
            'amp_resistance_genes': self.amp_resistance_genes,
            'membrane_genes': self.membrane_genes,
            'regulatory_genes': self.regulatory_genes,
            'hypermutable_regions': self.hypermutable_regions,
            'essential_genes': self.essential_genes,
            'amr_databases': self.amr_databases,
            'metadata': self.metadata,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "SpeciesConfig":
        config = cls(
            species=data.get('species', ''),
            strain=data.get('strain', ''),
            ncbi_taxid=data.get('ncbi_taxid', ''),
            genome_size=data.get('genome_size', 0),
            gc_content=data.get('gc_content', 0.0),
            chromosome_count=data.get('chromosome_count', 1),
            amp_resistance_genes=data.get('amp_resistance_genes', []),
            membrane_genes=data.get('membrane_genes', []),
            regulatory_genes=data.get('regulatory_genes', []),
            hypermutable_regions=data.get('hypermutable_regions', []),
            essential_genes=data.get('essential_genes', []),
            amr_databases=data.get('amr_databases', ['card', 'resfinder']),
            metadata=data.get('metadata', {}),
        )

        for gene_name, gene_data in data.get('resistance_genes', {}).items():
            config.resistance_genes[gene_name] = ResistanceGene(
                gene_name=gene_data.get('gene_name', gene_name),
                locus_tag=gene_data.get('locus_tag', ''),
                product=gene_data.get('product', ''),
                resistance_class=gene_data.get('resistance_class', ''),
                mechanism=gene_data.get('mechanism', ''),
                known_mutations=gene_data.get('known_mutations', []),
                notes=gene_data.get('notes', ''),
            )

        return config

    def save(self, filepath: str) -> None:
        """Save configuration to JSON file."""
        with open(filepath, 'w') as f:
            json.dump(self.to_dict(), f, indent=2)

    @classmethod
    def load(cls, filepath: str) -> "SpeciesConfig":
        """Load configuration from JSON file."""
        with open(filepath, 'r') as f:
            data = json.load(f)
        return cls.from_dict(data)

    def is_resistance_gene(self, gene_name: str) -> bool:
        """Check if gene is a known resistance gene."""
        gene_lower = gene_name.lower()
        return (
            gene_lower in [g.lower() for g in self.resistance_genes.keys()] or
            gene_lower in [g.lower() for g in self.amp_resistance_genes]
        )

    def get_gene_info(self, gene_name: str) -> Optional[ResistanceGene]:
        """Get resistance gene information."""
        for name, gene in self.resistance_genes.items():
            if name.lower() == gene_name.lower():
                return gene
            if gene.locus_tag.lower() == gene_name.lower():
                return gene
        return None

    def is_membrane_related(self, gene_name: str) -> bool:
        """Check if gene is membrane-related (relevant for AMP resistance)."""
        return gene_name.lower() in [g.lower() for g in self.membrane_genes]


# Pre-defined species configurations

def get_staphylococcus_aureus_config() -> SpeciesConfig:
    """Get S. aureus specific configuration."""
    config = SpeciesConfig(
        species="Staphylococcus aureus",
        ncbi_taxid="1280",
        genome_size=2800000,
        gc_content=0.33,
        chromosome_count=1,
    )

    # Known resistance genes
    config.resistance_genes = {
        'mecA': ResistanceGene(
            gene_name='mecA',
            product='Penicillin-binding protein 2a',
            resistance_class='beta-lactam',
            mechanism='target alteration',
            notes='MRSA marker',
        ),
        'vanA': ResistanceGene(
            gene_name='vanA',
            product='D-alanine--D-lactate ligase',
            resistance_class='glycopeptide',
            mechanism='target alteration',
        ),
        'dfrA': ResistanceGene(
            gene_name='dfrA',
            product='Dihydrofolate reductase',
            resistance_class='trimethoprim',
            mechanism='target bypass',
        ),
        'norA': ResistanceGene(
            gene_name='norA',
            product='Quinolone resistance protein',
            resistance_class='fluoroquinolone',
            mechanism='efflux',
        ),
        'mprF': ResistanceGene(
            gene_name='mprF',
            product='Phosphatidylglycerol lysyltransferase',
            resistance_class='AMP',
            mechanism='membrane modification',
            known_mutations=['S295L', 'L341F', 'T345A'],
            notes='Key gene for AMP resistance via membrane charge modification',
        ),
        'dltA': ResistanceGene(
            gene_name='dltA',
            product='D-alanine--D-alanyl carrier protein ligase',
            resistance_class='AMP',
            mechanism='membrane modification',
            notes='dlt operon - teichoic acid D-alanylation',
        ),
        'vraS': ResistanceGene(
            gene_name='vraS',
            product='Sensor histidine kinase VraS',
            resistance_class='cell wall stress',
            mechanism='regulatory',
            notes='VraSR two-component system',
        ),
        'graS': ResistanceGene(
            gene_name='graS',
            product='Sensor histidine kinase GraS',
            resistance_class='AMP',
            mechanism='regulatory',
            notes='GraSR system regulates mprF and dlt operon',
        ),
        'rpoB': ResistanceGene(
            gene_name='rpoB',
            product='RNA polymerase beta subunit',
            resistance_class='rifampicin',
            mechanism='target alteration',
            known_mutations=['H481Y', 'H481N', 'S486L'],
        ),
        'gyrA': ResistanceGene(
            gene_name='gyrA',
            product='DNA gyrase subunit A',
            resistance_class='fluoroquinolone',
            mechanism='target alteration',
            known_mutations=['S84L', 'E88K'],
        ),
    }

    # AMP resistance specific genes
    config.amp_resistance_genes = [
        'mprF', 'dltA', 'dltB', 'dltC', 'dltD',  # Membrane charge modification
        'graS', 'graR', 'vraS', 'vraR',  # Two-component systems
        'pgsA', 'cls', 'pssA',  # Phospholipid synthesis
        'liaF', 'liaS', 'liaR',  # Cell envelope stress response
    ]

    # Membrane-related genes
    config.membrane_genes = [
        'mprF', 'pgsA', 'cls', 'pssA', 'cdsA',  # Phospholipid synthesis
        'dltA', 'dltB', 'dltC', 'dltD',  # Teichoic acid modification
        'tagO', 'tarO', 'ltaS',  # Wall teichoic acid
        'atl', 'lytM',  # Autolysins
    ]

    # Regulatory genes
    config.regulatory_genes = [
        'agrA', 'agrB', 'agrC', 'agrD',  # Agr quorum sensing
        'sarA', 'sarS', 'sarR',  # Sar regulators
        'sigB', 'rsbU', 'rsbV', 'rsbW',  # Sigma B
        'codY', 'ccpA',  # Metabolic regulators
        'walK', 'walR',  # WalKR essential TCS
    ]

    return config


def get_escherichia_coli_config() -> SpeciesConfig:
    """Get E. coli specific configuration."""
    config = SpeciesConfig(
        species="Escherichia coli",
        ncbi_taxid="562",
        genome_size=4600000,
        gc_content=0.51,
        chromosome_count=1,
    )

    config.resistance_genes = {
        'blaTEM': ResistanceGene(
            gene_name='blaTEM',
            product='Beta-lactamase TEM',
            resistance_class='beta-lactam',
            mechanism='enzymatic inactivation',
        ),
        'blaCTX-M': ResistanceGene(
            gene_name='blaCTX-M',
            product='CTX-M beta-lactamase',
            resistance_class='beta-lactam',
            mechanism='enzymatic inactivation',
        ),
        'acrA': ResistanceGene(
            gene_name='acrA',
            product='Efflux pump membrane fusion protein',
            resistance_class='multidrug',
            mechanism='efflux',
        ),
        'acrB': ResistanceGene(
            gene_name='acrB',
            product='Multidrug efflux pump',
            resistance_class='multidrug',
            mechanism='efflux',
        ),
        'pmrA': ResistanceGene(
            gene_name='pmrA',
            product='Response regulator PmrA',
            resistance_class='AMP',
            mechanism='LPS modification',
            notes='Regulates lipid A modification',
        ),
        'pmrB': ResistanceGene(
            gene_name='pmrB',
            product='Sensor kinase PmrB',
            resistance_class='AMP',
            mechanism='LPS modification',
        ),
        'phoP': ResistanceGene(
            gene_name='phoP',
            product='Response regulator PhoP',
            resistance_class='AMP',
            mechanism='LPS modification',
        ),
        'arnT': ResistanceGene(
            gene_name='arnT',
            product='4-amino-4-deoxy-L-arabinose transferase',
            resistance_class='AMP',
            mechanism='LPS modification',
            notes='Adds L-Ara4N to lipid A',
        ),
    }

    config.amp_resistance_genes = [
        'pmrA', 'pmrB', 'pmrC', 'pmrD',  # PmrAB system
        'phoP', 'phoQ',  # PhoPQ system
        'arnA', 'arnB', 'arnC', 'arnT',  # L-Ara4N modification
        'eptA', 'eptB',  # Phosphoethanolamine addition
        'lpxM', 'lpxP',  # Lipid A modification
        'mlaA', 'mlaB', 'mlaC',  # Outer membrane lipid asymmetry
    ]

    config.membrane_genes = [
        'lpxA', 'lpxB', 'lpxC', 'lpxD',  # Lipid A biosynthesis
        'ompA', 'ompC', 'ompF',  # Outer membrane porins
        'tolC', 'acrA', 'acrB',  # Efflux system
        'mlaA', 'mlaB', 'mlaC', 'mlaD',  # Phospholipid transport
    ]

    config.regulatory_genes = [
        'marA', 'marR', 'marB',  # Mar regulon
        'soxS', 'soxR',  # SoxRS system
        'rob', 'acrR',  # Efflux regulators
        'cpxA', 'cpxR',  # Envelope stress response
        'rcsA', 'rcsB', 'rcsC',  # Rcs phosphorelay
    ]

    return config


def get_pseudomonas_aeruginosa_config() -> SpeciesConfig:
    """Get P. aeruginosa specific configuration."""
    config = SpeciesConfig(
        species="Pseudomonas aeruginosa",
        ncbi_taxid="287",
        genome_size=6300000,
        gc_content=0.66,
        chromosome_count=1,
    )

    config.resistance_genes = {
        'mexA': ResistanceGene(
            gene_name='mexA',
            product='Efflux pump membrane fusion protein',
            resistance_class='multidrug',
            mechanism='efflux',
        ),
        'mexB': ResistanceGene(
            gene_name='mexB',
            product='Multidrug efflux pump',
            resistance_class='multidrug',
            mechanism='efflux',
        ),
        'oprD': ResistanceGene(
            gene_name='oprD',
            product='Outer membrane porin D',
            resistance_class='carbapenem',
            mechanism='reduced permeability',
            notes='Loss of OprD causes imipenem resistance',
        ),
        'ampC': ResistanceGene(
            gene_name='ampC',
            product='Beta-lactamase AmpC',
            resistance_class='beta-lactam',
            mechanism='enzymatic inactivation',
        ),
        'parR': ResistanceGene(
            gene_name='parR',
            product='Response regulator ParR',
            resistance_class='AMP',
            mechanism='LPS modification',
        ),
        'parS': ResistanceGene(
            gene_name='parS',
            product='Sensor kinase ParS',
            resistance_class='AMP',
            mechanism='LPS modification',
        ),
    }

    config.amp_resistance_genes = [
        'parR', 'parS',  # ParRS system
        'pmrA', 'pmrB',  # PmrAB system
        'arnA', 'arnB', 'arnC', 'arnT',  # L-Ara4N modification
        'cprR', 'cprS',  # CprRS system
        'phoP', 'phoQ',  # PhoPQ system
    ]

    return config


# Species configuration registry
SPECIES_CONFIGS = {
    'staphylococcus_aureus': get_staphylococcus_aureus_config,
    's_aureus': get_staphylococcus_aureus_config,
    'saureus': get_staphylococcus_aureus_config,
    'escherichia_coli': get_escherichia_coli_config,
    'e_coli': get_escherichia_coli_config,
    'ecoli': get_escherichia_coli_config,
    'pseudomonas_aeruginosa': get_pseudomonas_aeruginosa_config,
    'p_aeruginosa': get_pseudomonas_aeruginosa_config,
    'paeruginosa': get_pseudomonas_aeruginosa_config,
}


def get_species_config(species: str) -> Optional[SpeciesConfig]:
    """
    Get species-specific configuration.

    Args:
        species: Species name (flexible matching)

    Returns:
        SpeciesConfig or None if not found
    """
    # Normalize species name
    species_key = species.lower().replace(' ', '_').replace('.', '')

    if species_key in SPECIES_CONFIGS:
        return SPECIES_CONFIGS[species_key]()

    # Try partial matching
    for key, config_fn in SPECIES_CONFIGS.items():
        if key in species_key or species_key in key:
            return config_fn()

    return None


def list_available_species() -> List[str]:
    """List available species configurations."""
    # Get unique species names
    seen = set()
    species_list = []

    for config_fn in SPECIES_CONFIGS.values():
        config = config_fn()
        if config.species not in seen:
            seen.add(config.species)
            species_list.append(config.species)

    return species_list
