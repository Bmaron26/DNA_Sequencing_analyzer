"""
AMR database integration for CARD, ResFinder, and other resistance databases.

Provides:
- Local database caching
- Query interface for known resistance mutations
- Cross-reference of variants with AMR databases
"""

import os
import json
import csv
import subprocess
from pathlib import Path
from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional, Set
from collections import defaultdict
import logging
import urllib.request
import gzip
import shutil

logger = logging.getLogger(__name__)


@dataclass
class AMRMutation:
    """Known AMR mutation from database."""
    gene: str
    mutation: str  # e.g., "S84L"
    position: int
    reference_aa: str
    mutant_aa: str
    drug_class: str
    antibiotics: List[str] = field(default_factory=list)
    resistance_level: str = ""  # high, moderate, low
    source_database: str = ""
    evidence_level: str = ""  # clinical, in_vitro, predicted
    pmid: str = ""  # PubMed reference
    notes: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return {
            'gene': self.gene,
            'mutation': self.mutation,
            'position': self.position,
            'reference_aa': self.reference_aa,
            'mutant_aa': self.mutant_aa,
            'drug_class': self.drug_class,
            'antibiotics': self.antibiotics,
            'resistance_level': self.resistance_level,
            'source_database': self.source_database,
            'evidence_level': self.evidence_level,
            'pmid': self.pmid,
            'notes': self.notes,
        }


@dataclass
class AMRHit:
    """Result of matching a variant against AMR database."""
    variant_id: str
    chromosome: str
    position: int
    gene: str
    amino_acid_change: str
    is_known_resistance: bool = False
    resistance_info: Optional[AMRMutation] = None
    is_resistance_gene: bool = False
    gene_function: str = ""
    resistance_class: str = ""
    mechanism: str = ""
    confidence: str = ""  # known, predicted, possible
    notes: str = ""

    def to_dict(self) -> Dict[str, Any]:
        result = {
            'variant_id': self.variant_id,
            'chromosome': self.chromosome,
            'position': self.position,
            'gene': self.gene,
            'amino_acid_change': self.amino_acid_change,
            'is_known_resistance': self.is_known_resistance,
            'is_resistance_gene': self.is_resistance_gene,
            'gene_function': self.gene_function,
            'resistance_class': self.resistance_class,
            'mechanism': self.mechanism,
            'confidence': self.confidence,
            'notes': self.notes,
        }
        if self.resistance_info:
            result['resistance_info'] = self.resistance_info.to_dict()
        return result


class AMRDatabase:
    """
    Interface to antimicrobial resistance databases.

    Supports:
    - CARD (Comprehensive Antibiotic Resistance Database)
    - ResFinder point mutations
    - Custom species-specific databases
    - Local caching for offline use
    """

    CARD_ONTOLOGY_URL = "https://card.mcmaster.ca/latest/ontology"
    RESFINDER_URL = "https://bitbucket.org/genomicepidemiology/resfinder_db/raw/master/"

    def __init__(self, cache_dir: Optional[str] = None):
        self.cache_dir = Path(cache_dir or os.path.expanduser("~/.bma_cache/amr"))
        self.cache_dir.mkdir(parents=True, exist_ok=True)

        # In-memory databases
        self.known_mutations: Dict[str, Dict[str, AMRMutation]] = defaultdict(dict)
        self.resistance_genes: Dict[str, Dict[str, Any]] = {}
        self.gene_drug_class: Dict[str, str] = {}

        # Load cached data if available
        self._load_cache()

    def _load_cache(self) -> None:
        """Load cached database from disk."""
        mutations_cache = self.cache_dir / "mutations.json"
        genes_cache = self.cache_dir / "genes.json"

        if mutations_cache.exists():
            try:
                with open(mutations_cache, 'r') as f:
                    data = json.load(f)
                    for gene, mutations in data.items():
                        for mut_str, mut_data in mutations.items():
                            self.known_mutations[gene][mut_str] = AMRMutation(
                                gene=mut_data['gene'],
                                mutation=mut_data['mutation'],
                                position=mut_data['position'],
                                reference_aa=mut_data['reference_aa'],
                                mutant_aa=mut_data['mutant_aa'],
                                drug_class=mut_data['drug_class'],
                                antibiotics=mut_data.get('antibiotics', []),
                                resistance_level=mut_data.get('resistance_level', ''),
                                source_database=mut_data.get('source_database', ''),
                                evidence_level=mut_data.get('evidence_level', ''),
                                pmid=mut_data.get('pmid', ''),
                                notes=mut_data.get('notes', ''),
                            )
                logger.info(f"Loaded {sum(len(m) for m in self.known_mutations.values())} cached mutations")
            except Exception as e:
                logger.warning(f"Could not load mutations cache: {e}")

        if genes_cache.exists():
            try:
                with open(genes_cache, 'r') as f:
                    self.resistance_genes = json.load(f)
                logger.info(f"Loaded {len(self.resistance_genes)} cached resistance genes")
            except Exception as e:
                logger.warning(f"Could not load genes cache: {e}")

    def _save_cache(self) -> None:
        """Save database to disk cache."""
        mutations_cache = self.cache_dir / "mutations.json"
        genes_cache = self.cache_dir / "genes.json"

        # Convert AMRMutation objects to dicts
        mutations_data = {}
        for gene, mutations in self.known_mutations.items():
            mutations_data[gene] = {
                mut_str: mut.to_dict() for mut_str, mut in mutations.items()
            }

        with open(mutations_cache, 'w') as f:
            json.dump(mutations_data, f, indent=2)

        with open(genes_cache, 'w') as f:
            json.dump(self.resistance_genes, f, indent=2)

        logger.info(f"Saved AMR database cache to {self.cache_dir}")

    def load_builtin_data(self) -> None:
        """Load built-in resistance data for common species."""
        # S. aureus mutations
        self._add_saureus_mutations()
        # E. coli mutations
        self._add_ecoli_mutations()
        # Common AMP resistance mutations
        self._add_amp_mutations()

        self._save_cache()

    def _add_saureus_mutations(self) -> None:
        """Add known S. aureus resistance mutations."""
        mutations = [
            # Fluoroquinolone resistance
            AMRMutation('gyrA', 'S84L', 84, 'S', 'L', 'fluoroquinolone',
                       ['ciprofloxacin', 'levofloxacin'], 'high', 'literature'),
            AMRMutation('gyrA', 'S84A', 84, 'S', 'A', 'fluoroquinolone',
                       ['ciprofloxacin'], 'high', 'literature'),
            AMRMutation('gyrA', 'E88K', 88, 'E', 'K', 'fluoroquinolone',
                       ['ciprofloxacin', 'levofloxacin'], 'moderate', 'literature'),
            AMRMutation('grlA', 'S80F', 80, 'S', 'F', 'fluoroquinolone',
                       ['ciprofloxacin'], 'high', 'literature'),
            AMRMutation('grlA', 'S80Y', 80, 'S', 'Y', 'fluoroquinolone',
                       ['ciprofloxacin'], 'moderate', 'literature'),

            # Rifampicin resistance
            AMRMutation('rpoB', 'H481Y', 481, 'H', 'Y', 'rifampicin',
                       ['rifampicin'], 'high', 'literature'),
            AMRMutation('rpoB', 'H481N', 481, 'H', 'N', 'rifampicin',
                       ['rifampicin'], 'high', 'literature'),
            AMRMutation('rpoB', 'S486L', 486, 'S', 'L', 'rifampicin',
                       ['rifampicin'], 'high', 'literature'),

            # Daptomycin/AMP resistance (mprF)
            AMRMutation('mprF', 'S295L', 295, 'S', 'L', 'AMP',
                       ['daptomycin', 'cationic AMPs'], 'high', 'literature',
                       notes='Gain-of-function, increased lysylPG'),
            AMRMutation('mprF', 'L341F', 341, 'L', 'F', 'AMP',
                       ['daptomycin', 'cationic AMPs'], 'high', 'literature'),
            AMRMutation('mprF', 'T345A', 345, 'T', 'A', 'AMP',
                       ['daptomycin', 'cationic AMPs'], 'moderate', 'literature'),
            AMRMutation('mprF', 'S337L', 337, 'S', 'L', 'AMP',
                       ['daptomycin'], 'high', 'literature'),

            # Vancomycin intermediate resistance
            AMRMutation('vraS', 'I5N', 5, 'I', 'N', 'glycopeptide',
                       ['vancomycin'], 'moderate', 'literature',
                       notes='VISA-associated'),
            AMRMutation('walK', 'G223D', 223, 'G', 'D', 'glycopeptide',
                       ['vancomycin'], 'moderate', 'literature'),

            # Mupirocin resistance
            AMRMutation('ileS', 'V588F', 588, 'V', 'F', 'mupirocin',
                       ['mupirocin'], 'high', 'literature'),

            # Linezolid resistance
            AMRMutation('rplC', 'G152D', 152, 'G', 'D', 'oxazolidinone',
                       ['linezolid'], 'high', 'literature'),
        ]

        for mut in mutations:
            self.known_mutations[mut.gene][mut.mutation] = mut

        # Add resistance genes
        self.resistance_genes.update({
            'mecA': {'class': 'beta-lactam', 'mechanism': 'target alteration'},
            'vanA': {'class': 'glycopeptide', 'mechanism': 'target alteration'},
            'vanB': {'class': 'glycopeptide', 'mechanism': 'target alteration'},
            'ermA': {'class': 'macrolide', 'mechanism': 'target modification'},
            'ermB': {'class': 'macrolide', 'mechanism': 'target modification'},
            'ermC': {'class': 'macrolide', 'mechanism': 'target modification'},
            'tetK': {'class': 'tetracycline', 'mechanism': 'efflux'},
            'tetM': {'class': 'tetracycline', 'mechanism': 'ribosomal protection'},
            'dfrA': {'class': 'trimethoprim', 'mechanism': 'target bypass'},
            'fusB': {'class': 'fusidic acid', 'mechanism': 'ribosomal protection'},
            'norA': {'class': 'fluoroquinolone', 'mechanism': 'efflux'},
            'mprF': {'class': 'AMP', 'mechanism': 'membrane modification'},
            'dltA': {'class': 'AMP', 'mechanism': 'membrane modification'},
        })

    def _add_ecoli_mutations(self) -> None:
        """Add known E. coli resistance mutations."""
        mutations = [
            # Fluoroquinolone resistance
            AMRMutation('gyrA', 'S83L', 83, 'S', 'L', 'fluoroquinolone',
                       ['ciprofloxacin', 'levofloxacin'], 'high', 'literature'),
            AMRMutation('gyrA', 'D87N', 87, 'D', 'N', 'fluoroquinolone',
                       ['ciprofloxacin'], 'high', 'literature'),
            AMRMutation('gyrA', 'D87G', 87, 'D', 'G', 'fluoroquinolone',
                       ['ciprofloxacin'], 'moderate', 'literature'),
            AMRMutation('parC', 'S80I', 80, 'S', 'I', 'fluoroquinolone',
                       ['ciprofloxacin'], 'moderate', 'literature'),

            # Colistin/AMP resistance
            AMRMutation('pmrA', 'L14R', 14, 'L', 'R', 'AMP',
                       ['colistin', 'polymyxin B'], 'high', 'literature'),
            AMRMutation('pmrB', 'L10R', 10, 'L', 'R', 'AMP',
                       ['colistin', 'polymyxin B'], 'high', 'literature'),
            AMRMutation('phoP', 'R81S', 81, 'R', 'S', 'AMP',
                       ['colistin'], 'moderate', 'literature'),

            # Rifampicin resistance
            AMRMutation('rpoB', 'S531L', 531, 'S', 'L', 'rifampicin',
                       ['rifampicin'], 'high', 'literature'),
            AMRMutation('rpoB', 'H526Y', 526, 'H', 'Y', 'rifampicin',
                       ['rifampicin'], 'high', 'literature'),
        ]

        for mut in mutations:
            self.known_mutations[mut.gene][mut.mutation] = mut

        # Add resistance genes
        self.resistance_genes.update({
            'blaTEM': {'class': 'beta-lactam', 'mechanism': 'enzymatic inactivation'},
            'blaSHV': {'class': 'beta-lactam', 'mechanism': 'enzymatic inactivation'},
            'blaCTX-M': {'class': 'beta-lactam', 'mechanism': 'enzymatic inactivation'},
            'blaKPC': {'class': 'carbapenem', 'mechanism': 'enzymatic inactivation'},
            'blaNDM': {'class': 'carbapenem', 'mechanism': 'enzymatic inactivation'},
            'mcr-1': {'class': 'AMP', 'mechanism': 'LPS modification'},
            'qnrA': {'class': 'fluoroquinolone', 'mechanism': 'target protection'},
            'qnrB': {'class': 'fluoroquinolone', 'mechanism': 'target protection'},
            'aac(6\')': {'class': 'aminoglycoside', 'mechanism': 'enzymatic modification'},
            'aph(3\')': {'class': 'aminoglycoside', 'mechanism': 'enzymatic modification'},
        })

    def _add_amp_mutations(self) -> None:
        """Add AMP-specific resistance mutations across species."""
        # These are general AMP resistance mutations
        amp_mutations = [
            # Common membrane modification genes
            AMRMutation('pgsA', 'any', 0, '', '', 'AMP',
                       ['cationic AMPs'], 'variable', 'literature',
                       notes='PG synthase - loss may affect membrane'),
            AMRMutation('cls', 'any', 0, '', '', 'AMP',
                       ['cationic AMPs'], 'variable', 'literature',
                       notes='Cardiolipin synthase'),
        ]

        # Don't add these to known_mutations since they're not specific
        # Just note that these genes are relevant

    def query_mutation(self, gene: str, mutation: str) -> Optional[AMRMutation]:
        """
        Query database for a specific mutation.

        Args:
            gene: Gene name
            mutation: Mutation string (e.g., "S84L")

        Returns:
            AMRMutation if found, None otherwise
        """
        gene_mutations = self.known_mutations.get(gene, {})
        return gene_mutations.get(mutation)

    def query_gene(self, gene: str) -> Optional[Dict[str, Any]]:
        """
        Query database for resistance gene information.

        Args:
            gene: Gene name

        Returns:
            Gene info dict if found, None otherwise
        """
        # Check exact match
        if gene in self.resistance_genes:
            return self.resistance_genes[gene]

        # Check case-insensitive
        gene_lower = gene.lower()
        for name, info in self.resistance_genes.items():
            if name.lower() == gene_lower:
                return info

        return None

    def is_resistance_gene(self, gene: str) -> bool:
        """Check if gene is a known resistance gene."""
        return self.query_gene(gene) is not None

    def get_gene_drug_class(self, gene: str) -> Optional[str]:
        """Get drug class for a resistance gene."""
        info = self.query_gene(gene)
        if info:
            return info.get('class')
        return None

    def annotate_variant(self, variant: Dict[str, Any],
                        species_config=None) -> AMRHit:
        """
        Annotate a variant with AMR information.

        Args:
            variant: Variant dictionary
            species_config: Optional species configuration

        Returns:
            AMRHit with resistance annotation
        """
        gene = variant.get('gene_name') or variant.get('locus_tag') or ''
        aa_change = variant.get('amino_acid_change', '')

        hit = AMRHit(
            variant_id=f"{variant.get('chromosome')}:{variant.get('position')}",
            chromosome=variant.get('chromosome', ''),
            position=variant.get('position', 0),
            gene=gene,
            amino_acid_change=aa_change,
        )

        # Check if known resistance gene
        gene_info = self.query_gene(gene)
        if gene_info:
            hit.is_resistance_gene = True
            hit.resistance_class = gene_info.get('class', '')
            hit.mechanism = gene_info.get('mechanism', '')
            hit.gene_function = f"Resistance gene: {hit.resistance_class}"

        # Check species-specific config
        if species_config:
            if species_config.is_resistance_gene(gene):
                hit.is_resistance_gene = True
                gene_detail = species_config.get_gene_info(gene)
                if gene_detail:
                    hit.resistance_class = gene_detail.resistance_class
                    hit.mechanism = gene_detail.mechanism

            if species_config.is_membrane_related(gene):
                hit.gene_function = hit.gene_function or "Membrane-related"
                hit.notes = "Membrane gene - may affect AMP susceptibility"

        # Check for known mutation
        if aa_change:
            # Parse amino acid change (e.g., "S84L" or "Ser84Leu")
            mutation_str = self._normalize_aa_change(aa_change)
            known_mut = self.query_mutation(gene, mutation_str)

            if known_mut:
                hit.is_known_resistance = True
                hit.resistance_info = known_mut
                hit.confidence = 'known'
                hit.resistance_class = known_mut.drug_class
                hit.notes = known_mut.notes
            elif hit.is_resistance_gene:
                # Mutation in resistance gene but not known mutation
                hit.confidence = 'possible'
                hit.notes = f"Novel mutation in {hit.resistance_class} resistance gene"

        return hit

    def _normalize_aa_change(self, aa_change: str) -> str:
        """Normalize amino acid change format (e.g., Ser84Leu -> S84L)."""
        # Three-letter to one-letter
        aa_map = {
            'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D',
            'Cys': 'C', 'Gln': 'Q', 'Glu': 'E', 'Gly': 'G',
            'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K',
            'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S',
            'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
            'Ter': '*', 'Stop': '*',
        }

        result = aa_change
        for three, one in aa_map.items():
            result = result.replace(three, one)

        # Remove common prefixes
        result = result.replace('p.', '')

        return result

    def batch_annotate(self, variants: List[Dict[str, Any]],
                      species_config=None) -> List[AMRHit]:
        """
        Annotate multiple variants.

        Args:
            variants: List of variant dictionaries
            species_config: Optional species configuration

        Returns:
            List of AMRHit annotations
        """
        return [self.annotate_variant(v, species_config) for v in variants]

    def get_summary(self, hits: List[AMRHit]) -> Dict[str, Any]:
        """Generate summary of AMR hits."""
        summary = {
            'total_variants': len(hits),
            'known_resistance_mutations': 0,
            'resistance_gene_mutations': 0,
            'by_drug_class': defaultdict(int),
            'by_mechanism': defaultdict(int),
            'known_mutations': [],
            'possible_resistance': [],
        }

        for hit in hits:
            if hit.is_known_resistance:
                summary['known_resistance_mutations'] += 1
                summary['known_mutations'].append(hit.to_dict())
            if hit.is_resistance_gene:
                summary['resistance_gene_mutations'] += 1
                if not hit.is_known_resistance:
                    summary['possible_resistance'].append(hit.to_dict())

            if hit.resistance_class:
                summary['by_drug_class'][hit.resistance_class] += 1
            if hit.mechanism:
                summary['by_mechanism'][hit.mechanism] += 1

        summary['by_drug_class'] = dict(summary['by_drug_class'])
        summary['by_mechanism'] = dict(summary['by_mechanism'])

        return summary
