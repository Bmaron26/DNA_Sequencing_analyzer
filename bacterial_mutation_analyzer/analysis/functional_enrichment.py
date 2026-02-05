"""
Functional enrichment analysis for bacterial mutations.
Includes COG categories, functional groups, and pathway-like analysis.
"""

import os
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set
from dataclasses import dataclass, field
from collections import defaultdict
import math

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats


# COG (Clusters of Orthologous Groups) functional categories
COG_CATEGORIES = {
    # Information storage and processing
    'J': ('Translation, ribosomal structure and biogenesis', 'Information'),
    'A': ('RNA processing and modification', 'Information'),
    'K': ('Transcription', 'Information'),
    'L': ('Replication, recombination and repair', 'Information'),
    'B': ('Chromatin structure and dynamics', 'Information'),

    # Cellular processes and signaling
    'D': ('Cell cycle control, cell division', 'Cellular'),
    'Y': ('Nuclear structure', 'Cellular'),
    'V': ('Defense mechanisms', 'Cellular'),
    'T': ('Signal transduction mechanisms', 'Cellular'),
    'M': ('Cell wall/membrane biogenesis', 'Cellular'),
    'N': ('Cell motility', 'Cellular'),
    'Z': ('Cytoskeleton', 'Cellular'),
    'W': ('Extracellular structures', 'Cellular'),
    'U': ('Intracellular trafficking and secretion', 'Cellular'),
    'O': ('Post-translational modification, protein turnover', 'Cellular'),

    # Metabolism
    'C': ('Energy production and conversion', 'Metabolism'),
    'G': ('Carbohydrate transport and metabolism', 'Metabolism'),
    'E': ('Amino acid transport and metabolism', 'Metabolism'),
    'F': ('Nucleotide transport and metabolism', 'Metabolism'),
    'H': ('Coenzyme transport and metabolism', 'Metabolism'),
    'I': ('Lipid transport and metabolism', 'Metabolism'),
    'P': ('Inorganic ion transport and metabolism', 'Metabolism'),
    'Q': ('Secondary metabolites biosynthesis', 'Metabolism'),

    # Poorly characterized
    'R': ('General function prediction only', 'Unknown'),
    'S': ('Function unknown', 'Unknown'),
    'X': ('Mobilome: prophages, transposons', 'Mobile'),
}

# S. aureus specific functional groups relevant to AMP resistance
SAUREUS_FUNCTIONAL_GROUPS = {
    'Cell membrane': [
        'membrane', 'lipid', 'phospholipid', 'fatty acid', 'mprF', 'dltA', 'dltB',
        'dltC', 'dltD', 'pgsA', 'cls', 'fmtA', 'lysC'
    ],
    'Cell wall': [
        'peptidoglycan', 'murein', 'teichoic', 'wall', 'pbp', 'murA', 'murB',
        'murC', 'murD', 'murE', 'murF', 'murG', 'ddl', 'femA', 'femB', 'femX'
    ],
    'Surface proteins': [
        'surface', 'adhesin', 'fibronectin', 'collagen', 'clumping', 'protein A',
        'spa', 'fnbA', 'fnbB', 'clfA', 'clfB', 'sdrC', 'sdrD', 'sdrE'
    ],
    'Transport': [
        'transporter', 'permease', 'ABC', 'efflux', 'pump', 'antiporter',
        'symporter', 'channel', 'import', 'export'
    ],
    'Two-component systems': [
        'histidine kinase', 'response regulator', 'sensor', 'two-component',
        'agrA', 'agrB', 'agrC', 'saeR', 'saeS', 'vraR', 'vraS', 'walR', 'walK',
        'graR', 'graS', 'arlR', 'arlS'
    ],
    'Stress response': [
        'stress', 'heat shock', 'cold shock', 'oxidative', 'osmotic', 'sigma',
        'sigB', 'rsbU', 'rsbV', 'rsbW', 'clpC', 'clpP', 'clpX', 'groEL', 'dnaK'
    ],
    'Virulence': [
        'toxin', 'hemolysin', 'leukocidin', 'enterotoxin', 'exotoxin', 'protease',
        'lipase', 'nuclease', 'hla', 'hlb', 'hlg', 'lukS', 'lukF', 'agrA', 'sarA'
    ],
    'DNA/RNA metabolism': [
        'replication', 'transcription', 'ribosom', 'polymerase', 'helicase',
        'gyrase', 'topoisomerase', 'dnaA', 'dnaB', 'dnaE', 'rpoA', 'rpoB', 'rpoC'
    ],
    'Regulation': [
        'regulator', 'repressor', 'activator', 'transcription factor',
        'sigma factor', 'anti-sigma'
    ],
    'Unknown': []
}


@dataclass
class EnrichmentResult:
    """Result of enrichment analysis for a functional category."""
    category: str
    category_name: str
    observed: int
    expected: float
    fold_enrichment: float
    p_value: float
    adjusted_p_value: float = 0.0
    genes: List[str] = field(default_factory=list)
    is_significant: bool = False


class FunctionalEnrichment:
    """
    Perform functional enrichment analysis on mutation data.
    """

    def __init__(self, annotation_file: Optional[str] = None):
        """
        Initialize enrichment analyzer.

        Args:
            annotation_file: GFF3 file with gene annotations
        """
        self.annotation_file = annotation_file
        self.gene_functions = {}  # gene -> functional groups
        self.gene_products = {}   # gene -> product description
        self.total_genes = 0

        if annotation_file:
            self._parse_annotation()

    def _parse_annotation(self):
        """Parse GFF3 annotation file to extract gene functions."""
        with open(self.annotation_file, 'r', encoding='utf-8', errors='replace') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue

                feature_type = parts[2]
                if feature_type not in ['gene', 'CDS']:
                    continue

                # Parse attributes
                attrs = {}
                for attr in parts[8].split(';'):
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        attrs[key] = value

                gene_id = attrs.get('locus_tag', attrs.get('ID', ''))
                gene_name = attrs.get('gene', attrs.get('Name', ''))
                product = attrs.get('product', '')

                if gene_id:
                    self.gene_products[gene_id] = product
                    self.gene_functions[gene_id] = self._classify_gene(
                        gene_name, product
                    )
                    self.total_genes += 1

    def _classify_gene(self, gene_name: str, product: str) -> List[str]:
        """Classify a gene into functional groups based on name and product."""
        groups = []
        search_text = f"{gene_name} {product}".lower()

        for group_name, keywords in SAUREUS_FUNCTIONAL_GROUPS.items():
            if group_name == 'Unknown':
                continue
            for keyword in keywords:
                if keyword.lower() in search_text:
                    groups.append(group_name)
                    break

        if not groups:
            groups.append('Unknown')

        return groups

    def classify_mutations(self, mutations: List[Dict]) -> Dict[str, List[Dict]]:
        """
        Classify mutations into functional groups.

        Args:
            mutations: List of mutation dictionaries

        Returns:
            Dict mapping functional group to list of mutations
        """
        classified = defaultdict(list)

        for mut in mutations:
            gene = mut.get('GENE', mut.get('gene_name', ''))
            locus_tag = mut.get('LOCUS_TAG', mut.get('locus_tag', ''))
            product = mut.get('PRODUCT', mut.get('product', ''))

            # Try to get classification from annotation
            gene_id = locus_tag or gene
            if gene_id in self.gene_functions:
                groups = self.gene_functions[gene_id]
            else:
                # Classify based on available info
                groups = self._classify_gene(gene, product)

            for group in groups:
                classified[group].append(mut)

        return dict(classified)

    def calculate_enrichment(self, mutations: List[Dict],
                            background_size: Optional[int] = None) -> List[EnrichmentResult]:
        """
        Calculate enrichment of functional categories in mutations.

        Args:
            mutations: List of mutation dictionaries
            background_size: Total number of genes in genome (for expected calculation)

        Returns:
            List of EnrichmentResult objects sorted by p-value
        """
        if background_size is None:
            background_size = self.total_genes or 2800  # Default S. aureus gene count

        # Count mutations per group
        classified = self.classify_mutations(mutations)
        total_mutations = len(mutations)

        # Count genes in background per group
        background_counts = defaultdict(int)
        for gene_id, groups in self.gene_functions.items():
            for group in groups:
                background_counts[group] += 1

        # Set default background counts if annotation not loaded
        if not background_counts:
            # Rough estimates for S. aureus
            background_counts = {
                'Cell membrane': 150,
                'Cell wall': 80,
                'Surface proteins': 50,
                'Transport': 200,
                'Two-component systems': 40,
                'Stress response': 60,
                'Virulence': 70,
                'DNA/RNA metabolism': 250,
                'Regulation': 150,
                'Unknown': 1700
            }

        results = []

        for group_name in SAUREUS_FUNCTIONAL_GROUPS.keys():
            observed = len(classified.get(group_name, []))
            background = background_counts.get(group_name, 100)

            # Expected number based on proportion in genome
            expected = total_mutations * (background / background_size)

            # Fold enrichment
            fold_enrichment = observed / expected if expected > 0 else 0

            # Hypergeometric test (Fisher's exact)
            # Using binomial approximation for simplicity
            if total_mutations > 0 and background > 0:
                p_value = stats.binom_test(
                    observed, total_mutations,
                    background / background_size,
                    alternative='greater'
                ) if observed > expected else stats.binom_test(
                    observed, total_mutations,
                    background / background_size,
                    alternative='less'
                )
            else:
                p_value = 1.0

            genes = [m.get('GENE', m.get('gene_name', ''))
                    for m in classified.get(group_name, [])]

            results.append(EnrichmentResult(
                category=group_name,
                category_name=group_name,
                observed=observed,
                expected=expected,
                fold_enrichment=fold_enrichment,
                p_value=p_value,
                genes=genes
            ))

        # Benjamini-Hochberg correction
        results.sort(key=lambda x: x.p_value)
        n = len(results)
        for i, result in enumerate(results):
            result.adjusted_p_value = min(result.p_value * n / (i + 1), 1.0)
            result.is_significant = result.adjusted_p_value < 0.05

        return results

    def plot_enrichment(self, results: List[EnrichmentResult],
                       output_path: str, title: str = '') -> str:
        """
        Create enrichment bar plot.

        Args:
            results: List of EnrichmentResult objects
            output_path: Where to save the plot
            title: Plot title

        Returns:
            Path to saved plot
        """
        fig, axes = plt.subplots(1, 2, figsize=(14, 8))

        # Sort by fold enrichment
        results_sorted = sorted(results, key=lambda x: x.fold_enrichment, reverse=True)
        results_sorted = [r for r in results_sorted if r.observed > 0]

        if not results_sorted:
            plt.close()
            return output_path

        categories = [r.category for r in results_sorted]
        fold_enrichments = [r.fold_enrichment for r in results_sorted]
        observed_counts = [r.observed for r in results_sorted]
        p_values = [r.adjusted_p_value for r in results_sorted]

        # Color by significance
        colors = ['#E41A1C' if p < 0.05 else '#377EB8' if p < 0.1 else '#999999'
                 for p in p_values]

        # Plot 1: Fold enrichment
        ax1 = axes[0]
        y_pos = range(len(categories))
        bars = ax1.barh(y_pos, fold_enrichments, color=colors)
        ax1.axvline(x=1, color='black', linestyle='--', alpha=0.5)
        ax1.set_yticks(y_pos)
        ax1.set_yticklabels(categories)
        ax1.set_xlabel('Fold Enrichment')
        ax1.set_title('Functional Category Enrichment')
        ax1.invert_yaxis()

        # Add count labels
        for i, (bar, count) in enumerate(zip(bars, observed_counts)):
            ax1.text(bar.get_width() + 0.05, bar.get_y() + bar.get_height()/2,
                    f'n={count}', va='center', fontsize=9)

        # Plot 2: Mutation counts
        ax2 = axes[1]
        ax2.barh(y_pos, observed_counts, color=colors)
        ax2.set_yticks(y_pos)
        ax2.set_yticklabels([])
        ax2.set_xlabel('Number of Mutations')
        ax2.set_title('Mutation Counts by Category')
        ax2.invert_yaxis()

        # Add legend
        legend_elements = [
            plt.Line2D([0], [0], color='#E41A1C', lw=4, label='p < 0.05'),
            plt.Line2D([0], [0], color='#377EB8', lw=4, label='p < 0.10'),
            plt.Line2D([0], [0], color='#999999', lw=4, label='Not significant')
        ]
        ax2.legend(handles=legend_elements, loc='lower right')

        if title:
            fig.suptitle(title, fontsize=14, fontweight='bold')

        plt.tight_layout()
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close()

        return output_path

    def plot_category_heatmap(self, samples_data: Dict[str, List[Dict]],
                             output_path: str, title: str = '') -> str:
        """
        Create heatmap of functional categories across samples.

        Args:
            samples_data: Dict mapping sample name to list of mutations
            output_path: Where to save the plot
            title: Plot title

        Returns:
            Path to saved plot
        """
        # Build matrix
        categories = list(SAUREUS_FUNCTIONAL_GROUPS.keys())
        samples = list(samples_data.keys())

        matrix = np.zeros((len(categories), len(samples)))

        for j, sample in enumerate(samples):
            classified = self.classify_mutations(samples_data[sample])
            for i, cat in enumerate(categories):
                matrix[i, j] = len(classified.get(cat, []))

        # Create heatmap
        fig, ax = plt.subplots(figsize=(max(12, len(samples) * 0.5), 10))

        sns.heatmap(matrix, xticklabels=samples, yticklabels=categories,
                   cmap='YlOrRd', annot=True, fmt='.0f', ax=ax,
                   cbar_kws={'label': 'Number of Mutations'})

        ax.set_xlabel('Sample')
        ax.set_ylabel('Functional Category')

        if title:
            ax.set_title(title, fontsize=14, fontweight='bold')

        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close()

        return output_path


def run_functional_analysis(results_dir: str, output_dir: str,
                           annotation_file: Optional[str] = None,
                           caller: str = 'bcftools') -> Dict[str, str]:
    """
    Run functional enrichment analysis on all samples.

    Args:
        results_dir: Directory containing sample results
        output_dir: Where to save analysis outputs
        annotation_file: GFF3 annotation file
        caller: Which variant caller results to use

    Returns:
        Dict of output file paths
    """
    os.makedirs(output_dir, exist_ok=True)

    analyzer = FunctionalEnrichment(annotation_file)

    # Load all sample mutations
    results_path = Path(results_dir)
    samples_data = {}
    all_mutations = []

    for sample_dir in results_path.iterdir():
        if not sample_dir.is_dir():
            continue

        sample_name = sample_dir.name

        # Find mutation file
        if caller == 'bcftools':
            mut_file = sample_dir / f"{sample_name}_mutations.csv"
        else:
            mut_file = sample_dir / f"{sample_name}_mutations_{caller}.csv"

        if not mut_file.exists():
            mut_file = sample_dir / f"{sample_name}_mutations.csv"
            if not mut_file.exists():
                continue

        df = pd.read_csv(mut_file)
        mutations = df.to_dict('records')
        samples_data[sample_name] = mutations
        all_mutations.extend(mutations)

    outputs = {}

    # 1. Overall enrichment analysis
    print("Running enrichment analysis...")
    results = analyzer.calculate_enrichment(all_mutations)

    # Save enrichment table
    enrichment_df = pd.DataFrame([{
        'category': r.category,
        'observed': r.observed,
        'expected': f'{r.expected:.1f}',
        'fold_enrichment': f'{r.fold_enrichment:.2f}',
        'p_value': f'{r.p_value:.4f}',
        'adjusted_p_value': f'{r.adjusted_p_value:.4f}',
        'significant': r.is_significant,
        'genes': ';'.join(set(r.genes))
    } for r in results])

    outputs['enrichment_table'] = os.path.join(output_dir,
                                               f"functional_enrichment_{caller}.csv")
    enrichment_df.to_csv(outputs['enrichment_table'], index=False)
    print(f"Created: {outputs['enrichment_table']}")

    # 2. Enrichment plot
    outputs['enrichment_plot'] = os.path.join(output_dir,
                                              f"enrichment_plot_{caller}.png")
    analyzer.plot_enrichment(results, outputs['enrichment_plot'],
                            title=f'Functional Enrichment ({caller})')
    print(f"Created: {outputs['enrichment_plot']}")

    # 3. Category heatmap across samples
    if len(samples_data) > 1:
        outputs['category_heatmap'] = os.path.join(output_dir,
                                                   f"category_heatmap_{caller}.png")
        analyzer.plot_category_heatmap(samples_data, outputs['category_heatmap'],
                                      title=f'Mutations by Functional Category ({caller})')
        print(f"Created: {outputs['category_heatmap']}")

    # 4. Per-sample classification
    classification_rows = []
    for sample_name, mutations in samples_data.items():
        classified = analyzer.classify_mutations(mutations)
        for category, muts in classified.items():
            for mut in muts:
                classification_rows.append({
                    'sample': sample_name,
                    'category': category,
                    'gene': mut.get('GENE', mut.get('gene_name', '')),
                    'locus_tag': mut.get('LOCUS_TAG', mut.get('locus_tag', '')),
                    'product': mut.get('PRODUCT', mut.get('product', '')),
                    'effect': mut.get('EFFECT', mut.get('effect', '')),
                    'aa_change': mut.get('AA_CHANGE', mut.get('amino_acid_change', ''))
                })

    if classification_rows:
        outputs['mutation_classification'] = os.path.join(
            output_dir, f"mutation_classification_{caller}.csv"
        )
        pd.DataFrame(classification_rows).to_csv(
            outputs['mutation_classification'], index=False
        )
        print(f"Created: {outputs['mutation_classification']}")

    return outputs
