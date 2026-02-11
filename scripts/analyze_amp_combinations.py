#!/usr/bin/env python3
"""
Comprehensive analysis of AMP resistance evolution mutations.

Analyzes:
1. Single AMP similarity (which AMPs target similar genes?)
2. Combination vs Parents (novel, parent1-like, parent2-like, shared)
3. Similarity networks (singles + combinations)
4. Convergent evolution analysis
5. Functional category enrichment
6. Replicate consistency score
7. Parent dominance index
8. Venn diagrams for combinations
9. Phylogenetic/clustering tree of treatments

Usage:
    python analyze_amp_combinations.py <singles_dir> <combinations_dir> <output_dir>
"""

import sys
import os
from pathlib import Path
from typing import Dict, List, Tuple, Set, Optional
from collections import defaultdict
import re

import pandas as pd
import numpy as np

# Standard color map for AMP treatments
EVO_COLOR_MAP = {
    'Mel': '#729ECEFF',
    'Pex': '#FF9E4AFF',
    'BmKn': '#67BF5CFF',
    'Puro': '#ED665DFF',
    'Pleu': '#AD8BC9FF',
    'Smp': '#A8786EFF',
    'NPSA': '#999999'
}

# Combination color map (blend of parents)
COMBO_COLOR_MAP = {
    'Mel_Pex': '#9C9F8CFF',
    'Mel_BmKn': '#6DA78DFF',
    'Mel_Pleu': '#8E92CFFF',
    'Mel_Puro': '#8186A6FF',
    'Mel_Smp': '#8D8B89FF',
    'Pex_BmKn': '#B3AA54FF',
    'Pex_Pleu': '#D692B4FF',
    'Pex_Puro': '#F69F54FF',
    'Pex_Smp': '#CCA05CFF',
    'BmKn_Pleu': '#8DB393FF',
    'BmKn_Puro': '#8D9A5DFF',
    'BmKn_Smp': '#8FA866FF',
    'Puro_Pleu': '#DD7193FF',
    'Puro_Smp': '#CB6F66FF',
    'Smp_Pleu': '#AA80ACFF',
}


def is_synonymous(effect: str) -> bool:
    """Check if a mutation effect is synonymous."""
    if pd.isna(effect):
        return False
    effect_lower = str(effect).lower()
    return any(term in effect_lower for term in ['synonymous', 'silent', 'syn_coding'])


def extract_treatment_from_sample(sample_name: str) -> str:
    """Extract treatment from sample name (e.g., 'Mel_Pex1' -> 'Mel_Pex')."""
    if '_' in sample_name:
        # Could be combination like 'Mel_Pex1' or 'BmKn_Pleu3'
        parts = sample_name.rsplit('_', 1)
        if len(parts) == 2:
            # Check if last part ends with digit
            if parts[1] and parts[1][-1].isdigit():
                # Remove trailing digit
                return re.sub(r'\d+$', '', sample_name)
    # Single like 'Mel1' -> 'Mel'
    return re.sub(r'\d+$', '', sample_name)


def get_parent_treatments(combo_name: str) -> Tuple[str, str]:
    """Get parent treatment names from combination (e.g., 'Mel_Pex' -> ('Mel', 'Pex'))."""
    parts = combo_name.split('_')
    if len(parts) == 2:
        return parts[0], parts[1]
    return combo_name, combo_name


def dice_similarity(set1: set, set2: set) -> float:
    """Calculate Dice similarity coefficient."""
    if len(set1) == 0 and len(set2) == 0:
        return 1.0
    if len(set1) == 0 or len(set2) == 0:
        return 0.0
    intersection = len(set1 & set2)
    return (2 * intersection) / (len(set1) + len(set2))


def jaccard_similarity(set1: set, set2: set) -> float:
    """Calculate Jaccard similarity coefficient."""
    if len(set1) == 0 and len(set2) == 0:
        return 1.0
    if len(set1) == 0 or len(set2) == 0:
        return 0.0
    intersection = len(set1 & set2)
    union = len(set1 | set2)
    return intersection / union


def load_mutations_from_directory(results_dir: Path) -> Dict[str, Set[str]]:
    """
    Load mutations from a results directory.
    Returns dict mapping sample_name -> set of mutated genes.
    """
    sample_genes = {}

    for sample_dir in results_dir.iterdir():
        if not sample_dir.is_dir():
            continue

        sample_name = sample_dir.name

        # Try freebayes results first
        mut_file = sample_dir / f"{sample_name}_mutations_freebayes.csv"
        if not mut_file.exists():
            mut_file = sample_dir / f"{sample_name}_mutations.csv"
        if not mut_file.exists():
            continue

        df = pd.read_csv(mut_file)

        # Detect gene column
        gene_col = None
        for col in ['gene_name', 'GENE', 'gene', 'locus_tag']:
            if col in df.columns:
                gene_col = col
                break

        # Detect effect column
        effect_col = None
        for col in ['effect', 'EFFECT']:
            if col in df.columns:
                effect_col = col
                break

        genes = set()
        for _, row in df.iterrows():
            # Skip synonymous
            if effect_col and is_synonymous(row.get(effect_col, '')):
                continue

            if gene_col:
                gene = row.get(gene_col, '')
                if not pd.isna(gene) and str(gene).strip():
                    genes.add(str(gene).strip())

        sample_genes[sample_name] = genes

    return sample_genes


def aggregate_by_treatment(sample_genes: Dict[str, Set[str]]) -> Dict[str, Set[str]]:
    """Aggregate genes across replicates for each treatment."""
    treatment_genes = defaultdict(set)

    for sample, genes in sample_genes.items():
        treatment = extract_treatment_from_sample(sample)
        treatment_genes[treatment].update(genes)

    return dict(treatment_genes)


def get_treatment_gene_counts(sample_genes: Dict[str, Set[str]]) -> Dict[str, Dict[str, int]]:
    """Count how many replicates have each gene mutated per treatment."""
    treatment_gene_counts = defaultdict(lambda: defaultdict(int))

    for sample, genes in sample_genes.items():
        treatment = extract_treatment_from_sample(sample)
        for gene in genes:
            treatment_gene_counts[treatment][gene] += 1

    return dict(treatment_gene_counts)


def analyze_combination_vs_parents(
    combo_genes: Set[str],
    parent1_genes: Set[str],
    parent2_genes: Set[str]
) -> Dict[str, Set[str]]:
    """
    Classify mutations in combination relative to parents.

    Returns dict with:
    - 'parent1_only': genes from parent1 but not parent2
    - 'parent2_only': genes from parent2 but not parent1
    - 'shared_parents': genes in both parents
    - 'novel': genes not in either parent
    """
    in_p1 = combo_genes & parent1_genes
    in_p2 = combo_genes & parent2_genes
    in_both = combo_genes & parent1_genes & parent2_genes
    novel = combo_genes - parent1_genes - parent2_genes

    return {
        'parent1_only': in_p1 - parent2_genes,
        'parent2_only': in_p2 - parent1_genes,
        'shared_parents': in_both,
        'novel': novel,
        'total': combo_genes
    }


def create_similarity_heatmap(
    treatment_genes: Dict[str, Set[str]],
    output_path: Path,
    title: str = "Treatment Mutation Similarity"
):
    """Create heatmap of pairwise treatment similarities."""
    import matplotlib.pyplot as plt
    import seaborn as sns

    treatments = sorted(treatment_genes.keys())
    n = len(treatments)
    sim_matrix = np.zeros((n, n))

    for i, t1 in enumerate(treatments):
        for j, t2 in enumerate(treatments):
            sim_matrix[i, j] = dice_similarity(treatment_genes[t1], treatment_genes[t2])

    df_sim = pd.DataFrame(sim_matrix, index=treatments, columns=treatments)

    # Create color mapping for treatments
    colors = []
    for t in treatments:
        if t in EVO_COLOR_MAP:
            colors.append(EVO_COLOR_MAP[t])
        elif t in COMBO_COLOR_MAP:
            colors.append(COMBO_COLOR_MAP[t])
        else:
            colors.append('#888888')

    fig, ax = plt.subplots(figsize=(12, 10))

    sns.heatmap(
        df_sim,
        annot=True,
        fmt='.2f',
        cmap='RdYlBu_r',
        vmin=0, vmax=1,
        square=True,
        linewidths=1,
        cbar_kws={'label': 'Dice Similarity'},
        ax=ax,
        annot_kws={'fontsize': 9}
    )

    ax.set_title(title, fontsize=14, fontweight='bold', pad=15)

    # Color tick labels
    for i, label in enumerate(ax.get_xticklabels()):
        label.set_color(colors[i])
        label.set_fontweight('bold')
    for i, label in enumerate(ax.get_yticklabels()):
        label.set_color(colors[i])
        label.set_fontweight('bold')

    plt.setp(ax.get_xticklabels(), rotation=45, ha='right', fontsize=10)
    plt.setp(ax.get_yticklabels(), rotation=0, fontsize=10)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()

    return df_sim


def create_combination_breakdown_plot(
    combo_analysis: Dict[str, Dict[str, Set[str]]],
    output_path: Path
):
    """Create stacked bar chart showing combination mutation breakdown."""
    import matplotlib.pyplot as plt

    combos = sorted(combo_analysis.keys())

    data = {
        'Parent 1 only': [],
        'Parent 2 only': [],
        'Shared (both parents)': [],
        'Novel': []
    }

    for combo in combos:
        analysis = combo_analysis[combo]
        total = len(analysis['total']) if len(analysis['total']) > 0 else 1
        data['Parent 1 only'].append(len(analysis['parent1_only']) / total * 100)
        data['Parent 2 only'].append(len(analysis['parent2_only']) / total * 100)
        data['Shared (both parents)'].append(len(analysis['shared_parents']) / total * 100)
        data['Novel'].append(len(analysis['novel']) / total * 100)

    fig, ax = plt.subplots(figsize=(14, 8))

    x = np.arange(len(combos))
    width = 0.6

    colors = ['#729ECEFF', '#FF9E4AFF', '#67BF5CFF', '#ED665DFF']

    bottom = np.zeros(len(combos))
    for i, (label, values) in enumerate(data.items()):
        ax.bar(x, values, width, label=label, bottom=bottom, color=colors[i])
        bottom += np.array(values)

    ax.set_ylabel('Percentage of Mutated Genes', fontsize=12)
    ax.set_xlabel('Combination Treatment', fontsize=12)
    ax.set_title('Mutation Origin in Combination Treatments\n(Compared to Parent Single Treatments)',
                 fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(combos, rotation=45, ha='right', fontsize=10)
    ax.legend(loc='upper right', fontsize=10)
    ax.set_ylim(0, 105)

    # Add counts on top
    for i, combo in enumerate(combos):
        total = len(combo_analysis[combo]['total'])
        ax.text(i, 102, f'n={total}', ha='center', fontsize=9)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()


def create_similarity_network(
    treatment_genes: Dict[str, Set[str]],
    output_path: Path,
    min_edge_weight: float = 0.1,
    title: str = "Treatment Similarity Network"
):
    """Create network visualization of treatment similarities."""
    import matplotlib.pyplot as plt
    import networkx as nx

    treatments = list(treatment_genes.keys())
    G = nx.Graph()

    # Add nodes
    for t in treatments:
        is_combo = '_' in t and not t.startswith('NPSA')
        G.add_node(t, is_combo=is_combo, n_genes=len(treatment_genes[t]))

    # Add edges
    for i, t1 in enumerate(treatments):
        for j, t2 in enumerate(treatments):
            if i < j:
                sim = dice_similarity(treatment_genes[t1], treatment_genes[t2])
                if sim >= min_edge_weight:
                    G.add_edge(t1, t2, weight=sim, distance=1.0 - sim + 0.1)

    # Layout
    if G.edges():
        pos = nx.spring_layout(G, k=2, iterations=100, seed=42, weight='distance')
    else:
        pos = nx.spring_layout(G, k=2, iterations=50, seed=42)

    fig, ax = plt.subplots(figsize=(16, 14))

    # Node colors and sizes
    node_colors = []
    node_sizes = []
    for n in G.nodes():
        if n in EVO_COLOR_MAP:
            node_colors.append(EVO_COLOR_MAP[n])
        elif n in COMBO_COLOR_MAP:
            node_colors.append(COMBO_COLOR_MAP[n])
        else:
            node_colors.append('#888888')

        # Size based on number of genes
        n_genes = G.nodes[n]['n_genes']
        node_sizes.append(500 + n_genes * 30)

    # Draw nodes
    singles = [n for n in G.nodes() if not G.nodes[n]['is_combo']]
    combos = [n for n in G.nodes() if G.nodes[n]['is_combo']]

    # Draw singles as circles
    nx.draw_networkx_nodes(
        G, pos,
        nodelist=singles,
        node_size=[node_sizes[list(G.nodes()).index(n)] for n in singles],
        node_color=[node_colors[list(G.nodes()).index(n)] for n in singles],
        edgecolors='black',
        linewidths=2,
        ax=ax
    )

    # Draw combos as squares (using scatter)
    if combos:
        combo_pos = np.array([pos[n] for n in combos])
        combo_colors = [node_colors[list(G.nodes()).index(n)] for n in combos]
        combo_sizes = [node_sizes[list(G.nodes()).index(n)] for n in combos]
        ax.scatter(combo_pos[:, 0], combo_pos[:, 1],
                   s=combo_sizes, c=combo_colors, marker='s',
                   edgecolors='black', linewidths=2, zorder=3)

    # Draw edges
    edges = G.edges(data=True)
    if edges:
        for (u, v, d) in edges:
            weight = d['weight']
            width = 0.5 + 5 * weight
            alpha = 0.2 + 0.6 * weight

            if weight > 0.5:
                color = '#2E7D32'
            elif weight > 0.3:
                color = '#FF8F00'
            else:
                color = '#BDBDBD'

            nx.draw_networkx_edges(
                G, pos,
                edgelist=[(u, v)],
                width=width,
                alpha=alpha,
                edge_color=color,
                ax=ax
            )

    # Labels
    nx.draw_networkx_labels(G, pos, font_size=10, font_weight='bold', ax=ax)

    ax.set_title(title, fontsize=16, fontweight='bold', pad=20)

    # Legend
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D

    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='gray',
               markeredgecolor='black', markersize=12, label='Single AMP (circle)'),
        Line2D([0], [0], marker='s', color='w', markerfacecolor='gray',
               markeredgecolor='black', markersize=12, label='Combination (square)'),
        Line2D([0], [0], color='#2E7D32', linewidth=4, label='High similarity (>0.5)'),
        Line2D([0], [0], color='#FF8F00', linewidth=3, label='Medium (0.3-0.5)'),
        Line2D([0], [0], color='#BDBDBD', linewidth=2, label='Low (<0.3)')
    ]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=10)

    ax.text(0.02, 0.02,
            'Node size = number of mutated genes\n'
            'Edge thickness = Dice similarity\n'
            'Circles = single AMPs, Squares = combinations',
            transform=ax.transAxes, fontsize=9, verticalalignment='bottom',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    ax.axis('off')
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()


def find_convergent_genes(treatment_gene_counts: Dict[str, Dict[str, int]],
                          min_treatments: int = 3) -> pd.DataFrame:
    """Find genes mutated across multiple treatments."""
    gene_treatments = defaultdict(list)

    for treatment, gene_counts in treatment_gene_counts.items():
        for gene, count in gene_counts.items():
            gene_treatments[gene].append((treatment, count))

    convergent = []
    for gene, treatment_list in gene_treatments.items():
        if len(treatment_list) >= min_treatments:
            convergent.append({
                'gene': gene,
                'n_treatments': len(treatment_list),
                'treatments': ', '.join([t for t, _ in treatment_list]),
                'total_occurrences': sum([c for _, c in treatment_list])
            })

    df = pd.DataFrame(convergent)
    if not df.empty:
        df = df.sort_values('n_treatments', ascending=False)
    return df


def create_convergent_genes_heatmap(
    treatment_gene_counts: Dict[str, Dict[str, int]],
    output_path: Path,
    top_n: int = 30
):
    """Create heatmap showing top convergent genes across treatments."""
    import matplotlib.pyplot as plt
    import seaborn as sns

    # Find top convergent genes
    gene_treatment_count = defaultdict(int)
    for treatment, gene_counts in treatment_gene_counts.items():
        for gene in gene_counts:
            gene_treatment_count[gene] += 1

    top_genes = sorted(gene_treatment_count.items(), key=lambda x: -x[1])[:top_n]
    top_gene_names = [g for g, _ in top_genes]

    treatments = sorted(treatment_gene_counts.keys())

    # Build matrix
    matrix = np.zeros((len(top_gene_names), len(treatments)))
    for i, gene in enumerate(top_gene_names):
        for j, treatment in enumerate(treatments):
            matrix[i, j] = treatment_gene_counts[treatment].get(gene, 0)

    df_matrix = pd.DataFrame(matrix, index=top_gene_names, columns=treatments)

    fig, ax = plt.subplots(figsize=(14, 12))

    sns.heatmap(
        df_matrix,
        cmap='YlOrRd',
        annot=True,
        fmt='.0f',
        linewidths=0.5,
        cbar_kws={'label': 'Replicate count'},
        ax=ax
    )

    ax.set_title(f'Top {top_n} Convergent Genes Across Treatments\n(Number of replicates with mutation)',
                 fontsize=14, fontweight='bold', pad=15)
    ax.set_xlabel('Treatment', fontsize=12)
    ax.set_ylabel('Gene', fontsize=12)

    plt.setp(ax.get_xticklabels(), rotation=45, ha='right', fontsize=10)
    plt.setp(ax.get_yticklabels(), rotation=0, fontsize=9)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()

    return df_matrix


# =============================================================================
# ADDITIONAL PUBLICATION-WORTHY ANALYSES
# =============================================================================

def load_mutations_with_details(results_dir: Path) -> Dict[str, pd.DataFrame]:
    """
    Load full mutation data from a results directory.
    Returns dict mapping sample_name -> DataFrame with all mutation details.
    """
    sample_mutations = {}

    for sample_dir in results_dir.iterdir():
        if not sample_dir.is_dir():
            continue

        sample_name = sample_dir.name

        # Try freebayes results first
        mut_file = sample_dir / f"{sample_name}_mutations_freebayes.csv"
        if not mut_file.exists():
            mut_file = sample_dir / f"{sample_name}_mutations.csv"
        if not mut_file.exists():
            continue

        df = pd.read_csv(mut_file)
        sample_mutations[sample_name] = df

    return sample_mutations


def categorize_gene_function(gene_name: str, product: str = "") -> str:
    """
    Categorize gene by functional category based on gene name and product.
    S. aureus focused categories.
    """
    if pd.isna(gene_name):
        gene_name = ""
    if pd.isna(product):
        product = ""

    gene_lower = str(gene_name).lower()
    product_lower = str(product).lower()
    combined = gene_lower + " " + product_lower

    # Membrane/Cell wall related
    membrane_terms = ['mpr', 'mre', 'lta', 'dlt', 'pbp', 'mur', 'fts', 'rod', 'cls',
                      'pgs', 'membrane', 'cell wall', 'peptidoglycan', 'lipoteichoic',
                      'transporter', 'permease', 'efflux', 'atpa', 'atpb', 'atpc', 'atpd',
                      'atpe', 'atpf', 'atpg', 'atph', 'atp synthase']

    # Transcription/Regulation
    regulation_terms = ['agr', 'sar', 'sig', 'rpo', 'graa', 'grab', 'wal', 'vra',
                        'regulator', 'transcription', 'sigma', 'two-component',
                        'sensor', 'kinase', 'response regulator']

    # Metabolism
    metabolism_terms = ['met', 'ald', 'pyk', 'pfk', 'eno', 'gap', 'ldh', 'pfl',
                        'dehydrogenase', 'kinase', 'synthase', 'metabolism',
                        'glycolysis', 'tca', 'oxidoreductase', 'reductase']

    # DNA/RNA
    dna_rna_terms = ['dna', 'rna', 'gyr', 'top', 'pol', 'lig', 'rec', 'mut',
                    'ribosom', 'rrna', 'trna', 'helicase', 'replication',
                    'transcription', 'translation', 'gyrase', 'topoisomerase']

    # Stress response
    stress_terms = ['clp', 'gro', 'dna', 'hsp', 'csp', 'htr', 'stress', 'heat shock',
                   'cold shock', 'chaperone', 'protease', 'oxidative']

    # Virulence
    virulence_terms = ['spa', 'fnb', 'clf', 'sdr', 'cap', 'hla', 'hlb', 'hld',
                      'luk', 'sec', 'sed', 'tst', 'virulence', 'toxin', 'adhesin',
                      'protein a', 'coagulase', 'hemolysin']

    # Check categories in order of specificity
    if any(term in combined for term in virulence_terms):
        return 'Virulence'
    if any(term in combined for term in membrane_terms):
        return 'Membrane/Cell Wall'
    if any(term in combined for term in regulation_terms):
        return 'Regulation'
    if any(term in combined for term in dna_rna_terms):
        return 'DNA/RNA'
    if any(term in combined for term in stress_terms):
        return 'Stress Response'
    if any(term in combined for term in metabolism_terms):
        return 'Metabolism'

    return 'Other/Unknown'


def analyze_functional_enrichment(
    sample_mutations: Dict[str, pd.DataFrame],
    output_dir: Path
) -> pd.DataFrame:
    """
    Analyze functional category enrichment across treatments.
    Creates pie charts and stacked bar charts.
    """
    import matplotlib.pyplot as plt

    # Aggregate by treatment
    treatment_categories = defaultdict(lambda: defaultdict(int))
    treatment_genes = defaultdict(set)

    for sample, df in sample_mutations.items():
        treatment = extract_treatment_from_sample(sample)

        # Detect columns
        gene_col = None
        for col in ['gene_name', 'GENE', 'gene', 'locus_tag']:
            if col in df.columns:
                gene_col = col
                break

        product_col = None
        for col in ['product', 'PRODUCT', 'annotation']:
            if col in df.columns:
                product_col = col
                break

        effect_col = None
        for col in ['effect', 'EFFECT']:
            if col in df.columns:
                effect_col = col
                break

        for _, row in df.iterrows():
            # Skip synonymous
            if effect_col and is_synonymous(row.get(effect_col, '')):
                continue

            if gene_col:
                gene = row.get(gene_col, '')
                if pd.isna(gene) or not str(gene).strip():
                    continue

                gene = str(gene).strip()
                if gene in treatment_genes[treatment]:
                    continue  # Already counted this gene

                treatment_genes[treatment].add(gene)
                product = row.get(product_col, '') if product_col else ''
                category = categorize_gene_function(gene, product)
                treatment_categories[treatment][category] += 1

    # Create DataFrame
    categories = ['Membrane/Cell Wall', 'Regulation', 'Metabolism', 'DNA/RNA',
                  'Stress Response', 'Virulence', 'Other/Unknown']

    data = []
    for treatment in sorted(treatment_categories.keys()):
        row = {'treatment': treatment}
        total = sum(treatment_categories[treatment].values())
        for cat in categories:
            count = treatment_categories[treatment].get(cat, 0)
            row[cat] = count
            row[f'{cat}_pct'] = count / total * 100 if total > 0 else 0
        row['total'] = total
        data.append(row)

    df_enrichment = pd.DataFrame(data)
    df_enrichment.to_csv(output_dir / 'functional_enrichment.csv', index=False)

    # Create stacked bar chart
    fig, ax = plt.subplots(figsize=(14, 8))

    treatments = df_enrichment['treatment'].tolist()
    x = np.arange(len(treatments))
    width = 0.7

    # Color palette for categories
    category_colors = {
        'Membrane/Cell Wall': '#E41A1C',
        'Regulation': '#377EB8',
        'Metabolism': '#4DAF4A',
        'DNA/RNA': '#984EA3',
        'Stress Response': '#FF7F00',
        'Virulence': '#A65628',
        'Other/Unknown': '#999999'
    }

    bottom = np.zeros(len(treatments))
    for cat in categories:
        pct_col = f'{cat}_pct'
        values = df_enrichment[pct_col].values
        ax.bar(x, values, width, label=cat, bottom=bottom, color=category_colors[cat])
        bottom += values

    ax.set_ylabel('Percentage of Mutated Genes', fontsize=12)
    ax.set_xlabel('Treatment', fontsize=12)
    ax.set_title('Functional Category Enrichment by Treatment\n(Non-synonymous mutations)',
                 fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(treatments, rotation=45, ha='right', fontsize=10)
    ax.legend(loc='upper right', bbox_to_anchor=(1.25, 1), fontsize=9)
    ax.set_ylim(0, 105)

    # Add gene counts on top
    for i, treatment in enumerate(treatments):
        total = df_enrichment[df_enrichment['treatment'] == treatment]['total'].values[0]
        ax.text(i, 102, f'n={total}', ha='center', fontsize=8)

    plt.tight_layout()
    plt.savefig(output_dir / 'functional_enrichment_stacked.png', dpi=150,
                bbox_inches='tight', facecolor='white')
    plt.close()

    # Create individual pie charts for singles
    singles = [t for t in treatments if '_' not in t or t.startswith('NPSA')]
    if singles:
        n_singles = len(singles)
        cols = min(3, n_singles)
        rows = (n_singles + cols - 1) // cols

        fig, axes = plt.subplots(rows, cols, figsize=(5*cols, 5*rows))
        if n_singles == 1:
            axes = np.array([axes])
        axes = axes.flatten()

        for idx, treatment in enumerate(singles):
            ax = axes[idx]
            row_data = df_enrichment[df_enrichment['treatment'] == treatment].iloc[0]

            sizes = [row_data.get(cat, 0) for cat in categories]
            colors = [category_colors[cat] for cat in categories]
            labels_with_pct = [f'{cat}\n({row_data.get(f"{cat}_pct", 0):.1f}%)'
                              if row_data.get(cat, 0) > 0 else '' for cat in categories]

            # Only show non-zero wedges
            non_zero = [(s, c, l) for s, c, l in zip(sizes, colors, labels_with_pct) if s > 0]
            if non_zero:
                sizes_nz, colors_nz, labels_nz = zip(*non_zero)
                ax.pie(sizes_nz, colors=colors_nz, labels=labels_nz,
                       autopct='', startangle=90, textprops={'fontsize': 8})

            color = EVO_COLOR_MAP.get(treatment, '#888888')
            ax.set_title(treatment, fontsize=12, fontweight='bold', color=color)

        # Hide unused axes
        for idx in range(len(singles), len(axes)):
            axes[idx].axis('off')

        plt.suptitle('Functional Categories by Single AMP Treatment',
                    fontsize=14, fontweight='bold', y=1.02)
        plt.tight_layout()
        plt.savefig(output_dir / 'functional_enrichment_pies_singles.png', dpi=150,
                   bbox_inches='tight', facecolor='white')
        plt.close()

    return df_enrichment


def calculate_replicate_consistency(
    sample_genes: Dict[str, Set[str]],
    output_dir: Path
) -> pd.DataFrame:
    """
    Calculate replicate consistency score for each treatment.
    Score = average pairwise Dice similarity within replicates.
    """
    import matplotlib.pyplot as plt

    # Group samples by treatment
    treatment_samples = defaultdict(list)
    for sample in sample_genes.keys():
        treatment = extract_treatment_from_sample(sample)
        treatment_samples[treatment].append(sample)

    consistency_data = []

    for treatment, samples in sorted(treatment_samples.items()):
        n_samples = len(samples)
        if n_samples < 2:
            consistency_data.append({
                'treatment': treatment,
                'n_replicates': n_samples,
                'consistency_score': np.nan,
                'std_dev': np.nan,
                'min_similarity': np.nan,
                'max_similarity': np.nan,
                'n_shared_genes': 0,
                'n_union_genes': len(sample_genes.get(samples[0], set())) if samples else 0
            })
            continue

        # Calculate all pairwise similarities
        similarities = []
        for i in range(n_samples):
            for j in range(i + 1, n_samples):
                sim = dice_similarity(sample_genes[samples[i]], sample_genes[samples[j]])
                similarities.append(sim)

        # Calculate shared and union genes
        all_gene_sets = [sample_genes[s] for s in samples]
        shared_genes = set.intersection(*all_gene_sets) if all_gene_sets else set()
        union_genes = set.union(*all_gene_sets) if all_gene_sets else set()

        consistency_data.append({
            'treatment': treatment,
            'n_replicates': n_samples,
            'consistency_score': np.mean(similarities),
            'std_dev': np.std(similarities),
            'min_similarity': np.min(similarities),
            'max_similarity': np.max(similarities),
            'n_shared_genes': len(shared_genes),
            'n_union_genes': len(union_genes)
        })

    df_consistency = pd.DataFrame(consistency_data)
    df_consistency.to_csv(output_dir / 'replicate_consistency.csv', index=False)

    # Create bar chart
    fig, ax = plt.subplots(figsize=(14, 8))

    valid_data = df_consistency[~df_consistency['consistency_score'].isna()]
    treatments = valid_data['treatment'].tolist()
    scores = valid_data['consistency_score'].tolist()
    std_devs = valid_data['std_dev'].tolist()

    x = np.arange(len(treatments))

    # Color bars by treatment type
    colors = []
    for t in treatments:
        if t in EVO_COLOR_MAP:
            colors.append(EVO_COLOR_MAP[t])
        elif t in COMBO_COLOR_MAP:
            colors.append(COMBO_COLOR_MAP[t])
        else:
            colors.append('#888888')

    bars = ax.bar(x, scores, width=0.6, color=colors, edgecolor='black', linewidth=1)
    ax.errorbar(x, scores, yerr=std_devs, fmt='none', color='black', capsize=5)

    ax.set_ylabel('Replicate Consistency Score\n(Mean Pairwise Dice Similarity)', fontsize=12)
    ax.set_xlabel('Treatment', fontsize=12)
    ax.set_title('Replicate Consistency Across Treatments\n(Higher = more consistent mutation profiles)',
                 fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(treatments, rotation=45, ha='right', fontsize=10)
    ax.set_ylim(0, 1.1)
    ax.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5, label='50% threshold')

    # Add score labels on bars
    for i, (bar, score) in enumerate(zip(bars, scores)):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.05,
               f'{score:.2f}', ha='center', fontsize=9)

    plt.tight_layout()
    plt.savefig(output_dir / 'replicate_consistency.png', dpi=150,
                bbox_inches='tight', facecolor='white')
    plt.close()

    return df_consistency


def calculate_parent_dominance_index(
    combo_treatment_genes: Dict[str, Set[str]],
    single_treatment_genes: Dict[str, Set[str]],
    output_dir: Path
) -> pd.DataFrame:
    """
    Calculate Parent Dominance Index for each combination.
    PDI = (|combo ∩ P1| - |combo ∩ P2|) / (|combo ∩ P1| + |combo ∩ P2|)
    Range: -1 (P2 dominant) to +1 (P1 dominant), 0 = balanced
    """
    import matplotlib.pyplot as plt

    dominance_data = []

    for combo, combo_genes in sorted(combo_treatment_genes.items()):
        p1, p2 = get_parent_treatments(combo)

        if p1 not in single_treatment_genes or p2 not in single_treatment_genes:
            continue

        p1_genes = single_treatment_genes[p1]
        p2_genes = single_treatment_genes[p2]

        # Overlap with each parent
        overlap_p1 = len(combo_genes & p1_genes)
        overlap_p2 = len(combo_genes & p2_genes)

        # Parent Dominance Index
        if overlap_p1 + overlap_p2 > 0:
            pdi = (overlap_p1 - overlap_p2) / (overlap_p1 + overlap_p2)
        else:
            pdi = 0

        # Classification
        if pdi > 0.3:
            classification = f'{p1}-dominant'
        elif pdi < -0.3:
            classification = f'{p2}-dominant'
        else:
            classification = 'Balanced'

        dominance_data.append({
            'combination': combo,
            'parent1': p1,
            'parent2': p2,
            'overlap_parent1': overlap_p1,
            'overlap_parent2': overlap_p2,
            'parent_dominance_index': pdi,
            'classification': classification,
            'total_combo_genes': len(combo_genes),
            'novel_genes': len(combo_genes - p1_genes - p2_genes)
        })

    df_dominance = pd.DataFrame(dominance_data)
    df_dominance.to_csv(output_dir / 'parent_dominance_index.csv', index=False)

    if df_dominance.empty:
        return df_dominance

    # Create visualization
    fig, ax = plt.subplots(figsize=(14, 8))

    combos = df_dominance['combination'].tolist()
    pdi_values = df_dominance['parent_dominance_index'].tolist()

    x = np.arange(len(combos))

    # Color by dominance direction
    colors = []
    for pdi in pdi_values:
        if pdi > 0.3:
            colors.append('#729ECEFF')  # Parent1 dominant (blue)
        elif pdi < -0.3:
            colors.append('#FF9E4AFF')  # Parent2 dominant (orange)
        else:
            colors.append('#67BF5CFF')  # Balanced (green)

    bars = ax.bar(x, pdi_values, width=0.6, color=colors, edgecolor='black', linewidth=1)

    ax.axhline(y=0, color='black', linewidth=1)
    ax.axhline(y=0.3, color='gray', linestyle='--', alpha=0.5)
    ax.axhline(y=-0.3, color='gray', linestyle='--', alpha=0.5)

    ax.set_ylabel('Parent Dominance Index\n(+1 = P1 dominant, -1 = P2 dominant)', fontsize=12)
    ax.set_xlabel('Combination Treatment', fontsize=12)
    ax.set_title('Parent Dominance Index\n(Which parent\'s mutation profile dominates?)',
                 fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(combos, rotation=45, ha='right', fontsize=10)
    ax.set_ylim(-1.1, 1.1)

    # Add parent labels
    for i, row in df_dominance.iterrows():
        combo = row['combination']
        p1, p2 = row['parent1'], row['parent2']
        pdi = row['parent_dominance_index']
        y_pos = pdi + 0.1 if pdi >= 0 else pdi - 0.15
        ax.text(i, y_pos, f'{p1}/{p2}', ha='center', fontsize=8, rotation=90)

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#729ECEFF', edgecolor='black', label='Parent 1 dominant (PDI > 0.3)'),
        Patch(facecolor='#67BF5CFF', edgecolor='black', label='Balanced (-0.3 ≤ PDI ≤ 0.3)'),
        Patch(facecolor='#FF9E4AFF', edgecolor='black', label='Parent 2 dominant (PDI < -0.3)')
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=10)

    plt.tight_layout()
    plt.savefig(output_dir / 'parent_dominance_index.png', dpi=150,
                bbox_inches='tight', facecolor='white')
    plt.close()

    return df_dominance


def create_venn_diagrams(
    combo_treatment_genes: Dict[str, Set[str]],
    single_treatment_genes: Dict[str, Set[str]],
    output_dir: Path
):
    """
    Create Venn diagrams for each combination showing overlap with parents.
    """
    import matplotlib.pyplot as plt

    try:
        from matplotlib_venn import venn3, venn3_circles
        has_venn = True
    except ImportError:
        has_venn = False
        print("   Note: matplotlib_venn not installed. Creating alternative overlap plots.")

    venn_dir = output_dir / 'venn_diagrams'
    venn_dir.mkdir(exist_ok=True)

    valid_combos = []
    for combo in sorted(combo_treatment_genes.keys()):
        p1, p2 = get_parent_treatments(combo)
        if p1 in single_treatment_genes and p2 in single_treatment_genes:
            valid_combos.append(combo)

    if not valid_combos:
        return

    # Create individual Venn diagrams
    for combo in valid_combos:
        p1, p2 = get_parent_treatments(combo)

        combo_genes = combo_treatment_genes[combo]
        p1_genes = single_treatment_genes[p1]
        p2_genes = single_treatment_genes[p2]

        fig, ax = plt.subplots(figsize=(10, 8))

        if has_venn:
            # Create Venn diagram
            venn = venn3(
                [combo_genes, p1_genes, p2_genes],
                set_labels=(combo, p1, p2),
                ax=ax
            )
            venn3_circles([combo_genes, p1_genes, p2_genes], ax=ax, linewidth=2)

            # Color the circles
            if venn.get_patch_by_id('100'):
                venn.get_patch_by_id('100').set_color(COMBO_COLOR_MAP.get(combo, '#888888'))
            if venn.get_patch_by_id('010'):
                venn.get_patch_by_id('010').set_color(EVO_COLOR_MAP.get(p1, '#888888'))
            if venn.get_patch_by_id('001'):
                venn.get_patch_by_id('001').set_color(EVO_COLOR_MAP.get(p2, '#888888'))
        else:
            # Alternative: bar chart showing overlaps
            only_combo = len(combo_genes - p1_genes - p2_genes)
            only_p1 = len(p1_genes - combo_genes - p2_genes)
            only_p2 = len(p2_genes - combo_genes - p1_genes)
            combo_p1 = len((combo_genes & p1_genes) - p2_genes)
            combo_p2 = len((combo_genes & p2_genes) - p1_genes)
            p1_p2 = len((p1_genes & p2_genes) - combo_genes)
            all_three = len(combo_genes & p1_genes & p2_genes)

            categories = [f'{combo}\nonly', f'{p1}\nonly', f'{p2}\nonly',
                         f'{combo}∩{p1}', f'{combo}∩{p2}', f'{p1}∩{p2}', 'All three']
            values = [only_combo, only_p1, only_p2, combo_p1, combo_p2, p1_p2, all_three]
            colors = [COMBO_COLOR_MAP.get(combo, '#888888'),
                     EVO_COLOR_MAP.get(p1, '#888888'),
                     EVO_COLOR_MAP.get(p2, '#888888'),
                     '#9999CC', '#CC9999', '#99CC99', '#666666']

            bars = ax.bar(categories, values, color=colors, edgecolor='black')
            ax.set_ylabel('Number of Genes', fontsize=12)
            ax.set_title(f'{combo} vs Parents ({p1}, {p2})\nGene Overlap',
                        fontsize=14, fontweight='bold')

            # Add value labels
            for bar, val in zip(bars, values):
                if val > 0:
                    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                           str(val), ha='center', fontsize=10)

        ax.set_title(f'{combo} Mutation Overlap\n(Combination vs Parent Treatments: {p1}, {p2})',
                    fontsize=14, fontweight='bold')

        plt.tight_layout()
        plt.savefig(venn_dir / f'venn_{combo}.png', dpi=150,
                   bbox_inches='tight', facecolor='white')
        plt.close()

    # Create summary figure with all Venn diagrams
    n_combos = len(valid_combos)
    cols = min(3, n_combos)
    rows = (n_combos + cols - 1) // cols

    fig, axes = plt.subplots(rows, cols, figsize=(6*cols, 5*rows))
    if n_combos == 1:
        axes = np.array([axes])
    axes = axes.flatten() if n_combos > 1 else [axes]

    for idx, combo in enumerate(valid_combos):
        ax = axes[idx]
        p1, p2 = get_parent_treatments(combo)

        combo_genes = combo_treatment_genes[combo]
        p1_genes = single_treatment_genes[p1]
        p2_genes = single_treatment_genes[p2]

        # Create simple bar representation
        only_combo = len(combo_genes - p1_genes - p2_genes)
        shared_p1 = len(combo_genes & p1_genes)
        shared_p2 = len(combo_genes & p2_genes)
        shared_both = len(combo_genes & p1_genes & p2_genes)

        categories = ['Novel', f'With\n{p1}', f'With\n{p2}', 'With\nBoth']
        values = [only_combo, shared_p1 - shared_both, shared_p2 - shared_both, shared_both]
        colors = ['#ED665DFF', EVO_COLOR_MAP.get(p1, '#888888'),
                 EVO_COLOR_MAP.get(p2, '#888888'), '#67BF5CFF']

        bars = ax.bar(categories, values, color=colors, edgecolor='black')

        for bar, val in zip(bars, values):
            if val > 0:
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                       str(val), ha='center', fontsize=9)

        ax.set_title(combo, fontsize=12, fontweight='bold',
                    color=COMBO_COLOR_MAP.get(combo, '#888888'))
        ax.set_ylabel('Genes', fontsize=10)

    # Hide unused axes
    for idx in range(len(valid_combos), len(axes)):
        axes[idx].axis('off')

    plt.suptitle('Combination Mutations: Overlap with Parent Treatments',
                fontsize=16, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(output_dir / 'venn_summary.png', dpi=150,
               bbox_inches='tight', facecolor='white')
    plt.close()


def create_phylogenetic_tree(
    treatment_genes: Dict[str, Set[str]],
    output_dir: Path
):
    """
    Create hierarchical clustering tree of treatments based on mutation similarity.
    """
    import matplotlib.pyplot as plt
    from scipy.cluster.hierarchy import dendrogram, linkage
    from scipy.spatial.distance import squareform

    treatments = sorted(treatment_genes.keys())
    n = len(treatments)

    if n < 3:
        print("   Not enough treatments for clustering tree")
        return

    # Build distance matrix (1 - Dice similarity)
    dist_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i != j:
                sim = dice_similarity(treatment_genes[treatments[i]],
                                     treatment_genes[treatments[j]])
                dist_matrix[i, j] = 1 - sim

    # Convert to condensed form for linkage
    condensed_dist = squareform(dist_matrix)

    # Perform hierarchical clustering
    linkage_matrix = linkage(condensed_dist, method='average')

    # Create figure
    fig, ax = plt.subplots(figsize=(14, 10))

    # Create dendrogram
    dendro = dendrogram(
        linkage_matrix,
        labels=treatments,
        ax=ax,
        leaf_rotation=45,
        leaf_font_size=11,
        color_threshold=0.7
    )

    # Color the labels
    xlbls = ax.get_xmajorticklabels()
    for lbl in xlbls:
        treatment = lbl.get_text()
        if treatment in EVO_COLOR_MAP:
            lbl.set_color(EVO_COLOR_MAP[treatment])
        elif treatment in COMBO_COLOR_MAP:
            lbl.set_color(COMBO_COLOR_MAP[treatment])
        lbl.set_fontweight('bold')

    ax.set_ylabel('Distance (1 - Dice Similarity)', fontsize=12)
    ax.set_xlabel('Treatment', fontsize=12)
    ax.set_title('Hierarchical Clustering of Treatments\n(Based on Shared Mutated Genes)',
                fontsize=14, fontweight='bold')

    # Add horizontal line at common threshold
    ax.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5, label='50% similarity')
    ax.axhline(y=0.7, color='red', linestyle='--', alpha=0.5, label='30% similarity')
    ax.legend(loc='upper right')

    plt.tight_layout()
    plt.savefig(output_dir / 'treatment_clustering_tree.png', dpi=150,
                bbox_inches='tight', facecolor='white')
    plt.close()

    # Also create circular dendrogram if possible
    try:
        fig, ax = plt.subplots(figsize=(12, 12), subplot_kw={'projection': 'polar'})

        # Create circular layout
        n_leaves = len(treatments)
        angles = np.linspace(0, 2*np.pi, n_leaves, endpoint=False)

        # Plot as polar scatter with labels
        for i, (angle, treatment) in enumerate(zip(angles, dendro['ivl'])):
            color = EVO_COLOR_MAP.get(treatment, COMBO_COLOR_MAP.get(treatment, '#888888'))
            ax.scatter(angle, 1, s=200, c=[color], edgecolors='black', zorder=3)
            ax.annotate(treatment, (angle, 1.15), ha='center', va='center',
                       fontsize=10, fontweight='bold', color=color)

        # Add connections based on clustering
        # This is simplified - just show the clustering order
        ax.set_ylim(0, 1.5)
        ax.set_title('Circular Treatment Clustering\n(Based on Mutation Similarity)',
                    fontsize=14, fontweight='bold', pad=20)
        ax.axis('off')

        plt.savefig(output_dir / 'treatment_clustering_circular.png', dpi=150,
                   bbox_inches='tight', facecolor='white')
        plt.close()
    except Exception:
        pass  # Skip circular plot if it fails


def main():
    if len(sys.argv) < 4:
        print("Usage: python analyze_amp_combinations.py <singles_dir> <combinations_dir> <output_dir>")
        print()
        print("Arguments:")
        print("  singles_dir:      Directory with single AMP results")
        print("  combinations_dir: Directory with combination results")
        print("  output_dir:       Output directory for analysis")
        print()
        print("Example:")
        print("  python analyze_amp_combinations.py results_single_amps results_combinations analysis_output")
        sys.exit(1)

    singles_dir = Path(sys.argv[1])
    combinations_dir = Path(sys.argv[2])
    output_dir = Path(sys.argv[3])

    output_dir.mkdir(parents=True, exist_ok=True)

    print("="*70)
    print("AMP COMBINATION MUTATION ANALYSIS")
    print("="*70)

    # Load single AMP data
    print("\n1. Loading single AMP mutations...")
    single_sample_genes = load_mutations_from_directory(singles_dir)
    print(f"   Loaded {len(single_sample_genes)} single AMP samples")

    # Load combination data
    print("\n2. Loading combination mutations...")
    combo_sample_genes = load_mutations_from_directory(combinations_dir)
    print(f"   Loaded {len(combo_sample_genes)} combination samples")

    # Aggregate by treatment
    print("\n3. Aggregating by treatment...")
    single_treatment_genes = aggregate_by_treatment(single_sample_genes)
    combo_treatment_genes = aggregate_by_treatment(combo_sample_genes)

    print(f"   Single treatments: {list(single_treatment_genes.keys())}")
    print(f"   Combinations: {list(combo_treatment_genes.keys())}")

    # All treatments combined
    all_treatment_genes = {**single_treatment_genes, **combo_treatment_genes}

    # Gene counts per treatment
    single_gene_counts = get_treatment_gene_counts(single_sample_genes)
    combo_gene_counts = get_treatment_gene_counts(combo_sample_genes)
    all_gene_counts = {**single_gene_counts, **combo_gene_counts}

    # ===== ANALYSIS 1: Single AMP Similarity =====
    print("\n" + "="*70)
    print("ANALYSIS 1: Single AMP Similarity")
    print("="*70)

    singles_sim_path = output_dir / 'singles_similarity_heatmap.png'
    singles_sim_df = create_similarity_heatmap(
        single_treatment_genes,
        singles_sim_path,
        title="Single AMP Treatment Similarity\n(Dice coefficient based on shared mutated genes)"
    )
    singles_sim_df.to_csv(output_dir / 'singles_similarity_matrix.csv')
    print(f"   Saved: {singles_sim_path}")

    # ===== ANALYSIS 2: Combination vs Parents =====
    print("\n" + "="*70)
    print("ANALYSIS 2: Combination vs Parents Analysis")
    print("="*70)

    combo_analysis = {}
    combo_results = []

    for combo, combo_genes in combo_treatment_genes.items():
        p1, p2 = get_parent_treatments(combo)

        if p1 not in single_treatment_genes or p2 not in single_treatment_genes:
            print(f"   Skipping {combo}: parents {p1}, {p2} not found in singles")
            continue

        analysis = analyze_combination_vs_parents(
            combo_genes,
            single_treatment_genes[p1],
            single_treatment_genes[p2]
        )
        combo_analysis[combo] = analysis

        total = len(analysis['total'])
        result = {
            'combination': combo,
            'parent1': p1,
            'parent2': p2,
            'total_genes': total,
            'parent1_only': len(analysis['parent1_only']),
            'parent2_only': len(analysis['parent2_only']),
            'shared_parents': len(analysis['shared_parents']),
            'novel': len(analysis['novel']),
            'parent1_pct': len(analysis['parent1_only']) / total * 100 if total > 0 else 0,
            'parent2_pct': len(analysis['parent2_only']) / total * 100 if total > 0 else 0,
            'shared_pct': len(analysis['shared_parents']) / total * 100 if total > 0 else 0,
            'novel_pct': len(analysis['novel']) / total * 100 if total > 0 else 0
        }
        combo_results.append(result)

        print(f"   {combo}: {total} genes total")
        print(f"      {p1}-only: {len(analysis['parent1_only'])} ({result['parent1_pct']:.1f}%)")
        print(f"      {p2}-only: {len(analysis['parent2_only'])} ({result['parent2_pct']:.1f}%)")
        print(f"      Shared: {len(analysis['shared_parents'])} ({result['shared_pct']:.1f}%)")
        print(f"      Novel: {len(analysis['novel'])} ({result['novel_pct']:.1f}%)")

    # Save combo analysis
    combo_df = pd.DataFrame(combo_results)
    combo_df.to_csv(output_dir / 'combination_vs_parents.csv', index=False)

    # Create breakdown plot
    if combo_analysis:
        breakdown_path = output_dir / 'combination_breakdown.png'
        create_combination_breakdown_plot(combo_analysis, breakdown_path)
        print(f"   Saved: {breakdown_path}")

    # ===== ANALYSIS 3: All Treatments Similarity Network =====
    print("\n" + "="*70)
    print("ANALYSIS 3: Treatment Similarity Network")
    print("="*70)

    # All treatments heatmap
    all_sim_path = output_dir / 'all_treatments_similarity_heatmap.png'
    all_sim_df = create_similarity_heatmap(
        all_treatment_genes,
        all_sim_path,
        title="All Treatments Similarity\n(Singles + Combinations)"
    )
    all_sim_df.to_csv(output_dir / 'all_treatments_similarity_matrix.csv')
    print(f"   Saved: {all_sim_path}")

    # Network
    network_path = output_dir / 'treatment_similarity_network.png'
    create_similarity_network(
        all_treatment_genes,
        network_path,
        min_edge_weight=0.15,
        title="Treatment Similarity Network\n(Circles = Singles, Squares = Combinations)"
    )
    print(f"   Saved: {network_path}")

    # ===== ANALYSIS 4: Convergent Evolution =====
    print("\n" + "="*70)
    print("ANALYSIS 4: Convergent Evolution")
    print("="*70)

    convergent_df = find_convergent_genes(all_gene_counts, min_treatments=3)
    convergent_df.to_csv(output_dir / 'convergent_genes.csv', index=False)
    print(f"   Found {len(convergent_df)} genes mutated in ≥3 treatments")

    if not convergent_df.empty:
        print("\n   Top convergent genes:")
        for _, row in convergent_df.head(10).iterrows():
            print(f"      {row['gene']}: {row['n_treatments']} treatments")

    # Convergent genes heatmap
    conv_heatmap_path = output_dir / 'convergent_genes_heatmap.png'
    create_convergent_genes_heatmap(all_gene_counts, conv_heatmap_path, top_n=25)
    print(f"   Saved: {conv_heatmap_path}")

    # ===== ANALYSIS 5: Functional Category Enrichment =====
    print("\n" + "="*70)
    print("ANALYSIS 5: Functional Category Enrichment")
    print("="*70)

    # Load full mutation data for functional analysis
    single_mutations = load_mutations_with_details(singles_dir)
    combo_mutations = load_mutations_with_details(combinations_dir)
    all_mutations = {**single_mutations, **combo_mutations}

    df_enrichment = analyze_functional_enrichment(all_mutations, output_dir)
    print(f"   Analyzed {len(df_enrichment)} treatments for functional enrichment")
    print(f"   Saved: functional_enrichment.csv, functional_enrichment_stacked.png")

    # ===== ANALYSIS 6: Replicate Consistency Score =====
    print("\n" + "="*70)
    print("ANALYSIS 6: Replicate Consistency Score")
    print("="*70)

    # Combine all samples for consistency analysis
    all_sample_genes = {**single_sample_genes, **combo_sample_genes}
    df_consistency = calculate_replicate_consistency(all_sample_genes, output_dir)

    print("\n   Replicate Consistency Scores:")
    for _, row in df_consistency.iterrows():
        if not np.isnan(row['consistency_score']):
            print(f"      {row['treatment']}: {row['consistency_score']:.3f} "
                  f"(n={row['n_replicates']} replicates)")
    print(f"   Saved: replicate_consistency.csv, replicate_consistency.png")

    # ===== ANALYSIS 7: Parent Dominance Index =====
    print("\n" + "="*70)
    print("ANALYSIS 7: Parent Dominance Index")
    print("="*70)

    df_dominance = calculate_parent_dominance_index(
        combo_treatment_genes,
        single_treatment_genes,
        output_dir
    )

    if not df_dominance.empty:
        print("\n   Parent Dominance Index:")
        for _, row in df_dominance.iterrows():
            print(f"      {row['combination']}: PDI={row['parent_dominance_index']:.3f} "
                  f"({row['classification']})")
    print(f"   Saved: parent_dominance_index.csv, parent_dominance_index.png")

    # ===== ANALYSIS 8: Venn Diagrams =====
    print("\n" + "="*70)
    print("ANALYSIS 8: Venn Diagrams (Combination vs Parents)")
    print("="*70)

    create_venn_diagrams(combo_treatment_genes, single_treatment_genes, output_dir)
    print(f"   Saved: venn_summary.png, venn_diagrams/ directory")

    # ===== ANALYSIS 9: Phylogenetic/Clustering Tree =====
    print("\n" + "="*70)
    print("ANALYSIS 9: Hierarchical Clustering Tree")
    print("="*70)

    create_phylogenetic_tree(all_treatment_genes, output_dir)
    print(f"   Saved: treatment_clustering_tree.png")

    # ===== SUMMARY =====
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)

    print(f"\nOutput files saved to: {output_dir}")
    print("\nFiles created:")
    print("\n  Core Analyses:")
    print("  - singles_similarity_heatmap.png       : Single AMP similarity heatmap")
    print("  - singles_similarity_matrix.csv        : Single AMP similarity data")
    print("  - combination_vs_parents.csv           : Combination breakdown data")
    print("  - combination_breakdown.png            : Stacked bar chart of mutation origins")
    print("  - all_treatments_similarity_heatmap.png: All treatments heatmap")
    print("  - all_treatments_similarity_matrix.csv : All treatments similarity data")
    print("  - treatment_similarity_network.png     : Network visualization")
    print("  - convergent_genes.csv                 : Convergent genes list")
    print("  - convergent_genes_heatmap.png         : Convergent genes heatmap")
    print("\n  Additional Publication-Worthy Analyses:")
    print("  - functional_enrichment.csv            : Functional category counts")
    print("  - functional_enrichment_stacked.png    : Functional enrichment stacked bars")
    print("  - functional_enrichment_pies_singles.png: Pie charts for single AMPs")
    print("  - replicate_consistency.csv            : Replicate consistency scores")
    print("  - replicate_consistency.png            : Consistency bar chart")
    print("  - parent_dominance_index.csv           : Parent dominance data")
    print("  - parent_dominance_index.png           : PDI visualization")
    print("  - venn_summary.png                     : Overview of all combination overlaps")
    print("  - venn_diagrams/                       : Individual Venn diagrams per combo")
    print("  - treatment_clustering_tree.png        : Hierarchical clustering dendrogram")


if __name__ == '__main__':
    main()
