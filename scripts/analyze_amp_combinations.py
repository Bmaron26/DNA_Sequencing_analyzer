#!/usr/bin/env python3
"""
Comprehensive analysis of AMP resistance evolution mutations.

Analyzes:
1. Single AMP similarity (which AMPs target similar genes?)
2. Combination vs Parents (novel, parent1-like, parent2-like, shared)
3. Similarity networks (singles + combinations)
4. Convergent evolution analysis
5. Functional category enrichment

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
        Patch(facecolor='gray', edgecolor='black', label='Single AMP (circle)'),
        Patch(facecolor='gray', edgecolor='black', label='Combination (square)', marker='s'),
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

    # ===== SUMMARY =====
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)

    print(f"\nOutput files saved to: {output_dir}")
    print("\nFiles created:")
    print("  - singles_similarity_heatmap.png     : Single AMP similarity")
    print("  - singles_similarity_matrix.csv      : Single AMP similarity matrix")
    print("  - combination_vs_parents.csv         : Combination breakdown data")
    print("  - combination_breakdown.png          : Stacked bar chart")
    print("  - all_treatments_similarity_heatmap.png : All treatments heatmap")
    print("  - all_treatments_similarity_matrix.csv  : All treatments similarity")
    print("  - treatment_similarity_network.png   : Network visualization")
    print("  - convergent_genes.csv               : Convergent genes list")
    print("  - convergent_genes_heatmap.png       : Convergent genes heatmap")


if __name__ == '__main__':
    main()
