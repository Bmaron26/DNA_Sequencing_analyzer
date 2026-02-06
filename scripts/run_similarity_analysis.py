#!/usr/bin/env python3
"""
Standalone script for mutation profile similarity analysis.

Calculates Dice similarity between mutation profiles and creates:
1. Sample-level clustered heatmap
2. Treatment-level heatmap
3. Similarity network visualization

Usage:
    python scripts/run_similarity_analysis.py <input_file> <output_dir>

Example:
    python scripts/run_similarity_analysis.py \
        results_single_amps/combined_analysis/freebayes/all_mutations_freebayes.csv \
        similarity_results
"""

import sys
import os
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from itertools import combinations
from collections import defaultdict

import pandas as pd
import numpy as np


def dice_similarity(set1: set, set2: set) -> float:
    """Calculate Dice similarity coefficient between two sets."""
    if len(set1) == 0 and len(set2) == 0:
        return 1.0
    if len(set1) == 0 or len(set2) == 0:
        return 0.0
    intersection = len(set1 & set2)
    return (2 * intersection) / (len(set1) + len(set2))


def load_long_format(filepath: Path) -> Dict[str, set]:
    """Load mutations from long format CSV and extract mutation sets per sample."""
    df = pd.read_csv(filepath)

    print(f"  Columns found: {list(df.columns)[:10]}...")

    # Detect column names
    if 'position' in df.columns:
        pos_col = 'position'
    elif 'POS' in df.columns:
        pos_col = 'POS'
    else:
        pos_col = 'pos'

    if 'reference' in df.columns:
        ref_col = 'reference'
    elif 'REF' in df.columns:
        ref_col = 'REF'
    else:
        ref_col = 'ref'

    if 'alternative' in df.columns:
        alt_col = 'alternative'
    elif 'ALT' in df.columns:
        alt_col = 'ALT'
    else:
        alt_col = 'alt'

    if 'sample' in df.columns:
        sample_col = 'sample'
    elif 'SAMPLE' in df.columns:
        sample_col = 'SAMPLE'
    else:
        raise ValueError("No sample column found in file")

    print(f"  Using columns: pos={pos_col}, ref={ref_col}, alt={alt_col}, sample={sample_col}")

    # Create variant keys and group by sample
    sample_mutations = {}
    for _, row in df.iterrows():
        sample = row[sample_col]
        if pd.isna(sample):
            continue

        pos = row[pos_col]
        ref = row[ref_col] if not pd.isna(row[ref_col]) else ''
        alt = row[alt_col] if not pd.isna(row[alt_col]) else ''

        key = f"{int(pos)}_{ref}_{alt}"

        if sample not in sample_mutations:
            sample_mutations[sample] = set()
        sample_mutations[sample].add(key)

    return sample_mutations


def calculate_pairwise_similarity(sample_mutations: Dict[str, set]) -> pd.DataFrame:
    """Calculate pairwise Dice similarity between all samples."""
    samples = sorted(sample_mutations.keys())
    n = len(samples)
    sim_matrix = np.zeros((n, n))

    for i, s1 in enumerate(samples):
        for j, s2 in enumerate(samples):
            if i == j:
                sim_matrix[i, j] = 1.0
            elif i < j:
                sim = dice_similarity(sample_mutations[s1], sample_mutations[s2])
                sim_matrix[i, j] = sim
                sim_matrix[j, i] = sim

    return pd.DataFrame(sim_matrix, index=samples, columns=samples)


def extract_treatment_from_sample(sample_name: str) -> str:
    """Extract treatment name from sample name (e.g., 'BmKn1_1' -> 'BmKn1')."""
    parts = sample_name.rsplit('_', 1)
    if len(parts) == 2 and parts[1].isdigit():
        return parts[0]
    return sample_name


def calculate_treatment_similarity(similarity_matrix: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Calculate average similarity within and between treatment groups."""
    samples = similarity_matrix.index.tolist()
    sample_to_treatment = {s: extract_treatment_from_sample(s) for s in samples}
    treatments = sorted(set(sample_to_treatment.values()))

    treatment_samples = {t: [s for s in samples if sample_to_treatment[s] == t] for t in treatments}

    n_treatments = len(treatments)
    treatment_sim = np.zeros((n_treatments, n_treatments))
    stats_data = []

    for i, t1 in enumerate(treatments):
        for j, t2 in enumerate(treatments):
            samples1 = treatment_samples[t1]
            samples2 = treatment_samples[t2]

            similarities = []
            for s1 in samples1:
                for s2 in samples2:
                    if s1 != s2:
                        similarities.append(similarity_matrix.loc[s1, s2])

            if similarities:
                mean_sim = np.mean(similarities)
                std_sim = np.std(similarities)
            else:
                mean_sim = 1.0 if t1 == t2 else 0.0
                std_sim = 0.0

            treatment_sim[i, j] = mean_sim

            if i <= j:
                stats_data.append({
                    'treatment1': t1,
                    'treatment2': t2,
                    'type': 'within' if t1 == t2 else 'between',
                    'mean_similarity': mean_sim,
                    'std_similarity': std_sim,
                    'n_comparisons': len(similarities)
                })

    treatment_sim_df = pd.DataFrame(treatment_sim, index=treatments, columns=treatments)
    stats_df = pd.DataFrame(stats_data)

    return treatment_sim_df, stats_df


def create_similarity_clustermap(similarity_matrix: pd.DataFrame, output_path: Path):
    """Create a clustermap visualization of sample similarities."""
    import matplotlib.pyplot as plt
    import seaborn as sns

    samples = similarity_matrix.index.tolist()
    treatments = [extract_treatment_from_sample(s) for s in samples]
    unique_treatments = sorted(set(treatments))

    colors = plt.cm.Set2(np.linspace(0, 1, len(unique_treatments)))
    treatment_colors = {t: colors[i] for i, t in enumerate(unique_treatments)}
    row_colors = [treatment_colors[extract_treatment_from_sample(s)] for s in samples]

    g = sns.clustermap(
        similarity_matrix,
        method='average',
        metric='euclidean',
        cmap='RdYlBu_r',
        vmin=0, vmax=1,
        row_colors=row_colors,
        col_colors=row_colors,
        figsize=(14, 12),
        dendrogram_ratio=(0.15, 0.15),
        cbar_pos=(0.02, 0.8, 0.03, 0.15),
        linewidths=0.5,
        xticklabels=True,
        yticklabels=True
    )

    g.fig.suptitle("Mutation Profile Similarity (Dice Coefficient)", fontsize=14, fontweight='bold', y=1.02)

    legend_handles = [plt.Rectangle((0, 0), 1, 1, facecolor=treatment_colors[t],
                                     edgecolor='black', linewidth=0.5)
                      for t in unique_treatments]
    g.ax_heatmap.legend(legend_handles, unique_treatments,
                        loc='upper left', bbox_to_anchor=(1.15, 1.0),
                        title='Treatment', frameon=True, fontsize=9)

    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right', fontsize=9)
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=9)

    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  Saved clustermap to {output_path}")


def create_treatment_heatmap(treatment_sim: pd.DataFrame, output_path: Path):
    """Create a heatmap of treatment-level similarities."""
    import matplotlib.pyplot as plt
    import seaborn as sns

    fig, ax = plt.subplots(figsize=(10, 8))

    sns.heatmap(
        treatment_sim,
        annot=True,
        fmt='.3f',
        cmap='RdYlBu_r',
        vmin=0, vmax=1,
        square=True,
        linewidths=1,
        cbar_kws={'label': 'Dice Similarity', 'shrink': 0.8},
        ax=ax,
        annot_kws={'fontsize': 11, 'fontweight': 'bold'}
    )

    ax.set_title("Treatment-Level Mutation Similarity", fontsize=14, fontweight='bold', pad=15)
    ax.set_xlabel('Treatment', fontsize=12)
    ax.set_ylabel('Treatment', fontsize=12)

    plt.setp(ax.get_xticklabels(), rotation=45, ha='right', fontsize=11)
    plt.setp(ax.get_yticklabels(), rotation=0, fontsize=11)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  Saved treatment heatmap to {output_path}")


def create_similarity_network(treatment_sim: pd.DataFrame, output_path: Path, min_edge_weight: float = 0.1):
    """Create a network visualization of treatment similarities."""
    import matplotlib.pyplot as plt
    import networkx as nx

    treatments = treatment_sim.index.tolist()
    G = nx.Graph()

    for t in treatments:
        within_sim = treatment_sim.loc[t, t]
        G.add_node(t, within_similarity=within_sim)

    for i, t1 in enumerate(treatments):
        for j, t2 in enumerate(treatments):
            if i < j:
                sim = treatment_sim.loc[t1, t2]
                if sim >= min_edge_weight:
                    G.add_edge(t1, t2, weight=sim)

    pos = nx.spring_layout(G, k=2, iterations=50, seed=42)

    fig, ax = plt.subplots(figsize=(12, 10))

    node_sizes = [1500 + 1000 * G.nodes[n].get('within_similarity', 0.5) for n in G.nodes()]
    colors = plt.cm.Set2(np.linspace(0, 1, len(treatments)))
    node_colors = {t: colors[i] for i, t in enumerate(treatments)}

    nx.draw_networkx_nodes(
        G, pos,
        node_size=node_sizes,
        node_color=[node_colors[n] for n in G.nodes()],
        edgecolors='black',
        linewidths=2,
        ax=ax
    )

    edges = G.edges(data=True)
    if edges:
        edge_weights = [d['weight'] for _, _, d in edges]
        max_weight = max(edge_weights) if edge_weights else 1

        for (u, v, d) in edges:
            weight = d['weight']
            width = 1 + 8 * (weight / max_weight)
            alpha = 0.3 + 0.7 * (weight / max_weight)

            if weight > 0.5:
                color = 'darkgreen'
            elif weight > 0.3:
                color = 'orange'
            else:
                color = 'lightgray'

            nx.draw_networkx_edges(G, pos, edgelist=[(u, v)], width=width, alpha=alpha, edge_color=color, ax=ax)

        edge_labels = {(u, v): f"{d['weight']:.2f}" for u, v, d in edges if d['weight'] > 0.15}
        nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=9, font_weight='bold', ax=ax)

    nx.draw_networkx_labels(G, pos, font_size=12, font_weight='bold', ax=ax)

    ax.set_title("AMP Mutation Similarity Network", fontsize=14, fontweight='bold', pad=15)

    legend_elements = [
        plt.Line2D([0], [0], color='darkgreen', linewidth=4, label='High similarity (>0.5)'),
        plt.Line2D([0], [0], color='orange', linewidth=3, label='Medium similarity (0.3-0.5)'),
        plt.Line2D([0], [0], color='lightgray', linewidth=2, label='Low similarity (<0.3)')
    ]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=10)

    ax.text(0.02, 0.02, 'Node size = within-treatment similarity\nEdge thickness = between-treatment similarity',
            transform=ax.transAxes, fontsize=9, verticalalignment='bottom',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    ax.axis('off')
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  Saved similarity network to {output_path}")


def main():
    if len(sys.argv) < 3:
        print("Usage: python run_similarity_analysis.py <input_file> <output_dir>")
        print()
        print("Input file: Long format CSV with 'sample', 'position', 'reference', 'alternative' columns")
        print()
        print("Example:")
        print("  python scripts/run_similarity_analysis.py \\")
        print("      results/summary_freebayes/all_mutations_freebayes.csv \\")
        print("      similarity_results")
        sys.exit(1)

    input_file = Path(sys.argv[1])
    output_dir = Path(sys.argv[2])

    if not input_file.exists():
        print(f"Error: Input file not found: {input_file}")
        sys.exit(1)

    output_dir.mkdir(parents=True, exist_ok=True)

    print("="*60)
    print("MUTATION PROFILE SIMILARITY ANALYSIS")
    print("="*60)
    print(f"Input: {input_file}")
    print(f"Output: {output_dir}")
    print()

    # Load data
    print("Loading mutations...")
    sample_mutations = load_long_format(input_file)
    print(f"  Found {len(sample_mutations)} samples")
    for sample, muts in list(sample_mutations.items())[:3]:
        print(f"    {sample}: {len(muts)} mutations")

    # Calculate pairwise similarities
    print("\nCalculating pairwise Dice similarities...")
    sim_matrix = calculate_pairwise_similarity(sample_mutations)

    sim_matrix_path = output_dir / 'sample_similarity_matrix.csv'
    sim_matrix.to_csv(sim_matrix_path)
    print(f"  Saved similarity matrix to {sim_matrix_path}")

    # Calculate treatment-level similarities
    print("\nCalculating treatment-level similarities...")
    treatment_sim, treatment_stats = calculate_treatment_similarity(sim_matrix)

    treatment_sim_path = output_dir / 'treatment_similarity_matrix.csv'
    treatment_sim.to_csv(treatment_sim_path)
    print(f"  Saved treatment similarity matrix to {treatment_sim_path}")

    treatment_stats_path = output_dir / 'treatment_similarity_stats.csv'
    treatment_stats.to_csv(treatment_stats_path, index=False)
    print(f"  Saved treatment statistics to {treatment_stats_path}")

    # Create visualizations
    print("\nCreating visualizations...")

    clustermap_path = output_dir / 'sample_similarity_clustermap.png'
    create_similarity_clustermap(sim_matrix, clustermap_path)

    treatment_heatmap_path = output_dir / 'treatment_similarity_heatmap.png'
    create_treatment_heatmap(treatment_sim, treatment_heatmap_path)

    network_path = output_dir / 'treatment_similarity_network.png'
    create_similarity_network(treatment_sim, network_path)

    # Print summary
    print()
    print("="*60)
    print("SIMILARITY ANALYSIS SUMMARY")
    print("="*60)

    within_stats = treatment_stats[treatment_stats['type'] == 'within']
    print("\nWithin-treatment similarity (replicate consistency):")
    for _, row in within_stats.iterrows():
        print(f"  {row['treatment1']}: {row['mean_similarity']:.3f} ± {row['std_similarity']:.3f}")

    between_stats = treatment_stats[treatment_stats['type'] == 'between'].copy()
    between_stats = between_stats.sort_values('mean_similarity', ascending=False)

    print("\nMost similar treatment pairs:")
    for _, row in between_stats.head(5).iterrows():
        print(f"  {row['treatment1']} - {row['treatment2']}: {row['mean_similarity']:.3f}")

    print("\nLeast similar treatment pairs:")
    for _, row in between_stats.tail(3).iterrows():
        print(f"  {row['treatment1']} - {row['treatment2']}: {row['mean_similarity']:.3f}")

    print()
    print("="*60)
    print("OUTPUT FILES")
    print("="*60)
    print(f"  {output_dir}/sample_similarity_matrix.csv")
    print(f"  {output_dir}/treatment_similarity_matrix.csv")
    print(f"  {output_dir}/treatment_similarity_stats.csv")
    print(f"  {output_dir}/sample_similarity_clustermap.png")
    print(f"  {output_dir}/treatment_similarity_heatmap.png")
    print(f"  {output_dir}/treatment_similarity_network.png")


if __name__ == '__main__':
    main()
