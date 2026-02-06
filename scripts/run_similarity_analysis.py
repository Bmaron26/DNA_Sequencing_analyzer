#!/usr/bin/env python3
"""
Standalone script for mutation profile similarity analysis.

Calculates Dice similarity between mutation profiles and creates:
1. Sample-level clustered heatmap
2. Treatment-level heatmap (aggregated across replicates)
3. Sample-level similarity network (individual replicates)
4. Treatment-level similarity network (aggregated)

Similarity is based on mutations in the SAME GENE (not exact position),
excluding synonymous mutations.

Usage:
    python scripts/run_similarity_analysis.py <input_file> <output_dir>

Example:
    python scripts/run_similarity_analysis.py \
        results/summary_freebayes/all_mutations_freebayes_filtered.csv \
        similarity_results
"""

import sys
import os
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
from itertools import combinations
from collections import defaultdict

import pandas as pd
import numpy as np
import re

# Standard color map for AMP treatments
EVO_COLOR_MAP = {
    'Mel': '#729ECEFF',
    'Pex': '#FF9E4AFF',
    'BmKn': '#67BF5CFF',
    'Puro': '#ED665DFF',
    'Pleu': '#AD8BC9FF',
    'Smp': '#A8786EFF'
}


def dice_similarity(set1: set, set2: set) -> float:
    """Calculate Dice similarity coefficient between two sets."""
    if len(set1) == 0 and len(set2) == 0:
        return 1.0
    if len(set1) == 0 or len(set2) == 0:
        return 0.0
    intersection = len(set1 & set2)
    return (2 * intersection) / (len(set1) + len(set2))


def is_synonymous(effect: str) -> bool:
    """Check if a mutation effect is synonymous (silent)."""
    if pd.isna(effect):
        return False
    effect_lower = str(effect).lower()
    synonymous_terms = ['synonymous', 'silent', 'syn_coding']
    return any(term in effect_lower for term in synonymous_terms)


def load_long_format(filepath: Path, use_genes: bool = True) -> Tuple[Dict[str, Set[str]], Dict[str, Set[str]]]:
    """
    Load mutations from long format CSV and extract mutation sets per sample.

    Args:
        filepath: Path to CSV file
        use_genes: If True, use gene names for similarity (excluding synonymous).
                   If False, use exact position+ref+alt.

    Returns:
        sample_genes: Dict mapping sample -> set of mutated genes (non-synonymous)
        sample_positions: Dict mapping sample -> set of position keys (for reference)
    """
    df = pd.read_csv(filepath)

    print(f"  Columns found: {list(df.columns)[:15]}...")
    print(f"  Total rows: {len(df)}")

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

    # Gene column
    if 'gene_name' in df.columns:
        gene_col = 'gene_name'
    elif 'GENE' in df.columns:
        gene_col = 'GENE'
    elif 'gene' in df.columns:
        gene_col = 'gene'
    elif 'locus_tag' in df.columns:
        gene_col = 'locus_tag'
    else:
        gene_col = None
        print("  WARNING: No gene column found, falling back to position-based similarity")

    # Effect column (to filter synonymous)
    if 'effect' in df.columns:
        effect_col = 'effect'
    elif 'EFFECT' in df.columns:
        effect_col = 'EFFECT'
    else:
        effect_col = None

    print(f"  Using columns: sample={sample_col}, gene={gene_col}, effect={effect_col}")

    # Create mutation sets per sample
    sample_genes = {}  # Gene-based (non-synonymous)
    sample_positions = {}  # Position-based (all mutations)

    synonymous_count = 0
    no_gene_count = 0

    for _, row in df.iterrows():
        sample = row[sample_col]
        if pd.isna(sample):
            continue

        # Initialize sample sets
        if sample not in sample_genes:
            sample_genes[sample] = set()
            sample_positions[sample] = set()

        # Position-based key
        pos = row[pos_col]
        ref = row[ref_col] if not pd.isna(row[ref_col]) else ''
        alt = row[alt_col] if not pd.isna(row[alt_col]) else ''
        pos_key = f"{int(pos)}_{ref}_{alt}"
        sample_positions[sample].add(pos_key)

        # Gene-based key (excluding synonymous)
        if effect_col and is_synonymous(row.get(effect_col, '')):
            synonymous_count += 1
            continue

        if gene_col and gene_col in df.columns:
            gene = row[gene_col]
            if not pd.isna(gene) and str(gene).strip():
                sample_genes[sample].add(str(gene).strip())
            else:
                no_gene_count += 1
        else:
            # Fallback to position if no gene column
            sample_genes[sample].add(pos_key)

    print(f"  Excluded {synonymous_count} synonymous mutations")
    print(f"  {no_gene_count} mutations had no gene annotation")

    return sample_genes, sample_positions


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
    """
    Extract treatment name from sample name.

    Examples:
        'BmKn1' -> 'BmKn'
        'BmKn2' -> 'BmKn'
        'Mel1' -> 'Mel'
        'Puro6' -> 'Puro'
        'BmKn1_1' -> 'BmKn'
    """
    # First handle underscore format (e.g., 'BmKn1_1' -> 'BmKn1')
    if '_' in sample_name:
        sample_name = sample_name.rsplit('_', 1)[0]

    # Now extract treatment by removing trailing digits
    # 'BmKn1' -> 'BmKn', 'Mel2' -> 'Mel'
    match = re.match(r'^([A-Za-z]+)', sample_name)
    if match:
        return match.group(1)
    return sample_name


def calculate_treatment_similarity(similarity_matrix: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Calculate average similarity within and between treatment groups.

    Within-treatment: similarity between replicates of the same treatment
                     (e.g., BmKn1 vs BmKn2, both are BmKn treatment)
    Between-treatment: similarity between samples of different treatments
                      (e.g., BmKn1 vs Mel1)

    Note: Self-comparisons (sample vs itself) are excluded.
    """
    samples = similarity_matrix.index.tolist()
    sample_to_treatment = {s: extract_treatment_from_sample(s) for s in samples}
    treatments = sorted(set(sample_to_treatment.values()))

    print(f"  Detected treatments: {treatments}")
    for t in treatments:
        samples_in_t = [s for s in samples if sample_to_treatment[s] == t]
        print(f"    {t}: {len(samples_in_t)} samples ({', '.join(samples_in_t[:3])}...)")

    treatment_samples = {t: [s for s in samples if sample_to_treatment[s] == t] for t in treatments}

    n_treatments = len(treatments)
    treatment_sim = np.zeros((n_treatments, n_treatments))
    stats_data = []

    for i, t1 in enumerate(treatments):
        for j, t2 in enumerate(treatments):
            samples1 = treatment_samples[t1]
            samples2 = treatment_samples[t2]

            similarities = []

            if t1 == t2:
                # Within-treatment: compare different replicates of same treatment
                # e.g., BmKn1 vs BmKn2, BmKn1 vs BmKn3, etc.
                for idx1, s1 in enumerate(samples1):
                    for idx2, s2 in enumerate(samples1):
                        if idx1 < idx2:  # Avoid self-comparison and duplicates
                            similarities.append(similarity_matrix.loc[s1, s2])
            else:
                # Between-treatment: compare all samples across treatments
                for s1 in samples1:
                    for s2 in samples2:
                        similarities.append(similarity_matrix.loc[s1, s2])

            if similarities:
                mean_sim = np.mean(similarities)
                std_sim = np.std(similarities)
            else:
                mean_sim = 0.0
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

    # Use the standard EVO color map
    row_colors = [EVO_COLOR_MAP.get(extract_treatment_from_sample(s), '#888888') for s in samples]

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

    # Create legend with EVO colors
    legend_handles = [plt.Rectangle((0, 0), 1, 1, facecolor=EVO_COLOR_MAP.get(t, '#888888'),
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
    """Create a heatmap of treatment-level similarities with EVO colors."""
    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib.patches import Patch

    fig, ax = plt.subplots(figsize=(10, 8))

    treatments = treatment_sim.index.tolist()

    sns.heatmap(
        treatment_sim,
        annot=True,
        fmt='.3f',
        cmap='RdYlBu_r',
        vmin=0, vmax=1,
        square=True,
        linewidths=2,
        cbar_kws={'label': 'Dice Similarity', 'shrink': 0.8},
        ax=ax,
        annot_kws={'fontsize': 12, 'fontweight': 'bold'}
    )

    ax.set_title("Treatment-Level Mutation Similarity\n(Mean Dice coefficient across replicates)",
                 fontsize=14, fontweight='bold', pad=15)
    ax.set_xlabel('Treatment', fontsize=12)
    ax.set_ylabel('Treatment', fontsize=12)

    # Color the tick labels using EVO colors
    for i, label in enumerate(ax.get_xticklabels()):
        treatment = treatments[i]
        label.set_color(EVO_COLOR_MAP.get(treatment, 'black'))
        label.set_fontweight('bold')

    for i, label in enumerate(ax.get_yticklabels()):
        treatment = treatments[i]
        label.set_color(EVO_COLOR_MAP.get(treatment, 'black'))
        label.set_fontweight('bold')

    plt.setp(ax.get_xticklabels(), rotation=45, ha='right', fontsize=12)
    plt.setp(ax.get_yticklabels(), rotation=0, fontsize=12)

    # Add color legend
    legend_patches = [Patch(facecolor=EVO_COLOR_MAP.get(t, '#888888'), edgecolor='black', label=t)
                      for t in treatments]
    ax.legend(handles=legend_patches, loc='upper left', bbox_to_anchor=(1.25, 1.0),
              title='Treatment', fontsize=10)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  Saved treatment heatmap to {output_path}")


def create_sample_network(similarity_matrix: pd.DataFrame, output_path: Path, min_edge_weight: float = 0.3):
    """
    Create a network visualization of sample (replicate) similarities.

    - Higher similarity = nodes closer together
    - Edge thickness = similarity strength
    - Node color = treatment (using EVO_COLOR_MAP)
    - Shows individual replicates
    """
    import matplotlib.pyplot as plt
    import networkx as nx

    samples = similarity_matrix.index.tolist()
    G = nx.Graph()

    # Add nodes with treatment info
    for s in samples:
        treatment = extract_treatment_from_sample(s)
        G.add_node(s, treatment=treatment)

    # Add edges with similarity as weight
    for i, s1 in enumerate(samples):
        for j, s2 in enumerate(samples):
            if i < j:
                sim = similarity_matrix.loc[s1, s2]
                if sim >= min_edge_weight:
                    G.add_edge(s1, s2, weight=sim, distance=1.0 - sim + 0.1)

    # Use spring layout with distance-based positioning
    if G.edges():
        pos = nx.spring_layout(
            G,
            k=1.5,
            iterations=200,
            seed=42,
            weight='distance'
        )
    else:
        pos = nx.spring_layout(G, k=2, iterations=50, seed=42)

    fig, ax = plt.subplots(figsize=(16, 14))

    # Node colors based on treatment
    node_colors = [EVO_COLOR_MAP.get(G.nodes[n]['treatment'], '#888888') for n in G.nodes()]

    # Draw nodes
    nx.draw_networkx_nodes(
        G, pos,
        node_size=800,
        node_color=node_colors,
        edgecolors='black',
        linewidths=1.5,
        ax=ax
    )

    # Draw edges
    edges = G.edges(data=True)
    if edges:
        edge_weights = [d['weight'] for _, _, d in edges]
        max_weight = max(edge_weights) if edge_weights else 1

        for (u, v, d) in edges:
            weight = d['weight']
            # Check if same treatment (within) or different (between)
            t1 = G.nodes[u]['treatment']
            t2 = G.nodes[v]['treatment']

            width = 0.5 + 4 * (weight / max_weight)
            alpha = 0.2 + 0.6 * (weight / max_weight)

            if t1 == t2:
                # Within-treatment edges: use treatment color
                color = EVO_COLOR_MAP.get(t1, '#888888')
            else:
                # Between-treatment edges: gray
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

    # Draw labels
    nx.draw_networkx_labels(G, pos, font_size=8, font_weight='bold', ax=ax)

    ax.set_title("Sample-Level Mutation Similarity Network\n(Gene-based, excluding synonymous)",
                 fontsize=16, fontweight='bold', pad=20)

    # Legend for treatments
    from matplotlib.patches import Patch
    treatments = sorted(set(extract_treatment_from_sample(s) for s in samples))
    legend_patches = [Patch(facecolor=EVO_COLOR_MAP.get(t, '#888888'),
                            edgecolor='black', label=t) for t in treatments]
    ax.legend(handles=legend_patches, loc='upper left', fontsize=10, title='Treatment')

    ax.text(0.02, 0.02,
            f'Showing edges with similarity >= {min_edge_weight}\n'
            'Colored edges = within-treatment\n'
            'Gray edges = between-treatment\n'
            'Distance = inverse of similarity (closer = more similar)',
            transform=ax.transAxes, fontsize=9, verticalalignment='bottom',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    ax.axis('off')
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  Saved sample network to {output_path}")


def create_similarity_network(treatment_sim: pd.DataFrame, output_path: Path, min_edge_weight: float = 0.1):
    """
    Create a network visualization of treatment similarities (aggregated).

    - Higher similarity = nodes closer together
    - Edge thickness = similarity strength
    - Node color = treatment (using EVO_COLOR_MAP)
    - No edge labels (cleaner visualization)
    """
    import matplotlib.pyplot as plt
    import networkx as nx

    treatments = treatment_sim.index.tolist()
    G = nx.Graph()

    # Add nodes
    for t in treatments:
        within_sim = treatment_sim.loc[t, t]
        G.add_node(t, within_similarity=within_sim)

    # Add edges with similarity as weight
    for i, t1 in enumerate(treatments):
        for j, t2 in enumerate(treatments):
            if i < j:
                sim = treatment_sim.loc[t1, t2]
                if sim >= min_edge_weight:
                    G.add_edge(t1, t2, weight=sim, distance=1.0 - sim + 0.1)

    # Use spring layout with distance-based positioning
    if G.edges():
        pos = nx.spring_layout(
            G,
            k=2,
            iterations=100,
            seed=42,
            weight='distance'
        )
    else:
        pos = nx.spring_layout(G, k=2, iterations=50, seed=42)

    fig, ax = plt.subplots(figsize=(12, 10))

    # Node sizes based on within-treatment similarity
    node_sizes = [2000 + 1500 * G.nodes[n].get('within_similarity', 0.5) for n in G.nodes()]

    # Use EVO color map
    node_colors = [EVO_COLOR_MAP.get(n, '#888888') for n in G.nodes()]

    # Draw nodes
    nx.draw_networkx_nodes(
        G, pos,
        node_size=node_sizes,
        node_color=node_colors,
        edgecolors='black',
        linewidths=2,
        ax=ax
    )

    # Draw edges with varying thickness based on similarity
    edges = G.edges(data=True)
    if edges:
        edge_weights = [d['weight'] for _, _, d in edges]
        max_weight = max(edge_weights) if edge_weights else 1

        for (u, v, d) in edges:
            weight = d['weight']
            width = 1 + 10 * (weight / max_weight)
            alpha = 0.3 + 0.6 * (weight / max_weight)

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

    # Draw labels
    nx.draw_networkx_labels(G, pos, font_size=14, font_weight='bold', ax=ax)

    ax.set_title("Treatment-Level Mutation Similarity Network\n(Aggregated across replicates)",
                 fontsize=16, fontweight='bold', pad=20)

    legend_elements = [
        plt.Line2D([0], [0], color='#2E7D32', linewidth=5, label='High similarity (>0.5)'),
        plt.Line2D([0], [0], color='#FF8F00', linewidth=3, label='Medium similarity (0.3-0.5)'),
        plt.Line2D([0], [0], color='#BDBDBD', linewidth=2, label='Low similarity (<0.3)')
    ]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=11)

    ax.text(0.02, 0.02,
            'Node size = within-treatment replicate consistency\n'
            'Edge thickness = between-treatment similarity\n'
            'Distance = inverse of similarity (closer = more similar)',
            transform=ax.transAxes, fontsize=10, verticalalignment='bottom',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    ax.axis('off')
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  Saved treatment network to {output_path}")


def main():
    if len(sys.argv) < 3:
        print("Usage: python run_similarity_analysis.py <input_file> <output_dir>")
        print()
        print("Input file: Long format CSV with 'sample', 'gene_name'/'locus_tag', 'effect' columns")
        print()
        print("Similarity is calculated based on mutations in the SAME GENE")
        print("(excluding synonymous mutations)")
        print()
        print("Example:")
        print("  python scripts/run_similarity_analysis.py \\")
        print("      results/summary_freebayes/all_mutations_freebayes_filtered.csv \\")
        print("      similarity_results")
        sys.exit(1)

    input_file = Path(sys.argv[1])
    output_dir = Path(sys.argv[2])

    if not input_file.exists():
        print(f"Error: Input file not found: {input_file}")
        sys.exit(1)

    output_dir.mkdir(parents=True, exist_ok=True)

    print("="*70)
    print("MUTATION PROFILE SIMILARITY ANALYSIS")
    print("="*70)
    print(f"Input: {input_file}")
    print(f"Output: {output_dir}")
    print()
    print("Similarity method: Gene-based (mutations in same gene = match)")
    print("                   Excluding synonymous mutations")
    print()

    # Load data - get gene-based mutation sets
    print("Loading mutations...")
    sample_genes, sample_positions = load_long_format(input_file, use_genes=True)

    print(f"\n  Found {len(sample_genes)} samples")
    print("\n  Sample mutation counts (non-synonymous genes):")
    for sample, genes in sorted(sample_genes.items())[:6]:
        print(f"    {sample}: {len(genes)} genes mutated")
    if len(sample_genes) > 6:
        print(f"    ... and {len(sample_genes) - 6} more samples")

    # Calculate pairwise similarities using GENES
    print("\n" + "="*70)
    print("CALCULATING GENE-BASED DICE SIMILARITIES")
    print("="*70)
    sim_matrix = calculate_pairwise_similarity(sample_genes)

    sim_matrix_path = output_dir / 'sample_similarity_matrix.csv'
    sim_matrix.to_csv(sim_matrix_path)
    print(f"  Saved sample similarity matrix to {sim_matrix_path}")

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
    print("\n" + "="*70)
    print("CREATING VISUALIZATIONS")
    print("="*70)

    # 1. Sample-level clustermap
    print("\n1. Sample-level clustered heatmap...")
    clustermap_path = output_dir / 'sample_similarity_clustermap.png'
    create_similarity_clustermap(sim_matrix, clustermap_path)

    # 2. Treatment-level heatmap
    print("\n2. Treatment-level heatmap...")
    treatment_heatmap_path = output_dir / 'treatment_similarity_heatmap.png'
    create_treatment_heatmap(treatment_sim, treatment_heatmap_path)

    # 3. Sample-level network (replicate level)
    print("\n3. Sample-level (replicate) network...")
    sample_network_path = output_dir / 'sample_similarity_network.png'
    create_sample_network(sim_matrix, sample_network_path, min_edge_weight=0.3)

    # 4. Treatment-level network (aggregated)
    print("\n4. Treatment-level (aggregated) network...")
    treatment_network_path = output_dir / 'treatment_similarity_network.png'
    create_similarity_network(treatment_sim, treatment_network_path, min_edge_weight=0.1)

    # Print summary
    print()
    print("="*70)
    print("SIMILARITY ANALYSIS SUMMARY")
    print("="*70)

    within_stats = treatment_stats[treatment_stats['type'] == 'within']
    print("\nWithin-treatment similarity (replicate consistency):")
    print("  (How consistently do replicates of the same AMP mutate the same genes?)")
    for _, row in within_stats.sort_values('mean_similarity', ascending=False).iterrows():
        n_comp = int(row['n_comparisons'])
        print(f"  {row['treatment1']:8s}: {row['mean_similarity']:.3f} ± {row['std_similarity']:.3f} (n={n_comp} pairs)")

    between_stats = treatment_stats[treatment_stats['type'] == 'between'].copy()
    between_stats = between_stats.sort_values('mean_similarity', ascending=False)

    print("\nMost similar treatment pairs (target similar genes/pathways):")
    for _, row in between_stats.head(5).iterrows():
        print(f"  {row['treatment1']:8s} - {row['treatment2']:8s}: {row['mean_similarity']:.3f}")

    print("\nLeast similar treatment pairs (target different genes/pathways):")
    for _, row in between_stats.tail(3).iterrows():
        print(f"  {row['treatment1']:8s} - {row['treatment2']:8s}: {row['mean_similarity']:.3f}")

    print()
    print("="*70)
    print("OUTPUT FILES")
    print("="*70)
    print(f"  {output_dir}/sample_similarity_matrix.csv      - Pairwise sample similarities")
    print(f"  {output_dir}/treatment_similarity_matrix.csv   - Treatment-level similarities")
    print(f"  {output_dir}/treatment_similarity_stats.csv    - Within/between statistics")
    print(f"  {output_dir}/sample_similarity_clustermap.png  - Clustered heatmap (samples)")
    print(f"  {output_dir}/treatment_similarity_heatmap.png  - Treatment heatmap")
    print(f"  {output_dir}/sample_similarity_network.png     - Replicate-level network")
    print(f"  {output_dir}/treatment_similarity_network.png  - Treatment-level network")


if __name__ == '__main__':
    main()
