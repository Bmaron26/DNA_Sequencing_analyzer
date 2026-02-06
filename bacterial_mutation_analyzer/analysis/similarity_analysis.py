"""
Mutation profile similarity analysis using Dice coefficient.

Analyzes similarity between mutation profiles across samples and treatments
to identify which AMPs induce similar mutations/target similar pathways.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from itertools import combinations
import warnings


def dice_similarity(set1: set, set2: set) -> float:
    """
    Calculate Dice similarity coefficient between two sets.

    Dice = 2|A∩B| / (|A| + |B|)

    Returns 0 if both sets are empty, 1 if identical.
    """
    if len(set1) == 0 and len(set2) == 0:
        return 1.0  # Both empty = identical
    if len(set1) == 0 or len(set2) == 0:
        return 0.0  # One empty, one not = no similarity

    intersection = len(set1 & set2)
    return (2 * intersection) / (len(set1) + len(set2))


def load_mutation_matrix(filepath: Path) -> Tuple[pd.DataFrame, Dict[str, set]]:
    """
    Load mutation matrix and extract mutation sets per sample.

    Returns:
        df: The original dataframe
        sample_mutations: Dict mapping sample name to set of variant keys
    """
    df = pd.read_csv(filepath)

    # Identify sample columns (exclude metadata columns)
    metadata_cols = ['position', 'ref', 'alt', 'gene', 'product', 'effect',
                     'variant_key', 'chromosome', 'CHROM', 'POS', 'REF', 'ALT',
                     'GENE', 'PRODUCT', 'EFFECT', 'locus_tag', 'feature_type',
                     'strand', 'gene_id', 'gene_name', 'variant_type', 'quality',
                     'depth', 'allele_frequency', 'reference', 'alternative']

    sample_cols = [c for c in df.columns if c.lower() not in [m.lower() for m in metadata_cols]]

    # Create variant key if not present
    if 'variant_key' not in df.columns:
        # Try to create from position/ref/alt
        if 'position' in df.columns:
            pos_col, ref_col, alt_col = 'position', 'ref', 'alt'
        elif 'POS' in df.columns:
            pos_col, ref_col, alt_col = 'POS', 'REF', 'ALT'
        else:
            # Use row index as key
            df['variant_key'] = df.index.astype(str)
            pos_col = None

        if pos_col:
            df['variant_key'] = df.apply(
                lambda r: f"{r[pos_col]}_{r.get(ref_col, '')}_{r.get(alt_col, '')}",
                axis=1
            )

    # Extract mutation sets per sample
    sample_mutations = {}
    for sample in sample_cols:
        # A mutation is present if the value is 1 (or truthy)
        mask = df[sample].fillna(0).astype(bool)
        mutations = set(df.loc[mask, 'variant_key'].tolist())
        sample_mutations[sample] = mutations

    return df, sample_mutations


def load_long_format(filepath: Path) -> Dict[str, set]:
    """
    Load mutations from long format CSV and extract mutation sets per sample.
    """
    df = pd.read_csv(filepath)

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
        raise ValueError("No sample column found in long format file")

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
    """
    Calculate pairwise Dice similarity between all samples.

    Returns:
        DataFrame with samples as both index and columns, values are Dice similarities
    """
    samples = sorted(sample_mutations.keys())
    n = len(samples)

    # Initialize similarity matrix
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

    Assumes format like 'BmKn1_1', 'Nisin_2', etc.
    Returns the part before the last underscore.
    """
    parts = sample_name.rsplit('_', 1)
    if len(parts) == 2 and parts[1].isdigit():
        return parts[0]
    return sample_name


def calculate_treatment_similarity(
    similarity_matrix: pd.DataFrame,
    sample_to_treatment: Optional[Dict[str, str]] = None
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Calculate average similarity within and between treatment groups.

    Returns:
        treatment_sim: Mean similarity between treatment pairs
        treatment_stats: DataFrame with within/between treatment statistics
    """
    samples = similarity_matrix.index.tolist()

    # Auto-detect treatment mapping if not provided
    if sample_to_treatment is None:
        sample_to_treatment = {s: extract_treatment_from_sample(s) for s in samples}

    treatments = sorted(set(sample_to_treatment.values()))

    # Group samples by treatment
    treatment_samples = {t: [s for s in samples if sample_to_treatment[s] == t]
                         for t in treatments}

    # Calculate mean similarity between each treatment pair
    n_treatments = len(treatments)
    treatment_sim = np.zeros((n_treatments, n_treatments))

    stats_data = []

    for i, t1 in enumerate(treatments):
        for j, t2 in enumerate(treatments):
            samples1 = treatment_samples[t1]
            samples2 = treatment_samples[t2]

            # Get all pairwise similarities between the groups
            similarities = []
            for s1 in samples1:
                for s2 in samples2:
                    if s1 != s2:  # Exclude self-comparisons
                        similarities.append(similarity_matrix.loc[s1, s2])

            if similarities:
                mean_sim = np.mean(similarities)
                std_sim = np.std(similarities)
            else:
                mean_sim = 1.0 if t1 == t2 else 0.0
                std_sim = 0.0

            treatment_sim[i, j] = mean_sim

            if i <= j:  # Store stats for unique pairs
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


def create_similarity_clustermap(
    similarity_matrix: pd.DataFrame,
    output_path: Path,
    title: str = "Mutation Profile Similarity (Dice Coefficient)",
    figsize: Tuple[int, int] = (14, 12),
    cmap: str = "RdYlBu_r"
) -> None:
    """
    Create a clustermap visualization of sample similarities.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns

    # Extract treatment for color coding
    samples = similarity_matrix.index.tolist()
    treatments = [extract_treatment_from_sample(s) for s in samples]
    unique_treatments = sorted(set(treatments))

    # Create color palette for treatments
    colors = plt.cm.Set2(np.linspace(0, 1, len(unique_treatments)))
    treatment_colors = {t: colors[i] for i, t in enumerate(unique_treatments)}
    row_colors = [treatment_colors[extract_treatment_from_sample(s)] for s in samples]

    # Create clustermap
    g = sns.clustermap(
        similarity_matrix,
        method='average',
        metric='euclidean',
        cmap=cmap,
        vmin=0, vmax=1,
        row_colors=row_colors,
        col_colors=row_colors,
        figsize=figsize,
        dendrogram_ratio=(0.15, 0.15),
        cbar_pos=(0.02, 0.8, 0.03, 0.15),
        linewidths=0.5,
        xticklabels=True,
        yticklabels=True
    )

    g.fig.suptitle(title, fontsize=14, fontweight='bold', y=1.02)

    # Add legend for treatments
    legend_handles = [plt.Rectangle((0, 0), 1, 1, facecolor=treatment_colors[t],
                                     edgecolor='black', linewidth=0.5)
                      for t in unique_treatments]
    g.ax_heatmap.legend(legend_handles, unique_treatments,
                        loc='upper left', bbox_to_anchor=(1.15, 1.0),
                        title='Treatment', frameon=True, fontsize=9)

    # Rotate labels
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right', fontsize=9)
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=9)

    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"  Saved clustermap to {output_path}")


def create_treatment_heatmap(
    treatment_sim: pd.DataFrame,
    output_path: Path,
    title: str = "Treatment-Level Mutation Similarity",
    figsize: Tuple[int, int] = (10, 8),
    cmap: str = "RdYlBu_r"
) -> None:
    """
    Create a heatmap of treatment-level similarities.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns

    fig, ax = plt.subplots(figsize=figsize)

    # Create heatmap with annotations
    sns.heatmap(
        treatment_sim,
        annot=True,
        fmt='.3f',
        cmap=cmap,
        vmin=0, vmax=1,
        square=True,
        linewidths=1,
        cbar_kws={'label': 'Dice Similarity', 'shrink': 0.8},
        ax=ax,
        annot_kws={'fontsize': 11, 'fontweight': 'bold'}
    )

    ax.set_title(title, fontsize=14, fontweight='bold', pad=15)
    ax.set_xlabel('Treatment', fontsize=12)
    ax.set_ylabel('Treatment', fontsize=12)

    plt.setp(ax.get_xticklabels(), rotation=45, ha='right', fontsize=11)
    plt.setp(ax.get_yticklabels(), rotation=0, fontsize=11)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"  Saved treatment heatmap to {output_path}")


def create_similarity_network(
    treatment_sim: pd.DataFrame,
    output_path: Path,
    title: str = "AMP Mutation Similarity Network",
    min_edge_weight: float = 0.1,
    figsize: Tuple[int, int] = (12, 10)
) -> None:
    """
    Create a network visualization of treatment similarities.

    Nodes = treatments, edges = similarity (weighted by Dice coefficient)
    """
    import matplotlib.pyplot as plt
    import networkx as nx

    treatments = treatment_sim.index.tolist()

    # Create graph
    G = nx.Graph()

    # Add nodes
    for t in treatments:
        # Within-treatment similarity as node attribute
        within_sim = treatment_sim.loc[t, t]
        G.add_node(t, within_similarity=within_sim)

    # Add edges (between-treatment similarities)
    for i, t1 in enumerate(treatments):
        for j, t2 in enumerate(treatments):
            if i < j:
                sim = treatment_sim.loc[t1, t2]
                if sim >= min_edge_weight:
                    G.add_edge(t1, t2, weight=sim)

    # Layout
    pos = nx.spring_layout(G, k=2, iterations=50, seed=42)

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Node sizes based on within-treatment similarity
    node_sizes = [1500 + 1000 * G.nodes[n].get('within_similarity', 0.5) for n in G.nodes()]

    # Node colors
    colors = plt.cm.Set2(np.linspace(0, 1, len(treatments)))
    node_colors = {t: colors[i] for i, t in enumerate(treatments)}

    # Draw nodes
    nx.draw_networkx_nodes(
        G, pos,
        node_size=node_sizes,
        node_color=[node_colors[n] for n in G.nodes()],
        edgecolors='black',
        linewidths=2,
        ax=ax
    )

    # Draw edges with width proportional to similarity
    edges = G.edges(data=True)
    if edges:
        edge_weights = [d['weight'] for _, _, d in edges]
        max_weight = max(edge_weights) if edge_weights else 1

        # Draw edges with varying thickness and color
        for (u, v, d) in edges:
            weight = d['weight']
            width = 1 + 8 * (weight / max_weight)
            alpha = 0.3 + 0.7 * (weight / max_weight)

            # Color based on similarity strength
            if weight > 0.5:
                color = 'darkgreen'
            elif weight > 0.3:
                color = 'orange'
            else:
                color = 'lightgray'

            nx.draw_networkx_edges(
                G, pos,
                edgelist=[(u, v)],
                width=width,
                alpha=alpha,
                edge_color=color,
                ax=ax
            )

        # Add edge labels
        edge_labels = {(u, v): f"{d['weight']:.2f}" for u, v, d in edges if d['weight'] > 0.15}
        nx.draw_networkx_edge_labels(
            G, pos,
            edge_labels=edge_labels,
            font_size=9,
            font_weight='bold',
            ax=ax
        )

    # Draw labels
    nx.draw_networkx_labels(
        G, pos,
        font_size=12,
        font_weight='bold',
        ax=ax
    )

    ax.set_title(title, fontsize=14, fontweight='bold', pad=15)

    # Add legend
    legend_elements = [
        plt.Line2D([0], [0], color='darkgreen', linewidth=4, label='High similarity (>0.5)'),
        plt.Line2D([0], [0], color='orange', linewidth=3, label='Medium similarity (0.3-0.5)'),
        plt.Line2D([0], [0], color='lightgray', linewidth=2, label='Low similarity (<0.3)')
    ]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=10)

    # Add note about node size
    ax.text(0.02, 0.02, 'Node size = within-treatment similarity\nEdge thickness = between-treatment similarity',
            transform=ax.transAxes, fontsize=9, verticalalignment='bottom',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    ax.axis('off')
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"  Saved similarity network to {output_path}")


def run_similarity_analysis(
    input_file: Path,
    output_dir: Path,
    file_format: str = 'auto',
    min_edge_weight: float = 0.1
) -> Dict:
    """
    Run complete similarity analysis pipeline.

    Args:
        input_file: Path to mutation matrix (matrix format) or all_mutations.csv (long format)
        output_dir: Directory to save results
        file_format: 'matrix', 'long', or 'auto' (auto-detect)
        min_edge_weight: Minimum similarity to show edge in network

    Returns:
        Dictionary with analysis results
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading mutations from {input_file}...")

    # Auto-detect format
    if file_format == 'auto':
        df_check = pd.read_csv(input_file, nrows=5)
        if 'sample' in df_check.columns.str.lower().tolist():
            file_format = 'long'
        else:
            file_format = 'matrix'
        print(f"  Detected format: {file_format}")

    # Load data
    if file_format == 'matrix':
        _, sample_mutations = load_mutation_matrix(input_file)
    else:
        sample_mutations = load_long_format(input_file)

    print(f"  Found {len(sample_mutations)} samples")
    for sample, muts in list(sample_mutations.items())[:3]:
        print(f"    {sample}: {len(muts)} mutations")

    # Calculate pairwise similarities
    print("\nCalculating pairwise Dice similarities...")
    sim_matrix = calculate_pairwise_similarity(sample_mutations)

    # Save similarity matrix
    sim_matrix_path = output_dir / 'sample_similarity_matrix.csv'
    sim_matrix.to_csv(sim_matrix_path)
    print(f"  Saved similarity matrix to {sim_matrix_path}")

    # Calculate treatment-level similarities
    print("\nCalculating treatment-level similarities...")
    treatment_sim, treatment_stats = calculate_treatment_similarity(sim_matrix)

    # Save treatment results
    treatment_sim_path = output_dir / 'treatment_similarity_matrix.csv'
    treatment_sim.to_csv(treatment_sim_path)
    print(f"  Saved treatment similarity matrix to {treatment_sim_path}")

    treatment_stats_path = output_dir / 'treatment_similarity_stats.csv'
    treatment_stats.to_csv(treatment_stats_path, index=False)
    print(f"  Saved treatment statistics to {treatment_stats_path}")

    # Create visualizations
    print("\nCreating visualizations...")

    # Sample-level clustermap
    clustermap_path = output_dir / 'sample_similarity_clustermap.png'
    create_similarity_clustermap(sim_matrix, clustermap_path)

    # Treatment-level heatmap
    treatment_heatmap_path = output_dir / 'treatment_similarity_heatmap.png'
    create_treatment_heatmap(treatment_sim, treatment_heatmap_path)

    # Similarity network
    network_path = output_dir / 'treatment_similarity_network.png'
    create_similarity_network(treatment_sim, network_path, min_edge_weight=min_edge_weight)

    # Summary statistics
    print("\n" + "="*60)
    print("SIMILARITY ANALYSIS SUMMARY")
    print("="*60)

    # Within-treatment similarities (replicate consistency)
    within_stats = treatment_stats[treatment_stats['type'] == 'within']
    print("\nWithin-treatment similarity (replicate consistency):")
    for _, row in within_stats.iterrows():
        print(f"  {row['treatment1']}: {row['mean_similarity']:.3f} ± {row['std_similarity']:.3f}")

    # Most similar treatment pairs
    between_stats = treatment_stats[treatment_stats['type'] == 'between'].copy()
    between_stats = between_stats.sort_values('mean_similarity', ascending=False)

    print("\nMost similar treatment pairs:")
    for _, row in between_stats.head(5).iterrows():
        print(f"  {row['treatment1']} - {row['treatment2']}: {row['mean_similarity']:.3f}")

    print("\nLeast similar treatment pairs:")
    for _, row in between_stats.tail(3).iterrows():
        print(f"  {row['treatment1']} - {row['treatment2']}: {row['mean_similarity']:.3f}")

    return {
        'sample_similarity': sim_matrix,
        'treatment_similarity': treatment_sim,
        'treatment_stats': treatment_stats,
        'sample_mutations': sample_mutations
    }


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Analyze mutation profile similarity between samples and treatments'
    )
    parser.add_argument('input_file', help='Mutation matrix CSV or long-format mutations CSV')
    parser.add_argument('output_dir', help='Output directory for results')
    parser.add_argument('--format', choices=['auto', 'matrix', 'long'], default='auto',
                        help='Input file format (default: auto-detect)')
    parser.add_argument('--min-edge-weight', type=float, default=0.1,
                        help='Minimum similarity for network edges (default: 0.1)')

    args = parser.parse_args()

    run_similarity_analysis(
        Path(args.input_file),
        Path(args.output_dir),
        file_format=args.format,
        min_edge_weight=args.min_edge_weight
    )
