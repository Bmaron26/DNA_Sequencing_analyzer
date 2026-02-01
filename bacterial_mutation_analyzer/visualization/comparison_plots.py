"""
Visualization module for multi-sample comparison results.

Generates:
- Convergent mutation heatmaps
- Group comparison plots
- Mutation overlap Venn diagrams
- Treatment-specific mutation profiles
"""

import os
from pathlib import Path
from typing import Dict, Any, List, Optional
from collections import defaultdict
import logging

logger = logging.getLogger(__name__)

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    import matplotlib.colors as mcolors
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False

try:
    import seaborn as sns
    SEABORN_AVAILABLE = True
except ImportError:
    SEABORN_AVAILABLE = False

try:
    import numpy as np
    NUMPY_AVAILABLE = True
except ImportError:
    NUMPY_AVAILABLE = False

try:
    import plotly.graph_objects as go
    import plotly.express as px
    from plotly.subplots import make_subplots
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False


class ComparisonVisualizer:
    """
    Generates visualizations for multi-sample comparison results.
    """

    # Color palettes
    GROUP_COLORS = [
        '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
        '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'
    ]

    IMPACT_COLORS = {
        'HIGH': '#d62728',
        'MODERATE': '#ff7f0e',
        'LOW': '#2ca02c',
        'MODIFIER': '#1f77b4',
    }

    def __init__(self, output_dir: str = "comparison_plots",
                 format: str = "png", dpi: int = 150):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.format = format
        self.dpi = dpi

        if SEABORN_AVAILABLE:
            sns.set_style("whitegrid")

    def generate_all(self, comparison_data: Dict[str, Any],
                    experiment_data: Optional[Dict] = None) -> List[str]:
        """
        Generate all comparison visualizations.

        Args:
            comparison_data: Output from MultiSampleComparison.generate_comparison_report()
            experiment_data: Optional experiment metadata

        Returns:
            List of generated file paths
        """
        generated = []

        if not MATPLOTLIB_AVAILABLE:
            logger.warning("Matplotlib not available, skipping visualizations")
            return generated

        # Mutation presence heatmap
        try:
            path = self.plot_mutation_heatmap(comparison_data)
            if path:
                generated.append(path)
        except Exception as e:
            logger.warning(f"Could not generate mutation heatmap: {e}")

        # Convergent mutations bar plot
        try:
            path = self.plot_convergent_mutations(comparison_data)
            if path:
                generated.append(path)
        except Exception as e:
            logger.warning(f"Could not generate convergent mutations plot: {e}")

        # Group comparison
        try:
            path = self.plot_group_comparison(comparison_data)
            if path:
                generated.append(path)
        except Exception as e:
            logger.warning(f"Could not generate group comparison: {e}")

        # Effect distribution by group
        try:
            path = self.plot_effect_by_group(comparison_data)
            if path:
                generated.append(path)
        except Exception as e:
            logger.warning(f"Could not generate effect distribution: {e}")

        # Interactive plots
        if PLOTLY_AVAILABLE:
            try:
                path = self.plot_interactive_comparison(comparison_data)
                if path:
                    generated.append(path)
            except Exception as e:
                logger.warning(f"Could not generate interactive comparison: {e}")

        return generated

    def plot_mutation_heatmap(self, data: Dict[str, Any]) -> Optional[str]:
        """
        Create heatmap of mutation presence across samples.
        """
        if not MATPLOTLIB_AVAILABLE or not NUMPY_AVAILABLE:
            return None

        convergent = data.get('convergent_mutations', [])
        if not convergent:
            return None

        # Collect samples and mutations
        all_samples = set()
        for cm in convergent:
            all_samples.update(cm.get('samples', []))

        samples = sorted(all_samples)
        mutations = []
        matrix = []

        for cm in convergent[:50]:  # Limit to top 50
            gene = cm.get('gene_name') or cm.get('locus_tag') or ''
            aa_change = cm.get('amino_acid_change', '')
            label = f"{gene} {aa_change}" if aa_change else gene

            mutations.append(label)

            row = []
            freqs = cm.get('frequencies', {})
            for sample in samples:
                freq = freqs.get(sample, 0)
                row.append(freq if freq > 0 else 0)
            matrix.append(row)

        if not matrix:
            return None

        matrix = np.array(matrix)

        # Create figure
        fig_height = max(8, len(mutations) * 0.3)
        fig_width = max(10, len(samples) * 0.5)

        fig, ax = plt.subplots(figsize=(fig_width, fig_height))

        if SEABORN_AVAILABLE:
            sns.heatmap(matrix, xticklabels=samples, yticklabels=mutations,
                       cmap='YlOrRd', ax=ax, cbar_kws={'label': 'Allele Frequency'})
        else:
            im = ax.imshow(matrix, cmap='YlOrRd', aspect='auto')
            ax.set_xticks(range(len(samples)))
            ax.set_xticklabels(samples, rotation=45, ha='right')
            ax.set_yticks(range(len(mutations)))
            ax.set_yticklabels(mutations)
            plt.colorbar(im, ax=ax, label='Allele Frequency')

        ax.set_title('Mutation Presence Across Samples')
        ax.set_xlabel('Sample')
        ax.set_ylabel('Mutation')

        plt.tight_layout()

        output_path = self.output_dir / f"mutation_heatmap.{self.format}"
        plt.savefig(output_path, dpi=self.dpi, bbox_inches='tight')
        plt.close()

        return str(output_path)

    def plot_convergent_mutations(self, data: Dict[str, Any]) -> Optional[str]:
        """
        Bar plot of convergent mutations by sample count and gene.
        """
        if not MATPLOTLIB_AVAILABLE:
            return None

        convergent = data.get('convergent_mutations', [])
        if not convergent:
            return None

        # Group by gene
        gene_counts = defaultdict(int)
        gene_samples = defaultdict(int)

        for cm in convergent:
            gene = cm.get('gene_name') or cm.get('locus_tag') or 'Unknown'
            gene_counts[gene] += 1
            gene_samples[gene] = max(gene_samples[gene], cm.get('sample_count', 0))

        # Sort by sample count then mutation count
        genes = sorted(gene_counts.keys(),
                      key=lambda g: (gene_samples[g], gene_counts[g]),
                      reverse=True)[:20]

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        # Mutations per gene
        counts = [gene_counts[g] for g in genes]
        colors = plt.cm.viridis([gene_samples[g] / max(gene_samples.values())
                                for g in genes])

        bars = ax1.barh(genes[::-1], counts[::-1], color=colors[::-1])
        ax1.set_xlabel('Number of Convergent Mutations')
        ax1.set_ylabel('Gene')
        ax1.set_title('Convergent Mutations by Gene')

        # Add colorbar for sample count
        sm = plt.cm.ScalarMappable(cmap='viridis',
                                   norm=plt.Normalize(1, max(gene_samples.values())))
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax1)
        cbar.set_label('Max Samples')

        # Sample count distribution
        sample_counts = [cm.get('sample_count', 0) for cm in convergent]
        ax2.hist(sample_counts, bins=range(1, max(sample_counts) + 2),
                edgecolor='white', alpha=0.7)
        ax2.set_xlabel('Number of Samples')
        ax2.set_ylabel('Number of Mutations')
        ax2.set_title('Distribution of Convergent Mutations by Sample Count')

        plt.tight_layout()

        output_path = self.output_dir / f"convergent_mutations.{self.format}"
        plt.savefig(output_path, dpi=self.dpi, bbox_inches='tight')
        plt.close()

        return str(output_path)

    def plot_group_comparison(self, data: Dict[str, Any]) -> Optional[str]:
        """
        Compare mutation counts and characteristics between groups.
        """
        if not MATPLOTLIB_AVAILABLE:
            return None

        summaries = data.get('group_summaries', {})
        if not summaries:
            return None

        groups = list(summaries.keys())
        n_groups = len(groups)

        if n_groups < 2:
            return None

        fig, axes = plt.subplots(2, 2, figsize=(12, 10))

        # 1. Total mutations per group
        ax1 = axes[0, 0]
        total_muts = [summaries[g].get('total_unique_mutations', 0) for g in groups]
        colors = [self.GROUP_COLORS[i % len(self.GROUP_COLORS)] for i in range(n_groups)]
        ax1.bar(groups, total_muts, color=colors)
        ax1.set_ylabel('Total Unique Mutations')
        ax1.set_title('Mutations by Group')
        ax1.tick_params(axis='x', rotation=45)

        # 2. Convergence rate per group
        ax2 = axes[0, 1]
        conv_rates = [summaries[g].get('convergence_rate', 0) * 100 for g in groups]
        ax2.bar(groups, conv_rates, color=colors)
        ax2.set_ylabel('Convergence Rate (%)')
        ax2.set_title('Within-Group Convergence Rate')
        ax2.tick_params(axis='x', rotation=45)

        # 3. Genes mutated per group
        ax3 = axes[1, 0]
        genes_mutated = [summaries[g].get('genes_mutated', 0) for g in groups]
        ax3.bar(groups, genes_mutated, color=colors)
        ax3.set_ylabel('Genes Mutated')
        ax3.set_title('Number of Genes with Mutations')
        ax3.tick_params(axis='x', rotation=45)

        # 4. Top genes comparison
        ax4 = axes[1, 1]

        # Get union of top genes across groups
        top_genes = set()
        for g in groups:
            top = list(summaries[g].get('top_mutated_genes', {}).keys())[:5]
            top_genes.update(top)

        top_genes = list(top_genes)[:10]

        if top_genes:
            x = np.arange(len(top_genes))
            width = 0.8 / n_groups

            for i, group in enumerate(groups):
                gene_counts = summaries[group].get('top_mutated_genes', {})
                values = [gene_counts.get(g, 0) for g in top_genes]
                ax4.bar(x + i * width, values, width, label=group,
                       color=self.GROUP_COLORS[i % len(self.GROUP_COLORS)])

            ax4.set_xticks(x + width * (n_groups - 1) / 2)
            ax4.set_xticklabels(top_genes, rotation=45, ha='right')
            ax4.set_ylabel('Mutation Count')
            ax4.set_title('Top Mutated Genes by Group')
            ax4.legend()

        plt.tight_layout()

        output_path = self.output_dir / f"group_comparison.{self.format}"
        plt.savefig(output_path, dpi=self.dpi, bbox_inches='tight')
        plt.close()

        return str(output_path)

    def plot_effect_by_group(self, data: Dict[str, Any]) -> Optional[str]:
        """
        Compare mutation effects between groups.
        """
        if not MATPLOTLIB_AVAILABLE:
            return None

        summaries = data.get('group_summaries', {})
        if not summaries:
            return None

        groups = list(summaries.keys())

        # Collect all effects
        all_effects = set()
        for g in groups:
            effects = summaries[g].get('effect_distribution', {})
            all_effects.update(effects.keys())

        if not all_effects:
            return None

        # Limit to top effects
        effect_totals = defaultdict(int)
        for g in groups:
            effects = summaries[g].get('effect_distribution', {})
            for e, c in effects.items():
                effect_totals[e] += c

        top_effects = sorted(effect_totals.keys(),
                            key=lambda e: effect_totals[e], reverse=True)[:10]

        fig, ax = plt.subplots(figsize=(12, 6))

        x = np.arange(len(top_effects))
        width = 0.8 / len(groups)

        for i, group in enumerate(groups):
            effects = summaries[group].get('effect_distribution', {})
            values = [effects.get(e, 0) for e in top_effects]
            ax.bar(x + i * width, values, width, label=group,
                  color=self.GROUP_COLORS[i % len(self.GROUP_COLORS)])

        ax.set_xticks(x + width * (len(groups) - 1) / 2)
        ax.set_xticklabels([e[:25] for e in top_effects], rotation=45, ha='right')
        ax.set_ylabel('Count')
        ax.set_title('Mutation Effects by Group')
        ax.legend()

        plt.tight_layout()

        output_path = self.output_dir / f"effect_by_group.{self.format}"
        plt.savefig(output_path, dpi=self.dpi, bbox_inches='tight')
        plt.close()

        return str(output_path)

    def plot_interactive_comparison(self, data: Dict[str, Any]) -> Optional[str]:
        """
        Create interactive comparison dashboard.
        """
        if not PLOTLY_AVAILABLE:
            return None

        summaries = data.get('group_summaries', {})
        convergent = data.get('convergent_mutations', [])

        if not summaries:
            return None

        groups = list(summaries.keys())

        fig = make_subplots(
            rows=2, cols=2,
            specs=[[{"type": "bar"}, {"type": "bar"}],
                   [{"type": "scatter"}, {"type": "heatmap"}]],
            subplot_titles=('Mutations by Group', 'Convergence Rate',
                          'Convergent Mutations', 'Group Similarity')
        )

        colors = self.GROUP_COLORS[:len(groups)]

        # 1. Total mutations
        total_muts = [summaries[g].get('total_unique_mutations', 0) for g in groups]
        fig.add_trace(go.Bar(x=groups, y=total_muts, marker_color=colors,
                            name='Total Mutations'), row=1, col=1)

        # 2. Convergence rate
        conv_rates = [summaries[g].get('convergence_rate', 0) * 100 for g in groups]
        fig.add_trace(go.Bar(x=groups, y=conv_rates, marker_color=colors,
                            name='Convergence %'), row=1, col=2)

        # 3. Convergent mutations scatter
        if convergent:
            for cm in convergent[:100]:
                groups_affected = list(set(cm.get('groups', [])))
                samples = cm.get('samples', [])
                gene = cm.get('gene_name') or cm.get('locus_tag') or ''

                fig.add_trace(go.Scatter(
                    x=[cm.get('sample_count', 0)],
                    y=[cm.get('convergence_score', 0)],
                    mode='markers',
                    marker=dict(size=10, opacity=0.6),
                    text=f"{gene}<br>Samples: {', '.join(samples[:5])}",
                    hoverinfo='text',
                    showlegend=False,
                ), row=2, col=1)

        # 4. Group similarity matrix (based on shared mutations)
        if 'pairwise_comparisons' in data:
            comparisons = data['pairwise_comparisons']
            similarity_matrix = np.zeros((len(groups), len(groups)))

            for comp in comparisons:
                g1, g2 = comp['groups']
                shared = comp['shared_count']
                total = shared + comp['unique_to_first'] + comp['unique_to_second']
                similarity = shared / max(total, 1)

                i1 = groups.index(g1)
                i2 = groups.index(g2)
                similarity_matrix[i1, i2] = similarity
                similarity_matrix[i2, i1] = similarity

            np.fill_diagonal(similarity_matrix, 1)

            fig.add_trace(go.Heatmap(
                z=similarity_matrix,
                x=groups,
                y=groups,
                colorscale='Blues',
                showscale=True,
            ), row=2, col=2)

        fig.update_layout(
            title='Multi-Sample Comparison Dashboard',
            height=800,
            showlegend=False,
        )

        output_path = self.output_dir / "comparison_dashboard.html"
        fig.write_html(str(output_path))

        return str(output_path)
