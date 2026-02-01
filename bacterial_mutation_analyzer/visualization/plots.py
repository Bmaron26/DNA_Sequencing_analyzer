"""
Visualization module for generating plots from analysis results.

Generates:
- Mutation spectrum plot
- Variant distribution across genome
- Coverage plot
- Quality metrics
- Interactive HTML reports
"""

import os
from pathlib import Path
from typing import Dict, Any, List, Optional
from collections import defaultdict
import logging

logger = logging.getLogger(__name__)

try:
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False

try:
    import seaborn as sns
    SEABORN_AVAILABLE = True
except ImportError:
    SEABORN_AVAILABLE = False

try:
    import plotly.graph_objects as go
    import plotly.express as px
    from plotly.subplots import make_subplots
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False

try:
    import numpy as np
    NUMPY_AVAILABLE = True
except ImportError:
    NUMPY_AVAILABLE = False


class VisualizationGenerator:
    """
    Generates visualizations from pipeline results.
    """

    # Mutation type colors
    MUTATION_COLORS = {
        'C>A': '#1E88E5',  # Blue
        'C>G': '#000000',  # Black
        'C>T': '#D81B60',  # Red
        'T>A': '#757575',  # Gray
        'T>C': '#43A047',  # Green
        'T>G': '#F4511E',  # Orange
        'A>C': '#757575',
        'A>G': '#43A047',
        'A>T': '#1E88E5',
        'G>A': '#D81B60',
        'G>C': '#000000',
        'G>T': '#F4511E',
    }

    # Effect impact colors
    IMPACT_COLORS = {
        'HIGH': '#D32F2F',
        'MODERATE': '#FFA000',
        'LOW': '#388E3C',
        'MODIFIER': '#1976D2',
    }

    def __init__(self, output_dir: str = "visualizations",
                 format: str = "png",
                 dpi: int = 150):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.format = format
        self.dpi = dpi

        if SEABORN_AVAILABLE:
            sns.set_style("whitegrid")
            sns.set_context("paper", font_scale=1.2)

    def generate_all(self, result: Dict[str, Any]) -> List[str]:
        """
        Generate all available visualizations.

        Args:
            result: Pipeline result dictionary

        Returns:
            List of generated file paths
        """
        generated = []

        variants = result.get('variants', [])

        if not variants:
            logger.warning("No variants to visualize")
            return generated

        # Generate static plots (matplotlib/seaborn)
        if MATPLOTLIB_AVAILABLE:
            try:
                path = self.plot_mutation_spectrum(variants)
                if path:
                    generated.append(path)
            except Exception as e:
                logger.warning(f"Could not generate mutation spectrum: {e}")

            try:
                path = self.plot_variant_types(variants)
                if path:
                    generated.append(path)
            except Exception as e:
                logger.warning(f"Could not generate variant types plot: {e}")

            try:
                path = self.plot_effect_distribution(variants)
                if path:
                    generated.append(path)
            except Exception as e:
                logger.warning(f"Could not generate effect distribution: {e}")

            try:
                path = self.plot_chromosome_distribution(variants)
                if path:
                    generated.append(path)
            except Exception as e:
                logger.warning(f"Could not generate chromosome distribution: {e}")

            try:
                path = self.plot_quality_metrics(result)
                if path:
                    generated.append(path)
            except Exception as e:
                logger.warning(f"Could not generate quality metrics: {e}")

        # Generate interactive plots (plotly)
        if PLOTLY_AVAILABLE:
            try:
                path = self.plot_interactive_genome_view(variants)
                if path:
                    generated.append(path)
            except Exception as e:
                logger.warning(f"Could not generate genome view: {e}")

            try:
                path = self.plot_interactive_summary(result)
                if path:
                    generated.append(path)
            except Exception as e:
                logger.warning(f"Could not generate interactive summary: {e}")

        return generated

    def plot_mutation_spectrum(self, variants: List[Dict]) -> Optional[str]:
        """
        Plot mutation spectrum (6 or 96 class).

        Shows the distribution of single nucleotide substitution types.
        """
        if not MATPLOTLIB_AVAILABLE:
            return None

        # Count mutation types
        snps = [v for v in variants if v.get('variant_type') == 'SNP']
        if not snps:
            return None

        mutation_counts = defaultdict(int)
        for v in snps:
            ref = v.get('reference', '').upper()
            alt = v.get('alternative', '').upper()
            if len(ref) == 1 and len(alt) == 1:
                # Normalize to pyrimidine reference
                if ref in 'AG':
                    ref = {'A': 'T', 'G': 'C'}[ref]
                    alt = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}[alt]
                mutation_counts[f"{ref}>{alt}"] += 1

        # Create plot
        fig, ax = plt.subplots(figsize=(10, 6))

        categories = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
        counts = [mutation_counts.get(cat, 0) for cat in categories]
        colors = [self.MUTATION_COLORS.get(cat, '#888888') for cat in categories]

        bars = ax.bar(categories, counts, color=colors, edgecolor='white', linewidth=1)

        ax.set_xlabel('Mutation Type')
        ax.set_ylabel('Count')
        ax.set_title('Mutation Spectrum')

        # Add count labels
        for bar, count in zip(bars, counts):
            if count > 0:
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                       str(count), ha='center', va='bottom', fontsize=10)

        plt.tight_layout()

        output_path = self.output_dir / f"mutation_spectrum.{self.format}"
        plt.savefig(output_path, dpi=self.dpi, bbox_inches='tight')
        plt.close()

        return str(output_path)

    def plot_variant_types(self, variants: List[Dict]) -> Optional[str]:
        """Plot distribution of variant types (SNP, insertion, deletion, etc.)."""
        if not MATPLOTLIB_AVAILABLE:
            return None

        type_counts = defaultdict(int)
        for v in variants:
            vtype = v.get('variant_type', 'unknown')
            type_counts[vtype] += 1

        if not type_counts:
            return None

        fig, ax = plt.subplots(figsize=(8, 6))

        types = list(type_counts.keys())
        counts = list(type_counts.values())

        colors = plt.cm.Set2(range(len(types)))
        wedges, texts, autotexts = ax.pie(
            counts,
            labels=types,
            autopct='%1.1f%%',
            colors=colors,
            explode=[0.02] * len(types),
            shadow=False,
            startangle=90
        )

        ax.set_title('Variant Type Distribution')

        # Add legend with counts
        legend_labels = [f"{t}: {c}" for t, c in zip(types, counts)]
        ax.legend(wedges, legend_labels, title="Types",
                 loc="center left", bbox_to_anchor=(1, 0, 0.5, 1))

        plt.tight_layout()

        output_path = self.output_dir / f"variant_types.{self.format}"
        plt.savefig(output_path, dpi=self.dpi, bbox_inches='tight')
        plt.close()

        return str(output_path)

    def plot_effect_distribution(self, variants: List[Dict]) -> Optional[str]:
        """Plot distribution of variant effects and their impact."""
        if not MATPLOTLIB_AVAILABLE:
            return None

        # Count by effect
        effect_counts = defaultdict(int)
        impact_counts = defaultdict(int)

        for v in variants:
            effect = v.get('effect', 'unknown')
            impact = v.get('effect_impact', 'unknown')
            effect_counts[effect] += 1
            impact_counts[impact] += 1

        if not effect_counts:
            return None

        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        # Effect distribution
        ax1 = axes[0]
        effects = sorted(effect_counts.items(), key=lambda x: x[1], reverse=True)[:10]
        effect_names = [e[0] for e in effects]
        effect_values = [e[1] for e in effects]

        bars = ax1.barh(effect_names[::-1], effect_values[::-1], color='steelblue')
        ax1.set_xlabel('Count')
        ax1.set_title('Top 10 Variant Effects')

        # Add count labels
        for bar, count in zip(bars, effect_values[::-1]):
            ax1.text(bar.get_width() + 0.5, bar.get_y() + bar.get_height()/2,
                    str(count), ha='left', va='center', fontsize=9)

        # Impact distribution
        ax2 = axes[1]
        impact_order = ['HIGH', 'MODERATE', 'LOW', 'MODIFIER']
        impacts = [i for i in impact_order if i in impact_counts]
        impact_values = [impact_counts[i] for i in impacts]
        colors = [self.IMPACT_COLORS.get(i, '#888888') for i in impacts]

        bars = ax2.bar(impacts, impact_values, color=colors, edgecolor='white')
        ax2.set_xlabel('Impact')
        ax2.set_ylabel('Count')
        ax2.set_title('Variant Impact Distribution')

        for bar, count in zip(bars, impact_values):
            if count > 0:
                ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                        str(count), ha='center', va='bottom', fontsize=10)

        plt.tight_layout()

        output_path = self.output_dir / f"effect_distribution.{self.format}"
        plt.savefig(output_path, dpi=self.dpi, bbox_inches='tight')
        plt.close()

        return str(output_path)

    def plot_chromosome_distribution(self, variants: List[Dict]) -> Optional[str]:
        """Plot variant distribution across chromosomes/contigs."""
        if not MATPLOTLIB_AVAILABLE or not NUMPY_AVAILABLE:
            return None

        # Group by chromosome
        by_chrom = defaultdict(list)
        for v in variants:
            chrom = v.get('chromosome', 'unknown')
            pos = v.get('position', 0)
            by_chrom[chrom].append(pos)

        if not by_chrom:
            return None

        # Sort chromosomes
        chroms = sorted(by_chrom.keys())

        fig, ax = plt.subplots(figsize=(12, max(6, len(chroms) * 0.5)))

        y_positions = []
        for i, chrom in enumerate(chroms):
            positions = by_chrom[chrom]
            y = [i] * len(positions)
            y_positions.append(i)
            ax.scatter(positions, y, alpha=0.6, s=10, label=chrom)

        ax.set_yticks(y_positions)
        ax.set_yticklabels(chroms)
        ax.set_xlabel('Position (bp)')
        ax.set_ylabel('Chromosome/Contig')
        ax.set_title('Variant Distribution Across Genome')

        plt.tight_layout()

        output_path = self.output_dir / f"chromosome_distribution.{self.format}"
        plt.savefig(output_path, dpi=self.dpi, bbox_inches='tight')
        plt.close()

        return str(output_path)

    def plot_quality_metrics(self, result: Dict[str, Any]) -> Optional[str]:
        """Plot quality metrics summary."""
        if not MATPLOTLIB_AVAILABLE:
            return None

        fig, axes = plt.subplots(2, 2, figsize=(12, 10))

        # QC metrics
        if result.get('qc_stats'):
            qc = result['qc_stats']
            ax = axes[0, 0]

            metrics = ['Q20 Rate', 'Q30 Rate', 'GC Content', 'Read Retention']
            values = [
                qc.get('q20_rate', 0) * 100,
                qc.get('q30_rate', 0) * 100,
                qc.get('gc_content', 0) * 100,
                qc.get('read_retention_rate', 0) * 100,
            ]

            colors = ['#2ecc71' if v >= 80 else '#f39c12' if v >= 60 else '#e74c3c'
                     for v in values]

            bars = ax.bar(metrics, values, color=colors)
            ax.axhline(y=80, color='gray', linestyle='--', alpha=0.5)
            ax.set_ylabel('Percentage (%)')
            ax.set_title('Quality Control Metrics')
            ax.set_ylim(0, 105)

            for bar, val in zip(bars, values):
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                       f'{val:.1f}%', ha='center', va='bottom', fontsize=9)
        else:
            axes[0, 0].text(0.5, 0.5, 'No QC data', ha='center', va='center',
                           transform=axes[0, 0].transAxes)

        # Alignment metrics
        if result.get('alignment_stats'):
            align = result['alignment_stats']
            ax = axes[0, 1]

            metrics = ['Mapping\nRate', 'Coverage\nBreadth', 'Uniformity']
            values = [
                align.get('mapping_rate', 0) * 100,
                align.get('coverage_breadth', 0) * 100,
                align.get('coverage_uniformity', 0) * 100,
            ]

            colors = ['#2ecc71' if v >= 90 else '#f39c12' if v >= 70 else '#e74c3c'
                     for v in values]

            bars = ax.bar(metrics, values, color=colors)
            ax.axhline(y=90, color='gray', linestyle='--', alpha=0.5)
            ax.set_ylabel('Percentage (%)')
            ax.set_title('Alignment Metrics')
            ax.set_ylim(0, 105)

            for bar, val in zip(bars, values):
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                       f'{val:.1f}%', ha='center', va='bottom', fontsize=9)
        else:
            axes[0, 1].text(0.5, 0.5, 'No alignment data', ha='center', va='center',
                           transform=axes[0, 1].transAxes)

        # Variant counts
        if result.get('variant_stats'):
            var = result['variant_stats']
            ax = axes[1, 0]

            categories = ['SNPs', 'Insertions', 'Deletions', 'Complex']
            counts = [
                var.get('snps', 0),
                var.get('insertions', 0),
                var.get('deletions', 0),
                var.get('complex_variants', 0) + var.get('mnps', 0),
            ]

            colors = plt.cm.Set2(range(4))
            bars = ax.bar(categories, counts, color=colors)
            ax.set_ylabel('Count')
            ax.set_title('Variant Counts by Type')

            for bar, count in zip(bars, counts):
                if count > 0:
                    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                           str(count), ha='center', va='bottom', fontsize=10)
        else:
            axes[1, 0].text(0.5, 0.5, 'No variant data', ha='center', va='center',
                           transform=axes[1, 0].transAxes)

        # Coverage histogram (if available)
        ax = axes[1, 1]
        if result.get('alignment_stats'):
            align = result['alignment_stats']
            mean_cov = align.get('mean_coverage', 0)
            median_cov = align.get('median_coverage', 0)

            ax.text(0.5, 0.6, f'Mean Coverage: {mean_cov:.1f}x',
                   ha='center', va='center', fontsize=14,
                   transform=ax.transAxes)
            ax.text(0.5, 0.4, f'Median Coverage: {median_cov:.1f}x',
                   ha='center', va='center', fontsize=14,
                   transform=ax.transAxes)
            ax.set_title('Coverage Summary')
            ax.axis('off')
        else:
            ax.text(0.5, 0.5, 'No coverage data', ha='center', va='center',
                   transform=ax.transAxes)
            ax.axis('off')

        plt.tight_layout()

        output_path = self.output_dir / f"quality_metrics.{self.format}"
        plt.savefig(output_path, dpi=self.dpi, bbox_inches='tight')
        plt.close()

        return str(output_path)

    def plot_interactive_genome_view(self, variants: List[Dict]) -> Optional[str]:
        """Create interactive genome-wide variant view using Plotly."""
        if not PLOTLY_AVAILABLE:
            return None

        # Prepare data
        chroms = []
        positions = []
        types = []
        effects = []
        hover_texts = []

        for v in variants:
            chroms.append(v.get('chromosome', 'unknown'))
            positions.append(v.get('position', 0))
            types.append(v.get('variant_type', 'unknown'))
            effects.append(v.get('effect', 'unknown'))

            hover = (f"Position: {v.get('chromosome', '')}:{v.get('position', '')}<br>"
                    f"Change: {v.get('reference', '')}>{v.get('alternative', '')}<br>"
                    f"Type: {v.get('variant_type', '')}<br>"
                    f"Effect: {v.get('effect', '')}<br>"
                    f"Gene: {v.get('gene_name', '') or v.get('locus_tag', 'N/A')}")
            hover_texts.append(hover)

        fig = go.Figure()

        # Add traces for each variant type
        for vtype in set(types):
            mask = [t == vtype for t in types]
            fig.add_trace(go.Scatter(
                x=[p for p, m in zip(positions, mask) if m],
                y=[c for c, m in zip(chroms, mask) if m],
                mode='markers',
                name=vtype,
                text=[h for h, m in zip(hover_texts, mask) if m],
                hoverinfo='text',
                marker=dict(size=8, opacity=0.7)
            ))

        fig.update_layout(
            title='Interactive Genome-Wide Variant View',
            xaxis_title='Position (bp)',
            yaxis_title='Chromosome/Contig',
            hovermode='closest',
            template='plotly_white',
            height=600,
            legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1)
        )

        output_path = self.output_dir / "genome_view.html"
        fig.write_html(str(output_path))

        return str(output_path)

    def plot_interactive_summary(self, result: Dict[str, Any]) -> Optional[str]:
        """Create interactive summary dashboard using Plotly."""
        if not PLOTLY_AVAILABLE:
            return None

        fig = make_subplots(
            rows=2, cols=2,
            specs=[[{"type": "domain"}, {"type": "bar"}],
                   [{"type": "bar"}, {"type": "bar"}]],
            subplot_titles=('Variant Types', 'Mutation Spectrum',
                          'Effect Impact', 'Quality Metrics')
        )

        variants = result.get('variants', [])

        # Pie chart for variant types
        type_counts = defaultdict(int)
        for v in variants:
            type_counts[v.get('variant_type', 'unknown')] += 1

        fig.add_trace(go.Pie(
            labels=list(type_counts.keys()),
            values=list(type_counts.values()),
            hole=0.4,
            name='Types'
        ), row=1, col=1)

        # Mutation spectrum
        snps = [v for v in variants if v.get('variant_type') == 'SNP']
        mut_counts = defaultdict(int)
        for v in snps:
            ref = v.get('reference', '').upper()
            alt = v.get('alternative', '').upper()
            if len(ref) == 1 and len(alt) == 1:
                if ref in 'AG':
                    ref = {'A': 'T', 'G': 'C'}[ref]
                    alt = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}[alt]
                mut_counts[f"{ref}>{alt}"] += 1

        categories = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
        colors = [self.MUTATION_COLORS.get(c, '#888888') for c in categories]

        fig.add_trace(go.Bar(
            x=categories,
            y=[mut_counts.get(c, 0) for c in categories],
            marker_color=colors,
            name='Mutations'
        ), row=1, col=2)

        # Effect impact
        impact_counts = defaultdict(int)
        for v in variants:
            impact_counts[v.get('effect_impact', 'unknown')] += 1

        impact_order = ['HIGH', 'MODERATE', 'LOW', 'MODIFIER']
        impact_colors = [self.IMPACT_COLORS.get(i, '#888888') for i in impact_order]

        fig.add_trace(go.Bar(
            x=impact_order,
            y=[impact_counts.get(i, 0) for i in impact_order],
            marker_color=impact_colors,
            name='Impact'
        ), row=2, col=1)

        # Quality metrics
        metrics = []
        values = []

        if result.get('qc_stats'):
            qc = result['qc_stats']
            metrics.extend(['Q30 Rate', 'GC Content'])
            values.extend([qc.get('q30_rate', 0) * 100, qc.get('gc_content', 0) * 100])

        if result.get('alignment_stats'):
            align = result['alignment_stats']
            metrics.extend(['Mapping Rate', 'Coverage'])
            values.extend([align.get('mapping_rate', 0) * 100,
                          min(align.get('coverage_breadth', 0) * 100, 100)])

        colors = ['#2ecc71' if v >= 80 else '#f39c12' if v >= 60 else '#e74c3c'
                 for v in values]

        fig.add_trace(go.Bar(
            x=metrics,
            y=values,
            marker_color=colors,
            name='Quality'
        ), row=2, col=2)

        fig.update_layout(
            title='Analysis Summary Dashboard',
            showlegend=False,
            height=700,
            template='plotly_white'
        )

        output_path = self.output_dir / "summary_dashboard.html"
        fig.write_html(str(output_path))

        return str(output_path)

    def generate_html_report(self, result: Dict[str, Any],
                            plots: List[str]) -> str:
        """Generate a comprehensive HTML report."""
        try:
            from jinja2 import Template
        except ImportError:
            logger.warning("Jinja2 not available for HTML report generation")
            return ""

        template = Template(HTML_REPORT_TEMPLATE)

        html_content = template.render(
            sample_name=result.get('sample_name', 'Unknown'),
            start_time=result.get('start_time', ''),
            end_time=result.get('end_time', ''),
            qc_stats=result.get('qc_stats', {}),
            alignment_stats=result.get('alignment_stats', {}),
            variant_stats=result.get('variant_stats', {}),
            annotation_summary=result.get('annotation_summary', {}),
            variants=result.get('variants', [])[:100],  # Limit to 100 for display
            total_variants=len(result.get('variants', [])),
            warnings=result.get('warnings', []),
            outliers=result.get('outliers', []),
            plots=[Path(p).name for p in plots if p],
        )

        output_path = self.output_dir / "analysis_report.html"
        with open(output_path, 'w') as f:
            f.write(html_content)

        return str(output_path)


HTML_REPORT_TEMPLATE = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Mutation Analysis Report - {{ sample_name }}</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; background: #f5f5f5; }
        .container { max-width: 1200px; margin: 0 auto; background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }
        h1 { color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }
        h2 { color: #34495e; margin-top: 30px; }
        .stats-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; margin: 20px 0; }
        .stat-card { background: #ecf0f1; padding: 15px; border-radius: 5px; text-align: center; }
        .stat-value { font-size: 24px; font-weight: bold; color: #2980b9; }
        .stat-label { color: #7f8c8d; font-size: 14px; }
        table { width: 100%; border-collapse: collapse; margin: 20px 0; }
        th, td { padding: 10px; text-align: left; border-bottom: 1px solid #ddd; }
        th { background: #3498db; color: white; }
        tr:hover { background: #f5f5f5; }
        .warning { background: #fff3cd; border-left: 4px solid #ffc107; padding: 10px; margin: 10px 0; }
        .high-impact { color: #c0392b; font-weight: bold; }
        .moderate-impact { color: #e67e22; }
        .low-impact { color: #27ae60; }
        .plots { display: grid; grid-template-columns: repeat(auto-fit, minmax(400px, 1fr)); gap: 20px; }
        .plot img { max-width: 100%; border: 1px solid #ddd; border-radius: 5px; }
    </style>
</head>
<body>
    <div class="container">
        <h1>Mutation Analysis Report</h1>
        <p><strong>Sample:</strong> {{ sample_name }}</p>
        <p><strong>Analysis Date:</strong> {{ start_time }}</p>

        {% if warnings %}
        <h2>Warnings</h2>
        {% for warning in warnings %}
        <div class="warning">{{ warning }}</div>
        {% endfor %}
        {% endif %}

        <h2>Quality Control</h2>
        <div class="stats-grid">
            <div class="stat-card">
                <div class="stat-value">{{ "{:,.0f}".format(qc_stats.get('total_reads', 0)) }}</div>
                <div class="stat-label">Total Reads</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{{ "{:.1%}".format(qc_stats.get('q30_rate', 0)) }}</div>
                <div class="stat-label">Q30 Rate</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{{ "{:.1%}".format(qc_stats.get('gc_content', 0)) }}</div>
                <div class="stat-label">GC Content</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{{ "{:.1f}".format(qc_stats.get('mean_length', 0)) }}</div>
                <div class="stat-label">Mean Read Length</div>
            </div>
        </div>

        <h2>Alignment</h2>
        <div class="stats-grid">
            <div class="stat-card">
                <div class="stat-value">{{ "{:.1%}".format(alignment_stats.get('mapping_rate', 0)) }}</div>
                <div class="stat-label">Mapping Rate</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{{ "{:.1f}x".format(alignment_stats.get('mean_coverage', 0)) }}</div>
                <div class="stat-label">Mean Coverage</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{{ "{:.1%}".format(alignment_stats.get('coverage_breadth', 0)) }}</div>
                <div class="stat-label">Coverage Breadth</div>
            </div>
        </div>

        <h2>Variants</h2>
        <div class="stats-grid">
            <div class="stat-card">
                <div class="stat-value">{{ variant_stats.get('total_variants', 0) }}</div>
                <div class="stat-label">Total Variants</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{{ variant_stats.get('snps', 0) }}</div>
                <div class="stat-label">SNPs</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{{ variant_stats.get('insertions', 0) }}</div>
                <div class="stat-label">Insertions</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{{ variant_stats.get('deletions', 0) }}</div>
                <div class="stat-label">Deletions</div>
            </div>
        </div>

        {% if annotation_summary %}
        <h2>Annotation Summary</h2>
        <div class="stats-grid">
            <div class="stat-card">
                <div class="stat-value">{{ annotation_summary.get('genes_affected', [])|length }}</div>
                <div class="stat-label">Genes Affected</div>
            </div>
            <div class="stat-card">
                <div class="stat-value high-impact">{{ annotation_summary.get('by_impact', {}).get('HIGH', 0) }}</div>
                <div class="stat-label">High Impact</div>
            </div>
            <div class="stat-card">
                <div class="stat-value moderate-impact">{{ annotation_summary.get('by_impact', {}).get('MODERATE', 0) }}</div>
                <div class="stat-label">Moderate Impact</div>
            </div>
            <div class="stat-card">
                <div class="stat-value low-impact">{{ annotation_summary.get('by_impact', {}).get('LOW', 0) }}</div>
                <div class="stat-label">Low Impact</div>
            </div>
        </div>
        {% endif %}

        <h2>Mutations (showing first 100 of {{ total_variants }})</h2>
        <table>
            <tr>
                <th>Position</th>
                <th>Change</th>
                <th>Type</th>
                <th>Gene</th>
                <th>Effect</th>
                <th>AA Change</th>
                <th>Frequency</th>
            </tr>
            {% for v in variants %}
            <tr>
                <td>{{ v.chromosome }}:{{ v.position }}</td>
                <td>{{ v.reference }}>{{ v.alternative }}</td>
                <td>{{ v.variant_type }}</td>
                <td>{{ v.gene_name or v.locus_tag or '-' }}</td>
                <td class="{{ 'high-impact' if v.effect_impact == 'HIGH' else 'moderate-impact' if v.effect_impact == 'MODERATE' else '' }}">{{ v.effect }}</td>
                <td>{{ v.amino_acid_change or '-' }}</td>
                <td>{{ "{:.1%}".format(v.allele_frequency) if v.allele_frequency else '-' }}</td>
            </tr>
            {% endfor %}
        </table>

        {% if plots %}
        <h2>Visualizations</h2>
        <div class="plots">
            {% for plot in plots %}
            <div class="plot">
                <img src="{{ plot }}" alt="{{ plot }}">
            </div>
            {% endfor %}
        </div>
        {% endif %}
    </div>
</body>
</html>
"""
