"""
Circular genome plots for bacterial mutation visualization.
Creates Circos-style plots showing mutations across the genome.
"""

import os
import json
import math
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Wedge, FancyArrowPatch
from matplotlib.collections import PatchCollection
import matplotlib.colors as mcolors


@dataclass
class GenomeFeature:
    """Represents a genomic feature for plotting."""
    name: str
    start: int
    end: int
    strand: str = '+'
    feature_type: str = 'gene'


@dataclass
class MutationMarker:
    """Represents a mutation to plot on the circle."""
    position: int
    sample: str
    variant_type: str  # SNP, insertion, deletion
    effect: str  # synonymous, nonsynonymous, intergenic, etc.
    gene: str = ''
    amino_acid_change: str = ''
    allele_frequency: float = 1.0
    locus_tag: str = ''


def extract_short_gene_name(gene: str, product: str = '', locus_tag: str = '') -> str:
    """
    Extract a short gene name from gene/product/locus_tag fields.

    Prioritizes:
    1. Standard gene names (e.g., agrC, clpX, mprF)
    2. Gene name from product description
    3. Short locus tag
    4. Truncated product name as last resort

    Args:
        gene: Gene name field
        product: Product description
        locus_tag: Locus tag

    Returns:
        Short gene name suitable for plotting
    """
    import re

    # Common gene name patterns (case-insensitive)
    gene_pattern = re.compile(r'\b([a-zA-Z]{2,4}[A-Z0-9]?)\b')

    # 1. If gene field looks like a standard gene name (2-5 chars, starts with lowercase)
    if gene and len(gene) <= 6 and re.match(r'^[a-z]{2,4}[A-Z0-9]?$', gene):
        return gene

    # 2. Check if gene field contains a gene name (e.g., "agrC" in "accessory gene regulator C")
    if gene:
        # Look for capitalized gene-like pattern at end of string
        match = re.search(r'\b([A-Z][a-z]{1,3}[A-Z0-9]?)\s*$', gene)
        if match:
            return match.group(1)
        # Look for standard gene name pattern
        match = re.search(r'\b([a-z]{2,4}[A-Z][0-9]?)\b', gene)
        if match:
            return match.group(1)

    # 3. Extract gene name from product description
    if product:
        # Common patterns: "... protein ClpX", "... subunit AgrC", "... MprF"
        # Look for capitalized gene name at end
        match = re.search(r'\b([A-Z][a-z]{1,3}[A-Z0-9]?)\s*$', product)
        if match:
            name = match.group(1)
            if len(name) >= 3:
                return name

        # Look for gene name pattern anywhere
        match = re.search(r'\b([A-Z][a-z]{1,3}[A-Z0-9])\b', product)
        if match:
            return match.group(1)

        # Try lowercase gene names (e.g., "rpoB", "gyrA")
        match = re.search(r'\b([a-z]{3,4}[A-Z][0-9]?)\b', product)
        if match:
            return match.group(1)

    # 4. Use locus tag if short enough
    if locus_tag:
        # Extract just the number part if it's like "SAOUHSC_00123"
        match = re.search(r'_(\d+)$', locus_tag)
        if match:
            return f"#{match.group(1)}"
        if len(locus_tag) <= 12:
            return locus_tag

    # 5. Last resort: truncate product name
    if product:
        # Get first meaningful word
        words = product.split()
        for word in words:
            if len(word) >= 3 and word[0].isupper():
                return word[:10]
        return product[:10] + "..."

    if gene:
        return gene[:10] + "..." if len(gene) > 10 else gene

    return "unknown"


class CircosPlot:
    """
    Create Circos-style circular genome plots.
    """

    # Color schemes for different mutation types
    EFFECT_COLORS = {
        'missense_variant': '#E41A1C',      # Red - non-synonymous
        'nonsynonymous': '#E41A1C',
        'non-synonymous': '#E41A1C',
        'missense': '#E41A1C',
        'synonymous_variant': '#377EB8',    # Blue - synonymous
        'synonymous': '#377EB8',
        'stop_gained': '#984EA3',           # Purple - stop
        'stop_lost': '#984EA3',
        'frameshift': '#FF7F00',            # Orange - frameshift
        'frameshift_variant': '#FF7F00',
        'inframe_insertion': '#FFFF33',     # Yellow
        'inframe_deletion': '#FFFF33',
        'intergenic': '#4DAF4A',            # Green - intergenic
        'intergenic_region': '#4DAF4A',
        'upstream': '#A65628',              # Brown
        'downstream': '#A65628',
        'unknown': '#999999',               # Gray
    }

    # Treatment colors (for different AMPs)
    TREATMENT_COLORS = [
        '#E41A1C', '#377EB8', '#4DAF4A', '#984EA3',
        '#FF7F00', '#FFFF33', '#A65628', '#F781BF'
    ]

    def __init__(self, genome_size: int, title: str = ''):
        """
        Initialize CircosPlot.

        Args:
            genome_size: Size of the genome in base pairs
            title: Plot title
        """
        self.genome_size = genome_size
        self.title = title
        self.tracks = []  # List of track data
        self.genes = []   # Gene annotations

    def add_gene_track(self, genes: List[GenomeFeature],
                       color: str = '#CCCCCC', label: str = 'Genes'):
        """Add a gene annotation track."""
        self.genes = genes

    def add_mutation_track(self, mutations: List[MutationMarker],
                          sample_name: str, color: str = None,
                          ring_index: int = 0):
        """
        Add a track of mutations for a sample.

        Args:
            mutations: List of MutationMarker objects
            sample_name: Name of the sample
            color: Base color for this track (if None, use effect colors)
            ring_index: Which ring to place this track on (0 = outermost)
        """
        self.tracks.append({
            'mutations': mutations,
            'sample_name': sample_name,
            'color': color,
            'ring_index': ring_index
        })

    def _position_to_angle(self, position: int) -> float:
        """Convert genomic position to angle in radians."""
        # Start at top (90 degrees) and go clockwise
        fraction = position / self.genome_size
        angle = (0.5 - fraction) * 2 * math.pi
        return angle

    def _get_effect_color(self, effect: str) -> str:
        """Get color for mutation effect."""
        effect_lower = effect.lower() if effect else 'unknown'
        for key, color in self.EFFECT_COLORS.items():
            if key in effect_lower:
                return color
        return self.EFFECT_COLORS['unknown']

    def plot(self, output_path: str, figsize: Tuple[int, int] = (12, 12),
             dpi: int = 150, show_genes: bool = True,
             show_legend: bool = True):
        """
        Generate the circular plot.

        Args:
            output_path: Path to save the plot
            figsize: Figure size in inches
            dpi: Resolution
            show_genes: Whether to show gene labels for mutated genes
            show_legend: Whether to show the legend
        """
        fig, ax = plt.subplots(figsize=figsize, subplot_kw={'aspect': 'equal'})
        ax.set_xlim(-1.5, 1.5)
        ax.set_ylim(-1.5, 1.5)
        ax.axis('off')

        # Draw genome backbone circle
        n_tracks = len(self.tracks)
        outer_radius = 1.0
        track_width = 0.08 if n_tracks <= 6 else 0.05
        gap = 0.02

        # Draw outer circle (genome backbone)
        backbone = plt.Circle((0, 0), outer_radius + 0.02,
                              fill=False, color='black', linewidth=2)
        ax.add_patch(backbone)

        # Draw position markers (0, 1MB, 2MB, etc.)
        self._draw_position_markers(ax, outer_radius + 0.05)

        # Draw each track (ring)
        mutated_genes = set()
        for track in self.tracks:
            ring_idx = track['ring_index']
            radius = outer_radius - (ring_idx * (track_width + gap))

            # Draw ring background
            ring = plt.Circle((0, 0), radius, fill=False,
                             color='#EEEEEE', linewidth=0.5)
            ax.add_patch(ring)

            # Plot mutations
            for mut in track['mutations']:
                angle = self._position_to_angle(mut.position)
                x = radius * math.cos(angle)
                y = radius * math.sin(angle)

                # Determine color
                if track['color']:
                    color = track['color']
                else:
                    color = self._get_effect_color(mut.effect)

                # Determine marker based on effect
                if 'synonymous' in (mut.effect or '').lower() and 'non' not in (mut.effect or '').lower():
                    marker = 'o'  # Circle for synonymous
                elif 'intergenic' in (mut.effect or '').lower():
                    marker = 's'  # Square for intergenic
                else:
                    marker = '*'  # Star for non-synonymous/other

                ax.plot(x, y, marker=marker, color=color, markersize=8,
                       markeredgecolor='black', markeredgewidth=0.5)

                # Track mutated genes for labeling
                if mut.gene:
                    mutated_genes.add((mut.gene, mut.position))

        # Add gene labels for mutated genes
        if show_genes and mutated_genes:
            self._add_gene_labels(ax, mutated_genes, outer_radius + 0.15)

        # Add title
        if self.title:
            ax.set_title(self.title, fontsize=14, fontweight='bold', pad=20)

        # Add legend
        if show_legend:
            self._add_legend(ax)

        # Add sample legend (track colors)
        if len(self.tracks) > 1:
            self._add_track_legend(ax)

        plt.tight_layout()
        plt.savefig(output_path, dpi=dpi, bbox_inches='tight',
                   facecolor='white', edgecolor='none')
        plt.close()

        return output_path

    def _draw_position_markers(self, ax, radius: float):
        """Draw genome position markers (0MB, 1MB, etc.)."""
        mb_size = 1_000_000
        n_markers = int(self.genome_size / mb_size) + 1

        for i in range(n_markers):
            pos = i * mb_size
            if pos > self.genome_size:
                break
            angle = self._position_to_angle(pos)
            x = radius * math.cos(angle)
            y = radius * math.sin(angle)

            # Draw tick
            x_inner = (radius - 0.03) * math.cos(angle)
            y_inner = (radius - 0.03) * math.sin(angle)
            ax.plot([x_inner, x], [y_inner, y], 'k-', linewidth=1)

            # Add label
            label = f"{i} MB" if i > 0 else "0 MB"
            ha = 'left' if x > 0 else 'right' if x < 0 else 'center'
            va = 'bottom' if y > 0 else 'top' if y < 0 else 'center'
            ax.text(x * 1.08, y * 1.08, label, fontsize=9,
                   ha=ha, va=va)

    def _add_gene_labels(self, ax, mutated_genes: set, radius: float):
        """Add labels for mutated genes with overlap avoidance."""
        # Sort genes by position
        genes_sorted = sorted(mutated_genes, key=lambda x: x[1])

        # Group nearby genes to avoid overlap
        min_angle_diff = 0.15  # Minimum angle difference between labels (radians)
        placed_angles = []

        for gene_name, position in genes_sorted:
            if not gene_name or gene_name == 'unknown':
                continue

            angle = self._position_to_angle(position)

            # Check if too close to already placed label
            too_close = False
            for placed_angle in placed_angles:
                if abs(angle - placed_angle) < min_angle_diff:
                    too_close = True
                    break

            if too_close:
                continue  # Skip this label to avoid overlap

            placed_angles.append(angle)

            # Calculate label position
            x = radius * math.cos(angle)
            y = radius * math.sin(angle)

            # Rotate text to be readable
            rotation = math.degrees(angle) - 90
            if rotation < -90 or rotation > 90:
                rotation += 180
                ha = 'right'
            else:
                ha = 'left'

            ax.text(x, y, gene_name, fontsize=7, fontstyle='italic',
                   rotation=rotation, ha=ha, va='center')

    def _add_legend(self, ax):
        """Add mutation type legend."""
        legend_elements = [
            plt.Line2D([0], [0], marker='*', color='w',
                      markerfacecolor='#E41A1C', markersize=10,
                      markeredgecolor='black', label='Non-synonymous SNP'),
            plt.Line2D([0], [0], marker='o', color='w',
                      markerfacecolor='#377EB8', markersize=8,
                      markeredgecolor='black', label='Synonymous SNP'),
            plt.Line2D([0], [0], marker='s', color='w',
                      markerfacecolor='#4DAF4A', markersize=8,
                      markeredgecolor='black', label='Intergenic SNP'),
            plt.Line2D([0], [0], marker='*', color='w',
                      markerfacecolor='#984EA3', markersize=10,
                      markeredgecolor='black', label='Stop gained/lost'),
            plt.Line2D([0], [0], marker='*', color='w',
                      markerfacecolor='#FF7F00', markersize=10,
                      markeredgecolor='black', label='Frameshift'),
        ]
        ax.legend(handles=legend_elements, loc='lower left',
                 fontsize=8, framealpha=0.9)

    def _add_track_legend(self, ax):
        """Add legend for sample tracks."""
        legend_elements = []
        for i, track in enumerate(self.tracks):
            color = track['color'] or self.TREATMENT_COLORS[i % len(self.TREATMENT_COLORS)]
            legend_elements.append(
                plt.Line2D([0], [0], marker='o', color='w',
                          markerfacecolor=color, markersize=8,
                          label=track['sample_name'])
            )

        ax.legend(handles=legend_elements, loc='lower right',
                 fontsize=8, framealpha=0.9, title='Samples')


def create_treatment_circos_plots(results_dir: str, output_dir: str,
                                  experiment_file: Optional[str] = None,
                                  genome_size: int = 2800000,
                                  caller: str = 'bcftools') -> List[str]:
    """
    Create one Circos plot per treatment group.

    Args:
        results_dir: Directory containing sample results
        output_dir: Where to save the plots
        experiment_file: CSV/JSON with sample->group mapping
        genome_size: Genome size in bp
        caller: Which variant caller results to use

    Returns:
        List of generated plot paths
    """
    os.makedirs(output_dir, exist_ok=True)

    # Load experiment info if provided
    sample_groups = {}
    if experiment_file:
        if experiment_file.endswith('.csv'):
            df = pd.read_csv(experiment_file)
            for _, row in df.iterrows():
                sample_groups[row['sample_id']] = row.get('group', 'Unknown')
        else:
            with open(experiment_file) as f:
                exp = json.load(f)
                for sample_id, sample in exp.get('samples', {}).items():
                    sample_groups[sample_id] = sample.get('group', 'Unknown')

    # Load all sample results
    results_path = Path(results_dir)
    samples_by_group = {}

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
            # Try alternative naming
            mut_file = sample_dir / f"{sample_name}_mutations.csv"
            if not mut_file.exists():
                continue

        # Determine group
        group = sample_groups.get(sample_name, 'Unknown')
        if group not in samples_by_group:
            samples_by_group[group] = []

        # Load mutations
        df = pd.read_csv(mut_file)
        mutations = []
        for _, row in df.iterrows():
            # Extract short gene name for display
            raw_gene = row.get('GENE', row.get('gene_name', ''))
            product = row.get('PRODUCT', row.get('product', ''))
            locus_tag = row.get('LOCUS_TAG', row.get('locus_tag', ''))
            short_gene = extract_short_gene_name(raw_gene, product, locus_tag)

            mut = MutationMarker(
                position=int(row.get('POS', row.get('position', 0))),
                sample=sample_name,
                variant_type=row.get('TYPE', row.get('variant_type', 'SNP')),
                effect=row.get('EFFECT', row.get('effect', '')),
                gene=short_gene,  # Use short gene name
                amino_acid_change=row.get('AA_CHANGE', row.get('amino_acid_change', '')),
                allele_frequency=float(row.get('FREQ', row.get('allele_frequency', 1.0))),
                locus_tag=locus_tag
            )
            mutations.append(mut)

        samples_by_group[group].append({
            'name': sample_name,
            'mutations': mutations
        })

    # Create one plot per group
    plot_paths = []
    colors = ['#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#FF7F00', '#A65628']

    for group_name, samples in samples_by_group.items():
        if not samples:
            continue

        plot = CircosPlot(genome_size=genome_size,
                         title=f"{group_name} - Mutations ({caller})")

        for i, sample in enumerate(samples):
            plot.add_mutation_track(
                sample['mutations'],
                sample_name=sample['name'],
                color=colors[i % len(colors)],
                ring_index=i
            )

        output_path = os.path.join(output_dir, f"circos_{group_name}_{caller}.png")
        plot.plot(output_path, show_legend=True)
        plot_paths.append(output_path)
        print(f"Created: {output_path}")

    return plot_paths


def create_combined_circos_plot(results_dir: str, output_path: str,
                                experiment_file: Optional[str] = None,
                                genome_size: int = 2800000,
                                caller: str = 'bcftools',
                                convergent_only: bool = True) -> str:
    """
    Create a single Circos plot with one ring per treatment.
    Shows only convergent mutations (appearing in multiple replicates).

    Args:
        results_dir: Directory containing sample results
        output_path: Where to save the plot
        experiment_file: CSV/JSON with sample->group mapping
        genome_size: Genome size in bp
        caller: Which variant caller results to use
        convergent_only: Only show mutations appearing in 2+ replicates

    Returns:
        Path to generated plot
    """
    # Load experiment info if provided
    sample_groups = {}
    if experiment_file:
        if experiment_file.endswith('.csv'):
            df = pd.read_csv(experiment_file)
            for _, row in df.iterrows():
                sample_groups[row['sample_id']] = row.get('group', 'Unknown')
        else:
            with open(experiment_file) as f:
                exp = json.load(f)
                for sample_id, sample in exp.get('samples', {}).items():
                    sample_groups[sample_id] = sample.get('group', 'Unknown')

    # Load all sample results and group mutations
    results_path = Path(results_dir)
    mutations_by_group = {}

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

        group = sample_groups.get(sample_name, 'Unknown')
        if group not in mutations_by_group:
            mutations_by_group[group] = {}

        df = pd.read_csv(mut_file)
        for _, row in df.iterrows():
            pos = int(row.get('POS', row.get('position', 0)))
            key = f"{pos}_{row.get('REF', '')}_{row.get('ALT', '')}"

            if key not in mutations_by_group[group]:
                # Extract short gene name
                raw_gene = row.get('GENE', row.get('gene_name', ''))
                product = row.get('PRODUCT', row.get('product', ''))
                locus_tag = row.get('LOCUS_TAG', row.get('locus_tag', ''))
                short_gene = extract_short_gene_name(raw_gene, product, locus_tag)

                mutations_by_group[group][key] = {
                    'position': pos,
                    'effect': row.get('EFFECT', row.get('effect', '')),
                    'gene': short_gene,  # Use short gene name
                    'samples': [],
                    'aa_change': row.get('AA_CHANGE', row.get('amino_acid_change', ''))
                }
            mutations_by_group[group][key]['samples'].append(sample_name)

    # Create plot with one ring per treatment
    plot = CircosPlot(genome_size=genome_size,
                     title=f"Mutations by Treatment ({caller})")

    colors = ['#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#FF7F00', '#A65628']

    for i, (group_name, mutations) in enumerate(sorted(mutations_by_group.items())):
        track_mutations = []

        for key, mut_data in mutations.items():
            # Filter for convergent if requested
            if convergent_only and len(mut_data['samples']) < 2:
                continue

            track_mutations.append(MutationMarker(
                position=mut_data['position'],
                sample=group_name,
                variant_type='SNP',
                effect=mut_data['effect'],
                gene=mut_data['gene'],
                amino_acid_change=mut_data['aa_change']
            ))

        if track_mutations:
            plot.add_mutation_track(
                track_mutations,
                sample_name=f"{group_name} (n={len(track_mutations)})",
                color=colors[i % len(colors)],
                ring_index=i
            )

    plot.plot(output_path, show_legend=True)
    print(f"Created: {output_path}")

    return output_path
