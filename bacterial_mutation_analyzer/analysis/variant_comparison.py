"""
Compare variants between different callers (bcftools vs freebayes).
"""

import os
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set
from dataclasses import dataclass, field

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
import seaborn as sns


@dataclass
class VariantKey:
    """Unique identifier for a variant."""
    chrom: str
    pos: int
    ref: str
    alt: str

    def __hash__(self):
        return hash((self.chrom, self.pos, self.ref, self.alt))

    def __eq__(self, other):
        return (self.chrom == other.chrom and self.pos == other.pos and
                self.ref == other.ref and self.alt == other.alt)

    def __str__(self):
        return f"{self.chrom}:{self.pos}:{self.ref}>{self.alt}"


@dataclass
class CallerComparison:
    """Results of comparing two variant callers."""
    caller1: str
    caller2: str
    shared_variants: List[VariantKey] = field(default_factory=list)
    caller1_only: List[VariantKey] = field(default_factory=list)
    caller2_only: List[VariantKey] = field(default_factory=list)
    concordance: float = 0.0
    jaccard_index: float = 0.0

    def summary(self) -> Dict:
        """Return summary statistics."""
        total_unique = len(self.shared_variants) + len(self.caller1_only) + len(self.caller2_only)
        return {
            'caller1': self.caller1,
            'caller2': self.caller2,
            'shared': len(self.shared_variants),
            f'{self.caller1}_only': len(self.caller1_only),
            f'{self.caller2}_only': len(self.caller2_only),
            'total_unique': total_unique,
            'concordance': self.concordance,
            'jaccard_index': self.jaccard_index
        }


class VariantCallerComparison:
    """
    Compare variants called by different tools.
    """

    def __init__(self, results_dir: str):
        """
        Initialize comparison.

        Args:
            results_dir: Directory containing sample results
        """
        self.results_dir = Path(results_dir)
        self.samples = {}
        self._load_samples()

    def _load_samples(self):
        """Load all sample data."""
        for sample_dir in self.results_dir.iterdir():
            if not sample_dir.is_dir():
                continue

            sample_name = sample_dir.name
            self.samples[sample_name] = {
                'dir': sample_dir,
                'callers': {}
            }

            # Look for different caller results
            for caller in ['bcftools', 'freebayes']:
                # Check for caller-specific mutation file
                mut_file = sample_dir / f"{sample_name}_mutations_{caller}.csv"
                if not mut_file.exists() and caller == 'bcftools':
                    # Try default name (bcftools is usually the default)
                    mut_file = sample_dir / f"{sample_name}_mutations.csv"

                if mut_file.exists():
                    self.samples[sample_name]['callers'][caller] = mut_file

    def _load_variants(self, csv_path: Path) -> Dict[VariantKey, Dict]:
        """Load variants from CSV file."""
        variants = {}
        df = pd.read_csv(csv_path)

        for _, row in df.iterrows():
            key = VariantKey(
                chrom=str(row.get('CHROM', row.get('chromosome', ''))),
                pos=int(row.get('POS', row.get('position', 0))),
                ref=str(row.get('REF', row.get('reference', ''))),
                alt=str(row.get('ALT', row.get('alternative', '')))
            )
            variants[key] = row.to_dict()

        return variants

    def compare_sample(self, sample_name: str,
                       caller1: str = 'bcftools',
                       caller2: str = 'freebayes') -> Optional[CallerComparison]:
        """
        Compare variants from two callers for a single sample.

        Args:
            sample_name: Name of sample
            caller1: First caller name
            caller2: Second caller name

        Returns:
            CallerComparison object or None if data missing
        """
        if sample_name not in self.samples:
            return None

        sample = self.samples[sample_name]

        if caller1 not in sample['callers'] or caller2 not in sample['callers']:
            return None

        vars1 = self._load_variants(sample['callers'][caller1])
        vars2 = self._load_variants(sample['callers'][caller2])

        keys1 = set(vars1.keys())
        keys2 = set(vars2.keys())

        shared = keys1 & keys2
        only1 = keys1 - keys2
        only2 = keys2 - keys1

        # Calculate concordance metrics
        total = len(keys1 | keys2)
        concordance = len(shared) / total if total > 0 else 0
        jaccard = len(shared) / total if total > 0 else 0

        return CallerComparison(
            caller1=caller1,
            caller2=caller2,
            shared_variants=list(shared),
            caller1_only=list(only1),
            caller2_only=list(only2),
            concordance=concordance,
            jaccard_index=jaccard
        )

    def compare_all_samples(self, caller1: str = 'bcftools',
                           caller2: str = 'freebayes') -> pd.DataFrame:
        """
        Compare all samples between two callers.

        Returns:
            DataFrame with comparison statistics per sample
        """
        results = []

        for sample_name in self.samples:
            comparison = self.compare_sample(sample_name, caller1, caller2)
            if comparison:
                summary = comparison.summary()
                summary['sample'] = sample_name
                results.append(summary)

        return pd.DataFrame(results)

    def plot_venn_diagram(self, sample_name: str, output_path: str,
                          caller1: str = 'bcftools',
                          caller2: str = 'freebayes') -> str:
        """
        Create Venn diagram for a single sample.

        Args:
            sample_name: Sample to plot
            output_path: Where to save the plot
            caller1: First caller
            caller2: Second caller

        Returns:
            Path to saved plot
        """
        comparison = self.compare_sample(sample_name, caller1, caller2)
        if not comparison:
            raise ValueError(f"Cannot compare {sample_name}: missing data")

        fig, ax = plt.subplots(figsize=(8, 8))

        venn2(
            subsets=(len(comparison.caller1_only),
                    len(comparison.caller2_only),
                    len(comparison.shared_variants)),
            set_labels=(caller1.upper(), caller2.upper()),
            ax=ax
        )

        ax.set_title(f'{sample_name}\nVariant Caller Comparison', fontsize=14)

        # Add stats text
        stats_text = (f"Concordance: {comparison.concordance:.1%}\n"
                     f"Jaccard Index: {comparison.jaccard_index:.2f}")
        ax.text(0.02, 0.02, stats_text, transform=ax.transAxes,
               fontsize=10, verticalalignment='bottom',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

        plt.tight_layout()
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close()

        return output_path

    def plot_concordance_summary(self, output_path: str,
                                 caller1: str = 'bcftools',
                                 caller2: str = 'freebayes') -> str:
        """
        Create summary plot of concordance across all samples.

        Args:
            output_path: Where to save the plot
            caller1: First caller
            caller2: Second caller

        Returns:
            Path to saved plot
        """
        df = self.compare_all_samples(caller1, caller2)
        if df.empty:
            raise ValueError("No samples to compare")

        fig, axes = plt.subplots(2, 2, figsize=(14, 12))

        # 1. Stacked bar chart of variants
        ax1 = axes[0, 0]
        df_sorted = df.sort_values('sample')
        x = range(len(df_sorted))

        ax1.bar(x, df_sorted['shared'], label='Shared', color='#4DAF4A')
        ax1.bar(x, df_sorted[f'{caller1}_only'], bottom=df_sorted['shared'],
               label=f'{caller1} only', color='#377EB8')
        ax1.bar(x, df_sorted[f'{caller2}_only'],
               bottom=df_sorted['shared'] + df_sorted[f'{caller1}_only'],
               label=f'{caller2} only', color='#E41A1C')

        ax1.set_xlabel('Sample')
        ax1.set_ylabel('Number of Variants')
        ax1.set_title('Variant Distribution by Caller')
        ax1.set_xticks(x)
        ax1.set_xticklabels(df_sorted['sample'], rotation=45, ha='right')
        ax1.legend()

        # 2. Concordance bar chart
        ax2 = axes[0, 1]
        colors = ['#4DAF4A' if c > 0.8 else '#FFFF33' if c > 0.5 else '#E41A1C'
                 for c in df_sorted['concordance']]
        ax2.bar(x, df_sorted['concordance'], color=colors)
        ax2.axhline(y=0.8, color='green', linestyle='--', alpha=0.5, label='80% threshold')
        ax2.axhline(y=0.5, color='orange', linestyle='--', alpha=0.5, label='50% threshold')
        ax2.set_xlabel('Sample')
        ax2.set_ylabel('Concordance')
        ax2.set_title('Variant Concordance (Jaccard Index)')
        ax2.set_xticks(x)
        ax2.set_xticklabels(df_sorted['sample'], rotation=45, ha='right')
        ax2.set_ylim(0, 1)
        ax2.legend()

        # 3. Scatter plot: caller1 vs caller2 variant counts
        ax3 = axes[1, 0]
        caller1_total = df['shared'] + df[f'{caller1}_only']
        caller2_total = df['shared'] + df[f'{caller2}_only']
        ax3.scatter(caller1_total, caller2_total, alpha=0.7, s=100)

        # Add diagonal line
        max_val = max(caller1_total.max(), caller2_total.max()) * 1.1
        ax3.plot([0, max_val], [0, max_val], 'k--', alpha=0.3, label='1:1 line')

        ax3.set_xlabel(f'{caller1.upper()} variant count')
        ax3.set_ylabel(f'{caller2.upper()} variant count')
        ax3.set_title('Total Variants: Caller Comparison')
        ax3.legend()

        # Add sample labels
        for i, row in df.iterrows():
            ax3.annotate(row['sample'],
                        (caller1_total.iloc[i], caller2_total.iloc[i]),
                        fontsize=7, alpha=0.7)

        # 4. Summary statistics
        ax4 = axes[1, 1]
        ax4.axis('off')

        summary_text = f"""
        Variant Caller Comparison Summary
        ═══════════════════════════════════

        Samples compared: {len(df)}

        {caller1.upper()} Statistics:
          • Total variants: {(df['shared'] + df[f'{caller1}_only']).sum():,}
          • Unique to {caller1}: {df[f'{caller1}_only'].sum():,}
          • Mean per sample: {(df['shared'] + df[f'{caller1}_only']).mean():.1f}

        {caller2.upper()} Statistics:
          • Total variants: {(df['shared'] + df[f'{caller2}_only']).sum():,}
          • Unique to {caller2}: {df[f'{caller2}_only'].sum():,}
          • Mean per sample: {(df['shared'] + df[f'{caller2}_only']).mean():.1f}

        Overlap Statistics:
          • Shared variants: {df['shared'].sum():,}
          • Mean concordance: {df['concordance'].mean():.1%}
          • Median concordance: {df['concordance'].median():.1%}
        """

        ax4.text(0.1, 0.9, summary_text, transform=ax4.transAxes,
                fontsize=11, verticalalignment='top', fontfamily='monospace',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

        plt.tight_layout()
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close()

        return output_path

    def export_comparison_table(self, output_path: str,
                                caller1: str = 'bcftools',
                                caller2: str = 'freebayes') -> str:
        """
        Export detailed comparison to CSV.

        Args:
            output_path: Where to save the CSV
            caller1: First caller
            caller2: Second caller

        Returns:
            Path to saved CSV
        """
        df = self.compare_all_samples(caller1, caller2)
        df.to_csv(output_path, index=False)
        return output_path

    def export_unique_variants(self, output_dir: str,
                               caller1: str = 'bcftools',
                               caller2: str = 'freebayes') -> Dict[str, str]:
        """
        Export variants unique to each caller for all samples.

        Args:
            output_dir: Directory to save files
            caller1: First caller
            caller2: Second caller

        Returns:
            Dict mapping sample names to output file paths
        """
        os.makedirs(output_dir, exist_ok=True)
        output_files = {}

        for sample_name in self.samples:
            sample = self.samples[sample_name]

            if caller1 not in sample['callers'] or caller2 not in sample['callers']:
                continue

            vars1 = self._load_variants(sample['callers'][caller1])
            vars2 = self._load_variants(sample['callers'][caller2])

            keys1 = set(vars1.keys())
            keys2 = set(vars2.keys())

            # Create comparison DataFrame
            rows = []

            # Shared variants
            for key in keys1 & keys2:
                row = vars1[key].copy()
                row['status'] = 'shared'
                row['caller1_present'] = True
                row['caller2_present'] = True
                rows.append(row)

            # Caller1 only
            for key in keys1 - keys2:
                row = vars1[key].copy()
                row['status'] = f'{caller1}_only'
                row['caller1_present'] = True
                row['caller2_present'] = False
                rows.append(row)

            # Caller2 only
            for key in keys2 - keys1:
                row = vars2[key].copy()
                row['status'] = f'{caller2}_only'
                row['caller1_present'] = False
                row['caller2_present'] = True
                rows.append(row)

            if rows:
                df = pd.DataFrame(rows)
                output_path = os.path.join(output_dir,
                                          f"{sample_name}_caller_comparison.csv")
                df.to_csv(output_path, index=False)
                output_files[sample_name] = output_path

        return output_files


def compare_callers(results_dir: str, output_dir: str,
                    caller1: str = 'bcftools',
                    caller2: str = 'freebayes') -> Dict[str, str]:
    """
    Convenience function to run full caller comparison.

    Args:
        results_dir: Directory with sample results
        output_dir: Where to save comparison outputs
        caller1: First caller
        caller2: Second caller

    Returns:
        Dict of output file paths
    """
    os.makedirs(output_dir, exist_ok=True)

    comparator = VariantCallerComparison(results_dir)

    outputs = {}

    # Summary table
    outputs['summary_csv'] = comparator.export_comparison_table(
        os.path.join(output_dir, f"caller_comparison_{caller1}_vs_{caller2}.csv"),
        caller1, caller2
    )
    print(f"Created: {outputs['summary_csv']}")

    # Concordance plot
    try:
        outputs['concordance_plot'] = comparator.plot_concordance_summary(
            os.path.join(output_dir, f"concordance_{caller1}_vs_{caller2}.png"),
            caller1, caller2
        )
        print(f"Created: {outputs['concordance_plot']}")
    except Exception as e:
        print(f"Could not create concordance plot: {e}")

    # Per-sample Venn diagrams (first 6 samples)
    venn_dir = os.path.join(output_dir, "venn_diagrams")
    os.makedirs(venn_dir, exist_ok=True)

    for i, sample_name in enumerate(list(comparator.samples.keys())[:6]):
        try:
            venn_path = comparator.plot_venn_diagram(
                sample_name,
                os.path.join(venn_dir, f"venn_{sample_name}.png"),
                caller1, caller2
            )
            print(f"Created: {venn_path}")
        except Exception as e:
            print(f"Could not create Venn for {sample_name}: {e}")

    # Export unique variants
    unique_dir = os.path.join(output_dir, "caller_unique_variants")
    outputs['unique_variants'] = comparator.export_unique_variants(
        unique_dir, caller1, caller2
    )
    print(f"Exported unique variants to: {unique_dir}")

    return outputs
