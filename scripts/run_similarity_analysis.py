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

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pathlib import Path
from bacterial_mutation_analyzer.analysis.similarity_analysis import run_similarity_analysis


def main():
    if len(sys.argv) < 3:
        print("Usage: python run_similarity_analysis.py <input_file> <output_dir>")
        print()
        print("Input file can be:")
        print("  - Long format (all_mutations.csv): one row per mutation with 'sample' column")
        print("  - Matrix format (mutation_matrix.csv): samples as columns, variants as rows")
        print()
        print("Example:")
        print("  python scripts/run_similarity_analysis.py \\")
        print("      results_single_amps/combined_analysis/freebayes/all_mutations_freebayes.csv \\")
        print("      similarity_results")
        sys.exit(1)

    input_file = Path(sys.argv[1])
    output_dir = Path(sys.argv[2])

    if not input_file.exists():
        print(f"Error: Input file not found: {input_file}")
        sys.exit(1)

    print("="*60)
    print("MUTATION PROFILE SIMILARITY ANALYSIS")
    print("="*60)
    print(f"Input: {input_file}")
    print(f"Output: {output_dir}")
    print()

    results = run_similarity_analysis(
        input_file,
        output_dir,
        file_format='auto',
        min_edge_weight=0.1
    )

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
