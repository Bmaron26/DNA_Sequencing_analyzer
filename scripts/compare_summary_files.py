#!/usr/bin/env python3
"""
Compare mutation summary files from different variant callers.
Usage: python compare_summary_files.py bcf_file.csv freebayes_file.csv [output_dir]
"""

import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict

def load_mutations(filepath):
    """Load mutations from CSV file and create variant keys."""
    df = pd.read_csv(filepath)

    # Detect column names (handle different naming conventions)
    pos_col = 'POS' if 'POS' in df.columns else 'position' if 'position' in df.columns else None
    ref_col = 'REF' if 'REF' in df.columns else 'reference' if 'reference' in df.columns else None
    alt_col = 'ALT' if 'ALT' in df.columns else 'alternative' if 'alternative' in df.columns else None

    if pos_col is None:
        # Try to find position-like column
        for col in df.columns:
            if 'pos' in col.lower():
                pos_col = col
                break

    print(f"  Columns: {list(df.columns)}")
    print(f"  Using: pos={pos_col}, ref={ref_col}, alt={alt_col}")

    variants = {}

    for _, row in df.iterrows():
        if pos_col and pos_col in df.columns:
            pos = row[pos_col]
            ref = row.get(ref_col, '') if ref_col else ''
            alt = row.get(alt_col, '') if alt_col else ''

            # Handle NaN
            if pd.isna(pos):
                continue

            key = f"{int(pos)}_{ref}_{alt}"
            variants[key] = row.to_dict()

    return variants, df

def compare_variants(bcf_variants, fb_variants):
    """Compare variants between two callers."""
    bcf_keys = set(bcf_variants.keys())
    fb_keys = set(fb_variants.keys())

    shared = bcf_keys & fb_keys
    bcf_only = bcf_keys - fb_keys
    fb_only = fb_keys - bcf_keys

    return shared, bcf_only, fb_only

def main():
    if len(sys.argv) < 3:
        print("Usage: python compare_summary_files.py bcf_file.csv freebayes_file.csv [output_dir]")
        sys.exit(1)

    bcf_file = sys.argv[1]
    fb_file = sys.argv[2]
    output_dir = sys.argv[3] if len(sys.argv) > 3 else "comparison_output"

    os.makedirs(output_dir, exist_ok=True)

    print(f"\n=== Variant Caller Comparison ===\n")
    print(f"BCFtools file: {bcf_file}")
    print(f"FreeBayes file: {fb_file}")
    print()

    # Load data
    print("Loading BCFtools data...")
    bcf_variants, bcf_df = load_mutations(bcf_file)
    print(f"  Loaded {len(bcf_variants)} unique variants\n")

    print("Loading FreeBayes data...")
    fb_variants, fb_df = load_mutations(fb_file)
    print(f"  Loaded {len(fb_variants)} unique variants\n")

    # Compare
    shared, bcf_only, fb_only = compare_variants(bcf_variants, fb_variants)

    # Print summary
    print("=" * 50)
    print("COMPARISON SUMMARY")
    print("=" * 50)
    print(f"Shared variants:      {len(shared):>6}")
    print(f"BCFtools only:        {len(bcf_only):>6}")
    print(f"FreeBayes only:       {len(fb_only):>6}")
    print(f"Total unique:         {len(shared) + len(bcf_only) + len(fb_only):>6}")
    print()

    if len(shared) + len(bcf_only) + len(fb_only) > 0:
        concordance = len(shared) / (len(shared) + len(bcf_only) + len(fb_only))
        print(f"Concordance (Jaccard): {concordance:.1%}")
    print()

    # Show BCFtools-only variants
    if bcf_only:
        print(f"\n--- BCFtools-only variants ({len(bcf_only)}) ---")
        for i, key in enumerate(list(bcf_only)[:10]):
            v = bcf_variants[key]
            gene = v.get('GENE', v.get('gene_name', v.get('gene', '-')))
            print(f"  {key}: {gene}")
        if len(bcf_only) > 10:
            print(f"  ... and {len(bcf_only) - 10} more")

    # Show FreeBayes-only variants
    if fb_only:
        print(f"\n--- FreeBayes-only variants ({len(fb_only)}) ---")
        for i, key in enumerate(list(fb_only)[:10]):
            v = fb_variants[key]
            gene = v.get('GENE', v.get('gene_name', v.get('gene', '-')))
            print(f"  {key}: {gene}")
        if len(fb_only) > 10:
            print(f"  ... and {len(fb_only) - 10} more")

    # Save detailed comparison
    print(f"\nSaving detailed comparison to {output_dir}/...")

    # Shared variants
    if shared:
        shared_data = []
        for key in shared:
            row = {'variant_key': key, 'status': 'shared'}
            row.update({f'bcf_{k}': v for k, v in bcf_variants[key].items()})
            shared_data.append(row)
        pd.DataFrame(shared_data).to_csv(f"{output_dir}/shared_variants.csv", index=False)
        print(f"  Saved: shared_variants.csv ({len(shared)} variants)")

    # BCFtools only
    if bcf_only:
        bcf_only_data = [{'variant_key': k, **bcf_variants[k]} for k in bcf_only]
        pd.DataFrame(bcf_only_data).to_csv(f"{output_dir}/bcftools_only.csv", index=False)
        print(f"  Saved: bcftools_only.csv ({len(bcf_only)} variants)")

    # FreeBayes only
    if fb_only:
        fb_only_data = [{'variant_key': k, **fb_variants[k]} for k in fb_only]
        pd.DataFrame(fb_only_data).to_csv(f"{output_dir}/freebayes_only.csv", index=False)
        print(f"  Saved: freebayes_only.csv ({len(fb_only)} variants)")

    # Create Venn diagram
    try:
        from matplotlib_venn import venn2

        fig, ax = plt.subplots(figsize=(10, 8))
        venn2(subsets=(len(bcf_only), len(fb_only), len(shared)),
              set_labels=('BCFtools', 'FreeBayes'), ax=ax)
        ax.set_title('Variant Caller Comparison', fontsize=14, fontweight='bold')

        # Add summary text
        summary_text = f"Concordance: {concordance:.1%}" if (len(shared) + len(bcf_only) + len(fb_only)) > 0 else ""
        ax.text(0.5, -0.1, summary_text, transform=ax.transAxes,
                ha='center', fontsize=12)

        plt.savefig(f"{output_dir}/venn_comparison.png", dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  Saved: venn_comparison.png")
    except ImportError:
        print("  (matplotlib-venn not available, skipping Venn diagram)")

    print(f"\n✓ Comparison complete!")
    print(f"  Results saved to: {output_dir}/")

if __name__ == "__main__":
    main()
