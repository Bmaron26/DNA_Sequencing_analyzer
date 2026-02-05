#!/usr/bin/env python3
"""
Compare mutation summary files from different variant callers.
Handles both matrix format (bcftools) and long format (freebayes).

Usage: python compare_summary_files.py bcf_matrix.csv freebayes_long.csv [output_dir]
"""

import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict

def load_matrix_format(filepath):
    """Load mutations from matrix format CSV (bcftools mutation_matrix.csv)."""
    df = pd.read_csv(filepath)

    # Matrix format has: mutation_id, chromosome, position, ref, alt, gene, effect, sample1, sample2, ...
    # Sample columns have allele frequencies (0 = not present)

    info_cols = ['mutation_id', 'chromosome', 'position', 'ref', 'alt', 'gene', 'effect']
    sample_cols = [c for c in df.columns if c not in info_cols]

    variants = {}
    for _, row in df.iterrows():
        pos = row['position']
        ref = row['ref']
        alt = row['alt']
        key = f"{int(pos)}_{ref}_{alt}"

        # Get which samples have this variant (non-zero value)
        samples_with_variant = []
        for sample in sample_cols:
            if row[sample] > 0:
                samples_with_variant.append(sample)

        variants[key] = {
            'position': pos,
            'ref': ref,
            'alt': alt,
            'gene': row['gene'],
            'effect': row['effect'],
            'samples': samples_with_variant,
            'sample_count': len(samples_with_variant)
        }

    return variants, df, sample_cols

def load_long_format(filepath):
    """Load mutations from long format CSV (freebayes all_mutations.csv)."""
    df = pd.read_csv(filepath)

    # Detect column names - check what columns exist
    print(f"  Columns found: {list(df.columns)[:10]}...")

    # Position column
    if 'POS' in df.columns:
        pos_col = 'POS'
    elif 'position' in df.columns:
        pos_col = 'position'
    else:
        pos_col = 'pos'

    # Reference column
    if 'REF' in df.columns:
        ref_col = 'REF'
    elif 'reference' in df.columns:
        ref_col = 'reference'
    else:
        ref_col = 'ref'

    # Alternative column
    if 'ALT' in df.columns:
        alt_col = 'ALT'
    elif 'alternative' in df.columns:
        alt_col = 'alternative'
    else:
        alt_col = 'alt'

    # Sample column
    if 'sample' in df.columns:
        sample_col = 'sample'
    elif 'SAMPLE' in df.columns:
        sample_col = 'SAMPLE'
    else:
        sample_col = None

    # Gene column
    if 'GENE' in df.columns:
        gene_col = 'GENE'
    elif 'gene_name' in df.columns:
        gene_col = 'gene_name'
    elif 'gene' in df.columns:
        gene_col = 'gene'
    else:
        gene_col = None

    # Effect column
    if 'EFFECT' in df.columns:
        effect_col = 'EFFECT'
    elif 'effect' in df.columns:
        effect_col = 'effect'
    else:
        effect_col = None

    print(f"  Using columns: pos={pos_col}, ref={ref_col}, alt={alt_col}, sample={sample_col}")

    variants = {}
    for _, row in df.iterrows():
        pos = row[pos_col]

        # Use direct indexing instead of .get() for pandas Series
        ref = row[ref_col] if ref_col in df.columns else ''
        alt = row[alt_col] if alt_col in df.columns else ''

        # Handle NaN values
        if pd.isna(ref):
            ref = ''
        if pd.isna(alt):
            alt = ''

        if pd.isna(pos):
            continue

        key = f"{int(pos)}_{ref}_{alt}"

        if key not in variants:
            # Get gene and effect with proper null handling
            gene = row[gene_col] if gene_col and gene_col in df.columns else ''
            effect = row[effect_col] if effect_col and effect_col in df.columns else ''
            if pd.isna(gene):
                gene = ''
            if pd.isna(effect):
                effect = ''

            variants[key] = {
                'position': pos,
                'ref': ref,
                'alt': alt,
                'gene': gene,
                'effect': effect,
                'samples': [],
                'sample_count': 0
            }

        if sample_col and sample_col in df.columns:
            sample_val = row[sample_col]
            if not pd.isna(sample_val):
                variants[key]['samples'].append(sample_val)
                variants[key]['sample_count'] = len(variants[key]['samples'])

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
        print("Usage: python compare_summary_files.py bcf_matrix.csv freebayes_long.csv [output_dir]")
        sys.exit(1)

    bcf_file = sys.argv[1]
    fb_file = sys.argv[2]
    output_dir = sys.argv[3] if len(sys.argv) > 3 else "caller_comparison_results"

    os.makedirs(output_dir, exist_ok=True)

    print(f"\n{'='*60}")
    print("VARIANT CALLER COMPARISON: BCFtools vs FreeBayes")
    print(f"{'='*60}\n")
    print(f"BCFtools file: {bcf_file}")
    print(f"FreeBayes file: {fb_file}")
    print()

    # Load BCFtools data (matrix format)
    print("Loading BCFtools data (matrix format)...")
    bcf_variants, bcf_df, sample_cols = load_matrix_format(bcf_file)
    print(f"  Loaded {len(bcf_variants)} unique variants")
    print(f"  Samples: {len(sample_cols)}\n")

    # Load FreeBayes data (long format)
    print("Loading FreeBayes data (long format)...")
    fb_variants, fb_df = load_long_format(fb_file)
    print(f"  Loaded {len(fb_variants)} unique variants\n")

    # Compare
    shared, bcf_only, fb_only = compare_variants(bcf_variants, fb_variants)

    # Print summary
    print("=" * 60)
    print("COMPARISON SUMMARY")
    print("=" * 60)
    print(f"{'Shared variants:':<30} {len(shared):>6}")
    print(f"{'BCFtools only:':<30} {len(bcf_only):>6}")
    print(f"{'FreeBayes only:':<30} {len(fb_only):>6}")
    print(f"{'Total unique:':<30} {len(shared) + len(bcf_only) + len(fb_only):>6}")
    print()

    total = len(shared) + len(bcf_only) + len(fb_only)
    if total > 0:
        concordance = len(shared) / total
        print(f"{'Concordance (Jaccard):':<30} {concordance:.1%}")
    print()

    # Show BCFtools-only variants
    if bcf_only:
        print(f"\n--- BCFtools-only variants ({len(bcf_only)}) ---")
        print("These variants were found by bcftools but NOT by freebayes:")
        for i, key in enumerate(sorted(bcf_only)[:15]):
            v = bcf_variants[key]
            gene = v.get('gene', '-')[:25]
            n_samples = v.get('sample_count', 0)
            print(f"  {v['position']:>10} {v['ref']}->{v['alt']:<5} {gene:<25} (n={n_samples})")
        if len(bcf_only) > 15:
            print(f"  ... and {len(bcf_only) - 15} more")

    # Show FreeBayes-only variants
    if fb_only:
        print(f"\n--- FreeBayes-only variants ({len(fb_only)}) ---")
        print("These variants were found by freebayes but NOT by bcftools:")
        for i, key in enumerate(sorted(fb_only)[:15]):
            v = fb_variants[key]
            gene = str(v.get('gene', '-'))[:25]
            n_samples = v.get('sample_count', 0)
            print(f"  {v['position']:>10} {v['ref']}->{v['alt']:<5} {gene:<25} (n={n_samples})")
        if len(fb_only) > 15:
            print(f"  ... and {len(fb_only) - 15} more")

    # Show shared variants
    if shared:
        print(f"\n--- Shared variants ({len(shared)}) ---")
        print("These variants were found by BOTH callers:")
        for i, key in enumerate(sorted(shared)[:10]):
            v = bcf_variants[key]
            gene = v.get('gene', '-')[:25]
            n_samples = v.get('sample_count', 0)
            print(f"  {v['position']:>10} {v['ref']}->{v['alt']:<5} {gene:<25} (n={n_samples})")
        if len(shared) > 10:
            print(f"  ... and {len(shared) - 10} more")

    # Save detailed comparison
    print(f"\n{'='*60}")
    print(f"Saving detailed results to {output_dir}/...")

    # Shared variants
    if shared:
        shared_data = []
        for key in shared:
            row = {
                'variant_key': key,
                'position': bcf_variants[key]['position'],
                'ref': bcf_variants[key]['ref'],
                'alt': bcf_variants[key]['alt'],
                'gene': bcf_variants[key]['gene'],
                'effect': bcf_variants[key]['effect'],
                'bcf_sample_count': bcf_variants[key]['sample_count'],
                'fb_sample_count': fb_variants[key]['sample_count'],
            }
            shared_data.append(row)
        pd.DataFrame(shared_data).to_csv(f"{output_dir}/shared_variants.csv", index=False)
        print(f"  ✓ shared_variants.csv ({len(shared)} variants)")

    # BCFtools only
    if bcf_only:
        bcf_only_data = [{
            'variant_key': k,
            'position': bcf_variants[k]['position'],
            'ref': bcf_variants[k]['ref'],
            'alt': bcf_variants[k]['alt'],
            'gene': bcf_variants[k]['gene'],
            'effect': bcf_variants[k]['effect'],
            'sample_count': bcf_variants[k]['sample_count'],
            'samples': ';'.join(bcf_variants[k]['samples'])
        } for k in bcf_only]
        pd.DataFrame(bcf_only_data).to_csv(f"{output_dir}/bcftools_only.csv", index=False)
        print(f"  ✓ bcftools_only.csv ({len(bcf_only)} variants)")

    # FreeBayes only
    if fb_only:
        fb_only_data = [{
            'variant_key': k,
            'position': fb_variants[k]['position'],
            'ref': fb_variants[k]['ref'],
            'alt': fb_variants[k]['alt'],
            'gene': fb_variants[k]['gene'],
            'effect': fb_variants[k]['effect'],
            'sample_count': fb_variants[k]['sample_count'],
            'samples': ';'.join(fb_variants[k]['samples'])
        } for k in fb_only]
        pd.DataFrame(fb_only_data).to_csv(f"{output_dir}/freebayes_only.csv", index=False)
        print(f"  ✓ freebayes_only.csv ({len(fb_only)} variants)")

    # Create Venn diagram
    try:
        from matplotlib_venn import venn2

        fig, ax = plt.subplots(figsize=(10, 8))
        v = venn2(subsets=(len(bcf_only), len(fb_only), len(shared)),
              set_labels=('BCFtools', 'FreeBayes'), ax=ax)

        # Style the diagram
        if v.get_patch_by_id('10'):
            v.get_patch_by_id('10').set_color('#377EB8')
            v.get_patch_by_id('10').set_alpha(0.7)
        if v.get_patch_by_id('01'):
            v.get_patch_by_id('01').set_color('#E41A1C')
            v.get_patch_by_id('01').set_alpha(0.7)
        if v.get_patch_by_id('11'):
            v.get_patch_by_id('11').set_color('#4DAF4A')
            v.get_patch_by_id('11').set_alpha(0.7)

        ax.set_title('Variant Caller Comparison\nBCFtools vs FreeBayes',
                    fontsize=14, fontweight='bold')

        # Add summary text
        if total > 0:
            summary_text = f"Concordance: {concordance:.1%}"
            ax.text(0.5, -0.1, summary_text, transform=ax.transAxes,
                    ha='center', fontsize=12, fontweight='bold')

        plt.savefig(f"{output_dir}/venn_comparison.png", dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  ✓ venn_comparison.png")
    except ImportError:
        print("  (matplotlib-venn not available, skipping Venn diagram)")
    except Exception as e:
        print(f"  (Could not create Venn: {e})")

    print(f"\n{'='*60}")
    print("✓ COMPARISON COMPLETE!")
    print(f"{'='*60}")
    print(f"\nResults saved to: {output_dir}/")

    # Interpretation
    print("\n--- INTERPRETATION ---")
    if len(bcf_only) == 0 and len(fb_only) == 0:
        print("Perfect concordance! Both callers found the same variants.")
    elif concordance > 0.9:
        print("Excellent concordance (>90%). Minor differences may be due to")
        print("different filtering parameters.")
    elif concordance > 0.7:
        print("Good concordance (70-90%). FreeBayes may have found additional")
        print("low-frequency variants in your pooled samples.")
    else:
        print("Moderate concordance. Review the unique variants from each caller")
        print("to understand the differences.")

    if len(fb_only) > len(bcf_only):
        print(f"\nFreeBayes found {len(fb_only)} additional variants not in bcftools.")
        print("These may be low-frequency variants in your pooled samples.")

if __name__ == "__main__":
    main()
