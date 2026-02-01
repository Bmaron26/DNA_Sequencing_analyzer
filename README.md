# Bacterial Mutation Analyzer

A comprehensive pipeline for whole genome sequencing (WGS) analysis of bacterial samples to identify mutations. Particularly useful for studying antimicrobial resistance evolution and bacterial adaptation under selective pressure (e.g., antimicrobial peptides).

## Features

- **Complete WGS Pipeline**: Quality control, alignment, variant calling, and annotation
- **Flexible Input**: Supports single-end and paired-end FASTQ files (gzipped or plain)
- **Multiple Annotation Formats**: GFF3 and GenBank annotation support
- **Protein Effect Prediction**: Identifies synonymous, missense, nonsense, and frameshift mutations
- **Outlier Detection**: Flags unusual patterns like mutation clusters and quality anomalies
- **Rich Output**: CSV/TSV mutation tables, VCF files, and interactive visualizations
- **Configurable**: YAML-based configuration for all pipeline parameters

## Installation

### Prerequisites

The following bioinformatics tools must be installed and available in your PATH:

**Required:**
- BWA (>= 0.7.17) - Read alignment
- SAMtools (>= 1.10) - BAM file manipulation
- BCFtools (>= 1.10) - Variant calling

**Optional (recommended):**
- fastp - Fast FASTQ preprocessing
- FreeBayes - Alternative variant caller
- bgzip/tabix - VCF compression and indexing

### Install the Python package

```bash
# Clone the repository
git clone https://github.com/your-repo/bacterial-mutation-analyzer.git
cd bacterial-mutation-analyzer

# Create virtual environment (recommended)
python -m venv venv
source venv/bin/activate  # Linux/Mac
# or: venv\Scripts\activate  # Windows

# Install the package
pip install -e .

# Or install dependencies only
pip install -r requirements.txt
```

## Quick Start

### Basic Usage

```bash
# Single-end reads
bma analyze -1 reads.fastq.gz -r reference.fasta -a annotation.gff -n my_sample

# Paired-end reads
bma analyze -1 R1.fastq.gz -2 R2.fastq.gz -r reference.fasta -a annotation.gff -o results

# Check available tools
bma check

# Create default configuration
bma init -o my_config.yaml
```

### Python API

```python
from bacterial_mutation_analyzer import PipelineRunner, PipelineConfig

# Create configuration
config = PipelineConfig()
config.variant_calling.min_depth = 20
config.visualization.enabled = True

# Run pipeline
runner = PipelineRunner(config)
result = runner.run(
    fastq_files=['R1.fastq.gz', 'R2.fastq.gz'],
    reference='reference.fasta',
    annotation='genes.gff',
    sample_name='evolved_strain_1',
    output_dir='results/strain_1'
)

# Access results
print(f"Total variants: {len(result.variants)}")
print(f"High-impact variants: {result.annotation_summary['by_impact'].get('HIGH', 0)}")

# Export mutations
import pandas as pd
df = pd.DataFrame(result.variants)
df.to_csv('mutations.csv', index=False)
```

## Pipeline Steps

### 1. Quality Control
- Adapter trimming
- Quality filtering (Q20/Q30)
- Length filtering
- Complexity filtering
- Generates QC statistics and reports

### 2. Alignment
- BWA-MEM alignment to reference genome
- BAM sorting and indexing
- Duplicate marking
- Generates alignment statistics

### 3. Variant Calling
- BCFtools mpileup + call
- Quality and depth filtering
- VCF output with annotations

### 4. Annotation
- Maps variants to genomic features
- Determines effect type (synonymous, missense, etc.)
- Calculates amino acid changes
- Identifies affected genes

### 5. Outlier Detection
- Quality metric anomalies
- Mutation clustering (potential hypermutation)
- Coverage uniformity issues
- Unusual Ti/Tv ratios

## Output Files

After running the pipeline, you'll find:

```
results/
├── sample_name/
│   ├── qc/
│   │   ├── sample_R1_filtered.fastq.gz
│   │   ├── sample_R2_filtered.fastq.gz
│   │   ├── sample_fastp.json
│   │   └── sample_qc_report.txt
│   ├── alignment/
│   │   ├── sample.sorted.bam
│   │   └── sample.sorted.bam.bai
│   ├── variants/
│   │   ├── sample.filtered.vcf.gz
│   │   └── sample.filtered.vcf.gz.tbi
│   ├── visualizations/
│   │   ├── mutation_spectrum.png
│   │   ├── variant_types.png
│   │   ├── effect_distribution.png
│   │   ├── genome_view.html
│   │   └── summary_dashboard.html
│   ├── sample_mutations.csv
│   ├── sample_mutations.tsv
│   └── sample_pipeline_result.json
```

## Mutation Output Format

The CSV/TSV output includes:

| Column | Description |
|--------|-------------|
| chromosome | Chromosome/contig name |
| position | Genomic position |
| reference | Reference allele |
| alternative | Alternative allele |
| variant_type | SNP, insertion, deletion |
| quality | Variant quality score |
| depth | Read depth |
| allele_frequency | Variant allele frequency |
| gene_name | Gene name |
| locus_tag | Locus tag |
| product | Gene product |
| effect | Variant effect (missense, synonymous, etc.) |
| effect_impact | Impact level (HIGH, MODERATE, LOW) |
| amino_acid_change | Amino acid change (e.g., A123V) |
| codon_change | Codon change (e.g., GCT>GTT) |

## Configuration

Create a configuration file with `bma init` or copy `examples/config.yaml`:

```yaml
# Key configuration options
qc:
  min_quality: 20
  min_length: 50
  adapter_trimming: true

variant_calling:
  min_depth: 10
  min_variant_frequency: 0.1
  ploidy: 1  # Bacteria are haploid

visualization:
  enabled: true
  interactive: true
```

See `examples/config.yaml` for all available options.

## Visualization Examples

The pipeline generates several visualizations:

1. **Mutation Spectrum**: 6-class substitution pattern
2. **Variant Types**: Pie chart of SNPs, insertions, deletions
3. **Effect Distribution**: Bar chart of variant effects and impacts
4. **Genome View**: Interactive variant positions across the genome
5. **Quality Metrics**: Summary of QC and alignment statistics
6. **Summary Dashboard**: Interactive HTML dashboard with all metrics

## Interpreting Results

### Effect Impact Levels

- **HIGH**: Likely disruptive (stop gained, frameshift, start lost)
- **MODERATE**: Possibly functional (missense, in-frame indel)
- **LOW**: Unlikely to affect protein (synonymous)
- **MODIFIER**: Non-coding or intergenic

### Quality Indicators

- **Ti/Tv Ratio**: Expected ~2.0-2.5 for real variants; low values may indicate artifacts
- **Q30 Rate**: Should be >80% for good quality data
- **Mapping Rate**: Should be >90% if using correct reference
- **Coverage**: Recommend ≥30x for reliable variant calling

### Outlier Flags

- **Mutation clusters**: May indicate hypermutation or recombination
- **Low mapping rate**: May indicate wrong reference or contamination
- **Unusual Ti/Tv**: May indicate sequencing artifacts

## Use Cases

### AMP Resistance Evolution

```python
# Compare mutations between ancestral and evolved strains
from bacterial_mutation_analyzer.analysis import MutationStatistics

# Load results for multiple strains
ancestral = load_result('results/ancestral/pipeline_result.json')
evolved_1 = load_result('results/evolved_1/pipeline_result.json')

# Compare mutation patterns
stats1 = MutationStatistics(ancestral['variants'])
stats2 = MutationStatistics(evolved_1['variants'])

comparison = stats1.compare_samples(stats2)
print(f"New mutations in evolved strain: {len(comparison['unique_to_sample2'])}")
```

### Batch Processing

```bash
#!/bin/bash
# Process multiple samples
for sample in sample1 sample2 sample3; do
    bma analyze \
        -1 fastq/${sample}_R1.fastq.gz \
        -2 fastq/${sample}_R2.fastq.gz \
        -r reference.fasta \
        -a annotation.gff \
        -n $sample \
        -o results
done
```

## Troubleshooting

### Common Issues

1. **Low mapping rate**: Check if using correct reference genome
2. **No variants found**: Check coverage depth and quality filtering
3. **Missing annotations**: Ensure GFF/GBK sequence IDs match reference FASTA

### Getting Help

```bash
# Check tool availability
bma check

# Verbose output
bma analyze -1 reads.fq.gz -r ref.fa -v --log-file analysis.log
```

## Citation

If you use this tool in your research, please cite:

```
Bacterial Mutation Analyzer: A WGS pipeline for bacterial mutation identification
[Your publication details]
```

## License

MIT License - see LICENSE file for details.

## Contributing

Contributions are welcome! Please submit issues and pull requests on GitHub.
