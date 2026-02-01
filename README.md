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

### New: Multi-Sample & AMR Features

- **Treatment Groups**: Define experimental groups (e.g., Melittin, Cecropin, Control)
- **Sample Origin Tracking**: Track single colony vs population (~10 colony) samples
- **Batch Processing**: Analyze multiple samples with one command
- **Convergent Evolution Detection**: Identify mutations appearing across parallel replicates
- **AMR Database Integration**: Cross-reference with CARD, ResFinder, and species-specific databases
- **Species-Specific Analysis**: Pre-configured for S. aureus, E. coli, P. aeruginosa
- **AMP Resistance Focus**: Highlight mutations in membrane modification and regulatory genes

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

### Interpreting Population Samples

When samples originate from multiple colonies (~10), low-frequency variants may represent:
- **Subpopulation mutations**: Present in some but not all colonies
- **Emerging resistance**: Early-stage adaptive mutations
- **Hitchhiker mutations**: Neutral mutations in a subset of cells

The tool flags population samples and provides frequency information to help interpretation.

### AMR Interpretation

The AMR annotation provides:
- **Known resistance mutations**: Validated mutations from CARD/ResFinder with clinical evidence
- **Resistance gene mutations**: Novel mutations in known resistance genes (possible emerging resistance)
- **AMP-related mutations**: Mutations in membrane/LPS modification genes relevant for AMP resistance

For S. aureus AMP resistance, key genes include:
- **mprF**: Lysyl-phosphatidylglycerol synthesis (positive membrane charge)
- **dltABCD operon**: D-alanylation of teichoic acids
- **graSR/vraSR**: Two-component regulatory systems
- **Membrane phospholipid genes**: cls, pgsA, pssA

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

### Batch Processing with Treatment Groups

For experiments with multiple samples and treatment groups, use the batch analysis feature:

```bash
# 1. Create a sample sheet (CSV)
# See examples/sample_sheet.csv for format

# 2. Run batch analysis with AMR annotation
bma batch -e samples.csv -r reference.fasta -a genes.gff --species "S. aureus" -o results
```

**Sample sheet format (CSV):**
```csv
sample_id,group,replicate,origin,colony_count,fastq_r1,fastq_r2,notes
Mel1,Melittin,1,population,10,fastq/Mel1_R1.fastq.gz,fastq/Mel1_R2.fastq.gz,Evolved with Melittin
Mel2,Melittin,2,population,10,fastq/Mel2_R1.fastq.gz,fastq/Mel2_R2.fastq.gz,Evolved with Melittin
Ctrl1,Control,1,population,10,fastq/Ctrl1_R1.fastq.gz,fastq/Ctrl1_R2.fastq.gz,No treatment
```

**Sample origin options:**
- `single_colony`: Sample from a single isolated colony
- `population`: Sample from pooled colonies (~10) - may contain subpopulation variants

### Compare Samples Across Groups

```bash
# Compare mutations across samples and identify convergent evolution
bma compare results/ -e samples.csv -o comparison_results

# Output includes:
# - convergent_mutations.csv: Mutations appearing in multiple samples
# - mutation_matrix.csv: Sample x mutation presence matrix
# - comparison_report.json: Full comparison statistics
```

### AMR Database Annotation

```bash
# Annotate variants with known resistance mutations
bma amr results/Mel1/Mel1_pipeline_result.json --species "S. aureus"

# List available species
bma species
```

**Supported species:**
- *Staphylococcus aureus* - includes mprF, dlt operon, vraSR, graSR for AMP resistance
- *Escherichia coli* - includes pmrAB, phoPQ, arn genes for colistin/AMP resistance
- *Pseudomonas aeruginosa* - includes parRS, cprRS systems

### Python API for Multi-Sample Analysis

```python
from bacterial_mutation_analyzer.experiment import Experiment, SampleMetadata, TreatmentGroup, SampleOrigin
from bacterial_mutation_analyzer.analysis import MultiSampleComparison
from bacterial_mutation_analyzer.amr import AMRAnnotator

# Define experiment
exp = Experiment(
    name="AMP_Evolution",
    species="Staphylococcus aureus",
    strain="ATCC 29213"
)

# Add treatment groups
exp.add_group(TreatmentGroup(
    name="Melittin",
    treatment_type="AMP",
    treatment_agent="Melittin",
    concentration="0.5x MIC"
))

# Add samples (6 parallel replicates, each from ~10 colonies)
for i in range(1, 7):
    exp.add_sample(SampleMetadata(
        sample_id=f"Mel{i}",
        group="Melittin",
        replicate=i,
        origin=SampleOrigin.POPULATION,
        colony_count=10,
        fastq_r1=f"fastq/Mel{i}_R1.fastq.gz",
        fastq_r2=f"fastq/Mel{i}_R2.fastq.gz",
    ))

# Run comparison analysis
comparison = MultiSampleComparison(exp)
comparison.load_all_samples("results/")

# Find convergent mutations (parallel evolution)
convergent = comparison.get_convergent_mutations(min_samples=2)
for cm in convergent:
    print(f"{cm.mutation.gene_name}: {cm.mutation.amino_acid_change} "
          f"in {cm.replicates_affected} samples (score: {cm.convergence_score})")

# AMR annotation
annotator = AMRAnnotator("S. aureus")
for sample_id, result in comparison.sample_results.items():
    report = annotator.annotate_variants(result['variants'], sample_id)
    print(f"{sample_id}: {report.amp_related_mutations} AMP-related mutations")
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
