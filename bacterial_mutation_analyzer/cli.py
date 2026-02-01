"""
Command-line interface for the bacterial mutation analyzer.
"""

import os
import sys
import logging
from pathlib import Path
from typing import Optional, List

import click
from rich.console import Console
from rich.table import Table
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TaskProgressColumn
from rich.panel import Panel
from rich.logging import RichHandler

from .config import Config, PipelineConfig, create_default_config
from .pipeline.runner import PipelineRunner, run_pipeline
from .experiment import Experiment, SampleMetadata, TreatmentGroup, SampleOrigin

console = Console()


def setup_logging(verbose: bool = False, log_file: Optional[str] = None):
    """Configure logging."""
    level = logging.DEBUG if verbose else logging.INFO

    handlers = [RichHandler(console=console, show_time=True, show_path=False)]

    if log_file:
        handlers.append(logging.FileHandler(log_file))

    logging.basicConfig(
        level=level,
        format="%(message)s",
        handlers=handlers
    )


@click.group()
@click.version_option(version="1.0.0", prog_name="Bacterial Mutation Analyzer")
def main():
    """
    Bacterial Mutation Analyzer - WGS pipeline for mutation identification.

    A comprehensive tool for analyzing bacterial whole genome sequencing data
    to identify mutations, particularly useful for studying antimicrobial
    resistance evolution.
    """
    pass


@main.command()
@click.option('-1', '--fastq1', required=True, type=click.Path(exists=True),
              help='Forward reads (R1) FASTQ file')
@click.option('-2', '--fastq2', type=click.Path(exists=True),
              help='Reverse reads (R2) FASTQ file (optional, for paired-end)')
@click.option('-r', '--reference', required=True, type=click.Path(exists=True),
              help='Reference genome FASTA file')
@click.option('-a', '--annotation', type=click.Path(exists=True),
              help='Annotation file (GFF3 or GenBank format)')
@click.option('-o', '--output', default='results', type=click.Path(),
              help='Output directory (default: results)')
@click.option('-n', '--name', default='sample',
              help='Sample name for output files (default: sample)')
@click.option('-c', '--config', type=click.Path(exists=True),
              help='Configuration file (YAML)')
@click.option('-t', '--threads', default=4, type=int,
              help='Number of threads (default: 4)')
@click.option('--min-quality', default=20, type=int,
              help='Minimum base quality (default: 20)')
@click.option('--min-depth', default=10, type=int,
              help='Minimum read depth for variants (default: 10)')
@click.option('--min-freq', default=0.1, type=float,
              help='Minimum variant frequency (default: 0.1)')
@click.option('--no-visualization', is_flag=True,
              help='Disable visualization generation')
@click.option('-v', '--verbose', is_flag=True,
              help='Enable verbose output')
@click.option('--log-file', type=click.Path(),
              help='Log file path')
def analyze(fastq1: str, fastq2: Optional[str], reference: str,
            annotation: Optional[str], output: str, name: str,
            config: Optional[str], threads: int, min_quality: int,
            min_depth: int, min_freq: float, no_visualization: bool,
            verbose: bool, log_file: Optional[str]):
    """
    Run the complete mutation analysis pipeline.

    This command performs:
    - Quality control and read filtering
    - Alignment to reference genome
    - Variant calling
    - Functional annotation
    - Outlier detection
    - Report generation

    Examples:

    \b
    # Single-end analysis
    bma analyze -1 reads.fastq.gz -r reference.fasta -a annotation.gff

    \b
    # Paired-end analysis with custom output
    bma analyze -1 R1.fastq.gz -2 R2.fastq.gz -r ref.fa -a genes.gff -o my_results -n sample1

    \b
    # Using a configuration file
    bma analyze -1 R1.fq.gz -2 R2.fq.gz -r ref.fa -c config.yaml
    """
    setup_logging(verbose, log_file)

    # Build FASTQ file list
    fastq_files = [fastq1]
    if fastq2:
        fastq_files.append(fastq2)

    # Load or create configuration
    if config:
        pipeline_config = PipelineConfig.from_yaml(config)
    else:
        pipeline_config = PipelineConfig()

    # Override config with CLI options
    pipeline_config.alignment.threads = threads
    pipeline_config.qc.min_quality = min_quality
    pipeline_config.variant_calling.min_depth = min_depth
    pipeline_config.variant_calling.min_variant_frequency = min_freq
    pipeline_config.filtering.min_depth = min_depth
    pipeline_config.visualization.enabled = not no_visualization

    # Print header
    console.print(Panel.fit(
        "[bold blue]Bacterial Mutation Analyzer[/bold blue]\n"
        "[dim]Whole Genome Sequencing Analysis Pipeline[/dim]",
        border_style="blue"
    ))

    # Print input summary
    _print_input_summary(fastq_files, reference, annotation, name)

    # Check tools
    tool_config = Config()
    is_ok, missing = tool_config.check_requirements()
    if not is_ok:
        console.print("[red]Missing required tools:[/red]")
        for tool in missing:
            console.print(f"  [red]✗[/red] {tool}")
        console.print("\nPlease install the missing tools and try again.")
        sys.exit(1)

    # Run pipeline
    try:
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TaskProgressColumn(),
            console=console,
        ) as progress:
            task = progress.add_task("[cyan]Running pipeline...", total=100)

            runner = PipelineRunner(pipeline_config)

            progress.update(task, advance=10, description="[cyan]Validating inputs...")
            progress.update(task, advance=10, description="[cyan]Quality control...")

            result = runner.run(
                fastq_files=fastq_files,
                reference=reference,
                annotation=annotation,
                sample_name=name,
                output_dir=output
            )

            progress.update(task, completed=100, description="[green]Complete!")

        # Print results summary
        _print_results_summary(result)

        # Generate outputs
        if result.status == "completed":
            _export_results(result, output, name, not no_visualization)

            console.print(f"\n[green]✓[/green] Analysis complete!")
            console.print(f"[dim]Results saved to: {output}[/dim]")

            if result.warnings:
                console.print("\n[yellow]Warnings:[/yellow]")
                for warning in result.warnings:
                    console.print(f"  [yellow]![/yellow] {warning}")

            if result.outliers:
                console.print("\n[yellow]Outliers detected:[/yellow]")
                for outlier in result.outliers:
                    console.print(f"  [yellow]![/yellow] {outlier['message']}")

    except Exception as e:
        console.print(f"\n[red]Error:[/red] {str(e)}")
        if verbose:
            console.print_exception()
        sys.exit(1)


@main.command()
@click.option('-o', '--output', default='config.yaml',
              help='Output configuration file path')
def init(output: str):
    """
    Create a default configuration file.

    This generates a YAML configuration file with all available options
    and their default values.
    """
    create_default_config(output)
    console.print(f"[green]✓[/green] Created configuration file: {output}")
    console.print("[dim]Edit this file to customize the pipeline settings.[/dim]")


@main.command()
def check():
    """
    Check availability of required tools.

    Verifies that all required bioinformatics tools are installed
    and accessible in the system PATH.
    """
    config = Config()

    console.print(Panel.fit(
        "[bold]Tool Availability Check[/bold]",
        border_style="blue"
    ))

    # Required tools
    table = Table(title="Required Tools", show_header=True)
    table.add_column("Tool", style="cyan")
    table.add_column("Status", justify="center")
    table.add_column("Description")

    for tool, desc in Config.REQUIRED_TOOLS.items():
        if tool in config.available_tools:
            status = "[green]✓[/green]"
        else:
            status = "[red]✗[/red]"
        table.add_row(tool, status, desc)

    console.print(table)

    # Optional tools
    table2 = Table(title="Optional Tools", show_header=True)
    table2.add_column("Tool", style="cyan")
    table2.add_column("Status", justify="center")
    table2.add_column("Description")

    for tool, desc in Config.OPTIONAL_TOOLS.items():
        if tool in config.available_tools:
            status = "[green]✓[/green]"
        else:
            status = "[dim]-[/dim]"
        table2.add_row(tool, status, desc)

    console.print(table2)

    is_ok, missing = config.check_requirements()
    if is_ok:
        console.print("\n[green]✓[/green] All required tools are available!")
    else:
        console.print(f"\n[red]✗[/red] Missing {len(missing)} required tool(s)")
        sys.exit(1)


@main.command()
@click.argument('result_file', type=click.Path(exists=True))
@click.option('-f', '--format', 'output_format', default='csv',
              type=click.Choice(['csv', 'tsv', 'json', 'excel']),
              help='Output format')
@click.option('-o', '--output', help='Output file path')
def export(result_file: str, output_format: str, output: Optional[str]):
    """
    Export results to various formats.

    Takes a pipeline result JSON file and exports the mutation data
    to CSV, TSV, JSON, or Excel format.
    """
    import json
    import pandas as pd

    with open(result_file, 'r') as f:
        result = json.load(f)

    variants = result.get('variants', [])
    if not variants:
        console.print("[yellow]No variants found in result file.[/yellow]")
        return

    df = pd.DataFrame(variants)

    if output is None:
        base_name = Path(result_file).stem
        output = f"{base_name}_mutations.{output_format}"

    if output_format == 'csv':
        df.to_csv(output, index=False)
    elif output_format == 'tsv':
        df.to_csv(output, index=False, sep='\t')
    elif output_format == 'json':
        df.to_json(output, orient='records', indent=2)
    elif output_format == 'excel':
        df.to_excel(output, index=False)

    console.print(f"[green]✓[/green] Exported {len(variants)} variants to {output}")


@main.command()
@click.argument('result_file', type=click.Path(exists=True))
@click.option('-o', '--output-dir', default='visualizations',
              help='Output directory for plots')
@click.option('--format', 'plot_format', default='png',
              type=click.Choice(['png', 'svg', 'pdf', 'html']),
              help='Plot format')
def visualize(result_file: str, output_dir: str, plot_format: str):
    """
    Generate visualizations from analysis results.

    Creates various plots including:
    - Mutation spectrum
    - Coverage plot
    - Variant distribution
    - Quality metrics
    """
    import json
    from .visualization.plots import VisualizationGenerator

    with open(result_file, 'r') as f:
        result = json.load(f)

    os.makedirs(output_dir, exist_ok=True)

    viz = VisualizationGenerator(output_dir=output_dir, format=plot_format)
    plots = viz.generate_all(result)

    console.print(f"[green]✓[/green] Generated {len(plots)} visualizations in {output_dir}")
    for plot in plots:
        console.print(f"  [dim]- {plot}[/dim]")


@main.command()
@click.option('-e', '--experiment', required=True, type=click.Path(exists=True),
              help='Experiment definition file (JSON or CSV)')
@click.option('-r', '--reference', required=True, type=click.Path(exists=True),
              help='Reference genome FASTA file')
@click.option('-a', '--annotation', type=click.Path(exists=True),
              help='Annotation file (GFF3 or GenBank format)')
@click.option('-o', '--output', default='results', type=click.Path(),
              help='Output directory')
@click.option('--species', default='',
              help='Species name for AMR analysis (e.g., "S. aureus")')
@click.option('-c', '--config', type=click.Path(exists=True),
              help='Pipeline configuration file (YAML)')
@click.option('-t', '--threads', default=4, type=int,
              help='Number of threads per sample')
@click.option('-v', '--verbose', is_flag=True,
              help='Enable verbose output')
def batch(experiment: str, reference: str, annotation: Optional[str],
          output: str, species: str, config: Optional[str],
          threads: int, verbose: bool):
    """
    Run batch analysis for multiple samples with treatment groups.

    Processes all samples defined in an experiment file and generates
    individual and comparison reports.

    Example experiment CSV format:
    sample_id,group,replicate,origin,colony_count,fastq_r1,fastq_r2

    \b
    Examples:
    # Using CSV sample sheet
    bma batch -e samples.csv -r reference.fasta -a genes.gff --species "S. aureus"

    # Using JSON experiment file
    bma batch -e experiment.json -r ref.fa -a genes.gff -o results
    """
    setup_logging(verbose)

    console.print(Panel.fit(
        "[bold blue]Bacterial Mutation Analyzer[/bold blue]\n"
        "[dim]Batch Analysis Mode[/dim]",
        border_style="blue"
    ))

    # Load experiment
    exp_path = Path(experiment)
    if exp_path.suffix == '.csv':
        exp = Experiment.from_csv(experiment, name=exp_path.stem)
    else:
        exp = Experiment.load(experiment)

    # Set reference and annotation
    exp.reference_genome = reference
    exp.annotation_file = annotation or ""
    if species:
        exp.species = species

    console.print(f"[cyan]Experiment:[/cyan] {exp.name}")
    console.print(f"[cyan]Species:[/cyan] {exp.species or 'Not specified'}")
    console.print(f"[cyan]Samples:[/cyan] {len(exp.samples)}")
    console.print(f"[cyan]Groups:[/cyan] {', '.join(exp.groups.keys())}")

    # Print group summary
    table = Table(title="Treatment Groups", show_header=True)
    table.add_column("Group", style="cyan")
    table.add_column("Samples", justify="right")
    table.add_column("Treatment")

    for group_name, group in exp.groups.items():
        sample_count = len(exp.get_samples_by_group(group_name))
        treatment = group.treatment_agent or group.treatment_type or "-"
        table.add_row(group_name, str(sample_count), treatment)

    console.print(table)

    # Load config
    if config:
        pipeline_config = PipelineConfig.from_yaml(config)
    else:
        pipeline_config = PipelineConfig()

    pipeline_config.alignment.threads = threads

    # Process samples
    os.makedirs(output, exist_ok=True)
    results = {}
    failed = []

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TaskProgressColumn(),
        console=console,
    ) as progress:
        task = progress.add_task(
            f"[cyan]Processing samples...",
            total=len(exp.samples)
        )

        for sample_id, sample in exp.samples.items():
            progress.update(task, description=f"[cyan]Processing {sample_id}...")

            # Build FASTQ list
            fastq_files = [sample.fastq_r1]
            if sample.fastq_r2:
                fastq_files.append(sample.fastq_r2)

            # Check files exist
            if not all(os.path.exists(f) for f in fastq_files):
                console.print(f"[yellow]Skipping {sample_id}: FASTQ files not found[/yellow]")
                failed.append(sample_id)
                progress.advance(task)
                continue

            try:
                runner = PipelineRunner(pipeline_config)
                result = runner.run(
                    fastq_files=fastq_files,
                    reference=reference,
                    annotation=annotation,
                    sample_name=sample_id,
                    output_dir=os.path.join(output, sample_id)
                )
                results[sample_id] = result

            except Exception as e:
                console.print(f"[red]Failed {sample_id}: {e}[/red]")
                failed.append(sample_id)

            progress.advance(task)

    console.print(f"\n[green]✓[/green] Processed {len(results)} samples")
    if failed:
        console.print(f"[yellow]![/yellow] Failed: {', '.join(failed)}")

    # Run comparison analysis
    if len(results) >= 2:
        console.print("\n[cyan]Running multi-sample comparison...[/cyan]")

        from .analysis.comparison import MultiSampleComparison

        comparison = MultiSampleComparison(exp)
        for sample_id, result in results.items():
            result_path = os.path.join(output, sample_id, f"{sample_id}_pipeline_result.json")
            if os.path.exists(result_path):
                comparison.load_sample_result(sample_id, result_path)

        # Generate comparison report
        comparison_path = os.path.join(output, "comparison_report.json")
        comparison.generate_comparison_report(comparison_path)

        # Export mutation matrix
        matrix_path = os.path.join(output, "mutation_matrix.csv")
        comparison.export_mutation_matrix(matrix_path)

        # Show convergent mutations
        convergent = comparison.get_convergent_mutations(min_samples=2)
        if convergent:
            console.print(f"\n[bold]Convergent Mutations (n={len(convergent)})[/bold]")
            table = Table(show_header=True)
            table.add_column("Gene")
            table.add_column("Mutation")
            table.add_column("Samples")
            table.add_column("Groups")
            table.add_column("Score", justify="right")

            for cm in convergent[:10]:  # Show top 10
                table.add_row(
                    cm.mutation.gene_name or cm.mutation.locus_tag,
                    cm.mutation.amino_acid_change or f"{cm.mutation.reference}>{cm.mutation.alternative}",
                    str(cm.replicates_affected),
                    ", ".join(cm.groups_affected),
                    f"{cm.convergence_score:.1f}"
                )

            console.print(table)

        console.print(f"\n[dim]Comparison report: {comparison_path}[/dim]")
        console.print(f"[dim]Mutation matrix: {matrix_path}[/dim]")

    # Run AMR annotation if species specified
    if species and results:
        console.print(f"\n[cyan]Running AMR annotation for {species}...[/cyan]")

        from .amr import AMRAnnotator

        annotator = AMRAnnotator(species)

        for sample_id, result in results.items():
            amr_report = annotator.annotate_variants(
                result.variants,
                sample_id=sample_id
            )

            amr_path = os.path.join(output, sample_id, f"{sample_id}_amr_report.json")
            amr_report.save(amr_path)

            if amr_report.known_resistance_mutations > 0:
                console.print(
                    f"  [yellow]{sample_id}[/yellow]: "
                    f"{amr_report.known_resistance_mutations} known resistance mutations"
                )

    console.print(f"\n[green]✓[/green] Batch analysis complete!")
    console.print(f"[dim]Results saved to: {output}[/dim]")


@main.command()
@click.option('-o', '--output', default='experiment.json',
              help='Output experiment file path')
@click.option('--name', default='experiment',
              help='Experiment name')
@click.option('--species', default='',
              help='Species name')
@click.option('--example', is_flag=True,
              help='Create example experiment with sample data')
def create_experiment(output: str, name: str, species: str, example: bool):
    """
    Create a new experiment definition file.

    Creates a JSON or CSV file to define your samples and treatment groups.
    """
    if example:
        from .experiment import create_example_experiment
        exp = create_example_experiment()
    else:
        exp = Experiment(name=name, species=species)

        # Add default control group
        exp.add_group(TreatmentGroup(
            name="Control",
            description="Untreated control",
            treatment_type="control"
        ))

    if output.endswith('.csv'):
        exp.to_csv(output)
    else:
        exp.save(output)

    console.print(f"[green]✓[/green] Created experiment file: {output}")

    if example:
        console.print("[dim]Example experiment created with Melittin and Cecropin groups.[/dim]")
    else:
        console.print("[dim]Edit this file to add your samples and groups.[/dim]")

    console.print("\n[bold]Next steps:[/bold]")
    console.print("1. Edit the experiment file to add your samples")
    console.print("2. Run: bma batch -e experiment.json -r reference.fasta -a annotation.gff")


@main.command()
@click.argument('result_file', type=click.Path(exists=True))
@click.option('--species', required=True,
              help='Species name (e.g., "S. aureus", "E. coli")')
@click.option('-o', '--output', help='Output AMR report path')
def amr(result_file: str, species: str, output: Optional[str]):
    """
    Annotate variants with AMR database information.

    Cross-references mutations with CARD, ResFinder, and species-specific
    resistance databases.

    \b
    Examples:
    bma amr results/sample1/sample1_pipeline_result.json --species "S. aureus"
    """
    import json
    from .amr import AMRAnnotator

    with open(result_file, 'r') as f:
        result = json.load(f)

    sample_id = result.get('sample_name', Path(result_file).stem)

    annotator = AMRAnnotator(species)
    report = annotator.annotate_variants(result.get('variants', []), sample_id)

    if output is None:
        output = str(Path(result_file).with_name(f"{sample_id}_amr_report.json"))

    report.save(output)

    # Print summary
    console.print(Panel.fit(
        f"[bold blue]AMR Analysis: {sample_id}[/bold blue]",
        border_style="blue"
    ))

    table = Table(show_header=False, box=None)
    table.add_column("Metric", style="cyan")
    table.add_column("Value", justify="right")

    table.add_row("Total variants", str(report.total_variants))
    table.add_row("Known resistance mutations", str(report.known_resistance_mutations))
    table.add_row("Mutations in resistance genes", str(report.mutations_in_resistance_genes))
    table.add_row("AMP-related mutations", str(report.amp_related_mutations))

    console.print(table)

    if report.known_hits:
        console.print("\n[bold]Known Resistance Mutations:[/bold]")
        for hit in report.known_hits[:10]:
            info = hit.get('resistance_info', {})
            console.print(
                f"  [yellow]•[/yellow] {hit['gene']} {hit['amino_acid_change']}: "
                f"{info.get('drug_class', 'Unknown')} resistance"
            )

    if report.amp_related_hits:
        console.print("\n[bold]AMP-Related Mutations:[/bold]")
        for hit in report.amp_related_hits[:10]:
            console.print(
                f"  [cyan]•[/cyan] {hit['gene']} {hit.get('amino_acid_change', '')}: "
                f"{hit.get('notes', '')}"
            )

    if report.warnings:
        console.print("\n[yellow]Warnings:[/yellow]")
        for warning in report.warnings:
            console.print(f"  [yellow]![/yellow] {warning}")

    console.print(f"\n[dim]Full report saved to: {output}[/dim]")


@main.command()
@click.argument('results_dir', type=click.Path(exists=True))
@click.option('-e', '--experiment', type=click.Path(exists=True),
              help='Experiment file for group information')
@click.option('-o', '--output', default='comparison',
              help='Output directory for comparison results')
@click.option('--min-samples', default=2, type=int,
              help='Minimum samples for convergent mutation (default: 2)')
def compare(results_dir: str, experiment: Optional[str], output: str,
            min_samples: int):
    """
    Compare mutations across multiple samples.

    Identifies convergent evolution and treatment-specific mutations.

    \b
    Examples:
    bma compare results/ -e experiment.json -o comparison_results
    """
    import json
    from .analysis.comparison import MultiSampleComparison

    # Load experiment if provided
    exp = None
    if experiment:
        exp_path = Path(experiment)
        if exp_path.suffix == '.csv':
            exp = Experiment.from_csv(experiment)
        else:
            exp = Experiment.load(experiment)

    comparison = MultiSampleComparison(exp)

    # Find and load results
    results_path = Path(results_dir)
    loaded = 0

    for sample_dir in results_path.iterdir():
        if not sample_dir.is_dir():
            continue

        result_file = sample_dir / f"{sample_dir.name}_pipeline_result.json"
        if result_file.exists():
            group = None
            if exp and sample_dir.name in exp.samples:
                group = exp.samples[sample_dir.name].group

            comparison.load_sample_result(sample_dir.name, str(result_file), group)
            loaded += 1

    if loaded == 0:
        console.print("[red]No sample results found in directory[/red]")
        return

    console.print(f"[green]✓[/green] Loaded {loaded} samples")

    os.makedirs(output, exist_ok=True)

    # Generate reports
    report_path = os.path.join(output, "comparison_report.json")
    comparison.generate_comparison_report(report_path)

    matrix_path = os.path.join(output, "mutation_matrix.csv")
    comparison.export_mutation_matrix(matrix_path)

    # Show convergent mutations
    convergent = comparison.get_convergent_mutations(min_samples=min_samples)

    console.print(f"\n[bold]Convergent Mutations (appearing in ≥{min_samples} samples)[/bold]")

    if convergent:
        table = Table(show_header=True)
        table.add_column("Gene")
        table.add_column("Change")
        table.add_column("Effect")
        table.add_column("Samples", justify="right")
        table.add_column("Groups")
        table.add_column("Interpretation")

        for cm in convergent[:15]:
            table.add_row(
                cm.mutation.gene_name or cm.mutation.locus_tag or "-",
                cm.mutation.amino_acid_change or
                    f"{cm.mutation.reference}>{cm.mutation.alternative}",
                cm.mutation.effect[:30] if cm.mutation.effect else "-",
                str(cm.replicates_affected),
                ", ".join(cm.groups_affected)[:20],
                cm.interpretation[:40] + "..." if len(cm.interpretation) > 40 else cm.interpretation
            )

        console.print(table)

        # Export convergent mutations
        convergent_path = os.path.join(output, "convergent_mutations.csv")
        import csv
        with open(convergent_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=[
                'gene', 'position', 'reference', 'alternative', 'amino_acid_change',
                'effect', 'sample_count', 'samples', 'groups', 'convergence_score',
                'interpretation'
            ])
            writer.writeheader()
            for cm in convergent:
                writer.writerow({
                    'gene': cm.mutation.gene_name or cm.mutation.locus_tag,
                    'position': cm.mutation.position,
                    'reference': cm.mutation.reference,
                    'alternative': cm.mutation.alternative,
                    'amino_acid_change': cm.mutation.amino_acid_change,
                    'effect': cm.mutation.effect,
                    'sample_count': cm.replicates_affected,
                    'samples': ';'.join(cm.mutation.samples),
                    'groups': ';'.join(cm.groups_affected),
                    'convergence_score': cm.convergence_score,
                    'interpretation': cm.interpretation,
                })

        console.print(f"\n[dim]Convergent mutations: {convergent_path}[/dim]")
    else:
        console.print("[dim]No convergent mutations found[/dim]")

    # Show group summaries
    summaries = comparison.get_all_group_summaries()
    if summaries:
        console.print("\n[bold]Group Summaries[/bold]")
        table = Table(show_header=True)
        table.add_column("Group", style="cyan")
        table.add_column("Samples", justify="right")
        table.add_column("Unique Mutations", justify="right")
        table.add_column("Convergence Rate", justify="right")
        table.add_column("Genes Mutated", justify="right")

        for group_name, summary in summaries.items():
            table.add_row(
                group_name,
                str(summary.get('sample_count', 0)),
                str(summary.get('total_unique_mutations', 0)),
                f"{summary.get('convergence_rate', 0):.1%}",
                str(summary.get('genes_mutated', 0))
            )

        console.print(table)

    console.print(f"\n[dim]Comparison report: {report_path}[/dim]")
    console.print(f"[dim]Mutation matrix: {matrix_path}[/dim]")


@main.command()
def species():
    """
    List available species configurations for AMR analysis.
    """
    from .amr.species_config import list_available_species, get_species_config

    console.print(Panel.fit(
        "[bold]Available Species Configurations[/bold]",
        border_style="blue"
    ))

    species_list = list_available_species()

    for sp in species_list:
        config = get_species_config(sp)
        console.print(f"\n[cyan]{sp}[/cyan]")
        console.print(f"  Resistance genes: {len(config.resistance_genes)}")
        console.print(f"  AMP-related genes: {len(config.amp_resistance_genes)}")
        console.print(f"  Membrane genes: {len(config.membrane_genes)}")

    console.print("\n[dim]Use --species flag with the species name in batch or amr commands[/dim]")


def _print_input_summary(fastq_files: List[str], reference: str,
                        annotation: Optional[str], sample_name: str):
    """Print input file summary."""
    table = Table(title="Input Summary", show_header=False, box=None)
    table.add_column("Parameter", style="cyan")
    table.add_column("Value")

    table.add_row("Sample Name", sample_name)
    table.add_row("Read Type", "Paired-end" if len(fastq_files) == 2 else "Single-end")

    for i, fq in enumerate(fastq_files, 1):
        table.add_row(f"FASTQ R{i}", Path(fq).name)

    table.add_row("Reference", Path(reference).name)

    if annotation:
        table.add_row("Annotation", Path(annotation).name)

    console.print(table)
    console.print()


def _print_results_summary(result):
    """Print analysis results summary."""
    console.print("\n[bold]Results Summary[/bold]")

    # QC Stats
    if result.qc_stats:
        qc = result.qc_stats
        table = Table(title="Quality Control", show_header=False, box=None)
        table.add_column("Metric", style="cyan")
        table.add_column("Value", justify="right")

        table.add_row("Total Reads", f"{qc.get('total_reads', 0):,}")
        table.add_row("Filtered Reads", f"{qc.get('reads_after_filter', 0):,}")
        table.add_row("Q30 Rate", f"{qc.get('q30_rate', 0):.1%}")
        table.add_row("GC Content", f"{qc.get('gc_content', 0):.1%}")

        console.print(table)

    # Alignment Stats
    if result.alignment_stats:
        align = result.alignment_stats
        table = Table(title="Alignment", show_header=False, box=None)
        table.add_column("Metric", style="cyan")
        table.add_column("Value", justify="right")

        table.add_row("Mapping Rate", f"{align.get('mapping_rate', 0):.1%}")
        table.add_row("Mean Coverage", f"{align.get('mean_coverage', 0):.1f}x")
        table.add_row("Coverage Breadth", f"{align.get('coverage_breadth', 0):.1%}")

        console.print(table)

    # Variant Stats
    if result.variant_stats:
        var = result.variant_stats
        table = Table(title="Variants", show_header=False, box=None)
        table.add_column("Metric", style="cyan")
        table.add_column("Value", justify="right")

        table.add_row("Total Variants", f"{var.get('total_variants', 0):,}")
        table.add_row("SNPs", f"{var.get('snps', 0):,}")
        table.add_row("Insertions", f"{var.get('insertions', 0):,}")
        table.add_row("Deletions", f"{var.get('deletions', 0):,}")
        table.add_row("Ti/Tv Ratio", f"{var.get('ti_tv_ratio', 0):.2f}")

        console.print(table)

    # Annotation Summary
    if result.annotation_summary:
        ann = result.annotation_summary
        table = Table(title="Annotation", show_header=False, box=None)
        table.add_column("Metric", style="cyan")
        table.add_column("Value", justify="right")

        table.add_row("Genes Affected", f"{len(ann.get('genes_affected', []))}")
        table.add_row("High Impact", f"{ann.get('by_impact', {}).get('HIGH', 0)}")
        table.add_row("Moderate Impact", f"{ann.get('by_impact', {}).get('MODERATE', 0)}")
        table.add_row("Low Impact", f"{ann.get('by_impact', {}).get('LOW', 0)}")

        console.print(table)


def _export_results(result, output_dir: str, sample_name: str, visualize: bool):
    """Export results to files."""
    import pandas as pd

    os.makedirs(output_dir, exist_ok=True)

    # Export mutations to CSV
    if result.variants:
        mutations_csv = os.path.join(output_dir, f"{sample_name}_mutations.csv")
        mutations_tsv = os.path.join(output_dir, f"{sample_name}_mutations.tsv")

        df = pd.DataFrame(result.variants)
        df.to_csv(mutations_csv, index=False)
        df.to_csv(mutations_tsv, index=False, sep='\t')

        console.print(f"[dim]Mutations exported to: {mutations_csv}[/dim]")

    # Generate visualizations
    if visualize and result.variants:
        try:
            from .visualization.plots import VisualizationGenerator

            viz_dir = os.path.join(output_dir, "visualizations")
            os.makedirs(viz_dir, exist_ok=True)

            viz = VisualizationGenerator(output_dir=viz_dir)
            plots = viz.generate_all(result.to_dict())

            if plots:
                console.print(f"[dim]Visualizations saved to: {viz_dir}[/dim]")
        except ImportError:
            console.print("[dim]Visualization libraries not available[/dim]")
        except Exception as e:
            console.print(f"[dim]Could not generate visualizations: {e}[/dim]")


if __name__ == '__main__':
    main()
