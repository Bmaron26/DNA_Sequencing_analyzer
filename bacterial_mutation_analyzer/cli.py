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
