import click
from rich.console import Console
from rich.panel import Panel
from pathlib import Path
from . import __version__
from .pipeline import run_pipeline
from .config import generate_config, validate_config
from .report import generate_report

console = Console()

@click.group()
@click.version_option(__version__)
def main():
    """ViralUnknown — Pipeline for detecting unknown viruses (Nanopore)"""
    pass


@main.command()
@click.option("--fastq-dir",    required=True,  help="Folder containing the Nanopore FASTQ files")
@click.option("--kraken2-db",   required=True,  help="KrakenDB directory")
@click.option("--output",       default="config.yaml", help="Output config file")
@click.option("--threads",      default=16,     type=int)
@click.option("--memory",       default=30,     type=int, help="Max RAM for SPAdes/Flye (Go)")
@click.option("--min-quality",  default=10,     type=int, help="Min quality score ")
@click.option("--min-len",      default=500,    type=int, help="Min reads length")
@click.option("--blast-db",     required=True,  help="BlastDB directory")
def init(fastq_dir, kraken2_db, output, threads, memory, min_quality, min_len,blast_db):
    """Generates the configuration file from a FASTQ folder."""
    console.print(Panel(f"[bold cyan]ViralUnknown v{__version__}[/bold cyan]\nConfiguration Generation"))
    generate_config(
        fastq_dir   = fastq_dir,
        kraken2_db  = kraken2_db,
        output      = output,
        threads     = threads,
        memory      = memory,
        min_quality = min_quality,
        min_len     = min_len,
        blast_db    = blast_db,
    )
    console.print(f"[green]✓[/green]  : Generated config  {output}")



@main.command()
@click.option("--config",   required=True,  help="YAML config file")
@click.option("--cores",    default=16,     type=int)
@click.option("--dryrun",   is_flag=True,   help="run test")
@click.option("--conda-frontend", default="conda")
@click.option("--latency-wait",   default=30, type=int)
def run(config, cores, dryrun, conda_frontend, latency_wait):
    """Run the Snakemake pipeline."""
    console.print(Panel(f"[bold cyan]ViralUnknown v{__version__}[/bold cyan]\nPipeline Launch"))
    validate_config(config)
    run_pipeline(
        config          = config,
        cores           = cores,
        dryrun          = dryrun,
        conda_frontend  = conda_frontend,
        latency_wait    = latency_wait,
    )


@main.command()
@click.option("--config",   required=True,  help="YAML config file")
@click.option("--output",   default="results/report.html", help="Output HTML file")
def report(config, output):
    """Generates an HTML report of the results."""
    console.print(Panel(f"[bold cyan]ViralUnknown v{__version__}[/bold cyan]\nGenerating the HTML report"))
    generate_report(config=config, output=output)
    console.print(f"[green]✓[/green] Generated report : {output}")

