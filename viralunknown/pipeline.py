import subprocess
import sys
from pathlib import Path
from rich.console import Console

console = Console()

WORKFLOW_DIR = Path(__file__).parent.parent / "workflow"


def run_pipeline(config, cores, dryrun, conda_frontend, latency_wait):
    snakefile = WORKFLOW_DIR / "Snakefile"

    cmd = [
        "snakemake",
        "--snakefile",      str(snakefile),
        "--configfile",     config,
        "--cores",          str(cores),
        "--use-conda",
        "--conda-frontend", conda_frontend,
        "--latency-wait",   str(latency_wait),
        "--rerun-incomplete",
    ]

    if dryrun:
        cmd.append("--dry-run")

    console.print(f"[bold]Commande :[/bold] {' '.join(cmd)}")

    result = subprocess.run(cmd)

    if result.returncode != 0:
        console.print("[red]✗The pipeline failed[/red]")
        sys.exit(result.returncode)
    else:
        console.print("[green]✓ Pipeline completed successfully [/green]")