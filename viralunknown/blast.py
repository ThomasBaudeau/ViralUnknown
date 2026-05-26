import subprocess
import sys
from pathlib import Path
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn

console = Console()


def build_blast_db(fasta, output, title="custom_db", dbtype="nucl"):
    """
    Build a BLAST nucleotide database from a FASTA file.
    Requires BLAST+ to be installed (makeblastdb).
    """
    fasta_path  = Path(fasta)
    output_path = Path(output)

    # Validation
    if not fasta_path.exists():
        console.print(f"[red]✗ FASTA file not found: {fasta}[/red]")
        sys.exit(1)

    # Compter les séquences
    n_seqs = sum(1 for line in fasta_path.read_text().splitlines()
                 if line.startswith(">"))
    console.print(f"  Sequences found   : [cyan]{n_seqs}[/cyan]")
    console.print(f"  Output database   : [cyan]{output}[/cyan]")
    console.print(f"  Database type     : [cyan]{dbtype}[/cyan]")

    # Créer le dossier de sortie si nécessaire
    output_path.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "makeblastdb",
        "-in",       str(fasta_path),
        "-dbtype",   dbtype,
        "-out",      str(output_path),
        "-title",    title,
        "-parse_seqids",   # permet de récupérer des séquences par ID
        "-blastdb_version", "5",
    ]

    console.print(f"\n[bold]Command:[/bold] {' '.join(cmd)}\n")

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        transient=True,
    ) as progress:
        progress.add_task("Building BLAST database...", total=None)
        result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        console.print(f"[red]✗ makeblastdb failed:[/red]\n{result.stderr}")
        sys.exit(result.returncode)

    # Vérifier les fichiers générés
    extensions = [".nhr", ".nin", ".nsq"] if dbtype == "nucl" else [".phr", ".pin", ".psq"]
    generated = [f"{output}{ext}" for ext in extensions
                 if Path(f"{output}{ext}").exists()]

    console.print(f"[green]✓ Database built successfully[/green]")
    console.print(f"  Files generated:")
    for f in generated:
        size = Path(f).stat().st_size / (1024 * 1024)
        console.print(f"    [dim]{f}[/dim] ({size:.1f} MB)")

    if result.stdout:
        console.print(f"\n[dim]{result.stdout}[/dim]")

    return str(output_path)