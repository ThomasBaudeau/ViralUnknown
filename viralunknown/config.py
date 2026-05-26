import os
import re
import yaml
from pathlib import Path
from rich.console import Console

console = Console()

def generate_config(fastq_dir, kraken2_db, output, threads, memory,
                    min_quality, min_len,blast_db):
    fastq_dir = Path(fastq_dir)
    fastq_files = sorted([
        f for f in fastq_dir.iterdir()
        if re.search(r"\.fastq(\.gz)?$", f.name)
    ])

    samples = {}
    for f in fastq_files:
        name = re.sub(r"\.fastq(\.gz)?$", "", f.name)
        samples[name] = {"reads": str(f.resolve())}
        console.print(f"  [cyan]→[/cyan] {name}")

    config = {
        "kraken2_db":             str(Path(kraken2_db).resolve()),
        "threads":                threads,
        "flye_memory_gb":         memory,
        "min_quality":            min_quality,
        "min_read_length_assembly": min_len,
        "mmseqs_min_seq_id":      0.80,
        "dust_score_threshold":   2,
        "max_masked_frac":        0.4,
        "samples":                samples,
        "blast_db":               str(Path(blast_db).resolve()),
        "blast_evalue":           "1e-5",
        "blast_max_hits":         5,
        "blast_min_pident":       80,
        "blast_min_qcov":         50,   
    }

    with open(output, "w") as f:
        yaml.dump(config, f, default_flow_style=False, sort_keys=False)

    console.print(f"[green]✓[/green] {len(samples)} detected samples")


def validate_config(config_path):
    with open(config_path) as f:
        config = yaml.safe_load(f)

    required = ["kraken2_db", "threads", "samples"]
    for key in required:
        if key not in config:
            raise ValueError(f"Missing key in the config : {key}")

    if not Path(config["kraken2_db"]).exists():
        raise FileNotFoundError(f"Kraken2 database not found  : {config['kraken2_db']}")

    for sample, paths in config["samples"].items():
        if not Path(paths["reads"]).exists():
            raise FileNotFoundError(f"File not found for {sample} : {paths['reads']}")

    console.print(f"[green]✓[/green] Valid config — {len(config['samples'])} samples")
    return config