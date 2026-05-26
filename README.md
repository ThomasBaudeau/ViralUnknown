
# ViralUnknown

**ViralUnknown** is a command-line bioinformatics pipeline designed for the detection of unknown or novel viruses in Nanopore long-read sequencing data from tumor samples or any biological material suspected to harbor uncharacterized viral sequences.

## Overview

The pipeline combines multiple filtering, assembly, and clustering strategies to progressively remove known biological signals (human reads, known organisms) and isolate candidate sequences that do not match any known reference — the most likely candidates for novel viral sequences.

## Pipeline steps

Raw Nanopore reads
│
├── Raw statistics (seqkit)
│
▼
rRNA + mitochondrial filtering (minimap2)
│
▼
Low-complexity filtering (dustmasker)
│  Removes reads with >40% masked bases
│  or consecutive runs >20 masked bases
│
▼
Taxonomic classification (Kraken2)
│
├──► Non-human reads (excludes cellular organisms, taxid 131567)
│
└──► Unclassified reads (taxid 0) ──► MetaFlye assembly
│
▼
Unmapped reads extraction
│
Short reads ──► aside
│
▼
Quality filter (Q ≥ 10, seqkit)
│
▼
MMseqs2 clustering (80% identity)
│
▼
Merge assembly + cluster representatives
│
▼
BLASTn against database
│
Filter: keep sequences with
no hit or pident < 80% AND
query coverage < 50%
│
▼
HTML report
## Key features

- **Multi-stage human decontamination** — removes human reads via Kraken2 taxonomic classification, with explicit exclusion of all cellular organisms to maximize sensitivity for viral signal
- **rRNA and mitochondrial filtering** — dedicated minimap2 alignment step against human rRNA (18S, 28S, 5.8S, 5S) and mitochondrial genome (chrM GRCh38) to eliminate common contaminants in RNA sequencing data
- **Low-complexity filtering** — dustmasker-based removal of repetitive and low-complexity sequences (poly-A tails, tandem repeats) that would generate false BLAST hits
- **De novo viral assembly** — MetaFlye with `--meta` mode optimized for low-coverage viral genomes in metagenomic context, with low minimum overlap to capture rare viruses
- **Residual read clustering** — MMseqs2 groups similar unassembled reads into clusters, reducing the BLAST workload from thousands of reads to a handful of representative sequences
- **Custom BLAST database builder** — build targeted databases from any FASTA file (e.g., all known CMV strains, specific virus family) for focused detection
- **Automated HTML report** — interactive Plotly report showing read progression funnels, Kraken2 species pie charts, and per-sample statistics
- **Fully reproducible** — all steps run in isolated conda environments managed by Snakemake

## Installation

```bash
git clone https://github.com/baudeau/viralunknown
cd viralunknown
conda activate snakemake
pip install -e .
```

## Quick start

```bash
# 1. Generate configuration from a folder of FASTQ files
viralunknown init \
    --fastq-dir  /path/to/nanopore/fastq/ \
    --kraken2-db /path/to/kraken2_db \
    --threads    16 \
    --memory     30

# 2. (Optional) Build a custom BLAST database
viralunknown makeblastdb \
    --fasta   references/cmv_strains.fasta \
    --output  databases/cmv/cmv_db \
    --title   "CMV 200 strains"

# 3. Run the pipeline
viralunknown run \
    --config config.yaml \
    --cores  16

# 4. Generate the HTML report
viralunknown report \
    --config config.yaml \
    --output results/report.html
```

## Commands

| Command | Description |
|---|---|
| `viralunknown init` | Generate `config.yaml` from a FASTQ directory |
| `viralunknown run` | Execute the full Snakemake pipeline |
| `viralunknown report` | Generate the interactive HTML report |
| `viralunknown makeblastdb` | Build a BLAST database from a FASTA file |

## Dependencies

All bioinformatics tools are installed automatically in isolated conda environments by Snakemake. The only manual requirement is a working conda/mamba installation and Snakemake ≥ 7.x.

| Tool | Version | Purpose |
|---|---|---|
| Kraken2 | 2.1.3 | Taxonomic classification |
| KrakenTools | 1.2 | Read extraction by taxid |
| Flye / MetaFlye | 2.9.5 | De novo metagenomic assembly |
| minimap2 | 2.28 | Read alignment |
| samtools | 1.21 | BAM processing |
| MMseqs2 | 15.6f452 | Sequence clustering |
| BLAST+ | 2.15.0 | Sequence similarity search |
| seqkit | 2.8.0 | FASTQ/FASTA manipulation |
| FastQC | latest | Quality control |
| MultiQC | 1.21 | QC report aggregation |
| dustmasker | (via BLAST+) | Low-complexity masking |

## Output structure
results/
├── stats/                    # Raw read statistics (seqkit)
├── rrna_mito_filter/         # Reads after rRNA/mito removal
├── dustmasker/               # Reads after low-complexity filtering
├── kraken2/                  # Kraken2 classification reports
├── non_human/                # Non-human classified reads
├── unclassified/             # Reads with no Kraken2 hit
├── metaflye/                 # MetaFlye assembly per sample
├── final/metaflye/           # Unmapped reads after assembly
├── quality_filter/           # High/low quality read separation
├── rrna_mito_filter/         # Final rRNA/mito removal
├── mmseqs/                   # Cluster representatives + summary
├── blast_input/              # Merged assembly + clusters FASTA
├── blast/                    # BLAST results + filtered candidates
├── fastqc/                   # FastQC reports
├── multiqc/                  # MultiQC aggregated report
└── report.html               # Interactive HTML summary report
## Configuration

Key parameters in `config.yaml`:

```yaml
kraken2_db:               "/path/to/kraken2_db"
rrna_mito_ref:            "references/rrna_mito/human_rrna_mito.fasta"
blast_db:                 "/path/to/blast_db"
threads:                  16
flye_memory_gb:           30
min_quality:              10      # Minimum Q score for final reads
min_read_length_assembly: 500     # Minimum read length for assembly
mmseqs_min_seq_id:        0.80    # MMseqs2 clustering identity threshold
blast_min_pident:         80      # Minimum % identity to exclude a hit
blast_min_qcov:           50      # Minimum % query coverage to exclude a hit
```

## License

GNU License 
