#!/bin/bash
set -euo pipefail

mkdir -p references/filter_db
cd references/filter_db

echo "============================================"
echo " Downloading reference sequences"
echo "============================================"

# ── Human rRNA ────────────────────────────────────────────────
echo ""
echo "[ 1/3 ] Human rRNA sequences..."

RRNA_ACC=(
    "NR_003286.4"   # 18S rRNA
    "NR_003287.4"   # 28S rRNA
    "NR_003285.2"   # 5.8S rRNA
    "X12811.1"   # 5S rRNA
    "U13369.1"      # 45S precursor (18S + ITS1 + 5.8S + ITS2 + 28S)
    "NC_012920.1"   # Human mitochondrial genome
)

for acc in "${RRNA_ACC[@]}"; do
    echo "  Downloading $acc..."
    wget -q -O ${acc}.fasta \
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${acc}&rettype=fasta&retmode=text"
    echo "  ✓ ${acc} ($(grep -c '^>' ${acc}.fasta) seq)"
    sleep 0.5
done

cat ${RRNA_ACC[@]/%/.fasta} > rrna_mito.fasta
echo "  → rrna_mito.fasta : $(grep -c '^>' rrna_mito.fasta) sequences"

# ── CMV (HCMV) ───────────────────────────────────────────────
echo ""
echo "[ 2/3 ] CMV (Human cytomegalovirus) genomes..."

CMV_ACC=(
    "NC_006273.2"   # Merlin — reference strain (wild-type)
    "X17403.1"      # AD169  — most used lab strain
    "GQ221974.1"    # Toledo — clinical strain
    "FJ616285.1"    # TB40/E — endothelial tropism strain
    "AC146905.1"    # Towne  — attenuated vaccine strain
    "GU179001.1"    # Merlin BAC (reconstructed wild-type)
)

for acc in "${CMV_ACC[@]}"; do
    echo "  Downloading $acc..."
    wget -q -O ${acc}.fasta \
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${acc}&rettype=fasta&retmode=text"
    echo "  ✓ ${acc} ($(grep -c '^>' ${acc}.fasta) seq)"
    sleep 0.5
done

cat ${CMV_ACC[@]/%/.fasta} > cmv_genomes.fasta
echo "  → cmv_genomes.fasta : $(grep -c '^>' cmv_genomes.fasta) sequences"

# ── Mycoplasma ────────────────────────────────────────────────
echo ""
echo "[ 3/3 ] Mycoplasma genomes..."

MYCO_ACC=(
    "NC_000912.1"   # M. pneumoniae M129 — respiratory pathogen
    "NC_000908.2"   # M. genitalium G37  — urogenital pathogen
    "NC_013511.1"   # M. hominis PG21    — urogenital commensal/pathogen
    "NC_014921.1"   # M. fermentans PG18 — immunocompromised patients
)

for acc in "${MYCO_ACC[@]}"; do
    echo "  Downloading $acc..."
    wget -q -O ${acc}.fasta \
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${acc}&rettype=fasta&retmode=text"
    echo "  ✓ ${acc} ($(grep -c '^>' ${acc}.fasta) seq)"
    sleep 0.5
done

cat ${MYCO_ACC[@]/%/.fasta} > mycoplasma_genomes.fasta
echo "  → mycoplasma_genomes.fasta : $(grep -c '^>' mycoplasma_genomes.fasta) sequences"

# ── Merge all ─────────────────────────────────────────────────
echo ""
echo "Merging all sequences..."
cat rrna_mito.fasta \
    cmv_genomes.fasta \
    mycoplasma_genomes.fasta \
    > filter_db_all.fasta

echo "Total sequences : $(grep -c '^>' filter_db_all.fasta)"

# ── Build BLAST database ──────────────────────────────────────
echo ""
echo "Building BLAST database..."
makeblastdb \
    -in            filter_db_all.fasta \
    -dbtype        nucl \
    -out           filter_db \
    -title         "Human rRNA + mito + CMV + Mycoplasma" \
    -parse_seqids \
    -blastdb_version 5

echo ""
echo "============================================"
echo " Database summary"
echo "============================================"
echo "  rRNA + mito   : $(grep -c '^>' rrna_mito.fasta) sequences"
echo "  CMV genomes   : $(grep -c '^>' cmv_genomes.fasta) sequences"
echo "  Mycoplasma    : $(grep -c '^>' mycoplasma_genomes.fasta) sequences"
echo "  TOTAL         : $(grep -c '^>' filter_db_all.fasta) sequences"
echo ""
echo "  BLAST db      : references/filter_db/filter_db"
echo "  minimap2 ref  : references/filter_db/filter_db_all.fasta"
echo ""
echo "✓ Done"