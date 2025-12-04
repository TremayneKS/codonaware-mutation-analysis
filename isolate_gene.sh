#!/usr/bin/env bash
set -euo pipefail

# Usage: blast_gene_vs_genomes.sh multi_genomes.fa gene.fa output_prefix

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 multi_genomes.fa gene.fa output_prefix" >&2
    exit 1
fi

MULTI_GENOMES=$1
GENE_FASTA=$2
OUT_PREFIX=$3

if [[ ! -f "$MULTI_GENOMES" ]]; then
    echo "Error: multi-genome FASTA not found: $MULTI_GENOMES" >&2
    exit 1
fi

if [[ ! -f "$GENE_FASTA" ]]; then
    echo "Error: gene FASTA not found: $GENE_FASTA" >&2
    exit 1
fi

if ! command -v makeblastdb >/dev/null 2>&1; then
    echo "Error: makeblastdb not found on PATH (install BLAST+)." >&2
    exit 1
fi

if ! command -v blastn >/dev/null 2>&1; then
    echo "Error: blastn not found on PATH (install BLAST+)." >&2
    exit 1
fi

DB_PREFIX="${OUT_PREFIX}_db"
OUT_TSV="${OUT_PREFIX}_blast.tsv"

echo "Building BLAST database from: $MULTI_GENOMES"
makeblastdb -in "$MULTI_GENOMES" -dbtype nucl -out "$DB_PREFIX" >/dev/null

echo "Running BLAST: $GENE_FASTA vs $MULTI_GENOMES"
blastn \
    -query "$GENE_FASTA" \
    -db "$DB_PREFIX" \
    -out "$OUT_TSV" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

echo "BLAST results written to: $OUT_TSV"

echo "Cleaning up BLAST database files..."
rm -f "${DB_PREFIX}".*
