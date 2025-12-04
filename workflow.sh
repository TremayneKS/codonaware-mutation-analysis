#!/usr/bin/env bash
set -euo pipefail

# Simple example workflow:
# 1) BLAST example_gene.fa against example_genomes.fa (isolate_gene.sh)
# 2) Extract best hit per genome + prepend reference gene (extract_gene.py)
# 3) Align all gene sequences (aligner.py)
# 4) Call codon-level mutations (call_mutations.py)

GENE_FASTA="example_gene.fa"
GENOMES_FASTA="example_genomes.fa"
PREFIX="example"

if [[ ! -f "$GENE_FASTA" ]]; then
    echo "Missing $GENE_FASTA" >&2
    exit 1
fi

if [[ ! -f "$GENOMES_FASTA" ]]; then
    echo "Missing $GENOMES_FASTA" >&2
    exit 1
fi

GENE_NAME="$(basename "$GENE_FASTA")"
GENE_NAME="${GENE_NAME%%.*}"

BLAST_TSV="${PREFIX}_blast.tsv"
EXTRACTED_FASTA="${PREFIX}_gene.fa"
ALIGNED_FASTA="aligned_${PREFIX}_gene.fa"

echo "[1/5] BLAST: $GENE_FASTA vs $GENOMES_FASTA"
./isolate_gene.sh "$GENOMES_FASTA" "$GENE_FASTA" "$PREFIX"

if [[ ! -f "$BLAST_TSV" ]]; then
    echo "Expected BLAST output not found: $BLAST_TSV" >&2
    exit 1
fi

echo "[2/5] Extracting gene regions from genomes"
python extract_gene.py "$BLAST_TSV" "$GENOMES_FASTA" "$GENE_FASTA" "$EXTRACTED_FASTA" "$GENE_NAME"

echo "[3/5] Aligning extracted sequences with MAFFT"
python aligner.py "$EXTRACTED_FASTA" "$ALIGNED_FASTA"

echo "[4/5] Calling mutations from aligned CDS"
python call_mutations.py "$ALIGNED_FASTA"
echo "[5/5] Analyzing mutation effects"
python analysis.py "${ALIGNED_FASTA%.fa}"
echo "Done."
echo "Outputs:"
echo "  Aligned FASTA: $ALIGNED_FASTA"
echo "  Mutation CSVs: ${ALIGNED_FASTA%.fa}_nt_aa_events.csv, ${ALIGNED_FASTA%.fa}_nt_position_counts.csv, ${ALIGNED_FASTA%.fa}_effect_counts.csv"
