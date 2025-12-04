import os
import sys

# Parse BLAST outfmt 6 and keep best hit per subject
def parse_blast_hits(blast_path):
    best_hits = {}
    with open(blast_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 12:
                continue
            qseqid = cols[0]
            sseqid = cols[1]
            try:
                length = int(cols[3])
                sstart = int(cols[8])
                send = int(cols[9])
                evalue = float(cols[10])
                bitscore = float(cols[11])
            except ValueError:
                continue
            hit = {
                "qseqid": qseqid,
                "sseqid": sseqid,
                "length": length,
                "sstart": sstart,
                "send": send,
                "evalue": evalue,
                "bitscore": bitscore,
            }
            current = best_hits.get(sseqid)
            if current is None:
                best_hits[sseqid] = hit
            else:
                if (
                    hit["bitscore"] > current["bitscore"]
                    or (hit["bitscore"] == current["bitscore"] and hit["length"] > current["length"])
                ):
                    best_hits[sseqid] = hit
    return best_hits


# Load multi-FASTA into dict id -> sequence (first token after '>')
def load_fasta_sequences(fasta_path):
    seqs = {}
    current_id = None
    chunks = []
    with open(fasta_path) as f:
        for line in f:
            if not line:
                continue
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_id is not None:
                    seqs[current_id] = "".join(chunks)
                header = line[1:].strip()
                current_id = header.split()[0]
                chunks = []
            else:
                if current_id is not None:
                    chunks.append(line.strip())
    if current_id is not None:
        seqs[current_id] = "".join(chunks)
    return seqs


# Load first FASTA record and return (header_without_gt, sequence)
def load_first_fasta_record(fasta_path):
    header = None
    chunks = []
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                if header is None:
                    header = line[1:].strip()
                    continue
                else:
                    break
            else:
                if header is not None:
                    chunks.append(line.strip())
    if header is None:
        raise SystemExit(f"No sequences found in FASTA: {fasta_path}")
    seq = "".join(chunks)
    if not seq:
        raise SystemExit(f"Empty sequence for first record in FASTA: {fasta_path}")
    return header, seq


# Reverse-complement a nucleotide sequence
def revcomp(seq):
    comp = str.maketrans("ACGTacgtnN", "TGCAtgcanN")
    return seq.translate(comp)[::-1]


# Extract hit regions from genome sequences, handling strand
def extract_regions(seqs, hits):
    extracted = {}
    for sseqid, hit in hits.items():
        if sseqid not in seqs:
            print(f"Warning: {sseqid} in BLAST hits but not in FASTA; skipping", file=sys.stderr)
            continue
        seq = seqs[sseqid]
        sstart = hit["sstart"]
        send = hit["send"]
        if sstart <= send:
            start = sstart - 1
            end = send
            subseq = seq[start:end]
            strand = "+"
        else:
            start = send - 1
            end = sstart
            subseq = seq[start:end]
            subseq = revcomp(subseq)
            strand = "-"
        if not subseq:
            print(f"Warning: empty subsequence for {sseqid}; coords {sstart}-{send}", file=sys.stderr)
            continue
        extracted[sseqid] = {
            "seq": subseq,
            "strand": strand,
            "start": min(sstart, send),
            "end": max(sstart, send),
            "bitscore": hit["bitscore"],
            "length": hit["length"],
            "evalue": hit["evalue"],
            "qseqid": hit["qseqid"],
        }
    return extracted


# Write reference gene first, then extracted hits to FASTA
def write_fasta_with_ref(extracted, out_path, ref_header, ref_seq, gene_name):
    with open(out_path, "w") as out:
        ref_id = ref_header.split()[0]
        out.write(f">{ref_id}|gene={gene_name}|role=reference\n")
        for i in range(0, len(ref_seq), 60):
            out.write(ref_seq[i:i+60] + "\n")
        for sseqid, info in extracted.items():
            header = (
                f">{sseqid}|gene={gene_name}|strand={info['strand']}"
                f"|coords={info['start']}-{info['end']}"
                f"|bitscore={info['bitscore']}"
                f"|length={info['length']}"
                f"|evalue={info['evalue']}"
            )
            out.write(header + "\n")
            seq = info["seq"]
            for i in range(0, len(seq), 60):
                out.write(seq[i:i+60] + "\n")


# CLI entrypoint: blast.tsv genomes.fa gene.fa output.fa [gene_name]
def main():
    if len(sys.argv) < 5 or len(sys.argv) > 6:
        print(
            f"Usage: {sys.argv[0]} blast_hits.tsv genomes.fa gene.fa output.fa [gene_name]",
            file=sys.stderr,
        )
        sys.exit(1)

    blast_file = sys.argv[1]
    genomes_fasta = sys.argv[2]
    gene_fasta = sys.argv[3]
    output_fasta = sys.argv[4]

    if not os.path.exists(blast_file):
        print(f"BLAST file not found: {blast_file}", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(genomes_fasta):
        print(f"Genomes FASTA not found: {genomes_fasta}", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(gene_fasta):
        print(f"Gene FASTA not found: {gene_fasta}", file=sys.stderr)
        sys.exit(1)

    ref_header, ref_seq = load_first_fasta_record(gene_fasta)
    default_gene_name = ref_header.split()[0]
    gene_name = sys.argv[5] if len(sys.argv) == 6 else default_gene_name

    print("Parsing BLAST hits...")
    hits = parse_blast_hits(blast_file)
    if not hits:
        print("No hits found in BLAST file (check outfmt / columns).", file=sys.stderr)
        sys.exit(1)
    print(f"Loaded {len(hits)} best hits (one per subject sequence).")

    print("Loading genome FASTA...")
    seqs = load_fasta_sequences(genomes_fasta)
    print(f"Loaded {len(seqs)} sequences from {genomes_fasta}.")

    print(f"Extracting {gene_name} regions...")
    extracted = extract_regions(seqs, hits)
    print(f"Extracted {len(extracted)} sequences.")

    print(f"Writing reference + extracted sequences to {output_fasta}...")
    write_fasta_with_ref(extracted, output_fasta, ref_header, ref_seq, gene_name=gene_name)
    print("Done.")


if __name__ == "__main__":
    main()
