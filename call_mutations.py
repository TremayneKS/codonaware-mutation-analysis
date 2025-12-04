import csv
import sys
from collections import Counter
import os

ALIGNED_FASTA = sys.argv[1]
BASE = os.path.splitext(ALIGNED_FASTA)[0]

DETAIL_CSV = f"{BASE}_nt_aa_events.csv"
SUMMARY_CSV = f"{BASE}_nt_position_counts.csv"
EVENT_EFFECT_SUMMARY_CSV = f"{BASE}_effect_counts.csv"


BASES_STRICT = set("ACGT")

CODON_TABLE = {
    "ATA": "I", "ATC": "I", "ATT": "I", "ATG": "M",
    "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
    "AAC": "N", "AAT": "N", "AAA": "K", "AAG": "K",
    "AGC": "S", "AGT": "S", "AGA": "R", "AGG": "R",
    "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",
    "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
    "CAC": "H", "CAT": "H", "CAA": "Q", "CAG": "Q",
    "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
    "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
    "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
    "GAC": "D", "GAT": "D", "GAA": "E", "GAG": "E",
    "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
    "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
    "TTC": "F", "TTT": "F", "TTA": "L", "TTG": "L",
    "TAC": "Y", "TAT": "Y", "TAA": "*", "TAG": "*",
    "TGC": "C", "TGT": "C", "TGA": "*", "TGG": "W",
}


# Extract a compact strain ID from a FASTA header
def parse_header_id(header_line: str) -> str:
    h = header_line.strip()
    if h.startswith(">"):
        h = h[1:]
    h = h.split()[0]
    if "|" in h:
        h = h.split("|")[0]
    return h


# Read aligned FASTA and return (ids, seqs)
def read_fasta_aligned(path: str):
    ids = []
    seqs = {}
    current_id = None
    current_parts = []
    try:
        with open(path) as f:
            for line in f:
                line = line.rstrip("\n")
                if not line:
                    continue
                if line.startswith(">"):
                    if current_id is not None:
                        seqs[current_id] = "".join(current_parts)
                    current_id = parse_header_id(line)
                    ids.append(current_id)
                    current_parts = []
                else:
                    current_parts.append(line.strip())
            if current_id is not None:
                seqs[current_id] = "".join(current_parts)
    except FileNotFoundError:
        print(f"Error: cannot open FASTA file '{path}'", file=sys.stderr)
        sys.exit(1)

    if not ids:
        print("Error: no sequences found in FASTA.", file=sys.stderr)
        sys.exit(1)

    lengths = {len(s) for s in seqs.values()}
    if len(lengths) != 1:
        print("Warning: sequences have different lengths; using min length.", file=sys.stderr)

    return ids, seqs


# Map alignment index to ungapped reference nt position (1-based)
def build_ref_nt_positions(ref_seq: str):
    L = len(ref_seq)
    ref_pos = [None] * L
    nt_pos = 0
    for i, r in enumerate(ref_seq):
        r = r.upper()
        if r in BASES_STRICT or r == "N":
            nt_pos += 1
            ref_pos[i] = nt_pos
        else:
            ref_pos[i] = None
    return ref_pos, nt_pos


# Find last reference nt position before a given alignment column
def find_anchor_nt(col_index: int, ref_pos):
    j = col_index - 1
    while j >= 0:
        if ref_pos[j] is not None:
            return ref_pos[j]
        j -= 1
    return 0


# Translate ungapped nucleotide CDS with standard bacterial code
def translate_ungapped(seq: str) -> str:
    dna = seq.upper().replace("-", "").replace("\n", "").strip()
    prot = []
    for i in range(0, len(dna), 3):
        codon = dna[i:i+3]
        if len(codon) < 3:
            break
        aa = CODON_TABLE.get(codon, "X")
        prot.append(aa)
        if aa == "*":
            break
    return "".join(prot)


# Classify amino-acid change as synonymous/missense/nonsense/stop_loss
def classify_aa_change(ref_aa: str, alt_aa: str):
    if ref_aa == "" or alt_aa == "":
        return "", ""
    if ref_aa == alt_aa:
        return "synonymous", ""
    if ref_aa != "*" and alt_aa == "*":
        return "nonsense", ""
    if ref_aa == "*" and alt_aa != "*":
        return "stop_loss", ""
    return "missense", ""


# Assign a rough functional effect label from event type and AA effect
def assign_functional_effect(event_type, aa_effect, codon_idx, length_nt, ref_prot_len):
    if aa_effect == "truncated":
        return "irrelevant_after_truncation"

    if aa_effect in ("frameshift_del", "frameshift_ins"):
        return "loss_of_function"

    if aa_effect == "nonsense":
        if ref_prot_len and codon_idx and codon_idx > 0:
            trunc_frac = (ref_prot_len - codon_idx + 1) / ref_prot_len
            if trunc_frac >= 0.10:
                return "loss_of_function"
            if trunc_frac >= 0.05:
                return "likely_damaging"
            return "possibly_tolerated"
        return "loss_of_function"

    if aa_effect in ("inframe_del", "inframe_ins"):
        aa_len = length_nt // 3 if length_nt else 0
        if aa_len >= 3:
            return "likely_damaging"
        if aa_len >= 1:
            return "possibly_tolerated"
        return "unknown"

    if aa_effect == "synonymous":
        return "benign"
    if aa_effect == "missense":
        return "unknown"
    if aa_effect == "stop_loss":
        return "unknown"

    if event_type in ("DEL", "INS") and not aa_effect:
        return "unknown"

    return "unknown"


# Call SNP/DEL/INS events for one strain vs reference and attach AA context
def analyse_strain_nt_and_aa_events(
    strain_id,
    strain_seq,
    ref_seq,
    ref_pos,
    ref_prot,
    strain_prot,
):
    events = []
    L = min(len(ref_seq), len(strain_seq))
    i = 0
    ref_prot_len = len(ref_prot)

    while i < L:
        r = ref_seq[i].upper()
        q = strain_seq[i].upper()

        if r in BASES_STRICT and q == "-":
            start_pos = ref_pos[i]
            if start_pos is None:
                i += 1
                continue
            length = 0
            j = i
            last_pos = start_pos
            while j < L:
                rj = ref_seq[j].upper()
                qj = strain_seq[j].upper()
                if not (rj in BASES_STRICT and qj == "-"):
                    break
                posj = ref_pos[j]
                if posj is not None:
                    last_pos = posj
                length += 1
                j += 1

            codon_idx = (start_pos - 1) // 3 + 1 if start_pos is not None else None

            if length % 3 != 0:
                aa_effect = "frameshift_del"
            else:
                aa_effect = "inframe_del"

            functional_effect = assign_functional_effect(
                event_type="DEL",
                aa_effect=aa_effect,
                codon_idx=codon_idx,
                length_nt=length,
                ref_prot_len=ref_prot_len,
            )

            events.append({
                "strain": strain_id,
                "event": "DEL",
                "start_pos": start_pos,
                "end_pos": last_pos,
                "ref": "",
                "alt": "-",
                "length": length,
                "aa_codon": codon_idx,
                "ref_aa": "",
                "alt_aa": "",
                "aa_mut": "",
                "aa_effect": aa_effect,
                "functional_effect": functional_effect,
            })
            i = j
            continue

        if r == "-" and q in BASES_STRICT:
            anchor_nt = find_anchor_nt(i, ref_pos)
            inserted = []
            length = 0
            j = i
            while j < L:
                rj = ref_seq[j].upper()
                qj = strain_seq[j].upper()
                if not (rj == "-" and qj in BASES_STRICT):
                    break
                inserted.append(qj)
                length += 1
                j += 1
            ins_seq = "".join(inserted)

            codon_idx = (anchor_nt // 3) + 1 if anchor_nt > 0 else 0

            if length % 3 != 0:
                aa_effect = "frameshift_ins"
            else:
                aa_effect = "inframe_ins"

            functional_effect = assign_functional_effect(
                event_type="INS",
                aa_effect=aa_effect,
                codon_idx=codon_idx,
                length_nt=length,
                ref_prot_len=ref_prot_len,
            )

            events.append({
                "strain": strain_id,
                "event": "INS",
                "start_pos": anchor_nt,
                "end_pos": anchor_nt,
                "ref": "-",
                "alt": ins_seq,
                "length": length,
                "aa_codon": codon_idx,
                "ref_aa": "",
                "alt_aa": "",
                "aa_mut": "",
                "aa_effect": aa_effect,
                "functional_effect": functional_effect,
            })
            i = j
            continue

        if r in BASES_STRICT and q in BASES_STRICT and r != q:
            pos = ref_pos[i]
            if pos is not None:
                codon_idx = (pos - 1) // 3 + 1
                ref_aa = ref_prot[codon_idx - 1] if 1 <= codon_idx <= len(ref_prot) else ""

                if codon_idx > len(strain_prot):
                    alt_aa = ""
                    aa_effect = "truncated"
                    aa_mut = f"truncated_at_{codon_idx}"
                else:
                    alt_aa = strain_prot[codon_idx - 1] if codon_idx >= 1 else ""
                    aa_effect, _ = classify_aa_change(ref_aa, alt_aa)
                    aa_mut = f"{ref_aa}{codon_idx}{alt_aa}" if ref_aa and alt_aa else ""

                functional_effect = assign_functional_effect(
                    event_type="SNP",
                    aa_effect=aa_effect,
                    codon_idx=codon_idx,
                    length_nt=1,
                    ref_prot_len=ref_prot_len,
                )

                events.append({
                    "strain": strain_id,
                    "event": "SNP",
                    "start_pos": pos,
                    "end_pos": pos,
                    "ref": r,
                    "alt": q,
                    "length": 1,
                    "aa_codon": codon_idx,
                    "ref_aa": ref_aa,
                    "alt_aa": alt_aa,
                    "aa_mut": aa_mut,
                    "aa_effect": aa_effect,
                    "functional_effect": functional_effect,
                })

        i += 1

    return events


# Write per-strain nucleotide and AA events to CSV
def write_detail_csv(events, path):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "Strain", "Event",
            "Start_Pos", "End_Pos", "Ref", "Alt", "Length",
            "AA_Codon_Index", "Ref_AA", "Alt_AA", "AA_Mutation",
            "AA_Effect", "Functional_Effect",
        ])
        for ev in events:
            w.writerow([
                ev["strain"],
                ev["event"],
                ev["start_pos"],
                ev["end_pos"],
                ev["ref"],
                ev["alt"],
                ev["length"],
                ev["aa_codon"],
                ev["ref_aa"],
                ev["alt_aa"],
                ev["aa_mut"],
                ev["aa_effect"],
                ev["functional_effect"],
            ])


# Build SNP/DEL/INS counts per nucleotide position
def build_position_counts(events):
    counts = {}
    for ev in events:
        ev_type = ev["event"]
        start_pos = ev["start_pos"]
        end_pos = ev["end_pos"]

        if ev_type == "SNP":
            pos = start_pos
            if pos not in counts:
                counts[pos] = {"snp": 0, "del": 0, "ins": 0}
            counts[pos]["snp"] += 1

        elif ev_type == "DEL":
            if start_pos is None or end_pos is None:
                continue
            for pos in range(start_pos, end_pos + 1):
                if pos not in counts:
                    counts[pos] = {"snp": 0, "del": 0, "ins": 0}
                counts[pos]["del"] += 1

        elif ev_type == "INS":
            pos = start_pos
            if pos not in counts:
                counts[pos] = {"snp": 0, "del": 0, "ins": 0}
            counts[pos]["ins"] += 1

    rows = []
    for pos in sorted(counts.keys()):
        c = counts[pos]
        total = c["snp"] + c["del"] + c["ins"]
        rows.append({
            "pos": pos,
            "snp": c["snp"],
            "del": c["del"],
            "ins": c["ins"],
            "total": total,
        })
    return rows


# Write per-position event counts to CSV
def write_summary_csv(rows, path):
    rows_sorted = sorted(rows, key=lambda r: (-r["total"], r["pos"]))
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "NT_Position",
            "SNP_Count", "Deletion_Count", "Insertion_Count", "Total_Count",
        ])
        for r in rows_sorted:
            w.writerow([
                r["pos"],
                r["snp"],
                r["del"],
                r["ins"],
                r["total"],
            ])


# Collapse events into counts by (Position, Event, AA_Effect, AA_Mutation)
def build_event_effect_counts(events):
    counts = Counter()
    for ev in events:
        pos = ev["start_pos"]
        event = ev["event"]
        effect = ev["aa_effect"]
        aa_mut = ev["aa_mut"]
        if pos is None or not event:
            continue
        key = (pos, event, effect, aa_mut)
        counts[key] += 1

    rows = []
    for (pos, event, effect, aa_mut), count in counts.items():
        rows.append({
            "pos": pos,
            "event": event,
            "effect": effect,
            "aa_mut": aa_mut,
            "count": count,
        })
    return rows


# Write event/effect/AA-mutation count table to CSV
def write_event_effect_summary_csv(rows, path):
    rows_sorted = sorted(
        rows,
        key=lambda r: (-r["count"], r["pos"], r["event"], r["effect"], r["aa_mut"]),
    )
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Position", "Event", "Effect", "AA_Mutation", "Count"])
        for r in rows_sorted:
            w.writerow([r["pos"], r["event"], r["effect"], r["aa_mut"], r["count"]])


# Run full pipeline: read FASTA, call events, and write all CSV outputs
def main():
    print(f"Reading aligned FASTA: {ALIGNED_FASTA}")
    ids, seqs = read_fasta_aligned(ALIGNED_FASTA)

    ref_id = ids[0]
    ref_seq = seqs[ref_id]
    print(f"Using reference (first sequence): {ref_id}")

    ref_pos, cds_len = build_ref_nt_positions(ref_seq)
    print(f"Reference CDS length (ungapped): {cds_len} nt")

    ref_prot = translate_ungapped(ref_seq)
    print(f"Reference protein length: {len(ref_prot)} aa")

    all_events = []

    for sid in ids[1:]:
        print(f"Analyzing strain: {sid}")
        strain_seq = seqs[sid]
        strain_prot = translate_ungapped(strain_seq)
        events = analyse_strain_nt_and_aa_events(
            strain_id=sid,
            strain_seq=strain_seq,
            ref_seq=ref_seq,
            ref_pos=ref_pos,
            ref_prot=ref_prot,
            strain_prot=strain_prot,
        )
        all_events.extend(events)

    print(f"Writing detailed events to: {DETAIL_CSV}")
    write_detail_csv(all_events, DETAIL_CSV)

    print("Building per-position counts...")
    rows = build_position_counts(all_events)
    print(f"Writing per-position summary to: {SUMMARY_CSV}")
    write_summary_csv(rows, SUMMARY_CSV)

    print("Building event/effect/AA-mutation counts...")
    ee_rows = build_event_effect_counts(all_events)
    print(f"Writing event/effect summary to: {EVENT_EFFECT_SUMMARY_CSV}")
    write_event_effect_summary_csv(ee_rows, EVENT_EFFECT_SUMMARY_CSV)

    print("Done.")


if __name__ == "__main__":
    main()
