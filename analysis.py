#!/usr/bin/env python3

import argparse
import os
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def find_fasta(base):
    for ext in ["fa", "fasta", "fna"]:
        p = f"{base}.{ext}"
        if os.path.exists(p):
            return p
    sys.exit(f"Missing FASTA for base '{base}'")


def count_fasta_seqs(path):
    c = 0
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                c += 1
    if c == 0:
        sys.exit("FASTA contains no sequences")
    return c


def load_effect_counts(base):
    path = f"{base}_effect_counts.csv"
    if not os.path.exists(path):
        sys.exit("Missing effect-counts CSV")
    df = pd.read_csv(path)
    df["Position"] = df["Position"].astype(int)
    df["Count"] = df["Count"].astype(int)
    df["Event"] = df["Event"].astype(str)
    df["Effect"] = df["Effect"].astype(str)
    df["AA_Mutation"] = df["AA_Mutation"].fillna("").astype(str)
    return df


def build_event_freq(df, n):
    df = df[df["Count"] > 0].copy()
    df["Frequency"] = df["Count"] / n * 100
    return df.sort_values("Position")


def build_pos_totals(df, n):
    df = (
        df.groupby("Position", as_index=False)["Count"]
        .sum()
        .rename(columns={"Count": "Event_Count"})
    )
    df["Frequency"] = df["Event_Count"] / n * 100
    return df.sort_values("Position")


def plot_freq(df_pos, df_events, out_path):
    if df_pos.empty or df_events.empty:
        return

    plt.figure(figsize=(11, 4))
    plt.plot(df_pos["Position"], df_pos["Frequency"], color="black", lw=1)

    syn = df_events["Effect"] == "synonymous"
    mis = df_events["Effect"] == "missense"
    non = df_events["Effect"] == "nonsense"
    stl = df_events["Effect"] == "stop_loss"
    trc = df_events["Effect"] == "truncated"
    fsd = (df_events["Event"] == "DEL") & (df_events["Effect"] == "frameshift_del")
    ifd = (df_events["Event"] == "DEL") & (df_events["Effect"] == "inframe_del")
    fsi = (df_events["Event"] == "INS") & (df_events["Effect"] == "frameshift_ins")
    ifi = (df_events["Event"] == "INS") & (df_events["Effect"] == "inframe_ins")

    known = syn | mis | non | stl | trc | fsd | ifd | fsi | ifi
    assert known.all(), f"Unmapped effects present:\n{df_events[~known]}"

    def scatter(mask, fc, ec, label):
        d = df_events[mask]
        if not d.empty:
            plt.scatter(d["Position"], d["Frequency"],
                        marker="o", s=30,
                        facecolors=fc, edgecolors=ec,
                        lw=0.6, label=label)

    scatter(syn, "none", "black", "synonymous")
    scatter(mis, "black", "black", "missense")
    scatter(non, "purple", "purple", "nonsense")
    scatter(stl, "blue", "blue", "stop-loss")
    scatter(trc, "none", "yellow", "post-trunc")
    scatter(fsd, "red", "red", "frameshift DEL")
    scatter(ifd, "none", "red", "inframe DEL")
    scatter(fsi, "green", "green", "frameshift INS")
    scatter(ifi, "none", "green", "inframe INS")

    plt.xlabel("Nucleotide position")
    plt.ylabel("Mutation frequency (%)")
    plt.legend(fontsize=6, frameon=False)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()


def plot_top(df_eff, out_path, n=20):
    df = df_eff.copy()

    def lab(r):
        if r["AA_Mutation"]:
            return f"{r['Position']}:{r['AA_Mutation']}"
        if r["Effect"]:
            return f"{r['Position']}:{r['Effect']}"
        return f"{r['Position']}:{r['Event']}"

    df["Label"] = df.apply(lab, axis=1)
    df = df.sort_values("Count", ascending=False).head(n)
    plt.figure(figsize=(max(8, n * 0.4), 4))
    plt.bar(df["Label"], df["Count"])
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("base")
    p.add_argument("--top-n", type=int, default=20)
    a = p.parse_args()
    fasta = find_fasta(a.base)
    total = count_fasta_seqs(fasta)
    strains = max(total - 1, 1)
    df = load_effect_counts(a.base)
    df_events = build_event_freq(df, strains)
    df_pos = build_pos_totals(df, strains)
    plot_freq(df_pos, df_events, f"{a.base}_pos_frequency_by_effect.png")
    plot_top(df, f"{a.base}_top_effects.png", n=a.top_n)


if __name__ == "__main__":
    main()
