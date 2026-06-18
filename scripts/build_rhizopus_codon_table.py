#!/usr/bin/env python

import argparse
import json
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from codon_count_functions import CodonFrequencyByAA, fastaFileToStrings


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Build Rhizopus delemar codon table from featureCounts output: "
            "top X%% expressed genes + ribosomal genes."
        )
    )
    p.add_argument(
        "--counts",
        required=True,
        help="featureCounts output file (e.g. Rd_featureCounts.txt)",
    )
    p.add_argument(
        "--cds-fasta",
        default="data/Rd_all_CDS.fa",
        help="FASTA of all CDS (e.g. from gffread). Default: data/Rd_all_CDS.fa",
    )
    p.add_argument(
        "--ribosomal-ids",
        default="data/Rd_ribosomal_gene_ids.txt",
        help="Text file with one ribosomal Geneid per line.",
    )
    p.add_argument(
        "--top-percent",
        type=float,
        default=5.0,
        help="Top X percent expressed genes to use (default: 5.0).",
    )
    p.add_argument(
        "--out-selected-cds",
        default="data/Rhizopus_delemar_top5pct_plus_ribo_CDS.fa",
        help="Output FASTA with CDS of selected genes.",
    )
    p.add_argument(
        "--out-codon-json",
        default="data/Rd_codon_freqs_by_AA.json",
        help="Output JSON with codon frequency table (for inspection/debug).",
    )
    return p.parse_args()


def load_expression_from_featureCounts(counts_path, top_percent):
    # featureCounts lines starting with "#" are comments
    df = pd.read_csv(counts_path, sep="\t", comment="#")
    if "Geneid" not in df.columns:
        raise ValueError(f"'Geneid' column not found in {counts_path}")

    # Identify expression columns:
    # featureCounts has: Geneid, Chr, Start, End, Strand, Length, then 1+ sample columns
    meta_cols = {"Geneid", "Chr", "Start", "End", "Strand", "Length"}
    expr_cols = [c for c in df.columns if c not in meta_cols]

    if not expr_cols:
        raise ValueError("No expression columns found in featureCounts file.")

    # For one sample: expr_cols will have a single column
    # For multiple samples: use mean across samples
    expr = df[expr_cols].mean(axis=1)

    df_expr = pd.DataFrame({
        "Geneid": df["Geneid"].astype(str),
        "expr": expr
    })

    # Sort by descending expression
    df_sorted = df_expr.sort_values("expr", ascending=False)

    n_total = len(df_sorted)
    n_top = max(1, int(round(top_percent / 100.0 * n_total)))
    top_df = df_sorted.head(n_top)

    top_gene_ids = set(top_df["Geneid"])

    print(f"Total genes in counts: {n_total}")
    print(f"Top {top_percent}% = {n_top} genes by expression")
    return top_gene_ids


def load_ribosomal_ids(path):
    ribos = set()
    p = Path(path)
    if not p.exists():
        print(f"WARNING: ribosomal ID file {path} not found; continuing with none.")
        return ribos

    with p.open() as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                ribos.add(line)
    print(f"Loaded {len(ribos)} ribosomal Geneids from {path}")
    return ribos


def gene_id_variants(gid):
    """
    Generate possible variants of an ID to help match between featureCounts
    Geneid and FASTA record IDs.

    Examples:
    - AB3477_000001-T1 -> ["AB3477_000001-T1", "AB3477_000001"]
    - AB3477_000001.1  -> ["AB3477_000001.1", "AB3477_000001"]
    """
    variants = {gid}
    if "-T" in gid:
        variants.add(gid.split("-T")[0])
    if "." in gid:
        variants.add(gid.split(".")[0])
    return variants


def select_cds_fasta(cds_fasta, selected_gene_ids, out_fasta):
    # Build a lookup of all variants of selected Geneids
    selected_variants = set()
    for gid in selected_gene_ids:
        selected_variants |= gene_id_variants(str(gid))

    records_in = list(SeqIO.parse(cds_fasta, "fasta"))
    records_out = []

    for rec in records_in:
        rec_id = rec.id
        rec_variants = gene_id_variants(rec_id)
        if selected_variants & rec_variants:
            records_out.append(rec)

    print(
        f"Selected {len(records_out)} CDS records from {len(records_in)} total in {cds_fasta}"
    )

    if not records_out:
        raise RuntimeError(
            "No CDS records selected. Check that Geneid in counts and FASTA headers match."
        )

    Path(out_fasta).parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records_out, out_fasta, "fasta")
    print(f"Wrote selected CDS FASTA to {out_fasta}")


def build_codon_table_from_fasta(selected_cds_fasta, out_json):
    seq_strings = fastaFileToStrings(selected_cds_fasta)
    codon_freqs = CodonFrequencyByAA(seq_strings, outputtype="freq")

    with open(out_json, "w") as f:
        json.dump(codon_freqs, f, indent=2)

    print(f"Wrote codon frequency table (by amino acid) to {out_json}")
    return codon_freqs


def main():
    args = parse_args()

    # 1) top X% expressed genes from featureCounts
    top_gene_ids = load_expression_from_featureCounts(args.counts, args.top_percent)

    # 2) ribosomal genes
    ribo_ids = load_ribosomal_ids(args.ribosomal_ids)

    # 3) combine
    selected_ids = top_gene_ids | ribo_ids
    print(f"Total selected Geneids (top + ribosomal): {len(selected_ids)}")

    # 4) subset CDS FASTA
    select_cds_fasta(args.cds_fasta, selected_ids, args.out_selected_cds)

    # 5) build codon usage table
    build_codon_table_from_fasta(args.out_selected_cds, args.out_codon_json)


if __name__ == "__main__":
    main()
