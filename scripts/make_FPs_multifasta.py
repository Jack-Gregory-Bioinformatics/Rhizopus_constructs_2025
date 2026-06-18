#!/usr/bin/env python

import glob
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def main():
    fp_dir = Path("data")
    out_dir = Path("alignment/FP_AA_alignments")
    out_dir.mkdir(parents=True, exist_ok=True)

    out_fasta = out_dir / "FPs_all_aa.fasta"

    # adjust this glob if your AA files have a different suffix
    fasta_paths = sorted(fp_dir.glob("*.fasta"))

    if not fasta_paths:
        raise SystemExit(f"No .fasta files found in {fp_dir}")

    records_out = []

    for fpath in fasta_paths:
        for rec in SeqIO.parse(str(fpath), "fasta"):
            # Use a clean ID: first token on the header line
            clean_id = rec.id.split()[0]
            # Optional: strip any trailing '*' stop codon from AA sequences
            seq_str = str(rec.seq).replace(" ", "").replace("\n", "").rstrip("*")

            new_rec = SeqRecord(
                rec.seq.__class__(seq_str),
                id=clean_id,
                description=""  # keep headers simple
            )
            records_out.append(new_rec)

    if not records_out:
        raise SystemExit("No sequences parsed from input FASTA files.")

    SeqIO.write(records_out, str(out_fasta), "fasta")
    print(f"Wrote {len(records_out)} FP AA sequences to {out_fasta}")

if __name__ == "__main__":
    main()
