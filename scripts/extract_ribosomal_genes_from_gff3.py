#!/usr/bin/env python3
import argparse

def parse_args():
    p = argparse.ArgumentParser(description="Extract ribosomal gene IDs from a GFF3 via GO:0003735")
    p.add_argument("--gff", required=True, help="Input GFF3 annotation")
    p.add_argument("--out", default="data/Rd_ribosomal_gene_ids.txt", help="Output gene list")
    return p.parse_args()

def main():
    args = parse_args()
    ribo_ids = set()

    with open(args.gff) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            
            feature_type = fields[2]
            attrs = fields[8]

            # Only check gene/mRNA/Transcript entries
            if feature_type not in ("gene", "mRNA", "transcript"):
                continue

            # Look for GO:0003735 anywhere in the attributes
            if "GO:0003735" in attrs:
                # Extract the ID=... field
                for part in attrs.split(";"):
                    if part.startswith("ID="):
                        gid = part.replace("ID=", "").strip()
                        ribo_ids.add(gid)
                        break

    print(f"Found {len(ribo_ids)} ribosomal genes with GO:0003735.")

    with open(args.out, "w") as out:
        for gid in sorted(ribo_ids):
            out.write(gid + "\n")

    print(f"Wrote {len(ribo_ids)} IDs to {args.out}")

if __name__ == "__main__":
    main()
