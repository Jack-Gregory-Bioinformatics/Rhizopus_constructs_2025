#!/usr/bin/env python

import argparse
from pathlib import Path

from Bio import SeqIO


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Standardise N- and/or C-terminal codon blocks for conserved amino-acid motifs\n"
            "across codon-optimised fluorescent proteins.\n\n"
            "You can:\n"
            "  - choose region: N, C, or both\n"
            "  - provide separate master GenBank + AA motif for N and C\n"
            "  - only sequences whose CDS actually has the AA motif at that end are changed."
        )
    )
    p.add_argument(
        "--designs-dir",
        default="designs",
        help="Directory containing FP GenBank designs (default: designs)",
    )
    p.add_argument(
        "--out-dir",
        default="designs_std",
        help="Output directory for standardised GenBank files (default: designs_std)",
    )
    p.add_argument(
        "--mode",
        choices=["N", "C", "both"],
        default="N",
        help="Which termini to standardise: N, C, or both (default: N).",
    )

    # N-terminus options
    p.add_argument(
        "--N-master-gbk",
        help="GenBank file for the master FP (source of N-terminal motif codons).",
    )
    p.add_argument(
        "--N-motif",
        help="N-terminal amino-acid motif to standardise (e.g. MVSKGEE).",
    )

    # C-terminus options
    p.add_argument(
        "--C-master-gbk",
        help="GenBank file for the master FP (source of C-terminal motif codons).",
    )
    p.add_argument(
        "--C-motif",
        help="C-terminal amino-acid motif to standardise (e.g. GGMDELYK).",
    )

    return p.parse_args()


def get_cds_feature(record):
    """Return the first CDS feature in a GenBank record."""
    for feat in record.features:
        if feat.type == "CDS":
            return feat
    raise ValueError(f"No CDS feature found in record {record.id}")


def extract_cds_nt(record, cds_feat):
    """Extract CDS nucleotides from record using CDS feature coordinates."""
    strand = cds_feat.location.strand

    # For FP constructs from DNAchisel, CDS should be on + strand
    if strand not in (1, None):
        # You can extend this if you ever have minus-strand constructs
        raise NotImplementedError(
            f"CDS on non-plus strand in {record.id}; script assumes +1 strand."
        )

    start = int(cds_feat.location.start)
    end = int(cds_feat.location.end)
    return record.seq[start:end], start, end



def extract_N_motif_nt(master_record, motif):
    """Extract nucleotide block encoding the N-terminal motif from master FP."""
    cds_feat = get_cds_feature(master_record)
    cds_nt, _, _ = extract_cds_nt(master_record, cds_feat)
    aa = cds_nt.translate(to_stop=False)

    if not str(aa).startswith(motif):
        raise ValueError(
            f"Master CDS does not start with motif {motif}.\n"
            f"Master translation starts with: {str(aa)[:20]}"
        )

    motif_nt_len = len(motif) * 3
    motif_nt = cds_nt[:motif_nt_len]

    print(
        f"N-terminal: master {master_record.id} -> using motif {motif} codon block ({motif_nt_len} nt)."
    )
    return motif_nt


def extract_C_motif_nt(master_record, motif):
    """Extract nucleotide block encoding the C-terminal motif from master FP."""
    cds_feat = get_cds_feature(master_record)
    cds_nt, _, _ = extract_cds_nt(master_record, cds_feat)
    aa = cds_nt.translate(to_stop=False)

    if not str(aa).endswith(motif):
        raise ValueError(
            f"Master CDS does not end with motif {motif}.\n"
            f"Master translation ends with: {str(aa)[-20:]}"
        )

    motif_nt_len = len(motif) * 3
    motif_nt = cds_nt[-motif_nt_len:]

    print(
        f"C-terminal: master {master_record.id} -> using motif {motif} codon block ({motif_nt_len} nt)."
    )
    return motif_nt


def standardise_record(record, N_motif, N_motif_nt, C_motif, C_motif_nt,
                       do_N, do_C):
    """
    Apply N- and/or C-terminal standardisation to a single record.
    Returns (modified_record, changed_flag).
    """
    changed = False
    cds_feat = get_cds_feature(record)
    cds_nt, cds_start, cds_end = extract_cds_nt(record, cds_feat)

    aa = cds_nt.translate(to_stop=False)
    cds_nt_new = cds_nt

    # N-terminal block
    if do_N and N_motif and N_motif_nt:
        N_nt_len = len(N_motif_nt)
        N_aa_len = N_nt_len // 3

        if str(aa).startswith(N_motif):
            print(f"  {record.id}: N-term motif {N_motif} found -> standardising.")
            if len(cds_nt_new) < N_nt_len:
                raise ValueError(
                    f"CDS in {record.id} shorter than N motif block ({len(cds_nt_new)} < {N_nt_len})."
                )
            cds_nt_new = N_motif_nt + cds_nt_new[N_nt_len:]
            changed = True
        else:
            # Not a target for N-term standardisation
            pass

    # Recompute aa if N changed (for the C-terminal check below)
    if changed and do_N:
        aa = cds_nt_new.translate(to_stop=False)

    # C-terminal block
    if do_C and C_motif and C_motif_nt:
        C_nt_len = len(C_motif_nt)
        C_aa_len = C_nt_len // 3

        if str(aa).endswith(C_motif):
            print(f"  {record.id}: C-term motif {C_motif} found -> standardising.")
            if len(cds_nt_new) < C_nt_len:
                raise ValueError(
                    f"CDS in {record.id} shorter than C motif block ({len(cds_nt_new)} < {C_nt_len})."
                )
            cds_nt_new = cds_nt_new[:-C_nt_len] + C_motif_nt
            changed = True
        else:
            # Not a target for C-term standardisation
            pass

    # If nothing was changed, just return the original
    if not changed:
        return record, False

    # Rebuild full sequence
    full_seq = record.seq
    new_full_seq = full_seq[:cds_start] + cds_nt_new + full_seq[cds_end:]

    record.seq = new_full_seq

    # Remove stale 'translation' qualifier so downstream tools recalc if needed
    for feat in record.features:
        if feat is cds_feat:
            if "translation" in feat.qualifiers:
                del feat.qualifiers["translation"]
            break

    # Final sanity check: translation still matches motifs where expected
    new_cds_nt_check = new_full_seq[cds_start:cds_end]
    new_aa = new_cds_nt_check.translate(to_stop=False)
    if do_N and N_motif and str(new_aa).startswith(N_motif) is False:
        raise RuntimeError(
            f"After N-term replacement, CDS in {record.id} no longer starts with {N_motif}."
        )
    if do_C and C_motif and str(new_aa).endswith(C_motif) is False:
        raise RuntimeError(
            f"After C-term replacement, CDS in {record.id} no longer ends with {C_motif}."
        )

    return record, True


def main():
    args = parse_args()

    designs_dir = Path(args.designs_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    do_N = args.mode in ("N", "both")
    do_C = args.mode in ("C", "both")

    # Validate args
    if do_N and (not args.N_master_gbk or not args.N_motif):
        raise SystemExit(
            "Mode includes N-terminal standardisation, but --N-master-gbk and/or --N-motif not provided."
        )
    if do_C and (not args.C_master_gbk or not args.C_motif):
        raise SystemExit(
            "Mode includes C-terminal standardisation, but --C-master-gbk and/or --C-motif not provided."
        )

    # Load master motifs
    N_motif_nt = None
    C_motif_nt = None

    if do_N:
        N_master_record = SeqIO.read(args.N_master_gbk, "genbank")
        N_motif_nt = extract_N_motif_nt(N_master_record, args.N_motif)

    if do_C:
        C_master_record = SeqIO.read(args.C_master_gbk, "genbank")
        C_motif_nt = extract_C_motif_nt(C_master_record, args.C_motif)

    gbk_paths = sorted(designs_dir.glob("*.gbk"))
    if not gbk_paths:
        raise SystemExit(f"No .gbk files found in {designs_dir}")

    print(f"Scanning {len(gbk_paths)} GenBank files in {designs_dir}...")

    changed_count = 0

    for gbk_path in gbk_paths:
        record = SeqIO.read(str(gbk_path), "genbank")

        new_record, changed = standardise_record(
            record,
            N_motif=args.N_motif if do_N else None,
            N_motif_nt=N_motif_nt,
            C_motif=args.C_motif if do_C else None,
            C_motif_nt=C_motif_nt,
            do_N=do_N,
            do_C=do_C,
        )

        if changed:
            changed_count += 1
            out_name = gbk_path.stem + "_std.gbk"
        else:
            # still output, just without "_std" suffix to keep everything together
            out_name = gbk_path.name

        out_path = out_dir / out_name
        SeqIO.write(new_record, str(out_path), "genbank")

    print("\nStandardisation complete.")
    print(f"  Modified {changed_count} FP constructs.")
    print(f"  Output written to: {out_dir}")


if __name__ == "__main__":
    main()
