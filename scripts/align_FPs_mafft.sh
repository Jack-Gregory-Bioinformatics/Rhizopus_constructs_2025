#!/usr/bin/env bash
set -euo pipefail

OUT_DIR="alignment/FP_AA_alignments"
IN_FASTA="${OUT_DIR}/FPs_all_aa.fasta"
ALN_FASTA="${OUT_DIR}/FPs_all_aa_mafft.aln.fasta"

if [ ! -s "${IN_FASTA}" ]; then
    echo "ERROR: ${IN_FASTA} not found or empty. Run make_FP_multifasta.py first."
    exit 1
fi

echo "Running MAFFT on ${IN_FASTA}..."
# load mafft if needed:
# module load mafft
mafft --auto --thread 8 "${IN_FASTA}" > "${ALN_FASTA}"

echo "Done."
echo "  Alignment: ${ALN_FASTA}"
