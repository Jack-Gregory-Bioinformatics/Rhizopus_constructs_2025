#!/usr/bin/env bash
set -euo pipefail

# Activate your environment before running, for example:
# mamba activate hisat2_gffread

# Resolve paths relative to the repository root
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Input/reference files
GENOME_INDEX="${REPO_DIR}/alignment/Rdel_genome/Rdel_hisat2_index"
GFF="${REPO_DIR}/alignment/Rdel_genome/Rdel_AB3477_EC_removed_variants_removed.gff3"

# Output directories
OUTDIR_ALIGN="${REPO_DIR}/alignment/results"
OUTDIR_COUNTS="${REPO_DIR}/alignment/results/counts"

THREADS=8

mkdir -p "$OUTDIR_ALIGN" "$OUTDIR_COUNTS"

# Raw read files
# Update these paths to point to your local or HPC raw data location.
RAW_DATA_DIR="/path/to/raw_data/R_del_RNA-seq"

R1="${RAW_DATA_DIR}/combined_wt_all_1.fastq"
R2="${RAW_DATA_DIR}/combined_wt_all_2.fastq"

BAM="${OUTDIR_ALIGN}/Rdel.sorted.bam"

echo "=== Aligning sample combined_wt_all ==="
hisat2 -x "${GENOME_INDEX}" \
     -1 "${R1}" -2 "${R2}" \
     -p "${THREADS}" \
| samtools view -bS - \
| samtools sort -o "${BAM}" -

samtools index "${BAM}"


# featureCounts over all BAMs together to get one matrix
# relies on having subread in the environment
echo "=== Running featureCounts ==="
featureCounts \
  -T "${THREADS}" \
  -p \
  -t exon \
  -g Parent \
  -a "${GFF}" \
  -o "${OUTDIR_COUNTS}/Rd_featureCounts.txt" \
  ${OUTDIR_ALIGN}/*.sorted.bam

echo "Done. Counts in ${OUTDIR_COUNTS}/Rd_featureCounts.txt"
