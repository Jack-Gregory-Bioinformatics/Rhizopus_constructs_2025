#!/usr/bin/env bash
set -euo pipefail

# mamba activate histat2_gffread

# Adjust these to your directory structure
GENOME_INDEX="/lustre/home/jbg209/Research_Project-T117972/Active_Projects/32_FP_R_delemar/Rhizopus_constructs_2025/alignment/Rdel_genome/Rdel_hisat2_index"
GFF="/lustre/home/jbg209/Research_Project-T117972/Active_Projects/32_FP_R_delemar/Rhizopus_constructs_2025/alignment/Rdel_genome/Rdel_AB3477_EC_removed_variants_removed.gff3"
OUTDIR_ALIGN="/lustre/home/jbg209/Research_Project-T117972/Active_Projects/32_FP_R_delemar/Rhizopus_constructs_2025/alignment/results"
OUTDIR_COUNTS="/lustre/home/jbg209/Research_Project-T117972/Active_Projects/32_FP_R_delemar/Rhizopus_constructs_2025/alignment/results/counts"
THREADS=8

mkdir -p "$OUTDIR_ALIGN" "$OUTDIR_COUNTS"


R1="/lustre/home/jbg209/Research_Project-T117972/raw_data/09_LR_Mucorales/R_del_RNA-seq/combined_wt_all_1.fastq"
R2="/lustre/home/jbg209/Research_Project-T117972/raw_data/09_LR_Mucorales/R_del_RNA-seq/combined_wt_all_2.fastq"

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
