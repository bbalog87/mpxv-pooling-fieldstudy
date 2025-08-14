#!/bin/bash

###############################################################################
# Pipeline: Mpox Variant & Coverage Summary
# Author: Julien Nguinkal
# Description:
#   Processes paired-end FASTQ files to:
#     - Align to MPXV reference
#     - Mark duplicates
#     - Compute coverage
#     - Call variants with iVar
#     - Generate consensus sequences
#     - Extract variant stats
###############################################################################

# ===================== PARSE ARGUMENTS =====================
usage() {
    echo -e "\nUsage: $0 -i <input_fastq_folder> -o <output_project_folder> -r <reference_fasta> -t <threads>\n"
    exit 1
}

while getopts ":i:o:r:t:h" opt; do
    case $opt in
        i) IN_DIR="$OPTARG" ;;
        o) OUT_DIR="$OPTARG" ;;
        r) REF="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        h|*) usage ;;
    esac
done

[[ -z "$IN_DIR" || -z "$OUT_DIR" || -z "$REF" || -z "$THREADS" ]] && usage

# ===================== SETUP =====================
mkdir -p "$OUT_DIR"
SUMMARY_FILE="$OUT_DIR/summary.tsv"
echo -e "Sample\tMeanDepth\tPctCovered10x\tTotalVariants\tVariantsQUAL20\tVariantsQUAL30\tPctQUAL20\tPctQUAL30\tTsTv" > "$SUMMARY_FILE"

bwa index "$REF"
samtools faidx "$REF"

# ===================== PER SAMPLE LOOP =====================
for R1 in "$IN_DIR"/*_R1.fastq.gz; do
    sample=$(basename "$R1" _R1.fastq.gz)
    R2="$IN_DIR/${sample}_R2.fastq.gz"
    SAMPLE_DIR="$OUT_DIR/sample_${sample}"
    mkdir -p "$SAMPLE_DIR"

    echo -e "\n[INFO] Processing $sample..."

    # === ALIGNMENT ===
    bwa mem -t "$THREADS" -R "@RG\tID:$sample\tSM:$sample\tPL:ILLUMINA" "$REF" "$R1" "$R2" |
        samtools view -b - | samtools sort -@ "$THREADS" -o "$SAMPLE_DIR/aligned.bam"

    # === DEDUPLICATION ===
    samtools view -b -f 1 -F 12 "$SAMPLE_DIR/aligned.bam" > "$SAMPLE_DIR/paired.bam"
    picard MarkDuplicates -I "$SAMPLE_DIR/paired.bam" \
                          -O "$SAMPLE_DIR/dedup.bam" \
                          -M "$SAMPLE_DIR/dup_metrics.txt" \
                          --REMOVE_DUPLICATES true \
                          --ASSUME_SORT_ORDER coordinate

    # === COVERAGE ===
    samtools depth -aa "$SAMPLE_DIR/dedup.bam" > "$SAMPLE_DIR/coverage.tsv"
    awk '{sum+=$3; if($3>=10) c10++} END {if (NR>0) print sum/NR, c10/NR*100; else print 0, 0}' "$SAMPLE_DIR/coverage.tsv" > "$SAMPLE_DIR/depth_stats.tsv"
    read mean cov10 <<< $(cat "$SAMPLE_DIR/depth_stats.tsv")

    # === CONSENSUS FROM READS ===
    samtools mpileup -aa -A -Q 0 -d 10000 -f "$REF" "$SAMPLE_DIR/dedup.bam" | \
        ivar consensus -p "$SAMPLE_DIR/consensus.fa" -m 10 -n N -t 0.05

    # === VARIANT CALLING (.tsv only, no .vcf redirection) ===
    samtools mpileup -aa -A -Q 20 -d 10000 -f "$REF" "$SAMPLE_DIR/dedup.bam" | \
        ivar variants -p "$SAMPLE_DIR/variants" -m 10 -q 10 -t 0.05

    VAR_TSV="$SAMPLE_DIR/variants.tsv"

    if [[ ! -s "$VAR_TSV" || $(wc -l < "$VAR_TSV") -le 1 ]]; then
        echo "[INFO] No variants found for $sample."
        echo -e "$sample\t$mean\t$cov10\t0\t0\t0\t0\t0\tNA" >> "$SUMMARY_FILE"
        continue
    fi



  # === VARIANT STATS ===
total_vars=$(tail -n +2 "$VAR_TSV" | wc -l)
vars_q20=$(awk 'NR>1 && $12 >= 20' "$VAR_TSV" | wc -l)
vars_q30=$(awk 'NR>1 && $12 >= 30' "$VAR_TSV" | wc -l)
pct_q20=$(awk -v a="$vars_q20" -v t="$total_vars" 'BEGIN{if (t==0) print 0; else printf "%.2f", (a/t)*100}')
pct_q30=$(awk -v a="$vars_q30" -v t="$total_vars" 'BEGIN{if (t==0) print 0; else printf "%.2f", (a/t)*100}')

# === Ts/Tv Classification (only valid ACGT SNVs) ===
ts=$(awk 'NR > 1 && $4 ~ /^[ACGT]$/ && $5 ~ /^[ACGT]$/ &&
            (($4 == "A" && $5 == "G") || ($4 == "G" && $5 == "A") ||
             ($4 == "C" && $5 == "T") || ($4 == "T" && $5 == "C"))' "$VAR_TSV" | wc -l)

tv=$(awk 'NR > 1 && $4 ~ /^[ACGT]$/ && $5 ~ /^[ACGT]$/ &&
            !(($4 == "A" && $5 == "G") || ($4 == "G" && $5 == "A") ||
              ($4 == "C" && $5 == "T") || ($4 == "T" && $5 == "C"))' "$VAR_TSV" | wc -l)

tstv=$(awk -v ts="$ts" -v tv="$tv" 'BEGIN { if (tv == 0) print "NA"; else printf "%.2f", ts/tv }')

# === SUMMARY OUTPUT ===
echo -e "$sample\t$mean\t$cov10\t$total_vars\t$vars_q20\t$vars_q30\t$pct_q20\t$pct_q30\t$tstv" >> "$SUMMARY_FILE"


done 

echo -e "\n✅ Pipeline complete. Summary written to: $SUMMARY_FILE"