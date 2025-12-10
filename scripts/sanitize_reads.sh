#!/bin/bash

# Usage: ./sanitize_reads.sh -i <input_fastq_folder> -o <clean_output_folder>

set -e

usage() {
  echo -e "USAGE: $0 -i <input_fastq_folder> -o <output_clean_folder>\n"
  echo "This script uses BBMap's repair.sh to clean and synchronize paired FASTQ reads."
  exit 1
}

# === Parse arguments ===
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -i) IN_DIR="$2"; shift ;;
        -o) OUT_DIR="$2"; shift ;;
        -h|--help) usage ;;
        *) echo "[ERROR] Unknown parameter passed: $1"; usage ;;
    esac
    shift
done

[[ -z "$IN_DIR" || -z "$OUT_DIR" ]] && usage

mkdir -p "$OUT_DIR"

# === Loop through R1 files ===
for R1 in "$IN_DIR"/*.nH_R1.fastq*; do
    sample=$(basename "$R1" | sed 's/\.nH_R1\.fastq.*$//')
    R2="$IN_DIR/${sample}.nH_R2.fastq"
    [[ ! -f "$R2" ]] && R2="$IN_DIR/${sample}.nH_R2.fastq.gz"
    [[ ! -f "$R2" ]] && echo "[WARN] R2 missing for $sample, skipping" && continue

    echo "[INFO] Sanitizing $sample"

    # Output names
    CLEAN_R1="$OUT_DIR/${sample}.nH_clean_R1.fastq"
    CLEAN_R2="$OUT_DIR/${sample}.nH_clean_R2.fastq"
    SINGLETONS="$OUT_DIR/${sample}.singleton.fastq"

    # Repair reads
    repair.sh in1="$R1" in2="$R2" out1="$CLEAN_R1" out2="$CLEAN_R2" outs="$SINGLETONS" overwrite=t
done

echo -e "\n All FASTQ files sanitized and stored in: $OUT_DIR"
