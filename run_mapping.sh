#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Single-end mapping with R1+R2 concatenated
#
# Outputs (per sample):
#   sample
#   raw_reads          # R1+R2
#   mapped_reads       # SE alignments to MPXV reference
#   mean_depth         # mean depth over genome
#   genome_ge_10x_pct  # % genome with depth ≥10x
#   Q20(%)             # from seqkit stats -a
#   Q30(%)             # from seqkit stats -a
#   mean_qual          # mean per-read avg PHRED (R1+R2)
#   median_qual        # median per-read avg PHRED (R1+R2)
###############################################################################

export LC_ALL=C
export LANG=C

THREADS=4
INPUT_DIR=""
REF_FA=""
OUTDIR="results-mpxv"

usage() {
    cat >&2 <<EOF
Usage: $0 -i <input_fastq_dir> -r <ref_fasta> [-o <outdir>] [-t <threads>]

  -i   Directory containing SAMPLE_R1.fastq(.gz) / SAMPLE_R2.fastq(.gz)
  -r   MPXV reference genome in FASTA (e.g. MPXV_ref.fasta)
  -o   Output directory (default: results-mpxv)
  -t   Number of threads for aligner/samtools (default: 4)
EOF
    exit 1
}

trap 'echo "[ERROR] Pipeline failed at line $LINENO. Check logs/ and seqkit/." >&2' ERR

# -------------------- parse arguments --------------------
while getopts "i:r:o:t:" opt; do
    case "$opt" in
        i) INPUT_DIR="$OPTARG" ;;
        r) REF_FA="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        *) usage ;;
    esac
done

[[ -z "$INPUT_DIR" || -z "$REF_FA" ]] && usage

# -------------------- choose aligner --------------------
ALIGNER_CMD=""
if command -v bwa-mem2 >/dev/null 2>&1; then
    ALIGNER_CMD="bwa-mem2"
    echo "[INFO] Using bwa-mem2 for alignment (R1+R2 concatenated, SE mode)."
elif command -v bwa >/dev/null 2>&1; then
    ALIGNER_CMD="bwa"
    echo "[INFO] Using bwa mem (R1+R2 concatenated, SE mode)."
else
    echo "[ERROR] Neither bwa-mem2 nor bwa found in PATH." >&2
    exit 1
fi

# -------------------- dependency checks --------------------
for cmd in seqkit samtools gawk; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "[ERROR] Required command '$cmd' not found in PATH: $cmd" >&2
        exit 1
    fi
done

# -------------------- prepare output structure --------------------
mkdir -p "$OUTDIR"/{bam,depth,summary,logs,seqkit}

SUMMARY_TSV="$OUTDIR/summary/mpxv_qc_for_manuscript.tsv"
echo -e "sample\traw_reads\tmapped_reads\tmean_depth\tgenome_ge_10x_pct\tQ20(%)\tQ30(%)\tmean_qual\tmedian_qual" > "$SUMMARY_TSV"

# -------------------- reference indexing --------------------
if [[ "$ALIGNER_CMD" == "bwa-mem2" ]]; then
    [[ -f "${REF_FA}.0123" ]] || bwa-mem2 index "$REF_FA"
else
    [[ -f "${REF_FA}.bwt"  ]] || bwa index "$REF_FA"
fi
[[ -f "${REF_FA}.fai" ]] || samtools faidx "$REF_FA"

# -------------------- main loop over samples --------------------
for R1 in "$INPUT_DIR"/*_R1.fastq*; do
    [[ -e "$R1" ]] || continue

    base=$(basename "$R1")

    if [[ "$base" == *_R1.fastq.gz ]]; then
        sample=${base%_R1.fastq.gz}
        ext=".fastq.gz"
    elif [[ "$base" == *_R1.fastq ]]; then
        sample=${base%_R1.fastq}
        ext=".fastq"
    else
        echo "[WARN] Cannot parse R1 filename: $base, skipping."
        continue
    fi

    R2="$INPUT_DIR/${sample}_R2$ext"
    if [[ ! -f "$R2" ]]; then
        echo "[WARN] R2 not found for sample $sample, skipping."
        continue
    fi

    echo "============================================================"
    echo "[INFO] Processing sample: $sample"
    echo "------------------------------------------------------------"

    ###########################################################################
    # 1) SEQKIT STATS (R1+R2) → raw_reads, Q20(%), Q30(%)
    ###########################################################################
    stats_file="$OUTDIR/seqkit/${sample}.stats.txt"
    seqkit stats -a "$R1" "$R2" > "$stats_file"

    # raw_reads = sum of num_seqs over R1+R2
    raw_reads=$(awk '
        NR==1 {
            for (i=1; i<=NF; i++) if ($i=="num_seqs") c_num=i
        }
        NR>1 {
            gsub(/,/, "", $c_num)
            sum += $c_num
        }
        END { print (sum+0) }
    ' "$stats_file")

    # Q20(%), Q30(%) – generic header-based extraction
    read Q20_pct Q30_pct < <(
        awk '
            NR==1 {
                for (i=1; i<=NF; i++) {
                    h=tolower($i)
                    if (index(h,"q20(")>0 && index(h,"%)")>0) c_q20=i
                    if (index(h,"q30(")>0 && index(h,"%)")>0) c_q30=i
                }
            }
            NR>1 {
                gsub(/,/, "", $c_q20)
                gsub(/,/, "", $c_q30)
                n    = $c_num + 0
                q20  = $c_q20 + 0.0
                q30  = $c_q30 + 0.0

                wQ20  += q20 * n
                wQ30  += q30 * n
                total += n
            }
            END {
                if (total==0) {
                    printf "0.00 0.00\n"
                } else {
                    printf "%.2f %.2f\n", wQ20/total, wQ30/total
                }
            }
        ' c_num="$(awk 'NR==1{for(i=1;i<=NF;i++)if($i=="num_seqs")print i}' "$stats_file")" "$stats_file"
    )

    ###########################################################################
    # 2) MEAN & MEDIAN PER-READ AVERAGE QUALITY (fx2tab --avg-qual)
    #
    # seqkit fx2tab --avg-qual prints a per-read mean PHRED.
    # We pool R1+R2, then compute:
    #   mean_qual   = mean of per-read avg-qual
    #   median_qual = median of per-read avg-qual
    ###########################################################################
    read mean_qual median_qual < <(
        seqkit fx2tab --avg-qual "$R1" "$R2" \
        | awk '
            NR>1 {
                q = $NF + 0.0   # avg-qual is last field
                total++
                sum += q
                qual[total] = q
            }
            END {
                if (total == 0) {
                    printf "0.00 0.00\n"
                    exit
                }
                mean = sum / total
                n = asort(qual)
                if (n % 2 == 1) {
                    med = qual[(n + 1) / 2]
                } else {
                    med = (qual[n/2] + qual[n/2 + 1]) / 2
                }
                printf "%.2f %.2f\n", mean, med
            }
        '
    )

    ###########################################################################
    # 3) MAPPING (R1+R2 concatenated, SE mode)
    ###########################################################################
    BAM="$OUTDIR/bam/${sample}.sorted.bam"
    LOG="$OUTDIR/logs/${sample}.${ALIGNER_CMD}.log"

    echo "[INFO] Mapping with $ALIGNER_CMD (R1+R2 concatenated)..."
    if [[ "$ALIGNER_CMD" == "bwa-mem2" ]]; then
        bwa-mem2 mem -t "$THREADS" "$REF_FA" <(cat "$R1" "$R2") \
            2> "$LOG" | samtools sort -@ "$THREADS" -o "$BAM"
    else
        bwa mem -t "$THREADS" "$REF_FA" <(cat "$R1" "$R2") \
            2> "$LOG" | samtools sort -@ "$THREADS" -o "$BAM"
    fi

    samtools index "$BAM"
    mapped_reads=$(samtools view -c -F 4 "$BAM")

    ###########################################################################
    # 4) DEPTH & % GENOME ≥ 10x
    ###########################################################################
    DEPTH_TSV="$OUTDIR/depth/${sample}.depth.tsv"
    echo "[INFO] Computing depth profile..."
    samtools depth -aa -d 0 "$BAM" > "$DEPTH_TSV"

    read mean_depth genome_ge_10x_pct < <(
        awk '
            {
                d_sum += $3
                if ($3 >= 10) c10++
                n++
            }
            END {
                if (n == 0) {
                    printf "0.00 0.00\n"
                } else {
                    mean_d = d_sum / n
                    pct10  = (c10 * 100.0) / n
                    printf "%.2f %.2f\n", mean_d, pct10
                }
            }
        ' "$DEPTH_TSV"
    )

    ###########################################################################
    # 5) WRITE SUMMARY LINE
    ###########################################################################
    printf "%s\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n" \
        "$sample" \
        "$raw_reads" \
        "$mapped_reads" \
        "$mean_depth" \
        "$genome_ge_10x_pct" \
        "$Q20_pct" \
        "$Q30_pct" \
        "$mean_qual" \
        "$median_qual" >> "$SUMMARY_TSV"

    echo "[INFO] Done: $sample"
done

echo "============================================================"
echo "[INFO] Summary table written to: $SUMMARY_TSV"
echo "============================================================"
