# MPXV-Pooling-VSP2-FieldStudy

This repository supports the manuscript:

**“Evaluating High-Plex Sample Pooling for the Illumina Viral Surveillance Panel v2 in Mpox Virus Genomic Detection under Resource-Constrained Conditions”**
It includes bioinformatics workflows, scripts, metadata, and processed results from an evaluation of sample multiplexing strategies during the 2024 Mpox outbreak in Burundi. The study compares sequencing performance across different pooling levels (3, 6, 9, and 15 samples per capture) using the Illumina VSP2 panel and iSeq 100 platform.

---

## 📂 Repository Contents

| Folder/File         | Description                                                                 |
|---------------------|-----------------------------------------------------------------------------|
| `data/`             | Sample metadata, Ct values, genome coverage, and variant tables             |
| `scripts/`          | Shell scripts and Python notebooks for preprocessing, mapping, and analysis |
| `results/`          | Processed outputs including variant metrics, assembly stats, and coverage   |
| `figures/`          | Publication-ready plots (PNG and PDF formats)                               |
| `MANUSCRIPT.md`     | Draft manuscript or link to final publication                               |
| `LICENSE`           | Licensing information (MIT by default)                                      |

---

## Bioinformatics Workflow Overview

This workflow processes Mpox VSP iSeq100 paired-end sequencing data across different pooling levels (3, 6, 9, and 15 samples per run). The aim is to evaluate the feasibility of cost-efficient genomic surveillance in low-resource contexts by analyzing genome coverage, variant detection, and consensus sequence generation.

##  Methods Overview

- **Preprocessing**: FastP v0.20.1 for quality filtering
- **Host Depletion**: Kraken2, SeqKit filtering
- **Alignment**: BWA-MEM v0.7.17
- **Variant Calling**: iVar, bcftools
- **Assembly**: SPAdes v3.15, QUAST
- **Taxonomic Classification**: KrakenUniq, Kaiju, Centrifuge

---

## 1. 🔍 Raw Read Quality Control  
**Tool:** `FastQC` (v0.11.9)  
**Input:** Raw paired-end FASTQ files  

**Purpose:**  
- Assess base quality distribution  
- Detect adapter contamination  
- Evaluate per-base sequence content  
- Identify technical issues prior to trimming  

---

## 2. Adapter Trimming & Quality Filtering  
**Tool:** `fastp` (v0.23.4)

**Key parameters:**  
- Automatic adapter detection  
- Sliding-window Q20 trimming  
- Minimum length: 50 bp  
- HTML + JSON QC reports

**Purpose:** Improve read quality, remove adapters, and filter out very short or low-quality reads.

---

## 3. 🔁 Post-trimming QC  
**Tool:** `FastQC` (on trimmed FASTQs)

**Purpose:** Confirm improvement after trimming and detect any remaining quality issues.

---

## 4. 📊 Consolidated QC Summary  
**Tool:** `MultiQC`  
**Setting:** `TZ=UTC` to avoid timestamp errors  

**Output:** One consolidated HTML report summarizing all FastQC and fastp metrics.

---

## 5.  Host-Read Removal  
**Tools:**  
- `Kraken2` (v2.1.2)  
- `src_scrubber` (v0.4.0)  

**Database:** Standard Kraken2 DB with human genome (GRCh38) and viral sequences  

**Process:**  
- Taxonomic classification with `kraken2`  
  - `--paired` mode  
  - `--confidence 0.1`  
  - `--report` and `--classified-out` flags  
- Host-read filtering using `src_scrubber`  
  - Filters out reads classified as human or ambiguous taxa  
  - Retains high-confidence viral reads  

**Output:**  
- Cleaned, host-depleted paired FASTQ files ready for downstream alignment  

---

## 6. 🧭 Alignment to Mpox Reference Genome  
**Tool:** `BWA-MEM2` (v2.2.1)  
**Reference:** MPXV Clade IIb (e.g., ON563414.3)

**Process:**  
- Alignment (`bwa-mem2 mem`)  
- Conversion & sorting (`samtools view`, `samtools sort`)  
- Indexing (`samtools index`)

**Outputs:**  
- Sorted BAM  
- BAM index  

---

## 7. 📏 Genome Coverage Estimation  
**Tool:** `samtools depth`

**Purpose:**  
- Compute average depth  
- Identify low-coverage genome segments  
- Generate depth profiles for comparing pooling levels  

---

## 8. 🧬 Variant Calling  
**Tools:**  
- `samtools mpileup`  
- `ivar variants` (v1.4.3)

**Parameters:**  
- MQ ≥ 20  
- BQ ≥ 20  
- Variant frequency ≥ 0.25  
- Minimum depth = 10  

**Output:** `.tsv` variant tables

**Purpose:** Detect SNPs and indels across pooling strategies.

---

## 9. 🧬 Consensus Sequence Generation  
**Tool:** `ivar consensus`

**Parameters:**  
- Base frequency threshold: 0.60  
- Depth cutoff: 10  
- Mask low-depth bases as "N"

**Output:** Consensus FASTA for each sample.

---

## 10. 📈 Variant Metrics & Summaries  
**Script:** `variant_summary.py`

**Purpose:**  
- Compile variant counts per sample  
- Distinguish high-confidence (QUAL ≥ 20) and very-high-confidence (QUAL ≥ 30) variants  
- Compare variant retention across pooling levels

---

## 11. 📂 Output Directory Structure  

```
results/
├── QC/
│   ├── raw_fastqc/
│   ├── trimmed_fastqc/
│   └── multiqc_report.html
├── trimmed_reads/
├── alignments/
│   ├── *.bam
│   └── *.bai
├── coverage/
│   └── depth.txt
├── variants/
│   └── *.tsv
├── consensus/
│   └── *.fa
├── summary/
│   └── variant_metrics.tsv
```

---

## 12. 🧪 Interpretation & Use Case  
This workflow supports:  
- Evaluating depth loss at increasing pooling levels  
- Assessing variant detection sensitivity  
- Quantifying completeness of consensus genomes  
- Supporting public health decisions on cost-efficient genomic surveillance  
- Reproducible benchmarking of sequencing strategies for field and mobile-lab settings


## Highlights

- MPXV enrichment remains robust in pools ≤6
- Larger pools reduce coverage, variant yield, and assembly completeness
- High Ct (low viral load) samples underperform in high-pooling conditions
- SPAdes assembly quality degrades significantly beyond 6-sample pools

---





